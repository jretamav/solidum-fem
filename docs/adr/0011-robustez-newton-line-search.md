# ADR 0011 — Robustez de Newton-Raphson: line search y telemetría de divergencia

- **Estado**: aceptado con enmienda (ver §"Enmiendas tras implementación")
- **Fecha**: 2026-05-18
- **Alcance**: `NonlinearSolver` y `NewtonNewmarkSolver` en `solidum/math/solvers/`; mensaje de `RuntimeError` cuando el solver diverge; catálogo de solvers.

## Contexto

La auditoría de la [fase A](../auditorias/solvers_robustez_fase_A.md) (11 tests en
`tests/test_solver_robustness.py`) confirmó que los solvers no lineales son más
robustos de lo que sugerían los huecos de validación. De los regímenes ejercitados,
solo uno produce un fallo no atribuible a una limitación física genuina:

**Newton oscila entre dos estados sin converger** cuando el incremento completo
de Newton aterriza fuera del pozo de convergencia. Reproducido en Drucker-Prager
perfectamente plástico (`H=0`) con carga ~50× la cohesión característica: el
residuo alterna entre ratios 2.1 y 440 durante las 30 iteraciones permitidas,
el paso se biseca, el rebote se repite, y el solver muere al alcanzar
`min_delta_lambda` con mensaje `"el solver ha divergido"`.

El patrón es canónico y bien estudiado en la literatura: el incremento `δU`
predicho por la jacobiana es matemáticamente correcto pero numéricamente
agresivo. La cura estándar es **escalar el incremento por un factor `α ∈ (0, 1]`
elegido para garantizar decrecimiento del residuo** — es lo que se conoce como
*line search* o *globalización* del Newton-Raphson.

Síntoma secundario: el mensaje de `RuntimeError` al divergir no distingue entre
modos de fallo cualitativamente distintos. Hoy se confunden en un único error:

1. Newton fuera del pozo de convergencia, rescatable con line search.
2. Singularidad real de `K_t` (problema físico mal planteado).
3. Carga > capacidad plástica del modelo (sin solución estática).
4. Tolerancia mal calibrada para la escala del problema.

Sin distinción tipada, todo se ve como "no convergió" — el diagnóstico recae
enteramente en leer el log de iteraciones.

## Decisión

Se introducen dos mejoras coordinadas: **line search Armijo con backtracking**
y **telemetría de divergencia tipada** en el `RuntimeError`.

### Line search Armijo en `NonlinearSolver` y `NewtonNewmarkSolver`

Tras computar el incremento de Newton `δU` y antes de aceptarlo, se busca
`α ∈ (0, 1]` que satisfaga la condición de Armijo (suficiente decrecimiento):

```
‖R(U_iter + α·δU)‖ ≤ (1 − c₁·α) · ‖R(U_iter)‖
```

con `c₁ = 1e-4` (valor canónico, Nocedal-Wright §3.1). Backtracking simple:
arranca con `α = 1.0`, y si la condición no se cumple, se reduce
`α ← ρ·α` con `ρ = 0.5` hasta encontrar uno válido o agotar 10 retrocesos.
Si se agotan, se acepta el último `α` y el comportamiento degrada al Newton
estándar (no se aborta — se confía en que el siguiente paso del solver
(bisección de Δλ, paso temporal, etc.) decida si rescatar o fallar).

El caso `α = 1.0` es el primer intento; en problemas donde Newton ya converge
hoy, el line search **no introduce iteraciones adicionales ni cambia el
resultado**: la condición de Armijo se satisface en el primer intento y se
sigue por el camino actual. El sobrecoste es **una evaluación adicional de
`F_int(U_iter + δU)` por iteración convergente**, lo que con tangente
consistente y problemas pequeños es marginal.

`ArcLengthSolver` **queda fuera del alcance**: el método de Crisfield ya
incorpora una forma implícita de globalización vía la ecuación cuadrática
de restricción, que selecciona la raíz menos abrupta. Aplicar Armijo sobre
los dos `δU` del arco (`du_R` y `du_t`) por separado rompería la geometría
del arco. Si emerge un caso real donde el arc-length oscila, se evaluará
entonces una variante específica.

#### API pública

Constructor de `NonlinearSolver` y `NewtonNewmarkSolver` admite un parámetro
opcional:

```python
def __init__(self, assembler, ..., line_search: bool = True):
```

Default `True` porque rescata casos sin sobrecoste perceptible cuando no
es necesario. El usuario puede desactivarlo con `line_search=False` para
diagnóstico o para reproducir comportamiento legacy. Sin parámetros de
configuración adicionales (`c₁`, `ρ`, `max_backtracks`) en la API pública —
viven como constantes en `solidum/constants.py` con los valores canónicos.

### Telemetría de divergencia tipada

Se introduce un módulo nuevo `solidum/math/solvers/diagnostics.py` con un tipo
de excepción discriminado:

```python
class SolverDivergedError(RuntimeError):
    """Base para fallos de convergencia con diagnóstico tipado."""

    def __init__(self, *, mode: str, last_residual: float, last_delta: float,
                 last_load_factor: float, n_bisections: int, hint: str):
        ...

class OscillatingNewtonError(SolverDivergedError):
    mode = "oscillating"
    hint = ("El residuo no decrece monótonamente entre iteraciones consecutivas. "
            "Probable rebote fuera del pozo de convergencia. Considera activar "
            "line_search=True (si no lo está) o reducir la carga.")

class SingularTangentError(SolverDivergedError):
    mode = "singular_tangent"
    hint = ("La rigidez tangente perdió rango. Probable bifurcación, punto "
            "límite o problema físico mal planteado. Considera ArcLengthSolver.")

class LoadExceedsCapacityError(SolverDivergedError):
    mode = "load_exceeds_capacity"
    hint = ("El residuo se estabiliza por encima de la tolerancia tras muchas "
            "iteraciones, sin oscilación ni divergencia. Probable carga "
            "superior a la capacidad plástica del modelo.")
```

Detección del modo dentro del solver:

- **`oscillating`**: el residuo de la iteración `k+1` es mayor que el de `k`
  durante ≥ 3 iteraciones consecutivas, **con line search activado** (sin
  line search, oscilación leve es esperable y no diagnóstica). Captura el
  caso del Drucker-Prager perfectamente plástico documentado.
- **`singular_tangent`**: el solver lineal devolvió `RuntimeError` por
  matriz singular, o Cholesky degradó a LU **y además** el residuo crece.
- **`load_exceeds_capacity`**: residuo estable (variación < 1% entre
  iteraciones) durante ≥ 5 iteraciones, sin cruzar la tolerancia. Heurística;
  vale más en plasticidad perfecta que en daño.
- **Modo `unknown`**: el solver agotó iteraciones sin caer en ningún
  patrón claro. Comportamiento por defecto, equivalente al `RuntimeError`
  actual pero con telemetría incluida.

El campo `hint` se imprime en el mensaje del error junto con las métricas
numéricas. No es prescriptivo (no instruye al usuario qué hacer); es
informativo y enlaza al diagnóstico probable.

### Actualización del catálogo de solvers

[`docs/catalogo_solvers.md`](../catalogo_solvers.md) se actualiza en el mismo
commit del ADR para reflejar:

- `NonlinearSolver` y `NewtonNewmarkSolver`: nuevo parámetro `line_search`.
- `NonlinearSolver`: la línea actual *"no puede atravesar puntos límite
  (snap-through) ni recorrer ramas con derivada negativa"* se matiza a
  *"con bisección adaptativa y cinemática no lineal (corotacional), puede
  atravesar puntos límite suaves; no atraviesa snap-back con `du/dλ < 0`,
  para eso se requiere `ArcLengthSolver`"*. Fundamentado por el test
  `test_snap_through_with_load_control` de la fase A.

## Consecuencias

**Inmediatas**:

- El régimen documentado en `test_drucker_prager_overload_oscillates`
  cambia de modo: en lugar de oscilar y morir con bisección agotada, el
  line search detecta el patrón y rescata con `α` reducido (o, si la carga
  realmente excede la capacidad, falla con `LoadExceedsCapacityError`
  informativo en lugar de `RuntimeError` genérico).
- Problemas donde Newton ya converge hoy no cambian de resultado ni de
  velocidad observable (Armijo se satisface con `α=1` en la primera
  comprobación).
- Diagnóstico de divergencia pasa de "leer 30 líneas de log" a "leer el
  tipo de excepción y su hint".

**Hacia adelante**:

- HHT-α ya consume el mismo patrón: `NewtonHHTSolver` (2026-05-18,
  ADR 0009 variante de fase 4) hereda el `line_search` opcional y la
  jerarquía de excepciones tipadas de `NewtonNewmarkSolver` sin tocar
  el módulo `solidum/math/solvers/diagnostics.py`. Otros solvers no
  lineales futuros (Riks generalizado, BFGS) consumirán el mismo
  patrón.
- El modo `oscillating` con line search activado y aun así divergiendo es
  señal de que el problema requiere una globalización más fuerte (trust
  region, Cauchy mixto). Marca natural para retomar.

**Migración**:

- API: `NonlinearSolver(...)` y `NewtonNewmarkSolver(...)` aceptan el nuevo
  kwarg con default. Código existente que no lo pase obtiene line search
  activado — cambio de comportamiento en casos patológicos, neutro en el
  resto.
- Tests existentes: ninguno cambia de resultado. Validado por la fase C
  (cuando se ejecute) corriendo la suite completa antes y después.
- Test `test_drucker_prager_overload_oscillates` se reformula: en lugar
  de esperar `RuntimeError`, espera `LoadExceedsCapacityError` o
  convergencia con line search (a determinar empíricamente).

## Alternativas consideradas

- **Damping fijo** (`α = 0.5` siempre). Más simple que Armijo (sin búsqueda),
  pero impone sobrecoste innecesario en problemas donde Newton convergería
  con `α=1`. Rechazado: peor relación coste/beneficio que Armijo, mismo
  esfuerzo de implementación.

- **Trust region (Levenberg-Marquardt-style)**. Globalización más fuerte que
  Armijo, controla la dirección además del paso. Costoso de implementar
  bien (requiere resolver `(K + λI)·δU = R` con `λ` adaptativo) y
  desproporcionado para el tipo de problemas que la fase A documenta.
  Postergado hasta que aparezca un caso donde Armijo no baste.

- **Newton-Krylov inexacto** (resolver el sistema lineal solo aproximadamente,
  reusar Krylov entre iteraciones). Orientado a problemas grandes donde
  factorizar `K_t` cada iteración es prohibitivo. Solidum hoy resuelve mallas
  pequeñas-medianas con factorización directa; over-engineering.

- **Aplicar line search también a `ArcLengthSolver`**. Considerado y
  descartado por la razón geométrica explicada en §"Decisión": el arc-length
  ya tiene una globalización implícita; Armijo sobre los dos `δU` del arco
  rompería la geometría. Si emerge necesidad real, se evaluará variante
  específica (line search sobre `dU_iter` total, no sobre las componentes).

- **Detectar oscilación sin line search, solo reportar mejor el error**.
  Detección útil pero pasiva. Con line search disponible, mejor *rescatar*
  que solo *reportar*. La detección queda igualmente, como criterio del
  modo `oscillating` cuando el rescate también falla.

## Deuda heredada

- **`ArcLengthSolver` sin globalización explícita**. La justificación es
  geométrica, no fundamental — si aparece un caso de oscilación en
  arc-length, este ADR no lo cubre. Se decidirá en un ADR posterior si
  emerge.

- **Detección heurística del modo `load_exceeds_capacity`**. La heurística
  ("residuo estable durante ≥5 iter") es pragmática, no rigurosa. Casos
  donde el residuo se estanca por otras razones (mal escalado severo,
  acoplamiento no diagonal) podrían misclasificarse. Aceptable porque el
  `hint` invita al usuario a verificar, no prescribe.

- **Constantes del line search en `solidum/constants.py`, no configurables**.
  `c₁ = 1e-4`, `ρ = 0.5`, `max_backtracks = 10` son canónicos en la
  literatura (Nocedal-Wright). Si un caso real requiere ajuste fino,
  promover a parámetros de constructor en una segunda iteración.

## Referencias

- Nocedal, J. & Wright, S. (2006). *Numerical Optimization*, 2nd ed., §3.1
  (line search methods), §3.5 (Armijo condition and backtracking).
- Dennis, J. E. & Schnabel, R. B. (1983). *Numerical Methods for
  Unconstrained Optimization and Nonlinear Equations*, §6.3.
- Crisfield, M. A. (1991). *Non-linear Finite Element Analysis of Solids
  and Structures*, vol. 1, cap. 9 — discusión del line search en MEF.
- ADR 0007 — Tolerancias de convergencia (criterio que consume este
  solver; line search no lo modifica).
- ADR 0009 — Análisis modal y dinámico (origen de `NewtonNewmarkSolver`
  afectado por este ADR).
- Auditoría fase A — [`docs/auditorias/solvers_robustez_fase_A.md`](../auditorias/solvers_robustez_fase_A.md).
- Grippo, L., Lampariello, F. & Lucidi, S. (1986). *A nonmonotone line
  search technique for Newton's method*, SIAM J. Numer. Anal. 23(4),
  707–716. — origen de la variante de descenso no monótono adoptada en
  la enmienda.
- de Borst, R. & Sluys, L. J. (1999). *Computational Methods in
  Non-linear Solid Mechanics*, TU Delft, §3.3-3.5 (Newton-Raphson,
  line-search, criterios de convergencia) y §3.7 (ejemplo del shallow
  truss von Mises). En [`docs/referencias/`](../referencias/). Tabla 3.1
  documenta la inutilidad del line-search con full Newton-Raphson y
  tangente consistente — respalda la enmienda del default.

## Enmiendas tras implementación

### Enmienda 1 (2026-05-18) — Default cambiado a `False`

**Decisión original**: ``line_search: bool = True`` por defecto, con el
argumento de que "rescata casos sin sobrecoste perceptible cuando no es
necesaria (α=1 satisface Armijo en el primer intento si Newton ya converge)".

**Hallazgo empírico**: la afirmación resultó falsa. Al activar el line
search por defecto, rompió **3 tests de daño 2D** que pasaban antes
(`test_damage_active_converges_in_few_iter`,
`test_damage_progresses_under_increasing_load`,
`test_unloading_no_damage_increment`). El patrón: en problemas con tangente
consistente cuasi-cuadrática y materiales con historia (daño activo,
plasticidad cerca de la transición elástico→postcrítico), F_int(U) no es
monótono respecto a U a lo largo de la dirección δU. El residuo a veces
**sube transitoriamente** durante el avance de Newton, que sin embargo
está yendo en la dirección correcta y converge en pocas iteraciones más.
Cualquier line search que exija "R no debe subir" (Armijo, descenso no
monótono GLL simplificado) **rechaza esos pasos correctos** y mete al
solver en un bucle de backtracking → α minúsculos → R apenas se mueve →
bisección de Δλ → repetición.

Probadas dos variantes:

1. **Armijo clásico** con `c₁=1e-4`: rechazo agresivo, oscilación severa.
2. **Descenso no monótono GLL** (`R_after ≤ R_before` sin término de
   decrecimiento): rechazo menos agresivo, pero el mismo patrón
   subsiste en daño 2D activo. Newton retrocede a α≈0.008 cada
   iteración, R baja marginalmente, repetición.

**Decisión enmendada**: ``line_search: bool = False`` por defecto. La
infraestructura del line search y de las excepciones tipadas se conserva
intacta — el usuario activa ``line_search=True`` cuando observa
oscilación en su problema (caso documentado: Drucker-Prager perfectamente
plástico con sobrecarga).

**Verificación empírica del beneficio cuando es necesario**: comparativa
en el caso `test_drucker_prager_overload_oscillates` (carga 50×
cohesión, ``H=0``):

| Modo | Resultado | Último ‖R‖ |
|---|---|---|
| ``line_search=False`` | diverge | 5.18e-1 |
| ``line_search=True``  | diverge (carga excede capacidad) | 5.89e-3 |

El line search baja el residuo dos órdenes de magnitud aunque no rescata
la convergencia (físicamente no hay solución estática para esta carga).
Es decir: el line search **sí es útil cuando hay oscilación**, pero
**dañino cuando Newton convergería sin él**. La separación que la
enmienda introduce ("activar solo cuando hace falta") refleja esta
realidad.

**Razón conceptual del fracaso de la afirmación original**: el ADR
asumía que la condición Armijo se satisface trivialmente con α=1 cuando
Newton converge cuadráticamente. Esto es cierto en problemas
**asintóticamente cuadráticos cerca de la solución**, pero falso en la
**fase inicial** del Newton para problemas no lineales con materiales
con historia: las primeras iteraciones del Newton aterrizan
deliberadamente en regiones de residuo creciente para acomodar la
evolución del state interno (kappa, eps_p), y la convergencia cuadrática
solo emerge tras 1-2 iteraciones de "exploración".

**Variante GLL adoptada en el helper**: se mantiene la condición
``R_after ≤ R_before`` (sin término Armijo de suficiente decrecimiento)
en el helper, porque Armijo puro era aún más agresivo. La variante GLL
es la implementación efectiva cuando el usuario activa ``line_search=True``.

**Sin cambio en el resto del ADR**: la decisión arquitectural
(infraestructura disponible, excepciones tipadas, separación de modos de
divergencia) se mantiene intacta. Solo cambia el default.

**Respaldo bibliográfico (de Borst & Sluys, TU Delft 1999, §3.4)**: el
hallazgo experimental coincide con la conclusión explícita del libro
de texto canónico de la materia. Tabla 3.1 (full Newton-Raphson con
plasticidad von Mises sobre el cilindro Lo-Scordelis): el line-search
con `ψ ∈ {0.4, 0.6, 0.8}` cuesta **102% del CPU** y solo ahorra 2
iteraciones (34 vs 36 sin line-search). Cita textual de p. 51:

> *"Line-searches are only then useful when no proper tangent stiffness
> relation is adopted. In case of a full Newton-Raphson method no
> savings in computer time are obtained when applying line-searches.
> ... line-searches can only enhance the performance of a full
> Newton-Raphson scheme in the very first iteration, in order to bring
> the solution within the radius of convergence of the Newton-Raphson
> method."*

Solidum usa **full Newton-Raphson** con **tangente consistente**
(materiales con `compute_state` que devuelve tangente algorítmica:
`Elastoplastic1D`, `VonMises2D`, `DruckerPrager2D`, `IsotropicDamage1D/2D`).
La conclusión del libro aplica directamente: el line-search es
contraproducente como default y solo aporta cuando el Newton inicial
está fuera del radio de convergencia (caso del Drucker-Prager
perfectamente plástico con sobrecarga, donde activarlo explícitamente
sí baja el residuo dos órdenes de magnitud).

Esto convierte la enmienda en una decisión bien fundada bibliográficamente,
no solo empíricamente.

