# ADR 0007 — Tolerancias de convergencia en solvers no lineales

- **Estado**: aceptado; enmendado el 2026-11-27 (piso de redondeo y calibración sin fuerzas, ver «Enmienda» al final)
- **Fecha**: 2026-05-13
- **Alcance**: `NonlinearSolver` y `ArcLengthSolver` en `solidum/math/solvers.py` (hoy `solidum/math/solvers/`, con Newmark, HHT y el corrector compartido del ADR 0015); constantes en `solidum/constants.py`; tests de convergencia en `tests/`.

## Contexto

Los dos solvers no lineales (`NonlinearSolver` con Newton-Raphson incremental y `ArcLengthSolver` con longitud de arco cilíndrica de Crisfield) usan hoy un criterio puramente relativo:

```
err_disp  = ‖δU‖     / (‖U_iter‖ + ZERO_TOL)
ref_force = max(‖F_ext_step‖, ‖F_int‖, ZERO_TOL)
err_force = ‖R[free]‖ / ref_force
error     = max(err_disp, err_force)
converged ↔ error < self.tol
```

con `self.tol = CONVERGENCE_TOL = 1e-5` y `ZERO_TOL = 1e-12`. El esquema funciona en los problemas actuales, pero arrastra **tres limitaciones estructurales** análogas a las que motivaron el ADR 0006:

1. **No hay tolerancia absoluta real.** `ZERO_TOL` solo actúa como piso del denominador (`ref_force`), no como término absoluto sobre la comparación. En regímenes donde la carga característica colapsa transitoriamente — primer sub-incremento con `delta_lambda·F_ext ≈ 0`, fases de carga balanceada, problemas con cargas que se cancelan — el ruido numérico hace que `err_force` salte erráticamente entre cero y valores grandes sin reflejar fallo real de equilibrio. La cura conceptual es exactamente la del ADR 0006: `‖R‖ ≤ atol + rtol · escala`, con `atol` en unidades físicas (fuerza para residuo, longitud para corrección de desplazamiento).

2. **Una sola `tol` para dos criterios físicamente distintos.** `err_disp` y `err_force` son ambos adimensionales pero miden cosas diferentes: ratio de "incrementabilidad" del iterado vs ratio de "equilibrio". En problemas casi-rígidos `err_disp` baja vertiginosamente y `err_force` baja despacio; en plástico perfecto puede ser al revés. Forzar el mismo umbral hace que el criterio fácil sature trivialmente y el difícil mande siempre. Los códigos comerciales lo separan: Abaqus `RTOL`/`CTOL`, OpenSees `tolForce`/`tolDispl`, Code_Aster `RESI_GLOB_RELA` con criterios independientes.

3. **Fórmula replicada en dos solvers.** `NonlinearSolver` y `ArcLengthSolver` repiten el cómputo del criterio a mano. Mismo problema que tenía `VonMises2D` antes del centralizado de `admissibility_tol`: la política vive en dos sitios, cambiarla obliga a recordar tocar ambos.

## Decisión

Se adopta el patrón **`atol + rtol · escala(estado)`** del ADR 0006, separado por criterio y centralizado en un único punto de la base de los solvers no lineales.

### Política del criterio

Convergencia simultánea de los dos criterios (AND, como hoy — no OR):

```
‖R[free_dofs]‖ ≤ atol_force + rtol_force · ref_force(U_iter, F_ext, F_int)
‖δU‖          ≤ atol_disp  + rtol_disp  · ref_disp(U_iter)
```

> **Nota 2026-09-30 (deuda #31).** `R[free_dofs]` se escribió pensando en condiciones de Dirichlet. La norma que entra al criterio es la del equilibrio reducido del ADR 0004, `‖TᵀR‖`, que con sólo Dirichlet coincide con `‖R[free_dofs]‖`; con restricciones multipunto la fuerza de cada esclavo pasa a sus maestros y `R[free_dofs]` no se anula en la solución. Los solvers estáticos la usan desde esta fecha, como ya hacía `NewtonNewmarkSolver`. La escala `max(‖F_ext‖, ‖F_int‖)` sigue usando `F_int` completo, que incluye las fuerzas de las restricciones —reacciones de Dirichlet y fuerzas de las MPC sobre esclavos y maestros—, igual que aquí se decidió para las reacciones; sólo afecta a la tolerancia (auditoría de las MPC del mismo día).

con las escalas evaluadas en el estado corriente de la iteración:

- `ref_force = max(‖F_ext_step‖, ‖F_int‖)`. Sin `ZERO_TOL` añadido — `atol_force` ya es el piso absoluto.
- `ref_disp  = ‖U_iter‖`. Sin `ZERO_TOL` añadido — `atol_disp` ya es el piso absoluto.

Las normas son L2 sobre DOFs libres (igual que hoy: prescritos contienen reacciones, no fallos de equilibrio).

### Constantes nuevas en `solidum/constants.py`

Todas adimensionales: las tolerancias relativas se aplican directamente al ratio, y las "absolutas" entran como **factores** que se multiplican por escalas autoderivadas del problema en su primera evaluación. Esto evita codificar valores con unidades en una constante global — el código deja de presuponer que el usuario trabaja en N/m vs MPa/mm.

```python
# Tolerancias relativas (rtol del patrón atol + rtol · escala):
CONVERGENCE_RTOL_FORCE = 1.0e-5
CONVERGENCE_RTOL_DISP  = 1.0e-5

# Factores absolutos adimensionales. El atol efectivo se autoderiva como
# `atol = factor · escala_inicial`, donde escala_inicial se computa una vez
# al iniciar solve() a partir del primer ensamblaje. Funcionan como piso
# de la comparación cuando ref(estado) colapsa transitoriamente.
CONVERGENCE_ATOL_FORCE_FACTOR = 1.0e-9
CONVERGENCE_ATOL_DISP_FACTOR  = 1.0e-9
```

`CONVERGENCE_TOL = 1e-5` y `ZERO_TOL` (en su rol de piso del denominador) se retiran del módulo `solvers.py` en el mismo commit que introduce el nuevo criterio. `ZERO_TOL` se conserva en `constants.py` por otros usos legítimos del proyecto (arc-length sign, etc.) pero deja de tener semántica de tolerancia de convergencia.

### Escalas autoderivadas

- `escala_force_inicial = ‖F_ext_global‖`, evaluada una sola vez al iniciar `solve()` antes del primer paso. En análisis con carga nula al inicio (control en desplazamiento puro), se autoderiva de la reacción del primer ensamblaje: `escala_force_inicial = max(‖F_ext‖, ‖F_int_initial‖)`.
- `escala_disp_inicial = escala_force_inicial / ‖K_initial‖_∞`, aproximación de orden 1 del desplazamiento característico esperado. Se computa también una vez tras el primer ensamblaje.

Estas escalas iniciales producen los `atol` efectivos:

```
atol_force = CONVERGENCE_ATOL_FORCE_FACTOR · escala_force_inicial
atol_disp  = CONVERGENCE_ATOL_DISP_FACTOR  · escala_disp_inicial
```

Ventaja: el código es **completamente invariante bajo cambio de unidades**, no solo en el término relativo. El usuario no tiene que ajustar nada al pasar de N/m a kN/mm.

> **Enmienda 2026-11-27.** La invariancia no se cumplía con control en desplazamiento desde el reposo: en el primer ensamblaje aún no se ha impuesto el desplazamiento, `F_int` es nula y la implementación usaba una unidad de fuerza (`max(‖F_ext‖, ‖F_int‖, 1)`). Ver «Enmienda» al final: la escala nula es ahora válida y el piso lo da un término de redondeo calculado con el modelo.

### Semántica de convergencia: AND, no OR

Convergencia simultánea de los dos criterios (AND). Es la elección canónica en código de investigación (Bathe, Crisfield, Owen-Hinton). OR es más permisivo y se ve en algunos códigos comerciales donde el vendor relajó el default para reducir tickets de soporte — incompatible con el rigor declarado en Reglas.md §0.

### API del solver: ruptura limpia, sin atajo legacy

Los constructores de `NonlinearSolver` y `ArcLengthSolver` reemplazan el parámetro `tol=` por uno solo de configuración:

```python
def __init__(self, assembler,
             convergence: Optional[ConvergenceCriterion] = None,
             ...)
```

Si `convergence` es `None`, se instancia con los defaults del proyecto. No se acepta `tol=` como atajo: la presencia de un único número escondería la separación entre fuerza y desplazamiento y entre `atol` y `rtol`, que es justamente lo que el ADR introduce. Los YAML y tests existentes se actualizan en el mismo commit. Coste de migración bajo (proyecto de un solo usuario).

### Punto único de la política: clase `ConvergenceCriterion`

Se introduce un módulo nuevo `solidum/math/convergence.py` con una clase pequeña que encapsula configuración + estado calibrado + evaluación:

```python
@dataclass
class ConvergenceState:
    converged: bool
    ratio_force: float   # ‖R‖ / tol_force, < 1 ⇒ ok
    ratio_disp: float    # ‖δU‖ / tol_disp,  < 1 ⇒ ok

class ConvergenceCriterion:
    def __init__(self,
                 rtol_force=CONVERGENCE_RTOL_FORCE,
                 rtol_disp =CONVERGENCE_RTOL_DISP,
                 atol_force_factor=CONVERGENCE_ATOL_FORCE_FACTOR,
                 atol_disp_factor =CONVERGENCE_ATOL_DISP_FACTOR):
        ...

    def calibrate(self, F_ext_global: np.ndarray,
                  F_int_initial: np.ndarray,
                  K_norm: float) -> None:
        """Computa atol_force y atol_disp a partir de las escalas
        del problema. Se invoca una vez al inicio de solve()."""
        ...

    def evaluate(self, residual_norm: float, ref_force: float,
                 delta_u_norm: float, u_norm: float) -> ConvergenceState:
        """Aplica la política atol + rtol · escala separada por criterio.
        Devuelve estado + ratios normalizados para logging."""
        ...
```

Razón para usar clase pequeña en lugar de función libre: la configuración (cuatro tolerancias) es estable durante toda la corrida pero los `atol` requieren una fase de calibración tras el primer ensamblaje. Una clase encapsula configuración + estado calibrado naturalmente; una función libre obligaría a pasar 6-8 argumentos en cada llamada y a externalizar la lógica de calibración.

Razón para no crear `NonlinearSolverBase`: los dos solvers comparten *solo* el criterio de convergencia. NR y arc-length tienen estructuras radicalmente distintas (incremental-iterativo vs predictor-corrector con ecuación cuadrática de restricción). Una clase base forzada sería abstracción especulativa, contraria a Reglas.md §1. Cuando aparezca un tercer solver no lineal y se vea qué comparten realmente, se evalúa extraer una base.

## Consecuencias

**Inmediatas**:

- Misma física en N, kN o MN converge con la misma precisión efectiva. Cubierto por test `TestConvergenceUnitInvariance` (a escribir).
- Régimen degenerado de carga inicial casi nula deja de oscilar espuriamente: el piso absoluto `atol_force` define cuándo se da por convergido un residuo numéricamente cero.
- Problemas casi-rígidos y plásticos perfectos se afinan independientemente vía `rtol_disp` vs `rtol_force` sin compromisos.
- `ArcLengthSolver` consume la misma función → cambiar la política llega automáticamente a ambos.

**Hacia adelante**:

- Solvers futuros (Newton-Krylov, Riks generalizado, predictor-corrector dinámico) heredan el patrón sin reescribir el criterio.
- Variantes del criterio (energético, residuo escalado por la diagonal de `K`, etc.) se introducen como implementaciones alternativas de `is_converged` en un dispatcher análogo al de `linalg`.
- Telemetría más rica: el ratio `residual/tol_force` por iteración da al usuario una métrica de "qué tan cerca del umbral" que no era visible con un único `err < tol`.

**Migración**:

- Ruptura limpia: tests existentes y configuraciones YAML que pasan `tol=` se actualizan en el mismo commit que introduce `ConvergenceCriterion`. No queda código legacy.
- Tests nuevos: invariancia bajo unidades (mismo problema en N/m, kN/mm, MN/m converge en el mismo número de iteraciones hasta paridad de bits gracias a `atol` autoderivado), comportamiento en régimen degenerado (control en desplazamiento con `F_ext = 0`), separación de criterios (problema casi-rígido donde solo `err_disp` baja debe converger; problema postcrítico donde solo `err_force` baja también).
- Manual de usuario y especificaciones YAML actualizadas con la nueva forma de pasar tolerancias.

## Alternativas consideradas

- **Mantener el esquema actual con un atol añadido al denominador** (`ref + atol`). Mejora marginal: protege del `0/0` pero sigue siendo conceptualmente "relativo con piso de referencia", no atol+rtol genuino. Rechazada por coherencia con ADR 0006.
- **Una sola tolerancia con criterio energético** (`δUᵀ·R ≤ tol·|F_ext·U_iter|`). Más elegante teóricamente, una sola escala. Descartada como default: oculta diagnósticos útiles (cuando δU es chico y R es grande, dice algo el operador; el energético los mezcla). Queda como variante futura introducible en `is_converged`.
- **Convergencia con OR semántica** (basta uno de los dos criterios). Más permisiva, estilo Code_Aster. Rechazada como default: en régimen postcrítico el desplazamiento puede ser microscópico aunque haya residuo grande, y viceversa; AND es la red de seguridad.
- **Tolerancias por DOF** (`atol` como vector de tamaño `n_free` con piso diferenciado para grados traslacionales vs rotacionales). Justificable físicamente (un giro de `1e-9 rad` y un desplazamiento de `1e-9 m` son cosas distintas), pero complica la API. Difere hasta que entre el primer problema donde la mezcla rotación-traslación cause un falso convergido visible.

## Deuda heredada

- **Sin enforcement** de que un solver no lineal futuro consuma `ConvergenceCriterion`. Si alguien escribe su propio criterio en un `RiksSolver` o variante, la convención no se puede imponer en construcción (idéntico a la deuda del ADR 0006).
- **`escala_disp_inicial` depende de `‖K_initial‖_∞`**, que requiere una norma matricial del primer ensamblaje. Para mallas grandes con `K` sparse esto es barato (norma del máximo absoluto por fila), pero queda como deuda evaluarlo cuidadosamente cuando lleguen problemas de varios millones de DOFs. Una aproximación más rápida si fuera necesario: tomar la diagonal máxima de `K_initial`. (2026-11-27: es lo que se hace, `stiffness_diag_scale`; el coste por iteración del piso de redondeo, medido, está en la enmienda.)

## Enmienda 2026-11-27 — piso de redondeo y calibración sin fuerzas

**Decisión del usuario** (Reglas §3: criterio de convergencia), tomada sobre la medición siguiente: un residuo del orden del error de redondeo cuenta como cero, en cualquier unidad. Verificación de `teorico-solidum` y revisión de `arquitecto-solidum` del 2026-11-27 (informes fuera del repositorio).

**Hallazgo** (al programar el contacto, ADR 0043). La escala de fuerza de la calibración era `max(‖F_ext‖, ‖F_int(U₀)‖, 1)`. Con control en desplazamiento desde el reposo las dos normas son nulas —el desplazamiento impuesto todavía no se ha aplicado en el primer ensamblaje— y el piso absoluto quedaba en `1e-9` *unidades de fuerza* del modelo. Medido:

- Un bloque que baja por una holgura de contacto (sin ninguna fuerza hasta tocar) se detiene en N con 30 y 100 pasos: su residuo de redondeo (1e-10 a 2e-9 N) no baja de 1e-9 N. Con el mismo problema en kN o en MN converge.
- Un voladizo con flecha impuesta acepta en MN un residuo 13 veces mayor que en N (tolerancia/fuerza en la punta: 2.2e-4 frente a 1.7e-5).
- Un voladizo de 3000 elementos se detiene en N: el residuo no baja del redondeo.
- Al continuar un estático con desplazamiento impuesto con un dinámico, la comprobación de equilibrio inicial avisaba de un desequilibrio inexistente (residuo relativo 1.0: el propio ruido frente a sí mismo).

**Decisión.** El criterio de fuerza queda

```
‖TᵀR‖ ≤ atol_force + rtol_force · max(‖F_ext‖, ‖F_int‖) + c_r · ε · ‖ |T|ᵀ · Σ |A|·|x| ‖
```

con los términos `(A, x)` del residuo de cada solver: `(K, U)` en Newton y en longitud de arco; `(K, u)`, `(C, u̇)` y `(M, ü)` en Newmark y HHT. Los valores absolutos se toman entrada a entrada, `T` es el operador de la reducción `u = T·u_libre + g` (ADR 0004), `ε` el épsilon de máquina y `c_r = roundoff_factor = 10` (`CONVERGENCE_ROUNDOFF_FACTOR`; clave YAML `convergence.roundoff_factor`). Calcular `A·x` en doble precisión deja en cada componente un error del orden de `ε·(|A||x|)_i` que ninguna iteración de Newton reduce (la resolución lineal sólo afecta a la velocidad de convergencia: el residuo se reevalúa en cada iteración); un residuo por debajo de ese piso no es un desequilibrio que el cálculo pueda resolver. Cada componente de la escala tiene las unidades de la componente del residuo —fuerza o momento—, así que el piso es tan invariante ante las unidades como el propio residuo, giros incluidos. En HHT-α los pesos de los términos del paso `n` y del `n+1` suman 1 para `α ∈ [−1/3, 0]`: basta la suma sin pesos en `n+1`. En la calibración se quita el `1`: sin fuerzas, `atol_force = 0` y la escala de desplazamiento es la norma de los valores impuestos `‖g‖` (`Assembler.imposed_values`). Las escalas nulas son válidas; con tolerancia nula (problema sin fuerzas ni desplazamientos) el paso converge si la magnitud es exactamente nula. La comprobación de equilibrio al continuar un estático con un dinámico descuenta el mismo piso, componente a componente, sin cambiar su umbral relativo (1e-6).

**Ruido medido** (`‖R‖` tras converger, iteraciones forzadas, frente a la escala; dos muestreos independientes, del autor y del teórico): 0.2–0.5 en voladizos y marcos con giros en m, mm y km (máximo 0.48), 0.07–0.13 en Hex8 elástico y plástico, hasta 0.31 en Newton-Newmark y Newton-HHT, 0.13–0.33 en contacto y hasta 0.45 en el bloque con holgura; en pasos sin fuerza cuyo residuo sale casi exacto, mucho menos. Con `ε·max|diag K|·‖U‖` (la primera versión de esta enmienda) el cociente bajaba a 1e-3 en m y a 1e-6 en km en el voladizo: la diagonal mayor mezcla rigideces de giro y de traslación, y el piso sobrestimaba el ruido hasta un millón de veces.

**Por qué `c_r = 10`.** Es empírico, no una cota rigurosa: el peor caso teórico es `c·m·ε·max|diag K|·‖U‖`, con `m` las entradas por fila (18 en Quad4, 81 en Hex8). Queda unas 20 veces por encima del mayor ruido medido (0.48). En los problemas cuyas fuerzas se resuelven bien es despreciable: en un voladizo de 10 elementos vale 3e-6 de la tolerancia, que es la misma relativa a la fuerza en N y en MN (1.7e-5).

**Coste del piso** (medido): un producto `|A|·|x|` por término y por iteración, sin copiar la matriz (núcleo de Numba sobre CSR): 29 ms con 31 millones de entradas, frente a 15 ms de leer la diagonal y 24 ms de un producto `K·x` de SciPy. Es del orden del 3 % de una iteración que reutiliza la factorización a 4e5 grados de libertad, y del 0.1 % de una que factoriza. En dinámica son hasta tres productos (`C` no tiene entradas sin amortiguamiento de Rayleigh).

**Contrato** (ruptura limpia, como en la decisión original; versión 0.1.0). El protocolo `NewtonProblem` (ADR 0015) gana el método obligatorio `roundoff_force(x, state)`, que el corrector pasa a `ConvergenceCriterion.evaluate(..., roundoff_force=)`, argumento de palabra clave también obligatorio, para que ningún solver vuelva en silencio a un piso que dependa de las unidades. `default_calibration_scales(F_ext, F_int, K, assembler)` recibe el ensamblador. Quién se rompe y cómo, todo ruidoso: un `NewtonProblem` externo sin el método → `TypeError` al entrar en `NewtonCorrector.run`, antes del primer ensamblaje, con los métodos que faltan; una subclase de `ConvergenceCriterion` con la firma anterior de `evaluate` → `TypeError`; quien llame a `default_calibration_scales` sin el ensamblador → `TypeError`. Las subclases de `ArcProblem` (módulos de usuario) heredan el método y no se rompen. `ConvergenceCriterion` valida que los cinco factores sean números finitos no negativos (el mensaje nombra la clave) y `ConvergenceState.floor_force` dice qué parte de la tolerancia aporta el piso.

**Alternativas descartadas.** (1) Escala de fuerza `‖K·g‖` (hipótesis inicial de la IA), la reacción de los desplazamientos impuestos con los nodos libres quietos: afloja 330 veces el criterio de un voladizo esbelto de 100 elementos (tolerancia/fuerza 5.7e-3 frente a 1.7e-5), porque la rigidez local del elemento de la punta es mucho mayor que la de la estructura. (2) `ε·max|diag K|·‖U‖`: barata, pero no invariante con giros (arriba). (3) `ε·‖diag(K)∘U‖`: el ruido llega a 10 veces esa escala en mm.

**Limitaciones.**

- *Límite de la doble precisión.* En un problema muy mal condicionado el piso supera a `rtol·fuerza` y se acepta el residuo en ese límite; el criterio de desplazamiento (AND) sigue exigiendo `rtol_disp`. En el voladizo de acero de 1 m y sección 10 × 10 mm con flecha impuesta (N-m): piso/(atol + rtol·ref) = 0.01 con 100 elementos, 0.47 con 300, 31 con 1000 y 1470 con 3000; con 3000 la tolerancia vale el 2.6 % de la fuerza en la punta y el residuo que queda, el 0.05 %. Es coherente con la cota `c_r·ε·max|diag K|/λ_min(K)` del error relativo en desplazamiento que el piso admite (2.9e-3 con 1000 elementos y 0.23 con 3000, *calculado*): con 1000 elementos o más la propia resolución lineal ya no da 1e-5.
- *Norma con unidades mezcladas.* El piso sigue al ruido en cualquier unidad, pero la tolerancia `rtol·max(‖F_ext‖, ‖F_int‖)` con la que se compara mezcla fuerzas y momentos en una sola norma (ver «Tolerancias por DOF» en las alternativas): en el voladizo, el cociente piso/tolerancia cambia 17 veces entre m y mm con 10 elementos y 206 veces con 100, porque en mm domina el momento del empotramiento (tolerancia/fuerza en la punta 1.0e-2 en mm frente a 1.7e-5 en m con 100 elementos).
- *Elementos geométricamente no lineales.* Los corotacionales (`Truss2DCorot`, `Frame2DEulerCorot`; previsiblemente `Truss3DCorot` y los cables) forman `x = X + u` y luego `x₂ − x₁`: su error de redondeo escala con las coordenadas `‖X‖`, no con `‖U‖`, y el piso no lo cubre. Medido por el teórico: ruido de 64 a 2.8e3 veces la escala en el origen (m, mm y km) y hasta 2.3e7 a 1e5 m de él; una armadura corotacional trasladada sin fuerzas se detiene (con este criterio y con el anterior: no es una regresión), y un marco con cargas reales se detiene a 1e5–1e6 m del origen (coordenadas UTM). Pendiente aparte: formar `(X₂ − X₁) + (u₂ − u₁)` y añadir al piso el término propio de estos elementos.
- Los umbrales de detección de divergencia del θ-método y de las diferencias centrales también tienen un `max(…, 1.0)`, pero no son criterios de convergencia y quedan fuera de esta enmienda.

**Pruebas.** `tests/test_convergence_criterion.py` (`TestRoundoffFloorUnit`, `TestForceFreeConvergence`): el bloque con holgura en N, kN y MN, la tolerancia del voladizo igual en N y en MN, el ruido frente a la escala en m, mm y km con giros, el voladizo de 3000 elementos, la etapa que parte de una traslación sin fuerzas, la continuación dinámica sin aviso falso, la fórmula del piso de cada problema de Newton, el rechazo de un problema incompleto y la validación de los factores. Las mutaciones están en `informes/prototipos_2026-11-27/convergencia/mutar_convergencia.py` (fuera del repositorio).

## Referencias

- ADR 0006 — Tolerancias del criterio de admisibilidad constitutiva (patrón base atol+rtol).
- Bathe, K.-J. (2014), *Finite Element Procedures*, §8.4.4 — criterios de convergencia en análisis incremental.
- Crisfield, M. A. (1991), *Non-linear Finite Element Analysis of Solids and Structures*, vol. 1, cap. 9.
- Abaqus Analysis User's Guide — *Convergence criteria for nonlinear problems* (`RTOL`, `CTOL`, `RTOLA`, `CTOLA`).
- ANSYS Mechanical APDL — `CNVTOL` (control de `MINREF` como piso absoluto).
- Code_Aster — `STAT_NON_LINE` / `CRITERE` (`RESI_GLOB_RELA`, `RESI_GLOB_MAXI`).
- SUNDIALS, *CVODES User Guide*, tolerancias mixtas atol+rtol.
