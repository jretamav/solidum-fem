# ADR 0007 — Tolerancias de convergencia en solvers no lineales

- **Estado**: aceptado
- **Fecha**: 2026-05-13
- **Alcance**: `NonlinearSolver` y `ArcLengthSolver` en `fenix/math/solvers.py`; constantes en `fenix/constants.py`; tests de convergencia en `tests/`.

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

con las escalas evaluadas en el estado corriente de la iteración:

- `ref_force = max(‖F_ext_step‖, ‖F_int‖)`. Sin `ZERO_TOL` añadido — `atol_force` ya es el piso absoluto.
- `ref_disp  = ‖U_iter‖`. Sin `ZERO_TOL` añadido — `atol_disp` ya es el piso absoluto.

Las normas son L2 sobre DOFs libres (igual que hoy: prescritos contienen reacciones, no fallos de equilibrio).

### Constantes nuevas en `fenix/constants.py`

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

Se introduce un módulo nuevo `fenix/math/convergence.py` con una clase pequeña que encapsula configuración + estado calibrado + evaluación:

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
- **`escala_disp_inicial` depende de `‖K_initial‖_∞`**, que requiere una norma matricial del primer ensamblaje. Para mallas grandes con `K` sparse esto es barato (norma del máximo absoluto por fila), pero queda como deuda evaluarlo cuidadosamente cuando lleguen problemas de varios millones de DOFs. Una aproximación más rápida si fuera necesario: tomar la diagonal máxima de `K_initial`.

## Referencias

- ADR 0006 — Tolerancias del criterio de admisibilidad constitutiva (patrón base atol+rtol).
- Bathe, K.-J. (2014), *Finite Element Procedures*, §8.4.4 — criterios de convergencia en análisis incremental.
- Crisfield, M. A. (1991), *Non-linear Finite Element Analysis of Solids and Structures*, vol. 1, cap. 9.
- Abaqus Analysis User's Guide — *Convergence criteria for nonlinear problems* (`RTOL`, `CTOL`, `RTOLA`, `CTOLA`).
- ANSYS Mechanical APDL — `CNVTOL` (control de `MINREF` como piso absoluto).
- Code_Aster — `STAT_NON_LINE` / `CRITERE` (`RESI_GLOB_RELA`, `RESI_GLOB_MAXI`).
- SUNDIALS, *CVODES User Guide*, tolerancias mixtas atol+rtol.
