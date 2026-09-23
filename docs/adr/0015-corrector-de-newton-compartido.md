# ADR 0015 — Corrector de Newton compartido: un bucle, cinco solvers

- **Estado**: aceptado
- **Fecha**: 2026-09-22
- **Alcance**: `solidum/math/solvers/corrector.py` (nuevo) y los cinco solvers iterativos: `NonlinearSolver`, `ArcLengthSolver`, `DissipationArcLengthSolver`, `NewtonNewmarkSolver` y `NewtonHHTSolver`. **No cambia** ningún resultado, ninguna formulación, la API pública de los solvers (constructores, `solve`, metadatos, excepciones tipadas) ni el parser YAML. Retira ayudantes privados (`_armijo_step`, `_solve_reduced`, `_armijo_step_dynamic`, `_solve`, `_make_linalg`) y expone `solver.corrector`.

## Contexto

Los cinco solvers no lineales ejecutaban el mismo bucle de corrección escrito cinco veces: ensamblar en el iterado, formar el residuo, calibrar y evaluar el criterio dual (ADR 0007), resolver con el despachador algebraico y su degradación Cholesky → LU (ADR 0003), congelar la factorización en el Newton modificado, aplicar el line search de descenso no monótono (ADR 0011), acumular el historial para clasificar la divergencia (ADR 0011) y comitear al converger. Dos variantes del line search (`_armijo_step` estático y `_armijo_step_dynamic`), tres gestiones distintas del backend, dos formas de clasificar la divergencia.

La auditoría global de esta mañana lo hizo visible: la corrección "un ensamblaje por iteración, convergencia evaluada con `(‖R(x_k)‖, ‖δ_{k−1}‖)`" hubo que aplicarla en los cinco, y la coherencia del estado committed con el iterado guardado también. Cada solver nuevo (ADR 0009 añadió dos) volvía a copiar el bucle, y cada cambio de semántica exigía cinco diffs y cinco oportunidades de divergir. Lo que viene después —modelo de cargas con historia, acoplamiento termomecánico, un backend algebraico multihilo— tocaría de nuevo los cinco.

La restricción (`Reglas.md` §1, memoria `feedback_consistencia_arquitectural`): centralizar cuando hay dos casos reales; aquí había cinco.

## Decisión

### 1. Un bucle, un protocolo

`NewtonCorrector.run(problem, x0)` ejecuta el bucle; lo propio de cada solver entra por un objeto `NewtonProblem` que el solver construye **por paso**:

| Método | Qué aporta el solver |
|---|---|
| `assemble(x)` | el único ensamblaje de la iteración; devuelve el estado (`(K, F_int)`) |
| `residual(x, state)`, `residual_norm(R)` | forma del residuo y norma que entra al criterio (DOF libres o `‖TᵀR‖`) |
| `calibration_scales`, `reference_force`, `x_norm` | escalas del criterio dual (ADR 0007) |
| `correction(x, state, R, solve)` | reducción, sistema tangente y corrección `dx`, resolviendo con el `solve` del corrector |
| `apply(x, dx, α)` | iterado avanzado y norma del incremento |
| `on_converged(x, state)` | commit del estado trial del ensamblaje convergido |

El iterado `x` y el estado son opacos para el corrector: `U` en el Newton incremental, `(U, λ, ΔU)` en el arc-length, `(u, u̇, ü, ü_libre)` en los dinámicos. Un objeto por paso, y no una función, porque el problema lleva datos del paso (carga, predictores de Newmark, longitud de arco, términos `α·X_n` de HHT) y estado del paso (`ΔU` acumulado, `λ` corriente) que sólo tienen sentido dentro de él; el corrector conserva lo que persiste entre pasos: el backend y su degradación a LU.

### 2. Lo que el corrector posee

- **Backend algebraico** (ADR 0003): despachador por propiedades declaradas (`is_symmetric`, `is_positive_definite`), degradación permanente a LU si Cholesky reporta no-positividad, Newton modificado con factorización congelada a partir de la iteración `freeze_tangent_after_iter` y descartada al cerrar el paso.
- **Line search** (ADR 0011, variante GLL): ensambla en cada `α` probado y devuelve el estado del punto aceptado, así que nunca hay reensamblaje redundante; sin line search aplica `α = 1` y ensambla una vez.
- **Telemetría** (ADR 0011): `CorrectorResult` con historial de residuos e incrementos, bandera de tangente singular, último `α`; `divergence_error(...)` devuelve la excepción tipada lista para lanzar. `CorrectionAborted` permite al problema abandonar el paso (raíces imaginarias del arco, `α ≈ 0` en disipación) sin marcarlo como tangente singular.
- **Dos matices de arranque**: `check_initial` (el Newton incremental no evalúa la iteración 0, porque el iterado inicial aún no incorpora el incremento de Dirichlet del paso; arc-length y dinámicos sí evalúan el iterado predicho) e `initial_delta_norm` (norma del incremento predictor para el criterio en desplazamiento de la iteración 0).

### 3. Lo que queda en cada solver

El control de paso y su adaptatividad: bisección de `Δλ` y aceleración (Newton incremental); predictor tangente, sentido de avance, `dl`/`τ`, modos y switching (arc-length); predictores de Newmark, carga del instante, caché `F_int_n` de HHT (dinámicos). La restricción del arc-length es un método `constraint` de `_ArcProblem` que la variante por disipación sobreescribe para su modo: la subclase aporta once líneas de física y ninguna de bucle.

### 4. Equivalencia

La migración no cambia ninguna operación ni su orden: cada solver ensambla, resuelve y actualiza lo mismo que antes. La suite completa (1 447 tests, incluidos los de convergencia cuadrática, line search, Newton modificado, fallback SPD → LU, modos de divergencia y los de fidelidad del estado committed) pasa sin cambiar ningún valor esperado. Las únicas diferencias observables: el arc-length registra sus iteraciones con el mismo formato que el Newton incremental (`Iteración k | lam=… | R/tol_F …`), y con line search activo los dinámicos ahorran un ensamblaje redundante por iteración.

## Consecuencias

**Coste del solver N+1**: un solver iterativo nuevo escribe su `NewtonProblem` (residuo, sistema tangente, actualización) y su control de paso; hereda backend, Newton modificado, line search, calibración, criterio y telemetría sin copiar nada. Un cambio de semántica del bucle (otro criterio de parada, otra globalización, otro backend) se hace una vez.

**Superficie**: `solver.corrector` (público: `solve`, `line_search_step`, `frozen_factor`, `is_positive_definite`) sustituye a los ayudantes privados; los tests de internos se reescribieron contra él. `NonlinearSolver.make_problem(F_ext, λ)` construye el problema de un paso para pruebas.

**Hacia adelante**: el modelo de cargas con historia entra por `residual` y `calibration_scales` de un problema, no por cinco bucles; un residuo acoplado (u, T) es otro `NewtonProblem`; un backend nuevo se enchufa en `NewtonCorrector.solve`.

## Alternativas descartadas

- **Clase base abstracta con método plantilla** (los solvers heredan `solve` y sobreescriben ganchos): acopla el control de paso al bucle y obliga a una jerarquía entre solvers que no comparten control de paso (arc-length y Newmark no tienen nada en común fuera del corrector). El objeto-problema separa las dos cosas.
- **Corrector como función con *callbacks*** sueltos: seis funciones por llamada y sin sitio para el estado del paso; el protocolo agrupa lo mismo con nombre.
- **Mantener el duplicado** con la auditoría como red: es lo que había; cinco diffs por cambio no escalan al acoplamiento ni al modelo de cargas.

## Diferido

- Line search en el arc-length (hoy sólo lo usan el Newton incremental y los dinámicos; el corrector ya lo admite, faltaría justificarlo con un caso).
- Un gancho por iteración (`step_callback` existe sólo por paso) para instrumentación externa.
