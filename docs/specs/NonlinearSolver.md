# NonlinearSolver — Análisis estático no lineal incremental-iterativo (Newton-Raphson con paso adaptativo)

> Spec **retroactiva**: el solver existe desde la primera fase del proyecto y está validado por tests de plasticidad, daño, corotacional, embedded. Esta spec lo documenta sin cambiar comportamiento — H-5.3 de la auditoría 2026-05-18.

---

## Especificación física

### 0. Descripción general

Análisis estático **no lineal** con **control de carga** (load-controlled) en una estructura discretizada por FEM. La no-linealidad puede provenir de la geometría (formulación corotacional, grandes desplazamientos pequeñas deformaciones) o del material (plasticidad, daño, embedded discontinuity). Esquema **incremental-iterativo**: la carga objetivo se aplica por incrementos $\Delta\lambda$, y en cada incremento se resuelve el equilibrio por **Newton-Raphson** con tangente reensamblada y paso bisectado al fallar.

### 1. Ecuación de equilibrio resuelta

Forma fuerte semidiscreta, residuo:

$$\mathbf R(\mathbf U, \lambda) \;=\; \lambda\,\mathbf F_{\text{ext}}^{\text{ref}} - \mathbf F_{\text{int}}(\mathbf U) \;=\; \mathbf 0$$

con $\lambda \in [0, 1]$ factor de carga aplicado a la carga de referencia $\mathbf F_{\text{ext}}^{\text{ref}}$, y $\mathbf F_{\text{int}}(\mathbf U)$ fuerzas internas no lineales (dependientes del estado actual de los puntos de Gauss).

### 2. Condiciones de contorno

Idénticas a `LinearSolver`: Dirichlet (homogéneo / no homogéneo) y MPC vía ADR 0004. Las cargas Dirichlet no homogéneas se aplican proporcionalmente al factor de carga (incremento $\Delta\mathbf g = \mathbf g \cdot \Delta\lambda$).

### 3. Salidas físicas

- $\mathbf U$ final converged en $\lambda = 1$.
- Historia disponible vía `step_callback(step, U, load_factor)` que el usuario puede pasar al solver.
- Diagnóstico de divergencia con excepciones tipadas (ADR 0011): `OscillatingNewtonError`, `SingularTangentError`, `LoadExceedsCapacityError`, `UnknownDivergenceError`, todas subclases de `SolverDivergedError`.

---

## Formulación numérica

### 4. Esquema operativo

```
λ = 0
Δλ = 1 / num_steps
mientras λ < 1:
    intentar avanzar λ → λ + Δλ:
        U_iter = U_current
        prepare_all_steps(U_current)              # hook ADR 0010 (embedded)
        para iter en range(max_iter):
            K_t, F_int = assemble_non_linear_system(U_iter)
            R = (λ+Δλ)·F_ext - F_int
            δU = K_t^{-1} · R                     # reducido ADR 0004
            U_iter, ‖R‖, F_int = line_search(...) # ADR 0011, opt-in
            si converge (ADR 0007 dual):
                commit_all_states()
                λ += Δλ
                si iter < threshold y Δλ < base:
                    Δλ ← min(Δλ · growth, base)   # adaptativo: acelerar
                break
        si no converge:
            Δλ /= 2                               # bisección
            si Δλ < min_delta_lambda:
                clasificar divergencia y raise    # ADR 0011
```

### 5. Predictor / corrector

**Predictor**: $\Delta\mathbf U^{(0)} = \mathbf 0$ — el incremento empieza desde el estado convergido del paso anterior.

**Corrector** (Newton iter $k$): $\mathbf K_t^{(k)}\,\delta\mathbf U^{(k)} = \mathbf R^{(k)}$, $\mathbf U^{(k+1)} = \mathbf U^{(k)} + \alpha\,\delta\mathbf U^{(k)}$, con $\alpha$ del line search (default $\alpha = 1$).

### 6. Criterio de convergencia (ADR 0007)

Patrón dual fuerza + desplazamiento, con escalas autoderivadas:

$$\lVert\mathbf R[\text{libres}]\rVert \le \text{atol}_F + \text{rtol}_F \cdot \max(\lVert\mathbf F_{\text{ext}}\rVert, \lVert\mathbf F_{\text{int}}\rVert)$$

$$\lVert\Delta\mathbf U\rVert \le \text{atol}_d + \text{rtol}_d \cdot \lVert\mathbf U\rVert$$

Semántica **AND** (ambos simultáneamente). Calibración: en el primer ensamblaje del solve, `force_scale = max(‖F_ext‖, ‖F_int‖, 1)` y `disp_scale = force_scale / max(|diag(K)|)`.

### 7. Imposición de Dirichlet y MPC

Idéntica a `LinearSolver` — ADR 0004. La carga Dirichlet no homogénea $\mathbf g$ se aplica proporcionalmente al factor de carga $\lambda$ a través del operador del Assembler (`U_current=U_iter, load_factor=next_load_factor` en la llamada a `reduce`).

### 8. Backend algebraico

Despacho automático ADR 0003 idéntico a `LinearSolver`. Fallback Cholesky→LU activado.

**Newton modificado (ADR 0003 fase 2)**: parámetro opcional `freeze_tangent_after_iter`. Si es `int N`, el solver factoriza fresca durante las primeras $N$ iteraciones del paso y reutiliza la factorización en las siguientes iteraciones del mismo paso (no entre pasos — la tangente cambia con $\mathbf U$). Default `None` = Newton estándar (factoriza cada iter).

### 9. Adaptatividad / control de paso

**Acelerado**: si un paso converge en $< $ `NEWTON_ADAPTIVE_GROWTH_ITER_THRESHOLD` iteraciones y $\Delta\lambda < 1/$ `num_steps`, $\Delta\lambda$ se multiplica por `NEWTON_ADAPTIVE_GROWTH_FACTOR` (cota: $1/$`num_steps`).

**Bisección**: si un paso no converge en `max_iter`, $\Delta\lambda \to \Delta\lambda / 2$ y se reintenta. Aborto si $\Delta\lambda < $ `min_delta_lambda`.

**Line search ADR 0011** (opt-in, default `False`): variante GLL (Grippo-Lampariello-Lucidi) — backtracking simple por descenso no monótono. Justificación: line search Wolfe puro rechaza pasos Newton correctos en regímenes plásticos/dañados con residuo transitorio creciente. Documentado en ADR 0011.

### 10. Caveats numéricos

- **Snap-through / snap-back**: NO soportado — `NonlinearSolver` con control de carga falla en puntos límite con $dU/d\lambda \to \infty$ o $d\lambda/dU < 0$. Usar `ArcLengthSolver`. Excepción: el solver atraviesa puntos límite *suaves* (sin cambio de signo en $d\lambda$) con bisección + cinemática corotacional, hallazgo empírico de la auditoría 2026-05-18.
- **Softening severo con penalty cohesivo stiff** (embedded): el solver no atraviesa la transición elástico→softening con $K_t$ casi singular en $\kappa_0$. Limitación específica al embedded; el daño bulk continuo (1D/2D) sí se atraviesa.
- **Plasticidad perfecta con line search off**: puede oscilar (residuo periódico, no diverge ni converge). Mitigación: activar `line_search=True`.
- **Tangente singular**: lanza `SingularTangentError` con el último estado.
- **Newton modificado**: con `freeze_tangent_after_iter` no nulo, la convergencia cuadrática se degrada a lineal en las iteraciones congeladas. Útil cuando el coste de factorizar domina (`K_t` grande, materiales con tangente cara de evaluar).

---

## Contrato de implementación

```yaml
name: NonlinearSolver
kind: solver
status: validated

interface:
  yaml_type: nonlinear
  output: SolveResult
  pipeline_kind: static

parameters:
  - { name: convergence,                 type: ConvergenceCriterion, required: false, default: "ConvergenceCriterion()",
      desc: "Política dual fuerza+desplazamiento (ADR 0007)" }
  - { name: max_iter,                    type: int,   required: false, default: 20,
      desc: "Iteraciones Newton máximas por paso antes de bisectar (keyword-only)" }
  - { name: num_steps,                   type: int,   required: false, default: 10,
      desc: "Número base de incrementos para alcanzar λ=1 (keyword-only)" }
  - { name: adaptive,                    type: bool,  required: false, default: true,
      desc: "Habilita acelerado + bisección de Δλ (keyword-only)" }
  - { name: min_delta_lambda,            type: float, required: false, default: NEWTON_DEFAULT_MIN_DELTA_LAMBDA,
      desc: "Cota inferior para Δλ; por debajo se aborta con SolverDivergedError (keyword-only)" }
  - { name: linear_algebra,              type: str,   required: false, default: "auto",
      desc: "'auto' (Cholesky→LU según simetría), 'cholesky', 'lu' (keyword-only)" }
  - { name: freeze_tangent_after_iter,   type: int|None, required: false, default: null,
      desc: "Newton modificado: factoriza fresca primeras N iter, congela después (keyword-only)" }
  - { name: line_search,                 type: bool,  required: false, default: false,
      desc: "Line search GLL ADR 0011 (keyword-only)" }

requirements:
  - "Modelo con no-linealidad bien definida (material con tangente válida en todo el rango de carga)"
  - "Cargas Dirichlet no homogéneas constantes en magnitud (se aplican proporcionalmente a λ)"
  - "Apoyos suficientes para descartar modos rígidos"

conventions:
  units: "heredadas del modelo (ADR 0008)"
  stability: "no garantizada — los puntos límite con d λ/dU<0 abortan vía bisección"

out_of_scope:
  - "Snap-through / snap-back con d λ/dU<0 ⇒ usar ArcLengthSolver"
  - "Análisis dinámico ⇒ usar NewtonNewmarkSolver, NewtonHHTSolver, etc."
  - "Buckling lineal (eigenproblema generalizado) ⇒ futuro BucklingSolver"

acceptance:
  verification:
    - name: recuperacion_caso_lineal
      setup: "modelo enteramente lineal (Elastic2D + Quad4 + sin corotacional)"
      expect: "resultado idéntico a LinearSolver; Newton converge en 1 iteración por paso"
      tol_rel: 1.0e-12

    - name: plasticidad_uniaxial_truss
      setup: "Truss2D con Elastoplastic1D, carga axial creciente hasta plástico"
      expect: "respuesta bilineal exacta (E hasta yield, E_t = E·H/(E+H) después)"
      tol_rel: 1.0e-6

    - name: corotacional_voladizo_grande
      setup: "viga Frame2DEulerCorot, carga lateral grande (rotación >30°)"
      expect: "δ obtenido coincide con benchmark Bathe-Bolourchi 1979"
      tol_rel: 1.0e-3

    - name: damage_bulk_softening
      setup: "Quad4 con IsotropicDamage2D, tracción monotónica hasta post-pico"
      expect: "rama softening trazada; convergencia paso a paso con bisección"
      tol_rel: 1.0e-3

    - name: divergencia_clasificada
      setup: "carga que excede la capacidad última (plasticidad perfecta con mecanismo)"
      expect: "raise OscillatingNewtonError o LoadExceedsCapacityError con métricas estructuradas (ADR 0011)"

    - name: paso_adaptativo_acelera_y_bisecta
      setup: "trayectoria con tramos elásticos rápidos y zona plástica lenta"
      expect: "Δλ crece en tramos fáciles y bisecta en zonas duras; alcanza λ=1"

  specific:
    - name: newton_modificado_factoriza_menos
      setup: "freeze_tangent_after_iter=2, 5 iter típicas por paso"
      expect: "factorizaciones por paso ≤ 3 en lugar de 5; resultado dentro de tolerancia"
      tol_rel: 1.0e-8

    - name: line_search_corrige_oscilacion_DP
      setup: "Drucker-Prager 2D perfectamente plástico con carga que provoca oscilación"
      expect: "line_search=True converge; line_search=False oscila o requiere bisección extra"

references:
  - "Crisfield M.A. (1991). Non-linear Finite Element Analysis of Solids and Structures, Vol. 1. Wiley. §1-§9."
  - "Bonet J., Wood R.D. (2008). Nonlinear Continuum Mechanics for Finite Element Analysis. Cambridge. §8-§9."
  - "de Borst R., Sluys L.J. (1999). Computational Methods in Non-linear Solid Mechanics. TU Delft. §3 (Newton-Raphson y variantes)."
  - "ADR 0003, ADR 0004, ADR 0006, ADR 0007, ADR 0011."
```

---

## Implementación

- **Archivo**: [solidum/math/solvers/nonlinear.py](../../solidum/math/solvers/nonlinear.py).
- **Clase**: `NonlinearSolver`, registrada vía `@SolverRegistry.register` con `PIPELINE_KIND = "static"`.
- **Argumentos keyword-only**: todos menos `assembler` y `convergence` son keyword-only (auditoría H-1.5, sesión 2026-05-19, commit `4e4ed54`).
- **Bisección y clasificación**: en `solve`, al fallar `Δλ < min_delta_lambda` se llama a `classify_divergence(residual_history, delta_history, singular_tangent_detected)` que devuelve la subclase apropiada de `SolverDivergedError` con métricas tipadas.
- **Entrypoint público**: `solidum.run_static(model, solver="nonlinear", ...)`.
- **Tests**: cobertura amplísima — todos los tests de plasticidad/daño/corotacional/embedded del pipeline estático lo ejercitan. Tests específicos del comportamiento adaptativo y de divergencia tipada en `tests/test_solver_robustness.py` y `tests/test_solver_diagnostics.py`.

---

## Diálogo

- **2026-05-19** · Spec creada retroactivamente para cerrar el hueco H-5.3. Solver anterior a la convención de specs. Comportamiento documenta el estado tras la sesión de saneamiento post-auditoría (commits `4e4ed54` para keyword-only, ADR 0011 ya integrado).
