# ArcLengthSolver — Análisis estático no lineal por longitud de arco cilíndrico (Crisfield)

> Spec **retroactiva**: el solver existe desde la fase de daño/softening y está validado por tests de snap-through, snap-back y daño post-pico. Esta spec lo documenta sin cambiar comportamiento — H-5.3 de la auditoría 2026-05-18.

---

## Especificación física

### 0. Descripción general

Análisis estático **no lineal** que traza curvas de equilibrio $\mathbf U(\lambda)$ atravesando puntos límite (snap-through, snap-back) donde $d\lambda/dU$ cambia de signo. A diferencia de `NonlinearSolver` (control de carga, $\lambda$ creciente impuesto), el arc-length **trata $\lambda$ como incógnita adicional** y restringe el avance combinado en $(U, \lambda)$ a un arco fijo $dl$ en el espacio aumentado. Esto permite seguir la respuesta postcrítica de estructuras con softening, pandeo geométrico, y daño con rama descendente.

Algoritmo: **Crisfield cilíndrico** (1981) — restricción cuadrática $\lVert\Delta\mathbf U\rVert^2 = dl^2$, $\lambda$ libre.

### 1. Ecuación de equilibrio resuelta

Forma fuerte:

$$\mathbf F_{\text{int}}(\mathbf U) \;=\; \lambda\,\mathbf F_{\text{ext}}^{\text{ref}}$$

con $\lambda$ incógnita, sujeta a la restricción cilíndrica de Crisfield en cada paso:

$$\lVert\Delta\mathbf U^{\text{paso}}\rVert^2 \;=\; dl^2$$

(versión **cilíndrica**: no incluye $\Delta\lambda^2$, equivalente a tomar $\psi = 0$ en el parámetro de escala carga-desplazamiento; la versión esférica $\lVert\Delta\mathbf U\rVert^2 + \psi^2\Delta\lambda^2\lVert\mathbf F_{\text{ext}}^{\text{ref}}\rVert^2 = dl^2$ no está implementada).

### 2. Condiciones de contorno

Idénticas a `NonlinearSolver`: Dirichlet (homogéneo / no homogéneo) y MPC vía ADR 0004. Las condiciones Dirichlet no homogéneas se aplican proporcionalmente a $\lambda$.

### 3. Salidas físicas

- $\mathbf U$ final en $\lambda \to \lambda_{\max}$ (o cuando se alcance `max_steps`).
- Historia disponible vía `step_callback(step, U, lambda)`.
- Trazas $\lambda(U_i)$ para visualizar snap-through/snap-back.

---

## Formulación numérica

### 4. Esquema operativo

```
λ = 0
dl = (sin fijar)
mientras λ < max_lambda y step < max_steps:
    # Predictor tangente
    K_t, F_int = assemble_non_linear_system(U_iter)
    δu_t = K_t^{-1} · F_ext^ref          # desplazamiento tangente unitario
    si dl sin fijar:                     # primer paso (§5)
        dl = Δl₁ = initial_dl  ó  initial_dlambda · ‖δu_t‖
    sign = sign(δU_prev · δu_t)          # evitar revertir dirección
    Δλ_pred = sign · dl / ‖δu_t‖
    ΔU_pred = Δλ_pred · δu_t

    # Corrector iterativo
    para iter en range(max_iter):
        K_t, F_int = assemble_non_linear_system(U_iter)
        R = (λ + Δλ)·F_ext^ref - F_int
        δu_R = K_t^{-1} · R               # corrección por residuo
        δu_t = K_t^{-1} · F_ext^ref       # tangente (refresh)

        # Restricción cuadrática de Crisfield
        # ‖ΔU_new + ddλ·δu_t‖² = dl²
        # ⇒ resolver ecuación cuadrática en ddλ, elegir raíz por menor ángulo
        ddλ = raiz_de_menor_angulo(...)
        ΔU += δu_R + ddλ · δu_t
        Δλ += ddλ
        U_iter = U_current + ΔU
        si converge (ADR 0007):
            si λ + Δλ > max_lambda:              # cruzó el objetivo (§10)
                no consolidar; llegada exacta:
                Newton puro con λ = max_lambda desde (U_current, λ),
                arrancando de la interpolación lineal del tramo recorrido
            commit_all_states()
            λ += Δλ
            ajustar_dl_adaptativamente(iter)
            break
    si no converge (o falla la llegada):
        dl /= 2
```

### 5. Predictor / corrector

**Predictor tangente**: $\delta\mathbf u_t = \mathbf K_t^{-1}\,\mathbf F_{\text{ext}}^{\text{ref}}$ (vector que mide la pendiente local del camino de equilibrio). $\Delta\lambda_{\text{pred}} = \pm dl / \lVert\delta\mathbf u_t\rVert$ con signo determinado por el producto escalar contra el incremento previo $\Delta\mathbf U^{\text{paso prev}}$ — así no se retrocede en el camino tras pasar un punto límite.

**Primer paso (2026-09-23).** $dl$ es una longitud en unidades de desplazamiento, y quien plantea el modelo no conoce de antemano sus desplazamientos. El primer paso se declara, por omisión, como **fracción de la carga de referencia** $\Delta\lambda_1$ (`initial_dlambda`, 0.1 si no se declara) y la longitud de arco se deriva del predictor elástico del primer paso (Crisfield 1991, cap. 9):

$$\Delta l_1 = \Delta\lambda_1\,\lVert\mathbf K_0^{-1}\,\mathbf F_{\text{ext}}^{\text{ref}}\rVert .$$

En régimen elástico el primer paso lleva exactamente $\Delta\lambda_1$, en cualquier sistema de unidades. `initial_dl` fija $\Delta l_1$ directamente (longitud explícita) y es excluyente con `initial_dlambda`. $\Delta l_1$ queda en `dl_reference` y es la referencia de `dl_max_factor` y del umbral de aborto. Una carga de referencia que no produce desplazamiento ($\lVert\delta\mathbf u_t\rVert = 0$) no da escala al primer paso y es un error explícito. La norma mezcla todos los grados de libertad; en marcos, giros y desplazamientos, cuyo peso relativo depende de las unidades (deuda #22): el primer paso es invariante, los siguientes heredan esa mezcla.

**Corrector** (Newton modificado para arc-length, Crisfield): cada iteración resuelve dos sistemas (residuo y tangente), combina con $dd\lambda$ que satisface la restricción cuadrática:

$$\lVert\Delta\mathbf U + dd\lambda\,\delta\mathbf u_t\rVert^2 = dl^2$$

Ecuación cuadrática en $dd\lambda$ con dos raíces. **Selección de raíz**: la que produce **menor ángulo** con el incremento acumulado $\Delta\mathbf U$ (evita pivotear el sentido del arco).

### 6. Criterio de convergencia (ADR 0007)

Idéntico a `NonlinearSolver` — patrón dual fuerza + desplazamiento con escala autoderivada. Compartido vía `ConvergenceCriterion`. La calibración del solver usa $\lVert\mathbf F_{\text{ext}}^{\text{ref}}\rVert$ como referencia de fuerza (no $\lambda\,\mathbf F_{\text{ext}}^{\text{ref}}$, porque $\lambda$ varía).

### 7. Imposición de Dirichlet y MPC

Idéntico al `NonlinearSolver` — `Assembler.reduce(...)` con `U_current` y `load_factor=lambda_iter`.

### 8. Backend algebraico

Régimen **postcrítico**: $\mathbf K_t$ puede ser indefinida (autovalores negativos tras bifurcación). El solver fuerza `is_positive_definite=False` al construir el `linalg`, así el despachador ADR 0003 elige LU automáticamente. Si el usuario fuerza `linear_algebra="cholesky"` desde YAML y $\mathbf K_t$ se vuelve indefinida en algún paso, hay fallback automático a LU con warning.

### 9. Adaptatividad — ajuste de $dl$

Política de auto-ajuste basada en iteraciones del paso:

- Converge en $<$ `dl_grow_iter_threshold` iteraciones ⇒ $dl \to dl \cdot $ `dl_grow_factor` (cota: $\Delta l_1 \cdot$ `dl_max_factor`).
- Converge en $>$ `dl_shrink_iter_threshold` iteraciones ⇒ $dl \to dl \cdot $ `dl_shrink_factor`.
- No converge en `max_iter` ⇒ $dl /= 2$ y reintentar el paso.
- Cota inferior: $dl < $ `ARCLENGTH_MIN_DL_FACTOR` $\cdot\,\Delta l_1$ ⇒ aborto.

### 10. Caveats numéricos

- **Cancelación o bifurcación verdadera**: en bifurcaciones simétricas múltiples, la selección de raíz por menor ángulo puede saltar a un camino paralelo. Para análisis de bifurcación auténtica usar perturbación + post-procesado modal.
- **Raíces imaginarias**: si el discriminante de la cuadrática es negativo (paso demasiado grande, mal condicionamiento severo), el solver aborta el paso y bisecta $dl$.
- **Llegada a $\lambda_{\max}$ (2026-09-23)**: todos los pasos llevan la restricción cilíndrica. Si uno converge con $\lambda > \lambda_{\max}$, su estado no se consolida y el paso se repite desde el último estado convergido con $\lambda = \lambda_{\max}$ fijo (Newton puro), arrancando de la interpolación lineal del tramo recién recorrido; como ese tramo cruza $\lambda_{\max}$, el equilibrio existe y está cerca. Si la llegada no converge, se biseca $dl$. Sustituye al "último paso" anterior, que se decidía extrapolando la tangente ($\lambda + \Delta\lambda_{\text{pred}} \ge \lambda_{\max}$) e imponía $\lambda_{\max}$ en control de carga sin recorrer la curva: podía saltar un punto límite, y un paso de arco que convergía por encima de $\lambda_{\max}$ terminaba el trazado ahí (medido: arco de von Mises, $\lambda_{\text{final}} = 0.398$ con $\lambda_{\max} = 0.39$).
- **Softening con penalty cohesivo stiff** (CST_Embedded2D): no se atraviesa la transición elástico→softening por la quasi-singularidad de $\mathbf K_t$ cerca de $\kappa_0$. Limitación específica al embedded.
- **Diagnóstico de bifurcación por inercia** (Sturm sequence): placeholder en `_negative_pivots()` — retorna `None` hasta que se implemente un backend LDLᵀ verdadero (Bunch-Kaufman) con conteo de pivots negativos.

---

## Contrato de implementación

```yaml
name: ArcLengthSolver
kind: solver
status: validated

interface:
  yaml_type: arclength
  output: SolveResult
  pipeline_kind: static

parameters:
  - { name: convergence,              type: ConvergenceCriterion, required: false, default: "ConvergenceCriterion()",
      desc: "Política dual fuerza+desplazamiento (ADR 0007), compartida con NonlinearSolver" }
  - { name: max_iter,                 type: int,   required: false, default: 20,
      desc: "Iteraciones Newton (corrector) máximas por paso" }
  - { name: max_lambda,               type: float, required: false, default: 1.0,
      desc: "Factor de carga objetivo" }
  - { name: initial_dlambda,          type: float, required: false, default: 0.1,
      desc: "Fracción de la carga de referencia del primer paso (adimensional, keyword-only); Δl₁ = Δλ₁·‖K₀⁻¹F_ref‖. Excluyente con initial_dl" }
  - { name: initial_dl,               type: float, required: false, default: null,
      desc: "Longitud de arco del primer paso en unidades de desplazamiento; sólo si se conoce la escala del problema. Excluyente con initial_dlambda" }
  - { name: max_steps,                type: int,   required: false, default: 100,
      desc: "Cota dura sobre el número de pasos" }
  - { name: dl_grow_factor,           type: float, required: false, default: 1.5 }
  - { name: dl_max_factor,            type: float, required: false, default: 5.0,
      desc: "dl máximo = Δl₁ · dl_max_factor" }
  - { name: dl_shrink_factor,         type: float, required: false, default: 0.6 }
  - { name: dl_grow_iter_threshold,   type: int,   required: false, default: 4 }
  - { name: dl_shrink_iter_threshold, type: int,   required: false, default: 8 }
  - { name: linear_algebra,           type: str,   required: false, default: "auto",
      desc: "'auto' (Pardiso si está instalado, si no LU: el solver declara la tangente potencialmente indefinida y no usa Cholesky), 'cholesky', 'pardiso', 'lu', 'iterative[:amg|jacobi|none]' (CG que pasa a MINRES al detectar curvatura negativa)" }

requirements:
  - "Modelo con potencial snap-through/snap-back (geometría inestable, softening del material)"
  - "Apoyos suficientes; sin modos rígidos"
  - "Carga de referencia que produzca desplazamiento (escala del primer paso)"

conventions:
  units: "heredadas del modelo (ADR 0008)"
  stability: "incondicional respecto a puntos límite; bifurcaciones complejas pueden saltar de rama"

out_of_scope:
  - "Variante esférica del arc-length (ψ ≠ 0) — sólo cilíndrica implementada"
  - "Análisis dinámico ⇒ usar NewmarkSolver y derivados"
  - "Continuación de bifurcaciones múltiples — requiere perturbación explícita por el usuario"
  - "Buckling lineal (eigenproblema) ⇒ futuro BucklingSolver"
  - "Sign-of-pivot tracking (Sturm sequence) — placeholder hasta LDLᵀ verdadero"

acceptance:
  verification:
    - name: barra_traccion_elastica_reproduce_caso_lineal
      setup: "Truss2D elástica, carga axial; arc-length avanza hasta λ=1"
      expect: "U(λ=1) = U_lineal exacto; equilibrio en cada paso"
      tol_rel: 1.0e-10

    - name: snap_through_lee_frame
      setup: "frame Lee con dos barras corotacionales en V, carga vertical en el vértice"
      expect: "se atraviesa el punto límite; respuesta U-λ coincide con Crisfield 1991 Vol.1 Fig. 9.5"
      tol_rel: 1.0e-3

    - name: snap_back_truss_inestable
      setup: "armadura con elemento corotacional cuya rama descendente requiere d λ/dU < 0"
      expect: "se atraviesa el snap-back; signo de Δλ cambia coherentemente"
      tol_rel: 1.0e-3

    - name: damage_softening_post_pico
      setup: "Quad4 con IsotropicDamage2D, tracción monotónica que entra a rama descendente"
      expect: "rama post-pico trazada; convergencia paso a paso"
      tol_rel: 1.0e-3

    - name: ajuste_dl_adaptativo
      setup: "trayectoria con tramos elásticos rápidos y zona post-pico lenta"
      expect: "dl crece en tramos fáciles, se reduce cerca del pico; converge globalmente"

    - name: ultimo_paso_cierra_en_target
      setup: "max_lambda = 0.7 sobre estructura elástica"
      expect: "λ_final = 0.7 exacto (no por overshoot)"
      tol_abs: 1.0e-12

    - name: primer_paso_fraccion_de_carga
      setup: "placa elástica; sin declarar el primer paso, y con initial_dlambda = 0.25"
      expect: "λ del primer paso = 0.1 y 0.25 respectivamente"
      tol_abs: 1.0e-12

    - name: trazado_invariante_ante_unidades
      setup: "misma placa J2 en m y en mm (longitudes ×1000, módulos /10⁶), primer paso por omisión"
      expect: "misma secuencia de λ, U_mm = 1000·U_m, dl_reference ×1000; con una longitud explícita igual en ambos, primeros λ distintos"
      tol_rel: 1.0e-8

    - name: llegada_exacta_dentro_del_tramo
      setup: "arco de von Mises h/L = 0.1, max_lambda = 0.39 (sobre el pico, 0.385), initial_dl = 0.2 y 1.0"
      expect: "λ_final = 0.39 exacto y equilibrio en ese estado (antes: 0.398 y 0.440)"
      tol_abs: 1.0e-12

  specific:
    - name: seleccion_raiz_no_revierte
      setup: "paso tras un snap-back con dos raíces de signos opuestos"
      expect: "se elige la raíz que continúa el camino (menor ángulo con incremento previo)"

    - name: fallback_cholesky_lu_indefinitud
      setup: "linear_algebra='cholesky' forzado, modelo con punto límite"
      expect: "warning + LU al detectar no-PD; análisis continúa"

references:
  - "Crisfield M.A. (1981). A fast incremental/iterative solution procedure that handles snap-through. Computers and Structures 13, 55-62."
  - "Crisfield M.A. (1991). Non-linear Finite Element Analysis of Solids and Structures, Vol. 1. Wiley. §9 (arc-length cilíndrico)."
  - "Riks E. (1979). An incremental approach to the solution of snapping and buckling problems. Int. J. Solids and Structures 15, 529-551."
  - "de Borst R., Sluys L.J. (1999). Computational Methods in Non-linear Solid Mechanics. TU Delft. §3.5."
  - "ADR 0003, ADR 0007."
```

---

## Implementación

- **Archivo**: [solidum/math/solvers/arclength.py](../../solidum/math/solvers/arclength.py).
- **Clase**: `ArcLengthSolver`, registrada vía `@SolverRegistry.register` con `PIPELINE_KIND = "static"`.
- **Restricción cuadrática y selección de raíz**: implementada en `solve`. Comparación de ángulos `theta1 = ΔU_iter · (ΔU_new + ddl1·δu_t)`, idem `theta2`, seleccionar el mayor → menor ángulo con la dirección previa.
- **`_make_linalg`**: fuerza `is_positive_definite=False` independientemente de la simetría del dominio (régimen postcrítico).
- **Entrypoint público**: `solidum.run(domain, solver=ArcLengthSolver(assembler, ...), F_applied=...)` o `solidum.run_yaml(path)`.
- **Tests**:
  - [tests/test_solver_robustness.py](../../tests/test_solver_robustness.py) (`TestArcLengthRobustness`: punto límite, snap-back, recuperación lineal), [tests/test_snap_through_corot.py](../../tests/test_snap_through_corot.py) (armadura de von Mises contra la solución cerrada) y [tests/test_integration.py](../../tests/test_integration.py) (`test_arclength_solver_elastoplastic`).
  - [tests/test_solver_robustness.py](../../tests/test_solver_robustness.py) · `test_arc_length_traverses_damage_softening`.
  - [tests/test_arclength_paso_inicial.py](../../tests/test_arclength_paso_inicial.py) — primer paso como fracción de carga, invariancia ante las unidades, declaración inválida, carga de referencia nula, clave YAML y llegada exacta a `max_lambda`.

---

## Diálogo

- **2026-05-19** · Spec creada retroactivamente para cerrar el hueco H-5.3. Solver anterior a la convención de specs. La spec recoge la implementación tal como está al cierre de la sesión de saneamiento post-auditoría.
- **2026-09-22** · Auditoría global: corrector con un ensamblaje por iteración y convergencia evaluada antes de resolver: el estado committed corresponde exactamente a `U_current` (antes se evaluaba R en U_k, se comiteaba el trial de U_k y se guardaba U_{k+1}: un iterado de desfase, del orden de la tolerancia). La parada por `max_steps` sin alcanzar `max_lambda` emite WARNING y se expone en `reached_max_lambda`, `lambda_final` y `steps_done`; `solidum.run` escala `F_applied` por λ_final (reacciones coherentes) y reporta `converged`/`num_steps` reales (antes siempre `True`/1).
- **2026-09-22** · ADR 0015: el bucle de Newton pasa al corrector compartido `NewtonCorrector`; este solver aporta su problema por paso (`_ArcProblem`: iterado `(U, λ, ΔU)`, dos resoluciones por iteración y restricción cilíndrica como método `constraint`; raíces imaginarias → `CorrectionAborted`; predictor tangente en `_tangent_predictor`; retirados `_solve` y `_make_linalg`) y conserva su control de paso. Sin cambio de formulación ni de resultados (suite completa sin tocar ningún valor esperado).
- **2026-09-23** · ADR 0017-0019: `linear_algebra` admite `pardiso` e `iterative` (CG → MINRES automático ante curvatura negativa, adecuado a la tangente indefinida del régimen postcrítico); al empezar `solve` se rechaza un mecanismo rígido (`MechanismError`). Dentro del trazado no se rechaza ningún sistema casi singular: cerca de un punto límite es parte del algoritmo.
- **2026-09-23** · Revisión documental: `tests/test_arclength.py` ya no existe; la cobertura está repartida en `test_solver_robustness.py`, `test_snap_through_corot.py` y `test_integration.py`.
- **2026-09-23** · **Cambio de formulación validado por el usuario** (deuda #26, destapada en el ejemplo 4 del manual de ejemplos). (1) Primer paso adimensional: `initial_dlambda` (0.1 por omisión) y $\Delta l_1 = \Delta\lambda_1\lVert\mathbf K_0^{-1}\mathbf F_{\text{ref}}\rVert$; `initial_dl` pasa a ser una longitud explícita opcional, excluyente, sin valor por omisión. Antes el valor por omisión era `initial_dl = 0.1`, una longitud: en el cilindro J2 del ejemplo 4 pedía $\Delta\lambda = 86$ en el primer paso, el predictor rebasaba `max_lambda` y el solver imponía $\lambda = 1.5$ en control de carga, por encima del colapso, sobre un equilibrio espurio del modelo discreto (0,98 m de desplazamiento). (2) Llegada exacta a `max_lambda` dentro del tramo recorrido, en lugar del paso final decidido por el predictor (§10). Los 84 tests existentes de los dos solvers de arco pasan sin cambiar ningún valor esperado; los nuevos fallan con el código anterior (9 de 11).
