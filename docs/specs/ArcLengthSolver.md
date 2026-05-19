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
dl = initial_dl
mientras λ < max_lambda y step < max_steps:
    # Predictor tangente
    K_t, F_int = assemble_non_linear_system(U_iter)
    δu_t = K_t^{-1} · F_ext^ref          # desplazamiento tangente unitario
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
            commit_all_states()
            λ += Δλ
            ajustar_dl_adaptativamente(iter)
            break
    si no converge:
        dl /= dl_shrink_factor
```

### 5. Predictor / corrector

**Predictor tangente**: $\delta\mathbf u_t = \mathbf K_t^{-1}\,\mathbf F_{\text{ext}}^{\text{ref}}$ (vector que mide la pendiente local del camino de equilibrio). $\Delta\lambda_{\text{pred}} = \pm dl / \lVert\delta\mathbf u_t\rVert$ con signo determinado por el producto escalar contra el incremento previo $\Delta\mathbf U^{\text{paso prev}}$ — así no se retrocede en el camino tras pasar un punto límite.

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

- Converge en $<$ `dl_grow_iter_threshold` iteraciones ⇒ $dl \to dl \cdot $ `dl_grow_factor` (cota: `initial_dl · dl_max_factor`).
- Converge en $>$ `dl_shrink_iter_threshold` iteraciones ⇒ $dl \to dl \cdot $ `dl_shrink_factor`.
- No converge en `max_iter` ⇒ $dl /= 2$ y reintentar el paso.
- Cota inferior: $dl < $ `ARCLENGTH_MIN_DL_FACTOR · initial_dl` ⇒ aborto.

### 10. Caveats numéricos

- **Cancelación o bifurcación verdadera**: en bifurcaciones simétricas múltiples, la selección de raíz por menor ángulo puede saltar a un camino paralelo. Para análisis de bifurcación auténtica usar perturbación + post-procesado modal.
- **Raíces imaginarias**: si el discriminante de la cuadrática es negativo (paso demasiado grande, mal condicionamiento severo), el solver aborta el paso y bisecta $dl$.
- **Último paso**: cuando $\lambda + \Delta\lambda_{\text{pred}} \ge \lambda_{\max}$, el solver fija $\Delta\lambda = \lambda_{\max} - \lambda$ y corrige solo en desplazamientos (Newton puro), sin restricción cuadrática. Permite cerrar exactamente en el target.
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
  - { name: initial_dl,               type: float, required: false, default: 0.1,
      desc: "Longitud de arco inicial" }
  - { name: max_steps,                type: int,   required: false, default: 100,
      desc: "Cota dura sobre el número de pasos" }
  - { name: dl_grow_factor,           type: float, required: false, default: 1.5 }
  - { name: dl_max_factor,            type: float, required: false, default: 5.0,
      desc: "dl máximo = initial_dl · dl_max_factor" }
  - { name: dl_shrink_factor,         type: float, required: false, default: 0.6 }
  - { name: dl_grow_iter_threshold,   type: int,   required: false, default: 4 }
  - { name: dl_shrink_iter_threshold, type: int,   required: false, default: 8 }
  - { name: linear_algebra,           type: str,   required: false, default: "auto",
      desc: "'auto' (LU para postcrítico), 'cholesky', 'lu'. En arc-length el default fuerza LU por indefinitud potencial" }

requirements:
  - "Modelo con potencial snap-through/snap-back (geometría inestable, softening del material)"
  - "Apoyos suficientes; sin modos rígidos"
  - "dl inicial razonable (orden de magnitud del desplazamiento característico del primer paso)"

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

- **Archivo**: [fenix/math/solvers/arclength.py](../../fenix/math/solvers/arclength.py).
- **Clase**: `ArcLengthSolver`, registrada vía `@SolverRegistry.register` con `PIPELINE_KIND = "static"`.
- **Restricción cuadrática y selección de raíz**: implementada en `solve`. Comparación de ángulos `theta1 = ΔU_iter · (ΔU_new + ddl1·δu_t)`, idem `theta2`, seleccionar el mayor → menor ángulo con la dirección previa.
- **`_make_linalg`**: fuerza `is_positive_definite=False` independientemente de la simetría del dominio (régimen postcrítico).
- **Entrypoint público**: `fenix.run_static(model, solver="arclength", ...)`.
- **Tests**:
  - [tests/test_arclength.py](../../tests/test_arclength.py) · snap-through Lee frame, snap-back truss, recuperación lineal.
  - [tests/test_solver_robustness.py](../../tests/test_solver_robustness.py) · `test_arc_length_traverses_damage_softening`.

---

## Diálogo

- **2026-05-19** · Spec creada retroactivamente para cerrar el hueco H-5.3. Solver anterior a la convención de specs. La spec recoge la implementación tal como está al cierre de la sesión de saneamiento post-auditoría.
