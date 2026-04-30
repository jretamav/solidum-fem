# Catálogo de solvers

> Referencia rápida de los métodos de solución implementados. Una entrada por solver. Para detalles del algoritmo → código fuente.

---

## LinearSolver — solución lineal directa en un paso

- **Propósito**: resolver `K · U = F` en un único `spsolve`.
- **Esquema**: ensamblaje único → eliminación directa de DOFs prescritos (ADR 0004) → `scipy.sparse.linalg.spsolve`.
- **Parámetros**: `linear_algebra` (default `auto`; ADR 0003 §4).
- **Cuándo usarlo**: problemas estrictamente lineales (todos los materiales con tangente constante, sin grandes desplazamientos, sin contacto).
- **No converge / no aplica**: cualquier no-linealidad material (daño, plasticidad) o geométrica (corotacional).
- **Archivo**: [fenix/math/solvers.py](fenix/math/solvers.py)

---

## NonlinearSolver — Newton-Raphson incremental con paso adaptativo

- **Propósito**: resolver problemas no lineales con control de carga, dividiendo la carga total en pasos y aplicando Newton-Raphson dentro de cada paso.
- **Esquema**:
  - Bucle externo: incrementar factor de carga `λ` desde 0 hasta 1 en `num_steps` pasos.
  - Bucle interno: Newton-Raphson con `K_t · ΔU = R = λ·F_ext − F_int`.
  - **Criterio de convergencia dual**: `max(err_disp, err_force) < tol`, con
    `err_disp = ‖ΔU‖ / (‖U‖ + ε)` y `err_force = ‖R‖ / max(‖λ·F_ext‖, ‖F_int‖, ε)`.
  - **Adaptatividad** (`adaptive=True`): si converge en < 5 iteraciones, agranda el siguiente paso (×1.5); si no converge, biseca (÷2). Falla si `Δλ < min_delta_lambda`.
  - Tras converger un paso → `assembler.commit_all_states()` (trial → committed en todos los elementos).
- **Parámetros**: `tol`, `max_iter`, `num_steps`, `adaptive`, `min_delta_lambda`, `linear_algebra`, `freeze_tangent_after_iter`.
- **Cuándo usarlo**: la opción por defecto para no-linealidad material o geométrica suave (sin snap-back/snap-through).
- **Limitación**: no puede atravesar puntos límite (snap-through) ni recorrer ramas con derivada negativa de la curva carga-desplazamiento. Para eso → `ArcLengthSolver`.
- **Referencia**: Crisfield, *Non-linear Finite Element Analysis of Solids and Structures*, vol. 1, cap. 9.
- **Archivo**: [fenix/math/solvers.py](fenix/math/solvers.py)

---

## ArcLengthSolver — método de longitud de arco cilíndrico (Crisfield)

- **Propósito**: trazar curvas de equilibrio con snap-through, snap-back o pérdida de unicidad de carga, controlando simultáneamente desplazamientos y factor de carga.
- **Esquema**:
  - **Predictor tangente**: `du_t = K_t⁻¹ · F_ext_ref`; `dλ = sign · dl / ‖du_t‖`.
    El `sign` se elige por proyección con el incremento del paso anterior (evitar regresar).
  - **Corrector iterativo**: en cada iteración resuelve dos sistemas (`du_R = K⁻¹·R` y `du_t = K⁻¹·F_ext`) y aplica la restricción cuadrática cilíndrica de Crisfield: `‖dU‖² = dl²`.
  - De las dos raíces se elige la que mantiene el ángulo positivo con el incremento previo.
  - **Caso especial — final_step**: si el predictor sobrepasaría `max_lambda`, fija λ exactamente a `max_lambda` y resuelve solo desplazamientos (Newton-Raphson puro).
  - **Auto-ajuste de `dl`**: < 4 iter → ampliar (×1.5, tope `5·dl_initial`); > 8 iter → reducir (×0.6); no converge → biseca (÷2).
  - Convergencia: mismo criterio dual que `NonlinearSolver`.
- **Parámetros**: `tol`, `max_iter`, `max_lambda`, `initial_dl`, `max_steps`, `dl_grow_factor`, `dl_max_factor`, `dl_shrink_factor`, `dl_grow_iter_threshold`, `dl_shrink_iter_threshold`.
- **Cuándo usarlo**: problemas con softening pronunciado (daño, post-pandeo, snap-through de cúpulas), o cuando `NonlinearSolver` diverge cerca de un punto límite.
- **Limitación**: más caro por paso (dos `spsolve` por iteración); requiere ajuste de `initial_dl` para problemas nuevos.
- **Referencia**: Crisfield, "A fast incremental/iterative solution procedure that handles snap-through" (Computers & Structures, 1981); Crisfield vol. 1, cap. 9.
- **Archivo**: [fenix/math/solvers.py](fenix/math/solvers.py)

---

## Cómo añadir un solver nuevo

`/fenix-new solver <Name>` — genera archivo en `fenix/math/solver_<snake>.py`, decorador `@SolverRegistry.register`, esqueleto de test.
Convención de interfaz: constructor recibe `assembler` + parámetros; método `solve(F_ext_global, step_callback=None)` devuelve `U_final` y compromete los estados internos vía `assembler.commit_all_states()` al converger cada paso.
Tras implementar, **añadir una entrada a este catálogo** siguiendo el formato de arriba.
