# Catálogo de solvers

> Referencia rápida de los métodos de solución implementados. Una entrada por solver. Para detalles del algoritmo → código fuente.

---

## LinearSolver — solución lineal directa en un paso

- **Propósito**: resolver `K · U = F` en un único `spsolve`.
- **Esquema**: ensamblaje único → eliminación directa de DOFs prescritos (ADR 0004) → `scipy.sparse.linalg.spsolve`.
- **Parámetros**: `linear_algebra` (default `auto`; ADR 0003 §4).
- **Cuándo usarlo**: problemas estrictamente lineales (todos los materiales con tangente constante, sin grandes desplazamientos, sin contacto).
- **No converge / no aplica**: cualquier no-linealidad material (daño, plasticidad) o geométrica (corotacional).
- **Archivo**: [fenix/math/solvers/linear.py](../fenix/math/solvers/linear.py)

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
- **Archivo**: [fenix/math/solvers/nonlinear.py](../fenix/math/solvers/nonlinear.py)

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
- **Archivo**: [fenix/math/solvers/arclength.py](../fenix/math/solvers/arclength.py)

---

## ModalSolver — análisis modal por autovalores generalizados (ADR 0009)

- **Propósito**: calcular frecuencias naturales `ω_n` y modos `φ_n` de vibración no amortiguados resolviendo `K · φ = ω² · M · φ`.
- **Esquema**:
  - Ensamblaje de `K` (rigidez lineal en `u = 0`) y `M` (masa consistente, ADR 0009 §1) reusando la topología COO cacheada del `Assembler`.
  - Reducción simultánea por Dirichlet `Assembler.reduce_pair(K, M)` → `K_red, M_red, T`.
  - `EigenSolver` envuelve `scipy.sparse.linalg.eigsh` con **shift-invert** en `sigma` (`sigma=0` → frecuencias más bajas) y `which="LM"`. ARPACK factoriza `(K_red − σ·M_red)` una sola vez y la reusa en todas las iteraciones Lanczos.
  - Expansión de modos al espacio completo: `φ = T · φ_red`. En DOFs prescritos por Dirichlet la componente es 0.
  - Salida en `ModalResult`: `frequencies_rad`, `frequencies_hz`, `periods`, `modes` (M-ortonormales), `n_modes`, `converged`.
- **Parámetros**: `n_modes` (obligatorio), `sigma` (default 0.0), `which` (default `"LM"`), `tolerance` (default 1.0e-9), `lumping` (default `"consistent"`; fase 1 solo admite ese valor), `linear_algebra` (reservado; ADR 0003).
- **Firma del solver**: `solve()` sin argumentos (no consume `F_applied`). El entrypoint público es `fenix.run_modal(domain, n_modes=N)` o YAML con `solver: type: ModalSolver`. `fenix.run_yaml` despacha automáticamente al detectar el tipo.
- **Cuándo usarlo**: caracterización dinámica de la estructura (frecuencias naturales, formas modales) antes de pasar a análisis transitorio o sísmico. Validación de la matriz de masa contra fórmulas analíticas (barra empotrada-libre, viga Bernoulli-Euler).
- **Limitaciones**: lineal (`K` evaluada en `u = 0`); sin amortiguamiento; sin masas modales efectivas ni factores de participación (fase 7 del ADR 0009); modos complejos no soportados.
- **Referencia**: Bathe, *Finite Element Procedures*, cap. 9 (masa) y cap. 11 (algoritmos de autovalor); Hughes, *The Finite Element Method*, cap. 7; ARPACK Users' Guide (Lehoucq–Sorensen–Yang).
- **Spec**: [docs/specs/ModalSolver.md](specs/ModalSolver.md)
- **Archivos**: [fenix/math/solvers/modal.py](../fenix/math/solvers/modal.py) (clase `ModalSolver`), [fenix/math/linalg/eigen.py](../fenix/math/linalg/eigen.py) (clase `EigenSolver`), [fenix/results.py](../fenix/results.py) (`ModalResult`).

---

## NewmarkSolver — análisis dinámico transitorio lineal (ADR 0009 fase 3)

- **Propósito**: integrar en el tiempo la ecuación de movimiento semidiscreta `M·ü + C·u̇ + K·u = F(t)` con amortiguamiento Rayleigh `C = α·M + β·K` y apoyos Dirichlet constantes.
- **Esquema**:
  - Familia Newmark-β con `(β, γ)` parametrizables; default `(1/4, 1/2)` — average acceleration, incondicionalmente estable, sin amortiguamiento numérico, error `O(Δt²)`.
  - **Predictores**: `ũ = u_n + Δt·u̇_n + (Δt²/2)(1−2β)·ü_n`, `ũ̇ = u̇_n + Δt·(1−γ)·ü_n`.
  - **Sistema efectivo**: `A_eff·ü_{n+1} = F_{n+1} − C·ũ̇ − K·ũ` con `A_eff = M + γΔt·C + βΔt²·K`.
  - **Correctores**: `u_{n+1} = ũ + βΔt²·ü_{n+1}`, `u̇_{n+1} = ũ̇ + γΔt·ü_{n+1}`.
  - Reducción por Dirichlet (operador `T`) aplicada a `(K, M, C)` simultáneamente; factorización única de `A_eff_red` reusada en todos los pasos.
  - Aceleración inicial consistente con la ecuación de movimiento en `t=0`.
- **Parámetros**: `t_end`, `dt` (obligatorios), `beta` (default 0.25), `gamma` (default 0.5), `rayleigh` (None | `{alpha, beta}` directos | `{xi1, omega1, xi2, omega2}` para calibración modal), `u0`, `u0_dot` (default ceros), `F_func` (callback Python `t → F(t)`; default vibración libre), `linear_algebra`, `lumping`.
- **Firma del solver**: `solve()` sin argumentos, retorna `TransientResult` con historiales `t_history`, `u_history`, `udot_history`, `uddot_history` (shape `(n_dof, n_steps + 1)`).
- **Cuándo usarlo**: respuesta a cargas dependientes del tiempo (sísmica time-history, impacto, vibración forzada), evolución transitoria desde condiciones iniciales arbitrarias, validación de respuesta libre/forzada de estructuras lineales.
- **Limitaciones**:
  - Solo lineal (`K`, `M` constantes). La extensión a no-linealidad (Newton dentro de cada paso) es fase 4 del ADR 0009.
  - Apoyos Dirichlet constantes en el tiempo. Multi-support excitation sísmica diferida.
  - Paso `Δt` fijo. Adaptativo diferido.
  - HHT-α y generalized-α (con amortiguamiento numérico controlado) no incluidos.
- **Referencia**: Newmark (1959), J. Eng. Mech. Div. ASCE; Bathe, *FEP* cap. 9.4; Hughes, *FEM* cap. 9; Chopra, *Dynamics of Structures* cap. 5 y 15.
- **Spec**: [docs/specs/NewmarkSolver.md](specs/NewmarkSolver.md)
- **Archivos**: [fenix/math/solvers/newmark.py](../fenix/math/solvers/newmark.py) (clase `NewmarkSolver`), [fenix/math/damping.py](../fenix/math/damping.py) (Rayleigh), [fenix/results.py](../fenix/results.py) (`TransientResult`).

---

## Cómo añadir un solver nuevo

`/fenix-new solver <Name>` — genera archivo en `fenix/math/solvers/<snake>.py`, decorador `@SolverRegistry.register`, esqueleto de test.

Convenciones de interfaz:

- **Solvers estáticos** (lineales, no lineales, arc-length): constructor recibe `assembler` + parámetros; método `solve(F_ext_global, step_callback=None) → U_final`. Comprometen los estados internos vía `assembler.commit_all_states()` al converger cada paso. Retornan el campo de desplazamientos completo.
- **Solvers de autovalor** (modal — ADR 0009 — y futuros pandeo lineal): constructor recibe `assembler` + parámetros; método `solve() → ModalResult` (u otro tipo de resultado específico). No consumen vector de cargas. El entrypoint `fenix.run_yaml` despacha por `isinstance` al entrypoint específico (`fenix.run_modal` para modal).

Tras implementar, **añadir una entrada a este catálogo** siguiendo el formato de arriba y promover la spec a `status: validated`.
