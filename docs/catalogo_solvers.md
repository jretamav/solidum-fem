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
- **Parámetros**: `convergence`, `max_iter`, `num_steps`, `adaptive`, `min_delta_lambda`, `linear_algebra`, `freeze_tangent_after_iter`, **`line_search`** (default `False`, ADR 0011).
- **Cuándo usarlo**: la opción por defecto para no-linealidad material o geométrica suave (sin snap-back).
- **Cuándo activar `line_search=True`**: cuando se observe oscilación del residuo entre iteraciones (síntoma clásico: ratios alternados sin descenso monótono). Caso documentado: plasticidad perfecta con carga cerca/sobre la capacidad. **Default `False`** porque en problemas con tangente consistente cuasi-cuadrática (daño activo, plasticidad estándar) el line search rechaza pasos correctos del Newton donde el residuo sube transitoriamente — ver ADR 0011 §"Enmiendas".
- **Diagnóstico al diverger**: lanza una subclase tipada de `RuntimeError` (`OscillatingNewtonError`, `SingularTangentError`, `LoadExceedsCapacityError`, `UnknownDivergenceError`) con métricas (último ‖R‖, último ‖δU‖, factor de carga, bisecciones consumidas) y `hint` textual. Ver `fenix/math/solvers/diagnostics.py`.
- **Limitación**: con bisección adaptativa y cinemática no lineal (corotacional) puede atravesar puntos límite suaves; no atraviesa snap-back con `du/dλ < 0`. Para eso → `ArcLengthSolver`. Hallazgo de la auditoría fase A: el solver es más capaz de lo asumido históricamente (ver [`docs/auditorias/solvers_robustez_fase_A.md`](auditorias/solvers_robustez_fase_A.md)).
- **Referencia**: Crisfield, *Non-linear Finite Element Analysis of Solids and Structures*, vol. 1, cap. 9. Line search: Grippo-Lampariello-Lucidi 1986 (variante de descenso no monótono).
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

## NewtonNewmarkSolver — análisis dinámico transitorio no lineal (ADR 0009 fase 4)

- **Propósito**: integrar en el tiempo `M·ü + C·u̇ + F_int(u) = F_ext(t)` con materiales no lineales (plasticidad, daño) o no linealidad geométrica. **Variante** de `NewmarkSolver` (Reglas §4): hereda predictores/correctores y reducción de Dirichlet; añade Newton-Raphson dentro de cada paso temporal.
- **Esquema**:
  - Predictores Newmark idénticos a la fase 3.
  - **Bucle interno de Newton-Raphson** sobre el residuo dinámico
    `R(ü_{n+1}) = F_ext(t_{n+1}) − F_int(u_{n+1}) − C·u̇_{n+1} − M·ü_{n+1}`,
    con jacobiano `J = M + γΔt·C + βΔt²·K_t` y `K_t` la rigidez tangente corriente.
  - **Convergencia dual** (ADR 0007): `R/(atol + rtol·F_ref) ≤ 1` y `δu/(atol + rtol·‖u‖) ≤ 1`.
  - **Commit de estado**: tras converger cada paso, `assembler.commit_all_states()` (trial → committed). Si no converge en `max_iter`, lanza `RuntimeError` y se descarta el state trial.
  - **Amortiguamiento Rayleigh constante en el tiempo**, calibrado con la **rigidez elástica de referencia** `K_0` (al inicio del análisis, `u = 0`). Convención estándar (Abaqus, ANSYS, OpenSees); evita acoplamiento ad-hoc entre disipación viscosa y plástica.
  - **Newton modificado opcional** (`freeze_tangent_after_iter`, ADR 0003 fase 2): factoriza fresco las primeras N iter de cada paso y reusa la factorización en las siguientes.
  - **Recuperación del caso lineal**: con materiales lineales (`K_t ≡ K_0`), el residuo se anula en 1 iter y el resultado coincide exactamente con `NewmarkSolver`. Validado en tests.
- **Parámetros**: heredados de `NewmarkSolver` (`t_end`, `dt`, `beta`, `gamma`, `rayleigh`, `u0`, `u0_dot`, `F_func`, `linear_algebra`, `lumping`) más `convergence` (`ConvergenceCriterion`, ADR 0007), `max_iter` (default 20), `freeze_tangent_after_iter` (default `None`), **`line_search`** (default `False`, ADR 0011 — mismo criterio que `NonlinearSolver`).
- **Diagnóstico al diverger**: lanza una subclase tipada de `RuntimeError` con métricas y `hint` (mismo patrón que `NonlinearSolver`, ADR 0011).
- **Firma del solver**: `solve()` sin argumentos, retorna `TransientResult` (mismo formato que `NewmarkSolver`).
- **Cuándo usarlo**: respuesta dinámica de estructuras con plasticidad transitoria, fatiga de bajo ciclo, impacto en hormigón friccional, vibraciones de marcos con disipación plástica.
- **Limitaciones**: igual que `NewmarkSolver` en lo lineal (apoyos constantes, paso fijo). Además: no resuelve snap-back en problemas con softening severo — combinar con `ArcLengthSolver` si emerge.
- **Despacho YAML**: `solver.type: NewtonNewmarkSolver`. Como es subclase de `NewmarkSolver`, `entry.run_yaml` lo detecta vía `isinstance` y despacha a `run_transient` automáticamente.
- **Referencia**: ADR 0009 §Fase 4; Hughes, *FEM* cap. 7-9 (Newton dentro de paso temporal); Crisfield vol. 2 cap. 24 (dinámica no lineal).
- **Spec**: [docs/specs/NewtonNewmarkSolver.md](specs/NewtonNewmarkSolver.md)
- **Archivo**: [fenix/math/solvers/newmark.py](../fenix/math/solvers/newmark.py) (clase `NewtonNewmarkSolver`, junto al padre).

---

## HHTSolver — análisis dinámico transitorio lineal con disipación numérica (Hilber-Hughes-Taylor α)

- **Propósito**: integrar `M·ü + C·u̇ + K·u = F(t)` con **amortiguamiento numérico controlado** que disipa selectivamente los modos de alta frecuencia espurios sin afectar los modos resueltos. **Variante** de `NewmarkSolver` (Reglas §4): hereda predictores/correctores, reducción Dirichlet y factorización única de `A_eff`.
- **Esquema**:
  - Ecuación de equilibrio temporal: `M·ü_{n+1} + (1+α)·C·u̇_{n+1} − α·C·u̇_n + (1+α)·K·u_{n+1} − α·K·u_n = (1+α)·F_{n+1} − α·F_n`.
  - **Auto-derivación de β, γ desde α** (Hilber 1977): `β = (1−α)²/4`, `γ = (1−2α)/2`. Preserva orden 2 y estabilidad incondicional.
  - **Sistema efectivo**: `A_eff = M + (1+α)·γΔt·C + (1+α)·βΔt²·K`. Constante en el tiempo → factorización única reutilizable (como Newmark).
  - **Radio espectral en alta frecuencia**: `ρ_∞ = (1+α)/(1−α)`. Con `α=0` → `ρ_∞=1` (sin disipación, recupera Newmark trapezoidal). Con `α=−0.05` → `ρ_∞≈0.905`. Con `α=−1/3` → `ρ_∞=0.5` (disipación máxima manteniendo orden 2).
  - **Disipación selectiva**: modos con `ωΔt > π` se atenúan significativamente; modos con `ωΔt « 1` casi intactos. Es lo que distingue HHT-α de subir γ Newmark (que disipa en todas las frecuencias y baja a orden 1).
- **Parámetros**: `t_end`, `dt` obligatorios. **`alpha`** (default `-0.05`, rango `[−1/3, 0]`). `beta`, `gamma` autoderivados desde `alpha` salvo override explícito (no recomendado fuera de experimentación). Resto heredado de `NewmarkSolver`: `rayleigh`, `u0`, `u0_dot`, `F_func`, `linear_algebra`, `lumping`.
- **Cuándo usarlo**: cuando la respuesta transitoria contenga modos de alta frecuencia espurios (mallas con modos altos no de interés, impacto, contacto, propagación con `ωΔt` cerca o por encima de π) que enmascaren la señal. Si quieres exactamente Newmark trapezoidal sin disipación, usa `NewmarkSolver`.
- **Default de `alpha`**: `-0.05` es el valor canónico (Hilber 1977; Abaqus, OpenSees, ANSYS): disipación >9% por ciclo en altas frecuencias, <0.5% en modos de interés. Si pides `HHTSolver` con `alpha=0.0`, es funcionalmente idéntico a `NewmarkSolver(beta=0.25, gamma=0.5)`.
- **Despacho YAML**: `solver.type: HHTSolver`. Como subclase de `NewmarkSolver`, `entry.run_yaml` lo detecta vía `isinstance` y despacha a `run_transient`.
- **Referencia**: Hilber, Hughes & Taylor (1977), *Earthquake Eng. Struct. Dyn.* 5(3) 283-292; Hughes *FEM* §9.3; Chopra *Dynamics of Structures* §15.4.
- **Spec**: [docs/specs/HHTSolver.md](specs/HHTSolver.md)
- **Archivo**: [fenix/math/solvers/newmark.py](../fenix/math/solvers/newmark.py) (clase `HHTSolver`).

---

## NewtonHHTSolver — análisis dinámico transitorio no lineal HHT-α

- **Propósito**: combinar HHT-α (disipación numérica) con Newton-Raphson interno para problemas con materiales no lineales (plasticidad, daño) o no linealidad geométrica. **Variante** de `NewtonNewmarkSolver` (Reglas §4) que aplica el esquema temporal HHT-α al residuo dinámico no lineal.
- **Esquema**:
  - Residuo dinámico: `R(ü_{n+1}) = (1+α)·F_{n+1} − α·F_n − [(1+α)·F_int(u_{n+1}) − α·F_int(u_n)] − [(1+α)·C·u̇_{n+1} − α·C·u̇_n] − M·ü_{n+1}`.
  - Jacobiano: `J = M + (1+α)·γΔt·C + (1+α)·βΔt²·K_t`.
  - Todo lo demás idéntico a `NewtonNewmarkSolver`: convergencia dual (ADR 0007), commit/rollback, Rayleigh con K_0, Newton modificado opcional, **`line_search`** opt-in (ADR 0011), telemetría tipada al diverger.
  - **Recuperación**: con materiales lineales, el residuo de Newton se anula en una iteración y el resultado coincide con `HHTSolver`. Validado en tests.
- **Parámetros**: `alpha` (default `-0.05`) + heredados de `NewtonNewmarkSolver` (`convergence`, `max_iter`, `freeze_tangent_after_iter`, `line_search`, …).
- **Cuándo usarlo**: respuesta dinámica de estructuras con plasticidad transitoria + presencia de modos altos espurios que el Newmark trapezoidal no atenuaría.
- **Despacho YAML**: `solver.type: NewtonHHTSolver`. Vía atributo `PIPELINE_KIND="transient"` heredado de `NewmarkSolver` → `run_transient`.
- **Spec**: [docs/specs/HHTSolver.md](specs/HHTSolver.md) (cubre ambas variantes lineal y no lineal).
- **Archivo**: [fenix/math/solvers/newmark.py](../fenix/math/solvers/newmark.py) (clase `NewtonHHTSolver`).

---

## CentralDifferenceSolver — integración explícita por diferencias centradas (ADR 0009 fase 5)

- **Propósito**: integrar `M·ü + C·u̇ + F_int(u) = F(t)` con esquema **explícito** de segundo orden — la aceleración se actualiza en cada paso por `M⁻¹·rhs` con `M` diagonal lumped, sin sistema lineal ni Newton interno. Cubre análisis lineal y no lineal con un solo solver (parámetro `nonlinear: bool`).
- **Esquema (leapfrog Belytschko-Liu-Moran)**:
  - Inicialización: `a_0 = M⁻¹·(F(0) − F_int(u_0) − C·v_0 − F_dir)`, `v_{1/2} = v_0 + (Δt/2)·a_0`.
  - Paso `n → n+1`: `u_{n+1} = u_n + Δt·v_{n+1/2}`; reevaluar `F_int(u_{n+1})`; `a_{n+1} = M⁻¹·(F_{n+1} − F_int − C·v_{n+1/2} − F_dir)`; `v_{n+1} = v_{n+1/2} + (Δt/2)·a_{n+1}`.
  - Modo lineal: `F_int = K·u` (K constante, ensamblada una vez). Modo no lineal: `Assembler.assemble_non_linear_system(u)` cada paso (descarta la tangente).
- **Estabilidad condicional CFL**: `Δt < 2/ω_max`. Para barras: `Δt < L_el/c_p = L_el·√(ρ/E)`. Si se excede, el solver detecta la explosión exponencial a posteriori y aborta con `RuntimeError` informativo. **El solver no calcula `Δt_crit` automáticamente** — responsabilidad del usuario en esta versión inicial (power iteration sobre `M⁻¹K` queda diferida).
- **Parámetros**: `t_end`, `dt`, `nonlinear` (default `False`), `rayleigh`, `u0`, `u0_dot`, `F_func`, `lumping` (default `"lumped"`; `"consistent"` rechazado), `divergence_threshold` (default `1e8`).
- **Cuándo usarlo**: wave propagation, impacto, dinámica de respuesta rápida — situaciones donde el `Δt` natural del problema cae por debajo del paso de estabilidad incondicional de Newmark. Para problemas con plasticidad transitoria de respuesta moderadamente rápida, sigue siendo competitivo frente a Newton-Newmark por evitar la factorización de la tangente cada iteración.
- **Rechazos defensivos**: `lumping="consistent"` → `ValueError` (la inversión de M consistente cada paso anula la ventaja); Frame3D con eje oblicuo a globales → `ValueError` (M_red no es estrictamente diagonal; ADR 0009 fase 2 documenta la bloque-diagonalidad inherente).
- **Despacho YAML**: `solver.type: CentralDifferenceSolver`. Vía atributo `PIPELINE_KIND="transient"` → `run_transient`. **Primer solver fuera de la jerarquía de Newmark**, lo que disparó la regla C de refactor (atributo declarativo en lugar de cadena de isinstance).
- **Spec**: [docs/specs/CentralDifferenceSolver.md](specs/CentralDifferenceSolver.md).
- **Archivo**: [fenix/math/solvers/central_difference.py](../fenix/math/solvers/central_difference.py).

---

## Cómo añadir un solver nuevo

`/fenix-new solver <Name>` — genera archivo en `fenix/math/solvers/<snake>.py`, decorador `@SolverRegistry.register`, esqueleto de test.

Convenciones de interfaz:

- **Solvers estáticos** (lineales, no lineales, arc-length): `PIPELINE_KIND = "static"`. Constructor recibe `assembler` + parámetros; método `solve(F_ext_global, step_callback=None) → U_final`. Comprometen los estados internos vía `assembler.commit_all_states()` al converger cada paso. Retornan el campo de desplazamientos completo.
- **Solvers modales / autovalor** (modal — ADR 0009 — y futuros pandeo lineal): `PIPELINE_KIND = "modal"`. Constructor recibe `assembler` + parámetros; método `solve() → ModalResult` (u otro tipo específico). No consumen vector de cargas. `run_yaml` despacha a `run_modal`.
- **Solvers transitorios** (Newmark, HHT, central differences, futuros harmonic/response spectrum): `PIPELINE_KIND = "transient"`. Constructor recibe `assembler`, `t_end`, `dt` + parámetros; método `solve() → TransientResult`. `run_yaml` despacha a `run_transient`.

El **dispatch en `run_yaml`** se hace por el atributo de clase `PIPELINE_KIND` (regla C, 2026-05-18). Solvers nuevos no clásicos (que no hereden de los existentes) sólo declaran su `PIPELINE_KIND` y quedan automáticamente cableados — no requieren tocar `entry.py`.

Tras implementar, **añadir una entrada a este catálogo** siguiendo el formato de arriba y promover la spec a `status: validated`.
