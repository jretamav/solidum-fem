# Análisis Dinámico

Fenix FEM expone tres análisis dinámicos cerrados por la fase 4 del ADR 0009: modal (autovalores generalizados), transitorio lineal (Newmark-$\beta$) y transitorio no lineal (Newton-Newmark). Todos consumen la matriz de masa consistente que cada elemento devuelve por `compute_mass_matrix(lumping="consistent")` y la densidad declarada en el material (`density`; obligatoria para cualquier análisis dinámico, ADR 0008).

## `ModalSolver` — análisis modal por autovalores generalizados

Resuelve el problema generalizado $\mathbf K \boldsymbol\phi = \omega^2 \mathbf M \boldsymbol\phi$ con $\mathbf K$ evaluada en $\mathbf u = 0$ (régimen lineal) y $\mathbf M$ ensamblada por cuadratura consistente sobre todos los elementos. La capa algebraica (`EigenSolver`) envuelve `scipy.sparse.linalg.eigsh` con *shift-invert* en $\sigma$ (ARPACK Lanczos), de modo que basta una sola factorización de $\mathbf K - \sigma\mathbf M$ para extraer las $n$ frecuencias más bajas.

**Parámetros**: `n_modes` (obligatorio), `sigma` (default 0.0; `sigma=0` extrae las frecuencias más bajas), `which` (default `"LM"`), `tolerance` (default `1.0e-9`), `lumping` (default `"consistent"`; la masa concentrada queda diferida a fase 2 del ADR 0009).

**Salida**: instancia de `ModalResult` con `frequencies_rad`, `frequencies_hz`, `periods`, `modes` (M-ortonormales). El método `ModalResult.free_vibration(M, u0, u0_dot, t)` reconstruye analíticamente la respuesta libre no amortiguada por superposición modal.

```yaml
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}

materials:
  - {id: 1, type: Elastic1D, E: 200.0e9, density: 7850.0}

elements:
  - {id: 1, type: Truss2D, material: 1, A: 0.01, nodes: [1, 2]}

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0}
  - {node_id: 2, uy: 0.0}

solver:
  type: ModalSolver
  n_modes: 3
  sigma: 0.0
```

```callout Nota
La densidad es *obligatoria* para todo material involucrado en análisis modal o transitorio. La omisión lanza `ValueError` (ADR 0008) listando los materiales afectados por nombre antes de iniciar el cálculo.
```

## `NewmarkSolver` — transitorio lineal Newmark-$\beta$

Integra en el tiempo $\mathbf M\ddot{\mathbf u} + \mathbf C\dot{\mathbf u} + \mathbf K\mathbf u = \mathbf F(t)$ por la familia Newmark-$\beta$. Con $(\beta, \gamma) = (1/4, 1/2)$ — *average acceleration*, default — es *incondicionalmente estable*, sin amortiguamiento numérico, con error $\mathcal O(\Delta t^2)$.

**Esquema operativo**:

- Predictor: $\tilde{\mathbf u} = \mathbf u_n + \Delta t\,\dot{\mathbf u}_n + (\Delta t^2/2)(1 - 2\beta)\ddot{\mathbf u}_n$ y análogo para $\tilde{\dot{\mathbf u}}$.
- Sistema efectivo en aceleración: $\mathbf A_{\text{eff}}\ddot{\mathbf u}_{n+1} = \mathbf F_{n+1} - \mathbf C\tilde{\dot{\mathbf u}} - \mathbf K\tilde{\mathbf u}$ con $\mathbf A_{\text{eff}} = \mathbf M + \gamma\Delta t\,\mathbf C + \beta\Delta t^2\,\mathbf K$. Factorización *única* al inicio del análisis y reusada en todos los pasos (`FactorizedSolver`, ADR 0003 fase 2).
- Corrector: $\mathbf u_{n+1} = \tilde{\mathbf u} + \beta\Delta t^2\,\ddot{\mathbf u}_{n+1}$, $\dot{\mathbf u}_{n+1} = \tilde{\dot{\mathbf u}} + \gamma\Delta t\,\ddot{\mathbf u}_{n+1}$.

**Amortiguamiento Rayleigh** $\mathbf C = \alpha\mathbf M + \beta\mathbf K$, declarado de dos formas:

- Directa: `rayleigh: {alpha: 0.1, beta: 1.0e-4}`.
- Calibración modal: `rayleigh: {xi1: 0.02, omega1: 50.0, xi2: 0.02, omega2: 200.0}` — Fenix resuelve los $(\alpha, \beta)$ que reproducen los amortiguamientos $\xi_1, \xi_2$ a las frecuencias $\omega_1, \omega_2$.

**Parámetros**: `t_end`, `dt` (obligatorios); `beta` (default 0.25), `gamma` (default 0.5), `rayleigh` (opcional), `u0`, `u0_dot` (condiciones iniciales, default ceros), `F_func` (callback Python $t \to \mathbf F(t)$; default vibración libre).

**Salida**: `TransientResult` con `t_history`, `u_history`, `udot_history`, `uddot_history` (forma $(n_{\text{dof}}, n_{\text{steps}} + 1)$) y los $\alpha_{\text{Rayleigh}}, \beta_{\text{Rayleigh}}$ efectivos.

```yaml
solver:
  type: NewmarkSolver
  t_end: 2.0
  dt: 0.005
  beta: 0.25
  gamma: 0.5
  rayleigh:
    xi1: 0.02
    omega1: 50.0
    xi2: 0.02
    omega2: 200.0
```

**Limitaciones**: lineal ($\mathbf K, \mathbf M$ constantes; para no-linealidad material o geométrica → `NewtonNewmarkSolver`). Apoyos Dirichlet constantes en el tiempo. Paso $\Delta t$ fijo (adaptativo diferido). HHT-$\alpha$ y generalized-$\alpha$ (con amortiguamiento numérico controlado) no incluidos.

## `NewtonNewmarkSolver` — transitorio no lineal Newton-Newmark

Subclase de `NewmarkSolver` (Reglas §4 — variante de componente existente) que añade un bucle Newton-Raphson dentro de cada paso temporal sobre el residuo dinámico:

$$\mathbf R = \mathbf F_{\text{ext}}(t) - \mathbf F_{\text{int}}(\mathbf u) - \mathbf C\dot{\mathbf u} - \mathbf M\ddot{\mathbf u}, \quad \mathbf J = \mathbf M + \gamma\Delta t\,\mathbf C + \beta\Delta t^2\,\mathbf K_t$$

**Convergencia dual** (ADR 0007): mismas tolerancias `rtol_force` / `rtol_disp` y factores `atol_*_factor` que los solvers no lineales estáticos. Tras converger cada paso se invoca `assembler.commit_all_states()` (trial → committed en todos los elementos); si agota `max_iter` se lanza `RuntimeError` y se descarta el estado trial.

**Amortiguamiento Rayleigh constante en el tiempo**, calibrado con la rigidez elástica de referencia $\mathbf K_0$ al inicio del análisis. Convención estándar (Abaqus, ANSYS, OpenSees) que evita acoplamiento ad-hoc entre disipación viscosa y plástica.

**Newton modificado opcional** (`freeze_tangent_after_iter`, ADR 0003 fase 2): factoriza fresco las primeras $N$ iter y reusa la factorización en las siguientes para abaratar pasos largos.

**Recuperación del caso lineal**: con materiales sin historia o no plastificados, el residuo se anula en 1 iter y el resultado coincide *a paridad de bits* con `NewmarkSolver`. Validado en tests.

**Parámetros**: heredados de `NewmarkSolver` más `max_iter` (default 20), bloque `convergence` (`rtol_force`, `rtol_disp`, `atol_force_factor`, `atol_disp_factor`; ADR 0007) y `freeze_tangent_after_iter` (default `None`). El despacho YAML es automático: como `NewtonNewmarkSolver` es subclase de `NewmarkSolver`, `fenix.run_yaml` lo detecta vía `isinstance` y enruta a `run_transient`.

```yaml
materials:
  - {id: 1, type: VonMises2D, E: 200.0e9, nu: 0.3, sigma_y: 250.0e6, H: 2.0e9, hypothesis: plane_strain, density: 7850.0}

solver:
  type: NewtonNewmarkSolver
  t_end: 0.5
  dt: 0.001
  beta: 0.25
  gamma: 0.5
  rayleigh:
    xi1: 0.02
    omega1: 100.0
    xi2: 0.02
    omega2: 500.0
  max_iter: 15
  convergence:
    rtol_force: 1.0e-6
    rtol_disp: 1.0e-6
```

**Cuándo usarlo**: respuesta dinámica de estructuras con plasticidad transitoria (sísmica con disipación plástica, impacto en hormigón friccional, fatiga de bajo ciclo), vibraciones de marcos elastoplásticos. Materiales con historia disponibles en dinámica: `Elastoplastic1D`, `VonMises2D` (*plane strain* y *plane stress*), `DruckerPrager2D`, `IsotropicDamage1D/2D`.

**Limitaciones**: igual que `NewmarkSolver` en lo lineal (apoyos constantes, paso fijo, sin HHT-$\alpha$). Además: no resuelve *snap-back* dinámico con softening severo — combinar con `ArcLengthSolver` si emerge.
