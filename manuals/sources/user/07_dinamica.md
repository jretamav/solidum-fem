# Análisis Dinámico

Solidum FEM expone el subsistema modal/dinámico/espectral **completo** (ADR 0009 cerrado en su totalidad 2026-05-18). Cubre los siete regímenes canónicos de la mecánica computacional clásica:

1. **Modal**: autovalores generalizados $\mathbf K \boldsymbol\phi = \omega^2\mathbf M\boldsymbol\phi$ (`ModalSolver`).
2. **Mass lumping**: matriz de masa diagonal vía HRZ canónico (sólidos) o nodal directo (frames). Disponible en todos los elementos vía `compute_mass_matrix(lumping="lumped")`.
3. **Transitorio implícito lineal**: Newmark-$\beta$ (`NewmarkSolver`) y variante HHT-$\alpha$ con disipación numérica (`HHTSolver`).
4. **Transitorio implícito no lineal**: Newton-Newmark (`NewtonNewmarkSolver`) y Newton-HHT (`NewtonHHTSolver`).
5. **Transitorio explícito**: diferencias centradas leapfrog (`CentralDifferenceSolver`).
6. **Frecuencia**: respuesta forzada armónica (`HarmonicSolver`).
7. **Espectral / sísmico**: combinación modal SRSS/CQC contra espectros normativos (`ResponseSpectrumSolver`).

Todos los solvers consumen la masa global que cada elemento devuelve por `compute_mass_matrix(lumping)` con `lumping ∈ {"consistent", "lumped"}` y la densidad declarada en el material (`density`; obligatoria para cualquier análisis dinámico, ADR 0008).

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
- Calibración modal: `rayleigh: {xi1: 0.02, omega1: 50.0, xi2: 0.02, omega2: 200.0}` — Solidum resuelve los $(\alpha, \beta)$ que reproducen los amortiguamientos $\xi_1, \xi_2$ a las frecuencias $\omega_1, \omega_2$.

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

**Parámetros**: heredados de `NewmarkSolver` más `max_iter` (default 20), bloque `convergence` (`rtol_force`, `rtol_disp`, `atol_force_factor`, `atol_disp_factor`; ADR 0007) y `freeze_tangent_after_iter` (default `None`). El despacho YAML es automático: como `NewtonNewmarkSolver` es subclase de `NewmarkSolver`, `solidum.run_yaml` lo detecta vía `isinstance` y enruta a `run_transient`.

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

**Limitaciones**: igual que `NewmarkSolver` en lo lineal (apoyos constantes, paso fijo). No resuelve *snap-back* dinámico con softening severo — combinar con `ArcLengthSolver` si emerge.

## `HHTSolver` y `NewtonHHTSolver` — variantes Hilber-Hughes-Taylor

Variantes de Newmark con **disipación numérica controlada** en altas frecuencias (Hilber, Hughes & Taylor 1977). Útiles cuando modos altos espurios contaminan la respuesta y conviene atenuarlos preservando segundo orden y estabilidad incondicional.

El parámetro $\alpha \in [-1/3, 0]$ controla la disipación; los coeficientes $(\beta, \gamma)$ se auto-derivan canónicamente:

$$\beta = \frac{(1-\alpha)^2}{4}, \quad \gamma = \frac{1-2\alpha}{2}$$

con radio espectral $\rho_\infty = (1+\alpha)/(1-\alpha)$. Recuperación exacta a Newmark cuando $\alpha = 0$.

**Esquema operativo**: idéntico a Newmark pero con factor $(1+\alpha)$ en los términos elásticos y de amortiguamiento del residuo dinámico, y $-\alpha$ multiplicando los del paso anterior. Variantes lineal (`HHTSolver`) y no lineal (`NewtonHHTSolver`) heredan toda la maquinaria de `NewmarkSolver` / `NewtonNewmarkSolver` (Reglas §4 sobre variantes).

```yaml
solver:
  type: HHTSolver           # o NewtonHHTSolver para materiales con historia
  t_end: 2.0
  dt: 0.005
  alpha: -0.05              # default; disipación ligera de altas frecuencias
```

## `CentralDifferenceSolver` — explícito por diferencias centradas (ADR 0009 fase 5)

Integración **explícita** de segundo orden por leapfrog Belytschko-Liu-Moran. Sin sistema lineal en cada paso — la única "inversión" es $\mathbf M^{-1}$, trivial con masa lumped diagonal. Apropiado para wave propagation, impacto, dinámica de respuesta rápida.

**Esquema operativo** (forma predictor-corrector):

$$\mathbf a_0 = \mathbf M^{-1}(\mathbf F(0) - \mathbf F_{\text{int}}(\mathbf u_0) - \mathbf C\mathbf v_0); \quad \mathbf v_{1/2} = \mathbf v_0 + \tfrac{\Delta t}{2}\mathbf a_0$$

Paso $n \to n+1$:

$$\mathbf u_{n+1} = \mathbf u_n + \Delta t\,\mathbf v_{n+1/2}; \quad \mathbf a_{n+1} = \mathbf M^{-1}(\mathbf F_{n+1} - \mathbf F_{\text{int}}(\mathbf u_{n+1}) - \mathbf C\mathbf v_{n+1/2})$$

$$\mathbf v_{n+1} = \mathbf v_{n+1/2} + \tfrac{\Delta t}{2}\mathbf a_{n+1}$$

Una sola clase cubre lineal y no lineal con parámetro `nonlinear`: con `False` usa $\mathbf F_{\text{int}} = \mathbf K\mathbf u$ con $\mathbf K$ constante; con `True` recalcula $\mathbf F_{\text{int}}$ cada paso vía `Assembler.assemble_non_linear_system`.

**Condición de estabilidad CFL**: $\Delta t < 2/\omega_{\max}$. Para barras: $\Delta t < L_{el}\sqrt{\rho/E}$. El solver detecta divergencia exponencial a posteriori y aborta con `RuntimeError` informativo si la condición se viola.

**Requisitos**: `lumping="lumped"` obligatorio (`"consistent"` rechazado con `ValueError` — el coste de invertir $\mathbf M$ consistente cada paso anula la ventaja del explícito). Frame3D con eje oblicuo a los ejes globales también rechazado (la rotación rompe la diagonalidad estricta — limitación documentada estándar, Cook-Malkus-Plesha §11.4).

```yaml
solver:
  type: CentralDifferenceSolver
  t_end: 2.5
  dt: 0.005                  # debe respetar CFL: Δt < L_el·√(ρ/E)
  lumping: lumped            # obligatorio
  nonlinear: false           # true para materiales con historia
  u0: [0.0, 0.0, 1.0, 0.0]
```

## `HarmonicSolver` — respuesta forzada armónica en frecuencia (ADR 0009 fase 6)

Análisis **lineal en estado estacionario**: dada una excitación armónica $\mathbf F(t) = \text{Re}\{\hat{\mathbf F}\,e^{i\omega t}\}$, calcula la amplitud compleja $\hat{\mathbf u}(\omega)$ de la respuesta estacionaria $\mathbf u(t) = \text{Re}\{\hat{\mathbf u}(\omega)\,e^{i\omega t}\}$ para un barrido de frecuencias.

**Ecuación resuelta** por frecuencia:

$$(-\omega^2\mathbf M + i\omega\mathbf C + \mathbf K)\,\hat{\mathbf u}(\omega) = \hat{\mathbf F}$$

Aritmética compleja directa con `scipy.sparse.linalg.spsolve`. Factorización LU compleja por frecuencia (no se cachea — la matriz cambia con $\omega$). Coste dominante: $O(n_\omega \cdot \text{coste\_factorizacion})$.

**Barrido configurable**:

- Vector explícito vía `omega` (lista de rad/s).
- Lineal: `omega_min`, `omega_max`, `n_omega`, `scale: linear`.
- Logarítmico: `scale: log` con `np.geomspace`.

**Cargas**: las `point_loads` del bloque YAML estándar se inyectan automáticamente como `F_amplitude` real con fase cero. Cargas con fase compleja requieren construcción desde código Python (YAML no soporta tipos complejos nativos).

**Salida**: `HarmonicResult` con `omega`, `u_complex` (shape $(n_{\text{dof}}, n_\omega)$), `F_complex`, y métodos `.amplitude()` ($|\hat{\mathbf u}|$) y `.phase()` ($\arg \hat{\mathbf u}$). Propiedad derivada `frequencies_hz`.

```yaml
point_loads:
  - {node_id: 2, ux: 1.0}    # amplitud real (fase 0)

solver:
  type: HarmonicSolver
  omega_min: 1.0
  omega_max: 10.0
  n_omega: 91
  scale: linear              # o "log" para FRFs sobre rangos amplios
  rayleigh:
    alpha: 0.5
    beta: 0.0
```

**Cuándo usarlo**: funciones de transferencia (FRFs), diagnóstico de resonancias antes de un transitorio caro, respuesta forzada por maquinaria rotante o cargas armónicas descompuestas.

**Caveats**: cerca de resonancia exacta sin amortiguamiento $\mathbf Z(\omega_n)$ es singular y el `spsolve` puede fallar — añadir Rayleigh.

## `ResponseSpectrumSolver` — análisis sísmico por combinación modal (ADR 0009 fase 7)

Análisis **sísmico clásico**: el suelo se caracteriza por un espectro de respuesta (en aceleración $S_a(T)$ o desplazamiento $S_d(T)$ vs período), no por un acelerograma temporal. El solver combina las respuestas máximas modales en una envolvente máxima — pierde la fase temporal entre modos por construcción, sólo conserva amplitudes envolventes.

**Esquema operativo**:

1. Análisis modal interno vía `ModalSolver(n_modes, sigma, lumping)`.
2. Factores de participación $\Gamma_n = \boldsymbol\phi_n^T\mathbf M\,\mathbf r$ con $\mathbf r$ el vector de excitación rígida (dirección sísmica).
3. Respuesta máxima por modo: $\mathbf u_n^{\max} = \Gamma_n\,\boldsymbol\phi_n\,S_d(\omega_n)$.
4. **Combinación SRSS** (Square Root of Sum of Squares, modos bien separados):

$$u^{\max}_k = \sqrt{\sum_n (u_{k,n}^{\max})^2}$$

5. **Combinación CQC** (Complete Quadratic Combination, Der Kiureghian 1980, para modos cercanos):

$$u^{\max}_k = \sqrt{\sum_i \sum_j \rho_{ij}\,u_{k,i}^{\max}\,u_{k,j}^{\max}}, \quad \rho_{ij} = \frac{8\xi^2(1+r) r^{3/2}}{(1-r^2)^2 + 4\xi^2 r(1+r)^2}$$

con $r = \omega_i/\omega_j$. Cuando $\xi \to 0$ o $r \to 0$: $\rho_{ij} \to \delta_{ij}$ y CQC degenera a SRSS.

**Interfaz del espectro**:

- Callable directo `spectrum(omega) -> S_d`.
- Dict tabulado: `{type: tabulated, kind: Sd|Sa, periods: [...], values: [...]}`. Interpolación lineal en período; clamp fuera de rango.
- Constante: `{type: constant_acceleration, Sa: ...}`.

**Dirección de excitación**:

- Vector explícito (`np.ndarray` shape $(n_{\text{dof}},)$).
- Atajo por DOF: `{dof_name: ux}` — el solver construye el vector unitario en ese DOF de todos los nodos.

```yaml
solver:
  type: ResponseSpectrumSolver
  n_modes: 5
  direction:
    dof_name: ux
  spectrum:
    type: tabulated
    kind: Sa                 # aceleración espectral; el solver convierte a Sd
    periods: [0.1, 0.5, 1.0, 2.0, 5.0]
    values:  [10.0, 10.0, 5.0, 1.0, 0.25]
  combination: SRSS          # o "CQC" para modos cercanos
  damping: 0.05              # ξ — sólo usado por CQC
```

**Salida**: `ResponseSpectrumResult` con `u_combined` (envolvente, positiva), `u_per_mode` (contribución modal con signo), `participation_factors` ($\Gamma_n$), `effective_masses` ($\Gamma_n^2$), `frequencies_rad`, `combination`, `damping`, `direction`. Método `.cumulative_effective_mass_ratio()` para verificar que la masa modal acumulada alcanza el objetivo normativo ($\geq 0.9$).

**Cuándo usarlo**: diseño sísmico contra espectros normativos (ASCE 7, Eurocódigo 8, NCh 433, NTC-DS México). Diagnóstico de qué modos dominan la respuesta en una dirección de excitación.

**Caveats**: precisión depende del número de modos incluidos (verifica `cumulative_effective_mass_ratio` ≥ 0.9). Excitación multi-direccional simultánea (CQC3, regla 100/30/30) requiere combinar varios análisis fuera del solver.

## Mass lumping (ADR 0009 fase 2)

Independiente de cualquier solver: todos los elementos del catálogo soportan `compute_mass_matrix(lumping="lumped")`. El esquema se elige por familia:

- **Sólidos isoparamétricos** (Tri3, Quad4, Tri6, Quad8, Quad9): **HRZ canónico** (Hinton-Rock-Zienkiewicz 1976) centralizado en `solidum.math.mass_lumping.lump_hrz`. Preserva masa total y proporcionalidad diagonal; evita masas negativas en elementos de orden superior.
- **Vigas y marcos** (Truss, Cable, Frame2D Euler/Timoshenko/EulerCorot, Frame3D): **lumping nodal directo** — $\rho AL/2$ traslacional + $\rho I L/2$ rotacional por nodo. Único esquema que sobrevive a la rotación local→global manteniendo $\mathbf M$ diagonal.

**Diagonalidad en globales**: garantizada en todos los elementos excepto Frame3D con eje oblicuo (queda **bloque-diagonal** por nodo — limitación documentada estándar; `CentralDifferenceSolver` rechaza este caso explícitamente).
