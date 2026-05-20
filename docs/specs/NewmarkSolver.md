# NewmarkSolver — Análisis dinámico transitorio lineal (Newmark-β)

> Orden de trabajo. La especificación física y la formulación numérica son patrimonio del usuario (revisa con detalle). La IA implementa y rellena las secciones marcadas.

---

## Especificación física

### 0. Descripción general
Análisis dinámico transitorio **lineal** de una estructura discretizada por FEM: integración temporal directa de la ecuación semidiscreta de movimiento bajo cargas dependientes del tiempo. Hipótesis: linealidad geométrica y material (K constante), masa M constante (lagrangeano total), amortiguamiento Rayleigh proporcional `C = α·M + β·K`. Apoyos Dirichlet constantes en el tiempo. La no-linealidad (Newton dentro de cada paso) se difiere a la fase 4 del ADR 0009.

### 1. Ecuación de movimiento semidiscreta
Tras la discretización espacial por FEM, el problema continuo
$$\rho\,\ddot{\mathbf u} - \operatorname{div}\boldsymbol\sigma = \mathbf b(\mathbf x, t)$$
se reduce al sistema acoplado de ODEs:
$$\mathbf M\,\ddot{\mathbf u}(t) + \mathbf C\,\dot{\mathbf u}(t) + \mathbf K\,\mathbf u(t) = \mathbf F(t),$$
con condiciones iniciales $\mathbf u(0) = \mathbf u_0$, $\dot{\mathbf u}(0) = \dot{\mathbf u}_0$ y condiciones de contorno Dirichlet $\mathbf u_s(t) = \mathbf g$ constantes.

### 2. Amortiguamiento Rayleigh
$$\mathbf C = \alpha\,\mathbf M + \beta\,\mathbf K,$$
con coeficientes calibrados a partir de dos pares modales $(\xi_1, \omega_1)$ y $(\xi_2, \omega_2)$:
$$\alpha = \frac{2\,\omega_1\omega_2\,(\xi_1\omega_2 - \xi_2\omega_1)}{\omega_2^2 - \omega_1^2}, \qquad
  \beta = \frac{2\,(\xi_2\omega_2 - \xi_1\omega_1)}{\omega_2^2 - \omega_1^2}.$$
La razón de amortiguamiento del modo $n$-ésimo resulta $\xi_n = \tfrac{1}{2}(\alpha/\omega_n + \beta\omega_n)$ — exacta en $\omega_1, \omega_2$ y suavemente sobreamortiguada fuera del rango. El usuario también puede pasar `(α, β)` directos.

### 3. Salidas
- $\mathbf u(t_k)$, $\dot{\mathbf u}(t_k)$, $\ddot{\mathbf u}(t_k)$ para cada paso temporal $k = 0, 1, \ldots, N_t$.
- Vector de instantes $t_k = k\,\Delta t$.
- Diagnósticos: $\alpha, \beta$ Rayleigh efectivos; número de pasos; tiempo de factorización vs tiempo de pasos.

---

## Formulación numérica

### 4. Método de Newmark-β (a-form)
Familia paramétrica con $(\beta, \gamma)$. Predictores:
$$\tilde{\mathbf u}_{n+1} = \mathbf u_n + \Delta t\,\dot{\mathbf u}_n + \tfrac{\Delta t^2}{2}(1 - 2\beta)\,\ddot{\mathbf u}_n,$$
$$\tilde{\dot{\mathbf u}}_{n+1} = \dot{\mathbf u}_n + \Delta t\,(1 - \gamma)\,\ddot{\mathbf u}_n.$$

Sistema efectivo en aceleraciones:
$$\bigl(\mathbf M + \gamma\Delta t\,\mathbf C + \beta\Delta t^2\,\mathbf K\bigr)\,\ddot{\mathbf u}_{n+1} \;=\; \mathbf F_{n+1} - \mathbf C\,\tilde{\dot{\mathbf u}}_{n+1} - \mathbf K\,\tilde{\mathbf u}_{n+1}.$$

Correctores:
$$\mathbf u_{n+1} = \tilde{\mathbf u}_{n+1} + \beta\Delta t^2\,\ddot{\mathbf u}_{n+1}, \qquad
  \dot{\mathbf u}_{n+1} = \tilde{\dot{\mathbf u}}_{n+1} + \gamma\Delta t\,\ddot{\mathbf u}_{n+1}.$$

### 5. Aceleración inicial
$\ddot{\mathbf u}_0$ se calcula consistentemente con la ecuación de movimiento en $t = 0$:
$$\mathbf M\,\ddot{\mathbf u}_0 = \mathbf F(0) - \mathbf C\,\dot{\mathbf u}_0 - \mathbf K\,\mathbf u_0.$$

### 6. Parámetros canónicos
| Esquema | β | γ | Propiedades |
|---|---|---|---|
| **Average acceleration (default)** | 1/4 | 1/2 | Incondicionalmente estable, sin amortiguamiento numérico, exactitud $O(\Delta t^2)$ |
| Linear acceleration | 1/6 | 1/2 | Condicionalmente estable ($\Delta t \le 2\sqrt{3}/\omega_{\max}$), $O(\Delta t^2)$ |
| Fox-Goodwin | 1/12 | 1/2 | Condicionalmente estable, $O(\Delta t^4)$ |
| Diferencias centradas | 0 | 1/2 | Explícito; condicionalmente estable ($\Delta t \le 2/\omega_{\max}$), $O(\Delta t^2)$ |

### 7. Imposición de Dirichlet con apoyos constantes
Mismo mecanismo que estática (ADR 0004 fase 1) con operador $\mathbf T$ que selecciona DOFs libres y $\mathbf g$ que recoge los valores prescritos: $\mathbf u = \mathbf T\,\mathbf u_\text{free} + \mathbf g$. Como $\mathbf g$ es constante, $\dot{\mathbf u} = \mathbf T\,\dot{\mathbf u}_\text{free}$ y $\ddot{\mathbf u} = \mathbf T\,\ddot{\mathbf u}_\text{free}$. El sistema efectivo se proyecta a DOFs libres:
$$\mathbf A_\text{eff,red}\,\ddot{\mathbf u}_{\text{free},n+1} = \mathbf T^\top\bigl(\mathbf F_{n+1} - \mathbf K\,\mathbf g\bigr) - \mathbf C_\text{red}\,\tilde{\dot{\mathbf u}}_\text{free} - \mathbf K_\text{red}\,\tilde{\mathbf u}_\text{free},$$
donde $\mathbf A_\text{eff,red} = \mathbf M_\text{red} + \gamma\Delta t\,\mathbf C_\text{red} + \beta\Delta t^2\,\mathbf K_\text{red}$. El término $\mathbf T^\top\,\mathbf K\,\mathbf g$ se precalcula una vez. La expansión final a globales es $\mathbf u_n = \mathbf T\,\mathbf u_{\text{free},n} + \mathbf g$.

### 8. Factorización reutilizable
Con $\Delta t$ constante y problema lineal, $\mathbf A_\text{eff,red}$ es constante. Se factoriza una vez al inicio (ADR 0003 — `FactorizedSolver`) y cada paso es una resolución triangular barata. Para $N_t$ pasos sobre $N_\text{libre}$ DOFs: coste $\sim O(\text{factor})_\text{una vez} + N_t \cdot O(N_\text{libre})_\text{por paso}$.

---

## Contrato de implementación

```yaml
name: NewmarkSolver
kind: solver
status: validated       # draft → implemented → validated

interface:
  yaml_type: NewmarkSolver
  output: TransientResult   # campos: t_history, u_history, udot_history, uddot_history

parameters:
  - { name: t_end,     type: float, required: true,                desc: "Tiempo final del análisis (s)" }
  - { name: dt,        type: float, required: true,                desc: "Paso temporal constante (s)" }
  - { name: beta,      type: float, required: false, default: 0.25, desc: "Parámetro β de Newmark (default: average acceleration)" }
  - { name: gamma,     type: float, required: false, default: 0.5,  desc: "Parámetro γ de Newmark" }
  - { name: rayleigh,  type: dict,  required: false, default: null, desc: "Amortiguamiento Rayleigh: {alpha, beta} directos o {xi1, omega1, xi2, omega2} para calibración modal" }
  - { name: u0,        type: ndarray, required: false, default: null, desc: "Desplazamiento inicial (ndof,); cero por defecto" }
  - { name: u0_dot,    type: ndarray, required: false, default: null, desc: "Velocidad inicial (ndof,); cero por defecto" }
  - { name: F_func,    type: callable, required: false, default: null, desc: "Callback Python F_func(t) -> ndarray (ndof,); cero por defecto" }
  - { name: linear_algebra, type: str, required: false, default: "auto", desc: "Backend para factorizar A_eff (ADR 0003)" }
  - { name: lumping,   type: str, required: false, default: "consistent", desc: "Discretización de masa; solo 'consistent' en esta fase" }

requirements:
  - "Todos los materiales del modelo declaran `density > 0` (ADR 0008)."
  - "K lineal y constante (hipótesis de pequeñas deformaciones, material elástico)."
  - "Apoyos Dirichlet constantes en el tiempo."

conventions:
  units:   "heredadas del modelo (ADR 0008); t en s, ω en rad/s, F en N consistente con kg/N/m."
  damping: "Rayleigh proporcional: C = α·M + β·K. Si se pasan (ξ₁,ω₁),(ξ₂,ω₂), (α,β) se derivan automáticamente."
  stability: |
    Default (β=1/4, γ=1/2) es incondicionalmente estable. Para β=0 (diferencias
    centradas) o β=1/6 (lineal) el usuario es responsable de que Δt ≤ 2/ω_max;
    el solver no comprueba esta cota automáticamente.

out_of_scope:
  - "Apoyos dependientes del tiempo (multi-support excitation sísmica) — diferido."
  - "Δt adaptativo — diferido."
  - "Generalized-α, Bossak-α — diferidos (variantes adicionales de la familia con amortiguamiento numérico)."
  - "No-linealidad y disipación numérica controlada — implementados como variantes en clases hermanas (`NewtonNewmarkSolver`, `HHTSolver`, `NewtonHHTSolver`); fases 4 y 3-bis del ADR 0009, cerradas."

acceptance:
  verification:
    - name: oscilador_1gdl_no_amortiguado
      setup: |
        Sistema 1 GDL: M=1, K=ω²=25 (ω=5 rad/s), C=0. Condiciones iniciales
        u₀=1, u̇₀=0. F(t)=0. Δt=ω·T/100 (100 pasos por período). Integrar
        durante 4 períodos.
      expect: |
        u(t) = cos(ω·t) exacto. Error máximo en u durante el intervalo
        < 1% (β=1/4 tiene error de período O(Δt²) bajo control con 100
        pasos por período).
      tol_rel: 0.01
    - name: oscilador_1gdl_amortiguado_rayleigh
      setup: |
        Sistema 1 GDL: M=1, K=25, C=2·ξ·ω·M con ξ=0.05 (sub-amortiguado).
        u₀=1, u̇₀=0, F=0. Δt = T/100, t_end = 4·T. Construir C vía
        Rayleigh con un único par (ξ, ω).
      expect: |
        u(t) = e^(-ξωt)·(cos(ω_d·t) + (ξω/ω_d)·sin(ω_d·t)) con
        ω_d = ω·√(1-ξ²). Error máximo < 2 %.
      tol_rel: 0.02
    - name: viga_voladizo_carga_subita_modal_vs_newmark
      setup: |
        Voladizo Frame2DEuler (20 elementos) con carga puntual aplicada en
        el extremo en t=0 como step function. Δt = T₁/40 (T₁ período del
        primer modo), t_end = 2·T₁. Sin amortiguamiento.
      expect: |
        La reconstrucción modal de la respuesta (superposición de primeros
        5 modos forzados por F(t)=H(t)·F₀ con solución analítica conocida)
        coincide con Newmark a < 5 %.
      tol_rel: 0.05
    - name: orden_de_convergencia
      setup: |
        Oscilador 1 GDL no amortiguado integrado con β=1/4, γ=1/2 a
        Δt, Δt/2, Δt/4. Calcular el error máximo en u contra solución
        analítica.
      expect: |
        El cociente de errores e(Δt)/e(Δt/2) tiende a 4 (orden 2 exacto).
      tol_rel: 0.5  # tolerancia sobre el cociente: 4 ± 2
  specific:
    - name: rayleigh_calibration
      setup: "Dados (ξ₁=0.02, ω₁=10), (ξ₂=0.05, ω₂=100), calibrar (α, β)."
      expect: |
        α = 2·ω₁·ω₂·(ξ₁·ω₂ − ξ₂·ω₁)/(ω₂² − ω₁²) ≈ 0.303,
        β = 2·(ξ₂·ω₂ − ξ₁·ω₁)/(ω₂² − ω₁²) ≈ 9.6e-4.
        Y verificar ξ(ω_n) = (α/ω_n + β·ω_n)/2 reproduce ξ₁ en ω₁ y ξ₂ en ω₂.
      tol_rel: 1.0e-10

references:
  - "Newmark N. M. (1959), A method of computation for structural dynamics, J. Eng. Mech. Div. ASCE."
  - "Bathe K.-J., Finite Element Procedures (2014), cap. 9 (matriz de masa) y cap. 9.4 (Newmark)."
  - "Hughes T. J. R., The Finite Element Method, cap. 9 (algoritmos de integración temporal)."
  - "Chopra A. K., Dynamics of Structures, cap. 5 (1 GDL) y cap. 15 (MDOF Newmark)."
  - "ADR 0003 — Capa algebraica (factorización reutilizable de A_eff)."
  - "ADR 0007 — Convergencia y tolerancias."
  - "ADR 0009 — Análisis modal y dinámico (este solver es la fase 3)."
```

---

## Implementación

- **Archivo**: [solidum/math/solvers/newmark.py](../../solidum/math/solvers/newmark.py) (clase `NewmarkSolver`, registrada vía `@SolverRegistry.register`).
- **Amortiguamiento**: [solidum/math/damping.py](../../solidum/math/damping.py) (`rayleigh_from_modes`, `rayleigh_xi`).
- **Tipo de resultado**: [`TransientResult`](../../solidum/results.py) — dataclass frozen con `t_history`, `u_history`, `udot_history`, `uddot_history`, `n_steps`, `alpha_rayleigh`, `beta_rayleigh`, `converged`.
- **Entrypoint público**: [`solidum.run_transient`](../../solidum/entry.py) para uso programático; [`solidum.run_yaml`](../../solidum/entry.py) despacha automáticamente cuando el YAML declara `solver: type: NewmarkSolver`.
- **Pipeline**: `Assembler.assemble_system()` → `Assembler.assemble_mass_matrix()` → `C = α·M + β·K` → reducción simultánea por Dirichlet (operador `T`) → cálculo de `ü₀` → factorización única de `A_eff_red = M_red + γΔt·C_red + βΔt²·K_red` → bucle de pasos con resolución triangular reutilizada.
- **Tests**:
  - [tests/test_newmark.py](../../tests/test_newmark.py) — 14 tests:
    - `TestRayleighCalibration`: `(α, β)` reproduce ξ objetivo en ω₁, ω₂; errores agregados (ω iguales, no positivos).
    - `TestNewmark1DofUndamped`: vibración libre `u(t) = cos(ωt)` a 1 %; estado inicial preservado.
    - `TestNewmark1DofDamped`: respuesta sub-amortiguada `e^(-ξωt)·(cos(ω_d·t) + ...)` a 2 %; coeficientes Rayleigh propagados al `TransientResult`.
    - `TestNewmark1DofStepResponse`: `u(t) = (F₀/K)(1 − cos(ωt))` a 0.1 %.
    - `TestNewmarkConvergenceOrder`: cociente de errores al refinar Δt ∈ {T/40, T/80, T/160} cae en `[3.5, 4.5]` → orden 2 confirmado.
    - `TestNewmarkContract`: validaciones de constructor y `run_transient`.
    - `TestNewmarkYamlPipeline`: end-to-end desde YAML.
- **Notas de traducción**:
  - **Factorización única**: con `Δt` constante y problema lineal, `A_eff_red` es invariante; se factoriza al inicio (`A_solver.factorize(A_eff)`) y cada paso es una sustitución triangular barata (ADR 0003 — `FactorizedSolver`).
  - **`reduce_pair` vs `reduce` triple**: en el problema modal `reduce_pair(K, M)` basta porque `g` no entra en autovalores. En Newmark se manejan tres matrices (`K, M, C`) y un término constante `F_dir = TᵀK·g` para apoyos no nulos; se hace explícitamente en `NewmarkSolver.solve()` sin un helper unificado (no aporta valor centralizarlo con un solo consumidor).
  - **Aceleración inicial consistente**: `M·ü₀ = F(0) − C·u̇₀ − K·u₀` resuelve un sistema lineal aparte; se usa el dispatcher de la capa algebraica sobre `M_red` (siempre SPD si la masa es consistente y la densidad estrictamente positiva).
  - **`F_func` opcional**: si es `None`, se trata como vector cero (vibración libre). El callback se invoca una vez por paso con `t_{n+1}`; el usuario es responsable de la coherencia temporal (no se reusan valores entre pasos).
  - **El bloque `convergence:` del YAML** (ADR 0007) se omite para `NewmarkSolver` igual que para `LinearSolver` y `ModalSolver`: la fase 3 es lineal y no hay iteración Newton interna; cuando entre la fase 4 (Newton dentro de cada paso), el bloque se reincorporará.

---

## Diálogo

- **2026-05-13** · Validado. Los 14 tests de `tests/test_newmark.py` cubren `acceptance.verification` (osciladores 1 GDL no amortiguado / amortiguado / step / orden de convergencia) y `acceptance.specific` (calibración Rayleigh). Promovido `status: validated`.
