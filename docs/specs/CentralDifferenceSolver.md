# CentralDifferenceSolver — integración explícita por diferencias centradas (ADR 0009 fase 5)

> Spec de **solver nuevo** (no variante). Cubre el integrador explícito de segundo orden para análisis dinámico transitorio lineal y no lineal. Aprovecha la masa lumped diagonal habilitada por la fase 2 del ADR 0009 para evitar resolver sistema lineal en cada paso.

---

## Especificación física

### 0. Descripción general

Análisis dinámico transitorio **explícito** de segundo orden. Resuelve `M·ü + C·u̇ + K·u = F(t)` (lineal) o `M·ü + C·u̇ + F_int(u) = F(t)` (no lineal) avanzando la aceleración en cada paso por inversión trivial de la masa diagonal lumped — sin sistema lineal, sin Newton interno. Apropiado para wave propagation, impacto, dinámica de respuesta rápida.

Hipótesis cinemáticas globales: las mismas que admite el bloque de elementos del modelo (corotacional, pequeñas deformaciones, etc.). Material lineal o no lineal con historia — en modo `nonlinear=True` el solver llama `Assembler.assemble_non_linear_system(u_n)` cada paso para reevaluar `F_int`.

### 1. Ecuación de equilibrio resuelta

$$\mathbf{M}\,\ddot{\mathbf{u}} + \mathbf{C}\,\dot{\mathbf{u}} + \mathbf{F}_\text{int}(\mathbf{u}) = \mathbf{F}_\text{ext}(t)$$

con $\mathbf{F}_\text{int}(\mathbf{u}) = \mathbf{K}\,\mathbf{u}$ en el caso lineal.

### 2. Condiciones iniciales / de contorno

- Apoyos Dirichlet constantes en el tiempo (eliminación directa, ADR 0004).
- $\mathbf{u}(0) = \mathbf{u}_0$, $\dot{\mathbf{u}}(0) = \dot{\mathbf{u}}_0$ (cero por defecto).
- Cargas externas $\mathbf{F}_\text{ext}(t)$ vía `F_func(t)`.
- Amortiguamiento Rayleigh $\mathbf{C} = \alpha\,\mathbf{M} + \beta\,\mathbf{K}$ (mismo contrato que `NewmarkSolver`).

MPC no soportado en esta fase — no aparece en uso real con dinámica explícita.

### 3. Salidas físicas

`TransientResult` con historiales:

- `u_history[ndof, n_steps+1]` desplazamientos.
- `udot_history[ndof, n_steps+1]` velocidades.
- `uddot_history[ndof, n_steps+1]` aceleraciones.
- `t_history[n_steps+1]` instantes de tiempo.

---

## Formulación numérica

### 4. Esquema operativo — leapfrog Belytschko-Liu-Moran

Inicialización:

$$\mathbf{a}_0 = \mathbf{M}^{-1}\big(\mathbf{F}(0) - \mathbf{F}_\text{int}(\mathbf{u}_0) - \mathbf{C}\,\mathbf{v}_0 - \mathbf{F}_\text{dir}\big)$$

$$\mathbf{v}_{1/2} = \mathbf{v}_0 + \tfrac{\Delta t}{2}\,\mathbf{a}_0$$

Paso $n \to n+1$:

$$\mathbf{u}_{n+1} = \mathbf{u}_n + \Delta t\,\mathbf{v}_{n+1/2}$$

$$\mathbf{F}_\text{int}^{n+1} = \mathbf{K}\,\mathbf{u}_{n+1} \quad\text{(lineal)} \quad\text{o}\quad \text{ensamble}(\mathbf{u}_{n+1}) \quad\text{(no lineal)}$$

$$\mathbf{a}_{n+1} = \mathbf{M}^{-1}\big(\mathbf{F}(t_{n+1}) - \mathbf{F}_\text{int}^{n+1} - \mathbf{C}\,\mathbf{v}_{n+1/2} - \mathbf{F}_\text{dir}\big)$$

$$\mathbf{v}_{n+1} = \mathbf{v}_{n+1/2} + \tfrac{\Delta t}{2}\,\mathbf{a}_{n+1} \quad\text{(salida de paso)}$$

$$\mathbf{v}_{n+3/2} = \mathbf{v}_{n+1} + \tfrac{\Delta t}{2}\,\mathbf{a}_{n+1} \quad\text{(para el paso siguiente)}$$

Sin Newton interno — la ecuación de aceleración es explícita en $\mathbf{F}_\text{int}^{n+1}$ que ya se evaluó en el estado actual.

### 5. Predictor / corrector

El esquema leapfrog **no tiene predictor/corrector** en el sentido implícito. La actualización $\mathbf{u}_{n+1}$ se calcula sólo con datos en $n$ y $n+1/2$; la velocidad $\mathbf{v}_{n+1}$ se obtiene del semipaso.

### 6. Criterio de convergencia

No aplica — esquema no iterativo. La convergencia "del problema" se materializa en la condición CFL y en la detección de divergencia.

### 7. Imposición de Dirichlet y MPC

Eliminación directa idéntica a `NewmarkSolver` (operador $\mathbf{T}$ + offset $\mathbf{g}$). MPC fuera de alcance en esta fase.

### 8. Backend algebraico

**Ninguno** — la única operación lineal es $\mathbf{M}^{-1} \cdot \mathbf{r}$ con $\mathbf{M}$ diagonal, implementada como $1/\text{diag}(\mathbf{M})$ y multiplicación elemento a elemento. El solver verifica al construir $\mathbf{M}_\text{red}$ que la diagonal es estrictamente positiva y que los off-diagonals son numéricamente nulos (tolerancia $10^{-12} \cdot \max|\text{diag}|$); si no, lanza `ValueError`.

### 9. Adaptatividad / control de paso

No adaptativo en esta fase — paso temporal constante $\Delta t$. La adaptatividad explícita (e.g. estimación de $\omega_\max$ por power iteration con actualización por paso) queda diferida hasta que aparezca un caso de uso con $\omega_\max$ variable significativamente en el tiempo.

### 10. Caveats numéricos

- **Estabilidad condicional**. La condición CFL exige $\Delta t < 2/\omega_\max$ donde $\omega_\max$ es la frecuencia natural máxima del sistema discreto. Para una barra axial elástica con elemento de longitud $L_e$ y velocidad de propagación $c_p = \sqrt{E/\rho}$: $\Delta t_\text{crit} \approx L_e/c_p$. Para sólidos 2D la fórmula análoga es $\Delta t_\text{crit} \approx h_e/c_d$ donde $c_d$ es la velocidad dilatacional. Si se excede, la solución diverge exponencialmente — el solver lo detecta a posteriori (umbral configurable `divergence_threshold`, default $10^8 \cdot |\mathbf{u}_0|$) y aborta con `RuntimeError` informativo.

- **Frame3D oblicuo**. La fase 2 del ADR 0009 documenta que $\mathbf{M}_\text{lumped}$ es bloque-diagonal (no estrictamente diagonal) en Frame3D cuando el eje del elemento no coincide con un eje global. Para central differences esto significaría invertir bloques 6×6 por nodo cada paso, **anulando** la ventaja del esquema explícito. El solver rechaza ese caso con `ValueError`.

- **Lumping obligatorio**. `lumping="consistent"` se rechaza al construir el solver — invertir una $\mathbf{M}$ consistente cada paso es caro y elimina la razón de ser del explícito. El usuario debe pasar `lumping="lumped"` (default).

- **Modos espurios de altas frecuencias**. El esquema no introduce amortiguamiento numérico (a diferencia de HHT-α); las altas frecuencias de la malla pueden contaminar la respuesta de bajas frecuencias en problemas de wave propagation. Mitigaciones: refinamiento de malla, lumping puro (mejor filtrado que consistent), o cambio a HHT-α implícito.

---

## Contrato de implementación

```yaml
name: CentralDifferenceSolver
kind: solver
status: validated

interface:
  yaml_type: CentralDifferenceSolver
  output: TransientResult

parameters:
  - { name: t_end, type: float, required: true, desc: "tiempo final (s)" }
  - { name: dt, type: float, required: true, desc: "paso temporal constante (s); debe satisfacer Δt < Δt_crit" }
  - { name: nonlinear, type: bool, required: false, desc: "default False; True recalcula F_int cada paso vía ensamble" }
  - { name: rayleigh, type: dict | null, required: false, desc: "C = α·M + β·K; mismo contrato que NewmarkSolver" }
  - { name: u0, type: ndarray | null, required: false, desc: "desplazamiento inicial global" }
  - { name: u0_dot, type: ndarray | null, required: false, desc: "velocidad inicial global" }
  - { name: F_func, type: callable | null, required: false, desc: "F(t) → ndarray global" }
  - { name: lumping, type: str, required: false, desc: "default 'lumped'; 'consistent' rechazado" }
  - { name: divergence_threshold, type: float, required: false, desc: "default 1e8; multiplicador sobre escala inicial para abortar" }

requirements:
  - "Materiales con `density` declarada (ADR 0008)."
  - "Todos los elementos del modelo deben soportar `compute_mass_matrix(lumping='lumped')` (ADR 0009 fase 2)."
  - "Frame3D con eje oblicuo a globales NO soportado (M_red no es diagonal pura)."

conventions:
  units:    "heredadas del modelo (ADR 0008)."
  stability: "Condicional — Δt < 2/ω_max (CFL). Para barras: Δt < L_el/c_p. Sin power iteration en esta fase — responsabilidad del usuario calcular o estimar Δt_crit."

out_of_scope:
  - "Estimación automática de Δt_crit (power iteration sobre M⁻¹K) — diferida."
  - "Paso adaptativo (Δt variable) — diferido."
  - "MPC (restricciones lineales multipunto) — diferido."
  - "Frame3D con eje oblicuo (rechazado por arquitectura del lumping)."

acceptance:
  verification:
    - name: free_vibration_1dof_undamped
      setup: "Truss2D 1 elemento con E=25, ρ=2, A=1, L=1 (lumped). u₀=1, v₀=0. Δt = T/200."
      expect: "u(t) = cos(ωt) con ω=5 rad/s, error < 2%"
      tol_rel: 0.02
    - name: step_response_1dof
      setup: "Mismo modelo con F₀=50·H(t)"
      expect: "u(t) = (F₀/K)(1 − cos(ωt)), error < 2% del régimen permanente"
      tol_rel: 0.02
    - name: matches_newmark_at_small_dt
      setup: "Mismo oscilador, Newmark β=1/4 γ=1/2 con masa lumped vs CentralDiff con Δt = T/500"
      expect: "u_Newmark − u_CD máximo < 1%"
      tol_rel: 0.01
  specific:
    - name: cfl_violation_diverges
      setup: "Δt = 3·Δt_crit"
      expect: "RuntimeError con mensaje CFL informativo"
      tol_rel: 0
    - name: cfl_threshold_matches_2_over_omega_max_1dof
      setup: "Oscilador 1-GDL (ω=5, Δt_crit=2/ω=0.4): integrar con Δt = 0.95·Δt_crit vs 1.05·Δt_crit"
      expect: "0.95·Δt_crit produce trayectoria acotada (30 ciclos); 1.05·Δt_crit lanza RuntimeError CFL"
      tol_rel: 0
    - name: cfl_threshold_matches_2_over_omega_max_2dof
      setup: |
        Cadena 2-DOF (3 trusses E=A=L=1, ρ=1, lumped) con K=[[2,-1],[-1,2]], M=I →
        autovalores cerrados ω² ∈ {1, 3}, ω_max=√3, Δt_crit=2/√3. Estado inicial
        excita modo antisimétrico (1, -1) — puro ω_max.
      expect: |
        (1) Δt = 0.95·Δt_crit estable.
        (2) Δt = 1.05·Δt_crit lanza RuntimeError CFL.
        (3) Δt = Δt_crit/10 reproduce u(t) = (cos(√3·t), -cos(√3·t)) con error < 2%.
      tol_rel: 0.02
      ref: "tests/test_central_difference.py::TestCentralDifferenceCFLAnalytic"
    - name: linear_mode_matches_nonlinear_mode_for_linear_material
      setup: "Mismo problema con material lineal puro, solver con nonlinear=False vs nonlinear=True"
      expect: "Coincidencia exacta a precisión de máquina"
      tol_rel: 1.0e-9
    - name: consistent_lumping_rejected
      setup: "Construir solver con lumping='consistent'"
      expect: "ValueError"
      tol_rel: 0
    - name: frame3d_oblique_rejected
      setup: "Frame3D con eje en dirección (1,1,1)/√3 e Iy ≠ Iz"
      expect: "ValueError 'no es estrictamente diagonal'"
      tol_rel: 0

references:
  - "Belytschko, T., Liu, W. K., Moran, B. (2014). *Nonlinear Finite Elements for Continua and Structures*, §6.2."
  - "Hughes, T. J. R. (2000). *The Finite Element Method*, §9.1."
  - "Crisfield, M. A. (1997). *Non-Linear Finite Element Analysis of Solids and Structures*, vol. 2, §24.6."
  - "Bathe, K.-J. (2014). *Finite Element Procedures*, §9.2 (estabilidad y CFL)."
```

---

## Implementación

- Archivo: [fenix/math/solvers/central_difference.py](../../fenix/math/solvers/central_difference.py)
- Clase: `CentralDifferenceSolver` (registrada en `SolverRegistry`)
- Atributo de despacho: `PIPELINE_KIND = "transient"` — `run_yaml` despacha a `run_transient` por el mecanismo declarativo de la regla C (aplicada 2026-05-18 al introducir este solver, primer caso fuera de la jerarquía de Newmark)
- Tests: [tests/test_central_difference.py](../../tests/test_central_difference.py) (10 tests verdes)

Notas de traducción: ninguna — la convención de signos del proyecto (Reglas §5) ya es la stress-resultant en el interior, sin reinterpretación específica de central differences.

---

## Diálogo

- *2026-05-18*: regla C aplicada al introducir este solver. Antes de él, `run_yaml::isinstance(NewmarkSolver)` capturaba HHT-α por subclasificación; CentralDifference no podía heredar de NewmarkSolver (flujo radicalmente distinto, sin sistema lineal en cada paso, estabilidad condicional). En lugar de añadir una tercera rama `isinstance`, se refactorizó el dispatch a `PIPELINE_KIND` declarativo — futuros solvers no clásicos (Harmonic, ResponseSpectrum, …) no requieren tocar `entry.py`.
