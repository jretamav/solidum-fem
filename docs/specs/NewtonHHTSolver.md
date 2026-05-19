# NewtonHHTSolver — Análisis dinámico transitorio no lineal por HHT-α + Newton

> Variante de [`NewtonNewmarkSolver`](NewtonNewmarkSolver.md) según la regla
> de variantes (Reglas.md §4). Aplica el esquema temporal HHT-α (Hilber-
> Hughes-Taylor 1977) al residuo dinámico no lineal: amortiguamiento numérico
> controlado en altas frecuencias dentro del bucle Newton de cada paso.
>
> Reusa todas las decisiones arquitecturales del ADR 0009 (esquema temporal,
> contrato `PIPELINE_KIND`, telemetría tipada, line search opt-in) y de
> [`HHTSolver`](HHTSolver.md) (formulación temporal, autoderivación de β/γ).
> **Sin ADR nuevo**. Documenta solo lo que cambia respecto a los dos padres.

---

## Qué cambia respecto a `NewtonNewmarkSolver`

`NewtonNewmarkSolver` integra `M·ü + C·u̇ + F_int(u) = F(t)` con Newmark-β
clásico (β=1/4, γ=1/2 default, sin amortiguamiento numérico). Cada paso
abre un bucle Newton porque `F_int(u)` es no lineal en `u` (plasticidad,
daño, corotacional, embedded).

`NewtonHHTSolver` reemplaza el residuo Newmark por el residuo HHT-α: las
**fuerzas internas, viscosas y externas** se evalúan en un instante
intermedio `t_{n+1-α}` mientras la **fuerza inercial** se mantiene en
`t_{n+1}`. El bucle Newton se ejecuta sobre este residuo modificado.

### Residuo dinámico HHT-α no lineal

$$\mathbf R(\ddot{\mathbf u}_{n+1}) \;=\; (1+\alpha)\,\mathbf F_{n+1} - \alpha\,\mathbf F_n
   - \big[(1+\alpha)\,\mathbf F_{\text{int}}(\mathbf u_{n+1}) - \alpha\,\mathbf F_{\text{int}}(\mathbf u_n)\big]
   - \big[(1+\alpha)\,\mathbf C\,\dot{\mathbf u}_{n+1} - \alpha\,\mathbf C\,\dot{\mathbf u}_n\big]
   - \mathbf M\,\ddot{\mathbf u}_{n+1}$$

El término no lineal $\mathbf F_{\text{int}}(\mathbf u_{n+1})$ depende del
estado actualizado por los correctores Newmark, así que el residuo es
implícito en $\ddot{\mathbf u}_{n+1}$ y requiere Newton.

### Jacobiano efectivo

$$\mathbf J \;=\; \mathbf M + (1+\alpha)\,\gamma\,\Delta t\,\mathbf C
                            + (1+\alpha)\,\beta\,\Delta t^2\,\mathbf K_t$$

donde $\mathbf K_t = \partial\mathbf F_{\text{int}}/\partial\mathbf u$ es
la tangente material/geométrica reensamblada en cada iteración Newton.

### Parámetros heredados de `HHTSolver`

- Rango admisible $\alpha \in [-1/3, 0]$. Valores fuera lanzan `ValueError`.
- $\beta$ y $\gamma$ auto-derivados desde $\alpha$:

$$\beta = \frac{(1-\alpha)^2}{4}, \qquad \gamma = \frac{1-2\alpha}{2}$$

- Override explícito de $\beta$/$\gamma$ permitido pero genera `RuntimeWarning`
  (combinaciones distintas pueden romper estabilidad incondicional u orden 2).
- Radio espectral en alta frecuencia: $\rho_\infty = (1+\alpha)/(1-\alpha)$.
  Disipación selectiva en modos con $\omega\Delta t > \pi$.

### Comportamiento que NO cambia

Todo lo heredado de `NewtonNewmarkSolver`:

- Convergencia dual fuerza + desplazamiento (ADR 0007) con autocalibración.
- Line search GLL opt-in (ADR 0011, default `False`).
- Diagnóstico tipado de divergencia (`SolverDivergedError` y subclases).
- Rayleigh con $\mathbf K_0$ constante (matriz de rigidez inicial).
- Commit/rollback del estado material por paso.
- Detección de oscilación por ventana móvil.
- Despacho algebraico ADR 0003 con fallback Cholesky→LU.
- Lumping de masa (consistent / lumped) ADR 0009 fase 2.

### Recuperación de `NewtonNewmark` con $\alpha = 0$

Con $\alpha = 0$: $\beta = 1/4$, $\gamma = 1/2$, el residuo HHT colapsa al
de Newmark estándar. Verificado en tests:
`tests/test_dynamic_nonlinear_hht.py::test_alpha_zero_recovers_newton_newmark`.

### Recuperación de `HHTSolver` con materiales lineales

Con materiales lineales ($\mathbf F_{\text{int}}(\mathbf u) = \mathbf K\,\mathbf u$),
el residuo de Newton se anula en una iteración y el solver reproduce
exactamente `HHTSolver`. Tangente constante = rigidez secante. Verificado:
`tests/test_dynamic_nonlinear_hht.py::test_linear_material_matches_hht_solver`.

---

## Contrato de implementación

```yaml
name: NewtonHHTSolver
kind: solver
status: validated

extends: NewtonNewmarkSolver       # variante por Reglas §4
also_references: HHTSolver         # formulación temporal heredada

interface:
  yaml_type: newton_hht
  output: TransientResult
  pipeline_kind: transient

parameters:
  # Solo se listan los parámetros propios o con default distinto al padre.
  # Resto heredados de NewtonNewmarkSolver.
  - { name: alpha, type: float, required: false, default: -0.05,
      desc: "Parámetro HHT-α ∈ [-1/3, 0]. α=0 recupera NewtonNewmark; α=-1/3 disipación máxima preservando orden 2" }
  - { name: beta,  type: float, required: false, default: null,
      desc: "Override del β autoderivado (no recomendado — emite RuntimeWarning)" }
  - { name: gamma, type: float, required: false, default: null,
      desc: "Override del γ autoderivado (no recomendado — emite RuntimeWarning)" }
  # heredados de NewtonNewmarkSolver:
  #   t_end, dt, convergence, max_iter, freeze_tangent_after_iter,
  #   line_search, rayleigh, u0, u0_dot, F_func, linear_algebra, lumping

requirements:
  - "Modelo dinámico con no-linealidad de material o geometría"
  - "M con density > 0 en todos los materiales activos (ADR 0008)"
  - "Δt razonable (incondicionalmente estable para sistemas lineales, no garantizado para fuertes no-linealidades)"

conventions:
  units: "heredadas del modelo (ADR 0008)"
  stability: "incondicionalmente estable en problemas lineales para α ∈ [-1/3, 0]; condicionalmente para no lineales severos"
  disipacion: "ρ_∞ = (1+α)/(1-α); selectiva en alta frecuencia"

out_of_scope:
  - "α fuera de [-1/3, 0] (pierde estabilidad u orden 2)"
  - "Generalized-α / Bossak (variantes con disipación distinta) — diferidos"
  - "Análisis estático ⇒ usar NonlinearSolver"

acceptance:
  verification:
    - name: alpha_cero_recupera_newton_newmark
      setup: "mismo problema con NewtonNewmarkSolver y NewtonHHTSolver(α=0)"
      expect: "respuesta idéntica paso a paso"
      tol_rel: 1.0e-12

    - name: material_lineal_recupera_hht
      setup: "mismo problema con HHTSolver y NewtonHHTSolver, ambos α=−0.1"
      expect: "respuesta idéntica; Newton converge en 1 iteración por paso"
      tol_rel: 1.0e-12

    - name: disipacion_modo_alta_frecuencia
      setup: "respuesta libre con modos altos espurios excitados; α=−0.1"
      expect: "amplitud de modos espurios decae ~ρ_∞ por periodo"
      tol_rel: 0.05

    - name: respuesta_plastica_dinamica
      setup: "barra Truss2D con Elastoplastic1D, carga dinámica que cruza el yield"
      expect: "estado plástico converged paso a paso; histéresis dinámica coherente"
      tol_rel: 1.0e-3

    - name: rechazo_alpha_fuera_de_rango
      setup: "NewtonHHTSolver(alpha=-0.5) o (alpha=0.1)"
      expect: "ValueError con mensaje claro"

    - name: warning_override_beta_gamma
      setup: "NewtonHHTSolver(alpha=-0.1, beta=0.3)"
      expect: "RuntimeWarning emitido al construir"

  specific:
    - name: tangente_singular_dispara_excepcion_tipada
      setup: "modelo con material que produce K_t singular en algún paso"
      expect: "raise SingularTangentError con métricas estructuradas"

references:
  - "Hilber H.M., Hughes T.J.R., Taylor R.L. (1977). Improved numerical dissipation for time integration algorithms in structural dynamics. Earthquake Eng. Struct. Dyn. 5, 283-292."
  - "Hughes T.J.R. (2000). The Finite Element Method. Dover. §9.3."
  - "Newmark N.M. (1959). A method of computation for structural dynamics. ASCE 85(EM3), 67-94."
  - "ADR 0009 (fase 4: HHT-α — formulación temporal en residuo no lineal)."
  - "ADR 0011 (excepciones tipadas, line search opt-in)."
  - "Spec padre: HHTSolver.md (formulación temporal lineal)."
  - "Spec padre: NewtonNewmarkSolver.md (Newton dentro de Newmark no lineal)."
```

---

## Implementación

- **Archivo**: [fenix/math/solvers/newmark.py](../../fenix/math/solvers/newmark.py) (clase `NewtonHHTSolver`).
- **Clase**: hereda de `NewtonNewmarkSolver`; añade el atributo `self.alpha` y reemplaza el cómputo del residuo y del jacobiano en `solve()`.
- **Auto-derivación** de β y γ vía `_hht_autoderive_beta_gamma(alpha)` (función helper compartida con `HHTSolver`).
- **Validación temprana** en `__init__`: `alpha ∈ [-1/3, 0]`, warning si override de β/γ.
- **Entrypoint público**: `fenix.run_transient(model, solver="newton_hht", ...)` (despacho por `PIPELINE_KIND = "transient"`).
- **Tests**:
  - [tests/test_dynamic_nonlinear_hht.py](../../tests/test_dynamic_nonlinear_hht.py) — incluye recuperación de `NewtonNewmark` con α=0, recuperación de `HHTSolver` con materiales lineales, respuesta plástica dinámica.

---

## Diálogo

- **2026-05-19** · Spec corta tipo extensión creada para cerrar el hueco H-5.3 (variante de `HHTSolver`/`NewtonNewmarkSolver` por Reglas §4). Sin ADR nuevo: reúsa ADR 0009 (formulación temporal HHT) y ADR 0011 (Newton + diagnóstico). El comportamiento no se modifica.
