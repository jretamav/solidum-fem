# HHTSolver — Análisis dinámico transitorio lineal por Hilber-Hughes-Taylor α

> Variante de [`NewmarkSolver`](NewmarkSolver.md) según la regla de variantes
> (Reglas.md §4). Documenta **solo lo que cambia** respecto al padre. Implementa
> el integrador HHT-α (1977) para introducir **amortiguamiento numérico
> controlado** en altas frecuencias preservando segundo orden de precisión y
> estabilidad incondicional. Reusa todas las decisiones arquitecturales del
> ADR 0009 — sin ADR nuevo.

---

## Qué cambia respecto a `NewmarkSolver`

`NewmarkSolver` integra `M·ü + C·u̇ + K·u = F(t)` evaluando todas las fuerzas
en el instante final `t_{n+1}` del paso. Con la elección β=1/4, γ=1/2 (regla
trapezoidal, default) el integrador es **no disipativo**: los modos de alta
frecuencia espurios introducidos por la discretización persisten
indefinidamente en la respuesta y pueden enmascarar la señal de interés.

HHT-α modifica la ecuación de equilibrio temporal para evaluar las **fuerzas
viscosas, elásticas y externas** en un instante intermedio `t_{n+1-α}`,
mientras la **fuerza inercial** se mantiene en `t_{n+1}`. Esto introduce
disipación numérica controlada por un único parámetro escalar `α ∈ [−1/3, 0]`,
sin sacrificar la precisión de segundo orden global ni la estabilidad
incondicional. Con α = 0 se recupera Newmark trapezoidal.

### Ecuación de movimiento HHT-α

$$\mathbf M\,\ddot{\mathbf u}_{n+1}
  + (1+\alpha)\,\mathbf C\,\dot{\mathbf u}_{n+1} - \alpha\,\mathbf C\,\dot{\mathbf u}_n
  + (1+\alpha)\,\mathbf K\,\mathbf u_{n+1} - \alpha\,\mathbf K\,\mathbf u_n
  \;=\; (1+\alpha)\,\mathbf F(t_{n+1}) - \alpha\,\mathbf F(t_n)$$

Para preservar el orden 2 y la estabilidad incondicional, los parámetros
Newmark se derivan automáticamente del α:

$$\beta = \frac{(1-\alpha)^2}{4}, \qquad \gamma = \frac{1-2\alpha}{2}$$

El usuario pasa solo `alpha`; β y γ se autoderivan salvo override explícito.

### Sistema efectivo

Sustituyendo los predictores y correctores Newmark estándar (idénticos a la
fase 3 con los β, γ autoderivados), el sistema efectivo a resolver en cada
paso es:

$$\mathbf A_\text{eff} \,\ddot{\mathbf u}_{n+1} \;=\; \mathbf b$$

con

$$\mathbf A_\text{eff} \;=\; \mathbf M + (1+\alpha)\,\gamma\Delta t\,\mathbf C + (1+\alpha)\,\beta\Delta t^2\,\mathbf K$$

$$\mathbf b \;=\; (1+\alpha)\,\mathbf F_{n+1} - \alpha\,\mathbf F_n
              - \mathbf F_\text{dir}
              - (1+\alpha)\,\mathbf C\,\tilde{\dot{\mathbf u}} + \alpha\,\mathbf C\,\dot{\mathbf u}_n
              - (1+\alpha)\,\mathbf K\,\tilde{\mathbf u} + \alpha\,\mathbf K\,\mathbf u_n$$

`A_eff` es constante (depende solo de M, C, K y de los coeficientes), por lo
que se **factoriza una sola vez** al inicio del análisis (igual que Newmark).
Cada paso temporal es una resolución triangular barata.

### Disipación numérica controlada

El radio espectral del operador de amplificación en alta frecuencia
(`Δt → ∞`) es

$$\rho_\infty = \frac{1+\alpha}{1-\alpha}$$

Para `α = 0`: `ρ_∞ = 1` (sin disipación, persistencia indefinida del modo).
Para `α = −0.05`: `ρ_∞ ≈ 0.905` (disipación ligera, modo decae 9.5% por
periodo de alta frecuencia). Para `α = −1/3` (límite teórico): `ρ_∞ = 0.5`
(disipación máxima manteniendo orden 2).

La curva de disipación es **selectiva**: amortigua significativamente solo los
modos con `ωΔt > π`, dejando casi intactos los modos de interés con
`ωΔt « 1`. Es lo que distingue a HHT-α de simplemente subir el γ Newmark
(que disipa en *todas* las frecuencias por igual y baja el orden a 1).

### Recuperación del comportamiento Newmark

Con `α = 0`, los coeficientes auto-derivados quedan `β = 1/4, γ = 1/2` y la
ecuación de movimiento HHT colapsa exactamente a la de Newmark trapezoidal.
`HHTSolver(alpha=0.0)` reproduce **bit a bit** los resultados de
`NewmarkSolver` con `beta=0.25, gamma=0.5`. Esto se valida en los acceptance.

### Variante no lineal

`NewtonHHTSolver` (subclase de `NewtonNewmarkSolver`) aplica el mismo cambio
de la ecuación de movimiento al residuo dinámico. El Newton interno y su
jacobiano se ven afectados análogamente:

$$\mathbf R(\ddot{\mathbf u}_{n+1}) \;=\; (1+\alpha)\,\mathbf F_{n+1} - \alpha\,\mathbf F_n
  - [(1+\alpha)\,\mathbf F_\text{int}(\mathbf u_{n+1}) - \alpha\,\mathbf F_\text{int}(\mathbf u_n)]
  - [(1+\alpha)\,\mathbf C\,\dot{\mathbf u}_{n+1} - \alpha\,\mathbf C\,\dot{\mathbf u}_n]
  - \mathbf M\,\ddot{\mathbf u}_{n+1}$$

$$\mathbf J \;=\; \mathbf M + (1+\alpha)\,\gamma\Delta t\,\mathbf C + (1+\alpha)\,\beta\Delta t^2\,\mathbf K_\text{t}$$

Todo lo demás (criterio de convergencia, line search opt-in del ADR 0011,
telemetría tipada de divergencia, commit/rollback de estado, amortiguamiento
Rayleigh constante en el tiempo) se hereda de `NewtonNewmarkSolver` sin
cambios.

---

## Contrato de implementación

```yaml
name: HHTSolver
kind: solver
status: validated
extends: NewmarkSolver       # variante según Reglas §4

interface:
  type_yaml: HHTSolver
  pipeline_kind: transient    # mismo que NewmarkSolver

parameters_extra_to_NewmarkSolver:
  - { name: alpha,  type: float, required: false, default: -0.05,
      desc: "Parámetro HHT en [-1/3, 0]. α=0 recupera Newmark trapezoidal.
             α=-0.05 (default) es disipación ligera estándar; α=-1/3 disipación máxima." }

parameters_overridden_with_autoderivation:
  - { name: beta,   type: float, required: false, default: "(1-alpha)^2/4",
      desc: "Auto-derivado del alpha si no se pasa; override explícito posible
             (no recomendado salvo experimentación)." }
  - { name: gamma,  type: float, required: false, default: "(1-2*alpha)/2",
      desc: "Auto-derivado del alpha si no se pasa." }

parameters_inherited:
  # Idénticos a NewmarkSolver: t_end, dt, rayleigh, u0, u0_dot, F_func,
  # linear_algebra, lumping. Ver docs/specs/NewmarkSolver.md.

acceptance:
  recovery:
    - name: alpha_zero_matches_NewmarkSolver
      setup: "Oscilador 1 GDL con HHTSolver(alpha=0.0) y NewmarkSolver(beta=0.25, gamma=0.5)"
      expect: "u_history coincide bit-a-bit (atol 1e-10)"

  damping:
    - name: high_frequency_dissipation
      setup: "Dos osciladores acoplados con ω₁=1 rad/s y ω₂=100 rad/s; condición inicial excita ambos modos"
      expect: "HHTSolver(alpha=-0.1) atenúa el modo alto en >50% tras 5 periodos del modo bajo,
               manteniendo amplitud del modo bajo dentro de 5% del valor inicial"

    - name: spectral_radius_at_high_frequency_limit
      setup: "Δt grande (>10 periodos del modo); HHTSolver con varios α"
      expect: "ρ_∞ observado = (1+α)/(1-α) dentro de tolerancia 5%"

    - name: amplification_matrix_spectral_radius_analytic
      setup: "Oscilador 1-GDL con HHTSolver(α∈{-0.05, -0.10, -1/3}); Δt tal que Ω=ω·Δt=2"
      expect: |
        Razón de contracción por paso medida en la simulación
        coincide con max|λ(A)|, donde A es la matriz de amplificación 3×3 HHT-α
        (Hughes 1987 §9.3, Bathe 2014 §9.5.4): D=1+β(1+α)Ω², A[0,0]=(1+αβΩ²)/D, …
        Tolerancia 1%.
      ref: "tests/test_hht.py::TestHHTNumericalDampingAnalytic"

  stability:
    - name: unconditional_stability_large_dt
      setup: "Sistema lineal con Δt = 100·T (T periodo natural), 10 pasos"
      expect: "u_history acotada (no diverge) para todo α ∈ [-1/3, 0]"

  nonlinear:
    - name: oscilador_plastico_con_disipacion_HHT
      setup: "Truss elastoplástico 1 GDL, vibración libre desde u_0 más allá del yield"
      expect: "NewtonHHTSolver(alpha=-0.1) reduce oscilación espuria post-yielding
               manteniendo decay físico de la disipación plástica"

references:
  - "Hilber, H. M., Hughes, T. J. R. & Taylor, R. L. (1977). Improved numerical
     dissipation for time integration algorithms in structural dynamics.
     Earthquake Eng. Struct. Dyn. 5(3), 283-292."
  - "Hughes, T.J.R. (2000). The Finite Element Method, §9.3 (HHT-α y generalized-α)."
  - "Chopra, A. K. (2017). Dynamics of Structures, §15.4 (integradores con disipación numérica)."
  - "ADR 0009 §Hoja de ruta — fila 3 (Newmark / HHT-α)."
  - "Spec padre: [docs/specs/NewmarkSolver.md](NewmarkSolver.md)."
  - "Spec hermana: [docs/specs/NewtonNewmarkSolver.md](NewtonNewmarkSolver.md) — padre del NewtonHHTSolver."
```

---

## Implementación

- **Archivo**: [solidum/math/solvers/newmark.py](../../solidum/math/solvers/newmark.py) — junto a `NewmarkSolver` y `NewtonNewmarkSolver` (lógicamente cercanos; comparten predictores Newmark, reducción Dirichlet y amortiguamiento Rayleigh).
- **Clases**:
  - `HHTSolver(NewmarkSolver)` — variante lineal con disipación numérica.
  - `NewtonHHTSolver(NewtonNewmarkSolver)` — variante no lineal con Newton interno + disipación numérica.
  Ambas registradas con `@SolverRegistry.register`.
- **Auto-derivación de β, γ**: en `__init__`, si el usuario no pasa β o γ explícitos, se computan desde α antes de invocar `super().__init__`. Esto asegura que el padre vea la β, γ ya correctas y la mayor parte del flujo herede sin cambios.
- **Solve overrideado**: `solve()` se reescribe en ambas subclases para usar las ecuaciones HHT-α (ver §"Sistema efectivo" y §"Variante no lineal" arriba). El estado anterior `F_int_n`, `F_n`, `u̇_n` necesario para el término `−α·X_n` se mantiene en variables locales del bucle.
- **Despacho YAML**: `solver: type: HHTSolver` (o `NewtonHHTSolver`). Como `HHTSolver` hereda de `NewmarkSolver` y `NewtonHHTSolver` hereda de `NewtonNewmarkSolver`, el `isinstance` en `entry.run_yaml` los detecta como `transient` automáticamente. **No dispara la regla C arquitectural** (que esperaría una rama no-clásica que NO sea subclase de las existentes).
- **Tests**: `tests/test_hht.py` nuevo con los acceptance arriba.

---

## Diálogo

- **2026-05-18** · Spec creada como **variante** de `NewmarkSolver` (y `NewtonHHTSolver` como variante de `NewtonNewmarkSolver`) aplicando la regla §4. Sin ADR nuevo: reusa la decisión del ADR 0009 fila 3 ("Transitorio lineal Newmark / HHT-α"). Sin entrada propia en el manual estructural — visible al usuario solo a través del `type` YAML y el parámetro `alpha`.

- **2026-05-18** · Decisión sobre el default de `alpha`: **−0.05** (disipación ligera estándar) en lugar de `0` (sin disipación). Razones físicas: (i) si el usuario pide `HHTSolver`, está pidiendo disipación — default `α=0` lo convierte en un alias trivial de `NewmarkSolver`, ruido en el catálogo; (ii) `−0.05` es el valor canónico de Hilber 1977 (artículo original) y de los códigos de referencia (Abaqus, OpenSees) por ser balance entre disipación útil (>9% por ciclo en altas frecuencias) y precisión preservada (<0.5% en modos de interés); (iii) no introduce regresiones — `HHTSolver` es un solver nuevo, nadie lo usa hoy. Diferente del default `line_search=False` del ADR 0011, que era *no* poder romper código existente; aquí no hay código existente que romper.

- **2026-05-18** · Decisión sobre auto-derivación de β, γ: por defecto sí, desde `alpha`. Override explícito posible pero no recomendado (combinaciones α, β, γ arbitrarias pueden perder estabilidad incondicional u orden 2). El parámetro de la API es **el** α, no la terna; β, γ son consecuencia matemática de α salvo experimentación expresa. Esto simplifica el uso típico y centraliza la responsabilidad: si el usuario solo conoce α, el solver hace el resto.
