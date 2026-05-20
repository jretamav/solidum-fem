# Evolución de la arquitectura

Este capítulo presenta la trayectoria prevista de la arquitectura de Solidum FEM. Su función es preservar la imagen de la arquitectura del programa a medida que crece: cada fase incorpora un nuevo tipo de análisis con sus contratos propios, sus variables de estado y, eventualmente, capas o subsistemas adicionales. La diferencia con los pendientes del capítulo 4 es de granularidad: aquellos enumeran componentes individuales; este capítulo describe los **estados sucesivos de la arquitectura** del programa, cada uno con su mapa de capas correspondiente.

Las fases posteriores a la actual están marcadas como pendientes y se materializarán mediante los Architecture Decision Records (ADR) anunciados en el capítulo 8. La estructura de cada fase se redefinirá en su ADR de apertura; lo aquí descrito constituye una hipótesis razonada de partida.

## Fase 1 (actual) — Estática lineal y no lineal

**Estado**: implementado.

**Capas activas**: las seis capas descritas en el capítulo 2 (entrada, inicialización, interpretación del caso, dominio, capa numérica, salida).

**Tipos de análisis**: estática lineal (`LinearSolver`) y estática no lineal incremental (`NonlinearSolver`, `ArcLengthSolver`).

**Variables de estado**: desplazamientos `U`. Variables internas del material (deformaciones plásticas, daño, historia) gestionadas mediante `ElementState` con semántica *trial* / comprometido.

**Matrices globales**: matriz de rigidez `K` (única).

**Contrato de elemento**: `K_local`, `f_internal`, `internal_forces(U)`. La matriz `B` y el residuo se calculan en cada iteración.

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K, solvers] - [Salida]
```

## Fase 2 — Análisis modal

**Estado**: implementado (ADR 0009 fase 1, sesión del 2026-05-13).

El análisis modal resuelve el problema generalizado de valores y vectores característicos `K·φ = ω²·M·φ`. Es un problema algebraico cuyas soluciones admiten **interpretación** como frecuencias propias y formas de vibración libre no amortiguada, pero el análisis en sí es estático: no hay paso de tiempo ni evolución temporal. Se acometió antes que la dinámica transitoria porque solo requiere la matriz de masa nueva, aísla esa pieza, y deja la infraestructura preparada para Newmark.

**Cambios respecto a la fase 1**:

- Aparece la matriz de masa `M`, ensamblada con la misma infraestructura que `K` (mismo patrón de sparsity, misma topología COO cacheada).
- Cada elemento implementa `compute_mass_matrix(lumping)` — masa consistente translacional invariante en barras y cables, axial + Hermitiana cúbica + inercia rotacional propia de sección en frames 2D/3D, integral `∫ρ·NᵀN·t dA` en sólidos 2D (Tri6 usa cuadratura específica `tri_6` por subintegración con `tri_3`).
- Nuevo subsistema en la capa algebraica: `EigenSolver` que envuelve ARPACK Lanczos con shift-invert (`scipy.sparse.linalg.eigsh`).
- Nueva clase de solver: `ModalSolver`, registrado en `SolverRegistry` (no se fragmenta el registro). Reutiliza Dirichlet por eliminación directa con un `reduce_pair(K, M)` que aplica `T` a ambas matrices simultáneamente.
- Nueva salida especializada: `ModalResult` (frozen) con `frequencies_rad/hz`, `periods`, `modes` M-ortonormales, y un método `free_vibration(M, u0, u0_dot, t)` que reconstruye analíticamente la respuesta libre sin amortiguamiento por superposición modal — sin necesidad de integrador temporal.
- `Material.density` (ADR 0008) ya disponía de la única propiedad física nueva necesaria.

**Topología efectiva**:

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K, M, eigensolver] - [Salida: ModalResult]
```

```callout Riesgo de arquitectura
Bajo. El análisis modal añade un nuevo tipo de problema algebraico (autovalor generalizado en lugar de sistema lineal) y una nueva familia de resultados, sin modificar contratos existentes ni la infraestructura de Dirichlet.
```

## Fase 3 — Dinámica estructural transitoria lineal

**Estado**: implementado (ADR 0009 fase 3, sesión del 2026-05-13).

Integración temporal directa de la ecuación de movimiento semidiscreta `M·ü + C·u̇ + K·u = F(t)` por el método de Newmark-β. Es el primer análisis con dependencia temporal real; aprovecha la matriz de masa introducida en la fase 2 y reutiliza las decisiones arquitecturales fijadas en el ADR 0009 (registro unificado, estado dinámico fuera de `Node`, amortiguamiento Rayleigh como entrada estándar).

**Cambios respecto a la fase 2**:

- Aparece la matriz de amortiguamiento `C = α·M + β·K` (Rayleigh proporcional). Los coeficientes pueden pasarse directos o calibrarse modalmente a partir de `(ξ₁, ω₁), (ξ₂, ω₂)` con `solidum.math.damping.rayleigh_from_modes`.
- Las variables de estado se amplían a `(u, u̇, ü)`. Coherente con el ADR 0009 §5, el estado dinámico **no contamina** la clase `Node`: vive en el nuevo `TransientResult` (frozen) como historiales `t_history`, `u_history`, `udot_history`, `uddot_history` indexados por paso temporal.
- Nuevo solver `NewmarkSolver` con esquema `(β, γ)` parametrizable (default 1/4, 1/2 — average acceleration, incondicionalmente estable, error `O(Δt²)`). Cargas dependientes del tiempo por callback Python (`F_func: t → F`).
- El subsistema algebraico recibe una matriz efectiva `A_eff = M + γΔt·C + βΔt²·K`. Con `Δt` constante y problema lineal, `A_eff` es invariante: se factoriza **una sola vez** al inicio (ADR 0003 — `FactorizedSolver`) y cada paso temporal se reduce a una sustitución triangular barata.
- Aceleración inicial calculada de forma consistente con la ecuación de movimiento: `M·ü₀ = F(0) − C·u̇₀ − K·u₀`.

**Topología efectiva**:

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K, M, C, A_eff, solver] - [Salida: TransientResult]
                                                                                  |
                                                                  [Integrador temporal Newmark]
```

```callout Riesgo de arquitectura
Bajo. La decisión 5 del ADR 0009 (estado dinámico fuera de `Node`, en `TransientResult`) preserva la pureza topológica/geométrica de `Node` y evita que análisis estáticos paguen el coste de campos que no usan. La semántica *trial* / comprometido del estado interno de elementos (ADR sobre `ElementState`) sigue intacta porque la fase 3 es lineal: no hay variables internas que evolucionen.
```

**Limitaciones declaradas en esta fase 3**:

- Solo lineal: `K, M` constantes. No-linealidad (Newton dentro de Newmark) entregada en la fase 4.
- Apoyos Dirichlet constantes en el tiempo. Multi-support excitation sísmica diferida.
- `Δt` fijo; adaptativo diferido.
- Esquemas con amortiguamiento numérico controlado (HHT-α) entregados en 2026-05-18 como variante de la familia Newmark.

## Fase 4 — Dinámica estructural transitoria no lineal

**Estado**: implementado (ADR 0009 fase 4, sesión del 2026-05-14).

**Cambios respecto a la fase 3**:

- Nueva clase `NewtonNewmarkSolver` como **variante** de `NewmarkSolver` (Reglas §4 — extensión documentada en spec corta, sin ADR nuevo porque reúsa todas las decisiones del ADR 0009): hereda predictores, correctores y reducción de Dirichlet; añade un bucle Newton-Raphson dentro de cada paso temporal.
- En cada paso, el residuo `R = F_ext(t_{n+1}) − F_int(u_{n+1}) − C·u̇_{n+1} − M·ü_{n+1}` se anula iterativamente con jacobiano `J = M + γΔt·C + βΔt²·K_t`, donde `K_t = ∂F_int/∂u` es la rigidez tangente corriente. Convergencia dual del ADR 0007 calibrada al primer paso.
- Reaprovecha íntegramente la maquinaria del `NonlinearSolver` estático: `ConvergenceCriterion`, `freeze_tangent_after_iter` (Newton modificado, ADR 0003 fase 2), `commit_all_states()` tras converger cada paso temporal — el estado interno de los elementos recupera plenamente su semántica *trial* / comprometido también en dinámica.
- Amortiguamiento Rayleigh `C = α·M + β·K_0` **constante en el tiempo**, calibrado con la rigidez elástica de referencia `K_0` al inicio del análisis. Convención estándar (Abaqus, ANSYS, OpenSees); evita acoplamiento ad-hoc entre disipación viscosa y plástica.
- **Recuperación del caso lineal**: con materiales sin historia (o plásticos no plastificados), `K_t ≡ K_0` y el residuo se anula en una iteración. El resultado coincide a paridad de bits con `NewmarkSolver`, validado en tests.
- Despacho YAML automático: como `NewtonNewmarkSolver` es subclase de `NewmarkSolver`, `entry.run_yaml` lo detecta vía `isinstance` y enruta a `run_transient` sin tocar el dispatcher.

**Topología efectiva**: igual que la fase 3, con el bucle Newton añadido dentro del bucle temporal y la rigidez tangente reensamblada en cada iteración.

```callout Riesgo de arquitectura
Bajo en la práctica: toda la infraestructura ya existía. El `NonlinearSolver` estático aportaba calibración de convergencia, fallback de positividad y Newton modificado; la fase 3 aportaba predictores/correctores Newmark y reducción de Dirichlet. La fase 4 es la composición exacta de ambas, formalizada como subclase para no duplicar código.
```

## Fase 4 (extensión 2026-05-18) — Cierre completo del subsistema modal/dinámico/espectral

**Estado**: implementado (ADR 0009 fases 2, 5, 6, 7 + variante HHT-α + reglas C y D de auditoría aplicadas).

Las cuatro fases pendientes del ADR 0009 se cierran en una secuencia obligada (cada una precondición de la siguiente) sin reabrir las decisiones arquitecturales tomadas en el ADR original:

**Mass lumping (fase 2)** — `compute_mass_matrix(lumping="lumped")` operativo en todos los elementos. Helper centralizado `solidum.math.mass_lumping.lump_hrz` con esquema HRZ canónico (Hinton-Rock-Zienkiewicz 1976) para sólidos isoparamétricos y fórmula nodal directa para frames. Diagonalidad estricta en globales para todos los elementos excepto Frame3D con eje oblicuo (queda bloque-diagonal — limitación documentada, Cook-Malkus-Plesha §11.4).

**HHT-α (variante)** — `HHTSolver` y `NewtonHHTSolver` como variantes de la familia Newmark (Reglas §4 — spec corta sin ADR nuevo). Disipación numérica controlada por `α ∈ [−1/3, 0]` con `(β, γ)` auto-derivados canónicamente.

**Central differences (fase 5)** — `CentralDifferenceSolver` con leapfrog explícito Belytschko-Liu-Moran. Una sola clase cubre lineal y no lineal con parámetro `nonlinear`. Requiere masa lumped (`lumping="consistent"` rechazado con `ValueError`). Estabilidad condicional CFL con detección a posteriori de divergencia exponencial.

**Harmonic (fase 6)** — `HarmonicSolver` resuelve `(−ω²M + iωC + K)·û = F̂` con aritmética compleja para un barrido en `ω`. Factorización LU compleja por frecuencia. Nuevo `HarmonicResult` con métodos `.amplitude()` y `.phase()`. Cuarto valor del literal `PIPELINE_KIND` (ver regla C abajo).

**Response spectrum (fase 7)** — `ResponseSpectrumSolver` combina respuestas máximas modales bajo un espectro `S_d(ω)` o `S_a(ω)` por SRSS (default) o CQC (Der Kiureghian 1980, para modos cercanos). Devuelve `ResponseSpectrumResult` con respuesta envolvente, contribución modal individual, factores de participación, masas efectivas y método `.cumulative_effective_mass_ratio()` para verificación normativa (≥ 0.9).

**Reglas C y D aplicadas durante esta extensión** (auditoría arquitectural 2026-05-13):

- **Regla C** — `entry.py::run_yaml` despacha por atributo de clase `PIPELINE_KIND ∈ {"static", "modal", "transient", "harmonic", "spectrum"}` en lugar de cadena de `isinstance`. Solvers no clásicos futuros no requieren tocar `entry.py`.
- **Regla D** — `solidum/math/modal_response.py` agrupa todos los cómputos sobre modos: `free_vibration` (movido desde `results.py`), `participation_factors`, `response_spectrum_srss`, `response_spectrum_cqc`, helpers de espectros (`spectrum_from_sa`, `spectrum_tabulated`). `ModalResult.free_vibration` queda como wrapper delgado preservando la API histórica. `results.py` vuelve a su propósito declarado: dataclasses inmutables, no algoritmos.

**Topología efectiva** tras la extensión:

```
[Entrada] - [Inicializacion] - [Interprete: PIPELINE_KIND] - [Dominio] - [Numerica: K, M, C, lump_hrz, modal_response] - [Salida]
                                                                                                                          |
                            [run, run_modal, run_transient, run_harmonic, run_response_spectrum]    →    [Result type específico]
```

```callout Riesgo de arquitectura
Bajo. Las cinco piezas son extensiones aditivas: nuevo módulo `mass_lumping.py`, nuevo módulo `modal_response.py`, cinco solvers nuevos en `solidum/math/solvers/`, dos `Result` nuevos (`HarmonicResult`, `ResponseSpectrumResult`). El despacho declarativo `PIPELINE_KIND` reemplaza la cadena de `isinstance` sin romper contratos — los solvers existentes heredan el nuevo atributo. La regla D centraliza algoritmos sin modificar la API pública de `ModalResult`.
```

## Fase 5 — Problema térmico estacionario y transitorio

**Estado**: [PENDIENTE: Implementación de la fase 5 — problema térmico.]

**Cambios respecto a las fases anteriores**:

- Aparece un campo escalar nuevo: la temperatura `T`. Esto obliga a generalizar el concepto de DOF: hasta ahora todos los DOF eran mecánicos (desplazamientos, rotaciones); ahora se incorporan DOF térmicos.
- Aparece la matriz de conductividad `K_t` (análogo térmico de la rigidez) y la matriz de capacidad `C_t` (análogo térmico de la masa).
- Se introduce una nueva familia de elementos térmicos, paralela a la mecánica. Cada elemento térmico declara `DOF_NAMES = ['T']` y consume materiales con propiedades térmicas (conductividad, capacidad, fuente).
- Aparece una nueva familia de materiales térmicos con interfaz propia: `compute_flux_and_tangent(grad_T, state)`.
- En transitorio, se reutilizan los integradores temporales de la fase 2, adaptados al sistema térmico (Crank-Nicolson, theta-método).

**Topología prevista**:

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K_m, M_m, K_t, C_t, solvers] - [Salida]
                                                  |
                              [DOF mecanicos] + [DOF termicos]
                              [Materiales mec.] + [Materiales term.]
                              [Elementos mec.]  + [Elementos term.]
```

```callout Riesgo de arquitectura
Medio. Es la primera fase que rompe la suposición implícita "todos los DOF son mecánicos" presente hoy en numerosos puntos del código. La generalización requiere atributo `DOF_KIND` o similar en cada DOF (mecánico, térmico, futuro presión, futuro químico); posible separación entre `MechanicalElement` y `ThermalElement` con base común `Element`, o un único `Element` polimórfico (decisión a tomar en el ADR de apertura); adaptación del `VtkExporter` para escribir campos escalares además de vectoriales.
```

## Fase 6 — Acoplamiento termo-mecánico

**Estado**: [PENDIENTE: Implementación de la fase 6 — acoplamiento termo-mecánico.]

**Cambios respecto a la fase 5**:

- Aparece la interacción entre los campos térmico y mecánico: dilatación térmica que induce deformaciones, disipación mecánica que genera calor.
- Dos estrategias distintas en términos de arquitectura, a decidir en el ADR de apertura:
  - *Estrategia desacoplada (staggered)*: en cada paso, se resuelve primero el problema térmico con los desplazamientos de la última iteración, después el problema mecánico con las temperaturas recién calculadas. Cada subsistema mantiene su propia matriz; el acoplamiento se realiza mediante intercambio de campos. Implementación de bajo riesgo de arquitectura pero convergencia limitada en problemas de fuerte interacción.
  - *Estrategia monolítica*: se ensambla un sistema único `K_acoplada` que combina los DOF mecánicos y térmicos. Convergencia robusta pero exige una nueva infraestructura de ensamblaje multifísico.
- Aparece una nueva familia de materiales acoplados con interfaz extendida: `compute_response(strain, T, grad_T, state)`.

**Topología prevista (estrategia monolítica)**:

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K_acoplada, M_acoplada, solvers] - [Salida]
                                                  |
                                       [Materiales termo-mec.]
                                       [Elementos termo-mec.]
                                       [Solver acoplado]
```

```callout Riesgo de arquitectura
Alto. Es la primera fase con auténtica multifísica. Requiere generalización del subsistema algebraico para manejar matrices con bloques no homogéneos (un bloque mecánico simétrico positivo definido y un bloque térmico no simétrico, en la formulación habitual); generalización del concepto de prueba de validación, ya que no basta una solución analítica monodisciplinaria y los problemas de referencia deben tener acoplamiento conocido; posible introducción de precondicionadores específicos para sistemas acoplados (precondicionamiento por bloques, complemento de Schur).
```

## Mantenimiento de este capítulo

Cada vez que se cierre la implementación de una fase y se acepte su ADR de apertura, esta sección se actualiza:

1. La fase pasa de "[PENDIENTE: ...]" a "implementado".
2. La hipótesis inicial sobre la arquitectura se sustituye por la descripción real, ajustada a las decisiones que se hayan tomado durante la implementación.
3. El diagrama de capas se redibuja conforme a la topología efectiva.
4. Las consecuencias no previstas (ajustes en capas anteriores, contratos modificados) se registran en una subsección "Lecciones de la fase N", de cara a fases posteriores.

De este modo, el capítulo se convierte en una crónica ordenada del crecimiento del programa: en cualquier momento futuro, abrir el manual permite reconstruir no solo el estado actual de Solidum FEM sino la trayectoria completa que lo ha llevado hasta él.
