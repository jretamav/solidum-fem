# Evolución de la arquitectura

Este capítulo presenta la trayectoria prevista de la arquitectura de Fenix FEM. Su función es preservar la imagen de la arquitectura del programa a medida que crece: cada fase incorpora un nuevo tipo de análisis con sus contratos propios, sus variables de estado y, eventualmente, capas o subsistemas adicionales. La diferencia con los pendientes del capítulo 4 es de granularidad: aquellos enumeran componentes individuales; este capítulo describe los **estados sucesivos de la arquitectura** del programa, cada uno con su mapa de capas correspondiente.

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

- Aparece la matriz de amortiguamiento `C = α·M + β·K` (Rayleigh proporcional). Los coeficientes pueden pasarse directos o calibrarse modalmente a partir de `(ξ₁, ω₁), (ξ₂, ω₂)` con `fenix.math.damping.rayleigh_from_modes`.
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

**Limitaciones declaradas en esta fase**:

- Solo lineal: `K, M` constantes. No-linealidad (Newton dentro de Newmark) diferida a la fase 4.
- Apoyos Dirichlet constantes en el tiempo. Multi-support excitation sísmica diferida.
- `Δt` fijo; adaptativo diferido.
- Esquemas con amortiguamiento numérico controlado (HHT-α, generalized-α) diferidos a subclase futura.

## Fase 4 — Dinámica estructural transitoria no lineal

**Estado**: [PENDIENTE: Implementación de la fase 4 — Newton-Raphson dentro de cada paso de Newmark.]

**Cambios respecto a la fase 3**:

- En cada paso temporal, el residuo `R = F(t_{n+1}) − F_int(u_{n+1}) − M·ü_{n+1} − C·u̇_{n+1}` se anula iterativamente con Newton-Raphson, reutilizando la infraestructura del `NonlinearSolver` estático (ADR 0007 — `ConvergenceCriterion` calibrado).
- La matriz tangente efectiva incorpora la rigidez dependiente del estado: `A_eff = M + γΔt·C + βΔt²·K_t(u)`. Se refactoriza cuando `K_t` cambia significativamente; con `K_t` constante (régimen lineal local) se reaprovecha la factorización del paso anterior.
- El estado interno de los elementos (`ElementState`) recupera plenamente su semántica *trial* / comprometido: se commit al converger cada paso temporal, no solo cada paso de carga.

**Topología prevista**: igual que la fase 3, con el bucle Newton añadido dentro del bucle temporal.

```callout Riesgo de arquitectura
Bajo. Toda la infraestructura ya existe — el `NonlinearSolver` estático tiene calibración de convergencia, fallback de positividad, y Newton modificado. La fase 4 es composición de las fases 1 y 3.
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

De este modo, el capítulo se convierte en una crónica ordenada del crecimiento del programa: en cualquier momento futuro, abrir el manual permite reconstruir no solo el estado actual de Fenix FEM sino la trayectoria completa que lo ha llevado hasta él.
