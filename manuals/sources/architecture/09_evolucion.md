# Evolución de la arquitectura

Este capítulo presenta la trayectoria prevista de la arquitectura de Solidum FEM. Su función es preservar la imagen de la arquitectura del programa a medida que crece: cada fase incorpora un nuevo tipo de análisis con sus contratos propios, sus variables de estado y, eventualmente, capas o subsistemas adicionales. La diferencia con los pendientes del capítulo 4 es de granularidad: aquellos enumeran componentes individuales; este capítulo describe los **estados sucesivos de la arquitectura** del programa, cada uno con su mapa de capas correspondiente.

Las fases aún no implementadas están marcadas como pendientes y se materializarán mediante los Architecture Decision Records (ADR) anunciados en el capítulo 8. La estructura de cada fase se redefinirá en su ADR de apertura; lo aquí descrito constituye una hipótesis razonada de partida.

## Fase 1 — Estática lineal y no lineal

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

**Estado**: implementada en su núcleo de conducción pura (Etapa 8, agosto de 2026). Quedan fuera la convección, la radiación y toda forma de no linealidad térmica.

**Cambios efectivamente introducidos**:

- Un campo escalar nuevo, la temperatura `T`, con elementos que declaran `DOF_NAMES = ['T']`.
- Las matrices de conductividad y de capacidad calorífica, análogos térmicos de la rigidez y de la masa.
- Una familia de elementos térmicos paralela a la mecánica (`Quad4Thermal`, `Hex8Thermal`) sobre base común `_ThermalSolid`, y una familia de materiales térmicos con registro e interfaz propios.
- Un solver de integración temporal de primer orden, `ThetaMethodSolver`, con tipo de resultado específico.

**Contraste entre lo previsto y lo ocurrido.** Esta fase se había estimado de riesgo *medio* por ser «la primera que rompe la suposición implícita de que todos los DOF son mecánicos, presente en numerosos puntos del código». Se anticipaban tres consecuencias: un atributo `DOF_KIND` para clasificar cada grado de libertad, una posible escisión del contrato `Element` en ramas mecánica y térmica, y la adaptación del exportador de resultados.

Ninguna de las dos primeras resultó necesaria. La suposición que se temía no estaba en el código: el ensamblador, la imposición de restricciones por eliminación y el despacho algebraico operan sobre índices de grado de libertad y sobre propiedades de la matriz, sin consultar en ningún punto el significado físico del campo. Las condiciones de frontera y las cargas nodales se aplican por *nombre* de grado de libertad, de modo que `T` funcionó por el mismo mecanismo genérico que `ux`. El régimen estacionario, en consecuencia, no requirió solver alguno: el solver lineal existente resolvió el problema de conducción sin una sola línea de modificación.

La tercera predicción sí se confirmó al cerrar la etapa: **el exportador VTK escribía únicamente campos mecánicos** —desplazamientos y rotaciones— y no exportaba el campo de temperatura. Era la limitación práctica más visible del subsistema térmico; desde el 22 de septiembre de 2026 el exportador escribe el campo nodal `Temperature` cuando algún nodo tiene el grado de libertad `T`.

La lección de arquitectura que deja la fase es que el coste real de la extensión no se concentró donde la previsión lo situaba. Lo caro no fue generalizar la infraestructura —que ya era agnóstica al campo, aunque nadie lo hubiera comprobado— sino la física propia del problema: la elección de los esquemas por defecto contra el principio del máximo de la ecuación de difusión, y la constatación de que la ecuación de conducción, al ser de primer orden en el tiempo, no admite reutilizar la familia de integradores de Newmark. La interfaz del material tampoco fue la anticipada `compute_flux_and_tangent(grad_T, state)`, sino `compute_flux(grad_T)` sin argumento de estado, porque la conducción de Fourier carece de variables internas con historia.

**Topología resultante**:

```
[Entrada] - [Inicializacion] - [Interprete] - [Dominio] - [Numerica: K_m, M_m, K_t, C_t, solvers] - [Salida]
                                                  |
                              [DOF mecanicos] + [DOF termicos]
                              [Materiales mec.] + [Materiales term.]
                              [Elementos mec.]  + [Elementos term.]
```

```callout Riesgo de arquitectura
Resultó **bajo**, frente al riesgo medio estimado. La infraestructura existente absorbió el campo escalar sin generalización previa, lo que constituye una verificación empírica de las decisiones de despacho algebraico e imposición de restricciones frente a un campo físico que no existía cuando se tomaron. La deuda que quedó era de salida, no de arquitectura (la exportación del campo de temperatura a VTK), y se saldó el 22 de septiembre de 2026.
```

## Fase 6 — Acoplamiento termo-mecánico

**Estado**: estrategia decidida en el ADR 0016 (22 de septiembre de 2026), sin implementar.

**Lo que ya existe.** Antes de decidir se verificó sobre el código que la infraestructura admite nodos con desplazamiento y temperatura a la vez: la numeración de ecuaciones itera los grados de libertad por nombre, el ensamblador reúne ambos campos en una sola matriz, la reducción acepta condiciones mixtas sobre un mismo nodo, masa y capacidad calorífica se ensamblan con el mismo operador y un análisis estático de dos campos se resuelve de extremo a extremo. La previsión que figuraba aquí —una nueva infraestructura de ensamblaje multifísico— resultó innecesaria.

**Lo que falta** es la física de los bloques fuera de la diagonal: la dilatación térmica (el campo de temperatura induce fuerzas mecánicas, efecto de primer orden) y el calor de deformación (efecto de segundo orden salvo en régimen adiabático).

**Decisiones del ADR 0016**:

- *Acoplamiento débil unidireccional primero*: se resuelve el problema térmico y su campo entra en el mecánico, sin retroalimentación. El monolítico queda como extensión aditiva cuando exista el caso real, sobre la misma infraestructura.
- *Deformación propia genérica* `ε₀`, con `σ = C:(ε − ε₀)`: la térmica es el caso particular `ε₀ = α·ΔT`, y la higroscópica —dominante en madera y bambú— entra por el mismo canal sin volver a tocar el contrato de material.
- `ε₀` vive **por punto de Gauss**, no en el material (compartido por todos los elementos de una familia de lote) ni en el elemento (varía dentro de él).
- *Sin familia de elementos acoplados*: el elemento mecánico interpola el campo de temperatura con sus propias funciones de forma, con mallas conformes.

```callout Riesgo de arquitectura
Bajo, tras la verificación del ADR 0016. El cambio de contrato se limita a un canal aditivo para `ε₀` con valor por defecto cero, que deja los materiales existentes intactos.
```

## Fase 7 — Rendimiento y robustez de la capa numérica

**Estado**: implementado (ADR 0014, 0015, 0017, 0018 y 0019, 22 y 23 de septiembre de 2026; capa 4 de la red de seguridad, 27 de octubre de 2026).

Ninguna de estas piezas añade un tipo de análisis: cambian cómo se ensambla, cómo se resuelve y qué se comprueba, sin cambiar un resultado.

- **Ensamblaje por lotes** (ADR 0014): la familia de elementos se deriva de los contratos declarativos y se evalúa en un kernel compilado, serie o paralelo, con el camino por elemento como referencia obligatoria; el estado interno se declara (`STATE_SCHEMA`) y se guarda en arreglos.
- **Corrector de Newton compartido** (ADR 0015): un solo bucle de corrección para todos los solvers iterativos; cada solver aporta su problema por paso y su control de paso.
- **Backends algebraicos** (ADR 0017 y 0018): Pardiso multihilo opcional, y un iterativo (CG/MINRES con multimalla algebraica sembrada con los modos de cuerpo rígido) sólo a petición.
- **Red de seguridad del análisis estático** (ADR 0019): mecanismos rígidos detectados antes de resolver y descritos en términos del modelo; mecanismos internos, solución del sistema y pérdida de cifras por mal condicionamiento comprobados después.

```
[Numerica] = [Ensamblaje por familias (Numba)] + [NewtonCorrector] + [Despacho: Cholesky | Pardiso | LU | iterativo] + [Red de seguridad]
```

```callout Riesgo de arquitectura
Bajo. Todo es interno a la capa numérica y se verificó contra el camino anterior a precisión de máquina; el contrato que gana un componente nuevo es opcional (`BATCH_KINEMATICS`, `BATCH_KERNEL`) salvo `STATE_SCHEMA`, una línea.
```

## Fase 8 — Módulos de usuario y frontera del programa principal

**Estado**: implementado (ADR 0020, 24 de septiembre de 2026; revisión del 29 de septiembre con el criterio de admisión de `Reglas.md` §1).

Las formulaciones no estándar salen del programa principal a módulos de usuario (`solidum/user/<nombre>/`), con carga explícita (`user_modules` en el YAML o `solidum.load_user_module`). El programa principal sólo ofrece piezas genéricas: registros con metadatos de familia, un lector YAML con un solo bucle sobre las familias, referencias declaradas por el elemento (`REFERENCE_KWARGS`), el gancho `prepare_step` y la cinemática pública. Dos módulos: `discontinuities` (discontinuidades interiores embebidas, ADR 0010) y `dissipation_control` (arco por disipación). Un test de pureza vigila que el programa principal no los nombre.

```
[Programa principal: solidum/]  <-- piezas genéricas --  [Modulos de usuario: solidum/user/<nombre>/]
```

```callout Riesgo de arquitectura
Bajo. La frontera es aditiva y la vigilan tres controles automáticos (estructural, léxico y en ejecución). Una pieza estándar que hoy sólo usa un módulo sigue en él hasta que la necesite el programa principal; así volvió la familia cohesiva el 10 de noviembre de 2026 (fase 12).
```

## Fase 9 — Secciones de marco

**Estado**: implementado (ADR 0021, 29 de septiembre de 2026).

Una cuarta familia de objetos con sección YAML propia, `SectionRegistry` (`sections`), sin estado: el estado de cada fibra vive en el elemento. `FiberSection` da a los marcos Euler-Bernoulli plasticidad distribuida con cualquier material 1D; los marcos aceptan un material elástico con propiedades geométricas o una sección con `n_ip` puntos de Gauss-Legendre.

```callout Riesgo de arquitectura
Bajo para la arquitectura, alto para los modelos existentes: un marco con material inelástico sin sección deja de construirse, porque antes daba flectores incorrectos en silencio.
```

## Fase 10 — Cargas variables y verificación por soluciones manufacturadas

**Estado**: implementado (ADR 0022 y 0023, 27 de octubre de 2026).

Las cargas distribuidas aceptan campos de la posición integrados sobre la geometría real, también desde el YAML con expresiones compiladas por un intérprete con lista blanca. El contrato declarativo de la celda de referencia (`REFERENCE_CELL`, `REFERENCE_NODES`) habilita un generador de mallas estructuradas (`solidum.mesh`) y la verificación por soluciones manufacturadas (`solidum.verification`), que recorre el registro: un elemento nuevo queda verificado sin escribir pruebas.

```callout Riesgo de arquitectura
Bajo. Extensiones aditivas de las cargas y de la entrada; la verificación es una herramienta, no una capa.
```

## Fase 11 — Análisis por etapas y continuidad

**Estado**: implementado (ADR 0024 a 0031, 9 de noviembre de 2026).

Hasta aquí cada análisis empezaba en `U = 0` y terminaba en sí mismo. Esta fase encadena análisis y hace legible cómo terminó cada uno:

- **Masas de nudo y grados de libertad sin masa** (ADR 0024): la masa de nudo es un dato del dominio; el modal condensa de forma exacta los grados de libertad sin masa y los transitorios implícitos los admiten con aceleración inicial consistente.
- **Terminación y registro** (ADR 0025): `termination` legible en todos los solvers incrementales y en el tiempo, excepciones tipadas, parada en un valor objetivo de la magnitud controlada (`control_target`) y una política de niveles del registro.
- **Etapas de carga** (ADR 0026): `solve(F, u0=…, constant_load=…)`, `solidum.run(..., previous=…)` y el bloque YAML `stages`.
- **Cargas con sus cargas equivalentes** (ADR 0027): `LoadVector` combina linealmente las cargas equivalentes por elemento con la carga.
- **Longitud característica del punto** (ADR 0028): viaja en el estado del punto de integración para regularizar el ablandamiento (`ManderConcrete1D`).
- **Continuidad en el tiempo** (ADR 0029), **reacciones dinámicas** (ADR 0030) y **análisis lineales sobre el estado de una etapa** (ADR 0031): un transitorio continúa un estático u otro transitorio; el modal, el armónico, el espectral y el lineal linealizan alrededor del estado de la etapa sin cambiarlo.

```
[Etapa 1: estatico] --previous--> [Etapa 2: estatico | transitorio | modal | armonico | espectral | lineal] --previous--> ...
                         (estado de los elementos + carga constante + cargas equivalentes)
```

```callout Riesgo de arquitectura
Medio, ya cerrado. Una etapa hereda el estado que guardan los elementos, así que cada solver debía arrancar de un estado ajeno; se resolvió con mezclas comunes (`TransientTerminationMixin`, `LinearizedStateMixin`) y con comprobación de equilibrio al empezar, sin tocar el contrato de elemento ni de material.
```

## Fase 12 — Interfaces, refuerzo embebido y restricciones generadas

**Estado**: implementado (ADR 0032 y 0033, 10 de noviembre de 2026).

- **Restricciones generadas** (ADR 0032): una familia de generadores (`ConstraintRegistry`, bloque YAML `constraints`) crea sobre el dominio mallado los nodos, elementos y MPC que su geometría pide; la localización de puntos invierte el mapeo isoparamétrico con lo que todo sólido declara. Primer generador: el refuerzo embebido, con barras de dos o tres nodos (`Truss2D3`, `Truss3D3`).
- **Materiales de interfaz** (ADR 0033): una quinta familia, `InterfaceMaterialRegistry` (`interface_materials`), de leyes tracción–desplazamiento relativo para superficies (grietas) y líneas (barras embebidas). Con `bond`, el refuerzo genera interfaces de adherencia-deslizamiento no registradas; la familia cohesiva del módulo `discontinuities` volvió al programa principal y se amplió con la ley de modo mixto.
- En el módulo `discontinuities`, el seguimiento de grieta (ADR 0010, fase E).

```
[Dominio mallado] -- generador de restricciones --> [nodos + elementos de barra + MPC | interfaces de adherencia]
[Materiales de interfaz] --> [interfaces de adherencia (linea)] + [discontinuidades embebidas (superficie)]
```

```callout Riesgo de arquitectura
Bajo. Los generadores producen objetos que el resto del programa ya sabía tratar (nodos, elementos, MPC); las interfaces generadas no se registran ni se exportan a VTK, y sólo exigieron dos piezas genéricas (`VTK_EXPORT` en el exportador y `YAML_KEYS` en el lector).
```

## Fase 13 — Fisuración distribuida del concreto

**Estado**: implementado (ampliación del ADR 0028 y ADR 0034, 11 de noviembre de 2026).

- **Banda de grieta** (`IsotropicDamageCrackBand2D/3D`): el daño isótropo regularizado por energía de fractura, con la medida de Rankine y un tope de daño propio. Los diez sólidos 2D y 3D siembran la longitud característica de cada punto (ampliación del ADR 0028, `CHARACTERISTIC_LENGTH_FACTOR`); un elemento que no la proporciona rechaza el material al construirse.
- **IMPL-EX y pseudo-tiempo** (ADR 0034): la integración implícita-explícita del daño, con tangente secante simétrica y definida positiva, y la pieza genérica por la que los solvers incrementales informan el incremento de su parámetro de control.
- Es el primer modelo de la familia de modelos de concreto; los siguientes (daño de tracción y compresión, grieta rotante, daño-plasticidad) reutilizarán la banda de grieta, la longitud característica y el pseudo-tiempo.

```
[Solver incremental] -- Δτ del intento --> [Ensamblador] -- r = Δτ(n+1)/Δτ(n) --> [Material IMPL-EX | batch_params]
[Elemento sólido] -- h por punto de Gauss --> [estado del punto] --> [Material con banda de grieta]
```

```callout Riesgo de arquitectura
Bajo. La longitud característica reutiliza la pieza del ADR 0028; el pseudo-tiempo es un dato más que el ensamblador pasa a los materiales que lo piden, sin cambiar la firma de `compute_state` ni la del kernel compilado. Un solver nuevo que avance por pasos debe declarar su pseudo-tiempo para admitir materiales IMPL-EX.
```

## Fase 14 — Mallas de Gmsh completas y mallado desde la geometría

**Estado**: implementado (ADR 0035, 11 de noviembre de 2026).

- **Lector de Gmsh** para todas las celdas de los sólidos 2D y 3D, de primer y segundo orden, con el orden de los nodos deducido de las coordenadas de referencia y la orientación corregida.
- **Mallado desde la geometría** en el YAML: el `.geo` con sus parámetros (tamaño de elemento), el orden y la forma de los elementos; Solidum llama a Gmsh y lee la malla. Un estudio de malla es cambiar un dato del YAML.
- **Barras incrustadas** desde un grupo físico de líneas con `element`.
- **Salida de resultados en el paquete** (ADR 0036): `python -m solidum modelo.yaml` escribe la historia de la respuesta en CSV y los VTK en los instantes elegidos.
- **Figuras de los resultados** (ADR 0039): `output.plots` y `python -m solidum.post` dibujan la malla, las deformadas y las curvas con los archivos de la salida.
- **Energías y magnitudes de elemento en la historia** (ADR 0038): trabajo externo, energía elástica y disipada; esfuerzo, deformación, variables internas y fuerzas internas de marcos por elemento o conjunto de elementos.
- **Mallas estructuradas por bloques y calidad de la malla** (ADR 0037): tramos con gradación, celdas excluidas, diagonales de los triángulos y fusión de bloques, desde el YAML con grupos por lado y por condición; resumen y avisos de calidad al leer o generar una malla.

```
[YAML: mesh.geometry + parameters] --> [Gmsh (opcional)] --> [.msh] --> [GmshParser: celda → REFERENCE_NODES del elemento] --> [Domain]
```

```callout Riesgo de arquitectura
Bajo. El lector reutiliza el contrato de la celda de referencia del ADR 0022 y la receta por grupo del ADR 0020; las mallas de primer orden existentes se leen igual.
```

## Mantenimiento de este capítulo

Cada vez que se cierre la implementación de una fase y se acepte su ADR de apertura, esta sección se actualiza:

1. La fase pasa de "[PENDIENTE: ...]" a "implementado".
2. La hipótesis inicial sobre la arquitectura se sustituye por la descripción real, ajustada a las decisiones que se hayan tomado durante la implementación.
3. El diagrama de capas se redibuja conforme a la topología efectiva.
4. Las consecuencias no previstas (ajustes en capas anteriores, contratos modificados) se registran en una subsección "Lecciones de la fase N", de cara a fases posteriores.

De este modo, el capítulo se convierte en una crónica ordenada del crecimiento del programa: en cualquier momento futuro, abrir el manual permite reconstruir no solo el estado actual de Solidum FEM sino la trayectoria completa que lo ha llevado hasta él.
