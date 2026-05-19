# Bloques funcionales

Este capítulo describe los grandes bloques del programa desde el punto de vista funcional: la responsabilidad de cada uno, los componentes disponibles en la versión actual y el espacio de extensión natural. No se trata de un manual de referencia: para la formulación detallada de cada componente debe consultarse `manuals/Reference_manual.pdf` y los catálogos en `docs/`.

A lo largo del capítulo se utilizan las denominaciones abreviadas habituales para las dimensiones del espacio físico: una dimensión (1D), dos dimensiones (2D) y tres dimensiones (3D). En el contexto de los materiales, la **notación de Voigt** representa el tensor simétrico de deformación o tensión como un vector: con tres componentes para el caso bidimensional (Voigt-3) y con seis componentes para el caso tridimensional (Voigt-6). El atributo `STRAIN_DIM` de cada material indica el tamaño de ese vector y, por tanto, la dimensión geométrica del problema sobre el que opera.

## Materiales

Un material en Fenix FEM es una ley constitutiva: dada una deformación (escalar para problemas 1D, en notación de Voigt-3 para problemas 2D, o en notación de Voigt-6 para problemas 3D, según el atributo `STRAIN_DIM`) y un estado interno, devuelve la tensión y el módulo tangente consistente para el solver. El material no posee información sobre el elemento que lo utiliza; el elemento no accede al interior del material. El acoplamiento se establece mediante el contrato declarativo `STRAIN_DIM` y la interfaz `compute_stress_and_tangent(strain, state)`.

Familias actualmente implementadas:

- *Elásticos lineales*: `Elastic1D` para deformación axial escalar; `Elastic2D` para problemas planos (tensión plana o deformación plana) en notación de Voigt 3.
- *Plásticos 1D*: `Elastoplastic1D` con criterio J2 axial, endurecimiento isótropo lineal y tangente algorítmica consistente.
- *Plásticos 2D — J2*: `VonMises2D` con criterio de fluencia J2 y algoritmo de retorno radial (*return mapping*) en notación de Voigt 3. Soporta dos hipótesis cinemáticas mutuamente excluyentes: deformación plana (Simó-Hughes §3.3) y tensión plana (*plane stress projected algorithm*, Simó-Hughes §3.4.1). Tangente algorítmica consistente cerrada en ambas.
- *Plásticos 2D — friccional*: `DruckerPrager2D` con criterio cohesivo-friccional (cono circular suave de Mohr-Coulomb), plasticidad no asociada por defecto y dos ramas cerradas de return (regular sobre el cono y al ápice). Tres calibraciones con Mohr-Coulomb (`plane_strain_matched`, `outer_cone`, `inner_cone`). `IS_SYMMETRIC = False` cuando `ψ ≠ φ`; el despachador algebraico (ADR 0003) degrada Cholesky→LU automáticamente.
- *Daño isótropo*: `IsotropicDamage1D` y `IsotropicDamage2D`, daño escalar con softening exponencial. Tangente algorítmica consistente cerrada en carga activa (negativa en 1D durante softening; asimétrica en 2D porque `(C_e·ε) ⊗ (M·ε)` no es simétrica). Secante en descarga.
- *Cables*: `CableMaterial1D`, ley axial con condición de tracción exclusiva (sin respuesta a compresión).

Espacio natural de extensión: cualquier ley constitutiva que pueda expresarse como una función `(strain, state) → (stress, tangent)` se incorpora sin modificación del resto del programa. Candidatos directos:

- [PENDIENTE: Modelo de Mohr-Coulomb 2D con esquinas (multi-superficie) para mayor fidelidad geotécnica respecto a Drucker-Prager.]
- [PENDIENTE: Endurecimiento cinemático (Bauschinger) en J2 para plasticidad cíclica y fatiga.]
- [PENDIENTE: Hiperelasticidad isótropa (Neo-Hooke, Mooney-Rivlin, Ogden).]
- [PENDIENTE: Viscoplasticidad de Perzyna y Duvaut-Lions.]
- [PENDIENTE: Hiperelasticidad incompresible con tratamiento mixto desplazamiento-presión.]
- [PENDIENTE: Drucker-Prager en tensión plana (la proyección con `σ_zz=0` acoplada a flujo dilatante queda fuera del alcance de la versión actual).]

En todos los casos basta declarar el `STRAIN_DIM` apropiado y el atributo `PRIMARY_STATE_VAR` para visualización.

## Elementos

Un elemento finito en Fenix FEM es una entidad geométrica que construye su matriz de gradientes `B`, su matriz de rigidez tangente y su vector de fuerzas internas a partir de los desplazamientos nodales y del material asociado. Cada elemento declara su `DOF_NAMES`, su `STRAIN_DIM` y el número de puntos de integración mediante `N_INTEGRATION_POINTS`. La clase base concentra la lógica común (registro de DOF, ensamblaje local→global, validación material↔elemento, gestión del `ElementState`) y cada subclase aporta exclusivamente lo físicamente específico (matriz `B`, esquema de integración, invocación al material).

El estado del programa por familia se presenta a continuación en forma de **matrices de cobertura**. Cada celda combina dos ejes ortogonales: el régimen geométrico (linealidad o no linealidad) y el régimen material (lineal o no lineal). La no linealidad material se obtiene mediante la asignación al elemento de un material no lineal del catálogo (plasticidad, daño); ningún elemento es intrínsecamente no lineal en el material. La marca ✓ indica combinación implementada y validada; la marca [Pendiente: ...] indica el hueco como candidato natural de extensión; el guion (—) indica combinación que no aplica físicamente o que no se contempla.

### Familia armadura (truss)

Elemento articulado que transmite exclusivamente esfuerzo axial. Material consumido: 1D escalar (`Elastic1D`, `Elastoplastic1D`, `IsotropicDamage1D`).

**Primer orden (interpolación lineal, dos nodos)**

[TABLA: Cobertura de la familia armadura — primer orden.]
| Régimen | 2D · material lineal | 2D · material no lineal | 3D · material lineal | 3D · material no lineal |
|---|---|---|---|---|
| Linealidad geométrica | ✓ `Truss2D` | ✓ `Truss2D` | ✓ `Truss3D` | ✓ `Truss3D` |
| No linealidad geométrica | ✓ `Truss2DCorot` | ✓ `Truss2DCorot` | ✓ `Truss3DCorot` | ✓ `Truss3DCorot` |

**Segundo orden (interpolación cuadrática, tres nodos)**

[TABLA: Cobertura de la familia armadura — segundo orden.]
| Régimen | 2D · material lineal | 2D · material no lineal | 3D · material lineal | 3D · material no lineal |
|---|---|---|---|---|
| Linealidad geométrica | [PENDIENTE: Armadura 2D de segundo orden con linealidad geométrica.] | [PENDIENTE: Armadura 2D de segundo orden con linealidad geométrica y material no lineal.] | [PENDIENTE: Armadura 3D de segundo orden con linealidad geométrica.] | [PENDIENTE: Armadura 3D de segundo orden con linealidad geométrica y material no lineal.] |
| No linealidad geométrica | [PENDIENTE: Armadura 2D de segundo orden corrotacional.] | [PENDIENTE: Armadura 2D de segundo orden corrotacional con material no lineal.] | [PENDIENTE: Armadura 3D de segundo orden corrotacional.] | [PENDIENTE: Armadura 3D de segundo orden corrotacional con material no lineal.] |

### Familia cable

Elemento de tracción exclusiva en formulación corrotacional. Material consumido: `CableMaterial1D` (sin respuesta a compresión).

**Primer orden (interpolación lineal, dos nodos)**

[TABLA: Cobertura de la familia cable — primer orden.]
| Régimen | 2D · material lineal a tracción | 2D · material no lineal a tracción | 3D · material lineal a tracción | 3D · material no lineal a tracción |
|---|---|---|---|---|
| Linealidad geométrica | — | — | — | — |
| No linealidad geométrica | ✓ `Cable2DCorot` | [PENDIENTE: Cable 2D corrotacional con ley constitutiva no lineal a tracción (relajación, fluencia, daño tensional).] | ✓ `Cable3DCorot` | [PENDIENTE: Cable 3D corrotacional con ley constitutiva no lineal a tracción.] |

La fila "Linealidad geom." se marca como no aplicable: un cable es, por su naturaleza física, un elemento de grandes desplazamientos (la rigidez axial es la única que aporta y depende de la configuración corriente).

**Segundo orden (interpolación cuadrática, tres nodos)**

[TABLA: Cobertura de la familia cable — segundo orden.]
| Régimen | 2D · material lineal | 2D · material no lineal | 3D · material lineal | 3D · material no lineal |
|---|---|---|---|---|
| No linealidad geométrica | [PENDIENTE: Cable 2D de segundo orden corrotacional.] | [PENDIENTE: Cable 2D de segundo orden corrotacional con material no lineal.] | [PENDIENTE: Cable 3D de segundo orden corrotacional.] | [PENDIENTE: Cable 3D de segundo orden corrotacional con material no lineal.] |

### Familia marco / viga

Elemento con flexión y, en su caso, cortante y torsión. Material consumido: 1D escalar. La no linealidad material se introduce mediante el módulo tangente `E_t` del material, que escala uniformemente la sección; **la plasticidad distribuida en la sección (formulación fibra a fibra) no está implementada** y se marca como pendiente independiente.

**Primer orden, formulación Euler-Bernoulli**

[TABLA: Cobertura de la familia marco — primer orden, formulación Euler-Bernoulli.]
| Régimen | 2D · material lineal | 2D · material no lineal (sección uniforme) | 3D · material lineal | 3D · material no lineal (sección uniforme) |
|---|---|---|---|---|
| Linealidad geométrica | ✓ `Frame2DEuler` | ✓ `Frame2DEuler` con `Elastoplastic1D` o `IsotropicDamage1D` | ✓ `Frame3D` | ✓ `Frame3D` con material 1D no lineal |
| No linealidad geométrica | ✓ `Frame2DEulerCorot` | ✓ `Frame2DEulerCorot` con material 1D no lineal | [PENDIENTE: Marco 3D Euler-Bernoulli corrotacional con grandes rotaciones.] | [PENDIENTE: Marco 3D Euler-Bernoulli corrotacional con material no lineal.] |

**Primer orden, formulación Timoshenko (con deformación por cortante)**

[TABLA: Cobertura de la familia marco — primer orden, formulación Timoshenko.]
| Régimen | 2D · material lineal | 2D · material no lineal (sección uniforme) | 3D · material lineal | 3D · material no lineal |
|---|---|---|---|---|
| Linealidad geométrica | ✓ `Frame2DTimoshenko` | ✓ `Frame2DTimoshenko` con material 1D no lineal | [PENDIENTE: Marco 3D Timoshenko con deformación por cortante.] | [PENDIENTE: Marco 3D Timoshenko con material no lineal.] |
| No linealidad geométrica | [PENDIENTE: Marco 2D Timoshenko corrotacional.] | [PENDIENTE: Marco 2D Timoshenko corrotacional con material no lineal.] | [PENDIENTE: Marco 3D Timoshenko corrotacional.] | [PENDIENTE: Marco 3D Timoshenko corrotacional con material no lineal.] |

**Segundo orden y formulaciones avanzadas**

- [PENDIENTE: Marco de segundo orden (interpolación cuadrática) en cualquiera de las formulaciones anteriores.]
- [PENDIENTE: Plasticidad distribuida fibra a fibra en marcos 2D y 3D, con discretización transversal de la sección.]
- [PENDIENTE: Marco con sección variable a lo largo del eje (sección no prismática).]

### Familia sólido 2D

Elemento de medio continuo bidimensional bajo hipótesis de tensión plana o deformación plana. Material consumido: 2D en notación de Voigt 3 (`Elastic2D`, `VonMises2D`, `IsotropicDamage2D`, `DruckerPrager2D`).

**Primer orden**

[TABLA: Cobertura de la familia sólido bidimensional — primer orden.]
| Topología | Material lineal | Material no lineal |
|---|---|---|
| Triangular lineal (3 nodos) | ✓ `Tri3` | ✓ `Tri3` con `VonMises2D`, `IsotropicDamage2D` o `DruckerPrager2D` |
| Cuadrilátero bilineal (4 nodos) | ✓ `Quad4` | ✓ `Quad4` con `VonMises2D`, `IsotropicDamage2D` o `DruckerPrager2D` |

La columna de no linealidad geométrica se omite en esta familia: los elementos sólidos 2D actuales no incluyen formulación corrotacional ni lagrangiana actualizada.

- [PENDIENTE: Sólido 2D con no linealidad geométrica (formulación lagrangiana actualizada o total).]

**Segundo orden**

[TABLA: Cobertura de la familia sólido bidimensional — segundo orden.]
| Topología | Material lineal | Material no lineal |
|---|---|---|
| Triangular cuadrático $P_2$ (6 nodos) | ✓ `Tri6` | ✓ `Tri6` con `VonMises2D`, `IsotropicDamage2D` o `DruckerPrager2D` |
| Cuadrilátero serendípito $Q_2^{\text{seren}}$ (8 nodos) | ✓ `Quad8` | ✓ `Quad8` con cualquier material 2D no lineal |
| Cuadrilátero lagrangiano $Q_2$ (9 nodos) | ✓ `Quad9` | ✓ `Quad9` con cualquier material 2D no lineal |

### Familia sólido 3D

- [PENDIENTE: Tetraedro lineal `Tet4` (cuatro nodos), primer orden.]
- [PENDIENTE: Tetraedro cuadrático `Tet10` (diez nodos), segundo orden.]
- [PENDIENTE: Hexaedro trilineal `Hex8` (ocho nodos), primer orden.]
- [PENDIENTE: Hexaedro serendípito `Hex20` (veinte nodos), segundo orden.]
- [PENDIENTE: Hexaedro lagrangiano `Hex27` (veintisiete nodos), segundo orden.]

### Familias adicionales no implementadas

- [PENDIENTE: Elementos de cáscara y lámina (formulación de Mindlin-Reissner para cáscaras gruesas; formulación de Kirchhoff-Love para cáscaras delgadas).]
- [PENDIENTE: Elementos mixtos desplazamiento-presión (Q1-P0, Taylor-Hood) para problemas casi incompresibles.]
- [PENDIENTE: Elementos de interfaz para modelado de fricción y contacto cohesivo entre superficies.]

### Contrato común y convención de signos

Todos los elementos tipo armadura, cable o marco implementan el contrato `internal_forces(U)` y devuelven los esfuerzos en la convención de signos del proyecto (ADR 0002).

## Solvers no lineales

El solver no lineal corresponde al nivel estratégico del cálculo: orquesta la subdivisión en pasos, las iteraciones internas y los criterios de convergencia. Es el componente que el usuario selecciona en el archivo YAML mediante el campo `solver.type`.

Disponibles en la versión actual:

- `LinearSolver` — un único paso, una única resolución del sistema `K · U = F`. Aplicable a problemas estrictamente lineales (rigidez constante, sin actualización de geometría, sin no linealidad material).
- `NonlinearSolver` — método de Newton-Raphson incremental con control de carga. Aplica fracciones crecientes de la carga total y, en cada paso, itera hasta la convergencia. Criterio de parada dual (desplazamiento y fuerza).
- `ArcLengthSolver` — método de Crisfield. Introduce una incógnita adicional (el factor de carga) y una restricción geométrica sobre la trayectoria en el espacio (U, λ), lo que permite recorrer ramas con derivada infinita o negativa. Imprescindible para problemas con *snap-back* o *snap-through*.
- `ModalSolver` (ADR 0009 fase 1) — análisis modal: problema generalizado de valores y vectores característicos `K · φ = ω² · M · φ`. No es análisis dinámico propiamente dicho (no hay paso de tiempo ni evolución temporal); sus soluciones admiten interpretación como frecuencias propias y formas de vibración libre no amortiguada. Internamente delega en `EigenSolver` (ARPACK Lanczos con shift-invert) sobre el sistema reducido por Dirichlet. Devuelve `ModalResult` con frecuencias, períodos y modos M-ortonormales; expone además `free_vibration(M, u0, u0_dot, t)` para reconstruir analíticamente la respuesta libre por superposición modal.
- `NewmarkSolver` (ADR 0009 fase 3) — análisis dinámico transitorio lineal: integración temporal directa de `M · ü + C · u̇ + K · u = F(t)` por el método de Newmark-β con `(β, γ)` parametrizables (default 1/4, 1/2 — *average acceleration*, incondicionalmente estable, error `O(Δt²)`). Amortiguamiento Rayleigh `C = α · M + β · K` con coeficientes directos o calibrados modalmente a partir de dos pares `(ξ, ω)`. Cargas dependientes del tiempo por callback Python. Devuelve `TransientResult` con historiales `(t, u, u̇, ü)`. La matriz efectiva `A_eff = M + γΔt · C + βΔt² · K` se factoriza una sola vez al inicio y se reutiliza en todos los pasos.
- `NewtonNewmarkSolver` (ADR 0009 fase 4) — análisis dinámico transitorio no lineal. **Variante** de `NewmarkSolver` (Reglas §4): hereda predictores, correctores y reducción de Dirichlet; añade un bucle Newton-Raphson dentro de cada paso temporal sobre el residuo dinámico `R = F_ext − F_int(u) − C · u̇ − M · ü` con jacobiano `J = M + γΔt · C + βΔt² · K_t` (rigidez tangente corriente). Convergencia dual (ADR 0007). Amortiguamiento Rayleigh constante en el tiempo, calibrado con `K_0` (rigidez elástica de referencia; convención estándar Abaqus/ANSYS/OpenSees). Recupera el caso lineal a paridad de bits cuando los materiales no tienen historia activa. Newton modificado opcional (`freeze_tangent_after_iter`, ADR 0003 fase 2).
- `HHTSolver` y `NewtonHHTSolver` (ADR 0009, variantes HHT-α) — disipación numérica controlada en altas frecuencias (Hilber-Hughes-Taylor 1977) con parámetro `α ∈ [−1/3, 0]` y `(β, γ)` auto-derivados. Subclases de `NewmarkSolver` y `NewtonNewmarkSolver` respectivamente (Reglas §4 — variantes).
- `CentralDifferenceSolver` (ADR 0009 fase 5) — integración explícita por diferencias centradas (leapfrog Belytschko-Liu-Moran). `M⁻¹` trivial sobre masa lumped diagonal. Lineal y no lineal en una sola clase con parámetro `nonlinear`. Estabilidad condicional CFL con detección a posteriori de divergencia.
- `HarmonicSolver` (ADR 0009 fase 6) — respuesta forzada armónica en el dominio de la frecuencia. Aritmética compleja con `scipy.sparse.linalg.spsolve` complejo; factorización LU compleja por frecuencia (no se cachea). Barrido lineal/logarítmico/explícito. Devuelve `HarmonicResult` con métodos `.amplitude()` y `.phase()`.
- `ResponseSpectrumSolver` (ADR 0009 fase 7) — análisis sísmico por combinación modal SRSS/CQC contra espectro de respuesta. Orquesta `ModalSolver` interno y delega los algoritmos en `fenix.math.modal_response` (centralización por regla D de auditoría). Devuelve `ResponseSpectrumResult` con respuesta envolvente, factores de participación y masas efectivas.

Espacio natural de extensión:

- [PENDIENTE: Newton-Raphson con búsqueda lineal (line-search) para mejorar la robustez ante no convergencia. **Parcialmente entregado por ADR 0011** — line search disponible como opt-in (`line_search=True`); default `False` por la justificación bibliográfica de Borst & Sluys 1999.]
- [PENDIENTE: Solver cuasi-Newton tipo BFGS para reducir el coste de actualización de la matriz tangente.]
- [PENDIENTE: Generalized-α, Bossak-α — esquemas adicionales de la familia con disipación numérica.]
- [PENDIENTE: Solver de relajación dinámica para búsqueda de configuraciones de equilibrio en estructuras flexibles.]
- [PENDIENTE: Excitación sísmica multi-direccional simultánea (CQC3, regla 100/30/30) y multi-support seismic excitation.]
- [PENDIENTE: Δt adaptativo para transitorios.]

La incorporación de un solver nuevo se realiza en `fenix/math/solvers/<nombre>.py` mediante el decorador `@SolverRegistry.register`. El **dispatch a entrypoints** en `fenix.entry.run_yaml` es declarativo (regla C de auditoría aplicada 2026-05-18): cada solver expone un atributo de clase `PIPELINE_KIND ∈ {"static", "modal", "transient", "harmonic", "spectrum"}` y `run_yaml` despacha por valor, sin tocarse cuando entran solvers nuevos no clásicos.

## Subsistema algebraico

Por debajo del solver no lineal, dentro de cada iteración, se requiere la resolución de un sistema lineal `K · x = b`. Esta es la responsabilidad del subsistema algebraico (`fenix/math/linalg/`, ADR 0003). Constituye un nivel táctico: no decide la estrategia del cálculo, sino que selecciona el algoritmo numérico adecuado para resolver cada sistema lineal.

Algoritmos disponibles:

- `LUSolver` — factorización LU mediante SuperLU. Algoritmo universal: aplicable a cualquier matriz no singular, simétrica o no, definida o no.
- `CholeskySolver` — factorización de Cholesky mediante CHOLMOD. Significativamente más rápido y con menor consumo de memoria, pero aplicable únicamente cuando `K` es simétrica y definida positiva.
- `EigenSolver` (ADR 0009) — algoritmo de autovalor generalizado simétrico para el problema `K · φ = ω² · M · φ`. Envuelve `scipy.sparse.linalg.eigsh` (ARPACK Lanczos con shift-invert centrado en `σ`). Vive en `fenix/math/linalg/eigen.py` y no comparte la interfaz `solve(K, b)` de los anteriores — su firma natural es `solve(K, M, n_modes) → (λ, φ)`.
- [PENDIENTE: `LDLTSolver` — factorización LDLᵀ para el caso simétrico no definido positivo, reservada para la fase 2 del ADR 0003.]
- [PENDIENTE: Algoritmos algebraicos directos paralelos (Pardiso, MUMPS) para problemas de gran tamaño.]
- [PENDIENTE: Algoritmos algebraicos iterativos con precondicionamiento (gradiente conjugado precondicionado, GMRES, multimalla algebraica).]

**Despachador.** La función `select_solver(props)` inspecciona las propiedades de la matriz (encapsuladas en `StiffnessProperties`: simetría, definición positiva) y selecciona el algoritmo adecuado de forma automática. El usuario puede forzar un algoritmo específico desde el archivo YAML mediante el campo opcional `linear_algebra` en calidad de herramienta de diagnóstico, no como decisión de modelado.

**Justificación de la separación en dos niveles.** El solver no lineal y el subsistema algebraico resuelven cuestiones distintas. El primero decide cómo recorrer la respuesta del sistema (incrementos, iteraciones, longitud de arco). El segundo decide con qué método resolver cada sistema lineal individual. Esta separación permite, por ejemplo, mejorar el rendimiento del paso elástico sin afectar al método de Newton-Raphson, o introducir longitud de arco sin reescribir el código del subsistema algebraico.

## Tipos de análisis

El tipo de análisis en Fenix FEM corresponde al problema físico que afronta el solver, no al solver en sí.

Implementados en la versión actual:

- *Estática lineal*. Desplazamientos pequeños, rigidez constante, una resolución directa. Combinación: `LinearSolver` con materiales lineales y elementos en formulación lineal.
- *Estática no lineal incremental*. No linealidad geométrica (corrotacional), material (plasticidad, daño) o ambas. Combinación: `NonlinearSolver` o `ArcLengthSolver` con materiales no lineales y/o elementos corrotacionales.
- *Análisis modal* (ADR 0009 fase 1). Problema algebraico de valores y vectores característicos generalizado `K · φ = ω² · M · φ`. Combinación: `ModalSolver` con materiales que declaren `density` y elementos que implementen `compute_mass_matrix`.
- *Dinámica estructural transitoria lineal* (ADR 0009 fase 3). Integración temporal directa de la ecuación de movimiento por Newmark-β, con amortiguamiento Rayleigh proporcional y cargas dependientes del tiempo por callback Python. Combinación: `NewmarkSolver` con material lineal, elementos en formulación lineal y `density` declarada.
- *Dinámica estructural transitoria no lineal* (ADR 0009 fase 4). Integración Newmark con Newton-Raphson dentro de cada paso sobre el residuo dinámico `R = F_ext − F_int(u) − C·u̇ − M·ü`. Combinación: `NewtonNewmarkSolver` con cualquier material con historia (`Elastoplastic1D`, `VonMises2D` en ambas hipótesis, `DruckerPrager2D`, `IsotropicDamage1D/2D`) y `density` declarada. Habilita sísmica con disipación plástica, impacto en hormigón friccional, fatiga de bajo ciclo y vibraciones de marcos elastoplásticos.
- *Dinámica estructural transitoria con disipación HHT-α* (ADR 0009, variante de fases 3 y 4). Mismas combinaciones que Newmark, con la familia HHT-α que controla la disipación numérica de modos altos espurios.
- *Dinámica estructural transitoria explícita* (ADR 0009 fase 5). Integración explícita por diferencias centradas. Combinación: `CentralDifferenceSolver` con masa lumped (`lumping="lumped"`) y la condición de estabilidad CFL `Δt < L_el · √(ρ/E)` respetada por el usuario. Apropiado para *wave propagation*, impacto y dinámica de respuesta rápida.
- *Respuesta forzada armónica en el dominio de la frecuencia* (ADR 0009 fase 6). Resolución directa del sistema complejo `(−ω²M + iωC + K)·û = F̂` por barrido en `ω`. Combinación: `HarmonicSolver` con amortiguamiento Rayleigh y excitación armónica.
- *Análisis sísmico por combinación modal espectral* (ADR 0009 fase 7). Cálculo modal interno + factores de participación + combinación SRSS o CQC contra un espectro de respuesta tabulado o callable. Combinación: `ResponseSpectrumSolver` con materiales lineales, `density` declarada y un espectro normativo o equivalente.

Previstos y no implementados (dirección del proyecto, sin compromiso de calendario):

- [PENDIENTE: Problema térmico estacionario lineal y no lineal (conductividad dependiente de la temperatura).]
- [PENDIENTE: Problema térmico transitorio con esquemas de integración temporal (Crank-Nicolson, theta-método).]
- [PENDIENTE: Acoplamiento termo-mecánico desacoplado (staggered) para casos en que la influencia mecánica sobre el campo térmico sea despreciable.]
- [PENDIENTE: Acoplamiento termo-mecánico monolítico para problemas con fuerte interacción bidireccional.]
- [PENDIENTE: Excitación sísmica multi-direccional simultánea (CQC3, 100/30/30) y multi-support para diferencias de movimiento entre apoyos.]

La estructura de capas y registros se ha diseñado para incorporar dichos análisis sin reformular las piezas existentes; el ADR correspondiente se redactará al inicio de cada fase de diseño.
