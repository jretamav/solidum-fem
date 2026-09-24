# Conceptos clave — test de no-caja-negra

> Inventario de los conceptos sostenedores del sistema. **Test operativo**: el usuario debe poder explicar cualquiera de estas entradas en términos generales en cualquier momento. Si una se vuelve opaca, pedir refresco.
>
> Cada entrada: 2–3 frases conceptuales. Para detalles → código.

---

## Plumbing de software (idioms Python que estructuran el proyecto)

### 1. Registry + decorador `@register`
Diccionario global por categoría (materiales, materiales térmicos, elementos, solvers) que mapea nombre → clase. Cada clase, al definirse, se autorregistra mediante un decorador que la añade al diccionario sin tocar listas centrales. Permite que el `YamlParser` instancie cualquier clase por su nombre sin conocerla en compile-time. Un registro que declara su sección YAML (`YAML_SECTION`) es una **familia de material** (`Registry.families()`, ADR 0020): el parser recorre todas las familias con el mismo código.

### 2. Autodiscover (`solidum/autodiscover.py`)
Al hacer `import solidum`, un solo recorrido (`pkgutil.iter_modules`) importa todos los módulos de `solidum/materials`, `solidum/elements`, `solidum/math`. Como los decoradores se ejecutan al importar el módulo, esto basta para poblar los Registry. Equivalente conceptual al `INCLUDE` automático en Fortran moderno, pero en runtime. No recorre `solidum/user/`: los módulos de usuario sólo se cargan cuando un modelo los pide (entrada 10d).

### 3. Contrato declarativo vía `ClassVar` (`STRAIN_DIM`, `DOF_NAMES`, `N_INTEGRATION_POINTS`, `PRIMARY_STATE_VAR`)
Cada elemento y material declara como atributos de clase qué espera y qué produce. La base abstracta los lee y se autoconfigura: registra DOFs, valida compatibilidad material↔elemento (mismo `STRAIN_DIM`), inicializa el `ElementState` con la forma correcta. Sustituye métodos `setup()` repetitivos en cada subclase.

### 4. Validación temprana en construcción
Al instanciar un elemento, el `__init__` de la base verifica que el material sea dimensionalmente compatible (un `Truss2D` con `Elastic2D` falla aquí, no a 200 iteraciones después con un crash críptico). Coste: una comprobación; beneficio: errores físicos imposibles de cometer sin verlos al construir el caso.

### 5. Generic `YamlParser` por introspección de `kwargs`
El parser no contiene un `if material_type == "Elastic1D":` por cada material. Mira el `__init__` de la clase del Registry, extrae sus kwargs, y pasa lo que el YAML provee. Añadir un material nuevo no toca el parser nunca. Tampoco una familia de material nueva ni una referencia nueva entre objetos: la sección YAML la declara el registro, y qué parámetros de un elemento son el id de un objeto de otra familia lo declara el elemento (`REFERENCE_KWARGS`, p. ej. `{"material": ThermalMaterialRegistry}` en los térmicos). Una sección de primer nivel desconocida es error, no se ignora.

### 6. Skill `/solidum-new` (`.claude/skills/solidum-new/SKILL.md`)
Skill versionada con el repo que la IA invoca cuando el usuario pide un material/elemento/solver nuevo. Genera el archivo en su carpeta canónica, con el decorador correcto y un test esqueleto. Cierra el ciclo: la arquitectura optimizada para extensión + la herramienta que materializa la extensión.

### 6b. Ensamblaje por lotes: familias y un solo kernel compilado (ADR 0014)
Los elementos que comparten clase, instancia de material, cuadratura y número de nodos forman una **familia**; cada familia se evalúa entera dentro de un único kernel Numba —en serie o repartido entre hilos— en vez de elemento por elemento, y su estado interno vive en arreglos. El resultado es idéntico al del camino por elemento, que sigue siendo el contrato obligatorio; un componente que no declara kernel simplemente sigue ese camino. Medido: hasta ×110 en el ensamblaje.

### 6c. Corrector de Newton compartido (ADR 0015)
Los cinco solvers iterativos no escriben su bucle de Newton: lo ejecuta `NewtonCorrector`, que posee el backend algebraico, el Newton modificado, el line search y el diagnóstico de divergencia. Cada solver aporta sólo su física por paso (residuo, sistema tangente, actualización) en un objeto `NewtonProblem`, y conserva su control de paso. Un solver nuevo hereda todo lo demás.

---

## Zona gris (clases base y semántica que el usuario debe reconocer)

### 7. `Element` base con boilerplate compartido
Toda la lógica común (registro de DOFs, validación de material, creación de `ElementState`, ensamblaje local→global, extracción de desplazamientos por nodo) vive en la base. Las subclases solo aportan **lo que es físicamente específico**: matriz B, integración, ley de comportamiento llamada al material.

### 8. `Material` base + `STRAIN_DIM`
La base define la interfaz `compute_stress_and_tangent(strain, state)`. El `STRAIN_DIM` (1 / 3 / 6) declara si la deformación de entrada es escalar, vector Voigt 2D o vector Voigt 3D. Es el contrato que permite enchufar materiales sin que el elemento conozca su interior.

### 9. `ElementState` con semántica trial/commit
Cada elemento tiene un objeto `ElementState` con dos copias de las variables internas: **trial** (en exploración durante la iteración) y **committed** (último paso convergido). El solver llama `commit_state()` solo cuando el paso converge — esto evita corromper el historial plástico/de daño con tentativas no convergidas.

### 10. `PRIMARY_STATE_VAR` (zona gris materiales ↔ exporter)
Cada material declara cuál de sus variables internas es la "principal" para visualización (`'damage'`, `'alpha'`, etc.). El `VtkExporter` la lee genéricamente sin saber de qué material proviene. Permite añadir materiales nuevos con visualización automática.

### 10b. `STATE_SCHEMA`: las variables internas declaradas (ADR 0014)
Cada material declara sus variables internas como `{nombre: forma}` (p. ej. `{"eps_p": (4,), "alpha": ()}`; `{}` si no tiene memoria). Con esa declaración el estado se guarda en arreglos por familia y el commit deja de copiar diccionarios. Es obligatorio en todo material del catálogo: lo exige el barrido de contratos.

### 10c. Familias paralelas de materiales
No todo "material" relaciona esfuerzo con deformación. Los **térmicos** (`ThermalMaterial`, Etapa 8) relacionan flujo de calor con gradiente de temperatura (declaran `FLUX_DIM`, no `STRAIN_DIM`, y no usan notación de Voigt). Cada familia tiene clase base, registro y sección YAML propios, para que el parser y los elementos no tengan que distinguir el tipo en cada uso. El programa principal tiene dos, mecánica y térmica; un módulo de usuario puede declarar la suya, como la **cohesiva** del módulo `discontinuities` (`CohesiveMaterial`, ADR 0010), que relaciona tracción con salto de desplazamiento.

### 10d. Programa principal y módulos de usuario (ADR 0020)
Dos términos, los de FEAP. El **programa principal** es Solidum estándar: elementos, materiales y solvers de elementos finitos clásicos. Un **módulo de usuario** es una formulación no estándar en su propia carpeta, `solidum/user/<nombre>/` (hoy `discontinuities`, la discontinuidad embebida), que el programa principal no importa ni nombra: el modelo la pide con `user_modules: [<nombre>]` en el YAML o `solidum.load_user_module("<nombre>")`. El módulo se conecta por piezas genéricas —registros y familias, `REFERENCE_KWARGS`, el gancho de inicio de paso `prepare_step` que los solvers llaman una vez por paso con el estado convergido—, y cambiar esas piezas afecta a todo módulo de usuario.

---

## Núcleo numérico (territorio del usuario, listado para completitud)

### 11. Numeración de DOFs en `Domain`
`Domain` recorre nodos y elementos: cada elemento declara qué DOFs necesita en cada nodo (vía `DOF_NAMES`); el dominio asigna índices globales. La numeración resultante define el tamaño y la topología de la matriz global.

### 12. Ensamblaje sparse con cache COO
El primer ensamblaje calcula los pares (i, j) de cada contribución elemental y los guarda. En las iteraciones siguientes, solo se reescriben los `data` sobre la misma topología — sin recomputar índices. Análogo al patrón "compute once, reuse" típico en FEM compilado, pero implementado sobre `scipy.sparse`.

### 13. Eliminación directa de Dirichlet (ADR 0004)
Las condiciones de frontera Dirichlet se imponen por eliminación directa del sistema, no por penalización. Toda restricción se expresa como afín `u_s = g_s + Σ α_si·u_mi` y se acumula en un `ConstraintSet` (`solidum/bc/constraints.py`). El `Assembler` produce el par `(T, g)` tal que `u = T·u_libre + g`; el solver ve siempre el sistema reducido `K_red·u_libre = F_red` con `K_red = TᵀKT`, `F_red = Tᵀ(F − K·g)`. La imposición es exacta a redondeo (no aproximada como en penalización), preserva la simetría y positividad de `K`, y abre la puerta a MPC lineales sin tocar el solver. Las reacciones se calculan post-hoc como `R = F_int(U) − F_applied` evaluado en los DOFs prescritos.

### 14. Criterio de convergencia dual
Newton-Raphson termina cuando **ambos** criterios están bajo tolerancia: norma del incremento de desplazamientos (relativa) y norma del residuo de fuerza (relativa). Cualquiera por separado puede dar falsos positivos en problemas con softening o con cargas casi nulas.

### 15. `ArcLengthSolver` (Crisfield)
Cuando el problema tiene snap-back / snap-through, Newton-Raphson con control de carga falla. Arc-length añade una incógnita (factor de carga λ) y una restricción geométrica sobre la trayectoria en el espacio (U, λ), permitiendo recorrer ramas con derivada infinita o negativa. Vive en `solidum/math/solvers/arclength.py`.

### 16. Cuadraturas de Gauss centralizadas (`solidum/math/integration.py`)
Tablas y reglas de cuadratura Gauss-Legendre 1D/2D/3D centralizadas. Cada elemento declara `N_INTEGRATION_POINTS` y consume estos puntos/pesos — sin duplicar tablas en cada subclase.

### 17. Return mapping (J2 plasticity, Von Mises 2D)
Algoritmo predictor-corrector: se asume paso elástico, si la tensión predictora viola el criterio de fluencia se proyecta de vuelta a la superficie de fluencia integrando el flujo plástico. La condición de consistencia (Kuhn-Tucker) se resuelve localmente en cada punto de Gauss; el módulo tangente consistente se deriva linealizando este algoritmo (no es el tangente continuo).

### 18. Formulación corotacional (Updated Lagrangian en `Frame2D`, `Truss`)
Se separa el movimiento de cuerpo rígido (rotación de la cuerda del elemento) del estiramiento/flexión local. La ley constitutiva trabaja en el sistema corrotado donde las deformaciones son pequeñas; las rotaciones grandes se manejan por la cinemática del marco. Permite usar materiales lineales con grandes desplazamientos.

### 19. Numba JIT en hot loops
Funciones críticas (ensamblaje elemento→global, return mapping interior) decoradas con `@njit` se compilan a código nativo en la primera llamada. La primera ejecución paga el coste de compilación; las siguientes corren a velocidad cercana a Fortran. Restricción: solo tipos numéricos primitivos y arrays NumPy — no objetos Python.

### 20. Capa algebraica vs. solver de análisis (ADR 0003)
Hay **dos capas de "solver"** y conviene no confundirlas:
- **Solver de análisis** (los 14 del catálogo: `LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`, `DissipationArcLengthSolver`, `IndirectDisplacementSolver`, `ModalSolver`, `NewmarkSolver`, `HHTSolver`, `NewtonNewmarkSolver`, `NewtonHHTSolver`, `CentralDifferenceSolver`, `HarmonicSolver`, `ResponseSpectrumSolver`, `ThetaMethodSolver`): orquesta la estrategia de paso, las iteraciones de Newton, los criterios de convergencia, la longitud de arco, la integración temporal o el barrido en frecuencia.
- **Capa algebraica** (`solidum/math/linalg/`): resuelve el sistema lineal `K·δU = R` (o `Z(ω)·û = F̂` en complejos, o `K·φ = ω²M·φ` en autovalor) que aparece dentro de cada iteración del solver de análisis. Tiene varios backends (Cholesky, Pardiso, LU, el iterativo CG/MINRES, ARPACK para autovalores) y un **despachador interno** que elige el adecuado según las propiedades del operador (simétrica, positiva definida, …).

El usuario solo ve la primera capa; la segunda es plumbing automático. Solo se expone el campo opcional `linear_algebra` en YAML como herramienta de diagnóstico — no como decisión de modelado.

### 21. Solver directo frente a solver iterativo (ADR 0017, 0018)
Un solver **directo** factoriza la matriz y resuelve por sustitución: robusto ante cualquier condicionamiento, pero sus factores tienen muchos más no nulos que la matriz (*relleno*), así que su tiempo y su memoria crecen más deprisa que el modelo. Un solver **iterativo** aproxima la solución hasta una tolerancia sin factorizar: memoria proporcional al modelo, pero sensible al condicionamiento. Solidum usa siempre un directo por defecto (Pardiso multihilo si está instalado) y ofrece el iterativo a petición para modelos grandes, igual que ANSYS y Abaqus.

### 22. Modos de cuerpo rígido derivados del nombre de los DOF (ADR 0018, 0019)
Los movimientos que un modelo sin apoyos puede hacer sin deformarse —traslaciones y giros, y el valor constante de un campo escalar— se calculan a partir de `DOF_NAMES` y las coordenadas, sin saber qué elementos hay. Son el núcleo exacto de `K` sin apoyos. Los usa el multimalla del solver iterativo (sin ellos converge peor que sin precondicionar) y la detección de mecanismos.

### 23. Red de seguridad del análisis estático (ADR 0019)
Un solver directo ante una matriz singular devuelve resultados absurdos sin avisar. Por eso el análisis estático comprueba, antes de resolver, que los apoyos impiden todo movimiento de sólido rígido —y si no, dice cuál queda libre en términos del modelo—, y después, en el caso lineal, que la solución está en equilibrio y que la matriz no tiene pivotes nulos. Dentro de un Newton no rechaza nada: cerca de un punto límite resolver un sistema casi singular es legítimo.

---

## Mantenimiento de este documento

La IA añade una entrada cuando introduce una pieza conceptual nueva (no cuando cambia la implementación de una existente). El usuario solicita refresco si alguna entrada se le vuelve opaca, o pide convertirla en explicación larga vía "explícame X".
