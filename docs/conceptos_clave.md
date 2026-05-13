# Conceptos clave — test de no-caja-negra

> Inventario de los conceptos sostenedores del sistema. **Test operativo**: el usuario debe poder explicar cualquiera de estas entradas en términos generales en cualquier momento. Si una se vuelve opaca, pedir refresco.
>
> Cada entrada: 2–3 frases conceptuales. Para detalles → código.

---

## Plumbing de software (idioms Python que estructuran el proyecto)

### 1. Registry + decorador `@register`
Diccionario global por categoría (materiales, elementos, solvers) que mapea nombre → clase. Cada clase, al definirse, se autorregistra mediante un decorador que la añade al diccionario sin tocar listas centrales. Permite que el `YamlParser` instancie cualquier clase por su nombre sin conocerla en compile-time.

### 2. Autodiscover (`fenix/autodiscover.py`)
Al hacer `import fenix`, un solo recorrido (`pkgutil.iter_modules`) importa todos los módulos de `fenix/materials`, `fenix/elements`, `fenix/math`. Como los decoradores se ejecutan al importar el módulo, esto basta para poblar los Registry. Equivalente conceptual al `INCLUDE` automático en Fortran moderno, pero en runtime.

### 3. Contrato declarativo vía `ClassVar` (`STRAIN_DIM`, `DOF_NAMES`, `N_INTEGRATION_POINTS`, `PRIMARY_STATE_VAR`)
Cada elemento y material declara como atributos de clase qué espera y qué produce. La base abstracta los lee y se autoconfigura: registra DOFs, valida compatibilidad material↔elemento (mismo `STRAIN_DIM`), inicializa el `ElementState` con la forma correcta. Sustituye métodos `setup()` repetitivos en cada subclase.

### 4. Validación temprana en construcción
Al instanciar un elemento, el `__init__` de la base verifica que el material sea dimensionalmente compatible (un `Truss2D` con `Elastic2D` falla aquí, no a 200 iteraciones después con un crash críptico). Coste: una comprobación; beneficio: errores físicos imposibles de cometer sin verlos al construir el caso.

### 5. Generic `YamlParser` por introspección de `kwargs`
El parser no contiene un `if material_type == "Elastic1D":` por cada material. Mira el `__init__` de la clase del Registry, extrae sus kwargs, y pasa lo que el YAML provee. Añadir un material nuevo no toca el parser nunca.

### 6. Skill `/fenix-new` (`.claude/skills/fenix-new/SKILL.md`)
Skill versionada con el repo que la IA invoca cuando el usuario pide un material/elemento/solver nuevo. Genera el archivo en su carpeta canónica, con el decorador correcto y un test esqueleto. Cierra el ciclo: la arquitectura optimizada para extensión + la herramienta que materializa la extensión.

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

---

## Núcleo numérico (territorio del usuario, listado para completitud)

### 11. Numeración de DOFs en `Domain`
`Domain` recorre nodos y elementos: cada elemento declara qué DOFs necesita en cada nodo (vía `DOF_NAMES`); el dominio asigna índices globales. La numeración resultante define el tamaño y la topología de la matriz global.

### 12. Ensamblaje sparse con cache COO
El primer ensamblaje calcula los pares (i, j) de cada contribución elemental y los guarda. En las iteraciones siguientes, solo se reescriben los `data` sobre la misma topología — sin recomputar índices. Análogo al patrón "compute once, reuse" típico en FEM compilado, pero implementado sobre `scipy.sparse`.

### 13. Eliminación directa de Dirichlet (ADR 0004)
Las condiciones de frontera Dirichlet se imponen por eliminación directa del sistema, no por penalización. Toda restricción se expresa como afín `u_s = g_s + Σ α_si·u_mi` y se acumula en un `ConstraintSet` (`fenix/bc/constraints.py`). El `Assembler` produce el par `(T, g)` tal que `u = T·u_libre + g`; el solver ve siempre el sistema reducido `K_red·u_libre = F_red` con `K_red = TᵀKT`, `F_red = Tᵀ(F − K·g)`. La imposición es exacta a redondeo (no aproximada como en penalización), preserva la simetría y positividad de `K`, y abre la puerta a MPC lineales sin tocar el solver. Las reacciones se calculan post-hoc como `R = F_int(U) − F_applied` evaluado en los DOFs prescritos.

### 14. Criterio de convergencia dual
Newton-Raphson termina cuando **ambos** criterios están bajo tolerancia: norma del incremento de desplazamientos (relativa) y norma del residuo de fuerza (relativa). Cualquiera por separado puede dar falsos positivos en problemas con softening o con cargas casi nulas.

### 15. `ArcLengthSolver` (Crisfield)
Cuando el problema tiene snap-back / snap-through, Newton-Raphson con control de carga falla. Arc-length añade una incógnita (factor de carga λ) y una restricción geométrica sobre la trayectoria en el espacio (U, λ), permitiendo recorrer ramas con derivada infinita o negativa. Vive en `fenix/math/solvers/arclength.py`.

### 16. Cuadraturas de Gauss centralizadas (`fenix/math/integration.py`)
Tablas y reglas de cuadratura Gauss-Legendre 1D/2D/3D centralizadas. Cada elemento declara `N_INTEGRATION_POINTS` y consume estos puntos/pesos — sin duplicar tablas en cada subclase.

### 17. Return mapping (J2 plasticity, Von Mises 2D)
Algoritmo predictor-corrector: se asume paso elástico, si la tensión predictora viola el criterio de fluencia se proyecta de vuelta a la superficie de fluencia integrando el flujo plástico. La condición de consistencia (Kuhn-Tucker) se resuelve localmente en cada punto de Gauss; el módulo tangente consistente se deriva linealizando este algoritmo (no es el tangente continuo).

### 18. Formulación corotacional (Updated Lagrangian en `Frame2D`, `Truss`)
Se separa el movimiento de cuerpo rígido (rotación de la cuerda del elemento) del estiramiento/flexión local. La ley constitutiva trabaja en el sistema corrotado donde las deformaciones son pequeñas; las rotaciones grandes se manejan por la cinemática del marco. Permite usar materiales lineales con grandes desplazamientos.

### 19. Numba JIT en hot loops
Funciones críticas (ensamblaje elemento→global, return mapping interior) decoradas con `@njit` se compilan a código nativo en la primera llamada. La primera ejecución paga el coste de compilación; las siguientes corren a velocidad cercana a Fortran. Restricción: solo tipos numéricos primitivos y arrays NumPy — no objetos Python.

### 20. Capa algebraica vs. solver no lineal (ADR 0003)
Hay **dos capas de "solver"** y conviene no confundirlas:
- **Solver no lineal** (`LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`): orquesta la estrategia de paso, iteraciones de Newton, criterios de convergencia, longitud de arco. Lo que el usuario elige en el YAML con `solver.type`.
- **Capa algebraica** (`fenix/math/linalg/`): resuelve el sistema lineal `K·δU = R` que aparece dentro de cada iteración del solver no lineal. Tiene varios backends (Cholesky, LU, …) y un **despachador interno** que elige el adecuado según las propiedades de `K` (simétrica, positiva definida, …).

El usuario solo ve la primera capa; la segunda es plumbing automático. Solo se expone el campo opcional `linear_algebra` en YAML como herramienta de diagnóstico — no como decisión de modelado.

---

## Mantenimiento de este documento

La IA añade una entrada cuando introduce una pieza conceptual nueva (no cuando cambia la implementación de una existente). El usuario solicita refresco si alguna entrada se le vuelve opaca, o pide convertirla en explicación larga vía "explícame X".
