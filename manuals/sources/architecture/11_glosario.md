# Glosario

Términos técnicos centrales de Solidum FEM, ordenados alfabéticamente. Cada entrada da una definición breve (de una a tres frases) y, cuando procede, remite al capítulo donde el concepto se desarrolla.

**Ajuste de paso (incremento de carga)**. División de la carga total en fracciones que se aplican sucesivamente en problemas no lineales. Permite que el método de Newton-Raphson converja a partir de una aproximación cercana a la solución del paso anterior.

**Ardura corrotacional (formulación)**. Descomposición del movimiento de un elemento en una rotación de cuerpo rígido más una deformación local pequeña. La ley constitutiva trabaja en la configuración local; la cinemática del elemento absorbe las grandes rotaciones. Aplicable a armaduras, cables y marcos en régimen geométricamente no lineal con deformación pequeña.

**Architecture Decision Record (ADR)**. Documento breve que registra una decisión de arquitectura junto con su contexto, alternativas consideradas y consecuencias. Reside en `docs/adr/000N-titulo.md`. Ver capítulo 8.

**Autodiscover (descubrimiento automático)**. Procedimiento que recorre las carpetas canónicas del proyecto e importa todos sus módulos al inicio del programa, lo que dispara la ejecución de los decoradores `@register` y, con ello, la inscripción de cada componente en su registro. No recorre `solidum/user/`: los módulos de usuario sólo se cargan cuando un caso los pide. Ver capítulo 5.

**Capa algebraica**. Subsistema en `solidum/math/linalg/` responsable de resolver el sistema lineal `K · x = b` que aparece dentro de cada iteración del solver no lineal. Está separada del solver no lineal y selecciona automáticamente el algoritmo de factorización adecuado mediante un despachador. Ver capítulo 4 y ADR 0003.

**Corrector de Newton**. Bucle de corrección compartido por los solvers iterativos (`NewtonCorrector`, ADR 0015): ensambla una vez en el iterado, evalúa el criterio dual con el residuo de ese ensamblaje y el incremento previo, resuelve con el backend algebraico (degradación a LU, factorización congelada opcional), aplica el line search opcional y comitea al converger. Cada solver le entrega un `NewtonProblem` con su residuo, su sistema tangente y su actualización, y conserva su propio control de paso. Ver capítulos 3 y 4.

**Catálogo**. Documento navegable que enumera todos los componentes implementados de una categoría (elementos, materiales, solvers) con sus parámetros, validez y referencias. Reside en `docs/catalogo_<categoría>.md`. Es índice y resumen; el detalle vive en las especificaciones.

**Contrato declarativo**. Mecanismo por el cual una clase declara, mediante atributos de clase como `STRAIN_DIM`, `DOF_NAMES`, `N_INTEGRATION_POINTS` o `PRIMARY_STATE_VAR`, las propiedades fijas que el sistema necesita conocer de ella. La clase base lee estos atributos y se autoconfigura, sustituyendo métodos de inicialización repetitivos en cada subclase. Ver capítulo 5.

**Despachador algebraico**. Función `select_solver(props, override)` que recibe las propiedades de la matriz de rigidez (simetría, definición positiva) y devuelve el algoritmo de factorización adecuado: Cholesky cuando aplica, LU en caso contrario. Ver ADR 0003.

**`DOF_NAMES`**. Atributo de clase de cada elemento que enumera los nombres de los grados de libertad activos en cada uno de sus nodos. La clase base lo lee para registrar los grados de libertad globales en el `Domain`.

**`Domain`**. Objeto central de la capa de dominio que mantiene la lista de nodos, la conectividad de los elementos, la numeración global de los grados de libertad y, tras la solución, el agregado `SolveResult`.

**Eliminación directa (de condiciones de Dirichlet)**. Estrategia de imposición de condiciones esenciales que retira los grados de libertad prescritos del sistema antes de resolverlo, en oposición a la penalización (que añade un valor grande a la diagonal) y a los multiplicadores de Lagrange (que amplían el sistema). Toda restricción se expresa en forma afín `u_s = g_s + Σ α_si · u_mi` y se acumula en un `ConstraintSet` (`solidum/bc/`). El ensamblador construye un operador disperso `T` y un vector `g` tales que `u = T · u_libre + g`, y entrega al solver el sistema reducido `K_red = TᵀKT`, `F_red = Tᵀ(F − K · g)`. La imposición es exacta a redondeo y preserva la simetría y la positividad definida de la matriz original. Ver ADR 0004.

**Elemento**. Entidad geométrica que construye su matriz de gradientes `B`, su matriz de rigidez tangente y su vector de fuerzas internas a partir de los desplazamientos nodales y del material asignado. Sus contratos declarativos están en `solidum/elements/`. Ver capítulo 4.

**`ElementState`**. Objeto que encapsula las dos copias de las variables internas de cada elemento: una copia *trial* (la que se explora durante las iteraciones del solver) y otra comprometida (la del último paso convergido). Ver capítulo 5.

**Ensamblaje con topología cacheada**. Estrategia de ensamblaje en la que la primera invocación calcula los pares de índices `(i, j)` de cada contribución elemental en formato Coordinate (COO), deriva la estructura CSR de la matriz y un mapa COO → CSR; las invocaciones siguientes reescriben únicamente el vector de datos y lo acumulan sobre la misma estructura, sin recalcular índices ni reordenar. Ver capítulo 5 y ADR 0014.

**Ensamblaje por lotes**. Evaluación de todos los puntos de Gauss de una familia de lote dentro de un único kernel compilado (`solid_family_kernel`, con variante paralela `solid_family_kernel_parallel` bit a bit equivalente), que recibe la cinemática del elemento y la constitutiva del material como funciones tipadas por su firma. El post-proceso por familia (`Family.gauss_state`) usa las mismas funciones. Coexiste con el camino por elemento, que sigue siendo el contrato obligatorio y la referencia física. Ver capítulo 5 y ADR 0014.

**Especificación (spec)**. Documento que describe un componente físico nuevo antes de su implementación: especificación física, formulación numérica, contrato YAML y criterios de aceptación. Reside en `docs/specs/<Nombre>.md`. Es a la vez orden de trabajo y referencia detallada del componente. Ver capítulo 7.

**Estado *trial* y estado comprometido**. Las dos copias de las variables internas de cada elemento. El estado *trial* se explora durante las iteraciones del solver no lineal; al converger un paso, el solver invoca `commit_state()` y promueve el *trial* a estado comprometido. Esta semántica es indispensable en problemas con plasticidad o daño para no contaminar el historial con tentativas posteriormente descartadas.

**Familia de elementos**. Conjunto de elementos que comparten la cinemática esencial: armaduras, cables, marcos (vigas), sólidos 2D, sólidos 3D, cáscaras. Ver capítulo 4.

**Familia de lote**. Conjunto de elementos de un dominio que comparten clase de elemento, instancia de material, regla de cuadratura y número de nodos, y que el ensamblador evalúa juntos en un kernel compilado. Su estado interno vive en un `FamilyState` (arreglos por punto de Gauss) del que `elem.state` es una vista (`BatchedElementState`). Ver capítulo 5 y ADR 0014.

**Familia de material**. Registro que declara la sección de primer nivel del YAML con sus objetos (`YAML_SECTION`): `materials` para los mecánicos, `thermal_materials` para los térmicos y, con el módulo de usuario `discontinuities` cargado, `cohesive_materials`. `Registry.families()` las enumera y el intérprete del YAML las recorre sin conocerlas por nombre; un elemento declara a qué familia pertenecen sus materiales con `REFERENCE_KWARGS`. Ver capítulo 5 y ADR 0020.

**Gradiente conjugado (CG) y MINRES**. Métodos iterativos de Krylov para sistemas lineales simétricos: CG para matrices definidas positivas, MINRES para indefinidas. El solver iterativo de Solidum pasa de uno a otro de forma automática al detectar curvatura negativa. Capítulo 4, ADR 0018.

**Grado de libertad (DOF)**. Cada incógnita escalar del problema discreto. En estática mecánica, los grados de libertad son los desplazamientos nodales y, en su caso, las rotaciones nodales. La numeración global se establece en el `Domain`.

**Longitud de arco (método de Crisfield)**. Estrategia de solución no lineal que introduce el factor de carga como incógnita adicional y añade una restricción geométrica sobre la trayectoria en el espacio (U, λ). Permite recorrer ramas con derivada infinita o negativa, imprescindible para problemas con *snap-back* o *snap-through*.

**Material**. Ley constitutiva que, dada una deformación y un estado interno, devuelve el esfuerzo y el módulo tangente consistente. Su contrato declarativo está en `solidum/materials/`. Ver capítulo 4.

**Mecanismo**. Modelo cuyos apoyos o cuya conectividad permiten un movimiento sin deformación: la matriz de rigidez es singular y el problema estático no tiene solución única. Si el movimiento es de sólido rígido del conjunto, Solidum lo detecta antes de resolver y lo describe; si es interno (una rótula de más), lo detecta tras resolver. Capítulo 5, ADR 0019.

**Método de Elementos Finitos (MEF)**. Técnica de discretización para resolver ecuaciones diferenciales parciales mediante la subdivisión del dominio en elementos sobre los que se interpolan los campos incógnita. Solidum FEM trabaja en aproximación de desplazamientos.

**Modos de cuerpo rígido**. Movimientos de un modelo sin apoyos que no producen deformación: traslaciones y giros infinitesimales, y el valor constante en un campo escalar. En Solidum se derivan del nombre de los grados de libertad y los usan el precondicionador multimalla y la detección de mecanismos. Capítulo 5.

**Módulo de usuario**. Formulación no estándar que vive fuera del programa principal, en su propia carpeta `solidum/user/<nombre>/`, a la manera de los elementos y materiales de usuario de FEAP. El programa principal no la importa ni la nombra: un caso la carga con `user_modules: [<nombre>]` en el YAML o con `solidum.load_user_module` en Python. El primero es `discontinuities` (discontinuidad interior embebida). Ver capítulo 5 y ADR 0020.

**Multimalla algebraico (AMG)**. Precondicionador que resuelve el sistema en una jerarquía de problemas cada vez más pequeños construidos a partir de la propia matriz. Con los modos de cuerpo rígido como casi-núcleo, el número de iteraciones apenas depende del tamaño del modelo. ADR 0018.

**Multipoint constraint (MPC)**. Restricción afín lineal que liga el desplazamiento de un grado de libertad esclavo a una combinación lineal de los desplazamientos de uno o varios grados de libertad maestros, posiblemente con un término independiente. Modelan apoyos en plano oblicuo, periodicidad de celda unitaria, uniones rígidas entre nodos y simetrías no alineadas con los ejes globales. Se declaran mediante `Domain.add_linear_constraint` o desde el bloque `linear_constraints` del archivo YAML, y se imponen por eliminación directa con la misma maquinaria que las condiciones de Dirichlet. Ver capítulo 5 y ADR 0004.

**Notación de Voigt**. Representación de un tensor simétrico de segundo orden como un vector. En 2D, los tres componentes son las dos componentes normales y la componente de cortante (Voigt-3); en 3D, los seis componentes son las tres normales y las tres de cortante (Voigt-6).

**Numba (compilación Just-In-Time, JIT)**. Decorador `@njit` que compila a código nativo, en la primera invocación, las funciones a las que se aplica. La primera ejecución asume el coste de compilación; las siguientes corren a velocidad cercana a Fortran. Restringido a tipos primitivos y arreglos NumPy.

**Pivote nulo**. Valor de la diagonal de la factorización numéricamente igual a cero: indica una matriz singular. Solidum lo considera nulo por debajo de 10⁻¹³ veces la mayor entrada de la matriz, el mismo criterio con que MKL Pardiso los perturba. ADR 0019.

**`PRIMARY_STATE_VAR`**. Atributo de clase de cada material que declara cuál de sus variables internas es la principal a efectos de visualización. El exportador VTK la lee genéricamente sin conocimiento del material que la origina.

**Programa principal**. Solidum estándar: los elementos, materiales y solvers de elementos finitos clásicos, es decir, todo `solidum/` salvo `solidum/user/`. No conoce los módulos de usuario; un test de pureza (`tests/test_main_program_purity.py`) lo vigila. Ver capítulo 5 y ADR 0020.

**Registro (registry)**. Diccionario global por categoría (`MaterialRegistry`, `ThermalMaterialRegistry`, `ElementRegistry`, `SolverRegistry`) que mapea el nombre de cada componente con su clase. Todos son subclases de la clase genérica `Registry`, que puede declarar además metadatos de familia (`YAML_SECTION`, `YAML_LABEL`) y el tipo de spec de sus componentes (`SPEC_KIND`). Cada clase se inscribe automáticamente en su registro mediante un decorador. Ver capítulo 5.

**Regla de la mano derecha (RHR)**. Convención de orientación tridimensional para ejes y vectores momento, adoptada universalmente en Solidum FEM para convenciones de signos de magnitudes vectoriales en 3D. Ver capítulo 6.

**Relleno (fill-in)**. Entradas no nulas que aparecen en los factores de una factorización directa donde la matriz tenía ceros. Crece con el tamaño del modelo y domina el tiempo y la memoria de un solver directo: de ×10 a ×42 la matriz original entre 4 000 y 53 000 grados de libertad. ADR 0017.

**Retorno radial (return mapping)**. Algoritmo predictor-corrector para integrar la ecuación constitutiva de plasticidad J2: se asume un paso elástico (predictor), se evalúa el criterio de fluencia y, si el esfuerzo predictor lo viola, se proyecta de vuelta a la superficie de fluencia (corrector). El módulo tangente consistente se obtiene linealizando el algoritmo discreto.

**Solver directo / solver iterativo**. El directo factoriza la matriz y resuelve por sustitución: robusto ante cualquier condicionamiento, pero su memoria crece más deprisa que el modelo. El iterativo aproxima la solución hasta una tolerancia sin factorizar: memoria proporcional al modelo, pero sensible al condicionamiento. Capítulo 4, ADR 0017 y 0018.

**`SolveResult`**. Agregado inmutable que el `Domain` construye al final de la solución. Contiene los desplazamientos globales, las cargas aplicadas, las reacciones y las fuerzas internas (N, V, M, T) por elemento. Es la interfaz pública para consumidores externos. Ver ADR 0002.

**Solver no lineal**. Componente que orquesta la solución de un problema no lineal: subdivisión en pasos, iteraciones internas, criterio de convergencia. En Solidum FEM están implementados `LinearSolver`, `NonlinearSolver` y `ArcLengthSolver`. Es nivel estratégico, distinto del subsistema algebraico. Ver capítulo 4.

**`STRAIN_DIM`**. Atributo de clase de cada material y elemento que declara la dimensión del vector de deformación que manejan: 1 (escalar axial), 3 (Voigt-3 bidimensional), 6 (Voigt-6 tridimensional). El sistema lo usa para validar la compatibilidad material↔elemento en el momento de construcción.

**Tensión plana / deformación plana**. Hipótesis de modelado para problemas bidimensionales. La tensión plana asume que la tensión perpendicular al plano es nula (placas delgadas); la deformación plana asume que la deformación perpendicular es nula (cuerpos extensos en una dimensión, secciones de tuberías largas).

**Validación temprana en construcción**. Comprobación que la clase base de elemento ejecuta cuando se le asigna un material: verifica que la dimensión de deformación del material coincida con la del elemento. Errores físicos de combinación se detectan en la construcción del caso, no doscientas iteraciones después.

**Visualization Toolkit (VTK)**. Conjunto de formatos de archivo y bibliotecas para representación de mallas y campos discretos en visualizadores de simulación numérica como ParaView. El `VtkExporter` de Solidum FEM produce los archivos correspondientes.

**YAML (YAML Ain't Markup Language)**. Formato de serialización legible por humanos en el que se describen los casos de Solidum FEM (nodos, materiales, elementos, condiciones de contorno, cargas, solver). Es contrato público del programa.
