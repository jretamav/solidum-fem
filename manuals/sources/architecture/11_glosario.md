# Glosario

Términos técnicos centrales de Solidum FEM, ordenados alfabéticamente. Cada entrada da una definición breve (de una a tres frases) y, cuando procede, remite al capítulo donde el concepto se desarrolla.

**Ajuste de paso (incremento de carga)**. División de la carga total en fracciones que se aplican sucesivamente en problemas no lineales. Permite que el método de Newton-Raphson converja a partir de una aproximación cercana a la solución del paso anterior.

**Ardura corrotacional (formulación)**. Descomposición del movimiento de un elemento en una rotación de cuerpo rígido más una deformación local pequeña. La ley constitutiva trabaja en la configuración local; la cinemática del elemento absorbe las grandes rotaciones. Aplicable a armaduras, cables y marcos en régimen geométricamente no lineal con deformación pequeña.

**Architecture Decision Record (ADR)**. Documento breve que registra una decisión de arquitectura junto con su contexto, alternativas consideradas y consecuencias. Reside en `docs/adr/000N-titulo.md`. Ver capítulo 8.

**Autodiscover (descubrimiento automático)**. Procedimiento que recorre las carpetas canónicas del proyecto e importa todos sus módulos al inicio del programa, lo que dispara la ejecución de los decoradores `@register` y, con ello, la inscripción de cada componente en su registro. Ver capítulo 5.

**Capa algebraica**. Subsistema en `solidum/math/linalg/` responsable de resolver el sistema lineal `K · x = b` que aparece dentro de cada iteración del solver no lineal. Está separada del solver no lineal y selecciona automáticamente el algoritmo de factorización adecuado mediante un despachador. Ver capítulo 4 y ADR 0003.

**Catálogo**. Documento navegable que enumera todos los componentes implementados de una categoría (elementos, materiales, solvers) con sus parámetros, validez y referencias. Reside en `docs/catalogo_<categoría>.md`. Es índice y resumen; el detalle vive en las especificaciones.

**Contrato declarativo**. Mecanismo por el cual una clase declara, mediante atributos de clase como `STRAIN_DIM`, `DOF_NAMES`, `N_INTEGRATION_POINTS` o `PRIMARY_STATE_VAR`, las propiedades fijas que el sistema necesita conocer de ella. La clase base lee estos atributos y se autoconfigura, sustituyendo métodos de inicialización repetitivos en cada subclase. Ver capítulo 5.

**Despachador algebraico**. Función `select_solver(props, override)` que recibe las propiedades de la matriz de rigidez (simetría, definición positiva) y devuelve el algoritmo de factorización adecuado: Cholesky cuando aplica, LU en caso contrario. Ver ADR 0003.

**`DOF_NAMES`**. Atributo de clase de cada elemento que enumera los nombres de los grados de libertad activos en cada uno de sus nodos. La clase base lo lee para registrar los grados de libertad globales en el `Domain`.

**`Domain`**. Objeto central de la capa de dominio que mantiene la lista de nodos, la conectividad de los elementos, la numeración global de los grados de libertad y, tras la solución, el agregado `SolveResult`.

**Eliminación directa (de condiciones de Dirichlet)**. Estrategia de imposición de condiciones esenciales que retira los grados de libertad prescritos del sistema antes de resolverlo, en oposición a la penalización (que añade un valor grande a la diagonal) y a los multiplicadores de Lagrange (que amplían el sistema). Toda restricción se expresa en forma afín `u_s = g_s + Σ α_si · u_mi` y se acumula en un `ConstraintSet` (`solidum/bc/`). El ensamblador construye un operador disperso `T` y un vector `g` tales que `u = T · u_libre + g`, y entrega al solver el sistema reducido `K_red = TᵀKT`, `F_red = Tᵀ(F − K · g)`. La imposición es exacta a redondeo y preserva la simetría y la positividad definida de la matriz original. Ver ADR 0004.

**Elemento**. Entidad geométrica que construye su matriz de gradientes `B`, su matriz de rigidez tangente y su vector de fuerzas internas a partir de los desplazamientos nodales y del material asignado. Sus contratos declarativos están en `solidum/elements/`. Ver capítulo 4.

**`ElementState`**. Objeto que encapsula las dos copias de las variables internas de cada elemento: una copia *trial* (la que se explora durante las iteraciones del solver) y otra comprometida (la del último paso convergido). Ver capítulo 5.

**Ensamblaje con caché Coordinate (COO)**. Estrategia de ensamblaje en la que la primera invocación calcula los pares de índices `(i, j)` de cada contribución elemental y los almacena en formato Coordinate; las invocaciones siguientes reescriben únicamente el vector de datos sobre la misma topología, sin recalcular índices. Ver capítulo 5.

**Especificación (spec)**. Documento que describe un componente físico nuevo antes de su implementación: especificación física, formulación numérica, contrato YAML y criterios de aceptación. Reside en `docs/specs/<Nombre>.md`. Es a la vez orden de trabajo y referencia detallada del componente. Ver capítulo 7.

**Estado *trial* y estado comprometido**. Las dos copias de las variables internas de cada elemento. El estado *trial* se explora durante las iteraciones del solver no lineal; al converger un paso, el solver invoca `commit_state()` y promueve el *trial* a estado comprometido. Esta semántica es indispensable en problemas con plasticidad o daño para no contaminar el historial con tentativas posteriormente descartadas.

**Familia de elementos**. Conjunto de elementos que comparten la cinemática esencial: armaduras, cables, marcos (vigas), sólidos 2D, sólidos 3D, cáscaras. Ver capítulo 4.

**Grado de libertad (DOF)**. Cada incógnita escalar del problema discreto. En estática mecánica, los grados de libertad son los desplazamientos nodales y, en su caso, las rotaciones nodales. La numeración global se establece en el `Domain`.

**Longitud de arco (método de Crisfield)**. Estrategia de solución no lineal que introduce el factor de carga como incógnita adicional y añade una restricción geométrica sobre la trayectoria en el espacio (U, λ). Permite recorrer ramas con derivada infinita o negativa, imprescindible para problemas con *snap-back* o *snap-through*.

**Material**. Ley constitutiva que, dada una deformación y un estado interno, devuelve la tensión y el módulo tangente consistente. Su contrato declarativo está en `solidum/materials/`. Ver capítulo 4.

**Método de Elementos Finitos (MEF)**. Técnica de discretización para resolver ecuaciones diferenciales parciales mediante la subdivisión del dominio en elementos sobre los que se interpolan los campos incógnita. Solidum FEM trabaja en aproximación de desplazamientos.

**Multipoint constraint (MPC)**. Restricción afín lineal que liga el desplazamiento de un grado de libertad esclavo a una combinación lineal de los desplazamientos de uno o varios grados de libertad maestros, posiblemente con un término independiente. Modelan apoyos en plano oblicuo, periodicidad de celda unitaria, uniones rígidas entre nodos y simetrías no alineadas con los ejes globales. Se declaran mediante `Domain.add_linear_constraint` o desde el bloque `linear_constraints` del archivo YAML, y se imponen por eliminación directa con la misma maquinaria que las condiciones de Dirichlet. Ver capítulo 5 y ADR 0004.

**Notación de Voigt**. Representación de un tensor simétrico de segundo orden como un vector. En 2D, los tres componentes son las dos componentes normales y la componente de cortante (Voigt-3); en 3D, los seis componentes son las tres normales y las tres de cortante (Voigt-6).

**Numba (compilación Just-In-Time, JIT)**. Decorador `@njit` que compila a código nativo, en la primera invocación, las funciones a las que se aplica. La primera ejecución asume el coste de compilación; las siguientes corren a velocidad cercana a Fortran. Restringido a tipos primitivos y arreglos NumPy.

**`PRIMARY_STATE_VAR`**. Atributo de clase de cada material que declara cuál de sus variables internas es la principal a efectos de visualización. El exportador VTK la lee genéricamente sin conocimiento del material que la origina.

**Registro (registry)**. Diccionario global por categoría (`MaterialRegistry`, `ElementRegistry`, `SolverRegistry`) que mapea el nombre de cada componente con su clase. Cada clase se inscribe automáticamente en su registro mediante un decorador. Ver capítulo 5.

**Regla de la mano derecha (RHR)**. Convención de orientación tridimensional para ejes y vectores momento, adoptada universalmente en Solidum FEM para convenciones de signos de magnitudes vectoriales en 3D. Ver capítulo 6.

**Retorno radial (return mapping)**. Algoritmo predictor-corrector para integrar la ecuación constitutiva de plasticidad J2: se asume un paso elástico (predictor), se evalúa el criterio de fluencia y, si la tensión predictora lo viola, se proyecta de vuelta a la superficie de fluencia (corrector). El módulo tangente consistente se obtiene linealizando el algoritmo discreto.

**`SolveResult`**. Agregado inmutable que el `Domain` construye al final de la solución. Contiene los desplazamientos globales, las cargas aplicadas, las reacciones y los esfuerzos internos por elemento. Es la interfaz pública para consumidores externos. Ver ADR 0002.

**Solver no lineal**. Componente que orquesta la solución de un problema no lineal: subdivisión en pasos, iteraciones internas, criterio de convergencia. En Solidum FEM están implementados `LinearSolver`, `NonlinearSolver` y `ArcLengthSolver`. Es nivel estratégico, distinto del subsistema algebraico. Ver capítulo 4.

**`STRAIN_DIM`**. Atributo de clase de cada material y elemento que declara la dimensión del vector de deformación que manejan: 1 (escalar axial), 3 (Voigt-3 bidimensional), 6 (Voigt-6 tridimensional). El sistema lo usa para validar la compatibilidad material↔elemento en el momento de construcción.

**Tensión plana / deformación plana**. Hipótesis de modelado para problemas bidimensionales. La tensión plana asume que la tensión perpendicular al plano es nula (placas delgadas); la deformación plana asume que la deformación perpendicular es nula (cuerpos extensos en una dimensión, secciones de tuberías largas).

**Validación temprana en construcción**. Comprobación que la clase base de elemento ejecuta cuando se le asigna un material: verifica que la dimensión de deformación del material coincida con la del elemento. Errores físicos de combinación se detectan en la construcción del caso, no doscientas iteraciones después.

**Visualization Toolkit (VTK)**. Conjunto de formatos de archivo y bibliotecas para representación de mallas y campos discretos en visualizadores de simulación numérica como ParaView. El `VtkExporter` de Solidum FEM produce los archivos correspondientes.

**YAML (YAML Ain't Markup Language)**. Formato de serialización legible por humanos en el que se describen los casos de Solidum FEM (nodos, materiales, elementos, condiciones de contorno, cargas, solver). Es contrato público del programa.
