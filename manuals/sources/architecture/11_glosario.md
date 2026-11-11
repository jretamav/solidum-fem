# Glosario

Términos técnicos centrales de Solidum FEM, ordenados alfabéticamente. Cada entrada da una definición breve (de una a tres frases) y, cuando procede, remite al capítulo donde el concepto se desarrolla.

**Adherencia-deslizamiento**. Transmisión de fuerza entre una barra de refuerzo y el sólido que la rodea a través de una ley entre el esfuerzo de adherencia y el deslizamiento relativo de la barra, en lugar de la unión perfecta. En Solidum la genera el refuerzo embebido con `bond`: una interfaz `BondSlipInterface` por tramo de barra, con una ley de interfaz de línea (`BondSlipModelCode`, `BondSlipMultilinear`). ADR 0032 y 0033.

**Ajuste de paso (incremento de carga)**. División de la carga total en fracciones que se aplican sucesivamente en problemas no lineales. Permite que el método de Newton-Raphson converja a partir de una aproximación cercana a la solución del paso anterior.

**Corrotacional (formulación)**. Descomposición del movimiento de un elemento en una rotación de cuerpo rígido más una deformación local pequeña. La ley constitutiva trabaja en la configuración local; la cinemática del elemento absorbe las grandes rotaciones. Aplicable a armaduras, cables y marcos en régimen geométricamente no lineal con deformación pequeña.

**Architecture Decision Record (ADR)**. Documento breve que registra una decisión de arquitectura junto con su contexto, alternativas consideradas y consecuencias. Reside en `docs/adr/000N-titulo.md`. Ver capítulo 8.

**`ArcProblem`**. Problema de Newton por paso del arc-length (iterado `(U, λ, ΔU)`, dos resoluciones por iteración y la restricción del paso como método `constraint`), exportado con nombre público en `solidum.math.solvers`. Una variante de arco que sólo cambia la restricción lo hereda y redefine `constraint`: así lo hacen `IndirectDisplacementSolver` en el programa principal y `DissipationArcLengthSolver` en el módulo de usuario `dissipation_control`. Ver capítulo 5, ADR 0015 (protocolo `NewtonProblem`) y ADR 0020.

**Autodiscover (descubrimiento automático)**. Procedimiento que recorre las carpetas canónicas del proyecto e importa todos sus módulos al inicio del programa, lo que dispara la ejecución de los decoradores `@register` y, con ello, la inscripción de cada componente en su registro. No recorre `solidum/user/`: los módulos de usuario sólo se cargan cuando un caso los pide. Ver capítulo 5.

**Capa algebraica**. Subsistema en `solidum/math/linalg/` responsable de resolver el sistema lineal `K · x = b` que aparece dentro de cada iteración del solver no lineal. Está separada del solver no lineal y selecciona automáticamente el algoritmo de factorización adecuado mediante un despachador. Ver capítulo 4 y ADR 0003.

**Corrector de Newton**. Bucle de corrección compartido por los solvers iterativos (`NewtonCorrector`, ADR 0015): ensambla una vez en el iterado, evalúa el criterio dual con el residuo de ese ensamblaje y el incremento previo, resuelve con el backend algebraico (degradación a LU, factorización congelada opcional), aplica el line search opcional y comitea al converger. Cada solver le entrega un `NewtonProblem` con su residuo, su sistema tangente y su actualización, y conserva su propio control de paso. Ver capítulos 3 y 4.

**Catálogo**. Documento navegable que enumera todos los componentes implementados de una categoría (elementos, materiales, solvers, restricciones generadas) con sus parámetros, validez y referencias. Reside en `docs/catalogo_<categoría>.md`. Es índice y resumen; el detalle vive en las especificaciones.

**Contrato declarativo**. Mecanismo por el cual una clase declara, mediante atributos de clase como `STRAIN_DIM`, `DOF_NAMES`, `N_INTEGRATION_POINTS`, `PRIMARY_STATE_VAR` o `IS_LINEAR_ELASTIC`, las propiedades fijas que el sistema necesita conocer de ella. La clase base lee estos atributos y se autoconfigura, sustituyendo métodos de inicialización repetitivos en cada subclase. Ver capítulo 5.

**Despachador algebraico**. Función `select_solver(props, override)` que recibe las propiedades de la matriz de rigidez (simetría, definición positiva) y devuelve el algoritmo adecuado: Cholesky cuando la matriz es simétrica definida positiva y CHOLMOD está instalado, Pardiso si está instalado (ADR 0017) y LU en otro caso; el iterativo sólo se usa a petición (ADR 0018). Ver ADR 0003.

**`DOF_NAMES`**. Atributo de clase de cada elemento que enumera los nombres de los grados de libertad activos en cada uno de sus nodos. La clase base lo lee para registrar los grados de libertad globales en el `Domain`.

**`Domain`**. Objeto central de la capa de dominio que mantiene la lista de nodos, la conectividad de los elementos, la numeración global de los grados de libertad y, tras la solución, el agregado `SolveResult`.

**Etapa (de análisis)**. Análisis que parte del estado que dejó otro: los desplazamientos, lo que guardan los elementos y la carga que ya actúa, que pasa como carga constante. Así se aplican la gravedad y después el empuje lateral, un transitorio tras un estático o un modal sobre una estructura cargada. Desde Python, `solidum.run(..., previous=resultado)`; desde el YAML, el bloque `stages`. ADR 0026, 0029 y 0031.

**Eliminación directa (de condiciones de Dirichlet)**. Estrategia de imposición de condiciones esenciales que retira los grados de libertad prescritos del sistema antes de resolverlo, en oposición a la penalización (que añade un valor grande a la diagonal) y a los multiplicadores de Lagrange (que amplían el sistema). Toda restricción se expresa en forma afín `u_s = g_s + Σ α_si · u_mi` y se acumula en un `ConstraintSet` (`solidum/bc/`). El ensamblador construye un operador disperso `T` y un vector `g` tales que `u = T · u_libre + g`, y entrega al solver el sistema reducido `K_red = TᵀKT`, `F_red = Tᵀ(F − K · g)`. La imposición es exacta a redondeo y preserva la simetría y la positividad definida de la matriz original. Ver ADR 0004.

**Elemento**. Entidad geométrica que construye su matriz de gradientes `B`, su matriz de rigidez tangente y su vector de fuerzas internas a partir de los desplazamientos nodales y del material asignado (o, en un marco, de su sección). Sus contratos declarativos están en `solidum/elements/`. `Element.materials()` enumera todos los materiales que usa —el suyo, o los de las fibras de su sección— para las comprobaciones agregadas de masa, peso propio y simetría (ADR 0021). Ver capítulo 4.

**`ElementState`**. Objeto que encapsula las dos copias de las variables internas de cada elemento: una copia *trial* (la que se explora durante las iteraciones del solver) y otra comprometida (la del último paso convergido). Ver capítulo 5.

**Ensamblaje con topología cacheada**. Estrategia de ensamblaje en la que la primera invocación calcula los pares de índices `(i, j)` de cada contribución elemental en formato Coordinate (COO), deriva la estructura CSR de la matriz y un mapa COO → CSR; las invocaciones siguientes reescriben únicamente el vector de datos y lo acumulan sobre la misma estructura, sin recalcular índices ni reordenar. Ver capítulo 5 y ADR 0014.

**Ensamblaje por lotes**. Evaluación de todos los puntos de Gauss de una familia de lote dentro de un único kernel compilado (`solid_family_kernel`, con variante paralela `solid_family_kernel_parallel` bit a bit equivalente), que recibe la cinemática del elemento y la constitutiva del material como funciones tipadas por su firma. El post-proceso por familia (`Family.gauss_state`) usa las mismas funciones. Coexiste con el camino por elemento, que sigue siendo el contrato obligatorio y la referencia física. Ver capítulo 5 y ADR 0014.

**Esfuerzo plano / deformación plana**. Hipótesis de modelado para problemas bidimensionales. El esfuerzo plano asume que el esfuerzo perpendicular al plano es nulo (placas delgadas); la deformación plana asume que la deformación perpendicular es nula (cuerpos extensos en una dimensión, secciones de tuberías largas).

**Especificación (spec)**. Documento que describe un componente físico nuevo antes de su implementación: especificación física, formulación numérica, contrato YAML y criterios de aceptación. Reside en `docs/specs/<Nombre>.md`. Es a la vez orden de trabajo y referencia detallada del componente. Ver capítulo 7.

**Estado *trial* y estado comprometido**. Las dos copias de las variables internas de cada elemento. El estado *trial* se explora durante las iteraciones del solver no lineal; al converger un paso, el solver invoca `commit_state()` y promueve el *trial* a estado comprometido. Esta semántica es indispensable en problemas con plasticidad o daño para no contaminar el historial con tentativas posteriormente descartadas.

**Familia de elementos**. Conjunto de elementos que comparten la cinemática esencial: armaduras, cables, marcos (vigas), sólidos 2D, sólidos 3D, cáscaras. Ver capítulo 4.

**Familia de lote**. Conjunto de elementos de un dominio que comparten clase de elemento, instancia de material, regla de cuadratura y número de nodos, y que el ensamblador evalúa juntos en un kernel compilado. Su estado interno vive en un `FamilyState` (arreglos por punto de Gauss) del que `elem.state` es una vista (`BatchedElementState`). Ver capítulo 5 y ADR 0014.

**Familia de material**. Registro que declara la sección de primer nivel del YAML con sus objetos (`YAML_SECTION`): `materials` para los mecánicos, `thermal_materials` para los térmicos e `interface_materials` para las leyes tracción–desplazamiento relativo de una interfaz (ADR 0033). `Registry.families()` las enumera, en su orden de definición, junto con las familias que no son de materiales (hoy, la de secciones de marco, `sections`), y el intérprete del YAML las recorre sin conocerlas por nombre; un elemento declara a qué familia pertenecen sus materiales o su sección con `REFERENCE_KWARGS`. Ver capítulo 5, ADR 0020 y ADR 0021.

**`from_yaml` (gancho)**. Método de clase opcional `from_yaml(data, resolve)` de una clase de familia con referencias anidadas a otras familias (una sección por fibras referencia materiales dentro de sus parches). El intérprete del YAML lo llama en lugar de `cls(**datos)` y le pasa `resolve(registro, id, contexto)`, que devuelve el objeto de una familia ya construida o lanza un error con el contexto de la referencia. Las familias se construyen en orden de definición. Ver capítulo 5 y ADR 0021.

**Gradiente conjugado (CG) y MINRES**. Métodos iterativos de Krylov para sistemas lineales simétricos: CG para matrices definidas positivas, MINRES para indefinidas. El solver iterativo de Solidum pasa de uno a otro de forma automática al detectar curvatura negativa. Capítulo 4, ADR 0018.

**Generador de restricciones (restricción generada)**. Objeto que, sobre un dominio ya mallado, crea los nodos, los elementos y las restricciones multipunto que pide una geometría, antes de numerar las ecuaciones (`ConstraintRegistry`, bloque YAML `constraints`, `solidum/bc/`). El primero es el refuerzo embebido. ADR 0032.

**Grado de libertad (DOF)**. Cada incógnita escalar del problema discreto. En estática mecánica, los grados de libertad son los desplazamientos nodales y, en su caso, las rotaciones nodales. La numeración global se establece en el `Domain`.

**Grado de libertad sin masa**. Grado de libertad en el que la matriz de masa es nula (masas sólo en algunos nudos, elementos sin densidad, giros sin inercia rotacional). El análisis modal los condensa estáticamente de forma exacta y los transitorios implícitos los admiten con una aceleración inicial consistente. ADR 0024.

**IMPL-EX (integración implícita-explícita)**. Forma de integrar un modelo de daño (Oliver, Huespe y Cante 2008): la variable histórica se extrapola de los dos últimos pasos confirmados para calcular el esfuerzo del paso, de modo que el material es lineal dentro del paso (tangente secante, simétrica y definida positiva), y se actualiza de forma implícita al confirmar el paso. Robusta en el ablandamiento a cambio de un error proporcional al paso. ADR 0034.

**`IS_LINEAR_ELASTIC`**. Atributo de clase de cada material (`False` por defecto) que declara que su respuesta es lineal y sin historia: `True` en `Elastic1D`, `Elastic2D`, `Elastic3D` y `Orthotropic2D`. Los marcos lo exigen en su camino de forma cerrada (`material` + propiedades geométricas) y rechazan al construir cualquier otro material; la no linealidad material de un marco entra por su sección. Ver capítulo 5 y ADR 0021.

**`LoadVector`**. Vector de cargas que lleva consigo las cargas equivalentes por elemento de su parte distribuida y las combina linealmente con la carga (`F_a ± F_b`, `a·F`). Así las fuerzas internas de los elementos cargados son correctas con varias cargas, casos o etapas. ADR 0027.

**Localización de puntos**. Búsqueda del elemento sólido que contiene un punto y de sus coordenadas naturales, por inversión del mapeo isoparamétrico con Newton a partir de lo que el elemento declara (`REFERENCE_CELL`, funciones de forma), sin código por tipo de elemento (`solidum/mesh/point_location.py`). La usan los generadores de restricciones. ADR 0032.

**Longitud característica (del punto de integración)**. Longitud `h` asociada a cada punto de integración con la que un material con ablandamiento regularizado por energía de fractura ajusta su rama descendente para disipar `G_f/h` por unidad de volumen. La siembra el elemento al construirse y viaja en el estado del punto; el material declara si la usa (`uses_characteristic_length`). ADR 0028.

**Longitud de arco (método de Crisfield)**. Estrategia de solución no lineal que introduce el factor de carga como incógnita adicional y añade una restricción geométrica sobre la trayectoria en el espacio (U, λ). Permite recorrer ramas con derivada infinita o negativa, imprescindible para problemas con *snap-back* o *snap-through*.

**Masa de nudo**. Masa concentrada declarada como dato del dominio en grados de libertad de nudos, traslacionales o rotacionales (`Domain.add_nodal_mass`, bloques YAML `nodal_masses*`). Es sólo inercia: no genera peso con `gravity`. ADR 0024.

**Material**. Ley constitutiva que, dada una deformación y un estado interno, devuelve el esfuerzo y el módulo tangente consistente. Su contrato declarativo está en `solidum/materials/`. Ver capítulo 4.

**Material de interfaz**. Ley entre la tracción sobre una interfaz y el desplazamiento relativo entre sus dos lados, en el sistema local de la interfaz (normales primero): de superficie para grietas y juntas, de línea para barras embebidas. Familia propia (`InterfaceMaterialRegistry`, sección YAML `interface_materials`, base `InterfaceMaterial`). ADR 0033.

**Mecanismo**. Modelo cuyos apoyos o cuya conectividad permiten un movimiento sin deformación: la matriz de rigidez es singular y el problema estático no tiene solución única. Si el movimiento es de sólido rígido del conjunto, Solidum lo detecta antes de resolver y lo describe; si es interno (una rótula de más), lo detecta tras resolver. Capítulo 5, ADR 0019.

**Método de Elementos Finitos (MEF)**. Técnica de discretización para resolver ecuaciones diferenciales parciales mediante la subdivisión del dominio en elementos sobre los que se interpolan los campos incógnita. Solidum FEM trabaja en aproximación de desplazamientos.

**Modos de cuerpo rígido**. Movimientos de un modelo sin apoyos que no producen deformación: traslaciones y giros infinitesimales, y el valor constante en un campo escalar. En Solidum se derivan del nombre de los grados de libertad y los usan el precondicionador multimalla y la detección de mecanismos. Capítulo 5.

**Módulo de usuario**. Formulación no estándar que vive fuera del programa principal, en su propia carpeta `solidum/user/<nombre>/`, a la manera de los elementos y materiales de usuario de FEAP. El programa principal no la importa ni la nombra: un caso la carga con `user_modules: [<nombre>]` en el YAML o con `solidum.load_user_module` en Python. El primero es `discontinuities` (discontinuidad interior embebida); el segundo, `dissipation_control` (arco por disipación de energía, un solver). Ver capítulo 5 y ADR 0020.

**Multimalla algebraico (AMG)**. Precondicionador que resuelve el sistema en una jerarquía de problemas cada vez más pequeños construidos a partir de la propia matriz. Con los modos de cuerpo rígido como casi-núcleo, el número de iteraciones apenas depende del tamaño del modelo. ADR 0018.

**Multipoint constraint (MPC)**. Restricción afín lineal que liga el desplazamiento de un grado de libertad esclavo a una combinación lineal de los desplazamientos de uno o varios grados de libertad maestros, posiblemente con un término independiente. Modelan apoyos en plano oblicuo, periodicidad de celda unitaria, uniones rígidas entre nodos y simetrías no alineadas con los ejes globales. Se declaran mediante `Domain.add_linear_constraint` o desde el bloque `linear_constraints` del archivo YAML, y se imponen por eliminación directa con la misma maquinaria que las condiciones de Dirichlet. Ver capítulo 5 y ADR 0004.

**Notación de Voigt**. Representación de un tensor simétrico de segundo orden como un vector. En 2D, los tres componentes son las dos componentes normales y la componente de cortante (Voigt-3); en 3D, los seis componentes son las tres normales y las tres de cortante (Voigt-6).

**Numba (compilación Just-In-Time, JIT)**. Decorador `@njit` que compila a código nativo, en la primera invocación, las funciones a las que se aplica. La primera ejecución asume el coste de compilación; las siguientes corren a velocidad cercana a Fortran. Restringido a tipos primitivos y arreglos NumPy.

**Pivote pequeño**. Valor de la diagonal de la factorización pequeño frente a la diagonal de la matriz (`|U_ii|/√(|K_rr|·|K_cc|) < 10⁻⁸`): puede delatar una matriz singular, pero su tamaño también depende del contraste de rigideces y del tamaño del modelo. Desde el 20 de octubre de 2026 Solidum lo usa sólo como filtro: el veredicto sobre un mecanismo interno lo da la energía por elemento del modo de menor energía. ADR 0019.

**Perturbación lineal**. Análisis lineal estático (`LinearSolver`) sobre el estado de una etapa, con la rigidez tangente de los elementos y sin cambiar ese estado: da casos superponibles sobre la estructura cargada (`SolveResult.perturbation`). ADR 0031.

**Pseudo-tiempo**. Parámetro de control de un solver que avanza por pasos (factor de carga, longitud de arco, magnitud del control indirecto, tiempo). El solver declara su incremento en cada intento de paso y el ensamblador entrega a los materiales IMPL-EX la razón con el incremento del último paso confirmado. ADR 0034.

**`PRIMARY_STATE_VAR`**. Atributo de clase de cada material que declara cuál de sus variables internas es la principal a efectos de visualización. El exportador VTK la lee genéricamente sin conocimiento del material que la origina.

**Programa principal**. Solidum estándar: los elementos, materiales y solvers de elementos finitos clásicos, es decir, todo `solidum/` salvo `solidum/user/`. Admite un componente según el criterio de admisión de Reglas.md §1 (textos de referencia y código general establecido). No conoce los módulos de usuario; el test de pureza (`tests/test_main_program_purity.py`) lo vigila con tres controles: estructural (ningún import de un módulo de usuario), léxico (ninguna mención de lo que un módulo registra, con los nombres derivados automáticamente) y en ejecución (`import solidum` no carga ningún módulo). Ver capítulo 5 y ADR 0020.

**Registro (registry)**. Diccionario global por categoría (`MaterialRegistry`, `ThermalMaterialRegistry`, `InterfaceMaterialRegistry`, `SectionRegistry`, `ElementRegistry`, `SolverRegistry`, `ConstraintRegistry`) que mapea el nombre de cada componente con su clase. Todos son subclases de la clase genérica `Registry` y declaran sus metadatos con atributos públicos de clase: `YAML_LABEL`, la etiqueta de un objeto del registro en los mensajes de error, que tiene todo registro (`"elemento"`, `"solver"`, `"material térmico"`, `"material de interfaz"`, `"restricción"`) y sustituye al antiguo atributo privado `_kind`, ya retirado; `SPEC_KIND`, el tipo de spec de sus componentes; y, sólo en las familias (de material o de secciones), `YAML_SECTION`. Cada clase se inscribe automáticamente en su registro mediante un decorador. Un módulo de usuario sólo puede añadir entradas: reasignar una existente hace fallar y deshacer su carga. Ver capítulo 5.

**Refuerzo embebido**. Barra de refuerzo con su propia malla de elementos de línea que atraviesa una malla de sólidos sin coincidir con sus nodos y queda atada al campo interpolado del sólido (adherencia perfecta) o unida a él por interfaces de adherencia-deslizamiento (`EmbeddedReinforcement`). ADR 0032.

**Regla de la mano derecha (RHR)**. Convención de orientación tridimensional para ejes y vectores momento, adoptada universalmente en Solidum FEM para convenciones de signos de magnitudes vectoriales en 3D. Ver capítulo 6.

**Relleno (fill-in)**. Entradas no nulas que aparecen en los factores de una factorización directa donde la matriz tenía ceros. Crece con el tamaño del modelo y domina el tiempo y la memoria de un solver directo: de ×10 a ×42 la matriz original entre 4 000 y 53 000 grados de libertad. ADR 0017.

**Retorno radial (return mapping)**. Algoritmo predictor-corrector para integrar la ecuación constitutiva de plasticidad J2: se asume un paso elástico (predictor), se evalúa el criterio de fluencia y, si el esfuerzo predictor lo viola, se proyecta de vuelta a la superficie de fluencia (corrector). El módulo tangente consistente se obtiene linealizando el algoritmo discreto.

**Seguimiento de grieta**. Procedimiento del módulo de usuario `discontinuities` que decide, entre pasos, por dónde continúa una grieta embebida de un elemento al siguiente, para que las particiones de elementos vecinos sean coherentes: global por región conexa (`GlobalCrackTracking`) o con una trayectoria prescrita (`PrescribedCrackPath`). ADR 0010, fase E.

**Sección de marco**. Ley constitutiva de un marco a nivel de sección transversal (ADR 0021): recibe las deformaciones generalizadas del eje en un punto de Gauss, `[ε₀, κ_z]` o `[ε₀, κ_y, κ_z]`, y devuelve las fuerzas internas `[N, M_z]` o `[N, M_y, M_z]` y su tangente. Familia propia (`SectionRegistry`, sección YAML `sections`, `kind: section`), con la base abstracta `Section` en `solidum/core/section.py`. No guarda estado: la historia de sus puntos materiales en cada punto de Gauss vive en el elemento. La consumen los marcos Euler-Bernoulli con su parámetro `section`. Ver capítulo 4.

**Sección por fibras (`FiberSection`)**. Sección de marco discretizada en fibras `(y, z, A, material)`, referidas al eje del elemento, cada una con su material 1D y su historia; las fuerzas internas son las sumas de `σ_i A_i` ponderadas por `1`, `z_i` y `−y_i`. Da plasticidad y daño distribuidos en la sección. Las fibras se generan con parches rectangulares y circulares, capas de barras y fibras sueltas. Ver capítulo 4, ADR 0021 y la spec `FiberSection`.

**Solver directo / solver iterativo**. El directo factoriza la matriz y resuelve por sustitución: robusto ante cualquier condicionamiento, pero su memoria crece más deprisa que el modelo. El iterativo aproxima la solución hasta una tolerancia sin factorizar: memoria proporcional al modelo, pero sensible al condicionamiento. Capítulo 4, ADR 0017 y 0018.

**`SolveResult`**. Agregado inmutable que el `Domain` construye al final de la solución. Contiene los desplazamientos globales, las cargas aplicadas, las reacciones y las fuerzas internas (N, V, M, T) por elemento. Es la interfaz pública para consumidores externos. Ver ADR 0002.

**Solver no lineal**. Componente que orquesta la solución de un problema no lineal: subdivisión en pasos, iteraciones internas, criterio de convergencia. En el programa principal de Solidum FEM, los estáticos no lineales son `NonlinearSolver`, `ArcLengthSolver` e `IndirectDisplacementSolver`, y los transitorios no lineales `NewtonNewmarkSolver` y `NewtonHHTSolver`; el módulo de usuario `dissipation_control` añade `DissipationArcLengthSolver`. Es nivel estratégico, distinto del subsistema algebraico. Ver capítulo 4.

**`STRAIN_DIM`**. Atributo de clase de cada material y elemento que declara la dimensión del vector de deformación que manejan: 1 (escalar axial), 3 (Voigt-3 bidimensional), 6 (Voigt-6 tridimensional). El sistema lo usa para validar la compatibilidad material↔elemento en el momento de construcción.

**Terminación**. Cómo terminó un análisis, legible en atributos y no sólo en el registro: `termination` (`max_lambda`, `control_target`, `max_steps`, `t_end` o `diverged`), `completed`, el número de pasos y el último estado convergido, también tras una excepción. ADR 0025.

**Validación temprana en construcción**. Comprobación que la clase base de elemento ejecuta cuando se le asigna un material: verifica que la dimensión de deformación del material coincida con la del elemento. Errores físicos de combinación se detectan en la construcción del caso, no doscientas iteraciones después.

**Visualization Toolkit (VTK)**. Conjunto de formatos de archivo y bibliotecas para representación de mallas y campos discretos en visualizadores de simulación numérica como ParaView. El `VtkExporter` de Solidum FEM produce los archivos correspondientes.

**YAML (YAML Ain't Markup Language)**. Formato de serialización legible por humanos en el que se describen los casos de Solidum FEM (nodos, materiales, elementos, condiciones de contorno, cargas, solver). Es contrato público del programa.
