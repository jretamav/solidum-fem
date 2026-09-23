# Mecanismos transversales

Los patrones que se describen a continuación estructuran el código de Solidum FEM sin pertenecer a una capa concreta. Aparecen indistintamente en materiales, elementos y solvers; su comprensión es necesaria para la lectura de cualquier zona del programa.

## Registros con auto-registro

Cada categoría de componente dispone de su propio diccionario global: `MaterialRegistry`, `ElementRegistry`, `SolverRegistry`. Cada clase concreta se inscribe en su diccionario mediante un decorador (`@MaterialRegistry.register` y análogos) que se ejecuta en el momento de definición de la clase. El intérprete del archivo YAML no requiere conocimiento *a priori* de los materiales o elementos existentes: consulta el registro por el nombre que aparece en el archivo. La incorporación de un componente nuevo no exige la modificación de listas centrales.

## Descubrimiento automático

El módulo `solidum/autodiscover.py` se invoca una sola vez durante la importación del paquete. Recorre, mediante `pkgutil.iter_modules`, las carpetas canónicas (`solidum/materials`, `solidum/elements`, `solidum/math`) e importa cada módulo. Como los decoradores `@register` se ejecutan en el momento de la importación, esta operación basta para poblar todos los registros sin enumeración manual. Este mecanismo es el equivalente conceptual al `INCLUDE` automático de Fortran moderno, ejecutado a tiempo de ejecución.

## Contratos declarativos mediante atributos de clase

En lugar de imponer métodos abstractos cuando el comportamiento puede inferirse del dato, Solidum FEM declara los contratos como atributos de clase: `STRAIN_DIM` (1, 3 o 6), `DOF_NAMES` (lista de nombres de DOF por nodo), `N_INTEGRATION_POINTS`, `PRIMARY_STATE_VAR`. La clase base los lee y se autoconfigura: registra los DOF, valida la compatibilidad material↔elemento e inicializa el `ElementState` con la forma correcta. Este mecanismo sustituye los métodos `setup()` repetitivos que cada subclase tendría que implementar.

## Validación temprana en construcción

Durante la instanciación de un elemento, su método `__init__` verifica que el material recibido sea dimensionalmente compatible. Un elemento `Truss2D` con material `Elastic2D` produce un error en la construcción del caso, no después de doscientas iteraciones en forma de fallo opaco durante la solución. El coste es una comprobación trivial; el beneficio consiste en que ciertos errores físicos se detectan de forma inmediata en el momento de su introducción.

## Intérprete genérico por introspección de argumentos

El componente `YamlParser` no contiene ramificaciones del tipo `if material_type == "Elastic1D"` por cada material. Inspecciona los argumentos del constructor de la clase recuperada del registro y le transmite los campos del archivo YAML. La incorporación de un material nuevo no requiere modificación del intérprete; es suficiente con que el constructor declare los parámetros con nombres compatibles con el archivo YAML.

## Estados *trial* y comprometido en `ElementState`

Cada elemento mantiene un objeto `ElementState` con dos copias de las variables internas: una copia *trial* (la que se explora durante las iteraciones del solver) y una copia comprometida (la que corresponde al último paso convergido). El solver invoca `commit_state()` solo cuando un paso converge, lo que evita la contaminación del historial plástico o de daño con tentativas posteriormente descartadas. Esta semántica es crítica para problemas con plasticidad o daño: en su ausencia, una iteración no convergida que casualmente quedase dentro de tolerancia contaminaría el historial de forma irreversible.

## Ensamblaje disperso con topología cacheada

El primer ensamblaje de la matriz global calcula los pares de índices (i, j) de cada contribución elemental en formato Coordinate (COO), deriva de ellos la estructura comprimida por filas (CSR) de la matriz y un mapa que lleva cada entrada COO a su posición CSR. En las iteraciones siguientes se reescribe únicamente el vector de datos sobre la misma topología (las matrices elementales de cada familia se escriben directamente en su bloque) y se reduce a CSR con un kernel compilado sobre el mapa inverso —las entradas COO que alimentan cada posición CSR, en orden estable—, serie o paralelo, sin recalcular índices ni reordenar: la conversión COO a CSR que `scipy` rehacía en cada ensamblaje (unos 115 ms para 8 000 hexaedros) desaparece. Se trata del patrón "calcular una vez, reutilizar" típico en programas de MEF compilados, implementado aquí sobre `scipy.sparse` (ADR 0014).

## Ensamblaje por lotes: familias y un único kernel compilado

Los elementos que comparten clase, instancia de material, regla de cuadratura y número de nodos forman una *familia de lote* (`solidum/math/batch/`). La clave se deriva de atributos que ya existen; un elemento nuevo cae en su propia familia sin tocar el ensamblador. Cada familia guarda su estado interno como arreglos (`FamilyState`: variables internas y esfuerzos, committed y trial, una fila por punto de Gauss) y se evalúa completa dentro de un único bucle compilado, `solid_family_kernel`, que recibe la cinemática del elemento (`BATCH_KINEMATICS`) y la constitutiva del material (`BATCH_KERNEL`) como funciones tipadas por su firma. Esa tipificación (`numba.types.FunctionType`) es lo que permite compilar el kernel una sola vez para todas las parejas elemento × material y reutilizarlo desde el caché en disco: un despachador tiene identidad propia por proceso y obligaría a recompilar por pareja en cada ejecución.

El camino por elemento (`compute_element_state`, `compute_state`) sigue siendo el contrato obligatorio y la referencia física. El camino por lotes es una capacidad que el componente declara: los elementos sin cinemática compilada (estructurales 1D, discontinuidad embebida) y los materiales sin esquema de estado se ensamblan como siempre dentro del mismo `Assembler`. Ambos caminos ejecutan las mismas funciones compiladas por punto de Gauss; el barrido de contratos exige que `K` y `F_int` coincidan a precisión de máquina en toda combinación registrada. El kernel recibe siempre coordenadas de referencia y desplazamientos y recomputa jacobiano y `B` en cada evaluación: nada derivado de la configuración deformada se cachea, de modo que el diseño vale para elementos lineales, corotacionales y formulaciones lagrangianas futuras. Medido sobre 10 000 Quad4 con plasticidad J2: 157 → 5,7 µs por elemento y ensamblaje con el kernel serie, y 1,4 µs con el paralelo (ADR 0014).

El kernel de familia existe en dos variantes con la misma aritmética: serie y paralela (`prange` sobre bloques de elementos, con espacio de trabajo por bloque y sin estado compartido entre hilos). Como cada elemento suma en el mismo orden, el resultado es bit a bit el mismo con cualquier número de hilos, cosa que el barrido de contratos exige. Una limitación de Numba condiciona el contrato de la cinemática: una excepción lanzada dentro de una región paralela se pierde en silencio, así que las cinemáticas por lotes no lanzan sino que devuelven un jacobiano no positivo, el kernel marca el punto y el ensamblador convierte la marca en la excepción de siempre, ahora con el identificador del elemento. El mismo protocolo rige la variante serie.

El post-proceso sigue el mismo camino. Cada familia deja en sus elementos una referencia a sí misma al adoptar sus estados, de modo que el exportador VTK y `gauss_states(domain, U)`, que sólo ven el dominio, evalúan por familia —un kernel que devuelve deformación y esfuerzo por punto de Gauss desde el estado committed, sin tocar el trial— y caen al `compute_gauss_state` de cada elemento sólo en los que no pertenecen a ninguna familia. La rigidez inicial (`assemble_system`, `compute_global_stiffness`) es una evaluación auxiliar y restaura el estado trial previo al salir.

## Eliminación directa de condiciones de frontera

Las condiciones de frontera, tanto Dirichlet nodales como restricciones multipunto lineales (MPC), se imponen por eliminación directa (ADR 0004). Toda restricción se expresa en forma afín `u_s = g_s + Σ α_si · u_mi` y se acumula en un `ConstraintSet` (subpaquete `solidum.bc`). El ensamblador construye un operador disperso `T` y un vector `g` tales que `u = T · u_libre + g`: en filas correspondientes a DOF libres, `T` es la matriz identidad rectangular; en filas esclavas, lleva los coeficientes `α_si` en las columnas de los maestros. El sistema entregado al subsistema algebraico es `K_red = TᵀKT`, `F_red = Tᵀ(F − K · g)`. La imposición es exacta a redondeo, preserva la simetría y la positividad definida de `K`, y deja el sistema reducido apto para solvers iterativos (CG, GMRES).

La forma afín cubre con la misma maquinaria los casos de empotramiento, asentamiento prescrito, apoyo en plano oblicuo, periodicidad de celda unitaria y unión rígida master-slave entre nodos. Las restricciones declaradas como cadenas (un esclavo cuyo maestro es a su vez esclavo de otra restricción) se resuelven por cierre transitivo en el momento de construir `T`; los ciclos y las redeclaraciones inconsistentes se detectan en validación temprana. La declaración se hace mediante `Domain.add_linear_constraint` o desde el bloque `linear_constraints` del archivo YAML.

## Compilación Just-In-Time mediante Numba

Las funciones críticas (cinemática de los elementos, núcleos de los algoritmos de retorno, el kernel de familia del ensamblaje por lotes) se decoran con `@njit(cache=True)` y se compilan a código nativo en su primera invocación; el binario se guarda en `__pycache__` y las ejecuciones siguientes lo cargan en centésimas de segundo. Lo que se compila es el procedimiento, no los datos: la configuración deformada, los parámetros y las variables internas entran como argumentos en cada llamada, así que un cambio de geometría o de estado nunca exige recompilar. La restricción consiste en que el código compilado solo admite tipos primitivos y arreglos NumPy, no objetos arbitrarios de Python; por eso los materiales exponen su estado como filas de arreglo (`STATE_SCHEMA`) y sus constantes como un vector de parámetros (`batch_params`).

## Variable principal de visualización

Cada material declara como atributo de clase la variable interna principal a efectos de visualización: por ejemplo `'damage'` en un modelo de daño o `'alpha'` en plasticidad con endurecimiento. El componente `VtkExporter` la lee de forma genérica sin conocimiento del material que la origina. Este mecanismo permite la incorporación de materiales nuevos con visualización automática, sin modificación del exportador.

## Cuadraturas centralizadas

Las tablas y reglas de cuadratura de Gauss-Legendre 1D, 2D y 3D residen centralizadas en `solidum/math/integration.py`. Cada elemento declara su `N_INTEGRATION_POINTS` y consume los puntos y pesos correspondientes, sin reproducción de tablas en cada subclase.

## Política unificada de tolerancias: patrón `atol + rtol · escala`

Toda comparación con significado físico o iterativo del proyecto — admisibilidad constitutiva, convergencia de Newton-Raphson, convergencia de arc-length, futuros criterios de cierre de gap en contacto o de inversión del jacobiano — sigue una misma fórmula estructural:

```
magnitud  ≤  atol  +  rtol · escala(estado)
```

donde `atol` es un piso absoluto en las unidades físicas del problema, `rtol` es una banda relativa adimensional y `escala(estado)` es la magnitud característica del criterio en el estado corriente. Esta forma única garantiza tres propiedades simultáneamente: invariancia bajo cambio de unidades (la `escala` se ajusta), adaptatividad al estado (la tolerancia crece con la evolución del problema, p. ej. con el endurecimiento plástico) y robustez en regímenes degenerados (cuando la escala colapsa transitoriamente, `atol` mantiene la comparación significativa).

El patrón vive centralizado por subsistema, no replicado: la admisibilidad constitutiva se aplica en `Material.is_admissible` y `Material.admissibility_tol` (ADR 0006); la convergencia de los solvers no lineales reside en la clase `ConvergenceCriterion` de `solidum/math/convergence.py` (ADR 0007). Cada material o solver declara su `escala` mediante un método ligero (por ejemplo `Material.admissibility_scale`, que devuelve la fluencia corriente o el umbral de daño), y la fórmula no se reproduce a mano en ningún sitio salvo cuando una restricción técnica (kernel compilado con Numba) obliga a precomputar la tolerancia fuera del kernel — y, aun así, vía el método centralizado.

Adicionalmente, los términos absolutos `atol` se autoderivan de las escalas del problema en su primer ensamblaje, no se codifican como constantes globales con unidades. Las constantes globales del proyecto que rigen esta política son adimensionales (`CONVERGENCE_RTOL_FORCE`, `CONVERGENCE_ATOL_FORCE_FACTOR`, `ADMISSIBILITY_TOL_REL`, etc.), lo que mantiene el código independiente del sistema de unidades elegido por el usuario.

La importancia de este mecanismo es estructural: la convergencia es lo que marca el éxito de un análisis numérico, y la diferencia entre un código robusto y uno frágil está en la disciplina con que se construyen estas comparaciones. Toda extensión futura que introduzca un nuevo criterio de comparación contra cero adopta este patrón.

## Modos de cuerpo rígido derivados del nombre de los grados de libertad

Los modos de cuerpo rígido de un modelo —traslaciones, giros infinitesimales y, para un campo escalar como la temperatura, el vector constante— se derivan **sólo del nombre de los grados de libertad** (`ux`, `uy`, `uz`, `rx`, `ry`, `rz`, `T`) y de las coordenadas nodales, sin información del tipo de elemento (`solidum/math/linalg/nullspace.py`). Por eso cubren cualquier combinación del catálogo, incluidos los dominios mixtos, y un elemento nuevo los hereda sin hacer nada. Son el núcleo exacto de la matriz de rigidez de un modelo sin apoyos, verificado en todas las familias.

Los consumen dos piezas: el precondicionador multimalla del solver iterativo, que sin ellos converge peor que sin precondicionar (ADR 0018), y la detección de mecanismos del análisis estático (ADR 0019). El ensamblador los calcula una vez y los cachea mientras no cambie la topología.

## Red de seguridad del análisis estático

El usuario de Solidum no elige el solver algebraico ni tiene por qué conocer cómo falla. Un solver directo ante una matriz singular devuelve resultados absurdos sin avisar (medido: desplazamientos de miles de kilómetros). El análisis estático interpone tres capas (ADR 0019, `solidum/math/solvers/model_checks.py`):

1. **Antes de resolver**: los modos de cuerpo rígido que las restricciones no impiden son un **mecanismo**. Se calculan como el núcleo de la matriz que evalúa cada modo sobre las restricciones —incluidas las lineales, que entran por el operador de reducción del ADR 0004— y se describen en términos del modelo: *"traslación en la dirección x"*, *"giro alrededor del eje z que pasa por (0, 0.5)"*, *"el campo 'T' no tiene ningún valor prescrito"*. No se aplica a análisis modales o dinámicos, donde un modelo libre es legítimo.
2. **Después de resolver** (estático lineal): la solución debe satisfacer el equilibrio, `‖F − K·u‖ ≤ 10⁻⁸·‖F‖`.
3. **Después de resolver** (estático lineal): la factorización no debe tener pivotes numéricamente nulos, con el mismo criterio en SuperLU y en Pardiso. Delata un mecanismo interno aunque ninguna carga lo active.

Dentro de un Newton no se rechaza nada: cerca de un punto límite resolver un sistema casi singular es parte del algoritmo, y la garantía la da el propio Newton, que exige equilibrio real. Allí el residuo del sistema lineal sólo se usa para diagnosticar mejor la divergencia.

El principio que el mecanismo materializa es general: **un problema mal planteado lo detecta el sistema y lo explica en el lenguaje del modelo**, no lo diagnostica el usuario a partir de un síntoma numérico.

