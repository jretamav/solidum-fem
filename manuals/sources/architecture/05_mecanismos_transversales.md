# Mecanismos transversales

Los patrones que se describen a continuación estructuran el código de Fenix FEM sin pertenecer a una capa concreta. Aparecen indistintamente en materiales, elementos y solvers; su comprensión es necesaria para la lectura de cualquier zona del programa.

## Registros con auto-registro

Cada categoría de componente dispone de su propio diccionario global: `MaterialRegistry`, `ElementRegistry`, `SolverRegistry`. Cada clase concreta se inscribe en su diccionario mediante un decorador (`@MaterialRegistry.register` y análogos) que se ejecuta en el momento de definición de la clase. El intérprete del archivo YAML no requiere conocimiento *a priori* de los materiales o elementos existentes: consulta el registro por el nombre que aparece en el archivo. La incorporación de un componente nuevo no exige la modificación de listas centrales.

## Descubrimiento automático

El módulo `fenix/autodiscover.py` se invoca una sola vez durante la importación del paquete. Recorre, mediante `pkgutil.iter_modules`, las carpetas canónicas (`fenix/materials`, `fenix/elements`, `fenix/math`) e importa cada módulo. Como los decoradores `@register` se ejecutan en el momento de la importación, esta operación basta para poblar todos los registros sin enumeración manual. Este mecanismo es el equivalente conceptual al `INCLUDE` automático de Fortran moderno, ejecutado a tiempo de ejecución.

## Contratos declarativos mediante atributos de clase

En lugar de imponer métodos abstractos cuando el comportamiento puede inferirse del dato, Fenix FEM declara los contratos como atributos de clase: `STRAIN_DIM` (1, 3 o 6), `DOF_NAMES` (lista de nombres de DOF por nodo), `N_INTEGRATION_POINTS`, `PRIMARY_STATE_VAR`. La clase base los lee y se autoconfigura: registra los DOF, valida la compatibilidad material↔elemento e inicializa el `ElementState` con la forma correcta. Este mecanismo sustituye los métodos `setup()` repetitivos que cada subclase tendría que implementar.

## Validación temprana en construcción

Durante la instanciación de un elemento, su método `__init__` verifica que el material recibido sea dimensionalmente compatible. Un elemento `Truss2D` con material `Elastic2D` produce un error en la construcción del caso, no después de doscientas iteraciones en forma de fallo opaco durante la solución. El coste es una comprobación trivial; el beneficio consiste en que ciertos errores físicos se detectan de forma inmediata en el momento de su introducción.

## Intérprete genérico por introspección de argumentos

El componente `YamlParser` no contiene ramificaciones del tipo `if material_type == "Elastic1D"` por cada material. Inspecciona los argumentos del constructor de la clase recuperada del registro y le transmite los campos del archivo YAML. La incorporación de un material nuevo no requiere modificación del intérprete; es suficiente con que el constructor declare los parámetros con nombres compatibles con el archivo YAML.

## Estados *trial* y comprometido en `ElementState`

Cada elemento mantiene un objeto `ElementState` con dos copias de las variables internas: una copia *trial* (la que se explora durante las iteraciones del solver) y una copia comprometida (la que corresponde al último paso convergido). El solver invoca `commit_state()` solo cuando un paso converge, lo que evita la contaminación del historial plástico o de daño con tentativas posteriormente descartadas. Esta semántica es crítica para problemas con plasticidad o daño: en su ausencia, una iteración no convergida que casualmente quedase dentro de tolerancia contaminaría el historial de forma irreversible.

## Ensamblaje disperso con caché en formato Coordinate

El primer ensamblaje de la matriz global calcula los pares de índices (i, j) de cada contribución elemental y los almacena en formato Coordinate (COO). En las iteraciones siguientes, se reescribe únicamente el vector de datos sobre la misma topología, sin recalculo de índices. Se trata del patrón "calcular una vez, reutilizar" típico en programas de MEF compilados, implementado aquí sobre `scipy.sparse`.

## Eliminación directa de condiciones de Dirichlet

Las condiciones de contorno de tipo Dirichlet se imponen por eliminación directa (ADR 0004). Toda restricción se expresa en forma afín `u_s = g_s + Σ α_si · u_mi` y se acumula en un `ConstraintSet` (subpaquete `fenix.bc`). El ensamblador construye un operador disperso `T` que selecciona los DOF libres y un vector `g` con los términos prescritos, de modo que `u = T · u_libre + g`. El sistema entregado al subsistema algebraico es `K_red = TᵀKT`, `F_red = Tᵀ(F − K · g)`. La imposición es exacta a redondeo, preserva la simetría y la positividad definida de `K`, y deja el sistema reducido apto para solvers iterativos (CG, GMRES). La forma afín cubre además, sin cambios en el solver, los casos de asentamiento prescrito y, en una fase posterior, las restricciones multipunto lineales (uniones rígidas, periodicidad, simetrías oblicuas).

## Compilación Just-In-Time mediante Numba

Las funciones críticas (ensamblaje elemento→global, núcleo del algoritmo de retorno radial) se decoran con `@njit` y se compilan a código nativo en su primera invocación. La primera ejecución asume el coste de compilación; las invocaciones siguientes se ejecutan a velocidad cercana a la de un programa Fortran compilado. La restricción consiste en que el código compilado solo admite tipos primitivos y arreglos NumPy, no objetos arbitrarios de Python.

## Variable principal de visualización

Cada material declara como atributo de clase la variable interna principal a efectos de visualización: por ejemplo `'damage'` en un modelo de daño o `'alpha'` en plasticidad con endurecimiento. El componente `VtkExporter` la lee de forma genérica sin conocimiento del material que la origina. Este mecanismo permite la incorporación de materiales nuevos con visualización automática, sin modificación del exportador.

## Cuadraturas centralizadas

Las tablas y reglas de cuadratura de Gauss-Legendre 1D, 2D y 3D residen centralizadas en `fenix/math/integration.py`. Cada elemento declara su `N_INTEGRATION_POINTS` y consume los puntos y pesos correspondientes, sin reproducción de tablas en cada subclase.
