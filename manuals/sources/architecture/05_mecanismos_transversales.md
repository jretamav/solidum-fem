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

## Ensamblaje disperso con caché en formato Coordinate

El primer ensamblaje de la matriz global calcula los pares de índices (i, j) de cada contribución elemental y los almacena en formato Coordinate (COO). En las iteraciones siguientes, se reescribe únicamente el vector de datos sobre la misma topología, sin recalculo de índices. Se trata del patrón "calcular una vez, reutilizar" típico en programas de MEF compilados, implementado aquí sobre `scipy.sparse`.

## Eliminación directa de condiciones de frontera

Las condiciones de frontera, tanto Dirichlet nodales como restricciones multipunto lineales (MPC), se imponen por eliminación directa (ADR 0004). Toda restricción se expresa en forma afín `u_s = g_s + Σ α_si · u_mi` y se acumula en un `ConstraintSet` (subpaquete `solidum.bc`). El ensamblador construye un operador disperso `T` y un vector `g` tales que `u = T · u_libre + g`: en filas correspondientes a DOF libres, `T` es la matriz identidad rectangular; en filas esclavas, lleva los coeficientes `α_si` en las columnas de los maestros. El sistema entregado al subsistema algebraico es `K_red = TᵀKT`, `F_red = Tᵀ(F − K · g)`. La imposición es exacta a redondeo, preserva la simetría y la positividad definida de `K`, y deja el sistema reducido apto para solvers iterativos (CG, GMRES).

La forma afín cubre con la misma maquinaria los casos de empotramiento, asentamiento prescrito, apoyo en plano oblicuo, periodicidad de celda unitaria y unión rígida master-slave entre nodos. Las restricciones declaradas como cadenas (un esclavo cuyo maestro es a su vez esclavo de otra restricción) se resuelven por cierre transitivo en el momento de construir `T`; los ciclos y las redeclaraciones inconsistentes se detectan en validación temprana. La declaración se hace mediante `Domain.add_linear_constraint` o desde el bloque `linear_constraints` del archivo YAML.

## Compilación Just-In-Time mediante Numba

Las funciones críticas (ensamblaje elemento→global, núcleo del algoritmo de retorno radial) se decoran con `@njit` y se compilan a código nativo en su primera invocación. La primera ejecución asume el coste de compilación; las invocaciones siguientes se ejecutan a velocidad cercana a la de un programa Fortran compilado. La restricción consiste en que el código compilado solo admite tipos primitivos y arreglos NumPy, no objetos arbitrarios de Python.

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
