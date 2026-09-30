# Diagnóstico de Problemas

- **`YamlValidationError: parámetro 'X' no aceptado por 'Y'.`**
  El parser detecta un parámetro que el constructor del elemento no admite. La salida lista los parámetros aceptados; revisar el catálogo o la spec del componente.

- **`MechanismError: El modelo no está suficientemente apoyado para un análisis estático…`**
  Los apoyos no impiden que el modelo se mueva como sólido rígido, así que el problema estático no tiene solución única. Solidum lo detecta **antes de resolver** y dice exactamente qué movimiento queda libre, por ejemplo *"traslación en la dirección x"* o *"giro alrededor del eje z que pasa por (0, 0.5)"*. Añadir apoyos (o restricciones lineales) que impidan esos movimientos. Ojo con los rodillos alineados: impedir solo `uy` a lo largo de una recta vertical deja libre el giro alrededor de cualquier punto de esa recta. En un modelo térmico estacionario, el mensaje equivalente es que el campo `T` no tiene ningún valor prescrito. En análisis modal o dinámico un modelo libre sí es válido y no se comprueba.

- **`IllPosedSystemError: Resultado descartado…`**
  El modelo está bien apoyado en conjunto, pero la solución no es fiable: o no satisface el equilibrio, o la matriz tiene pivotes numéricamente nulos. Casi siempre es un **mecanismo interno**: una rótula o articulación de más, una barra o elemento suelto, o dos partes unidas por un solo nodo que pueden girar entre sí. Más raramente, rigideces muy desproporcionadas entre partes del modelo (penalizaciones, unidades mezcladas). Revisar la conectividad. Antes de esta comprobación, un modelo así devolvía desplazamientos absurdos sin ningún aviso. Si el mensaje dice que *el solver directo no pudo factorizar la matriz*, falta rigidez por completo en algún grado de libertad: un nudo al que sólo llega una barra en línea recta, un grado de libertad que ningún elemento usa, un giro sin rigidez.

- **`IllPosedSystemError: Análisis descartado (NonlinearSolver): la matriz tangente inicial (U = 0)…`**
  Lo mismo en un análisis no lineal (`NonlinearSolver`, `ArcLengthSolver`, `IndirectDisplacementSolver`): la matriz tangente inicial es singular, así que el modelo tiene un mecanismo interno, aunque la carga no lo active. Si ningún elemento es geométricamente no lineal, la rigidez del modelo no puede crecer al deformarse (la no linealidad del material sólo la reduce) y el mecanismo lo es durante todo el análisis: por eso se rechaza antes de iterar. Para localizarlo, hacer un análisis modal del mismo modelo con densidad en los materiales: cada mecanismo aparece como una frecuencia nula y su forma modal muestra qué parte se mueve.

- **Aviso: *la matriz tangente inicial (U = 0) tiene N pivote(s) numéricamente nulo(s)…***
  El mismo caso con elementos geométricamente no lineales (corrotacionales, cables) o materiales unilaterales. La rigidez geométrica, o la de un cable al tensarse, puede estabilizar el mecanismo al deformarse —una red de cables sin pretensado es el ejemplo típico—, así que el análisis sigue. Si el modelo no es de ese tipo, localizar el mecanismo con un análisis modal.

- **Aviso: *con el solver algebraico iterativo (linear_algebra: iterative) no se comprueba si el modelo tiene un mecanismo interno…***
  El solver iterativo no factoriza la matriz y no puede delatar un mecanismo interno que la carga no active: la solución cumpliría el equilibrio, pero con una parte arbitraria del mecanismo. El aviso sale una vez por análisis estático. Para comprobarlo, un análisis modal con densidad o una ejecución con `linear_algebra: auto`.

- **`SingularTangentError` en un análisis no lineal**
  La rigidez tangente perdió rango: un punto límite o bifurcación (usar `ArcLengthSolver`), o un mecanismo (ver los puntos anteriores). Si el mensaje añade *"el sistema tangente no se pudo resolver con precisión"*, la causa más probable es un mecanismo.

- **`MemoryError: Memoria insuficiente para factorizar…`**
  El modelo es demasiado grande para el solver directo, cuya memoria crece más deprisa que el modelo. Usar `linear_algebra: iterative` en el bloque `solver` (ver el anexo de la capa algebraica).

- **Cables y materiales con ablandamiento.**
  En modelos con cables o materiales con softening pronunciado, un elemento que se destensa o se degrada puede dejar un mecanismo: comprobar que existen caminos de carga alternativos.

- **`Newton-Raphson no convergió / bisecando incremento.`**
  El paso de carga es demasiado grande. Aumentar `num_steps`, asegurar `adaptive: true`. Verificar tipeo de los parámetros constitutivos ($\sigma_y$, $\kappa_0$, etc.). Si el problema tiene snap-through o softening, cambiar a `ArcLengthSolver`.

- **`No se importó ningún elemento (Solid2D).`**
  El parser de Gmsh no encontró superficies bidimensionales. En Gmsh, crear una *Plane Surface* y generar la malla 2D antes de exportar.

- **`Cable totalmente destensado — KT = 0.`**
  Pretensar el cable o garantizar que la estructura tiene rigidez residual por otros elementos. Considerar también pasos de carga finos para capturar la transición tensado → destensado.
