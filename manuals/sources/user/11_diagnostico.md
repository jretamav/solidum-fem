# Diagnóstico de Problemas

- **`YamlValidationError: parámetro 'X' no aceptado por 'Y'.`**
  El parser detecta un parámetro que el constructor del elemento no admite. La salida lista los parámetros aceptados; revisar el catálogo o la spec del componente.

- **`MechanismError: El modelo no está suficientemente apoyado para un análisis estático…`**
  Los apoyos no impiden que el modelo se mueva como sólido rígido, así que el problema estático no tiene solución única. Solidum lo detecta **antes de resolver** y dice exactamente qué movimiento queda libre, por ejemplo *"traslación en la dirección x"* o *"giro alrededor del eje z que pasa por (0, 0.5)"*. Añadir apoyos (o restricciones lineales) que impidan esos movimientos. Ojo con los rodillos alineados: impedir solo `uy` a lo largo de una recta vertical deja libre el giro alrededor de cualquier punto de esa recta. En un modelo térmico estacionario, el mensaje equivalente es que el campo `T` no tiene ningún valor prescrito. En análisis modal o dinámico un modelo libre sí es válido y no se comprueba.

- **`IllPosedSystemError: Resultado descartado…`**
  El modelo está bien apoyado en conjunto, pero la solución no es fiable: o tiene un **mecanismo interno**, o el solver no consiguió resolver el sistema (*error hacia atrás normado* por encima de 10⁻¹⁰, algo que no debería pasar con un modelo bien planteado). Un mecanismo interno es una parte que puede moverse sin deformar ningún elemento: una rótula o articulación de más, una barra o elemento suelto, o dos partes unidas por un solo nodo que pueden girar entre sí. El mensaje dice qué elementos se mueven y, si lo hacen como un sólido rígido, cómo: *"la parte formada por los elementos 1025–2048 puede moverse sin deformarse: giro alrededor del eje z que pasa por (1, 1)"*. Si el modelo es correcto pero tiene rigideces muy desproporcionadas (zonas rígidas, un material casi incompresible), no es un mecanismo: su modo de menor energía deforma algún elemento y no se rechaza por esto. Más raramente, rigideces muy desproporcionadas entre partes del modelo (penalizaciones, unidades mezcladas). Revisar la conectividad. Antes de esta comprobación, un modelo así devolvía desplazamientos absurdos sin ningún aviso. Si el mensaje dice que *el solver directo no pudo factorizar la matriz*, falta rigidez por completo en algún grado de libertad: un nudo al que sólo llega una barra en línea recta, un grado de libertad que ningún elemento usa, un giro sin rigidez.

- **`IllPosedSystemError: Análisis descartado (NonlinearSolver): la matriz tangente inicial (U = 0)…`**
  Lo mismo en un análisis no lineal (`NonlinearSolver`, `ArcLengthSolver`, `IndirectDisplacementSolver`): la matriz tangente inicial es singular, así que el modelo tiene un mecanismo interno, aunque la carga no lo active, y el mensaje dice qué elementos se mueven. Si ningún elemento es geométricamente no lineal, la rigidez del modelo no puede crecer al deformarse (la no linealidad del material sólo la reduce) y el mecanismo lo es durante todo el análisis: por eso se rechaza antes de iterar. Para ver la forma del mecanismo, hacer un análisis modal del mismo modelo con densidad en los materiales: cada mecanismo aparece como una frecuencia nula y su forma modal muestra qué parte se mueve.

- **Aviso:** *la matriz tangente inicial (U = 0) es singular: el modelo tiene N mecanismo(s) interno(s)…*
  El mismo caso con elementos geométricamente no lineales (corrotacionales, cables) o materiales unilaterales. La rigidez geométrica, o la de un cable al tensarse, puede estabilizar el mecanismo al deformarse —una red de cables sin pretensado es el ejemplo típico—, así que el análisis sigue. Si el modelo no es de ese tipo, revisar la conectividad de los elementos que nombra el aviso.

- **Aviso:** *sistema mal condicionado (κ̃ ≈ …): la solución puede haber perdido hasta N de sus unas 16 cifras significativas…*
  El análisis lineal se resolvió, pero el sistema está tan mal condicionado que la solución puede tener pocas cifras correctas. El ordenador trabaja con unas 16 cifras; con un número de condición κ, la solución conserva unas $16 - \log_{10}\kappa$, y en la práctica una o dos más. Solidum lo estima en todo análisis lineal estático y térmico estacionario (queda en `result.condition_estimate`) y avisa por encima de $10^{12}$, cuando la solución baja a unas 5 cifras. Causas habituales: rigideces muy desproporcionadas (zonas rígidas de rigidez ficticia enorme, penalizaciones), un coeficiente de Poisson muy cercano a 0.5, o una barra dividida en demasiados elementos: en una viga κ crece con la cuarta potencia del número de elementos, y un voladizo de 10 000 elementos de Euler-Bernoulli da una flecha errónea en un 50 %. Remedios: rigideces de zona rígida unas $10^3$–$10^4$ veces la del resto (no $10^8$), una restricción multipunto en vez de una zona rígida, o menos elementos por barra (un elemento de Euler-Bernoulli por tramo ya da los nodos exactos con cargas nodales). Con un material casi incompresible el remedio sería una formulación mixta, que Solidum no implementa todavía: mientras tanto, conviene un $\nu$ no más allá de 0.499. No depende de las unidades.

- **Aviso:** *con el solver algebraico iterativo (linear_algebra: iterative) no se comprueba si el modelo tiene un mecanismo interno…*
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
