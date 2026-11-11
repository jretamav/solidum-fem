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

- **`StepConstraintError` en `ArcLengthSolver` (o en las variantes del arco)**
  La restricción del paso no tiene solución aunque se reduzca la longitud del paso: con el arco cilíndrico, las raíces de su ecuación son imaginarias («Raíces imaginarias» en el detalle) porque el arco no corta la trayectoria de equilibrio. Ocurre en un retroceso (snap-back) por localización o tras un punto límite agudo. Usar `IndirectDisplacementSolver`, con una magnitud que crezca de forma monótona (el alargamiento de la zona que se ablanda, la apertura de la grieta), o un primer paso menor. El detalle del mensaje dice también el último paso convergido; el mismo dato queda en `solver.lambda_final` y `solver.U_final`.

- **Aviso:** *trazado detenido en λ=… < max_lambda=… tras agotar max_steps*
  El trazado no terminó: el resultado es el del último paso convergido y `solver.termination` vale `"max_steps"` (`SolveResult.converged` es `False`). Aumentar `max_steps` o el tamaño del paso.

- **`TimeIntegrationInstabilityError` en un análisis en el tiempo**
  El paso de tiempo supera el límite de estabilidad del esquema: diferencia central por encima de la condición CFL ($\Delta t < 2/\omega_{\max}$) o θ-método con $\theta < 1/2$. La solución crece sin límite y el análisis se detiene. Reducir `dt` o usar un esquema incondicionalmente estable (Newmark con $\beta = 1/4$, $\gamma = 1/2$; HHT-α; $\theta \ge 1/2$). El historial hasta el último paso estable está en `exc.partial_result`.

- **`MasslessDofsError` en un transitorio**
  Hay grados de libertad libres sin masa (sin densidad en los elementos que los comparten ni masa de nudo, o sin inercia rotacional), y el mensaje los nombra («nodo 2 (ux)»). Los integradores implícitos los admiten (aceleración inicial consistente, ADR 0024 §6 bis): el error sólo aparece si además forman un mecanismo (`massless_mechanism`: sin masa y sin rigidez que los fije), y se corrige revisando la conectividad y los apoyos. La diferencia central (`explicit_integration`) exige masa en todos: dar masa a esos grados de libertad (`density` o `add_nodal_mass`) o usar un integrador implícito.

- **Aviso:** *la carga aplicada es un np.ndarray simple, sin sus cargas equivalentes por elemento…*
  Las fuerzas internas de los elementos con carga distribuida pueden salir mal: la carga perdió sus cargas equivalentes al convertirla con `np.asarray`/`np.array` o al sumarla con `+=` sobre un arreglo de ceros, y se ensambló más de una carga de cuerpo. Combinar las cargas que devuelven `assemble_self_weight` y `assemble_body_load` con `+` y escalares (`F = F_a + λ·F_b`). Lo mismo avisa el historial de fuerzas internas de un transitorio cuyo `F_func` devuelve un `np.ndarray` simple.

- **`ValueError`: *la longitud característica h = … supera el límite 2·G_f·E/ft²***
  Con la energía de fractura `G_f`, la rama de ablandamiento en tracción termina en $\varepsilon_{tu} = 2G_f/(f_t h)$. Si el punto de integración es tan largo que $\varepsilon_{tu}$ quedaría por debajo de la deformación de agrietamiento $f_t/E$, la rama tendría retroceso en el propio material. El elemento lo rechaza al construirse y dice cuál es. Usar elementos más cortos.

- **`ValueError`: *ManderConcrete1D con G_f necesita la longitud característica del punto de integración***
  El material regularizado se usa en un elemento que no le entrega la longitud de sus puntos. La entregan las barras y los marcos con sección por fibras. Con otro elemento, dar la rama descendente en deformaciones (`eps_tu`, `eps_cu`, `eps_sp`).

- **Un transitorio con concreto deja de converger en la vibración libre.**
  El cierre de una grieta cambia bruscamente la rigidez de la fibra (de la secante agrietada a la de compresión) y el Newton puede rebotar entre las dos ramas. Activar `line_search: true` en el solver, o reducir `dt`.

- **Aviso:** *la etapa empieza fuera de equilibrio (‖F_c − F_int(u0)‖/‖F‖ = …)*
  Una etapa (`solidum.run(..., previous=…)` o `solve(F, u0=…, constant_load=…)`) no parte del estado convergido de la anterior: otro `u0`, otra carga constante o un `Assembler` distinto, que no guarda el estado interno de los elementos. Continuar con `previous=resultado` y el mismo `Assembler`. La variante *el modelo conserva el estado interno de un análisis anterior, pero esta corrida empieza en U = 0* avisa de un análisis nuevo sobre un modelo ya cargado: para continuarlo, `previous`; para empezar de cero, construir el modelo de nuevo.

- **`ValueError`: *con 'previous' pase también el 'assembler' de la etapa anterior…***
  El `Assembler` guarda el estado interno de los elementos (plasticidad, daño, grietas); uno nuevo partiría del estado virgen. Pasar a cada etapa el mismo `assembler`.

- **`ValueError`: *… no admite etapas de carga (u0, constant_load)…***
  Sólo los solvers incrementales estáticos (`NonlinearSolver`, `ArcLengthSolver`, `IndirectDisplacementSolver`) continúan una etapa con carga constante. Un `LinearSolver` sobre una etapa es una perturbación lineal (ADR 0031) y los transitorios continúan el estado en el tiempo (ADR 0029).

- **`ValueError`: *'previous' es una perturbación lineal…***
  El resultado de un `LinearSolver` sobre una etapa no cambia el estado de los elementos, así que no se puede continuar desde él: la etapa siguiente debe tomar como `previous` la etapa sobre la que se calculó la perturbación.

- **`ValueError`: *una etapa estática no puede continuar un análisis en el tiempo…***
  Con `stages`, una etapa estática sólo puede seguir a otra estática: qué carga se mantiene y qué pasa con la velocidad al pasar de un análisis en el tiempo a uno estático son decisiones físicas todavía pendientes. Las etapas modal, armónica y espectral sí pueden seguir a cualquiera (linealizan alrededor de su estado).

- **`YamlValidationError` en `stages`** (*stages: con etapas, el solver y las cargas se declaran en cada etapa…*, *stages[i]: falta el bloque 'solver'*)
  Con `stages`, el `solver` y las cargas (`point_loads…`, `gravity`, `body_force`, `tractions…`, `thermal_loads`) van dentro de cada etapa (`loads` y `solver`); declararlos además en el nivel raíz es un error. Cada etapa admite sólo `name`, `loads` y `solver`.

- **Aviso:** *ModalSolver: N autovalor(es) negativo(s)… el estado es inestable en esos modos*
  Un modal sobre el estado de una etapa (ADR 0031) usa su rigidez tangente; un autovalor negativo indica que ese estado es inestable en ese modo (por ejemplo, una carga axial por encima de la de pandeo). Los modos quedan en `ModalResult.unstable_modes`, con su autovalor en `eigenvalues` y frecuencia y periodo `NaN`. Un `ResponseSpectrumSolver` sobre ese estado se detiene con `ValueError` (*el estado sobre el que se calcula tiene modos inestables*): revisar la carga de la etapa.

- **`ValueError` de `EmbeddedReinforcement`** (refuerzo embebido, bloque `constraints` o `solidum.embed_reinforcement`)
  - *el punto … de la barra no está dentro de ningún elemento anfitrión* / *la barra sale de los anfitriones cerca de …* / *el punto medio … de un tramo no está dentro de ningún anfitrión*: la barra debe quedar entera dentro de los sólidos anfitriones (o de los `hosts` indicados); revisar sus `points` o ampliar `hosts`.
  - *el dominio no tiene elementos de sólido 2D/3D donde embeber la barra* / *los elementos … no son sólidos isoparamétricos*: la barra se embebe en sólidos con funciones de forma (Tri3, Quad4, Hex8…); la llamada va después de crear esos elementos.
  - *… no es un elemento 1D de N nodos con traslaciones …* / *… tiene N nodos; la barra admite elementos de línea de dos o tres nodos*: `element` debe ser una barra (`Truss2D`, `Truss3D`, `Truss2D3`, `Truss3D3`…), de la dimensión del sólido.
  - *Bond: law debe ser un material de interfaz de línea (GEOMETRY='line')*, *Bond: diameter debe ser > 0*, *bond sin 'law'* o *la ley de adherencia … no admite N componentes*: la ley de `bond` debe ser una ley de interfaz de **línea** (`BondSlipModelCode`, `BondSlipMultilinear`) declarada en `interface_materials`, con `diameter` > 0; una ley de grieta (de superficie) no sirve.
  - *bond_forces: con adherencia-deslizamiento no hay restricciones de embebido…* / *bond_state: la barra tiene adherencia perfecta*: con `bond` las fuerzas de adherencia se leen con `EmbeddedBar.bond_state`; sin `bond`, con `bond_forces`.

- **Aviso:** *CST_Embedded2D: el Newton local del salto no convergió en N iteraciones…*
  El salto de un elemento con discontinuidad embebida no se pudo despejar en el paso: el cálculo sigue con el residuo local, que el Newton global ve como desequilibrio. Si se repite, reducir el paso. Con la ley de modo mixto `CohesiveDamageMixedMode`, la causa típica es un elemento demasiado grande para el ablandamiento en cortante: el sistema local sólo tiene solución única si $\beta^2 |T'|\,h < G$ (sección «Materiales de interfaz» del capítulo «Catálogo de Modelos Constitutivos»); refinar la malla en la zona de la grieta.

- **Niveles del registro.**
  Con el nivel por omisión (`INFO`) cada solver incremental escribe una línea por paso convergido. El detalle por iteración de Newton sale con `solidum.set_log_level("DEBUG")`. Un `WARNING` de bisección («sin convergencia (tangente singular); bisección…») es un evento que el solver resolvió solo; un `ERROR` significa siempre que el análisis se detuvo.

- **Consola de Windows (cp1252) o ASCII.**
  Si la consola no tiene un carácter de un mensaje (letras griegas, subíndices, flechas), Solidum lo escribe en texto plano: `Δλ₁` sale como `Deltalambda_1` y `⇒` como `=>`. No es un error. Con una consola UTF-8 (`PYTHONIOENCODING=utf-8`) los mensajes salen tal cual.

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
