# ADR 0019 — Red de seguridad del análisis estático: mecanismos, equilibrio y pivotes nulos

- **Estado**: aceptado; ampliado el 2026-09-30 (§Cierres del 2026-09-30) y el 2026-10-20 (§Ampliación del 2026-10-20: mecanismos internos por la energía del modo)
- **Fecha**: 2026-09-23
- **Alcance**: `solidum/math/solvers/model_checks.py` (nuevo); `MechanismError` e `IllPosedSystemError` en `diagnostics.py`; los tres solvers estáticos (`LinearSolver`, `NonlinearSolver`, `ArcLengthSolver` y su variante por disipación); telemetría del residuo lineal en el corrector (ADR 0015); `n_zero_pivots` en SuperLU y Pardiso; mensaje de memoria insuficiente en los backends directos; `rigid_body_basis` y `Assembler.rigid_basis` (reutilizan los modos del ADR 0018). **No cambia** ninguna formulación ni ningún resultado de un modelo bien planteado; cambia qué pasa con los mal planteados.

## Contexto

El usuario planteó el supuesto que define este ADR: el usuario de Solidum no tiene por qué dominar los problemas numéricos de la solución del sistema algebraico, y sería arriesgado pedirle que elija el solver según el tipo y el tamaño del problema. ¿Qué red de seguridad tiene Solidum?

La primera mitad de la respuesta ya existía: el usuario **no elige**. Desde el ADR 0003 el despachador decide por las propiedades que declaran materiales y elementos, siempre un solver directo, y `linear_algebra` es una herramienta de diagnóstico; el iterativo (ADR 0018) nunca se activa solo. Dentro de un Newton, además, el equilibrio se verifica siempre con el residuo real.

La segunda mitad se midió, con el error de modelado más común de un usuario no experto: un **mecanismo** (apoyos insuficientes). El resultado fue un agujero en el camino por defecto:

| Backend | Análisis estático lineal | Análisis no lineal |
|---|---|---|
| SuperLU | **devuelve en silencio \|u\| = 5·10⁶ m** | `UnknownDivergenceError` ("modo no clasificado") |
| Pardiso (el default desde el ADR 0017) | **devuelve en silencio \|u\| = 7·10⁴ m** | `UnknownDivergenceError` |
| Iterativo | falla con la causa | `SingularTangentError` |

Un solver directo ante una matriz singular no avisa. SuperLU arrastra un pivote de 10⁻¹⁶ y Pardiso lo **perturba por diseño** y sigue adelante. ANSYS y Abaqus sí lo detectan: avisan de singularidad numérica o pivote nulo e indican el nodo y el grado de libertad.

## Decisión: tres capas

### Capa 1 — mecanismos rígidos, antes de resolver

`ensure_statically_restrained(assembler)` se ejecuta al empezar los tres solvers estáticos. Reutiliza los modos de cuerpo rígido del ADR 0018 (derivados del nombre de los DOF, núcleo exacto de `K` sin apoyos) y comprueba cuáles dejan libres las restricciones.

Un modo `r` es compatible con las restricciones homogéneas si `r = T·r_libre` (ADR 0004), es decir, si se anula `S·c = (B − T·B_libre)[esclavos]·c`. El **núcleo de `S`** son las combinaciones de modos que ningún apoyo impide, incluidas las que no son un modo puro. Un pasador único deja libre el giro **alrededor del pasador**, no alrededor del centroide: esa rotación es una combinación de giro y traslación, y el núcleo de `S` la encuentra. Las restricciones lineales (MPC) cuentan como apoyo sin tratamiento especial, porque entran por `T`.

Cada movimiento libre se describe en términos del modelo, no de pivotes, y se lanza `MechanismError` (un `ValueError`: es un error de modelado):

```
El modelo no está suficientemente apoyado para un análisis estático (LinearSolver):
puede moverse como sólido rígido sin deformarse [...]. Movimientos libres:
  - giro alrededor del eje z que pasa por (0, 0.5)
  - traslación en la dirección x
Añade apoyos (o restricciones lineales) que impidan esos movimientos. [...]
```

Un movimiento `u(x) = t + w × (x − c)` se describe como traslación si `|w|·L` es despreciable frente a `|t|`. Si no, como giro alrededor del eje `w` que pasa por `p = c + (w × t)/|w|²`, el punto donde el movimiento es paralelo al eje; en 3D se añade el deslizamiento a lo largo del eje si existe. Los campos escalares (`T`) se tratan por bloques separados, porque las restricciones no los acoplan con los mecánicos: *"el campo 'T' no tiene ningún valor prescrito"*.

**Solo en análisis estáticos.** En uno modal o dinámico un modelo libre es un problema bien planteado (la masa o la capacidad lo regularizan) y la comprobación no se hace.

El detector encontró más de lo que yo le había descrito al usuario. En el caso "rodillos en `uy` a lo largo de `x = 0`" hay **dos** movimientos libres: la traslación en x y el giro alrededor de cualquier punto de esa recta, que mueve esos nodos en horizontal.

### Capa 2 — equilibrio, después de resolver (estático lineal)

`check_linear_solution` exige `‖F − K·u‖ ≤ EQUILIBRIUM_RTOL·‖F‖` con `EQUILIBRIUM_RTOL = 1e-8`. Medido: un sistema bien planteado queda en ~10⁻¹³ con cualquier directo y en ≤ 10⁻¹⁰ con el iterativo; un mecanismo cargado da ~1. Deja cinco órdenes de margen a cada lado. Atrapa lo que la capa 1 no ve: **mecanismos internos** (una rótula de más, dos partes unidas por un solo nodo) y matrices tan mal condicionadas que la solución carece de sentido. Lanza `IllPosedSystemError`.

### Capa 3 — pivotes numéricamente nulos (estático lineal)

> *Desde el 2026-10-20 la capa 3 no decide por el tamaño del pivote: un pivote pequeño frente a la diagonal sólo dispara la búsqueda del modo de menor energía, y decide su energía por elemento (§Ampliación del 2026-10-20). Lo que sigue es la decisión original.*

Un mecanismo interno **sin carga que lo active** da un sistema singular pero consistente: el residuo sale perfecto (10⁻¹⁴) y la capa 2 no lo ve. Lo delatan los pivotes. Ambos directos exponen ahora `n_zero_pivots` con **el mismo criterio**, `|pivote| < ZERO_PIVOT_RTOL·max|K|` con `ZERO_PIVOT_RTOL = 1e-13`, que es el umbral con el que MKL Pardiso decide perturbar un pivote (`iparm(10) = 13`):

- **Pardiso**: `iparm(14)`, pivotes perturbados.
- **SuperLU**: la diagonal de `U`. SciPy sólo la expone copiando el factor, así que se calcula de forma perezosa, una vez por factorización y sólo cuando la pide el análisis estático lineal.

Medido: pivote mínimo relativo 3·10⁻³ en un modelo bien planteado y 1,6·10⁻¹⁶ con la rótula interna.

### Por qué las capas 2 y 3 solo en el estático lineal

Dentro de un Newton o de un arc-length, resolver un sistema casi singular es **legítimo**: cerca de un punto límite la tangente lo es por definición, y el arc-length la usa. Rechazarlo rompería el algoritmo. Ahí la red es el propio Newton, que exige equilibrio real. Lo que se añade es solo **diagnóstico**: el corrector registra el residuo relativo de cada resolución con factorización fresca (con la tangente congelada del Newton modificado no aplica) y, si el paso fracasa con un residuo lineal por encima de `EQUILIBRIUM_RTOL`, la divergencia se clasifica como tangente singular y el mensaje lo dice. El algoritmo no cambia.

### Memoria insuficiente

Los backends directos traducen un fallo de memoria (MemoryError de SuperLU, código −2 de la MKL) a un `MemoryError` con la salida: *"usa `linear_algebra: iterative`"*. No se cambia de backend automáticamente: pasar al iterativo sin que el usuario lo sepa le introduciría justamente los riesgos de convergencia que esta red intenta evitarle.

## Resultado

Mecanismo rígido (rodillos alineados, pasador único): `MechanismError` **antes de resolver**, con los seis backends × solvers probados (SuperLU, Pardiso, iterativo; lineal, no lineal, arc-length).

Mecanismo interno (dos bloques unidos por un nodo):

| | SuperLU | Pardiso | Iterativo |
|---|---|---|---|
| Cargado, lineal | `IllPosedSystemError` (residuo) | `IllPosedSystemError` (pivote) | `IterativeNotConvergedError` |
| Sin carga que lo active, lineal | `IllPosedSystemError` (pivote) | `IllPosedSystemError` (pivote) | **devuelve una solución** |
| Cargado, no lineal | `SingularTangentError` (antes: "no clasificado") | `SingularTangentError` | `SingularTangentError` |

*Desde el 2026-09-30* cambian dos filas: con el iterativo, la solución sin carga que active el mecanismo llega con un aviso, y con los directos el no lineal rechaza el mecanismo interno de este modelo (geométricamente lineal) en la tangente inicial, con `IllPosedSystemError`, antes de iterar (con el iterativo, aviso y después `SingularTangentError`). Ver §Cierres del 2026-09-30.

**Cero falsos positivos**: la suite completa (1 479 pruebas con Pardiso y 1 471 con SuperLU, que ejercita su comprobación de pivotes) pasa sin cambios. Las tres capas se ejecutan en todos los análisis estáticos de la suite: marcos, armaduras, cables corrotacionales, discontinuidades embebidas, restricciones lineales, térmico estacionario, 3D.

**Coste** (Hex8 40³, 2·10⁵ DOF, frente a 17 s de factorizar y resolver con Pardiso): capa 1, 0,42 s (2,5 %); capas 2 y 3, 30 ms. La capa 1 bajó de 0,85 s al corregir dos cosas: una SVD completa que construía una matriz densa de (restricciones × restricciones), unos 200 MB con 5 000 restricciones, sustituida por QR económica + SVD del factor `R` de 6×6; y un recorrido por nodos repetido por cada nombre de DOF.

## Cierres del 2026-09-30

Tres huecos de la lista original dejaban pasar sin aviso un **mecanismo interno que la carga no activa**, el caso que sólo delatan los pivotes (el residuo sale perfecto). El caso de referencia es el último del ejemplo 11 del manual de ejemplos: la armadura Pratt sin la diagonal del segundo panel, cargada sólo en horizontal.

### 1. Cholesky cuenta sus pivotes nulos (deuda #28)

`CholeskyFactorized.n_zero_pivots` cuenta los pivotes `d_i` de `LDLᵀ` por debajo de `ZERO_PIVOT_RTOL·max|K|`. Para una matriz simétrica y la misma ordenación, los pivotes de `LDLᵀ` son los de la eliminación gaussiana (`d_i = U_ii`), así que la regla es la de SuperLU; los valores concretos dependen de la ordenación de cada backend, pero no el número de pivotes nulos de una matriz singular. Pardiso compara con la matriz escalada, no con `max|K|` (ver «Lo que queda sin cubrir»). Se leen de la diagonal del factor sin copiarlo: `Factor.D()` en scikit-sparse 0.4 (`L_ii²` si el factor es `LLᵀ`) y la vista `CholeskyFactor.factor` en la 0.5 (`L_ii²` o `D_ii` según `is_ll`).

Al revisar la dependencia apareció un problema mayor: **scikit-sparse 0.5 (2025) rehízo su API** —`cholesky(A)` devuelve ahora las matrices `(R, p)`, no un factor—, y el extra `cholmod` (`scikit-sparse>=0.4`) instalaría hoy la 0.5, con la que el backend fallaba en la primera resolución. El backend admite ahora las dos versiones. Con la 0.5, un pivote diminuto no se detecta al factorizar sino al resolver, y sólo si su estimación de condición es exactamente 0; esa comprobación se adelanta a la factorización para que el fallback a LU actúe en el mismo sitio con las dos versiones, y se silencia su aviso *Matrix is nearly singular*, que salta con tangentes legítimamente mal condicionadas y no distingue una matriz singular.

scikit-sparse no se puede instalar en la máquina Windows de desarrollo. La lógica se prueba con un `sksparse.cholmod` simulado que imita las dos API sobre una `LDLᵀ` densa, y la prueba real (`skipUnless`) queda para un job nuevo del CI de Linux, `test-cholmod`, que compila scikit-sparse contra `libsuitesparse-dev` y falla si CHOLMOD no queda activo. Es un job aparte, y no parte del de backends opcionales, porque **el despachador prefiere Cholesky a Pardiso** para toda matriz simétrica definida positiva (Cholesky → Pardiso → LU, ADR 0017): instalarlo allí sacaría a Pardiso de la mayor parte de la suite.

### 2. Tangente inicial de los solvers estáticos no lineales (decisión del usuario)

`NonlinearSolver`, `ArcLengthSolver`, `IndirectDisplacementSolver` y el arco por disipación del módulo de usuario `dissipation_control` inspeccionan los pivotes nulos de la **tangente en `U = 0`**. Lo que se hace con ellos lo decidió el usuario:

- **Modelo geométricamente lineal** ⇒ `IllPosedSystemError`. En pequeños desplazamientos la matriz `B` es constante; si `v` es un mecanismo de `K(0)`, la energía `∫(Bv)ᵀ D₀ (Bv)` es nula y, con `D₀` definida positiva, `Bv = 0` en todo punto, así que `K(U)·v = 0` para todo `U` y `F_int(U + s·v) = F_int(U)`: el mecanismo lo es durante todo el análisis, porque la no linealidad del material sólo reduce la rigidez. Si la carga no lo activa, el Newton convergía con una parte arbitraria del mecanismo; si lo activa, no hay equilibrio.
- **Con elementos geométricamente no lineales** ⇒ sólo un aviso en el log. La rigidez geométrica crece con las fuerzas internas y puede estabilizar el mecanismo al deformarse (una red de cables sin pretensado).

La distinción la da un **contrato declarativo nuevo**, `Element.GEOMETRICALLY_NONLINEAR` (`False` por omisión; `True` en `Truss2DCorot`, `Truss3DCorot`, `Cable2DCorot`, `Cable3DCorot` y `Frame2DEulerCorot`). El barrido de contrato de elementos lo exige booleano y lo contrasta con la cinemática: con material elástico, la rigidez de un elemento geométricamente lineal no cambia tras un giro rígido pequeño y la de uno corrotacional sí. El elemento embebido del módulo `discontinuities` es geométricamente lineal y hereda el valor por omisión.

**Refinamiento de la IA, confirmado por el usuario el 2026-09-30**: la demostración anterior exige `D₀` definida positiva, y un material unilateral (`Material.IS_UNILATERAL`, hoy `CableMaterial1D`) tiene rigidez nula sin deformar y la gana al traccionarse. Un marco geométricamente lineal con una sección por fibras que mezcle tendones de cable puede, en teoría, ganar rigidez respecto de `U = 0`. Esos modelos también reciben sólo el aviso. En la práctica el caso es marginal: si las fibras bilaterales no bastan para que `K(0)` sea regular, suele ser exactamente singular y el Newton no puede ni empezar.

**Reutiliza la factorización** que el solver ya hace en `U = 0` (la primera iteración del Newton o el predictor del primer arco), en vez de factorizar dos veces. El corrector compartido (ADR 0015) gana una pieza genérica, `NewtonCorrector.inspect_next_factorization(inspector)`: entrega su próxima factorización fresca al inspector antes de usarla, o `None` si la factorización directa no se completó por un pivote exactamente nulo. `model_checks.check_initial_tangent(assembler, corrector, solver)` arma el inspector al empezar el análisis, justo después de `ensure_statically_restrained`, y es pública: el solver del módulo de usuario la llama igual. Dos sitios que traducían todo `RuntimeError` en "tangente singular, bisecar" dejan pasar `IllPosedSystemError`: el bucle del corrector y el predictor del arco. Medido: el mismo número de factorizaciones con la comprobación que sin ella. Su coste es el de contar los pivotes una vez por análisis: nulo con Pardiso (4 µs, `iparm(14)`) y una copia del factor con SuperLU (2-3 % de una factorización: 30 ms de 1,05 s a 6 084 DOF, 95 ms de 4,65 s a 13 872).

Los solvers dinámicos y el modal no hacen la comprobación: ahí un modelo libre o con mecanismos es legítimo, porque la masa lo regulariza. Las tangentes posteriores a `U = 0` tampoco se inspeccionan, por la razón de siempre: cerca de un punto límite una tangente casi singular es el algoritmo.

### 3. El backend iterativo lo dice

Un método de Krylov no factoriza y no tiene pivotes. Distinguir un mecanismo (autovalor exactamente nulo) de un modelo sólo mal condicionado exigiría estimar el menor autovalor con una precisión relativa de ~10⁻¹⁴, que un método de Krylov no da de forma fiable. En lugar de fingir que lo detecta, el hueco deja de ser silencioso: `IterativeFactorized.n_zero_pivots = None` declara que no hay pivotes, y `model_checks.warn_if_pivots_unavailable` avisa en el log, una vez por análisis estático (`LinearSolver` al factorizar, los no lineales en la tangente inicial), de que un mecanismo interno no activado no se comprueba y de cómo comprobarlo: un análisis modal con densidad o el solver directo.

### Y una factorización imposible

Al probar un nudo al que sólo llega una barra recta apareció otro caso sin explicar: una fila de `K` vacía hace que SuperLU lance *Factor is exactly singular* y que Pardiso rechace la matriz por filas vacías, y `LinearSolver` dejaba salir ese `RuntimeError` del backend. Ahora es `IllPosedSystemError` con el mensaje del backend adjunto, y en la tangente inicial de los no lineales, el mismo veredicto que un pivote nulo.

### Resultado

| Mecanismo interno sin carga que lo active | SuperLU | Pardiso | Cholesky | Iterativo |
|---|---|---|---|---|
| `LinearSolver` | error (pivote) | error (pivote) | error (pivote) | **aviso** y solución |
| No lineal, geométricamente lineal | **error** (tangente inicial) | **error** | **error** | **aviso** y solución |
| No lineal, corrotacional | **aviso**; el análisis sigue | **aviso**; sigue | **aviso**; sigue | **aviso**; sigue |

En negrita lo que cambia (antes: solución sin aviso). La suite completa pasa sin tocar ningún modelo de test existente. El único test afectado comprobaba el comportamiento anterior sobre un modelo que de verdad es un mecanismo interno (dos bloques unidos por un nudo, cargado): pasa de `SingularTangentError` tras agotar las bisecciones a `IllPosedSystemError` antes de iterar. El diagnóstico de tangente singular del corrector conserva su prueba con la armadura corrotacional cargada.

## Ampliación del 2026-10-20 — mecanismos internos por la energía del modo (decisión del usuario)

### Lo que se encontró

La revisión pendiente del criterio de pivote nulo (la dependencia de las unidades, en «Lo que queda sin cubrir») proponía medir cada pivote frente a la diagonal de su fila y su columna, `|U_ii|/√(|K_rr|·|K_cc|)`, que no cambia al escalar un tipo de grado de libertad. Al medirlo apareció un defecto mayor, que afecta también al criterio que había: **el pivote de un mecanismo real no es cero, es ruido de redondeo, y crece con el tamaño del modelo**. Dos mecanismos que la carga no activa, medidos con SuperLU (cociente frente a la diagonal):

| Mecanismo | DOF | Pivote mínimo |
|---|---|---|
| Dos placas Quad4 unidas por un nodo | 28 / 4 288 / 66 304 / 263 680 | 1.2e-15 / 6.9e-12 / 1.4e-12 / 2.4e-11 |
| Dos bloques Hex8 unidos por una arista | 126 / 7 590 / 54 180 | 1.3e-14 / 2.4e-13 / 7.8e-12 |

Con el umbral de entonces (`1e-13·max|K|`) SuperLU dejaba de verlos desde unos 4 000 DOF, y Pardiso, el backend por omisión, desde unos 30 000 (28 560 DOF en el Hex8, 66 304 en el Quad4: su conteo de pivotes perturbados, `iparm(14)`, era 0). El mecanismo pasaba en silencio. En el otro lado, los modelos bien planteados con contrastes de rigidez grandes también dan pivotes pequeños: el mínimo de las 1 658 factorizaciones que la red de seguridad inspecciona en la suite es 3.5e-8 (el pórtico con zonas rígidas ×1e6), un material casi incompresible da ≈ 1 − 2ν (2e-7 con ν = 0.4999999) y unas zonas rígidas ×c dan ~1/c (3.5e-10 con SuperLU y 1.8e-11 con una LDLᵀ tipo CHOLMOD con c = 1e8). **Las dos poblaciones se solapan hacia 1e-11: ningún umbral de pivote las separa.** Un umbral de 1e-9 detectaría todos los mecanismos medidos a cambio de rechazar contrastes de 1e7–1e8 (MSC Nastran, con `PARAM,MAXRATIO = 1e7` por omisión, es aún más estricto); se descartó a favor de lo que sigue.

### Decisión

Un pivote pequeño **no decide**: sólo dispara la búsqueda del modo, y decide la **energía de deformación de ese modo en cada elemento**, que es la definición de mecanismo (un movimiento con `B·v = 0` en todos los elementos).

1. **Filtro** (`n_small_pivots`, constante `SUSPECT_PIVOT_RTOL = 1e-8`): los tres backends directos exponen sus pivotes pequeños frente a la diagonal. SuperLU con la diagonal de `U` y sus dos permutaciones (`LU = Pr·K·Pc`: el pivote `i` está en la fila `argsort(perm_r)[i]` y la columna `argsort(perm_c)[i]`); Cholesky con los pivotes de `LDLᵀ` y su permutación (`Factor.P()` en scikit-sparse 0.4, `CholeskyFactor.perm` en la 0.5); Pardiso con `pardiso_getdiag`, habilitado al crear su handle con `pardisoinit` (los mismos parámetros que la MKL usa por omisión) más `iparm(56) = 1`. Pardiso devuelve los pivotes de su matriz escalada y emparejada (MC64), de diagonal unitaria, así que ya están relativos a la diagonal. Medido sobre 90 000 DOF: factorizar 1.60 s sin la opción y 1.63 s con ella, leer la diagonal 0.2 ms, la misma solución a 3e-14. Sin pivotes pequeños no se hace nada más, y un modelo normal no paga nada.
2. **Modo de menor energía** del problema escalado `K·v = λ·diag(K)·v` (λ no depende de las unidades): iteración inversa en bloque con la factorización existente, `p = n_small_pivots + 2` vectores, tres iteraciones, D-ortonormalización y Rayleigh-Ritz. La separación entre el autovalor de un mecanismo (~1e-16) y el siguiente (≥ 1e-3 en los modelos medidos) hace que converja en una o dos.
3. **Veredicto por elementos** (`MECHANISM_ENERGY_RTOL = 1e-12`, valor fijado por el usuario): un modo con `λ < 1e3·MECHANISM_ENERGY_RTOL` es un mecanismo si todos los elementos que mueve tienen energía relativa `r_e = v_eᵀ·K_e·v_e / v_eᵀ·D_e·v_e < 1e-12`, con `D_e` la diagonal ensamblada en los DOF del elemento (invariante a las unidades, y no nula para una barra que se mueve en transversal, cuya `K_e` sí tiene diagonal nula en esa dirección). Medido:

   | Caso | Pivote mínimo / diag | λ mínimo | r_e máximo |
   |---|---|---|---|
   | Mecanismos, de 17 a 263 680 DOF | 3e-16 a 2.4e-11 | ≤ 1.4e-16 | **≤ 8.4e-17** |
   | Zonas rígidas, c = 1e6 / 1e8 / 1e10 | 3.5e-8 / 3.5e-10 / 3.5e-12 | 3.1e-10 / 3.1e-12 / 3.1e-14 | **6.8e-2** en los tres |
   | ν = 0.4999999 / 0.49999999999 | 2.1e-7 / 3.0e-11 | 1.4e-9 / 1.4e-13 | **1e-7 / 1e-11** |
   | Marco 3D ordinario | 2.8e-2 | 1.4e-2 | 1e-2 |

   La energía por elemento separa las dos poblaciones por al menos cuatro órdenes (quince en las zonas rígidas) y no depende del tamaño del modelo: el redondeo de `v_eᵀ·K_e·v_e` es local al elemento, y el error del modo entra al cuadrado.
4. **Localización.** Los modos que mueven los mismos elementos se agrupan en un mecanismo de `dimension` movimientos independientes (`InternalMechanism`), y se describen en términos del modelo: los elementos que se mueven y, si lo hacen como un sólido rígido (ajuste por mínimos cuadrados de `u = t + w × (x − c)`), el movimiento en una base legible: primero las traslaciones puras, después los giros sin la traslación que ya está libre, y tres giros con un punto fijo común como «giro alrededor de cualquier eje que pasa por P (como una rótula esférica)». Un modo que sólo mueve un campo escalar es «una parte no unida a ningún valor prescrito del campo 'T'». Ejemplos reales: *«la parte formada por los elementos 1025–2048 puede moverse sin deformarse: giro alrededor del eje z que pasa por (1, 1)»*; *«… de 2 formas independientes: traslación en la dirección y; giro alrededor del eje z que pasa por (2, 0.5)»*. Era la extensión anotada en «Lo que queda sin cubrir» (nodo y DOF concretos), resuelta por el modo en vez de por la permutación del pivote.
5. **Dónde actúa.** `LinearSolver` llama a `ensure_no_internal_mechanisms` una vez por factorización, antes de resolver (antes la capa 3 se miraba tras cada `solve`); la tangente inicial de los no lineales, en el inspector del corrector, que ahora recibe también la matriz factorizada (`inspector(factor, A)`). La política de la tangente inicial no cambia (error si nada puede ganar rigidez, aviso si sí). `IllPosedSystemError` gana `mechanisms` (las descripciones) y `zero_pivots` pasa a contar los movimientos independientes. Si el modo no pudiera calcularse (solución no finita), se rechaza como antes, por los pivotes.
6. **Coste.** Filtro: nulo con Pardiso, una copia de la diagonal de `U` con SuperLU (ya se hacía). Comprobación, sólo cuando el filtro salta: `3·p` resoluciones (entre el 1 % y el 6 % de la factorización con SuperLU, hasta el ~50 % con Pardiso, que factoriza muy deprisa) y una pasada por los elementos con `compute_global_stiffness` (1–4 s a 28 000–90 000 DOF en Python; sólo cuando hay un modo candidato). Los modelos bien planteados de la suite casi nunca la disparan; los extremos (zonas rígidas, ν muy cerca de 1/2) la pagan y pasan.

### Resultado

`tests/test_internal_mechanisms.py`: mecanismos que el criterio anterior dejaba pasar (Quad4 de 4 288 DOF, Hex8 de 4 104) rechazados y descritos igual con SuperLU y Pardiso; rótula esférica, traslación más giro, parte térmica sin valor prescrito y tangente inicial; zonas rígidas ×1e6, 1e8 y 1e10 en m y en mm, y ν = 0.4999999, sin mecanismo; el filtro salta en m y en mm; sin pivotes pequeños no se busca el modo; Pardiso detecta, por la diagonal de su factor, un mecanismo de 28 560 DOF que su `iparm(14)` no veía. La suite completa pasa sin tocar ningún modelo de test (2 709 con Pardiso, 2 701 sólo con SuperLU); los resultados de los modelos bien planteados no cambian más allá del redondeo (≤ 4e-15, por la inicialización explícita del handle de Pardiso). El capítulo 11 del manual de ejemplos cuenta la comprobación nueva.

## Lo que queda sin cubrir

- **Iterativo con un mecanismo interno sin carga que lo active**: sigue sin detectarse (un método de Krylov no tiene pivotes que inspeccionar), pero desde el 2026-09-30 se avisa una vez por análisis. El iterativo solo se activa a petición, así que el usuario no experto no llega a él por defecto.
- **Modelos que pueden ganar rigidez al deformarse** (corrotacionales, materiales unilaterales): un pivote nulo en la tangente inicial es sólo un aviso, porque la rigidez geométrica puede estabilizarlo; si no lo hace, el Newton fracasa después como tangente singular.
- **Mecanismos que aparecen durante el análisis** (un material que se ablanda hasta rigidez nula, una rótula plástica que completa un mecanismo): no son errores de planteamiento sino respuesta física (colapso), y no se rechazan; el Newton los reporta como `LoadExceedsCapacityError` o `SingularTangentError`.
- ~~**CHOLMOD real**~~: cubierto el 2026-09-30. La prueba real corre en el job `test-cholmod` del CI de Linux; en su primera ejecución (scikit-sparse 0.5.0) pasaron la verificación de activación y la suite completa por CHOLMOD.
- **Mal condicionamiento sin singularidad** (κ de 10¹⁰ a 10¹³): un directo es estable hacia atrás, así que el residuo sale pequeño aunque el error de la solución pueda ser grande. No hay estimador de condición. Es la advertencia que ANSYS y Abaqus hacen con penalizaciones y rigideces desproporcionadas.
- **Cambio de backend por tamaño**: deliberadamente no automático (ver §Memoria).
- ~~**Nodo y DOF concretos** de un mecanismo interno~~: resuelto el 2026-10-20 por el modo de menor energía, que dice qué elementos se mueven y, si lo hacen como un sólido rígido, cómo. Un mecanismo que no es rígido (articulaciones de más, como el panel sin diagonal de una armadura) se describe nombrando todos los elementos que mueve; agruparlos en subpartes rígidas lo haría más preciso.
- ~~**El criterio de pivote nulo de SuperLU y Cholesky depende de las unidades**~~: resuelto el 2026-10-20; el filtro mide el pivote frente a la diagonal y el veredicto es la energía del modo, que no depende de las unidades. Quedan dos restos menores, porque sólo afectan al filtro: el pivoteo parcial de SuperLU puede elegir otra secuencia de pivotes al cambiar de unidades (el pórtico con zonas rígidas ×1e8 da 6 pivotes pequeños en m y 2 en mm; el filtro salta en los dos), y el escalado MC64 de Pardiso no es exactamente invariante (×10 en el pivote mínimo de ese pórtico entre m y mm). Texto original (revisión adversaria del 2026-09-30): `ZERO_PIVOT_RTOL·max|K|` se aplica a la matriz sin escalar: es invariante a un escalado uniforme (`E × 10³` no cambia nada, medido), pero no a un cambio de unidades que escale distinto traslaciones y giros. En un pórtico, pasar de m a mm baja el menor pivote relativo unas 10⁴ veces (3,1·10⁻⁵ frente a 1,9·10⁻⁹). Pardiso no lo sufre: con `mtype = 11` y los `iparm` por omisión, la MKL escala y empareja la matriz (MC64) y compara con `10⁻¹³·‖A₂‖∞` de la matriz ya escalada. Consecuencia medida: un pórtico con zonas rígidas de 0,3 m (`A·c`, `I·c`) en mm da con `c = 10⁷` dos pivotes nulos con SuperLU, así que el no lineal lanza `IllPosedSystemError` y `LinearSolver` rechaza el modelo, y ninguno con Pardiso. Una `LDLᵀ` sin pivoteo, como la de CHOLMOD (emulada con SuperLU en modo simétrico, no con CHOLMOD real), cuenta ya dos con `c = 10⁶`, y siete con `c = 10⁷`. Escalada por su diagonal, esa matriz tiene un pivote mínimo relativo de 3,5·10⁻⁹. En unidades SI el criterio de pivotes no rechaza ese pórtico ni con `c = 10⁸`. Los mecanismos internos medidos (sólidos 2D y 3D unidos por un nudo o una arista, rótula por MPC, placa atada por una sola MPC) dan entre 10⁻¹⁶ y 10⁻¹⁴ con cualquiera de los dos criterios. Propuesta: medir el pivote respecto de la diagonal de su fila y su columna, `|U_ii|/√(|K_rr|·|K_cc|)` (`d_i/K_pp` en Cholesky), que es el cociente que reportan ANSYS y Abaqus.
- **Falso positivo de la capa 2 con zonas rígidas** (preexistente, medido el 2026-09-30). `‖F − K·u‖ ≤ 10⁻⁸·‖F‖` rechaza con SuperLU y con Pardiso el pórtico anterior desde `c = 10⁶` en m (desde `c = 10⁵` en mm), aunque su matriz escalada tiene un condicionamiento de ~10⁷–10⁸ y la solución es exacta a efectos prácticos. Un solver estable hacia atrás garantiza un residuo del orden de `ε·‖K‖·‖u‖`, no de `ε·‖F‖`, y con rigideces muy desproporcionadas el primero es mucho mayor. Propuesta: el error hacia atrás normado, `‖F − K·u‖/(‖K‖·‖u‖ + ‖F‖)`. Decisión pendiente del usuario. *Medido además el 2026-10-20*: un material casi incompresible, ν = 0.4999999, también se rechaza (residuo relativo 3e-8 frente a 1e-8), aunque la capa 3 lo deja pasar.

## Alternativas descartadas

- **Rechazar pivotes nulos o residuos grandes también dentro del Newton**: rompería el arc-length cerca de los puntos límite, donde resolver con una tangente casi singular es el algoritmo.
- **Sólo la capa 2 (equilibrio)**: no ve los mecanismos consistentes. Medido: residuo de 10⁻¹⁴ con un mecanismo con cargas equilibradas.
- **Sólo pivotes**: no dice al usuario *qué* está mal. La capa 1 describe el movimiento libre, que es lo que un usuario no experto puede corregir.
- **Estimador de condición** (`onenormest` sobre la factorización): cuesta varias resoluciones adicionales y el umbral de "demasiado mal condicionado" depende del problema. Se deja para cuando haya un caso real.
- **Pasar automáticamente al iterativo ante falta de memoria**: ver §Memoria.
- **Factorizar aparte la tangente inicial** (2026-09-30): una comprobación autónoma, al estilo de la capa 1, que ensamblara y factorizara `K(0)` antes de iterar. Duplica la operación que se lleva el 97 % del tiempo de un análisis, cuando el solver ya factoriza esa misma matriz en su primera iteración; el inspector del corrector la reutiliza sin coste de factorización.
- **Rechazar también los modelos corrotacionales con la tangente inicial singular** (2026-09-30): convertiría en error una red de cables sin pretensado, que la rigidez geométrica estabiliza al deformarse. Decisión del usuario: aviso.
- **Un umbral de pivote más holgado, 1e-9 frente a la diagonal** (2026-10-20): detectaría todos los mecanismos medidos, pero el pivote no separa las dos poblaciones (el de un mecanismo crece con el tamaño del modelo y se solapa con el de un contraste de rigidez legítimo hacia 1e-11), así que rechazaría zonas rígidas de 1e7–1e8 y aun así dejaría un margen de sólo 40 veces sobre el mecanismo más grande medido. Decisión del usuario: decidir por la energía del modo.
- **Estimar el menor autovalor con el backend iterativo** (2026-09-30): Lanczos o LOBPCG sobre `K` para distinguir un mecanismo. La frontera entre un autovalor nulo y uno de 10⁻¹⁴·‖K‖ no se resuelve de forma fiable con un método de Krylov, y un detector poco fiable con apariencia de fiable es peor que un aviso honesto.
