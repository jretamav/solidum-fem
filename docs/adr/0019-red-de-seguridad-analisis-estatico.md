# ADR 0019 — Red de seguridad del análisis estático: mecanismos, equilibrio y pivotes nulos

- **Estado**: aceptado
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

**Cero falsos positivos**: la suite completa (1 479 pruebas con Pardiso y 1 471 con SuperLU, que ejercita su comprobación de pivotes) pasa sin cambios. Las tres capas se ejecutan en todos los análisis estáticos de la suite: marcos, armaduras, cables corrotacionales, discontinuidades embebidas, restricciones lineales, térmico estacionario, 3D.

**Coste** (Hex8 40³, 2·10⁵ DOF, frente a 17 s de factorizar y resolver con Pardiso): capa 1, 0,42 s (2,5 %); capas 2 y 3, 30 ms. La capa 1 bajó de 0,85 s al corregir dos cosas: una SVD completa que construía una matriz densa de (restricciones × restricciones), unos 200 MB con 5 000 restricciones, sustituida por QR económica + SVD del factor `R` de 6×6; y un recorrido por nodos repetido por cada nombre de DOF.

## Lo que queda sin cubrir

- **Iterativo con un mecanismo interno sin carga que lo active**: un método de Krylov no tiene pivotes que inspeccionar, y sobre un sistema singular consistente converge a una solución, en la práctica razonable pero de un modelo inestable. El iterativo solo se activa a petición, así que el usuario no experto no llega a él por defecto.
- **Backend Cholesky (extra opcional `cholmod`)**: `CholeskyFactorized` no expone `n_zero_pivots`, así que la capa 3 no actúa con él. Hallado el 2026-09-24 al verificar el ejemplo 11 del manual de ejemplos; sin comprobar con CHOLMOD real. Deuda #28 en [STATUS](../STATUS.md).
- **Mal condicionamiento sin singularidad** (κ de 10¹⁰ a 10¹³): un directo es estable hacia atrás, así que el residuo sale pequeño aunque el error de la solución pueda ser grande. No hay estimador de condición. Es la advertencia que ANSYS y Abaqus hacen con penalizaciones y rigideces desproporcionadas.
- **Cambio de backend por tamaño**: deliberadamente no automático (ver §Memoria).
- **Nodo y DOF concretos** de un mecanismo interno: ANSYS y Abaqus los reportan a partir del pivote. Aquí se dice qué pasa, no dónde, para los internos; para los rígidos sí se dice el movimiento. Extensión posible: mapear el pivote nulo a su DOF a través de la permutación de la factorización.

## Alternativas descartadas

- **Rechazar pivotes nulos o residuos grandes también dentro del Newton**: rompería el arc-length cerca de los puntos límite, donde resolver con una tangente casi singular es el algoritmo.
- **Sólo la capa 2 (equilibrio)**: no ve los mecanismos consistentes. Medido: residuo de 10⁻¹⁴ con un mecanismo con cargas equilibradas.
- **Sólo pivotes**: no dice al usuario *qué* está mal. La capa 1 describe el movimiento libre, que es lo que un usuario no experto puede corregir.
- **Estimador de condición** (`onenormest` sobre la factorización): cuesta varias resoluciones adicionales y el umbral de "demasiado mal condicionado" depende del problema. Se deja para cuando haya un caso real.
- **Pasar automáticamente al iterativo ante falta de memoria**: ver §Memoria.
