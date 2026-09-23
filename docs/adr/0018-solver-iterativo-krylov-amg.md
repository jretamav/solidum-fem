# ADR 0018 — Solver algebraico iterativo: CG / MINRES con AMG y modos de cuerpo rígido

- **Estado**: aceptado. La política de tolerancia del §4 —la parte matemática del ADR; el resto es plumbing— **validada por el usuario el 2026-09-23**.
- **Fecha**: 2026-09-23
- **Alcance**: `solidum/math/linalg/iterative.py` y `nullspace.py` (nuevos); `StiffnessProperties.near_nullspace` y el despachador (ADR 0003); `Assembler.near_nullspace`; el corrector de Newton (ADR 0015) y los sitios que construyen propiedades de matriz (`LinearSolver`, familia Newmark, `ThetaMethodSolver`); validación YAML de `linear_algebra`; `pyproject.toml` (extra `iterative`); CI (job con backends opcionales); `examples/benchmarks/bench_linear_algebra.py`. **No cambia** el backend por defecto, ninguna formulación ni ningún contrato de elemento o material.

## Contexto

El ADR 0017 mostró que la factorización se lleva el 97 % del tiempo de un análisis y añadió Pardiso. El 2026-09-23 se evaluó un solver iterativo y se **difirió** (opción I del ROADMAP) con tres obstáculos: CG sólo sirve para matrices SPD; `pyamg`, el precondicionador que escala, no instalaba en la máquina del usuario (Python 3.14, Windows); y el Newton inexacto sería una decisión de formulación.

El usuario objetó que ANSYS y Abaqus ofrecen varios solvers algebraicos, no uno, y pidió el iterativo "perfectamente implementado". Revisados, los tres obstáculos no aguantaban:

1. **SPD**: se resuelve con una familia (CG para SPD, MINRES para simétricas indefinidas, solver directo para no simétricas), que es lo que hacen los programas comerciales.
2. **`pyamg`**: la máquina sí tenía Visual Studio 2022 con el toolset C++. pip fallaba porque la carpeta de `vswhere.exe` no está en `PATH`: `vcvarsall.bat` imprime un aviso en ANSI que desalinea la salida UTF-16 que setuptools decodifica, y el entorno del compilador llega vacío. Con esa carpeta en `PATH`, `pyamg` 5.3.0 compila e instala sobre Python 3.14.
3. **Newton inexacto**: la objeción estaba **sobredimensionada**. El Newton inexacto aparece cuando se *relaja* la tolerancia del sistema lineal para ahorrar iteraciones. Con una tolerancia fija y estricta (ANSYS usa 1e-8 por defecto en PCG) el iterativo devuelve la misma corrección que el directo y el método de Newton no cambia; sólo cambia cómo se resuelve cada sistema, igual que con Pardiso.

La distinción que el iterativo aporta no es de velocidad sino de **memoria**: un directo guarda los factores, cuyo relleno crece con la talla; un iterativo guarda la matriz, el precondicionador y unos pocos vectores.

## Decisión

### 1. Opción explícita, nunca automática

`linear_algebra: iterative` en el bloque `solver` del YAML (o `linear_algebra="iterative"` en Python), con sufijo opcional de precondicionador: `iterative:amg`, `iterative:jacobi`, `iterative:none`; sin sufijo, AMG si `pyamg` está instalado y ninguno si no (con aviso único). El despachador **nunca** lo elige solo: por debajo de ~10⁵ DOF el directo multihilo es más rápido (§Mediciones), y el directo es robusto ante cualquier condicionamiento. Es la misma división que ANSYS (`SPARSE` por defecto, `PCG` a petición) y Abaqus (*direct* por defecto, *iterative* a petición). El parser YAML valida el valor al leer (`is_valid_override`) y rechaza `iterative:amg` si `pyamg` no está instalado.

### 2. Un método por tipo de matriz

- **CG precondicionado** para simétricas definidas positivas. Implementado en el módulo y no con `scipy.sparse.linalg.cg`, para detectar la **curvatura negativa** (`pᵀKp ≤ 0`, o `rᵀz ≤ 0` en el precondicionador): es la prueba de que `K` no es definida positiva, la misma información que da Cholesky al fallar. Al detectarla, el backend pasa a MINRES para esa matriz en vez de devolver basura.
- **MINRES** para simétricas indefinidas (arc-length cerca de un punto límite, ablandamiento), precondicionado con Jacobi sobre `|diag K|`. MINRES exige un precondicionador definido positivo y un AMG construido sobre una matriz indefinida no lo garantiza.
- Si la diagonal tiene alguna entrada `≤ 0`, la matriz no puede ser SPD (`d_i = e_iᵀKe_i`): se va directamente a MINRES sin intentar construir AMG, que invierte la diagonal y fallaría.
- Si CG agota las iteraciones sin exhibir curvatura negativa, se intenta MINRES antes de rendirse: sobre una matriz indefinida CG puede estancarse sin llegar a mostrarla.
- Las matrices **no simétricas** (Drucker-Prager no asociado, cargas seguidoras) van al solver directo automático con un aviso único. CG y MINRES exigen simetría, y la teoría de la agregación suavizada es la de operadores SPD.

### 3. AMG con los modos de cuerpo rígido como casi-núcleo

El multimalla por agregación suavizada construye sus niveles gruesos a partir de los vectores que el operador casi anula. En un problema escalar basta el vector constante (el default de `pyamg`); en elasticidad el núcleo de `K` sin apoyos son los **modos de cuerpo rígido**: 3 en el plano, 6 en el espacio.

`rigid_body_modes(domain)` (`solidum/math/linalg/nullspace.py`) los deriva **sólo del nombre de los DOF**, sin información del tipo de elemento: una traslación por cada `ux`/`uy`/`uz` presente; la rotación infinitesimal `u = e_a × (x − c)` (con `r_a = 1` si existe el giro) alrededor de cada eje que deja rastro en el modelo; y un modo constante por cada DOF no mecánico (`T`). En un dominio mixto los bloques no se mezclan. Las coordenadas se centran en el centroide y las rotaciones se escalan por la dimensión característica, para que el resultado no dependa de las unidades.

Verificado: `K·B = 0` a precisión de máquina (≤ 2,2·10⁻¹⁶ relativo) en modelos sin apoyos de **Quad4, Hex8, Truss2D, Truss3D, Frame2DEuler, Frame2DTimoshenko, Frame3D, Quad4Thermal y un dominio mixto Quad4 + Quad4Thermal**. Son el núcleo exacto, no una aproximación.

El efecto es decisivo (`Hex8 40³`, 2·10⁵ DOF, medido antes de implementar):

| Precondicionador | Iteraciones | Tiempo (setup + CG) |
|---|---:|---:|
| ninguno | 538 | 12,4 s |
| AMG sin modos rígidos | 112 | 26,8 s |
| **AMG con modos rígidos** | **11** | **9,3 s** |

**Sin los modos, AMG es peor que no precondicionar.** Con ellos, el número de iteraciones apenas depende de la talla (10, 12, 15 y 14 iteraciones de 3 600 a 2·10⁵ DOF), que es la propiedad que define a un multimalla bien construido.

Dos ajustes de `pyamg`, ambos por medición:

- **Nivel grueso con `splu`**, no con la pseudoinversa densa por defecto. Con seis modos cada agregado aporta seis incógnitas gruesas y la jerarquía puede detenerse con miles de incógnitas en el último nivel (2 058 en `Hex8 20³`); su pseudoinversa es una SVD densa O(n³) que costaba ~3 s frente a ~0,1 s.
- **Construcción completa en `factorize`**: `pyamg` factoriza el nivel grueso en la primera aplicación; se fuerza con una aplicación de calentamiento para que todo el coste quede donde el contrato del ADR 0003 lo pone (el Newton modificado reutiliza la "factorización") y no aparezca escondido en la primera resolución. Detectado por una anomalía: 12 iteraciones a 26k DOF tardaban 3,4 s y 15 a 86k tardaban 1,6 s.

### 4. Criterio de parada — *validado por el usuario (2026-09-23)*

```
‖b − K·x‖  ≤  ITERATIVE_RTOL · ‖b‖,      ITERATIVE_RTOL = 1e-10
```

- **Sobre el residuo verdadero**, recalculado al terminar, no sobre el residuo recursivo de CG, que deriva del verdadero por redondeo. Si el recursivo converge y el verdadero no cumple, CG se reinicia desde el iterado con el residuo recalculado (*residual replacement*, hasta `ITERATIVE_MAX_RESTARTS = 3`).
- **Cinco órdenes por debajo** de la tolerancia del Newton (`CONVERGENCE_RTOL_FORCE = 1e-5`, ADR 0007). Dentro del Newton equivale a un término de forzado constante `η = 1e-10`; la teoría del Newton inexacto (Dembo, Eisenstat y Steihaug, 1982) da `‖R_{k+1}‖ ≲ η·‖R_k‖ + C·‖R_k‖²`, y con `η = 1e-10` el término cuadrático domina en todo el rango útil. **No es Newton inexacto**: la tolerancia no se relaja con el residuo. Verificado: un análisis J2 no lineal da el mismo campo `U` que con el directo a `9·10⁻¹⁷` relativo.
- **En un análisis lineal** no hay criterio exterior que corrija, y el error de la solución está acotado por `κ(K)·rtol`. La cota es pesimista: en las mallas medidas el error relativo frente al directo fue de 1·10⁻¹¹ a 6·10⁻¹¹. Pero en un problema muy mal condicionado (láminas, penalizaciones rígidas) la cota puede ser relevante; es la razón de fondo por la que los manuales de ANSYS y Abaqus desaconsejan sus iterativos en esos casos, y este ADR hereda la advertencia.
- `ITERATIVE_MAX_ITER = 10 000` por resolución. Agotarlo es un fallo explícito (§5).

La elección de `1e-10` frente al `1e-8` de ANSYS es deliberadamente conservadora: el coste de dos órdenes más de tolerancia con AMG son unas pocas iteraciones (la convergencia es geométrica), y a cambio el iterativo no introduce una fuente de error por encima del orden del directo en ningún caso práctico.

**Validado por el usuario el 2026-09-23**, con tres decisiones asociadas:

- **Se mantiene `1e-10`** y no el `1e-8` de ANSYS. Medido con AMG: pasar de 10⁻⁸ a 10⁻¹⁰ cuesta unas 3 iteraciones (9/12/11 → 12/15/14 a 26k/86k/202k DOF), y a cambio el iterativo no aporta un error propio por encima del orden del directo. En un código de investigación prima que el resultado no dependa del solver algebraico elegido.
- **No se expone en el YAML.** Pensando en el usuario no experto en numérica (ADR 0019), una perilla de tolerancia invita a relajarla para "que converja". Queda configurable desde Python (`IterativeSolver(rtol=...)`); se expondrá si un caso real lo pide (`Reglas.md` §1).
- **Sin forzado de Eisenstat-Walker** (Newton inexacto propiamente dicho): cambiaría la sucesión de iterados y hoy no hay caso que lo justifique.

**Caveat heredado, registrado como mejora**: `‖b‖` es una norma euclídea que mezcla unidades (fuerzas y momentos en marcos; fuerzas y flujos de calor en un dominio mixto), así que el criterio relativo depende de la elección de unidades. Es la misma limitación del criterio del Newton (ADR 0007) y debe resolverse en ambos a la vez, midiendo el residuo escalado por la diagonal de `K` (deuda #22 de `docs/STATUS.md`).

### 5. Nunca falla en silencio

Si no se alcanza la tolerancia, se lanza `IterativeNotConvergedError` con el método, el precondicionador, las iteraciones, el residuo relativo alcanzado, la tolerancia, el tamaño y las salidas (directo, o AMG si no se usaba). **Nunca** se devuelve una solución que no cumple la tolerancia.

Es `RuntimeError`, así que los solvers la tratan como cualquier fallo del sistema lineal: el Newton abandona el paso. Pero el corrector (ADR 0015) la captura **antes** que el `RuntimeError` genérico, que hasta ahora se reportaba siempre como "matriz singular". Guarda el mensaje en `CorrectorResult.linear_solver_note` y lo añade a la excepción tipada final. Así, un iterativo que no converge no se diagnostica como "la rigidez tangente perdió rango", que sería falso.

### 6. Transporte de los modos: un proveedor perezoso

El backend sólo ve la matriz reducida; los modos los conoce el ensamblador (coordenadas, nombres de DOF, restricciones). `StiffnessProperties` gana un campo `near_nullspace: Callable[[], ndarray] | None` (fuera de la igualdad) y `Assembler.near_nullspace()` devuelve los modos **restringidos a los DOF libres** en el orden de las columnas de `T`. Con restricciones afines basta la restricción, porque el valor en un esclavo lo reconstruye `T` a partir de sus maestros. El resultado se cachea con la huella de topología y restricciones.

Es perezoso: sólo AMG lo invoca, así que el resto de backends no pagan su cálculo. Lo pasan el corrector (y con él los cinco solvers iterativos), `LinearSolver`, los seis sitios de la familia Newmark y `ThetaMethodSolver`. El despachador construye con propiedades sólo los backends que lo declaran (`ACCEPTS_PROPERTIES = True`).

### 7. Dependencia y CI

`pip install solidum-fem[iterative]` instala `pyamg` (wheels para Python 3.10–3.12; en 3.14/Windows se compila, ver anexo de la capa algebraica). Sin él, el backend funciona sin AMG.

El CI gana un job, `test-optional-backends` (Linux, Python 3.12, `pip install -e .[dev,fast,iterative]`), que **comprueba que Pardiso y `pyamg` están activos antes de correr la suite**: sin esa comprobación, un fallo silencioso de instalación dejaría saltados los tests de los backends y el job saldría verde sin probar nada. El job original sigue sin dependencias opcionales, que es la garantía del camino base.

## Mediciones

`examples/benchmarks/bench_linear_algebra.py`, `Hex8 n³` empotrada en `x = 0` (sistema reducido, apoyos reales), Windows, 16 hilos:

| Malla | DOF | SuperLU | Pardiso | Iterativo (setup + CG) | Iter. | Δ relativa frente al directo |
|---|---:|---:|---:|---:|---:|---:|
| Hex8 10³ | 3 630 | 0,55 s | **0,08 s** | 0,12 + 0,04 = 0,15 s | 10 | 2,2·10⁻¹¹ |
| Hex8 20³ | 26 460 | 20,1 s | **0,70 s** | 1,01 + 0,41 = 1,43 s | 12 | 1,4·10⁻¹¹ |
| Hex8 30³ | 86 490 | — | 4,29 s | 2,94 + 1,61 = 4,55 s | 15 | 5,3·10⁻¹¹ |
| Hex8 40³ | 201 720 | — | 17,6 s | 6,70 + 3,48 = **10,2 s** | 14 | 5,7·10⁻¹¹ |

Memoria pico (Pardiso: `iparm` 15–17 de la MKL; iterativo: `K` más el pico de `factorize` + `solve` medido con `tracemalloc`):

| DOF | `K` | Pico Pardiso | Pico iterativo |
|---:|---:|---:|---:|
| 26 460 | 22 MB | 313 MB (×14) | ~92 MB (×4,2) |
| 86 490 | 75 MB | 1 461 MB (×20) | ~312 MB (×4,2) |
| 201 720 | 178 MB | 4 434 MB (×25) | **~740 MB (×4,2)** |

- **Tiempo**: el cruce está en torno a 10⁵ DOF; a 2·10⁵ el iterativo gana ×1,7. Por debajo, Pardiso es más rápido: el iterativo es para modelos grandes, no para acelerar los medianos.
- **Memoria**: el iterativo ocupa un múltiplo **constante** de `K` (×4,2: `K`, una jerarquía AMG de tamaño comparable —complejidad de operador 1,19— y temporales de construcción); el directo, un múltiplo **creciente** (×14 → ×25). A 2·10⁵ DOF la diferencia es ×6 y aumenta con la talla. *Corrige una estimación previa de la sesión, "×25 menos memoria", que comparaba con `K` sola y olvidaba la jerarquía AMG.*
- **En un Newton**, el setup de AMG se paga por iteración como la factorización del directo. Con el Newton modificado (`freeze_tangent_after_iter`) se reutiliza, igual que la factorización.

## Consecuencias

**Coste del componente N+1**: un precondicionador nuevo (ILU, AMG por bloques, un AMG de otra biblioteca) es una función que devuelve `apply_M`; un método nuevo (GMRES para no simétricas) es otra rama de `IterativeFactorized`. Ni el despachador ni los solvers cambian.

**Superficie**: `IterativeSolver`, `IterativeNotConvergedError` y `rigid_body_modes` se exportan desde `solidum.math.linalg`; `Assembler.near_nullspace()` es público; `linear_algebra` admite `iterative[:amg|jacobi|none]`.

**Lo que este ADR NO hace**:

- **No cambia el default.** El despachador sigue eligiendo Cholesky → Pardiso → LU.
- **No implementa Newton inexacto** (tolerancia ligada al residuo, Eisenstat–Walker). Ahorraría iteraciones de Krylov lejos de la solución, pero cambia la sucesión de iterados: es una decisión de formulación del usuario, no de implementación.
- **No resuelve iterativamente las no simétricas.** GMRES/BiCGSTAB con un precondicionador adecuado es la extensión natural si aparece un caso real grande con plasticidad no asociada.
- **No usa agregación por bloques** (los DOF de un nodo agregados juntos, formato BSR). Con apoyos, el sistema reducido pierde la estructura por bloques y la agregación por incógnita ya da 10-15 iteraciones; queda como mejora si un caso lo pide.
- **No da inercia**: `n_negative_pivots` es `None`. La deuda #7 sigue abierta.

## Alternativas descartadas

- **`scipy.sparse.linalg.cg`**: no detecta la curvatura negativa; sobre una matriz indefinida puede estancarse o devolver una solución sin avisar.
- **AMG sin modos rígidos**: medido, peor que no precondicionar (112 iteraciones frente a 11).
- **Jacobi como precondicionador por defecto**: medido, iguala a no precondicionar en elasticidad.
- **Pseudoinversa densa en el nivel grueso** (default de `pyamg`): O(n³) cuando la jerarquía se detiene con miles de incógnitas.
- **Iterativo en la regla automática**: pierde en tiempo por debajo de 10⁵ DOF, que es hoy la totalidad de los modelos del proyecto.
- **Entorno aparte con Python 3.12** para tener `pyamg` con wheels: innecesario una vez diagnosticado el fallo de pip con Visual Studio.
