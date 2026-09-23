# ADR 0014 — Ensamblaje por lotes: familias derivadas de los contratos, estado por arreglos y un único kernel compilado

- **Estado**: aceptado
- **Fecha**: 2026-09-22
- **Alcance**: `Assembler` (`solidum/math/assembly.py`), nuevo paquete `solidum/math/batch/`, contratos opcionales en `Element` y `Material`/`ThermalMaterial`, kernels puntuales de los diez sólidos isoparamétricos, de los dos elementos térmicos y de los once materiales de sólidos y térmico. **No toca** el contrato obligatorio de elementos y materiales (`compute_element_state`, `compute_state`) ni la formulación de ningún componente. El cierre del mismo día (§9) añade, como consecuencias colaterales, el post-proceso por familia en el exportador VTK, el parámetro opt-in `record_internal_forces` de los solvers transitorios (y el campo `TransientResult.element_forces_history`) y el aislamiento del estado trial en la rigidez inicial. Documento de diseño previo: [`docs/propuestas/2026-09-22-vectorizacion-por-lotes.pdf`](../propuestas/2026-09-22-vectorizacion-por-lotes.pdf). Cierra las deudas técnicas #19, #15 y la parte de estado de #18 de [`STATUS.md`](../STATUS.md).

## Contexto

El ensamblaje recorría los elementos uno a uno en Python y, dentro de cada uno, los puntos de Gauss uno a uno. Un perfil sobre 10 000 Quad4 con J2 (propuesta, §1) mostró que más del 90 % del tiempo era sobrecarga de despacho —llamadas a funciones, arreglos pequeños, diccionarios de estado— y no aritmética: 120–170 µs por elemento y por ensamblaje, más 0,4 s por commit copiando en profundidad 40 000 diccionarios. El estado interno como diccionarios ocupaba 462 B por punto de Gauss, más que la propia matriz global.

Los kernels `@njit` existentes (cinemática, return mapping) aceleraban lo que ocurría dentro de ellos, no el hecho de invocarlos 80 000 veces desde Python. La solución conocida desde los supercomputadores vectoriales de los años ochenta (Hughes, Ferencz y Hallquist 1987; bloques `nblock` de la interfaz VUMAT de Abaqus; `MatrixFree` de deal.II; libCEED) es invertir el bucle: agrupar los elementos que comparten receta y evaluar todos sus puntos de Gauss dentro de un único bucle compilado.

La restricción del proyecto (`Reglas.md` §1, memoria `feedback_no_particularizar_al_presente`) es que la solución no fije hipótesis que el elemento, material o cinemática N+1 tengan que deshacer: nada de listas de tipos en el ensamblador, nada de cachés que presupongan geometría fija, nada que obligue a escribir dos veces la física.

## Decisión

### 1. La familia se deriva de los contratos, no de listas

Una *familia* es el conjunto de elementos del dominio con la misma clase de elemento, la misma **instancia** de material, la misma regla de cuadratura y el mismo número de nodos (`solidum/math/batch/family.py::family_key`). Todo se lee de atributos que ya existen. Un elemento nuevo cae en su propia familia sin tocar el ensamblador; dos instancias del mismo material son dos familias porque sus parámetros entran al kernel como constantes.

### 2. El camino por elemento sigue siendo el contrato obligatorio; el camino por lotes es una capacidad declarada

`compute_element_state` y `compute_state` no cambian y siguen siendo la referencia física. El camino por lotes lo declaran los componentes con atributos opcionales:

| Contrato | Quién | Qué |
|---|---|---|
| `BATCH_KINEMATICS` | elemento | función `@njit` con firma `KIN_SIG`: `detJ = kin(pt, coords, B)`, rellena `B` en un punto natural a partir de las coordenadas de referencia |
| `BATCH_SCALE` | elemento | nombre del atributo escalar que multiplica `dV` (`"thickness"` en los planos) |
| `STATE_SCHEMA` | material | `{nombre: forma}` de las variables internas; **obligatorio** en todo el catálogo (`{}` sin historia) |
| `BATCH_KERNEL` / `batch_kernel()` | material | función `@njit` con firma `MAT_SIG`: `C_t = mat(strain, S_old, S_new, params, C, sigma, C_out, flag)` |
| `batch_params()`, `batch_matrix()`, `initial_state()`, `batch_report()` | material | constantes de la familia, estado inicial por punto de Gauss y reporte de incidencias |

Un elemento sin `BATCH_KINEMATICS`, o cuyo material no declara esquema y kernel, sigue el camino por elemento dentro del mismo `Assembler` (dominios mixtos). Hoy quedan ahí los diez elementos estructurales 1D (fase 4 de la propuesta, diferida: pocos puntos de Gauss, ganancia marginal) y `CST_Embedded2D` (su `DiscontinuityState` y su Newton local no encajan en filas regulares).

### 3. Un único kernel de familia, compilado una vez y cacheado en disco

La cinemática y la constitutiva entran al kernel de familia como **funciones de primera clase tipadas por su firma** (`numba.types.FunctionType`), no como despachadores. Un despachador tiene identidad propia por proceso: pasarlo como argumento obligaría a compilar el kernel para cada pareja elemento × material en cada ejecución (medido: 3 s por pareja, sin caché). Una función tipada por su firma es un tipo estructural, así que `solid_family_kernel` (`solidum/math/batch/kernels.py`) se compila una sola vez para todas las parejas y el caché de Numba lo reutiliza entre procesos (0,14 s de carga). Los kernels puntuales se compilan perezosamente la primera vez que una familia los usa y también se cachean.

El kernel de familia no contiene física: recorre elementos y puntos de Gauss, llama a `kin` y `mat`, acumula `K_e` y `F_e` con bucles explícitos sobre el índice contiguo (para matrices de `3×8` a `6×81` la llamada a BLAS cuesta más que la operación) y vuelca los resultados por trozos. Existe en dos variantes con la misma firma y la misma aritmética, serie y paralela (§9).

### 4. Los dos caminos ejecutan las mismas funciones compiladas

Para que la equivalencia sea demostrable y no supuesta, las cinemáticas del camino por elemento se reescribieron como envoltorios de los mismos núcleos compilados que usa el camino por lotes (`_kin2d_core`, `_kin3d_core`, `_grad2d_core`, `_grad3d_core`; las envolturas `_*_from_dN` añaden el `ValueError` del camino por elemento, ver §9), y los materiales con historia exponen adaptadores que llaman a sus return mapping existentes. Los return mapping J2 (2D plane strain y 3D) se reescribieron sin asignaciones de memoria (`_j2_plane_strain_core`, `_j2_3d_core`) **con el mismo orden de operaciones** que la versión con arreglos: verificado bit a bit sobre 40 000 estados aleatorios en ramas elástica y plástica. El daño isótropo 2D/3D movió su física a un núcleo `@njit` que consumen ambos caminos.

Lo único que difiere entre caminos es el orden de suma en `B·u` (BLAS frente a bucle explícito) y en `Bᵀ·C·B`: `K` y `F_int` coinciden a `1e-14` relativo y el estado interno a `1e-13`, sobre las 44 combinaciones registradas elemento × material (`tests/test_batch_assembly.py`).

### 5. El estado interno se declara y se almacena por arreglos

`FamilyState` guarda variables internas y esfuerzos, committed y trial, como arreglos `(N·n_gp, n_state)` y `(N·n_gp, n_sigma)`; el commit es una copia de arreglo (no un intercambio de referencias: tras el commit, `vars_trial` y `vars` coinciden, como en `ElementState`). Cada elemento de la familia recibe un `BatchedElementState`, subclase de `ElementState` cuyas listas `vars`, `vars_trial`, `stresses`, `stresses_trial` son vistas que construyen el diccionario bajo demanda a partir de la fila y del esquema, y que escriben la fila al asignar. `compute_gauss_state`, el exportador VTK, `commit_state` por elemento y los tests que leen el estado siguen funcionando sin cambios. `Assembler.invalidate()` devuelve `ElementState` clásicos con el contenido actual, sin perder historia. Los elementos térmicos (sin estado) no reciben vista.

Memoria: 80 B por punto de Gauss en J2 2D (dos filas de 5 dobles, dos esfuerzos de 3 dobles) frente a 462 B con diccionarios, y sin duplicación transitoria en el commit.

### 6. Nada derivado de `x = X + u` se cachea

El kernel recibe siempre coordenadas de referencia y desplazamientos y recomputa jacobiano, inversa y `B` en cada evaluación (principio 4 de la propuesta). Cachear `B` por punto de Gauss duplicaría la memoria por elemento (78 KB por Hex20) y fijaría una hipótesis de geometría fija que una formulación lagrangiana futura no cumple. Recomputarla cuesta menos que el return mapping.

### 7. Trozos con presupuesto de memoria y mapa COO → CSR

Las `K_e` se escriben directamente en el bloque de la familia del vector COO (vista `(n, n_dof, n_dof)` sobre `data`, §9); el único temporal del trozo es `F_e`, acotado con `BATCH_MEMORY_BUDGET_BYTES` (64 MB por defecto; `Assembler(..., batch_memory_budget=…)`): `n_c = presupuesto / (8·n_dof)`. En la práctica un trozo cubre la familia entera y el límite sólo actúa en mallas enormes. El pico no crece con la malla.

El vector COO se ordena por familia (bloques contiguos) y se reutiliza entre ensamblajes, y la conversión a CSR —que `scipy` rehacía en cada ensamblaje ordenando y sumando duplicados, 115 ms para 8 000 Hex8— se sustituye por un mapa inverso COO → CSR calculado una vez (fuentes COO de cada posición CSR, en orden estable) y una reducción compilada, serie o paralela (§9), con el mismo orden de suma que `np.bincount`. Las filas y columnas COO se descartan tras construir el mapa. Aplica también al camino por elemento y a la matriz de masa.

### 8. Activación

`Assembler(domain, batch=None)` lee `BATCH_ASSEMBLY_DEFAULT = True`; `batch=False` fuerza el camino por elemento en todo el modelo (diagnóstico, comparaciones). `solidum.run` y el parser YAML no exponen la opción: el camino por lotes es equivalente por contrato y no hay razón física para elegir.

### 9. Cierre (mismo día): kernel paralelo, post-proceso por familia y señalización sin excepciones

Tres consecuencias del diseño anterior, implementadas en la misma sesión una vez medido el kernel serie:

**Kernel paralelo.** `solid_family_kernel_parallel` reparte los elementos de cada trozo en bloques de `PARALLEL_BLOCK = 32` con `prange`; cada bloque asigna su espacio de trabajo y cada elemento escribe sólo sus filas, así que no hay estado compartido y el orden de las sumas dentro de un elemento es el del kernel serie: el resultado es **bit a bit el mismo** con cualquier número de hilos (verificado sobre el barrido completo). La variante con espacios de trabajo por hilo (`numba.get_thread_id`) es algo más rápida pero **no se cachea en disco** (usa punteros dinámicos), lo que rompería el arranque; se descartó. `Assembler(parallel=None)` lee `BATCH_PARALLEL_DEFAULT = True`; el número de hilos lo gobierna Numba (`NUMBA_NUM_THREADS`, `numba.set_num_threads`). La reducción COO → CSR pasa de `np.bincount` a un kernel compilado sobre el mapa inverso cacheado (fuentes COO de cada posición CSR, en orden estable: misma suma que `bincount`), serie o paralelo según el ensamblador, y las `K_e` se escriben directamente en el bloque de la familia del vector COO (una vista, sin copia ni buffer aparte).

**Señalización sin excepciones.** Una excepción lanzada dentro de una región `prange` de Numba **se pierde**: la iteración aborta en silencio (medido: el elemento queda con `K = 0` y el bucle continúa). Por eso el contrato `KIN_SIG` cambia: la cinemática por lotes **no lanza**, devuelve `det J ≤ 0` cuando el jacobiano degenera (mismo chequeo relativo de Hadamard) y no escribe `B`. Los núcleos se dividen en `_grad2d_core`/`_kin2d_core` y `_grad3d_core`/`_kin3d_core` (no lanzan; los consumen las `BATCH_KINEMATICS`) y las envolturas `_*_from_dN` que lanzan `ValueError` (camino por elemento, sin cambios). El kernel marca `flags[fila] = −1` y salta el punto; `Family.evaluate` lanza el `ValueError` con la clase y el **id del elemento**, cosa que el camino por elemento nunca dio. Mismo protocolo en las variantes serie y paralela, de modo que se comportan igual.

**Post-proceso por familia.** `Family.gauss_state(U)` evalúa `compute_gauss_state` de toda la familia en un kernel (`GAUSS_SIG`, serie o paralelo) desde el estado committed y sin tocar el trial; devuelve `strain`/`stress` (o `grad_T`/`flux` con `q = −k·∇T` en las térmicas) con el eje de elemento delante, más `points_global = N(ξ_g)·X_e` a partir de `BATCH_SHAPE_FUNCTIONS` (o de `_SHAPE_FN`/`_shape_functions` en las bases que ya las tienen). Cada familia deja en sus elementos la referencia `_batch_family` al adoptar sus estados, así que `solidum.math.batch.postprocess.gauss_states(domain, U)` y el exportador VTK —que sólo ven el dominio— evalúan por familia y caen a `compute_gauss_state` en los elementos sueltos. `Material.batch_out_of_plane_stress(sigma, S)` es la versión por arreglos de `out_of_plane_stress` (vectorizada en los cuatro materiales 2D con plane strain; por defecto punto a punto). Medido: exportación VTK 245 → 38 µs/elemento (Quad4) y 331 → 44 (Hex8); `gauss_states` ×13 y ×18.

**Aislamiento del trial.** Como corolario, `Assembler.assemble_system` (y `Element.compute_global_stiffness`) ya no dejan el trial evaluado en `u = 0`: la evaluación es auxiliar y el trial previo se restaura (deuda #18).

## Consecuencias

**Rendimiento** (`examples/benchmarks/bench_assembly.py`, Python 3.14, NumPy 2.3.3, Numba 0.65, Windows, 16 hilos, 2026-09-22):

| Malla | Por elemento | Por lotes (serie) | Por lotes (paralelo) | Aceleración serie / paralelo | Commit por elemento | Commit por lotes | max Δ K / Δ F_int |
|---|---:|---:|---:|---:|---:|---:|---:|
| Quad4 100×100 + VonMises2D (10 000) | 1 528 ms (153 µs/elem) | 56,6 ms (5,7 µs/elem) | 13,8 ms (1,4 µs/elem) | ×27 / ×110 | 414 ms | 0,4 ms | 3,9e-16 / 6,7e-16 |
| Hex8 20³ + VonMises3D (8 000) | 2 810 ms (351 µs/elem) | 213 ms (26,6 µs/elem) | 33,6 ms (4,2 µs/elem) | ×13 / ×84 | 628 ms | 0,8 ms | 1,0e-15 / 7,0e-16 |

Serie y paralelo son bit a bit iguales. Se cumplen los umbrales de aceptación de la propuesta (≥ 20 en Quad4, ≥ 10 en Hex8). Tras esta ganancia el ensamblaje deja de dominar: la factorización dispersa pasa a ser el término principal de cada iteración de Newton, lo que motiva el Cholesky del ADR 0003 (o un backend multihilo) en los modelos grandes.

**Coste del componente N+1**: un material nuevo declara `STATE_SCHEMA` (una línea, obligatoria) y, si quiere acelerarse, un adaptador de diez líneas que llama a su return mapping; un sólido sobre una base existente hereda el kernel de la base; un elemento con cinemática nueva escribe un `BATCH_KINEMATICS` con las mismas funciones `@njit` que su `compute_element_state`. Nada de esto es necesario para que el componente funcione.

**Semántica**: `elem.state` puede ser una vista (`BatchedElementState`); lee y escribe como siempre. Cambiar el material o el espesor de un elemento existente tras el primer ensamblaje exige `Assembler.invalidate()` (las constantes de la familia se leen al construirla); añadir elementos o renumerar se detecta solo.

**Primer arranque**: `import solidum` compila el kernel de familia la primera vez en una máquina (2–3 s) y lo cachea en `__pycache__`; cada pareja elemento × material compila sus kernels puntuales en su primer uso (0,5–1 s) y también los cachea. En CI cada trabajo paga esto una vez.

## Alternativas descartadas

- **Cachear `B` por punto de Gauss** — duplica memoria y fija geometría (§6).
- **Un kernel por elemento (`BATCH_KERNEL` en cada clase)**, como planteaba la propuesta: obligaba a escribir el bucle de integración en cada elemento y a recompilarlo por material. Las funciones tipadas lo hacen innecesario: un solo bucle, un solo binario.
- **Despachadores como argumentos** o **clausuras sobre despachadores**: no se cachean entre procesos (medido: 3 s y 2,8 s por pareja en cada ejecución).
- **`jitclass`, C++ con pybind11, procesos en paralelo**: descartados en la propuesta (§9) por romper los contratos declarativos, el modelo de colaboración o por multiplicar el despacho en vez de eliminarlo.

## Diferido

- Fase 4 de la propuesta (elementos 1D): la ganancia es marginal con uno o dos puntos de Gauss; se retoma si un benchmark de marcos grandes lo justifica.
- Adaptadores sin asignaciones para Drucker-Prager y daño (hoy asignan arreglos por punto de Gauss dentro del return mapping; correctos, más lentos que J2).
- En el exportador VTK queda por arreglos lo que era por elemento; el bucle sobre nodos (diccionarios de DOF) y la escritura de `meshio` son ahora el resto del tiempo (~35 µs por elemento en total).
- Backend algebraico multihilo (Pardiso, CHOLMOD): con el ensamblaje a 1–4 µs por elemento, la factorización es el término dominante en 3D a partir de ~10⁵ grados de libertad.
