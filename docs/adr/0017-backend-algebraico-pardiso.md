# ADR 0017 — Backend algebraico multihilo: la factorización era el cuello de botella

- **Estado**: aceptado
- **Fecha**: 2026-09-22
- **Alcance**: `solidum/math/linalg/pardiso.py` (nuevo), el despachador del ADR 0003 (`dispatcher.py`), `solidum/math/linalg/__init__.py`, `pyproject.toml` (extras `fast` y `cholmod`) y `examples/benchmarks/bench_linear_algebra.py` (nuevo). **No cambia** ninguna formulación, ningún contrato de elemento o material, ni la API pública de los solvers. Los resultados son los mismos a precisión de máquina.

## Contexto

El ADR 0014 dejó el ensamblaje por lotes en 1,4 µs/elemento (Quad4+J2) y 4,2 µs (Hex8+J2), entre ×84 y ×110 más rápido que el camino por elemento. Cerrado aquello, la pregunta era qué quedaba en la misma línea de optimización, y la respuesta honesta exigía **medir**, no razonar: una lista de mejoras escrita de memoria arrastraba tanto trabajo de rendimiento como de mantenibilidad, mezclados.

El perfilado de extremo a extremo de un análisis no lineal completo (Hex8 15³, 12 288 DOF, `VonMises3D`, 5 pasos de carga) dio el reparto real:

| Fase | Tiempo | Fracción |
|---|---|---|
| **Factorización (`gstrf`, SuperLU)** | **29,8 s** | **97 %** |
| Ensamblaje tangente (15 llamadas) | 0,33 s | 1 % |
| Reducción, residuo, commit, resto | ~0,5 s | 2 % |

**El ADR 0014 optimizó ×110 una fase que hoy representa el 1 % del tiempo.** Fue el trabajo correcto entonces —el ensamblaje sí dominaba— pero está amortizado: seguir puliéndolo (los diferidos de aquel ADR: elementos 1D por lotes, adaptadores sin asignaciones de Drucker-Prager y daño) produciría mejoras invisibles.

La causa de fondo del coste de la factorización es el **relleno** (*fill-in*). Medido sobre `Hex8 n³` con SuperLU:

| Malla | DOF | nnz(K) | nnz(L+U)/nnz(K) | `splu` |
|---|---:|---:|---:|---:|
| Hex8 10³ | 3 993 | 262 877 | ×10,0 | 0,47 s |
| Hex8 15³ | 12 288 | 867 600 | ×15,0 | 3,24 s |
| Hex8 20³ | 27 783 | 2 023 214 | ×23,5 | 18,32 s |
| Hex8 25³ | 52 728 | 3 920 073 | ×41,9 | 110,10 s |

El relleno crece con la talla y arrastra tiempo **y** memoria: los 3,9 M de no-ceros de `K` en la malla mayor se convierten en 164 M en `L+U`. El coste se multiplica por ~5,5 cada vez que el problema se duplica.

El propio ADR 0014 anotó esto como diferido, con una estimación que la medición corrige: *"la factorización es el término dominante en 3D a partir de ~10⁵ grados de libertad"*. **Domina desde ~4 000.**

### Lo que se descartó por medición, no por criterio

Antes de añadir una dependencia se probó el atajo: cambiar el reordenamiento de SuperLU vía `permc_spec`. Sobre `Hex8 20³`, las cuatro opciones de SciPy dan:

| `permc_spec` | nnz(L+U) | fill-in | tiempo |
|---|---:|---:|---:|
| **COLAMD** (default) | **50,3 M** | **×24,8** | **19,1 s** |
| MMD_AT_PLUS_A | 63,3 M | ×31,3 | 191,3 s |
| MMD_ATA | 65,9 M | ×32,6 | 34,2 s |
| NATURAL | 73,6 M | ×36,4 | 37,5 s |

**El default ya es el mejor de los disponibles.** No hay mejora gratis en configuración; SciPy no ofrece disección anidada (METIS) ni factorización multihilo. Hace falta una biblioteca externa.

## Decisión

### 1. `PardisoSolver` como backend opcional, hermano de `CholeskySolver`

Se añade un backend sobre **Intel MKL Pardiso** (vía `pypardiso`), con el mismo patrón que Cholesky: el módulo importa su dependencia y lanza `ImportError` si falta; `solidum.math.linalg.__init__` y el despachador lo absorben. Sin la dependencia, el sistema se comporta exactamente como antes.

Pardiso aporta las dos cosas que SciPy no tiene: **factorización multihilo** y **reordenamiento por disección anidada**, que es lo que ataca el relleno.

Medido (Windows, 16 hilos, `examples/benchmarks/bench_linear_algebra.py`):

| Malla | DOF | SuperLU | Pardiso | Aceleración | Δ relativa |
|---|---:|---:|---:|---:|---:|
| Hex8 10³ | 3 993 | 0,47 s | 0,10 s | ×4,7 | 6,4e-13 |
| Hex8 15³ | 12 288 | 3,24 s | 0,27 s | ×12,0 | 1,8e-13 |
| Hex8 20³ | 27 783 | 18,32 s | 0,78 s | ×23,5 | 1,1e-13 |
| Hex8 25³ | 52 728 | 110,10 s | 1,97 s | ×56,0 | 6,4e-14 |

Sobre el **análisis no lineal completo** (5 pasos, 10 factorizaciones), con el campo `U` final idéntico a `1e-16` relativo:

| Malla | LU | Pardiso | Aceleración |
|---|---:|---:|---:|
| Hex8 15³ | 32,3 s | 5,8 s | ×5,6 |
| Hex8 20³ | 212,0 s | 9,8 s | ×21,7 |

**La aceleración crece con el tamaño**, que es la propiedad que importa: el backend vale más cuanto más grande es el problema, justo donde hoy el proyecto se atasca.

### 2. Pardiso ocupa el lugar de **LU**, no el de Cholesky

En la regla automática del despachador, Pardiso entra donde caería `LUSolver`: mismo dominio de aplicación (real, sin hipótesis de simetría ni positividad) y misma solución. La preferencia por Cholesky en matrices SPD se mantiene cuando `scikit-sparse` está instalado, porque explotar la simetría sigue siendo la mejor opción disponible para ese caso.

El orden resultante: **Cholesky** (SPD, si está) → **Pardiso** (si está) → **LU** (siempre). El `override` explícito del usuario (`linear_algebra: lu`) se respeta como antes: sigue siendo herramienta de diagnóstico, no decisión de modelado.

Consecuencia deliberada: con Pardiso instalado, el aviso *"scikit-sparse no está instalado, degradando a LU"* deja de emitirse, porque ya no se degrada a LU sino a algo mejor.

### 3. Un handle de MKL por proceso, y una factorización obsoleta falla de forma ruidosa

Construir un `PyPardisoSolver` localiza la DLL de la MKL recorriendo el sistema de archivos con `glob`: **2,5 s por instancia**, medido. Con una instancia por factorización, un Newton de diez iteraciones gastaba **25 s buscando la biblioteca frente a 7 s factorizando** — el perfilado lo mostró como `glob._iterdir` con 750 000 llamadas. El handle se crea una sola vez por proceso y se comparte.

Eso introduce un riesgo que se blinda explícitamente: el handle guarda **una sola** factorización, así que factorizar de nuevo invalida la anterior. Si alguien mantiene dos vivas a la vez (p. ej. `M` y `A` en un solver dinámico) y resuelve con la vieja, obtendría **la solución del otro sistema, sin error**. `PardisoFactorized` lleva un testigo y lanza `RuntimeError` en ese caso, nombrando la salida (refactorizar, o usar `lu`). Un fallo ruidoso es preferible a un resultado silenciosamente equivocado.

### 4. El sistema vacío `0×0` se atiende en el backend

Un modelo de un elemento controlado por desplazamiento deja el sistema reducido sin incógnitas. SuperLU acepta el caso y devuelve un vector vacío; la MKL lo rechaza como singular. Lo destapó la suite: **13 tests de Drucker-Prager y daño** empezaron a fallar con `SingularTangentError` al activar el backend, y el diagnóstico —matriz `0×0`, no materiales rotos— exigió capturar la matriz que llegaba al backend. Se resuelve devolviendo una factorización vacía, para que Pardiso sea sustituto equivalente de LU también en el borde.

### 5. Sin conteo de inercia: la deuda #7 **no** queda cerrada

`pypardiso` usa por defecto `mtype = 11` (real no simétrica), que recibe la matriz completa en CSR y **no reporta inercia**. El tipo que sí la reporta, `mtype = -2` (real simétrica indefinida), espera **sólo el triángulo superior**: pasarle la matriz completa produce una violación de acceso en la MKL —no una excepción de Python— como se comprobó al evaluarlo.

Por tanto `PardisoFactorized.n_negative_pivots` es `None`, y el `sign-of-pivot tracking` exacto que pide la deuda #7 (LDLᵀ de Bunch-Kaufman para el caso cohesivo+embedded del `DissipationArcLengthSolver`) **sigue abierto**. El `ldlt.py` del ADR 0003 ya nombraba `pypardiso` como candidato para ello; este ADR confirma que la biblioteca sirve, pero que el trabajo es otro: construir el triángulo superior, validar la simetría y exponer `iparm[22]`. No se hace aquí porque el objetivo era el tiempo de factorización, no el diagnóstico de bifurcaciones.

### 6. Dependencia declarada como extra

`pip install solidum-fem[fast]` instala `pypardiso`, que trae wheels para Windows, Linux y macOS y no requiere compilador. Se añade además el extra `cholmod` para `scikit-sparse`, que en Windows sí exige Visual Studio (`conda install -c conda-forge scikit-sparse` es la vía práctica) — distinción que antes sólo vivía en el texto de un warning.

El CI **no** instala ninguno de los dos: es la garantía de que el camino sin dependencias opcionales sigue siendo correcto. La evidencia es el **CI en verde** (Linux × 3 Pythons, sin `pypardiso`) sobre el commit de este ADR.

> **Corrección (2026-09-23, ADR 0018).** La versión original de este párrafo afirmaba además una verificación local "con `pypardiso` bloqueado: 1 453 pasan, 9 skipped, idéntico a con él instalado". Esa verificación **no fue válida**: el bloqueador usaba `find_module`/`load_module`, que Python 3.12+ ya no consulta en `sys.meta_path`, así que no bloqueaba nada — de ahí que el recuento saliera idéntico cuando debían saltarse los 6 tests del backend. Con un bloqueador por `find_spec`, la suite sin `pypardiso` da el recuento esperado (los 6 tests de Pardiso saltados). La conclusión se sostenía por el CI; la afirmación local era falsa. Desde el ADR 0018, además, el CI tiene un job que instala los backends opcionales y comprueba que están activos antes de correr la suite.
>
> El benchmark `bench_linear_algebra.py` se rehízo en el ADR 0018: ahora mide sistemas **con apoyos reales** en vez de `K + 10⁶·I` regularizada, y añade el backend iterativo y la memoria. Las cifras de la tabla de arriba corresponden a la versión regularizada; con apoyos son del mismo orden (Hex8 20³: 0,70 s con Pardiso frente a 0,78 s).

## Consecuencias

**Rendimiento**: los análisis grandes dejan de estar limitados por la factorización monohilo. Hex8 25³ (52 728 DOF) pasa de 110 s a 2 s por factorización; el análisis no lineal completo de Hex8 20³, de 212 s a 9,8 s.

**Equivalencia**: la solución es la misma —`1e-13` relativo en el sistema lineal, `1e-16` en el campo `U` final de un análisis no lineal—. No cambia la sucesión de iterados del Newton ni el criterio de convergencia: es la misma matemática resuelta más rápido.

**Coste del backend N+1**: el ADR 0003 ya tenía la estructura correcta. Añadir Pardiso fueron un módulo nuevo, tres líneas en el despachador y dos en el `__init__`; ningún solver, elemento o material se enteró. Un backend futuro (MUMPS, CHOLMOD por METIS, un iterativo) entra igual.

**Lo que este ADR NO hace**, y conviene no confundir:

- **No cambia el método de Newton.** El Newton modificado (`freeze_tangent_after_iter`) sigue desactivado por defecto. Medido, `freeze=1` reduce Hex8 15³ de 35,2 s a 17,7 s, pero congelar la tangente cambia la convergencia de cuadrática a lineal, y cuánto cuesta eso depende del régimen (plasticidad severa, softening, rotaciones grandes). Esa es una decisión de formulación, no de implementación: queda para un estudio con barrido sobre regímenes contrastados.
- **No cierra la deuda #7** (ver §5).
- **No introduce un solver iterativo.** Un CG precondicionado evitaría la factorización por completo y con ella el límite de memoria que impone el relleno; es el salto cualitativo siguiente y merece componente propio con su spec.

## Alternativas descartadas

- **Otro `permc_spec` de SciPy.** Medido: COLAMD, el default, ya es el mejor de los cuatro (§Contexto). Descartado por datos.
- **`scikit-sparse` (CHOLMOD) como backend principal.** Es mejor que Pardiso para SPD —explota la simetría—, pero en Windows no tiene wheels y exige Visual Studio; el intento de instalación falló con `Unable to find a compatible Visual Studio installation`. Sigue soportado y con preferencia en la regla automática cuando está presente: no se sustituye, se complementa.
- **MUMPS (`python-mumps`).** Equivalente en prestaciones y con LDLᵀ real, pero de instalación difícil en Windows —el mismo problema que CHOLMOD—. Candidato si alguna vez se necesita la inercia (deuda #7) y Pardiso con `mtype = -2` resulta impracticable.
- **Instancia de `PyPardisoSolver` por factorización.** Es la forma obvia y cuesta 2,5 s de `glob` cada vez (§3). Descartada por medición.
- **Hacer Pardiso obligatorio.** Arrastraría la MKL (cientos de MB) a toda instalación, incluidos los modelos pequeños donde no aporta, y ataría el proyecto a un binario propietario. Como extra opcional, el camino base sigue siendo SciPy puro.
