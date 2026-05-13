# ADR 0003 — Capa algebraica `fenix.math.linalg` y despachador de solver para `K·x = y`

**Fecha**: 2026-04-30
**Estado**: Aceptado (2026-04-30)

## Resumen ejecutivo

Fenix introduce una capa nueva, `fenix.math.linalg`, que abstrae la resolución del sistema algebraico `K·x = y` detrás de una interfaz uniforme con varios backends (Cholesky, LDLᵀ, LU, iterativos, eigensolver). Un **despachador interno** los selecciona automáticamente a partir de propiedades declarativas de `K` (simetría, positividad, talla) — el usuario nunca elige solver algebraico desde el YAML como decisión obligatoria. Esta capa es el cimiento sobre el que crecerán los próximos tipos de análisis (modal, pandeo, dinámica implícita, térmico-mecánico) sin romper los solvers existentes ni la API pública.

## Contexto

Hoy los tres solvers de `fenix/math/solvers.py` (`LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`) llaman directamente a `scipy.sparse.linalg.spsolve` (SuperLU, factorización LU general). Eso tiene dos problemas:

**1. Coste algebraico subóptimo.** SuperLU no detecta ni explota simetría ni positividad. Para una `K` simétrica y positiva definida — el caso de la mayoría del catálogo actual — pagamos del orden de **2× tiempo y 2× memoria** respecto a un Cholesky sparse. El sobrecoste compone con la talla del problema.

**2. Cuello de botella arquitectural.** Los próximos tipos de análisis que el proyecto va a necesitar (ver hoja de ruta abajo) no se reducen a `K·x = y`: requieren factorización reusable, productos `K^{-1}·v` repetidos, eigensolvers, sistemas con shift `(K − σM)`. Si cada solver no lineal llama a `spsolve` directamente, añadir esas capacidades implica tocar todos los solvers en cada extensión. Contrario a Reglas.md §1.

**Estado del arte.** Todos los códigos FEM de referencia (ANSYS, Abaqus, NASTRAN, Code_Aster, OpenSees, FEAP) resuelven este dilema con la misma arquitectura: **un dispatcher central con backends algebraicos especializados**. Difieren en filosofía:

- **ANSYS / Abaqus**: el dispatcher es interno, default `auto`, el usuario solo elige "directo vs iterativo" en casos especiales.
- **OpenSees / FEAP**: el dispatcher es explícito, el usuario elige el backend desde el archivo de entrada (cultura investigación).
- **Backends industriales**: Pardiso (Intel MKL) y MUMPS dominan; cubren LU + LDLᵀ + Cholesky en una sola librería con auto-detección. CHOLMOD (open-source, solo SPD) es la alternativa ligera.

FEAP en particular es ilustrativo porque su `psolve.f` despacha sobre `stype` (entero: 0 diagonal, ≤−1 directo, ≥1 iterativo), incluye soporte mixto simétrico/no-simétrico (primeras `neqs` ecuaciones simétricas, resto no), y aprovecha la factorización LDLᵀ para reportar el número de pivotes negativos — diagnóstico gratuito de bifurcaciones (Sturm sequence).

**Filosofía elegida para Fenix.** Por el modelo de colaboración de Reglas.md §2 — el usuario no escribe código ni decide plumbing — Fenix adopta la línea ANSYS/Abaqus: dispatcher interno con auto-detección, default `auto`, override YAML solo como herramienta de diagnóstico. La elección del solver algebraico **no es decisión de modelado**, es consecuencia de qué materiales/elementos/análisis se hayan declarado.

## Decisión

### 1. Capa nueva `fenix.math.linalg`

Módulo nuevo con la siguiente estructura:

```
fenix/math/linalg/
    __init__.py
    base.py            # Protocol LinearAlgebraSolver, dataclass StiffnessProperties
    lu.py              # LUSolver (envuelve spsolve actual)
    cholesky.py        # CholeskySolver (scikit-sparse / CHOLMOD)
    ldlt.py            # LDLTSolver (diferido; placeholder)
    iterative.py       # CGSolver, GMRESSolver (diferido)
    eigen.py           # EigenSolver (ARPACK; diferido para modal/pandeo)
    dispatcher.py      # select_solver(props, override) -> LinearAlgebraSolver
    cache.py           # FactorizedSolver: factoriza una vez, resuelve N veces
```

Interfaz uniforme:

```python
class LinearAlgebraSolver(Protocol):
    def solve(self, K: sp.csr_matrix, b: np.ndarray) -> np.ndarray: ...
    def factorize(self, K: sp.csr_matrix) -> "FactorizedSolver": ...

class FactorizedSolver(Protocol):
    def solve(self, b: np.ndarray) -> np.ndarray: ...        # reutiliza factorización
    def shifted_solve(self, b: np.ndarray, sigma: float, M) -> np.ndarray: ...  # opcional
```

`factorize` es la pieza estratégica: separa la fase cara (factorización numérica) de la barata (sustitución). Lo necesitan Newton modificado, dinámica implícita con paso fijo, y los eigensolvers tipo Lanczos con shift-invert.

### 2. Propiedades declarativas de `K`

Por analogía con `STRAIN_DIM`, `DOF_NAMES`, etc., los flags algebraicos son **declarativos** y se agregan automáticamente:

```python
@dataclass(frozen=True)
class StiffnessProperties:
    is_symmetric: bool
    is_positive_definite: bool
    size: int
```

Origen de cada flag:
- `is_symmetric`: AND lógico sobre `material.is_symmetric` (default `True`; se sobrescribe en plasticidad no asociada, daño anisótropo con back-stress, etc.) y `element.preserves_symmetry` (default `True`; falso para follower loads y formulaciones corrotacionales con rotaciones finitas no simetrizadas).
- `is_positive_definite`: arranca en `True` para `LinearSolver` y `NonlinearSolver`; arranca en `False` para `ArcLengthSolver` (régimen postcrítico). Se degrada automáticamente a `False` si el backend reporta pérdida de positividad (Cholesky aborta, LDLᵀ encuentra pivote negativo).
- `size`: `domain.total_dofs`.

Ningún flag viene del YAML.

### 3. Despachador interno

```python
def select_solver(props: StiffnessProperties,
                  override: str | None = None) -> LinearAlgebraSolver:
    if override is not None:
        return _REGISTRY[override]()              # diagnóstico/benchmark

    if props.size > LARGE_PROBLEM_THRESHOLD:
        return CGSolver() if props.is_positive_definite else GMRESSolver()

    if props.is_symmetric and props.is_positive_definite:
        return CholeskySolver()
    if props.is_symmetric:
        return LDLTSolver()                       # cuando esté implementado; fallback LU mientras tanto
    return LUSolver()
```

Los tres solvers no lineales actuales pasan a llamar a `select_solver(...)` una vez por paso y reutilizan la instancia. **Cero cambio en su lógica de iteración o convergencia**: la migración es localizada en las cuatro líneas que hoy hacen `spla.spsolve(K, R)`.

### 4. Override desde YAML — solo diagnóstico

```yaml
solver:
  type: linear
  linear_algebra: auto    # auto | cholesky | ldlt | lu | cg | gmres
```

`auto` es el default y el 99% de los YAMLs no escribirán este campo. Sirve para tres cosas: (a) benchmarks comparativos, (b) forzar LU si se sospecha bug en la auto-detección de simetría tras un material nuevo, (c) reproducibilidad bit-a-bit de un caso de regresión. **No se documenta como elección de modelado** — se documenta en la sección "Diagnóstico" del manual.

### 5. Fallback automático

Si `CholeskySolver.solve` o `CholeskySolver.factorize` lanza excepción de no-positividad (CHOLMOD: `CholmodNotPositiveDefiniteError`), el solver no lineal:

1. Registra warning con el paso e iteración donde ocurrió.
2. Marca `props.is_positive_definite = False` para el resto del análisis.
3. Reintenta la iteración con `LDLTSolver` (o `LUSolver` mientras LDLᵀ no exista).

Esto convierte un fallo abrupto en una transición silenciosa al cruzar un punto de bifurcación — comportamiento esperado en pandeo y snap-through.

## Hoja de ruta — qué tipos de análisis futuros habilita esta capa

Esta es la parte estratégica del ADR: cada tipo de análisis que el proyecto va a necesitar **se reduce a una capacidad de la capa algebraica**. Diseñar la interfaz hoy con esas capacidades en mente significa que añadir el análisis futuro no requerirá refactor de la capa.

| Análisis futuro | Operación algebraica clave | Capacidad necesaria | Fase ADR |
|---|---|---|---|
| Estático lineal (hoy) | `K·u = f` | `solve(K, b)` | 1 |
| Estático no lineal Newton (hoy) | `K_t·δu = R` repetido | `solve(K, b)` | 1 |
| Estático no lineal Newton modificado | `K_t·δu = R` con `K_t` congelada | `factorize(K)` + `solve(b)` repetido | 2 |
| Arc-length / pandeo postcrítico (hoy) | `K_t` indefinida simétrica | `LDLTSolver` con conteo de pivotes negativos (Sturm) | 2 |
| **Análisis modal** | `K·φ = λ·M·φ` | `EigenSolver` ARPACK + `factorize(K)` para shift-invert | 3 |
| **Pandeo lineal** | `K·φ = λ·K_geom·φ` (autovalor generalizado) | `EigenSolver` + `factorize(K − σ·K_geom)` | 3 |
| **Dinámica implícita** (Newmark, HHT-α) | `[K + a·C + b·M]·δu = R` | `factorize(K_eff)` reutilizable durante muchos pasos si `Δt` constante | 4 |
| **Dinámica explícita** | `M·a = F`, `M` lumped diagonal | Diagonal solve (ya existe en spirit) | 4 |
| **Térmico estacionario** | `K_T·T = Q`, simétrica SPD | `CholeskySolver` (mismo backend que estático) | gratis con fase 1 |
| **Termo-mecánico acoplado monolítico** | `[K_uu K_uT; K_Tu K_TT]·x = R`, no simétrica | `LUSolver` (mismo backend que plasticidad no asociada) | gratis con fase 1 |
| **Termo-mecánico staggered** | Resolver `K_T` y `K_uu` por separado en el mismo paso | Dos instancias independientes de la capa algebraica | gratis con fase 1 |
| **Mallas grandes 3D** (≳ 10⁵–10⁶ DOF) | `K·u = f` con factorización inviable en RAM | `CGSolver` / `GMRESSolver` precondicionado | 5 |
| **Optimización topológica / sensibilidades** | `K·λ = ∂f/∂u` (sistema adjunto) repetido con misma `K` | `factorize(K)` reutilizable | 2 |

**Lectura clave**: la fase 1 (Cholesky + LU + fallback) ya habilita gratis el térmico estacionario y los acoplamientos monolíticos/staggered. La fase 2 (factorize reusable + LDLᵀ) habilita Newton modificado, optimización y mejor diagnóstico de pandeo. La fase 3 (eigen) abre toda la familia modal/pandeo. La fase 4 (dinámica) reusa lo de fase 2. La fase 5 (iterativo) escala a problemas grandes.

## Fases de implementación

Cada fase es un commit cerrado con tests, no bloquea las siguientes y deja Fenix en estado funcional.

### Fase 1 — Cimientos (esta es la orden de trabajo inmediata si apruebas el ADR)

1. Crear `fenix/math/linalg/` con `base.py`, `lu.py`, `dispatcher.py`. `LUSolver` envuelve el `spla.spsolve` actual.
2. Añadir flags `is_symmetric` y `preserves_symmetry` a `Material` y `Element` base con default `True`.
3. Reemplazar las cuatro llamadas `spla.spsolve(...)` en `solvers.py` por `self._linalg.solve(...)`. El despachador se instancia al comienzo de `solve()`.
4. Añadir `CholeskySolver` con import lazy de `scikit-sparse`. Si el import falla, el despachador degrada silenciosamente a `LUSolver` con warning. Fenix sigue funcionando sin la dependencia.
5. Lógica de fallback automático SPD→LU (sección 5 del ADR).
6. Override `linear_algebra` en el parser YAML.
7. Tests:
   - Equivalencia bit-a-bit `LU` (nuevo) vs `spsolve` (viejo) en todos los benchmarks actuales — debe ser idéntico.
   - Equivalencia `LU` vs `Cholesky` en problema elástico SPD a tolerancia de redondeo.
   - Override YAML: forzar `linear_algebra: lu` en problema SPD y verificar coincidencia con `cholesky`.
   - Fallback: problema artificialmente indefinido, verificar que el despachador degrada con warning sin abortar.

### Fase 2 — Factorización reusable + LDLᵀ

8. Implementar `FactorizedSolver` para `LUSolver` (vía `scipy.sparse.linalg.splu`) y para `CholeskySolver` (CHOLMOD ya lo soporta nativo).
9. Flag `freeze_tangent_after_iter` en `NonlinearSolver` (Newton modificado).
10. `LDLTSolver`: envoltorio sobre `pypardiso` o sobre LDLᵀ propio basado en la factorización simbólica de CHOLMOD. Decisión deferida al inicio de fase 2 — ver alternativas en el último apartado.
11. Sturm sequence: tras factorizar con LDLᵀ, exponer `n_negative_pivots` en el resultado del solver. `ArcLengthSolver` lo registra en `SolveResult.diagnostics` para detectar bifurcaciones automáticamente.

### Fase 3 — Eigensolver (modal y pandeo lineal)

12. `EigenSolver` envolviendo `scipy.sparse.linalg.eigsh` (Lanczos para simétrica) y `eigs` (Arnoldi para general). Soporte de shift-invert con `factorize(K − σ·M)`.
13. Solver nuevo `ModalSolver` en `fenix/math/solvers.py` con su propia entrada en el catálogo.
14. Solver nuevo `LinearBucklingSolver` (autovalor generalizado `K·φ = λ·K_geom·φ`).

### Fase 4 — Dinámica implícita

15. `DynamicImplicitSolver` (Newmark-β o HHT-α) que reusa `FactorizedSolver` cuando `Δt` es constante, refactoriza en cambio de paso adaptativo.
16. `DynamicExplicitSolver` con masa diagonal lumped — no necesita factorización; consume la capa algebraica solo trivialmente.

### Fase 5 — Escala (iterativos)

17. `CGSolver` y `GMRESSolver` con precondicionador ILU. Threshold `LARGE_PROBLEM_THRESHOLD` configurable (default ~5·10⁴ DOF en 3D).
18. Heurística de selección automática iterativo↔directo basada en talla y disponibilidad de memoria.

## Consecuencias

**Positivas**

- Coste algebraico ~2× menor en lineal elástico y régimen no lineal estable. Beneficio inmediato sin cambiar API ni YAML.
- Cada análisis futuro mapea a una capacidad de la capa algebraica ya prevista. Añadir modal/pandeo/dinámica no requiere refactor.
- Reglas.md §1 cumplido: añadir un solver algebraico nuevo (paralelo, GPU, externo) es agregar una clase y registrar una regla.
- La elección queda interna y declarativa: nadie tiene que saber álgebra lineal numérica para usar Fenix.
- El override YAML da una ventana de diagnóstico sin convertir la elección en decisión obligatoria.

**Negativas / costes**

- Dependencia nueva: `scikit-sparse` (CHOLMOD nativo). En Windows requiere conda-forge, no instala con `pip` puro. Mitigación: `CholeskySolver` es opcional; si la dependencia falta, Fenix sigue funcionando con `LUSolver`.
- Cada material y elemento debe declarar su flag `is_symmetric` / `preserves_symmetry`. Coste pequeño y unidireccional.
- Branch nuevo en el solver no lineal (asumir SPD, caer a LU al primer fallo). Acotado y testeable.

**Alternativas consideradas**

- *Solver único (status quo)*: rechazado por coste numérico y por bloquear todas las extensiones futuras del catálogo.
- *Elegir solver desde el YAML como decisión obligatoria*: rechazado. Información redundante con material+elemento+análisis; invita a errores físicos sutiles cuando se cambia un material.
- *Detectar simetría inspeccionando `K - Kᵀ` numéricamente*: rechazado. O(nnz) por iteración, depende de tolerancias de redondeo. La vía declarativa es exacta y barata.
- *Empezar con Pardiso (Intel MKL) en lugar de CHOLMOD*: Pardiso cubre LU+LDLᵀ+Cholesky en un binario y es el estándar industrial. **Diferido a fase 2** como decisión: en fase 1 priorizamos stack open-source puro y dependencia ligera; cuando llegue LDLᵀ reabriremos esta puerta. Backend candidato fase 2: `pypardiso` (free, no requiere licencia MKL Pardiso para uso académico) o `MUMPS` vía `python-mumps`.
- *Empezar con MUMPS*: equivalente a Pardiso pero más pesado de instalar en Windows. Mismo razonamiento que arriba.

**Implicaciones para la documentación**

- `arquitectura.md`: añadir `fenix.math.linalg` como subsistema y diagrama de flujo `Solver no lineal → Despachador → Backend algebraico`.
- `conceptos_clave.md`: entrada nueva *"Capa algebraica vs. solver no lineal"* explicando la separación entre las dos capas.
- `catalogo_solvers.md`: no cambia en fase 1 (sigue catalogando `LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`). En fase 3 añade `ModalSolver` y `LinearBucklingSolver`. En fase 4 añade los dinámicos.
- Manual de usuario: sección nueva *"Diagnóstico y override de solver algebraico"* documentando `linear_algebra: <name>` como herramienta de soporte, **no** como elección de modelado.

## Paralelismo: alcance de este ADR y diferimientos

En FEM existen tres niveles de paralelismo, independientes entre sí:

1. **Algebraico** — dentro del backend que resuelve `K·x = y`. Lo aportan las librerías nativas (CHOLMOD, Pardiso, MUMPS, MKL/OpenBLAS) vía OpenMP. No requiere código en Fenix.
2. **Ensamblaje** — cálculo concurrente de `K_e` y `F_int_e` por elemento, con resolución del conflicto de escritura sobre entradas globales compartidas (coloreado de grafo, reducción por hilo, atómicos).
3. **Descomposición de dominio** — partición de la malla entre procesos MPI, cada uno ensambla y resuelve su porción, comunicación solo en interfaces.

**Este ADR cubre exclusivamente el nivel 1**, y lo hace gratis: las cinco fases lo absorben automáticamente porque cada backend (CHOLMOD desde fase 1, Pardiso/MUMPS desde fase 2, ARPACK desde fase 3, iterativos desde fase 5) ya es multihilo cuando se compila contra BLAS paralelo. Sin escribir código adicional, Fenix pasa a usar todos los cores disponibles para la operación más cara (la factorización) en problemas medianos típicos: ~4×–8× speedup en una estación de 8 cores.

**Los niveles 2 y 3 quedan diferidos a ADRs futuros**, en este orden:

- **ADR 0004 (futuro)** — ensamblaje paralelo con coloreado de grafo. Se activa cuando el ensamblaje pase a dominar el coste por paso (> 50% del tiempo). Heurística diagnóstica: instrumentar tiempos `assemble` vs `solve` en `SolveResult.diagnostics`; mientras `solve` domine, este ADR es prematuro.
- **ADR 0005 (futuro lejano)** — descomposición de dominio MPI. Se activa solo si Fenix necesita resolver problemas que no caben en RAM de una máquina (> 10⁶–10⁷ DOF en 3D). Es un esfuerzo grande (parfeap es ~30% del tamaño del FEAP serie). Para el rango de investigación en mecánica de sólidos sobre estación de trabajo individual, Pardiso/MUMPS multihilo en una máquina de 16–32 cores suele ser suficiente.

Coherente con Reglas.md §1: no abstraer hacia paralelismo especulativo sin al menos dos casos reales del dominio que lo pidan. Hoy, la fase 1 ya te da la mejora más barata posible sin abrir frentes nuevos.

## Criterio de cierre del ADR

Se cierra como `Aceptado` cuando el usuario valida el resumen ejecutivo y la hoja de ruta. La fase 1 se ejecuta inmediatamente tras la aceptación bajo este ADR (no requiere spec adicional, porque no es un componente físico nuevo — es plumbing arquitectural cubierto por Reglas.md §3, columna IA-decide). Las fases 2–5 se ejecutan cuando el usuario dé la orden de trabajo correspondiente al activar el análisis que las consume.
