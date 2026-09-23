# Anexo — Capa algebraica

> Esta sección es referencia técnica del subsistema interno que resuelve el sistema lineal `K·x = b` que aparece en el corazón de cada solver no lineal. La elección del backend es **automática**; el usuario no decide. Esta sección documenta cómo se decide y cuáles son los puntos de override para diagnóstico.
>
> Diseño completo: ver ADR 0003 en `docs/adr/0003-despachador-de-solver-algebraico.md`.

## Dos capas que se llaman ambas «solver»

Conviene distinguirlas:

- **Solver no lineal** — `LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`. Orquesta la estrategia de paso, las iteraciones de Newton, el control de longitud de arco, los criterios de convergencia. Es lo que el usuario elige en el YAML con `solver.type`.
- **Capa algebraica** — `solidum.math.linalg`. Resuelve el sistema lineal `K·δU = R` que aparece dentro de cada iteración del solver no lineal. Tiene varios *backends* numéricos y un despachador interno que elige el adecuado según las propiedades de `K`. Es plumbing automático: el usuario no la ve salvo que pida diagnóstico.

## Backends disponibles y criterio de selección

El despachador (`solidum.math.linalg.dispatcher.select_solver`) elige el backend en función de tres flags declarativos sobre `K`: simetría, positividad y talla.

| Régimen de `K` | Backend elegido | Notas |
|---|---|---|
| Simétrica + positiva definida | **Cholesky** (CHOLMOD) | ~2× más rápido y ~2× menos memoria que LU. Requiere la dependencia opcional `scikit-sparse`. Si falta, se pasa al siguiente disponible. |
| Simétrica indefinida | **Pardiso** si está; si no, LU (placeholder LDLᵀ) | Ideal para LDLᵀ con conteo de pivotes negativos (Sturm sequence) en pandeo y snap-through. La interfaz está reservada; mientras no haya backend nativo se resuelve sin pérdida de corrección, solo sin diagnóstico de bifurcación. |
| No simétrica | **Pardiso** si está; si no, **LU general** (SuperLU) | Backend universal. Cubre plasticidad no asociada, follower loads, contacto con fricción. |

El orden de preferencia es **Cholesky → Pardiso → LU**: Cholesky sólo aplica al caso SPD, pero ahí explota la simetría y es la mejor opción; Pardiso cubre el mismo terreno que LU —cualquier matriz real— y lo sustituye cuando está instalado (ADR 0017), porque es multihilo y reordena por disección anidada. LU (SuperLU, en SciPy) es el suelo que siempre está disponible.

**Origen de los flags**, todos derivables del modelo:
- `is_symmetric`: AND lógico de `material.IS_SYMMETRIC` y `element.PRESERVES_SYMMETRY` sobre todos los componentes del dominio. Defaults `True`.
- `is_positive_definite`: arranca en `True` para `LinearSolver` y `NonlinearSolver`; en `False` para `ArcLengthSolver`. Se degrada automáticamente si Cholesky reporta no-positividad.
- `size`: número total de DOFs.

Ningún flag se pide al usuario.

**`EigenSolver` (ADR 0009, problema generalizado)**. La capa algebraica incluye además `solidum.math.linalg.eigen.EigenSolver`, que resuelve el problema generalizado simétrico `K · φ = ω² · M · φ` envolviendo `scipy.sparse.linalg.eigsh` (ARPACK Lanczos con shift-invert centrado en `σ`). No comparte el `Protocol` `LinearAlgebraSolver` de los backends de `K·x = b` porque su firma natural es `solve(K, M, n_modes) → (λ, φ)`. El `ModalSolver` lo invoca para el análisis modal; los cálculos internos de shift-invert generan una factorización de `(K − σ · M)` reutilizada en cada iteración Lanczos, así que la palanca de optimización es la misma — Cholesky enchufado en lugar de SuperLU bajaría 2× el coste del modal en problemas SPD (pendiente, ver memoria de cierre).

## Fallback automático SPD → LU

Si en algún paso `K` deja de ser positiva definida (paso a régimen postcrítico, daño con reblandecimiento), Cholesky lanza `CholeskyNotPositiveDefiniteError`. Los tres solvers no lineales lo capturan y reinstancian el backend con LU general para el resto del análisis. La transición es silenciosa — solo se imprime un mensaje informativo en stdout.

Esto significa que **forzar Cholesky desde YAML como diagnóstico nunca rompe el análisis**: si la matriz se vuelve incompatible, el sistema se autocorrige.

## Factorización reusable y Newton modificado

La interfaz expone `solver.factorize(K)` que retorna un objeto `FactorizedSolver` con `solve(b)` reutilizable. Habilita:

- **Newton modificado**: el `NonlinearSolver` admite el parámetro `freeze_tangent_after_iter: int` que congela la factorización tras N iteraciones del paso. Útil cuando la tangente cambia poco y el coste dominante es la factorización.
- **Dinámica implícita lineal** (ADR 0009 fase 3): el `NewmarkSolver` factoriza una sola vez la matriz efectiva `A_eff = M + γΔt · C + βΔt² · K` al inicio del análisis y reutiliza la factorización en cada paso temporal. Con `Δt` constante y problema lineal, los cientos o miles de pasos del análisis se reducen a sustituciones triangulares baratas. Es el mismo `FactorizedSolver` del despachador.
- **Análisis futuros** (pandeo lineal, dinámica implícita no lineal con factorización congelada entre pasos cuando la tangente cambia poco) reutilizan la misma interfaz sin extensiones adicionales.

Ejemplo YAML:

```yaml
solver:
  type: NonlinearSolver
  num_steps: 10
  tol: 1.0e-8
  freeze_tangent_after_iter: 2
```

## Override desde YAML — herramienta de diagnóstico

Solo en caso de necesidad de *diagnóstico* o *benchmarking*, el bloque `solver` admite el campo opcional `linear_algebra`:

```yaml
solver:
  type: LinearSolver
  linear_algebra: auto         # default; el despachador decide
  # linear_algebra: cholesky   # forzar Cholesky (requiere scikit-sparse)
  # linear_algebra: pardiso    # forzar Pardiso (requiere pypardiso)
  # linear_algebra: lu         # forzar LU (siempre disponible)
  # linear_algebra: ldlt       # placeholder; degrada a LU con warning
```

**Cuándo usar el override:**

- Comparar tiempos entre Cholesky y LU en un mismo problema.
- Forzar LU si se sospecha un bug en la auto-detección de simetría tras añadir un material nuevo.
- Reproducir un caso de regresión bit-a-bit con un backend específico.

**Cuándo NO usarlo**: como parte del setup habitual de un caso de análisis. No es decisión de modelado; el default `auto` cubre todos los regímenes correctamente.

## Backends opcionales: por qué importan

En un análisis grande, **la factorización del sistema es el cuello de botella**, no el ensamblaje. Medido sobre un análisis no lineal completo de `Hex8 15³` (12 288 grados de libertad): el 97 % del tiempo se va en factorizar y el 1 % en ensamblar (ADR 0017). La razón de fondo es el **relleno** (*fill-in*): al factorizar, `L+U` tiene muchos más no-ceros que `K`, y la proporción empeora con la talla — de ×10 a 4 000 grados de libertad a ×42 a 53 000 —, llevándose por delante tiempo y memoria.

Los backends opcionales atacan exactamente eso. Sin ellos Solidum funciona igual, con SuperLU; con ellos los modelos grandes dejan de ser inviables.

### Pardiso — recomendado (`pip install solidum-fem[fast]`)

```bash
pip install solidum-fem[fast]
```

Intel MKL Pardiso: factorización **multihilo** con reordenamiento por **disección anidada**. Trae binarios para Windows, Linux y macOS, así que no necesita compilador. Medido (16 hilos), resolviendo el sistema de `Hex8 n³`:

| Malla | Grados de libertad | SuperLU | Pardiso | Aceleración |
|---|---:|---:|---:|---:|
| Hex8 10³ | 3 993 | 0,47 s | 0,10 s | ×4,7 |
| Hex8 15³ | 12 288 | 3,24 s | 0,27 s | ×12 |
| Hex8 20³ | 27 783 | 18,32 s | 0,78 s | ×23,5 |
| Hex8 25³ | 52 728 | 110,10 s | 1,97 s | ×56 |

La aceleración **crece con el tamaño**, que es donde hace falta. Sobre el análisis no lineal completo de `Hex8 20³` (5 pasos): 212 s → 9,8 s. La solución es la misma a precisión de máquina (`1e-16` relativo en el campo de desplazamientos).

### Cholesky / CHOLMOD (`scikit-sparse`)

```bash
conda install -c conda-forge scikit-sparse
```

Sólo aplica a matrices simétricas y positivas definidas, pero ahí explota la simetría y es la mejor opción. En Windows requiere Visual Studio para compilar, de ahí que la vía práctica sea `conda` y no `pip`.

Instalado cualquiera de los dos, el despachador empieza a usarlo automáticamente — sin tocar YAML, sin recompilar, sin reescribir specs.
