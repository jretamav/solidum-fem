# Anexo — Capa algebraica

> Esta sección es referencia técnica del subsistema interno que resuelve el sistema lineal `K·x = b` que aparece en el corazón de cada solver no lineal. La elección del backend es **automática**; el usuario no decide. Esta sección documenta cómo se decide y cuáles son los puntos de override para diagnóstico.
>
> Diseño completo: ver ADR 0003 en `docs/adr/0003-despachador-de-solver-algebraico.md`.

## Dos capas que se llaman ambas «solver»

Conviene distinguirlas:

- **Solver no lineal** — `LinearSolver`, `NonlinearSolver`, `ArcLengthSolver`. Orquesta la estrategia de paso, las iteraciones de Newton, el control de longitud de arco, los criterios de convergencia. Es lo que el usuario elige en el YAML con `solver.type`.
- **Capa algebraica** — `fenix.math.linalg`. Resuelve el sistema lineal `K·δU = R` que aparece dentro de cada iteración del solver no lineal. Tiene varios *backends* numéricos y un despachador interno que elige el adecuado según las propiedades de `K`. Es plumbing automático: el usuario no la ve salvo que pida diagnóstico.

## Backends disponibles y criterio de selección

El despachador (`fenix.math.linalg.dispatcher.select_solver`) elige el backend en función de tres flags declarativos sobre `K`: simetría, positividad y talla.

| Régimen de `K` | Backend elegido | Notas |
|---|---|---|
| Simétrica + positiva definida | **Cholesky** (CHOLMOD) | ~2× más rápido y ~2× menos memoria que LU. Requiere la dependencia opcional `scikit-sparse`. Si falta, degrada a LU con un *warning* una sola vez. |
| Simétrica indefinida | LU general (placeholder LDLᵀ) | Ideal para LDLᵀ con conteo de pivotes negativos (Sturm sequence) en pandeo y snap-through. La interfaz está reservada; mientras no haya backend nativo, se usa LU sin pérdida de corrección, solo sin diagnóstico de bifurcación. |
| No simétrica | **LU general** (SuperLU) | Backend universal. Cubre plasticidad no asociada, follower loads, contacto con fricción. |

**Origen de los flags**, todos derivables del modelo:
- `is_symmetric`: AND lógico de `material.IS_SYMMETRIC` y `element.PRESERVES_SYMMETRY` sobre todos los componentes del dominio. Defaults `True`.
- `is_positive_definite`: arranca en `True` para `LinearSolver` y `NonlinearSolver`; en `False` para `ArcLengthSolver`. Se degrada automáticamente si Cholesky reporta no-positividad.
- `size`: número total de DOFs.

Ningún flag se pide al usuario.

**`EigenSolver` (ADR 0009, problema generalizado)**. La capa algebraica incluye además `fenix.math.linalg.eigen.EigenSolver`, que resuelve el problema generalizado simétrico `K · φ = ω² · M · φ` envolviendo `scipy.sparse.linalg.eigsh` (ARPACK Lanczos con shift-invert centrado en `σ`). No comparte el `Protocol` `LinearAlgebraSolver` de los backends de `K·x = b` porque su firma natural es `solve(K, M, n_modes) → (λ, φ)`. El `ModalSolver` lo invoca para el análisis modal; los cálculos internos de shift-invert generan una factorización de `(K − σ · M)` reutilizada en cada iteración Lanczos, así que la palanca de optimización es la misma — Cholesky enchufado en lugar de SuperLU bajaría 2× el coste del modal en problemas SPD (pendiente, ver memoria de cierre).

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
  # linear_algebra: lu         # forzar LU (siempre disponible)
  # linear_algebra: ldlt       # placeholder; degrada a LU con warning
```

**Cuándo usar el override:**

- Comparar tiempos entre Cholesky y LU en un mismo problema.
- Forzar LU si se sospecha un bug en la auto-detección de simetría tras añadir un material nuevo.
- Reproducir un caso de regresión bit-a-bit con un backend específico.

**Cuándo NO usarlo**: como parte del setup habitual de un caso de análisis. No es decisión de modelado; el default `auto` cubre todos los regímenes correctamente.

## Activación de Cholesky (opcional)

```bash
conda install -c conda-forge scikit-sparse
```

Una vez instalado, los problemas estáticos lineales y los Newton estables empiezan a usar Cholesky automáticamente — sin tocar YAML, sin recompilar, sin reescribir specs. Si la dependencia no está disponible, Fenix funciona idéntico con LU.
