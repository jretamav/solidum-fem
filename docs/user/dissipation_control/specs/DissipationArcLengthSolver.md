# DissipationArcLengthSolver — Arc-length con restricción de disipación (Gutiérrez 2004)

> **Módulo de usuario `dissipation_control`** ([ADR 0020](../../../adr/0020-modulos-de-usuario.md)): código en [`solidum/user/dissipation_control/`](../../../../solidum/user/dissipation_control/), tests en [`tests/user/dissipation_control/`](../../../../tests/user/dissipation_control/). `import solidum` no lo carga: en YAML, `user_modules: [dissipation_control]`; en Python, `solidum.load_user_module("dissipation_control")`. Estuvo en el programa principal (`docs/specs/`, `solidum/math/solvers/`) del 2026-05-19 al 2026-09-29.
>
> **Clasificación** (criterio de admisión, Reglas.md §1, 2026-09-29): **no estándar**. La condición (ii) no se cumple: ningún código general establecido implementa la restricción de Gutiérrez de forma nativa. DIANA ofrece como restricciones del arco el *Spherical Path* y el *Updated Normal Plane*, y su método «energy based» es una regla de tamaño de paso (fija el incremento de carga con el producto de los incrementos del primer predictor), no esta restricción; los `TYPE` del `PILOTAGE` de Code_Aster (`DDL_IMPO`, `LONG_ARC`, `ANA_LIM`, `DEFORMATION`, `PRED_ELAS`, `SAUT_IMPO`, `SAUT_LONG_ARC`) no la incluyen; en Abaqus sólo existen implementaciones de investigación mediante subrutinas de usuario.
>
> **Spec corta tipo extensión** (Reglas.md §4). Variante de [`ArcLengthSolver`](../../../specs/ArcLengthSolver.md) (programa principal) que **solo cambia la restricción del paso**: pasa de cilíndrica ($\lVert\Delta\mathbf U\rVert^2 = dl^2$) a basada en **disipación de energía** ($g(\Delta\mathbf U, \Delta\lambda) = \tau$). Reusa todo el plumbing del padre (predictor tangente, corrector Newton, adaptatividad de paso, Dirichlet/MPC, backend algebraico) a través de su parte pública: la clase `ArcLengthSolver` y su problema de Newton por paso `solidum.math.solvers.ArcProblem`. Sin ADR nuevo: reusa el contrato arquitectural ya aceptado del arc-length.
>
> Motivación: el `ArcLengthSolver` cilíndrico no atraviesa la transición elástico→softening con penalty stiff (`K_e` de la discontinuidad embebida, ADR 0010 fase 4). La restricción cuadrática queda casi degenerada cerca del pico — el discriminante de la cuadrática tiende a cero, el solver bisecta indefinidamente y aborta. *(Nota 2026-09-29: ese diagnóstico es anterior a la corrección de la tangente cohesiva del 2026-09-23, commit `597f3d9`, que presentaba a los solvers una tangente de signo contrario en casi toda la rama de ablandamiento, y no se ha vuelto a comprobar; con control por desplazamiento, `NonlinearSolver` sí atraviesa hoy esa rama. Ver el caveat de [`CST_Embedded2D`](../../discontinuities/specs/CST_Embedded2D.md).)* La restricción de disipación es **lineal en $(\Delta\mathbf U, \Delta\lambda)$** y se activa **proporcionalmente a la energía consumida** por la grieta cohesiva, navegando la rama post-pico de forma estable. Esta variante apuntaba a desbloquear el benchmark de fractura faithful Van Vliet ([`CST_Embedded2D`](../../discontinuities/specs/CST_Embedded2D.md) §"Out of scope"), el balance energético end-to-end de $G_F$ (hueco abierto en [[project_validacion_fase_b]]), y validación cuantitativa de cohesivos en pipeline real. Esos benchmarks usan componentes de **otro módulo de usuario**, `discontinuities` (elemento embebido y ley cohesiva): un modelo que los combine con este solver declara los dos, `user_modules: [discontinuities, dissipation_control]`.

---

## Especificación física

### 0. Descripción general

Solver estático no lineal con **traza por disipación**: la magnitud que se mantiene constante de paso a paso ya no es la longitud euclídea $\lVert\Delta\mathbf U\rVert$ sino la **disipación incremental de energía** $\Delta E_d$ asociada al paso. Adecuado para:

- Fractura cohesiva (Mode-I y, en el futuro, mixto) con penalty `K_e` stiff.
- Daño con softening severo donde la rama post-pico tiene pendiente casi vertical en el plano $(U, \lambda)$.
- Localización de deformación plástica con bandas estrechas.

Régimen elástico (sin disipación): la restricción degenera a $g = 0$. El solver detecta este caso y **conmuta automáticamente** a arc-length cilíndrico durante el régimen elástico, retornando a control por disipación cuando la disipación del paso $\Delta E_d$ (§1) supera un umbral relativo (§5). La detección es global, por energía: el solver no consulta las variables internas del material (daño, deformación plástica, apertura de la grieta). El alcance real en plasticidad es un punto abierto (§5).

### 1. Ecuación de equilibrio resuelta

Idéntica al padre:

$$\mathbf F_{\text{int}}(\mathbf U) = \lambda\,\mathbf F_{\text{ext}}^{\text{ref}}$$

con $\lambda$ incógnita. La **restricción del paso** cambia respecto al padre.

#### Restricción de disipación (Gutiérrez 2004, Verhoosel et al. 2009)

La energía disipada por paso, en problemas conservativos en el componente elástico, es:

$$\Delta E_d = \frac{1}{2}\left(\lambda_n\,\mathbf F_{\text{ext}}^{\text{ref}T}\,\Delta\mathbf U \;-\; \Delta\lambda\,\mathbf F_{\text{ext}}^{\text{ref}T}\,\mathbf U_n\right)$$

donde $(\mathbf U_n, \lambda_n)$ son los valores committed al inicio del paso y $(\Delta\mathbf U, \Delta\lambda)$ los incrementos del paso. La restricción es:

$$g(\Delta\mathbf U, \Delta\lambda) \;\equiv\; \frac{1}{2}\left(\lambda_n\,\mathbf F_{\text{ext}}^{\text{ref}T}\,\Delta\mathbf U \;-\; \mathbf U_n^T\,\mathbf F_{\text{ext}}^{\text{ref}}\,\Delta\lambda\right) \;=\; \tau$$

con $\tau$ el **incremento de disipación objetivo** (parámetro del solver, papel análogo al `initial_dl` del padre).

**Linealidad**: $g$ es lineal en $(\Delta\mathbf U, \Delta\lambda)$. Por eso el corrector tiene **una sola raíz**, no dos como en la cuadrática cilíndrica — eliminando la selección por menor ángulo y la patología de raíces imaginarias en pasos grandes.

### 2. Condiciones de contorno

Idénticas al padre (`ArcLengthSolver` §2). Dirichlet homogéneo/no homogéneo y MPC vía ADR 0004.

### 3. Salidas físicas

Idénticas al padre: `SolveResult` con $\mathbf U$ final y `step_callback(step, U, λ)` para la historia, con la misma firma que en el padre. El solver **no expone** la disipación del paso ni la acumulada $\sum_n \Delta E_d^{(n)}$: la calcula tras cada paso convergido para decidir el modo (§5) y sólo la escribe en el registro (`ΔE_d=…` en el log del solver). Quien la necesite —por ejemplo, para contrastar **$G_F$ de extremo a extremo** con $\int_{\Gamma_d} t\cdot d[\![u]\!] \approx \sum_n \Delta E_d^{(n)}$— la reconstruye desde la historia $(\mathbf U_n, \lambda_n)$ del callback con la fórmula de $\Delta E_d$ de §1, como hace el test 5 de [`test_dctl_dissipation_arclength.py`](../../../../tests/user/dissipation_control/test_dctl_dissipation_arclength.py).

---

## Formulación numérica — solo lo que cambia respecto al padre

### 4. Esquema operativo

Estructura idéntica al padre con tres diferencias: el predictor y la restricción del corrector dependen del modo, y tras cada paso convergido se decide el modo del siguiente (§5). El esquema describe el `solve` del código (`dissipation_arclength.py`):

```
λ = 0;  U = 0;  ΔU_prev = 0
τ = initial_tau          # análogo a initial_dl, pero en unidades de energía
modo = "cylindrical"      # siempre arranca así: con λ_n = 0, U_n = 0 es α = 0
dl = (sin fijar)
mientras λ < max_lambda y step < max_steps:
    step += 1
    prepare_all_steps(U)                  # gancho de inicio de paso (ADR 0020, P5)

    # Predictor tangente (el del padre)
    δu_t = K_t^{-1} · F_ext^ref
    si K_t singular:                      # el predictor no puede resolver
        si dl fijado: dl /= 2
        τ /= 2                            # se dividen los dos, sin comprobar cotas
        continuar                         # reintento del paso
    si dl sin fijar: dl = dl_ref = Δl₁     # primer paso del padre (initial_dlambda / initial_dl)
    sign = sign(ΔU_prev · δu_t)           # regla del padre

    modo_paso = modo
    si modo_paso == "dissipation":
        α = ½·(λ_n·F_ext^ref·δu_t − F_ext^ref·U_n)
        escala = |½·λ_n·F_ext^ref·δu_t|
        si (escala > 0 y |α| < 1e-6·escala) o |α| < ZERO_TOL:
            modo_paso = modo = "cylindrical"   # régimen lineal: α no es significativo
    si modo_paso == "cylindrical":
        Δλ_pred = sign · dl / (‖δu_t‖ + ZERO_TOL)
    sino:
        Δλ_pred = sign · |τ/α|            # signo del padre, no el de τ/α (punto abierto, §5)
    ΔU_pred = Δλ_pred · δu_t

    # Corrector iterativo (NewtonCorrector compartido, ADR 0015)
    para iter en range(max_iter):
        K_t, F_int = assemble(...)
        R = (λ + Δλ)·F_ext^ref − F_int
        δu_R = K_t^{-1} · R
        δu_t = K_t^{-1} · F_ext^ref
        si modo_paso == "cylindrical":
            # Restricción cuadrática, idéntica al padre (ArcLengthSolver §5)
            ddλ = raiz_de_menor_angulo(...)
        sino:
            # Restricción lineal de Gutiérrez: g(ΔU + δu_R + ddλ·δu_t, Δλ + ddλ) = τ
            num = τ − ½·(λ_n·F_ext^ref·(ΔU+δu_R) − Δλ·F_ext^ref·U_n)
            den = ½·(λ_n·F_ext^ref·δu_t − F_ext^ref·U_n)     # α con el δu_t de la iteración
            si |den| < ZERO_TOL: abortar el paso              # CorrectionAborted → bisección
            ddλ = num / den
        ΔU += δu_R + ddλ·δu_t
        Δλ += ddλ
        si converge: break

    si converge y λ + Δλ > max_lambda:    # llegada exacta del padre (ArcLengthSolver._land)
        repetir desde (U_n, λ_n) con λ = max_lambda fijo
        si converge: consolidar y continuar     # sin decidir modo ni adaptar dl/τ
        sino: tratar como paso fallido
    si converge:
        ΔE_d = ½·(λ_n·F_ext^ref·ΔU − Δλ·F_ext^ref·U_n)       # paso completo, tras corregir
        umbral = dissipation_threshold · max(‖F_ext^ref‖·‖U_n+1‖, 1)
        si modo == "cylindrical" y ΔE_d > umbral:
            modo = "dissipation";  τ = max(ΔE_d, initial_tau)
        sino si modo == "dissipation" y ΔE_d < umbral:
            modo = "cylindrical"
        adaptar dl (modo cilíndrico) o τ (modo disipación) según el modo del SIGUIENTE paso (§9)
        consolidar: U += ΔU;  λ += Δλ;  ΔU_prev = ΔU
    sino:                                 # bisección según el modo del paso fallido
        si modo_paso == "cylindrical":
            dl *= 0.5;  si dl < ARCLENGTH_MIN_DL_FACTOR·dl_ref: abortar (RuntimeError)
        sino:
            τ *= 0.5;   si τ < ARCLENGTH_MIN_DL_FACTOR·initial_tau: abortar (RuntimeError)
```

El paso fallido reduce τ a la mitad (`τ *= 0.5`); `tau_shrink_factor` sólo actúa sobre pasos que convergen lentamente (§9).

### 5. Predictor / corrector — punto clave del cambio

#### Switching automático cilíndrico ↔ disipación

El solver mantiene un atributo `modo` con dos estados:

- `"cylindrical"`: arc-length cilíndrico estándar; restricción $\lVert\Delta\mathbf U\rVert^2 = dl^2$. **Régimen elástico**, sin disipación activa.
- `"dissipation"`: arc-length por energía; restricción $g = \tau$. **Régimen disipativo**.

**Criterio de transición** (inspirado en Verhoosel et al. 2009 §3.2; es el que implementa el código). El único indicador es la disipación del paso $\Delta E_d$ de §1, evaluada con el paso completo tras la corrección; el solver no consulta las variables internas del material (daño, deformación plástica, apertura). Tras cada paso convergido que no es la llegada a `max_lambda`:

$$\varepsilon_{\text{diss}} = \texttt{dissipation\_threshold}\cdot\max\bigl(\lVert\mathbf F_{\text{ext}}^{\text{ref}}\rVert\,\lVert\mathbf U_{n+1}\rVert,\ 1\bigr)$$

- En `"cylindrical"`: si $\Delta E_d > \varepsilon_{\text{diss}}$, el siguiente paso entra en `"dissipation"` y $\tau \leftarrow \max(\Delta E_d, \texttt{initial\_tau})$, para dar continuidad al primer paso disipativo.
- En `"dissipation"`: si $\Delta E_d < \varepsilon_{\text{diss}}$ (disipación nula o negativa, p. ej. descarga), el siguiente paso vuelve a `"cylindrical"`.

Además, **en el predictor** de un paso en `"dissipation"`, con $\alpha = \tfrac12(\lambda_n\,\mathbf F_{\text{ext}}^{\text{ref}}\cdot\delta\mathbf u_t - \mathbf F_{\text{ext}}^{\text{ref}}\cdot\mathbf U_n)$ (denominador de la restricción): si $|\alpha| < 10^{-6}\cdot\tfrac12|\lambda_n\,\mathbf F_{\text{ext}}^{\text{ref}}\cdot\delta\mathbf u_t|$, o $|\alpha|$ es menor que `ZERO_TOL`, el paso se hace en cilíndrico y el modo pasa a `"cylindrical"`. En un problema lineal monótono $\alpha \equiv 0$ exactamente ($\mathbf F\cdot\mathbf U_n = \lambda_n\,\mathbf F\cdot\delta\mathbf u_t$), así que esta detección es la matemática del régimen sin disipación, no un error. Tras ese paso vuelve a aplicarse el criterio de $\Delta E_d$, que puede devolver el modo a `"dissipation"`. Dentro del corrector, $|\alpha| < $ `ZERO_TOL` con el $\delta\mathbf u_t$ de la iteración aborta el paso (`CorrectionAborted`), que se biseca.

No hay criterio por número de iteraciones: la spec describía hasta el 2026-09-29 una vuelta a cilíndrico con un `tau_relax_iter_threshold` que el código no tiene.

#### Puntos físicos abiertos (pendientes del usuario, 2026-09-29)

Registrados al trasladar el solver al módulo; la spec describe el código tal como está y **no los resuelve**:

1. **Signo del predictor en modo de disipación.** El código usa $\Delta\lambda_\text{pred} = \text{sign}\cdot|\tau/\alpha|$, con el signo del producto escalar del padre, en lugar de $\Delta\lambda_\text{pred} = \tau/\alpha$, cuyo signo fijaría la propia restricción. En el predictor $g = \alpha\,\Delta\lambda_\text{pred}$, así que cuando $\text{sign} \neq \operatorname{sgn}(\alpha)$ el predictor arranca en $g = -\tau$ y es el corrector quien lo lleva a $+\tau$.
2. **Piso absoluto de la escala del umbral.** $\max(\lVert\mathbf F_{\text{ext}}^{\text{ref}}\rVert\,\lVert\mathbf U\rVert,\ 1)$ tiene un piso de 1 en unidades de energía del modelo: el umbral $\varepsilon_{\text{diss}}$ depende del sistema de unidades cuando $\lVert\mathbf F\rVert\,\lVert\mathbf U\rVert < 1$.
3. **Alcance en plasticidad.** §0 y `requirements` declaran el solver apto para plasticidad, pero la $\Delta E_d$ de §1 es la expresión de daño con descarga secante; el código no tiene otra.

#### Orientación del predictor (sin seguimiento de inercia)

El sentido del predictor sale de la regla del padre: el signo del producto escalar del predictor tangente con el incremento del paso anterior, $\text{sign}(\Delta\mathbf U_{\text{prev}}\cdot\delta\mathbf u_t)$ (Crisfield 1991, cap. 9), tanto en modo cilíndrico como en modo de disipación (§4). El solver **no sigue la inercia de $\mathbf K_t$** (número de pivotes negativos), y ninguna decisión del paso depende de ella.

Un seguimiento de inercia exigiría una factorización LDLᵀ de Bunch-Kaufman que exponga el número de pivotes negativos (ley de inercia de Sylvester): permitiría detectar el paso por un punto límite y orientar el predictor por el cambio de signo de un pivote, no sólo por el producto escalar. Es la deuda técnica #7 de [STATUS.md](../../../STATUS.md). El caso patológico de bifurcación simétrica, en el que el "menor ángulo" elige una rama paralela (caveat documentado en `ArcLengthSolver` §10), queda por tanto sin resolver aquí.

*Historia*: hasta el 2026-09-29 el solver tenía un método `_negative_pivots()` que aproximaba la paridad del número de pivotes negativos con el signo del determinante de $\mathbf K_t$ (densificando la matriz), y la spec afirmaba que orientaba el predictor. Ningún `solve` lo llamaba; se retiró junto con el placeholder homónimo de `ArcLengthSolver`, sin cambio de resultados.

### 6. Criterio de convergencia (ADR 0007)

Idéntico al padre — patrón dual fuerza + desplazamiento, escala autoderivada.

### 7. Imposición de Dirichlet y MPC

Idéntico al padre. La restricción de disipación opera sobre DOFs libres después de la reducción Dirichlet/MPC.

### 8. Backend algebraico

Idéntico al padre (`is_positive_definite=False`, despachador LU). Sin adiciones: el solver no pide al backend ningún dato de inercia. Un seguimiento de inercia (§5) requeriría LDLᵀ de Bunch-Kaufman con conteo de pivotes negativos; hoy las factorizaciones exponen `n_negative_pivots` sólo en Cholesky (0 por construcción) y lo dejan en `None` en LU, Pardiso y el iterativo.

### 9. Adaptatividad — ajuste de $\tau$

Análoga a la adaptatividad de $dl$ del padre, con $\tau$ en lugar de $dl$:

- Converge en $<$ `tau_grow_iter_threshold` iteraciones ⇒ $\tau \to \tau \cdot $ `tau_grow_factor` (cota: `initial_tau · tau_max_factor`).
- Converge en $>$ `tau_shrink_iter_threshold` iteraciones ⇒ $\tau \to \tau \cdot $ `tau_shrink_factor`.
- No converge en `max_iter` (o el corrector aborta, o falla la llegada a `max_lambda`) ⇒ $\tau \to \tau/2$ (`τ *= 0.5`, no `tau_shrink_factor`) y reintentar.
- Cota inferior: $\tau < $ `ARCLENGTH_MIN_DL_FACTOR · initial_tau` ⇒ aborto (el mismo factor que acota $dl$ en el padre, `solidum/constants.py`).
- Tangente singular en el predictor ⇒ se dividen a la vez $dl$ (si ya está fijado) y $\tau$ entre 2, sea cual sea el modo, y se reintenta el paso; aquí no se comprueba la cota inferior.
- Al pasar de cilíndrico a disipación, $\tau \leftarrow \max(\Delta E_d, \texttt{initial\_tau})$ (§5).

Las iteraciones que se comparan con los umbrales son las resoluciones del corrector en el paso (`n_solves`). La adaptatividad se aplica según el modo del **siguiente** paso, ya decidido: en `"cylindrical"` (régimen elástico inicial o vuelta a cilíndrico) opera sobre $dl$, idéntica al padre; en `"dissipation"`, sobre $\tau$.

### 10. Caveats numéricos

Heredados del padre (`ArcLengthSolver` §10) más los específicos:

- **Régimen disipativo no estrictamente monótono**: si el material vuelve transitoriamente a régimen elástico (descarga local), $\Delta E_d$ puede ser negativo. El solver tolera este caso pero conmuta a `"cylindrical"` si la disipación neta del paso cae por debajo del umbral.
- **Inicialización de $\tau$**: debe ser del orden de la disipación esperada por paso, es decir, una fracción de la disipación total esperada en la rama post-pico (en fractura cohesiva, típicamente una fracción pequeña de $G_F \cdot l_d$ por elemento agrietado). Demasiado pequeño ⇒ pasos numerosos; demasiado grande ⇒ saltos sobre la rama post-pico.
- **Disipación de origen no físico**: si el modelo numérico introduce disipación espuria (algorithmic damping en transitorio, regularización no consistente), $\Delta E_d$ contiene esa componente. **Este solver es estático**: si el usuario invoca `DissipationArcLengthSolver` sobre un modelo con material no disipativo (elástico puro), el modo permanece en `"cylindrical"` indefinidamente — comportamiento equivalente al padre. No es un caveat operacional, es la definición del switching.
- **MPC con offset $\mathbf g$ no homogéneo**: la restricción de disipación se evalúa sobre DOFs libres reducidos; si las MPC introducen offset, $g$ se evalúa sobre el sistema reducido y la disipación medida es consistente con la energía aportada por las MPC (tratamiento idéntico al residuo en el padre).

---

## Contrato de implementación

```yaml
name: DissipationArcLengthSolver
kind: solver
status: implemented      # draft → implemented → validated. Validación parcial — ver §"Implementación".

parent_spec: ArcLengthSolver     # programa principal: docs/specs/ArcLengthSolver.md

interface:
  yaml_type: dissipation_arclength
  output: SolveResult
  pipeline_kind: static

parameters:
  # Heredados del padre (todos con misma semántica):
  - { name: convergence,              type: ConvergenceCriterion, required: false, default: "ConvergenceCriterion()" }
  - { name: max_iter,                 type: int,   required: false, default: 20 }
  - { name: max_lambda,               type: float, required: false, default: 1.0 }
  - { name: max_steps,                type: int,   required: false, default: 100 }
  - { name: initial_dlambda,          type: float, required: false, default: 0.1,
      desc: "Heredado: fracción de la carga de referencia del primer paso (modo cilíndrico); Δl₁ = Δλ₁·‖K₀⁻¹F_ref‖. Excluyente con initial_dl" }
  - { name: initial_dl,               type: float, required: false, default: null,
      desc: "Heredado: longitud de arco explícita del primer paso, en unidades de desplazamiento. Excluyente con initial_dlambda" }
  - { name: dl_grow_factor,           type: float, required: false, default: 1.5 }
  - { name: dl_max_factor,            type: float, required: false, default: 5.0 }
  - { name: dl_shrink_factor,         type: float, required: false, default: 0.6 }
  - { name: dl_grow_iter_threshold,   type: int,   required: false, default: 4 }
  - { name: dl_shrink_iter_threshold, type: int,   required: false, default: 8 }
  - { name: linear_algebra,           type: str,   required: false, default: "auto",
      desc: "Heredado de ArcLengthSolver: 'auto' (Cholesky → Pardiso → LU según simetría y dependencias instaladas), 'cholesky', 'pardiso', 'lu', 'iterative[:amg|jacobi|none]' (ADR 0003, 0017, 0018)" }

  # Específicos de esta variante:
  - { name: initial_tau,              type: float, required: true,
      desc: "Incremento de disipación objetivo (unidades de energía: F·U). Calibrar como fracción de la disipación total esperada en la rama post-pico, del orden de E_d,total/n_pasos_post_pico (en fractura cohesiva, E_d,total ≈ G_F·l_d por elemento agrietado)" }
  - { name: tau_grow_factor,          type: float, required: false, default: 1.5 }
  - { name: tau_max_factor,           type: float, required: false, default: 5.0 }
  - { name: tau_shrink_factor,        type: float, required: false, default: 0.6 }
  - { name: tau_grow_iter_threshold,  type: int,   required: false, default: 4 }
  - { name: tau_shrink_iter_threshold,type: int,   required: false, default: 8 }
  - { name: dissipation_threshold,    type: float, required: false, default: 1.0e-12,
      desc: "Umbral relativo bajo el cual ΔE_d se considera 0 (mantiene modo cilíndrico)" }

requirements:
  - "Modelo con material disipativo (plasticidad, daño, cohesivo) — en modelo puramente elástico el solver opera equivalente a ArcLengthSolver cilíndrico"
  - "F_ext^ref no nulo (la restricción de disipación se evalúa sobre F_ext^ref·U)"
  - "Heredados del padre: apoyos suficientes; sin modos rígidos; dl_inicial razonable"
  - "Módulo de usuario cargado: user_modules: [dissipation_control] en el YAML o solidum.load_user_module('dissipation_control')"

conventions:
  units: "heredadas del modelo (ADR 0008); τ en unidades de energía consistentes"
  stability: "incondicional respecto a puntos límite y rama post-pico de softening severo. Ante bifurcaciones, la misma que el padre: el sentido del predictor sale del producto escalar con el incremento previo (sin seguimiento de inercia, §5)"

out_of_scope:
  - "Bifurcaciones — el solver no sigue la inercia de K_t (sin sign-of-pivot tracking, deuda #7); una bifurcación, simple o múltiple, requiere perturbación explícita"
  - "Análisis dinámico — la formulación de disipación es válida solo en régimen cuasi-estático conservativo"
  - "Modo II y mixto I-II en cohesivo — heredado del scope ADR 0010 fases F+"
  - "3D — coherente con el scope global del proyecto en 2026-05"

acceptance:
  verification:
    - name: recuperacion_del_arclength_cilindrico_en_elastico
      setup: "Modelo elástico puro (sin disipación); ejecutar DissipationArcLengthSolver y ArcLengthSolver con los mismos parámetros y misma malla"
      expect: "El modo permanece 'cylindrical'; trayectoria U(λ) coincide a paridad de bits con el padre"
      tol_rel: 1.0e-12

    - name: snap_through_lee_frame
      setup: "Mismo benchmark que ArcLengthSolver §acceptance — frame de Lee, dos barras corotacionales"
      expect: "Atraviesa el punto límite; U-λ coincide con Crisfield 1991 Fig. 9.5"
      tol_rel: 1.0e-3

    - name: damage_softening_post_pico
      setup: "Quad4 con IsotropicDamage2D, tracción que entra en rama descendente"
      expect: "Modo conmuta a 'dissipation' al detectar Δω>0; rama post-pico trazada; coincide con ArcLengthSolver del padre dentro de tolerancia (ambos resuelven el mismo problema)"
      tol_rel: 1.0e-3

    # Los tres casos siguientes usan componentes de OTRO módulo de usuario,
    # `discontinuities` (CST_Embedded2D + CohesiveDamageIsotropic): el modelo
    # declara user_modules: [discontinuities, dissipation_control] y su test
    # irá en tests/user/dissipation_control/, cargando los dos módulos.
    - name: cohesive_embedded_softening_postpeak       # benchmark CRÍTICO — desbloquea fase 4 ADR 0010
      setup: "Requiere el módulo de usuario discontinuities. Probeta rectangular con CST_Embedded2D + CohesiveDamageIsotropic; tracción uniaxial monotónica que activa Rankine y entra en rama post-pico con penalty K_e stiff (1e13)"
      expect: "El solver atraviesa la transición elástico→softening (régimen donde el ArcLengthSolver cilíndrico aborta); cierra el ciclo de softening; daño cohesivo committed alcanza saturación"
      tol_rel: 1.0e-3

    - name: balance_energetico_GF
      setup: "Requiere el módulo de usuario discontinuities. Mismo modelo cohesivo+embedded del test anterior; integrar ΣΔE_d sobre toda la historia hasta saturación (ω→1) y comparar con G_F · longitud_grieta_abierta del material"
      expect: "Σ ΔE_d ≈ G_F · l_d a tolerancia de discretización temporal (1% típico)"
      tol_rel: 1.0e-2

    - name: convergencia_temporal_tau
      setup: "Requiere el módulo de usuario discontinuities. Reducir initial_tau en factor 2 sucesivamente sobre el benchmark cohesivo, observar U(λ_final)"
      expect: "Convergencia primer orden en τ (esquema lineal en restricción): error → 0 como O(τ)"
      tol_rel: 0.2

  specific:
    - name: switching_cilindrico_a_disipacion_en_pico
      setup: "Modelo con daño que cruza σ_t0 en el paso N"
      expect: "Modos: 'cylindrical' en pasos 1..N-1; 'dissipation' desde N+1 en adelante; sin oscilación de modo en pasos contiguos"

    - name: sign_of_pivot_tracking_bifurcacion_simple     # no implementado: requiere LDLᵀ (deuda #7), §5
      setup: "Modelo con un punto límite simple (snap-back); registrar el número de pivots negativos antes y después del paso del límite"
      expect: "Conteo de pivots negativos cambia por 1 al atravesar el límite; predictor mantiene la dirección física"

    - name: fallback_a_cilindrico_en_descarga_global
      setup: "Modelo con grieta saturada (ω≈1) y descarga: λ decrece, F_ext invertida"
      expect: "Modo conmuta de vuelta a 'cylindrical' al detectar ΔE_d < umbral"

references:
  - "Gutiérrez M.A. (2004). Energy release control for numerical simulations of failure in quasi-brittle solids. *Communications in Numerical Methods in Engineering* 20(1), 19-29. (Restricción lineal de disipación; fórmula central de esta spec.)"
  - "Verhoosel C.V., Remmers J.J.C., Gutiérrez M.A. (2009). A dissipation-based arc-length method for robust simulation of brittle and ductile failure. *International Journal for Numerical Methods in Engineering* 77(9), 1290-1321. (Versión robusta con switching y aplicación a cohesivos.)"
  - "Crisfield M.A. (1991). *Non-linear Finite Element Analysis of Solids and Structures*, Vol. 1, Wiley. §9.4 (alternative constraints to cylindrical/spherical)."
  - "de Borst R., Sluys L.J. (1999). *Computational Methods in Non-linear Solid Mechanics*, TU Delft. §3.5 (controlled descent)."
  - "ADR 0010 §Hoja de ruta (fase 4 diferida por solver)."
  - "ADR 0020 (módulos de usuario) y criterio de admisión de Reglas.md §1."
  - "Spec padre: docs/specs/ArcLengthSolver.md (programa principal)."
```

---

## Implementación

- **Archivo**: [`solidum/user/dissipation_control/dissipation_arclength.py`](../../../../solidum/user/dissipation_control/dissipation_arclength.py) (hasta el 2026-09-29, `solidum/math/solvers/dissipation_arclength.py`).
- **Clase**: `DissipationArcLengthSolver(ArcLengthSolver)` — subclase. Reusa predictor, corrector, ajuste de paso, manejo de Dirichlet/MPC y backend del padre. Override completo de `solve()` con bifurcación por `self._mode`. Su problema por paso, `_DissipationArcProblem`, hereda de `solidum.math.solvers.ArcProblem` (el problema de Newton del arco que el programa principal exporta como punto de extensión) y sólo redefine `constraint`.
- **Atributo `self._mode`**: `"cylindrical" | "dissipation"`. Transición tras commit según umbral `dissipation_threshold·max(‖F_ref‖·‖U‖, 1)` (§5; el piso 1 es un punto abierto).
- **Orientación del predictor**: la del padre, `sign(ΔU_prev · du_t)`, en los dos modos. Sin seguimiento de inercia (§5); el método `_negative_pivots()` (signo del `slogdet` de `K_t`, que ningún `solve` llamaba) se retiró el 2026-09-29.
- **Llegada a `max_lambda`** (desde 2026-09-23, la del padre): todos los pasos llevan su restricción (cilíndrica o de disipación); el que converge por encima de `max_lambda` no se consolida y se repite con λ = `max_lambda` fijo desde el estado anterior (`ArcLengthSolver._land`). Si la llegada falla se biseca dl o τ según el modo del paso. Sustituye a la salvaguarda anterior contra el `final_step` prematuro, que bisecaba cuando el predictor rebasaba `max_lambda` más de 3 veces lo que faltaba y, por debajo de ese factor, imponía `max_lambda` en control de carga.
- **Detección de α≈0 en el predictor**: `|α| < 1e-6·½·|λ_n·F·du_t|` (umbral relativo) o `|α| < ZERO_TOL` revierte a cilíndrico — para problemas lineales monotónicos α≡0 exactamente, no es un error sino la matemática del régimen sin disipación física. En el corrector, `|α| < ZERO_TOL` aborta el paso (`CorrectionAborted`) y se biseca τ.
- **Bisección**: según el modo del paso fallido, `dl *= 0.5` o `τ *= 0.5`, con aborto (`RuntimeError`) por debajo de `ARCLENGTH_MIN_DL_FACTOR` veces `dl_ref` o `initial_tau`. Si el predictor encuentra la tangente singular, divide a la vez `dl` y `τ` y reintenta.
- **Entrypoint público**: registrado vía `@SolverRegistry.register` al cargar el módulo; en YAML, `user_modules: [dissipation_control]` y `solver: {type: DissipationArcLengthSolver, ...}`; en Python, `solidum.load_user_module("dissipation_control")` o `from solidum.user.dissipation_control import DissipationArcLengthSolver`. Sin el módulo cargado, el YAML falla con un error que sugiere declararlo y lista los módulos disponibles.

### Tests (suite [tests/user/dissipation_control/](../../../../tests/user/dissipation_control/))

[`test_dctl_dissipation_arclength.py`](../../../../tests/user/dissipation_control/test_dctl_dissipation_arclength.py) (hasta el 2026-09-29, `tests/test_dissipation_arclength.py`) cubre:

- ✅ Recuperación cilíndrica en régimen elástico puro (regresión a paridad de bits con `ArcLengthSolver`).
- ✅ Daño 1D con softening pronunciado (rama post-pico atravesada).
- ✅ Switching automático cilíndrico↔disipación observable.
- ✅ Daño 2D continuo (Quad4 + IsotropicDamage2D plane strain).
- ✅ Balance energético cualitativo (ΣΔE_d > 0 con daño activo) sobre una barra de daño continuo 1D, sin grieta ni cohesivo.

Los dos tests del `_negative_pivots` se retiraron con el método. [`test_dctl_contratos.py`](../../../../tests/user/dissipation_control/test_dctl_contratos.py) cubre el contrato de solver del barrido del programa principal, el primer paso heredado, `initial_tau`, el gancho de inicio de paso (ADR 0020, P5) y la carga desde YAML con y sin `user_modules`; [`test_dctl_specs.py`](../../../../tests/user/dissipation_control/test_dctl_specs.py), esta spec frente a la clase registrada.

### Validación parcial — qué falta para status `validated`

No están cubiertos: `snap_through_lee_frame`, los tres casos con componentes del módulo de usuario `discontinuities` (`cohesive_embedded_softening_postpeak`, `balance_energetico_GF`, `convergencia_temporal_tau`), `sign_of_pivot_tracking_bifurcacion_simple` (capacidad no implementada, §5) y `fallback_a_cilindrico_en_descarga_global`.

El acceptance `cohesive_embedded_softening_postpeak` **no está cubierto** en esta primera implementación. Diagnóstico:

La activación discreta del cohesivo en `prepare_step` (Rankine onset al cruzar `σ_t0`) introduce un salto en `F_int` y `K_t` entre el paso pre-activación y el paso post-activación. El predictor por `sign = sign(dot(δU_prev, du_t))` puede heredar una orientación inconsistente cuando el Newton del paso de activación reorienta `dU_iter` para acomodar el nuevo `K_t` con el cohesivo activo. El paso siguiente entonces predice hacia atrás (descarga) en lugar de continuar la rama post-pico.

Resolverlo limpiamente requiere **una de las dos vías siguientes**, ambas más allá del scope incremental de esta spec:

1. **Sign-of-pivot tracking con LDLᵀ Bunch-Kaufman**. Permitiría detectar la inercia exacta de `K_t` y orientar el predictor por cambio de signo de pivot en lugar del producto escalar `δU·du_t`. Esto desacopla la orientación del predictor de la cinemática perturbada por la activación. Requiere extender el backend algebraico (ADR 0003 fase 2) con un solver LDLᵀ verdadero.
2. **Control por CMOD/CTOD del salto** mediante MPC sobre el DOF interno del cohesivo, en lugar de la restricción de disipación global. Cambia la naturaleza del solver — sería una spec separada (`CMODControlSolver`) en lugar de variante de arc-length. *Nota 2026-09-29*: desde el 2026-09-23 el programa principal tiene el control indirecto de desplazamiento, [`IndirectDisplacementSolver`](../../../specs/IndirectDisplacementSolver.md) (de Borst 1987), que controla una combinación de desplazamientos nodales (p. ej. una apertura medida entre nodos), no el DOF interno del salto.

*Nota 2026-09-29*: este diagnóstico (2026-05-19) es anterior a la corrección de la tangente cohesiva del 2026-09-23 (`597f3d9`) y no se ha repetido; con la tangente corregida el comportamiento de este solver en el caso cohesivo+embebido no está comprobado.

Recomendación: dejar esta validación como deuda técnica priorizada en STATUS.md §"Deuda técnica priorizada" hasta que el LDLᵀ Bunch-Kaufman entre en otro contexto (lo amerita por más razones, no solo esta).

---

## Diálogo

- **2026-05-19 (tarde)** · Spec promovida de `draft` a `implemented`. 7/7 tests verde sobre los casos de daño continuo (1D y 2D); el caso crítico cohesivo+embedded reveló una limitación más profunda del `sign` por producto escalar tras activación discreta del Rankine, que requiere LDLᵀ Bunch-Kaufman real o control CMOD/CTOD — ambos fuera del scope incremental de esta variante. Documentado en §"Implementación" como deuda técnica. La spec mantiene su valor: provee la arquitectura del switching cilíndrico↔disipación y la fórmula lineal de Gutiérrez funcional para todos los problemas con softening continuo (daño bulk 1D y 2D), que es donde más casos prácticos de Solidum se concentran.

- **2026-05-19 (mañana)** · Spec creada como variante de `ArcLengthSolver` para desbloquear ADR 0010 fase 4 (embedded discontinuity con softening severo) y validación cuantitativa de $G_F$ en pipeline real. Motivada por el cierre de la campaña de validación externa Tandas 12-14 ([[project_validacion_externa_tandas_12_14]]) que identificó esta pieza como **el primer componente nuevo priorizado** tras la campaña. Sin ADR nuevo: reusa decisiones arquitecturales del padre (estructura del bucle predictor/corrector, manejo de Dirichlet/MPC, ADR 0003 algebraico, ADR 0007 convergencia, ADR 0008 unidades). El cambio se circunscribe a la **restricción del paso** (cuadrática → lineal) más mejoras secundarias (switching automático cilíndrico↔disipación, sign-of-pivot tracking).

- **Puntos abiertos para el usuario** (físicos/de alcance, no de plumbing — ver `[[feedback_filtrar_puntos_abiertos]]`):
  1. ¿La restricción de disipación de Gutiérrez 2004 es la fórmula que quieres, o prefieres una alternativa (p. ej. constraint sobre `K_t`-norma de Crisfield 1991 §9.4)? La de Gutiérrez es la más extendida en literatura cohesiva.
  2. ¿El switching automático cilíndrico↔disipación es deseable, o prefieres dos solvers separados (`CylindricalArcLengthSolver` y `DissipationArcLengthSolver` puros) que el usuario elige según el caso? El switching es estándar en Verhoosel et al. 2009; dos solvers separados son más predecibles pero exigen que el usuario diagnostique cuándo cambiar.
  3. ¿Sign-of-pivot tracking se aborda en esta spec o se difiere a una spec independiente (`PivotTracker`)? Está acoplado al backend algebraico (ADR 0003) — si el LDLᵀ Bunch-Kaufman se introduce como helper compartido, otros solvers se beneficiarían. Mi propuesta es **abordarlo aquí en versión "signo del determinante LU"** (aproximación válida para indefinitud simple) y dejar el LDLᵀ verdadero como deuda técnica del proyecto. *(2026-09-29: la versión por signo del determinante se implementó pero nunca se conectó al predictor, y se retiró; ver §5.)*
- **2026-09-22** · Auditoría global: hereda la reestructuración del corrector y los metadatos de trazado de `ArcLengthSolver` (ensamblaje único por iteración, estado committed coherente, `reached_max_lambda`).
- **2026-09-22** · ADR 0015: el bucle de Newton pasa al corrector compartido `NewtonCorrector`; este solver aporta su problema por paso (`_DissipationArcProblem`: sobreescribe `constraint` con la restricción lineal de Gutiérrez para el modo de disipación y hereda los modos cilíndrico y de cierre exacto; `α ≈ 0` → `CorrectionAborted`) y conserva su control de paso. Sin cambio de formulación ni de resultados (suite completa sin tocar ningún valor esperado).
- **2026-09-23** · ADR 0017-0019: hereda de `ArcLengthSolver` los valores nuevos de `linear_algebra` y la comprobación de mecanismo rígido al empezar `solve`.
- **2026-09-23** · Hereda el cambio de formulación del padre validado por el usuario: primer paso como fracción de carga (`initial_dlambda`, 0.1 por omisión; `initial_dl` explícito y excluyente) y llegada exacta a `max_lambda` dentro del tramo recorrido, que sustituye a la salvaguarda del `final_step`. Ver el Diálogo de [`ArcLengthSolver`](../../../specs/ArcLengthSolver.md). Los tests de esta variante pasan sin cambiar ningún valor esperado (todos declaran `initial_dl` explícito).
- **2026-09-29** · **Sale del programa principal al módulo de usuario `dissipation_control`** (ADR 0020), por el criterio de admisión de Reglas.md §1: ningún código general establecido implementa la restricción de Gutiérrez de forma nativa (ver la clasificación al principio). Spec, código y tests pasan a `docs/user/dissipation_control/specs/`, `solidum/user/dissipation_control/` y `tests/user/dissipation_control/`; el YAML necesita `user_modules: [dissipation_control]`. El programa principal gana una sola pieza, genérica: el nombre público `solidum.math.solvers.ArcProblem` para el problema de Newton por paso del arco, que ya extendía `IndirectDisplacementSolver`. En la misma revisión se retira `_negative_pivots()` —que ningún `solve` llamaba, en este solver ni en el padre— con sus dos tests, y la spec pasa a describir lo que hace el código: el predictor se orienta por el producto escalar con el incremento previo, sin seguimiento de inercia (§5, §8); el callback no expone la disipación acumulada (§3); el retorno a cilíndrico usa el umbral de $\Delta E_d$ y la detección de $\alpha \approx 0$, no un `tau_relax_iter_threshold` (§5). Los tres acceptance con `CST_Embedded2D` declaran que requieren el módulo `discontinuities`. Sin cambio de formulación ni de resultados.
- **2026-09-29** · Segunda revisión documental del día, contra `solidum/user/dissipation_control/dissipation_arclength.py` (sin cambiar el código ni la formulación). (a) §3 ya no promete que el callback exponga $\sum\Delta E_d$: dice que no lo hace y cómo se reconstruye. (b) §5 describe el criterio real de vuelta a cilíndrico ($\Delta E_d < \varepsilon_{\text{diss}}$ tras el commit y detección de $\alpha \approx 0$ en el predictor) y retira el de `tau_relax_iter_threshold`, que no existe; recoge también que al entrar en disipación $\tau \leftarrow \max(\Delta E_d, \texttt{initial\_tau})$, no $\Delta E_d$, y §0 que la detección es por energía, no por variables internas del material. (c) El pseudocódigo de §4 pasa a describir el `solve`: al fallar un paso $\tau$ se divide entre 2 (`τ *= 0.5`), no entre `shrink_factor`; con tangente singular en el predictor se dividen a la vez $dl$ y $\tau$; se añaden la detección de $\alpha \approx 0$, el aborto del corrector con $\alpha \approx 0$ y la llegada exacta a `max_lambda`, que no decide modo ni adapta el paso. §9 e Implementación, en consonancia. Quedan registrados en §5, **sin resolver y pendientes del usuario**, los tres puntos físicos abiertos: el signo del predictor en modo de disipación ($\text{sign}\cdot|\tau/\alpha|$ frente a $\tau/\alpha$), el piso absoluto 1 de la escala del umbral y el alcance en plasticidad. Notas nuevas en la Motivación y en la Validación parcial: los diagnósticos con penalización cohesiva rígida son anteriores a la corrección de la tangente cohesiva (`597f3d9`) y no se han repetido.
