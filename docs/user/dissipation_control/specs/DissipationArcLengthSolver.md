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

Solver estático no lineal con **traza por disipación**: la magnitud que se mantiene constante de paso a paso ya no es la longitud euclídea $\lVert\Delta\mathbf U\rVert$ sino la **disipación incremental de energía** $\Delta E_d$ asociada al paso. Pensado para:

- Fractura cohesiva (Mode-I y, en el futuro, mixto) con penalty `K_e` stiff.
- Daño con softening severo donde la rama post-pico tiene pendiente casi vertical en el plano $(U, \lambda)$.

**Alcance (2026-09-29): daño con descarga secante** —continuo o cohesivo—, es decir, materiales cuya fuerza interna es $\mathbf F_{\text{int}} = \mathbf K_s(\text{daño})\,\mathbf U$ y que descargan hacia el origen. Es la hipótesis de la que sale la $\Delta E_d$ de §1, y el código no tiene otra. **Plasticidad, fuera del alcance**: con deformación plástica la descarga no es secante y esa expresión no mide la disipación. La variante para plasticidad de Verhoosel et al. (2009), que descuenta la energía de la descarga elástica, queda como ampliación futura.

**Estado de uso medido (2026-09-30)**: desde el cambio de formulación del 2026-09-30 (§5: tangente de carga en el predictor y arranque del modo de disipación por $\alpha_L$) la rama de ablandamiento se resuelve de verdad con la restricción de disipación, y el solver sigue el retroceso exacto de la barra con daño localizado (acceptance `retroceso_barra_dano_localizado`: 172 puntos post-pico, distancia normalizada máxima $3.2\cdot10^{-4}$, $u_{\min}/u_{\text{pico}} = 0.616376$). Hasta el 2026-09-29 el predictor del modo de disipación resultaba degenerado ($\alpha \approx 0$ con la tangente secante del estado consolidado) y el solver se comportaba en la práctica como el arco cilíndrico del padre (deuda #27 de STATUS.md).

Régimen elástico (sin disipación): la restricción degenera a $g = 0$. El solver detecta este caso y **conmuta automáticamente** a arc-length cilíndrico durante el régimen elástico, y pasa a control por disipación cuando la **tangente de carga** del último paso convergido se separa de la secante (§5). La detección es global, sobre matrices y vectores ensamblados: el solver no consulta las variables internas del material (daño, apertura de la grieta).

### 1. Ecuación de equilibrio resuelta

Idéntica al padre:

$$\mathbf F_{\text{int}}(\mathbf U) = \lambda\,\mathbf F_{\text{ext}}^{\text{ref}}$$

con $\lambda$ incógnita. La **restricción del paso** cambia respecto al padre.

#### Restricción de disipación (Gutiérrez 2004, Verhoosel et al. 2009)

Hipótesis: daño con descarga secante, $\mathbf F_{\text{int}} = \mathbf K_s\,\mathbf U$ con $\mathbf K_s$ función sólo del daño. La energía elástica almacenada es $\tfrac12\mathbf U^T\mathbf K_s\mathbf U$ y la tasa de disipación es la liberación de energía por degradación de $\mathbf K_s$:

$$\dot E_d = -\tfrac12\,\mathbf U^T\dot{\mathbf K}_s\,\mathbf U = \tfrac12\left(\lambda\,\mathbf F_{\text{ext}}^{\text{ref}T}\,\dot{\mathbf U} - \dot\lambda\,\mathbf F_{\text{ext}}^{\text{ref}T}\,\mathbf U\right),$$

donde la segunda igualdad sale de derivar el equilibrio $\mathbf K_s\mathbf U = \lambda\,\mathbf F_{\text{ext}}^{\text{ref}}$. Integrada en el paso con el estado inicial (Gutiérrez 2004), la energía disipada por paso es:

$$\Delta E_d = \frac{1}{2}\left(\lambda_n\,\mathbf F_{\text{ext}}^{\text{ref}T}\,\Delta\mathbf U \;-\; \Delta\lambda\,\mathbf F_{\text{ext}}^{\text{ref}T}\,\mathbf U_n\right)$$

donde $(\mathbf U_n, \lambda_n)$ son los valores committed al inicio del paso y $(\Delta\mathbf U, \Delta\lambda)$ los incrementos del paso. Con deformación plástica la descarga no es secante y esta expresión no mide la disipación (§0). La restricción es:

$$g(\Delta\mathbf U, \Delta\lambda) \;\equiv\; \frac{1}{2}\left(\lambda_n\,\mathbf F_{\text{ext}}^{\text{ref}T}\,\Delta\mathbf U \;-\; \mathbf U_n^T\,\mathbf F_{\text{ext}}^{\text{ref}}\,\Delta\lambda\right) \;=\; \tau$$

con $\tau$ el **incremento de disipación objetivo** (parámetro del solver, papel análogo al `initial_dl` del padre).

**Linealidad**: $g$ es lineal en $(\Delta\mathbf U, \Delta\lambda)$. Por eso el corrector tiene **una sola raíz**, no dos como en la cuadrática cilíndrica — eliminando la selección por menor ángulo y la patología de raíces imaginarias en pasos grandes.

### 2. Condiciones de contorno

Idénticas al padre (`ArcLengthSolver` §2). Dirichlet homogéneo/no homogéneo y MPC vía ADR 0004.

### 3. Salidas físicas

Idénticas al padre: `SolveResult` con $\mathbf U$ final y `step_callback(step, U, λ)` para la historia, con la misma firma que en el padre. El solver **no expone** la disipación del paso ni la acumulada $\sum_n \Delta E_d^{(n)}$: la calcula tras cada paso convergido para decidir la vuelta a cilíndrico (§5) y sólo la escribe en el registro (`ΔE_d=…` en el log del solver). Quien la necesite —por ejemplo, para contrastar **$G_F$ de extremo a extremo** con $\int_{\Gamma_d} t\cdot d[\![u]\!] \approx \sum_n \Delta E_d^{(n)}$— la reconstruye desde la historia $(\mathbf U_n, \lambda_n)$ del callback con la fórmula de $\Delta E_d$ de §1, como hace el test 5 de [`test_dctl_dissipation_arclength.py`](../../../../tests/user/dissipation_control/test_dctl_dissipation_arclength.py).

---

## Formulación numérica — solo lo que cambia respecto al padre

### 4. Esquema operativo

Estructura idéntica al padre con cuatro diferencias: el predictor y la restricción del corrector dependen del modo, el arranque del modo de disipación se decide en el predictor con la tangente de carga, y tras cada paso convergido se decide la vuelta a cilíndrico (§5). El esquema describe el `solve` del código (`dissipation_arclength.py`):

```
λ = 0;  U = 0;  ΔU_prev = 0
τ = initial_tau          # análogo a initial_dl, pero en unidades de energía
modo = "cylindrical"      # siempre arranca así: con λ_n = 0, U_n = 0 es α = 0
dl = (sin fijar)
K_L = (ninguna)           # tangente de carga: la del último iterado convergido
mientras λ < max_lambda y step < max_steps:
    step += 1
    prepare_all_steps(U)                  # gancho de inicio de paso (ADR 0020, P5)

    # Predictor tangente (el del padre), ensamblado en el estado consolidado U_n
    δu_t = K_t(U_n)^{-1} · F_ext^ref      # con daño y descarga secante, K_t(U_n) = K_s
    si K_t singular:                      # el predictor no puede resolver
        si dl fijado: dl /= 2
        τ /= 2                            # se dividen los dos, sin comprobar cotas
        continuar                         # reintento del paso
    si dl sin fijar: dl = dl_ref = Δl₁     # primer paso del padre (initial_dlambda / initial_dl)
    sign = sign(ΔU_prev · δu_t)           # regla del padre (sólo modo cilíndrico)

    # Tangente de carga (desde 2026-09-30)
    δu_L = K_L^{-1} · F_ext^ref           # ninguna si aún no hay paso convergido o K_L es singular

    # (b) Arranque del modo de disipación por α_L (desde 2026-09-30)
    si modo == "cylindrical" y δu_L existe:
        α_L = ½·(λ_n·F_ext^ref·δu_L − F_ext^ref·U_n)
        si |α_L| > 1e-3 · |½·λ_n·F_ext^ref·δu_L|:   # la tangente de carga se separa de la secante
            modo = "dissipation";  τ = initial_tau

    modo_paso = modo
    δu = δu_t
    si modo_paso == "dissipation":
        δu = δu_L si existe                # (a) predictor sobre la tangente de carga
        α = ½·(λ_n·F_ext^ref·δu − F_ext^ref·U_n)
        escala = |½·λ_n·F_ext^ref·δu|
        si (escala > 0 y |α| < 1e-6·escala) o α_nulo(α):
            modo_paso = modo = "cylindrical"   # α no significativo (§5)
        # α_nulo(α): |α| ≤ ZERO_TOL·½·(|λ_n·F·δu| + |F·U_n|), relativa (2026-09-30)
    si modo_paso == "cylindrical":
        Δλ_pred = sign · dl / (‖δu_t‖ + ZERO_TOL)
        δu = δu_t
    sino:
        Δλ_pred = τ/α                     # la restricción fija valor y signo (desde 2026-09-29)
    ΔU_pred = Δλ_pred · δu

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
            si α_nulo(den): abortar el paso                   # CorrectionAborted → bisección
            ddλ = num / den
        ΔU += δu_R + ddλ·δu_t
        Δλ += ddλ
        si converge: break

    si converge y λ + Δλ > max_lambda:    # llegada exacta del padre (ArcLengthSolver._land)
        repetir desde (U_n, λ_n) con λ = max_lambda fijo
        si converge: consolidar y continuar     # sin decidir modo ni adaptar dl/τ
        sino: tratar como paso fallido
    si converge:
        K_L = K_t del ensamblaje en el iterado convergido   # prueba de carga respecto del estado anterior
        ΔE_d = ½·(λ_n·F_ext^ref·ΔU − Δλ·F_ext^ref·U_n)       # paso completo, tras corregir
        umbral = dissipation_threshold · ‖F_ext^ref‖·‖U_n+1‖    # relativo, sin piso (desde 2026-09-29)
        si modo == "dissipation" y ΔE_d < umbral:
            modo = "cylindrical"          # el arranque ya no se decide aquí (desde 2026-09-30)
        adaptar dl (modo cilíndrico) o τ (modo disipación) según el modo del SIGUIENTE paso (§9)
        consolidar: U += ΔU;  λ += Δλ;  ΔU_prev = ΔU
    sino:                                 # bisección según el modo del paso fallido
        si modo_paso == "cylindrical":
            dl *= 0.5;  si dl < ARCLENGTH_MIN_DL_FACTOR·dl_ref: abortar (RuntimeError)
        sino:
            τ *= 0.5;   si τ < ARCLENGTH_MIN_DL_FACTOR·initial_tau: abortar (RuntimeError)
```

El paso fallido reduce τ a la mitad (`τ *= 0.5`); `tau_shrink_factor` sólo actúa sobre pasos que convergen lentamente (§9). La tangente de carga se guarda sólo tras un paso con restricción que converge: un paso fallido conserva la del último paso convergido, y la llegada a `max_lambda` termina el trazado.

Coste: una resolución lineal más por paso (la de $\delta\mathbf u_L$, con el mismo backend y la misma reducción de Dirichlet/MPC que el predictor del padre) y una matriz tangente guardada entre pasos. El predictor del padre, en el estado consolidado, se mantiene: calibra el criterio de convergencia en el primer paso, detecta la tangente singular y da la dirección del modo cilíndrico.

### 5. Predictor / corrector — punto clave del cambio

#### Switching automático cilíndrico ↔ disipación

El solver mantiene un atributo `modo` con dos estados:

- `"cylindrical"`: arc-length cilíndrico estándar; restricción $\lVert\Delta\mathbf U\rVert^2 = dl^2$. **Régimen elástico**, sin disipación activa.
- `"dissipation"`: arc-length por energía; restricción $g = \tau$. **Régimen disipativo**.

#### Tangente de carga (formulación adoptada, decisión del usuario del 2026-09-30)

En el predictor de un paso intervienen dos tangentes:

- $\mathbf K_s$: la que ensambla el predictor del padre en el estado consolidado $\mathbf U_n$. Ahí los materiales de daño con descarga secante no están en carga activa ($\varepsilon_{eq} = \kappa$, sin incremento) y devuelven la rigidez secante. Como $\mathbf F_{\text{int}}(\mathbf U_n) = \mathbf K_s\,\mathbf U_n = \lambda_n\mathbf F + \mathbf r_n$, con $\delta\mathbf u_t = \mathbf K_s^{-1}\mathbf F$ el denominador de la restricción es sólo el residuo de equilibrio, $\alpha = -\tfrac12\,\mathbf F^T\mathbf K_s^{-1}\mathbf r_n$ (§ «Detección de $\alpha \approx 0$»).
- $\mathbf K_L$, **tangente de carga**: la del **último iterado convergido del paso anterior**, la que tenía el Newton al converger. Está ensamblada en $\mathbf U_n$ con la prueba de carga de Kuhn-Tucker respecto del estado consolidado anterior ($\kappa_{n-1}$), así que los elementos que dañaban en ese paso aportan su tangente consistente. El solver la toma del corrector compartido: es la matriz del `(K, F_int)` del ensamblaje en el iterado convergido (`CorrectorResult.state`), sin ensamblaje adicional.

Con $\delta\mathbf u_L = \mathbf K_L^{-1}\mathbf F$ y $\mathbf U_n = \mathbf K_s^{-1}(\lambda_n\mathbf F + \mathbf r_n)$:

$$\alpha_L = \tfrac12\left(\lambda_n\,\mathbf F\cdot\delta\mathbf u_L - \mathbf F\cdot\mathbf U_n\right) = \tfrac12\,\lambda_n\,\mathbf F^T\left(\mathbf K_L^{-1} - \mathbf K_s^{-1}\right)\mathbf F \;-\; \tfrac12\,\mathbf F^T\mathbf K_s^{-1}\mathbf r_n .$$

El primer término mide el crecimiento del daño (se anula si $\mathbf K_L = \mathbf K_s$); el segundo es el ruido del equilibrio. En la barra con retroceso, en el paso siguiente al inicio del daño, $\alpha_L \approx -0.42$ J frente a una energía elástica $\tfrac12\mathbf F\cdot\mathbf U_n \approx 1.35$ J; es negativo, como corresponde a un retroceso.

**(a) Predictor del modo de disipación**: $\Delta\mathbf U = \Delta\lambda\,\delta\mathbf u_L$ y $\Delta\lambda_\text{pred} = \tau/\alpha_L$ (siguiente apartado). Si $\mathbf K_L$ no existe aún o es singular, el paso usa $\delta\mathbf u_t$ como antes del 2026-09-30, con lo que $\alpha \approx 0$ lo devuelve a cilíndrico.

**(b) Arranque del modo de disipación**: en modo `"cylindrical"`, antes de cada paso con $\mathbf K_L$ disponible,

$$|\alpha_L| > 10^{-3}\cdot\tfrac12\left|\lambda_n\,\mathbf F\cdot\delta\mathbf u_L\right| \;\Rightarrow\; \text{modo} \leftarrow \texttt{"dissipation"},\quad \tau \leftarrow \texttt{initial\_tau},$$

es decir, cuando la tangente de carga se separa de la secante (constante `ALPHA_SWITCH_RTOL` del módulo). Umbral: en régimen sin disipación $\mathbf K_L = \mathbf K_s$ y $|\alpha_L|/\text{escala} = |\mathbf F^T\mathbf K_s^{-1}\mathbf r_n| / |\lambda_n\,\mathbf F^T\mathbf K_s^{-1}\mathbf F| \le (\lVert\mathbf r_n\rVert/\lVert\lambda_n\mathbf F\rVert)/\cos\theta$, con $\theta$ el ángulo entre $\mathbf F$ y $\mathbf K_s^{-1}\mathbf F$: del orden de la tolerancia relativa de fuerza del criterio de convergencia ($10^{-5}$ por omisión). Medido en la barra con retroceso: $\approx 10^{-15}$ en la rama elástica y $\approx 10^{-6}$ con la secante tras el inicio del daño. $10^{-3}$ deja dos órdenes de margen, y con daño en crecimiento la razón es $O(0.1\text{–}1)$ (0.45 en el paso que conmuta en la barra). Si el usuario afloja el criterio de convergencia hasta $\text{rtol} \gtrsim 10^{-4}$, el margen desaparece. Como $\mathbf K_L$ es la del paso anterior, el cambio de modo llega en el paso **siguiente** al inicio del daño: si el daño empieza en el paso convergido $N$, los pasos $1..N$ son cilíndricos y los siguientes, de disipación (acceptance `switching_cilindrico_a_disipacion_en_pico`).

Sustituye al criterio anterior, $\Delta E_d > \varepsilon_{\text{diss}}$ tras un paso cilíndrico, que en un pico anguloso no se alcanza: en la barra con daño localizado el arco cilíndrico se atascaba en la esquina con pasos de $\sim 10^{-10}$ m cuya disipación ($\le 1.1\cdot10^{-12}$ J) quedaba bajo el umbral ($\approx 5.3\cdot10^{-12}$ J). Con $\alpha_L$, un paso con daño en crecimiento, por pequeño que sea, basta: la tangente del elemento débil salta de $E$ a $-\alpha_d\,\kappa_{0w}\,E = -0.475\,E$ ($\alpha_d$, parámetro de ablandamiento del material). Al entrar, $\tau$ vuelve a `initial_tau` (antes, $\max(\Delta E_d, \texttt{initial\_tau})$ con la $\Delta E_d$ del paso que conmutaba).

**Alcance del criterio**: $\alpha_L \neq 0$ también sin disipación cuando la tangente no es la secante por otra causa: con no linealidad geométrica ($\mathbf K_L = \mathbf K_t(\mathbf U_n)$ y $\mathbf U_n \neq \lambda_n\mathbf K_t^{-1}\mathbf F$) y con plasticidad (tangente elastoplástica frente a descarga elástica). El criterio sólo vale dentro del alcance de §0: **pequeñas deformaciones y descarga secante**. No se ha podido confirmar en el texto de Verhoosel et al. (2009) qué tangente usa su predictor; queda como comprobación bibliográfica pendiente.

#### Vuelta a cilíndrico (sin cambios el 2026-09-30)

- **En el predictor** de un paso en `"dissipation"`: $|\alpha| < 10^{-6}\cdot$escala (constante `ALPHA_LINEAR_RTOL`) o $\alpha$ nulo a precisión de máquina (guarda relativa, abajo) ⇒ el paso se hace en cilíndrico y el modo pasa a `"cylindrical"` (apartado «Detección de $\alpha \approx 0$»). Con la tangente de carga ocurre cuando $\mathbf K_L$ vuelve a ser la secante: el último paso no hizo crecer el daño (descarga) o el daño saturó en `DAMAGE_MAX`, donde el material devuelve la secante. Con el umbral de arranque forma una **histéresis** ($10^{-3}$ para entrar, $10^{-6}$ para salir) que evita que el modo oscile entre pasos contiguos.
- **Tras cada paso convergido** en `"dissipation"` que no es la llegada a `max_lambda`: $\Delta E_d < \varepsilon_{\text{diss}}$ (disipación nula o negativa) ⇒ el siguiente paso vuelve a `"cylindrical"`, con

$$\varepsilon_{\text{diss}} = \texttt{dissipation\_threshold}\cdot\lVert\mathbf F_{\text{ext}}^{\text{ref}}\rVert\,\lVert\mathbf U_{n+1}\rVert .$$

Umbral **relativo** a la escala energética del estado, sin piso absoluto (desde el 2026-09-29). Hasta entonces era $\texttt{dissipation\_threshold}\cdot\max(\lVert\mathbf F\rVert\,\lVert\mathbf U\rVert,\ 1)$: ese 1 está en unidades de energía del modelo y hacía que el umbral dependiera del sistema de unidades cuando $\lVert\mathbf F\rVert\,\lVert\mathbf U\rVert < 1$ (ADR 0006). Un paso de disipación converge con $\Delta E_d = \tau$, así que este criterio sólo actúa si $\tau$ cae por debajo de $\varepsilon_{\text{diss}}$; en la práctica la vuelta la decide el predictor.

#### Predictor en modo de disipación (desde el 2026-09-29)

Con el predictor tangente $\Delta\mathbf U = \Delta\lambda\,\delta\mathbf u$ —desde el 2026-09-30, $\delta\mathbf u = \delta\mathbf u_L = \mathbf K_L^{-1}\mathbf F_{\text{ext}}^{\text{ref}}$, la tangente de carga— la restricción lineal da $g = \alpha\,\Delta\lambda$ con

$$\alpha = \tfrac12\left(\lambda_n\,\mathbf F_{\text{ext}}^{\text{ref}}\cdot\delta\mathbf u - \mathbf F_{\text{ext}}^{\text{ref}}\cdot\mathbf U_n\right),$$

así que exigir $g = \tau > 0$ (disipación positiva) fija el valor **y el signo** del incremento de carga, sin regla de signo (Gutiérrez 2004; Verhoosel et al. 2009):

$$\Delta\lambda_\text{pred} = \tau/\alpha .$$

Con $\alpha < 0$ el paso reduce la carga, que es lo que exige un retroceso (*snap-back*). Hasta el 2026-09-29 el código usaba $\text{sign}\cdot|\tau/\alpha|$ con el signo del producto escalar del padre: cuando $\text{sign} \neq \operatorname{sgn}(\alpha)$ el predictor arrancaba en $g = -\tau$ (disipación negativa) y el corrector tenía que darle la vuelta. En el prototipo de diagnóstico del 2026-09-29 (tangente de carga en el predictor, arriba) el signo del padre discrepaba del de $\alpha$ en 57 de 141 predictores del retroceso, y con la regla antigua el trazado se perdía.

#### Detección de $\alpha \approx 0$ en el predictor

**En el predictor** de un paso en `"dissipation"`, si $|\alpha| < 10^{-6}\cdot\tfrac12|\lambda_n\,\mathbf F_{\text{ext}}^{\text{ref}}\cdot\delta\mathbf u|$, o $\alpha$ es nulo a precisión de máquina, el paso se hace en cilíndrico y el modo pasa a `"cylindrical"`. Tras ese paso vuelve a aplicarse el criterio de arranque por $\alpha_L$, que puede devolver el modo a `"dissipation"`. Dentro del corrector, $\alpha$ nulo con el $\delta\mathbf u_t$ de la iteración (la tangente del iterado, en carga) aborta el paso (`CorrectionAborted`), que se biseca.

**Guarda de $\alpha$ nulo, relativa** (revisión adversaria del 2026-09-30): $|\alpha| \le \texttt{ZERO\_TOL}\cdot\tfrac12\big(|\lambda_n\,\mathbf F\cdot\delta\mathbf u| + |\mathbf F\cdot\mathbf U_n|\big)$, es decir, relativa a los dos términos de los que $\alpha$ es diferencia; salta cuando ambos son nulos ($\lambda_n = 0$, $\mathbf U_n = \mathbf 0$). Hasta entonces se comparaba $|\alpha|$ con `ZERO_TOL` $= 10^{-12}$ en unidades de energía: con energías del orden de $10^{-12}$ (una probeta de micras en SI; medido con la barra de validación en unidades de longitud y fuerza $10^6$ veces mayores) el $\alpha$ de un paso de ablandamiento genuino ($\approx 4\cdot10^{-13}$) caía por debajo, el paso volvía a cilíndrico en la esquina del pico y el solver abortaba en el paso 24. Con la guarda relativa el trazado coincide con el de SI (1.5e-10); lo fija `test_el_trazado_no_depende_de_las_unidades`. Queda una operación no covariante en el camino, el `+ ZERO_TOL` del $\Delta\lambda$ del predictor cilíndrico, heredado de la fórmula de `ArcLengthSolver`: sólo mueve el punto de partida de Newton (un 1 % con $\lVert\mathbf K^{-1}\mathbf F\rVert \approx 10^{-10}$).

$\alpha$ es nulo en estas situaciones, y la detección no las distingue:

- **Problema lineal monótono**: $\mathbf U_n = \lambda_n\mathbf K^{-1}\mathbf F$ ⇒ $\alpha \equiv 0$. Es la matemática del régimen sin disipación.
- **Tangente de carga igual a la secante**: el último paso no hizo crecer el daño (régimen elástico, descarga) o el daño saturó en `DAMAGE_MAX`. Entonces $\alpha_L = -\tfrac12\,\mathbf F^T\mathbf K_s^{-1}\mathbf r_n$, sólo el residuo. Es físicamente correcto: moverse sobre la secante (descarga y recarga elásticas) no disipa, y la restricción no puede fijar $\Delta\lambda$ en esa dirección. En la barra con retroceso, al saturar el daño del elemento débil el solver vuelve así a cilíndrico.
- **Tangente secante en el predictor** (medido 2026-09-29, **resuelto el 2026-09-30** con la tangente de carga). Con $\delta\mathbf u_t = \mathbf K_s^{-1}\mathbf F$ del estado consolidado resultaba
  $$\alpha = \tfrac12\left(\lambda_n\,\mathbf F^T\mathbf K_s^{-1}\mathbf F - \mathbf F^T\mathbf U_n\right) = -\tfrac12\,\mathbf F^T\mathbf K_s^{-1}\mathbf r_n$$
  incluso con el daño creciendo, y la detección lo tomaba por régimen lineal. Medido entonces en los tests del módulo: con `IsotropicDamage1D` (barra de un elemento) y con `IsotropicDamage2D` (Quad4), ningún problema de Newton se construía en modo de disipación (500 de 500 y 300 de 300 en cilíndrico), con $|\alpha|/\text{escala} \le 7\cdot10^{-14}$; el modo `"dissipation"` sólo se veía en el callback, entre la conmutación tras el commit y el predictor siguiente, que lo revertía. En la barra con retroceso, $|\alpha|/\text{escala} \approx 10^{-6}$, el nivel de la tolerancia de equilibrio. Hoy sólo ocurre si $\mathbf K_L$ es singular y el paso recurre a $\delta\mathbf u_t$.

No hay criterio por número de iteraciones: la spec describía hasta el 2026-09-29 una vuelta a cilíndrico con un `tau_relax_iter_threshold` que el código no tiene.

#### Puntos físicos resueltos (2026-09-29 y 2026-09-30)

Los tres registrados al trasladar el solver al módulo quedaron resueltos el 2026-09-29 con visto bueno del usuario: predictor $\Delta\lambda = \tau/\alpha$ (arriba), umbral relativo sin piso absoluto (arriba) y alcance restringido a daño con descarga secante (§0). Los dos que dejó abiertos el test de retroceso (`test_dctl_retroceso.py`) quedaron resueltos el **2026-09-30 por decisión del usuario**, con la formulación de «Tangente de carga» (arriba):

1. **Tangente del predictor en modo de disipación** → (a), la tangente de carga $\mathbf K_L$ del último iterado convergido. Se descartó la alternativa del predictor secante por extrapolación del paso anterior ($\Delta\mathbf U^0 = s\,\Delta\mathbf U_{n-1}$, $\Delta\lambda^0 = s\,\Delta\lambda_{n-1}$, con $s = \tau/g_n(\Delta\mathbf U_{n-1}, \Delta\lambda_{n-1})$): no sirve en el primer paso tras el pico, porque el paso anterior está todavía sobre la rama elástica y da otra vez $\alpha \approx 10^{-6}\cdot$escala.
2. **Arranque del modo de disipación en un pico anguloso** → (b), conmutación por $|\alpha_L| > 10^{-3}\cdot$escala. Ningún umbral de $\Delta E_d$ por paso detecta la esquina con fiabilidad, porque los pasos del cilíndrico atascado apenas disipan ($\Delta E_d/\Delta W_{\text{ext}} \approx 4\cdot10^{-7}$ en la barra).

Resultado, igual al del prototipo de diagnóstico: la barra con daño localizado sigue el retroceso con las tolerancias del test de `IndirectDisplacementSolver` (172 puntos post-pico, distancia máxima normalizada $3.2\cdot10^{-4}$, $u_{\min}/u_{\text{pico}} = 0.616376$, el exacto; sólo el elemento débil dañado), con el primer paso del caso de validación ($\Delta\lambda_1 = 0.02$). Límites que quedan (§10): con un primer paso grande el paso cilíndrico que cruza el pico salta a otra rama; el criterio (b) sólo vale dentro del alcance de §0.

#### Orientación del predictor (sin seguimiento de inercia)

En modo cilíndrico el sentido del predictor sale de la regla del padre: el signo del producto escalar del predictor tangente con el incremento del paso anterior, $\text{sign}(\Delta\mathbf U_{\text{prev}}\cdot\delta\mathbf u_t)$ (Crisfield 1991, cap. 9). En modo de disipación lo fija la restricción ($\Delta\lambda = \tau/\alpha$, desde el 2026-09-29; con $\alpha = \alpha_L$ de la tangente de carga desde el 2026-09-30). El solver **no sigue la inercia de $\mathbf K_t$** (número de pivotes negativos), y ninguna decisión del paso depende de ella.

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
- Al pasar de cilíndrico a disipación (en el predictor, por $\alpha_L$), $\tau \leftarrow \texttt{initial\_tau}$ (§5). Hasta el 2026-09-29, cuando la conmutación se decidía tras el commit, era $\max(\Delta E_d, \texttt{initial\_tau})$.

Las iteraciones que se comparan con los umbrales son las resoluciones del corrector en el paso (`n_solves`). La adaptatividad se aplica según el modo del **siguiente** paso, ya decidido: en `"cylindrical"` (régimen elástico inicial o vuelta a cilíndrico) opera sobre $dl$, idéntica al padre; en `"dissipation"`, sobre $\tau$.

### 10. Caveats numéricos

Heredados del padre (`ArcLengthSolver` §10) más los específicos:

- **Régimen disipativo no estrictamente monótono**: si el material vuelve transitoriamente a régimen elástico (descarga local), $\Delta E_d$ puede ser negativo. El solver tolera este caso pero conmuta a `"cylindrical"` si la disipación neta del paso cae por debajo del umbral, o en el predictor si la tangente de carga vuelve a ser la secante ($\alpha \approx 0$).
- **Primer paso grande: salto de rama** (medido 2026-09-30). El primer paso con daño es siempre cilíndrico: el arranque del modo de disipación necesita un paso convergido con daño en crecimiento (§5). Si ese paso cilíndrico es grande, puede converger a otro equilibrio. En la barra con daño localizado, con $\Delta\lambda_1 = 0.1$ los pasos crecen (×1.5 hasta $5\,\Delta l_1$) y el que cruza el pico lleva el predictor elástico, de deformación uniforme, por encima del umbral de los diez elementos a la vez (el margen entre $\kappa_{0w}$ y $\kappa_0$ es del 5 %): Newton converge a un equilibrio con los diez elementos ablandándose en serie ($u/u_{\text{pico}} = 1.313$, $F/F_{\text{pico}} = 0.923$), y desde ahí el solver sigue esa rama, ya en modo de disipación, sin retroceso. No es un defecto del método sino la **no unicidad del problema incremental con ablandamiento**: a una misma fuerza cada elemento puede cargar o descargar, y el arco no tiene cómo elegir la rama localizada. Remedio: un primer paso lo bastante pequeño para que el paso que cruza el pico no rebase el umbral de los elementos que deben descargarse (con $\Delta\lambda_1 = 0.02$ la barra sigue el retroceso). El test `test_limite_conocido_con_un_primer_paso_grande` lo deja fijado.
- **Fin de la rama de ablandamiento (saturación del daño)** (medido 2026-09-30). $g$ es la regla del trapecio del trabajo externo menos el cambio de energía elástica, exacta sólo si el camino del paso es lineal en $(\mathbf U, \lambda)$. El paso de disipación que cruza la saturación (`DAMAGE_MAX`) tiene por delante menos disipación que $\tau$ y cumple $g = \tau$ acreditando trabajo a lo largo de la rama de rigidez residual $(1 - \texttt{DAMAGE\_MAX})\,E$: en la barra salta de $\varepsilon_w \approx 10.4\,\kappa_{0w}$ a $\approx 46.5\,\kappa_{0w}$, sobre la curva exacta. Después la tangente de carga es la secante y el solver vuelve a cilíndrico. El punto es un equilibrio, pero el tramo intermedio no se traza. En la franja 2D de 10×2 Quad4 con retroceso (revisión adversaria del 2026-09-30, fuera de la suite) ese paso acredita $7.8\cdot10^{-3}$ J más de lo que disipa el material (0.57 % de la disipación total).
- **$\sum\Delta E_d$ como medida de la energía disipada** (medido 2026-09-30). En los pasos de disipación $\Delta E_d = \tau$ coincide con la disipación del material hasta el error del trapecio (1.4e-6 relativo en una franja 2D de 6×2 Quad4 sin retroceso). Los pasos cilíndricos, en cambio, no controlan la disipación: el que cruza el pico, crecido hasta $5\,\Delta l_1$, recorta con su cuerda la esquina del pico, y en esa franja $\sum\Delta E_d$ queda $5.2\cdot10^{-3}$ J por debajo (3.7 % de lo disipado en ese paso, 0.18 % al final del trazado). Sumado al paso de saturación (arriba), $\sum\Delta E_d$ de todo el trazado es una estimación, no un balance exacto; la disipación exacta es la del material, $\sum_{gp} V_{gp}\,D(\kappa_{gp})$.
- **Inicialización de $\tau$**: debe ser del orden de la disipación esperada por paso, es decir, una fracción de la disipación total esperada en la rama post-pico (en fractura cohesiva, típicamente una fracción pequeña de $G_F \cdot l_d$ por elemento agrietado). Demasiado pequeño ⇒ pasos numerosos; demasiado grande ⇒ saltos sobre la rama post-pico.
- **Criterio de convergencia laxo**: el umbral de arranque $10^{-3}$ de §5 supone que el ruido de equilibrio en $\alpha_L$ (del orden de `rtol_force`) queda dos órdenes por debajo. Con `rtol_force` $\gtrsim 10^{-4}$ el margen desaparece y el modo de disipación podría arrancar sin daño en crecimiento.
- **Disipación de origen no físico**: si el modelo numérico introduce disipación espuria (algorithmic damping en transitorio, regularización no consistente), $\Delta E_d$ contiene esa componente. **Este solver es estático**: si el usuario invoca `DissipationArcLengthSolver` sobre un modelo con material no disipativo (elástico puro, pequeñas deformaciones), el modo permanece en `"cylindrical"` indefinidamente — comportamiento equivalente al padre. No es un caveat operacional, es la definición del switching. Con no linealidad geométrica, en cambio, $\alpha_L \neq 0$ sin disipación y el modo de disipación arrancaría (fuera del alcance, §0 y §5).
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
      desc: "Umbral relativo: tras un paso en modo de disipación, ΔE_d < dissipation_threshold·‖F_ref‖·‖U‖ se considera 0 y devuelve el modo cilíndrico. Sin piso absoluto desde 2026-09-29. Desde 2026-09-30 no decide el arranque del modo de disipación (criterio por α_L con la tangente de carga, §5)" }

requirements:
  - "Modelo con material de daño con descarga secante (daño continuo o cohesivo) — en modelo puramente elástico el solver opera equivalente a ArcLengthSolver cilíndrico"
  - "F_ext^ref no nulo (la restricción de disipación se evalúa sobre F_ext^ref·U)"
  - "Heredados del padre: apoyos suficientes; sin modos rígidos; dl_inicial razonable"
  - "Módulo de usuario cargado: user_modules: [dissipation_control] en el YAML o solidum.load_user_module('dissipation_control')"

conventions:
  units: "heredadas del modelo (ADR 0008); τ en unidades de energía consistentes"
  stability: "incondicional respecto a puntos límite y rama post-pico de softening severo. Ante bifurcaciones, la misma que el padre: el sentido del predictor sale del producto escalar con el incremento previo (sin seguimiento de inercia, §5)"

out_of_scope:
  - "Plasticidad — la ΔE_d de §1 supone descarga secante; la variante para plasticidad de Verhoosel et al. (2009) queda como ampliación futura"
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

    - name: retroceso_barra_dano_localizado             # validado 2026-09-30 (tangente de carga, §5); xfail(strict) del 2026-09-29 al 2026-09-30
      setup: "Barra de 10 Truss2D con IsotropicDamage1D, elemento central con κ₀ un 5 % menor (caso de IndirectDisplacementSolver); initial_dlambda 0.02, initial_tau 1e-3"
      expect: "Puntos post-pico sobre la curva exacta (distancia normalizada < 1e-3, más de 100 puntos) y desplazamiento del extremo que retrocede hasta el mínimo exacto 0.616·u_pico (±2e-3); sólo el elemento débil dañado; la rama de ablandamiento, resuelta en modo de disipación. Medido: 172 puntos, distancia máxima 3.2e-4, u_min/u_pico = 0.616376. Con initial_dlambda 0.1 el paso cilíndrico que cruza el pico salta a la rama con los diez elementos ablandándose (límite conocido, §10)"
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
    - name: switching_cilindrico_a_disipacion_en_pico     # cubierto 2026-09-30 con el modo de los problemas de Newton
      setup: "Modelo con daño que empieza en el paso convergido N"
      expect: "Problemas de Newton 'cylindrical' en pasos 1..N (el N, cilíndrico por construcción: el arranque necesita la tangente de carga de un paso con daño en crecimiento); 'dissipation' desde N+1 en adelante mientras el daño crece; sin oscilación de modo en pasos contiguos"

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
- **Atributo `self._mode`**: `"cylindrical" | "dissipation"`. Arranque del modo de disipación en el predictor, por `|α_L| > ALPHA_SWITCH_RTOL·½·|λ_n·F·du_L|` con `ALPHA_SWITCH_RTOL = 1e-3` (§5, desde el 2026-09-30; antes, tras el commit por `ΔE_d > dissipation_threshold·‖F_ref‖·‖U‖`). Vuelta a cilíndrico en el predictor por α≈0 o tras el commit por `ΔE_d` bajo el umbral relativo (§5; sin el piso `max(·, 1)` desde el 2026-09-29).
- **Tangente de carga** (desde el 2026-09-30): atributo `self._K_load`, la matriz del `(K, F_int)` que devuelve el corrector compartido en `CorrectorResult.state` tras cada paso con restricción que converge; `_load_tangent_direction` resuelve `du_L = K_L⁻¹·F_ref` con la misma reducción de Dirichlet/MPC (`assembler.reduce`/`expand`) y el mismo backend (`corrector.solve`) que el predictor del padre, y devuelve `None` si aún no hay `K_L` o es singular. Todo por la API pública del programa principal; el programa principal no gana piezas nuevas.
- **Predictor**: en modo cilíndrico, el del padre, orientado por `sign(ΔU_prev · du_t)`; en modo de disipación, `ΔU = Δλ·du_L` y `Δλ = τ/α_L`, con valor y signo fijados por la restricción (§5, desde el 2026-09-29 `τ/α`, antes `sign·|τ/α|`; con la tangente de carga desde el 2026-09-30). Sin seguimiento de inercia (§5); el método `_negative_pivots()` (signo del `slogdet` de `K_t`, que ningún `solve` llamaba) se retiró el 2026-09-29.
- **Llegada a `max_lambda`** (desde 2026-09-23, la del padre): todos los pasos llevan su restricción (cilíndrica o de disipación); el que converge por encima de `max_lambda` no se consolida y se repite con λ = `max_lambda` fijo desde el estado anterior (`ArcLengthSolver._land`). Si la llegada falla se biseca dl o τ según el modo del paso. Sustituye a la salvaguarda anterior contra el `final_step` prematuro, que bisecaba cuando el predictor rebasaba `max_lambda` más de 3 veces lo que faltaba y, por debajo de ese factor, imponía `max_lambda` en control de carga.
- **Detección de α≈0 en el predictor**: `|α| < ALPHA_LINEAR_RTOL·½·|λ_n·F·du|` (`ALPHA_LINEAR_RTOL = 1e-6`, umbral relativo) o α nulo a precisión de máquina (`_alpha_is_null`: `|α| ≤ ZERO_TOL·½·(|λ_n·F·du| + |F·U_n|)`, relativa desde el 2026-09-30; antes `|α| < ZERO_TOL` en unidades de energía, §5) revierte a cilíndrico. En problemas lineales monótonos α≡0 exactamente, y con la tangente de carga también cuando ésta vuelve a ser la secante (descarga, daño saturado). Hasta el 2026-09-29 ocurría además con el daño creciendo, porque el predictor usaba la tangente del estado consolidado, que es la secante (§5), y la restricción de disipación no llegaba a usarse. En el corrector, α nulo con la misma guarda aborta el paso (`CorrectionAborted`) y se biseca τ.
- **Bisección**: según el modo del paso fallido, `dl *= 0.5` o `τ *= 0.5`, con aborto (`RuntimeError`) por debajo de `ARCLENGTH_MIN_DL_FACTOR` veces `dl_ref` o `initial_tau`. Si el predictor encuentra la tangente singular, divide a la vez `dl` y `τ` y reintenta.
- **Entrypoint público**: registrado vía `@SolverRegistry.register` al cargar el módulo; en YAML, `user_modules: [dissipation_control]` y `solver: {type: DissipationArcLengthSolver, ...}`; en Python, `solidum.load_user_module("dissipation_control")` o `from solidum.user.dissipation_control import DissipationArcLengthSolver`. Sin el módulo cargado, el YAML falla con un error que sugiere declararlo y lista los módulos disponibles.

### Tests (suite [tests/user/dissipation_control/](../../../../tests/user/dissipation_control/))

[`test_dctl_dissipation_arclength.py`](../../../../tests/user/dissipation_control/test_dctl_dissipation_arclength.py) (hasta el 2026-09-29, `tests/test_dissipation_arclength.py`) cubre:

- ✅ Recuperación cilíndrica en régimen elástico puro (regresión a paridad de bits con `ArcLengthSolver`; todos los problemas de Newton cilíndricos).
- ✅ Daño 1D con softening pronunciado (rama post-pico atravesada, con pasos post-pico en modo de disipación).
- ✅ Switching automático cilíndrico→disipación (acceptance `switching_cilindrico_a_disipacion_en_pico`): problemas de Newton cilíndricos hasta el paso en que empieza el daño, de disipación desde el siguiente, sin oscilación.
- ✅ Daño 2D continuo (Quad4 + IsotropicDamage2D plane strain), con pasos post-pico en modo de disipación.
- ✅ Balance energético sobre una barra de daño continuo 1D, sin grieta ni cohesivo: ΣΔE_d > 0, cada paso de disipación con ΔE_d > 0 y, desde el 2026-09-30, la ΣΔE_d de los pasos en modo de disipación igual a la disipación analítica del material entre sus estados inicial y final, $E_d(\kappa) = W(\kappa) - \Psi(\kappa)$ (tolerancia derivada 1e-4; medido 1.3e-10).
- ✅ Daño 2D localizado con retroceso (revisión adversaria del 2026-09-30): franja de 4×2 Quad4 con `IsotropicDamage2D` en tensión plana (ν = 0.2, medida energética por omisión) y la columna central un 5 % más débil. El daño se localiza en la columna, el extremo retrocede, los problemas de Newton son cilíndricos hasta el paso de inicio del daño y de disipación desde el siguiente (67 pasos), y la ΣΔE_d de esos pasos coincide con la disipación del material en los puntos de Gauss, $\sum V_{gp}\,[D(\kappa_{\text{fin}}) - D(\kappa_{\text{ini}})]$ con $D(\kappa) = \int \tfrac12 E\kappa^2 d'(\kappa)\,d\kappa$, dentro de una cota sin márgenes ajustados: residuo de equilibrio en los extremos más el error del trapecio estimado por duplicación de paso (medido 5.9e-4 J sobre 3.23 J, frente a una estimación del trapecio de 5.3e-4 J). Con el solver anterior al 2026-09-30 falla. Fuera de la suite, la misma franja con 10 columnas sigue la curva de `IndirectDisplacementSolver` (distancia máxima 7.4e-6, mismo $u_{\min}/u_{\text{pico}} = 0.61494$).

Qué restricción resolvió cada paso lo comprueba un espía del corrector compartido ([`_dctl_espia.py`](../../../../tests/user/dissipation_control/_dctl_espia.py)), que anota el modo de cada problema de Newton. *Nota 2026-09-29*: hasta entonces los tests 2-5 sólo observaban el modo `"dissipation"` en el callback, sin que ninguno ejercitara la restricción de disipación: el predictor siguiente revertía el modo a cilíndrico por $\alpha \approx 0$ (§5), y todos los problemas de Newton de esos trazados eran cilíndricos. Desde el 2026-09-30 esos trazados resuelven la rama de ablandamiento en modo de disipación (495 de 500 pasos en la barra 1D; 291 de 300 en el Quad4) y los tests pasan sin cambiar sus valores esperados.

[`test_dctl_retroceso.py`](../../../../tests/user/dissipation_control/test_dctl_retroceso.py) cubre `retroceso_barra_dano_localizado`: la barra de 10 elementos con daño localizado del test de `IndirectDisplacementSolver`, con sus tolerancias: puntos post-pico sobre la curva exacta, retroceso hasta el mínimo exacto y sólo el elemento débil dañado. Del 2026-09-29 al 2026-09-30 los dos primeros iban con `xfail(strict=True)` (el arco cilíndrico se atascaba en el pico sin conmutar y abortaba, con $u_{\min}/u_{\text{pico}} = 1.00001$); desde el 2026-09-30 pasan sin marca (172 puntos, distancia máxima $3.2\cdot10^{-4}$, $u_{\min}/u_{\text{pico}} = 0.616376$). Dos tests más: la rama de ablandamiento del elemento débil sale entera de problemas de Newton en modo de disipación (cilíndricos antes y en el paso de inicio del daño, y de nuevo tras la saturación) y el límite conocido con `initial_dlambda = 0.1` (§10): el paso cilíndrico que cruza el pico salta a la rama con los diez elementos ablandándose, que el solver sigue después en modo de disipación sobre su curva analítica. De la revisión adversaria del 2026-09-30: todos los puntos con daño, también los de la rama de rigidez residual tras la saturación (hasta $\varepsilon_w \approx 230\,\kappa_{0w}$, que el filtro $\varepsilon_w \le 60\,\kappa_{0w}$ del test de la curva deja fuera), sobre la solución exacta evaluada en su propia $\varepsilon_w$ (medido 7.7e-6); y el trazado no depende de las unidades: con el modelo en mm (energías ×1e3) y en unidades $10^6$ veces mayores (energías ×1e-12) la secuencia de intentos es idéntica y los puntos coinciden con los de SI salvo redondeo.

Los dos tests del `_negative_pivots` se retiraron con el método. [`test_dctl_contratos.py`](../../../../tests/user/dissipation_control/test_dctl_contratos.py) cubre el contrato de solver del barrido del programa principal, el primer paso heredado, `initial_tau`, el gancho de inicio de paso (ADR 0020, P5) y la carga desde YAML con y sin `user_modules`; [`test_dctl_specs.py`](../../../../tests/user/dissipation_control/test_dctl_specs.py), esta spec frente a la clase registrada.

### Validación parcial — qué falta para status `validated`

Cubiertos desde el 2026-09-30: `retroceso_barra_dano_localizado` y `switching_cilindrico_a_disipacion_en_pico`. No están cubiertos: `snap_through_lee_frame` (con no linealidad geométrica el criterio de arranque por $\alpha_L$ queda fuera de su alcance, §5), los tres casos con componentes del módulo de usuario `discontinuities` (`cohesive_embedded_softening_postpeak`, `balance_energetico_GF`, `convergencia_temporal_tau`), `sign_of_pivot_tracking_bifurcacion_simple` (capacidad no implementada, §5) y `fallback_a_cilindrico_en_descarga_global`.

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
- **2026-09-29** · **Cambio de formulación con visto bueno del usuario** (tercera revisión del día). (1) Predictor del modo de disipación $\Delta\lambda = \tau/\alpha$: la restricción lineal fija valor y signo; sustituye a $\text{sign}\cdot|\tau/\alpha|$, que con $\text{sign} \neq \operatorname{sgn}(\alpha)$ arrancaba en $g = -\tau$. (2) Umbral de conmutación relativo, $\texttt{dissipation\_threshold}\cdot\lVert\mathbf F\rVert\,\lVert\mathbf U\rVert$, sin el piso $\max(\cdot, 1)$ en unidades de energía. (3) Alcance restringido a daño con descarga secante (§0, §1 con la derivación de $\Delta E_d$, `requirements`, `out_of_scope`); la variante para plasticidad de Verhoosel et al. (2009) queda como ampliación futura. Los tres puntos abiertos de la revisión anterior quedan así resueltos. Los tests existentes del módulo pasan sin cambios (19 y 6 subtests).
  La hipótesis de que (1) resolvería la deuda #27 **no se cumple**. Nuevo test `test_dctl_retroceso.py` (acceptance `retroceso_barra_dano_localizado`, barra de 10 elementos con daño localizado, tolerancias del test de `IndirectDisplacementSolver`) en `xfail(strict=True)`. Diagnóstico con registro por paso (§5): el predictor $\tau/\alpha$ no llega a ejecutarse, porque con la tangente secante del estado consolidado $\alpha = -\tfrac12\mathbf F^T\mathbf K_s^{-1}\mathbf r_n \approx 0$ y el paso vuelve a cilíndrico; en los tests existentes, con daño 1D y 2D, ningún problema de Newton se construye en modo de disipación. Además, en el pico anguloso de la barra el arco cilíndrico se atasca sin dar un paso con $\Delta E_d$ sobre el umbral: el solver no conmuta y aborta ($u_{\min}/u_{\text{pico}} = 1.00001$). El cambio (2) es neutro en este caso ($\lVert\mathbf F\rVert\,\lVert\mathbf U\rVert \approx 5.3 > 1$: el piso no actuaba). Un prototipo de diagnóstico, fuera del repositorio, con la tangente del último iterado convergido en el predictor y la conmutación por $\alpha_L$, sigue el retroceso con las tolerancias del test; con esa tangente, la regla de signo antigua pierde el trazado ($u_{\min}/u_{\text{pico}} = -1.0$) y la nueva no: (1) es necesario pero no suficiente. Quedan dos decisiones de formulación para el usuario (§5, puntos abiertos 1 y 2).
- **2026-09-30** · **Cambio de formulación por decisión del usuario**: resueltos los puntos abiertos 1 y 2 de §5. (a) **Tangente de carga en el predictor** del modo de disipación: $\delta\mathbf u_L = \mathbf K_L^{-1}\mathbf F$ con $\mathbf K_L$ la tangente del último iterado convergido del paso anterior (prueba de carga respecto de $\kappa_{n-1}$), en lugar de la del estado consolidado, que con daño y descarga secante es la secante y da $\alpha = -\tfrac12\mathbf F^T\mathbf K_s^{-1}\mathbf r_n \approx 0$; entonces $\alpha_L = \tfrac12\lambda_n\,\mathbf F^T(\mathbf K_L^{-1} - \mathbf K_s^{-1})\mathbf F$ (más el ruido del residuo) mide el crecimiento del daño y en el retroceso es negativo, así que $\Delta\lambda = \tau/\alpha_L < 0$. (b) **Arranque del modo de disipación por $\alpha_L$**: conmuta cuando $|\alpha_L| > 10^{-3}\cdot\tfrac12|\lambda_n\,\mathbf F\cdot\delta\mathbf u_L|$, en lugar de exigir $\Delta E_d > \varepsilon_{\text{diss}}$ tras un paso cilíndrico; al entrar, $\tau \leftarrow$ `initial_tau`. La vuelta a cilíndrico no cambia ($\alpha \approx 0$ en el predictor, $\Delta E_d < \varepsilon_{\text{diss}}$ tras el commit); `dissipation_threshold` sólo interviene en ella. No se añadió la condición $\Delta E_d \ge 0$ al arranque: el prototipo no la usaba, y $\alpha_L \neq 0$ exige que el paso anterior acabara con daño en crecimiento, que disipa. Implementación sólo en el módulo, con API pública del programa principal (`CorrectorResult.state`, `corrector.solve`, `assembler.reduce`/`expand`); el programa principal no cambia. Resultados, iguales a los del prototipo: la barra con daño localizado sigue el retroceso exacto (172 puntos post-pico, distancia máxima $3.2\cdot10^{-4}$, $u_{\min}/u_{\text{pico}} = 0.616376$, sólo el elemento débil dañado) y los dos tests de `test_dctl_retroceso.py` pasan sin `xfail` con las tolerancias del test de `IndirectDisplacementSolver`; el acceptance `retroceso_barra_dano_localizado` queda validado (deuda #27). Tests nuevos: un espía del corrector (`_dctl_espia.py`) comprueba el modo de los problemas de Newton, no sólo el `_mode` del callback; en la barra, la rama de ablandamiento sale entera en modo de disipación y el solver vuelve a cilíndrico al saturar el daño; en los tests existentes, la rama post-pico se resuelve ahora en modo de disipación (495 de 500 pasos en la barra 1D, 291 de 300 en el Quad4) y pasan sin cambiar ningún valor esperado; el switching queda cubierto con el modo de los problemas de Newton (acceptance `switching_cilindrico_a_disipacion_en_pico`); balance energético cuantitativo en la barra 1D: la $\sum\Delta E_d$ de los pasos de disipación coincide con la disipación analítica del material (1.3e-10 relativo, tolerancia derivada 1e-4). Límites documentados en §10: con un primer paso grande ($\Delta\lambda_1 = 0.1$) el paso cilíndrico que cruza el pico salta a la rama con los diez elementos ablandándose —no unicidad del problema incremental, fijada como caso conocido en un test—; el paso que cruza la saturación del daño salta por la rama residual; el criterio (b) sólo vale con pequeñas deformaciones y descarga secante, y supone `rtol_force` $\ll 10^{-3}$. Pendiente: confirmar en el texto de Verhoosel et al. (2009) qué tangente usa su predictor. `status` sigue en `implemented` (quedan acceptance sin cubrir).
- **2026-09-30** · Revisión adversaria de (a) y (b), sin cambio de formulación. Comprobado: $\mathbf K_L$ es la tangente del iterado convergido del paso anterior (en la barra coincide con la tangente analítica de carga del elemento débil a 1e-15 en los 399 predictores) y no toca el estado de los elementos (no se ensambla); tras un paso fallido, también con fallos inyectados en pasos cilíndricos y de disipación, el intento siguiente usa la misma $\mathbf K_L$ y el trazado sigue sobre la curva exacta; el paso que rebasa `max_lambda` no la actualiza, porque no se consolida. En la rama elástica $|\alpha_L|/\text{escala} \le 1.6\cdot10^{-15}$. Corregido: la guarda de $\alpha$ nulo del predictor y del corrector comparaba $|\alpha|$ con `ZERO_TOL` en unidades de energía y rompía el trazado con energías de ~1e-12; pasa a ser relativa (§5, «Guarda de $\alpha$ nulo»). Tests nuevos: invariancia frente a las unidades y puntos de la rama residual en `test_dctl_retroceso.py`; daño 2D localizado con retroceso y balance de disipación en `test_dctl_dissipation_arclength.py` (Tests, arriba). Los valores esperados de los tests existentes no cambian.
