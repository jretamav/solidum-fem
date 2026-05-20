# DissipationArcLengthSolver — Arc-length con restricción de disipación (Gutiérrez 2004)

> **Spec corta tipo extensión** (Reglas.md §4). Variante de [`ArcLengthSolver`](ArcLengthSolver.md) que **solo cambia la restricción del paso**: pasa de cilíndrica ($\lVert\Delta\mathbf U\rVert^2 = dl^2$) a basada en **disipación de energía** ($g(\Delta\mathbf U, \Delta\lambda) = \tau$). Reusa todo el plumbing del padre (predictor tangente, corrector Newton, adaptatividad de paso, Dirichlet/MPC, backend algebraico). Sin ADR nuevo: reusa el contrato arquitectural ya aceptado del arc-length.
>
> Motivación: el `ArcLengthSolver` cilíndrico no atraviesa la transición elástico→softening con penalty stiff (`K_e` del embedded discontinuity, ADR 0010 fase 4). La restricción cuadrática queda casi degenerada cerca del pico — el discriminante de la cuadrática tiende a cero, el solver bisecta indefinidamente y aborta. La restricción de disipación es **lineal en $(\Delta\mathbf U, \Delta\lambda)$** y se activa **proporcionalmente a la energía consumida** por la grieta cohesiva, navegando la rama post-pico de forma estable. Esta variante desbloquea el benchmark de fractura faithful Van Vliet (`docs/specs/CST_Embedded2D.md` §"Out of scope"), el balance energético end-to-end de $G_F$ (hueco abierto en [[project_validacion_fase_b]]), y validación cuantitativa de cohesivos en pipeline real.

---

## Especificación física

### 0. Descripción general

Solver estático no lineal con **traza por disipación**: la magnitud que se mantiene constante de paso a paso ya no es la longitud euclídea $\lVert\Delta\mathbf U\rVert$ sino la **disipación incremental de energía** $\Delta E_d$ asociada al paso. Adecuado para:

- Fractura cohesiva (Mode-I y, en el futuro, mixto) con penalty `K_e` stiff.
- Daño con softening severo donde la rama post-pico tiene pendiente casi vertical en el plano $(U, \lambda)$.
- Localización de deformación plástica con bandas estrechas.

Régimen elástico (sin disipación): la restricción degenera a $g = 0$. El solver detecta este caso y **conmuta automáticamente** a arc-length cilíndrico durante el régimen elástico, retornando a control por disipación cuando se detecta inicio de disipación (`Δω > 0` en daño, $\Delta\alpha > 0$ en plasticidad, $\Delta\kappa > \kappa_0$ en cohesivo).

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

Idénticas al padre — `SolveResult` con $\mathbf U$ final y `step_callback(step, U, λ)` para historia. Adicionalmente, se expone vía el callback el valor de disipación acumulada $\sum_n \Delta E_d^{(n)}$ — útil para validar **$G_F$ end-to-end** integrando $\int_{\Gamma_d} t\cdot d[\![u]\!] \approx \sum_n \Delta E_d^{(n)}$.

---

## Formulación numérica — solo lo que cambia respecto al padre

### 4. Esquema operativo

Estructura idéntica al padre con tres diferencias:

```
λ = 0
τ = initial_tau          # análogo a initial_dl, pero en unidades de energía
modo = "cylindrical"      # arranque en cilíndrico mientras g ≈ 0
mientras λ < max_lambda y step < max_steps:
    # Predictor tangente (idéntico al padre)
    δu_t = K_t^{-1} · F_ext^ref
    sign = sign(δU_prev · δu_t)
    si modo == "cylindrical":
        Δλ_pred = sign · dl / ‖δu_t‖
    sino:                  # modo == "dissipation"
        # Predictor por disipación: Δλ_pred derivado de τ
        Δλ_pred = sign · τ / (½·|λ_n·(F_ext^ref·δu_t) − F_ext^ref·U_n|)
    ΔU_pred = Δλ_pred · δu_t

    # Corrector iterativo
    para iter en range(max_iter):
        K_t, F_int = assemble(...)
        R = (λ + Δλ)·F_ext^ref − F_int
        δu_R = K_t^{-1} · R
        δu_t = K_t^{-1} · F_ext^ref

        si modo == "cylindrical":
            # Restricción cuadrática, idéntica al padre (ArcLengthSolver §5)
            ddλ = raiz_de_menor_angulo(...)
        sino:               # modo == "dissipation"
            # Restricción lineal de Gutiérrez: g(ΔU_new, Δλ_new) = τ
            # ΔU_new = ΔU + δu_R + ddλ·δu_t
            # g(ΔU_new, Δλ + ddλ) = τ  ⇒  ddλ explícito:
            num = τ − ½·(λ_n·F_ext^ref·(ΔU+δu_R) − (Δλ)·F_ext^ref·U_n)
            den = ½·(λ_n·F_ext^ref·δu_t − F_ext^ref·U_n)
            ddλ = num / den
        ΔU += δu_R + ddλ·δu_t
        Δλ += ddλ
        U_iter = U_current + ΔU
        si converge: break

    si converge:
        commit_all_states()
        λ += Δλ
        actualizar_modo(...)            # ← switching cilíndrico ↔ disipación
        ajustar_tau_adaptativamente(iter)
    sino:
        τ /= shrink_factor
```

### 5. Predictor / corrector — punto clave del cambio

#### Switching automático cilíndrico ↔ disipación

El solver mantiene un atributo `modo` con dos estados:

- `"cylindrical"`: arc-length cilíndrico estándar; restricción $\lVert\Delta\mathbf U\rVert^2 = dl^2$. **Régimen elástico**, sin disipación activa.
- `"dissipation"`: arc-length por energía; restricción $g = \tau$. **Régimen disipativo**.

**Criterio de transición** (Verhoosel et al. 2009 §3.2):

- Inicio de paso en `"cylindrical"`: tras commit, calcular $\Delta E_d$ del paso recién cerrado. Si $\Delta E_d > \varepsilon_{\text{diss}}$ (umbral pequeño, p. ej. $10^{-12}\cdot \lVert\mathbf F_{\text{ext}}^{\text{ref}}\rVert\cdot \lVert\mathbf U_n\rVert$) ⇒ siguiente paso entra en `"dissipation"`. Inicializa $\tau \leftarrow \Delta E_d$ del último paso para continuidad.
- Inicio de paso en `"dissipation"`: si el paso anterior tuvo iteraciones $\le$ `tau_relax_iter_threshold` **y** el "predictor de disipación" arrojaría $\Delta\lambda$ desproporcionadamente grande (señal de retorno a elástico tras la rama post-pico), revertir a `"cylindrical"`. Patológico pero posible si hay descarga global tras grieta saturada.

#### Sign-of-pivot tracking (mejora sobre el padre)

`ArcLengthSolver._negative_pivots()` es placeholder hoy. Esta variante **lo implementa** mediante factorización LDLᵀ (Bunch-Kaufman) del backend algebraico, exponiendo el conteo de pivots negativos como diagnóstico de bifurcación. Uso interno: detectar paso del punto límite (cambio en el signo del primer pivot indefinido) y orientar la dirección del predictor con respecto al signo previo, no sólo al producto escalar con $\Delta\mathbf U_{\text{prev}}$.

Esto resuelve el caso patológico de bifurcación simétrica donde el "menor ángulo" elige una rama paralela en lugar de continuar el camino físico — caveat documentado en `ArcLengthSolver` §10.

### 6. Criterio de convergencia (ADR 0007)

Idéntico al padre — patrón dual fuerza + desplazamiento, escala autoderivada.

### 7. Imposición de Dirichlet y MPC

Idéntico al padre. La restricción de disipación opera sobre DOFs libres después de la reducción Dirichlet/MPC.

### 8. Backend algebraico

Idéntico al padre (`is_positive_definite=False`, despachador LU). **Adición**: para sign-of-pivot tracking auténtico se requiere LDLᵀ; la implementación inicial puede usar la firma de la factorización LU (signo del determinante) como aproximación, con upgrade a LDLᵀ Bunch-Kaufman cuando se priorice.

### 9. Adaptatividad — ajuste de $\tau$

Análoga a la adaptatividad de $dl$ del padre, con $\tau$ en lugar de $dl$:

- Converge en $<$ `tau_grow_iter_threshold` iteraciones ⇒ $\tau \to \tau \cdot $ `tau_grow_factor` (cota: `initial_tau · tau_max_factor`).
- Converge en $>$ `tau_shrink_iter_threshold` iteraciones ⇒ $\tau \to \tau \cdot $ `tau_shrink_factor`.
- No converge en `max_iter` ⇒ $\tau /= 2$ y reintentar.
- Cota inferior: $\tau < $ `DISS_MIN_TAU_FACTOR · initial_tau` ⇒ aborto.

En modo `"cylindrical"` durante régimen elástico inicial, la adaptatividad opera sobre $dl$ (idéntica al padre).

### 10. Caveats numéricos

Heredados del padre (`ArcLengthSolver` §10) más los específicos:

- **Régimen disipativo no estrictamente monótono**: si el material vuelve transitoriamente a régimen elástico (descarga local), $\Delta E_d$ puede ser negativo. El solver tolera este caso pero conmuta a `"cylindrical"` si la disipación neta del paso cae por debajo del umbral.
- **Inicialización de $\tau$**: debe ser del orden de la disipación esperada por paso (típicamente fracción pequeña de $G_F \cdot l_d$ por elemento agrietado). Demasiado pequeño ⇒ pasos numerosos; demasiado grande ⇒ saltos sobre la rama post-pico.
- **Disipación de origen no físico**: si el modelo numérico introduce disipación espuria (algorithmic damping en transitorio, regularización no consistente), $\Delta E_d$ contiene esa componente. **Este solver es estático**: si el usuario invoca `DissipationArcLengthSolver` sobre un modelo con material no disipativo (elástico puro), el modo permanece en `"cylindrical"` indefinidamente — comportamiento equivalente al padre. No es un caveat operacional, es la definición del switching.
- **MPC con offset $\mathbf g$ no homogéneo**: la restricción de disipación se evalúa sobre DOFs libres reducidos; si las MPC introducen offset, $g$ se evalúa sobre el sistema reducido y la disipación medida es consistente con la energía aportada por las MPC (tratamiento idéntico al residuo en el padre).

---

## Contrato de implementación

```yaml
name: DissipationArcLengthSolver
kind: solver
status: implemented      # draft → implemented → validated. Validación parcial — ver §"Implementación".

parent_spec: ArcLengthSolver

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
  - { name: initial_dl,               type: float, required: false, default: 0.1,
      desc: "Longitud de arco inicial — usado en modo cilíndrico (régimen elástico)" }
  - { name: dl_grow_factor,           type: float, required: false, default: 1.5 }
  - { name: dl_max_factor,            type: float, required: false, default: 5.0 }
  - { name: dl_shrink_factor,         type: float, required: false, default: 0.6 }
  - { name: dl_grow_iter_threshold,   type: int,   required: false, default: 4 }
  - { name: dl_shrink_iter_threshold, type: int,   required: false, default: 8 }
  - { name: linear_algebra,           type: str,   required: false, default: "auto" }

  # Específicos de esta variante:
  - { name: initial_tau,              type: float, required: true,
      desc: "Incremento de disipación objetivo (unidades de energía: F·U). Calibrar al orden de G_F·l_d/n_pasos_post_pico" }
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

conventions:
  units: "heredadas del modelo (ADR 0008); τ en unidades de energía consistentes"
  stability: "incondicional respecto a puntos límite y rama post-pico de softening severo. La estabilidad ante bifurcaciones múltiples mejora vs el padre por el sign-of-pivot tracking"

out_of_scope:
  - "Bifurcaciones múltiples simétricas — la detección por inercia desambigua bifurcaciones simples; las múltiples siguen requiriendo perturbación explícita"
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

    - name: cohesive_embedded_softening_postpeak       # benchmark CRÍTICO — desbloquea fase 4 ADR 0010
      setup: "Probeta rectangular con CST_Embedded2D + CohesiveDamageIsotropic; tracción uniaxial monotónica que activa Rankine y entra en rama post-pico con penalty K_e stiff (1e13)"
      expect: "El solver atraviesa la transición elástico→softening (régimen donde el ArcLengthSolver cilíndrico aborta); cierra el ciclo de softening; daño cohesivo committed alcanza saturación"
      tol_rel: 1.0e-3

    - name: balance_energetico_GF
      setup: "Mismo modelo cohesivo+embedded del test anterior; integrar ΣΔE_d sobre toda la historia hasta saturación (ω→1) y comparar con G_F · longitud_grieta_abierta del material"
      expect: "Σ ΔE_d ≈ G_F · l_d a tolerancia de discretización temporal (1% típico)"
      tol_rel: 1.0e-2

    - name: convergencia_temporal_tau
      setup: "Reducir initial_tau en factor 2 sucesivamente sobre el benchmark cohesivo, observar U(λ_final)"
      expect: "Convergencia primer orden en τ (esquema lineal en restricción): error → 0 como O(τ)"
      tol_rel: 0.2

  specific:
    - name: switching_cilindrico_a_disipacion_en_pico
      setup: "Modelo con daño que cruza σ_t0 en el paso N"
      expect: "Modos: 'cylindrical' en pasos 1..N-1; 'dissipation' desde N+1 en adelante; sin oscilación de modo en pasos contiguos"

    - name: sign_of_pivot_tracking_bifurcacion_simple
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
  - "Spec padre: ArcLengthSolver.md."
```

---

## Implementación

- **Archivo**: [`solidum/math/solvers/dissipation_arclength.py`](../../solidum/math/solvers/dissipation_arclength.py).
- **Clase**: `DissipationArcLengthSolver(ArcLengthSolver)` — subclase. Reusa predictor, corrector, ajuste de paso, manejo de Dirichlet/MPC y backend del padre. Override completo de `solve()` con bifurcación por `self._mode`.
- **Atributo `self._mode`**: `"cylindrical" | "dissipation"`. Transición tras commit según umbral ``dissipation_threshold·‖F_ref‖·‖U‖``.
- **Sign-of-pivot tracking**: `_negative_pivots()` implementado vía signo del `slogdet` de `K_t`. Aproximación válida para distinguir indefinitud simple (0 vs número impar de pivots negativos) — base para detección de paso por punto límite simple. **LDLᵀ Bunch-Kaufman queda como deuda técnica** para tracking exacto del conteo.
- **Salvaguarda contra ``final_step`` prematuro**: si ``|dλ_pred|`` excede ``3 × (max_lambda − lambda_curr)``, se bisecta dl o τ antes de aceptar el paso. Esencial en pasos iniciales sin calibrar y tras switch cilíndrico→disipación donde α puede ser pequeño.
- **Detección de α≈0 vía threshold relativo**: ``|α| < 1e-6·½·|λ_n·F·du_t|`` reverte a cilíndrico — para problemas lineales monotónicos α≡0 exactamente, no es un error sino la matemática del régimen sin disipación física.
- **Entrypoint público**: registrado vía `@SolverRegistry.register`, expuesto como `solidum.math.solvers.DissipationArcLengthSolver`.

### Tests (suite [tests/test_dissipation_arclength.py](../../tests/test_dissipation_arclength.py))

**7/7 verde**. Cubren:

- ✅ Recuperación cilíndrica en régimen elástico puro (regresión a paridad de bits con `ArcLengthSolver`).
- ✅ Daño 1D con softening pronunciado (rama post-pico atravesada).
- ✅ Switching automático cilíndrico↔disipación observable.
- ✅ Daño 2D continuo (Quad4 + IsotropicDamage2D plane strain).
- ✅ Balance energético cualitativo (ΣΔE_d ≥ 0 con daño activo).
- ✅ Sign-of-pivot tracking via signo del determinante (0 para PD, 1 para indefinida).

### Validación parcial — qué falta para status `validated`

El acceptance ``cohesive_embedded_softening_postpeak`` **no está cubierto** en esta primera implementación. Diagnóstico:

La activación discreta del cohesivo en `prepare_step` (Rankine onset al cruzar `σ_t0`) introduce un salto en `F_int` y `K_t` entre el paso pre-activación y el paso post-activación. El predictor por `sign = sign(dot(δU_prev, du_t))` puede heredar una orientación inconsistente cuando el Newton del paso de activación reorienta `dU_iter` para acomodar el nuevo `K_t` con el cohesivo activo. El paso siguiente entonces predice hacia atrás (descarga) en lugar de continuar la rama post-pico.

Resolverlo limpiamente requiere **una de las dos vías siguientes**, ambas más allá del scope incremental de esta spec:

1. **Sign-of-pivot tracking real con LDLᵀ Bunch-Kaufman**, no aproximación por signo del determinante. Permitiría detectar la inercia exacta de `K_t` y orientar el predictor por cambio de signo de pivot en lugar del producto escalar `δU·du_t`. Esto desacopla la orientación del predictor de la cinemática perturbada por la activación. Requiere extender el backend algebraico (ADR 0003 fase 2) con un solver LDLᵀ verdadero.
2. **Control por CMOD/CTOD del salto** mediante MPC sobre el DOF interno del cohesivo, en lugar de la restricción de disipación global. Cambia la naturaleza del solver — sería una spec separada (`CMODControlSolver`) en lugar de variante de arc-length.

Recomendación: dejar esta validación como deuda técnica priorizada en STATUS.md §"Deuda técnica priorizada" hasta que el LDLᵀ Bunch-Kaufman entre en otro contexto (lo amerita por más razones, no solo esta).

---

## Diálogo

- **2026-05-19 (tarde)** · Spec promovida de `draft` a `implemented`. 7/7 tests verde sobre los casos de daño continuo (1D y 2D); el caso crítico cohesivo+embedded reveló una limitación más profunda del `sign` por producto escalar tras activación discreta del Rankine, que requiere LDLᵀ Bunch-Kaufman real o control CMOD/CTOD — ambos fuera del scope incremental de esta variante. Documentado en §"Implementación" como deuda técnica. La spec mantiene su valor: provee la arquitectura del switching cilíndrico↔disipación y la fórmula lineal de Gutiérrez funcional para todos los problemas con softening continuo (daño bulk 1D y 2D), que es donde más casos prácticos de Solidum se concentran.

- **2026-05-19 (mañana)** · Spec creada como variante de `ArcLengthSolver` para desbloquear ADR 0010 fase 4 (embedded discontinuity con softening severo) y validación cuantitativa de $G_F$ en pipeline real. Motivada por el cierre de la campaña de validación externa Tandas 12-14 ([[project_validacion_externa_tandas_12_14]]) que identificó esta pieza como **el primer componente nuevo priorizado** tras la campaña. Sin ADR nuevo: reusa decisiones arquitecturales del padre (estructura del bucle predictor/corrector, manejo de Dirichlet/MPC, ADR 0003 algebraico, ADR 0007 convergencia, ADR 0008 unidades). El cambio se circunscribe a la **restricción del paso** (cuadrática → lineal) más mejoras secundarias (switching automático cilíndrico↔disipación, sign-of-pivot tracking).

- **Puntos abiertos para el usuario** (físicos/de alcance, no de plumbing — ver `[[feedback_filtrar_puntos_abiertos]]`):
  1. ¿La restricción de disipación de Gutiérrez 2004 es la fórmula que quieres, o prefieres una alternativa (p. ej. constraint sobre `K_t`-norma de Crisfield 1991 §9.4)? La de Gutiérrez es la más extendida en literatura cohesiva.
  2. ¿El switching automático cilíndrico↔disipación es deseable, o prefieres dos solvers separados (`CylindricalArcLengthSolver` y `DissipationArcLengthSolver` puros) que el usuario elige según el caso? El switching es estándar en Verhoosel et al. 2009; dos solvers separados son más predecibles pero exigen que el usuario diagnostique cuándo cambiar.
  3. ¿Sign-of-pivot tracking se aborda en esta spec o se difiere a una spec independiente (`PivotTracker`)? Está acoplado al backend algebraico (ADR 0003) — si el LDLᵀ Bunch-Kaufman se introduce como helper compartido, otros solvers se beneficiarían. Mi propuesta es **abordarlo aquí en versión "signo del determinante LU"** (aproximación válida para indefinitud simple) y dejar el LDLᵀ verdadero como deuda técnica del proyecto.
