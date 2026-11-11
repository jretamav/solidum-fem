# COHESIVE DAMAGE MIXED MODE — ley cohesiva de modo mixto con apertura efectiva y daño escalar

> Orden de trabajo pre-redactada por la IA el 2026-11-10 con las decisiones físicas del usuario del mismo día: (1) ley de modo mixto con apertura efectiva de Camacho y Ortiz (1996), con energía de fractura independiente del modo (`G_II = G_I`, resistencia a cortante `β·σ_t0`); la energía de modo II mayor, con fricción y dilatancia, queda para la ley elastoplástica normal–tangencial que viene después; (2) tres formas de ablandamiento: lineal, exponencial y bilineal del fib Model Code 2010. Los datos que respaldan las decisiones están en los informes de búsqueda del 2026-11-10 (fuera del repositorio). El usuario revisa la física y aprueba antes de que la IA escriba código.

---

## Especificación física

### 0. Descripción general

Ley tracción–salto de **daño escalar** para una superficie de discontinuidad (grieta) en **modo mixto I–II**, para materiales cuasi-frágiles (concreto, roca, mortero). Generaliza `CohesiveDamageIsotropic` (modo I, la grieta desliza libre) acoplando el deslizamiento al daño mediante una **apertura efectiva** con un peso `β` sobre el deslizamiento (Camacho y Ortiz 1996; Ortiz y Pandolfi 1999):

- en tracción pura se comporta igual que la ley de modo I: resistencia `σ_t0`, energía `G_F`;
- en cortante puro tiene resistencia `τ_c = β·σ_t0` y la misma energía `G_F`;
- en cualquier trayectoria proporcional disipa exactamente `G_F` (energía independiente del modo);
- con `β = 0` es exactamente `CohesiveDamageIsotropic` (acoplamiento nulo, deslizamiento libre).

Material de interfaz de superficie (ADR 0033, `GEOMETRY = "surface"`), con salto de 2 componentes `(n, s)` en 2D o 3 componentes `(n, s1, s2)` en 3D. Lo puede usar cualquier elemento con una superficie de discontinuidad; hoy, el elemento embebido del módulo de usuario `discontinuities`.

**Por qué `G_II = G_I`** (decisión del usuario, 2026-11-10). Las mediciones de ensayos dominados por cortante dan `G_II` de 8 a 25 veces `G_I` (Bažant y Pfeiffer 1986: 24; Reinhardt y Xu 2000: `G_IIF = 2058 N/m`; Swartz et al. 1988: 8–10), pero ninguno es modo II puro y la cifra incluye la trabazón de agregados bajo compresión. En las simulaciones de ensayos de modo mixto de concreto (Arrea–Ingraffea, Gálvez et al. 1998, probetas de doble entalla), las curvas son prácticamente insensibles a `G_II` entre 1 y 1000 veces `G_I` porque la grieta propaga en modo I (Gálvez, Cendón y Planas 2002); García-Álvarez, Gettu y Carol (2012) usan `G_II/G_I = 1`. La energía de modo II grande es física de fricción y dilatancia, que corresponde a la ley elastoplástica posterior.

### 1. Descomposición del salto

Sin descomposición elástica-plástica (daño puro). En el sistema local de `Γ_d`:

$$\llbracket\mathbf u\rrbracket = \delta_n\,\mathbf n + \boldsymbol\delta_s, \qquad \boldsymbol\delta_s = \delta_s\,\mathbf s \ \text{(2D)} \quad\text{o}\quad \delta_{s1}\,\mathbf s_1 + \delta_{s2}\,\mathbf s_2 \ \text{(3D)}$$

con `δ_n > 0` apertura y `‖δ_s‖` deslizamiento. En 3D la ley es isótropa en el plano tangente: sólo depende de `‖δ_s‖`.

### 2. Ley elástica del salto (sin daño)

$$\mathbf t = K_e\,\mathbf W\,\llbracket\mathbf u\rrbracket, \qquad \mathbf W = \operatorname{diag}(1,\ \beta^2[,\ \beta^2])$$

Rigidez `K_e` en dirección normal y `β²·K_e` en dirección tangencial. `K_e` es una rigidez de penalización (misma guía que la ley de modo I: `K_e ≈ 10·E/ℓ_c`, declarada por el usuario, sin valor por omisión).

### 3. Apertura efectiva y criterio de activación

$$\delta_\text{eff} = \sqrt{\langle\delta_n\rangle^2 + \beta^2\,\|\boldsymbol\delta_s\|^2}, \qquad f = \delta_\text{eff} - \kappa \le 0$$

con `⟨·⟩` la parte positiva (Macaulay). Con la grieta comprimida (`δ_n < 0`) sólo el deslizamiento cuenta: `δ_eff = β‖δ_s‖`.

**Criterio de iniciación equivalente en tracciones.** En la rama elástica (`δ_n ≥ 0`), `t_n = K_e δ_n` y `t_s = β² K_e δ_s`, de modo que `δ_eff = κ_0 = σ_t0/K_e` equivale a

$$\left(\frac{t_n}{\sigma_{t0}}\right)^2 + \left(\frac{\|\mathbf t_s\|}{\beta\,\sigma_{t0}}\right)^2 = 1$$

el criterio cuadrático de iniciación de Camanho y Dávila (2002, ec. 27) con resistencia a cortante `S = β·σ_t0`. La tracción efectiva conjugada de `δ_eff` es `t_eff = √(t_n² + ‖t_s‖²/β²)` para `δ_n ≥ 0`.

### 4. Variable histórica y evolución (Kuhn-Tucker)

$$\kappa_{n+1} = \max\bigl(\kappa_n,\ \delta_{\text{eff},\,n+1}\bigr), \qquad \kappa \ge \kappa_0 = \frac{\sigma_{t0}}{K_e}$$

- **Carga activa** (`δ_eff > κ_n`): `κ` crece y el daño aumenta.
- **Descarga y recarga** (`δ_eff ≤ κ_n`): `κ` congelado; secante al origen con rigidez reducida.
- **Compresión** (`δ_n < 0`): la componente normal recupera la rigidez inicial sin daño, `t_n = K_e·δ_n`, que impide la interpenetración (como la ley de modo I, Alfaiate et al. 2002, p. 667). El deslizamiento **sí** puede hacer crecer el daño bajo compresión (`δ_eff = β‖δ_s‖`), y la tracción tangencial se degrada con él: **no hay fricción**. Bajo compresión la resistencia a cortante no aumenta; ver §12.
- Irreversibilidad: `κ` y `ω` monótonos no decrecientes.

Condiciones de Kuhn-Tucker: `κ̇ ≥ 0`, `f ≤ 0`, `κ̇·f = 0`.

### 5. Tracción y envolvente de ablandamiento

Con la secante `S(κ) = T(κ)/κ`:

$$t_n = \begin{cases} S(\kappa)\,\delta_n & \delta_n \ge 0 \\ K_e\,\delta_n & \delta_n < 0 \end{cases}, \qquad \mathbf t_s = \beta^2\,S(\kappa)\,\boldsymbol\delta_s$$

En carga activa (`κ = δ_eff`) la tracción efectiva sigue la envolvente: `t_eff = T(δ_eff)`. La envolvente `T(κ)` es la de la ley de modo I, ahora con tres formas. Para `κ ≤ κ_0`, `T = K_e·κ` (rama elástica). Para `κ > κ_0` se mide la apertura de ablandamiento `w = κ − κ_0` y la energía de ablandamiento

$$H = G_F - \tfrac12\,\sigma_{t0}\,\kappa_0$$

de modo que el área total bajo `T(κ)` (rama elástica más ablandamiento) es exactamente `G_F` para las tres formas:

- **Lineal**: `T = σ_t0·(1 − w/w_c')`, con `w_c' = 2H/σ_t0`. Coincide exactamente con la lineal de `CohesiveDamageIsotropic` (termina en `κ = 2G_F/σ_t0`).
- **Exponencial**: `T = σ_t0·exp(−σ_t0·w/H)`. Coincide exactamente con la exponencial de `CohesiveDamageIsotropic`.
- **Bilineal** (fib Model Code 2010, §5.1.5.2, en su forma general con dos parámetros de forma: la razón de esfuerzo en el quiebre `r` y el factor de apertura final `c`):

  $$T = \begin{cases} \sigma_{t0}\left[1 - (1 - r)\,\dfrac{w}{w_1}\right] & 0 \le w \le w_1 \\[2mm] r\,\sigma_{t0}\,\dfrac{w_c' - w}{w_c' - w_1} & w_1 < w \le w_c' \\[2mm] 0 & w > w_c' \end{cases}, \qquad w_c' = \frac{c\,H}{\sigma_{t0}}, \quad w_1 = \frac{(2 - r\,c)\,H}{\sigma_{t0}}$$

  `w_1` sale de imponer que el área de ablandamiento valga `H`: `½(1+r)σ_t0 w_1 + ½ r σ_t0 (w_c' − w_1) = ½σ_t0 (w_1 + r w_c') = H`. Valores por omisión, los del Model Code 2010: `r = 0.2`, `c = 5`, que dan `w_1 = H/σ_t0` y las ramas `σ_t0(1 − 0.8 w/w_1)` y `σ_t0(0.25 − 0.05 w/w_1)`, las del Model Code con `G_F` sustituida por `H`. Otra calibración usual es la de Petersson (1981): `r = 1/3`, `c = 3.6` (`w_1 = 0.8 H/σ_t0`).

Con `K_e` de penalización, `κ_0 ≪ w` y `H ≈ G_F`: el desplazamiento por `κ_0` es la corrección que hace exacta la energía, como ya hace la forma exponencial de la ley de modo I.

**Energía en trayectorias proporcionales.** Con `δ = λ·(a, b)` y carga monótona (`κ = δ_eff`):

$$\int \mathbf t\cdot d\boldsymbol\delta = \int S\,(\delta_n\,d\delta_n + \beta^2\,\delta_s\,d\delta_s) = \int S\,\delta_\text{eff}\,d\delta_\text{eff} = \int T(\delta_\text{eff})\,d\delta_\text{eff} = G_F$$

para cualquier dirección: la energía no depende del modo. En trayectorias no proporcionales la energía disipada depende de la trayectoria (como en toda ley de daño escalar con apertura efectiva).

### 6. Variables internas

- `kappa : float ≥ κ_0` — máxima apertura efectiva alcanzada.
- `damage : float ∈ [0, 1]` — `ω = 1 − S(κ)/K_e`. **PRIMARY_STATE_VAR** (como la ley de modo I; se exporta al post-proceso).

Mismo esquema de estado que `CohesiveDamageIsotropic`: un elemento puede cambiar de una ley a la otra sin cambiar su manejo del estado.

### 7. Sistema local y notación

El material recibe `[[u]] = (δ_n, δ_s)` o `(δ_n, δ_s1, δ_s2)` y devuelve `t` en los mismos ejes locales; la transformación a ejes globales es del elemento. `δ_n > 0` ⇔ apertura. Signo de `t_s`: el de `δ_s` (la tracción tangencial se opone al deslizamiento relativo de la cara opuesta, como en la ley de adherencia).

---

## Formulación numérica

### 8. Algoritmo (explícito, sin Newton local)

Dado `κ_n` y el salto de prueba:

1. `δ_n⁺ = max(δ_n, 0)`, `δ_eff = √(δ_n⁺² + β²‖δ_s‖²)`.
2. `carga = δ_eff > κ_n`; `κ = max(κ_n, δ_eff)`.
3. `T(κ)`, `T'(κ)` de la envolvente (§5); `S = T/κ`.
4. Tracción por §5. La secante se evalúa como `T/κ`, nunca como `(1 − ω)·K_e` (con `K_e` de penalización, `1 − ω ~ κ_0/κ` perdería cifras).
5. Tangente por §9.

### 9. Tangente algorítmica consistente

Sea `v = (δ_n⁺, β² δ_s)` (el gradiente `∂δ_eff/∂δ = v/δ_eff`) y `S' = dS/dκ = (T'·κ − T)/κ²`.

- **Rama elástica, descarga y recarga** (`δ_eff ≤ κ_n` o `κ ≤ κ_0`): secante
  $$\mathbf D = S\,\mathbf W, \qquad D_{nn} = K_e \text{ si } \delta_n < 0$$
- **Carga activa** (`δ_eff > κ_n`, `κ > κ_0`):
  $$\mathbf D = S\,\mathbf W + \frac{S'}{\delta_\text{eff}}\,\mathbf v\otimes\mathbf v, \qquad \text{con } D_{nn} = K_e,\ D_{ns} = D_{sn} = 0 \text{ si } \delta_n < 0$$
  (con `δ_n < 0`, `v_n = 0` y el término de rango 1 sólo actúa en el bloque tangencial).

Casos particulares de control: en modo I puro con carga, `D = diag(T', β² S)`; con la grieta totalmente abierta (lineal y bilineal, `κ ≥ w_c'`), `T = T' = 0` y `D = 0` en apertura. Sin tope de daño ni rigidez residual, por la misma razón que la ley de modo I (su spec, §9).

### 10. Simetría

`W` es diagonal y el término de carga es un producto `v⊗v`: **`D` es simétrica** en todas las ramas, también con la grieta comprimida. La razón física: la ley deriva de un potencial secante `ψ = ½ S(κ)·(δ⁺ᵀ W δ⁺)` más la penalización de cierre. `IS_SYMMETRIC = True`; el despachador algebraico (ADR 0003) puede seguir usando la factorización simétrica. (La ley de modo I preveía que el modo mixto rompería la simetría; con la apertura efectiva no ocurre. Sí la romperá la ley con dilatancia.)

### 11. Decisión de régimen durante el Newton global

Como la ley de modo I: carga o descarga se deciden con `κ_n` confirmado del paso anterior; las transiciones son discontinuas en la tangente (al cruzar `κ_0`, `κ_n` y `δ_n = 0`).

### 12. Advertencias numéricas y límites del modelo

- **Penalización.** `K_e` en normal y `β²·K_e` en tangencial; con `β` pequeño la rigidez tangencial inicial es pequeña. Misma guía que la ley de modo I; sin valor por omisión.
- **Sin fricción.** Con la grieta cerrada y comprimida, el cortante se degrada con el deslizamiento igual que con la grieta abierta, y la resistencia a cortante no crece con la compresión. En problemas gobernados por cortante con compresión (interfaces, ensayos de empuje directo, cortante en vigas con trabazón) la ley **subestima** el cortante transmitido; para eso es la ley elastoplástica normal–tangencial con fricción y dilatancia (siguiente bloque).
- **Acoplamiento por el daño.** El deslizamiento bajo compresión hace crecer `κ`, y al reabrirse la grieta la resistencia a tracción ya está reducida: es la consecuencia de un daño escalar único (Camacho–Ortiz), no un error.
- **Elemento embebido con cinemática KOS.** La ley transmite cortante; con la cinemática KOS el elemento puede presentar bloqueo de esfuerzos con la grieta muy abierta y esfuerzo paralelo grande (medido el 2026-11-10). Es una propiedad del elemento, no del material; la decisión KOS frente a SKON se toma con el ensayo de referencia de modo mixto.
- **Objetividad de malla** en el ablandamiento: a nivel del elemento, como la ley de modo I.

---

## Contrato de implementación

```yaml
name: CohesiveDamageMixedMode
kind: interface_material        # familia de materiales de interfaz, ADR 0033
status: draft

interface:
  geometry: surface             # superficie de discontinuidad: componentes (n, s) o (n, s1, s2)
  jump_dim: null                # 2 (2D) o 3 (3D)
  primary_state_var: damage
  is_symmetric: true            # potencial secante: D simétrica en todas las ramas (§10)

admission:                      # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - de_borst_crisfield: "cap. 6, modelos de zona cohesiva, y cap. 13, interfaces y discontinuidades (el tratamiento del modo mixto se confirmó por título)"
  general_codes:
    - abaqus: "elementos cohesivos: desplazamiento efectivo (Camanho–Dávila), evolución del daño independiente del modo (opción por omisión), ablandamiento lineal, exponencial o tabular"
    - ansys: "CZM bilineal de modo mixto con salto efectivo ponderado por β (constante C6, Alfano–Crisfield 2001)"
    - code_aster: "CZM_TAC_MIX: ablandamiento acoplado en los tres modos con κ = máx ‖δ‖ y contacto unilateral"

parameters:
  - { name: sigma_t0,  type: float, required: true,  desc: "Resistencia a tracción σ_t0 [Pa]" }
  - { name: G_f,       type: float, required: true,  desc: "Energía de fractura G_F [N/m], la misma en modo I, en modo II y en modo mixto" }
  - { name: K_e,       type: float, required: true,  desc: "Rigidez de penalización normal K_e [Pa/m]; la tangencial es β²·K_e" }
  - { name: beta,      type: float, required: true,  desc: "Peso del deslizamiento en la apertura efectiva (≥ 0); resistencia a cortante τ_c = β·σ_t0. Sin valor por omisión: las calibraciones publicadas de concreto dan τ_c/σ_t0 entre 1.0 y 2.6. β = 0 recupera CohesiveDamageIsotropic" }
  - { name: softening, type: str,   required: true,  desc: "'linear' | 'exponential' | 'bilinear'" }
  - { name: knee_stress_ratio,   type: float, required: false, desc: "Sólo 'bilinear': r = esfuerzo en el quiebre / σ_t0. Por omisión 0.2 (Model Code 2010)" }
  - { name: final_opening_factor, type: float, required: false, desc: "Sólo 'bilinear': c = w_c'·σ_t0/H. Por omisión 5.0 (Model Code 2010)" }

signature:
  compute_traction: "([[u]]: ndarray(2,) o (3,), state_vars=None) -> (t, D: ndarray(n,n), state_vars')"
  jump_kind: "sistema local: (δ_n, δ_s) o (δ_n, δ_s1, δ_s2)"

state_schema:
  kappa: "float ≥ κ_0 — máxima apertura efectiva alcanzada"
  damage: "float ∈ [0, 1] — 1 − S(κ)/K_e"

conventions:
  sign: "δ_n > 0 ⇔ apertura; t_s con el signo de δ_s"
  frame: "ejes locales de Γ_d; transformación a globales por el elemento consumidor"
  units: "[σ_t0] = Pa, [G_F] = N/m, [K_e] = Pa/m, β, r, c adimensionales"

validity:
  - "σ_t0 > 0, G_F > 0, K_e > 0, β ≥ 0"
  - "K_e > σ_t0²/(2·G_F)  (H > 0)"
  - "bilineal: 0 < r < 1, c > 0, r·c < 2 (w_1 > 0) y c·(1 + r) > 2 (w_1 < w_c')"
  - "knee_stress_ratio y final_opening_factor sólo con softening='bilinear' (con otra forma se rechazan)"
  - "salto de 2 o 3 componentes"

out_of_scope:
  - "G_II ≠ G_I (criterios de Benzeggagh–Kenane o de ley de potencia): η sin datos verificados para concreto; posible ley aparte"
  - "Fricción y dilatancia en la grieta: ley elastoplástica normal–tangencial (siguiente bloque)"
  - "Resistencia a cortante dependiente de la compresión"
  - "Degradación cíclica más allá de la descarga secante; fatiga; dependencia de la velocidad"
  - "Función auxiliar para obtener σ_t0 y G_F del Model Code 2010 a partir de f_cm (posible después)"

numerical_caveats:
  - "Penalización K_e (normal) y β²·K_e (tangencial), declaradas por el usuario"
  - "Secante evaluada como T(κ)/κ, sin pasar por 1 − ω (cancelación con K_e de penalización)"
  - "Sin tope de daño ni rigidez residual; tangente consistente en toda la rama"
  - "Sin fricción: subestima el cortante transmitido con la grieta comprimida (§12)"

acceptance:
  verification:
    - name: respuesta_elastica_bajo_umbral
      setup: "saltos con δ_eff < κ_0 en direcciones θ ∈ {0°, 30°, 60°, 90°} del plano (δ_n, β·δ_s), y con δ_n < 0"
      expect: "t = K_e·W·δ exacto (t_n = K_e δ_n con δ_n < 0); ω = 0; D = K_e·W"
      tol_rel: 1.0e-12
    - name: superficie_de_iniciacion
      setup: "carga radial hasta δ_eff = κ_0 en las mismas direcciones"
      expect: "(t_n/σ_t0)² + (t_s/(β σ_t0))² = 1; cortante puro |t_s| = β·σ_t0; tracción pura t_n = σ_t0"
      tol_rel: 1.0e-10
    - name: energia_independiente_del_modo
      setup: "carga monótona proporcional en θ ∈ {0°, 30°, 60°, 90°} hasta la apertura total (exponencial: hasta 30·G_F/σ_t0), para las tres formas, en 2D y en 3D"
      expect: "∫ t·dδ = G_F en todas las direcciones y formas"
      tol_rel: 1.0e-3
    - name: bilineal_model_code_2010
      setup: "softening='bilinear' con r y c por omisión; tracción pura monótona"
      expect: "T(κ_0 + w_1) = 0.2·σ_t0 con w_1 = H/σ_t0; T = 0 desde κ_0 + 5H/σ_t0; ramas σ_t0(1 − 0.8 w/w_1) y σ_t0(0.25 − 0.05 w/w_1); la calibración de Petersson (r = 1/3, c = 3.6) también integra G_F"
      tol_rel: 1.0e-10
    - name: limite_beta_cero
      setup: "β = 0, softening lineal y exponencial; 200 historias aleatorias de salto (apertura, cierre, deslizamiento) aplicadas a esta ley y a CohesiveDamageIsotropic"
      expect: "tracciones, tangentes y estado idénticos a CohesiveDamageIsotropic"
      tol_rel: 1.0e-12
    - name: descarga_y_recarga_secantes
      setup: "carga en modo mixto hasta κ ∈ (κ_0, w_c/2), descarga al origen y recarga en otra dirección por debajo de κ, después por encima"
      expect: "descarga lineal al origen con secante S(κ)·W; sin daño nuevo mientras δ_eff ≤ κ; al superarlo, κ = δ_eff y t_eff = T(δ_eff)"
      tol_rel: 1.0e-10

  specific:
    - name: tangente_simetrica
      setup: "todas las ramas: elástica, carga en apertura, carga con δ_n < 0 y deslizamiento, descarga, grieta totalmente abierta; 2D y 3D"
      expect: "‖D − Dᵀ‖ / ‖D‖ < 1e-14"
      tol_rel: 1.0e-14
    - name: tangente_diferencia_finita
      setup: "K_e = 1e15 (penalización real), las tres formas, κ/w_ref ∈ {0.01, 0.05, 0.3, 0.7, 0.95}, direcciones de carga θ ∈ {0°, 30°, 60°, 90°} y con δ_n < 0; diferencia centrada relativa"
      expect: "D igual a la diferencia finita en carga activa y en descarga"
      tol_rel: 1.0e-6
    - name: compresion_sin_dano_normal
      setup: "δ_n < 0 con δ_s = 0, con la grieta virgen, parcial y totalmente abierta"
      expect: "t_n = K_e δ_n, D_nn = K_e; κ y ω no cambian; t_n continua en δ_n = 0"
      tol_rel: 1.0e-12
    - name: deslizamiento_bajo_compresion_dana_sin_friccion
      setup: "δ_n < 0 fijo y deslizamiento creciente hasta la apertura total"
      expect: "δ_eff = β|δ_s|; |t_s| sigue β·T(β|δ_s|) (sin dependencia de δ_n) hasta 0; t_n = K_e δ_n"
      tol_rel: 1.0e-10
    - name: disipacion_no_negativa
      setup: "1000 historias aleatorias de salto en 2D y 3D, las tres formas"
      expect: "energía disipada acumulada no negativa y no decreciente paso a paso (trabajo menos energía secante almacenada)"
      tol_rel: 1.0e-10
    - name: isotropia_en_el_plano_tangente
      setup: "3D: mismo δ_n y ‖δ_s‖ con δ_s girado en el plano tangente"
      expect: "t_n igual y t_s girado con δ_s; D girada en el bloque tangencial"
      tol_rel: 1.0e-12
    - name: rechazos
      setup: "parámetros fuera de validity; parámetros de la bilineal con otra forma; salto de 1 o 4 componentes"
      expect: "ValueError con mensaje claro"
    - name: elemento_embebido_tangente_condensada
      setup: "CST_Embedded2D (módulo discontinuities) con esta ley, grieta abierta y deslizando en carga activa"
      expect: "tangente condensada del elemento igual a la diferencia finita de sus fuerzas internas"
      tol_rel: 1.0e-5
    - name: elemento_embebido_limite_modo_I
      setup: "tracción uniaxial con CST_Embedded2D y β = 0 hasta la separación completa"
      expect: "la misma respuesta que con CohesiveDamageIsotropic (solución exacta del ensayo existente)"
      tol_rel: 1.0e-10

  arch:
    - name: registro_de_interfaz
      expect: "Registro vía InterfaceMaterialRegistry; sección YAML interface_materials; YAML = Python"
    - name: is_symmetric_attribute
      expect: "CohesiveDamageMixedMode.IS_SYMMETRIC == True"

references:
  - "Camacho, G.T., Ortiz, M. (1996). Computational modelling of impact damage in brittle materials. Int. J. Solids Struct. 33(20-22), 2899-2938. DOI 10.1016/0020-7683(95)00255-3. Apertura efectiva con peso β; la forma δ_eff = √(⟨δ_n⟩² + β²δ_s²) y su uso en compresión (β‖δ_s‖) se verificaron en fuente secundaria (Liu 2020, CMAME 363:112913, §4), no en el original."
  - "Ortiz, M., Pandolfi, A. (1999). Finite-deformation irreversible cohesive elements for three-dimensional crack-propagation analysis. IJNME 44(9), 1267-1282."
  - "Camanho, P.P., Dávila, C.G. (2002). Mixed-mode decohesion finite elements for the simulation of delamination in composite materials. NASA/TM-2002-211737. Ec. (27): criterio cuadrático de iniciación; ec. (28): desplazamiento efectivo de Abaqus; ec. (35): criterio BK (fuera de alcance)."
  - "fib (2013). fib Model Code for Concrete Structures 2010. Ernst & Sohn. §5.1.5.2: relación bilineal esfuerzo–apertura de grieta (verificada sólo en fuente secundaria; confirmar con el texto)."
  - "Petersson, P.-E. (1981). Crack growth and development of fracture zones in plain concrete and similar materials. Report TVBM-1006, Lund Institute of Technology. Calibración bilineal r = 1/3, c = 3.6."
  - "Gálvez, J.C., Cendón, D.A., Planas, J. (2002). Influence of shear parameters on mixed-mode fracture of concrete. Int. J. Fract. 118(2), 163-189. Insensibilidad a G_II en ensayos de modo mixto."
  - "García-Álvarez, V.O., Gettu, R., Carol, I. (2012). Analysis of mixed-mode fracture in concrete using interface elements and a cohesive crack model. Sadhana 37(1), 187-205. G_II/G_I = 1."
  - "Bažant, Z.P., Pfeiffer, P.A. (1986). Shear fracture tests of concrete. Mater. Struct. 19(2), 111-121. Reinhardt, H.W., Xu, S. (2000). Int. J. Fract. 105(2), 107-125. Mediciones de G_II."
  - "Alfaiate, J., Wells, G.N., Sluys, L.J. (2002). Eng. Fract. Mech. 69, 661-686. Penalización y cierre (p. 667), como en la ley de modo I."
  - "Spec CohesiveDamageIsotropic — ley de modo I que esta generaliza (límite β = 0)."
  - "ADR 0033 — Familia de materiales de interfaz."
```

---

## Implementación

*Rellena la IA tras programar.*

- Archivo: —
- Clase: —
- Tests: —
- Notas de traducción: —

---

## Diálogo

- **2026-11-10** · Spec pre-redactada por la IA con las decisiones físicas del usuario del mismo día: ley de Camacho–Ortiz con energía independiente del modo («tu recomendación») y tres formas de ablandamiento, lineal, exponencial y bilineal del Model Code 2010 («si»). Propuesta de la IA incluida para aprobación: `β` obligatorio, sin valor por omisión.
- **2026-11-10** · Spec aprobada por el usuario («si»), con `β` obligatorio. La IA arranca la implementación.
