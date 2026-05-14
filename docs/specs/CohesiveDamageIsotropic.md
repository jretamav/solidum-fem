# COHESIVE DAMAGE ISOTROPIC — material cohesivo traction-jump con daño escalar Modo-I

> Orden de trabajo. Pre-redactada por la IA como borrador para validación del usuario (autor de la formulación, Retama 2010). La física y la formulación numérica son traducción directa del Cap. 3 de la tesis; la mecánica de implementación (Voigt, registry, signatura, tests) sigue las convenciones del proyecto. El usuario revisa, corrige donde proceda, y aprueba antes de que la IA toque código.

---

## Especificación física

### 0. Descripción general

Material cohesivo **traction-jump** para la superficie de discontinuidad `Γ_d` interior a un elemento finito sólido 2D. Modelo de **daño escalar isótropo** con activación tipo Rankine (Modo-I), ablandamiento gobernado por la energía de fractura `G_F`, y bulk vecino que descarga elásticamente sin disipación adicional.

Distinto en dominio de los materiales continuos `Elastic2D`, `VonMises2D`, `IsotropicDamage2D`:

- **No opera sobre puntos de Gauss del bulk.** Opera sobre el salto `[[u]] ∈ ℝ²` definido sobre `Γ_d` y devuelve tracciones `t ∈ ℝ²` definidas sobre la misma superficie.
- **Unidades**: parámetros físicos son `σ_t0` (Pa, resistencia a tracción) y `G_F` (N/m, energía de fractura) — no `E` y `ν`.
- **Es la primera implementación** de la familia `CohesiveMaterial` introducida por ADR 0010.

Pensado para hormigón, roca, cerámicos en régimen cuasi-frágil bajo tracción dominante. Diseñado específicamente para alimentar el elemento `CST_Embedded2D` (fase 2 del ADR 0010), pero la abstracción es reutilizable por futuros elementos de interfaz cohesiva clásica (CZM) si entran al proyecto.

### 1. Descomposición del salto

Sin descomposición elástica-plástica (modelo de daño puro, no plasticidad cohesiva). El salto total `[[u]]` se proyecta sobre el frame local `(n, s)` de `Γ_d`:

$$\llbracket\mathbf u\rrbracket = \llbracket u_n\rrbracket\,\mathbf n + \llbracket u_s\rrbracket\,\mathbf s$$

donde `n` es la normal a `Γ_d` (fijada al activarse el elemento, perpendicular a `σ_I` por criterio Rankine) y `s` es la tangente.

**Modo-I (fase 1)**: sólo `⟨[[u_n]]⟩` (parte positiva del salto normal) gobierna la disipación. La componente tangencial `[[u_s]]` desliza sin restricción y sin disipación. Modo mixto I-II queda diferido a la fase G del ADR 0010.

### 2. Ley elástica del salto (sin daño)

$$\mathbf t = (1 - \omega)\,\mathbf T^{el}\,\llbracket\mathbf u\rrbracket$$

con `T^el` la rigidez elástica del salto:

$$\mathbf T^{el} = K_e\,(\mathbf n \otimes \mathbf n)$$

es decir: rigidez `K_e` en dirección normal, cero en dirección tangencial (consecuencia explícita de Modo-I).

`K_e` es una rigidez de penalty — debe escogerse suficientemente grande para que la fase elástica de la grieta sea geométricamente despreciable (`κ_0 = σ_t0/K_e ≪` escala del problema), pero no tan grande que perjudique el condicionamiento del sistema. Guía típica: `K_e ≈ E_bulk / ℓ_c` con `ℓ_c` una longitud característica del elemento. Ver §12.

### 3. Variable equivalente y función de fluencia

$$\llbracket u\rrbracket_\text{eq} = \langle\llbracket u_n\rrbracket\rangle$$

donde `⟨·⟩` son los brackets de McAuley (parte positiva). La función de fluencia (Eq. 3.13 de Retama 2010):

$$f(\llbracket\mathbf u\rrbracket, \kappa) = \llbracket u\rrbracket_\text{eq} - \kappa \leq 0$$

donde `κ` es el historial del salto equivalente máximo.

**Interpretación**: la grieta sólo "siente" la tracción (y por tanto sólo daña) cuando se abre (`[[u_n]] > 0`). En compresión (`[[u_n]] < 0`) o deslizamiento puro, `f = -κ < 0`, no hay evolución de daño.

### 4. Variable histórica y evolución (Kuhn-Tucker)

$$\kappa_{n+1} = \max\bigl(\kappa_n,\ \langle\llbracket u_n\rrbracket_{n+1}\rangle\bigr), \qquad \kappa_0 = \frac{\sigma_{t0}}{K_e}$$

Codifica:

- **Carga activa** (`⟨[[u_n]]⟩ > κ_n`): `κ` crece; el daño puede aumentar.
- **Descarga / recarga elástica** (`⟨[[u_n]]⟩ ≤ κ_n`): `κ` queda congelado; daño no cambia (secante respecto al origen con rigidez reducida `(1−ω)·K_e`).
- **Compresión** (`[[u_n]] < 0`): `⟨[[u_n]]⟩ = 0 ≤ κ_n`, `κ` congelado. El daño no crece. *Caveat de cierre*: la rigidez efectiva en compresión sigue siendo `(1−ω)·K_e`, no se recupera el contacto rígido — limitación documentada en §12.
- **Irreversibilidad**: `κ` monótono no decreciente; `ω(κ)` monótono no decreciente.

Condiciones de Kuhn-Tucker estándar (Eq. 3.15):

$$\dot\omega \geq 0, \qquad f \leq 0, \qquad \dot\omega\,f = 0$$

### 5. Ley de evolución del daño — softening lineal y exponencial

El daño `ω(κ)` se define implícitamente para que la tracción en carga monótona siga una curva prescrita `T_\text{soft}(κ)`:

$$t_n(\kappa) = (1 - \omega(\kappa))\,K_e\,\kappa \;\stackrel{!}{=}\; T_\text{soft}(\kappa) \;\Rightarrow\; \omega(\kappa) = 1 - \frac{T_\text{soft}(\kappa)}{K_e\,\kappa}$$

con `κ ≥ κ_0`. Para `κ ≤ κ_0`, `ω = 0` por definición.

**Softening lineal**:

$$T_\text{soft}(\kappa) = \sigma_{t0}\left(1 - \frac{\kappa - \kappa_0}{w_c - \kappa_0}\right), \qquad w_c = \frac{2\,G_F}{\sigma_{t0}}$$

válida en `κ ∈ [κ_0, w_c]`. Para `κ ≥ w_c`: `T_\text{soft} = 0`, `ω = 1` (saturada a `DAMAGE_MAX`).

**Softening exponencial**:

$$T_\text{soft}(\kappa) = \sigma_{t0}\,\exp\left(-\frac{\sigma_{t0}\,(\kappa - \kappa_0)}{H}\right), \qquad H = G_F - \frac{\sigma_{t0}\,\kappa_0}{2}$$

asintótica a cero; `ω(κ) → 1` cuando `κ → ∞`, saturada a `DAMAGE_MAX` en la práctica. `H` se deriva imponiendo que la integral total de tracción a lo largo de la curva sea `G_F` (consistencia energética).

La validez de la rama exponencial exige `G_F > σ_{t0}·κ_0/2`, es decir `K_e > σ_{t0}²/(2·G_F)` — condición geométrica trivial cuando `K_e` se elige como penalty.

### 6. Variables internas

- `κ : float` — historial del salto equivalente máximo. **PRIMARY_STATE_VAR alternativa** (más fundamental, controla todo el resto).
- `ω : float ∈ [0, DAMAGE_MAX]` — daño escalar. **PRIMARY_STATE_VAR adoptada** (más interpretable físicamente; se exporta al post-proceso).
- `ω` es función determinista de `κ` (ec. §5); se cachea en el estado para evitar recálculo y para diagnóstico.

### 7. Frame local y notación

Convención del proyecto:

- `n`: normal a `Γ_d`, unitaria. Se fija al activar el elemento (perpendicular a `σ_I`).
- `s`: tangente, ortogonal a `n`, orientada por regla de la mano derecha con `z` saliendo del plano (consistente con Reglas.md §5).
- `[[u]] = [[u_n]]·n + [[u_s]]·s` con `[[u_n]] = [[u]]·n`, `[[u_s]] = [[u]]·s`.
- `t = t_n·n + t_s·s` con `t_n = K_e·(1−ω)·[[u_n]]` en Modo-I, `t_s = 0`.

El material recibe `[[u]]` y devuelve `t` en **ejes locales `(n, s)` de `Γ_d`**. La transformación a ejes globales es responsabilidad del elemento que lo consume, **no** del material. Mismo patrón que en `VonMises2D` (recibe `ε` en Voigt, no se ocupa de transformaciones de elemento).

---

## Formulación numérica

### 8. Algoritmo en un paso (explícito, sin Newton local)

Dado `κ_n` y `[[u]]_{n+1}` (en ejes locales):

1. Calcular `[[u]]_eq = ⟨[[u_n]]⟩ = max(0, [[u]]·n_local)`. *(En el frame local, `n_local = (1, 0)` y `[[u_n]] = [[u]][0]`. Se mantiene la notación abstracta por claridad física.)*
2. Si `[[u]]_eq > κ_n` y `[[u]]_eq > κ_0`: **carga activa**, `κ_{n+1} = [[u]]_eq`. Si no: **descarga / no daño**, `κ_{n+1} = κ_n`.
3. Aplicar `ω(κ_{n+1})` por la fórmula §5; saturar a `DAMAGE_MAX`.
4. Tracción: `t_n = (1 − ω)·K_e·[[u_n]]`, `t_s = 0`.
5. Tangente: ver §9.

El algoritmo es **explícito**: no requiere Newton local. La actualización de `κ` y `ω` es directa. Mismo patrón que `IsotropicDamage2D`.

### 9. Tangente algorítmica consistente

$$\mathbf T^\text{alg} = \frac{\partial\,\mathbf t_{n+1}}{\partial\,\llbracket\mathbf u\rrbracket_{n+1}}$$

Derivando `t = (1−ω)·K_e·[[u_n]]·n` en el frame local con `n = (1,0)`:

- **Sin daño activo** (`κ ≤ κ_0`):

$$\mathbf T^\text{alg} = K_e\,(\mathbf n \otimes \mathbf n) = K_e\begin{pmatrix}1 & 0 \\ 0 & 0\end{pmatrix}$$

- **Descarga** (`[[u]]_eq ≤ κ_n`):

$$\mathbf T^\text{alg} = (1-\omega)\,K_e\,(\mathbf n \otimes \mathbf n)$$

(secante reducida; `∂κ/∂[[u]] = 0` anula el segundo término).

- **Carga activa** (`[[u]]_eq > κ_n`, `κ > κ_0`, `ω < DAMAGE_MAX`):

$$\mathbf T^\text{alg} = (1-\omega)\,K_e\,(\mathbf n \otimes \mathbf n) - \frac{d\omega}{d\kappa}\,K_e\,\llbracket u_n\rrbracket\,(\mathbf n \otimes \mathbf n)$$

$$= \left[(1-\omega) - \frac{d\omega}{d\kappa}\,\llbracket u_n\rrbracket\right]\,K_e\,(\mathbf n \otimes \mathbf n)$$

- **Saturado** (`ω = DAMAGE_MAX`):

$$\mathbf T^\text{alg} = (1-\text{DAMAGE\_MAX})\,K_e\,(\mathbf n \otimes \mathbf n)$$

(`dω/dκ` efectivamente cero al saturar; mismo corte que `IsotropicDamage2D`).

Las derivadas `dω/dκ` cerradas:

- **Lineal**, en `κ ∈ [κ_0, w_c]`:
  $$\frac{d\omega}{d\kappa} = \frac{\sigma_{t0}}{K_e\,\kappa^2}\left(\frac{w_c}{w_c - \kappa_0}\right) - \frac{\sigma_{t0}}{K_e\,\kappa\,(w_c - \kappa_0)} = \frac{\sigma_{t0}\,(w_c - \kappa)}{K_e\,\kappa^2\,(w_c - \kappa_0)} \cdot \frac{1}{1}\;.\;.$$
  *(Forma explícita: simplificar tras derivar `ω(κ) = 1 − T_\text{soft}(κ)/(K_e·κ)`. El usuario revisa el álgebra durante la implementación; la spec exige que la forma cerrada se documente en el docstring de la clase.)*

- **Exponencial**, en `κ > κ_0`:
  $$\frac{d\omega}{d\kappa} = \frac{T_\text{soft}(\kappa)}{K_e\,\kappa}\left(\frac{1}{\kappa} + \frac{\sigma_{t0}}{H}\right) = (1 - \omega(\kappa))\left(\frac{1}{\kappa} + \frac{\sigma_{t0}}{H}\right)$$

  *(Análoga a la forma de `IsotropicDamage2D` con `α ↔ σ_{t0}/H`.)*

### 10. Simetría de la tangente

A diferencia de `IsotropicDamage2D` (donde el producto externo `(C_e·ε) ⊗ (M·ε)` rompe la simetría), aquí ambos vectores que entran al rank-1 colapsan a `n`:

$$\mathbf T^\text{alg} = \alpha(\kappa, \omega, \llbracket u_n\rrbracket)\,K_e\,(\mathbf n \otimes \mathbf n)$$

con `α` escalar. `n ⊗ n` es simétrica, por tanto `T^alg` es **simétrica** en Modo-I.

Consecuencia: `IS_SYMMETRIC = True` para este material. El despachador algebraico (ADR 0003) puede usar Cholesky/LDLᵀ si el resto del dominio también es simétrico. **Modo mixto (fase G) reintroducirá asimetría** — se evaluará entonces.

### 11. Decisión de régimen durante el Newton global

Dentro de una iteración del Newton global, el material recibe `[[u]]_trial` y decide carga/descarga comparando con `κ_n` (estado committed del paso anterior). Mismo patrón que `IsotropicDamage2D` §9. La transición loading↔unloading es discontinua en la tangente; en problemas patológicos cerca del umbral, considerar arc-length (`ArcLengthSolver` ya disponible).

### 12. Caveats numéricos

- **Elección de `K_e` (penalty)**. Si `K_e` muy pequeño: `κ_0 = σ_{t0}/K_e` no es despreciable y la respuesta elástica deforma significativamente la curva (la grieta "se abre" antes del pico). Si `K_e` muy grande: condicionamiento del sistema empeora (`K_e ≫ E_bulk·L_e` introduce diferencia de magnitudes en `K_{ũũ}`). Regla práctica: `K_e ≈ 10·E_bulk/ℓ_c` con `ℓ_c` la longitud característica del elemento (su lado más corto). El usuario lo declara explícitamente en YAML; no hay default automático.
- **Cierre en compresión**. Cuando `[[u_n]] < 0` el modelo devuelve `t_n = (1−ω)·K_e·[[u_n]]` (compresión reducida). Físicamente, la grieta cerrada debería transmitir compresión a rigidez completa `K_e` (contacto unilateral). En fase 1 esta limitación se acepta tal cual — los benchmarks de la tesis (Van Vliet, viga SEN) son traccionados puros. Un modelo con corrección de cierre se contempla como **deuda explícita** para casos con cambio de signo en la grieta.
- **Saturación numérica**. `DAMAGE_MAX < 1` evita rigidez exactamente nula (singularidad del Newton). Constante del proyecto `fenix.constants.DAMAGE_MAX` (mismo valor que `IsotropicDamage2D`).
- **Sensibilidad de malla** en régimen de softening. Sin regularización, el patrón de fallo localiza en función del tamaño y orientación de elementos. La integración con el elemento `CST_Embedded2D` (fase 2 del ADR 0010) y con el `l_d = (A/h)·cos(θ−α)` de Cap. 6 da objetividad parcial respecto a la longitud de la grieta — pero la dirección de propagación sigue dependiendo de la dirección de tensiones principales en el elemento. Mitigación a nivel del *elemento*, no del material.
- **Validación energética**. Por construcción `∫_0^{w_c} t(w)·dw = G_F`. Esto se verifica en los tests (§acceptance). Si la curva implementada no integra a `G_F`, hay error en `w_c` o en `H`.

---

## Contrato de implementación

```yaml
name: CohesiveDamageIsotropic
kind: cohesive_material         # nueva familia, ADR 0010
status: draft

interface:
  jump_dim: 2                   # [[u]] ∈ ℝ² en 2D (componentes n, s en frame local)
  primary_state_var: damage     # ω, escalar exportable al post-proceso
  is_symmetric: true            # tangente simétrica en Modo-I (ver §10)

parameters:
  - { name: sigma_t0,           type: float, required: true,  desc: "Resistencia a tracción [Pa]; iniciación del daño" }
  - { name: G_f,                type: float, required: true,  desc: "Energía de fractura [N/m]; área bajo la curva t–[[u_n]]" }
  - { name: K_e,                type: float, required: true,  desc: "Rigidez elástica del salto en dirección normal [Pa/m] (penalty); ver §12" }
  - { name: softening,          type: str,   required: true,  desc: "Forma de la curva de softening: 'linear' | 'exponential'" }

signature:
  compute_traction: "([[u]]: ndarray(2,), state_vars=None) -> (t: ndarray(2,), T_tan: ndarray(2,2), state_vars')"
  jump_kind: "frame local — [[u]] = ([[u_n]], [[u_s]]) en ejes (n, s) de Γ_d"

state_schema:
  kappa: "float ≥ κ_0 — historial del salto equivalente máximo"
  damage: "float ∈ [0, DAMAGE_MAX] — variable de daño escalar"

conventions:
  sign: "[[u_n]] > 0 ⇔ apertura (tracción); el modelo daña sólo en apertura (Modo-I)"
  frame: "ejes locales (n, s) de Γ_d; transformación a globales por el elemento consumidor"
  damage_irreversibility: "κ monótono no decreciente vía max(κ_old, ⟨[[u_n]]⟩); ω(κ) monótono creciente"
  units: "[σ_t0] = Pa, [G_F] = N/m, [K_e] = Pa/m (asumiendo SI; el sistema no convierte unidades)"

validity:
  - "σ_t0 > 0, G_F > 0, K_e > 0"
  - "K_e > σ_t0² / (2·G_F)  — consistencia energética para softening exponencial (§5)"
  - "softening ∈ {'linear', 'exponential'}"

out_of_scope:
  - "Modo mixto I-II (cizalla acoplada al daño); diferido a fase G del ADR 0010"
  - "Contacto / cierre unilateral en compresión; el modelo reduce rigidez al estar dañado, no recupera contacto rígido"
  - "Anisotropía del daño (tensor de daño D vs escalar ω)"
  - "Acoplamiento viscoso (rate-dependent); rate-independent puro en fase 1"
  - "Fatiga / acumulación cíclica del daño"
  - "Regularización para mesh-objectivity en propagación; se aborda a nivel de elemento (l_d correcto, Cap. 6) y no de material"

numerical_caveats:
  - "Penalty K_e: balance entre κ_0 despreciable (K_e grande) y condicionamiento (K_e moderado). Guía: K_e ≈ 10·E_bulk/ℓ_c. Declarado por el usuario en YAML, sin default automático."
  - "Cierre en compresión: rigidez reducida (1−ω)·K_e, no recuperación de contacto rígido. Aceptable para benchmarks puramente traccionados; deuda explícita para casos con cambio de signo."
  - "Saturación ω = DAMAGE_MAX evita singularidad numérica; la tangente en saturación es secante pura sin corrección consistente."
  - "Régimen elastic↔softening discontinuo en la tangente al cruzar κ_0; régimen loading↔unloading discontinuo al cruzar κ_n committed. Mismo patrón que IsotropicDamage2D."
  - "Para softening lineal con κ → w_c, ω → 1 con derivada finita: la transición a ω=DAMAGE_MAX puede dejar un residuo pequeño en la curva (corte numérico). Aceptable si DAMAGE_MAX está cerca de 1."

acceptance:
  verification:
    - name: respuesta_elastica_bajo_umbral
      setup: "carga monótona en [[u_n]] ∈ [0, κ_0); [[u_s]] = 0"
      expect: "ω = 0; t_n = K_e·[[u_n]] exacto; t_s = 0; T_tan = K_e·(n⊗n)"
      tol_rel: 1.0e-12

    - name: inicio_dano_en_umbral
      setup: "[[u_n]] cruza κ_0 = σ_t0/K_e por encima por primera vez"
      expect: "ω salta de 0 a un valor positivo continuo en κ_0 (ω(κ_0) = 0+); t_n(κ_0) = σ_t0 exacto"
      tol_rel: 1.0e-10

    - name: energia_fractura_softening_lineal
      setup: "carga monótona [[u_n]] de 0 a w_c = 2·G_F/σ_t0 con softening='linear', K_e tal que κ_0 ≪ w_c"
      expect: "∫_0^{w_c} t_n([[u_n]])·d[[u_n]] = G_F  (cuadratura trapezoidal con ≥1000 puntos)"
      tol_rel: 1.0e-3
      note: "La tolerancia 1e-3 absorbe el error de cuadratura; la consistencia analítica es exacta por construcción de w_c."

    - name: energia_fractura_softening_exponencial
      setup: "carga monótona [[u_n]] de 0 a 20·G_F/σ_t0 con softening='exponential', K_e tal que κ_0 ≪ G_F/σ_t0"
      expect: "∫ t_n·d[[u_n]] ≈ G_F·(1 − ε_truncamiento), con ε_truncamiento < 1e-6 por la cola exponencial cortada"
      tol_rel: 1.0e-3

    - name: descarga_secante
      setup: "carga hasta κ ∈ (κ_0, w_c/2), después descarga [[u_n]] → 0"
      expect: "rama de descarga lineal con pendiente (1−ω(κ))·K_e; al volver a [[u_n]]=0, t_n = 0 y ω no ha cambiado"
      tol_rel: 1.0e-10

    - name: recarga_hasta_mismo_kappa_sin_dano_extra
      setup: "carga → descarga → recarga hasta [[u_n]] = κ (mismo κ histórico)"
      expect: "ω no cambia (κ_{n+1} = max(κ_n, [[u_n]]) = κ_n)"
      tol_rel: 1.0e-12

    - name: recarga_mas_alla_continua_danando
      setup: "carga → descarga → recarga hasta [[u_n]] > κ histórico"
      expect: "ω evoluciona desde el valor previo hacia ω([[u_n]]_nuevo) sin saltos"
      tol_rel: 1.0e-10

  specific:
    - name: tangente_simetrica_en_carga_y_descarga
      setup: "estados de carga activa, descarga, sin daño, saturado"
      expect: "T_tan = T_tan^T (simetría exacta en Modo-I); ‖T_tan − T_tan^T‖_F / ‖T_tan‖_F < 1e-14"
      tol_rel: 1.0e-14

    - name: tangente_diferencia_finita_consistente
      setup: "carga activa con [[u_n]] ∈ (κ_0, w_c/2); perturbar [[u]] con ε = 1e-7"
      expect: "‖T_tan_FD − T_tan_analitica‖ / ‖T_tan‖ < 1e-5"
      tol_rel: 1.0e-5

    - name: tangencial_no_dana_ni_transmite
      setup: "[[u_s]] arbitrario, [[u_n]] = 0"
      expect: "t_s = 0; t_n = 0; ω no cambia"
      tol_rel: 1.0e-14

    - name: compresion_no_dana
      setup: "[[u_n]] < 0 puro (penetración)"
      expect: "⟨[[u_n]]⟩ = 0 ⇒ f = −κ ≤ 0 ⇒ ω no cambia. t_n = (1−ω_previo)·K_e·[[u_n]] (caveat de cierre, §12)"
      tol_rel: 1.0e-12

    - name: saturacion_en_w_c_softening_lineal
      setup: "[[u_n]] = w_c·1.01 con softening='linear'"
      expect: "ω = DAMAGE_MAX (saturado); t_n ≈ (1−DAMAGE_MAX)·K_e·[[u_n]] (residuo numérico controlado)"
      tol_rel: 1.0e-10

    - name: degeneracion_a_elasticidad_intacta
      setup: "σ_t0 → ∞ (10^20 Pa) o cualquier [[u_n]] tan pequeño que jamás active"
      expect: "ω = 0, t_n = K_e·[[u_n]], T_tan = K_e·(n⊗n) (material intacto puro)"
      tol_rel: 1.0e-12

  arch:
    - name: registry_kind_es_cohesive
      expect: "Registro vía CohesiveMaterialRegistry, no MaterialRegistry; kind='cohesive_material' en YAML"

    - name: is_symmetric_attribute
      expect: "CohesiveDamageIsotropic.IS_SYMMETRIC == True"

references:
  - "Retama Velasco, J. (2010). Formulation and Approximation to Problems in Solids by Embedded Discontinuity Models. Tesis Doctoral, UNAM. **Cap. 3 — Discrete damage models**. Fórmulas §3.1 (energía libre, ec. 3.1–3.4); §3.1.1 (tensor tangente, ec. 3.7–3.12); §3.1.2 (función de fluencia, ec. 3.13–3.14); §3.1.3 (Kuhn-Tucker, ec. 3.15–3.16)."
  - "Hillerborg, A., Modéer, M., Petersson, P.-E. (1976). Analysis of crack formation and crack growth in concrete by means of fracture mechanics and finite elements. *Cement and Concrete Research* 6, 773-782. **Origen del modelo cohesivo con energía de fractura G_F**."
  - "Barenblatt, G.I. (1962). The mathematical theory of equilibrium cracks in brittle fracture. *Advances in Applied Mechanics* 7, 55-129."
  - "Simó, J.C., Ju, J.W. (1987). Strain- and stress-based continuum damage models — I. Formulation. *Int. J. Solids Struct.* 23, 821-840."
  - "ADR 0010 — Discontinuidades interiores embebidas (familia CohesiveMaterial, hoja de ruta)."
  - "Spec [IsotropicDamage2D](IsotropicDamage2D.md) — patrón análogo de daño isótropo escalar para el caso continuo."
```

---

## Implementación

*Rellena la IA tras programar.*

- Archivo: —
- Clase: —
- Tests:
  - —
- Notas de traducción: —

---

## Diálogo

- **2026-05-13** · Spec pre-redactada por la IA como **orden de trabajo** sobre el Cap. 3 de Retama 2010 (tesis del usuario). La física es traducción directa de la tesis; la mecánica del proyecto (Voigt, registry, signatura, tests) sigue convenciones de `IsotropicDamage2D`. Pendiente de validación del usuario.
- **Punto abierto** · Decisión sobre **PRIMARY_STATE_VAR** (`damage` vs `kappa`). En el borrador se adoptó `damage` por interpretabilidad física (mismo criterio que `IsotropicDamage2D`); si el usuario prefiere `kappa` por ser la variable más fundamental, cambia trivialmente.
- **Punto abierto** · La fórmula cerrada de `dω/dκ` para softening lineal se dejó simplificable durante la implementación. Sale directa de derivar `ω(κ) = 1 − [σ_t0·(1 − (κ−κ_0)/(w_c−κ_0))]/(K_e·κ)` con regla del cociente — el usuario revisa el álgebra al cerrar el docstring de la clase.
- **Punto abierto** · `K_e` queda como parámetro obligatorio sin default automático. Alternativa considerada: derivar `K_e` de `E_bulk` y `ℓ_c` automáticamente al construir el elemento. Se descartó para no acoplar bulk y cohesivo en el contrato del material — el bulk es responsabilidad del elemento, el cohesivo se construye independiente y se acopla en el elemento. Si el usuario prefiere lo contrario, se revisa.
