# DRUCKER-PRAGER 2D — plasticidad friccional con endurecimiento isótropo lineal

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general

Modelo elastoplástico **friccional** independiente de la velocidad. Suelos, hormigón, granulares, materiales cuasi-frágiles cohesivo-friccionales. Suaviza la pirámide hexagonal de Mohr-Coulomb a un **cono circular** en el espacio de esfuerzos principales, eliminando las aristas que complican el return mapping.

Diferencia clave con J2:

- **Dependencia de la presión**. La resistencia al corte aumenta con la presión de confinamiento (componente friccional).
- **Plasticidad no asociada**. La dilatancia volumétrica real del flujo plástico ($\psi$) es típicamente menor que la fricción del criterio ($\phi$). Forzar $\psi = \phi$ (asociada) sobreestima la dilatancia y el volumen post-fluencia.
- **Retorno al ápice**. El criterio define un cono con vértice; estados de prueba en zona puramente hidrostática extrema retornan al vértice (zona singular en la dirección de flujo).

Hipótesis cinemática soportada: **`plane_strain`** únicamente en esta spec.

### 1. Descomposición aditiva y elasticidad

Pequeñas deformaciones: $\boldsymbol\varepsilon = \boldsymbol\varepsilon^e + \boldsymbol\varepsilon^p$. Ley elástica isótropa: $\boldsymbol\sigma = \mathbf C_e (\boldsymbol\varepsilon - \boldsymbol\varepsilon^p)$.

### 2. Criterio de fluencia (Drucker-Prager)

$$f(\boldsymbol\sigma, \alpha) = \sqrt{J_2} + \eta_f\,I_1 - k(\alpha) \;\leq\; 0$$

con:

- $J_2 = \tfrac{1}{2}\,\mathbf s : \mathbf s$, segundo invariante del desviador.
- $I_1 = \operatorname{tr}(\boldsymbol\sigma) = \sigma_{xx} + \sigma_{yy} + \sigma_{zz}$, primer invariante.
- $\eta_f$, coeficiente friccional del cono (adimensional, $\geq 0$).
- $k(\alpha) = k_0 + H_k\,\alpha$, radio del cono (el parámetro de cohesión del criterio: $k_0$ sale de $c_0$ según la calibración, §3) con endurecimiento isótropo lineal en el multiplicador plástico acumulado $\alpha$ (§5).

La superficie es un cono con eje en la línea hidrostática ($\sigma_{xx} = \sigma_{yy} = \sigma_{zz}$), vértice en $I_1 = k/\eta_f$, abertura controlada por $\eta_f$.

### 3. Calibración con Mohr-Coulomb

Los parámetros $(\eta_f, k_0)$ se derivan de los parámetros físicos $(c_0, \phi)$ — cohesión y ángulo de fricción interna — según la variante del cono elegida:

| variant | $\eta_f$ | $k_0$ |
|---|---|---|
| `plane_strain_matched` (default) | $\tan(\phi)/\sqrt{9 + 12\tan^2(\phi)}$ | $3\,c_0/\sqrt{9 + 12\tan^2(\phi)}$ |
| `outer_cone` | $2\sin(\phi)/[\sqrt 3\,(3 - \sin(\phi))]$ | $6\,c_0\cos(\phi)/[\sqrt 3\,(3 - \sin(\phi))]$ |
| `inner_cone` | $2\sin(\phi)/[\sqrt 3\,(3 + \sin(\phi))]$ | $6\,c_0\cos(\phi)/[\sqrt 3\,(3 + \sin(\phi))]$ |

`plane_strain_matched` es la variante que coincide exactamente con Mohr-Coulomb en condiciones de plane strain (Drucker, 1953). Las otras dos circunscriben (outer) o inscriben (inner) el hexágono de MC en el plano π. Como esta spec solo soporta plane strain, `plane_strain_matched` es el default natural.

### 4. Regla de flujo (no asociada por defecto)

Potencial plástico $g$ con coeficiente de **dilatancia** $\eta_g$ (función del ángulo de dilatancia $\psi$, calibrado con la misma fórmula que $\eta_f$ pero usando $\psi$):

$$g(\boldsymbol\sigma) = \sqrt{J_2} + \eta_g\,I_1$$

$$\dot{\boldsymbol\varepsilon}^p = \dot\gamma\,\frac{\partial g}{\partial \boldsymbol\sigma} = \dot\gamma\,\left(\frac{\mathbf s}{2\sqrt{J_2}} + \eta_g\,\mathbf I\right)$$

Descomposición:
- **Parte desviadora**: $\dot{\mathbf e}^p = \dot\gamma\,\mathbf s/(2\sqrt{J_2})$ — análoga a J2.
- **Parte volumétrica**: $\operatorname{tr}(\dot{\boldsymbol\varepsilon}^p) = 3\dot\gamma\,\eta_g$ — **dilatación plástica positiva** si $\eta_g > 0$.

Si $\psi = \phi$ (i.e. $\eta_g = \eta_f$) el modelo se vuelve **asociado**; tangente algorítmica simétrica. Para $\psi < \phi$ (típico en suelos: $\psi \approx \phi/2$ o incluso $\psi = 0$, "non-dilatant"), el modelo es **no asociado** y la tangente algorítmica es **asimétrica**.

### 5. Endurecimiento isótropo lineal

$$\dot\alpha = \dot\gamma$$

con $\alpha$ el **multiplicador plástico acumulado** $\gamma$ (convención del proyecto para Drucker-Prager). La superficie crece manteniendo eje y abertura, solo $k$ aumenta linealmente.

**Parámetro del usuario (desde el 2026-09-29): la parametrización de de Souza Neto, Perić y Owen (2008, §8.3).** La cohesión se endurece con la deformación plástica equivalente, $c(\bar\varepsilon^p) = c_0 + H\,\bar\varepsilon^p$, con $\bar\varepsilon^p = \xi\,\gamma$ y $k = \xi\,c$, donde $\xi = k_0/c_0$ es el factor de la calibración (§3). El parámetro `H` que escribe el usuario es $H = H_c = dc/d\bar\varepsilon^p$. El return mapping trabaja con el módulo del radio del cono respecto del multiplicador acumulado, $H_k = dk/d\alpha$, que es

$$H_k = \xi^2\,H_c .$$

Endurecimiento por fricción/dilatancia (cambio de $\eta_f$ o $\eta_g$ con $\alpha$) queda fuera de esta spec — escenario que no necesitamos hoy y que añade no-linealidad fuerte al return mapping.

### 6. Condiciones de Kuhn-Tucker

$$\dot\gamma \geq 0,\qquad f \leq 0,\qquad \dot\gamma\,f = 0$$

### 7. Variables internas

| Símbolo | Tipo | Significado |
|---|---|---|
| $\boldsymbol\varepsilon^p$ | tensor 4 componentes `[xx, yy, zz, xy_tensorial]` | deformación plástica |
| $\alpha$ | escalar $\geq 0$ | multiplicador plástico acumulado (= acumulado de $\Delta\gamma$; no es una deformación plástica equivalente, §5) |

`alpha` es la variable principal exportable (`PRIMARY_STATE_VAR = 'alpha'`).

---

## Formulación numérica

### 8. Esquema implícito — Backward Euler

Dado $\{\boldsymbol\varepsilon^p_n, \alpha_n\}$ y $\boldsymbol\varepsilon_{n+1}$, predictor elástico:

$$\boldsymbol\varepsilon^\text{trial} = \boldsymbol\varepsilon_{n+1} - \boldsymbol\varepsilon^p_n, \qquad \mathbf s^\text{trial} = 2G\,\mathbf e^\text{trial}_\text{dev},\qquad p^\text{trial} = K\,\operatorname{tr}(\boldsymbol\varepsilon^\text{trial})$$

con $I_1^\text{trial} = 3\,p^\text{trial}$, $\sqrt{J_2^\text{trial}} = \|\mathbf s^\text{trial}\|/\sqrt{2}$.

Función de fluencia trial:

$$f^\text{trial} = \sqrt{J_2^\text{trial}} + \eta_f\,I_1^\text{trial} - k(\alpha_n)$$

Si $f^\text{trial} \leq \text{tol}$ → paso elástico. Si no, intentar **return regular**.

### 9. Return regular (cone surface)

La dirección desviadora $\hat{\mathbf n} = \mathbf s^\text{trial}/(2\sqrt{J_2^\text{trial}}) = \mathbf s^{n+1}/(2\sqrt{J_2^{n+1}})$ es invariante bajo carga radial ($\hat{\mathbf n} = \partial\sqrt{J_2}/\partial\boldsymbol\sigma$, de norma tensorial $1/\sqrt 2$; es la misma que usa la tangente de §12). Las actualizaciones son lineales en $\Delta\gamma$:

$$\sqrt{J_2^{n+1}} = \sqrt{J_2^\text{trial}} - G\,\Delta\gamma$$
$$p^{n+1} = p^\text{trial} - 3K\,\eta_g\,\Delta\gamma, \qquad I_1^{n+1} = 3\,p^{n+1}$$

La consistencia $f(\sigma_{n+1}, \alpha_{n+1}) = 0$ con $\alpha_{n+1} = \alpha_n + \Delta\gamma$ y $k_{n+1} = k_0 + H_k(\alpha_n + \Delta\gamma)$ da:

$$\boxed{\;\Delta\gamma = \frac{f^\text{trial}}{G + 9K\,\eta_f\,\eta_g + H_k}\;}$$

**Forma cerrada** — sin Newton local. Cada término del denominador:
- $G$: cambio en $\sqrt{J_2}$ por la parte desviadora del flujo.
- $9K\,\eta_f\,\eta_g$: cambio en $\eta_f I_1$ por la parte volumétrica dilatante del flujo (multiplicado por $\eta_f$ del criterio y $\eta_g$ del potencial).
- $H_k$: endurecimiento del radio del cono, $dk/d\alpha$ (§5).

**Actualizaciones**:

$$\mathbf s^{n+1} = \frac{\sqrt{J_2^{n+1}}}{\sqrt{J_2^\text{trial}}}\,\mathbf s^\text{trial}$$

$$\boldsymbol\sigma^{n+1} = \mathbf s^{n+1} + p^{n+1}\,\mathbf I$$

$$\boldsymbol\varepsilon^p_{n+1} = \boldsymbol\varepsilon^p_n + \Delta\gamma\,\left(\frac{\mathbf s^{n+1}}{2\sqrt{J_2^{n+1}}} + \eta_g\,\mathbf I\right)$$

**Condición de validez del return regular**: $\sqrt{J_2^{n+1}} \geq 0$. Si $\Delta\gamma_\text{reg} > \sqrt{J_2^\text{trial}}/G$, el desviador "se pasa de cero" y el algoritmo regular es inválido — toca **return al ápice**.

### 10. Return al ápice

Cuando el estado de prueba cae en zona puramente hidrostática extrema (p_trial muy traccionante respecto a la cohesión), la proyección al cono cae en su vértice, donde $\mathbf s = 0$ y la dirección desviadora $\hat{\mathbf n}$ no está definida.

En el ápice: $\boldsymbol\sigma^{n+1} = \bigl(k_{n+1}/(3\eta_f)\bigr)\,\mathbf I$ (puramente hidrostático; $I_1^{n+1} = k_{n+1}/\eta_f$ es el vértice del cono de §2), $\sqrt{J_2^{n+1}} = 0$.

Por conservación volumétrica:

$$p^{n+1} = p^\text{trial} - 3K\,\eta_g\,\Delta\gamma = \frac{k_{n+1}}{3\,\eta_f}$$

con $k_{n+1} = k_0 + H_k(\alpha_n + \Delta\gamma)$. Multiplicando por $3\eta_f$:

$$3\eta_f\,p^\text{trial} - 9K\,\eta_f\,\eta_g\,\Delta\gamma = k(\alpha_n) + H_k\,\Delta\gamma$$

y, con $I_1^\text{trial} = 3\,p^\text{trial}$, despejando:

$$\boxed{\;\Delta\gamma_\text{apex} = \frac{I_1^\text{trial}\,\eta_f - k(\alpha_n)}{9K\,\eta_f\,\eta_g + H_k}\;}$$

Es la fórmula regular sin la parte desviadora, porque la rama de apex es puramente volumétrica: el numerador es $f^\text{trial} - \sqrt{J_2^\text{trial}}$ y el denominador pierde el término $G$.

**Actualizaciones (apex)**:

$$\mathbf s^{n+1} = 0,\qquad p^{n+1} = p^\text{trial} - 3K\,\eta_g\,\Delta\gamma_\text{apex} \;\equiv\; \frac{k_{n+1}}{3\,\eta_f}, \qquad p^\text{trial} = K\operatorname{tr}(\boldsymbol\varepsilon_{n+1} - \boldsymbol\varepsilon^p_n)$$

(el código evalúa directamente $p^{n+1} = k_{n+1}/(3\eta_f)$.)

$$\boldsymbol\sigma^{n+1} = p^{n+1}\,\mathbf I$$

$$\boldsymbol\varepsilon^p_{n+1} = \boldsymbol\varepsilon^p_n + \Delta\gamma_\text{apex}\,\eta_g\,\mathbf I + \mathbf e^\text{trial}_\text{dev}$$

con $\mathbf e^\text{trial}_\text{dev}$ el desviador de la deformación de prueba de §8 (la parte desviadora plástica absorbe completamente el desviador de prueba, ya que $\mathbf s^{n+1} = 0$).

### 11. Detección del régimen de retorno

Algoritmo de despacho:

1. Calcular $\Delta\gamma_\text{reg}$ por la fórmula regular.
2. Si $\sqrt{J_2^\text{trial}} - G\,\Delta\gamma_\text{reg} \geq 0$ → return regular.
3. Si no → return al ápice con $\Delta\gamma_\text{apex}$.

El criterio es geométricamente: el predictor cae más allá del ápice si la proyección desviadora exigiría reducir $\sqrt{J_2}$ por debajo de cero.

### 12. Tangente algorítmica consistente

**Return regular** ($\sqrt{J_2^{n+1}} > 0$). Linealización del algoritmo de §9, forma compacta (de Souza Neto §8.3.4) en Voigt 3 plane strain; es la que implementa el código y la verifica `tests/test_consistent_tangent.py` por diferencias finitas:

$$\mathbf C^\text{alg} = K\,\mathbf v\otimes\mathbf v + 2G\,(1 - \beta)\,\mathbf I_\text{dev} + 4G\,\beta\,\hat{\mathbf n}\otimes\hat{\mathbf n} - \frac{1}{A}\,\mathbf b_g \otimes \mathbf b_f$$

con:

- $\mathbf v = [1, 1, 0]$ (gradiente de $I_1$ respecto de la deformación en Voigt 3 *engineering*).
- $\mathbf I_\text{dev}$: bloque en plano $3\times 3$ del proyector desviador, $\begin{bmatrix} 2/3 & -1/3 & 0 \\ -1/3 & 2/3 & 0 \\ 0 & 0 & 1/2 \end{bmatrix}$ (el $1/2$ del cortante lleva $\gamma_{xy}$ a $\varepsilon_{xy}$ tensorial).
- $\hat{\mathbf n} = [s_{xx}, s_{yy}, s_{xy}]^\text{trial}/(2\sqrt{J_2^\text{trial}})$, con el cortante tensorial (sin factor 2); coincide con $\mathbf s^{n+1}/(2\sqrt{J_2^{n+1}})$ (§9).
- $\beta = G\,\Delta\gamma/\sqrt{J_2^\text{trial}}$.
- $\mathbf b_f = 2G\,\hat{\mathbf n} + 3K\,\eta_f\,\mathbf v$, $\mathbf b_g = 2G\,\hat{\mathbf n} + 3K\,\eta_g\,\mathbf v$.
- $A = G + 9K\,\eta_f\,\eta_g + H_k$ (denominador de $\Delta\gamma$, §9).

El término $+4G\beta\,\hat{\mathbf n}\otimes\hat{\mathbf n}$ recoge el cambio de la dirección de flujo con $\boldsymbol\varepsilon$ a través de $\mathbf s^\text{trial}$ (no sólo el de su magnitud); sin él la tangente cierra mal en las componentes de cortante alineadas con el flujo (lo detecta la diferencia finita). Con $\mathbf N = \sqrt 2\,\hat{\mathbf n}$ (normal unitaria del desviador) es $2G\beta\,\mathbf N\otimes\mathbf N$, y el conjunto es la forma de de Souza Neto et al. (2008, §8.3) con $\eta = 3\eta_f$, $\bar\eta = 3\eta_g$ (su criterio usa $p$, éste $I_1 = 3p$) y $\xi^2 H_c = H_k$. **Asimétrica si $\eta_f \neq \eta_g$** (no asociada): $\mathbf b_g \otimes \mathbf b_f$ sólo es simétrico con $\mathbf b_g = \mathbf b_f$.

**Return al ápice** ($\sqrt{J_2^{n+1}} = 0$):

$$\mathbf C^\text{alg} = K\,\Big(1 - \frac{3K\eta_f\eta_g}{3K\eta_f\eta_g + H_k/3}\Big)\,\mathbf v\otimes\mathbf v$$

— rigidez puramente volumétrica reducida (el material en el ápice se comporta como un fluido con compresibilidad finita; sin rigidez desviadora alguna), con $\mathbf v = [1, 1, 0]$ como en el return regular.

Forma equivalente más limpia: $\mathbf C^\text{alg}_\text{apex} = \dfrac{K\,H_k}{9K\eta_f\eta_g + H_k}\,\mathbf v\otimes\mathbf v$ (con $H_k=0$, el ápice se vuelve completamente plástico volumétrico y la tangente colapsa).

**Piso de la tangente en el ápice.** Si $K H_k/(9K\eta_f\eta_g + H_k) < 10^{-6} K$, el código devuelve $10^{-6} K\,\mathbf v\otimes\mathbf v$ en lugar de la tangente exacta, para no dejar la matriz global singular mientras el punto siga en el ápice. Actúa en **todo retorno al ápice con $H_k = 0$** (perfectamente plástico, el caso más común, donde la tangente exacta es nula) y, en general, siempre que $H_k/(9K\eta_f\eta_g + H_k) < 10^{-6}$; con $\psi = 0$ el retorno al ápice no tiene solución y el kernel lanza `ValueError` antes de llegar aquí (Diálogo, 2026-09-22). Es una regularización de la tangente, no la tangente consistente: no cambia $\boldsymbol\sigma$ ni el estado convergido, sólo la convergencia del Newton global.

### 13. Tolerancia de fluencia

ADR 0006: `admissibility_scale = k(α) = k_0 + H_k·α` (radio corriente del cono). Tiene unidades de esfuerzo, igual que $f$.

---

## Contrato de implementación

```yaml
name: DruckerPrager2D
kind: material
status: validated

interface:
  strain_dim: 3
  primary_state_var: alpha
  is_symmetric: false      # tangente asimétrica si no asociada (ψ ≠ φ)

admission:                # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - de_souza_neto: "§8.3 (Drucker-Prager: retorno al cono y al ápice, tangentes)"
  general_codes:
    - abaqus: "*DRUCKER PRAGER (lineal) en CPE*"
    - code_aster: "DRUCK_PRAGER / DRUCK_PRAG_N_A"
    - opensees: "nDMaterial DruckerPrager"

parameters:
  - { name: E,          type: float, required: true,  desc: "Módulo de Young" }
  - { name: nu,         type: float, required: true,  desc: "Coeficiente de Poisson" }
  - { name: cohesion,   type: float, required: true,  desc: "Cohesión c_0 inicial (esfuerzo)" }
  - { name: phi_deg,    type: float, required: true,  desc: "Ángulo de fricción interna en grados; típico suelos 20-40°, hormigón 30-37°" }
  - { name: psi_deg,    type: float, required: false, default: null, desc: "Ángulo de dilatancia en grados. Default = phi_deg (asociada). Típico no asociada: ψ < φ, e.g. 0 o φ/2." }
  - { name: H,          type: float, required: false, default: 0.0, desc: "Módulo de endurecimiento isótropo lineal de la cohesión H = dc/dε̄ᵖ (≥0), parametrización de de Souza Neto et al. (2008, §8.3); internamente H_k = ξ²·H con ξ = k_0/c_0 (§5). H=0 ⇒ perfectamente plástico." }
  - { name: hypothesis, type: str,   required: false, default: "plane_strain", desc: "plane_strain (único soportado)" }
  - { name: variant,    type: str,   required: false, default: "plane_strain_matched", desc: "plane_strain_matched | outer_cone | inner_cone — calibración con Mohr-Coulomb" }
  - { name: density,    type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(3,), state_vars=None) -> (σ: ndarray(3,), C_alg: ndarray(3,3), state_vars')"
  strain_kind: "plano Voigt — [ε_xx, ε_yy, γ_xy] con γ_xy = 2·ε_xy"

state_schema:
  eps_p:   "ndarray(4,) — [ε^p_xx, ε^p_yy, ε^p_zz, ε^p_xy_tensorial]"
  alpha:   "float ≥ 0 — multiplicador plástico acumulado"

conventions:
  sign: "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "γ_xy = 2·ε_xy en deformación; ε^p_xy en tensorial (sin factor 2)"
  associated_flow: "ψ = φ ⇒ asociada (tangente simétrica); ψ ≠ φ ⇒ no asociada"

validity:
  - "E > 0, c_0 > 0, 0 ≤ φ < 90°, 0 ≤ ψ ≤ φ"
  - "ν ∈ (-1, 0.5)"
  - "H ≥ 0"
  - "hypothesis == 'plane_strain' (único soportado en esta spec)"
  - "variant ∈ {'plane_strain_matched', 'outer_cone', 'inner_cone'}"

out_of_scope:
  - "plane_stress — la proyección de Drucker-Prager con σ_zz=0 acoplada a regla de flujo dilatante es notoriamente delicada; difiere hasta caso real"
  - "endurecimiento por cambio de φ y/o ψ con α — no lineal fuerte en el return mapping"
  - "endurecimiento cinemático (back-stress hidrostático y desviador)"
  - "cap (modelo cap-cone) para suelos compactables"
  - "ablandamiento (H<0) — requiere regularización"
  - "grandes deformaciones"

numerical_caveats:
  - "Tangente algorítmica asimétrica si ψ ≠ φ → IS_SYMMETRIC=False, el despachador algebraico usa LU (ADR 0003)."
  - "Modos de carga muy traccionantes pueden caer en la rama de ápice. La detección automática del régimen (regular vs apex) está en el algoritmo; si en el futuro aparece oscilación entre ramas en pasos cerca de la transición, considerar smoothing con vertex regularization."
  - "Con H_k = 0 (perfectamente plástico), en la rama de ápice la tangente volumétrica exacta se anula → singularidad. El código la sustituye por el piso K·1e-6·1⊗1 en todo retorno al ápice con H_k/(9Kη_fη_g + H_k) < 1e-6 (§12): regularización de la tangente, sin efecto sobre σ; típicamente sale del ápice en el siguiente paso. Con ψ = 0 el retorno al ápice no tiene solución y el kernel lanza ValueError."
  - "El return regular puede generar α creciente con incrementos pequeños cerca de la transición regular↔apex. Newton global converge igual; si en problema concreto se detecta oscilación, usar arc-length."

acceptance:
  unit:
    - name: parameter_validation
      expect: "constructor rechaza variant inválida, φ ≥ 90°, ψ > φ, hypothesis distinto de plane_strain"

    - name: calibration_matched_consistency
      expect: "para variant='plane_strain_matched' con ψ = φ, asociada → η_f = η_g; con ψ = 0, η_g = 0 (sin dilatación)"

    - name: paso_elastico_below_yield
      expect: "ε pequeño; σ = C_e·ε, eps_p y alpha intactos"
      tol_rel: 1.0e-10

    - name: regular_return_under_shear
      setup: "ε cortante puro suficiente para fluir (ε_xy = γ/2, sin parte hidrostática)"
      expect: "return regular activo, α > 0, √J₂ del estado final = k(α) − η_f·I₁ ≥ 0 (f = 0, cumple criterio)"
      tol_rel: 1.0e-8

    - name: apex_return_under_pure_tension
      setup: "ε puramente hidrostática traccionante (ε_xx = ε_yy > 0, sin cortante)"
      expect: "return al ápice; σ_xx = σ_yy = k/(3η_f) tras converger; √J₂ = 0"
      tol_rel: 1.0e-8

    - name: associated_tangent_is_symmetric
      setup: "ψ = φ (asociada), estado plástico activo"
      expect: "C_alg simétrica (‖C_alg − C_alg^T‖ < tol)"
      tol_abs: 1.0e-8

    - name: nonassociated_tangent_is_asymmetric
      setup: "ψ = 0, φ = 30°, estado plástico activo no degenerado"
      expect: "C_alg − C_alg^T ≠ 0 con norma significativa"

    - name: tangent_matches_finite_difference
      setup: "estado plástico regular; comparar C_alg cerrada con diferencia finita centrada"
      expect: "error relativo < 1e-5"
      tol_rel: 1.0e-5

    - name: alpha_monotonic_under_increasing_load
      expect: "bajo carga creciente α no decrece"

    - name: unloading_recovers_C_e
      setup: "cargar a plástico, descargar elásticamente"
      expect: "tangente = C_e en el paso de descarga, α se conserva"
      tol_rel: 1.0e-10

    - name: unit_invariance_MPa_vs_Pa
      expect: "mismo path adimensional en (MPa,mm,N) y (Pa,m,N) ⇒ α idéntico"
      tol_rel: 1.0e-12

  integration:
    - name: quad4_pipeline_elastico
      setup: "Quad4 unitario plane strain, carga muy por debajo del yield"
      expect: "respuesta elástica pura, α=0"
      tol_rel: 1.0e-8

    - name: quad4_pipeline_plastico_converges
      setup: "Quad4 unitario plane strain, carga que produce plasticidad activa, 4 pasos, max_iter=8"
      expect: "Newton converge en cada paso, α > 0 al final"

    - name: quad4_yield_onset_confined_tension          # DP-1 (2026-05-19)
      setup: "Quad4 unitario plane strain con desplazamientos prescritos en los 4 nodos para ε_xx=k·ε_y, ε_yy=0, γ_xy=0"
      expect: "ε_xx=0.5·ε_y elástico (σ_xx coincide con (λ+2μ)·ε exacto, α=0); ε_xx=0.99·ε_y aún α=0; ε_xx=1.5·ε_y plástico con α>0 idéntico en los 4 Gauss points"
      tol_rel: 1.0e-8

    - name: quad4_flow_rule_invariant                   # DP-2 (2026-05-19)
      setup: "Quad4 unitario plane strain, carga monótona en rama regular (5 niveles 1.2·ε_y → 2.6·ε_y), ψ ≠ φ"
      expect: "invariante cinemática tr(ε_p) = 3·η_g·α en cada Gauss point y en cada nivel"
      tol_rel: 1.0e-8

    - name: quad4_apex_under_biaxial_tension            # DP-2bis (2026-05-19)
      setup: "Quad4 unitario plane strain, ε_xx = ε_yy ~ 1e-2 (predictor hidrostático extremo)"
      expect: "estado final cae en ápice: σ_xx ≈ σ_yy = k(α)/(3·η_f), σ_xy ≈ 0, en cada Gauss point"
      tol_rel: 1.0e-6

references:
  - "Drucker D.C., Prager W. (1952). Soil mechanics and plastic analysis or limit design. Quart. Appl. Math. 10, 157-165."
  - "Drucker D.C. (1953). Coulomb friction, plasticity, and limit loads. J. Appl. Mech. 21, 71-74. (calibración plane_strain_matched)"
  - "Simó J.C., Hughes T.J.R. (1998). Computational Inelasticity. Springer. §6 (cone plasticity, apex return)."
  - "de Souza Neto E.A., Perić D., Owen D.R.J. (2008). Computational Methods for Plasticity. Wiley. §8 (Drucker-Prager, tangentes algorítmicas detalladas)."
  - "Chen W.F., Han D.J. (1988). Plasticity for Structural Engineers. Springer. §7 (cohesivo-friccionales)."
```

---

## Implementación

- **Archivo**: [solidum/materials/drucker_prager_2d.py](../../solidum/materials/drucker_prager_2d.py).
- **Clase**: `DruckerPrager2D`, registrada vía `@MaterialRegistry.register`. `IS_SYMMETRIC = False` declarativo (la tangente puede ser asimétrica; el despachador agrega el flag sobre todos los materiales del dominio).
- **Kernel Numba** `_compute_drucker_prager_plane_strain`:
  1. Predictor elástico desviador-volumétrico (paralelo a J2 plane strain, $\varepsilon^p_{zz}$ libre).
  2. Check $f^\text{trial}$.
  3. Si elástico: devolver σ_trial restaurando contribución hidrostática y volumétrica.
  4. Si no: intentar $\Delta\gamma_\text{reg}$ cerrado; comprobar $\sqrt{J_2^{n+1}} \geq 0$.
  5. Si regular válido: actualizar σ, $\boldsymbol\varepsilon^p$, $\alpha$, tangente regular (de Souza Neto §8.3.4).
  6. Si regular inválido: rama de ápice — $\Delta\gamma_\text{apex}$ cerrado, σ puramente hidrostático, tangente volumétrica reducida.
- **Tangente regular** (§12): `C_alg = K·v⊗v + 2G·(1−β)·I_dev + 4G·β·n̂⊗n̂ − (1/A)·b_g⊗b_f` con `b = 2G·n̂ + 3K·η·v`.
- **Piso de la tangente en el ápice** (§12): `K_apex_eff < K·1e-6` ⇒ `K·1e-6`; actúa en todo retorno al ápice con `H = 0`.
- **Cálculo de $(\eta_f, k_0, \eta_g)$**: en `__init__` desde $(c_0, \phi, \psi, \text{variant})$ — sin coste runtime. Se valida $\psi \leq \phi$ (físicamente sensato).
- **`admissibility_scale`**: $k(\alpha) = k_0 + H_k \cdot \alpha$, $H_k = \xi^2 H$.
- **Tests**:
  - [tests/test_materials_unit.py](../../tests/test_materials_unit.py) · `TestDruckerPrager2D`.
  - [tests/test_solid_2d_drucker_prager.py](../../tests/test_solid_2d_drucker_prager.py).
  - [tests/test_consistent_tangent.py](../../tests/test_consistent_tangent.py) · `TestDruckerPrager2DConsistentTangent` — tangente regular y de ápice frente a diferencias finitas.

---

## Diálogo

- **2026-05-14** · Spec creada. Decisiones clave del alcance MVP:
  - Solo `plane_strain` en esta entrega; `plane_stress` con cono en σ_zz=0 es notoriamente delicado y de uso geotécnico raro — se difiere.
  - Parámetros físicos `(c_0, φ, ψ)` en lugar de adimensionales `(k_0, η_f, η_g)` — más natural para el usuario de geotecnia. Calibración interna por `variant`.
  - **No asociada por defecto** (ψ es parámetro independiente con default = φ que la fuerza a asociada solo si el usuario no especifica). Plasticidad asociada en suelos sobreestima dilatancia.
  - Endurecimiento solo en cohesión (lineal); fricción/dilatancia constantes.
  - `IS_SYMMETRIC = False` declarativo a nivel material (independiente del caso particular de uso; conservador y consistente con el patrón del proyecto: el despachador agrega).
- **2026-05-14** · Sin ADR. El modelo encaja en el patrón establecido por `VonMises2D` (kernel Numba, despacho hypothesis aunque solo soporte una, semántica trial/commit) y no rompe contratos transversales.
- **2026-09-22** · Auditoría global: retorno al ápice con ψ = 0: sin componente volumétrica del flujo no existe solución consistente con σ = C_e:(ε − ε^p) (con H = 0 se clampaba `A_apex` a 1e-30 y salía Δγ ≈ 1e30; con H > 0 se devolvía el σ del ápice sin la ε^p volumétrica que lo justifica). Ahora el kernel lanza `ValueError` al entrar en la rama del ápice con η_g = 0 y el constructor avisa al construir con ψ = 0. Aclarado que el endurecimiento actúa sobre `k(α) = k₀ + H·α` (radio del cono), no sobre la cohesión `c` (dk/dc vale 0.83–1.2 según la calibración). Nuevo `out_of_plane_stress(σ, state)`: σ_zz de plane strain desde ε^p_zz (helper `materials/_plane_strain.py`), validado contra `DruckerPrager3D` con ε_zz = 0.
- **2026-09-22** · Ensamblaje por lotes (ADR 0014): declara `STATE_SCHEMA` y un adaptador `@njit` con firma `MAT_SIG` que llama al mismo return mapping; la tolerancia de admisibilidad (ADR 0006) se evalúa dentro del kernel con las mismas operaciones que `admissibility_tol`. Formulación sin cambios; equivalencia con `compute_state` exacta en esfuerzo, tangente y estado (`tests/test_batch_kernels.py`).
- **2026-09-22** · Cierre del ADR 0014 (§9): `batch_out_of_plane_stress(sigma, S)`, versión por arreglos de `out_of_plane_stress` para el post-proceso por familia (exportador VTK); misma fórmula, `ε^p_zz` leída de la fila de estado cuando aplica. Verificado archivo VTK idéntico por ambos caminos.
- **2026-09-29** · Revisión documental (sin cambio de código ni de resultados): (1) $\alpha$ es el multiplicador plástico acumulado $\gamma$, no una deformación plástica equivalente, y $H$ endurece el radio del cono $k$ respecto de $\gamma$; conversión desde la parametrización $c(\bar\varepsilon^p)$ de de Souza Neto et al. (2008), $H_k = \xi^2 H_c$ con $\xi = k_0/c_0$ (§5). Sustituye a la cota «dk/dc vale 0.83–1.2» de la entrada anterior, que sólo vale para $\phi = 30°$. (2) El piso $K\cdot 10^{-6}$ de la tangente en el ápice actúa en todo retorno al ápice con $H_k = 0$, no sólo con $\psi = 0$, caso que desde el 2026-09-22 lanza `ValueError` (§12 y `numerical_caveats`).
- **2026-09-29** · Erratas de la spec frente al código, que es correcto (tangente regular y de ápice verificadas por diferencias finitas en `tests/test_consistent_tangent.py`; ápice contra `k/(3η_f)` en `tests/test_materials_unit.py`). Sólo cambia el texto; ni el código ni la formulación implementada cambian. (a) §10: el esfuerzo en el ápice es $\boldsymbol\sigma^{n+1} = (k_{n+1}/(3\eta_f))\,\mathbf I$, no $(k_{n+1}/\eta_f)\,\mathbf I$; es la presión $p = k/(3\eta_f)$ de la línea siguiente, del código (`p_new = k_new/(3η_f)`) y de la spec 3D. (b) §10: la primera expresión del $\Delta\gamma_\text{apex}$ encuadrado, con numerador $\eta_f(p^\text{trial} - k)$ y denominador $3K\eta_f\eta_g + H_k$, no era igual a la final; queda sólo la correcta, $(I_1^\text{trial}\eta_f - k)/(9K\eta_f\eta_g + H_k)$, con el paso intermedio del despeje. (c) §12: la forma compacta de la tangente regular usaba $\mathbf b = G\hat{\mathbf n} + 3K\eta\,\mathbf 1$ con $\hat{\mathbf n} = \mathbf s/(2\sqrt{J_2})$ y omitía el término $+4G\beta\,\hat{\mathbf n}\otimes\hat{\mathbf n}$; el código usa $2G\hat{\mathbf n}$ y lleva ese término, como documenta la spec `DruckerPrager3D` §14. Se retira además la forma «por linealización» previa, con términos $\Theta$, $\mathbf a_f$, $\mathbf a_g$ y «corrección» sin definir. Erratas menores encontradas al hacerlo, en el mismo sentido: la actualización del ápice escribía $p^{n+1} = K(\operatorname{tr}\boldsymbol\varepsilon_{n+1} - 3\eta_g\Delta\gamma)$ sin restar $\boldsymbol\varepsilon^p_n$ (ahora $p^\text{trial} - 3K\eta_g\Delta\gamma$ con $p^\text{trial}$ de §8); el desviador que absorbe $\boldsymbol\varepsilon^p$ en el ápice es $\mathbf e^\text{trial}_\text{dev}$ de §8; §9 definía $\hat{\mathbf n} = \mathbf s/\sqrt{2J_2}$ (normal unitaria) y §12 $\mathbf s/(2\sqrt{J_2})$, ahora ambas la segunda, como el código; las tangentes de ápice usan $\mathbf v = [1,1,0]$ en Voigt 3, y el `expect` de `regular_return_under_shear` decía $\sqrt{J_2} = \sqrt{2/3}\,k - \eta_f I_1$, cuando el criterio de §2 da $\sqrt{J_2} = k - \eta_f I_1$ (el test comprueba $\alpha > 0$ y el cortante, no ese valor).
- **2026-09-29** · Cambio del significado del parámetro `H`, aprobado por el usuario (hallazgo clas-6 de la auditoría): pasa a ser $H_c = dc/d\bar\varepsilon^p$, la parametrización de de Souza Neto et al. (2008, §8.3), y el material calcula $H_k = \xi^2 H_c$ con $\xi = k_0/c_0$. Antes `H` era directamente $H_k$: un $H_c$ tomado de un ensayo o de la literatura daba un endurecimiento distinto por el factor $\xi^2$, sin aviso. Con endurecimiento lineal la física es la misma; cambia sólo el número que se escribe. El return mapping, la tangente y los resultados con `H = 0` no cambian. Blindado en `tests/test_drucker_prager_hardening.py` ($k(\gamma) = \xi\,c(\xi\gamma)$ en las tres calibraciones y estado final sobre la superficie endurecida con $H_k$); los tests con `H > 0` pasan a usar $H_k$ en sus expresiones analíticas.
