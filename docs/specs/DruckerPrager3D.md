# DRUCKER-PRAGER 3D — plasticidad friccional con endurecimiento isótropo lineal

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA propone, ejecuta y rellena **implementación** + **diálogo**.
>
> **Pre-redactada por la IA** sobre la base de `DruckerPrager2D.md` (mismo modelo físico) y `VonMises3D.md` (mismo patrón Voigt 6D). Segunda entrega de la sub-etapa **A.bis** (materiales 3D no lineales).

---

## Especificación física

### 0. Descripción general

Modelo elastoplástico **friccional** independiente de la velocidad. Suelos, hormigón, granulares, materiales cuasi-frágiles cohesivo-friccionales en régimen 3D completo. Suaviza la pirámide hexagonal de Mohr-Coulomb a un **cono circular** en el espacio de tensiones principales, eliminando las aristas que complican el return mapping. Formulado íntegramente en Voigt 6D del proyecto (ADR 0012).

Diferencias clave con J2 3D (`VonMises3D`):

- **Dependencia de la presión**. La resistencia al corte aumenta con la presión de confinamiento (componente friccional).
- **Plasticidad no asociada**. La dilatancia volumétrica real del flujo plástico ($\psi$) es típicamente menor que la fricción del criterio ($\phi$). Forzar $\psi = \phi$ (asociada) sobreestima la dilatancia.
- **Retorno al ápice**. El criterio define un cono con vértice; estados de prueba en zona puramente hidrostática extrema retornan al vértice (zona singular en la dirección de flujo).

Comparado con `DruckerPrager2D`:

- **Sin variantes de hipótesis** (igual que VM3D): en 3D todas las componentes son activas.
- **El catálogo de calibraciones cambia**: la variante 2D `plane_strain_matched` no tiene análogo en 3D (Mohr-Coulomb en 3D es una pirámide hexagonal sin coincidencia circular global). Las calibraciones 3D disponibles son `outer_cone` (circunscribe MC en el meridiano de compresión, **default**) e `inner_cone` (inscribe MC en el meridiano de extensión).

### 1. Descomposición aditiva y elasticidad

Pequeñas deformaciones: $\boldsymbol\varepsilon = \boldsymbol\varepsilon^e + \boldsymbol\varepsilon^p$. Ley elástica isótropa con la matriz $\mathbf C_e$ 6×6 de `Elastic3D`:

$$\boldsymbol\sigma = \mathbf C_e\,(\boldsymbol\varepsilon - \boldsymbol\varepsilon^p)$$

equivalente desviador + presión:

$$\boldsymbol\sigma = \mathbf s + p\,\mathbf I, \qquad \mathbf s = 2G\,\mathbf e^e_\text{dev}, \qquad p = K\,\operatorname{tr}(\boldsymbol\varepsilon - \boldsymbol\varepsilon^p)$$

A diferencia de J2, $\operatorname{tr}(\boldsymbol\varepsilon^p) \neq 0$ en general (dilatancia plástica si $\eta_g > 0$), por lo que $p$ depende de la deformación plástica acumulada vía $\operatorname{tr}(\boldsymbol\varepsilon^p)$.

### 2. Criterio de fluencia (Drucker-Prager)

$$f(\boldsymbol\sigma, \alpha) = \sqrt{J_2} + \eta_f\,I_1 - k(\alpha) \;\leq\; 0$$

con:

- $J_2 = \tfrac{1}{2}\,\mathbf s : \mathbf s$, segundo invariante del desviador.
- $I_1 = \operatorname{tr}(\boldsymbol\sigma) = \sigma_{xx} + \sigma_{yy} + \sigma_{zz}$, primer invariante.
- $\eta_f$, coeficiente friccional del cono (adimensional, $\geq 0$).
- $k(\alpha) = k_0 + H_k\,\alpha$, parámetro cohesivo con endurecimiento isótropo lineal.

La superficie es un **cono circular** con eje en la línea hidrostática ($\sigma_{xx} = \sigma_{yy} = \sigma_{zz}$), vértice en $I_1 = k/\eta_f$, abertura controlada por $\eta_f$.

**Cómputo de $J_2$ en Voigt 6D del proyecto** (componentes cortantes tensoriales sin factor 2):

$$J_2 = \tfrac{1}{2}\,\Big(s_{xx}^2 + s_{yy}^2 + s_{zz}^2 + 2\,(s_{xy}^2 + s_{yz}^2 + s_{xz}^2)\Big)$$

El factor 2 contabiliza las dos posiciones simétricas del tensor.

### 3. Calibración con Mohr-Coulomb (3D)

Los parámetros $(\eta_f, k_0)$ se derivan de los parámetros físicos $(c_0, \phi)$ — cohesión y ángulo de fricción interna — según la variante del cono elegida. En 3D, **el cono no puede coincidir con la pirámide hexagonal de MC en todo el plano π** (las aristas de MC quedan fuera o dentro del círculo), por lo que las calibraciones eligen un meridiano de referencia:

| variant | $\eta_f$ | $k_0$ | Comportamiento |
|---|---|---|---|
| `outer_cone` (default) | $2\sin(\phi)/[\sqrt 3\,(3 - \sin(\phi))]$ | $6\,c_0\cos(\phi)/[\sqrt 3\,(3 - \sin(\phi))]$ | Circunscribe MC en el meridiano de compresión triaxial ($\theta_L = -30°$). Más conservador en compresión, menos conservador en extensión. |
| `inner_cone` | $2\sin(\phi)/[\sqrt 3\,(3 + \sin(\phi))]$ | $6\,c_0\cos(\phi)/[\sqrt 3\,(3 + \sin(\phi))]$ | Inscribe MC en el meridiano de extensión triaxial ($\theta_L = +30°$). Cota inferior segura para diseño con MC. |

El default `outer_cone` se elige porque la compresión triaxial es el modo dominante en problemas geotécnicos típicos (cimentaciones, muros de contención, taludes). Si la aplicación prioriza la extensión (e.g. excavaciones), `inner_cone` es la elección segura.

**Caveat fundamental**: el cono no captura efectos de ángulo de Lode (la pirámide hexagonal de MC sí). Para problemas donde el modo de tensión rota fuera de los meridianos puros, considerar Mohr-Coulomb 3D directo (out-of-scope) o variantes con dependencia explícita del Lode angle (e.g. Matsuoka-Nakai, también out-of-scope).

**No se incluye `plane_strain_matched`** (variante 2D de DP2D): es un artefacto matemático que coincide con MC exactamente solo cuando $\varepsilon_{zz} = 0$ está impuesto cinemáticamente; no tiene análogo 3D porque MC en 3D depende del Lode angle activo.

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

con $\alpha$ la **deformación plástica equivalente** (convención del proyecto para Drucker-Prager). La superficie crece manteniendo eje y abertura, solo $k$ aumenta linealmente. Endurecimiento por fricción/dilatancia (cambio de $\eta_f$ o $\eta_g$ con $\alpha$) queda fuera de esta spec.

### 6. Condiciones de Kuhn-Tucker

$$\dot\gamma \geq 0,\qquad f \leq 0,\qquad \dot\gamma\,f = 0$$

### 7. Variables internas

| Símbolo | Tipo | Significado |
|---|---|---|
| $\boldsymbol\varepsilon^p$ | tensor 6 componentes `[xx, yy, zz, xy, yz, xz]` (cortantes tensoriales) | deformación plástica |
| $\alpha$ | escalar $\geq 0$ | deformación plástica equivalente (acumulado de $\Delta\gamma$) |

`alpha` es la variable principal exportable (`PRIMARY_STATE_VAR = 'alpha'`). $\boldsymbol\varepsilon^p$ **no es incompresible** en este modelo (a diferencia de J2): $\operatorname{tr}(\boldsymbol\varepsilon^p) = 3\eta_g\,\alpha$ en cualquier estado, lo que es un invariante cinemático verificable por test.

### 8. Notación Voigt 6D

Convención del proyecto (ADR 0012):

$$\boldsymbol\varepsilon = [\varepsilon_{xx},\,\varepsilon_{yy},\,\varepsilon_{zz},\,\gamma_{xy},\,\gamma_{yz},\,\gamma_{xz}]^\top, \qquad \gamma_{ij} = 2\,\varepsilon_{ij}$$

**Almacenamiento interno de $\boldsymbol\varepsilon^p$**: 6 componentes con cortantes **tensoriales sin factor 2**, idéntica convención a `VonMises3D` (consistencia entre materiales 3D no lineales del proyecto).

### 9. Relación con `DruckerPrager2D` plane strain

`DruckerPrager2D` plane strain con variante `outer_cone` (o `inner_cone`) es la restricción de este modelo bajo $\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0$ con $\varepsilon^p_{zz} \neq 0$ libre (dilatancia volumétrica). Esta relación se verifica como test cruzado en §13.

**No se puede establecer la equivalencia con la variante `plane_strain_matched` de DP2D**: esa calibración es 2D-only por construcción. El test cruzado usa `outer_cone` (o `inner_cone`) que es la única variante compartida.

---

## Formulación numérica

### 10. Esquema implícito — Backward Euler

Dado $\{\boldsymbol\varepsilon^p_n, \alpha_n\}$ y $\boldsymbol\varepsilon_{n+1}$, predictor elástico. Conversión engineering→tensorial dividiendo los cortantes por 2 (única vez al construir el desviador):

$$\boldsymbol\varepsilon^\text{tens}_{n+1} = [\varepsilon_{xx},\,\varepsilon_{yy},\,\varepsilon_{zz},\,\tfrac{\gamma_{xy}}{2},\,\tfrac{\gamma_{yz}}{2},\,\tfrac{\gamma_{xz}}{2}]$$

Deformación elástica de prueba (incluye toda la $\boldsymbol\varepsilon^p_n$, incluida su parte volumétrica por dilatancia):

$$\boldsymbol\varepsilon^{e,\text{trial}}_{n+1} = \boldsymbol\varepsilon^\text{tens}_{n+1} - \boldsymbol\varepsilon^p_n$$

Descomposición volumétrica-desviadora:

$$\operatorname{tr}(\boldsymbol\varepsilon^{e,\text{trial}}_{n+1}) = (\varepsilon_{xx} + \varepsilon_{yy} + \varepsilon_{zz}) - (\varepsilon^p_{xx,n} + \varepsilon^p_{yy,n} + \varepsilon^p_{zz,n})$$

$$\mathbf e^\text{trial} = \operatorname{dev}(\boldsymbol\varepsilon^{e,\text{trial}}_{n+1}), \qquad \mathbf s^\text{trial} = 2G\,\mathbf e^\text{trial}, \qquad p^\text{trial} = K\,\operatorname{tr}(\boldsymbol\varepsilon^{e,\text{trial}}_{n+1})$$

con $I_1^\text{trial} = 3\,p^\text{trial}$ y $\sqrt{J_2^\text{trial}} = \|\mathbf s^\text{trial}\|/\sqrt{2}$ (norma Frobenius en Voigt 6D tensorial).

Función de fluencia trial:

$$f^\text{trial} = \sqrt{J_2^\text{trial}} + \eta_f\,I_1^\text{trial} - k(\alpha_n)$$

Si $f^\text{trial} \leq \text{tol}$ → paso elástico. Si no → return mapping.

### 11. Return regular (cone surface)

La dirección desviadora $\hat{\mathbf n} = \mathbf s^\text{trial}/(2\sqrt{J_2^\text{trial}})$ es invariante bajo carga radial. Las actualizaciones son lineales en $\Delta\gamma$:

$$\sqrt{J_2^{n+1}} = \sqrt{J_2^\text{trial}} - G\,\Delta\gamma$$
$$p^{n+1} = p^\text{trial} - 3K\,\eta_g\,\Delta\gamma, \qquad I_1^{n+1} = 3\,p^{n+1}$$

La consistencia $f(\boldsymbol\sigma_{n+1}, \alpha_{n+1}) = 0$ con $\alpha_{n+1} = \alpha_n + \Delta\gamma$ y $k_{n+1} = k_0 + H_k(\alpha_n + \Delta\gamma)$ da:

$$\boxed{\;\Delta\gamma = \frac{f^\text{trial}}{G + 9K\,\eta_f\,\eta_g + H_k}\;}$$

**Forma cerrada** — sin Newton local. Idéntica al 2D plane strain en estructura algebraica.

**Actualizaciones**:

$$\mathbf s^{n+1} = \frac{\sqrt{J_2^{n+1}}}{\sqrt{J_2^\text{trial}}}\,\mathbf s^\text{trial}, \qquad \boldsymbol\sigma^{n+1} = \mathbf s^{n+1} + p^{n+1}\,\mathbf I$$

$$\boldsymbol\varepsilon^p_{n+1} = \boldsymbol\varepsilon^p_n + \Delta\gamma\,\left(\frac{\mathbf s^{n+1}}{2\sqrt{J_2^{n+1}}} + \eta_g\,\mathbf I\right)$$

donde $\mathbf I = [1,1,1,0,0,0]$ en Voigt 6D (parte volumétrica) y la componente desviadora de $\Delta\boldsymbol\varepsilon^p$ se almacena en convención tensorial.

**Condición de validez**: $\sqrt{J_2^{n+1}} \geq 0$. Si $\Delta\gamma_\text{reg} > \sqrt{J_2^\text{trial}}/G$, el desviador "se pasa de cero" y el algoritmo regular es inválido → toca **return al ápice**.

### 12. Return al ápice

Cuando el estado de prueba cae en zona puramente hidrostática extrema (predictor muy traccionante respecto a la cohesión), la proyección al cono cae en su vértice, donde $\mathbf s = 0$ y la dirección desviadora $\hat{\mathbf n}$ no está definida.

En el ápice: $\boldsymbol\sigma^{n+1} = (k_{n+1}/(3\eta_f))\,\mathbf I$ (puramente hidrostático), $\sqrt{J_2^{n+1}} = 0$.

$$\boxed{\;\Delta\gamma_\text{apex} = \frac{I_1^\text{trial}\,\eta_f - k(\alpha_n)}{9K\,\eta_f\,\eta_g + H_k}\;}$$

(Equivalente a la fórmula regular sin el término $G$, porque la rama de apex es puramente volumétrica.)

**Actualizaciones (apex)**:

$$\mathbf s^{n+1} = 0,\qquad p^{n+1} = \frac{k_{n+1}}{3\,\eta_f}, \qquad \boldsymbol\sigma^{n+1} = p^{n+1}\,\mathbf I$$

$$\boldsymbol\varepsilon^p_{n+1} = \boldsymbol\varepsilon^p_n + \mathbf e^\text{trial} + \Delta\gamma_\text{apex}\,\eta_g\,\mathbf I$$

(El desviador trial se absorbe completamente; la parte volumétrica del flujo añade $3\eta_g\Delta\gamma$ a la traza, distribuida $\eta_g\Delta\gamma$ por componente diagonal.)

### 13. Detección del régimen de retorno

Algoritmo de despacho (idéntico al 2D):

1. Calcular $\Delta\gamma_\text{reg}$ por la fórmula regular.
2. Si $\sqrt{J_2^\text{trial}} - G\,\Delta\gamma_\text{reg} \geq 0$ → return regular.
3. Si no → return al ápice con $\Delta\gamma_\text{apex}$.

### 14. Tangente algorítmica consistente

**Return regular** ($\sqrt{J_2^{n+1}} > 0$), forma compacta (de Souza Neto §8.3.4 extendida a 6D):

$$\mathbf C^\text{alg} = K\,\mathbf v\otimes\mathbf v + 2G\,\Big(1 - \tfrac{G\Delta\gamma}{\sqrt{J_2^\text{trial}}}\Big)\,\mathbf I_\text{dev} + 4G\,\beta\,\hat{\mathbf n}\otimes\hat{\mathbf n} - \frac{1}{A}\,\mathbf b_g \otimes \mathbf b_f$$

con:

- $\mathbf v = [1, 1, 1, 0, 0, 0]$ (gradiente de $I_1$ en Voigt 6D engineering input).
- $\mathbf I_\text{dev}$: proyector desviador 6×6 con $1/2$ en los cortantes (mapea engineering $\gamma_{ij}$ a tensorial $\varepsilon_{ij}$).
- $\hat{\mathbf n} = \mathbf s^\text{trial}/(2\sqrt{J_2^\text{trial}})$ en Voigt 6D tensorial.
- $\beta = G\,\Delta\gamma/\sqrt{J_2^\text{trial}}$.
- $\mathbf b_f = 2G\,\hat{\mathbf n} + 3K\,\eta_f\,\mathbf v$, $\mathbf b_g = 2G\,\hat{\mathbf n} + 3K\,\eta_g\,\mathbf v$.
- $A = G + 9K\,\eta_f\,\eta_g + H_k$.

**Asimétrica si $\eta_f \neq \eta_g$** (no asociada). El término $+4G\beta\,\hat{\mathbf n}\otimes\hat{\mathbf n}$ recoge el cambio de dirección de flujo con $\boldsymbol\varepsilon$ a través de $\mathbf s^\text{trial}$ (no solo el cambio de magnitud); sin él la tangente cierra mal en componentes de cortante alineados con el flujo (verificable por diferencia finita, idéntico al hallazgo del DP2D).

**Return al ápice** ($\sqrt{J_2^{n+1}} = 0$):

$$\mathbf C^\text{alg}_\text{apex} = \frac{K\,H_k}{9K\eta_f\eta_g + H_k}\,\mathbf v\otimes\mathbf v$$

Rigidez puramente volumétrica reducida; sin rigidez desviadora alguna. Con $H_k = 0$ y $\eta_g > 0$ la rigidez tangente colapsa por completo — el código devuelve un escalado mínimo de $K$ (idéntico patrón al 2D) para evitar matriz globalmente singular hasta que el siguiente paso salga del ápice.

### 15. Tolerancia de fluencia

ADR 0006: `admissibility_scale = k(α) = k_0 + H_k·α` (escala cohesiva corriente del cono). Tiene unidades de esfuerzo, igual que $f$.

---

## Contrato de implementación

```yaml
name: DruckerPrager3D
kind: material
status: validated        # draft → implemented → validated

interface:
  strain_dim: 6
  primary_state_var: alpha
  is_symmetric: false      # tangente asimétrica si no asociada (ψ ≠ φ);
                           # override por instancia a True cuando ψ = φ (mismo
                           # patrón que DP2D — el despachador algebraico ADR
                           # 0003 elige Cholesky/LDLᵀ en problemas asociados).

parameters:
  - { name: E,          type: float, required: true,  desc: "Módulo de Young (>0)" }
  - { name: nu,         type: float, required: true,  desc: "Coeficiente de Poisson ∈ (-1, 0.5)" }
  - { name: cohesion,   type: float, required: true,  desc: "Cohesión c_0 inicial (esfuerzo, >0)" }
  - { name: phi_deg,    type: float, required: true,  desc: "Ángulo de fricción interna en grados; típico suelos 20-40°, hormigón 30-37°" }
  - { name: psi_deg,    type: float, required: false, default: null,
      desc: "Ángulo de dilatancia en grados. Default = phi_deg (asociada). Típico no asociada: ψ < φ, e.g. 0 o φ/2." }
  - { name: H,          type: float, required: false, default: 0.0,
      desc: "Endurecimiento isótropo lineal en cohesión H_k (≥0). H=0 ⇒ perfectamente plástico." }
  - { name: variant,    type: str,   required: false, default: "outer_cone",
      desc: "outer_cone (default, circunscribe MC en compresión) | inner_cone (inscribe MC en extensión)" }
  - { name: density,    type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(6,), state_vars=None) -> (σ: ndarray(6,), C_alg: ndarray(6,6), state_vars')"
  strain_kind: "3D Voigt — [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

state_schema:
  eps_p: "ndarray(6,) — [ε^p_xx, ε^p_yy, ε^p_zz, ε^p_xy, ε^p_yz, ε^p_xz] (cortantes tensoriales, sin factor 2). tr(ε^p) = 3·η_g·α invariante cinemático."
  alpha: "float ≥ 0 — multiplicador plástico acumulado (= deformación plástica equivalente)"

conventions:
  sign:  "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "[xx, yy, zz, xy, yz, xz] (ADR 0012); γ_ij = 2·ε_ij engineering input;
          ε^p_ij almacenada tensorial; σ_ij tensorial sin factor"
  associated_flow: "ψ = φ ⇒ asociada (tangente simétrica); ψ ≠ φ ⇒ no asociada (tangente asimétrica)"

validity:
  - "E > 0, c_0 > 0, 0 ≤ φ < 90°, 0 ≤ ψ ≤ φ"
  - "ν ∈ (-1, 0.5) estricto"
  - "H ≥ 0"
  - "variant ∈ {'outer_cone', 'inner_cone'}"

out_of_scope:
  - "Mohr-Coulomb 3D con aristas (pirámide hexagonal): el cono circular DP no captura efectos de ángulo de Lode"
  - "Variantes con dependencia explícita del Lode angle (Matsuoka-Nakai, Lade-Duncan): formulaciones distintas"
  - "Calibración 'plane_strain_matched': sin análogo 3D (la coincidencia DP↔MC en plane strain es 2D-only)"
  - "Endurecimiento por cambio de φ y/o ψ con α (no lineal fuerte en return mapping)"
  - "Endurecimiento cinemático (back-stress hidrostático y desviador)"
  - "Cap (modelo cap-cone) para suelos compactables"
  - "Ablandamiento (H<0) — requiere regularización"
  - "Grandes deformaciones"

numerical_caveats:
  - "Tangente algorítmica asimétrica si ψ ≠ φ → IS_SYMMETRIC=False, el despachador algebraico usa LU (ADR 0003). Cuando ψ = φ se vuelve simétrica y el override por instancia activa Cholesky/LDLᵀ."
  - "Modos de carga muy traccionantes pueden caer en la rama de ápice. La detección automática del régimen (regular vs apex) está en el algoritmo."
  - "Con ψ = 0 (sin dilatancia) y H_k = 0, en la rama de ápice la tangente volumétrica se anula → singularidad. El código devuelve un valor de respaldo (rigidez volumétrica reducida no nula, K·1e-6) para evitar matriz singular; típicamente sale del ápice en el siguiente paso."
  - "Locking volumétrico en Hex8/Tet4 cuando ψ > 0 (flujo dilatante) combinado con ν moderado-alto. Mismo régimen documentado para Elastic3D/VonMises3D; mitigaciones B-bar/F-bar diferidas."
  - "El cono circular DP es independiente del ángulo de Lode. Para trayectorias de carga donde el modo rota fuera del meridiano de calibración (compresión triaxial → cortante → extensión triaxial), la diferencia DP vs MC puede ser hasta el 15-20% en la frontera. Es limitación inherente del modelo, no del algoritmo."

acceptance:
  unit:
    - name: parameter_validation
      setup: "constructor con variant inválida, φ ≥ 90°, ψ > φ, E≤0, c_0≤0, H<0, ν∉(-1,0.5)"
      expect: "ValueError con mensaje claro en cada caso"

    - name: calibration_outer_inner_consistency
      setup: "comparar η_f de outer_cone vs inner_cone para φ=30°"
      expect: "outer_cone η_f > inner_cone η_f (el outer es geométricamente mayor); ambos > 0; k_0 análogo"

    - name: paso_elastico_below_yield
      setup: "ε pequeño 6D, todas las componentes activas"
      expect: "σ = C_e·ε; eps_p y alpha intactos"
      tol_rel: 1.0e-10

    - name: regular_return_under_pure_shear
      setup: "ε cortante puro suficiente para fluir (γ_xy creciente, otras componentes 0)"
      expect: "return regular activo, α > 0; estado final cumple f ≤ tol"
      tol_rel: 1.0e-8

    - name: apex_return_under_hydrostatic_tension
      setup: "ε hidrostática traccionante (ε_xx = ε_yy = ε_zz > 0 suficientemente grande, sin cortante)"
      expect: "return al ápice; σ_xx = σ_yy = σ_zz = k(α)/(3η_f); √J₂ ≈ 0"
      tol_rel: 1.0e-8

    - name: tr_eps_p_invariant
      setup: "carga monotónica en rama regular (5 niveles), ψ ≠ φ (no asociada con η_g > 0)"
      expect: "tr(ε^p) = 3·η_g·α en cada paso plástico"
      tol_rel: 1.0e-10

    - name: associated_tangent_is_symmetric
      setup: "ψ = φ (asociada), estado plástico activo en rama regular"
      expect: "C_alg simétrica (‖C_alg − C_alg^T‖_∞ ≤ tol)"
      tol_abs: 1.0e-8

    - name: nonassociated_tangent_is_asymmetric
      setup: "ψ = 0, φ = 30°, estado plástico activo no degenerado"
      expect: "‖C_alg − C_alg^T‖_∞ con magnitud significativa relativa a ‖C_alg‖_∞"

    - name: tangent_matches_finite_difference
      setup: "estado plástico regular; comparar C_alg cerrada con diferencia finita centrada"
      expect: "error relativo < 1e-4 en todas las componentes"
      tol_rel: 1.0e-4

    - name: alpha_monotonic
      setup: "carga creciente con plasticidad sostenida"
      expect: "α no decrece nunca"

    - name: unloading_recovers_C_e
      setup: "cargar a plástico, descargar elásticamente"
      expect: "tangente = C_e en el paso de descarga; α se conserva"
      tol_rel: 1.0e-10

    - name: unit_invariance_MPa_vs_Pa
      setup: "mismo path adimensional con (MPa,mm,N) y (Pa,m,N)"
      expect: "α y eps_p (adimensionales) idénticos; σ/E idéntico"
      tol_rel: 1.0e-12

  cross_consistency:
    - name: equivalencia_plane_strain_outer_cone
      setup: "DP3D con variant='outer_cone' y ε = (ε_xx, ε_yy, 0, γ_xy, 0, 0);
             DP2D con variant='outer_cone' (plane_strain) y ε = (ε_xx, ε_yy, γ_xy);
             path multi-paso con plasticidad activa en rama regular"
      expect: "σ_xx, σ_yy, σ_xy y α coinciden entre 3D y 2D PS;
              componentes plásticas planas coinciden (xx, yy, zz, xy_tens);
              ε^p_yz = ε^p_xz = 0 en VM3D (no se activan)"
      tol_rel: 1.0e-10

  integration:
    - name: hex8_pipeline_elastico
      setup: "Hex8 unitario, carga muy por debajo del yield"
      expect: "respuesta elástica pura, α=0 en todos los Gauss points"
      tol_rel: 1.0e-8

    - name: hex8_pipeline_plastico_converges
      setup: "Hex8 unitario con confinamiento + cortante; carga plástica significativa; 4 pasos, max_iter=8"
      expect: "Newton converge en cada paso; α > 0 al final; tangente asimétrica activa LU"

    - name: hex8_apex_under_hydrostatic
      setup: "Hex8 unitario con ε hidrostática traccionante (todos los nodos del cubo desplazados radialmente desde el centro)"
      expect: "estado final cae en ápice en los 8 puntos de Gauss: σ_xx ≈ σ_yy ≈ σ_zz, √J₂ ≈ 0"
      tol_rel: 1.0e-6

    - name: tet4_basic_pipeline
      setup: "Tet4 con confinamiento + cortante, plasticidad activa"
      expect: "convergencia + α > 0 + tr(ε^p) = 3·η_g·α invariante en el único Gauss point"

references:
  - "Drucker D.C., Prager W. (1952). Soil mechanics and plastic analysis or limit design. Quart. Appl. Math. 10, 157-165."
  - "Simó J.C., Hughes T.J.R. (1998). Computational Inelasticity. Springer. §6 (cone plasticity, apex return)."
  - "de Souza Neto E.A., Perić D., Owen D.R.J. (2008). Computational Methods for Plasticity. Wiley. §8 (Drucker-Prager 3D, calibraciones outer/inner, tangentes algorítmicas detalladas)."
  - "Chen W.F., Han D.J. (1988). Plasticity for Structural Engineers. Springer. §7 (cohesivo-friccionales 3D)."
  - "ADR 0012 — Sólidos 3D: convención Voigt 6D del proyecto."
  - "ADR 0003 — Despachador algebraico: tangente asimétrica fuerza LU global."
  - "ADR 0006 — Tolerancias adimensionales: patrón atol + rtol·escala con escala física k(α)."
```

---

## Implementación

- **Archivo**: [solidum/materials/drucker_prager_3d.py](../../solidum/materials/drucker_prager_3d.py).
- **Clase**: `DruckerPrager3D`, registrada vía `@MaterialRegistry.register`. `IS_SYMMETRIC = False` a nivel de clase; override por instancia a `True` cuando `ψ = φ` (asociada).
- **Kernel Numba** `_compute_drucker_prager_3d`: predictor elástico desviador-volumétrico en Voigt 6D tensorial (sustracción de $\boldsymbol\varepsilon^p_n$ antes del desviador — *necesario* en DP por la dilatancia plástica), detección automática regular/apex, tangente cerrada para cada rama.
- **Calibración**: reusa la función `_calibrate_drucker_prager` de `drucker_prager_2d.py` por composición (mismas fórmulas para outer/inner; `plane_strain_matched` queda filtrada por el validador de variant antes de llegar a la función compartida).
- **State schema**: `eps_p` ndarray(6,) `[xx, yy, zz, xy, yz, xz]` con cortantes tensoriales, `alpha` escalar. **No es incompresible**: `tr(eps_p) = 3·η_g·α` por construcción del flujo.
- **`admissibility_scale`**: `k(α) = k_0 + H·α` (cohesión efectiva corriente).
- **Tests**:
  - [tests/test_materials_unit.py](../../tests/test_materials_unit.py) · `TestDruckerPrager3D` (17 tests unitarios) + `TestDruckerPrager3DvsPlaneStrain` (1 cruzado contra DP2D outer_cone plane_strain).
  - [tests/test_solid_3d_plasticity.py](../../tests/test_solid_3d_plasticity.py) · 4 tests de integración: `TestHex8DruckerPrager3DElastic` (1 — respuesta elástica), `TestHex8DruckerPrager3DRegularPlastic` (1 — cortante puro + tangente asimétrica activa LU + invariante $\operatorname{tr}(\boldsymbol\varepsilon^p) = 3\eta_g\alpha$), `TestHex8DruckerPrager3DApex` (1 — expansión hidrostática a 8 GP en rama ápice), `TestTet4DruckerPrager3DBasic` (1 — smoke test Tet4).

- **Notas de traducción**:
  - Convención mixta engineering/tensorial idéntica a VM3D: entrada `strain` engineering (`γ_ij = 2·ε_ij`), salida `sigma` con cortantes tensoriales, estado `eps_p` con cortantes tensoriales. Conversión engineering→tensorial una sola vez al construir el desviador (`strain[3..5]/2`).
  - El predictor construye el desviador como $\mathbf s^\text{trial} = 2G\,\operatorname{dev}(\boldsymbol\varepsilon^\text{tens} - \boldsymbol\varepsilon^p_n)$, **no** como $2G\,(\operatorname{dev}(\boldsymbol\varepsilon^\text{tens}) - \boldsymbol\varepsilon^p_n)$. A diferencia de J2, $\boldsymbol\varepsilon^p$ en DP tiene parte volumétrica no nula (dilatancia plástica si $\eta_g > 0$), así que la sustracción debe hacerse *antes* del operador desviador para que la parte volumétrica del flujo se contabilice correctamente.
  - Tangente regular usa $\mathbf I_\text{dev}$ 6×6 con $1/2$ en los cortantes (mapea engineering input a tensorial output), idéntica estructura al I_dev de VM3D.
  - Tangente apex: rigidez puramente volumétrica reducida con escalado mínimo $K \cdot 10^{-6}$ cuando $H_k = 0$ y $\eta_g$ pequeño (mismo patrón que DP2D para evitar matriz globalmente singular en el paso degenerado).
  - Validación de parámetros en constructor: `E > 0`, `c_0 > 0`, `0 ≤ φ < 90°`, `0 ≤ ψ ≤ φ`, `H ≥ 0`, `ν ∈ (-1, 0.5)`, `variant ∈ {outer_cone, inner_cone}`, `density ≥ 0`. Mensajes en español; rechazo explícito de `plane_strain_matched` con explicación de por qué es 2D-only.

---

## Diálogo

- **2026-05-21** · Spec pre-redactada por la IA al avanzar la sub-etapa **A.bis** (materiales 3D no lineales) tras cerrar `VonMises3D` con 824 tests verdes. Base: `DruckerPrager2D.md` (modelo físico, calibraciones, two-branch return mapping) + `VonMises3D.md` (estructura Voigt 6D, patrón cross-consistency).

- **2026-05-21** · **Decisión de plumbing (decide la IA, Reglas §3)**: kernels Numba separados `_compute_drucker_prager_3d` (este material) vs `_compute_drucker_prager_plane_strain` (existente en DP2D). **No** se extrae un `_DPCore` compartido en esta entrega, por las mismas razones que motivaron mantener separados los kernels de VM2D plane strain y VM3D: Numba *nopython* desfavorece arrays de tamaño variable, los kernels son cortos (~80 líneas cada uno), y el sub-patrón compartido (descomposición desviadora, branch detection, tangente cerrada) se documenta cruzadamente entre las specs.

- **2026-05-21** · **Decisión de alcance de variantes**: `DruckerPrager3D` ofrece **`outer_cone`** (default) e **`inner_cone`**, omitiendo `plane_strain_matched`. Razón: la variante 2D `plane_strain_matched` es un artefacto matemático que reproduce MC exactamente solo bajo $\varepsilon_{zz} = 0$ cinemático impuesto; no tiene análogo 3D porque MC en 3D depende del Lode angle activo. Default `outer_cone` se elige por ser el modo dominante en problemas geotécnicos (cimentaciones, muros). Si la aplicación es extension-dominated (excavaciones), el usuario puede elegir `inner_cone` explícitamente.

- **2026-05-21** · **Decisión de tipo de tests cruzados**: `cross_consistency` compara DP3D (`variant='outer_cone'`) contra DP2D (`variant='outer_cone'`, plane_strain), no contra `plane_strain_matched`. Razón: solo las dos calibraciones cono (outer/inner) existen en ambos modelos. La equivalencia matemática es exacta porque las fórmulas de calibración no dependen de la dimensión — solo de $(c_0, \phi, \psi)$ — y el cono 3D con $\varepsilon_{zz}=\gamma_{yz}=\gamma_{xz}=0$ impuesto se reduce a la formulación 2D plane strain en las componentes activas.

- **2026-05-21** · **Decisión de plumbing**: el predictor elástico construye el desviador como $\mathbf s^\text{trial} = 2G\,\operatorname{dev}(\boldsymbol\varepsilon^\text{tens} - \boldsymbol\varepsilon^p_n)$, **no** como $2G\,(\operatorname{dev}(\boldsymbol\varepsilon^\text{tens}) - \boldsymbol\varepsilon^p_n)$. Diferencia material: a diferencia de J2, $\boldsymbol\varepsilon^p$ en DP tiene parte volumétrica no nula ($\operatorname{tr}(\boldsymbol\varepsilon^p) = 3\eta_g\alpha$), así que la sustracción debe hacerse *antes* del desviador para que la parte volumétrica del flujo plástico se contabilice correctamente. Esta es la lección documentada en el comentario del kernel de DP2D (línea 39–42) y se aplica desde el inicio en DP3D.

- **2026-05-21** · **Implementación completa, status `validated`**. Kernel + clase + 17 tests unitarios + 1 test cruzado + 4 tests de integración (Hex8 elástico, Hex8 rama regular con tangente asimétrica + invariante $\operatorname{tr}(\boldsymbol\varepsilon^p) = 3\eta_g\alpha$, Hex8 rama ápice bajo expansión hidrostática, Tet4 smoke test). **24 tests añadidos a la suite (824 → 848 verdes; suite global sin regresiones).** Todos los tests pasaron a la primera, incluyendo el cruzado DP3D ↔ DP2D outer_cone plane_strain a 10 decimales — señal fuerte de que la convención Voigt 6D, la calibración compartida y la sustracción de $\boldsymbol\varepsilon^p_n$ antes del desviador son correctas. La integración Hex8 con cortante puro y tangente asimétrica ($\psi \neq \phi$) confirma que el despachador algebraico (ADR 0003) selecciona LU correctamente vía el flag declarativo `IS_SYMMETRIC = False`. Validación contra solución analítica cerrada (cilindro con falla DP, esfera hueca bajo presión interna, etc.) **diferida a la campaña 3D consolidada** que se hará al final de la sub-etapa A.bis (tras completar `IsotropicDamage3D`).
