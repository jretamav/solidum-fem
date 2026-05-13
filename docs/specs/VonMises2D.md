# VON MISES 2D — plasticidad J2 con endurecimiento isótropo lineal

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> Spec **retroactiva con ampliación**: la hipótesis `plane_strain` está implementada y testeada desde sesiones anteriores; esta spec la documenta y añade `plane_stress` como hipótesis a implementar en esta sesión.

---

## Especificación física

### 0. Descripción general

Modelo elastoplástico independiente de la velocidad ("rate-independent") basado en el criterio **J2 / Von Mises** con regla de flujo **asociada** y **endurecimiento isótropo lineal**. Captura plasticidad de metales dúctiles en régimen de pequeñas deformaciones, sin efectos cinemáticos (Bauschinger), térmicos ni viscosos.

Soporta dos hipótesis cinemáticas 2D mutuamente excluyentes seleccionadas en construcción:

- **`plane_strain`** — $\varepsilon_{zz} = \varepsilon_{xz} = \varepsilon_{yz} = 0$. Cuerpos prismáticos largos cargados transversalmente (presa, túnel, cilindro largo).
- **`plane_stress`** — $\sigma_{zz} = \sigma_{xz} = \sigma_{yz} = 0$. Láminas o membranas delgadas en su plano (chapa traccionada, placa fina).

### 1. Descomposición aditiva

En pequeñas deformaciones, el tensor de deformación se descompone aditivamente:

$$\boldsymbol\varepsilon = \boldsymbol\varepsilon^e + \boldsymbol\varepsilon^p$$

La parte elástica $\boldsymbol\varepsilon^e$ gobierna el esfuerzo vía la ley elástica isótropa; la parte plástica $\boldsymbol\varepsilon^p$ es **histórica** y evoluciona según la regla de flujo.

### 2. Ley elástica

$$\boldsymbol\sigma = \mathbf C_e : \boldsymbol\varepsilon^e = \mathbf C_e : (\boldsymbol\varepsilon - \boldsymbol\varepsilon^p)$$

con $\mathbf C_e$ tensor de rigidez elástica isótropa parametrizado por $(E, \nu)$ o equivalentemente $(K, G)$:

$$K = \frac{E}{3(1-2\nu)}, \qquad G = \frac{E}{2(1+\nu)}$$

### 3. Superficie de fluencia (J2)

$$f(\boldsymbol\sigma, \alpha) = \|\mathbf s\| - \sqrt{\tfrac{2}{3}}\,\big(\sigma_y + H\,\alpha\big) \;\leq\; 0$$

donde $\mathbf s = \boldsymbol\sigma - \tfrac{1}{3}\operatorname{tr}(\boldsymbol\sigma)\mathbf I$ es el tensor desviador, $\|\mathbf s\| = \sqrt{\mathbf s:\mathbf s}$, $\sigma_y > 0$ es la fluencia inicial y $H \geq 0$ el módulo de endurecimiento isótropo lineal. El factor $\sqrt{2/3}$ alinea $\sigma_y$ con la fluencia uniaxial estándar.

### 4. Regla de flujo (asociada)

$$\dot{\boldsymbol\varepsilon}^p = \dot\gamma \,\frac{\partial f}{\partial \boldsymbol\sigma} = \dot\gamma \,\mathbf N, \qquad \mathbf N = \frac{\mathbf s}{\|\mathbf s\|}$$

El flujo es puramente desviador ($\operatorname{tr}(\mathbf N) = 0$) — **plasticidad incompresible**, propiedad central de J2.

### 5. Endurecimiento isótropo lineal

$$\dot\alpha = \sqrt{\tfrac{2}{3}}\,\dot\gamma$$

$\alpha$ es la **deformación plástica acumulada equivalente** (escalar, monótona no decreciente). La superficie de fluencia se expande uniformemente con $\alpha$ manteniendo su centro en el origen del espacio de tensiones desviadoras.

### 6. Condiciones de Kuhn-Tucker

$$\dot\gamma \geq 0, \qquad f \leq 0, \qquad \dot\gamma\,f = 0$$

Garantizan que $\dot\gamma > 0$ solo si el estado intenta salir de la superficie, y que el estado real siempre cumple $f \leq 0$.

### 7. Variables internas

| Símbolo | Tipo | Significado |
|---|---|---|
| $\boldsymbol\varepsilon^p$ | tensor (4 componentes Voigt extendido en plane strain, ver §8) | deformación plástica |
| $\alpha$ | escalar | deformación plástica acumulada equivalente |

`alpha` es la variable principal exportable (`PRIMARY_STATE_VAR = 'alpha'`); `eps_p` es interna al algoritmo de return mapping.

### 8. Notación Voigt y manejo de $\varepsilon^p_{zz}$

Convención del proyecto para 2D: $\boldsymbol\varepsilon = [\varepsilon_{xx},\,\varepsilon_{yy},\,\gamma_{xy}]$ con $\gamma_{xy} = 2\varepsilon_{xy}$. Pero la **descomposición desviadora** requiere todas las componentes 3D, así que internamente el material maneja $\boldsymbol\varepsilon^p = [\varepsilon^p_{xx},\,\varepsilon^p_{yy},\,\varepsilon^p_{zz},\,\varepsilon^p_{xy}]$ (4 componentes, $\varepsilon^p_{xy}$ tensorial sin factor 2).

- **Plane strain**: $\varepsilon_{zz} = 0$ se impone *en deformación total*. La componente $\varepsilon^p_{zz}$ es no nula y evoluciona libremente (la incompresibilidad plástica acopla $\varepsilon^p_{zz} = -(\varepsilon^p_{xx} + \varepsilon^p_{yy})$). Por simetría, $\sigma_{zz}$ resulta no nulo pero no se expone en el API público 2D — vive solo en el cómputo interno.
- **Plane stress**: $\sigma_{zz} = 0$ se impone *en esfuerzo*. La componente $\varepsilon_{zz}$ es no nula (extensión transversal de Poisson + Poisson plástico), y aparece como variable adicional en el problema local de retorno (ver §10).

---

## Formulación numérica

### 9. Esquema temporal — Backward Euler implícito

Dado el estado convergido $\{\boldsymbol\varepsilon^p_n, \alpha_n\}$ al paso $n$ y la deformación total prescrita $\boldsymbol\varepsilon_{n+1}$, resolver para $\{\boldsymbol\sigma_{n+1}, \boldsymbol\varepsilon^p_{n+1}, \alpha_{n+1}\}$ con discretización implícita:

$$\boldsymbol\varepsilon^p_{n+1} = \boldsymbol\varepsilon^p_n + \Delta\gamma\,\mathbf N_{n+1}, \qquad \alpha_{n+1} = \alpha_n + \sqrt{\tfrac{2}{3}}\,\Delta\gamma$$

más $f(\boldsymbol\sigma_{n+1}, \alpha_{n+1}) \leq 0$ con la condición de complementariedad.

### 10. Return mapping plane strain — **algoritmo radial cerrado** (existente)

Descomposición volumétrica-desviadora del estado de prueba:

$$\boldsymbol\varepsilon^{\text{trial}}_{n+1} = \boldsymbol\varepsilon_{n+1} - \boldsymbol\varepsilon^p_n, \qquad \mathbf s^{\text{trial}} = 2G\,\mathbf e^{\text{trial}}, \qquad p^{\text{trial}} = K\,\operatorname{tr}(\boldsymbol\varepsilon_{n+1})$$

donde $\mathbf e^{\text{trial}}$ es la parte desviadora de $\boldsymbol\varepsilon^{\text{trial}}_{n+1}$. Función de fluencia trial:

$$f^{\text{trial}} = \|\mathbf s^{\text{trial}}\| - \sqrt{\tfrac{2}{3}}\,(\sigma_y + H\,\alpha_n)$$

Si $f^{\text{trial}} \leq \text{tol}$ → paso elástico, $\boldsymbol\sigma_{n+1} = \mathbf C_e:\boldsymbol\varepsilon_{n+1}$, estado intacto.

Si no, el flujo es radial (la dirección no cambia entre trial y final por isotropía + asociado en J2):

$$\Delta\gamma = \frac{f^{\text{trial}}}{2G + \tfrac{2}{3}H} \qquad\text{(forma cerrada)}$$

$$\mathbf N = \frac{\mathbf s^{\text{trial}}}{\|\mathbf s^{\text{trial}}\|}, \quad \mathbf s_{n+1} = \mathbf s^{\text{trial}} - 2G\,\Delta\gamma\,\mathbf N, \quad \boldsymbol\sigma_{n+1} = \mathbf s_{n+1} + p^{\text{trial}}\,\mathbf I$$

**Tangente consistente algorítmica** (derivada del algoritmo, no de la ley continua):

$$\mathbf C^{\text{alg}} = K\,\mathbf 1\otimes\mathbf 1 + 2G(1-\beta)\,\mathbf I_{\text{dev}} - 2G\,\bar\gamma\,\mathbf N\otimes\mathbf N$$

con $\beta = \dfrac{2G\,\Delta\gamma}{\|\mathbf s^{\text{trial}}\|}$, $\bar\gamma = \dfrac{1}{1 + H/(3G)} - \beta$.

### 11. Return mapping plane stress — **algoritmo proyectado** (a implementar)

Referencia: Simó & Hughes 1998, *Computational Inelasticity*, §3.4.1 ("Plane stress projected algorithm").

La restricción $\sigma_{zz} = 0$ se impone **eliminando** $\varepsilon_{zz}$ del problema mediante un operador de proyección, no como ecuación extra. Resultado: el problema local sigue siendo escalar en $\Delta\gamma$, sin nested Newton.

**Matriz elástica plane stress** (3×3, base Voigt $[\varepsilon_{xx},\varepsilon_{yy},\gamma_{xy}]$):

$$\mathbf C^{ps}_e = \frac{E}{1-\nu^2}\begin{bmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & (1-\nu)/2 \end{bmatrix}$$

**Predictor elástico** en el subespacio plane stress:

$$\boldsymbol\sigma^{\text{trial}} = \mathbf C^{ps}_e\,(\boldsymbol\varepsilon_{n+1} - \mathbf P_\varepsilon\,\boldsymbol\varepsilon^p_n)$$

con $\mathbf P_\varepsilon$ extrayendo las componentes plane $[\varepsilon^p_{xx},\varepsilon^p_{yy},2\varepsilon^p_{xy}]$ del tensor 4D.

**Operador $\mathbf P$ — norma desviadora bajo $\sigma_{zz}=0$** (Simó-Hughes Box 3.1):

$$\mathbf P = \frac{1}{3}\begin{bmatrix} 2 & -1 & 0 \\ -1 & 2 & 0 \\ 0 & 0 & 6 \end{bmatrix}$$

Esta matriz codifica el producto escalar desviador en el subespacio plane stress; reemplaza $\|\mathbf s\|^2$ por $\boldsymbol\sigma^\top \mathbf P\,\boldsymbol\sigma$ con todas las componentes 3D ya eliminadas analíticamente.

**Función de fluencia plane stress proyectada**:

$$\bar f(\boldsymbol\sigma, \alpha) = \tfrac{1}{2}\,\boldsymbol\sigma^\top \mathbf P\,\boldsymbol\sigma - \tfrac{1}{3}\big(\sigma_y + H\,\alpha\big)^2 \;\leq\; 0$$

**Esquema de retorno**:

$$\boldsymbol\sigma_{n+1} = \big[\mathbf C^{ps}_e\big]^{-1}\,\big(\boldsymbol\sigma^{\text{trial}}\big) - \Delta\gamma\,\mathbf P\,\boldsymbol\sigma_{n+1}$$

reordenando:

$$\boldsymbol\sigma_{n+1} = \mathbf A(\Delta\gamma)^{-1}\,\boldsymbol\sigma^{\text{trial}}, \qquad \mathbf A(\Delta\gamma) = \mathbf I + \Delta\gamma\,\mathbf C^{ps}_e\,\mathbf P$$

La matriz $\mathbf A$ es 3×3, diagonalizable en los autovectores ortonormales de $\mathbf C^{ps}_e\mathbf P$. Los autovalores son cerrados:

$$\mu_1 = \frac{E}{3(1-\nu)}, \quad \mu_2 = 2G, \quad \mu_3 = 2G$$

con autovectores $\mathbf v_1 = \tfrac{1}{\sqrt 2}[1,1,0]$ (hidrostático plano), $\mathbf v_2 = \tfrac{1}{\sqrt 2}[1,-1,0]$ (desviador plano), $\mathbf v_3 = [0,0,1]$ (cortante). Diagonalizan simultáneamente $\mathbf C^{ps}_e$ y $\mathbf P$, con $\mathbf v_i^\top\mathbf P\,\mathbf v_i \in \{1/3, 1, 2\}$.

**Actualización de la deformación plástica acumulada equivalente**. La regla $\dot\alpha = \sqrt{2/3}\,\dot\gamma$ usada en plane strain **no** es válida aquí: en plane stress projected $\dot\gamma$ tiene unidades $[1/\text{esfuerzo}]$ (porque $\dot{\boldsymbol\varepsilon}^p = \dot\gamma\,\mathbf P\boldsymbol\sigma$ con $\mathbf P\boldsymbol\sigma$ en unidades de esfuerzo). La fórmula físicamente correcta y dimensionalmente consistente, obtenida de $\dot\alpha = \sqrt{2/3}\,\|\dot{\boldsymbol\varepsilon}^p\|_\text{Frob}$:

$$\alpha_{n+1} = \alpha_n + \Delta\gamma \,\sqrt{\tfrac{2}{3}\,\boldsymbol\sigma_{n+1}^\top\mathbf P\,\boldsymbol\sigma_{n+1}}$$

**Newton local** sobre $\Delta\gamma$. En base autovectores $\sigma_i(\Delta\gamma) = a_i / (1+\Delta\gamma\,\mu_i)$, así $\boldsymbol\sigma^\top\mathbf P\,\boldsymbol\sigma$ se expresa cerrado:

$$\boldsymbol\sigma^\top\mathbf P\,\boldsymbol\sigma(\Delta\gamma) = \frac{a_1^2}{3(1+\Delta\gamma\mu_1)^2} + \frac{a_2^2}{(1+\Delta\gamma\mu_2)^2} + \frac{2\,a_3^2}{(1+\Delta\gamma\mu_3)^2}$$

donde $a_i$ son las proyecciones de $\boldsymbol\sigma^{\text{trial}}$ en los autovectores. La ecuación residual:

$$\bar f(\Delta\gamma) = \tfrac{1}{2}\,\boldsymbol\sigma^\top\mathbf P\,\boldsymbol\sigma(\Delta\gamma) - \tfrac{1}{3}\big(\sigma_y + H\,\alpha(\Delta\gamma)\big)^2 = 0$$

con $\alpha(\Delta\gamma) = \alpha_n + \Delta\gamma\sqrt{(2/3)\,\boldsymbol\sigma^\top\mathbf P\,\boldsymbol\sigma}$ es monótona decreciente; tangente cerrada con regla de la cadena. **Newton local converge** típicamente en 3-6 iteraciones desde $\Delta\gamma = 0$. Coste: O(1) por punto de Gauss, comparable a plane strain.

**Tangente consistente plane stress**. Sea $\mathbf D = \mathbf A^{-1}\mathbf C^{ps}_e$, $\mathbf q = \boldsymbol\sigma_{n+1}^\top\mathbf P\,\mathbf D\,\mathbf P\,\boldsymbol\sigma_{n+1}$, $w = \sqrt{(2/3)\boldsymbol\sigma_{n+1}^\top\mathbf P\,\boldsymbol\sigma_{n+1}}$, $m = (2/3)\,R\,H$, $n = 2\Delta\gamma/(3w)$. De la condición de consistencia ($\dot{\bar f}=0$) con la regla correcta de $\alpha$:

$$\mathbf C^{\text{alg}}_{ps} = \mathbf D - \frac{(\mathbf D\,\mathbf P\boldsymbol\sigma_{n+1})\,(\boldsymbol\sigma_{n+1}^\top\mathbf P\,\mathbf D)}{q + m\,w / (1 - m\,n)}$$

Si $H=0$ (plasticidad perfecta) el término $m$ se anula y el denominador se reduce a $q$. Forma cerrada — sin diferenciación numérica.

**Predictor elástico con $\boldsymbol\varepsilon^p \neq 0$**. Tanto en plane strain como en plane stress, el predictor debe construirse como $\boldsymbol\sigma^{\text{trial}} = \mathbf C_e\,(\boldsymbol\varepsilon - \boldsymbol\varepsilon^p_n)$ (o equivalente desviador + presión), **no** como $\mathbf C_e\,\boldsymbol\varepsilon$, que ignoraría la deformación plástica acumulada. Importante en descargas elásticas desde estado plástico y al reevaluar el material con `compute_gauss_state(U_final)` tras un análisis converged.

### 12. Tolerancia de fluencia

Vía ADR 0006 (`Material.admissibility_tol`), patrón atol + rtol·escala:

- **Plane strain**: `admissibility_scale = √(2/3)(σ_y + H·α)` (norma desviadora característica en la frontera).
- **Plane stress**: `admissibility_scale = (σ_y + H·α)²/3` (escala de $\bar f$ proyectada; los dos términos son del mismo orden).

---

## Contrato de implementación

```yaml
name: VonMises2D
kind: material
status: validated

interface:
  strain_dim: 3
  primary_state_var: alpha   # deformación plástica acumulada equivalente

parameters:
  - { name: E,          type: float, required: true,  desc: "Módulo de Young" }
  - { name: nu,         type: float, required: true,  desc: "Coeficiente de Poisson" }
  - { name: sigma_y,    type: float, required: true,  desc: "Tensión de fluencia inicial" }
  - { name: H,          type: float, required: false, default: 0.0, desc: "Módulo de endurecimiento isótropo lineal (≥0). H=0 ⇒ plasticidad perfecta" }
  - { name: hypothesis, type: str,   required: false, default: "plane_strain", desc: "plane_strain | plane_stress" }
  - { name: density,    type: float, required: false, default: null, desc: "Densidad (kg/m³). Opcional al construir; obligatoria solo si se ensambla peso propio o masa (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(3,), state_vars=None) -> (σ: ndarray(3,), C_alg: ndarray(3,3), state_vars')"
  strain_kind: "plano Voigt — [ε_xx, ε_yy, γ_xy] con γ_xy = 2·ε_xy"

state_schema:
  eps_p:   "ndarray(4,) — [ε^p_xx, ε^p_yy, ε^p_zz, ε^p_xy] (xy tensorial, sin factor 2)"
  alpha:   "float ≥ 0 — deformación plástica acumulada equivalente"

conventions:
  sign: "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "γ_xy = 2·ε_xy en deformación total; ε^p_xy almacenada en componente tensorial (sin factor 2)"
  hypothesis_dispatch: "selección por kwarg al construir; mutuamente excluyentes"

validity:
  - "E > 0, σ_y > 0, H ≥ 0"
  - "ν ∈ (-1, 0.5) en plane_strain; ν ∈ (-1, 0.5) en plane_stress (ν=0.5 no admitido por singularidad)"
  - "|ε| ≲ 1e-2 (régimen de pequeñas deformaciones)"
  - "hypothesis ∈ {'plane_strain', 'plane_stress'}"

out_of_scope:
  - "endurecimiento cinemático (back-stress); usar variante futura VonMises2DKinematic"
  - "endurecimiento isótropo no lineal (Voce, exponencial)"
  - "endurecimiento por ablandamiento (H<0); requiere regularización (longitud característica, no-local)"
  - "regla de flujo no asociada"
  - "acoplamiento con daño; usar modelo plástico-dañado dedicado"
  - "viscoplasticidad y dependencia de la velocidad de deformación"
  - "grandes deformaciones (descomposición multiplicativa F = F^e·F^p)"
  - "anisotropía (Hill, Barlat) — el modelo es estrictamente J2 isótropo"

numerical_caveats:
  - "Plane strain: incompresibilidad plástica genera locking volumétrico en elementos de bajo orden (Quad4, Tri3) cuando ν→0.5 o el campo plástico domina. Mitigaciones futuras: B-bar, F-bar, elementos mixtos u-p. No incluidas en esta spec."
  - "Plane stress: el Newton local sobre Δγ asume H ≥ 0; con H = 0 (perfectamente plástico) la ecuación residual sigue siendo monótona pero su pendiente al final del paso plástico tiende a 0 — convergencia más lenta. Mitigación: bisección de respaldo si Newton local diverge."
  - "Para ν cercano a 0.5 en plane strain, K → ∞ y el predictor elástico pierde precisión numérica. Recomendado ν ≤ 0.49 en metales (típico ν=0.3)."

acceptance:
  # Cubierto por tests/test_materials_unit.py · TestVonMises2D y TestVonMises2DPlaneStress
  unit_plane_strain:
    - name: paso_elastico_bajo_fluencia
      expect: "ε=[1e-5,0,0]: σ=C_e·ε exacto, alpha=0, tangente=C_e"
      tol_rel: 1.0e-10
    - name: fluencia_traccion_uniaxial_confinada
      expect: "fluencia activa con ε=[0.01,-0.01,0]: alpha>0 tras un solo paso"
    - name: monotonicidad_alpha
      expect: "α no decrece bajo carga creciente con eps_yy=-ε; α_final > 0"
    - name: no_further_yield_on_frontier
      expect: "repetir ε en frontera de fluencia no incrementa α (return mapping consistente)"
      tol_rel: 1.0e-10
    - name: hypothesis_validation
      expect: "constructor rechaza hypothesis distinto de 'plane_strain'/'plane_stress' con ValueError"

  unit_plane_stress:
    - name: paso_elastico_bajo_fluencia
      expect: "ε pequeño: σ=C_e^ps·ε, tangente=C_e^ps, alpha=0, eps_p=0"
      tol_rel: 1.0e-10
    - name: degeneracion_a_1d_traccion_pura_elastico
      expect: "ε=(eps,-ν·eps,0): σ_xx=E·eps, σ_yy≈0 (regimen elástico, no E/(1-ν²))"
      tol_rel: 1.0e-4
    - name: yield_uniaxial_at_sigma_y
      expect: "primer yield en σ_xx=σ_y bajo tracción uniaxial pura (no E/(1-ν²)·ε)"
      tol_rel: 1.0e-2
    - name: traccion_biaxial_isotropa
      expect: "fluencia biaxial isótropa exactamente en σ_xx=σ_yy=σ_y (||s||=σ·√(2/3) iguala √(2/3)σ_y); ε_yield=σ_y(1-ν)/E"
      tol_rel: 1.0e-2
    - name: incompresibilidad_plastica
      expect: "tr(ε^p) = ε^p_xx + ε^p_yy + ε^p_zz = 0 tras cualquier paso plástico"
      tol_abs: 1.0e-12
    - name: alpha_monotonic
      expect: "α no decrece bajo carga monotónica creciente"
    - name: descarga_recupera_C_e_ps
      expect: "tras estado plástico, paso siguiente con menor ε produce tangente=C_e^ps y α intacto"
      tol_rel: 1.0e-10
    - name: tangente_simetrica
      expect: "C_alg simétrica (J2 asociado ⇒ IS_SYMMETRIC=True)"
      tol_abs: 1.0e-8
    - name: unit_invariance_MPa_vs_Pa
      expect: "mismo path adimensional en (MPa,mm,N) y (Pa,m,N) ⇒ α y ε^p idénticos"
      tol_rel: 1.0e-12

  # Cubierto por tests/test_solid_2d_plasticity.py
  integration:
    - name: pipeline_plane_strain_elastico
      setup: "Quad4 unitario, confinamiento uy=0; tracción horizontal por debajo del yield"
      expect: "u_x = ε_xx_target exacto; α=0"
      tol_rel: 1.0e-8

    - name: pipeline_plane_strain_plastico_vs_material_directo
      setup: "mismo setup, carga 1.5×σ_xx_yield → régimen plástico; integrar 10 pasos"
      expect: "σ_xx FEM (gauss avg) coincide con material standalone en ε_xx_final bajo trayectoria monótona; α coincide"
      tol_rel: 1.0e-3

    - name: pipeline_plane_stress_elastico_uniaxial
      setup: "Quad4 unitario, lado izq apoyado en ux; nodo (0,0) y (1,0) en uy; nodo (1,1) libre en y. Tracción horizontal por debajo del yield"
      expect: "u_x = ε_xx_target; contracción transversal u_y(0,1) = -ν·ε_xx (Poisson); α=0"
      tol_rel: 1.0e-6

    - name: pipeline_plane_stress_plastico_vs_1d
      setup: "mismo setup, carga 1.2σ_y para entrar plástico; comparar contra Elastoplastic1D bajo mismo path ε_xx"
      expect: "σ_xx FEM = σ_1D, σ_yy FEM ≈ 0 (uniaxial puro emergente); α_FEM = α_1D"
      tol_rel: 1.0e-3

    - name: pipeline_converges_few_iter_plane_strain
      setup: "carga plástica significativa, 4 incrementos, max_iter=8"
      expect: "solver converge en todos los pasos sin agotar max_iter; alpha final > 0"

    - name: pipeline_converges_few_iter_plane_stress
      setup: "tracción uniaxial pura plane stress 1.5σ_y, 4 incrementos, max_iter=8"
      expect: "solver converge en todos los pasos; alpha final > 0"

references:
  - "Simó J.C., Hughes T.J.R. (1998). Computational Inelasticity. Springer. §3.3 (plane strain return mapping), §3.4 (plane stress projected algorithm), Box 3.1 (operador P)."
  - "Crisfield M.A. (1991). Non-linear Finite Element Analysis of Solids and Structures, Vol. 1. Wiley. Cap. 6 (J2 plasticity)."
  - "Chen W.F., Han D.J. (1988). Plasticity for Structural Engineers. Springer. Cap. 5 (cilindro a presión interna, sol. cerrada Hill)."
  - "de Souza Neto E.A., Perić D., Owen D.R.J. (2008). Computational Methods for Plasticity. Wiley. §9.4 (plane stress projection alternativa)."
```

---

## Implementación

- **Archivo**: [fenix/materials/von_mises_2d.py](../../fenix/materials/von_mises_2d.py).
- **Clase**: `VonMises2D`, registrada vía `@MaterialRegistry.register`. Despacho por `hypothesis` en construcción (sin coste runtime).
- **Kernels Numba**:
  - `_compute_j2_plane_strain`: descomposición volumétrica-desviadora 3D ($\varepsilon^p_{zz}$ libre), return mapping radial cerrado, tangente algorítmica consistente cerrada. Predictor elástico construido como $\mathbf s_\text{trial} + p\mathbf I$ con $\mathbf s_\text{trial} = 2G(\mathbf e^\text{dev} - \boldsymbol\varepsilon^p_n)$ — respeta $\boldsymbol\varepsilon^p_n$ en descargas/reevaluaciones.
  - `_compute_j2_plane_stress`: predictor $\boldsymbol\sigma^\text{trial} = \mathbf C^{ps}_e(\boldsymbol\varepsilon - \boldsymbol\varepsilon^p_n)$, función de fluencia $\bar f$ proyectada, Newton local sobre $\Delta\gamma$ con tangente cerrada y máximo `_PLANE_STRESS_MAX_LOCAL_ITER = 25` iteraciones. Construcción de $\mathbf A^{-1}$ y σ_new vía base de autovectores conocida ($\mathbf v_1, \mathbf v_2, \mathbf v_3$); incompresibilidad plástica cierra $\varepsilon^p_{zz} = -(\varepsilon^p_{xx} + \varepsilon^p_{yy})$.
- **State schema**: `eps_p` ndarray(4,) `[xx, yy, zz, xy_tensorial]`, `alpha` escalar. Las cuatro componentes de `eps_p` son obligatorias en ambas hipótesis para preservar invariantes (incompresibilidad, normas tensoriales).
- **`admissibility_scale`**:
  - plane_strain: `√(2/3)·(σ_y + H·α)`
  - plane_stress: `(σ_y + H·α)²/3`
- **Tests**:
  - [tests/test_materials_unit.py](../../tests/test_materials_unit.py) · `TestVonMises2D` (5 tests plane strain pre-existentes, no regresión) + `TestVonMises2DPlaneStress` (9 tests plane stress).
  - [tests/test_solid_2d_plasticity.py](../../tests/test_solid_2d_plasticity.py) · 6 tests de integración del pipeline Quad4 + VonMises2D + NonlinearSolver en ambas hipótesis.
  - [tests/test_density_self_weight.py](../../tests/test_density_self_weight.py) · `test_vonmises2d_density` — densidad ADR 0008.
- **Notas de traducción**:
  - El kernel plane stress trabaja con $\boldsymbol\sigma$ en Voigt mixto: las componentes $\boldsymbol\sigma = [\sigma_{xx},\sigma_{yy},\sigma_{xy}]$ son las "engineering" mientras que $\boldsymbol\varepsilon$ usa $\gamma_{xy} = 2\varepsilon_{xy}$. El operador $\mathbf P$ en esta convención tiene la forma con `6` en la entrada (3,3); $\mathbf P\boldsymbol\sigma$ devuelve la componente cortante como $2\sigma_{xy}$ (convención engineering), que se convierte a tensorial dividiendo por 2 antes de acumular en $\varepsilon^p_{xy}$.
  - La actualización $\alpha_{n+1} = \alpha_n + \Delta\gamma\sqrt{(2/3)\boldsymbol\sigma^\top\mathbf P\boldsymbol\sigma}$ se evalúa **con $\boldsymbol\sigma$ final** del paso plástico (tras converger Newton local). Durante las iteraciones se usa la misma fórmula con el $\boldsymbol\sigma$ corriente.
  - Validación de parámetros en constructor: `E>0`, `σ_y>0`, `H≥0`, `ν∈(-1,0.5)`, `hypothesis ∈ {plane_strain, plane_stress}`. Mensajes en español describiendo la violación.

### Bug fixes detectados al implementar

- **Predictor elástico plane strain ignoraba $\boldsymbol\varepsilon^p$**. Versión previa devolvía $\boldsymbol\sigma = \mathbf C_e\boldsymbol\varepsilon$ cuando $f^\text{trial} \leq 0$ — correcto solo si $\boldsymbol\varepsilon^p = 0$. En descargas elásticas desde estado plástico, o al reevaluar el material vía `compute_gauss_state(U_final)` tras un análisis converged, devolvía σ inconsistente con el estado interno. Corregido a $\boldsymbol\sigma = \mathbf s_\text{trial} + p\mathbf I$. Detectado por el test de integración `test_plastic_regime_matches_standalone_material`: σ FEM reportado coincidía con el valor de equilibrio (F/A) en Newton (donde sí entra a return mapping), pero `compute_gauss_state` post-converged daba el valor incorrecto del predictor elástico, divergiendo del material standalone.

---

## Diálogo

- **2026-05-14** · Spec creada retroactivamente sobre material existente (J2 plane strain) y extendida con J2 plane stress. Decisión de algoritmo plane stress: **proyección de Simó-Hughes §3.4.1** en lugar de nested Newton sobre $\varepsilon_{zz}$. Justificación: forma cerrada del problema local (Newton local escalar sobre $\Delta\gamma$ con tangente cerrada, 3-6 iteraciones, O(1) por Gauss); Numba-compatible (kernel plano); tangente consistente cerrada vía operador $\mathbf P$. Nested Newton implicaría ~5-10× el coste por anidación de return mapping completo dentro de cada iteración externa sobre $\sigma_{zz}=0$, justificable solo con endurecimiento no lineal complejo (Voce, plasticidad mixta) que no es el caso. Registrado aquí, no en ADR — afecta a un material, no a contratos arquitecturales transversales.
- **2026-05-14** · Corrección durante la implementación: la regla heurística $\alpha_{n+1} = \alpha_n + \sqrt{2/3}\,\Delta\gamma$ heredada de plane strain rompe la invariancia bajo cambio de unidades en plane stress, porque allí $\Delta\gamma$ tiene unidades $[1/\text{esfuerzo}]$ (la regla de flujo $\dot{\boldsymbol\varepsilon}^p = \dot\gamma\,\mathbf P\boldsymbol\sigma$ tiene $\mathbf P\boldsymbol\sigma$ con dimensión de esfuerzo). La fórmula físicamente correcta es $\alpha_{n+1} = \alpha_n + \Delta\gamma\,\sqrt{(2/3)\boldsymbol\sigma^\top\mathbf P\boldsymbol\sigma}$ — adimensional. Detectado por el test unitario `test_unit_invariance_MPa_vs_Pa` que daba un factor 1e6 entre los dos sistemas. Implica recalcular la tangente consistente con $\partial\alpha/\partial\Delta\gamma \neq \sqrt{2/3}$ (ver §11).
- **2026-05-14** · Corrección plane strain: el predictor elástico devolvía $\boldsymbol\sigma = \mathbf C_e\boldsymbol\varepsilon$, ignorando $\boldsymbol\varepsilon^p_n$. Bug latente preexistente sin tests que lo cubrieran (los unitarios solo probaban primer paso con $\boldsymbol\varepsilon^p = 0$). Detectado al validar el pipeline FEM: `compute_gauss_state(U_final)` reportaba σ del predictor elástico (incorrecto) en lugar del σ post return-mapping. Corregido a $\boldsymbol\sigma = \mathbf s_\text{trial} + p\mathbf I$. La rama plane stress ya lo hacía bien desde el principio (predictor incluía la sustracción de $\boldsymbol\varepsilon^p$).
- **2026-05-14** · Validación numérica de plane strain queda implícita en los tests unitarios existentes; los tests de integración `tests/test_solid_2d_plasticity.py` abren por primera vez el pipeline completo Quad4 + VonMises2D + NonlinearSolver en ambas hipótesis — hasta hoy no había cobertura de integración sólido 2D + material plástico.
- **2026-05-14** · No se introduce ADR. La spec no rompe contratos: extiende un material existente con una rama nueva del despacho por `hypothesis`. Sí se actualiza `docs/catalogo_materiales.md` al validar.
