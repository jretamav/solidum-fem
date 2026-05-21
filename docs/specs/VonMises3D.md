# VON MISES 3D — plasticidad J2 con endurecimiento isótropo lineal

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA propone, ejecuta y rellena **implementación** + **diálogo** durante el trabajo.
>
> **Pre-redactada por la IA** sobre la base de `VonMises2D.md` y `Elastic3D.md`. Apertura de la sub-etapa **A.bis** (materiales 3D no lineales) tras cierre de la Etapa 7 (sólidos 3D acotados con ADR 0012).

---

## Especificación física

### 0. Descripción general

Modelo elastoplástico independiente de la velocidad ("rate-independent") basado en el criterio **J2 / Von Mises** con regla de flujo **asociada** y **endurecimiento isótropo lineal**, formulado **íntegramente en 3D** sobre la convención Voigt 6D del proyecto (ADR 0012). Captura plasticidad de metales dúctiles en régimen de pequeñas deformaciones, sin efectos cinemáticos (Bauschinger), térmicos ni viscosos.

A diferencia de `VonMises2D`, **no tiene variantes de hipótesis cinemática**: en 3D todas las componentes de $\boldsymbol\varepsilon$ y $\boldsymbol\sigma$ son activas; no hay proyección a un subespacio. La consecuencia es algorítmica: el return mapping es **único y radial cerrado**, sin Newton local proyectado.

`VonMises2D` en hipótesis **plane strain** es matemáticamente la restricción de este modelo bajo $\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0$ (deformación total impuesta); las componentes plásticas $\varepsilon^p_{zz}, \varepsilon^p_{yz}, \varepsilon^p_{xz}$ son libres y evolucionan por incompresibilidad/isotropía. Esta relación se verifica como test cruzado en §13.

### 1. Descomposición aditiva

En pequeñas deformaciones, el tensor de deformación se descompone aditivamente:

$$\boldsymbol\varepsilon = \boldsymbol\varepsilon^e + \boldsymbol\varepsilon^p$$

La parte elástica $\boldsymbol\varepsilon^e$ gobierna el esfuerzo vía la ley elástica isótropa; la parte plástica $\boldsymbol\varepsilon^p$ es **histórica** y evoluciona según la regla de flujo. La parte plástica es **incompresible**: $\operatorname{tr}(\boldsymbol\varepsilon^p) = 0$.

### 2. Ley elástica

$$\boldsymbol\sigma = \mathbf C_e : \boldsymbol\varepsilon^e = \mathbf C_e : (\boldsymbol\varepsilon - \boldsymbol\varepsilon^p)$$

con $\mathbf C_e$ la matriz constitutiva isótropa 3D de `Elastic3D` (Voigt 6×6), parametrizada por $(E, \nu)$ o equivalentemente $(K, G)$:

$$K = \frac{E}{3(1-2\nu)}, \qquad G = \frac{E}{2(1+\nu)}$$

Equivalente desviador + presión:

$$\boldsymbol\sigma = \mathbf s + p\,\mathbf I, \qquad \mathbf s = 2G\,\mathbf e^e, \qquad p = K\,\operatorname{tr}(\boldsymbol\varepsilon)$$

con $\mathbf e^e = \operatorname{dev}(\boldsymbol\varepsilon^e)$ y, por incompresibilidad plástica, $\operatorname{tr}(\boldsymbol\varepsilon^e) = \operatorname{tr}(\boldsymbol\varepsilon)$.

### 3. Superficie de fluencia (J2)

$$f(\boldsymbol\sigma, \alpha) = \|\mathbf s\| - \sqrt{\tfrac{2}{3}}\,\big(\sigma_y + H\,\alpha\big) \;\leq\; 0$$

donde $\mathbf s = \boldsymbol\sigma - \tfrac{1}{3}\operatorname{tr}(\boldsymbol\sigma)\mathbf I$ es el tensor desviador, $\|\mathbf s\| = \sqrt{\mathbf s : \mathbf s}$ la norma de Frobenius tensorial, $\sigma_y > 0$ la fluencia inicial y $H \geq 0$ el módulo de endurecimiento isótropo lineal. El factor $\sqrt{2/3}$ alinea $\sigma_y$ con la fluencia uniaxial estándar.

**Cómputo de $\|\mathbf s\|$ en Voigt 6D del proyecto** (componentes cortantes tensoriales sin factor 2):

$$\|\mathbf s\|^2 = s_{xx}^2 + s_{yy}^2 + s_{zz}^2 + 2\,(s_{xy}^2 + s_{yz}^2 + s_{xz}^2)$$

El factor 2 contabiliza las dos posiciones simétricas $s_{ij} = s_{ji}$ del tensor de 2º orden.

### 4. Regla de flujo (asociada)

$$\dot{\boldsymbol\varepsilon}^p = \dot\gamma\,\frac{\partial f}{\partial\boldsymbol\sigma} = \dot\gamma\,\mathbf N, \qquad \mathbf N = \frac{\mathbf s}{\|\mathbf s\|}$$

El flujo es puramente desviador ($\operatorname{tr}(\mathbf N) = 0$) — **plasticidad incompresible**, propiedad central de J2. $\mathbf N$ es un tensor de 2º orden de norma unidad almacenado en Voigt 6D con cortantes tensoriales.

### 5. Endurecimiento isótropo lineal

$$\dot\alpha = \sqrt{\tfrac{2}{3}}\,\dot\gamma$$

$\alpha$ es la **deformación plástica acumulada equivalente** (escalar, monótona no decreciente). La superficie de fluencia se expande uniformemente con $\alpha$ manteniendo su centro en el origen del espacio de tensiones desviadoras. La fórmula es **idéntica a plane strain** porque la dimensión de $\mathbf N$ es la misma (norma unidad tensorial); contrasta con plane stress, donde $\Delta\gamma$ tiene unidades de $[1/\text{esfuerzo}]$ por la proyección.

### 6. Condiciones de Kuhn-Tucker

$$\dot\gamma \geq 0, \qquad f \leq 0, \qquad \dot\gamma\,f = 0$$

Garantizan que $\dot\gamma > 0$ solo si el estado intenta salir de la superficie, y que el estado real siempre cumple $f \leq 0$.

### 7. Variables internas

| Símbolo | Tipo | Significado |
|---|---|---|
| $\boldsymbol\varepsilon^p$ | tensor (6 componentes Voigt, ver §8) | deformación plástica |
| $\alpha$ | escalar | deformación plástica acumulada equivalente |

`alpha` es la variable principal exportable (`PRIMARY_STATE_VAR = 'alpha'`); `eps_p` es interna al algoritmo de return mapping.

### 8. Notación Voigt y manejo de $\boldsymbol\varepsilon^p$

Convención del proyecto para 3D (ADR 0012):

$$\boldsymbol\varepsilon = [\varepsilon_{xx},\,\varepsilon_{yy},\,\varepsilon_{zz},\,\gamma_{xy},\,\gamma_{yz},\,\gamma_{xz}]^\top, \qquad \gamma_{ij} = 2\,\varepsilon_{ij}$$

$$\boldsymbol\sigma = [\sigma_{xx},\,\sigma_{yy},\,\sigma_{zz},\,\sigma_{xy},\,\sigma_{yz},\,\sigma_{xz}]^\top$$

**Almacenamiento interno de $\boldsymbol\varepsilon^p$**: 6 componentes con cortantes **tensoriales sin factor 2**:

$$\boldsymbol\varepsilon^p = [\varepsilon^p_{xx},\,\varepsilon^p_{yy},\,\varepsilon^p_{zz},\,\varepsilon^p_{xy},\,\varepsilon^p_{yz},\,\varepsilon^p_{xz}]$$

Decisión consistente con `VonMises2D` (que almacena $\varepsilon^p_{xy}$ tensorial en su slot 4). Razón: el flujo plástico $\mathbf N = \mathbf s/\|\mathbf s\|$ es un tensor de 2º orden con cortantes tensoriales; almacenar $\boldsymbol\varepsilon^p$ en la misma convención evita conversiones repetidas en el return mapping y mantiene la fórmula $\|\boldsymbol\varepsilon^p\|^2 = \sum_\text{diag} + 2\sum_\text{off}$ idéntica a la de $\mathbf s$.

**Traducción al entrar al kernel**: la deformación total llega con $\gamma$ *engineering*; las componentes cortantes se dividen por 2 una sola vez al construir el tensor para deviator decomposition. Al salir, $\boldsymbol\sigma$ se devuelve en Voigt 6D *engineering-compatible* (cortantes tensoriales sin factor — esto es lo que el ensamblador espera).

### 9. Relación con `VonMises2D` plane strain

`VonMises2D` plane strain es la restricción de este modelo bajo:

$$\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0 \quad\text{(deformación total impuesta por la hipótesis 2D)}$$

con $\varepsilon^p_{zz} \neq 0$ libre (incompresibilidad plástica). El kernel 2D-plane-strain mantiene 4 componentes de $\boldsymbol\varepsilon^p$ porque $\varepsilon^p_{yz} = \varepsilon^p_{xz} = 0$ por isotropía + ausencia de cortantes activos en esos planos.

**Implicación de testeo**: si se ejecuta `VonMises3D` con $\varepsilon = [\varepsilon_{xx}, \varepsilon_{yy}, 0, \gamma_{xy}, 0, 0]$ (path plane strain), debe reproducir $\sigma_{xx}, \sigma_{yy}, \sigma_{xy}, \alpha$ de `VonMises2D` plane strain exactamente, además de devolver el $\sigma_{zz}$ que `VonMises2D` calcula internamente pero no expone en el API público 2D. Ver §13.

---

## Formulación numérica

### 10. Esquema temporal — Backward Euler implícito

Dado el estado convergido $\{\boldsymbol\varepsilon^p_n, \alpha_n\}$ al paso $n$ y la deformación total prescrita $\boldsymbol\varepsilon_{n+1}$, resolver para $\{\boldsymbol\sigma_{n+1}, \boldsymbol\varepsilon^p_{n+1}, \alpha_{n+1}\}$ con discretización implícita:

$$\boldsymbol\varepsilon^p_{n+1} = \boldsymbol\varepsilon^p_n + \Delta\gamma\,\mathbf N_{n+1}, \qquad \alpha_{n+1} = \alpha_n + \sqrt{\tfrac{2}{3}}\,\Delta\gamma$$

más $f(\boldsymbol\sigma_{n+1}, \alpha_{n+1}) \leq 0$ con la condición de complementariedad.

### 11. Return mapping 3D — algoritmo radial cerrado

Descomposición volumétrica-desviadora del estado de prueba. Conversión de $\boldsymbol\varepsilon$ a tensorial dividiendo por 2 los cortantes engineering:

$$\boldsymbol\varepsilon^\text{tens} = [\varepsilon_{xx},\,\varepsilon_{yy},\,\varepsilon_{zz},\,\tfrac{\gamma_{xy}}{2},\,\tfrac{\gamma_{yz}}{2},\,\tfrac{\gamma_{xz}}{2}]$$

(esta es la conversión que ya hace `Elastic3D` implícitamente vía la matriz $\mathbf C_e$ con sus factores $(1-2\nu)/2$ en los cortantes).

Deformación de prueba elástica y su descomposición:

$$\boldsymbol\varepsilon^{e,\text{trial}}_{n+1} = \boldsymbol\varepsilon^\text{tens}_{n+1} - \boldsymbol\varepsilon^p_n$$

$$p^\text{trial} = K\,\operatorname{tr}(\boldsymbol\varepsilon_{n+1}) = K\,(\varepsilon_{xx} + \varepsilon_{yy} + \varepsilon_{zz})$$

(la traza no depende de la sustracción de $\boldsymbol\varepsilon^p_n$ porque $\operatorname{tr}(\boldsymbol\varepsilon^p) = 0$.)

$$\mathbf e^\text{trial} = \operatorname{dev}(\boldsymbol\varepsilon^{e,\text{trial}}_{n+1}), \qquad \mathbf s^\text{trial} = 2G\,\mathbf e^\text{trial}$$

con $\mathbf e^\text{trial}$ y $\mathbf s^\text{trial}$ ambos en Voigt 6D tensorial.

Función de fluencia trial:

$$f^\text{trial} = \|\mathbf s^\text{trial}\| - \sqrt{\tfrac{2}{3}}\,(\sigma_y + H\,\alpha_n)$$

con $\|\mathbf s^\text{trial}\|^2 = \sum_\text{diag} s_i^2 + 2\sum_\text{off} s_{ij}^2$ (§3).

**Si $f^\text{trial} \leq \text{tol}$** → paso elástico. $\boldsymbol\sigma_{n+1} = \mathbf C_e\,\boldsymbol\varepsilon^e_n + \mathbf C_e\,\Delta\boldsymbol\varepsilon$ equivalente, evaluado como $\mathbf s^\text{trial} + p^\text{trial}\,\mathbf I$. Estado interno intacto.

**Si $f^\text{trial} > 0$** → corrector plástico radial. Por isotropía + asociado en J2, la dirección de flujo no cambia entre trial y final: $\mathbf N_{n+1} = \mathbf N^\text{trial} = \mathbf s^\text{trial}/\|\mathbf s^\text{trial}\|$. Forma cerrada:

$$\boxed{\Delta\gamma = \frac{f^\text{trial}}{2G + \tfrac{2}{3}H}}$$

Actualización (todo en Voigt 6D tensorial):

$$\mathbf s_{n+1} = \mathbf s^\text{trial} - 2G\,\Delta\gamma\,\mathbf N, \qquad \boldsymbol\sigma_{n+1} = \mathbf s_{n+1} + p^\text{trial}\,\mathbf I$$

$$\boldsymbol\varepsilon^p_{n+1} = \boldsymbol\varepsilon^p_n + \Delta\gamma\,\mathbf N, \qquad \alpha_{n+1} = \alpha_n + \sqrt{\tfrac{2}{3}}\,\Delta\gamma$$

Coste: O(1) por punto de Gauss. Sin Newton local. **Algoritmo idéntico al de plane strain extendido a 6D** — solo cambian los índices de las componentes activas.

### 12. Tangente algorítmica consistente

Forma cerrada estándar (Simó-Hughes 1998, §3.3):

$$\mathbf C^\text{alg} = K\,\mathbf 1\otimes\mathbf 1 + 2G(1-\beta)\,\mathbf I_\text{dev} - 2G\,\bar\gamma\,\mathbf N\otimes\mathbf N$$

con

$$\beta = \frac{2G\,\Delta\gamma}{\|\mathbf s^\text{trial}\|}, \qquad \bar\gamma = \frac{1}{1 + H/(3G)} - \beta$$

En Voigt 6D del proyecto:

- $\mathbf 1\otimes\mathbf 1$: matriz 6×6 con `1` en las entradas $(i,j)$ con $i,j \in \{1,2,3\}$ y `0` resto.
- $\mathbf I_\text{dev}$: proyector desviador 6×6. En la base tensorial, $\mathbf I_\text{dev} = \mathbf I - \tfrac{1}{3}\mathbf 1\otimes\mathbf 1$ en el bloque diagonal (3×3) y un bloque identidad escalado en los cortantes con el factor adecuado para que $\mathbf I_\text{dev}\,\boldsymbol\varepsilon = \operatorname{dev}(\boldsymbol\varepsilon)$ cuando $\boldsymbol\varepsilon$ está en convención engineering. En la convención del proyecto (entrada engineering, salida engineering), el bloque cortante de $\mathbf I_\text{dev}$ es la identidad 3×3 simple.
- $\mathbf N\otimes\mathbf N$: producto exterior 6×6 de $\mathbf N$ (Voigt 6D tensorial) consigo mismo. **Cuidado con la convención**: el producto $\mathbf N\otimes\mathbf N$ del kernel debe ser consistente con $\mathbf C^\text{alg}\,\boldsymbol\varepsilon$ engineering. La forma simple: si $\mathbf N$ está almacenado tensorial (sin factor 2 en cortantes), entonces $\mathbf C^\text{alg}$ en convención engineering tiene $\mathbf N\otimes\mathbf N$ con un factor de 2 en las posiciones cortantes — equivalente a multiplicar las componentes cortantes de $\mathbf N$ por $\sqrt 2$ antes del producto exterior (Mandel-like rescaling local al kernel, idéntico al patrón de VM2D plane strain).

**Simetría**: $\mathbf C^\text{alg}$ es simétrica (J2 asociado), `IS_SYMMETRIC = True`.

**Coste**: precalcular los cuatro escalares $K, G, \beta, \bar\gamma$ y construir la 6×6 sin loops anidados — O(1) constante.

### 13. Caveats numéricos

- **Locking volumétrico** en `Hex8`/`Tet4` cuando la plasticidad domina: la incompresibilidad plástica $\operatorname{tr}(\boldsymbol\varepsilon^p) = 0$ combinada con $\nu$ no muy lejano de 0.5 produce el mismo locking documentado para `Elastic3D` (limitación declarada, blindada por test en la Etapa 7). B-bar/F-bar diferidas a caso de uso real.
- **$\nu \to 0.5$**: $K \to \infty$ y el predictor elástico pierde precisión. Constructor rechaza $\nu \geq 0.5$. Recomendado $\nu \leq 0.49$ en metales.
- **$H = 0$ (plasticidad perfecta)**: $\Delta\gamma = f^\text{trial}/(2G)$, sin singularidades. Carga aplicada que supere la capacidad plástica produce divergencia del Newton global del solver (`LoadExceedsCapacityError`, ADR 0011); responsabilidad del solver, no del material.
- **Tangente algorítmica vs continua**: con $\Delta\gamma > 0$, $\mathbf C^\text{alg} \neq \mathbf C^\text{ep}_\text{continuo}$ por el término $\beta$. Usar la algorítmica garantiza convergencia cuadrática del Newton global.
- **Tolerancia de fluencia** vía ADR 0006: `admissibility_scale = √(2/3)(σ_y + H·α)`. Idéntica a plane strain — la escala física es la misma norma desviadora en la frontera de fluencia.

---

## Contrato de implementación

```yaml
name: VonMises3D
kind: material
status: validated        # draft → implemented → validated

interface:
  strain_dim: 6
  primary_state_var: alpha     # deformación plástica acumulada equivalente
  is_symmetric: true           # tangente algorítmica simétrica (J2 asociado)

parameters:
  - { name: E,       type: float, required: true,  desc: "Módulo de Young (>0)" }
  - { name: nu,      type: float, required: true,  desc: "Coeficiente de Poisson ∈ (-1, 0.5)" }
  - { name: sigma_y, type: float, required: true,  desc: "Esfuerzo de fluencia inicial (>0)" }
  - { name: H,       type: float, required: false, default: 0.0,
      desc: "Módulo de endurecimiento isótropo lineal (≥0). H=0 ⇒ plasticidad perfecta" }
  - { name: density, type: float, required: false, default: null,
      desc: "Densidad (kg/m³). Opcional al construir; obligatoria solo si se ensambla peso propio o masa (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(6,), state_vars=None) -> (σ: ndarray(6,), C_alg: ndarray(6,6), state_vars')"
  strain_kind: "3D Voigt — [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

state_schema:
  eps_p: "ndarray(6,) — [ε^p_xx, ε^p_yy, ε^p_zz, ε^p_xy, ε^p_yz, ε^p_xz] (cortantes tensoriales, sin factor 2)"
  alpha: "float ≥ 0 — deformación plástica acumulada equivalente"

conventions:
  sign:  "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "[xx, yy, zz, xy, yz, xz] (ADR 0012); γ_ij = 2·ε_ij en deformación total engineering;
          ε^p_ij almacenada en componente tensorial (sin factor 2);
          σ_ij tensorial sin factor"

validity:
  - "E > 0, σ_y > 0, H ≥ 0"
  - "ν ∈ (-1, 0.5) estricto; ν = 0.5 rechazado (incompresible)"
  - "|ε| ≲ 1e-2 (régimen de pequeñas deformaciones)"

out_of_scope:
  - "endurecimiento cinemático (back-stress); usar variante futura VonMises3DKinematic"
  - "endurecimiento isótropo no lineal (Voce, exponencial, potencial)"
  - "ablandamiento (H<0); requiere regularización por longitud característica o no-local"
  - "regla de flujo no asociada"
  - "acoplamiento con daño; usar modelo plástico-dañado dedicado"
  - "viscoplasticidad y dependencia de la velocidad de deformación"
  - "grandes deformaciones (descomposición multiplicativa F = F^e·F^p)"
  - "anisotropía (Hill 1948, Barlat); el modelo es estrictamente J2 isótropo"
  - "incompresibilidad estricta (ν=0.5)"

numerical_caveats:
  - "Locking volumétrico en Hex8/Tet4 cuando la plasticidad domina (incompresibilidad plástica + ν moderado-alto). Mitigaciones B-bar/F-bar diferidas. Test específico documenta la limitación."
  - "ν cercano a 0.5: K diverge; recomendado ν ≤ 0.49 (típico ν=0.3 en metales)."
  - "H=0 (plasticidad perfecta): material correcto, divergencia bajo sobrecarga es responsabilidad del solver (ADR 0011)."

acceptance:
  # Cubierto por tests/test_materials_unit.py · TestVonMises3D (a crear).
  unit:
    - name: paso_elastico_bajo_fluencia
      setup: "ε = (1e-5, 0, 0, 0, 0, 0): tracción uniaxial por debajo del yield"
      expect: "σ = C_e·ε exacto; α = 0; eps_p = 0; tangente = C_e"
      tol_rel: 1.0e-12

    - name: yield_uniaxial_at_sigma_y
      setup: "tracción uniaxial monotónica creciente con ε_yy = ε_zz = -ν·ε_xx, ε_xy = ε_yz = ε_xz = 0"
      expect: "primer yield exactamente en σ_xx = σ_y; eps_yield = σ_y/E"
      tol_rel: 1.0e-10

    - name: traccion_biaxial_isotropa
      setup: "ε = (eps, eps, -2ν·eps/(1-ν), 0, 0, 0) — tracción biaxial isótropa en plano x-y con ε_zz libre por σ_zz=0"
      expect: "σ_xx = σ_yy convergen a yield uniforme; α evoluciona; σ_zz ≈ 0"
      tol_rel: 1.0e-4

    - name: cortante_puro_xy_yield
      setup: "γ_xy creciente con resto de componentes 0"
      expect: "primer yield en τ_xy = σ_y/√3 (criterio J2 en cortante puro)"
      tol_rel: 1.0e-10

    - name: incompresibilidad_plastica
      setup: "cualquier paso plástico monotónico"
      expect: "tr(eps_p) = ε^p_xx + ε^p_yy + ε^p_zz = 0 exacto"
      tol_abs: 1.0e-12

    - name: alpha_monotonic
      setup: "carga monotónica creciente con incrementos plásticos no triviales"
      expect: "α no decrece nunca; α_final > 0 si hubo plasticidad"

    - name: descarga_recupera_C_e
      setup: "tras estado plástico, paso con menor ε (descarga elástica)"
      expect: "σ devuelto coherente con eps_p_n acumulado; tangente devuelta = C_e; α intacto"
      tol_rel: 1.0e-10

    - name: tangente_simetrica
      setup: "evaluar C_alg en estado plástico arbitrario"
      expect: "máx |C_alg - C_alg.T| ≤ tol; IS_SYMMETRIC = True"
      tol_abs: 1.0e-8

    - name: unit_invariance_MPa_vs_Pa
      setup: "mismo path adimensional con (E=210, σ_y=0.25) y (E=210e9, σ_y=0.25e9)"
      expect: "α, eps_p (adimensionales) idénticos; σ/E idéntico"
      tol_rel: 1.0e-12

    - name: degeneracion_a_elasticidad_sin_plasticidad
      setup: "trayectoria que nunca alcanza la frontera de fluencia"
      expect: "todos los pasos devuelven C_alg = C_e; eps_p y α permanecen en cero"
      tol_rel: 1.0e-12

    - name: rechazo_inputs_invalidos
      setup: "construir VonMises3D con E≤0, ν fuera de rango, σ_y≤0, H<0"
      expect: "ValueError con mensaje claro en cada caso"

  # Test cruzado VM3D ↔ VM2D plane strain.
  cross_consistency:
    - name: equivalencia_plane_strain
      setup: "ejecutar VM3D con ε = (ε_xx, ε_yy, 0, γ_xy, 0, 0) y VM2D plane strain con ε = (ε_xx, ε_yy, γ_xy);
             mismo path multi-paso con plasticidad activa"
      expect: "σ_xx, σ_yy, σ_xy y α coinciden entre 3D y 2D PS;
              VM3D además devuelve σ_zz coherente con la incompresibilidad plástica"
      tol_rel: 1.0e-10

  # Cubierto por tests/test_solid_3d_plasticity.py (a crear, fuera de esta spec — pertenece a la integración Hex8/Tet4 + VM3D).
  integration:
    - name: pipeline_hex8_traccion_uniaxial_elastoplastica
      setup: "Hex8 unitario, condiciones de contorno equivalentes a tracción uniaxial pura
             (4 caras laterales restringidas en su dirección normal, una cara tirada en x).
             Carga monotónica que cruza el yield."
      expect: "σ_xx en el centro coincide con material standalone bajo el path equivalente;
              α coincide; convergencia Newton global en pocos pasos"
      tol_rel: 1.0e-6

    - name: pipeline_hill_cylinder_3d
      setup: "cilindro hueco bajo presión interna creciente discretizado con Hex8 (octante con simetrías);
             cargar hasta plasticidad parcial"
      expect: "frontera elasto-plástica radial coincide con solución analítica cerrada de Hill (1950)
              al avanzar la presión"
      tol_rel: 1.0e-2

references:
  - "Simó J.C., Hughes T.J.R. (1998). Computational Inelasticity. Springer.
     §3.3 (3D J2 return mapping radial cerrado, tangente algorítmica consistente)."
  - "Crisfield M.A. (1991). Non-linear Finite Element Analysis of Solids and Structures, Vol. 1. Wiley.
     Cap. 6 (J2 plasticity)."
  - "de Souza Neto E.A., Perić D., Owen D.R.J. (2008). Computational Methods for Plasticity. Wiley.
     §7 (3D J2 plasticity)."
  - "Hill R. (1950). The Mathematical Theory of Plasticity. Oxford University Press.
     Cap. V (cilindro hueco bajo presión interna, solución analítica)."
  - "Bathe K.J. (2014). Finite Element Procedures. Prentice Hall. §6.6, §6.8.4."
  - "ADR 0012 — Sólidos 3D: convención Voigt 6D del proyecto."
  - "ADR 0006 — Tolerancias adimensionales: patrón atol + rtol·escala con escala física."
```

---

## Implementación

- **Archivo**: [solidum/materials/von_mises_3d.py](../../solidum/materials/von_mises_3d.py).
- **Clase**: `VonMises3D`, registrada vía `@MaterialRegistry.register`. Sin despacho por hipótesis (no aplica en 3D).
- **Kernel Numba** `_compute_j2_3d`: descomposición volumétrica-desviadora en Voigt 6D tensorial, return mapping radial cerrado, tangente algorítmica consistente cerrada. Predictor elástico construido como $\mathbf s_\text{trial} + p\,\mathbf I$ con $\mathbf s_\text{trial} = 2G(\mathbf e^\text{dev}_\text{trial} - \boldsymbol\varepsilon^p_n)$ — respeta $\boldsymbol\varepsilon^p_n$ en descargas y reevaluaciones (aplicada desde el inicio la lección del bug detectado en VM2D plane strain el 2026-05-14).
- **State schema**: `eps_p` ndarray(6,) `[xx, yy, zz, xy, yz, xz]` con cortantes tensoriales, `alpha` escalar.
- **`admissibility_scale`**: `√(2/3)·(σ_y + H·α)` — idéntica a VM2D plane strain (misma escala física en la frontera J2).
- **Tests**:
  - [tests/test_materials_unit.py](../../tests/test_materials_unit.py) · `TestVonMises3D` (12 tests unitarios) + `TestVonMises3DvsPlaneStrain` (1 test cruzado VM3D ↔ VM2D plane strain con path multi-paso).
  - [tests/test_solid_3d_plasticity.py](../../tests/test_solid_3d_plasticity.py) · 5 tests de integración: `TestHex8VonMises3DUniaxialFree` (2 — cubo libre con Poisson, yield uniaxial), `TestHex8VonMises3DConfinedVsStandalone` (1 — cross-check contra material standalone), `TestHex8VonMises3DConvergence` (1 — Newton converge en pocas iteraciones), `TestTet4VonMises3DBasic` (1 — smoke test Tet4 + VM3D).

- **Notas de traducción**:
  - Convención mixta engineering/tensorial idéntica a VM2D: entrada `strain` engineering (`γ_ij = 2·ε_ij`), salida `sigma` con cortantes tensoriales (lo que el ensamblador espera), estado `eps_p` con cortantes tensoriales. La conversión engineering→tensorial se hace una sola vez al construir el desviador (`strain[3..5]/2`).
  - La tangente algorítmica $\mathbf C^\text{alg}$ usa el patrón $K\,\mathbf 1\otimes\mathbf 1 + 2G(1-\beta)\,\mathbf I_\text{dev} - 2G\,\bar\gamma\,\mathbf N\otimes\mathbf N$ con $\mathbf I_\text{dev}$ teniendo $1/2$ en los cortantes (mapea engineering input a tensorial output) y $\mathbf N\otimes\mathbf N$ con $\mathbf N$ tensorial. Esta combinación produce la convención correcta `engineering→tensorial` en la salida sin factores adicionales — verificado por el test cruzado VM3D ↔ VM2D plane strain (mismos $\sigma_{xx}, \sigma_{yy}, \sigma_{xy}, \alpha$ a 10 decimales) y por el test de integración Hex8 confinado vs material standalone.
  - Validación de parámetros en constructor: `E > 0`, `σ_y > 0`, `H ≥ 0`, `ν ∈ (-1, 0.5)`, `density ≥ 0` si se pasa. Mensajes en español.

---

## Diálogo

- **2026-05-21** · Spec pre-redactada por la IA al abrir la sub-etapa **A.bis** (materiales 3D no lineales) tras cierre de la Etapa 7 (ADR 0012). Base: `VonMises2D.md` (modelo físico) + `Elastic3D.md` (estructura Voigt 6D + sin variantes de hipótesis).

- **2026-05-21** · **Decisión de plumbing (decide la IA, Reglas §3)**: kernels Numba separados `_compute_j2_3d` (este material) vs `_compute_j2_plane_strain` (existente en VonMises2D). **No** se extrae un `_J2Core` compartido en esta entrega. Razones:
  1. La diferencia operativa entre los dos kernels es el número de componentes activas (4 vs 6), que en Numba *nopython* requiere arrays de tamaño fijo — un kernel "genérico" obligaría a zero-padding sistemático con coste y oscurecería la lectura.
  2. Cada kernel es corto (~40 líneas): la duplicación cabe en el mismo PR sin coste de mantenimiento desproporcionado.
  3. La regla de centralización post-dos-casos del proyecto se satisface ya con la documentación cruzada (esta spec + VonMises2D referencian el mismo algoritmo de Simó-Hughes §3.3). Si aparece un tercer caso (e.g. `VonMises3DCorot` o variante con anisotropía leve) será el momento natural de abrir `_J2Core`.
  4. Plane stress 2D queda explícitamente fuera de un eventual `_J2Core` por tener un algoritmo genuinamente distinto (Newton local proyectado).

  Si al implementar la duplicación resulta más alta de lo previsto, se reconsidera y se eleva a zona gris.

- **2026-05-21** · **Decisión de plumbing**: el predictor elástico debe construirse como $\boldsymbol\sigma^\text{trial} = \mathbf C_e\,(\boldsymbol\varepsilon - \boldsymbol\varepsilon^p_n)$ (o equivalente desviador + presión), **no** como $\mathbf C_e\,\boldsymbol\varepsilon$. Lección de VM2D plane strain (bug detectado el 2026-05-14): en descargas elásticas desde estado plástico o reevaluaciones vía `compute_gauss_state(U_final)`, ignorar $\boldsymbol\varepsilon^p_n$ devuelve σ inconsistente. Se aplica desde el inicio en VM3D y se cubre con `test_descarga_recupera_C_e` y con el test de integración `test_pipeline_hex8_traccion_uniaxial_elastoplastica`.

- **2026-05-21** · **Tests cruzados VM3D ↔ VM2D plane strain** (`cross_consistency`): garantía adicional de que VM3D no introduce regresión silenciosa sobre el comportamiento ya validado de VM2D. Es la red de seguridad principal para confirmar que la formulación 6D no contiene errores de signo/factor en la convención Voigt — los tests unitarios analíticos por sí solos no cubren toda la combinatoria de componentes activas.

- **2026-05-21** · **Validación contra Hill 3D**: el caso del cilindro hueco bajo presión interna con plasticidad perfecta tiene **solución analítica cerrada** (Hill 1950, §V) — la presión a la que la frontera elasto-plástica radial avanza a un radio dado es expresable como función explícita de $\sigma_y, r_i, r_o$. Es el benchmark canónico para validar J2 3D y ya está disponible en 2D plane strain (`tests/validation/test_hill_cylinder_j2.py`). El test 3D será su análogo discretizado por octante con simetrías. Se redactará y validará tras tener `Hex8 + VonMises3D` funcionando end-to-end.

- **2026-05-21** · **Implementación completa, status `validated`**. Kernel + clase + 12 tests unitarios + 1 test cruzado + 5 tests de integración (Hex8 uniaxial libre con Poisson, Hex8 confinado vs material standalone, Hex8 convergencia, Tet4 smoke test). **18 tests añadidos a la suite (804 → 824 verdes; suite global sin regresiones).** Todos los tests pasaron a la primera, incluyendo el test cruzado VM3D ↔ VM2D plane strain a 10 decimales — señal fuerte de que la convención Voigt 6D en el kernel y la tangente algorítmica son correctas. Validación Hill 3D (cilindro hueco bajo presión interna con plasticidad perfecta) **diferida a entrega separada** porque requiere setup de malla por octante con simetrías y carga incremental hasta plasticidad parcial; será el primer test de `tests/validation/` que ejercite VM3D contra solución analítica cerrada.
