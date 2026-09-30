# ISOTROPIC DAMAGE 3D — daño escalar isótropo con ablandamiento exponencial

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA propone, ejecuta y rellena **implementación** + **diálogo**.
>
> **Pre-redactada por la IA** sobre la base de `IsotropicDamage2D.md` (mismo modelo físico) y `Elastic3D.md` / `VonMises3D.md` (patrón Voigt 6D). Tercera y última entrega de la sub-etapa **A.bis** (materiales 3D no lineales). Desde el 2026-09-29 la **deformación equivalente es configurable** entre tres medidas de la literatura (ver §3 y *Diálogo*).

---

## Especificación física

### 0. Descripción general

Modelo de **daño isótropo escalar** en pequeñas deformaciones para sólidos 3D, formulado íntegramente en Voigt 6D del proyecto (ADR 0012). Captura degradación progresiva e isótropa de la rigidez en respuesta a la deformación máxima histórica, con ablandamiento exponencial. Pensado para hormigón, cerámicos y materiales cuasi frágiles donde la microfisuración distribuida puede modelarse mediante un escalar. Con las medidas de Mazars o de von Mises modificado (§3) el inicio del daño distingue tracción de compresión.

Diferencias respecto a `IsotropicDamage2D`:

- **Sin variantes de hipótesis** (igual que VM3D, DP3D, Elastic3D): en 3D todas las componentes son activas.
- **Algoritmo idéntico al 2D extendido a 6D**: misma ley de evolución exponencial, mismo régimen carga/descarga por $\kappa$, mismas tres medidas de deformación equivalente, misma estructura de tangente algorítmica consistente. **Sin Newton local** (el algoritmo es explícito).
- **`IsotropicDamage2D`** es la restricción exacta de este modelo en las dos hipótesis y con las tres medidas: en deformación plana bajo $\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0$; en tensión plana con $\varepsilon_{zz} = -\nu/(1-\nu)\,(\varepsilon_{xx} + \varepsilon_{yy})$ y $\gamma_{yz} = \gamma_{xz} = 0$ (y $\sigma_{zz} = 0$ resulta).

### 1. Cinemática

Pequeñas deformaciones, descomposición elástica pura: no hay deformación plástica (sin plasticidad). El daño se manifiesta como degradación de la rigidez efectiva. Sin variables internas tensoriales — solo el historial escalar $\kappa$ y el daño $d$.

### 2. Esfuerzo nominal vs esfuerzo efectivo

$$\boldsymbol\sigma = (1 - d)\,\boldsymbol\sigma^\text{eff}, \qquad \boldsymbol\sigma^\text{eff} = \mathbf C_e\,\boldsymbol\varepsilon$$

con $d \in [0, d_\text{max}]$ la variable de daño isótropa y $\mathbf C_e$ el tensor constitutivo elástico isótropo 3D 6×6 de `Elastic3D` (sin variantes de hipótesis en 3D).

El esfuerzo efectivo $\boldsymbol\sigma^\text{eff}$ es el que el material intacto soportaría con la deformación corriente. El esfuerzo nominal $\boldsymbol\sigma$ es lo que efectivamente transmite considerando la degradación.

### 3. Deformación equivalente

Cantidad escalar que cuantifica la "intensidad" de la deformación y gobierna el daño. Se elige con el parámetro `equivalent_strain` entre tres medidas estándar, en notación Voigt 6D $[\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, \gamma_{xy}, \gamma_{yz}, \gamma_{xz}]$ con $\gamma_{ij} = 2\varepsilon_{ij}$. Las tres están **calibradas en tracción uniaxial de esfuerzo**: $\varepsilon_\text{eq} = \varepsilon$, de modo que el daño empieza en $\sigma = E\,\kappa_0$ y $\kappa_0 = f_t/E$. Se diferencian en compresión y en estados multiaxiales.

1. **Energética** (`energy`, por defecto; Simó y Ju 1987):
   $$\varepsilon_\text{eq} = \sqrt{\frac{\boldsymbol\varepsilon^\top\mathbf C_e\,\boldsymbol\varepsilon}{E}} = \sqrt{\frac{\boldsymbol\sigma^\text{eff}:\boldsymbol\varepsilon}{E}}$$
   Simétrica en tracción y compresión (en compresión uniaxial de esfuerzo, $\varepsilon_\text{eq} = |\varepsilon|$).

2. **Mazars** (`mazars`; Mazars 1986):
   $$\varepsilon_\text{eq} = \sqrt{\textstyle\sum_{i=1}^{3}\langle\varepsilon_i\rangle_+^2}$$
   con $\varepsilon_i$ las deformaciones principales. Sólo las extensiones dañan: en compresión uniaxial de esfuerzo daña la extensión lateral, $\varepsilon_\text{eq} = \sqrt2\,\nu\,|\varepsilon|$; en compresión hidrostática no hay daño.

3. **von Mises modificado** (`modified_von_mises`; de Vree, Brekelmans y van Gils 1995), con $k = f_c/f_t$:
   $$\varepsilon_\text{eq} = \frac{k-1}{2k(1-2\nu)}\,I_1 + \frac{1}{2k}\sqrt{\left(\frac{k-1}{1-2\nu}\right)^2 I_1^2 + \frac{12k}{(1+\nu)^2}\,J_2}$$
   con $I_1 = \varepsilon_{xx} + \varepsilon_{yy} + \varepsilon_{zz}$ y $J_2 = \tfrac16\left[(\varepsilon_{xx}-\varepsilon_{yy})^2 + (\varepsilon_{yy}-\varepsilon_{zz})^2 + (\varepsilon_{zz}-\varepsilon_{xx})^2\right] + \tfrac14(\gamma_{xy}^2 + \gamma_{yz}^2 + \gamma_{xz}^2)$. En compresión uniaxial de esfuerzo, $\varepsilon_\text{eq} = |\varepsilon|/k$; con $k > 1$ no hay daño en compresión hidrostática.

Hasta el 2026-09-29 la medida era la norma del tensor, $\varepsilon_\text{eq} = \sqrt{\boldsymbol\varepsilon:\boldsymbol\varepsilon} = \sqrt{\boldsymbol\varepsilon^\top\mathbf M\boldsymbol\varepsilon}$ con $\mathbf M = \operatorname{diag}(1,1,1,\tfrac12,\tfrac12,\tfrac12)$: convención propia, no estándar y mal calibrada (en tracción uniaxial de esfuerzo daba $\varepsilon_\text{eq} = \varepsilon\sqrt{1+2\nu^2}$, así que $\kappa_0 \ne f_t/E$). Se retiró (ver *Diálogo*).

### 4. Variable histórica y condición de carga (Kuhn-Tucker)

$$\kappa_{n+1} = \max(\kappa_n, \varepsilon_\text{eq}(\boldsymbol\varepsilon_{n+1})), \qquad \kappa_0 \text{ inicial}$$

Esto codifica:

- **Carga activa** ($\dot\kappa > 0$): $\kappa$ crece con $\varepsilon_\text{eq}$; el daño puede aumentar.
- **Descarga / recarga** ($\dot\kappa = 0$): $\kappa$ queda congelado; el daño no cambia.
- **Irreversibilidad**: $\kappa$ es monótono no decreciente.

### 5. Ley de evolución del daño

$$d(\kappa) = \begin{cases}
0 & \kappa \leq \kappa_0 \\
1 - \dfrac{\kappa_0}{\kappa}\,e^{-\alpha(\kappa - \kappa_0)} & \kappa > \kappa_0
\end{cases}$$

con saturación numérica $d \leq d_\text{max} = $ `DAMAGE_MAX` (evita rigidez nula).

Propiedades idénticas al 2D: $d(\kappa_0) = 0$ (continuidad); $d'(\kappa) > 0$ para $\kappa > \kappa_0$ (irreversibilidad); $\lim_{\kappa \to \infty} d = 1$ (asintótico); $\alpha$ controla la velocidad de degradación. **Misma función de ablandamiento** que en 1D y 2D — centralizada en `solidum.materials._softening.evaluate_exponential_damage`.

---

## Formulación numérica

### 6. Algoritmo en un paso

Dado $\{\kappa_n\}$ y $\boldsymbol\varepsilon_{n+1}$:

1. $\boldsymbol\sigma^\text{eff} = \mathbf C_e\,\boldsymbol\varepsilon_{n+1}$ y $\varepsilon_\text{eq}$ con la medida elegida (§3), junto con su derivada. Las deformaciones principales de Mazars salen de `numpy.linalg.eigh` compilado con Numba.
2. Si $\varepsilon_\text{eq} > \kappa_n$: **carga activa**, $\kappa_{n+1} = \varepsilon_\text{eq}$. Si no: **descarga**, $\kappa_{n+1} = \kappa_n$.
3. Aplicar ley de daño con $\kappa_{n+1}$; saturar a $d_\text{max}$.
4. $\boldsymbol\sigma_{n+1} = (1 - d)\,\mathbf C_e\,\boldsymbol\varepsilon_{n+1}$.
5. Tangente: **secante** en descarga / sin daño / saturación; **algorítmica consistente** en carga activa con daño no saturado (ver §7).

El algoritmo es **explícito**: no requiere Newton local. La actualización de $\kappa$ y $d$ es directa, idéntica al 2D.

### 7. Tangente algorítmica consistente

$$\mathbf C^\text{alg} = \frac{\partial \boldsymbol\sigma_{n+1}}{\partial \boldsymbol\varepsilon_{n+1}} = (1-d)\,\mathbf C_e - \boldsymbol\sigma^\text{eff} \otimes \frac{\partial d}{\partial \boldsymbol\varepsilon},
\qquad
\frac{\partial d}{\partial \boldsymbol\varepsilon} = \frac{\partial d}{\partial \kappa}\,\frac{\partial \kappa}{\partial \boldsymbol\varepsilon}$$

con $\boldsymbol\sigma^\text{eff} = \mathbf C_e\,\boldsymbol\varepsilon$ (vector 6×1), $\partial d/\partial\kappa = (1 - d)(1/\kappa + \alpha)$ (idéntica al 2D) y $\partial\kappa/\partial\boldsymbol\varepsilon = \partial\varepsilon_\text{eq}/\partial\boldsymbol\varepsilon$ en carga activa, $\mathbf 0$ en descarga. Derivadas respecto de la deformación en Voigt *engineering*:

- **Energética**: $\dfrac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon} = \dfrac{\mathbf C_e\,\boldsymbol\varepsilon}{E\,\varepsilon_\text{eq}}$.
- **Mazars**: $\dfrac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon} = \dfrac{1}{\varepsilon_\text{eq}}\displaystyle\sum_i\langle\varepsilon_i\rangle_+\,\frac{\partial\varepsilon_i}{\partial\boldsymbol\varepsilon}$, con $\partial\varepsilon_i/\partial\boldsymbol\varepsilon$ el proyector $\mathbf n_i\otimes\mathbf n_i$ en Voigt *engineering* (componente de $\gamma_{ab}$: $n_a n_b$). La suma es la parte positiva del tensor, $\boldsymbol\varepsilon_+ = \sum_i\langle\varepsilon_i\rangle_+\,\mathbf n_i\otimes\mathbf n_i$, escrita en Voigt con sus componentes tensoriales: $[(\varepsilon_+)_{xx}, (\varepsilon_+)_{yy}, (\varepsilon_+)_{zz}, (\varepsilon_+)_{xy}, (\varepsilon_+)_{yz}, (\varepsilon_+)_{xz}]/\varepsilon_\text{eq}$. Bien definida con deformaciones principales repetidas: $\sum_i\langle\lambda_i\rangle_+^2$ es una función simétrica $C^1$ de los valores principales y los proyectores de un autoespacio repetido llevan el mismo peso.
- **von Mises modificado**: con $S = \sqrt{b^2 I_1^2 + c\,J_2}$, $b = (k-1)/(1-2\nu)$, $c = 12k/(1+\nu)^2$:
  $$\frac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon} = \frac{b + b^2 I_1/S}{2k}\,[1, 1, 1, 0, 0, 0] + \frac{c}{4kS}\,[e_{xx}, e_{yy}, e_{zz}, \tfrac12\gamma_{xy}, \tfrac12\gamma_{yz}, \tfrac12\gamma_{xz}]$$
  con $e$ la parte desviadora de la deformación.

**Tangente final**:

$$\mathbf C^\text{alg} = \begin{cases}
\mathbf C_e & \kappa \leq \kappa_0 \;(\text{sin daño}) \\
(1-d)\,\mathbf C_e & \text{descarga} \\
(1-d)\,\mathbf C_e \;-\; (1-d)\left(\dfrac1\kappa + \alpha\right)(\mathbf C_e\boldsymbol\varepsilon)\otimes\dfrac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon} & \text{carga activa, } d < d_\text{max} \\
(1-d_\text{max})\,\mathbf C_e & d = d_\text{max} \;(\text{saturación})
\end{cases}$$

### 8. Simetría de la tangente

- **Energética**: el término de carga activa es $\propto(\mathbf C_e\boldsymbol\varepsilon)\otimes(\mathbf C_e\boldsymbol\varepsilon)$: tangente consistente **simétrica** (en la implementación, bit a bit).
- **Mazars y von Mises modificado**: $\mathbf C_e\boldsymbol\varepsilon$ y $\partial\varepsilon_\text{eq}/\partial\boldsymbol\varepsilon$ no son proporcionales en general: tangente **no simétrica**.

Consecuencia (ADR 0003): `IS_SYMMETRIC` es `False` en la clase (cota conservadora) y la **instancia** lo pone a `True` con la medida energética; `domain_is_symmetric` lo lee de la instancia. Con Mazars o von Mises modificado el despachador elige LU (o Pardiso). Con la energética puede elegir Cholesky si está disponible; como la tangente de ablandamiento puede no ser definida positiva, el corrector degrada a LU en cuanto Cholesky la rechaza.

### 9. Decisión de régimen durante el Newton

Idéntica al 2D: la transición loading↔unloading se decide con $\kappa_n$ committed (no con $\kappa$ trial dentro del Newton). Estándar y evita oscilaciones del Newton entre ramas en pasos con ambigüedad. Si el paso converge, $\kappa$ se compromete y la decisión es consistente con la trayectoria.

---

## Contrato de implementación

```yaml
name: IsotropicDamage3D
kind: material
status: validated        # draft → implemented → validated

interface:
  strain_dim: 6
  primary_state_var: damage
  is_symmetric: false      # cota de clase; la instancia es true con equivalent_strain='energy' (ver §8)

admission:                # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - de_borst_crisfield: "§6.2 (daño isótropo escalar con ablandamiento exponencial y medidas de deformación equivalente)"
    - de_souza_neto: "cap. 12 (mecánica del daño)"
  general_codes:
    - code_aster: "MAZARS, ENDO_ISOT_BETON (daño escalar isótropo; MAZARS con la deformación equivalente de Mazars)"

parameters:
  - { name: E,       type: float, required: true,  desc: "Módulo de Young intacto (>0)" }
  - { name: nu,      type: float, required: true,  desc: "Coeficiente de Poisson ∈ (-1, 0.5)" }
  - { name: kappa_0, type: float, required: true,  desc: "Umbral de deformación equivalente; daño nulo mientras κ ≤ κ_0 (>0). Con las tres medidas κ_0 = f_t/E" }
  - { name: alpha,   type: float, required: true,  desc: "Velocidad de degradación exponencial (>0). α mayor → ablandamiento más abrupto" }
  - { name: equivalent_strain, type: str, required: false, default: "energy", desc: "energy (Simó-Ju 1987) | mazars (Mazars 1986) | modified_von_mises (de Vree et al. 1995)" }
  - { name: k,       type: float, required: false, default: null, desc: "f_c/f_t (≥ 1); obligatorio con modified_von_mises y no admitido con las otras medidas" }
  - { name: density, type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(6,), state_vars=None) -> (σ: ndarray(6,), C_tan: ndarray(6,6), state_vars')"
  strain_kind: "3D Voigt — [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

state_schema:
  kappa:  "float ≥ κ_0 — máxima deformación equivalente histórica"
  damage: "float ∈ [0, DAMAGE_MAX] — variable de daño escalar isótropa"

conventions:
  sign:  "σ > 0 ⇔ tracción; energy no distingue tracción de compresión; mazars sólo daña con extensiones principales; modified_von_mises daña en compresión uniaxial con ε_eq = |ε|/k"
  voigt: "γ_ij = 2·ε_ij en deformación; en ∂ε_eq/∂ε la componente de γ_ab es la derivada respecto de la deformación engineering"
  equivalent_strain: "energy: √(εᵀ·C_e·ε/E); mazars: √(Σ⟨ε_i⟩₊²); modified_von_mises: (k−1)/(2k(1−2ν))·I1 + 1/(2k)·√(((k−1)/(1−2ν))²·I1² + 12k/(1+ν)²·J2). Las tres dan ε_eq = ε en tracción uniaxial de esfuerzo"
  damage_irreversibility: "κ monótono no decreciente vía max(κ_old, ε_eq); d(κ) creciente"

validity:
  - "E > 0, κ_0 > 0, α > 0"
  - "ν ∈ (-1, 0.5) estricto"
  - "equivalent_strain ∈ {'energy', 'mazars', 'modified_von_mises'}"
  - "k ≥ 1 finito, sólo con equivalent_strain = 'modified_von_mises' (obligatorio ahí); con k < 1 no dañaría en tracción triaxial y sí en compresión hidrostática (decisión del usuario, 2026-09-30)"

out_of_scope:
  - "modelo de Mazars completo con dos daños (d_t, d_c) ponderados por α_t, α_c; aquí sólo su medida de deformación equivalente con un único daño escalar"
  - "anisotropía del daño (tensor de daño D vs escalar)"
  - "regularización para evitar mesh-dependency en regime de ablandamiento (longitud característica, no-local, gradient)"
  - "acoplamiento con plasticidad"
  - "fatiga / acumulación cíclica del daño"
  - "restitución de rigidez al cerrar fisuras (efecto unilateral)"

numerical_caveats:
  - "IS_SYMMETRIC por instancia: true con energy (tangente consistente simétrica), false con mazars y modified_von_mises (LU/Pardiso, ADR 0003). Con energy y Cholesky disponible, la tangente de ablandamiento puede no ser definida positiva: el corrector degrada a LU al primer rechazo."
  - "Transición loading↔unloading discontinua en la tangente. Newton converge antes de oscilar entre ramas en casos típicos. Para problemas patológicos cerca del umbral, considerar arc-length."
  - "Mazars: ε_eq es C¹ pero no C² donde una deformación principal cruza cero (entra o sale una extensión). ∂ε_eq/∂ε y la tangente consistente son continuas en ese cruce; lo que salta es su derivada (la curvatura de ⟨ε_i⟩₊²), y la tangente sigue siendo Lipschitz, así que el Newton conserva la convergencia cuadrática. Con deformaciones principales repetidas la derivada está bien definida. Las principales salen de numpy.linalg.eigh compilado (LAPACK) por punto de Gauss."
  - "Saturación d=DAMAGE_MAX corta la corrección consistente; el algoritmo devuelve tangente secante en ese caso para evitar términos espurios."
  - "Sin regularización, la solución es mesh-dependent en régimen de ablandamiento (localización en una banda de elementos). Es limitación inherente del modelo continuo de daño isótropo escalar; out-of-scope esta entrega."

acceptance:
  unit:
    - name: paso_elastico_bajo_umbral
      setup: "ε pequeño 6D con todas las componentes activas, ε_eq < κ_0"
      expect: "σ = C_e·ε exacto; d=0; C_tan = C_e (sin corrección rank-1)"
      tol_rel: 1.0e-12

    - name: damage_evolution_exponential
      setup: "ε grande tal que ε_eq >> κ_0"
      expect: "d cumple la ley exponencial; κ = ε_eq de la medida (energy, ε = [ε_xx, 0, …]: ε_eq = ε_xx·√((1−ν)/[(1+ν)(1−2ν)]))"
      tol_rel: 1.0e-10

    - name: eps_eq_matches_independent_reference
      setup: "25 estados aleatorios por medida; referencia sobre el tensor 3D (σ:ε con Lamé, autovalores de numpy, J2 = ½·e:e)"
      expect: "κ tras un paso desde κ = 0 coincide con ε_eq de referencia"
      tol_rel: 1.0e-12

    - name: uniaxial_stress_calibration
      setup: "tracción y compresión uniaxial de esfuerzo (ε = [ε, −νε, −νε, 0, 0, 0])"
      expect: "tracción: ε_eq = ε con las tres medidas (κ_0 = f_t/E); compresión: energy |ε|, mazars √2·ν·|ε|, modified_von_mises |ε|/k"
      tol_rel: 1.0e-13

    - name: hydrostatic_compression
      expect: "mazars y modified_von_mises (k > 1): ε_eq = 0; energy: ε_eq = e·√(3/(1−2ν))"

    - name: rotational_invariance
      expect: "ε_eq(R·ε·Rᵀ) = ε_eq(ε) con las tres medidas"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_on_unloading
      setup: "carga hasta κ > κ_0; descarga con ε de menor ε_eq"
      expect: "C_tan = (1-d)·C_e exacto (segundo término se anula con ∂κ/∂ε = 0)"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_below_threshold
      setup: "ε pequeño, κ_old = κ_0 (estado inicial intacto)"
      expect: "C_tan = C_e (sin daño activo)"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_at_saturation
      setup: "carga hasta κ enorme tal que d alcanza DAMAGE_MAX"
      expect: "C_tan = (1-DAMAGE_MAX)·C_e (el término consistente se corta)"
      tol_rel: 1.0e-12

    - name: tangent_symmetry_in_loading
      setup: "carga activa con estado de deformación no degenerado (ε_xx ≠ ε_yy ≠ ε_zz, cortantes activos)"
      expect: "energy: C_tan = C_tanᵀ bit a bit; mazars y modified_von_mises: ‖C_tan − C_tanᵀ‖_F/‖C_tan‖_F ≥ 0.01"

    - name: tangent_finite_difference_consistency
      setup: "carga activa desde κ = κ_0 con cada medida; estados general, todo extensión y dos principales repetidas; diferencia centrada con h = 1e-10 sobre |ε| ~ 1e-4"
      expect: "|C_FD − C_alg| ≤ 1e-6·max|C_FD| (truncamiento ~1e-12 y redondeo ~1e-10 relativos)"
      tol_rel: 1.0e-6

    - name: kappa_monotonic
      setup: "trayectoria multi-paso con carga, descarga, recarga"
      expect: "κ es monótono no decreciente en todos los pasos"

    - name: descarga_recarga_irreversibilidad
      setup: "cargar a κ_load, descargar elásticamente, recargar"
      expect: "d permanece constante durante toda la descarga; al recargar, d sigue creciendo desde el valor de descarga"

    - name: unit_invariance_MPa_vs_Pa
      setup: "mismo path adimensional con (MPa,mm,N) y (Pa,m,N)"
      expect: "d, κ y σ/E idénticos (κ_0 es una deformación, ε es adimensional → invariancia exacta)"
      tol_rel: 1.0e-12

    - name: rechazo_inputs_invalidos
      setup: "construir con E≤0, κ_0≤0, α≤0, ν∉(-1, 0.5), density<0, equivalent_strain desconocida, modified_von_mises sin k, k con otra medida, k ≤ 0 o no finito"
      expect: "ValueError con mensaje claro en cada caso"

    - name: is_symmetric_attribute
      expect: "IsotropicDamage3D.IS_SYMMETRIC == False en la clase; instancia True con energy, False con mazars y modified_von_mises"

    - name: batch_kernel_matches_compute_state
      expect: "el kernel por lotes reproduce bit a bit σ, C_tan y la fila de estado de compute_state con las tres medidas"

  cross_consistency:
    - name: equivalencia_plane_strain
      setup: "Damage3D con ε = (ε_xx, ε_yy, 0, γ_xy, 0, 0); Damage2D plane_strain con ε = (ε_xx, ε_yy, γ_xy);
             path multi-paso con carga, descarga, recarga"
      expect: "σ_xx, σ_yy, σ_xy, κ, d coinciden exactamente entre 3D y 2D PS;
              σ_zz en 3D no nulo (consistente con C_e plane_strain interno);
              componentes 3D ausentes en 2D PS (σ_yz, σ_xz) permanecen nulas"
      tol_rel: 1.0e-12

    - name: restriccion_2d_con_las_tres_medidas
      setup: "tests/test_damage_2d_3d_consistency.py: cuatro direcciones, tensión plana (ε_zz = −ν/(1−ν)·(ε_xx+ε_yy)) y deformación plana, energy/mazars/modified_von_mises"
      expect: "κ, d y σ en el plano iguales entre 2D y 3D; σ_zz 3D = 0 en tensión plana"
      tol_rel: 1.0e-12

  integration:
    - name: hex8_elastico
      setup: "Hex8 unitario con carga uniforme muy por debajo del umbral"
      expect: "d = 0 en todos los Gauss points; σ = C_e·ε exacto"
      tol_rel: 1.0e-8

    - name: hex8_damage_newton_converges_in_few_iter
      setup: "Hex8 unitario, carga uniaxial hasta régimen plenamente dañado (κ >> κ_0), 8 pasos"
      expect: "Newton converge en ≤6 iteraciones por paso (evidencia de convergencia cuadrática con tangente consistente)"

    - name: tet4_basic_pipeline
      setup: "Tet4 con confinamiento + tracción uniaxial, carga que activa daño"
      expect: "convergencia, d > 0 en el único Gauss point, κ ≥ κ_0"

    - name: uniaxial_confinado_curva_analitica
      setup: "tests/validation/test_uniaxial_softening_damage_3d.py: Hex8 confinado, u_x prescrito, las tres medidas"
      expect: "σ_xx y d coinciden con la curva analítica con el ε_eq de cada medida bajo confinamiento"
      tol_rel: 1.0e-10

    - name: yaml_modified_von_mises
      setup: "Hex8 en tracción uniaxial de esfuerzo leído desde YAML con equivalent_strain: modified_von_mises y k"
      expect: "reacción (1−d)·E·δ analítica (ε_eq = δ)"
      tol_rel: 1.0e-8

references:
  - "Lemaitre J., Chaboche J.-L. (1990). Mechanics of Solid Materials. Cambridge UP. Cap. 7 (continuum damage mechanics)."
  - "Simó J.C., Ju J.W. (1987). Strain- and stress-based continuum damage models — I. Formulation. Int. J. Solids Struct. 23, 821-840 (deformación equivalente energética; tangente consistente para daño isótropo)."
  - "Mazars J. (1986). A description of micro- and macroscale damage of concrete structures. Eng. Fract. Mech. 25, 729-737 (deformación equivalente de Mazars)."
  - "de Vree J.H.P., Brekelmans W.A.M., van Gils M.A.J. (1995). Comparison of nonlocal approaches in continuum damage mechanics. Comput. Struct. 55, 581-588 (deformación equivalente de von Mises modificado)."
  - "de Borst R., Crisfield M.A., Remmers J.J.C., Verhoosel C.V. (2012). Nonlinear Finite Element Analysis of Solids and Structures, 2.ª ed. Wiley. Cap. 6 (daño)."
  - "Oliver J. (1989). A consistent characteristic length for smeared cracking models. IJNME 28, 461-474."
  - "de Souza Neto E.A., Perić D., Owen D.R.J. (2008). Computational Methods for Plasticity. Wiley. §12 (damage models)."
  - "ADR 0012 — Sólidos 3D: convención Voigt 6D del proyecto."
  - "ADR 0003 — Despachador algebraico: simetría declarada por instancia; tangente asimétrica fuerza LU global."
  - "ADR 0006 — Tolerancias adimensionales: patrón atol + rtol·escala con escala física E·κ_0."
```

---

## Implementación

- **Archivo**: [solidum/materials/damage_3d.py](../../solidum/materials/damage_3d.py); medidas de deformación equivalente en [solidum/materials/_equivalent_strain.py](../../solidum/materials/_equivalent_strain.py), compartidas con `IsotropicDamage2D`.
- **Clase**: `IsotropicDamage3D`, registrada vía `@MaterialRegistry.register`. `IS_SYMMETRIC = False` como ClassVar; el constructor pone `self.IS_SYMMETRIC = True` con `equivalent_strain='energy'`.
- **Composición**:
  - `Elastic3D` como base elástica (acceso a $\mathbf C_e$ 6×6 vía `self.elastic_base.C`).
  - `evaluate_exponential_damage` de `solidum.materials._softening` para la ley de softening (centralizada, compartida con `IsotropicDamage1D` e `IsotropicDamage2D`).
  - `equivalent_strain_3d` de `solidum.materials._equivalent_strain`: devuelve $(\varepsilon_\text{eq}, \mathbf v, s)$ con $\partial\varepsilon_\text{eq}/\partial\boldsymbol\varepsilon = s\,\mathbf v$. Con la energética $\mathbf v = \boldsymbol\sigma^\text{eff}$ y $s = 1/(E\varepsilon_\text{eq})$, de modo que la tangente sale simétrica bit a bit. Mazars con `numpy.linalg.eigh` (atajo exacto cuando todas las principales son ≥ 0 o ≤ 0).
- **Estado**: dict `{'kappa': float, 'damage': float}`. `PRIMARY_STATE_VAR = 'damage'`.
- **Algoritmo `compute_state`**: la física vive en el núcleo `@njit` `_compute_damage_3d`, que consumen `compute_state` y el adaptador por lotes `_damage_3d_batch` (ADR 0014). Estructura idéntica a `IsotropicDamage2D` en Voigt 6D.
- **Camino por lotes**: `batch_params = [E, κ_0, α, tol_abs, tol_rel, DAMAGE_MAX, ν, código de la medida, k]` (códigos 0 energy, 1 mazars, 2 modified_von_mises; `k = 0` sin uso fuera de von Mises modificado). Cacheable en disco (verificado: kernel de familia y núcleos del material cargan del caché de Numba sin recompilar).
- **`admissibility_scale`**: `E·κ_0` constante respecto al estado (idéntico al patrón 1D/2D — el criterio se mide contra `κ_0`, no contra `κ` corriente).
- **Tests**:
  - [tests/test_materials_unit.py](../../tests/test_materials_unit.py) · `TestIsotropicDamage3D` + `TestIsotropicDamage3DvsPlaneStrain` (cruzado contra Damage2D plane_strain).
  - [tests/test_damage_equivalent_strain.py](../../tests/test_damage_equivalent_strain.py): las tres medidas (valores frente a referencia independiente, calibración, compresión hidrostática, invariancia, tangente por diferencias finitas, lotes, simetría, constructor, YAML).
  - [tests/test_damage_2d_3d_consistency.py](../../tests/test_damage_2d_3d_consistency.py): el 2D es la restricción exacta del 3D con las tres medidas.
  - [tests/test_solid_3d_plasticity.py](../../tests/test_solid_3d_plasticity.py) · 3 tests de integración: `TestHex8IsotropicDamage3DElastic` (1 — sin daño), `TestHex8IsotropicDamage3DActiveLoad` (1 — daño activo), `TestTet4IsotropicDamage3DBasic` (1 — smoke Tet4).
  - [tests/validation/test_uniaxial_softening_damage_3d.py](../../tests/validation/test_uniaxial_softening_damage_3d.py): curva σ-ε analítica bajo confinamiento con las tres medidas.

- **Notas de traducción**:
  - Convención engineering input (γ_ij): en $J_2$ los cortantes entran como $\gamma^2/4 = \varepsilon_{ij}^2$ y en el tensor de Mazars como $\gamma/2$; la derivada respecto de $\gamma_{ab}$ es la componente tensorial correspondiente del gradiente.
  - Validación de parámetros en constructor: `E > 0`, `κ_0 > 0`, `α > 0`, `ν ∈ (-1, 0.5)`, `density ≥ 0` si se pasa, `equivalent_strain` y `k`. Mensajes en español con el nombre del material y del parámetro.

---

## Diálogo

- **2026-05-21** · Spec pre-redactada por la IA al cerrar el segundo material de A.bis (`DruckerPrager3D`). Base: `IsotropicDamage2D.md` (modelo físico, tangente algorítmica consistente, semánticas carga/descarga, cap DAMAGE_MAX) + `VonMises3D.md`/`DruckerPrager3D.md` (estructura Voigt 6D, patrón cross-consistency).

- **2026-05-21** · **Decisión de plumbing (decide la IA, Reglas §3)**: reuso de la función centralizada `evaluate_exponential_damage` de `solidum.materials._softening`. Esta función ya está compartida entre `IsotropicDamage1D` e `IsotropicDamage2D` (auditoría H-3.4 cerrada el 2026-05-19); añadir `IsotropicDamage3D` como tercer consumidor satisface la regla de centralización del proyecto. Sin duplicación de la ley de softening.

- **2026-05-21** · **Decisión de plumbing**: el material usa `Elastic3D` internamente como base elástica (paralelo a cómo `IsotropicDamage2D` usa `Elastic2D`). Acceso a la matriz constitutiva 6×6 vía `self.elastic_base.C` — idéntico patrón al 2D, sin duplicación de la construcción de $\mathbf C_e$.

- **2026-05-21** · **Decisión de plumbing**: la matriz $\mathbf M$ del cálculo de $\varepsilon_\text{eq}$ y de la derivada $\partial\kappa/\partial\boldsymbol\varepsilon$ se materializa implícitamente en código (multiplicaciones por $1/2$ sobre los slots cortantes) en lugar de allocar una 6×6 — ahorro mínimo de memoria pero código más legible y consistente con el patrón seguido en `damage_2d.py`. *(Obsoleta desde el 2026-09-29: la norma del tensor se retiró.)*

- **2026-05-21** · **Cross-consistency** contra `IsotropicDamage2D` plane_strain, no plane_stress. La equivalencia matemática es exacta bajo $\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0$ impuesto (plane strain): la matriz $\mathbf C_e$ 3D restringida a ε_zz=0 input se reduce a la matriz $\mathbf C_e^\text{plane strain}$ 2D, y $\varepsilon_\text{eq}^\text{3D}$ se reduce a $\varepsilon_\text{eq}^\text{2D}$ bajo la misma restricción. Plane stress 2D usa un $\mathbf C_e$ proyectado distinto y por tanto no es comparable directamente con 3D bajo restricción cinemática. *(Desde el 2026-09-29 también la tensión plana es la restricción exacta del 3D, con el $\varepsilon_{zz}$ que anula $\sigma_{zz}$; ver `IsotropicDamage2D.md`.)*

- **2026-05-21** · **Sin variantes de hipótesis** en 3D (igual que VM3D, DP3D, Elastic3D). Constructor más simple que el 2D.

- **2026-05-21** · **No se incluye consistent tangent para `IsotropicDamage1D`** (deuda gemela del 2D): fuera del alcance de esta entrega A.bis. La fórmula 1D es trivial reducción de la 3D si se prioriza después.

- **2026-05-21** · **Implementación completa, status `validated`**. Clase + 12 tests unitarios + 1 test cruzado + 3 tests de integración (Hex8 elástico, Hex8 con daño activo + tangente asimétrica, Tet4 smoke). **16 tests añadidos a la suite (848 → 866 verdes; suite global sin regresiones).** Todos los tests pasaron a la primera, incluyendo:
  - Cruzado Damage3D ↔ Damage2D plane_strain a **14 decimales en `d` y `κ`** y 10 decimales en σ — la equivalencia bajo $\varepsilon_{zz} = \gamma_{yz} = \gamma_{xz} = 0$ es exacta a precisión máquina (misma ley centralizada, mismo $\mathbf C_e$ plane_strain bajo restricción cinemática).
  - Diferencia finita de la tangente algorítmica vs cerrada con error rel. < 1e-5 en las 36 componentes 6×6 en carga activa.
  - Asimetría significativa de la tangente en estados no degenerados (norma rel. ≥ 0.01).
  - Invariancia bajo cambio de unidades (E en MPa vs Pa) — `d`, `κ` idénticos a 12 decimales; σ escala linealmente con E.

  **Con esta entrega cierra la sub-etapa A.bis** (materiales 3D no lineales): VonMises3D + DruckerPrager3D + IsotropicDamage3D. Validación contra benchmark publicado (e.g. ensayo tensión-softening con valor cuantitativo de carga máxima en barra cargada axialmente con softening exponencial calibrado por G_f efectivo) **diferida a la campaña 3D consolidada** que se hará al cierre de A.bis junto con Hill 3D para VM3D y esfera hueca para DP3D.
- **2026-09-22** · Ensamblaje por lotes (ADR 0014): la física de `compute_state` se movió a un núcleo `@njit` (`_compute_damage_3d`) que consumen el camino por elemento y el adaptador por lotes; `evaluate_exponential_damage` se compila y recibe el cap como argumento. Sin cambio de formulación; `STATE_SCHEMA = {kappa, damage}` con estado inicial `κ = κ₀`.
- **2026-09-29** · **Cambio de formulación aprobado por el usuario** (segunda mitad del hallazgo clas-3 de la auditoría de pureza; la primera, la $\varepsilon_{zz}$ del 2D en tensión plana, en `IsotropicDamage2D.md`): se retira la norma del tensor $\sqrt{\boldsymbol\varepsilon:\boldsymbol\varepsilon}$, convención propia no estándar y mal calibrada (en tracción uniaxial de esfuerzo $\varepsilon_\text{eq} = \varepsilon\sqrt{1+2\nu^2}$, así que $\kappa_0 \ne f_t/E$). La deformación equivalente pasa a ser configurable con `equivalent_strain`: **energética** por defecto (Simó y Ju 1987; simétrica en tracción y compresión, como la anterior), y **Mazars** y **von Mises modificado** (este con `k = f_c/f_t` obligatorio) como opciones que distinguen tracción de compresión. Las tres dan $\varepsilon_\text{eq} = \varepsilon$ en tracción uniaxial ($\kappa_0 = f_t/E$). Ley de daño, variable histórica, `DAMAGE_MAX`, esfuerzo y ramas de la tangente sin cambios. Con la energética la tangente consistente es simétrica y la instancia declara `IS_SYMMETRIC = True` (§8). Las tres medidas viven en `solidum/materials/_equivalent_strain.py`, compartido con el 2D; el camino por lotes recibe la medida como código numérico en `batch_params`. El benchmark de validación uniaxial confinado se recalculó con el $\varepsilon_\text{eq}$ de cada medida bajo confinamiento (energética: $\varepsilon_{xx}\sqrt{(1-\nu)/[(1+\nu)(1-2\nu)]}$).
