# ISOTROPIC DAMAGE 2D — daño escalar con ablandamiento exponencial

> Orden de trabajo. El usuario escribe **especificación física**, **formulación** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> Spec **retroactiva con ampliación**: el modelo está implementado desde sesiones anteriores con **tangente secante**. Esta spec lo documenta y añade la **tangente algorítmica consistente** para recuperar convergencia cuadrática del Newton global. Desde el 2026-09-29 la **deformación equivalente es configurable** entre tres medidas de la literatura (ver §3 y *Diálogo*).

---

## Especificación física

### 0. Descripción general

Modelo de **daño isótropo escalar** en pequeñas deformaciones para sólidos 2D. Captura degradación progresiva e isótropa de la rigidez en respuesta a la deformación máxima histórica, con ablandamiento exponencial. Pensado para hormigón, cerámicos y materiales cuasi frágiles donde la microfisuración distribuida puede modelarse mediante un escalar de daño en lugar de microestructura explícita. Con las medidas de Mazars o de von Mises modificado (§3) el inicio del daño distingue tracción de compresión.

### 1. Cinemática

Pequeñas deformaciones, descomposición elástica pura: la deformación plástica no existe en este modelo (sin plasticidad). El daño se manifiesta como degradación de la rigidez efectiva.

### 2. Esfuerzo nominal vs esfuerzo efectivo

$$\boldsymbol\sigma = (1 - d)\,\boldsymbol\sigma^\text{eff}, \qquad \boldsymbol\sigma^\text{eff} = \mathbf C_e\,\boldsymbol\varepsilon$$

con $d \in [0, d_\text{max}]$ la variable de daño isótropa y $\mathbf C_e$ el tensor constitutivo elástico isótropo (2D, plane stress o plane strain según `hypothesis` heredada del `Elastic2D` interno).

El esfuerzo efectivo $\boldsymbol\sigma^\text{eff}$ es el que el material intacto soportaría con la deformación corriente. El esfuerzo nominal $\boldsymbol\sigma$ es lo que efectivamente transmite considerando la degradación.

### 3. Deformación equivalente

Cantidad escalar que cuantifica la "intensidad" de la deformación y gobierna el daño. Se elige con el parámetro `equivalent_strain` entre tres medidas estándar. Todas se evalúan sobre el **estado 3D asociado** al estado plano, en notación Voigt $[\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]$ con $\gamma_{xy} = 2\varepsilon_{xy}$:

$$\boldsymbol\varepsilon^\text{3D} = [\varepsilon_{xx},\ \varepsilon_{yy},\ \varepsilon_{zz},\ \gamma_{xy},\ 0,\ 0],
\qquad
\varepsilon_{zz} = c_{zz}\,(\varepsilon_{xx} + \varepsilon_{yy}),\quad
c_{zz} = \begin{cases} -\dfrac{\nu}{1-\nu} & \text{tensión plana} \\[4pt] 0 & \text{deformación plana} \end{cases}$$

En tensión plana $\varepsilon_{zz}$ no depende del daño: con daño escalar isótropo, $\sigma_{zz} = (1-d)\,(\lambda\,\mathrm{tr}\,\boldsymbol\varepsilon + 2\mu\,\varepsilon_{zz}) = 0$ para cualquier $d$. La dirección $z$ es principal ($\gamma_{xz} = \gamma_{yz} = 0$). Con cualquiera de las tres medidas el 2D es la **restricción exacta** de `IsotropicDamage3D` en las dos hipótesis.

Las tres medidas están **calibradas en tracción uniaxial de esfuerzo**: $\varepsilon_\text{eq} = \varepsilon$, de modo que el daño empieza en $\sigma = E\,\kappa_0$ y $\kappa_0 = f_t/E$. Se diferencian en compresión y en estados multiaxiales.

1. **Energética** (`energy`, por defecto; Simó y Ju 1987):
   $$\varepsilon_\text{eq} = \sqrt{\frac{\boldsymbol\varepsilon^\top\mathbf C_e\,\boldsymbol\varepsilon}{E}}$$
   con $\mathbf C_e$ la matriz elástica Voigt del propio material (la de tensión plana o la de deformación plana). Con la de tensión plana, $\boldsymbol\varepsilon^\top\mathbf C_\text{ps}\boldsymbol\varepsilon = \boldsymbol\sigma:\boldsymbol\varepsilon$ del 3D con $\sigma_{zz} = 0$; en deformación plana el término $\sigma_{zz}\varepsilon_{zz}$ es nulo. Por eso esta medida no necesita tratar $\varepsilon_{zz}$ aparte. Simétrica en tracción y compresión (en compresión uniaxial de esfuerzo, $\varepsilon_\text{eq} = |\varepsilon|$).

2. **Mazars** (`mazars`; Mazars 1986):
   $$\varepsilon_\text{eq} = \sqrt{\textstyle\sum_{i=1}^{3}\langle\varepsilon_i\rangle_+^2}$$
   con $\varepsilon_i$ las deformaciones principales del tensor 3D asociado: las dos del plano y $\varepsilon_{zz}$. Sólo las extensiones dañan: en compresión uniaxial de esfuerzo daña la extensión lateral, $\varepsilon_\text{eq} = \sqrt2\,\nu\,|\varepsilon|$.

3. **von Mises modificado** (`modified_von_mises`; de Vree, Brekelmans y van Gils 1995), con $k = f_c/f_t$:
   $$\varepsilon_\text{eq} = \frac{k-1}{2k(1-2\nu)}\,I_1 + \frac{1}{2k}\sqrt{\left(\frac{k-1}{1-2\nu}\right)^2 I_1^2 + \frac{12k}{(1+\nu)^2}\,J_2}$$
   con $I_1 = \varepsilon_{xx} + \varepsilon_{yy} + \varepsilon_{zz}$ y $J_2 = \tfrac16\left[(\varepsilon_{xx}-\varepsilon_{yy})^2 + (\varepsilon_{yy}-\varepsilon_{zz})^2 + (\varepsilon_{zz}-\varepsilon_{xx})^2\right] + \tfrac14\gamma_{xy}^2$. En compresión uniaxial de esfuerzo, $\varepsilon_\text{eq} = |\varepsilon|/k$.

Hasta el 2026-09-29 la medida era la norma del tensor, $\varepsilon_\text{eq} = \sqrt{\boldsymbol\varepsilon:\boldsymbol\varepsilon}$: convención propia, no estándar y mal calibrada (en tracción uniaxial de esfuerzo daba $\varepsilon_\text{eq} = \varepsilon\sqrt{1+2\nu^2}$, así que $\kappa_0 \ne f_t/E$). Se retiró (ver *Diálogo*).

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

con saturación numérica $d \leq d_\text{max} = $ `DAMAGE_MAX` (`solidum.constants.DAMAGE_MAX`, evita rigidez nula).

Propiedades:
- $d(\kappa_0) = 0$ (continuidad).
- $d'(\kappa) > 0$ para $\kappa > \kappa_0$ (irreversibilidad).
- $\lim_{\kappa \to \infty} d = 1$ (asintótico).
- $\alpha$ controla la velocidad de degradación (ablandamiento más abrupto si $\alpha$ grande).

---

## Formulación numérica

### 6. Algoritmo en un paso

Dado $\{\kappa_n\}$ y $\boldsymbol\varepsilon_{n+1}$:

1. $\boldsymbol\sigma^\text{eff} = \mathbf C_e\,\boldsymbol\varepsilon_{n+1}$ y $\varepsilon_\text{eq}$ con la medida elegida (§3), junto con su derivada $\partial\varepsilon_\text{eq}/\partial\boldsymbol\varepsilon$.
2. Si $\varepsilon_\text{eq} > \kappa_n$: **carga activa**, $\kappa_{n+1} = \varepsilon_\text{eq}$. Si no: **descarga**, $\kappa_{n+1} = \kappa_n$.
3. Aplicar ley de daño con $\kappa_{n+1}$; saturar a $d_\text{max}$.
4. $\boldsymbol\sigma_{n+1} = (1 - d)\,\mathbf C_e\,\boldsymbol\varepsilon_{n+1}$.
5. Tangente: **secante** en descarga / sin daño / saturación; **algorítmica consistente** en carga activa con daño no saturado (ver §7).

El algoritmo es **explícito**: no requiere Newton local. La actualización de $\kappa$ y $d$ es directa. Las deformaciones principales de Mazars salen en forma cerrada: las dos del plano del bloque $2\times2$ y la tercera es $\varepsilon_{zz}$.

### 7. Tangente algorítmica consistente (ampliación)

$$\mathbf C^\text{alg} = \frac{\partial \boldsymbol\sigma_{n+1}}{\partial \boldsymbol\varepsilon_{n+1}}$$

Derivando $\boldsymbol\sigma = (1-d)\,\mathbf C_e\,\boldsymbol\varepsilon$:

$$\mathbf C^\text{alg} = (1-d)\,\mathbf C_e - \boldsymbol\sigma^\text{eff} \otimes \frac{\partial d}{\partial \boldsymbol\varepsilon},
\qquad
\frac{\partial d}{\partial \boldsymbol\varepsilon} = \frac{\partial d}{\partial \kappa}\,\frac{\partial \kappa}{\partial \boldsymbol\varepsilon}$$

donde:

- $\dfrac{\partial d}{\partial \kappa} = (1 - d)\left(\dfrac{1}{\kappa} + \alpha\right)$ para la ley exponencial. (Demostración: derivando $d = 1 - (\kappa_0/\kappa)e^{-\alpha(\kappa-\kappa_0)}$ y reagrupando vía $1 - d = (\kappa_0/\kappa)e^{-\alpha(\kappa-\kappa_0)}$.)
- $\dfrac{\partial \kappa}{\partial \boldsymbol\varepsilon} = \dfrac{\partial \varepsilon_\text{eq}}{\partial \boldsymbol\varepsilon}$ en **carga activa** y $\mathbf 0$ en **descarga**.

Derivadas de la deformación equivalente respecto de la deformación en Voigt *engineering*. En las medidas que usan el estado 3D asociado, la regla de la cadena por $\varepsilon_{zz} = c_{zz}(\varepsilon_{xx} + \varepsilon_{yy})$ añade $c_{zz}\,\partial\varepsilon_\text{eq}/\partial\varepsilon_{zz}$ a las componentes $xx$ e $yy$:

- **Energética**: $\dfrac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon} = \dfrac{\mathbf C_e\,\boldsymbol\varepsilon}{E\,\varepsilon_\text{eq}}$ (la cadena de $\varepsilon_{zz}$ ya está en $\mathbf C_\text{ps}$).
- **Mazars**: $\dfrac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon} = \dfrac{1}{\varepsilon_\text{eq}}\displaystyle\sum_i \langle\varepsilon_i\rangle_+\,\frac{\partial\varepsilon_i}{\partial\boldsymbol\varepsilon}$, con $\partial\varepsilon_i/\partial\boldsymbol\varepsilon$ el proyector $\mathbf n_i\otimes\mathbf n_i$ en Voigt *engineering* (componente de $\gamma_{ab}$: $n_a n_b$). La suma es la parte positiva del tensor, $\boldsymbol\varepsilon_+ = \sum_i\langle\varepsilon_i\rangle_+\mathbf n_i\otimes\mathbf n_i$: en 2D, $[(\varepsilon_+)_{xx} + c_{zz}\langle\varepsilon_{zz}\rangle_+,\ (\varepsilon_+)_{yy} + c_{zz}\langle\varepsilon_{zz}\rangle_+,\ (\varepsilon_+)_{xy}]/\varepsilon_\text{eq}$. Está bien definida con deformaciones principales repetidas, porque $\sum_i\langle\lambda_i\rangle_+^2$ es una función simétrica $C^1$ de los valores principales.
- **von Mises modificado**: con $S = \sqrt{b^2 I_1^2 + c\,J_2}$, $b = (k-1)/(1-2\nu)$, $c = 12k/(1+\nu)^2$:
  $$\frac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon^\text{3D}} = \frac{b + b^2 I_1/S}{2k}\,\frac{\partial I_1}{\partial\boldsymbol\varepsilon^\text{3D}} + \frac{c}{4kS}\,\frac{\partial J_2}{\partial\boldsymbol\varepsilon^\text{3D}},
  \qquad \frac{\partial I_1}{\partial\boldsymbol\varepsilon^\text{3D}} = [1, 1, 1, 0, 0, 0],\quad
  \frac{\partial J_2}{\partial\boldsymbol\varepsilon^\text{3D}} = [e_{xx}, e_{yy}, e_{zz}, \tfrac12\gamma_{xy}, \tfrac12\gamma_{yz}, \tfrac12\gamma_{xz}]$$
  ($e$ la parte desviadora) y la cadena de $\varepsilon_{zz}$ para pasar al 2D.

**Tangente final**:

$$\mathbf C^\text{alg} = \begin{cases}
\mathbf C_e & \kappa \leq \kappa_0 \;(\text{sin daño}) \\
(1-d)\,\mathbf C_e & \text{descarga} \\
(1-d)\,\mathbf C_e \;-\; (1-d)\left(\dfrac1\kappa + \alpha\right)(\mathbf C_e\boldsymbol\varepsilon) \otimes \dfrac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon} & \text{carga activa, } d < d_\text{max} \\
(1-d_\text{max})\,\mathbf C_e & d = d_\text{max} \;(\text{saturación})
\end{cases}$$

### 8. Simetría de la tangente

- **Energética**: $\partial\varepsilon_\text{eq}/\partial\boldsymbol\varepsilon \parallel \mathbf C_e\boldsymbol\varepsilon$, así que el término de carga activa es $\propto(\mathbf C_e\boldsymbol\varepsilon)\otimes(\mathbf C_e\boldsymbol\varepsilon)$: la tangente consistente es **simétrica** (en la implementación, bit a bit).
- **Mazars y von Mises modificado**: $\mathbf C_e\boldsymbol\varepsilon$ y $\partial\varepsilon_\text{eq}/\partial\boldsymbol\varepsilon$ no son proporcionales en general: tangente **no simétrica**.

Consecuencia algebraica (ADR 0003): `IS_SYMMETRIC` es `False` en la clase (cota conservadora) y la **instancia** lo pone a `True` con la medida energética; `domain_is_symmetric` lo lee de la instancia (mismo patrón que `DruckerPrager2D` asociado). Con Mazars o von Mises modificado el despachador elige LU (o Pardiso). Con la energética puede elegir Cholesky si está disponible; como la tangente de ablandamiento puede no ser definida positiva, el corrector degrada a LU en cuanto Cholesky la rechaza.

### 9. Decisión de régimen durante el Newton

Durante una iteración del Newton global, el material recibe un $\boldsymbol\varepsilon$ trial y debe decidir carga/descarga comparando con $\kappa_n$ (el estado committed del paso anterior). La transición entre ramas (loading ↔ unloading) es **discontinua en la tangente**: al cruzar $\varepsilon_\text{eq} = \kappa_n$ el segundo término aparece/desaparece bruscamente. En la práctica esto no genera problemas porque el Newton converge antes de oscilar entre ramas; en problemas patológicos puede requerir continuation o arc-length.

---

## Contrato de implementación

```yaml
name: IsotropicDamage2D
kind: material
status: validated

interface:
  strain_dim: 3
  primary_state_var: damage
  is_symmetric: false      # cota de clase; la instancia es true con equivalent_strain='energy' (ver §8)

admission:                # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - de_borst_crisfield: "§6.2 (daño isótropo escalar con ablandamiento exponencial y medidas de deformación equivalente)"
    - de_souza_neto: "cap. 12 (mecánica del daño)"
  general_codes:
    - code_aster: "MAZARS, ENDO_ISOT_BETON (daño escalar isótropo; MAZARS con la deformación equivalente de Mazars)"

parameters:
  - { name: E,          type: float, required: true,  desc: "Módulo de Young intacto" }
  - { name: nu,         type: float, required: true,  desc: "Coeficiente de Poisson" }
  - { name: kappa_0,    type: float, required: true,  desc: "Umbral de deformación equivalente; daño nulo mientras κ ≤ κ_0. Con las tres medidas κ_0 = f_t/E" }
  - { name: alpha,      type: float, required: true,  desc: "Velocidad de degradación exponencial (>0)" }
  - { name: hypothesis, type: str,   required: false, default: "plane_stress", desc: "plane_stress | plane_strain (delegado a Elastic2D interno)" }
  - { name: equivalent_strain, type: str, required: false, default: "energy", desc: "energy (Simó-Ju 1987) | mazars (Mazars 1986) | modified_von_mises (de Vree et al. 1995)" }
  - { name: k,          type: float, required: false, default: null, desc: "f_c/f_t (>0); obligatorio con modified_von_mises y no admitido con las otras medidas" }
  - { name: density,    type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(3,), state_vars=None) -> (σ: ndarray(3,), C_tan: ndarray(3,3), state_vars')"
  strain_kind: "plano Voigt — [ε_xx, ε_yy, γ_xy] con γ_xy = 2·ε_xy"

state_schema:
  kappa:  "float ≥ κ_0 — máxima deformación equivalente histórica"
  damage: "float ∈ [0, DAMAGE_MAX] — variable de daño escalar isótropa"

conventions:
  sign: "σ > 0 ⇔ tracción; energy no distingue tracción de compresión; mazars sólo daña con extensiones principales; modified_von_mises daña en compresión uniaxial con ε_eq = |ε|/k"
  voigt: "γ_xy = 2·ε_xy en deformación; ε_eq sobre el estado 3D [ε_xx, ε_yy, ε_zz, γ_xy, 0, 0] con ε_zz = c_zz·(ε_xx + ε_yy): c_zz = −ν/(1−ν) en tensión plana, 0 en deformación plana"
  equivalent_strain: "energy: √(εᵀ·C_e·ε/E) con la C_e plana del material; mazars: √(Σ⟨ε_i⟩₊²) sobre las tres principales; modified_von_mises: (k−1)/(2k(1−2ν))·I1 + 1/(2k)·√(((k−1)/(1−2ν))²·I1² + 12k/(1+ν)²·J2). Las tres dan ε_eq = ε en tracción uniaxial de esfuerzo"
  damage_irreversibility: "κ monótono no decreciente vía max(κ_old, ε_eq); d(κ) creciente"

validity:
  - "E > 0, κ_0 > 0, α > 0"
  - "ν ∈ (-1, 0.5) (delegación a Elastic2D)"
  - "hypothesis ∈ {'plane_stress', 'plane_strain'}"
  - "equivalent_strain ∈ {'energy', 'mazars', 'modified_von_mises'}"
  - "k > 0 finito, sólo con equivalent_strain = 'modified_von_mises' (obligatorio ahí)"

out_of_scope:
  - "modelo de Mazars completo con dos daños (d_t, d_c) ponderados por α_t, α_c; aquí sólo su medida de deformación equivalente con un único daño escalar"
  - "anisotropía del daño (tensor de daño D vs escalar)"
  - "regularización para evitar mesh-dependency en regime de ablandamiento (longitud característica, no-local, gradient)"
  - "acoplamiento con plasticidad"
  - "fatiga / acumulación cíclica del daño"
  - "restitución de rigidez al cerrar fisuras (efecto unilateral) — el daño escalar degrada igual en tracción y compresión una vez activado"

numerical_caveats:
  - "IS_SYMMETRIC por instancia: true con energy (tangente consistente simétrica), false con mazars y modified_von_mises (LU/Pardiso, ADR 0003). Con energy y Cholesky disponible, la tangente de ablandamiento puede no ser definida positiva: el corrector degrada a LU al primer rechazo."
  - "Transición loading↔unloading es discontinua en la tangente. Newton converge antes de oscilar entre ramas en casos típicos. Para problemas patológicos cerca del umbral, considerar arc-length."
  - "Mazars: ε_eq es C¹ pero no C² donde una deformación principal (incluida ε_zz en esfuerzo plano) cruza cero, es decir, donde entra o sale una extensión. ∂ε_eq/∂ε y la tangente consistente son continuas en ese cruce; lo que salta es su derivada (la curvatura de ⟨ε_i⟩₊²), y la tangente sigue siendo Lipschitz, así que el Newton conserva la convergencia cuadrática. Con deformaciones principales repetidas la derivada está bien definida."
  - "Saturación d=DAMAGE_MAX corta la corrección consistente (∂d/∂κ deja de ser efectivamente (1-d)(1/κ+α)); el algoritmo devuelve tangente secante en ese caso para evitar términos espurios."
  - "Sin regularización, la solución es mesh-dependent en régimen de ablandamiento (localización en una banda de elementos). Para benchmark cuantitativo se requiere mismo refinamiento de malla; para comparaciones cualitativas (presencia de banda, dirección) basta con regularización de gradiente o longitud característica (fuera de scope)."

acceptance:
  # tests/test_materials_unit.py · TestIsotropicDamage2D, TestIsotropicDamage2DConsistentTangent;
  # tests/test_damage_equivalent_strain.py (tres medidas); tests/test_damage_2d_3d_consistency.py;
  # tests/test_solid_2d_damage.py (integración).
  unit_existing_no_regression:
    - name: secant_response_below_threshold
      expect: "κ ≤ κ_0 ⇒ d=0, σ = C_e·ε, tangente = C_e"
      tol_rel: 1.0e-12

    - name: damage_evolution_exponential
      expect: "para κ > κ_0, d cumple ley exponencial con ε_eq de la medida (energy: ε_eq = ε_xx/√(1−ν²) para ε = [ε_xx, 0, 0] en tensión plana)"
      tol_rel: 1.0e-10

  unit_equivalent_strain:
    - name: eps_eq_matches_independent_reference
      setup: "25 estados aleatorios por medida e hipótesis; referencia sobre el tensor 3D (σ:ε con Lamé, autovalores de numpy, J2 = ½·e:e)"
      expect: "κ tras un paso desde κ = 0 coincide con ε_eq de referencia"
      tol_rel: 1.0e-12

    - name: uniaxial_stress_calibration
      setup: "tracción y compresión uniaxial de esfuerzo en tensión plana (ε = [ε, −νε, 0])"
      expect: "tracción: ε_eq = ε con las tres medidas (κ_0 = f_t/E); compresión: energy |ε|, mazars √2·ν·|ε|, modified_von_mises |ε|/k"
      tol_rel: 1.0e-13

    - name: damage_onset_at_kappa_0
      expect: "d = 0 con ε_xx = 0,999·κ_0 y d > 0 con ε_xx = 1,001·κ_0 en tracción uniaxial de esfuerzo, con las tres medidas"

    - name: rotational_invariance
      expect: "ε_eq invariante frente a giros en el plano con las tres medidas"
      tol_rel: 1.0e-12

    - name: constructor_validation
      expect: "ValueError con el nombre del material si equivalent_strain es desconocida, si falta k con modified_von_mises, si se da k con otra medida o si k ≤ 0 / no finito"

  unit_consistent_tangent:
    - name: tangent_equals_secant_on_unloading
      setup: "carga hasta κ>κ_0, descarga con ε de menor ε_eq"
      expect: "C_tan = (1-d)·C_e exacto (segundo término se anula con ∂κ/∂ε = 0)"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_below_threshold
      setup: "ε pequeño, κ_old = κ_0"
      expect: "C_tan = C_e (sin daño activo, ni secante con d=0 ni corrección)"
      tol_rel: 1.0e-12

    - name: tangent_equals_secant_at_saturation
      setup: "carga hasta κ enorme tal que d alcanza DAMAGE_MAX"
      expect: "C_tan = (1-DAMAGE_MAX)·C_e (el término consistente se corta porque ∂d/∂κ efectivamente 0 al saturar)"
      tol_rel: 1.0e-12

    - name: tangent_consistent_term_on_loading
      setup: "carga activa con κ > κ_0 y d < DAMAGE_MAX"
      expect: "C_tan ≠ (1-d)·C_e — incluye el término rank-1 −(1-d)(1/κ+α)·(C_e·ε) ⊗ ∂ε_eq/∂ε (§7)"
      tol_rel: 1.0e-9

    - name: tangent_symmetry_in_loading
      setup: "carga activa con estado no degenerado (ε_xx ≠ ε_yy, γ_xy ≠ 0)"
      expect: "energy: C_tan = C_tanᵀ bit a bit; mazars y modified_von_mises: ‖C_tan − C_tanᵀ‖_F/‖C_tan‖_F > 1e-3"

    - name: tangent_finite_difference_consistency
      setup: "carga activa desde κ = κ_0, con cada medida, en tensión y deformación plana; estados general, todo extensión y principales del plano repetidas; diferencia centrada con h = 1e-10 sobre |ε| ~ 1e-4"
      expect: "|C_FD − C_tan| ≤ 1e-6·max|C_FD| (truncamiento ~1e-12 y redondeo ~1e-10 relativos)"
      tol_rel: 1.0e-6

  unit_arch:
    - name: is_symmetric_attribute
      expect: "IsotropicDamage2D.IS_SYMMETRIC == False en la clase; instancia True con energy, False con mazars y modified_von_mises"

    - name: batch_kernel_matches_compute_state
      expect: "el kernel por lotes reproduce bit a bit σ, C_tan y la fila de estado de compute_state con las tres medidas (carga, descarga, recarga, compresión)"

  cross_consistency:
    - name: restriction_of_isotropic_damage_3d
      setup: "historia de carga, descarga y recarga en cuatro direcciones; tensión plana (ε_zz = c_zz·(ε_xx+ε_yy)) y deformación plana (ε_zz = 0); las tres medidas"
      expect: "κ, d y σ en el plano iguales a IsotropicDamage3D con el estado 3D asociado; σ_zz 3D = 0 en tensión plana"
      tol_rel: 1.0e-12

  integration:
    - name: quad4_damage_newton_converges_in_few_iter
      setup: "Quad4 unitario, tensión plana, tracción uniaxial de esfuerzo por desplazamiento prescrito δ = 4·κ_0 en 5 pasos, sin bisección; las tres medidas"
      expect: "Newton converge en ≤6 iteraciones por paso; d y σ_xx = (1−d)·E·δ analíticos (ε_eq = δ)"
      tol_rel: 1.0e-8

    - name: yaml_equivalent_strain
      setup: "Quad4 con mazars y con la medida por defecto, Hex8 (IsotropicDamage3D) con modified_von_mises y k, leídos desde YAML"
      expect: "parámetros leídos; reacción (1−d)·E·δ analítica; modified_von_mises sin k rechazado al leer"

references:
  - "Lemaitre J., Chaboche J.-L. (1990). Mechanics of Solid Materials. Cambridge UP. Cap. 7 (continuum damage mechanics)."
  - "Oliver J. (1989). A consistent characteristic length for smeared cracking models. IJNME 28, 461-474."
  - "Simó J.C., Ju J.W. (1987). Strain- and stress-based continuum damage models — I. Formulation. Int. J. Solids Struct. 23, 821-840 (deformación equivalente energética; tangente consistente para daño isótropo)."
  - "Mazars J. (1986). A description of micro- and macroscale damage of concrete structures. Eng. Fract. Mech. 25, 729-737 (deformación equivalente de Mazars)."
  - "de Vree J.H.P., Brekelmans W.A.M., van Gils M.A.J. (1995). Comparison of nonlocal approaches in continuum damage mechanics. Comput. Struct. 55, 581-588 (deformación equivalente de von Mises modificado)."
  - "de Borst R., Crisfield M.A., Remmers J.J.C., Verhoosel C.V. (2012). Nonlinear Finite Element Analysis of Solids and Structures, 2.ª ed. Wiley. Cap. 6 (daño)."
  - "de Souza Neto E.A., Perić D., Owen D.R.J. (2008). Computational Methods for Plasticity. Wiley. §12 (damage models)."
```

---

## Implementación

- **Archivo**: [solidum/materials/damage_2d.py](../../solidum/materials/damage_2d.py); medidas de deformación equivalente en [solidum/materials/_equivalent_strain.py](../../solidum/materials/_equivalent_strain.py), compartidas con `IsotropicDamage3D`.
- **Clase**: `IsotropicDamage2D`, registrada vía `@MaterialRegistry.register`. `IS_SYMMETRIC = False` como ClassVar; el constructor pone `self.IS_SYMMETRIC = True` con `equivalent_strain='energy'`.
- **Estado**: dict `{'kappa': float, 'damage': float}`. `PRIMARY_STATE_VAR = 'damage'`.
- **Constructor**: valida `equivalent_strain` y `k` con `parse_equivalent_strain` (mensajes con el nombre del material y del parámetro) y guarda `equivalent_strain`, `k` (`None` salvo en von Mises modificado) y `c_zz`.
- **Algoritmo `compute_state`**: la física vive en el núcleo `@njit` `_compute_damage_2d`, que consumen `compute_state` y el adaptador por lotes (ADR 0014):
  1. $\boldsymbol\sigma^\text{eff} = \mathbf C_e\boldsymbol\varepsilon$ y `equivalent_strain_2d` → $(\varepsilon_\text{eq}, \mathbf v, s)$ con $\partial\varepsilon_\text{eq}/\partial\boldsymbol\varepsilon = s\,\mathbf v$. Con la energética $\mathbf v = \boldsymbol\sigma^\text{eff}$ y $s = 1/(E\varepsilon_\text{eq})$: el producto externo $(\boldsymbol\sigma^\text{eff})\otimes\mathbf v$ es simétrico bit a bit.
  2. Kuhn-Tucker con $\kappa_n$, ley exponencial, saturación a `DAMAGE_MAX`.
  3. $\boldsymbol\sigma = (1-d)\mathbf C_e\boldsymbol\varepsilon$ y tangente según §7.
- **Mazars en 2D**: forma cerrada del bloque $2\times2$ del plano (casos todo extensión, todo contracción y mixto, este último con $(\mathbf A - e_2\mathbf I)/(e_1 - e_2)$ evaluado sin cancelación) más $\varepsilon_{zz}$. **von Mises modificado en 2D**: sobre el estado 3D asociado y regla de la cadena.
- **Camino por lotes**: `batch_params = [E, κ_0, α, tol_abs, tol_rel, DAMAGE_MAX, c_zz, ν, código de la medida, k]` (códigos 0 energy, 1 mazars, 2 modified_von_mises; `k = 0` sin uso fuera de von Mises modificado). Sigue siendo cacheable en disco (verificado: el kernel de familia y los núcleos del material cargan del caché de Numba sin recompilar).
- **`admissibility_scale`**: `E·κ_0` (umbral inicial — la escala no evoluciona con el estado porque el criterio se mide contra $\kappa_0$, no contra $\kappa$).
- **Tests**:
  - [tests/test_materials_unit.py](../../tests/test_materials_unit.py) · `TestIsotropicDamage2D` y `TestIsotropicDamage2DConsistentTangent`.
  - [tests/test_damage_equivalent_strain.py](../../tests/test_damage_equivalent_strain.py): las tres medidas (valores, calibración, tangente por diferencias finitas, lotes, simetría, constructor, YAML).
  - [tests/test_damage_2d_3d_consistency.py](../../tests/test_damage_2d_3d_consistency.py): restricción exacta del 3D con las tres medidas.
  - [tests/test_solid_2d_damage.py](../../tests/test_solid_2d_damage.py): integración Quad4 + NonlinearSolver.

---

## Diálogo

- **2026-05-14** · Spec creada retroactivamente sobre material existente (tangente secante) y extendida con tangente algorítmica consistente. Sin ADR: extiende un material existente añadiendo un término a la tangente; no rompe contratos arquitecturales transversales. Sí cambia `IS_SYMMETRIC` para este material, lo cual el despachador algebraico ya soporta vía agregación sobre el dominio (no requiere refactor).
- **2026-05-14** · Decisión: la transición loading↔unloading se decide con $\kappa_n$ committed (no con $\kappa$ trial dentro del Newton). Esto es lo estándar y evita oscilaciones del Newton entre ramas en pasos con ambigüedad. Si el paso converge, $\kappa$ se compromete y la decisión es consistente con la trayectoria.
- **2026-05-14** · `IsotropicDamage1D` queda con tangente secante. Misma deuda gemela conceptual, pero fuera del alcance pactado (A2 era solo 2D). Si el usuario la prioriza luego, replicar es trivial: misma fórmula reducida a escalar.
- **2026-09-22** · Auditoría global: nuevo `out_of_plane_stress` delegado al elástico base (σ = (1−d)·C_e·ε escala todas las componentes, así que σ_zz = ν(σ_xx + σ_yy) sigue valiendo).
- **2026-09-22** · Ensamblaje por lotes (ADR 0014): la física de `compute_state` se movió a un núcleo `@njit` (`_compute_damage_2d`) que consumen el camino por elemento y el adaptador por lotes; `evaluate_exponential_damage` se compila y recibe el cap como argumento. Sin cambio de formulación; `STATE_SCHEMA = {kappa, damage}` con estado inicial `κ = κ₀`.
- **2026-09-22** · Cierre del ADR 0014 (§9): `batch_out_of_plane_stress(sigma, S)`, versión por arreglos de `out_of_plane_stress` para el post-proceso por familia (exportador VTK); misma fórmula, `ε^p_zz` leída de la fila de estado cuando aplica. Verificado archivo VTK idéntico por ambos caminos.
- **2026-09-29** · Cambio de formulación aprobado por el usuario (hallazgo clas-3 de la auditoría de pureza): en tensión plana $\varepsilon_\text{eq}$ incluye $\varepsilon_{zz} = -\nu/(1-\nu)\,(\varepsilon_{xx}+\varepsilon_{yy})$, que con daño escalar isótropo no depende de $d$. Antes se omitía y, con el mismo estado de esfuerzo, el 2D retrasaba el inicio del daño respecto de `IsotropicDamage3D` (con $E = 30$ GPa, $\nu = 0{,}2$, $\kappa_0 = 10^{-4}$: $d = 0$ en 2D frente a $d = 1{,}4\,\%$ en 3D). La tangente consistente incluye el término de $\varepsilon_{zz}$ en $\partial\varepsilon_\text{eq}/\partial\boldsymbol\varepsilon$ (verificada por diferencias finitas). Blindado en `tests/test_damage_2d_3d_consistency.py`: el 2D reproduce al 3D en toda una historia de carga, descarga y recarga, en tensión y en deformación plana, y el test falla si se retira el término. Queda abierta la otra mitad del hallazgo: si la norma del tensor es la medida adecuada o conviene una de la literatura (Mazars, von Mises modificado); pendiente del usuario.
- **2026-09-29** · Segunda mitad del hallazgo, **cambio aprobado por el usuario**: se retira la norma del tensor $\sqrt{\boldsymbol\varepsilon:\boldsymbol\varepsilon}$, convención propia no estándar y mal calibrada (en tracción uniaxial de esfuerzo $\varepsilon_\text{eq} = \varepsilon\sqrt{1+2\nu^2}$, así que $\kappa_0 \ne f_t/E$). La deformación equivalente pasa a ser configurable con `equivalent_strain`: **energética** por defecto (Simó y Ju 1987, la referencia que ya citaba la spec; simétrica en tracción y compresión, como la anterior), y **Mazars** y **von Mises modificado** (este con `k = f_c/f_t` obligatorio) como opciones que distinguen tracción de compresión. Las tres dan $\varepsilon_\text{eq} = \varepsilon$ en tracción uniaxial ($\kappa_0 = f_t/E$). Ley de daño, variable histórica, `DAMAGE_MAX`, esfuerzo y ramas de la tangente sin cambios; `IsotropicDamage1D` tampoco cambia ($|\varepsilon|$ es la restricción 1D de las tres). Consecuencias: (i) con la energética la tangente consistente es simétrica y la instancia declara `IS_SYMMETRIC = True` (§8), lo que habilita Cholesky cuando está instalado, con la degradación a LU ya existente si la tangente de ablandamiento no es definida positiva; (ii) las medidas viven en `solidum/materials/_equivalent_strain.py`, compartido con el 3D (dos casos reales, Reglas §1); (iii) el camino por lotes recibe la medida como código numérico en `batch_params`. Tests existentes recalculados con la medida nueva a partir de las fórmulas; `tests/test_solid_2d_damage.py` pasa a desplazamiento prescrito, porque con fuerza por encima del pico ($\sigma \approx E\kappa_0$) sólo había equilibrio con el daño saturado y el paso adaptativo ocultaba las iteraciones.
