# ISOTROPIC DAMAGE CRACK BAND 3D — daño escalar regularizado por banda de grieta, con integración implícita o IMPL-EX

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> **Variante de componente existente** (Reglas §4): subclase de `IsotropicDamage3D` ([spec](IsotropicDamage3D.md), validada), que no se modifica. Es la versión 3D de [IsotropicDamageCrackBand2D](IsotropicDamageCrackBand2D.md): **toda la física, las ecuaciones, los requisitos de interfaz y las preguntas al usuario (Q1–Q5) están allí**. Aquí sólo lo que difiere.
>
> **Estado: aprobada por el usuario (2026-11-11, junto con la 2D) y validada**: bloque 1 (banda de grieta) y bloque 2 (IMPL-EX, [ADR 0034](../adr/0034-pseudo-tiempo-de-los-solvers-incrementales.md)); ver la spec 2D.

---

## Especificación física

### 0. Descripción general

Daño isótropo escalar 3D con ablandamiento exponencial regularizado por banda de grieta y con integración implícita o IMPL-EX. `STRAIN_DIM = 6`, Voigt del proyecto $[\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, \gamma_{xy}, \gamma_{yz}, \gamma_{xz}]$ con $\gamma_{ij} = 2\varepsilon_{ij}$ (ADR 0012). Sin variantes de hipótesis.

### 1. Lo que es igual a la variante 2D

- Ley exponencial, $\kappa_0 = f_t/E$, regularización $\alpha(h) = 1/\bigl(G_f/(h f_t) - \kappa_0/2\bigr)$, límite $h < h_\text{max} = 2EG_f/f_t^2$ y energía $D(\kappa)$ (2D, §1–§2).
- Tope propio `damage_max` (2D, §3, Q5).
- Integración implícita e IMPL-EX, simetría y caveats (2D, §7–§10).
- Requisitos de interfaz: siembra de $h$ por punto de integración y razón de pseudo-tiempo de IMPL-EX (2D, §6).

### 2. Longitud característica (2D, §4 y Q2)

Convenio de Abaqus en 3D: $h = \sqrt[3]{V_e}$ en elementos lineales (`Hex8`, `Tet4`) y $h = \tfrac12\sqrt[3]{V_e}$ en cuadráticos (`Hex20`, `Hex27`, `Tet10`), con $V_e$ el volumen del elemento; el mismo $h$ en todos sus puntos de integración. Esos cinco elementos deben sembrarlo. Si la prueba de objetividad con `Tet4` se aparta más de un 5 % de $G_f\,A$, se decide su corrección igual que para `Tri3` en 2D.

**Resultado (2026-11-11)**: con $\sqrt[3]{V_e}$ la barra de `Tet4` (seis tetraedros por cubo) disipaba $1.81\,G_fA$; `Tet4` usa $h = \sqrt[3]{6V_e}$, el análogo del $\sqrt{2A_e}$ de Rots (el cubo que forman seis tetraedros), y disipa $0.999\,G_fA$, igual que `Hex8`.

**Cuadráticos (2026-11-11)**: con $\tfrac12\sqrt[3]{V_e}$ la barra disipaba $1.97\,G_fA$ (`Hex20`, `Hex27`) y $3.23\,G_fA$ (`Tet10`): la deformación se localiza en el elemento entero. Con la regla del 5 % (usuario, «si»; spec 2D §4) usan el factor del lineal de su celda, $\sqrt[3]{V_e}$ (`Hex20`, `Hex27`) y $\sqrt[3]{6V_e}$ (`Tet10`), y disipan $0.998$–$0.999\,G_fA$.

### 3. Medida de Rankine (2D, §5 y Q1)

$$\varepsilon_\text{eq} = \frac{\langle\bar\sigma_I\rangle_+}{E}, \qquad \bar{\boldsymbol\sigma} = \mathbf C_e\,\boldsymbol\varepsilon,$$

con $\bar\sigma_I$ el mayor valor propio del tensor $3\times3$ de esfuerzo efectivo y $\mathbf n_I$ su vector propio. En Voigt de esfuerzo,

$$\mathbf p = [n_x^2,\ n_y^2,\ n_z^2,\ 2n_xn_y,\ 2n_yn_z,\ 2n_xn_z], \qquad \frac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon} = \frac1E\,\mathbf C_e\,\mathbf p .$$

Con la variante 2D como restricción exacta de ésta en tensión plana y deformación plana (prueba `restriction_of_3d` de la spec 2D). Valores propios repetidos: subgradiente, como en 2D.

---

## Contrato de implementación

```yaml
name: IsotropicDamageCrackBand3D
kind: material
status: validated          # bloques 1 y 2 (2026-11-11)

interface:
  strain_dim: 6
  primary_state_var: damage
  is_symmetric: false      # cota de clase; la instancia es true con equivalent_strain='energy' o con integration='impl-ex'
  batch_kernel: true       # estado [κ, d, h, κ_prev, e]
  uses_pseudo_time: "true con integration='impl-ex' (ADR 0034)"

admission:                # Reglas §1; variante de un material admitido (IsotropicDamage3D); ver la spec 2D
  textbooks:
    - de_borst_crisfield: "§6.2, daño isótropo escalar (el padre); la banda de grieta en el texto, por confirmar"
    - de_souza_neto: "cap. 12, mecánica del daño (el padre); la banda de grieta en el texto, por confirmar"
  general_codes:
    - abaqus: "Concrete Damaged Plasticity: tracción regularizada por energía de fractura con la longitud característica del elemento"
    - diana: "Total Strain Crack Model con ancho de banda de grieta h"
    - code_aster: "ENDO_ISOT_BETON y MAZARS; METHODE='IMPLEX' de STAT_NON_LINE (doc. R5.03.81)"
    - opensees: "ASDConcrete3D: -autoRegularization y -implex"

parameters:
  - { name: E,          type: float, required: true,  desc: "Módulo de Young intacto [Pa]" }
  - { name: nu,         type: float, required: true,  desc: "Coeficiente de Poisson" }
  - { name: f_t,        type: float, required: true,  desc: "Resistencia a tracción [Pa]; κ_0 = f_t/E" }
  - { name: G_f,        type: float, required: true,  desc: "Energía de fractura en modo I [N/m]" }
  - { name: equivalent_strain, type: str, required: false, default: "rankine", desc: "rankine | energy | mazars | modified_von_mises (Q1 de la spec 2D)" }
  - { name: k,          type: float, required: false, default: null, desc: "f_c/f_t, sólo y obligatorio con modified_von_mises" }
  - { name: integration, type: str, required: false, default: "implicit", desc: "implicit | impl-ex (§8 de la spec 2D; ADR 0034)" }
  - { name: damage_max, type: float, required: false, default: 0.999999, desc: "Tope del daño de la variante (Q5 de la spec 2D)" }
  - { name: density,    type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(6,), state_vars) -> (σ: ndarray(6,), C_tan: ndarray(6,6), state_vars')"
  initial_point_state: "(h: float) -> state_vars; ValueError si h no es finito y positivo o si h ≥ 2·E·G_f/f_t²"
  strain_kind: "Voigt 6D — [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

state_schema:
  kappa:      "float ≥ κ_0 — máxima deformación equivalente histórica"
  damage:     "float ∈ [0, damage_max] — daño con que se calcula el esfuerzo"
  h:          "float > 0 — longitud característica del punto de integración"
  kappa_prev: "float — κ del paso confirmado anterior (IMPL-EX)"
  damage_error: "float ≥ 0 — indicador |d(κ) − d̃| de IMPL-EX; 0 en la implícita"

conventions:
  sign: "σ > 0 ⇔ tracción"
  voigt: "orden del proyecto (ADR 0012), γ_ij = 2·ε_ij"
  units: "E, f_t en Pa; G_f en N/m; h en m"

validity:
  - "E > 0, f_t > 0, G_f > 0, ν ∈ (−1, 0.5)"
  - "0 < h < 2·E·G_f/f_t² en cada punto"
  - "damage_max ∈ (0, 1)"
  - "integration = 'impl-ex' sólo con un solver que declare su pseudo-tiempo (ADR 0034)"

out_of_scope:
  - "lo mismo que la variante 2D"

acceptance:
  material_point:
    - name: same_as_2d
      setup: "los casos material_point e impl_ex de la spec 2D, en Voigt 6D"
      expect: "mismos criterios y tolerancias"
    - name: rankine_3d_values
      setup: "tracción uniaxial, triaxial igual, cortante puro, compresión uniaxial, biaxial y triaxial"
      expect: "ε_eq = ε en tracción uniaxial; inicio en σ = f_t en tracción triaxial igual y en τ = f_t en cortante puro; sin daño en compresión"
      tol_rel: 1.0e-12
    - name: consistent_tangent_fd
      setup: "carga activa, las cuatro medidas, estados generales con valores propios distintos"
      expect: "|C_FD − C_tan| ≤ 1e-6·max|C_FD|"
      tol_rel: 1.0e-6
  elements:
    - name: characteristic_length_seeded
      setup: "Hex8, Tet4, Hex20, Hex27, Tet10, malla regular y distorsionada"
      expect: "h de cada punto = ∛V_e (Hex8, Hex20, Hex27) o ∛(6V_e) (Tet4, Tet10)"
      tol_rel: 1.0e-12
  mesh_objectivity:
    - name: bar_in_tension_3d
      setup: "barra en tracción con debilitamiento en una capa; 1, 2, 4 y 8 elementos por unidad de longitud; Hex8 y Tet4"
      expect: "energía disipada = G_f·A_grieta ±1 % (Hex8 y, con ∛(6V_e), Tet4)"
      tol_rel: 1.0e-2

references:
  - "Las de IsotropicDamageCrackBand2D y las de IsotropicDamage3D."
```

---

## Implementación

- Archivo: `solidum/materials/damage_crack_band_3d.py` (común con la 2D en `solidum/materials/_crack_band.py`; Rankine 3D en `solidum/materials/_equivalent_strain.py`).
- Clase: `IsotropicDamageCrackBand3D(CrackBandMixin, IsotropicDamage3D)`, con el núcleo compilado del padre (`_compute_damage_3d`).
- Tests:
  - `tests/test_damage_crack_band.py` (punto material 3D, Rankine 3D, tangente por diferencias finitas, restricción exacta de la 2D, kernel por lotes; siembra de $h$ en `Hex8`, `Tet4`, `Hex20`, `Hex27`, `Tet10`).
  - `tests/validation/test_crack_band_objectivity.py` (barra 3D con `Hex8` y `Tet4`).
- Notas de traducción: ninguna.

---

## Diálogo

- **2026-11-11** · Spec creada por la IA junto con [IsotropicDamageCrackBand2D](IsotropicDamageCrackBand2D.md); las preguntas al usuario (Q1–Q5) están en esa spec y valen para las dos.
- **2026-11-11** · Aprobada por el usuario con la 2D («si»). Bloque 1 implementado (sin IMPL-EX); `Tet4` con $\sqrt[3]{6V_e}$ por la prueba de objetividad (§2). Hallazgos en el Diálogo de la spec 2D.
- **2026-11-11** · Bloque 2: IMPL-EX con la misma formulación que la 2D (§8 de esa spec; ADR 0034). Pruebas en 3D: tangente secante exacta y simétrica definida positiva, kernel por lotes bit a bit (implícita e IMPL-EX, razón variable y «sin paso») en `tests/test_damage_crack_band.py` y `tests/test_damage_crack_band_implex.py`. La prueba estructural (Van Vliet) es 2D.
- **2026-11-11** · Factor de $h$ de los cuadráticos corregido con la prueba de objetividad (§2; Diálogo de la spec 2D).
