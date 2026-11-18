# TENSION COMPRESSION DAMAGE 3D — daño de tracción y de compresión con cierre de grieta, banda de grieta e IMPL-EX

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> Versión 3D de [TensionCompressionDamage2D](TensionCompressionDamage2D.md): **toda la física, las ecuaciones y las preguntas al usuario (Q1–Q6) están allí**. Aquí sólo lo que difiere.
>
> **Estado: spec completa, junto con la 2D: Q1–Q4 y Q6 decididas el 2026-11-17 (ley (a) de Faria, Oliver y Cervera); Q5 el 2026-11-18 (programa principal).**

---

## Especificación física

### 0. Descripción general

Dos daños escalares, $d^+$ de tracción y $d^-$ de compresión, con cierre de grieta, en 3D. `STRAIN_DIM = 6`, Voigt del proyecto $[\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, \gamma_{xy}, \gamma_{yz}, \gamma_{xz}]$ con $\gamma_{ij} = 2\varepsilon_{ij}$ (ADR 0012). Sin variantes de hipótesis.

### 1. Lo que es igual a la 2D

- Ley de esfuerzo (a), $\boldsymbol\sigma = (1-d^-)\,\bar{\boldsymbol\sigma} - (d^+-d^-)\,\bar{\boldsymbol\sigma}^+$, comportamiento uniaxial y caveat de energía (2D, §1: con $\nu \ne 0$ no deriva de un potencial; los ciclos se midieron también en 3D).
- Daño de tracción con Rankine sobre el esfuerzo efectivo y banda de grieta (2D, §2); daño de compresión opcional con Drucker-Prager y la ley parabólica de Feenstra (2D, §3).
- Variables internas (2D, §4), tangente (2D, §7, ya escrita en forma tensorial), integración implícita e IMPL-EX (2D, §8), simetría y caveats (2D, §9–§10).

### 2. Lo que difiere

- **Esfuerzo** (2D, §6): $\bar{\boldsymbol\sigma} = \mathbf C_e:\boldsymbol\varepsilon$ con la elasticidad 3D y descomposición espectral 3×3 de $\bar{\boldsymbol\sigma}$; $\bar{\boldsymbol\sigma}^+ = \sum_i\langle\bar\sigma_i\rangle_+\,\mathbf n_i\otimes\mathbf n_i$ con los tres autovalores.
- **Tangente con daño fijo** (2D, §7): $\mathbf Q^+$ con los tres términos $H(\bar\sigma_i)\,\mathbf M_{ii}\otimes\mathbf M_{ii}$ y los tres pares $\tfrac12 c_{ij}\,\mathbf S_{ij}\otimes\mathbf S_{ij}$, con su límite en autovalores iguales (dobles y triples); en Voigt 3D, factor 2 en las tres columnas de cortante.
- **Rankine**: $\bar\sigma_I$, el mayor autovalor del esfuerzo efectivo 3D.
- **Longitud característica**: la siembran los sólidos 3D (`Hex8`, `Hex20`, `Hex27`: $\sqrt[3]{V}$; `Tet4`, `Tet10`: $\sqrt[3]{6V}$).

---

## Contrato de implementación

```yaml
name: TensionCompressionDamage3D
kind: material
status: draft

interface:
  strain_dim: 6
  primary_state_var: damage_t
  is_symmetric: false      # siempre, también con IMPL-EX (ley (a), 2D §9)
  state_schema: {kappa_t: [], kappa_c: [], damage_t: [], damage_c: [], h: [], kappa_t_prev: [], kappa_c_prev: [], damage_error: []}
  batch_kernel: true
  uses_pseudo_time: "true con integration='impl-ex' (ADR 0034)"

admission:                # la de la 2D (Q5: programa principal, usuario, 2026-11-18)
  textbooks:
    - de_souza_neto: "cap. 12, §12.6.1–12.6.2, cierre de grieta con partición espectral del esfuerzo: la familia, con otra forma (un solo D y hD; texto leído el 2026-11-18; 2D, Q5)"
  general_codes:
    - opensees: "ASDConcrete3D (la misma partición del esfuerzo efectivo; criterios de daño distintos)"
    - code_aster: "ENDO_ORTH_BETON; ENDO_ISOT_BETON"
    - abaqus: "Concrete Damaged Plasticity (un solo daño, sin partición espectral)"

parameters:
  - { name: E,          type: float, required: true,  desc: "Módulo de Young intacto [Pa]" }
  - { name: nu,         type: float, required: true,  desc: "Coeficiente de Poisson, en [0, 0.5)" }
  - { name: f_t,        type: float, required: true,  desc: "Resistencia a tracción [Pa]" }
  - { name: G_f,        type: float, required: true,  desc: "Energía de fractura en tracción [N/m]" }
  - { name: f_c,        type: float, required: false, default: null, desc: "Resistencia a compresión uniaxial [Pa]; sin ella d⁻ ≡ 0" }
  - { name: G_c,        type: float, required: false, default: null, desc: "Energía de aplastamiento [N/m]; obligatoria con f_c" }
  - { name: beta,       type: float, required: false, default: 1.16, desc: "f_bc/f_c" }
  - { name: integration, type: str,  required: false, default: "implicit", desc: "implicit | impl-ex" }
  - { name: damage_max, type: float, required: false, default: 0.999999, desc: "Tope de los dos daños" }
  - { name: density,    type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

state_variables:
  - { name: kappa_t,      type: float, shape: [], desc: "como la 2D" }
  - { name: kappa_c,      type: float, shape: [], desc: "como la 2D" }
  - { name: damage_t,     type: float, shape: [], desc: "como la 2D" }
  - { name: damage_c,     type: float, shape: [], desc: "como la 2D" }
  - { name: h,            type: float, shape: [], desc: "como la 2D" }
  - { name: kappa_t_prev, type: float, shape: [], desc: "como la 2D" }
  - { name: kappa_c_prev, type: float, shape: [], desc: "como la 2D" }
  - { name: damage_error, type: float, shape: [], desc: "como la 2D" }

conventions:
  voigt: "[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz], γ_ij = 2·ε_ij (ADR 0012)"
  units: "como la 2D"
  sign: "σ > 0 ⇔ tracción; f_c > 0"

validity:
  - "como la 2D: pequeñas deformaciones; E > 0, 0 ≤ ν < 0.5, f_t > 0, G_f > 0; con compresión, f_c > f_t y G_c > 0"
  - "0 < h < 2·E·G_f/f_t² en cada punto"
  - "sólo en sólidos 3D que proporcionan h"

out_of_scope:
  - "los de la 2D"

acceptance:
  material_point:
    - name: same_as_2d
      setup: "los casos de punto material de la 2D en 3D (los de esfuerzo plano, con el estado uniaxial 3D)"
      expect: "mismos resultados y tolerancias"
    - name: tangent_fd_3d
      setup: "daño fijo con autovalores distintos, dobles y triples del mismo signo; implícita en carga activa con autovalores distintos (doble o triple negativo con d⁻); se excluyen la arista de Rankine y los autovalores nulos (2D, §10)"
      expect: "|C_FD − C_tan| ≤ 1e-6·max|C_FD|; con ν = 0 y daño fijo, simétrica"
      tol_rel: 1.0e-6
    - name: batch_kernel_matches
      expect: "kernel por lotes bit a bit con compute_state"
  impl_ex:
    - name: same_as_2d
      expect: "el bloque impl_ex de la 2D en 3D"
  mesh_objectivity:
    - name: bar_in_tension_3d
      setup: "barra con una capa debilitada, Hex8, Tet4, Hex20, Hex27, Tet10"
      expect: "D/(G_f·A) = 1 ± 1 %"
      tol_rel: 1.0e-2
    - name: bar_in_compression_3d
      setup: "la barra de compresión de la 2D (L = 0.2 m, h ≤ 40 mm) con los sólidos 3D"
      expect: "trabajo sobre los elementos de la banda desde su pico = G_c·A ± 2 % e igual entre mallas"
      tol_rel: 2.0e-2

references:
  - "Las de la 2D."
```

---

## Implementación

*Pendiente de la aprobación de la spec (con la 2D).*

---

## Diálogo

- **2026-11-11** · Spec creada por la IA junto con [TensionCompressionDamage2D](TensionCompressionDamage2D.md); las preguntas al usuario (Q1–Q6) están en esa spec y valen para las dos.
- **2026-11-17** · Q1 decidida por el usuario: ley (a) de Faria, Oliver y Cervera. §1–§2 y el contrato reescritos para (a) junto con la 2D (la versión anterior describía la forma complementaria (d), descartada). Después de la revisión de `teorico-solidum` (H13): `state_variables`, `validity` y `out_of_scope` explícitos, `bar_in_compression_3d`, bloque `impl_ex` y `tangent_fd_3d` sin los puntos no diferenciables.
- **2026-11-18** · Q5 decidida por el usuario: programa principal, como la 2D (Diálogo de esa spec).
