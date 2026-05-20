# Elastic3D — Elasticidad lineal isótropa 3D

> Material elástico isótropo en 3D, espejo natural de `Elastic2D` con la convención Voigt 6D fijada en ADR 0012.

---

## Especificación física

### 0. Descripción general

Modelo elástico lineal isótropo tridimensional, parametrizado por $(E, \nu)$. A diferencia de `Elastic2D`, **no tiene variantes de hipótesis cinemática**: en 3D no aplican `plane_stress` ni `plane_strain` — toda la deformación es activa. Independiente de la velocidad, sin variables internas. Es el material elástico canónico para todos los sólidos 3D del catálogo.

### 1. Descomposición de la deformación

No aplica — toda la deformación es elástica:

$$\boldsymbol\varepsilon = \boldsymbol\varepsilon^e$$

### 2. Ley elástica

$$\boldsymbol\sigma = \mathbf C_e\,\boldsymbol\varepsilon$$

con $\boldsymbol\varepsilon$ y $\boldsymbol\sigma$ en **notación Voigt 6D del proyecto** (ADR 0012, `Reglas.md §5`):

$$\boldsymbol\varepsilon = [\varepsilon_{xx},\ \varepsilon_{yy},\ \varepsilon_{zz},\ \gamma_{xy},\ \gamma_{yz},\ \gamma_{xz}]^\top, \qquad \gamma_{ij} = 2\,\varepsilon_{ij}$$

$$\boldsymbol\sigma = [\sigma_{xx},\ \sigma_{yy},\ \sigma_{zz},\ \sigma_{xy},\ \sigma_{yz},\ \sigma_{xz}]^\top$$

**Matriz constitutiva elástica 6×6**:

$$\mathbf C_e = \frac{E}{(1+\nu)(1-2\nu)}\begin{bmatrix}
1-\nu & \nu & \nu & 0 & 0 & 0 \\
\nu & 1-\nu & \nu & 0 & 0 & 0 \\
\nu & \nu & 1-\nu & 0 & 0 & 0 \\
0 & 0 & 0 & \tfrac{1-2\nu}{2} & 0 & 0 \\
0 & 0 & 0 & 0 & \tfrac{1-2\nu}{2} & 0 \\
0 & 0 & 0 & 0 & 0 & \tfrac{1-2\nu}{2}
\end{bmatrix}$$

Tangente constante $\mathbf C_e$ en todo régimen. Simétrica y positiva definida en el rango admisible de $\nu$.

El factor $\tfrac{1-2\nu}{2}$ en los cortantes (no $1-2\nu$) es consecuencia del uso de $\gamma_{ij}$ *engineering*: $\sigma_{ij} = G\,\gamma_{ij}$ con $G = \tfrac{E}{2(1+\nu)} = \tfrac{E}{(1+\nu)(1-2\nu)}\cdot\tfrac{1-2\nu}{2}$.

### 3-6. Superficie de fluencia, flujo, endurecimiento, Kuhn-Tucker

No aplican — el material es indefinidamente elástico.

### 7. Variables internas

Ninguna — `state_vars = None` se acepta y se devuelve intacto. Sin `PRIMARY_STATE_VAR`.

### 8. Notación Voigt

Orden y factores fijados por **ADR 0012**:

- Componentes diagonales primero: `[xx, yy, zz]`.
- Bloque cortante después: `[xy, yz, xz]` (extensión natural del orden 2D del proyecto `[xx, yy, xy]`).
- Deformaciones angulares *engineering*: $\gamma_{ij} = 2\,\varepsilon_{ij}$.
- Esfuerzos tensoriales sin factor.

**Compatibilidad con códigos externos**: ABAQUS usa `[11, 22, 33, 12, 13, 23]`. La diferencia con la convención del proyecto es la permutación del bloque cortante (`yz ↔ xz`). Si se importan datos desde ABAQUS, aplicar permutación 5↔6 una vez en el preprocesador.

---

## Formulación numérica

### 9. Esquema temporal

No aplica — sin historia.

### 10. Return mapping / actualización

No aplica. Evaluación reducida a producto matriz-vector:

$$\boldsymbol\sigma_{n+1} = \mathbf C_e\,\boldsymbol\varepsilon_{n+1}, \qquad \mathbf C^{\text{alg}} = \mathbf C_e$$

La matriz $\mathbf C_e$ se precalcula y guarda en el constructor — coste de `compute_state` es un producto 6×6 · 6 sin alocación adicional.

### 11. Tangente algorítmica consistente

$$\mathbf C^{\text{alg}} = \mathbf C_e$$

Simétrica (`IS_SYMMETRIC = True`). Positiva definida para $\nu \in (-1, 0.5)$ estricto; diverge en $\nu \to 0.5$ (incompresible — requiere formulación mixta u-p, no soportada).

### 12. Caveats numéricos

- **Incompresibilidad $\nu \to 0.5$**: la matriz constitutiva diverge ($1-2\nu \to 0$). El constructor rechaza $\nu \ge 0.5$ con `ValueError`. Recomendado $\nu \le 0.49$ para metales y polímeros casi-incompresibles.
- **Locking volumétrico** en `Hex8`/`Tet4` con $\nu \to 0.5$: declarado limitación arquitectural en la spec de cada elemento, sin mitigación implementada (B-bar/F-bar diferidas a caso de uso real, política idéntica al 2D).
- Para $\nu < 0$ el material es admisible termodinámicamente (auxético) pero infrecuente en la práctica; el constructor lo permite hasta $\nu > -1$.

---

## Contrato de implementación

```yaml
name: Elastic3D
kind: material
status: validated

interface:
  strain_dim: 6
  primary_state_var: null
  is_symmetric: true

parameters:
  - { name: E,       type: float, required: true,  desc: "Módulo de Young (>0)" }
  - { name: nu,      type: float, required: true,  desc: "Coeficiente de Poisson ∈ (-1, 0.5)" }
  - { name: density, type: float, required: false, default: null,
      desc: "Densidad (kg/m³). Opcional al construir; obligatoria solo si se ensambla peso propio o masa (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(6,), state_vars=None) -> (σ: ndarray(6,), C_e: ndarray(6,6), state_vars)"
  strain_kind: "3D Voigt — [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

state_variables: []

conventions:
  sign: "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "[xx, yy, zz, xy, yz, xz] (ADR 0012); γ_ij = 2·ε_ij en deformación; σ_ij tensorial sin factor"
  hypothesis_dispatch: "no aplica en 3D"

validity:
  - "E > 0"
  - "ν ∈ (-1, 0.5) estricto; ν = 0.5 rechazado (incompresible)"
  - "pequeñas deformaciones (|ε| ≲ 1e-2)"

out_of_scope:
  - "anisotropía (ortótropo, monoclínico, totalmente anisótropo)"
  - "plasticidad ⇒ VonMises3D, DruckerPrager3D (sub-etapa posterior)"
  - "daño ⇒ IsotropicDamage3D (sub-etapa posterior)"
  - "incompresibilidad estricta (ν=0.5) ⇒ requiere formulación mixta u-p"
  - "grandes deformaciones"

acceptance:
  verification:
    - name: simetria_C_e
      setup: "construir Elastic3D y comparar C_e con su traspuesta"
      expect: "C_e == C_e.T exacto componente a componente"
      tol_abs: 1.0e-14
    - name: positividad_definida
      setup: "autovalores de C_e con ν ∈ {0, 0.3, 0.49}"
      expect: "todos los autovalores estrictamente positivos"
      tol_abs: 0.0
    - name: traccion_uniaxial
      setup: "ε = (ε_0, -ν·ε_0, -ν·ε_0, 0, 0, 0)"
      expect: "σ = (E·ε_0, 0, 0, 0, 0, 0) exacto"
      tol_rel: 1.0e-12
    - name: cortante_puro_xy
      setup: "ε = (0, 0, 0, γ_0, 0, 0) con γ_0 = 1e-3"
      expect: "σ = (0, 0, 0, G·γ_0, 0, 0) con G = E/[2(1+ν)]"
      tol_rel: 1.0e-12
    - name: cortante_puro_yz
      setup: "ε = (0, 0, 0, 0, γ_0, 0)"
      expect: "σ = (0, 0, 0, 0, G·γ_0, 0)"
      tol_rel: 1.0e-12
    - name: cortante_puro_xz
      setup: "ε = (0, 0, 0, 0, 0, γ_0)"
      expect: "σ = (0, 0, 0, 0, 0, G·γ_0)"
      tol_rel: 1.0e-12
    - name: compresion_hidrostatica
      setup: "ε = (ε_0, ε_0, ε_0, 0, 0, 0) con ε_0 = -1e-3 (compresión isótropa)"
      expect: "σ = (K·3ε_0, K·3ε_0, K·3ε_0, 0, 0, 0) con K = E/[3(1-2ν)] (módulo de compresibilidad)"
      tol_rel: 1.0e-12
    - name: rechazo_inputs_invalidos
      setup: "construir Elastic3D con E≤0, ν fuera de rango (-1, 0.5), density<0"
      expect: "ValueError con mensaje claro en cada caso"

references:
  - "Timoshenko S., Goodier J.N. (1970). Theory of Elasticity. McGraw-Hill. §6 (3D elasticity)."
  - "Bathe K.J. (2014). Finite Element Procedures. Prentice Hall. §6.6."
  - "Cook R.D., Malkus D.S., Plesha M.E., Witt R.J. (2002). Concepts and Applications of FEA. Wiley. §6.8."
  - "ADR 0012 — Sólidos 3D: convención Voigt 6D y cierre del contrato internal_forces."
```

---

## Implementación

- **Archivo**: [solidum/materials/elastic_3d.py](../../solidum/materials/elastic_3d.py).
- **Clase**: `Elastic3D`, registrada vía `@MaterialRegistry.register`.
- **Despacho**: matriz `C` precalculada en el constructor. `compute_state` reduce a `self.C @ strain` sin ramificaciones — coste mínimo.
- **Validación temprana** en constructor: `E > 0`, `ν ∈ (-1, 0.5)`, `density ≥ 0` si se pasa. Mensajes en español.
- **Tests**: [tests/test_solid_3d.py](../../tests/test_solid_3d.py) — clase `TestElastic3D` (11 tests: simetría C, positividad definida, tracción uniaxial, los tres cortantes puros, compresión hidrostática con K = E/[3(1-2ν)], rechazos de inputs inválidos, density opcional).
- **Validación de integración**: `tests/validation/test_cube_lame_3d.py` consume Elastic3D end-to-end vía Hex8 y Tet4 (tracción uniaxial + compresión hidrostática con solución analítica exacta).

---

## Diálogo

- **2026-05-19** · Spec creada como pieza inicial de la Etapa 7 (sólidos 3D, alcance acotado Hex8 + Tet4 + Elastic3D). Convención Voigt 6D fijada en ADR 0012. Sin variantes de hipótesis (no aplica plane_stress/plane_strain en 3D).
