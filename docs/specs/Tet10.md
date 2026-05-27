# Tet10 — tetraedro cuadrático isoparamétrico (10 nodos)

> Elemento sólido 3D de orden cuadrático sobre el simplex tetraédrico. Análogo 3D del ``Tri6`` 2D. Reproduce campos cuadráticos completos en (ξ, η, ζ). Mitiga drásticamente el shear locking y el locking volumétrico del ``Tet4`` (CST 3D). Sub-fase 3 de A.ter; subclase de ``_HigherOrderSolid3D``.

---

## Especificación física

### 0. Descripción general

Tetraedro tridimensional de segundo orden con **10 nodos**: 4 vértices + 6 medios de arista. Aproximación cuadrática $C^0$ en coordenadas baricéntricas $L_1, L_2, L_3, L_4$ con $\sum_i L_i = 1$. Formulación isoparamétrica. Régimen de pequeños desplazamientos y deformaciones.

### 1-5. Cinemática, ecuación constitutiva, equilibrio

Idénticas al ``Hex20``/``Hex27`` — la diferencia con los hexaedros está exclusivamente en la geometría de referencia (simplex tetraédrico vs cubo).

---

## Formulación numérica (FEM)

### 6. Discretización — convención VTK_QUADRATIC_TETRA

Los 10 nodos sobre el tetraedro de referencia con vértices $(0,0,0)$, $(1,0,0)$, $(0,1,0)$, $(0,0,1)$:

| Nodo  | Tipo            | Coords naturales $(\xi, \eta, \zeta)$ |
|-------|-----------------|---------------------------------------|
| 0     | vértice (idéntico Tet4) | $(0, 0, 0)$                  |
| 1     | vértice (idéntico Tet4) | $(1, 0, 0)$                  |
| 2     | vértice (idéntico Tet4) | $(0, 1, 0)$                  |
| 3     | vértice (idéntico Tet4) | $(0, 0, 1)$                  |
| 4     | medio arista 0-1 | $(0.5, 0, 0)$                       |
| 5     | medio arista 1-2 | $(0.5, 0.5, 0)$                     |
| 6     | medio arista 2-0 | $(0, 0.5, 0)$                       |
| 7     | medio arista 0-3 | $(0, 0, 0.5)$                       |
| 8     | medio arista 1-3 | $(0.5, 0, 0.5)$                     |
| 9     | medio arista 2-3 | $(0, 0.5, 0.5)$                     |

Mapeo isoparamétrico desde el tetraedro de referencia con $L_1 = 1 - \xi - \eta - \zeta$, $L_2 = \xi$, $L_3 = \eta$, $L_4 = \zeta$.

### 7. Funciones de forma

**Vértices** ($i \in \{1, 2, 3, 4\}$ — corresponde a los nodos locales $\{0, 1, 2, 3\}$):

$$N_i = L_i\,(2\,L_i - 1).$$

**Medios de arista** entre vértices $i, j$:

$$N_{ij} = 4\,L_i\,L_j.$$

Las seis aristas se etiquetan en el orden VTK_QUADRATIC_TETRA: $(1,2), (2,3), (3,1), (1,4), (2,4), (3,4)$ — nodos locales $(0,1), (1,2), (2,0), (0,3), (1,3), (2,3)$. Garantizan partición de la unidad $\sum N_i = 1$ y propiedad delta de Kronecker.

**Espacio reproducido**: el ``Tet10`` reproduce **completamente** el espacio cuadrático en barycentric coords — todos los polinomios de grado total $\le 2$ en $(\xi, \eta, \zeta)$. 10 monomios, idéntico al número de nodos.

### 8. Matriz $\mathbf B$

Tamaño **6×30**, ensamblada con la misma plantilla 6×3 del ``Hex8``/``Hex20``/``Hex27`` por nodo. Derivadas en coords naturales obtenidas con regla del producto: para vértices $\partial N_i / \partial x_k = (4 L_i - 1)\,\partial L_i / \partial x_k$; para medios $\partial N_{ij}/\partial x_k = 4(L_j\,\partial L_i/\partial x_k + L_i\,\partial L_j/\partial x_k)$.

### 9-10. Rigidez y fuerzas internas

Implementadas en la base ``_HigherOrderSolid3D``. Tamaño 30×30 para $\mathbf K_e$.

### 11. Cuadratura

Por defecto **Stroud 4 puntos** (``tet_4``): orden de exactitud 2. Para Tet10 con Jacobiano constante el integrando de $\mathbf B^\top \mathbf D\, \mathbf B$ es a lo sumo de grado 2 (B es lineal en barycentric), así que la regla es **exacta** para K.

Para problemas con Jacobiano variable (geometría distorsionada, mid-edges fuera del centro real) el integrando puede subir de grado; en esos casos el usuario puede seleccionar `tet_15` (Keast 15 puntos, orden 5) vía el parámetro `quadrature`.

**Cuadratura de masa**: el ``Tet10`` declara ``_MASS_QUADRATURE = "tet_15"`` (Keast 15 puntos, orden 5) **independientemente** de la cuadratura de K. El integrando $\rho\,\mathbf N^\top \mathbf N$ es de grado 4 (cuadrático × cuadrático) y la regla de 4 puntos lo subintegra; ``tet_15`` lo integra exactamente. Esto garantiza que la masa total sea correcta y que análisis modal/transitorio con `Tet10` no sufran errores de cuadratura.

### 12. Cargas distribuidas consistentes

**Body load**: ``∫ N^T b dV`` con la cuadratura del elemento. Reparto no uniforme entre vértices y medios; suma global ``b·V_e`` exacta para ``b`` uniforme y geometría sin distorsión severa.

**Face traction**: cada cara del ``Tet10`` tiene **6 nodos** (3 vértices + 3 medios de arista) y se aproxima por funciones cuadráticas Tri6 sobre el triángulo de referencia $(0,0)$-$(1,0)$-$(0,1)$. Cuadratura ``tri_3`` (3 puntos, orden 2) — exacta para tracción uniforme sobre cara plana.

| Cara | Vértices Tet10 | Medios Tet10 (en orden Tri6) | Normal saliente   |
|------|----------------|-------------------------------|-------------------|
| 0    | (1, 2, 3)      | (5, 9, 8)                     | opuesta al nodo 0 |
| 1    | (0, 3, 2)      | (7, 9, 6)                     | opuesta al nodo 1 |
| 2    | (0, 1, 3)      | (4, 8, 7)                     | opuesta al nodo 2 |
| 3    | (0, 2, 1)      | (6, 5, 4)                     | opuesta al nodo 3 |

El orden dentro de cada cara reproduce la convención de ``_N_tri6`` del paquete 2D: vértices `(0, 1, 2)` antihorarios, luego medios `(0-1, 1-2, 2-0)`.

### 13. Salida por punto de Gauss

`compute_gauss_state(U)` con 4 puntos por defecto. Sin `internal_forces` (ADR 0012).

---

## Limitaciones declaradas

### Locking volumétrico atenuado pero presente

`Tet10` con cuadratura ``tet_4`` y material casi-incompresible (`Elastic3D` con $\nu \to 0.5$, o `VonMises3D` plástico desarrollado) sufre **locking volumétrico** mucho menor que el `Tet4` (cuyo shear locking severo es la motivación principal del `Tet10`), pero todavía presente. Sin mitigación implementada. Blindado por test (mismo patrón que `Hex20`/`Hex27`).

### Distorsión severa

Para Tet10 con mid-edges desplazados fuera del centro real de la arista (mid-edges curvos para representar geometría curva), el Jacobiano deja de ser constante y la cuadratura ``tet_4`` puede subintegrar K. Recomendación: usar ``tet_15`` en mallas con curvatura severa.

---

## Contrato de implementación

```yaml
name: Tet10
kind: element
status: validated

interface:
  dof_names: [ux, uy, uz]
  n_nodes: 10
  strain_dim: 6
  n_integration_points: 4       # default Stroud tet_4

parameters:
  - { name: quadrature, type: str, required: false, default: "tet_4",
      desc: "Cuadratura tetraédrica desde QuadratureRegistry. Default
             tet_4 (Stroud 4 puntos, orden 2, suficiente para K con J
             constante). Alternativa tet_15 (Keast 15 puntos, orden 5)
             para geometrías distorsionadas." }

material_contract:
  signature: "compute_state(ε, state) -> (σ, C_tangent, state')"
  strain_kind: "Voigt 3D [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

conventions:
  sign: "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "[xx, yy, zz, xy, yz, xz] con γ_ij = 2·ε_ij (ADR 0012)"
  node_orientation: "VTK_QUADRATIC_TETRA; 0-3 vértices (Tet4 con volumen
                     positivo), 4-9 medios de arista en orden (0-1, 1-2,
                     2-0, 0-3, 1-3, 2-3)"
  face_numbering: "ADR 0012 — 4 caras triangulares con normal saliente,
                   paritarias con Tet4 en los vértices; cada cara añade
                   3 medios de arista en orden Tri6"

validity:
  - "pequeños desplazamientos y deformaciones"
  - "geometría sin distorsión severa (det J > tol en todos los Gauss)"
  - "compatible con materiales con STRAIN_DIM = 6"

out_of_scope:
  - "grandes deformaciones, formulaciones corotacionales"
  - "internal_forces (ADR 0012)"
  - "mid-edges severamente curvos (usar tet_15 si es necesario)"

acceptance:
  verification:
    - name: patch_test_lineal
      setup: "campo lineal u_i = a_ij·x_j en los 10 nodos del contorno"
      expect: "ε constante = sym(a) en todos los Gauss"
      tol_rel: 1.0e-12
    - name: patch_test_cuadratico
      setup: "campo u_x = c·x² en los 10 nodos"
      expect: "ε_xx(x) = 2c·x exacto en todos los Gauss"
      tol_rel: 1.0e-12
    - name: modos_de_cuerpo_rigido
      setup: "K_e proyectado sobre 6 modos rígidos 3D"
      expect: "exactamente 6 autovalores nulos"
      tol_abs: 1.0e-9
    - name: simetria_K_e
      setup: "K_e con Elastic3D"
      expect: "K_e == K_e.T exacto"
      tol_abs: 1.0e-12
  specific:
    - name: cubo_lame_3d_uniaxial
      setup: "malla 5 tetraedros del cubo unitario, refinados a Tet10,
              con tracción uniaxial; el campo analítico es lineal y
              Tet10 lo reproduce exacto"
      expect: "u_x = p·x/E etc. exacto a precisión máquina"
      tol_rel: 1.0e-10
    - name: macneal_3d_tet10_dominates_tet4
      setup: "cantilever de Tet10 con la misma malla que Tet4 sufre
              dramáticamente menos shear locking"
      expect: "Tet10 alcanza > 80% u_EB con ~12 tetraedros (Tet4 << 20%)"
    - name: body_load_uniforme
      setup: "Tet10 ref con b = (0, 0, -g)"
      expect: "Σf_z = -g·V_e = -g/6"
      tol_rel: 1.0e-12
    - name: face_traction_balance
      setup: "tracción uniforme sobre una cara del Tet10"
      expect: "Σf = A_face · t̄ exacto; nodos fuera de la cara reciben 0"
      tol_rel: 1.0e-12
    - name: mass_consistent_total_with_tet15
      setup: "Tet10 cubo ref con ρ=7850, masa consistente con tet_15"
      expect: "Σ M = 3·ρ·V_e = 3·ρ/6"
      tol_rel: 1.0e-9
    - name: mass_lumped_HRZ_positivo
      setup: "masa lumped HRZ del Tet10"
      expect: "diagonal estricta, todas las entradas positivas,
              conservación de masa total"
      tol_abs: 1.0e-12

references:
  - "Bathe K.J. (2014). Finite Element Procedures. §5.3.3 (tetraedro 3D)."
  - "Cook R.D., Malkus D.S., Plesha M.E., Witt R.J. (2002). Concepts and
    Applications of FEA. §6.5."
  - "Keast P. (1986). Moderate-degree tetrahedral quadrature formulas.
    CMAME 55(3), 339-348 (regla tet_15 de masa)."
  - "Stroud A.H. (1971). Approximate calculation of multiple integrals
    (regla tet_4 de K)."
  - "ADR 0012 — Sólidos 3D: convención Voigt 6D y caras."
```

---

## Implementación

- **Archivo**: [solidum/elements/solid_3d/tet10.py](../../solidum/elements/solid_3d/tet10.py).
- **Clase**: `Tet10` (registrada vía `@ElementRegistry.register`); subclase de [`_HigherOrderSolid3D`](../../solidum/elements/solid_3d/_shared.py).
- **Funciones núcleo**: `_N_tet10(xi, eta, zeta)`, `_dN_tet10(xi, eta, zeta)` en [solidum/elements/solid_3d/_shared.py](../../solidum/elements/solid_3d/_shared.py). Implementadas en numpy puro vía coordenadas baricéntricas $L_i$.
- **Cuadraturas registradas**: `tet_4` (Stroud, orden 2) y `tet_15` (Keast, orden 5) en [solidum/math/integration.py](../../solidum/math/integration.py). El `tet_4` cubre la integración de K; el `tet_15` se usa exclusivamente para la masa consistente (declarada en `Tet10._MASS_QUADRATURE`).
- **Tests**:
  - `tests/test_solid_3d_higher_order.py` — `TestTet10Element` (patch lineal + cuadrático, simetría K, body load, face traction, masa consistente con tet_15 y lumped HRZ, gauss state).
  - `tests/test_rigid_body_modes.py` — `TestRBMTet10` (3 trans + 3 rot + rank-deficiency = 6).
  - `tests/validation/test_cube_lame_3d.py` — extensión con malla Tet10 del cubo unitario; uniaxial e hidrostática exactas a precisión máquina.

---

## Diálogo

- **2026-05-27** · Spec inicial. Sub-fase 3 de A.ter — cierra el catálogo cuadrático 3D con la rama tetraédrica. Subclase de `_HigherOrderSolid3D` igual que `Hex20`/`Hex27`, pero con cara triangular (`Tri6` shape functions importadas del paquete 2D, cuadratura `tri_3`) y cuadratura volumétrica tetraédrica nueva (`tet_4` + `tet_15`).
- **2026-05-27** · Implementado y validado. Patch test lineal y cuadrático exactos a precisión máquina, 6 modos rígidos, body load suma exacta `b·V_e`, face traction balance exacto sobre la cara triangular, masa consistente con `tet_15` exacta a precisión máquina (conservación `3·ρ·V`), lumped HRZ con todas las entradas positivas. Spec a `status: validated`.
