# Hex27 — hexaedro Lagrangiano 3D triquadrático (27 nodos)

> Elemento sólido 3D de orden cuadrático con interpolación **Lagrangiana completa** — producto tensorial de Lagrange cuadrático en las tres direcciones. Reproduce **todos** los polinomios triquadráticos ``ξ^a η^b ζ^c`` con ``a, b, c ∈ {0, 1, 2}``, incluyendo los puramente triquadráticos (``ξ²η²ζ²`` y combinaciones) que faltan en el serendípito ``Hex20``. Sub-fase 2 de A.ter; dispara la centralización en ``_HigherOrderSolid3D`` (regla de los dos casos reales).

---

## Especificación física

### 0. Descripción general

Hexaedro tridimensional de segundo orden con **27 nodos**: 8 vértices + 12 medios de arista + 6 centros de cara + 1 centro del cuerpo. Aproximación Lagrangiana $C^0$ completa en coordenadas naturales $(\xi, \eta, \zeta) \in [-1, 1]^3$. Formulación isoparamétrica. Régimen de pequeños desplazamientos y deformaciones.

### 1-5. Cinemática, ecuación constitutiva, equilibrio

Idénticas al ``Hex20`` — la diferencia entre serendipity y Lagrange está exclusivamente en el espacio de aproximación discreto, no en la formulación variacional.

---

## Formulación numérica (FEM)

### 6. Discretización — convención VTK_TRIQUADRATIC_HEXAHEDRON

Los 27 nodos en el cubo de referencia $[-1, 1]^3$:

| Rango  | Tipo            | Posiciones $(\xi_i, \eta_i, \zeta_i)$ |
|--------|-----------------|---------------------------------------|
| 0-7    | vértices        | idénticos al `Hex8`/`Hex20`           |
| 8-19   | medios de arista | idénticos al `Hex20`                  |
| 20     | centro de cara  | $(-1, 0, 0)$ — cara 5 (−ξ)            |
| 21     | centro de cara  | $(+1, 0, 0)$ — cara 3 (+ξ)            |
| 22     | centro de cara  | $(0, -1, 0)$ — cara 2 (−η)            |
| 23     | centro de cara  | $(0, +1, 0)$ — cara 4 (+η)            |
| 24     | centro de cara  | $(0, 0, -1)$ — cara 0 (−ζ)            |
| 25     | centro de cara  | $(0, 0, +1)$ — cara 1 (+ζ)            |
| 26     | centro del cuerpo | $(0, 0, 0)$                         |

### 7. Funciones de forma — producto tensorial de Lagrange cuadrático

Para cada nodo se asigna una terna de índices $(i, j, k) \in \{-1, 0, +1\}^3$ correspondiente a su posición natural. Las funciones de forma son productos tensoriales:

$$N_n(\xi, \eta, \zeta) = L_{i_n}(\xi)\,L_{j_n}(\eta)\,L_{k_n}(\zeta),$$

con los tres polinomios de Lagrange cuadrático 1D:

$$L_{-1}(x) = \tfrac{1}{2}\,x(x - 1),\qquad L_0(x) = 1 - x^2,\qquad L_{+1}(x) = \tfrac{1}{2}\,x(x + 1).$$

Cumplen $L_{-1}(-1) = L_0(0) = L_{+1}(+1) = 1$ y cero en los demás. Garantiza propiedad delta de Kronecker $N_m(\xi_n, \eta_n, \zeta_n) = \delta_{mn}$ y partición de la unidad $\sum_n N_n = 1$.

**Espacio reproducido**: el ``Hex27`` reproduce **completamente** el espacio triquadrático $Q_2 = \mathrm{span}\{\xi^a \eta^b \zeta^c : a, b, c \in \{0, 1, 2\}\}$ — 27 monomios, idéntico al número de nodos. La diferencia con el ``Hex20`` serendípito está en los términos puramente triquadráticos ``ξ²η²``, ``η²ζ²``, ``ξ²ζ²``, ``ξ²η²ζ``, ``ξ²ηζ²``, ``ξη²ζ²``, ``ξ²η²ζ²`` que ``Hex20`` no captura.

### 8. Matriz $\mathbf B$

Tamaño **6×81**, ensamblada con la misma plantilla 6×3 del ``Hex8``/``Hex20`` por nodo. Derivadas obtenidas vía regla del producto sobre los polinomios de Lagrange:

$$\frac{\partial N_n}{\partial \xi} = L'_{i_n}(\xi)\,L_{j_n}(\eta)\,L_{k_n}(\zeta),$$

con

$$L'_{-1}(x) = x - \tfrac{1}{2},\qquad L'_0(x) = -2x,\qquad L'_{+1}(x) = x + \tfrac{1}{2}.$$

### 9-10. Rigidez y fuerzas internas

Idénticas a ``Hex8``/``Hex20`` en estructura. El tamaño es 81×81.

### 11. Cuadratura

Por defecto **Gauss–Legendre 3×3×3** (27 puntos): orden de exactitud 5 en cada dirección. Suficiente para integrar exactamente $\mathbf B^\top \mathbf D\, \mathbf B$ del ``Hex27`` sobre Jacobiano constante (el integrando es a lo sumo de grado 3 en cada dirección — derivada cuadrático×cuadrático) y la masa consistente $\rho \mathbf N^\top \mathbf N$ (grado 4 en cada dirección).

**Aviso sobre integración reducida**: el esquema $2 \times 2 \times 2$ (`hex_2x2x2`, 8 puntos) está disponible pero **introduce muchos modos espurios** sobre un elemento aislado: verificación numérica directa del proyecto cuenta **27 modos de hourglass** (= 81 DOFs − 8 Gauss × 6 strain − 6 rigid; jump espectral limpio entre eigs[32] ≈ 1e-10 y eigs[33] ≈ 7e3). Es muchísimo peor que en el `Hex20` reducido (6 hourglass) — el espacio extra del Lagrangiano completo es subintegrado severamente por 2×2×2. La cita clásica de "1 modo bubble" (Cook-Malkus-Plesha-Witt §11.6) corresponde al **único modo que persiste tras ensamblar** en mallas estructuradas; en un elemento aislado el conteo real es el de arriba. Sin estabilización implementada. Recomendación operativa estricta: usar el default $3 \times 3 \times 3$ con `Hex27`.

### 12. Cargas distribuidas consistentes

**Fuerza de cuerpo**: integración con cuadratura del elemento, reparto no uniforme entre vértices, medios, centros de cara y centro del cuerpo. Suma global $\mathbf b\cdot V_e$ exacta para $\mathbf b$ uniforme.

**Tracción de superficie**: cada cara del ``Hex27`` tiene **9 nodos** (4 vértices + 4 medios de arista + 1 centro de cara) que se aproximan por funciones Lagrangianas Quad9 sobre el parámetro 2D $(s, t) \in [-1, 1]^2$.

| Cara | Vértices  | Medios de arista | Centro de cara | Normal saliente |
|------|-----------|------------------|----------------|-----------------|
| 0    | (0,3,2,1) | (11, 10, 9, 8)   | 24             | −ζ              |
| 1    | (4,5,6,7) | (12, 13, 14, 15) | 25             | +ζ              |
| 2    | (0,1,5,4) | (8, 17, 12, 16)  | 22             | −η              |
| 3    | (1,2,6,5) | (9, 18, 13, 17)  | 21             | +ξ              |
| 4    | (2,3,7,6) | (10, 19, 14, 18) | 23             | +η              |
| 5    | (3,0,4,7) | (11, 16, 15, 19) | 20             | −ξ              |

El orden dentro de cada cara reproduce la convención de ``_N_quad9`` del paquete 2D: primero los 4 vértices, después los 4 medios en orden cíclico, y finalmente el nodo central. Cuadratura 2D sobre la cara: Gauss $3 \times 3$ (exacta para tracción constante sobre cara plana).

### 13. Salida por punto de Gauss

`compute_gauss_state(U)` devuelve dict con $\boldsymbol\varepsilon$ y $\boldsymbol\sigma$ en cada uno de los 27 puntos de Gauss del elemento, junto con coordenadas naturales y globales. Sin `internal_forces` (ADR 0012).

---

## Limitaciones declaradas

### Coste computacional vs ``Hex20``

`Hex27` tiene **81 DOFs por elemento** vs 60 del `Hex20` (35% más). El espacio extra que captura (términos puramente triquadráticos) **rara vez gobierna la solución** en problemas físicos estándar — el `Hex20` serendípito es la elección por defecto en la mayoría de aplicaciones (Bathe FEP §5.3.2). El `Hex27` se prefiere cuando:

- La geometría tiene curvatura severa y se mapea con elementos isoparamétricos cuadráticos: el nodo central de cuerpo ayuda a representar la curvatura interior.
- La sensibilidad al locking es crítica y se prefiere intentar el Q_2 completo antes que recurrir a B-bar/F-bar.
- Se busca convergencia óptima $O(h^{p+1})$ con $p = 2$ en problemas donde el contenido espectral cubre el espacio completo.

Para problemas con flexión simple y geometría rectilínea, `Hex20` y `Hex27` dan resultados prácticamente idénticos a un coste menor para el primero.

### Locking volumétrico con $\nu \to 0.5$

Comparable al ``Hex20``: el espacio extra triquadrático añade flexibilidad sin eliminar la incompresibilidad — formulación mixta u-p o B-bar/F-bar son necesarias para eliminarlo. Declarado como limitación, blindado por test cuantitativo.

### Hourglass con integración reducida

`Hex27` integrado con `hex_2x2x2` tiene **27 modos de hourglass espurios por elemento aislado** (vs 6 del `Hex20` reducido y 12 del `Hex8` reducido). El espacio Lagrangiano completo es subintegrado severamente por 8 puntos. En mallas estructuradas la mayoría no propagan tras ensamblar, pero el conteo a nivel de elemento es alto. Sin estabilización.

---

## Contrato de implementación

```yaml
name: Hex27
kind: element
status: validated

interface:
  dof_names: [ux, uy, uz]
  n_nodes: 27
  strain_dim: 6
  n_integration_points: 27      # default Gauss 3×3×3

parameters:
  - { name: quadrature, type: str, required: false, default: "hex_3x3x3",
      desc: "Regla de cuadratura desde QuadratureRegistry; hex_2x2x2
             reducida con 1 modo de hourglass espurio (bubble) por
             elemento aislado, sin estabilización" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, C_tangent, state')"
  strain_kind: "Voigt 3D [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

conventions:
  sign: "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "[xx, yy, zz, xy, yz, xz] con γ_ij = 2·ε_ij (ADR 0012)"
  node_orientation: "VTK_TRIQUADRATIC_HEXAHEDRON; 0-7 vértices (Hex8), 8-19
                     medios de arista (Hex20), 20-25 centros de cara en orden
                     (-x, +x, -y, +y, -z, +z), 26 centro del cuerpo"
  face_numbering: "ADR 0012 — 6 caras con normal saliente, paritarias con
                   Hex8/Hex20 en los vértices y medios; cada cara añade el
                   nodo de centro correspondiente (índice 20-25)"

validity:
  - "pequeños desplazamientos y deformaciones"
  - "geometría sin distorsión severa (det J > tol en todos los puntos de Gauss)"
  - "compatible con materiales con STRAIN_DIM = 6 (Elastic3D, VonMises3D,
     DruckerPrager3D, IsotropicDamage3D)"

out_of_scope:
  - "grandes deformaciones, formulaciones corotacionales"
  - "estabilización contra hourglass para integración reducida"
  - "tracción variable sobre cara"
  - "internal_forces (ADR 0012)"
  - "Tet10 cuadrático (sub-fase 3)"

acceptance:
  verification:
    - name: patch_test_lineal
      setup: "campo lineal u_i = a_ij·x_j en los 27 nodos del contorno"
      expect: "ε constante = sym(a) en todos los Gauss; σ uniforme"
      tol_rel: 1.0e-12
    - name: patch_test_triquadratico
      setup: "campo u_x = c·x²·y²·z² impuesto en los 27 nodos"
      expect: "ε_xx(x,y,z) = 2c·x·y²·z² exacto en todos los Gauss (el espacio
              Q_2 lagrangiano completo lo contiene; Hex20 serendípito NO)"
      tol_rel: 1.0e-12
    - name: modos_de_cuerpo_rigido
      setup: "evaluar K_e y proyectar sobre 6 modos rígidos 3D"
      expect: "K_e · u_rigid ≈ 0; exactamente 6 autovalores nulos"
      tol_abs: 1.0e-9
    - name: simetria_K_e
      setup: "K_e con Elastic3D"
      expect: "K_e == K_e.T exacto"
      tol_abs: 1.0e-12
  specific:
    - name: cubo_lame_3d
      setup: "1 Hex27 en cubo unitario con tracción uniaxial y BCs mínimas"
      expect: "u_x = p·x/E etc. exacto a precisión máquina (campo lineal)"
      tol_rel: 1.0e-10
    - name: body_load_uniforme
      setup: "Hex27 cubo unitario con b = (0, 0, -g)"
      expect: "Σf_z = -g·V_e = -g"
      tol_rel: 1.0e-12
    - name: face_traction_balance
      setup: "tracción (p, 0, 0) sobre cara +ξ del cubo unitario"
      expect: "Σf_x = p·A_cara = p; reparto consistente Lagrange Q_2 (vértices,
              medios y centro de cara con sus pesos respectivos)"
      tol_rel: 1.0e-12
    - name: reduced_integration_introduces_hourglass
      setup: "Hex27 con cuadratura reducida hex_2x2x2"
      expect: "81 - 8·6 = 33 modos cero totales (6 rígidos + 27 hourglass).
              Verificación con jump espectral limpio entre eigs[32] ≈ 1e-10
              y eigs[33] ≈ O(10^3) (escala del módulo de cortante × volumen)"
    - name: mass_consistent_total
      setup: "Hex27 cubo unitario con ρ=7850, masa consistente"
      expect: "Σ M = 3·ρ·V = 23550"
      tol_rel: 1.0e-9
    - name: mass_lumped_total_HRZ
      setup: "Hex27 cubo unitario con ρ=7850, masa lumped HRZ"
      expect: "diag(M) suma = 3·ρ·V; matriz diagonal; todas las entradas
              positivas (HRZ canónico)"
      tol_abs: 1.0e-12

references:
  - "Bathe K.J. (2014). Finite Element Procedures. §5.3.2 (Lagrangiano 3D)."
  - "Cook R.D., Malkus D.S., Plesha M.E., Witt R.J. (2002). Concepts and
    Applications of FEA. §6.6, §11.6 (Q2/E27 con 2×2×2 reducida)."
  - "Zienkiewicz O.C., Taylor R.L. (2005). The Finite Element Method, vol. 1,
    §8.6 (familia Lagrangiana 3D)."
  - "ADR 0012 — Sólidos 3D: convención Voigt 6D y caras."
```

---

## Implementación

- **Archivo**: [solidum/elements/solid_3d/hex27.py](../../solidum/elements/solid_3d/hex27.py).
- **Clase**: `Hex27` (registrada vía `@ElementRegistry.register`); subclase de [`_HigherOrderSolid3D`](../../solidum/elements/solid_3d/_shared.py).
- **Funciones núcleo**: `_N_hex27(xi, eta, zeta)`, `_dN_hex27(xi, eta, zeta)` en [solidum/elements/solid_3d/_shared.py](../../solidum/elements/solid_3d/_shared.py). Implementadas en numpy puro vía productos tensoriales de los polinomios Lagrange 1D `_L_lagrange_quad` y `_dL_lagrange_quad`.
- **Centralización (sub-fase 2)**: `_HigherOrderSolid3D` abstrae el bucle de Gauss para `K`, `F_int`, `gauss_state`, `body_load`, `face_traction` y masa consistente. Comparte con `Hex20` la misma kinematics y la misma maquinaria; sólo cambian las funciones de forma 3D, la cuadratura, las funciones de forma 2D de cara y los `FACE_NODES`.
- **Tests**:
  - `tests/test_solid_3d_higher_order.py` — `TestHex27Element` (patch lineal, patch **triquadrático**, simetría K, jacobiano degenerado, body load, face traction balance, masa consistente y lumped HRZ, reduced integration 1 hourglass).
  - `tests/test_rigid_body_modes.py` — `TestRBMHex27` (3 trans + 3 rot + rank-deficiency = 6).
  - `tests/validation/test_cube_lame_3d.py` — `test_uniaxial_traction_hex27_exact`, `test_hydrostatic_compression_hex27`.
  - `tests/validation/test_macneal_beam_3d.py` — `test_macneal_3d_hex27_*` (convergencia espejo del Hex20).
  - `tests/test_volumetric_locking_3d.py` — `TestHex27VolumetricLocking` (mitigación comparable al Hex20).

---

## Diálogo

- **2026-05-27** · Spec inicial. Sub-fase 2 de A.ter — el segundo caso real (Hex27) dispara la centralización en `_HigherOrderSolid3D` (regla de los dos casos reales aplicada en el momento canónico, igual que el `_HigherOrderSolid2D` cuando entró el `Quad9` después del `Quad8`).
- **2026-05-27** · Implementado y validado. Patch test lineal, bilineal mixto (`u_x = c·xy`) y triquadrático completo (`u_x = c·x²y²z²`) exactos a precisión máquina — esta última es la capacidad distintiva del Hex27 respecto al Hex20 serendípito. K simétrica, 6 modos rígidos, masas consistente y lumped HRZ con todas las entradas positivas, balance de cargas distribuidas correcto. Spec a `status: validated`.
