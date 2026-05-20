# Hex8 — hexaedro trilineal isoparamétrico 3D

> Elemento sólido 3D de primer orden, espejo natural de `Quad4`. Pieza central de la Etapa 7. Convención Voigt 6D y caras fijadas por ADR 0012.

---

## Especificación física

### 0. Descripción general

Elemento sólido tridimensional de primer orden, hexaedro de ocho nodos. Aproximación **trilineal** $C^0$ del campo de desplazamientos en coordenadas naturales $(\xi, \eta, \zeta) \in [-1, 1]^3$. Formulación isoparamétrica: la geometría y los desplazamientos comparten las mismas funciones de forma. Régimen de pequeños desplazamientos y deformaciones.

### 1. Cinemática de desplazamientos

$$u_x(\xi, \eta, \zeta) = \sum_{i=1}^{8} N_i(\xi, \eta, \zeta)\, u_x^{(i)}, \qquad u_y = \sum N_i\, u_y^{(i)}, \qquad u_z = \sum N_i\, u_z^{(i)}.$$

### 2. Cinemática de deformaciones

Tensor infinitesimal en notación Voigt 6D del proyecto (ADR 0012):

$$\boldsymbol\varepsilon = [\varepsilon_{xx},\ \varepsilon_{yy},\ \varepsilon_{zz},\ \gamma_{xy},\ \gamma_{yz},\ \gamma_{xz}]^\top$$

con

$$\varepsilon_{xx} = \tfrac{\partial u_x}{\partial x},\quad \varepsilon_{yy} = \tfrac{\partial u_y}{\partial y},\quad \varepsilon_{zz} = \tfrac{\partial u_z}{\partial z},$$

$$\gamma_{xy} = \tfrac{\partial u_x}{\partial y} + \tfrac{\partial u_y}{\partial x},\quad \gamma_{yz} = \tfrac{\partial u_y}{\partial z} + \tfrac{\partial u_z}{\partial y},\quad \gamma_{xz} = \tfrac{\partial u_x}{\partial z} + \tfrac{\partial u_z}{\partial x}.$$

### 3. Ecuación constitutiva

Delegada al material 3D (`STRAIN_DIM = 6`). El elemento llama `material.compute_state(ε, state) → (σ, C_tangent, state')` por punto de Gauss. En esta etapa solo `Elastic3D` está disponible; el contrato deja la puerta abierta a futuros `VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`.

### 4. Equilibrio — forma fuerte

$$-\nabla\cdot\boldsymbol\sigma = \mathbf b \quad \text{en } \Omega, \qquad \boldsymbol\sigma\cdot\mathbf n = \bar{\mathbf t} \quad \text{en } \partial\Omega_N, \qquad \mathbf u = \bar{\mathbf u} \quad \text{en } \partial\Omega_D.$$

### 5. Equilibrio — forma débil

$$\int_\Omega \delta\boldsymbol\varepsilon^\top \boldsymbol\sigma\, dV = \int_\Omega \delta\mathbf u^\top \mathbf b\, dV + \int_{\partial\Omega_N} \delta\mathbf u^\top \bar{\mathbf t}\, dS, \qquad \forall\, \delta\mathbf u \in V_0.$$

---

## Formulación numérica (FEM)

### 6. Discretización

Ocho nodos en orden estándar de hexaedro VTK (compatible con la mayoría de pre-procesadores). En el cubo de referencia $[-1, 1]^3$:

| Nodo local | $(\xi_i, \eta_i, \zeta_i)$ |
|---|---|
| 0 | $(-1, -1, -1)$ |
| 1 | $(+1, -1, -1)$ |
| 2 | $(+1, +1, -1)$ |
| 3 | $(-1, +1, -1)$ |
| 4 | $(-1, -1, +1)$ |
| 5 | $(+1, -1, +1)$ |
| 6 | $(+1, +1, +1)$ |
| 7 | $(-1, +1, +1)$ |

Cara inferior $\zeta = -1$ recorrida `(0, 1, 2, 3)` antihorario visto desde dentro; cara superior $\zeta = +1$ recorrida `(4, 5, 6, 7)` antihorario visto desde fuera. Mapeo isoparamétrico desde el cubo de referencia:

$$x(\xi, \eta, \zeta) = \sum_{i=0}^{7} N_i\, x^{(i)}, \quad y = \sum N_i\, y^{(i)}, \quad z = \sum N_i\, z^{(i)}.$$

Jacobiano $\mathbf J = \partial(x, y, z) / \partial(\xi, \eta, \zeta)$ es una matriz 3×3 con determinante $\det\mathbf J$. El código aborta con `ValueError` si $\det\mathbf J \le \text{tol}$ (elemento degenerado o nodos en orden invertido).

### 7. Funciones de forma

Producto tensorial trilineal:

$$N_i(\xi, \eta, \zeta) = \tfrac{1}{8}(1 + \xi_i\,\xi)(1 + \eta_i\,\eta)(1 + \zeta_i\,\zeta), \qquad (\xi_i, \eta_i, \zeta_i) \in \{-1, +1\}^3.$$

### 8. Matriz deformación-desplazamiento

Las derivadas en globales se obtienen vía $\partial N_i / \partial \mathbf x = \mathbf J^{-1}\,\partial N_i / \partial \boldsymbol\xi$. La matriz $\mathbf B$ tiene tamaño **6×24** y se ensambla por nodo $i \in \{0,\dots,7\}$:

$$\mathbf B_i = \begin{bmatrix}
\partial N_i / \partial x & 0 & 0 \\
0 & \partial N_i / \partial y & 0 \\
0 & 0 & \partial N_i / \partial z \\
\partial N_i / \partial y & \partial N_i / \partial x & 0 \\
0 & \partial N_i / \partial z & \partial N_i / \partial y \\
\partial N_i / \partial z & 0 & \partial N_i / \partial x
\end{bmatrix}, \qquad \boldsymbol\varepsilon = \mathbf B\,\mathbf u_e.$$

El orden de las filas reproduce exactamente el orden Voigt del proyecto (ADR 0012). El orden de las columnas por nodo es `[ux, uy, uz]`, agrupado por nodo $i$ para `i ∈ {0..7}`.

### 9. Rigidez elemental

$$\mathbf K_e = \int_\Omega \mathbf B^\top\,\mathbf D\,\mathbf B\, dV = \sum_{p=1}^{n_g} \mathbf B(\boldsymbol\xi_p)^\top\,\mathbf D_p\,\mathbf B(\boldsymbol\xi_p)\,\det\mathbf J_p\, w_p.$$

A diferencia del 2D, no hay parámetro de espesor — el "espesor" está incluido en el volumen integrado.

### 10. Fuerzas internas

$$\mathbf F_{\text{int}} = \sum_{p=1}^{n_g} \mathbf B(\boldsymbol\xi_p)^\top\,\boldsymbol\sigma_p\,\det\mathbf J_p\, w_p.$$

### 11. Cuadratura

Por defecto **Gauss–Legendre 2×2×2** (ocho puntos): orden de exactitud 3 en cada dirección — integra exactamente $\mathbf K_e$ para Jacobiano constante (geometría de paralelepípedo) y captura el comportamiento dominante en geometrías irregulares moderadas. Soporta esquemas alternativos vía el parámetro `quadrature` desde `QuadratureRegistry` (e.g. `hex_3x3x3` para problemas con materiales no lineales severos).

**Aviso sobre integración reducida**: el esquema $1 \times 1 \times 1$ (un solo punto central) está disponible pero el código emite advertencia explícita; un único punto no detecta los modos de hourglass (energía nula con desplazamientos alternados) que en 3D son **doce** modos espurios (no cuatro como en 2D). No se implementa estabilización tipo Flanagan-Belytschko. Recomendación operativa: usar el default $2 \times 2 \times 2$ salvo necesidad explícita.

### 12. Cargas distribuidas consistentes

**Fuerza de cuerpo $\mathbf b$** (e.g. peso propio $\mathbf b = (0, 0, -\rho g)$). Vector nodal equivalente:

$$\mathbf f_e^{b} = \int_{\Omega_e} \mathbf N^\top\,\mathbf b\, dV = \sum_{p=1}^{n_g} \mathbf N(\boldsymbol\xi_p)^\top\,\mathbf b\,\det\mathbf J_p\, w_p,$$

donde $\mathbf N$ es la matriz $3 \times 24$ de funciones de forma. Se reutiliza la cuadratura del elemento (Gauss $2 \times 2 \times 2$ por defecto). Exacto para $\mathbf b$ uniforme y geometría de paralelepípedo.

**Tracción de superficie $\bar{\mathbf t}$** sobre una cara $\Gamma_e \subset \partial\Omega_N$. Caras numeradas con normal saliente (ADR 0012):

| Cara | Nodos | Normal saliente |
|------|-------|-----------------|
| 0 | (0, 3, 2, 1) | −ζ (inferior) |
| 1 | (4, 5, 6, 7) | +ζ (superior) |
| 2 | (0, 1, 5, 4) | −η (frontal) |
| 3 | (1, 2, 6, 5) | +ξ (derecha) |
| 4 | (2, 3, 7, 6) | +η (trasera) |
| 5 | (3, 0, 4, 7) | −ξ (izquierda) |

$$\mathbf f_e^{t} = \int_{\Gamma_e} \mathbf N^\top\,\bar{\mathbf t}\, dS.$$

Cuadratura 2D sobre la cara: Gauss $2 \times 2$ (cuatro puntos por cara). Exacto para tracción constante sobre cara plana. $\bar{\mathbf t}$ se especifica en **coordenadas globales** $(t_x, t_y, t_z)$; presiones normales se obtienen multiplicando previamente por la normal exterior de la cara. Tracción variable a lo largo de la cara no soportada en este paso.

### 13. Salida por punto de Gauss

`compute_gauss_state(U)` devuelve un dict con $\boldsymbol\varepsilon$ y $\boldsymbol\sigma$ en cada uno de los $n_g$ puntos de Gauss del elemento (8 por defecto), junto con coordenadas naturales y globales. Habilita post-proceso fino (mapas no promediados, suavizado nodal, recovery superconvergente futuro). Por ADR 0012, `internal_forces` devuelve `None` (sólidos no tienen "esfuerzo seccional").

---

## Limitaciones declaradas

### Locking volumétrico con $\nu \to 0.5$

`Hex8` con integración completa $2\times 2\times 2$ y material casi-incompresible (`Elastic3D` con $\nu \to 0.5$, o un futuro `VonMises3D` en régimen plástico desarrollado) sufre **locking volumétrico**: la rigidez crece artificialmente porque el espacio de funciones trilineal no contiene desplazamientos puramente isocóricos suficientes. El elemento se queda atrapado en respuestas erróneamente rígidas.

Tratamiento en este proyecto: declarar limitación, **blindar con test** (`test_volumetric_locking_3d.py`, espejo del análogo 2D) que documente cuantitativamente el colapso, **no implementar mitigación** (B-bar, F-bar, mixed u-p). Política idéntica al `Quad4` (ver `STATUS.md` "Limitaciones declaradas").

### Hourglass con integración reducida

`Hex8` integrado con un solo punto central tiene **12 modos de hourglass** (energía nula con desplazamientos alternados). Reducción disponible vía `quadrature` con warning; estabilización Flanagan-Belytschko no implementada.

### Mass lumping en orientación arbitraria

`compute_mass_matrix(lumping="lumped")` usa HRZ canónico (`solidum/math/mass_lumping.py::lump_hrz`). Diagonal en ejes locales del elemento. En sólidos 3D la matriz lumped traslacional es estrictamente diagonal porque `m_t·I₃` es invariante bajo SO(3) — sin el caveat de `Frame3D` que sí aparece en su bloque rotacional (no hay rotaciones nodales en sólidos).

---

## Contrato de implementación

```yaml
name: Hex8
kind: element
status: validated

interface:
  dof_names: [ux, uy, uz]
  n_nodes: 8
  strain_dim: 6                 # Voigt 3D [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]
  n_integration_points: 8       # default Gauss 2×2×2

parameters:
  - { name: quadrature, type: str, required: false, default: "hex_2x2x2",
      desc: "Regla de cuadratura desde QuadratureRegistry (hex_2x2x2 default, hex_3x3x3 alternativa)" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, C_tangent, state')"
  strain_kind: "Voigt 3D [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

conventions:
  sign: "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "[xx, yy, zz, xy, yz, xz] con γ_ij = 2·ε_ij (ADR 0012)"
  node_orientation: "convención VTK_HEXAHEDRON; nodos 0-3 cara inferior antihorario visto desde dentro, 4-7 cara superior antihorario visto desde fuera; asegura det J > 0"
  face_numbering: "ADR 0012 — 6 caras con normal saliente: 0 (−ζ), 1 (+ζ), 2 (−η), 3 (+ξ), 4 (+η), 5 (−ξ)"

validity:
  - "pequeños desplazamientos y deformaciones"
  - "geometría sin distorsión severa (det J > tol en todos los puntos de Gauss)"
  - "compatible con materiales con STRAIN_DIM = 6 (en esta etapa: Elastic3D)"

out_of_scope:
  - "grandes deformaciones, formulaciones corotacionales o lagrangianas totales"
  - "estabilización contra hourglass para integración reducida"
  - "elementos de alto orden (Hex20, Hex27 — sub-etapa posterior)"
  - "tracción variable sobre cara o presión normal"
  - "internal_forces (ADR 0012: sólidos exponen compute_gauss_state)"

acceptance:
  verification:
    - name: patch_test_3d
      setup: "8 Hex8 alrededor de un nodo interior (cubo 2×2×2 partido en
             ocho octantes) con campo lineal u_i = a_ij·x_j impuesto en
             los 26 nodos del contorno; nodo interior libre"
      expect: "nodo interior adopta el campo lineal exacto y ε constante e
              igual a 0.5·(a_ij + a_ji) en todos los puntos de Gauss; σ
              uniforme"
      tol_rel: 1.0e-10
    - name: modos_de_cuerpo_rigido
      setup: "evaluar K_e sobre un Hex8 cualquiera y proyectar sobre los 6
             modos de cuerpo rígido (3 traslaciones + 3 rotaciones)"
      expect: "K_e · u_rigid ≈ 0 en cada uno de los 6 modos"
      tol_abs: 1.0e-9
    - name: simetria_K_e
      setup: "evaluar K_e con Elastic3D"
      expect: "K_e == K_e.T exacto"
      tol_abs: 1.0e-12
  specific:
    - name: cubo_lame_3d
      setup: "cubo unitario con Elastic3D; tracción uniaxial σ_xx = p en
             cara +ξ; cara opuesta restringida en ux; restricciones
             mínimas para evitar movimiento rígido"
      expect: "u_x(L, y, z) = p·L/E; u_y = -ν·p·y/E; u_z = -ν·p·z/E
              exactos en todos los nodos"
      tol_rel: 1.0e-10
    - name: traccion_uniaxial_sobre_cara
      setup: "Hex8 con tracción constante (p, 0, 0) sobre cara 3 (+ξ)"
      expect: "Σf_x = p·A_cara, Σf_y = Σf_z = 0; reparto p·A_cara/4 a cada
              nodo de la cara para Hex8 regular"
      tol_rel: 1.0e-12
    - name: body_load_uniforme
      setup: "Hex8 regular y distorsionado con b uniforme"
      expect: "Σf = b·V_e (invariante de geometría); por simetría en Hex8
              regular cada nodo recibe b·V_e / 8"
      tol_rel: 1.0e-10
    - name: jacobiano_degenerado_abortado
      setup: "Hex8 con cuatro nodos coplanares de la cara superior sobre
             la inferior (det J ≈ 0)"
      expect: "ValueError en _compute_kinematics"
    - name: locking_volumetrico_documentado
      setup: "viga esbelta empotrada-libre con Hex8 plane strain forzado
             (cara z restringida); ν ∈ {0.3, 0.49, 0.4999}"
      expect: "deflexión vs Frame3D decrece monotonamente con ν → 0.5
              (locking documentado, no mitigado)"
      tol_rel: "documentar valores, sin tolerancia de aceptación"

references:
  - "Bathe K.J. (2014). Finite Element Procedures. §6.6 (sólidos 3D isoparamétricos)."
  - "Cook R.D., Malkus D.S., Plesha M.E., Witt R.J. (2002). Concepts and Applications of FEA. §6.8."
  - "Zienkiewicz O.C., Taylor R.L. (2005). The Finite Element Method, vol. 1, §6.7."
  - "ADR 0012 — Sólidos 3D: convención Voigt 6D y caras."
```

---

## Implementación

- **Archivo**: [solidum/elements/solid_3d/hex8.py](../../solidum/elements/solid_3d/hex8.py).
- **Clase**: `Hex8` (registrada vía `@ElementRegistry.register`).
- **Funciones núcleo**: `_compute_kinematics_hex8(xi, eta, zeta, coords)`, `_det_jacobian_hex8(...)` y `_shape_functions_hex8(...)` con `@njit` (Numba) en [solidum/elements/solid_3d/_shared.py](../../solidum/elements/solid_3d/_shared.py). Integrandos comunes 3D (`_compute_integrands_3d`) y expansor masa traslacional (`_expand_scalar_mass_3d`) también en `_shared`.
- **Tests**:
  - [tests/test_solid_3d.py](../../tests/test_solid_3d.py) — clase `TestHex8Element` (10 tests: dimensiones/DOFs, simetría K, patch tracción uniaxial, jacobiano negativo abortado, body load, face traction (balance + nodos activos + índice fuera de rango), masa consistente y lumped HRZ, gauss state).
  - [tests/test_rigid_body_modes.py](../../tests/test_rigid_body_modes.py) — `TestRBMHex8` (5 tests: traslación, 3 rotaciones independientes en (x, y, z), rank-deficiency = 6 modos rígidos).
  - [tests/validation/test_cube_lame_3d.py](../../tests/validation/test_cube_lame_3d.py) — `test_uniaxial_traction_hex8_exact`, `test_hydrostatic_compression_hex8` (solución exacta a precisión máquina).
  - [tests/validation/test_macneal_beam_3d.py](../../tests/validation/test_macneal_beam_3d.py) — shear locking documentado en malla coarse + convergencia h monótona.
  - [tests/test_volumetric_locking_3d.py](../../tests/test_volumetric_locking_3d.py) — locking volumétrico declarado y blindado por test.

---

## Diálogo

- **2026-05-19** · Spec inicial. Elemento espejo natural de `Quad4`; convención Voigt 6D y numeración de caras fijadas por ADR 0012. Sin variantes corotacionales ni de alto orden en esta etapa. Locking volumétrico declarado como limitación con test que lo blinda; sin mitigación implementada (política idéntica al 2D).
