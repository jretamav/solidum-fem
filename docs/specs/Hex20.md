# Hex20 — hexaedro serendípito 3D de orden 2 (20 nodos)

> Elemento sólido 3D de orden cuadrático con interpolación serendípita.
> Reproduce campos cuadráticos completos en 3D (todos los polinomios de grado total ≤ 2 más algunos de orden 3 sin los términos triquadráticos). Es al `Hex8` lo que el `Quad8` al `Quad4`: mucho más preciso en flexión y geometrías curvas con el mismo número de elementos. Primer elemento de la sub-etapa A.ter.

---

## Especificación física

### 0. Descripción general

Hexaedro tridimensional de segundo orden con **20 nodos**: 8 vértices + 12 medios de arista, sin nodos de cara ni centroide. Aproximación serendípita $C^0$ en coordenadas naturales $(\xi, \eta, \zeta) \in [-1, 1]^3$. Formulación isoparamétrica: la geometría y los desplazamientos comparten las mismas funciones de forma. Régimen de pequeños desplazamientos y deformaciones.

### 1. Cinemática de desplazamientos

$$u_x(\xi, \eta, \zeta) = \sum_{i=0}^{19} N_i(\xi, \eta, \zeta)\, u_x^{(i)}, \qquad u_y = \sum N_i\, u_y^{(i)}, \qquad u_z = \sum N_i\, u_z^{(i)}.$$

### 2. Cinemática de deformaciones

Tensor infinitesimal en notación Voigt 6D del proyecto (ADR 0012):

$$\boldsymbol\varepsilon = [\varepsilon_{xx},\ \varepsilon_{yy},\ \varepsilon_{zz},\ \gamma_{xy},\ \gamma_{yz},\ \gamma_{xz}]^\top$$

con $\gamma_{ij} = 2\,\varepsilon_{ij}$ *engineering*. Idéntica al `Hex8`.

### 3. Ecuación constitutiva

Delegada al material 3D (`STRAIN_DIM = 6`). El elemento llama `material.compute_state(ε, state) → (σ, C_tangent, state')` por punto de Gauss. Compatible con todo el catálogo 3D actual: `Elastic3D`, `VonMises3D`, `DruckerPrager3D`, `IsotropicDamage3D`.

### 4-5. Equilibrio

Forma fuerte y débil idénticas al `Hex8` — la diferencia entre `Hex8` y `Hex20` está exclusivamente en el espacio de aproximación discreto.

---

## Formulación numérica (FEM)

### 6. Discretización

Veinte nodos en orden VTK_QUADRATIC_HEXAHEDRON (estándar gmsh/VTK/Abaqus). En el cubo de referencia $[-1, 1]^3$:

| Nodo local | $(\xi_i, \eta_i, \zeta_i)$ | Tipo |
|---|---|---|
| 0  | $(-1, -1, -1)$ | vértice |
| 1  | $(+1, -1, -1)$ | vértice |
| 2  | $(+1, +1, -1)$ | vértice |
| 3  | $(-1, +1, -1)$ | vértice |
| 4  | $(-1, -1, +1)$ | vértice |
| 5  | $(+1, -1, +1)$ | vértice |
| 6  | $(+1, +1, +1)$ | vértice |
| 7  | $(-1, +1, +1)$ | vértice |
| 8  | $(0, -1, -1)$ | medio arista 0-1 (eje ξ, cara inferior) |
| 9  | $(+1, 0, -1)$ | medio arista 1-2 (eje η, cara inferior) |
| 10 | $(0, +1, -1)$ | medio arista 2-3 (eje ξ, cara inferior) |
| 11 | $(-1, 0, -1)$ | medio arista 3-0 (eje η, cara inferior) |
| 12 | $(0, -1, +1)$ | medio arista 4-5 (eje ξ, cara superior) |
| 13 | $(+1, 0, +1)$ | medio arista 5-6 (eje η, cara superior) |
| 14 | $(0, +1, +1)$ | medio arista 6-7 (eje ξ, cara superior) |
| 15 | $(-1, 0, +1)$ | medio arista 7-4 (eje η, cara superior) |
| 16 | $(-1, -1, 0)$ | medio arista 0-4 (eje ζ, vertical) |
| 17 | $(+1, -1, 0)$ | medio arista 1-5 (eje ζ, vertical) |
| 18 | $(+1, +1, 0)$ | medio arista 2-6 (eje ζ, vertical) |
| 19 | $(-1, +1, 0)$ | medio arista 3-7 (eje ζ, vertical) |

Mapeo isoparamétrico desde el cubo de referencia:

$$x(\xi, \eta, \zeta) = \sum_{i=0}^{19} N_i\, x^{(i)},\quad y = \sum N_i\, y^{(i)},\quad z = \sum N_i\, z^{(i)}.$$

Jacobiano $\mathbf J = \partial(x, y, z) / \partial(\xi, \eta, \zeta)$ es una matriz 3×3 con determinante $\det\mathbf J$. El código aborta con `ValueError` si $\det\mathbf J \le \text{tol}$ en cualquier punto de Gauss (elemento severamente distorsionado o nodos en orden invertido).

### 7. Funciones de forma serendípitas

**Vértices** ($i \in \{0,\dots,7\}$, $(\xi_i, \eta_i, \zeta_i) \in \{\pm 1\}^3$):

$$N_i = \frac{1}{8}(1 + \xi_i\xi)(1 + \eta_i\eta)(1 + \zeta_i\zeta)(\xi_i\xi + \eta_i\eta + \zeta_i\zeta - 2).$$

El multiplicador $(\xi_i\xi + \eta_i\eta + \zeta_i\zeta - 2)$ vale $1$ en el propio vértice y $0$ en los demás vértices y en todos los medios de arista — garantía $\delta_{ij}$ del esquema.

**Medios de arista paralelos a $\xi$** (nodos 8, 10, 12, 14; $\xi_i = 0$, $\eta_i, \zeta_i \in \{\pm 1\}$):

$$N_i = \frac{1}{4}(1 - \xi^2)(1 + \eta_i\eta)(1 + \zeta_i\zeta).$$

**Medios de arista paralelos a $\eta$** (nodos 9, 11, 13, 15; $\eta_i = 0$):

$$N_i = \frac{1}{4}(1 + \xi_i\xi)(1 - \eta^2)(1 + \zeta_i\zeta).$$

**Medios de arista paralelos a $\zeta$** (nodos 16, 17, 18, 19; $\zeta_i = 0$):

$$N_i = \frac{1}{4}(1 + \xi_i\xi)(1 + \eta_i\eta)(1 - \zeta^2).$$

Espacio reproducido: $\{1, \xi, \eta, \zeta, \xi^2, \eta^2, \zeta^2, \xi\eta, \eta\zeta, \xi\zeta, \xi^2\eta, \xi\eta^2, \eta^2\zeta, \eta\zeta^2, \xi^2\zeta, \xi\zeta^2, \xi\eta\zeta, \xi^2\eta\zeta, \xi\eta^2\zeta, \xi\eta\zeta^2\}$ — 20 monomios, omitiendo los puramente cúbicos ($\xi^3$, etc.) y el término $\xi^2\eta^2\zeta^2$.

### 8. Matriz deformación-desplazamiento

Las derivadas en globales se obtienen vía $\partial N_i / \partial \mathbf x = \mathbf J^{-1}\,\partial N_i / \partial \boldsymbol\xi$. La matriz $\mathbf B$ tiene tamaño **6×60** y se ensambla por nodo $i \in \{0,\dots,19\}$ con la misma plantilla 6×3 del `Hex8`:

$$\mathbf B_i = \begin{bmatrix}
\partial N_i / \partial x & 0 & 0 \\
0 & \partial N_i / \partial y & 0 \\
0 & 0 & \partial N_i / \partial z \\
\partial N_i / \partial y & \partial N_i / \partial x & 0 \\
0 & \partial N_i / \partial z & \partial N_i / \partial y \\
\partial N_i / \partial z & 0 & \partial N_i / \partial x
\end{bmatrix}, \qquad \boldsymbol\varepsilon = \mathbf B\,\mathbf u_e.$$

El orden de las filas reproduce exactamente el orden Voigt del proyecto (ADR 0012). El orden de las columnas por nodo es `[ux, uy, uz]`, agrupado por nodo $i$ para `i ∈ {0..19}`.

**Derivadas explícitas** (cerradas vía regla del producto):

Vértices ($i \in \{0,\dots,7\}$, con $u = \xi_i\xi$, $v = \eta_i\eta$, $w = \zeta_i\zeta$):

$$\frac{\partial N_i}{\partial \xi} = \frac{\xi_i}{8}(1 + v)(1 + w)(2u + v + w - 1),$$

$$\frac{\partial N_i}{\partial \eta} = \frac{\eta_i}{8}(1 + u)(1 + w)(u + 2v + w - 1),$$

$$\frac{\partial N_i}{\partial \zeta} = \frac{\zeta_i}{8}(1 + u)(1 + v)(u + v + 2w - 1).$$

Medio paralelo a $\xi$ (nodos 8, 10, 12, 14):

$$\frac{\partial N_i}{\partial \xi} = -\frac{\xi}{2}(1 + \eta_i\eta)(1 + \zeta_i\zeta),\quad \frac{\partial N_i}{\partial \eta} = \frac{\eta_i}{4}(1 - \xi^2)(1 + \zeta_i\zeta),\quad \frac{\partial N_i}{\partial \zeta} = \frac{\zeta_i}{4}(1 - \xi^2)(1 + \eta_i\eta).$$

Medio paralelo a $\eta$ (nodos 9, 11, 13, 15):

$$\frac{\partial N_i}{\partial \xi} = \frac{\xi_i}{4}(1 - \eta^2)(1 + \zeta_i\zeta),\quad \frac{\partial N_i}{\partial \eta} = -\frac{\eta}{2}(1 + \xi_i\xi)(1 + \zeta_i\zeta),\quad \frac{\partial N_i}{\partial \zeta} = \frac{\zeta_i}{4}(1 + \xi_i\xi)(1 - \eta^2).$$

Medio paralelo a $\zeta$ (nodos 16, 17, 18, 19):

$$\frac{\partial N_i}{\partial \xi} = \frac{\xi_i}{4}(1 + \eta_i\eta)(1 - \zeta^2),\quad \frac{\partial N_i}{\partial \eta} = \frac{\eta_i}{4}(1 + \xi_i\xi)(1 - \zeta^2),\quad \frac{\partial N_i}{\partial \zeta} = -\frac{\zeta}{2}(1 + \xi_i\xi)(1 + \eta_i\eta).$$

### 9. Rigidez elemental

$$\mathbf K_e = \int_\Omega \mathbf B^\top\,\mathbf D\,\mathbf B\, dV = \sum_{p=1}^{n_g} \mathbf B(\boldsymbol\xi_p)^\top\,\mathbf D_p\,\mathbf B(\boldsymbol\xi_p)\,\det\mathbf J_p\, w_p.$$

Sin parámetro de espesor — el volumen ya está incluido en $\det\mathbf J_p\cdot w_p$.

### 10. Fuerzas internas

$$\mathbf F_{\text{int}} = \sum_{p=1}^{n_g} \mathbf B(\boldsymbol\xi_p)^\top\,\boldsymbol\sigma_p\,\det\mathbf J_p\, w_p.$$

### 11. Cuadratura

Por defecto **Gauss–Legendre 3×3×3** (27 puntos): orden de exactitud 5 en cada dirección. Suficiente para integrar exactamente $\mathbf B^\top \mathbf D\, \mathbf B$ del `Hex20` sobre Jacobiano constante (el integrando es a lo sumo de grado 3 en cada dirección) y para integrar exactamente la masa consistente $\rho \mathbf N^\top \mathbf N$ (grado 4 en cada dirección).

**Aviso sobre integración reducida**: el esquema $2 \times 2 \times 2$ (`hex_2x2x2`, 8 puntos) está disponible y reduce el coste computacional, pero introduce **6 modos de hourglass espurios** sobre un elemento aislado (verificación numérica directa del proyecto; Belytschko-Liu-Moran §8.4.4). Estos modos **no propagan** a través de mallas ensambladas con condiciones de Dirichlet razonables — la rigidez global ensamblada queda rank-sufficient en la mayoría de casos prácticos, pero un elemento aislado o una capa con BC pobres puede mostrar hourglass. El código permite seleccionar la cuadratura reducida vía el parámetro `quadrature` pero **no estabiliza** los modos espurios. Recomendación operativa: usar el default $3 \times 3 \times 3$ salvo necesidad experimental. La matriz de masa se integra con $3 \times 3 \times 3$ independientemente de la cuadratura elegida para $\mathbf K$, para evitar subintegrarla.

### 12. Cargas distribuidas consistentes

**Fuerza de cuerpo $\mathbf b$**. Vector nodal equivalente:

$$\mathbf f_e^{b} = \int_{\Omega_e} \mathbf N^\top\,\mathbf b\, dV = \sum_{p=1}^{n_g} \mathbf N(\boldsymbol\xi_p)^\top\,\mathbf b\,\det\mathbf J_p\, w_p,$$

donde $\mathbf N$ es la matriz $3 \times 60$ de funciones de forma. Se reutiliza la cuadratura del elemento. Para $\mathbf b$ uniforme y geometría regular el resultado es exacto. **Para `Hex20` el reparto no es uniforme entre nodos**: la fórmula serendípita asigna pesos distintos a vértices y medios (los medios reciben más carga total que los vértices, con suma global $\rho\cdot V_e$ correcta).

**Tracción de superficie $\bar{\mathbf t}$** sobre una cara. Cada cara del `Hex20` tiene **8 nodos** (4 vértices + 4 medios de arista) y se aproxima por funciones serendípitas Quad8 sobre el parámetro 2D $(s, t) \in [-1, 1]^2$. La numeración de caras conserva la convención del `Hex8` (ADR 0012) extendida con los medios:

| Cara | Vértices  | Medios de arista          | Normal saliente |
|------|-----------|---------------------------|-----------------|
| 0    | (0,3,2,1) | (11, 10, 9, 8)            | −ζ (inferior)   |
| 1    | (4,5,6,7) | (12, 13, 14, 15)          | +ζ (superior)   |
| 2    | (0,1,5,4) | (8, 17, 12, 16)           | −η (frontal)    |
| 3    | (1,2,6,5) | (9, 18, 13, 17)           | +ξ (derecha)    |
| 4    | (2,3,7,6) | (10, 19, 14, 18)          | +η (trasera)    |
| 5    | (3,0,4,7) | (11, 16, 15, 19)          | −ξ (izquierda)  |

El orden dentro de cada cara reproduce la convención VTK_QUADRATIC_HEXAHEDRON: primero los 4 vértices, luego los 4 medios en el mismo orden cíclico (medio entre vértice $k$ y vértice $k+1$).

$$\mathbf f_e^{t} = \int_{\Gamma_e} \mathbf N^\top\,\bar{\mathbf t}\, dS.$$

Cuadratura 2D sobre la cara: Gauss $3 \times 3$ (nueve puntos) — exacta para tracción constante sobre cara plana (integrando degree 2 en $(s,t)$). $\bar{\mathbf t}$ se especifica en **coordenadas globales** $(t_x, t_y, t_z)$; presiones normales se obtienen multiplicando previamente por la normal exterior de la cara. Tracción variable sobre la cara no soportada en este paso. **Para `Hex20` el reparto entre vértices y medios no es trivial** (la integral cuadrática de Quad8 sobre cara plana asigna fracciones específicas: los vértices reciben $-1/12$ de la fuerza total cada uno y los medios $4/12$ cada uno, con suma global exacta — ver test de balance).

### 13. Salida por punto de Gauss

`compute_gauss_state(U)` devuelve un dict con $\boldsymbol\varepsilon$ y $\boldsymbol\sigma$ en cada uno de los $n_g$ puntos de Gauss del elemento (27 por defecto), junto con coordenadas naturales y globales. Habilita post-proceso fino (mapas no promediados, suavizado nodal, recovery superconvergente futuro). Por ADR 0012, `internal_forces` devuelve `None` (sólidos no tienen fuerzas internas seccionales discretas).

---

## Limitaciones declaradas

### Locking volumétrico con $\nu \to 0.5$

`Hex20` con integración completa $3\times 3\times 3$ y material casi-incompresible (`Elastic3D` con $\nu \to 0.5$, o `VonMises3D` en régimen plástico desarrollado) sufre **locking volumétrico** menos severo que el `Hex8` pero todavía presente: el espacio serendípito triquadrático tampoco contiene desplazamientos puramente isocóricos suficientes en mallas distorsionadas (Cook-Malkus-Plesha §13.6). 

Tratamiento en este proyecto: declarar limitación, **blindar con test cuantitativo** (`test_volumetric_locking_hex20.py`, espejo del análogo `Hex8`) que documente la mitigación parcial respecto al lineal, **no implementar B-bar/F-bar/mixed u-p**. Política idéntica a `Quad4`/`Quad8` y `Hex8` (ver `STATUS.md` "Limitaciones declaradas").

### Hourglass con integración reducida

`Hex20` integrado con `hex_2x2x2` tiene **6 modos de hourglass espurios** por elemento aislado (vs 12 del `Hex8` reducido). En mallas ensambladas con Dirichlet razonable los modos suelen no propagar — la rigidez global queda rank-sufficient — pero un elemento aislado o una capa con BC pobres lo mostraría. Reducción disponible vía `quadrature`; estabilización no implementada. Modo recomendado: full integration $3\times 3\times 3$.

### Mass lumping

`compute_mass_matrix(lumping="lumped")` usa HRZ canónico (`solidum/math/mass_lumping.py::lump_hrz`). Para `Hex20` (con nodos de medio de arista) HRZ es la única opción razonable: row-sum produciría masas **negativas** en los vértices (los pesos de Gauss son negativos para los vértices en la cuadratura serendípita Lobatto equivalente), inadmisible para análisis dinámico (Bathe FEP §9.2.4). HRZ preserva masa total escalando la diagonal consistente por un factor común.

---

## Contrato de implementación

```yaml
name: Hex20
kind: element
status: validated

interface:
  dof_names: [ux, uy, uz]
  n_nodes: 20
  strain_dim: 6                 # Voigt 3D [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]
  n_integration_points: 27      # default Gauss 3×3×3

parameters:
  - { name: quadrature, type: str, required: false, default: "hex_3x3x3",
      desc: "Regla de cuadratura desde QuadratureRegistry (hex_3x3x3 default,
             hex_2x2x2 reducida con 4 modos de hourglass espurios)" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, C_tangent, state')"
  strain_kind: "Voigt 3D [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

conventions:
  sign: "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "[xx, yy, zz, xy, yz, xz] con γ_ij = 2·ε_ij (ADR 0012)"
  node_orientation: "convención VTK_QUADRATIC_HEXAHEDRON; vértices 0-7 idénticos
                     al Hex8; medios 8-11 bottom edges, 12-15 top edges, 16-19
                     vertical edges; asegura det J > 0 en todos los Gauss"
  face_numbering: "ADR 0012 — 6 caras con normal saliente, idéntica a Hex8 en
                   los vértices; cada cara añade 4 nodos medios en el mismo
                   orden cíclico (Quad8 serendipity)"

validity:
  - "pequeños desplazamientos y deformaciones"
  - "geometría sin distorsión severa (det J > tol en todos los puntos de Gauss)"
  - "compatible con materiales con STRAIN_DIM = 6 (Elastic3D, VonMises3D,
     DruckerPrager3D, IsotropicDamage3D)"

out_of_scope:
  - "grandes deformaciones, formulaciones corotacionales o lagrangianas totales"
  - "estabilización contra hourglass para integración reducida hex_2x2x2"
  - "tracción variable sobre cara o presión normal con gradiente"
  - "internal_forces (ADR 0012: sólidos exponen compute_gauss_state)"
  - "Hex27 lagrangiano (sub-fase 2 de A.ter)"

acceptance:
  verification:
    - name: patch_test_lineal
      setup: "Hex20 con campo lineal u_i = a_ij·x_j impuesto en los 20 nodos
             del contorno; sin nodos interiores libres (todos están en la
             frontera de un elemento aislado, pero al menos los medios no
             están fijados en su posición referencial)"
      expect: "ε constante e igual a sym(a) en todos los puntos de Gauss; σ
              uniforme"
      tol_rel: 1.0e-12
    - name: patch_test_cuadratico
      setup: "Hex20 con campo u_x = c·x², u_y = u_z = 0 impuesto en los 20 nodos
             del contorno (todos)"
      expect: "ε_xx(x) = 2c·x exacto en todos los puntos de Gauss"
      tol_rel: 1.0e-12
    - name: modos_de_cuerpo_rigido
      setup: "K_e del Hex20 proyectado sobre los 6 modos de cuerpo rígido"
      expect: "K_e · u_rigid ≈ 0 (autovalores ≈ 0: exactamente 6)"
      tol_abs: 1.0e-9
    - name: simetria_K_e
      setup: "evaluar K_e con Elastic3D"
      expect: "K_e == K_e.T exacto"
      tol_abs: 1.0e-12
  specific:
    - name: cubo_lame_3d
      setup: "cubo unitario con un Hex20; tracción uniaxial σ_xx = p en cara +ξ;
             cara opuesta restringida en ux; restricciones mínimas para evitar
             movimiento rígido"
      expect: "u_x(L, y, z) = p·L/E; u_y = −ν·p·y/E; u_z = −ν·p·z/E exactos en
              todos los nodos"
      tol_rel: 1.0e-10
    - name: body_load_uniforme
      setup: "Hex20 regular cubo unitario con b = (0, 0, −g)"
      expect: "Σf_z = −g·V_e = −g; reparto no uniforme entre vértices y medios"
      tol_rel: 1.0e-12
    - name: face_traction_balance
      setup: "Hex20 cubo unitario con tracción (p, 0, 0) sobre cara +ξ (face 3)"
      expect: "Σf_x = p·A_cara = p·1; vértices reciben −p/12 cada uno y medios
              4p/12 cada uno (suma = 4·(−1/12) + 4·(4/12) = 12/12 = 1, OK)"
      tol_rel: 1.0e-12
    - name: jacobiano_degenerado_abortado
      setup: "Hex20 con un nodo medio desplazado fuera del rango cuasi-recto"
      expect: "ValueError en _compute_kinematics_hex20 cuando det J ≤ tol"
    - name: convergencia_cantilever_macneal
      setup: "cantilever esbelto L/h = 30 con 1×1 sección + N elementos a lo
             largo; flecha vs Euler-Bernoulli"
      expect: "Hex20 con 1×1×6 alcanza >0.9 · u_EB (Hex8 con 1×1×6 solo ≈ 0.5);
              convergencia h monotónica al refinar"
      tol_rel: "documentar — sin tolerancia estricta"
    - name: locking_volumetrico_atenuado
      setup: "cantilever esbelto con Hex20; ν ∈ {0.3, 0.4999}"
      expect: "ratio u_tip(ν=0.4999) / u_tip(ν=0.3) significativamente mayor
              que el del Hex8 (Hex20 < 0.6 → Hex20 > 0.8 esperado al menos)"
      tol_rel: "documentar valores; cota relativa Hex20 vs Hex8"
    - name: mass_consistent_total
      setup: "Hex20 cubo unitario con ρ = 7850, masa consistente"
      expect: "Σ M (todas las entradas) = 3·ρ·V = 3·7850·1 (conservación 3D)"
      tol_rel: 1.0e-9
    - name: mass_lumped_total_HRZ
      setup: "Hex20 cubo unitario con ρ = 7850, masa lumped HRZ"
      expect: "diag(M) suma = 3·ρ·V; matriz estrictamente diagonal; todas las
              entradas estrictamente positivas (HRZ canónico, no row-sum)"
      tol_abs: 1.0e-12

references:
  - "Bathe K.J. (2014). Finite Element Procedures. §5.3.2 (serendipity 3D)."
  - "Cook R.D., Malkus D.S., Plesha M.E., Witt R.J. (2002). Concepts and
    Applications of FEA. §6.5, §13.6 (locking volumétrico Hex20)."
  - "Zienkiewicz O.C., Taylor R.L. (2005). The Finite Element Method, vol. 1,
    §8.7 (familia serendípita 3D)."
  - "ADR 0012 — Sólidos 3D: convención Voigt 6D y caras."
```

---

## Implementación

- **Archivo**: [solidum/elements/solid_3d/hex20.py](../../solidum/elements/solid_3d/hex20.py).
- **Clase**: `Hex20` (registrada vía `@ElementRegistry.register`); subclase de [`_HigherOrderSolid3D`](../../solidum/elements/solid_3d/_shared.py) desde **sub-fase 2** de A.ter — la entrada del `Hex27` disparó la centralización (regla de los dos casos reales). El cuerpo de `hex20.py` quedó como pura declaración de atributos: funciones de forma, cuadraturas, FACE_NODES y funciones de cara Quad8 (`_N_quad8`/`_dN_quad8` del paquete 2D).
- **Funciones núcleo**: `_N_hex20(xi, eta, zeta)`, `_dN_hex20(xi, eta, zeta)` en [solidum/elements/solid_3d/_shared.py](../../solidum/elements/solid_3d/_shared.py). Numpy puro (no `@njit`) por la complejidad de la rama por nodo serendipity, paritario con `_N_quad8`/`_dN_quad8` del 2D. El kinematics genérico ``_kinematics_higher_order_3d`` (también en `_shared.py`) es reutilizado por la base con cualquier par ``(shape_fn, grad_fn)`` 3D.
- **Tests**:
  - `tests/test_solid_3d_higher_order.py` — clase `TestHex20Element` (dimensiones/DOFs, patch lineal, patch cuadrático, simetría K, jacobiano degenerado abortado, body load suma, face traction suma + reparto vértice/medio, masa consistente y lumped HRZ, gauss state).
  - `tests/test_rigid_body_modes.py` — `TestRBMHex20` (3 trans + 3 rot + rank-deficiency = 6).
  - `tests/validation/test_cube_lame_3d.py` (extender) — `test_uniaxial_traction_hex20_exact`.
  - `tests/validation/test_macneal_beam_3d.py` (extender) — convergencia Hex20 1×1 sección vs Hex8 1×1.
  - `tests/test_volumetric_locking_3d.py` (extender) — `TestHex20VolumetricLocking` documentando mitigación parcial.

---

## Diálogo

- **2026-05-27** · Spec inicial. Sub-fase 1 de la sub-etapa A.ter (sólidos 3D cuadráticos). Hex20 standalone — sin base interna compartida en esta sub-fase; el segundo caso (Hex27) entra en sub-fase 2 y dispara la centralización en `_HigherOrderSolid3D` (regla de los dos casos reales aplicada).
- **2026-05-27** · Implementado y validado. Patch test lineal y cuadrático exactos a precisión máquina, 6 modos rígidos exactos, cubo Lamé 3D uniaxial e hidrostático exactos. Demostración cuantitativa de la mitigación del shear locking sobre MacNeal 3D: Hex20 6×1×1 alcanza 97% de la flecha analítica EB, frente a < 55% de Hex8 12×1×1. Mitigación parcial del locking volumétrico documentada: ratio u(ν=0.4999)/u(ν=0.3) ≈ 0.80 (vs < 0.6 del Hex8). Conteo de modos espurios con cuadratura reducida `hex_2x2x2`: 12 = 6 rígidos + 6 hourglass por elemento aislado (corregido desde el "4" inicial tras verificación numérica directa; coincide con Belytschko-Liu-Moran §8.4.4). Spec actualizada a `status: validated`.
- **2026-05-27** · Refactor a la base `_HigherOrderSolid3D` con la entrada del `Hex27` (sub-fase 2). El cuerpo del archivo `hex20.py` queda como declaración de atributos (`_SHAPE_FN`, `_GRAD_FN`, `_DEFAULT_QUADRATURE`, `_MASS_QUADRATURE`, `_FACE_N_FN`, `_FACE_DN_FN`, `_FACE_QUADRATURE`, `FACE_NODES`). La maquinaria de cálculo (bucle de Gauss para K, gauss state, body load, face traction, masa consistente, HRZ) ahora vive en la base — paritaria con la centralización 2D `_HigherOrderSolid2D` (Quad8/Quad9/Tri6). Los 26 tests previos pasan sin modificaciones.
