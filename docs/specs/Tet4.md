# Tet4 — tetraedro lineal de deformación constante 3D (CST 3D)

> Elemento sólido 3D de primer orden, espejo natural de `Tri3`. Pieza de la Etapa 7. Convención Voigt 6D y numeración de caras por ADR 0012.

---

## Especificación física

### 0. Descripción general

Elemento sólido tridimensional de primer orden, tetraedro de cuatro nodos. Aproximación **lineal** $C^0$ del campo de desplazamientos en coordenadas baricéntricas. Análogo 3D del `Tri3`: la matriz $\mathbf B$ es **constante sobre el elemento**, por lo que $\boldsymbol\varepsilon$ y $\boldsymbol\sigma$ son uniformes dentro de cada `Tet4` (*constant strain tetrahedron* en la literatura). Régimen de pequeños desplazamientos y deformaciones.

### 1. Cinemática de desplazamientos

$$u_x = \sum_{i=1}^{4} N_i\, u_x^{(i)}, \qquad u_y = \sum N_i\, u_y^{(i)}, \qquad u_z = \sum N_i\, u_z^{(i)}.$$

### 2. Cinemática de deformaciones

Notación Voigt 6D del proyecto (ADR 0012):

$$\boldsymbol\varepsilon = [\varepsilon_{xx},\ \varepsilon_{yy},\ \varepsilon_{zz},\ \gamma_{xy},\ \gamma_{yz},\ \gamma_{xz}]^\top$$

Por la linealidad de $N_i$, las derivadas $\partial N_i / \partial \mathbf x$ son constantes y $\boldsymbol\varepsilon$ es uniforme en el elemento. Idéntico patrón conceptual al `Tri3` en 2D.

### 3. Ecuación constitutiva

Delegada al material 3D (`STRAIN_DIM = 6`). En esta etapa solo `Elastic3D` está disponible.

### 4. Equilibrio — forma fuerte

Idéntica a la del `Hex8`. El principio variacional no depende de la forma del elemento.

### 5. Equilibrio — forma débil

Idéntica a la del `Hex8`.

---

## Formulación numérica (FEM)

### 6. Discretización

Cuatro nodos en orden estándar de tetraedro VTK. En el tetraedro de referencia con vértices:

| Nodo local | $(\xi_i, \eta_i, \zeta_i)$ |
|---|---|
| 0 | $(0, 0, 0)$ |
| 1 | $(1, 0, 0)$ |
| 2 | $(0, 1, 0)$ |
| 3 | $(0, 0, 1)$ |

El orden de nodos asegura que el vector $(\mathbf x_1 - \mathbf x_0) \times (\mathbf x_2 - \mathbf x_0)$ apunte hacia $\mathbf x_3$ (volumen orientado positivo). Mapeo isoparamétrico desde el tetraedro de referencia; jacobiano $\mathbf J$ constante en todo el elemento; $\det\mathbf J = 6\,V_e$ con $V_e$ el volumen del tetraedro. El código aborta con `ValueError` si $\det\mathbf J \le \text{tol}$ (nodos coplanares o en orden invertido).

### 7. Funciones de forma

Coordenadas baricéntricas / tetraédricas:

$$N_1 = 1 - \xi - \eta - \zeta, \qquad N_2 = \xi, \qquad N_3 = \eta, \qquad N_4 = \zeta.$$

Sus derivadas en coordenadas naturales son constantes:

$$\partial N_i / \partial \xi = (-1, 1, 0, 0),\quad \partial N_i / \partial \eta = (-1, 0, 1, 0),\quad \partial N_i / \partial \zeta = (-1, 0, 0, 1).$$

### 8. Matriz deformación-desplazamiento

$\mathbf B$ tiene tamaño $6 \times 12$ y, una vez transformadas las derivadas a globales vía $\mathbf J^{-1}$, sus entradas son **constantes** sobre el elemento:

$$\mathbf B_i = \begin{bmatrix}
\partial N_i / \partial x & 0 & 0 \\
0 & \partial N_i / \partial y & 0 \\
0 & 0 & \partial N_i / \partial z \\
\partial N_i / \partial y & \partial N_i / \partial x & 0 \\
0 & \partial N_i / \partial z & \partial N_i / \partial y \\
\partial N_i / \partial z & 0 & \partial N_i / \partial x
\end{bmatrix}, \qquad \boldsymbol\varepsilon = \mathbf B\,\mathbf u_e.$$

Mismo patrón filas/columnas que `Hex8` (ver §8 de su spec).

### 9. Rigidez elemental

$$\mathbf K_e = \mathbf B^\top\,\mathbf D\,\mathbf B\, V_e,$$

con $V_e = \tfrac{1}{6}\,\det\mathbf J$ el volumen del tetraedro. La integración es **exacta con un único punto central** porque $\mathbf B$ y $\boldsymbol\sigma$ son constantes.

### 10. Fuerzas internas

$$\mathbf F_{\text{int}} = \mathbf B^\top\,\boldsymbol\sigma\, V_e.$$

### 11. Cuadratura

Un punto central (baricentro del tetraedro) con peso $w = 1/6$ sobre el tetraedro de referencia. Es exacto para formas constantes de $\mathbf B$ y $\boldsymbol\sigma$.

### 12. Cargas distribuidas consistentes

**Fuerza de cuerpo $\mathbf b$**. Con $N_i$ lineales y $\mathbf b$ uniforme la integral es exacta:

$$\mathbf f_e^{b} = \int_{\Omega_e} \mathbf N^\top\,\mathbf b\, dV \;\Longrightarrow\; \text{cada nodo recibe } \tfrac{1}{4}\,\mathbf b\,V_e.$$

**Tracción de superficie $\bar{\mathbf t}$** sobre una cara triangular. Caras numeradas con normal saliente (ADR 0012):

| Cara | Nodos | Nodo opuesto |
|------|-------|--------------|
| 0 | (1, 2, 3) | 0 |
| 1 | (0, 3, 2) | 1 |
| 2 | (0, 1, 3) | 2 |
| 3 | (0, 2, 1) | 3 |

Sobre una cara triangular con tracción constante el reparto es exacto:

$$\mathbf f_e^{t} = \tfrac{1}{3}\,\bar{\mathbf t}\,A_\text{cara} \text{ por nodo del triángulo}, \qquad 0 \text{ al nodo opuesto}.$$

Cuadratura 1 punto centroide de la cara (suficiente por linealidad de las funciones de forma sobre la cara). $\bar{\mathbf t}$ se especifica en coordenadas globales $(t_x, t_y, t_z)$. Tracción variable o presión normal no soportadas en este paso.

### 13. Salida por punto de Gauss

`compute_gauss_state(U)` devuelve la misma estructura que `Hex8` pero con un único punto (centroide). Para `Tet4` $\boldsymbol\varepsilon$ y $\boldsymbol\sigma$ son constantes sobre el elemento, así que el dato del punto coincide con el promedio del elemento. Por ADR 0012, `internal_forces` devuelve `None`.

---

## Limitación crítica — *shear locking*

El `Tet4` es notoriamente rígido en flexión y dominado por cortante, exactamente como su análogo `Tri3` en 2D: cuatro DOFs de desplazamiento por dirección no bastan para representar el modo de flexión puro de una viga discretizada en tetraedros, por lo que el elemento responde con cortante espurio y subestima sistemáticamente la deflexión. La severidad crece con la relación $L/h$ del problema.

**Recomendación operativa**: usar `Hex8` por defecto en mallas hexaédricas o mixtas; reservar `Tet4` para zonas de transición de malla donde la generación automática de hexaedros falla y se requiere malla tetraédrica de baja regularidad. Para problemas dominados por flexión se necesitan tetraedros cuadráticos `Tet10` (sub-etapa posterior). Locking volumétrico con $\nu \to 0.5$ es **aún peor** que en `Hex8` por la misma razón: el espacio lineal es muy pobre.

---

## Contrato de implementación

```yaml
name: Tet4
kind: element
status: validated

interface:
  dof_names: [ux, uy, uz]
  n_nodes: 4
  strain_dim: 6                 # Voigt 3D [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]
  n_integration_points: 1

parameters: []

material_contract:
  signature: "compute_state(ε, state) -> (σ, C_tangent, state')"
  strain_kind: "Voigt 3D [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz] con γ_ij = 2·ε_ij"

conventions:
  sign: "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "[xx, yy, zz, xy, yz, xz] con γ_ij = 2·ε_ij (ADR 0012)"
  node_orientation: "convención VTK_TETRA; nodos 0-2 base con normal hacia 3 (volumen orientado positivo); asegura det J > 0"
  face_numbering: "ADR 0012 — 4 caras con normal saliente; cara i opuesta al nodo i"

validity:
  - "pequeños desplazamientos y deformaciones"
  - "campos de tensión que varían suavemente sobre el elemento"
  - "compatible con materiales con STRAIN_DIM = 6 (en esta etapa: Elastic3D)"

out_of_scope:
  - "flexión dominante (shear locking severo)"
  - "casi-incompresibilidad (ν → 0.5 — locking volumétrico aún peor que Hex8)"
  - "elementos de alto orden (Tet10 — sub-etapa posterior)"
  - "grandes deformaciones"
  - "internal_forces (ADR 0012: sólidos exponen compute_gauss_state)"

acceptance:
  verification:
    - name: patch_test_3d
      setup: "malla tetraédrica con al menos un nodo interior (e.g.
             división de un cubo en 5 tetraedros con un nodo central
             libre); campo lineal u_i = a_ij·x_j impuesto en el contorno"
      expect: "nodo interior adopta el campo lineal exacto y ε constante
              e igual a 0.5·(a_ij + a_ji) en cada elemento; σ uniforme"
      tol_rel: 1.0e-10
    - name: modos_de_cuerpo_rigido
      setup: "evaluar K_e sobre un Tet4 cualquiera y proyectar sobre los 6
             modos de cuerpo rígido"
      expect: "K_e · u_rigid ≈ 0 en cada uno de los 6 modos"
      tol_abs: 1.0e-9
    - name: simetria_K_e
      setup: "evaluar K_e con Elastic3D"
      expect: "K_e == K_e.T exacto"
      tol_abs: 1.0e-12
  specific:
    - name: parche_deformacion_constante
      setup: "malla tetraédrica regular con campo lineal y Elastic3D"
      expect: "ε constante (CST 3D reproduce campos lineales por construcción)"
      tol_rel: 1.0e-12
    - name: jacobiano_degenerado_abortado
      setup: "cuatro nodos coplanares"
      expect: "ValueError en _compute_kinematics_tet4"
    - name: body_load_uniforme
      setup: "Tet4 con b uniforme"
      expect: "cada nodo recibe (1/4)·b·V_e; Σ = b·V_e"
      tol_rel: 1.0e-12
    - name: traccion_uniforme_cara
      setup: "Tet4 con tracción constante en una cara"
      expect: "Σf = t̄·A_cara; reparto A_cara/3·t̄ a cada uno de los 3
              nodos de la cara, 0 al nodo opuesto"
      tol_rel: 1.0e-12
    - name: cubo_lame_3d_malla_tetraedrica
      setup: "cubo unitario discretizado con tetraedros (e.g. 5 o 6 Tet4
             por cubo); tracción uniaxial σ_xx = p en cara +ξ"
      expect: "u_x(L) ≈ p·L/E con error decreciente al refinar; orden
              de convergencia O(h) en norma energética"
      tol_rel: "documentar convergencia, sin tolerancia absoluta"

references:
  - "Cook R.D., Malkus D.S., Plesha M.E., Witt R.J. (2002). Concepts and Applications of FEA. §6.8 (CST 3D)."
  - "Zienkiewicz O.C., The Finite Element Method, vol. 1, §6.7 (limitations of linear tets)."
  - "Bathe K.J. (2014). Finite Element Procedures. §6.6.2."
  - "ADR 0012 — Sólidos 3D: convención Voigt 6D y caras."
```

---

## Implementación

- **Archivo**: [fenix/elements/solid_3d/tet4.py](../../fenix/elements/solid_3d/tet4.py).
- **Clase**: `Tet4` (registrada vía `@ElementRegistry.register`).
- **Función núcleo**: `_compute_kinematics_tet4(coords)` con `@njit` en [fenix/elements/solid_3d/_shared.py](../../fenix/elements/solid_3d/_shared.py). Reutiliza `_compute_integrands_3d` (también en `_shared`, compartido con Hex8) con peso `1/6` (volumen del tet de referencia).
- **Masa consistente analítica**: ``M_ij = ρ·V_e·(1 + δ_ij)/20`` (matriz escalar 4×4, expandida 3D vía `_expand_scalar_mass_3d`). Evita cuadratura insuficiente en el centroide.
- **Tests**:
  - [tests/test_solid_3d.py](../../tests/test_solid_3d.py) — clase `TestTet4Element` (10 tests: dimensiones/DOFs, simetría exacta de K (no numérica, por fórmula cerrada), volumen, patch tracción uniaxial, jacobiano degenerado abortado, body load, face traction (cara opuesta a nodo recibe cero, reparto A/3 a los 3 nodos restantes), masa consistente y lumped HRZ (reparto equitativo ρV/4 por DOF)).
  - [tests/test_rigid_body_modes.py](../../tests/test_rigid_body_modes.py) — `TestRBMTet4` (5 tests: traslación, 3 rotaciones, rank-deficiency = 6).
  - [tests/validation/test_cube_lame_3d.py](../../tests/validation/test_cube_lame_3d.py) — `test_uniaxial_traction_tet4_mesh_exact` (cubo dividido en 5 tetraedros con tracción uniaxial; solución exacta nodal).

---

## Diálogo

- **2026-05-19** · Spec inicial. Elemento espejo natural de `Tri3`; convención Voigt 6D y numeración de caras fijadas por ADR 0012. Sin elementos de alto orden en esta etapa. Shear locking y locking volumétrico declarados como limitaciones; `Hex8` es la opción preferente cuando la malla lo permite.
