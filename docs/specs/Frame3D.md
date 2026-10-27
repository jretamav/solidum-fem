# ELEMENTO MARCO/VIGA 3D EULER-BERNOULLI

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Elemento 1D inmerso en el espacio tridimensional que modela una **viga esbelta 3D**. Dos nodos rígidamente conectados, 6 DOFs por nodo ($u_x, u_y, u_z, r_x, r_y, r_z$). Captura seis acciones internas:

- **Axial** (fuerza normal $N$ a lo largo del eje local $x$).
- **Cortante en dos planos** ($V_y$ en el plano $xy$, $V_z$ en el plano $xz$).
- **Flexión en dos planos** ($M_z$ alrededor del eje local $z$, $M_y$ alrededor del eje local $y$).
- **Torsión** ($M_x$ alrededor del eje longitudinal).

Régimen de **linealidad geométrica** (pequeños desplazamientos y rotaciones). Hipótesis de Euler-Bernoulli: secciones planas perpendiculares al eje deformado (sin cizalladura transversal).

### 1. Hipótesis
- Secciones planas y perpendiculares al eje neutro deformado (igual que Frame2DEuler).
- Flexiones en los dos planos **desacopladas** entre sí (ejes principales de inercia).
- Torsión de Saint-Venant pura (sin alabeo). Válida para secciones macizas o cerradas; abiertas de pared delgada requerirían formulación con alabeo (fuera de alcance).
- Camino con material: material isótropo 1D elástico lineal; $G = E/[2(1+\nu)]$. Camino con sección: la torsión es elástica y desacoplada, con el $GJ$ que declara la sección ([FiberSection](FiberSection.md) §1).

### 2. Sistema de ejes locales
Tres ejes ortonormales en cada elemento:
- **$\hat{\mathbf x}$ local**: dirección del eje de la barra ($\hat{\mathbf x} = (\mathbf x_2 - \mathbf x_1)/L$).
- **$\hat{\mathbf y}$ local** y **$\hat{\mathbf z}$ local**: ejes principales de inercia de la sección, ortogonales al eje.

La orientación de $\hat{\mathbf y}, \hat{\mathbf z}$ no está determinada solo por la dirección del eje — hace falta un **vector de referencia** proporcionado por el usuario (o un default). Es la diferencia estructural respecto al caso 2D, donde el plano era único.

### 3. Medidas de deformación / fuerzas internas
- Axial: $\varepsilon_{axial} = du/ds$, fuerza axial $N = EA\,\varepsilon_{axial}$ (o $N = A\,\sigma$ con $\sigma$ del material).
- Flexión en plano $xy$: curvatura $\kappa_z = d^2v/ds^2$, momento $M_z = EI_z\,\kappa_z$.
- Flexión en plano $xz$: curvatura $\kappa_y = -d^2w/ds^2$ (con $w' = -\theta_y$ por la regla de la mano derecha), momento $M_y = EI_y\,\kappa_y$.
- Torsión: $d\theta_x/ds$, momento torsor $M_x = GJ\,d\theta_x/ds$.

**Constitutiva: dos caminos excluyentes** ([ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md)):

- **Con `material` + `A`, `Iy`, `Iz`, `J`** (+ `nu`): las leyes anteriores en forma cerrada. Sólo materiales **elásticos lineales** (`IS_LINEAR_ELASTIC = True`); un material inelástico se rechaza al construir con un mensaje que remite a la sección por fibras.
- **Con `section`** (+ `n_ip`): $N$, $M_y$, $M_z$ los da la ley de sección, que integra el esfuerzo $\sigma_{xx}$ de cada fibra con su material 1D ([FiberSection](FiberSection.md) §2-§4); la torsión es $T = GJ\,\theta_x'$ con el `GJ` de la sección, obligatorio.

Hasta el 2026-09-29 el camino con material admitía materiales inelásticos y escalaba toda la matriz local con la $E_t$ axial, torsión incluida ($G = E_t/[2(1+\nu)]$); esa ley desapareció (ver *Diálogo*).

### 4. Equilibrio — forma débil
PTV. El elemento calcula solo el lado interno:
$$\delta W_{\text{int}} = \int_0^L \!\left(N\,\delta\varepsilon + M_z\,\delta\kappa_z + M_y\,\delta\kappa_y + M_x\,\delta\theta_x'\right)\,ds.$$

---

## Formulación numérica (FEM)

### 5. Discretización
Dos nodos, 6 DOFs por nodo:
$$\mathbf u_e = [u_x^{(1)}, u_y^{(1)}, u_z^{(1)}, r_x^{(1)}, r_y^{(1)}, r_z^{(1)},\; u_x^{(2)}, u_y^{(2)}, u_z^{(2)}, r_x^{(2)}, r_y^{(2)}, r_z^{(2)}]^\top.$$

Vector de 12 componentes → matrices $12\times 12$.

### 6. Construcción de ejes locales
Dados $\mathbf X_1, \mathbf X_2$ y un vector de referencia $\mathbf v_{\text{ref}}$ proporcionado:

1. $\hat{\mathbf x} = (\mathbf X_2 - \mathbf X_1)/L$.
2. Proyectar $\mathbf v_{\text{ref}}$ al plano perpendicular a $\hat{\mathbf x}$:
   $$\tilde{\mathbf y} = \mathbf v_{\text{ref}} - (\mathbf v_{\text{ref}} \cdot \hat{\mathbf x})\,\hat{\mathbf x}.$$
3. $\hat{\mathbf y} = \tilde{\mathbf y}/\|\tilde{\mathbf y}\|$ (falla si $\|\tilde{\mathbf y}\| \approx 0$: $\mathbf v_{\text{ref}}$ paralelo al eje).
4. $\hat{\mathbf z} = \hat{\mathbf x} \times \hat{\mathbf y}$.

**Default de $\mathbf v_{\text{ref}}$** (cuando el usuario no lo proporciona): $[0, 0, 1]$ (eje $z$ global como "vertical"). Si el eje del elemento es cercanamente paralelo a $z$ global (viga vertical), se usa $[1, 0, 0]$ como fallback y se emite una advertencia.

### 7. Matriz de transformación global → local
$$\boldsymbol{\lambda} = \begin{bmatrix} \hat{\mathbf x}^\top \\ \hat{\mathbf y}^\top \\ \hat{\mathbf z}^\top \end{bmatrix} \quad (3\times 3,\ \text{ortogonal}).$$

$$\mathbf T = \text{bloques diag}(\boldsymbol{\lambda}, \boldsymbol{\lambda}, \boldsymbol{\lambda}, \boldsymbol{\lambda}) \quad (12\times 12).$$

Aplicada a los 4 bloques de 3 DOFs: desplazamientos nodo 1, rotaciones nodo 1, desplazamientos nodo 2, rotaciones nodo 2. Las rotaciones infinitesimales en 3D transforman como vectores.

$$\mathbf u_{\text{local}} = \mathbf T\,\mathbf u_e.$$

### 8. Matriz de rigidez local
Desacoplada en cuatro contribuciones: axial, torsión, flexión $xy$ (plano con $u_y$ y $r_z$, inercia $I_z$), flexión $xz$ (plano con $u_z$ y $r_y$, inercia $I_y$).

Definiciones:
$$a_z = \tfrac{12EI_z}{L^3},\ b_z = \tfrac{6EI_z}{L^2},\ c_z = \tfrac{4EI_z}{L},\ d_z = \tfrac{2EI_z}{L},$$
$$a_y = \tfrac{12EI_y}{L^3},\ b_y = \tfrac{6EI_y}{L^2},\ c_y = \tfrac{4EI_y}{L},\ d_y = \tfrac{2EI_y}{L}.$$

La rigidez $\mathbf K_{\text{local}}$ es simétrica, con los bloques desacoplados (0s fuera de su plano correspondiente). Los signos del bloque de flexión $xz$ se invierten respecto al $xy$ porque $u_z$ y $r_y$ forman una pareja dextrógira opuesta (consecuencia del producto vectorial $\hat{\mathbf x}\times\hat{\mathbf y}=\hat{\mathbf z}$).

Forma concreta en términos de los 12 DOFs locales (la matriz completa vive en la implementación).

### 9. Fuerzas internas
$$\mathbf F_{\text{int,local}} = \mathbf K_{\text{local}}\,\mathbf u_{\text{local}},$$
con la componente axial escrita con el $\sigma$ que devuelve el material, $F_1 = -\sigma A$, $F_7 = +\sigma A$; con material elástico lineal (el único admitido en este camino) coincide con la fila axial del producto.

Camino con sección: con $\mathbf e_g = [\varepsilon_0,\ \kappa_y,\ \kappa_z]^\top(x_g) = \mathbf B_g\,\mathbf u_{\text{local}}$, donde la fila de $\kappa_z$ es $[N_1'', N_2'', N_3'', N_4'']$ sobre $(v_1, \theta_{z1}, v_2, \theta_{z2})$ y la de $\kappa_y$ es $[-N_1'', N_2'', -N_3'', N_4'']$ sobre $(w_1, \theta_{y1}, w_2, \theta_{y2})$ ([FiberSection](FiberSection.md) §5),
$$\mathbf F_{\text{int,local}} = \sum_g w_g\,\tfrac{L}{2}\,\mathbf B_g^\top\,\mathbf s_g, \qquad \mathbf K_{\text{local}} = \sum_g w_g\,\tfrac{L}{2}\,\mathbf B_g^\top\,\mathbf k_{s,g}\,\mathbf B_g,$$
más la torsión elástica: $GJ/L$ en el bloque $(\theta_{x1}, \theta_{x2})$ de $\mathbf K_{\text{local}}$ y $\mp T$ con $T = GJ\,(\theta_{x2} - \theta_{x1})/L$ en $\mathbf F_{\text{int,local}}$. Con fibras elásticas y `n_ip` ≥ 2 reproduce la matriz de §8 exactamente.

### 10. Ensamblaje global
$$\mathbf K_{\text{global}} = \mathbf T^\top\,\mathbf K_{\text{local}}\,\mathbf T, \qquad \mathbf F_{\text{int}} = \mathbf T^\top\,\mathbf F_{\text{int,local}}.$$

### 11. Cuadratura
Camino con material: no aplica (forma cerrada; integración analítica de los polinomios de Hermite). Camino con sección: Gauss-Legendre de `n_ip` puntos a lo largo del elemento (3 por defecto, mínimo 2), con `N_INTEGRATION_POINTS = n_ip` como atributo de instancia.

### 12. Matriz de masa
Cuatro bloques desacoplados en ejes locales, como la rigidez. La flexión **no lleva inercia rotatoria** de la sección (hipótesis de Euler-Bernoulli); la torsión sí lleva la inercia del giro de la sección alrededor del eje, $\rho J_p$ con $J_p = I_y + I_z$ (momento polar geométrico, no la constante de Saint-Venant $J$, que rige la rigidez). Con sección, $\rho A = \sum_i\rho_iA_i$ y $\rho J_p = \rho I_y + \rho I_z$ ([FiberSection](FiberSection.md) §6).

**Consistente** (`lumping="consistent"`, por defecto):

- Axial $(u_1, u_2)$: $\tfrac{\rho A L}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$.
- Torsión $(\theta_{x1}, \theta_{x2})$: $\tfrac{\rho J_p L}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$.
- Flexión en el plano $xy$, $(v_1, \theta_{z1}, v_2, \theta_{z2})$ con $v' = \theta_z$: la hermítica de [Frame2DEuler](Frame2DEuler.md) §13,
$$\mathbf M_{xy} = \frac{\rho A L}{420}\begin{bmatrix}
156 & 22L & 54 & -13L \\
22L & 4L^2 & 13L & -3L^2 \\
54 & 13L & 156 & -22L \\
-13L & -3L^2 & -22L & 4L^2
\end{bmatrix}.$$
- Flexión en el plano $xz$, $(w_1, \theta_{y1}, w_2, \theta_{y2})$: por $w' = -\theta_y$, los términos de acoplamiento traslación-giro cambian de signo, $\mathbf M_{xz} = \mathbf S\,\mathbf M_{xy}\,\mathbf S$ con $\mathbf S = \operatorname{diag}(1, -1, 1, -1)$:
$$\mathbf M_{xz} = \frac{\rho A L}{420}\begin{bmatrix}
156 & -22L & 54 & 13L \\
-22L & 4L^2 & -13L & -3L^2 \\
54 & -13L & 156 & 22L \\
13L & -3L^2 & 22L & 4L^2
\end{bmatrix}.$$

Reproduce exactamente la energía cinética de un giro rígido de velocidad angular $\boldsymbol\omega$ alrededor de un extremo: $\tfrac12\,\rho A\,\tfrac{L^3}{3}\,|\boldsymbol\omega_\perp|^2 + \tfrac12\,\rho J_p L\,\omega_\parallel^2$, con $\omega_\parallel$ la componente sobre el eje de la barra.

**Concentrada** (`lumping="lumped"`): por nodo, $\rho A L/2$ en $u_x, u_y, u_z$, $\rho J_p L/2$ en $\theta_x$ y la HRZ de la flexión, $\rho A L^3/78$, en $\theta_y$ y $\theta_z$. Diagonal en locales; en globales sólo lo es si el eje del elemento coincide con un eje global: con eje oblicuo el bloque rotacional $\operatorname{diag}(\rho J_p L/2,\ \rho AL^3/78,\ \rho AL^3/78)$ no es escalar y la matriz queda diagonal por bloques por nodo (`CentralDifferenceSolver` la rechaza).

En globales, $\mathbf M_{\text{global}} = \mathbf T^\top\,\mathbf M_{\text{local}}\,\mathbf T$. Hasta el 2026-09-29 la flexión sumaba la inercia rotatoria $\tfrac{\rho I L}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$ en $\theta_y$ ($I_y$) y $\theta_z$ ($I_z$), y la concentrada usaba $\rho I_y L/2$ y $\rho I_z L/2$ (ver *Diálogo*); la torsión no cambia.

---

## Contrato de implementación

```yaml
name: Frame3D
kind: element
status: validated

interface:
  dof_names: [ux, uy, uz, rx, ry, rz]
  n_nodes: 2
  strain_dim: 1
  n_integration_points: 1   # camino con material; con section, N_INTEGRATION_POINTS = n_ip (atributo de instancia)

admission:                # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - bathe: "cap. 5 (elementos de viga)"
    - zienkiewicz_taylor: "vol. 2 (6.ª ed.), cap. 10: barras (formulación base, material elástico)"
  general_codes:
    - opensees: "elasticBeamColumn (3D)"
    - abaqus: "B33"

parameters:
  - { name: A,  type: float, required: false, desc: "Área de la sección transversal. Obligatoria con material; no se declara con section" }
  - { name: Iy, type: float, required: false, desc: "Momento de inercia alrededor del eje local y. Obligatorio con material; no se declara con section" }
  - { name: Iz, type: float, required: false, desc: "Momento de inercia alrededor del eje local z. Obligatorio con material; no se declara con section" }
  - { name: J,  type: float, required: false, desc: "Constante torsional de Saint-Venant. Obligatoria con material; con section la torsión usa el GJ de la sección" }
  - { name: nu, type: float, required: false, default: 0.3, desc: "Coeficiente de Poisson (si el material no lo expone). Sólo con material" }
  - { name: ref_vector, type: list, required: false, default: "[0, 0, 1]",
      desc: "Vector 3D de referencia que fija la orientación de los ejes locales y, z de la sección. Default: eje z global." }
  - { name: section, type: "id de sections", required: false, desc: "Sección de marco (FiberSection) con GJ declarado. Excluyente con material, A, Iy, Iz, J y nu (ADR 0021)" }
  - { name: n_ip, type: int, required: false, default: 3, desc: "Puntos de Gauss-Legendre a lo largo del elemento, ≥ 2. Sólo con section" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: "axial escalar (ε = (u_2 − u_1)/L en sistema local)"
  nonlinearity_model: "camino con material: sólo materiales elásticos lineales (IS_LINEAR_ELASTIC = True; un material inelástico se rechaza al construir). La no linealidad material entra por section (FiberSection): cada fibra evalúa su material 1D con ε = ε₀ − y·κ_z + z·κ_y; torsión elástica desacoplada con el GJ de la sección"
  poisson_source: "material.nu si existe; en su defecto parámetro `nu` del elemento (sólo camino con material)"

conventions:
  sign: "σ > 0 ⇔ ε > 0 (elongación); momentos siguen regla de la mano derecha respecto a los ejes locales"
  local_axes: "x_local = eje de la barra; y_local, z_local = ejes principales de inercia definidos por ref_vector"
  node_orientation: "eje local del nodo 1 al nodo 2"
  configuration: "inicial fija — L, T, ejes locales no se actualizan"

validity:
  - "vigas esbeltas (L/h ≳ 10 en ambos planos de flexión)"
  - "pequeños desplazamientos y pequeñas rotaciones en el espacio"
  - "|ε_axial| ≲ 1e-2"
  - "sección con ejes principales de inercia bien definidos por ref_vector"

out_of_scope:
  - grandes rotaciones espaciales (requeriría Frame3DCorot; pendiente)
  - torsión con alabeo restringido (secciones abiertas de pared delgada)
  - material inelástico con material + A, Iy, Iz, J (se rechaza al construir; la plasticidad y el daño entran por section, FiberSection)
  - degradación de la torsión y su interacción con N, My, Mz en el camino con sección (torsión elástica desacoplada, FiberSection §7)
  - pandeo por bifurcación
  - deformación por cortante transversal significativa (Timoshenko 3D, fuera de alcance)
  - "torsor distribuido (la fuerza de cuerpo, centrada en el eje, no produce torsión)"

numerical_caveats:
  - "Si ref_vector es paralelo al eje de la barra, la proyección al plano perpendicular se anula y el constructor aborta con mensaje claro — el usuario debe elegir otro vector."
  - "Para barras verticales (eje casi paralelo a z global) con ref_vector default, el constructor usa [1, 0, 0] como fallback y emite advertencia por consola."

acceptance:
  - name: soluciones manufacturadas (ADR 0022)                   # 2026-10-27
    setup: "línea recta inclinada; solución exacta impuesta en un extremo y sus fuerzas internas como cargas nodales en el otro; carga de línea como fuerza de cuerpo b(x); línea orientada al azar; u*, v*, w* y θx* lineal (no hay torsor distribuido); n = 4, 8, 16, 32"
    expect: "con campos polinómicos (carga de grado ≤ 4), error nodal de redondeo con 1 y 3 elementos (exactitud nodal, Tong 1969); con campos suaves, orden ≥ teórico − 0.1 por componente: u 2, v y w 4, θy y θz 3 en L²; ε 1, κy y κz 2 en energía; torsión exacta"
    tol_abs: 0.1
  - name: respuesta axial pura
    setup: "voladizo alineado con eje x, carga F axial en extremo libre"
    expect: "u_x = F·L/(E·A); demás componentes nulas"
    tol_rel: 1.0e-10

  - name: flexión en plano xy (coherencia con Frame2D)
    setup: "voladizo alineado con eje x, carga P en dirección y, ref_vector = [0,0,1]"
    expect: "v_y = P·L³/(3·E·Iz), θ_z = P·L²/(2·E·Iz)"
    tol_rel: 1.0e-10

  - name: flexión en plano xz
    setup: "voladizo alineado con eje x, carga P en dirección z, ref_vector = [0,1,0]"
    expect: "v_z = P·L³/(3·E·Iy), θ_y = -P·L²/(2·E·Iy) (signo por mano derecha)"
    tol_rel: 1.0e-10

  - name: torsión pura
    setup: "voladizo con momento T aplicado alrededor del eje de la barra en el extremo libre"
    expect: "θ_x = T·L/(G·J); demás componentes nulas"
    tol_rel: 1.0e-10

  - name: simetría de K
    setup: "viga en orientación arbitraria con ref_vector genérico"
    expect: "K_global = K_globalᵀ"
    tol_abs: 1.0e-8

  - name: guarda de material elástico lineal            # añadido 2026-09-29 (ADR 0021)
    setup: "material Elastoplastic1D o IsotropicDamage1D con A, Iy, Iz, J"
    expect: "ValueError al construir, con un mensaje que nombra IS_LINEAR_ELASTIC y remite a FiberSection con 'section'; una sección sin GJ también se rechaza"

  - name: camino con sección                            # añadido 2026-09-29 (ADR 0021)
    setup: "section = FiberSection; casos de acceptance de docs/specs/FiberSection.md"
    expect: "reproducción elástica con n_ip = 2 (elemento oblicuo, sección rectangular y mixta), M–κ en los dos ejes, interacción N–M, sección circular a 45°, signos"

  - name: matrices de masa consistente y concentrada   # añadido 2026-09-29
    setup: "elemento oblicuo con ref_vector genérico, caminos con material y con sección; M_local = T·M·Tᵀ"
    expect: "consistente = axial + torsión ρJ_pL/6·[[2,1],[1,2]] (J_p = I_y + I_z) + hermítica en xy y, en xz, la hermítica con los acoplamientos traslación-giro de signo opuesto; sin inercia rotatoria de flexión. Concentrada = diag(ρAL/2 ×3, ρJ_pL/2, ρAL³/78, ρAL³/78) por nodo"
    tol_rel: 1.0e-12

  - name: energía cinética de movimientos rígidos        # añadido 2026-09-29
    setup: "traslación v; giro ω oblicuo alrededor de un extremo"
    expect: "½u̇ᵀMu̇ = ρAL|v|²/2 (consistente y concentrada); ½·ρA·L³/3·|ω⊥|² + ½·ρJ_p·L·ω∥² (consistente)"
    tol_rel: 1.0e-12

references:
  - "Przemieniecki J.S., Theory of Matrix Structural Analysis, Tabla 11.3"
  - "Cook R.D., Concepts and Applications of FEA, §2.8"
  - "Bathe K.-J., Finite Element Procedures, cap. 5"
  - "Hinton E., Rock T., Zienkiewicz O.C. (1976). A note on mass lumping and related processes in the finite element method. Earthquake Eng. Struct. Dyn. 4, 245-249 (masa concentrada HRZ)"
```

---

## Implementación

- **Archivo**: [solidum/elements/frame3d.py](../../solidum/elements/frame3d.py) — archivo propio, sin dependencia de ningún otro elemento.
- **Clase**: `Frame3D(FrameSectionMixin, Element)`. Métodos privados `_build_local_frame` (longitud, matriz $\boldsymbol{\lambda}$ 3×3 y $\mathbf T$ 12×12), `_build_local_stiffness` (matriz $\mathbf K_{\text{local}}$ 12×12 con bloques axial, torsión y flexión en ambos planos, camino con material), `_build_section_B` y `_section_element_state` (camino con sección), `_build_local_mass`, `_build_local_mass_lumped` y `_mass_properties` (masa, §12). Comparte con los marcos planos la constitutiva genérica de `FrameSectionMixin` ([solidum/elements/frame/_section.py](../../solidum/elements/frame/_section.py), ADR 0021) y los bloques de masa de la flexión de [solidum/elements/frame/_shared.py](../../solidum/elements/frame/_shared.py) (`hermite_consistent_mass`, `hrz_rotational_mass`). `REFERENCE_KWARGS = {"material": MaterialRegistry, "section": SectionRegistry}`.
- **Tests**: [`tests/test_line_elements.py`](../../tests/test_line_elements.py) (campo interior, energía igual a la de la rigidez, carga como campo) y [`tests/validation/test_mms_frames.py`](../../tests/validation/test_mms_frames.py) (soluciones manufacturadas); [tests/test_frame3d.py](../../tests/test_frame3d.py) · `TestFrame3DAcceptance` — cubre los 5 criterios + validación de `ref_vector` paralelo al eje + registro:
  - `test_acceptance_respuesta_axial_pura`
  - `test_acceptance_flexion_xy_coincide_con_frame2d`
  - `test_acceptance_flexion_xz`
  - `test_acceptance_torsion_pura`
  - `test_acceptance_simetria_K`
  - `test_rechazo_ref_vector_paralelo`
  - `test_registro_en_registry`
  - Criterio 6 (guarda): [tests/test_fiber_section.py](../../tests/test_fiber_section.py) · `TestGuardas` (`test_material_inelastico_remite_a_la_seccion`, `test_exactamente_un_camino`, `test_frame3d_exige_GJ`).
  - Criterio 7 (camino con sección): [tests/validation/test_fiber_section_validation.py](../../tests/validation/test_fiber_section_validation.py) (casos 1, 2, 3, 6 y 7 de la spec FiberSection) y, en `tests/test_fiber_section.py`, `TestTangente.test_frame3d_con_torsion`, `TestMarcosConSeccion.test_torsion_elastica_desacoplada`, `TestMasa`, `TestPesoPropio` y `TestYaml.test_frame3d_torsion_con_GJ_del_yaml`.
  - Criterios 8 y 9 (masa): [tests/test_frame_mass.py](../../tests/test_frame_mass.py) · `TestFrame3D` (consistente en los dos planos y torsión, signo del plano $xz$, HRZ, $I_y$ e $I_z$ sólo en la torsión), `TestCaminoConSeccion.test_frame3d_con_rho_A_y_rho_Jp_de_las_fibras`, `test_traslacion_rigida_3d`, `test_giro_rigido_3d_consistente` y `test_giro_rigido_por_plano_local_3d`. Diagonalidad por bloques en orientación oblicua: [tests/test_mass_lumping.py](../../tests/test_mass_lumping.py) · `TestLumpedFrame3D`.
- **Notas de traducción**:
  - La matriz $\boldsymbol{\lambda}$ 3×3 se construye como $[\hat{\mathbf x}^\top; \hat{\mathbf y}^\top; \hat{\mathbf z}^\top]$ — las filas son los ejes locales expresados en coordenadas globales. $\mathbf T$ 12×12 es un `blkdiag` de cuatro copias de $\boldsymbol{\lambda}$ (traslaciones y rotaciones transforman como vectores en 3D infinitesimal).
  - El default de `ref_vector` es `[0, 0, 1]`; el fallback `[1, 0, 0]` se activa cuando $|\hat{\mathbf x} \cdot \hat{\mathbf z}_{\text{global}}| > 0.99$ (barra cercanamente vertical). Ambas rutas emiten advertencia por consola la primera vez que se usan con default.
  - La matriz $\mathbf K_{\text{local}}$ se construye rellenando entradas no nulas en sus 4 bloques desacoplados (axial, torsión, flexión $xy$, flexión $xz$). Los signos de la flexión $xz$ están invertidos respecto a la flexión $xy$ porque la pareja $(u_z, r_y)$ tiene orientación dextrógira opuesta a $(u_y, r_z)$: este punto se verifica explícitamente en `test_acceptance_flexion_xz`.
  - `compute_internal_forces` devuelve cortantes, torsión y momentos en ambos extremos separados — útil para post-proceso estructural (diagramas de fuerzas internas). Con sección añade `section_deformations` y `section_forces` `(n_ip, 3)` = $[\varepsilon_0, \kappa_y, \kappa_z]$ y $[N, M_y, M_z]$ en los puntos de Gauss; como `internal_forces`, evalúa sobre una instantánea del trial y la restaura.
  - Camino con sección: sin `nu` (`self.nu = None`); la torsión la da `section.GJ`, y una sección sin `GJ` se rechaza al construir. Masa con $\rho A$ de la sección y $\rho J_p = \rho I_y + \rho I_z$ para la torsión, en las mismas funciones de masa consistente y concentrada; peso propio con $\rho A$ y fuerza de cuerpo con el área bruta $\sum A_i$.
  - Masa (§12): `_mass_properties` devuelve $(\rho A, \rho J_p)$; el plano $xz$ se obtiene de la hermítica del plano $xy$ con $\mathbf S = \operatorname{diag}(1, -1, 1, -1)$ (`_XZ_SIGNS`), el mismo cambio de signo que la rigidez.

---

## Diálogo

- **2026-04-21** · Elemento creado como componente totalmente autónomo: no hereda de `Frame2DEuler` ni de ninguna otra viga. La maquinaria cinemática 3D se implementa íntegra dentro de `frame3d.py`, coherente con la filosofía aplicada a cables.
- **2026-04-21** · Decisión sobre `ref_vector` vs "tercer nodo de orientación". Se optó por `ref_vector` (vector 3D proyectado al plano perpendicular al eje) por ser más compacto en YAML y no requerir nodos fantasma en la malla. El default `[0, 0, 1]` cubre la mayoría de casos de ingeniería estructural (estructuras terrestres con gravedad en $-z$); el fallback para barras verticales usa `[1, 0, 0]` con advertencia explícita.
- **2026-04-21** · Verificación cruzada con `Frame2DEuler`: el criterio 2 aplica una viga alineada con el eje $x$ global bajo carga en $y$; los desplazamientos resultantes deben coincidir con la solución 2D clásica $v = PL^3/(3EI_z)$, $\theta = PL^2/(2EI_z)$. El test lo verifica a 8 decimales. Esto blinda la coherencia dimensional: el 3D se reduce al 2D sin error acumulado cuando el problema es efectivamente plano.
- **2026-09-22** · Auditoría global: `internal_forces(U, equivalent_load=None)`: cuando el `Assembler` ensambló peso propio o fuerza de cuerpo, `solidum.run` pasa el vector nodal equivalente del elemento (escalado por λ_final) y las fuerzas internas de extremo son `F_int − f_eq` en ejes locales. Antes se devolvía `K·u`, que con carga distribuida son fuerzas nodales, no fuerzas internas: un voladizo bajo peso propio daba V_i = qL/2 y M_i = −5qL²/12 en vez de qL y −qL²/2, contradiciendo a las reacciones.
- **2026-09-29** · Sección por fibras, etapa 1 ([ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md), [FiberSection](FiberSection.md)). El modelo no lineal anterior —la $E_t$ de la deformación axial escalando toda la matriz local, con $G = E_t/[2(1+\nu)]$— desaparece: con fluencia axial anulaba flectores y torsor, con ablandamiento les cambiaba el signo. Dos caminos excluyentes: `material` + `A`, `Iy`, `Iz`, `J` (+ `nu`), la forma cerrada, sin cambio de resultados, sólo con materiales `IS_LINEAR_ELASTIC = True`; y `section` + `n_ip`, con $\mathbf e = [\varepsilon_0, \kappa_y, \kappa_z]$ integrado en Gauss-Legendre y torsión elástica desacoplada con el `GJ` de la sección (obligatorio). Contrato: parámetros `section` y `n_ip`; `A`, `Iy`, `Iz`, `J` obligatorios sólo con material; criterios 6 (guarda) y 7 (camino con sección). Corrección de notación en §3, sin efecto en código ni resultados: $\kappa_y = -d^2w/ds^2$ (antes $+d^2w/ds^2$), la definición coherente con $M_y = EI_y\,\kappa_y$, con $M_y > 0$ ⇒ tracción en $z > 0$ (Reglas.md §5) y con la matriz de §8, que ya usaba $w' = -\theta_y$.
- **2026-09-29** · Masa de flexión sin inercia rotatoria (hallazgo clas-4 de la auditoría; **cambio físico aprobado por el usuario**). La consistente sumaba $\tfrac{\rho I L}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$ en $\theta_y$ (con $I_y$) y en $\theta_z$ (con $I_z$): cinemática de Timoshenko, mientras que Euler-Bernoulli desprecia la inercia rotatoria por hipótesis. Se retira, por los dos caminos: flexión con la hermítica clásica en los dos planos (en $xz$, con $w' = -\theta_y$, los acoplamientos traslación-giro de signo opuesto) y concentrada HRZ $\rho A L^3/78$ en $\theta_y$ y $\theta_z$ en lugar de $\rho I_y L/2$ y $\rho I_z L/2$. La torsión no cambia ($\rho J_p L/6$ y $\rho J_p L/2$, $J_p = I_y + I_z$). Consecuencia: $I_y$ e $I_z$ ya sólo entran en la masa a través de $J_p$. Contrato: §12 y criterios 8 y 9. Los bloques de flexión se comparten ahora con los marcos planos (`_shared.py`).
- **2026-10-27** · Cargas variables y verificación por soluciones manufacturadas (ADR 0022, fase 4). `compute_body_load` acepta un campo `b(x)`, proyectado a ejes locales y repartido con las funciones de la rigidez (axial lineal, Hermite en los dos planos con `w' = −θy`); centrada en el eje, la carga no produce torsión. Campo interior (`line_field`: `u`, `v`, `w`, `θx`, `θy`, `θz`, `ε`, `θx'`, `κy`, `κz`) con `section_rigidities` (`EA`, `GJ`, `EIy`, `EIz`). Nueva aceptación «soluciones manufacturadas»: exactitud nodal, órdenes 2, 4, 3, 1 y 2 y torsión exacta (sin torsor distribuido, decisión del usuario, `θx` lineal). Se corrige además una nota desactualizada de `out_of_scope`: la fuerza de cuerpo ya la repartía `compute_body_load`.
