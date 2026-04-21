# ELEMENTO PÓRTICO/VIGA 3D EULER-BERNOULLI

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Elemento 1D inmerso en el espacio tridimensional que modela una **viga esbelta 3D**. Dos nodos rígidamente conectados, 6 DOFs por nodo ($u_x, u_y, u_z, r_x, r_y, r_z$). Captura seis acciones internas:

- **Axial** (esfuerzo normal $N$ a lo largo del eje local $x$).
- **Cortante en dos planos** ($V_y$ en el plano $xy$, $V_z$ en el plano $xz$).
- **Flexión en dos planos** ($M_z$ alrededor del eje local $z$, $M_y$ alrededor del eje local $y$).
- **Torsión** ($M_x$ alrededor del eje longitudinal).

Régimen de **linealidad geométrica** (pequeños desplazamientos y rotaciones). Hipótesis de Euler-Bernoulli: secciones planas perpendiculares al eje deformado (sin cizalladura transversal).

### 1. Hipótesis
- Secciones planas y perpendiculares al eje neutro deformado (igual que Frame2DEuler).
- Flexiones en los dos planos **desacopladas** entre sí (ejes principales de inercia).
- Torsión de Saint-Venant pura (sin alabeo). Válida para secciones macizas o cerradas; abiertas de pared delgada requerirían formulación con alabeo (fuera de alcance).
- Material isótropo 1D; $G = E/[2(1+\nu)]$.

### 2. Sistema de ejes locales
Tres ejes ortonormales en cada elemento:
- **$\hat{\mathbf x}$ local**: dirección del eje de la barra ($\hat{\mathbf x} = (\mathbf x_2 - \mathbf x_1)/L$).
- **$\hat{\mathbf y}$ local** y **$\hat{\mathbf z}$ local**: ejes principales de inercia de la sección, ortogonales al eje.

La orientación de $\hat{\mathbf y}, \hat{\mathbf z}$ no está determinada solo por la dirección del eje — hace falta un **vector de referencia** proporcionado por el usuario (o un default). Es la diferencia estructural respecto al caso 2D, donde el plano era único.

### 3. Medidas de deformación / acciones internas
- Axial: $\varepsilon_{axial} = du/ds$, tensión $N = EA\,\varepsilon_{axial}$ (o $N = A\,\sigma$ con $\sigma$ del material).
- Flexión en plano $xy$: curvatura $\kappa_z = d^2v/ds^2$, momento $M_z = EI_z\,\kappa_z$.
- Flexión en plano $xz$: curvatura $\kappa_y = d^2w/ds^2$, momento $M_y = EI_y\,\kappa_y$.
- Torsión: $d\theta_x/ds$, momento torsor $M_x = GJ\,d\theta_x/ds$.

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
con la componente axial corregida por el $\sigma$ que devuelve el material: $F_1 = -\sigma A$, $F_7 = +\sigma A$. El resto de componentes viene del producto lineal (permite materiales no-lineales en la rama axial).

### 10. Ensamblaje global
$$\mathbf K_{\text{global}} = \mathbf T^\top\,\mathbf K_{\text{local}}\,\mathbf T, \qquad \mathbf F_{\text{int}} = \mathbf T^\top\,\mathbf F_{\text{int,local}}.$$

### 11. Cuadratura
No aplica (forma cerrada; integración analítica de los polinomios de Hermite).

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
  n_integration_points: 1

parameters:
  - { name: A,  type: float, required: true, desc: "Área de la sección transversal" }
  - { name: Iy, type: float, required: true, desc: "Momento de inercia alrededor del eje local y" }
  - { name: Iz, type: float, required: true, desc: "Momento de inercia alrededor del eje local z" }
  - { name: J,  type: float, required: true, desc: "Constante torsional de Saint-Venant" }
  - { name: nu, type: float, required: false, default: 0.3, desc: "Coeficiente de Poisson (si el material no lo expone)" }
  - { name: ref_vector, type: list, required: false, default: "[0, 0, 1]",
      desc: "Vector 3D de referencia que fija la orientación de los ejes locales y, z de la sección. Default: eje z global." }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: "axial escalar (ε = (u_2 − u_1)/L en sistema local)"
  nonlinearity_model: "E_tangent escala toda la matriz local; G se recalcula como E_t/[2(1+ν)]"
  poisson_source: "material.nu si existe; en su defecto parámetro `nu` del elemento"

conventions:
  sign: "tracción positiva; momentos siguen regla de la mano derecha respecto a los ejes locales"
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
  - plasticidad distribuida en la sección (fibras)
  - pandeo por bifurcación
  - deformación por cortante transversal significativa (Timoshenko 3D, fuera de alcance)
  - "fuerzas de cuerpo distribuidas: el usuario reparte carga a los nodos"

numerical_caveats:
  - "Si ref_vector es paralelo al eje de la barra, la proyección al plano perpendicular se anula y el constructor aborta con mensaje claro — el usuario debe elegir otro vector."
  - "Para barras verticales (eje casi paralelo a z global) con ref_vector default, el constructor usa [1, 0, 0] como fallback y emite advertencia por consola."

acceptance:
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

references:
  - "Przemieniecki J.S., Theory of Matrix Structural Analysis, Tabla 11.3"
  - "Cook R.D., Concepts and Applications of FEA, §2.8"
  - "Bathe K.-J., Finite Element Procedures, cap. 5"
```

---

## Implementación

- **Archivo**: [fenix/elements/frame3d.py](../../fenix/elements/frame3d.py) — archivo propio, sin dependencia de ningún otro elemento.
- **Clase**: `Frame3D` — hereda directamente de `Element`. Métodos privados `_build_local_frame` (longitud, matriz $\boldsymbol{\lambda}$ 3×3 y $\mathbf T$ 12×12) y `_build_local_stiffness` (matriz $\mathbf K_{\text{local}}$ 12×12 con bloques axial, torsión y flexión en ambos planos).
- **Tests**: [tests/test_frame3d.py](../../tests/test_frame3d.py) · `TestFrame3DAcceptance` — cubre los 5 criterios + validación de `ref_vector` paralelo al eje + registro:
  - `test_acceptance_respuesta_axial_pura`
  - `test_acceptance_flexion_xy_coincide_con_frame2d`
  - `test_acceptance_flexion_xz`
  - `test_acceptance_torsion_pura`
  - `test_acceptance_simetria_K`
  - `test_rechazo_ref_vector_paralelo`
  - `test_registro_en_registry`
- **Notas de traducción**:
  - La matriz $\boldsymbol{\lambda}$ 3×3 se construye como $[\hat{\mathbf x}^\top; \hat{\mathbf y}^\top; \hat{\mathbf z}^\top]$ — las filas son los ejes locales expresados en coordenadas globales. $\mathbf T$ 12×12 es un `blkdiag` de cuatro copias de $\boldsymbol{\lambda}$ (traslaciones y rotaciones transforman como vectores en 3D infinitesimal).
  - El default de `ref_vector` es `[0, 0, 1]`; el fallback `[1, 0, 0]` se activa cuando $|\hat{\mathbf x} \cdot \hat{\mathbf z}_{\text{global}}| > 0.99$ (barra cercanamente vertical). Ambas rutas emiten advertencia por consola la primera vez que se usan con default.
  - La matriz $\mathbf K_{\text{local}}$ se construye rellenando entradas no nulas en sus 4 bloques desacoplados (axial, torsión, flexión $xy$, flexión $xz$). Los signos de la flexión $xz$ están invertidos respecto a la flexión $xy$ porque la pareja $(u_z, r_y)$ tiene orientación dextrógira opuesta a $(u_y, r_z)$: este punto se verifica explícitamente en `test_acceptance_flexion_xz`.
  - `compute_internal_forces` devuelve cortantes, torsión y momentos en ambos extremos separados — útil para post-proceso estructural (diagramas de esfuerzos).

---

## Diálogo

- **2026-04-21** · Elemento creado como componente totalmente autónomo: no hereda de `Frame2DEuler` ni de ninguna otra viga. La maquinaria cinemática 3D se implementa íntegra dentro de `frame3d.py`, coherente con la filosofía aplicada a cables.
- **2026-04-21** · Decisión sobre `ref_vector` vs "tercer nodo de orientación". Se optó por `ref_vector` (vector 3D proyectado al plano perpendicular al eje) por ser más compacto en YAML y no requerir nodos fantasma en la malla. El default `[0, 0, 1]` cubre la mayoría de casos de ingeniería estructural (estructuras terrestres con gravedad en $-z$); el fallback para barras verticales usa `[1, 0, 0]` con advertencia explícita.
- **2026-04-21** · Verificación cruzada con `Frame2DEuler`: el criterio 2 aplica una viga alineada con el eje $x$ global bajo carga en $y$; los desplazamientos resultantes deben coincidir con la solución 2D clásica $v = PL^3/(3EI_z)$, $\theta = PL^2/(2EI_z)$. El test lo verifica a 8 decimales. Esto blinda la coherencia dimensional: el 3D se reduce al 2D sin error acumulado cuando el problema es efectivamente plano.
