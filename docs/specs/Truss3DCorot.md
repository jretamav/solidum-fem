# ARMADURA 3D COROTACIONAL — grandes desplazamientos, pequeñas deformaciones

> Orden de trabajo. Usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Barra articulada 3D equivalente a `Truss3D` pero en régimen de **grandes desplazamientos y rotaciones** con **pequeña deformación axial**. Updated Lagrangian / corotacional: el eje de la barra se redefine en cada iteración sobre la configuración corriente y la deformación se mide sobre la longitud actual. La rigidez tangente suma material + geométrica.

### 1. Cinemática
Configuraciones:
- Referencia: $\mathbf X_1, \mathbf X_2 \in \mathbb{R}^3$, longitud $L_0 = \|\mathbf X_2 - \mathbf X_1\|$.
- Corriente: $\mathbf x_i = \mathbf X_i + \mathbf u_i$, longitud $l = \|\mathbf x_2 - \mathbf x_1\|$.

Vector unitario del eje corriente:
$$\hat{\mathbf e} = \frac{\mathbf x_2 - \mathbf x_1}{l} = (c_x, c_y, c_z), \qquad c_x^2 + c_y^2 + c_z^2 = 1.$$

### 2. Medida de deformación
Ingenieril corotacional:
$$\varepsilon = \frac{l - L_0}{L_0}.$$

### 3. Ecuación constitutiva
Material 1D con `STRAIN_DIM = 1`: $\sigma = \sigma(\varepsilon)$ (el elemento no hardcodea la ley; la provee el material).

### 4. Equilibrio — forma débil incremental
PTV estándar. El elemento calcula solo el lado interno:
$$\delta W_{\text{int}} = \delta\mathbf u_e^\top\,\mathbf F_{\text{int}}, \qquad \mathbf K_T = \frac{\partial\mathbf F_{\text{int}}}{\partial\mathbf u_e} = \mathbf K_M + \mathbf K_G.$$

---

## Formulación numérica (FEM)

### 6. Discretización
Dos nodos, 3 DOFs por nodo ($u_x, u_y, u_z$). Vector elemental de 6 componentes:
$$\mathbf u_e = [u_x^{(1)}, u_y^{(1)}, u_z^{(1)}, u_x^{(2)}, u_y^{(2)}, u_z^{(2)}]^\top.$$

### 7. Vector dirección actual
$$\mathbf d = [-c_x,\; -c_y,\; -c_z,\; c_x,\; c_y,\; c_z]^\top,$$
con los cosenos evaluados en configuración **corriente**. $\|\mathbf d\|^2 = 2$.

### 8. Matriz B corotacional
$$\mathbf B = \frac{1}{L_0}\,\mathbf d^\top, \qquad \varepsilon = \mathbf B\,\mathbf u_e.$$

### 9. Rigidez tangente
**Material**:
$$\mathbf K_M = \frac{E_t\,A}{L_0}\,\mathbf d\,\mathbf d^\top \quad (6\times 6,\ \text{rango 1}).$$

**Geométrica**: en 3D la dirección perpendicular al eje es un **plano** (no un vector único). Se usa el proyector perpendicular $\mathbf P_3 = \mathbf I_3 - \hat{\mathbf e}\hat{\mathbf e}^\top$ (matriz $3\times 3$, rango 2). La rigidez geométrica en el elemento de 6 DOFs:
$$\mathbf K_G = \frac{N}{l}\,\begin{bmatrix}\mathbf P_3 & -\mathbf P_3 \\ -\mathbf P_3 & \mathbf P_3\end{bmatrix} \quad (6\times 6,\ \text{rango 2}), \qquad N = \sigma A.$$

### 10. Fuerzas internas
$$\mathbf F_{\text{int}} = N\,\mathbf d.$$

### 11. Cuadratura
No aplica (forma cerrada).

---

## Contrato de implementación

```yaml
name: Truss3DCorot
kind: element
status: validated

interface:
  dof_names: [ux, uy, uz]
  n_nodes: 2
  strain_dim: 1
  n_integration_points: 1

parameters:
  - { name: A, type: float, required: true, desc: Área de la sección transversal }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: axial escalar (deformación ingenieril corotacional)

conventions:
  sign: "σ > 0 ⇔ ε > 0 (elongación)"
  voigt: "N/A (escalar)"
  node_orientation: "libre"
  configuration: "Updated Lagrangian — l, c_x, c_y, c_z se recalculan en cada evaluación"

validity:
  - "|ε| ≲ 1e-2"
  - rotaciones y desplazamientos de cualquier magnitud en el espacio
  - cargas exclusivamente axiales

out_of_scope:
  - flexión, cortante, momentos
  - pandeo por bifurcación
  - plasticidad con grandes deformaciones
  - "fuerzas de cuerpo distribuidas: el usuario reparte masa a los nodos"

acceptance:
  - name: equivalencia con Truss3D en régimen lineal
    setup: "barra 3D con orientación arbitraria, carga axial pequeña (ε ~ 1e-6)"
    expect: "K_T y F_int coinciden con Truss3D lineal a tolerancia rtol=1e-4"
    tol_rel: 1.0e-4

  - name: invariancia bajo rotación rígida 3D
    setup: "aplicar una rotación rígida arbitraria (p. ej. 30° sobre eje z, luego 45° sobre eje x)"
    expect: "F_int = 0 y σ = 0"
    tol_abs: 1.0e-10

  - name: rigidez geométrica transversa (plano perpendicular)
    setup: "barra en tensión con N > 0; aplicar K_G sobre dos vectores perpendiculares al eje linealmente independientes"
    expect: "K_G·v_perp = (N/l)·v_perp_struct, donde v_perp_struct es la extensión del vector transverso al espacio de 6 DOFs con signos opuestos en los dos nodos"
    tol_rel: 1.0e-10

references:
  - "Crisfield M.A., Non-linear Finite Element Analysis of Solids and Structures, vol.1, §3.3"
  - "Belytschko T., Nonlinear Finite Elements for Continua and Structures, §4.5"
```

---

## Implementación

- **Archivo**: [fenix/elements/truss.py](../../fenix/elements/truss.py)
- **Clase**: `Truss3DCorot` — hereda de `Truss3D` (reutiliza `__init__`, tolerancia a nodos 2D/3D, metadatos de registro) y sobrescribe `compute_element_state` y `compute_internal_forces` para evaluar longitud y cosenos directores en configuración corriente.
- **Tests**: [tests/test_truss.py](../../tests/test_truss.py) · `TestTruss3DCorot` — tres tests, uno por cada `acceptance`:
  - `test_acceptance_linear_limit_matches_truss3d` (criterio 1).
  - `test_acceptance_rigid_body_rotation` (criterio 2 — dos rotaciones rígidas distintas, verifica σ=0 y F_int=0 en ambas).
  - `test_acceptance_geometric_stiffness_transverse_plane` (criterio 3 — aplica $\mathbf K_G$ a **dos** direcciones linealmente independientes del plano perpendicular al eje).
- **Notas de traducción**:
  - La clave respecto a 2D es el proyector $\mathbf P_3 = \mathbf I - \hat{\mathbf e}\hat{\mathbf e}^\top$ (3×3, rango 2): el plano perpendicular al eje es bidimensional. $\mathbf K_G$ resulta entonces de rango 2, no rango 1.
  - Método estático `_perpendicular_projector(cx, cy, cz)` extrae esta operación para legibilidad y posible reuso.
  - La longitud corriente $l$ y los cosenos $c_x, c_y, c_z$ se recalculan en cada evaluación (Updated Lagrangian); no se cachean entre llamadas del solver.
  - Autovalores de $\mathbf K_G$ sobre modos transversos (rotación de la barra alrededor de su centro): $2N/l$ sobre el vector sin normalizar — igual que en 2D, el factor 2 viene de que el vector de 6 componentes tiene norma $\sqrt 2$. Sobre el versor el autovalor es $N/l$.

---

## Diálogo

- **2026-04-21** · Implementación tras validar el diseño spec-first. Se optó por herencia `Truss3DCorot(Truss3D)` por paralelismo con `Truss2DCorot(Truss2D)` y reutilización del constructor. La diferencia estructural con el caso 2D es el proyector $\mathbf P$: en 2D tenía rango 1 y generaba un $\mathbf K_G$ de rango 1 (vector $\mathbf n$ único); en 3D tiene rango 2 y $\mathbf K_G$ rango 2 (plano perpendicular al eje). Los tests del criterio 3 verifican explícitamente **dos** direcciones transversas linealmente independientes, capturando este cambio dimensional.
- **2026-04-21** · Test de invariancia bajo rotación rígida 3D: se verifican dos rotaciones distintas (plano x-y y plano x-z) para descartar coincidencias accidentales en ejes canónicos.
