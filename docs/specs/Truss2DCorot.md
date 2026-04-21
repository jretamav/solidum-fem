# ARMADURA 2D COROTACIONAL — grandes desplazamientos, pequeñas deformaciones

> Orden de trabajo. Usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Barra articulada 2D equivalente a `Truss2D` pero en régimen de **grandes desplazamientos y rotaciones** con **pequeña deformación axial**. La cinemática se actualiza en configuración corriente (corotacional / Updated Lagrangian): el eje de la barra gira con el movimiento y la deformación se mide sobre la longitud actual. Rigidez tangente = rigidez material + rigidez geométrica debida al esfuerzo axial.

### 1. Cinemática
Configuraciones:
- Referencia: nodos $\mathbf X_1, \mathbf X_2$, longitud $L_0 = \|\mathbf X_2 - \mathbf X_1\|$.
- Corriente: $\mathbf x_i = \mathbf X_i + \mathbf u_i$, longitud $l = \|\mathbf x_2 - \mathbf x_1\|$.

Cosenos directores en configuración corriente:
$$c_\theta = \frac{x_2 - x_1}{l}, \qquad s_\theta = \frac{y_2 - y_1}{l}.$$

### 2. Medida de deformación
Ingenieril corotacional (válida para pequeña deformación aunque la rotación sea grande):
$$\varepsilon = \frac{l - L_0}{L_0}.$$

### 3. Ecuación constitutiva
Material 1D: $\sigma = E\,\varepsilon$ (elástico lineal); el contrato admite cualquier material con `STRAIN_DIM = 1`.

### 4. Equilibrio — forma débil incremental
Principio de trabajos virtuales completo:
$$\underbrace{\int_V \delta\varepsilon\, \sigma\, dV}_{\delta W_{\text{int}}} \;=\; \underbrace{\int_V \delta\mathbf u^\top \mathbf b\, dV + \int_{\partial V_t} \delta\mathbf u^\top \bar{\mathbf t}\, dA}_{\delta W_{\text{ext}}}.$$

El elemento calcula **solo el lado interno**:
$$\delta W_{\text{int}} = \delta\mathbf u_e^\top \mathbf F_{\text{int}},$$
con la rigidez tangente consistente separada en parte material y geométrica (ver §9). Las cargas externas (nodales y, si las hubiera, de cuerpo o tracción) las ensambla el solver a partir del input del usuario; este elemento **no** produce vector nodal equivalente de fuerzas de cuerpo distribuidas a lo largo del eje.

---

## Formulación numérica (FEM)

### 6. Discretización
Dos nodos, 2 DOFs por nodo ($u_x, u_y$). Vector elemental $\mathbf u_e = [u_x^{(1)}, u_y^{(1)}, u_x^{(2)}, u_y^{(2)}]^\top$.

### 7. Vector dirección actual
$$\mathbf d = [-c_\theta,\; -s_\theta,\; c_\theta,\; s_\theta]^\top, \qquad \mathbf n = [-s_\theta,\; c_\theta,\; s_\theta,\; -c_\theta]^\top,$$
con $c_\theta, s_\theta$ evaluados en la configuración corriente.

### 8. Matriz B corotacional
$$\mathbf B = \frac{1}{L_0}\,\mathbf d^\top, \qquad \varepsilon = \mathbf B\, \mathbf u_e.$$

### 9. Rigidez tangente
$$\mathbf K_T = \mathbf K_M + \mathbf K_G,$$
$$\mathbf K_M = \frac{E A}{L_0}\, \mathbf d\,\mathbf d^\top,$$
$$\mathbf K_G = \frac{N}{l}\, \mathbf n\,\mathbf n^\top, \qquad N = \sigma A.$$

### 10. Fuerzas internas
$$\mathbf F_{\text{int}} = N\, \mathbf d.$$

### 11. Cuadratura
No aplica (forma cerrada; un único punto conceptual en el eje).

---

## Contrato de implementación

```yaml
name: Truss2DCorot
kind: element
status: validated

interface:
  dof_names: [ux, uy]
  n_nodes: 2
  strain_dim: 1
  n_integration_points: 1

parameters:
  - { name: A, type: float, required: true, desc: Área de la sección transversal }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: axial escalar (deformación ingenieril corotacional)

conventions:
  sign: "tracción positiva (ε > 0 ⇔ σ > 0)"
  voigt: "N/A (escalar)"
  node_orientation: "libre; K_T y F_int invariantes bajo permutación"
  configuration: "Updated Lagrangian — c_θ, s_θ, l se recalculan en cada evaluación"

validity:
  - "|ε| ≲ 1e-2 (pequeña deformación axial)"
  - rotaciones y desplazamientos de cualquier magnitud
  - cargas exclusivamente axiales

out_of_scope:
  - flexión, cortante, momentos
  - pandeo por bifurcación (captura de snap-through requiere arc-length)
  - plasticidad con grandes deformaciones
  - "fuerzas de cuerpo distribuidas (peso propio sobre la barra): el elemento no calcula vector nodal equivalente; el usuario reparte masa a los nodos en el input"

acceptance:
  - name: equivalencia con Truss2D en régimen lineal
    setup: "barra axial empotrada, carga F pequeña (ε ~ 1e-6)"
    expect: "u(L) = F·L/(E·A) con tol_rel 1e-8"
    tol_rel: 1.0e-8

  - name: invariancia bajo rotación de cuerpo rígido
    setup: "aplicar rotación rígida de 45° a ambos nodos (sin deformar)"
    expect: "F_int = 0 y σ = 0"
    tol_abs: 1.0e-10

  - name: rigidez geométrica bajo tracción
    setup: "barra traccionada con N > 0; verificar autovalor transverso"
    expect: "autovalor de K_T en dirección n = N/l"
    tol_rel: 1.0e-10

references:
  - "Crisfield M.A., Non-linear Finite Element Analysis of Solids and Structures, vol.1, §3.3"
  - "Belytschko T., Nonlinear Finite Elements for Continua and Structures, §4.5"
```

---

## Implementación

- **Archivo**: [fenix/elements/structural.py](../../fenix/elements/structural.py)
- **Clase**: `Truss2DCorot` — hereda de `Truss2D` (reutiliza `__init__`, `ClassVar`s y metadatos de registro) y sobrescribe `compute_element_state` y `compute_internal_forces` para evaluar cosenos directores y longitud en configuración corriente.
- **Tests**: [tests/test_structural.py](../../tests/test_structural.py) · `TestTruss2DCorot` — tres tests, uno por cada `acceptance`:
  - `test_acceptance_linear_limit_matches_truss2d` (criterio 1).
  - `test_acceptance_rigid_body_rotation` (criterio 2).
  - `test_acceptance_geometric_stiffness_under_traction` (criterio 3).
- **Notas de traducción**:
  - La longitud corriente $l$ y los cosenos $c_\theta, s_\theta$ se recalculan en cada evaluación (no se cachean entre llamadas del solver, como exige Updated Lagrangian).
  - El vector $\mathbf d$ de §7 se construye con signos $[-c_\theta, -s_\theta, c_\theta, s_\theta]$; el transverso $\mathbf n$ como $[-s_\theta, c_\theta, s_\theta, -c_\theta]$. Ambos con norma $\sqrt 2$ — el factor se absorbe en los coeficientes $EA/L_0$ y $N/l$.
  - La rigidez tangente se devuelve ya sumada $\mathbf K_M + \mathbf K_G$; el solver no distingue las dos contribuciones.
  - El test 3 confronta $\mathbf K_T$ con la descomposición esperada y además verifica $\mathbf K_G \mathbf n = (2N/l)\mathbf n$ (autovalor sobre el vector sin normalizar, que es $N/l$ sobre el versor — así queda la ambigüedad resuelta explícitamente).

---

## Diálogo

- **2026-04-21** · Implementación inicial tras validar el diseño spec-first. Se optó por herencia `Truss2DCorot(Truss2D)` para reutilizar el constructor (L0, cosenos iniciales) y los `ClassVar`s del registro; los métodos que dependen de configuración corriente se sobrescriben. La bandera `large_strains` que tenía `Truss2D` previamente se elimina en el mismo commit — el régimen corotacional vive ahora como elemento independiente.
- **2026-04-21** · Normalización del autovalor transverso (criterio 3): la fórmula del libro da $N/l$ sobre el **versor** $\hat{\mathbf n} = \mathbf n/\sqrt 2$; sobre el vector $\mathbf n$ directo (norma $\sqrt 2$) el autovalor es $2N/l$. El test verifica la segunda forma para evitar raíces cuadradas y la nota de traducción documenta la equivalencia.
