# ELEMENTO FINITO SÓLIDO UNIDIMENSIONAL EN 3D

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.

---

## Especificación física

### 0. Descripción general
Elemento sólido 1D de primer orden, inmerso en el espacio tridimensional. Dos nodos articulados; transmite únicamente esfuerzo axial. Aproximación $C^0$ del desplazamiento. Régimen estrictamente de **linealidad geométrica**: pequeños desplazamientos, pequeñas rotaciones, pequeñas deformaciones.

### 1. Cinemática de desplazamientos
Interpolación lineal del desplazamiento axial a lo largo del eje $s \in [0, L]$:
$$u(s) = a_0 + a_1 s.$$

### 2. Cinemática de deformaciones
Única componente (axial, en el sistema local de la barra):
$$\varepsilon_{ss}(s) = \frac{du}{ds}.$$

### 3. Ecuación constitutiva
$$\sigma_{ss} = E\,\varepsilon_{ss}.$$
El elemento admite cualquier material 1D con `STRAIN_DIM = 1` (elástico lineal, plástico 1D, daño 1D).

### 4. Equilibrio — forma fuerte
$$-\frac{d}{ds}\!\left(EA\,\frac{du}{ds}\right) = b(s), \quad s \in (0, L),$$
con condiciones de contorno Dirichlet o Neumann ($\sigma A = \bar F$) en los extremos.

### 5. Equilibrio — forma débil
$$\int_0^L \frac{d\,\delta u}{ds}\,EA\,\frac{du}{ds}\,ds \;=\; \int_0^L \delta u\,b\,ds \;+\; \bigl[\delta u\,\bar F\bigr]_{\partial\Omega}, \quad \forall\,\delta u \in V_0.$$

El elemento calcula **solo el lado interno**; las cargas externas las ensambla el solver a partir del input del usuario.

---

## Formulación numérica (FEM)

### 6. Discretización
Dos nodos en $s = 0$ y $s = L$. Cosenos directores del eje en coordenadas globales:
$$c_x = \frac{x_2 - x_1}{L}, \quad c_y = \frac{y_2 - y_1}{L}, \quad c_z = \frac{z_2 - z_1}{L},$$
con $L = \|\mathbf X_2 - \mathbf X_1\|$ (configuración inicial, fija).

### 7. Funciones de forma
$$N_1(s) = 1 - \tfrac{s}{L}, \qquad N_2(s) = \tfrac{s}{L}.$$

### 8. Matriz deformación-desplazamiento
Local (1 DOF axial por nodo): $\mathbf B_{\text{loc}} = \tfrac{1}{L}[-1,\; 1]$.
Global con $\mathbf u_e = [u_x^{(1)}, u_y^{(1)}, u_z^{(1)}, u_x^{(2)}, u_y^{(2)}, u_z^{(2)}]^\top$:
$$\mathbf B = \tfrac{1}{L}\,[-c_x,\; -c_y,\; -c_z,\; c_x,\; c_y,\; c_z], \qquad \varepsilon = \mathbf B\,\mathbf u_e.$$

### 9. Rigidez elemental
Integración analítica exacta ($EA$ constante, $\mathbf B$ constante):
$$\mathbf K_e = \int_0^L \mathbf B^\top EA\,\mathbf B\,ds = \frac{EA}{L}\,\mathbf d\,\mathbf d^\top, \qquad \mathbf d = [-c_x,\; -c_y,\; -c_z,\; c_x,\; c_y,\; c_z]^\top.$$

$\mathbf K_e$ es una matriz $6\times 6$ de rango 1.

### 10. Fuerzas internas
$$\mathbf F_{\text{int}} = \sigma A\,\mathbf d = N\,\mathbf d.$$

### 11. Cuadratura
No aplica (forma cerrada, $\mathbf B$ y $\sigma$ uniformes por elemento). El campo `n_integration_points: 1` se refiere al número de estados materiales conceptuales (uno por elemento), no a puntos de Gauss.

---

## Contrato de implementación

```yaml
name: Truss3D
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
  strain_kind: axial escalar

conventions:
  sign: "tracción positiva (ε > 0 ⇔ σ > 0)"
  voigt: "N/A (escalar)"
  node_orientation: "libre; K_e y F_int invariantes bajo permutación"
  configuration: "inicial fija — L, c_x, c_y, c_z no se actualizan con los desplazamientos"

validity:
  - "|ε| ≲ 1e-2"
  - pequeños desplazamientos y rotaciones
  - cargas exclusivamente axiales

out_of_scope:
  - flexión, cortante, momentos
  - grandes desplazamientos o rotaciones (para ese régimen usar `Truss3DCorot`)
  - pandeo
  - "fuerzas de cuerpo distribuidas sobre el eje: el usuario reparte masa a los nodos"

acceptance:
  - name: barra axial bajo carga puntual (3D)
    setup: "empotrada en s=0, carga F axial en s=L, orientación arbitraria en el espacio"
    expect: "u(L) = F·L / (E·A) proyectado sobre el eje"
    tol_rel: 1.0e-10

references:
  - "Bathe K.-J., Finite Element Procedures, §4.2.1"
  - "Cook R.D., Concepts and Applications of FEA, §2.3"
```

---

## Implementación

- **Archivo**: [fenix/elements/truss.py](../../fenix/elements/truss.py)
- **Clase**: `Truss3D` — hereda directamente de `Element` (sin relación de herencia con `Truss2D` ni con `Truss2DCorot`).
- **Tests**: [tests/test_truss.py](../../tests/test_truss.py) · `TestTruss3D` — verifica `L0`, cosenos directores $(c_x, c_y, c_z)$, entradas de $\mathbf K_e$ (incluyendo simetría, dimensiones $6\times 6$) y evaluación de $\mathbf F_{\text{int}}$ sobre desplazamientos conocidos con geometría de prueba $L=13$, $(c_x,c_y,c_z)=(3/13, 4/13, 12/13)$.
- **Notas de traducción**:
  - `L0`, `cx`, `cy`, `cz` se calculan una vez en `__init__` sobre la configuración inicial y no se actualizan — coherente con el régimen geométricamente lineal declarado.
  - El constructor tolera nodos con 2 ó 3 coordenadas (completa con $z=0$ si falta), para facilitar transición de problemas planos embebidos.
  - El contrato con el material es el estándar del proyecto: `material.compute_state(ε, state) → (σ, Eₜ, state')`.
  - Para grandes desplazamientos/rotaciones en 3D usar `Truss3DCorot`, que hereda de esta clase.

---

## Diálogo

- **2026-04-21** · Apertura de spec retroactiva. La clase `Truss3D` preexistía al flujo spec-first; la spec se crea documentando la formulación ya implementada y estableciendo explícitamente su alcance: régimen de linealidad geométrica, sin bandera ni variante corotacional. Los tests unitarios existentes cubren el criterio de `acceptance` (vía los valores de $\mathbf K_e = (EA/L)\mathbf d\mathbf d^\top$ que implican algebraicamente $u(L)=FL/(EA)$). Promovido `status: validated`.
