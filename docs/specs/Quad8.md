# Quad8 — cuadrilátero serendípito 2D de orden 2

> Elemento sólido 2D de orden cuadrático con interpolación serendípita.
> Reproduce exactamente campos hasta $Q_2^{seren}$ (todos los polinomios de grado total ≤ 2 más algunos de orden 3 sin el término $\xi^2\eta^2$). Mucho más preciso que `Quad4` en flexión y geometrías curvas con el mismo grado de malla.

---

## Especificación física

### 0. Descripción general
Cuadrilátero plano de 8 nodos: 4 vértices + 4 nodos en los puntos medios de los bordes. Pequeñas deformaciones, isoparamétrico.

### 1. Cinemática de desplazamientos
$$u_x(\xi, \eta) = \sum_{i=1}^{8} N_i(\xi, \eta)\, u_x^{(i)}, \qquad u_y(\xi, \eta) = \sum_{i=1}^{8} N_i(\xi, \eta)\, u_y^{(i)}.$$

### 2. Cinemática de deformaciones
Voigt 2D $[\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]$, idéntica al `Quad4`.

### 3. Ecuación constitutiva
Material 2D (`STRAIN_DIM = 3`).

### 4-5. Equilibrio
Idénticas a `Quad4`.

---

## Formulación numérica

### 6. Discretización
Nodos en orden estándar: 0-3 vértices antihorarios desde $(-1, -1)$; 4 medio del borde 0-1, 5 medio del 1-2, 6 medio del 2-3, 7 medio del 3-0. Mapeo isoparamétrico.

### 7. Funciones de forma serendípitas
Vértices ($i = 0, 1, 2, 3$ con $(\xi_i, \eta_i) \in \{\pm 1\}^2$):
$$N_i = \tfrac{1}{4}(1 + \xi_i\xi)(1 + \eta_i\eta)(\xi_i\xi + \eta_i\eta - 1).$$

Medios de borde:
- Bordes con $\xi = 0$ (nodos 4, 6): $N_i = \tfrac{1}{2}(1 - \xi^2)(1 + \eta_i\eta)$ con $\eta_i = \mp 1$.
- Bordes con $\eta = 0$ (nodos 5, 7): $N_i = \tfrac{1}{2}(1 + \xi_i\xi)(1 - \eta^2)$ con $\xi_i = \pm 1$.

### 8. Matriz $\mathbf B$
Ensamblaje estándar a partir de $\partial N_i / \partial x = \mathbf J^{-1} \partial N_i / \partial \xi$. Tamaño $3 \times 16$.

### 9-10. Rigidez y fuerzas internas
$$\mathbf K_e = \int \mathbf B^\top \mathbf D\, \mathbf B\, t\, dA, \qquad \mathbf F_{\text{int}}^e = \int \mathbf B^\top \boldsymbol\sigma\, t\, dA.$$
Cuadratura Gauss $3 \times 3$ por defecto.

### 11. Cuadratura
Gauss–Legendre $3 \times 3$ (9 puntos) — exacta para polinomios hasta grado 5, suficiente para integrar $\mathbf B^\top \mathbf D \mathbf B$ del Quad8 sobre Jacobiano constante. La opción `2x2` queda disponible (reducida) pero introduce modos espurios; emite advertencia.

### 12. Cargas distribuidas consistentes
- **Body load**: $\int N^\top \mathbf b\, t\, dA$ con la cuadratura del elemento (3×3).
- **Edge traction** uniforme: en un borde recto con N cuadráticas la integral analítica reparte $L/6$ a cada nodo extremo y $4L/6$ al nodo medio (regla 1-4-1, Simpson). Bordes:

| Edge | Nodos extremos | Nodo medio |
|------|----------------|------------|
| 0    | (0, 1) | 4 |
| 1    | (1, 2) | 5 |
| 2    | (2, 3) | 6 |
| 3    | (3, 0) | 7 |

### 13. Salida por punto de Gauss
`compute_gauss_state(U)` con la misma estructura que `Quad4`: 9 puntos por defecto.

---

## Contrato de implementación

```yaml
name: Quad8
kind: element
status: validated

interface:
  dof_names: [ux, uy]
  n_nodes: 8
  strain_dim: 3
  n_integration_points: 9       # default Gauss 3×3

parameters:
  - { name: thickness, type: float, required: false, default: 1.0 }
  - { name: quadrature, type: tuple, required: false, default: "3x3" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, C_tangent, state')"
  strain_kind: "Voigt 2D [ε_xx, ε_yy, γ_xy]"

conventions:
  sign: "tensión positiva = tracción"
  voigt: "[xx, yy, xy] con γ_xy = 2·ε_xy"
  node_orientation: "vértices antihorarios; medios en el orden de los bordes"

validity:
  - "pequeños desplazamientos y deformaciones"

out_of_scope:
  - "grandes deformaciones"
  - "estabilización contra modos espurios bajo integración reducida"

acceptance:
  verification:
    - name: patch_lineal
      setup: "campo u = a₀ + a₁·x + a₂·y impuesto en el contorno"
      expect: "ε constante en todos los Gauss e igual al gradiente analítico"
      tol_rel: 1.0e-10
    - name: patch_cuadratico
      setup: "campo u_x = x²·1e-3 impuesto en los 8 nodos"
      expect: "ε_xx(x, y) = 2e-3·x exacto en todos los puntos de Gauss"
      tol_rel: 1.0e-10
  specific:
    - name: cargas consistentes — body load uniforme
      setup: "Σf = b·A·t (regular y distorsionado)"
      tol_rel: 1.0e-10
    - name: cargas consistentes — tracción uniforme en un borde
      setup: "Σf = t̄·L·t y reparto 1/6, 4/6, 1/6 en (vértice, medio, vértice)"
      tol_rel: 1.0e-12

references:
  - "Bathe K.-J., Finite Element Procedures, §5.3"
  - "Cook, Malkus, Plesha, Witt, Concepts and Applications of FEA, §6"
  - "Zienkiewicz O.C., The Finite Element Method, vol. 1, §8"
```

---

## Implementación

- **Archivo**: [fenix/elements/solid_2d.py](../../fenix/elements/solid_2d.py) · clase `Quad8`.
- **Tests**: [tests/test_higher_order_solid_2d.py](../../tests/test_higher_order_solid_2d.py).
