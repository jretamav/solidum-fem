# Quad4 — cuadrilátero bilineal isoparamétrico 2D

> Documentación retroactiva del elemento ya implementado. Cubre la formulación, el contrato de implementación y la trazabilidad con tests existentes.

---

## Especificación física

### 0. Descripción general
Elemento sólido bidimensional plano de primer orden, cuadrilátero de cuatro nodos. Aproximación bilineal $C^0$ del campo de desplazamientos en coordenadas naturales $(\xi, \eta) \in [-1, 1]^2$. Formulación isoparamétrica: la geometría y los desplazamientos comparten las mismas funciones de forma. Régimen de pequeños desplazamientos y deformaciones.

### 1. Cinemática de desplazamientos
$$u_x(\xi, \eta) = \sum_{i=1}^{4} N_i(\xi, \eta)\, u_x^{(i)}, \qquad u_y(\xi, \eta) = \sum_{i=1}^{4} N_i(\xi, \eta)\, u_y^{(i)}.$$

### 2. Cinemática de deformaciones
Tensor infinitesimal en notación Voigt 2D $[\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]$:
$$\varepsilon_{xx} = \frac{\partial u_x}{\partial x}, \quad \varepsilon_{yy} = \frac{\partial u_y}{\partial y}, \quad \gamma_{xy} = \frac{\partial u_x}{\partial y} + \frac{\partial u_y}{\partial x}.$$

### 3. Ecuación constitutiva
Delegada al material 2D (`STRAIN_DIM = 3`). El elemento llama `material.compute_state(ε, state) → (σ, C_tangent, state')` por punto de Gauss y soporta cualquier ley constitutiva 2D (`Elastic2D`, `IsotropicDamage2D`, `VonMises2D`).

### 4. Equilibrio — forma fuerte
$$-\nabla \cdot \boldsymbol{\sigma} = \mathbf{b} \quad \text{en } \Omega, \qquad \boldsymbol{\sigma} \cdot \mathbf{n} = \bar{\mathbf{t}} \quad \text{en } \partial\Omega_N, \qquad \mathbf{u} = \bar{\mathbf{u}} \quad \text{en } \partial\Omega_D.$$

### 5. Equilibrio — forma débil
$$\int_\Omega \delta \boldsymbol{\varepsilon}^\top \boldsymbol{\sigma}\, dV = \int_\Omega \delta \mathbf{u}^\top \mathbf{b}\, dV + \int_{\partial\Omega_N} \delta \mathbf{u}^\top \bar{\mathbf{t}}\, dS, \qquad \forall\, \delta \mathbf{u} \in V_0.$$

---

## Formulación numérica (FEM)

### 6. Discretización
Cuatro nodos en orden antihorario sobre el plano $(x, y)$. Mapeo isoparamétrico desde el cuadrado de referencia $[-1, 1]^2$:
$$x(\xi, \eta) = \sum_{i=1}^{4} N_i(\xi, \eta)\, x^{(i)}, \qquad y(\xi, \eta) = \sum_{i=1}^{4} N_i(\xi, \eta)\, y^{(i)}.$$
Jacobiano $\mathbf{J} = \partial(x, y) / \partial(\xi, \eta)$ y su determinante $\det \mathbf{J}$ entran en el cambio de variable de la integración. El código aborta con error si $\det \mathbf{J} \le \text{tol}$ (elemento degenerado o nodos en orden invertido).

### 7. Funciones de forma
Producto tensorial bilineal:
$$N_i(\xi, \eta) = \tfrac{1}{4}(1 + \xi_i\,\xi)(1 + \eta_i\,\eta), \qquad (\xi_i, \eta_i) = (\pm 1, \pm 1).$$
Convención de numeración (antihorario): $(\xi_1, \eta_1) = (-1, -1)$, $(\xi_2, \eta_2) = (+1, -1)$, $(\xi_3, \eta_3) = (+1, +1)$, $(\xi_4, \eta_4) = (-1, +1)$.

### 8. Matriz deformación-desplazamiento
Las derivadas en globales se obtienen vía $\partial N_i / \partial x = \mathbf{J}^{-1}\, \partial N_i / \partial \xi$. La matriz $\mathbf{B}$ de tamaño $3 \times 8$ se ensambla por nodo:
$$\mathbf{B}_i = \begin{bmatrix} \partial N_i / \partial x & 0 \\ 0 & \partial N_i / \partial y \\ \partial N_i / \partial y & \partial N_i / \partial x \end{bmatrix}, \qquad \boldsymbol{\varepsilon} = \mathbf{B}\, \mathbf{u}_e.$$

### 9. Rigidez elemental
$$\mathbf{K}_e = \int_\Omega \mathbf{B}^\top\, \mathbf{D}\, \mathbf{B}\, t\, dA = \sum_{p=1}^{n_g} \mathbf{B}(\xi_p, \eta_p)^\top\, \mathbf{D}_p\, \mathbf{B}(\xi_p, \eta_p)\, \det \mathbf{J}_p\, w_p\, t,$$
donde $t$ es el espesor (`thickness`) y $\mathbf{D}_p$ es el tensor constitutivo tangente devuelto por el material en cada punto.

### 10. Fuerzas internas
$$\mathbf{F}_{\text{int}} = \sum_{p=1}^{n_g} \mathbf{B}(\xi_p, \eta_p)^\top\, \boldsymbol{\sigma}_p\, \det \mathbf{J}_p\, w_p\, t.$$

### 11. Cuadratura
Por defecto Gauss–Legendre $2 \times 2$ (cuatro puntos): orden de exactitud 3 — integra exactamente $\mathbf{K}_e$ para Jacobiano constante (geometría de paralelogramo) y captura el comportamiento dominante en geometrías irregulares moderadas. Soporta esquemas alternativos vía el parámetro `quadrature` (lista nombrada en `QuadratureRegistry`).

**Aviso sobre integración reducida**: el esquema $1 \times 1$ está disponible pero el código emite advertencia explícita; un único punto central no detecta los modos de hourglass (energía nula con desplazamientos de signo alternado), que pueden contaminar la respuesta sin estabilización adicional. No se implementa estabilización tipo Flanagan-Belytschko.

---

## Contrato de implementación

```yaml
name: Quad4
kind: element
status: validated

interface:
  dof_names: [ux, uy]
  n_nodes: 4
  strain_dim: 3                 # Voigt 2D [ε_xx, ε_yy, γ_xy]
  n_integration_points: 4       # default Gauss 2×2

parameters:
  - { name: thickness, type: float, required: false, default: 1.0, desc: Espesor para estado plano }
  - { name: quadrature, type: tuple, required: false, default: "2x2", desc: Regla de cuadratura desde QuadratureRegistry }

material_contract:
  signature: "compute_state(ε, state) -> (σ, C_tangent, state')"
  strain_kind: "Voigt 2D [ε_xx, ε_yy, γ_xy]"

conventions:
  sign: "tensión positiva = tracción"
  voigt: "[xx, yy, xy] con γ_xy = 2·ε_xy (engineering shear strain)"
  node_orientation: "antihorario (asegura det J > 0)"

validity:
  - "pequeños desplazamientos y deformaciones"
  - "geometría sin distorsión severa (det J > tol en todos los puntos de Gauss)"
  - "compatible con materiales Elastic2D, IsotropicDamage2D, VonMises2D (este último solo plane_strain)"

out_of_scope:
  - "grandes deformaciones, formulaciones corotacionales o totales lagrangianas"
  - "estabilización contra hourglass para integración reducida"
  - "elementos de alto orden (Q8, Q9)"

acceptance:
  - name: parche de tracción uniaxial
    setup: "malla regular con carga axial uniforme y material Elastic2D"
    expect: "deformación constante ε_xx = σ_aplicada / E_efectivo en todo el dominio"
    tol_rel: 1.0e-12
  - name: jacobiano degenerado abortado
    setup: "elemento con nodos colineales o en orden inverso"
    expect: "ValueError en _compute_kinematics"

references:
  - "Bathe K.-J., Finite Element Procedures, §5.3 (isoparametric elements)"
  - "Cook, Malkus, Plesha, Witt, Concepts and Applications of FEA, §6 (plane stress/strain)"
```

---

## Implementación

- **Archivo**: [fenix/elements/solid_2d.py](../../fenix/elements/solid_2d.py)
- **Clase**: `Quad4` (registrada vía `@ElementRegistry.register`)
- **Funciones núcleo**: `_compute_kinematics(xi, eta, coords)` y `_compute_integrands(B, C, σ, detJ, w, t)`, ambas decoradas con `@njit` (Numba) para evaluar $\mathbf{B}$, $\det \mathbf{J}$ y los integrandos a velocidad cercana a Fortran.
- **Tests**: [tests/test_solid_2d.py](../../tests/test_solid_2d.py) — verifica el ensamblaje de $\mathbf{B}$, simetría de $\mathbf{K}_e$ en geometrías regulares, conservación de modos de cuerpo rígido (translación y rotación pura producen $\mathbf{F}_{\text{int}} = 0$), y el patch test de tracción uniaxial.

---

## Diálogo

- **2026-04-30** · Spec retroactiva. `Quad4` precedía al protocolo spec-first; la spec se creó documentando la formulación ya implementada y verificada. Promovido `status: validated`.
