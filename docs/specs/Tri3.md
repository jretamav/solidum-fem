# Tri3 — triángulo lineal de deformación constante (CST) 2D

> Documentación retroactiva del elemento ya implementado.

---

## Especificación física

### 0. Descripción general
Elemento sólido bidimensional plano de primer orden, triangular de tres nodos. Aproximación lineal $C^0$ del campo de desplazamientos en coordenadas naturales $(\xi, \eta)$. Conocido en la literatura como elemento *constant strain triangle* (CST): la matriz $\mathbf{B}$ es constante sobre el elemento, por lo que la deformación $\boldsymbol{\varepsilon}$ y la tensión $\boldsymbol{\sigma}$ son uniformes dentro de cada Tri3. Régimen de pequeños desplazamientos y deformaciones.

### 1. Cinemática de desplazamientos
$$u_x(\xi, \eta) = \sum_{i=1}^{3} N_i(\xi, \eta)\, u_x^{(i)}, \qquad u_y(\xi, \eta) = \sum_{i=1}^{3} N_i(\xi, \eta)\, u_y^{(i)}.$$

### 2. Cinemática de deformaciones
Notación Voigt 2D $[\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]$, idéntica a la de `Quad4`. Por la linealidad de $N_i$, las derivadas $\partial N_i / \partial x$ son constantes y $\boldsymbol{\varepsilon}$ es uniforme en el elemento.

### 3. Ecuación constitutiva
Delegada al material 2D (`STRAIN_DIM = 3`). Soporta `Elastic2D`, `IsotropicDamage2D`, `VonMises2D` (este último solo `plane_strain`).

### 4. Equilibrio — forma fuerte
$$-\nabla \cdot \boldsymbol{\sigma} = \mathbf{b} \quad \text{en } \Omega, \quad \boldsymbol{\sigma} \cdot \mathbf{n} = \bar{\mathbf{t}} \text{ en } \partial\Omega_N, \quad \mathbf{u} = \bar{\mathbf{u}} \text{ en } \partial\Omega_D.$$

### 5. Equilibrio — forma débil
Idéntica a la del `Quad4`; el principio variacional no depende de la forma del elemento.

---

## Formulación numérica (FEM)

### 6. Discretización
Tres nodos en orden antihorario sobre el plano $(x, y)$. Mapeo isoparamétrico desde el triángulo de referencia con vértices en $(0, 0)$, $(1, 0)$, $(0, 1)$. Jacobiano $\mathbf{J} = \partial(x, y) / \partial(\xi, \eta)$ es constante sobre el elemento; el código aborta con error si $\det \mathbf{J} \le \text{tol}$ (nodos colineales o invertidos).

### 7. Funciones de forma
Coordenadas baricéntricas:
$$N_1 = 1 - \xi - \eta, \qquad N_2 = \xi, \qquad N_3 = \eta.$$
Sus derivadas en coordenadas naturales son constantes:
$$\partial N_1 / \partial \xi = -1, \quad \partial N_2 / \partial \xi = 1, \quad \partial N_3 / \partial \xi = 0;$$
$$\partial N_1 / \partial \eta = -1, \quad \partial N_2 / \partial \eta = 0, \quad \partial N_3 / \partial \eta = 1.$$

### 8. Matriz deformación-desplazamiento
$\mathbf{B}$ tiene tamaño $3 \times 6$ y, una vez transformadas las derivadas a globales vía $\mathbf{J}^{-1}$, sus entradas son **constantes** sobre el elemento:
$$\mathbf{B}_i = \begin{bmatrix} \partial N_i / \partial x & 0 \\ 0 & \partial N_i / \partial y \\ \partial N_i / \partial y & \partial N_i / \partial x \end{bmatrix}.$$

### 9. Rigidez elemental
$$\mathbf{K}_e = \mathbf{B}^\top\, \mathbf{D}\, \mathbf{B}\, A_e\, t,$$
con $A_e = \tfrac{1}{2}\, \det \mathbf{J}$ el área del triángulo y $t$ el espesor. La integración es exacta con un único punto central porque $\mathbf{B}$ es constante.

### 10. Fuerzas internas
$$\mathbf{F}_{\text{int}} = \mathbf{B}^\top\, \boldsymbol{\sigma}\, A_e\, t.$$

### 11. Cuadratura
Un punto central (en el baricentro) con peso $w = 1/2$ sobre el triángulo de referencia. Es exacto para formas constantes de $\mathbf{B}$ y $\boldsymbol{\sigma}$.

---

## Limitación crítica — *shear locking*

El Tri3 es notoriamente rígido en flexión: tres DOFs de desplazamiento no bastan para representar el modo de flexión puro de una viga discretizada en triángulos, por lo que el elemento responde con cortante espurio (*shear locking*) y subestima sistemáticamente la deflexión. La severidad crece con la relación $L/h$ del problema. Recomendación operativa: usar `Quad4` por defecto; reservar `Tri3` para zonas de transición de malla en regiones donde el campo de tensiones varía suavemente.

---

## Contrato de implementación

```yaml
name: Tri3
kind: element
status: validated

interface:
  dof_names: [ux, uy]
  n_nodes: 3
  strain_dim: 3                 # Voigt 2D [ε_xx, ε_yy, γ_xy]
  n_integration_points: 1

parameters:
  - { name: thickness, type: float, required: false, default: 1.0, desc: Espesor para estado plano }

material_contract:
  signature: "compute_state(ε, state) -> (σ, C_tangent, state')"
  strain_kind: "Voigt 2D [ε_xx, ε_yy, γ_xy]"

conventions:
  sign: "tensión positiva = tracción"
  voigt: "[xx, yy, xy] con γ_xy = 2·ε_xy"
  node_orientation: "antihorario (asegura det J > 0)"

validity:
  - "pequeños desplazamientos y deformaciones"
  - "campos de tensión que varían suavemente sobre el elemento"

out_of_scope:
  - "flexión dominante (shear locking severo)"
  - "elementos de alto orden (Tri6)"
  - "grandes deformaciones"

acceptance:
  - name: parche de tracción uniaxial
    setup: "malla triangular con carga axial uniforme y material Elastic2D"
    expect: "deformación constante en todo el dominio"
    tol_rel: 1.0e-12
  - name: jacobiano degenerado abortado
    setup: "tres nodos colineales"
    expect: "ValueError en _compute_kinematics_tri3"

references:
  - "Cook, Malkus, Plesha, Witt, Concepts and Applications of FEA, §3.6 (CST)"
  - "Zienkiewicz O.C., The Finite Element Method, vol. 1, §5 (limitations of CST)"
```

---

## Implementación

- **Archivo**: [fenix/elements/solid_2d.py](../../fenix/elements/solid_2d.py)
- **Clase**: `Tri3` (registrada vía `@ElementRegistry.register`)
- **Función núcleo**: `_compute_kinematics_tri3(coords)`, decorada con `@njit`. Reutiliza `_compute_integrands` con peso $w = 0{.}5$.
- **Tests**: [tests/test_solid_2d.py](../../tests/test_solid_2d.py) — incluye el patch test sobre malla triangular y verifica que el modo de cuerpo rígido produce fuerzas internas nulas.

---

## Diálogo

- **2026-04-30** · Spec retroactiva. `Tri3` precedía al protocolo spec-first; la spec se creó documentando la formulación ya implementada. Promovido `status: validated`.
