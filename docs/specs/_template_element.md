# <Nombre> — <descripción corta>

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.

---

## Especificación física

### 0. Descripción general
*Tipo de elemento, dimensión, grado, hipótesis cinemáticas globales (pequeñas deformaciones, corotacional, gran deformación…).*

### 1. Cinemática de desplazamientos
*Campo $u(\xi)$ aproximado: orden polinómico, continuidad ($C^0$, $C^1$…).*

### 2. Cinemática de deformaciones
*Componentes activas y su definición a partir de $u$. Notación Voigt si aplica.*

### 3. Ecuación constitutiva
*Relación $\sigma$–$\varepsilon$ usada (o referencia al material). Variables internas si aplica.*

### 4. Equilibrio — forma fuerte
*PDE de gobierno con condiciones de contorno.*

### 5. Equilibrio — forma débil
*Forma variacional integrada por partes; espacios de prueba y test.*

---

## Formulación numérica (FEM)

### 6. Discretización
*Geometría del elemento, sistema local vs global, transformaciones (cosenos directores, jacobiano).*

### 7. Funciones de forma
*$N_i(\xi)$ explícitas o referencia al esquema (Lagrange, serendípito, jerárquico…).*

### 8. Matriz deformación-desplazamiento
*$\mathbf B(\xi)$ local y, si aplica, transformación a globales.*

### 9. Rigidez elemental
*$\mathbf K_e = \int \mathbf B^\top \mathbf D\, \mathbf B\, dV$ — analítica o cuadratura.*

### 10. Fuerzas internas
*$\mathbf F_{\text{int}} = \int \mathbf B^\top \boldsymbol\sigma\, dV$.*

### 11. Cuadratura
*Esquema (Gauss–Legendre 2×2, etc.), nº de puntos, justificación (orden de exactitud, locking, hourglassing).*

---

## Contrato de implementación

```yaml
name: <Nombre>
kind: element
status: draft            # draft → implemented → validated

interface:
  dof_names: [ux, uy]    # DOFs por nodo
  n_nodes: 0
  strain_dim: 0          # 1 axial · 3 plano (ε_xx, ε_yy, γ_xy) · 6 3D
  n_integration_points: 0

parameters:
  - { name: , type: , required: true, desc:  }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: <axial escalar | plano Voigt | 3D Voigt>

conventions:
  sign: ""
  voigt: ""              # orden de componentes, factor en cortante
  node_orientation: ""   # antihorario, libre, etc.

validity:
  - ""

out_of_scope:
  - ""

acceptance:
  # Bloque de verificación obligatorio para todo elemento nuevo.
  # Sigue la disciplina V&V de NAFEMS (verification = el código resuelve
  # bien las ecuaciones; validation = las ecuaciones modelan bien la física).
  verification:
    - name: patch_test_macneal_harder
      setup: "malla con nodos interiores y elementos distorsionados; campo
             lineal u = a₀ + a₁·x + a₂·y, v = b₀ + b₁·x + b₂·y impuesto en
             los nodos del contorno"
      expect: "nodos interiores adoptan el campo lineal y ε es constante e
              igual a (a₁, b₂, a₂+b₁) en todos los puntos de Gauss"
      tol_rel: 1.0e-10
    - name: nafems_LE_<n>
      setup: "<benchmark NAFEMS aplicable: LE1 elasticidad plana, LE10 placa
             gruesa, T1 térmico…>"
      expect: "<magnitud y punto de referencia tabulados en el benchmark>"
      tol_rel: <según tabla del benchmark>
  # Tests adicionales específicos del elemento (degeneración, modos de
  # cuerpo rígido, particularidades de la formulación).
  specific:
    - name: ""
      setup: ""
      expect: ""
      tol_rel: 1.0e-10

references:
  - ""
```

---

## Implementación

*Rellena la IA tras programar.*

- Archivo: —
- Clase: —
- Tests:
  - —
- Notas de traducción: —

---

## Diálogo

*Preguntas, aclaraciones y hallazgos durante la implementación. Entradas fechadas.*

- *(vacío)*
