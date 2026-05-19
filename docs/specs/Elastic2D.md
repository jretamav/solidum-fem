# Elastic2D — Elasticidad lineal isótropa 2D (plane stress / plane strain)

> Spec **retroactiva**: el material existe desde el origen del subsistema 2D y está validado por tests preexistentes (todos los pipelines sólido 2D elástico). Esta spec lo documenta sin cambiar comportamiento — H-5.3 de la auditoría 2026-05-18.

---

## Especificación física

### 0. Descripción general

Modelo elástico lineal isótropo bidimensional, parametrizado por $(E, \nu)$. Soporta dos hipótesis cinemáticas 2D mutuamente excluyentes seleccionadas en construcción:

- **`plane_stress`** (default) — $\sigma_{zz} = \sigma_{xz} = \sigma_{yz} = 0$. Láminas o membranas delgadas en su plano (chapa traccionada, placa fina).
- **`plane_strain`** — $\varepsilon_{zz} = \varepsilon_{xz} = \varepsilon_{yz} = 0$. Cuerpos prismáticos largos cargados transversalmente (presa, túnel, cilindro largo).

Independiente de la velocidad, sin variables históricas. Es el material elástico canónico para todos los sólidos 2D del catálogo.

### 1. Descomposición de la deformación

No aplica — toda la deformación es elástica:

$$\boldsymbol\varepsilon = \boldsymbol\varepsilon^e$$

### 2. Ley elástica

$$\boldsymbol\sigma = \mathbf C_e\,\boldsymbol\varepsilon$$

con $\boldsymbol\varepsilon = [\varepsilon_{xx},\,\varepsilon_{yy},\,\gamma_{xy}]^\top$ en notación Voigt del proyecto ($\gamma_{xy} = 2\varepsilon_{xy}$) y $\boldsymbol\sigma = [\sigma_{xx},\,\sigma_{yy},\,\sigma_{xy}]^\top$.

**Matriz constitutiva plane stress** (3×3):

$$\mathbf C^{ps}_e = \frac{E}{1-\nu^2}\begin{bmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & (1-\nu)/2 \end{bmatrix}$$

**Matriz constitutiva plane strain** (3×3):

$$\mathbf C^{pe}_e = \frac{E}{(1+\nu)(1-2\nu)}\begin{bmatrix} 1-\nu & \nu & 0 \\ \nu & 1-\nu & 0 \\ 0 & 0 & (1-2\nu)/2 \end{bmatrix}$$

Tangente constante $\mathbf C_e$ en todo régimen. Simétrica y positiva definida en el rango admisible de $\nu$.

### 3. Superficie de fluencia / criterio de daño

No aplica — el material es indefinidamente elástico.

### 4-6. Flujo, endurecimiento, Kuhn-Tucker

No aplican.

### 7. Variables internas

Ninguna — `state_vars = None` se acepta y se devuelve intacto. Sin `PRIMARY_STATE_VAR`.

### 8. Notación Voigt

Convención del proyecto: $\boldsymbol\varepsilon = [\varepsilon_{xx},\,\varepsilon_{yy},\,\gamma_{xy}]$ con $\gamma_{xy} = 2\varepsilon_{xy}$ (deformación angular *engineering*). Los esfuerzos son los componentes tensoriales reales $[\sigma_{xx},\,\sigma_{yy},\,\sigma_{xy}]$.

---

## Formulación numérica

### 9. Esquema temporal

No aplica — sin historia.

### 10. Return mapping / actualización

No aplica. Evaluación reducida a producto matriz-vector:

$$\boldsymbol\sigma_{n+1} = \mathbf C_e\,\boldsymbol\varepsilon_{n+1}, \qquad \mathbf C^{\text{alg}} = \mathbf C_e$$

La matriz $\mathbf C_e$ se precalcula y guarda en el constructor — coste de `compute_state` es un producto 3×3 · 3 sin alocación adicional.

### 11. Tangente algorítmica consistente

$$\mathbf C^{\text{alg}} = \mathbf C_e$$

Simétrica (`IS_SYMMETRIC = True` por default). En plane stress la matriz es PD para $\nu \in (-1, 0.5)$; en plane strain también, con singularidad en $\nu \to 0.5$ (incompresible — requiere formulación mixta, no soportada).

### 12. Caveats numéricos

- **Plane strain $\nu \to 0.5$**: la matriz constitutiva diverge ($1-2\nu \to 0$). Recomendado $\nu \le 0.49$ para metales (típico $\nu=0.3$). El constructor rechaza $\nu \ge 0.5$ con `ValueError`.
- **Locking volumétrico** en plane strain con elementos de bajo orden (Quad4, Tri3) cuando $\nu \to 0.5$: declarado limitación arquitectural, sin mitigación implementada (B-bar, F-bar diferidas).
- Para $\nu < 0$ el material es admisible termodinámicamente (auxético) pero infrecuente en la práctica; el constructor lo permite hasta $\nu > -1$.

---

## Contrato de implementación

```yaml
name: Elastic2D
kind: material
status: validated

interface:
  strain_dim: 3
  primary_state_var: null
  is_symmetric: true

parameters:
  - { name: E,          type: float, required: true,  desc: "Módulo de Young (>0)" }
  - { name: nu,         type: float, required: true,  desc: "Coeficiente de Poisson ∈ (-1, 0.5)" }
  - { name: hypothesis, type: str,   required: false, default: "plane_stress",
      desc: "plane_stress | plane_strain" }
  - { name: density,    type: float, required: false, default: null,
      desc: "Densidad (kg/m³). Opcional al construir; obligatoria solo si se ensambla peso propio o masa (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(3,), state_vars=None) -> (σ: ndarray(3,), C_e: ndarray(3,3), state_vars)"
  strain_kind: "plano Voigt — [ε_xx, ε_yy, γ_xy] con γ_xy = 2·ε_xy"

state_variables: []

conventions:
  sign: "σ > 0 ⇔ tracción (Reglas §5)"
  voigt: "γ_xy = 2·ε_xy en deformación; σ_xy tensorial sin factor"
  hypothesis_dispatch: "selección por kwarg al construir; matriz C precalculada"

validity:
  - "E > 0"
  - "ν ∈ (-1, 0.5) estricto; ν = 0.5 rechazado (incompresible)"
  - "hypothesis ∈ {'plane_stress', 'plane_strain'}"
  - "pequeñas deformaciones (|ε| ≲ 1e-2)"

out_of_scope:
  - "anisotropía (ortótropo, monoclínico)"
  - "plasticidad ⇒ VonMises2D, DruckerPrager2D"
  - "daño ⇒ IsotropicDamage2D"
  - "incompresibilidad estricta (ν=0.5) ⇒ requiere formulación mixta u-p"
  - "grandes deformaciones"

acceptance:
  # Cubierto por tests existentes de elementos sólidos 2D que usan
  # Elastic2D como material por default.
  verification:
    - name: ley_de_hooke_3x3_exacta
      setup: "evaluar compute_state con ε arbitrario; comparar σ = C_e·ε en ambas hipótesis"
      expect: "σ = C_e·ε exacto componente a componente; tangente = C_e"
      tol_rel: 1.0e-14
    - name: simetria_C_e
      setup: "construir C_e en ambas hipótesis"
      expect: "C_e == C_e.T exacto"
      tol_abs: 1.0e-14
    - name: positividad_definida
      setup: "autovalores de C_e en ambas hipótesis con ν ∈ {0, 0.3, 0.49}"
      expect: "todos los autovalores estrictamente positivos"
      tol_abs: 0.0
    - name: plane_stress_traccion_uniaxial
      setup: "ε = (ε_0, -ν·ε_0, 0); plane_stress"
      expect: "σ_xx = E·ε_0, σ_yy = 0, σ_xy = 0 (uniaxial puro emergente)"
      tol_rel: 1.0e-12
    - name: plane_strain_confinamiento
      setup: "ε = (ε_0, 0, 0); plane_strain"
      expect: "σ_yy = ν/(1-ν) · σ_xx (reacción de confinamiento clásica)"
      tol_rel: 1.0e-12
    - name: rechazo_inputs_invalidos
      setup: "construir Elastic2D con E≤0, ν fuera de rango, hypothesis desconocida, density<0"
      expect: "ValueError con mensaje claro en cada caso"

references:
  - "Timoshenko S., Goodier J.N. (1970). Theory of Elasticity. McGraw-Hill. §3-4 (plane stress / plane strain)."
  - "Cook R.D., Malkus D.S., Plesha M.E., Witt R.J. (2002). Concepts and Applications of Finite Element Analysis. Wiley. §3.6."
  - "Bathe K.J. (2014). Finite Element Procedures. Prentice Hall. §4.2."
```

---

## Implementación

- **Archivo**: [fenix/materials/elastic_2d.py](../../fenix/materials/elastic_2d.py).
- **Clase**: `Elastic2D`, registrada vía `@MaterialRegistry.register`.
- **Despacho por hipótesis**: la matriz `C` se precalcula en el constructor según `hypothesis`. `compute_state` solo hace `self.C @ strain` sin ramificaciones — coste mínimo.
- **Validación temprana** en constructor: `E > 0`, `ν ∈ (-1, 0.5)`, `hypothesis ∈ {plane_stress, plane_strain}`, `density ≥ 0` si se pasa. Mensajes en español.
- **Tests**: cobertura implícita amplísima vía todos los pipelines de sólidos 2D (Quad4, Quad8, Quad9, Tri3, Tri6) que usan Elastic2D como material por default (`tests/test_solid_2d_*.py`, decenas de tests).

---

## Diálogo

- **2026-05-19** · Spec creada retroactivamente para cerrar el hueco H-5.3 de la auditoría 2026-05-18. Material anterior a la convención de specs validadas; no se modifica código. Aceptación cumplida por cobertura indirecta en pipelines existentes.
