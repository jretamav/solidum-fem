# CABLE 2D COROTACIONAL — tensión unilateral, grandes rotaciones

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Elemento 1D inmerso en el plano 2D que modela un **cable**: dos nodos articulados que solo transmiten **esfuerzo axial de tensión**. En compresión o estado sin tensar, el elemento no transmite nada — se destensa. Cinemática **corotacional** (Updated Lagrangian): el eje rota con el movimiento, la longitud se mide sobre la configuración corriente.

La propiedad de unilateralidad vive en el **material**; el elemento es la maquinaria cinemática que aplica al material la deformación correcta. Para obtener el comportamiento propio de cable se empareja este elemento con un material unilateral como `CableMaterial1D`. Con un material elástico clásico, el elemento se comporta como una barra articulada corotacional (útil para pruebas, no para modelado realista de cables).

### 1. Cinemática
Configuraciones:
- Referencia: $\mathbf X_1, \mathbf X_2 \in \mathbb{R}^2$, longitud $L_0 = \|\mathbf X_2 - \mathbf X_1\|$.
- Corriente: $\mathbf x_i = \mathbf X_i + \mathbf u_i$, longitud $l = \|\mathbf x_2 - \mathbf x_1\|$.

Cosenos directores del eje **corriente**:
$$c_\theta = \frac{x_2 - x_1}{l}, \qquad s_\theta = \frac{y_2 - y_1}{l}.$$

### 2. Medida de deformación
Ingenieril corotacional:
$$\varepsilon = \frac{l - L_0}{L_0}.$$

### 3. Ecuación constitutiva
Delegada al material con `STRAIN_DIM = 1`. El elemento invoca `material.compute_state(ε, state)` y recibe $(\sigma, E_t, \text{nuevo estado})$ sin inspeccionar su naturaleza. La unilateralidad (si la hay) es responsabilidad del material — el elemento nunca filtra $\sigma$ ni $\varepsilon$ por signo.

### 4. Equilibrio — forma débil incremental
Principio de trabajos virtuales, lado interno:
$$\delta W_{\text{int}} = \int_V \delta\varepsilon\,\sigma\,dV = \delta\mathbf u_e^\top\,\mathbf F_{\text{int}}, \qquad \mathbf K_T = \frac{\partial\mathbf F_{\text{int}}}{\partial\mathbf u_e} = \mathbf K_M + \mathbf K_G.$$

Las cargas externas las ensambla el solver. El elemento no calcula vector nodal equivalente de fuerzas de cuerpo.

---

## Formulación numérica (FEM)

### 6. Discretización
Dos nodos, 2 DOFs por nodo ($u_x, u_y$). Vector elemental:
$$\mathbf u_e = [u_x^{(1)},\; u_y^{(1)},\; u_x^{(2)},\; u_y^{(2)}]^\top.$$

### 7. Vectores de dirección (configuración corriente)
$$\mathbf d = [-c_\theta,\; -s_\theta,\; c_\theta,\; s_\theta]^\top, \qquad \mathbf n = [-s_\theta,\; c_\theta,\; s_\theta,\; -c_\theta]^\top.$$

$\mathbf d$ apunta a lo largo del eje corriente (norma $\sqrt 2$); $\mathbf n$ al perpendicular (norma $\sqrt 2$). Son ortogonales: $\mathbf d^\top\mathbf n = 0$.

### 8. Matriz B corotacional
$$\mathbf B = \frac{1}{L_0}\,\mathbf d^\top, \qquad \varepsilon = \mathbf B\,\mathbf u_e.$$

### 9. Rigidez tangente
$$\mathbf K_T = \mathbf K_M + \mathbf K_G,$$
$$\mathbf K_M = \frac{E_t\,A}{L_0}\,\mathbf d\,\mathbf d^\top, \qquad \mathbf K_G = \frac{N}{l}\,\mathbf n\,\mathbf n^\top, \qquad N = \sigma A.$$

**Caso destensado** (material unilateral con $\varepsilon \leq 0$): el material devuelve $\sigma = 0$, $E_t = 0$, luego $N = 0$ y **ambas contribuciones se anulan** — $\mathbf K_T = \mathbf 0$. El elemento aporta cero al sistema global, como físicamente corresponde a un cable flojo.

### 10. Fuerzas internas
$$\mathbf F_{\text{int}} = N\,\mathbf d.$$

Cuando el cable está destensado ($N = 0$), $\mathbf F_{\text{int}} = \mathbf 0$.

### 11. Cuadratura
No aplica (forma cerrada).

---

## Contrato de implementación

```yaml
name: Cable2DCorot
kind: element
status: validated

interface:
  dof_names: [ux, uy]
  n_nodes: 2
  strain_dim: 1
  n_integration_points: 1

parameters:
  - { name: A, type: float, required: true, desc: "Área de la sección transversal del cable" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: axial escalar (deformación ingenieril corotacional)
  expected_behaviour: "material unilateral (p. ej. CableMaterial1D) para obtener respuesta física de cable; cualquier material con STRAIN_DIM=1 es técnicamente aceptado pero producirá comportamiento de barra"

conventions:
  sign: "σ > 0 ⇔ ε > 0 (elongación); la responsabilidad de anular σ en compresión es del material"
  voigt: "N/A (escalar)"
  node_orientation: "libre; K_T y F_int invariantes bajo permutación"
  configuration: "Updated Lagrangian — l, c_θ, s_θ se recalculan en cada evaluación"

validity:
  - "|ε| ≲ 1e-2 en régimen tensado"
  - rotaciones y desplazamientos de cualquier magnitud
  - cargas axiales (las transversales se transmiten por el destensado-retensado del cable)

out_of_scope:
  - flexión, cortante, momentos
  - dinámica con inercia distribuida (sin matriz de masa)
  - pandeo (irrelevante: el cable no transmite compresión)
  - "fuerzas de cuerpo distribuidas: el usuario reparte masa/peso a los nodos"
  - pretensión mediante parámetro del elemento (modelar con longitud inicial ajustada)

numerical_caveats:
  - "Un cable completamente destensado aporta K_T = 0 al sistema global. Si otros elementos no garantizan rigidez suficiente en los DOFs compartidos, la matriz tangente global puede ser singular o mal condicionada."
  - "En transiciones tensado ↔ destensado (ε cruza 0), Newton-Raphson puede oscilar por el salto discontinuo en E_t. Usar arc-length o pasos de carga finos si el problema atraviesa destensiones."

acceptance:
  - name: cable tensado coincide con barra corotacional
    setup: "cable horizontal con material lineal (Elastic1D, E=1000), u_x nodo 2 = 1e-3 ⇒ ε = 1e-3 > 0"
    expect: "σ = E·ε, E_t = E, K_T = K_M + K_G con valores de barra corotacional"
    tol_rel: 1.0e-10

  - name: cable destensado produce aportación nula
    setup: "cable horizontal con CableMaterial1D, u_x nodo 2 = -1e-3 ⇒ ε < 0"
    expect: "σ = 0, F_int = 0, K_T = matriz cero"
    tol_abs: 1.0e-12

  - name: invariancia bajo rotación rígida
    setup: "cable inicialmente horizontal, rotación rígida de 45° de ambos nodos"
    expect: "σ = 0 y F_int = 0 (independiente del material)"
    tol_abs: 1.0e-10

  - name: cruce por cero
    setup: "u_e tal que ε = 0 exactamente, material unilateral"
    expect: "σ = 0, E_t = 0, F_int = 0, K_T = 0"
    tol_abs: 1.0e-12

references:
  - "Crisfield M.A., Non-linear Finite Element Analysis of Solids and Structures, vol.1, §3.3"
  - "Irvine H.M., Cable Structures, MIT Press, §1.2"
```

---

## Implementación

- **Archivo**: [solidum/elements/cable.py](../../solidum/elements/cable.py) — archivo propio, no referencia a los módulos de armadura ni a ningún otro elemento.
- **Clase**: `Cable2DCorot` — hereda directamente de `Element` (no de `Truss2DCorot` ni de ningún elemento de armadura). La maquinaria cinemática corotacional se reimplementa íntegra dentro de la clase.
- **Tests**: [tests/test_cable_elements.py](../../tests/test_cable_elements.py) · `TestCable2DCorot` — los cuatro criterios de `acceptance` más un test de registro:
  - `test_acceptance_cable_tensado_como_barra_corot` (criterio 1).
  - `test_acceptance_cable_destensado_aportacion_nula` (criterio 2).
  - `test_acceptance_rotacion_rigida` (criterio 3).
  - `test_acceptance_cruce_por_cero` (criterio 4).
  - `test_registro_en_registry` — autodiscover.
- **Notas de traducción**:
  - La clase duplica deliberadamente la lógica de `Truss2DCorot`. No hay herencia. Es el precio de la independencia: si mañana se toca la armadura corotacional, este cable no se mueve.
  - `_current_geometry(u_e)` devuelve `(l, c_θ, s_θ)` — método privado para aislar la única operación que depende de configuración corriente.
  - En estado destensado, $\mathbf K_M = \mathbf 0$ (porque $E_t = 0$ viene del material) y $\mathbf K_G = \mathbf 0$ (porque $N = 0$); el elemento aporta literalmente ceros al ensamblaje, como debe ser físicamente un cable flojo. Los tests lo verifican con tolerancia absoluta.
  - El elemento no valida que el material sea unilateral: es responsabilidad del usuario emparejar con `CableMaterial1D`. Emparejar con `Elastic1D` produce respuesta de barra corotacional y se usa internamente en el test del criterio 1 para comparar contra valores analíticos.

---

## Diálogo

- **2026-04-21** · Elemento creado como componente totalmente independiente de `Truss2DCorot` por decisión explícita del usuario: la maquinaria cinemática se duplica dentro de la clase de cable para desacoplar su evolución futura del elemento de armadura. El coste es ~40 líneas duplicadas; la ganancia es que ninguna refactorización de las armaduras puede contaminar silenciosamente el comportamiento de los cables.
- **2026-04-21** · Validación de material unilateral: se optó por **no** imponerla en `__init__`. Motivo: el contrato elemento↔material del proyecto es genérico sobre `STRAIN_DIM=1`; introducir una verificación específica crearía un acoplamiento entre la clase del elemento y la identidad del material. Se documenta en `material_contract.expected_behaviour` del YAML para orientar al usuario.
