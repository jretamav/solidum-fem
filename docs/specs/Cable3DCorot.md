# CABLE 3D COROTACIONAL — tensión unilateral, grandes rotaciones en el espacio

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Elemento 1D inmerso en el espacio 3D que modela un **cable**: dos nodos articulados que solo transmiten **esfuerzo axial de tensión**. En compresión o estado sin tensar, el elemento no transmite nada. Cinemática **corotacional** (Updated Lagrangian): el eje rota libremente en el espacio, la longitud se mide sobre la configuración corriente.

La unilateralidad vive en el **material**; el elemento es la maquinaria cinemática 3D. Para obtener respuesta física de cable, emparejar con `CableMaterial1D`. Con un material elástico clásico, el elemento se comporta como una barra corotacional 3D.

### 1. Cinemática
- Referencia: $\mathbf X_1, \mathbf X_2 \in \mathbb{R}^3$, longitud $L_0 = \|\mathbf X_2 - \mathbf X_1\|$.
- Corriente: $\mathbf x_i = \mathbf X_i + \mathbf u_i$, longitud $l = \|\mathbf x_2 - \mathbf x_1\|$.

Vector unitario del eje corriente:
$$\hat{\mathbf e} = \frac{\mathbf x_2 - \mathbf x_1}{l} = (c_x, c_y, c_z), \qquad c_x^2 + c_y^2 + c_z^2 = 1.$$

### 2. Medida de deformación
Ingenieril corotacional:
$$\varepsilon = \frac{l - L_0}{L_0}.$$

### 3. Ecuación constitutiva
Delegada al material con `STRAIN_DIM = 1`. El elemento llama `material.compute_state(ε, state)` y recibe $(\sigma, E_t, \text{nuevo estado})$ sin inspeccionar su naturaleza. La unilateralidad (si la hay) es responsabilidad del material.

### 4. Equilibrio — forma débil incremental
PTV estándar; el elemento calcula solo el lado interno:
$$\delta W_{\text{int}} = \delta\mathbf u_e^\top\,\mathbf F_{\text{int}}, \qquad \mathbf K_T = \frac{\partial\mathbf F_{\text{int}}}{\partial\mathbf u_e} = \mathbf K_M + \mathbf K_G.$$

---

## Formulación numérica (FEM)

### 6. Discretización
Dos nodos, 3 DOFs por nodo ($u_x, u_y, u_z$). Vector elemental de 6 componentes:
$$\mathbf u_e = [u_x^{(1)}, u_y^{(1)}, u_z^{(1)}, u_x^{(2)}, u_y^{(2)}, u_z^{(2)}]^\top.$$

### 7. Vector dirección actual
$$\mathbf d = [-c_x,\; -c_y,\; -c_z,\; c_x,\; c_y,\; c_z]^\top, \qquad \|\mathbf d\|^2 = 2,$$
con los cosenos evaluados en configuración **corriente**.

### 8. Matriz B corotacional
$$\mathbf B = \frac{1}{L_0}\,\mathbf d^\top, \qquad \varepsilon = \mathbf B\,\mathbf u_e.$$

### 9. Rigidez tangente
**Material**:
$$\mathbf K_M = \frac{E_t\,A}{L_0}\,\mathbf d\,\mathbf d^\top \quad (6\times 6,\ \text{rango 1}).$$

**Geométrica**: en 3D la dirección perpendicular al eje es un plano bidimensional. Se usa el proyector perpendicular $\mathbf P_3 = \mathbf I_3 - \hat{\mathbf e}\hat{\mathbf e}^\top$ ($3\times 3$, rango 2):
$$\mathbf K_G = \frac{N}{l}\,\begin{bmatrix}\mathbf P_3 & -\mathbf P_3 \\ -\mathbf P_3 & \mathbf P_3\end{bmatrix} \quad (6\times 6,\ \text{rango 2}), \qquad N = \sigma A.$$

**Caso destensado** (material unilateral con $\varepsilon \leq 0$): el material devuelve $\sigma = 0$, $E_t = 0$, luego $N = 0$ y **ambas contribuciones se anulan** — $\mathbf K_T = \mathbf 0$. El elemento aporta cero al sistema global.

### 10. Fuerzas internas
$$\mathbf F_{\text{int}} = N\,\mathbf d.$$

### 11. Cuadratura
No aplica (forma cerrada).

---

## Contrato de implementación

```yaml
name: Cable3DCorot
kind: element
status: validated

interface:
  dof_names: [ux, uy, uz]
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
  configuration: "Updated Lagrangian — l, c_x, c_y, c_z se recalculan en cada evaluación"

validity:
  - "|ε| ≲ 1e-2 en régimen tensado"
  - rotaciones y desplazamientos de cualquier magnitud en el espacio
  - cargas axiales (las transversales se transmiten por destensado-retensado del cable)

out_of_scope:
  - flexión, cortante, momentos, torsión
  - dinámica con inercia distribuida (sin matriz de masa)
  - pandeo (irrelevante: el cable no transmite compresión)
  - "fuerzas de cuerpo distribuidas: el usuario reparte masa/peso a los nodos"
  - pretensión mediante parámetro del elemento

numerical_caveats:
  - "Cable destensado ⇒ K_T = 0 en los DOFs del elemento. Si otros elementos no garantizan rigidez en esos DOFs, la matriz tangente global puede ser singular."
  - "En transiciones tensado ↔ destensado (ε cruza 0), Newton-Raphson puede oscilar. Usar arc-length o pasos finos si el problema atraviesa destensiones."

acceptance:
  - name: cable tensado coincide con barra corotacional 3D
    setup: "cable con orientación arbitraria, material lineal (Elastic1D), desplazamiento axial pequeño (ε ~ 1e-6)"
    expect: "K_T y F_int coinciden con Truss3DCorot a tolerancia rtol = 1e-4"
    tol_rel: 1.0e-4

  - name: cable destensado produce aportación nula
    setup: "cable sobre eje x con CableMaterial1D, u_x nodo 2 = -1e-3 ⇒ ε < 0"
    expect: "σ = 0, F_int = 0, K_T = matriz cero"
    tol_abs: 1.0e-12

  - name: invariancia bajo rotación rígida 3D
    setup: "cable inicial sobre eje x, rotación rígida de 45° en plano x-z"
    expect: "σ = 0 y F_int = 0 (independiente del material)"
    tol_abs: 1.0e-10

  - name: cruce por cero
    setup: "u_e = 0 con material unilateral"
    expect: "σ = 0, E_t = 0, F_int = 0, K_T = 0"
    tol_abs: 1.0e-12

references:
  - "Crisfield M.A., Non-linear Finite Element Analysis of Solids and Structures, vol.1, §3.3"
  - "Belytschko T., Nonlinear Finite Elements for Continua and Structures, §4.5"
  - "Irvine H.M., Cable Structures, MIT Press, §1.2"
```

---

## Implementación

- **Archivo**: [fenix/elements/cable.py](../../fenix/elements/cable.py) — comparte archivo con `Cable2DCorot` por convención temática del proyecto; no comparte código.
- **Clase**: `Cable3DCorot` — hereda directamente de `Element` (no de `Cable2DCorot`, no de `Truss3DCorot`, no de ningún otro elemento). Maquinaria cinemática 3D autónoma.
- **Tests**: [tests/test_cable_elements.py](../../tests/test_cable_elements.py) · `TestCable3DCorot` — los 4 criterios de `acceptance` más un test de registro:
  - `test_acceptance_cable_tensado_como_barra_corot_3d` (criterio 1).
  - `test_acceptance_cable_destensado_aportacion_nula` (criterio 2).
  - `test_acceptance_rotacion_rigida_3d` (criterio 3).
  - `test_acceptance_cruce_por_cero` (criterio 4).
  - `test_registro_en_registry` — autodiscover.
- **Notas de traducción**:
  - La clase duplica la maquinaria corotacional 3D de `Truss3DCorot`. Duplicación deliberada por la misma razón que en 2D: desacoplar la evolución futura de armaduras y cables.
  - El proyector $\mathbf P_3 = \mathbf I_3 - \hat{\mathbf e}\hat{\mathbf e}^\top$ se construye con `fenix.math.geometry.perpendicular_projector(ê)`, helper compartido con `Truss3DCorot` (operación geométrica pura, ortogonal al hecho de que el material sea unilateral o no).
  - `_current_geometry` tolera nodos con 2 ó 3 coordenadas (completa con $z=0$ si falta), heredando la flexibilidad de `Truss3D` para problemas embebidos.
  - En estado destensado, $\mathbf K_M = \mathbf 0$ y $\mathbf K_G = \mathbf 0$ simultáneamente; el elemento aporta literalmente ceros al ensamblaje.
  - El test del criterio 1 compara directamente contra `Truss3DCorot` con material `Elastic1D`: si ambos producen $\mathbf K$ y $\mathbf F_{\text{int}}$ iguales a tolerancia, la cinemática del cable y de la barra son consistentes entre sí en el régimen tensado.

---

## Diálogo

- **2026-04-21** · Elemento creado como pieza totalmente autónoma por decisión del usuario: no hereda de `Cable2DCorot` ni de `Truss3DCorot`. Compartir archivo con `Cable2DCorot` (`fenix/elements/cable.py`) es convención temática del proyecto (elementos del mismo dominio conviven en un archivo, como las 4 armaduras y los 2 frames 2D en `structural.py`), no implica dependencia de código.
- **2026-04-21** · La única diferencia estructural con `Cable2DCorot` es el proyector 3D: en 2D la perpendicular es un vector único ($\mathbf n$); en 3D es un plano (dimensión 2). Esto se refleja en el rango de $\mathbf K_G$: 1 en 2D, 2 en 3D.
