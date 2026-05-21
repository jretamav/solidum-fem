# ELEMENTO MARCO/VIGA 2D EULER-BERNOULLI

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Elemento 1D inmerso en el plano 2D que modela una **viga esbelta** según la teoría de **Euler-Bernoulli**. Dos nodos rígidamente conectados (transmiten momento), 3 grados de libertad por nodo ($u_x, u_y, r_z$). Captura tres acciones internas: **fuerza axial**, **cortante transverso** y **momento flector**. Régimen de **linealidad geométrica** (pequeños desplazamientos y rotaciones).

### 1. Hipótesis cinemática de Euler-Bernoulli
- Las secciones transversales se mantienen **planas** tras la deformación.
- Las secciones permanecen **perpendiculares al eje neutro deformado** (cizalladura transversal nula: $\gamma_{xy} = 0$).

Consecuencia: la rotación de la sección es igual a la pendiente de la elástica:
$$\theta(s) = \frac{dv}{ds}.$$

Válido solo para vigas **esbeltas** ($L/h \gtrsim 10$). En vigas peraltadas la cizalladura no es despreciable → `Frame2DTimoshenko`.

### 2. Campos de desplazamiento
- Axial: $u(s)$ a lo largo del eje local $s \in [0, L]$.
- Transverso: $v(s)$ perpendicular al eje local.
- Rotación de la sección: $\theta(s) = dv/ds$ (derivada de la flecha).

### 3. Medidas de deformación
- Axial: $\varepsilon_{axial}(s) = du/ds$ (constante para interpolación lineal en $u$).
- Curvatura: $\kappa(s) = d^2v/ds^2$.

La deformación axial en un punto genérico de la sección es $\varepsilon(s, y) = \varepsilon_{axial}(s) - y\,\kappa(s)$; al integrar sobre la sección se obtienen axial ($N$) y momento flector ($M$) desacoplados.

### 4. Ecuaciones constitutivas
- Axial: $N = EA\,\varepsilon_{axial}$ (o $N = A\,\sigma$ con $\sigma$ del material).
- Flexión: $M = EI\,\kappa$.

El material delega la respuesta axial; la rigidez a flexión escala con el **módulo tangente** $E_t$ devuelto por el material (aproximación: toda la sección entra en régimen no-elástico simultáneamente — no hay plasticidad distribuida).

### 5. Equilibrio — forma débil
Principio de trabajos virtuales. El elemento calcula solo el lado interno:
$$\delta W_{\text{int}} = \int_0^L (N\,\delta\varepsilon_{axial} + M\,\delta\kappa)\,ds = \delta\mathbf u_e^\top\,\mathbf F_{\text{int}}.$$

---

## Formulación numérica (FEM)

### 6. Discretización
Dos nodos ($s = 0$ y $s = L$), 3 DOFs por nodo en coordenadas globales:
$$\mathbf u_e = [u_x^{(1)}, u_y^{(1)}, r_z^{(1)}, u_x^{(2)}, u_y^{(2)}, r_z^{(2)}]^\top.$$

### 7. Sistema local y transformación
Cosenos directores del eje (configuración inicial, fija): $c = (x_2 - x_1)/L$, $s_\theta = (y_2 - y_1)/L$. Matriz de transformación ortogonal $6\times 6$:
$$\mathbf T = \begin{bmatrix}
 c & s_\theta & 0 & 0 & 0 & 0 \\
-s_\theta & c & 0 & 0 & 0 & 0 \\
 0 & 0 & 1 & 0 & 0 & 0 \\
 0 & 0 & 0 & c & s_\theta & 0 \\
 0 & 0 & 0 & -s_\theta & c & 0 \\
 0 & 0 & 0 & 0 & 0 & 1
\end{bmatrix}.$$

Las rotaciones $r_z$ son invariantes ante rotación del sistema (escalares en el plano). $\mathbf u_{\text{local}} = \mathbf T\,\mathbf u_e$.

### 8. Funciones de forma
- Axial: lineal, $N_1^u(s) = 1 - s/L$, $N_2^u(s) = s/L$.
- Transverso (Hermite cúbicos, cuatro funciones que interpolan $v$ y $\theta$):
$$N_1^v = 1 - 3\xi^2 + 2\xi^3, \quad N_1^\theta = L(\xi - 2\xi^2 + \xi^3),$$
$$N_2^v = 3\xi^2 - 2\xi^3, \quad N_2^\theta = L(-\xi^2 + \xi^3),$$
con $\xi = s/L$.

### 9. Rigidez local
En el sistema local, la matriz $6\times 6$ tiene forma cerrada por integración analítica:

$$\mathbf K_{\text{local}} = \begin{bmatrix}
 \tfrac{EA}{L} & 0 & 0 & -\tfrac{EA}{L} & 0 & 0 \\
 0 & \tfrac{12EI}{L^3} & \tfrac{6EI}{L^2} & 0 & -\tfrac{12EI}{L^3} & \tfrac{6EI}{L^2} \\
 0 & \tfrac{6EI}{L^2} & \tfrac{4EI}{L} & 0 & -\tfrac{6EI}{L^2} & \tfrac{2EI}{L} \\
-\tfrac{EA}{L} & 0 & 0 & \tfrac{EA}{L} & 0 & 0 \\
 0 & -\tfrac{12EI}{L^3} & -\tfrac{6EI}{L^2} & 0 & \tfrac{12EI}{L^3} & -\tfrac{6EI}{L^2} \\
 0 & \tfrac{6EI}{L^2} & \tfrac{2EI}{L} & 0 & -\tfrac{6EI}{L^2} & \tfrac{4EI}{L}
\end{bmatrix}.$$

En no-linealidad material, $E \to E_t(\varepsilon_{axial})$ devuelto por el material — escala **todos** los términos de la matriz por igual. Esta es una aproximación: supone que la plasticidad/daño afecta uniformemente a la sección.

### 10. Fuerzas internas
En el sistema local:
$$\mathbf F_{\text{int,local}} = \mathbf K_{\text{local}}\,\mathbf u_{\text{local}}.$$

La componente axial se **corrige** con el esfuerzo axial verdadero del material: $F_1 = -\sigma A$, $F_4 = +\sigma A$. Esto permite que el material aporte $\sigma$ no-lineal sin que la formulación de flexión interfiera.

### 11. Ensamblaje global
$$\mathbf K_{\text{global}} = \mathbf T^\top\,\mathbf K_{\text{local}}\,\mathbf T, \qquad \mathbf F_{\text{int}} = \mathbf T^\top\,\mathbf F_{\text{int,local}}.$$

### 12. Cuadratura
No aplica (forma cerrada para la matriz; integración analítica de los polinomios de Hermite). El campo `n_integration_points: 1` se refiere al único estado material conceptual del elemento (la deformación axial media).

---

## Contrato de implementación

```yaml
name: Frame2DEuler
kind: element
status: validated

interface:
  dof_names: [ux, uy, rz]
  n_nodes: 2
  strain_dim: 1
  n_integration_points: 1

parameters:
  - { name: A, type: float, required: true, desc: "Área de la sección transversal" }
  - { name: I, type: float, required: true, desc: "Momento de inercia respecto al eje perpendicular al plano (z)" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: "axial escalar (deformación media axial ε = (u_2 - u_1)/L en sistema local)"
  nonlinearity_model: "E_tangent escala toda la matriz local — no hay distinción axial vs flexión en régimen no-elástico"

conventions:
  sign: "σ > 0 ⇔ ε > 0 (elongación) en el eje axial"
  rotation_sign: "r_z positivo antihorario (convención matemática estándar)"
  node_orientation: "el eje local va del nodo 1 al nodo 2; permutar nodos invierte el signo de las fuerzas internas"
  configuration: "inicial fija — L, c, s, T no se actualizan con los desplazamientos"

validity:
  - "vigas esbeltas: L/h ≳ 10 (h = canto de la sección)"
  - pequeños desplazamientos y pequeñas rotaciones
  - "|ε_axial| ≲ 1e-2"

out_of_scope:
  - vigas peraltadas con deformación por cortante significativa (usar Frame2DTimoshenko)
  - grandes rotaciones o desplazamientos (no hay variante corotacional en el catálogo actual)
  - plasticidad distribuida en la sección (fibras); el modelo escala rigidez global por E_t
  - pandeo (requiere rigidez geométrica, no implementada en esta viga lineal)
  - "fuerzas de cuerpo distribuidas: el usuario reparte carga equivalente a los nodos"

acceptance:
  - name: respuesta axial pura
    setup: "viga horizontal, empotrada en extremo izquierdo, carga axial F en extremo derecho"
    expect: "u_x (extremo libre) = F·L/(E·A); momentos y cortantes nulos"
    tol_rel: 1.0e-10

  - name: voladizo con carga transversal (flecha analítica)
    setup: "viga horizontal, empotrada en extremo izquierdo, carga vertical P en extremo derecho"
    expect: "v (extremo libre) = P·L³/(3·E·I); rotación θ = P·L²/(2·E·I)"
    tol_rel: 1.0e-10

  - name: simetría de K
    setup: "cualquier configuración"
    expect: "K_global = K_globalᵀ (elemento conservativo en régimen elástico)"
    tol_abs: 1.0e-12

references:
  - "Bathe K.-J., Finite Element Procedures, §4.2 (formulación isoparamétrica) y cap. 5 (vigas)"
  - "Cook R.D., Concepts and Applications of FEA, §2.7 (elementos de viga)"
```

---

## Implementación

- **Archivo**: [solidum/elements/frame/euler.py](../../solidum/elements/frame/euler.py) — submódulo del paquete `solidum/elements/frame/` que aloja también `Frame2DTimoshenko` y `Frame2DEulerCorot`.
- **Clase**: `Frame2DEuler` — hereda directamente de `Element`. La construcción de la matriz de transformación $\mathbf T$ se delega a `build_geometry_2d` en [solidum/elements/frame/_shared.py](../../solidum/elements/frame/_shared.py), compartida con `Frame2DTimoshenko` (la versión corotacional reconstruye $\mathbf T$ desde `alpha0` por su lógica propia).
- **Tests**: [tests/test_frame.py](../../tests/test_frame.py) · `TestFrame2DEulerAcceptance` — los tres criterios físicos y el registro:
  - `test_acceptance_respuesta_axial_pura` (criterio 1) — resuelve el voladizo con solver, verifica $u_x(L) = FL/(EA)$.
  - `test_acceptance_voladizo_carga_transversal` (criterio 2) — carga transversal $P$, verifica flecha $v = PL^3/(3EI)$ y rotación $\theta = PL^2/(2EI)$.
  - `test_acceptance_simetria_K` (criterio 3) — viga oblicua (c=0.6, s=0.8) para evitar ejes canónicos, verifica $\mathbf K = \mathbf K^\top$.
  - `test_registro_en_registry` — autodiscover.
- **Notas de traducción**:
  - La matriz $\mathbf T$ se calcula una sola vez en `__init__` sobre la configuración inicial y se guarda como atributo — coherente con el régimen geométricamente lineal.
  - En `compute_element_state`, la componente axial de $\mathbf F_{\text{int,local}}$ se sobrescribe tras el producto $\mathbf K_{\text{local}}\mathbf u_{\text{local}}$ con el valor $\sigma A$ devuelto por el material. Esto permite usar materiales no-lineales en la rama axial sin que la parte lineal de flexión interfiera.
  - La aproximación del modelo no-lineal es que $E_t$ escala **toda** la matriz local: no hay distinción entre rigidez axial y de flexión en régimen no-elástico. Documentado en `material_contract.nonlinearity_model`.

---

## Diálogo

- **2026-04-21** · Elemento movido a archivo propio `solidum/elements/frame.py` y desacoplado del helper compartido `_frame_geometry`. El método estático `_build_geometry` reimplementa la construcción de $\mathbf T$ dentro de la clase. `Frame2DTimoshenko` sigue en `structural.py` usando su propio acceso al helper compartido; si mañana también se independiza, duplicará o internalizará su propia versión.
- **2026-04-21** · Tests de aceptación cubren la flecha analítica del voladizo a 8 decimales ($v = PL^3/(3EI)$) usando el solver completo del proyecto — validan no solo la cinemática del elemento sino también su integración con ensamblador y solver.
- **2026-05-13** · `frame.py` se parte en paquete `solidum/elements/frame/{euler,timoshenko,euler_corot,_shared}.py`. La duplicación de `_build_geometry` entre Euler y Timoshenko se elimina: la construcción de $\mathbf T$ vuelve a vivir en un único lugar (`build_geometry_2d` en `_shared.py`), pero esta vez sin función helper externa flotante — pertenece al paquete `frame/`. Frame2DEulerCorot mantiene su propia reconstrucción de $\mathbf T$ desde `alpha0` porque su cinemática corotacional no se beneficia del helper común.
