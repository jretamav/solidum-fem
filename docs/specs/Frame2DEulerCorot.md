# ELEMENTO MARCO/VIGA 2D EULER-BERNOULLI COROTACIONAL

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Viga 2D esbelta según Euler-Bernoulli en formulación **corotacional** (Updated Lagrangian). Dos nodos rígidamente conectados, 3 DOFs por nodo ($u_x, u_y, r_z$). Captura grandes desplazamientos y grandes rotaciones del elemento como un todo con **pequeñas deformaciones axiales** ($|\varepsilon_{axial}| \lesssim 1\%$) y rotaciones nodales deformacionales moderadas.

Abre dominios que la viga lineal no captura: **pandeo por flexión**, **post-pandeo**, **snap-through** de arcos planos, brazos flexibles esbeltos.

### 1. Cinemática corotacional
Se separa el movimiento del elemento en:

1. **Rotación rígida** $\alpha_e$: el ángulo que el eje de la barra ha girado desde la configuración inicial.
2. **Rotación deformacional de cada sección** $\bar\theta_i$: la rotación nodal menos la rotación rígida (es lo que genera momento flector).

Sean:
$$\alpha_0 = \arctan2(Y_2 - Y_1,\; X_2 - X_1) \quad \text{(ángulo del eje en config. inicial, fijo)},$$
$$\alpha = \arctan2(y_2 - y_1,\; x_2 - x_1) \quad \text{(ángulo corriente)},$$
$$\alpha_e = \alpha - \alpha_0 \quad \text{(rotación rígida del elemento)}.$$

Rotaciones deformacionales:
$$\bar\theta_1 = \theta_1 - \alpha_e, \qquad \bar\theta_2 = \theta_2 - \alpha_e,$$
donde $\theta_i = r_z^{(i)}$ son las rotaciones totales de los nodos.

### 2. Medidas de deformación
- Axial: $\varepsilon_{axial} = (l - L_0)/L_0$ con $l$ longitud corriente.
- Curvatura efectiva: $\kappa = (\bar\theta_2 - \bar\theta_1)/L_0$ (constante si se usa interpolación cúbica).

### 3. Ecuaciones constitutivas
- Axial: delegado al material ($\sigma, E_t$).
- Flexión: $M = E_t I\,\kappa$ (el mismo $E_t$ del material escala la rigidez de flexión, como en `Frame2DEuler`).

### 4. Equilibrio — forma débil
PTV en configuración corriente. El elemento calcula solo el lado interno:
$$\delta W_{\text{int}} = \delta\mathbf u_e^\top\,\mathbf F_{\text{int}}, \qquad \mathbf K_T = \mathbf K_{\text{material}} + \mathbf K_{\sigma}.$$

---

## Formulación numérica (FEM)

### 5. DOFs locales deformacionales
Se reduce la descripción del estado deformacional del elemento a **3 magnitudes escalares**:
$$\mathbf q_l = [\,u_l,\; \bar\theta_1,\; \bar\theta_2\,]^\top, \qquad u_l = l - L_0.$$

### 6. Rigidez local
En el sistema de coordenadas corrotado (eje de la barra corriente), la rigidez local $3\times 3$:
$$\mathbf K_{ll} = \begin{bmatrix}
\tfrac{E_t A}{L_0} & 0 & 0 \\
0 & \tfrac{4 E_t I}{L_0} & \tfrac{2 E_t I}{L_0} \\
0 & \tfrac{2 E_t I}{L_0} & \tfrac{4 E_t I}{L_0}
\end{bmatrix}.$$

Fuerzas locales $\mathbf F_l = \mathbf K_{ll}\,\mathbf q_l = [N,\; M_1,\; M_2]^\top$, con $N$ sustituido por $\sigma A$ directamente del material para permitir respuestas no-lineales.

### 7. Matriz de transformación $\mathbf B$ (3×6)
Relación entre variaciones de $\mathbf q_l$ y de los DOFs globales $\delta\mathbf u_e$:
$$\mathbf B = \begin{bmatrix}
-c & -s & 0 & c & s & 0 \\
-s/l & c/l & 1 & s/l & -c/l & 0 \\
-s/l & c/l & 0 & s/l & -c/l & 1
\end{bmatrix}, \qquad c = \cos\alpha,\; s = \sin\alpha.$$

$c, s$ y $l$ se evalúan en **configuración corriente** — esta es la esencia corotacional.

Fuerzas internas globales:
$$\mathbf F_{\text{int}} = \mathbf B^\top\,\mathbf F_l.$$

### 8. Rigidez tangente
Se descompone en contribución material más geométrica:

$$\mathbf K_T = \mathbf K_{\text{material}} + \mathbf K_\sigma,$$
$$\mathbf K_{\text{material}} = \mathbf B^\top\,\mathbf K_{ll}\,\mathbf B,$$
$$\mathbf K_\sigma = \frac{N}{l}\,\mathbf z\,\mathbf z^\top + \frac{M_1 + M_2}{l^2}\,\bigl(\mathbf r\,\mathbf z^\top + \mathbf z\,\mathbf r^\top\bigr),$$

con vectores geométricos de 6 componentes:
$$\mathbf r = [-c,\;-s,\;0,\; c,\; s,\;0]^\top \quad \text{(dirección axial)},$$
$$\mathbf z = [-s,\; c,\;0,\; s,\;-c,\;0]^\top \quad \text{(perpendicular)}.$$

La parte $(N/l)\mathbf z\mathbf z^\top$ es análoga al $\mathbf K_G$ de `Truss2DCorot`: captura la rigidización por tensión (y el ablandamiento precrítico por compresión, base del pandeo). La parte con momentos acopla rotaciones con traslaciones — crítica para problemas de flexión con desplazamientos significativos.

### 9. Cuadratura
No aplica (forma cerrada con Hermite cúbicos).

---

## Contrato de implementación

```yaml
name: Frame2DEulerCorot
kind: element
status: validated

interface:
  dof_names: [ux, uy, rz]
  n_nodes: 2
  strain_dim: 1
  n_integration_points: 1

parameters:
  - { name: A, type: float, required: true, desc: "Área de la sección transversal" }
  - { name: I, type: float, required: true, desc: "Momento de inercia respecto al eje perpendicular al plano" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: "axial escalar corotacional (ε = (l − L₀)/L₀)"
  nonlinearity_model: "E_tangent escala K_material local; el K_σ geométrico depende solo de N, M₁, M₂ corrientes"

conventions:
  sign: "σ > 0 ⇔ ε > 0 (elongación); r_z positivo antihorario"
  node_orientation: "eje local del nodo 1 al nodo 2"
  configuration: "Updated Lagrangian — l, c, s, α se recalculan en cada evaluación"

validity:
  - "|ε_axial| ≲ 1e-2 (pequeña deformación axial)"
  - rotaciones rígidas del elemento de cualquier magnitud, pero acumuladas en |α_e| < π entre commits (rotaciones continuas > π requerirían tracking de vueltas — fuera de alcance)
  - rotaciones deformacionales de las secciones moderadas (|θ̄_i| ≲ 30°); para valores mayores la interpolación Hermite cúbica pierde precisión
  - vigas esbeltas (L/h ≳ 10)

out_of_scope:
  - rotaciones continuas del eje > π sin pasar por commit (pendulums; aliasing en atan2)
  - deformaciones axiales grandes (requiere medida logarítmica)
  - plasticidad distribuida en la sección (fibras)
  - cortante transverso significativo (usar versión Timoshenko, pendiente)
  - "fuerzas de cuerpo distribuidas: el usuario reparte carga a los nodos"

numerical_caveats:
  - "El ángulo α_e se unwrappea a (-π, π] por llamada. Si el problema tiene saltos de rotación > π entre iteraciones consecutivas del solver, la formulación produce resultados incorrectos. Usar pasos de carga finos."
  - "En problemas con pandeo, la rigidez tangente se hace singular al cruzar el punto crítico. Usar `ArcLengthSolver` para seguir la rama post-crítica."

acceptance:
  - name: límite lineal coincide con Frame2DEuler
    setup: "voladizo horizontal, cargas pequeñas (ε ~ 1e-8, v/L ~ 1e-6)"
    expect: "K_T y F_int coinciden con Frame2DEuler lineal a tolerancia rtol=1e-4"
    tol_rel: 1.0e-4

  - name: invariancia bajo rotación rígida
    setup: "viga inicialmente horizontal, rotada rígidamente 45° (ambos nodos)"
    expect: "F_int = 0, σ = 0 (ningún esfuerzo interno bajo rotación rígida pura)"
    tol_abs: 1.0e-10

  - name: coherencia tangente/fuerza interna (chequeo por diferencias finitas)
    setup: "configuración arbitraria no trivial; K_T analítica vs derivada numérica de F_int"
    expect: "‖K_T_analítica − K_T_numérica‖ / ‖K_T_analítica‖ ≲ 1e-5"
    tol_rel: 1.0e-5

  - name: simetría de K
    setup: "configuración arbitraria"
    expect: "K_T = K_Tᵀ"
    tol_abs: 1.0e-6

references:
  - "Crisfield M.A., Non-linear Finite Element Analysis of Solids and Structures, vol.1, §7.3 (corotational beams)"
  - "Belytschko T., Nonlinear Finite Elements for Continua and Structures, §4.11"
  - "Wriggers P., Nonlinear Finite Element Methods, cap. 4"
```

---

## Implementación

- **Archivo**: [fenix/elements/frame.py](../../fenix/elements/frame.py) — convive con `Frame2DEuler` y `Frame2DTimoshenko` por familia temática, sin heredar de ninguno.
- **Clase**: `Frame2DEulerCorot` — hereda directamente de `Element`. Métodos internos `_current_geometry` (longitud, cosenos, ángulo corriente y rigid-body $\alpha_e$) y `_unwrap` (reducción a $(-\pi, \pi]$).
- **Tests**: [tests/test_frame.py](../../tests/test_frame.py) · `TestFrame2DEulerCorotAcceptance` — 4 criterios + registro:
  - `test_acceptance_limite_lineal_coincide_con_frame2deuler`
  - `test_acceptance_invariancia_rotacion_rigida`
  - `test_acceptance_coherencia_tangente_fuerza_interna` — **finite-difference check** de $\mathbf K_T$ contra derivada numérica de $\mathbf F_{\text{int}}$, red de seguridad clave en esta formulación.
  - `test_acceptance_simetria_K`
  - `test_registro_en_registry`
- **Notas de traducción**:
  - La formulación se basa en 3 DOFs deformacionales $[u, \bar\theta_1, \bar\theta_2]$ y una matriz $\mathbf B$ $3\times 6$ que los relaciona con los DOFs globales. Evaluada en configuración corriente.
  - **Signo de la contribución de momentos en $\mathbf K_\sigma$**: tras derivar $\partial(\mathbf z/l)/\partial\mathbf u_e$ analíticamente se obtiene $-(1/l^2)(\mathbf r\mathbf z^\top + \mathbf z\mathbf r^\top)$. El signo es **negativo** — un error común es ponerlo positivo. El test de diferencias finitas lo cazó en la primera implementación.
  - `compute_internal_forces` proyecta $\mathbf F_{\text{int}}$ al sistema local corriente para reportar axial/cortante separados, utilizando la submatriz $2\times 2$ de rotación del eje corriente.
  - `_unwrap` lleva los ángulos al intervalo $(-\pi, \pi]$; la formulación es válida mientras la rotación acumulada por paso sea menor que $\pi$.

---

## Diálogo

- **2026-04-21** · Implementación siguiendo Crisfield §7.3. Convive en `frame.py` con las vigas lineales por familia temática (todas son vigas 2D), sin herencia ni helpers compartidos.
- **2026-04-21** · La primera versión tenía el signo equivocado en la parte de momentos de $\mathbf K_\sigma$ ($+$ en lugar de $-$). El test de coherencia tangente/fuerza interna por diferencias finitas lo detectó inmediatamente: $\|\mathbf K_{T,\text{analítica}} - \mathbf K_{T,\text{FD}}\|_{\text{rel}}$ del orden de 1e-1 en lugar de 1e-5. Tras derivar $\partial(\mathbf z/l)/\partial\mathbf u_e$ cuidadosamente (ver notas) el signo correcto es negativo. Este episodio justifica retrospectivamente incluir un test FD como criterio de aceptación en toda formulación no-lineal geométrica futura.
- **2026-04-21** · El test del límite lineal usa $P=10$ N para obtener $v/L \sim 10^{-4}$: suficientemente pequeño para estar en régimen lineal pero suficientemente grande para que el residuo del solver converja por debajo de `tol=1e-8` (cargas menores producen desplazamientos del orden del ruido numérico y el Newton oscila en la tolerancia relativa).
