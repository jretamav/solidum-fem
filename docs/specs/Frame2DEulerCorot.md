# ELEMENTO MARCO/VIGA 2D EULER-BERNOULLI COROTACIONAL

> Orden de trabajo. El usuario escribe **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.

---

## Especificación física

### 0. Descripción general
Viga 2D esbelta según Euler-Bernoulli en formulación **corotacional** (Updated Lagrangian). Dos nodos rígidamente conectados, 3 DOFs por nodo ($u_x, u_y, r_z$). Captura grandes desplazamientos y grandes rotaciones del elemento como un todo con **pequeñas deformaciones axiales** ($|\varepsilon_{axial}| \lesssim 1$ %) y rotaciones nodales deformacionales moderadas.

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
La ley del elemento vive en el sistema básico $[\Delta l,\ \bar\theta_1,\ \bar\theta_2]$ con fuerzas básicas $[N,\ M_1,\ M_2]$. Dos caminos excluyentes, elegidos al construir ([ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md)):

- **Con `material` + `A`, `I`** (forma cerrada): $N = EA\,\varepsilon_{axial}$ y $M = EI\,\kappa$. Sólo materiales **elásticos lineales** (`IS_LINEAR_ELASTIC = True`); un material inelástico se rechaza al construir con un mensaje que remite a la sección por fibras.
- **Con `section`** (+ `n_ip`): la ley de sección ([FiberSection](FiberSection.md)) se integra en los puntos de Gauss del elemento en el sistema básico corrotado (§6).

Hasta el 2026-09-29 el camino con material admitía materiales inelásticos y escalaba también la rigidez de flexión con la $E_t$ axial, como `Frame2DEuler`; esa ley desapareció (ver *Diálogo*).

### 4. Equilibrio — forma débil
PTV en configuración corriente. El elemento calcula solo el lado interno:
$$\delta W_{\text{int}} = \delta\mathbf u_e^\top\,\mathbf F_{\text{int}}, \qquad \mathbf K_T = \mathbf K_{\text{material}} + \mathbf K_{\sigma}.$$

---

## Formulación numérica (FEM)

### 5. DOFs locales deformacionales
Se reduce la descripción del estado deformacional del elemento a **3 magnitudes escalares**:
$$\mathbf q_l = [\,u_l,\; \bar\theta_1,\; \bar\theta_2\,]^\top, \qquad u_l = l - L_0.$$

### 6. Rigidez local
En el sistema de coordenadas corrotado (eje de la barra corriente), camino con material elástico lineal, la rigidez local $3\times 3$:
$$\mathbf K_{ll} = \begin{bmatrix}
\tfrac{E A}{L_0} & 0 & 0 \\
0 & \tfrac{4 E I}{L_0} & \tfrac{2 E I}{L_0} \\
0 & \tfrac{2 E I}{L_0} & \tfrac{4 E I}{L_0}
\end{bmatrix}.$$

Fuerzas locales $\mathbf F_l = \mathbf K_{ll}\,\mathbf q_l = [N,\; M_1,\; M_2]^\top$, con $N = \sigma A$ del material (con material elástico lineal coincide con la fila axial del producto).

Camino con sección: en cada punto de Gauss $x_g = \xi_g L_0$,
$$\mathbf e_g = \begin{bmatrix}\varepsilon_0\\ \kappa_z\end{bmatrix} = \mathbf B_g\,\mathbf q_l, \qquad \mathbf B_g = \frac{1}{L_0}\begin{bmatrix}1 & 0 & 0\\ 0 & 6\xi_g - 4 & 6\xi_g - 2\end{bmatrix},$$
y $\mathbf F_l = \sum_g w_g L_0\,\mathbf B_g^\top\,\mathbf s_g$, $\mathbf K_{ll} = \sum_g w_g L_0\,\mathbf B_g^\top\,\mathbf k_{s,g}\,\mathbf B_g$ (pesos de Gauss-Legendre sobre $(0, 1)$). Con fibras elásticas y `n_ip` ≥ 2 reproduce la matriz anterior exactamente. La cinemática corrotacional de §7-§8, incluida $\mathbf K_\sigma$ con $N$, $M_1$ y $M_2$, es la misma en los dos caminos.

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
$$\mathbf K_\sigma = \frac{N}{l}\,\mathbf z\,\mathbf z^\top - \frac{M_1 + M_2}{l^2}\,\bigl(\mathbf r\,\mathbf z^\top + \mathbf z\,\mathbf r^\top\bigr),$$

con vectores geométricos de 6 componentes:
$$\mathbf r = [-c,\;-s,\;0,\; c,\; s,\;0]^\top \quad \text{(dirección axial)},$$
$$\mathbf z = [-s,\; c,\;0,\; s,\;-c,\;0]^\top \quad \text{(perpendicular)}.$$

El signo **menos** del término de momentos va con este $\mathbf z$ (con $\mathbf z \to -\mathbf z$, como en Crisfield, el término cambia de signo y $\mathbf z\mathbf z^\top$ no). Errata corregida el 2026-11-27: la spec decía «+»; el código (`euler_corot.py`) siempre usó «−», y la tangente frente a diferencias finitas de $\mathbf F_{\text{int}}$ lo confirma (1.8e-10 con «−», 1.1e-3 con «+»; `teorico-solidum`, `informes/teorico/consulta_perturbacion_corot_2026-11-27.md`).

La parte $(N/l)\mathbf z\mathbf z^\top$ es análoga al $\mathbf K_G$ de `Truss2DCorot`: captura la rigidización por tracción (y el ablandamiento precrítico por compresión, base del pandeo). La parte con momentos acopla rotaciones con traslaciones — crítica para problemas de flexión con desplazamientos significativos.

### 9. Cuadratura
Camino con material: no aplica (forma cerrada con Hermite cúbicos). Camino con sección: Gauss-Legendre de `n_ip` puntos (3 por defecto, mínimo 2), con `N_INTEGRATION_POINTS = n_ip` como atributo de instancia.

### 10. Matriz de masa
La de [Frame2DEuler](Frame2DEuler.md) §13, evaluada en la configuración de referencia ($L_0$, $\alpha_0$; Lagrangeano total): **sin inercia rotatoria** de la sección (Euler-Bernoulli), consistente $\tfrac{\rho A L_0}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$ axial más la hermítica $\tfrac{\rho A L_0}{420}\mathbf H$ en $(v_1, \theta_1, v_2, \theta_2)$, y concentrada HRZ $\rho A L_0/2$ traslacional y $\rho A L_0^3/78$ rotacional por nodo. Depende sólo de $\rho A$ ($\sum_i \rho_i A_i$ con sección). Hasta el 2026-09-29 llevaba además la inercia rotatoria $\rho I$ (ver *Diálogo*).

### 11. Carga por unidad de longitud y fuerzas internas a lo largo del elemento
Las de [Frame2DEuler](Frame2DEuler.md) §14–§15 ([ADR 0044](../adr/0044-cargas-de-linea-y-diagramas-de-fuerzas-internas.md)), con dos precisiones propias de la cinemática corrotacional:

- **Carga muerta en la configuración de referencia.** El vector consistente se calcula una vez con $\alpha_0$ y $L_0$, como la fuerza de cuerpo: una carga en ejes locales se refiere a los ejes **iniciales** del elemento y no gira con él. Con giros grandes del elemento, los momentos nodales de la carga ($\pm q_\perp L_0^2/12$) quedan con la componente transversal de la configuración inicial; el error de la solución es de orden 2 en el tamaño del elemento (medido por el teórico el 2026-11-27: órdenes 1.97–2.00 en un voladizo con giros de 0.6 a 1.3 rad, y la misma solución límite que con un vector consistente con la cinemática corrotacional, que no converge mejor y exige una rigidez de carga). Es lo que hacen OpenSees en 2D (`CorotCrdTransf2d`: las cargas del elemento se suman en los ejes locales iniciales; leído en el código, no documentado) y Code_Aster (`FORCE_POUTRE` sólo admite `FIXE_CSTE`, carga sobre la geometría inicial); Abaqus, en cambio, trata las cargas de viga en ejes locales (`P1`, `P2`) como seguidoras, con rigidez de carga, y las globales como muertas (`informes/bibliotecario/busqueda_cargas_linea_2026-11-27.md`, P2).
- **Diagramas en la cuerda corriente.** Las fuerzas de extremo se proyectan sobre los ejes de la cuerda corriente ($\alpha$, como `internal_forces`) y la carga, fija en el espacio, también. Con $\xi = s_0/L_0 \in [0, 1]$ la coordenada material y $l$ la longitud corriente de la cuerda, la carga se integra por unidad de longitud de referencia y los brazos se miden sobre la cuerda corriente:
$$Q(\xi) = L_0\!\int_0^\xi \tilde q'\,d\xi', \qquad N = N_i - Q_{x'}, \quad V = V_i + Q_{y'}, \quad M = M_i + \xi\,l\,V_i + L_0\,l\!\int_0^\xi (\xi - \xi')\,\tilde q_{y'}\,d\xi'.$$
  $N$ y $V$ cierran siempre en $\xi = 1$, y $M$ cierra exactamente con carga uniforme a cualquier giro y alargamiento. Con carga variable queda, para una carga lineal, $\Delta M = \Delta q_\perp L_0^2\,(\lambda\cos\varphi - 1)/60$ ($\Delta q_\perp$: variación de la componente transversal inicial entre los extremos; $\varphi$: giro de la cuerda; $\lambda = l/L_0$), de segundo orden en $\varphi$ y $O(L_0^3)$ por elemento: se anula al refinar.
- **Límite dentro del elemento.** El diagrama es estática sobre la cuerda: no incluye el momento de la fuerza axial por la separación entre el eje deformado y la cuerda ($N\,v_l$), y $N$ y $V$ se proyectan sobre la cuerda, no sobre la normal de cada sección. Es exacto en los nodos (salvo el cierre anterior) y converge al refinar; con efectos de segundo orden dentro de la barra hacen falta varios elementos por barra.
- **Diagramas de una perturbación lineal** (ADR 0031; decisión del usuario del 2026-11-27, opción (a) de la consulta al teórico). El resultado es el total linealizado alrededor del estado base $u_b$: fuerzas de extremo $F_\text{int}(u_b) + K_T(u_b)\,\Delta U - f$ proyectadas con los ejes de la cuerda **del estado base**. El diagrama es la linealización de la estática del tramo de arriba: se expresa en esos mismos ejes, con la longitud $l_b$ de la cuerda del estado base ($s = \xi\,l_b$), y al flector se le suma el término geométrico de primer orden
$$\Delta M_\text{geo}(\xi) = \Delta l\,\Bigl[\xi\,V_{i,b} + L_0\!\int_0^\xi (\xi - \xi')\,\tilde q_{y',b}\,d\xi'\Bigr] + l_b\,\Delta\alpha\,\Bigl[\xi\,N_{i,b} - L_0\!\int_0^\xi (\xi - \xi')\,\tilde q_{x',b}\,d\xi'\Bigr],$$
  con $\Delta l = \mathbf e_{x',b}\cdot\Delta\mathbf d$, $\Delta\alpha = \mathbf e_{y',b}\cdot\Delta\mathbf d / l_b$ y $\Delta\mathbf d = \Delta\mathbf u_j - \Delta\mathbf u_i$ (traslaciones nodales del incremento), y $N_{i,b}$, $V_{i,b}$ y $\tilde q_b$ las fuerzas de extremo y la carga **del estado base**. El primer sumando es el alargamiento de los brazos, $(\Delta l/l_b)\,[M_b(\xi) - M_{i,b}]$; el segundo, el efecto P-Δ del tramo. Sale de linealizar $M(\xi) = M_i + \xi\,(\mathbf d\times\mathbf p_i)_z + L_0\int_0^\xi(\xi-\xi')\,(\mathbf d\times\mathbf q)_z\,d\xi'$ ($\mathbf d$: cuerda corriente; $\mathbf p_i$: fuerza del nodo $i$ sobre el elemento; $\mathbf p_i\cdot\mathbf e_{x'} = -N_i$, $\mathbf p_i\cdot\mathbf e_{y'} = V_i$ en la convención pública): $\delta\mathbf d\times\mathbf p_{i,b} = \Delta l\,V_{i,b} + l_b\,\Delta\alpha\,N_{i,b}$ y $\delta\mathbf d\times\mathbf q_b = \Delta l\,q_{y',b} - l_b\,\Delta\alpha\,q_{x',b}$. $N$ y $V$ no llevan término geométrico: son las componentes sobre la cuerda del estado base, exactas al primer orden en esos ejes. Así el diagrama cierra con las fuerzas de extremo, es lineal en $\Delta U$ (superponible, como pide el ADR 0031) y su flector difiere del de un análisis no lineal en segundo orden (*calculado* por el teórico: 0.22 % frente al 2.3 % sin el término, con $\Delta\alpha \approx 0.03$; `informes/teorico/consulta_perturbacion_corot_2026-11-27.md`). $N$ y $V$, en cambio, difieren del no lineal en primer orden, $N_{i,b}\,\Delta\alpha$ y $V_{i,b}\,\Delta\alpha$, porque éste los da sobre la cuerda que gira; si se reportan así es una segunda pregunta al usuario, pendiente. Con carga variable en el estado base persiste el cierre de $M$ del punto anterior (de segundo orden). Un elemento geométricamente lineal no tiene término geométrico y sus ejes no dependen del estado: su diagrama de una perturbación es el de siempre.

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
  n_integration_points: 1   # camino con material; con section, N_INTEGRATION_POINTS = n_ip (atributo de instancia)

admission:                # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - de_borst_crisfield: "§9.3, elementos corrotacionales (formulación base, material elástico; 1.ª ed.: Crisfield 1991, §7.3)"
  general_codes:
    - opensees: "elasticBeamColumn + geomTransf Corotational"

parameters:
  - { name: A, type: float, required: false, desc: "Área de la sección transversal. Obligatoria con material; no se declara con section" }
  - { name: I, type: float, required: false, desc: "Momento de inercia respecto al eje perpendicular al plano. Obligatorio con material; no se declara con section" }
  - { name: section, type: "id de sections", required: false, desc: "Sección de marco (FiberSection), evaluada en el sistema básico corrotado. Excluyente con material, A e I (ADR 0021)" }
  - { name: n_ip, type: int, required: false, default: 3, desc: "Puntos de Gauss-Legendre a lo largo del elemento, ≥ 2. Sólo con section" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: "axial escalar corotacional (ε = (l − L₀)/L₀)"
  nonlinearity_model: "camino con material: sólo materiales elásticos lineales (IS_LINEAR_ELASTIC = True; un material inelástico se rechaza al construir). La no linealidad material entra por section (FiberSection) en el sistema básico; el K_σ geométrico depende sólo de N, M₁, M₂ corrientes en los dos caminos"

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
  - material inelástico con material + A, I (se rechaza al construir; la plasticidad y el daño entran por section, FiberSection)
  - cortante transverso significativo (usar versión Timoshenko, pendiente)
  - "fuerzas de cuerpo sobre la geometría corriente (se reparten en la de referencia, constantes o como campo b(x))"

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

  - name: Bathe-Bolourchi — cuarto de círculo                 # añadido 2026-05-19
    setup: "voladizo recto, 10 elementos, momento puro M = (π/2)·EI/L en la punta, 8 pasos"
    expect: "tip se mueve a (R·sin(π/2)-L, R·(1-cos(π/2))) con R=EI/M; θ_tip = π/2"
    tol_abs: 0.005   # sobre L (cuerda ~0.2% en 10 elementos rectos)

  - name: Bathe-Bolourchi — medio círculo                     # añadido 2026-05-19
    setup: "voladizo recto, 20 elementos, momento puro M = π·EI/L en la punta, 20 pasos"
    expect: "tip a (-L, 2L/π) y θ_tip = π — la viga se enrosca en semicírculo cerrado"
    tol_abs: 0.01    # sobre L

  - name: guarda de material elástico lineal                  # añadido 2026-09-29 (ADR 0021)
    setup: "material Elastoplastic1D o IsotropicDamage1D con A e I"
    expect: "ValueError al construir, con un mensaje que nombra IS_LINEAR_ELASTIC y remite a FiberSection con 'section'"

  - name: camino con sección, giros grandes                   # añadido 2026-09-29 (ADR 0021)
    setup: "section = FiberSection elástica frente a material + A, I integrados; giros de sólido rígido grandes y arco de Bathe-Bolourchi (caso 8 de FiberSection)"
    expect: "K_T y F_int iguales; misma solución del arco"
    tol_rel: 1.0e-10

  - name: matrices de masa                                    # añadido 2026-09-29
    setup: "elemento oblicuo, caminos con material y con sección"
    expect: "las de Frame2DEuler §13 en la configuración de referencia: consistente axial + hermítica sin inercia rotatoria, concentrada HRZ diag(ρAL₀/2, ρAL₀/2, ρAL₀³/78); ½u̇ᵀMu̇ exacta en traslación y, con la consistente, en giro rígido alrededor de un extremo"
    tol_rel: 1.0e-12

references:
  - "Crisfield M.A., Non-linear Finite Element Analysis of Solids and Structures, vol.1, §7.3 (corotational beams)"
  - "Belytschko T., Nonlinear Finite Elements for Continua and Structures, §4.11"
  - "Wriggers P., Nonlinear Finite Element Methods, cap. 4"
```

---

## Implementación

- **Archivo**: [solidum/elements/frame/euler_corot.py](../../solidum/elements/frame/euler_corot.py) — convive con `Frame2DEuler` y `Frame2DTimoshenko` en el paquete `solidum/elements/frame/`, sin heredar de ninguno. Comparte con ellos los helpers libres de [_shared.py](../../solidum/elements/frame/_shared.py) (carga de cuerpo, masa consistente, traducción a `ElementForces`), pero no `build_geometry_2d`: reconstruye `T` desde `alpha0` por su lógica corotacional propia.
- **Clase**: `Frame2DEulerCorot(FrameSectionMixin, Element)`. Métodos internos `_current_geometry` (longitud, cosenos, ángulo corriente y rigid-body $\alpha_e$), `_unwrap` (reducción a $(-\pi, \pi]$), `_basic_deformations` y `_basic_response` (fuerzas básicas y $\mathbf K_{ll}$ por el camino con material o por la sección). La constitutiva común de los marcos vive en `FrameSectionMixin` ([solidum/elements/frame/_section.py](../../solidum/elements/frame/_section.py), ADR 0021). `REFERENCE_KWARGS = {"material": MaterialRegistry, "section": SectionRegistry}`.
- **Tests**: [tests/test_frame.py](../../tests/test_frame.py) · `TestFrame2DEulerCorotAcceptance` — 4 criterios + registro:
  - `test_acceptance_limite_lineal_coincide_con_frame2deuler`
  - `test_acceptance_invariancia_rotacion_rigida`
  - `test_acceptance_coherencia_tangente_fuerza_interna` — **finite-difference check** de $\mathbf K_T$ contra derivada numérica de $\mathbf F_{\text{int}}$, red de seguridad clave en esta formulación.
  - `test_acceptance_simetria_K`
  - `test_registro_en_registry`
  - Bathe-Bolourchi (criterios 5 y 6): `TestFrame2DEulerCorotBatheBolourchi` en el mismo archivo.
  - Guarda (criterio 7): [tests/test_fiber_section.py](../../tests/test_fiber_section.py) · `TestGuardas`.
  - Camino con sección (criterio 8): [tests/validation/test_fiber_section_validation.py](../../tests/validation/test_fiber_section_validation.py) · `test_caso1_reproduccion_elastica_marcos_2d`, `test_caso8_corotacional_elemento_giro_grande`, `test_caso8_corotacional_arco_de_Bathe_Bolourchi`; tangente frente a diferencias finitas con plasticidad y giro de sólido rígido grande en `tests/test_fiber_section.py` · `TestTangente.test_frame2d_euler_corot_con_giro_grande`.
  - Masa (criterio 9): [tests/test_frame_mass.py](../../tests/test_frame_mass.py) · `TestMarco2DEulerBernoulli` (parametrizado con los dos marcos, incluye `test_corrotacional_igual_al_lineal`), `TestCaminoConSeccion`, `test_traslacion_rigida_2d` y `test_giro_rigido_2d_consistente`.
- **Notas de traducción**:
  - La formulación se basa en 3 DOFs deformacionales $[u, \bar\theta_1, \bar\theta_2]$ y una matriz $\mathbf B$ $3\times 6$ que los relaciona con los DOFs globales. Evaluada en configuración corriente.
  - **Signo de la contribución de momentos en $\mathbf K_\sigma$**: tras derivar $\partial(\mathbf z/l)/\partial\mathbf u_e$ analíticamente se obtiene $-(1/l^2)(\mathbf r\mathbf z^\top + \mathbf z\mathbf r^\top)$. El signo es **negativo** — un error común es ponerlo positivo. El test de diferencias finitas lo cazó en la primera implementación.
  - `compute_internal_forces` proyecta $\mathbf F_{\text{int}}$ al sistema local corriente para reportar axial/cortante separados, utilizando la submatriz $2\times 2$ de rotación del eje corriente.
  - `_unwrap` lleva los ángulos al intervalo $(-\pi, \pi]$; la formulación es válida mientras la rotación acumulada por paso sea menor que $\pi$.
  - Camino con sección: $\mathbf B_g$ (2×3) del sistema básico se precalcula al construir (`_build_basic_B`); `_basic_response` devuelve $[N, M_1, M_2]$ y $\mathbf K_{ll}$ integrados, y el resto de `compute_element_state` (matriz $\mathbf B$ 3×6, $\mathbf K_\sigma$) no distingue el camino. `internal_forces` y `compute_internal_forces` evalúan sobre una instantánea del trial y la restauran; el diccionario legado añade `section_deformations` y `section_forces` `(n_ip, 2)`. Masa y peso propio en la configuración de referencia con $\rho A$ de la sección.
  - Masa (§10): las mismas funciones de `_shared.py` que `Frame2DEuler` (`_frame2d_consistent_mass_local`, `_frame2d_lumped_mass_local`) con la densidad lineal de `_linear_density` y la $\mathbf T$ de referencia reconstruida desde `alpha0`.

---

## Diálogo

- **2026-04-21** · Implementación siguiendo Crisfield §7.3. Convive en `frame.py` con las vigas lineales por familia temática (todas son vigas 2D), sin herencia ni helpers compartidos.
- **2026-05-13** · `frame.py` se parte en paquete `solidum/elements/frame/`. EulerCorot vive ahora en `euler_corot.py`, sin herencia con las otras vigas 2D. Pasan a compartirse los helpers libres de carga de cuerpo, masa y traducción de fuerzas a través de `_shared.py`; la lógica corotacional propia (reconstrucción de `T` desde `alpha0`, cinemática actual) se mantiene íntegra dentro de la clase.
- **2026-04-21** · La primera versión tenía el signo equivocado en la parte de momentos de $\mathbf K_\sigma$ ($+$ en lugar de $-$). El test de coherencia tangente/fuerza interna por diferencias finitas lo detectó inmediatamente: $\|\mathbf K_{T,\text{analítica}} - \mathbf K_{T,\text{FD}}\|_{\text{rel}}$ del orden de 1e-1 en lugar de 1e-5. Tras derivar $\partial(\mathbf z/l)/\partial\mathbf u_e$ cuidadosamente (ver notas) el signo correcto es negativo. Este episodio justifica retrospectivamente incluir un test FD como criterio de aceptación en toda formulación no-lineal geométrica futura.
- **2026-04-21** · El test del límite lineal usa $P=10$ N para obtener $v/L \sim 10^{-4}$: suficientemente pequeño para estar en régimen lineal pero suficientemente grande para que el residuo del solver converja por debajo de `tol=1e-8` (cargas menores producen desplazamientos del orden del ruido numérico y el Newton oscila en la tolerancia relativa).
- **2026-09-22** · Auditoría global: `internal_forces(U, equivalent_load=None)`: cuando el `Assembler` ensambló peso propio o fuerza de cuerpo, `solidum.run` pasa el vector nodal equivalente del elemento (escalado por λ_final) y las fuerzas internas de extremo son `F_int − f_eq` en ejes locales. Antes se devolvía `K·u`, que con carga distribuida son fuerzas nodales, no fuerzas internas: un voladizo bajo peso propio daba V_i = qL/2 y M_i = −5qL²/12 en vez de qL y −qL²/2, contradiciendo a las reacciones.
- **2026-09-29** · Sección por fibras, etapa 1 ([ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md), [FiberSection](FiberSection.md)). El modelo no lineal anterior —$E_t$ de la deformación axial escalando toda $\mathbf K_{ll}$— desaparece: con fluencia axial anulaba el flector y con ablandamiento le cambiaba el signo. Dos caminos excluyentes: `material` + `A`, `I` (forma cerrada, sin cambio de resultados, sólo materiales con `IS_LINEAR_ELASTIC = True`) y `section` + `n_ip` (la sección se integra en el sistema básico corrotado; la cinemática corrotacional y $\mathbf K_\sigma$ no cambian). Contrato: parámetros `section` y `n_ip`; `A` e `I` obligatorios sólo con material; criterios 7 (guarda) y 8 (camino con sección, verificado a 1e-10 con giros de hasta 150° y en el arco de Bathe-Bolourchi hasta θ = π/2).
- **2026-09-29** · Masa sin inercia rotatoria (hallazgo clas-4 de la auditoría; **cambio físico aprobado por el usuario**), igual que en [Frame2DEuler](Frame2DEuler.md): se retira el bloque $\tfrac{\rho I L_0}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$ de los giros (cinemática de Timoshenko; Euler-Bernoulli desprecia la inercia rotatoria) y la concentrada pasa de $\rho I L_0/2$ a la HRZ $\rho A L_0^3/78$ por giro. Por los dos caminos; la masa ya no depende de $I$. Contrato: §10 y criterio 9.
- **2026-10-27** · `compute_body_load` acepta un campo `b(x)` (ADR 0022, fase 4), repartido en la geometría de referencia con las funciones de Hermite, como el `Frame2DEuler`; declara `REFERENCE_CELL = "line"` y su campo interior es el linealizado en la referencia (no entra en la verificación por soluciones manufacturadas, que es lineal). Se corrige una nota desactualizada de `out_of_scope`: la fuerza de cuerpo ya la repartía `compute_body_load`.
- **2026-11-27** · **Validación de §11 como usuario** (`validador-solidum`, versión corta aprobada por el usuario, «ok la 1»; reporte `informes/reportes/reporte_cargas_linea_3d_corot_2026-11-27.pdf`, corridas en `informes/cargas_linea_yaml/corot/`). Veredicto: **funciona**. Marco portal con carga pequeña: la diferencia con `Frame2DEuler` es proporcional a la carga (1.5·10⁻³, 10⁻⁴ y 10⁻⁵ en las reacciones con factores 1, 0.1 y 0.01). Voladizo con gran flecha bajo carga vertical uniforme muerta ($qL^3/EI = 10$, flecha 0.70 $L$, giro 60°) frente a la elástica resuelta aparte por dos métodos que coinciden a 10⁻¹²: con 4, 8, 16 y 32 elementos el error de la punta baja con orden 1.94–2.02 y queda por debajo de 0.1 % desde 16; el momento del empotramiento es $q\sum w_n x_n$ sobre la deformada (error ≤ 3.3·10⁻¹⁰) y el diagrama cierra con él y en los nudos. `q_local` es muerta: con giros de 51 a 60° coincide con la carga global equivalente (diferencia 0 en la historia, 2·10⁻¹⁰ en los diagramas). Observaciones de software: la columna `s` del archivo de diagramas es $\xi\,l$ sobre la cuerda corriente (falta decirlo en el manual); la figura `deformed` de los marcos no colorea por `field` ni dibuja la configuración inicial.
- **2026-11-27** · **Diagrama en una perturbación lineal sobre este elemento con carga repartida** (consulta a `teorico-solidum`, `informes/teorico/consulta_perturbacion_corot_2026-11-27.md`). **Decisión del usuario (2026-11-27, «si»): opción (a)**, formulación en §11. Hoy las fuerzas de extremo de la perturbación van en los ejes de la cuerda del estado base y el diagrama usa los de $u_b + \Delta U$; además al flector le falta un término geométrico de primer orden, $\Delta M_\text{geo}(\xi) = (\Delta l/l_b)\,[M_b(\xi) - M_{i,b}] + l_b\,\Delta\alpha\,[\xi N_{i,b} - L_0\int_0^\xi(\xi-\xi')\,q_{x',b}\,d\xi']$ (alargamiento de los brazos y efecto P-Δ del tramo). Recomendación del teórico: diagrama en los ejes y la longitud de la cuerda del estado base más $\Delta M_\text{geo}$ (cierre a 8·10⁻¹¹; error de $M$ frente al no lineal de 2.3 % a 0.22 % con $\Delta\alpha \approx 0.03$; lineal en $\Delta U$ y sin cambiar `element_forces`). Aparte, y también del usuario: si $N$ y $V$ de una perturbación siguen en los ejes del estado base o pasan a la linealización con la cuerda que gira. Las armaduras corrotacionales reciben la misma corrección de ejes (su $N$ no lleva término geométrico; en `Truss2DCorot` con carga transversal el diagrama de una perturbación no cerraba: −0.111 con $\Delta\alpha = -0.0109$, *medido* por el teórico).
- **2026-11-28** · **Implementación de la opción (a) verificada** por `teorico-solidum` (`informes/teorico/verificacion_teorica_2026-11-28_perturbacion.md`, «sí al commit»): coincide con su deducción independiente a 10⁻¹⁰ en cinco combinaciones de carga; cierre a redondeo con carga uniforme en el estado base; orden 2 del flector frente al no lineal (también con desplazamiento lateral libre, sección por fibras, base en dos etapas y perturbación sólo por desplazamiento impuesto); superponible; los elementos lineales, bit a bit sin cambio (*medido*). Atendidos sus hallazgos menores: caso con barra axialmente blanda en la prueba de orden (protege la forma interior del alargamiento de los brazos: el mutante $\xi^2 V_i$ cerraba en los extremos y pasaba), casos sin carga en el estado base y `Truss3DCorot` en la prueba de cierre, error si el estado base perdió las muestras de su carga, y el matiz del cierre con carga variable (de segundo orden, 1.2·10⁻⁴ con carga trapecial) en el catálogo y el ADR 0044. Mutaciones: 6/6 (agente principal) y 9/10 + la corregida (teórico).
