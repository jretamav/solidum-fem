# SECCIÓN POR FIBRAS PARA MARCOS EULER-BERNOULLI

> Orden de trabajo. La especificación física y la formulación numérica las decidió el usuario el 2026-09-29, una cuestión por vez (ver *Diálogo*); la arquitectura que las aloja está en el [ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md).

---

## Especificación física

### 0. Descripción general
Ley de sección para marcos Euler-Bernoulli en la que las fuerzas internas de la sección ($N$, $M_z$ en 2D; $N$, $M_y$, $M_z$ en 3D) se obtienen **integrando el esfuerzo axial $\sigma_{xx}$ de cada fibra** sobre la sección transversal. Cada fibra tiene su propio material 1D del catálogo y su propia historia, de modo que la plastificación o el daño avanzan fibra a fibra desde las fibras extremas hacia el eje neutro: **plasticidad distribuida en la sección**. Una misma sección puede mezclar materiales (concreto y acero). La consumen `Frame2DEuler`, `Frame2DEulerCorot` y `Frame3D` (etapa 1); `Frame2DTimoshenko` queda para la etapa 2.

Sustituye la ley anterior de los marcos con material inelástico, que escalaba toda la rigidez local con la $E_t$ de la deformación axial: con fluencia axial el flector se anulaba y con ablandamiento cambiaba de signo (auditoría 2026-09-29).

### 1. Hipótesis
- **Secciones planas y perpendiculares al eje deformado** (Euler-Bernoulli): $\gamma_{xy} = \gamma_{xz} = 0$ en la cinemática.
- **Fibras uniaxiales**: cada fibra sólo tiene $\varepsilon_{xx}$ y $\sigma_{xx}$; no hay esfuerzos cortantes en las fibras.
- **Torsión elástica y desacoplada** (3D): $T = GJ\,\theta_x'$, con $GJ$ constante declarado en la sección. No se degrada cuando las fibras fluyen ni interactúa con $N$, $M_y$, $M_z$.
- **Pequeñas deformaciones** en la sección. En `Frame2DEulerCorot` los giros de sólido rígido pueden ser grandes; la sección trabaja en el sistema local corrotado.
- **Sin regularización del ablandamiento** (ver §8).

### 2. Cinemática de la sección
Coordenadas $(y, z)$ de la fibra $i$ en los ejes locales del elemento, medidas desde el **eje del elemento** (la línea de nodos), no desde el centroide: una sección excéntrica acopla $N$ y $M$ de forma natural. Con la convención de Reglas §5 ($M_z > 0$ ⇒ tracción en $y < 0$; $M_y > 0$ ⇒ tracción en $z > 0$):

- 2D: $\ \varepsilon_i = \varepsilon_0 - y_i\,\kappa_z$, con $\mathbf e = [\varepsilon_0,\ \kappa_z]^\top$.
- 3D: $\ \varepsilon_i = \varepsilon_0 - y_i\,\kappa_z + z_i\,\kappa_y$, con $\mathbf e = [\varepsilon_0,\ \kappa_y,\ \kappa_z]^\top$.

### 3. Constitutiva por fibra
$(\sigma_i,\ E_{t,i},\ \text{estado}_i') = \text{material}_i.\text{compute\_state}(\varepsilon_i,\ \text{estado}_i)$, con cualquier material de `STRAIN_DIM = 1`.

### 4. Fuerzas internas y tangente de la sección
2D:
$$\mathbf s = \begin{bmatrix}N\\ M_z\end{bmatrix} = \sum_i \sigma_i A_i \begin{bmatrix}1\\ -y_i\end{bmatrix},\qquad
\mathbf k_s = \sum_i E_{t,i} A_i \begin{bmatrix}1 & -y_i\\ -y_i & y_i^2\end{bmatrix}.$$

3D:
$$\mathbf s = \begin{bmatrix}N\\ M_y\\ M_z\end{bmatrix} = \sum_i \sigma_i A_i \begin{bmatrix}1\\ z_i\\ -y_i\end{bmatrix},\qquad
\mathbf k_s = \sum_i E_{t,i} A_i \begin{bmatrix}1 & z_i & -y_i\\ z_i & z_i^2 & -y_i z_i\\ -y_i & -y_i z_i & y_i^2\end{bmatrix}.$$

Son las integrales $N = \int\sigma\,dA$, $M_y = \int z\,\sigma\,dA$, $M_z = -\int y\,\sigma\,dA$ de Reglas §5 aproximadas por la regla del punto medio sobre las fibras. $\mathbf k_s = \partial\mathbf s/\partial\mathbf e$ es exacta respecto de esa aproximación y simétrica para cualquier valor de $E_{t,i}$. En régimen elástico con sección simétrica respecto del eje, $\mathbf k_s = \mathrm{diag}(EA,\ EI)$ (2D) o $\mathrm{diag}(EA,\ EI_y,\ EI_z)$ (3D).

### 5. Elemento en desplazamientos
Interpolación de los marcos actuales: axial lineal y transversal cúbica de Hermite. En el sistema local, con $x \in [0, L]$:

- $\varepsilon_0 = (u_2 - u_1)/L$, constante en el elemento.
- $\kappa_z(x) = N_1''\,v_1 + N_2''\,\theta_{z1} + N_3''\,v_2 + N_4''\,\theta_{z2}$, con $N_1'' = (12x/L - 6)/L^2$, $N_2'' = (6x/L - 4)/L$, $N_3'' = (6 - 12x/L)/L^2$, $N_4'' = (6x/L - 2)/L$.
- 3D, plano $xz$ ($w' = -\theta_y$ por la regla de la mano derecha): $\kappa_y(x) = -w'' = -N_1''\,w_1 + N_2''\,\theta_{y1} - N_3''\,w_2 + N_4''\,\theta_{y2}$.
- 3D, torsión: $\theta_x' = (\theta_{x2} - \theta_{x1})/L$, $T = GJ\,\theta_x'$.

Con $\mathbf e(x_g) = \mathbf B_g\,\mathbf u_{\text{local}}$ en los puntos de Gauss-Legendre $x_g$ de pesos $w_g$:
$$\mathbf F_{\text{int,local}} = \sum_{g=1}^{n_{ip}} w_g\,\tfrac{L}{2}\,\mathbf B_g^\top\,\mathbf s_g\ (+\ \text{torsión}),\qquad
\mathbf K_{\text{local}} = \sum_{g=1}^{n_{ip}} w_g\,\tfrac{L}{2}\,\mathbf B_g^\top\,\mathbf k_{s,g}\,\mathbf B_g\ (+\ GJ/L\ \text{en los DOF de torsión}).$$

- **Número de puntos**: $n_{ip} = 3$ por defecto, mínimo 2. Con 2 puntos la rigidez elástica se integra exactamente (el integrando es cuadrático) y se reproduce el marco elástico en forma cerrada; con 1 punto el modo de flexión antisimétrica no tiene energía.
- **Corrotacional** (`Frame2DEulerCorot`): la sección trabaja en el sistema básico $\mathbf v = [\Delta l,\ \bar\theta_1,\ \bar\theta_2]$ del elemento corrotado. $\varepsilon_0 = \Delta l/L_0$, $\kappa_z(x) = [(6x/L_0 - 4)\,\bar\theta_1 + (6x/L_0 - 2)\,\bar\theta_2]/L_0$; las fuerzas básicas $\mathbf q = [N,\ M_1,\ M_2]$ y su tangente se integran igual y se transforman con la cinemática corrotacional existente (incluida su rigidez geométrica con $N$, $M_1$, $M_2$).

### 6. Masa y peso propio
$\rho A = \sum_i \rho_i A_i$, $\rho I_z = \sum_i \rho_i y_i^2 A_i$, $\rho I_y = \sum_i \rho_i z_i^2 A_i$, con las mismas fórmulas de masa consistente y concentrada de los marcos elásticos (en 3D, la inercia torsional con $\rho J_p = \rho I_y + \rho I_z$). El peso propio usa $\rho A$. Se desprecia el acoplamiento traslación-rotación de una distribución de masa excéntrica ($\sum\rho_i y_i A_i \ne 0$).

### 7. Límites de la formulación
- **Curvatura lineal y $\varepsilon_0$ constante** dentro del elemento: una rótula plástica se reparte sobre todo el elemento y hacen falta varios elementos por miembro cerca de las rótulas. Cuando el eje neutro se desplaza, $N$ en los puntos de Gauss sólo cumple el equilibrio axial en promedio.
- **Torsión no degradada** y sin interacción con la flexión: del lado no conservador cuando la torsión es relevante.
- **Sin deformación por cortante** (etapa 2: `Frame2DTimoshenko`).

### 8. Ablandamiento (sin regularización)
Con fibras en ablandamiento, el daño se localiza en uno o dos puntos de Gauss de un elemento y la energía disipada es proporcional a $l_c = w_g L/2$: la respuesta posterior al pico depende de la malla y de $n_{ip}$. Como el daño continuo de los sólidos, que tampoco está regularizado. La regularización genérica (longitud característica que cada elemento entrega a todo material con ablandamiento) es deuda registrada en STATUS ([deuda #29](../STATUS.md)).

---

## Contrato

```yaml
name: FiberSection
kind: section
status: validated

interface:
  deformations_2d: [eps0, kappa_z]
  forces_2d: [N, Mz]
  deformations_3d: [eps0, kappa_y, kappa_z]
  forces_3d: [N, My, Mz]
  torsion: "elástica desacoplada, T = GJ·θx' (sólo 3D)"
  consumers: [Frame2DEuler, Frame2DEulerCorot, Frame3D]

admission:
  textbooks:
    - bathe: "vigas isoparamétricas con integración numérica de σ sobre la sección para no linealidad material"
    - belytschko_liu_moran: "vigas y láminas con cuadratura a través del espesor para plasticidad"
  general_codes:
    - opensees: "dispBeamColumn con section Fiber (patch, layer, fiber); -GJ para la torsión 3D"
    - abaqus: "vigas B2x/B3x con integración numérica de la sección"

parameters:
  - { name: patches, type: list, required: false, desc: "Parches rectangulares (y, z, ny, nz, material), circulares (centro, radios, n_rad, n_circ, material), capas de barras (extremos, n, área, material)" }
  - { name: fibers, type: list, required: false, desc: "Fibras sueltas (y, z, A, material)" }
  - { name: GJ, type: float, required: false, desc: "Rigidez torsional elástica; obligatoria si la sección la usa un Frame3D" }

element_parameters:
  - { name: section, type: "id de sections", desc: "Excluyente con material + propiedades geométricas" }
  - { name: n_ip, type: int, default: 3, min: 2, desc: "Puntos de Gauss-Legendre a lo largo del elemento" }

acceptance:
  - name: reproducción elástica
    setup: "fibras elásticas; los tres marcos con n_ip = 2 frente a su versión en forma cerrada con EA y EI (EIy, EIz) integrados por las fibras"
    expect: "K y F_int iguales; EI de n capas = b h³/12 · (1 − 1/n²)"
    tol_rel: 1.0e-10

  - name: momento–curvatura, sección rectangular elastoplástica perfecta
    setup: "Elastoplastic1D con H = 0, N = 0, curvatura uniforme impuesta; 2D y los dos ejes de Frame3D"
    expect: "M_y = σ_y b h²/6, M_p = σ_y b h²/4, κ_y = 2σ_y/(E h), M(κ) = M_p [1 − (κ_y/κ)²/3] para κ ≥ κ_y"
    tol_rel: "del orden del error de discretización de las fibras, O(1/n²)"

  - name: interacción axial–flector
    setup: "N impuesta, curvatura creciente"
    expect: "M → M_p [1 − (N/N_p)²], N_p = σ_y b h"

  - name: colapso de un voladizo con carga en el extremo
    setup: "mallas de 4, 8, 16 y 32 elementos, control por desplazamiento"
    expect: "carga máxima decreciente con la malla y convergente a P_c = M_p/L"

  - name: colapso de una viga biempotrada con carga puntual central
    setup: "mallas crecientes, control por desplazamiento"
    expect: "carga máxima decreciente con la malla y convergente a P_c = 8 M_p/L (rótulas en empotramientos y centro)"

  - name: sección circular en flexión biaxial a 45°
    setup: "Frame3D, sección circular maciza de diámetro d, curvatura impuesta a 45° de los ejes"
    expect: "misma curva M–κ que en flexión uniaxial; M_y = π σ_y d³/32, M_p = σ_y d³/6"

  - name: convenciones de signo
    setup: "curvaturas positivas en cada eje"
    expect: "Mz > 0 tracciona las fibras con y < 0; My > 0 las fibras con z > 0 (Reglas §5)"

  - name: corrotacional elástico con giros grandes
    setup: "Frame2DEulerCorot con fibras elásticas frente a la versión con A e I integrados"
    expect: "misma respuesta"
    tol_rel: 1.0e-10

references:
  - "Spacone E., Filippou F.C., Taucer F.F. (1996). Fibre beam–column model for non-linear analysis of R/C frames: Part I. Formulation. Earthquake Engineering and Structural Dynamics 25:711–725."
  - "Bathe K.-J. (2014). Finite Element Procedures, 2.ª ed., cap. 6 (elementos de viga con no linealidad material)."
  - "Neal B.G. (1977). The Plastic Methods of Structural Analysis, 3.ª ed. (M–κ de la sección rectangular, interacción N–M, cargas de colapso)."
```

---

## Implementación

- **Archivo**: [solidum/sections/fiber.py](../../solidum/sections/fiber.py) — la sección y sus primitivas de discretización; el `__init__.py` del paquete [solidum/sections/](../../solidum/sections/) las reexporta (`from solidum.sections import FiberSection, RectangularPatch, CircularPatch, BarLayer, Fiber`). El paquete `solidum/sections/` lo recorre el descubrimiento automático como `materials/` y `elements/`.
- **Clase**: `FiberSection`, registrada en `SectionRegistry` (familia YAML `sections`, `SPEC_KIND = "section"`, [solidum/registry.py](../../solidum/registry.py)) y subclase de la base abstracta `Section` ([solidum/core/section.py](../../solidum/core/section.py)), que fija el contrato que consumen los marcos: `response_2d`, `response_3d`, `strains_2d`, `strains_3d`, `initial_states`, `materials`, `mass_properties`, `total_area` y el atributo `GJ`. Primitivas: `RectangularPatch`, `CircularPatch`, `BarLayer`, `Fiber`.
- **Consumidores**: la constitutiva común de los marcos vive en [solidum/elements/frame/_section.py](../../solidum/elements/frame/_section.py) (`FrameSectionMixin`: validación de los dos caminos, guarda de material elástico lineal, estado por punto de Gauss, integración `K = Σ w_g Bᵀ k_s B`, `F = Σ w_g Bᵀ s` y evaluación sin huella en el trial). Cada marco aporta su cinemática: [euler.py](../../solidum/elements/frame/euler.py) (`B_g` 2×6), [euler_corot.py](../../solidum/elements/frame/euler_corot.py) (`B_g` 2×3 del sistema básico) y [frame3d.py](../../solidum/elements/frame3d.py) (`B_g` 3×12 más la torsión).
- **Piezas genéricas que añade** ([ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md) §2 y §4): el gancho `from_yaml(data, resolve)` del lector YAML ([solidum/utils/yaml_parser.py](../../solidum/utils/yaml_parser.py), `_parse_families` y `_resolver`), `Material.IS_LINEAR_ELASTIC` y `Element.materials()`.
- **Tests**:
  - [tests/validation/test_fiber_section_validation.py](../../tests/validation/test_fiber_section_validation.py) — los ocho casos de `acceptance`, escritos a ciegas desde esta spec (43 tests):
    1. reproducción elástica — `test_caso1_capas_en_el_centro_y_EI_de_n_capas`, `test_caso1_rigidez_de_seccion_elastica`, `test_caso1_reproduccion_elastica_marcos_2d` y `test_caso1_reproduccion_elastica_frame3d` (sección rectangular y mixta concreto + acero; `n_ip = 2` y 3; elementos oblicuos);
    2. momento–curvatura — `test_caso2_momento_curvatura_seccion` y `test_caso2_momento_curvatura_elemento` (2D y los dos ejes de `Frame3D`, con descarga elástica que prueba la historia de las fibras);
    3. interacción axial–flector — `test_caso3_interaccion_axial_flector` (cinco combinaciones de eje y `N/N_p`, en tracción y en compresión);
    4. colapso del voladizo — `test_caso4_colapso_voladizo_carga_en_el_extremo`;
    5. colapso de la viga biempotrada — `test_caso5_colapso_viga_biempotrada_carga_central`;
    6. sección circular a 45° — `test_caso6_disposicion_de_la_seccion_circular`, `test_caso6_momento_curvatura_circular_uniaxial_y_45` y `test_caso6_frame3d_flexion_biaxial_45`;
    7. convenciones de signo — `test_caso7_signos_deformaciones_de_fibra` y `test_caso7_fibra_excentrica_fuerzas_internas_y_tangente`;
    8. corrotacional con giros grandes — `test_caso8_corotacional_elemento_giro_grande` (α = 60°, 150°, −120°) y `test_caso8_corotacional_arco_de_Bathe_Bolourchi` (voladizo enrollado hasta θ = π/2).
  - [tests/test_fiber_section.py](../../tests/test_fiber_section.py) — pruebas unitarias de la implementación (116 tests): construcción y validaciones de la sección y de las primitivas, cinemática y signos, tangente frente a diferencias finitas en régimen inelástico (sección mixta con daño en ablandamiento y acero con endurecimiento; los tres marcos, el corrotacional con giro grande y `Frame3D` con torsión), guardas de los marcos, estado por punto de Gauss, masa y peso propio, simetría agregada (ADR 0003), YAML de extremo a extremo y exportación VTK.
- **Notas de traducción**:
  - **Forma matricial de §4.** Con `W` la matriz de filas $\mathbf w_i^\top = [1,\ -y_i]$ (2D) o $[1,\ z_i,\ -y_i]$ (3D): $\boldsymbol\varepsilon = \mathbf W\mathbf e$, $\mathbf s = \mathbf W^\top(\boldsymbol\sigma\circ\mathbf A)$ y $\mathbf k_s = \mathbf W^\top\operatorname{diag}(\mathbf E_t\circ\mathbf A)\,\mathbf W$. `W` se construye una vez al crear la sección. $\mathbf k_s$ se simetriza con $\tfrac12(\mathbf k_s + \mathbf k_s^\top)$ para eliminar el redondeo del producto; la matriz local del elemento, igual.
  - **Sin traducción de signos.** La cinemática de §2 es la de Reglas.md §5 tal cual; la API pública 2D (`internal_forces`) sigue exponiendo `V` con el signo de viga estructural a través de la misma traducción de siempre (`_frame2d_forces_from_local`), que no depende del camino.
  - **Sección sin estado.** `initial_states()` devuelve una lista nueva de `None` (estado virgen de cada material) y `response_*` devuelve la lista de estados trial sin tocar la recibida. En el elemento, `state.vars[g]` guarda la lista de estados de las fibras del punto `g` y `state.stresses[g]` las fuerzas internas $\mathbf s_g$; `N_INTEGRATION_POINTS = n_ip` es atributo de instancia (la clase conserva `1`, el valor del camino con material).
  - **Evaluación agrupada por material** ([ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md) §6): las fibras se agrupan por identidad del material y se evalúan con un bucle escalar sobre `compute_state`. Los marcos no van por el camino por lotes del ADR 0014; se vectorizará si la medición lo justifica.
  - **Cuadratura.** `gauss_legendre_unit(n)` da los puntos $\xi_g \in (0, 1)$ y pesos que suman 1, de modo que $\int_0^L f\,dx \approx L\sum_g w_g f(\xi_g L)$: es la regla de §5 con $w_g^{(0,1)} = w_g^{(-1,1)}/2$.
  - **Corrotacional.** $\mathbf B_g$ (2×3) del sistema básico: $\varepsilon_0 = \Delta l/L_0$, $\kappa_z = [(6\xi - 4)\,\bar\theta_1 + (6\xi - 2)\,\bar\theta_2]/L_0$. Las fuerzas básicas $\mathbf q = [N, M_1, M_2]$ y su tangente entran en la cinemática corrotacional existente sin cambios, incluida la rigidez geométrica.
  - **Frame3D.** La fila de $\kappa_y$ de $\mathbf B_g$ es $[-N_1'',\ N_2'',\ -N_3'',\ N_4'']$ sobre $(w_1, \theta_{y1}, w_2, \theta_{y2})$; la torsión suma $GJ/L$ en los DOF $\theta_x$ y $T = GJ\,(\theta_{x2} - \theta_{x1})/L$. Una sección sin `GJ` se rechaza al construir el elemento.
  - **CircularPatch.** División radial en `n_rad` anillos de **igual espesor** $\Delta r = (r_\text{ext} - r_\text{int})/n_\text{rad}$ (como el `patch circ` de OpenSees) y angular en `n_circ` sectores iguales. Área exacta del sector $A = \tfrac{\Delta\theta}{2}(r_2^2 - r_1^2)$ y centroide a $r_c = \tfrac23\,\tfrac{\sin\alpha}{\alpha}\,\tfrac{r_2^3 - r_1^3}{r_2^2 - r_1^2}$ del centro sobre la bisectriz ($\alpha = \Delta\theta/2$); ángulos desde $+y$ hacia $+z$, primer sector en 0. Con `n_circ` múltiplo de 8 la sección discreta es invariante ante un giro de 45°, que es lo que hace exacto el caso 6.
  - **Masa y peso propio** (§6). `mass_properties()` devuelve `rhoA`, `rhoIy`, `rhoIz` respecto del eje del elemento y lanza `ValueError` con los materiales sin `density`; los marcos llaman a las mismas funciones de masa consistente y concentrada del camino con material con $(\rho, A, I) = (1, \rho A, \rho I_z)$, y `Frame3D` con $\rho J_p = \rho I_y + \rho I_z$. La fuerza de cuerpo uniforme usa el área bruta $\sum A_i$ (`total_area()`); el peso propio, $\rho A$.
  - **Post-proceso.** `internal_forces` y `compute_internal_forces` evalúan el elemento sobre una instantánea del trial y la restauran: no son un iterado del paso. El diccionario legado añade `section_deformations` y `section_forces`, de forma `(n_ip, 2)` o `(n_ip, 3)`. En VTK, `Internal_State` promedia la `PRIMARY_STATE_VAR` sobre fibras y puntos de Gauss.
  - **YAML.** Cada parche declara `shape` ∈ {`rectangle`, `circle`, `layer`} y sus parámetros son los del constructor de la primitiva (se validan por introspección: una clave desconocida es error). Los errores nombran el contexto, p. ej. `sections[0] (id=7), patches[0].material: referencia a material inexistente (9)`.

---

## Resultados de validación

Datos comunes: $E = 200$ GPa, $\sigma_y = 250$ MPa, $b = 0.3$ m, $h = 0.5$ m, $L = 3$ m. Valores medidos el 2026-09-29.

| Caso | Resultado | Referencia |
|---|---|---|
| 1. Reproducción elástica | $K$ y $F_{int}$ de los tres marcos iguales a la forma cerrada a 1e-10 con `n_ip` = 2 y 3; $EI$ de 10 capas igual a $Ebh^3/12\,(1 - 1/n^2)$ a redondeo | forma cerrada |
| 2. M–κ rectangular, 200 capas | $\max\lvert\Delta M\rvert/M_p$ = 1.7e-5 en los tres ejes (cota $\tfrac53/n^2$ = 4.2e-5); $M$ con todas las fibras plastificadas = $M_p$ a 2e-15 | Neal (1977) |
| 3. Interacción N–M | $\max\lvert\Delta M\rvert/M_p$ = 1.7e-5 (cota $\tfrac83/n^2$ = 6.7e-5); en $\kappa = 100\,\kappa_y$ y $N/N_p$ = ±0.3, $M/M_p$ = 0.90995 frente a la forma cerrada 0.90997 y al límite plástico $1 - (N/N_p)^2$ = 0.91 (con 0.6: 0.63995, 0.63997 y 0.64) | Neal (1977) |
| 4. Voladizo, $P_{max}/P_c$ | 1.0710, 1.0351, 1.0176, 1.0088 con 4, 8, 16, 32 elementos (cota del modelo discreto 1.0767, 1.0369, 1.0181, 1.0090) | $P_c = M_p/L$ |
| 5. Biempotrada, $P_{max}/P_c$ | 1.3475, 1.1496, 1.0710, 1.0351 con 4, 8, 16, 32 elementos (cota 1.3983, 1.1661, 1.0767, 1.0369) | $P_c = 8M_p/L$ |
| 6. Circular, 20 × 64 sectores | $\max\lvert\Delta M\rvert/M_p$ = 7.2e-4 (cota 6.0e-3); a 45°, $M_y = M_z$ y $\lvert\mathbf M\rvert$ igual a la uniaxial a 1e-9 $M_p$, en la sección y en `Frame3D` | $M_y = \pi\sigma_y d^3/32$, $M_p = \sigma_y d^3/6$ |
| 7. Signos | $M_z > 0$ tracciona $y < 0$; $M_y > 0$ tracciona $z > 0$; fibra excéntrica con $\mathbf s$ y $\mathbf k_s$ exactos | Reglas.md §5 |
| 8. Corrotacional elástico | $K$ y $F_{int}$ iguales a 1e-10 con giros de 60°, 150° y −120°; arco de Bathe-Bolourchi hasta θ = π/2 con $\lVert\Delta U\rVert/\lVert U\rVert \le$ 1e-10 | versión con $A$ e $I$ integrados |

Los casos 4 y 5 sobreestiman la carga de colapso en primer orden del tamaño del elemento (spec §7): el elemento en desplazamientos evalúa la sección en los puntos de Gauss, nunca en el extremo, y la rótula se reparte hacia el interior del elemento, con un brazo efectivo menor. La cota del modelo discreto la deriva el test por análisis de límite exacto del elemento raíz ($P_\infty/P_c = 1.398/(1.398 - 0.398/n)$ en el voladizo). En la biempotrada cada cuarto de vano es un voladizo de $n/4$ elementos, por eso sus cifras con $n$ elementos repiten las del voladizo con $n/4$.

---

## Diálogo

- **2026-09-29** · Decisiones del usuario, una por vez: (1) sección por fibras en lugar de restringir los marcos a material elástico; (2) elemento en desplazamientos, coherente con Reglas §0, no en fuerzas; (3) etapa 1 con los tres marcos Euler-Bernoulli, Timoshenko en una etapa 2 con su propia cinemática; (4) Gauss-Legendre, 3 puntos por defecto y mínimo 2; (5) torsión elástica desacoplada con GJ declarado en la sección; (6) sin regularización del ablandamiento en esta etapa, con la regularización genérica registrada como deuda; (7) los ocho casos de validación del contrato.
- **2026-09-29** · Implementada y validada: los ocho casos de `acceptance` pasan (43 tests en `tests/validation/test_fiber_section_validation.py`, escritos a ciegas desde esta spec, más 116 pruebas unitarias en `tests/test_fiber_section.py`); `status: validated`. Cifras en *Resultados de validación*. Dos observaciones para el usuario, ninguna es un defecto: (a) los casos 4 y 5 confirman §7 con números — la carga de colapso converge por arriba y en primer orden del tamaño del elemento, así que con pocos elementos por miembro la sobreestimación no es despreciable (biempotrada con 4 elementos en el vano: +35 %; con 32: +3.5 %); (b) la tolerancia del caso 2 (`tol_rel` "del orden de O(1/n²)") quedó fijada en el test como la cota $\tfrac53 M_p/n^2$ derivada para la regla del punto medio, y la del caso 3 como $\tfrac83 M_p/n^2$; los casos 4 y 5 se comprueban contra la cota de análisis de límite del modelo discreto y contra la predicción de la estática del elemento raíz con la M–κ analítica. Los cuatro marcos rechazan ahora un material inelástico en el camino con `material` + propiedades geométricas (`IS_LINEAR_ELASTIC`), con un mensaje que remite a esta sección.
