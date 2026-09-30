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

La deformación axial en un punto genérico de la sección es $\varepsilon(s, y) = \varepsilon_{axial}(s) - y\,\kappa(s)$; al integrar sobre la sección elástica y simétrica respecto del eje se obtienen axial ($N$) y momento flector ($M$) desacoplados.

### 4. Ecuaciones constitutivas
Dos caminos excluyentes, elegidos al construir el elemento ([ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md)):

- **Con `material` + `A`, `I`** (forma cerrada): $N = EA\,\varepsilon_{axial}$, $M = EI\,\kappa$. Sólo admite materiales **elásticos lineales** (`IS_LINEAR_ELASTIC = True`); un material inelástico se rechaza al construir con un mensaje que remite a la sección por fibras.
- **Con `section`** (+ `n_ip`): $N$ y $M$ son los de la ley de sección, que integra el esfuerzo $\sigma_{xx}$ de cada fibra con su propio material 1D ([FiberSection](FiberSection.md) §2-§4). Plasticidad y daño distribuidos en la sección.

Hasta el 2026-09-29 el camino con material admitía materiales inelásticos y escalaba toda la rigidez local con la $E_t$ de la deformación axial: con fluencia axial el flector se anulaba y con ablandamiento cambiaba de signo. Esa ley desapareció (ver *Diálogo*).

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

Es la forma cerrada del camino con material elástico lineal. En el camino con sección la matriz local se integra numéricamente, $\mathbf K_{\text{local}} = \sum_g w_g\,\tfrac{L}{2}\,\mathbf B_g^\top\,\mathbf k_{s,g}\,\mathbf B_g$ con $\mathbf e_g = [\varepsilon_0,\ \kappa_z(x_g)]^\top = \mathbf B_g\,\mathbf u_{\text{local}}$, las mismas funciones de forma del §8 y la tangente de la sección $\mathbf k_s$ ([FiberSection](FiberSection.md) §5); con fibras elásticas y `n_ip` ≥ 2 reproduce esta matriz exactamente.

### 10. Fuerzas internas
En el sistema local, camino con material:
$$\mathbf F_{\text{int,local}} = \mathbf K_{\text{local}}\,\mathbf u_{\text{local}}.$$

La componente axial se escribe con el esfuerzo que devuelve el material, $F_1 = -\sigma A$, $F_4 = +\sigma A$; con material elástico lineal coincide con la fila axial del producto. Camino con sección: $\mathbf F_{\text{int,local}} = \sum_g w_g\,\tfrac{L}{2}\,\mathbf B_g^\top\,\mathbf s_g$, con $\mathbf s_g = [N,\ M_z]$ de la sección en cada punto de Gauss.

### 11. Ensamblaje global
$$\mathbf K_{\text{global}} = \mathbf T^\top\,\mathbf K_{\text{local}}\,\mathbf T, \qquad \mathbf F_{\text{int}} = \mathbf T^\top\,\mathbf F_{\text{int,local}}.$$

### 12. Cuadratura
Camino con material: no aplica (forma cerrada; integración analítica de los polinomios de Hermite). El campo `n_integration_points: 1` se refiere al único estado material conceptual del elemento (la deformación axial media).

Camino con sección: Gauss-Legendre de `n_ip` puntos a lo largo del elemento (3 por defecto, mínimo 2). El elemento fija `N_INTEGRATION_POINTS = n_ip` como atributo de instancia; el atributo de clase sigue en 1.

### 13. Matriz de masa
**Sin inercia rotatoria de la sección.** La hipótesis de Euler-Bernoulli ($\theta = dv/ds$, §1) la desprecia: la energía cinética es la de la traslación del eje, $T = \tfrac12\int_0^L \rho A\,(\dot u^2 + \dot v^2)\,ds$. La masa depende sólo de la densidad lineal $\rho A$ ($\rho\cdot A$ con material; $\sum_i \rho_i A_i$ con sección, [FiberSection](FiberSection.md) §6), no de $I$.

**Consistente** (`lumping="consistent"`, por defecto), $\int_0^L \rho A\,\mathbf N^\top\mathbf N\,ds$ con las funciones de forma del §8:
$$\mathbf M_{\text{axial}} = \frac{\rho A L}{6}\begin{bmatrix}2 & 1\\ 1 & 2\end{bmatrix}\ \text{en}\ (u_1, u_2), \qquad
\mathbf M_{\text{flex}} = \frac{\rho A L}{420}\begin{bmatrix}
156 & 22L & 54 & -13L \\
22L & 4L^2 & 13L & -3L^2 \\
54 & 13L & 156 & -22L \\
-13L & -3L^2 & -22L & 4L^2
\end{bmatrix}\ \text{en}\ (v_1, \theta_1, v_2, \theta_2).$$

Reproduce exactamente la energía cinética de los movimientos rígidos, que la interpolación contiene: $\rho A L\,|\mathbf v|^2/2$ en una traslación y $\tfrac12\,\rho A\,\tfrac{L^3}{3}\,\omega^2$ en un giro alrededor de un extremo. Con ella las frecuencias discretas son cotas superiores de las de Euler-Bernoulli (Rayleigh-Ritz conforme).

**Concentrada** (`lumping="lumped"`), HRZ (Hinton, Rock y Zienkiewicz 1976): la diagonal consistente escalada para que los términos traslacionales de cada dirección sumen $\rho A L$. En la flexión la diagonal es $\tfrac{\rho A L}{420}(156, 4L^2, 156, 4L^2)$ y el factor $420/312$:
$$m_u = m_v = \frac{\rho A L}{2}, \qquad m_\theta = \frac{\rho A L^3}{78} \quad \text{por nodo}.$$
El bloque traslacional por nodo es $m\,\mathbf I_2$, invariante ante la rotación: la matriz es diagonal también en globales, requisito de `CentralDifferenceSolver`.

En globales, $\mathbf M_{\text{global}} = \mathbf T^\top\,\mathbf M_{\text{local}}\,\mathbf T$.

Hasta el 2026-09-29 la consistente sumaba $\tfrac{\rho I L}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$ en $(\theta_1, \theta_2)$ y la concentrada usaba $\rho I L/2$ por giro: inercia rotatoria de la cinemática de Timoshenko, no de Euler-Bernoulli (ver *Diálogo*).

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
  n_integration_points: 1   # camino con material; con section, N_INTEGRATION_POINTS = n_ip (atributo de instancia)

admission:                # criterio de admisión al programa principal (Reglas.md §1)
  textbooks:
    - zienkiewicz_taylor: "vol. 2 (6.ª ed.), cap. 10: barras de Euler-Bernoulli (formulación base, material elástico)"
    - bathe: "cap. 4-5 (elementos de viga)"
  general_codes:
    - opensees: "elasticBeamColumn (2D)"
    - abaqus: "B23"

parameters:
  - { name: A, type: float, required: false, desc: "Área de la sección transversal. Obligatoria con material; no se declara con section" }
  - { name: I, type: float, required: false, desc: "Momento de inercia respecto al eje perpendicular al plano (z). Obligatorio con material; no se declara con section" }
  - { name: section, type: "id de sections", required: false, desc: "Sección de marco (FiberSection). Excluyente con material, A e I (ADR 0021)" }
  - { name: n_ip, type: int, required: false, default: 3, desc: "Puntos de Gauss-Legendre a lo largo del elemento, ≥ 2. Sólo con section" }

material_contract:
  signature: "compute_state(ε, state) -> (σ, E_tangent, state')"
  strain_kind: "axial escalar (deformación media axial ε = (u_2 - u_1)/L en sistema local)"
  nonlinearity_model: "camino con material: sólo materiales elásticos lineales (IS_LINEAR_ELASTIC = True; un material inelástico se rechaza al construir). La no linealidad material entra por section (FiberSection): cada fibra evalúa su material 1D con ε = ε₀ − y·κ_z"

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
  - grandes rotaciones o desplazamientos (usar Frame2DEulerCorot)
  - material inelástico con material + A, I (se rechaza al construir; la plasticidad y el daño entran por section, FiberSection)
  - pandeo (requiere rigidez geométrica, no implementada en esta viga lineal)
  - "fuerzas de cuerpo no uniformes (trapezoidales, puntuales en el vano): sólo b uniforme vía compute_body_load"

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

  - name: guarda de material elástico lineal            # añadido 2026-09-29 (ADR 0021)
    setup: "material Elastoplastic1D o IsotropicDamage1D con A e I"
    expect: "ValueError al construir, con un mensaje que nombra IS_LINEAR_ELASTIC y remite a FiberSection con 'section'"

  - name: camino con sección                            # añadido 2026-09-29 (ADR 0021)
    setup: "section = FiberSection; casos de acceptance de docs/specs/FiberSection.md"
    expect: "reproducción elástica con n_ip = 2, momento–curvatura, interacción N–M, colapso del voladizo y de la viga biempotrada, signos"

  - name: matrices de masa consistente y concentrada   # añadido 2026-09-29
    setup: "elemento oblicuo, caminos con material y con sección (ρA = Σρ_i A_i); M_local = T·M·Tᵀ"
    expect: "consistente = axial ρAL/6·[[2,1],[1,2]] + hermítica ρAL/420·H, sin inercia rotatoria (independiente de I); concentrada HRZ = diag(ρAL/2, ρAL/2, ρAL³/78) por nodo, diagonal también en globales"
    tol_rel: 1.0e-12

  - name: energía cinética de movimientos rígidos        # añadido 2026-09-29
    setup: "traslación v y giro ω alrededor de un extremo"
    expect: "½u̇ᵀMu̇ = ρAL|v|²/2 (consistente y concentrada) y ½·ρA·L³/3·ω² (consistente)"
    tol_rel: 1.0e-12

references:
  - "Bathe K.-J., Finite Element Procedures, §4.2 (formulación isoparamétrica) y cap. 5 (vigas)"
  - "Cook R.D., Concepts and Applications of FEA, §2.7 (elementos de viga)"
  - "Hinton E., Rock T., Zienkiewicz O.C. (1976). A note on mass lumping and related processes in the finite element method. Earthquake Eng. Struct. Dyn. 4, 245-249 (masa concentrada HRZ)"
```

---

## Implementación

- **Archivo**: [solidum/elements/frame/euler.py](../../solidum/elements/frame/euler.py) — submódulo del paquete `solidum/elements/frame/` que aloja también `Frame2DTimoshenko` y `Frame2DEulerCorot`.
- **Clase**: `Frame2DEuler(FrameSectionMixin, Element)`. La construcción de la matriz de transformación $\mathbf T$ se delega a `build_geometry_2d` en [solidum/elements/frame/_shared.py](../../solidum/elements/frame/_shared.py), compartida con `Frame2DTimoshenko` (la versión corotacional reconstruye $\mathbf T$ desde `alpha0` por su lógica propia). La constitutiva de los dos caminos (validación de argumentos, guarda de material elástico lineal, estado por punto de Gauss, integración de la sección) vive en `FrameSectionMixin`, [solidum/elements/frame/_section.py](../../solidum/elements/frame/_section.py), compartida con `Frame2DEulerCorot` y `Frame3D` (ADR 0021). `REFERENCE_KWARGS = {"material": MaterialRegistry, "section": SectionRegistry}`.
- **Tests**: [tests/test_frame.py](../../tests/test_frame.py) · `TestFrame2DEulerAcceptance` — los tres criterios físicos y el registro:
  - `test_acceptance_respuesta_axial_pura` (criterio 1) — resuelve el voladizo con solver, verifica $u_x(L) = FL/(EA)$.
  - `test_acceptance_voladizo_carga_transversal` (criterio 2) — carga transversal $P$, verifica flecha $v = PL^3/(3EI)$ y rotación $\theta = PL^2/(2EI)$.
  - `test_acceptance_simetria_K` (criterio 3) — viga oblicua (c=0.6, s=0.8) para evitar ejes canónicos, verifica $\mathbf K = \mathbf K^\top$.
  - `test_registro_en_registry` — autodiscover.
  - Criterio 4 (guarda): [tests/test_fiber_section.py](../../tests/test_fiber_section.py) · `TestGuardas` (`test_material_inelastico_remite_a_la_seccion`, `test_exactamente_un_camino`, `test_n_ip_invalido`, `test_n_ip_solo_con_seccion`).
  - Criterio 5 (camino con sección): [tests/validation/test_fiber_section_validation.py](../../tests/validation/test_fiber_section_validation.py) (casos 1-5 y 7 de la spec FiberSection) y `TestMarcosConSeccion`, `TestTangente`, `TestMasa`, `TestPesoPropio`, `TestYaml` de `tests/test_fiber_section.py`.
  - Criterios 6 y 7 (masa): [tests/test_frame_mass.py](../../tests/test_frame_mass.py) · `TestMarco2DEulerBernoulli` (consistente, HRZ, independencia de $I$), `TestCaminoConSeccion`, `test_traslacion_rigida_2d`, `test_giro_rigido_2d_consistente` y `test_frecuencias_biapoyada_cota_superior_de_euler_bernoulli`; la hermítica se contrasta además con la integral $\int\rho A\,\mathbf N^\top\mathbf N$ (`test_hermitica_es_la_integral_de_las_funciones_de_hermite`). Diagonalidad y masa total de la concentrada en [tests/test_mass_lumping.py](../../tests/test_mass_lumping.py) · `TestLumpedFrame2D`.
- **Notas de traducción**:
  - La matriz $\mathbf T$ se calcula una sola vez en `__init__` sobre la configuración inicial y se guarda como atributo — coherente con el régimen geométricamente lineal.
  - Camino con material: la componente axial de $\mathbf F_{\text{int,local}}$ se escribe con el $\sigma A$ que devuelve el material; con material elástico lineal (el único admitido) coincide con la fila axial de $\mathbf K_{\text{local}}\mathbf u_{\text{local}}$.
  - Camino con sección: $\mathbf B_g$ (2×6) se precalcula al construir (`_build_section_B`) con $N_1''\ldots N_4''$ de `hermite_second_derivatives`; los pesos ya llevan la longitud ($w_g L$ con la regla sobre $(0, 1)$). `internal_forces` y `compute_internal_forces` evalúan sobre una instantánea del trial y la restauran; el diccionario legado añade `section_deformations` y `section_forces` `(n_ip, 2)`. Masa con $\rho A$ de la sección; peso propio con $\rho A$ y fuerza de cuerpo con el área bruta $\sum A_i$.
  - Masa (§13): `_frame2d_consistent_mass_local` y `_frame2d_lumped_mass_local` de `_shared.py`, que reciben la densidad lineal $\rho A$ (`_linear_density`). Los bloques de la flexión (`hermite_consistent_mass` y la masa rotacional HRZ `hrz_rotational_mass`) los comparten los tres marcos Euler-Bernoulli; el argumento opcional `rotary_inertia` $= \rho I$ sólo lo pasa `Frame2DTimoshenko`.

---

## Diálogo

- **2026-04-21** · Elemento movido a archivo propio `solidum/elements/frame.py` y desacoplado del helper compartido `_frame_geometry`. El método estático `_build_geometry` reimplementa la construcción de $\mathbf T$ dentro de la clase. `Frame2DTimoshenko` sigue en `structural.py` usando su propio acceso al helper compartido; si mañana también se independiza, duplicará o internalizará su propia versión.
- **2026-04-21** · Tests de aceptación cubren la flecha analítica del voladizo a 8 decimales ($v = PL^3/(3EI)$) usando el solver completo del proyecto — validan no solo la cinemática del elemento sino también su integración con ensamblador y solver.
- **2026-05-13** · `frame.py` se parte en paquete `solidum/elements/frame/{euler,timoshenko,euler_corot,_shared}.py`. La duplicación de `_build_geometry` entre Euler y Timoshenko se elimina: la construcción de $\mathbf T$ vuelve a vivir en un único lugar (`build_geometry_2d` en `_shared.py`), pero esta vez sin función helper externa flotante — pertenece al paquete `frame/`. Frame2DEulerCorot mantiene su propia reconstrucción de $\mathbf T$ desde `alpha0` porque su cinemática corotacional no se beneficia del helper común.
- **2026-09-22** · Auditoría global: `internal_forces(U, equivalent_load=None)`: cuando el `Assembler` ensambló peso propio o fuerza de cuerpo, `solidum.run` pasa el vector nodal equivalente del elemento (escalado por λ_final) y las fuerzas internas de extremo son `F_int − f_eq` en ejes locales. Antes se devolvía `K·u`, que con carga distribuida son fuerzas nodales, no fuerzas internas: un voladizo bajo peso propio daba V_i = qL/2 y M_i = −5qL²/12 en vez de qL y −qL²/2, contradiciendo a las reacciones.
- **2026-09-29** · Sección por fibras, etapa 1 ([ADR 0021](../adr/0021-secciones-de-marco-y-seccion-por-fibras.md), [FiberSection](FiberSection.md)). La auditoría de pureza verificó que el modelo no lineal anterior —evaluar el material con la deformación axial y escalar **toda** la matriz local con su $E_t$— daba flectores incorrectos en silencio: nulos con fluencia axial, de signo cambiado con ablandamiento. Desaparece. El elemento tiene ahora dos caminos excluyentes: `material` + `A`, `I`, la forma cerrada de siempre, sin cambio de resultados, que sólo admite materiales con `IS_LINEAR_ELASTIC = True` (un material inelástico se rechaza al construir y el mensaje remite a la sección por fibras); y `section` + `n_ip` (3 por defecto, mínimo 2), elemento en desplazamientos que integra la sección en Gauss-Legendre. Contrato: parámetros `section` y `n_ip`; `A` e `I` pasan a ser obligatorios sólo con material; criterios de aceptación 4 (guarda) y 5 (camino con sección). Un modelo que antes corría con material inelástico ahora falla al construir: antes daba resultados incorrectos.
- **2026-09-29** · Masa sin inercia rotatoria (hallazgo clas-4 de la auditoría; **cambio físico aprobado por el usuario**). La consistente sumaba $\tfrac{\rho I L}{6}\left[\begin{smallmatrix}2&1\\1&2\end{smallmatrix}\right]$ en $(\theta_1, \theta_2)$: la energía cinética $\tfrac12\rho I\dot\theta^2$ de un giro de sección interpolado linealmente e independiente de la flecha, que es cinemática de Timoshenko; la hipótesis de Euler-Bernoulli desprecia la inercia rotatoria. Se retira, por los dos caminos: la consistente es la clásica (axial lineal + hermítica) y la concentrada pasa de $\rho I L/2$ por giro (su justificación por OpenSees/SAP2000 no se sostenía) a la HRZ, $\rho AL/2$ traslacional y $m_\theta = \rho A L^3/78$ (§13). Consecuencias: la masa ya no depende de $I$; las frecuencias de marcos esbeltos cambian poco, porque la inercia rotatoria escala con $(r/L)^2$, y con la consistente vuelven a ser cotas superiores de las de Euler-Bernoulli. Contrato: §13 y criterios 6 y 7. `Frame2DTimoshenko` conserva su masa con inercia rotatoria.
