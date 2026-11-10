# GlobalCrackTracking — seguimiento global de la trayectoria de grieta

> Orden de trabajo. Redactada por la IA el 2026-11-10 a partir de las decisiones del usuario sobre su formulación (Retama 2010), de los informes del `formulador` y del `teorico_fem`, y de la búsqueda bibliográfica sobre seguimiento de grieta. Aprobada por el usuario el 2026-11-10; implementada y validada el mismo día, con los ajustes del Diálogo (ε, prioridad de la propagación, regla del vértice en la punta, tolerancias de energía).

---

## Especificación física

### 0. Descripción general

Al inicio de cada paso, el algoritmo decide qué elementos `CST_Embedded2D` se agrietan y con qué partición, es decir, qué nodo queda solo a un lado de la grieta. El objetivo es que las grietas elementales formen trayectorias continuas y coherentes entre vecinos. El enriquecimiento sigue siendo elemental (KOS): cada elemento conserva su propio salto. Lo que se impone es que dos vecinos asignen siempre el mismo lado a cada nodo que comparten. Sin esa coherencia, una grieta recta en una malla no alineada se bloquea: la fuerza queda en 1.2 a 4.3 veces el pico, o en 3–24 % del pico con la ley de modo I, en vez de caer a cero (informes del 2026-11-10).

Decisiones del usuario que fija esta spec (2026-11-10):
1. El nodo solitario lo dicta la trayectoria, no una regla elemental.
2. Seguimiento **global**: un problema escalar auxiliar, tipo conducción anisótropa (Oliver et al. 2004a, 2004b; en concreto reforzado 2D, Oliver et al. 2008). La trayectoria es una isolínea y el lado de cada nodo es el signo del campo.
3. Si el nodo aislado por la trayectoria da $\cos\psi \le 0$, la trayectoria pasa por ese vértice: el nodo cambia de lado.

La formulación del elemento no cambia: $l_d = (A/h)\cos\psi$, de la ec. (8) de la deducción del usuario (tesis, ec. 6.24), con $\cos\psi > 0$ en todo elemento cortado gracias a la decisión 3.

### 1. Campo de direcciones

En cada elemento $e$ de la región de seguimiento (los `CST_Embedded2D` del modelo) hay una tangente unitaria $\mathbf t_e$ de la trayectoria:
- en un elemento agrietado, la tangente $\mathbf s_e$ fijada al activarse, que no se reorienta;
- en un elemento intacto, la perpendicular a la dirección del esfuerzo principal mayor $\sigma_I$ del estado convergido, el mismo criterio de Rankine de la activación.

Es un campo de líneas: $\mathbf t_e$ y $-\mathbf t_e$ son equivalentes.

### 2. Problema auxiliar de cada grieta

Cada grieta $k$ tiene un campo escalar $\theta_k$ lineal por elemento (los nodos de la malla) que minimiza

$$\Pi(\theta_k) = \tfrac12\int_\Omega (\mathbf t\cdot\nabla\theta_k)^2\,d\Omega + \tfrac12\,\varepsilon\int_\Omega |\nabla\theta_k|^2\,d\Omega .$$

El primer término hace a $\theta_k$ constante a lo largo de las líneas del campo de direcciones, de modo que sus isolíneas son trayectorias. Equivale a una conducción con conductividad $\mathbf t\otimes\mathbf t$. El segundo término regulariza la matriz, que sin él es singular (Riccardi et al. 2017). $\varepsilon$ es adimensional: $\varepsilon = 10^{-6}$. El regularizador difunde $\theta_k$ a través de las líneas del campo una distancia del orden de $\sqrt{\varepsilon}\,L$, independiente de $h$; con la propuesta inicial ($10^{-4}$) la isolínea se apartaba de una trayectoria exacta 0.18·h en la malla más fina (Diálogo, 2026-11-10). La prueba de aceptación comprueba que la trayectoria apenas cambia por debajo de $10^{-6}$.

**Condiciones de Dirichlet: la grieta ya formada fija su trayectoria.** En los nodos de los elementos ya agrietados de $k$, $\theta_k$ vale la distancia con signo a la grieta: el signo es el lado ya asignado y la magnitud es la distancia a la recta de la grieta en ese elemento, con un mínimo de $10^{-3}h_e$. Al nacer la grieta, la raíz es el elemento inicial, con la recta por su centroide y normal de Rankine. Así la parte formada no se mueve cuando rotan las direcciones principales, y el campo solo sigue la dirección actual por delante de las puntas. Este es el remedio de Riccardi et al. (2017) para la pérdida de continuidad cuando rotan las direcciones y no se actualiza la raíz.

### 3. Trayectoria y partición

La trayectoria de $k$ es la isolínea $\theta_k = 0$, y el lado de cada nodo respecto de $k$ es el signo de $\theta_k$. Un elemento está cortado por $k$ si sus nodos tienen los dos signos. La partición es coherente por construcción, porque el lado de cada nodo es único. La posición de la recta dentro del elemento no entra en las ecuaciones: solo cuentan el nodo solitario y $l_d$ (informes del 2026-11-10). Se guarda para el posproceso.

### 4. Activación (al inicio del paso, con el estado convergido)

Los elementos intactos con $\sigma_I > \sigma_{t0}$ se recorren en orden de $\sigma_I/\sigma_{t0}$ decreciente:

1. **Propagación.** Si el elemento está en una punta de una grieta $k$, se activa en $k$. Está en una punta si $\theta_k$ cambia de signo en él y comparte con un elemento agrietado de $k$ un lado que la isolínea cruza. La partición sale de $\theta_k$. $\mathbf n$ es la normal de Rankine, orientada hacia el nodo solitario. Si $\cos\psi \le 0$, se aplica la regla del vértice (§5).
2. **Nucleación.** Si el elemento no comparte ningún nodo con un elemento agrietado, nace una grieta nueva. La normal es la de Rankine y la partición la de la recta por el centroide. Si el nodo que esa recta aísla da $\cos\psi \le 0$, la recta pasa al otro lado del nodo intermedio y aísla el del extremo opuesto: de las dos particiones con esa normal, al menos una tiene $\cos\psi > 0$ (comprobado en 200 000 triángulos al azar). Es la regla del vértice aplicada a la raíz, donde no hay vecinos que puedan repartirse el corte.
3. **Exclusión.** Si comparte un nodo con un elemento agrietado pero no está en una punta, no se activa en ese paso. Dos grietas distintas no nacen en elementos contiguos: la separación mínima es de un anillo de elementos.

Tras cada activación se resuelve de nuevo el campo de la grieta afectada y se sigue con la lista. Una grieta puede avanzar varios elementos en un paso. **La propagación tiene prioridad sobre la nucleación**: se recorren primero las puntas hasta que ninguna puede avanzar; sólo entonces nace una grieta (la del elemento más solicitado que pueda nacer) y se vuelve a propagar. Sin esta prioridad, un elemento dos posiciones por delante de la punta, también por encima de $\sigma_{t0}$ en el mismo paso, hacía nacer otra grieta sobre la misma trayectoria (Diálogo, 2026-11-10). Crece desde la raíz en los dos sentidos y termina en el contorno. Un elemento en las puntas de dos grietas se asigna a la de menor índice: la unión de grietas y las dos discontinuidades por elemento quedan fuera de alcance.

La actualización al inicio del paso es la estándar. Parvaneh y Foster (2016) muestran que con pasos pequeños rinde como la actualización dentro de las iteraciones, más cara.

### 5. Regla del vértice (decisión 3)

Si el nodo aislado $i^*$ de un elemento cortado por $k$ da $\cos\psi = \mathbf n\cdot\mathbf n_e \le 0$, se cambia el signo de $\theta_k$ en $i^*$: pasa a ser una condición de Dirichlet de magnitud $10^{-3}h_e$. Después se recalcula la partición de todos los elementos que comparten $i^*$. La regla se aplica a toda la trayectoria prevista de $k$ (los elementos intactos cortados por la isolínea y unidos a la grieta por lados cruzados), no sólo a la punta: así la trayectoria queda corregida antes de que la grieta llegue. **En un elemento de la punta el nodo aislado está siempre en el lado por el que entra la grieta, ya fijado por ella, y no puede cambiar de lado**; cambia entonces el nodo libre del elemento que deje $\cos\psi > 0$, y la trayectoria pasa al otro lado de ese vértice. La normal $\mathbf n$ se orienta hacia el nodo aislado respecto del punto medio de los dos cruces de la isolínea con los lados del elemento. La trayectoria pasa entonces por ese vértice, la astilla deja de estar cortada y los vecinos se reparten de forma coherente. La desviación es, como mucho, la altura de la astilla (mediana: 10 % del elemento). Cada nodo cambia de lado como mucho una vez por paso, para evitar el ciclo observado en un caso con la regla ingenua. Si después queda algún $\cos\psi \le 0$, ese elemento no se activa en el paso y se emite un aviso.

Datos (2026-11-10, elemento real, 5 mallas Delaunay con $\cos\psi < 0$ y 3 con $\cos\psi = 0$):

| Tratamiento | Mallas en que la franja se separa | Energía / $G_F W t$ |
|---|---|---|
| $l_d$ con su signo | 3 de 5 | se bloquea: 5.5–6.1 |
| $\lvert\cos\psi\rvert$ (código anterior) | 5 de 5 | 1.004–1.021 |
| no activar ese elemento | 0 de 3 | 3.1–7.0 |
| **regla del vértice** | **8 de 8** | **1.0014–1.0015** |

---

## Formulación numérica

### 6. Discretización del problema auxiliar
- Elementos lineales sobre la malla de los `CST_Embedded2D`, uno por grieta.
- Matriz $\mathbf K_k = \sum_e A_e\,\nabla N^\top(\mathbf t_e\mathbf t_e^\top + \varepsilon\mathbf I)\nabla N$, simétrica y definida positiva con las condiciones de Dirichlet.
- Resolución directa dispersa. El costo es un sistema escalar por grieta en cada paso con activaciones; los pasos sin activaciones no resuelven nada.

### 7. Interacción con el elemento
El seguimiento crea el `DiscontinuityState` de cada elemento que activa: normal, nodo solitario y $l_d = (A/h)\cos\psi > 0$. El criterio de Rankine y la irreversibilidad no cambian. El elemento deja de elegir su nodo solitario por su cuenta.

---

## Contrato de implementación

```yaml
name: GlobalCrackTracking
kind: crack_tracker
status: validated

interface:
  region: "todos los CST_Embedded2D del dominio (por omisión) o los ids dados"
  hook: "una vez por paso, desde prepare_step de los elementos de la región (gancho P5 del ADR 0020), antes de su propia activación"

parameters:
  - { name: epsilon, type: float, required: false, default: 1.0e-6, desc: "Regularización adimensional del problema auxiliar" }

validity:
  - una discontinuidad por elemento
  - grietas en 2D con elementos CST_Embedded2D
  - actualización al inicio del paso (pasos de carga moderados cerca de la propagación)

out_of_scope:
  - unión y ramificación de grietas (dos grietas que se encuentran sin ser colineales reparten de forma distinta su lado común: el seguimiento lo avisa, side_conflicts)
  - reorientación de la parte ya formada de una grieta
  - 3D

acceptance:
  verification:
    - name: trayectoria_con_campo_conocido
      setup: "sólo el seguimiento, con un campo de direcciones prescrito cuyas trayectorias exactas son arcos de círculo; mallas Delaunay h = 1/16 … 1/64"
      expect: "sin ruido, la desviación máxima de la isolínea respecto del arco es ≤ 0.02·h y disminuye con h; con ruido de 5° por elemento la isolínea sigue continua de borde a borde (la comparación con el seguimiento local se descartó: Diálogo, 2026-11-10)"
    - name: sensibilidad_a_epsilon
      setup: "lo mismo con ε = 1e-8, 1e-6, 1e-2"
      expect: "la trayectoria cambia menos de 0.1·h entre 1e-8 y 1e-6; con 1e-2 se aparta del arco más de h"
    - name: rotacion_de_direcciones
      setup: "media grieta formada; se gira 30° el campo de direcciones"
      expect: "la parte formada no cambia de partición; la nueva sigue la dirección nueva desde la punta, sin conflictos de lado"
  specific:
    - name: franjas_no_alineadas
      setup: "franja 0.1×0.05 en tracción, plano débil en x = 0.0537 (σ_t0 = σ_0(0.85 + 0.3·y/W + 3·|x − x_0|/L)), banda embebida de ±0.02 alrededor; mallas estructurada y sesgada 16×8 y Delaunay 24×12 (dos semillas); elemento real, NonlinearSolver, 400 pasos"
      expect: "0 conflictos de lado; una grieta que separa la franja (F_fin < 1e-3·F_pico, daño 1 en todos sus elementos); Σl_d/W ∈ [0.99, 1.012]; trabajo/(G_F·W·t) ∈ [0.995, 1.012] (el exceso es Σl_d/W − 1, efecto de las normales de Rankine: Diálogo); 0 fallos del Newton local"
    - name: astillas
      setup: "las 5 primeras mallas Delaunay 16×8 en que la recta x = x_c deja un cos ψ < −0.01, con PrescribedCrackPath x = x_c (regla del vértice) y embebidos los elementos cortados y su anillo"
      expect: "todas se separan, 0 conflictos, 0 avisos y 0 fallos del Newton local; trabajo/(G_F·W·t) ∈ [0.999, 1.012]"
    - name: dos_grietas
      setup: "franja con dos planos débiles (x = 0.03 y 0.07) y resistencia creciente con la altura, cargada hasta poco después de la nucleación"
      expect: "dos grietas con raíces en los dos planos, 0 conflictos de lado"

references:
  - "Oliver J., Huespe A. E., Samaniego E., Chaves E. W. V. (2004). Continuum approach to the numerical simulation of material failure in concrete. Int. J. Numer. Anal. Meth. Geomech. 28(7-8):609-632. doi:10.1002/nag.365"
  - "Oliver J., Huespe A. E. (2004). Continuum approach to material failure in strong discontinuity settings. CMAME 193(30-32):3195-3220. doi:10.1016/j.cma.2003.07.013"
  - "Oliver J., Linero D. L., Huespe A. E., Manzoli O. L. (2008). Two-dimensional modeling of material failure in reinforced concrete by means of a continuum strong discontinuity approach. CMAME 197(5):332-348. doi:10.1016/j.cma.2007.05.017"
  - "Riccardi F., Kishta E., Richard B. (2017). A step-by-step global crack-tracking approach in E-FEM simulations of quasi-brittle materials. Eng. Fract. Mech. 170:44-58. doi:10.1016/j.engfracmech.2016.11.032"
  - "Parvaneh S. M., Foster C. D. (2016). On numerical aspects of different updating schedules for tracking fracture path in strain localization modeling. Eng. Fract. Mech. 152:26-57. doi:10.1016/j.engfracmech.2015.11.011"
  - "Retama Velasco J. (2010). Tesis doctoral, UNAM, cap. 6 (l_d)."
```

---

## Implementación

- Archivo: `solidum/user/discontinuities/crack_tracking.py`
- Clase: `GlobalCrackTracking` (registrada en `CrackTrackerRegistry`, `SPEC_KIND = "crack_tracker"`, sección YAML `crack_trackers`); variante `PrescribedCrackPath` (spec propia).
- Funciones puras del problema auxiliar: `triangle_gradients`, `tracking_matrix`, `solve_tracking_field`, `partition`, `root_partition`, `isoline_segments`, `polyline_signed_distance`.
- Conexión con el elemento: `CST_Embedded2D(..., crack_tracker=None)` (referencia YAML `crack_tracker: <id>`). Sin referencia, cada región conexa de elementos embebidos recibe un `GlobalCrackTracking` con los valores por omisión (`default_tracker`, agrupación por nodos compartidos). `CST_Embedded2D.prepare_step` delega en `tracker.prepare_step(elemento, U)`: la primera llamada de cada ronda (o con otro campo) ejecuta `update`, que es idempotente (itera hasta que nada cambia). El elemento activa con `activate(normal, solitary_node, point, crack)`, que rechaza `l_d ≤ 0`.
- Esfuerzos principales vectorizados con la `B` del CST y `batch_matrix()` del bulk (iguales a `CST_Embedded2D.principal_stress`, prueba); si ningún intacto supera σ_t0, el paso no resuelve nada.
- Diagnóstico: `unactivated_warnings` (punta con cos ψ ≤ 0 sin nodo libre), `side_conflict_warnings` y `side_conflicts()` (grietas que se encuentran), `cracks` (`Crack`: raíz, elementos, valores de Dirichlet por nodo, `side(nodo)`), `field(k)`.
- Constantes en `constants.py`: `TRACKING_EPSILON = 1e-6`, `TRACKING_MIN_DISTANCE = 1e-3` (relativa a $h_e = \sqrt{2A_e}$).
- Tests: `tests/user/discontinuities/test_disc_crack_tracking.py`.

---

## Diálogo

- **2026-11-10** · Decisiones del usuario 1–3 (§0). Borrador para su aprobación. Lo que la IA propone y el usuario confirma o corrige:
  - el campo de direcciones (§1);
  - la regularización y las condiciones de Dirichlet sobre la grieta formada (§2);
  - la nucleación y la exclusión de un anillo (§4);
  - la actualización al inicio del paso (§4).
- **2026-11-10** · **Aprobada por el usuario** («si») junto con los cambios a `CST_Embedded2D`. Implementada y validada el mismo día. Hallazgos y ajustes (los físicos, para que el usuario los revise):
  - **ε = 10⁻⁶ en vez de 10⁻⁴.** Arcos exactos, Delaunay: con 10⁻⁴ la desviación de la isolínea no baja al refinar (2.9e-3, 1.8e-3, 2.9e-3 con h = 1/16, 1/32, 1/64; 0.18·h en la fina); con 10⁻⁶ baja con h (0.0145·h, 0.0079·h, 0.0101·h). Entre 10⁻⁶ y 10⁻⁸ la trayectoria apenas cambia.
  - **Con ruido en las direcciones, el global no es más exacto que el local.** Con 5° de ruido por elemento, la desviación media del global respecto del arco es 1.9 a 3.6 veces la del seguimiento local, para cualquier ε (h = 1/16 a 1/64). Explicación: el ruido equivale, en promedio, a una difusión transversal de orden σ² ≈ 8·10⁻³, que en trayectorias curvas desplaza la isolínea como una ε grande. Esto contradice el prototipo del 2026-11-10 que sostuvo la recomendación de la decisión 2, que medía otra cosa: el seguimiento local se detenía en vértices en algunas mallas. La ventaja que sí se confirma es la coherencia: 0 conflictos y separación completa en todas las mallas no alineadas, donde el reparto elemental se bloqueaba. El criterio «≤ la mitad del local» se sustituyó por la continuidad de la isolínea.
  - **Propagación antes que nucleación** (§4). Sin ella, en una franja completamente embebida nacían 3 a 17 grietas paralelas en el paso del pico y Newton divergía.
  - **Paso de carga.** Con la activación al inicio del paso, todo elemento que supera σ_t0 en el paso del pico se activa a la vez. Con pasos de 8 % de σ_t0 aparecían grietas pequeñas de 1 a 5 elementos junto a la principal; con pasos de 2 a 4 % desaparecen. Es la condición «pasos de carga moderados» de `validity`.
  - **Energía con normales de Rankine.** Con pasos finos el balance es exacto: trabajo / (G_F·Σl_d·t) = 1.0001–1.0002. Σl_d/W vale 0.999–1.000 en mallas estructuradas, pero 1.0025–1.0087 en Delaunay: cada elemento usa su propia normal de Rankine, que gira entre 2.5° y 4° de un vecino a otro, y entonces la suma telescópica de $l_d$ no da exactamente W. Con la normal de la trayectoria (prototipo del 2026-11-10) daba 1.0000. Las tolerancias de energía de las franjas pasan a [0.995, 1.012]. Queda para el usuario si la normal debe ser la de Rankine (tesis) o la de la trayectoria; va ligado a la decisión del modo mixto.
  - **Regla del vértice en la punta** (§5): el nodo aislado de un elemento de la punta ya está fijado por la grieta; cambia de lado el nodo libre. Las 5 mallas con astillas se separan con 0 fallos.
  - **Unión de grietas** (fuera de alcance): dos grietas nacidas en los dos bordes del mismo plano se unen bien si son colineales (malla estructurada, energía 1.0006). Si no lo son (Delaunay), pasan una junto a otra y reparten de forma distinta su lado común; el seguimiento lo avisa.
  - **Banda delgada de grieta prescrita.** En el ensayo de Van Vliet de la tesis, la banda de un solo elemento de espesor no deja seguir la isolínea cuando las direcciones principales giran cerca de la entalla. La grieta se detenía y nacían otras. Ese modelo prescribe la grieta, así que usa la variante [`PrescribedCrackPath`](PrescribedCrackPath.md), que reproduce exactamente el resultado anterior.
- **2026-11-10** · Decisión 4 del usuario: **se mantiene la normal de Rankine** (la de la tesis). Con compresión lateral y deslizamiento, la normal de la trayectoria cambia la energía ≤ 0.3 % (spec `CST_Embedded2D`, Diálogo).
