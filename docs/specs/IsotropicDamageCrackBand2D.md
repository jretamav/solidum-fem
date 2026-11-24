# ISOTROPIC DAMAGE CRACK BAND 2D — daño escalar regularizado por banda de grieta, con integración implícita o IMPL-EX

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> **Variante de componente existente** (Reglas §4): subclase de `IsotropicDamage2D` ([spec](IsotropicDamage2D.md), validada), que no se modifica. Esta spec documenta **sólo lo que cambia**: (i) la velocidad de ablandamiento $\alpha$ deja de ser un dato y se calcula en cada punto de integración a partir de la energía de fractura $G_f$ y de la longitud característica $h$ del elemento (banda de grieta, Bažant y Oh 1983), para que la energía disipada por unidad de área de grieta sea $G_f$ con cualquier malla; (ii) una medida de deformación equivalente más, la de **Rankine**; (iii) una segunda forma de integrar el daño, **IMPL-EX** (Oliver, Huespe y Cante 2008); (iv) un tope de daño propio. Todo lo demás (esfuerzo nominal y efectivo, estado 3D asociado, variable histórica, Kuhn-Tucker, ley exponencial, medidas existentes, tangente consistente) es el del padre.
>
> Uso previsto: agrietamiento de estructuras de concreto con fisuración distribuida, en 2D y 3D (la variante 3D está en [IsotropicDamageCrackBand3D](IsotropicDamageCrackBand3D.md)).
>
> **Estado: aprobada por el usuario (2026-11-11, «si», con las cinco recomendaciones Q1–Q5) y validada.** Bloque 1 (banda de grieta, medida de Rankine, tope propio y siembra de $h$ en los sólidos 2D y 3D) y bloque 2 (integración IMPL-EX, §8, con la pieza genérica del pseudo-tiempo, [ADR 0034](../adr/0034-pseudo-tiempo-de-los-solvers-incrementales.md), y la prueba estructural de Van Vliet) implementados el mismo día.
>
> **Bloque 3 (2026-11-17): ancho de banda según la orientación de la grieta — estudiado y cerrado sin programar** (usuario, 2026-11-23; §4, «Ancho según la orientación»; Q6–Q11 del Diálogo). El estudio se conserva: con la grieta libre ningún ancho quita la dependencia de la malla, porque lo que queda viene de la cinemática de la banda, no de la fórmula de $h$. Límites medidos del modelo en §10.

---

## Especificación física

### 0. Descripción general

Daño isótropo escalar en pequeñas deformaciones para sólidos 2D (tensión plana o deformación plana), con ablandamiento exponencial **regularizado por banda de grieta**: la rama de ablandamiento de cada punto de integración se ajusta a su longitud característica $h$ para que el material disipe $G_f/h$ por unidad de volumen. Cuando la deformación se localiza en una fila de elementos de ancho $h$, la energía disipada por unidad de área de grieta es $G_f$, independiente del tamaño de los elementos. Es la forma estándar de quitar la dependencia de malla de la energía en los modelos de fisuración distribuida (Bažant y Oh 1983; Rots 1988; Oliver 1989).

`STRAIN_DIM = 3`, Voigt del proyecto $[\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]$ con $\gamma_{xy} = 2\varepsilon_{xy}$.

### 1. Lo que no cambia respecto al padre

$$\boldsymbol\sigma = (1-d)\,\bar{\boldsymbol\sigma}, \qquad \bar{\boldsymbol\sigma} = \mathbf C_e\,\boldsymbol\varepsilon, \qquad \kappa_{n+1} = \max\bigl(\kappa_n,\ \varepsilon_\text{eq}(\boldsymbol\varepsilon_{n+1})\bigr),$$

$$d(\kappa) = \begin{cases} 0 & \kappa \le \kappa_0 \\[2pt] 1 - \dfrac{\kappa_0}{\kappa}\,e^{-\alpha(\kappa-\kappa_0)} & \kappa > \kappa_0 \end{cases}$$

con $\bar{\boldsymbol\sigma}$ el esfuerzo efectivo y las medidas de deformación equivalente evaluadas sobre el estado 3D asociado (padre, §3). En tracción uniaxial de esfuerzo todas las medidas dan $\varepsilon_\text{eq} = \varepsilon$, de modo que

$$\kappa_0 = \frac{f_t}{E}, \qquad \sigma = (1-d)\,E\,\kappa = f_t\,e^{-\alpha(\kappa-\kappa_0)} \quad (\kappa > \kappa_0).$$

El parámetro de entrada pasa a ser la resistencia a tracción $f_t$ (más natural para concreto que $\kappa_0$); $\kappa_0$ se deriva.

### 2. Regularización por banda de grieta

**Energía por unidad de volumen** en tracción uniaxial monótona hasta la separación completa (la rama elástica más la de ablandamiento; al final el material no almacena energía, así que toda es disipada):

$$g = \int_0^\infty \sigma\,d\varepsilon = \tfrac12\,f_t\,\kappa_0 + \int_{\kappa_0}^{\infty} f_t\,e^{-\alpha(\kappa-\kappa_0)}\,d\kappa = \tfrac12\,f_t\,\kappa_0 + \frac{f_t}{\alpha}.$$

**Condición de banda de grieta**: $g = G_f/h$. Despejando,

$$\boxed{\;\alpha(h) = \frac{1}{\dfrac{G_f}{h\,f_t} - \dfrac{\kappa_0}{2}}\;}$$

que exige $G_f/h > \tfrac12 f_t\kappa_0$, es decir,

$$\boxed{\;h < h_\text{max} = \frac{2\,E\,G_f}{f_t^{2}} = 2\,\ell_\text{ch}\;}, \qquad \ell_\text{ch} = \frac{E\,G_f}{f_t^{2}} \ \text{(Hillerborg, Modéer y Petersson 1976)}.$$

Con $h \ge h_\text{max}$ la rama de ablandamiento tendría retroceso en el propio material (la energía elástica almacenada en el pico ya supera $G_f/h$): el material lo rechaza con `ValueError` al recibir $h$ (`initial_point_state`), con $h$, $h_\text{max}$ y los parámetros en el mensaje, igual que `ManderConcrete1D`. Para concreto usual ($E = 30$ GPa, $f_t = 3$ MPa, $G_f = 100$ N/m): $\kappa_0 = 10^{-4}$, $\ell_\text{ch} = 0.333$ m, $h_\text{max} = 0.667$ m.

**Energía disipada en un estado intermedio** (carga monótona hasta $\kappa$; la parte almacenada es $\tfrac12\sigma\kappa$, recuperable al descargar por la secante):

$$D(\kappa) = \tfrac12 f_t\kappa_0 + \frac{f_t}{\alpha}\Bigl(1 - e^{-\alpha(\kappa-\kappa_0)}\Bigr) - \tfrac12\,f_t\,e^{-\alpha(\kappa-\kappa_0)}\,\kappa, \qquad D(\infty) = g = \frac{G_f}{h}.$$

Es la fórmula de referencia de las pruebas de energía en el punto material.

### 3. Tope del daño (pregunta Q5)

El padre satura el daño en $d_\text{max} = $ `DAMAGE_MAX` $= 0.999$. Con banda de grieta ese tope **no es inocuo**: se alcanza a mitad del ablandamiento y, a partir de ahí, el punto transmite $(1-d_\text{max})\,E\,\varepsilon$, que **crece** con la apertura. Cálculo exacto en tracción uniaxial con el concreto usual de arriba ($D_\text{cap}$ = energía disipada al llegar al tope; $w = \kappa_\text{cap}\,h$ la apertura equivalente):

| $h$ [mm] | $\alpha$ | $d_\text{max} = 0.999$: $w$ al tope [mm] | $D_\text{cap}\,h/G_f$ | $d_\text{max} = 1-10^{-6}$: $D_\text{cap}\,h/G_f$ |
|---|---|---|---|---|
| 5 | 151 | 0.067 | 0.732 | 0.998 |
| 10 | 305 | 0.083 | 0.816 | 0.999 |
| 25 | 779 | 0.104 | 0.894 | 0.999 |
| 50 | 1622 | 0.120 | 0.935 | 1.000 |
| 100 | 3529 | 0.133 | 0.962 | 1.000 |
| 400 | 30000 | 0.118 | 0.994 | 1.000 |

Con $d_\text{max} = 0.999$ y $h = 5$ mm se pierde el 27 % de $G_f$, y con una apertura de 0.1 mm la banda transmite todavía $0.20\,f_t$ (con $h = 10$ mm, $0.10\,f_t$): la energía deja de ser objetiva respecto a la malla, que es justo lo que la variante pretende. **Propuesta**: tope propio de la variante, parámetro `damage_max` con valor por omisión $1 - 10^{-6}$ (pierde menos del 0.25 % de $G_f$ para $h \ge 5$ mm; la rigidez residual $10^{-6}E$ evita la singularidad exacta de la matriz). El padre conserva `DAMAGE_MAX`. Mismo criterio que la ley cohesiva `CohesiveDamageIsotropic`, que no trunca el daño por la misma razón.

### 4. Longitud característica $h$ (pregunta Q2)

$h$ es el ancho de la banda en que se localiza la deformación cuando la grieta atraviesa el elemento. Opciones:

- **(a) Convenio de Abaqus** («línea típica a través del elemento»): $h = \sqrt{A_e}$ en elementos lineales (`Tri3`, `Quad4`) y $h = \tfrac12\sqrt{A_e}$ en cuadráticos (`Tri6`, `Quad8`, `Quad9`), con $A_e$ el área del elemento; el mismo $h$ en todos sus puntos de integración.
- **(b) Triángulos lineales de Rots (1988)**: $h = \sqrt{2A_e}$ en `Tri3`; una grieta que cruza una malla de triángulos atraviesa una banda más ancha que $\sqrt{A_e}$.
- **(c) Longitud proyectada de Oliver (1989)**: depende de la dirección de la grieta al iniciarse el daño; es la más precisa para bandas oblicuas a la malla, pero exige fijar esa dirección en el punto (una variable interna más) y cambia de valor según el ángulo.

**Recomendación**: (a), y que la prueba de objetividad con mallas de triángulos mida $D/(G_f\,A_\text{grieta})$; si en `Tri3` se aparta más de un 5 % de 1, pasar a (b) en triángulos lineales. (c) queda fuera de alcance.

**Decisión y resultado (2026-11-11).** Aprobada (a) con la regla de los triángulos. La prueba de objetividad (barra con una columna de cuadrados debilitada, partidos en dos triángulos) dio con $\sqrt{A_e}$ en `Tri3` $D/(G_f A) = 1.41$ (la banda la forman los dos triángulos de cada cuadrado: su ancho es $\sqrt{2A_e}$), así que `Tri3` usa **(b) $h = \sqrt{2A_e}$** y da 1.00, igual que `Quad4`. En 3D lo mismo con `Tet4` (seis tetraedros por cubo): $\sqrt[3]{V_e}$ daba 1.81 y $\sqrt[3]{6V_e}$ da 1.00. Los demás elementos usaban el convenio de Abaqus.

**Cuadráticos (2026-11-11, tras el diagnóstico de la viga SEN).** La misma prueba ampliada a los cuadráticos, con el ½ de Abaqus, disipaba $1.97\,G_fA$ (`Quad8`, `Quad9`, `Hex20`, `Hex27`), $2.68\,G_fA$ (`Tri6`) y $3.23\,G_fA$ (`Tet10`), igual con todas las mallas: con el esfuerzo uniforme de la barra la deformación se localiza en el elemento entero. Aplicada la regla del 5 % (usuario, «si»): cada cuadrático usa el ancho del elemento lineal de su celda, $c\,\Omega_e^{1/d}$ (`CHARACTERISTIC_LENGTH_FACTOR`): `Quad4`/`Quad8`/`Quad9` 1, `Tri3`/`Tri6` $\sqrt2$; `Hex8`/`Hex20`/`Hex27` 1, `Tet4`/`Tet10` $6^{1/3}$.

**Factor de localización de los cuadráticos (2026-11-12).** En flexión los cuadráticos disipaban sólo $0.69$–$0.80\,G_fA$ (viga de Petersson, objetivos entre mallas; `Quad4` 0.94). Diagnóstico por punto de Gauss: el elemento lineal tiene deformación constante a través de la banda y se ablanda entero; el cuadrático, más flexible, tiene deformación lineal dentro del elemento y, con gradiente, el ablandamiento se localiza en una parte de sus puntos: en `Quad8` en dos de sus tres filas de Gauss $3\times3$ en todos los elementos de la banda (fracción $13/18 = 0.72$), en `Tri6` en dos de sus tres puntos o en los tres (fracción $2/3$ o 1, según el triángulo). Es el resultado de Jirásek y Bauer (2012, §4: en 1D, una fila de dos con 2 puntos, dos (a veces una) de tres con 3 puntos) y de Slobbe, Hendriks y Rots (2013). Por eso
$$h = \alpha\,c\,\Omega_e^{1/d},$$
con $\alpha$ la fracción del elemento que se ablanda, declarada por el elemento según su regla (`CRACK_BAND_LOCALIZATION`):

| Elemento | Regla | $\alpha$ |
|---|---|---|
| `Tri3`, `Quad4`, `Tet4`, `Hex8` | cualquiera | 1 |
| `Quad8`, `Quad9` | Gauss $3\times3$ | $13/18$ |
| `Quad8`, `Quad9` | Gauss $2\times2$ | $1/2$ |
| `Tri6`, `Hex20`, `Hex27`, `Tet10` | cualquiera | 1, sin factor establecido |

Comprobación en una viga entallada a flexión con malla estructurada alineada ($D/G_fA$ en $\delta = 2$ mm): `Quad4` 1.015; `Quad8` 0.773 sin $\alpha$ y 1.029 con $13/18$. En `Tri6` la fracción varía de un triángulo a otro: sin factor 0.79–0.87 (disipa de menos, del lado seguro) y con el $2/3$ de la literatura 1.17–1.28 (de más); queda con $\alpha = 1$ (usuario). En 3D no hay un factor medido ni publicado: $\alpha = 1$. Con un cuadrático Solidum avisa una vez y describe la limitación, sin recomendar un orden de elemento (usuario, 2026-11-12): la corrección sólo es fiable con la banda alineada con la malla (Jirásek y Bauer: con la malla inclinada $13/18$ da una respuesta demasiado dúctil). Con esfuerzo uniforme en el elemento (la barra de la prueba de objetividad) el cuadrático se ablanda entero y $\alpha < 1$ sobrestima la energía: $G_fA/\alpha = 1.385\,G_fA$ en `Quad8` y `Quad9`.

#### Ancho según la orientación de la grieta (bloque 3, 2026-11-17; estudio cerrado sin programar el 2026-11-23)

> Lo que sigue es el estudio del bloque 3 tal como se hizo. No se implementó: Solidum usa el $h$ de arriba. El usuario cerró el bloque el 2026-11-23 (Diálogo) con las preguntas Q6–Q11 sin decidir.

**Problema.** $h = \alpha\,c\,\Omega_e^{1/d}$ es el ancho de la banda sólo cuando la banda sigue las líneas de la malla. Una grieta inclinada respecto de la malla cruza una **escalera** de elementos: la banda mínima que separa el cuerpo son los elementos que la grieta corta (dos elementos que sólo comparten un nodo no separan nada), y su área por unidad de longitud de grieta es mayor que el lado del elemento. Con $h$ fijo la energía disipada depende del ángulo entre la grieta y la malla. Jirásek y Bauer (2012, §5.3 y §7) lo cuantifican: el criterio de la raíz del área puede equivaler a equivocar $G_f$ en un 50 % o más.

**Candidatos.** Los dos usan la dirección $\mathbf n_c$ normal a la grieta (Q7) y dan un solo $h$ por elemento; los dos reproducen el $h$ actual cuando la grieta sigue los lados de los elementos en mallas estructuradas. No cuando sigue las diagonales de una malla de triángulos: ahí $h_P = \Delta/\sqrt2$ y el actual es $\Delta$ (el 0.682 de `Tri3` a 45° en la tabla).

- **(P) Proyección**, la «Pm» de Jirásek y Bauer (2012, §5.3.2) cuando la dirección se toma en el centro del elemento; ellos la atribuyen a Červenka y Pukl (1994), Červenka, Pukl y Eligehausen (1995) y Bažant y Planas (1998), y la fórmula del rectángulo a Bažant (1985) (atribuciones por verificar con el bibliotecario). El ancho del elemento medido en la dirección $\mathbf n_c$,
  $$h_P = \max_a\,(\mathbf x_a\cdot\mathbf n_c) - \min_a\,(\mathbf x_a\cdot\mathbf n_c),$$
  sobre los nodos de vértice $\mathbf x_a$. Sentido geométrico: el área de la escalera de elementos que cruza una grieta recta, por unidad de longitud de grieta, es el ancho proyectado de un elemento. En un cuadrado de lado $\Delta$ con la grieta a un ángulo $a$ de sus lados, $h_P = \Delta\,(\cos a + \sin a)$.
- **(O2) Segundo método de Oliver** (Oliver 1989; Jirásek y Bauer 2012, ec. 26):
  $$h_{O2} = \frac{A_e}{l_d}, \qquad l_d = \oint_{\Gamma_e} \varphi\,(\mathbf n_c\cdot\mathbf n_e)\,ds = \int_{\Omega_e} \nabla\varphi\cdot\mathbf n_c\,d\Omega,$$
  con $\varphi = 1$ en los vértices de un lado de la recta que pasa por el centro del elemento perpendicular a $\mathbf n_c$, $\varphi = 0$ en los del otro, interpolada con las funciones de forma de los vértices. $l_d$ es la cuerda de la grieta que pasa por el centro, proyectada. En la CST es el $l_d = (A_e/h_s)\cos(\theta - \alpha)$ de la tesis (cap. 6; el de `CST_Embedded2D`). En el cuadrado, $h_{O2} = \Delta/\cos a$.

Coinciden a 0° y a 45° en cuadrados; entre medias $h_{O2} < h_P$ (a 22.5°, $1.082\,\Delta$ frente a $1.307\,\Delta$).

**Geometría (agente teórico, 2026-11-17, `informes/teorico/verificacion_teorica_2026-11-17_banda_orientacion.md`).** Una recta con normal $\mathbf n$ corta un elemento convexo si y sólo si su posición cae en un intervalo de longitud $h_P$; promediando sobre la posición de la recta, el área de los elementos cortados por unidad de longitud de grieta es el promedio de $h_P$ (pesado por el área). Por eso, si cada elemento cortado se ablanda entero y disipa $(G_f/h)\,A_e$, **con $h = h_P$ la energía por unidad de longitud de grieta es exactamente $G_f$**, en cualquier malla de elementos convexos, 2D o 3D: $h_P$ es el área del elemento entre la cuerda **media** de las rectas que lo cortan. En O2, $l_d = \Delta\cos a$ (la distancia proyectada entre los puntos medios de los dos lados que cruza la recta por el centro) es mayor que la cuerda media, $\Delta/(\cos a + \sin a)$, y $h_{O2}$ resulta menor que el ancho real de la banda. Con rectas aleatorias en mallas estructuradas, $\sum (A_e/h)/\ell$: P 0.9997–1.0003; O2 en cuadrados 1.141, 1.207 y 1.183 a 10°, 22.5° y 30°.

**Contrapeso: sensibilidad a la dirección.** P depende de $\mathbf n_c$ a primer orden cerca de la alineación y O2 a segundo orden: con la banda alineada con la malla y $\mathbf n_c$ desviada 10° (por ejemplo, por el cortante de los `Quad4` que se agrietan en flexión, Jirásek y Bauer Fig. 20b), $h_P = 1.158\,\Delta$ (energía −14 %) y $h_{O2} = 1.015\,\Delta$ (−1.5 %). Jirásek y Bauer lo observan en sus corridas (§6.1 con la dirección por punto de Gauss; §6.4, banda en el eje de simetría: P y Pm son los que más se apartan y O2 = O1m el que menos). Con la dirección exacta y la malla de cuadrados sólo probaron 0° y 45°, donde coinciden. Jirásek y Bauer descartan además el factor de orientación de Červenka et al. (1995), $\gamma = 1$ a 1.5, que reduciría la energía hasta un 33 % (§6.3, §7).

**Medición (prototipo, 2026-11-17; `informes/prototipos_2026-11-17/banda_orientacion/`, no versionado).** Barra en tracción de 200 × 100 mm, espesor 10, celdas de 5 mm, malla estructurada girada un ángulo $a$ respecto de la barra; $E = 30\,000$ MPa, $\nu = 0$, $f_t = 3$ MPa, $G_f = 0.1$ N/mm, integración implícita, 800 pasos hasta 0.4 mm (carga final < 0.4 % del pico). La grieta se fuerza perpendicular a la carga: sólo los elementos que corta la recta $x = 102.5$ mm pueden dañarse (los demás son elásticos, para que la grieta no empiece en los bordes escalonados), y $\mathbf n_c = \mathbf e_x$ se impone. $D$ es la energía disipada al final; $\ell$, la longitud exacta de la recta dentro de la malla. Valores de $D/(G_f\,t\,\ell)$; el exacto es 1:

| $a$ | `Quad4` actual | `Quad4` P | `Quad4` O2 | `Tri3` actual | `Tri3` P | `Tri3` O2 |
|---|---|---|---|---|---|---|
| 0° | 0.998 | 0.998 | 0.998 | 0.998 | 0.998 | 0.998 |
| 10° | 1.141 | 0.989 | 1.124 | 0.959 | 0.974 | 0.945 |
| 22.5° | 1.264 | 0.977 | 1.172 | 0.900 | 0.974 | 0.831 |
| 30° | 1.364 | 1.009 | 1.188 | 0.868 | 1.002 | 1.185 |
| 45° | 1.341 | 0.962 | 0.962 | 0.682 | 0.964 | 0.964 |

En todos los casos se rompieron todos o casi todos los elementos cortados (por ejemplo, 26 de 27 en `Quad4` a 22.5°) y ninguno fuera de ellos. **Lectura: la prueba es casi circular.** Con la banda forzada a los elementos cortados y la dirección impuesta, $D$ es la identidad geométrica del párrafo anterior ($\sum A_e/h_e$ sobre los elementos rotos entre $\ell$: 0.968 frente a 0.977 medido en `Quad4` 22.5° P); la corrida sólo añade que cada elemento cortado se ablanda entero. La dispersión de ±4 % es de los bordes. **Con la banda libre** (franja debilitada de ±30 mm, dirección impuesta; corridas del agente teórico, sólo indicativas): `Quad4` P 0.907 y 0.869 a 10° y 22.5°, O2 1.032 y 1.044, actual 1.048 y 1.128; a 30°, con P se formaron dos bandas (1.69). La banda libre rompe menos elementos que la escalera recta, y P disipa de menos. 

**Viga entallada con la grieta libre y la dirección estimada (prototipo, 2026-11-17; `informes/prototipos_2026-11-17/banda_orientacion/viga_entallada.py`).** La viga de `notched_beam_bending` (flexión en cuatro puntos, 410 × 100 × 50 mm, entalla de 10 × 50 mm, concreto de Petersson, IMPL-EX con 1000 pasos hasta δ = 2 mm) con la dirección estimada como propone Q7 (deformación principal en el centro del elemento con el último paso confirmado sin daño, fija desde que el elemento se daña). Mallas: `S` estructurada alineada (la de la prueba existente) y cuatro no estructuradas de gmsh: `U` (10 mm, Frontal-Delaunay), `Ua` (10 mm, otro algoritmo), `U75` (7.5 mm) y `U5` (5 mm). $D/(G_f\,A_\text{lig})$ en δ = 2 mm (carga final < 5 % del pico):

| Malla | `Quad4` actual | `Quad4` P | `Quad4` O2 | `Tri3` actual | `Tri3` P | `Tri3` O2 |
|---|---|---|---|---|---|---|
| `S` | 1.015 | 1.015 | 1.015 | 1.038 | 1.004 | 1.074 |
| `U` | 1.110 | 1.011 | 1.091 | 1.277 | 1.222 | 1.260 |
| `Ua` | 1.454 | 1.207 | 1.419 | 1.240 | 1.145 | 1.217 |
| `U75` | 1.079 | 0.945 | 1.046 | 1.285 | 1.248 | 1.321 |
| `U5` | 1.019 | 0.904 | 1.003 | 0.924 | 0.886 | 0.943 |
| error medio / máximo | 16 % / 45 % | **11 % / 25 %** | 15 % / 42 % | | | |

(Error medio y máximo de $|D/G_fA - 1|$ sobre las diez mallas de primer orden.) Segundo orden, serie aparte: `Quad8` `S` 1.029 / 1.028 / 1.026 y `U` 1.220 / 1.069 / 1.038; `Tri6` 0.79–0.82 con los tres (sin factor $\alpha$, como se sabía).

Lectura: (1) la proyección tiene el menor error medio y la menor dispersión entre mallas, y siempre da la menor energía de las tres; (2) **ningún método quita la dependencia de la malla**: la dispersión entre mallas (0.89–1.25 con P) la causan la grieta que se desvía o se ensancha a dos elementos, el daño parcial de los vecinos y daño espurio bajo un punto de carga (mapas `dano_viga_*.png`), que no dependen del $h$ de los elementos rotos; lo mismo hallan Slobbe, Hendriks y Rots (2013) con la proyección en mallas giradas (0.5–1.9 de $G_f$); (3) **la proyección falla donde avisó el agente teórico**: en `U5` (gmsh dio una malla casi regular) la grieta sube por una columna fuera del eje de simetría, la dirección estimada se tuerce por el cortante y P quita un 10 % de energía (O2 y el actual, ≈ 1.0).

**Literatura y códigos** (`informes/bibliotecario/busqueda_banda_orientacion_2026-11-17.md`). La proyección es la opción «Govindjee» de DIANA (Theory Manual §37.12.1; también ATENA, multiplicada por el factor $\gamma$ de Červenka); Abaqus, ANSYS, Code_Aster y OpenSees usan un $h$ geométrico sin orientación. El segundo método de Oliver no está en ningún código de la lista. Ninguno de los ocho textos de Reglas §1 trata un ancho según la orientación (no verificado del todo: de Borst et al. 2012, cap. 6, no se pudo leer). Criterio de admisión: la proyección cumple (ii) y no (i), así que el criterio no decide y decide el usuario; O2 no cumple ninguna de las dos, así que iría a un módulo de usuario. Oliver (1989) no se pudo leer; su antecedente (Oliver y Fernández Roure 1987) confirma la fórmula de la longitud equivalente.

**Propuesta que se estudió** (cada punto, una pregunta del Diálogo; cerradas sin decidir el 2026-11-23):

- **Q6. Método**: P, en el programa principal (cumple la condición de los códigos: DIANA), y descartar O2 (no está en ningún código y en la viga no mejoró al $h$ actual). Con su límite declarado: no quita la dependencia de la malla y puede quitar ~10 % de energía cuando la banda sigue las líneas de la malla fuera de un eje de simetría.
- **Q7. Dirección y momento**: $\mathbf n_c$ = dirección de la deformación principal mayor **en el plano** ($\tan 2\theta = \gamma_{xy}/(\varepsilon_{xx} - \varepsilon_{yy})$, con $\gamma$ de Voigt), evaluada en el centro del elemento ($\boldsymbol\xi = \mathbf 0$ en cuadriláteros; baricentro en triángulos), con la deformación del último paso confirmado sin daño en el elemento, y fija desde el paso en que el primero de sus puntos pasa el umbral (se confirma con el paso: un paso repetido tras una bisección da el mismo $h$); el mismo $h$ en todos los puntos del elemento. Con daño isótropo el esfuerzo efectivo es coaxial con la deformación (en el plano, $\sigma_1 - \sigma_2 = E(\varepsilon_1 - \varepsilon_2)/(1+\nu)$), así que da igual tomar la dirección del esfuerzo principal. Jirásek y Bauer: tomar la dirección en cada punto de Gauss por separado da peores resultados que tomarla en el centro (§6.1). Oliver (1989) la obtiene de un análisis de localización; con daño isótropo y Rankine coincide con la deformación principal mayor **en tracción biaxial con $\sigma_2/\sigma_1 \ge \nu$** (Jirásek y Bauer §5.3.3, con Jirásek y Horák 2010); fuera de ese caso (tracción uniaxial con $\nu > 0$, u otras medidas de deformación equivalente) la regla es una aproximación. **Falta una regla para deformaciones principales iguales**, donde la dirección no está definida y $h_P$ de un cuadrado va de $\Delta$ a $\sqrt2\,\Delta$ (decisión del usuario; propuesta: con diferencia relativa de las principales menor que un umbral, el ancho medio del elemento, perímetro/π en 2D, 1.273 Δ en el cuadrado).
- **Q8. Por omisión**: el método de Q6 para todos los sólidos, con un parámetro `crack_band_width = "projection"` u `"oliver"` (por omisión, según Q6) | `"element"` (el $h$ actual, para comparar con códigos que usan el tamaño del elemento). Cambian los resultados de las mallas en que la grieta no sigue la malla (Van Vliet con `Tri3` no estructurada); con mallas estructuradas y la grieta alineada, no.
- **Q9. Cuadráticos**: $h = \alpha\,h_P$ con la proyección de los nodos de vértice y el $\alpha$ actual. La limitación de los cuadráticos con la banda inclinada (Jirásek y Bauer, §6.5) no la resuelve la orientación; el aviso sigue igual.
- **Q10. Límite $h_{\max}$**: hoy se comprueba al construir. Con $h$ según la orientación se comprueba al construir con el mayor ancho posible del elemento (la mayor distancia entre dos vértices), que es una cota exacta de $h_P$ en cualquier dirección (con $\alpha$, $\alpha$ por esa distancia): una malla que pasa esa comprobación nunca falla a mitad del análisis. Efecto (decisión del usuario, Reglas §3): el tamaño de elemento admisible baja por $\sqrt2$ en `Quad4` y `Tri3` estructurados y por $\sqrt3$ en `Hex8` y `Tet4` (con $h_{\max} = 667$ mm, $\Delta < 471$ y $< 385$ mm), y un elemento alargado se rechaza por su diagonal aunque la grieta corra a lo largo. Con O2 la cota es otra.
- **Q11. Aceptación**: (i) verificación geométrica: la barra girada de la tabla con la banda forzada y $\nu = 0$, `Quad4` y `Tri3` con las dos diagonales, `Hex8`/`Tet4` girados alrededor de $z$ en los dos sentidos y con una normal genérica, p. ej. $(1, -1, 1)/\sqrt3$: $D/(G_f\,A) = 1 \pm 0.05$; (ii) $h$ y $\mathbf n_c$ contra fórmulas cerradas en elementos girados y distorsionados, y una prueba de cortante puro (dirección a 45°) que blinda el factor $\gamma = 2\varepsilon$; (iii) validación (agente validador): viga entallada con la malla girada y con malla no estructurada frente a la malla alineada, energía a ±5 %; con el método por omisión, volver a correr `notched_beam_bending` y `van_vliet_tension`.

**3D.** P vale tal cual en 3D: el ancho del elemento en la dirección $\mathbf n_c$ (de la deformación principal mayor), con los vértices del hexaedro o del tetraedro. Con la malla alineada da el $h$ actual (`Hex8` $\sqrt[3]{V_e}$, `Tet4` de la partición en seis tetraedros $\sqrt[3]{6V_e}$). Ver [IsotropicDamageCrackBand3D](IsotropicDamageCrackBand3D.md), §2.

### 5. Medida de Rankine (pregunta Q1)

Se añade a las tres existentes (`energy`, `mazars`, `modified_von_mises`):

$$\varepsilon_\text{eq} = \frac{\langle \bar\sigma_I\rangle_+}{E},$$

con $\bar\sigma_I$ el mayor esfuerzo principal **efectivo** del estado 3D asociado, $\bar{\boldsymbol\sigma}^\text{3D} = \mathbf C^\text{3D}_e\,\boldsymbol\varepsilon^\text{3D}$ y $\boldsymbol\varepsilon^\text{3D} = [\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, \gamma_{xy}, 0, 0]$, $\varepsilon_{zz} = c_{zz}(\varepsilon_{xx}+\varepsilon_{yy})$ (padre, §3). En tensión plana $\bar\sigma_{zz} = 0$; en deformación plana $\bar\sigma_{zz} = \nu(\bar\sigma_{xx} + \bar\sigma_{yy})$. Así

$$\bar\sigma_I = \max\Bigl(\tfrac12(\bar\sigma_{xx}+\bar\sigma_{yy}) + \sqrt{\tfrac14(\bar\sigma_{xx}-\bar\sigma_{yy})^2 + \bar\sigma_{xy}^2},\ \bar\sigma_{zz}\Bigr).$$

Calibración: tracción uniaxial de esfuerzo $\Rightarrow \varepsilon_\text{eq} = \varepsilon$ ($\kappa_0 = f_t/E$, como las demás). Tracción biaxial igual: el daño empieza en $\sigma = f_t$ (criterio de esfuerzo máximo). Cortante puro: en $\tau = f_t$. **Compresión uniaxial o biaxial: $\varepsilon_\text{eq} = 0$, nunca daña.**

**Derivada** (para la tangente consistente). Con $\mathbf n_I$ la dirección principal de $\bar\sigma_I$ y, en Voigt de esfuerzo, $\mathbf p = [n_x^2,\ n_y^2,\ 2n_xn_y]$ (porque $\bar\sigma_I = \mathbf n_I\cdot\bar{\boldsymbol\sigma}\,\mathbf n_I = n_x^2\bar\sigma_{xx} + n_y^2\bar\sigma_{yy} + 2n_xn_y\bar\sigma_{xy}$):

$$\frac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon} = \frac{1}{E}\,\mathbf C_e^{\top}\mathbf p = \frac{1}{E}\,\mathbf C_e\,\mathbf p \quad (\bar\sigma_I \text{ en el plano}),
\qquad
\frac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon} = \frac{\lambda}{E}\,[1,\ 1,\ 0] \quad (\bar\sigma_I = \bar\sigma_{zz},\ \text{deformación plana}),$$

con $\mathbf C_e$ la matriz plana del material y $\lambda$ el primer parámetro de Lamé. La tangente de carga activa es la del padre con esta derivada:

$$\mathbf C^\text{alg} = (1-d)\,\mathbf C_e - (1-d)\Bigl(\frac1\kappa + \alpha\Bigr)\,\bar{\boldsymbol\sigma}\otimes\frac{\partial\varepsilon_\text{eq}}{\partial\boldsymbol\varepsilon}.$$

No es simétrica en general ($\bar{\boldsymbol\sigma}$ y $\mathbf C_e\mathbf p$ sólo son paralelos en estados uniaxiales de esfuerzo efectivo). Con dos esfuerzos principales efectivos iguales y positivos, $\bar\sigma_I$ tiene una arista: se toma cualquiera de las dos direcciones (subgradiente); el Newton puede perder la convergencia cuadrática en ese estado, sin consecuencias en la respuesta.

**Qué hace cada medida en compresión uniaxial de esfuerzo** (con banda de grieta, todo daño ablanda con la energía $G_f$ de tracción):

| Medida | Inicio del daño en compresión uniaxial | Con $\nu = 0.2$, $f_t = 3$ MPa |
|---|---|---|
| `energy` | $\lvert\sigma\rvert = f_t$ | 3 MPa: el concreto se rompe en compresión igual que en tracción |
| `mazars` | $\lvert\sigma\rvert = f_t/(\sqrt2\,\nu)$ | 10.6 MPa ($\approx 3.5 f_t$), con la energía de tracción |
| `modified_von_mises` | $\lvert\sigma\rvert = k\,f_t = f_c$ | $f_c$, pero ablanda con $G_f$ de tracción: compresión demasiado frágil |
| `rankine` | nunca | elástico en compresión |

### 6. Requisitos de interfaz (arquitectura; la decide la IA, ADR nuevo para la pieza de IMPL-EX)

1. **Longitud característica por punto de integración.** La pieza genérica ya existe (`Material.uses_characteristic_length`, `Material.initial_point_state(h)`, `Element._seed_characteristic_lengths`), pero hoy sólo la usan las barras y los marcos (con `ManderConcrete1D`). Los elementos sólidos 2D `Tri3`, `Quad4`, `Tri6`, `Quad8` y `Quad9` deben sembrar $h$ en cada punto de integración según §4. El material guarda $h$ en su estado (variable `h`) y calcula $\alpha(h)$ en cada llamada. Un elemento que no siembre $h$ con este material produce un error claro al construir el modelo, no un resultado silencioso.
2. **Razón de incrementos de pseudo-tiempo para IMPL-EX** (§8): $r = \Delta\tau_{n+1}/\Delta\tau_n$. Pieza genérica del [ADR 0034](../adr/0034-pseudo-tiempo-de-los-solvers-incrementales.md): el solver declara el incremento de su parámetro de control en cada intento de paso (`Assembler.prepare_all_steps(U, pseudo_time_increment=Δτ)` o, si lo conoce tras el predictor, `Assembler.set_pseudo_time_increment(Δτ)`); el ensamblador calcula $r$ respecto del último paso confirmado y la entrega al material (`set_pseudo_time_ratio`) y a los parámetros de su familia por lotes. Fuera de un paso (análisis lineales, post-proceso, commit final) el material recibe «sin paso» y congela su daño en el confirmado. Un solver que abre un intento sin declarar su pseudo-tiempo se rechaza al confirmar el primer paso, con un error que nombra los solvers que sí lo declaran.

---

## Formulación numérica

### 7. Integración implícita (`integration = 'implicit'`, por omisión)

La del padre (§6–§9 de su spec) con $\alpha = \alpha(h)$ del punto, el tope `damage_max` de la variante y la medida elegida (incluida Rankine, §5). Tangente consistente; se pasa a la secante en descarga, sin daño y al llegar al tope.

### 8. Integración IMPL-EX (`integration = 'impl-ex'`)

Oliver, Huespe y Cante (2008): la variable histórica se **extrapola** con los dos últimos pasos convergidos para calcular el esfuerzo, y se **actualiza implícitamente** para la extrapolación siguiente. Dentro del paso el material es lineal (daño congelado), así que la tangente es la secante, simétrica y definida positiva.

**Pseudo-tiempo** $\tau$: el parámetro de control del solver (factor de carga, desplazamiento prescrito, longitud de arco o la magnitud del control indirecto). Su incremento en cada paso se conoce antes de resolverlo; $r = \Delta\tau_{n+1}/\Delta\tau_n$.

Dado el estado convergido $(\kappa_n, \kappa_{n-1}, \tilde d_n)$ y $\boldsymbol\varepsilon_{n+1}$:

1. **Extrapolación**: $\tilde\kappa_{n+1} = \kappa_n + r\,(\kappa_n - \kappa_{n-1})$.
2. **Daño del paso** (irreversible): $\tilde d_{n+1} = \max\bigl(\tilde d_n,\ d(\tilde\kappa_{n+1})\bigr)$, con $d$ la ley de §1 con $\alpha(h)$ y el tope de §3.
3. **Esfuerzo y tangente**: $\boldsymbol\sigma_{n+1} = (1-\tilde d_{n+1})\,\mathbf C_e\,\boldsymbol\varepsilon_{n+1}$, $\ \mathbf C^\text{tan} = (1-\tilde d_{n+1})\,\mathbf C_e$ (exacta: $\tilde d_{n+1}$ no depende de $\boldsymbol\varepsilon_{n+1}$).
4. **Actualización implícita** (en el mismo `compute_state`, para el estado de prueba): $\kappa_{n+1} = \max\bigl(\kappa_n,\ \varepsilon_\text{eq}(\boldsymbol\varepsilon_{n+1})\bigr)$, $\ \kappa^\text{prev}_{n+1} = \kappa_n$. Al confirmarse el paso quedan como estado convergido.

Consecuencias:

- Si no hay otras no linealidades, el problema global de cada paso es lineal y el Newton converge en una iteración (más la de comprobación).
- **Disipación no negativa**: $\mathcal D_{n+1} = \psi_0(\boldsymbol\varepsilon_{n+1})\,(\tilde d_{n+1} - \tilde d_n) \ge 0$ por el paso 2. El $\max$ del paso 2 evita una «curación» espuria cuando una extrapolación anterior se pasó y el punto descarga; en Oliver et al. (2008) la monotonía se razona sobre la variable extrapolada (por confirmar con el artículo si la imponen de la misma forma).
- **Retraso de un paso**: el daño entra un paso después de que $\varepsilon_\text{eq}$ supere $\kappa_0$ (al inicio $\kappa_n = \kappa_{n-1} = \kappa_0$). **Orden**: la extrapolación lineal es exacta si $\kappa(\tau)$ es lineal y su error local es $O(\Delta\tau^2)$ donde $\kappa(\tau)$ es suave; en un quiebre de $\kappa(\tau)$ (inicio del daño, descarga, recarga) el error es $O(\Delta\tau)$. Como $\kappa$ se actualiza de forma implícita, el error no se acumula. En la práctica manda el quiebre del inicio del daño: cerca del pico hacen falta pasos pequeños (prueba de Van Vliet: el pico se pasa un 7.7 % con 0.5 µm y un 0.8 % con 0.0625 µm).
- **Bisección**: repetir un paso con otro $\Delta\tau_{n+1}$ sólo cambia $r$; el estado convergido no cambia, así que el resultado es idempotente para el mismo intento.
- **Fuera de un paso** ($r$ «sin paso»): $\tilde d = \tilde d_n$ y el estado no cambia. Reevaluar el estado confirmado en el desplazamiento convergido reproduce el esfuerzo del paso (post-proceso, análisis lineales sobre el estado de una etapa, commit final de `solidum.run`).
- **Indicador de error** por punto: $e_{n+1} = \lvert d(\kappa_{n+1}) - \tilde d_{n+1}\rvert$, en la variable de estado `damage_error` (post-proceso por punto o promedio por elemento). El control adaptativo del paso con este indicador (Oliver et al. 2008, §5) queda fuera de alcance.
- **Newton**: sin otras no linealidades el problema de cada paso es lineal; el corrector hace una resolución y la de comprobación (corrección nula), dos por paso.

### 9. Simetría

- Implícito: como el padre (simétrica sólo con `energy`); con `rankine`, `mazars` y `modified_von_mises`, no simétrica.
- IMPL-EX: siempre simétrica y definida positiva mientras $\tilde d < 1$ (lo garantiza el tope). La instancia declara `IS_SYMMETRIC = True` con IMPL-EX.

### 10. Caveats numéricos

- Con la implícita, el ablandamiento hace la tangente indefinida y el Newton puede fallar cerca del pico o con varias bandas compitiendo; es la razón de ofrecer IMPL-EX.
- Con IMPL-EX la robustez se paga con exactitud: el resultado depende del tamaño del paso (primer orden). La prueba de aceptación lo mide.
- **Sesgo de malla** (medido, §4): la banda tiende a seguir las líneas de la malla; con grietas oblicuas a la malla se ensancha o zigzaguea y la energía se aparta de $G_f$. Viga entallada a flexión con la grieta libre (prototipo con este material, 2026-11-17), diez mallas `Quad4` y `Tri3` (una alineada y cuatro no estructuradas de cada uno): $D/(G_f A)$ = 0.92–1.45, error medio 16 % y máximo 45 %; con la malla alineada, 1.015 (`Quad4`) y 1.038 (`Tri3`). Ningún ancho según la orientación de la grieta lo quita (§4, bloque 3): lo que queda es cinemático (elementos que la grieta cruza y no se rompen) y, en triángulos, el daño isótropo con estados mixtos (disipa también la energía de cortante).
- **Estructuras con cortante: la biela comprimida se degrada** (medido, viga SEN de Schlangen, 2026-11-11, `informes/reportes/reporte_viga_schlangen_banda_2026-11-11.pdf` y [TensionCompressionDamage2D](TensionCompressionDamage2D.md) §1). El daño escalar quita también la rigidez de la biela comprimida que cruza la banda: con mallas finas la banda baja casi vertical desde la entalla, en lugar de ir en diagonal hacia la placa de carga, y disipa 1.54 y 1.98·$G_f t L$ con $h/4$ y $h/8$. No es objetiva. El daño sólo de tracción del modelo 2 da la grieta diagonal con las mismas mallas.
- **Modo mixto** (medido por el validador, panel en L de Winkler, 2026-11-18, `informes/reportes/reporte_modelo2_panel_L_2026-11-18.pdf`; Quad4 en mallas libres de 10, 5 y 2.5 mm). Pico objetivo (2.7 % entre mallas); en el ablandamiento, la carga en $\delta$ = 1 mm es 309 / 317 / 397 N (26 % entre mallas) y la energía hasta $\delta$ = 2 mm difiere 7.7 %.
- Transferencia espuria de esfuerzos: los elementos vecinos de la banda, dañados en parte, siguen transmitiendo esfuerzo a través de ella en estados mixtos. Los rotos ($d \to 1$) no transmiten nada: 0–0.3 MPa en el panel en L, ≈ 0.02 MPa calculado (`informes/teorico/verificacion_teorica_2026-11-18_cortante_grieta_abierta.md`).
- Rankine: arista con esfuerzos principales efectivos iguales (§5).

---

## Contrato de implementación

```yaml
name: IsotropicDamageCrackBand2D
kind: material
status: validated          # bloques 1 y 2 (2026-11-11)

interface:
  strain_dim: 3
  primary_state_var: damage
  is_symmetric: false      # cota de clase; la instancia es true con equivalent_strain='energy' o con integration='impl-ex'
  batch_kernel: true       # estado [κ, d, h, κ_prev, e]; integración y razón de pseudo-tiempo en batch_params
  uses_pseudo_time: "true con integration='impl-ex' (ADR 0034)"

admission:                # Reglas §1; variante de un material admitido (IsotropicDamage2D)
  textbooks:
    - de_borst_crisfield: "§6.2, daño isótropo escalar con ablandamiento exponencial (el padre); que el texto trate la regularización por banda de grieta está por confirmar con el texto"
    - de_souza_neto: "cap. 12, mecánica del daño (el padre); la banda de grieta en el texto, por confirmar"
  general_codes:
    - abaqus: "Concrete Damaged Plasticity y Brittle Cracking: ablandamiento en tracción por energía de fractura con la longitud característica del elemento (línea típica a través del elemento; la mitad en elementos cuadráticos)"
    - diana: "Total Strain Crack Model: ablandamiento por energía de fractura con el ancho de banda de grieta h del elemento"
    - code_aster: "ENDO_ISOT_BETON y MAZARS (daño escalar del padre); METHODE='IMPLEX' de STAT_NON_LINE para leyes de daño (doc. R5.03.81)"
    - opensees: "ASDConcrete3D: regularización por energía de fractura (-autoRegularization) e integración IMPL-EX (-implex)"

parameters:
  - { name: E,          type: float, required: true,  desc: "Módulo de Young intacto [Pa]" }
  - { name: nu,         type: float, required: true,  desc: "Coeficiente de Poisson" }
  - { name: f_t,        type: float, required: true,  desc: "Resistencia a tracción [Pa]; κ_0 = f_t/E" }
  - { name: G_f,        type: float, required: true,  desc: "Energía de fractura en modo I [N/m]; el material disipa G_f/h por unidad de volumen" }
  - { name: hypothesis, type: str,   required: false, default: "plane_stress", desc: "plane_stress | plane_strain (como el padre)" }
  - { name: equivalent_strain, type: str, required: false, default: "rankine", desc: "rankine (nueva) | energy | mazars | modified_von_mises (Q1)" }
  - { name: k,          type: float, required: false, default: null, desc: "f_c/f_t, sólo y obligatorio con modified_von_mises (como el padre)" }
  - { name: integration, type: str, required: false, default: "implicit", desc: "implicit | impl-ex (Oliver, Huespe y Cante 2008, §8); impl-ex requiere un solver que declare su pseudo-tiempo (ADR 0034) y pasos pequeños cerca del pico" }
  - { name: damage_max, type: float, required: false, default: 0.999999, desc: "Tope del daño de la variante, en (0, 1); ver §3 y Q5" }
  - { name: density,    type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

signature:
  compute_state: "(ε: ndarray(3,), state_vars) -> (σ: ndarray(3,), C_tan: ndarray(3,3), state_vars')"
  initial_point_state: "(h: float) -> state_vars; ValueError si h no es finito y positivo o si h ≥ h_max = 2·E·G_f/f_t²"
  strain_kind: "plano Voigt — [ε_xx, ε_yy, γ_xy] con γ_xy = 2·ε_xy"

state_schema:
  kappa:      "float ≥ κ_0 — máxima deformación equivalente histórica"
  damage:     "float ∈ [0, damage_max] — daño con que se calcula el esfuerzo"
  h:          "float > 0 — longitud característica del punto de integración, fijada por el elemento (0 sin sembrar: sólo admitido en régimen elástico)"
  kappa_prev: "float — κ del paso confirmado anterior (extrapolación de IMPL-EX; en la implícita, el κ confirmado)"
  damage_error: "float ≥ 0 — indicador |d(κ) − d̃| de IMPL-EX; 0 en la implícita"

conventions:
  sign: "σ > 0 ⇔ tracción"
  voigt: "γ_xy = 2·ε_xy; ε_eq sobre el estado 3D asociado (padre §3)"
  units: "E, f_t en Pa; G_f en N/m; h en m (cualquier sistema coherente)"
  regularization: "α(h) = 1/(G_f/(h·f_t) − κ_0/2), h < h_max = 2·E·G_f/f_t²"

validity:
  - "E > 0, f_t > 0, G_f > 0, ν ∈ (−1, 0.5)"
  - "0 < h < 2·E·G_f/f_t² en cada punto (se comprueba al recibir h)"
  - "damage_max ∈ (0, 1)"
  - "integration = 'impl-ex' sólo con un solver que declare su pseudo-tiempo (§6.2, ADR 0034)"
  - "pequeñas deformaciones; agrietamiento dominado por tracción"
  - "sólo en elementos que proporcionan h (sólidos 2D y 3D, barras, marcos); otro elemento lo rechaza al construirse"

out_of_scope:
  - "daño separado de tracción y compresión; energía de fractura en compresión"
  - "cierre de grietas con recuperación de la rigidez (efecto unilateral)"
  - "h según la dirección de la grieta (proyección, segundo método de Oliver): estudiado en el bloque 3 (§4) y cerrado sin programar (usuario, 2026-11-23)"
  - "control adaptativo del paso con el indicador de error de IMPL-EX"
  - "anisotropía del daño; grieta fija o rotante (siguiente paso posible: grieta rotante de deformación total)"
  - "regularización no local o de gradiente"

numerical_caveats:
  - "Implícita: tangente indefinida en ablandamiento; posibles fallos del Newton cerca del pico."
  - "IMPL-EX: error de primer orden en los quiebres de κ(τ) (inicio del daño) y de segundo donde κ(τ) es suave; el daño entra un paso tarde; pico sobrepasado si el paso es grueso."
  - "Sesgo de la banda por la malla y transferencia espuria de esfuerzos (§10)."
  - "Rankine: arista con esfuerzos principales efectivos iguales; subgradiente."

acceptance:
  material_point:
    - name: alpha_closed_form
      setup: "h ∈ {5, 25, 100, 400} mm con E = 30 GPa, f_t = 3 MPa, G_f = 100 N/m"
      expect: "α del estado = 1/(G_f/(h·f_t) − κ_0/2)"
      tol_rel: 1.0e-14
    - name: peak_and_onset
      setup: "tracción uniaxial de esfuerzo monótona, cada medida"
      expect: "d = 0 hasta ε = κ_0 = f_t/E; σ máximo = f_t"
      tol_rel: 1.0e-12
    - name: dissipated_energy_per_volume
      setup: "tracción uniaxial monótona hasta d = damage_max, por deformación, pasos finos; h ∈ {5, 25, 100, 400} mm"
      expect: "∫σ dε − ½σε = D(κ) de §2 en estados intermedios; al final D = (G_f/h)·(fracción del tope, §3), > 0.997·G_f/h con damage_max = 1 − 1e-6"
      tol_rel: 1.0e-6
    - name: h_max_rejected
      setup: "initial_point_state con h = h_max y con h = 1.01·h_max; h ≤ 0, NaN"
      expect: "ValueError con h, h_max y los parámetros en el mensaje"
    - name: rankine_values
      setup: "tensión plana y deformación plana: tracción uniaxial, tracción biaxial igual, cortante puro, compresión uniaxial y biaxial"
      expect: "ε_eq = ε; inicio en σ = f_t; inicio en τ = f_t; compresión: d = 0 con |ε| = 100·κ_0"
      tol_rel: 1.0e-12
    - name: consistent_tangent_fd
      setup: "carga activa, cada medida (incluida rankine), tensión y deformación plana, estados generales; diferencia centrada"
      expect: "|C_FD − C_tan| ≤ 1e-6·max|C_FD| (implícita)"
      tol_rel: 1.0e-6
    - name: restriction_of_3d
      setup: "historia de carga, descarga y recarga; las cuatro medidas; ambas hipótesis; mismo h"
      expect: "κ, d y σ en el plano iguales a IsotropicDamageCrackBand3D con el estado 3D asociado"
      tol_rel: 1.0e-12
    - name: batch_kernel_matches
      expect: "el kernel por lotes reproduce bit a bit σ, C_tan y la fila de estado (incluida h) de compute_state, implícita e IMPL-EX"
  impl_ex:
    - name: secant_tangent_exact
      expect: "C_tan = (1 − d̃)·C_e exacta, simétrica y definida positiva"
      tol_rel: 1.0e-14
    - name: first_order_convergence
      setup: "camino monótono de tracción uniaxial desde κ_0 con pasos Δ, Δ/2, Δ/4, Δ/8: lineal en el pseudo-tiempo (quiebre en el inicio del daño) y cuadrático (κ(τ) suave)"
      expect: "el error máximo de σ frente a la implícita se reduce a la mitad al reducir el paso a la mitad en el lineal (razón 2 ± 15 %) y a la cuarta parte en el cuadrático (4 ± 15 %)"
    - name: no_negative_dissipation
      setup: "carga, descarga parcial y recarga con extrapolación que se pasa"
      expect: "d̃ nunca decrece; disipación por paso ≥ 0"
    - name: bisection_idempotent
      setup: "mismo intento de paso dos veces tras una bisección"
      expect: "estado de prueba idéntico"
  elements:
    - name: characteristic_length_seeded
      setup: "Tri3, Quad4, Tri6, Quad8, Quad9 con una malla regular y otra distorsionada"
      expect: "h de cada punto = α·√A_e (Quad4; Quad8 y Quad9 con α = 13/18 en Gauss 3×3 y ½ en 2×2) o √(2A_e) (Tri3, Tri6), decisiones de §4"
      tol_rel: 1.0e-12
    - name: localization_warning
      expect: "un aviso por elemento y regla con los cuadráticos: con factor (Quad8, Quad9) o sin factor establecido (Tri6); ninguno con los lineales"
    - name: missing_seed_error
      expect: "un elemento que no siembre h con este material produce un error claro al construir el modelo"
  mesh_objectivity:
    - name: bar_in_tension
      setup: "barra o placa en tracción con un debilitamiento en una fila; 1, 2, 4 y 8 elementos por unidad de longitud en la dirección de la carga; Quad4 y Tri3 estructurados; implícita e IMPL-EX"
      expect: "energía disipada al final = G_f·A_grieta y la curva fuerza–alargamiento posterior al pico igual en todas las mallas, ±1 % (Quad4, Tri3, Tri6); Quad8 y Quad9 con α = 13/18, G_f·A/α (esfuerzo uniforme: el elemento se ablanda entero), objetivo entre mallas"
      tol_rel: 1.0e-2
    - name: notched_beam_bending
      setup: "viga entallada a flexión en cuatro puntos, H = 100 mm, malla estructurada alineada de celdas de 10 mm, entalla de la celda central hasta H/2; IMPL-EX, 1000 pasos hasta δ = 2 mm; Quad4, Quad8, Tri3, Tri6"
      expect: "D/(G_f·A_ligamento) = 1 ± 6 % con Quad4 y Tri3; Quad8 (α = 13/18) disipa lo mismo que Quad4 ±5 % y su pico ±8 %; Tri6 sin factor disipa entre 0.65 y 0.9 veces lo de Tri3"
      tol_rel: 5.0e-2
  structural:
    - name: van_vliet_tension
      setup: "ensayo de tracción de Van Vliet de la tesis (§8.1; modelo de tests/user/discontinuities/_van_vliet.py) con fisuración distribuida, dos mallas; implícita e IMPL-EX con dos tamaños de paso"
      expect: "energía disipada final frente a G_f·área del ligamento ±3 %; carga máxima entre mallas ±3 %; IMPL-EX con el paso fino a ±2 % de la implícita; balance de energía; la curva de la tesis sólo como referencia"
      tol_rel: 3.0e-2
    - name: schlangen_demo
      setup: "opcional, demostración: viga de Schlangen sobre la malla de FEAP del usuario"
      expect: "sin criterio cuantitativo; respuesta global, trayectoria de la banda y robustez del cálculo"

references:
  - "Bažant Z.P., Oh B.H. (1983). Crack band theory for fracture of concrete. Matériaux et Constructions 16, 155–177."
  - "Hillerborg A., Modéer M., Petersson P.-E. (1976). Analysis of crack formation and crack growth in concrete by means of fracture mechanics and finite elements. Cem. Concr. Res. 6, 773–781."
  - "Rots J.G. (1988). Computational modeling of concrete fracture. Tesis doctoral, TU Delft."
  - "Oliver J. (1989). A consistent characteristic length for smeared cracking models. IJNME 28, 461–474. doi:10.1002/nme.1620280214 (no leído; bloque 3)"
  - "Oliver J., Fernández Roure M. (1987). Rev. Int. Métodos Numéricos para Cálculo y Diseño en Ingeniería 3(4), 343–368 (longitud característica y longitud equivalente de grieta; antecedente de Oliver 1989)."
  - "Govindjee S., Kay G.J., Simo J.C. (1995). Anisotropic modelling and numerical simulation of brittle damage in concrete. IJNME 38, 3611–3633. doi:10.1002/nme.1620382105 (proyección; informe UCB/SEMM-94/18 leído)"
  - "Červenka V., Pukl R., Eligehausen R. (1995). Mesh sensitivity effects in smeared finite element analysis of concrete fracture. FraMCoS-2, 1387–1396 (proyección con factor de orientación γ)."
  - "DIANA FEA, Theory Manual, cap. 37, §37.12.1 Crack Bandwidth (Rots; Govindjee): https://manuals.dianafea.com/d102/Theory/Theorych37.html"
  - "Jirásek M., Bauer M. (2012). Numerical aspects of the crack band approach. Computers & Structures 110–111, 60–78. doi:10.1016/j.compstruc.2012.06.006"
  - "Slobbe A.T., Hendriks M.A.N., Rots J.G. (2013). Systematic assessment of directional mesh bias with periodic boundary conditions: Applied to the crack band model. Eng. Fract. Mech. 109, 186–208. doi:10.1016/j.engfracmech.2013.06.005"
  - "Oliver J., Huespe A.E., Cante J.C. (2008). An implicit/explicit integration scheme to increase computability of non-linear material and contact/friction problems. CMAME 197, 1865–1889."
  - "Code_Aster, doc. R5.03.81, méthode IMPLEX: https://ericca.uqtr.ca/fr14.4/man_r/r5/r5.03.81.pdf"
  - "OpenSees, ASDConcrete3D: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/ASDConcrete3D.html"
  - "Las del padre (Simó y Ju 1987; Mazars 1986; de Vree et al. 1995; de Borst et al. 2012; de Souza Neto et al. 2008)."
```

---

## Implementación

- Archivo: `solidum/materials/damage_crack_band_2d.py` (lo común con la 3D en `solidum/materials/_crack_band.py`; Rankine en `solidum/materials/_equivalent_strain.py`, código `RANKINE`, que los padres no admiten).
- Clase: `IsotropicDamageCrackBand2D(CrackBandMixin, IsotropicDamage2D)`; reutiliza el núcleo compilado del padre (`_compute_damage_2d`) con $\alpha(h)$ del punto y el tope `damage_max`.
- Siembra de $h$: pieza genérica del elemento base (`Element.CHARACTERISTIC_LENGTH_FACTOR`, `Element.CRACK_BAND_LOCALIZATION`, `Element._seed_geometric_characteristic_length`): $h = \alpha\,c\,\Omega_e^{1/d}$ con $\Omega_e = \sum_g w_g \det J_g$ de la propia cinemática compilada del elemento y $\alpha$ según la regla de cuadratura del elemento (§4), con un aviso único por elemento y regla en los cuadráticos; un elemento sin `CHARACTERISTIC_LENGTH_FACTOR` ni `SEEDS_CHARACTERISTIC_LENGTH` rechaza el material al construirse.
- Tests:
  - `tests/test_damage_crack_band.py` (punto material 2D y 3D: α(h), h_max, pico, energía $D(\kappa)$, tope, Rankine, tangente por diferencias finitas, restricción exacta del 3D, kernel por lotes bit a bit; siembra de $h$ en los diez sólidos, mallas regulares y distorsionadas, con el factor de localización por regla y su aviso; elemento sin $h$).
  - `tests/validation/test_crack_band_objectivity.py` (barra en tracción con los diez sólidos, mallas del generador estructurado; viga entallada a flexión con `Quad4`, `Quad8`, `Tri3` y `Tri6`).
- Notas de traducción: la energía de fractura se da por unidad de área ($N/m$ en SI); con unidades N–mm, $G_f$ en N/mm.

---

## Diálogo

- **2026-11-11** · Spec creada por la IA a pedido del usuario («¿qué otra alternativa podríamos implementar para simular el agrietamiento del concreto… smeared cracking?»; aprobó con «si» preparar esta opción). Clasificación: la banda de grieta es estándar (Abaqus, DIANA, Code_Aster, OpenSees; textos del padre, la sección concreta por confirmar). Arquitectura (IA): subclase del padre validado, que no cambia; los elementos sólidos 2D/3D siembran $h$ con la pieza genérica existente; la razón de pseudo-tiempo de IMPL-EX es una pieza genérica nueva con ADR propio.

**Preguntas al usuario** (una decisión física por pregunta; cada una con recomendación):

- **Q1. Medida de deformación equivalente para el concreto.** Con banda de grieta todo daño ablanda con la energía de tracción $G_f$, así que la medida decide qué estados «agrietan» (tabla de §5): `energy` rompe el concreto en compresión a $f_t$; `mazars` a $\approx 3.5 f_t$ (con $\nu = 0.2$); `modified_von_mises` a $f_c$, pero con un ablandamiento de compresión demasiado frágil; `rankine` nunca daña en compresión. **Recomendación: `rankine` por omisión**, para agrietamiento dominado por tracción, que es el uso previsto; el comportamiento no lineal en compresión queda fuera de esta variante (grieta rotante de deformación total o daño de tracción y compresión, más adelante). Las otras tres siguen disponibles.
- **Q2. Longitud característica $h$.** **Recomendación: convenio de Abaqus** ($\sqrt{A_e}$ en lineales, $\tfrac12\sqrt{A_e}$ en cuadráticos) y decidir los triángulos lineales con la prueba de objetividad: si en `Tri3` la energía se aparta más de un 5 % de $G_f\,A$, usar $\sqrt{2A_e}$ (Rots 1988).
- **Q3. Clasificación de IMPL-EX.** Está implementado en dos códigos generales de la lista (Code_Aster, `METHODE='IMPLEX'`; OpenSees, `ASDConcrete3D -implex`), pero no consta en los textos de referencia de Reglas §1: el criterio no decide y la clasificación es del usuario. **Recomendación: programa principal**, como opción de integración de un material admitido (no cambia la ley, sólo cómo se integra en el tiempo); la alternativa es llevar la variante IMPL-EX a un módulo de usuario.
- **Q4. Aceptación.** Casos y tolerancias del bloque `acceptance`: energía exacta en el punto material; objetividad de malla ±1 % en una barra (Quad4); Van Vliet ±3 % entre mallas y frente a $G_f\,A$; IMPL-EX de primer orden y a ±2 % de la implícita con paso fino. Sin comparación punto a punto con curvas experimentales. **Recomendación: aprobarla así.**
- **Q5. Tope del daño.** Con el `DAMAGE_MAX = 0.999` del padre se pierde entre el 0.6 % ($h = 400$ mm) y el 27 % ($h = 5$ mm) de $G_f$, y la banda abierta 0.1 mm transmite hasta $0.2\,f_t$ (§3). **Recomendación: tope propio de la variante, `damage_max` $= 1 - 10^{-6}$ por omisión** (pérdida < 0.25 % para $h \ge 5$ mm).

- **2026-11-11** · **Usuario: «si»** a la spec con las cinco recomendaciones: Q1 Rankine por omisión (las otras tres disponibles); Q2 convenio de Abaqus y $\sqrt{2A}$ en `Tri3` si la energía se aparta más de un 5 %; Q3 IMPL-EX en el programa principal; Q4 aceptación tal cual; Q5 `damage_max` propio $= 1 - 10^{-6}$. La IA divide el trabajo en dos bloques: el **bloque 1** (este) sin IMPL-EX; el **bloque 2** añade IMPL-EX con su pieza genérica (razón de pseudo-tiempo, ADR nuevo), `integration`, `kappa_prev` y la prueba de Van Vliet.
- **2026-11-11** · **Bloque 1 implementado. Hallazgos para el usuario:**
  - **Q2 resuelta con la prueba**: en la barra con una columna debilitada, `Tri3` con $\sqrt{A_e}$ disipaba $1.41\,G_fA$ (la banda son los dos triángulos de cada cuadrado) y con $\sqrt{2A_e}$ disipa $0.999\,G_fA$, igual que `Quad4`. Aplicada la regla. Lo mismo en 3D con `Tet4` (seis por cubo): $\sqrt[3]{V_e}$ daba $1.81\,G_fA$; con $\sqrt[3]{6V_e}$, el análogo de la regla de Rots, $0.999\,G_fA$, igual que `Hex8` (la spec 3D pedía decidir `Tet4` «igual que `Tri3`»).
  - **Objetividad**: con $\nu = 0$, $D/(G_fA) = 0.9992$, 0.9985, 0.9974 y 0.9955 con 5, 10, 20 y 40 elementos a lo largo de la barra (`Quad4`; idéntico en `Tri3`, `Hex8` y `Tet4`), y la carga en $u = 0.05$ mm difiere menos del 0.4 % entre mallas.
  - **Con $\nu = 0.2$ la barra disipa $(1-\nu^2)\,G_fA = 0.96\,G_fA$ en todas las mallas.** No es error de discretización: las columnas vecinas se descargan e impiden la contracción lateral de la banda ($\varepsilon_{yy} \approx 0$ en ella), Rankine ve $\varepsilon_\text{eq} = \varepsilon_{xx}/(1-\nu^2)$ y el trabajo de $\sigma_{xx}$ en la banda es $(1-\nu^2)\,G_f/h$. La calibración de §2 supone tracción uniaxial de esfuerzo; una banda confinada lateralmente disipa ese 4 % menos. Queda como dato; corregirlo (p. ej. calibrar con la banda en deformación plana) sería una decisión física del usuario.
  - **Van Vliet con integración implícita**: con la malla de la tesis en `Tri3` (1200 elementos, paso de 1 µm) el Newton se detiene a λ ≈ 0.16 (desplazamiento de la rótula 93 µm), con $D \approx 0.75\,G_fA$ disipado: es justo la falta de robustez de la implícita en ablandamiento que motiva IMPL-EX. La prueba estructural `van_vliet_tension` pasa al bloque 2.
- **2026-11-11** · **Usuario: «si»** a dejar el efecto de Poisson documentado y sin corregir: con $\nu = 0.2$ la banda confinada lateralmente disipa $(1-\nu^2)\,G_f = 0.96\,G_f$ por unidad de área en todas las mallas (la objetividad se conserva). Una corrección de la calibración sólo sería exacta para un estado de esfuerzos de la banda, y el 4 % está muy por debajo de la incertidumbre usual de $G_F$. El efecto es propio del daño isótropo; la grieta rotante (modelo 3 de la familia de concreto) no lo tiene.
- **2026-11-11** · **Bloque 2 implementado: IMPL-EX** (ADR 0034). Resultados:
  - **Punto material**: tangente $(1-\tilde d)\,\mathbf C_e$ exacta, simétrica y definida positiva; disipación $\psi_0\,\Delta\tilde d \ge 0$ con carga, descarga y recarga; kernel por lotes bit a bit con la implícita y con IMPL-EX, con razón variable y «sin paso». **Orden**: con un camino lineal en el pseudo-tiempo que arranca en $\kappa_0$ la razón de errores al dividir el paso es 2.0 (primer orden, por el quiebre del inicio del daño); con un camino cuadrático (suave) es 4.0 (segundo orden). Lo de la spec («primer orden») vale en los quiebres; se corrigió §8 y el caso de aceptación.
  - **Global** (barra de la prueba de objetividad, $\nu = 0$): con IMPL-EX el corrector hace dos resoluciones por paso (una y la de comprobación) en los 400 pasos; la implícita, hasta 26 ensamblajes en el pico (bisecciones). Con paso de 0.125 µm la energía disipada es $0.9989$–$0.9996\,G_fA$ con 5, 10 y 20 elementos y la carga en $u = 0.05$ mm coincide con la implícita a 0.3 %; con 0.5 µm la extrapolación pasa el pico y disipa $1.01$–$1.03\,G_fA$.
  - **Van Vliet** (todo el concreto agrietable, `Tri3`, dos mallas de 1200 y 4396 elementos): implícita, pico 26.49 / 26.95 kN (1.7 %) y energía en $u = 0.12$ mm 0.966 / 0.956 J (1.1 %); al final ($u = 1.2$ mm, carga 0.36 % del pico) $D = 1.378$ J $= 0.981\,(1-\nu^2)\,G_fA$. IMPL-EX con 0.0625 µm: pico 26.72 kN (+0.9 %), energía +0.2 %; con 0.5 µm, pico +7.7 %. La implícita ya no se detiene: con las placas de la prueba de la embebida (dos capas de 8 mm) la carga puntual las flexionaba y el concreto junto a ellas llegaba a $\sigma_I$ 1.7 veces mayor que en el cuello, así que la probeta se rompía junto a las placas y el Newton fallaba ahí (el «λ ≈ 0.16» del bloque 1). Con placas de dos capas de 25 mm (supuesto propio, la tesis no da el espesor) la grieta se forma en el cuello. Prueba: `tests/validation/test_crack_band_van_vliet.py`.
  - **Salvaguarda $\tilde d = \max(\tilde d_n, d(\tilde\kappa))$**: implementada como en §8; el artículo de Oliver, Huespe y Cante (2008) no la impone; se conserva por decisión del usuario (ver la entrada siguiente del Diálogo).
- **2026-11-11** · **Contraste con Oliver, Huespe y Cante (2008), CMAME 197:1865–1889, Tabla 3 y §3.6** (artículo conseguido por el usuario; copia local en `docs/referencias/`, no versionada). Coinciden: la extrapolación $\tilde\alpha_{n+1} = \alpha_n + (\Delta t_{n+1}/\Delta t_n)\,\Delta\alpha_n$ (ec. 35), el esfuerzo $\tilde{\boldsymbol\sigma} = (1-\tilde d)\,\mathbf C^e:\boldsymbol\varepsilon$ y la tangente secante (ecs. 38–39), y el inicio del daño con un paso de retraso (error en los puntos no suaves de $\alpha(t)$, §3.6). Difieren en dos detalles: (1) el artículo linealiza la ley con la pendiente del paso anterior, $\tilde q = q_n + H_n\,\Delta\tilde\alpha$ y $\tilde d = 1 - \tilde q/\tilde\alpha$ (ec. 37), mientras Solidum evalúa la ley exponencial exacta en $\tilde\kappa$ (ambas de primer orden; con la exponencial convexa la linealización daña algo más y con pasos grandes puede dar $\tilde q < 0$); (2) el artículo no impone $\tilde d \ge \tilde d_n$ (en ablandamiento su $\tilde d \ge d_n$ implícito, pero puede bajar del extrapolado anterior); Solidum conserva la salvaguarda $\max(\tilde d_n, d(\tilde\kappa))$, que evita la «curación» aparente entre pasos. **Usuario: «si» a conservar la versión de Solidum y documentar las diferencias.** Para el futuro control adaptativo del paso: estimador de la ec. 78, $\Delta t_{n+1}^2 = \xi\,\alpha_\text{ref}\,\Delta t_n^2\,\min_x 1/|\Delta\lambda_n(x) - \Delta\lambda_{n-1}(x)|$.
- **2026-11-11** · **Longitud característica de los cuadráticos corregida** (usuario, «si», tras el diagnóstico de la viga SEN, `informes/` no versionado). El diagnóstico (Tri6, Tri3 no estructurada y un prototipo de daño sólo de tracción) mostró que la banda vertical de las mallas finas la causa el daño isótropo, que degrada también la biela comprimida en la zona de cortante, y no el elemento ni la orientación de la malla; de paso halló que `Tri6` con el $\tfrac12\sqrt{A_e}$ de Abaqus disipaba $2.68\,G_fA$ en la barra. La prueba de objetividad ampliada a los seis cuadráticos dio con ½: $1.97\,G_fA$ (`Quad8`, `Quad9`, `Hex20`, `Hex27`), $2.68$ (`Tri6`) y $3.23$ (`Tet10`), en todas las mallas. La regla Q2 del 5 % da a cada cuadrático el factor del elemento lineal de su celda (1, $\sqrt2$ o $6^{1/3}$) y la barra disipa $0.997$–$0.999\,G_fA$ con todos (§4). Sin cambio en el material. Siguiente: spec del modelo 2 (daño de tracción y compresión) antes del bloque 3.
- **2026-11-12** · **Factor de localización de los cuadráticos** (usuario: «aplícalo»; Tri6 sin factor: «SI»). La viga de Petersson con `Quad8` y `Tri6` (`informes/reportes/reporte_viga_petersson_segundo_orden_2026-11-12.pdf`, no versionado) fue objetiva entre mallas pero disipó sólo $0.69$–$0.80\,G_fA$ en $\delta = 8$ mm (`Quad4` 0.94). El usuario lo atribuyó a la regularización, no al elemento: la formulación en desplazamientos gana flexibilidad al subir el orden, y la rigidez del primer orden es la que obliga a su banda a ablandarse entera. El diagnóstico por punto de Gauss lo confirmó (`Quad8`: dos de tres filas de Gauss en todos los elementos de la banda, $13/18$; `Tri6`: $2/3$ o 1) y coincide con Jirásek y Bauer (2012, §4) y Slobbe, Hendriks y Rots (2013) (búsqueda del bibliotecario en `informes/bibliotecario/busqueda_banda_cuadraticos_2026-11-12.md`; PDFs en `docs/referencias/`). Aplicado $h = \alpha\,c\,\Omega_e^{1/d}$ (§4) con $\alpha = 13/18$ en `Quad8`/`Quad9` 3×3 y ½ en 2×2. En `Tri6` el $2/3$ de la literatura sobrecorregía en la viga entallada (1.17–1.28 frente a 0.79–0.87 sin él): el usuario eligió dejarlo sin factor, del lado seguro. Limitaciones documentadas: $\alpha$ sólo es fiable con la banda alineada con la malla, y con esfuerzo uniforme en el elemento sobrestima la energía ($1/\alpha$). Siguiente: la $h$ por proyección del segundo método de Oliver (1989), $h = A_e/l_d$, para la orientación en todos los elementos (bloque 3).
- **2026-11-17** · **El aviso de los cuadráticos ya no recomienda elementos lineales** (usuario, 2026-11-12: el aviso se había escrito sin consultarlo y no estaba de acuerdo; lo que Solidum recomienda a quien lo usa es decisión del usuario, Reglas §3). El aviso, el manual de usuario, el catálogo, la spec 3D y el ADR 0028 describen la limitación sin recomendar un orden: con `Quad8`/`Quad9` la corrección es fiable con la banda alineada con la malla, con la banda inclinada da una respuesta demasiado dúctil y con esfuerzo uniforme disipa $G_f/\alpha$; `Tri6` y los cuadráticos 3D, sin factor, disipan de menos en flexión. Prueba: `test_aviso_describe_la_limitacion_sin_recomendar_un_orden` (falla con el aviso anterior). Sin cambio en la física ni en los resultados.
- **2026-11-17** · **Bloque 3: ancho de banda según la orientación de la grieta (propuesta de la IA, en revisión; §4, «Ancho según la orientación»).** El plan aprobado el 2026-11-11 era el segundo método de Oliver (1989), $h = A_e/l_d$. Al preparar la spec, la IA midió con un prototipo (barra con la malla girada) los dos candidatos que recomiendan Jirásek y Bauer (2012), que sólo los compararon a 0° y 45°, donde coinciden en cuadrados. En ángulos intermedios, con la banda y la dirección impuestas, la proyección disipa $G_f$ a ±4 % y el segundo método de Oliver se aparta hasta +19 % (`Quad4`) y −17 % (`Tri3`) (tabla de §4). La revisión del agente teórico (`informes/teorico/verificacion_teorica_2026-11-17_banda_orientacion.md`) mostró que esa medición es casi una identidad geométrica, que con la banda libre la proyección disipa de menos (0.87–0.91) y que es más sensible que O2 al error de dirección (−14 % frente a −1.5 % con 10° de error cerca de la alineación); corregidas la §4 y las preguntas. Preguntas, en este orden (una por vez):
  - **Q6. Método.** P (exacta en la geometría) u O2 (aprobada el 2026-11-11; menos sensible al error de dirección). El usuario aprobó (2026-11-17) un prototipo de la viga entallada con grieta libre y dirección estimada, en primer orden y, aparte, segundo orden: P tuvo el menor error medio (11 % frente a 16 % del actual y 15 % de O2) y la menor dispersión, pero ninguno quita la dependencia de la malla, y P quitó 10 % de energía con una malla casi regular y la banda fuera del eje. Literatura: P es la opción «Govindjee» de DIANA; O2 no está en ningún código (criterio de admisión: P lo decide el usuario; O2 iría a un módulo de usuario). **Recomendación: P en el programa principal y descartar O2**; el valor por omisión, en Q8.
  - **Q7. Dirección y momento.** **Recomendación**: deformación principal mayor en el centro del elemento, en el paso en que empieza su daño, fija después; un $h$ por elemento.
  - **Q8. Por omisión.** **Recomendación**: `crack_band_width = "projection"` por omisión en todos los sólidos; `"element"` disponible.
  - **Q9. Cuadráticos.** **Recomendación**: $\alpha\,h_P$ con los vértices; el aviso sigue igual.
  - **Q10. Límite $h_{\max}$.** **Recomendación**: comprobarlo al construir con la mayor distancia entre vértices, para que nunca falle a mitad del análisis.
  - **Q11. Aceptación.** **Recomendación**: la barra girada ±5 % (`Quad4`, `Tri3`, `Hex8`, `Tet4`), $h_P$ y $\mathbf n_c$ contra fórmulas cerradas, y la validación con la viga entallada en malla girada y no estructurada.
- **2026-11-17** · **Q6 en espera; línea de investigación abierta** (usuario: «hasta el momento sólo hemos probado lo que proponen otros autores… ¿qué propones?»). Propuesta de la IA, fuera del programa principal (no estándar, ADR 0020): **banda de grieta con trayectoria seguida y cuerda exacta**. El seguimiento global del módulo `discontinuities` traza la grieta; sólo se ablandan los elementos que cruza; cada uno usa $h_e = A_e/\ell_e$, con $\ell_e$ la longitud del segmento real de la grieta dentro del elemento, de modo que $D = G_f\sum_e \ell_e = G_f\,L$ por construcción (generaliza el $l_d$ de la tesis, la cuerda por el centro, a la cuerda de la grieta real). Motivo: en la viga entallada el error restante viene de dónde se forma la banda (desvío, ensanchamiento, vecinos parcialmente dañados, dirección), que ninguna estimación local de $h$ corrige. Primer paso aprobado: búsqueda de novedad (`bibliotecario-solidum`); si es nueva, spec propia en un módulo de usuario.
- **2026-11-17** · **Prototipo de la banda con trayectoria prescrita** (`informes/prototipos_2026-11-17/banda_orientacion/viga_trayectoria.py`; el usuario: «hazlo»). Viga entallada, 10 mallas (Quad4 y Tri3; alineadas y no estructuradas), grieta prescrita vertical desde la entalla, daño sólo en los elementos que cruza; cuatro anchos: actual, P, O2 y cuerda exacta $h = \sum A/\sum\ell$ (elementos de cuerda muy corta agrupados con un vecino de banda para que la banda quede conectada por lados, el caso «crossover» de Burnett y Schreyer). Error medio / máximo de $D/(G_f A)$: actual 14 % / 36 %, P 11 % / 33 %, O2 11 % / 32 %, **cuerda 9 % / 23 %**. La cuerda es la mejor, pero lejos del «exacto por construcción», porque fallan sus dos premisas: (a) no todos los elementos cruzados se rompen (la grieta encuentra cómo separarse sin alguno: disipa de menos, Quad4 `Ua` 0.84); (b) en los triángulos de una banda en zigzag la grieta abre con deformación mixta y el daño isótropo disipa también la energía de cortante (Tri3 `Ua` 1.15 con todos los cruzados rotos; no es error de IMPL-EX: con 4000 pasos, 1.145). **Lectura: en mallas no estructuradas el error restante es cinemático y constitutivo, no de la fórmula de $h$.**
- **2026-11-17** · **Banda restringida con cuerda exacta y daño sólo de tracción** (usuario: «si»; `informes/prototipos_2026-11-17/banda_orientacion/viga_trayectoria_m2.py`, `resultados_viga_trayectoria_m2.csv`). Mismas 10 mallas, trayectoria y anchos que el prototipo anterior, con dos materiales en la banda: el daño isótropo de Solidum y la partición de Faria, Oliver y Cervera (1998) con $d^- = 0$ (modelo 2, opción (a) de su Q1), mismo criterio de Rankine, misma ley y la misma banda de grieta. Nuevo: balance de energía por punto de integración ($D_g = \sum\tfrac12(\sigma_n+\sigma_{n+1})\Delta\varepsilon\,dV - \tfrac12\sigma\varepsilon\,dV$; su suma coincide con la $D$ global hasta el redondeo) frente a la energía nominal $G_f\,dV/h_g$, que separa el **déficit** de los puntos no rotos del **exceso** de los rotos. Error medio / máximo de $D/(G_f A)$ con la cuerda: isótropo 9.0 % / 22.7 %, **tracción 6.1 % / 14.9 %**. (1) **El exceso de los triángulos lo causa el daño isótropo** (medido): con daño isótropo, los Tri3 rotos disipan 1.07–1.62 veces su energía nominal (`U`, `Ua`, `U75`) y el indicador del estado de esfuerzos al dañarse, $r = \sum\psi_0\,\Delta d / \sum(\bar\sigma_1^2/2E)\,\Delta d$ (1 en tracción uniaxial, $2(1+\nu)$ en cortante puro con Rankine), lo explica (1.11–1.51); con daño sólo de tracción el exceso desaparece (rotos 0.92–1.03 de su nominal; Tri3 `U` 1.227 → 0.995). (2) **El déficit de los elementos cruzados que no se rompen no depende del material** (medido): Quad4 `Ua` 0.135 con los dos, Tri3 `Ua` 0.072, Quad4 `U75` 0.053, Tri3 `U5` 0.038 con tracción; es cinemático. (3) Con daño de tracción, P y O2 tienen menor error medio que la cuerda (4.5 y 4.4 %) porque en varias mallas asignan más energía nominal que $G_f L$ (hasta +24 % P y +36 % O2) y así compensan el déficit; la cuerda da el objetivo exacto y su error es sólo el déficit cinemático más 1–8 % de energía que los rotos aún no han disipado al final del ensayo (*hipótesis*: la carga se detiene en $\delta = 2$ mm con 1–3 % del pico). El daño de tracción necesita Newton completo (la partición espectral hace el paso no lineal aun con IMPL-EX); con la tangente congelada se detuvo en el paso 377.
- **2026-11-23** · **Bloque 3 cerrado sin programar; límites medidos del modelo** (usuario: «la 1» al elegir el siguiente frente, cerrar la banda y el modelo 2 tal como están con sus límites escritos, y «si» al plan). Q6–Q11 quedan sin decidir y el estudio de §4 se conserva como tal; el `h` de Solidum sigue siendo $\alpha\,c\,\Omega_e^{1/d}$. Razón (medida en §4 y en las entradas del 2026-11-17): con la grieta libre ningún ancho quita la dependencia de la malla (error medio 16 % actual, 11 % proyección, 15 % O2); con la trayectoria prescrita y la cuerda exacta queda un 9 % (6 % con daño sólo de tracción), y lo que queda es cinemático. Jirásek y Zimmermann (1998, p. 277; leído) atribuyen la transferencia espuria a que la interpolación continua no reproduce el salto de una grieta que se abre, y dan como remedio de fondo la discontinuidad dentro del elemento. Se cierra también la línea de la banda con trayectoria seguida y cuerda exacta (módulo de usuario, 2026-11-17). Del contrato se quitan el parámetro `crack_band_width` y la aceptación `rotated_mesh_bar`, que eran propuesta. En §10 se añaden los límites medidos: sesgo de malla con cifras, la viga SEN (biela comprimida degradada, no objetiva) y el panel en L (ablandamiento a 26 % entre mallas). Sin cambio de código ni de resultados.
