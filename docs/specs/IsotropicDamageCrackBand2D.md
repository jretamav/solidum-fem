# ISOTROPIC DAMAGE CRACK BAND 2D — daño escalar regularizado por banda de grieta, con integración implícita o IMPL-EX

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> **Variante de componente existente** (Reglas §4): subclase de `IsotropicDamage2D` ([spec](IsotropicDamage2D.md), validada), que no se modifica. Esta spec documenta **sólo lo que cambia**: (i) la velocidad de ablandamiento $\alpha$ deja de ser un dato y se calcula en cada punto de integración a partir de la energía de fractura $G_f$ y de la longitud característica $h$ del elemento (banda de grieta, Bažant y Oh 1983), para que la energía disipada por unidad de área de grieta sea $G_f$ con cualquier malla; (ii) una medida de deformación equivalente más, la de **Rankine**; (iii) una segunda forma de integrar el daño, **IMPL-EX** (Oliver, Huespe y Cante 2008); (iv) un tope de daño propio. Todo lo demás (esfuerzo nominal y efectivo, estado 3D asociado, variable histórica, Kuhn-Tucker, ley exponencial, medidas existentes, tangente consistente) es el del padre.
>
> Uso previsto: agrietamiento de estructuras de concreto con fisuración distribuida, en 2D y 3D (la variante 3D está en [IsotropicDamageCrackBand3D](IsotropicDamageCrackBand3D.md)).
>
> **Estado: aprobada por el usuario (2026-11-11, «si», con las cinco recomendaciones Q1–Q5) y validada.** Bloque 1 (banda de grieta, medida de Rankine, tope propio y siembra de $h$ en los sólidos 2D y 3D) y bloque 2 (integración IMPL-EX, §8, con la pieza genérica del pseudo-tiempo, [ADR 0034](../adr/0034-pseudo-tiempo-de-los-solvers-incrementales.md), y la prueba estructural de Van Vliet) implementados el mismo día.

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
- Sesgo de malla: la banda tiende a seguir las líneas de la malla; con grietas oblicuas a la malla puede ensancharse (dos filas de elementos) y disipar más que $G_f$.
- Transferencia espuria de esfuerzos: los elementos vecinos a la banda pueden seguir transmitiendo esfuerzo a través de ella en estados mixtos.
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
  - "longitud proyectada de Oliver (1989) y cualquier h que dependa de la dirección de la grieta"
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
  - "Oliver J. (1989). A consistent characteristic length for smeared cracking models. IJNME 28, 461–474."
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
- **2026-11-12** · **Factor de localización de los cuadráticos** (usuario: «aplícalo»; Tri6 sin factor: «SI»). La viga de Petersson con `Quad8` y `Tri6` (`informes/reporte_viga_petersson_segundo_orden_2026-11-12.pdf`, no versionado) fue objetiva entre mallas pero disipó sólo $0.69$–$0.80\,G_fA$ en $\delta = 8$ mm (`Quad4` 0.94). El usuario lo atribuyó a la regularización, no al elemento: la formulación en desplazamientos gana flexibilidad al subir el orden, y la rigidez del primer orden es la que obliga a su banda a ablandarse entera. El diagnóstico por punto de Gauss lo confirmó (`Quad8`: dos de tres filas de Gauss en todos los elementos de la banda, $13/18$; `Tri6`: $2/3$ o 1) y coincide con Jirásek y Bauer (2012, §4) y Slobbe, Hendriks y Rots (2013) (búsqueda del bibliotecario en `informes/busqueda_banda_cuadraticos_2026-11-12.md`; PDFs en `docs/referencias/`). Aplicado $h = \alpha\,c\,\Omega_e^{1/d}$ (§4) con $\alpha = 13/18$ en `Quad8`/`Quad9` 3×3 y ½ en 2×2. En `Tri6` el $2/3$ de la literatura sobrecorregía en la viga entallada (1.17–1.28 frente a 0.79–0.87 sin él): el usuario eligió dejarlo sin factor, del lado seguro. Limitaciones documentadas: $\alpha$ sólo es fiable con la banda alineada con la malla, y con esfuerzo uniforme en el elemento sobrestima la energía ($1/\alpha$). Siguiente: la $h$ por proyección del segundo método de Oliver (1989), $h = A_e/l_d$, para la orientación en todos los elementos (bloque 3).
- **2026-11-17** · **El aviso de los cuadráticos ya no recomienda elementos lineales** (usuario, 2026-11-12: el aviso se había escrito sin consultarlo y no estaba de acuerdo; lo que Solidum recomienda a quien lo usa es decisión del usuario, Reglas §3). El aviso, el manual de usuario, el catálogo, la spec 3D y el ADR 0028 describen la limitación sin recomendar un orden: con `Quad8`/`Quad9` la corrección es fiable con la banda alineada con la malla, con la banda inclinada da una respuesta demasiado dúctil y con esfuerzo uniforme disipa $G_f/\alpha$; `Tri6` y los cuadráticos 3D, sin factor, disipan de menos en flexión. Prueba: `test_aviso_describe_la_limitacion_sin_recomendar_un_orden` (falla con el aviso anterior). Sin cambio en la física ni en los resultados.
