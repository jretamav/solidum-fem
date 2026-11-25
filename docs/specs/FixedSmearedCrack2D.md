# FIXED SMEARED CRACK 2D — grieta distribuida fija con retención de cortante y banda de grieta

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> **Componente nuevo** (modelo 3 de la familia de modelos de concreto, [ROADMAP](../ROADMAP.md)): grieta distribuida **fija** con descomposición de la deformación en la del concreto entre grietas y la de la grieta (Rots 1988; de Borst y Nauta 1985; de Borst, Crisfield, Remmers y Verhoosel 2012, §6.8.1). La grieta recuerda su dirección y transmite cortante según una ley de retención. Reutiliza piezas de los modelos 1 y 2: longitud característica $h$ sembrada por el elemento ([IsotropicDamageCrackBand2D](IsotropicDamageCrackBand2D.md) §4), ley de Hordijk en la apertura (§2b de esa spec) e integración IMPL-EX ([ADR 0034](../adr/0034-pseudo-tiempo-de-los-solvers-incrementales.md)).
>
> **Objetivo** (usuario, 2026-11-24, «si, vamos con la opción 1»): el modelo estándar de grieta distribuida que ofrecen DIANA, Abaqus y ANSYS, para estructuras de concreto (sobre todo reforzado, con grietas repartidas) y para comparar con esos programas y con la discontinuidad embebida, como hizo Rots (1988). **No** es un remedio para la bajada de la viga SEN: en la viga de Arrea–Ingraffea todas las grietas distribuidas quedan demasiado rígidas tras el pico, y cuando alguna coincide con el ensayo es por bloqueo de esfuerzos, no por trabazón de agregados (Rots 1988, pp. 85–89, *leído*; `informes/bibliotecario/busqueda_modelo3_2026-11-24.md`). La SEN entra sólo como observación. Criterio de parada: objetividad entre mallas en la escalera de validación, como el modelo 2.
>
> **Estado: borrador (2026-11-24).** Q1 (sólo grieta fija), Q2 (retención ligada a la apertura, sin parámetro) y Q3 (aplastamiento opcional con las piezas del modelo 2) decididas. Q4 (aceptación) aprobada. Revisada por `teorico-solidum` el 2026-11-24, con sus correcciones aplicadas y H1 (límite de $h$) decidida: **lista para programar el bloque 1** (tracción y cortante). **Bloque 1 programado el 2026-11-24** (sección «Implementación» y Diálogo). **Bloque 2 (aplastamiento) redactado el 2026-11-24** (H4: §4, §8b, §9, §10 y contrato) y revisado por `teorico-solidum` («sí, se puede programar», correcciones aplicadas); Q5 decidida. **Bloque 2 programado y verificado el 2026-11-24** («Implementación» y Diálogo). **Validación del 2026-11-25: no validado; limitación documentada por decisión del usuario** (§10, primer caveat, y Diálogo): en vigas con entalla en modo I las grietas fijas nacen inclinadas y el modelo sobrestima el pico (11–24 %) y la rama de ablandamiento. `status: implemented`; la escalera no sigue (panel en L y SEN sin correr).

---

## Especificación física

### 0. Descripción general

Grieta distribuida en pequeñas deformaciones para sólidos 2D (esfuerzo plano o deformación plana). El material es elástico isótropo hasta que el esfuerzo principal mayor alcanza $f_t$; en ese momento se forma una grieta cuya **dirección queda fija** (normal $\mathbf n$ = dirección principal mayor del esfuerzo en el instante del agrietamiento). Puede formarse después una **segunda grieta perpendicular** a la primera (dirección $\mathbf t$), como en DIANA, Abaqus y ANSYS (Q1).

`STRAIN_DIM = 3`, Voigt del proyecto $[\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]$ con $\gamma_{xy} = 2\varepsilon_{xy}$. Convención: $\sigma > 0$ ⇔ tracción.

### 1. Descomposición de la deformación (formulación elegida por la IA, 2026-11-24)

$$\boldsymbol\varepsilon = \boldsymbol\varepsilon^{co} + \boldsymbol\varepsilon^{cr}, \qquad \boldsymbol\varepsilon^{co} = \mathbf C_e^{-1}:\boldsymbol\sigma, \qquad \boldsymbol\varepsilon^{cr} = \mathbf N\,\mathbf e^{cr},$$

con $\mathbf C_e$ la elasticidad isótropa ($E$, $\nu$) del concreto entre grietas, $\mathbf e^{cr} = [e_{nn}, e_{tt}, \gamma_{nt}]$ las deformaciones de la grieta en sus ejes locales (apertura normal a $\mathbf n$, apertura normal a $\mathbf t$ de la segunda grieta y deslizamiento, con $\gamma_{nt}$ de ingeniería) y $\mathbf N$ la matriz de transformación de los ejes de la grieta a los globales (de Borst et al. 2012, ecs. 6.127–6.135). Con $\mathbf n = (\cos\theta, \sin\theta)$, $\mathbf t = (-\sin\theta, \cos\theta)$, $c = \cos\theta$, $s = \sin\theta$ y el Voigt del proyecto:

$$\mathbf N = \begin{bmatrix} c^2 & s^2 & -cs \\ s^2 & c^2 & cs \\ 2cs & -2cs & c^2 - s^2 \end{bmatrix},$$

de modo que $\mathbf N\,\mathbf e^{cr}$ es el Voigt de $e_{nn}\,\mathbf n\otimes\mathbf n + e_{tt}\,\mathbf t\otimes\mathbf t + \tfrac12\gamma_{nt}(\mathbf n\otimes\mathbf t + \mathbf t\otimes\mathbf n)$, la tracción en la grieta es $\mathbf t^{cr} = \mathbf N^\top\boldsymbol\sigma = [\sigma_{nn}, \sigma_{tt}, \tau_{nt}]$ y $\boldsymbol\sigma:\boldsymbol\varepsilon^{cr} = \sigma_{nn}e_{nn} + \sigma_{tt}e_{tt} + \tau_{nt}\gamma_{nt}$ (calculado por el teórico). Se escribe aquí porque el libro mide el ángulo de dos maneras (p. 202, del eje $x$ al eje $s$; p. 205, de $x$ a la normal) y copiarlo puede cambiar el signo de la columna de cortante. Las dos grietas comparten $\gamma_{nt}$ (el tensor $\mathbf n\otimes\mathbf t + \mathbf t\otimes\mathbf n$ es el mismo): tres columnas bastan. $\mathbf t^{cr}$ se relaciona con $\mathbf e^{cr}$ por la ley de la grieta (§2–§3).

**Ángulo de la grieta:** $\theta = \tfrac12\operatorname{atan2}(2\sigma_{xy},\ \sigma_{xx} - \sigma_{yy}) \in (-\pi/2, \pi/2]$, la dirección principal mayor del esfuerzo de prueba en el paso de inicio (§8), como en el modelo 2.

**Banda de grieta.** La grieta se reparte en el ancho $h$ del elemento: apertura $w = h\,e_{nn}$ y deslizamiento $s = h\,\gamma_{nt}$. La ley de la grieta se escribe en $w$ y $s$, de modo que es una ley cohesiva repartida en la banda.

**Por qué esta formulación y no la de deformación total de DIANA** (decisión de software de la IA; el usuario retiró el 2026-11-24 la restricción sobre modelos por reglas en sólidos y pidió «como sea mejor para Solidum y el modelado del concreto»):
- El concreto entre grietas conserva su Poisson: abrir una grieta no contrae lateralmente el material vecino (Rots 1988 conserva el Poisson completo, nota de la p. 38, *leído* por el bibliotecario). **Corrección del 2026-11-25:** en esfuerzo plano esto no distingue las dos formulaciones. La de deformación total con el Poisson reducido con el daño, $\nu_{nt} = \nu_0 E_n/E_0$ (DIANA §46.5.4.1; Slobbe 2015, ecs. 2.9 y 2.12), tiene la flexibilidad $\bigl[\begin{smallmatrix}1/E_n & -\nu_0/E_0\\ -\nu_0/E_0 & 1/E_t\end{smallmatrix}\bigr]$, que es la de esta spec con $1/E_n = 1/E_0 + c_n/h$: el «Poisson reducido» es el Poisson entero del concreto entre grietas escrito de otra forma (*calculado* por el teórico con SymPy, `informes/teorico/verificacion_teorica_2026-11-25_modelo3_slobbe.md` §2). La elección de esta forma se sostiene por los puntos siguientes, no por éste. Deformación plana: no comprobado.
- La energía disipada es exactamente $G_f/h$ por unidad de volumen en tracción uniaxial, también en deformación plana (en el modelo 1 queda el factor $1-\nu^2$, documentado sin corregir).
- La ley de la grieta en $w$ es la ley cohesiva del material: la de Hordijk entra tal cual.
- Es la base natural del modelo 4 (fractura con plasticidad, $\boldsymbol\varepsilon = \boldsymbol\varepsilon^e + \boldsymbol\varepsilon^p + \boldsymbol\varepsilon^{cr}$, como ATENA).
- Costo: un Newton local en el punto material (dos incógnitas en 2D, las aperturas; el cortante sale en forma cerrada, §8). Con IMPL-EX la ley de la grieta queda secante y conocida dentro del paso, y el esfuerzo sale en forma cerrada. Con aplastamiento, hasta tres incógnitas acopladas y sin forma cerrada en IMPL-EX (§8b, §9).

**Deformación plana.** La grieta sólo está en el plano; $\varepsilon_{zz} = 0$ y $\sigma_{zz} = \nu(\sigma_{xx} + \sigma_{yy})$ sale de la elasticidad del concreto entre grietas. Una grieta normal a $z$ queda fuera de alcance.

### 2. Ley normal de la grieta (apertura)

- **Inicio:** la grieta se forma cuando el esfuerzo principal mayor alcanza $f_t$ (criterio de Rankine sobre el esfuerzo, que antes de agrietarse es el del concreto entre grietas: el elástico, o el de §4 si hay aplastamiento). La segunda grieta, cuando $\sigma_{tt}$ alcanza $f_t$. El aplastamiento no cambia el inicio de la primera grieta: $\boldsymbol\sigma = \bar{\boldsymbol\sigma} - d^-\bar{\boldsymbol\sigma}^-$ tiene las direcciones principales de $\bar{\boldsymbol\sigma}$ y los esfuerzos principales $\bar\sigma_i$ (si $\bar\sigma_i \ge 0$) o $(1-d^-)\,\bar\sigma_i$ (si $\bar\sigma_i < 0$), así que su principal mayor positivo es $\bar\sigma_I$ y $\theta$ es la misma (deducción). Sí puede adelantar la segunda: $\sigma_{tt} = \bar\sigma_{tt} - d^-\bar\sigma^-_{tt} \ge \bar\sigma_{tt}$ cuando los ejes de la grieta no son principales.
- **Ablandamiento:** $\sigma_{nn} = f(w)$ con la ley elegida en `softening`: exponencial en la apertura, $f = f_t\,e^{-f_t w/G_f}$ (la «exponential» de DIANA), o la de Hordijk, $f = f_t\,g(w/w_c)$ con $w_c = 5.136\,G_f/f_t$ (modelos 1 y 2, §2b). Por omisión la exponencial, como decidió el usuario para los modelos 1 y 2 (Q13 del 2026-11-24). Nota: la «exponencial» de este modelo es exponencial en $w$; la del modelo 1 es exponencial en la deformación total con $\alpha(h)$. Encierran la misma $G_f$ y difieren en la forma, más cuanto mayor es $h$ (calculado por el teórico; diferencia máxima de $\sigma(w)$ y razón de pendientes iniciales modelo 1 / modelo 3): $h$ = 5 mm, 0.13 % de $f_t$ y 0.993; 25 mm, 0.65 % y 0.964; 100 mm, 2.6 % y 0.870; 300 mm, 7.4 % y 0.690 ($E$ = 30 GPa, $f_t$ = 3 MPa, $G_f$ = 100 N/m). Coinciden cuando $h \to 0$.
- **Descarga y recarga:** secante hacia el origen de la grieta, $\sigma_{nn} = f(w_\text{max})\,w/w_\text{max}$, con $w_\text{max}$ la apertura máxima histórica.
- **Cierre:** con $w = 0$ la grieta está cerrada y transmite toda la compresión normal ($e_{nn} \ge 0$).
- **Límite del tamaño del elemento** (H1, decidida por el usuario el 2026-11-24): $h < E/((1+\nu)\lvert f'(0)\rvert)$, es decir $h < \ell_{ch}/(1+\nu)$ con la exponencial en $w$ ($\ell_{ch} = EG_f/f_t^2$) y $h < 0.738\,\ell_{ch}/(1+\nu)$ con la de Hordijk; por encima, `ValueError` al sembrar $h$, con mensaje claro. Es la condición necesaria y suficiente para que el problema local tenga una sola solución también con las dos direcciones en carga (§8), en esfuerzo plano y en deformación plana. Con una sola dirección bastaría $h < E/\lvert f'(0)\rvert$ (el límite de §2b del modelo 1); con concreto usual ($E$ = 30 GPa, $f_t$ = 3 MPa, $G_f$ = 100 N/m, $\nu$ = 0.2) el límite queda en 278 mm (exponencial) y 205 mm (Hordijk), frente a 333 y 246 mm con una dirección. El modelo 1 con su exponencial admite $h < 2\,\ell_{ch}$.

### 3. Ley de cortante de la grieta (Q2, decidida por el usuario el 2026-11-24: ligada a la apertura, sin parámetro)

La grieta resiste el deslizamiento con la mitad de su rigidez normal secante:

$$\tau_{nt} = k_s\,s, \qquad k_s = \frac{k_n}{2} = \frac{f(w_\text{max})}{2\,w_\text{max}}, \qquad k_n = \frac{f(w_\text{max})}{w_\text{max}},$$

con $s = h\,\gamma_{nt}$ el deslizamiento y $f$ la ley de ablandamiento de §2. En deformaciones de la grieta, $\tau_{nt} = D_{II}\,\gamma_{nt}$ con $D_{II} = h\,k_s$.

- **Equivalencia con la opción «por daño» de DIANA, con una grieta y tracción normal uniaxial** (calculado, y comprobado por el teórico). Si el módulo de cortante reducido es $G_\text{red} = E_\text{red}/(2(1+\nu_\text{red}))$ con $\nu_\text{red} = \nu E_\text{red}/E$ y $E_\text{red} = \sigma_{nn}/\varepsilon_{nn}$ el módulo secante de la banda, $\varepsilon_{nn} = \sigma_{nn}/E + w/h$, entonces $1/G_\text{red} - 1/G = 2/E_\text{red} - 2/E = 2w/(h\,\sigma_{nn})$, y el deslizamiento de la grieta, $s = h\,\tau\,(1/G_\text{red} - 1/G)$, da exactamente $\tau/s = \sigma_{nn}/(2w)$. Precisiones: el manual de DIANA 10.4 dice que se reducen rigideces y Poisson y que $G = E/(2(1+\nu))$ (ec. 46.103), pero no muestra una ecuación propia de «damage based»; Slobbe (2015, ec. 2.14), con el mismo modelo, da $1/G_\text{red} - 1/G_0 = 2/E_\text{min} - 2/E_0$, que con una grieta es exactamente esta ley (*calculado* por el teórico, 2026-11-25). **Corrección del 2026-11-25:** con $\sigma_{tt} \ne 0$ la coincidencia se mantiene (aquí se había supuesto lo contrario como *hipótesis*). Con dos grietas, la suma en serie es elección de esta spec, coherente con de Borst et al. ec. 6.134; Slobbe usa $E_\text{min}$, la más dañada, y queda más rígido a cortante (pico del punto material 1.324 frente a 1.265·$f_t$ con la grieta a 30°, *calculado*).
- **Objetividad** (calculado): la ley $\tau$–$s$ no contiene $h$; la rigidez al deslizamiento de la grieta es la misma con cualquier malla. Con una fracción fija $\beta$ del cortante elástico (descartada) sería $\tau/s = \beta G/((1-\beta)\,h)$: el doble con elementos de la mitad de tamaño.
- **Evolución:** el cortante transmitido cae a cero conforme la grieta se abre ($k_s \propto f(w)/w$). Antes de agrietarse no hay deslizamiento; se implementa con la flexibilidad $c_s = 1/k_s = 2\,w_\text{max}/f(w_\text{max})$, que vale 0 en $w_\text{max} = 0$ (sin infinitos).
- **Descarga y cierre:** $k_s$ depende sólo de la historia ($w_\text{max}$); al cerrarse la grieta conserva ese valor, sin recuperar rigidez a cortante (como el modelo 2, que tampoco la recupera; Cervera y Tesei 2017 fuera de alcance).
- **Dos grietas:** cada una desliza con su propia ley; las flexibilidades se suman (en serie).
- Sin tope de resistencia a cortante (el de la MCFT que usa ATENA queda fuera de alcance).

Referencia de la familia: de Borst et al. (2012), ec. 6.117 ($\beta = 1-\omega$, retención ligada al daño); DIANA «damage based»; ATENA usa la misma forma con un factor ajustable, $\tau/s = s_F\,\sigma_{nn}(w)/w$ con $s_F = 20$ por omisión (ec. 2.47; aquí $s_F = \tfrac12$ fijo).

### 4. Compresión (Q3, decidida por el usuario el 2026-11-24: aplastamiento opcional con las piezas del modelo 2)

Opcional: con $f_c$ y $G_c$ dados se activa; sin ellos el concreto entre grietas es elástico en compresión. El daño de compresión actúa sobre el **concreto entre grietas**; la tracción ya la llevan las grietas:

$$\boldsymbol\sigma = \bar{\boldsymbol\sigma} - d^-\,\bar{\boldsymbol\sigma}^-, \qquad \bar{\boldsymbol\sigma} = \mathbf C_e:\boldsymbol\varepsilon^{co},$$

con $\bar{\boldsymbol\sigma}^-$ la parte negativa (espectral) del esfuerzo efectivo del estado 3D asociado (esfuerzo plano: $\bar\sigma_{zz} = 0$; deformación plana: $\bar\sigma_{zz} = \lambda\,\operatorname{tr}\boldsymbol\varepsilon^{co}$), y $d^-$ el del modelo 2 ([TensionCompressionDamage2D](TensionCompressionDamage2D.md) §3): medida de Drucker-Prager sobre $\bar{\boldsymbol\sigma}^-$ con $\beta = f_{bc}/f_c$ (1.16 por omisión), parábola de Feenstra con la energía de aplastamiento $G_c/h$, sin el límite de DIANA, y el mismo tope. Es el modelo 2 con $d^+ = 0$ aplicado al concreto entre grietas.

**Evolución de $d^-$** (bloque 2): $\kappa^-_{n+1} = \max\bigl(\kappa^-_n,\ \varepsilon^-_\text{eq}(\bar{\boldsymbol\sigma}^-)\bigr)$ con $\bar{\boldsymbol\sigma} = \mathbf C_e:\boldsymbol\varepsilon^{co}$ y $\boldsymbol\varepsilon^{co} = \boldsymbol\varepsilon - \mathbf N\mathbf e$. Como el modelo 2 es explícito en la deformación, $d^-$ es una función explícita de las deformaciones de la grieta $\mathbf e$: no es una incógnita más del problema local (§8b), pero lo vuelve no lineal y lo acopla. $d^- = 0$ mientras $\kappa^- \le f_c/(3E)$ (tramo lineal de la parábola de Feenstra). En deformación plana $\boldsymbol\varepsilon^{co}_{zz} = 0$ (la grieta está en el plano), así que $\bar\sigma_{zz} = \lambda\,(\varepsilon^{co}_{xx} + \varepsilon^{co}_{yy})$ y $\sigma_{zz} = \bar\sigma_{zz} - d^-\langle\bar\sigma_{zz}\rangle_-$, como en el modelo 2.

- **Grieta cerrada:** la compresión normal pasa por el concreto entre grietas y puede aplastarlo.
- **Compresión paralela a una grieta abierta: se conserva completa.** Abrir la grieta no contrae el concreto vecino (§1). $\sigma_{tt}$ sólo depende de la apertura a través del Poisson con la tracción que todavía transmite la grieta: $\sigma_{tt} = E\,\varepsilon_t + \nu\,\sigma_{nn}$ en esfuerzo plano y $\sigma_{tt} = E\,\varepsilon_t/(1-\nu^2) + \nu\,\sigma_{nn}/(1-\nu)$ en deformación plana (calculado por el teórico; con $\varepsilon_t = -10^{-3}$, de −29.4 a −30.0 MPa entre $w = 0$ y $w = w_c$), y tiende a la compresión elástica completa cuando la grieta se abre del todo. Es la diferencia de fondo con el modelo 2, que la pierde por Poisson a partir de un umbral que depende de $h$ (caveat de su §1).
- **Aplastamiento paralelo a una grieta abierta** (deducción, esfuerzo plano, grieta en $\mathbf n$ abierta y ejes de la grieta principales): $\bar\sigma_{nn} = \sigma_{nn} \ge 0$ no entra en $\bar{\boldsymbol\sigma}^-$, que es uniaxial, $\bar{\boldsymbol\sigma}^- = \bar\sigma_{tt}\,\mathbf t\otimes\mathbf t$ con $\bar\sigma_{tt} = E\,\varepsilon_t + \nu\,\sigma_{nn}$. En carga monótona, $-\sigma_{tt} = \sigma_c(\kappa^-)$ con $\kappa^- = \lvert\bar\sigma_{tt}\rvert/E$: la curva uniaxial de Feenstra en $\varepsilon_t + \nu\sigma_{nn}/E$, con el pico en $f_c$ **con cualquier apertura** (sin reducción de Vecchio–Collins, Q3) y el área de la rama descendente $G_c/h$. En esos ejes $\sigma_{nn} = \bar\sigma_{nn}$ no depende de $d^-$: el problema normal de la grieta es el del bloque 1. Comprobado por el teórico con el problema local completo ($-\sigma_{tt} = \sigma_c$ a $1.2\cdot10^{-12} f_c$, pico $1.00000\,f_c$ con $w$ = 0, 0.1, 0.5 y 1·$w_c$). En deformación plana, $\bar\sigma_{zz} = \nu(\bar\sigma_{nn} + \bar\sigma_{tt}) < 0$ confina: $\bar{\boldsymbol\sigma}^- = \operatorname{diag}(0, \bar\sigma_{tt}, \bar\sigma_{zz})$ con $\bar\sigma_{zz}/\bar\sigma_{tt} = \rho = \nu\,(1 - \sigma_{nn}/\lvert\bar\sigma_{tt}\rvert)$ y el pico es $f_c/\varphi(\rho)$, con $\varphi = E\,\varepsilon^-_\text{eq}/\lvert\bar\sigma_{tt}\rvert$ del estado $(0, -1, -\rho)$ (calculado por el teórico, H7): 1.1397·$f_c$ sin tracción en la grieta ($\nu$ = 0.2, $\beta$ = 1.16; el 1.14·$f_c$ del modelo 2, caveat de su §10) y 1.1371, 1.1317 y 1.1127·$f_c$ con $\sigma_{nn}/\lvert\bar\sigma_{tt}\rvert$ = 0.02, 0.06 y 0.2.
- **Sin acoplamiento entre el aplastamiento y la grieta** (herencia del modelo 2, §5 de su spec): un concreto aplastado conserva su resistencia $f_t$ y su ley de ablandamiento en tracción; la grieta se forma y se abre con la misma ley con cualquier $d^-$.
- **Sin reducción por agrietamiento transversal** (Vecchio–Collins; DIANA ecs. 46.139–46.141, ATENA ec. 2.78). Descartada en Q3 porque su fórmula usa la deformación de tracción transversal, que en la banda es $w/h$: dependería de la malla. Se puede añadir después en una forma con $w$ y la separación entre grietas si una validación con paneles reforzados lo pide.
- Con $d^- > 0$ y $\nu \ne 0$ el concreto entre grietas hereda el caveat del modelo 2: su ley no deriva de un potencial (trabajo no nulo en ciclos cerrados con daño fijo).
- **Lo que el aplastamiento cambia en la formulación numérica** (H4 del teórico, `informes/teorico/verificacion_teorica_2026-11-24_modelo3_spec.md` §6, calculado por él en esfuerzo plano; tratamiento en §8b y §9):
  - $d^-$ evoluciona dentro del problema local, como función explícita de $\mathbf e$ (arriba).
  - El cortante deja de desacoplarse del problema normal cuando los ejes de la grieta no son principales de $\bar{\boldsymbol\sigma}$: hasta tres incógnitas y jacobiano no simétrico.
  - **El límite de $h$ de §2 deja de garantizar una sola solución.** La rigidez del concreto normal a una grieta que puede abrir baja hasta cerca de $(1-d^-)\,C_{11}$ cuando la normal queda próxima a la dirección aplastada (mínimo $K_{nn}/C_{11}$ con $\sigma_{nn} > 0$: 0.51, 0.10 y 0.011 con $d^-$ = 0.5, 0.9 y 0.99, en 20 000 estados al azar); con $d^- = 0.9$ haría falta $h < 0.1\,\ell_{ch}$ (cota pesimista), y ningún $h$ basta cuando $d^- \to 1$. Hay además un segundo mecanismo, la grieta y el aplastamiento ablandando en el mismo paso (§8b). No ocurre en el caso típico, la biela comprimida paralela a la grieta (viñeta anterior: $K_{nn} = C_{11}$); en caminos al azar el teórico no encontró ningún estado con varias soluciones, sólo en estados construidos con $d^- \ge 0.8$–0.99 (§8b). En esos estados Solidum elige la solución por continuidad, sin límite ni aviso nuevos (Q5, decidida por el usuario).
  - IMPL-EX pierde la forma cerrada: con $\tilde d^-$ congelado, $\bar{\boldsymbol\sigma}^-$ no es lineal en $\mathbf e$ (las direcciones principales de $\bar{\boldsymbol\sigma}$ giran con $\mathbf e$ cuando sus esfuerzos principales tienen signos distintos).
  - Se pierden el potencial χ de §6, el trabajo nulo en ciclos con la historia fija y la simetría de la secante; por encima de $d^* = 2\sqrt{1-\nu^2}/(1+\sqrt{1-\nu^2})$ (0.9898 con $\nu$ = 0.2 en esfuerzo plano; §9) la parte simétrica de la tangente del concreto a daño fijo deja de ser definida positiva (valor propio mínimo de su parte simétrica, $+0.208E$, $+0.042E$ y $-0.0002E$ con $d^-$ = 0.5, 0.9 y 0.99; herencia del modelo 2).

### 5. Alternativas estudiadas (Q1)

| Variante | Qué hace | Admisión (Reglas §1) | Lo que dice la literatura |
|---|---|---|---|
| **Fija (adoptada)** | Dirección congelada al agrietarse; segunda grieta perpendicular; cortante según la ley de retención | de Borst et al. 2012 §6.8.1 + DIANA, Abaqus, ANSYS | Con retención significativa, respuesta rígida y bloqueo por esfuerzo inclinado; con retención despreciable ablanda (Rots 1988, pp. 86–89, Tabla 6.1) |
| Rotante | La grieta gira con la deformación principal; sin cortante en la grieta | de Borst et al. 2012 ecs. 6.122–6.126 + sólo DIANA | Converge mejor, más sensible a la orientación de la malla (Rots pp. 81, 87–88); bloqueo en la viga de Arrea–Ingraffea (Jirásek y Zimmermann 1998); su β equivalente se hace negativo (Rots p. 40) |
| Fija y rotante a elegir | Las dos | — | — |

Razón de la decisión: la fija es lo nuevo respecto al modelo 2, que ya gira con las direcciones principales y no tiene cortante en la grieta; es la que traen por omisión los códigos generales y la que se usa en concreto reforzado.

---

## Formulación numérica

### 6. Las dos direcciones de grieta y su ley secante

Al formarse la primera grieta (§2) quedan fijas las dos direcciones $\mathbf n$ y $\mathbf t$; la segunda grieta es la dirección $\mathbf t$, que se comporta igual pero con $w_{\text{max},t} = 0$ hasta que $\sigma_{tt}$ alcanza $f_t$. Con la historia fija, cada dirección $i \in \{n, t\}$ tiene la flexibilidad (por unidad de longitud de banda)

$$c_i = \frac{w_{\text{max},i}}{f(w_{\text{max},i})} \quad (c_i = 0 \text{ si } w_{\text{max},i} = 0), \qquad c_s = 2\,(c_n + c_t)$$

(§3: cada grieta desliza con $k_s = k_n/2$ y las dos en serie). $c_i$ se calcula siempre con el tope ya aplicado (abajo): con Hordijk y $w \ge w_c$, $f = 0$ y $w/f$ dividiría entre cero; lo mismo con la $\tilde w_\text{max}$ extrapolada de IMPL-EX. La ley es unilateral en la apertura (grieta cerrada con $\sigma_{ii} \le 0$) y lineal en el cortante:

$$e_{ii} = \frac{c_i}{h}\,\langle\sigma_{ii}\rangle_+, \qquad \gamma_{nt} = \frac{c_s}{h}\,\tau_{nt}.$$

**Potencial a historia fija** (calculado). Con la compresión elástica la ley deriva de la energía complementaria

$$\chi(\boldsymbol\sigma) = \tfrac12\,\boldsymbol\sigma:\mathbf C_e^{-1}:\boldsymbol\sigma + \frac{1}{2h}\Bigl(c_n\langle\sigma_{nn}\rangle_+^2 + c_t\langle\sigma_{tt}\rangle_+^2 + c_s\,\tau_{nt}^2\Bigr), \qquad \boldsymbol\varepsilon = \partial\chi/\partial\boldsymbol\sigma,$$

estrictamente convexa, coerciva y de clase $C^1$: la relación $\boldsymbol\varepsilon \mapsto \boldsymbol\sigma$ es única, continua, lineal a trozos y positivamente homogénea; la secante es simétrica y definida positiva y el trabajo en un ciclo cerrado con la historia fija es nulo, **también con $\nu \ne 0$** (a diferencia del modelo 2). Comprobado por el teórico: de las cuatro combinaciones abierta/cerrada, exactamente una es consistente en 20 000 de 20 000 estados al azar; trabajo en ciclos con $\nu = 0.2$ ≤ 7.5·10⁻⁷ de la energía máxima (error de cuadratura), y 0.3–0.67 si se recupera el cortante al cerrar la grieta: **no recuperar la rigidez a cortante al cerrar es lo que conserva el potencial**. Con aplastamiento ($d^- > 0$) se pierde esta propiedad (§4).

**Tope de la flexibilidad** (`damage_max`, como los modelos 1 y 2). Daño uniaxial equivalente de la banda en la dirección $i$: $\omega_i = (E\,w_{\text{max},i}/h)/(f(w_{\text{max},i}) + E\,w_{\text{max},i}/h)$, que es el $d$ de los modelos 1 y 2 en tracción uniaxial. Cuando $\omega_i$ llega al tope, $c_i$ deja de crecer: $c_i \le c_\text{max} = h\,\omega_\text{max}/((1-\omega_\text{max})\,E)$, y la grieta conserva una rigidez residual. Equivale a sustituir la envolvente por $\max(f(w),\ w/c_\text{max})$, que sigue siendo convexa; $c_s \le 4\,c_\text{max}$. $\omega_n$ es la variable principal del post-proceso ($\omega_i$ se define con $E$ en las dos hipótesis; es el $d$ uniaxial de los modelos 1 y 2 en esfuerzo plano).

**Energía disipada hasta el tope** (referencia exacta de `tension_energy`, calculada por el teórico): $D_\text{tope} = \frac1h\Bigl[\int_0^{w_\text{tope}} f\,dw - \tfrac12 f(w_\text{tope})\,w_\text{tope}\Bigr]$, con $w_\text{tope}$ la raíz de $w = c_\text{max}\,f(w)$. Con `damage_max` = 1 − 10⁻⁶, $D_\text{tope}\,h/G_f$ = 0.99757 / 0.99934 / 0.99979 (exponencial, $h$ = 5 / 25 / 100 mm) y 0.99912 / 0.99982 / 0.99996 (Hordijk).

### 7. Ley de evolución y disipación

En carga la apertura alcanza su máximo histórico y la grieta sigue la ley de ablandamiento: $w_{\text{max},i} \leftarrow \max(w_{\text{max},i},\, h\,e_{ii})$, y si $w_i = w_{\text{max},i}$ crece, $\sigma_{ii} = f(w_i)$. Las flexibilidades $c_i$ y $c_s$ sólo crecen ($dc/dw = (f - w f')/f^2 > 0$ porque $f' \le 0$).

Disipación por unidad de volumen: $\dot{\mathcal D} = \frac{1}{2h}\bigl(\langle\sigma_{nn}\rangle_+^2\,\dot c_n + \langle\sigma_{tt}\rangle_+^2\,\dot c_t + \tau_{nt}^2\,\dot c_s\bigr) \ge 0$. En tracción pura hasta la separación completa (calculado): $\int\frac{f^2}{2h}\,dc = \int_0^\infty\frac{f - wf'}{2h}\,dw = \frac{G_f}{h}$, exacta en esfuerzo plano y en deformación plana (sólo entra $\sigma_{nn} = f(w)$; con Hordijk, 1.000000000000000·$G_f$, calculado por el teórico). Con cortante, la apertura libera además la energía de cortante almacenada en la grieta ($\tfrac{1}{2h}\tau^2\,\Delta c_s$). **En modo mixto la grieta disipa más que $G_f$** (calculado por el teórico): con deslizamiento proporcional a la apertura, $s = \lambda w$, disipa $(1 + \lambda^2/2)\,G_f$ por unidad de área, sin depender de $h$. En la validación de modo mixto, $D/(G_f A) > 1$ no es sólo bloqueo.

Al agrietarse no hay salto de rigidez a cortante: $c_s = 2w/f(w) \to 0$ cuando $w \to 0$ (con $\beta$ constante el módulo de cortante salta de $G$ a $\beta G$, lo que empeora la convergencia, Rots 1988 p. 87).

### 8. Integración implícita (bloque 1, sin aplastamiento)

Dado $\boldsymbol\varepsilon_{n+1}$ y la historia confirmada ($w_{\text{max},i}$, dirección de la grieta), se resuelve el equilibrio entre el concreto entre grietas y las grietas, $\mathbf N^\top\boldsymbol\sigma(\boldsymbol\varepsilon_{n+1} - \mathbf N\mathbf e) = \mathbf t^{cr}(\mathbf e)$, con las deformaciones de la grieta $\mathbf e = [e_{nn}, e_{tt}, \gamma_{nt}]$ como incógnitas.

**El cortante se desacopla** (calculado por el teórico): $\mathbf N^\top\mathbf C_e\mathbf N$ es diagonal por bloques en los ejes de la grieta, en las dos hipótesis. Por eso el problema local es:

1. **Problema normal de dos incógnitas** $(w_n, w_t)$: $\mathbf A\,(\boldsymbol\varepsilon_\text{loc} - \mathbf w/h) = \mathbf F(\mathbf w)$, con $\boldsymbol\varepsilon_\text{loc} = (\varepsilon_{nn}, \varepsilon_{tt})$ la deformación total en los ejes de la grieta, $\mathbf A = \begin{bmatrix} C_{11} & C_{12} \\ C_{12} & C_{11}\end{bmatrix}$ el bloque normal de $\mathbf C_e$ en esos ejes y $F_i(w_i)$ la ley de cada dirección: sin grieta o cerrada ($w_i = 0$, tramo vertical hasta $f_t$ o con $\sigma_{ii} \le 0$), en descarga ($F_i = w_i/c_i$, hasta $w_{\text{max},i}$) o en carga ($F_i = f(w_i)$, $w_i > w_{\text{max},i}$).
2. **Cortante en forma cerrada**, con la $w_\text{max}$ ya actualizada: $\tau_{nt} = \gamma_{nt}^\text{tot}/(1/G + c_s/h)$, con $\gamma_{nt}^\text{tot}$ la deformación angular total en los ejes de la grieta. Si $c_s = 0$, $\tau_{nt} = G\,\gamma_{nt}^\text{tot}$.

**Unicidad** (calculado por el teórico). Cada $F_i$ es un grafo monótono por tramos (vertical en $w = 0$, secante hasta $w_\text{max}$, envolvente convexa con pendiente $\ge f'(0)$). El operador $\mathbf F(\mathbf w) + \mathbf A\,\mathbf w/h$ es estrictamente monótono si y sólo si $(C_{11} - C_{12})/h > \lvert f'(0)\rvert$, y $C_{11} - C_{12} = E/(1+\nu)$ en las dos hipótesis: **hace falta $h < E/((1+\nu)\lvert f'(0)\rvert)$**, es decir $\ell_{ch}/(1+\nu)$ con la exponencial y $0.738\,\ell_{ch}/(1+\nu)$ con Hordijk. Es el modo antisimétrico: una grieta abre y la otra descarga. Entre ese valor y el límite de una dirección de §2 el teórico midió hasta cinco soluciones con las dos direcciones cerca del pico (doble inicio casi equibiaxial); con mallas usuales ($h \ll \ell_{ch}$) no ocurre. Decidido por el usuario (H1): ese intervalo se rechaza con el límite de §2.

**Algoritmo.** Bajo ese límite el operador es monótono y la rama de cada dirección se decide **evaluando el residuo en los quiebres** ($w = 0$ y $w = w_\text{max}$), no por los signos tras converger; en la rama de carga la función es creciente y convexa, y Newton desde el lado derecho baja de forma monótona (el argumento de §2b del modelo 1), con salvaguarda de bisección dentro del intervalo de la rama. Con dos direcciones en carga, Newton de 2×2 protegido. Convergencia con residuo ≤ tolerancia·$f_t$. (El Newton con conjunto activo literal, que cambia de estado tras converger, no convergió en 0.4–3.2 % de los estados al azar con $h \ge 0.6\,h_\text{lím}$; medido por el teórico, H2.) Lo programado (bloque 1) se describe en «Implementación»: el problema de dos incógnitas se resuelve anidado, con un problema escalar monótono en cada nivel, y la decisión en los quiebres usa la tolerancia de admisibilidad del ADR 0006 (una dirección carga si el esfuerzo de prueba supera la envolvente en $w_\text{max}$ en más de `tol`, y queda cerrada con $\sigma_{ii} \le$ `tol`).

**Inicio de la grieta:** sin grieta, el esfuerzo de prueba es elástico; si su esfuerzo principal mayor supera $f_t$, la grieta se forma con el ángulo $\theta$ de §1 calculado con ese esfuerzo de prueba y se resuelve el paso con la grieta en carga (error de primer orden en el paso de inicio, como en Rots 1988). En esos ejes $\tau_{nt} = 0$ y, por el desacople, sigue siéndolo en el paso de inicio; la respuesta es continua al cruzar el umbral. En deformación plana $\sigma_{zz} = \nu(\sigma_1 + \sigma_2)$ no supera $\sigma_1 > 0$. Un intento de paso rechazado no fija $\theta$: el ángulo sólo se confirma con el paso.

**Tangente consistente:** $\mathbf C_\text{tan} = \mathbf C_b - \mathbf C_b\mathbf N_a\,(\mathbf N_a^\top\mathbf C_b\mathbf N_a + \mathbf D_\text{tan})^{-1}\mathbf N_a^\top\mathbf C_b$, con $\mathbf C_b$ la tangente del concreto entre grietas ($\mathbf C_e$ sin aplastamiento), $\mathbf N_a$ las columnas de $\mathbf N$ de las direcciones abiertas (una dirección sin grieta o cerrada tiene $e_{ii} = 0$ y su columna se quita; también la de cortante si $c_s = 0$, porque la forma en rigidez no admite $D = \infty$) y $\mathbf D_\text{tan} = \partial\mathbf t^{cr}/\partial\mathbf e$: en carga, $h f'(w) < 0$ en la diagonal normal y el acoplamiento $\partial\tau_{nt}/\partial e_{ii} = -\tau_{nt}\,h\,(dc_s/dw_i)/c_s$. **No es simétrica** en carga activa con cortante (asimetría hasta 0.32·max|C| con dos grietas en carga, medida por el teórico; la fórmula coincide con diferencias finitas a 10⁻⁹–10⁻¹¹); en carga normal pura y a historia fija sí. `IS_SYMMETRIC = False`. En el paso de inicio $\theta$ depende de $\boldsymbol\varepsilon$ y la fórmula no incluye $\partial\theta/\partial\boldsymbol\varepsilon$: ahí la tangente no es la consistente.

**Deformación plana:** $\mathbf C_e$ es la elasticidad de deformación plana del plano; $\sigma_{zz} = \nu(\sigma_{xx}+\sigma_{yy})$ sin aplastamiento, o la del concreto entre grietas con $d^-$ (§4).

### 8b. Integración implícita con aplastamiento (bloque 2)

Con $f_c$ y $G_c$ el concreto entre grietas sigue la ley explícita de §4, $\boldsymbol\sigma_b(\boldsymbol\varepsilon^{co}) = \bar{\boldsymbol\sigma} - d^-(\kappa^-)\,\bar{\boldsymbol\sigma}^-$, con $\bar{\boldsymbol\sigma} = \mathbf C_e:\boldsymbol\varepsilon^{co}$ y $\kappa^- = \max\bigl(\kappa^-_n,\ \varepsilon^-_\text{eq}(\bar{\boldsymbol\sigma}^-)\bigr)$. Sin $f_c$ y $G_c$, todo es el bloque 1 (§8).

**Problema local.** Las incógnitas siguen siendo las deformaciones de la grieta $\mathbf e = [e_{nn}, e_{tt}, \gamma_{nt}]$ (aperturas $w_i = h\,e_{ii}$, deslizamiento $s = h\,\gamma_{nt}$), con $\boldsymbol\varepsilon^{co} = \boldsymbol\varepsilon_{n+1} - \mathbf N\mathbf e$, $\boldsymbol\sigma = \boldsymbol\sigma_b(\boldsymbol\varepsilon^{co})$ y $[\sigma_{nn}, \sigma_{tt}, \tau_{nt}] = \mathbf N^\top\boldsymbol\sigma$. Es el problema de §8 escrito sin suponer el desacople:

- **Dirección normal $i \in \{n, t\}$**, complementariedad: $w_i \ge 0$, $\ g_i(w_i) - \sigma_{ii} \ge 0$, $\ w_i\,\bigl(g_i(w_i) - \sigma_{ii}\bigr) = 0$, con $g_i$ la ley de §2 con la historia confirmada: $g_i(w) = w/c_i$ hasta $w_{\text{max},i}$ y $f(w)$ con el tope después; en una dirección sin grieta, $g_i(0) = f_t$ (cerrada hasta $f_t$) y la envolvente desde ahí.
- **Cortante:** $s = c_s\,\tau_{nt}$ ($s = h\,\gamma_{nt}$), con $c_s = 2(c_n + c_t)$ de las aperturas máximas actualizadas. Residuo $(c_s\,\tau_{nt} - s)/(c_s + h/G)$: en unidades de esfuerzo y bien condicionado desde $c_s = 0$ (da $s = 0$ sin dividir) hasta $4\,c_\text{max}$ (H2 de la revisión del bloque 2).
- $d^-$ no es incógnita: es función explícita de $\mathbf e$ a través de $\boldsymbol\varepsilon^{co}$ (§4). $\kappa^-_{n+1}$ sale de la solución final, no del esfuerzo de prueba.

**Algoritmo.**

1. **Predictor sin aplastamiento nuevo.** Si $d^-_n = 0$ (es decir, $\kappa^-_n \le f_c/(3E)$), se resuelve el problema del bloque 1 (concreto elástico, §8). Si en esa solución $d^-(\kappa^-_{n+1}) = 0$, es la solución del paso. Así el bloque 2 reproduce exactamente el bloque 1 mientras el concreto no se aplasta.
2. **Con aplastamiento** ($d^-_n > 0$, o el predictor da $d^-(\kappa^-_{n+1}) > 0$): Newton sobre las componentes activas de $\mathbf e$ (hasta tres). La complementariedad se escribe con la función mínimo, $\min\bigl(k\,w_i,\ g_i(w_i) - \sigma_{ii}\bigr) = 0$ con $k = C_{11}/h$ para igualar unidades (Newton semiliso: un conjunto activo reevaluado en cada iteración). En los quiebres ($w = 0$, $w = w_\text{max}$, el tope) las decisiones usan la tolerancia del ADR 0006, como el bloque 1 («Implementación»). Arranque: las deformaciones de la grieta del paso confirmado, es decir, la solución se elige por continuidad con el paso anterior (Q5, decidida por el usuario: sin límite ni aviso nuevos).
   - **Cascada** (H1 de la revisión del bloque 2, *medido* por el teórico): (a) Newton de paso completo desde el paso confirmado; (b) si no converge, el mismo con búsqueda lineal sobre la norma del residuo; (c) si tampoco, desde $\mathbf e = 0$; sólo entonces la bisección del paso (abajo). La búsqueda lineal sola, como primer intento, empeora la robustez con pasos grandes: la norma del residuo de la función mínimo no es diferenciable y la dirección de Newton no siempre la hace bajar; con una grieta ancha que descarga en un paso grande se queda en $w \approx 1$ mm con residuo de 2–50 kPa (fallos en hasta 0.56 % de los pasos con $\Delta\varepsilon \approx 1$–$2\cdot10^{-3}$, frente a 0.09 % con paso completo y 0 con la cascada; 0 en 1.24 millones de pasos con $\Delta\varepsilon \approx 1$–$3\cdot10^{-4}$). El arranque desde el paso confirmado es esencial: desde $\mathbf e = 0$ fallan 0.02–3.3 % de los pasos, así que sólo sirve de último recurso.
   - Jacobiano: $\mathbf N_a^\top\mathbf C_b\mathbf N_a + \mathbf D_\text{tan}$, con $\mathbf N_a$ y $\mathbf D_\text{tan}$ los de §8 y $\mathbf C_b = \partial\boldsymbol\sigma_b/\partial\boldsymbol\varepsilon^{co}$ la tangente algorítmica del modelo 2 con $d^+ = 0$ ([TensionCompressionDamage2D](TensionCompressionDamage2D.md) §7): $\mathbf C_b = \mathbf C^d - \bar{\boldsymbol\sigma}^-\otimes\partial d^-/\partial\boldsymbol\varepsilon^{co}$, con $\mathbf C^d = \mathbf C_e - d^-(\mathbb I - \mathbf Q^+)\,\mathbf C_e$ y, en deformación plana, el término de $\bar\sigma_{zz}$ en $\partial d^-/\partial\boldsymbol\varepsilon^{co}$; el término de $d^-$ sólo en carga de compresión ($\varepsilon^-_\text{eq} > \kappa^-_n$ y $d^- > 0$).
   - Convergencia: el residuo de cada ecuación, en unidades de esfuerzo, ≤ tolerancia·$f_t$.
   - **Si la cascada no converge**, el punto devuelve un esfuerzo no finito, el corrector abandona el intento sin resolver con él y el solver repite el paso con un incremento menor (bisección); el aviso de la bisección lleva el motivo («el residuo no es finito: algún punto material no devolvió un esfuerzo válido»), no «tangente singular».

**Unicidad.** Sin aplastamiento ($d^- = 0$) vale la de §8: única bajo el límite de $h$ de §2. Con $d^- > 0$ no está garantizada (§4), por dos mecanismos (H3 de la revisión del bloque 2): (1) a daño fijo, una grieta en carga sobre concreto aplastado cuya rigidez proyectada sobre la grieta baja por debajo de $h\,\lvert f'\rvert$; (2) la grieta ablanda y el aplastamiento crece en el mismo paso: dos ablandamientos en serie. En la biela paralela a la grieta (ejes de la grieta principales) el problema normal es el del bloque 1 y la unicidad sí vale (§4). **Medido por el teórico:** en 1.9 millones de pasos de caminos al azar con continuidad (aplastamiento hasta $d^-$ = 0.99, una y dos grietas, cortante, giro de direcciones, inversión de carga, $h$ = 10–100 mm), ningún estado con más de una solución ni dependencia del arranque, y el predictor nunca coexistió con otra solución; también 0 en 72 caminos dirigidos (biela a 20°–70° de la normal, tracción a 5°–30° de la dirección aplastada). Los estados con varias soluciones (hasta 7) sólo aparecieron construidos a propósito y con concreto casi sin capacidad: 0 con $d^- \le 0.95$ y $h \le 50$ mm por el mecanismo (1), y 0.03–0.17 % de los estados construidos desde $d^-$ = 0.9 ($h$ = 50 mm) o 0.8 ($h$ = 100 mm) por el (2); son órdenes de magnitud del muestreo, no probabilidades. La cota «$h < 0.1\,\ell_{ch}$ con $d^- = 0.9$» de §4 es pesimista. Tratamiento (Q5, decidida por el usuario): selección por continuidad, sin límite ni aviso nuevos.

**Tangente consistente.** La fórmula de §8 con $\mathbf C_b$ la tangente algorítmica de arriba:

$$\mathbf C_\text{tan} = \mathbf C_b - \mathbf C_b\mathbf N_a\,(\mathbf N_a^\top\mathbf C_b\mathbf N_a + \mathbf D_\text{tan})^{-1}\mathbf N_a^\top\mathbf C_b,$$

que sale de linealizar $\boldsymbol\sigma = \boldsymbol\sigma_b(\boldsymbol\varepsilon - \mathbf N\mathbf e)$ y $\mathbf N_a^\top\boldsymbol\sigma = \mathbf t^{cr}_a(\mathbf e)$ y vale también con $\mathbf C_b$ no simétrica. Sin direcciones activas, $\mathbf C_\text{tan} = \mathbf C_b$ (el modelo 2 con $d^+ = 0$). No simétrica. Hereda del modelo 2 los quiebres de $\mathbf Q^+$: $H(0) = 0$ determinista con un esfuerzo principal efectivo nulo (§7 de su spec).

**Inicio de la grieta con aplastamiento.** El esfuerzo de prueba es $\boldsymbol\sigma_b(\boldsymbol\varepsilon_{n+1})$ con $\mathbf e = 0$ (el modelo 2 con $d^+ = 0$). Si su principal mayor supera $f_t$, la grieta se forma con el $\theta$ de §1, que es el de $\bar{\boldsymbol\sigma}$ (§2), y el paso se resuelve con el algoritmo de arriba. Como en el bloque 1, un intento rechazado no fija $\theta$ ni $\kappa^-$.

### 9. Integración IMPL-EX

Como en los modelos 1 y 2 ([ADR 0034](../adr/0034-pseudo-tiempo-de-los-solvers-incrementales.md)): la apertura máxima se extrapola con los dos últimos pasos confirmados, $\tilde w_{\text{max},i} = \max\bigl(\tilde w_{\text{max},i}^{(n)},\ w_{\text{max},i,n} + r\,(w_{\text{max},i,n} - w_{\text{max},i,n-1})\bigr)$, $r = \Delta\tau_{n+1}/\Delta\tau_n$; con ella las flexibilidades $\tilde c_i$, $\tilde c_s$ quedan conocidas en el paso y el esfuerzo sale de la ley secante de §6 en forma cerrada: el problema normal de 2×2 en sus a lo sumo cuatro combinaciones de grietas abiertas o cerradas (se toma la única que cumple los signos, por la convexidad de χ; con $\sigma_{ii} \le$ `tol`, la tolerancia de admisibilidad del ADR 0006, la dirección se toma cerrada, lo que incluye el empate $\sigma_{ii} = 0$ exacto, como el $H(0) = 0$ del modelo 2) y el cortante en forma cerrada. Tangente: la secante, simétrica y definida positiva; con las dos direcciones en el tope su valor propio mínimo es el del cortante con $c_s = 4\,c_\text{max}$, $(1-\omega_\text{max})\,E/(4\,\omega_\text{max})$ (2.5·10⁻⁷·$E$ con `damage_max` por omisión; medido por el teórico, verificación del bloque 1). La historia se actualiza con la integración implícita de §8 sobre la deformación del paso, para la extrapolación siguiente; indicador de error `damage_error` $= \max_i\lvert\omega_i(w_{\text{max},i}) - \omega_i(\tilde w_{\text{max},i})\rvert$. Una grieta nueva se detecta en esa etapa implícita y entra en la extrapolación desde el paso siguiente (error de primer orden en el inicio, como el modelo 1). La etapa implícita es el problema de §8, así que su límite de $h$ y su algoritmo valen también aquí.

**Con aplastamiento (bloque 2).** $\kappa^-$ se extrapola como en el modelo 2: $\tilde\kappa^-_{n+1} = \kappa^-_n + r\,(\kappa^-_n - \kappa^-_{n-1})$ y $\tilde d^-_{n+1} = \max\bigl(\tilde d^-_n,\ d^-(\tilde\kappa^-_{n+1})\bigr)$. Con $\tilde d^-$, $\tilde c_i$ y $\tilde c_s$ congelados, el esfuerzo sale del problema local de §8b con la historia fija: el concreto entre grietas $\bar{\boldsymbol\sigma} - \tilde d^-\bar{\boldsymbol\sigma}^-$ y la ley secante de cada grieta (sin rama de envolvente).
- Si $\tilde d^- = 0$, el concreto es elástico y vale la forma cerrada de arriba (bloque 1).
- Si $\tilde d^- > 0$ **no hay forma cerrada** (§4): el problema sigue siendo no lineal por la partición espectral, aunque la historia esté fija. Se resuelve con el Newton de §8b; la tangente es la consistente de ese problema, la fórmula de §8b con $\mathbf C_b = \mathbf C^d$ (la del modelo 2 a daño fijo) y $\mathbf D$ la secante. **No es simétrica** con $\nu \ne 0$ (como el modelo 2, §9 de su spec), y el Newton global itera dentro del paso como en el modelo 2 (§8 de su spec).
- **Unicidad con la historia fija:** la secante de la grieta es monótona y el concreto a daño fijo $\tilde d^-$ es estrictamente monótono mientras la parte simétrica de $\mathbf C^d$ es definida positiva, así que la solución es única. Eso vale hasta el umbral (calculado por el teórico, H5)
  $$d^* = \frac{2\sqrt{1-\nu^2}}{1 + \sqrt{1-\nu^2}}$$
  en esfuerzo plano, y la misma fórmula con $\nu/(1-\nu)$ en lugar de $\nu$ en deformación plana: con $\nu$ = 0.2, $d^*$ = 0.98980 y 0.98387 (0.99177 y 0.98765 con $\nu$ = 0.18). Por encima se pierde la garantía (herencia del modelo 2), no la unicidad observada: 0 estados con varias soluciones con $\tilde d^-$ hasta 0.9999 (*medido* por el teórico, 1500 estados por caso). Desde $\mathbf e = 0$ el Newton falla en 0.5–1.5 % de los estados con $\tilde d^- = 0.9999$, así que la cascada de §8b se usa también aquí.
- Indicador de error: `damage_error` $= \max\bigl(\max_i\lvert\omega_i(w_{\text{max},i}) - \omega_i(\tilde w_{\text{max},i})\rvert,\ \lvert d^-(\kappa^-_{n+1}) - \tilde d^-\rvert\bigr)$.
- La historia ($w_{\text{max},i}$ y $\kappa^-$) se actualiza con el problema implícito de §8b sobre la deformación del paso.

### 10. Caveats

- **Sobrestimación en vigas con entalla en modo I: las grietas nacen inclinadas** (validación del 2026-11-25; limitación documentada por decisión del usuario, Diálogo).
  - *Medido*: viga de Petersson, Quad4 h = 10 / 5 / 2.5 mm, IMPL-EX: pico 868.6 / 928.4 / 908.0 N (6.6 % entre mallas) frente a 746–750 N de los modelos 1 y 2; $D(8\text{ mm}) = 2.43 / 2.43 / 2.91\,G_F A$ frente a 0.94–0.95. Viga de Körmeling y Reinhardt de Slobbe (2015), Hordijk: Quad4 1606 N (+15 %) y Quad8 3×3, el elemento de Slobbe, 1540 N (+11.5 %) frente a 1381–1392 N de Slobbe y 1394 N del modelo 1; a 0.3 mm de flecha, 1052 frente a 757 N con Quad8 (+39 %). Ensayo de Van Vliet (tracción casi uniforme): objetivo entre mallas, pero $D(1.2\text{ mm}) = 1.14\,G_f A$ y subiendo (grieta principal 1.02, segundas grietas 0.05, deslizamiento 0.07). Reportes: `informes/reportes/reporte_modelo3_vanvliet_petersson_2026-11-24.pdf`, `informes/modelo3_slobbe_yaml/resumen.md`.
  - *Causa*: al nacer la grieta junto a la punta de la banda, el estado de esfuerzos es casi equibiaxial ($\sigma_{II}/\sigma_I$ mediana 0.74; en el campo de modo I la diferencia sobre el eje es sólo el término $T$), y la dirección principal queda mal condicionada: un cortante del 3–11 % de $f_t$ la gira 19–35° (*calculado*). $\theta$ sigue bien esa dirección (mediana 0.2–1.2°, *medido* con Quad4, Quad8 y Tri3, implícita e IMPL-EX, con el mismo pico): no lo causan el elemento, la integración ni la asignación del ángulo, y no se quita al refinar. La grieta fija conserva la inclinación; la apertura posterior de la banda la carga en modo mixto, abre la segunda grieta y hace trabajar la retención de §3: con vecinos rígidos y la grieta a 20–30°, el punto material da pico +11–25 % y $1.9$–$2.5\,G_f$ (*calculado*, `informes/teorico/verificacion_teorica_2026-11-25_modelo3_petersson.md`).
  - *Publicado* (*leído*): Rots et al. 1985 (HERON 30(1), pp. 17–18), en esta misma viga de Petersson: «off-centre Gauss-points have non-vertical cracks»; usan $\beta = 0.001$ porque con retención alta la resistencia residual queda «considerably overestimated». Rots 1988, p. 81: «fixed cracks perform adequately only in conjunction with zero or almost zero β»; Rots y Blaauwendraad 1989 (HERON 34(1), §5.2): bloqueo de esfuerzos por rotación del esfuerzo principal en dirección inclinada. Slobbe 2015 lo compensa con un ancho de banda empírico $l_p(\theta)$, que llama «a quantitative compensation for the numerically obtained incorrect energy dissipation» (p. 68). En un modelo de articulación, esa compensación cierra la mitad o dos tercios de la brecha, y la proyección geométrica de Govindjee sola cierra poco (*calculado*, `informes/teorico/verificacion_teorica_2026-11-25_modelo3_slobbe.md`). No se adopta ninguna: no corrigen la causa.
  - Remedios publicados que cambiarían la física (fuera de esta spec): retención casi nula (contradice Q2), grieta rotante o multidireccional con ángulo umbral (fuera de alcance por Q1), «rotating to fixed» de DIANA (§46.5.12), ancho de banda dependiente de $\theta$; y de modelado: puntos de integración sobre el eje de la banda (Rots 1988, p. 74), banda predefinida. Búsqueda: `informes/bibliotecario/busqueda_grieta_fija_modo_I_2026-11-25.md`.
- **Reconstrucción de esfuerzo en direcciones inclinadas.** Con grietas fijas ortogonales, el esfuerzo principal puede superar $f_t$ en una dirección distinta de $\mathbf n$ y $\mathbf t$ sin que se forme una grieta nueva; con retención de cortante es la causa del bloqueo «por esfuerzo inclinado» (Rots 1988, pp. 24, 81–82, 86). Se registra en una prueba y se reporta en la validación.
- **Bloqueo por compatibilidad entre elementos** (Rots 1988, §4.2.3): cinemático, como en los modelos 1 y 2; ningún cambio de material lo quita.
- **Sesgo direccional de la malla** (Rots 1988, pp. 87–88; Jirásek y Zimmermann 1998, p. 280).
- Límite de $h$ de §2; tope de la flexibilidad de §6.
- IMPL-EX: error de primer orden en el inicio de la grieta; pasos finos cerca del pico.
- Descarga secante: sin deformación residual al cerrar la grieta (de Borst et al. 2012, p. 206).
- **Con aplastamiento** (bloque 2): sin unicidad garantizada del problema local sobre concreto aplastado con grietas inclinadas respecto de la compresión (§4, §8b, Q5); sin potencial ni trabajo nulo en ciclos (§4); IMPL-EX con tangente no simétrica y unicidad garantizada sólo hasta $\tilde d^- = d^*$ (0.9898 con $\nu$ = 0.2 en esfuerzo plano, §9). Heredados del modelo 2 ([TensionCompressionDamage2D](TensionCompressionDamage2D.md) §10): pico mayor que $f_c$ en deformación plana uniaxial y con confinamiento (cono de Drucker-Prager abierto, sin tapa); energía total en compresión mayor que $G_c A$ por la rama ascendente; compresión de una barra con $\nu \ne 0$ sin localización clara; con Quad8 y Quad9 la $h$ lleva el factor 13/18 también en $\kappa_u$ (pendiente del usuario en el modelo 2, H1); elementos grandes sin límite en compresión.

---

## Contrato de implementación

```yaml
name: FixedSmearedCrack2D
kind: material
status: implemented      # draft → implemented → validated. No validado: limitación documentada en vigas con entalla en modo I (§10, 2026-11-25)

interface:
  strain_dim: 3
  primary_state_var: "damage_n"
  is_symmetric: false    # tangente consistente no simétrica en carga activa (§8, §8b); la de IMPL-EX es simétrica sin aplastamiento y no simétrica con él (§9)
  state_schema: {cracked: [], theta: [], w_max_n: [], w_max_t: [], w_n: [], w_t: [], slip: [], damage_n: [], damage_t: [], kappa_c: [], damage_c: [], h: [], w_max_n_prev: [], w_max_t_prev: [], kappa_c_prev: [], damage_error: []}
  batch_kernel: true
  uses_pseudo_time: "true con integration='impl-ex' (ADR 0034)"

admission:
  textbooks:
    - de_borst_crisfield: "§6.8.1 «Elasticity-based smeared crack models», pp. 201–206: grieta fija (ecs. 6.112–6.121, retención β en 6.115 y 6.117) y descomposición multidireccional (ecs. 6.127–6.135)"
  general_codes:
    - diana: "Total Strain Crack Model, crack orientation «fixed» (Theory 10.4 §46.5)"
    - abaqus: "Concrete smeared cracking (User's Guide §23.6.1, Theory §4.5.1) y brittle cracking (§23.6.2, Theory §4.5.3): grietas fijas ortogonales con retención de cortante"
    - ansys: "SOLID65, TB,CONCR: hasta tres grietas ortogonales con coeficientes de transferencia de cortante β_t y β_c"

parameters:
  - { name: E, type: float, required: true, desc: "módulo de Young del concreto entre grietas [Pa]" }
  - { name: nu, type: float, required: true, desc: "relación de Poisson del concreto entre grietas" }
  - { name: f_t, type: float, required: true, desc: "resistencia a tracción [Pa]" }
  - { name: G_f, type: float, required: true, desc: "energía de fractura [N/m]" }
  - { name: hypothesis, type: str, required: false, desc: "plane_stress (por omisión, como los demás materiales 2D) | plane_strain" }
  - { name: softening, type: str, required: false, desc: "exponential (en w, por omisión) | hordijk" }
  - { name: f_c, type: float, required: false, desc: "resistencia a compresión [Pa]; con G_c activa el aplastamiento (§4)" }
  - { name: G_c, type: float, required: false, desc: "energía de aplastamiento [N/m]; va con f_c" }
  - { name: beta, type: float, required: false, desc: "f_bc/f_c de Drucker-Prager (1.16 por omisión, como el modelo 2)" }
  - { name: damage_max, type: float, required: false, desc: "tope del daño uniaxial equivalente de cada banda (el de los modelos 1 y 2)" }
  - { name: integration, type: str, required: false, desc: "implicit (por omisión) | impl-ex (como los modelos 1 y 2)" }
  - { name: density, type: float, required: false, desc: "densidad [kg/m³]" }

state_variables:
  - { name: cracked, type: float, shape: [], desc: "1 si ya se formó la primera grieta (direcciones fijas), 0 si no" }
  - { name: theta, type: float, shape: [], desc: "ángulo de la normal n de la primera grieta, antihorario desde +x [rad]" }
  - { name: w_max_n, type: float, shape: [], desc: "apertura máxima histórica de la grieta normal a n [m]" }
  - { name: w_max_t, type: float, shape: [], desc: "apertura máxima histórica de la grieta normal a t [m]" }
  - { name: w_max_n_prev, type: float, shape: [], desc: "w_max_n del paso confirmado anterior (IMPL-EX)" }
  - { name: w_max_t_prev, type: float, shape: [], desc: "w_max_t del paso confirmado anterior (IMPL-EX)" }
  - { name: w_n, type: float, shape: [], desc: "apertura actual de la grieta normal a n [m] (post-proceso)" }
  - { name: w_t, type: float, shape: [], desc: "apertura actual de la grieta normal a t [m] (post-proceso)" }
  - { name: slip, type: float, shape: [], desc: "deslizamiento actual de las grietas, h·γ_nt [m] (post-proceso)" }
  - { name: damage_n, type: float, shape: [], desc: "daño uniaxial equivalente de la banda en n, ω_n (§6); variable principal" }
  - { name: damage_t, type: float, shape: [], desc: "ídem en t" }
  - { name: kappa_c, type: float, shape: [], desc: "máxima deformación equivalente de compresión (§4)" }
  - { name: damage_c, type: float, shape: [], desc: "daño de compresión d⁻ del concreto entre grietas (§4)" }
  - { name: kappa_c_prev, type: float, shape: [], desc: "kappa_c del paso confirmado anterior (IMPL-EX)" }
  - { name: h, type: float, shape: [], desc: "longitud característica del punto, sembrada por el elemento" }
  - { name: damage_error, type: float, shape: [], desc: "indicador de IMPL-EX (§9)" }

conventions:
  voigt: "[ε_xx, ε_yy, γ_xy], γ_xy = 2ε_xy; σ = [σ_xx, σ_yy, σ_xy]"
  units: "SI coherente: E, f_t en Pa; G_f en N/m; h en m"
  sign: "σ > 0 tracción; w ≥ 0 apertura; ángulo de la grieta antihorario desde +x"

validity:
  - "pequeñas deformaciones; carga monótona o descargas moderadas (descarga secante)"

out_of_scope:
  - "grieta rotante (Q1, 2026-11-24)"
  - "grietas no ortogonales (multidireccional con ángulo umbral, Rots 1988)"
  - "grieta normal a z en deformación plana"
  - "factor ajustable s_F, β constante o tope de resistencia a cortante (Q2)"
  - "reducción de la resistencia a compresión por agrietamiento transversal, Vecchio–Collins (Q3)"
  - "recuperación de la rigidez a cortante al cerrar la grieta"
  - "deformaciones plásticas (modelo 4)"
  - "3D (spec aparte, FixedSmearedCrack3D, con tres grietas ortogonales)"

acceptance:
  material_point:
    - name: elastic_before_cracking
      setup: "estados con esfuerzo principal mayor < f_t, esfuerzo plano y deformación plana"
      expect: "σ = C_e:ε; tangente C_e"
      tol_rel: 1.0e-12
    - name: uniaxial_tension_curve
      setup: "tracción uniaxial de esfuerzo en esfuerzo plano, grieta a 0° y a 30°, leyes exponential y hordijk"
      expect: "σ = f(w) con ε = σ/E + w/h; contracción lateral −ν·σ/E (sólo la del concreto entre grietas: abrir la grieta no contrae); la grieta conserva el ángulo inicial"
      tol_rel: 1.0e-10
    - name: plane_strain_uniaxial
      setup: "deformación plana, tracción con σ_yy = 0"
      expect: "σ_zz = ν·σ_xx; ε_xx = (1 − ν²)·σ/E + w/h; la misma σ–w que en esfuerzo plano"
      tol_rel: 1.0e-10
    - name: tension_energy
      setup: "tracción uniaxial monótona hasta el tope; h ∈ {5, 25, 100} mm; exponential y hordijk; esfuerzo plano y deformación plana"
      expect: "energía disipada = D_tope de §6 (en deformación plana sin el factor 1 − ν² del modelo 1); ≥ 2·10⁴ pasos; control negativo: con la ley de otra G_f la prueba falla"
      tol_rel: 1.0e-6
    - name: mixed_mode_energy
      setup: "apertura y deslizamiento proporcionales, s = λ·w con λ ∈ {0.5, 1, 2}, hasta el tope; h ∈ {5, 25, 100} mm; G_f y h ≠ 1"
      expect: "energía disipada por unidad de área = (1 + λ²/2)·G_f (§7, con el descuento del tope), igual con las tres h (protege c_s en la disipación; Reglas §6)"
      tol_rel: 1.0e-5
    - name: dissipation_non_negative
      setup: "caminos al azar (semilla fija) con apertura, deslizamiento, descarga, cierre y reapertura, en las dos direcciones; sin aplastamiento (con d⁻ > 0 no hay potencial, §4)"
      expect: "disipación de cada paso ≥ −1e-12·energía elástica (§7)"
    - name: shear_retention_law
      setup: "grieta abierta hasta w_max y después deslizamiento con la apertura fija; w_max en el ablandamiento y en descarga"
      expect: "τ = f(w_max)/(2·w_max)·s; igual al de la opción «por daño» de DIANA en deformación total, G_red = E_red/(2(1 + ν·E_red/E)) (§3)"
      tol_rel: 1.0e-12
    - name: shear_objectivity
      setup: "el caso anterior con h ∈ {5, 25, 100} mm y la misma apertura w_max"
      expect: "la misma curva τ–s (h no entra en la ley de la grieta)"
      tol_rel: 1.0e-12
    - name: fixed_direction_and_inclined_stress
      setup: "grieta formada a +θ₀ y a −θ₀ (θ₀ ≠ 0, ejercita la columna de cortante de N); después se giran los ejes principales de la carga 30° y 60° con la grieta abierta"
      expect: "θ no cambia; τ_nt = γ_nt/(1/G + c_s/h) (§8) y el tensor completo con N, en forma cerrada; se registra el esfuerzo principal en dirección inclinada, que puede superar f_t sin grieta nueva (caveat de §10)"
      tol_rel: 1.0e-10
    - name: second_crack_orthogonal
      setup: "tracción biaxial: primero en n hasta ablandar, después en t"
      expect: "la segunda grieta se forma a θ₀ + 90° cuando σ_tt = f_t; cada dirección sigue su ley; flexibilidad a cortante c_s = 2(c_n + c_t); energía de las dos grietas 2·G_f/h menos el tope"
      tol_rel: 1.0e-10
    - name: cap_reached
      setup: "apertura más allá de w_tope, exponential y hordijk (con w ≥ w_c)"
      expect: "σ_nn = h·e_nn/c_max, c_s ≤ 4·c_max; sin división entre cero con hordijk y w ≥ w_c"
      tol_rel: 1.0e-12
    - name: local_problem_robustness
      setup: "estados al azar (semilla fija) cerca de los quiebres w = 0 y w = w_max, una y dos direcciones, h hasta 0.95 del límite de §8; y el doble inicio casi equibiaxial"
      expect: "residuo ≤ tolerancia en todos; la misma solución desde arranques distintos; el doble inicio con resultado determinista"
    - name: implicit_step_idempotent
      setup: "el mismo intento de paso dos veces en la integración implícita, incluido el paso de inicio de la grieta tras un intento rechazado; con y sin aplastamiento"
      expect: "estado de prueba idéntico; ni θ ni κ⁻ quedan fijados por el intento rechazado"
    - name: unloading_closure_reopening
      setup: "apertura hasta w_max, descarga, compresión y reapertura"
      expect: "descarga secante σ = f(w_max)·w/w_max; cerrada: rigidez E y Poisson ν en compresión; reapertura por la secante hasta w_max y después la envolvente; w_max nunca decrece"
      tol_rel: 1.0e-12
    - name: parallel_compression_open_crack
      setup: "esfuerzo plano, grieta abierta en n con aperturas w ∈ {0, 0.5, 1, 2}·w_c (o hasta el tope) y compresión ε_t paralela"
      expect: "σ_tt − ν·σ_nn = E·ε_t (esfuerzo plano) y σ_tt = E·ε_t/(1 − ν²) + ν·σ_nn/(1 − ν) (deformación plana), con cualquier apertura (§4; contraste con el caveat del modelo 2)"
      tol_rel: 1.0e-10
    - name: zero_work_fixed_history
      setup: "ν = 0.2, ciclos cerrados genéricos de deformación con la historia fija (w_max sin crecer), con grietas abiertas y cerradas, sin aplastamiento"
      expect: "∮σ:dε = 0 frente a la energía elástica máxima del ciclo (la ley deriva de χ, §6); control negativo: recuperar el cortante al cerrar la grieta debe dar ≥ 0.1"
      tol_rel: 1.0e-6
    - name: tangent_fd
      setup: "implícita en carga activa normal, normal con cortante, dos grietas, descarga y cerrada; esfuerzo plano y deformación plana; lejos de los cambios de estado y fuera del paso de inicio (§8). Bloque 2 (§8b): aplastamiento en carga sin grieta, con grieta cerrada, con grieta abierta paralela a la compresión y con grieta inclinada y cortante; en deformación plana con σ̄_zz < 0"
      expect: "|C_FD − C_tan| ≤ 1e-6·max|C_FD| (1e-5 en los casos con grieta abierta y aplastamiento en carga: p99 medido 1.5e-6); asimetría > 0 en carga con cortante (falla si se olvida ∂τ/∂e) y con aplastamiento; a historia fija y sin aplastamiento, simétrica; en deformación plana con aplastamiento, falla si se quita el término de σ̄_zz (error medido 1.8·max|C|)"
      tol_rel: 1.0e-6
    - name: compression_unchanged_below_crushing
      setup: "con f_c y G_c, caminos al azar (semilla fija) con grietas y compresión tal que κ⁻ ≤ f_c/(3E) en todo el camino"
      expect: "σ, C_tan y estado idénticos bit a bit a los del mismo caso sin f_c ni G_c (predictor de §8b), salvo kappa_c y kappa_c_prev, que con f_c se actualizan aunque no haya daño"
    - name: predictor_independence
      setup: "pasos al azar con aplastamiento incipiente (κ⁻ cerca de f_c/(3E)), resueltos con y sin el predictor de §8b"
      expect: "la misma solución (medido por el teórico: siempre)"
      tol_rel: 1.0e-10
    - name: crushing_parallel_to_open_crack
      setup: "esfuerzo plano; grieta en n abierta hasta w ∈ {0.1, 0.5, 1}·w_c (o hasta el tope) con la apertura fija, impuesta con ε_n(ε_t) = w/h + (σ_nn − ν·σ̄_tt)/E (fijar ε_n cierra la grieta); después compresión ε_t paralela hasta pasar el pico de aplastamiento; h ∈ {10, 50} mm"
      expect: "−σ_tt = σ_c(κ⁻) con κ⁻ = |E·ε_t + ν·σ_nn|/E (§4), tol 1e-10 (medido 1e-12): pico f_c con cualquier apertura; σ_nn = f(w) la del bloque 1; área de la rama descendente G_c/h con tolerancia propia 1e-3 (cuadratura)"
      tol_rel: 1.0e-10
    - name: crushing_parallel_plane_strain
      setup: "el caso anterior en deformación plana"
      expect: "pico = f_c/φ(ρ), ρ = ν·(1 − σ_nn/|σ̄_tt|) (§4): 1.1397·f_c con σ_nn = 0 y ν = 0.2, β = 1.16 (la tolerancia la limita el muestreo del camino)"
      tol_rel: 1.0e-6
    - name: crack_initiation_after_crushing
      setup: "compresión en una dirección hasta d⁻ ≈ 0.5, descarga y tracción en la dirección perpendicular y con ejes girados ±θ₀"
      expect: "la grieta se forma cuando σ̄_I = f_t, con θ = θ(σ̄) (§2)"
      tol_rel: 1.0e-10
    - name: closed_crack_crushing
      setup: "grieta abierta en n, descarga hasta cerrarla y compresión normal a la grieta hasta pasar el pico de aplastamiento"
      expect: "la misma curva que compression_like_model2 en compresión uniaxial; w_n = 0 y w_max,n sin cambio"
      tol_rel: 1.0e-10
    - name: crushing_then_cracking
      setup: "compresión uniaxial en x hasta d⁻ ≈ 0.5 en la rama descendente, descarga a cero y tracción en x"
      expect: "la grieta se forma con σ_xx = f_t y θ = 0 (tol 1e-10); d⁻ sin cambio en la tracción; energía disipada = D_tope de §6 con la tolerancia de tension_energy (1e-6; sin acoplamiento entre los dos daños, §4)"
      tol_rel: 1.0e-10
    - name: local_problem_robustness_crushing
      setup: "caminos al azar (semilla fija) con aplastamiento hasta d⁻ = 0.9, una y dos grietas, cortante, giro de las direcciones principales e inversión de carga; h ∈ {10, 25, 50} mm; esfuerzo plano y deformación plana; ≥ 60 pasos por tramo; y caminos entre vértices con grietas anchas que descargan (ε hasta f_t/E + w_c/h) en 3 pasos por tramo, h ∈ {10, 50, 100} mm, que sí necesitan la cascada de §8b (fallan sin ella o arrancando desde x = 0)"
      expect: "el Newton de §8b converge en todos los pasos de los caminos (arranque desde el paso confirmado, con la cascada); se informa el número de estados con más de una solución encontrada desde arranques distintos (dato, sin criterio)"
    - name: local_newton_failure_flagged
      setup: "un punto material de prueba que devuelve un esfuerzo no finito con pasos grandes, también en un elemento con todos sus grados de libertad prescritos; y FixedSmearedCrack2D con una deformación no finita"
      expect: "el material devuelve σ no finito sin lanzar; el corrector abandona el intento con el motivo «el residuo no es finito» (no «tangente singular»), aunque el NaN sólo llegue a filas de reacciones; el solver repite el paso más corto y el análisis termina"
    - name: confirmed_state_reevaluated_coupled
      setup: "caminos con grieta y aplastamiento en carga; tras cada paso, el estado confirmado se vuelve a evaluar con la misma deformación"
      expect: "w_max y κ⁻ no crecen y la tangente es la de descarga (la de una deformación apenas menor) a 1e-4; falla sin la banda δ o sin la tolerancia de κ⁻"
    - name: element_open_crack_and_crushing
      setup: "un Quad4 de 100 mm con tracción en x y compresión en y proporcionales (deformación uniforme), 200 pasos, implícita e IMPL-EX"
      expect: "sin bisecciones, ≤ 6 iteraciones por paso (8 con IMPL-EX); estado y reacciones iguales al punto material a 1e-8"
    - name: compression_like_model2
      setup: "con f_c y G_c, sin grieta: compresión uniaxial, biaxial igual e hidrostática; h ∈ {10, 50, 100} mm; esfuerzo plano y deformación plana; implícita e IMPL-EX"
      expect: "las mismas curvas y daños que TensionCompressionDamage2D con d⁺ = 0 y la misma integración (pico f_c, pico biaxial β·f_c, área de la rama descendente G_c/h, sin daño hidrostático)"
      tol_rel: 1.0e-10
    - name: snapback_limit
      setup: "h justo por encima y justo por debajo de ℓ_ch/(1+ν) (exponential) y de 0.738·ℓ_ch/(1+ν) (hordijk), esfuerzo plano y deformación plana"
      expect: "ValueError con mensaje claro por encima; se admite por debajo (§2, H1)"
    - name: units_invariance
      setup: "el mismo caso en Pa–m y en MPa–mm, con todos los parámetros físicos ≠ 1 (E, ν, f_t, G_f, f_c, G_c, β, h)"
      expect: "esfuerzos, aperturas, daños y energías iguales tras convertir unidades (Reglas §6)"
      tol_rel: 1.0e-12
    - name: batch_kernel_matches
      expect: "el kernel por lotes reproduce bit a bit σ, C_tan y la fila de estado de compute_state, implícita e IMPL-EX"
  impl_ex:
    - name: secant_tangent
      expect: "C_tan = secante de §6 con la historia extrapolada, igual a diferencias finitas con esa historia congelada; simétrica y definida positiva. Con d̃⁻ > 0 (§9): igual a diferencias finitas con la historia congelada, lejos de los quiebres; asimetría > 0"
      tol_rel: 1.0e-6
    - name: first_order_convergence
      setup: "como los modelos 1 y 2, en tracción con y sin cortante y en compresión con aplastamiento (con y sin grieta abierta paralela), pasos Δ, Δ/2, Δ/4, Δ/8"
      expect: "por forma del camino, como el modelo 2: con descarga y recarga el error máximo de σ frente a la implícita se reduce a la mitad al reducir el paso a la mitad (razón 2 ± 15 %); en compresión monótona (también con la grieta abierta paralela y σ_nn fija) la extrapolación de κ⁻ es exacta desde el segundo paso y el orden es 2 (razón 4)"
    - name: history_monotone
      setup: "carga, descarga parcial y recarga con extrapolación que se pasa, en tracción y en compresión con aplastamiento"
      expect: "w̃_max y d̃⁻ nunca decrecen; damage_error incluye |d⁻(κ⁻) − d̃⁻| (§9)"
    - name: bisection_idempotent
      expect: "el mismo intento de paso dos veces da el mismo estado de prueba"
    - name: outside_step_frozen
      expect: "compute_state fuera de un paso no cambia la historia"
  mesh_objectivity:
    - name: bar_in_tension
      setup: "barra con una columna debilitada, Tri3, Quad4, Tri6, Quad8, Quad9, 5–20 elementos; esfuerzo plano y deformación plana; ν = 0 y ν = 0.2; implícita e IMPL-EX"
      expect: "D/(G_f·A) = 1 ± 1 % (Quad8 y Quad9 con el factor de localización de los modelos 1 y 2) y carga posterior al pico igual entre mallas ± 1 %; con ν = 0.2 en deformación plana también 1 (sin el factor 1 − ν²)"
      tol_rel: 1.0e-2
    - name: bar_in_compression
      setup: "la barra del modelo 2 en compresión con f_c y G_c"
      expect: "trabajo de la banda desde su pico = G_c·A ± 2 % e igual entre mallas"
      tol_rel: 2.0e-2
  structural:
    - name: van_vliet_tension
      setup: "el ensayo de Van Vliet de los modelos 1 y 2"
      expect: "energía frente a G_f·A ± 3 %"
      tol_rel: 3.0e-2
    - name: petersson_beam
      setup: "viga de Petersson (flexión en tres puntos, modo I), tres mallas, como los modelos 1 y 2"
      expect: "carga máxima entre mallas ± 3 %; D(8 mm)/(G_f·A) igual entre mallas ± 5 %"
      tol_rel: 5.0e-2
    - name: l_panel
      setup: "panel en L de Winkler (modo mixto), tres mallas, como el modelo 2"
      expect: "carga máxima entre mallas ± 3 %; cola y energía: se reportan frente a los modelos 1 y 2 (observación, por el bloqueo de §10)"
      tol_rel: 3.0e-2
    - name: sen_beam
      setup: "viga SEN de Schlangen como el modelo 2"
      expect: "observación, sin criterio: se reporta frente al modelo 2 y al ensayo"

references:
  - "Rots, J. G. (1988). Computational modeling of concrete fracture. Tesis doctoral, TU Delft."
  - "de Borst, R., Crisfield, M. A., Remmers, J. J. C. y Verhoosel, C. V. (2012). Non-linear Finite Element Analysis of Solids and Structures, 2.ª ed., §6.8.1. Wiley."
  - "de Borst, R. y Nauta, P. (1985). Non-orthogonal cracks in a smeared finite element model. Engineering Computations 2(1), 35–46."
  - "Rots, J. G., Nauta, P., Kusters, G. M. A. y Blaauwendraad, J. (1985). Smeared crack approach and fracture localization in concrete. HERON 30(1)."
  - "Kostense, N. W., Yang, Y., Hendriks, M. A. N. y Rots, J. G. (2025). Smeared cracking and discrete crack dilatancy applied to tension-shear model problems. FraMCoS-12. DOI 10.21012/fc12.1384"
  - "Jirásek, M. y Zimmermann, T. (1998). Rotating crack model with transition to scalar damage. J. Eng. Mech. 124(3), 277–284. DOI 10.1061/(ASCE)0733-9399(1998)124:3(277)"
  - "Feenstra, P. H. (1993). Computational aspects of biaxial stress in plain and reinforced concrete. Tesis doctoral, TU Delft (compresión, por el modelo 2)."
  - "Faria, R., Oliver, J. y Cervera, M. (1998). A strain-based plastic viscous-damage model for massive concrete structures. IJSS 35(14), 1533–1558 (medida de compresión, por el modelo 2)."
  - "DIANA FEA, Theory Manual 10.4, §46.5 Total Strain Crack Model (retención de cortante «damage based», §46.5.4.1 y §46.5.7)."
  - "Červenka, V., Jendele, L. y Červenka, J. (2025). ATENA Program Documentation — Theory, ec. 2.47 (factor de cortante s_F)."
```

---

## Implementación

Bloque 1 (tracción y cortante, sin aplastamiento), 2026-11-24.

- Archivo: `solidum/materials/fixed_smeared_crack_2d.py` (núcleo Numba por punto y adaptador por lotes, como los modelos 1 y 2; reutiliza `_softening.py` y `CrackBandMixin`).
- Clase: `FixedSmearedCrack2D` (programa principal; `f_c`, `G_c` y `beta` se rechazan con un mensaje hasta el bloque 2).
- Tests: `tests/test_fixed_smeared_crack.py` (aceptación `material_point` sin `compression_like_model2`, e `impl_ex`; tangente global ensamblada, implícita e IMPL-EX) y `tests/validation/test_crack_band_objectivity.py` (`bar_in_tension`: Quad4, Tri3, Tri6, Quad8, Quad9; ν = 0.2 en esfuerzo plano y deformación plana; IMPL-EX con las dos leyes).
- Notas de traducción:
  - Ejes: el código usa `T = N⁻¹` (`ε_loc = T·ε`, `σ = Tᵀ·σ_loc`); las pruebas usan la `N` de §1, escrita aparte.
  - Problema local (§8), anidado, con un problema escalar monótono en cada nivel: dentro, `w_n` dada `w_t` (secante en forma cerrada; envolvente con Newton desde la derecha; el tope en forma cerrada); fuera, `w_t` como raíz de la función reducida `R(w_t)`, que es creciente (minimización parcial de una función estrictamente convexa), acotada entre la raíz con `w_n` congelada en `w_n*(0)` y la raíz con `w_n = 0`, con Newton protegido por bisección. La rama de cada dirección se decide una vez en los quiebres y no se reclasifica después. Equivale al «Newton de 2×2 protegido» de §8 (misma solución, única bajo el límite de `h`).
  - Tolerancia en los quiebres (ADR 0006): carga si `b > (1/c + a)·w_max + tol`; cerrada si `b ≤ tol` (`f_t + tol` sin grieta). Sin ella, al empezar cada paso el estado confirmado (`w = w_max`) caía en la carga o en la descarga según el redondeo; en la barra de Quad8 con 20 columnas eso mezclaba tangentes de signo contrario, rompía la simetría del predictor y Newton se iba a la rama localizada dentro del elemento y ciclaba (medido; Diálogo del 2026-11-24). La historia crece sólo en la carga sobre la envolvente. Efecto sobre la respuesta (verificación del bloque 1, O1): en la banda de tolerancia la ley es la secante a historia fija, que deriva de χ (sin disipación espuria); el error de esfuerzo es menor que `tol` en las direcciones con grieta y, en el inicio de una dirección sin grieta, está acotado por `tol·a/(a − |f'(0)|)` ($a = C_{11}/h$), ≈ `tol` con mallas usuales y hasta 4.3·`tol` medido con $h = 0.95\,h_\text{lím}$.
  - IMPL-EX, ley secante (§9): se resuelve con el mismo algoritmo anidado (lineal a trozos), no enumerando las cuatro combinaciones; coincide con la forma cerrada a 6·10⁻¹⁴ (verificación del bloque 1).
  - Cobertura de `bar_in_tension`: los cinco sólidos con ν = 0 en esfuerzo plano y con ν = 0.2 en deformación plana, Quad4 con ν = 0.2 en esfuerzo plano, IMPL-EX con Quad4 (exponencial y Hordijk en esfuerzo plano; exponencial con ν = 0.2 en deformación plana). La barra está alineada con la malla (θ = 0); la banda inclinada respecto de la malla queda para la escalera de validación.
  - IMPL-EX: la salvaguarda `w̃_max ≥ w̃_max^(n)` de §9 se aplica sobre `ω` (`damage_n`, `damage_t` guardan el `ω̃` del esfuerzo), monótona también en el redondeo; `c̃ = c(ω̃)`. `w_n`, `w_t`, `slip` son los de la ley secante del paso (los del esfuerzo); fuera de un paso el estado no cambia.
  - `IS_SYMMETRIC`: falso en la clase (tangente implícita) y verdadero en la instancia con IMPL-EX (secante, simetrizada).
  - `σ_zz` de deformación plana: `ν·(σ_xx + σ_yy)` (concreto entre grietas elástico en el bloque 1).
  - `kappa_c`, `damage_c`, `kappa_c_prev` quedan en 0 sin aplastamiento.

Bloque 2 (aplastamiento), 2026-11-24.

- Archivos: el mismo material (`_fsc_crush_point`, `_implicit_point_crush`, `_frozen_point_crush` y el problema acoplado `_coupled_*`); la validación de `f_c`, `G_c`, `beta` y los parámetros del núcleo, compartidos con el modelo 2 (`parse_compression`, `compression_kernel_params` en `_tension_compression.py`); el corrector abandona un intento con residuo no finito (`solidum/math/solvers/corrector.py`, `nonfinite_residual`).
- Tests: `tests/test_fixed_smeared_crack_crushing.py` (aceptación del bloque 2 en el punto material), `tests/test_corrector_nonfinite_residual.py` (residuo no finito → bisección; falla sin el cambio del corrector, comprobado) y en `tests/validation/test_crack_band_objectivity.py` las barras en compresión del modelo 3 (`bar_in_compression` con Quad4, Tri3, Tri6, y Quad8 y Quad9 con el factor; igual al modelo 2 paso a paso, implícita e IMPL-EX).
- Notas de traducción:
  - Sin `f_c` el núcleo es el del bloque 1 sin tocar (`_fsc_b1_point`): todos sus resultados quedan iguales bit a bit.
  - `κ⁻` empieza en $f_c/(3E)$, como en el modelo 2 (IMPL-EX no crea daño de compresión en un punto que la implícita deja elástico).
  - Predictor de §8b: el bloque 1 tal cual, y su esfuerzo es el efectivo del concreto entre grietas (elástico), con el que se evalúa $d^-(\kappa^-)$.
  - Problema acoplado: incógnitas $x = [w_n, w_t, s]$ siempre (una dirección rígida de IMPL-EX con $\tilde c = 0$ tiene la fila $k\,w_i = 0$). La función mínimo usa $k = C_{11}/h$. El valor de $g_i$ es continuo (secante hasta $w_\text{max}$, envolvente después); la pendiente y el crecimiento de la historia usan la banda $w_\text{max} + \delta$, $\delta = $ `tol`$/(1/c + k)$, la misma que el bloque 1, así que dentro de la banda la dirección descarga.
  - Convergencia con max$|r_i| \le$ `tol` (la de admisibilidad, en esfuerzo) y después un paso más de pulido, que se acepta si no aumenta el residuo: deja la solución cerca del redondeo, como el bloque 1. Sin él, la tolerancia contaminaba las diferencias finitas de la tangente (1.2·10⁻⁵ medido).
  - Carga de $\kappa^-$ en el problema acoplado con tolerancia: $\varepsilon^-_\text{eq} > \kappa^-_n + $ `tol_c`$/E$ (la de admisibilidad del ADR 0006, como la banda $\delta$). Al volver a evaluar el estado confirmado, el Newton acoplado reproduce $\varepsilon^-_\text{eq}$ sólo hasta el redondeo; sin la tolerancia $\kappa^-$ «crecía» y la tangente era la de carga en el 5–7 % de esos pasos (H1 de la verificación; el mismo tipo de defecto que la barra de Quad8 del bloque 1). Sin grieta no hace falta ($\varepsilon^-_\text{eq}$ se recalcula idéntico) y no se usa, para que sea el modelo 2 bit a bit.
  - Corrector: el intento se abandona si el residuo completo no es finito, no sólo $\lVert\mathbf T^\top\mathbf R\rVert$: el punto que falla puede estar en un elemento con todos sus grados de libertad prescritos (H2 de la verificación).
  - Tangente: $\mathbf C_\text{tan} = \mathbf C_b + \mathbf C_b\mathbf N\,\mathbf J^{-1}\mathbf R_\varepsilon/h$ por la función implícita sobre $r(x, \varepsilon) = 0$, con el jacobiano $\mathbf J$ del Newton en la solución; equivale a la fórmula de §8b. Sistemas de 3×3 por eliminación con pivote parcial (sin excepciones dentro de Numba).
  - IMPL-EX con aplastamiento: `damage_c` guarda el $\tilde d^-$ del esfuerzo, como `damage_n`, `damage_t` guardan $\tilde\omega$; `IS_SYMMETRIC` es verdadero en la instancia sólo con IMPL-EX y sin aplastamiento.
  - $\sigma_{zz}$ de deformación plana con $d^- > 0$: $\bar{\boldsymbol\sigma}$ se recupera de $\boldsymbol\sigma$ y de $d^-$ (mismas direcciones principales), como en el modelo 2.
- Medido al programar: la compresión sin grieta coincide con el modelo 2 a 5·10⁻¹⁶; el aplastamiento paralelo a una grieta abierta llega a $f_c$ con la apertura fija; la tangente coincide con diferencias finitas a ≤ 6·10⁻⁸ en 1674 estados con aplastamiento (implícita e IMPL-EX, esfuerzo plano y deformación plana). Con IMPL-EX y 1000 pasos la barra en compresión (igual en el modelo 2) sobrepasa el pico 0.15–0.36 % y, con Tri3 y 5 columnas, disipa 0.63·$G_c A$ desde el pico: error del método en la estructura, para la validación.

---

## Diálogo

- **2026-11-24** · **Inicio del modelo 3** (usuario: «empecemos con el modelo 3 de la familia de concreto»). Búsqueda de `bibliotecario-solidum`: `informes/bibliotecario/busqueda_modelo3_2026-11-24.md`. Admisión: la fija cumple con de Borst et al. 2012 §6.8.1 y tres códigos (DIANA, Abaqus, ANSYS); la rotante, con el mismo texto y sólo DIANA. Dato dicho al usuario: en la viga de Arrea–Ingraffea (Rots 1988, Fig. 6.6, lectura gráfica) la fija con β = 0.05 se queda en el 90 % del pico, la fija con β variable no ablanda, la rotante y la fija con β ≈ 0 se estancan en el 45 % y el 25 % del pico; «all smeared crack results are too stiff in the post-peak regime» (p. 86). **Usuario: «si, vamos con la opción 1»**: seguir con el objetivo redefinido del encabezado.
- **2026-11-24** · **Formulación (IA):** descomposición de la deformación en forma secante (§1), con sus razones. Ley de ablandamiento por omisión: la exponencial, como en los modelos 1 y 2; Hordijk opcional.
- **Q1. Dirección de la grieta. Decidida: sólo grieta fija** (usuario, 2026-11-24, «si, vamos con la opción 1»), con segunda grieta perpendicular. Descartadas: fija y rotante a elegir (la rotante se parece a lo que ya hace el modelo 2) y sólo rotante (sin cortante en la grieta, más sensible a la orientación de la malla, un solo código).
- **Q2. Ley de retención de cortante. Decidida: ligada a la apertura, sin parámetro** (usuario, 2026-11-24, «si, vamos con la opción 1»): $\tau/s = \sigma_{nn}(w)/(2w)$, equivalente exacta de la «por daño» de DIANA y objetiva por construcción (§3). Descartadas: la misma forma con factor ajustable $s_F$ (ATENA; se volvería un parámetro de calibración; se puede añadir después sin cambiar lo demás), fracción fija $\beta$ (Rots, ANSYS, Abaqus; no objetiva y, con $\beta = 0.05$, la viga de Arrea–Ingraffea no ablanda, Rots 1988 p. 86) y caída a cero en $w = d_\text{agg}/2$ (DIANA; la que más sobrestima el cortante en Kostense et al. 2025).
- **Q3. Compresión. Decidida: aplastamiento opcional con las piezas del modelo 2** (usuario, 2026-11-24, «si, vamos con la opción 1»): con $f_c$ y $G_c$, Drucker-Prager sobre la parte negativa del esfuerzo efectivo del concreto entre grietas y parábola de Feenstra con $G_c/h$ (§4); la compresión paralela a una grieta abierta se conserva completa. Descartadas: añadir la reducción de Vecchio–Collins por agrietamiento transversal (usa la deformación $w/h$ de la banda: dependería de la malla; se puede añadir después con $w$ y la separación entre grietas) y compresión siempre elástica (una viga reforzada que falle por aplastamiento daría un resultado falso). Se programa en un segundo bloque, como en el modelo 2.
- **2026-11-24** · **Formulación numérica, variables internas y aceptación redactadas por la IA** (§6–§10 y contrato). Plan: bloque 1, tracción y cortante (§1–§3, §6–§9), implícita e IMPL-EX, pruebas del punto material y barra en tracción; bloque 2, aplastamiento (§4) y barra en compresión; después, la escalera de validación (Van Vliet, Petersson, panel en L, SEN como observación) y `status: validated`. 3D en una spec aparte.
- **Q4. Aceptación. Decidida: aprobada como está** (usuario, 2026-11-24, «si, vamos con la opción 1»): objetividad estricta en el punto material, las barras, Van Vliet y Petersson; en el panel en L sólo la carga máxima es criterio, y la cola y la energía se reportan como observación (bloqueo de §10); la viga SEN, sólo observación. Descartado: hacer de la cola del panel en L un criterio (± 10 % entre mallas), que lo más probable es que no se cumpla por una causa cinemática.
- **2026-11-24** · **Revisión de la spec por `teorico-solidum`** (`informes/teorico/verificacion_teorica_2026-11-24_modelo3_spec.md`; guiones en `informes/prototipos_2026-11-24/modelo3_spec/`). Veredicto: el bloque 1 se puede programar tras corregir H1–H5; el bloque 2 no está listo (H4). La formulación no tiene problema de fondo; equivalencia de §3, potencial χ, energía $G_f/h$, tangente y forma cerrada de IMPL-EX comprobadas; ningún hallazgo contradice Q1–Q4. **Aplicado por la IA:** H2 (§8 reescrito: problema normal de dos incógnitas con el cortante en forma cerrada, rama decidida en los quiebres, Newton protegido), H3 (`parallel_compression_open_crack`: σ_tt − ν·σ_nn = E·ε_t; el texto de §4 matizado), H4 (pendiente del bloque 2 anotado en §4), H5 (matriz N y θ explícitos en §1, columnas quitadas en la tangente, tope antes de w/f), H6 (precisión de la equivalencia con DIANA), H7 (energía en modo mixto (1 + λ²/2)·G_f en §7), H8 (tolerancias y seis pruebas nuevas en el contrato) y H9 (cifras de la diferencia entre las dos «exponential» en §2). **H1, pendiente del usuario** (umbral del `ValueError`, Reglas §3).
- **H1. Límite del tamaño del elemento. Decidida: bajar el límite** (usuario, 2026-11-24, «si, vamos con la opción 1»): `ValueError` con $h \ge \ell_{ch}/(1+\nu)$ (exponencial) o $0.738\,\ell_{ch}/(1+\nu)$ (Hordijk), la condición necesaria y suficiente de unicidad del problema local con dos direcciones en carga (§2, §8). Descartado: conservar el límite de una dirección y elegir con una regla fija entre las hasta cinco soluciones medidas por el teórico.
- **2026-11-24** · **Spec lista para programar el bloque 1** (tracción y cortante, implícita e IMPL-EX). El bloque 2 (aplastamiento) espera la redacción de H4 y otra revisión del teórico.
- **2026-11-24** · **Bloque 1 programado** (sección «Implementación»). Pruebas de aceptación `material_point` (salvo `compression_like_model2`, del bloque 2) e `impl_ex` verdes, y `bar_in_tension` con los cinco sólidos 2D: energía G_f·A (Quad8 y Quad9, G_f·A/α como los modelos 1 y 2) también con ν = 0.2 en deformación plana, sin el factor 1 − ν² del modelo 1 (*medido*). Caveat de §10 registrado en `fixed_direction_and_inclined_stress`: con la grieta recién pasada el pico (w_max = 0.02·G_f/f_t) y la carga girada 60°, el esfuerzo principal llega a 1.36·f_t a unos 45° de la grieta sin que se forme otra (30°: 1.006·f_t) (*medido*). **Hallazgo al programar** (*medido*): en la barra de Quad8 con 20 columnas Newton ciclaba en el ablandamiento. En esa barra la solución uniforme dentro del elemento es un equilibrio inestable con los dos modelos (la tangente global tiene 5 a 9 valores propios negativos); el modelo 1 la conserva porque compara la misma magnitud calculada al decidir carga o descarga. El modelo 3 resolvía de nuevo el estado confirmado y la decisión en el quiebre `w = w_max` dependía del redondeo: al empezar cada paso entre 1 y 16 de los 36 puntos de la banda recibían la tangente de ablandamiento y los demás la secante, y el predictor dejaba de ser simétrico. Corregido con la tolerancia de admisibilidad del ADR 0006 en los quiebres (sin efecto físico: desplaza la decisión en `tol` ≈ 10⁻¹⁰·f_t); con una imperfección geométrica de 10⁻⁹ la barra de Quad8 ya se comporta como con el modelo 1. Lo protegen la prueba `test_estado_confirmado_vuelto_a_evaluar_no_carga` y las barras de Quad8, Quad9 y Tri6.
- **2026-11-24** · **Verificación del bloque 1 por `teorico-solidum`: sí al commit** (`informes/teorico/verificacion_teorica_2026-11-24_modelo3_bloque1.md`; guiones en `informes/prototipos_2026-11-24/modelo3_bloque1/teorico/`). Con su prototipo anterior al código (enumeración de ramas con Newton de 2×2, otro método que el anidado): σ a 5·10⁻¹³ en 4783 estados al azar (dos leyes, dos hipótesis, una y dos grietas, carga, descarga, cierre, tope, inicio, h hasta 0.95 del límite), tangente a 6.5·10⁻¹⁵ de su fórmula y a 6.4·10⁻¹⁰ de diferencias finitas, IMPL-EX y kernel por lotes al redondeo; energía hasta el tope a ≤ 2.2·10⁻⁸ de D_tope, modo mixto (1 + λ²/2)·G_f con las dos leyes, ciclos con historia fija → 0 como Δ². Las dos decisiones de implementación: (a) el problema anidado da la misma solución y es robusto hasta 0.999999 del límite de h (media de 5 iteraciones fuera); (b) la tolerancia en los quiebres no cambia la física más allá de unas pocas veces `tol`, es Kuhn–Tucker con tolerancia y respeta χ. Nada que decida el usuario. **Aplicado por la IA:** O1 (cota del efecto de la tolerancia, en «Implementación»), O2 (valor propio mínimo de la secante de IMPL-EX con el tope, §9), O3 (prueba del intercambio n↔t, que recorre el algoritmo por la otra ruta), O4 (barra con ν = 0.2 en deformación plana en los cinco sólidos, IMPL-EX también en deformación plana, energía de la implícita con Hordijk; cobertura anotada en «Implementación»), O5 (modo mixto también con Hordijk). Campaña de mutaciones: 20/20 detectadas (`informes/prototipos_2026-11-24/modelo3_bloque1/mutaciones_resultados.md`).
- **2026-11-24** · **Bloque 2 redactado por la IA (H4)**, con el plan aprobado por el usuario («si»): §2 (el aplastamiento no cambia el inicio de la primera grieta), §4 (evolución de $d^-$ dentro del problema local, aplastamiento paralelo a una grieta abierta con pico $f_c$, sin acoplamiento entre los daños, lo que cambia en la formulación numérica), §8b (problema local con complementariedad, predictor que reproduce el bloque 1 sin aplastamiento, Newton semiliso con búsqueda lineal, tangente, inicio), §9 (IMPL-EX sin forma cerrada y con tangente no simétrica con $\tilde d^- > 0$), §10 y seis pruebas nuevas en el contrato. Pendiente: revisión de `teorico-solidum` antes de programar y, después, Q5 al usuario.
- **2026-11-24** · **Revisión del bloque 2 por `teorico-solidum`: sí, se puede programar** tras aplicar H1 y H4 (`informes/teorico/verificacion_teorica_2026-11-24_modelo3_bloque2_spec.md`; guiones en `informes/prototipos_2026-11-24/modelo3_bloque2_spec/`). Formulación de §8b, predictor (exacto), tangente (C_b a 1.4·10⁻⁸ de diferencias finitas; sin el término de σ̄_zz, error 1.8·max|C|) e IMPL-EX correctos; ningún hallazgo contradice Q1–Q4 ni H1. En 1.9 millones de pasos de caminos al azar, ningún estado con varias soluciones. **Aplicado por la IA:** H1 (cascada del Newton local: paso completo, búsqueda lineal, arranque en e = 0, bisección), H2 (residuo de cortante en unidades de esfuerzo), H3 (dos mecanismos de pérdida de unicidad y frecuencias medidas, §4 y §8b), H4 (contrato: tolerancias y planteamientos de cinco pruebas, cinco pruebas nuevas), H5 (umbral exacto d* de §9), H6 (rectificación de «lineal a trozos» también en §1 del modelo 2) y H7 (pico en deformación plana, §4). Sólo Q5 es del usuario.
- **Q5. Concreto aplastado con una grieta inclinada: qué hace Solidum si el problema local puede tener varias soluciones. Decidida: la propuesta de la IA** (usuario, 2026-11-24, «si»; también la recomendación del teórico): sin límite ni aviso nuevos; solución por continuidad con el paso confirmado (§8b); si la cascada del Newton local no converge, el paso se repite más corto y queda un aviso único; el validador informa cuántas veces se activa ese aviso en la escalera de validación.
- **2026-11-24** · **Bloque 2 programado** (sección «Implementación»). Pruebas de aceptación del bloque 2 verdes (84 en el punto material, barras en compresión). Pieza genérica nueva en el programa principal: el corrector abandona un intento con residuo no finito y el solver repite el paso, con su motivo (antes se habría registrado como «tangente singular»). Hallazgos al programar (*medido*): sin el paso de pulido del Newton acoplado, la tolerancia de convergencia contaminaba las diferencias finitas de la tangente (1.2·10⁻⁵); con IMPL-EX y 1000 pasos la barra en compresión del modelo 2, y por tanto la del modelo 3, sobrepasa el pico 0.15–0.36 % y con Tri3 en la malla más gruesa disipa 0.63·$G_c A$ desde el pico (queda para la validación).
- **2026-11-24** · **Verificación del bloque 2 por `teorico-solidum`: sí al commit** (`informes/teorico/verificacion_teorica_2026-11-24_modelo3_bloque2.md`; guiones en `informes/prototipos_2026-11-24/modelo3_bloque2/teorico/`). Con su prototipo anterior al código: σ, estado y tangente al redondeo en unos 480 000 pasos (σ ≤ 3.2·10⁻¹³·$f_t$, tangente ≤ 6.5·10⁻¹² de su fórmula), sin NaN; sin `f_c`, el bloque 1 bit a bit en 28 800 pasos; el modelo 2 bit a bit tras el refactor; la banda $\delta$ y el paso de pulido no cambian la física más allá de la tolerancia (0.77·`tol`). Mutaciones: 34 de 37 válidas detectadas, todas las de física y tangente. Nada que decida el usuario. **Aplicado por la IA:** H1 (tolerancia de la carga de $\kappa^-$ en el problema acoplado), H2 (el corrector mira el residuo completo), H3 (pruebas `confirmed_state_reevaluated_coupled` y de caminos que necesitan la cascada; comprobado que fallan sin la banda $\delta$, sin la tolerancia de $\kappa^-$ y sin la cascada) y H4 (tolerancias al nivel medido, prueba a nivel de elemento, contrato al día). Con $d^- > 0$ el límite de $h$ ya no garantiza una sola solución (§4, §8b); ningún límite de $h$ lo arregla cuando $d^- \to 1$. Propuesta de la IA: sin límite ni aviso nuevos; la solución se elige por continuidad (Newton desde el paso confirmado); si el Newton local no converge, el paso se repite más corto y queda un aviso; la validación cuenta cuántas veces pasa (`local_problem_robustness_crushing`). Alternativas: avisar en cada punto donde la unicidad no está garantizada (muchos avisos sin acción posible para el usuario) o impedir que una grieta se cargue sobre concreto muy aplastado (cambia la física).
- **2026-11-25** · **Escalera de validación, primera tanda (Van Vliet y Petersson) con `validador-solidum`** (usuario: «si» a dos tandas; `informes/reportes/reporte_modelo3_vanvliet_petersson_2026-11-24.pdf`). Van Vliet **a medias**: objetivo entre mallas (pico −0.93 %, energía −1.44 %), pero $D(1.2\text{ mm}) = 1.14\,G_f A$ y subiendo; la malla fina se detuvo en 0.62 mm (Newton oscila entre dos estados). Petersson **no validado**: pico 868.6 / 928.4 / 908.0 N (6.6 % entre mallas; modelos 1 y 2: 746–750 N), $D(8) = 2.43$–$2.91\,G_F A$. 0 avisos de residuo no finito (Q5). La segunda tanda (panel en L, SEN) no se corrió.
- **2026-11-25** · **Diagnóstico** (usuario: «si»; `informes/modelo3_petersson_yaml/diagnostico/`): Quad4 y Tri3, implícita e IMPL-EX: $\theta$ sigue la dirección principal (mediana 0.7–1.2°), mismos ángulos (|θ| medio 29–33°) y mismo pico con las dos integraciones; Quad4 con implícita no converge tras el pico. **Estudio de escritorio** (usuario: «ok»): `teorico-solidum` (causa: estado casi equibiaxial al nacer, §10) y `bibliotecario-solidum` (Rots et al. 1985 describen el fenómeno en esta viga; DIANA, Abaqus y ANSYS fijan la dirección igual que Solidum). **Reproducción de Slobbe 2015** (usuario: «ok»; `informes/modelo3_slobbe_yaml/`): el modelo 1 coincide con Slobbe (+0.2 %); el modelo 3 no (+11.5 % con Quad8 3×3, +39 % a 0.3 mm). **Separación de causas** por `teorico-solidum` (usuario: «si»): la deformación total de Slobbe es, en esfuerzo plano, la formulación de esta spec (correcciones en §1 y §3); lo que lo separa es su ancho de banda empírico $l_p(\theta)$, una compensación, que no cierra toda la brecha. Dato dicho una vez al usuario, sin recomendarlo porque reabre Q1: DIANA ofrece «rotating to fixed».
- **Q6. Qué hacer con la sobrestimación en modo I. Decidida: limitación documentada** (usuario, 2026-11-25, «si»; recomendación de la IA y del teórico): la causa medida y publicada en §10 y en el catálogo, correcciones de §1 y §3, `status: implemented` (no `validated`), sin segunda tanda de validación y sin código nuevo. Descartadas: adoptar el ancho $l_p(\theta)$ de Slobbe (compensación empírica que no cierra la brecha y probablemente no está en un código general, *hipótesis*) y sólo la proyección de Govindjee (a 0.3 mm la brecha pasa de +42 % a +28 %).
