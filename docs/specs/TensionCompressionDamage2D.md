# TENSION COMPRESSION DAMAGE 2D — daño de tracción y de compresión con cierre de grieta, banda de grieta e IMPL-EX

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> **Componente nuevo** (modelo 2 de la familia de modelos de concreto): dos daños escalares, $d^+$ de tracción y $d^-$ de compresión, con **efecto unilateral**: una grieta cerrada recupera la rigidez en compresión. Cambia la ley de esfuerzo respecto al daño isótropo, así que no es una variante de [IsotropicDamageCrackBand2D](IsotropicDamageCrackBand2D.md); reutiliza sus piezas: banda de grieta $\alpha(h)$ (§2 de esa spec), medida de Rankine (§5), tope del daño (§3), siembra de $h$ en los sólidos (§4) e integración IMPL-EX con el pseudo-tiempo de los solvers ([ADR 0034](../adr/0034-pseudo-tiempo-de-los-solvers-incrementales.md)).
>
> **Motivación (diagnóstico de la viga SEN, 2026-11-11).** Con `IsotropicDamageCrackBand2D` la viga SEN de Schlangen no es objetiva con mallas finas: la banda baja casi vertical desde la entalla y disipa 1.5–2 veces $G_f$. El diagnóstico descartó el elemento (`Tri6`) y la orientación de la malla (malla no estructurada) y lo atribuyó al daño isótropo: en la zona de cortante el elemento dañado pierde también la rigidez de la biela comprimida que cruza la banda. Un prototipo que daña sólo la parte de tracción dio, con las mismas mallas, la grieta diagonal, la misma energía en $h/4$ y $h/8$ ($0.85\,G_f t L$) y diferencias entre mallas de 1.2 kN.
>
> Uso previsto: estructuras de concreto con agrietamiento por tracción y, si se dan $f_c$ y $G_c$, aplastamiento en compresión, en 2D (esta spec) y 3D ([TensionCompressionDamage3D](TensionCompressionDamage3D.md)).
>
> **Estado: borrador para revisión del usuario. Q1 (ley (a) de Faria, Oliver y Cervera), Q2 (aplastamiento opcional), Q3 (sin deformaciones plásticas), Q4 (sólo Rankine) y Q6 (aceptación) decididas el 2026-11-17; el límite de DIANA en compresión, no adoptado. Q5 (clasificación) decidida el 2026-11-18: programa principal, con el criterio (i) cumplido por la familia y no por la forma. Spec completa.**

---

## Especificación física

### 0. Descripción general

Daño escalar doble en pequeñas deformaciones para sólidos 2D (esfuerzo plano o deformación plana): $d^+ \in [0, 1)$ degrada la parte de tracción del esfuerzo y $d^- \in [0, 1)$ la de compresión. Sin deformaciones plásticas: la descarga es secante hacia el origen. Cada daño ablanda con su propia energía regularizada por banda de grieta: $G_f/h$ en tracción y, si se activa la compresión, $G_c/h$ en la rama descendente de compresión.

`STRAIN_DIM = 3`, Voigt del proyecto $[\varepsilon_{xx}, \varepsilon_{yy}, \gamma_{xy}]$ con $\gamma_{xy} = 2\varepsilon_{xy}$. Convención: $\sigma > 0$ ⇔ tracción.

### 1. Ley de esfuerzo (Q1, decidida por el usuario el 2026-11-17: (a))

Partición espectral del esfuerzo efectivo de Faria, Oliver y Cervera (1998), la misma partición que usa OpenSees `ASDConcrete3D` (sus criterios de daño son otros: superficies con $I_1$, $\sqrt{J_2}$ y $\sigma_\text{max}$ de Petracca et al. 2022):

$$\boxed{\;\boldsymbol\sigma = (1-d^+)\,\bar{\boldsymbol\sigma}^+ + (1-d^-)\,\bar{\boldsymbol\sigma}^- = (1-d^-)\,\bar{\boldsymbol\sigma} - (d^+-d^-)\,\bar{\boldsymbol\sigma}^+\;}$$

$$\bar{\boldsymbol\sigma} = \mathbf C_e:\boldsymbol\varepsilon, \qquad \bar{\boldsymbol\sigma}^+ = \sum_i\langle\bar\sigma_i\rangle_+\,\mathbf n_i\otimes\mathbf n_i, \qquad \bar{\boldsymbol\sigma}^- = \bar{\boldsymbol\sigma} - \bar{\boldsymbol\sigma}^+,$$

con $\bar\sigma_i$, $\mathbf n_i$ los esfuerzos y las direcciones principales del esfuerzo efectivo **del estado 3D asociado**:

- **Esfuerzo plano**: $\bar\sigma_{zz} = 0$, es decir $\bar{\boldsymbol\sigma} = \mathbf C_\text{ep}:\boldsymbol\varepsilon$ con la elasticidad de esfuerzo plano y $\varepsilon_{zz} = -\tfrac{\nu}{1-\nu}(\varepsilon_{xx}+\varepsilon_{yy})$; la partición es la de 2×2 del plano y el esfuerzo nominal cumple también $\sigma_{zz} = 0$.
- **Deformación plana**: $\varepsilon_{zz} = 0$; $\bar\sigma_{zz} = \lambda(\varepsilon_{xx}+\varepsilon_{yy})$ es un esfuerzo principal (eje $z$) y entra en la partición, de modo que $\sigma_{zz} = (1-d^\pm)\,\bar\sigma_{zz}$ según su signo.

Los dos casos son la restricción exacta de la ley 3D (§11).

**Comportamiento uniaxial** (calculado, §11):

- Tracción (esfuerzo plano y 3D): $\sigma = (1-d^+)\,E\,\varepsilon$ **exacto**, igual que el daño isótropo, así que la calibración de la banda de grieta del modelo 1 sirve tal cual. La contracción lateral es $\varepsilon_\perp = -\nu\,\varepsilon$: la deformación de apertura de la grieta también contrae lateralmente, rasgo de los modelos formulados en el esfuerzo efectivo. En deformación plana con el lateral libre, $\sigma_{xx} = (1-d^+)\,E\,\varepsilon/(1-\nu^2)$, $\sigma_{zz} = \nu\,\sigma_{xx}$, y la energía disipada es $(1-\nu^2)\,G_f/h$ con el pico en $f_t$: el mismo rasgo del modelo 1 con Rankine, documentado allí sin corregir (2026-11-11).
- Compresión con la grieta cerrada ($d^+ \to 1$, $d^- = 0$): rigidez $E$ y Poisson $\nu$, **recuperación completa**.

**Energía.** Con los daños fijos, $\boldsymbol\sigma(\boldsymbol\varepsilon)$ es continua, lineal a trozos y positivamente homogénea de grado 1: la descarga es secante hacia el origen y devuelve $\tfrac12\boldsymbol\sigma:\boldsymbol\varepsilon$. Con $\nu = 0$ ($\bar{\boldsymbol\sigma} = E\,\boldsymbol\varepsilon$) la ley deriva del potencial $\psi = \tfrac E2\bigl[(1-d^+)\,\boldsymbol\varepsilon^+:\boldsymbol\varepsilon^+ + (1-d^-)\,\boldsymbol\varepsilon^-:\boldsymbol\varepsilon^-\bigr]$. **Con $\nu \ne 0$ no deriva de ningún potencial.** Faria, Oliver y Cervera la presentan como derivada de $\psi = (1-d^+)\psi_0^+ + (1-d^-)\psi_0^-$, $\psi_0^\pm = \tfrac12\bar{\boldsymbol\sigma}^\pm:\mathbf C_e^{-1}:\bar{\boldsymbol\sigma}$ (sus ecs. 3–41), pero el paso de su ec. (38) a la (39) usa $(\partial\bar{\boldsymbol\sigma}^+/\partial\boldsymbol\varepsilon):\boldsymbol\varepsilon = \bar{\boldsymbol\sigma}^+$ con la simetría mayor de $\partial\bar{\boldsymbol\sigma}^+/\partial\boldsymbol\varepsilon$, que sólo se cumple con $\nu = 0$; su Remark 2.7 admite que la tangente no es simétrica (leído en la fuente; calculado: $\partial\psi/\partial\boldsymbol\varepsilon$ difiere de $\boldsymbol\sigma$ hasta 135 % con $\nu = 0.2$). Carol y Willam (1996) tratan el problema en la clase general de los modelos de recuperación de rigidez. La tangente con daño fijo no es simétrica y en un ciclo cerrado de deformación con los daños fijos y direcciones principales que giran el trabajo neto no es nulo, de los dos signos. Medido con $\nu = 0.2$ (400 ciclos al azar por caso, $d^+ \in \{0.5, 0.9\}$, $d^- \in \{0, 0.3\}$, 3D y esfuerzo plano): mediana de 2–6 % de la energía elástica máxima del ciclo, percentil 90 de 6–19 %, máximo de 15–48 %, y la mitad de los ciclos crea energía. En carga monótona el efecto es pequeño: residuo del balance de energía de 0.6–0.9 % en la viga SEN (tabla de abajo). Es un caveat documentado: el modelo no está pensado para ciclos sobre material dañado (no tiene deformaciones plásticas, §5), y la prueba `closed_cycle_work` lo mide. Los daños nunca decrecen (§2, §3) y en tracción uniaxial monótona (esfuerzo plano y 3D) el material disipa exactamente $G_f/h$ por unidad de volumen.

**Compresión paralela a una grieta abierta** (calculado por el teórico, `informes/verificacion_teorica_2026-11-17_modelo2_spec.md`, H1). Con $d^+ \to 1$, el esfuerzo paralelo a la grieta es $\sigma_{tt} = \tfrac{E}{1-\nu^2}(\varepsilon_t + \nu\,\varepsilon_n)$ mientras sea negativo: la deformación de apertura $\varepsilon_n$ entra por Poisson en el esfuerzo efectivo y anula la compresión paralela cuando $\varepsilon_n \ge |\varepsilon_t|/\nu$ (esfuerzo plano; $(1-\nu)|\varepsilon_t|/\nu$ en deformación plana). Como $\varepsilon_n = w/h$, el umbral depende de la malla: con una apertura $w = 0.05$ mm, $h = 10$ mm borra una compresión paralela de hasta $|\varepsilon_t| = 10^{-3}$ (≈ 32 MPa) y $h = 100$ mm sólo hasta $10^{-4}$ (≈ 3 MPa); además $\bar{\boldsymbol\sigma}^-$ se reduce y $d^-$ se retrasa en esas zonas. Es el mismo mecanismo que evita el bloqueo de esfuerzos en la viga SEN, visto del otro lado: una biela comprimida paralela a grietas abiertas pierde capacidad, más cuanto más fina la malla. Caveat documentado; la prueba `parallel_compression_open_crack` lo registra.

#### Alternativas estudiadas (Q1)

Cálculo con daño fijo, 3D, $\nu = 0.2$, $d^+ = 0.9$ (`informes/prototipos_2026-11-11/modelo2/` y `informes/prototipos_2026-11-17/modelo2_ley_a/`):

| Partición | Tracción uniaxial $\sigma/(E\varepsilon)$ (exacto 0.100) | Grieta cerrada en compresión (perfecto 1) | ¿Deriva de un potencial? | Códigos |
|---|---|---|---|---|
| **(a) Esfuerzo efectivo (adoptada)**: $\boldsymbol\sigma = (1-d^+)\bar{\boldsymbol\sigma}^+ + (1-d^-)\bar{\boldsymbol\sigma}^-$, $\bar{\boldsymbol\sigma} = \mathbf C_e:\boldsymbol\varepsilon$ (Faria, Oliver y Cervera 1998) | 0.100 | 1.000 | **No**, con $\nu \ne 0$: tangente con daño fijo asimétrica (hasta 30 % del mayor coeficiente); trabajo en ciclos cerrados con daño fijo de los dos signos (mediana 2–6 %, máximo 15–48 % de la energía elástica; Carol y Willam 1996). Con $\nu = 0$, sí | OpenSees `ASDConcrete3D` |
| (b) Deformación: $\psi = (1-d^+)\psi_0^+(\boldsymbol\varepsilon) + (1-d^-)\psi_0^-(\boldsymbol\varepsilon)$ (Ju 1989) | 0.109 (descalibra la banda) | 0.834 | Sí | Code_Aster `ENDO_ISOT_BETON` (un daño) |
| (c) Complementaria con la traza dañada (Ladevèze y Lemaitre 1984) | 0.100 | 1.000 | Sí, pero **no convexa** con daño alto: con $d^+ = 0.9$ hay dos esfuerzos para la misma deformación | HYPLAS (de Souza Neto et al. 2008, §12.6.1, ec. 12.99: un solo $D$ y $d^- = hD$) |
| (d) Complementaria con Poisson intacto (descartada, ver abajo) | 0.100 | 1.000 | Sí, convexa; ciclos con trabajo nulo | — |

**Viga SEN con (a) y (d)** (malla de FEAP del usuario subdividida, Tri3, IMPL-EX con paso de CMSD de 0.125 µm, carga del ensayo, sólo $d^+$; prototipos en `informes/prototipos_2026-11-11/modelo2/` y `feap_tesis/diagnostico/`; la fila del daño isótropo es la del modelo 1):

| Ley | Malla | $P_\text{max}$ [kN] | $P$ en CMSD 0.05 / 0.10 / 0.20 mm [kN] | $D/(G_f t L)$ | $x_{40}$ [mm] | Residuo del balance |
|---|---|---|---|---|---|---|
| Isótropo (modelo 1) | $h/4$ · $h/8$ | 37.6 · 36.5 | 25.3 / 2.2 / 1.0 · 24.8 / 4.1 / 1.2 | 1.54 · 1.98 | 226 · 225 (vertical) | 0.14 % |
| (a) Faria et al. | $h/4$ · $h/8$ | 39.8 · 39.2 | 17.0 / 3.7 / 1.5 · 15.8 / 3.7 / 1.4 | 0.85 · 0.86 | 248 · 246 (diagonal) | 0.9 % · 0.6 % |
| (d) Complementaria | $h/2$ · $h/4$ · $h/8$ | 40.8 · 39.7 · 39.4 | 22.5 / 8.3 / 2.8 · 25.0 / 12.4 / 7.5 · 16.6 / 8.0 / 5.2 | 0.84 · 1.01 · 0.86 | 249 · 244 · 243 (diagonal) | 0.01 % · 0.01 % · 0.07 % |

($x_{40}$: posición de la banda a media altura; la grieta física va en diagonal hacia la placa de carga, 240–250 mm. Residuo: máximo de $(W_\text{ext} - E_\text{el} - D)/W_\text{ext}$.)

- Las dos forman la grieta diagonal y llegan al mismo pico (39–41 kN), que es lo que el daño isótropo no lograba.
- **(a) es objetiva**: entre $h/4$ y $h/8$ la curva entera difiere a lo sumo 1.2 kN y la energía por unidad de longitud de grieta es la misma. Su balance de energía no cierra exactamente (0.6–0.9 %), por la parte del esfuerzo que no deriva del potencial.
- **(d) cierra el balance de energía** (≤ 0.07 %), pero **no es objetiva** después del pico: entre mallas la carga difiere hasta 8 kN en CMSD = 0.05 mm, la cola no tiende a cero (2.8 a 7.5 kN en 0.2 mm, sin orden con la malla) y aparecen ramas de daño secundarias. Causa: **bloqueo de esfuerzos** (*stress locking*, Rots 1988). En los puntos ya totalmente dañados de $h/4$ el esfuerzo paralelo a la grieta es de −8.4 MPa en promedio (mínimo −33 MPa): la banda inclinada de triángulos de deformación constante no reproduce una apertura pura, aparecen deformaciones paralelas a la grieta espurias y (d), que conserva intacta la rigidez paralela a una grieta abierta (físicamente correcto), las convierte en esfuerzo. En (a) el mismo dato es −0.05 MPa: la apertura de la grieta entra en el esfuerzo efectivo como $\lambda\,\varepsilon_\text{apertura}$ y anula la compresión paralela. **Es justo la parte de (a) que no deriva de un potencial la que evita el bloqueo.**

En resumen: con triángulos lineales, la ley consistente (d) deja ver el bloqueo que (a) disimula. Los remedios del bloqueo son cinemáticos (seguimiento de la grieta, Cervera y Chiumenti 2006; elementos mixtos deformación–desplazamiento, Cervera, Chiumenti y Codina 2010; discontinuidades embebidas) o una transición a daño escalar con la apertura (Jirásek y Zimmermann 1998), y ninguno es objeto de esta spec.

#### Forma complementaria (d), descartada

Energía complementaria a daño fijo, con $\boldsymbol\sigma^\pm$ la partición espectral del esfuerzo **nominal** (Mazars, Berthaud y Ramtani 1990; es la parte elástica del modelo de La Borderie 1991):

$$\chi(\boldsymbol\sigma, d^+, d^-) = \frac{\boldsymbol\sigma^+:\boldsymbol\sigma^+}{2E\,(1-d^+)} + \frac{\boldsymbol\sigma^-:\boldsymbol\sigma^-}{2E\,(1-d^-)} + \frac{\nu}{2E}\Bigl(\boldsymbol\sigma:\boldsymbol\sigma - (\operatorname{tr}\boldsymbol\sigma)^2\Bigr), \qquad \boldsymbol\varepsilon = \frac{\partial\chi}{\partial\boldsymbol\sigma}.$$

Es convexa con $\nu < 0.5$ y $d^\pm < 1$, así que la relación $\boldsymbol\varepsilon \mapsto \boldsymbol\sigma$ es única, $\psi = \tfrac12\boldsymbol\sigma:\boldsymbol\varepsilon$ y $\boldsymbol\sigma = \partial\psi/\partial\boldsymbol\varepsilon$ exactamente; disipación $Y^+\dot d^+ + Y^-\dot d^- \ge 0$ con $Y^\pm = \boldsymbol\sigma^\pm:\boldsymbol\sigma^\pm/(2E(1-d^\pm)^2)$; trabajo nulo en ciclos cerrados y tangente con daño fijo simétrica y definida positiva. El esfuerzo sale en ejes principales de $\boldsymbol\varepsilon$ probando los patrones de signo de los esfuerzos principales (`informes/prototipos_2026-11-11/modelo2/run_compl.py`). Tracción uniaxial exacta con $\varepsilon_\perp = -\nu(1-d^+)\,\varepsilon$; cierre completo. **Se descartó** porque con triángulos lineales no es objetiva (bloqueo de esfuerzos, tabla de la viga SEN); usarla exigiría un remedio cinemático (elementos mixtos de Cervera, Chiumenti y Codina 2010; discontinuidad embebida) o la transición a daño escalar de Jirásek y Zimmermann (1998).

### 2. Daño de tracción $d^+$

Igual que el modelo 1 con la medida de Rankine sobre el esfuerzo **efectivo** $\bar{\boldsymbol\sigma} = \mathbf C_e:\boldsymbol\varepsilon$ del estado 3D asociado ([IsotropicDamageCrackBand2D](IsotropicDamageCrackBand2D.md) §5) (Q4, decidida por el usuario el 2026-11-17):

$$\varepsilon^+_\text{eq} = \frac{\langle\bar\sigma_I\rangle_+}{E}, \qquad \kappa^+_{n+1} = \max(\kappa^+_n,\ \varepsilon^+_\text{eq}), \qquad
d^+(\kappa^+) = 1 - \frac{\kappa_0}{\kappa^+}\,e^{-\alpha(h)\,(\kappa^+-\kappa_0)}\ \ (\kappa^+ > \kappa_0),$$

con $\kappa_0 = f_t/E$, $\alpha(h) = 1/(G_f/(h f_t) - \kappa_0/2)$, $h < h_\text{max} = 2EG_f/f_t^2$ y el tope `damage_max` (modelo 1, §2–§3). Como la tracción uniaxial es exacta (§1), el material disipa $G_f/h$ por unidad de volumen y la prueba de objetividad de la barra vale sin cambios. En compresión $\varepsilon^+_\text{eq} = 0$: la compresión nunca hace crecer $d^+$.

### 3. Daño de compresión $d^-$ (Q2, decidida por el usuario el 2026-11-17)

Opcional: con $f_c$ y $G_c$ dados se activa; sin ellos $d^- \equiv 0$ (compresión elástica, como el prototipo de la viga).

**Medida equivalente** (Drucker-Prager sobre la parte negativa del esfuerzo efectivo, Faria, Oliver y Cervera 1998), normalizada para que en compresión uniaxial valga $|\bar\sigma|/E$:

$$\varepsilon^-_\text{eq} = \frac{3\,\bigl(K\,\bar\sigma^-_\text{oct} + \bar\tau^-_\text{oct}\bigr)}{(\sqrt2 - K)\,E}, \qquad K = \frac{\sqrt2\,(\beta - 1)}{2\beta - 1},$$

con $\bar\sigma^-_\text{oct} = \tfrac13\operatorname{tr}\bar{\boldsymbol\sigma}^-$, $\bar\tau^-_\text{oct} = \sqrt{\tfrac23 J_2(\bar{\boldsymbol\sigma}^-)}$ y $\beta = f_{bc}/f_c$ la razón entre las resistencias biaxial y uniaxial (1.16 por omisión, Kupfer). $J_2$ es siempre el del tensor 3D de $\bar{\boldsymbol\sigma}^-$, también en esfuerzo plano ($\bar\sigma_{zz} = 0$ cuenta en el desviador) y en deformación plana ($\bar\sigma^-_{zz}$ cuenta en la traza y en $J_2$). Compresión biaxial igual: el pico es $\beta f_c$. Compresión hidrostática: $\varepsilon^-_\text{eq}\,E/|\bar\sigma| = -3K/(\sqrt2-K) < 0$ (−0.414 con $\beta = 1.16$), así que nunca daña, rasgo de Drucker-Prager.

Faria et al. escriben la medida como $\bar\tau^- = \sqrt{\sqrt3\,(K\bar\sigma^-_\text{oct} + \bar\tau^-_\text{oct})}$ (su ec. 14b, con el mismo $K$, ec. 18); la forma lineal de arriba es una renormalización monótona con el mismo umbral y la misma superficie (leído en la fuente por el teórico). Su medida de tracción es la norma energética (ec. 14a); aquí se usa Rankine (Q4).

**Ley uniaxial** (parabólica de Feenstra 1993, la de DIANA, con la energía de aplastamiento $G_c$ regularizada por $h$), en la variable histórica $\kappa^- = \max \varepsilon^-_\text{eq}$:

$$\sigma_c(\kappa) = \begin{cases}
E\kappa & \kappa \le \kappa_{c/3} = \dfrac{f_c}{3E} \\[6pt]
\dfrac{f_c}{3}\Bigl(1 + 4\,\xi - 2\,\xi^2\Bigr),\ \ \xi = \dfrac{\kappa-\kappa_{c/3}}{\kappa_c-\kappa_{c/3}} & \kappa_{c/3} < \kappa \le \kappa_c = \dfrac{5f_c}{3E} \\[6pt]
f_c\Bigl(1 - \Bigl(\dfrac{\kappa-\kappa_c}{\kappa_u-\kappa_c}\Bigr)^2\Bigr) & \kappa_c < \kappa \le \kappa_u = \kappa_c + \dfrac{3\,G_c}{2\,h\,f_c} \\[6pt]
0 & \kappa > \kappa_u
\end{cases}
\qquad d^-(\kappa) = 1 - \frac{\sigma_c(\kappa)}{E\,\kappa},$$

con el mismo tope `damage_max`. La rama ascendente empieza con pendiente $E$ y llega a $f_c$ con pendiente nula; el área de la rama descendente es $\tfrac23 f_c(\kappa_u - \kappa_c) = G_c/h$, la regularización por banda de grieta de la compresión. La rama ascendente no se regulariza (es comportamiento del continuo antes del pico). $d^-$ crece con $\kappa^-$ en todo el intervalo, porque la secante $\sigma_c/\kappa$ decrece (calculado: $d(\sigma_c/E\kappa)/d\kappa = -3E\xi(2\xi+1)/(f_c(4\xi+1)^2) < 0$ en la rama ascendente).

**Elementos grandes (decidido por el usuario el 2026-11-17: sin límite).** DIANA limita la rama descendente: $\kappa_u = \max\bigl(\kappa_c + \tfrac{3G_c}{2hf_c},\ 2.5\,\kappa_c\bigr)$ (manual de teoría de DIANA 10.3, ecs. 42.111–42.115, en valor absoluto; leído por el teórico). El límite actúa cuando $h > 0.6\,E\,G_c/f_c^2$ (≈ 0.30 m con $E$ = 30 GPa, $f_c$ = 30 MPa, $G_c$ = 15 kN/m); por encima, el elemento disipa más que $G_c/h$. **No se adopta ese límite ni se rechaza $h$ grande**: con la fórmula de arriba la energía de aplastamiento es exactamente $G_c/h$ para cualquier $h$, a costa de una rama descendente casi vertical con $h$ grande (más pasos o retroceso de la estructura). Descartados: el límite de DIANA (los elementos grandes disiparían más que $G_c$ y la energía volvería a depender de la malla) y rechazar $h$ por encima del límite (obligaría a refinar).

### 4. Variables internas

| Variable | Significado |
|---|---|
| `kappa_t`, `kappa_c` | máxima deformación equivalente histórica de tracción y de compresión |
| `damage_t`, `damage_c` | daños $d^+$, $d^-$ con que se calcula el esfuerzo (`PRIMARY_STATE_VAR = damage_t`) |
| `h` | longitud característica del punto, sembrada por el elemento |
| `kappa_t_prev`, `kappa_c_prev` | $\kappa$ del paso confirmado anterior (extrapolación de IMPL-EX) |
| `damage_error` | indicador de IMPL-EX, $\max_\pm \lvert d^\pm(\kappa^\pm) - \tilde d^\pm\rvert$ |

### 5. Fuera de alcance (Q3, decidida por el usuario el 2026-11-17)

Sin deformaciones plásticas (descarga secante): la deformación residual del concreto aplastado y la histéresis en ciclos quedan para el modelo 4 (daño-plasticidad). Sin acoplamiento entre los dos daños (un concreto aplastado conserva su rigidez a tracción si $d^+ = 0$, como en Faria, Oliver y Cervera 1998; Abaqus CDP, por omisión, hace lo contrario: con $w_t = 0$ la tracción no se recupera tras aplastar). Sin recuperación parcial de la rigidez a cortante al cerrar la grieta (Cervera y Tesei 2017). Sin efectos de velocidad.

---

## Formulación numérica

### 6. Esfuerzo a partir de la deformación (forma explícita)

1. Esfuerzo efectivo del estado asociado (§1): $\bar{\boldsymbol\sigma} = \mathbf C_\text{ep}:\boldsymbol\varepsilon$ en esfuerzo plano; en deformación plana, $\bar{\boldsymbol\sigma}$ del plano con la elasticidad de deformación plana y $\bar\sigma_{zz} = \lambda(\varepsilon_{xx}+\varepsilon_{yy})$.
2. Esfuerzos y direcciones principales del plano en forma cerrada: $\bar\sigma_{1,2} = \tfrac12(\bar\sigma_{xx}+\bar\sigma_{yy}) \pm \sqrt{\tfrac14(\bar\sigma_{xx}-\bar\sigma_{yy})^2 + \bar\sigma_{xy}^2}$, $\theta = \tfrac12\operatorname{atan2}(2\bar\sigma_{xy},\ \bar\sigma_{xx}-\bar\sigma_{yy})$; en deformación plana, además $\bar\sigma_{zz}$ (eje $z$).
3. $\bar{\boldsymbol\sigma}^+ = \sum_i\langle\bar\sigma_i\rangle_+\,\mathbf n_i\otimes\mathbf n_i$.
4. $\boldsymbol\sigma = (1-d^-)\,\bar{\boldsymbol\sigma} - (d^+-d^-)\,\bar{\boldsymbol\sigma}^+$.

Sin iteración local. La medida de Rankine (§2) usa $\bar\sigma_I = \max(\bar\sigma_1, \bar\sigma_2)$: en deformación plana $\bar\sigma_{zz} = \nu(\bar\sigma_1+\bar\sigma_2)$ no lo supera cuando $\bar\sigma_I > 0$ (con $\nu \le 0.5$; comprobado, §11).

### 7. Tangente

**Con daño fijo** (IMPL-EX y descarga):

$$\mathbf C^{d} = (1-d^-)\,\mathbf C_e - (d^+-d^-)\,\mathbf Q^+:\mathbf C_e, \qquad
\mathbf Q^+ = \frac{\partial\bar{\boldsymbol\sigma}^+}{\partial\bar{\boldsymbol\sigma}} = \sum_i H(\bar\sigma_i)\,\mathbf M_{ii}\otimes\mathbf M_{ii} + \sum_{i<j}\tfrac12\,c_{ij}\,\mathbf S_{ij}\otimes\mathbf S_{ij},$$

con $\mathbf M_{ij} = \mathbf n_i\otimes\mathbf n_j$, $\mathbf S_{ij} = \mathbf M_{ij} + \mathbf M_{ji}$, $H$ la función escalón y $c_{ij} = (\langle\bar\sigma_i\rangle_+ - \langle\bar\sigma_j\rangle_+)/(\bar\sigma_i - \bar\sigma_j)$, con el límite $\tfrac12(H_i + H_j)$ si $\bar\sigma_i = \bar\sigma_j$, y $H(0) = 0$ de forma determinista (un esfuerzo principal efectivo nulo cuenta como no positivo; el esfuerzo no cambia, sólo evita que la tangente alterne entre iteraciones). En Voigt, las columnas de cortante de $\mathbf Q^+$ llevan el factor 2 (la variación de $\bar\sigma_{ij}$ y $\bar\sigma_{ji}$ a la vez). En deformación plana, para $\mathbf Q^+$ y $\mathbf C^d$ basta el bloque del plano: $\mathbf Q^+$ no acopla $\bar\sigma_{zz}$ con las componentes del plano. **No es simétrica** con $\nu \ne 0$, porque $\mathbf Q^+$ y $\mathbf C_e$ no conmutan; con $\nu = 0$ sí lo es. Igual a diferencias finitas a $1.2\times10^{-8}$ (§11).

**Implícita** (carga activa): se añaden los términos del daño,

$$\mathbf C^\text{alg} = \mathbf C^{d} - \bar{\boldsymbol\sigma}^+\otimes\frac{\partial d^+}{\partial\boldsymbol\varepsilon} - \bar{\boldsymbol\sigma}^-\otimes\frac{\partial d^-}{\partial\boldsymbol\varepsilon},
\qquad \frac{\partial d^+}{\partial\boldsymbol\varepsilon} = \frac{d^{+\prime}(\kappa^+)}{E}\,\mathbf C_e:(\mathbf n_I\otimes\mathbf n_I),
\qquad \frac{\partial d^-}{\partial\boldsymbol\varepsilon} = d^{-\prime}(\kappa^-)\,\frac{\partial\varepsilon^-_\text{eq}}{\partial\bar{\boldsymbol\sigma}^-}:(\mathbb I - \mathbf Q^+):\mathbf C_e,$$

la de $d^+$ como en el modelo 1 (§5 de su spec, Rankine) y la de $d^-$ con la medida de §3. No simétrica. **En deformación plana, $\partial d^-/\partial\boldsymbol\varepsilon$ no se reduce al bloque del plano**: $\varepsilon^-_\text{eq}$ depende de $\bar\sigma^-_{zz}$ (traza y $J_2$ del tensor 3D) y $\partial\bar\sigma_{zz}/\partial\boldsymbol\varepsilon = \lambda\,[1, 1, 0]$, así que se suma $\dfrac{\partial\varepsilon^-_\text{eq}}{\partial\bar\sigma^-_{zz}}\,\bigl(1 - H(\bar\sigma_{zz})\bigr)\,\lambda\,[1, 1, 0]$ (sin ese término, error del 8.6 % frente a diferencias finitas; con él, $2.8\times10^{-9}$; calculado por el teórico).

### 8. Integración implícita e IMPL-EX

Parámetro `integration` como en el modelo 1 (§7–§8 de su spec), aplicado a los dos daños: en IMPL-EX, $\tilde\kappa^\pm_{n+1} = \kappa^\pm_n + r\,(\kappa^\pm_n - \kappa^\pm_{n-1})$, $\tilde d^\pm_{n+1} = \max(\tilde d^\pm_n, d^\pm(\tilde\kappa^\pm_{n+1}))$, esfuerzo por §6 con los daños congelados y actualización implícita de $\kappa^\pm$ para el paso siguiente. A diferencia del daño isótropo, dentro del paso el material **no** es lineal (la partición cambia con la deformación), así que el Newton itera con $\mathbf C^{d}$: hasta tres resoluciones por paso en la viga SEN. Congelar la tangente después de la primera iteración no basta: en la viga a flexión con la banda restringida se detuvo en el paso 377 de 1000, y con Newton completo terminó (medido, §11).

### 9. Simetría

No simétrica en las dos integraciones (`IS_SYMMETRIC = False`), también con IMPL-EX: la tangente con daño fijo ya es asimétrica con $\nu \ne 0$ (§7).

### 10. Caveats numéricos

- **Sin potencial con $\nu \ne 0$** (§1): trabajo de los dos signos en ciclos cerrados con daño fijo; pequeño en carga monótona. Documentado y medido por `closed_cycle_work`.
- Contracción lateral $\varepsilon_\perp = -\nu\,\varepsilon$ con la grieta abierta, y su reverso: la compresión paralela a una grieta abierta se pierde con un umbral que depende de $h$ (§1). Registrado por `parallel_compression_open_crack`.
- Esfuerzos principales efectivos iguales: con $H_i = H_j$ el cociente $c_{ij}$ vale exactamente $H_i$ y se toma así, sin tolerancia (§Implementación); $H(0) = 0$. La función no es diferenciable en la arista de Rankine (doble autovalor máximo positivo con $d^+$ activo; error de 24 % frente a diferencias finitas) ni con un esfuerzo principal efectivo exactamente nulo en $\mathbf C^d$ (75 %); las barras a tracción uniaxial trabajan sobre ese quiebre ($\bar\sigma_\text{lateral} = 0$ ± redondeo), lo que no afecta al esfuerzo y puede costar iteraciones.
- Energía de tracción en deformación plana: $(1-\nu^2)\,G_f/h$ (§1), como el modelo 1.
- Con daños altos, el tope `damage_max` mantiene la rigidez positiva en tracción.
- Compresión con $h$ grande: $\kappa_u - \kappa_c$ pequeño (rama descendente casi vertical); retroceso a nivel de la estructura posible; sin límite, por decisión del usuario (§3).
- **Compresión con Quad8 y Quad9** (verificación del 2026-11-18, H1): la $h$ que siembran lleva el factor de localización $\alpha = 13/18$ de la tracción (decisión del usuario del 2026-11-12) y el material la usa también en $\kappa_u$. Con esfuerzo uniforme el elemento entero se ablanda y la banda entrega $G_c A/\alpha = 1.385\,G_c A$, objetivo entre mallas (medido; prueba `test_modelo2_barra_en_compresion_cuadraticos_con_factor_disipan_1_sobre_alfa`), igual que en tracción. Si el factor debe aplicarse también a la compresión está **pendiente de decisión del usuario**; hoy se aplica.
- **Energía total en compresión**: el trabajo de la banda desde su pico es objetivo ($G_c A$), pero el trabajo externo total incluye la energía de la banda antes del pico, proporcional a $h$: en la barra de `bar_in_compression`, 1.116, 1.081 y 1.064·$G_c A$ con $h$ = 40, 20 y 10 mm (calculado en forma cerrada y medido por el teórico, 2026-11-18). Es consecuencia de la definición de §3 (área de la rama descendente, como DIANA), no un defecto.
- **Deformación plana uniaxial** (lateral libre): $\bar\sigma_{zz} = \nu(\bar\sigma_{xx}+\bar\sigma_{yy})$ confina y el pico es $f_c/(E\,\varepsilon^-_\text{eq})$ del estado $(-1, 0, -\nu)$: 1.14·$f_c$ con $\nu = 0.2$ y $\beta = 1.16$ (calculado), rasgo de Drucker-Prager sobre $\bar{\boldsymbol\sigma}^-$.
- **Confinamiento lateral**: el cono de Drucker-Prager es abierto; con un confinamiento lateral $c$ veces el esfuerzo axial el pico crece (1.18·$f_c$ con $c$ = 0.1 y 1.84·$f_c$ con $c$ = 0.3, $\beta = 1.21$) y desde $c = (\sqrt2 - K)/(\sqrt2 + 2K)$ (0.707 con $\beta = 1.16$) la compresión no daña nunca (calculado por el teórico, 2026-11-18). Sin tapa en compresión: el aplastamiento con confinamiento alto queda fuera, como la plasticidad (modelo 4).
- IMPL-EX: error de primer orden en el inicio del daño (modelo 1); pasos finos cerca del pico; Newton completo dentro del paso (§8).
- Bloqueo de esfuerzos con triángulos: (a) lo atenúa (esfuerzo paralelo a la grieta en los puntos rotos de la viga SEN: −0.05 MPa frente a −8.4 MPa con (d)), porque la apertura entra en el esfuerzo efectivo y anula la compresión paralela (§1).
- Aviso de la literatura, sin comprobar: He et al. (2019, `docs/referencias/`) sostienen que la partición espectral puede no converger con cargas biaxiales o de inversión por el conflicto entre las superficies de tracción y de compresión, y usan el signo de la traza; vigilarlo en la validación.

### 11. Prototipos

- **Ley (a), punto material**: `informes/prototipos_2026-11-17/modelo2_ley_a/check_ley_a.py` (3D, esfuerzo plano y deformación plana). Tangente con daño fijo frente a diferencias finitas: $1.2\times10^{-8}$ en 2000 estados, asimetría hasta 30 %, nula con $\nu = 0$; tracción uniaxial $\sigma/(E\varepsilon) = 1 - d^+$ a $10^{-12}$ con $\varepsilon_\perp = -\nu\varepsilon$; cierre completo; ciclo radial: trabajo $= \tfrac12\boldsymbol\sigma:\boldsymbol\varepsilon$ exacto; ciclos cerrados con $\nu = 0$: $5.5\times10^{-8}$ (error de la cuadratura del trabajo); esfuerzo plano $=$ 3D con $\bar\sigma_{zz} = 0$ a $2.7\times10^{-14}$; deformación plana: $\bar\sigma_{zz} \le \bar\sigma_I$ cuando $\bar\sigma_I > 0$. Ciclos con $\nu = 0.2$: `ciclos.py` y `resultados_ciclos.txt`.
- **Ley (a), estructuras** (sólo $d^+$, esfuerzo plano, IMPL-EX): viga SEN, `informes/prototipos_2026-11-11/feap_tesis/diagnostico/run_diag.py` (casos `tc_L2`, `tc_L3`); viga a flexión en cuatro puntos con la banda restringida a una trayectoria prescrita, `informes/prototipos_2026-11-17/banda_orientacion/viga_trayectoria_m2.py` (los elementos rotos disipan 0.92–1.03 veces su energía nominal; con el daño isótropo, hasta 1.62).
- **Forma (d), descartada**: `informes/prototipos_2026-11-11/modelo2/run_compl.py`.
- **Verificación teórica independiente de la spec** (`teorico-solidum`, 2026-11-17): `informes/verificacion_teorica_2026-11-17_modelo2_spec.md`, guiones en `informes/prototipos_2026-11-17/modelo2_teorico/`. Ley en esfuerzo plano y deformación plana, $\mathbf C^d$, tangente implícita, medida de compresión con numerador 3, parábola de Feenstra y monotonía de $d^-$, $\bar\sigma_{zz} \le \bar\sigma_I$ y cifras de los ciclos deducidas por separado: coinciden. Las pruebas detectan las roturas a propósito (factor 2 de cortante, $(\mathbb I - \mathbf Q^+)$, término de $\bar\sigma_{zz}$, numerador $\sqrt3$).

---

## Contrato de implementación

```yaml
name: TensionCompressionDamage2D
kind: material
status: draft

interface:
  strain_dim: 3
  primary_state_var: damage_t
  is_symmetric: false      # siempre, también con IMPL-EX (ley (a), §9)
  state_schema: {kappa_t: [], kappa_c: [], damage_t: [], damage_c: [], h: [], kappa_t_prev: [], kappa_c_prev: [], damage_error: []}
  batch_kernel: true
  uses_pseudo_time: "true con integration='impl-ex' (ADR 0034)"

admission:                # Reglas §1 (Q5: el criterio no decide por la forma; el usuario clasifica: programa principal, 2026-11-18)
  textbooks:
    - de_souza_neto: "cap. 12, §12.6.1 (pp. 504–510, ecs. 12.83–12.117, Box 12.7) y §12.6.2 (pp. 510–511): cierre de grieta con la partición espectral tracción/compresión del tensor de esfuerzos (ecs. 12.93–12.97). La familia, con otra forma: la de Ladevèze y Lemaitre (forma (c) de §1), un solo daño D con la compresión dañada por hD (0 ≤ h ≤ 1), en flexibilidad y con la traza también partida (ec. 12.99). Texto leído el 2026-11-18"
  general_codes:
    - opensees: "ASDConcrete3D: la misma partición espectral del esfuerzo efectivo (ley (a)), daños de tracción y compresión con regularización por energía de fractura e IMPL-EX; sus criterios de daño son otros (Petracca et al. 2022)"
    - code_aster: "ENDO_ORTH_BETON: daño de tracción (tensor) y de compresión (escalar) con cierre de grieta, partición en deformación y potencial; ENDO_ISOT_BETON: un daño con partición de la deformación (MAZARS no recupera la rigidez al cerrar en 3D)"
    - abaqus: "Concrete Damaged Plasticity: un solo daño sobre el esfuerzo efectivo con recuperación de la rigidez al cerrar la grieta (pesos w_t, w_c), sin partición espectral"
    - diana: "Total Strain Crack Model: compresión parabólica de Feenstra con energía G_c y ancho de banda h"

parameters:
  - { name: E,          type: float, required: true,  desc: "Módulo de Young intacto [Pa]" }
  - { name: nu,         type: float, required: true,  desc: "Coeficiente de Poisson, en [0, 0.5)" }
  - { name: f_t,        type: float, required: true,  desc: "Resistencia a tracción [Pa]; κ_0 = f_t/E" }
  - { name: G_f,        type: float, required: true,  desc: "Energía de fractura en tracción [N/m]" }
  - { name: f_c,        type: float, required: false, default: null, desc: "Resistencia a compresión uniaxial [Pa], positiva; sin ella d⁻ ≡ 0" }
  - { name: G_c,        type: float, required: false, default: null, desc: "Energía de aplastamiento [N/m] (rama descendente de compresión); obligatoria con f_c" }
  - { name: beta,       type: float, required: false, default: 1.16, desc: "f_bc/f_c, razón de resistencias biaxial y uniaxial en compresión; β ≥ 1 (con β < 1 la compresión hidrostática dañaría)" }
  - { name: hypothesis, type: str,   required: false, default: "plane_stress", desc: "plane_stress | plane_strain" }
  - { name: integration, type: str,  required: false, default: "implicit", desc: "implicit | impl-ex (modelo 1, §8)" }
  - { name: damage_max, type: float, required: false, default: 0.999999, desc: "Tope de los dos daños, en (0, 1)" }
  - { name: density,    type: float, required: false, default: null, desc: "Densidad (ADR 0008)" }

state_variables:
  - { name: kappa_t,      type: float, shape: [], desc: "máxima ε⁺_eq histórica (Rankine)" }
  - { name: kappa_c,      type: float, shape: [], desc: "máxima ε⁻_eq histórica (Drucker-Prager sobre σ̄⁻)" }
  - { name: damage_t,     type: float, shape: [], desc: "d⁺ ∈ [0, damage_max]" }
  - { name: damage_c,     type: float, shape: [], desc: "d⁻ ∈ [0, damage_max]" }
  - { name: h,            type: float, shape: [], desc: "longitud característica del punto, sembrada por el elemento" }
  - { name: kappa_t_prev, type: float, shape: [], desc: "κ⁺ del paso confirmado anterior (IMPL-EX)" }
  - { name: kappa_c_prev, type: float, shape: [], desc: "κ⁻ del paso confirmado anterior (IMPL-EX)" }
  - { name: damage_error, type: float, shape: [], desc: "indicador de IMPL-EX; 0 en la implícita" }

conventions:
  voigt: "[ε_xx, ε_yy, γ_xy], γ_xy = 2·ε_xy; σ = [σ_xx, σ_yy, σ_xy]"
  units: "E, f_t, f_c en Pa; G_f, G_c en N/m; h en m (cualquier sistema coherente)"
  sign: "σ > 0 ⇔ tracción; f_c > 0"

validity:
  - "pequeñas deformaciones; E > 0, 0 ≤ ν < 0.5, f_t > 0, G_f > 0; con compresión, f_c > f_t y G_c > 0"
  - "0 < h < 2·E·G_f/f_t² en cada punto (como el modelo 1)"
  - "sólo en elementos que proporcionan h (sólidos 2D y 3D)"

out_of_scope:
  - "deformaciones plásticas, histéresis y deformación residual (modelo 4)"
  - "acoplamiento entre d⁺ y d⁻; recuperación parcial del cortante al cerrar la grieta"
  - "anisotropía del daño (grieta fija o rotante: modelo 3)"
  - "control adaptativo del paso de IMPL-EX"

acceptance:
  material_point:
    - name: uniaxial_tension_exact
      setup: "tracción uniaxial de esfuerzo en esfuerzo plano (y en 3D), d⁺ ∈ {0, 0.5, 0.9}"
      expect: "σ = (1 − d⁺)·E·ε; ε_⊥ = −ν·ε; las demás componentes del esfuerzo nulas"
      tol_rel: 1.0e-12
    - name: plane_strain_uniaxial
      setup: "deformación plana, tracción con σ_yy = 0, d⁺ ∈ {0, 0.5, 0.9} (valor independiente del 3D)"
      expect: "σ_xx = (1 − d⁺)·E·ε/(1 − ν²); σ_zz = ν·σ_xx"
      tol_rel: 1.0e-12
    - name: pure_shear_closed_form
      setup: "cortante puro γ en esfuerzo plano con daños fijos, p. ej. (d⁺, d⁻) = (0.8, 0.25), y con daño activo"
      expect: "σ = G·γ·[−(d⁺ − d⁻)/2, −(d⁺ − d⁻)/2, 1 − (d⁺ + d⁻)/2]; d⁺ empieza en τ = f_t (blinda γ = 2ε en el esfuerzo)"
      tol_rel: 1.0e-12
    - name: crack_closure
      setup: "tracción hasta d⁺ = 0.9 y después compresión uniaxial"
      expect: "rigidez en compresión E y Poisson ν (recuperación completa)"
      tol_rel: 1.0e-12
    - name: parallel_compression_open_crack
      setup: "esfuerzo plano, d⁺ = damage_max, d⁻ = 0, ε_t = −10⁻³ paralela a la grieta y apertura ε_n ∈ {0, 1, 2, 4, 5}·10⁻³"
      expect: "registra el caveat de §1: σ_tt = E/(1 − ν²)·(ε_t + ν·ε_n) mientras sea negativa y ≈ 0 (orden 1 − damage_max) desde ε_n = |ε_t|/ν"
      tol_rel: 1.0e-5
    - name: radial_unloading_secant
      setup: "deformaciones al azar con daños fijos; carga y descarga por la misma recta"
      expect: "σ(s·ε) = s·σ(ε) para s > 0; el trabajo devuelto en la descarga es ½·σ:ε"
      tol_rel: 1.0e-12
    - name: potential_when_nu_zero
      setup: "ν = 0, ciclos cerrados genéricos de deformación con daños fijos (d⁺, d⁻) ∈ {(0.9, 0), (0.9, 0.3)}; ≥ 16001 puntos de cuadratura por ciclo, o ψ de §1 evaluado en forma cerrada"
      expect: "∮σ:dε = 0 frente a la energía elástica máxima del ciclo (con ν = 0 la ley deriva de un potencial)"
      tol_rel: 1.0e-6
    - name: closed_cycle_work
      setup: "ν = 0.2, conjunto fijo (semilla) de ciclos cerrados genéricos con daños fijos, esfuerzo plano y 3D"
      expect: "registra el defecto documentado en §1: trabajo de los dos signos, mediana de |∮σ:dε|/ψ_max entre 1 y 10 %; los valores del conjunto, iguales a los de referencia"
      tol_rel: 1.0e-8
    - name: tension_energy
      setup: "tracción uniaxial monótona hasta d⁺ = damage_max en esfuerzo plano y 3D; h ∈ {5, 25, 100, 400} mm; en deformación plana con el lateral libre se espera (1 − ν²)·G_f/h (§1)"
      expect: "energía disipada = D(κ) del modelo 1 (§2); al final > 0.997·G_f/h"
      tol_rel: 1.0e-6
    - name: compression_curve
      setup: "compresión uniaxial monótona con f_c y G_c; h ∈ {10, 50, 100} mm"
      expect: "σ = −σ_c(κ) de §3 con |σ + σ_c(κ)| ≤ tol·f_c (escala f_c: σ → 0 cerca de κ_u); pico f_c en κ_c = 5f_c/(3E); área de la rama descendente G_c/h, descontado analíticamente el efecto del tope damage_max"
      tol_rel: 1.0e-8
    - name: biaxial_compression_peak
      setup: "compresión biaxial igual"
      expect: "pico β·f_c"
      tol_rel: 1.0e-8
    - name: hydrostatic_compression_no_damage
      setup: "compresión hidrostática creciente (3D) con f_c y G_c"
      expect: "ε⁻_eq < 0 y d⁻ = 0 en toda la historia (rasgo de Drucker-Prager, §3)"
    - name: damage_irreversible
      setup: "historias de carga, descarga, cierre y reapertura"
      expect: "d± nunca decrecen; la compresión no hace crecer d⁺"
    - name: tangent_fd
      setup: "daño fijo con esfuerzos principales efectivos distintos y dobles del mismo signo; implícita en carga activa con d⁺ o d⁻ (distintos; doble negativo con d⁻); esfuerzo plano y deformación plana, ésta con d⁻ activo y σ̄_zz de los dos signos; se excluyen la arista de Rankine y los esfuerzos principales efectivos nulos (§10)"
      expect: "|C_FD − C_tan| ≤ 1e-6·max|C_FD|; con ν = 0 y daño fijo, simétrica"
      tol_rel: 1.0e-6
    - name: restriction_of_3d
      setup: "historias con d⁺ y d⁻ activos y σ̄_zz de los dos signos, en esfuerzo plano (estado 3D con σ̄_zz = 0) y deformación plana (ε_zz = 0)"
      expect: "σ y variables internas iguales a TensionCompressionDamage3D con el estado 3D asociado"
      tol_rel: 1.0e-12
    - name: units_invariance
      setup: "el mismo caso en Pa–m y en MPa–mm, con todos los parámetros físicos ≠ 1 (E, ν, f_t, G_f, f_c, G_c, β, h)"
      expect: "esfuerzos, daños y energías iguales tras convertir unidades (Reglas §6)"
      tol_rel: 1.0e-12
    - name: batch_kernel_matches
      expect: "el kernel por lotes reproduce bit a bit σ, C_tan y la fila de estado de compute_state, implícita e IMPL-EX"
  impl_ex:
    - name: fixed_damage_tangent
      expect: "C_tan = C^d(d̃⁺, d̃⁻) de §7, igual a diferencias finitas con los daños extrapolados congelados"
      tol_rel: 1.0e-6
    - name: first_order_convergence
      setup: "como el modelo 1, en tracción (d⁺) y en compresión (d⁻), con pasos Δ, Δ/2, Δ/4, Δ/8"
      expect: "el error máximo de σ frente a la implícita se reduce a la mitad al reducir el paso a la mitad en el camino lineal (razón 2 ± 15 %) y a la cuarta parte en el cuadrático (4 ± 15 %); en compresión, a la cuarta parte en el camino lineal desde el estado virgen (d⁻ empieza con pendiente nula) y en el cuadrático, y a la mitad con descarga y recarga (medido el 2026-11-18, §Implementación)"
    - name: damage_monotone
      setup: "carga, descarga parcial y recarga con extrapolación que se pasa, en tracción y en compresión"
      expect: "d̃⁺ y d̃⁻ nunca decrecen"
    - name: bisection_idempotent
      setup: "mismo intento de paso dos veces tras una bisección"
      expect: "estado de prueba idéntico"
    - name: outside_step_frozen
      setup: "compute_state fuera de un paso (sin pseudo-tiempo)"
      expect: "κ± y d̃± sin cambio, como el modelo 1"
  mesh_objectivity:
    - name: bar_in_tension
      setup: "barra con una columna debilitada (la de la banda de grieta), los diez sólidos, 5–20 elementos; los 2D en esfuerzo plano; como el modelo 1: ν = 0, y además ν = 0.2 con Quad4; IMPL-EX con paso de 0.125 µm"
      expect: "D/(G_f·A) = 1 ± 1 % (Quad8 y Quad9: G_f·A/α = 18/13, factor de localización del modelo 1; con ν = 0.2: 1 − ν²) y carga posterior al pico (u = 0.05 mm) igual entre mallas ± 1 %"
      tol_rel: 1.0e-2
    - name: bar_in_compression
      setup: "la misma barra en compresión con f_c y G_c, de longitud fija L = 0.2 m con 5–20 elementos (h ≤ 40 mm: sin retroceso, du/dκ > 0 en toda la rama), control por desplazamiento"
      expect: "trabajo sobre los elementos de la banda desde su pico (∫σ dε en la banda) = G_c·A ± 2 % e igual entre mallas"
      tol_rel: 2.0e-2
  structural:
    - name: sen_beam
      setup: "viga SEN de Schlangen sobre la malla de FEAP del usuario y sus subdivisiones h/2, h/4, h/8 (Tri3), carga del ensayo, control indirecto del CMSD, IMPL-EX"
      expect: "banda diagonal hacia la placa de carga; D/(G_f·t·L) igual entre mallas ± 5 %; carga máxima entre mallas ± 3 %; residuo del balance de energía ≤ 1 % (caveat de la ley (a), §1; medido 0.6–0.9 % con sólo d⁺); sin comparación punto a punto con la curva experimental"
      tol_rel: 5.0e-2
    - name: van_vliet_tension
      setup: "el ensayo de Van Vliet del modelo 1"
      expect: "mismos criterios que el modelo 1 (energía frente a G_f·A ± 3 %)"
      tol_rel: 3.0e-2

references:
  - "Mazars J., Berthaud Y., Ramtani S. (1990). The unilateral behaviour of damaged concrete. Engineering Fracture Mechanics 35(4–5), 629–635. DOI 10.1016/0013-7944(90)90145-7"
  - "La Borderie C. (1991). Phénomènes unilatéraux dans un matériau endommageable : modélisation et application à l'analyse de structures en béton. Tesis doctoral, Université Paris 6."
  - "Ladevèze P., Lemaitre J. (1984). Damage effective stress in quasi unilateral conditions. 16th IUTAM Congress, Lyngby."
  - "Faria R., Oliver J., Cervera M. (1998). A strain-based plastic viscous-damage model for massive concrete structures. IJSS 35(14), 1533–1558. DOI 10.1016/S0020-7683(97)00119-4"
  - "Cervera M., Oliver J., Faria R. (1995). Seismic evaluation of concrete dams via continuum damage models. EESD 24(9), 1225–1245. DOI 10.1002/eqe.4290240905"
  - "Wu J.Y., Li J., Faria R. (2006). An energy release rate-based plastic-damage model for concrete. IJSS 43(3–4), 583–612. DOI 10.1016/j.ijsolstr.2005.05.038"
  - "Carol I., Willam K. (1996). Spurious energy dissipation/generation in stiffness recovery models for elastic degradation and damage. IJSS 33(20–22), 2939–2957. DOI 10.1016/0020-7683(95)00254-5"
  - "Ju J.W. (1989). On energy-based coupled elastoplastic damage theories: constitutive modeling and computational aspects. IJSS 25(7), 803–833. DOI 10.1016/0020-7683(89)90015-2"
  - "Cervera M., Tesei C. (2017). An energy-equivalent d+/d− damage model with enhanced microcrack closure-reopening capabilities for cohesive-frictional materials. Materials 10(4), 433. DOI 10.3390/ma10040433"
  - "Feenstra P.H. (1993). Computational aspects of biaxial stress in plain and reinforced concrete. Tesis doctoral, TU Delft."
  - "Code_Aster, doc. R7.01.04, loi ENDO_ISOT_BETON."
  - "OpenSees, ASDConcrete3D: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/ASDConcrete3D.html"
  - "Las del modelo 1 (Bažant y Oh 1983; Hillerborg et al. 1976; Rots 1988; Oliver, Huespe y Cante 2008)."
```

---

## Implementación

Plan: bloque 1, ley de §1 con $d^+$ (Rankine y banda de grieta), §6–§9 en 2D y 3D, implícita e IMPL-EX, pruebas del punto material y barra en tracción; bloque 2, daño de compresión de §3 con la barra en compresión; después, validación (viga SEN, Van Vliet) y `status: validated`.

**Bloque 1 (2026-11-18).** Archivos:

- `solidum/materials/_spectral_split.py`: $\bar{\boldsymbol\sigma}^+$ y $\mathbf Q^+$ (§7) en forma cerrada en 2D y con `eigh` en 3D; devuelve también $\bar\sigma_I$ y $\mathbf n_I\otimes\mathbf n_I$ para Rankine, sin una segunda descomposición. Esfuerzos principales iguales (§10): con $H_i = H_j$ se toma $c_{ij} = H_i$, que es el valor exacto del cociente, así que no hace falta tolerancia; con $H_i \ne H_j$ los dos esfuerzos tienen signo distinto y $\bar\sigma_i - \bar\sigma_j > 0$. Con todos los esfuerzos principales positivos (o no positivos), $\bar{\boldsymbol\sigma}^+ = \bar{\boldsymbol\sigma}$ y $\mathbf Q^+ = \mathbb I$ (o cero) exactamente.
- `solidum/materials/_tension_compression.py`: estado por punto y actualización de $d^+$ (implícita e IMPL-EX), común a 2D y 3D; reutiliza de `_crack_band.py` $\alpha(h)$, la ley exponencial con el tope, `implex_update`, la validación y la siembra de $h$ (`CrackBandMixin`).
- `solidum/materials/damage_tension_compression_2d.py` y `_3d.py`: núcleo por punto (§6–§8), adaptador por lotes (ADR 0014) y clase. En deformación plana el núcleo sólo necesita el plano (§6); $\sigma_{zz}$ para el post-proceso (`out_of_plane_stress`) se recupera de $\boldsymbol\sigma$ y de los daños, porque $\boldsymbol\sigma$ y $\bar{\boldsymbol\sigma}$ son coaxiales con $\sigma_i = (1-d^\pm)\,\bar\sigma_i$.
- Sólo evoluciona $d^+$: $d^-$ se arrastra del estado sin cambio (0 en un análisis) y la compresión es elástica. Los parámetros `f_c`, `G_c` y `beta` entran con el bloque 2.
- Caché de Numba: los núcleos de los dos materiales se compilan con `_spectral_split` y `_tension_compression` dentro; un cambio en esos módulos no invalida su caché (mismo caveat que el modelo 1 con `_equivalent_strain`).

Pruebas: `tests/test_tension_compression_damage.py` (aceptación `material_point` e `impl_ex` sin compresión; el trabajo en ciclos con $\nu = 0.2$ se compara ciclo a ciclo con el prototipo independiente, `informes/prototipos_2026-11-18/modelo2_bloque1/ref_ciclos.py`) y `tests/validation/test_crack_band_objectivity.py` (`bar_in_tension` con los diez sólidos, implícita e IMPL-EX). Mutaciones detectadas por las pruebas: `informes/prototipos_2026-11-18/modelo2_bloque1/mutaciones.py`.

**Bloque 2 (2026-11-18): daño de compresión (§3).**

- `_tension_compression.py`: medida de Drucker-Prager sobre $\bar{\boldsymbol\sigma}^-$ del tensor 3D y su gradiente $\mathbf g = \partial\varepsilon^-_\text{eq}/\partial\bar{\boldsymbol\sigma}^- = (K\,\mathbf I + \mathbf s/\bar\tau_\text{oct})/((\sqrt2-K)E)$ (`compression_measure`); parábola de Feenstra y $d^{-\prime}(\kappa) = (\sigma_c - \kappa\,\sigma_c')/(E\kappa^2)$ (`compression_damage`), con el umbral de admisibilidad del ADR 0006 a escala $f_c/3$ y el tope; actualización implícita e IMPL-EX con la misma estructura que $d^+$ (`compression_damage_update`); $\kappa_u$ con la $h$ del punto (`compression_ultimate`). Validación de `f_c`, `G_c` (van juntos, $f_c > f_t$, $G_c > 0$) y `beta`.
- Núcleos 2D y 3D: con `f_c` se calcula $\bar{\boldsymbol\sigma}^- = \bar{\boldsymbol\sigma} - \bar{\boldsymbol\sigma}^+$, en deformación plana con $\bar\sigma^-_{zz} = \min(\bar\sigma_{zz}, 0)$, y en carga activa se resta $\bar{\boldsymbol\sigma}^-\otimes\partial d^-/\partial\boldsymbol\varepsilon$, con $\partial d^-/\partial\boldsymbol\varepsilon = d^{-\prime}\,\mathbf g_v(\mathbb I - \mathbf Q^+)\mathbf C_e$ ($\mathbf g_v$ con el factor 2 de Voigt en el cortante) más, en deformación plana con $\bar\sigma_{zz} \le 0$, $d^{-\prime} g_{zz}\,\nu\,(\mathbf C_e[0,:] + \mathbf C_e[1,:])$ (que es $\lambda[1,1,0]$). Indicador de IMPL-EX: el mayor de los dos daños.
- **Estado inicial**: $\kappa^-_0 = f_c/(3E)$, fin del tramo elástico, como $\kappa^+_0 = \kappa_0$: así IMPL-EX nunca crea daño en un punto que la implícita deja elástico (con $\kappa^-_0 = 0$ la extrapolación del crecimiento elástico podría dar $\tilde d^- > 0$, irreversible, a un punto que se detiene justo bajo el umbral). La implícita no depende de este valor (diferencia medida: 0). Sin `f_c`, $\kappa^- = 0$ y $d^-$ se arrastra sin cambio (compresión elástica, como el bloque 1).
- **`beta ≥ 1`** (validación): con $\beta < 1$, $K < 0$ y la compresión hidrostática dañaría, contra la superficie de §3. El gradiente $\mathbf s/\bar\tau_\text{oct}$ es singular en el estado hidrostático con cualquier $\beta$; con $\beta \ge 1$ ahí nunca hay carga activa ($\varepsilon^-_\text{eq} \le 0 < \kappa^-$), así que no se usa. No hace falta cota superior ($\sqrt2 - K > 0$ para todo $\beta \ge 1$).
- Sin $h$ sembrada, el daño de compresión avisa igual que el de tracción.

Pruebas: `tests/test_tension_compression_damage.py`, sección «Bloque 2» (parámetros, `compression_curve` con la parábola escrita aparte en la prueba y el área de la rama descendente por cuadratura de Gauss con el tope descontado en forma cerrada, `biaxial_compression_peak` con β ∈ {1, 1.16, 1.3}, `hydrostatic_compression_no_damage`, la medida contra una evaluación independiente con numpy, `damage_irreversible` con los dos daños, `tangent_fd` con $d^-$ en carga en las dos ramas, con $d^+$ fijo o también en carga, con $\bar\sigma_{zz}$ de los dos signos y con esfuerzos principales negativos dobles, `restriction_of_3d`, `units_invariance`, `batch_kernel_matches` e IMPL-EX en compresión); `tests/validation/test_crack_band_objectivity.py`, `bar_in_compression` con Quad4, Tri3, Tri6 (5, 10, 20 columnas), Hex8 y Tet4 (5, 10). Mutaciones: `informes/prototipos_2026-11-18/modelo2_bloque2/mutaciones.py`.

Medido (2026-11-18): curva uniaxial igual a la parábola a $2\times10^{-15} f_c$; pico biaxial $\beta f_c$ exacto; compresión hidrostática sin daño; barra en compresión, trabajo de la banda desde su pico 0.998–1.000·$G_c A$ y pico $0.95 f_c A$ con Quad4, Tri3, Tri6, Hex8, Tet4, Hex20, Tet10 (y Hex27, fuera de la suite por tiempo), y $G_c A/\alpha = 1.385\,G_c A$ con Quad8 y Quad9 (§10). **IMPL-EX en compresión, camino lineal desde el estado virgen: segundo orden** (razón 4.0 frente a 2 en tracción): la extrapolación es exacta desde el segundo paso y el único error es el del primero, $(2f_c/3)\,\xi_1^2$, porque $d^-$ empieza con pendiente nula ($\sigma_c' = E$ en $\kappa_c/3$). No es una propiedad general de la compresión: **con descarga y recarga vuelve el primer orden** (razón 1.98–1.99, por el quiebre de $\kappa(t)$ al recargar dentro de un paso). La prueba de `first_order_convergence` en compresión exige 4 ± 15 % en el camino lineal desde el estado virgen y en el cuadrático, y 2 ± 15 % con descarga y recarga; otra prueba con pasos desiguales comprueba que la extrapolación con la razón $r$ es exacta desde el segundo paso y el valor del indicador.

Verificación independiente de `teorico-solidum` (`informes/verificacion_teorica_2026-11-18_modelo2_bloque2.md`): deducción y código coinciden al redondeo ($10^{-13}$ en 3600 estados por hipótesis); sin errores de formulación ni de implementación. Hallazgos de pruebas aplicados (H2–H7): tangente implícita en descarga con $d^- > 0$, IMPL-EX con razón variable y el indicador, aviso sin $h$ en 3D, descarga y recarga, una sola expresión de la tolerancia $f_c/3$ (bit a bit con los lotes), ramas contadas en la tangente y negativos dobles en las dos ramas y en deformación plana. Mutaciones: las del bloque 2 más las cinco que la verificación encontró sin detectar.

---

## Diálogo

- **2026-11-11** · Spec creada por la IA tras el «si» del usuario a «corregir $h$ de los cuadráticos y después la spec del modelo 2 antes del bloque 3». La hoja de ruta aprobada nombraba Faria, Oliver y Cervera (1998) / Mazars. Al preparar la spec la IA midió que la partición del esfuerzo efectivo de Faria et al. (a) no deriva de un potencial con $\nu \ne 0$, buscó una forma termodinámicamente consistente con tracción exacta y cierre completo, halló que la complementaria de Ladevèze y Lemaitre con la traza dañada (c) no es convexa y que la de Mazars, Berthaud y Ramtani (1990) (d) sí lo es, y probó (d) en la viga SEN con tres mallas: cierra el balance de energía, pero el bloqueo de esfuerzos la hace no objetiva (§1). Por eso la recomendación de Q1 es (a), con su caveat documentado.

- **2026-11-17** · **Dato nuevo para Q1** (prototipo de la banda de grieta con trayectoria prescrita y cuerda exacta, viga entallada a flexión en cuatro puntos, 10 mallas Quad4/Tri3; detalle en el Diálogo de [IsotropicDamageCrackBand2D](IsotropicDamageCrackBand2D.md) de la misma fecha). Con la partición (a) y $d^- = 0$, los elementos rotos disipan 0.92–1.03 veces la energía de la calibración uniaxial; con el daño isótropo, hasta 1.62 en las bandas de triángulos, por el daño que crece con estados mixtos de cortante y compresión. Error medio de la energía con la cuerda: 9.0 % (isótropo) → 6.1 % (a). Lo que queda es cinemático (elementos cruzados que no se rompen) e igual con los dos materiales.

- **2026-11-17** · **Q1: el usuario aprueba (a)** («sigamos con lo que tú recomiendas, sobre el modelo 2»). Al preparar la reescritura de la spec para (a), la IA comprobó la ley (`informes/prototipos_2026-11-17/modelo2_ley_a/`): tracción uniaxial exacta y lateral $\varepsilon_\perp = -\nu\varepsilon$ (no $-\nu(1-d^+)\varepsilon$ como en (d)); cierre completo; descarga secante exacta; tangente con daño fijo $(1-d^-)\mathbf C_e - (d^+-d^-)\,\mathbf Q^+\mathbf C_e$ igual a diferencias finitas a $10^{-8}$, asimétrica hasta 30 %; con $\nu = 0$ deriva de un potencial (trabajo nulo en ciclos); esfuerzo plano = restricción exacta del 3D. **Pero el trabajo en ciclos cerrados de deformación con daño fijo es mayor que el ±6 % citado en Q1**: con $\nu = 0.2$ y 400 ciclos al azar por caso, mediana 2–6 % de la energía elástica máxima del ciclo, percentil 90 de 6–19 % y máximo de 15–48 %, la mitad de las veces creando energía. En carga monótona sigue siendo pequeño (residuo de 0.6–0.9 % en la viga SEN). La cifra corregida no cambia la recomendación (el defecto aparece en ciclos con daño fijo, fuera del uso de este modelo sin deformaciones plásticas; (d) no es objetiva con triángulos): se informó al usuario y **(a) queda adoptada**; el defecto se documenta con estas cifras y una prueba de aceptación lo mide.

- **2026-11-17** · Q2, Q3, Q4 y Q6 decididas por el usuario con la recomendación (lista de abajo); Q5 la decide la IA con el criterio de admisión, con la evidencia que reúne `bibliotecario-solidum`. Al repasar §3 la IA halló y corrigió un error de la spec: la medida de compresión llevaba $\sqrt3$ en el numerador y debe llevar 3; con $\sqrt3$ la compresión uniaxial daba $\varepsilon^-_\text{eq} = |\bar\sigma|/(\sqrt3E)$ y el pico habría salido en $\sqrt3 f_c$ (calculado; el pico biaxial $\beta f_c$ sale bien con las dos).

- **2026-11-17** · **Revisión de la spec por `teorico-solidum`** (`informes/verificacion_teorica_2026-11-17_modelo2_spec.md`): consistente, sin errores de fórmula que cambien resultados. La IA aplicó sus correcciones de precisión: §7, en deformación plana $\partial d^-/\partial\boldsymbol\varepsilon$ lleva el término de $\bar\sigma_{zz}$ (H3); §3, $J_2$ del tensor 3D, compresión hidrostática sin daño y la forma original de Faria (H7, H11); §1, la deducción de Faria et al. desde $\psi$ sólo vale con $\nu = 0$ y OpenSees comparte la partición, no los criterios de daño (H6, H7); §1 y §10, energía $(1-\nu^2)G_f/h$ en deformación plana (H12) y compresión paralela a una grieta abierta con umbral que depende de $h$ (H1, documentada como caveat con la prueba `parallel_compression_open_crack`); §7 y §10, $H(0) = 0$ y los puntos no diferenciables (H5); `acceptance`: definiciones precisas de `bar_in_compression` (trabajo sobre la banda desde su pico, $L$ = 0.2 m), `tangent_fd`, `potential_when_nu_zero` y `compression_curve`, y pruebas nuevas (deformación plana con valor independiente, cortante puro, hidrostática, invariancia de unidades, bloque `impl_ex`) sin cambiar ninguna tolerancia aprobada (H4, H5, H8–H11). H2 (límite de DIANA para la rama descendente con elementos grandes): **el usuario decide no adoptarlo** (energía exacta con cualquier $h$; §3).

- **2026-11-17** · **Evidencia del criterio de admisión** (`bibliotecario-solidum`, `informes/busqueda_admision_modelo2_2026-11-17.md`): (ii) se cumple con la misma partición (OpenSees `ASDConcrete3D`, leído; sus criterios de daño son otros); (i) sólo por familia (de Souza Neto et al. 2008, §12.6, otra forma, por índice y código HYPLAS); la forma (a) no está en los ocho textos. El criterio no decide: Q5 pasa al usuario, con recomendación. La IA corrigió el bloque `admission` (Code_Aster: `ENDO_ORTH_BETON` en lugar de `MAZARS`, que no recupera la rigidez en 3D; Abaqus CDP: un solo daño sin partición espectral) y anotó dos datos: Abaqus CDP acopla por omisión la tracción al aplastamiento (§5) y He et al. (2019) advierten de problemas de convergencia de la partición espectral (§10). DOI de las referencias de la spec verificados en Crossref.

- **2026-11-18** · **Q5: el usuario clasifica el modelo en el programa principal** (opción recomendada). El criterio (ii) se cumple con la misma partición (OpenSees `ASDConcrete3D`) y el (i) por la familia (de Souza Neto et al. 2008, §12.6); la nota va al catálogo al validar. La lectura de pp. 504–512 del libro queda como refuerzo opcional, no como condición. Spec completa; la implementación empieza por el bloque 1 (§Implementación).

- **2026-11-18** · **de Souza Neto et al. (2008), §12.6, leído** (el PDF está en la biblioteca digital del usuario, `Mi_Biblioteca_Digital/Books/Constitutive models/Computational methods for plasticity. 9780470694527.pdf`; pp. 504–512 del libro = pp. 529–537 del PDF). Confirma lo que el bibliotecario había deducido del código HYPLAS: §12.6.1 presenta el modelo de Ladevèze (1983), Ladevèze y Lemaitre (1984) y Lemaitre (1996), con un solo daño $D$ y un factor de cierre $h \in [0, 1]$ ($h \approx 0.2$ típico; $h = 0$, cierre completo; $h = 1$, sin cierre). La deformación se expresa en función del esfuerzo nominal partido espectralmente en tracción y compresión, $\boldsymbol\sigma = \boldsymbol\sigma^+ + \boldsymbol\sigma^-$ (ecs. 12.93–12.97, la misma partición que §1 aplica al esfuerzo efectivo), con la traza también partida (ec. 12.99); el esfuerzo se obtiene por Newton en esfuerzos principales (Box 12.7). Es la forma (c) de la tabla de §1 con $d^+ = D$ y $d^- = hD$. §12.6.2 extiende el cierre a la evolución del daño dúctil de Lemaitre (tasa de liberación $Y$ partida, ec. 12.118), con plasticidad, fuera del alcance del concreto. No cita a Faria et al. ni a Mazars. Conclusión: el criterio (i) se cumple por la familia (daño unilateral con partición espectral del esfuerzo), no por la forma (a); la clasificación de Q5 no cambia.

- **2026-11-18** · **Bloque 1 implementado y verificado por `teorico-solidum`** (`informes/verificacion_teorica_2026-11-18_modelo2_bloque1.md`): deducción independiente y código coinciden al redondeo; sin errores de formulación ni de implementación; listo para commit. La IA aplicó sus hallazgos, todos de pruebas o de redacción: pruebas nuevas de la tangente implícita con $d^- \ne 0$ (H1), de $H(0) = 0$ (H2), de la tangente global ensamblada con partición de signo mixto, implícita e IMPL-EX (H3), de autovalores dobles girados (H4), tolerancia absoluta de `damage_error` (H5) y deformación plana en simetría y unidades (P4); `math.hypot` en la forma cerrada 2D (P2); la frase de §10 sobre la tolerancia de $c_{ij}$ (P1). **Aclaración de la aceptación `bar_in_tension` (H6)**: el texto aprobado decía «los diez sólidos… 1 ± 1 %»; las pruebas siguen, como el modelo 1, las decisiones del usuario del 2026-11-11 y 2026-11-12 (ν = 0; Quad8/Quad9 disipan $G_f A/\alpha$ = 18/13 en la barra de esfuerzo uniforme; con ν = 0.2, $1-\nu^2$), y con IMPL-EX el paso de 0.125 µm del modelo 1 (con 1 µm IMPL-EX disipa de más, igual en los dos modelos).

- **2026-11-18** · **Escalón «barra» de la validación (`validador-solidum`, `informes/reporte_modelo2_barra_2026-11-18.pdf`): validado con observaciones.** Objetivo entre mallas (Quad4, 5/10/20 columnas): $D/(G_f A)$ = 0.997–0.999 con ν = 0 y 0.958–0.961 con ν = 0.2 ($1-\nu^2$), carga en $u$ = 0.05 mm igual a 0.3 %; en tracción coincide con el modelo 1 (por construcción). **Cierre de grieta**: tras abrir la grieta hasta $d^+$ = 0.9995, la barra descarga por la secante y en compresión recupera exactamente $EA/L$ (el modelo 1 se queda en 0.0049·$EA/L$); el daño y $D$ no cambian al comprimir. Observaciones, iguales en el modelo 1: con IMPL-EX y pasos de 1 µm el pico se rebasa un 13–15 %, se daña toda la barra y la energía disipada crece al refinar (1.036, 1.062, 1.094·$G_f A$; hipótesis: el error depende de $\Delta u/h$); con 0.125 µm, 0.999–1.000. La implícita biseca una vez tras el pico. Seis huecos de Solidum como usuario (etapas con desplazamiento impuesto, historia con etapas, figuras), en `ESTADO.md`. Se reescribió el `setup`/`expect` para decirlo; no cambia ningún criterio ni tolerancia.

- **2026-11-18** · **Bloque 2 (daño de compresión) implementado** (§Implementación). Sin decisiones físicas nuevas: todo sale de §3 y §7. Dos precisiones de la IA, que se le dicen al usuario: (1) `beta ≥ 1` como validación (con $\beta < 1$ la compresión hidrostática dañaría, contra lo que dice §3); (2) IMPL-EX en compresión, camino lineal desde el estado virgen: segundo orden (medido: razón 4.0), porque $d^-$ empieza con pendiente nula; con descarga y recarga, primer orden como en tracción. Barra en compresión: trabajo de la banda 0.998–1.000·$G_c A$.

- **2026-11-18** · **Verificación independiente del bloque 2** (`teorico-solidum`, `informes/verificacion_teorica_2026-11-18_modelo2_bloque2.md`): consistente, listo para commit; deducción y código coinciden al redondeo. La IA aplicó los hallazgos de pruebas y redacción (H2–H7) y documentó en §10 los rasgos que calculó el teórico (energía total en compresión, deformación plana uniaxial, confinamiento). **H1 queda para el usuario**: Quad8 y Quad9 aplican en compresión el factor $\alpha = 13/18$ de la tracción y entregan $G_c A/\alpha$ en la barra de esfuerzo uniforme; documentado y protegido por una prueba, como en tracción.

**Preguntas al usuario** (una decisión física por pregunta; cada una con recomendación):

- **Q1. Ley de esfuerzo (partición tracción/compresión). Decidida: (a), Faria, Oliver y Cervera (1998)** (usuario, 2026-11-17; entradas del Diálogo de esa fecha). Razones: la única objetiva en la viga SEN con mallas finas, la ley estándar de la familia (OpenSees `ASDConcrete3D`) y la que suprime el exceso de energía de las bandas de triángulos; su defecto (no deriva de un potencial con $\nu \ne 0$) se documenta en §1 y lo mide `closed_cycle_work`. (d) es exacta en energía pero no es objetiva con triángulos lineales (bloqueo de esfuerzos); (b) descalibra la banda y cierra mal la grieta; (c) no es convexa.
- **Q2. Compresión. Decidida: aplastamiento opcional** (usuario, 2026-11-17): Drucker-Prager sobre la parte negativa del esfuerzo efectivo (Faria et al. 1998, con $\beta = f_{bc}/f_c$) y ley parabólica de Feenstra (DIANA) con $G_c/h$; sin $f_c$ y $G_c$ la compresión es elástica. Descartadas: compresión siempre elástica (no representa el aplastamiento) y la ley exponencial de Faria et al. (parámetros de ajuste, sin regularización por $G_c$: dependería de la malla).
- **Q3. Deformaciones plásticas. Decidida: no en este modelo** (usuario, 2026-11-17): descarga secante, sin deformación residual en descargas ni ciclos; la plasticidad del concreto aplastado va al modelo 4.
- **Q4. Criterio de tracción. Decidida: sólo Rankine sobre el esfuerzo efectivo, como el modelo 1** (usuario, 2026-11-17): la compresión no hace crecer $d^+$ (de ella se encarga $d^-$) y la calibración de la banda es exacta. Descartado ofrecer las otras medidas del modelo 1 (energía, Mazars, von Mises modificado): hacen crecer $d^+$ también con compresión y duplicarían el efecto de $d^-$.
- **Q5. Clasificación. Decidida: programa principal** (usuario, 2026-11-18, con la recomendación). En el catálogo, la nota «criterio (i) cumplido por la familia (de Souza Neto et al. 2008, §12.6), no por la forma (a)». Evidencia de `bibliotecario-solidum` (`informes/busqueda_admision_modelo2_2026-11-17.md`): el criterio (ii) se cumple con la misma partición del esfuerzo (OpenSees `ASDConcrete3D`, leído); el (i) sólo **por familia**: de Souza Neto et al. (2008), §12.6.1–12.6.2, trata el cierre de grieta en daño elástico, pero con otra forma (un solo daño y un factor de cierre $h$); la forma (a) no aparece en ninguno de los ocho textos (está en artículos: Faria et al. 1998, Cervera y Chiumenti 2006, Saloustros et al. 2018). **Recomendación: programa principal**, con la nota en el catálogo de que (i) se cumple por la familia y no por la forma; es el modelo estándar de concreto de la familia y el paso natural del modelo 1. La alternativa es un módulo de usuario. Afianzada el 2026-11-18 con la lectura de §12.6 de de Souza Neto (Diálogo de esa fecha); sigue sin leer §7.7 de de Borst et al. (2012).
- **Q6. Aceptación.** Casos y tolerancias del bloque `acceptance` (punto material exacto, descarga secante, potencial con $\nu = 0$, registro del trabajo en ciclos con $\nu = 0.2$, objetividad en barras de tracción y compresión, viga SEN entre mallas ±5 % en energía y ±3 % en carga máxima con residuo del balance ≤ 1 %, Van Vliet como el modelo 1; sin comparación punto a punto con curvas experimentales). **Decidida: aprobada así** (usuario, 2026-11-17).
