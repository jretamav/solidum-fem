# ManderConcrete1D — concreto uniaxial de Mander con descarga de daño-plasticidad

> Orden de trabajo. La especificación física y la formulación numérica las revisa el usuario con detalle; la IA implementa y rellena las secciones marcadas.

---

## Especificación física

### 0. Descripción general
Ley uniaxial de concreto (`STRAIN_DIM = 1`) en pequeñas deformaciones, independiente de la velocidad de deformación. Sirve a cualquier elemento que use un material 1D: barras de armadura y fibras de una `FiberSection` (ADR 0021) en los elementos de marco. Admite carga monótona, cíclica y dinámica, en análisis estáticos (incluidas las etapas de carga, ADR 0026) y en el tiempo, desde Python y desde el YAML.

**Clase de modelo**: fenomenológico (macroscópico), con estructura teórica de **daño-plasticidad 1D** y de nivel material. Los datos experimentales vienen de Mander, Priestley y Park (1988): la envolvente de compresión, la deformación plástica y la degradación de la tracción por compresión previa. La forma de la descarga y de la recarga la da el marco de daño-plasticidad, no las reglas curvas de Mander.

Decisiones del usuario (2026-11-09):

1. **Compresión**: envolvente de Mander con la ecuación de Popovics (§2); el confinamiento se calcula fuera del material (§6).
2. **Descarga y recarga en compresión**: lineales, con rigidez $(1-d)\,E$, hacia la deformación plástica de Mander (punto focal, ecs. 31–34); sin lazos en los ciclos parciales (§3).
3. **Tracción**: lineal hasta $f_t$ y rama descendente lineal hasta cero (§4). Con historia: descarga secante hacia el origen desplazado, recuperación de la rigidez de compresión al cerrarse la grieta, y degradación de la tracción por compresión previa según Mander (ecs. 42–45).
4. **Regularización opcional por energía de fractura** de las ramas descendentes, en tracción y en compresión (§5).
5. **Confinamiento**: ecuación de Mander para presiones laterales iguales y aproximación de Chang y Mander (1994) para presiones distintas (§6).

**Diferencias con Mander (1988)**, todas deliberadas:
- La descarga y la recarga son rectas (decisión 2), no las curvas de las ecs. 35–41 y 46–52.
- La tracción tiene rama descendente; Mander la anula al agrietarse (ec. 30b). Con $\varepsilon_{tu} = \varepsilon_t$ se recupera su caída brusca.
- La condición de la ec. 42 está corregida: el texto dice «si $\varepsilon_{pl} < \varepsilon_{cc}$ entonces $f_t = 0$»; debe ser $\varepsilon_{pl} \ge \varepsilon_{cc}$, como muestra la propia ecuación.
- Donde Mander escribe $\varepsilon_{cc}$ (ecs. 32, 33, 42) se usa la deformación de pico del material, $\varepsilon_c$: $\varepsilon_{cc}$ en el núcleo confinado, $\varepsilon_{c0}$ en el recubrimiento.

### 1. Convención de signos
$\sigma > 0$ y $\varepsilon > 0$ en tracción (Reglas §5). Los parámetros de compresión se dan como **magnitudes positivas**. En lo que sigue, $e = -\varepsilon \ge 0$ es la deformación de compresión y $s \ge 0$ el esfuerzo de compresión, con $\sigma = -s$.

### 2. Envolvente de compresión
Con $x = e/\varepsilon_c$, $E_{sec} = f_c/\varepsilon_c$ y $r = E/(E - E_{sec})$ (Mander 1988, ecs. 3–8):

$$s_{env}(e) = \begin{cases} \dfrac{f_c\, x\, r}{r - 1 + x^{r}} & 0 \le e \le \varepsilon_{cu} \\[2mm] s_{env}(\varepsilon_{cu})\,\dfrac{\varepsilon_{sp} - e}{\varepsilon_{sp} - \varepsilon_{cu}} & \varepsilon_{cu} < e \le \varepsilon_{sp} \\[2mm] 0 & e > \varepsilon_{sp} \end{cases}$$

- $f_c$ y $\varepsilon_c$ son los del material: $f'_{cc}$ y $\varepsilon_{cc}$ en el núcleo confinado, $f'_{c}$ y $\varepsilon_{c0}$ en el recubrimiento.
- **Núcleo confinado**: $\varepsilon_{cu}$ es la deformación última (rotura del primer estribo, §6) y, por omisión, $\varepsilon_{sp} = \varepsilon_{cu}$: el esfuerzo cae a cero al alcanzarla.
- **Recubrimiento** (Mander 1988, p. 1808): la curva de Popovics hasta $\varepsilon_{cu} = 2\varepsilon_{c0}$ y una recta hasta cero en la deformación de desprendimiento $\varepsilon_{sp}$.
- La pendiente inicial es $E$. Se exige $E > E_{sec}$ ($r > 1$). Mander toma $E = 5000\sqrt{f'_{c}}$ MPa (ec. 7); aquí $E$ es un parámetro, y la función auxiliar del §6 da ese valor por omisión.
- La secante $s_{env}(e)/e$ decrece con $e$, así que $s_{env}(e) \le E\,e$ en toda la envolvente. De esto dependen las propiedades del §3.

### 3. Deformación plástica, daño y descarga–recarga en compresión
Con $e_m$ la máxima compresión alcanzada y $s_m = s_{env}(e_m)$, la deformación plástica es la del punto focal de Mander (ecs. 31–34, p. 1814):

$$\varepsilon_a = a\sqrt{e_m\,\varepsilon_c},\qquad a = \max\!\left(\frac{\varepsilon_c}{\varepsilon_c + e_m},\ \frac{0.09\,e_m}{\varepsilon_c}\right),\qquad e_p = e_m - \frac{(e_m + \varepsilon_a)\,s_m}{s_m + E\,\varepsilon_a},$$

con $e_p = 0$ si $e_m = 0$. Mander la calibró con los ensayos de Karsan y Jirsa (1969) y de Sinha et al. (1964) en concreto simple, y con los de Mander et al. (1984) en concreto confinado.

**Forma de daño-plasticidad.** Entre $e_p$ y $e_m$ el material responde elásticamente con la rigidez dañada:

$$s = (1 - d)\,E\,(e - e_p),\qquad 1 - d = \frac{E_r}{E},\qquad E_r = \frac{s_m}{e_m - e_p} = \frac{s_m + E\,\varepsilon_a}{e_m + \varepsilon_a}.$$

La descarga y la recarga siguen la misma recta, y al pasar de $e_m$ el material vuelve a la envolvente. $d$ y $e_p$ sólo cambian sobre la envolvente ($e = e_m$, $\dot e > 0$).

**Propiedades por construcción** (consecuencia de $s_{env} \le E\,e$):
- $0 \le e_p \le e_m$;
- $E_r \le E$: la descarga nunca es más rígida que la carga inicial, es decir $0 \le d < 1$ mientras $s_m > 0$.

**Energía.** La energía libre de la parte de compresión es $\psi_c = \tfrac12 (1-d)\,E\,(e - e_p)^2 = \tfrac12\, s\,(e - e_p)$, y la disipación es $D = s\,\dot e_p + \tfrac12 E (e - e_p)^2\,\dot d$. Se comprobó con un prototipo (fuera del repositorio) que el trabajo realizado menos la energía almacenada no baja de cero en 1 200 trayectorias aleatorias de compresión, descarga, tracción, agrietamiento y aplastamiento. La violación mayor fue de $6\cdot10^{-8}\,f_c\varepsilon_c$, del orden del error de integración.

**Aplastamiento.** Si $s_m = 0$ (más allá de $\varepsilon_{sp}$), $e_p = e_m$ y el factor de tracción del §4 se anula: la fibra no resiste nada en adelante, ni en compresión ni en tracción.

### 4. Tracción y cierre de la grieta
La deformación de tracción se mide desde el origen desplazado: $t = \varepsilon + e_p$. Para $t < 0$ la fibra está en compresión (§3); para $t \ge 0$, en tracción.

**Envolvente sin daño** (decisión del 2026-11-09, como `Concrete02` de OpenSees), con $\varepsilon_t = f_t/E$:

$$\sigma^0_{t}(t) = \begin{cases} E\,t & 0 \le t \le \varepsilon_t \\[1mm] f_t\,\dfrac{\varepsilon_{tu} - t}{\varepsilon_{tu} - \varepsilon_t} & \varepsilon_t < t \le \varepsilon_{tu} \\[2mm] 0 & t > \varepsilon_{tu} \end{cases}$$

$f_t = 0$: sin tracción. $\varepsilon_{tu} = \varepsilon_t$: caída brusca (la de Mander, ec. 30b). Para representar el concreto entre grietas de un elemento reforzado, valores típicos de $\varepsilon_{tu}$ van de 5 a 20 veces $\varepsilon_t$. Con energía de fractura, $\varepsilon_{tu}$ sale del §5.

**Degradación por compresión previa** (Mander 1988, ecs. 42–45, p. 1815, con la errata corregida):

$$\phi = \max\!\left(0,\ 1 - \frac{e_p}{\varepsilon_c}\right),\qquad \sigma_{t,env}(t) = \phi\,\sigma^0_t(t).$$

Equivale a Mander: resistencia $\phi f_t$, rigidez $E_t = \phi E$ y la misma deformación de agrietamiento $\varepsilon_t = f_t/E$. Valores para concreto sin confinar de 30 MPa:

| $\eta = e_m/\varepsilon_c$ | 0.5 | 1.0 | 1.5 | 2.0 | 2.21 |
|---|---|---|---|---|---|
| $\phi$ | 0.96 | 0.78 | 0.50 | 0.16 | 0 |

**Historia.** Con $t_m$ la máxima deformación de tracción alcanzada:
- Mientras $t_m \le \varepsilon_t$ (sin agrietar), la respuesta es elástica: $\sigma = \phi E t$.
- Agrietada ($t_m > \varepsilon_t$): para $0 \le t \le t_m$, **descarga secante** hacia el origen desplazado, $\sigma = \sigma_{t,env}(t_m)\, t/t_m$; para $t > t_m$, la envolvente.
- **Cierre de la grieta**: para $t < 0$ rige la compresión del §3, con su rigidez $(1-d)E$, que no depende del daño de tracción. Es la recuperación total de la rigidez de compresión (como el concreto de Abaqus con $w_c = 1$).

**Por qué no crea energía.** La energía almacenada en tracción es $\psi_t = \tfrac12\,\sigma\,t$. $\phi$ y $e_p$ sólo cambian sobre la envolvente de compresión, donde $t < 0$ y $\psi_t = 0$. $t_m$ conserva su valor cuando el origen se desplaza.

### 5. Regularización por energía de fractura (opcional)
Una rama descendente concentra la deformación en un punto de integración, y la energía que disipa es el área bajo la curva por la longitud de ese punto: sin regularizar, el resultado depende de la malla. Para cada rama descendente se da la deformación (§2, §4) **o** la energía; con la energía, la rama se ajusta a la **longitud característica $h$** del punto de integración, que proporciona el elemento:
- barra de armadura: $h = L$;
- elemento de marco con `FiberSection`: la longitud tributaria del punto de Gauss, $h = w_i\,L$, con $w_i$ el peso en $[0, 1]$ (Coleman y Spacone 2001).

**Tracción** (Hillerborg et al. 1976; Bažant y Oh 1983). Con ablandamiento lineal, el área bajo $\sigma^0_t$ es $\tfrac12 f_t\,\varepsilon_{tu}$, y se iguala a $G_f/h$:

$$\varepsilon_{tu} = \frac{2\,G_f}{f_t\,h}.$$

Exige $\varepsilon_{tu} > \varepsilon_t$, es decir $h < 2\,G_f\,E/f_t^2$. Si no se cumple, la rama tendría retroceso (*snap-back*) en el material: Solidum lo rechaza al construir el modelo y pide elementos más cortos. Ejemplo: $f_t = 3$ MPa, $E = 27.4$ GPa, $G_f = 140$ N/m da $h < 0.85$ m. Con la degradación del §4, la energía disipada es $\phi\,G_f$.

**Compresión** (Coleman y Spacone 2001). La parte de la envolvente posterior al pico se estira en el eje de deformaciones:

$$s_{env}(e) = s^{ref}_{env}\!\left(\varepsilon_c + \frac{e - \varepsilon_c}{\lambda}\right)\quad (e > \varepsilon_c),\qquad \lambda = \frac{G_{fc}}{h\,A_{post}},\qquad A_{post} = \int_{\varepsilon_c}^{\varepsilon_{sp}} s^{ref}_{env}\,de,$$

donde $s^{ref}_{env}$ es la envolvente del §2 sin regularizar. El área posterior al pico vale entonces $G_{fc}/h$, y $\varepsilon_{cu}$ y $\varepsilon_{sp}$ se escalan como $\varepsilon_c + \lambda(\cdot - \varepsilon_c)$. La rama ascendente no cambia. La deformación plástica del §3 se evalúa con la envolvente regularizada; $s_{env} \le E\,e$ se conserva, y con ella las propiedades del §3 (el prototipo lo comprobó con $\lambda = 0.3$ y $\lambda = 3$).

**Valor de $G_f$.** Si no se tiene el dato del ensayo, una función auxiliar da el del Model Code 2010: $G_f = 73\, f_{cm}^{0.18}$ N/m, con $f_{cm}$ en MPa. $G_{fc}$ no tiene valor por omisión.

### 6. Confinamiento (función auxiliar, fuera del material)
Calcula $f'_{cc}$, $\varepsilon_{cc}$ y $\varepsilon_{cu}$ del núcleo a partir del concreto sin confinar ($f'_{c}$, $\varepsilon_{c0} = 0.002$ por omisión) y del refuerzo transversal. El material sólo recibe el resultado.

**Presión lateral efectiva** (Mander 1988, ecs. 9–28): $f'_l = k_e\,f_l$, con
- sección circular: $f_l = \tfrac12\,\rho_s\,f_{yh}$ (ec. 18), $k_e = (1 - s'/2d_s)^2/(1 - \rho_{cc})$ con aros (ec. 14) o $(1 - s'/2d_s)/(1 - \rho_{cc})$ con espiral (ec. 15);
- sección rectangular: $f_{lx} = \rho_x f_{yh}$ y $f_{ly} = \rho_y f_{yh}$ (ecs. 25–26), con $k_e$ de la ec. 22 (anchos libres $w'_i$ entre barras longitudinales, separación libre $s'$, dimensiones del núcleo $b_c$ y $d_c$).

$k_e$ puede darse directamente o calcularse con esas ecuaciones.

**Resistencia confinada**, con $K = f'_{cc}/f'_{c}$:
- **Presiones iguales** (sección circular): Mander, ec. 29,
$$K = -1.254 + 2.254\sqrt{1 + 7.94\,f'_l/f'_c} - 2\,f'_l/f'_c.$$
- **Sección rectangular** (presiones $f'_{l1} \le f'_{l2}$, iguales o no): Chang y Mander (1994), ecs. 3-70a–e, p. 3-32, con $q_i = f'_{li}/f'_c$, $\bar x = (q_1 + q_2)/2$ y $r = q_1/q_2$:
$$K = 1 + A\,\bar x\left(0.1 + \frac{0.9}{1 + B\,\bar x}\right),\quad A = 6.8886 - (0.6069 + 17.275\,r)\,e^{-4.989\,r},\quad B = \frac{4.5}{\frac{5}{A}\left(0.9849 - 0.6306\,e^{-3.8939\,r}\right) - 0.1} - 5.$$
  Aproxima el ábaco de Mander (Fig. 4, superficie de cinco parámetros de Willam y Warnke). Con presiones iguales difiere de la ec. 29 en menos de 0.003 en $K$; se usa también en ese caso para que la función no tenga saltos. Reproduce el ejemplo numérico de Mander (p. 1812: $f'_c = 30$ MPa, $f'_{ly} = 2.7$ MPa, $f'_{lx} = 5.1$ MPa): $K = 1.647$, frente a 1.65 leído en el ábaco. Validez: $q_2 \le 0.3$; por encima, aviso.

**Deformaciones**: $\varepsilon_{cc} = \varepsilon_{c0}\,[1 + 5(K - 1)]$ (ec. 5). Deformación última por balance de energía, en la aproximación de Paulay y Priestley (1992, §3.2.2): $\varepsilon_{cu} = 0.004 + 1.4\,\rho_s\,f_{yh}\,\varepsilon_{su}/f'_{cc}$, con $\rho_s = \rho_x + \rho_y$ en la sección rectangular y $\varepsilon_{su}$ la deformación del acero transversal a la resistencia máxima.

**Módulo**: $E = 5000\sqrt{f'_c}$ MPa (ec. 7) por omisión.

### 7. Variables internas
| Variable | Significado |
|---|---|
| `eps_cmax` | $e_m \ge 0$: máxima deformación de compresión alcanzada (`PRIMARY_STATE_VAR`). |
| `eps_tmax` | $t_m \ge 0$: máxima deformación de tracción alcanzada, medida desde $-e_p$. |

$e_p$, $d$ y $\phi$ se deducen de $e_m$; el aplastamiento, de $s_m = 0$. La longitud característica $h$ no es una variable interna: es un dato fijo del punto de integración.

---

## Formulación numérica

### 9. Actualización del estado
Explícita, sin ecuación local que resolver. Dado el estado comprometido $(e_m, t_m)$ y la deformación total $\varepsilon_{n+1}$:
1. $e = -\varepsilon_{n+1}$; si $e > e_m$, la fibra está en la envolvente de compresión y $e_m^{trial} = e$.
2. Con $e_m^{trial}$ se calculan $e_p$, $d$ y $\phi$, y $t = \varepsilon_{n+1} + e_p$.
3. $t < 0$: compresión (envolvente o recta de daño). $t \ge 0$: tracción (envolvente si $t > t_m$, con $t_m^{trial} = t$; elástica o secante si no).

### 11. Tangente consistente
La derivada de esa función por ramas respecto de $\varepsilon$, con el estado comprometido:
- envolvente de compresión: $ds_{env}/de$, que en la curva de Popovics es
$$\frac{ds_{env}}{de} = \frac{f_c\, r\,(r-1)\,(1 - x^{r})}{\varepsilon_c\,(r - 1 + x^{r})^2}$$
(negativa tras el pico); $1/\lambda$ veces el valor de referencia en la rama regularizada; $-s_{env}(\varepsilon_{cu})/(\varepsilon_{sp} - \varepsilon_{cu})$ en la rama lineal;
- recta de daño en compresión: $(1-d)\,E$;
- tracción: $\phi E$ elástica, $-\phi f_t/(\varepsilon_{tu} - \varepsilon_t)$ en la rama descendente, $\sigma_{t,env}(t_m)/t_m$ en la secante;
- $0$ en las zonas sin esfuerzo (aplastada, agrietada por completo).

Sobre la envolvente, $e_p$, $d$ y $\phi$ cambian, pero el esfuerzo sólo depende de $e$, así que la tangente es la de la envolvente. En los quiebres (picos, $\varepsilon_{cu}$, $\varepsilon_{sp}$, inicio de la descarga, cierre de la grieta) se toma el valor lateral de la rama en la que se está. Es simétrica (escalar).

### 12. Caveats numéricos
- **Ablandamiento sin regularizar**: tras los picos de tracción y de compresión, la deformación se localiza y la respuesta depende de la malla (§5). Con $G_f$ y $G_{fc}$ deja de depender.
- **Tangente negativa**: en las ramas descendentes. Para seguir la respuesta estructural después del pico hace falta control de desplazamiento o de longitud de arco.
- **Caídas a cero** en $\varepsilon_{sp}$ (si $\varepsilon_{sp} = \varepsilon_{cu}$) y en $\varepsilon_{tu} = \varepsilon_t$: son discontinuidades del esfuerzo, y el Newton bisecará al cruzarlas.
- **Quiebre al cerrarse la grieta**: la rigidez salta de la de tracción a $(1-d)E$; es continuo en esfuerzo.

---

## Contrato de implementación

```yaml
name: ManderConcrete1D
kind: material
status: validated        # draft → implemented → validated

interface:
  strain_dim: 1
  primary_state_var: eps_cmax
  is_symmetric: true
  state_schema: {eps_cmax: (), eps_tmax: (), h: ()}   # h: longitud característica del punto, fija (ADR 0028)

admission:                # criterio de admisión (Reglas §1); clasificado estándar por el usuario (2026-11-09)
  textbooks:
    - paulay_priestley: "cap. 3, §3.2.2 Confined Concrete, pp. 98–104: f'cc, εcc y εcu de Mander et al. (1988); la forma exacta de la curva, en Mander (1988) y Popovics (1973)"
  general_codes:
    - opensees: "uniaxialMaterial Concrete04 (envolvente de Popovics/Mander, descarga y recarga lineales) y Concrete02 (tracción con rama lineal)"
    - abaqus: "Concrete Damaged Plasticity: daño-plasticidad con recuperación de la rigidez de compresión al cerrarse la grieta (w_c = 1) y regularización de la tracción por energía de fractura"

parameters:
  - { name: fc,      type: float, required: true,  desc: "Resistencia de pico en compresión, magnitud > 0 (f'cc en el núcleo, f'c en el recubrimiento)" }
  - { name: eps_c,   type: float, required: true,  desc: "Deformación de pico, magnitud > 0 (εcc o εc0)" }
  - { name: E,       type: float, required: true,  desc: "Módulo inicial > fc/eps_c (Mander: 5000·√f'c MPa)" }
  - { name: eps_cu,  type: float, required: true,  desc: "Fin de la curva de Popovics, magnitud ≥ eps_c (εcu en el núcleo, 2·εc0 en el recubrimiento)" }
  - { name: eps_sp,  type: float, required: false, default: null, desc: "Deformación en que el esfuerzo llega a cero por la rama lineal (≥ eps_cu); por omisión eps_cu: caída a cero en eps_cu" }
  - { name: ft,      type: float, required: false, default: 0.0,  desc: "Resistencia a tracción ≥ 0; 0 = sin tracción" }
  - { name: eps_tu,  type: float, required: false, default: null, desc: "Deformación de tracción en que el esfuerzo llega a cero (≥ ft/E). Si ft > 0 se da eps_tu o G_f, no ambos" }
  - { name: G_f,     type: float, required: false, default: null, desc: "Energía de fractura en tracción (N/m); regulariza la rama descendente con la longitud característica del punto de integración" }
  - { name: G_fc,    type: float, required: false, default: null, desc: "Energía de fractura en compresión (N/m); regulariza la rama posterior al pico" }
  - { name: density, type: float, required: false, default: null, desc: "Densidad (kg/m³), ADR 0008" }
  - { name: confinement, type: dict, required: false, default: null, desc: "Argumentos de mander_confined_concrete (§6); calcula fc, eps_c, eps_cu y E del núcleo. Excluye fc, eps_c y eps_cu" }

state_variables:
  - { name: eps_cmax, type: float, shape: (), desc: "Máxima deformación de compresión alcanzada, magnitud ≥ 0" }
  - { name: eps_tmax, type: float, shape: (), desc: "Máxima deformación de tracción alcanzada, medida desde la deformación plástica de compresión" }
  - { name: h, type: float, shape: (), desc: "Longitud característica del punto de integración, fija; 0 si el elemento no la proporcionó (ADR 0028)" }

conventions:
  sign: "σ > 0 y ε > 0 en tracción; parámetros de compresión en magnitud positiva"
  units: "[fc] = [ft] = [E] = Pa; [G_f] = [G_fc] = N/m; deformaciones adimensionales"

validity:
  - "pequeñas deformaciones; estado uniaxial de esfuerzo; independiente de la velocidad de deformación"
  - "carga monótona, cíclica y dinámica; análisis estático (con etapas) y en el tiempo"
  - "confinamiento por aros, espirales y estribos rectangulares con presiones laterales iguales o distintas (f'l/f'c ≤ 0.3)"

out_of_scope:
  - "curvas de descarga y recarga y lazos de los ciclos parciales (Mander 1988, ecs. 35–41 y 46–52; Chang y Mander 1994): si un caso real los pide, entran como otro modelo completo de la familia uniaxial, sin tocar éste"
  - "degradación por ciclos repetidos a la misma amplitud"
  - "efecto de la velocidad de deformación (Mander 1988, último apartado)"
  - "confinamiento por encamisado de acero o de polímeros reforzados con fibras"
  - "concreto 2D/3D: futura familia de daño-plasticidad tipo Lee-Fenves, que compartirá los parámetros"
  - "dependencia de la temperatura"

acceptance:
  verification:
    - name: puntos_notables
      setup: "Núcleo (fc = 34.5 MPa, eps_c = 0.0045, eps_cu = 0.015, E = 25.0 GPa) y recubrimiento (fc = 25 MPa, eps_c = 0.002, eps_cu = 0.004, eps_sp = 0.006, ft = 2.5 MPa, eps_tu = 10·ε_t)"
      expect: "σ(−eps_c) = −fc; σ(−eps_cu) igual a la fórmula de Popovics; σ(−eps_sp) = 0; σ(ε_t) = ft y σ(eps_tu) = 0; pendiente inicial E"
      tol_rel: 1.0e-10
    - name: deformacion_plastica_mander
      setup: "Descarga desde la envolvente en η = e_m/eps_c ∈ {0.25, 0.5, 1, 2, 3}, en núcleo y recubrimiento"
      expect: "e_p igual a las ecs. 31–34 de Mander evaluadas de forma independiente; 0 ≤ e_p ≤ e_m; E_r ≤ E"
      tol_rel: 1.0e-12
    - name: descarga_recarga_lineal
      setup: "Compresión hasta e_m en la rama ascendente y en la descendente, descarga hasta e_p, recarga hasta pasar e_m"
      expect: "Descarga y recarga por la misma recta de pendiente (1 − d)·E; variables internas constantes mientras e < e_m; vuelta a la envolvente en e_m"
      tol_rel: 1.0e-12
    - name: traccion_con_historia
      setup: "Compresión previa hasta η ∈ {0.5, 1, 1.5, 2, 2.5} y después tracción; tracción hasta agrietar, descarga y compresión"
      expect: "Resistencia φ·ft y rigidez φ·E con φ = max(0, 1 − e_p/eps_c) (0 para e_p ≥ eps_c); descarga secante hacia −e_p; al cerrarse la grieta la rigidez de compresión es (1 − d)·E, igual que sin agrietar"
      tol_rel: 1.0e-12
    - name: regularizacion_traccion
      setup: "Tracción hasta la rotura completa con G_f dado y h ∈ {0.05, 0.1, 0.5} m"
      expect: "Energía disipada por unidad de volumen = G_f/h; error al construir si h ≥ 2·G_f·E/ft²"
      tol_rel: 1.0e-10
    - name: regularizacion_compresion
      setup: "Compresión hasta eps_sp regularizado con G_fc dado y h ∈ {0.05, 0.1, 0.5} m"
      expect: "Área de la envolvente posterior al pico = G_fc/h; rama ascendente sin cambio"
      tol_rel: 1.0e-10
    - name: confinamiento
      setup: "Función auxiliar: sección circular y rectangular, con k_e dado y calculado"
      expect: "Ec. 29 de Mander con presiones iguales; Chang y Mander con presiones distintas, evaluados de forma independiente; diferencia entre ambas < 0.003 en K con presiones iguales y q ≤ 0.3; ejemplo numérico de Mander (p. 1812): K = 1.65 ± 0.01; εcc de la ec. 5"
      tol_rel: 1.0e-12
    - name: tangente_consistente
      setup: "Barridos de deformación por cada rama (envolvente antes y después del pico, regularizada o no, rama lineal, recta de daño, tracción elástica, rama descendente, secante agrietada), salvo en los quiebres"
      expect: "La tangente coincide con diferencias finitas centradas (patrón de tests/test_material_tangents_fd.py)"
      tol_rel: 1.0e-6
    - name: disipacion_no_negativa
      setup: "Trayectorias aleatorias de deformación que cruzan compresión, descarga, tracción, agrietamiento, cierre y aplastamiento, con y sin regularización"
      expect: "Trabajo realizado − energía almacenada ≥ 0 en todo instante, integrando el trabajo de forma exacta por tramos"
  specific:
    - name: familia_uniaxial
      setup: "Prueba común aplicada automáticamente a todo material 1D registrado (sin darlo de alta): ciclos cerrados de deformación y tangente por diferencias finitas"
      expect: "Energía disipada ≥ 0 en ciclos cerrados; tangente consistente fuera de los quiebres"
  validation:
    - name: momento_curvatura
      setup: "Sección rectangular de concreto reforzado (núcleo confinado con la función auxiliar, recubrimiento y barras de Elastoplastic1D en una FiberSection), sin carga axial y con ella"
      expect: "Momento de primera fluencia del acero, momento máximo y curvatura en εcu iguales a los de una integración independiente de la sección (franjas finas, equilibrio por bisección del eje neutro, sin código de Solidum)"
      tol_rel: 2.0e-2
    - name: voladizo_concreto_reforzado
      setup: "Voladizo de concreto reforzado sin carga axial y con carga axial constante (etapa previa, ADR 0026) y empuje lateral, 16 elementos de Frame2DEuler con la sección anterior"
      expect: "Hasta el pico, la fuerza a desplazamiento dado coincide con la integración de la curvatura de la solución independiente"
      tol_rel: 2.0e-2
    - name: objetividad_de_malla
      setup: "(a) Barra de concreto simple a tracción con un elemento debilitado, 2, 4, 8 y 16 elementos, con G_f; (b) el voladizo anterior con N = −0.2·f'c·A_g después del pico, con 8 y 16 elementos, eps_tu en tracción y G_fc en compresión"
      expect: "(a) curva fuerza-desplazamiento analítica de la banda de fisuración y energía disipada A·G_f·(1 − σ/f_t) con cualquier malla (sin regularizar, proporcional al tamaño del elemento); (b) fuerza a desplazamiento dado después del pico con diferencia ≤ 5 % entre 8 y 16 elementos (decisión del usuario, 2026-11-09). Con 4 elementos, documentado y no exigido: la malla es más gruesa que la rótula plástica"
      tol_rel: 5.0e-2
    - name: respuesta_ciclica_y_dinamica
      setup: "El voladizo anterior con desplazamiento lateral cíclico de amplitud creciente (IndirectDisplacementSolver por etapas) y con un pulso de carga lateral seguido de vibración libre (NewmarkSolver)"
      expect: "Ciclos completos sin pérdida de convergencia, con fuerza residual tras los ciclos plásticos y picos simétricos; en el transitorio, la energía disipada por el material (trabajo de las fuerzas internas, F_ext − M·ü, menos la energía almacenada en las fibras) no negativa y no decreciente"
      tol_rel: 1.0e-2
    - name: yaml_equivalente
      setup: "El voladizo desde el YAML, con el confinamiento calculado por la función auxiliar desde el propio YAML"
      expect: "Resultados idénticos a los de Python"
      tol_rel: 1.0e-12

references:
  - "Mander J. B., Priestley M. J. N., Park R. (1988). Theoretical stress-strain model for confined concrete. J. Struct. Eng. 114(8), 1804-1826. doi:10.1061/(ASCE)0733-9445(1988)114:8(1804)."
  - "Popovics S. (1973). A numerical approach to the complete stress-strain curve of concrete. Cem. Concr. Res. 3(5), 583-599. doi:10.1016/0008-8846(73)90096-3."
  - "Chang G. A., Mander J. B. (1994). Seismic energy based fatigue damage analysis of bridge columns: Part I – Evaluation of seismic capacity. Technical Report NCEER-94-0006, National Center for Earthquake Engineering Research, Buffalo. Ecs. 3-70a–e, p. 3-32."
  - "Karsan I. D., Jirsa J. O. (1969). Behavior of concrete under compressive loadings. J. Struct. Div. ASCE 95(12), 2543-2564. doi:10.1061/JSDEAG.0002424. (Datos de calibración de la deformación plástica de Mander.)"
  - "Hillerborg A., Modéer M., Petersson P.-E. (1976). Analysis of crack formation and crack growth in concrete by means of fracture mechanics and finite elements. Cem. Concr. Res. 6(6), 773-781. doi:10.1016/0008-8846(76)90007-7."
  - "Bažant Z. P., Oh B. H. (1983). Crack band theory for fracture of concrete. Matériaux et Constructions 16(3), 155-177. doi:10.1007/BF02486267."
  - "Coleman J., Spacone E. (2001). Localization issues in force-based frame elements. J. Struct. Eng. 127(11), 1257-1265. doi:10.1061/(ASCE)0733-9445(2001)127:11(1257)."
  - "Lee J., Fenves G. L. (1998). Plastic-damage model for cyclic loading of concrete structures. J. Eng. Mech. 124(8), 892-900. doi:10.1061/(ASCE)0733-9399(1998)124:8(892). (Marco de daño-plasticidad del concreto de Abaqus.)"
  - "fib (2013). fib Model Code for Concrete Structures 2010. Ernst & Sohn. doi:10.1002/9783433604090. (G_f = 73·f_cm^0.18.)"
  - "Paulay T., Priestley M. J. N. (1992). Seismic Design of Reinforced Concrete and Masonry Buildings. Wiley, New York. ISBN 0-471-54915-0. Cap. 3, §3.2.2."
  - "ADR 0021 (secciones por fibras), ADR 0026 (etapas de carga)."
```

---

## Implementación

- Archivo: `solidum/materials/concrete_1d.py` (material) y `solidum/materials/concrete_confinement.py` (confinamiento, `k_e`, módulo de Mander, energía de fractura del Model Code).
- Clase: `ManderConcrete1D` (registrada en `MaterialRegistry`).
- Pieza general nueva: la longitud característica del punto de integración ([ADR 0028](../adr/0028-longitud-caracteristica-del-punto-de-integracion.md)): `Material.uses_characteristic_length`, `Material.initial_point_state(h)`, `Element._seed_characteristic_lengths`, `FiberSection.initial_states(characteristic_length=h)`; la siembran las barras (`h = L`) y los marcos con sección (`h = w_g·L`).
- Tests:
  - `tests/test_mander_concrete_1d.py` (verificación: 20 casos).
  - `tests/test_uniaxial_family.py` (prueba común de todo material 1D: trabajo desde el estado virgen ≥ 0, ciclos cerrados con respuesta repetida ≥ 0 —postulado de Ilyushin— y tangente por diferencias finitas fuera de los quiebres).
  - `tests/validation/test_mander_concrete_validation.py` (validación).
  - Muestra en `MUESTRAS` de `tests/test_material_contract_sweep.py` (con tracción: sin ella no hay régimen elástico para ε > 0).
- Notas de implementación:
  - Actualización explícita, sin ecuación local. `e_p`, `d` y `φ` se recalculan de `eps_cmax` en cada llamada (no se guardan).
  - `A_post` (área de la rama posterior al pico sin regularizar) se integra una vez al construir, con 64 puntos de Gauss-Legendre.
  - `compression_envelope(e, h)`, `plastic_strain(e_m, s_m)` y `tension_softening_end(h)` son públicos (los usan las pruebas y el post-proceso de energía).
- Notas de las pruebas:
  - La disipación no negativa se integra con la regla del trapecio en pasos finos (tolerancia `1e-6·fc·eps_c`), no de forma exacta por tramos.
  - El balance de energía del transitorio con `F_int = F_ext − M·ü` se cumple de forma idéntica con la regla del trapecio de Newmark (aceleración media), así que no prueba nada: se sustituyó por la disipación del material, trabajo interno menos energía almacenada en las fibras, en tres instantes.

---

## Diálogo

- **2026-11-09** · Spec redactada por la IA con las tres decisiones físicas del usuario (compresión de Mander con la forma de Popovics y confinamiento fuera del material; tracción lineal con rama descendente lineal; descarga de Karsan y Jirsa) y su clasificación como estándar. Pendiente de la revisión del usuario, en particular: (a) la cota $E_r \le E$ de la descarga; (b) la deformación de tracción medida desde $-e_p$ y la descarga en tracción por la secante; (c) el confinamiento uniforme en la función auxiliar; (d) la rama lineal de $e_p$ para $\eta \ge 2$.
- **2026-11-09** · Referencias verificadas por el bibliotecario (`verificacion_bib_2026-11-09.md`, fuera de Git). Paulay y Priestley (1992) presenta el modelo de Mander en su §3.2.2 (pp. 98–104): sirve como texto de la lista de admisión. Karsan y Jirsa: páginas 2543–2564. La expresión de $e_p$ se confirmó en el código de `Concrete04` de OpenSees, que la usa sólo para $\eta < 2$ y una rama lineal para $\eta \ge 2$ (punto (d)).
- **2026-11-09** · **Reformulación.** El usuario pidió que el material sea una capacidad general de un programa de elementos finitos (carga monótona, cíclica y dinámica; YAML), no a la medida de un caso, y sin parches: el borrador anterior combinaba Karsan y Jirsa, la rama lineal de Taucer et al. (1991) que usa OpenSees y una cota propia ($E_r \le E$), y creaba energía en $\eta = 2$. Decisiones del usuario: (1) forma de daño-plasticidad 1D con la deformación plástica de Mander (ecs. 31–34), que elimina los puntos (a) y (d); (2) descarga y recarga lineales, sin lazos (las reglas curvas de Mander sólo están definidas para ciclos completos y dieron problemas de convergencia en fibras; Chang y Mander tiene unas 15 variables y errores de energía de hasta 47 % por ciclo); (3) tracción con descarga secante, recuperación al cerrar la grieta y degradación de Mander (ecs. 42–45, con la errata corregida), que responde al punto (b); (4) regularización opcional por energía de fractura; (5) confinamiento con presiones distintas por Chang y Mander (1994), que responde al punto (c). El módulo $E_c = 5000\sqrt{f'_c}$ MPa sí está en Mander (1988), ec. 7, p. 1808. Puntos para la revisión del usuario: (i) la regularización de la compresión estirando la rama posterior al pico (§5): Coleman y Spacone la hicieron con la rama lineal de Kent y Park; estirarla conserva la forma de Mander; (ii) $h = w_i L$ en los elementos de marco de Solidum, que son en desplazamientos (Coleman y Spacone lo plantean para elementos en fuerzas): lo decide la prueba `objetividad_de_malla`, y si no se cumple se trae al usuario; (iii) $t_m$ conserva su valor cuando el origen de la tracción se desplaza por compresión posterior (§4).
- **2026-11-09** · **Implementación y hallazgos.** Las pruebas de verificación pasan. Hallazgos de la validación:
  - **Momento–curvatura**: coincide con la integración independiente en < 0.1 % (con y sin carga axial). Una primera discrepancia era del programa de prueba (Newton sin salvaguarda en ε₀ cuando ∂N/∂ε₀ < 0 por el aplastamiento del recubrimiento), no del material.
  - **Voladizo**: hasta F = 0.92·F_pico, error ≤ 0.8 % con 8 elementos. Cerca del pico, con carga axial, el error de δ a fuerza dada converge con la malla (23.6 %, 8 %, 2.5 % con 8, 16 y 32 elementos): es del elemento en desplazamientos en la meseta tras la fluencia, no del material. La prueba compara la fuerza a desplazamiento dado, bien condicionada, con 16 elementos (con 8 da 2.2 %).
  - **Barra a tracción con `G_f`**: objetividad exacta (2, 4, 8 y 16 elementos dan la curva analítica de la banda de fisuración a precisión de máquina; energía `A·G_f·(1 − σ/f_t)`).
  - **Rama de tracción en concreto reforzado**: con acero de poco endurecimiento (1 %), la pérdida de la rigidez entre grietas tras la fluencia hace bajar un poco el momento de la sección y la deformación se localiza: con 16 elementos el control de desplazamiento no pasa (retroceso). Corrige la afirmación de la decisión de tracción del 2026-11-09 («en secciones de concreto reforzado no localiza»): también localiza si el acero endurece poco.
  - **`G_f` en concreto reforzado**: regularizar con `G_f` la rama que representa la rigidez entre grietas hace depender de la malla la respuesta antes del pico (−5.4 % con 4 elementos frente a 16). Recomendación de uso, en el catálogo: en concreto reforzado, `eps_tu` en tracción y `G_fc` en compresión; `G_f` para el concreto simple.
  - **Objetividad del voladizo después del pico** (punto (ii)), diferencia máxima de la fuerza a desplazamiento dado frente a 16 elementos: sin regularizar, 16 % con 4 elementos y 7.8 % con 8 (y 16 elementos no terminan); con `eps_tu` en tracción y `G_fc` en compresión, 6.6 % con 4 y 3.1 % con 8. Con `h = L` en vez de `w_g·L`: 5.8 % y 3.1 %. La elección de `h` apenas influye; lo que queda con 4 elementos es error de discretización. **Pendiente de la decisión del usuario** el criterio de la prueba con 4 elementos.
  - **Transitorio**: el cierre de las grietas es un quiebre de la rigidez y el Newton puro rebota entre ramas (con dt = 5 ms falla en t = 0.63 s); con búsqueda lineal o con dt = 2.5 ms pasa. La prueba usa `line_search=True`; recomendado en el catálogo y en el manual.
- **2026-11-09** · **Decisión del usuario (punto (ii))**: la prueba de objetividad del voladizo exige ≤ 5 % entre 8 y 16 elementos (medido 1.9 % en la configuración de la prueba); el caso de 4 elementos queda documentado y no se exige (7.1 %; 16 % sin regularizar): es error de discretización del elemento en desplazamientos, más largo que la rótula plástica, no de la regularización. **El usuario pide anotar el elemento de marco en fuerzas para implementarlo** (ROADMAP, sección por fibras, etapa 3): es el remedio para mallas gruesas y el elemento para el que Coleman y Spacone (2001) plantearon esta regularización. Con todas las pruebas de `acceptance` verdes, `status: validated`.
