# ISOTROPIC DAMAGE CRACK BAND 2D — daño escalar regularizado por banda de grieta, con integración implícita o IMPL-EX

> Orden de trabajo. El usuario revisa **especificación física**, **formulación numérica** y **contrato**; la IA rellena **implementación** y responde en **diálogo**.
>
> **Variante de componente existente** (Reglas §4): subclase de `IsotropicDamage2D` ([spec](IsotropicDamage2D.md), validada), que no se modifica. Esta spec documenta **sólo lo que cambia**: (i) la velocidad de ablandamiento $\alpha$ deja de ser un dato y se calcula en cada punto de integración a partir de la energía de fractura $G_f$ y de la longitud característica $h$ del elemento (banda de grieta, Bažant y Oh 1983), para que la energía disipada por unidad de área de grieta sea $G_f$ con cualquier malla; (ii) una medida de deformación equivalente más, la de **Rankine**; (iii) una segunda forma de integrar el daño, **IMPL-EX** (Oliver, Huespe y Cante 2008); (iv) un tope de daño propio. Todo lo demás (esfuerzo nominal y efectivo, estado 3D asociado, variable histórica, Kuhn-Tucker, ley exponencial, medidas existentes, tangente consistente) es el del padre.
>
> Uso previsto: agrietamiento de estructuras de concreto con fisuración distribuida, en 2D y 3D (la variante 3D está en [IsotropicDamageCrackBand3D](IsotropicDamageCrackBand3D.md)).
>
> **Estado: aprobada por el usuario (2026-11-11, «si», con las cinco recomendaciones Q1–Q5). Bloque 1 implementado** (banda de grieta, medida de Rankine, tope propio y siembra de $h$ en los sólidos 2D y 3D). **IMPL-EX (§8) entra en el bloque 2**, con su ADR y la prueba estructural de Van Vliet; hasta entonces el material sólo se integra de forma implícita y no admite el parámetro `integration`.

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

**Decisión y resultado (2026-11-11).** Aprobada (a) con la regla de los triángulos. La prueba de objetividad (barra con una columna de cuadrados debilitada, partidos en dos triángulos) dio con $\sqrt{A_e}$ en `Tri3` $D/(G_f A) = 1.41$ (la banda la forman los dos triángulos de cada cuadrado: su ancho es $\sqrt{2A_e}$), así que `Tri3` usa **(b) $h = \sqrt{2A_e}$** y da 1.00, igual que `Quad4`. En 3D lo mismo con `Tet4` (seis tetraedros por cubo): $\sqrt[3]{V_e}$ daba 1.81 y $\sqrt[3]{6V_e}$ da 1.00. Los demás elementos usan el convenio de Abaqus. Factor $h = c\,\Omega_e^{1/d}$ declarado por cada elemento (`CHARACTERISTIC_LENGTH_FACTOR`): `Quad4` 1, `Tri3` $\sqrt2$, `Tri6`/`Quad8`/`Quad9` ½; `Hex8` 1, `Tet4` $6^{1/3}$, `Hex20`/`Hex27`/`Tet10` ½.

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
2. **Razón de incrementos de pseudo-tiempo para IMPL-EX** (§8): $r = \Delta\tau_{n+1}/\Delta\tau_n$, que hoy no llega al material. El solver la da al inicio de cada intento de paso (junto al gancho `prepare_all_steps` del ADR 0020, P5), con $\Delta\tau_n$ el incremento del último paso convergido. Un solver que no declare su pseudo-tiempo no admite materiales IMPL-EX: error al empezar el análisis.

---

## Formulación numérica

### 7. Integración implícita (`integration = 'implicit'`, por omisión)

La del padre (§6–§9 de su spec) con $\alpha = \alpha(h)$ del punto, el tope `damage_max` de la variante y la medida elegida (incluida Rankine, §5). Tangente consistente; se pasa a la secante en descarga, sin daño y al llegar al tope.

### 8. Integración IMPL-EX (`integration = 'impl-ex'`) — bloque 2, pendiente

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
- **Retraso de un paso**: el daño entra un paso después de que $\varepsilon_\text{eq}$ supere $\kappa_0$ (al inicio $\kappa_n = \kappa_{n-1} = \kappa_0$). El error de la extrapolación es de primer orden en $\Delta\tau$; cerca del pico hacen falta pasos pequeños.
- **Bisección**: repetir un paso con otro $\Delta\tau_{n+1}$ sólo cambia $r$; el estado convergido no cambia, así que el resultado es idempotente para el mismo intento.
- **Indicador de error** por punto: $e_{n+1} = \lvert d(\kappa_{n+1}) - \tilde d_{n+1}\rvert$, el máximo del modelo se informa por paso. El control adaptativo del paso con este indicador (Oliver et al. 2008, §5) queda fuera de alcance.

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
status: implemented        # bloque 1; validated al cerrar el bloque 2 (IMPL-EX y Van Vliet)

interface:
  strain_dim: 3
  primary_state_var: damage
  is_symmetric: false      # cota de clase; la instancia es true con equivalent_strain='energy' (y, en el bloque 2, con integration='impl-ex')
  batch_kernel: true       # h viaja en la fila de estado [κ, d, h]

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
  # integration: "implicit | impl-ex" entra en el bloque 2 (§8); hoy, sólo implícita.
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
  # kappa_prev (κ del paso convergido anterior) entra con IMPL-EX en el bloque 2.

conventions:
  sign: "σ > 0 ⇔ tracción"
  voigt: "γ_xy = 2·ε_xy; ε_eq sobre el estado 3D asociado (padre §3)"
  units: "E, f_t en Pa; G_f en N/m; h en m (cualquier sistema coherente)"
  regularization: "α(h) = 1/(G_f/(h·f_t) − κ_0/2), h < h_max = 2·E·G_f/f_t²"

validity:
  - "E > 0, f_t > 0, G_f > 0, ν ∈ (−1, 0.5)"
  - "0 < h < 2·E·G_f/f_t² en cada punto (se comprueba al recibir h)"
  - "damage_max ∈ (0, 1)"
  - "bloque 2: integration = 'impl-ex' sólo con un solver que declare su pseudo-tiempo (§6.2)"
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
  - "IMPL-EX: error de primer orden en el paso; el daño entra un paso tarde."
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
      setup: "camino monótono de tracción uniaxial con pasos Δ, Δ/2, Δ/4, Δ/8"
      expect: "el error máximo de σ frente a la implícita se reduce a la mitad al reducir el paso a la mitad (razón entre 1.7 y 2.3)"
    - name: no_negative_dissipation
      setup: "carga, descarga parcial y recarga con extrapolación que se pasa"
      expect: "d̃ nunca decrece; disipación por paso ≥ 0"
    - name: bisection_idempotent
      setup: "mismo intento de paso dos veces tras una bisección"
      expect: "estado de prueba idéntico"
  elements:
    - name: characteristic_length_seeded
      setup: "Tri3, Quad4, Tri6, Quad8, Quad9 con una malla regular y otra distorsionada"
      expect: "h de cada punto = √A_e (Quad4), √(2A_e) (Tri3, decisión de §4) o ½·√A_e (cuadráticos)"
      tol_rel: 1.0e-12
    - name: missing_seed_error
      expect: "un elemento que no siembre h con este material produce un error claro al construir el modelo"
  mesh_objectivity:
    - name: bar_in_tension
      setup: "barra o placa en tracción con un debilitamiento en una fila; 1, 2, 4 y 8 elementos por unidad de longitud en la dirección de la carga; Quad4 y Tri3 estructurados; implícita e IMPL-EX"
      expect: "energía disipada al final = G_f·A_grieta y la curva fuerza–alargamiento posterior al pico igual en todas las mallas, ±1 % (Quad4); en Tri3 se informa D/(G_f·A_grieta) para decidir Q2"
      tol_rel: 1.0e-2
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
  - "Oliver J., Huespe A.E., Cante J.C. (2008). An implicit/explicit integration scheme to increase computability of non-linear material and contact/friction problems. CMAME 197, 1865–1889."
  - "Code_Aster, doc. R5.03.81, méthode IMPLEX: https://ericca.uqtr.ca/fr14.4/man_r/r5/r5.03.81.pdf"
  - "OpenSees, ASDConcrete3D: https://opensees.github.io/OpenSeesDocumentation/user/manual/material/ndMaterials/ASDConcrete3D.html"
  - "Las del padre (Simó y Ju 1987; Mazars 1986; de Vree et al. 1995; de Borst et al. 2012; de Souza Neto et al. 2008)."
```

---

## Implementación

- Archivo: `solidum/materials/damage_crack_band_2d.py` (lo común con la 3D en `solidum/materials/_crack_band.py`; Rankine en `solidum/materials/_equivalent_strain.py`, código `RANKINE`, que los padres no admiten).
- Clase: `IsotropicDamageCrackBand2D(CrackBandMixin, IsotropicDamage2D)`; reutiliza el núcleo compilado del padre (`_compute_damage_2d`) con $\alpha(h)$ del punto y el tope `damage_max`.
- Siembra de $h$: pieza genérica del elemento base (`Element.CHARACTERISTIC_LENGTH_FACTOR`, `Element._seed_geometric_characteristic_length`): $h = c\,\Omega_e^{1/d}$ con $\Omega_e = \sum_g w_g \det J_g$ de la propia cinemática compilada del elemento; un elemento sin `CHARACTERISTIC_LENGTH_FACTOR` ni `SEEDS_CHARACTERISTIC_LENGTH` rechaza el material al construirse.
- Tests:
  - `tests/test_damage_crack_band.py` (punto material 2D y 3D: α(h), h_max, pico, energía $D(\kappa)$, tope, Rankine, tangente por diferencias finitas, restricción exacta del 3D, kernel por lotes bit a bit; siembra de $h$ en los diez sólidos, mallas regulares y distorsionadas; elemento sin $h$).
  - `tests/validation/test_crack_band_objectivity.py` (barra en tracción con Quad4, Tri3, Hex8 y Tet4).
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
