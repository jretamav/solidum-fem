# CONTACT PENALTY — contacto unilateral por penalización con fricción de Coulomb opcional

> Orden de trabajo pre-redactada por la IA el 2026-11-26 dentro del bloque de **contacto acotado** aprobado por el usuario el mismo día: superficies con pares conocidos y desplazamientos pequeños, ley de interfaz unilateral sin fricción con Coulomb opcional, sobre la familia de materiales de interfaz (ADR 0033). La consume el elemento de interfaz de espesor nulo (spec `Interface2D` / `Interface3D`, bloque hermano). El usuario revisa la física y aprueba antes de que la IA escriba código; las decisiones abiertas están en el **Diálogo**.

---

## Especificación física

### 0. Descripción general

Ley tracción–holgura para el **contacto entre dos superficies** con pares de puntos conocidos y desplazamientos pequeños. Material de interfaz de superficie (ADR 0033, `GEOMETRY = "surface"`), en 2D `(n, s)` y en 3D `(n, s1, s2)` (`JUMP_DIM = None`).

- **Normal**: condición unilateral de Signorini (no interpenetración, presión de contacto sólo de compresión, sin adhesión), impuesta por **penalización** con rigidez `K_n`.
- **Tangencial**: sin fricción (`μ = 0`, la superficie desliza libre) o **fricción de Coulomb** con coeficiente `μ`, regularizada con una adherencia elástica de rigidez `K_t` (formulación elastoplástica de la fricción: la adherencia es la «parte elástica» del deslizamiento y el deslizamiento por fricción la «parte plástica»).

No hay daño, ni cohesión, ni dilatancia: una superficie que se separa no transmite nada.

### 1. Holgura y desplazamiento relativo

En el sistema local de la interfaz, el elemento entrega a la ley

$$\mathbf g = \begin{bmatrix} g_n \\ \mathbf g_s\end{bmatrix}, \qquad g_n = g_0 + \Delta u_n, \qquad \mathbf g_s = \Delta\mathbf u_s$$

- `g_n`: **holgura normal actual**; `g_0` es la holgura inicial del par, que mide el elemento sobre la geometría de referencia (0 en una interfaz de espesor nulo con nodos coincidentes; definición en la spec del elemento). `g_n > 0` separación; `g_n < 0` interpenetración. La ley declara `ACCEPTS_INITIAL_GAP = True` (ADR 0043 §4): las leyes que no lo declaran (cohesivas, adherencia) sólo se admiten con `g_0 = 0`.
- `g_s`: desplazamiento relativo tangencial (1 componente en 2D, 2 en 3D), medido desde la configuración de referencia.

El deslizamiento se descompone en una parte de adherencia (elástica, recuperable) y una de deslizamiento por fricción (irreversible):

$$\mathbf g_s = \mathbf g_s^{e} + \mathbf g_s^{p}$$

### 2. Ley normal: contacto unilateral por penalización

$$t_n = K_n\,\langle g_n\rangle_-, \qquad \langle x\rangle_- = \min(x,\,0), \qquad p_N = -t_n = K_n\,\langle -g_n\rangle_+ \ \ge 0$$

con `p_N` la **presión de contacto**. En el límite `K_n → ∞` se recuperan las condiciones de Signorini (Kuhn–Tucker del contacto):

$$g_n \ge 0, \qquad p_N \ge 0, \qquad p_N\, g_n = 0$$

Con `K_n` finito la condición `g_n ≥ 0` se viola en `g_n = −p_N/K_n`: la interpenetración es el error de la penalización, inversamente proporcional a `K_n`.

### 3. Ley tangencial: fricción de Coulomb regularizada

$$\mathbf t_s = K_t\,(\mathbf g_s - \mathbf g_s^{p})$$

$$\Phi(\mathbf t_s, p_N) = \|\mathbf t_s\| - \mu\,p_N \le 0$$

$$\dot{\mathbf g}_s^{p} = \dot\lambda\,\frac{\mathbf t_s}{\|\mathbf t_s\|}, \qquad \dot\lambda \ge 0,\quad \Phi \le 0,\quad \dot\lambda\,\Phi = 0$$

- **Adherencia** (`Φ < 0`): la superficie se deforma elásticamente en el plano tangente con rigidez `K_t`; `K_t → ∞` es la adherencia perfecta de Coulomb.
- **Deslizamiento** (`Φ = 0`): `‖t_s‖ = μ·p_N`, en la dirección del deslizamiento.
- La regla de flujo no tiene componente normal: **sin dilatancia**. Como el límite `μ·p_N` depende de la holgura normal, la ley no deriva de un potencial: es **no asociada** y su tangente no es simétrica (§10).
- En 3D la ley es isótropa en el plano tangente: sólo depende de `‖t_s‖`.
- Con `μ = 0`: `t_s = 0` siempre, `K_t` no interviene.

### 4. Separación y recontacto

Con la superficie abierta (`g_n > 0`), `p_N = 0` y el límite de fricción es nulo: `t_s = 0` y todo el desplazamiento tangencial es deslizamiento, `g_s^p = g_s`. Al volver a cerrarse, la adherencia empieza desde la posición de recontacto (el «punto de adherencia» se reinicia). No es una regla aparte: sale del algoritmo de §8 con `μ·p_N = 0`.

Como la presión baja de forma continua hasta cero antes de abrirse, la tracción tangencial también llega a cero de forma continua (desliza mientras la presión cae): no hay saltos en `t`.

### 5. Energía

$$\psi = \tfrac12\,K_n\,\langle g_n\rangle_-^2 + \tfrac12\,K_t\,\|\mathbf g_s - \mathbf g_s^{p}\|^2, \qquad \mathcal D = \mathbf t_s\cdot\dot{\mathbf g}_s^{p} = \mu\,p_N\,\|\dot{\mathbf g}_s^{p}\| \ \ge 0$$

`ψ` es la energía almacenada en la penalización y la adherencia (recuperable); `𝒟` la potencia disipada por la fricción. En un ciclo tangencial de amplitud `±a` bajo presión constante, la energía disipada es `4·μ·p_N·(a − μ·p_N/K_t)`.

### 6. Variables internas

| Variable | Tipo | Significado |
|---|---|---|
| `slip` | vector, 1 (2D) o 2 (3D) componentes | deslizamiento por fricción `g_s^p` (confirmado) |
| `slip_accumulated` | escalar ≥ 0 | `∫‖ġ_s^p‖ dt` con la superficie cerrada: longitud recorrida deslizando en contacto (lo que se mueve abierta no cuenta) |
| `contact_pressure` | escalar ≥ 0 | `p_N` (informativa) |
| `status` | entero (guardado como real) | 0 abierta, 1 cerrada sin deslizar (adherencia; sin fricción, cerrada), 2 cerrada deslizando con fricción (informativa). Es un código (`STATE_CODES`, ADR 0043 §7): el VTK por celda no lo promedia; se lee por punto y en la historia |

`PRIMARY_STATE_VAR = 'contact_pressure'`.

### 7. Sistema local y notación

`g` y `t` en los ejes locales `(n, s)` o `(n, s1, s2)` de la interfaz, normales primero (ADR 0033); el paso a globales es del elemento. `n` apunta del lado A al lado B del par, de modo que `g_n > 0` es separación.

`t` es la **tracción sobre la cara A** (cuya normal exterior es `n`), conjugada de `g` en `δW = ∫δg·t dΓ`. Signos: `t_n > 0` sería tracción (nunca ocurre: `t_n ≤ 0`); en compresión `t = −p_N·n` sobre A, que coincide con la convención de la presión del ADR 0023 (`p > 0` compresión). `t_s = K_t(g_s − g_s^p)` tiene el **mismo sentido que el deslizamiento de B respecto de A** (B arrastra a A); la tracción sobre B es `−t` y es la que se opone a ese deslizamiento (revisión del teórico, 2026-11-26).

---

## Formulación numérica

### 8. Algoritmo (retorno radial cerrado, sin Newton local)

Dado `g_{n+1}` y el estado confirmado `g^p_{s,n}`:

1. Presión: `p_N = K_n·max(−g_n, 0)`.
2. Prueba elástica: `t_s^tr = K_t·(g_s − g^p_{s,n})`, `Φ^tr = ‖t_s^tr‖ − μ·p_N`.
3. Si `Φ^tr ≤ 0` → **adherencia**: `t_s = t_s^tr`, `g^p_s = g^p_{s,n}`.
4. Si `Φ^tr > 0` → **deslizamiento**: `m = t_s^tr/‖t_s^tr‖`, `t_s = μ·p_N·m`, `g^p_s = g^p_{s,n} + (Φ^tr/K_t)·m`.

Euler implícito con el límite `μ·p_N` evaluado en `n+1`. Las ecuaciones discretas de Euler implícito tienen **solución cerrada** porque el límite no depende del deslizamiento (Coulomb sin endurecimiento ni dilatancia): no hace falta Newton local. Lo exacto es esa solución discreta, no la integración de la ley en tasas: con trayectoria recta desde el reposo coincide con la integración fina (≤ 1e-15), pero si la dirección del deslizamiento gira dentro del paso (3D) o la presión varía, el error es de primer orden en el paso (*calculado* por el teórico: un giro de 90° en un paso da 8.7 % en `t`; con presión variable en 2D, 38 % en la disipación `μ·p_N·Δλ` del paso).

### 9. Tangente algorítmica consistente

`T = ∂t_{n+1}/∂g_{n+1}`, por ramas (`c = 1` si la superficie está cerrada, 0 si abierta; ver §11 para `g_n = 0`):

- **Abierta**: `T = 0`.
- **Cerrada, adherencia**: `T = diag(K_n, K_t·I)`.
- **Cerrada, deslizamiento, 3D**:

$$\frac{\partial t_n}{\partial g_n} = K_n, \qquad \frac{\partial \mathbf t_s}{\partial g_n} = -\mu K_n\,\mathbf m, \qquad \frac{\partial \mathbf t_s}{\partial \mathbf g_s} = \frac{\mu\,p_N\,K_t}{\|\mathbf t_s^{tr}\|}\,(\mathbf I - \mathbf m\otimes\mathbf m), \qquad \frac{\partial t_n}{\partial\mathbf g_s} = \mathbf 0$$

- **Cerrada, deslizamiento, 2D** (`m = ±1`): `∂t_s/∂g_n = −μ K_n m`, `∂t_s/∂g_s = 0`.

### 10. Simetría

- `μ = 0`: `T = diag(c·K_n, 0)`, simétrica.
- `μ > 0`: en deslizamiento, `∂t_s/∂g_n = −μK_n m ≠ ∂t_n/∂g_s = 0`: **no simétrica** (ley no asociada). `IS_SYMMETRIC` se declara por instancia: `True` sólo si `μ = 0` (como `IsotropicDamage2D`, que la fija según la medida de deformación). El despachador del ADR 0003 elige entonces LU. Simetrizar la tangente perdería la convergencia cuadrática del Newton: **decidido** usar la consistente no simétrica siempre que `μ > 0` (Diálogo Q3).

### 11. Decisión de régimen durante el Newton global

Abierta/cerrada y adherencia/deslizamiento se deciden en cada iteración desde la prueba (no desde el estado confirmado), como en la plasticidad. Las transiciones son no suaves: el Newton puede oscilar entre dos estados cerca de un cambio (abrir/cerrar, adherir/deslizar). Remedios estándar: búsqueda lineal, pasos más pequeños (bisección automática del paso).

**Holgura nula** (`g_n = 0`, el caso de partida de una interfaz de espesor nulo): `t_n = 0` en las dos ramas, pero la tangente no está definida. Se toma la rama **cerrada** (`g_n ≤ 0` ⇒ `c = 1`) para que la primera iteración tenga rigidez normal; con la rama abierta, un cuerpo que sólo se apoya en el contacto sería un mecanismo en la primera iteración y la matriz sería singular (decidido por la IA, Diálogo Q4: no cambia ningún resultado convergido). Con `g^p` sin cambios y `t_s^tr = 0` se toma adherencia (`Φ^tr = 0 ≤ 0`): rigidez lateral `K_t` en la primera iteración. Es una **rigidez de arranque**, no la derivada (con `p_N = 0` la derivada tangencial real es 0; la normal por la izquierda es `[K_n, −μK_n m]` y por la derecha 0, así que la rama cerrada es una elección lateral legítima); el Newton la corrige en la iteración siguiente.

### 12. Advertencias numéricas y límites del modelo

- **Penalización `K_n`**: la interpenetración es `p_N/K_n`; demasiado grande degrada el condicionamiento. Obligatorio, sin valor por omisión; guía del manual `K_n ≈ 10·E/h` (E y h del sólido vecino); con `K_n` fijo el error de la penalización no baja al refinar (Diálogo Q1, con los valores de los códigos generales).
- **Adherencia `K_t`**: el deslizamiento antes del umbral es `μ·p_N/K_t`. Obligatoria con `μ > 0`, constante, sin valor por omisión; guía del manual «empezar con `K_t = K_n`» (Diálogo Q2).
- Una superficie abierta no aporta rigidez: un cuerpo sostenido sólo por contactos que se abren es un mecanismo (el solver lo diagnostica como siempre).
- **Energías de la historia**: con fricción o con holgura inicial, `½U·F_int` deja de ser la energía almacenada (*medido* por el arquitecto: disipación de −0.5 en un resorte unilateral con `g_0 = 1`). La ley entrega `ψ` de §5 por `stored_energy(g, estado)` y el elemento suma `Σ w·ψ` (ADR 0043 §4–§5); `D = W − E` es entonces la energía disipada por fricción. `SECANT_UNLOADING = False`.
- **Ley unilateral** (`IS_UNILATERAL = True`, ADR 0019): una interfaz abierta en `U = 0` (holgura inicial positiva) no aporta rigidez; la red de seguridad avisa en lugar de rechazar el modelo, como con los cables (Diálogo Q6).
- **Sólo análisis no lineales**: en el análisis lineal y en el modal, un modelo con esta ley es error (decisión del usuario, Diálogo Q4 de `Interface2D4`). Implementado en todos los análisis de rigidez fija (también la perturbación sobre una etapa, el armónico, el espectral y los transitorios lineales): `ensure_no_unilateral_interfaces`, ADR 0043 «Consecuencias».
- Desplazamientos pequeños: los pares, la normal y los ejes tangentes son los de la configuración de referencia y no cambian durante el análisis. Válido mientras el deslizamiento sea pequeño frente al tamaño de los elementos de la interfaz.

---

## Contrato de implementación

```yaml
name: ContactPenalty
kind: interface_material        # familia de materiales de interfaz, ADR 0033
status: validated

interface:
  geometry: surface             # superficie de contacto: componentes (n, s) o (n, s1, s2)
  jump_dim: null                # 2 (2D) o 3 (3D)
  primary_state_var: contact_pressure
  is_symmetric: false           # por instancia: True sólo con mu = 0 (§10)
  accepts_initial_gap: true     # g_n = g_0 + Δu_n (ADR 0043 §4)
  is_unilateral: true           # sin rigidez abierta (ADR 0019)
  secant_unloading: false       # energía almacenada por stored_energy (§5, ADR 0043 §5)

admission:                      # criterio de admisión al programa principal (Reglas.md §1); evidencia en informes/bibliotecario/busqueda_contacto_2026-11-26.md
  textbooks:
    - zienkiewicz_taylor: "vol. 2, 7.ª ed., cap. 8: contacto nodo a nodo en desplazamientos pequeños §8.2 (pp. 237–242), penalización §8.2.2.3 (ec. 8.16), fricción de Coulomb §8.4.3 (pp. 257–263); juntas con Coulomb §4.11.2"
    - belytschko_liu_moran: "2.ª ed., cap. 10: Coulomb §10.3.2, adherencia elástica de asperezas (Curnier) §10.3.3, penalización en elastostática §10.5.6"
    - de_borst_crisfield: "§7.1, bloque con Coulomb y tangente no simétrica (ecs. 7.14–7.21)"
  general_codes:
    - abaqus: "GAPUNI (holgura inicial, contacto «hard»); fricción de Coulomb por penalización con deslizamiento elástico"
    - ansys: "CONTA178: nodo a nodo por penalización (FKN, FKS), holgura inicial GAP, Coulomb"
    - code_aster: "DEFI_CONTACT, método PENALISATION (E_N, E_T), U4.44.11 §3.5.1.1"
    - opensees: "zeroLengthContact2D/3D (Kn, Kt, μ)"

parameters:
  - { name: K_n, type: float, required: true,  desc: "Rigidez de penalización normal [Pa/m]. Guía en el Diálogo Q1" }
  - { name: mu,  type: float, required: false, desc: "Coeficiente de fricción de Coulomb (≥ 0). Sin él, contacto sin fricción (Diálogo Q5)" }
  - { name: K_t, type: float, required: false, desc: "Rigidez de adherencia tangencial [Pa/m]; obligatoria si mu > 0, se rechaza si mu = 0. Guía en el Diálogo Q2" }

signature:
  compute_traction: "(g: ndarray(2,) o (3,), state_vars=None) -> (t, T: ndarray(n,n), state_vars')"
  jump_kind: "sistema local: (g_n, g_s) o (g_n, g_s1, g_s2); g_n incluye la holgura inicial g_0 (la aporta el elemento)"

state_schema:
  slip: "ndarray(d−1) — deslizamiento por fricción g_s^p"
  slip_accumulated: "float ≥ 0 — ∫‖ġ_s^p‖"
  contact_pressure: "float ≥ 0 — p_N"
  status: "int — 0 abierto, 1 adherencia, 2 deslizamiento"

conventions:
  sign: "g_n > 0 separación, g_n < 0 interpenetración; t_n ≤ 0 (compresión); p_N = −t_n ≥ 0"
  frame: "ejes locales de la interfaz, n del lado A al lado B; transformación a globales por el elemento"
  units: "[K_n] = [K_t] = Pa/m, μ adimensional"

validity:
  - "K_n > 0; mu ≥ 0; K_t > 0 si mu > 0"
  - "salto de 2 o 3 componentes"
  - "desplazamientos y deslizamientos pequeños frente al tamaño de los elementos de la interfaz"

out_of_scope:
  - "Multiplicadores de Lagrange y Lagrangiano aumentado (imposición exacta): cambian el conjunto de restricciones entre iteraciones (pendientes del ADR 0040)"
  - "Deslizamiento grande, búsqueda de pares, contacto segmento a segmento o mortar"
  - "Fricción con coeficientes estático y dinámico distintos, dependiente de la velocidad o de la presión; fricción anisótropa"
  - "Dilatancia, adhesión, cohesión (para eso, las leyes cohesivas)"
  - "Contacto térmico (conductancia de la interfaz)"

numerical_caveats:
  - "Penalización K_n: interpenetración p_N/K_n; condicionamiento"
  - "Adherencia K_t: deslizamiento previo al umbral μ·p_N/K_t"
  - "Tangente no simétrica con fricción (LU)"
  - "Transiciones no suaves: posible oscilación del Newton entre estados"
  - "Holgura nula: rama cerrada (§11)"

acceptance:
  verification:
    - name: unilateral_normal
      setup: "g_n ∈ {−2, −1e-6, 0, 1e-6, 3}·(escala) con g_s = 0, K_n ≠ 1"
      expect: "t_n = K_n·g_n si g_n < 0; t = 0 si g_n > 0; t_n = 0 y T_nn = K_n en g_n = 0"
      tol_rel: 1.0e-14
    - name: adherencia_y_umbral_de_coulomb
      setup: "presión constante p_N (g_n < 0 fijo) y g_s monótono creciente; K_n, K_t, μ ≠ 1"
      expect: "t_s = K_t·g_s hasta ‖t_s‖ = μ·p_N; después ‖t_s‖ = μ·p_N y slip = g_s − μ·p_N/K_t"
      tol_rel: 1.0e-12
    - name: ciclo_de_friccion
      setup: "presión constante, g_s cíclico de amplitud ±a > μ·p_N/K_t, varios ciclos; medir desde el segundo ciclo (el primer cuarto disipa μ·p_N·(a − μ·p_N/K_t)); trabajo integrado con pasos que caen en los quiebres de la ley"
      expect: "energía disipada por ciclo estacionario 4·μ·p_N·(a − μ·p_N/K_t); trabajo = ψ + D; D ≥ 0"
      tol_rel: 1.0e-10
    - name: paso_unico_que_cruza_ramas
      setup: "un solo paso de adherencia a deslizamiento y de abierta a cerrada"
      expect: "t y estado iguales a la solución cerrada de §8"
      tol_rel: 1.0e-12
    - name: presion_decreciente_hasta_abrir
      setup: "deslizamiento sostenido mientras g_n crece hasta abrir"
      expect: "‖t_s‖ = μ·p_N baja de forma continua a 0; abierta, t = 0 y slip = g_s"
      tol_rel: 1.0e-12
    - name: energia_almacenada
      setup: "stored_energy frente al trabajo integrado exactamente en trayectorias sin deslizamiento; holgura g_0 > 0 en U = 0"
      expect: "ψ = trabajo; ψ(g_0 > 0, U = 0) = 0"
      tol_rel: 1.0e-12
    - name: separacion_y_recontacto
      setup: "cerrar con deslizamiento, abrir (g_n > 0), desplazar tangencialmente, volver a cerrar"
      expect: "t = 0 abierta; al cerrar, t_s = 0 en la posición de recontacto (punto de adherencia reiniciado)"
      tol_rel: 1.0e-12
  specific:
    - name: isotropia_tangencial_3d
      setup: "trayectoria circular de g_s bajo presión constante; rotación de los ejes (s1, s2) del plano tangente"
      expect: "‖t_s‖ = μ·p_N en deslizamiento con t_s paralelo a la prueba; resultado invariante ante la rotación de (s1, s2)"
      tol_rel: 1.0e-12
    - name: sin_friccion
      setup: "mu = 0 (o ausente), g_s arbitrario"
      expect: "t_s = 0, T = diag(c·K_n, 0); IS_SYMMETRIC True; K_t rechazado"
      tol_rel: 1.0e-14
    - name: tangente_diferencia_finita_consistente
      setup: "puntos en cada rama (abierta, adherencia, deslizamiento) en 2D y 3D, lejos de las transiciones"
      expect: "T = diferencias centradas de t"
      tol_rel: 1.0e-6
    - name: cambio_de_direccion_3d
      setup: "3D: deslizar en s1 y después en s2 (90°); subdividir el segundo tramo en 1, 2, 4… pasos"
      expect: "la solución converge a la integración fina con orden 1"
      tol_rel: 1.0e-2
    - name: prueba_y_confirmacion
      setup: "repetir la misma llamada desde el mismo estado confirmado"
      expect: "mismo t, misma tangente, mismo estado de prueba (status, slip_accumulated incluidos)"
      tol_rel: 0.0
    - name: escalado_con_los_parametros
      setup: "multiplicar K_n, K_t por α y las holguras por 1/α; multiplicar μ por β"
      expect: "t invariante con α; t_s en deslizamiento escala con β (Reglas §6: parámetros ≠ 1)"
      tol_rel: 1.0e-12
  arch:
    - name: registro_de_interfaz
      expect: "Registro vía InterfaceMaterialRegistry; sección YAML interface_materials"
    - name: is_symmetric_por_instancia
      expect: "ContactPenalty(K_n=…).IS_SYMMETRIC is True; con mu > 0, False"

references:
  - "Zienkiewicz, O. C., Taylor, R. L., Fox, D. D. (2014). The Finite Element Method for Solid and Structural Mechanics, 7.ª ed., cap. 8 (Treatment of Constraints: Contact and Tied Interfaces). Butterworth-Heinemann. DOI 10.1016/B978-1-85617-634-7.00008-9"
  - "Belytschko, T., Liu, W. K., Moran, B., Elkhodary, K. I. (2014). Nonlinear Finite Elements for Continua and Structures, 2.ª ed., cap. 10. Wiley."
  - "de Borst, R., Crisfield, M. A., Remmers, J. J. C., Verhoosel, C. V. (2012). Non-Linear Finite Element Analysis of Solids and Structures, 2.ª ed., §7.1. Wiley. DOI 10.1002/9781118375938"
  - "Wriggers, P. (2006). Computational Contact Mechanics, 2.ª ed. Springer. DOI 10.1007/978-3-540-32609-0"
  - "Curnier, A. (1984). A theory of friction. Int. J. Solids Struct. 20(7), 637–647. DOI 10.1016/0020-7683(84)90021-0"
  - "Laursen, T. A., Simo, J. C. (1993). Algorithmic symmetrization of Coulomb frictional problems using augmented Lagrangians. CMAME 108, 133–146. DOI 10.1016/0045-7825(93)90157-S"
  - "ADR 0033 — Familia de materiales de interfaz; ADR 0043 — Elementos de interfaz de espesor nulo."
```

---

## Implementación

- Archivo: `solidum/materials/interface/contact_penalty.py`; contrato de la familia ampliado en `solidum/core/interface_material.py` (`IS_UNILATERAL`, `ACCEPTS_INITIAL_GAP`, `STATE_SCHEMA`, `STATE_CODES`, `stored_energy`; ADR 0043 §4).
- Clase: `ContactPenalty` (`IS_SYMMETRIC` por instancia, `True` sólo con `mu = 0`; `STATE_CODES = ("status",)`).
- Tests:
  - `tests/test_contact_penalty.py`: todos los casos de `acceptance`; 12/12 mutaciones detectadas (`informes/prototipos_2026-11-26/contacto/mutar_contacto.py`).
  - `tests/test_interface_material_contract.py`: `STATE_SCHEMA` de cada ley de interfaz frente al estado que crea.
  - `tests/test_contact_requires_nonlinear.py`: los análisis de rigidez fija rechazan la ley.
  - Elemento y energía: `tests/test_interface_elements.py`, `test_interface_energy.py`, `test_interface_output.py`.
- Verificación: `informes/teorico/verificacion_teorica_2026-11-26_contacto.md` (ley y tangente frente a una referencia independiente en 16 000 puntos, diferencia 0; *calculado*).
- Notas de traducción: ninguna; la convención es la del proyecto (`t` sobre la cara A, `p > 0` compresión, ADR 0023).

---

## Diálogo

*Preguntas, aclaraciones y hallazgos. Entradas fechadas.*

- **2026-11-26 — Decisiones abiertas para el usuario** (se preguntan una a la vez, después de la revisión de `teorico-solidum`):
  - **Q1 — DECIDIDO por el usuario el 2026-11-26 («si»): `K_n` obligatorio, sin valor por omisión, con la guía `K_n ≈ 10·E/h` (E y h del sólido vecino) en el manual.** Razón dada: igual que `K_e` de las cohesivas, el valor que fija el error queda a la vista, no exige reglas nuevas (qué E con dos materiales, con daño o con un lado casi rígido) y deja abierto añadir después un cálculo automático cuando no se declare, sin cambiar modelos existentes. Con K_n fijo, el error de la penalización no baja al refinar (en Hertz ≈ 1 % con h = a/10). Pregunta: `K_n`: ¿obligatorio sin valor por omisión (como `K_e` de las cohesivas), con una guía en el manual? Lo que hacen los códigos (*leído*, informe del bibliotecario): ANSYS, factor FKN = 1 (rango 0.1–10) sobre una rigidez que calcula del módulo y del tamaño del elemento; Abaqus (v2025), 10 veces la rigidez representativa del elemento subyacente; Code_Aster, sin valor por omisión («del orden del mayor módulo de Young» por una longitud, y subir hasta resultados estables; en juntas `E/e`); OpenSees y DIANA, lo da el usuario. Recomendación de la IA: obligatorio, con la guía `K_n ≈ 10·E/h` (E y h del sólido adyacente), como `K_e`.
  - **Q2 — DECIDIDO por el usuario el 2026-11-26 («1»): `K_t` obligatorio con `μ > 0`, constante, sin valor por omisión; guía en el manual «empezar con `K_t = K_n`»** (que `μ·p/K_t` sea pequeño frente al deslizamiento esperado). Pregunta y datos: Abaqus fija el deslizamiento elástico admisible (0.005 de una longitud característica), con lo que su `K_t` efectivo crece con la presión; ANSYS lo hace proporcional a `μ·K_n`; Code_Aster (juntas), `K_t = K_n` por omisión; OpenSees, del usuario. Recomendación de la IA: obligatorio con `μ > 0`, constante (Code_Aster, OpenSees, Goodman), con una guía en el manual.
  - **Q3 — DECIDIDO por el usuario el 2026-11-26 («si»): tangente consistente no simétrica siempre que `μ > 0`.** Pregunta: tangente con fricción: ¿consistente no simétrica (convergencia cuadrática, LU) o simetrizada? Abaqus usa la no simétrica automáticamente con `μ > 0.2`; ANSYS la recomienda (a veces obligatoria) con `μ > 0.2`; OpenSees y DIANA, no simétrica. Recomendación de la IA: consistente no simétrica siempre que `μ > 0`.
  - **Q4 — Decidido por la IA (2026-11-26): rama cerrada en `g_n = 0`.** No cambia ningún resultado convergido (la fuerza es 0 en las dos ramas), sólo la rigidez de arranque: es numérica (Reglas §3), se le dijo al usuario en una línea. Pregunta original: holgura nula: ¿rama cerrada o abierta?
  - **Q5 — DECIDIDO por el usuario el 2026-11-26 («1»): sin `mu` declarado, contacto sin fricción (`μ = 0`)**, como Abaqus y ANSYS. Pregunta: sin `mu` declarado, ¿contacto sin fricción o `mu` obligatorio?
  - **Q6 — Resuelto por la regla existente (2026-11-26): aviso.** El ADR 0019 ya fija que con materiales unilaterales un pivote nulo en la tangente inicial es aviso (si no se estabiliza, el Newton falla después como tangente singular); la ley declara `IS_UNILATERAL` y se aplica esa regla, como con los cables. Se le dijo al usuario en una línea.
  - Integración del elemento (nodal o de Gauss), definición de `g_0` y de la normal: en la spec del elemento.
- **2026-11-26 — Revisión de la spec por `teorico-solidum`** (`informes/teorico/revision_spec_contacto_2026-11-26.md`, veredicto «cambios antes de programar»): ley, tangente (diferencias centradas en las seis ramas, ≤ 2.7e-12), energía y disipación por ciclo correctas (*calculado*). Aplicado: sentido de `t_s` (§7), «solución cerrada de Euler implícito» en lugar de «exacto» (§8), rigidez de arranque en `g_n = 0` (§11), pruebas nuevas y ajuste de `ciclo_de_friccion`. Opinión del teórico en las decisiones abiertas: coincide con las recomendaciones de la IA en Q1–Q6. Hallazgo para el elemento y la salida: `D = W − E` con el trapecio del ADR 0038 da una disipación espuria cuando un paso cruza un quiebre de la ley (cierre de holgura en un paso: `D/ψ = 0.50`; en 10 pasos, 5e-3; *calculado*); es del integrador del trabajo, no de la ley. Su signo depende del quiebre (verificación del 2026-11-26): positiva al cerrar una holgura, negativa al pasar de adherencia a deslizamiento (la disipación de ese paso sale la mitad de la real).
- **2026-11-26 — Revisión del plan por `arquitecto-solidum`** (`informes/arquitecto/revision_plan_contacto_2026-11-26.md`): la holgura inicial la calcula el elemento y la ley la recibe dentro de `g_n`; contrato ampliado de la familia (`ACCEPTS_INITIAL_GAP`, `IS_UNILATERAL`, `STATE_SCHEMA`, `stored_energy`) y energía almacenada por elemento: ADR 0043. La historia sólo exporta variables escalares: `slip` (vector) queda en el estado y se exportan `contact_pressure`, `slip_accumulated` y `status`.
- **2026-11-27 — Validación (`validador-solidum`)**: escalera de la spec `Interface2D4` §7b en 2D (parcela, despegue, fricción, Hertz) y con etapas (despegue, fricción 2D de ida y vuelta, fricción 3D con trayectoria en L, con `Interface3D8`); resultados en `Interface2D4` §7b; `status: validated`. Viga SEN con la placa apoyada: funciona sin fricción; con μ = 0.5 y la placa libre en horizontal, el análisis se detiene en la cola. Es un problema de cómo se representa la placa, no de la ley (*hipótesis*); queda abierto en el Diálogo de `Interface2D4`.
