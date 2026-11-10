# BondSlipModelCode — adherencia-deslizamiento de barras corrugadas (Model Code 2010 / Eligehausen et al. 1983)

> Orden de trabajo. El usuario escribe la **especificación física**, la **formulación numérica** y el **contrato**. La IA rellena **implementación** y responde en **diálogo** durante el trabajo.
>
> Material de interfaz de **línea** (ADR 0033): relaciona la tracción con el desplazamiento relativo entre una barra embebida y el sólido que la rodea, en el sistema local de la barra: componentes normales primero y deslizamiento después, `(n, s)` en 2D y `(n1, n2, s)` en 3D.

---

## Especificación física

### 0. Descripción general
Ley local de adherencia-deslizamiento de una barra corrugada en concreto, con tres partes. Todas son decisiones del usuario del 2026-11-10:
- **Tangencial (deslizamiento `s`)**: envolvente de Eligehausen, Popov y Bertero (1983, §5.3.1), adoptada por el Model Code 2010 (Tabla 6.1-1). Rama ascendente potencial, meseta, descenso lineal y rama residual de fricción. Antes de la rama potencial hay un tramo lineal con rigidez inicial `k0`.
- **Normal**: elástica lineal y desacoplada del deslizamiento (como DIANA, §44.5).
- **Descarga y recarga**: geometría de Eligehausen et al. (1983, Fig. 5.3b), sin degradación y con fricción constante igual a la residual `τ_f` (regla de DIANA, R1). La degradación por energía de Eligehausen et al. (R3) queda para un módulo de usuario.

Es una ley **local**: no depende de la deformación de la barra ni de la presión transversal. Esos efectos, que el Model Code trata en §6.1.1.3 y Murcia-Delso y Shing en su ley de 2015, quedan fuera.

### 1. Componentes
- `s` — deslizamiento: desplazamiento de la barra respecto del concreto a lo largo de su eje `t`; `s > 0` cuando la barra avanza en `+t` respecto del concreto.
- `Δu_n` (`Δu_n1`, `Δu_n2` en 3D) — desplazamiento relativo normal al eje.
- Tracciones: `τ` (esfuerzo de adherencia, sobre el perímetro de la barra) y `t_n` (sobre su ancho proyectado). El ancho sobre el que actúa cada una lo pone el elemento (spec `EmbeddedReinforcement`).

### 2. Envolvente monótona (Eligehausen et al. 1983, §5.3.1; MC2010, ecs. 6.1-1 a 6.1-4)
Con `a = |s|`:

$$
\tau_{env}(a) = \begin{cases}
k_0\,a & 0 \le a \le s_0 \\
\tau_{max}\,(a/s_1)^{\alpha} & s_0 < a \le s_1 \\
\tau_{max} & s_1 < a \le s_2 \\
\tau_{max} - (\tau_{max} - \tau_f)\,\dfrac{a - s_2}{s_3 - s_2} & s_2 < a \le s_3 \\
\tau_f & a > s_3
\end{cases}
\qquad
s_0 = s_1\left(\frac{k_0\,s_1}{\tau_{max}}\right)^{1/(\alpha - 1)}
$$

`s0` es el punto donde la recta `k0·a` corta la curva potencial (continuidad de τ). Con α < 1, la rigidez de la rama potencial en `a → 0` es infinita; el tramo lineal la acota (decisión del usuario del 2026-11-10, pregunta 3). Se exige `k0 > τ_max/s1` para que el corte ocurra antes de `s1`.

La misma ley vale si la barra se jala o se empuja (Eligehausen et al. 1983, §5.3.1): `τ(−s) = −τ(s)` en carga monótona.

### 3. Componente normal (DIANA §44.5)
$$t_{n} = k_n\,\Delta u_n$$
en cada componente normal, en tracción y en compresión por igual, independiente de `s`. Una barra recta sólo transmite fuerza normal en los quiebres de la polilínea (`N·Δθ`). La presión radial de cuña es autoequilibrada en el perímetro: el hendimiento entra por los parámetros de la envolvente (columnas de hendimiento del Model Code), no por `t_n`.

### 4. Descarga y recarga (R1: Eligehausen et al. 1983, Fig. 5.3b, sin degradación)
Variables de historia: el deslizamiento extremo de cada lado y su tracción, `(s⁺, τ_E⁺)` con `s⁺ ≥ 0` y `(s⁻, τ_E⁻)` con `s⁻ ≤ 0`, ambos nulos en estado virgen.
- Recta de descarga y recarga de cada lado: `ℓ⁺(s) = τ_E⁺ + k0 (s − s⁺)` (sólo si `s⁺ > 0`) y `ℓ⁻(s) = τ_E⁻ + k0 (s − s⁻)` (sólo si `s⁻ < 0`).
- Fricción: `τ_fr = min(τ_f, τ_env(max(s⁺, −s⁻)))`. Es nula en estado virgen y nunca supera la mayor tracción alcanzada.
- Cotas:
$$U(s) = \max\!\big(\tau_{fr},\ \min(\tau_{env}(s)\,[s > 0],\ \ell^+(s))\big), \qquad
L(s) = -\max\!\big(\tau_{fr},\ \min(\tau_{env}(-s)\,[s < 0],\ -\ell^-(s))\big)$$
- Paso: predictor `τ* = τ_n + k0 (s − s_n)` desde el estado confirmado y `τ = min(max(τ*, L(s)), U(s))`. Si `s` supera un extremo, se actualiza ese extremo con la tracción resultante.

El recorrido de la Fig. 5.3b sale de estas reglas:
- **Descarga desde E**: rama rígida con `k0` hasta `−τ_f` (camino EFG).
- **Fricción**: a `−τ_f` hasta cortar la envolvente del lado opuesto, que parte del origen (camino GHI). Después sigue esa envolvente.
- **Al invertir otra vez**: rama rígida hasta `+τ_f`, y fricción hasta cortar la recta de descarga anterior (punto N).
- **Recarga**: con `k0` por esa recta hasta la envolvente.
- **Ciclos completos entre extremos fijos**: lazos rígido-plásticos sin degradación.

### 5. Variables internas
`s_pos`, `tau_pos` (extremo positivo y su tracción), `s_neg`, `tau_neg` (extremo negativo), `s`, `tau` (último estado confirmado, para el predictor), `slip_max = max(s⁺, −s⁻)` (`PRIMARY_STATE_VAR`, para el post-proceso).

---

## Formulación numérica

### 6. Tangente
`∂τ/∂s = k0` si `L < τ* < U`. Si una cota está activa, vale la pendiente de esa cota: 0 en la fricción, `τ_env'(|s|)` en la envolvente (con su signo en el ablandamiento, negativo, y 0 en la meseta y la residual), `k0` en una recta de descarga. Normal: `∂t_n/∂Δu_n = k_n`. Sin acoplamiento: la tangente es diagonal y **simétrica** (`IS_SYMMETRIC = True`).

### 7. Función auxiliar del Model Code 2010
`BondSlipModelCode.model_code_2010(fcm, bond, failure, ...)` construye la ley con los parámetros de la Tabla 6.1-1. `fcm` es la resistencia media a compresión en cilindro, en MPa en las fórmulas, y `c` la separación libre entre corrugas.

| | Arrancamiento, buena | Arrancamiento, otras | Hendimiento, buena, sin estribos | Hendimiento, buena, con estribos | Hendimiento, otras, sin estribos | Hendimiento, otras, con estribos |
|---|---|---|---|---|---|---|
| τ_max (MPa) | 2.5√fcm | 1.25√fcm | 7.0(fcm/25)^0.25 | 8.0(fcm/25)^0.25 | 5.0(fcm/25)^0.25 | 5.5(fcm/25)^0.25 |
| s1 (mm) | 1.0 | 1.8 | s(τ_max) | s(τ_max) | s(τ_max) | s(τ_max) |
| s2 (mm) | 2.0 | 3.6 | s1 | s1 | s1 | s1 |
| s3 | c | c | 1.2·s1 | 0.5·c | 1.2·s1 | 0.5·c |
| α | 0.4 | 0.4 | 0.4 | 0.4 | 0.4 | 0.4 |
| τ_f | 0.40τ_max | 0.40τ_max | 0 | 0.40τ_max | 0 | 0.40τ_max |

- **`s(τ_max)` en hendimiento** es el deslizamiento en que la curva de arrancamiento de la misma condición alcanza el `τ_max` de hendimiento: `s1 = s1_PO·(τ_max,SP/τ_max,PO)^(1/α)`.
- **Valores por omisión** (decisiones del usuario del 2026-11-10): `k0 = 180 N/mm³` = 1.8·10¹¹ Pa/m (Eligehausen et al. 1983, §5.3.4) y `k_n = 2·E_c/d_b` (Murcia-Delso y Shing 2016). Este último exige `E_c` y `d_b` si no se da `k_n`.
- **Fuentes de la tabla**: secundarias, sin lectura del texto del Model Code. Son Bastos (2022), Tabla 8, y la librería abierta de la fib `structuralcodes`, que coinciden. Discrepancia sin resolver: Tang y Cheng (2020) dan s2 = 3.0 mm para arrancamiento con buena adherencia; se adopta 2.0 mm, de las fuentes que coinciden. Ver el Diálogo.
- **YAML**: `type: BondSlipModelCode` con los parámetros explícitos, o con un bloque `model_code_2010: {fcm, bond, failure, stirrups, clear_rib_spacing, ...}` que llama a la función auxiliar.

---

## Contrato de implementación

```yaml
name: BondSlipModelCode
kind: interface_material
status: validated

interface:
  geometry: line                 # barra embebida: (n, s) en 2D, (n1, n2, s) en 3D
  jump_dim: null                 # 2 en 2D y 3 en 3D
  primary_state_var: slip_max
  is_symmetric: true

admission:
  textbooks:
    - de_borst_crisfield: "§13.1 elementos de interfaz (marco); §6.8.2 p. 207 sólo menciona la adherencia-deslizamiento: el criterio no decide y lo decidió el usuario"
  general_codes:
    - diana: "BONDSL 6, ley fib Model Code 2010 con tramo lineal inicial, normal lineal y descarga hasta ±τ_f"

parameters:
  - { name: tau_max, type: float, required: true,  desc: "resistencia máxima de adherencia τ_max [Pa]" }
  - { name: s1,      type: float, required: true,  desc: "deslizamiento al alcanzar τ_max [m]" }
  - { name: s2,      type: float, required: true,  desc: "fin de la meseta [m], s2 ≥ s1" }
  - { name: s3,      type: float, required: true,  desc: "inicio de la rama residual [m], s3 ≥ s2" }
  - { name: alpha,   type: float, required: true,  desc: "exponente de la rama ascendente, 0 < α < 1" }
  - { name: tau_f,   type: float, required: true,  desc: "resistencia residual de fricción [Pa], 0 ≤ τ_f ≤ τ_max" }
  - { name: k0,      type: float, required: true,  desc: "rigidez inicial y de descarga [Pa/m], k0 > τ_max/s1" }
  - { name: k_n,     type: float, required: true,  desc: "rigidez normal [Pa/m], > 0" }

conventions:
  sign: "s > 0 ⇔ la barra avanza en +t respecto del concreto; τ con el signo de s en carga monótona"
  frame: "local de la barra, normales primero; el paso a globales y los anchos son del elemento"
  units: "τ, t_n en Pa; s, Δu_n en m; k0, k_n en Pa/m (SI; el sistema no convierte unidades)"

validity:
  - barras corrugadas en concreto confinado (arrancamiento) o con los parámetros de hendimiento del Model Code
  - acero sin fluir (la Tabla 6.1-1 está limitada a ε_s < ε_y)
  - carga monótona o con pocos ciclos (sin degradación cíclica)

out_of_scope:
  - degradación cíclica de la envolvente y de la fricción (Eligehausen et al. 1983, R3) → módulo de usuario
  - efecto de la fluencia de la barra y de la presión transversal sobre τ (MC2010 §6.1.1.3; Murcia-Delso y Shing 2015)
  - acoplamiento normal-tangencial, dilatancia, separación por hendimiento (Murcia-Delso y Shing 2015/2016)
  - fatiga, fluencia lenta, velocidad de carga

acceptance:
  verification:
    - name: envolvente
      setup: "carga monótona en s > 0 y s < 0 con los valores de Eligehausen et al. (τ1 = 13.5 MPa, s1 = 1, s2 = 3, s3 = 10.5 mm, α = 0.4, τ_f = 5 MPa, k0 = 180 N/mm³)"
      expect: "τ igual a §2 en puntos de cada rama y τ(−s) = −τ(s); continuidad en s0, s1, s2 y s3"
      tol_rel: 1.0e-12
    - name: ciclo_figura_5_3b
      setup: "historia 0 → 1.0 → −1.0 → 1.5 mm con los mismos valores"
      expect: "descarga a −τ_f en s = 1 − (13.5 + 5)/180 = 0.8972 mm; fricción hasta la envolvente negativa en |s| = (5/13.5)^2.5 = 0.0835 mm; envolvente hasta −13.5 MPa en −1 mm; rama rígida a +5 MPa en −0.8972 mm; fricción hasta el punto N, s = 1 − 8.5/180 = 0.9528 mm; recarga con k0 hasta 13.5 MPa en 1.0 mm y envolvente hasta 1.5 mm"
      tol_rel: 1.0e-12
    - name: tangente_diferencias_finitas
      setup: "estados en cada rama (lineal inicial, potencial, meseta, ablandamiento, residual, descarga, fricción, recarga), 2D y 3D"
      expect: "tangente consistente con diferencias finitas centradas"
      tol_rel: 1.0e-6
    - name: disipacion_no_negativa
      setup: "1000 historias aleatorias de deslizamiento"
      expect: "trabajo acumulado ∫τ ds − τ²/(2·k0) ≥ 0 en todo instante"
      tol_abs: 1.0e-9
    - name: normal_desacoplada
      setup: "Δu_n arbitrario con cualquier s e historia"
      expect: "t_n = k_n·Δu_n, sin término cruzado en la tangente"
      tol_rel: 1.0e-14
  specific:
    - name: model_code_2010
      setup: "función auxiliar con fcm = 30 MPa, c = 10 mm, d_b = 16 mm, E_c = 30 GPa en las seis columnas"
      expect: "parámetros de la tabla de §7; k0 = 1.8e11 Pa/m y k_n = 2·E_c/d_b por omisión"
      tol_rel: 1.0e-12
    - name: rechazos
      setup: "k0 ≤ τ_max/s1, s2 < s1, s3 < s2, α ∉ (0, 1), τ_f > τ_max, k_n ≤ 0"
      expect: "ValueError al construir con mensaje que nombra el parámetro"
  arch:
    - name: registro_de_interfaz
      expect: "en InterfaceMaterialRegistry; GEOMETRY = 'line'; accepts('line', 2) y accepts('line', 3); no accepts('surface', 2)"

references:
  - "Eligehausen R., Popov E. P., Bertero V. V. (1983). Local bond stress-slip relationships of deformed bars under generalized excitations. Report UCB/EERC-83/23, University of California, Berkeley. §5.3 (pp. 56–62), Fig. 5.3."
  - "fib (2013). fib Model Code for Concrete Structures 2010. Wiley. §6.1.1, Tabla 6.1-1 (leída por fuentes secundarias: Bastos 2022; fib structuralcodes)."
  - "Murcia-Delso J., Shing P. B. (2016). Elastoplastic dilatant interface model for cyclic bond-slip behavior of reinforcing bars. J. Eng. Mech. 142(2), 04015082 (k_n = 2·E_c/d_b)."
  - "DIANA FEA, Material Library §26.3 / Theory §44.5, bond-slip interface (normal lineal, descarga)."
```

---

## Implementación

- **Archivo**: [solidum/materials/interface/bond_slip.py](../../solidum/materials/interface/bond_slip.py) · clase `BondSlipModelCode` (registrada en `InterfaceMaterialRegistry`), sobre la base común `_BondSlipLaw`, que implementa la regla R1 (§4); cada ley aporta sólo `envelope(a) → (τ, dτ/da)`.
- **Función auxiliar**: `BondSlipModelCode.model_code_2010(fcm, bond=, failure=, stirrups=, clear_rib_spacing=, diameter=, E_c=, k0=, k_n=)`, en SI. En el YAML, el bloque `model_code_2010: {...}` (declarado en `YAML_KEYS` y traducido por `from_yaml`).
- **Tests**: [tests/test_bond_slip_laws.py](../../tests/test_bond_slip_laws.py): envolvente y simetría, ciclo de la Fig. 5.3b con los valores exactos de `acceptance`, tangente por diferencias finitas en todas las ramas (2D y 3D), tangente k0 en el estado virgen, disipación no negativa en 1000 historias aleatorias, normal desacoplada, las seis columnas del Model Code 2010, rechazos, registro y geometría. El uso dentro del refuerzo, en [tests/test_embedded_bond_slip.py](../../tests/test_embedded_bond_slip.py).
- **Notas de traducción**:
  - Las cotas `U`, `L` actúan sólo si el predictor las **rebasa**, y la fricción sólo si supera **estrictamente** la envolvente o la recta. En el estado virgen con `s = 0` las dos cotas y la fricción valen cero; con un empate resuelto hacia la fricción, la tangente salía nula y la interfaz no tenía rigidez al empezar. Lo detectó la red de seguridad del ADR 0019 como mecanismo y lo cubre `test_tangente_inicial_es_k0`.
  - `density = 0.0` (clase base de la familia): las leyes de interfaz no tienen masa (ADR 0008).

---

## Diálogo

- **2026-11-10** · Decisiones del usuario, con los datos en `ESTADO.md` y en los informes `busqueda_leyes_adherencia_2026-11-10.md` y `busqueda_refuerzo_acoplamiento_2026-11-10.md`:
  - (2) Forma Model Code / Eligehausen con parámetros explícitos, función auxiliar del Model Code 2010 y ley multilineal; Murcia-Delso y Shing, más adelante, como módulo de usuario.
  - (3) Rigidez inicial `k0`, 180 N/mm³ por omisión. Con k0 ≥ 180 el resultado no cambia (±0.01 %); con 83/55/36 N/mm³ el deslizamiento al 30 % de la carga sube 5/19/52 %.
  - (4) Normal elástica lineal desacoplada, sobre `d_b`, con `k_n = 2·E_c/d_b` por omisión.
  - (5) Descarga R1, sin degradación. En el 2.º ciclo a ±0.25/±1/±3 mm disipa 4.1/18.5/58.5 N/mm, frente a 0.6/4.0/13.5 de Eligehausen completo (R3).
- **2026-11-10** · Clasificación: el criterio de admisión no decide la adherencia por los textos. De Borst et al. 2012 §6.8.2, p. 207, sólo menciona el «bond-slip behaviour», y Paulay–Priestley la trata como criterio de diseño. Los códigos sí la tienen: DIANA y Code_Aster. El usuario decidió ponerla en el programa principal (2026-11-10).
- **2026-11-10** · Discrepancia de la Tabla 6.1-1 sin resolver: s2 = 2.0 mm (Bastos 2022, `structuralcodes`, Croce et al. 2022) frente a 3.0 mm (Tang y Cheng 2020, que coincide con Eligehausen et al. 1983). Se adopta 2.0 mm. Conviene confirmarlo con el texto del Model Code cuando esté disponible.
- **2026-11-10** · Spec **aprobada por el usuario** («si»).
- **2026-11-10** · Implementada y validada: todos los casos de `acceptance` en verde. Promovido `status: validated`.
