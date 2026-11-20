# CST_Embedded2D — triángulo CST 2D con discontinuidad interior embebida (KOS)

> Orden de trabajo. Pre-redactada por la IA como borrador para validación del usuario (autor de la formulación, Retama 2010). La cinemática KOS y la fórmula `l_d = (A/h)·cos(θ−α)` son traducción directa de los Caps. 2, 5, 6 y 7 de la tesis; la mecánica de implementación (Voigt, registry, signatura, contrato `Element`, condensación) sigue las convenciones del proyecto. El usuario revisa, corrige donde proceda, y aprueba antes de que la IA toque código.

---

## Especificación física

### 0. Descripción general

Elemento finito sólido 2D triangular de tres nodos (CST padre) **enriquecido** con una **discontinuidad interior embebida** `Γ_d` que cruza el dominio del elemento cuando se cumple un criterio de activación. Es la primera implementación de la familia de elementos con DOFs enriquecidos elementales introducida por **ADR 0010**, y la **fase 2** de la hoja de ruta de discontinuidades embebidas.

Estructura del elemento:

- **Estado intacto** (`discontinuity_state is None`): se comporta exactamente como un [`Tri3`](../../../specs/Tri3.md) estándar — CST con un único punto de Gauss, B constante, integración exacta del bulk con peso `A_e·t`.
- **Estado agrietado** (`discontinuity_state is not None`): el campo de desplazamientos del bulk se enriquece con un salto `[[u]] ∈ ℝ²` localizado en `Γ_d`. El bulk vecino descarga elásticamente; la disipación se concentra en `Γ_d` y la gobierna un material cohesivo (material de interfaz de superficie del programa principal, `InterfaceMaterial`, ADR 0033; antes `CohesiveMaterial`, ADR 0010 §1).

Aproximación adoptada: **discrete approach** (Retama 2010), **KOS** (Kinematically Optimal Symmetric) según la clasificación de Jirásek (2000). La normal `n` se fija al activar el elemento (perpendicular a `σ_I`, criterio de Rankine) y no se reorienta después. **El nodo solitario lo dicta la trayectoria de la grieta**, coherente entre vecinos, mediante el seguimiento de grieta de la spec [`CrackTracking`](CrackTracking.md) (decisiones del usuario del 2026-11-10; desde el 2026-11-19, prolongación recta del último tramo de la grieta). Antes lo elegía cada elemento con la recta por su centroide, y en mallas no alineadas eso bloqueaba la grieta. Modo-I puro (consistente con `CohesiveDamageIsotropic` fase 1); modo mixto y KSON quedan diferidos (fases G y H del ADR 0010).

Régimen: pequeños desplazamientos y deformaciones, cuasi-estático. Pensado para hormigón, roca, cerámicos cuasi-frágiles. Bulk **elástico** por construcción de la formulación (la disipación está localizada en `Γ_d`); materiales no lineales del bulk quedan `out_of_scope` en fase 1.

### 1. Cinemática de desplazamientos (formulación KOS)

Campo total de desplazamientos en el elemento:

$$\mathbf u(\mathbf x) = \underbrace{\sum_{i=1}^{3} N_i(\mathbf x)\,\mathbf d_i}_{\mathbf u^{\text{std}}(\mathbf x)} \;+\; \underbrace{\bigl(H_{\Gamma_d}(\mathbf x) - \varphi(\mathbf x)\bigr)\,\mathbf G \,\llbracket\mathbf u\rrbracket}_{\mathbf u^{\text{enr}}(\mathbf x)}$$

donde:

- `N_i(x)` son las funciones de forma estándar del CST (coordenadas baricéntricas; ver [`Tri3`](../../../specs/Tri3.md) §7).
- `d_i ∈ ℝ²` son los desplazamientos nodales estándar (DOFs globales del modelo, 2 por nodo).
- `[[u]] ∈ ℝ²` es el salto de desplazamientos a través de `Γ_d`, expresado en el **frame local** `(n, s)` de la discontinuidad. Es un **DOF enriquecido elemental** — vive en el elemento, no se ensambla al sistema global (ADR 0010 §2).
- `H_{Γ_d}(x)` es la función escalón de Heaviside relativa a `Γ_d` (0 en `Ω⁻`, 1 en `Ω⁺`). `Ω⁺` es el subdominio del lado al que apunta `n`.
- `φ(x)` es una función de soporte continua que vale `1` en los nodos de `Ω⁺` y `0` en los de `Ω⁻`. Para CST con `Γ_d` cortando dos lados, `φ(x)` es la combinación lineal `Σ N_i^+(x)` sobre los nodos del lado positivo (típicamente un único nodo solitario, ver §6).
- `G` es la matriz que transforma el salto del frame local al sistema cartesiano global:

$$\mathbf G = \begin{bmatrix} n_x & s_x \\ n_y & s_y \end{bmatrix}$$

La construcción `(H_{Γ_d} - φ)` garantiza que `u^enr` sea **continuo en el contorno del elemento** (porque `H_{Γ_d} = φ` en los nodos) y que el salto correcto `[[u]] = G^T·[[u]]_global` ocurra únicamente sobre `Γ_d` interior. Es la elección **kinematically optimal**: el campo discontinuo está bien localizado y el salto no se transmite a elementos vecinos (a diferencia de XFEM, donde el salto sí cruza fronteras inter-elemento).

### 2. Cinemática de deformaciones

Deformación en el bulk (parte regular del campo):

$$\boldsymbol\varepsilon^{\text{bulk}}(\mathbf x) = \mathbf B_{\text{std}}\,\mathbf d \;-\; \mathbf B^\varphi(\mathbf x)\,\mathbf G\,\llbracket\mathbf u\rrbracket$$

donde `B_std` es la matriz `B` del CST estándar (Voigt 2D, 3×6, constante, ver [`Tri3`](../../../specs/Tri3.md) §8) y `B^φ(x)` es la matriz Voigt construida con `∇φ`:

$$\mathbf B^\varphi = \begin{bmatrix} \partial\varphi/\partial x & 0 \\ 0 & \partial\varphi/\partial y \\ \partial\varphi/\partial y & \partial\varphi/\partial x \end{bmatrix}$$

Como `φ = Σ N_i^+`, las derivadas `∂φ/∂x`, `∂φ/∂y` son **constantes sobre el elemento** (CST lineal). Por tanto `B^φ` es constante.

Salto cinemático en `Γ_d` (en frame local):

$$\llbracket\mathbf u(\mathbf x_{\Gamma_d})\rrbracket = \mathbf G^\top\,\bigl[\,\mathbf u(\mathbf x^+) - \mathbf u(\mathbf x^-)\,\bigr] \;=\; \llbracket\mathbf u\rrbracket \qquad (\text{i.e. los propios DOFs enriquecidos})$$

por construcción de la cinemática KOS — el salto en `Γ_d` coincide exactamente con la variable enriquecida, sin contribución de `d` (continuidad del campo estándar atravesando `Γ_d`).

### 3. Ecuación constitutiva

**Bulk** (en todo el dominio del elemento `Ω_e`):
$$\boldsymbol\sigma(\mathbf x) = \mathbf C_e \,\boldsymbol\varepsilon^{\text{bulk}}(\mathbf x)$$
con `C_e` la matriz constitutiva elástica isótropa (parámetros `E`, `ν`, hipótesis plane stress / plane strain del material `Elastic2D`).

**Discontinuidad** (sobre `Γ_d`):
$$\mathbf t(\mathbf x_{\Gamma_d}) = \mathbf T_{\text{coh}}\bigl(\llbracket\mathbf u\rrbracket, \text{estado cohesivo}\bigr)$$
con `T_coh` el operador del material cohesivo (en frame local; ver las specs [`CohesiveDamageIsotropic`](../../../specs/CohesiveDamageIsotropic.md) y [`CohesiveDamageMixedMode`](../../../specs/CohesiveDamageMixedMode.md)).

### 4. Equilibrio — forma fuerte

Balance estático en cada subdominio (`Ω⁻`, `Ω⁺`):

$$\nabla\cdot\boldsymbol\sigma + \mathbf b = \mathbf 0 \quad \text{en } \Omega \setminus \Gamma_d .$$

**La forma débil KOS (§5) no impone la continuidad de tracción $\boldsymbol\sigma\cdot\mathbf n = \mathbf t$ punto a punto** (corrección del 2026-11-10, verificada simbólicamente por el `formulador` y el `teorico_fem`). Con $\psi$ el ángulo entre $\mathbf n$ y la normal $\mathbf n_e$ al lado opuesto al nodo solitario, la ecuación del salto da

$$t_n = \frac{A}{h\,l_d}\,(\cos\psi\,\sigma_{nn} + \sin\psi\,\sigma_{ns}), \qquad t_s = \frac{A}{h\,l_d}\,(\cos\psi\,\sigma_{ns} + \sin\psi\,\sigma_{ss}),$$

y con $l_d = (A/h)\cos\psi$ (§11): $t_n - \sigma_{nn} = \tan\psi\,\sigma_{ns}$ y $t_s - \sigma_{ns} = \tan\psi\,\sigma_{ss}$. En la activación, con $\mathbf n$ principal, $t_n = \sigma_I$ y $t_s = \sigma_{II}\tan\psi$. Sólo con $\psi = 0$ la continuidad es completa.

### 5. Equilibrio — forma débil (residuales acoplados)

Aplicando el principio de trabajos virtuales con variaciones independientes `δd` y `δ[[u]]`, se obtienen **dos residuales acoplados**:

$$\mathbf R^d = \int_{\Omega_e} \mathbf B_{\text{std}}^\top\,\boldsymbol\sigma\,d\Omega \;-\; \mathbf f^{\text{ext}}_d = \mathbf 0$$

$$\mathbf R^{\llbracket u\rrbracket} = -\int_{\Omega_e} \mathbf G^\top\,(\mathbf B^\varphi)^\top\,\boldsymbol\sigma\,d\Omega \;+\; l_d\,\mathbf t = \mathbf 0$$

con `l_d` la longitud efectiva de `Γ_d` dentro del elemento (ver §11). El primer residual se ensambla al sistema global; el segundo se resuelve **localmente** dentro del elemento (ver §10, condensación estática).

**Cargas exteriores.** Sólo entran en `f^ext_d`, el residual de los DOFs regulares; el salto no recibe carga exterior. Una tracción o presión sobre un lado del elemento (`tractions` y `tractions_by_group` del YAML; `compute_edge_traction`) da el mismo vector nodal que el `Tri3` padre, $\int_{\Gamma} \mathbf N^\top \bar{\mathbf t}\,t\,ds$, esté o no agrietado el elemento y aunque `Γ_d` cruce el lado cargado (2026-11-20; antes el elemento no admitía cargas de borde).

---

## Formulación numérica (FEM)

### 6. Discretización del CST padre

Idéntica a [`Tri3`](../../../specs/Tri3.md): tres nodos en orden antihorario, mapeo isoparamétrico desde el triángulo de referencia con vértices `(0,0)`, `(1,0)`, `(0,1)`. Jacobiano constante; aborta con error si `det J ≤ tol`.

**Identificación del nodo solitario.** Cuando `Γ_d` corta dos lados del triángulo, exactamente un nodo queda aislado en `Ω⁺` (lado al que apunta `n`); los otros dos quedan en `Ω⁻`. El nodo solitario es el **único** que tiene `φ = 1`; los otros dos tienen `φ = 0`. **El lado de cada nodo lo da la trayectoria** (signo de los valores nodales de la grieta en [`CrackTracking`](CrackTracking.md)), igual para todos los elementos que comparten el nodo. La regla anterior (signo de `(x_i − x_c)·n` respecto del centroide de cada elemento) sólo sobrevive para la raíz de una grieta nueva.

### 7. Funciones de forma + función `φ`

`N_i(x)` idénticas a [`Tri3`](../../../specs/Tri3.md) §7. Para `φ`:

$$\varphi(\mathbf x) = N_{i^*}(\mathbf x)$$

donde `i*` es el índice del nodo solitario. Es lineal y `C^0` continua, vale 1 en `i*` y 0 en los otros dos. Sus derivadas son `(∂N_{i*}/∂x, ∂N_{i*}/∂y)`, **constantes** sobre el elemento.

### 8. Matrices B estándar + B^φ

`B_std` es la matriz B clásica del CST (constante, 3×6); ver [`Tri3`](../../../specs/Tri3.md) §8.

`B^φ` es la matriz Voigt 3×2 construida con `(∂φ/∂x, ∂φ/∂y)`:

$$\mathbf B^\varphi = \begin{bmatrix} \partial\varphi/\partial x & 0 \\ 0 & \partial\varphi/\partial y \\ \partial\varphi/\partial y & \partial\varphi/\partial x \end{bmatrix}$$

Por construcción, `B^φ` coincide con la **submatriz B_std correspondiente al nodo solitario** (las columnas `2·i*` y `2·i*+1` de `B_std`). Eso permite implementación trivial sin recomputar derivadas.

### 9. Sistema acoplado — matrices elementales

Linealización del residual acoplado conduce al sistema 8×8 local (6 DOFs estándar + 2 enriquecidos):

$$\begin{bmatrix} \mathbf K_{dd} & \mathbf K_{d\llbracket u\rrbracket} \\ \mathbf K_{\llbracket u\rrbracket d} & \mathbf K_{\llbracket u\rrbracket\llbracket u\rrbracket} \end{bmatrix} \begin{bmatrix} \Delta\mathbf d \\ \Delta\llbracket\mathbf u\rrbracket \end{bmatrix} = \begin{bmatrix} -\mathbf R^d \\ -\mathbf R^{\llbracket u\rrbracket} \end{bmatrix}$$

con (todas las cantidades evaluadas en el único punto de Gauss central, peso `A_e·t`):

- `K_dd = B_std^T · C_e · B_std · A_e · t` (6×6) — rigidez estándar del CST.
- `K_{d[[u]]} = -B_std^T · C_e · B^φ · G · A_e · t` (6×2) — acoplamiento.
- `K_{[[u]] d} = -G^T · (B^φ)^T · C_e · B_std · A_e · t` (2×6) — traspuesta de `K_{d[[u]]}` en KOS (simétrica por construcción).
- `K_{[[u]][[u]]} = G^T · (B^φ)^T · C_e · B^φ · G · A_e · t + l_d · T_coh` (2×2) — rigidez del salto, suma de contribución del bulk (descarga elástica del modo enriquecido) y de la rigidez tangente cohesiva `T_coh ≡ ∂t/∂[[u]]` proyectada sobre `Γ_d`.

Convención de signos: el `-` en `K_{d[[u]]}` y `K_{[[u]]d}` viene del signo del término enriquecido en `B^enr = -B^φ·G` (recordar `(H - φ)` con `∇H = δ_{Γ_d}·n`).

### 10. Condensación estática local

El sistema 8×8 se condensa eliminando los DOFs enriquecidos `Δ[[u]]` antes de devolver `K_e` al ensamblador (ADR 0010 §3):

$$\Delta\llbracket\mathbf u\rrbracket = -\mathbf K_{\llbracket u\rrbracket\llbracket u\rrbracket}^{-1}\,\bigl(\mathbf R^{\llbracket u\rrbracket} + \mathbf K_{\llbracket u\rrbracket d}\,\Delta\mathbf d\bigr)$$

Sustituyendo en la primera ecuación:

$$\boxed{\;\mathbf K_{\text{cond}} = \mathbf K_{dd} - \mathbf K_{d\llbracket u\rrbracket}\,\mathbf K_{\llbracket u\rrbracket\llbracket u\rrbracket}^{-1}\,\mathbf K_{\llbracket u\rrbracket d}\;}$$

$$\boxed{\;\mathbf F_{\text{int}}^{\text{cond}} = \mathbf F_{\text{int}}^d - \mathbf K_{d\llbracket u\rrbracket}\,\mathbf K_{\llbracket u\rrbracket\llbracket u\rrbracket}^{-1}\,\mathbf R^{\llbracket u\rrbracket}\;}$$

donde `F_int^d = B_std^T · σ · A_e · t` es la fuerza interna estándar antes de condensar.

**Inversión de `K_{[[u]][[u]]}` (2×2)**: cerrada analíticamente. El despachador algebraico (ADR 0003) no participa al ser una matriz elemental densa de 2×2.

**Sin recuperación posterior** (corregido 2026-09-29): el salto se resuelve con un Newton local dentro de `compute_element_state()` para cada `d` (Retama 2010, §7.2), y la condensación conserva el residuo local. `commit_state()` sólo copia `trial → committed`.

### 11. Activación: criterio Rankine + dirección `n` + `l_d` del Cap. 6

**Criterio de activación**. Con cada campo convergido (al final de cada paso, antes de confirmarlo, y al inicio del siguiente; nunca dentro de las iteraciones del Newton — ver §14), el elemento intacto evalúa su esfuerzo principal mayor `σ_I` en el centroide:

$$\sigma_I = \frac{\sigma_{xx} + \sigma_{yy}}{2} + \sqrt{\Bigl(\frac{\sigma_{xx} - \sigma_{yy}}{2}\Bigr)^2 + \sigma_{xy}^2}$$

Si `σ_I > σ_t0` (resistencia a tracción del material cohesivo), se activa:

- `n` = dirección principal de `σ_I` (vector unitario; criterio de Rankine): la grieta es perpendicular a ella.
- `s` = ortogonal a `n` orientada por regla de la mano derecha con `+z` saliendo del plano (consistente con Reglas.md §5).
- El nodo solitario `i*` lo da la trayectoria (§6). La posición de `Γ_d` dentro del elemento no entra en las ecuaciones: sólo `i*` y `l_d` (informes del 2026-11-10). La «superficie de equilibrio» de la ec. (8) de la deducción del usuario (la recta de normal `n` cuya cuerda mide `l_d`) queda como indicador de posproceso; existe en todo triángulo con ángulo mínimo ≥ arcsin(1/3) ≈ 19.47°.
- **Regla del vértice** (decisión del usuario, 2026-11-10): si la partición por la trayectoria da `cos ψ ≤ 0`, un nodo cambia de lado y la trayectoria pasa por ese vértice; en la punta de la grieta el nodo aislado ya está fijado y cambia el tercer nodo ([`CrackTracking`](CrackTracking.md) §5). Así todo elemento activado tiene `cos ψ > 0` y `l_d > 0`.

**Longitud efectiva `l_d` — Cap. 6 de Retama (2010)**. Fórmula cerrada:

$$\boxed{\;l_d = \frac{A_e}{h}\,\cos(\theta - \alpha)\;}$$

donde:
- `A_e` es el área del triángulo.
- `h` es la altura desde el nodo solitario `i*` al lado opuesto.
- `θ` es el ángulo de `Γ_d` respecto a `x` global (i.e. `θ = atan2(s_y, s_x)` con `s = (s_x, s_y)` la tangente a `Γ_d`).
- `α` es el ángulo del lado opuesto al nodo solitario respecto a `x` global.

Esta fórmula es **decisión cerrada en el ADR 0010** (§6); no se relitiga aquí. **Se toma con su signo** (ya no `|cos(θ−α)|`): la regla del vértice garantiza `cos ψ > 0`. Es válida para todo θ: la identidad (6.21) de la tesis sólo vale con `cos 2θ > 0`, pero su error de signo se compensa con el de la expresión de `σ_I` con raíz, y (6.24) es correcta en general (verificado simbólicamente, 2026-11-10). La alternativa ingenua `l_d = A_e/h` (independiente de `θ`) produce stress locking cuando la grieta no es paralela al lado opuesto; ver fase 3 del ADR 0010 para test numérico.

### 12. Estado de la discontinuidad — `DiscontinuityState`

Dataclass paralela a `ElementState` (ADR 0010 §4), específica de la discontinuidad. Vive en `solidum/user/discontinuities/discontinuity_state.py`:

```python
@dataclass
class DiscontinuityState:
    normal: np.ndarray        # n ∈ ℝ², unitario, fijado al activar
    tangent: np.ndarray       # s ∈ ℝ², s = (−n_y, n_x)
    centroid: np.ndarray      # punto por el que pasa Γ_d
    solitary_node: int        # índice (0, 1 o 2) del nodo en Ω⁺
    l_d: float                # (A_e/h)·cos(θ − α)
    jump_committed: np.ndarray   # [[u]] al cierre del último paso convergido
    jump_trial: np.ndarray       # [[u]] dentro de la iteración Newton
    cohesive_state_committed: dict  # estado interno del material de interfaz (committed)
    cohesive_state_trial: dict      # estado interno del material de interfaz (trial)
```

Semántica **trial / commit** consistente con `ElementState`: `jump_trial` y `cohesive_state_trial` evolucionan en cada iteración del Newton global; al converger el paso, `commit_state()` copia `trial → committed` (el salto ya lo resolvió el Newton local, §10).

### 13. Cuadratura

Un único punto de Gauss central, peso `A_e·t` — idéntico a `Tri3`. Las matrices `B_std`, `B^φ` y la deformación del bulk son constantes; la integración es exacta. El término cohesivo `l_d·t` se evalúa directamente sobre `Γ_d` (no requiere cuadratura porque `t` es constante a lo largo de `Γ_d` por la cinemática KOS — el salto es constante).

### 14. Activación dentro del paso — interacción con el solver

La activación se evalúa **sólo con campos convergidos**, nunca dentro de las iteraciones del Newton: ahí `σ_I` oscila por el predictor lineal y puede cruzar y descruzar el umbral varias veces antes de converger (*chattering*, ADR 0010 §5).

**Al final de cada paso, antes de confirmarlo** (decisión del usuario del 2026-11-18; ADR 0041): cuando el Newton del paso converge, el seguimiento de grieta evalúa el criterio con ese campo. Si activa algún elemento, el paso **no se confirma**: se repite desde el mismo estado confirmado del paso anterior, con las grietas nuevas (salto inicial nulo). Se repite hasta que el seguimiento no activa ningún elemento nuevo (**cascada**: al abrirse una grieta se carga el vecino, que puede entrar en el mismo paso). Así, al confirmar un paso, ningún elemento intacto que las reglas del seguimiento permiten activar (punta de una grieta o nucleación admisible) queda por encima de su resistencia. Los que el seguimiento excluye (contigüidad, `max_cracks` alcanzado, punta con cos ψ ≤ 0 sin nodo libre: `CrackTracking` §4–§5) pueden quedar intactos por encima de `σ_t0`. La evaluación al converger con repetición del paso es la del programa de la tesis (FEAP), con la que Solidum reproduce su salida (Diálogo, 2026-11-11); la diferencia es que allí entra un elemento por punta en cada evaluación y aquí todos los que el seguimiento admite.

**Arranque del salto al repetir.** El Newton local del salto arranca de su último valor de prueba, también cuando el paso se repite (y tras una bisección, como antes de este cambio): la ronda repetida hereda el salto de la descartada. En pasos regulares el efecto es de 1e-15; en pasos críticos ese arranque decide si el Newton converge (`informes/verificacion_teorica_2026-11-18_activacion_precommit.md`, hallazgo 1). Devolverlo al confirmado al repetir cambia qué equilibrio se alcanza en esos pasos y se medirá aparte.

**Intento que fracasa.** Si el intento repetido no converge, o supera el máximo de revisiones (ADR 0041), el solver lo trata como cualquier intento fallido (bisección donde la hay; sin ella, el análisis se detiene). Las grietas activadas durante ese intento **se conservan** (irreversibilidad), aunque su normal y su partición vengan de un campo que no se confirma; el paso reducido parte de ellas. Es lo que hacía la variante D del prototipo. Con la ley penalizada el efecto es pequeño: el elemento agrietado se comporta elásticamente hasta que `t_n` alcanza `σ_t0`. Caso medido en que importa (`teorico-solidum`, 2026-11-18): franja con `max_cracks = 1` y 400 pasos; la revisión activa un elemento, la ronda repetida no converge (ciclo de periodo 3) y las bisecciones no lo arreglan porque la activación persiste (la regla anterior se detenía en el mismo punto).

**Al inicio de cada paso** sigue la evaluación con el campo confirmado del paso anterior (gancho P5 del ADR 0020, `prepare_step(U_committed)`). Con la revisión al final de cada paso normalmente no encuentra nada nuevo; actúa en el primer paso de un análisis o de una etapa, cuyo campo de partida no pasó por una revisión.

**Implementación**: dos ganchos genéricos del programa principal. `prepare_step(U_committed)` (P5) lo invocan los solvers una vez por intento de paso, antes del primer Newton; `revise_step(U_trial) -> bool` (P9, ADR 0041) lo invoca el corrector compartido al converger, antes de confirmar, y si devuelve `True` repite el intento. La base hace los dos no-op; `CST_Embedded2D` los sobreescribe y delega en su seguimiento de grieta. Una vez activado, el elemento queda agrietado **irreversiblemente** (no revierte si el siguiente paso lo descarga), lo que garantiza que la cascada termina.

---

## Caveats numéricos

- **Invertibilidad de `K_{[[u]][[u]]}` en ablandamiento**. La tangente cohesiva es `diag(dT_soft/dκ, 0)`: negativa en el ablandamiento, nula con la grieta totalmente abierta y nula siempre en la dirección tangencial (Modo-I). La matriz 2×2 sigue siendo invertible por el término del bulk `G^T·(B^φ)^T·C_e·B^φ·G·vol`, definido positivo, mientras el elemento no sea tan grande que el ablandamiento lo supere (aproximadamente `l_d·|dT_soft/dκ| < E/h`). Por eso el cohesivo no necesita ninguna rigidez residual; la que tuvo hasta el 2026-09-23 (`(1 − DAMAGE_MAX)·K_e`) daba a la tangente condensada errores del 26 al 45 % en el ablandamiento ([`TestCondensedTangentInSoftening`](../../../../tests/user/discontinuities/test_disc_cst_embedded.py)).
- **Mecanismo de deslizamiento**. Con rigidez de modo II nula (Retama 2010, p. 67), una grieta que atraviesa todo el sólido deja deslizar sin resistencia una parte respecto de la otra en la dirección de la grieta. Las condiciones de apoyo tienen que impedir ese movimiento; si no, el Newton global deriva sin converger en desplazamientos aunque el residuo sea nulo.
- **Solvers que atraviesan el ablandamiento con penalización cohesiva rígida** (nota del 2026-09-29, recogida de la limitación que se retiró ese día de la spec de `NonlinearSolver`). `NonlinearSolver` con control por desplazamiento atraviesa el ablandamiento hasta la separación completa con $K_e = 10^{15}$, frente a la solución exacta ([`test_disc_embedded_uniaxial_softening.py`](../../../../tests/user/discontinuities/test_disc_embedded_uniaxial_softening.py)); la limitación antigua («no atraviesa la transición elástico→ablandamiento con $\mathbf K_t$ casi singular en $\kappa_0$») venía de la tangente cohesiva de signo erróneo, corregida el 2026-09-23 (commit `597f3d9`). El fallo registrado de `ArcLengthSolver` con esa misma penalización es anterior a esa corrección y **no se ha vuelto a comprobar**: puede que ya no exista, pero ningún test lo confirma.
- **Activación dentro del paso** (§14). Ningún elemento que el seguimiento puede activar queda confirmado intacto por encima de su resistencia: desaparece el sobrepaso de un paso que tenía la activación al inicio del paso siguiente (Van Vliet: +0.9 kN en el pico con pasos de 1 µm, Diálogo del 2026-09-29). El sesgo cambia de signo: el pico confirmado queda acotado por arriba y muestreado, por debajo del real hasta ≈ (Δδ_exc/δ_pico)·m/(1 − m), con Δδ_exc lo que el paso rebasa el alargamiento del pico y m = |f′|·L/E (probeta de la aceptación: −1.2 % lineal, −2.8 % exponencial; *calculado* por `teorico-solidum`); converge al pico al refinar el paso. También depende del tamaño del paso cuántos elementos entran a la vez: la primera revisión usa el campo sin la grieta nueva, y todos los elementos que lo superan con ese campo se activan en esa revisión, con su normal; con pasos grandes cerca del pico pueden entrar juntos elementos que con pasos finos entrarían uno tras otro (Diálogo de `CrackTracking`, 2026-11-10). La regla nueva no corrige eso. Cada revisión repite el Newton del paso: un paso con *n* revisiones cuesta *n* + 1 Newton.
- **Bulk descarga elástica**. Por la cinemática KOS, el bulk evalúa `ε^bulk = B_std·d − B^φ·G·[[u]]` — descarga conforme `[[u]]` crece. Esto es físicamente correcto (la disipación va a `Γ_d`) y consistente con la discrete approach de Retama 2010. El bulk **no acumula daño** propio (out_of_scope fase 1).
- **`Γ_d` paralelo a un lado**. Caso degenerado: `cos(θ−α) = 1`, `l_d = A_e/h`. La fórmula del Cap. 6 sigue siendo correcta, pero numéricamente el nodo solitario debe identificarse con tolerancia para evitar oscilación entre dos nodos casi en el mismo plano.
- **`Γ_d` casi vertical en CST con jacobiano malacondicionado**. Sin patología nueva: hereda los caveats del CST padre. El test `det J > tol` de `Tri3` cubre el caso.
- **Errores de la tracción KOS y bloqueo de esfuerzos con estados biaxiales oblicuos** (2026-11-10). Por §4, la tracción difiere de `σ·n` en `tanψ·σ_ns` (normal) y `tanψ·σ_ss` (tangencial). Con la ley de modo I puro (`t_s = 0`) la segunda ecuación del salto impone `σ_ns = −tanψ·σ_ss` en el bulk, aparece un deslizamiento espurio (`ũ_s/ũ_n → sinψ cosψ/(1 + sin²ψ)`) y, en un elemento aislado con `ε_yy = 0.5·ε_xx`, al abrir del todo quedan esfuerzos de bulk `σ_xx/σ_t0` = 0.89, 3.9 y 9.4 para `ψ` = 15°, 30° y 45°. **Decisión 4 del usuario (2026-11-10): se mantiene la KOS con t_s = 0 y normal de Rankine.** Ese bloqueo es del elemento aislado con deformación impuesta en sus tres nodos, que no puede redistribuir; en una malla con particiones coherentes no aparece (datos en el Diálogo, 2026-11-10).
- **Stress locking**. La fase 3 del ADR 0010 valida numéricamente que `l_d = (A_e/h)·cos(θ−α)` evita el locking que aparece con `l_d = A_e/h` ingenuo. Test cubierto en `tests/user/discontinuities/test_disc_ld_chapter_6_validation.py`: con `Γ_d` oblicua respecto al lado opuesto al nodo solitario (`cos(θ−α) = √3/2` en el setup), la versión correcta produce una apertura `[[u_n]]` ≈ 5.45 % mayor que la versión ingenua para una misma carga, en rama de softening exponencial activo. La versión ingenua "infla" la rigidez cohesiva (`l_d · T_coh`) por un factor `1/cos(θ−α) > 1`, restringe la apertura y el elemento se comporta como si fuera más rígido — eso es el *stress locking*. Cuando `Γ_d` es paralela al lado opuesto (`cos(θ−α) = ±1`), las dos fórmulas coinciden exactamente y no hay locking.

---

## Contrato de implementación

```yaml
name: CST_Embedded2D
kind: element
status: validated

interface:
  dof_names: [ux, uy]           # 2 DOFs globales por nodo (los enriquecidos NO son globales)
  n_nodes: 3
  strain_dim: 3                 # Voigt 2D [ε_xx, ε_yy, γ_xy]
  n_integration_points: 1       # CST con cuadratura central

parameters:
  - { name: thickness,         type: float, required: false, default: 1.0,
      desc: "Espesor para estado plano (idéntico a Tri3)." }
  - { name: cohesive_material, type: str,   required: true,
      desc: "Id (en interface_materials, InterfaceMaterialRegistry) del material de interfaz de superficie que gobierna Γ_d cuando se activa." }
  - { name: crack_tracker,     type: str,   required: false,
      desc: "Seguimiento de grieta (CrackTrackerRegistry: CrackTracking o PrescribedCrackPath) compartido por la región. Sin él, un CrackTracking por región conexa de elementos embebidos." }

material_contract:
  bulk:
    signature: "compute_state(ε, state) -> (σ, C_tangent, state')"
    strain_kind: "Voigt 2D [ε_xx, ε_yy, γ_xy]"
    accepted: ["Elastic2D"]      # fase 1: bulk elástico; otros → out_of_scope
  cohesive:
    signature: "compute_traction([[u]], state) -> (t, T_tan, state')"
    jump_kind: "frame local — [[u]] = ([[u_n]], [[u_s]]) en ejes (n, s) de Γ_d"
    accepted: ["CohesiveDamageIsotropic", "CohesiveDamageMixedMode"]  # cualquier InterfaceMaterial de superficie que admita 2 componentes

conventions:
  sign: "tracción interna positiva; convención stress-resultant del proyecto (Reglas.md §5)."
  voigt: "[xx, yy, xy] con γ_xy = 2·ε_xy."
  node_orientation: "antihorario (det J > 0)."
  frame_local_discontinuity: "n perpendicular a σ_I al activarse; s = (−n_y, n_x) (RHR con +z saliendo del plano)."
  activation: "irreversible; evaluada con el campo convergido de cada paso antes de confirmarlo (el paso se repite con las grietas nuevas hasta que no entra ninguna; ADR 0041) y al inicio del paso con el campo confirmado (ADR 0020, P5); nunca dentro de las iteraciones del Newton."

validity:
  - "Pequeños desplazamientos y deformaciones (no large strain)."
  - "Bulk elástico (Elastic2D) y material de interfaz de superficie con dos componentes (modo I, CohesiveDamageIsotropic, o modo mixto, CohesiveDamageMixedMode)."
  - "Γ_d cortando exactamente dos lados del triángulo (configuración geométrica estándar)."

out_of_scope:
  - "Bulk no lineal (J2, daño continuo, Drucker-Prager): la disipación debe localizarse en Γ_d, no co-existir con disipación de bulk."
  - "Cinemática que transmita cortante sin bloqueo de esfuerzos con la ley de modo mixto (SKON, no simétrica): con CohesiveDamageMixedMode el elemento sigue siendo KOS, que puede bloquear esfuerzos con la grieta muy abierta y esfuerzo paralelo grande; la elección se decidía con un ensayo de referencia de modo mixto (decisión del usuario del 2026-11-10); tras las vigas de Gálvez el usuario decidió el 2026-11-11 seguir sólo con la KOS simétrica de la tesis (Diálogo)."
  - "Reorientación de n tras la activación."
  - "Unión y ramificación de grietas (una discontinuidad por elemento; ver CrackTracking)."
  - "Múltiples discontinuidades por elemento."
  - "Contacto con fricción en la grieta cerrada (el cohesivo sólo recupera la penalización normal K_e en compresión)."
  - "Elementos de orden superior con discontinuidad embebida (Tri6_Embedded, Quad8_Embedded; fase J)."
  - "3D (fase I)."

acceptance:
  verification:
    - name: estado_intacto_reproduce_Tri3
      setup: "elemento sin activar (discontinuity_state=None), bulk Elastic2D; aplicar campo lineal de desplazamientos."
      expect: "K_e y F_int_e idénticos (bit-exact en doble) a los de un Tri3 con el mismo bulk."
      tol_rel: 1.0e-14

    - name: patch_test_intacto
      setup: "patch MacNeal-Harder con 4 CST_Embedded2D distorsionados, ninguno activado, σ < σ_t0 en todo el dominio."
      expect: "campo lineal exacto en nodos interiores; ε constante por elemento."
      tol_rel: 1.0e-10

    - name: activacion_rankine_correcta
      setup: "tracción uniaxial pura en x sobre un único CST con σ_t0 conocido; incrementar carga hasta cruzar σ_t0."
      expect: "activación en el paso donde σ_xx > σ_t0; n = (1, 0) y s = (0, 1) con error angular < 1e-10."
      tol_rel: 1.0e-10

    - name: bulk_descarga_elasticamente_post_activacion
      setup: "activar elemento, después fijar [[u]] e imponer Δd; verificar que la deformación del bulk responde como elástica con C_e."
      expect: "σ_bulk = C_e · (B_std·d − B^φ·G·[[u]]); coincide con C_e aplicado a la ε neta."
      tol_rel: 1.0e-12

    - name: condensacion_simetrica_KOS
      setup: "elemento activado con cohesivo en carga activa simétrica (Modo-I puro)."
      expect: "K_cond simétrica: ‖K_cond − K_cond^T‖_F / ‖K_cond‖_F < 1e-12."
      tol_rel: 1.0e-12

    - name: ld_formula_cerrada
      setup: "varias orientaciones de Γ_d (θ ∈ {0, π/6, π/4, π/3, π/2}) sobre un triángulo conocido."
      expect: "l_d coincide con (A_e/h)·cos(θ−α) calculado analíticamente."
      tol_rel: 1.0e-14

    - name: recuperacion_local_del_salto
      setup: "carga monótona hasta varios pasos post-activación; tras cada commit, [[u]]_committed debe satisfacer R^{[[u]]} = 0."
      expect: "‖R^{[[u]]}‖ < 1e-8 tras commit_state."
      tol_rel: 1.0e-8

    - name: activacion_dentro_del_paso_sin_sobrepaso
      setup: "probeta de dos CST en tracción uniaxial hasta la separación (test_disc_embedded_uniaxial_softening.py), 75 pasos de desplazamiento impuesto, softening lineal y exponencial; el paso que cruza el pico (paso 10) lleva el alargamiento 6.7 % (1/15) más allá del de pico: con el elemento intacto, σ = 1.067·σ_t0."
      expect: "ningún estado confirmado supera σ_t0·W·t (F ≤ F_pico·(1 + 1e-12)); el paso 10 se confirma agrietado tras una revisión, con F = 0.988288·F_pico (lineal) y 0.972121·F_pico (exponencial), sobre la rama de ablandamiento exacta de la ley penalizada (κ = w + κ₀), tol_rel 1e-6. Con la activación al inicio del paso siguiente ese estado se confirmaba intacto con F = 1.0667·F_pico (la prueba falla con ese comportamiento). Valores calculados por teorico-solidum con un modelo 1D propio (informes/verificacion_teorica_2026-11-18_activacion_dentro_del_paso.md)."
      tol_rel: 1.0e-6

    - name: traccion_de_borde_como_Tri3
      setup: "triángulo distorsionado con espesor 0.37; tracción constante, presión y campo t(x, n) en los tres lados; el mismo elemento antes y después de activarlo; viga SEN sin placa inferior con tractions_by_group (pressure) en la cara x = 0.23–0.25."
      expect: "el vector nodal es el de Tri3 bit a bit e igual antes y después de la activación; con presión uniforme, p·L·t/2 en cada nodo del lado según −n y cero en el tercero; en la viga SEN, p·t·L/2 por extremo de cada lado (las fuerzas que el prototipo del 2026-11-19 sumaba a mano) y resultante 10/11 de la carga hacia arriba."
      tol_rel: 1.0e-12

  specific:
    - name: bulk_no_lineal_rechazado
      setup: "instanciar CST_Embedded2D con bulk_material=IsotropicDamage2D."
      expect: "ValueError clara en construcción, mensaje listando los bulks aceptados (fase 1: Elastic2D)."

    - name: jump_dim_validado
      setup: "instanciar con cohesive_material cuyo JUMP_DIM=3 (futuro 3D)."
      expect: "ValueError con mensaje sobre incompatibilidad dimensional (elemento 2D requiere JUMP_DIM=2)."

    - name: activacion_no_chattering_dentro_de_newton
      setup: "elemento cerca del umbral; correr iteraciones Newton sin invocar prepare_step."
      expect: "estado discontinuity_state no cambia durante las iteraciones; sólo prepare_step modifica activación."

    - name: identificacion_nodo_solitario
      setup: "construir CST con coordenadas tales que solo un nodo cae en Ω⁺ para una n dada."
      expect: "DiscontinuityState.solitary_node identifica correctamente el índice; los otros dos tienen φ=0."

    - name: ld_invariante_bajo_reflejo_de_n
      setup: "activar con n y comparar con activar con −n (mismo problema físico)."
      expect: "l_d, K_cond, F_int^cond idénticos (la convención de n positivo está bien definida — apunta hacia el nodo solitario)."
      tol_rel: 1.0e-14

    - name: tres_triangulos_traccion_hasta_separacion
      setup: "rectángulo (ψ = 0°), agudo (ψ = 24°) y obtuso (ángulo mínimo 17°, ψ = 73°), ε_xx impuesta creciente (ε_yy = γ_xy = 0), partición de la raíz con la regla del vértice"
      expect: "t_n ≥ −1e-6·σ_t0 en todo el proceso (la regla anterior con |cos| daba −400·σ_t0 en el obtuso); l_d > 0; 0 fallos del Newton local; rectángulo y agudo se separan con disipación = G_F·l_d·t (error ≤ 1e-3 y ≤ 3e-3); el obtuso se bloquea (t_n final 0.88·σ_t0): bloqueo de esfuerzos KOS con ψ grande en un elemento aislado (en una malla no aparece: decisión 4)"

  arch:
    - name: registro_en_ElementRegistry
      expect: "CST_Embedded2D registrado con kind='element'."

    - name: STRAIN_DIM_y_DOF_NAMES_correctos
      expect: "STRAIN_DIM=3, DOF_NAMES=['ux','uy'], N_NODES=3, N_INTEGRATION_POINTS=1."

    - name: prepare_step_es_no_op_en_Element_base
      expect: "Element.prepare_step existe y no hace nada por defecto; CST_Embedded2D lo sobreescribe."

    - name: revise_step_es_no_op_en_Element_base
      expect: "Element.revise_step existe y devuelve False por defecto; CST_Embedded2D lo sobreescribe y devuelve True sólo en la revisión que activa algún elemento, medido por el número de elementos agrietados antes y después (un cambio de lado de un nodo sin activación no cuenta); idempotente con el mismo campo."

references:
  - "Retama Velasco, J. (2010). Formulation and Approximation to Problems in Solids by Embedded Discontinuity Models. Tesis Doctoral, UNAM. **Caps. 2 (cinemática KOS), 5 (formulación variacional discreta), 6 (l_d correcto), 7 (condensación estática)**."
  - "Jirásek, M. (2000). Comparative study on finite elements with embedded discontinuities. CMAME 188, 307-330. — Clasificación SOS/KOS/SKON."
  - "Oliver, J., Huespe, A.E., Sánchez, P.J. (2006). A comparative study on finite elements for capturing strong discontinuities. CMAME 195, 4732-4752."
  - "Linder, C., Armero, F. (2007). Finite elements with embedded strong discontinuities. IJNME 72, 1391-1433."
  - "ADR 0010 — Discontinuidades interiores embebidas (hoja de ruta; fase 2 = este elemento)."
  - "Specs [CohesiveDamageIsotropic](../../../specs/CohesiveDamageIsotropic.md) y [CohesiveDamageMixedMode](../../../specs/CohesiveDamageMixedMode.md) — materiales de interfaz consumidos en Γ_d."
  - "Spec [Tri3](../../../specs/Tri3.md) — CST padre del elemento; mismo bulk en estado intacto."
```

---

## Implementación

- Archivo: `solidum/user/discontinuities/embedded_cst.py`
- Clase: `CST_Embedded2D` (registrada en `ElementRegistry`)
- Estado de la discontinuidad: `solidum/user/discontinuities/discontinuity_state.py` (`DiscontinuityState`)
- Ganchos de paso en clase base: `Element.prepare_step(U_committed)` (P5, no-op por defecto) y `Element.revise_step(U_trial) -> bool` (P9, ADR 0041, `False` por defecto)
- Integración en solvers:
  - `NonlinearSolver`: invoca `self.assembler.prepare_all_steps(U_current)` al inicio de cada paso, antes del Newton.
  - `ArcLengthSolver` y sus variantes: idem, antes del predictor.
  - Revisión del paso convergido: la hace el corrector compartido (`NewtonCorrector.run`) en todos los solvers que lo usan, con `assembler.revise_all_steps(U)` antes de confirmar; si algún elemento cambió, repite el intento desde su iterado inicial.
- Parser YAML: el material se declara en la sección `interface_materials:` del programa principal (hasta el 2026-11-10, `cohesive_materials:` del módulo) + campo `cohesive_material:` por elemento.
- Tests: `tests/user/discontinuities/test_disc_cst_embedded.py` — 24 tests cubriendo `acceptance` completo (verification, specific, arch) + integración end-to-end del YAML, y la tracción de borde (`TestEdgeTraction`, `test_presion_por_grupo_en_la_viga_sen`, 2026-11-20). El aislamiento del trial del salto frente a `compute_global_stiffness` (deuda #18) se prueba en `tests/user/discontinuities/test_disc_contratos.py::TestTrialIsolation`.
- Notas de traducción:
  - Estado intacto (`discontinuity_state is None`) → bit-exact con `Tri3` (mismo `B`, mismo `compute_integrands`). La cinemática del CST se importa de la API pública del programa principal, `solidum.elements.solid_2d` (`compute_kinematics_tri3`, `compute_integrands`; ADR 0020, P6), no del módulo interno `_shared`.
  - Estado agrietado: Newton local sobre `[[u]]` hasta `R^{[[u]]} = 0`, después condensación estática local cerrada (`K_jj` es 2×2). Tolerancia interna: `_LOCAL_JUMP_RTOL = 1e-10`, `_LOCAL_JUMP_MAX_ITER = 30`.
  - `B^φ` se obtiene como la submatriz `B_std[:, 2·i*:2·i*+2]` correspondiente al nodo solitario (consistente con `φ = N_{i*}`).
  - Activación: la decide el seguimiento de grieta (`crack_tracker`, spec [`CrackTracking`](CrackTracking.md)); el elemento expone `activate(normal, solitary_node, point, crack)`, que rechaza `l_d ≤ 0`, y `principal_stress(U)`. `n` apunta hacia el nodo solitario. `_activate(n, coords)` aplica la partición de una raíz (recta por el centroide, o la del otro lado del nodo intermedio si aquella da `cos ψ ≤ 0`).
  - `_compute_ld` usaba `|cos(θ − α)|` hasta el 2026-11-10; desde entonces `l_d` lleva su signo, positivo por la regla del vértice. La invariancia bajo `n → −n` se mantiene porque `n` apunta siempre hacia el nodo solitario.
  - `compute_gauss_state` extendido: en estado agrietado añade clave `'discontinuity': {normal, tangent, centroid, solitary_node, l_d, jump, traction, damage}` para post-proceso VTK.
  - Bulks aceptados en fase 1: `_ACCEPTED_BULKS = ('Elastic2D',)`. Restricción declarativa, ampliación futura es aditiva.

---

## Diálogo

- **2026-05-18** · Spec pre-redactada por la IA como **orden de trabajo** sobre los Caps. 2, 5, 6 y 7 de Retama 2010. La cinemática KOS, la condensación estática y el `l_d = (A/h)·cos(θ−α)` son traducción directa de la tesis; las decisiones arquitecturales del ADR 0010 (familia paralela, condensación local, irreversibilidad, activación al principio del paso) están fijadas. Pendiente de validación del usuario.

- **2026-05-18** · Usuario delega las decisiones a la IA. Cierre con justificación arquitectural:
  - **A — Bulk restringido a `Elastic2D` en fase 1.** Adoptado. Razón: la *discrete approach* del Cap. 3-7 de la tesis presupone bulk elástico (la disipación se localiza en `Γ_d`); mezclar plasticidad continua del bulk con cohesivo de `Γ_d` es una formulación distinta no cubierta por el ADR 0010. Validación cruzada en `__init__` con mensaje listando bulks aceptados. Apertura futura a otros materiales es aditiva.
  - **B — Hook `prepare_step(U_committed)` en `Element` base, no-op por defecto.** Adoptado. Lo prescribe el ADR 0010 §5 explícitamente; aditivo (no rompe llamadores) y testeable directamente. Alternativa de meter activación en `compute_element_state` con guarda `_step_id` descartada por ser frágil y mezclar conceptos (preparación vs cálculo).
  - **C — `INTERNAL_DOF_NAMES` no se añade al contrato declarativo todavía.** Resistir abstracción especulativa (Reglas.md §1: justificar con dos casos reales). Cuando entre el segundo elemento enriquecido (Quad4_Embedded o Tet4_Embedded), se centraliza la convención.
  - **D — Post-proceso vía `compute_gauss_state(U)` extendido.** Adoptado. En estado intacto se comporta como `Tri3`; en estado agrietado añade la clave `'discontinuity': {'normal', 'jump', 'traction', 'damage', 'l_d'}`. Consistente con el patrón actual del `Tri3` y aditivo para `VtkExporter`. Un método separado `crack_state()` duplicaría el patrón sin beneficio.
  - **E — Constructor `__init__(element_id, nodes, material, cohesive_material, thickness=1.0)`.** Mantengo el nombre `material` para el bulk (consistente con `Element.__init__`, que ya lo recibe así); `cohesive_material` es el parámetro nuevo. El YAML mapea `bulk_material: <name>` → `material=...` y `cohesive_material: <name>` → `cohesive_material=...`.
  - **F — Tolerancia `1e-14` (bit-exact)** en `estado_intacto_reproduce_Tri3`. Si en la implementación emerge un desliz numérico documentado (por canceladas algebraicas tipo `G·G^T = I`), se relaja con razón explícita en el Diálogo. Empezar estricta es mejor que pre-laxar.
- **2026-05-18** · Status → en draft hasta que el código + tests cierren; se promoverá a `validated` cuando los `acceptance` pasen. La IA arranca implementación.
- **2026-09-22** · Auditoría global: tolerancia del jacobiano ahora **relativa** (`det J ≤ JACOBIAN_RTOL · Π‖fila_i(J)‖`, cota de Hadamard, adimensional): la absoluta `1e-10` sobre `det J` rechazaba mallas finas en metros sin distorsión.
- **2026-09-22** · Deuda #18 (estado): `compute_global_stiffness` restaura al salir el trial del material y el del salto (`jump_trial`, `cohesive_state_trial`), que el Newton local escribía al evaluar en `u = 0` (`tests/test_trial_isolation.py`).
- **2026-09-29** · Correcciones documentales: §10 y §12 describían una recuperación del salto en `commit_state()` que el código no hace (Newton local de la §7.2); §14 ya no es punto abierto; se retira el parámetro `activation_criterion` del contrato, que el constructor no acepta (la activación es siempre Rankine).
- **2026-09-29** · **Validación: ensayo de Van Vliet con grieta prescrita** (Retama 2010, §8.1), `tests/user/discontinuities/test_disc_van_vliet.py`, figura en `docs/user/discontinuities/validacion/van_vliet_fase1.png`. Probeta en hueso de perro con carga excéntrica; banda de `CST_Embedded2D` a media altura; `NonlinearSolver` con desplazamiento impuesto en pasos de 0.5 µm (convergidos: 0.25 y 0.125 µm cambian el pico < 0.1 kN). Frente a la fig. 8.3 de la tesis (valores leídos): exponencial, pico 26.6 kN (tesis 27.0) a 7.1 µm, P(20/50/100 µm) = 17.8 / 8.3 / 3.3 kN (19.1 / 8.3 / 3.1); lineal, pico 27.4 kN (27.4), 22.0 / 8.5 / 3.0 kN (22.4 / 8.5 / 2.4). Criterio de aceptación validado por el usuario: pico ±0.8 kN, δ en el pico entre 5 y 9 µm, P(20) ±1.5, P(50) ±1.0, P(100) ±0.7 kN. `K_e` = 1e13 o 1e14 da lo mismo. Diferencias abiertas: (a) más allá de ~100 µm la carga queda alta (exponencial 1.5 kN a 200 µm frente a 0.8): las mitades giran sobre el borde de la banda opuesto a la carga, ese elemento queda comprimido (σ_yy ≈ −31 MPa) y nunca se activa, y la cola depende de la malla (1.19 kN con banda más fina); (b) con pasos de 1 µm el pico sube ~0.9 kN por el sobrepaso de la activación con el estado del paso anterior (fase 2 del plan de 2026-09-24). Supuestos no fijados por la tesis: acero E = 210 GPa, rótulas sin desplazamiento lateral (reacción horizontal nula).
- **2026-09-29** · Correcciones documentales tras endurecer la frontera del programa principal (ADR 0020): (1) la entrada del 2026-09-22 sobre la deuda #18 cita `tests/test_trial_isolation.py`; el aislamiento del trial del salto se prueba ahora en `tests/user/discontinuities/test_disc_contratos.py::TestTrialIsolation` (sección Implementación). (2) La cinemática que reutiliza el estado intacto es la API pública P6 del programa principal, `solidum.elements.solid_2d` (`compute_kinematics_tri3`, `compute_integrands`, sin guion bajo). (3) Nuevo caveat con lo que sigue siendo cierto de la limitación retirada de la spec de `NonlinearSolver`: este solver, con control por desplazamiento, atraviesa el ablandamiento hasta la separación; el fallo de `ArcLengthSolver` con penalización rígida es anterior a `597f3d9` y no se ha vuelto a comprobar. Sin cambio de formulación.
- **2026-11-10** · **Revisión de la formulación** con dos informes independientes (`formulador` y `teorico_fem`, excluidos de Git). El código implementaba exactamente la KOS de la tesis; los defectos estaban en la regla del nodo solitario (elemental e incoherente entre vecinos: bloqueo en mallas no alineadas) y en el `|cos|` (grieta que abre en compresión). **Decisiones del usuario**: (1) el nodo solitario lo dicta la trayectoria; (2) seguimiento global ([`GlobalCrackTracking`](GlobalCrackTracking.md)); (3) regla del vértice con `cos ψ ≤ 0`; `l_d` con su signo. §4 corregida: KOS no impone `σ·n = t` punto a punto. Pendiente del usuario: la componente tangencial (modo mixto).
- **2026-11-10** · Implementado el seguimiento (`crack_tracking.py`). `l_d` con su signo y `activate` lo exige positivo. Aceptación `tres_triangulos_traccion_hasta_separacion` medida: el rectángulo (ψ = 0°) disipa exactamente G_F·l_d·t, el agudo (ψ = 24°) 0.18 % de más, y el obtuso (ψ = 73°, la única partición con cos ψ > 0) se bloquea con t_n final 0.88·σ_t0, sin compresión espuria. Es el bloqueo de esfuerzos KOS con ψ grande que señaló el `teorico_fem` y entra en la decisión pendiente del modo mixto. En una malla la regla del vértice evita estas particiones cuando hay vecinos que puedan repartirse el corte (5/5 mallas con astillas se separan).
- **2026-11-10** · **Decisión 4 del usuario: se mantiene la formulación de la tesis** (KOS simétrica, deslizamiento libre t_s = 0, normal de Rankine), sin cambios de código. Datos (prototipos fuera del repositorio): placa 0.1×0.05 traccionada hasta separarse con compresión lateral constante p = 0, 0.5 y 1.0·σ_t0 (malla estructurada y dos Delaunay): KOS con t_s = 0 da F_fin/F_máx < 1e-3 y trabajo/(G_F·W·t) 1.0004–1.0068; SKON (Jirásek 2000, misma cinemática con t = σ·n, matriz no simétrica) da lo mismo (F_fin ≤ 1e-3, 0.998–1.008); el cortante penalizado y la normal de la trayectoria cambian la energía ≤ 0.3 %. Grieta abierta al 60 % de w_c y después deslizada: con t_s = 0, |F_y| máx = 0.0006–0.002·F_pico con KOS y con SKON (sin grieta sería ≈ 5·F_pico), es decir, sin bloqueo; con cortante penalizado (1 − ω)·K_e·[[u_s]], 0.65·F_pico y de 4 a 11 grietas secundarias, la transferencia de cortante de esa ley. El modo mixto (una ley cohesiva con componente tangencial, p. ej. trabazón de agregados en concreto) queda como ampliación de la ley, no como corrección: interesa al usuario para modelar estructuras de concreto (fase G del ADR 0010).
- **2026-11-10** · El elemento admite la ley de modo mixto del programa principal `CohesiveDamageMixedMode` (spec propia, decisiones del usuario del 2026-11-10: apertura efectiva de Camacho–Ortiz, energía independiente del modo; la cinemática sigue siendo KOS, decisión 1 del modo II). Sin cambio de código en el elemento. Con una ley que transmite cortante aparece una **condición de tamaño en cortante**: el sistema local del salto es definido positivo sólo si `β²·|T'|·h < G` (además de `|T'|·h < E'` en la normal); medido en `tests/user/discontinuities/test_disc_mixed_mode.py`, detalle en §12 de la spec de la ley.
- **2026-11-11** · **Decisiones del usuario** tras las vigas de Gálvez (la KOS se bloquea tras el pico en modo mixto; el prototipo SKON diverge): (1) las mejoras se buscan dentro de la KOS simétrica de la tesis (reglas del elemento, seguimiento, Newton local, malla); por ahora no se proponen KSON, SKON ni saltos lineales de otros autores; (2) los elementos de interfaz no sustituyen a la discontinuidad embebida, cuyo objetivo es modelar la grieta sin ellos. Fase 2 en pausa.
- **2026-11-11** · **Verificación contra el FEAP de la tesis** (modelo de la viga de Schlangen, §8.2, `Iviga4p2.dat`, fuera del repositorio; prototipos en `informes/`, fuera de Git). Sobre ese modelo exacto Solidum reproduce la salida de FEAP: P_max 40.989 frente a 40.999 kN, los mismos 19 elementos agrietados en el mismo orden, normales iguales a ±1°, ΔP medio 0.024 kN (máx. 0.23) en 168 pasos. El elemento KOS, `l_d`, la condensación y la ley no tienen error. Para reproducirla hacen falta: ablandamiento exponencial (la salida de FEAP imprime «Linear»), ley inicialmente rígida (κ₀ = 1e-15), una sola grieta por centroides que avanza por sus dos puntas, un elemento por punta y paso, y activación dentro del paso (Solidum activa al inicio del siguiente). Estas reglas no están en Solidum: el usuario decidió no incorporarlas por ahora. Con fuerzas en razón 1:10 y continuidad de la grieta, la respuesta es objetiva a partir de h/2; con la malla de FEAP (~10 mm) las variantes nuevas no avanzan. Formulación **en pausa** por decisión del usuario.
- **2026-11-18** · **Decisión del usuario: se levanta la exclusión de los saltos lineales** fijada el 2026-11-11 («por ahora no se proponen KSON, SKON ni saltos lineales de otros autores»). Motivo: los límites que quedan en el elemento son cinemáticos; el salto constante no deja girar una cara respecto a la otra (cola de Van Vliet: el elemento de la articulación, comprimido, nunca se activa). Siguen vigentes: la exclusión de KSON y SKON y que los elementos de interfaz no sustituyen a la embebida. Primer paso: búsqueda bibliográfica de si ya existe una KOS simétrica con salto lineal condensado por elemento en triángulos lineales con ley cohesiva (Linder–Armero 2007; Alfaiate–Simone–Sluys 2003, con salto continuo y grados de libertad globales; Dias-da-Costa et al. 2009). Fase 2 reabierta en esa dirección.
- **2026-11-18** · **Decisión del usuario: activación dentro del paso** (paso 1 del plan para la viga SEN, aprobado con «si»). Se lleva a Solidum la variante D del prototipo del 2026-11-11: el elemento de la tesis sin cambios, activación dentro del paso con cascada (§14) y una sola grieta con trayectoria continua (`GlobalCrackTracking`, parámetro `max_cracks`). Sustituye la parte de la entrada del 2026-11-11 que decía que estas reglas no se incorporaban. **Ley cohesiva: la de Solidum basta** (*medido*, prototipo `informes/prototipos_2026-11-18/embebida_paso1/ley_solidum.py`): con `CohesiveDamageIsotropic` exponencial y `K_e = 1e14`, en lugar de la ley inicialmente rígida de la tesis, la variante D en la malla h/2 da un pico de 44.70 kN (prototipo 44.69) y P(CMSD = 0.05/0.10/0.15/0.20 mm) = 19.14/4.55/1.49/1.02 kN (19.13/4.55/1.55/0.92), con 48 elementos agrietados (49). No hace falta una ley nueva. La activación dentro del paso es una pieza genérica del programa principal (gancho P9, ADR 0041) que también quita el sobrepaso de un paso de la activación al inicio del paso siguiente (aceptación `activacion_dentro_del_paso_sin_sobrepaso`).
- **2026-11-18** · **Implementado** (ADR 0041, gancho P9; `GlobalCrackTracking.max_cracks`). Revisiones previas: `teorico-solidum` (lista con cambios de redacción, aplicados; `informes/verificacion_teorica_2026-11-18_activacion_dentro_del_paso.md`) y `arquitecto-solidum` (aprobado con cambios, aplicados; `informes/revision_arquitectura_2026-11-18_adr0041.md`). **Solidum reproduce la variante D** sin los arreglos del prototipo (`informes/prototipos_2026-11-18/embebida_paso1/reproduce_D.py`; *medido*): h/2, pico 44.70 kN (prototipo 44.69), diferencia entre curvas 0.33 kN (0.7 % del pico) sobre una rejilla de 400 puntos y 1.58 kN (3.5 %) en un único paso (el mismo suceso, desplazado medio paso; percentil 99: 0.34 %); h/4, 43.45 kN (43.44), 0.80 kN (1.8 %, en la cola). Entre h/2 y h/4, 1.82 kN (4.1 %) en CMSD = 0.037 mm, como en el prototipo. Con la malla de FEAP coincide a 0.1 % hasta CMSD = 0.099 mm, donde el prototipo saltaba sin equilibrio (12 elementos de golpe, residuo de energía 72 %); Solidum no da ese salto: biseca hasta el mínimo y se detiene (esa malla ya estaba descartada el 2026-11-11). Las pruebas existentes de la embebida pasan sin cambiar tolerancias. **Pendientes del usuario**: (1) si las grietas activadas en un intento que después se abandona (bisección) persisten, como hoy y como en el prototipo, o se revocan (hallazgo 3 del teórico y 2 del arquitecto; nunca ocurrió en la viga SEN); (2) qué cuenta de resoluciones usa la adaptación del paso tras una cascada (hoy, la de la última ronda; ADR 0041 §6).
- **2026-11-18** · **Verificación de la implementación** (`teorico-solidum`, commit con cambios de texto, aplicados; `informes/verificacion_teorica_2026-11-18_activacion_precommit.md`): valores de aceptación a 1e-10, emulación independiente a 1e-15, cascada genuina sin activables por encima de σ_t0 y energía = G_F·t·Σl_d. **Penalización en la viga SEN**: con `K_e = 1e14` la cola de h/2 alterna entre dos respuestas según el tamaño del paso (hasta 2.7 % en área); con `K_e = 1e15` no, y Solidum coincide con el prototipo D (pico 44.689 frente a 44.686 kN, máx |ΔP| = 0.004 kN hasta CMSD 0.14 mm). **Validación como usuario** (`validador-solidum`, `informes/reporte_embebida_sen_2026-11-18.pdf`, mallas libres de Gmsh, `K_e = 1e14`): **a medias**. h = 10 mm y 2.5 mm funcionan (picos 44.04 y 43.53 kN; 2.5 mm frente a la reproducción h/4 del prototipo: 2.1 %), pero h = 5 mm se detiene: **la grieta da la vuelta alrededor de un nodo** (cinco elementos seguidos aíslan el mismo nodo, cuando la trayectoria corre casi sobre una línea de la malla regular) y la punta deja de avanzar; los elementos del camino real superan σ_t0 hasta 6.8 veces sin quedar cortados por la trayectoria. *Medido*: con la activación al inicio del paso (regla anterior) ocurre lo mismo; es un defecto previo de `GlobalCrackTracking`, no de este cambio (`informes/prototipos_2026-11-18/embebida_paso1/diag_h05*.py`).
- **2026-11-20** · **Cargas de borde** (para validar el seguimiento con la viga SEN cargada como presión, decisión del usuario de ese día en el Diálogo de [`CrackTracking`](CrackTracking.md)): el elemento admite `tractions` y `tractions_by_group` (`compute_edge_traction`), con el vector nodal del `Tri3` padre; la carga entra sólo en `f^ext_d` (§5), como ya decía la forma débil. Sin cambio de formulación. Prueba nueva `TestEdgeTraction` y `test_presion_por_grupo_en_la_viga_sen` (aceptación `traccion_de_borde_como_Tri3`): las cinco fallan sin el método (comprobado). Hasta ese día la presión de la viga SEN se añadía a mano como fuerzas nodales en los guiones de los prototipos.
