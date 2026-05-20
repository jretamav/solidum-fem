# CST_Embedded2D — triángulo CST 2D con discontinuidad interior embebida (KOS)

> Orden de trabajo. Pre-redactada por la IA como borrador para validación del usuario (autor de la formulación, Retama 2010). La cinemática KOS y la fórmula `l_d = (A/h)·cos(θ−α)` son traducción directa de los Caps. 2, 5, 6 y 7 de la tesis; la mecánica de implementación (Voigt, registry, signatura, contrato `Element`, condensación) sigue las convenciones del proyecto. El usuario revisa, corrige donde proceda, y aprueba antes de que la IA toque código.

---

## Especificación física

### 0. Descripción general

Elemento finito sólido 2D triangular de tres nodos (CST padre) **enriquecido** con una **discontinuidad interior embebida** `Γ_d` que cruza el dominio del elemento cuando se cumple un criterio de activación. Es la primera implementación de la familia de elementos con DOFs enriquecidos elementales introducida por **ADR 0010**, y la **fase 2** de la hoja de ruta de discontinuidades embebidas.

Estructura del elemento:

- **Estado intacto** (`discontinuity_state is None`): se comporta exactamente como un [`Tri3`](Tri3.md) estándar — CST con un único punto de Gauss, B constante, integración exacta del bulk con peso `A_e·t`.
- **Estado agrietado** (`discontinuity_state is not None`): el campo de desplazamientos del bulk se enriquece con un salto `[[u]] ∈ ℝ²` localizado en `Γ_d`. El bulk vecino descarga elásticamente; la disipación se concentra en `Γ_d` y la gobierna un material cohesivo (`CohesiveMaterial`, ADR 0010 §1).

Aproximación adoptada: **discrete approach** (Retama 2010), **KOS** (Kinematically Optimal Symmetric) según la clasificación de Jirásek (2000). Tracking trivial: la normal `n` y la posición de `Γ_d` se fijan al activar el elemento (perpendicular a `σ_I` por criterio Rankine; `Γ_d` pasa por el centroide) y no se reorientan después. Modo-I puro (consistente con `CohesiveDamageIsotropic` fase 1); modo mixto y KSON quedan diferidos (fases G y H del ADR 0010).

Régimen: pequeños desplazamientos y deformaciones, cuasi-estático. Pensado para hormigón, roca, cerámicos cuasi-frágiles. Bulk **elástico** por construcción de la formulación (la disipación está localizada en `Γ_d`); materiales no lineales del bulk quedan `out_of_scope` en fase 1.

### 1. Cinemática de desplazamientos (formulación KOS)

Campo total de desplazamientos en el elemento:

$$\mathbf u(\mathbf x) = \underbrace{\sum_{i=1}^{3} N_i(\mathbf x)\,\mathbf d_i}_{\mathbf u^{\text{std}}(\mathbf x)} \;+\; \underbrace{\bigl(H_{\Gamma_d}(\mathbf x) - \varphi(\mathbf x)\bigr)\,\mathbf G \,\llbracket\mathbf u\rrbracket}_{\mathbf u^{\text{enr}}(\mathbf x)}$$

donde:

- `N_i(x)` son las funciones de forma estándar del CST (coordenadas baricéntricas; ver [`Tri3`](Tri3.md) §7).
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

donde `B_std` es la matriz `B` del CST estándar (Voigt 2D, 3×6, constante, ver [`Tri3`](Tri3.md) §8) y `B^φ(x)` es la matriz Voigt construida con `∇φ`:

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
con `T_coh` el operador del material cohesivo (en frame local; ver spec [`CohesiveDamageIsotropic`](CohesiveDamageIsotropic.md)).

### 4. Equilibrio — forma fuerte

Balance estático en cada subdominio (`Ω⁻`, `Ω⁺`) más continuidad de tracciones a través de `Γ_d`:

$$\nabla\cdot\boldsymbol\sigma + \mathbf b = \mathbf 0 \quad \text{en } \Omega \setminus \Gamma_d, \qquad \boldsymbol\sigma\cdot\mathbf n = \mathbf t \quad \text{en } \Gamma_d.$$

### 5. Equilibrio — forma débil (residuales acoplados)

Aplicando el principio de trabajos virtuales con variaciones independientes `δd` y `δ[[u]]`, se obtienen **dos residuales acoplados**:

$$\mathbf R^d = \int_{\Omega_e} \mathbf B_{\text{std}}^\top\,\boldsymbol\sigma\,d\Omega \;-\; \mathbf f^{\text{ext}}_d = \mathbf 0$$

$$\mathbf R^{\llbracket u\rrbracket} = -\int_{\Omega_e} \mathbf G^\top\,(\mathbf B^\varphi)^\top\,\boldsymbol\sigma\,d\Omega \;+\; l_d\,\mathbf t = \mathbf 0$$

con `l_d` la longitud efectiva de `Γ_d` dentro del elemento (ver §11). El primer residual se ensambla al sistema global; el segundo se resuelve **localmente** dentro del elemento (ver §10, condensación estática).

---

## Formulación numérica (FEM)

### 6. Discretización del CST padre

Idéntica a [`Tri3`](Tri3.md): tres nodos en orden antihorario, mapeo isoparamétrico desde el triángulo de referencia con vértices `(0,0)`, `(1,0)`, `(0,1)`. Jacobiano constante; aborta con error si `det J ≤ tol`.

**Identificación del nodo solitario.** Cuando `Γ_d` corta dos lados del triángulo, exactamente un nodo queda aislado en `Ω⁺` (lado al que apunta `n`); los otros dos quedan en `Ω⁻`. El nodo solitario es el **único** que tiene `φ = 1`; los otros dos tienen `φ = 0`. La identificación se hace al activar el elemento, comparando el producto `(x_i − x_c)·n` para cada nodo (`x_c` = centroide), donde el signo del producto determina el lado.

### 7. Funciones de forma + función `φ`

`N_i(x)` idénticas a [`Tri3`](Tri3.md) §7. Para `φ`:

$$\varphi(\mathbf x) = N_{i^*}(\mathbf x)$$

donde `i*` es el índice del nodo solitario. Es lineal y `C^0` continua, vale 1 en `i*` y 0 en los otros dos. Sus derivadas son `(∂N_{i*}/∂x, ∂N_{i*}/∂y)`, **constantes** sobre el elemento.

### 8. Matrices B estándar + B^φ

`B_std` es la matriz B clásica del CST (constante, 3×6); ver [`Tri3`](Tri3.md) §8.

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

**Recuperación local** tras converger el Newton global: `[[u]]_committed += Δ[[u]]` con la fórmula de arriba aplicada a `Δd = d_n+1 − d_n`. Esta recuperación se hace en `commit_state()`, no en `compute_element_state()`.

### 11. Activación: criterio Rankine + dirección `n` + `l_d` del Cap. 6

**Criterio de activación**. Al inicio de cada paso de carga (no dentro del Newton — ver §14), el elemento intacto evalúa su tensión principal mayor `σ_I` en el centroide:

$$\sigma_I = \frac{\sigma_{xx} + \sigma_{yy}}{2} + \sqrt{\Bigl(\frac{\sigma_{xx} - \sigma_{yy}}{2}\Bigr)^2 + \sigma_{xy}^2}$$

Si `σ_I > σ_t0` (resistencia a tracción del material cohesivo), se activa:

- `n` = autovector unitario de `σ_I` (perpendicular a la dirección de tracción máxima, criterio Rankine).
- `s` = ortogonal a `n` orientada por regla de la mano derecha con `+z` saliendo del plano (consistente con Reglas.md §5).
- `Γ_d` pasa por el centroide del CST (decisión cerrada: la posición no es DOF; se fija en el activación, igual que `n`).
- Identificación del nodo solitario `i*` (§6).

**Longitud efectiva `l_d` — Cap. 6 de Retama (2010)**. Fórmula cerrada:

$$\boxed{\;l_d = \frac{A_e}{h}\,\cos(\theta - \alpha)\;}$$

donde:
- `A_e` es el área del triángulo.
- `h` es la altura desde el nodo solitario `i*` al lado opuesto.
- `θ` es el ángulo de `Γ_d` respecto a `x` global (i.e. `θ = atan2(s_y, s_x)` con `s = (s_x, s_y)` la tangente a `Γ_d`).
- `α` es el ángulo del lado opuesto al nodo solitario respecto a `x` global.

Esta fórmula es **decisión cerrada en el ADR 0010** (§6); no se relitiga aquí. La alternativa ingenua `l_d = A_e/h` (independiente de `θ`) produce stress locking cuando la grieta no es paralela al lado opuesto; ver fase 3 del ADR 0010 para test numérico.

### 12. Estado de la discontinuidad — `DiscontinuityState`

Dataclass paralela a `ElementState` (ADR 0010 §4), específica de la discontinuidad. Vive en `solidum/core/discontinuity_state.py`:

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
    cohesive_state_committed: dict  # estado interno del CohesiveMaterial (committed)
    cohesive_state_trial: dict      # estado interno del CohesiveMaterial (trial)
```

Semántica **trial / commit** consistente con `ElementState`: `jump_trial` y `cohesive_state_trial` evolucionan en cada iteración del Newton global; al converger el paso, `commit_state()` copia `trial → committed` y aplica la recuperación local de `Δ[[u]]` (§10).

### 13. Cuadratura

Un único punto de Gauss central, peso `A_e·t` — idéntico a `Tri3`. Las matrices `B_std`, `B^φ` y la deformación del bulk son constantes; la integración es exacta. El término cohesivo `l_d·t` se evalúa directamente sobre `Γ_d` (no requiere cuadratura porque `t` es constante a lo largo de `Γ_d` por la cinemática KOS — el salto es constante).

### 14. Activación al principio del paso — interacción con el solver

ADR 0010 §5 fija que la activación se evalúa **al inicio de cada paso de carga, no dentro del Newton**. Razón: dentro del Newton, `σ_I` oscila por el predictor lineal y puede cruzar/descruzar el umbral varias veces antes de converger — *chattering*.

**Implementación propuesta** (a validar — punto abierto en §Diálogo): añadir al contrato de `Element` un método opcional `prepare_step(U_committed)` que el `NonlinearSolver`/`ArcLengthSolver` invoca **una vez por paso** antes del primer Newton. La base lo hace no-op (la mayoría de elementos no necesitan preparación); `CST_Embedded2D` lo sobreescribe para evaluar la activación con el estado convergido del paso anterior. Una vez activado, el elemento queda agrietado **irreversiblemente** (no revierte si el siguiente paso lo descarga).

---

## Caveats numéricos

- **Singularidad de `K_{[[u]][[u]]}` en Modo-I puro saturado**. Cuando `ω → 1` en el cohesivo, la rigidez tangente cohesiva `T_coh` se aproxima a la diagonal `[(1−DAMAGE_MAX)·K_e, 0]` (la componente tangencial es cero por construcción del Modo-I). La suma con la contribución del bulk `G^T·(B^φ)^T·C_e·B^φ·G` introduce rigidez en la dirección tangencial vía el bulk; la matriz 2×2 sigue siendo invertible. Verificable analíticamente.
- **Activación al principio del paso**. Si el incremento es muy grande, `σ_I` puede saltar muy por encima de `σ_t0` antes de la activación, generando un Δ[[u]] grande en el primer Newton. Mitigación: pasos suficientemente pequeños o `ArcLengthSolver` cerca del pico.
- **Bulk descarga elástica**. Por la cinemática KOS, el bulk evalúa `ε^bulk = B_std·d − B^φ·G·[[u]]` — descarga conforme `[[u]]` crece. Esto es físicamente correcto (la disipación va a `Γ_d`) y consistente con la discrete approach de Retama 2010. El bulk **no acumula daño** propio (out_of_scope fase 1).
- **`Γ_d` paralelo a un lado**. Caso degenerado: `cos(θ−α) = 1`, `l_d = A_e/h`. La fórmula del Cap. 6 sigue siendo correcta, pero numéricamente el nodo solitario debe identificarse con tolerancia para evitar oscilación entre dos nodos casi en el mismo plano.
- **`Γ_d` casi vertical en CST con jacobiano malacondicionado**. Sin patología nueva: hereda los caveats del CST padre. El test `det J > tol` de `Tri3` cubre el caso.
- **Stress locking**. La fase 3 del ADR 0010 valida numéricamente que `l_d = (A_e/h)·cos(θ−α)` evita el locking que aparece con `l_d = A_e/h` ingenuo. Test cubierto en `tests/test_ld_chapter_6_validation.py`: con `Γ_d` oblicua respecto al lado opuesto al nodo solitario (`cos(θ−α) = √3/2` en el setup), la versión correcta produce una apertura `[[u_n]]` ≈ 5.45 % mayor que la versión ingenua para una misma carga, en rama de softening exponencial activo. La versión ingenua "infla" la rigidez cohesiva (`l_d · T_coh`) por un factor `1/cos(θ−α) > 1`, restringe la apertura y el elemento se comporta como si fuera más rígido — eso es el *stress locking*. Cuando `Γ_d` es paralela al lado opuesto (`cos(θ−α) = ±1`), las dos fórmulas coinciden exactamente y no hay locking.

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
      desc: "Nombre (en CohesiveMaterialRegistry) del material cohesivo que gobierna Γ_d cuando se activa." }
  - { name: activation_criterion, type: str, required: false, default: rankine,
      desc: "Criterio de activación: 'rankine' (única opción fase 1)." }

material_contract:
  bulk:
    signature: "compute_state(ε, state) -> (σ, C_tangent, state')"
    strain_kind: "Voigt 2D [ε_xx, ε_yy, γ_xy]"
    accepted: ["Elastic2D"]      # fase 1: bulk elástico; otros → out_of_scope
  cohesive:
    signature: "compute_traction([[u]], state) -> (t, T_tan, state')"
    jump_kind: "frame local — [[u]] = ([[u_n]], [[u_s]]) en ejes (n, s) de Γ_d"
    accepted: ["CohesiveDamageIsotropic"]  # cualquier CohesiveMaterial con JUMP_DIM=2

conventions:
  sign: "tracción interna positiva; convención stress-resultant del proyecto (Reglas.md §5)."
  voigt: "[xx, yy, xy] con γ_xy = 2·ε_xy."
  node_orientation: "antihorario (det J > 0)."
  frame_local_discontinuity: "n perpendicular a σ_I al activarse; s = (−n_y, n_x) (RHR con +z saliendo del plano)."
  activation: "irreversible; evaluada al principio de cada paso con el estado convergido previo."

validity:
  - "Pequeños desplazamientos y deformaciones (no large strain)."
  - "Bulk elástico (Elastic2D) y cohesivo Modo-I (JUMP_DIM = 2)."
  - "Γ_d cortando exactamente dos lados del triángulo (configuración geométrica estándar)."

out_of_scope:
  - "Bulk no lineal (J2, daño continuo, Drucker-Prager): la disipación debe localizarse en Γ_d, no co-existir con disipación de bulk."
  - "Modo mixto I-II en el cohesivo (fase G del ADR 0010)."
  - "Reorientación de n tras la activación (tracking no trivial, fase F)."
  - "Múltiples discontinuidades por elemento."
  - "Contacto unilateral en compresión (cierre rígido)."
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

  arch:
    - name: registro_en_ElementRegistry
      expect: "CST_Embedded2D registrado con kind='element'."

    - name: STRAIN_DIM_y_DOF_NAMES_correctos
      expect: "STRAIN_DIM=3, DOF_NAMES=['ux','uy'], N_NODES=3, N_INTEGRATION_POINTS=1."

    - name: prepare_step_es_no_op_en_Element_base
      expect: "Element.prepare_step existe y no hace nada por defecto; CST_Embedded2D lo sobreescribe."

references:
  - "Retama Velasco, J. (2010). Formulation and Approximation to Problems in Solids by Embedded Discontinuity Models. Tesis Doctoral, UNAM. **Caps. 2 (cinemática KOS), 5 (formulación variacional discreta), 6 (l_d correcto), 7 (condensación estática)**."
  - "Jirásek, M. (2000). Comparative study on finite elements with embedded discontinuities. CMAME 188, 307-330. — Clasificación SOS/KOS/SKON."
  - "Oliver, J., Huespe, A.E., Sánchez, P.J. (2006). A comparative study on finite elements for capturing strong discontinuities. CMAME 195, 4732-4752."
  - "Linder, C., Armero, F. (2007). Finite elements with embedded strong discontinuities. IJNME 72, 1391-1433."
  - "ADR 0010 — Discontinuidades interiores embebidas (hoja de ruta; fase 2 = este elemento)."
  - "Spec [CohesiveDamageIsotropic](CohesiveDamageIsotropic.md) — material cohesivo consumido en Γ_d."
  - "Spec [Tri3](Tri3.md) — CST padre del elemento; mismo bulk en estado intacto."
```

---

## Implementación

- Archivo: `solidum/elements/solid_2d/embedded_cst.py`
- Clase: `CST_Embedded2D` (registrada en `ElementRegistry`)
- Estado de la discontinuidad: `solidum/core/discontinuity_state.py` (`DiscontinuityState`)
- Hook de paso en clase base: `Element.prepare_step(U_committed)` (no-op por defecto)
- Integración en solvers:
  - `NonlinearSolver`: invoca `self.assembler.prepare_all_steps(U_current)` al inicio de cada paso, antes del Newton.
  - `ArcLengthSolver`: idem, antes del predictor.
- Parser YAML: nueva sección `cohesive_materials:` + campo `cohesive_material:` por elemento.
- Tests: `tests/test_cst_embedded.py` — 24 tests cubriendo `acceptance` completo (verification, specific, arch) + integración end-to-end del YAML.
- Notas de traducción:
  - Estado intacto (`discontinuity_state is None`) → bit-exact con `Tri3` (mismo `B`, mismo `_compute_integrands`).
  - Estado agrietado: Newton local sobre `[[u]]` hasta `R^{[[u]]} = 0`, después condensación estática local cerrada (`K_jj` es 2×2). Tolerancia interna: `_LOCAL_JUMP_RTOL = 1e-10`, `_LOCAL_JUMP_MAX_ITER = 30`.
  - `B^φ` se obtiene como la submatriz `B_std[:, 2·i*:2·i*+2]` correspondiente al nodo solitario (consistente con `φ = N_{i*}`).
  - Activación: `n` apunta hacia el nodo solitario (convención fija); si el autovector de `σ_I` deja 2 nodos en lado positivo, se invierte el signo de `n`.
  - `_compute_ld` usa `|cos(θ − α)|` (valor absoluto): la longitud es positiva y simétrica bajo `n → −n` (test `ld_invariante_bajo_reflejo_de_n`).
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
