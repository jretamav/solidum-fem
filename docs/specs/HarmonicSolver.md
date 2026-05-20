# HarmonicSolver — respuesta forzada armónica en el dominio de la frecuencia (ADR 0009 fase 6)

> Spec de **solver nuevo**. Resuelve `(−ω²M + iωC + K)·û = F̂` para un barrido de frecuencias y devuelve la amplitud compleja `û(ω)`. Análisis lineal en frecuencia, sin transitorio.

---

## Especificación física

### 0. Descripción general

Análisis **lineal** de respuesta forzada armónica en estado estacionario. La excitación `F(t) = Re{F̂·e^{iωt}}` produce una respuesta `u(t) = Re{û(ω)·e^{iωt}}` también armónica con la misma frecuencia. El solver calcula `û(ω)` para un barrido de frecuencias. Útil para:

- Funciones de transferencia (FRFs) entrada→salida.
- Diagnóstico de resonancias previo a análisis transitorio caro.
- Respuesta forzada en estructuras con cargas periódicas (maquinaria rotante, sismo armónico, viento turbulento descompuesto en armónicos).

Hipótesis cinemáticas y material: lineales (`K` y `M` constantes). No se modelan grandes desplazamientos, plasticidad ni daño — quien quiera no linealidad debe usar `NewtonNewmarkSolver` con excitación armónica.

### 1. Ecuación de equilibrio resuelta

$$\big(-\omega^2\mathbf{M} + i\omega\mathbf{C} + \mathbf{K}\big)\,\hat{\mathbf{u}}(\omega) = \hat{\mathbf{F}}$$

para cada `ω` del barrido. `C = α·M + β·K` Rayleigh.

### 2. Condiciones iniciales / de contorno

- Apoyos Dirichlet constantes (eliminación directa, ADR 0004).
- `F̂` ∈ ℂ^ndof. Real puro si la carga lleva fase cero; complejo si lleva fase.
- **No** hay condiciones iniciales — el estado estacionario es independiente de ellas (el transitorio se descarta por definición).
- MPC no soportado en esta fase.

### 3. Salidas físicas

`HarmonicResult` (inmutable):

- `omega[n_omega]` frecuencias evaluadas (rad/s).
- `frequencies_hz` propiedad derivada.
- `u_complex[ndof, n_omega]` amplitud compleja del desplazamiento.
- `F_complex[ndof]` amplitud compleja de la carga (informativa).
- Métodos auxiliares: `.amplitude()` → `|û|`, `.phase()` → `arg(û)`.

---

## Formulación numérica

### 4. Esquema operativo

```
Ensamblar K, M.
C = α·M + β·K (Rayleigh).
Reducir por Dirichlet: K_red, M_red, C_red, F_red.
Para cada ω en el barrido:
    Z(ω) = -ω²·M_red + iω·C_red + K_red   [matriz compleja]
    û_red(ω) = Z(ω)⁻¹ · F_red             [spsolve complejo]
    û(ω) = T · û_red(ω) [expand a globales]
```

### 5. Predictor / corrector

No aplica — un solve directo por frecuencia.

### 6. Criterio de convergencia

No aplica — esquema no iterativo. La precisión depende del backend
algebraico (`scipy.sparse.linalg.spsolve` LU complejo).

### 7. Imposición de Dirichlet y MPC

Eliminación directa idéntica a `NewmarkSolver`. En DOFs prescritos
`û = 0` (el apoyo no oscila). MPC fuera de alcance.

### 8. Backend algebraico

Factorización LU compleja por frecuencia vía
`scipy.sparse.linalg.spsolve(Z.tocsc(), F_red)`. **No hay caché**: `Z`
cambia con `ω`. Coste dominante del análisis es `O(n_omega ·
coste_factorizacion(Z))`.

Posible optimización futura (fuera de esta fase): superposición modal
— ortogonalizar la respuesta sobre los primeros `n_modes` modos del
problema generalizado, lo que reduce el barrido a `n_omega` operaciones
escalares. Justifica el dispositivo cuando `n_modes ≪ n_dof` y el
operador es simétrico.

### 9. Adaptatividad / control de paso

Espaciado del barrido controlado por el usuario:

- `scale="linear"`: `np.linspace(omega_min, omega_max, n_omega)`.
- `scale="log"`: `np.geomspace(omega_min, omega_max, n_omega)` —
  apropiado para FRFs sobre rangos de varias décadas.
- Vector explícito vía `omega=...` — usuarios con criterios
  específicos (densificar alrededor de resonancias conocidas, etc.).

No hay refinamiento adaptativo de la malla de frecuencia.

### 10. Caveats numéricos

- **Resonancia exacta sin amortiguamiento**: `Z(ω_n)` es singular —
  `spsolve` puede fallar o devolver valores enormes. Mitigación:
  evitar `ω` coincidente con un modo, o añadir amortiguamiento
  Rayleigh.
- **Mal condicionamiento cerca de resonancia**: la factorización LU
  pierde precisión a medida que `Z` se aproxima a la singularidad. Es
  intrínseco al método directo; superposición modal es la alternativa
  estable.
- **Cargas con fase**: `F̂` compleja requiere construcción desde código
  Python; YAML no soporta tipos complejos nativos. Si el usuario YAML
  necesita fase, debe construir el dominio programáticamente.
- **No analiza estabilidad** del estado estacionario (no aplica al
  problema lineal — el estado estacionario siempre existe y es único
  para `C` definida positiva y `K` regular).

---

## Contrato de implementación

```yaml
name: HarmonicSolver
kind: solver
status: validated

interface:
  yaml_type: HarmonicSolver
  output: HarmonicResult

parameters:
  - { name: omega, type: array-like | null, required: false, desc: "vector explícito (rad/s); excluye omega_min/omega_max/n_omega/scale" }
  - { name: omega_min, type: float, required: false, desc: "límite inferior del barrido (rad/s)" }
  - { name: omega_max, type: float, required: false, desc: "límite superior del barrido (rad/s)" }
  - { name: n_omega, type: int, required: false, desc: "número de puntos del barrido, default 100" }
  - { name: scale, type: str, required: false, desc: "'linear' (default) o 'log'" }
  - { name: F_amplitude, type: ndarray | null, required: false, desc: "amplitud compleja shape (ndof,); default ceros; YAML inyecta desde point_loads como real" }
  - { name: rayleigh, type: dict | null, required: false, desc: "C = α·M + β·K; sin amortiguamiento ⇒ pico singular en resonancia" }
  - { name: lumping, type: str, required: false, desc: "'consistent' (default) o 'lumped'" }

requirements:
  - "Materiales con `density` declarada (ADR 0008)."
  - "Problema lineal — `K` y `M` constantes en el barrido."

conventions:
  units:    "heredadas del modelo (ADR 0008). ω en rad/s — coherente con frecuencias de ModalSolver."
  stability: "Lineal en frecuencia — estable para C definida positiva. Mal condicionado en ω → ω_n sin amortiguamiento."

out_of_scope:
  - "Superposición modal — diferida; el barrido directo es más simple y suficiente para los tamaños actuales del catálogo."
  - "MPC (restricciones lineales multipunto) — diferido."
  - "No linealidad (geométrica o material) — usar NewtonNewmarkSolver con excitación armónica."
  - "Cargas con fase declaradas desde YAML — usar construcción desde código."

acceptance:
  verification:
    - name: transfer_function_no_damping
      setup: "1 GDL con K=25, M=1; barrido en ω lejos de resonancia"
      expect: "û_real(ω) = F̂/(K − ω²M); û_imag ≈ 0"
      tol_rel: 1.0e-10
    - name: transfer_function_rayleigh_damping
      setup: "1 GDL con α-mass Rayleigh"
      expect: "û = F̂/(K − ω²M + iωc) complejo"
      tol_rel: 1.0e-9
  specific:
    - name: resonance_peak_location
      setup: "1 GDL con ξ=2% Rayleigh; barrido fino alrededor de ω_n"
      expect: "pico en ω = ω_n·√(1 − 2ξ²)"
      tol_rel: 0.01
    - name: matches_newmark_steady_state
      setup: "Newmark con F(t)=F₀cos(ωt), suficiente t_end para descartar transitorio (8τ)"
      expect: "amplitud Newmark estacionaria = |û_harmonic|"
      tol_rel: 0.02
    - name: log_sweep_geometric_spacing
      setup: "scale='log', omega_min=1, omega_max=1000, n_omega=4"
      expect: "omega = geomspace(1, 1000, 4) = [1, 10, 100, 1000]"
      tol_rel: 1.0e-12

references:
  - "Bathe, K.-J. (2014). *Finite Element Procedures*, §9.5."
  - "Clough, R. W., Penzien, J. (1993). *Dynamics of Structures*, §4.6."
  - "Chopra, A. K. (2017). *Dynamics of Structures*, cap. 4."
```

---

## Implementación

- Archivo: [solidum/math/solvers/harmonic.py](../../solidum/math/solvers/harmonic.py)
- Clase: `HarmonicSolver` (registrada en `SolverRegistry`)
- Atributo de despacho: `PIPELINE_KIND = "harmonic"` → `solidum.entry.run_harmonic` → `HarmonicResult`.
- Tests: [tests/test_harmonic.py](../../tests/test_harmonic.py) (11 tests verdes).
- Resultado: [`HarmonicResult`](../../solidum/results.py) — dataclass inmutable con `omega`, `u_complex`, `F_complex`, propiedades `frequencies_hz`, `n_omega`, métodos `amplitude()` y `phase()`.

---

## Diálogo

- *2026-05-18*: nuevo `PIPELINE_KIND = "harmonic"` añadido al despacho declarativo de `run_yaml` (regla C, ya aplicada en D2). Tercer valor del literal junto a `"modal"`, `"transient"`, `"static"`. El parser YAML inyecta automáticamente las `point_loads` del bloque estándar como `F_amplitude` real cuando el usuario no la especifica explícitamente — convención coherente con el uso de `external_forces` en análisis estáticos.
