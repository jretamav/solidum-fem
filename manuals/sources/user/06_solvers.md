# Esquemas de Solución

Este capítulo cubre los **solvers estáticos** del catálogo: lineal, no lineal (Newton-Raphson) y arc-length (Crisfield cilíndrico). Los solvers de análisis modal, dinámico (Newmark, HHT-α, central differences) y de análisis en frecuencia / espectral (Harmonic, Response Spectrum) se documentan en el siguiente capítulo, "Análisis Dinámico".

## `LinearSolver` — solución directa lineal

Resuelve $\mathbf K \cdot \mathbf U = \mathbf F$ en un único `spsolve`. Las condiciones de Dirichlet se imponen por eliminación directa (ADR 0004): el sistema entregado al backend algebraico es ya el reducido $\mathbf K_{\text{red}} \cdot \mathbf U_{\text{libre}} = \mathbf F_{\text{red}}$.

- **Cuándo usarlo**: problemas estrictamente lineales (todos los materiales con tangente constante, sin grandes desplazamientos, sin contacto).
- **Parámetros**: `linear_algebra` (default `auto`; ADR 0003 §4).
- **No aplica**: cualquier no-linealidad material (daño, plasticidad, cable) o geométrica (corotacional). Producirá resultados inconsistentes.

```yaml
solver:
  type: LinearSolver
```

## `NonlinearSolver` — Newton-Raphson incremental con paso adaptativo

Bucle externo de incrementos del factor de carga $\lambda \in [0, 1]$; bucle interno Newton-Raphson sobre $\mathbf K_t \cdot \Delta\mathbf U = \mathbf R = \lambda\,\mathbf F_{\text{ext}} - \mathbf F_{\text{int}}$.

**Criterio de convergencia dual**:

$$\max(\text{err}_{\text{disp}}, \text{err}_{\text{force}}) < \text{tol}$$

con $\text{err}_{\text{disp}} = \lVert\Delta\mathbf U\rVert / (\lVert\mathbf U\rVert + \epsilon)$ y $\text{err}_{\text{force}} = \lVert\mathbf R\rVert / \max(\lVert\lambda\,\mathbf F_{\text{ext}}\rVert, \lVert\mathbf F_{\text{int}}\rVert, \epsilon)$.

**Adaptatividad** (`adaptive: true`):

- Si converge en menos de 5 iteraciones, agranda el siguiente paso ($\times 1.5$).
- Si no converge, biseca el paso ($\div 2$). Falla definitivamente si $\Delta\lambda < \texttt{min\_delta\_lambda}$.

```yaml
solver:
  type: NonlinearSolver
  num_steps: 20
  max_iter: 15
  adaptive: true
  convergence:
    rtol_force: 1.0e-5
    rtol_disp: 1.0e-5
```

**Parámetros**: `convergence` (bloque con `rtol_force`, `rtol_disp`, `atol_force_factor`, `atol_disp_factor`; ver ADR 0007), `max_iter`, `num_steps`, `adaptive`, `min_delta_lambda`, `linear_algebra`, `freeze_tangent_after_iter`.

**Limitación**: no puede atravesar puntos límite (snap-through) ni recorrer ramas con derivada negativa de la curva carga-desplazamiento. Para esos casos → `ArcLengthSolver`.

## `ArcLengthSolver` — longitud de arco cilíndrico (Crisfield)

Traza curvas de equilibrio con snap-through, snap-back o pérdida de unicidad de carga, controlando simultáneamente desplazamientos y factor de carga $\lambda$.

**Esquema**:

- **Predictor tangente**: $\Delta\mathbf U_t = \mathbf K_t^{-1} \cdot \mathbf F_{\text{ext}}^{\text{ref}}$; $\Delta\lambda = \pm\,dl / \lVert\Delta\mathbf U_t\rVert$. Signo elegido por proyección con el incremento previo.
- **Corrector iterativo**: en cada iteración resuelve dos sistemas ($\mathbf K \cdot \Delta\mathbf U_R = \mathbf R$ y $\mathbf K \cdot \Delta\mathbf U_t = \mathbf F_{\text{ext}}$) y aplica la restricción cuadrática cilíndrica $\lVert\delta\mathbf U\rVert^2 = dl^2$.
- **Auto-ajuste de `dl`**: se agranda en convergencias rápidas, se reduce en lentas, biseca al fallar.
- **Caso especial `final_step`**: si el predictor sobrepasaría `max_lambda`, fija $\lambda$ a ese valor y resuelve solo desplazamientos (Newton-Raphson puro).

```yaml
solver:
  type: ArcLengthSolver
  max_iter: 15
  max_lambda: 1.0
  initial_dl: 0.05
  max_steps: 200
  convergence:
    rtol_force: 1.0e-5
    rtol_disp: 1.0e-5
```

**Cuándo usarlo**: problemas con softening pronunciado (daño, post-pandeo, snap-through de cúpulas), o cuando `NonlinearSolver` diverge cerca de un punto límite. Más caro por paso (dos `spsolve` por iteración) pero indispensable en estos casos.

**Referencia**: Crisfield, "A fast incremental/iterative solution procedure that handles snap-through" (*Computers & Structures*, 1981).
