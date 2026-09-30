# Esquemas de Solución

Este capítulo cubre los **solvers estáticos** del catálogo: lineal, no lineal (Newton-Raphson), arc-length (Crisfield cilíndrico) y control indirecto de desplazamiento; al final, el arc-length por disipación de energía del módulo de usuario `dissipation_control`. Los solvers de análisis modal, dinámico (Newmark, HHT-α, central differences) y de análisis en frecuencia / espectral (Harmonic, Response Spectrum) se documentan en el siguiente capítulo, "Análisis Dinámico".

## `LinearSolver` — solución directa lineal

Resuelve $\mathbf K \cdot \mathbf U = \mathbf F$ con una sola factorización del sistema (el backend lo elige Solidum: Cholesky, Pardiso o LU según la matriz y lo instalado). Las condiciones de Dirichlet se imponen por eliminación directa (ADR 0004): el sistema entregado al backend algebraico es ya el reducido $\mathbf K_{\text{red}} \cdot \mathbf U_{\text{libre}} = \mathbf F_{\text{red}}$.

- **Cuándo usarlo**: problemas estrictamente lineales (todos los materiales con tangente constante, sin grandes desplazamientos, sin contacto).
- **Parámetros**: `linear_algebra` (default `auto`; ver *Solver* en el capítulo del YAML).
- **Comprobaciones automáticas** (ADR 0019): antes de resolver, que el modelo esté suficientemente apoyado; después, que la solución cumpla el equilibrio y que la matriz no sea singular. Si algo falla, el análisis se detiene con una explicación en vez de devolver resultados absurdos (ver *Diagnóstico de Problemas*). Con `linear_algebra: iterative` la matriz no se factoriza y no se puede comprobar si es singular: Solidum lo avisa en el log, una vez por análisis.
- **No aplica**: cualquier no-linealidad material (daño, plasticidad, cable) o geométrica (corotacional). Producirá resultados inconsistentes.

```yaml
solver:
  type: LinearSolver
```

## `NonlinearSolver` — Newton-Raphson incremental con paso adaptativo

Bucle externo de incrementos del factor de carga $\lambda \in [0, 1]$; bucle interno Newton-Raphson sobre $\mathbf K_t \cdot \Delta\mathbf U = \mathbf R = \lambda\,\mathbf F_{\text{ext}} - \mathbf F_{\text{int}}$.

**Criterio de convergencia dual** (se exigen los dos a la vez):

$$\lVert\mathbf T^\top\mathbf R\rVert \le \texttt{atol\_force} + \texttt{rtol\_force}\cdot\max(\lVert\mathbf F_{\text{ext}}\rVert, \lVert\mathbf F_{\text{int}}\rVert), \qquad \lVert\Delta\mathbf U\rVert \le \texttt{atol\_disp} + \texttt{rtol\_disp}\cdot\lVert\mathbf U\rVert$$

donde $\mathbf T^\top\mathbf R$ es el residuo de equilibrio en los grados de libertad libres: las reacciones de los apoyos no cuentan y, con restricciones multipunto, el residuo de cada esclavo se suma a sus maestros multiplicado por su coeficiente. Las tolerancias se configuran en el bloque `convergence` (capítulo «Archivo de Entrada YAML»).

**Desplazamientos prescritos**: los valores de Dirichlet no nulos y el término independiente $g$ de las restricciones multipunto se aplican proporcionalmente al factor de carga, como $\lambda\,g$, igual que las cargas.

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

Traza curvas de equilibrio con puntos límite de carga (snap-through, colapso), controlando simultáneamente desplazamientos y factor de carga $\lambda$. **No sigue los retrocesos (snap-back) por localización**, en los que el desplazamiento también disminuye tras el pico: para ellos, `IndirectDisplacementSolver`.

**Esquema**:

- **Predictor tangente**: $\Delta\mathbf U_t = \mathbf K_t^{-1} \cdot \mathbf F_{\text{ext}}^{\text{ref}}$; $\Delta\lambda = \pm\,dl / \lVert\Delta\mathbf U_t\rVert$. Signo elegido por proyección con el incremento previo.
- **Corrector iterativo**: en cada iteración resuelve dos sistemas ($\mathbf K \cdot \Delta\mathbf U_R = \mathbf R$ y $\mathbf K \cdot \Delta\mathbf U_t = \mathbf F_{\text{ext}}$) y aplica la restricción cuadrática cilíndrica $\lVert\delta\mathbf U\rVert^2 = dl^2$.
- **Primer paso**: se declara con `initial_dlambda`, como **fracción de la carga de referencia** (0.1 si no se declara, es decir, el 10 % de la carga). El solver lo convierte en longitud de arco con una resolución elástica, $dl_1 = \Delta\lambda_1\,\lVert\mathbf K^{-1}\mathbf F_{\text{ext}}^{\text{ref}}\rVert$, así que no hace falta conocer de antemano los desplazamientos del modelo y el valor significa lo mismo en cualquier sistema de unidades. La alternativa `initial_dl` fija la longitud de arco directamente, en unidades de desplazamiento; sólo tiene sentido si se conoce la escala del problema, y no puede declararse junto con `initial_dlambda`.
- **Métrica del arco**: $\lVert\delta\mathbf U\rVert$ se mide sobre todos los grados de libertad, también los esclavos de las restricciones multipunto y los prescritos. Un desplazamiento prescrito grande frente a la deformación consume casi toda la longitud de arco y el paso se comporta como un control por carga. Un modelo con nodos atados recorre una secuencia de pasos distinta de la del modelo con los nodos fundidos: todos los puntos son de equilibrio, y con plasticidad o daño los resultados difieren dentro del error de discretización de la trayectoria de carga, como al cambiar el tamaño de paso.
- **Auto-ajuste de `dl`**: se agranda en convergencias rápidas (hasta `dl_max_factor` veces la del primer paso), se reduce en lentas, biseca al fallar.
- **Llegada a `max_lambda`**: todos los pasos siguen la curva. Si uno la cruza, se repite desde el estado anterior hasta llegar exactamente a `max_lambda`. El solver nunca impone la carga máxima sin haber recorrido la curva hasta ella; si hay un colapso antes, el trazado se detiene por `max_steps` y lo avisa.

```yaml
solver:
  type: ArcLengthSolver
  max_iter: 15
  max_lambda: 1.0
  initial_dlambda: 0.05      # primer paso: 5 % de la carga de referencia
  max_steps: 200
  convergence:
    rtol_force: 1.0e-5
    rtol_disp: 1.0e-5
```

**Cuándo usarlo**: problemas con softening pronunciado (daño, post-pandeo, snap-through de cúpulas), o cuando `NonlinearSolver` diverge cerca de un punto límite. Más caro por paso (dos `spsolve` por iteración) pero indispensable en estos casos.

**Referencia**: Crisfield, "A fast incremental/iterative solution procedure that handles snap-through" (*Computers & Structures*, 1981).

## `IndirectDisplacementSolver` — control indirecto de desplazamiento

Variante de `ArcLengthSolver` que controla, en vez de la longitud de arco de todos los grados de libertad, el incremento de **una magnitud elegida**: el alargamiento de la zona que se ablanda o la apertura de una grieta (CMOD), como se controla en laboratorio un ensayo de fractura. La carga $\lambda$ queda libre y puede bajar, así que sigue los **retrocesos** (snap-back) por localización, que `ArcLengthSolver` no sigue (tampoco el arc-length por disipación del módulo de usuario `dissipation_control`).

$$\mathbf c^\top \Delta\mathbf U = \Delta l$$

- **`control`** (obligatorio): la magnitud, como lista de términos `{node, dof, coef}`. El alargamiento entre los nodos 6 y 7 en $x$ es `[{node: 7, dof: ux, coef: 1.0}, {node: 6, dof: ux, coef: -1.0}]`. Los grados de libertad deben ser libres, y la magnitud tiene que crecer al cargar; si no, el solver lo explica con un error.
- **Primer paso**: `initial_dlambda`, fracción de la carga de referencia, como en `ArcLengthSolver`.
- **Sin crecimiento del paso por omisión** (`dl_max_factor: 1.0`). Con ablandamiento local, un paso que sobrepasa el pico más que la imperfección que localiza el daño puede converger a otro equilibrio, con más zonas ablandadas de las que el camino físico produciría. Esto le pasa a cualquier solver, también al control por desplazamiento; la defensa es un paso pequeño cerca del pico.

```yaml
solver:
  type: IndirectDisplacementSolver
  initial_dlambda: 0.02
  max_lambda: 1.2
  max_steps: 300
  control:
    - {node: 7, dof: ux, coef: 1.0}
    - {node: 6, dof: ux, coef: -1.0}
```

**Referencia**: R. de Borst, "Computation of post-bifurcation and post-failure behavior of strain-softening solids" (*Computers & Structures*, 1987).

## Módulo de usuario `dissipation_control`: arc-length por disipación de energía

No es un solver estándar según el criterio de admisión del proyecto: ningún código general establecido implementa esta restricción de forma nativa. Por eso no forma parte del programa principal: vive en el **módulo de usuario** `dissipation_control` (`solidum/user/dissipation_control/`, ADR 0020). `import solidum` no lo carga; el YAML que lo use debe declararlo:

```yaml
user_modules: [dissipation_control]
```

Sin esa línea, la lectura se detiene con un error que nombra el solver desconocido y lista los módulos de usuario disponibles. En Python, `solidum.load_user_module("dissipation_control")` o `from solidum.user.dissipation_control import DissipationArcLengthSolver`.

### `DissipationArcLengthSolver` — arc-length por disipación de energía

Variante de `ArcLengthSolver` que controla la **disipación incremental de energía** por paso en lugar de la longitud de arco euclídea. Es la herramienta indicada cuando la curva $\mathbf U$–$\lambda$ tiene tramos casi verticales, donde la restricción cilíndrica se vuelve mal condicionada.

**Diferencia esencial**: la restricción de paso

$$g(\Delta\mathbf U, \Delta\lambda) = \tfrac{1}{2}\left(\lambda_n \mathbf F \cdot \Delta\mathbf U - \Delta\lambda\, \mathbf F \cdot \mathbf U_n\right) = \tau$$

es **lineal** en $(\Delta\mathbf U, \Delta\lambda)$, frente a la cuadrática del cilíndrico. El corrector tiene por tanto una sola raíz: desaparecen la selección de raíz por menor ángulo y la patología de raíces imaginarias.

- **Switching automático**: arranca en modo cilíndrico (la restricción por disipación es idénticamente nula mientras $\lambda_n = \mathbf U_n = 0$) y conmuta a disipación en cuanto detecta disipación neta sobre el umbral. La conmutación es automática en ambos sentidos.
- **Parámetros**: todos los de `ArcLengthSolver`, más `initial_tau` (obligatorio) y la familia `tau_grow_factor`, `tau_max_factor`, `tau_shrink_factor`, `tau_grow_iter_threshold`, `tau_shrink_iter_threshold`, `dissipation_threshold`.
- **Cuándo usarlo**: softening de daño continuo (1D/2D) con rama post-pico pronunciada, donde el cilíndrico converge mal.
- **Limitación medida**: tampoco sigue el retroceso de una barra con daño localizado. Su predictor degenera en el régimen de descarga del material y el paso que cruza el pico puede saltar a otro equilibrio (deuda #27). Para retrocesos, `IndirectDisplacementSolver`.
- **Limitación conocida**: **no** resuelve el caso de discontinuidad embebida (`CST_Embedded2D`, módulo de usuario `discontinuities`) con penalty cohesivo rígido. La activación discreta del criterio de Rankine produce un salto en las fuerzas internas que la regla de signo del predictor (producto escalar con el paso anterior) no maneja. Es una limitación documentada del solver, no un error de configuración. El diagnóstico es anterior a la corrección de la tangente cohesiva del 2026-09-23 y no se ha vuelto a comprobar; con desplazamiento impuesto, `NonlinearSolver` sí recorre hoy esa rama. Un modelo que combine los dos módulos los declara juntos: `user_modules: [discontinuities, dissipation_control]`.
- **Puntos abiertos**: la spec del solver registra tres cuestiones físicas pendientes de decisión (el signo del predictor en modo de disipación, el piso de la escala del umbral de conmutación y el alcance en plasticidad). Mientras no se cierren, conviene limitar su uso al daño continuo validado.

```yaml
user_modules: [dissipation_control]
solver:
  type: DissipationArcLengthSolver
  max_lambda: 1.0
  initial_dlambda: 0.05
  initial_tau: 1.0e-4
  max_steps: 300
```

**Referencia**: Gutiérrez, "Energy release control for numerical simulations of failure in quasi-brittle solids" (*Communications in Numerical Methods in Engineering*, 2004); switching a la Verhoosel et al. (2009).

La formulación completa está en la spec del solver, en el apéndice «Módulos de usuario» del Manual de Referencia.
