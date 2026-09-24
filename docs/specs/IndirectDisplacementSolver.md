# IndirectDisplacementSolver — control indirecto de desplazamiento (de Borst 1987)

> **Spec corta tipo extensión** (Reglas.md §4). Variante de [`ArcLengthSolver`](ArcLengthSolver.md) que **sólo cambia la restricción del paso**: en lugar de fijar la norma del incremento de todos los grados de libertad ($\lVert\Delta\mathbf U\rVert = \Delta l$), fija el incremento de una **combinación elegida** de ellos ($\mathbf c^\top\Delta\mathbf U = \Delta l$). Reusa el predictor tangente, el corrector de Newton compartido (ADR 0015), el primer paso adimensional (deuda #26), la llegada exacta a `max_lambda`, Dirichlet/MPC y el backend algebraico. Sin ADR nuevo.
>
> Motivación (medida el 2026-09-23 al preparar el ejemplo 5 del manual de ejemplos): en una barra con daño que se localiza en un elemento, el retroceso (*snap-back*) empieza en el propio pico. Ahí el incremento de la parte que se descarga apunta casi al revés que el del paso anterior, y ni el `ArcLengthSolver` cilíndrico (su selección de raíz por menor ángulo elige la rama equivocada o no converge) ni el `DissipationArcLengthSolver` (ver deuda #27) siguen la curva. Controlar el alargamiento del elemento que se daña, como un ensayo de fractura controlado por la apertura de la boca de fisura (CMOD), la sigue exacta. Es también el control que el ADR 0010 dejó previsto para la fase 4 (Van Vliet).

---

## Especificación física

### 0. Descripción general

Solver estático no lineal que traza la curva de equilibrio $\mathbf F_\text{int}(\mathbf U) = \lambda\,\mathbf F_\text{ref}$ avanzando en cada paso un incremento fijo $\Delta l$ de una **magnitud cinemática elegida por el usuario**: una combinación lineal de grados de libertad, típicamente el alargamiento de la zona que se ablanda o la apertura de una grieta. La carga $\lambda$ es incógnita y puede bajar, así que la curva se sigue a través de puntos límite y de retrocesos, **siempre que la magnitud controlada crezca monótonamente** a lo largo del camino físico.

### 1. Restricción del paso

$$
g(\Delta\mathbf U) = \mathbf c^\top\Delta\mathbf U - \Delta l = 0,
\qquad
\mathbf c^\top\mathbf U = \sum_k c_k\,u_{(n_k,\,d_k)}
$$

con $\mathbf c$ el vector de coeficientes sobre los grados de libertad controlados. Para el alargamiento entre los nodos $a$ y $b$ en $x$: $c_{b,x} = +1$, $c_{a,x} = -1$.

### 2. Predictor

Con $\mathbf{du}_t = \mathbf K_t^{-1}\mathbf F_\text{ref}$:

$$
\Delta\lambda = \frac{\Delta l}{\mathbf c^\top\mathbf{du}_t}, \qquad \Delta\mathbf U = \Delta\lambda\,\mathbf{du}_t .
$$

**Sin regla de signo**: el signo de $\Delta\lambda$ sale de la propia restricción (negativo en un retroceso, donde la tangente hace bajar la carga mientras la magnitud controlada sigue creciendo). El arco cilíndrico, en cambio, necesita decidir el sentido con el producto escalar con el incremento anterior.

### 3. Corrector

En cada iteración, con $\mathbf{du}_R = \mathbf K^{-1}\mathbf R$ y $\mathbf{du}_t = \mathbf K^{-1}\mathbf F_\text{ref}$:

$$
\delta\lambda = \frac{\Delta l - \mathbf c^\top(\Delta\mathbf U + \mathbf{du}_R)}{\mathbf c^\top\mathbf{du}_t},
\qquad
\delta\mathbf U = \mathbf{du}_R + \delta\lambda\,\mathbf{du}_t .
$$

Es **lineal** en $\delta\lambda$: no hay raíces que elegir ni discriminantes negativos. Tras converger, $\mathbf c^\top\Delta\mathbf U = \Delta l$ se cumple a precisión de máquina.

### 4. Primer paso y tamaño de paso

- **Primer paso adimensional** (deuda #26): `initial_dlambda` (0.1 por omisión) es la fracción de la carga de referencia; $\Delta l_1 = \Delta\lambda_1\,\mathbf c^\top\mathbf K_0^{-1}\mathbf F_\text{ref}$. `initial_dl` fija $\Delta l_1$ explícito, en las unidades de la magnitud controlada.
- **Sin crecimiento por omisión** (`dl_max_factor = 1.0`, frente a 5.0 del padre): la magnitud avanza a ritmo constante, como en un ensayo a velocidad de CMOD constante. La razón es física, ver §5.
- Se conserva la reducción por iteraciones y la bisección del padre.

### 5. No unicidad del problema incremental con ablandamiento

Con ablandamiento local, el problema de un paso (equilibrio con la historia $\kappa$ del paso anterior) admite varios equilibrios: además del físico, en que sólo se ablanda la zona debilitada, otros en que también se ablandan zonas que en el camino continuo se habrían descargado. El Newton converge a uno u otro según la estimación inicial. **Un paso que sobrepasa el pico más que la imperfección que localiza el daño puede converger a otro equilibrio.** Medido en la barra de 10 elementos con uno debilitado un 5 %: con pasos del 2 % o el 5 % de la carga el trazado es exacto; con el 10 %, o con el crecimiento ×5 del padre, se dañan los diez elementos. El control por desplazamiento del extremo (`NonlinearSolver`) hace lo mismo con pasos grandes: no es una propiedad del solver sino del problema discreto. De ahí el valor por omisión sin crecimiento y la recomendación de pasos menores que el margen de la imperfección cerca del pico.

---

## Contrato de implementación

```yaml
name: IndirectDisplacementSolver
kind: solver
status: validated

parent_spec: ArcLengthSolver

interface:
  yaml_type: IndirectDisplacementSolver
  output: SolveResult
  pipeline_kind: static

parameters:
  # Específico de esta variante:
  - { name: control,                  type: list,  required: true,
      desc: "Combinación controlada Σ coef·u(nodo, dof): lista de (nodo, dof, coef) o de {node, dof, coef} (coef = 1 si se omite); nodo = Node o id. Grados de libertad libres; la combinación debe crecer con la carga de referencia" }
  # Heredados del padre, con la semántica de la magnitud controlada:
  - { name: convergence,              type: ConvergenceCriterion, required: false, default: "ConvergenceCriterion()" }
  - { name: max_iter,                 type: int,   required: false, default: 20 }
  - { name: max_lambda,               type: float, required: false, default: 1.0 }
  - { name: max_steps,                type: int,   required: false, default: 100 }
  - { name: initial_dlambda,          type: float, required: false, default: 0.1,
      desc: "Fracción de la carga de referencia del primer paso; Δl₁ = Δλ₁·cᵀK₀⁻¹F_ref. Excluyente con initial_dl" }
  - { name: initial_dl,               type: float, required: false, default: null,
      desc: "Incremento explícito de la magnitud controlada en el primer paso, en sus unidades. Excluyente con initial_dlambda" }
  - { name: dl_grow_factor,           type: float, required: false, default: 1.5 }
  - { name: dl_max_factor,            type: float, required: false, default: 1.0,
      desc: "Sin crecimiento por omisión (ver §5)" }
  - { name: dl_shrink_factor,         type: float, required: false, default: 0.6 }
  - { name: dl_grow_iter_threshold,   type: int,   required: false, default: 4 }
  - { name: dl_shrink_iter_threshold, type: int,   required: false, default: 8 }
  - { name: linear_algebra,           type: str,   required: false, default: "auto" }

requirements:
  - "La magnitud controlada crece monótonamente a lo largo del camino físico (apertura de una grieta, alargamiento de la zona que se ablanda). Si decrece con la carga de referencia al empezar, el solver lo rechaza con un error que lo explica"
  - "Los grados de libertad controlados son libres (sin apoyo ni restricción multipunto). Para controlar un nodo respecto de un apoyo fijo basta con el nodo solo"
  - "Heredados del padre: F_ref no nulo; modelo sin mecanismos (ADR 0019)"

conventions:
  units: "Δl en las unidades de la magnitud controlada (longitud para desplazamientos con coeficientes adimensionales)"
  sign: "la magnitud controlada crece Δl > 0 en cada paso; λ puede crecer o decrecer"

out_of_scope:
  - "Elección automática de la magnitud controlada (control local adaptativo, p. ej. el mayor alargamiento de la malla en cada paso)"
  - "Magnitudes que cambian de signo de crecimiento a lo largo del camino (descarga global): el control exige monotonía"
  - "Análisis dinámico (régimen cuasiestático)"

acceptance:
  verification:
    - name: retroceso_barra_con_dano_localizado
      setup: "Barra de 10 Truss2D con IsotropicDamage1D (α·κ₀w = 0.475), elemento central con κ₀ un 5 % menor; control del alargamiento del elemento débil, initial_dlambda = 0.02"
      expect: "Sólo se daña el elemento débil; los puntos tras el pico están sobre la curva exacta (distancia normalizada < 1e-3); el desplazamiento del extremo retrocede hasta el mínimo exacto 0.616·u_pico"
      tol_rel: 2.0e-3

    - name: rama_descendente_sin_retroceso
      setup: "La misma barra con 2 elementos ((n − 1)·α·κ₀w < 1)"
      expect: "Rama elástica exacta y rama descendente sobre la curva exacta"
      tol_rel: 2.0e-3

    - name: restriccion_satisfecha
      setup: "Barra de 10 elementos del primer caso"
      expect: "El alargamiento controlado crece exactamente Δl por paso"
      tol_rel: 1.0e-8

  specific:
    - name: primer_paso_como_fraccion_de_carga
      setup: "Barra elástica, control del desplazamiento del extremo, initial_dlambda = 0.25"
      expect: "λ₁ = 0.25 exacto; Δλ = 0.25 en cada paso; solución lineal exacta en λ = 1"
    - name: errores_claros
      setup: "control vacío o mal formado, nodo o dof inexistentes, dof restringido, términos que se anulan, magnitud que decrece con la carga"
      expect: "ValueError con un mensaje que explica el problema"
    - name: desde_yaml
      setup: "solver: {type: IndirectDisplacementSolver, control: [{node, dof, coef}, …]}"
      expect: "Se construye y resuelve"

references:
  - "de Borst R. (1987). Computation of post-bifurcation and post-failure behavior of strain-softening solids. *Computers & Structures* 25(2), 211-224. (Control indirecto de desplazamiento.)"
  - "Crisfield M.A. (1991). *Non-linear Finite Element Analysis of Solids and Structures*, Vol. 1, Wiley, cap. 9. (Métodos de continuación y restricciones alternativas.)"
  - "de Borst R., Sluys L.J. (1999). *Computational Methods in Non-linear Solid Mechanics*, TU Delft, §3."
  - "Spec padre: ArcLengthSolver.md. ADR 0015 (corrector de Newton compartido)."
```

---

## Implementación

- Archivo: `solidum/math/solvers/indirect_displacement.py` — `_IndirectProblem` (restricción lineal sobre el `_ArcProblem` del padre) e `IndirectDisplacementSolver` (registrado en `SolverRegistry`, `PIPELINE_KIND = "static"`).
- El padre expone cuatro puntos de extensión que usa esta variante: `_control_measure` (magnitud que fija la restricción: norma en el padre, $\mathbf c^\top\mathbf{du}$ aquí), `_predictor_dlambda`, `_step_problem` y `_on_solve_start` (construye $\mathbf c$ sobre la numeración de ecuaciones). Sin cambio de comportamiento del padre.
- Tests: [`tests/test_indirect_displacement.py`](../../tests/test_indirect_displacement.py) — todos los `acceptance`.

---

## Diálogo

- **2026-09-23** · Formulación propuesta por la IA y validada por el usuario al descubrir, preparando el ejemplo 5 del manual de ejemplos, que ninguno de los dos solvers de arco seguía el retroceso de una barra con daño localizado. Diagnóstico registrado como deuda #27 (arco por disipación) y en §5 (no unicidad del problema incremental).
