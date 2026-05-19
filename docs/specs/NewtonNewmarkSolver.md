# NewtonNewmarkSolver — Análisis dinámico transitorio no lineal (Newmark-β + Newton)

> Variante de [`NewmarkSolver`](NewmarkSolver.md) según la regla de variantes
> (Reglas.md §4). Documenta **solo lo que cambia** respecto al padre. Implementa
> la **fase 4 del ADR 0009**; reusa todas las decisiones arquitecturales del
> ADR padre — sin ADR nuevo.

---

## Qué cambia respecto a `NewmarkSolver`

`NewmarkSolver` integra `M·ü + C·u̇ + K·u = F(t)` asumiendo **K constante**:
factoriza `A_eff = M + γΔt·C + βΔt²·K` una vez y cada paso temporal es una
resolución triangular barata.

`NewtonNewmarkSolver` resuelve el mismo problema cuando **la rigidez no es
constante**: hay materiales con historia (plasticidad, daño) o no linealidad
geométrica. En cada paso temporal se introduce un bucle de **Newton-Raphson**
sobre el residuo dinámico, hasta cumplir el criterio de convergencia.

### Residuo dinámico en cada paso

En el instante $t_{n+1}$ las incógnitas son $(\mathbf u_{n+1}, \dot{\mathbf u}_{n+1}, \ddot{\mathbf u}_{n+1})$
ligadas por los correctores Newmark:

$$\mathbf u_{n+1} = \tilde{\mathbf u}_{n+1} + \beta\Delta t^2\,\ddot{\mathbf u}_{n+1}, \qquad
  \dot{\mathbf u}_{n+1} = \tilde{\dot{\mathbf u}}_{n+1} + \gamma\Delta t\,\ddot{\mathbf u}_{n+1}$$

(con predictores idénticos a la fase 3). El residuo dinámico es:

$$\mathbf R(\ddot{\mathbf u}_{n+1}) \;=\; \mathbf F_\text{ext}(t_{n+1}) \;-\; \mathbf F_\text{int}(\mathbf u_{n+1}) \;-\; \mathbf C\,\dot{\mathbf u}_{n+1} \;-\; \mathbf M\,\ddot{\mathbf u}_{n+1} \;-\; \mathbf F_\text{dir}$$

donde $\mathbf F_\text{int}(\mathbf u)$ es el vector de fuerzas internas no lineales
(plástico, dañado, etc.) que el `Assembler.assemble_non_linear_system` ya
calcula para los solvers estáticos (ADR 0003 + ADR 0006 + ADR 0008).

### Jacobiano dinámico tangente

Derivando $\mathbf R$ respecto a $\ddot{\mathbf u}_{n+1}$ usando los correctores:

$$\mathbf J \;=\; -\frac{\partial \mathbf R}{\partial \ddot{\mathbf u}_{n+1}} \;=\; \mathbf M \;+\; \gamma\Delta t\,\mathbf C \;+\; \beta\Delta t^2\,\mathbf K_\text{t}$$

donde $\mathbf K_\text{t} = \partial \mathbf F_\text{int}/\partial \mathbf u$ es la
rigidez tangente del estado corriente. La iteración es:

$$\mathbf J\,\delta \ddot{\mathbf u} = \mathbf R, \qquad \ddot{\mathbf u}_{n+1}^{(k+1)} = \ddot{\mathbf u}_{n+1}^{(k)} + \delta \ddot{\mathbf u}$$

Los desplazamientos y velocidades se actualizan con los mismos correctores
Newmark tras cada update de $\ddot{\mathbf u}_{n+1}$.

### Amortiguamiento Rayleigh — heredado, constante en el tiempo

$\mathbf C = \alpha\,\mathbf M + \beta_R\,\mathbf K_0$ con $\mathbf K_0$ la
rigidez **elástica de referencia** (la primera ensamblada en $t = 0$, $\mathbf u = 0$).
No se actualiza con la plasticidad: el Rayleigh es modelo aproximado de
amortiguamiento estructural, no consecuencia de la cinemática elastoplástica.
Mantenerlo constante evita un acoplamiento ad-hoc que ningún código de
referencia (Abaqus, ANSYS, OpenSees) realiza por defecto.

### Reducción y commit de estado

- **Reducción por Dirichlet** heredada de la fase 3 (ADR 0004): mismo operador $\mathbf T$ y $\mathbf g$ constantes; el residuo y el jacobiano se proyectan a DOFs libres antes de cada solve.
- **Commit de variables internas**: al converger el paso, se llama `assembler.commit_all_states()` — los estados $(\boldsymbol\varepsilon^p, \alpha, \kappa, d, \ldots)$ pasan de *trial* a *committed*. Mismo patrón que `NonlinearSolver` y `ArcLengthSolver`. Si el paso no converge, los estados trial se descartan (no se llama commit) y se reporta error.

### Convergencia (ADR 0007)

Criterio dual fuerza + desplazamiento, calibrado al primer paso temporal con
la escala del ensamblaje inicial:

- $\|\mathbf R\|/(\text{atol} + \text{rtol}\,\|\mathbf F\|)\,\leq 1$
- $\|\delta \mathbf u\|/(\text{atol} + \text{rtol}\,\|\mathbf u\|)\,\leq 1$

Idéntico al `NonlinearSolver` estático. La calibración se hace una vez al
inicio del análisis; las tolerancias son adimensionales (invariantes bajo
cambio de unidades) gracias a la calibración con `force_scale`/`disp_scale`.

### Factorización por iteración

A diferencia de `NewmarkSolver` (que factoriza $\mathbf J$ una sola vez al
inicio), aquí $\mathbf K_\text{t}$ cambia con cada iteración, así que $\mathbf J$
debe re-factorizarse. Se mantiene el patrón de **Newton modificado opcional**
(ADR 0003 fase 2): si `freeze_tangent_after_iter` es `int`, se factoriza fresco
las primeras N iteraciones del paso y se reusa la factorización en las
siguientes. `None` ⇒ Newton estándar (re-factoriza cada iteración).

### Recuperación del comportamiento lineal

Si todos los materiales del modelo son lineales, $\mathbf K_\text{t} \equiv \mathbf K_0$
constante y el residuo se anula en una iteración de Newton (Δx exacto en el
primer solve). Por tanto **`NewtonNewmarkSolver` con materiales lineales
reproduce exactamente `NewmarkSolver`** (sin discretización extra ni error
adicional). Esto se valida en los acceptance.

---

## Contrato de implementación

```yaml
name: NewtonNewmarkSolver
kind: solver
status: validated
extends: NewmarkSolver       # variante según Reglas §4

interface:
  type_yaml: NewtonNewmarkSolver
  pipeline_kind: transient    # mismo que NewmarkSolver

parameters_extra_to_NewmarkSolver:
  - { name: convergence,                type: ConvergenceCriterion, required: false, desc: "Política ADR 0007. Default: ConvergenceCriterion() con rtols de fenix.constants." }
  - { name: max_iter,                   type: int,                 required: false, default: 20, desc: "Máximo iteraciones Newton por paso temporal." }
  - { name: freeze_tangent_after_iter,  type: int or None,          required: false, default: null, desc: "Si int, Newton modificado: factoriza fresco las primeras N iter y reusa después (ADR 0003 fase 2)." }

parameters_inherited:
  # Idénticos a NewmarkSolver: t_end, dt, beta, gamma, rayleigh, u0, u0_dot,
  # F_func, linear_algebra, lumping. Ver docs/specs/NewmarkSolver.md.

state_handling:
  commit: "assembler.commit_all_states() tras converger cada paso"
  rollback: "estado trial descartado si el paso no converge → RuntimeError"

acceptance:
  recovery_linear:
    - name: linear_material_matches_NewmarkSolver
      setup: "Mismo modelo y carga que un caso analítico de NewmarkSolver (oscilador 1 GDL, viga), pero con NewtonNewmarkSolver"
      expect: "u_history idéntico (rtol 1e-10) al de NewmarkSolver — convergencia en 1-2 iter por paso"

  nonlinear:
    - name: oscilador_1gdl_elastoplastico_libre
      setup: "1 nodo, Truss2D con Elastoplastic1D, masa concentrada, condición inicial u_0 más allá del yield, vibración libre (F=0)"
      expect: "u(t) refleja decaimiento por disipación plástica; α (deformación plástica acumulada) > 0 al final"

    - name: viga_voladizo_elastoplastica_pulso
      setup: "Frame2DEuler con varios elementos, masa distribuida, pulso de carga concentrada en el extremo libre"
      expect: "respuesta no lineal con redistribución plástica; converge en pocas iter por paso (max ≤6 típicamente)"

  convergence:
    - name: newton_converges_in_few_iter
      setup: "paso plástico activo en un benchmark"
      expect: "número de iteraciones Newton por paso ≤ max_iter; raramente excede 8 con tangente algorítmica consistente"

references:
  - "ADR 0009 §Fase 4 — Transitorio no lineal Newmark + Newton"
  - "Hughes T.J.R. (2000). The Finite Element Method, cap. 7-9 (Newmark, Newton dentro de paso temporal)"
  - "Crisfield M.A. (1991). Non-Linear Finite Element Analysis of Solids and Structures, vol. 2 cap. 24 (dinámica no lineal)"
  - "Spec padre: [docs/specs/NewmarkSolver.md](NewmarkSolver.md)"
```

---

## Implementación

- **Archivo**: [fenix/math/solvers/newmark.py](../../fenix/math/solvers/newmark.py) — junto al padre.
- **Clase**: `NewtonNewmarkSolver(NewmarkSolver)`, registrada con `@SolverRegistry.register`. Hereda `__init__` con kwargs adicionales `convergence`, `max_iter`, `freeze_tangent_after_iter`. Sobrescribe `solve()`.
- **Flujo de `solve()`**:
  1. Ensamblar $\mathbf M$, ensamblar sistema no lineal inicial para obtener $\mathbf K_0$ y $\mathbf F_\text{int,0}$.
  2. Calibrar amortiguamiento Rayleigh con $\mathbf K_0$ (constante en el tiempo).
  3. Reducir por Dirichlet (mismo $\mathbf T$, $\mathbf g$, $\mathbf F_\text{dir}$ que la fase 3).
  4. Aceleración inicial consistente: $\mathbf M\,\ddot{\mathbf u}_0 = \mathbf F_\text{ext}(0) - \mathbf F_\text{int}(\mathbf u_0) - \mathbf C\,\dot{\mathbf u}_0 - \mathbf F_\text{dir}$.
  5. Bucle temporal con Newton interno: predictor → iteración $(\mathbf R, \mathbf J = \mathbf M + \gamma\Delta t\,\mathbf C + \beta\Delta t^2\,\mathbf K_\text{t})$, factoriza/solve, actualiza correctores; converger; commit.
  6. Volcar historial; reportar éxito.
- **Refactor del padre**: la fase 3 (`NewmarkSolver`) sigue intacta. La subclase no extrae métodos del padre por ahora (el flujo no lineal es lo suficientemente distinto como para que la herencia sea estructural, no compartición de pasos internos). Cuando emergió la variante HHT-α (2026-05-18, `NewtonHHTSolver` como subclase de `NewtonNewmarkSolver`) la herencia estructural funcionó tal cual — los métodos `solve()` se sobrescriben íntegramente en cada subclase sin compartir pasos internos. Si en el futuro entran generalized-α u otras variantes y se observa que comparten *este* flujo no lineal específico, entonces se refactorizará al patrón con `_solve_step` virtual.
- **Despacho YAML**: `solver: type: NewtonNewmarkSolver` activado vía `SolverRegistry` y `run_yaml`. Como en `NewmarkSolver`, `PIPELINE_KIND = "transient"` declarativo: si futura regla disparo C reemplaza el `isinstance` en `entry.run_yaml`, no se rompe.
- **Tests**:
  - [tests/test_newmark_nonlinear.py](../../tests/test_newmark_nonlinear.py) nuevo: acceptance.
  - Cobertura de regresión de los tests de `NewmarkSolver` permanece intacta (la subclase no toca el padre).

---

## Diálogo

- **2026-05-14** · Spec creada como **variante** de `NewmarkSolver` aplicando la nueva regla §4 (registrada en el mismo commit que esta spec). Sin ADR nuevo: reusa las decisiones de ADR 0009 fase 4 (residuo, jacobiano, conmutar Newton dentro del paso). Sin entrada propia en el manual estructural — la subclase es invisible al usuario API (excepto por el `type` en YAML y los parámetros extra).
- **2026-05-14** · Decisión: **subclase** sobre `NewmarkSolver` en lugar de bandera `nonlinear=True` en la misma clase. Razones: (i) más explícito en el código y en el YAML (`type` distinto para análisis distinto); (ii) evita condicionales que ensucian el flujo lineal de la fase 3; (iii) deja sitio para variantes adicionales (HHT-α no lineal, generalized-α) por composición. El ADR 0009 fila "Fase 4" decía "`NewmarkSolver` (mismo, con `nonlinear=True`)" — interpretación textualista vs. arquitectural: la regla §4 nueva legitima la elección de subclase, y el patrón conserva el espíritu (mismo análisis, mismas decisiones, mismo entrypoint) sin la implementación literal.
- **2026-05-14** · Amortiguamiento Rayleigh queda fijado con $\mathbf K_0$ (rigidez elástica de referencia) y constante en el tiempo. Alternativa rechazada: Rayleigh con $\mathbf K_\text{t}$ corriente — no es la convención estándar y crea coupling artificial entre la disipación viscosa y la plástica. Documentado en sección "Amortiguamiento Rayleigh — heredado".
