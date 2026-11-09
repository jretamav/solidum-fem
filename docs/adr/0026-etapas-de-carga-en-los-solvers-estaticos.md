# ADR 0026 — Etapas de carga en los solvers estáticos

- **Estado**: aceptado (2026-11-09)
- **Fecha**: 2026-11-09
- **Alcance**:
  - `NonlinearSolver.solve`, `ArcLengthSolver.solve` (y `IndirectDisplacementSolver`, que lo hereda): parámetros `u0` y `constant_load`;
  - `Assembler.reduce` (`g_start`), `Assembler.commit_all_states` (`Domain.committed_history`), `Assembler.equivalent_loads_version`;
  - `solidum/math/solvers/stages.py`;
  - `solidum.run(..., previous=...)` y `SolveResult` (`F_constant`, `F_variable`, `load_factor`, `equivalent_loads`).
- **No incluye**: `DissipationArcLengthSolver` (módulo de usuario), que rechaza las etapas con un error explicativo (§6).

## Contexto

Un análisis de empuje lateral correcto aplica primero la gravedad, que se queda, y después empuja con una carga lateral que crece. Sin eso no hay fuerzas axiales en las columnas, ni efecto P-Δ con la carga vertical, ni patrones de carga que se recalculan en cada paso desde el estado vigente. La necesidad es general: cualquier análisis con cargas permanentes, como el peso propio con una sobrecarga creciente, las secuencias constructivas o el pretensado.

Medido sobre el código anterior:
- Los solvers estáticos arrancaban siempre en `U = 0`, con λ = 0, y su residuo no tenía término constante: era `λ·F_ref − F_int`.
- Una segunda corrida tras aplicar la gravedad la perdía: la reacción vertical caía a 6e-6 N cuando debía ser 5.1e7 N.
- Tras una primera etapa que plastificaba, el estado interno de los elementos se conservaba pero `U` volvía a cero. Con `U = 0` la fuerza interna era 7 veces la de primera fluencia, en un modelo que debía estar en reposo.

Con carga proporcional, el elemento corrotacional ya resolvía el efecto P-Δ con orden 2. Sólo faltaba escalonar las cargas.

## Decisión

### 1. Una etapa: estado de partida y carga constante

`solve(F_ref, step_callback=None, *, u0=None, constant_load=None)` en `NonlinearSolver` y en la familia del arco:

$$\mathbf R(\mathbf U, \lambda) = \mathbf F_c + \lambda\,\mathbf F_\text{ref} - \mathbf F_\text{int}(\mathbf U),\qquad \lambda = 0 \text{ al empezar}.$$

- **`u0`**: los desplazamientos de partida, los de la etapa anterior. El estado interno de los elementos lo guarda el ensamblador (las familias por lotes, ADR 0014) y los elementos: se usa el mismo `Assembler`.
- **`constant_load`**: la carga ya aplicada, que actúa completa durante toda la etapa. Se genera con las herramientas de siempre: cargas puntuales, `assemble_self_weight`, `assemble_body_load` y cargas distribuidas.
- La semántica de λ, `max_lambda`, `control_target` (ADR 0025) y la terminación son las de una corrida, medidas en la etapa. `control_target` es un valor absoluto de la magnitud controlada, así que continúa desde donde quedó.
- El criterio de convergencia se calibra con la carga de la etapa completa, $\mathbf F_c + \mathbf F_\text{ref}$.

Sin `u0` ni `constant_load` todo es exactamente como antes.

### 2. Valores impuestos

Los desplazamientos prescritos y el término independiente de las restricciones multipunto crecen desde lo que valían en `u0` hasta el valor del modelo:

$$\mathbf g(\lambda) = \mathbf g_0 + \lambda\,(\mathbf g - \mathbf g_0),\qquad \mathbf g_0 = \mathbf U_0 - \mathbf T\,\mathbf U_0[\text{libres}].$$

- Un apoyo que no cambia entre etapas queda constante; uno que cambia crece con la etapa.
- Sin estado inicial, $\mathbf g_0 = \mathbf 0$ y se recupera el $\lambda\,\mathbf g$ de siempre (ADR 0004).
- `Assembler.reduce` recibe `g_start`. En la dirección tangente $\partial\mathbf U/\partial\lambda$, la parte impuesta es $\mathbf g - \mathbf g_0$.

### 3. Equilibrio al empezar

Cada etapa comprueba que empieza en equilibrio: $\lVert\mathbf T^\top(\mathbf F_c - \mathbf F_\text{int}(\mathbf U_0))\rVert$ con el criterio de fuerza del solver (ADR 0007), calibrado en una copia. El residuo relativo queda en `solver.initial_residual`, y si no hay equilibrio se avisa con un `WARNING` que explica cómo continuar.

La misma comprobación delata el caso medido en el contexto. `commit_all_states` marca el dominio con `committed_history`, y una corrida que empieza en `U = 0` sobre un modelo con historia avisa si ese estado no está en equilibrio ahí. Un modelo elástico, con `F_int(0) = 0`, no avisa.

### 4. `solidum.run(..., previous=resultado)`

La forma habitual de continuar es pasar el `SolveResult` de la etapa anterior:

```python
asm = Assembler(dom)
g = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm), F_applied=F_gravedad)
e = solidum.run(dom, assembler=asm, F_applied=F_lateral, previous=g,
                solver=IndirectDisplacementSolver(asm, control=[(tope, "ux", 1.0)],
                                                  max_lambda=1e9, control_target=0.09))
```

- `u0` es `previous.U` y `constant_load` es `previous.F_applied`: lo aplicado sigue aplicado, como en los pasos de Abaqus y ANSYS o con `loadConst` de OpenSees.
- El resultado separa la carga constante (`F_constant`) de la que creció (`F_variable = F_applied − F_constant`, `load_factor` = λ final). Las reacciones se calculan con la total.
- Un solver sin etapas (`LinearSolver`) da un error explicativo.

### 5. Cargas distribuidas en las fuerzas internas

Las fuerzas internas de extremo restan las cargas equivalentes de la carga distribuida (`F_int − f_eq`). En una etapa:

$$\mathbf f_\text{eq} = \mathbf f_\text{eq}^{\,c} + \lambda\,\mathbf f_\text{eq}^{\,\text{ref}}.$$

- La parte constante es la de la etapa anterior, guardada en `SolveResult.equivalent_loads`.
- La parte de la carga que crece es la del registro del ensamblador **sólo si se ensambló después** de la etapa anterior. Lo dice `Assembler.equivalent_loads_version`, que aumenta con cada carga de cuerpo ensamblada.
- Así, el peso propio de la etapa 1 entra con factor 1 y una fuerza de cuerpo lateral de la etapa 2, con λ.

### 6. Fuera de este ADR

- **`DissipationArcLengthSolver`** (módulo de usuario) rechaza `u0` y `constant_load` con `NotImplementedError`. Su restricción usa el trabajo de la carga de referencia, y con una carga constante la energía disipada del paso cambia. Eso es física por decidir, no un detalle de implementación.
- **Etapas desde el YAML**: no hay sintaxis todavía (deuda #36).
- **Dos cargas de cuerpo en una misma etapa.** El registro del ensamblador sólo guarda las cargas equivalentes de la última ensamblada, así que `element_forces` descuenta sólo esa. Es un defecto previo, también en análisis de una etapa (deuda #35). Lo resuelve un modelo de cargas que lleve sus cargas equivalentes (mejora 4 de la lista del 2026-09-22).

## Verificación (criterios de aceptación)

- **P-Δ por etapas.** Voladizo corrotacional de 20 fibras con P = 0.3·P_cr y H = 1e4 N:
  - la etapa de gravedad seguida de la lateral (con Newton o con el arco) da la flecha de una sola etapa proporcional a 1e-13;
  - contra la viga-columna con acortamiento axial (δ = 2.045488e-4 m), el error es 3.42e-4 con 16 elementos y el orden observado 2.00.
- **Interacción N–M.** Acero perfectamente plástico con n = 0.3:
  - el exceso sobre M_p(1 − n²)/L es +11.7 %, +5.9 % y +3.0 % con 4, 8 y 16 elementos (+1.5 % con 32);
  - es de primer orden y su extrapolación de Richardson da la carga límite (−0.04 %);
  - la reacción vertical suma P durante todo el empuje.
- **Coherencia elástica.** Las etapas reproducen la superposición a 1e-10.
- **Continuación.** Continuar el empuje de 0.03 a 0.06 m da la misma curva que empujar de una vez hasta 0.06 m, a 1e-9 en la fuerza.
- **Equilibrio inicial.** Tras una etapa plástica, el residuo inicial es de 1e-15. Empezar en `U = 0` con el estado plástico avisa.
- **Compatibilidad.** Con apoyos con desplazamiento impuesto, el oráculo del nodo fundido de las restricciones multipunto, la red de seguridad (ADR 0019) y las cargas distribuidas constantes y crecientes en las fuerzas internas.

## Alternativas consideradas

- **La carga de la etapa nueva como incremento sobre la anterior**, con $\mathbf F(\lambda) = \mathbf F_\text{prev} + \lambda(\mathbf F_\text{nueva} - \mathbf F_\text{prev})$: generaliza, pero oculta qué carga es constante y obliga a reescribir toda la carga en cada etapa.
- **Que el solver recuerde su último estado y continúe solo**: acopla la continuación a reutilizar el mismo objeto y no sirve para cambiar de solver entre etapas (Newton para la gravedad, control indirecto para el empuje).
- **Un objeto de carga con sus cargas equivalentes**: es la solución general del §5 y de la deuda #35, pero cambia la API de todas las cargas. Queda para el modelo de cargas.

## Referencias

- de Borst R., Crisfield M. A., Remmers J. J. C., Verhoosel C. V. (2012), *Nonlinear Finite Element Analysis of Solids and Structures*, 2.ª ed., Wiley, cap. 4: procedimientos incrementales-iterativos.
- OpenSees, `loadConst`; Abaqus, pasos de análisis (*step*) con cargas propagadas; ANSYS, pasos de carga.
- ADR 0004 (Dirichlet y restricciones multipunto), ADR 0007 (convergencia), ADR 0014 (estado por familias), ADR 0019 (red de seguridad), ADR 0025 (terminación y control_target).
