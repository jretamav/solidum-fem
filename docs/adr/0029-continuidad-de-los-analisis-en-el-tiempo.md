# ADR 0029 — Continuidad de los análisis en el tiempo

- **Estado**: aceptado (2026-11-09)
- **Fecha**: 2026-11-09
- **Alcance**:
  - `TransientTerminationMixin.continue_from(previous)` y sus piezas (`t_start`, `constant_load`, `_external_load`, `_number_of_steps`, `_check_initial_equilibrium`, `_carried_acceleration`), comunes a `NewmarkSolver`, `HHTSolver`, `NewtonNewmarkSolver`, `NewtonHHTSolver`, `CentralDifferenceSolver` y `ThetaMethodSolver`;
  - `run_transient(..., previous=)` y `run_thermal_transient(..., previous=)`;
  - el campo `F_constant` de `TransientResult` y `ThermalTransientResult`.
- **Extiende**: el ADR 0026 (etapas de carga en los solvers estáticos) a los análisis en el tiempo.

## Contexto

Las etapas del ADR 0026 encadenan análisis estáticos. Un programa general de elementos finitos encadena también análisis en el tiempo:
- una carga permanente aplicada de forma estática y, sobre ese estado, una carga dinámica;
- un transitorio que se continúa (por tramos, para cambiar el paso de tiempo, o tras un fallo, desde el último paso convergido);
- un transitorio térmico desde el estado estacionario o desde otro transitorio.

Hasta ahora sólo se podía pasar `u0` a mano y repetir la carga permanente dentro de `F_func`. Eso perdía sus cargas distribuidas por elemento si se sumaba como arreglo simple (ADR 0027), no comprobaba el equilibrio del estado inicial y obligaba a desplazar el tiempo a mano.

## Decisión

Un análisis en el tiempo **continúa** otro con `continue_from(previous)`, o con `previous=` en `run_transient` y `run_thermal_transient`. El estado interno de los elementos es el que guarda el ensamblador (el mismo, como en el ADR 0026).

| Análisis anterior | Estado inicial | Instante inicial | Carga |
|---|---|---|---|
| Estático (`SolveResult`) | sus desplazamientos (o temperaturas), en reposo | 0 | su carga aplicada, **constante**, más `F_func(t)` |
| Transitorio del mismo tipo | su estado final: `u`, `u̇` y `ü` (o `T`) | su instante final | su carga constante, más `F_func(t)` |

1. **Tiempo absoluto.** `t_history` empieza en el instante inicial, `F_func` y `dirichlet_func` se evalúan en tiempo absoluto, y `t_end` es el instante final absoluto. El θ-método conserva `n_steps` como número de pasos.
2. **Carga constante con sus cargas equivalentes.** Se suma a `F_func(t)` como `LoadVector` (ADR 0027): las fuerzas internas de los elementos con carga distribuida siguen siendo correctas. El resultado la guarda en `F_constant`, y un análisis que lo continúe la hereda.
3. **Equilibrio del estado inicial.** Al continuar un estático se comprueba `Tᵀ·F_c − F_int,red(u₀) ≈ 0`. Si no se cumple, se avisa: la diferencia actuaría como una fuerza súbita.
4. **La aceleración se hereda.** Al continuar un transitorio, la aceleración inicial es la final del anterior, no la del equilibrio instantáneo. Con HHT-α no coinciden, porque la ecuación del método promedia con α. Así **parar y continuar es igual que no parar**: medido en los cinco integradores mecánicos, con plasticidad en las fibras, y en el θ-método.
5. **Excluyentes.** Las condiciones iniciales del constructor (`u0`, `u0_dot`, `T_initial`) y `previous` se excluyen. Un resultado parcial (análisis detenido) se admite con un aviso.

## Consecuencias

- La gravedad estática seguida de una excitación dinámica es un análisis de dos líneas, con las fuerzas internas correctas desde el primer instante.
- Un transitorio que se detiene puede continuarse desde su último paso convergido con otro paso de tiempo.
- **Pendiente (decisión física del usuario, fase 3)**: el amortiguamiento de Rayleigh proporcional a la rigidez usa la tangente en el instante inicial de cada análisis (`K` en `u₀`). Al continuar un transitorio con `β ≠ 0` y material no lineal, la matriz `C` cambia respecto de la de una sola corrida.
- Fuera de este ADR: un estático que continúe un transitorio (qué carga se mantiene y qué pasa con la velocidad son decisiones físicas) y una perturbación lineal sobre el estado de una etapa (fase 3).
