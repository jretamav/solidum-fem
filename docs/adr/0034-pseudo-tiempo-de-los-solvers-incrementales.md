# ADR 0034 — Pseudo-tiempo de los solvers incrementales

- **Estado**: aceptado (2026-11-11)
- **Fecha**: 2026-11-11
- **Alcance**:
  - `solidum/math/assembly.py`: `Assembler.prepare_all_steps(U, pseudo_time_increment=None)`, `Assembler.set_pseudo_time_increment(Δτ)`, `Assembler.pseudo_time_increment` y el cierre del intento en `commit_all_states`;
  - `solidum/core/material.py`: `Material.uses_pseudo_time` y `Material.set_pseudo_time_ratio(r)`;
  - los solvers que avanzan por pasos: `NonlinearSolver`, `ArcLengthSolver` y sus variantes (`IndirectDisplacementSolver`), `NewtonNewmarkSolver`, `NewtonHHTSolver`, `CentralDifferenceSolver` y el del módulo de usuario `dissipation_control`;
  - primer consumidor: la integración IMPL-EX de `IsotropicDamageCrackBand2D/3D` ([spec](../specs/IsotropicDamageCrackBand2D.md), §8).
- **Extiende**: el ADR 0020, P5 (gancho de inicio de paso), y el ADR 0014 (camino por lotes).

## Contexto

La integración IMPL-EX (Oliver, Huespe y Cante 2008), aprobada por el usuario para el programa principal el 2026-11-11, extrapola la variable histórica de cada punto con los dos últimos pasos convergidos: $\tilde\kappa_{n+1} = \kappa_n + r\,(\kappa_n - \kappa_{n-1})$, con $r = \Delta\tau_{n+1}/\Delta\tau_n$ la razón entre el incremento del parámetro de control del paso en curso y el del último paso confirmado. Ese dato es del solver, no del punto, y hasta ahora no llegaba al material: `compute_state(strain, state_vars)` sólo recibe la deformación y la historia, y el camino por lotes calcula `batch_params()` una vez, al construir la familia.

El incremento sólo se conoce en momentos distintos según el solver: el control de carga lo sabe antes de empezar el paso; la longitud de arco, después del predictor (que necesita la tangente en el estado convergido). Y un paso puede repetirse con otro incremento tras una bisección.

## Decisión

1. **El ensamblador lleva el intento de paso.** `prepare_all_steps(U, pseudo_time_increment=None)` abre un intento (además del gancho de los elementos del ADR 0020). El solver declara el incremento de su parámetro de control, $\Delta\tau_{n+1} > 0$, en esa llamada si ya lo conoce, o después con `set_pseudo_time_increment(Δτ)`, siempre antes del corrector. `commit_all_states` cierra el intento: el incremento del paso confirmado pasa a ser la referencia $\Delta\tau_n$ del siguiente.
2. **El material recibe la razón, no el incremento.** El ensamblador calcula $r = \Delta\tau_{n+1}/\Delta\tau_n$ ($r = 1$ en el primer paso) y la entrega con `set_pseudo_time_ratio(r)` a los materiales de los elementos que declaran `uses_pseudo_time`. Fuera de un paso (antes del primero, entre el commit y el siguiente intento, y en los análisis lineales y el post-proceso) entrega `None`, y el material congela su daño en el valor confirmado: reevaluar el estado confirmado reproduce el esfuerzo del paso, y el commit final de `solidum.run` no cambia nada. La razón vive en la instancia del material: no es historia, y `compute_state` sigue siendo puro respecto de la historia.
3. **Camino por lotes.** Al cambiar la razón, el ensamblador refresca en sitio los `mat_params` de las familias de esos materiales (el kernel compilado los lee de ahí); la razón viaja como un parámetro más.
4. **Qué es el pseudo-tiempo de cada solver.**
   - `NonlinearSolver`: el incremento del factor de carga $\Delta\lambda$.
   - `ArcLengthSolver` y sus variantes: la longitud de arco del paso (en el control indirecto, el incremento de la magnitud controlada), tras el predictor.
   - `NewtonNewmarkSolver`, `NewtonHHTSolver`: el paso de tiempo del intento, incluidos los subpasos.
   - `CentralDifferenceSolver`: $\Delta t$ en cada instante.
   - `DissipationArcLengthSolver` (módulo `dissipation_control`): la magnitud del modo activo, $dl$ en el arco cilíndrico y $\tau$ en el de disipación.
5. **Bisección.** Repetir el intento con otro incremento sólo cambia $r$ respecto de la misma referencia confirmada: el resultado de un intento es el mismo la primera vez y tras una bisección.
6. **Rechazo.** Si un solver abre un intento con `prepare_all_steps` y lo confirma sin haber declarado su incremento, y el modelo tiene materiales con `uses_pseudo_time`, `commit_all_states` lanza `ValueError` nombrando los materiales y los solvers que sí lo declaran. Ocurre al confirmar el primer paso, antes de que el daño evolucione con una razón que no existe.

## Alternativas descartadas

- **La razón en el estado de cada punto.** Habría que escribirla en todas las filas de estado en cada intento, y ensuciaría el esquema de historia (`STATE_SCHEMA`) con un dato que no es historia.
- **Un argumento más en `compute_state` y en la firma del kernel por lotes.** Rompería el contrato de todos los materiales y la firma compilada `MAT_SIG` para servir a uno.
- **Una variable global del paquete.** Frágil con dos modelos a la vez o con análisis encadenados.

## Consecuencias

- Un solver nuevo que avance por pasos debe declarar su pseudo-tiempo para admitir materiales IMPL-EX; si no lo hace, el modelo lo rechaza en el primer commit con un mensaje claro. La plantilla de solver de la skill `solidum-new` lo recuerda.
- La referencia $\Delta\tau_n$ se conserva entre análisis sobre el mismo ensamblador (continuidad, ADR 0029). Si el análisis siguiente usa otro parámetro de control, la razón de su primer paso compara magnitudes distintas; lo mismo pasa en el solver de disipación en el primer paso tras un cambio de modo. El efecto es un paso con la extrapolación mal escalada.
- En la llegada exacta a `max_lambda` del arco, el paso final usa la razón de su intento.
- `Material.compute_state` documenta la excepción: la instancia guarda, además de los parámetros, la razón del intento en curso.
