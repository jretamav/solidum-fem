# ADR 0025 — Terminación de los solvers estáticos incrementales y registro por omisión

- **Estado**: aceptado (2026-11-09)
- **Fecha**: 2026-11-09
- **Alcance**:
  - `ArcLengthSolver`, `IndirectDisplacementSolver` (parámetro `control_target`) y `NonlinearSolver`: el estado de terminación y las excepciones tipadas;
  - `DissipationArcLengthSolver` del módulo de usuario `dissipation_control`, por el contrato de su padre;
  - `NewtonCorrector`, `SolveResult.termination` y `solidum.run`;
  - el registro de `solidum/logging.py`.
- **Amplía**: ADR 0005 (logging configurable) con la política de niveles, y ADR 0011 (diagnóstico tipado de divergencia) con la familia del arco.

## Contexto

Los solvers de trayectoria (arco cilíndrico y control indirecto de desplazamiento) son los del post-pico y del ablandamiento, que es donde más se falla. Quien los usaba desde Python no sabía, sin leer el registro, si el análisis terminó ni por qué se detuvo. Medido sobre el código anterior:

1. **No había parada por desplazamiento.** Un análisis de empuje lateral o un ensayo controlado se detienen en un desplazamiento, no en un factor de carga; en una meseta plástica λ casi no cambia y `max_lambda` no sirve. El `IndirectDisplacementSolver` sólo se detenía por `max_lambda` o `max_steps`, y el valor de retorno del `step_callback` se ignoraba.
2. **Un trazado que agotaba `max_steps` regresaba con normalidad.** `solve()` devolvía el último `U` convergido; sólo quedaban un aviso y `reached_max_lambda = False`.
3. **El fracaso del arco era un `RuntimeError` genérico** («Arc-Length fracasó irreparablemente»), fuera de la jerarquía `SolverDivergedError` del ADR 0011. El diagnóstico real, como «Raíces imaginarias», sólo quedaba en mensajes `ERROR` del registro.
4. **Estado inconsistente.** `reached_max_lambda`, `lambda_final` y `steps_done` se fijaban sólo al terminar con normalidad. Tras un fracaso, un solver reutilizado mostraba los de la corrida anterior.
5. **Registro.**
   - Un empuje de 30 pasos escribía 163 líneas `INFO`, 101 de ellas de iteraciones de Newton.
   - Una tangente singular que el arco recuperaba bisecando se escribía como `ERROR`.
   - Con la consola en cp1252 (Git Bash en Windows), los mensajes con letras griegas, subíndices o flechas producían trazas `--- Logging error ---`, aunque el análisis terminara bien.

## Decisión

### 1. Estado de terminación

Tras `solve()`, `NonlinearSolver` y la familia del arco dejan estos atributos:

| Atributo | Significado |
|---|---|
| `termination` | Motivo de la parada (`Termination`, enumeración de cadenas en `solidum.math.solvers.diagnostics`): `"max_lambda"`, `"control_target"`, `"max_steps"` o `"diverged"`; `None` antes de trazar o si `solve` falló antes del primer paso (modelo mal planteado). |
| `completed` | `True` si llegó a donde se pidió (`max_lambda` u objetivo). |
| `lambda_final`, `steps_done`, `n_bisections`, `U_final` | Hasta dónde llegó: factor de carga, desplazamientos y número de pasos del último paso convergido, y bisecciones del trazado. |

- El estado se fija **también cuando `solve` termina con una excepción**.
- `solve` lo **reinicia al empezar**, antes de cualquier comprobación que pueda lanzar.
- `solidum.run` copia `termination` en `SolveResult.termination` y toma `converged` de `completed`, así que llegar al objetivo cuenta como convergido.

**Agotar `max_steps` no lanza excepción**: trazar «hasta donde se pueda» es un uso legítimo. Pero ya no se parece a un éxito: `termination = "max_steps"`, `completed = False`, `SolveResult.converged = False` y un `WARNING`. Se descartó una opción para que lance: el atributo basta para distinguirlo sin interpretar textos, y una opción duplicaría el contrato.

### 2. Excepciones tipadas

Un trazado que fracasa (un paso no converge ni con la longitud mínima) lanza una subclase de `SolverDivergedError` con el diagnóstico de la causa y el último estado convergido en el `Detalle`:

- **`StepConstraintError`** (nueva): la restricción del paso no tiene solución. Son las raíces imaginarias del arco cilíndrico o una magnitud controlada que deja de responder a la carga. El `NewtonCorrector` guarda el motivo con que el problema abortó (`CorrectorResult.abort_reason`).
- **La clasificación del ADR 0011** en el resto de los casos: tangente singular, oscilación, carga sobre la capacidad o modo no clasificado.
- **`SingularTangentError`** si la tangente del predictor es singular en el estado inicial. Antes el bucle repetía el mismo predictor hasta agotar `max_steps` y lo informaba como un trazado incompleto.

Como `SolverDivergedError` hereda de `RuntimeError`, el código que capturaba `RuntimeError` sigue funcionando. Para distinguir el fracaso del arco de un veredicto de la red de seguridad (`IllPosedSystemError`, también `RuntimeError`), se captura `SolverDivergedError`.

### 3. Parada en un valor objetivo de la magnitud controlada

`IndirectDisplacementSolver(..., control_target=c)` detiene el trazado cuando $\mathbf c^\top\mathbf U$ alcanza $c > 0$.

- La restricción del paso fija exactamente $\mathbf c^\top\Delta\mathbf U = \Delta l$ en cada iteración, así que basta recortar el último paso a lo que falta. La llegada es exacta, a precisión de máquina, con cualquier tamaño de paso.
- Un objetivo que cae en medio de un paso reduce ese paso en vez de sobrepasarlo.
- Deja `termination = "control_target"` y `reached_control_target = True`; la magnitud final queda en `control_final`.
- Puede combinarse con `max_lambda`: se detiene en lo que llegue primero.

El mecanismo vive en la base como puntos de extensión que una variante redefine: `_target_reached`, `_limit_step` y `_step_summary`. No hay casos especiales.

### 4. Política de niveles del registro (amplía el ADR 0005)

| Nivel | Contenido |
|---|---|
| `DEBUG` | Detalle por iteración de Newton (línea `Iteración`) y por intento de paso. |
| `INFO` (por omisión) | Progreso: a lo sumo **una línea por paso convergido** (`[PASO n] λ = … \| dl = … \| k iter`), más el inicio y un resumen final. |
| `WARNING` | Lo que el solver recupera por sí solo, con su motivo: bisección por tangente singular, raíces imaginarias o falta de convergencia. También el trazado detenido por `max_steps` y los avisos de siempre. |
| `ERROR` | Sólo lo que detiene el análisis: una línea antes de la excepción tipada. |

El `NewtonCorrector` registra sus propios fracasos en `DEBUG`, porque sólo el solver que lo llama sabe si el paso se recupera.

### 5. Consola con cualquier codificación

El manejador de consola de `solidum` adapta cada mensaje a la codificación del flujo antes de escribirlo. Los caracteres que la codificación no tiene se escriben en texto plano:
- letras griegas por su nombre (`Δλ` → `Deltalambda`);
- subíndices y superíndices (`₁` → `_1`, `²` → `^2`);
- flechas y relaciones (`⇒` → `=>`, `≤` → `<=`);
- letras acentuadas sin acento, si hace falta (ASCII).

Con UTF-8 nada cambia. Se prefirió a reescribir los mensajes en ASCII: los textos siguen siendo legibles en una consola moderna, y la garantía no depende de que cada mensaje futuro respete una regla.

## Consecuencias

- Quien usa los solvers desde Python sabe cómo terminó un análisis por `termination`, sin leer el registro; un fracaso llega como excepción tipada con su causa.
- Un análisis de empuje lateral se detiene en el desplazamiento pedido con llegada exacta.
- El registro por omisión de un análisis de muchos pasos tiene una línea por paso, y un `ERROR` significa siempre que el análisis se detuvo.
- Ningún resultado numérico cambia. Los tests existentes siguen verdes; los nuevos están en `tests/test_trace_termination.py` y `tests/test_solver_logging.py`.
- Código que dependía del texto exacto de los mensajes (por ejemplo «Arc-Length fracasó irreparablemente», o las líneas de iteración en `INFO`) debe usar `termination` y las excepciones tipadas, o el nivel `DEBUG`.

## Alternativas consideradas

- **Detener el trazado con el valor de retorno de `step_callback`.** Da una parada, pero no una llegada exacta; el objetivo como parámetro sí.
- **Excepción al agotar `max_steps`.** Rompe el uso legítimo de «trazar hasta donde se pueda» (§1).
- **Mensajes sólo en ASCII.** Pierde legibilidad y obliga a vigilar cada mensaje nuevo (§5).

## Referencias

- de Borst R. (1987), *Computation of post-bifurcation and post-failure behavior of strain-softening solids*, Computers & Structures 25(2): control indirecto de desplazamiento.
- Crisfield M. A. (1991), *Non-linear Finite Element Analysis of Solids and Structures*, Vol. 1, cap. 9: arco cilíndrico y raíces de la restricción.
- ADR 0005 (logging), ADR 0011 (diagnóstico tipado), ADR 0015 (corrector compartido), ADR 0019 (red de seguridad).
