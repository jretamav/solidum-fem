# ADR 0036 — Salida de resultados del YAML: historia de la respuesta y frecuencia de VTK

- **Estado**: aceptado (2026-11-11)
- **Fecha**: 2026-11-11
- **Alcance**:
  - `solidum/utils/output_writer.py` (nuevo): `OutputWriter`, `validate_output_block`;
  - `solidum/entry.py`: `run_yaml(path, write_output=False)`, `YamlModel.solve(write_output=False)`;
  - `solidum/__main__.py` (nuevo): `python -m solidum modelo.yaml`;
  - `examples/ejecutar_yaml.py`: envoltura de `python -m solidum`;
  - `solidum/utils/yaml_parser.py`: validación de `output.results.frequency` y `output.history`.

## Contexto

Para validar los modelos de concreto por el camino normal de un usuario (geometría `.geo` + YAML, ADR 0035) faltaban dos salidas estándar en los programas de elementos finitos (la *history output* de Abaqus, las gráficas de FEAP):

- la **curva de respuesta** paso a paso (p. ej., carga–flecha): `SolveResult` sólo guarda el estado final;
- los **resultados VTK en instantes elegidos**: la salida era «todos los pasos» o «sólo el último»; en la viga de Petersson, todos son 4000 archivos.

Además, la escritura de archivos vivía en el script `examples/ejecutar_yaml.py`, con un despacho propio reducido (sólo las cargas nodales: sin etapas, peso propio ni cargas térmicas), distinto del de `solidum.run_yaml`.

## Decisión

1. **El escritor vive en el paquete.** `OutputWriter` se engancha al `step_callback` del análisis (el mismo de `run`, `run_transient`, `run_thermal_transient` y las etapas) y escribe junto al YAML lo que declara su bloque `output`: el modelo sin deformar, los VTK por paso con su `.pvd`, la historia, la salida de texto estilo FEAP (trasladada del script) y `<nombre>_final.vtu`.
2. **Escribir es opcional y explícito.** `run_yaml(path, write_output=True)`, `YamlModel.solve(write_output=True)` y `python -m solidum modelo.yaml`. Por omisión no se escriben archivos: las pruebas que resuelven los ejemplos no ensucian el repositorio. `examples/ejecutar_yaml.py` pasa a ser una envoltura de `python -m solidum` y usa el despacho completo.
3. **Frecuencia de VTK**: `last_step` (por omisión), `all_steps`, `{every: N}` (cada N pasos convergidos) o `{at: [f1, f2, ...]}` (el primer paso convergido con factor de carga, o tiempo, ≥ cada valor).
4. **Historia**: `output.history.columns`, cada columna con un nombre y una magnitud: `displacement` (media, mínimo o máximo de un grado de libertad en un nodo, una lista de nodos, un grupo físico o los nodos en una coordenada) o `reaction` (suma de las reacciones de un grado de libertad prescrito), con un factor `scale` opcional (p. ej. −1 para el signo de la carga). Se escribe un CSV fila a fila (`step`, `factor` y las columnas), de modo que un análisis interrumpido conserva lo calculado. La reacción de un paso es `F_int(U) − factor·F` en los grados de libertad prescritos (ADR 0004 §3); se admite sólo en análisis estáticos sin etapas, y el selector debe tener ese grado de libertad prescrito en todos sus nodos.
5. **Validación** de forma en la lectura del YAML (frecuencia, columnas, selectores), con todos los errores juntos; la resolución de nodos y grupos, al crear el escritor.

## Alternativas descartadas

- **Guardar la historia en `SolveResult`.** Obligaría a decidir qué cantidades guarda todo análisis; la historia depende de lo que el usuario quiere observar.
- **Escribir siempre que haya `output`.** Las pruebas que resuelven todos los ejemplos escribirían archivos en `examples/`.
- **Mantener la escritura en el script.** Duplicaba el despacho y dejaba fuera etapas, peso propio y cargas térmicas.

## Consecuencias

- Un estudio de validación es un `.geo`, un YAML y `python -m solidum`, con la curva en CSV y las deformadas en VTK.
- La reacción con restricciones multipunto de apoyo (rodillo oblicuo) no se registra en la historia: el selector exige el grado de libertad prescrito.
- Los archivos por paso pasan a llamarse `<nombre>_NNNN.vtu` (cuatro cifras).
