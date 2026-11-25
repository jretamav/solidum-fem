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

## Ampliación (2026-11-25): análisis detenido, variables internas y puntos de integración

Huecos que los reportes de validación tuvieron que rodear con guiones propios (`informes/modelo3_petersson_yaml/correr.py` guardaba el estado de los puntos en un `.npz`).

6. **Un análisis que se detiene conserva lo calculado.** Si el análisis lanza un error tras al menos un paso convergido (no converge, mecanismo, sistema mal planteado), `YamlModel.solve(write_output=True)` cierra la salida con el último paso convergido —`<nombre>_final.vtu`, `.pvd`, texto y figuras; la historia ya estaba escrita— y deja seguir el error original (`OutputWriter.finish_stopped`; un fallo al cerrar se registra y no lo tapa). `python -m solidum` lo explica y termina con el código 3. Sin ningún paso convergido no se escribe nada.
7. **El VTK por celda lleva todas las variables internas escalares** de los materiales, con su nombre (`damage_c`, `w_max_n`, `cracked`...), promediadas sobre los puntos materiales de la celda (por arreglos en las familias del ADR 0014); `NaN` en las celdas cuyo material no la tiene. `Internal_State` (la principal) se conserva. Un nombre que comparten dos materiales lleva en cada celda el valor del suyo, con el significado de su material.
8. **Ángulos (`STATE_ANGLES`).** Contrato declarativo nuevo de `Material` (por omisión `()`): las variables internas que son la dirección de un eje (la θ de la grieta fija). El VTK por celda no las promedia —la media de las direcciones de puntos distintos, o de un punto sin grieta, no es una dirección—.
9. **Archivo por puntos de integración** (`output.results.integration_points: true`, por omisión `false`): `<nombre>_puntos_NNNN.vtu` en los mismos instantes que el VTK por celda, `<nombre>_puntos_final.vtu` y `<nombre>_puntos.pvd`. Una nube de celdas `vertex` en la posición de referencia de cada punto de los sólidos mecánicos con `Displacements` (`N(ξ)·u_e`), `Element_Id`, `Point`, las seis componentes de esfuerzo, `Von_Mises` y todas las variables internas escalares, también los ángulos. Los valores son los del estado con que el análisis equilibró el paso (el confirmado; en el análisis lineal, el de prueba), los mismos que lee la historia (ADR 0038 §3), sin volver a evaluar el material.
10. **Nombres de campos validados.** Un nombre de `nodal_results` o `element_results` que no es un campo del VTK de ese modelo es un error al crear el escritor, con la lista de los que hay (antes se ignoraba sin aviso); si es un ángulo, el mensaje remite al archivo por puntos.

Alternativas descartadas: promediar los ángulos por celda como ejes (½·atan2 de las sumas de sen 2θ y cos 2θ): los puntos sin grieta guardan θ = 0, que no es una dirección, y haría falta otra declaración de cuándo θ vale; escribir siempre el archivo por puntos (duplica el tamaño de la salida de todos los modelos).
