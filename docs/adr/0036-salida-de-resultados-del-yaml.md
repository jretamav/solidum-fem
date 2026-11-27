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

## Ampliación (2026-11-27): historia de un análisis por etapas

Hueco del reporte de validación del contacto (`informes/reportes/reporte_contacto_2026-11-27.pdf`): con etapas, el escritor se creaba como si el análisis no fuera estático y sin las cargas de cada etapa, así que rechazaba las reacciones y las energías; la fuerza del despegue salió de la tracción de la interfaz y la comprobación de la disipación por fricción necesitó una corrida aparte en una sola etapa. Además, las fuerzas internas de un marco con carga distribuida en una etapa se calculaban sin descontar su carga equivalente, sin aviso.

11. **Columna `stage`.** Con etapas, la historia lleva el número de la etapa (desde 1) entre `step` y `factor`; `factor` es el λ de la etapa. Es numérica para que la lean las figuras (`curve`) y los lectores de CSV; los nombres están en el YAML y en `StagedResult.names`. `stage` pasa a ser un nombre de columna reservado, como `step` y `factor`.
12. **La carga de cada etapa.** `OutputWriter.begin_stage(número, F_constante, F_referencia)`, que el despacho de etapas llama antes de cada etapa con pasos: en una estática, la carga del paso es `F_c + λ·F_ref` (ADR 0026), con la que se calculan la reacción, la fuerza externa del VTK y de la salida de texto y, con las cargas equivalentes de cada parte (ADR 0027), las fuerzas internas de un marco con carga distribuida; en una etapa en el tiempo, `None, None`. Las de perturbación lineal, modal, armónica y espectral no tienen pasos y no la llaman.
13. **Reacciones y energías con etapas**, si ninguna etapa es en el tiempo (el escritor lo decide al crearse, con la clase del solver de cada etapa). El trabajo externo se acumula de una etapa a otra con el mismo trapecio de `F_int·ΔU`, que no depende de cómo se reparta la carga (ADR 0038). La reacción de un grado de libertad se admite si lo prescribe el nivel raíz o alguna etapa (ADR 0026 §7); en cada paso se suma sobre los nodos del selector que lo tienen prescrito en ese momento, y en una etapa en la que está libre en todos la celda queda vacía: `F_int − F` sería sólo el residuo del equilibrio, no la reacción de un apoyo. Una celda vacía es «no aplica» (las figuras la leen como `NaN`, un hueco en la curva; `ResultFiles.history` las convierte, corregido tras la validación: antes ninguna curva se dibujaba si la historia tenía una celda vacía).
14. **Fuerzas internas de marco con etapas en el tiempo.** El escritor no conoce las cargas de una etapa en el tiempo, y la carga equivalente constante de una etapa estática sigue actuando en ella: si alguna etapa es en el tiempo y alguna tiene carga distribuida, una columna de fuerza interna de marco es un error al crear la salida. Antes se escribía sin descontar la carga equivalente, sin aviso.

Validado como usuario (`informes/reportes/reporte_contacto_etapas_2026-11-27.pdf`, veredicto «funciona»): despegue con carga y descarga (reacción con error 1e-10 de la fuerza máxima, W exacto, |D| ≤ 5e-16·W), fricción 2D de ida y vuelta (fuerza de fricción 1e-15 de μpA, D exacto en todos los pasos) y fricción 3D con trayectoria en L (retorno radial discreto a 4e-10·μpA; frente a la solución continua, error de primer orden 1.2/0.30/0.08 % con 20/80/320 pasos).

Verificado con valores exactos (`tests/test_stage_boundary_conditions.py`, 23 pruebas): fricción con cierre, deslizamiento y retorno (fuerza μpA y pA, `D = μpA·Σ deslizamiento`), barra con carga, desplazamiento impuesto y vuelta (`R = k·u − P`, `W = E`, `D = 0`), y voladizo con peso propio y carga de punta en dos etapas (`M = −qL²/2 − λPL`); 20/20 mutaciones (`informes/prototipos_2026-11-27/etapas_apoyos/mutar_etapas.py`). Revisión con `/code-review` (2026-11-27): atendidos la carga del paso recalculada por nodo en la salida de texto, la validación al leer, los apoyos restaurados en cada `solve()`, la reacción vacía donde no hay apoyo y el rechazo de las fuerzas de marco en etapas en el tiempo.

Sigue pendiente (anterior a esta ampliación): con etapas, el `factor` del `.pvd` vuelve a empezar en cada etapa (instantes repetidos para ParaView y las figuras), `frequency: {at: [...]}` y `deformed: {at: [...]}` sólo alcanzan la primera; en las etapas en el tiempo la fuerza externa del VTK y de la salida de texto es cero; la figura de malla toma los apoyos del VTK final, que con etapas son los de la última (un grado de libertad que volvió a 0 aparece como apoyo fijo); una figura que falla impide dibujar las siguientes (validador, 2026-11-27).
