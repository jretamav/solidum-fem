# ADR 0039 — Figuras de los resultados desde el YAML

- **Estado**: aceptado (2026-11-11)
- **Fecha**: 2026-11-11
- **Alcance**:
  - `solidum/post/` (nuevo): `plots.py` (`make_plots`, `validate_plots_block`, `ResultFiles`) y `__main__.py` (`python -m solidum.post modelo.yaml`);
  - `solidum/utils/output_writer.py`: valida `output.plots`, dibuja las figuras al terminar y escribe siempre el `.pvd` cuando hay VTK por instante;
  - `pyproject.toml`: grupo opcional `plot` (matplotlib), también en `dev`.

## Contexto

Con los ADR 0035 a 0038, un estudio de validación se declara entero en el YAML: la geometría o la malla, el análisis y la salida (curva de respuesta, energías, daño, VTK en los instantes elegidos). Faltaban las **figuras**: la malla, la deformada coloreada por un campo y las curvas. Se hacían con scripts propios de la IA fuera de Solidum, lo que el usuario no quiere (2026-11-11): Solidum debe usarse como lo usaría un usuario final, y el código externo sólo debe cubrir lo que Solidum no hace. Los programas generales traen su propio post-proceso gráfico (Abaqus/Viewer, el visor de ANSYS, las gráficas de FEAP).

## Decisión

1. **Las figuras se dibujan con los archivos que escribe Solidum**: el `.pvd` como índice de factor → VTK por instante, `<nombre>_final.vtu` y la historia CSV. No se usa el estado en memoria, así que las mismas figuras salen al terminar `python -m solidum modelo.yaml` o, sin volver a calcular, con `python -m solidum.post modelo.yaml`. El `.pvd` se escribe siempre que haya al menos un VTK por instante (antes, sólo con dos o más).
2. **`output.plots`**, una lista de figuras:
   - `mesh`: la malla sin deformar, con los apoyos y los desplazamientos impuestos (un grado de libertad restringido que se movió) y las cargas nodales;
   - `deformed`: la deformada en los instantes escritos más cercanos a los factores pedidos (`final`, el último), coloreada por un campo del VTK (por celda, por nodo, o `displacement`), un panel por instante con la misma escala de color; escala de la deformada automática (el mayor desplazamiento dibujado mide el 5 % del tamaño) o dada; ventana `zoom`; rótulo de cada panel con una columna de la historia;
   - `curve`: columnas de la historia, también de otros modelos (`compare`), para los estudios de malla.
3. **Colores** de la guía de visualización del proyecto, que es la paleta validada para daltonismo: curvas en el orden categórico fijo, con leyenda desde dos series; campos positivos (daño) con una rampa azul de un solo tono cuyo cero se funde con el fondo; campos con signo con la divergente azul–gris–rojo centrada en cero.
4. **Un fallo al dibujar no invalida el análisis**: al terminar `python -m solidum`, el error se registra y los resultados quedan escritos; las figuras se regeneran con `python -m solidum.post` tras corregir el YAML.
5. **Alcance**: mallas 2D (sólidos planos y elementos de línea en el plano), en PNG, PDF o SVG. Las 3D, no todavía (error claro).

## Alternativas descartadas

- **Dibujar desde el estado en memoria al terminar.** Obligaría a resolver de nuevo para cambiar una figura.
- **Delegar en ParaView.** El usuario no maneja ParaView; las figuras de un reporte deben salir del mismo flujo que el cálculo.
- **Dibujar en cada paso.** Las figuras se piden en unos pocos instantes; los VTK ya se escriben sólo en esos instantes (ADR 0036).

## Consecuencias

- Un reporte de validación sale de Solidum: geometría, malla, deformadas, curvas y comparación entre mallas. Fuera de Solidum queda redactar e interpretar.
- matplotlib pasa a ser dependencia opcional (`plot`) y de desarrollo; sin ella, `output.plots` da un error claro y el análisis no se ve afectado.
- Pendiente: figuras 3D, campos nodales suavizados (hoy cada celda lleva un valor) y figuras por etapa.
