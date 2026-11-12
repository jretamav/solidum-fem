# ADR 0037 — Mallas estructuradas por bloques y calidad de la malla

- **Estado**: aceptado (2026-11-11)
- **Fecha**: 2026-11-11
- **Alcance**:
  - `solidum/mesh/structured.py`: `axis_lines`, y en `structured_mesh` los argumentos `lines`, `exclude`, `diagonal` y `merge_tol`;
  - `solidum/mesh/structured_blocks.py` (nuevo): `mesh: {structured: ...}` del YAML;
  - `solidum/mesh/quality.py` (nuevo): `mesh_quality`, `log_mesh_quality`;
  - `solidum/utils/expressions.py`: `compile_predicate` (condiciones);
  - `solidum/utils/yaml_parser.py`: validación y construcción de `mesh.structured`; resumen de calidad al leer una malla.

## Contexto

El generador de mallas estructuradas (ADR 0022) malla un único bloque uniforme, con los triángulos siempre por la misma diagonal; se usaba sólo desde Python, en pruebas y estudios de convergencia. Los problemas de validación de la familia de concreto (viga de Petersson, viga con entalla, panel en L) son geometrías de bloques rectos con una entalla o una ranura, que piden afinar la malla hacia la zona de fisuración y alinear una línea de la rejilla con la punta de la entalla y con el punto de carga. Con Gmsh (ADR 0035) se pueden mallar, pero una malla estructurada da control exacto de las líneas, una banda recta por construcción y la misma malla en primer y segundo orden; los programas generales la ofrecen (bloques de FEAP, mallado mapeado de ANSYS, líneas transfinitas con progresión de Gmsh).

Tampoco había ninguna medida de la calidad de los elementos: una malla distorsionada se notaba sólo por sus resultados o, si llegaba a invertir un elemento, por el error del ensamblaje.

## Decisión

1. **Líneas de la rejilla por tramos.** `axis_lines(segments, origin)` da las líneas de un eje como tramos consecutivos `(length, n, ratio)`: `n` celdas en progresión geométrica cuyo último tamaño es `ratio` veces el primero (la razón de sesgo de Abaqus). `structured_mesh(..., lines=...)` malla sobre esas líneas: la referencia sigue siendo `[0, 1]^d` por índice de celda y la caja es lineal por tramos en ella, de modo que los nodos intermedios de un elemento cuadrático quedan en el centro de sus bordes. Sin `lines`, el comportamiento anterior, sin cambios.
2. **Celdas excluidas.** `exclude(x)` se evalúa en el centro físico de cada celda; si es verdadero, la celda no se malla y no se crea ninguno de sus nodos. Las entallas y ranuras se alinean con la rejilla mediante los tramos.
3. **Partición en triángulos**: `diagonal` = `right` (la de siempre), `left`, `alternate` (en tablero de ajedrez) o `crossed` (cuatro triángulos alrededor del centro de la celda). Sólo en celdas `tri`.
4. **Varios bloques.** Con `merge_tol`, un nodo nuevo a esa distancia de un nodo que ya estaba en el dominio se funde con él. Dos bloques forman un dominio conforme si comparten las líneas de la rejilla en el lado común (el usuario da los mismos tramos).
5. **YAML**: `mesh: {structured: ...}` con un bloque (las claves directamente) o `blocks: [...]`. Cada bloque declara sus ejes, `origin`, `mapping` (expresiones de `x, y, z`), `exclude`, `diagonal`, y su receta de elemento con las mismas claves que un grupo físico de Gmsh (`element`, `material`, `thickness`, `quadrature`, parámetros del elemento), con `mesh_element`, `mesh_material`... por omisión. Al nivel de `structured`, `groups` y `merge_tol` (por omisión, 10⁻⁶ veces la celda más pequeña).
6. **Grupos**: los mismos objetos que da una malla de Gmsh (`physical_groups`, nodos; `physical_group_cells`, bordes o caras de contorno), para que `boundary_conditions_by_group`, `point_loads_by_group`, `tractions_by_group` y la historia de la salida (ADR 0036) funcionen igual. Se crean los lados de cada bloque (`xmin`... con un solo bloque sin nombre; `<nombre>.xmin`... en otro caso) y los de `groups`, definidos por una condición sobre la posición: sus nodos, y los bordes o caras del contorno libre con todos sus nodos en el grupo.
7. **Condiciones** en el compilador seguro de expresiones (ADR 0023): `compile_predicate` añade `< <= > >=`, `and`, `or` y `not`, con la misma lista blanca; la igualdad entre reales no se admite y en su lugar está `near(a, b[, tol])`, con tolerancia por omisión de 10⁻⁶ veces la diagonal de la caja que envuelve la malla (del bloque, en `exclude`). Las expresiones de carga siguen sin admitir comparaciones.
8. **Calidad de la malla.** `mesh_quality(domain)` mide cada sólido 2D y 3D: la razón de jacobianos (`max det J / min det J` en los nodos del elemento, con su propia cinemática compilada, la de `BATCH_KINEMATICS`; `inf` si `det J ≤ 0` en algún nodo), la razón de aspecto (arista mayor / menor entre vértices) y, en 2D, los ángulos de los vértices. Al leer una malla de Gmsh o generar una por bloques desde el YAML, el registro resume la malla (INFO) y avisa (WARNING) de los elementos que pasan los límites de aviso por omisión de la comprobación de forma de ANSYS: razón de jacobianos 30, aspecto 20, ángulo máximo 165° en triángulos y 155° en cuadriláteros. Es un diagnóstico y no detiene el análisis; un jacobiano no positivo lo sigue deteniendo el ensamblaje.

## Alternativas descartadas

- **Gradación como función de mapeo** (`mapping` que concentra la malla). No alinea las líneas con la entalla ni con el punto de carga, y el usuario tendría que escribir la función.
- **Bloques con fusión obligatoria por topología** (lados declarados como comunes). Obliga a describir la conectividad entre bloques; la fusión por posición es la de los generadores habituales y cubre los mismos casos.
- **Una entalla de ancho nulo** (nodos sin fundir a lo largo de una línea). Las entallas de las referencias tienen ancho; con una ranura de una o dos columnas de celdas el caso queda cubierto sin otra opción.
- **Comparaciones también en las cargas.** Una carga discontinua por zonas se declara por grupos; mantener la carga como función aritmética evita errores de integración en los saltos.
- **Calidad sólo en el post-proceso VTK.** El aviso debe llegar antes de resolver.

## Consecuencias

- Un estudio de validación con geometría de bloques es un YAML sin archivos externos: la malla, sus grupos, el análisis y la salida.
- Los elementos con celda `line` también se mallan por bloques en el YAML (un eje y, si hace falta, `mapping` a 2D o 3D).
- Un bloque con lados no alineados con otro (tramos distintos en el lado común) da nodos colgantes sin aviso; la prueba de la parcela o los resultados lo revelan. Se deja así hasta que haga falta un control.
- La primera llamada a `mesh_quality` compila un núcleo pequeño de Numba (queda en caché).
