# ADR 0035 — Mallado desde la geometría con Gmsh, y lector de Gmsh para todas las celdas

- **Estado**: aceptado (2026-11-11)
- **Fecha**: 2026-11-11
- **Alcance**:
  - `solidum/utils/gmsh_parser.py`: `GmshParser` importa todas las celdas de los sólidos de Solidum, en 2D y 3D;
  - `solidum/mesh/gmsh_geometry.py` (nuevo): `mesh_from_geometry`, `GeometryMeshOptions`, `validate_geometry_block`;
  - `solidum/utils/yaml_parser.py`: la clave `mesh` admite un `.geo` o un bloque con `geometry`; `mesh_quadrature` deja de imponer `2x2`;
  - `pyproject.toml`: grupo opcional `mesh` (`gmsh`).
- **Extiende**: el ADR 0022 (contrato `REFERENCE_CELL`/`REFERENCE_NODES`), el ADR 0020, P7 (receta de elemento por grupo físico) y el ADR 0023 (cargas por grupo físico).
- **Salda**: la deuda #18 de `docs/STATUS.md`.

## Contexto

Al validar la banda de grieta en la viga de Petersson (2026-11-11) las mallas se generaron con la API de Gmsh desde scripts que construían el modelo a mano, fuera del YAML y de `solidum.run`, y sin los archivos de resultados de Solidum. El camino normal de un usuario —un `.msh` declarado en el YAML— tenía dos huecos: el lector sólo importaba `quad` y `triangle` en 2D (ni elementos cuadráticos ni 3D, deuda #18), y el `.msh` había que generarlo por fuera, de modo que un estudio de malla (h = 10, 5, 2.5 mm) o el paso a segundo orden obligaba a regenerar archivos a mano. Dos ejemplos de `examples/` se omiten en las pruebas por esa razón.

## Decisión

1. **Lector de Gmsh para todas las celdas.** `triangle`, `triangle6`, `quad`, `quad8`, `quad9`, `tetra`, `tetra10`, `hexahedron`, `hexahedron20` y `hexahedron27` dan por omisión `Tri3`, `Tri6`, `Quad4`, `Quad8`, `Quad9`, `Tet4`, `Tet10`, `Hex8`, `Hex20` y `Hex27`. La dimensión de la malla es la mayor de sus celdas; sus nodos llevan 2 ó 3 coordenadas.
2. **Orden de nodos por coordenadas de referencia, sin tablas por tipo.** Una tabla única da, para cada celda de meshio, las coordenadas de referencia de sus nodos (las de `gmsh.model.mesh.getElementProperties` con la permutación `_gmsh_to_meshio_order` de meshio; una prueba las contrasta con Gmsh). El lector casa esas coordenadas con las `REFERENCE_NODES` del elemento destino (ADR 0022), que usan las mismas celdas de referencia. Un elemento nuevo con `REFERENCE_NODES` entra sin tocar el lector; uno sin ellas recibe los nodos en el orden de la celda (comportamiento anterior).
3. **Orientación.** Un sólido con jacobiano negativo en el centro de su geometría lineal (celda horaria en 2D, invertida en 3D) se reordena reflejando su celda de referencia (ξ ↔ η), con la misma correspondencia por coordenadas; uno degenerado es un error.
4. **Celdas de dimensión menor.** Líneas en 2D y caras en 3D son contorno (cargas y condiciones por grupo, ADR 0023). Sólo son elementos si su grupo físico nombra un `element` en el propio grupo (no el heredado de `mesh_element`): así entran barras incrustadas en un sólido.
5. **Cuadratura.** Sin `mesh_quadrature`, cada elemento usa su regla por omisión; antes se imponía `2x2`, que no sirve a `Tri6` ni a los sólidos 3D.
6. **Mallado desde la geometría.** `mesh` admite la ruta de un `.geo` o un bloque `{geometry, parameters, order, element_shape, serendipity, output}`. Solidum llama a Gmsh: fija los parámetros del `.geo` (`gmsh.parser.setNumber` antes de `gmsh.merge`; el `.geo` los declara con `DefineConstant` y su valor por omisión), recombina en cuadriláteros si se pide, malla en la dimensión de la geometría, eleva a orden 2 (serendípito por omisión: `Quad8`, `Hex20`), escribe el `.msh` junto al YAML (o en `output`) y lo lee con el lector del punto 1. Un parámetro que no aparece en el `.geo` es un error (errata). El bloque se valida antes de abrir Gmsh.
7. **Gmsh opcional.** Se importa sólo al mallar; sin él, el mensaje indica instalarlo (`pip install gmsh`, grupo `mesh`) o declarar el `.msh`. meshio sigue siendo el lector de `.msh`, sin dependencia nueva obligatoria.

## Alternativas descartadas

- **Tablas de permutación escritas a mano por tipo de celda.** Son la fuente clásica de errores silenciosos (Hex20 y Tet10 difieren entre Gmsh, meshio y VTK) y obligan a tocar el lector en cada elemento nuevo.
- **Leer el `.msh` con la API de Gmsh en lugar de meshio.** Haría obligatoria la dependencia de Gmsh para cualquier malla; meshio ya estaba y basta con conocer su orden.
- **Generar la malla en el preprocesador del usuario siempre.** Deja el estudio de malla y el segundo orden fuera del modelo, sin reproducibilidad.

## Consecuencias

- Un estudio de malla es cambiar `parameters` en el YAML; el paso a segundo orden, `order: 2`.
- Las pruebas que mallan con Gmsh se omiten si Gmsh no está instalado (la integración continua no lo instala todavía; las del lector usan `.msh` versionados en `tests/data/gmsh/` y corren siempre).
- Las mallas de primer orden ya existentes se leen igual (mismos nodos, elementos y orden); sólo cambia que una celda horaria se reorienta y que `mesh_quadrature` ausente deja la regla del elemento (la misma `2x2` en `Quad4`).
