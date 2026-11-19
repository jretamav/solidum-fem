# ADR 0042 — El elemento conoce su dominio (P10) y el contorno libre es una función del programa principal

- **Estado**: aceptado (2026-11-19), tras la revisión del plan por `arquitecto-solidum` (cambios antes de programar, aplicados; `informes/revision_arquitectura_2026-11-19_paso_b.md`).
- **Fecha**: 2026-11-19
- **Alcance**:
  - `solidum/core/element.py`: propiedad `Element.domain` (referencia débil `_domain_ref`, atributo de clase `None` por omisión) y `Element.__getstate__`, que no serializa la referencia;
  - `solidum/core/domain.py`: `add_element` enlaza y rechaza un elemento de otro dominio vivo; `replace_element` (nuevo) sustituye y desenlaza el anterior; `Domain.__setstate__` vuelve a enlazar al deserializar o copiar los elementos que no tienen dominio vivo (una copia superficial comparte los elementos con el original, que los conserva); `Domain.topology_version`, contador de `add_element` y `replace_element`;
  - `solidum/math/assembly.py`: la huella de topología del ensamblador incluye `topology_version`, así que un elemento sustituido se ve también por lotes (antes la huella era el número de elementos y el camino por lotes seguía con la matriz del elemento anterior, *medido* por el arquitecto);
  - `solidum/mesh/boundary.py` (nuevo): `free_boundary(elements, dim, dofs=None)`, exportada en `solidum.mesh`; la usan `solidum/mesh/structured_blocks.py` (grupos de las mallas estructuradas, que tenía su propia copia) y el seguimiento de grieta del módulo `discontinuities`;
  - `CST_Embedded2D` declara `EDGE_NODES`;
  - pruebas: `tests/test_domain_element_link.py` (incluye un control sintáctico: en `solidum/` sólo `Domain` cambia `domain.elements`), `tests/test_free_boundary.py`; las 8 escrituras directas en `domain.elements[...]` de cuatro archivos de prueba pasan a `add_element`/`replace_element`.
- **Extiende**: el ADR 0020 (pieza nueva **P10** del contrato con los módulos de usuario).

## Contexto

La regla de trayectoria de la grieta que decidió el usuario el 2026-11-19 (spec `CrackTracking`) hace salir la grieta por un lado del **borde libre** del dominio cuando el elemento de la punta lo toca. El seguimiento conoce sólo los elementos de su región (los `CST_Embedded2D` que lo referencian); el borde libre depende de **todos** los elementos del modelo: el lado entre el concreto embebido y una placa de apoyo no es libre, aunque la placa no sea de la región. Ni el elemento, ni sus nodos, ni el lector YAML le daban el dominio al seguimiento; el único gancho posterior a crear los elementos era `apply(domain)` de las restricciones (ADR 0032).

El programa principal ya calculaba el contorno libre (lados o caras de un solo elemento) dentro de `structured_blocks._predicate_groups`; el seguimiento habría sido el segundo consumidor de la misma lógica.

## Decisión

1. **P10: `Element.domain`.** Cada elemento guarda una referencia débil al dominio que lo contiene; la ponen sólo `Domain.add_element` y `Domain.replace_element`. La propiedad devuelve el dominio si sigue vivo **y** el elemento sigue siendo el suyo con ese id (`domain.elements.get(id) is self`); si no, `None`. Así un elemento colado a mano en `elements`, o sustituido, no apunta a un dominio que ya no lo tiene.
2. **Un solo camino de mutación.** `add_element` y `replace_element` son los únicos que cambian `domain.elements`; un elemento no puede estar en dos dominios vivos (`ValueError`). Un control sintáctico de las pruebas lo vigila en `solidum/` (escrituras `x.elements[...] = ...`, `del`, `pop`, `update`, `setdefault`, `clear`, `popitem`).
3. **Serialización y copias.** La referencia no se serializa (`Element.__getstate__`); `Domain.__setstate__` vuelve a enlazar sus elementos. Medido por el arquitecto antes del cambio: con la referencia tal cual, `pickle` de un dominio fallaba (`TypeError`) y `deepcopy` dejaba la copia apuntando al original.
4. **Contorno libre centralizado.** `solidum.mesh.free_boundary(elements, dim, dofs=None)`: lados (`EDGE_NODES`, dim 2) o caras (`FACE_NODES`, dim 3) que declara un solo elemento; el filtro `dofs` deja contar sólo los elementos con esos grados de libertad (los térmicos no tapan el contorno mecánico). Los elementos que no declaran lados (barras, marcos, resortes, interfaces) no cuentan.
5. **`CST_Embedded2D` declara `EDGE_NODES`.** Sin ellos el lado común de dos regiones embebidas, o de la región y un sólido, parecería libre. Declararlos no activa cargas por lado: la carga por grupo exige además `compute_edge_traction` (`solidum/utils/yaml_loads.py`), que el elemento no tiene.
6. **Huella de topología.** `Domain.topology_version` aumenta con cada `add_element` y `replace_element`. La usan el ensamblador (su huella de topología, que ya no confunde un elemento sustituido con el anterior) y el seguimiento, que guarda el borde libre con la referencia débil al dominio y esa huella: añadir una placa entre dos análisis cambia la decisión siguiente (prueba). La caché no se serializa con el seguimiento. Un elemento sustituido en una región de seguimiento reemplaza al anterior como miembro (mismo id).
7. **Uso perezoso.** El seguimiento sólo mira el dominio cuando una punta toca un lado que no comparte con otro elemento de su región; un elemento sin dominio da entonces un error claro (`RuntimeError`, «no pertenece a ningún dominio»). `PrescribedCrackPath` no lo necesita nunca.

## Justificación como pieza del programa principal (Reglas §4)

Un módulo de usuario no puede saber en qué dominio está un elemento: la pieza tiene que estar en el programa principal. Es genérica: la necesita cualquier componente que decida algo con la vista del modelo entero desde un elemento; el contacto, segundo consumidor previsto, tendrá que conocer las superficies libres. El contorno libre ya lo usaba el programa principal (mallas estructuradas): el seguimiento es su segundo consumidor y la función pasa a ser pública.

## Alternativas descartadas

- **Pasar el dominio al seguimiento desde el lector YAML.** No cubre los modelos construidos en Python ni los elementos añadidos después.
- **Un gancho tipo `apply(domain)` para los seguimientos** (como las restricciones, ADR 0032). Otro gancho de ciclo de vida con un solo consumidor; P10 es más pequeño y sirve a cualquier componente.
- **Contar sólo los elementos de tres nodos** (el prototipo). Una placa Quad4 dejaba libre el lado concreto–placa.
- **El borde libre de la región embebida y un caso especial para los miembros.** Con dos seguimientos vecinos el lado común parecía libre (hallazgo del arquitecto).
- **Un `WeakKeyDictionary` lateral elemento → dominio.** 110 B por elemento frente a 0.3–0.7 B del atributo (medido por el arquitecto), y la misma invalidación que resolver.

## Consecuencias

- **A favor**: el seguimiento ve el modelo entero sin que el programa principal lo nombre; el contorno libre tiene una sola implementación; `pickle` y `deepcopy` de un dominio siguen funcionando y la copia queda enlazada con la copia.
- **En contra**: escribir directamente en `domain.elements` deja de ser una vía válida (las 8 escrituras de las pruebas, migradas; el control sintáctico vigila `solidum/`). Un lado de un CST no coincide con el de un elemento cuadrático vecino (malla no conforme): una placa Tri6 o Quad8 pegada no tapa el lado (spec `CrackTracking` §3).
