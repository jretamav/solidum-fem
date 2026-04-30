# Crecimiento del sistema

Los capítulos anteriores describen el estado actual del programa. Este capítulo describe el procedimiento por el que se incorporan piezas nuevas — el aspecto del proyecto que la arquitectura optimiza explícitamente.

## Tres puntos de extensión canónicos

| Componente a incorporar | Ubicación | Mecanismo de registro | Modifica otros archivos |
|---|---|---|---|
| Material | `fenix/materials/<snake>.py` | `@MaterialRegistry.register` | No |
| Elemento | `fenix/elements/<snake>.py` | `@ElementRegistry.register` | No |
| Solver | `fenix/math/solver_<snake>.py` | `@SolverRegistry.register` | No |

La columna "modifica otros archivos" es deliberadamente "No". Si la incorporación de un componente exigiera la modificación del intérprete, el ensamblador o la inicialización, dicha situación constituiría un síntoma de degradación de la arquitectura y procedería un refactor con su correspondiente Architecture Decision Record (ADR).

## Skill `/fenix-new`

Para reducir también la fricción ergonómica de la incorporación, el repositorio incluye una *skill* versionada (`.claude/skills/fenix-new/`) que el asistente invoca cuando el usuario solicita un material, elemento o solver nuevo. La invocación es `/fenix-new material|element|solver <Name>` y produce automáticamente:

- El archivo en su carpeta canónica con la plantilla del tipo correspondiente.
- El decorador `@register` correctamente posicionado.
- Una prueba esqueleto en `tests/` con la estructura habitual (configuración, expectativas, validación contra solución analítica si aplica).

La *skill* no sustituye al razonamiento físico ni a la formulación; constituye estrictamente un mecanismo de generación de plantillas. Cierra el ciclo entre la arquitectura (optimizada para extensión) y la herramienta operativa que la materializa.

## Ciclo especificación → código → catálogo

Cuando el componente a incorporar es físico (no infraestructura), el flujo es:

1. *Apertura de la especificación.* El usuario abre una especificación en `docs/specs/<Nombre>.md` a partir de la plantilla `_template_<kind>.md`. La especificación contiene: especificación física, formulación numérica y contrato YAML (`interface`, `parameters`, `acceptance`, `references`).
2. *Bloqueo en ausencia de especificación completa.* El asistente no escribe código hasta que la especificación esté completa. Si falta información, formula consultas en la sección *Diálogo* de la propia especificación. Esto convierte la especificación, sucesivamente, en orden de trabajo y en referencia detallada del componente.
3. *Implementación.* El asistente implementa el componente y añade pruebas conforme a los criterios de `acceptance`.
4. *Cierre.* Tras la validación (todas las pruebas verdes contra los criterios), el asistente actualiza `status: validated` en la especificación e inserta una entrada en el catálogo correspondiente (`docs/catalogo_<elementos|materiales|solvers>.md`).

La relación es: la especificación constituye la orden de trabajo y la referencia detallada por componente; el catálogo constituye el índice navegable del conjunto. El manual de referencia (`Reference_manual.pdf`) se genera de forma automática a partir de las especificaciones.

## Architecture Decision Records: memoria de arquitectura persistente

Para los cambios que afectan al diseño del programa (no a la implementación interna de un módulo), el asistente crea un Architecture Decision Record en `docs/adr/000N-titulo.md` antes del *commit* o de forma simultánea al mismo. Los ADR son persistentes y consultables; sustituyen a las explicaciones efímeras del *chat* a fin de que el sistema permanezca comprensible meses o años después. Los ADR vigentes se enumeran en el siguiente capítulo.

## Mantenimiento de los documentos vivos

Los archivos `docs/arquitectura.md`, `docs/conceptos_clave.md` y este manual se mantienen sincronizados con el estado del código. El asistente los actualiza cuando se modifica la topología del sistema (introducción de una capa nueva, una carpeta canónica nueva o una responsabilidad transversal nueva), no cuando se modifica la implementación interna de un módulo. Si tras un refactor de tamaño medio cualquiera de estos documentos no refleja la realidad del código, dicha discrepancia constituye un defecto del propio refactor.
