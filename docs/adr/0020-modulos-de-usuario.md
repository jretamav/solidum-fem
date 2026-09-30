# ADR 0020 — Módulos de usuario: las formulaciones no estándar fuera del programa principal

- **Estado**: aceptado (2026-09-24)
- **Fecha**: 2026-09-24
- **Alcance**: nueva carpeta `solidum/user/`; la discontinuidad interior embebida (ADR 0010) pasa a `solidum/user/discontinuities/`; registros (`solidum/registry.py`), lector YAML (`solidum/utils/yaml_parser.py`), ruta gmsh (`solidum/utils/gmsh_parser.py`), autodescubrimiento, herramienta de specs (`solidum/tools/spec.py`), camino por lotes (`element_is_batchable`), reexports de `solidum/__init__.py`, CI, Reglas.md §4, CLAUDE.md, skill `solidum-new`, manual de referencia, `paper.md`. **No cambia** ninguna formulación ni ningún resultado.

## Vocabulario

Dos términos, los de FEAP:

- **Programa principal**: Solidum estándar. Elementos, materiales y solvers de elementos finitos clásicos. Es todo `solidum/` salvo `solidum/user/`.
- **Módulo de usuario**: una formulación no estándar, en su propia carpeta `solidum/user/<nombre>/`. El programa principal no la conoce hasta que un modelo la pide. Es la analogía de los elementos de usuario (`elmtNN`), materiales de usuario (`umati`) y macros de usuario (`umacr`) de FEAP.

## Contexto

El usuario, al decidir completar su formulación de discontinuidades interiores (2026-09-24): *"este elemento no es un elemento finito estándar por lo que su implementación en solidum no debería contaminar su kernel principal. Velo como algo parecido a lo que hacen otros programas, por ejemplo FEAP, en el que se tienen formulaciones FEM estándar pero si el usuario necesita implementar algo nuevo, esto lo hace de tal forma que no contamina el programa principal"*.

Hoy la contamina. Inventario (2026-09-24):

- Fuera de su propio archivo: `core/discontinuity_state.py` y `core/cohesive_material.py` (150 líneas), `CohesiveMaterialRegistry` en `registry.py`, 22 referencias en `yaml_parser.py` (una sección `cohesive_materials` y una clave `cohesive_material` como casos especiales), las constantes `EMBEDDED_LOCAL_JUMP_*`, reexports en la raíz, `autodiscover.py`, `tools/spec.py` y dos capítulos cableados en el manual de referencia.
- El registro por decorador ya permite registrar clases sin tocar el programa principal. Lo que no permite es que un módulo declare su propia **familia de material con sección YAML**, sus **referencias cruzadas**, su **tipo de spec** o su **capítulo de manual**: ahí están todos los casos especiales.
- Lo mismo le pasa hoy a la familia **térmica**: el lector YAML repite tres veces el mismo bloque (mecánica, térmica, cohesiva), los elementos térmicos resuelven su material con una heurística y no pueden salir de gmsh. Ese segundo caso real justifica las piezas genéricas de abajo con independencia de la embebida (Reglas §1).

## Decisión

### 1. Ubicación y carga

- La formulación completa (elemento `CST_Embedded2D`, familia cohesiva y su registro, `DiscontinuityState`, sus constantes y, en adelante, el seguimiento de trayectoria) vive en `solidum/user/discontinuities/`. Una sola instalación: viaja con Solidum.
- **Carga explícita, nunca automática.** `import solidum` no importa nada de `solidum/user/`. Un modelo la pide:
  - YAML: `user_modules: [discontinuities]`, procesado al principio de la lectura, antes de validar.
  - Python: `solidum.load_user_module("discontinuities")`, o importar el módulo directamente.
  - Si un YAML usa un tipo del módulo sin declararlo, el error lo dice y lista los módulos disponibles. Así un modelo no funciona en una PC y falla en otra según lo que se haya importado antes.

### 2. Lo que el programa principal ofrece a un módulo de usuario

Cada pieza se justifica por casos reales ajenos a la embebida.

| # | Pieza | Casos reales que la justifican |
|---|---|---|
| P1 | `Registry` público (hoy `_BaseRegistry`) con metadatos de familia: `YAML_SECTION`, `YAML_LABEL`, `SPEC_KIND`; `Registry.families()`. `__init_subclass__` crea `_items` (hoy olvidarlo comparte el dict de la base en silencio) y rechaza secciones o tipos de spec duplicados. *Precisado el 2026-09-29*: `YAML_LABEL` es la etiqueta de los mensajes de **todo** registro, no sólo de las familias; el atributo privado `_kind` desapareció | Familias mecánica y térmica |
| P2 | Lector YAML genérico: un solo bucle sobre `Registry.families()` en vez de tres copias. Claves de primer nivel desconocidas pasan a ser error (hoy se ignoran en silencio) | Mecánica y térmica; la térmica gana la validación de parámetros que hoy sólo tiene la mecánica |
| P3 | `Element.REFERENCE_KWARGS`: qué parámetros del elemento son referencias a otra familia, p. ej. `{"material": MaterialRegistry}`. El lector valida y resuelve las referencias sin heurística | Elementos térmicos (`{"material": ThermalMaterialRegistry}`) |
| P4 | `solidum.load_user_module(nombre)` y la clave YAML `user_modules`. *Precisado el 2026-09-29*: la carga es **atómica** y **sólo añade** entradas a los registros; si el módulo reasigna una entrada existente (p. ej. registra otra clase como `Elastic2D`), lanza `ValueError` y deshace la carga entera (registros, familias, tipos de spec y `sys.modules`) | Carga de cualquier módulo |
| P5 | Gancho de inicio de paso `Element.prepare_step(U_committed)`, **ya existe**. Se documenta como contrato genérico: idempotente para el mismo `U_committed`, porque tras una bisección el solver lo repite. Un test del programa principal, con una subclase de `Element` sin registrar, comprueba que lo llaman los cinco solvers que deben | Equivale a la tarea de inicio de paso de un elemento de usuario de FEAP |
| P6 | Funciones de cinemática públicas (`compute_kinematics_tri3`, `compute_integrands`). *Precisado el 2026-09-29*: su ruta pública es `solidum.elements.solid_2d` (`from solidum.elements.solid_2d import compute_kinematics_tri3, compute_integrands`), sin guion bajo; el módulo interno `solid_2d/_shared.py` no forma parte del contrato | API para elementos de usuario; hoy las usan Tri3 y Quad4 |
| P7 | Ruta gmsh con tipo de elemento, referencias y parámetros por grupo físico | Los elementos térmicos tampoco pueden salir de gmsh hoy |
| P8 | Contrato del **elemento con estado propio** (historia fuera de `ElementState`): sobrescribe `commit_state` y `compute_global_stiffness`, y no va por el camino por lotes. `element_is_batchable` lo garantiza: devuelve `False` si la clase sobrescribe `compute_element_state` o `commit_state` respecto de la que declaró `BATCH_KINEMATICS` | Hoy es una trampa silenciosa para cualquier elemento que herede de un sólido por lotes |

**Queda fuera, a propósito:**

- El servicio a nivel de modelo para el seguimiento de trayectoria (Cervera et al. 2010). Se diseña con su spec en la fase 3 del plan de la embebida. Diseñarlo ahora sería especulativo; la forma candidata es un método de instancia `Element.step_services()`, con el servicio declarado y configurado en YAML como una familia (P1-P3).
- Una marca `REQUIRES_STEP_PREPARATION` validada por los solvers: el módulo puede protegerse solo (rechaza una evaluación con desplazamientos sin haber pasado por `prepare_step`).
- Un bus de eventos, la carga automática por *entry points* y un estado declarado por tamaño al estilo `nh1/nh3` de FEAP: sin un segundo caso real.

### 3. API pública, YAML y tests

- `from solidum import CST_Embedded2D` (y la ley cohesiva y su registro) **deja de funcionar, sin periodo de transición**: la versión es beta y ningún ejemplo lo usa. El nuevo import es `from solidum.user.discontinuities import CST_Embedded2D`.
- Los YAML que usan la embebida añaden `user_modules: [discontinuities]`.
- Los tests del módulo pasan a `tests/user/discontinuities/`, con nombres de archivo propios (los tests se importan entre hermanos por nombre). Los barridos de contrato del programa principal recorren sólo las clases del programa principal: el módulo se carga en la misma sesión de pytest y sus clases quedan registradas.
- Specs del módulo en `docs/user/discontinuities/specs/`; entrada propia «Módulos de usuario» en los catálogos; apéndice «Módulos de usuario» en el manual de referencia.
- El ADR 0010 se conserva como registro histórico de la formulación, con una nota que remite a éste para su ubicación.
- Reglas.md §4, CLAUDE.md y la skill `solidum-new` incorporan el módulo de usuario como destino, para que el siguiente componente no estándar no acabe en el programa principal.
- `paper.md` (en inglés; desde el 2026-09-24 sólo en local, fuera del repositorio) deja de presentar la embebida como parte del programa y presenta Solidum como extensible por módulos de usuario, con la embebida como primero. El texto lo valida el usuario.

### 4. Verificación de que el programa principal queda limpio

`tests/test_main_program_purity.py`, con tres controles (redacción al día del 2026-09-29; la original tenía dos y listaba a mano los nombres de `discontinuities`):

1. **Estructural** (análisis sintáctico con `ast`): ningún archivo del programa principal —todo `solidum/` salvo los subpaquetes `solidum/user/<nombre>/`— importa un módulo de usuario de ninguna forma: absoluta (`import solidum.user.<x>`, `from solidum.user.<x> import …`, `from solidum.user import <x>`), relativa (`from ..user import x`), perezosa (dentro de una función), por cadena (`importlib.import_module("solidum.user." + n)`, `__import__`, `sys.modules[...]`) o con `load_user_module` y nombre literal. De `solidum.user` sólo se admite el cargador genérico (`load_user_module`, `available_user_modules`). Cada forma de escape se prueba contra el detector, para que el control no quede en verde por un detector ciego.
2. **Léxico con nombres derivados**: ningún archivo del programa principal menciona lo que aporta un módulo. Los nombres no se mantienen a mano: se derivan cargando cada módulo instalado en un proceso limpio (las entradas que añade a los registros, las familias que declara con su `YAML_SECTION` y su `SPEC_KIND`, las referencias nuevas de sus elementos, su `__all__` y su ruta de importación). Un módulo nuevo no se da de alta en el test; `NOMBRES_ADICIONALES` es sólo un complemento opcional para constantes que la derivación no ve.
3. **En ejecución**: en un subproceso limpio, tras `import solidum`, ningún `solidum.user.*` en `sys.modules`, toda clase registrada pertenece al programa principal y sólo existen las familias `materials` y `thermal_materials`.

La otra dirección la vigila, desde el 2026-09-29, `tests/test_user_module_public_imports.py`: un módulo de usuario (y sus tests) sólo importa rutas y nombres públicos del programa principal —nada con guion bajo— y sus subclases de `Registry` no declaran atributos `_*` (`_items` lo crea la base). Así el programa principal puede reorganizar sus internos sin romper un módulo.

## Migración

Siete pasos, con la suite en verde en cada uno. Los pasos 1 a 5 generalizan el programa principal **con la embebida todavía dentro**, así que el paso 6 se reduce a mover archivos.

| Paso | Contenido |
|---|---|
| 1 | P1 y P2: familias genéricas; la cohesiva sigue en su sitio como tercera familia |
| 2 | P3: desaparece el caso especial `cohesive_material` del lector |
| 3 | P6, P8 y `tools/spec.py` por `SPEC_KIND` |
| 4 | P4, con un módulo de juguete probado en subproceso |
| 5 | P7, ruta gmsh |
| 6 | Test de P5 en el programa principal; mover los tests y después el código a `solidum/user/discontinuities/`; retirar reexports y la línea de autodescubrimiento |
| 7 | Test de pureza, documentos, protocolo, manuales y `paper.md` |

Estimación: dos o tres sesiones.

## Alternativas descartadas

- **Dejarla donde está**: contradice la directiva del usuario y hace que cada componente no estándar siguiente (seguimiento de trayectoria, otra ley cohesiva, XFEM) entre también al programa principal.
  > *Nota 2026-09-29*: desde el **criterio de admisión** de Reglas.md §1 (una formulación es estándar si aparece en los textos de referencia del área y está implementada en un código general establecido), XFEM y las leyes cohesivas tracción-salto genéricas con elementos de interfaz son **estándar**: Belytschko-Liu-Moran y de Borst-Crisfield-Remmers-Verhoosel las tratan y Abaqus las implementa. Ejemplos correctos de no estándar son la discontinuidad embebida con cinemática KOS de Retama (2010), el seguimiento de trayectoria de grieta y el control por liberación de energía de Gutiérrez (2004). La familia cohesiva sigue en `discontinuities` porque su único consumidor es el elemento embebido; pasará al programa principal, sin física nueva, cuando entre un elemento de interfaz estándar.
- **Paquete Python aparte** (`solidum_embedded`, instalado por separado): obliga a instalar dos paquetes en cada PC y su nombre sugiere un programa distinto de Solidum. El usuario lo rechazó.
- **Repositorio aparte**: dos historias que sincronizar por Drive, sin beneficio con un solo autor.
- **Carga automática al importar Solidum**: haría depender el comportamiento de lo instalado o importado antes y el YAML no diría qué formulación necesita.

## Consecuencias

- **A favor**: el programa principal no conoce la embebida; un módulo de usuario nuevo no toca el lector YAML, gmsh ni la herramienta de specs; desaparecen tres copias de código; la familia térmica gana validación de parámetros y referencias, y puede salir de gmsh; la embebida se convierte en el primer ejemplo real de cómo se escribe un módulo de usuario.
- **En contra**: `from solidum import CST_Embedded2D` deja de funcionar; los YAML de la embebida necesitan `user_modules`; P1-P8 pasan a ser contrato estable del programa principal y cambiarlos afecta a todo módulo de usuario.

## Resultado

Migración hecha el 2026-09-24, con la suite en verde en cada paso y sin cambio de formulación ni de resultados. Los pasos 6 y 7 se partieron en dos commits.

| Paso | Commit | Contenido | Suite (pasan / skipped) |
|---|---|---|---|
| 1 | `9f49c6e` | P1 y P2: `Registry` público con familias; un solo bucle en el lector YAML; sección de primer nivel desconocida = error | 1568 / 9 |
| 2 | `47565a1` | P3: `Element.REFERENCE_KWARGS`; fuera la heurística de material y el caso especial `cohesive_material` | 1570 / 9 |
| 3 | `1f69c47` | P6 (cinemática pública), P8 (`element_is_batchable`) y specs localizadas por `SPEC_KIND` | 1572 / 9 |
| 4 | `953e211` | P4: `solidum/user/`, `load_user_module`, `available_user_modules` y clave `user_modules`, probados con un módulo de juguete en subproceso | 1574 / 9 |
| 5 | `ad95ad4` | P7: ruta gmsh con `mesh_element` y `element` por grupo físico | 1581 / 9 |
| 6a | `cdd8125` | Test de P5 (`tests/test_prepare_step_hook.py`); barridos de contrato restringidos a las clases del programa principal | 1589 / 9 |
| 6b | `67fcef7` | La embebida pasa a `solidum/user/discontinuities/`; fuera reexports, autodescubrimiento, registro y constantes del programa principal; tests a `tests/user/discontinuities/test_disc_*.py`, specs a `docs/user/discontinuities/specs/` | 1610 / 10 |
| 7a | `eca7ec9` | Test de pureza `tests/test_main_program_purity.py`; docstrings del programa principal sin nombrar la embebida | 1613 / 10 |

El resto del paso 7 (Reglas.md, CLAUDE.md, skill `solidum-new`, catálogos, manuales y `paper.md`) acompaña a esta sección.

Diferencias con el texto de la decisión, sin cambiarla:

- **P5**: el gancho lo invocan **seis** solvers, no cinco: `NonlinearSolver`, `ArcLengthSolver`, `IndirectDisplacementSolver`, `DissipationArcLengthSolver`, `NewtonNewmarkSolver` y `NewtonHHTSolver`. `LinearSolver` no lo invoca.
- **§4, test de pureza**: además de los nombres listados, vigila la ruta `solidum.user.discontinuities` y la palabra «embebida»/«embedded», y comprueba que tras `import solidum` sólo existen las familias `materials` y `thermal_materials`.

Recuento tras la migración. Programa principal: 22 elementos, 13 materiales mecánicos, 1 material térmico y 14 solvers. Módulo de usuario `discontinuities`: 1 elemento (`CST_Embedded2D`) y 1 material cohesivo (`CohesiveDamageIsotropic`, familia `cohesive_materials`). Suite completa: 1613 pasan / 10 skipped.

### 2026-09-29 — segundo módulo y endurecimiento de la frontera

Auditoría de pureza del programa principal, tras fijar el usuario el **criterio de admisión** (Reglas.md §1) y delegar por completo la arquitectura (Reglas.md §2-§3). Sin cambio de formulación ni de resultados.

- **Segundo módulo de usuario, `dissipation_control`**: `DissipationArcLengthSolver` (Gutiérrez 2004, con conmutación cilíndrico↔disipación a la manera de Verhoosel, Remmers y Gutiérrez 2009) sale del programa principal por el criterio de admisión: no cumple la condición (ii), porque ningún código general establecido implementa esa restricción del arco de forma nativa. Código en `solidum/user/dissipation_control/`, spec en `docs/user/dissipation_control/specs/`, tests en `tests/user/dissipation_control/` (`test_dctl_*`), entrada en la sección «Módulos de usuario — dissipation_control» del catálogo de solvers. Se carga con `user_modules: [dissipation_control]`; `from solidum.math.solvers import DissipationArcLengthSolver` deja de funcionar. Es el primer módulo que aporta un solver, como un macro de usuario (`umacr`) de FEAP.
- **Pieza genérica nueva, `solidum.math.solvers.ArcProblem`**: el problema de Newton por paso del arc-length, exportado con nombre público. Una variante de arco que sólo cambia la restricción del paso lo extiende y redefine su método `constraint`; hay dos casos reales, `IndirectDisplacementSolver` en el programa principal y `DissipationArcLengthSolver` en su módulo.
- **Pista de `user_modules` también para solvers**: un tipo de solver desconocido en el YAML sugiere declarar `user_modules` y lista los módulos disponibles, como ya hacían los tipos de elemento y las secciones de familia.
- **Retirada de `_negative_pivots`**: la aproximación por signo del determinante de `DissipationArcLengthSolver` y el marcador vacío de `ArcLengthSolver`; ningún `solve` los llamaba. El seguimiento de la inercia de `K_t` sigue esperando a un LDLᵀ verdadero (deuda #7 de STATUS).
- **Endurecimiento de la frontera**: P4 atómica y de sólo añadir (ver §2); test de pureza con tres controles y nombres derivados en vez de la lista a mano (ver §4; la palabra «embebida»/«embedded» ya no se vigila: los nombres vigilados son los derivados más `NOMBRES_ADICIONALES`); test nuevo `tests/test_user_module_public_imports.py` para la dirección contraria; ruta pública de P6 (`solidum.elements.solid_2d`); `YAML_LABEL` como etiqueta de todo registro, sin `_kind`.
- **Specs del programa principal con bloque `admission`**: cada spec declara en su contrato los textos (`textbooks`) y códigos (`general_codes`) que la admiten, con los identificadores canónicos de `solidum/tools/spec.py` (`ADMISSION_TEXTBOOKS`, `ADMISSION_GENERAL_CODES`).

Cifras al día, sin reescribir las de arriba (históricas): el gancho P5 lo invocan **cinco** solvers del programa principal (`NonlinearSolver`, `ArcLengthSolver`, `IndirectDisplacementSolver`, `NewtonNewmarkSolver`, `NewtonHHTSolver`) más `DissipationArcLengthSolver` en su módulo; tras una bisección lo repiten los estáticos que bisecan (los tres primeros y el del módulo), mientras que los dos dinámicos abortan el paso que no converge. Programa principal: 22 elementos, 13 materiales mecánicos, 1 material térmico y **13 solvers**. Módulos de usuario: **2**, `discontinuities` (1 elemento, 1 material cohesivo) y `dissipation_control` (1 solver).
