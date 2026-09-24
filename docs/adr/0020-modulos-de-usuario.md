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
| P1 | `Registry` público (hoy `_BaseRegistry`) con metadatos de familia: `YAML_SECTION`, `YAML_LABEL`, `SPEC_KIND`; `Registry.families()`. `__init_subclass__` crea `_items` (hoy olvidarlo comparte el dict de la base en silencio) y rechaza secciones o tipos de spec duplicados | Familias mecánica y térmica |
| P2 | Lector YAML genérico: un solo bucle sobre `Registry.families()` en vez de tres copias. Claves de primer nivel desconocidas pasan a ser error (hoy se ignoran en silencio) | Mecánica y térmica; la térmica gana la validación de parámetros que hoy sólo tiene la mecánica |
| P3 | `Element.REFERENCE_KWARGS`: qué parámetros del elemento son referencias a otra familia, p. ej. `{"material": MaterialRegistry}`. El lector valida y resuelve las referencias sin heurística | Elementos térmicos (`{"material": ThermalMaterialRegistry}`) |
| P4 | `solidum.load_user_module(nombre)` y la clave YAML `user_modules` | Puerta de carga de cualquier módulo |
| P5 | Gancho de inicio de paso `Element.prepare_step(U_committed)`, **ya existe**. Se documenta como contrato genérico: idempotente para el mismo `U_committed`, porque tras una bisección el solver lo repite. Un test del programa principal, con una subclase de `Element` sin registrar, comprueba que lo llaman los cinco solvers que deben | Equivale a la tarea de inicio de paso de un elemento de usuario de FEAP |
| P6 | Funciones de cinemática públicas (`compute_kinematics_tri3`, `compute_integrands`) | API para elementos de usuario; hoy las usan Tri3 y Quad4 |
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
- `paper.md` (público, en inglés) deja de presentar la embebida como parte del programa y presenta Solidum como extensible por módulos de usuario, con la embebida como primero. El texto lo valida el usuario.

### 4. Verificación de que el programa principal queda limpio

`tests/test_main_program_purity.py`:

1. Ningún archivo de `solidum/` fuera de `solidum/user/` importa `solidum.user` ni menciona los nombres del módulo (`CST_Embedded2D`, `CohesiveMaterial`, `CohesiveDamageIsotropic`, `DiscontinuityState`, `EMBEDDED_LOCAL_JUMP`, `cohesive_material`).
2. En un subproceso limpio, tras `import solidum`: ningún `solidum.user.*` en `sys.modules`, y toda clase registrada pertenece al programa principal.

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
- **Paquete Python aparte** (`solidum_embedded`, instalado por separado): obliga a instalar dos paquetes en cada PC y su nombre sugiere un programa distinto de Solidum. El usuario lo rechazó.
- **Repositorio aparte**: dos historias que sincronizar por Drive, sin beneficio con un solo autor.
- **Carga automática al importar Solidum**: haría depender el comportamiento de lo instalado o importado antes y el YAML no diría qué formulación necesita.

## Consecuencias

- **A favor**: el programa principal no conoce la embebida; un módulo de usuario nuevo no toca el lector YAML, gmsh ni la herramienta de specs; desaparecen tres copias de código; la familia térmica gana validación de parámetros y referencias, y puede salir de gmsh; la embebida se convierte en el primer ejemplo real de cómo se escribe un módulo de usuario.
- **En contra**: `from solidum import CST_Embedded2D` deja de funcionar; los YAML de la embebida necesitan `user_modules`; P1-P8 pasan a ser contrato estable del programa principal y cambiarlos afecta a todo módulo de usuario.
