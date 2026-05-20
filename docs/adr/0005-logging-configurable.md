# ADR 0005 — Logging configurable en lugar de `print` directo

**Fecha**: 2026-04-30
**Estado**: Aceptado (2026-04-30)

## Resumen ejecutivo

Solidum sustituye los ~36 `print(...)` diseminados por solvers, parsers, exportadores y elementos por un logger jerárquico estándar (`logging.getLogger("solidum")`). El comportamiento visible por defecto se preserva: los mensajes siguen apareciendo en stdout con el mismo contenido. Lo que cambia es que ahora son **silenciables, redirigibles y filtrables por nivel** sin parchear código fuente. Habilita los frontends GUI y los scripts batch paramétricos previstos en la visión kernel + frontends sin acumular más deuda con cada solver nuevo.

## Contexto

### Estado actual

`solidum/math/solvers.py`, `solidum/utils/yaml_parser.py`, `solidum/utils/gmsh_parser.py`, `solidum/utils/vtk_exporter.py`, `solidum/elements/frame.py`, `solidum/elements/frame3d.py` y `solidum/core/domain.py` emiten progreso, advertencias y errores con `print(...)` directo a stdout. Total: 36 ocurrencias.

### Problema

1. **No silenciable.** Un script que corra cien casos paramétricos llena stdout con el log de cien `--- INICIANDO SOLVER LINEAL ---`. No hay mecanismo para callar la salida sin tocar el código fuente.
2. **No redirigible.** Un frontend GUI quiere capturar el log en un panel propio, no en stdout. Hoy requiere monkey-patching de `print` o redirección de `sys.stdout`.
3. **Sin niveles.** Mensajes informativos y errores se mezclan en el mismo flujo. No se puede filtrar "solo advertencias" o "solo errores".
4. **Compone deuda.** Cada solver nuevo (Riks, dinámica explícita, modal) añade más `print`. Mantener el patrón equivocado encarece todas las extensiones futuras.

### Por qué actuar ahora

A diferencia de otras mejoras arquitecturales en backlog (capa Analysis, sistema de cargas, partición del NonlinearSolver), el coste de logging compone con cada componente N+1: cada solver nuevo añade más prints. Saldarla pronto evita que la deuda crezca; saldarla más tarde implica también convertir los prints introducidos entre tanto. Decisión coherente con Reglas.md §1.

## Decisión

### 1. Módulo `solidum.logging` — logger jerárquico estándar

Un módulo público fino que:

- Define el logger raíz `"solidum"` con un `StreamHandler` por defecto a stdout y formato corto (`[NIVEL] mensaje`).
- Expone `get_logger(name)` para que cada submódulo pida su logger hijo (`solidum.solvers`, `solidum.parsers.yaml`, etc.). La jerarquía permite filtrar por subsistema (silenciar solo el solver, p. ej.) sin mecanismos ad hoc.
- Expone `set_log_level(level)` como helper público de conveniencia para usuarios finales: `solidum.set_log_level("WARNING")` silencia el progreso normal y deja solo advertencias y errores.
- **No llama a `logging.basicConfig()`**: configurar solo el logger `"solidum"` para no contaminar la configuración de logging de aplicaciones que embeban Solidum (frontends GUI).

### 2. Nivel asignado a cada mensaje

| Categoría | Nivel | Ejemplos |
|---|---|---|
| Progreso normal del solver | `INFO` | `--- INICIANDO SOLVER LINEAL ---`, `[PASO 3] Factor: 0.6`, iteración por iteración |
| Convergencia alcanzada / acelerando incremento | `INFO` | `-> CONVERGENCIA ALCANZADA` |
| Lectura de modelos / parsers | `INFO` | `Leyendo modelo desde: foo.yaml`, `Nodos importados: 312` |
| Advertencias de validación | `WARNING` | integración reducida, hourglass, advertencias de elementos |
| Cholesky degrada a LU | `WARNING` | fallback automático, no es error |
| Matriz singular, raíces imaginarias en arc-length | `ERROR` | divergencia que aborta el paso |

Default `INFO`: preserva el comportamiento visible actual sin que el usuario cambie nada.

### 3. API pública

```python
from solidum import set_log_level
set_log_level("WARNING")  # silencia el progreso, deja avisos y errores
```

Acceso programático completo:

```python
import logging
logging.getLogger("solidum.solvers").setLevel(logging.DEBUG)
```

Re-exportada `set_log_level` y `get_logger` desde `solidum/__init__.py`.

### 4. Migración

Sustitución mecánica de cada `print(...)` por `logger.info(...)` / `logger.warning(...)` / `logger.error(...)` con el logger del subsistema correspondiente. Los mensajes mantienen su contenido textual exacto. Total: 7 archivos modificados, una decisión arquitectural única.

## Consecuencias

### Positivas

- **Frontends GUI viables** sin parchear stdout. El logger admite handlers personalizados (`addHandler(GuiPanelHandler())`).
- **Scripts batch silenciables**. `set_log_level("ERROR")` para correr cien casos sin ruido.
- **Solvers nuevos heredan la disciplina**. Cualquier solver futuro pide `get_logger(__name__)` y emite con niveles, sin que la decisión vuelva a discutirse.
- **Sin sorpresas de comportamiento**. Default visible idéntico al actual.

### Negativas

- **Una capa más entre el código y stdout**. Trivial, pero existe.
- **Tests que capturen stdout no ven el log por defecto**. El logger escribe en stdout vía StreamHandler — sigue siendo capturable, pero los tests que usen `capsys` deben adaptarse si hubiera asserts sobre el contenido. No hay tests de ese tipo en la suite actual.

### Neutras

- **Dependencia de `logging` de stdlib**. Cero coste, ya disponible.

## Alternativas consideradas

1. **Mantener `print` y filtrar con un flag global `VERBOSE`**. Rechazada: no soluciona la redirección a frontends, no soporta niveles, sigue componiendo deuda.
2. **Librería externa (`loguru`, `structlog`)**. Rechazada: aporta features que no necesitamos hoy y añade dependencia. `logging` de stdlib es suficiente y es lo que cualquier frontend o script Python sabe consumir nativamente.
3. **No hacer nada hasta que llegue el primer frontend**. Rechazada: la deuda compone con cada solver. Saldar ahora es coste fijo bajo; saldar después es coste creciente.

## Referencias

- Reglas.md §1 — coste del componente N+1, deuda que compone vs deuda estática.
- Visión kernel + frontends del usuario.
- Python `logging` de la stdlib.
