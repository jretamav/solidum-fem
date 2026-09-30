"""Módulos de usuario (ADR 0020).

Formulaciones no estándar que viven fuera del programa principal, como los
elementos, materiales y macros de usuario de FEAP (``elmtNN``, ``umati``,
``umacr``). Cada módulo es un subpaquete ``solidum/user/<nombre>/`` que se
conecta a las piezas genéricas del programa principal —registros, familias de
material con su sección YAML, referencias declaradas por el elemento, gancho
de inicio de paso— y que el programa principal no importa nunca. Lo carga el
modelo que lo pide:

- en YAML, ``user_modules: [<nombre>]``;
- en Python, ``solidum.load_user_module("<nombre>")``.

La carga es explícita a propósito: el YAML dice qué formulación necesita, y un
modelo no funciona en una sesión y falla en otra según lo que se haya
importado antes.

Un módulo de usuario sólo **añade** entradas a los registros; nunca sustituye
una existente. Si lo intenta (p. ej. registra otra clase como ``Elastic2D``),
la carga falla y se deshace: cambiaría en silencio el comportamiento del
programa principal para todos los modelos posteriores de la sesión.
"""
import importlib
import pkgutil
import sys


def available_user_modules() -> list:
    """Nombres de los módulos de usuario instalados (subpaquetes públicos)."""
    return sorted(name for _, name, ispkg in pkgutil.iter_modules(__path__)
                  if ispkg and not name.startswith("_"))


def load_user_module(name: str):
    """Carga el módulo de usuario ``name`` y lo devuelve.

    Importa recursivamente sus submódulos, como hace el programa principal
    con los suyos, para que sus decoradores registren elementos, materiales y
    familias. Cargarlo dos veces no hace nada.

    La carga es atómica: si falla —por un error del propio módulo o porque
    reasigna una entrada ya registrada (``ValueError``)— los registros y
    ``sys.modules`` vuelven a como estaban, y un segundo intento vuelve a
    fallar en vez de devolver un módulo a medio registrar."""
    disponibles = available_user_modules()
    if name not in disponibles:
        raise ValueError(
            f"Módulo de usuario '{name}' desconocido. Disponibles: {disponibles}."
        )
    from solidum.autodiscover import _discover_package
    package = f"{__name__}.{name}"
    snapshot = _snapshot()
    try:
        _discover_package(package)
        module = importlib.import_module(package)
        reassigned = _reassigned(snapshot)
        if reassigned:
            raise ValueError(
                f"Módulo de usuario '{name}': reasigna entradas ya registradas "
                f"({'; '.join(reassigned)}). Un módulo de usuario sólo puede "
                f"añadir entradas a los registros, nunca sustituirlas: da a sus "
                f"clases nombres propios. La carga se ha deshecho."
            )
    except BaseException:
        _restore(snapshot)
        raise
    return module


# ---------------------------------------------------------------------------
# Instantánea de los registros para la carga atómica
# ---------------------------------------------------------------------------

def _registries() -> list:
    """Todos los registros vivos: subclases de ``Registry`` a cualquier
    profundidad, incluidas las familias que declaran los módulos de usuario."""
    from solidum.registry import Registry
    found, pending = [], list(Registry.__subclasses__())
    while pending:
        reg = pending.pop()
        if reg not in found:
            found.append(reg)
            pending.extend(reg.__subclasses__())
    return found


def _snapshot() -> dict:
    """Entradas de cada registro, familias, tipos de spec y módulos de
    usuario importados, para detectar reasignaciones y deshacer la carga."""
    from solidum.registry import QuadratureRegistry, Registry
    return {
        "items": {reg: dict(reg._items) for reg in _registries()},
        "quadratures": dict(QuadratureRegistry._rules),
        "families": list(Registry._families),
        "spec_kinds": dict(Registry._spec_kinds),
        "modules": {m for m in sys.modules if m.startswith(f"{__name__}.")},
    }


def _qualified(obj) -> str:
    return (f"{getattr(obj, '__module__', '?')}."
            f"{getattr(obj, '__qualname__', type(obj).__name__)}")


def _reassigned(snapshot: dict) -> list:
    """Entradas de la instantánea que ahora apuntan a otro objeto."""
    from solidum.registry import QuadratureRegistry
    found = []
    for reg, items in snapshot["items"].items():
        for key, old in items.items():
            new = reg._items.get(key, old)
            if new is not old:
                found.append(f"{reg.__name__} '{key}': {_qualified(old)} → {_qualified(new)}")
    for key, old in snapshot["quadratures"].items():
        if QuadratureRegistry._rules.get(key, old) is not old:
            found.append(f"QuadratureRegistry '{key}'")
    return found


def _restore(snapshot: dict) -> None:
    """Devuelve los registros y ``sys.modules`` al estado de la instantánea."""
    from solidum.registry import QuadratureRegistry, Registry
    for reg in _registries():
        reg._items.clear()
        reg._items.update(snapshot["items"].get(reg, {}))
    QuadratureRegistry._rules.clear()
    QuadratureRegistry._rules.update(snapshot["quadratures"])
    Registry._families[:] = snapshot["families"]
    Registry._spec_kinds.clear()
    Registry._spec_kinds.update(snapshot["spec_kinds"])
    parent = sys.modules[__name__]
    for m in [m for m in sys.modules
              if m.startswith(f"{__name__}.") and m not in snapshot["modules"]]:
        del sys.modules[m]
        if m.count(".") == __name__.count(".") + 1:
            parent.__dict__.pop(m.rsplit(".", 1)[1], None)
