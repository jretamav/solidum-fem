"""Tests del descubrimiento automático de componentes (``solidum.autodiscover``).

Regresión de la auditoría H-1.1: el autodiscover usaba ``pkgutil.iter_modules``
no recursivo, así que módulos dentro de subpaquetes (p. ej.
``solidum/elements/solid_2d/embedded_cst.py``) sólo se registraban si alguien
los re-exportaba manualmente desde un ``__init__.py``. El fix migra a
``pkgutil.walk_packages`` (recursivo); este test asegura que toda clase con
decorador ``@<Registry>.register`` en el árbol ``solidum/`` aparece en el
registry correspondiente tras ``_initialize_registries()``.
"""
from __future__ import annotations

import ast
import unittest
from pathlib import Path

import solidum  # dispara autodiscover

from solidum.registry import (
    CohesiveMaterialRegistry,
    ElementRegistry,
    MaterialRegistry,
    SolverRegistry,
)


_FENIX_ROOT = Path(solidum.__file__).resolve().parent

# Registries cuyos miembros esperamos descubrir desde el árbol. ``Quadrature``
# se llena por llamadas explícitas en `math/integration.py` (no por decorador
# de clase) — fuera del alcance de este test.
_REGISTRY_BY_NAME = {
    "ElementRegistry": ElementRegistry,
    "MaterialRegistry": MaterialRegistry,
    "CohesiveMaterialRegistry": CohesiveMaterialRegistry,
    "SolverRegistry": SolverRegistry,
}


def _decorator_targets_registry(decorator: ast.expr) -> str | None:
    """Si el decorador es ``@<X>Registry.register`` (con o sin llamada),
    devuelve ``<X>Registry``; si no, ``None``."""
    target = decorator.func if isinstance(decorator, ast.Call) else decorator
    if isinstance(target, ast.Attribute) and target.attr == "register":
        if isinstance(target.value, ast.Name):
            return target.value.id
    return None


def _scan_decorated_classes(root: Path) -> dict[str, set[str]]:
    """Recorre ``root`` y devuelve ``{registry_name: {class_name, ...}}``
    para todas las clases con ``@<X>Registry.register`` en su lista de
    decoradores. Salta archivos privados (``_*.py``)."""
    found: dict[str, set[str]] = {name: set() for name in _REGISTRY_BY_NAME}
    for py_file in root.rglob("*.py"):
        if any(p.startswith("_") for p in py_file.relative_to(root).parts):
            continue
        try:
            tree = ast.parse(py_file.read_text(encoding="utf-8"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            for deco in node.decorator_list:
                target = _decorator_targets_registry(deco)
                if target in found:
                    found[target].add(node.name)
    return found


class TestAutodiscoverCoversAllDecoratedClasses(unittest.TestCase):
    """Toda clase decorada con ``@<X>Registry.register`` debe estar en el
    registry correspondiente tras ``_initialize_registries()``.

    Si este test falla con clases declaradas pero no registradas, hay un
    hueco en el descubrimiento — típicamente un subpaquete cuyo módulo no
    se está alcanzando recursivamente (regresión de H-1.1).
    """

    def test_all_decorated_classes_are_registered(self):
        declared = _scan_decorated_classes(_FENIX_ROOT)

        # Sanity check: encontramos al menos los componentes esperables.
        # Si el escáner devuelve vacío, el AST no detectó nada y el test
        # no estaría validando nada.
        total_found = sum(len(s) for s in declared.values())
        self.assertGreater(
            total_found, 30,
            f"AST scan found only {total_found} decorated classes — sospecha "
            f"que el escáner no funciona o el árbol está vacío.",
        )

        missing: dict[str, set[str]] = {}
        for registry_name, class_names in declared.items():
            registry = _REGISTRY_BY_NAME[registry_name]
            registered = set(registry.names())
            absent = class_names - registered
            if absent:
                missing[registry_name] = absent

        if missing:
            lines = [
                "Clases decoradas con @<X>Registry.register que NO aparecen "
                "en el registry tras inicialización:",
            ]
            for reg, classes in missing.items():
                lines.append(f"  {reg}: {sorted(classes)}")
            lines.append(
                "Sugiere un hueco en el descubrimiento recursivo "
                "(autodiscover.py) o un alias que el AST no atrapa."
            )
            self.fail("\n".join(lines))


if __name__ == "__main__":
    unittest.main()
