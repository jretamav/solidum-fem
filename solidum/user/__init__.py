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
"""
import importlib
import pkgutil


def available_user_modules() -> list:
    """Nombres de los módulos de usuario instalados (subpaquetes públicos)."""
    return sorted(name for _, name, ispkg in pkgutil.iter_modules(__path__)
                  if ispkg and not name.startswith("_"))


def load_user_module(name: str):
    """Carga el módulo de usuario ``name`` y lo devuelve.

    Importa recursivamente sus submódulos, como hace el programa principal
    con los suyos, para que sus decoradores registren elementos, materiales y
    familias. Cargarlo dos veces no hace nada."""
    disponibles = available_user_modules()
    if name not in disponibles:
        raise ValueError(
            f"Módulo de usuario '{name}' desconocido. Disponibles: {disponibles}."
        )
    from solidum.autodiscover import _discover_package
    package = f"{__name__}.{name}"
    _discover_package(package)
    return importlib.import_module(package)
