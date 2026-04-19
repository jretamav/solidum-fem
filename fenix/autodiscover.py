"""Descubrimiento automático de materiales, elementos y solvers.

Importa todos los submódulos de los paquetes registrables. Los decoradores
`@MaterialRegistry.register`, `@ElementRegistry.register`, etc. ejecutan
el registro como efecto secundario de la importación, así no hay que
mantener listas manuales de imports.
"""
import importlib
import pkgutil


def _discover_package(package_name: str) -> None:
    """Importa todos los submódulos directos de un paquete."""
    package = importlib.import_module(package_name)
    if not hasattr(package, "__path__"):
        return
    for _, modname, _ in pkgutil.iter_modules(package.__path__):
        importlib.import_module(f"{package_name}.{modname}")


def initialize() -> None:
    """Carga todos los materiales, elementos, solvers y cuadraturas registrables."""
    _discover_package("fenix.materials")
    _discover_package("fenix.elements")
    importlib.import_module("fenix.math.solvers")
    importlib.import_module("fenix.math.integration")
