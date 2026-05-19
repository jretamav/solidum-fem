"""Descubrimiento automático de materiales, elementos y solvers.

Importa **recursivamente** todos los submódulos de los paquetes
registrables. Los decoradores `@MaterialRegistry.register`,
`@ElementRegistry.register`, etc. ejecutan el registro como efecto
secundario de la importación, así no hay que mantener listas manuales
de imports.

El paseo es recursivo (``pkgutil.walk_packages``) para que componentes
que viven dentro de subpaquetes (p. ej.
``fenix/elements/solid_2d/embedded_cst.py``) se registren sin depender
de re-exports manuales en el ``__init__.py`` del subpaquete o de
``fenix/__init__.py``. Se saltan los módulos cuyo último componente
empieza con ``_`` (convención del proyecto para módulos privados
internos como ``_shared``); no tienen registros y no aportan al
descubrimiento.
"""
import importlib
import pkgutil


def _is_external_dependency_missing(exc: ModuleNotFoundError) -> bool:
    """True si el ``ModuleNotFoundError`` viene de una librería externa.

    Módulos opt-in de Fenix (p. ej. ``fenix.math.linalg.cholesky`` requiere
    ``scikit-sparse``) lanzan ``ModuleNotFoundError`` al importarse si la
    dependencia no está instalada — caso legítimo y absorbido por el
    ``__init__`` del subpaquete. No queremos que esto rompa el
    descubrimiento del resto. Pero un ``ModuleNotFoundError`` con
    ``name='fenix.algo_que_no_existe'`` SÍ debe propagar (bug real).
    """
    return exc.name is not None and not exc.name.startswith("fenix")


def _discover_package(package_name: str) -> None:
    """Importa recursivamente todos los submódulos públicos de un paquete."""
    package = importlib.import_module(package_name)
    if not hasattr(package, "__path__"):
        return
    for _, modname, _ in pkgutil.walk_packages(
        package.__path__, prefix=f"{package_name}.",
    ):
        # Saltar módulos privados (p. ej. fenix.elements.solid_2d._shared).
        if any(part.startswith("_") for part in modname.split(".")):
            continue
        try:
            importlib.import_module(modname)
        except ModuleNotFoundError as exc:
            if _is_external_dependency_missing(exc):
                # Dependencia externa opcional ausente; el subpaquete lo
                # absorbe (ver fenix.math.linalg.cholesky para el patrón).
                continue
            raise


def initialize() -> None:
    """Carga todos los materiales, elementos, solvers y cuadraturas registrables.

    `fenix.math` se descubre como paquete completo: cubre `solvers/`,
    `integration.py`, `assembly.py` y cualquier solver futuro que el
    usuario añada con su propio decorador, incluso si vive en un
    subpaquete anidado.
    """
    _discover_package("fenix.materials")
    _discover_package("fenix.cohesive_materials")
    _discover_package("fenix.elements")
    _discover_package("fenix.math")
