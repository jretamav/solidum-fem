"""Registros del módulo de discontinuidades (ADR 0010, ADR 0020).

La ley cohesiva que gobierna la superficie de discontinuidad pertenece a la
familia de materiales de interfaz del programa principal
(``InterfaceMaterialRegistry``, sección YAML ``interface_materials``,
ADR 0033): pasó allí cuando el programa principal tuvo un elemento de
interfaz propio (Reglas.md §4). Aquí queda sólo lo propio del módulo.

Sólo declara atributos públicos de :class:`solidum.registry.Registry`: el
almacén de clases lo crea la base al definir la subclase.
"""
from __future__ import annotations

from solidum.registry import Registry


class CrackTrackerRegistry(Registry):
    """Algoritmos de seguimiento de grieta: deciden qué elementos embebidos se
    agrietan y con qué partición (spec ``GlobalCrackTracking``).

    Familia con sección YAML propia (``crack_trackers``) para que un modelo
    pueda fijar sus parámetros y referenciarlo desde los elementos
    (``crack_tracker: <id>``, ADR 0020, P3). Sin referencia, cada región
    conexa de elementos embebidos recibe uno con los valores por omisión.
    """
    SPEC_KIND = "crack_tracker"
    YAML_SECTION = "crack_trackers"
    YAML_LABEL = "seguimiento de grieta"
