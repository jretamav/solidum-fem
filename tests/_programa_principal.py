"""Filtro compartido por los barridos de contrato del programa principal.

Los módulos de usuario (``solidum/user/``, ADR 0020) se cargan en la misma
sesión de pytest cuando corren sus propios tests, y sus clases quedan en los
registros globales. Cada módulo trae su propio barrido; los barridos de
``tests/`` cubren sólo el programa principal, y todos filtran con esta
función para que el criterio sea uno solo.
"""
from __future__ import annotations


def nombres_programa_principal(registry) -> list:
    """Nombres registrados en ``registry`` cuya clase pertenece al programa
    principal (su módulo no está bajo ``solidum.user.``)."""
    return [n for n in registry.names()
            if not registry.get(n).__module__.startswith("solidum.user.")]
