"""Elementos de interfaz de espesor nulo (ADR 0043).

``InterfaceElement`` es la base pública de los elementos de espesor nulo: un
módulo de usuario la hereda por esta ruta (ADR 0020). ``LocalJumpElement`` es
la base común con la adherencia barra–sólido (``BondSlipInterface``): el bucle
por punto, la masa nula y el registro genérico por punto
(``InterfacePointState``)."""
from solidum.elements.interface._shared import (
    INTEGRATIONS,
    InterfaceElement,
    InterfacePointState,
    LocalJumpElement,
)
from solidum.elements.interface.interface_2d import Interface2D4, Interface2D6
from solidum.elements.interface.interface_3d import Interface3D6, Interface3D8

__all__ = ["INTEGRATIONS", "InterfaceElement", "InterfacePointState", "LocalJumpElement",
           "Interface2D4", "Interface2D6", "Interface3D6", "Interface3D8"]
