"""``Interface3D6`` e ``Interface3D8`` — elementos de interfaz de espesor nulo
en 3D con caras triangulares o cuadriláteras (specs ``Interface3D6`` e
``Interface3D8``, extensiones de ``Interface2D4``; ADR 0043)."""
from __future__ import annotations

from solidum.elements.interface._shared import InterfaceElement
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Interface3D6(InterfaceElement):
    """Interfaz 3D de 6 nodos: cara A = (1, 2, 3), cara B = (4, 5, 6) con 4
    frente a 1, 5 frente a 2 y 6 frente a 3 (orden de COH3D6)."""

    DOF_NAMES = ["ux", "uy", "uz"]
    JUMP_DIM = 3
    FACE_CELL = "tri"
    N_FACE_NODES = 3
    A_NODES = (0, 1, 2)
    B_NODES = (3, 4, 5)
    N_INTEGRATION_POINTS = 3
    VTK_CELL_TYPE = "wedge"


@ElementRegistry.register
class Interface3D8(InterfaceElement):
    """Interfaz 3D de 8 nodos: cara A = (1, 2, 3, 4), cara B = (5, 6, 7, 8)
    con 5 frente a 1, etc. (orden de COH3D8)."""

    DOF_NAMES = ["ux", "uy", "uz"]
    JUMP_DIM = 3
    FACE_CELL = "quad"
    N_FACE_NODES = 4
    A_NODES = (0, 1, 2, 3)
    B_NODES = (4, 5, 6, 7)
    N_INTEGRATION_POINTS = 4
    VTK_CELL_TYPE = "hexahedron"
