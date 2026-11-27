"""``Interface2D4`` e ``Interface2D6`` — elementos de interfaz de espesor nulo
en 2D con caras lineales de 2 nodos o cuadráticas de 3 (specs
``docs/specs/Interface2D4.md`` e ``Interface2D6.md``, ADR 0043)."""
from __future__ import annotations

from solidum.elements.interface._shared import InterfaceElement
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Interface2D4(InterfaceElement):
    """Interfaz 2D de 4 nodos: cara A = (1, 2), cara B = (3, 4) con 4 frente a
    1 y 3 frente a 2 (orden de COH2D4; antihorario si B está del lado de
    ``n = e_z × s``)."""

    DOF_NAMES = ["ux", "uy"]
    JUMP_DIM = 2
    FACE_CELL = "line"
    N_FACE_NODES = 2
    A_NODES = (0, 1)
    B_NODES = (3, 2)
    N_INTEGRATION_POINTS = 2
    VTK_CELL_TYPE = "quad"


@ElementRegistry.register
class Interface2D6(InterfaceElement):
    """Interfaz 2D de 6 nodos con caras cuadráticas (posiblemente curvas):
    esquinas primero, como los sólidos cuadráticos. Cara A = (1, 2, 5) y cara
    B = (3, 4, 6), con 4 frente a 1, 3 frente a 2 y 6 frente a 5; 5 y 6 son
    los intermedios. Es el orden de la celda VTK 30 (``quad6``)."""

    DOF_NAMES = ["ux", "uy"]
    JUMP_DIM = 2
    FACE_CELL = "line3"
    N_FACE_NODES = 3
    A_NODES = (0, 1, 4)
    B_NODES = (3, 2, 5)
    N_INTEGRATION_POINTS = 3
    VTK_CELL_TYPE = "quad6"
