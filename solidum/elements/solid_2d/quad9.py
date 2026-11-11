"""``Quad9`` — cuadrilátero Lagrangiano 2D de orden 2 (9 nodos)."""
from solidum.elements.solid_2d._shared import (
    _HigherOrderSolid2D,
    _batch_kin_quad9,
    _dN_quad9,
    _N_quad9,
)
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Quad9(_HigherOrderSolid2D):
    """Cuadrilátero Lagrangiano 2D de orden 2 (9 nodos)."""
    N_INTEGRATION_POINTS = 9
    BATCH_KINEMATICS = _batch_kin_quad9
    _SHAPE_FN = staticmethod(_N_quad9)
    _GRAD_FN = staticmethod(_dN_quad9)
    _DEFAULT_QUADRATURE = "3x3"
    _MASS_QUADRATURE = "3x3"
    _LOAD_QUADRATURE = "3x3"

    REFERENCE_CELL = "quad"
    # Banda de grieta: √A, como Quad4, no el ½·√A de Abaqus para cuadráticos:
    # la deformación se localiza en el elemento entero y con ½·√A la prueba de
    # objetividad disipaba 1.97·G_f·A (usuario, 2026-11-11, regla Q2).
    CHARACTERISTIC_LENGTH_FACTOR = 1.0
    REFERENCE_NODES = ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0),
                       (0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0),
                       (0.0, 0.0))

    EDGE_NODES = (
        (0, 4, 1),
        (1, 5, 2),
        (2, 6, 3),
        (3, 7, 0),
    )
