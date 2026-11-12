"""``Quad8`` — cuadrilátero serendípito 2D de orden 2 (8 nodos)."""
from solidum.elements.solid_2d._shared import (
    _HigherOrderSolid2D,
    _batch_kin_quad8,
    _dN_quad8,
    _N_quad8,
)
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Quad8(_HigherOrderSolid2D):
    """Cuadrilátero serendípito 2D de orden 2 (8 nodos).

    Reproduce campos cuadráticos exactamente; sin shear locking severo en
    flexión. Default: Gauss 3×3.
    """
    N_INTEGRATION_POINTS = 9
    BATCH_KINEMATICS = _batch_kin_quad8
    _SHAPE_FN = staticmethod(_N_quad8)
    _GRAD_FN = staticmethod(_dN_quad8)
    _DEFAULT_QUADRATURE = "3x3"
    _MASS_QUADRATURE = "3x3"
    _LOAD_QUADRATURE = "3x3"

    REFERENCE_CELL = "quad"
    # Banda de grieta: h = α·√A. √A es el ancho de una banda de un elemento,
    # como en Quad4. α es la fracción del elemento que se ablanda: con
    # gradiente, la banda se localiza en dos de las tres filas de puntos de
    # Gauss 3×3 (13/18; una de dos con 2×2, ½). Jirásek y Bauer 2012, §4;
    # medido en la viga de Petersson (0.72); usuario, 2026-11-12. Con esfuerzo
    # uniforme el elemento se ablanda entero y α sobrestima la energía (1/α).
    CHARACTERISTIC_LENGTH_FACTOR = 1.0
    CRACK_BAND_LOCALIZATION = {"3x3": 13.0 / 18.0, "2x2": 0.5}
    REFERENCE_NODES = ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0),
                       (0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0))

    EDGE_NODES = (
        (0, 4, 1),
        (1, 5, 2),
        (2, 6, 3),
        (3, 7, 0),
    )
