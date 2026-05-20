"""``Quad8`` — cuadrilátero serendípito 2D de orden 2 (8 nodos)."""
import numpy as np

from solidum.elements.solid_2d._shared import (
    _HigherOrderSolid2D,
    _dN_quad8,
    _N_quad8,
    _quadratic_edge_traction,
)
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Quad8(_HigherOrderSolid2D):
    """Cuadrilátero serendípito 2D de orden 2 (8 nodos).

    Reproduce campos cuadráticos exactamente; sin shear locking severo en
    flexión. Default: Gauss 3×3.
    """
    N_INTEGRATION_POINTS = 9
    _SHAPE_FN = staticmethod(_N_quad8)
    _GRAD_FN = staticmethod(_dN_quad8)
    _DEFAULT_QUADRATURE = "3x3"

    EDGE_NODES = (
        (0, 4, 1),
        (1, 5, 2),
        (2, 6, 3),
        (3, 7, 0),
    )

    def compute_edge_traction(self, edge: int, t_vec: np.ndarray) -> np.ndarray:
        if edge not in (0, 1, 2, 3):
            raise ValueError(f"edge={edge} fuera de rango para Quad8 (0..3).")
        return _quadratic_edge_traction(self, edge, t_vec, n_dofs=16)
