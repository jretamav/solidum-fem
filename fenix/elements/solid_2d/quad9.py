"""``Quad9`` — cuadrilátero Lagrangiano 2D de orden 2 (9 nodos)."""
import numpy as np

from fenix.elements.solid_2d._shared import (
    _HigherOrderSolid2D,
    _dN_quad9,
    _N_quad9,
    _quadratic_edge_traction,
)
from fenix.registry import ElementRegistry


@ElementRegistry.register
class Quad9(_HigherOrderSolid2D):
    """Cuadrilátero Lagrangiano 2D de orden 2 (9 nodos)."""
    N_INTEGRATION_POINTS = 9
    _SHAPE_FN = staticmethod(_N_quad9)
    _GRAD_FN = staticmethod(_dN_quad9)
    _DEFAULT_QUADRATURE = "3x3"

    EDGE_NODES = (
        (0, 4, 1),
        (1, 5, 2),
        (2, 6, 3),
        (3, 7, 0),
    )

    def compute_edge_traction(self, edge: int, t_vec: np.ndarray) -> np.ndarray:
        if edge not in (0, 1, 2, 3):
            raise ValueError(f"edge={edge} fuera de rango para Quad9 (0..3).")
        return _quadratic_edge_traction(self, edge, t_vec, n_dofs=18)
