"""``Tri6`` — triángulo 2D de 6 nodos (cuadrático completo P₂)."""
import numpy as np

from solidum.elements.solid_2d._shared import (
    _HigherOrderSolid2D,
    _dN_tri6,
    _N_tri6,
    _quadratic_edge_traction,
)
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Tri6(_HigherOrderSolid2D):
    """Triángulo 2D de 6 nodos (cuadrático completo P₂).

    Cura el shear locking severo del Tri3; reproduce campos cuadráticos
    exactamente. Cuadratura del elemento ``tri_3`` (3 puntos en los
    puntos medios). Masa consistente con cuadratura ``tri_6`` (Dunavant,
    6 puntos, orden 4) porque la del elemento subintegra el producto
    cuadrático×cuadrático y dejaría modos nulos espurios en M.

    Caveat (auditoría H-2.7): la integración con ``tri_3`` (orden 2)
    es exacta para ``Bᵀ·C·B·detJ`` cuando el Jacobiano es constante
    (malla con lados rectos). En **mallas curvilíneas o distorsionadas**
    el integrando es polinomial de orden superior a 2 y ``tri_3``
    subintegra ligeramente — efecto típicamente beneficioso para
    aliviar locking volumétrico, pero puede introducir modos espurios
    sutiles en regímenes patológicos. Si la integración exacta de la
    rigidez es prioridad, pasar ``quadrature=tri_6_points``
    explícitamente al constructor (heredado de
    :class:`_HigherOrderSolid2D`).
    """
    N_INTEGRATION_POINTS = 3
    _SHAPE_FN = staticmethod(_N_tri6)
    _GRAD_FN = staticmethod(_dN_tri6)
    _DEFAULT_QUADRATURE = "tri_3"
    _MASS_QUADRATURE = "tri_6"

    EDGE_NODES = (
        (0, 3, 1),
        (1, 4, 2),
        (2, 5, 0),
    )

    def compute_edge_traction(self, edge: int, t_vec: np.ndarray) -> np.ndarray:
        if edge not in (0, 1, 2):
            raise ValueError(f"edge={edge} fuera de rango para Tri6 (0..2).")
        return _quadratic_edge_traction(self, edge, t_vec, n_dofs=12)
