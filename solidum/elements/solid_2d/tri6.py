"""``Tri6`` — triángulo 2D de 6 nodos (cuadrático completo P₂)."""
from solidum.elements.solid_2d._shared import (
    _HigherOrderSolid2D,
    _batch_kin_tri6,
    _dN_tri6,
    _N_tri6,
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
    BATCH_KINEMATICS = _batch_kin_tri6
    _SHAPE_FN = staticmethod(_N_tri6)
    _GRAD_FN = staticmethod(_dN_tri6)
    _QUADRATURE_FAMILY = "tri"
    _DEFAULT_QUADRATURE = "tri_3"
    _MASS_QUADRATURE = "tri_6"
    _LOAD_QUADRATURE = "tri_6"

    REFERENCE_CELL = "tri"
    # Banda de grieta: √(2A), el ancho de una banda de un elemento, como en
    # Tri3. Sin factor de localización (α = 1, con aviso): con gradiente unos
    # triángulos de la banda se ablandan en dos de sus tres puntos y otros
    # enteros, así que la fracción no es constante. El ⅔ de la literatura
    # (Slobbe, Hendriks y Rots 2013; DIANA 10.5 §93.1.7.4) disipaba 1.17-1.28
    # veces G_f·A en una viga entallada a flexión; sin él, 0.79-0.87 (del lado
    # seguro). Decisión del usuario, 2026-11-12.
    CHARACTERISTIC_LENGTH_FACTOR = 2.0 ** 0.5
    CRACK_BAND_LOCALIZATION: dict = {}
    REFERENCE_NODES = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0),
                       (0.5, 0.0), (0.5, 0.5), (0.0, 0.5))

    EDGE_NODES = (
        (0, 3, 1),
        (1, 4, 2),
        (2, 5, 0),
    )
