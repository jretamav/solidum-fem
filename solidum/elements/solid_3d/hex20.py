"""``Hex20`` — hexaedro serendípito 3D de orden 2 (ADR 0012, sub-etapa A.ter).

20 nodos: 8 vértices + 12 medios de arista. Sin nodos de cara ni
centroide. Convención VTK_QUADRATIC_HEXAHEDRON. Cuadratura por defecto
Gauss 3×3×3 (27 puntos). Subclase de :class:`_HigherOrderSolid3D` que
comparte el bucle de Gauss con el ``Hex27`` y, en sub-fase 3, con el
``Tet10``.

Spec: ``docs/specs/Hex20.md``.
"""
from solidum.elements.solid_2d._shared import _dN_quad8, _N_quad8
from solidum.elements.solid_3d._shared import (
    _HigherOrderSolid3D,
    _dN_hex20,
    _N_hex20,
)
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Hex20(_HigherOrderSolid3D):
    """Hexaedro serendípito 3D de 20 nodos (orden 2).

    Reproduce campos cuadráticos completos en 3D (todos los polinomios
    de grado total ≤ 2 más algunos de orden superior, sin los puramente
    triquadráticos). Análogo 3D del ``Quad8`` 2D.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        20 nodos en orden VTK_QUADRATIC_HEXAHEDRON: 0-7 vértices
        (idéntico Hex8), 8-11 medios de aristas de la cara inferior,
        12-15 medios de aristas de la cara superior, 16-19 medios de
        aristas verticales.
    material : Material
        Material 3D (``STRAIN_DIM = 6``).
    quadrature : str, optional
        Nombre de cuadratura desde ``QuadratureRegistry``. Default
        ``"hex_3x3x3"`` (27 puntos, full integration para serendipity
        3D). Alternativa ``"hex_2x2x2"`` (8 puntos, reducida — introduce
        6 modos de hourglass espurios por elemento aislado, sin
        estabilización implementada).

    Notes
    -----
    Locking volumétrico con ν → 0.5 atenuado respecto al ``Hex8`` pero
    aún presente; sin mitigación implementada. Hourglass con integración
    reducida ``hex_2x2x2`` (6 modos espurios por elemento aislado)
    declarado y blindado por test.
    """

    N_INTEGRATION_POINTS = 27  # default Gauss 3×3×3

    _SHAPE_FN = staticmethod(_N_hex20)
    _GRAD_FN = staticmethod(_dN_hex20)
    _DEFAULT_QUADRATURE = "hex_3x3x3"
    # Mass quadrature fijada a 3×3×3 para que la masa siempre quede exacta
    # incluso si el usuario elige reducida para K.
    _MASS_QUADRATURE = "hex_3x3x3"

    # Funciones de forma de cara: Quad8 serendípito 8 nodos (4 vértices +
    # 4 medios de arista). Cuadratura 3×3 (exacta para tracción uniforme
    # sobre cara plana con N cuadráticos).
    _FACE_N_FN = staticmethod(_N_quad8)
    _FACE_DN_FN = staticmethod(_dN_quad8)
    _FACE_QUADRATURE = "3x3"

    # ADR 0012 — 6 caras con normal saliente, paritarias con Hex8 en los
    # vértices. Cada cara añade 4 nodos medios en el mismo orden cíclico
    # que la convención VTK_QUADRATIC_HEXAHEDRON.
    FACE_NODES = (
        (0, 3, 2, 1, 11, 10,  9,  8),  # 0: -ζ (inferior)
        (4, 5, 6, 7, 12, 13, 14, 15),  # 1: +ζ (superior)
        (0, 1, 5, 4,  8, 17, 12, 16),  # 2: -η (frontal)
        (1, 2, 6, 5,  9, 18, 13, 17),  # 3: +ξ (derecha)
        (2, 3, 7, 6, 10, 19, 14, 18),  # 4: +η (trasera)
        (3, 0, 4, 7, 11, 16, 15, 19),  # 5: -ξ (izquierda)
    )
