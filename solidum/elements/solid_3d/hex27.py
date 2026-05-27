"""``Hex27`` — hexaedro Lagrangiano 3D triquadrático (ADR 0012, sub-etapa A.ter).

27 nodos: 8 vértices + 12 medios de arista + 6 centros de cara + 1
centro del cuerpo. Convención VTK_TRIQUADRATIC_HEXAHEDRON. Cuadratura
por defecto Gauss 3×3×3 (27 puntos). Subclase de
:class:`_HigherOrderSolid3D`.

Spec: ``docs/specs/Hex27.md``.
"""
from solidum.elements.solid_2d._shared import _dN_quad9, _N_quad9
from solidum.elements.solid_3d._shared import (
    _HigherOrderSolid3D,
    _dN_hex27,
    _N_hex27,
)
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Hex27(_HigherOrderSolid3D):
    """Hexaedro Lagrangiano 3D triquadrático (27 nodos).

    Producto tensorial completo de Lagrange cuadráticos en las tres
    direcciones naturales. Reproduce **todos** los polinomios hasta
    grado total 6 que caen en el espacio triquadrático
    (``ξ^a η^b ζ^c`` con ``a, b, c ∈ {0, 1, 2}``), incluyendo los
    términos puramente triquadráticos ``ξ²η²ζ²`` que faltan en el
    serendípito ``Hex20``. Análogo 3D del ``Quad9`` 2D.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        27 nodos en orden VTK_TRIQUADRATIC_HEXAHEDRON: 0-7 vértices
        (idéntico Hex8), 8-19 medios de arista (idéntico Hex20), 20-25
        centros de cara en el orden ``(-x, +x, -y, +y, -z, +z)``, 26
        centro del cuerpo.
    material : Material
        Material 3D (``STRAIN_DIM = 6``).
    quadrature : str, optional
        Nombre de cuadratura desde ``QuadratureRegistry``. Default
        ``"hex_3x3x3"`` (27 puntos, full integration para Hex27).
        Alternativa ``"hex_2x2x2"`` (8 puntos, reducida — introduce 27
        modos de hourglass espurios por elemento aislado, sin
        estabilización implementada; la mayoría no propaga en mallas
        estructuradas pero el `Hex27` reducido es la combinación más
        problemática del catálogo 3D).

    Notes
    -----
    Locking volumétrico con ν → 0.5 atenuado respecto al ``Hex8`` y
    comparable al ``Hex20``; el espacio triquadrático añade flexibilidad
    pero no elimina la incompresibilidad sin formulación mixta u-p.
    El espacio extra de ``Hex27`` vs ``Hex20`` mejora la convergencia
    en geometrías con curvatura severa pero el coste computacional es
    mayor (60 → 81 DOFs por elemento, 30% más por dirección).
    """

    N_INTEGRATION_POINTS = 27  # default Gauss 3×3×3

    _SHAPE_FN = staticmethod(_N_hex27)
    _GRAD_FN = staticmethod(_dN_hex27)
    _DEFAULT_QUADRATURE = "hex_3x3x3"
    # Mass quadrature fijada a 3×3×3: integra exactamente el producto
    # triquadrático×triquadrático (grado 4 en cada dirección) sobre
    # Jacobiano constante.
    _MASS_QUADRATURE = "hex_3x3x3"

    # Funciones de forma de cara: Quad9 Lagrangiano 9 nodos (4 vértices
    # + 4 medios de arista + 1 centro de cara). Cuadratura 3×3 (exacta
    # para tracción uniforme sobre cara plana con N cuadráticos).
    _FACE_N_FN = staticmethod(_N_quad9)
    _FACE_DN_FN = staticmethod(_dN_quad9)
    _FACE_QUADRATURE = "3x3"

    # ADR 0012 — 6 caras con normal saliente, paritarias con Hex8/Hex20 en
    # los vértices y medios. Cada cara añade el nodo de centro de cara
    # correspondiente (índice 20-25 según la cara).
    FACE_NODES = (
        (0, 3, 2, 1, 11, 10,  9,  8, 24),  # 0: -ζ (inferior), centro 24
        (4, 5, 6, 7, 12, 13, 14, 15, 25),  # 1: +ζ (superior), centro 25
        (0, 1, 5, 4,  8, 17, 12, 16, 22),  # 2: -η (frontal),  centro 22
        (1, 2, 6, 5,  9, 18, 13, 17, 21),  # 3: +ξ (derecha),  centro 21
        (2, 3, 7, 6, 10, 19, 14, 18, 23),  # 4: +η (trasera),  centro 23
        (3, 0, 4, 7, 11, 16, 15, 19, 20),  # 5: -ξ (izquierda), centro 20
    )
