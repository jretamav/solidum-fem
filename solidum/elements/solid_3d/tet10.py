"""``Tet10`` — tetraedro cuadrático isoparamétrico (ADR 0012, sub-etapa A.ter).

10 nodos: 4 vértices + 6 medios de arista. Convención VTK_QUADRATIC_TETRA.
Cuadratura por defecto Stroud 4 puntos (``tet_4``), orden 2 exacto —
suficiente para K con Jacobiano constante. Masa consistente integrada
con cuadratura ``tet_15`` (Keast 15 puntos, orden 5) independientemente
de la cuadratura de K, para que la masa quede exacta.

Subclase de :class:`_HigherOrderSolid3D` (la misma base que Hex20/Hex27);
sólo cambian funciones de forma, cuadratura y formas/caras del elemento.

Spec: ``docs/specs/Tet10.md``.
"""
from solidum.elements.solid_2d._shared import _dN_tri6, _N_tri6
from solidum.elements.solid_3d._shared import (
    _HigherOrderSolid3D,
    _dN_tet10,
    _N_tet10,
)
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Tet10(_HigherOrderSolid3D):
    """Tetraedro cuadrático 3D de 10 nodos.

    Reproduce campos cuadráticos completos en 3D sobre el simplex
    tetraédrico. Análogo 3D del ``Tri6`` 2D. Adecuado para mallas no
    estructuradas (donde el ``Hex20``/``Hex27`` requeriría mapeo
    estructurado).

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        10 nodos en orden VTK_QUADRATIC_TETRA: 0-3 vértices (idéntico
        Tet4: orden VTK_TETRA con volumen positivo), 4-9 medios de
        aristas:
            4: arista 0-1, 5: arista 1-2, 6: arista 2-0,
            7: arista 0-3, 8: arista 1-3, 9: arista 2-3.
    material : Material
        Material 3D (``STRAIN_DIM = 6``).
    quadrature : str, optional
        Nombre de cuadratura desde ``QuadratureRegistry``. Default
        ``"tet_4"`` (Stroud 4 puntos, orden 2 — suficiente para K con
        Jacobiano constante). Alternativa ``"tet_15"`` (Keast 15 puntos,
        orden 5) para problemas donde el Jacobiano varía (geometría
        distorsionada o materiales no lineales severos).

    Notes
    -----
    El ``Tet10`` con cuadratura ``tet_4`` integra la K exacta pero
    subintegra la masa consistente (integrando de grado 4). El elemento
    fuerza ``_MASS_QUADRATURE = "tet_15"`` para que la masa sea siempre
    exacta independientemente de la cuadratura de K.
    """

    N_INTEGRATION_POINTS = 4  # default Stroud 4 puntos

    _SHAPE_FN = staticmethod(_N_tet10)
    _GRAD_FN = staticmethod(_dN_tet10)
    _DEFAULT_QUADRATURE = "tet_4"
    # Mass quadrature fijada a tet_15: integra exactamente el producto
    # cuadrático×cuadrático (grado 4) sobre el tetraedro.
    _MASS_QUADRATURE = "tet_15"

    # Funciones de forma de cara: Tri6 cuadrático 6 nodos (3 vértices +
    # 3 medios de arista). Cuadratura tri_3 (3 puntos, orden 2) — exacta
    # para tracción uniforme sobre cara plana con N cuadráticos: la
    # integral de N_i (grado 2) sobre el triángulo se integra
    # exactamente.
    _FACE_N_FN = staticmethod(_N_tri6)
    _FACE_DN_FN = staticmethod(_dN_tri6)
    _FACE_QUADRATURE = "tri_3"

    # ADR 0012 — 4 caras triangulares con normal saliente, paritarias
    # con Tet4 en los vértices. Cada cara añade 3 medios de arista en
    # orden Tri6 (vértices 0, 1, 2 antihorarios; medios entre 0-1, 1-2,
    # 2-0). Cara i opuesta al vértice i.
    FACE_NODES = (
        # Face 0: opuesta al vértice 0, vértices Tet10 (1, 2, 3).
        # Tri6 medios: medio 0-1 (Tri6) = arista (1,2) Tet10 = nodo 5;
        #              medio 1-2 (Tri6) = arista (2,3) Tet10 = nodo 9;
        #              medio 2-0 (Tri6) = arista (3,1) Tet10 = nodo 8.
        (1, 2, 3, 5, 9, 8),
        # Face 1: opuesta al vértice 1, vértices Tet10 (0, 3, 2).
        # medios: 0-3 = nodo 7; 3-2 = nodo 9; 2-0 = nodo 6.
        (0, 3, 2, 7, 9, 6),
        # Face 2: opuesta al vértice 2, vértices Tet10 (0, 1, 3).
        # medios: 0-1 = nodo 4; 1-3 = nodo 8; 3-0 = nodo 7.
        (0, 1, 3, 4, 8, 7),
        # Face 3: opuesta al vértice 3, vértices Tet10 (0, 2, 1).
        # medios: 0-2 = nodo 6; 2-1 = nodo 5; 1-0 = nodo 4.
        (0, 2, 1, 6, 5, 4),
    )
