# solidum_fem/solidum/elements/thermal/hex8_thermal.py
"""``Hex8Thermal`` — hexaedro trilineal de conducción 3D (Etapa 8).

Spec: ``docs/specs/Hex8Thermal.md``. Hermano 3D del ``Quad4Thermal``, con
el que comparte toda la formulación vía ``_ThermalSolid``; aquí sólo vive
la geometría trilineal y el flujo de frontera por cara.
"""
from typing import ClassVar, List

import numpy as np

from solidum.core.node import Node
from solidum.core.thermal_material import ThermalMaterial
from solidum.elements.solid_2d._shared import _dN_quad4, _shape_functions_quad4
from solidum.elements.solid_3d._shared import (
    _batch_grad_hex8,
    _compute_gradient_kinematics_hex8,
    _shape_functions_hex8,
)
from solidum.elements.thermal._shared import _ThermalSolid
from solidum.math.distributed_loads import face_load
from solidum.registry import ElementRegistry, QuadratureRegistry


@ElementRegistry.register
class Hex8Thermal(_ThermalSolid):
    """Hexaedro trilineal isoparamétrico de conducción de calor 3D.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        8 nodos en orden VTK_HEXAHEDRON (volumen positivo).
    material : ThermalMaterial
        Material térmico 3D (``FLUX_DIM = 3``), p. ej. ``ThermalConduction``.
    quadrature : str, optional
        Clave de ``QuadratureRegistry``. Default ``"hex_2x2x2"``.

    Notes
    -----
    **Sin `thickness`**: el volumen sale de la geometría, igual que en los
    sólidos 3D mecánicos. Es la diferencia visible más inmediata respecto
    al hermano 2D.

    La cuadratura reducida ``"hex_1x1x1"`` está disponible pero **no se
    recomienda**: deja ``rango(K_e) = 3`` frente a los 7 necesarios, con
    **4 modos hourglass** de temperatura (frente al único del
    ``Quad4Thermal`` reducido). Sin estabilización.
    """

    FLUX_DIM: ClassVar[int] = 3
    N_NODES: ClassVar[int] = 8
    DEFAULT_QUADRATURE: ClassVar[str] = "hex_2x2x2"
    QUADRATURE_FAMILY: ClassVar[str] = "hex"
    CAPACITY_QUADRATURE: ClassVar[str] = "hex_2x2x2"
    N_INTEGRATION_POINTS: ClassVar[int] = 8  # default; la instancia lo ajusta
    # Camino por lotes (ADR 0014): B térmica = ∂N/∂x (3×8); thickness = 1.
    BATCH_KINEMATICS = _batch_grad_hex8
    BATCH_SCALE = "thickness"
    # Celda y nodos de referencia, los del Hex8 mecánico (ADR 0022).
    REFERENCE_CELL: ClassVar[str] = "hex"
    REFERENCE_NODES: ClassVar[tuple] = (
        (-1.0, -1.0, -1.0), (1.0, -1.0, -1.0), (1.0, 1.0, -1.0), (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0), (1.0, -1.0, 1.0), (1.0, 1.0, 1.0), (-1.0, 1.0, 1.0),
    )
    # Caras bilineales para el flujo de frontera, como en el Hex8 mecánico.
    _FACE_N_FN = staticmethod(_shape_functions_quad4)
    _FACE_DN_FN = staticmethod(_dN_quad4)
    _FACE_QUADRATURE: ClassVar[str] = "2x2"

    # Caras con normal saliente, paritarias con las del Hex8 mecánico
    # (ADR 0012): 0(−ζ) 1(+ζ) 2(−η) 3(+ξ) 4(+η) 5(−ξ).
    FACE_NODES: ClassVar[tuple] = (
        (0, 3, 2, 1),  # 0: −ζ (inferior)
        (4, 5, 6, 7),  # 1: +ζ (superior)
        (0, 1, 5, 4),  # 2: −η (frontal)
        (1, 2, 6, 5),  # 3: +ξ (derecha)
        (2, 3, 7, 6),  # 4: +η (trasera)
        (3, 0, 4, 7),  # 5: −ξ (izquierda)
    )

    def __init__(
        self,
        element_id: int,
        nodes: List[Node],
        material: ThermalMaterial,
        quadrature: str | None = None,
    ):
        # thickness = 1.0 heredado de la base: en 3D el volumen es
        # geométrico y el factor no debe intervenir.
        super().__init__(element_id, nodes, material, quadrature=quadrature)

    # ------------------------------------------------------------------
    # Contrato geométrico
    # ------------------------------------------------------------------

    def _shape_functions(self, point) -> np.ndarray:
        xi, eta, zeta = point
        return _shape_functions_hex8(xi, eta, zeta)

    def _gradient(self, point, coords):
        xi, eta, zeta = point
        return _compute_gradient_kinematics_hex8(xi, eta, zeta, coords)

    # ------------------------------------------------------------------
    # Flujo de frontera
    # ------------------------------------------------------------------

    def compute_face_flux(self, face: int, q_bar) -> np.ndarray:
        """Vector nodal del flujo impuesto sobre una cara (Neumann).

        Integra ``-∫_Γ q̄ Nᵀ dΓ`` sobre la cara cuadrilateral indicada, con
        su geometría bilineal real (alabeada o no) y cuadratura Gauss 2×2
        (:func:`solidum.math.distributed_loads.face_load`). Exacto para
        flujo constante sobre una cara plana y, sobre cualquier cara, para
        ``q̄ = g·n`` con ``g`` constante, porque ``g·n dA`` es
        ``g·(∂x/∂s × ∂x/∂r) ds dr``, polinómico.

        Parameters
        ----------
        face : int
            Índice según ``FACE_NODES`` (0..5), numeración con normal
            saliente paritaria con el ``Hex8`` mecánico.
        q_bar : float o callable
            Flujo normal impuesto [W/m²]: constante o campo ``q̄(x, n)`` de
            la posición y de la normal exterior unitaria de la cara (ADR
            0022). **Positivo = SALIENTE** del dominio (enfriamiento). Un
            flujo entrante se declara negativo.

        Returns
        -------
        np.ndarray
            Vector de 8 componentes, no nulo sólo en los 4 nodos de la
            cara. Para ``q_bar`` uniforme sobre cara plana de área ``A``,
            reparte ``-q̄·A/4`` a cada uno.
        """
        if face not in range(len(self.FACE_NODES)):
            raise ValueError(
                f"Hex8Thermal(id={self.id}): face={face} fuera de rango; "
                f"esperado 0..{len(self.FACE_NODES) - 1}."
            )

        local = list(self.FACE_NODES[face])
        pts, wts = QuadratureRegistry.get(self._FACE_QUADRATURE)
        N = np.array([self._FACE_N_FN(s, r) for s, r in pts], dtype=np.float64)
        dN = np.array([self._FACE_DN_FN(s, r) for s, r in pts], dtype=np.float64)
        X = self.get_coordinate_matrix(ndim=3)[local]
        f = np.zeros(self.N_NODES)
        f[local] = -face_load(N, dN, np.asarray(wts, dtype=np.float64), X, q_bar,
                              n_components=1)
        return f
