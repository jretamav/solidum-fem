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
from solidum.elements.solid_3d._shared import (
    _compute_gradient_kinematics_hex8,
    _shape_functions_hex8,
)
from solidum.elements.thermal._shared import _ThermalSolid
from solidum.registry import ElementRegistry


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
    N_INTEGRATION_POINTS: ClassVar[int] = 8  # default; la instancia lo ajusta

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

    def compute_face_flux(self, face: int, q_bar: float) -> np.ndarray:
        """Vector nodal del flujo impuesto sobre una cara (Neumann).

        Integra ``-∫_Γ q̄ Nᵀ dΓ`` sobre la cara cuadrilateral indicada con
        cuadratura Gauss 2×2 y funciones de forma bilineales. Exacto para
        cara plana y flujo constante.

        Parameters
        ----------
        face : int
            Índice según ``FACE_NODES`` (0..5), numeración con normal
            saliente paritaria con el ``Hex8`` mecánico.
        q_bar : float
            Flujo normal impuesto [W/m²]. **Positivo = SALIENTE** del
            dominio (enfriamiento). Un flujo entrante se declara negativo.

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

        face_local = self.FACE_NODES[face]
        X = np.array(
            [self.nodes[i].coordinates[:3] for i in face_local],
            dtype=np.float64,
        )

        gp = 1.0 / np.sqrt(3.0)
        s_pts = (-gp, gp, gp, -gp)
        t_pts = (-gp, -gp, gp, gp)
        s_sign = (-1.0, 1.0, 1.0, -1.0)
        t_sign = (-1.0, -1.0, 1.0, 1.0)

        f = np.zeros(self.N_NODES)
        for sp, tp in zip(s_pts, t_pts):
            N = np.zeros(4)
            dN_ds = np.zeros(4)
            dN_dt = np.zeros(4)
            for i in range(4):
                N[i] = 0.25 * (1.0 + s_sign[i] * sp) * (1.0 + t_sign[i] * tp)
                dN_ds[i] = 0.25 * s_sign[i] * (1.0 + t_sign[i] * tp)
                dN_dt[i] = 0.25 * t_sign[i] * (1.0 + s_sign[i] * sp)

            g_s = dN_ds @ X
            g_t = dN_dt @ X
            dA = float(np.linalg.norm(np.cross(g_s, g_t)))

            for i, node_local in enumerate(face_local):
                f[node_local] += -N[i] * float(q_bar) * dA

        return f
