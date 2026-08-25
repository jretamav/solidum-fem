# solidum_fem/solidum/elements/thermal/quad4_thermal.py
"""``Quad4Thermal`` — cuadrilátero bilineal de conducción 2D (Etapa 8).

Spec: ``docs/specs/Quad4Thermal.md``. Misma geometría, funciones de forma y
cuadratura que el ``Quad4`` mecánico; lo que cambia es el campo interpolado
—la temperatura, un DOF escalar por nodo— y la ecuación discretizada.
"""
from typing import ClassVar, List

import numpy as np

from solidum.core.node import Node
from solidum.core.thermal_material import ThermalMaterial
from solidum.elements.solid_2d._shared import (
    _compute_gradient_kinematics_quad4,
    _shape_functions_quad4,
)
from solidum.elements.thermal._shared import _ThermalSolid
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Quad4Thermal(_ThermalSolid):
    """Cuadrilátero bilineal isoparamétrico de conducción de calor 2D.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        4 nodos en sentido antihorario (asegura ``det J > 0``).
    material : ThermalMaterial
        Material térmico 2D (``FLUX_DIM = 2``), p. ej. ``ThermalConduction``.
    thickness : float, optional
        Espesor de la rebanada plana [m]. Default 1.0, que reproduce el
        planteamiento "por unidad de profundidad" habitual en la literatura
        térmica. El flujo total y la capacidad escalan con él.
    quadrature : str, optional
        Clave de ``QuadratureRegistry``. Default ``"2x2"``.

    Notes
    -----
    Es una clase distinta del ``Quad4`` mecánico, no un modo de operación
    suyo: ``DOF_NAMES`` es atributo de clase y el registro depende de poder
    leerlo sin construir el objeto. Es también la separación que adoptan
    Abaqus (``CPE4`` / ``DC2D4``) y ANSYS (``PLANE182`` / ``PLANE55``).

    La cuadratura reducida ``"1x1"`` está disponible pero **no se
    recomienda**: deja ``rango(K_e) = 2`` frente a los 3 necesarios, con un
    modo *hourglass* de temperatura no detectado. Sin estabilización.
    """

    FLUX_DIM: ClassVar[int] = 2
    N_NODES: ClassVar[int] = 4
    DEFAULT_QUADRATURE: ClassVar[str] = "2x2"
    N_INTEGRATION_POINTS: ClassVar[int] = 4  # default; la instancia lo ajusta

    # Bordes del cuadrilátero, paritarios con los del Quad4 mecánico.
    EDGE_NODES: ClassVar[tuple] = ((0, 1), (1, 2), (2, 3), (3, 0))

    def __init__(
        self,
        element_id: int,
        nodes: List[Node],
        material: ThermalMaterial,
        thickness: float = 1.0,
        quadrature: str | None = None,
    ):
        if thickness <= 0.0:
            raise ValueError(
                f"Quad4Thermal(id={element_id}): thickness={thickness} debe ser "
                "estrictamente positivo."
            )
        self.thickness = float(thickness)
        super().__init__(element_id, nodes, material, quadrature=quadrature)

    # ------------------------------------------------------------------
    # Contrato geométrico
    # ------------------------------------------------------------------

    def _shape_functions(self, point) -> np.ndarray:
        xi, eta = point
        return _shape_functions_quad4(xi, eta)

    def _gradient(self, point, coords):
        xi, eta = point
        return _compute_gradient_kinematics_quad4(xi, eta, coords)

    # ------------------------------------------------------------------
    # Flujo de frontera
    # ------------------------------------------------------------------

    def compute_edge_flux(self, edge: int, q_bar: float) -> np.ndarray:
        """Vector nodal del flujo impuesto sobre un borde (Neumann).

        Parameters
        ----------
        edge : int
            Índice según ``EDGE_NODES``: 0=(n0,n1), 1=(n1,n2), 2=(n2,n3),
            3=(n3,n0). Numeración paritaria con el ``Quad4`` mecánico.
        q_bar : float
            Flujo normal impuesto [W/m²]. **Positivo = SALIENTE** del
            dominio (enfriamiento), coherente con la normal exterior y con
            ``q = -k·∇T``. Un flujo entrante se declara negativo.

        Returns
        -------
        np.ndarray
            Vector de 4 componentes. Para ``q_bar`` uniforme sobre un borde
            recto de longitud ``L``, reparte ``-q_bar·L·t/2`` a cada uno de
            los dos nodos del borde y cero a los demás.

        Notes
        -----
        El signo negativo del retorno viene de la forma débil: el término
        de frontera aparece como ``-∫ w q̄ dΓ``, de modo que un flujo
        saliente positivo **extrae** energía del sistema.

        Un borde sin condición declarada queda con ``q̄ = 0`` — frontera
        adiabática, el default natural de Neumann homogéneo.
        """
        if not (0 <= edge < len(self.EDGE_NODES)):
            raise ValueError(
                f"Quad4Thermal(id={self.id}): edge={edge} fuera de rango; "
                f"esperado 0..{len(self.EDGE_NODES) - 1}."
            )

        a, c = self.EDGE_NODES[edge]
        pa = np.asarray(self.nodes[a].coordinates[:2], dtype=float)
        pc = np.asarray(self.nodes[c].coordinates[:2], dtype=float)
        L = float(np.linalg.norm(pc - pa))

        f = np.zeros(self.N_NODES)
        contrib = -0.5 * L * float(q_bar) * self.thickness
        f[a] = contrib
        f[c] = contrib
        return f
