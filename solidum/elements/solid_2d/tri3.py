"""``Tri3`` — triángulo lineal 2D de deformación constante (CST)."""
from typing import List

import numpy as np

from solidum.core.element import Element, validate_lumping_kwarg
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.elements.solid_2d._shared import (
    _PlaneSolidLoads,
    _batch_kin_tri3,
    compute_integrands,
    compute_kinematics_tri3,
    _expand_scalar_mass,
    _shape_functions_tri3,
)
from solidum.math.mass_lumping import lump_hrz
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Tri3(_PlaneSolidLoads, Element):
    """
    Elemento triangular lineal 2D de deformación constante (CST).

    1 punto de integración central. Propenso a shear locking; preferir Quad4
    salvo en transiciones de malla complejas.
    """
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 3
    N_INTEGRATION_POINTS = 1
    BATCH_KINEMATICS = _batch_kin_tri3
    BATCH_SCALE = "thickness"
    BATCH_SHAPE_FUNCTIONS = staticmethod(_shape_functions_tri3)
    REFERENCE_CELL = "tri"
    REFERENCE_NODES = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 thickness: float = 1.0):
        self.thickness = thickness
        # Regla de un punto (centroide, peso = área del triángulo de
        # referencia). El elemento la usa implícitamente; se materializa para
        # que el camino por lotes la lea como en el resto de sólidos.
        self.points = [(1.0 / 3.0, 1.0 / 3.0)]
        self.weights = [0.5]
        self.quadrature_key = "tri_1"
        super().__init__(element_id, nodes, material)

    def compute_element_state(self, u_e: np.ndarray):
        coords = self.get_coordinate_matrix(ndim=2)

        B, detJ = compute_kinematics_tri3(coords)
        strain = B @ u_e

        sigma, C_tangent, new_state = self.material.compute_state(strain, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        # 1 punto central, peso 0.5 (área del triángulo en coordenadas naturales)
        K_e, F_int_e = compute_integrands(B, C_tangent, sigma, detJ, 0.5, self.thickness)
        return K_e, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """σ y ε del único punto central — coincide con el promedio del elemento."""
        gs = self.compute_gauss_state(U_global)
        return {'stress': gs['stress'][0], 'strain': gs['strain'][0]}

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        """Estado en el punto de Gauss (1 punto central) del Tri3.

        Devuelve la misma estructura que ``Quad4.compute_gauss_state``;
        para Tri3 ε y σ son uniformes sobre el elemento (CST).
        """
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=2)

        B, _ = compute_kinematics_tri3(coords)
        strain = B @ u_e
        stress, _, _ = self.material.compute_state(strain, self.state.vars[0])

        # Punto central en coordenadas baricéntricas (1/3, 1/3): el centroide.
        centroid = coords.mean(axis=0)

        return {
            'points_natural': np.array([[1.0 / 3.0, 1.0 / 3.0]]),
            'points_global': centroid.reshape(1, 2),
            'strain': strain.reshape(1, 3),
            'stress': stress.reshape(1, 3),
        }

    # ------------------------------------------------------------------
    # Cargas distribuidas consistentes: compute_body_load y
    # compute_edge_traction de _PlaneSolidLoads (ADR 0022). La regla de
    # 3 puntos integra exactamente N_i·b con b lineal; la del elemento (el
    # centroide) sólo con b constante.
    # ------------------------------------------------------------------

    EDGE_NODES = ((0, 1), (1, 2), (2, 0))
    _LOAD_QUADRATURE = "tri_3"

    def _element_area(self) -> float:
        coords = self.get_coordinate_matrix(ndim=2)
        x = coords[:, 0]; y = coords[:, 1]
        return 0.5 * abs(
            (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])
        )

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa del Tri3 (CST) por fórmula analítica exacta (ADR 0009).

        **Consistente** (default): para funciones de forma lineales en
        coordenadas baricéntricas, ``∫N_i·N_j dA = A·(1+δ_ij)/12``. Resulta

            M_scalar = (ρ·t·A/12) · [[2, 1, 1],
                                     [1, 2, 1],
                                     [1, 1, 2]]

        La cuadratura del propio elemento (1 punto central) sería
        insuficiente (integraría exactamente solo polinomios lineales);
        por eso se usa la fórmula cerrada.

        **Lumped** (fase 2): HRZ canónico. Para Tri3 coincide con row-sum
        por simetría — ``ρAt/3`` en cada DOF traslacional (masa total ``ρAt``
        repartida equitativamente entre los 3 nodos).

        Returns
        -------
        np.ndarray, shape (6, 6)
        """
        m_total = self.material.density * self._element_area() * self.thickness
        coef = m_total / 12.0
        M_s = coef * np.array([[2.0, 1.0, 1.0],
                               [1.0, 2.0, 1.0],
                               [1.0, 1.0, 2.0]])
        M_consistent = _expand_scalar_mass(M_s)
        validate_lumping_kwarg(lumping, type(self).__name__)
        if lumping == "lumped":
            return lump_hrz(M_consistent, total_mass=m_total, n_translational_dirs=2)
        return M_consistent
