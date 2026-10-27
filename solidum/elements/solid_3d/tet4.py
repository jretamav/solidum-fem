"""``Tet4`` — tetraedro lineal 3D de deformación constante (CST 3D, ADR 0012)."""
from typing import List

import numpy as np

from solidum.core.element import Element, validate_lumping_kwarg
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.elements.solid_2d._shared import _dN_tri3, _shape_functions_tri3
from solidum.elements.solid_3d._shared import (
    _SolidLoads3D,
    _batch_kin_tet4,
    _compute_integrands_3d,
    _compute_kinematics_tet4,
    _expand_scalar_mass_3d,
    _shape_functions_tet4,
)
from solidum.math.mass_lumping import lump_hrz
from solidum.registry import ElementRegistry


@ElementRegistry.register
class Tet4(_SolidLoads3D, Element):
    """Elemento tetraédrico lineal 3D de deformación constante (CST 3D).

    1 punto de integración central. Análogo 3D del Tri3. Severo shear locking
    y locking volumétrico con ν → 0.5; preferir Hex8 en mallas hexaédricas.
    """

    DOF_NAMES = ['ux', 'uy', 'uz']
    STRAIN_DIM = 6
    N_INTEGRATION_POINTS = 1
    BATCH_KINEMATICS = _batch_kin_tet4
    BATCH_SHAPE_FUNCTIONS = staticmethod(_shape_functions_tet4)
    REFERENCE_CELL = "tet"
    REFERENCE_NODES = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))

    # ADR 0012 — 4 caras con normal saliente. Cara i opuesta al nodo i.
    FACE_NODES = (
        (1, 2, 3),  # 0: opuesta al nodo 0
        (0, 3, 2),  # 1: opuesta al nodo 1
        (0, 1, 3),  # 2: opuesta al nodo 2
        (0, 2, 1),  # 3: opuesta al nodo 3
    )

    def __init__(self, element_id: int, nodes: List[Node], material: Material):
        # Regla de un punto (baricentro, peso = volumen del tetraedro de
        # referencia); materializada para el camino por lotes.
        self.points = [(0.25, 0.25, 0.25)]
        self.weights = [1.0 / 6.0]
        self.quadrature_key = "tet_1"
        super().__init__(element_id, nodes, material)

    # ------------------------------------------------------------------
    # Contrato principal
    # ------------------------------------------------------------------

    def compute_element_state(self, u_e: np.ndarray):
        coords = self.get_coordinate_matrix(ndim=3)
        B, detJ = _compute_kinematics_tet4(coords)
        strain = B @ u_e

        sigma, C_tangent, new_state = self.material.compute_state(strain, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        # 1 punto central, peso 1/6 (= volumen del tet de referencia).
        K_e, F_int_e = _compute_integrands_3d(B, C_tangent, sigma, detJ, 1.0 / 6.0)
        return K_e, F_int_e

    # ------------------------------------------------------------------
    # Salida por punto de Gauss
    # ------------------------------------------------------------------

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """σ y ε del único punto central — coincide con el promedio del elemento."""
        gs = self.compute_gauss_state(U_global)
        return {'stress': gs['stress'][0], 'strain': gs['strain'][0]}

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        """Estado en el punto de Gauss (1 punto baricéntrico) del Tet4.

        ε y σ son uniformes sobre el elemento (CST 3D).
        """
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=3)

        B, _ = _compute_kinematics_tet4(coords)
        strain = B @ u_e
        stress, _, _ = self.material.compute_state(strain, self.state.vars[0])

        centroid = coords.mean(axis=0)

        return {
            'points_natural': np.array([[0.25, 0.25, 0.25]]),
            'points_global': centroid.reshape(1, 3),
            'strain': strain.reshape(1, 6),
            'stress': stress.reshape(1, 6),
        }

    # ------------------------------------------------------------------
    # Cargas distribuidas consistentes
    # ------------------------------------------------------------------

    def _element_volume(self) -> float:
        coords = self.get_coordinate_matrix(ndim=3)
        v1 = coords[1] - coords[0]
        v2 = coords[2] - coords[0]
        v3 = coords[3] - coords[0]
        return abs(float(np.dot(v1, np.cross(v2, v3)))) / 6.0

    # compute_body_load y compute_face_traction: de _SolidLoads3D (ADR 0022).
    # La regla de 4 puntos integra exactamente N_i·b con b lineal (con b
    # constante, b·V_e/4 a cada nodo); las caras, con la de 3 puntos.
    _LOAD_QUADRATURE = "tet_4"
    _FACE_N_FN = staticmethod(_shape_functions_tri3)
    _FACE_DN_FN = staticmethod(_dN_tri3)
    _FACE_QUADRATURE = "tri_3"

    # ------------------------------------------------------------------
    # Masa elemental (ADR 0009)
    # ------------------------------------------------------------------

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa del Tet4 (CST 3D) por fórmula analítica exacta.

        **Consistente** (default): para funciones de forma lineales en
        coordenadas baricéntricas sobre el tetraedro,
        ``∫N_i·N_j dV = V_e·(1+δ_ij)/20``. Resulta el patrón canónico

            M_scalar = (ρ·V_e/20) · (I + 1·1^T)

        donde ``I`` es la identidad 4×4 y ``1`` un vector de unos. Cada
        elemento diagonal vale ``2·ρ·V_e/20`` y cada off-diagonal ``ρ·V_e/20``.

        **Lumped**: HRZ canónico. Para Tet4 coincide con row-sum por simetría
        — ``ρV_e/4`` en cada DOF traslacional.

        Returns
        -------
        np.ndarray, shape (12, 12)
        """
        V = self._element_volume()
        m_total = self.material.density * V
        coef = m_total / 20.0
        M_s = coef * (np.eye(4) + np.ones((4, 4)))
        M_consistent = _expand_scalar_mass_3d(M_s)
        validate_lumping_kwarg(lumping, type(self).__name__)
        if lumping == "lumped":
            return lump_hrz(M_consistent, total_mass=m_total, n_translational_dirs=3)
        return M_consistent
