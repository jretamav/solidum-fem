"""``Tet4`` — tetraedro lineal 3D de deformación constante (CST 3D, ADR 0012)."""
from typing import List

import numpy as np

from fenix.core.element import Element, validate_lumping_kwarg
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.elements.solid_3d._shared import (
    _compute_integrands_3d,
    _compute_kinematics_tet4,
    _expand_scalar_mass_3d,
)
from fenix.math.mass_lumping import lump_hrz
from fenix.registry import ElementRegistry


@ElementRegistry.register
class Tet4(Element):
    """Elemento tetraédrico lineal 3D de deformación constante (CST 3D).

    1 punto de integración central. Análogo 3D del Tri3. Severo shear locking
    y locking volumétrico con ν → 0.5; preferir Hex8 en mallas hexaédricas.
    """

    DOF_NAMES = ['ux', 'uy', 'uz']
    STRAIN_DIM = 6
    N_INTEGRATION_POINTS = 1

    # ADR 0012 — 4 caras con normal saliente. Cara i opuesta al nodo i.
    FACE_NODES = (
        (1, 2, 3),  # 0: opuesta al nodo 0
        (0, 3, 2),  # 1: opuesta al nodo 1
        (0, 1, 3),  # 2: opuesta al nodo 2
        (0, 2, 1),  # 3: opuesta al nodo 3
    )

    def __init__(self, element_id: int, nodes: List[Node], material: Material):
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

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector de cargas nodales consistente con una fuerza de cuerpo uniforme.

        Para Tet4 con N lineal y b uniforme la integral es exacta:
        cada nodo recibe ``b · V_e / 4``.

        Returns
        -------
        ndarray (12,)
            Cargas equivalentes en los DOFs locales.
        """
        b = np.asarray(b, dtype=np.float64).reshape(3)
        V = self._element_volume()
        share = V / 4.0
        f = np.zeros(12)
        for i in range(4):
            f[3 * i]     = share * b[0]
            f[3 * i + 1] = share * b[1]
            f[3 * i + 2] = share * b[2]
        return f

    def compute_face_traction(self, face: int, t_vec: np.ndarray) -> np.ndarray:
        """Vector de cargas nodales consistente con una tracción uniforme en una cara.

        Para Tet4 con N lineal sobre la cara triangular y t̄ uniforme el
        reparto es exacto: ``A_cara · t̄ / 3`` a cada nodo del triángulo,
        cero al nodo opuesto.

        Parameters
        ----------
        face : int
            Índice de cara según ``FACE_NODES`` (0..3). Cara i opuesta al nodo i.
        t_vec : ndarray (3,)
            Tracción uniforme en globales sobre la cara.

        Returns
        -------
        ndarray (12,)
            Cargas equivalentes en los DOFs locales.
        """
        if face not in range(4):
            raise ValueError(f"face={face} fuera de rango para Tet4 (0..3).")
        t_vec = np.asarray(t_vec, dtype=np.float64).reshape(3)

        face_local = self.FACE_NODES[face]
        X = np.array(
            [self.nodes[i].coordinates[:3] for i in face_local],
            dtype=np.float64,
        )  # 3×3
        v1 = X[1] - X[0]
        v2 = X[2] - X[0]
        A = 0.5 * float(np.linalg.norm(np.cross(v1, v2)))

        share = A / 3.0
        f = np.zeros(12)
        for node_local in face_local:
            base = 3 * node_local
            f[base]     = share * t_vec[0]
            f[base + 1] = share * t_vec[1]
            f[base + 2] = share * t_vec[2]
        return f

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
