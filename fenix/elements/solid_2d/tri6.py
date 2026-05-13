"""``Tri6`` — triángulo 2D de 6 nodos (cuadrático completo P₂)."""
from typing import List

import numpy as np

from fenix.core.element import Element
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.elements.solid_2d._shared import (
    _compute_integrands,
    _dN_tri6,
    _expand_scalar_mass,
    _kinematics_higher_order,
    _N_tri6,
    _quadratic_edge_traction,
)
from fenix.registry import ElementRegistry, QuadratureRegistry


@ElementRegistry.register
class Tri6(Element):
    """Triángulo 2D de 6 nodos (cuadrático completo P₂).

    Cura el shear locking severo del Tri3; reproduce campos cuadráticos
    exactamente. Cuadratura `tri_3` (3 puntos en los puntos medios).
    """
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 3
    N_INTEGRATION_POINTS = 3

    EDGE_NODES = (
        (0, 3, 1),
        (1, 4, 2),
        (2, 5, 0),
    )

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 thickness: float = 1.0, quadrature: tuple = None):
        if quadrature is None:
            self.points, self.weights = QuadratureRegistry.get("tri_3")
        else:
            self.points, self.weights = quadrature
        self.thickness = thickness
        self.N_INTEGRATION_POINTS = len(self.points)
        super().__init__(element_id, nodes, material)

    def compute_element_state(self, u_e: np.ndarray):
        K_e = np.zeros((12, 12))
        F_int_e = np.zeros(12)
        coords = self.get_coordinate_matrix(ndim=2)
        for idx, ((xi, eta), w) in enumerate(zip(self.points, self.weights)):
            B, detJ = _kinematics_higher_order(_dN_tri6, xi, eta, coords, 6)
            strain = B @ u_e
            sigma, C, new_state = self.material.compute_state(strain, self.state.vars[idx])
            self.state.vars_trial[idx] = new_state
            self.state.stresses_trial[idx] = sigma
            K_contrib, F_contrib = _compute_integrands(B, C, sigma, detJ, w, self.thickness)
            K_e += K_contrib
            F_int_e += F_contrib
        return K_e, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        gs = self.compute_gauss_state(U_global)
        return {'stress': gs['stress'].mean(axis=0), 'strain': gs['strain'].mean(axis=0)}

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=2)
        n_g = self.N_INTEGRATION_POINTS
        nat = np.asarray(self.points, dtype=np.float64).reshape(n_g, 2)
        glb = np.zeros((n_g, 2))
        eps = np.zeros((n_g, 3))
        sig = np.zeros((n_g, 3))
        for idx, (xi, eta) in enumerate(self.points):
            N = _N_tri6(xi, eta)
            glb[idx] = N @ coords
            B, _ = _kinematics_higher_order(_dN_tri6, xi, eta, coords, 6)
            strain = B @ u_e
            stress, _, _ = self.material.compute_state(strain, self.state.vars[idx])
            eps[idx] = strain
            sig[idx] = stress
        return {'points_natural': nat, 'points_global': glb,
                'strain': eps, 'stress': sig}

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        b = np.asarray(b, dtype=np.float64).reshape(2)
        coords = self.get_coordinate_matrix(ndim=2)
        f = np.zeros(12)
        for (xi, eta), w in zip(self.points, self.weights):
            N = _N_tri6(xi, eta)
            _, detJ = _kinematics_higher_order(_dN_tri6, xi, eta, coords, 6)
            factor = detJ * w * self.thickness
            for i in range(6):
                f[2 * i]     += N[i] * b[0] * factor
                f[2 * i + 1] += N[i] * b[1] * factor
        return f

    def compute_edge_traction(self, edge: int, t_vec: np.ndarray) -> np.ndarray:
        if edge not in (0, 1, 2):
            raise ValueError(f"edge={edge} fuera de rango para Tri6 (0..2).")
        return _quadratic_edge_traction(self, edge, t_vec, n_dofs=12)

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa consistente del Tri6 con cuadratura específica ``tri_6``
        (Dunavant 6 puntos, orden 4) — independiente de la cuadratura del
        elemento (ADR 0009).

        La cuadratura ``tri_3`` del elemento integra exactamente la rigidez
        (producto lineal×lineal = orden 2) pero **subintegra la masa**
        (producto cuadrático×cuadrático = orden 4). La subintegración deja
        modos nulos espurios en M (semi-definida positiva en lugar de
        definida), que rompen el problema generalizado ``K·φ = ω²·M·φ``.
        Usar una cuadratura específica para masa es práctica estándar
        cuando rigidez y masa requieren órdenes distintos.

        Returns
        -------
        np.ndarray, shape (12, 12)
        """
        if lumping != "consistent":
            raise NotImplementedError(
                f"Tri6.compute_mass_matrix: lumping='{lumping}' no "
                f"implementado. Fase 1 (ADR 0009) solo admite 'consistent'."
            )
        rho = self.material.density
        coords = self.get_coordinate_matrix(ndim=2)
        # Cuadratura específica para masa (orden 4 exacto), no la del elemento.
        mass_points, mass_weights = QuadratureRegistry.get("tri_6")
        M_s = np.zeros((6, 6))
        for (xi, eta), w in zip(mass_points, mass_weights):
            N = _N_tri6(xi, eta)
            _, detJ = _kinematics_higher_order(_dN_tri6, xi, eta, coords, 6)
            M_s += rho * np.outer(N, N) * (detJ * w * self.thickness)
        return _expand_scalar_mass(M_s)
