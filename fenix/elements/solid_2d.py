# fenix_fem/fenix/elements/solid_2d.py
import numpy as np
from typing import List
from fenix.core.element import Element
from fenix.core.node import Node
from fenix.core.material import Material
from fenix.registry import ElementRegistry, QuadratureRegistry
from fenix.constants import ZERO_JACOBIAN_TOL
from numba import njit


@njit
def _compute_kinematics(xi, eta, coords):
    # Preasignar memoria en Numba es la forma 100% segura
    dN_dxi = np.zeros((2, 4), dtype=np.float64)
    dN_dxi[0, 0] = -(1.0 - eta) / 4.0
    dN_dxi[0, 1] =  (1.0 - eta) / 4.0
    dN_dxi[0, 2] =  (1.0 + eta) / 4.0
    dN_dxi[0, 3] = -(1.0 + eta) / 4.0
    dN_dxi[1, 0] = -(1.0 - xi) / 4.0
    dN_dxi[1, 1] = -(1.0 + xi) / 4.0
    dN_dxi[1, 2] =  (1.0 + xi) / 4.0
    dN_dxi[1, 3] =  (1.0 - xi) / 4.0

    # Jacobiano
    J = np.dot(dN_dxi, coords)
    detJ = J[0,0] * J[1,1] - J[0,1] * J[1,0]

    if detJ <= ZERO_JACOBIAN_TOL:
        raise ValueError("Jacobiano negativo o cero detectado en elemento Quad4. Revisa la conectividad o distorsion.")

    # Inversa del Jacobiano segura
    invJ = np.zeros((2, 2), dtype=np.float64)
    invJ[0, 0] =  J[1, 1] / detJ
    invJ[0, 1] = -J[0, 1] / detJ
    invJ[1, 0] = -J[1, 0] / detJ
    invJ[1, 1] =  J[0, 0] / detJ

    # Derivadas respecto a coordenadas globales (x, y)
    dN_dx = np.dot(invJ, dN_dxi)

    # Ensamblaje de la matriz B (Deformación-Desplazamiento)
    B = np.zeros((3, 8), dtype=np.float64)
    for i in range(4):
        B[0, 2*i]     = dN_dx[0, i]
        B[1, 2*i + 1] = dN_dx[1, i]
        B[2, 2*i]     = dN_dx[1, i]
        B[2, 2*i + 1] = dN_dx[0, i]

    return B, detJ


@njit
def _compute_integrands(B, C_alg, sigma, detJ, weight, thickness):
    dV = detJ * weight * thickness

    # Descomponer el triple producto tensorial evita errores de Numba
    temp = np.dot(B.T, C_alg)
    K_contrib = np.dot(temp, B) * dV

    F_contrib = np.dot(B.T, sigma) * dV
    return K_contrib, F_contrib


@njit
def _compute_kinematics_tri3(coords):
    # Derivadas de funciones de forma analíticas para Tri3 (xi, eta)
    dN_dxi = np.zeros((2, 3), dtype=np.float64)
    dN_dxi[0, 0] = -1.0; dN_dxi[0, 1] = 1.0; dN_dxi[0, 2] = 0.0
    dN_dxi[1, 0] = -1.0; dN_dxi[1, 1] = 0.0; dN_dxi[1, 2] = 1.0

    J = np.dot(dN_dxi, coords)
    detJ = J[0,0] * J[1,1] - J[0,1] * J[1,0]

    if detJ <= ZERO_JACOBIAN_TOL:
        raise ValueError("Jacobiano negativo o cero detectado en elemento Tri3.")

    invJ = np.zeros((2, 2), dtype=np.float64)
    invJ[0, 0] =  J[1, 1] / detJ
    invJ[0, 1] = -J[0, 1] / detJ
    invJ[1, 0] = -J[1, 0] / detJ
    invJ[1, 1] =  J[0, 0] / detJ

    dN_dx = np.dot(invJ, dN_dxi)

    B = np.zeros((3, 6), dtype=np.float64)
    for i in range(3):
        B[0, 2*i]     = dN_dx[0, i]
        B[1, 2*i + 1] = dN_dx[1, i]
        B[2, 2*i]     = dN_dx[1, i]
        B[2, 2*i + 1] = dN_dx[0, i]

    return B, detJ


@ElementRegistry.register
class Quad4(Element):
    """
    Elemento cuadrilátero bilineal 2D isoparamétrico con integración de Gauss.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        4 nodos en sentido antihorario.
    material : Material
        Material 2D (STRAIN_DIM=3).
    thickness : float, optional
        Espesor para estado plano. Default 1.0.
    quadrature : tuple, optional
        Tupla (puntos, pesos). Si es None, usa Gauss 2×2.

    Notes
    -----
    Integración reducida ("1×1") es propensa a hourglassing.
    """
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 3
    N_INTEGRATION_POINTS = 4   # default Gauss 2×2; el constructor lo sobreescribe a nivel de instancia si se pasa otra cuadratura.

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 thickness: float = 1.0, quadrature: tuple = None):
        if quadrature is None:
            self.points, self.weights = QuadratureRegistry.get("2x2")
        else:
            self.points, self.weights = quadrature

        self.thickness = thickness
        # ClassVar override: la base lee este atributo para dimensionar el ElementState
        self.N_INTEGRATION_POINTS = len(self.points)

        super().__init__(element_id, nodes, material)

    def compute_element_state(self, u_e: np.ndarray):
        """Calcula K_tangente local y Fuerzas Internas dadas las deformaciones actuales."""
        K_e = np.zeros((8, 8))
        F_int_e = np.zeros(8)

        coords = self.get_coordinate_matrix(ndim=2)

        for idx, (p, w) in enumerate(zip(self.points, self.weights)):
            xi, eta = p
            B, detJ = _compute_kinematics(xi, eta, coords)
            strain = B @ u_e

            sigma, C_tangent, new_state = self.material.compute_state(strain, self.state.vars[idx])
            self.state.vars_trial[idx] = new_state
            self.state.stresses_trial[idx] = sigma

            K_contrib, F_contrib = _compute_integrands(B, C_tangent, sigma, detJ, w, self.thickness)
            K_e += K_contrib
            F_int_e += F_contrib

        return K_e, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Utilidad para reportes de post-procesamiento."""
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=2)

        avg_stress = np.zeros(3)
        avg_strain = np.zeros(3)
        for idx, p in enumerate(self.points):
            xi, eta = p
            B, _ = _compute_kinematics(xi, eta, coords)
            strain = B @ u_e
            sigma, _, _ = self.material.compute_state(strain, self.state.vars[idx])
            avg_strain += strain
            avg_stress += sigma

        if self.N_INTEGRATION_POINTS == 0:
            return {'stress': np.zeros(3), 'strain': np.zeros(3)}
        n = self.N_INTEGRATION_POINTS
        return {'stress': avg_stress / n, 'strain': avg_strain / n}


@ElementRegistry.register
class Tri3(Element):
    """
    Elemento triangular lineal 2D de deformación constante (CST).

    1 punto de integración central. Propenso a shear locking; preferir Quad4
    salvo en transiciones de malla complejas.
    """
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 3
    N_INTEGRATION_POINTS = 1

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 thickness: float = 1.0):
        self.thickness = thickness
        super().__init__(element_id, nodes, material)

    def compute_element_state(self, u_e: np.ndarray):
        coords = self.get_coordinate_matrix(ndim=2)

        B, detJ = _compute_kinematics_tri3(coords)
        strain = B @ u_e

        sigma, C_tangent, new_state = self.material.compute_state(strain, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        # 1 punto central, peso 0.5 (área del triángulo en coordenadas naturales)
        K_e, F_int_e = _compute_integrands(B, C_tangent, sigma, detJ, 0.5, self.thickness)
        return K_e, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=2)

        B, _ = _compute_kinematics_tri3(coords)
        strain = B @ u_e
        sigma, _, _ = self.material.compute_state(strain, self.state.vars[0])
        return {'stress': sigma, 'strain': strain}
