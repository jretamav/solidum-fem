# fenix_fem/fenix/elements/solid_2d.py
import numpy as np
from typing import List
from fenix.core.element import Element
from fenix.core.node import Node
from fenix.core.material import Material
from fenix.math.integration import GaussQuadrature
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
    
    if detJ <= 0.0:
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
    
    if detJ <= 0.0:
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

class Quad4(Element):
    def __init__(self, element_id: int, nodes: List[Node], material: Material, thickness: float = 1.0, quadrature: str = "2x2"):
        super().__init__(element_id, nodes)
        self.material = material
        self.thickness = thickness
        self.quadrature = quadrature if quadrature else "2x2"
        
        # Puntos de integración
        if self.quadrature == "2x2":
            self.points, self.weights = GaussQuadrature.get_points_2d_2x2()
        else:
            self.points, self.weights = GaussQuadrature.get_points_2d_2x2() # Fallback por seguridad
            
        self.num_ip = len(self.points)
        # Variables internas dependientes dinámicamente de la regla de integración
        self.state_vars = [None] * self.num_ip
        self.state_vars_trial = [None] * self.num_ip
        self.stresses = [np.zeros(3) for _ in range(self.num_ip)]
        self.stresses_trial = [np.zeros(3) for _ in range(self.num_ip)]
        
        for node in self.nodes:
            node.add_dof('ux')
            node.add_dof('uy')

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'ux', 'uy', 'ux', 'uy', 'ux', 'uy']

    def commit_state(self):
        """Guarda permanentemente las variables internas una vez que el paso de carga converge."""
        self.state_vars = self.state_vars_trial.copy()
        self.stresses = self.stresses_trial.copy()

    def compute_element_state(self, u_e: np.ndarray):
        """Calcula K_tangente local y Fuerzas Internas dadas las deformaciones actuales."""
        K_e = np.zeros((8, 8))
        F_int_e = np.zeros(8)
        
        coords = self.get_coordinate_matrix(ndim=2)
        
        for idx, (p, w) in enumerate(zip(self.points, self.weights)):
            xi, eta = p
            B, detJ = _compute_kinematics(xi, eta, coords)
            strain = B @ u_e
            
            # Llamada constitutiva con memoria histórica
            sigma, C_tangent, new_state = self.material.compute_state(strain, self.state_vars[idx])
            self.state_vars_trial[idx] = new_state
            self.stresses_trial[idx] = sigma
            
            K_contrib, F_contrib = _compute_integrands(B, C_tangent, sigma, detJ, w, self.thickness)
            K_e += K_contrib
            F_int_e += F_contrib
            
        return K_e, F_int_e

    def compute_global_stiffness(self) -> np.ndarray:
        K_e, _ = self.compute_element_state(np.zeros(8))
        return K_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Utilidad para reportes de post-procesamiento."""
        u_e = np.zeros(8)
        for i, node in enumerate(self.nodes):
            u_e[2*i] = U_global[node.dofs['ux']]
            u_e[2*i+1] = U_global[node.dofs['uy']]
            
        coords = self.get_coordinate_matrix(ndim=2)
        
        avg_stress = np.zeros(3)
        avg_strain = np.zeros(3)
        for idx, p in enumerate(self.points):
            xi, eta = p
            B, _ = _compute_kinematics(xi, eta, coords)
            strain = B @ u_e
            sigma, _, _ = self.material.compute_state(strain, self.state_vars[idx])
            avg_strain += strain
            avg_stress += sigma
            
        if self.num_ip == 0: return {'stress': np.zeros(3), 'strain': np.zeros(3)}
        return {'stress': avg_stress / self.num_ip, 'strain': avg_strain / self.num_ip}


class Tri3(Element):
    def __init__(self, element_id: int, nodes: List[Node], material: Material, thickness: float = 1.0):
        super().__init__(element_id, nodes)
        self.material = material
        self.thickness = thickness
        
        self.state_vars = [None]
        self.state_vars_trial = [None]
        self.stresses = [np.zeros(3)]
        self.stresses_trial = [np.zeros(3)]
        
        for node in self.nodes:
            node.add_dof('ux')
            node.add_dof('uy')

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'ux', 'uy', 'ux', 'uy']

    def commit_state(self):
        self.state_vars = self.state_vars_trial.copy()
        self.stresses = self.stresses_trial.copy()

    def compute_element_state(self, u_e: np.ndarray):
        coords = self.get_coordinate_matrix(ndim=2)
        
        B, detJ = _compute_kinematics_tri3(coords)
        strain = B @ u_e
        
        sigma, C_tangent, new_state = self.material.compute_state(strain, self.state_vars[0])
        self.state_vars_trial[0] = new_state
        self.stresses_trial[0] = sigma
        
        # 1 punto de integración central con peso 0.5 (Área del triángulo en coord. naturales)
        K_e, F_int_e = _compute_integrands(B, C_tangent, sigma, detJ, 0.5, self.thickness)
        return K_e, F_int_e

    def compute_global_stiffness(self) -> np.ndarray:
        K_e, _ = self.compute_element_state(np.zeros(6))
        return K_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Utilidad para reportes de post-procesamiento."""
        u_e = np.zeros(6)
        for i, node in enumerate(self.nodes):
            u_e[2*i] = U_global[node.dofs['ux']]
            u_e[2*i+1] = U_global[node.dofs['uy']]
            
        coords = self.get_coordinate_matrix(ndim=2)
        
        B, _ = _compute_kinematics_tri3(coords)
        strain = B @ u_e
        sigma, _, _ = self.material.compute_state(strain, self.state_vars[0])
        return {'stress': sigma, 'strain': strain}
