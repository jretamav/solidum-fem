# fenix_fem/fenix/elements/structural.py
import math
import numpy as np
from typing import List
from fenix.core.element import Element
from fenix.core.node import Node
from fenix.core.material import Material

class Truss2D(Element):
    """
    Elemento estructural de armadura 2D adaptado para análisis no lineal.
    """
    def __init__(self, element_id: int, nodes: List[Node], material: Material, A: float):
        super().__init__(element_id, nodes)
        if len(nodes) != 2:
            raise ValueError("El elemento Truss2D requiere exactamente 2 nodos.")
        
        self.material = material
        self.A = A
        
        # Un punto de integración (una sola variable de estado por barra)
        self.state_vars = [None]
        self.state_vars_trial = [None]
        
        for node in self.nodes:
            node.add_dof('ux')
            node.add_dof('uy')

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'ux', 'uy']

    def commit_state(self):
        """Fija las variables de estado cuando el paso converge."""
        self.state_vars = self.state_vars_trial.copy()

    def compute_element_state(self, u_e: np.ndarray):
        """Calcula matriz tangente y fuerzas internas dadas las deformaciones."""
        x1, y1 = self.nodes[0].coordinates
        x2, y2 = self.nodes[1].coordinates
        dx, dy = x2 - x1, y2 - y1
        L = math.sqrt(dx**2 + dy**2)
        c, s = dx / L, dy / L
        
        # Matriz cinemática B
        B = np.array([-c, -s, c, s]) / L
        
        # Deformación actual eps = B * u_e
        epsilon = np.dot(B, u_e)
        
        # Respuesta constitutiva
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state_vars[0])
        self.state_vars_trial[0] = new_state
        
        # Rigidez tangente local K_e = B^T * E_t * B * A * L
        coef = (E_t * self.A) / L
        K_e = coef * np.array([
            [ c**2,   c*s,  -c**2,  -c*s],
            [ c*s,    s**2, -c*s,   -s**2],
            [-c**2,  -c*s,   c**2,   c*s],
            [-c*s,   -s**2,  c*s,    s**2]
        ])
        
        # Fuerza interna F_int = B^T * sigma * A * L
        F_int_e = (sigma * self.A) * np.array([-c, -s, c, s])
        
        return K_e, F_int_e

    def compute_global_stiffness(self) -> np.ndarray:
        """Retrocompatibilidad para análisis lineal puro."""
        K_e, _ = self.compute_element_state(np.zeros(4))
        return K_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Utilidad para reportes de post-procesamiento."""
        u_e = np.array([
            U_global[self.nodes[0].dofs['ux']],
            U_global[self.nodes[0].dofs['uy']],
            U_global[self.nodes[1].dofs['ux']],
            U_global[self.nodes[1].dofs['uy']]
        ])
        
        x1, y1 = self.nodes[0].coordinates
        x2, y2 = self.nodes[1].coordinates
        L = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        c, s = (x2 - x1) / L, (y2 - y1) / L
        B = np.array([-c, -s, c, s]) / L
        
        epsilon = np.dot(B, u_e)
        sigma, _, _ = self.material.compute_state(epsilon, self.state_vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}
