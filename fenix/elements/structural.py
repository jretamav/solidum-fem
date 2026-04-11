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

        # Optimización: calcular geometría inicial una sola vez
        coords1 = self.nodes[0].coordinates
        coords2 = self.nodes[1].coordinates
        x1, y1 = coords1[0], coords1[1]
        x2, y2 = coords2[0], coords2[1]
        
        dx, dy = x2 - x1, y2 - y1
        self.L0 = math.sqrt(dx**2 + dy**2)
        if self.L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")
        self.c = dx / self.L0
        self.s = dy / self.L0

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'ux', 'uy']

    def commit_state(self):
        """Fija las variables de estado cuando el paso converge."""
        self.state_vars = self.state_vars_trial.copy()

    def compute_element_state(self, u_e: np.ndarray):
        """Calcula matriz tangente y fuerzas internas dadas las deformaciones."""
        c, s, L = self.c, self.s, self.L0
        
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
        
        c, s, L = self.c, self.s, self.L0
        B = np.array([-c, -s, c, s]) / L
        
        epsilon = np.dot(B, u_e)
        sigma, _, _ = self.material.compute_state(epsilon, self.state_vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}

class Truss3D(Element):
    """
    Elemento estructural de armadura 3D adaptado para análisis no lineal.
    """
    def __init__(self, element_id: int, nodes: List[Node], material: Material, A: float):
        super().__init__(element_id, nodes)
        if len(nodes) != 2:
            raise ValueError("El elemento Truss3D requiere exactamente 2 nodos.")
        
        self.material = material
        self.A = A
        
        # Un punto de integración (una sola variable de estado por barra)
        self.state_vars = [None]
        self.state_vars_trial = [None]
        
        for node in self.nodes:
            node.add_dof('ux')
            node.add_dof('uy')
            node.add_dof('uz')

        # Optimización: calcular geometría inicial una vez
        coords1 = self.nodes[0].coordinates
        coords2 = self.nodes[1].coordinates
        
        # Soporte para nodos que puedan tener 2 o 3 coordenadas en su array por error/defecto
        x1, y1 = coords1[0], coords1[1]
        z1 = coords1[2] if len(coords1) > 2 else 0.0
        x2, y2 = coords2[0], coords2[1]
        z2 = coords2[2] if len(coords2) > 2 else 0.0

        dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
        self.L0 = math.sqrt(dx**2 + dy**2 + dz**2)
        
        if self.L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")
            
        self.cx = dx / self.L0
        self.cy = dy / self.L0
        self.cz = dz / self.L0

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'uz', 'ux', 'uy', 'uz']

    def commit_state(self):
        """Fija las variables de estado cuando el paso converge."""
        self.state_vars = self.state_vars_trial.copy()

    def compute_element_state(self, u_e: np.ndarray):
        """Calcula matriz tangente y fuerzas internas dadas las deformaciones."""
        cx, cy, cz, L = self.cx, self.cy, self.cz, self.L0
        
        # Vector dirección
        B_dir = np.array([-cx, -cy, -cz, cx, cy, cz])
        
        # Matriz cinemática B
        B = B_dir / L
        
        # Deformación actual eps = B * u_e
        epsilon = np.dot(B, u_e)
        
        # Respuesta constitutiva
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state_vars[0])
        self.state_vars_trial[0] = new_state
        
        # Rigidez tangente local K_e = B^T * E_t * B * A * L
        coef = (E_t * self.A) / L
        K_e = coef * np.outer(B_dir, B_dir)
        
        # Fuerza interna F_int = B^T * sigma * A * L
        F_int_e = (sigma * self.A) * B_dir
        
        return K_e, F_int_e

    def compute_global_stiffness(self) -> np.ndarray:
        """Retrocompatibilidad para análisis lineal puro."""
        K_e, _ = self.compute_element_state(np.zeros(6))
        return K_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Utilidad para reportes de post-procesamiento."""
        u_e = np.array([
            U_global[self.nodes[0].dofs['ux']],
            U_global[self.nodes[0].dofs['uy']],
            U_global[self.nodes[0].dofs['uz']],
            U_global[self.nodes[1].dofs['ux']],
            U_global[self.nodes[1].dofs['uy']],
            U_global[self.nodes[1].dofs['uz']]
        ])
        
        cx, cy, cz, L = self.cx, self.cy, self.cz, self.L0
        B = np.array([-cx, -cy, -cz, cx, cy, cz]) / L
        
        epsilon = np.dot(B, u_e)
        sigma, _, _ = self.material.compute_state(epsilon, self.state_vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}

class Frame2DEuler(Element):
    """
    Elemento de viga 2D basado en la teoría de Euler-Bernoulli (vigas esbeltas).
    Desprecia la deformación por cortante transversal.
    """
    def __init__(self, element_id: int, nodes: List[Node], material: Material, A: float, I: float):
        super().__init__(element_id, nodes)
        if len(nodes) != 2:
            raise ValueError("El elemento Frame2DEuler requiere exactamente 2 nodos.")
        
        self.material = material
        self.A = A
        self.I = I
        
        for node in self.nodes:
            node.add_dof('ux')
            node.add_dof('uy')
            node.add_dof('rz')  # Grado de libertad rotacional

        coords1 = self.nodes[0].coordinates
        coords2 = self.nodes[1].coordinates
        x1, y1 = coords1[0], coords1[1]
        x2, y2 = coords2[0], coords2[1]
        
        dx, dy = x2 - x1, y2 - y1
        self.L0 = math.sqrt(dx**2 + dy**2)
        if self.L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")
            
        c = dx / self.L0
        s = dy / self.L0
        
        # Matriz de transformación T (Local a Global 6x6)
        self.T = np.array([
            [ c,  s,  0,  0,  0,  0],
            [-s,  c,  0,  0,  0,  0],
            [ 0,  0,  1,  0,  0,  0],
            [ 0,  0,  0,  c,  s,  0],
            [ 0,  0,  0, -s,  c,  0],
            [ 0,  0,  0,  0,  0,  1]
        ])

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'rz', 'ux', 'uy', 'rz']

    def compute_element_state(self, u_e: np.ndarray):
        # Extracción segura del módulo elástico, compatible con la API de Fenix
        if hasattr(self.material, 'E'):
            E = self.material.E
        else:
            # Evaluar en deformación nula para obtener el módulo tangente inicial
            _, E, _ = self.material.compute_state(0.0)
            
        L = self.L0
        EA_L = E * self.A / L
        EI_L = E * self.I / L
        EI_L2 = EI_L / L
        EI_L3 = EI_L2 / L
        
        # Matriz de Rigidez Local K'
        K_local = np.array([
            [ EA_L,        0,        0, -EA_L,        0,        0],
            [    0, 12*EI_L3,  6*EI_L2,     0,-12*EI_L3,  6*EI_L2],
            [    0,  6*EI_L2,   4*EI_L,     0, -6*EI_L2,   2*EI_L],
            [-EA_L,        0,        0,  EA_L,        0,        0],
            [    0,-12*EI_L3, -6*EI_L2,     0, 12*EI_L3, -6*EI_L2],
            [    0,  6*EI_L2,   2*EI_L,     0, -6*EI_L2,   4*EI_L]
        ])
        
        K_global = self.T.T @ K_local @ self.T
        F_int_e = K_global @ u_e
        
        return K_global, F_int_e

    def compute_global_stiffness(self) -> np.ndarray:
        K_e, _ = self.compute_element_state(np.zeros(6))
        return K_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = np.array([
            U_global[self.nodes[0].dofs['ux']], U_global[self.nodes[0].dofs['uy']], U_global[self.nodes[0].dofs['rz']],
            U_global[self.nodes[1].dofs['ux']], U_global[self.nodes[1].dofs['uy']], U_global[self.nodes[1].dofs['rz']]
        ])
        K_global, F_int = self.compute_element_state(u_e)
        F_local = self.T @ F_int
        return {
            'axial_force': F_local[0], 
            'shear_force': F_local[1], 
            'moment_i': F_local[2], 
            'moment_j': F_local[5]
        }

class Frame2DTimoshenko(Element):
    """
    Elemento de viga 2D basado en la teoría de Timoshenko (vigas gruesas/cortas).
    Incluye la deformación por cortante en su formulación matricial.
    """
    def __init__(self, element_id: int, nodes: List[Node], material: Material, A: float, I: float, As: float, nu: float = 0.3):
        super().__init__(element_id, nodes)
        if len(nodes) != 2:
            raise ValueError("El elemento Frame2DTimoshenko requiere exactamente 2 nodos.")
        
        self.material = material
        self.A = A
        self.I = I
        self.As = As
        # Obtenemos Poisson del material o utilizamos el provisto por el usuario
        self.nu = getattr(material, 'nu', nu)
        
        for node in self.nodes:
            node.add_dof('ux')
            node.add_dof('uy')
            node.add_dof('rz')

        coords1 = self.nodes[0].coordinates
        coords2 = self.nodes[1].coordinates
        x1, y1 = coords1[0], coords1[1]
        x2, y2 = coords2[0], coords2[1]
        
        dx, dy = x2 - x1, y2 - y1
        self.L0 = math.sqrt(dx**2 + dy**2)
        if self.L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")
            
        c = dx / self.L0
        s = dy / self.L0
        
        self.T = np.array([
            [ c,  s,  0,  0,  0,  0],
            [-s,  c,  0,  0,  0,  0],
            [ 0,  0,  1,  0,  0,  0],
            [ 0,  0,  0,  c,  s,  0],
            [ 0,  0,  0, -s,  c,  0],
            [ 0,  0,  0,  0,  0,  1]
        ])

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'rz', 'ux', 'uy', 'rz']

    def compute_element_state(self, u_e: np.ndarray):
        if hasattr(self.material, 'E'):
            E = self.material.E
        else:
            _, E, _ = self.material.compute_state(0.0)
            
        G = E / (2.0 * (1.0 + self.nu))
        L = self.L0
        
        # Factor de rigidez al cortante de Timoshenko
        Phi = (12.0 * E * self.I) / (G * self.As * (L**2))
        
        EA_L = E * self.A / L
        EI_L = E * self.I / L
        
        a = 12 * EI_L / (L**2 * (1 + Phi))
        b = 6 * EI_L / (L * (1 + Phi))
        c = (4 + Phi) * EI_L / (1 + Phi)
        d = (2 - Phi) * EI_L / (1 + Phi)
        
        K_local = np.array([
            [ EA_L,  0,  0, -EA_L,  0,  0],
            [    0,  a,  b,     0, -a,  b],
            [    0,  b,  c,     0, -b,  d],
            [-EA_L,  0,  0,  EA_L,  0,  0],
            [    0, -a, -b,     0,  a, -b],
            [    0,  b,  d,     0, -b,  c]
        ])
        
        K_global = self.T.T @ K_local @ self.T
        F_int_e = K_global @ u_e
        
        return K_global, F_int_e

    def compute_global_stiffness(self) -> np.ndarray:
        K_e, _ = self.compute_element_state(np.zeros(6))
        return K_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = np.array([
            U_global[self.nodes[0].dofs['ux']], U_global[self.nodes[0].dofs['uy']], U_global[self.nodes[0].dofs['rz']],
            U_global[self.nodes[1].dofs['ux']], U_global[self.nodes[1].dofs['uy']], U_global[self.nodes[1].dofs['rz']]
        ])
        K_global, F_int = self.compute_element_state(u_e)
        F_local = self.T @ F_int
        return {
            'axial_force': F_local[0], 
            'shear_force': F_local[1], 
            'moment_i': F_local[2], 
            'moment_j': F_local[5]
        }
