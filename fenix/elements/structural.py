import math
import numpy as np
from typing import List
from fenix.core.element import Element
from fenix.core.node import Node
from fenix.core.material import Material
from fenix.core.element_state import ElementState

class Truss2D(Element):
    """
    Elemento estructural de armadura (biela/tensor) 2D para análisis no lineal.
    
    Formulación de elemento de barra articulada que solo soporta cargas axiales 
    (tracción/compresión). Integra compatibilidad con deformaciones no lineales 
    evaluando la respuesta constitutiva en su punto central.
    
    Parameters
    ----------
    element_id : int
        Identificador único numérico del elemento.
    nodes : List[Node]
        Lista de exactamente 2 objetos `Node` que definen los extremos de la barra.
    material : Material
        Instancia de la ley constitutiva del material.
    A : float
        Área de la sección transversal de la barra.
    """
    def __init__(self, element_id: int, nodes: List[Node], material: Material, A: float):
        super().__init__(element_id, nodes)
        if len(nodes) != 2:
            raise ValueError("El elemento Truss2D requiere exactamente 2 nodos.")
        
        self.material = material
        self.A = A
        
        self.state = ElementState(1, init_stress=0.0)
        
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

    @property
    def state_vars(self):
        return self.state.vars

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'ux', 'uy']

    def commit_state(self):
        self.state.commit()

    def compute_element_state(self, u_e: np.ndarray):
        """Calcula matriz tangente y fuerzas internas dadas las deformaciones."""
        c, s, L = self.c, self.s, self.L0
        
        # Matriz cinemática B
        B = np.array([-c, -s, c, s]) / L
        
        # Deformación actual eps = B * u_e
        epsilon = np.dot(B, u_e)
        
        # Respuesta constitutiva
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma
        
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
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}

class Truss3D(Element):
    """
    Elemento estructural de armadura (biela/tensor) espacial 3D para análisis no lineal.
    
    Extensión tridimensional del elemento `Truss2D`. Soporta componentes de 
    desplazamiento en los ejes X, Y y Z.
    
    Parameters
    ----------
    element_id : int
        Identificador único numérico del elemento.
    nodes : List[Node]
        Lista de exactamente 2 objetos `Node` (con coordenadas tridimensionales).
    material : Material
        Instancia de la ley constitutiva del material.
    A : float
        Área de la sección transversal de la barra.
    """
    def __init__(self, element_id: int, nodes: List[Node], material: Material, A: float):
        super().__init__(element_id, nodes)
        if len(nodes) != 2:
            raise ValueError("El elemento Truss3D requiere exactamente 2 nodos.")
        
        self.material = material
        self.A = A
        
        self.state = ElementState(1, init_stress=0.0)
        
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

    @property
    def state_vars(self):
        return self.state.vars

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'uz', 'ux', 'uy', 'uz']

    def commit_state(self):
        self.state.commit()

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
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma
        
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
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}

class Frame2DEuler(Element):
    """
    Elemento de pórtico/viga 2D basado en la teoría de Euler-Bernoulli.
    
    Formulado para vigas esbeltas donde la deformación por cortante transversal 
    se considera despreciable. Transmite esfuerzos axiales, cortantes y momentos flectores.
    
    Parameters
    ----------
    element_id : int
        Identificador único numérico del elemento.
    nodes : List[Node]
        Lista de exactamente 2 objetos `Node`.
    material : Material
        Instancia de la ley constitutiva del material.
    A : float
        Área de la sección transversal del pórtico.
    I : float
        Momento de inercia de la sección transversal respecto al eje de flexión (Z).
    """
    # ... (previous content) ...
    """
    Notes
    -----
    **Limitación de No-Linealidad:** La no-linealidad de este elemento se evalúa
    únicamente a partir de la deformación axial. El módulo tangente resultante (E_t)
    se usa para escalar tanto la rigidez axial como la de flexión. Esto no captura
    la plastificación progresiva de la sección debida a momentos flectores (rótulas plásticas).
    """
    def __init__(self, element_id: int, nodes: List[Node], material: Material, A: float, I: float):
        super().__init__(element_id, nodes)
        if len(nodes) != 2:
            raise ValueError("El elemento Frame2DEuler requiere exactamente 2 nodos.")
        
        self.material = material
        self.A = A
        self.I = I
        
        self.state = ElementState(1, init_stress=0.0)
        
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

    @property
    def state_vars(self):
        return self.state.vars

    def commit_state(self):
        self.state.commit()

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'rz', 'ux', 'uy', 'rz']

    def compute_element_state(self, u_e: np.ndarray):
        L = self.L0
        u_local = self.T @ u_e
        
        # Deformación axial
        epsilon = (u_local[3] - u_local[0]) / L
        
        # Respuesta constitutiva (evaluada en el centroide)
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma
        
        EA_L = E_t * self.A / L
        EI_L = E_t * self.I / L
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
        
        # Fuerzas internas locales (usando sigma real para componente axial)
        F_int_local = K_local @ u_local
        F_int_local[0] = -sigma * self.A
        F_int_local[3] =  sigma * self.A
        
        K_global = self.T.T @ K_local @ self.T
        F_int_e = self.T.T @ F_int_local
        
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
        
        u_local = self.T @ u_e
        epsilon = (u_local[3] - u_local[0]) / self.L0
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        
        return {
            'axial_force': F_local[0], 
            'shear_force': F_local[1], 
            'moment_i': F_local[2], 
            'moment_j': F_local[5],
            'stress': sigma,
            'strain': epsilon
        }

class Frame2DTimoshenko(Element):
    """
    Elemento de pórtico/viga 2D basado en la teoría de Timoshenko.
    
    Formulado para vigas gruesas, cortas o peraltadas. Incluye explícitamente 
    la contribución de la deformación por cortante transversal en su matriz de rigidez.
    
    Parameters
    ----------
    element_id : int
        Identificador único numérico del elemento.
    nodes : List[Node]
        Lista de exactamente 2 objetos `Node`.
    material : Material
        Instancia de la ley constitutiva del material.
    A : float
        Área total de la sección transversal del pórtico.
    I : float
        Momento de inercia de la sección transversal respecto al eje de flexión (Z).
    As : float
        Área efectiva de cortante (e.g., 5/6 del área para una sección rectangular).
    nu : float, optional
        Relación de Poisson utilizada para calcular el módulo de rigidez al cortante (G).
        Si no se provee, intenta extraerlo del `material`.
        
    Notes
    -----
    Previene automáticamente el bloqueo por cortante (*shear locking*) mediante 
    su formulación analítica exacta con factores de forma.
    """
    # ... (previous content) ...
    """
    Notes
    -----
    **Limitación de No-Linealidad:** Al igual que el elemento de Euler, la no-linealidad
    se evalúa a partir de la deformación axial. El módulo tangente (E_t) resultante
    escala toda la matriz de rigidez, lo que es una simplificación del comportamiento real de flexión no-lineal.
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
        
        self.state = ElementState(1, init_stress=0.0)
        
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

    @property
    def state_vars(self):
        return self.state.vars

    def commit_state(self):
        self.state.commit()

    def get_dofs(self) -> List[str]:
        return ['ux', 'uy', 'rz', 'ux', 'uy', 'rz']

    def compute_element_state(self, u_e: np.ndarray):
        L = self.L0
        u_local = self.T @ u_e
        
        # Deformación axial
        epsilon = (u_local[3] - u_local[0]) / L
        
        # Respuesta constitutiva (evaluada en el centroide)
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma
            
        G = E_t / (2.0 * (1.0 + self.nu))
        
        # Factor de rigidez al cortante de Timoshenko
        Phi = (12.0 * E_t * self.I) / (G * self.As * (L**2))
        
        EA_L = E_t * self.A / L
        EI_L = E_t * self.I / L
        
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
        
        F_int_local = K_local @ u_local
        F_int_local[0] = -sigma * self.A
        F_int_local[3] =  sigma * self.A
        
        K_global = self.T.T @ K_local @ self.T
        F_int_e = self.T.T @ F_int_local
        
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
        
        u_local = self.T @ u_e
        epsilon = (u_local[3] - u_local[0]) / self.L0
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        
        return {
            'axial_force': F_local[0], 
            'shear_force': F_local[1], 
            'moment_i': F_local[2], 
            'moment_j': F_local[5],
            'stress': sigma,
            'strain': epsilon
        }
