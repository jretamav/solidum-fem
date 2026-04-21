import math
import numpy as np
from typing import List
from fenix.core.element import Element
from fenix.core.node import Node
from fenix.core.material import Material
from fenix.registry import ElementRegistry


@ElementRegistry.register
class Truss2D(Element):
    """
    Elemento de armadura 2D, primer orden, formulación lineal infinitesimal.

    Barra articulada que transmite únicamente esfuerzo axial. Régimen válido:
    pequeños desplazamientos, pequeñas rotaciones, |ε| ≲ 1e-2. Para grandes
    rotaciones con pequeña deformación usar Truss2DCorot.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos.
    material : Material
        Material 1D (STRAIN_DIM=1).
    A : float
        Área de la sección transversal.
    """
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 A: float):
        if len(nodes) != 2:
            raise ValueError("El elemento Truss2D requiere exactamente 2 nodos.")

        self.A = A
        super().__init__(element_id, nodes, material)

        coords1 = self.nodes[0].coordinates
        coords2 = self.nodes[1].coordinates
        dx, dy = coords2[0] - coords1[0], coords2[1] - coords1[1]
        self.L0 = math.sqrt(dx**2 + dy**2)
        if self.L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")
        self.c = dx / self.L0
        self.s = dy / self.L0

    def compute_element_state(self, u_e: np.ndarray):
        c, s, L = self.c, self.s, self.L0
        B = np.array([-c, -s, c, s]) / L
        epsilon = np.dot(B, u_e)

        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        d = np.array([-c, -s, c, s])
        K_e = ((E_t * self.A) / L) * np.outer(d, d)
        F_int_e = (sigma * self.A) * d
        return K_e, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Utilidad para reportes de post-procesamiento."""
        u_e = self.get_local_displacements(U_global)
        c, s, L = self.c, self.s, self.L0
        B = np.array([-c, -s, c, s]) / L
        epsilon = np.dot(B, u_e)

        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}


@ElementRegistry.register
class Truss2DCorot(Truss2D):
    """
    Armadura 2D corotacional (Updated Lagrangian).

    Hereda identidad estructural de Truss2D (mismos DOFs, parámetros, material)
    y redefine la evaluación de estado para capturar grandes desplazamientos y
    rotaciones manteniendo la hipótesis de pequeña deformación (|ε| ≲ 1e-2).

    Rigidez tangente = K_M (material) + K_G (geométrica); las fuerzas internas
    se proyectan sobre la dirección corriente del eje de la barra.
    """

    def _current_geometry(self, u_e: np.ndarray):
        """Posiciones, longitud y cosenos directores en la configuración corriente."""
        x1 = self.nodes[0].coordinates[0] + u_e[0]
        y1 = self.nodes[0].coordinates[1] + u_e[1]
        x2 = self.nodes[1].coordinates[0] + u_e[2]
        y2 = self.nodes[1].coordinates[1] + u_e[3]
        dx, dy = x2 - x1, y2 - y1
        l = math.sqrt(dx**2 + dy**2)
        if l == 0.0:
            raise ValueError("Longitud corriente nula en Truss2DCorot.")
        return l, dx / l, dy / l

    def compute_element_state(self, u_e: np.ndarray):
        l, c_t, s_t = self._current_geometry(u_e)

        epsilon = (l - self.L0) / self.L0
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        N = sigma * self.A
        d = np.array([-c_t, -s_t, c_t, s_t])
        n = np.array([-s_t, c_t, s_t, -c_t])

        K_M = ((E_t * self.A) / self.L0) * np.outer(d, d)
        K_G = (N / l) * np.outer(n, n)
        F_int_e = N * d

        return K_M + K_G, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        l, _, _ = self._current_geometry(u_e)
        epsilon = (l - self.L0) / self.L0
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}


@ElementRegistry.register
class Truss3D(Element):
    """
    Elemento de armadura 3D, primer orden, formulación lineal infinitesimal.

    Barra articulada en el espacio que transmite únicamente esfuerzo axial.
    Dos nodos, tres DOFs de desplazamiento por nodo. Régimen válido: pequeños
    desplazamientos, pequeñas rotaciones, |ε| ≲ 1e-2. La configuración de
    referencia es fija; no existe variante corotacional en este elemento.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos. Acepta nodos con 2 ó 3 coordenadas
        (completa con z=0 si la tercera falta).
    material : Material
        Material 1D (STRAIN_DIM=1).
    A : float
        Área de la sección transversal.
    """
    DOF_NAMES = ['ux', 'uy', 'uz']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1

    def __init__(self, element_id: int, nodes: List[Node], material: Material, A: float):
        if len(nodes) != 2:
            raise ValueError("El elemento Truss3D requiere exactamente 2 nodos.")

        self.A = A
        super().__init__(element_id, nodes, material)

        coords1 = self.nodes[0].coordinates
        coords2 = self.nodes[1].coordinates
        # Soporte para nodos con 2 ó 3 coordenadas
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

    def compute_element_state(self, u_e: np.ndarray):
        cx, cy, cz, L = self.cx, self.cy, self.cz, self.L0
        B_dir = np.array([-cx, -cy, -cz, cx, cy, cz])
        B = B_dir / L
        epsilon = np.dot(B, u_e)

        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        coef = (E_t * self.A) / L
        K_e = coef * np.outer(B_dir, B_dir)
        F_int_e = (sigma * self.A) * B_dir

        return K_e, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        cx, cy, cz, L = self.cx, self.cy, self.cz, self.L0
        B = np.array([-cx, -cy, -cz, cx, cy, cz]) / L

        epsilon = np.dot(B, u_e)
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}


def _frame_geometry(nodes):
    """Devuelve (L0, c, s, T) para un par de nodos 2D que comparten DOFs ['ux','uy','rz']."""
    coords1 = nodes[0].coordinates
    coords2 = nodes[1].coordinates
    dx, dy = coords2[0] - coords1[0], coords2[1] - coords1[1]
    L0 = math.sqrt(dx**2 + dy**2)
    if L0 == 0.0:
        raise ValueError("La longitud del elemento no puede ser cero.")
    c = dx / L0
    s = dy / L0
    T = np.array([
        [ c,  s,  0,  0,  0,  0],
        [-s,  c,  0,  0,  0,  0],
        [ 0,  0,  1,  0,  0,  0],
        [ 0,  0,  0,  c,  s,  0],
        [ 0,  0,  0, -s,  c,  0],
        [ 0,  0,  0,  0,  0,  1]
    ])
    return L0, c, s, T


@ElementRegistry.register
class Frame2DEuler(Element):
    """
    Elemento de pórtico/viga 2D basado en Euler-Bernoulli.

    Para vigas esbeltas. Transmite axial, cortante y momento flector.

    Notes
    -----
    La no-linealidad se evalúa solo a partir de la deformación axial (E_t escala
    rigidez axial y de flexión por igual). No captura rótulas plásticas.
    """
    DOF_NAMES = ['ux', 'uy', 'rz']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 A: float, I: float):
        if len(nodes) != 2:
            raise ValueError("El elemento Frame2DEuler requiere exactamente 2 nodos.")

        self.A = A
        self.I = I
        super().__init__(element_id, nodes, material)
        self.L0, _, _, self.T = _frame_geometry(self.nodes)

    def compute_element_state(self, u_e: np.ndarray):
        L = self.L0
        u_local = self.T @ u_e

        epsilon = (u_local[3] - u_local[0]) / L
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        EA_L = E_t * self.A / L
        EI_L = E_t * self.I / L
        EI_L2 = EI_L / L
        EI_L3 = EI_L2 / L

        K_local = np.array([
            [ EA_L,        0,        0, -EA_L,        0,        0],
            [    0, 12*EI_L3,  6*EI_L2,     0,-12*EI_L3,  6*EI_L2],
            [    0,  6*EI_L2,   4*EI_L,     0, -6*EI_L2,   2*EI_L],
            [-EA_L,        0,        0,  EA_L,        0,        0],
            [    0,-12*EI_L3, -6*EI_L2,     0, 12*EI_L3, -6*EI_L2],
            [    0,  6*EI_L2,   2*EI_L,     0, -6*EI_L2,   4*EI_L]
        ])

        F_int_local = K_local @ u_local
        # Componente axial usa sigma directo del material (no la aproximación lineal)
        F_int_local[0] = -sigma * self.A
        F_int_local[3] =  sigma * self.A

        K_global = self.T.T @ K_local @ self.T
        F_int_e = self.T.T @ F_int_local
        return K_global, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        _, F_int = self.compute_element_state(u_e)
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


@ElementRegistry.register
class Frame2DTimoshenko(Element):
    """
    Elemento de pórtico/viga 2D basado en Timoshenko.

    Incluye explícitamente la deformación por cortante transversal. Apropiado
    para vigas gruesas, cortas o peraltadas. Previene shear locking.

    Parameters
    ----------
    A, I : float
        Área y momento de inercia respecto al eje Z.
    As : float
        Área efectiva de cortante.
    nu : float, optional
        Poisson para calcular G. Si el material expone `nu`, se usa ese valor.

    Notes
    -----
    Misma limitación que Euler en la no-linealidad: E_t escala toda la rigidez.
    """
    DOF_NAMES = ['ux', 'uy', 'rz']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 A: float, I: float, As: float, nu: float = 0.3):
        if len(nodes) != 2:
            raise ValueError("El elemento Frame2DTimoshenko requiere exactamente 2 nodos.")

        self.A = A
        self.I = I
        self.As = As
        super().__init__(element_id, nodes, material)

        # Obtener nu del material si está disponible; advertir si usamos default
        if hasattr(material, 'nu'):
            self.nu = material.nu
        else:
            if nu == 0.3:
                print(f"  [!] ADVERTENCIA Frame2DTimoshenko (id={element_id}): el material no expone 'nu'. "
                      f"Se usará nu={nu} (valor por defecto). Especifique 'nu' en el YAML del elemento si esto es incorrecto.")
            self.nu = nu

        self.L0, _, _, self.T = _frame_geometry(self.nodes)

    def compute_element_state(self, u_e: np.ndarray):
        L = self.L0
        u_local = self.T @ u_e

        epsilon = (u_local[3] - u_local[0]) / L
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        G = E_t / (2.0 * (1.0 + self.nu))
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

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        _, F_int = self.compute_element_state(u_e)
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
