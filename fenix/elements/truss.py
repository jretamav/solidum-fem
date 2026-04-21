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


@ElementRegistry.register
class Truss3DCorot(Truss3D):
    """
    Armadura 3D corotacional (Updated Lagrangian).

    Hereda identidad estructural de Truss3D (mismos DOFs, parámetro A,
    contrato con el material) y redefine la evaluación de estado para
    capturar grandes desplazamientos y rotaciones en el espacio, manteniendo
    la hipótesis de pequeña deformación (|ε| ≲ 1e-2).

    Rigidez tangente = K_M (material, rango 1) + K_G (geométrica, rango 2,
    proyector perpendicular al eje corriente); fuerzas internas proyectadas
    sobre la dirección actual de la barra.
    """

    def _current_geometry(self, u_e: np.ndarray):
        """Longitud y cosenos directores en configuración corriente."""
        c1 = self.nodes[0].coordinates
        c2 = self.nodes[1].coordinates
        z1 = c1[2] if len(c1) > 2 else 0.0
        z2 = c2[2] if len(c2) > 2 else 0.0
        x1 = c1[0] + u_e[0]; y1 = c1[1] + u_e[1]; zn1 = z1 + u_e[2]
        x2 = c2[0] + u_e[3]; y2 = c2[1] + u_e[4]; zn2 = z2 + u_e[5]
        dx, dy, dz = x2 - x1, y2 - y1, zn2 - zn1
        l = math.sqrt(dx**2 + dy**2 + dz**2)
        if l == 0.0:
            raise ValueError("Longitud corriente nula en Truss3DCorot.")
        return l, dx / l, dy / l, dz / l

    @staticmethod
    def _perpendicular_projector(cx: float, cy: float, cz: float) -> np.ndarray:
        """Proyector 3x3 al plano perpendicular al eje: P = I − ê·êᵀ."""
        e = np.array([cx, cy, cz])
        return np.eye(3) - np.outer(e, e)

    def compute_element_state(self, u_e: np.ndarray):
        l, cx, cy, cz = self._current_geometry(u_e)

        epsilon = (l - self.L0) / self.L0
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        N = sigma * self.A
        d = np.array([-cx, -cy, -cz, cx, cy, cz])
        P = self._perpendicular_projector(cx, cy, cz)

        K_M = ((E_t * self.A) / self.L0) * np.outer(d, d)

        K_G = np.zeros((6, 6))
        K_G[:3, :3] = P
        K_G[3:, 3:] = P
        K_G[:3, 3:] = -P
        K_G[3:, :3] = -P
        K_G *= (N / l)

        F_int_e = N * d
        return K_M + K_G, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        l, _, _, _ = self._current_geometry(u_e)
        epsilon = (l - self.L0) / self.L0
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        N = self.A * sigma
        return {'axial_force': N, 'stress': sigma, 'strain': epsilon}

