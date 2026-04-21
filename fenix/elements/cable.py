# fenix_fem/fenix/elements/cable.py
"""Elementos de cable.

Implementación **independiente** de los elementos de armadura (`Truss*`).
La duplicación de la maquinaria cinemática respecto a `Truss2DCorot` es
deliberada: si la formulación de armaduras cambia, los cables no se ven
arrastrados, y viceversa. Ver docs/specs/Cable2DCorot.md.
"""
import math
from typing import List

import numpy as np

from fenix.core.element import Element
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.registry import ElementRegistry


@ElementRegistry.register
class Cable2DCorot(Element):
    """Cable 2D corotacional — tensión unilateral, grandes rotaciones.

    Dos nodos articulados en el plano, 2 DOFs por nodo (ux, uy). La
    cinemática corotacional (Updated Lagrangian) recalcula longitud y
    cosenos directores en configuración corriente en cada evaluación.

    La propiedad de unilateralidad no vive en esta clase — el elemento
    delega la relación constitutiva al material. Para obtener respuesta
    físicamente correcta de cable, emparejar con un material unilateral
    como `CableMaterial1D`.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos.
    material : Material
        Material 1D (STRAIN_DIM=1). Para comportamiento de cable genuino,
        debe ser unilateral (p. ej. CableMaterial1D).
    A : float
        Área de la sección transversal del cable.
    """

    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 A: float):
        if len(nodes) != 2:
            raise ValueError("El elemento Cable2DCorot requiere exactamente 2 nodos.")

        self.A = A
        super().__init__(element_id, nodes, material)

        coords1 = self.nodes[0].coordinates
        coords2 = self.nodes[1].coordinates
        dx0, dy0 = coords2[0] - coords1[0], coords2[1] - coords1[1]
        self.L0 = math.sqrt(dx0**2 + dy0**2)
        if self.L0 == 0.0:
            raise ValueError("La longitud del cable no puede ser cero.")

    def _current_geometry(self, u_e: np.ndarray):
        """Longitud y cosenos directores en configuración corriente."""
        x1 = self.nodes[0].coordinates[0] + u_e[0]
        y1 = self.nodes[0].coordinates[1] + u_e[1]
        x2 = self.nodes[1].coordinates[0] + u_e[2]
        y2 = self.nodes[1].coordinates[1] + u_e[3]
        dx, dy = x2 - x1, y2 - y1
        l = math.sqrt(dx**2 + dy**2)
        if l == 0.0:
            raise ValueError("Longitud corriente nula en Cable2DCorot.")
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
