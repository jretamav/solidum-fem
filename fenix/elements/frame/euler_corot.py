"""``Frame2DEulerCorot`` — viga 2D Euler-Bernoulli corotacional
(Updated Lagrangian, grandes desplazamientos y rotaciones rígidas).
"""
import math
from typing import List

import numpy as np

from fenix.core.element import Element
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.elements.frame._shared import (
    _frame2d_consistent_body_load,
    _frame2d_consistent_mass_local,
    _frame2d_forces_from_local,
)
from fenix.registry import ElementRegistry
from fenix.results import ElementForces


@ElementRegistry.register
class Frame2DEulerCorot(Element):
    """Marco/viga 2D Euler-Bernoulli corotacional (Updated Lagrangian).

    Dos nodos, 3 DOFs por nodo (ux, uy, rz). Captura grandes desplazamientos
    y grandes rotaciones rígidas del elemento con pequeñas deformaciones
    axiales. Rotaciones deformacionales nodales moderadas (|θ̄| ≲ 30°).

    Implementación autónoma: hereda directamente de Element, sin relación
    con Frame2DEuler ni con ningún otro elemento. La construcción de la
    cinemática corotacional se hace íntegramente dentro de la clase.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos.
    material : Material
        Material 1D (STRAIN_DIM=1).
    A : float
        Área de la sección transversal.
    I : float
        Momento de inercia respecto al eje perpendicular al plano.
    """

    DOF_NAMES = ['ux', 'uy', 'rz']
    STRAIN_DIM = 1
    N_INTEGRATION_POINTS = 1

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 A: float, I: float):
        if len(nodes) != 2:
            raise ValueError("El elemento Frame2DEulerCorot requiere exactamente 2 nodos.")

        self.A = A
        self.I = I
        super().__init__(element_id, nodes, material)

        c1 = self.nodes[0].coordinates
        c2 = self.nodes[1].coordinates
        dx0, dy0 = c2[0] - c1[0], c2[1] - c1[1]
        self.L0 = math.sqrt(dx0**2 + dy0**2)
        if self.L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")
        self.alpha0 = math.atan2(dy0, dx0)

    @staticmethod
    def _unwrap(angle: float) -> float:
        """Lleva un ángulo al intervalo (-π, π]."""
        return math.atan2(math.sin(angle), math.cos(angle))

    def _current_geometry(self, u_e):
        c1 = self.nodes[0].coordinates
        c2 = self.nodes[1].coordinates
        x1 = c1[0] + u_e[0]; y1 = c1[1] + u_e[1]
        x2 = c2[0] + u_e[3]; y2 = c2[1] + u_e[4]
        dx, dy = x2 - x1, y2 - y1
        l = math.sqrt(dx**2 + dy**2)
        if l == 0.0:
            raise ValueError("Longitud corriente nula en Frame2DEulerCorot.")
        alpha = math.atan2(dy, dx)
        alpha_e = self._unwrap(alpha - self.alpha0)
        return l, dx / l, dy / l, alpha_e

    def compute_element_state(self, u_e):
        l, c, s, alpha_e = self._current_geometry(u_e)

        u_elong = l - self.L0
        theta_bar_1 = self._unwrap(u_e[2] - alpha_e)
        theta_bar_2 = self._unwrap(u_e[5] - alpha_e)

        epsilon = u_elong / self.L0
        sigma, E_t, new_state = self.material.compute_state(epsilon, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        EA_L0 = E_t * self.A / self.L0
        EI_L0 = E_t * self.I / self.L0
        K_ll = np.array([
            [EA_L0,        0,        0],
            [    0, 4*EI_L0, 2*EI_L0],
            [    0, 2*EI_L0, 4*EI_L0],
        ])

        N  = sigma * self.A
        M1 = K_ll[1, 1] * theta_bar_1 + K_ll[1, 2] * theta_bar_2
        M2 = K_ll[2, 1] * theta_bar_1 + K_ll[2, 2] * theta_bar_2
        F_l = np.array([N, M1, M2])

        B = np.array([
            [  -c,   -s, 0,    c,    s, 0],
            [-s/l,  c/l, 1,  s/l, -c/l, 0],
            [-s/l,  c/l, 0,  s/l, -c/l, 1],
        ])

        F_int_e = B.T @ F_l

        K_material = B.T @ K_ll @ B

        r = np.array([  -c,   -s, 0,    c,    s, 0])
        z = np.array([  -s,    c, 0,    s,   -c, 0])
        K_sigma = (N / l) * np.outer(z, z) \
                - ((M1 + M2) / (l**2)) * (np.outer(r, z) + np.outer(z, r))

        K_T = K_material + K_sigma
        return K_T, F_int_e

    def compute_internal_forces(self, U_global):
        u_e = self.get_local_displacements(U_global)
        l, c, s, _ = self._current_geometry(u_e)
        epsilon = (l - self.L0) / self.L0
        sigma, _, _ = self.material.compute_state(epsilon, self.state.vars[0])
        N = self.A * sigma

        _, F_int = self.compute_element_state(u_e)
        # Proyección a sistema local corotacional para separar axial/cortante
        lam2 = np.array([[c, s], [-s, c]])
        F1 = lam2 @ F_int[0:2]
        F2 = lam2 @ F_int[3:5]
        return {
            'axial_force': F1[0],
            'shear_force': F1[1],
            'moment_i': F_int[2],
            'moment_j': F_int[5],
            'stress': sigma,
            'strain': epsilon,
        }

    def internal_forces(self, U_global: np.ndarray) -> ElementForces:
        """API pública (ADR 0002): N, V, M en nodos i, j, convención §5.

        Usa los cosenos directores de la configuración corriente (corotacional);
        el ``state`` interno ya está comprometido tras ``solve()``.
        """
        u_e = self.get_local_displacements(U_global)
        _, c, s, _ = self._current_geometry(u_e)
        _, F_int = self.compute_element_state(u_e)
        lam2 = np.array([[c, s], [-s, c]])
        F_local = np.empty(6)
        F_local[0:2] = lam2 @ F_int[0:2]
        F_local[2]   = F_int[2]
        F_local[3:5] = lam2 @ F_int[3:5]
        F_local[5]   = F_int[5]
        return _frame2d_forces_from_local(F_local)

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector nodal consistente con peso propio uniforme, evaluado en la
        configuración de referencia.

        Para grandes rotaciones con cargas conservadoras (gravedad), la
        diferencia frente a la integración sobre la geometría corriente es
        de segundo orden y se acepta como aproximación estándar. Ver
        :func:`_frame2d_consistent_body_load`.
        """
        T = self._reference_transform()
        return _frame2d_consistent_body_load(b, self.A, self.L0, T)

    def _reference_transform(self) -> np.ndarray:
        """Matriz de transformación 6×6 en la configuración de referencia."""
        c, s = math.cos(self.alpha0), math.sin(self.alpha0)
        return np.array([
            [ c,  s, 0, 0, 0, 0],
            [-s,  c, 0, 0, 0, 0],
            [ 0,  0, 1, 0, 0, 0],
            [ 0,  0, 0, c, s, 0],
            [ 0,  0, 0,-s, c, 0],
            [ 0,  0, 0, 0, 0, 1],
        ])

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa consistente del Frame2DEulerCorot, evaluada en la configuración
        de referencia (Lagrangeano total, ADR 0009 §1). Análoga a Frame2DEuler
        pero la matriz de transformación se reconstruye desde ``alpha0`` ya que
        el corotacional no almacena ``self.T`` precalculada.
        """
        if lumping != "consistent":
            raise NotImplementedError(
                f"Frame2DEulerCorot.compute_mass_matrix: lumping='{lumping}' "
                f"no implementado. Fase 1 (ADR 0009) solo admite 'consistent'."
            )
        M_local = _frame2d_consistent_mass_local(
            self.material.density, self.A, self.I, self.L0
        )
        T = self._reference_transform()
        return T.T @ M_local @ T
