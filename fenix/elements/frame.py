# fenix_fem/fenix/elements/frame.py
"""Elementos de marco/viga 2D.

Agrupa `Frame2DEuler` (Euler-Bernoulli, vigas esbeltas) y
`Frame2DTimoshenko` (Timoshenko, vigas peraltadas). Ambas clases son
**independientes**: heredan directamente de `Element` y cada una replica
internamente su propia construcción de la matriz de transformación — no
se heredan entre sí ni comparten helpers. Conviven en este archivo por
convención temática (familia "marco/viga 2D"), no por dependencia de
código.

Ver docs/specs/Frame2DEuler.md y docs/specs/Frame2DTimoshenko.md.
"""
import math
from typing import List

import numpy as np

from fenix.core.element import Element
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.registry import ElementRegistry
from fenix.results import ElementForces


def _frame2d_forces_from_local(F_local: np.ndarray) -> ElementForces:
    """Traduce un vector F_local (6,) en convención stress-resultant interna a
    ``ElementForces`` en convención de viga estructural (Reglas.md §5).

    F_local tiene semántica ``K·u = F_ext``: fuerzas/momentos externos sobre los
    nodos del elemento en ejes locales, layout ``[Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]``.

    Mapeo §5:
      - Nodo i (cara −x): N = −F_local[0], V = +F_local[1], M = −F_local[2]
      - Nodo j (cara +x): N = +F_local[3], V = −F_local[4], M = +F_local[5]

    Verificado con cantilever (carga en punta → M hogging) y cantilever con
    momento en punta (M constante sagging).
    """
    return ElementForces(
        kind="frame2d",
        components={
            "N": np.array([-F_local[0],  F_local[3]]),
            "V": np.array([ F_local[1], -F_local[4]]),
            "M": np.array([-F_local[2],  F_local[5]]),
        },
    )


@ElementRegistry.register
class Frame2DEuler(Element):
    """Marco/viga 2D basado en Euler-Bernoulli.

    Dos nodos rígidamente conectados, 3 DOFs por nodo (ux, uy, rz). Transmite
    esfuerzo axial, cortante y momento flector. Régimen de linealidad
    geométrica y régimen de validez Euler-Bernoulli (vigas esbeltas,
    L/h ≳ 10).

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos.
    material : Material
        Material 1D (STRAIN_DIM=1). E_tangent escala toda la matriz local
        (aproximación: no hay plasticidad distribuida en la sección).
    A : float
        Área de la sección transversal.
    I : float
        Momento de inercia respecto al eje perpendicular al plano (z).
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

        self.L0, self.c, self.s, self.T = self._build_geometry(self.nodes)

    @staticmethod
    def _build_geometry(nodes: List[Node]):
        """Calcula longitud, cosenos directores y matriz T 6×6 (una sola vez)."""
        coords1 = nodes[0].coordinates
        coords2 = nodes[1].coordinates
        dx, dy = coords2[0] - coords1[0], coords2[1] - coords1[1]
        L0 = math.sqrt(dx**2 + dy**2)
        if L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")
        c = dx / L0
        s = dy / L0
        T = np.array([
            [ c,  s, 0, 0, 0, 0],
            [-s,  c, 0, 0, 0, 0],
            [ 0,  0, 1, 0, 0, 0],
            [ 0,  0, 0, c, s, 0],
            [ 0,  0, 0,-s, c, 0],
            [ 0,  0, 0, 0, 0, 1],
        ])
        return L0, c, s, T

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
            [ EA_L,         0,         0, -EA_L,         0,         0],
            [    0,  12*EI_L3,   6*EI_L2,     0, -12*EI_L3,   6*EI_L2],
            [    0,   6*EI_L2,    4*EI_L,     0,  -6*EI_L2,    2*EI_L],
            [-EA_L,         0,         0,  EA_L,         0,         0],
            [    0, -12*EI_L3,  -6*EI_L2,     0,  12*EI_L3,  -6*EI_L2],
            [    0,   6*EI_L2,    2*EI_L,     0,  -6*EI_L2,    4*EI_L],
        ])

        F_int_local = K_local @ u_local
        # Componente axial: usar sigma directo del material (permite materiales no-lineales)
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
            'strain': epsilon,
        }

    def internal_forces(self, U_global: np.ndarray) -> ElementForces:
        """API pública (ADR 0002): N, V, M en nodos i, j, convención §5."""
        u_e = self.get_local_displacements(U_global)
        _, F_int = self.compute_element_state(u_e)
        F_local = self.T @ F_int
        return _frame2d_forces_from_local(F_local)


@ElementRegistry.register
class Frame2DTimoshenko(Element):
    """Marco/viga 2D basado en Timoshenko.

    Dos nodos rígidamente conectados, 3 DOFs por nodo (ux, uy, rz). Incluye
    deformación por cortante transversal; apropiado para vigas peraltadas
    o cortas (L/h ≲ 10). El factor Φ = 12·E·I/(G·A_s·L²) corrige la rigidez
    para evitar shear locking en el límite esbelto.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        Exactamente 2 nodos extremos.
    material : Material
        Material 1D (STRAIN_DIM=1). Si expone atributo `nu`, se usa ese
        valor para el módulo de corte; en su defecto se usa el parámetro
        nu del elemento.
    A : float
        Área de la sección transversal.
    I : float
        Momento de inercia respecto al eje perpendicular al plano (z).
    As : float
        Área efectiva de cortante.
    nu : float, optional
        Coeficiente de Poisson si el material no lo expone. Default: 0.3.
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

        # Fuente del Poisson: material.nu si existe; si no, parámetro del elemento.
        if hasattr(material, 'nu'):
            self.nu = material.nu
        else:
            if nu == 0.3:
                print(f"  [!] ADVERTENCIA Frame2DTimoshenko (id={element_id}): "
                      f"el material no expone 'nu'. Se usará nu={nu} (default). "
                      f"Especifique 'nu' en el YAML del elemento si esto es incorrecto.")
            self.nu = nu

        self.L0, self.c, self.s, self.T = self._build_geometry(self.nodes)

    @staticmethod
    def _build_geometry(nodes: List[Node]):
        """Calcula longitud, cosenos directores y matriz T 6×6 (una sola vez)."""
        coords1 = nodes[0].coordinates
        coords2 = nodes[1].coordinates
        dx, dy = coords2[0] - coords1[0], coords2[1] - coords1[1]
        L0 = math.sqrt(dx**2 + dy**2)
        if L0 == 0.0:
            raise ValueError("La longitud del elemento no puede ser cero.")
        c = dx / L0
        s = dy / L0
        T = np.array([
            [ c,  s, 0, 0, 0, 0],
            [-s,  c, 0, 0, 0, 0],
            [ 0,  0, 1, 0, 0, 0],
            [ 0,  0, 0, c, s, 0],
            [ 0,  0, 0,-s, c, 0],
            [ 0,  0, 0, 0, 0, 1],
        ])
        return L0, c, s, T

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
        c_coef = (4 + Phi) * EI_L / (1 + Phi)
        d_coef = (2 - Phi) * EI_L / (1 + Phi)

        K_local = np.array([
            [ EA_L,  0,       0, -EA_L,  0,       0],
            [    0,  a,       b,     0, -a,       b],
            [    0,  b,  c_coef,     0, -b,  d_coef],
            [-EA_L,  0,       0,  EA_L,  0,       0],
            [    0, -a,      -b,     0,  a,      -b],
            [    0,  b,  d_coef,     0, -b,  c_coef],
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
            'strain': epsilon,
        }

    def internal_forces(self, U_global: np.ndarray) -> ElementForces:
        """API pública (ADR 0002): N, V, M en nodos i, j, convención §5."""
        u_e = self.get_local_displacements(U_global)
        _, F_int = self.compute_element_state(u_e)
        F_local = self.T @ F_int
        return _frame2d_forces_from_local(F_local)



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
