# fenix_fem/fenix/elements/frame.py
"""Elementos de pórtico/viga 2D.

Agrupa `Frame2DEuler` (Euler-Bernoulli, vigas esbeltas) y
`Frame2DTimoshenko` (Timoshenko, vigas peraltadas). Ambas clases son
**independientes**: heredan directamente de `Element` y cada una replica
internamente su propia construcción de la matriz de transformación — no
se heredan entre sí ni comparten helpers. Conviven en este archivo por
convención temática (familia "pórtico/viga 2D"), no por dependencia de
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


@ElementRegistry.register
class Frame2DEuler(Element):
    """Pórtico/viga 2D basado en Euler-Bernoulli.

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


@ElementRegistry.register
class Frame2DTimoshenko(Element):
    """Pórtico/viga 2D basado en Timoshenko.

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
