"""``Frame2DEuler`` — viga 2D basada en Euler-Bernoulli (vigas esbeltas).
"""
from typing import List

import numpy as np

from solidum.core.element import Element
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.elements.frame._shared import (
    _frame2d_consistent_body_load,
    _frame2d_consistent_mass_local,
    _frame2d_forces_from_local,
    _frame2d_lumped_mass_local,
    build_geometry_2d,
)
from solidum.registry import ElementRegistry
from solidum.results import ElementForces


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

        self.L0, self.c, self.s, self.T = build_geometry_2d(self.nodes)

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
        # Componente axial: usar sigma directo del material — permite
        # materiales no lineales en el axial uniforme (Elastoplastic1D,
        # daño 1D). Los términos transversal/flexional siguen siendo
        # ``K_local @ u_local`` con ``E_tangent`` aplicado globalmente:
        # **no captura plasticidad por flexión** (las fibras inferiores
        # en tracción no plastifican independientemente). La plasticidad
        # distribuida en la sección requiere ``FiberSection`` (auditoría
        # H-2.5; ver `project_pendiente_fiber_section.md` en memoria).
        F_int_local[0] = -sigma * self.A
        F_int_local[3] =  sigma * self.A

        K_global = self.T.T @ K_local @ self.T
        F_int_e = self.T.T @ F_int_local
        return K_global, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        _, F_int = self.compute_element_state(u_e)
        # ``compute_element_state`` aplica T.T para rotar F_int_local → F_int (globales).
        # Aquí lo rotamos de vuelta a locales para extraer N/V/M. La doble rotación
        # ``self.T @ (self.T.T @ F_int_local)`` es identidad porque T es ortogonal
        # (rotación pura) — pagamos un producto 6×6·6 cosmético por mantener el
        # contrato uniforme de ``compute_element_state`` devolviendo en globales
        # (auditoría H-2.8, severidad bajo).
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

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector nodal consistente con peso propio uniforme. Ver
        :func:`_frame2d_consistent_body_load`.
        """
        return _frame2d_consistent_body_load(b, self.A, self.L0, self.T)

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Matriz de masa del Frame2DEuler en ejes globales (ADR 0009).

        **Consistente** (default, fase 1): axial lineal + flexional
        Hermitiana cúbica + inercia rotacional de sección ``ρI``. Ver
        :func:`_frame2d_consistent_mass_local`.

        **Lumped** (fase 2): nodal directo ``ρAL/2`` traslacional + ``ρIL/2``
        rotacional por nodo. Ver :func:`_frame2d_lumped_mass_local`.

        Ambos casos pasan por ``T^T · M_local · T``. Para el lumped, T es
        ortogonal por bloques compatibles con la diagonalidad: el resultado
        global también es diagonal.

        Returns
        -------
        np.ndarray, shape (6, 6)
        """
        rho = self.material.density
        if lumping == "lumped":
            M_local = _frame2d_lumped_mass_local(rho, self.A, self.I, self.L0)
        elif lumping == "consistent":
            M_local = _frame2d_consistent_mass_local(rho, self.A, self.I, self.L0)
        else:
            raise NotImplementedError(
                f"Frame2DEuler.compute_mass_matrix: lumping='{lumping}' no "
                f"soportado. Valores admitidos: 'consistent', 'lumped'."
            )
        return self.T.T @ M_local @ self.T
