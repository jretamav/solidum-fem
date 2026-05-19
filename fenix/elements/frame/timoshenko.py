"""``Frame2DTimoshenko`` — viga 2D Timoshenko con deformación por cortante.
"""
from typing import List

import numpy as np

from fenix.core.element import Element
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.elements.frame._shared import (
    _frame2d_consistent_body_load,
    _frame2d_consistent_mass_local,
    _frame2d_forces_from_local,
    _frame2d_lumped_mass_local,
    _log,
    build_geometry_2d,
)
from fenix.registry import ElementRegistry
from fenix.results import ElementForces


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
                _log.warning(
                    f"Frame2DTimoshenko (id={element_id}): el material no expone 'nu'. "
                    f"Se usará nu={nu} (default). Especifique 'nu' en el YAML del elemento si esto es incorrecto."
                )
            self.nu = nu

        self.L0, self.c, self.s, self.T = build_geometry_2d(self.nodes)

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
        # Materiales no lineales sólo capturan el axial uniforme; los
        # términos transversales/flexionales se quedan en ``K_local @ u_local``
        # con ``E_tangent`` global — sin plasticidad por flexión hasta que
        # entre ``FiberSection`` (auditoría H-2.5, ver memoria
        # ``project_pendiente_fiber_section.md``).
        F_int_local[0] = -sigma * self.A
        F_int_local[3] =  sigma * self.A

        K_global = self.T.T @ K_local @ self.T
        F_int_e = self.T.T @ F_int_local
        return K_global, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        _, F_int = self.compute_element_state(u_e)
        # Doble rotación ``T @ (T.T @ F_int_local)`` = identidad (T ortogonal);
        # cosmético por uniformidad del contrato. Ver H-2.8.
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
        """Vector nodal consistente con peso propio uniforme. Para Timoshenko
        con interpolación Hermite-blended, la fórmula coincide con la de
        Euler-Bernoulli en el caso de carga uniforme. Ver
        :func:`_frame2d_consistent_body_load`.
        """
        return _frame2d_consistent_body_load(b, self.A, self.L0, self.T)

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Matriz de masa del Frame2DTimoshenko en ejes globales (ADR 0009).

        **Consistente**: idéntica en forma a la de Frame2DEuler — la
        densidad lineal ``ρA`` no depende de ``As`` ni de ``Φ``; la versión
        Hermite-blended introduce términos ``O(Φ²)`` despreciables. Ver
        :func:`_frame2d_consistent_mass_local`.

        **Lumped** (fase 2): nodal directo ``ρAL/2`` traslacional + ``ρIL/2``
        rotacional por nodo. Ver :func:`_frame2d_lumped_mass_local`.
        """
        rho = self.material.density
        if lumping == "lumped":
            M_local = _frame2d_lumped_mass_local(rho, self.A, self.I, self.L0)
        elif lumping == "consistent":
            M_local = _frame2d_consistent_mass_local(rho, self.A, self.I, self.L0)
        else:
            raise NotImplementedError(
                f"Frame2DTimoshenko.compute_mass_matrix: lumping='{lumping}' "
                f"no soportado. Valores admitidos: 'consistent', 'lumped'."
            )
        return self.T.T @ M_local @ self.T
