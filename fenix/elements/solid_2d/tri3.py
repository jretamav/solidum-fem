"""``Tri3`` — triángulo lineal 2D de deformación constante (CST)."""
from typing import List

import numpy as np

from fenix.core.element import Element, validate_lumping_kwarg
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.elements.solid_2d._shared import (
    _compute_integrands,
    _compute_kinematics_tri3,
    _expand_scalar_mass,
)
from fenix.math.mass_lumping import lump_hrz
from fenix.registry import ElementRegistry


@ElementRegistry.register
class Tri3(Element):
    """
    Elemento triangular lineal 2D de deformación constante (CST).

    1 punto de integración central. Propenso a shear locking; preferir Quad4
    salvo en transiciones de malla complejas.
    """
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 3
    N_INTEGRATION_POINTS = 1

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 thickness: float = 1.0):
        self.thickness = thickness
        super().__init__(element_id, nodes, material)

    def compute_element_state(self, u_e: np.ndarray):
        coords = self.get_coordinate_matrix(ndim=2)

        B, detJ = _compute_kinematics_tri3(coords)
        strain = B @ u_e

        sigma, C_tangent, new_state = self.material.compute_state(strain, self.state.vars[0])
        self.state.vars_trial[0] = new_state
        self.state.stresses_trial[0] = sigma

        # 1 punto central, peso 0.5 (área del triángulo en coordenadas naturales)
        K_e, F_int_e = _compute_integrands(B, C_tangent, sigma, detJ, 0.5, self.thickness)
        return K_e, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """σ y ε del único punto central — coincide con el promedio del elemento."""
        gs = self.compute_gauss_state(U_global)
        return {'stress': gs['stress'][0], 'strain': gs['strain'][0]}

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        """Estado en el punto de Gauss (1 punto central) del Tri3.

        Devuelve la misma estructura que ``Quad4.compute_gauss_state``;
        para Tri3 ε y σ son uniformes sobre el elemento (CST).
        """
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=2)

        B, _ = _compute_kinematics_tri3(coords)
        strain = B @ u_e
        stress, _, _ = self.material.compute_state(strain, self.state.vars[0])

        # Punto central en coordenadas baricéntricas (1/3, 1/3): el centroide.
        centroid = coords.mean(axis=0)

        return {
            'points_natural': np.array([[1.0 / 3.0, 1.0 / 3.0]]),
            'points_global': centroid.reshape(1, 2),
            'strain': strain.reshape(1, 3),
            'stress': stress.reshape(1, 3),
        }

    # ------------------------------------------------------------------
    # Cargas distribuidas consistentes
    # ------------------------------------------------------------------

    EDGE_NODES = ((0, 1), (1, 2), (2, 0))

    def _element_area(self) -> float:
        coords = self.get_coordinate_matrix(ndim=2)
        x = coords[:, 0]; y = coords[:, 1]
        return 0.5 * abs(
            (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])
        )

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector de cargas nodales consistente con una fuerza de cuerpo uniforme.

        Para Tri3 con N lineal y b uniforme la integral es exacta:
        cada nodo recibe ``b · A_e · t / 3``.

        Parameters
        ----------
        b : ndarray (2,)
            Fuerza de cuerpo por unidad de volumen en globales.

        Returns
        -------
        ndarray (6,)
            Cargas equivalentes ordenadas ``[fx_1, fy_1, ..., fx_3, fy_3]``.
        """
        b = np.asarray(b, dtype=np.float64).reshape(2)
        A = self._element_area()
        share = A * self.thickness / 3.0
        f = np.zeros(6)
        for i in range(3):
            f[2 * i]     = share * b[0]
            f[2 * i + 1] = share * b[1]
        return f

    def compute_edge_traction(self, edge: int, t_vec: np.ndarray) -> np.ndarray:
        """Vector de cargas nodales consistente con una tracción uniforme en un borde.

        Parameters
        ----------
        edge : int
            Índice del borde según ``EDGE_NODES``: 0=(n0,n1), 1=(n1,n2), 2=(n2,n0).
        t_vec : ndarray (2,)
            Tracción uniforme en globales sobre el borde.

        Returns
        -------
        ndarray (6,)
            Cargas equivalentes en los DOFs locales del elemento.
        """
        if edge not in (0, 1, 2):
            raise ValueError(f"edge={edge} fuera de rango para Tri3 (0..2).")
        t_vec = np.asarray(t_vec, dtype=np.float64).reshape(2)
        a, c = self.EDGE_NODES[edge]
        x_a = np.asarray(self.nodes[a].coordinates[:2], dtype=np.float64)
        x_c = np.asarray(self.nodes[c].coordinates[:2], dtype=np.float64)
        L = float(np.linalg.norm(x_c - x_a))
        f = np.zeros(6)
        f[2 * a]     = 0.5 * L * t_vec[0] * self.thickness
        f[2 * a + 1] = 0.5 * L * t_vec[1] * self.thickness
        f[2 * c]     = 0.5 * L * t_vec[0] * self.thickness
        f[2 * c + 1] = 0.5 * L * t_vec[1] * self.thickness
        return f

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa del Tri3 (CST) por fórmula analítica exacta (ADR 0009).

        **Consistente** (default): para funciones de forma lineales en
        coordenadas baricéntricas, ``∫N_i·N_j dA = A·(1+δ_ij)/12``. Resulta

            M_scalar = (ρ·t·A/12) · [[2, 1, 1],
                                     [1, 2, 1],
                                     [1, 1, 2]]

        La cuadratura del propio elemento (1 punto central) sería
        insuficiente (integraría exactamente solo polinomios lineales);
        por eso se usa la fórmula cerrada.

        **Lumped** (fase 2): HRZ canónico. Para Tri3 coincide con row-sum
        por simetría — ``ρAt/3`` en cada DOF traslacional (masa total ``ρAt``
        repartida equitativamente entre los 3 nodos).

        Returns
        -------
        np.ndarray, shape (6, 6)
        """
        m_total = self.material.density * self._element_area() * self.thickness
        coef = m_total / 12.0
        M_s = coef * np.array([[2.0, 1.0, 1.0],
                               [1.0, 2.0, 1.0],
                               [1.0, 1.0, 2.0]])
        M_consistent = _expand_scalar_mass(M_s)
        validate_lumping_kwarg(lumping, type(self).__name__)
        if lumping == "lumped":
            return lump_hrz(M_consistent, total_mass=m_total, n_translational_dirs=2)
        return M_consistent
