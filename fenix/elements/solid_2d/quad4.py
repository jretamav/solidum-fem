"""``Quad4`` — cuadrilátero bilineal 2D isoparamétrico."""
from typing import List

import numpy as np

from fenix.core.element import Element
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.elements.solid_2d._shared import (
    _compute_integrands,
    _compute_kinematics,
    _det_jacobian_quad4,
    _expand_scalar_mass,
    _shape_functions_quad4,
)
from fenix.math.mass_lumping import lump_hrz
from fenix.registry import ElementRegistry, QuadratureRegistry


@ElementRegistry.register
class Quad4(Element):
    """
    Elemento cuadrilátero bilineal 2D isoparamétrico con integración de Gauss.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        4 nodos en sentido antihorario.
    material : Material
        Material 2D (STRAIN_DIM=3).
    thickness : float, optional
        Espesor para estado plano. Default 1.0.
    quadrature : tuple, optional
        Tupla (puntos, pesos). Si es None, usa Gauss 2×2.

    Notes
    -----
    Integración reducida ("1×1") es propensa a hourglassing.
    """
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 3
    N_INTEGRATION_POINTS = 4   # default Gauss 2×2; el constructor lo sobreescribe a nivel de instancia si se pasa otra cuadratura.

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 thickness: float = 1.0, quadrature: tuple = None):
        if quadrature is None:
            self.points, self.weights = QuadratureRegistry.get("2x2")
        else:
            self.points, self.weights = quadrature

        self.thickness = thickness
        # ClassVar override: la base lee este atributo para dimensionar el ElementState
        self.N_INTEGRATION_POINTS = len(self.points)

        super().__init__(element_id, nodes, material)

    def compute_element_state(self, u_e: np.ndarray):
        """Calcula K_tangente local y Fuerzas Internas dadas las deformaciones actuales."""
        K_e = np.zeros((8, 8))
        F_int_e = np.zeros(8)

        coords = self.get_coordinate_matrix(ndim=2)

        for idx, (p, w) in enumerate(zip(self.points, self.weights)):
            xi, eta = p
            B, detJ = _compute_kinematics(xi, eta, coords)
            strain = B @ u_e

            sigma, C_tangent, new_state = self.material.compute_state(strain, self.state.vars[idx])
            self.state.vars_trial[idx] = new_state
            self.state.stresses_trial[idx] = sigma

            K_contrib, F_contrib = _compute_integrands(B, C_tangent, sigma, detJ, w, self.thickness)
            K_e += K_contrib
            F_int_e += F_contrib

        return K_e, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Promedios de σ y ε sobre el elemento (utilidad rápida).

        Para visualizar el campo no promediado úsese ``compute_gauss_state``.
        """
        gs = self.compute_gauss_state(U_global)
        if gs['stress'].shape[0] == 0:
            return {'stress': np.zeros(3), 'strain': np.zeros(3)}
        return {
            'stress': gs['stress'].mean(axis=0),
            'strain': gs['strain'].mean(axis=0),
        }

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        """Estado en cada punto de Gauss del elemento.

        Útil para post-proceso fino (mapas de esfuerzos, suavizado nodal,
        extrapolación de Barlow). Reevalúa el material en cada punto a partir
        de las deformaciones derivadas de ``U_global``.

        Returns
        -------
        dict con claves:
        - ``points_natural`` : ndarray (n_g, 2) con (ξ, η).
        - ``points_global``  : ndarray (n_g, 2) con (x, y) en el sistema global.
        - ``strain`` : ndarray (n_g, 3) — Voigt 2D [ε_xx, ε_yy, γ_xy].
        - ``stress`` : ndarray (n_g, 3) — Voigt 2D [σ_xx, σ_yy, σ_xy].
        """
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=2)
        n_g = self.N_INTEGRATION_POINTS

        nat = np.asarray(self.points, dtype=np.float64).reshape(n_g, 2)
        glb = np.zeros((n_g, 2))
        eps = np.zeros((n_g, 3))
        sig = np.zeros((n_g, 3))

        for idx, (xi, eta) in enumerate(self.points):
            N = _shape_functions_quad4(xi, eta)
            glb[idx] = N @ coords
            B, _ = _compute_kinematics(xi, eta, coords)
            strain = B @ u_e
            stress, _, _ = self.material.compute_state(strain, self.state.vars[idx])
            eps[idx] = strain
            sig[idx] = stress

        return {
            'points_natural': nat,
            'points_global': glb,
            'strain': eps,
            'stress': sig,
        }

    # ------------------------------------------------------------------
    # Cargas distribuidas consistentes
    # ------------------------------------------------------------------

    EDGE_NODES = ((0, 1), (1, 2), (2, 3), (3, 0))

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector de cargas nodales consistente con una fuerza de cuerpo uniforme.

        Integra ``∫ N^T b · t dA`` sobre el elemento usando la cuadratura del
        propio elemento. ``b`` se da en coordenadas globales ``(b_x, b_y)``.

        Parameters
        ----------
        b : ndarray (2,)
            Fuerza de cuerpo por unidad de volumen (e.g. ``ρg`` para peso
            propio: ``b = (0, -ρg)``).

        Returns
        -------
        ndarray (8,)
            Cargas equivalentes en los DOFs locales del elemento, ordenadas
            ``[fx_1, fy_1, ..., fx_4, fy_4]``.
        """
        b = np.asarray(b, dtype=np.float64).reshape(2)
        coords = self.get_coordinate_matrix(ndim=2)
        f = np.zeros(8)
        for (xi, eta), w in zip(self.points, self.weights):
            N = _shape_functions_quad4(xi, eta)
            detJ = _det_jacobian_quad4(xi, eta, coords)
            if detJ <= 0.0:
                raise ValueError(f"Jacobiano no positivo en Quad4 id={self.id}.")
            factor = detJ * w * self.thickness
            for i in range(4):
                f[2 * i]     += N[i] * b[0] * factor
                f[2 * i + 1] += N[i] * b[1] * factor
        return f

    def compute_edge_traction(self, edge: int, t_vec: np.ndarray) -> np.ndarray:
        """Vector de cargas nodales consistente con una tracción uniforme en un borde.

        Integra ``∫_Γ N^T t̄ · t dS`` sobre el borde indicado. ``t_vec`` es la
        tracción ``(t_x, t_y)`` en coordenadas globales y se asume **constante
        a lo largo del borde**. Para presión normal o tracción variable, el
        usuario debe descomponer/discretizar manualmente.

        Parameters
        ----------
        edge : int
            Índice del borde según ``EDGE_NODES``: 0=(n0,n1), 1=(n1,n2),
            2=(n2,n3), 3=(n3,n0).
        t_vec : ndarray (2,)
            Tracción uniforme en globales sobre el borde.

        Returns
        -------
        ndarray (8,)
            Cargas equivalentes en los DOFs locales del elemento.
        """
        if edge not in (0, 1, 2, 3):
            raise ValueError(f"edge={edge} fuera de rango para Quad4 (0..3).")
        t_vec = np.asarray(t_vec, dtype=np.float64).reshape(2)
        a, c = self.EDGE_NODES[edge]
        x_a = np.asarray(self.nodes[a].coordinates[:2], dtype=np.float64)
        x_c = np.asarray(self.nodes[c].coordinates[:2], dtype=np.float64)
        L = float(np.linalg.norm(x_c - x_a))
        # Borde recto entre dos nodos: el reparto consistente para tracción
        # uniforme y N lineales es L/2 en cada nodo, exacto sin cuadratura.
        f = np.zeros(8)
        f[2 * a]     = 0.5 * L * t_vec[0] * self.thickness
        f[2 * a + 1] = 0.5 * L * t_vec[1] * self.thickness
        f[2 * c]     = 0.5 * L * t_vec[0] * self.thickness
        f[2 * c + 1] = 0.5 * L * t_vec[1] * self.thickness
        return f

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa del Quad4 por cuadratura del propio elemento (ADR 0009).

        **Consistente** (default): ``M_e = ∫ρ·N^T·N·t·dA``. Para Quad4
        bilineal con cuadratura 2×2 la integral es exacta. Si se construyó
        con cuadratura 1×1 (subintegración para evitar locking), la masa
        queda subintegrada y la frecuencia se desvía respecto a la 2×2.

        **Lumped** (fase 2): HRZ canónico (Hinton-Rock-Zienkiewicz 1976)
        aplicado a la diagonal consistente. Para Quad4 coincide con
        row-sum por simetría — ``ρ·V_e/4`` en cada DOF traslacional.

        Returns
        -------
        np.ndarray, shape (8, 8)
        """
        rho = self.material.density
        coords = self.get_coordinate_matrix(ndim=2)
        M_s = np.zeros((4, 4))
        m_total = 0.0
        for (xi, eta), w in zip(self.points, self.weights):
            N = _shape_functions_quad4(xi, eta)
            detJ = _det_jacobian_quad4(xi, eta, coords)
            if detJ <= 0.0:
                raise ValueError(f"Jacobiano no positivo en Quad4 id={self.id}.")
            weight = detJ * w * self.thickness
            M_s += rho * np.outer(N, N) * weight
            m_total += rho * weight
        M_consistent = _expand_scalar_mass(M_s)
        if lumping == "lumped":
            return lump_hrz(M_consistent, total_mass=m_total, n_translational_dirs=2)
        if lumping != "consistent":
            raise NotImplementedError(
                f"Quad4.compute_mass_matrix: lumping='{lumping}' no "
                f"soportado. Valores admitidos: 'consistent', 'lumped'."
            )
        return M_consistent
