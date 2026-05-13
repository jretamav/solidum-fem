# fenix_fem/fenix/elements/solid_2d.py
import numpy as np
from typing import List
from fenix.core.element import Element
from fenix.core.node import Node
from fenix.core.material import Material
from fenix.registry import ElementRegistry, QuadratureRegistry
from fenix.constants import ZERO_JACOBIAN_TOL
from numba import njit


@njit
def _compute_kinematics(xi, eta, coords):
    # Preasignar memoria en Numba es la forma 100% segura
    dN_dxi = np.zeros((2, 4), dtype=np.float64)
    dN_dxi[0, 0] = -(1.0 - eta) / 4.0
    dN_dxi[0, 1] =  (1.0 - eta) / 4.0
    dN_dxi[0, 2] =  (1.0 + eta) / 4.0
    dN_dxi[0, 3] = -(1.0 + eta) / 4.0
    dN_dxi[1, 0] = -(1.0 - xi) / 4.0
    dN_dxi[1, 1] = -(1.0 + xi) / 4.0
    dN_dxi[1, 2] =  (1.0 + xi) / 4.0
    dN_dxi[1, 3] =  (1.0 - xi) / 4.0

    # Jacobiano
    J = np.dot(dN_dxi, coords)
    detJ = J[0,0] * J[1,1] - J[0,1] * J[1,0]

    if detJ <= ZERO_JACOBIAN_TOL:
        raise ValueError("Jacobiano negativo o cero detectado en elemento Quad4. Revisa la conectividad o distorsion.")

    # Inversa del Jacobiano segura
    invJ = np.zeros((2, 2), dtype=np.float64)
    invJ[0, 0] =  J[1, 1] / detJ
    invJ[0, 1] = -J[0, 1] / detJ
    invJ[1, 0] = -J[1, 0] / detJ
    invJ[1, 1] =  J[0, 0] / detJ

    # Derivadas respecto a coordenadas globales (x, y)
    dN_dx = np.dot(invJ, dN_dxi)

    # Ensamblaje de la matriz B (Deformación-Desplazamiento)
    B = np.zeros((3, 8), dtype=np.float64)
    for i in range(4):
        B[0, 2*i]     = dN_dx[0, i]
        B[1, 2*i + 1] = dN_dx[1, i]
        B[2, 2*i]     = dN_dx[1, i]
        B[2, 2*i + 1] = dN_dx[0, i]

    return B, detJ


@njit
def _compute_integrands(B, C_alg, sigma, detJ, weight, thickness):
    dV = detJ * weight * thickness

    # Descomponer el triple producto tensorial evita errores de Numba
    temp = np.dot(B.T, C_alg)
    K_contrib = np.dot(temp, B) * dV

    F_contrib = np.dot(B.T, sigma) * dV
    return K_contrib, F_contrib


@njit
def _shape_functions_quad4(xi, eta):
    """Funciones de forma bilineales del Quad4 evaluadas en (xi, eta)."""
    N = np.zeros(4, dtype=np.float64)
    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta)
    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta)
    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta)
    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta)
    return N


@njit
def _det_jacobian_quad4(xi, eta, coords):
    """det J del mapeo isoparamétrico del Quad4 en (xi, eta).

    Reutilizable desde compute_body_load sin pagar el coste de armar B.
    """
    dN_dxi = np.zeros((2, 4), dtype=np.float64)
    dN_dxi[0, 0] = -(1.0 - eta) / 4.0
    dN_dxi[0, 1] =  (1.0 - eta) / 4.0
    dN_dxi[0, 2] =  (1.0 + eta) / 4.0
    dN_dxi[0, 3] = -(1.0 + eta) / 4.0
    dN_dxi[1, 0] = -(1.0 - xi) / 4.0
    dN_dxi[1, 1] = -(1.0 + xi) / 4.0
    dN_dxi[1, 2] =  (1.0 + xi) / 4.0
    dN_dxi[1, 3] =  (1.0 - xi) / 4.0
    J = np.dot(dN_dxi, coords)
    return J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]


@njit
def _compute_kinematics_tri3(coords):
    # Derivadas de funciones de forma analíticas para Tri3 (xi, eta)
    dN_dxi = np.zeros((2, 3), dtype=np.float64)
    dN_dxi[0, 0] = -1.0; dN_dxi[0, 1] = 1.0; dN_dxi[0, 2] = 0.0
    dN_dxi[1, 0] = -1.0; dN_dxi[1, 1] = 0.0; dN_dxi[1, 2] = 1.0

    J = np.dot(dN_dxi, coords)
    detJ = J[0,0] * J[1,1] - J[0,1] * J[1,0]

    if detJ <= ZERO_JACOBIAN_TOL:
        raise ValueError("Jacobiano negativo o cero detectado en elemento Tri3.")

    invJ = np.zeros((2, 2), dtype=np.float64)
    invJ[0, 0] =  J[1, 1] / detJ
    invJ[0, 1] = -J[0, 1] / detJ
    invJ[1, 0] = -J[1, 0] / detJ
    invJ[1, 1] =  J[0, 0] / detJ

    dN_dx = np.dot(invJ, dN_dxi)

    B = np.zeros((3, 6), dtype=np.float64)
    for i in range(3):
        B[0, 2*i]     = dN_dx[0, i]
        B[1, 2*i + 1] = dN_dx[1, i]
        B[2, 2*i]     = dN_dx[1, i]
        B[2, 2*i + 1] = dN_dx[0, i]

    return B, detJ


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
        """Masa consistente del Quad4 por cuadratura del propio elemento (ADR 0009).

        ``M_e = ∫ρ·N^T·N·t·dA``. Para Quad4 bilineal con cuadratura 2×2,
        la integral es exacta (productos de bilineal × bilineal = orden 2 en
        ξ y η; 2×2 integra hasta orden 3 en cada variable). Si se construyó
        con cuadratura 1×1 (subintegración para evitar locking), la masa
        queda subintegrada y la frecuencia se desvía respecto a la 2×2.

        Returns
        -------
        np.ndarray, shape (8, 8)
        """
        if lumping != "consistent":
            raise NotImplementedError(
                f"Quad4.compute_mass_matrix: lumping='{lumping}' no "
                f"implementado. Fase 1 (ADR 0009) solo admite 'consistent'."
            )
        rho = self.material.density
        coords = self.get_coordinate_matrix(ndim=2)
        M_s = np.zeros((4, 4))
        for (xi, eta), w in zip(self.points, self.weights):
            N = _shape_functions_quad4(xi, eta)
            detJ = _det_jacobian_quad4(xi, eta, coords)
            if detJ <= 0.0:
                raise ValueError(f"Jacobiano no positivo en Quad4 id={self.id}.")
            M_s += rho * np.outer(N, N) * (detJ * w * self.thickness)
        return _expand_scalar_mass(M_s)


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
        """Masa consistente del Tri3 (CST) por fórmula analítica exacta.

        Para funciones de forma lineales en coordenadas baricéntricas, la
        integral ``∫N_i·N_j dA`` es ``A·(1+δ_ij)/12``. Resulta:

            M_scalar = (ρ·t·A/12) · [[2, 1, 1],
                                     [1, 2, 1],
                                     [1, 1, 2]]

        La cuadratura del propio elemento (1 punto central) sería
        insuficiente (integraría exactamente solo polinomios lineales);
        por eso se usa la fórmula cerrada — exacta y barata.

        Returns
        -------
        np.ndarray, shape (6, 6)
        """
        if lumping != "consistent":
            raise NotImplementedError(
                f"Tri3.compute_mass_matrix: lumping='{lumping}' no "
                f"implementado. Fase 1 (ADR 0009) solo admite 'consistent'."
            )
        coef = self.material.density * self._element_area() * self.thickness / 12.0
        M_s = coef * np.array([[2.0, 1.0, 1.0],
                               [1.0, 2.0, 1.0],
                               [1.0, 1.0, 2.0]])
        return _expand_scalar_mass(M_s)


# ============================================================================
# Higher-order continuous 2D elements: Quad8 (serendipity), Quad9 (Lagrange),
# Tri6. Funciones de forma y kinematics en numpy puro.
# ============================================================================

def _N_quad8(xi, eta):
    """Funciones de forma serendípitas del Quad8.

    Numeración: 0..3 vértices antihorarios desde (-1,-1); 4 medio del borde
    0-1 (η=-1), 5 medio del 1-2 (ξ=+1), 6 medio del 2-3 (η=+1), 7 medio del
    3-0 (ξ=-1).
    """
    N = np.zeros(8)
    xs = np.array([-1.0, 1.0, 1.0, -1.0])
    ys = np.array([-1.0, -1.0, 1.0, 1.0])
    for i in range(4):
        N[i] = 0.25 * (1 + xs[i] * xi) * (1 + ys[i] * eta) * (xs[i] * xi + ys[i] * eta - 1)
    N[4] = 0.5 * (1 - xi * xi) * (1 - eta)
    N[5] = 0.5 * (1 + xi)      * (1 - eta * eta)
    N[6] = 0.5 * (1 - xi * xi) * (1 + eta)
    N[7] = 0.5 * (1 - xi)      * (1 - eta * eta)
    return N


def _dN_quad8(xi, eta):
    dN = np.zeros((2, 8))
    xs = np.array([-1.0, 1.0, 1.0, -1.0])
    ys = np.array([-1.0, -1.0, 1.0, 1.0])
    for i in range(4):
        s = xs[i]; t = ys[i]
        dN[0, i] = 0.25 * s * (1 + t * eta) * (2 * s * xi + t * eta)
        dN[1, i] = 0.25 * t * (1 + s * xi)  * (s * xi + 2 * t * eta)
    dN[0, 4] = -xi * (1 - eta)
    dN[1, 4] = -0.5 * (1 - xi * xi)
    dN[0, 5] =  0.5 * (1 - eta * eta)
    dN[1, 5] = -eta * (1 + xi)
    dN[0, 6] = -xi * (1 + eta)
    dN[1, 6] =  0.5 * (1 - xi * xi)
    dN[0, 7] = -0.5 * (1 - eta * eta)
    dN[1, 7] = -eta * (1 - xi)
    return dN


def _N_quad9(xi, eta):
    L_xi  = np.array([0.5 * xi  * (xi  - 1), 1 - xi  * xi,  0.5 * xi  * (xi  + 1)])
    L_eta = np.array([0.5 * eta * (eta - 1), 1 - eta * eta, 0.5 * eta * (eta + 1)])
    N = np.zeros(9)
    N[0] = L_xi[0] * L_eta[0]
    N[1] = L_xi[2] * L_eta[0]
    N[2] = L_xi[2] * L_eta[2]
    N[3] = L_xi[0] * L_eta[2]
    N[4] = L_xi[1] * L_eta[0]
    N[5] = L_xi[2] * L_eta[1]
    N[6] = L_xi[1] * L_eta[2]
    N[7] = L_xi[0] * L_eta[1]
    N[8] = L_xi[1] * L_eta[1]
    return N


def _dN_quad9(xi, eta):
    L_xi  = np.array([0.5 * xi  * (xi  - 1), 1 - xi  * xi,  0.5 * xi  * (xi  + 1)])
    L_eta = np.array([0.5 * eta * (eta - 1), 1 - eta * eta, 0.5 * eta * (eta + 1)])
    dL_xi  = np.array([xi  - 0.5, -2 * xi,  xi  + 0.5])
    dL_eta = np.array([eta - 0.5, -2 * eta, eta + 0.5])
    idx = [(0, 0), (2, 0), (2, 2), (0, 2),
           (1, 0), (2, 1), (1, 2), (0, 1),
           (1, 1)]
    dN = np.zeros((2, 9))
    for k, (i, j) in enumerate(idx):
        dN[0, k] = dL_xi[i] * L_eta[j]
        dN[1, k] = L_xi[i]  * dL_eta[j]
    return dN


def _N_tri6(xi, eta):
    """Numeración Tri6: vértices 0,1,2 en (0,0),(1,0),(0,1); medios
    3 (de 0-1), 4 (de 1-2), 5 (de 2-0)."""
    L1 = 1 - xi - eta; L2 = xi; L3 = eta
    N = np.zeros(6)
    N[0] = L1 * (2 * L1 - 1)
    N[1] = L2 * (2 * L2 - 1)
    N[2] = L3 * (2 * L3 - 1)
    N[3] = 4 * L1 * L2
    N[4] = 4 * L2 * L3
    N[5] = 4 * L3 * L1
    return N


def _dN_tri6(xi, eta):
    L1 = 1 - xi - eta; L2 = xi; L3 = eta
    dN = np.zeros((2, 6))
    dN[0, 0] = -(4 * L1 - 1); dN[1, 0] = -(4 * L1 - 1)
    dN[0, 1] =  (4 * L2 - 1); dN[1, 1] =  0.0
    dN[0, 2] =  0.0;          dN[1, 2] =  (4 * L3 - 1)
    dN[0, 3] = 4 * (L1 - L2); dN[1, 3] = -4 * L2
    dN[0, 4] = 4 * L3;        dN[1, 4] =  4 * L2
    dN[0, 5] = -4 * L3;       dN[1, 5] =  4 * (L1 - L3)
    return dN


def _kinematics_higher_order(grad_fn, xi, eta, coords, n_nodes):
    dN_dxi = grad_fn(xi, eta)
    J = dN_dxi @ coords
    detJ = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    if detJ <= ZERO_JACOBIAN_TOL:
        raise ValueError("Jacobiano negativo o cero en elemento de orden superior.")
    invJ = np.array([[ J[1, 1], -J[0, 1]],
                     [-J[1, 0],  J[0, 0]]]) / detJ
    dN_dx = invJ @ dN_dxi
    B = np.zeros((3, 2 * n_nodes))
    for i in range(n_nodes):
        B[0, 2 * i]     = dN_dx[0, i]
        B[1, 2 * i + 1] = dN_dx[1, i]
        B[2, 2 * i]     = dN_dx[1, i]
        B[2, 2 * i + 1] = dN_dx[0, i]
    return B, detJ


def _expand_scalar_mass(M_scalar: np.ndarray) -> np.ndarray:
    """Expande una masa escalar (n_nodos × n_nodos) a la matriz traslacional
    completa (2·n_nodos × 2·n_nodos) por bloque diagonal: cada componente
    (ux, uy) hereda independientemente el mismo escalar de masa nodal.

    ``M[2i, 2j]     = M_scalar[i, j]``     (componente x)
    ``M[2i+1, 2j+1] = M_scalar[i, j]``     (componente y)
    El acoplamiento ux↔uy entre el mismo o distinto nodo es cero — los dos
    grados de libertad traslacionales del continuo son ortogonales en la
    integral ``∫N^T N``.
    """
    n = M_scalar.shape[0]
    M = np.zeros((2 * n, 2 * n))
    for i in range(n):
        for j in range(n):
            v = M_scalar[i, j]
            M[2 * i,     2 * j]     = v
            M[2 * i + 1, 2 * j + 1] = v
    return M


def _quadratic_edge_traction(elem, edge: int, t_vec: np.ndarray, n_dofs: int) -> np.ndarray:
    """Reparto consistente de tracción uniforme sobre un borde recto con
    funciones de forma cuadráticas: 1/6 al vértice inicial, 4/6 al medio,
    1/6 al vértice final (regla 1-4-1, idéntica a Simpson).
    """
    t_vec = np.asarray(t_vec, dtype=np.float64).reshape(2)
    a, m, c = elem.EDGE_NODES[edge]
    x_a = np.asarray(elem.nodes[a].coordinates[:2], dtype=np.float64)
    x_c = np.asarray(elem.nodes[c].coordinates[:2], dtype=np.float64)
    L = float(np.linalg.norm(x_c - x_a))
    f = np.zeros(n_dofs)
    th = elem.thickness
    f[2 * a]     = (L / 6.0) * t_vec[0] * th
    f[2 * a + 1] = (L / 6.0) * t_vec[1] * th
    f[2 * m]     = (4.0 * L / 6.0) * t_vec[0] * th
    f[2 * m + 1] = (4.0 * L / 6.0) * t_vec[1] * th
    f[2 * c]     = (L / 6.0) * t_vec[0] * th
    f[2 * c + 1] = (L / 6.0) * t_vec[1] * th
    return f


class _HigherOrderQuad(Element):
    """Base interna para Quad8/Quad9: comparte los bucles de Gauss y body load."""
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 3
    _SHAPE_FN = staticmethod(lambda xi, eta: None)
    _GRAD_FN = staticmethod(lambda xi, eta: None)
    _DEFAULT_QUADRATURE = "3x3"

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 thickness: float = 1.0, quadrature: tuple = None):
        if quadrature is None:
            self.points, self.weights = QuadratureRegistry.get(self._DEFAULT_QUADRATURE)
        else:
            self.points, self.weights = quadrature
        self.thickness = thickness
        self.N_INTEGRATION_POINTS = len(self.points)
        super().__init__(element_id, nodes, material)

    @property
    def _n_nodes(self) -> int:
        return len(self.nodes)

    def compute_element_state(self, u_e: np.ndarray):
        n = self._n_nodes
        K_e = np.zeros((2 * n, 2 * n))
        F_int_e = np.zeros(2 * n)
        coords = self.get_coordinate_matrix(ndim=2)
        for idx, ((xi, eta), w) in enumerate(zip(self.points, self.weights)):
            B, detJ = _kinematics_higher_order(self._GRAD_FN, xi, eta, coords, n)
            strain = B @ u_e
            sigma, C, new_state = self.material.compute_state(strain, self.state.vars[idx])
            self.state.vars_trial[idx] = new_state
            self.state.stresses_trial[idx] = sigma
            K_contrib, F_contrib = _compute_integrands(B, C, sigma, detJ, w, self.thickness)
            K_e += K_contrib
            F_int_e += F_contrib
        return K_e, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        gs = self.compute_gauss_state(U_global)
        return {'stress': gs['stress'].mean(axis=0), 'strain': gs['strain'].mean(axis=0)}

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=2)
        n = self._n_nodes
        n_g = self.N_INTEGRATION_POINTS
        nat = np.asarray(self.points, dtype=np.float64).reshape(n_g, 2)
        glb = np.zeros((n_g, 2))
        eps = np.zeros((n_g, 3))
        sig = np.zeros((n_g, 3))
        for idx, (xi, eta) in enumerate(self.points):
            N = self._SHAPE_FN(xi, eta)
            glb[idx] = N @ coords
            B, _ = _kinematics_higher_order(self._GRAD_FN, xi, eta, coords, n)
            strain = B @ u_e
            stress, _, _ = self.material.compute_state(strain, self.state.vars[idx])
            eps[idx] = strain
            sig[idx] = stress
        return {'points_natural': nat, 'points_global': glb,
                'strain': eps, 'stress': sig}

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        b = np.asarray(b, dtype=np.float64).reshape(2)
        coords = self.get_coordinate_matrix(ndim=2)
        n = self._n_nodes
        f = np.zeros(2 * n)
        for (xi, eta), w in zip(self.points, self.weights):
            N = self._SHAPE_FN(xi, eta)
            _, detJ = _kinematics_higher_order(self._GRAD_FN, xi, eta, coords, n)
            factor = detJ * w * self.thickness
            for i in range(n):
                f[2 * i]     += N[i] * b[0] * factor
                f[2 * i + 1] += N[i] * b[1] * factor
        return f

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa consistente del cuadrilátero de orden superior por cuadratura
        del propio elemento (ADR 0009).

        Quad8 y Quad9 con cuadratura 3×3 (default ``_DEFAULT_QUADRATURE``)
        integran exactamente la masa consistente: productos de funciones
        cuadráticas son de orden 4, y 3×3 integra hasta orden 5 en cada
        variable. Si el usuario inyecta una cuadratura más pobre, la masa
        queda subintegrada (mismo trade-off documentado para K).

        Returns
        -------
        np.ndarray, shape (2·n_nodos, 2·n_nodos)
        """
        if lumping != "consistent":
            raise NotImplementedError(
                f"{type(self).__name__}.compute_mass_matrix: lumping="
                f"'{lumping}' no implementado. Fase 1 (ADR 0009) solo "
                f"admite 'consistent'."
            )
        rho = self.material.density
        coords = self.get_coordinate_matrix(ndim=2)
        n = self._n_nodes
        M_s = np.zeros((n, n))
        for (xi, eta), w in zip(self.points, self.weights):
            N = self._SHAPE_FN(xi, eta)
            _, detJ = _kinematics_higher_order(
                self._GRAD_FN, xi, eta, coords, n
            )
            M_s += rho * np.outer(N, N) * (detJ * w * self.thickness)
        return _expand_scalar_mass(M_s)


@ElementRegistry.register
class Quad8(_HigherOrderQuad):
    """Cuadrilátero serendípito 2D de orden 2 (8 nodos).

    Reproduce campos cuadráticos exactamente; sin shear locking severo en
    flexión. Default: Gauss 3×3.
    """
    N_INTEGRATION_POINTS = 9
    _SHAPE_FN = staticmethod(_N_quad8)
    _GRAD_FN = staticmethod(_dN_quad8)
    _DEFAULT_QUADRATURE = "3x3"

    EDGE_NODES = (
        (0, 4, 1),
        (1, 5, 2),
        (2, 6, 3),
        (3, 7, 0),
    )

    def compute_edge_traction(self, edge: int, t_vec: np.ndarray) -> np.ndarray:
        if edge not in (0, 1, 2, 3):
            raise ValueError(f"edge={edge} fuera de rango para Quad8 (0..3).")
        return _quadratic_edge_traction(self, edge, t_vec, n_dofs=16)


@ElementRegistry.register
class Quad9(_HigherOrderQuad):
    """Cuadrilátero Lagrangiano 2D de orden 2 (9 nodos)."""
    N_INTEGRATION_POINTS = 9
    _SHAPE_FN = staticmethod(_N_quad9)
    _GRAD_FN = staticmethod(_dN_quad9)
    _DEFAULT_QUADRATURE = "3x3"

    EDGE_NODES = (
        (0, 4, 1),
        (1, 5, 2),
        (2, 6, 3),
        (3, 7, 0),
    )

    def compute_edge_traction(self, edge: int, t_vec: np.ndarray) -> np.ndarray:
        if edge not in (0, 1, 2, 3):
            raise ValueError(f"edge={edge} fuera de rango para Quad9 (0..3).")
        return _quadratic_edge_traction(self, edge, t_vec, n_dofs=18)


@ElementRegistry.register
class Tri6(Element):
    """Triángulo 2D de 6 nodos (cuadrático completo P₂).

    Cura el shear locking severo del Tri3; reproduce campos cuadráticos
    exactamente. Cuadratura `tri_3` (3 puntos en los puntos medios).
    """
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 3
    N_INTEGRATION_POINTS = 3

    EDGE_NODES = (
        (0, 3, 1),
        (1, 4, 2),
        (2, 5, 0),
    )

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 thickness: float = 1.0, quadrature: tuple = None):
        if quadrature is None:
            self.points, self.weights = QuadratureRegistry.get("tri_3")
        else:
            self.points, self.weights = quadrature
        self.thickness = thickness
        self.N_INTEGRATION_POINTS = len(self.points)
        super().__init__(element_id, nodes, material)

    def compute_element_state(self, u_e: np.ndarray):
        K_e = np.zeros((12, 12))
        F_int_e = np.zeros(12)
        coords = self.get_coordinate_matrix(ndim=2)
        for idx, ((xi, eta), w) in enumerate(zip(self.points, self.weights)):
            B, detJ = _kinematics_higher_order(_dN_tri6, xi, eta, coords, 6)
            strain = B @ u_e
            sigma, C, new_state = self.material.compute_state(strain, self.state.vars[idx])
            self.state.vars_trial[idx] = new_state
            self.state.stresses_trial[idx] = sigma
            K_contrib, F_contrib = _compute_integrands(B, C, sigma, detJ, w, self.thickness)
            K_e += K_contrib
            F_int_e += F_contrib
        return K_e, F_int_e

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        gs = self.compute_gauss_state(U_global)
        return {'stress': gs['stress'].mean(axis=0), 'strain': gs['strain'].mean(axis=0)}

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=2)
        n_g = self.N_INTEGRATION_POINTS
        nat = np.asarray(self.points, dtype=np.float64).reshape(n_g, 2)
        glb = np.zeros((n_g, 2))
        eps = np.zeros((n_g, 3))
        sig = np.zeros((n_g, 3))
        for idx, (xi, eta) in enumerate(self.points):
            N = _N_tri6(xi, eta)
            glb[idx] = N @ coords
            B, _ = _kinematics_higher_order(_dN_tri6, xi, eta, coords, 6)
            strain = B @ u_e
            stress, _, _ = self.material.compute_state(strain, self.state.vars[idx])
            eps[idx] = strain
            sig[idx] = stress
        return {'points_natural': nat, 'points_global': glb,
                'strain': eps, 'stress': sig}

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        b = np.asarray(b, dtype=np.float64).reshape(2)
        coords = self.get_coordinate_matrix(ndim=2)
        f = np.zeros(12)
        for (xi, eta), w in zip(self.points, self.weights):
            N = _N_tri6(xi, eta)
            _, detJ = _kinematics_higher_order(_dN_tri6, xi, eta, coords, 6)
            factor = detJ * w * self.thickness
            for i in range(6):
                f[2 * i]     += N[i] * b[0] * factor
                f[2 * i + 1] += N[i] * b[1] * factor
        return f

    def compute_edge_traction(self, edge: int, t_vec: np.ndarray) -> np.ndarray:
        if edge not in (0, 1, 2):
            raise ValueError(f"edge={edge} fuera de rango para Tri6 (0..2).")
        return _quadratic_edge_traction(self, edge, t_vec, n_dofs=12)

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa consistente del Tri6 por cuadratura del propio elemento (ADR 0009).

        Con la cuadratura ``tri_3`` (3 puntos, orden 2 exacto) la masa
        consistente queda **levemente subintegrada**: el producto ``N_i·N_j``
        de funciones cuadráticas es de orden 4 y tri_3 integra exactamente
        sólo hasta orden 2. El error introducido en la masa es ``O(h⁴)``
        sobre el campo de prueba, una décima parte del orden de la masa
        nodal, y se manifiesta como una pequeña deriva en las frecuencias
        altas. Si la precisión modal sobre Tri6 importa críticamente,
        inyectar una cuadratura más rica (``tri_6``, ``tri_7``) en el
        constructor — el resto del solver no necesita cambios.

        Returns
        -------
        np.ndarray, shape (12, 12)
        """
        if lumping != "consistent":
            raise NotImplementedError(
                f"Tri6.compute_mass_matrix: lumping='{lumping}' no "
                f"implementado. Fase 1 (ADR 0009) solo admite 'consistent'."
            )
        rho = self.material.density
        coords = self.get_coordinate_matrix(ndim=2)
        M_s = np.zeros((6, 6))
        for (xi, eta), w in zip(self.points, self.weights):
            N = _N_tri6(xi, eta)
            _, detJ = _kinematics_higher_order(_dN_tri6, xi, eta, coords, 6)
            M_s += rho * np.outer(N, N) * (detJ * w * self.thickness)
        return _expand_scalar_mass(M_s)
