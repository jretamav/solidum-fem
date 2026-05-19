"""Utilidades compartidas por los elementos del paquete ``fenix.elements.solid_2d``.

Cuatro grupos:

1. **Kinematics Numba** para los elementos lineales (Quad4, Tri3):
   :func:`_compute_kinematics`, :func:`_compute_kinematics_tri3`,
   :func:`_compute_integrands`, :func:`_shape_functions_quad4`,
   :func:`_det_jacobian_quad4`.
2. **Funciones de forma y derivadas en numpy puro** para los elementos
   de orden superior (Quad8, Quad9, Tri6): :func:`_N_quad8`,
   :func:`_dN_quad8`, :func:`_N_quad9`, :func:`_dN_quad9`,
   :func:`_N_tri6`, :func:`_dN_tri6`.
3. **Kinematics genérico de orden superior** y reparto consistente de
   tracciones cuadráticas: :func:`_kinematics_higher_order`,
   :func:`_quadratic_edge_traction`. Más el expansor escalar→bloque para
   masa traslacional: :func:`_expand_scalar_mass`.
4. **Base interna** :class:`_HigherOrderSolid2D` que comparten Quad8,
   Quad9 y Tri6 (mismo bucle de Gauss, distinto par
   ``_SHAPE_FN``/``_GRAD_FN``, y cuadratura específica para masa cuando
   la del elemento subintegra el producto cuadrático×cuadrático).

Todos los símbolos llevan guion bajo: son privados al paquete. Los tres
o cuatro que tests existentes importan se reexportan desde
``__init__.py`` para preservar compatibilidad.
"""
from typing import List

import numpy as np
from numba import njit

from fenix.constants import ZERO_JACOBIAN_TOL
from fenix.core.element import Element, validate_lumping_kwarg
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.math.mass_lumping import lump_hrz
from fenix.registry import QuadratureRegistry


# ---------------------------------------------------------------------------
# Kinematics Numba — elementos lineales Quad4 y Tri3.
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Funciones de forma y derivadas — elementos de orden superior (numpy puro).
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Base interna compartida por todos los elementos sólidos 2D isoparamétricos
# de orden superior (cuadrático o más): Quad8, Quad9, Tri6.
# ---------------------------------------------------------------------------

class _HigherOrderSolid2D(Element):
    """Base interna para sólidos 2D isoparamétricos de orden superior.

    Comparte los bucles de Gauss para K, F_int, gauss state, body load y
    masa consistente. Cada subclase declara las funciones de forma
    (``_SHAPE_FN``, ``_GRAD_FN``), la cuadratura por defecto del elemento
    (``_DEFAULT_QUADRATURE``) y, opcionalmente, una cuadratura
    independiente para la matriz de masa (``_MASS_QUADRATURE``). El
    número de nodos se infiere de ``len(self.nodes)``: la base no
    presupone si la geometría es cuadrilátera o triangular.
    """
    DOF_NAMES = ['ux', 'uy']
    STRAIN_DIM = 3
    _SHAPE_FN = staticmethod(lambda xi, eta: None)
    _GRAD_FN = staticmethod(lambda xi, eta: None)
    _DEFAULT_QUADRATURE = "3x3"
    # Cuadratura específica para masa cuando la del elemento subintegra el
    # producto cuadrático×cuadrático (orden 4). ``None`` = usa la del
    # elemento (caso Quad8/Quad9 con 3×3, orden 5 — masa exacta).
    _MASS_QUADRATURE: str | None = None

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
        """API legada: devuelve ``{stress, strain}`` promediados sobre
        los puntos de Gauss del elemento (consumida por el VTK exporter
        y scripts de post-proceso libre).

        **No es el contrato del ADR 0002** — ese es :meth:`internal_forces`
        que devuelve ``ElementForces`` para barras/vigas/cables y ``None``
        para sólidos (ver docstring base en :class:`Element` para la
        deuda técnica documentada). Para acceso por punto de integración
        usar :meth:`compute_gauss_state`, que esta función agrega.
        """
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
        """Masa del elemento sólido 2D de orden superior (ADR 0009).

        **Consistente** (default): cuando la cuadratura del elemento integra
        exactamente el producto cuadrático×cuadrático (Quad8/Quad9 con 3×3,
        orden 5), la masa se integra con esa misma cuadratura — caso
        ``_MASS_QUADRATURE = None``. Cuando subintegra (Tri6 con ``tri_3``,
        orden 2, frente al producto de orden 4 de la masa), la subclase
        declara ``_MASS_QUADRATURE`` con una regla específica
        suficientemente exacta (e.g. ``tri_6``, Dunavant orden 4), evitando
        modos nulos espurios en M.

        **Lumped** (fase 2): HRZ canónico aplicado a la diagonal
        consistente. En Tri6/Quad8/Quad9 (con nodos intermedios o
        centroide) el HRZ preserva masa total y mantiene la
        proporcionalidad de la diagonal; row-sum daría masas negativas en
        los vértices y debe evitarse (Bathe FEP §9.2.4).

        Returns
        -------
        np.ndarray, shape (2·n_nodos, 2·n_nodos)
        """
        rho = self.material.density
        coords = self.get_coordinate_matrix(ndim=2)
        n = self._n_nodes
        if self._MASS_QUADRATURE is None:
            mass_points, mass_weights = self.points, self.weights
        else:
            mass_points, mass_weights = QuadratureRegistry.get(self._MASS_QUADRATURE)
        M_s = np.zeros((n, n))
        m_total = 0.0
        for (xi, eta), w in zip(mass_points, mass_weights):
            N = self._SHAPE_FN(xi, eta)
            _, detJ = _kinematics_higher_order(
                self._GRAD_FN, xi, eta, coords, n
            )
            weight = detJ * w * self.thickness
            M_s += rho * np.outer(N, N) * weight
            m_total += rho * weight
        M_consistent = _expand_scalar_mass(M_s)
        validate_lumping_kwarg(lumping, type(self).__name__)
        if lumping == "lumped":
            return lump_hrz(M_consistent, total_mass=m_total, n_translational_dirs=2)
        return M_consistent
