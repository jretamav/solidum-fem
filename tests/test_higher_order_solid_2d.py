"""Tests de los elementos sólidos 2D de orden superior: Quad8, Quad9, Tri6.

Cubre:
- Forma y partición de la unidad de las funciones de forma.
- Identidad nodal: N_i en el nodo i = 1, en los demás = 0.
- Derivadas: Σ ∂N_i/∂ξ = 0 y Σ ∂N_i/∂η = 0 (rigid body translation).
- Patch test cuadrático: con un único elemento sometido a u = (x², 0)
  impuesto en sus nodos, ε_xx en cada Gauss = 2x — sólo posible si el
  elemento reproduce campos cuadráticos exactamente. Tri3/Quad4 fallan
  este test; Quad8/Quad9/Tri6 lo pasan.
- Cargas distribuidas: peso propio (Σ = b·A·t) y tracción uniforme con
  reparto consistente 1/6, 4/6, 1/6.
- Patch físico uniaxial sobre un Quad8 con tracción de borde: verifica
  σ_xx = p exacto en todos los puntos de Gauss.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.material import Material
from fenix.core.node import Node
from fenix.elements.solid_2d import Quad8, Quad9, Tri6
from fenix.elements.solid_2d import _N_quad8, _N_quad9, _N_tri6, _dN_quad8, _dN_quad9, _dN_tri6


class _Elastic2D(Material):
    STRAIN_DIM = 3

    def __init__(self, E=1000.0, nu=0.3):
        self.E = E; self.nu = nu
        fac = E / (1.0 - nu * nu)
        self.C = fac * np.array([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0],
        ])

    def compute_state(self, strain, sv=None):
        return self.C @ strain, self.C, sv


# Coordenadas de referencia de los nodos de cada elemento, en su orden canónico
_QUAD8_NAT = np.array([
    [-1, -1], [1, -1], [1, 1], [-1, 1],
    [0, -1], [1, 0], [0, 1], [-1, 0],
])
_QUAD9_NAT = np.vstack([_QUAD8_NAT, [[0, 0]]])
_TRI6_NAT = np.array([
    [0, 0], [1, 0], [0, 1],
    [0.5, 0], [0.5, 0.5], [0, 0.5],
])


class TestShapeFunctions(unittest.TestCase):
    def test_quad8_partition_of_unity_and_kronecker(self):
        for xi, eta in [(-0.7, 0.3), (0.0, 0.0), (0.5, -0.5), (1.0, 1.0)]:
            N = _N_quad8(xi, eta)
            self.assertAlmostEqual(N.sum(), 1.0, places=12,
                                   msg=f"Σ N ≠ 1 en ({xi},{eta})")
        for k, (xi_k, eta_k) in enumerate(_QUAD8_NAT):
            N = _N_quad8(float(xi_k), float(eta_k))
            for j in range(8):
                expected = 1.0 if j == k else 0.0
                self.assertAlmostEqual(N[j], expected, places=12,
                                       msg=f"Quad8 N[{j}] en nodo {k}")

    def test_quad9_partition_of_unity_and_kronecker(self):
        for xi, eta in [(-0.7, 0.3), (0.0, 0.0), (0.5, -0.5)]:
            N = _N_quad9(xi, eta)
            self.assertAlmostEqual(N.sum(), 1.0, places=12)
        for k, (xi_k, eta_k) in enumerate(_QUAD9_NAT):
            N = _N_quad9(float(xi_k), float(eta_k))
            for j in range(9):
                expected = 1.0 if j == k else 0.0
                self.assertAlmostEqual(N[j], expected, places=12)

    def test_tri6_partition_of_unity_and_kronecker(self):
        for xi, eta in [(0.2, 0.3), (1.0/3.0, 1.0/3.0), (0.5, 0.0), (0.0, 0.5)]:
            N = _N_tri6(xi, eta)
            self.assertAlmostEqual(N.sum(), 1.0, places=12)
        for k, (xi_k, eta_k) in enumerate(_TRI6_NAT):
            N = _N_tri6(float(xi_k), float(eta_k))
            for j in range(6):
                expected = 1.0 if j == k else 0.0
                self.assertAlmostEqual(N[j], expected, places=12)

    def test_derivatives_sum_to_zero(self):
        # Σ ∂N_i/∂ξ = 0 y Σ ∂N_i/∂η = 0 (consistencia de cuerpo rígido)
        for fn in (_dN_quad8, _dN_quad9, _dN_tri6):
            for xi, eta in [(0.0, 0.0), (0.3, -0.2), (0.4, 0.4)]:
                d = fn(xi, eta)
                self.assertAlmostEqual(d[0].sum(), 0.0, places=12)
                self.assertAlmostEqual(d[1].sum(), 0.0, places=12)


def _make_element(cls, nat_coords, mat, thickness=1.0):
    nodes = [Node(i + 1, list(c)) for i, c in enumerate(nat_coords)]
    for k, nd in enumerate(nodes):
        nd.add_dof('ux'); nd.add_dof('uy')
        nd.dofs['ux'] = 2 * k
        nd.dofs['uy'] = 2 * k + 1
    return cls(1, nodes, mat, thickness=thickness), nodes


def _u_field_quadratic(coords):
    """U construido a partir de u_x = x², u_y = 0 evaluado en cada nodo."""
    n = len(coords)
    U = np.zeros(2 * n)
    for k, c in enumerate(coords):
        x = c[0]
        U[2 * k] = x * x
        U[2 * k + 1] = 0.0
    return U


class TestQuadraticPatchTest(unittest.TestCase):
    """ε_xx(x) = 2x exacto en todos los Gauss para los 3 elementos."""

    def setUp(self):
        self.mat = _Elastic2D()

    def _check_patch(self, cls, coords):
        elem, _ = _make_element(cls, coords, self.mat)
        U = _u_field_quadratic(coords)
        gs = elem.compute_gauss_state(U)
        for k, p in enumerate(gs['points_global']):
            x = p[0]
            self.assertAlmostEqual(gs['strain'][k, 0], 2 * x, places=10,
                                   msg=f"{cls.__name__} ε_xx en Gauss {k}")
            self.assertAlmostEqual(gs['strain'][k, 1], 0.0, places=10)
            self.assertAlmostEqual(gs['strain'][k, 2], 0.0, places=10)

    def test_quad8(self):
        coords = [(0, 0), (2, 0), (2, 2), (0, 2),
                  (1, 0), (2, 1), (1, 2), (0, 1)]
        self._check_patch(Quad8, coords)

    def test_quad9(self):
        coords = [(0, 0), (2, 0), (2, 2), (0, 2),
                  (1, 0), (2, 1), (1, 2), (0, 1),
                  (1, 1)]
        self._check_patch(Quad9, coords)

    def test_tri6(self):
        coords = [(0, 0), (2, 0), (0, 2),
                  (1, 0), (1, 1), (0, 1)]
        self._check_patch(Tri6, coords)


class TestConsistentLoads(unittest.TestCase):
    def setUp(self):
        self.mat = _Elastic2D()

    def test_quad8_body_load_sum(self):
        coords = [(0, 0), (3, 0), (3, 2), (0, 2),
                  (1.5, 0), (3, 1), (1.5, 2), (0, 1)]
        elem, _ = _make_element(Quad8, coords, self.mat, thickness=0.5)
        b = np.array([0.0, -10.0])
        f = elem.compute_body_load(b)
        A = 6.0
        self.assertAlmostEqual(f[1::2].sum(), b[1] * A * 0.5, places=10)
        self.assertAlmostEqual(f[0::2].sum(), 0.0, places=10)

    def test_quad8_edge_traction_165(self):
        # Borde 0: nodos (0, 4, 1), longitud 3
        coords = [(0, 0), (3, 0), (3, 2), (0, 2),
                  (1.5, 0), (3, 1), (1.5, 2), (0, 1)]
        elem, _ = _make_element(Quad8, coords, self.mat, thickness=2.0)
        t = np.array([5.0, 0.0])
        f = elem.compute_edge_traction(0, t)
        L = 3.0
        # Σ = t · L · thickness
        self.assertAlmostEqual(f[0::2].sum(), t[0] * L * 2.0, places=12)
        # Vértice n0 (idx 0) = L/6 · t · thick; medio n4 (idx 4) = 4L/6 · t · thick
        self.assertAlmostEqual(f[0],     (L / 6.0) * t[0] * 2.0, places=12)
        self.assertAlmostEqual(f[2 * 4], (4 * L / 6.0) * t[0] * 2.0, places=12)
        self.assertAlmostEqual(f[2 * 1], (L / 6.0) * t[0] * 2.0, places=12)
        # Otros nodos: cero
        for j in (2, 3, 5, 6, 7):
            self.assertAlmostEqual(f[2 * j], 0.0, places=12)

    def test_tri6_body_load_sum(self):
        coords = [(0, 0), (3, 0), (0, 4),
                  (1.5, 0), (1.5, 2), (0, 2)]
        elem, _ = _make_element(Tri6, coords, self.mat, thickness=1.0)
        b = np.array([0.0, -2.0])
        f = elem.compute_body_load(b)
        A = 0.5 * 3 * 4
        self.assertAlmostEqual(f[1::2].sum(), b[1] * A, places=10)

    def test_tri6_edge_traction_165(self):
        coords = [(0, 0), (3, 0), (0, 4),
                  (1.5, 0), (1.5, 2), (0, 2)]
        elem, _ = _make_element(Tri6, coords, self.mat, thickness=1.0)
        # Borde 0: nodos (0, 3, 1), longitud 3
        t = np.array([1.0, 0.0])
        f = elem.compute_edge_traction(0, t)
        L = 3.0
        self.assertAlmostEqual(f[0::2].sum(), t[0] * L, places=12)
        self.assertAlmostEqual(f[2 * 3], (4 * L / 6.0) * t[0], places=12)


class TestQuad8UniaxialPatch(unittest.TestCase):
    """Patch físico: un Quad8 con tracción de borde reproduce σ_xx = p exacto."""

    def test_uniaxial(self):
        E = 1000.0; nu = 0.3; p = 5.0; L = 2.0; H = 1.0; t = 1.0
        mat = _Elastic2D(E, nu)
        coords = [(0, 0), (L, 0), (L, H), (0, H),
                  (L/2, 0), (L, H/2), (L/2, H), (0, H/2)]
        elem, nodes = _make_element(Quad8, coords, mat, thickness=t)

        K_e = elem.compute_global_stiffness()
        F = elem.compute_edge_traction(1, np.array([p, 0.0]))  # borde derecho

        # Restricciones: ux = 0 en los nodos del borde izquierdo (n0, n7, n3);
        # uy = 0 sólo en n0 para fijar el RBM en y.
        # DOFs (16): 2*idx + 0/1
        fixed = [2*0, 2*0+1, 2*7, 2*3]
        free = [k for k in range(16) if k not in fixed]
        K_ff = K_e[np.ix_(free, free)]
        u_f = np.linalg.solve(K_ff, F[free])
        U = np.zeros(16)
        U[free] = u_f

        gs = elem.compute_gauss_state(U)
        for k in range(gs['stress'].shape[0]):
            self.assertAlmostEqual(gs['stress'][k, 0], p,    places=8)
            self.assertAlmostEqual(gs['stress'][k, 1], 0.0,  places=8)
            self.assertAlmostEqual(gs['stress'][k, 2], 0.0,  places=8)


if __name__ == '__main__':
    unittest.main()
