"""Tests de la API por puntos de Gauss en Quad4 y Tri3.

Cubre:
- Forma y semántica del dict devuelto por ``compute_gauss_state``.
- Coordenadas globales: el mapeo isoparamétrico debe llevar los puntos
  naturales de Gauss al lugar correcto en el sistema global.
- Campo lineal ⇒ ε constante e idéntica en todos los puntos de Gauss.
- ``compute_internal_forces`` sigue devolviendo el promedio.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.material import Material
from solidum.core.node import Node
from solidum.elements.solid_2d import Quad4, Tri3


class _Elastic2D(Material):
    STRAIN_DIM = 3

    def __init__(self, E=1000.0, nu=0.3, density: float = 0.0):
        self.E = E
        self.nu = nu
        fac = E / (1.0 - nu * nu)
        self.C = fac * np.array([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0],
        ])

    def compute_state(self, strain, sv=None):
        return self.C @ strain, self.C, sv


def _quad4(coords, mat):
    nodes = [Node(i + 1, list(c)) for i, c in enumerate(coords)]
    for k, nd in enumerate(nodes):
        nd.add_dof('ux'); nd.add_dof('uy')
        nd.dofs['ux'] = 2 * k
        nd.dofs['uy'] = 2 * k + 1
    return Quad4(1, nodes, mat), nodes


def _tri3(coords, mat):
    nodes = [Node(i + 1, list(c)) for i, c in enumerate(coords)]
    for k, nd in enumerate(nodes):
        nd.add_dof('ux'); nd.add_dof('uy')
        nd.dofs['ux'] = 2 * k
        nd.dofs['uy'] = 2 * k + 1
    return Tri3(1, nodes, mat), nodes


class TestQuad4GaussState(unittest.TestCase):
    def setUp(self):
        self.mat = _Elastic2D()

    def test_shape_and_global_points(self):
        # Cuadrado [0,2] x [0,2]: los puntos de Gauss caen en (1±1/√3, 1±1/√3)
        coords = [(0, 0), (2, 0), (2, 2), (0, 2)]
        elem, _ = _quad4(coords, self.mat)
        gs = elem.compute_gauss_state(np.zeros(8))

        self.assertEqual(gs['points_natural'].shape, (4, 2))
        self.assertEqual(gs['points_global'].shape, (4, 2))
        self.assertEqual(gs['strain'].shape, (4, 3))
        self.assertEqual(gs['stress'].shape, (4, 3))

        d = 1.0 / np.sqrt(3.0)
        expected_global = np.array([
            [1 - d, 1 - d],
            [1 + d, 1 - d],
            [1 + d, 1 + d],
            [1 - d, 1 + d],
        ])
        self.assertTrue(np.allclose(gs['points_global'], expected_global, atol=1e-12))

    def test_uniform_strain_field(self):
        # u(x,y) = 1e-3·x, v=0  ⇒  ε_xx = 1e-3, resto 0  en TODOS los Gauss
        coords = [(0, 0), (2, 0), (2, 2), (0, 2)]
        elem, _ = _quad4(coords, self.mat)
        U = np.array([0.0, 0.0,  2e-3, 0.0,  2e-3, 0.0,  0.0, 0.0])
        gs = elem.compute_gauss_state(U)
        for i in range(4):
            self.assertAlmostEqual(gs['strain'][i, 0], 1e-3, places=12)
            self.assertAlmostEqual(gs['strain'][i, 1], 0.0,  places=12)
            self.assertAlmostEqual(gs['strain'][i, 2], 0.0,  places=12)

    def test_compute_internal_forces_matches_average(self):
        coords = [(0, 0), (2, 0), (2, 2), (0, 2)]
        elem, _ = _quad4(coords, self.mat)
        U = np.array([0.0, 0.0,  2e-3, 0.0,  2e-3, 0.0,  0.0, 0.0])
        gs = elem.compute_gauss_state(U)
        avg = elem.compute_internal_forces(U)
        self.assertTrue(np.allclose(avg['strain'], gs['strain'].mean(axis=0)))
        self.assertTrue(np.allclose(avg['stress'], gs['stress'].mean(axis=0)))


class TestTri3GaussState(unittest.TestCase):
    def setUp(self):
        self.mat = _Elastic2D()

    def test_shape_and_centroid(self):
        coords = [(0, 0), (3, 0), (0, 3)]
        elem, _ = _tri3(coords, self.mat)
        gs = elem.compute_gauss_state(np.zeros(6))
        self.assertEqual(gs['points_natural'].shape, (1, 2))
        self.assertEqual(gs['points_global'].shape, (1, 2))
        self.assertTrue(np.allclose(gs['points_global'][0], [1.0, 1.0]))

    def test_constant_strain(self):
        coords = [(0, 0), (2, 0), (0, 2)]
        elem, _ = _tri3(coords, self.mat)
        # ε_yy = 0.05, resto 0  →  v(x,y) = 0.05·y
        U = np.array([0.0, 0.0,  0.0, 0.0,  0.0, 0.10])
        gs = elem.compute_gauss_state(U)
        self.assertAlmostEqual(gs['strain'][0, 0], 0.0,  places=12)
        self.assertAlmostEqual(gs['strain'][0, 1], 0.05, places=12)
        self.assertAlmostEqual(gs['strain'][0, 2], 0.0,  places=12)


if __name__ == '__main__':
    unittest.main()
