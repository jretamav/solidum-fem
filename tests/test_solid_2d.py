import unittest
import numpy as np
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.node import Node
from solidum.core.material import Material
from solidum.elements.solid_2d import Quad4, Tri3

class DummyMaterial2D(Material):
    """Material elástico lineal simple para pruebas en 2D."""
    def __init__(self, E: float, nu: float):
        self.E = E
        self.nu = nu
        # Matriz constitutiva para tensión plana
        fac = E / (1.0 - nu**2)
        self.C = fac * np.array([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0]
        ])

    def compute_state(self, strain: np.ndarray, state_var=None):
        stress = self.C @ strain
        return stress, self.C, state_var

class TestQuad4(unittest.TestCase):
    def setUp(self):
        # Elemento cuadrado de 2x2 centrado en el origen
        self.node1 = Node(1, [-1.0, -1.0])
        self.node2 = Node(2, [ 1.0, -1.0])
        self.node3 = Node(3, [ 1.0,  1.0])
        self.node4 = Node(4, [-1.0,  1.0])
        
        for i, node in enumerate([self.node1, self.node2, self.node3, self.node4]):
            node.add_dof('ux')
            node.add_dof('uy')
            node.dofs['ux'] = 2*i
            node.dofs['uy'] = 2*i+1
            
        self.material = DummyMaterial2D(E=1000.0, nu=0.3)
        self.element = Quad4(1, [self.node1, self.node2, self.node3, self.node4], self.material, thickness=1.0)

    def test_initialization(self):
        # Solo verificamos atributos básicos de inicialización
        self.assertEqual(self.element.thickness, 1.0)
        self.assertEqual(len(self.element.nodes), 4)

    def test_compute_global_stiffness(self):
        K_e = self.element.compute_global_stiffness()
        self.assertEqual(K_e.shape, (8, 8))
        self.assertTrue(np.allclose(K_e, K_e.T)) # Simétrica

    def test_compute_internal_forces(self):
        # Aplicamos una deformación uniforme eps_x = 0.1, eps_y = 0, gamma_xy = 0
        # u(x,y) = 0.1 * x
        U_global = np.array([
            -0.1, 0.0,  # Nodo 1 (-1, -1)
             0.1, 0.0,  # Nodo 2 (1, -1)
             0.1, 0.0,  # Nodo 3 (1, 1)
            -0.1, 0.0   # Nodo 4 (-1, 1)
        ])
        
        results = self.element.compute_internal_forces(U_global)
        
        expected_strain = np.array([0.1, 0.0, 0.0])
        expected_stress = self.material.C @ expected_strain
        
        self.assertTrue(np.allclose(results['strain'], expected_strain))
        self.assertTrue(np.allclose(results['stress'], expected_stress))

class TestTri3(unittest.TestCase):
    def setUp(self):
        # Triángulo rectángulo en el origen de base b=2, altura h=2
        self.node1 = Node(1, [0.0, 0.0])
        self.node2 = Node(2, [2.0, 0.0])
        self.node3 = Node(3, [0.0, 2.0])
        
        for i, node in enumerate([self.node1, self.node2, self.node3]):
            node.add_dof('ux')
            node.add_dof('uy')
            node.dofs['ux'] = 2*i
            node.dofs['uy'] = 2*i+1
            
        self.material = DummyMaterial2D(E=1000.0, nu=0.3)
        self.element = Tri3(1, [self.node1, self.node2, self.node3], self.material, thickness=1.0)

    def test_initialization(self):
        self.assertEqual(self.element.thickness, 1.0)
        self.assertEqual(len(self.element.nodes), 3)

    def test_compute_global_stiffness(self):
        K_e = self.element.compute_global_stiffness()
        self.assertEqual(K_e.shape, (6, 6))
        self.assertTrue(np.allclose(K_e, K_e.T))

    def test_compute_internal_forces(self):
        # Deformación uniforme eps_y = 0.05
        # u(x,y) = 0, v(x,y) = 0.05 * y
        U_global = np.array([
            0.0, 0.0,    # Nodo 1 (0, 0)
            0.0, 0.0,    # Nodo 2 (2, 0)
            0.0, 0.1     # Nodo 3 (0, 2) -> 0.05 * 2 = 0.1
        ])
        
        results = self.element.compute_internal_forces(U_global)
        
        expected_strain = np.array([0.0, 0.05, 0.0])
        expected_stress = self.material.C @ expected_strain
        
        self.assertTrue(np.allclose(results['strain'], expected_strain))
        self.assertTrue(np.allclose(results['stress'], expected_stress))

if __name__ == '__main__':
    unittest.main()
