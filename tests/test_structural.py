import unittest
import numpy as np
import sys
import os

# Ensure the parent directory is in the path to import fenix
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.node import Node
from fenix.core.material import Material
from fenix.elements.structural import Truss2D, Truss3D

class DummyMaterial(Material):
    """Material elástico lineal simple para pruebas."""
    def __init__(self, E: float):
        self.E = E

    def compute_state(self, strain: float, state_var=None):
        stress = self.E * strain
        return stress, self.E, state_var

class TestTruss2D(unittest.TestCase):
    def setUp(self):
        self.node1 = Node(1, [0.0, 0.0])
        self.node2 = Node(2, [3.0, 4.0]) # Longitud = 5.0, c = 0.6, s = 0.8
        
        # Simulamos que los nodos ya tienen sus grados de libertad asignados con un índice
        self.node1.add_dof('ux')
        self.node1.add_dof('uy')
        self.node2.add_dof('ux')
        self.node2.add_dof('uy')
        
        self.node1.dofs['ux'] = 0
        self.node1.dofs['uy'] = 1
        self.node2.dofs['ux'] = 2
        self.node2.dofs['uy'] = 3
        
        self.material = DummyMaterial(E=1000.0)
        self.A = 2.0
        self.element = Truss2D(1, [self.node1, self.node2], self.material, self.A)

    def test_initialization(self):
        self.assertEqual(self.element.L0, 5.0)
        self.assertAlmostEqual(self.element.c, 0.6)
        self.assertAlmostEqual(self.element.s, 0.8)

    def test_compute_global_stiffness(self):
        K_e = self.element.compute_global_stiffness()
        
        # coef = E * A / L = 1000 * 2 / 5 = 400
        # K_e[0, 0] = coef * c^2 = 400 * 0.36 = 144
        # K_e[0, 1] = coef * c * s = 400 * 0.48 = 192
        self.assertAlmostEqual(K_e[0, 0], 144.0)
        self.assertAlmostEqual(K_e[0, 1], 192.0)
        self.assertAlmostEqual(K_e[2, 2], 144.0)
        self.assertAlmostEqual(K_e[3, 3], 256.0) # 400 * 0.64
        
        # Propiedades de la matriz
        self.assertTrue(np.allclose(K_e, K_e.T)) # Simétrica

    def test_compute_element_state(self):
        # Desplazamiento que causa un alargamiento de 0.5 unidades
        # vector_dir = [0.6, 0.8], move node2 along vector_dir by 0.5
        u_e = np.array([0.0, 0.0, 0.6 * 0.5, 0.8 * 0.5])
        
        K_e, F_int_e = self.element.compute_element_state(u_e)
        
        # epsilon = du / L = 0.5 / 5.0 = 0.1
        # sigma = E * epsilon = 1000 * 0.1 = 100
        # N = A * sigma = 2 * 100 = 200
        # F_int_e = N * [-c, -s, c, s]
        self.assertAlmostEqual(F_int_e[2], 200.0 * 0.6)
        self.assertAlmostEqual(F_int_e[3], 200.0 * 0.8)
        self.assertAlmostEqual(F_int_e[0], -200.0 * 0.6)

    def test_compute_internal_forces(self):
        U_global = np.array([0.0, 0.0, 0.3, 0.4]) # u_e equivalent
        results = self.element.compute_internal_forces(U_global)
        
        # epsilon = 0.5 / 5.0 = 0.1
        self.assertAlmostEqual(results['strain'], 0.1)
        self.assertAlmostEqual(results['stress'], 100.0)
        self.assertAlmostEqual(results['axial_force'], 200.0)

class TestTruss3D(unittest.TestCase):
    def setUp(self):
        self.node1 = Node(1, [0.0, 0.0, 0.0])
        self.node2 = Node(2, [3.0, 4.0, 12.0]) # Longitud = sqrt(9 + 16 + 144) = 13.0
        
        # Simulamos grados de libertad
        for node in [self.node1, self.node2]:
            node.add_dof('ux')
            node.add_dof('uy')
            node.add_dof('uz')
            
        self.node1.dofs['ux'] = 0
        self.node1.dofs['uy'] = 1
        self.node1.dofs['uz'] = 2
        self.node2.dofs['ux'] = 3
        self.node2.dofs['uy'] = 4
        self.node2.dofs['uz'] = 5
        
        self.material = DummyMaterial(E=2600.0) # E/L = 200
        self.A = 1.0
        self.element = Truss3D(1, [self.node1, self.node2], self.material, self.A)

    def test_initialization(self):
        self.assertEqual(self.element.L0, 13.0)
        self.assertAlmostEqual(self.element.cx, 3.0/13.0)
        self.assertAlmostEqual(self.element.cz, 12.0/13.0)
        
    def test_compute_global_stiffness(self):
        K_e = self.element.compute_global_stiffness()
        
        # coef = E * A / L = 2600 * 1 / 13 = 200
        # K_e[0, 0] = coef * cx^2 = 200 * (9/169) = 1800/169 = 10.6508...
        self.assertAlmostEqual(K_e[0, 0], 200 * (3/13)**2)
        self.assertAlmostEqual(K_e[2, 2], 200 * (12/13)**2)
        
        self.assertEqual(K_e.shape, (6, 6))
        self.assertTrue(np.allclose(K_e, K_e.T)) # Simétrica

    def test_compute_internal_forces(self):
        # Desplazamiento que alarga el elemento
        # dx, dy, dz a lo largo del elemento por 1.3 unidades
        U_global = np.array([0.0, 0.0, 0.0, 0.3, 0.4, 1.2]) 
        
        results = self.element.compute_internal_forces(U_global)
        
        # epsilon = dl / L = 1.3 / 13.0 = 0.1
        self.assertAlmostEqual(results['strain'], 0.1)
        self.assertAlmostEqual(results['stress'], 260.0) # E * eps = 2600 * 0.1
        self.assertAlmostEqual(results['axial_force'], 260.0) # A * stress = 1 * 260

if __name__ == '__main__':
    unittest.main()
