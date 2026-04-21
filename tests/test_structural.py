import math
import unittest
import numpy as np
import sys
import os

# Ensure the parent directory is in the path to import fenix
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.node import Node
from fenix.core.material import Material
from fenix.elements.structural import Truss2D, Truss2DCorot, Truss3D, Truss3DCorot

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

class TestTruss2DCorot(unittest.TestCase):
    """Tests que materializan los tres criterios de aceptación de la spec
    docs/specs/Truss2DCorot.md.
    """

    def _build(self, coords1, coords2, E=1000.0, A=2.0):
        n1 = Node(1, list(coords1))
        n2 = Node(2, list(coords2))
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy')
            node.dofs['ux'] = 2 * k
            node.dofs['uy'] = 2 * k + 1
        mat = DummyMaterial(E=E)
        return Truss2DCorot(1, [n1, n2], mat, A)

    def test_acceptance_linear_limit_matches_truss2d(self):
        """Criterio 1: en régimen infinitesimal reproduce Truss2D lineal."""
        corot = self._build([0.0, 0.0], [3.0, 4.0])  # L0=5, c=0.6, s=0.8
        linear = Truss2D(2, corot.nodes, corot.material, corot.A)
        u_tiny = np.array([0.0, 0.0, 0.6 * 1e-6, 0.8 * 1e-6])

        K_c, F_c = corot.compute_element_state(u_tiny)
        K_l, F_l = linear.compute_element_state(u_tiny)

        self.assertTrue(np.allclose(K_c, K_l, rtol=1e-4))
        self.assertTrue(np.allclose(F_c, F_l, rtol=1e-4))

    def test_acceptance_rigid_body_rotation(self):
        """Criterio 2: rotación rígida ⇒ σ=0 y F_int=0."""
        corot = self._build([0.0, 0.0], [1.0, 0.0])  # L0=1, horizontal
        theta = math.pi / 4  # 45°
        c, s = math.cos(theta), math.sin(theta)
        # Nodos rotados rígidamente (nodo 1 en origen, nodo 2 en (c, s))
        u_e = np.array([0.0, 0.0, c - 1.0, s - 0.0])

        _, F_int = corot.compute_element_state(u_e)

        self.assertAlmostEqual(corot.state.stresses_trial[0], 0.0, places=12)
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-12))

    def test_acceptance_geometric_stiffness_under_traction(self):
        """Criterio 3: K_G = (N/l) · n·nᵀ con n perpendicular al eje."""
        corot = self._build([0.0, 0.0], [1.0, 0.0], E=1000.0, A=1.0)
        # Estiramos axialmente: u_x2 = 0.001 ⇒ ε = 0.001 ⇒ σ = 1 ⇒ N = 1
        u_e = np.array([0.0, 0.0, 1e-3, 0.0])

        K_T, _ = corot.compute_element_state(u_e)

        # Geometría corriente tras estirar
        l = 1.0 + 1e-3
        c_t, s_t = 1.0, 0.0  # barra sigue horizontal
        N = 1.0
        E_t, A, L0 = 1000.0, 1.0, 1.0
        d = np.array([-c_t, -s_t, c_t, s_t])
        n = np.array([-s_t, c_t, s_t, -c_t])

        K_M_expected = (E_t * A / L0) * np.outer(d, d)
        K_G_expected = (N / l) * np.outer(n, n)

        self.assertTrue(np.allclose(K_T, K_M_expected + K_G_expected, atol=1e-10))
        # Verificación aislada de K_G: K_T - K_M debe ser rango-1 en n
        K_G_actual = K_T - K_M_expected
        # K_G @ n = (N/l) * n * (n@n) = (N/l) * n * 2
        self.assertTrue(np.allclose(K_G_actual @ n, (2 * N / l) * n, atol=1e-10))


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

class TestTruss3DCorot(unittest.TestCase):
    """Tests que materializan los tres criterios de aceptación de
    docs/specs/Truss3DCorot.md.
    """

    def _build(self, coords1, coords2, E=1000.0, A=1.0):
        n1 = Node(1, list(coords1))
        n2 = Node(2, list(coords2))
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy'); node.add_dof('uz')
            node.dofs['ux'] = 3 * k
            node.dofs['uy'] = 3 * k + 1
            node.dofs['uz'] = 3 * k + 2
        mat = DummyMaterial(E=E)
        return Truss3DCorot(1, [n1, n2], mat, A)

    def test_acceptance_linear_limit_matches_truss3d(self):
        """Criterio 1: en régimen infinitesimal reproduce Truss3D lineal."""
        corot = self._build([0.0, 0.0, 0.0], [3.0, 4.0, 12.0])  # L0=13
        linear = Truss3D(2, corot.nodes, corot.material, corot.A)
        # Desplazamiento axial pequeño del nodo 2 en la dirección del eje
        axis = np.array([3.0, 4.0, 12.0]) / 13.0
        u_tiny = np.concatenate([np.zeros(3), 1e-6 * axis])

        K_c, F_c = corot.compute_element_state(u_tiny)
        K_l, F_l = linear.compute_element_state(u_tiny)

        self.assertTrue(np.allclose(K_c, K_l, rtol=1e-4))
        self.assertTrue(np.allclose(F_c, F_l, rtol=1e-4))

    def test_acceptance_rigid_body_rotation(self):
        """Criterio 2: rotación rígida 3D ⇒ σ=0 y F_int=0."""
        corot = self._build([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])  # barra sobre eje x, L0=1
        # Rotación rígida 90° alrededor del eje z: (1,0,0) -> (0,1,0)
        # Nodo 1 queda en origen, nodo 2 va a (0,1,0)
        u_e = np.array([0.0, 0.0, 0.0, -1.0, 1.0, 0.0])

        _, F_int = corot.compute_element_state(u_e)

        self.assertAlmostEqual(corot.state.stresses_trial[0], 0.0, places=12)
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-12))

        # Segunda rotación compuesta: 45° alrededor del eje y sobre la barra ya rotada.
        # Posición esperada del nodo 2: R_y(45°) @ (0,1,0) = (0,1,0) (eje de rotación)
        # — degenerado, probemos rotando la barra original 45° en plano x-z:
        corot2 = self._build([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        theta = math.pi / 4
        c, s = math.cos(theta), math.sin(theta)
        # R_y(theta) @ (1,0,0) = (c, 0, -s); desplazamiento = posición final − inicial
        u_e2 = np.array([0.0, 0.0, 0.0, c - 1.0, 0.0, -s])
        _, F2 = corot2.compute_element_state(u_e2)
        self.assertAlmostEqual(corot2.state.stresses_trial[0], 0.0, places=12)
        self.assertTrue(np.allclose(F2, 0.0, atol=1e-12))

    def test_acceptance_geometric_stiffness_transverse_plane(self):
        """Criterio 3: K_G rigidiza el plano perpendicular al eje con factor N/l.

        En 3D la perpendicular al eje es un plano; hay dos direcciones
        linealmente independientes que deben producir la misma respuesta.
        """
        corot = self._build([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], E=1000.0, A=1.0)
        # Estirar axialmente: u_x node 2 = 0.001 ⇒ ε = 0.001 ⇒ σ = 1 ⇒ N = 1
        u_e = np.array([0.0, 0.0, 0.0, 1e-3, 0.0, 0.0])
        K_T, _ = corot.compute_element_state(u_e)

        # Geometría corriente tras estirar
        l = 1.0 + 1e-3
        N = 1.0
        cx, cy, cz = 1.0, 0.0, 0.0  # eje sigue sobre x
        E_t, A, L0 = 1000.0, 1.0, 1.0
        d = np.array([-cx, -cy, -cz, cx, cy, cz])
        P = np.eye(3) - np.outer([cx, cy, cz], [cx, cy, cz])

        K_M_expected = (E_t * A / L0) * np.outer(d, d)
        K_G_expected = np.zeros((6, 6))
        K_G_expected[:3, :3] = P;  K_G_expected[3:, 3:] = P
        K_G_expected[:3, 3:] = -P; K_G_expected[3:, :3] = -P
        K_G_expected *= (N / l)

        self.assertTrue(np.allclose(K_T, K_M_expected + K_G_expected, atol=1e-10))

        # Autovalor sobre modos transversos puros (rotación alrededor del centro).
        # v_y, v_z son ortogonales al eje y ortogonales entre sí → dos direcciones
        # del plano perpendicular. K_G · v = 2·(N/l) · v (factor 2 por la norma √2
        # del vector de 6 componentes, análogo al caso 2D).
        K_G_actual = K_T - K_M_expected
        v_y = np.array([0.0, -0.5, 0.0, 0.0, 0.5, 0.0])
        v_z = np.array([0.0, 0.0, -0.5, 0.0, 0.0, 0.5])
        self.assertTrue(np.allclose(K_G_actual @ v_y, (2 * N / l) * v_y, atol=1e-10))
        self.assertTrue(np.allclose(K_G_actual @ v_z, (2 * N / l) * v_z, atol=1e-10))


if __name__ == '__main__':
    unittest.main()
