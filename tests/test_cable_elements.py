"""Tests de los elementos de cable.

Materializan los criterios de `acceptance` declarados en
docs/specs/Cable2DCorot.md.
"""
import math
import unittest

import numpy as np

import fenix  # dispara autodiscover
from fenix.core.node import Node
from fenix.elements.cable import Cable2DCorot, Cable3DCorot
from fenix.elements.truss import Truss2D, Truss3D, Truss2DCorot, Truss3DCorot
from fenix.materials.cable_1d import CableMaterial1D
from fenix.materials.elastic import Elastic1D


class TestCable2DCorot(unittest.TestCase):

    def _build(self, coords1, coords2, material, A=1.0):
        n1 = Node(1, list(coords1))
        n2 = Node(2, list(coords2))
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy')
            node.dofs['ux'] = 2 * k
            node.dofs['uy'] = 2 * k + 1
        return Cable2DCorot(1, [n1, n2], material, A)

    def test_acceptance_cable_tensado_como_barra_corot(self):
        """Criterio 1: tensado con material lineal coincide con barra corot."""
        E, A = 1000.0, 1.0
        mat = Elastic1D(E=E)
        cable = self._build([0.0, 0.0], [1.0, 0.0], mat, A=A)
        u_e = np.array([0.0, 0.0, 1e-3, 0.0])  # ε = 1e-3 > 0 → tensado

        K_T, F_int = cable.compute_element_state(u_e)

        # Valores esperados analíticos
        l = 1.0 + 1e-3
        N_expected = E * 1e-3 * A  # σ·A
        d = np.array([-1.0, 0.0, 1.0, 0.0])
        n = np.array([0.0, 1.0, 0.0, -1.0])
        K_M_expected = (E * A / 1.0) * np.outer(d, d)
        K_G_expected = (N_expected / l) * np.outer(n, n)

        self.assertTrue(np.allclose(K_T, K_M_expected + K_G_expected, atol=1e-10))
        self.assertTrue(np.allclose(F_int, N_expected * d, atol=1e-10))

    def test_acceptance_cable_destensado_aportacion_nula(self):
        """Criterio 2: con material unilateral + ε < 0 ⇒ K_T = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0], [1.0, 0.0], mat, A=1.0)
        u_e = np.array([0.0, 0.0, -1e-3, 0.0])  # ε < 0 → destensado

        K_T, F_int = cable.compute_element_state(u_e)

        self.assertTrue(np.allclose(K_T, 0.0, atol=1e-12))
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-12))
        self.assertEqual(cable.state.stresses_trial[0], 0.0)

    def test_acceptance_rotacion_rigida(self):
        """Criterio 3: rotación rígida ⇒ ε = 0 ⇒ σ = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0], [1.0, 0.0], mat, A=1.0)
        theta = math.pi / 4
        c, s = math.cos(theta), math.sin(theta)
        u_e = np.array([0.0, 0.0, c - 1.0, s])  # nodo 2 en (c, s)

        _, F_int = cable.compute_element_state(u_e)

        self.assertAlmostEqual(cable.state.stresses_trial[0], 0.0, places=12)
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-10))

    def test_acceptance_cruce_por_cero(self):
        """Criterio 4: ε = 0 exacto ⇒ σ = 0, E_t = 0, K_T = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0], [1.0, 0.0], mat, A=1.0)
        u_e = np.zeros(4)  # ε = 0

        K_T, F_int = cable.compute_element_state(u_e)

        self.assertEqual(cable.state.stresses_trial[0], 0.0)
        self.assertTrue(np.allclose(K_T, 0.0, atol=1e-12))
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-12))

    def test_registro_en_registry(self):
        """El elemento queda registrado vía autodiscover."""
        from fenix.registry import ElementRegistry
        self.assertIn('Cable2DCorot', ElementRegistry._items)


class TestCable3DCorot(unittest.TestCase):

    def _build(self, coords1, coords2, material, A=1.0, cls=Cable3DCorot):
        n1 = Node(1, list(coords1))
        n2 = Node(2, list(coords2))
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy'); node.add_dof('uz')
            node.dofs['ux'] = 3 * k
            node.dofs['uy'] = 3 * k + 1
            node.dofs['uz'] = 3 * k + 2
        return cls(1, [n1, n2], material, A)

    def test_acceptance_cable_tensado_como_barra_corot_3d(self):
        """Criterio 1: tensado con material lineal coincide con Truss3DCorot."""
        E, A = 1000.0, 1.0
        cable = self._build([0.0, 0.0, 0.0], [3.0, 4.0, 12.0], Elastic1D(E=E), A=A)
        truss = self._build([0.0, 0.0, 0.0], [3.0, 4.0, 12.0], Elastic1D(E=E), A=A,
                            cls=Truss3DCorot)
        axis = np.array([3.0, 4.0, 12.0]) / 13.0
        u_tiny = np.concatenate([np.zeros(3), 1e-6 * axis])

        K_c, F_c = cable.compute_element_state(u_tiny)
        K_t, F_t = truss.compute_element_state(u_tiny)

        self.assertTrue(np.allclose(K_c, K_t, rtol=1e-4))
        self.assertTrue(np.allclose(F_c, F_t, rtol=1e-4))

    def test_acceptance_cable_destensado_aportacion_nula(self):
        """Criterio 2: con material unilateral + ε < 0 ⇒ K_T = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], mat, A=1.0)
        u_e = np.array([0.0, 0.0, 0.0, -1e-3, 0.0, 0.0])  # ε < 0

        K_T, F_int = cable.compute_element_state(u_e)

        self.assertTrue(np.allclose(K_T, 0.0, atol=1e-12))
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-12))
        self.assertEqual(cable.state.stresses_trial[0], 0.0)

    def test_acceptance_rotacion_rigida_3d(self):
        """Criterio 3: rotación rígida 3D ⇒ ε = 0 ⇒ σ = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], mat, A=1.0)
        theta = math.pi / 4
        c, s = math.cos(theta), math.sin(theta)
        # R_y(theta) @ (1,0,0) = (c, 0, -s); u_2 = final - inicial
        u_e = np.array([0.0, 0.0, 0.0, c - 1.0, 0.0, -s])

        _, F_int = cable.compute_element_state(u_e)

        self.assertAlmostEqual(cable.state.stresses_trial[0], 0.0, places=12)
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-10))

    def test_acceptance_cruce_por_cero(self):
        """Criterio 4: u_e = 0 ⇒ σ = 0, E_t = 0, K_T = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], mat, A=1.0)
        u_e = np.zeros(6)

        K_T, F_int = cable.compute_element_state(u_e)

        self.assertEqual(cable.state.stresses_trial[0], 0.0)
        self.assertTrue(np.allclose(K_T, 0.0, atol=1e-12))
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-12))

    def test_registro_en_registry(self):
        from fenix.registry import ElementRegistry
        self.assertIn('Cable3DCorot', ElementRegistry._items)


class TestUnilateralMaterialValidation(unittest.TestCase):
    """Validación declarativa material.IS_UNILATERAL vs element.ACCEPTS_UNILATERAL.

    Un material cuya rigidez tangente puede colapsar a cero (cables, gaps)
    solo es admisible en elementos preparados para manejarlo. La base
    Element atrapa la combinación incompatible en construcción.
    """

    def _make_nodes_2d(self):
        n1 = Node(1, [0.0, 0.0]); n2 = Node(2, [1.0, 0.0])
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy')
            node.dofs['ux'] = 2 * k
            node.dofs['uy'] = 2 * k + 1
        return n1, n2

    def _make_nodes_3d(self):
        n1 = Node(1, [0.0, 0.0, 0.0]); n2 = Node(2, [1.0, 0.0, 0.0])
        for k, node in enumerate((n1, n2)):
            for d in ('ux', 'uy', 'uz'):
                node.add_dof(d)
            node.dofs['ux'] = 3 * k
            node.dofs['uy'] = 3 * k + 1
            node.dofs['uz'] = 3 * k + 2
        return n1, n2

    def test_truss2d_rechaza_cable_material(self):
        n1, n2 = self._make_nodes_2d()
        mat = CableMaterial1D(E=1000.0)
        with self.assertRaisesRegex(ValueError, "unilateral"):
            Truss2D(1, [n1, n2], mat, A=1.0)

    def test_truss3d_rechaza_cable_material(self):
        n1, n2 = self._make_nodes_3d()
        mat = CableMaterial1D(E=1000.0)
        with self.assertRaisesRegex(ValueError, "unilateral"):
            Truss3D(1, [n1, n2], mat, A=1.0)

    def test_truss_corot_tambien_rechaza_cable_material(self):
        """Truss*Corot tampoco lo acepta — el material unilateral se reserva
        a Cable*Corot, donde la formulación está diseñada explícitamente
        para él."""
        n1, n2 = self._make_nodes_2d()
        with self.assertRaisesRegex(ValueError, "unilateral"):
            Truss2DCorot(1, [n1, n2], CableMaterial1D(E=1000.0), A=1.0)
        n3, n4 = self._make_nodes_3d()
        with self.assertRaisesRegex(ValueError, "unilateral"):
            Truss3DCorot(1, [n3, n4], CableMaterial1D(E=1000.0), A=1.0)

    def test_cable_corot_acepta_cable_material(self):
        """Combinación canónica — no debe lanzar."""
        n1, n2 = self._make_nodes_2d()
        Cable2DCorot(1, [n1, n2], CableMaterial1D(E=1000.0), A=1.0)
        n3, n4 = self._make_nodes_3d()
        Cable3DCorot(2, [n3, n4], CableMaterial1D(E=1000.0), A=1.0)

    def test_truss_sigue_aceptando_materiales_bilaterales(self):
        """Regresión: la nueva validación no rompe el contrato con Elastic1D."""
        n1, n2 = self._make_nodes_2d()
        Truss2D(1, [n1, n2], Elastic1D(E=1000.0), A=1.0)
        n3, n4 = self._make_nodes_3d()
        Truss3D(1, [n3, n4], Elastic1D(E=1000.0), A=1.0)


if __name__ == '__main__':
    unittest.main()
