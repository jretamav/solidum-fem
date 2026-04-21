"""Tests de los elementos de cable.

Materializan los criterios de `acceptance` declarados en
docs/specs/Cable2DCorot.md.
"""
import math
import unittest

import numpy as np

import fenix  # dispara autodiscover
from fenix.core.node import Node
from fenix.elements.cable import Cable2DCorot
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


if __name__ == '__main__':
    unittest.main()
