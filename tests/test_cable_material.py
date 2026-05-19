"""Tests del material de cable 1D.

Materializan los criterios de `acceptance` declarados en
docs/specs/CableMaterial1D.md.
"""
import unittest

import fenix  # dispara autodiscover → registra CableMaterial1D
from fenix.materials.cable_1d import CableMaterial1D


class TestCableMaterial1D(unittest.TestCase):

    def setUp(self):
        self.E = 210e9
        self.mat = CableMaterial1D(E=self.E)

    def test_acceptance_tension_pura(self):
        """ε > 0 ⇒ σ = E·ε, E_t = E."""
        eps = 1.0e-4
        sigma, Et, _ = self.mat.compute_state(eps)
        self.assertAlmostEqual(sigma, self.E * eps, delta=1e-12 * self.E * eps)
        self.assertEqual(Et, self.E)

    def test_acceptance_compresion_pura(self):
        """ε < 0 ⇒ σ = 0, E_t = 0."""
        sigma, Et, _ = self.mat.compute_state(-1.0e-4)
        self.assertEqual(sigma, 0.0)
        self.assertEqual(Et, 0.0)

    def test_acceptance_punto_de_corte(self):
        """ε = 0 ⇒ σ = 0, E_t = 0 (convención: régimen compresivo)."""
        sigma, Et, _ = self.mat.compute_state(0.0)
        self.assertEqual(sigma, 0.0)
        self.assertEqual(Et, 0.0)

    def test_acceptance_state_vars_passthrough(self):
        """state_vars se devuelve intacto en ambos regímenes."""
        for eps in (-1e-3, 0.0, 1e-3):
            sentinel = {'dummy': 42, 'lista': [1, 2, 3]}
            _, _, returned = self.mat.compute_state(eps, sentinel)
            self.assertIs(returned, sentinel)

    def test_registro_en_registry(self):
        """El material queda registrado vía autodiscover."""
        from fenix.registry import MaterialRegistry
        self.assertIn('CableMaterial1D', MaterialRegistry._items)
        instance = MaterialRegistry.create('CableMaterial1D', E=1000.0, density=0.0)
        self.assertIsInstance(instance, CableMaterial1D)

    def test_validacion_E_positivo(self):
        """E ≤ 0 se rechaza al construir."""
        with self.assertRaises(ValueError):
            CableMaterial1D(E=0.0)
        with self.assertRaises(ValueError):
            CableMaterial1D(E=-1.0)

    def test_switching_exacto_continuidad_de_sigma(self):
        """σ es continuo en ε = 0; E_t tiene un salto E → 0 (H-3.6).

        La unilateralidad del cable implica una discontinuidad **en la
        tangente** al cruzar de tracción a compresión, pero el esfuerzo
        debe permanecer continuo: σ(0⁺) = σ(0) = σ(0⁻) = 0. Verifica que
        la convención de asignar ε = 0 al régimen compresivo no introduce
        rigidez espuria ni esfuerzo residual.
        """
        eps_small = 1.0e-12  # mucho menor que cualquier strain físico
        sigma_pos, Et_pos, _ = self.mat.compute_state(+eps_small)
        sigma_zero, Et_zero, _ = self.mat.compute_state(0.0)
        sigma_neg, Et_neg, _ = self.mat.compute_state(-eps_small)

        # σ continuo en cero (los tres valores ≈ 0)
        self.assertAlmostEqual(sigma_pos, 0.0, delta=self.E * 10 * eps_small)
        self.assertEqual(sigma_zero, 0.0)
        self.assertEqual(sigma_neg, 0.0)

        # E_t con el salto esperado: E en 0⁺, 0 en 0 y 0⁻
        self.assertEqual(Et_pos, self.E)
        self.assertEqual(Et_zero, 0.0)
        self.assertEqual(Et_neg, 0.0)

    def test_switching_exacto_sin_histeresis(self):
        """El material es memoryless: el switching no depende del histórico.

        Carga monotónica positiva → cero → negativa → cero → positiva.
        Cada evaluación es independiente; no debe haber acumulación de
        estado interno que falsifique la respuesta unilateral.
        """
        trayectoria = [1e-4, 5e-5, 0.0, -1e-4, -5e-5, 0.0, 1e-4]
        for eps in trayectoria:
            sigma, Et, vars_out = self.mat.compute_state(eps)
            if eps > 0:
                self.assertAlmostEqual(sigma, self.E * eps,
                                       delta=1e-12 * self.E * abs(eps))
                self.assertEqual(Et, self.E)
            else:
                self.assertEqual(sigma, 0.0)
                self.assertEqual(Et, 0.0)
            self.assertIsNone(vars_out)


if __name__ == '__main__':
    unittest.main()
