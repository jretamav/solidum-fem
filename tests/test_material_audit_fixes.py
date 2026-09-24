"""Regresiones del lote de materiales de la auditoría 2026-09-22.

- Drucker-Prager con ψ = 0 y H = 0: el retorno al ápice no tiene solución
  y se rechaza con ``ValueError`` (antes: α ≈ 1e30 y σ inconsistente).
- Cohesivo lineal: ``w_c > κ_0`` obligatorio (misma condición que el
  exponencial).
- J2 plane stress: el Newton local converge también con predictores muy
  lejanos (warm start + tolerancia relativa a ``R_curr``).
- Timoshenko con ``E_t = 0`` no divide por cero.
"""
import logging
import unittest

import numpy as np

import solidum  # noqa: F401
from solidum.core.domain import Domain
from solidum.elements.frame import Frame2DTimoshenko
from solidum.materials.drucker_prager_2d import DruckerPrager2D
from solidum.materials.drucker_prager_3d import DruckerPrager3D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.materials.von_mises_2d import VonMises2D


class TestDruckerPragerApex(unittest.TestCase):

    def test_psi_zero_h_zero_apex_raises(self):
        with self.assertLogs("solidum.materials", level=logging.WARNING):
            m2 = DruckerPrager2D(E=200.0, nu=0.3, cohesion=1.0, phi_deg=30.0, psi_deg=0.0, H=0.0)
        with self.assertRaises(ValueError):
            m2.compute_state(np.array([0.03, 0.03, 0.0]))
        m3 = DruckerPrager3D(E=200.0, nu=0.3, cohesion=1.0, phi_deg=30.0, psi_deg=0.0, H=0.0)
        with self.assertRaises(ValueError):
            m3.compute_state(np.array([0.03, 0.03, 0.03, 0.0, 0.0, 0.0]))

    def test_psi_zero_h_zero_still_works_in_compression(self):
        m2 = DruckerPrager2D(E=200.0, nu=0.3, cohesion=1.0, phi_deg=30.0, psi_deg=0.0, H=0.0)
        sig, C, st = m2.compute_state(np.array([-0.05, -0.01, 0.02]))
        self.assertTrue(np.all(np.isfinite(sig)))
        self.assertLess(st["alpha"], 1.0)

    def test_psi_zero_apex_raises_even_with_hardening(self):
        # Con ψ = 0 y H > 0 la versión anterior devolvía el σ del ápice sin la
        # ε^p volumétrica que lo justifica (σ ≠ C_e:(ε − ε^p)).
        m2 = DruckerPrager2D(E=200.0, nu=0.3, cohesion=1.0, phi_deg=30.0, psi_deg=0.0, H=10.0)
        with self.assertRaises(ValueError):
            m2.compute_state(np.array([0.03, 0.03, 0.0]))

    def test_apex_with_dilatancy_stays_consistent(self):
        # Con ψ > 0 el ápice sí admite solución, también con H = 0: σ deriva
        # del estado interno vía σ = C_e:(ε − ε^p).
        m2 = DruckerPrager2D(E=200.0, nu=0.3, cohesion=1.0, phi_deg=30.0, psi_deg=15.0, H=0.0)
        eps = np.array([0.03, 0.03, 0.0])
        sig, _, st = m2.compute_state(eps)
        ep = st["eps_p"]
        # Relacion elastica 3D con a = eps - eps_p (a_zz = -eps_p_zz en plane
        # strain): sigma_ij = lambda tr(a) delta_ij + 2G a_ij. La matriz C_e
        # 3x3 de plane strain no sirve porque ignora eps_p_zz.
        lam = m2.K - 2.0 * m2.G / 3.0
        a = np.array([eps[0] - ep[0], eps[1] - ep[1], -ep[2], 0.5 * eps[2] - ep[3]])
        tr = a[0] + a[1] + a[2]
        sig_chk = np.array([lam * tr + 2 * m2.G * a[0], lam * tr + 2 * m2.G * a[1], 2 * m2.G * a[3]])
        np.testing.assert_allclose(sig, sig_chk, rtol=1e-8, atol=1e-10)
        self.assertLess(st["alpha"], 1.0)


class TestJ2PlaneStressLocalNewton(unittest.TestCase):

    def _vm(self, H):
        return VonMises2D(E=200e9, nu=0.3, sigma_y=250e6, H=H, hypothesis="plane_stress")

    def _on_surface(self, mat, sig, alpha):
        P_sig = np.array([(2 * sig[0] - sig[1]) / 3, (-sig[0] + 2 * sig[1]) / 3, 2 * sig[2]])
        f_bar = 0.5 * sig @ P_sig - (mat.sigma_y + mat.H * alpha) ** 2 / 3
        return abs(f_bar) / ((mat.sigma_y + mat.H * alpha) ** 2 / 3)

    def test_far_predictors_converge_h_zero_and_h_positive(self):
        for H in (0.0, 20e9):
            mat = self._vm(H)
            for scale in (1e-3, 1e-2, 1.0, 5.0):   # ε_trial hasta 5 (σ_trial/R ~ 4·10³)
                eps = np.array([scale, -0.2 * scale, 0.7 * scale])
                sig, C, st = mat.compute_state(eps)
                if st["alpha"] > 0.0:   # rama plástica: σ sobre la superficie
                    self.assertLess(self._on_surface(mat, sig, st["alpha"]), 1e-9, (H, scale))
                elif scale >= 1e-2:
                    self.fail(f"se esperaba plastificar con ε = {scale}")
                self.assertFalse(mat._local_newton_warned)
                # Tangente simétrica y finita.
                self.assertTrue(np.all(np.isfinite(C)))
                np.testing.assert_allclose(C, C.T, rtol=1e-10, atol=1e-6 * np.abs(C).max())

    def test_exhaustion_is_logged_not_silent(self):
        mat = self._vm(0.0)
        # Presupuesto artificialmente bajo para forzar el agotamiento.
        import solidum.materials.von_mises_2d as vm
        original = vm._PLANE_STRESS_MAX_LOCAL_ITER
        vm._PLANE_STRESS_MAX_LOCAL_ITER = 1
        try:
            with self.assertLogs("solidum.materials", level=logging.WARNING):
                mat.compute_state(np.array([0.5, 0.0, 0.0]))
        finally:
            vm._PLANE_STRESS_MAX_LOCAL_ITER = original


class TestTimoshenkoPerfectPlasticity(unittest.TestCase):

    def test_zero_tangent_modulus_does_not_divide_by_zero(self):
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0]); n2 = dom.add_node(2, [2.0, 0.0])
        mat = Elastoplastic1D(E=200e9, sigma_y=250e6, H=0.0)
        el = Frame2DTimoshenko(1, [n1, n2], mat, A=1e-2, I=1e-4, As=8e-3, nu=0.3)
        dom.add_element(el); dom.generate_equation_numbers()
        u = np.zeros(6)
        u[3] = 0.01     # ε = 0.005 ≫ ε_y = 1.25e-3 → E_t = 0
        K, F = el.compute_element_state(u)
        self.assertTrue(np.all(np.isfinite(K)))
        self.assertAlmostEqual(F[3], 250e6 * 1e-2, places=3)


if __name__ == "__main__":
    unittest.main()
