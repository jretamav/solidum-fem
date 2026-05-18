"""Tests del integrador HHT-α (Hilber-Hughes-Taylor 1977).

Variante de Newmark con disipación numérica controlada por ``α ∈ [−1/3, 0]``.
Cubre los acceptance de ``docs/specs/HHTSolver.md``:

1. **Recovery** — ``HHTSolver(alpha=0.0)`` reproduce ``NewmarkSolver`` con
   ``β=1/4, γ=1/2`` bit-a-bit en oscilador 1 GDL.
2. **Disipación de altas frecuencias** — sistema con dos modos bien separados;
   ``alpha=-0.1`` atenúa el modo alto sin tocar el bajo.
3. **Estabilidad incondicional** — paso grande no destruye la solución.
4. **Auto-derivación de β, γ** — con `alpha` solo, β y γ son los canónicos.
5. **Validación del alpha range** — `ValueError` fuera de [−1/3, 0].
6. **No lineal** — `NewtonHHTSolver` con material elastoplástico, recovery
   lineal (α=0) y disipación funcional.
"""
import math
import unittest

import numpy as np

import fenix  # autodiscover
from fenix.core.domain import Domain
from fenix.elements.truss import Truss2D
from fenix.materials.elastic import Elastic1D
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.solvers import (
    HHTSolver,
    NewmarkSolver,
    NewtonHHTSolver,
    NewtonNewmarkSolver,
)
from fenix.math.solvers.newmark import _hht_autoderive_beta_gamma


# Mismos parámetros que test_newmark.py para 1 GDL:
#   K_red = EA/L = 25, M_red = ρAL·2/6 = 1, ω = 5 rad/s.
E_1DOF = 25.0
RHO_1DOF = 3.0
A_1DOF = 1.0
L_1DOF = 1.0
OMEGA_1DOF = 5.0


def _build_1dof_oscillator(material=None) -> Domain:
    """Truss2D 1 elemento con K_red=25, M_red=1."""
    if material is None:
        material = Elastic1D(E=E_1DOF, density=RHO_1DOF)
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L_1DOF, 0.0])
    dom.add_element(Truss2D(1, [n1, n2], material, A=A_1DOF))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers()
    return dom


def _free_dof(dom: Domain) -> int:
    return dom.nodes[2].dofs["ux"]


# =============================================================================
# Auto-derivación de β, γ
# =============================================================================

class TestHHTAutoDerivation(unittest.TestCase):
    """Verifica las fórmulas canónicas β=(1−α)²/4, γ=(1−2α)/2 (Hilber 1977)."""

    def test_alpha_zero_recovers_newmark_trapezoidal(self):
        beta, gamma = _hht_autoderive_beta_gamma(0.0)
        self.assertAlmostEqual(beta, 0.25, places=12)
        self.assertAlmostEqual(gamma, 0.5, places=12)

    def test_alpha_minus_one_third_yields_canonical(self):
        beta, gamma = _hht_autoderive_beta_gamma(-1.0 / 3.0)
        # β = (1 + 1/3)² / 4 = (4/3)² / 4 = 16/9 / 4 = 4/9 ≈ 0.4444
        # γ = (1 + 2/3) / 2 = 5/6 ≈ 0.8333
        self.assertAlmostEqual(beta, 4.0 / 9.0, places=12)
        self.assertAlmostEqual(gamma, 5.0 / 6.0, places=12)


# =============================================================================
# Validación del rango de α
# =============================================================================

class TestHHTAlphaRangeValidation(unittest.TestCase):
    """α fuera de [−1/3, 0] pierde estabilidad incondicional u orden 2."""

    def test_alpha_positive_rejected(self):
        dom = _build_1dof_oscillator()
        with self.assertRaises(ValueError):
            HHTSolver(Assembler(dom), t_end=1.0, dt=0.01, alpha=0.1)

    def test_alpha_below_minus_one_third_rejected(self):
        dom = _build_1dof_oscillator()
        with self.assertRaises(ValueError):
            HHTSolver(Assembler(dom), t_end=1.0, dt=0.01, alpha=-0.5)

    def test_alpha_boundaries_accepted(self):
        dom = _build_1dof_oscillator()
        # Ambos límites del rango deben aceptarse sin error.
        HHTSolver(Assembler(dom), t_end=1.0, dt=0.01, alpha=0.0)
        HHTSolver(Assembler(dom), t_end=1.0, dt=0.01, alpha=-1.0 / 3.0)

    def test_newton_hht_validates_alpha(self):
        dom = _build_1dof_oscillator()
        with self.assertRaises(ValueError):
            NewtonHHTSolver(Assembler(dom), t_end=1.0, dt=0.01, alpha=0.5)


# =============================================================================
# Recovery: α=0 coincide con Newmark trapezoidal
# =============================================================================

class TestHHTRecoversNewmark(unittest.TestCase):
    """Con α=0, HHT colapsa exactamente a Newmark trapezoidal (β=1/4, γ=1/2)."""

    def test_undamped_free_vibration_matches_newmark(self):
        """Oscilador 1 GDL no amortiguado en vibración libre."""
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 50.0
        t_end = 2.0 * T

        # Newmark trapezoidal de referencia.
        dom_ref = _build_1dof_oscillator()
        ndof = dom_ref.total_dofs
        u0_ref = np.zeros(ndof); u0_ref[_free_dof(dom_ref)] = 0.1
        res_ref = NewmarkSolver(
            Assembler(dom_ref), t_end=t_end, dt=dt,
            beta=0.25, gamma=0.5, u0=u0_ref,
        ).solve()

        # HHTSolver con α=0.
        dom_hht = _build_1dof_oscillator()
        u0_hht = np.zeros(dom_hht.total_dofs); u0_hht[_free_dof(dom_hht)] = 0.1
        res_hht = HHTSolver(
            Assembler(dom_hht), t_end=t_end, dt=dt, alpha=0.0, u0=u0_hht,
        ).solve()

        np.testing.assert_allclose(res_hht.u_history, res_ref.u_history, atol=1.0e-12)
        np.testing.assert_allclose(res_hht.udot_history, res_ref.udot_history, atol=1.0e-11)

    def test_step_response_matches_newmark(self):
        """Step response: F constante a partir de t=0."""
        F0 = 1.0
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 60.0
        t_end = 1.5 * T

        dof = _free_dof(_build_1dof_oscillator())

        def F_func(t):
            v = np.zeros(4)
            v[dof] = F0
            return v

        # Newmark trapezoidal.
        dom_ref = _build_1dof_oscillator()
        res_ref = NewmarkSolver(
            Assembler(dom_ref), t_end=t_end, dt=dt,
            beta=0.25, gamma=0.5, F_func=F_func,
        ).solve()

        # HHT con α=0.
        dom_hht = _build_1dof_oscillator()
        res_hht = HHTSolver(
            Assembler(dom_hht), t_end=t_end, dt=dt, alpha=0.0, F_func=F_func,
        ).solve()

        np.testing.assert_allclose(res_hht.u_history, res_ref.u_history, atol=1.0e-12)

    def test_with_rayleigh_matches_newmark(self):
        """Con amortiguamiento Rayleigh, α=0 sigue coincidiendo con Newmark."""
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 80.0
        t_end = 3.0 * T

        rayleigh_cfg = {"alpha": 0.5, "beta": 0.001}

        dom_ref = _build_1dof_oscillator()
        u0 = np.zeros(dom_ref.total_dofs); u0[_free_dof(dom_ref)] = 0.1
        res_ref = NewmarkSolver(
            Assembler(dom_ref), t_end=t_end, dt=dt,
            beta=0.25, gamma=0.5, u0=u0, rayleigh=rayleigh_cfg,
        ).solve()

        dom_hht = _build_1dof_oscillator()
        u0_hht = np.zeros(dom_hht.total_dofs); u0_hht[_free_dof(dom_hht)] = 0.1
        res_hht = HHTSolver(
            Assembler(dom_hht), t_end=t_end, dt=dt, alpha=0.0,
            u0=u0_hht, rayleigh=rayleigh_cfg,
        ).solve()

        np.testing.assert_allclose(res_hht.u_history, res_ref.u_history, atol=1.0e-11)


# =============================================================================
# Disipación de altas frecuencias
# =============================================================================

class TestHHTHighFrequencyDissipation(unittest.TestCase):
    """α<0 disipa modos con ωΔt grande sin afectar modos resueltos.

    Estrategia: oscilador 1 GDL con Δt sobredimensionado de forma que ωΔt > 1
    (modo mal resuelto, "alta frecuencia" desde la óptica del integrador).
    Con α=0 (Newmark trapezoidal) el modo persiste con amplitud completa;
    con α=−0.1 el modo se atenúa progresivamente.
    """

    def test_high_frequency_mode_attenuated_by_negative_alpha(self):
        """Modo con ωΔt ~ 1.6 (mal resuelto): α<0 lo disipa, α=0 no."""
        # Con Δt = T/4 → ωΔt = 2π/4 ≈ 1.57, modo mal resuelto.
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 4.0
        n_periods = 20
        t_end = n_periods * T

        # Caso α=0 (Newmark trapezoidal sin disipación).
        dom_0 = _build_1dof_oscillator()
        u0_0 = np.zeros(dom_0.total_dofs); u0_0[_free_dof(dom_0)] = 1.0
        res_0 = HHTSolver(
            Assembler(dom_0), t_end=t_end, dt=dt, alpha=0.0, u0=u0_0,
        ).solve()

        # Caso α=−0.1.
        dom_1 = _build_1dof_oscillator()
        u0_1 = np.zeros(dom_1.total_dofs); u0_1[_free_dof(dom_1)] = 1.0
        res_1 = HHTSolver(
            Assembler(dom_1), t_end=t_end, dt=dt, alpha=-0.1, u0=u0_1,
        ).solve()

        dof = _free_dof(dom_0)
        amp_0_final = float(np.abs(res_0.u_history[dof, :]).max())
        amp_1_final = float(np.abs(res_1.u_history[dof, -10:]).max())

        # α=0: amplitud no decae (trapezoidal preserva energía).
        # α=−0.1: amplitud decae significativamente tras 20 "periodos" mal resueltos.
        self.assertGreater(amp_0_final, 0.9,
            f"α=0 debería preservar amplitud (≥0.9), obtuvo {amp_0_final:.4f}")
        self.assertLess(amp_1_final, 0.5 * amp_0_final,
            f"α=−0.1 debería atenuar al menos 50%; "
            f"amp_α=0={amp_0_final:.4f}, amp_α=−0.1={amp_1_final:.4f}")

    def test_well_resolved_mode_essentially_preserved_by_alpha(self):
        """Modo bien resuelto (Δt = T/100): α=−0.05 NO degrada apenas la amplitud.

        Verifica que la disipación HHT es selectiva: no toca modos resueltos.
        """
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 100.0  # bien resuelto
        n_periods = 5
        t_end = n_periods * T

        dom = _build_1dof_oscillator()
        u0 = np.zeros(dom.total_dofs); u0[_free_dof(dom)] = 1.0
        res = HHTSolver(
            Assembler(dom), t_end=t_end, dt=dt, alpha=-0.05, u0=u0,
        ).solve()

        dof = _free_dof(dom)
        amp_final = float(np.abs(res.u_history[dof, -20:]).max())
        # Tras 5 periodos bien resueltos, amplitud > 90% del valor inicial.
        self.assertGreater(amp_final, 0.9,
            f"Modo bien resuelto degradado en exceso: amp_final={amp_final:.4f}")


# =============================================================================
# Estabilidad incondicional
# =============================================================================

class TestHHTUnconditionalStability(unittest.TestCase):
    """Δt arbitrariamente grande no produce divergencia para α ∈ [−1/3, 0]."""

    def test_large_dt_does_not_diverge(self):
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = 50.0 * T  # 50 periodos por paso (totalmente sobredimensionado)
        t_end = 10.0 * dt

        for alpha in (0.0, -0.05, -0.15, -1.0 / 3.0):
            with self.subTest(alpha=alpha):
                dom = _build_1dof_oscillator()
                u0 = np.zeros(dom.total_dofs); u0[_free_dof(dom)] = 1.0
                res = HHTSolver(
                    Assembler(dom), t_end=t_end, dt=dt, alpha=alpha, u0=u0,
                ).solve()
                # Solución acotada: no debe haber NaN ni Inf, ni amplitud >> 1.
                max_abs = float(np.abs(res.u_history).max())
                self.assertTrue(np.isfinite(max_abs),
                    f"α={alpha}: solución no finita")
                self.assertLess(max_abs, 10.0,
                    f"α={alpha}: amplitud diverge a {max_abs:.4f}")


# =============================================================================
# Variante no lineal — NewtonHHTSolver
# =============================================================================

class TestNewtonHHT(unittest.TestCase):
    """NewtonHHTSolver — recuperación del caso lineal y régimen plástico."""

    def test_linear_material_recovers_HHT_lineal(self):
        """Con material elástico, NewtonHHT = HHT (residuo se anula en 1 iter)."""
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 50.0
        t_end = 2.0 * T

        # HHT lineal de referencia.
        dom_ref = _build_1dof_oscillator(Elastic1D(E=E_1DOF, density=RHO_1DOF))
        u0_ref = np.zeros(dom_ref.total_dofs); u0_ref[_free_dof(dom_ref)] = 0.1
        res_ref = HHTSolver(
            Assembler(dom_ref), t_end=t_end, dt=dt, alpha=-0.05, u0=u0_ref,
        ).solve()

        # NewtonHHT con material elástico (vía Elastoplastic1D con yield enorme).
        mat_quasi_elastic = Elastoplastic1D(
            E=E_1DOF, sigma_y=1.0e6, H=0.0, density=RHO_1DOF,
        )
        dom_nl = _build_1dof_oscillator(mat_quasi_elastic)
        u0_nl = np.zeros(dom_nl.total_dofs); u0_nl[_free_dof(dom_nl)] = 0.1
        res_nl = NewtonHHTSolver(
            Assembler(dom_nl), t_end=t_end, dt=dt, alpha=-0.05, u0=u0_nl,
            convergence=ConvergenceCriterion(rtol_force=1.0e-10, rtol_disp=1.0e-10),
            max_iter=5,
        ).solve()

        np.testing.assert_allclose(res_nl.u_history, res_ref.u_history, atol=1.0e-9)

    def test_alpha_zero_recovers_NewtonNewmark(self):
        """Con α=0, NewtonHHT coincide bit-a-bit con NewtonNewmarkSolver."""
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 60.0
        t_end = 1.5 * T

        # Material plástico activo: yield bajo para que sí plastifique.
        E_loc = E_1DOF
        rho_loc = RHO_1DOF
        F0 = 1.0
        dof_idx = _free_dof(_build_1dof_oscillator())

        def F_func(t):
            v = np.zeros(4)
            v[dof_idx] = F0
            return v

        # NewtonNewmark trapezoidal.
        mat_a = Elastoplastic1D(E=E_loc, sigma_y=0.3, H=0.5, density=rho_loc)
        dom_a = _build_1dof_oscillator(mat_a)
        res_a = NewtonNewmarkSolver(
            Assembler(dom_a), t_end=t_end, dt=dt,
            beta=0.25, gamma=0.5, F_func=F_func,
            convergence=ConvergenceCriterion(rtol_force=1.0e-10, rtol_disp=1.0e-10),
        ).solve()

        # NewtonHHT con α=0.
        mat_b = Elastoplastic1D(E=E_loc, sigma_y=0.3, H=0.5, density=rho_loc)
        dom_b = _build_1dof_oscillator(mat_b)
        res_b = NewtonHHTSolver(
            Assembler(dom_b), t_end=t_end, dt=dt, alpha=0.0, F_func=F_func,
            convergence=ConvergenceCriterion(rtol_force=1.0e-10, rtol_disp=1.0e-10),
        ).solve()

        np.testing.assert_allclose(res_b.u_history, res_a.u_history, atol=1.0e-10)
        # Plasticidad activada en ambos.
        elem_a = list(dom_a.elements.values())[0]
        elem_b = list(dom_b.elements.values())[0]
        self.assertGreater(elem_a.state.vars[0]['alpha'], 0.0)
        self.assertAlmostEqual(
            elem_a.state.vars[0]['alpha'], elem_b.state.vars[0]['alpha'],
            places=10,
        )

    def test_plastic_oscillator_converges_with_alpha_negative(self):
        """Oscilador 1 GDL elastoplástico con α<0: converge y plastifica."""
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 50.0
        t_end = 2.0 * T

        mat = Elastoplastic1D(E=E_1DOF, sigma_y=0.2, H=0.0, density=RHO_1DOF)
        dom = _build_1dof_oscillator(mat)
        u0 = np.zeros(dom.total_dofs); u0[_free_dof(dom)] = 0.05
        res = NewtonHHTSolver(
            Assembler(dom), t_end=t_end, dt=dt, alpha=-0.1, u0=u0,
            convergence=ConvergenceCriterion(rtol_force=1.0e-7, rtol_disp=1.0e-7),
        ).solve()

        self.assertTrue(res.converged)
        # La plasticidad se ha activado.
        elem = list(dom.elements.values())[0]
        self.assertGreater(elem.state.vars[0]['alpha'], 0.0,
            "El test no ejerció su régimen objetivo: no plastificó.")


if __name__ == '__main__':
    unittest.main()
