"""Tests del análisis dinámico transitorio Newmark — fase 3 del ADR 0009.

Cubre:

1. **Calibración Rayleigh** — `rayleigh_from_modes` reproduce ξ objetivo
   en ω₁ y ω₂; `rayleigh_xi` evalúa correctamente la fórmula.
2. **Oscilador 1 GDL no amortiguado** contra solución analítica
   `u(t) = u₀·cos(ωt) + (u̇₀/ω)·sin(ωt)`.
3. **Oscilador 1 GDL con Rayleigh** contra respuesta sub-amortiguada
   `u(t) = e^(-ξωt) · (cos(ω_d·t) + (ξω/ω_d)·sin(ω_d·t))`.
4. **Step response 1 GDL** con `F(t) = F₀·H(t)`: `u(t) = (F₀/K)·(1 − cos(ωt))`.
5. **Orden de convergencia** O(Δt²) para β=1/4, γ=1/2.
6. **Pipeline YAML** end-to-end con `solver: type: NewmarkSolver`.
"""
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import fenix  # autodiscover
from fenix.core.domain import Domain
from fenix.elements.truss import Truss2D
from fenix.entry import run_transient, run_yaml
from fenix.materials.elastic import Elastic1D
from fenix.math.assembly import Assembler
from fenix.math.damping import rayleigh_from_modes, rayleigh_xi
from fenix.math.solvers import NewmarkSolver
from fenix.results import TransientResult


# ---------------------------------------------------------------------------
# Helper: oscilador 1 GDL como un Truss2D axial empotrado-libre
# ---------------------------------------------------------------------------
# Parámetros escogidos para que el problema reducido tenga K=25, M=1:
#   K_red = EA/L = 25  →  A=1, E=25, L=1
#   M_red = ρ·A·L · 2/6 = ρ/3 = 1  →  ρ=3
# Frecuencia: ω = √(K/M) = 5 rad/s.

E_1DOF = 25.0
RHO_1DOF = 3.0
A_1DOF = 1.0
L_1DOF = 1.0
OMEGA_1DOF = 5.0  # rad/s — verificación analítica abajo


def _build_1dof_oscillator() -> Domain:
    """Truss2D 1 elemento con K_red=25, M_red=1 sobre el único DOF libre ux₂."""
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L_1DOF, 0.0])
    mat = Elastic1D(E=E_1DOF, density=RHO_1DOF)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=A_1DOF))
    # Empotramiento en n1; en n2 solo ux libre.
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers()
    return dom


def _free_dof(dom: Domain) -> int:
    return dom.nodes[2].dofs["ux"]


class TestRayleighCalibration(unittest.TestCase):
    """`rayleigh_from_modes` y `rayleigh_xi` — fórmulas analíticas exactas."""

    def test_reproduces_target_xi_at_calibration_frequencies(self):
        alpha, beta = rayleigh_from_modes(
            xi1=0.02, omega1=10.0, xi2=0.05, omega2=100.0
        )
        self.assertAlmostEqual(rayleigh_xi(alpha, beta, 10.0), 0.02, delta=1e-12)
        self.assertAlmostEqual(rayleigh_xi(alpha, beta, 100.0), 0.05, delta=1e-12)

    def test_equal_omegas_raises(self):
        with self.assertRaisesRegex(ValueError, "indeterminado"):
            rayleigh_from_modes(0.02, 10.0, 0.05, 10.0)

    def test_nonpositive_omega_raises(self):
        with self.assertRaisesRegex(ValueError, "positivas"):
            rayleigh_from_modes(0.02, 0.0, 0.05, 100.0)


class TestNewmark1DofUndamped(unittest.TestCase):
    """Vibración libre sin amortiguamiento — `u(t) = cos(ωt)` exacto."""

    def test_pure_cosine(self):
        dom = _build_1dof_oscillator()
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 100.0
        t_end = 4.0 * T

        ndof = dom.total_dofs
        dof = _free_dof(dom)
        u0 = np.zeros(ndof); u0[dof] = 1.0
        u0_dot = np.zeros(ndof)

        result = run_transient(dom, t_end=t_end, dt=dt, u0=u0, u0_dot=u0_dot)

        self.assertIsInstance(result, TransientResult)
        self.assertTrue(result.converged)
        expected = np.cos(OMEGA_1DOF * result.t_history)
        np.testing.assert_allclose(
            result.u_history[dof, :], expected, atol=0.01,
            err_msg="1 GDL no amortiguado fuera de tolerancia 1%.",
        )

    def test_initial_state_in_history(self):
        dom = _build_1dof_oscillator()
        T = 2.0 * math.pi / OMEGA_1DOF
        ndof = dom.total_dofs
        dof = _free_dof(dom)
        u0 = np.zeros(ndof); u0[dof] = 0.7
        u0_dot = np.zeros(ndof); u0_dot[dof] = 2.0
        result = run_transient(dom, t_end=T, dt=T / 50.0, u0=u0, u0_dot=u0_dot)
        np.testing.assert_allclose(result.u_history[:, 0], u0, atol=1e-12)
        np.testing.assert_allclose(result.udot_history[:, 0], u0_dot, atol=1e-12)


class TestNewmark1DofDamped(unittest.TestCase):
    """Vibración sub-amortiguada con Rayleigh α=0.5, β=0 (⇒ C=0.5, ξ=0.05).

    Solución analítica con u₀=1, u̇₀=0:
      u(t) = e^(-ξωt) · (cos(ω_d·t) + (ξω/ω_d)·sin(ω_d·t))
      ω_d = ω · √(1 − ξ²)
    """

    def test_decaying_oscillation(self):
        dom = _build_1dof_oscillator()
        xi = 0.05
        omega = OMEGA_1DOF
        omega_d = omega * math.sqrt(1.0 - xi * xi)
        T = 2.0 * math.pi / omega
        dt = T / 100.0
        t_end = 4.0 * T

        ndof = dom.total_dofs
        dof = _free_dof(dom)
        u0 = np.zeros(ndof); u0[dof] = 1.0
        u0_dot = np.zeros(ndof)

        # C = α·M + β·K. Con M=1, K=25: α=0.5, β=0 ⇒ C=0.5 ⇒ ξ = C/(2√(KM)) = 0.05.
        result = run_transient(
            dom, t_end=t_end, dt=dt, u0=u0, u0_dot=u0_dot,
            rayleigh={"alpha": 0.5, "beta": 0.0},
        )

        t = result.t_history
        expected = np.exp(-xi * omega * t) * (
            np.cos(omega_d * t) + (xi * omega / omega_d) * np.sin(omega_d * t)
        )
        np.testing.assert_allclose(
            result.u_history[dof, :], expected, atol=0.02,
            err_msg="1 GDL amortiguado fuera de tolerancia 2%.",
        )

    def test_rayleigh_coefficients_propagated_to_result(self):
        dom = _build_1dof_oscillator()
        result = run_transient(
            dom, t_end=1.0, dt=0.01,
            rayleigh={"alpha": 0.3, "beta": 0.01},
        )
        self.assertEqual(result.alpha_rayleigh, 0.3)
        self.assertEqual(result.beta_rayleigh, 0.01)


class TestNewmark1DofStepResponse(unittest.TestCase):
    """Carga súbita `F(t) = F₀·H(t)` sobre 1 GDL: `u(t) = (F₀/K)(1 − cos(ωt))`."""

    def test_step_response(self):
        dom = _build_1dof_oscillator()
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 100.0
        t_end = 2.0 * T

        ndof = dom.total_dofs
        dof = _free_dof(dom)
        F0 = 1.0
        F_vec = np.zeros(ndof); F_vec[dof] = F0
        def F_func(t):  # noqa: ARG001 — F constante (step) tras t=0
            return F_vec

        result = run_transient(dom, t_end=t_end, dt=dt, F_func=F_func)

        t = result.t_history
        K = E_1DOF * A_1DOF / L_1DOF  # 25.0
        expected = (F0 / K) * (1.0 - np.cos(OMEGA_1DOF * t))
        np.testing.assert_allclose(
            result.u_history[dof, :], expected, atol=0.001,
            err_msg="Step response fuera de tolerancia 0.1%.",
        )


class TestNewmarkConvergenceOrder(unittest.TestCase):
    """Esquema β=1/4, γ=1/2 tiene error O(Δt²). El cociente de errores
    al refinar Δt en factor 2 debe tender a 4."""

    def test_second_order_convergence(self):
        T = 2.0 * math.pi / OMEGA_1DOF
        t_end = 2.0 * T

        errors = []
        for n_per_T in (40, 80, 160):
            dom = _build_1dof_oscillator()
            ndof = dom.total_dofs
            dof = _free_dof(dom)
            u0 = np.zeros(ndof); u0[dof] = 1.0
            u0_dot = np.zeros(ndof)
            dt = T / n_per_T
            result = run_transient(dom, t_end=t_end, dt=dt,
                                     u0=u0, u0_dot=u0_dot)
            expected = np.cos(OMEGA_1DOF * result.t_history)
            err = float(np.max(np.abs(result.u_history[dof, :] - expected)))
            errors.append(err)

        # Cocientes consecutivos deben tender a 4 (orden 2 exacto).
        # Tolerancia amplia (3.5–4.5) — el régimen asintótico se alcanza con dt pequeño.
        r1 = errors[0] / errors[1]
        r2 = errors[1] / errors[2]
        self.assertGreater(r1, 3.5, f"Ratio dt→dt/2 = {r1:.3f}, esperado ≈ 4.")
        self.assertLess(r1, 4.5, f"Ratio dt→dt/2 = {r1:.3f}, esperado ≈ 4.")
        self.assertGreater(r2, 3.5, f"Ratio dt/2→dt/4 = {r2:.3f}, esperado ≈ 4.")
        self.assertLess(r2, 4.5, f"Ratio dt/2→dt/4 = {r2:.3f}, esperado ≈ 4.")


class TestNewmarkContract(unittest.TestCase):
    """Validaciones de contrato y errores agregados."""

    def test_runtransient_requires_solver_or_dt(self):
        dom = _build_1dof_oscillator()
        with self.assertRaisesRegex(ValueError, "solver.*t_end.*dt"):
            run_transient(dom)

    def test_negative_dt_raises(self):
        dom = _build_1dof_oscillator()
        asm = Assembler(dom)
        with self.assertRaisesRegex(ValueError, "dt"):
            NewmarkSolver(asm, t_end=1.0, dt=-0.01)

    def test_negative_t_end_raises(self):
        dom = _build_1dof_oscillator()
        asm = Assembler(dom)
        with self.assertRaisesRegex(ValueError, "t_end"):
            NewmarkSolver(asm, t_end=-1.0, dt=0.01)

    def test_invalid_rayleigh_dict_raises(self):
        dom = _build_1dof_oscillator()
        with self.assertRaisesRegex(ValueError, "rayleigh"):
            run_transient(dom, t_end=1.0, dt=0.01,
                           rayleigh={"foo": 1.0, "bar": 2.0})


class TestNewmarkYamlPipeline(unittest.TestCase):
    """`solver: type: NewmarkSolver` desde YAML."""

    def _write_yaml(self, content: str) -> Path:
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False,
                                          encoding="utf-8")
        f.write(content); f.close()
        return Path(f.name)

    def test_yaml_undamped_1dof(self):
        """Pipeline YAML completo del oscilador 1 GDL no amortiguado."""
        yaml_content = f"""
nodes:
  - {{id: 1, coords: [0.0, 0.0]}}
  - {{id: 2, coords: [{L_1DOF}, 0.0]}}

materials:
  - {{id: 1, type: Elastic1D, E: {E_1DOF}, density: {RHO_1DOF}}}

elements:
  - {{id: 1, type: Truss2D, material: 1, nodes: [1, 2], A: {A_1DOF}}}

boundary_conditions_by_node:
  - {{node_id: 1, ux: 0.0, uy: 0.0}}
  - {{node_id: 2, uy: 0.0}}

solver:
  type: NewmarkSolver
  t_end: {4.0 * 2.0 * math.pi / OMEGA_1DOF}
  dt: {2.0 * math.pi / OMEGA_1DOF / 50.0}
"""
        path = self._write_yaml(yaml_content)
        try:
            result = run_yaml(str(path))
        finally:
            path.unlink()

        self.assertIsInstance(result, TransientResult)
        # Vibración libre desde reposo (u0=u0_dot=0, F=0) — todo cero.
        np.testing.assert_allclose(result.u_history, 0.0, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
