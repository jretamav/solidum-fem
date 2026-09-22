"""Regresiones del lote dinámico de la auditoría 2026-09-22.

- Newton-Newmark / Newton-HHT: ü₀ consistente con apoyos prescritos no
  nulos (``g`` se impone antes de ensamblar ``F_int(u₀)``).
- Newton-Newmark / Newton-HHT con restricciones lineales (MPC con
  maestros): la incógnita reducida se reconstruye con ``T`` y el esclavo
  sigue al maestro.
- ``number_of_steps``: ``t_end = n·dt`` no añade un paso espurio.
- ``LUSolver.solve`` singular lanza ``RuntimeError`` en vez de ``NaN``.
- Newmark con DOF libre sin masa lanza ``ValueError`` accionable.
- Modal free-free con ``sigma = 0`` funciona para cualquier malla.
"""
import math
import unittest

import numpy as np
import scipy.sparse as sp

import solidum  # noqa: F401 — autodiscover
from solidum.core.domain import Domain
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.linalg import LUSolver
from solidum.math.solvers import (
    HHTSolver,
    ModalSolver,
    NewmarkSolver,
    NewtonHHTSolver,
    NewtonNewmarkSolver,
)
from solidum.math.solvers._shared import number_of_steps

E, RHO, A, L = 25.0, 2.0, 1.0, 1.0
OMEGA = 5.0  # K_red = 25, M_red (lumped) = 1


def _oscillator(g0=0.0):
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    dom.add_element(Truss2D(1, [n1, n2], Elastic1D(E=E, density=RHO), A=A))
    n1.fix_dof("ux", g0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, n2.dofs["ux"]


def _two_parallel_bars_periodic():
    """Dos barras paralelas independientes con periodicidad ``u4 = u2``."""
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0]); n2 = dom.add_node(2, [L, 0.0])
    n3 = dom.add_node(3, [0.0, 1.0]); n4 = dom.add_node(4, [L, 1.0])
    mat = Elastic1D(E=E, density=RHO)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
    dom.add_element(Truss2D(2, [n3, n4], mat, A=A))
    for n in (n1, n3):
        n.fix_dof("ux", 0.0); n.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0); n4.fix_dof("uy", 0.0)
    dom.add_linear_constraint(slave=(4, "ux"), masters=[(2, "ux")], coefficients=[1.0])
    dom.generate_equation_numbers(verbose=False)
    return dom, n2.dofs["ux"], n4.dofs["ux"]


class TestPrescribedSupportInitialAcceleration(unittest.TestCase):

    def _run(self, cls, dom, dt, t_end, u0):
        kwargs = dict(t_end=t_end, dt=dt, u0=u0, lumping="lumped")
        return cls(Assembler(dom), **kwargs).solve()

    def test_newton_newmark_matches_linear_with_nonzero_support(self):
        T = 2.0 * math.pi / OMEGA
        dt, t_end, g0 = T / 200.0, T, 0.5
        dom_l, dof = _oscillator(g0)
        u0 = np.zeros(dom_l.total_dofs); u0[dof] = 1.0
        # u0 del usuario NO impone g en el apoyo: el solver debe proyectarlo.
        res_l = self._run(NewmarkSolver, dom_l, dt, t_end, u0)
        dom_n, _ = _oscillator(g0)
        res_n = self._run(NewtonNewmarkSolver, dom_n, dt, t_end, u0)
        np.testing.assert_allclose(res_n.u_history, res_l.u_history, rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(res_n.uddot_history[:, 0], res_l.uddot_history[:, 0],
                                   rtol=1e-10, atol=1e-12)
        # ü₀ analítico: K(g₀ − u₀)/M = 25·(0.5 − 1)/1 = −12.5
        self.assertAlmostEqual(res_n.uddot_history[dof, 0], -12.5, places=9)
        # La trayectoria analítica u = g₀ + (u₀ − g₀)cos ωt.
        u_an = g0 + (1.0 - g0) * np.cos(OMEGA * res_n.t_history)
        self.assertLess(np.max(np.abs(res_n.u_history[dof] - u_an)), 5e-3)

    def test_newton_hht_matches_linear_with_nonzero_support(self):
        T = 2.0 * math.pi / OMEGA
        dt, t_end, g0 = T / 200.0, T, 0.5
        dom_l, dof = _oscillator(g0)
        u0 = np.zeros(dom_l.total_dofs); u0[dof] = 1.0
        res_l = HHTSolver(Assembler(dom_l), t_end=t_end, dt=dt, u0=u0, alpha=-0.05,
                          lumping="lumped").solve()
        dom_n, _ = _oscillator(g0)
        res_n = NewtonHHTSolver(Assembler(dom_n), t_end=t_end, dt=dt, u0=u0, alpha=-0.05,
                                lumping="lumped").solve()
        np.testing.assert_allclose(res_n.u_history, res_l.u_history, rtol=1e-9, atol=1e-12)


class TestLinearConstraintsInNewtonDynamics(unittest.TestCase):

    def test_newton_newmark_and_hht_honour_periodicity(self):
        T = 2.0 * math.pi / OMEGA
        dt, t_end = T / 200.0, T
        dom_l, d2, d4 = _two_parallel_bars_periodic()
        u0 = np.zeros(dom_l.total_dofs); u0[d2] = 1.0; u0[d4] = 1.0
        ref = NewmarkSolver(Assembler(dom_l), t_end=t_end, dt=dt, u0=u0, lumping="lumped").solve()
        for cls, extra in ((NewtonNewmarkSolver, {}), (NewtonHHTSolver, {"alpha": -0.05})):
            dom, d2, d4 = _two_parallel_bars_periodic()
            res = cls(Assembler(dom), t_end=t_end, dt=dt, u0=u0, lumping="lumped", **extra).solve()
            self.assertLess(float(np.max(np.abs(res.u_history[d4] - res.u_history[d2]))), 1e-12,
                            f"{cls.__name__}: el esclavo no sigue al maestro")
            if cls is NewtonNewmarkSolver:
                np.testing.assert_allclose(res.u_history, ref.u_history, rtol=1e-9, atol=1e-12)
            # Analítico: dos masas acopladas → mismo oscilador ω = 5.
            u_an = np.cos(OMEGA * res.t_history)
            self.assertLess(float(np.max(np.abs(res.u_history[d2] - u_an))), 5e-3)


class TestStepCountAndMassDiagnostics(unittest.TestCase):

    def test_number_of_steps_no_spurious_extra_step(self):
        self.assertEqual(number_of_steps(3 * 0.1, 0.1), 3)
        self.assertEqual(number_of_steps(1.0, 0.3), 4)
        self.assertEqual(number_of_steps(7 * 0.7, 0.7), 7)
        for n in range(1, 60):
            for dt in (0.1, 0.7, 0.01, 1e-3):
                self.assertEqual(number_of_steps(n * dt, dt), n, (n, dt))

    def test_newmark_final_time_is_t_end(self):
        dom, _ = _oscillator()
        res = NewmarkSolver(Assembler(dom), t_end=3 * 0.1, dt=0.1, lumping="lumped").solve()
        self.assertEqual(res.n_steps, 3)
        self.assertAlmostEqual(res.t_history[-1], 0.3, places=12)

    def test_lu_solver_singular_raises(self):
        K = sp.csr_matrix(np.array([[1.0, 0.0], [0.0, 0.0]]))
        with self.assertRaises(RuntimeError):
            LUSolver().solve(K, np.array([1.0, 1.0]))

    def test_massless_free_dof_raises_value_error(self):
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0]); n2 = dom.add_node(2, [1.0, 0.0])
        n3 = dom.add_node(3, [2.0, 0.0])
        dom.add_element(Truss2D(1, [n1, n2], Elastic1D(E=E, density=RHO), A=A))
        dom.add_element(Truss2D(2, [n2, n3], Elastic1D(E=E, density=0.0), A=A))
        n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
        n2.fix_dof("uy", 0.0); n3.fix_dof("uy", 0.0)
        dom.generate_equation_numbers(verbose=False)
        with self.assertRaises(ValueError):
            NewmarkSolver(Assembler(dom), t_end=0.1, dt=0.01, lumping="lumped").solve()
        with self.assertRaises(ValueError):
            NewtonNewmarkSolver(Assembler(dom), t_end=0.1, dt=0.01, lumping="lumped").solve()


class TestModalFreeFree(unittest.TestCase):

    def test_sigma_zero_free_free_any_mesh(self):
        for n_el in (1, 2, 3, 4, 10):
            dom = Domain()
            for i in range(n_el + 1):
                dom.add_node(i + 1, [i * L / n_el, 0.0])
            mat = Elastic1D(E=E, density=RHO)
            for e in range(n_el):
                dom.add_element(Truss2D(e + 1, [dom.nodes[e + 1], dom.nodes[e + 2]], mat, A=A))
            for nd in dom.nodes.values():
                nd.fix_dof("uy", 0.0)
            dom.generate_equation_numbers(verbose=False)
            n_modes = min(2, n_el)
            res = ModalSolver(Assembler(dom), n_modes=n_modes, lumping="lumped").solve()
            self.assertTrue(res.converged)
            self.assertAlmostEqual(res.frequencies_rad[0], 0.0, places=6, msg=f"n_el={n_el}")
            if n_modes > 1:
                self.assertGreater(res.frequencies_rad[1], 0.0)


if __name__ == "__main__":
    unittest.main()
