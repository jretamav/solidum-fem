"""Fuerzas internas por paso registradas durante el análisis transitorio (deuda #15).

Con ``record_internal_forces=True`` los solvers transitorios evalúan
``elem.internal_forces(u_k)`` en cada paso **con el estado interno de ese
paso** y lo guardan en ``TransientResult.element_forces_history``;
``internal_forces_history(domain)`` lo devuelve tal cual. Sin registro,
la reconstrucción lazy usa el estado final y no es fiel con historia.

Verificación de fidelidad: la fuerza registrada en el paso ``k`` debe
coincidir con la que da un análisis independiente truncado en ``t_k``
(cuyo estado final es, por construcción, el del paso ``k``).
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401
from solidum.core.domain import Domain
from solidum.elements.truss import Truss2D
from solidum.entry import run_transient
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.assembly import Assembler
from solidum.math.solvers import (
    CentralDifferenceSolver, HHTSolver, NewmarkSolver, NewtonHHTSolver, NewtonNewmarkSolver,
)

E, RHO, A, L, OMEGA = 25.0, 2.0, 1.0, 1.0, 5.0
T = 2.0 * math.pi / OMEGA


def _oscillator(material):
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    dom.add_element(Truss2D(1, [n1, n2], material, A=A))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, n2.dofs["ux"]


def _plastic():
    return Elastoplastic1D(E=E, sigma_y=5.0, H=5.0, density=RHO)


def _N(history, k):
    return float(history[1][k].components["N"][0])


class TestRegistroFiel(unittest.TestCase):

    def test_newton_newmark_plastico_coincide_con_analisis_truncado(self):
        dt = T / 200.0
        n_steps = 400
        dom, dof = _oscillator(_plastic())
        v0 = np.zeros(dom.total_dofs); v0[dof] = 6.0
        res = NewtonNewmarkSolver(
            Assembler(dom), t_end=n_steps * dt, dt=dt, u0_dot=v0, lumping="lumped",
            record_internal_forces=True,
        ).solve()
        self.assertIsNotNone(res.element_forces_history)
        hist = res.internal_forces_history(dom)
        self.assertEqual(len(hist[1]), n_steps + 1)
        self.assertGreater(dom.elements[1].state.vars[0]["alpha"], 0.0)

        for k in (n_steps // 3, 2 * n_steps // 3, n_steps):
            dom_k, _ = _oscillator(_plastic())
            res_k = NewtonNewmarkSolver(
                Assembler(dom_k), t_end=k * dt, dt=dt, u0_dot=v0, lumping="lumped",
            ).solve()
            np.testing.assert_allclose(res_k.u_history[dof, -1], res.u_history[dof, k], rtol=1e-12)
            N_k = float(dom_k.elements[1].internal_forces(res_k.u_history[:, -1]).components["N"][0])
            self.assertAlmostEqual(_N(hist, k), N_k, delta=1e-10 * max(1.0, abs(N_k)))

        # La reconstrucción lazy (estado final) NO reproduce los pasos
        # intermedios en régimen plástico: es la deuda que cierra el registro.
        lazy = NewtonNewmarkSolver(
            Assembler(_oscillator(_plastic())[0]), t_end=n_steps * dt, dt=dt, u0_dot=v0,
            lumping="lumped",
        ).solve()
        dom_lazy, _ = _oscillator(_plastic())
        res_lazy = NewtonNewmarkSolver(
            Assembler(dom_lazy), t_end=n_steps * dt, dt=dt, u0_dot=v0, lumping="lumped",
        ).solve()
        self.assertIsNone(res_lazy.element_forces_history)
        with self.assertLogs("solidum.results", level="WARNING"):
            hist_lazy = res_lazy.internal_forces_history(dom_lazy)
        diff = max(abs(_N(hist_lazy, k) - _N(hist, k)) for k in range(n_steps + 1))
        self.assertGreater(diff, 1e-3)
        del lazy

    def test_sin_historia_registro_y_lazy_coinciden_en_todos_los_solvers(self):
        dt = T / 100.0
        for name, make in (
            ("Newmark", lambda asm: NewmarkSolver(asm, t_end=T, dt=dt, record_internal_forces=True)),
            ("HHT", lambda asm: HHTSolver(asm, t_end=T, dt=dt, record_internal_forces=True)),
            ("NewtonNewmark", lambda asm: NewtonNewmarkSolver(asm, t_end=T, dt=dt, record_internal_forces=True)),
            ("NewtonHHT", lambda asm: NewtonHHTSolver(asm, t_end=T, dt=dt, record_internal_forces=True)),
            ("CentralDifference", lambda asm: CentralDifferenceSolver(
                asm, t_end=T, dt=dt, record_internal_forces=True)),
            ("CentralDifference-nl", lambda asm: CentralDifferenceSolver(
                asm, t_end=T, dt=dt, nonlinear=True, record_internal_forces=True)),
        ):
            with self.subTest(solver=name):
                dom, dof = _oscillator(Elastic1D(E=E, density=RHO))
                u0 = np.zeros(dom.total_dofs); u0[dof] = 0.1
                solver = make(Assembler(dom))
                solver.u0 = u0
                res = solver.solve()
                rec = res.internal_forces_history(dom)
                lazy = type(res)(**{**res.__dict__, "element_forces_history": None}
                                 ).internal_forces_history(dom)
                N_rec = np.array([_N(rec, k) for k in range(res.n_steps + 1)])
                N_lazy = np.array([_N(lazy, k) for k in range(res.n_steps + 1)])
                np.testing.assert_array_equal(N_rec, N_lazy)
                np.testing.assert_allclose(N_rec, E * A * res.u_history[dof, :] / L, rtol=1e-10)

    def test_run_transient_pasa_el_kwarg(self):
        dom, dof = _oscillator(Elastic1D(E=E, density=RHO))
        u0 = np.zeros(dom.total_dofs); u0[dof] = 0.1
        res = run_transient(dom, t_end=T, dt=T / 50.0, u0=u0, record_internal_forces=True)
        self.assertIsNotNone(res.element_forces_history)
        self.assertEqual(len(res.element_forces_history[1]), res.n_steps + 1)
        res2 = run_transient(dom, t_end=T, dt=T / 50.0, u0=u0)
        self.assertIsNone(res2.element_forces_history)


if __name__ == '__main__':
    unittest.main()
