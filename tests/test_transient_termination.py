"""Terminación de los análisis en el tiempo (ADR 0025, ampliado a los
transitorios).

Mismo contrato que los solvers estáticos incrementales:

- Al llegar al instante final: ``termination = "t_end"``, ``completed``,
  ``steps_done``, ``t_final`` y ``U_final`` coherentes con el resultado, que
  lleva ``termination = "t_end"`` y ``converged = True``.
- Un análisis que se detiene (Newton sin convergencia en Newmark/HHT no
  lineales, integración inestable en la diferencia central o en el θ-método)
  registra una línea ``ERROR`` y lanza una ``SolverDivergedError`` tipada con
  el historial hasta el último paso convergido en ``partial_result`` (también
  en el solver), sin columnas sin calcular. El estado del solver queda en
  ``"diverged"``.
- Un grado de libertad libre sin masa en un transitorio da
  ``MasslessDofsError`` que lo nombra (``"nodo 2 (ux)"``), con la causa
  ``initial_acceleration`` en los implícitos y ``explicit_integration`` en la
  diferencia central, en vez de un ``ValueError`` genérico.
- Los entrypoints de los análisis no estáticos son API pública del paquete.
"""
import logging
import unittest

import numpy as np

import solidum
from solidum import Assembler, Domain
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.solvers import (
    CentralDifferenceSolver,
    HHTSolver,
    NewmarkSolver,
    NewtonHHTSolver,
    NewtonNewmarkSolver,
)
from solidum.math.solvers.diagnostics import (
    MasslessDofsError,
    SolverDivergedError,
    Termination,
    TimeIntegrationInstabilityError,
)

E, A, RHO = 2.0e11, 1.0e-3, 7850.0


def bar(n=3, material=None, density=RHO):
    """Barra horizontal de ``n`` elementos, fija en el nodo 1; sólo ``ux``
    libre. Devuelve dominio y carga en el extremo."""
    material = material if material is not None else Elastic1D(E=E, density=density)
    dom = Domain()
    nodes = [dom.add_node(k + 1, [float(k), 0.0]) for k in range(n + 1)]
    for k in range(n):
        dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]], material, A))
    nodes[0].fix_dof("ux", 0.0)
    for node in nodes:
        node.fix_dof("uy", 0.0)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    F[nodes[-1].dofs["ux"]] = 1.0e4
    return dom, F


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.ERROR)
        self.addCleanup(logging.disable, logging.NOTSET)


class TestLlegaAlFinal(_Quiet):

    def test_todos_los_integradores(self):
        for cls, kw in ((NewmarkSolver, {}), (HHTSolver, {}), (NewtonNewmarkSolver, {}),
                        (NewtonHHTSolver, {}), (CentralDifferenceSolver, {"lumping": "lumped"})):
            with self.subTest(solver=cls.__name__):
                dom, F = bar()
                solver = cls(Assembler(dom), t_end=1.0e-3, dt=1.0e-5,
                             F_func=lambda t, F=F: F, **kw)
                self.assertIsNone(solver.termination)
                res = solver.solve()
                self.assertEqual(solver.termination, Termination.T_END)
                self.assertTrue(solver.completed)
                self.assertEqual(solver.steps_done, 100)
                self.assertAlmostEqual(solver.t_final, 1.0e-3, delta=1e-15)
                np.testing.assert_array_equal(solver.U_final, res.u_history[:, -1])
                self.assertIsNone(solver.partial_result)
                self.assertEqual(res.termination, "t_end")
                self.assertTrue(res.converged)
                self.assertEqual(res.n_steps, 100)


class TestSeDetiene(_Quiet):

    def _newton_failing_at(self, cls, t_fail):
        """Newton que deja de converger a partir de ``t_fail`` (sin
        iteraciones disponibles): la función de carga lo provoca."""
        dom, F = bar(material=Elastoplastic1D(E=E, sigma_y=2.5e8, H=1.0e9, density=RHO))
        solver = cls(Assembler(dom), t_end=1.0e-3, dt=1.0e-5,
                     F_func=lambda t: F * (1.0 + 100.0 * t))

        def F_func(t):
            if t > t_fail:
                solver.corrector.max_iter = 0
            return F * (1.0 + 100.0 * t)

        solver.F_func = F_func
        return solver

    def test_newton_sin_convergencia(self):
        for cls in (NewtonNewmarkSolver, NewtonHHTSolver):
            with self.subTest(solver=cls.__name__):
                solver = self._newton_failing_at(cls, 0.5e-3)
                logging.disable(logging.NOTSET)
                with self.assertLogs("solidum", level=logging.ERROR) as logs, \
                        self.assertRaises(SolverDivergedError) as info:
                    solver.solve()
                self.assertEqual(sum(r.levelno == logging.ERROR for r in logs.records), 1)
                parcial = info.exception.partial_result
                self.assertIs(parcial, solver.partial_result)
                self.assertEqual(solver.termination, Termination.DIVERGED)
                self.assertFalse(solver.completed)
                self.assertEqual(solver.steps_done, 50)
                self.assertAlmostEqual(solver.t_final, 0.5e-3, delta=1e-15)
                self.assertEqual(parcial.termination, "diverged")
                self.assertFalse(parcial.converged)
                self.assertEqual(parcial.u_history.shape[1], 51)
                self.assertEqual(parcial.n_steps, 50)
                np.testing.assert_array_equal(parcial.u_history[:, -1], solver.U_final)
                self.assertGreater(np.max(np.abs(parcial.u_history[:, -1])), 0.0)
                self.assertGreater(info.exception.last_load_factor, solver.t_final)

    def test_diferencia_central_inestable(self):
        # Elementos de 1 m de acero: Δt_crit ≈ L/c ≈ 2e-4 s; Δt = 2.5·Δt_crit.
        dom, F = bar(n=10)
        solver = CentralDifferenceSolver(Assembler(dom), t_end=5.0e-2, dt=5.0e-4,
                                         F_func=lambda t: F, lumping="lumped")
        with self.assertRaises(TimeIntegrationInstabilityError) as info:
            solver.solve()
        exc = info.exception
        self.assertIsInstance(exc, RuntimeError)
        self.assertRegex(str(exc), "divergencia.*CFL")
        parcial = exc.partial_result
        self.assertEqual(solver.termination, Termination.DIVERGED)
        self.assertTrue(np.all(np.isfinite(parcial.u_history)))
        self.assertEqual(parcial.u_history.shape[1], solver.steps_done + 1)
        self.assertLess(solver.steps_done, 100)

    def test_un_solver_reutilizado_reinicia_su_estado(self):
        solver = self._newton_failing_at(NewtonNewmarkSolver, 0.2e-3)
        with self.assertRaises(SolverDivergedError):
            solver.solve()
        dom, F = bar()
        solver.assembler = Assembler(dom)
        solver.corrector.max_iter = 20
        solver.F_func = lambda t: F
        solver.solve()
        self.assertEqual(solver.termination, Termination.T_END)
        self.assertIsNone(solver.partial_result)


class TestSinMasa(_Quiet):

    def _model(self):
        """Barra sin densidad con masa de nudo sólo en el extremo: el nodo 2
        no tiene masa en ``ux``."""
        dom, F = bar(n=2, density=0.0)
        dom.add_nodal_mass(3, ux=50.0)
        return dom, F

    def test_implicitos(self):
        for cls in (NewmarkSolver, HHTSolver, NewtonNewmarkSolver, NewtonHHTSolver):
            with self.subTest(solver=cls.__name__):
                dom, F = self._model()
                solver = cls(Assembler(dom), t_end=1.0e-3, dt=1.0e-4, F_func=lambda t: F)
                with self.assertRaises(MasslessDofsError) as info:
                    solver.solve()
                self.assertEqual(info.exception.cause, "initial_acceleration")
                self.assertEqual(info.exception.massless_dofs, ["nodo 2 (ux)"])
                self.assertIn("deuda #34", str(info.exception))
                self.assertIn("add_nodal_mass", str(info.exception))

    def test_diferencia_central(self):
        dom, F = self._model()
        solver = CentralDifferenceSolver(Assembler(dom), t_end=1.0e-3, dt=1.0e-6,
                                         F_func=lambda t: F, lumping="lumped")
        with self.assertRaises(MasslessDofsError) as info:
            solver.solve()
        self.assertEqual(info.exception.cause, "explicit_integration")
        self.assertEqual(info.exception.massless_dofs, ["nodo 2 (ux)"])


class TestApiPublica(unittest.TestCase):

    def test_entrypoints_exportados(self):
        for name in ("run_modal", "run_transient", "run_harmonic", "run_response_spectrum",
                     "run_thermal_transient", "ModalResult", "TransientResult",
                     "HarmonicResult", "ResponseSpectrumResult", "LoadVector"):
            self.assertTrue(hasattr(solidum, name), name)


if __name__ == "__main__":
    unittest.main()
