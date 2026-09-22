"""Regresiones del lote estático de la auditoría 2026-09-22.

- Un ensamblaje por iteración de Newton (antes dos).
- Arc-length: estado committed coherente con ``U_current`` (antes iba un
  iterado por detrás); terminación por ``max_steps`` reflejada en
  ``SolveResult.converged``; ``F_applied`` y reacciones escaladas por el
  factor de carga alcanzado; ``num_steps`` real.
"""
import logging
import unittest

import numpy as np

import solidum
from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.von_mises_2d import VonMises2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import ArcLengthSolver, NonlinearSolver


def _two_quads(material, thickness=0.5):
    dom = Domain()
    for i, (x, y) in enumerate([(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)], 1):
        dom.add_node(i, [float(x), float(y)])
    dom.add_element(Quad4(1, [dom.nodes[i] for i in (1, 2, 5, 4)], material, thickness=thickness))
    dom.add_element(Quad4(2, [dom.nodes[i] for i in (2, 3, 6, 5)], material, thickness=thickness))
    dom.generate_equation_numbers()
    for k in (1, 4):
        dom.nodes[k].fix_dof('ux'); dom.nodes[k].fix_dof('uy')
    F = np.zeros(dom.total_dofs)
    F[dom.nodes[3].dofs['ux']] = 8.0
    F[dom.nodes[6].dofs['ux']] = 8.0
    return dom, F


class _Instrumented:
    """Envuelve el Assembler para contar ensamblajes y comparar el U del
    último ensamblaje antes de cada commit con el U que el solver guarda."""

    def __init__(self, asm):
        self.asm = asm
        self.n_assemblies = 0
        self._last_U = None
        self.assembled_at_commit = []
        orig_assemble = asm.assemble_non_linear_system
        orig_commit = asm.commit_all_states

        def assemble(U):
            self.n_assemblies += 1
            self._last_U = np.array(U, copy=True)
            return orig_assemble(U)

        def commit():
            self.assembled_at_commit.append(np.array(self._last_U, copy=True))
            orig_commit()

        asm.assemble_non_linear_system = assemble
        asm.commit_all_states = commit


class TestOneAssemblyPerIteration(unittest.TestCase):

    def test_newton_assembles_once_per_solve_plus_one_per_step(self):
        dom, F = _two_quads(VonMises2D(E=1000.0, nu=0.3, sigma_y=5.0, H=50.0))
        asm = Assembler(dom)
        inst = _Instrumented(asm)
        solves = {"n": 0}
        solver = NonlinearSolver(asm, num_steps=5, adaptive=False)
        orig = solver._solve_reduced

        def counted(K, R, iteration=0):
            solves["n"] += 1
            return orig(K, R, iteration=iteration)
        solver._solve_reduced = counted
        solver.solve(F)
        # Un ensamblaje inicial por paso + uno por cada resolución.
        self.assertEqual(inst.n_assemblies, 5 + solves["n"])
        self.assertGreater(solves["n"], 5)


class TestCommittedStateConsistency(unittest.TestCase):

    def _check(self, make_solver):
        dom, F = _two_quads(VonMises2D(E=1000.0, nu=0.3, sigma_y=5.0, H=50.0))
        asm = Assembler(dom)
        inst = _Instrumented(asm)
        stored = []
        make_solver(asm).solve(F, step_callback=lambda s, U, lam: stored.append(np.array(U, copy=True)))
        self.assertEqual(len(stored), len(inst.assembled_at_commit))
        for U_stored, U_assembled in zip(stored, inst.assembled_at_commit):
            self.assertLess(float(np.max(np.abs(U_stored - U_assembled))), 1e-15)

    def test_newton(self):
        self._check(lambda asm: NonlinearSolver(asm, num_steps=6))

    def test_arclength(self):
        self._check(lambda asm: ArcLengthSolver(asm, initial_dl=0.02, max_lambda=1.0, max_steps=200))


class TestArcLengthResultMetadata(unittest.TestCase):

    def test_stopped_by_max_steps_is_reported(self):
        dom, F = _two_quads(Elastic2D(E=1000.0, nu=0.3))
        asm = Assembler(dom)
        solver = ArcLengthSolver(asm, initial_dl=1e-4, max_steps=3, max_lambda=1.0)
        with self.assertLogs("solidum.solvers", level=logging.WARNING) as cm:
            res = solidum.run(dom, assembler=asm, solver=solver, F_applied=F)
        self.assertTrue(any("max_steps" in m for m in cm.output))
        self.assertFalse(res.converged)
        self.assertFalse(solver.reached_max_lambda)
        self.assertEqual(res.num_steps, 3)
        self.assertLess(solver.lambda_final, 0.5)
        # Cargas realmente aplicadas = λ_final · F_ref ⇒ reacciones coherentes.
        np.testing.assert_allclose(res.F_applied, F * solver.lambda_final)
        self.assertAlmostEqual(float(np.sum(res.R)), -float(np.sum(res.F_applied)), places=8)

    def test_max_lambda_two_scales_loads_and_reactions(self):
        dom, F = _two_quads(Elastic2D(E=1000.0, nu=0.3))
        asm = Assembler(dom)
        solver = ArcLengthSolver(asm, initial_dl=0.05, max_steps=500, max_lambda=2.0)
        res = solidum.run(dom, assembler=asm, solver=solver, F_applied=F)
        self.assertTrue(res.converged)
        self.assertAlmostEqual(solver.lambda_final, 2.0, places=9)
        np.testing.assert_allclose(res.F_applied, 2.0 * F)
        # Equilibrio global en x: ΣR_x = −ΣF_x = −32.
        rx = sum(v['ux'] for v in res.reactions_by_node.values())
        self.assertAlmostEqual(rx, -32.0, places=6)
        self.assertEqual(res.num_steps, solver.steps_done)
        self.assertGreater(res.num_steps, 1)

    def test_newton_reports_converged_steps(self):
        dom, F = _two_quads(Elastic2D(E=1000.0, nu=0.3))
        asm = Assembler(dom)
        res = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=4), F_applied=F)
        self.assertTrue(res.converged)
        self.assertEqual(res.num_steps, 4)
        np.testing.assert_allclose(res.F_applied, F)


if __name__ == "__main__":
    unittest.main()
