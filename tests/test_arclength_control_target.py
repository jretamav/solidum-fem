"""Parada del arco cilíndrico en el valor de una magnitud (``control`` y
``control_target`` de ``ArcLengthSolver``; ADR 0025).

En el arco cilíndrico ``control`` sólo mide: la restricción del paso sigue
siendo la norma del incremento. El paso que cruza el objetivo se repite con
la restricción del control indirecto. Criterios:

1. Arco de von Mises trazado más allá del pico (salto dinámico que el control
   por carga no sigue) hasta un desplazamiento del vértice: el último punto
   coincide con el objetivo a precisión de máquina con cualquier longitud de
   arco, está en equilibrio, y es el mismo punto que da el control indirecto
   con el mismo objetivo.
2. Voladizo con plasticidad en las fibras empujado hasta un desplazamiento de
   cabeza: el mismo punto que el control indirecto.
3. ``max_lambda`` antes que el objetivo: se detiene en ``max_lambda``.
4. En una etapa el objetivo es absoluto: parte del valor en ``u0``.
5. Desde el YAML, y errores de uso claros.
"""
import logging
import textwrap
import unittest

import numpy as np

import solidum
from solidum import Assembler, load_yaml
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import ArcLengthSolver, IndirectDisplacementSolver, NonlinearSolver
from solidum.math.solvers.diagnostics import Termination
from test_arclength_paso_inicial import _arco_von_mises
from test_trace_termination import F_REF, cantilever

EPS = np.finfo(float).eps


def _criterion():
    return ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10)


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


def _apex(dom):
    return dom.nodes[2]


class TestArcoDeVonMises(_Quiet):
    """Criterio 1."""

    TARGET = -1.5          # vértice 1.5·h por debajo de su posición inicial

    def _trace(self, solver_cls, **kw):
        dom, F = _arco_von_mises()
        asm = Assembler(dom)
        apex = _apex(dom)
        solver = solver_cls(asm, convergence=_criterion(), control=[(apex.id, "uy", 1.0)],
                            control_target=self.TARGET, max_lambda=10.0, max_steps=400,
                            max_iter=30, **kw)
        U = solver.solve(F)
        return dom, asm, F, solver, U

    def test_llegada_exacta_y_mismo_punto_que_el_control_indirecto(self):
        _, _, _, ref, U_ref = self._trace(IndirectDisplacementSolver, initial_dl=0.05)
        self.assertEqual(ref.termination, Termination.CONTROL_TARGET)
        for dl in (0.07, 0.2, 0.45):
            with self.subTest(initial_dl=dl):
                dom, asm, F, solver, U = self._trace(ArcLengthSolver, initial_dl=dl)
                apex = _apex(dom)
                self.assertLessEqual(abs(U[apex.dofs["uy"]] - self.TARGET),
                                     8 * EPS * abs(self.TARGET))
                self.assertEqual(solver.termination, Termination.CONTROL_TARGET)
                self.assertTrue(solver.reached_control_target)
                self.assertTrue(solver.completed)
                self.assertFalse(solver.reached_max_lambda)
                self.assertEqual(solver.control_final, U[apex.dofs["uy"]])
                # Pasado el salto, la carga se invierte para sostener el vértice.
                self.assertLess(solver.lambda_final, 0.0)
                # Equilibrio del último punto.
                _, F_int = asm.assemble_non_linear_system(U)
                free = asm.constraint_set.free_dofs(dom.total_dofs)
                R = solver.lambda_final * F - F_int
                self.assertLess(np.linalg.norm(R[free]), 1e-8 * np.linalg.norm(F))
                # El mismo punto de equilibrio que el control indirecto.
                self.assertAlmostEqual(solver.lambda_final / ref.lambda_final, 1.0, delta=1e-8)
                np.testing.assert_allclose(U, U_ref, rtol=0.0, atol=1e-9)

    def test_max_lambda_antes_que_el_objetivo(self):
        # Con esta carga de referencia el pico está en λ ≈ 3.7e-3: max_lambda
        # por debajo de él se alcanza antes que el desplazamiento objetivo.
        dom, F = _arco_von_mises()
        solver = ArcLengthSolver(Assembler(dom), convergence=_criterion(),
                                 control=[(2, "uy", 1.0)], control_target=self.TARGET,
                                 max_lambda=2.0e-3, initial_dl=0.2, max_steps=400)
        U = solver.solve(F)
        self.assertEqual(solver.termination, Termination.MAX_LAMBDA)
        self.assertFalse(solver.reached_control_target)
        self.assertAlmostEqual(solver.lambda_final, 2.0e-3, places=14)
        self.assertGreater(U[dom.nodes[2].dofs["uy"]], self.TARGET)


class TestVoladizoPlastico(_Quiet):
    """Criterio 2."""

    def test_igual_que_el_control_indirecto(self):
        target = 0.03
        results = {}
        for cls, dl in ((IndirectDisplacementSolver, 2.0e-3), (ArcLengthSolver, 3.3e-3)):
            dom, tip, F = cantilever("plastic")
            solver = cls(Assembler(dom), convergence=_criterion(), control=[(tip.id, "ux", 1.0)],
                         control_target=target, initial_dl=dl, max_lambda=1.0e9,
                         max_steps=500)
            U = solver.solve(F)
            self.assertEqual(solver.termination, Termination.CONTROL_TARGET)
            self.assertLessEqual(abs(U[tip.dofs["ux"]] - target), 8 * EPS * target)
            results[cls] = solver.lambda_final
        # Más allá de la fluencia (λ·F_REF por encima del momento elástico).
        self.assertAlmostEqual(results[ArcLengthSolver] / results[IndirectDisplacementSolver],
                               1.0, delta=1e-7)


class TestEtapa(_Quiet):
    """Criterio 4."""

    def test_objetivo_absoluto_desde_u0(self):
        dom, tip, F = cantilever("elastic")
        asm = Assembler(dom)
        r1 = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=2),
                         F_applied=0.5 * F)
        u1 = r1.U[tip.dofs["ux"]]
        solver = ArcLengthSolver(asm, convergence=_criterion(), control=[(tip.id, "ux", 1.0)],
                                 control_target=2.0 * u1, initial_dl=0.37 * u1,
                                 max_lambda=1.0e9, max_steps=100)
        r2 = solidum.run(dom, assembler=asm, solver=solver, F_applied=F, previous=r1)
        self.assertEqual(r2.termination, "control_target")
        self.assertTrue(r2.converged)
        self.assertLessEqual(abs(r2.U[tip.dofs["ux"]] - 2.0 * u1), 8 * EPS * u1)
        # Elástico lineal: el desplazamiento adicional u1 lo da media carga.
        self.assertAlmostEqual(r2.load_factor, 0.5, delta=1e-9)


class TestYamlYErrores(_Quiet):
    """Criterio 5."""

    def test_yaml(self):
        import tempfile, os
        text = textwrap.dedent(f"""\
            nodes:
              - {{id: 1, coords: [0.0, 0.0]}}
              - {{id: 2, coords: [0.0, 1.5]}}
              - {{id: 3, coords: [0.0, 3.0]}}
            materials:
              - {{id: 1, type: Elastoplastic1D, E: 2.0e11, sigma_y: 2.5e8, H: 2.0e9}}
            sections:
              - id: 1
                type: FiberSection
                patches:
                  - {{shape: rectangle, material: 1, y: [-0.25, 0.25], z: [-0.15, 0.15], ny: 20, nz: 1}}
            elements:
              - {{id: 1, type: Frame2DEuler, nodes: [1, 2], section: 1}}
              - {{id: 2, type: Frame2DEuler, nodes: [2, 3], section: 1}}
            boundary_conditions_by_node:
              - {{node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}}
            point_loads:
              - {{node_id: 3, ux: {F_REF}}}
            solver:
              type: ArcLengthSolver
              control:
                - {{node: 3, dof: ux, coef: 1.0}}
              control_target: 0.04
              initial_dl: 0.006
              max_lambda: 1.0e9
              max_steps: 500
            """)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "modelo.yaml")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text)
            model = load_yaml(path)
            res = model.solve()
        self.assertEqual(res.termination, "control_target")
        self.assertLessEqual(abs(res.U[model.domain.nodes[3].dofs["ux"]] - 0.04), 8 * EPS * 0.04)

    def test_errores(self):
        dom, tip, F = cantilever("elastic")
        asm = Assembler(dom)
        with self.assertRaisesRegex(ValueError, "control_target necesita 'control'"):
            ArcLengthSolver(asm, control_target=0.01)
        with self.assertRaisesRegex(ValueError, "debe ser un número finito"):
            ArcLengthSolver(asm, control=[(tip.id, "ux", 1.0)], control_target=float("nan"))
        with self.assertRaisesRegex(ValueError, "lista no vacía"):
            ArcLengthSolver(asm, control=[], control_target=0.01)
        base = dom.nodes[1]
        with self.assertRaisesRegex(ValueError, "ArcLengthSolver: el grado de libertad 'ux' del "
                                                "nodo 1 está prescrito"):
            ArcLengthSolver(asm, control=[(base.id, "ux", 1.0)], control_target=0.01).solve(F)

    def test_sin_objetivo_control_solo_mide(self):
        """``control`` sin objetivo no cambia el trazado: mismos puntos que sin él."""
        paths = []
        for kw in ({}, {"control": [(5, "ux", 1.0)]}):
            dom, tip, F = cantilever("plastic")
            solver = ArcLengthSolver(Assembler(dom), max_lambda=1.2, initial_dl=4.0e-3,
                                     max_steps=200, **kw)
            lams = []
            solver.solve(F, step_callback=lambda k, U, lam: lams.append(lam))
            paths.append(np.array(lams))
        np.testing.assert_array_equal(paths[0], paths[1])
        self.assertEqual(solver.termination, Termination.MAX_LAMBDA)
        self.assertNotEqual(solver.control_final, 0.0)


if __name__ == "__main__":
    unittest.main()
