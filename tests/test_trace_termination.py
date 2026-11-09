"""Terminación de los solvers estáticos incrementales (ADR 0025).

Parada en un desplazamiento objetivo (``IndirectDisplacementSolver``,
``control_target``):

- Voladizo elástico empujado hasta un desplazamiento objetivo: el último
  punto coincide con el objetivo a precisión de máquina, con cualquier tamaño
  de paso, no sólo cuando el objetivo es múltiplo del paso. Un objetivo que se
  alcanza en medio de un paso reduce ese paso en vez de sobrepasarlo.
- ``termination``, ``reached_control_target``, ``completed`` y
  ``SolveResult.converged`` coherentes con haber llegado al objetivo.

Estado legible y excepción tipada del arco:

- Un trazado que agota ``max_steps`` y otro que fracasa por paso
  mínimo dejan, ambos, un estado legible y distinguible sin interpretar
  mensajes (``termination``).
- El fracaso del arco en la barra con daño localizado lanza una
  subclase de ``SolverDivergedError`` cuyo diagnóstico menciona la causa.
- Un solver reutilizado tras un fracaso no muestra los valores de la
  corrida anterior.
- Ninguna prueba existente cambia de resultado (la suite completa).
"""
import logging
import math
import os
import tempfile
import unittest
from unittest import mock

import numpy as np

import solidum
from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler
from solidum.elements.truss import Truss2D
from solidum.materials.damage_1d import IsotropicDamage1D
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.solvers import ArcLengthSolver, IndirectDisplacementSolver, NonlinearSolver
from solidum.math.solvers.arclength import ArcLengthSolver as _Arc
from solidum.math.solvers.diagnostics import (
    MechanismError,
    SingularTangentError,
    SolverDivergedError,
    StepConstraintError,
    Termination,
)
from solidum.sections import FiberSection, RectangularPatch

E, L, B, H = 2.0e11, 3.0, 0.30, 0.50
I_SEC = B * H ** 3 / 12.0
F_REF = 1.0e6
EPS = np.finfo(float).eps


def cantilever(material="elastic", n_el=4, hardening=2.0e9):
    """Columna vertical en voladizo con carga lateral de referencia en la
    cabeza. ``material="plastic"``: sección de 20 fibras de acero bilineal
    (por omisión, endurecimiento del 1 %; ``hardening=0`` es perfectamente
    plástico)."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, [0.0, i * L / n_el]) for i in range(n_el + 1)]
    for k in range(n_el):
        if material == "elastic":
            elem = Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]],
                                material=Elastic1D(E=E), A=B * H, I=I_SEC)
        else:
            steel = Elastoplastic1D(E=E, sigma_y=250.0e6, H=hardening)
            section = FiberSection(patches=[RectangularPatch(
                steel, y=(-H / 2, H / 2), z=(-B / 2, B / 2), ny=20, nz=1)])
            elem = Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], section=section, n_ip=3)
        dom.add_element(elem)
    for dof in ("ux", "uy", "rz"):
        nodes[0].fix_dof(dof, 0.0)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    F[nodes[-1].dofs["ux"]] = F_REF
    return dom, nodes[-1], F


def damaged_bar(n=10):
    """Barra con daño en ablandamiento y el elemento central un 5 % más débil
    (el modelo de ``test_indirect_displacement.py``): tiene retroceso, y el
    arco cilíndrico no lo pasa. Carga de referencia = carga pico."""
    E_bar, A, length, k0, alpha, weak = 30.0e9, 0.01, 1.0, 1.0e-4, 5000.0, 0.95
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * length / n, 0.0]) for i in range(n + 1)]
    for k in range(n):
        kappa = weak * k0 if k == n // 2 else k0
        dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]],
                                IsotropicDamage1D(E=E_bar, kappa_0=kappa, alpha=alpha), A=A))
    nodes[0].fix_dof("ux", 0.0)
    for node in nodes:
        node.fix_dof("uy", 0.0)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    F[nodes[-1].dofs["ux"]] = A * E_bar * weak * k0
    return dom, F


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


class TestParadaEnDesplazamientoObjetivo(_Quiet):
    """Llegada exacta con ``control_target``."""

    def _push(self, target, *, model="elastic", **kw):
        dom, tip, F = cantilever(model)
        kw.setdefault("max_lambda", 1.0e9)
        solver = IndirectDisplacementSolver(
            Assembler(dom), control=[(tip.id, "ux", 1.0)], control_target=target, **kw)
        path = []
        U = solver.solve(F, step_callback=lambda k, U, lam: path.append((U[tip.dofs["ux"]], lam)))
        return solver, U, tip, np.array(path)

    def test_llegada_exacta_con_cualquier_paso(self):
        target = 7.5e-3
        delta_unit = F_REF * L ** 3 / (3.0 * E * I_SEC)      # flecha con λ = 1
        for dl in (1.1e-3, 2.0e-3, 3.7e-3, 1.0e-2):
            with self.subTest(initial_dl=dl):
                solver, U, tip, path = self._push(target, initial_dl=dl)
                self.assertLessEqual(abs(U[tip.dofs["ux"]] - target), 4 * EPS * target)
                self.assertEqual(solver.termination, Termination.CONTROL_TARGET)
                self.assertEqual(solver.termination, "control_target")
                self.assertTrue(solver.reached_control_target)
                self.assertTrue(solver.completed)
                self.assertFalse(solver.reached_max_lambda)
                self.assertEqual(solver.control_final, U[tip.dofs["ux"]])
                np.testing.assert_array_equal(solver.U_final, U)
                # Un paso por cada dl entero y el último recortado a lo que falta.
                self.assertEqual(solver.steps_done, math.ceil(target / dl - 1e-9))
                increments = np.diff(np.concatenate([[0.0], path[:, 0]]))
                self.assertTrue(np.all(increments <= dl * (1 + 1e-12)))
                self.assertAlmostEqual(path[-1, 0], target, delta=4 * EPS * target)
                # Elástico: el factor de carga es el de la flecha objetivo.
                self.assertAlmostEqual(solver.lambda_final * delta_unit / target, 1.0, delta=1e-9)

    def test_objetivo_multiplo_del_paso_sin_paso_espurio(self):
        solver, U, tip, path = self._push(7.5e-3, initial_dl=2.5e-3)
        self.assertEqual(solver.steps_done, 3)
        self.assertLessEqual(abs(U[tip.dofs["ux"]] - 7.5e-3), 4 * EPS * 7.5e-3)

    def test_primer_paso_adimensional(self):
        solver, U, tip, _ = self._push(9.0e-3, initial_dlambda=0.3)
        self.assertLessEqual(abs(U[tip.dofs["ux"]] - 9.0e-3), 4 * EPS * 9.0e-3)
        self.assertEqual(solver.termination, Termination.CONTROL_TARGET)

    def test_meseta_plastica(self):
        """En la meseta plástica λ casi no cambia y ``max_lambda`` no sirve
        para detenerse; el objetivo sí, con llegada exacta."""
        target = 6.0e-2
        solver, U, tip, path = self._push(target, model="plastic", initial_dl=4.5e-3, max_steps=60)
        self.assertEqual(solver.termination, Termination.CONTROL_TARGET)
        self.assertLessEqual(abs(U[tip.dofs["ux"]] - target), 4 * EPS * target)
        # Meseta: en los últimos cuatro pasos λ varía menos de un 5 %.
        lam = path[:, 1]
        self.assertLess((lam[-1] - lam[-5]) / lam[-1], 0.05)

    def test_max_lambda_antes_que_el_objetivo(self):
        solver, U, tip, _ = self._push(1.0, initial_dl=2.0e-3, max_lambda=0.5)
        self.assertEqual(solver.termination, Termination.MAX_LAMBDA)
        self.assertTrue(solver.reached_max_lambda)
        self.assertFalse(solver.reached_control_target)
        self.assertAlmostEqual(solver.lambda_final, 0.5, delta=1e-12)

    def test_objetivo_invalido(self):
        dom, tip, _ = cantilever()
        for bad in (0.0, -1.0e-3, math.nan, math.inf):
            with self.subTest(control_target=bad):
                with self.assertRaisesRegex(ValueError, "control_target"):
                    IndirectDisplacementSolver(Assembler(dom), control=[(tip.id, "ux", 1.0)],
                                               control_target=bad)

    def test_solidum_run_converged_y_termination(self):
        dom, tip, F = cantilever()
        solver = IndirectDisplacementSolver(Assembler(dom), control=[(tip.id, "ux", 1.0)],
                                            control_target=5.0e-3, max_lambda=1.0e9,
                                            initial_dl=2.0e-3)
        res = solidum.run(dom, assembler=solver.assembler, solver=solver, F_applied=F)
        self.assertTrue(res.converged)
        self.assertEqual(res.termination, "control_target")
        self.assertLessEqual(abs(res.U[tip.dofs["ux"]] - 5.0e-3), 4 * EPS * 5.0e-3)

    def test_yaml(self):
        n = 4
        nodes = "\n".join(f"  - {{id: {i + 1}, coords: [0.0, {i * L / n!r}]}}" for i in range(n + 1))
        elems = "\n".join(f"  - {{id: {k + 1}, type: Frame2DEuler, material: 1, nodes: [{k + 1}, {k + 2}], "
                          f"A: {B * H!r}, I: {I_SEC!r}}}" for k in range(n))
        text = (f"nodes:\n{nodes}\nmaterials:\n  - {{id: 1, type: Elastic1D, E: {E!r}}}\n"
                f"elements:\n{elems}\nboundary_conditions:\n  - {{node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}}\n"
                f"point_loads:\n  - {{node_id: {n + 1}, ux: {F_REF!r}}}\n"
                "solver:\n  type: IndirectDisplacementSolver\n  max_lambda: 1.0e9\n"
                "  initial_dl: 2.0e-3\n  control_target: 7.0e-3\n"
                f"  control:\n    - {{node: {n + 1}, dof: ux}}\n")
        fd, path = tempfile.mkstemp(suffix=".yaml")
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(text)
        self.addCleanup(os.remove, path)
        res = solidum.run_yaml(path)
        self.assertEqual(res.termination, "control_target")
        self.assertTrue(res.converged)
        # Elástico: la carga final es la que da la flecha objetivo,
        # F = 3EI·δ/L³, y la reacción horizontal la equilibra.
        F_exact = 3.0 * E * I_SEC * 7.0e-3 / L ** 3
        self.assertAlmostEqual(float(np.sum(res.F_applied)) / F_exact, 1.0, delta=1e-9)
        self.assertAlmostEqual(res.reactions_by_node[1]["ux"] / -F_exact, 1.0, delta=1e-9)


class TestEstadoYExcepcionDelArco(_Quiet):
    """Estado legible, excepción tipada y reinicio en cada ``solve``."""

    def test_max_steps_y_fracaso_son_distinguibles(self):
        dom, tip, F = cantilever()
        exhausted = IndirectDisplacementSolver(Assembler(dom), control=[(tip.id, "ux", 1.0)],
                                               max_lambda=1.0e9, max_steps=3, initial_dl=1.0e-3)
        U = exhausted.solve(F)
        self.assertEqual(exhausted.termination, Termination.MAX_STEPS)
        self.assertFalse(exhausted.completed)
        self.assertEqual(exhausted.steps_done, 3)
        np.testing.assert_array_equal(exhausted.U_final, U)
        res = solidum.run(dom, assembler=Assembler(dom),
                          solver=IndirectDisplacementSolver(
                              Assembler(dom), control=[(tip.id, "ux", 1.0)], max_lambda=1.0e9,
                              max_steps=3, initial_dl=1.0e-3),
                          F_applied=F)
        self.assertFalse(res.converged)
        self.assertEqual(res.termination, "max_steps")

        dom_b, F_b = damaged_bar()
        failed = ArcLengthSolver(Assembler(dom_b), max_iter=30, max_lambda=1.2,
                                 max_steps=400, initial_dlambda=0.02)
        with self.assertRaises(SolverDivergedError):
            failed.solve(F_b)
        self.assertEqual(failed.termination, Termination.DIVERGED)
        self.assertFalse(failed.completed)
        self.assertNotEqual(exhausted.termination, failed.termination)

    def test_fracaso_del_arco_tipado_con_su_causa(self):
        dom, F = damaged_bar()
        solver = ArcLengthSolver(Assembler(dom), max_iter=30, max_lambda=1.2,
                                 max_steps=400, initial_dlambda=0.02)
        with self.assertRaises(StepConstraintError) as ctx:
            solver.solve(F)
        err = ctx.exception
        self.assertIsInstance(err, SolverDivergedError)
        self.assertIsInstance(err, RuntimeError)          # compatibilidad
        self.assertEqual(err.mode, "step_constraint")
        self.assertIn("Raíces imaginarias", str(err))
        self.assertGreater(err.n_bisections, 0)

    def test_solver_reutilizado_no_conserva_la_corrida_anterior(self):
        dom, F = damaged_bar()
        solver = ArcLengthSolver(Assembler(dom), max_iter=30, max_lambda=0.5,
                                 max_steps=400, initial_dlambda=0.02)
        solver.solve(F)
        first = (solver.termination, solver.lambda_final, solver.steps_done)
        self.assertEqual(first[0], Termination.MAX_LAMBDA)

        solver.max_lambda = 1.2
        seen = []
        with self.assertRaises(StepConstraintError):
            solver.solve(F, step_callback=lambda k, U, lam: seen.append((lam, U.copy())))
        self.assertEqual(solver.termination, Termination.DIVERGED)
        self.assertFalse(solver.reached_max_lambda)
        # Hasta dónde llegó la corrida que fracasó, no la anterior.
        self.assertEqual(solver.steps_done, len(seen))
        self.assertEqual(solver.lambda_final, seen[-1][0])
        np.testing.assert_array_equal(solver.U_final, seen[-1][1])
        self.assertNotEqual((solver.lambda_final, solver.steps_done), first[1:])

    def test_reinicio_aunque_falle_antes_del_primer_paso(self):
        dom, tip, F = cantilever()
        solver = IndirectDisplacementSolver(Assembler(dom), control=[(tip.id, "ux", 1.0)],
                                            max_lambda=1.0e9, max_steps=3, initial_dl=1.0e-3)
        solver.solve(F)
        self.assertEqual(solver.steps_done, 3)
        base = dom.nodes[1]
        for dof in ("ux", "uy", "rz"):
            base.boundary_conditions.pop(dof)     # sin apoyos: mecanismo
        with self.assertRaises(MechanismError):
            solver.solve(F)
        self.assertIsNone(solver.termination)
        self.assertEqual((solver.steps_done, solver.lambda_final, solver.n_bisections), (0, 0.0, 0))
        self.assertIsNone(solver.U_final)

    def test_tangente_singular_del_estado_inicial(self):
        """Sin paso previo no hay longitud que reducir: antes el bucle repetía
        el mismo predictor hasta agotar ``max_steps`` y lo informaba como un
        trazado incompleto; ahora es un fracaso tipado."""
        dom, tip, F = cantilever()
        solver = ArcLengthSolver(Assembler(dom), max_steps=50)
        with mock.patch.object(_Arc, "_tangent_predictor", return_value=None):
            with self.assertRaises(SingularTangentError):
                solver.solve(F)
        self.assertEqual(solver.termination, Termination.DIVERGED)
        self.assertEqual(solver.steps_done, 0)


class TestNonlinearSolver(_Quiet):
    """El mismo contrato en el Newton incremental."""

    def test_completo(self):
        dom, tip, F = cantilever("plastic")
        solver = NonlinearSolver(Assembler(dom), num_steps=5)
        res = solidum.run(dom, assembler=solver.assembler, solver=solver, F_applied=0.9 * F)
        self.assertEqual(solver.termination, Termination.MAX_LAMBDA)
        self.assertTrue(solver.completed)
        self.assertEqual(res.termination, "max_lambda")
        np.testing.assert_array_equal(solver.U_final, res.U)

    def test_fracaso_deja_el_ultimo_paso_convergido(self):
        """Acero perfectamente plástico y carga por encima de la capacidad de
        la sección (F_p = M_p/L ≈ 1.56 MN): el control de carga no la alcanza."""
        dom, tip, F = cantilever("plastic", hardening=0.0)
        solver = NonlinearSolver(Assembler(dom), num_steps=10, max_iter=15, min_delta_lambda=1e-4)
        seen = []
        with self.assertRaises(SolverDivergedError):
            solver.solve(5.0 * F, step_callback=lambda k, U, lam: seen.append((lam, U.copy())))
        self.assertEqual(solver.termination, Termination.DIVERGED)
        self.assertFalse(solver.completed)
        self.assertEqual(solver.lambda_final, seen[-1][0])
        self.assertEqual(solver.steps_done, len(seen))
        np.testing.assert_array_equal(solver.U_final, seen[-1][1])
        self.assertGreater(solver.n_bisections, 0)


if __name__ == "__main__":
    unittest.main()
