"""Continuidad de los análisis en el tiempo (ADR 0029).

Un análisis en el tiempo puede continuar otro: un estático (la gravedad,
cuya carga se mantiene) o un transitorio del mismo tipo (desde su estado e
instante finales). Criterios:

- **Estático → dinámico.** Con la gravedad estática primero y un pulso
  lateral después, el transitorio parte en reposo del estado estático, sin
  aceleración espuria, la carga axial se mantiene (el acortamiento no
  cambia en un modelo lineal) y el resultado es el mismo que pasar a mano
  ``u0`` y la gravedad dentro de ``F_func``. Las fuerzas internas del primer
  instante coinciden con las del estático, también con peso propio (cargas
  equivalentes, ADR 0027).
- **Parar y continuar es igual que no parar**, en los cinco integradores
  mecánicos (también con plasticidad en las fibras) y en el θ-método.
- Avisos y errores: estado inicial fuera de equilibrio, ``t_end`` anterior
  al instante inicial, condiciones iniciales repetidas, tipo de análisis.
"""
import dataclasses
import logging
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

import solidum
from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import (
    CentralDifferenceSolver,
    HHTSolver,
    LinearSolver,
    NewmarkSolver,
    NewtonHHTSolver,
    NewtonNewmarkSolver,
)
from solidum.math.solvers.theta_method import ThetaMethodSolver
from solidum.sections import FiberSection, RectangularPatch

from test_theta_method import barra_2d, dofs_fila

E, RHO, L, B, H = 2.0e11, 7850.0, 3.0, 0.2, 0.3
A, I = B * H, B * H ** 3 / 12.0
P_GRAV = 2.0e5


def tight():
    return ConvergenceCriterion(rtol_force=1e-12, rtol_disp=1e-12)


def column(n=6, *, plastic=False):
    """Voladizo vertical de ``n`` elementos: elástico con material, o con
    sección por fibras de acero elastoplástico."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, [0.0, i * L / n]) for i in range(n + 1)]
    if plastic:
        steel = Elastoplastic1D(E=E, sigma_y=2.5e8, H=2.0e9, density=RHO)
        section = FiberSection(patches=[RectangularPatch(
            steel, y=(-H / 2, H / 2), z=(-B / 2, B / 2), ny=16, nz=1)])
        make = lambda k: Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], section=section, n_ip=3)
    else:
        mat = Elastic1D(E=E, density=RHO)
        make = lambda k: Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], mat, A=A, I=I)
    for k in range(n):
        dom.add_element(make(k))
    for dof in ("ux", "uy", "rz"):
        nodes[0].fix_dof(dof, 0.0)
    dom.generate_equation_numbers()
    return dom, nodes[-1]


def pulse(dom, top, F0, t_p=0.05):
    dof = top.dofs["ux"]

    def F_func(t):
        F = np.zeros(dom.total_dofs)
        if 0.0 <= t < t_p:
            F[dof] = F0 * math.sin(math.pi * t / t_p)
        return F
    return F_func


def gravity(dom, top):
    F = np.zeros(dom.total_dofs)
    F[top.dofs["uy"]] = -P_GRAV
    return F


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.INFO)
        self.addCleanup(logging.disable, logging.NOTSET)


class TestEstaticoADinamico(_Quiet):

    def test_parte_en_reposo_y_la_gravedad_se_mantiene(self):
        dom, top = column()
        asm = Assembler(dom)
        r0 = solidum.run(dom, assembler=asm, F_applied=gravity(dom, top))
        solver = NewmarkSolver(asm, t_end=0.2, dt=2e-3, F_func=pulse(dom, top, 5e4))
        res = solidum.run_transient(dom, assembler=asm, solver=solver, previous=r0)
        np.testing.assert_array_equal(res.u_history[:, 0], r0.U)
        np.testing.assert_allclose(res.uddot_history[:, 0], 0.0, atol=1e-6)
        uy = res.u_history[top.dofs["uy"]]
        np.testing.assert_allclose(uy, r0.U[top.dofs["uy"]], rtol=1e-11)
        self.assertGreater(np.abs(res.u_history[top.dofs["ux"]]).max(), 1e-4)
        np.testing.assert_array_equal(res.F_constant, r0.F_applied)
        self.assertEqual(res.t_history[0], 0.0)

        # Igual que pasar a mano u0 y la gravedad dentro de F_func.
        dom2, top2 = column()
        F_g = gravity(dom2, top2)
        p = pulse(dom2, top2, 5e4)
        manual = NewmarkSolver(Assembler(dom2), t_end=0.2, dt=2e-3, u0=r0.U,
                               F_func=lambda t: F_g + p(t)).solve()
        np.testing.assert_allclose(res.u_history, manual.u_history, rtol=0,
                                   atol=1e-12 * np.abs(manual.u_history).max())

    def test_newton_y_fuerzas_internas_con_peso_propio(self):
        """El peso propio (carga distribuida) se mantiene con sus cargas
        equivalentes: la fuerza axial del primer instante es la estática y
        varía a lo largo de la columna."""
        dom, top = column()
        asm = Assembler(dom)
        F_sw = asm.assemble_self_weight(np.array([0.0, -9.81]))
        r0 = solidum.run(dom, assembler=asm, F_applied=F_sw)
        solver = NewtonNewmarkSolver(asm, t_end=0.05, dt=5e-3, convergence=tight(),
                                     F_func=pulse(dom, top, 2e4),
                                     record_internal_forces=True)
        with self.assertNoLogs("solidum", level=logging.WARNING):
            res = solidum.run_transient(dom, assembler=asm, solver=solver, previous=r0)
        for eid, ef in r0.element_forces.items():
            first = res.element_forces_history[eid][0]
            np.testing.assert_allclose(first.components["N"], ef.components["N"],
                                       rtol=1e-9, atol=1e-6)
        N_base = r0.element_forces[1].components["N"][0]
        N_top = r0.element_forces[6].components["N"][1]
        self.assertAlmostEqual(N_base, -RHO * 9.81 * A * L, delta=1e-6 * abs(N_base))
        self.assertAlmostEqual(N_top, 0.0, delta=1e-6 * abs(N_base))

    def test_estado_fuera_de_equilibrio_avisa(self):
        dom, top = column()
        asm = Assembler(dom)
        r0 = solidum.run(dom, assembler=asm, F_applied=gravity(dom, top))
        wrong = dataclasses.replace(r0, F_applied=2.0 * r0.F_applied)
        solver = NewmarkSolver(asm, t_end=0.01, dt=2e-3)
        logging.disable(logging.NOTSET)
        with self.assertLogs("solidum", level=logging.WARNING) as cm:
            solidum.run_transient(dom, assembler=asm, solver=solver, previous=wrong)
        self.assertTrue(any("no está en equilibrio" in m for m in cm.output))


def _split_vs_single(make_solver, plastic, t_mid, t_end, F0):
    """Historias de una corrida de 0 a ``t_end`` y de dos, 0 → ``t_mid`` y
    ``t_mid`` → ``t_end``, en dominios gemelos."""
    dom, top = column(plastic=plastic)
    single = solidum.run_transient(dom, assembler=Assembler(dom),
                                   solver=make_solver(Assembler(dom), t_end, pulse(dom, top, F0)))
    dom2, top2 = column(plastic=plastic)
    asm2 = Assembler(dom2)
    first = solidum.run_transient(dom2, assembler=asm2,
                                  solver=make_solver(asm2, t_mid, pulse(dom2, top2, F0)))
    second = solidum.run_transient(dom2, assembler=asm2,
                                   solver=make_solver(asm2, t_end, pulse(dom2, top2, F0)),
                                   previous=first)
    return single, first, second


class TestPararYContinuar(_Quiet):

    CASES = {
        "Newmark": (lambda asm, t, F: NewmarkSolver(asm, t_end=t, dt=2e-3, F_func=F), False),
        "HHT": (lambda asm, t, F: HHTSolver(asm, t_end=t, dt=2e-3, F_func=F), False),
        "NewtonNewmark": (lambda asm, t, F: NewtonNewmarkSolver(
            asm, t_end=t, dt=2e-3, F_func=F, convergence=tight()), True),
        "NewtonHHT": (lambda asm, t, F: NewtonHHTSolver(
            asm, t_end=t, dt=2e-3, F_func=F, convergence=tight()), True),
        "CentralDifference": (lambda asm, t, F: CentralDifferenceSolver(
            asm, t_end=t, dt=2e-5, F_func=F, nonlinear=True), True),
    }

    def test_igual_que_una_sola_corrida(self):
        for name, (make, plastic) in self.CASES.items():
            t_mid, t_end = (0.004, 0.008) if name == "CentralDifference" else (0.06, 0.12)
            F0 = 2.8e5 if plastic else 5e4
            with self.subTest(integrador=name):
                single, first, second = _split_vs_single(make, plastic, t_mid, t_end, F0)
                self.assertAlmostEqual(second.t_history[0], t_mid, delta=1e-12)
                self.assertAlmostEqual(second.t_history[-1], t_end, delta=1e-12)
                joined = np.hstack([first.u_history, second.u_history[:, 1:]])
                self.assertEqual(joined.shape, single.u_history.shape)
                scale = np.abs(single.u_history).max()
                np.testing.assert_allclose(joined, single.u_history, rtol=0,
                                           atol=1e-8 * scale)
                v = np.hstack([first.udot_history, second.udot_history[:, 1:]])
                np.testing.assert_allclose(v, single.udot_history, rtol=0,
                                           atol=1e-8 * np.abs(single.udot_history).max())

    def test_con_plasticidad_el_pulso_plastifica(self):
        dom, top = column(plastic=True)
        asm = Assembler(dom)
        solidum.run_transient(dom, assembler=asm, solver=NewtonNewmarkSolver(
            asm, t_end=0.06, dt=2e-3, F_func=pulse(dom, top, 2.8e5), convergence=tight()))
        alpha = max(st["alpha"] for el in dom.elements.values()
                    for g in el.state.vars for st in g if st)
        self.assertGreater(alpha, 0.0)


class TestTermico(_Quiet):

    def test_parar_y_continuar(self):
        def model():
            dom, nid = barra_2d(8, L=0.3, T_izq=100.0, T_der=0.0)
            return dom, Assembler(dom)
        dom, asm = model()
        single = solidum.run_thermal_transient(dom, assembler=asm, dt=5.0, n_steps=40,
                                               T_initial=0.0)
        dom2, asm2 = model()
        first = solidum.run_thermal_transient(dom2, assembler=asm2, dt=5.0, n_steps=15,
                                              T_initial=0.0)
        second = solidum.run_thermal_transient(dom2, assembler=asm2, dt=5.0, n_steps=25,
                                               previous=first)
        self.assertAlmostEqual(second.t_history[0], 75.0, delta=1e-12)
        joined = np.hstack([first.T_history, second.T_history[:, 1:]])
        np.testing.assert_allclose(joined, single.T_history, rtol=0, atol=1e-10)

    def test_desde_el_estacionario(self):
        """Del estacionario, en reposo térmico: sin cambio de condiciones, el
        campo no se mueve."""
        dom, nid = barra_2d(8, L=0.3, T_izq=100.0, T_der=0.0)
        asm = Assembler(dom)
        r0 = solidum.run(dom, assembler=asm)
        res = solidum.run_thermal_transient(dom, assembler=asm, dt=5.0, n_steps=10,
                                            previous=r0)
        for k in range(res.T_history.shape[1]):
            np.testing.assert_allclose(res.T_history[:, k], r0.U, rtol=0, atol=1e-9)


class TestErrores(_Quiet):

    def setUp(self):
        super().setUp()
        self.dom, self.top = column()
        self.asm = Assembler(self.dom)
        self.r0 = solidum.run(self.dom, assembler=self.asm, F_applied=gravity(self.dom, self.top))

    def test_t_end_anterior_al_inicio(self):
        first = solidum.run_transient(self.dom, assembler=self.asm,
                                      solver=NewmarkSolver(self.asm, t_end=0.02, dt=2e-3),
                                      previous=self.r0)
        with self.assertRaisesRegex(ValueError, "instante final absoluto"):
            solidum.run_transient(self.dom, assembler=self.asm,
                                  solver=NewmarkSolver(self.asm, t_end=0.01, dt=2e-3),
                                  previous=first)

    def test_condiciones_iniciales_repetidas(self):
        solver = NewmarkSolver(self.asm, t_end=0.02, dt=2e-3, u0=np.zeros(self.dom.total_dofs))
        with self.assertRaisesRegex(ValueError, "se excluyen"):
            solidum.run_transient(self.dom, assembler=self.asm, solver=solver,
                                  previous=self.r0)

    def test_termico_sin_temperatura_inicial(self):
        dom, nid = barra_2d(4, L=0.3, T_izq=100.0, T_der=0.0)
        with self.assertRaisesRegex(ValueError, "falta la temperatura inicial"):
            ThetaMethodSolver(Assembler(dom), dt=5.0, n_steps=3).solve()

    def test_tipo_de_analisis(self):
        dom, nid = barra_2d(4, L=0.3, T_izq=100.0, T_der=0.0)
        asm = Assembler(dom)
        first = solidum.run_thermal_transient(dom, assembler=asm, dt=5.0, n_steps=3,
                                              T_initial=0.0)
        with self.assertRaisesRegex(TypeError, "mismo tipo"):
            NewmarkSolver(self.asm, t_end=0.02, dt=2e-3).continue_from(first)


if __name__ == "__main__":
    unittest.main()
