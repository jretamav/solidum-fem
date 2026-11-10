"""Paso adaptativo de Newton-Newmark y Newton-HHT (decisión del usuario del
2026-11-09): el paso cuyo Newton no converge se repite en dos mitades, hasta
``max_subdivisions`` niveles; el historial queda en la malla del usuario.

Criterios (columna con plasticidad en fibras y un pulso lateral):

1. Un análisis que sin subdivisión se detiene (350 kN, Δt = 2 ms, en
   t = 0.082 s) termina; el historial tiene los instantes de la malla pedida,
   los pasos anteriores al primer subdividido son idénticos a los de la
   corrida sin subdivisión, y el máximo difiere < 0.5 % del de una corrida
   con Δt = 1 ms uniforme.
2. Sin pasos difíciles no se subdivide nada y el resultado no cambia.
3. ``max_subdivisions = 0`` se detiene como antes.
4. Si ni subdividiendo converge, el resultado parcial termina en el último
   subpaso convergido (fuera de la malla), que es el estado de los elementos,
   y un análisis que lo continúa (ADR 0029) arranca de él.
5. Newton-HHT igual; el parámetro se valida y se admite desde el YAML.
"""
import logging
import textwrap
import unittest

import numpy as np

import solidum
from solidum import Assembler
from solidum.math.solvers import NewtonHHTSolver, NewtonNewmarkSolver
from solidum.math.solvers.diagnostics import SolverDivergedError
from test_transient_continuation import column, pulse, tight

DT, T_END = 2e-3, 0.2


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


def run(F0, cls=NewtonNewmarkSolver, dt=DT, **kw):
    dom, top = column(plastic=True)
    asm = Assembler(dom)
    solver = cls(asm, t_end=T_END, dt=dt, F_func=pulse(dom, top, F0), convergence=tight(), **kw)
    return dom, top, asm, solver


class TestPasoDificil(_Quiet):
    """Criterios 1 y 3."""

    def test_termina_en_la_malla_pedida(self):
        _, top, _, s0 = run(3.5e5, max_subdivisions=0)
        with self.assertRaises(SolverDivergedError) as ctx:
            s0.solve()
        partial = ctx.exception.partial_result
        self.assertAlmostEqual(s0.t_final, 0.082, delta=1e-12)

        _, top, _, s = run(3.5e5)
        res = s.solve()
        self.assertEqual(res.termination, "t_end")
        self.assertGreater(s.n_subdivisions, 0)
        np.testing.assert_allclose(res.t_history, np.linspace(0.0, T_END, 101), rtol=0, atol=1e-15)
        n = partial.t_history.size
        np.testing.assert_array_equal(res.u_history[:, :n], partial.u_history)

        _, top1, _, s1 = run(3.5e5, dt=1e-3)
        ref = s1.solve()
        i = top.dofs["ux"]
        self.assertAlmostEqual(np.max(np.abs(res.u_history[i])) / np.max(np.abs(ref.u_history[i])),
                               1.0, delta=5e-3)


class TestSinPasosDificiles(_Quiet):
    """Criterio 2."""

    def test_no_subdivide_ni_cambia(self):
        _, _, _, a = run(1.0e5)
        _, _, _, b = run(1.0e5, max_subdivisions=0)
        ra, rb = a.solve(), b.solve()
        self.assertEqual(a.n_subdivisions, 0)
        np.testing.assert_array_equal(ra.u_history, rb.u_history)


class TestFracasoFinal(_Quiet):
    """Criterio 4."""

    def test_parcial_termina_en_el_subpaso(self):
        dom, top, asm, s = run(1.5e6, max_subdivisions=2)
        with self.assertRaises(SolverDivergedError) as ctx:
            s.solve()
        partial = ctx.exception.partial_result
        t_last = partial.t_history[-1]
        # Fuera de la malla de 2 ms: el último subpaso convergido.
        self.assertGreater(abs(t_last / DT - round(t_last / DT)), 1e-6)
        self.assertAlmostEqual(s.t_final, t_last, delta=1e-15)
        np.testing.assert_array_equal(s.U_final, partial.u_history[:, -1])
        self.assertEqual(partial.R_history.shape[1], partial.t_history.size)
        # Continuarlo con un paso menor arranca de ese estado.
        follow = NewtonNewmarkSolver(asm, t_end=t_last + 4 * 2.5e-4, dt=2.5e-4,
                                     F_func=pulse(dom, top, 1.5e6), convergence=tight())
        res = solidum.run_transient(dom, assembler=asm, solver=follow, previous=partial)
        self.assertAlmostEqual(res.t_history[0], t_last, delta=1e-15)
        np.testing.assert_array_equal(res.u_history[:, 0], partial.u_history[:, -1])


class TestHHTYParametro(_Quiet):
    """Criterio 5."""

    def test_hht(self):
        _, _, _, s0 = run(3.5e5, cls=NewtonHHTSolver, max_subdivisions=0)
        try:
            s0.solve()
            failed = False
        except SolverDivergedError:
            failed = True
        _, _, _, s = run(3.5e5, cls=NewtonHHTSolver)
        res = s.solve()
        self.assertEqual(res.termination, "t_end")
        if failed:
            self.assertGreater(s.n_subdivisions, 0)

    def test_parametro(self):
        dom, top = column(plastic=True)
        with self.assertRaisesRegex(ValueError, "max_subdivisions"):
            NewtonNewmarkSolver(Assembler(dom), t_end=0.1, dt=1e-3, max_subdivisions=-1)
        import os
        import tempfile
        text = textwrap.dedent("""\
            nodes:
              - {id: 1, coords: [0.0, 0.0]}
              - {id: 2, coords: [0.0, 3.0]}
            materials:
              - {id: 1, type: Elastic1D, E: 2.0e11, density: 7850.0}
            elements:
              - {id: 1, type: Frame2DEuler, nodes: [1, 2], material: 1, A: 0.06, I: 4.5e-4}
            boundary_conditions_by_node:
              - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
            point_loads:
              - {node_id: 2, ux: 1.0e3}
            solver: {type: NewtonHHTSolver, t_end: 0.01, dt: 0.002, max_subdivisions: 3}
            """)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "m.yaml")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text)
            model = solidum.load_yaml(path)
            solver = model.parser.get_solver(model.assembler)
        self.assertEqual(solver.max_subdivisions, 3)


if __name__ == "__main__":
    unittest.main()
