"""Progreso de los análisis en el tiempo (``step_callback``) y nivel de
registro desde el YAML (``output.log_level``).

Criterios:

1. Los cinco integradores mecánicos y el θ-método llaman
   ``step_callback(step, U, t)`` una vez por paso convergido, con el campo y
   el instante que quedan en el historial; nunca en el estado inicial.
2. ``run_transient`` y ``run_thermal_transient`` lo pasan durante la llamada
   y restituyen el del solver; ``run_yaml`` lo lleva a los transitorios y a
   las etapas en el tiempo, y avisa si el análisis no tiene pasos (modal).
3. ``output.log_level`` cambia el nivel durante el análisis y lo restituye
   después; un nivel desconocido es un error de validación.
"""
import logging
import os
import tempfile
import textwrap
import unittest

import numpy as np

import solidum
from solidum import Assembler, load_yaml, run_yaml
from solidum.math.solvers import (
    CentralDifferenceSolver,
    HHTSolver,
    NewmarkSolver,
    NewtonHHTSolver,
    NewtonNewmarkSolver,
)
from solidum.math.solvers.theta_method import ThetaMethodSolver
from solidum.utils.yaml_parser import YamlValidationError
from test_theta_method import barra_2d
from test_transient_continuation import column, pulse, tight

_MECHANICAL = {
    "Newmark": lambda asm, F: NewmarkSolver(asm, t_end=0.02, dt=2e-3, F_func=F),
    "HHT": lambda asm, F: HHTSolver(asm, t_end=0.02, dt=2e-3, F_func=F),
    "NewtonNewmark": lambda asm, F: NewtonNewmarkSolver(asm, t_end=0.02, dt=2e-3, F_func=F,
                                                        convergence=tight()),
    "NewtonHHT": lambda asm, F: NewtonHHTSolver(asm, t_end=0.02, dt=2e-3, F_func=F,
                                                convergence=tight()),
    "CentralDifference": lambda asm, F: CentralDifferenceSolver(asm, t_end=2e-4, dt=2e-5,
                                                                F_func=F),
}

_YAML_MODEL = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [0.0, 1.5]}
  - {id: 3, coords: [0.0, 3.0]}
materials:
  - {id: 1, type: Elastic1D, E: 2.0e11, density: 7850.0}
elements:
  - {id: 1, type: Frame2DEuler, nodes: [1, 2], material: 1, A: 0.06, I: 4.5e-4}
  - {id: 2, type: Frame2DEuler, nodes: [2, 3], material: 1, A: 0.06, I: 4.5e-4}
boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
"""


class _Recorder:
    def __init__(self):
        self.calls = []

    def __call__(self, step, U, t):
        self.calls.append((step, np.array(U, copy=True), t))


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.INFO)
        self.addCleanup(logging.disable, logging.NOTSET)


def _write(tmp, text):
    path = os.path.join(tmp, "modelo.yaml")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(textwrap.dedent(text))
    return path


class TestCallbackPorPaso(_Quiet):
    """Criterios 1 y 2."""

    def test_integradores_mecanicos(self):
        for name, make in _MECHANICAL.items():
            with self.subTest(solver=name):
                dom, top = column()
                asm = Assembler(dom)
                solver = make(asm, pulse(dom, top, 5.0e4, t_p=0.01))
                rec = _Recorder()
                res = solidum.run_transient(dom, assembler=asm, solver=solver,
                                            step_callback=rec)
                self.assertEqual([c[0] for c in rec.calls], list(range(1, res.n_steps + 1)))
                np.testing.assert_array_equal([c[2] for c in rec.calls], res.t_history[1:])
                for step, U, _ in rec.calls:
                    np.testing.assert_array_equal(U, res.u_history[:, step])
                self.assertIsNone(solver.step_callback)

    def test_theta_method(self):
        dom, _ = barra_2d(6, 1.0, T_izq=100.0)
        solver = ThetaMethodSolver(Assembler(dom), dt=5.0, n_steps=4, T_initial=20.0)
        rec = _Recorder()
        res = solidum.run_thermal_transient(dom, solver=solver, step_callback=rec)
        self.assertEqual([c[0] for c in rec.calls], [1, 2, 3, 4])
        np.testing.assert_array_equal([c[2] for c in rec.calls], res.t_history[1:])
        for step, T, _ in rec.calls:
            np.testing.assert_array_equal(T, res.T_history[:, step])

    def test_yaml_transitorio_etapas_y_modal(self):
        transient = _YAML_MODEL + """\
point_loads:
  - {node_id: 3, ux: 1.0e4}
solver: {type: NewmarkSolver, t_end: 0.01, dt: 0.002}
"""
        staged = _YAML_MODEL + """\
stages:
  - name: estatica
    loads:
      point_loads:
        - {node_id: 3, uy: -1.0e4}
    solver: {type: NonlinearSolver, num_steps: 3}
  - name: dinamica
    loads:
      point_loads:
        - {node_id: 3, ux: 1.0e4}
    solver: {type: NewmarkSolver, t_end: 0.01, dt: 0.002}
"""
        modal = _YAML_MODEL + "solver: {type: ModalSolver, n_modes: 2}\n"
        with tempfile.TemporaryDirectory() as tmp:
            rec = _Recorder()
            run_yaml(_write(tmp, transient), step_callback=rec)
            self.assertEqual([c[0] for c in rec.calls], [1, 2, 3, 4, 5])
            self.assertAlmostEqual(rec.calls[-1][2], 0.01, delta=1e-15)

            rec = _Recorder()
            res = run_yaml(_write(tmp, staged), step_callback=rec)
            # Tres pasos estáticos (λ) y cinco en el tiempo (t).
            self.assertEqual([c[0] for c in rec.calls], [1, 2, 3, 1, 2, 3, 4, 5])
            np.testing.assert_array_equal([c[2] for c in rec.calls[3:]],
                                          res["dinamica"].t_history[1:])

            logging.disable(logging.NOTSET)
            with self.assertLogs("solidum", level="WARNING") as logs:
                run_yaml(_write(tmp, modal), step_callback=_Recorder())
            self.assertTrue(any("step_callback ignorado" in m for m in logs.output))


class TestNivelDeRegistro(unittest.TestCase):
    """Criterio 3."""

    def test_nivel_durante_el_analisis_y_restitucion(self):
        root = logging.getLogger("solidum")
        before = root.level
        base = _YAML_MODEL + """\
point_loads:
  - {node_id: 3, ux: 1.0e4}
solver: {type: LinearSolver}
"""

        def levels_seen(text):
            """Niveles de los mensajes emitidos durante el análisis (la
            lectura del archivo, anterior, no cuenta)."""
            seen = []

            class Spy(logging.Handler):
                def emit(self, record):
                    seen.append(record.levelno)

            with tempfile.TemporaryDirectory() as tmp:
                model = load_yaml(_write(tmp, text))
            spy = Spy(level=logging.DEBUG)
            root.addHandler(spy)
            try:
                model.solve()
            finally:
                root.removeHandler(spy)
            return seen

        root.setLevel(logging.INFO)
        try:
            self.assertIn(logging.INFO, levels_seen(base))
            quiet = levels_seen(base + "output: {log_level: WARNING}\n")
            self.assertFalse([lvl for lvl in quiet if lvl < logging.WARNING])
            self.assertEqual(root.level, logging.INFO)
        finally:
            root.setLevel(before)

    def test_nivel_desconocido(self):
        text = _YAML_MODEL + "solver: {type: LinearSolver}\noutput: {log_level: MUCHO}\n"
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(YamlValidationError, "output.log_level: 'MUCHO'"):
                load_yaml(_write(tmp, text))


if __name__ == "__main__":
    unittest.main()
