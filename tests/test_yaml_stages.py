"""Etapas de carga desde el YAML (bloque ``stages``; ADR 0026 y 0029, deuda #36).

Cada etapa declara sus cargas (``loads``, los mismos bloques del nivel raíz)
y su ``solver``, y continúa la anterior. Criterios: el YAML da exactamente
lo mismo que el análisis por etapas hecho desde Python (estático → estático
con plasticidad, y estático → transitorio con las cargas de la etapa como
escalón); ``run_yaml`` devuelve un ``StagedResult``; los errores de forma se
informan al leer.
"""
import logging
import os
import tempfile
import textwrap
import unittest

import numpy as np

import solidum
from solidum import Assembler, Domain, StagedResult
from solidum.elements.frame import Frame2DEuler
from solidum.entry import run_yaml
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.solvers import IndirectDisplacementSolver, NewmarkSolver, NonlinearSolver
from solidum.sections import FiberSection, RectangularPatch
from solidum.utils.yaml_parser import YamlValidationError

E, SY, HH, RHO = 2.0e11, 2.5e8, 2.0e9, 7850.0
B, H, L, N = 0.2, 0.3, 3.0, 3

_MODEL = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [0.0, 1.0]}
  - {id: 3, coords: [0.0, 2.0]}
  - {id: 4, coords: [0.0, 3.0]}
materials:
  - {id: 1, type: Elastoplastic1D, E: 2.0e11, sigma_y: 2.5e8, H: 2.0e9, density: 7850.0}
sections:
  - id: 1
    type: FiberSection
    patches:
      - {shape: rectangle, material: 1, y: [-0.15, 0.15], z: [-0.1, 0.1], ny: 16, nz: 1}
elements:
  - {id: 1, type: Frame2DEuler, nodes: [1, 2], section: 1}
  - {id: 2, type: Frame2DEuler, nodes: [2, 3], section: 1}
  - {id: 3, type: Frame2DEuler, nodes: [3, 4], section: 1}
boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
"""

_STATIC_STAGES = _MODEL + """\
stages:
  - name: gravedad
    loads:
      point_loads:
        - {node_id: 4, uy: -1.0e6}
    solver: {type: NonlinearSolver, num_steps: 4}
  - name: empuje
    loads:
      point_loads:
        - {node_id: 4, ux: 1.0e5}
    solver:
      type: IndirectDisplacementSolver
      control:
        - {node: 4, dof: ux, coef: 1.0}
      initial_dl: 0.004
      control_target: 0.06
      max_lambda: 1.0e9
      max_steps: 500
"""

_DYNAMIC_STAGES = _MODEL + """\
stages:
  - name: gravedad
    loads:
      point_loads:
        - {node_id: 4, uy: -1.0e5}
    solver: {type: NonlinearSolver, num_steps: 2}
  - name: escalon
    loads:
      point_loads:
        - {node_id: 4, ux: 5.0e4}
    solver: {type: NewmarkSolver, t_end: 0.1, dt: 0.002}
"""


def _write(tmp, text):
    path = os.path.join(tmp, "modelo.yaml")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(textwrap.dedent(text))
    return path


def _python_model():
    steel = Elastoplastic1D(E=E, sigma_y=SY, H=HH, density=RHO)
    section = FiberSection(patches=[RectangularPatch(steel, y=(-H / 2, H / 2),
                                                     z=(-B / 2, B / 2), ny=16, nz=1)])
    dom = Domain()
    nodes = [dom.add_node(i + 1, [0.0, float(i)]) for i in range(N + 1)]
    for k in range(N):
        dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], section=section))
    for dof in ("ux", "uy", "rz"):
        nodes[0].fix_dof(dof, 0.0)
    dom.generate_equation_numbers()
    return dom, nodes[-1]


def _load(dom, node, dof, value):
    F = np.zeros(dom.total_dofs)
    F[node.dofs[dof]] = value
    return F


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


class TestEtapasEstaticas(_Quiet):

    def test_igual_que_python(self):
        with tempfile.TemporaryDirectory() as tmp:
            res = run_yaml(_write(tmp, _STATIC_STAGES))
        self.assertIsInstance(res, StagedResult)
        self.assertEqual(res.names, ("gravedad", "empuje"))
        self.assertIs(res["empuje"], res.final)
        self.assertEqual(len(res), 2)

        dom, top = _python_model()
        asm = Assembler(dom)
        r1 = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=4),
                         F_applied=_load(dom, top, "uy", -1.0e6))
        s2 = IndirectDisplacementSolver(asm, control=[(top.id, "ux", 1.0)], initial_dl=0.004,
                                        control_target=0.06, max_lambda=1e9, max_steps=500)
        r2 = solidum.run(dom, assembler=asm, solver=s2,
                         F_applied=_load(dom, top, "ux", 1.0e5), previous=r1)
        np.testing.assert_allclose(res["gravedad"].U, r1.U, rtol=0, atol=1e-14)
        np.testing.assert_allclose(res.final.U, r2.U, rtol=0, atol=1e-12)
        self.assertAlmostEqual(res.final.load_factor, r2.load_factor, delta=1e-12)
        # La carga axial sigue completa en la segunda etapa, que plastifica.
        self.assertAlmostEqual(res.final.reactions_by_node[1]["uy"], 1.0e6, delta=1e-3)
        self.assertEqual(res.final.termination, "control_target")


class TestEtapaEnElTiempo(_Quiet):

    def test_estatico_y_escalon_dinamico(self):
        with tempfile.TemporaryDirectory() as tmp:
            res = run_yaml(_write(tmp, _DYNAMIC_STAGES))
        dom, top = _python_model()
        asm = Assembler(dom)
        r0 = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=2),
                         F_applied=_load(dom, top, "uy", -1.0e5))
        F_step = _load(dom, top, "ux", 5.0e4)
        r1 = solidum.run_transient(dom, assembler=asm, previous=r0, solver=NewmarkSolver(
            asm, t_end=0.1, dt=0.002, F_func=lambda t: F_step))
        dyn = res["escalon"]
        np.testing.assert_allclose(dyn.u_history, r1.u_history, rtol=0, atol=1e-14)
        np.testing.assert_array_equal(dyn.u_history[:, 0], res["gravedad"].U)
        np.testing.assert_array_equal(dyn.F_constant, res["gravedad"].F_applied)
        self.assertGreater(np.abs(dyn.u_history[top.dofs["ux"]]).max(), 1e-4)


class TestErrores(_Quiet):

    def _parse_errors(self, text):
        from solidum.utils.yaml_parser import YamlParser
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(YamlValidationError) as cm:
                YamlParser(_write(tmp, text)).parse()
        return str(cm.exception)

    def test_solver_y_cargas_en_la_raiz(self):
        msg = self._parse_errors(_STATIC_STAGES + "solver: {type: LinearSolver}\n"
                                 "point_loads:\n  - {node_id: 4, ux: 1.0}\n")
        self.assertIn("se declaran en cada etapa", msg)
        self.assertIn("point_loads", msg)

    def test_formas_invalidas(self):
        text = _MODEL + textwrap.dedent("""\
            stages:
              - name: a
                loads: {cargas: 1}
                solver: {type: LinearSolver}
              - name: a
                solver: {type: NoExiste}
              - name: b
                extra: 1
            """)
        msg = self._parse_errors(text)
        self.assertIn("stages[0].loads: bloques desconocidos ['cargas']", msg)
        self.assertIn("nombre de etapa repetido", msg)
        self.assertIn("stages[1].solver: tipo desconocido 'NoExiste'", msg)
        self.assertIn("stages[2]: claves desconocidas ['extra']", msg)
        self.assertIn("stages[2]: falta el bloque 'solver'", msg)

    def test_modal_en_una_etapa(self):
        text = _MODEL + textwrap.dedent("""\
            stages:
              - name: gravedad
                loads:
                  point_loads:
                    - {node_id: 4, uy: -1.0e5}
                solver: {type: LinearSolver}
              - name: modos
                solver: {type: ModalSolver, n_modes: 2}
            """)
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "decisión física"):
                run_yaml(_write(tmp, text))

    def test_estatico_tras_transitorio(self):
        text = _MODEL + textwrap.dedent("""\
            stages:
              - name: dinamico
                solver: {type: NewmarkSolver, t_end: 0.01, dt: 0.002}
              - name: estatico
                loads:
                  point_loads:
                    - {node_id: 4, ux: 1.0e3}
                solver: {type: LinearSolver}
            """)
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "no puede continuar un análisis en el tiempo"):
                run_yaml(_write(tmp, text))


if __name__ == "__main__":
    unittest.main()
