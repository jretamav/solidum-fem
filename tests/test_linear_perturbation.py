"""``LinearSolver`` sobre una etapa: perturbación lineal (ADR 0031, decisión
del usuario del 2026-11-09). ``K_t(u_b)·ΔU = ΔF``, ``U = u_b + ΔU``; las
reacciones y las fuerzas internas son lineales en ``ΔU``; el estado no cambia.

Criterios:

1. Modelo elástico lineal: la perturbación da lo mismo que un análisis lineal
   con todas las cargas (desplazamientos, reacciones, fuerzas internas).
2. P-Δ: columna esbelta comprimida por P = 0.5·P_cr; una carga lateral H da
   en la cabeza ``H/k(P)`` con la rigidez lateral exacta, a 1e-3.
3. Material no lineal (fibras que fluyeron): dos casos laterales son
   superponibles en desplazamientos, reacciones y fuerzas internas; el
   equilibrio global se cumple con las reacciones lineales; una etapa
   posterior da lo mismo con o sin las perturbaciones.
4. Una perturbación no se puede continuar; en el YAML la etapa siguiente
   continúa la anterior.
"""
import logging
import math
import os
import tempfile
import textwrap
import unittest

import numpy as np

import solidum
from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler
from solidum.materials.elastic import Elastic1D
from solidum.math.solvers import IndirectDisplacementSolver, LinearSolver, NonlinearSolver
from test_linearized_stage import I_S, P_CR, axial_stage, k_lateral, steel_cantilever

E, L = 2.0e11, 3.0


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


def frame():
    """Marco elástico de un nivel con una diagonal."""
    dom = Domain()
    mat = Elastic1D(E=E)
    pts = {1: (0.0, 0.0), 2: (0.0, 3.0), 3: (5.0, 3.0), 4: (5.0, 0.0)}
    n = {k: dom.add_node(k, list(v)) for k, v in pts.items()}
    for eid, (a, b) in enumerate([(1, 2), (2, 3), (3, 4), (1, 3)], start=1):
        dom.add_element(Frame2DEuler(eid, [n[a], n[b]], material=mat, A=0.02, I=3e-4))
    for nid in (1, 4):
        for d in ("ux", "uy", "rz"):
            n[nid].fix_dof(d, 0.0)
    dom.generate_equation_numbers()
    return dom, n


def load(dom, items):
    F = np.zeros(dom.total_dofs)
    for node, dof, value in items:
        F[node.dofs[dof]] += value
    return F


class TestModeloLineal(_Quiet):
    """Criterio 1."""

    def test_igual_que_un_lineal_con_todo(self):
        dom, n = frame()
        Fg = load(dom, [(n[2], "uy", -2e5), (n[3], "uy", -2e5)])
        Fh = load(dom, [(n[2], "ux", 5e4)])
        asm = Assembler(dom)
        base = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=2),
                           F_applied=Fg)
        pert = solidum.run(dom, assembler=asm, solver=LinearSolver(asm), F_applied=Fh,
                           previous=base)
        self.assertTrue(pert.perturbation)
        dom2, n2 = frame()
        ref = solidum.run(dom2, F_applied=load(dom2, [(n2[2], "uy", -2e5), (n2[3], "uy", -2e5),
                                                     (n2[2], "ux", 5e4)]))
        scale = np.max(np.abs(ref.U))
        np.testing.assert_allclose(pert.U, ref.U, rtol=0, atol=1e-10 * scale)
        np.testing.assert_allclose(pert.R, ref.R, rtol=0, atol=1e-8 * np.max(np.abs(ref.R)))
        for eid, ef in ref.element_forces.items():
            for comp, val in ef.components.items():
                np.testing.assert_allclose(pert.element_forces[eid].components[comp], val,
                                           rtol=0, atol=1e-6 * max(1.0, np.max(np.abs(val))))


class TestPDelta(_Quiet):
    """Criterio 2."""

    def test_carga_lateral_sobre_columna_comprimida(self):
        dom, top, asm, base = axial_stage(0.5 * P_CR)
        H = 1.0e3
        F = np.zeros(dom.total_dofs)
        F[top.dofs["ux"]] = H
        pert = solidum.run(dom, assembler=asm, solver=LinearSolver(asm), F_applied=F,
                           previous=base)
        delta = pert.U[top.dofs["ux"]] - base.U[top.dofs["ux"]]
        self.assertAlmostEqual(delta / (H / k_lateral(0.5 * P_CR)), 1.0, delta=1e-3)
        # P-Δ: la compresión amplifica el desplazamiento de la columna sin
        # carga axial (≈ 1/(1 − P/P_cr) = 2).
        amplification = delta / (H * L ** 3 / (3 * E * I_S))
        self.assertGreater(amplification, 1.9)
        self.assertLess(amplification, 2.1)
        # Equilibrio global con las reacciones lineales.
        self.assertAlmostEqual(pert.reactions_by_node[1]["ux"] / -H, 1.0, delta=1e-9)
        self.assertAlmostEqual(pert.reactions_by_node[1]["uy"] / (0.5 * P_CR), 1.0, delta=1e-9)


class TestMaterialNoLineal(_Quiet):
    """Criterio 3."""

    def _base(self):
        dom, top = steel_cantilever()
        asm = Assembler(dom)
        F = np.zeros(dom.total_dofs)
        F[top.dofs["ux"]] = 1.0e5
        s = IndirectDisplacementSolver(asm, control=[(top.id, "ux", 1.0)], control_target=0.04,
                                       initial_dl=0.004, max_lambda=1e9, max_steps=200)
        base = solidum.run(dom, assembler=asm, solver=s, F_applied=F)
        return dom, top, asm, base

    def test_superposicion_y_estado_intacto(self):
        dom, top, asm, base = self._base()
        H1 = np.zeros(dom.total_dofs)
        H1[top.dofs["ux"]] = 2.0e4
        H2 = np.zeros(dom.total_dofs)
        H2[top.dofs["uy"]] = -3.0e4
        p1 = solidum.run(dom, assembler=asm, solver=LinearSolver(asm), F_applied=H1, previous=base)
        p2 = solidum.run(dom, assembler=asm, solver=LinearSolver(asm), F_applied=H2, previous=base)
        p12 = solidum.run(dom, assembler=asm, solver=LinearSolver(asm), F_applied=H1 + H2,
                          previous=base)
        scale = np.max(np.abs(p12.U - base.U))
        np.testing.assert_allclose(p12.U - base.U, (p1.U - base.U) + (p2.U - base.U),
                                   rtol=0, atol=1e-10 * scale)
        np.testing.assert_allclose(p12.R - base.R, (p1.R - base.R) + (p2.R - base.R),
                                   rtol=0, atol=1e-6)
        for eid in base.element_forces:
            for comp in ("N", "V", "M"):
                b = base.element_forces[eid].components[comp]
                np.testing.assert_allclose(
                    p12.element_forces[eid].components[comp] - b,
                    (p1.element_forces[eid].components[comp] - b)
                    + (p2.element_forces[eid].components[comp] - b), rtol=0, atol=1e-5)
        # Equilibrio global en x con las reacciones lineales.
        self.assertAlmostEqual(p12.reactions_by_node[1]["ux"],
                               -(base.F_applied[top.dofs["ux"]] + 2.0e4), delta=1e-4)
        # El estado no cambió: seguir empujando da lo mismo que sin ellas.
        dom2, top2, asm2, base2 = self._base()
        F = np.zeros(dom.total_dofs)
        F[top.dofs["ux"]] = 1.0e5
        push = lambda d, t, a, prev: solidum.run(d, assembler=a, solver=IndirectDisplacementSolver(
            a, control=[(t.id, "ux", 1.0)], control_target=0.06, initial_dl=0.004,
            max_lambda=1e9, max_steps=200), F_applied=F, previous=prev)
        np.testing.assert_allclose(push(dom, top, asm, base).U, push(dom2, top2, asm2, base2).U,
                                   rtol=0, atol=1e-12)


class TestNoSeContinua(_Quiet):
    """Criterio 4."""

    def test_guardas(self):
        dom, n = frame()
        asm = Assembler(dom)
        base = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=1),
                           F_applied=load(dom, [(n[2], "uy", -1e5)]))
        pert = solidum.run(dom, assembler=asm, solver=LinearSolver(asm),
                           F_applied=load(dom, [(n[2], "ux", 1e4)]), previous=base)
        for call in (
            lambda: solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=1),
                                F_applied=load(dom, [(n[2], "ux", 1e4)]), previous=pert),
            lambda: solidum.run_modal(dom, assembler=asm, n_modes=1, previous=pert),
        ):
            with self.assertRaisesRegex(ValueError, "perturbación lineal"):
                call()

    def test_yaml(self):
        model = textwrap.dedent("""\
            nodes:
              - {id: 1, coords: [0.0, 0.0]}
              - {id: 2, coords: [0.0, 1.5]}
              - {id: 3, coords: [0.0, 3.0]}
            materials:
              - {id: 1, type: Elastoplastic1D, E: 2.0e11, sigma_y: 2.5e8, H: 2.0e9}
            sections:
              - id: 1
                type: FiberSection
                patches:
                  - {shape: rectangle, material: 1, y: [-0.15, 0.15], z: [-0.1, 0.1], ny: 16, nz: 1}
            elements:
              - {id: 1, type: Frame2DEuler, nodes: [1, 2], section: 1}
              - {id: 2, type: Frame2DEuler, nodes: [2, 3], section: 1}
            boundary_conditions_by_node:
              - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
            """)
        gravity = """\
              - name: gravedad
                loads:
                  point_loads:
                    - {node_id: 3, uy: -5.0e5}
                solver: {type: NonlinearSolver, num_steps: 2}
            """
        linear = """\
              - name: caso_viento
                loads:
                  point_loads:
                    - {node_id: 3, ux: 2.0e4}
                solver: {type: LinearSolver}
            """
        push = """\
              - name: empuje
                loads:
                  point_loads:
                    - {node_id: 3, ux: 1.0e5}
                solver: {type: NonlinearSolver, num_steps: 4}
            """
        with tempfile.TemporaryDirectory() as tmp:
            def run(stages):
                path = os.path.join(tmp, "m.yaml")
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(model + "stages:\n" + "".join(
                        textwrap.indent(textwrap.dedent(st), "  ") for st in stages))
                return solidum.run_yaml(path)
            with_lin = run([gravity, linear, push])
            without = run([gravity, push])
        self.assertTrue(with_lin["caso_viento"].perturbation)
        np.testing.assert_allclose(with_lin["empuje"].U, without["empuje"].U, rtol=0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
