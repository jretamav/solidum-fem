"""Primer paso adimensional y llegada exacta a ``max_lambda`` del arc-length.

Cambio de formulación validado por el usuario el 2026-09-23 (deuda #26):

1. **Primer paso como fracción de carga.** ``initial_dlambda`` = Δλ₁ (0.1 por
   omisión) y el solver deriva la longitud de arco con el predictor
   elástico, Δl₁ = Δλ₁·‖K⁻¹·F_ref‖ (Crisfield 1991, cap. 9). Antes el valor
   por omisión era ``initial_dl = 0.1``, una longitud: el mismo número
   significaba cosas distintas según las unidades y la rigidez del modelo.
2. **Llegada exacta dentro del tramo recorrido.** Un paso de arco que
   converge por encima de ``max_lambda`` no se consolida; se repite con
   λ = ``max_lambda`` desde el estado anterior. Antes ese paso terminaba el
   trazado con λ > ``max_lambda`` (medido: arco de von Mises, 0.398 con
   ``max_lambda = 0.39``), y el paso final se decidía extrapolando la
   tangente, en control de carga.
"""
from __future__ import annotations

import math
import os
import tempfile
import unittest

import numpy as np

from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4
from solidum.elements.truss import Truss2DCorot
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.von_mises_2d import VonMises2D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import ArcLengthSolver, DissipationArcLengthSolver
from solidum.utils.yaml_parser import YamlParser


def _placa(material_factory, escala: float = 1.0):
    """Dos Quad4 (2 m × 1 m × 0.01 m en SI) traccionados en x.

    ``escala = 1000`` es el mismo modelo en mm, N y MPa: longitudes × 1000,
    módulos y esfuerzos / 1000², fuerzas iguales. Devuelve el dominio, la
    carga de referencia y el DOF ``ux`` de un nodo cargado.
    """
    dom = Domain()
    pts = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]
    nodos = [dom.add_node(i + 1, [x * escala, y * escala]) for i, (x, y) in enumerate(pts)]
    mat = material_factory(escala)
    t = 0.01 * escala
    dom.add_element(Quad4(1, [nodos[0], nodos[1], nodos[4], nodos[3]], mat, thickness=t))
    dom.add_element(Quad4(2, [nodos[1], nodos[2], nodos[5], nodos[4]], mat, thickness=t))
    for n in (nodos[0], nodos[3]):
        n.fix_dof("ux", 0.0)
    nodos[0].fix_dof("uy", 0.0)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    for n in (nodos[2], nodos[5]):
        F[n.dofs["ux"]] = 1.5e6                       # 3 MN en total
    return dom, F, nodos[2].dofs["ux"]


def _j2(escala):
    s2 = escala ** 2
    return VonMises2D(E=200e9 / s2, nu=0.3, sigma_y=250e6 / s2, H=2e9 / s2,
                      hypothesis="plane_stress")


def _elastico(escala):
    return Elastic2D(E=200e9 / escala ** 2, nu=0.3, hypothesis="plane_stress")


def _trazar(dom, F, dof, **kwargs):
    solver = ArcLengthSolver(Assembler(dom), max_steps=kwargs.pop("max_steps", 200), **kwargs)
    hist = []
    solver.solve(F, step_callback=lambda k, U, lam: hist.append((lam, U[dof])))
    lam, u = (np.array(v) for v in zip(*hist))
    return solver, lam, u


def _arco_von_mises():
    """Arco de dos barras de von Mises (h/L = 0.1): pico en λ ≈ 0.385 con
    la carga de referencia de ``test_snap_through_corot.py``."""
    L, h, E, A = 10.0, 1.0, 1.0e5, 1.0
    dom = Domain()
    n1, n2, n3 = (dom.add_node(1, [0.0, 0.0]), dom.add_node(2, [L, h]),
                  dom.add_node(3, [2 * L, 0.0]))
    mat = Elastic1D(E=E)
    dom.add_element(Truss2DCorot(1, [n1, n2], mat, A=A))
    dom.add_element(Truss2DCorot(2, [n2, n3], mat, A=A))
    for n in (n1, n3):
        n.fix_dof("ux", 0.0)
        n.fix_dof("uy", 0.0)
    n2.fix_dof("ux", 0.0)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["uy"]] = -10.0 * E * A * (h / L) ** 2
    return dom, F


class TestPrimerPasoAdimensional(unittest.TestCase):

    def test_por_omision_el_primer_paso_lleva_el_10_por_ciento(self):
        dom, F, dof = _placa(_elastico)
        _, lam, _ = _trazar(dom, F, dof, max_lambda=1.0)
        self.assertAlmostEqual(lam[0], 0.1, places=12)

    def test_primer_paso_es_la_fraccion_declarada(self):
        dom, F, dof = _placa(_elastico)
        _, lam, _ = _trazar(dom, F, dof, max_lambda=1.0, initial_dlambda=0.25)
        self.assertAlmostEqual(lam[0], 0.25, places=12)

    def test_trazado_invariante_ante_las_unidades(self):
        """El mismo modelo J2 en m y en mm da la misma secuencia de λ con el
        primer paso por omisión; con una longitud fija, no."""
        dom_m, F_m, dof_m = _placa(_j2, 1.0)
        dom_mm, F_mm, dof_mm = _placa(_j2, 1000.0)
        s_m, lam_m, u_m = _trazar(dom_m, F_m, dof_m, max_lambda=1.0)
        s_mm, lam_mm, u_mm = _trazar(dom_mm, F_mm, dof_mm, max_lambda=1.0)
        self.assertEqual(len(lam_m), len(lam_mm))
        np.testing.assert_allclose(lam_mm, lam_m, rtol=1e-9)
        np.testing.assert_allclose(u_mm, 1000.0 * u_m, rtol=1e-8)
        self.assertAlmostEqual(s_mm.dl_reference / s_m.dl_reference, 1000.0, places=6)
        # Contraste: la misma longitud explícita no significa lo mismo en m y
        # en mm (modelos nuevos: el estado plástico no se reutiliza).
        dom_m, F_m, dof_m = _placa(_j2, 1.0)
        dom_mm, F_mm, dof_mm = _placa(_j2, 1000.0)
        _, lam_m_dl, _ = _trazar(dom_m, F_m, dof_m, max_lambda=1.0, initial_dl=1e-4, max_steps=1)
        _, lam_mm_dl, _ = _trazar(dom_mm, F_mm, dof_mm, max_lambda=1.0, initial_dl=1e-4, max_steps=1)
        self.assertGreater(lam_m_dl[0] / lam_mm_dl[0], 100.0)

    def test_initial_dl_explicito_se_respeta(self):
        dom, F, dof = _placa(_elastico)
        solver, _, _ = _trazar(dom, F, dof, max_lambda=1.0, initial_dl=2e-4)
        self.assertEqual(solver.dl_reference, 2e-4)

    def test_declaracion_invalida_del_primer_paso(self):
        dom, _, _ = _placa(_elastico)
        asm = Assembler(dom)
        with self.assertRaises(ValueError):
            ArcLengthSolver(asm, initial_dl=0.1, initial_dlambda=0.1)
        with self.assertRaises(ValueError):
            ArcLengthSolver(asm, initial_dlambda=0.0)
        with self.assertRaises(ValueError):
            ArcLengthSolver(asm, initial_dl=-1.0)
        with self.assertRaises(ValueError):
            DissipationArcLengthSolver(asm, initial_tau=0.01, initial_dl=0.1,
                                       initial_dlambda=0.1)

    def test_carga_de_referencia_nula_es_error_claro(self):
        dom, F, _ = _placa(_elastico)
        with self.assertRaisesRegex(ValueError, "no produce desplazamiento"):
            ArcLengthSolver(Assembler(dom)).solve(np.zeros_like(F))

    def test_initial_dlambda_desde_yaml(self):
        yaml_txt = (
            "nodes:\n  - {id: 1, coords: [0.0, 0.0]}\n  - {id: 2, coords: [1.0, 0.0]}\n"
            "materials:\n  - {id: 1, type: Elastic1D, E: 1000.0}\n"
            "elements:\n  - {id: 1, type: Truss2D, material: 1, A: 1.0, nodes: [1, 2]}\n"
            "boundary_conditions_by_node:\n  - {node_id: 1, ux: 0.0, uy: 0.0}\n"
            "  - {node_id: 2, uy: 0.0}\n"
            "point_loads_by_node:\n  - {node_id: 2, ux: 10.0}\n"
            "solver: {type: ArcLengthSolver, initial_dlambda: 0.2}\n"
        )
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "m.yaml")
            with open(path, "w", encoding="utf-8") as f:
                f.write(yaml_txt)
            parser = YamlParser(path)
            dom = parser.parse()
            dom.generate_equation_numbers()
            solver = parser.get_solver(Assembler(dom))
        self.assertEqual(solver.initial_dlambda, 0.2)
        self.assertIsNone(solver.initial_dl)


class TestLlegadaExactaAMaxLambda(unittest.TestCase):

    def _arco(self, max_lambda, initial_dl):
        dom, F = _arco_von_mises()
        asm = Assembler(dom)
        solver = ArcLengthSolver(
            asm, convergence=ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8),
            max_lambda=max_lambda, initial_dl=initial_dl, max_steps=200, max_iter=30,
        )
        lams = []
        U = solver.solve(F, step_callback=lambda k, U, lam: lams.append(lam))
        return dom, asm, F, solver, U, lams

    def test_paso_que_cruza_max_lambda_aterriza_exacto(self):
        """Antes: λ_final = 0.398 con max_lambda = 0.39 (el paso de arco que
        cruzaba terminaba el trazado)."""
        for initial_dl in (0.2, 1.0):
            with self.subTest(initial_dl=initial_dl):
                dom, asm, F, solver, U, lams = self._arco(0.39, initial_dl)
                self.assertTrue(solver.reached_max_lambda)
                self.assertAlmostEqual(solver.lambda_final, 0.39, places=12)
                self.assertAlmostEqual(lams[-1], 0.39, places=12)
                # Equilibrio en el estado final con λ = max_lambda.
                _, F_int = asm.assemble_non_linear_system(U)
                free = asm.constraint_set.free_dofs(dom.total_dofs)
                R = 0.39 * F - F_int
                self.assertLess(np.linalg.norm(R[free]), 1e-6 * np.linalg.norm(F))

    def test_elastico_cierra_en_target(self):
        dom, F, dof = _placa(_elastico)
        solver, lam, _ = _trazar(dom, F, dof, max_lambda=0.7)
        self.assertAlmostEqual(lam[-1], 0.7, places=12)
        self.assertTrue(np.all(lam <= 0.7 + 1e-12))
        self.assertTrue(solver.reached_max_lambda)


if __name__ == "__main__":
    unittest.main()
