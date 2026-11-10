"""Modal, armónico y espectral sobre el estado de una etapa (ADR 0031,
decisión del usuario del 2026-11-09): rigidez tangente del estado.

Criterios:

1. Columna sin masa en voladizo (Frame2DEulerCorot, 20 elementos) con una
   masa en la cabeza y carga axial de compresión P aplicada en una etapa:
   la frecuencia es la del oscilador con la rigidez lateral exacta de la
   columna comprimida, ``k(P) = P·κ/(tan κL − κL)``, ``κ = √(P/EI)``, a 1e-3
   (P = 0.5·P_cr; columna esbelta, P/EA ≈ 1e-4: la fórmula supone la columna
   inextensible). Medido: el error baja de forma cuadrática con la malla, de
   1.2e-3 con 10 elementos a 2.3e-4 con 40. En el armónico,
   ``û = F/(k(P) − Ω²m)``; en el espectral, la cabeza ``Sₐ·m/k(P)``.
2. Por encima de la carga de pandeo (1.2·P_cr) el autovalor es negativo:
   se reporta con su signo, con frecuencia NaN y un aviso; el espectral se
   rechaza.
3. Acero que fluyó: la tangente del estado es la de descarga (elástica) y el
   periodo no cambia; concreto agrietado: el periodo se alarga.
4. El análisis no cambia el estado: una etapa posterior da lo mismo.
5. Con ``previous`` hace falta el ensamblador de la etapa.
"""
import logging
import math
import unittest

import numpy as np

import solidum
from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler, Frame2DEulerCorot
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.solvers import (
    HarmonicSolver,
    IndirectDisplacementSolver,
    ModalSolver,
    NonlinearSolver,
    ResponseSpectrumSolver,
)
from solidum.sections import FiberSection, RectangularPatch

E, L, B, H = 2.0e11, 3.0, 0.2, 0.3
A, I = B * H, B * H ** 3 / 12
M_TIP = 4.0e3
# Columna esbelta del criterio 1 (L/r ≈ 104): P/EA ≈ 1e-4, despreciable frente
# a la tolerancia. La fórmula clásica supone la columna inextensible; con una
# robusta (P/EA ≈ 1e-3) el acortamiento axial cambia k en el mismo orden.
H_SLENDER = 0.1
A_S, I_S = B * H_SLENDER, B * H_SLENDER ** 3 / 12
P_CR = math.pi ** 2 * E * I_S / (4 * L ** 2)


def column(n=20):
    dom = Domain()
    mat = Elastic1D(E=E, density=0.0)
    nodes = [dom.add_node(i + 1, [0.0, i * L / n]) for i in range(n + 1)]
    for k in range(n):
        dom.add_element(Frame2DEulerCorot(k + 1, [nodes[k], nodes[k + 1]], material=mat,
                                          A=A_S, I=I_S))
    for d in ("ux", "uy", "rz"):
        nodes[0].fix_dof(d, 0.0)
    dom.add_nodal_mass(nodes[-1].id, ux=M_TIP)
    dom.generate_equation_numbers()
    return dom, nodes[-1]


def k_lateral(P):
    """Rigidez lateral exacta de la cabeza de un voladizo comprimido por P."""
    kappa = math.sqrt(P / (E * I_S))
    return P * kappa / (math.tan(kappa * L) - kappa * L)


def axial_stage(P, n=20):
    dom, top = column(n)
    asm = Assembler(dom)
    F = np.zeros(dom.total_dofs)
    F[top.dofs["uy"]] = -P
    r = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=5), F_applied=F)
    return dom, top, asm, r


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


class TestColumnaComprimida(_Quiet):
    """Criterio 1."""

    def test_modal(self):
        dom, top, asm, r = axial_stage(0.5 * P_CR)
        elastic = solidum.run_modal(dom, assembler=asm, n_modes=1)
        staged = solidum.run_modal(dom, assembler=asm, n_modes=1, previous=r)
        w_exact = math.sqrt(k_lateral(0.5 * P_CR) / M_TIP)
        self.assertAlmostEqual(elastic.frequencies_rad[0] / math.sqrt(3 * E * I_S / L ** 3 / M_TIP),
                               1.0, delta=1e-9)
        self.assertAlmostEqual(staged.frequencies_rad[0] / w_exact, 1.0, delta=1e-3)
        self.assertEqual(staged.unstable_modes.size, 0)

    def test_armonico(self):
        dom, top, asm, r = axial_stage(0.5 * P_CR)
        k = k_lateral(0.5 * P_CR)
        Omega = np.array([0.3, 0.6]) * math.sqrt(k / M_TIP)
        F = np.zeros(dom.total_dofs, dtype=complex)
        F[top.dofs["ux"]] = 1.0e3
        res = solidum.run_harmonic(dom, assembler=asm, solver=HarmonicSolver(
            asm, omega=Omega, F_amplitude=F), previous=r)
        # El error de k (4e-4 con 20 elementos, cuadrático con la malla) se
        # amplifica por 1/(1 − Ω²/ω²) = 1.56 en Ω = 0.6·ω.
        np.testing.assert_allclose(res.u_complex[top.dofs["ux"]].real,
                                   1.0e3 / (k - Omega ** 2 * M_TIP), rtol=2e-3)

    def test_espectral(self):
        dom, top, asm, r = axial_stage(0.5 * P_CR)
        Sa = 2.0
        res = solidum.run_response_spectrum(dom, assembler=asm, solver=ResponseSpectrumSolver(
            asm, n_modes=1, direction={"dof_name": "ux"},
            spectrum={"type": "constant_acceleration", "Sa": Sa}), previous=r)
        self.assertAlmostEqual(res.u_combined[top.dofs["ux"]]
                               / (Sa * M_TIP / k_lateral(0.5 * P_CR)), 1.0, delta=1e-3)
        # La masa no cambia: el cortante basal es m·Sₐ.
        self.assertAlmostEqual(res.total_reaction("ux") / (M_TIP * Sa), 1.0, delta=1e-9)


class TestInestable(unittest.TestCase):
    """Criterio 2."""

    def test_autovalor_negativo(self):
        logging.disable(logging.NOTSET)
        dom, top, asm, r = axial_stage(1.2 * P_CR)
        with self.assertLogs("solidum", level="WARNING") as logs:
            res = solidum.run_modal(dom, assembler=asm, n_modes=1, previous=r)
        self.assertTrue(any("autovalor(es) negativo(s)" in m for m in logs.output))
        self.assertLess(res.eigenvalues[0], 0.0)
        self.assertTrue(math.isnan(res.frequencies_rad[0]))
        self.assertTrue(math.isnan(res.periods[0]))
        np.testing.assert_array_equal(res.unstable_modes, [0])
        with self.assertRaisesRegex(ValueError, "modos inestables"):
            solidum.run_response_spectrum(dom, assembler=asm, solver=ResponseSpectrumSolver(
                asm, n_modes=1, direction={"dof_name": "ux"},
                spectrum={"type": "constant_acceleration", "Sa": 1.0}), previous=r)


def steel_cantilever():
    steel = Elastoplastic1D(E=E, sigma_y=2.5e8, H=2.0e9, density=7850.0)
    sec = FiberSection(patches=[RectangularPatch(steel, y=(-H / 2, H / 2), z=(-B / 2, B / 2),
                                                 ny=20, nz=1)])
    dom = Domain()
    nodes = [dom.add_node(i + 1, [0.0, i * L / 6]) for i in range(7)]
    for k in range(6):
        dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], section=sec, n_ip=3))
    for d in ("ux", "uy", "rz"):
        nodes[0].fix_dof(d, 0.0)
    dom.add_nodal_mass(nodes[-1].id, ux=5000.0, uy=5000.0)
    dom.generate_equation_numbers()
    return dom, nodes[-1]


class TestMaterialNoLineal(_Quiet):
    """Criterios 3 y 4."""

    def _push(self, dom, top, asm, target):
        F = np.zeros(dom.total_dofs)
        F[top.dofs["ux"]] = 1.0e5
        s = IndirectDisplacementSolver(asm, control=[(top.id, "ux", 1.0)], control_target=target,
                                       initial_dl=0.004, max_lambda=1e9, max_steps=200)
        return solidum.run(dom, assembler=asm, solver=s, F_applied=F)

    def test_acero_descarga_elastica_y_estado_intacto(self):
        dom, top = steel_cantilever()
        asm = Assembler(dom)
        T0 = solidum.run_modal(dom, assembler=asm, n_modes=1).periods[0]
        r = self._push(dom, top, asm, 0.04)
        self.assertGreater(abs(r.load_factor), 0.0)
        T1 = solidum.run_modal(dom, assembler=asm, n_modes=1, previous=r).periods[0]
        self.assertAlmostEqual(T1 / T0, 1.0, delta=1e-9)
        # El modal no cambió el estado: seguir empujando da lo mismo.
        dom2, top2 = steel_cantilever()
        asm2 = Assembler(dom2)
        r2 = self._push(dom2, top2, asm2, 0.04)
        a = self._push_from(dom, top, asm, r, 0.06)
        b = self._push_from(dom2, top2, asm2, r2, 0.06)
        np.testing.assert_allclose(a.U, b.U, rtol=0, atol=1e-12)

    def _push_from(self, dom, top, asm, previous, target):
        F = np.zeros(dom.total_dofs)
        F[top.dofs["ux"]] = 1.0e5
        s = IndirectDisplacementSolver(asm, control=[(top.id, "ux", 1.0)], control_target=target,
                                       initial_dl=0.004, max_lambda=1e9, max_steps=200)
        return solidum.run(dom, assembler=asm, solver=s, F_applied=F, previous=previous)

    def test_concreto_agrietado(self):
        import sys
        from pathlib import Path
        sys.path.insert(0, str(Path(__file__).resolve().parent / "validation"))
        from test_mander_concrete_validation import cantilever, make_section
        dom, top, F_g, F_h = cantilever(8, make_section(ny=30, density=True), N_axial=-6.0e5)
        dom.add_nodal_mass(top.id, ux=2.0e4, uy=2.0e4)
        asm = Assembler(dom)
        T0 = solidum.run_modal(dom, assembler=asm, n_modes=1).periods[0]
        prev = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=4),
                           F_applied=F_g)
        s = IndirectDisplacementSolver(asm, control=[(top.id, "ux", 1.0)], control_target=0.02,
                                       initial_dl=1e-3, max_lambda=1e9, max_steps=400)
        prev = solidum.run(dom, assembler=asm, solver=s, F_applied=F_h, previous=prev)
        T1 = solidum.run_modal(dom, assembler=asm, n_modes=1, previous=prev).periods[0]
        self.assertGreater(T1 / T0, 1.3)


class TestEnsambladorDeLaEtapa(_Quiet):
    """Criterio 5."""

    def test_hace_falta(self):
        dom, top, asm, r = axial_stage(0.3 * P_CR, n=4)
        with self.assertRaisesRegex(ValueError, "pase también el 'assembler'"):
            solidum.run_modal(dom, n_modes=1, previous=r)
        # Con el solver (que trae su ensamblador) basta.
        res = solidum.run_modal(dom, solver=ModalSolver(asm, n_modes=1), previous=r)
        self.assertEqual(res.n_modes, 1)


if __name__ == "__main__":
    unittest.main()
