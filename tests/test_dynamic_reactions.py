"""Reacciones de los apoyos en los análisis dinámicos (decisión del usuario del
2026-11-09): la ecuación de movimiento en las filas de apoyo,
``R = F_int + M·ü + C·u̇ − F``, con la parte de amortiguamiento ``C·u̇`` aparte.

Criterios:

1. Transitorio: la suma de las reacciones en ``x`` equilibra la estructura
   entera, ``Σₓ R = Σₓ (M·ü + C·u̇) − Σₓ F`` (``Σₓ K·u = 0`` por traslación
   rígida), a 1e-10 con Newmark amortiguado y con la diferencia central sin
   amortiguar, y a la tolerancia del Newton con Newton-Newmark.
2. Edificio de cortante (columnas sin masa, uniones multipunto): el cortante
   basal es ``−k₁·(u₁ + a₁·u̇₁)`` con la historia del problema 3×3, y sin la
   parte de amortiguamiento ``−k₁·u₁``, a 1e-9.
3. Armónico: voladizo sin masa con masa en la cabeza y amortiguamiento
   proporcional a la rigidez: ``R̂ = −F̂·(k + iΩa₁k)/(k − Ω²m + iΩa₁k)`` y su
   parte de amortiguamiento ``−iΩa₁k·û``, a 1e-10.
4. Espectral: el cortante basal de cada modo es ``−Mₙ*·Sₐ`` y el combinado
   ``√(Σ (Mₙ*·Sₐ)²)`` (SRSS), a 1e-9; ``total_reaction`` suma antes de combinar.
"""
import logging
import math
import unittest

import numpy as np

from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler
from solidum.materials.elastic import Elastic1D
from solidum.math.solvers import (
    CentralDifferenceSolver,
    HarmonicSolver,
    NewmarkSolver,
    NewtonNewmarkSolver,
    ResponseSpectrumSolver,
)
from test_nodal_mass import E as E_SB, shear_building, shear_building_3x3
from test_transient_continuation import tight
from test_transient_massless_dofs import HEIGHTS, INERTIAS, MASSES, _oracle, _rayleigh

E, RHO, L, B, H, N = 2.0e11, 7850.0, 3.0, 0.3, 0.5, 6
A, I = B * H, B * H ** 3 / 12


def column(density=RHO, tip_mass=2000.0):
    dom = Domain()
    mat = Elastic1D(E=E, density=density)
    nodes = [dom.add_node(i + 1, [0.0, i * L / N]) for i in range(N + 1)]
    for k in range(N):
        dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], material=mat, A=A, I=I))
    for d in ("ux", "uy", "rz"):
        nodes[0].fix_dof(d, 0.0)
    if tip_mass:
        dom.add_nodal_mass(nodes[-1].id, ux=tip_mass, uy=tip_mass)
    dom.generate_equation_numbers()
    return dom, nodes


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


class TestEquilibrioGlobal(_Quiet):
    """Criterio 1."""

    def _check(self, make, rtol, rayleigh, lumping="consistent"):
        dom, nodes = column()
        top = nodes[-1].dofs["ux"]

        def F_func(t):
            F = np.zeros(dom.total_dofs)
            F[top] = 1.0e5 * math.sin(2 * math.pi * t / 0.02) if t < 0.01 else 0.0
            return F
        asm = Assembler(dom)
        res = make(asm, F_func, rayleigh).solve()
        M = asm.assemble_mass_matrix(lumping=lumping).toarray()
        a0, a1 = res.alpha_rayleigh, res.beta_rayleigh
        asm.assemble_system()
        C = a0 * M + a1 * asm.K_global.toarray()
        ux = [n.dofs["ux"] for n in nodes]
        Fh = np.column_stack([F_func(t) for t in res.t_history])
        expected = ((M @ res.uddot_history + C @ res.udot_history)[ux].sum(axis=0)
                    - Fh[ux].sum(axis=0))
        V = res.total_reaction("ux")
        np.testing.assert_allclose(V, expected, rtol=0.0, atol=rtol * np.max(np.abs(V)))
        self.assertEqual(set(res.reactions_by_node()), {nodes[0].id})
        self.assertEqual(set(res.reactions_by_node()[nodes[0].id]), {"ux", "uy", "rz"})
        return res

    def test_newmark_amortiguado(self):
        res = self._check(lambda asm, F, r: NewmarkSolver(asm, t_end=0.05, dt=2e-4, F_func=F,
                                                          rayleigh=r),
                          1e-10, {"alpha": 2.0, "beta": 5.0e-5})
        # Sin amortiguamiento la reacción difiere justo en su parte C·u̇.
        np.testing.assert_allclose(
            res.total_reaction("ux") - res.total_reaction("ux", damping=False),
            res.R_damping_history[[i for i, (_, d) in enumerate(res.reaction_labels)
                                   if d == "ux"]].sum(axis=0), rtol=0.0, atol=1e-6)

    def test_diferencia_central_sin_amortiguar(self):
        self._check(lambda asm, F, r: CentralDifferenceSolver(asm, t_end=5e-3, dt=2e-6,
                                                              F_func=F, lumping="lumped"),
                    1e-10, None, lumping="lumped")

    def test_newton_newmark(self):
        self._check(lambda asm, F, r: NewtonNewmarkSolver(asm, t_end=0.02, dt=2e-4, F_func=F,
                                                          rayleigh=r, convergence=tight()),
                    1e-8, {"alpha": 2.0, "beta": 5.0e-5})


class TestEdificioDeCortante(_Quiet):
    """Criterio 2."""

    def test_cortante_basal(self):
        dt, t_end, P = 0.01, 1.0, 1.0e5
        dom, floors = shear_building(HEIGHTS, MASSES, INERTIAS)
        dom.generate_equation_numbers()
        roof = floors[-1].dofs["ux"]

        def F_func(t):
            F = np.zeros(dom.total_dofs)
            F[roof] = P
            return F
        a0, a1 = _rayleigh(0.05)
        res = NewmarkSolver(Assembler(dom), t_end=t_end, dt=dt, F_func=F_func,
                            rayleigh={"alpha": a0, "beta": a1}).solve()
        K3, M3 = shear_building_3x3(HEIGHTS, MASSES, INERTIAS)
        u_ref, _ = _oracle(0.0, 0.25, 0.5, a0 * M3 + a1 * K3, dt, res.n_steps)
        k1 = 2.0 * 12.0 * E_SB * INERTIAS[0] / HEIGHTS[0] ** 3
        u1 = res.u_history[floors[0].dofs["ux"]]
        v1 = res.udot_history[floors[0].dofs["ux"]]
        np.testing.assert_allclose(u1, u_ref[0], rtol=0.0, atol=1e-9 * np.max(np.abs(u_ref)))
        V = res.total_reaction("ux")
        scale = k1 * np.max(np.abs(u_ref[0]))
        np.testing.assert_allclose(V, -k1 * (u1 + a1 * v1), rtol=0.0, atol=1e-9 * scale)
        np.testing.assert_allclose(res.total_reaction("ux", damping=False), -k1 * u1,
                                   rtol=0.0, atol=1e-9 * scale)


class TestArmonico(_Quiet):
    """Criterio 3."""

    def test_oscilador_con_amortiguamiento_de_rigidez(self):
        m, a1, F0 = 2000.0, 2.0e-3, 1.0e4
        dom, nodes = column(density=0.0, tip_mass=m)
        k = 3.0 * E * I / L ** 3
        wn = math.sqrt(k / m)
        omega = np.array([0.3, 0.9, 1.0, 1.4]) * wn
        F = np.zeros(dom.total_dofs, dtype=complex)
        F[nodes[-1].dofs["ux"]] = F0
        res = HarmonicSolver(Assembler(dom), omega=omega, F_amplitude=F,
                             rayleigh={"alpha": 0.0, "beta": a1}).solve()
        u = res.u_complex[nodes[-1].dofs["ux"]]
        den = k - omega ** 2 * m + 1j * omega * a1 * k
        np.testing.assert_allclose(u, F0 / den, rtol=1e-10)
        R = res.total_reaction("ux")
        np.testing.assert_allclose(R, -F0 * (k + 1j * omega * a1 * k) / den, rtol=1e-10)
        R_damp = R - res.total_reaction("ux", damping=False)
        np.testing.assert_allclose(R_damp, -1j * omega * a1 * k * u, rtol=1e-10)


class TestEspectral(_Quiet):
    """Criterio 4."""

    def test_cortante_basal_por_modo(self):
        Sa = 2.45
        for combination in ("SRSS", "CQC"):
            with self.subTest(combination=combination):
                dom, floors = shear_building(HEIGHTS, MASSES, INERTIAS)
                res = ResponseSpectrumSolver(
                    Assembler(dom), n_modes=3, direction={"dof_name": "ux"},
                    spectrum={"type": "constant_acceleration", "Sa": Sa},
                    combination=combination, damping=0.05).solve()
                rows = [i for i, (_, d) in enumerate(res.reaction_labels) if d == "ux"]
                V_modes = res.R_per_mode[rows].sum(axis=0)
                # Vₙ = −Γₙ·(Σ m·φₙ)·Sₐ = −Mₙ*·Sₐ, sea cual sea el signo del modo.
                np.testing.assert_allclose(V_modes, -res.effective_masses * Sa, rtol=1e-9)
                if combination == "SRSS":
                    expected = math.sqrt(float(np.sum((res.effective_masses * Sa) ** 2)))
                    self.assertAlmostEqual(res.total_reaction("ux") / expected, 1.0, delta=1e-9)
                # Combinar antes de sumar daría otra cosa (la suma de máximos).
                self.assertLessEqual(res.total_reaction("ux"),
                                     float(np.sum(res.R_combined[rows])) * (1 + 1e-12))


if __name__ == "__main__":
    unittest.main()
