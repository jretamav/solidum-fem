"""Tests de ``IndirectDisplacementSolver`` (control indirecto de desplazamiento).

Referencia
----------
de Borst, R. (1987). Computation of post-bifurcation and post-failure
behavior of strain-softening solids. *Computers & Structures* 25(2), 211-224.
Crisfield, M. A. (1991). *Non-linear Finite Element Analysis of Solids and
Structures*, vol. 1, cap. 9.

Caso de validación
------------------
Barra de longitud ``L`` y sección ``A`` en ``n`` elementos ``Truss2D`` con
``IsotropicDamage1D`` (ablandamiento exponencial), empotrada en ``x = 0`` y
cargada en el extremo. Un elemento tiene un umbral ``κ₀`` un 5 % menor: el
daño se localiza en él y el resto se descarga elásticamente. El esfuerzo es
uniforme, así que la respuesta exacta se escribe con la deformación ``ε_w`` del
elemento débil como parámetro::

    σ(ε_w) = (1 − d(ε_w))·E·ε_w,     d = 1 − (κ₀w/ε_w)·exp(−α(ε_w − κ₀w))
    F = A·σ,                          u = h·ε_w + (L − h)·σ/E,   h = L/n

(con el tope ``DAMAGE_MAX`` del material sobre ``d``). Hay retroceso
(*snap-back*) si ``(n − 1)·α·κ₀w > 1``: con ``α·κ₀w = 0.475`` la barra de 10
elementos lo tiene y la de 2 no. El control de carga se detiene en el pico, el
de desplazamiento del extremo no pasa el retroceso y el arco cilíndrico
tampoco (medido 2026-09-23); el control del alargamiento del elemento débil
lo sigue exacto.
"""
from __future__ import annotations

import logging
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.constants import DAMAGE_MAX
from solidum.core.domain import Domain
from solidum.elements.truss import Truss2D
from solidum.materials.damage_1d import IsotropicDamage1D
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.solvers import IndirectDisplacementSolver
from solidum.utils.yaml_parser import YamlParser

E, A, L = 30.0e9, 0.01, 1.0
K0, ALPHA, WEAK = 1.0e-4, 5000.0, 0.95
K0W = WEAK * K0
F_PEAK = A * E * K0W


def _bar(n: int):
    """Barra con el elemento central debilitado; carga de referencia = pico."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * L / n, 0.0]) for i in range(n + 1)]
    for k in range(n):
        kappa_0 = K0W if k == n // 2 else K0
        mat = IsotropicDamage1D(E=E, kappa_0=kappa_0, alpha=ALPHA)
        dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]], mat, A=A))
    nodes[0].fix_dof('ux', 0.0)
    for node in nodes:
        node.fix_dof('uy', 0.0)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[nodes[-1].dofs['ux']] = F_PEAK
    weak = dom.elements[n // 2 + 1]
    return dom, Assembler(dom), F, nodes[-1], weak


def _exact_curve(n: int):
    """(u, F) exactos en la rama de ablandamiento, densamente muestreados."""
    h = L / n
    ew = K0W * np.concatenate([np.linspace(1.0, 3.0, 20000), np.geomspace(3.0, 60.0, 20000)])
    d = np.minimum(1.0 - (K0W / ew) * np.exp(-ALPHA * (ew - K0W)), DAMAGE_MAX)
    sigma = (1.0 - d) * E * ew
    return h * ew + (L - h) * sigma / E, A * sigma


def _trace(n: int, **kw):
    dom, asm, F, end, weak = _bar(n)
    a, b = weak.nodes
    kw.setdefault('initial_dlambda', 0.02)
    solver = IndirectDisplacementSolver(
        asm, max_iter=30, max_lambda=1.2, max_steps=kw.pop('max_steps', 400),
        control=[(b.id, 'ux', 1.0), (a.id, 'ux', -1.0)], **kw)
    hist = []

    def cb(step, U, lam):
        damaged = sum(1 for el in dom.elements.values() if el.state.vars[0]['damage'] > 0.0)
        hist.append((U[end.dofs['ux']], lam * F_PEAK, damaged,
                     U[b.dofs['ux']] - U[a.dofs['ux']]))

    solver.solve(F, step_callback=cb)
    return solver, np.array(hist)


class _Quiet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        logging.disable(logging.CRITICAL)

    @classmethod
    def tearDownClass(cls):
        logging.disable(logging.NOTSET)


class TestSnapBackBar(_Quiet):
    """Barra de 10 elementos: retroceso con u_min = 0.616·u_pico."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.solver, cls.H = _trace(10)
        cls.u_ex, cls.F_ex = _exact_curve(10)

    def test_only_the_weak_element_damages(self):
        self.assertEqual(int(self.H[:, 2].max()), 1)

    def test_post_peak_points_on_exact_curve(self):
        u_pk = self.u_ex[0]
        post = self.H[(self.H[:, 2] >= 1) & (self.H[:, 1] > 1.5 * self.F_ex.min())]
        self.assertGreater(len(post), 100)
        dist = [np.min(np.hypot((self.u_ex - u) / u_pk, (self.F_ex - f) / F_PEAK))
                for u, f in post[:, :2]]
        self.assertLess(max(dist), 1e-3)

    def test_displacement_retreats_to_exact_minimum(self):
        u_pk = self.u_ex[0]
        u_min_exact = self.u_ex[self.F_ex < 0.999 * F_PEAK].min()
        u_min = self.H[self.H[:, 2] >= 1, 0].min()
        self.assertLess(u_min, 0.7 * u_pk)                 # retrocede de verdad
        self.assertAlmostEqual(u_min / u_pk, u_min_exact / u_pk, delta=2e-3)

    def test_controlled_quantity_grows_by_dl_each_step(self):
        """La restricción: el alargamiento del elemento débil crece Δl por paso."""
        dl = self.solver.dl_reference
        steps = np.diff(np.concatenate([[0.0], self.H[:, 3]]))
        np.testing.assert_allclose(steps, dl, rtol=1e-8)


class TestNoSnapBackBar(_Quiet):
    """Barra de 2 elementos: rama descendente sin retroceso."""

    def test_matches_exact_curve(self):
        _, H = _trace(2)
        u_ex, F_ex = _exact_curve(2)
        u_pk = u_ex[0]
        self.assertEqual(int(H[:, 2].max()), 1)
        post = H[(H[:, 2] >= 1) & (H[:, 1] > 1.5 * F_ex.min())]
        dist = [np.min(np.hypot((u_ex - u) / u_pk, (F_ex - f) / F_PEAK)) for u, f in post[:, :2]]
        self.assertLess(max(dist), 2e-3)
        # Rama elástica exacta antes del pico.
        pre = H[H[:, 2] == 0]
        np.testing.assert_allclose(pre[:, 1], E * A * pre[:, 0] / L, rtol=1e-10)


class TestLinearAndFirstStep(_Quiet):

    def _elastic_bar(self):
        dom = Domain()
        n1, n2, n3 = (dom.add_node(i + 1, [0.5 * i, 0.0]) for i in range(3))
        for k, (a, b) in enumerate(((n1, n2), (n2, n3))):
            dom.add_element(Truss2D(k + 1, [a, b], Elastic1D(E=E), A=A))
        n1.fix_dof('ux', 0.0)
        for node in (n1, n2, n3):
            node.fix_dof('uy', 0.0)
        dom.generate_equation_numbers(verbose=False)
        F = np.zeros(dom.total_dofs)
        F[n3.dofs['ux']] = 1.0e5
        return dom, Assembler(dom), F, n2, n3

    def test_first_step_is_a_load_fraction(self):
        dom, asm, F, n2, n3 = self._elastic_bar()
        solver = IndirectDisplacementSolver(
            asm, control=[(n3, 'ux', 1.0)], initial_dlambda=0.25, max_lambda=1.0)
        lams = []
        U = solver.solve(F, step_callback=lambda k, U, lam: lams.append(lam))
        self.assertAlmostEqual(lams[0], 0.25, places=12)
        np.testing.assert_allclose(np.diff([0.0] + lams), 0.25, rtol=1e-10)
        # Solución lineal exacta en λ = 1.
        self.assertAlmostEqual(U[n3.dofs['ux']], 1.0e5 * L / (E * A), delta=1e-12 * L)

    def test_decreasing_quantity_is_a_clear_error(self):
        dom, asm, F, n2, n3 = self._elastic_bar()
        solver = IndirectDisplacementSolver(asm, control=[(n3, 'ux', -1.0)])
        with self.assertRaisesRegex(ValueError, "no crece"):
            solver.solve(F)


class TestValidation(_Quiet):

    def _asm(self):
        dom, asm, F, end, weak = _bar(2)
        return dom, asm, F, end

    def test_malformed_control(self):
        dom, asm, F, end = self._asm()
        for bad in ([], {'node': 1, 'dof': 'ux'}, [(1, 'ux')], [(1, 'ux', 0.0)],
                    [{'node': 1, 'dir': 'ux'}], [(1, 3, 1.0)]):
            with self.subTest(control=bad), self.assertRaises(ValueError):
                IndirectDisplacementSolver(asm, control=bad)

    def test_unknown_node_or_dof(self):
        dom, asm, F, end = self._asm()
        for bad, pattern in (([(99, 'ux', 1.0)], "no existe"),
                             ([(end.id, 'uz', 1.0)], "no tiene")):
            with self.subTest(control=bad):
                with self.assertRaisesRegex(ValueError, pattern):
                    IndirectDisplacementSolver(asm, control=bad).solve(F)

    def test_constrained_dof_is_rejected(self):
        dom, asm, F, end = self._asm()
        with self.assertRaisesRegex(ValueError, "restringido"):
            IndirectDisplacementSolver(asm, control=[(1, 'ux', 1.0)]).solve(F)

    def test_cancelling_terms(self):
        dom, asm, F, end = self._asm()
        with self.assertRaisesRegex(ValueError, "se anulan"):
            IndirectDisplacementSolver(
                asm, control=[(end.id, 'ux', 1.0), (end.id, 'ux', -1.0)]).solve(F)

    def test_from_yaml(self):
        yaml_txt = (
            "nodes:\n  - {id: 1, coords: [0.0, 0.0]}\n  - {id: 2, coords: [0.5, 0.0]}\n"
            "  - {id: 3, coords: [1.0, 0.0]}\n"
            "materials:\n  - {id: 1, type: Elastic1D, E: 1000.0}\n"
            "elements:\n  - {id: 1, type: Truss2D, material: 1, A: 1.0, nodes: [1, 2]}\n"
            "  - {id: 2, type: Truss2D, material: 1, A: 1.0, nodes: [2, 3]}\n"
            "boundary_conditions_by_node:\n  - {node_id: 1, ux: 0.0, uy: 0.0}\n"
            "  - {node_id: 2, uy: 0.0}\n  - {node_id: 3, uy: 0.0}\n"
            "point_loads_by_node:\n  - {node_id: 3, ux: 10.0}\n"
            "solver:\n  type: IndirectDisplacementSolver\n  initial_dlambda: 0.5\n"
            "  control:\n    - {node: 3, dof: ux, coef: 1.0}\n    - {node: 2, dof: ux, coef: -1.0}\n"
        )
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "m.yaml")
            with open(path, "w", encoding="utf-8") as f:
                f.write(yaml_txt)
            parser = YamlParser(path)
            dom = parser.parse()
            dom.generate_equation_numbers()
            asm = Assembler(dom)
            solver = parser.get_solver(asm)
            U = solver.solve(parser.get_external_forces())
        self.assertIsInstance(solver, IndirectDisplacementSolver)
        self.assertAlmostEqual(U[dom.get_node(3).dofs['ux']], 10.0 * 1.0 / 1000.0, places=12)


if __name__ == '__main__':
    unittest.main()
