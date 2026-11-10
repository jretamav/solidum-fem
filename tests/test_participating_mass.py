"""Masa de referencia de la razón de masa efectiva acumulada (decisión del
usuario del 2026-11-09, ADR 0009 §8): la masa que puede participar,
``bᵀ·M_LL⁻¹·b`` con ``b = Tᵀ·M·r``, la suma de las masas efectivas de todos
los modos.

Criterios:

1. Con todos los modos la razón acumulada vale 1 (a 1e-10) en un marco con
   masa consistente y concentrada y mallas de 2 y 6 elementos por barra, en
   el modal y en el espectral. Con la masa total del modelo (``rᵀ·M·r``,
   la referencia anterior) no llegaba a 1 (0.875 con 2 elementos por barra y
   masa concentrada).
2. Con masa concentrada es la masa de los grados de libertad libres.
3. Una masa de nudo en un apoyo no cuenta: la razón llega a 1 (antes, 0.333
   en el caso de la auditoría con dos tercios de la masa en el apoyo).
4. Con grados de libertad sin masa y uniones multipunto (edificio de
   cortante) es la suma de las masas de piso.
"""
import logging
import unittest

import numpy as np

from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.math.solvers import ModalSolver, ResponseSpectrumSolver
from test_nodal_mass import shear_building
from test_transient_massless_dofs import HEIGHTS, INERTIAS, MASSES

E, RHO, A, I = 2.0e11, 7850.0, 0.15, 3.1e-3


def portal(n_el):
    """Marco de un nivel: columnas de 3 m y viga de 6 m, ``n_el`` elementos por barra."""
    dom = Domain()
    mat = Elastic1D(E=E, density=RHO)
    counter = {"n": 0, "e": 0}

    def node(x, y):
        counter["n"] += 1
        return dom.add_node(counter["n"], [x, y])

    def bar(a, b):
        (xa, ya), (xb, yb) = a.coordinates, b.coordinates
        pts = [a] + [node(xa + (xb - xa) * j / n_el, ya + (yb - ya) * j / n_el)
                     for j in range(1, n_el)] + [b]
        for p, q in zip(pts[:-1], pts[1:]):
            counter["e"] += 1
            dom.add_element(Frame2DEuler(counter["e"], [p, q], material=mat, A=A, I=I))

    n1, n2, n3, n4 = node(0.0, 0.0), node(0.0, 3.0), node(6.0, 3.0), node(6.0, 0.0)
    bar(n1, n2), bar(n2, n3), bar(n3, n4)
    for n in (n1, n4):
        for d in ("ux", "uy", "rz"):
            n.fix_dof(d, 0.0)
    dom.generate_equation_numbers()
    return dom


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


class TestTodosLosModos(_Quiet):
    """Criterios 1 y 2."""

    def test_razon_igual_a_uno(self):
        for n_el in (2, 6):
            for lumping in ("consistent", "lumped"):
                with self.subTest(n_el=n_el, lumping=lumping):
                    dom = portal(n_el)
                    n_free = len(Assembler(dom).constraint_set.free_dofs(dom.total_dofs))
                    modal = ModalSolver(Assembler(dom), n_modes=n_free, lumping=lumping).solve()
                    for d in ("ux", "uy"):
                        ratio = modal.cumulative_effective_mass_ratio(d)
                        self.assertAlmostEqual(ratio[-1], 1.0, delta=1e-10)
                        self.assertTrue(np.all(np.diff(ratio) >= -1e-14))
                    spec = ResponseSpectrumSolver(
                        Assembler(dom), n_modes=n_free, direction={"dof_name": "ux"},
                        spectrum=lambda w: 1.0, lumping=lumping).solve()
                    self.assertAlmostEqual(spec.cumulative_effective_mass_ratio()[-1], 1.0,
                                           delta=1e-10)
                    self.assertAlmostEqual(spec.participating_mass / modal.participating_mass["ux"],
                                           1.0, delta=1e-12)
                    self.assertAlmostEqual(spec.total_mass / modal.total_mass["ux"], 1.0,
                                           delta=1e-12)
                    if lumping == "lumped":
                        M = Assembler(dom).assemble_mass_matrix("lumped").toarray()
                        free = Assembler(dom).constraint_set.free_dofs(dom.total_dofs)
                        r = np.zeros(dom.total_dofs)
                        for n in dom.nodes.values():
                            r[n.dofs["ux"]] = 1.0
                        m_free = float(r[free] @ M[np.ix_(free, free)] @ r[free])
                        self.assertAlmostEqual(modal.participating_mass["ux"] / m_free, 1.0,
                                               delta=1e-12)
                        if n_el == 2:
                            # La referencia anterior: 0.875 con todos los modos.
                            old = float(np.sum(modal.effective_masses["ux"])) / modal.total_mass["ux"]
                            self.assertAlmostEqual(old, 0.875, delta=1e-3)


class TestMasaEnUnApoyo(_Quiet):
    """Criterio 3: barra de dos elementos sin densidad, masa de nudo en el
    extremo libre y el doble en el apoyo."""

    def test_razon_igual_a_uno(self):
        dom = Domain()
        nodes = [dom.add_node(i + 1, [float(i), 0.0]) for i in range(3)]
        mat = Elastic1D(E=E, density=0.0)
        for k in range(2):
            dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]], mat, A=1e-3))
        nodes[0].fix_dof("ux", 0.0)
        for n in nodes:
            n.fix_dof("uy", 0.0)
        dom.add_nodal_mass(1, ux=200.0)
        dom.add_nodal_mass(2, ux=50.0)
        dom.add_nodal_mass(3, ux=50.0)
        dom.generate_equation_numbers()
        res = ResponseSpectrumSolver(Assembler(dom), n_modes=2, direction={"dof_name": "ux"},
                                     spectrum=lambda w: 1.0).solve()
        self.assertAlmostEqual(res.total_mass, 300.0, delta=1e-9)
        self.assertAlmostEqual(res.participating_mass, 100.0, delta=1e-9)
        self.assertAlmostEqual(res.cumulative_effective_mass_ratio()[-1], 1.0, delta=1e-12)


class TestEdificioDeCortante(_Quiet):
    """Criterio 4."""

    def test_masas_de_piso(self):
        dom, _ = shear_building(HEIGHTS, MASSES, INERTIAS)
        modal = ModalSolver(Assembler(dom), n_modes=3).solve()
        self.assertAlmostEqual(modal.participating_mass["ux"] / sum(MASSES), 1.0, delta=1e-12)
        self.assertAlmostEqual(modal.cumulative_effective_mass_ratio("ux")[-1], 1.0, delta=1e-10)
        self.assertEqual(modal.participating_mass["uy"], 0.0)
        np.testing.assert_array_equal(modal.cumulative_effective_mass_ratio("uy"), 0.0)


if __name__ == "__main__":
    unittest.main()
