"""Cargas que llevan sus cargas equivalentes por elemento (ADR 0027).

Las fuerzas internas de extremo de un elemento con carga distribuida son
``F_int − f_eq``. Antes, el ensamblador recordaba sólo la última carga de
cuerpo ensamblada: dos cargas de cuerpo en una corrida, o dos casos de carga
ensamblados con el mismo ensamblador, daban fuerzas internas erróneas sin
aviso (deuda #35). Ahora cada carga lleva las suyas y la aritmética de las
cargas las combina.

Criterios:

- Álgebra: ``F1 ± F2``, ``a·F``, ``F/a`` y ``−F`` combinan linealmente las
  cargas equivalentes; un ``np.ndarray`` simple (cargas puntuales) no aporta;
  el resto de operaciones da un ``np.ndarray`` simple. Copia y serialización
  las conservan.
- Dos cargas de cuerpo en una corrida (peso propio y una fuerza de cuerpo
  lateral): las fuerzas internas coinciden con las de sumar a mano las cargas
  equivalentes de ambas (≤ 1e-10), y la fuerza axial varía a lo largo de la
  columna por el peso.
- Dos casos de carga ensamblados con el mismo ensamblador: cada caso descuenta
  sólo las suyas (igual que con un ensamblador propio, ≤ 1e-10).
- Un ``np.ndarray`` simple con varias cargas de cuerpo ensambladas avisa.
- Transitorio: con ``F(t) = F_peso + a(t)·F_lateral``, el historial de fuerzas
  internas descuenta en cada instante ``f_eq,peso + a(t)·f_eq,lateral``; con un
  ``np.ndarray`` simple avisa.
"""
import copy
import logging
import pickle
import unittest

import numpy as np

import solidum
from solidum import Assembler, Domain, LoadVector
from solidum.core.load_vector import EquivalentLoads
from solidum.elements.frame import Frame2DEuler
from solidum.materials.elastic import Elastic1D
from solidum.math.solvers import LinearSolver, NewmarkSolver
from solidum.results import build_solve_result

E, A, I, RHO = 2.0e11, 0.02, 8.0e-5, 7850.0
G = np.array([0.0, -9.81])
B_LAT = np.array([2.0e4, 0.0])


def column(n=4, length=3.0):
    dom = Domain()
    nodes = [dom.add_node(k + 1, [0.0, k * length / n]) for k in range(n + 1)]
    mat = Elastic1D(E=E, density=RHO)
    for k in range(n):
        dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], material=mat, A=A, I=I))
    for dof in ("ux", "uy", "rz"):
        nodes[0].fix_dof(dof, 0.0)
    dom.generate_equation_numbers()
    return dom, nodes[-1]


def forces(result, names=("N", "V", "M")):
    return {eid: np.concatenate([np.atleast_1d(np.asarray(ef.components[k], dtype=float))
                                 for k in names])
            for eid, ef in result.element_forces.items()}


class TestAlgebra(unittest.TestCase):

    def setUp(self):
        dom, self.top = column()
        self.asm = Assembler(dom)
        self.dom = dom
        self.F1 = self.asm.assemble_self_weight(G)
        self.F2 = self.asm.assemble_body_load(B_LAT)
        self.eq1 = self.F1.equivalent_loads.as_dict()
        self.eq2 = self.F2.equivalent_loads.as_dict()

    def check(self, F, a, b):
        self.assertIsInstance(F, LoadVector)
        got = F.equivalent_loads.as_dict()
        for eid in self.eq1:
            np.testing.assert_allclose(got[eid], a * self.eq1[eid] + b * self.eq2[eid],
                                       rtol=1e-14, atol=1e-12)

    def test_combinaciones_lineales(self):
        self.check(self.F1 + self.F2, 1.0, 1.0)
        self.check(self.F1 - self.F2, 1.0, -1.0)
        self.check(2.5 * self.F1, 2.5, 0.0)
        self.check(self.F1 * np.float64(0.5) + self.F2 / 4.0, 0.5, 0.25)
        self.check(-self.F2, 0.0, -1.0)
        self.check(1.2 * self.F1 + 1.6 * self.F2, 1.2, 1.6)

    def test_cargas_puntuales_no_aportan(self):
        P = np.zeros(self.dom.total_dofs)
        P[self.top.dofs["ux"]] = 1.0e3
        self.check(self.F1 + P, 1.0, 0.0)
        self.check(P + self.F1, 1.0, 0.0)
        F = self.F1.copy()
        F += self.F2
        self.check(F, 1.0, 1.0)
        F[self.top.dofs["ux"]] += 5.0
        self.check(F, 1.0, 1.0)
        np.testing.assert_allclose(np.asarray(F), np.asarray(self.F1) + np.asarray(self.F2)
                                   + 5.0 * (np.arange(F.size) == self.top.dofs["ux"]))

    def test_lo_que_no_es_una_carga_es_un_ndarray_simple(self):
        for value in (self.F1 * self.F2, np.abs(self.F1), self.F1[2:5],
                      np.asarray(self.F1), self.F1 > 0):
            self.assertNotIsInstance(value, LoadVector)
        self.assertIsInstance(float(np.linalg.norm(self.F1)), float)
        self.assertNotIsInstance(self.F1 @ self.F2, LoadVector)

    def test_copia_y_serializacion(self):
        for F in (copy.copy(self.F1 + self.F2), copy.deepcopy(self.F1 + self.F2),
                  pickle.loads(pickle.dumps(self.F1 + self.F2))):
            self.check(F, 1.0, 1.0)

    def test_desde_un_diccionario(self):
        F = solidum.core.load_vector.as_load_vector(np.asarray(self.F1), self.eq1)
        self.check(F + self.F2, 1.0, 1.0)


class TestDosCargasDeCuerpo(unittest.TestCase):
    """Deuda #35: peso propio y fuerza de cuerpo lateral en una corrida."""

    def reference(self, factors):
        """Fuerzas internas con las cargas equivalentes sumadas a mano."""
        dom, _ = column()
        asm = Assembler(dom)
        parts = [asm.assemble_self_weight(G), asm.assemble_body_load(B_LAT)]
        F = sum(np.asarray(f) * a for f, a in zip(parts, factors))
        eq = {}
        for f, a in zip(parts, factors):
            for eid, v in f.equivalent_loads.as_dict().items():
                eq[eid] = eq.get(eid, 0.0) + a * v
        U = LinearSolver(asm).solve(F)
        return forces(build_solve_result(dom, asm, U, np.asarray(F), equivalent_loads=eq))

    def test_igual_a_sumar_a_mano(self):
        dom, _ = column()
        asm = Assembler(dom)
        F = 1.2 * asm.assemble_self_weight(G) + 1.6 * asm.assemble_body_load(B_LAT)
        got = forces(solidum.run(dom, assembler=asm, F_applied=F))
        ref = self.reference((1.2, 1.6))
        for eid in ref:
            np.testing.assert_allclose(got[eid], ref[eid], rtol=1e-10,
                                       atol=1e-10 * np.max(np.abs(ref[eid])))

    def test_la_fuerza_axial_varia_con_el_peso(self):
        """En la columna vertical, N crece hacia la base por el peso propio:
        en cada elemento, N_i − N_j es el peso del elemento."""
        dom, _ = column()
        asm = Assembler(dom)
        F = asm.assemble_self_weight(G) + asm.assemble_body_load(B_LAT)
        res = solidum.run(dom, assembler=asm, F_applied=F)
        w_el = RHO * A * 9.81 * 3.0 / 4
        for ef in res.element_forces.values():
            N = np.asarray(ef.components["N"], dtype=float)
            self.assertAlmostEqual(N[-1] - N[0], w_el, delta=1e-9 * w_el)

    def test_dos_casos_con_el_mismo_ensamblador(self):
        """Cada caso descuenta sólo sus cargas, aunque se haya ensamblado otro
        después (antes descontaba las del último ensamblado)."""
        dom, _ = column()
        asm = Assembler(dom)
        F_a = asm.assemble_self_weight(G)
        F_b = asm.assemble_body_load(B_LAT)
        got_a = forces(solidum.run(dom, assembler=asm, F_applied=F_a))
        got_b = forces(solidum.run(dom, assembler=asm, F_applied=F_b))
        for case, factors in ((got_a, (1.0, 0.0)), (got_b, (0.0, 1.0))):
            ref = self.reference(factors)
            for eid in ref:
                np.testing.assert_allclose(case[eid], ref[eid], rtol=1e-10,
                                           atol=1e-10 * max(np.max(np.abs(ref[eid])), 1.0))

    def test_ndarray_simple_con_varias_cargas_avisa(self):
        dom, _ = column()
        asm = Assembler(dom)
        F = np.asarray(asm.assemble_self_weight(G)) + np.asarray(asm.assemble_body_load(B_LAT))
        with self.assertLogs("solidum", level=logging.WARNING) as logs:
            solidum.run(dom, assembler=asm, F_applied=F)
        self.assertTrue(any("np.ndarray simple" in m for m in logs.output))


class TestTransitorio(unittest.TestCase):

    def test_historial_descuenta_las_cargas_del_instante(self):
        dom, _ = column()
        asm = Assembler(dom)
        F_w = asm.assemble_self_weight(G)
        F_l = asm.assemble_body_load(B_LAT)

        def a(t):
            return np.sin(40.0 * t)

        solver = NewmarkSolver(asm, t_end=0.02, dt=1.0e-3, F_func=lambda t: F_w + a(t) * F_l,
                               record_internal_forces=True)
        res = solver.solve()
        eq_w, eq_l = F_w.equivalent_loads, F_l.equivalent_loads
        for k in (0, 7, 20):
            U, t = res.u_history[:, k], res.t_history[k]
            for eid, elem in dom.elements.items():
                f_eq = eq_w.get(eid) + a(t) * eq_l.get(eid)
                ref = elem.internal_forces(U, equivalent_load=f_eq)
                got = res.element_forces_history[eid][k]
                for name in ("N", "V", "M"):
                    np.testing.assert_allclose(
                        np.asarray(got.components[name], dtype=float),
                        np.asarray(ref.components[name], dtype=float), rtol=1e-12, atol=1e-9)

    def test_ndarray_simple_avisa(self):
        dom, _ = column()
        asm = Assembler(dom)
        F_w = np.asarray(asm.assemble_self_weight(G))
        solver = NewmarkSolver(asm, t_end=0.005, dt=1.0e-3, F_func=lambda t: F_w,
                               record_internal_forces=True)
        with self.assertLogs("solidum", level=logging.WARNING) as logs:
            solver.solve()
        self.assertEqual(sum("np.ndarray simple" in m for m in logs.output), 1)


class TestEquivalentLoads(unittest.TestCase):

    def test_disposiciones_distintas_se_fusionan(self):
        a = EquivalentLoads.from_dict({1: np.array([1.0, 2.0]), 2: np.array([3.0, 4.0])})
        b = EquivalentLoads.from_dict({2: np.array([1.0, 1.0]), 3: np.array([5.0, 6.0])})
        c = a.plus(b, 2.0).as_dict()
        np.testing.assert_allclose(c[1], [1.0, 2.0])
        np.testing.assert_allclose(c[2], [5.0, 6.0])
        np.testing.assert_allclose(c[3], [10.0, 12.0])
        with self.assertRaisesRegex(ValueError, "tamaños"):
            a.plus(EquivalentLoads.from_dict({1: np.array([1.0, 2.0, 3.0])}))


if __name__ == "__main__":
    unittest.main()
