"""Regresiones del lote de contrato/API de la auditoría 2026-09-22.

- ``internal_forces`` de marcos, armaduras y cables devuelve las fuerzas
  internas de extremo ``F_int − f_eq`` cuando hay carga distribuida
  (peso propio), coherentes con las reacciones.
- Cachés del ``Assembler``: un apoyo añadido tras el primer ``solve`` se
  honra; ``LinearSolver.invalidate_cache`` invalida también el ensamblador.
- ``solidum.run`` no enmascara ``TypeError`` internos ni ejecuta dos veces.
- Modal en modelos diminutos (1 GDL) y con ``n_modes = n_free``.
"""
import unittest

import numpy as np

import solidum
from solidum.core.domain import Domain
from solidum.elements.frame import Frame2DEuler, Frame2DTimoshenko
from solidum.elements.frame3d import Frame3D
from solidum.elements.solid_2d import Quad4
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver, ModalSolver, NonlinearSolver

E, RHO, A, I, L, G = 2.1e11, 7850.0, 1e-3, 8.33e-6, 2.0, 9.81


def _cantilever(cls, **extra):
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0]); n2 = dom.add_node(2, [L, 0.0])
    mat = Elastic1D(E=E, density=RHO)
    dom.add_element(cls(1, [n1, n2], mat, A=A, I=I, **extra))
    for d in ("ux", "uy", "rz"):
        n1.fix_dof(d, 0.0)
    dom.generate_equation_numbers()
    return dom


class TestInternalForcesWithDistributedLoad(unittest.TestCase):

    def _check_frame2d(self, cls, **extra):
        dom = _cantilever(cls, **extra)
        asm = Assembler(dom)
        F = asm.assemble_self_weight([0.0, -G])
        q = RHO * A * G                       # carga por unidad de longitud
        res = solidum.run(dom, assembler=asm, F_applied=F)
        ef = res.element_forces[1]
        # Convención §5: V_i = +qL, V_j = 0; M_i = −qL²/2 (hogging), M_j = 0.
        np.testing.assert_allclose(ef.components["V"], [q * L, 0.0], atol=1e-6 * q * L)
        np.testing.assert_allclose(ef.components["M"], [-q * L * L / 2.0, 0.0], atol=1e-6 * q * L * L)
        # Coherencia con las reacciones en el empotramiento.
        self.assertAlmostEqual(res.reactions_by_node[1]["uy"], q * L, delta=1e-6 * q * L)
        self.assertAlmostEqual(res.reactions_by_node[1]["rz"], q * L * L / 2.0, delta=1e-6 * q * L * L)

    def test_euler(self):
        self._check_frame2d(Frame2DEuler)

    def test_timoshenko(self):
        self._check_frame2d(Frame2DTimoshenko, As=0.8 * A, nu=0.3)

    def test_frame3d(self):
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0, 0.0]); n2 = dom.add_node(2, [L, 0.0, 0.0])
        mat = Elastic1D(E=E, density=RHO)
        dom.add_element(Frame3D(1, [n1, n2], mat, A=A, Iy=I, Iz=I, J=2 * I))
        for d in ("ux", "uy", "uz", "rx", "ry", "rz"):
            n1.fix_dof(d, 0.0)
        dom.generate_equation_numbers()
        asm = Assembler(dom)
        F = asm.assemble_self_weight([0.0, -G, 0.0])
        q = RHO * A * G
        res = solidum.run(dom, assembler=asm, F_applied=F)
        ef = res.element_forces[1]
        # Con el ref_vector por defecto de Frame3D, para un miembro sobre +x
        # global la carga en -y global cae en el eje local z (misma pareja
        # Vz/My que el test de carga en punta de test_frame3d_internal_forces).
        # Nodo i en la cara -x: Vz_i = +qL (como la reaccion), My_i = -qL^2/2.
        np.testing.assert_allclose(ef.components["Vz"], [q * L, 0.0], atol=1e-6 * q * L)
        np.testing.assert_allclose(ef.components["My"], [-q * L * L / 2.0, 0.0], atol=1e-6 * q * L * L)
        np.testing.assert_allclose(ef.components["Vy"], [0.0, 0.0], atol=1e-9 * q * L)
        np.testing.assert_allclose(ef.components["Mz"], [0.0, 0.0], atol=1e-9 * q * L * L)

    def test_truss_hanging_bar(self):
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0]); n2 = dom.add_node(2, [0.0, -L])
        mat = Elastic1D(E=E, density=RHO)
        dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
        n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0); n2.fix_dof("ux", 0.0)
        dom.generate_equation_numbers()
        asm = Assembler(dom)
        F = asm.assemble_self_weight([0.0, -G])
        res = solidum.run(dom, assembler=asm, F_applied=F)
        W = RHO * A * L * G
        # Barra colgada: N_i = W (todo el peso) y N_j = 0 (extremo libre).
        np.testing.assert_allclose(res.element_forces[1].components["N"], [W, 0.0], atol=1e-6 * W)

    def test_no_distributed_load_unchanged(self):
        dom = _cantilever(Frame2DEuler)
        F = np.zeros(dom.total_dofs); F[dom.nodes[2].dofs["uy"]] = -1.0
        res = solidum.run(dom, F_applied=F)
        np.testing.assert_allclose(res.element_forces[1].components["V"], [1.0, 1.0], atol=1e-9)
        np.testing.assert_allclose(res.element_forces[1].components["M"], [-L, 0.0], atol=1e-9)


class TestAssemblerCaches(unittest.TestCase):

    def _plate(self):
        dom = Domain()
        for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
            dom.add_node(i, [float(x), float(y)])
        dom.add_element(Quad4(1, [dom.nodes[i] for i in (1, 2, 3, 4)], Elastic2D(E=1000.0, nu=0.3, density=1.0)))
        dom.generate_equation_numbers()
        for k in (1, 4):
            dom.nodes[k].fix_dof("ux"); dom.nodes[k].fix_dof("uy")
        F = np.zeros(dom.total_dofs)
        F[dom.nodes[2].dofs["ux"]] = 10.0; F[dom.nodes[3].dofs["ux"]] = 10.0
        return dom, F

    def test_new_support_after_first_solve_is_honoured(self):
        dom, F = self._plate()
        asm = Assembler(dom); s = LinearSolver(asm)
        U1 = s.solve(F)
        self.assertNotEqual(U1[dom.nodes[3].dofs["uy"]], 0.0)
        dom.nodes[3].fix_dof("uy", 0.0)
        s.invalidate_cache()
        U2 = s.solve(F)
        self.assertEqual(U2[dom.nodes[3].dofs["uy"]], 0.0)
        # Y también sin invalidar el solver: el Newton lee el constraint_set fresco.
        dom.nodes[2].fix_dof("uy", 0.0)
        U3 = NonlinearSolver(asm, num_steps=1).solve(F)
        self.assertEqual(U3[dom.nodes[2].dofs["uy"]], 0.0)

    def test_new_element_after_first_assembly_is_included(self):
        dom, F = self._plate()
        asm = Assembler(dom)
        asm.assemble_system()
        n_entries = asm._total_entries
        dom.add_node(5, [2.0, 0.0]); dom.add_node(6, [2.0, 1.0])
        dom.add_element(Quad4(2, [dom.nodes[i] for i in (2, 5, 6, 3)], Elastic2D(E=1000.0, nu=0.3, density=1.0)))
        dom.generate_equation_numbers()
        asm.assemble_system()
        self.assertGreater(asm._total_entries, n_entries)
        self.assertEqual(asm.K_global.shape[0], dom.total_dofs)


class TestInvokeSolve(unittest.TestCase):

    def test_internal_type_error_propagates_and_runs_once(self):
        calls = {"n": 0}

        class _Broken:
            PIPELINE_KIND = "static"

            def solve(self, F, step_callback=None):
                calls["n"] += 1
                raise TypeError("fallo interno genuino")

        dom, F = TestAssemblerCaches()._plate()
        with self.assertRaises(TypeError):
            solidum.run(dom, solver=_Broken(), F_applied=F, step_callback=lambda *a: None)
        self.assertEqual(calls["n"], 1)


class TestModalSmallModels(unittest.TestCase):

    def test_single_dof_and_all_modes(self):
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0]); n2 = dom.add_node(2, [1.0, 0.0])
        dom.add_element(Truss2D(1, [n1, n2], Elastic1D(E=25.0, density=2.0), A=1.0))
        n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0); n2.fix_dof("uy", 0.0)
        dom.generate_equation_numbers()
        res = ModalSolver(Assembler(dom), n_modes=1, lumping="lumped").solve()
        self.assertAlmostEqual(res.frequencies_rad[0], 5.0, places=10)
        # M-ortonormalidad del modo denso.
        M = Assembler(dom).assemble_mass_matrix(lumping="lumped")
        phi = res.modes[:, 0]
        self.assertAlmostEqual(float(phi @ (M @ phi)), 1.0, places=10)


if __name__ == "__main__":
    unittest.main()
