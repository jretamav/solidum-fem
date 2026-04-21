"""Tests del elemento Frame2DTimoshenko.

Materializan los criterios de `acceptance` declarados en
docs/specs/Frame2DTimoshenko.md.
"""
import unittest

import numpy as np

import fenix
from fenix.core.domain import Domain
from fenix.elements.frame_timoshenko import Frame2DTimoshenko
from fenix.materials.elastic import Elastic1D
from fenix.math.assembly import Assembler
from fenix.math.solvers import NonlinearSolver


class TestFrame2DTimoshenkoAcceptance(unittest.TestCase):

    def _cantilever(self, L=1.0, E=210e9, A=1e-3, I=1e-6, As=0.8e-3, nu=0.3):
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [L, 0.0])
        mat = Elastic1D(E=E)
        elem = Frame2DTimoshenko(1, [n1, n2], mat, A=A, I=I, As=As, nu=nu)
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n1.fix_dof('rz', 0.0)
        dom.add_element(elem)
        dom.generate_equation_numbers(verbose=False)
        return dom, n1, n2, elem

    def test_acceptance_convergencia_euler_en_viga_esbelta(self):
        """Criterio 1: viga muy esbelta (L/h grande) ⇒ Φ → 0 ⇒ flecha ≈ PL³/(3EI)."""
        L, E, I = 100.0, 210e9, 1e-6
        # Área y A_s consistentes con una sección muy esbelta
        dom, n1, n2, elem = self._cantilever(L=L, E=E, A=1e-3, I=I, As=0.8e-3)
        P = 10.0
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs['uy']] = P

        solver = NonlinearSolver(Assembler(dom), tol=1e-10, num_steps=1)
        U = solver.solve(F_ext)

        v_euler = P * L**3 / (3 * E * I)
        # Tolerancia holgada: Timoshenko difiere de Euler por el término de cortante
        self.assertAlmostEqual(U[n2.dofs['uy']] / v_euler, 1.0, delta=1e-3)

    def test_acceptance_respuesta_axial_pura(self):
        """Criterio 2: carga axial F ⇒ u_x = F·L/(E·A), sin acoplamiento."""
        L, E, A = 2.0, 210e9, 1e-3
        dom, n1, n2, elem = self._cantilever(L=L, E=E, A=A, I=1e-6, As=0.8e-3)
        F = 1000.0
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs['ux']] = F

        solver = NonlinearSolver(Assembler(dom), tol=1e-10, num_steps=1)
        U = solver.solve(F_ext)

        self.assertAlmostEqual(U[n2.dofs['ux']], F * L / (E * A), places=10)
        self.assertAlmostEqual(U[n2.dofs['uy']], 0.0, places=12)
        self.assertAlmostEqual(U[n2.dofs['rz']], 0.0, places=12)

    def test_acceptance_simetria_K(self):
        """Criterio 3: K_global = K_global.T en viga oblicua."""
        n1 = fenix.core.node.Node(1, [0.0, 0.0])
        n2 = fenix.core.node.Node(2, [3.0, 4.0])  # L=5, c=0.6, s=0.8
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy'); node.add_dof('rz')
            node.dofs['ux'] = 3 * k
            node.dofs['uy'] = 3 * k + 1
            node.dofs['rz'] = 3 * k + 2
        mat = Elastic1D(E=210e9)
        elem = Frame2DTimoshenko(1, [n1, n2], mat, A=1e-3, I=1e-6, As=0.8e-3, nu=0.3)

        u_e = np.array([0.0, 0.0, 0.0, 1e-4, 2e-4, 3e-4])
        K, _ = elem.compute_element_state(u_e)

        self.assertTrue(np.allclose(K, K.T, atol=1e-6))

    def test_registro_en_registry(self):
        from fenix.registry import ElementRegistry
        self.assertIn('Frame2DTimoshenko', ElementRegistry._items)


if __name__ == '__main__':
    unittest.main()
