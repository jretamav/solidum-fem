"""Tests de los elementos de marco/viga 2D.

Materializan los criterios de `acceptance` declarados en
docs/specs/Frame2DEuler.md y docs/specs/Frame2DTimoshenko.md.
"""
import math
import unittest

import numpy as np

import fenix  # autodiscover
from fenix.core.domain import Domain
from fenix.elements.frame import Frame2DEuler, Frame2DTimoshenko, Frame2DEulerCorot
from fenix.materials.elastic import Elastic1D
from fenix.math.assembly import Assembler
from fenix.math.solvers import NonlinearSolver


class TestFrame2DEulerAcceptance(unittest.TestCase):

    def _cantilever(self, L=2.0, E=210e9, A=1e-3, I=1e-6):
        """Voladizo horizontal: nodo 1 empotrado, nodo 2 libre."""
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [L, 0.0])
        mat = Elastic1D(E=E)
        elem = Frame2DEuler(1, [n1, n2], mat, A=A, I=I)
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n1.fix_dof('rz', 0.0)
        dom.add_element(elem)
        dom.generate_equation_numbers(verbose=False)
        return dom, n1, n2, elem, E, A, I, L

    def test_acceptance_respuesta_axial_pura(self):
        """Criterio 1: carga axial F ⇒ u_x(L) = F·L/(E·A)."""
        dom, n1, n2, elem, E, A, I, L = self._cantilever()
        F = 1000.0
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs['ux']] = F

        solver = NonlinearSolver(Assembler(dom), tol=1e-10, num_steps=1)
        U = solver.solve(F_ext)

        self.assertAlmostEqual(U[n2.dofs['ux']], F * L / (E * A), places=10)
        self.assertAlmostEqual(U[n2.dofs['uy']], 0.0, places=12)
        self.assertAlmostEqual(U[n2.dofs['rz']], 0.0, places=12)

    def test_acceptance_voladizo_carga_transversal(self):
        """Criterio 2: carga transversal P ⇒ v(L) = P·L³/(3EI), θ(L) = P·L²/(2EI)."""
        dom, n1, n2, elem, E, A, I, L = self._cantilever()
        P = 500.0
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs['uy']] = P

        solver = NonlinearSolver(Assembler(dom), tol=1e-10, num_steps=1)
        U = solver.solve(F_ext)

        v_expected = P * L**3 / (3 * E * I)
        theta_expected = P * L**2 / (2 * E * I)
        self.assertAlmostEqual(U[n2.dofs['uy']], v_expected, places=8)
        self.assertAlmostEqual(U[n2.dofs['rz']], theta_expected, places=8)
        self.assertAlmostEqual(U[n2.dofs['ux']], 0.0, places=12)

    def test_acceptance_simetria_K(self):
        """Criterio 3: K_global = K_global.T en cualquier configuración (régimen elástico)."""
        # Frame oblicuo para no caer en casos triviales c=1,s=0 o c=0,s=1
        n1 = fenix.core.node.Node(1, [0.0, 0.0])
        n2 = fenix.core.node.Node(2, [3.0, 4.0])  # L = 5
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy'); node.add_dof('rz')
            node.dofs['ux'] = 3 * k
            node.dofs['uy'] = 3 * k + 1
            node.dofs['rz'] = 3 * k + 2
        mat = Elastic1D(E=210e9)
        elem = Frame2DEuler(1, [n1, n2], mat, A=1e-3, I=1e-6)

        u_e = np.array([0.0, 0.0, 0.0, 1e-4, 2e-4, 3e-4])
        K, _ = elem.compute_element_state(u_e)

        self.assertTrue(np.allclose(K, K.T, atol=1e-6))

    def test_registro_en_registry(self):
        from fenix.registry import ElementRegistry
        self.assertIn('Frame2DEuler', ElementRegistry._items)


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


class TestFrame2DEulerCorotAcceptance(unittest.TestCase):
    """Criterios de aceptación de docs/specs/Frame2DEulerCorot.md."""

    def _build_standalone(self, coords1, coords2, E=210e9, A=1e-3, I=1e-6):
        n1 = fenix.core.node.Node(1, list(coords1))
        n2 = fenix.core.node.Node(2, list(coords2))
        for k, node in enumerate((n1, n2)):
            for i, dof in enumerate(('ux', 'uy', 'rz')):
                node.add_dof(dof)
                node.dofs[dof] = 3 * k + i
        mat = Elastic1D(E=E)
        return Frame2DEulerCorot(1, [n1, n2], mat, A=A, I=I)

    def _cantilever_domain(self, L=2.0, E=210e9, A=1e-3, I=1e-6):
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [L, 0.0])
        mat = Elastic1D(E=E)
        elem = Frame2DEulerCorot(1, [n1, n2], mat, A=A, I=I)
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n1.fix_dof('rz', 0.0)
        dom.add_element(elem)
        dom.generate_equation_numbers(verbose=False)
        return dom, n1, n2, elem

    def test_acceptance_limite_lineal_coincide_con_frame2deuler(self):
        """Criterio 1: cargas pequeñas → solución coincide con Frame2DEuler."""
        L, E, A, I = 2.0, 210e9, 1e-3, 1e-6
        dom, n1, n2, _ = self._cantilever_domain(L=L, E=E, A=A, I=I)
        # P elegido para v/L ~ 1e-4 (régimen claramente lineal pero medible)
        P = 10.0
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs['uy']] = P
        U = NonlinearSolver(Assembler(dom), tol=1e-8, num_steps=1).solve(F_ext)
        v_corot = U[n2.dofs['uy']]
        v_linear = P * L**3 / (3 * E * I)
        self.assertAlmostEqual(v_corot / v_linear, 1.0, delta=1e-3)

    def test_acceptance_invariancia_rotacion_rigida(self):
        """Criterio 2: rotación rígida 45° ⇒ σ = 0, F_int = 0."""
        elem = self._build_standalone([0.0, 0.0], [1.0, 0.0])
        # Rotamos rígidamente: nodo 1 en origen, nodo 2 en (cos45°, sin45°)
        # Y las rotaciones nodales deben ser el mismo ángulo α_e = 45°
        theta = math.pi / 4
        c, s = math.cos(theta), math.sin(theta)
        u_e = np.array([
            0.0, 0.0, theta,           # nodo 1: no se mueve, gira 45°
            c - 1.0, s, theta,         # nodo 2: se mueve al nuevo punto, gira 45°
        ])
        K, F_int = elem.compute_element_state(u_e)
        self.assertAlmostEqual(elem.state.stresses_trial[0], 0.0, places=10)
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-10),
                        f"F_int no nulo bajo rotación rígida: {F_int}")

    def test_acceptance_coherencia_tangente_fuerza_interna(self):
        """Criterio 3: K_T analítica coincide con la derivada numérica de F_int."""
        elem = self._build_standalone([0.0, 0.0], [2.0, 0.0])
        # Configuración no trivial
        u_e = np.array([0.01, -0.005, 0.02, -0.008, 0.015, -0.03])

        K_analytical, F0 = elem.compute_element_state(u_e)

        # Jacobiano numérico por diferencias centrales
        h = 1e-7
        K_numerical = np.zeros((6, 6))
        for j in range(6):
            u_plus = u_e.copy(); u_plus[j] += h
            u_minus = u_e.copy(); u_minus[j] -= h
            _, F_plus = elem.compute_element_state(u_plus)
            _, F_minus = elem.compute_element_state(u_minus)
            K_numerical[:, j] = (F_plus - F_minus) / (2 * h)

        # Tolerancia relativa modesta — la fórmula analítica vs FD casan a ~1e-5
        rel_error = np.linalg.norm(K_analytical - K_numerical) / np.linalg.norm(K_analytical)
        self.assertLess(rel_error, 1e-5,
                        f"K_T analítica inconsistente con F_int: error relativo {rel_error:.2e}")

    def test_acceptance_simetria_K(self):
        """Criterio 4: K_T simétrica en configuración arbitraria."""
        elem = self._build_standalone([0.0, 0.0], [3.0, 4.0])  # L=5, diagonal
        u_e = np.array([0.01, -0.005, 0.02, -0.015, 0.008, 0.012])
        K, _ = elem.compute_element_state(u_e)
        self.assertTrue(np.allclose(K, K.T, atol=1e-6))

    def test_registro_en_registry(self):
        from fenix.registry import ElementRegistry
        self.assertIn('Frame2DEulerCorot', ElementRegistry._items)


if __name__ == '__main__':
    unittest.main()
