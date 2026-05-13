"""Tests del elemento Frame3D.

Materializan los criterios de `acceptance` declarados en
docs/specs/Frame3D.md.
"""
import unittest

import numpy as np

import fenix  # autodiscover
from fenix.core.domain import Domain
from fenix.elements.frame3d import Frame3D
from fenix.materials.elastic import Elastic1D
from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.solvers import NonlinearSolver


class TestFrame3DAcceptance(unittest.TestCase):

    def _cantilever(self, axis_coord=(1.0, 0.0, 0.0), L=2.0, E=210e9,
                    A=1e-3, Iy=1e-6, Iz=1e-6, J=2e-6, nu=0.3,
                    ref_vector=None):
        """Voladizo con nodo 1 empotrado en origen, nodo 2 en L·axis."""
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0, 0.0])
        n2 = dom.add_node(2, [L * axis_coord[0], L * axis_coord[1], L * axis_coord[2]])
        mat = Elastic1D(E=E)
        elem = Frame3D(1, [n1, n2], mat, A=A, Iy=Iy, Iz=Iz, J=J,
                       nu=nu, ref_vector=ref_vector)
        # Empotramiento completo en nodo 1
        for dof in ('ux', 'uy', 'uz', 'rx', 'ry', 'rz'):
            n1.fix_dof(dof, 0.0)
        dom.add_element(elem)
        dom.generate_equation_numbers(verbose=False)
        return dom, n1, n2, elem

    def test_acceptance_respuesta_axial_pura(self):
        """Criterio 1: carga axial F ⇒ u_x = F·L/(E·A)."""
        L, E, A = 2.0, 210e9, 1e-3
        dom, n1, n2, elem = self._cantilever(L=L, E=E, A=A)
        F = 1000.0
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs['ux']] = F

        U = NonlinearSolver(Assembler(dom), convergence=ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10), num_steps=1).solve(F_ext)

        self.assertAlmostEqual(U[n2.dofs['ux']], F * L / (E * A), places=10)
        # Demás DOFs del nodo libre nulos
        for other in ('uy', 'uz', 'rx', 'ry', 'rz'):
            self.assertAlmostEqual(U[n2.dofs[other]], 0.0, places=10)

    def test_acceptance_flexion_xy_coincide_con_frame2d(self):
        """Criterio 2: flexión en plano xy ⇒ v_y = P·L³/(3·E·Iz)."""
        L, E, Iz = 2.0, 210e9, 1e-6
        dom, n1, n2, elem = self._cantilever(L=L, E=E, Iz=Iz,
                                              ref_vector=[0.0, 0.0, 1.0])
        P = 500.0
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs['uy']] = P

        U = NonlinearSolver(Assembler(dom), convergence=ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10), num_steps=1).solve(F_ext)

        v_expected = P * L**3 / (3 * E * Iz)
        theta_expected = P * L**2 / (2 * E * Iz)
        self.assertAlmostEqual(U[n2.dofs['uy']], v_expected, places=8)
        self.assertAlmostEqual(U[n2.dofs['rz']], theta_expected, places=8)
        # Sin deformación en las demás direcciones
        self.assertAlmostEqual(U[n2.dofs['ux']], 0.0, places=10)
        self.assertAlmostEqual(U[n2.dofs['uz']], 0.0, places=10)
        self.assertAlmostEqual(U[n2.dofs['rx']], 0.0, places=10)
        self.assertAlmostEqual(U[n2.dofs['ry']], 0.0, places=10)

    def test_acceptance_flexion_xz(self):
        """Criterio 3: flexión en plano xz ⇒ v_z = P·L³/(3·E·Iy)."""
        L, E, Iy = 2.0, 210e9, 1e-6
        # ref_vector apuntando a y global para que eje z_local = z global
        dom, n1, n2, elem = self._cantilever(L=L, E=E, Iy=Iy,
                                              ref_vector=[0.0, 1.0, 0.0])
        P = 500.0
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs['uz']] = P

        U = NonlinearSolver(Assembler(dom), convergence=ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10), num_steps=1).solve(F_ext)

        v_expected = P * L**3 / (3 * E * Iy)
        self.assertAlmostEqual(U[n2.dofs['uz']], v_expected, places=8)
        # La rotación asociada es alrededor del eje y global con signo
        # consistente con la mano derecha (|θ_y| = PL²/(2EIy))
        self.assertAlmostEqual(abs(U[n2.dofs['ry']]), P * L**2 / (2 * E * Iy),
                               places=8)
        # Sin deformación en las demás
        self.assertAlmostEqual(U[n2.dofs['ux']], 0.0, places=10)
        self.assertAlmostEqual(U[n2.dofs['uy']], 0.0, places=10)
        self.assertAlmostEqual(U[n2.dofs['rx']], 0.0, places=10)
        self.assertAlmostEqual(U[n2.dofs['rz']], 0.0, places=10)

    def test_acceptance_torsion_pura(self):
        """Criterio 4: momento T ⇒ θ_x = T·L/(G·J)."""
        L, E, J, nu = 2.0, 210e9, 2e-6, 0.3
        dom, n1, n2, elem = self._cantilever(L=L, E=E, J=J, nu=nu)
        T_moment = 100.0
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs['rx']] = T_moment

        U = NonlinearSolver(Assembler(dom), convergence=ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10), num_steps=1).solve(F_ext)

        G = E / (2 * (1 + nu))
        theta_expected = T_moment * L / (G * J)
        self.assertAlmostEqual(U[n2.dofs['rx']], theta_expected, places=10)
        # Demás DOFs nulos
        for other in ('ux', 'uy', 'uz', 'ry', 'rz'):
            self.assertAlmostEqual(U[n2.dofs[other]], 0.0, places=10)

    def test_acceptance_simetria_K(self):
        """Criterio 5: K_global = K_global.T en orientación arbitraria."""
        n1 = fenix.core.node.Node(1, [0.0, 0.0, 0.0])
        n2 = fenix.core.node.Node(2, [3.0, 4.0, 12.0])  # dirección arbitraria
        for k, node in enumerate((n1, n2)):
            for i, dof in enumerate(('ux', 'uy', 'uz', 'rx', 'ry', 'rz')):
                node.add_dof(dof)
                node.dofs[dof] = 6 * k + i
        mat = Elastic1D(E=210e9)
        elem = Frame3D(1, [n1, n2], mat, A=1e-3, Iy=1e-6, Iz=2e-6, J=1.5e-6,
                       nu=0.3, ref_vector=[0.0, 0.0, 1.0])

        u_e = np.arange(12, dtype=float) * 1e-5
        K, _ = elem.compute_element_state(u_e)

        self.assertTrue(np.allclose(K, K.T, atol=1e-6))

    def test_rechazo_ref_vector_paralelo(self):
        """ref_vector paralelo al eje de la barra debe rechazarse."""
        with self.assertRaises(ValueError):
            self._cantilever(axis_coord=(1.0, 0.0, 0.0),
                             ref_vector=[1.0, 0.0, 0.0])

    def test_registro_en_registry(self):
        from fenix.registry import ElementRegistry
        self.assertIn('Frame3D', ElementRegistry._items)


if __name__ == '__main__':
    unittest.main()
