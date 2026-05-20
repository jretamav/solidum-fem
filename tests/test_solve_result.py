"""Test end-to-end de build_solve_result (ADR 0002, paso 4).

Monta un cantilever 2D de una sola barra Frame2DEuler, resuelve linealmente
con carga tip, y verifica que SolveResult contiene:
- U con los valores analíticos.
- R con las reacciones correctas en el nodo fijo.
- reactions_by_node consistente con R.
- element_forces con ElementForces en convención §5.
"""

import os
import sys
import unittest

import numpy as np
import scipy.sparse.linalg as spla

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.frame import Frame2DEuler
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver
from solidum.results import build_solve_result


class LinearElastic1D:
    STRAIN_DIM = 1

    def __init__(self, E, nu=0.3, density: float = 0.0):
        self.E = E
        self.nu = nu

    def compute_state(self, strain, state_var=None):
        return self.E * strain, self.E, state_var


class TestBuildSolveResult(unittest.TestCase):
    def setUp(self):
        """Cantilever: nodo 1 empotrado en origen, nodo 2 libre en (1,0).
        Carga tip F=-1 en uy del nodo 2. L=1, EI=1, EA=1.
        """
        self.domain = Domain()
        n1 = self.domain.add_node(1, [0.0, 0.0])
        n2 = self.domain.add_node(2, [1.0, 0.0])

        mat = LinearElastic1D(E=1.0)
        elem = Frame2DEuler(1, [n1, n2], mat, A=1.0, I=1.0)
        self.domain.add_element(elem)

        n1.fix_dof('ux', 0.0)
        n1.fix_dof('uy', 0.0)
        n1.fix_dof('rz', 0.0)

        self.domain.generate_equation_numbers()

        self.assembler = Assembler(self.domain)
        self.F_applied = np.zeros(self.domain.total_dofs)
        self.F_applied[n2.dofs['uy']] = -1.0

        solver = LinearSolver(self.assembler)
        self.U = solver.solve(self.F_applied.copy())

    def test_result_assembled(self):
        res = build_solve_result(self.domain, self.assembler, self.U, self.F_applied)

        # U: desplazamientos no nulos solo en nodo 2.
        self.assertAlmostEqual(self.U[0], 0.0, places=6)  # n1.ux
        self.assertAlmostEqual(self.U[1], 0.0, places=6)  # n1.uy
        self.assertAlmostEqual(self.U[2], 0.0, places=6)  # n1.rz

        # Reacciones: equilibrio de cantilever con carga -1 en y y brazo 1.
        # En nodo 1: Fy_reacción=+1, Mz_reacción=+1 (CCW para equilibrar).
        self.assertAlmostEqual(res.reactions_by_node[1]['uy'], 1.0, places=4)
        self.assertAlmostEqual(res.reactions_by_node[1]['rz'], 1.0, places=4)
        self.assertAlmostEqual(res.reactions_by_node[1]['ux'], 0.0, places=4)

        # R global = ceros en DOFs libres del nodo 2.
        for dof in ('ux', 'uy', 'rz'):
            n2_idx = self.domain.get_node(2).dofs[dof]
            self.assertAlmostEqual(res.R[n2_idx], 0.0, places=6)

        # element_forces: Frame2D en convención §5, cantilever → M hogging.
        ef = res.element_forces[1]
        self.assertEqual(ef.kind, "frame2d")
        np.testing.assert_allclose(ef.components["V"], [1.0, 1.0], atol=1e-4)
        np.testing.assert_allclose(ef.components["M"], [-1.0, 0.0], atol=1e-4)
        np.testing.assert_allclose(ef.components["N"], [0.0, 0.0], atol=1e-6)

        # Defaults de estado.
        self.assertTrue(res.converged)
        self.assertEqual(res.num_steps, 1)

    def test_last_result_attribute_initially_none(self):
        self.assertIsNone(Domain().last_result)


if __name__ == "__main__":
    unittest.main()
