"""Verificación de signos §5 en ``internal_forces`` de Frame3D.

Convención stress-resultant / RHR pura. Cantilever alineado con +x, nodo i
fijo en origen, nodo j libre a L=1. Cuatro cargas unitarias independientes
aplicadas en j; comparación con los valores analíticos derivados de §5.
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.node import Node
from fenix.elements.frame3d import Frame3D


class LinearElastic1D:
    STRAIN_DIM = 1

    def __init__(self, E, nu=0.3, density: float = 0.0):
        self.E = E
        self.nu = nu

    def compute_state(self, strain, state_var=None):
        return self.E * strain, self.E, state_var


def _cantilever():
    ni = Node(1, [0.0, 0.0, 0.0])
    nj = Node(2, [1.0, 0.0, 0.0])
    dof_names = ['ux', 'uy', 'uz', 'rx', 'ry', 'rz']
    for n in (ni, nj):
        for d in dof_names:
            n.add_dof(d)
    ni.dofs = {d: i for i, d in enumerate(dof_names)}
    nj.dofs = {d: i + 6 for i, d in enumerate(dof_names)}
    mat = LinearElastic1D(E=1.0, nu=0.0, density=0.0)  # G = 1/2, fijo
    # ref_vector=(0,1,0) fuerza y_local=global_y, z_local=global_z.
    # Con el default (0,0,1) el frame se construye rotado y las cargas globales
    # no coinciden con las locales, complicando la interpretación del test.
    elem = Frame3D(1, [ni, nj], mat, A=1.0, Iy=1.0, Iz=1.0, J=1.0,
                   ref_vector=[0.0, 1.0, 0.0])
    return elem


def _solve(elem, F_ext_j):
    """Resuelve con nodo i fijo (6 DOFs), nodo j libre; F_ext_j shape (6,)."""
    K = elem.compute_global_stiffness()
    free = list(range(6, 12))
    K_ff = K[np.ix_(free, free)]
    u_f = np.linalg.solve(K_ff, F_ext_j)
    U = np.zeros(12)
    U[6:12] = u_f
    return U


class TestFrame3DSigns(unittest.TestCase):
    def setUp(self):
        self.elem = _cantilever()

    def test_axial_tension(self):
        """F=+1 en +x en j → N=[+1,+1], resto=0."""
        U = _solve(self.elem, np.array([1., 0., 0., 0., 0., 0.]))
        ef = self.elem.internal_forces(U)
        np.testing.assert_allclose(ef.components["N"],  [1., 1.], atol=1e-9)
        np.testing.assert_allclose(ef.components["Vy"], [0., 0.], atol=1e-9)
        np.testing.assert_allclose(ef.components["Vz"], [0., 0.], atol=1e-9)
        np.testing.assert_allclose(ef.components["T"],  [0., 0.], atol=1e-9)
        np.testing.assert_allclose(ef.components["My"], [0., 0.], atol=1e-9)
        np.testing.assert_allclose(ef.components["Mz"], [0., 0.], atol=1e-9)

    def test_tip_load_y(self):
        """F=-1 en -y en j → Vy=[-1,-1], Mz=[-1,0] (hogging §5 A: Mz<0)."""
        U = _solve(self.elem, np.array([0., -1., 0., 0., 0., 0.]))
        ef = self.elem.internal_forces(U)
        np.testing.assert_allclose(ef.components["Vy"], [-1., -1.], atol=1e-9)
        np.testing.assert_allclose(ef.components["Mz"], [-1.,  0.], atol=1e-9)
        np.testing.assert_allclose(ef.components["Vz"], [0., 0.], atol=1e-9)
        np.testing.assert_allclose(ef.components["My"], [0., 0.], atol=1e-9)

    def test_tip_load_z(self):
        """F=-1 en -z en j → Vz=[-1,-1], My=[+1,0]."""
        U = _solve(self.elem, np.array([0., 0., -1., 0., 0., 0.]))
        ef = self.elem.internal_forces(U)
        np.testing.assert_allclose(ef.components["Vz"], [-1., -1.], atol=1e-9)
        np.testing.assert_allclose(ef.components["My"], [ 1.,  0.], atol=1e-9)
        np.testing.assert_allclose(ef.components["Vy"], [0., 0.], atol=1e-9)
        np.testing.assert_allclose(ef.components["Mz"], [0., 0.], atol=1e-9)

    def test_pure_torsion(self):
        """Mx=+1 en +x en j → T=[+1,+1]."""
        U = _solve(self.elem, np.array([0., 0., 0., 1., 0., 0.]))
        ef = self.elem.internal_forces(U)
        np.testing.assert_allclose(ef.components["T"],  [1., 1.], atol=1e-9)
        np.testing.assert_allclose(ef.components["N"],  [0., 0.], atol=1e-9)
        np.testing.assert_allclose(ef.components["Vy"], [0., 0.], atol=1e-9)
        np.testing.assert_allclose(ef.components["Vz"], [0., 0.], atol=1e-9)

    def test_kind(self):
        U = _solve(self.elem, np.array([1., 0., 0., 0., 0., 0.]))
        ef = self.elem.internal_forces(U)
        self.assertEqual(ef.kind, "frame3d")
        self.assertEqual(set(ef.components), {"N", "Vy", "Vz", "T", "My", "Mz"})


if __name__ == "__main__":
    unittest.main()
