"""Verificación de signos §5 en ``internal_forces`` de frames 2D.

Dos casos analíticos elementales:
1. Cantilever con carga en punta hacia abajo → M[i] hogging (negativo),
   V constante positivo, M[j]=0.
2. Cantilever con momento CCW en punta → M constante y positivo (sagging),
   V=0.

Se prueban las tres formulaciones (Euler, Timoshenko, Corot) con mismo setup.
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.node import Node
from fenix.elements.frame import Frame2DEuler, Frame2DTimoshenko, Frame2DEulerCorot


class LinearElastic1D:
    STRAIN_DIM = 1

    def __init__(self, E, nu=0.3):
        self.E = E
        self.nu = nu

    def compute_state(self, strain, state_var=None):
        return self.E * strain, self.E, state_var


def _setup_cantilever(element_cls, L=1.0):
    """Cantilever horizontal: nodo i en origen (fijo), nodo j a L en +x (libre)."""
    ni = Node(1, [0.0, 0.0])
    nj = Node(2, [L, 0.0])
    for n in (ni, nj):
        n.add_dof('ux'); n.add_dof('uy'); n.add_dof('rz')
    ni.dofs = {'ux': 0, 'uy': 1, 'rz': 2}
    nj.dofs = {'ux': 3, 'uy': 4, 'rz': 5}

    mat = LinearElastic1D(E=1.0)
    if element_cls is Frame2DTimoshenko:
        elem = element_cls(1, [ni, nj], mat, A=1.0, I=1.0, As=1e6)  # As grande → Timoshenko ≈ Euler
    else:
        elem = element_cls(1, [ni, nj], mat, A=1.0, I=1.0)
    return elem, ni, nj


def _solve_cantilever(elem, F_ext_j):
    """Resuelve el cantilever aplicando F_ext en los DOFs de j. Devuelve U global."""
    K = elem.compute_global_stiffness()
    # Nodo i fijo: DOFs 0,1,2 restringidos; libres: 3,4,5.
    K_ff = K[np.ix_([3, 4, 5], [3, 4, 5])]
    u_f = np.linalg.solve(K_ff, F_ext_j)
    U = np.zeros(6)
    U[3:6] = u_f
    return U


class TestCantileverTipLoad(unittest.TestCase):
    """Carga F=1 en -y aplicada en el nodo j del cantilever. L=1, EI=1.

    Resultado analítico (convención §5):
        N = [0, 0]  (sin axial)
        V = [+1, +1]  (constante)
        M = [-1, 0]  (hogging en la raíz, cero en la punta)
    """

    def _check(self, elem_cls, tol=1e-9):
        elem, _, _ = _setup_cantilever(elem_cls)
        U = _solve_cantilever(elem, F_ext_j=np.array([0.0, -1.0, 0.0]))
        ef = elem.internal_forces(U)

        self.assertEqual(ef.kind, "frame2d")
        np.testing.assert_allclose(ef.components["N"], [0.0, 0.0], atol=tol)
        np.testing.assert_allclose(ef.components["V"], [1.0, 1.0], atol=tol)
        np.testing.assert_allclose(ef.components["M"], [-1.0, 0.0], atol=tol)

    def test_euler(self):
        self._check(Frame2DEuler)

    def test_timoshenko(self):
        self._check(Frame2DTimoshenko, tol=1e-4)

    def test_corot_small_strain(self):
        # Corot: la cinemática (l-L0)/L0 introduce un falso N = O(u²) si se
        # resuelve linealmente con carga que genera desplazamientos grandes.
        # Se usa carga pequeña para que el segundo orden quede bajo tolerancia.
        elem, _, _ = _setup_cantilever(Frame2DEulerCorot)
        scale = 1e-4
        U = _solve_cantilever(elem, F_ext_j=np.array([0.0, -scale, 0.0]))
        ef = elem.internal_forces(U)
        np.testing.assert_allclose(ef.components["N"], [0.0, 0.0], atol=scale**2 * 10)
        np.testing.assert_allclose(ef.components["V"], [scale, scale], rtol=1e-4)
        np.testing.assert_allclose(ef.components["M"], [-scale, 0.0], atol=scale * 1e-3)


class TestCantileverTipMoment(unittest.TestCase):
    """Momento Mₐ=+1 (CCW, +z) aplicado en el nodo j. L=1, EI=1.

    Resultado analítico (convención §5):
        N = [0, 0]
        V = [0, 0]
        M = [+1, +1]  (constante, sagging positivo)
    """

    def _check(self, elem_cls, tol=1e-9):
        elem, _, _ = _setup_cantilever(elem_cls)
        U = _solve_cantilever(elem, F_ext_j=np.array([0.0, 0.0, 1.0]))
        ef = elem.internal_forces(U)

        np.testing.assert_allclose(ef.components["N"], [0.0, 0.0], atol=tol)
        np.testing.assert_allclose(ef.components["V"], [0.0, 0.0], atol=tol)
        np.testing.assert_allclose(ef.components["M"], [1.0, 1.0], atol=tol)

    def test_euler(self):
        self._check(Frame2DEuler)

    def test_timoshenko(self):
        self._check(Frame2DTimoshenko, tol=1e-4)

    def test_corot_small_strain(self):
        # Ver comentario en TestCantileverTipLoad.test_corot_small_strain.
        elem, _, _ = _setup_cantilever(Frame2DEulerCorot)
        scale = 1e-4
        U = _solve_cantilever(elem, F_ext_j=np.array([0.0, 0.0, scale]))
        ef = elem.internal_forces(U)
        np.testing.assert_allclose(ef.components["N"], [0.0, 0.0], atol=scale**2 * 10)
        np.testing.assert_allclose(ef.components["V"], [0.0, 0.0], atol=scale * 1e-4)
        np.testing.assert_allclose(ef.components["M"], [scale, scale], rtol=1e-4)


class TestAxialTension(unittest.TestCase):
    """Tracción pura F=+1 en +x aplicada en j. Debe dar N=[+1,+1], V=M=0."""

    def _check(self, elem_cls, tol=1e-9):
        elem, _, _ = _setup_cantilever(elem_cls)
        U = _solve_cantilever(elem, F_ext_j=np.array([1.0, 0.0, 0.0]))
        ef = elem.internal_forces(U)
        np.testing.assert_allclose(ef.components["N"], [1.0, 1.0], atol=tol)
        np.testing.assert_allclose(ef.components["V"], [0.0, 0.0], atol=tol)
        np.testing.assert_allclose(ef.components["M"], [0.0, 0.0], atol=tol)

    def test_euler(self):
        self._check(Frame2DEuler)

    def test_timoshenko(self):
        self._check(Frame2DTimoshenko, tol=1e-6)

    def test_corot_small_strain(self):
        self._check(Frame2DEulerCorot, tol=1e-6)


if __name__ == "__main__":
    unittest.main()
