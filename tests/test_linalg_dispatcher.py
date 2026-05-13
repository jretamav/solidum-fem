"""Tests de la capa algebraica fenix.math.linalg (ADR 0003 fase 1).

Cubre:
1. Equivalencia bit-a-bit del LUSolver nuevo con el spsolve directo previo.
2. Equivalencia LU vs Cholesky en problema SPD (si scikit-sparse está disponible).
3. Despachador automático: SPD → Cholesky (si disponible), no-SPD → LU.
4. Override desde YAML (linear_algebra) llega al solver y se respeta.
5. Fallback automático SPD→LU cuando Cholesky reporta no-positividad.
"""
import os
import sys
import unittest

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.math.linalg import LUSolver, StiffnessProperties, select_solver
from fenix.math.linalg.dispatcher import _HAS_CHOLESKY


def _build_spd_matrix(n: int = 50, seed: int = 0) -> sp.csr_matrix:
    """Construye una K dispersa SPD: tridiagonal con dominancia diagonal."""
    rng = np.random.default_rng(seed)
    main = 4.0 + rng.uniform(0.0, 1.0, size=n)
    off = -1.0 + rng.uniform(-0.1, 0.1, size=n - 1)
    K = sp.diags([off, main, off], offsets=[-1, 0, 1], format="csr")
    return K


def _build_nonsymmetric_matrix(n: int = 50) -> sp.csr_matrix:
    """K no simétrica: tridiagonal con coeficientes asimétricos."""
    main = np.full(n, 4.0)
    upper = np.full(n - 1, -2.0)
    lower = np.full(n - 1, -1.0)  # ≠ upper → no simétrica
    return sp.diags([lower, main, upper], offsets=[-1, 0, 1], format="csr")


class TestLUSolverEquivalence(unittest.TestCase):
    """LUSolver debe reproducir exactamente el comportamiento previo (spsolve)."""

    def test_lu_matches_raw_spsolve(self):
        K = _build_spd_matrix(n=80)
        b = np.linspace(1.0, 2.0, 80)

        x_old = spla.spsolve(K, b)
        x_new = LUSolver().solve(K, b)

        np.testing.assert_array_equal(x_old, x_new)


@unittest.skipUnless(_HAS_CHOLESKY, "scikit-sparse no instalado; omitiendo Cholesky")
class TestCholeskyVsLU(unittest.TestCase):
    """En problemas SPD, Cholesky y LU deben coincidir a tolerancia de redondeo."""

    def test_solutions_agree(self):
        from fenix.math.linalg import CholeskySolver  # type: ignore[attr-defined]

        K = _build_spd_matrix(n=120, seed=42)
        b = np.cos(np.linspace(0.0, 6.28, 120))

        x_lu = LUSolver().solve(K, b)
        x_chol = CholeskySolver().solve(K, b)

        np.testing.assert_allclose(x_chol, x_lu, rtol=1e-10, atol=1e-12)


class TestDispatcher(unittest.TestCase):
    """select_solver elige el backend correcto según las propiedades de K."""

    def test_dispatcher_picks_cholesky_for_spd(self):
        if not _HAS_CHOLESKY:
            self.skipTest("scikit-sparse no instalado")
        props = StiffnessProperties(is_symmetric=True, is_positive_definite=True, size=100)
        solver = select_solver(props)
        self.assertEqual(solver.name, "cholesky")

    def test_dispatcher_picks_lu_for_nonsymmetric(self):
        props = StiffnessProperties(is_symmetric=False, is_positive_definite=True, size=100)
        solver = select_solver(props)
        self.assertEqual(solver.name, "lu")

    def test_dispatcher_picks_lu_for_indefinite(self):
        props = StiffnessProperties(is_symmetric=True, is_positive_definite=False, size=100)
        solver = select_solver(props)
        # En fase 1 todavía no existe LDLᵀ, así que indefinida-simétrica → LU.
        self.assertEqual(solver.name, "lu")

    def test_override_forces_lu_even_when_spd(self):
        props = StiffnessProperties(is_symmetric=True, is_positive_definite=True, size=100)
        solver = select_solver(props, override="lu")
        self.assertEqual(solver.name, "lu")

    def test_override_auto_is_equivalent_to_none(self):
        props = StiffnessProperties(is_symmetric=False, is_positive_definite=False, size=100)
        a = select_solver(props, override=None)
        b = select_solver(props, override="auto")
        self.assertEqual(type(a), type(b))

    def test_override_unknown_raises(self):
        props = StiffnessProperties(is_symmetric=True, is_positive_definite=True, size=100)
        with self.assertRaises(ValueError):
            select_solver(props, override="pardiso_xyz")


@unittest.skipUnless(_HAS_CHOLESKY, "scikit-sparse no instalado; omitiendo Cholesky")
class TestCholeskyFallback(unittest.TestCase):
    """Cholesky aborta con CholeskyNotPositiveDefiniteError ante K no-SPD."""

    def test_cholesky_aborts_on_indefinite(self):
        from fenix.math.linalg import CholeskySolver  # type: ignore[attr-defined]
        from fenix.math.linalg.cholesky import CholeskyNotPositiveDefiniteError

        # K simétrica indefinida: diagonal con un autovalor negativo claro.
        K = sp.diags([1.0, 2.0, -1.0, 3.0], format="csr")
        b = np.array([1.0, 1.0, 1.0, 1.0])

        with self.assertRaises(CholeskyNotPositiveDefiniteError):
            CholeskySolver().solve(K, b)


class TestLinearSolverYAMLOverride(unittest.TestCase):
    """El campo YAML 'linear_algebra' llega al solver y se respeta."""

    def test_override_propagates_to_linear_solver(self):
        from fenix.elements.truss import Truss2D
        from fenix.core.domain import Domain
        from fenix.core.material import Material
        from fenix.math.assembly import Assembler
        from fenix.math.solvers import LinearSolver

        class _ElasticMat(Material):
            STRAIN_DIM = 1
            def __init__(self, E, density: float = 0.0):
                self.E = E
            def compute_state(self, strain, state_vars=None):
                return self.E * strain, self.E, state_vars

        domain = Domain()
        n1 = domain.add_node(1, [0.0, 0.0])
        n2 = domain.add_node(2, [1.0, 0.0])
        mat = _ElasticMat(E=1000.0)
        domain.add_element(Truss2D(1, [n1, n2], mat, A=1.0))
        n1.fix_dof("ux", 0.0)
        n1.fix_dof("uy", 0.0)
        n2.fix_dof("uy", 0.0)
        domain.generate_equation_numbers()

        F = np.zeros(domain.total_dofs)
        F[n2.dofs["ux"]] = 10.0

        for choice in ("auto", "lu"):
            solver = LinearSolver(Assembler(domain), linear_algebra=choice)
            U = solver.solve(F)
            # u = F·L/(E·A) = 10·1/(1000·1) = 0.01
            self.assertAlmostEqual(U[n2.dofs["ux"]], 0.01, places=10)


if __name__ == "__main__":
    unittest.main()
