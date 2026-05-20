"""Tests de la fase 2 del ADR 0003: factorize() reusable + Newton modificado.

Cubre:
1. ``LUSolver.factorize`` retorna un objeto reutilizable equivalente a solve repetido.
2. ``CholeskySolver.factorize`` ídem (si scikit-sparse está disponible).
3. ``CholeskyFactorized.n_negative_pivots == 0`` por construcción.
4. ``LUFactorized.n_negative_pivots is None`` (no expone Sturm).
5. ``LDLTSolver`` está registrado pero el despachador degrada a LU con warning.
6. Override ``linear_algebra: ldlt`` desde YAML degrada a LU con warning.
7. Newton modificado en ``NonlinearSolver`` produce el mismo resultado físico.
"""
import os
import sys
import unittest
import warnings

import numpy as np
import scipy.sparse as sp

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # autodiscover
from solidum.math.linalg import LUSolver, StiffnessProperties, select_solver
from solidum.math.linalg.dispatcher import _HAS_CHOLESKY


def _build_spd_matrix(n: int = 50, seed: int = 0) -> sp.csr_matrix:
    rng = np.random.default_rng(seed)
    main = 4.0 + rng.uniform(0.0, 1.0, size=n)
    off = -1.0 + rng.uniform(-0.1, 0.1, size=n - 1)
    return sp.diags([off, main, off], offsets=[-1, 0, 1], format="csr")


class TestLUFactorized(unittest.TestCase):

    def test_factorize_then_solve_matches_direct_solve(self):
        K = _build_spd_matrix(n=80)
        b1 = np.linspace(1.0, 2.0, 80)
        b2 = np.cos(np.linspace(0, np.pi, 80))

        lu = LUSolver()
        factor = lu.factorize(K)

        x1_direct = lu.solve(K, b1)
        x2_direct = lu.solve(K, b2)
        x1_reused = factor.solve(b1)
        x2_reused = factor.solve(b2)

        np.testing.assert_allclose(x1_reused, x1_direct, rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(x2_reused, x2_direct, rtol=1e-12, atol=1e-14)

    def test_lu_does_not_expose_sturm_count(self):
        # SuperLU usa pivoteo parcial por filas: no preserva la inercia.
        K = _build_spd_matrix(n=10)
        factor = LUSolver().factorize(K)
        self.assertIsNone(factor.n_negative_pivots)


@unittest.skipUnless(_HAS_CHOLESKY, "scikit-sparse no instalado; omitiendo Cholesky")
class TestCholeskyFactorized(unittest.TestCase):

    def test_factorize_then_solve_matches_lu(self):
        from solidum.math.linalg import CholeskySolver  # type: ignore[attr-defined]

        K = _build_spd_matrix(n=80, seed=11)
        b = np.sin(np.linspace(0, np.pi, 80))

        chol_factor = CholeskySolver().factorize(K)
        x_chol = chol_factor.solve(b)
        x_lu = LUSolver().solve(K, b)

        np.testing.assert_allclose(x_chol, x_lu, rtol=1e-10, atol=1e-12)

    def test_cholesky_factor_reports_zero_negative_pivots(self):
        from solidum.math.linalg import CholeskySolver  # type: ignore[attr-defined]
        K = _build_spd_matrix(n=10)
        factor = CholeskySolver().factorize(K)
        self.assertEqual(factor.n_negative_pivots, 0)


class TestLDLTPlaceholder(unittest.TestCase):
    """LDLᵀ está registrado pero degrada a LU con warning (fase 2)."""

    def setUp(self):
        # Reset del flag de warning una-vez para que se vuelva a emitir.
        from solidum.math.linalg.ldlt import LDLTSolver
        LDLTSolver._warning_emitted = False

    def test_override_ldlt_emits_warning_and_returns_lu(self):
        props = StiffnessProperties(is_symmetric=True, is_positive_definite=False, size=20)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            solver = select_solver(props, override="ldlt")
        self.assertEqual(solver.name, "lu")
        self.assertTrue(any("LDLTSolver" in str(w.message) for w in caught),
                        f"Esperaba warning sobre LDLTSolver, recibí: {[str(w.message) for w in caught]}")

    def test_ldlt_in_registry(self):
        from solidum.math.linalg.dispatcher import _REGISTRY
        self.assertIn("ldlt", _REGISTRY)


class TestNewtonModificado(unittest.TestCase):
    """Newton modificado: factoriza una vez por paso y reusa la factorización.

    Debe producir el mismo resultado físico que Newton estándar (puede tomar
    más iteraciones en pasos donde K_t cambia mucho, pero converge al mismo
    equilibrio).
    """

    def _cantilever_axial(self, F: float = 1000.0):
        from solidum.core.domain import Domain
        from solidum.elements.frame import Frame2DEuler
        from solidum.materials.elastic import Elastic1D

        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [2.0, 0.0])
        mat = Elastic1D(E=210e9)
        dom.add_element(Frame2DEuler(1, [n1, n2], mat, A=1e-3, I=1e-6))
        n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0); n1.fix_dof("rz", 0.0)
        dom.generate_equation_numbers(verbose=False)
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs["ux"]] = F
        return dom, n2, F_ext

    def test_modified_newton_matches_standard(self):
        from solidum.math.assembly import Assembler
        from solidum.math.convergence import ConvergenceCriterion
        from solidum.math.solvers import NonlinearSolver

        conv = lambda: ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10)
        dom_a, n2_a, F_a = self._cantilever_axial()
        std = NonlinearSolver(Assembler(dom_a), convergence=conv(), num_steps=2)
        U_std = std.solve(F_a)

        dom_b, n2_b, F_b = self._cantilever_axial()
        mod = NonlinearSolver(Assembler(dom_b), convergence=conv(), num_steps=2,
                              freeze_tangent_after_iter=1)
        U_mod = mod.solve(F_b)

        np.testing.assert_allclose(U_std, U_mod, rtol=1e-8, atol=1e-12)

    def test_modified_newton_invalidates_factor_between_steps(self):
        """Tras solve(), la factorización congelada debe haberse liberado."""
        from solidum.math.assembly import Assembler
        from solidum.math.convergence import ConvergenceCriterion
        from solidum.math.solvers import NonlinearSolver

        dom, n2, F = self._cantilever_axial()
        solver = NonlinearSolver(Assembler(dom),
                                 convergence=ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10),
                                 num_steps=2,
                                 freeze_tangent_after_iter=1)
        solver.solve(F)
        self.assertIsNone(solver._frozen_factor)


class TestSPDFallback(unittest.TestCase):
    """Cobertura del fallback automático SPD→LU sin depender de scikit-sparse.

    Inyecta un backend que siempre lanza CholeskyNotPositiveDefiniteError, y
    verifica que el solver lo captura y resuelve correctamente con LU.
    """

    def _cantilever_axial(self):
        from solidum.core.domain import Domain
        from solidum.elements.frame import Frame2DEuler
        from solidum.materials.elastic import Elastic1D

        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [2.0, 0.0])
        mat = Elastic1D(E=210e9)
        dom.add_element(Frame2DEuler(1, [n1, n2], mat, A=1e-3, I=1e-6))
        n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0); n1.fix_dof("rz", 0.0)
        dom.generate_equation_numbers(verbose=False)
        F_ext = np.zeros(dom.total_dofs)
        F_ext[n2.dofs["ux"]] = 1000.0
        return dom, n2, F_ext

    def _make_failing_cholesky(self):
        """Backend simulado que aborta como si CHOLMOD detectara no-positividad."""
        from solidum.math import solvers as solvers_module
        from solidum.math.linalg.lu import LUSolver

        class _FailingCholesky:
            name = "cholesky"
            def solve(self, K, b):
                raise solvers_module.CholeskyNotPositiveDefiniteError("simulated")
            def factorize(self, K):
                raise solvers_module.CholeskyNotPositiveDefiniteError("simulated")

        return _FailingCholesky, LUSolver

    def test_linear_solver_falls_back_to_lu(self):
        from solidum.math import solvers as solvers_module
        from solidum.math.solvers import linear as linear_module
        from solidum.math.assembly import Assembler

        FailingCholesky, _ = self._make_failing_cholesky()

        # ``LinearSolver`` resuelve ``select_solver`` por su import propio en
        # ``solidum.math.solvers.linear``; el monkey-patch debe apuntar a ese
        # submódulo, no al paquete.
        original = linear_module.select_solver
        linear_module.select_solver = lambda props, override=None: FailingCholesky()
        try:
            dom, n2, F = self._cantilever_axial()
            U = solvers_module.LinearSolver(Assembler(dom)).solve(F)
            expected = 1000.0 * 2.0 / (210e9 * 1e-3)
            # Tras eliminación directa (ADR 0004) la imposición es exacta a redondeo.
            np.testing.assert_allclose(U[n2.dofs["ux"]], expected, rtol=1e-12)
        finally:
            linear_module.select_solver = original

    def test_nonlinear_solver_falls_back_to_lu(self):
        from solidum.math import solvers as solvers_module
        from solidum.math.solvers import nonlinear as nonlinear_module
        from solidum.math.assembly import Assembler
        from solidum.math.linalg.lu import LUSolver as RealLU

        FailingCholesky, _ = self._make_failing_cholesky()

        # La primera llamada al despachador devuelve el cholesky simulado;
        # tras el fallback (is_pd=False) la segunda devuelve LU real.
        call_count = {"n": 0}
        def _patched(props, override=None):
            call_count["n"] += 1
            if call_count["n"] == 1:
                return FailingCholesky()
            return RealLU()

        # ``NonlinearSolver`` lo importa en su propio submódulo (ver
        # ``test_linear_solver_falls_back_to_lu`` para la lógica).
        original = nonlinear_module.select_solver
        nonlinear_module.select_solver = _patched
        try:
            dom, n2, F = self._cantilever_axial()
            from solidum.math.convergence import ConvergenceCriterion
            conv = ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10)
            U = solvers_module.NonlinearSolver(Assembler(dom), convergence=conv, num_steps=1).solve(F)
            expected = 1000.0 * 2.0 / (210e9 * 1e-3)
            # Tras eliminación directa (ADR 0004) la imposición es exacta a redondeo.
            np.testing.assert_allclose(U[n2.dofs["ux"]], expected, rtol=1e-12)
        finally:
            nonlinear_module.select_solver = original


if __name__ == "__main__":
    unittest.main()
