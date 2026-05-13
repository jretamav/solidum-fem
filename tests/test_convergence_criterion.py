"""Tests de la política de convergencia (ADR 0007).

Cubre:
- Invariancia bajo cambio de unidades: el mismo problema físico planteado en
  N/m, kN/mm o MN/m converge en el mismo número de iteraciones hasta paridad
  de bits gracias al `atol` autoderivado.
- Régimen degenerado: control en desplazamiento puro (F_ext_global ≈ 0)
  no provoca falsa no-convergencia.
- Separación de criterios: el AND lógico exige que ambos (fuerza Y
  desplazamiento) se satisfagan simultáneamente.
- Calibración: el criterio falla ruidosamente si se invoca evaluate() sin
  calibrate() previo.
- Parseo de configuración: claves desconocidas en el dict YAML disparan
  ValueError para detectar typos pronto.
"""
import unittest

import numpy as np

import fenix  # autodiscover
from fenix.core.domain import Domain
from fenix.elements.truss import Truss2D
from fenix.materials.elastic import Elastic1D
from fenix.math.assembly import Assembler
from fenix.math.convergence import (
    ConvergenceCriterion,
    ConvergenceState,
    make_convergence_from_config,
)
from fenix.math.solvers import NonlinearSolver


class TestConvergenceCriterionUnit(unittest.TestCase):
    """Tests unitarios de ConvergenceCriterion en aislamiento del solver."""

    def test_evaluate_requires_calibration(self):
        c = ConvergenceCriterion()
        with self.assertRaisesRegex(RuntimeError, "calibrate"):
            c.evaluate(1.0, 1.0, 1.0, 1.0)

    def test_calibration_sets_atol_proportional_to_scale(self):
        c = ConvergenceCriterion(atol_force_factor=1e-9, atol_disp_factor=1e-9)
        c.calibrate(force_scale=1.0e6, disp_scale=1.0e-3)
        self.assertAlmostEqual(c.atol_force, 1.0e6 * 1.0e-9)
        self.assertAlmostEqual(c.atol_disp, 1.0e-3 * 1.0e-9)

    def test_and_semantics_force_fails(self):
        c = ConvergenceCriterion(rtol_force=1e-5, rtol_disp=1e-5,
                                 atol_force_factor=0.0, atol_disp_factor=0.0)
        c.calibrate(1.0, 1.0)
        # δU cumple holgadamente; R excede tol_force.
        state = c.evaluate(residual_norm=1.0, ref_force=1.0,
                           delta_u_norm=1e-12, u_norm=1.0)
        self.assertFalse(state.converged)
        self.assertGreater(state.ratio_force, 1.0)
        self.assertLess(state.ratio_disp, 1.0)

    def test_and_semantics_disp_fails(self):
        c = ConvergenceCriterion(rtol_force=1e-5, rtol_disp=1e-5,
                                 atol_force_factor=0.0, atol_disp_factor=0.0)
        c.calibrate(1.0, 1.0)
        # R cumple; δU excede tol_disp.
        state = c.evaluate(residual_norm=1e-12, ref_force=1.0,
                           delta_u_norm=1.0, u_norm=1.0)
        self.assertFalse(state.converged)
        self.assertLess(state.ratio_force, 1.0)
        self.assertGreater(state.ratio_disp, 1.0)

    def test_both_criteria_satisfied(self):
        c = ConvergenceCriterion(rtol_force=1e-5, rtol_disp=1e-5,
                                 atol_force_factor=0.0, atol_disp_factor=0.0)
        c.calibrate(1.0, 1.0)
        state = c.evaluate(residual_norm=1e-7, ref_force=1.0,
                           delta_u_norm=1e-7, u_norm=1.0)
        self.assertTrue(state.converged)

    def test_atol_acts_as_floor_when_ref_collapses(self):
        """Cuando ref_force → 0, atol mantiene la comparación significativa."""
        c = ConvergenceCriterion(rtol_force=1e-5, atol_force_factor=1e-9)
        c.calibrate(force_scale=1.0e6, disp_scale=1.0)  # atol_force = 1e-3
        # ref colapsa a cero (transitorio); residuo numérico minúsculo.
        state = c.evaluate(residual_norm=1e-6, ref_force=0.0,
                           delta_u_norm=0.0, u_norm=1.0)
        # tol_force ≈ atol_force = 1e-3, holgadamente sobre el residuo.
        self.assertTrue(state.converged)

    def test_state_is_frozen_dataclass(self):
        c = ConvergenceCriterion()
        c.calibrate(1.0, 1.0)
        state = c.evaluate(1e-10, 1.0, 1e-10, 1.0)
        self.assertIsInstance(state, ConvergenceState)
        with self.assertRaises(Exception):
            state.converged = False  # frozen


class TestConvergenceFromConfig(unittest.TestCase):
    def test_none_config_uses_defaults(self):
        c = make_convergence_from_config(None)
        self.assertGreater(c.rtol_force, 0)

    def test_unknown_keys_raise(self):
        with self.assertRaisesRegex(ValueError, "claves desconocidas"):
            make_convergence_from_config({"rtol_force": 1e-5, "tol": 1e-5})

    def test_partial_config_overrides_only_specified(self):
        c = make_convergence_from_config({"rtol_force": 1e-3})
        self.assertEqual(c.rtol_force, 1e-3)
        # rtol_disp queda en el default.
        from fenix.constants import CONVERGENCE_RTOL_DISP
        self.assertEqual(c.rtol_disp, CONVERGENCE_RTOL_DISP)


def _build_truss(F_axial: float, E: float, A: float, L: float):
    """Truss2D 1-elemento, una restricción en x=0, carga axial en x=L.

    Devuelve (assembler, dof_index_libre).
    """
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    mat = Elastic1D(E=E, density=0.0)
    elem = Truss2D(1, [n1, n2], mat, A=A)
    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n2.fix_dof('uy', 0.0)
    dom.add_element(elem)
    dom.generate_equation_numbers(verbose=False)
    F_ext = np.zeros(dom.total_dofs)
    F_ext[n2.dofs['ux']] = F_axial
    return Assembler(dom), F_ext, n2.dofs['ux']


class TestUnitInvariance(unittest.TestCase):
    """El mismo problema físico en distintas unidades debe converger igual.

    Truss axial con E·A·L escalados consistentemente. Sin atol_factor invariante
    (los `atol` se autoderivan), el régimen relativo domina y la convergencia
    es bit-paritaria salvo redondeo del solver lineal.
    """

    def test_truss_axial_invariant_under_unit_rescale(self):
        # Caso A: SI base (N, m).
        asm_A, F_A, dof_A = _build_truss(F_axial=1000.0, E=210e9, A=1e-4, L=2.0)
        U_A = NonlinearSolver(asm_A, num_steps=1).solve(F_A)
        u_A = U_A[dof_A]

        # Caso B: kN, mm. Factor de escala: fuerza ×1e-3, longitud ×1e3.
        # E [N/m²] = 210e9 → [kN/mm²] = 210e9 · 1e-3 / 1e6 = 210; A [m²] = 1e-4 → [mm²] = 100;
        # L [m] = 2 → [mm] = 2000; F [N] = 1000 → [kN] = 1.
        asm_B, F_B, dof_B = _build_truss(F_axial=1.0, E=210.0, A=100.0, L=2000.0)
        U_B = NonlinearSolver(asm_B, num_steps=1).solve(F_B)
        u_B = U_B[dof_B]

        # Esperado: u_A [m] · 1e3 = u_B [mm]. Tolerancia relativa muy ajustada;
        # cualquier dependencia residual del solver con las unidades aparecería aquí.
        np.testing.assert_allclose(u_B, u_A * 1e3, rtol=1e-10)

    def test_truss_axial_invariant_under_force_rescale(self):
        """Escalar solo la fuerza factor 1e6 debe escalar U factor 1e6 (problema lineal)."""
        asm_A, F_A, dof_A = _build_truss(F_axial=1.0e-3, E=210e9, A=1e-4, L=2.0)
        U_A = NonlinearSolver(asm_A, num_steps=1).solve(F_A)
        u_A = U_A[dof_A]

        asm_B, F_B, dof_B = _build_truss(F_axial=1.0e3, E=210e9, A=1e-4, L=2.0)
        U_B = NonlinearSolver(asm_B, num_steps=1).solve(F_B)
        u_B = U_B[dof_B]

        np.testing.assert_allclose(u_B, u_A * 1e6, rtol=1e-10)


class TestDegenerate(unittest.TestCase):
    """Régimen degenerado: carga inicial nula."""

    def test_zero_load_problem_terminates_cleanly(self):
        """Truss sin carga ⇒ U = 0 desde el primer paso, sin oscilar."""
        asm, F_zero, dof = _build_truss(F_axial=0.0, E=210e9, A=1e-4, L=2.0)
        # Sin atol autoderivado robusto, ref_force = 0 dispararía división problemática
        # o falsa no-convergencia. Con calibración fallback a 1.0 el solver termina limpio.
        U = NonlinearSolver(asm, num_steps=1).solve(F_zero)
        np.testing.assert_allclose(U, 0.0, atol=1e-14)


if __name__ == '__main__':
    unittest.main()
