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
import dataclasses
import unittest

import numpy as np

import solidum  # autodiscover
from solidum.core.domain import Domain
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.convergence import (
    ConvergenceCriterion,
    ConvergenceState,
    make_convergence_from_config,
)
from solidum.math.solvers import NonlinearSolver


class TestConvergenceCriterionUnit(unittest.TestCase):
    """Tests unitarios de ConvergenceCriterion en aislamiento del solver."""

    def test_evaluate_requires_calibration(self):
        c = ConvergenceCriterion()
        with self.assertRaisesRegex(RuntimeError, "calibrate"):
            c.evaluate(1.0, 1.0, 1.0, 1.0, roundoff_force=0.0, roundoff_anchor=1.0)

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
                           delta_u_norm=1e-12, u_norm=1.0, roundoff_force=0.0, roundoff_anchor=1.0)
        self.assertFalse(state.converged)
        self.assertGreater(state.ratio_force, 1.0)
        self.assertLess(state.ratio_disp, 1.0)

    def test_and_semantics_disp_fails(self):
        c = ConvergenceCriterion(rtol_force=1e-5, rtol_disp=1e-5,
                                 atol_force_factor=0.0, atol_disp_factor=0.0)
        c.calibrate(1.0, 1.0)
        # R cumple; δU excede tol_disp.
        state = c.evaluate(residual_norm=1e-12, ref_force=1.0,
                           delta_u_norm=1.0, u_norm=1.0, roundoff_force=0.0, roundoff_anchor=1.0)
        self.assertFalse(state.converged)
        self.assertLess(state.ratio_force, 1.0)
        self.assertGreater(state.ratio_disp, 1.0)

    def test_both_criteria_satisfied(self):
        c = ConvergenceCriterion(rtol_force=1e-5, rtol_disp=1e-5,
                                 atol_force_factor=0.0, atol_disp_factor=0.0)
        c.calibrate(1.0, 1.0)
        state = c.evaluate(residual_norm=1e-7, ref_force=1.0,
                           delta_u_norm=1e-7, u_norm=1.0, roundoff_force=0.0, roundoff_anchor=1.0)
        self.assertTrue(state.converged)

    def test_atol_acts_as_floor_when_ref_collapses(self):
        """Cuando ref_force → 0, atol mantiene la comparación significativa."""
        c = ConvergenceCriterion(rtol_force=1e-5, atol_force_factor=1e-9)
        c.calibrate(force_scale=1.0e6, disp_scale=1.0)  # atol_force = 1e-3
        # ref colapsa a cero (transitorio); residuo numérico minúsculo.
        state = c.evaluate(residual_norm=1e-6, ref_force=0.0,
                           delta_u_norm=0.0, u_norm=1.0, roundoff_force=0.0, roundoff_anchor=1.0)
        # tol_force ≈ atol_force = 1e-3, holgadamente sobre el residuo.
        self.assertTrue(state.converged)

    def test_state_is_frozen_dataclass(self):
        c = ConvergenceCriterion()
        c.calibrate(1.0, 1.0)
        state = c.evaluate(1e-10, 1.0, 1e-10, 1.0, roundoff_force=0.0, roundoff_anchor=1.0)
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
        from solidum.constants import CONVERGENCE_RTOL_DISP
        self.assertEqual(c.rtol_disp, CONVERGENCE_RTOL_DISP)


def _build_truss(F_axial: float, E: float, A: float, L: float):
    """Truss2D 1-elemento, una restricción en x=0, carga axial en x=L.

    Devuelve (assembler, dof_index_libre).
    """
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    mat = Elastic1D(E=E)
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
        # Sin fuerzas ni valores impuestos las escalas de calibración son nulas
        # (enmienda del 2026-11-27): todo vale cero y el criterio lo acepta
        # (0 ≤ 0) sin divisiones por cero.
        U = NonlinearSolver(asm, num_steps=1).solve(F_zero)
        np.testing.assert_allclose(U, 0.0, atol=1e-14)


# ---------------------------------------------------------------------------
# Piso de redondeo y calibración sin fuerzas (ADR 0007, enmienda 2026-11-27)
# ---------------------------------------------------------------------------

class TestRoundoffFloorUnit(unittest.TestCase):
    """El piso de redondeo y las escalas nulas, en aislamiento del solver."""

    def test_roundoff_floor_adds_factor_times_scale_to_force_tolerance(self):
        c = ConvergenceCriterion(rtol_force=1e-5, atol_force_factor=0.0, roundoff_factor=7.0)
        c.calibrate(0.0, 0.0)
        state = c.evaluate(residual_norm=0.0, ref_force=2.0, delta_u_norm=0.0,
                           u_norm=1.0, roundoff_force=3.0e-8, roundoff_anchor=1.0)
        self.assertEqual(state.tol_force, 1e-5 * 2.0 + 7.0 * 3.0e-8)
        self.assertEqual(state.floor_force, 7.0 * 3.0e-8)

    def test_zero_scales_are_valid_and_an_all_zero_state_converges(self):
        c = ConvergenceCriterion()
        c.calibrate(0.0, 0.0)
        state = c.evaluate(0.0, 0.0, 0.0, 0.0, roundoff_force=0.0, roundoff_anchor=1.0)
        self.assertTrue(state.converged)
        self.assertEqual((state.ratio_force, state.ratio_disp), (0.0, 0.0))

    def test_nonzero_residual_with_zero_tolerance_does_not_converge(self):
        c = ConvergenceCriterion()
        c.calibrate(0.0, 0.0)
        state = c.evaluate(1e-30, 0.0, 0.0, 0.0, roundoff_force=0.0, roundoff_anchor=1.0)
        self.assertFalse(state.converged)
        self.assertEqual(state.ratio_force, float("inf"))

    def test_factors_that_are_not_finite_nonnegative_numbers_raise_naming_the_key(self):
        for key, value in (("roundoff_factor", -1.0), ("rtol_force", float("nan")),
                           ("atol_disp_factor", float("inf")), ("roundoff_factor", "diez"),
                           ("rtol_disp", None), ("atol_force_factor", True)):
            with self.subTest(key=key, value=value), self.assertRaisesRegex(ValueError, key):
                ConvergenceCriterion(**{key: value})

    def test_roundoff_factor_is_read_from_yaml_config(self):
        c = make_convergence_from_config({"roundoff_factor": 100.0})
        self.assertEqual(c.roundoff_factor, 100.0)

    def test_roundoff_scale_is_eps_times_norm_of_abs_operator_times_abs_vector(self):
        """``ε·‖|T|ᵀ·Σ|A|·|x|‖``: valores absolutos entrada a entrada (sin
        cancelaciones), proyectado con ``|T|`` a las componentes de ``TᵀR``."""
        import scipy.sparse as sp
        from solidum.math.convergence import roundoff_scale
        eps = np.finfo(float).eps
        K = sp.csr_matrix(np.array([[2.0, -1.0, 0.0], [-1.0, 4.0, -3.0], [0.0, -3.0, 5.0]]))
        M = sp.diags([1.0, 2.0, 1.0]).tocsr()
        x = np.array([1.0, -2.0, 0.5])
        absKx = np.array([2.0 + 2.0, 1.0 + 8.0 + 1.5, 6.0 + 2.5])
        self.assertAlmostEqual(roundoff_scale((K, x)) / (eps * np.linalg.norm(absKx)), 1.0,
                               places=14)
        # Varios términos se suman componente a componente; un vector nulo no aporta.
        self.assertEqual(roundoff_scale((K, x), (M, np.zeros(3))), roundoff_scale((K, x)))
        a = np.array([0.0, 1.0, 0.0])
        self.assertAlmostEqual(
            roundoff_scale((K, x), (M, a)) / (eps * np.linalg.norm(absKx + [0.0, 2.0, 0.0])),
            1.0, places=14)
        # Una matriz densa da lo mismo que la dispersa.
        self.assertAlmostEqual(roundoff_scale((K.toarray(), x)) / roundoff_scale((K, x)), 1.0,
                               places=14)
        # |T|ᵀ: el segundo grado de libertad prescrito y el tercero esclavo del
        # primero con coeficiente −2 (u₃ = −2·u₁).
        T = sp.csr_matrix(np.array([[1.0], [0.0], [-2.0]]))
        self.assertAlmostEqual(roundoff_scale((K, x), T=T) / (eps * (4.0 + 2.0 * 8.5)), 1.0,
                               places=14)
        self.assertEqual(roundoff_scale((K, np.zeros(3))), 0.0)

    def test_roundoff_scale_has_the_units_of_each_residual_component(self):
        """Con giros: escalar las longitudes por λ multiplica las filas de
        fuerza de ``|K||U|`` por 1 y las de momento por λ, igual que el
        residuo; una escala ``max|diag K|·‖U‖`` no lo hace."""
        import scipy.sparse as sp
        from solidum.math.convergence import roundoff_scale
        # Resorte de traslación k (F/L) y de giro c (F·L): u en L, θ adimensional.
        k, c, u, th, lam = 3.0e6, 2.0e3, 1.0e-3, 2.0e-3, 1.0e3
        K_m = sp.diags([k, c]).tocsr()
        K_mm = sp.diags([k / lam, c * lam]).tocsr()
        f_m = roundoff_scale((K_m, np.array([u, th])))
        f_mm = roundoff_scale((K_mm, np.array([u * lam, th])))
        eps = np.finfo(float).eps
        self.assertAlmostEqual(f_m / (eps * np.hypot(k * u, c * th)), 1.0, places=14)
        self.assertAlmostEqual(f_mm / (eps * np.hypot(k * u, c * th * lam)), 1.0, places=14)

    def test_each_newton_problem_floors_with_the_operators_of_its_residual(self):
        """Estática y arc-length: ``(K, U)``; Newmark/HHT: ``(K, u)``,
        ``(C, u̇)`` y ``(M, ü)``, los tres términos de su residuo."""
        import scipy.sparse as sp
        from solidum.math.solvers.arclength import _ArcProblem
        from solidum.math.solvers.newmark import _DynamicNewtonProblem
        from solidum.math.solvers.nonlinear import _IncrementalProblem
        eps = np.finfo(float).eps
        K = sp.diags([3.0, 7.0]).tocsr()
        C = sp.diags([0.5, 2.0]).tocsr()
        M = sp.diags([11.0, 1.0]).tocsr()
        u, v, a = np.array([3.0, 4.0]), np.array([0.0, 2.0]), np.array([1.0, 0.0])
        state = (K, np.zeros(2))
        # El primer grado de libertad prescrito: TᵀR es la segunda componente.
        T = sp.csr_matrix(np.array([[0.0], [1.0]]))
        # Sólo se usan los operadores, T y el iterado: se evita el constructor.
        static = object.__new__(_IncrementalProblem)
        static._T = T
        self.assertEqual(static.roundoff_force(u, state), eps * 7.0 * 4.0)
        arc = object.__new__(_ArcProblem)
        arc._T = T
        self.assertEqual(arc.roundoff_force((u, 0.3, np.zeros(2)), state), eps * 7.0 * 4.0)
        dyn = object.__new__(_DynamicNewtonProblem)
        dyn.M, dyn.C, dyn.T_op = M, C, T
        self.assertEqual(dyn.roundoff_force((u, v, a, None), state),
                         eps * (7.0 * 4.0 + 2.0 * 2.0 + 1.0 * 0.0))
        dyn.T_op = sp.identity(2, format="csr")
        self.assertAlmostEqual(
            dyn.roundoff_force((u, v, a, None), state)
            / (eps * np.hypot(3.0 * 3.0 + 0.0 + 11.0 * 1.0, 7.0 * 4.0 + 2.0 * 2.0 + 0.0)),
            1.0, places=14)

    def test_an_incomplete_newton_problem_is_rejected_before_assembling(self):
        from solidum.math.solvers.corrector import NewtonCorrector

        class OldProblem:
            def assemble(self, x):
                raise AssertionError("no debe ensamblar")
            residual = residual_norm = calibration_scales = reference_force = assemble
            x_norm = correction = apply = on_converged = assemble

        corrector = NewtonCorrector(ConvergenceCriterion(), max_iter=5, is_symmetric=True)
        with self.assertRaisesRegex(TypeError, r"roundoff_force.*ADR 0007"):
            corrector.run(OldProblem(), np.zeros(2))

    def test_calibration_without_forces_has_zero_force_scale_and_imposed_disp_scale(self):
        import scipy.sparse as sp
        from solidum.math.solvers.corrector import default_calibration_scales

        class _Imposed:
            def __init__(self, g):
                self.g = np.asarray(g, dtype=float)

            def imposed_values(self):
                return self.g

        K = sp.diags([4.0, 8.0, 2.0]).tocsr()
        zero = np.zeros(3)
        self.assertEqual(default_calibration_scales(zero, zero, K, _Imposed([0.0, 3.0, 4.0])),
                         (0.0, 5.0))
        # Con fuerzas, la calibración de siempre (ADR 0007).
        F = np.array([0.0, 6.0, 8.0])
        self.assertEqual(default_calibration_scales(F, zero, K, _Imposed(zero)), (10.0, 10.0 / 8.0))

    def test_assembler_imposed_values_are_the_prescribed_displacements(self):
        asm = _gap_block(3.0e10)
        g = asm.imposed_values()
        d = asm.domain
        self.assertEqual(g[d.nodes[5].dofs["uy"]], -2.0e-3)
        self.assertAlmostEqual(np.linalg.norm(g), 2.0e-3 * np.sqrt(2.0), places=15)


class _RecordingCriterion(ConvergenceCriterion):
    """Criterio que guarda cada evaluación ``(estado, roundoff_force)``."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.records: list[tuple[ConvergenceState, float]] = []

    def evaluate(self, *args, **kwargs):
        state = super().evaluate(*args, **kwargs)
        self.records.append((state, kwargs["roundoff_force"]))
        return state


class _ExtraIterations(ConvergenceCriterion):
    """Criterio que, al cumplirse, obliga a iterar hasta la 5.ª evaluación del
    paso: las evaluaciones forzadas miden el residuo que deja el redondeo.
    Guarda ``(estado, roundoff_force, forzada)``."""

    def __init__(self):
        super().__init__()
        self.records: list[tuple[ConvergenceState, float, bool]] = []
        self._n = 0

    def evaluate(self, *args, **kwargs):
        state = super().evaluate(*args, **kwargs)
        self._n += 1
        forced = state.converged and self._n % 5 != 0
        self.records.append((state, kwargs["roundoff_force"], forced or
                             (state.converged and self._n % 5 == 0)))
        if forced:
            return dataclasses.replace(state, converged=False)
        if state.converged:
            self._n = 0
        return state


def _gap_block(E: float, h: float = 0.1, g0: float = 1.0e-3, t: float = 1.0, *,
               density: float | None = None, impose: bool = True):
    """Bloque elástico de h × h sobre una base fija a través de una interfaz
    de contacto con holgura ``g0``; con ``impose`` la cara superior baja
    ``2·g0`` impuesta. Hasta cerrar la holgura el bloque se mueve sin ninguna
    fuerza. ``E`` en la unidad de fuerza elegida por m² (longitudes en m)."""
    from solidum import ContactPenalty, Elastic2D, Interface2D4, Quad4
    d = Domain()
    for nid, x in ((1, [0.0, 0.0]), (2, [h, 0.0]), (3, [h, g0]), (4, [0.0, g0]),
                   (5, [h, g0 + h]), (6, [0.0, g0 + h])):
        d.add_node(nid, x)
    n = d.nodes
    d.add_element(Interface2D4(1, [n[1], n[2], n[3], n[4]], ContactPenalty(K_n=10 * E / h),
                               thickness=t))
    d.add_element(Quad4(2, [n[4], n[3], n[5], n[6]], Elastic2D(E=E, nu=0.2, density=density),
                        thickness=t))
    for nid in (1, 2):
        n[nid].fix_dof("ux"); n[nid].fix_dof("uy")
    for nid in (3, 4, 5, 6):
        n[nid].fix_dof("ux", 0.0)
    if impose:
        for nid in (5, 6):
            n[nid].fix_dof("uy", -2.0 * g0)
    d.generate_equation_numbers(verbose=False)
    return Assembler(d)


def _cantilever(E: float, n_el: int, L: float = 1.0, b: float = 0.01, h: float = 0.01,
                delta: float = 0.01, density: float | None = None):
    """Voladizo Euler-Bernoulli con flecha impuesta ``delta`` en la punta.
    Devuelve el ensamblador y la fuerza en la punta ``P = 3EI·δ/L³``."""
    from solidum import Frame2DEuler
    d = Domain()
    for i in range(n_el + 1):
        d.add_node(i + 1, [L * i / n_el, 0.0])
    n = d.nodes
    A, I = b * h, b * h ** 3 / 12.0
    for i in range(n_el):
        d.add_element(Frame2DEuler(i + 1, [n[i + 1], n[i + 2]], Elastic1D(E=E, density=density),
                                   A=A, I=I))
    for dof in Frame2DEuler.DOF_NAMES:
        n[1].fix_dof(dof)
    n[n_el + 1].fix_dof("uy", delta)
    d.generate_equation_numbers(verbose=False)
    return Assembler(d), 3.0 * E * I * delta / L ** 3


class TestForceFreeConvergence(unittest.TestCase):
    """Sin fuerzas, la convergencia no depende de la unidad de fuerza.

    Antes de la enmienda el piso del criterio de fuerza era 1e-9 *unidades
    de fuerza*: el bloque que baja por la holgura se detenía en N (su residuo
    de redondeo, ~2e-9 N, no bajaba de 1e-9 N) y convergía en kN y MN; el
    voladizo con flecha impuesta aceptaba en MN un residuo 13 veces mayor
    que en N; el voladizo de 3000 elementos se detenía en N.
    """

    def test_gap_closing_without_force_converges_alike_in_N_kN_and_MN(self):
        U = {}
        for unit, E in (("N", 3.0e10), ("kN", 3.0e7), ("MN", 3.0e4)):
            asm = _gap_block(E)
            U[unit] = NonlinearSolver(asm, num_steps=30).solve(np.zeros(asm.domain.total_dofs))
        # Longitudes en m en los tres casos: los desplazamientos coinciden.
        np.testing.assert_allclose(U["kN"], U["N"], rtol=1e-9, atol=1e-15)
        np.testing.assert_allclose(U["MN"], U["N"], rtol=1e-9, atol=1e-15)
        # El bloque tocó: la cara inferior bajó más que la holgura.
        n3_uy = asm.domain.nodes[3].dofs["uy"]
        self.assertLess(U["N"][n3_uy], -1.0e-3)

    def test_stage_starting_from_a_force_free_translation_is_in_equilibrium(self):
        """Una etapa que parte del bloque trasladado dentro de la holgura (sin
        fuerzas) empieza en equilibrio: su residuo es sólo redondeo, en N y en
        MN, y no se avisa de un desequilibrio."""
        import logging
        for unit, E in (("N", 3.0e10), ("MN", 3.0e4)):
            asm = _gap_block(E)
            d = asm.domain
            u0 = np.zeros(d.total_dofs)
            for nid in (3, 4, 5, 6):
                u0[d.nodes[nid].dofs["uy"]] = -0.5e-3
            with self.subTest(unit=unit), self.assertNoLogs("solidum", level=logging.WARNING):
                NonlinearSolver(asm, num_steps=10).solve(np.zeros(d.total_dofs), u0=u0)

    def test_free_flight_towards_contact_converges_alike_in_N_kN_and_MN(self):
        """Newton-Newmark: el bloque cae con velocidad inicial por la holgura
        sin ninguna fuerza (escalas de calibración nulas). Prueba de cobertura,
        no de regresión: con el criterio anterior también convergía, porque una
        traslación idéntica de todos los nodos da deformación exactamente nula
        en doble precisión (medido). La fórmula del piso de Newmark la protege
        ``test_each_newton_problem_floors_with_the_operators_of_its_residual``."""
        from solidum.math.solvers import NewtonNewmarkSolver
        v0, t_end, dt = 1.0, 0.9e-3, 2.5e-5   # toca a 1 ms: sólo vuelo libre
        u_end = {}
        for unit, E, rho in (("N", 3.0e10, 2400.0), ("kN", 3.0e7, 2.4), ("MN", 3.0e4, 2.4e-3)):
            asm = _gap_block(E, density=rho, impose=False)
            d = asm.domain
            v = np.zeros(d.total_dofs)
            for nid in (3, 4, 5, 6):
                v[d.nodes[nid].dofs["uy"]] = -v0
            res = NewtonNewmarkSolver(asm, t_end=t_end, dt=dt, u0_dot=v).solve()
            self.assertTrue(res.converged, unit)
            u_end[unit] = res.u_history[:, -1]
        np.testing.assert_allclose(u_end["kN"], u_end["N"], rtol=1e-9, atol=1e-15)
        np.testing.assert_allclose(u_end["MN"], u_end["N"], rtol=1e-9, atol=1e-15)
        # Vuelo libre: traslación rígida u = −v0·t.
        uy3 = u_end["N"][asm.domain.nodes[3].dofs["uy"]]
        self.assertAlmostEqual(uy3 / (-v0 * t_end), 1.0, places=9)

    def test_force_tolerance_relative_to_tip_force_is_the_same_in_N_and_MN(self):
        ratios = {}
        for unit, E in (("N", 2.0e11), ("MN", 2.0e5)):
            asm, P = _cantilever(E, n_el=10)
            crit = _RecordingCriterion()
            NonlinearSolver(asm, convergence=crit, num_steps=5).solve(np.zeros(asm.domain.total_dofs))
            ratios[unit] = crit.records[-1][0].tol_force / P
        self.assertAlmostEqual(ratios["MN"] / ratios["N"], 1.0, places=9)
        # Y es la de rtol (1e-5 de la fuerza corriente): el piso no afloja.
        self.assertLess(ratios["N"], 2.0e-5)

    def test_roundoff_floor_is_negligible_when_forces_are_resolvable(self):
        asm, P = _cantilever(2.0e11, n_el=10)
        crit = _RecordingCriterion()
        NonlinearSolver(asm, convergence=crit, num_steps=5).solve(np.zeros(asm.domain.total_dofs))
        state, roundoff = crit.records[-1]
        # Medido: 3e-6 de la tolerancia (con ε·max|diag K|·‖U‖ era 4e-4: la
        # rigidez axial EA/h, la mayor de la diagonal, no es la que actúa).
        self.assertLess(crit.roundoff_factor * roundoff / state.tol_force, 1.0e-4)

    def test_roundoff_scale_tracks_the_residual_noise_with_rotations_in_m_mm_and_km(self):
        """Marco con giros: tras converger, el residuo que queda (redondeo)
        está entre 0.05 y 1 veces la escala del piso en m, mm y km (medido:
        0.18–0.32). Con ``max|diag K|·‖U‖`` el cociente bajaba a 1e-3 en m y a
        1e-6 en km: el piso sobrestimaba el ruido hasta un millón de veces,
        porque la diagonal mayor mezcla rigideces de giro y de traslación."""
        for unit, E, L, b, h, delta in (("m", 2.0e11, 1.0, 0.01, 0.01, 0.01),
                                        ("mm", 2.0e5, 1000.0, 10.0, 10.0, 10.0),
                                        ("km", 2.0e2, 1.0e-3, 1.0e-5, 1.0e-5, 1.0e-5)):
            asm, _ = _cantilever(E, n_el=10, L=L, b=b, h=h, delta=delta)
            crit = _ExtraIterations()
            NonlinearSolver(asm, convergence=crit, num_steps=2, max_iter=10).solve(
                np.zeros(asm.domain.total_dofs))
            # Las evaluaciones forzadas tras la convergencia: residuo = redondeo.
            noise = [st.err_force / r for st, r, forced in crit.records if forced and r > 0]
            with self.subTest(unit=unit):
                self.assertGreater(len(noise), 4)
                self.assertGreater(min(noise), 0.05)
                self.assertLess(max(noise), 1.0)

    def test_dynamic_continuation_of_a_force_free_static_does_not_warn(self):
        """Un estático con flecha impuesta (sin fuerzas) seguido de un
        dinámico: el estado inicial está en equilibrio y no se avisa (antes,
        residuo relativo 1.0: el propio ruido frente a sí mismo)."""
        import logging
        import solidum
        from solidum.math.solvers import NewmarkSolver, NewtonNewmarkSolver
        for static in ("lineal", "no lineal"):
            for Dyn in (NewmarkSolver, NewtonNewmarkSolver):
                asm, _ = _cantilever(2.0e11, n_el=10, density=7850.0)
                d = asm.domain
                kw = {} if static == "lineal" else {"solver": NonlinearSolver(asm, num_steps=2)}
                r0 = solidum.run(d, assembler=asm, **kw)
                with self.subTest(static=static, dynamic=Dyn.__name__), \
                        self.assertNoLogs("solidum", level=logging.WARNING):
                    solidum.run_transient(d, assembler=asm,
                                          solver=Dyn(asm, t_end=0.002, dt=1e-3), previous=r0)

    def test_residual_at_roundoff_level_converges_in_a_very_refined_cantilever(self):
        asm, P = _cantilever(2.0e11, n_el=3000)
        crit = _RecordingCriterion()
        U = NonlinearSolver(asm, convergence=crit, num_steps=5).solve(np.zeros(asm.domain.total_dofs))
        # La flecha impuesta se alcanzó y el giro de la punta es el de la viga
        # de Euler-Bernoulli con carga en la punta, θ = 3δ/(2L).
        tip = asm.domain.nodes[3001]
        self.assertAlmostEqual(U[tip.dofs["uy"]], 0.01, places=12)
        self.assertAlmostEqual(U[tip.dofs["rz"]] / (1.5 * 0.01), 1.0, places=6)


# ---------------------------------------------------------------------------
# Ancla del piso de redondeo (ADR 0007, segunda enmienda del 2026-11-27)
# ---------------------------------------------------------------------------

def _plastic_cantilever(n_el: int = 4):
    """Voladizo vertical de 3 m con sección de 20 fibras de acero
    perfectamente plástico (σ_y = 250 MPa) y carga de referencia lateral de
    1 MN en la punta. La capacidad lateral medida es ≈ 1.68 MN."""
    from solidum import Frame2DEuler
    from solidum.materials.plastic_1d import Elastoplastic1D
    from solidum.sections import FiberSection, RectangularPatch
    L, B, H = 3.0, 0.30, 0.50
    d = Domain()
    nodes = [d.add_node(i + 1, [0.0, i * L / n_el]) for i in range(n_el + 1)]
    for k in range(n_el):
        steel = Elastoplastic1D(E=2.0e11, sigma_y=250.0e6, H=0.0)
        section = FiberSection(patches=[RectangularPatch(
            steel, y=(-H / 2, H / 2), z=(-B / 2, B / 2), ny=20, nz=1)])
        d.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], section=section, n_ip=3))
    for dof in ("ux", "uy", "rz"):
        nodes[0].fix_dof(dof, 0.0)
    d.generate_equation_numbers(verbose=False)
    F = np.zeros(d.total_dofs)
    F[nodes[-1].dofs["ux"]] = 1.0e6
    return d, F, nodes[-1]


class _GuidedToMechanism:
    """Problema de un paso cuyo Newton recorre una secuencia fija de
    iterados, sin depender del redondeo de ninguna factorización: (1) la
    corrección de Newton real; (2) un salto al mecanismo de colapso, con la
    rótula en la base y el resto en movimiento de sólido rígido (giro
    ``2**50``, representable exactamente); (3) una corrección de 1 m en la
    punta, despreciable frente a ``rtol_disp·‖U‖``."""

    def __init__(self, problem, U_mechanism, tip_dof):
        self._problem = problem
        self._U_mechanism = U_mechanism
        self._tip_dof = tip_dof
        self._k = 0

    def __getattr__(self, name):
        return getattr(self._problem, name)

    def correction(self, U, state, R, solve):
        self._k += 1
        if self._k == 1:
            return self._problem.correction(U, state, R, solve)
        if self._k == 2:
            return self._U_mechanism - U
        dU = np.zeros_like(U)
        dU[self._tip_dof] = 1.0
        return dU


class _ScriptedProblem:
    """Problema de dos grados de libertad cuyo Newton aplica una lista fija
    de correcciones; residuo constante ``residual`` y piso nulo. Con
    ``revise_once``, la primera convergencia se revisa y el intento se repite
    desde ``x0`` (ADR 0041)."""

    def __init__(self, corrections, residual=1.0, revise_once=False):
        self._corrections = [np.asarray(c, dtype=float) for c in corrections]
        self._k = 0
        self._residual = float(residual)
        self._pending_revision = revise_once

    def assemble(self, x):
        return None

    def residual(self, x, state):
        return np.array([self._residual, 0.0])

    def residual_norm(self, R):
        return float(np.linalg.norm(R))

    def calibration_scales(self, x, state):
        return 1.0, 1.0

    def reference_force(self, x, state):
        return 1.0

    def roundoff_force(self, x, state):
        return 0.0

    def x_norm(self, x):
        return float(np.linalg.norm(x))

    def correction(self, x, state, R, solve):
        dx = self._corrections[self._k]
        self._k += 1
        return dx

    def apply(self, x, dx, alpha):
        return x + alpha * dx, float(np.linalg.norm(alpha * dx))

    def on_converged(self, x, state):
        pass

    def revise(self, x):
        changed, self._pending_revision = self._pending_revision, False
        return changed


class _AnchorRecorder(ConvergenceCriterion):
    """Criterio que guarda ``(u_norm, delta_u_norm, roundoff_anchor)`` de
    cada evaluación."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.calls: list[tuple[float, float, float]] = []

    def evaluate(self, residual_norm, ref_force, delta_u_norm, u_norm, **kwargs):
        self.calls.append((float(u_norm), float(delta_u_norm), kwargs["roundoff_anchor"]))
        return super().evaluate(residual_norm, ref_force, delta_u_norm, u_norm, **kwargs)


class TestRoundoffFloorAnchor(unittest.TestCase):
    """El piso de redondeo no crece con un iterado que se aleja más de ``G``
    veces del ancla del intento: un Newton que diverge no puede ganar
    tolerancia con el tamaño de su propio iterado."""

    def test_a_newton_driven_into_a_plastic_mechanism_is_not_accepted_by_the_roundoff_floor(self):
        """Carga de 2.0 MN, por encima de la capacidad: no hay equilibrio. En
        el mecanismo ‖U‖ ~ 5e15 m, el residuo (acotado por σ_y) es ~3e7 N y
        el piso sin ancla ~6e9 N: el criterio sin ancla lo aceptaba (medido;
        en el CI, la misma aceptación en el colapso de una corrida real)."""
        d, F, tip = _plastic_cantilever()
        solver = NonlinearSolver(Assembler(d), num_steps=10, max_iter=3)
        problem = solver.make_problem(5.0 * F, 0.40)
        theta = 2.0 ** 50
        U_mechanism = np.zeros(d.total_dofs)
        for node in d.nodes.values():
            y = node.coordinates[1]
            if y > 0.0:
                U_mechanism[node.dofs["ux"]] = theta * y
                U_mechanism[node.dofs["rz"]] = -theta
        guided = _GuidedToMechanism(problem, U_mechanism, tip.dofs["ux"])
        res = solver.corrector.run(guided, np.zeros(d.total_dofs), check_initial=False)
        self.assertFalse(res.converged)
        # La prueba actúa: sin ancla, el piso de ese iterado cubriría el
        # residuo y el criterio de desplazamiento ya se cumple.
        conv = res.last_conv
        unanchored_floor = solver.convergence.roundoff_factor * problem.roundoff_force(
            res.x, res.state)
        self.assertGreater(unanchored_floor, 100.0 * res.last_residual)
        self.assertLess(conv.ratio_disp, 1.0)
        self.assertGreater(conv.ratio_force, 1.0e4)

    def test_the_floor_of_an_iterate_beyond_G_times_the_anchor_is_that_of_norm_G_a(self):
        c = ConvergenceCriterion(rtol_force=0.0, atol_force_factor=0.0, roundoff_factor=10.0,
                                 roundoff_anchor_factor=3.0)
        c.calibrate(0.0, 0.0)
        rf, a = 2.0e-6, 0.5
        for u, s in ((0.0, 1.0), (1.0, 1.0), (1.5, 1.0), (15.0, 0.1), (1.5e12, 1.0e-12)):
            with self.subTest(u_norm=u):
                state = c.evaluate(0.0, 0.0, 0.0, u, roundoff_force=rf, roundoff_anchor=a)
                self.assertAlmostEqual(state.floor_force / (10.0 * rf * s), 1.0, places=12)
                self.assertEqual(state.tol_force, state.floor_force)
        # G = 0 quita el piso en cuanto el iterado no es nulo.
        c0 = ConvergenceCriterion(roundoff_anchor_factor=0.0)
        c0.calibrate(0.0, 0.0)
        self.assertEqual(c0.evaluate(0.0, 0.0, 0.0, 1.0, roundoff_force=rf,
                                     roundoff_anchor=a).floor_force, 0.0)
        self.assertEqual(c0.evaluate(0.0, 0.0, 0.0, 0.0, roundoff_force=rf,
                                     roundoff_anchor=a).floor_force, 10.0 * rf)

    def test_by_default_the_anchor_factor_is_the_newton_kantorovich_bound_3(self):
        self.assertEqual(ConvergenceCriterion().roundoff_anchor_factor, 3.0)

    def test_the_anchor_factor_is_validated_and_read_from_yaml(self):
        for value in (-1.0, float("nan"), float("inf"), "tres", True):
            with self.subTest(value=value), \
                    self.assertRaisesRegex(ValueError, "roundoff_anchor_factor"):
                ConvergenceCriterion(roundoff_anchor_factor=value)
        c = make_convergence_from_config({"roundoff_anchor_factor": 10.0})
        self.assertEqual(c.roundoff_anchor_factor, 10.0)

    def test_the_corrector_anchors_each_attempt_at_its_first_corrected_iterate(self):
        """``a = ‖x_1‖ + ‖δx_0‖``, fijo en el resto del intento aunque el
        iterado se dispare; el iterado 0 (predictor) se evalúa sin tope."""
        from solidum.math.solvers.corrector import NewtonCorrector
        crit = _AnchorRecorder()
        corrector = NewtonCorrector(crit, max_iter=3, is_symmetric=True)
        corrector.run(_ScriptedProblem([[3.0, 4.0], [0.0, 100.0], [0.0, 1.0e6]]),
                      np.zeros(2), check_initial=False)
        self.assertEqual([call[2] for call in crit.calls], [10.0, 10.0, 10.0])
        # Otro intento, otra ancla: con predictor, el iterado 0 va sin ancla
        # y la primera corrección fija la del intento.
        crit.calls.clear()
        corrector.run(_ScriptedProblem([[0.0, 2.0], [0.0, 1.0e6], [0.0, 1.0e6]]),
                      np.array([1.0, 0.0]), check_initial=True, initial_delta_norm=2.0)
        a = float(np.hypot(1.0, 2.0)) + 2.0
        self.assertEqual([call[2] for call in crit.calls], [None, a, a, a])

    def test_an_evaluation_without_anchor_keeps_the_full_floor_with_any_G(self):
        """El predictor y el estado inicial de una etapa no tienen iterados
        de Newton que anclar: ``roundoff_anchor=None`` no recorta el piso,
        ni con ``G < 1``."""
        for G in (0.0, 0.5, 3.0):
            c = ConvergenceCriterion(rtol_force=0.0, atol_force_factor=0.0, roundoff_factor=10.0,
                                     roundoff_anchor_factor=G)
            c.calibrate(0.0, 0.0)
            with self.subTest(G=G):
                state = c.evaluate(0.0, 0.0, 0.0, 1.0e12, roundoff_force=2.0e-6,
                                   roundoff_anchor=None)
                self.assertEqual(state.floor_force, 10.0 * 2.0e-6)

    def test_a_revision_of_the_converged_step_anchors_the_repeated_attempt_afresh(self):
        from solidum.math.solvers.corrector import NewtonCorrector
        crit = _AnchorRecorder(rtol_disp=10.0)
        corrector = NewtonCorrector(crit, max_iter=3, is_symmetric=True)
        res = corrector.run(_ScriptedProblem([[3.0, 4.0], [0.0, 7.0]], residual=0.0,
                                             revise_once=True),
                            np.zeros(2), check_initial=False)
        self.assertTrue(res.converged)
        self.assertEqual(res.n_revisions, 1)
        self.assertEqual([call[2] for call in crit.calls], [10.0, 14.0])


if __name__ == '__main__':
    unittest.main()
