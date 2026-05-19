"""Internos del `NonlinearSolver` — convergencia cuadrática, line search y
paso adaptativo.

Tres clases de tests que cubren propiedades del `NonlinearSolver` no
ejercitadas por los tests de robustez de `test_solver_robustness.py`
(mapeo de modos de fallo) ni por los tests indirectos vía plasticidad/daño
(que validan el resultado físico, no el comportamiento del solver):

1. **Convergencia cuadrática del Newton-Raphson** (Simo & Hughes §3.4,
   Belytschko §6.4.10). Verifica que con tangente consistente la sucesión
   de residuos satisface ``r_{n+1} ≤ C·r_n²`` en el régimen asintótico.
   Diagnóstico estándar para detectar bugs sutiles en `compute_state`
   (tangente que difiere ligeramente de ``∂σ/∂ε``).

2. **Line search del ADR 0011**. El default es `line_search=False`; el
   helper de descenso no monótono (Grippo-Lampariello-Lucidi) implementado
   en la fase C del ADR no estaba ejercitado por ningún test. Este bloque
   cubre el smoke test del code path activado + verificación directa de
   que `_armijo_step` hace backtracking cuando el paso completo aumenta
   el residuo.

3. **Paso adaptativo**: bisección automática cuando un intento no
   converge en `max_iter` iteraciones; comportamiento sin adaptivo
   (excepción tipada). Verifica el contrato declarado en el docstring
   del solver.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.domain import Domain
from fenix.elements.truss import Truss2D, Truss2DCorot
from fenix.materials.damage_1d import IsotropicDamage1D
from fenix.materials.elastic import Elastic1D
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.solvers import NonlinearSolver
from fenix.math.solvers.diagnostics import SolverDivergedError


class _ResidualRecorder:
    """Proxy de `ConvergenceCriterion` que captura `(residual_norm, delta_u_norm)`
    en cada llamada a `evaluate()`. Delega el resto al criterio interno.

    El solver llama a `convergence.evaluate(...)` una vez por iteración
    Newton (todas las iteraciones de todos los pasos), por lo que las
    listas resultantes son la secuencia cronológica completa del análisis.
    """

    def __init__(self, inner: ConvergenceCriterion):
        self._inner = inner
        self.residuals: list[float] = []
        self.deltas: list[float] = []

    @property
    def is_calibrated(self):
        return self._inner.is_calibrated

    @property
    def atol_force(self):
        return self._inner.atol_force

    @property
    def atol_disp(self):
        return self._inner.atol_disp

    def calibrate(self, *args, **kwargs):
        return self._inner.calibrate(*args, **kwargs)

    def evaluate(self, residual_norm, ref_force, delta_u_norm, u_norm):
        state = self._inner.evaluate(residual_norm, ref_force, delta_u_norm, u_norm)
        self.residuals.append(float(residual_norm))
        self.deltas.append(float(delta_u_norm))
        return state


def _build_plastic_truss(sigma_y: float = 100.0, H: float = 1.0e3,
                         E: float = 1.0e5, L: float = 1.0, A: float = 1.0):
    """Truss elastoplástico 1D estándar con hardening lineal pequeño.

    `H` pequeño (≈ E/100) garantiza tangente consistente bien condicionada
    y plasticidad clara al cargar por encima del yield. La curva σ-ε es
    bilineal — Newton converge exacto en 2 iter (1 predictor elástico + 1
    corrección plástica).
    """
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    mat = Elastoplastic1D(E=E, sigma_y=sigma_y, H=H)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n2.fix_dof('uy', 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), n2.dofs['ux']


def _build_shallow_von_mises_truss(L: float = 10.0, h: float = 1.0,
                                    E: float = 1.0e5, A: float = 1.0):
    """Shallow von Mises 2-bar (h/L=0.1), corotacional puro elástico.

    Geometría con non-linealidad geométrica suave en la rama pre-pico.
    Tangente consistente = K_M + K_G; Newton con carga sub-pico itera
    ~4-6 veces con convergencia cuadrática asintótica clara.
    """
    dom = Domain()
    n_left = dom.add_node(1, [0.0, 0.0])
    n_apex = dom.add_node(2, [L, h])
    n_right = dom.add_node(3, [2.0 * L, 0.0])
    mat = Elastic1D(E=E)
    dom.add_element(Truss2DCorot(1, [n_left, n_apex], mat, A=A))
    dom.add_element(Truss2DCorot(2, [n_apex, n_right], mat, A=A))
    n_left.fix_dof('ux', 0.0); n_left.fix_dof('uy', 0.0)
    n_right.fix_dof('ux', 0.0); n_right.fix_dof('uy', 0.0)
    n_apex.fix_dof('ux', 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), n_apex.dofs['uy']


def _build_damage_bar_1d(E: float = 1.0e5, A: float = 1.0, L: float = 1.0,
                          kappa_0: float = 1.0e-3, alpha: float = 200.0):
    """Barra empotrada-libre con IsotropicDamage1D — softening exponencial.

    Pico de carga elástico ``F_pico = E·A·κ_0`` (default = 100 N). El
    softening abrupto (`alpha=200`) hace que cargas moderadamente por
    encima del pico requieran muchas iteraciones (Newton oscila si la
    tangente es negativa) o impongan bisecciones del paso.
    """
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    mat = IsotropicDamage1D(E=E, kappa_0=kappa_0, alpha=alpha)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n2.fix_dof('uy', 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), n2.dofs['ux']


# =============================================================================
# 1. Convergencia cuadrática del Newton-Raphson
# =============================================================================

class TestNewtonQuadraticConvergence(unittest.TestCase):
    """Convergencia cuadrática asintótica del Newton con tangente consistente.

    Diagnóstico canónico (Simo & Hughes 1998, Belytschko et al. 2014) de la
    calidad de la tangente devuelta por `compute_state`. Una tangente que
    no es exactamente `∂σ/∂ε` produce convergencia lineal (cada iteración
    gana un dígito) en vez de cuadrática (cada iteración duplica los
    dígitos correctos) — bug silencioso que no rompe la suite pero
    degrada el rendimiento y puede hacer el solver no robusto.
    """

    def test_shallow_von_mises_shows_quadratic_jump_in_asymptotic_regime(self):
        """En las últimas iteraciones, alguna transición debe mostrar caída
        de al menos 2 órdenes de magnitud en el residuo (``r_{n+1} ≤ r_n²``
        cuando ``r_n ≤ 1``).

        El shallow von Mises 2-bar con corotacional puro tiene tangente
        consistente `K_M + K_G` (Cap. 9 Crisfield Vol.1). Cargado bajo el
        punto límite, Newton itera 4-6 veces con régimen asintótico claro.
        Convergencia lineal (tangente errónea) cae típicamente ~0.5 por
        iter; convergencia cuadrática cae cuadráticamente cuando el
        residuo entra en régimen pequeño. El test pide al menos UNA
        transición cuadrática limpia.

        Por qué no uso Elastoplastic1D: su curva σ-ε es bilineal exacta,
        Newton converge en 2 iter sin régimen asintótico medible.
        """
        dom, asm, dof_uy = _build_shallow_von_mises_truss()
        # Carga sub-pico (P_pico ≈ 38 N para h/L=0.1 con E·A=1e5): el
        # control de carga llega sin necesidad de arc-length.
        F = np.zeros(dom.total_dofs)
        F[dof_uy] = -20.0  # ~0.5·P_pico

        recorder = _ResidualRecorder(
            ConvergenceCriterion(rtol_force=1e-12, rtol_disp=1e-12)
        )
        solver = NonlinearSolver(
            asm, convergence=recorder, num_steps=1, max_iter=20,
            adaptive=False,
        )
        solver.solve(F)

        residuals = [r for r in recorder.residuals if r > 0.0]
        self.assertGreaterEqual(len(residuals), 3,
            f"Newton convergió en {len(residuals)} iter; insuficiente "
            "para evaluar régimen asintótico (esperadas ≥3).")

        # Debe ser monótonamente decreciente desde alguna iteración en
        # adelante (después de la iteración elástica inicial podría no
        # serlo perfectamente, pero el régimen asintótico sí lo es).
        last_three = residuals[-3:]
        self.assertTrue(
            last_three[1] < last_three[0] and last_three[2] < last_three[1],
            f"Los 3 últimos residuos no son decrecientes: {last_three}.",
        )

        # Test cuadrático: al menos UNA transición consecutiva en las
        # últimas iteraciones debe mostrar caída de ≥ 2 órdenes de
        # magnitud (log10 ratio ≤ -2). Convergencia lineal pura nunca
        # alcanza esa caída.
        found_quadratic = False
        for r_n, r_n1 in zip(residuals[:-1], residuals[1:]):
            if r_n > 0 and r_n1 > 0:
                drop = np.log10(r_n1) - np.log10(r_n)
                if drop <= -2.0:
                    found_quadratic = True
                    break
        self.assertTrue(found_quadratic,
            f"Ninguna transición muestra caída cuadrática (Δlog10 ≤ -2): "
            f"residuos={['%.3e' % r for r in residuals]}. Posible bug en "
            "la tangente consistente de Elastoplastic1D.")


# =============================================================================
# 2. Line search del ADR 0011
# =============================================================================

class TestLineSearch(unittest.TestCase):
    """Code path del line search opt-in (`line_search=True`)."""

    def test_line_search_preserves_answer_on_well_behaved_problem(self):
        """En un problema donde Newton avanza sin oscilación, activar
        `line_search` no debe cambiar la respuesta convergida.

        Sin este test, el code path `line_search=True` puede romperse
        silenciosamente sin que ningún test lo detecte (el default es
        `False` y ningún otro test lo activa).
        """
        dom_a, asm_a, dof_a = _build_plastic_truss(sigma_y=100.0, H=1.0e3)
        dom_b, asm_b, dof_b = _build_plastic_truss(sigma_y=100.0, H=1.0e3)
        F_a = np.zeros(dom_a.total_dofs); F_a[dof_a] = 150.0
        F_b = np.zeros(dom_b.total_dofs); F_b[dof_b] = 150.0

        solver_off = NonlinearSolver(asm_a, num_steps=4, line_search=False)
        solver_on = NonlinearSolver(asm_b, num_steps=4, line_search=True)

        U_off = solver_off.solve(F_a)
        U_on = solver_on.solve(F_b)

        np.testing.assert_allclose(U_off, U_on, rtol=1e-8, atol=1e-12,
            err_msg="line_search=True dio respuesta distinta a line_search=False "
                    "en un problema bien comportado.")

    def test_armijo_step_backtracks_when_full_step_increases_residual(self):
        """`_armijo_step` con `line_search=True` reduce α cuando ``α=1``
        aumenta el residuo.

        Construye manualmente un escenario: truss elástico 1D con
        `delta_U` deliberadamente exagerado (10× el incremento del
        Newton correcto). El paso completo lleva el residuo a ~9× su
        valor inicial; el backtracking debe encontrar ``α < 1`` que
        baje el residuo.
        """
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [1.0, 0.0])
        mat = Elastic1D(E=1.0e5)
        dom.add_element(Truss2D(1, [n1, n2], mat, A=1.0))
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n2.fix_dof('uy', 0.0)
        dom.generate_equation_numbers(verbose=False)
        asm = Assembler(dom)

        ndof = dom.total_dofs
        free_dofs = asm.constraint_set.free_dofs(ndof)

        F_ext_step = np.zeros(ndof)
        F_ext_step[n2.dofs['ux']] = 100.0  # solución exacta: u_x = 1e-3
        U_iter = np.zeros(ndof)
        delta_U = np.zeros(ndof)
        delta_U[n2.dofs['ux']] = 1.0e-2  # 10× el incremento correcto

        R_norm_current = float(np.linalg.norm(F_ext_step[free_dofs]))

        solver = NonlinearSolver(asm, line_search=True)
        alpha, R_after, _ = solver._armijo_step(
            U_iter, delta_U, R_norm_current, F_ext_step, free_dofs,
        )

        self.assertLess(alpha, 1.0,
            f"Line search no reaccionó al paso exagerado: α={alpha}. "
            "Esperado backtracking activo.")
        self.assertLessEqual(R_after, R_norm_current,
            f"Backtracking no logró bajar el residuo: α={alpha}, "
            f"R_after={R_after:.3e}, R_current={R_norm_current:.3e}.")

    def test_armijo_step_disabled_takes_full_step(self):
        """Con `line_search=False`, `_armijo_step` siempre devuelve α=1
        sin importar el efecto sobre el residuo (semántica pre-ADR 0011).
        """
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [1.0, 0.0])
        mat = Elastic1D(E=1.0e5)
        dom.add_element(Truss2D(1, [n1, n2], mat, A=1.0))
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n2.fix_dof('uy', 0.0)
        dom.generate_equation_numbers(verbose=False)
        asm = Assembler(dom)
        ndof = dom.total_dofs
        free_dofs = asm.constraint_set.free_dofs(ndof)
        F_ext_step = np.zeros(ndof); F_ext_step[n2.dofs['ux']] = 100.0
        delta_U = np.zeros(ndof); delta_U[n2.dofs['ux']] = 1.0e-2  # 10× exagerado

        solver = NonlinearSolver(asm, line_search=False)
        alpha, _, _ = solver._armijo_step(
            np.zeros(ndof), delta_U, 100.0, F_ext_step, free_dofs,
        )
        self.assertEqual(alpha, 1.0,
            f"Con line_search=False, α debe ser exactamente 1.0; obtuvo {alpha}.")


# =============================================================================
# 3. Paso adaptativo
# =============================================================================

class TestAdaptiveStepControl(unittest.TestCase):
    """Lógica del paso adaptativo: bisección automática al fallar, fallo
    limpio cuando el adaptivo está desactivado."""

    def test_step_bisects_when_iter_budget_exhausted(self):
        """Con `max_iter` deliberadamente bajo, el primer intento no
        converge → el adaptivo biseca el incremento y reintenta hasta
        terminar el análisis.

        Patrón: shallow von Mises 2-bar cargado a 0.7·P_pico con
        `num_steps=1` (incremento completo) y `max_iter=2` (presupuesto
        de iteración insuficiente). Newton no cabe en el primer intento;
        el solver biseca, reintenta con `δλ=0.5`, sigue. Al final el
        análisis llega a ``λ=1`` con más pasos convergidos que el
        `num_steps` inicial.
        """
        dom, asm, dof_uy = _build_shallow_von_mises_truss()
        F = np.zeros(dom.total_dofs); F[dof_uy] = -25.0  # ~0.7·P_pico

        trace: list[tuple[int, float]] = []

        def cb(step, _U, lam):
            trace.append((step, float(lam)))

        solver = NonlinearSolver(
            asm, num_steps=1, max_iter=2, adaptive=True,
            min_delta_lambda=1.0e-6,
        )
        solver.solve(F, step_callback=cb)

        self.assertGreater(len(trace), 1,
            f"Sin bisección activa: {trace}. Con num_steps=1 y max_iter=2 "
            "se esperaba ≥2 pasos convergidos tras bisecar.")
        self.assertAlmostEqual(trace[-1][1], 1.0, places=6,
            msg=f"Análisis no completó (λ_final={trace[-1][1]:.4f}); trace={trace}.")
        # Step counter monótono creciente (los intentos fallidos no
        # aparecen en el callback, pero el contador del solver sigue
        # avanzando — si el primer step del callback es > 1, es
        # evidencia directa de bisección al primer intento).
        steps_only = [s for s, _ in trace]
        self.assertEqual(steps_only, sorted(steps_only),
            f"Step counter no es monótono creciente: {steps_only}.")

    def test_non_adaptive_raises_typed_exception_on_first_failure(self):
        """Con `adaptive=False`, un primer intento fallido lanza una
        subclase de `SolverDivergedError` (ADR 0011) sin bisecar.

        Mismo problema del test anterior con `adaptive=False`. El Newton
        no converge en `max_iter=2` iteraciones; sin bisección, el solver
        clasifica el modo de divergencia y lanza la excepción tipada.
        """
        dom, asm, dof_uy = _build_shallow_von_mises_truss()
        F = np.zeros(dom.total_dofs); F[dof_uy] = -25.0

        solver = NonlinearSolver(
            asm, num_steps=1, max_iter=2, adaptive=False,
        )
        with self.assertRaises(SolverDivergedError) as ctx:
            solver.solve(F)

        # La excepción tipada debe llevar el contador de bisecciones a
        # cero (porque adaptive=False) y un residuo real > 0.
        self.assertEqual(ctx.exception.n_bisections, 0,
            f"Con adaptive=False, n_bisections debe ser 0; "
            f"obtuvo {ctx.exception.n_bisections}.")
        self.assertGreater(ctx.exception.last_residual, 0.0)


if __name__ == '__main__':
    unittest.main()
