# fenix_fem/fenix/math/solvers/nonlinear.py
"""``NonlinearSolver`` — Newton-Raphson incremental con paso adaptativo.
"""
from typing import Optional

import numpy as np

from fenix.math.convergence import (
    ConvergenceCriterion,
    stiffness_diag_scale,
)
from fenix.math.linalg import StiffnessProperties, select_solver
from fenix.math.solvers._shared import (
    CholeskyNotPositiveDefiniteError,
    _log,
    domain_is_symmetric,
)
from fenix.registry import SolverRegistry


@SolverRegistry.register
class NonlinearSolver:
    """Solucionador incremental-iterativo de Newton-Raphson con paso adaptativo."""
    def __init__(self, assembler, convergence: Optional[ConvergenceCriterion] = None, max_iter=20, num_steps=10, adaptive=True, min_delta_lambda=1e-5, linear_algebra: str = "auto", freeze_tangent_after_iter: int | None = None):
        self.assembler = assembler
        # Política de convergencia (ADR 0007). El criterio se calibra una vez
        # al inicio de solve() con la escala del primer ensamblaje.
        self.convergence = convergence if convergence is not None else ConvergenceCriterion()
        self.max_iter = max_iter
        self.num_steps = num_steps
        self.adaptive = adaptive
        self.min_delta_lambda = min_delta_lambda
        self.linear_algebra = linear_algebra
        # Newton modificado (ADR 0003 fase 2): si es int, se factoriza fresca
        # las primeras N iteraciones del paso y se reusa la factorización en
        # las siguientes. None ⇒ Newton estándar (factoriza cada iteración).
        self.freeze_tangent_after_iter = freeze_tangent_after_iter
        # Se inicializa al comienzo de solve(); se degrada a LU si Cholesky aborta.
        self._linalg = None
        self._is_pd = True
        self._frozen_factor = None  # FactorizedSolver | None

    def _make_linalg(self, ndof: int):
        props = StiffnessProperties(
            is_symmetric=domain_is_symmetric(self.assembler.domain),
            is_positive_definite=self._is_pd,
            size=ndof,
        )
        return select_solver(props, override=self.linear_algebra)

    def _solve_reduced(self, K_red, R_red, iteration: int = 0):
        """Resuelve K_red·δU_red = R_red con fallback SPD→LU y Newton modificado.

        Si ``freeze_tangent_after_iter`` está activo, factoriza fresca solo
        durante las primeras N iteraciones del paso y reusa la factorización
        cacheada en las siguientes (ADR 0003 §5 + fase 2).
        """
        if self.freeze_tangent_after_iter is None:
            try:
                return self._linalg.solve(K_red, R_red)
            except CholeskyNotPositiveDefiniteError:
                return self._fallback_to_lu_and_solve(K_red, R_red)

        threshold = self.freeze_tangent_after_iter
        if iteration < threshold or self._frozen_factor is None:
            try:
                self._frozen_factor = self._linalg.factorize(K_red)
            except CholeskyNotPositiveDefiniteError:
                self._is_pd = False
                self._linalg = self._make_linalg(K_red.shape[0])
                _log.warning("Cholesky reportó no-positividad. Degradando a LU para el resto del análisis.")
                self._frozen_factor = self._linalg.factorize(K_red)
        return self._frozen_factor.solve(R_red)

    def _fallback_to_lu_and_solve(self, K_red, R_red):
        _log.warning("Cholesky reportó no-positividad. Degradando a LU para el resto del análisis.")
        self._is_pd = False
        self._linalg = self._make_linalg(K_red.shape[0])
        return self._linalg.solve(K_red, R_red)

    def solve(self, F_ext_global: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U_current = np.zeros(ndof)

        _log.info("--- INICIANDO SOLVER NO LINEAL (CONTROL DE PASO ADAPTATIVO) ---")

        # Tras reducción el solver algebraico opera sobre n_libre, no sobre ndof.
        cs = self.assembler.constraint_set
        n_free = ndof - len(cs)
        free_dofs = cs.free_dofs(ndof)
        self._linalg = self._make_linalg(n_free)

        load_factor = 0.0
        target_load = 1.0
        delta_lambda = 1.0 / self.num_steps
        step = 0

        while load_factor < target_load - 1e-9:
            step += 1

            if load_factor + delta_lambda > target_load:
                delta_lambda = target_load - load_factor

            next_load_factor = load_factor + delta_lambda
            _log.info(f"[PASO {step}] Intentando Factor de Carga: {next_load_factor:.4f} (Incremento: {delta_lambda:.4f})")

            F_ext_step = F_ext_global * next_load_factor
            U_iter = U_current.copy()
            converged = False

            for iteration in range(self.max_iter):
                K_global, F_int_global = self.assembler.assemble_non_linear_system(U_iter)
                R = F_ext_step - F_int_global

                K_red, R_red, T_op, g_inc = self.assembler.reduce(
                    K_global, R, U_current=U_iter, load_factor=next_load_factor
                )

                # Calibración del criterio en el primer ensamblaje de la corrida
                # (ADR 0007). Las escalas se derivan del estado inicial real:
                # ‖F_ext_global‖ para fuerza, ‖F_ext‖/max|diag(K)| para
                # desplazamiento. Si ‖F_ext_global‖ = 0 (control en desplazamiento
                # puro), usar la reacción interna del primer paso o fallback 1.0.
                if not self.convergence.is_calibrated:
                    force_scale = max(
                        np.linalg.norm(F_ext_global),
                        np.linalg.norm(F_int_global),
                        1.0,
                    )
                    K_diag = stiffness_diag_scale(K_global)
                    disp_scale = force_scale / K_diag
                    self.convergence.calibrate(force_scale, disp_scale)

                try:
                    delta_U_red = self._solve_reduced(K_red, R_red, iteration=iteration)
                except RuntimeError:
                    _log.error("Matriz Singular detectada.")
                    break

                delta_U = self.assembler.expand(delta_U_red, T_op, g_inc)
                U_iter += delta_U

                # Criterio dual fuerza + desplazamiento (ADR 0007). El residuo se
                # evalúa en DOFs libres (en prescritos R contiene la reacción
                # física, no un fallo de equilibrio).
                ref_force = max(np.linalg.norm(F_ext_step), np.linalg.norm(F_int_global))
                state = self.convergence.evaluate(
                    residual_norm=np.linalg.norm(R[free_dofs]),
                    ref_force=ref_force,
                    delta_u_norm=np.linalg.norm(delta_U),
                    u_norm=np.linalg.norm(U_iter),
                )
                _log.info(
                    f"  Iteración {iteration+1:2d} | "
                    f"R/tol_F: {state.ratio_force:.4e} | "
                    f"dU/tol_d: {state.ratio_disp:.4e}"
                )

                if state.converged:
                    _log.info("  -> CONVERGENCIA ALCANZADA.")
                    self.assembler.commit_all_states()

                    U_current = U_iter
                    load_factor = next_load_factor
                    converged = True

                    if self.adaptive and iteration < 4 and delta_lambda < (1.0 / self.num_steps):
                        delta_lambda = min(delta_lambda * 1.5, 1.0 / self.num_steps)
                        _log.info(f"  -> Acelerando el próximo incremento a {delta_lambda:.4f}")

                    if step_callback:
                        step_callback(step, U_current, load_factor)

                    break

            # Cierre de paso (converja o no): la K_t cambia con U_current,
            # así que la factorización congelada deja de ser válida para el
            # siguiente paso (Newton modificado, ADR 0003 fase 2).
            self._frozen_factor = None

            if not converged:
                if self.adaptive:
                    delta_lambda /= 2.0
                    _log.warning(f"NO CONVERGIÓ. Bisección: reduciendo incremento a {delta_lambda:.4f}")
                    if delta_lambda < self.min_delta_lambda:
                        raise RuntimeError(f"El incremento de carga ({delta_lambda:.2e}) cayó por debajo del mínimo ({self.min_delta_lambda:.2e}). El solver ha divergido.")
                else:
                    raise RuntimeError(f"El solucionador no convergió en el paso {step}.")

        return U_current
