# fenix_fem/fenix/math/solvers/nonlinear.py
"""``NonlinearSolver`` — Newton-Raphson incremental con paso adaptativo.
"""
from __future__ import annotations

import numpy as np

from fenix.constants import (
    LINE_SEARCH_C1,
    LINE_SEARCH_MAX_BACKTRACKS,
    LINE_SEARCH_RHO,
)
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
from fenix.math.solvers.diagnostics import (
    SolverDivergedError,
    UnknownDivergenceError,
    classify_divergence,
)
from fenix.registry import SolverRegistry


@SolverRegistry.register
class NonlinearSolver:
    """Solucionador incremental-iterativo de Newton-Raphson con paso adaptativo."""

    PIPELINE_KIND = "static"

    def __init__(self, assembler, convergence: ConvergenceCriterion | None = None, max_iter=20, num_steps=10, adaptive=True, min_delta_lambda=1e-5, linear_algebra: str = "auto", freeze_tangent_after_iter: int | None = None, line_search: bool = False):
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
        # Line search por descenso no monótono (ADR 0011). Default ``False``:
        # el line search corrige el modo de oscilación documentado en
        # Drucker-Prager perfectamente plástico, pero es contraproducente
        # en regímenes donde Newton avanza correctamente a través de
        # transitorios de residuo creciente (daño con tangente consistente,
        # plasticidad cerca de la rama postcrítica). Activar explícitamente
        # cuando se observe oscilación.
        self.line_search = bool(line_search)
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

    def _armijo_step(self, U_iter: np.ndarray, delta_U: np.ndarray,
                     R_norm_current: float, F_ext_step: np.ndarray,
                     free_dofs):
        """Line search por descenso no monótono (ADR 0011, variante GLL).

        Devuelve ``(α, R_after_norm, F_int_after)`` donde ``α ∈ (0, 1]`` es
        el primer factor que satisface ``‖R(U + α·δU)‖ ≤ ‖R(U)‖`` (la
        condición Grippo-Lampariello-Lucidi 1986 simplificada — descenso
        no monótono). Más permisiva que Armijo puro: acepta ``α = 1``
        cuando Newton baja el residuo, **sin** exigir suficiente
        decrecimiento al modo Wolfe. Solo hace backtracking cuando el paso
        completo de Newton produce ``R`` mayor.

        Justificación: Armijo puro con ``c₁ > 0`` puede rechazar pasos de
        Newton correctos en problemas FEM no lineales (daño activo,
        plasticidad cerca de la rama postcrítica) donde el residuo a veces
        sube transitoriamente antes de converger. La condición de
        descenso no monótono preserva la velocidad cuadrática del Newton
        estándar y solo interviene cuando hay un rebote claro.

        Si el backtracking se agota sin encontrar un α que baje el
        residuo, devuelve los valores del último α probado y delega al
        control externo (oscillation, bisección del paso).

        Tras esta llamada, el ``state.vars_trial`` de los elementos
        corresponde al ensamblaje en ``U_iter + α·δU``: el bucle exterior
        evalúa la convergencia coherentemente y ``commit_all_states()``
        promueve el trial coherente con el U efectivamente avanzado.

        Cuando ``self.line_search=False`` ensambla una vez con α=1 para
        que el state trial quede coherente con el U avanzado (semántica
        equivalente al código pre-ADR 0011).
        """
        rho = LINE_SEARCH_RHO
        max_bt = LINE_SEARCH_MAX_BACKTRACKS
        # `LINE_SEARCH_C1` no se usa en esta variante (GLL no exige
        # suficiente decrecimiento). Se conserva en constants.py por
        # legibilidad y por si se introduce una variante Wolfe en futuro.

        if not self.line_search:
            U_trial = U_iter + delta_U
            _, F_int_trial = self.assembler.assemble_non_linear_system(U_trial)
            R_trial = F_ext_step - F_int_trial
            R_trial_norm = float(np.linalg.norm(R_trial[free_dofs]))
            return 1.0, R_trial_norm, F_int_trial

        alpha = 1.0
        F_int_trial = None
        R_trial_norm = float("inf")

        for _ in range(max_bt + 1):
            U_trial = U_iter + alpha * delta_U
            _, F_int_trial = self.assembler.assemble_non_linear_system(U_trial)
            R_trial = F_ext_step - F_int_trial
            R_trial_norm = float(np.linalg.norm(R_trial[free_dofs]))

            if R_trial_norm <= R_norm_current:
                return alpha, R_trial_norm, F_int_trial
            alpha *= rho

        return alpha, R_trial_norm, F_int_trial

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
        n_bisections = 0  # contador de bisecciones globales del paso (ADR 0011)

        while load_factor < target_load - 1e-9:
            step += 1

            if load_factor + delta_lambda > target_load:
                delta_lambda = target_load - load_factor

            next_load_factor = load_factor + delta_lambda
            _log.info(f"[PASO {step}] Intentando Factor de Carga: {next_load_factor:.4f} (Incremento: {delta_lambda:.4f})")

            F_ext_step = F_ext_global * next_load_factor
            U_iter = U_current.copy()
            converged = False

            # ADR 0010 §5: hook de preparación de paso. Evaluado con el estado
            # convergido del paso anterior para evitar chattering por
            # predictores lineales dentro del Newton. No-op para elementos
            # estándar; lo usan elementos con discontinuidad embebida
            # (CST_Embedded2D) para chequear activación.
            self.assembler.prepare_all_steps(U_current)

            # Historial de residuos del paso para clasificación de divergencia
            # (ADR 0011). Se resetea al inicio de cada intento de paso.
            residual_history: list[float] = []
            delta_history: list[float] = []
            singular_tangent_seen = False
            last_alpha = 1.0
            last_residual = float("inf")
            last_delta = 0.0

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
                    singular_tangent_seen = True
                    break

                delta_U = self.assembler.expand(delta_U_red, T_op, g_inc)

                # Line search Armijo con backtracking (ADR 0011). Devuelve también
                # el ‖R‖ y F_int post-actualización, garantizando que el state
                # trial queda coherente con el U avanzado para la evaluación de
                # convergencia y el commit posterior.
                R_norm_before = float(np.linalg.norm(R[free_dofs]))
                alpha, R_norm_after, F_int_after = self._armijo_step(
                    U_iter, delta_U, R_norm_before, F_ext_step, free_dofs,
                )
                last_alpha = alpha
                U_iter = U_iter + alpha * delta_U

                # Criterio dual fuerza + desplazamiento (ADR 0007), evaluado en
                # el U_iter actualizado (R_norm_after corresponde al state trial
                # coherente con este U). Residuo en DOFs libres.
                ref_force = max(np.linalg.norm(F_ext_step), np.linalg.norm(F_int_after))
                delta_U_norm = float(np.linalg.norm(alpha * delta_U))
                state = self.convergence.evaluate(
                    residual_norm=R_norm_after,
                    ref_force=ref_force,
                    delta_u_norm=delta_U_norm,
                    u_norm=np.linalg.norm(U_iter),
                )
                residual_history.append(R_norm_after)
                delta_history.append(delta_U_norm)
                last_residual = R_norm_after
                last_delta = delta_U_norm

                alpha_tag = f" | α={alpha:.3f}" if alpha < 1.0 else ""
                _log.info(
                    f"  Iteración {iteration+1:2d} | "
                    f"R/tol_F: {state.ratio_force:.4e} | "
                    f"dU/tol_d: {state.ratio_disp:.4e}{alpha_tag}"
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
                    n_bisections += 1
                    _log.warning(f"NO CONVERGIÓ. Bisección: reduciendo incremento a {delta_lambda:.4f}")
                    if delta_lambda < self.min_delta_lambda:
                        # Clasificar el modo y lanzar excepción tipada (ADR 0011).
                        err_cls = classify_divergence(
                            residual_history, delta_history,
                            singular_tangent_detected=singular_tangent_seen,
                        )
                        raise err_cls(
                            last_residual=last_residual,
                            last_delta=last_delta,
                            last_load_factor=next_load_factor,
                            n_bisections=n_bisections,
                            extra_message=(
                                f"Δλ={delta_lambda:.2e} < min={self.min_delta_lambda:.2e}; "
                                f"line_search={'on' if self.line_search else 'off'}; "
                                f"último α={last_alpha:.3f}."
                            ),
                        )
                else:
                    err_cls = classify_divergence(
                        residual_history, delta_history,
                        singular_tangent_detected=singular_tangent_seen,
                    )
                    raise err_cls(
                        last_residual=last_residual,
                        last_delta=last_delta,
                        last_load_factor=next_load_factor,
                        n_bisections=n_bisections,
                        extra_message=(
                            f"Paso {step} no convergió en {self.max_iter} iter; "
                            f"adaptive=False; line_search={'on' if self.line_search else 'off'}."
                        ),
                    )

        return U_current
