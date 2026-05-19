"""``DissipationArcLengthSolver`` — arc-length con restricción basada en disipación.

Variante de :class:`ArcLengthSolver` que cambia la restricción del paso de la
cuadrática cilíndrica (``‖ΔU‖² = dl²``) a la lineal por disipación de energía
(``g(ΔU, Δλ) = τ``) de Gutiérrez 2004 / Verhoosel et al. 2009. Reusa el
plumbing del padre (predictor tangente, corrector Newton, manejo de
Dirichlet/MPC, backend algebraico) y solo redefine la lógica de paso.

Está diseñado para atravesar:

- Régimen elástico inicial (modo `"cylindrical"` heredado del padre).
- Transición elástico→softening con penalty stiff ``K_e`` del embedded
  discontinuity (ADR 0010 fase 4), donde el cilíndrico falla por
  cuasi-degeneración del discriminante de la cuadrática cerca del pico.
- Rama post-pico estable de daño/cohesivos (modo `"dissipation"`).
- Descarga global con grieta saturada (volver a `"cylindrical"`).

El switching cilíndrico↔disipación es **automático** según se detecte
disipación neta por encima de un umbral relativo. La spec completa en
``docs/specs/DissipationArcLengthSolver.md``.
"""
from __future__ import annotations

import numpy as np

from fenix.constants import ARCLENGTH_MIN_DL_FACTOR, ZERO_TOL
from fenix.math.convergence import (
    ConvergenceCriterion,
    stiffness_diag_scale,
)
from fenix.math.solvers._shared import _log
from fenix.math.solvers.arclength import ArcLengthSolver
from fenix.registry import SolverRegistry


@SolverRegistry.register
class DissipationArcLengthSolver(ArcLengthSolver):
    """Arc-length con restricción de disipación (Gutiérrez 2004).

    Subclase de :class:`ArcLengthSolver`. Reusa los parámetros y la
    infraestructura del padre, y añade los específicos de la restricción
    de disipación: ``initial_tau`` y sus factores adaptativos.

    Parameters
    ----------
    initial_tau : float
        Incremento de disipación objetivo (unidades de energía: F·U).
        Calibrar al orden de ``G_F · l_d / n_pasos_post_pico``.
    tau_grow_factor, tau_max_factor, tau_shrink_factor : float
        Análogos a ``dl_*_factor`` del padre pero sobre τ.
    tau_grow_iter_threshold, tau_shrink_iter_threshold : int
        Análogos a ``dl_*_iter_threshold`` del padre pero sobre τ.
    dissipation_threshold : float
        Umbral relativo (respecto a ``‖F_ref‖·‖U‖``) bajo el cual la
        disipación neta del paso se considera 0 — el solver permanece
        o vuelve a modo cilíndrico.
    **arclength_kwargs
        Resto de parámetros, idénticos a :class:`ArcLengthSolver`.

    Notes
    -----
    La fórmula central (Gutiérrez 2004 ec. 13, Verhoosel et al. 2009 ec. 11):

    .. math::
        g(\\Delta\\mathbf U, \\Delta\\lambda) = \\tfrac{1}{2}\\left(
            \\lambda_n\\,\\mathbf F^T_{\\text{ref}}\\,\\Delta\\mathbf U
            - \\Delta\\lambda\\,\\mathbf F^T_{\\text{ref}}\\,\\mathbf U_n
        \\right) = \\tau

    es **lineal** en ``(ΔU, Δλ)`` — el corrector resuelve una ecuación
    lineal en ``ddλ`` (no cuadrática como el padre), eliminando la
    selección de raíz por menor ángulo y la patología de raíces
    imaginarias en pasos grandes.

    En el primer paso (``λ_n = 0``, ``U_n = 0``) el denominador α de la
    restricción es 0; por eso el solver siempre **arranca en modo
    cilíndrico** y solo conmuta a disipación tras detectar
    ``ΔE_d > umbral`` en un paso committed.
    """

    PIPELINE_KIND = "static"

    def __init__(self, assembler, convergence: ConvergenceCriterion | None = None,
                 max_iter: int = 20, max_lambda: float = 1.0,
                 initial_dl: float = 0.1, max_steps: int = 100,
                 dl_grow_factor: float = 1.5, dl_max_factor: float = 5.0,
                 dl_shrink_factor: float = 0.6,
                 dl_grow_iter_threshold: int = 4,
                 dl_shrink_iter_threshold: int = 8,
                 linear_algebra: str = "auto",
                 *,
                 initial_tau: float,
                 tau_grow_factor: float = 1.5, tau_max_factor: float = 5.0,
                 tau_shrink_factor: float = 0.6,
                 tau_grow_iter_threshold: int = 4,
                 tau_shrink_iter_threshold: int = 8,
                 dissipation_threshold: float = 1.0e-12):
        super().__init__(
            assembler, convergence=convergence, max_iter=max_iter,
            max_lambda=max_lambda, initial_dl=initial_dl, max_steps=max_steps,
            dl_grow_factor=dl_grow_factor, dl_max_factor=dl_max_factor,
            dl_shrink_factor=dl_shrink_factor,
            dl_grow_iter_threshold=dl_grow_iter_threshold,
            dl_shrink_iter_threshold=dl_shrink_iter_threshold,
            linear_algebra=linear_algebra,
        )
        if initial_tau <= 0.0:
            raise ValueError(
                f"DissipationArcLengthSolver: initial_tau={initial_tau} debe ser > 0."
            )
        self.initial_tau = float(initial_tau)
        self.tau_grow_factor = tau_grow_factor
        self.tau_max_factor = tau_max_factor
        self.tau_shrink_factor = tau_shrink_factor
        self.tau_grow_iter_threshold = tau_grow_iter_threshold
        self.tau_shrink_iter_threshold = tau_shrink_iter_threshold
        self.dissipation_threshold = dissipation_threshold
        # Estado: modo inicial cilíndrico; siempre arranca así porque la
        # restricción de disipación es degenerada en (λ=0, U=0).
        self._mode: str = "cylindrical"
        self._tau: float = self.initial_tau

    @staticmethod
    def _dissipation_increment(lambda_n: float, U_n: np.ndarray,
                                F_ref: np.ndarray, dU_step: np.ndarray,
                                dlambda_step: float) -> float:
        """Disipación incremental ``ΔE_d`` del paso (Gutiérrez 2004 ec. 13).

        Forma escalar de la restricción ``g(ΔU, Δλ)``:

        ``ΔE_d = ½ · (λ_n·F_ref·ΔU − Δλ·F_ref·U_n)``

        Positiva en pasos disipativos (carga + softening); cercana a 0 en
        régimen elástico puro.
        """
        return 0.5 * (
            lambda_n * float(F_ref @ dU_step)
            - dlambda_step * float(F_ref @ U_n)
        )

    def _negative_pivots(self, K) -> int | None:
        """Aproximación inicial al sign-of-pivot tracking vía signo del
        determinante.

        Por la ley de inercia de Sylvester, el número de pivots negativos
        en una factorización LDLᵀ es invariante bajo congruencia, pero LU
        no preserva inercia. El signo del **determinante** sí distingue
        entre "número par" y "número impar" de pivots negativos, lo cual
        basta para detectar el paso por un punto límite simple
        (transición ``0 → 1 → 0 → ...``).

        Para tracking exacto se requiere LDLᵀ Bunch-Kaufman (deuda técnica
        documentada en la spec).

        Returns
        -------
        int | None
            ``0`` si el determinante es positivo (asume PD inicial → 0
            pivots negativos), ``1`` si el determinante es negativo
            (número impar de pivots negativos), ``None`` si la matriz es
            singular o el cómputo falla.
        """
        try:
            K_dense = K.toarray() if hasattr(K, "toarray") else np.asarray(K)
            sign, _ = np.linalg.slogdet(K_dense)
        except Exception:
            return None
        if sign > 0:
            return 0
        if sign < 0:
            return 1
        return None

    def solve(self, F_ext_ref: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U_current = np.zeros(ndof)
        lambda_curr = 0.0
        step = 0
        dl = self.dl
        # Reset del estado adaptativo al inicio de una corrida nueva.
        self._tau = self.initial_tau
        self._mode = "cylindrical"

        delta_U_step = np.zeros(ndof)

        _log.info(
            "--- INICIANDO DISSIPATION-ARCLENGTH SOLVER (Gutiérrez 2004) ---"
        )

        cs = self.assembler.constraint_set
        n_free = ndof - len(cs)
        self._linalg = self._make_linalg(n_free)

        while lambda_curr < self.max_lambda and step < self.max_steps:
            step += 1
            _log.info(
                f"[PASO {step}] modo={self._mode!s} dl={dl:.4e} tau={self._tau:.4e}"
            )

            U_iter = U_current.copy()
            lambda_iter = lambda_curr
            converged = False

            # ADR 0010 §5: hook de preparación de paso.
            self.assembler.prepare_all_steps(U_current)

            # --- 1. PREDICTOR ---
            K_global, F_int_global = self.assembler.assemble_non_linear_system(U_iter)

            if not self.convergence.is_calibrated:
                force_scale = max(
                    np.linalg.norm(F_ext_ref),
                    np.linalg.norm(F_int_global),
                    1.0,
                )
                K_diag = stiffness_diag_scale(K_global)
                disp_scale = force_scale / K_diag
                self.convergence.calibrate(force_scale, disp_scale)

            K_t_red, F_t_red, T_t, g_t = self.assembler.reduce(
                K_global, F_ext_ref.copy(),
            )
            try:
                du_t_red = self._solve(K_t_red, F_t_red)
            except RuntimeError:
                _log.error("Matriz singular en predictor. Bisección...")
                dl /= 2.0
                self._tau /= 2.0
                continue
            du_t = self.assembler.expand(du_t_red, T_t, g_t)

            # Sentido del avance.
            sign = 1.0
            if step > 1 and np.dot(delta_U_step, du_t) < 0:
                sign = -1.0

            # dlambda según modo activo.
            mode_predictor = self._mode
            if mode_predictor == "dissipation":
                # α = ½·(λ_n·F·du_t − F·U_n)
                alpha = 0.5 * (
                    lambda_curr * float(F_ext_ref @ du_t)
                    - float(F_ext_ref @ U_current)
                )
                # Para problemas lineales monotónicos, α ≡ 0 exactamente
                # (U_n = λ_n·K^-1·F_ref ⇒ F·U_n = λ_n·F·K^-1·F = λ_n·F·du_t).
                # Por eso necesitamos un threshold **relativo** a la escala
                # energética del paso: si |α| ≪ ½·λ_n·F·du_t (= energía
                # elástica del estado escalada), el modo disipación no es
                # físicamente significativo y se revierte a cilíndrico.
                alpha_scale = abs(0.5 * lambda_curr * float(F_ext_ref @ du_t))
                if alpha_scale > 0 and abs(alpha) < 1.0e-6 * alpha_scale:
                    _log.info(
                        f"  α/escala = {abs(alpha)/alpha_scale:.2e} ≪ 1 "
                        f"⇒ régimen lineal: cilíndrico temporal."
                    )
                    mode_predictor = "cylindrical"
                    # Resetear el modo del solver también para evitar
                    # re-entrar en disipación en el siguiente paso sin
                    # disipación física genuina.
                    self._mode = "cylindrical"
                elif abs(alpha) < ZERO_TOL:
                    _log.info("  α ≈ 0 estrictamente → cilíndrico temporal.")
                    mode_predictor = "cylindrical"
                    self._mode = "cylindrical"

            if mode_predictor == "cylindrical":
                dlambda = sign * dl / (np.linalg.norm(du_t) + ZERO_TOL)
            else:  # dissipation
                dlambda = sign * abs(self._tau / alpha)

            # Cierre exacto en max_lambda.
            final_step = (
                sign > 0 and lambda_curr + dlambda >= self.max_lambda - ZERO_TOL
            )
            # Salvaguarda contra ``final_step`` prematuro en problemas con
            # softening severo: si ``dλ_pred`` excede ``max_lambda`` por un
            # factor importante, el Newton interno hereda un punto inicial
            # lejano que rebota en la singularidad del pico. Bisectar la
            # longitud característica (dl o τ) en lugar de aceptar el paso
            # final salva ese caso.
            #
            # Casos típicos donde se dispara:
            #   - Pasos 1-3 con dl/τ aún sin calibrar a la rigidez real.
            #   - Inmediatamente tras un switch cilíndrico→disipación, donde
            #     α puede ser pequeño y τ/α salta a ≫ Δλ del último paso
            #     cilíndrico.
            overshoot_factor = 3.0
            remaining = self.max_lambda - lambda_curr
            if final_step and abs(dlambda) > overshoot_factor * remaining:
                _log.info(
                    f"  Predictor excede max_lambda en {abs(dlambda)/remaining:.1f}×. "
                    f"Bisecando longitud de paso (modo={mode_predictor})."
                )
                if mode_predictor == "cylindrical":
                    dl *= 0.5
                else:
                    self._tau *= 0.5
                continue
            if final_step:
                dlambda = self.max_lambda - lambda_curr

            lambda_iter += dlambda
            dU_iter = dlambda * du_t
            U_iter += dU_iter

            # --- 2. CORRECTOR ITERATIVO ---
            for iteration in range(self.max_iter):
                K_global, F_int_global = self.assembler.assemble_non_linear_system(U_iter)
                R = lambda_iter * F_ext_ref - F_int_global

                K_t_red, F_t_red, T_t, g_t = self.assembler.reduce(
                    K_global, F_ext_ref.copy(),
                )
                K_red, R_red, T_R, g_R = self.assembler.reduce(
                    K_global, R, U_current=U_iter, load_factor=lambda_iter,
                )
                try:
                    du_R_red = self._solve(K_red, R_red)
                    du_t_red = self._solve(K_t_red, F_t_red)
                except RuntimeError:
                    _log.error("Matriz singular en corrector.")
                    break
                du_R = self.assembler.expand(du_R_red, T_R, g_R)
                du_t = self.assembler.expand(du_t_red, T_t, g_t)

                if final_step:
                    # Cierre exacto: lambda fijo, corrección puramente Newton.
                    ddlambda = 0.0
                    dU_update = du_R
                elif mode_predictor == "cylindrical":
                    # Restricción cuadrática (idéntica al padre).
                    dU_new = dU_iter + du_R
                    a = np.dot(du_t, du_t)
                    b = 2.0 * np.dot(dU_new, du_t)
                    c = np.dot(dU_new, dU_new) - dl**2
                    det = b**2 - 4.0 * a * c
                    if det < 0:
                        _log.error("Raíces imaginarias en cilíndrico.")
                        break
                    ddl1 = (-b + np.sqrt(det)) / (2.0 * a)
                    ddl2 = (-b - np.sqrt(det)) / (2.0 * a)
                    theta1 = np.dot(dU_iter, dU_new + ddl1 * du_t)
                    theta2 = np.dot(dU_iter, dU_new + ddl2 * du_t)
                    ddlambda = ddl1 if theta1 > theta2 else ddl2
                    dU_update = du_R + ddlambda * du_t
                    dU_iter = dU_new + ddlambda * du_t
                else:  # mode_predictor == "dissipation"
                    # Restricción lineal Gutiérrez en ddλ:
                    #   g(ΔU + du_R + ddλ·du_t, Δλ_pre + ddλ) = τ
                    # ⇒ ddλ = (τ − g_partial) / α
                    d_lambda_pre = lambda_iter - lambda_curr
                    alpha = 0.5 * (
                        lambda_curr * float(F_ext_ref @ du_t)
                        - float(F_ext_ref @ U_current)
                    )
                    g_partial = 0.5 * (
                        lambda_curr * float(F_ext_ref @ (dU_iter + du_R))
                        - d_lambda_pre * float(F_ext_ref @ U_current)
                    )
                    if abs(alpha) < ZERO_TOL:
                        _log.error("α ≈ 0 en corrector de disipación. Aborto del paso.")
                        break
                    ddlambda = (self._tau - g_partial) / alpha
                    dU_update = du_R + ddlambda * du_t
                    dU_iter = dU_iter + dU_update

                lambda_iter += ddlambda
                if final_step:
                    dU_iter = dU_iter + dU_update
                U_iter = U_current + dU_iter

                ref_force = max(
                    np.linalg.norm(F_ext_ref) * abs(lambda_iter),
                    np.linalg.norm(F_int_global),
                )
                state = self.convergence.evaluate(
                    residual_norm=np.linalg.norm(R[cs.free_dofs(ndof)]),
                    ref_force=ref_force,
                    delta_u_norm=np.linalg.norm(dU_update),
                    u_norm=np.linalg.norm(U_iter),
                )
                _log.info(
                    f"  Iter. {iteration+1:2d} | lam={lambda_iter:.4f} | "
                    f"R/tol_F: {state.ratio_force:.4e} | "
                    f"dU/tol_d: {state.ratio_disp:.4e}"
                )

                if state.converged:
                    # Disipación del paso completo (post-corrección).
                    dU_step_total = U_iter - U_current
                    dlambda_step_total = lambda_iter - lambda_curr
                    dE_d = self._dissipation_increment(
                        lambda_curr, U_current, F_ext_ref,
                        dU_step_total, dlambda_step_total,
                    )
                    _log.info(
                        f"  -> CONV. lam={lambda_iter:.4f} ΔE_d={dE_d:.4e}"
                    )

                    self.assembler.commit_all_states()

                    # Switching cilíndrico ↔ disipación tras commit.
                    ref_energy = max(
                        np.linalg.norm(F_ext_ref) * np.linalg.norm(U_iter),
                        1.0,
                    )
                    dE_threshold = self.dissipation_threshold * ref_energy
                    if self._mode == "cylindrical":
                        if dE_d > dE_threshold:
                            _log.info(
                                f"  Switch cilíndrico→disipación "
                                f"(ΔE_d={dE_d:.3e} > {dE_threshold:.3e})"
                            )
                            self._mode = "dissipation"
                            # Inicializa τ con la disipación recién observada
                            # para continuidad del primer paso disipativo.
                            self._tau = max(dE_d, self.initial_tau)
                    else:  # _mode == "dissipation"
                        if dE_d < dE_threshold:
                            _log.info(
                                f"  Switch disipación→cilíndrico "
                                f"(ΔE_d={dE_d:.3e} < {dE_threshold:.3e})"
                            )
                            self._mode = "cylindrical"

                    # Adaptatividad según el modo del **siguiente** paso.
                    if self._mode == "cylindrical":
                        if iteration < self.dl_grow_iter_threshold:
                            dl = min(
                                dl * self.dl_grow_factor,
                                self.dl * self.dl_max_factor,
                            )
                        elif iteration > self.dl_shrink_iter_threshold:
                            dl *= self.dl_shrink_factor
                    else:
                        if iteration < self.tau_grow_iter_threshold:
                            self._tau = min(
                                self._tau * self.tau_grow_factor,
                                self.initial_tau * self.tau_max_factor,
                            )
                        elif iteration > self.tau_shrink_iter_threshold:
                            self._tau *= self.tau_shrink_factor

                    U_current = U_iter
                    lambda_curr = lambda_iter
                    delta_U_step = dU_iter
                    converged = True

                    if step_callback is not None:
                        step_callback(step, U_current, lambda_curr)
                    break

            if not converged:
                # Bisección por modo del paso fallido.
                if mode_predictor == "cylindrical":
                    dl *= 0.5
                    _log.warning(
                        f"Bisección cilíndrica: dl → {dl:.4e}"
                    )
                    if dl < ARCLENGTH_MIN_DL_FACTOR * self.dl:
                        raise RuntimeError(
                            "DissipationArcLengthSolver: dl bajó del umbral "
                            f"({ARCLENGTH_MIN_DL_FACTOR * self.dl:.4e}). "
                            "Aborto irreparable."
                        )
                else:
                    self._tau *= 0.5
                    _log.warning(
                        f"Bisección disipación: τ → {self._tau:.4e}"
                    )
                    if self._tau < ARCLENGTH_MIN_DL_FACTOR * self.initial_tau:
                        raise RuntimeError(
                            "DissipationArcLengthSolver: τ bajó del umbral "
                            f"({ARCLENGTH_MIN_DL_FACTOR * self.initial_tau:.4e}). "
                            "Aborto irreparable."
                        )

        return U_current
