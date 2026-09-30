"""``DissipationArcLengthSolver`` — arc-length con restricción basada en disipación.

Variante de :class:`ArcLengthSolver` que cambia la restricción del paso de la
cuadrática cilíndrica (``‖ΔU‖² = dl²``) a la lineal por disipación de energía
(``g(ΔU, Δλ) = τ``) de Gutiérrez 2004 / Verhoosel et al. 2009. Reusa el
plumbing del padre (predictor tangente, corrector Newton, manejo de
Dirichlet/MPC, backend algebraico) y solo redefine la lógica de paso.

Está diseñado para atravesar:

- Régimen elástico inicial (modo `"cylindrical"` heredado del padre).
- Transición elástico→ablandamiento con una rigidez tangente casi singular
  en el pico, donde el cilíndrico falla por cuasi-degeneración del
  discriminante de la cuadrática.
- Rama post-pico estable de ablandamiento (modo `"dissipation"`), incluidos
  los retrocesos (*snap-back*): la restricción fija el signo de Δλ.
- Descarga global una vez agotado el ablandamiento (volver a
  `"cylindrical"`).

El switching cilíndrico↔disipación es **automático**. Desde el 2026-09-30 el
predictor del modo de disipación y el arranque de ese modo usan la **tangente
de carga** ``K_L`` —la del último iterado convergido del paso anterior— en
lugar de la tangente del estado consolidado, que en daño con descarga secante
es la secante y anula el denominador de la restricción. La spec completa en
``docs/user/dissipation_control/specs/DissipationArcLengthSolver.md`` (§5).

Módulo de usuario ``dissipation_control`` (ADR 0020): se conecta al programa
principal por su parte pública —:class:`ArcLengthSolver` y el problema de
Newton por paso :class:`ArcProblem`, cuyo método ``constraint`` redefine— y
por los métodos que la clase de arco ofrece a sus variantes (predictor
tangente, primer paso, llegada a ``max_lambda`` y cierre del trazado). La
tangente de carga la toma del resultado público del corrector compartido
(``CorrectorResult.state``, el ``(K, F_int)`` del ensamblaje en el iterado
convergido) y la resuelve con el backend del corrector (``corrector.solve``)
y la reducción de Dirichlet/MPC del ensamblador (``reduce``/``expand``), las
mismas piezas que usa el predictor del padre.
"""
from __future__ import annotations

import numpy as np

from solidum.constants import ARCLENGTH_MIN_DL_FACTOR, ZERO_TOL
from solidum.logging import get_logger
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import ArcLengthSolver, ArcProblem
from solidum.math.solvers.corrector import CorrectionAborted
from solidum.math.solvers.model_checks import (
    check_initial_tangent,
    ensure_statically_restrained,
)
from solidum.registry import SolverRegistry

# Mismo logger que los solvers del programa principal ("solvers"): los filtros
# y la captura de registros no distinguen la procedencia del solver.
_log = get_logger("solvers")

# Umbral relativo de ARRANQUE del modo de disipación (decisión del usuario del
# 2026-09-30, spec §5): se conmuta de cilíndrico a disipación cuando
#   |α_L| > ALPHA_SWITCH_RTOL · ½·|λ_n·F_ref·δu_L|,   δu_L = K_L⁻¹·F_ref,
# es decir, cuando la tangente de carga se separa de la secante. Suelo de
# ruido: en régimen sin disipación K_L = K_s y α_L = −½·F·K⁻¹·r_n, con r_n el
# residuo de equilibrio del paso anterior, así que |α_L|/escala ≤
# (‖r_n‖/‖λ_n·F‖)/cos θ, con θ el ángulo entre F y K⁻¹·F: del orden de la
# tolerancia relativa de fuerza del criterio de convergencia
# (CONVERGENCE_RTOL_FORCE = 1e-5 por omisión; medido en la barra con
# retroceso: ≈ 1e-15 en la rama elástica, ≈ 1e-6 con la secante tras el
# inicio del daño). 1e-3 deja dos órdenes de margen; con daño en
# crecimiento la razón es O(0.1-1) (0.45 en el paso que conmuta en la barra:
# α_L = −0.42 J frente a una escala de 0.93 J). Si el usuario afloja el
# criterio de convergencia hasta rtol ≳ 1e-4, el margen desaparece.
ALPHA_SWITCH_RTOL = 1.0e-3

# Umbral relativo de PERMANENCIA (anterior al 2026-09-30, sin cambios): en un
# paso en modo de disipación, |α| < ALPHA_LINEAR_RTOL·escala se toma por
# régimen sin disipación y el paso vuelve a cilíndrico. Junto con el umbral de
# arranque forma una histéresis (entrar exige 1e-3, salir 1e-6) que evita que
# el modo oscile entre pasos contiguos.
ALPHA_LINEAR_RTOL = 1.0e-6


def _alpha_is_null(alpha: float, lambda_n: float, F_du: float, F_U_n: float) -> bool:
    """``α = ½·(λ_n·F·δu − F·U_n)`` es nulo a precisión de máquina.

    Guarda contra la división por cero en ``Δλ = τ/α`` y ``ddλ = (τ − g)/α``,
    relativa a los dos términos de los que α es diferencia: ``|α| ≤ ZERO_TOL·
    ½·(|λ_n·F·δu| + |F·U_n|)``. Adimensional (revisión del 2026-09-30): hasta
    entonces se comparaba ``|α|`` con ``ZERO_TOL`` en unidades de energía, y en
    un modelo con energías de ~1e-12 (una probeta de micras en SI, o la barra
    de validación con longitudes y fuerzas en unidades 1e6 veces mayores) el
    α de un paso de ablandamiento genuino caía por debajo y el solver volvía a
    cilíndrico en la esquina del pico. Con ambos términos nulos (λ_n = 0 y
    U_n = 0) la guarda salta: α = 0 exactamente.
    """
    return abs(alpha) <= ZERO_TOL * 0.5 * (abs(lambda_n * F_du) + abs(F_U_n))


class _DissipationArcProblem(ArcProblem):
    """``ArcProblem`` con el modo ``"dissipation"``: restricción lineal de
    Gutiérrez en ``ddλ`` (ver :class:`DissipationArcLengthSolver`). Los
    modos ``"newton"`` y ``"cylindrical"`` se heredan del padre."""

    def __init__(self, *args, tau: float, **kwargs):
        super().__init__(*args, **kwargs)
        self.tau = float(tau)

    def constraint(self, x, du_R, du_t):
        if self.mode != "dissipation":
            return super().constraint(x, du_R, du_t)
        _, lambda_iter, dU_iter = x
        F_ext_ref = self.F_ext_ref
        # Restricción lineal Gutiérrez en ddλ:
        #   g(ΔU + du_R + ddλ·du_t, Δλ_pre + ddλ) = τ
        # ⇒ ddλ = (τ − g_partial) / α
        d_lambda_pre = lambda_iter - self.lambda_curr
        F_du_t = float(F_ext_ref @ du_t)
        F_U_n = float(F_ext_ref @ self.U_current)
        alpha = 0.5 * (self.lambda_curr * F_du_t - F_U_n)
        g_partial = 0.5 * (
            self.lambda_curr * float(F_ext_ref @ (dU_iter + du_R))
            - d_lambda_pre * F_U_n
        )
        if _alpha_is_null(alpha, self.lambda_curr, F_du_t, F_U_n):
            raise CorrectionAborted("α ≈ 0 en corrector de disipación. Aborto del paso.")
        ddlambda = (self.tau - g_partial) / alpha
        return du_R + ddlambda * du_t, ddlambda


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
        Calibrar como fracción de la disipación total esperada en la rama
        post-pico: del orden de ``E_d,total / n_pasos_post_pico``.
    tau_grow_factor, tau_max_factor, tau_shrink_factor : float
        Análogos a ``dl_*_factor`` del padre pero sobre τ.
    tau_grow_iter_threshold, tau_shrink_iter_threshold : int
        Análogos a ``dl_*_iter_threshold`` del padre pero sobre τ.
    dissipation_threshold : float
        Umbral relativo (respecto a ``‖F_ref‖·‖U‖``) bajo el cual la
        disipación neta de un paso en modo de disipación se considera 0 —
        el solver vuelve a modo cilíndrico. Desde el 2026-09-30 no decide el
        arranque del modo de disipación (ver Notes).
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
    cilíndrico**.

    **Tangente de carga** (2026-09-30). Con el predictor ``ΔU = Δλ·δu``,
    ``g = α·Δλ`` con ``α = ½·(λ_n·F·δu − F·U_n)``. Si ``δu`` se calcula con
    la tangente del estado consolidado ``U_n``, el daño con descarga secante
    está en descarga y devuelve la secante ``K_s``; como
    ``K_s·U_n = λ_n·F + r_n``, resulta ``α = −½·F·K_s⁻¹·r_n ≈ 0``. El solver
    usa en su lugar ``δu_L = K_L⁻¹·F`` con ``K_L`` la tangente del último
    iterado convergido del paso anterior (prueba de carga de Kuhn-Tucker
    respecto del estado consolidado anterior):

    .. math::
        \\alpha_L = \\tfrac12\\,\\lambda_n\\,\\mathbf F^T
            (\\mathbf K_L^{-1} - \\mathbf K_s^{-1})\\,\\mathbf F

    mide el crecimiento del daño. Con ella:

    - en modo de disipación, el predictor es ``Δλ = τ/α_L`` (negativo en un
      retroceso) con ``ΔU = Δλ·δu_L``;
    - en modo cilíndrico, el paso siguiente pasa a disipación cuando
      ``|α_L| > ALPHA_SWITCH_RTOL·½·|λ_n·F·δu_L|``: la tangente de carga se
      separa de la secante. Sustituye al criterio anterior
      (``ΔE_d > dissipation_threshold·‖F‖‖U‖`` tras un paso cilíndrico), que
      en un pico anguloso no se alcanza porque el arco cilíndrico se atasca
      con pasos que apenas disipan.

    Vale dentro del alcance de la spec: pequeñas deformaciones y descarga
    secante. Con no linealidad geométrica ``α_L ≠ 0`` sin disipación.
    """

    PIPELINE_KIND = "static"

    def __init__(self, assembler, convergence: ConvergenceCriterion | None = None,
                 max_iter: int = 20, max_lambda: float = 1.0,
                 initial_dl: float | None = None, max_steps: int = 100,
                 dl_grow_factor: float = 1.5, dl_max_factor: float = 5.0,
                 dl_shrink_factor: float = 0.6,
                 dl_grow_iter_threshold: int = 4,
                 dl_shrink_iter_threshold: int = 8,
                 linear_algebra: str = "auto",
                 *,
                 initial_tau: float,
                 initial_dlambda: float | None = None,
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
            initial_dlambda=initial_dlambda,
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
        # Tangente de carga K_L: la del último iterado convergido del paso
        # anterior (ver Notes). None hasta el primer paso convergido.
        self._K_load = None

    def _load_tangent_direction(self, F_ext_ref: np.ndarray) -> np.ndarray | None:
        """``δu_L = K_L⁻¹·F_ref`` con la tangente de carga ``K_L``.

        Misma reducción de Dirichlet/MPC y mismo backend que el predictor
        tangente del padre (``ArcLengthSolver._tangent_predictor``), con la
        matriz del último iterado convergido en lugar de la ensamblada en el
        estado consolidado. Devuelve ``None`` si aún no hay paso convergido o
        si ``K_L`` es singular: el paso sigue entonces con la tangente del
        estado consolidado, como antes del 2026-09-30.
        """
        if self._K_load is None:
            return None
        K_red, F_red, T_t, g_t = self.assembler.reduce(self._K_load, F_ext_ref.copy())
        try:
            du_red = self.corrector.solve(K_red, F_red)
        except RuntimeError:
            _log.warning("  Tangente de carga K_L singular: predictor con la tangente "
                         "del estado consolidado.")
            return None
        return self.assembler.expand(du_red, T_t, g_t)

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

    def solve(self, F_ext_ref: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U_current = np.zeros(ndof)
        lambda_curr = 0.0
        step = 0
        dl = None          # Δl₁: se fija con el predictor del primer paso (padre)
        dl_ref = None
        self.dl_reference = None
        # Reset del estado adaptativo al inicio de una corrida nueva.
        self._tau = self.initial_tau
        self._mode = "cylindrical"
        self._K_load = None

        delta_U_step = np.zeros(ndof)

        _log.info(
            "--- INICIANDO DISSIPATION-ARCLENGTH SOLVER (Gutiérrez 2004) ---"
        )
        # Red de seguridad (ADR 0019): mecanismos rígidos antes de resolver y
        # pivotes nulos de la tangente inicial, que factoriza el predictor
        # del primer paso. Ver ArcLengthSolver.solve.
        ensure_statically_restrained(self.assembler, type(self).__name__)
        check_initial_tangent(self.assembler, self.corrector, type(self).__name__)

        cs = self.assembler.constraint_set
        free_dofs = cs.free_dofs(ndof)
        steps_done = 0

        while lambda_curr < self.max_lambda and step < self.max_steps:
            step += 1

            # Gancho de inicio de paso (ADR 0020, P5). Evaluado con el estado
            # convergido del paso anterior, como en ArcLengthSolver.solve.
            self.assembler.prepare_all_steps(U_current)

            # --- 1. PREDICTOR ---
            pred = self._tangent_predictor(U_current, F_ext_ref, delta_U_step, step)
            if pred is None:
                _log.error("Matriz singular en predictor. Bisección...")
                if dl is not None:
                    dl /= 2.0
                self._tau /= 2.0
                continue
            du_t, sign = pred
            if dl is None:
                dl = dl_ref = self._first_dl(du_t)
            _log.info(
                f"[PASO {step}] modo={self._mode!s} dl={dl:.4e} tau={self._tau:.4e}"
            )

            # Tangente de carga (2026-09-30, ver Notes de la clase): δu_L =
            # K_L⁻¹·F con K_L la tangente del último iterado convergido. La
            # tangente ensamblada por el predictor del padre en U_n es, con
            # daño y descarga secante, la secante K_s, con la que
            # α = −½·F·K_s⁻¹·r_n ≈ 0; con K_L, α_L mide el crecimiento del daño.
            du_L = self._load_tangent_direction(F_ext_ref)

            # Arranque del modo de disipación por α_L (criterio (b), spec §5):
            # la tangente de carga se separa de la secante. Sustituye al umbral
            # de ΔE_d tras un paso cilíndrico, que no se alcanza en un pico
            # anguloso (el arco cilíndrico se atasca con pasos que apenas
            # disipan).
            if self._mode == "cylindrical" and du_L is not None:
                alpha_L = 0.5 * (
                    lambda_curr * float(F_ext_ref @ du_L)
                    - float(F_ext_ref @ U_current)
                )
                scale_L = abs(0.5 * lambda_curr * float(F_ext_ref @ du_L))
                if scale_L > 0 and abs(alpha_L) > ALPHA_SWITCH_RTOL * scale_L:
                    _log.info(
                        f"  Switch cilíndrico→disipación "
                        f"(|α_L|/escala = {abs(alpha_L) / scale_L:.3e} > "
                        f"{ALPHA_SWITCH_RTOL:g})"
                    )
                    self._mode = "dissipation"
                    self._tau = self.initial_tau

            # dlambda según modo activo.
            mode_predictor = self._mode
            du_pred = du_t
            if mode_predictor == "dissipation":
                # Predictor sobre la tangente de carga (criterio (a)); sin K_L
                # utilizable, la del estado consolidado, que da α ≈ 0 y
                # devuelve el paso a cilíndrico abajo.
                if du_L is not None:
                    du_pred = du_L
                # α = ½·(λ_n·F·δu − F·U_n)
                alpha = 0.5 * (
                    lambda_curr * float(F_ext_ref @ du_pred)
                    - float(F_ext_ref @ U_current)
                )
                # Para problemas lineales monotónicos, α ≡ 0 exactamente
                # (U_n = λ_n·K^-1·F_ref ⇒ F·U_n = λ_n·F·K^-1·F = λ_n·F·δu).
                # Por eso necesitamos un threshold **relativo** a la escala
                # energética del paso: si |α| ≪ ½·λ_n·F·δu (= energía
                # elástica del estado escalada), el modo disipación no es
                # físicamente significativo y se revierte a cilíndrico.
                F_du_pred = float(F_ext_ref @ du_pred)
                alpha_scale = abs(0.5 * lambda_curr * F_du_pred)
                if alpha_scale > 0 and abs(alpha) < ALPHA_LINEAR_RTOL * alpha_scale:
                    _log.info(
                        f"  α/escala = {abs(alpha)/alpha_scale:.2e} ≪ 1 "
                        f"⇒ régimen lineal: cilíndrico temporal."
                    )
                    mode_predictor = "cylindrical"
                    # Resetear el modo del solver también para evitar
                    # re-entrar en disipación en el siguiente paso sin
                    # disipación física genuina.
                    self._mode = "cylindrical"
                elif _alpha_is_null(alpha, lambda_curr, F_du_pred,
                                    float(F_ext_ref @ U_current)):
                    _log.info("  α ≈ 0 estrictamente → cilíndrico temporal.")
                    mode_predictor = "cylindrical"
                    self._mode = "cylindrical"

            if mode_predictor == "cylindrical":
                dlambda = sign * dl / (np.linalg.norm(du_t) + ZERO_TOL)
            else:  # dissipation
                # La restricción es lineal: con ΔU = Δλ·δu, g = Δλ·α. Exigir
                # g = τ > 0 (disipación positiva) fija el valor Y el signo de
                # Δλ, sin regla de signo (Gutiérrez 2004; Verhoosel et al.
                # 2009). Con α < 0 el paso reduce la carga: es lo que permite
                # seguir un retroceso. Hasta el 2026-09-29 el signo venía de la
                # regla del arco cilíndrico (sign·|τ/α|) y, cuando discrepaba
                # del de α, el predictor arrancaba con g = −τ.
                dlambda = self._tau / alpha

            # Llegada a max_lambda (2026-09-23): sin paso final decidido por el
            # predictor. Todos los pasos llevan su restricción; el que cruza
            # max_lambda no se consolida y se repite con llegada exacta dentro
            # del tramo recorrido (ArcLengthSolver._land). Sustituye a la
            # salvaguarda anterior, que bisecaba cuando el predictor rebasaba
            # max_lambda más de 3 veces lo que faltaba y, por debajo de ese
            # factor, fijaba λ = max_lambda en control de carga.
            dU_iter = dlambda * du_pred
            x0 = (U_current + dU_iter, lambda_curr + dlambda, dU_iter)

            # --- 2. CORRECTOR ITERATIVO ---
            # Mismo corrector compartido que el padre; sólo cambia la
            # restricción del modo activo.
            problem = _DissipationArcProblem(
                self.assembler, U_current, lambda_curr, F_ext_ref, free_dofs,
                mode=mode_predictor, dl=dl, tau=self._tau, max_lambda=self.max_lambda,
            )
            res = self.corrector.run(
                problem, x0, check_initial=True,
                initial_delta_norm=float(np.linalg.norm(dU_iter)),
            )

            if res.converged and problem.overshoot:
                landed = self._land(U_current, lambda_curr, res.x, F_ext_ref, free_dofs)
                if landed.converged:
                    U_current, lambda_curr, delta_U_step = landed.x
                    _log.info(f"  -> CONV. lam={lambda_curr:.4f} (llegada a max_lambda)")
                    steps_done += 1
                    if step_callback is not None:
                        step_callback(step, U_current, lambda_curr)
                    continue
                res = landed          # llegada fallida: bisección del modo del paso

            if res.converged:
                U_iter, lambda_iter, dU_iter = res.x
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

                # Tangente de carga del paso siguiente: la del ensamblaje en el
                # iterado convergido (CorrectorResult.state = (K, F_int)),
                # evaluada con la prueba de carga respecto del estado
                # consolidado anterior. El ensamblador crea una matriz nueva
                # en cada ensamblaje, así que guardarla no la expone a
                # sobrescrituras.
                self._K_load = res.state[0]

                # Vuelta a cilíndrico tras commit (sin cambios): umbral relativo
                # a la escala energética del estado, ‖F_ref‖·‖U‖, sin piso
                # absoluto (el max(·, 1) anterior, hasta el 2026-09-29, estaba
                # en unidades de energía; ADR 0006). El arranque del modo de
                # disipación ya no se decide aquí con ΔE_d > umbral, sino en el
                # predictor con α_L (2026-09-30).
                ref_energy = np.linalg.norm(F_ext_ref) * np.linalg.norm(U_iter)
                dE_threshold = self.dissipation_threshold * ref_energy
                if self._mode == "dissipation" and dE_d < dE_threshold:
                    _log.info(
                        f"  Switch disipación→cilíndrico "
                        f"(ΔE_d={dE_d:.3e} < {dE_threshold:.3e})"
                    )
                    self._mode = "cylindrical"

                # Adaptatividad según el modo del **siguiente** paso.
                n_solves = res.n_solves
                if self._mode == "cylindrical":
                    if n_solves < self.dl_grow_iter_threshold:
                        dl = min(
                            dl * self.dl_grow_factor,
                            dl_ref * self.dl_max_factor,
                        )
                    elif n_solves > self.dl_shrink_iter_threshold:
                        dl *= self.dl_shrink_factor
                else:
                    if n_solves < self.tau_grow_iter_threshold:
                        self._tau = min(
                            self._tau * self.tau_grow_factor,
                            self.initial_tau * self.tau_max_factor,
                        )
                    elif n_solves > self.tau_shrink_iter_threshold:
                        self._tau *= self.tau_shrink_factor

                U_current = U_iter
                lambda_curr = lambda_iter
                delta_U_step = dU_iter
                steps_done += 1

                if step_callback is not None:
                    step_callback(step, U_current, lambda_curr)
                continue

            # Bisección por modo del paso fallido.
            if mode_predictor == "cylindrical":
                dl *= 0.5
                _log.warning(
                    f"Bisección cilíndrica: dl → {dl:.4e}"
                )
                if dl < ARCLENGTH_MIN_DL_FACTOR * dl_ref:
                    raise RuntimeError(
                        "DissipationArcLengthSolver: dl bajó del umbral "
                        f"({ARCLENGTH_MIN_DL_FACTOR * dl_ref:.4e}). "
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

        self._finish(lambda_curr, steps_done)
        return U_current
