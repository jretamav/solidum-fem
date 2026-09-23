# solidum_fem/solidum/math/solvers/arclength.py
"""``ArcLengthSolver`` — método de longitud de arco cilíndrico (Crisfield).

Estructura del corrector (un ensamblaje por iteración, estado coherente)
------------------------------------------------------------------------
Cada iteración del corrector ensambla una sola vez en el iterado corriente
``U_k`` y evalúa la convergencia **antes** de resolver, con el par
``(‖R(U_k)‖, ‖δU_{k−1}‖)``. Al converger, el estado trial de los elementos
es exactamente el del ensamblaje en ``U_k`` y se comitea junto con
``U_current = U_k`` y ``λ_curr = λ_k``. La versión anterior evaluaba el
residuo en ``U_k`` pero almacenaba ``U_{k+1}``: el estado committed iba
un iterado por detrás del desplazamiento guardado (auditoría 2026-09-22).

El bucle es el :class:`~solidum.math.solvers.corrector.NewtonCorrector`
compartido (ADR 0015). Lo propio del arc-length —predictor tangente,
las dos resoluciones por iteración (``du_R`` y ``du_t``) y la restricción
que fija ``ddλ``— vive en :class:`_ArcProblem`; la variante por
disipación sólo redefine la restricción.

Fin del trazado
---------------
El bucle termina al alcanzar ``max_lambda`` o al agotar ``max_steps``. En
el segundo caso el solver **no** lanza excepción —trazar "hasta donde se
pueda" es un uso legítimo del arc-length— pero lo registra con un
``WARNING`` y lo expone en ``reached_max_lambda``/``lambda_final`` para
que ``solidum.run`` pueda reflejarlo en ``SolveResult.converged`` y escalar
``F_applied`` por el factor de carga realmente alcanzado.
"""
from __future__ import annotations

import numpy as np

from solidum.constants import ARCLENGTH_MIN_DL_FACTOR, ZERO_TOL
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers._shared import _log, domain_is_symmetric
from solidum.math.solvers.corrector import (
    CorrectionAborted,
    NewtonCorrector,
    default_calibration_scales,
)
from solidum.registry import SolverRegistry


class _ArcProblem:
    """Paso del arc-length en el protocolo ``NewtonProblem``.

    El iterado ``x`` es ``(U_iter, λ_iter, dU_iter)`` con ``dU_iter`` el
    incremento acumulado del paso (``U_iter = U_current + dU_iter``); el
    estado es ``(K, F_int)``. La corrección ``dx`` es ``(dU_update, ddλ)``.

    ``mode`` fija la restricción que determina ``ddλ``:

    - ``"newton"`` — último paso: ``λ`` fijo, corrección de Newton pura;
    - ``"cylindrical"`` — cuadrática de Crisfield ``‖ΔU‖² = dl²`` con
      selección de raíz por menor ángulo con el incremento previo.
    """

    def __init__(self, assembler, U_current, lambda_curr, F_ext_ref, free_dofs,
                 *, mode: str, dl: float):
        self.assembler = assembler
        self.U_current = U_current
        self.lambda_curr = float(lambda_curr)
        self.F_ext_ref = F_ext_ref
        self.free_dofs = free_dofs
        self.mode = mode
        self.dl = float(dl)

    # --- protocolo ----------------------------------------------------------

    def assemble(self, x):
        return self.assembler.assemble_non_linear_system(x[0])

    def residual(self, x, state):
        return x[1] * self.F_ext_ref - state[1]

    def residual_norm(self, R):
        return float(np.linalg.norm(R[self.free_dofs]))

    def calibration_scales(self, x, state):
        # La carga de referencia es la escala natural del arc-length (λ
        # varía, F_ext_ref es fijo).
        return default_calibration_scales(self.F_ext_ref, state[1], state[0])

    def reference_force(self, x, state):
        return max(float(np.linalg.norm(self.F_ext_ref)) * abs(x[1]),
                   float(np.linalg.norm(state[1])))

    def x_norm(self, x):
        return float(np.linalg.norm(x[0]))

    def log_context(self, x):
        return f" | lam={x[1]:.4f}"

    def correction(self, x, state, R, solve):
        U_iter, lambda_iter, dU_iter = x
        K_global = state[0]
        K_t_red, F_t_red, T_t, g_t = self.assembler.reduce(K_global, self.F_ext_ref.copy())
        K_red, R_red, T_R, g_R = self.assembler.reduce(
            K_global, R, U_current=U_iter, load_factor=lambda_iter,
        )
        du_R = self.assembler.expand(solve(K_red, R_red), T_R, g_R)
        du_t = self.assembler.expand(solve(K_t_red, F_t_red), T_t, g_t)
        return self.constraint(x, du_R, du_t)

    def apply(self, x, dx, alpha):
        U_iter, lambda_iter, dU_iter = x
        dU_update, ddlambda = dx
        dU_new = dU_iter + alpha * dU_update
        return ((self.U_current + dU_new, lambda_iter + alpha * ddlambda, dU_new),
                float(np.linalg.norm(alpha * dU_update)))

    def on_converged(self, x, state):
        # Estado trial del ensamblaje en U_iter: coherente con lo guardado.
        self.assembler.commit_all_states()

    # --- restricción --------------------------------------------------------

    def constraint(self, x, du_R, du_t):
        """``(dU_update, ddλ)`` según ``mode``."""
        _, _, dU_iter = x
        if self.mode == "newton":
            return du_R, 0.0
        if self.mode == "cylindrical":
            dU_new = dU_iter + du_R
            a = np.dot(du_t, du_t)
            b = 2.0 * np.dot(dU_new, du_t)
            c = np.dot(dU_new, dU_new) - self.dl ** 2
            det = b ** 2 - 4.0 * a * c
            if det < 0:
                raise CorrectionAborted("Raíces imaginarias. La solución diverge del arco.")
            ddl1 = (-b + np.sqrt(det)) / (2.0 * a)
            ddl2 = (-b - np.sqrt(det)) / (2.0 * a)
            # Raíz que produce el menor ángulo con el incremento previo.
            theta1 = np.dot(dU_iter, dU_new + ddl1 * du_t)
            theta2 = np.dot(dU_iter, dU_new + ddl2 * du_t)
            ddlambda = ddl1 if theta1 > theta2 else ddl2
            return du_R + ddlambda * du_t, ddlambda
        raise ValueError(f"_ArcProblem: modo desconocido {self.mode!r}.")


@SolverRegistry.register
class ArcLengthSolver:
    """
    Solucionador no lineal con Método de Longitud de Arco Cilíndrico (Crisfield).
    Permite trazar curvas de equilibrio con fenómenos de snap-through y snap-back
    variando simultáneamente los desplazamientos y la carga externa.

    Atributos de estado tras ``solve``
    ----------------------------------
    lambda_final : float
        Factor de carga del último paso convergido.
    reached_max_lambda : bool
        ``True`` si el trazado alcanzó ``max_lambda``; ``False`` si se detuvo
        antes por agotar ``max_steps``.
    steps_done : int
        Número de pasos convergidos.
    """

    PIPELINE_KIND = "static"

    def __init__(self, assembler, convergence: ConvergenceCriterion | None = None, max_iter=20, max_lambda=1.0, initial_dl=0.1, max_steps=100,
                 dl_grow_factor=1.5, dl_max_factor=5.0, dl_shrink_factor=0.6,
                 dl_grow_iter_threshold=4, dl_shrink_iter_threshold=8,
                 linear_algebra: str = "auto"):
        self.assembler = assembler
        # Política de convergencia (ADR 0007). Compartida con NonlinearSolver:
        # cambiar la política aquí llega automáticamente al arc-length.
        self.convergence = convergence if convergence is not None else ConvergenceCriterion()
        self.max_iter = max_iter
        self.max_lambda = max_lambda
        self.dl = initial_dl
        self.max_steps = max_steps
        # Factores de auto-ajuste de la longitud de arco:
        #   Si converge en < dl_grow_iter_threshold iter → ampliar dl × dl_grow_factor (max: initial_dl × dl_max_factor)
        #   Si converge en > dl_shrink_iter_threshold iter → reducir dl × dl_shrink_factor
        self.dl_grow_factor = dl_grow_factor
        self.dl_max_factor = dl_max_factor
        self.dl_shrink_factor = dl_shrink_factor
        self.dl_grow_iter_threshold = dl_grow_iter_threshold
        self.dl_shrink_iter_threshold = dl_shrink_iter_threshold
        self.linear_algebra = linear_algebra
        # Corrector compartido (ADR 0015). Régimen postcrítico: K_t puede
        # ser indefinida → no asumir SPD. Si el usuario fuerza
        # ``linear_algebra: cholesky`` y K_t se vuelve indefinida, el
        # corrector degrada a LU para que el override no rompa el análisis.
        self.corrector = NewtonCorrector(
            self.convergence,
            max_iter=self.max_iter,
            is_symmetric=domain_is_symmetric(assembler.domain),
            is_positive_definite=False,
            linear_algebra=self.linear_algebra,
        )
        # Estado del último trazado (ver docstring de la clase).
        self.lambda_final: float = 0.0
        self.reached_max_lambda: bool = False
        self.steps_done: int = 0

    def _negative_pivots(self, K) -> int | None:
        """Diagnóstico de bifurcación vía Sturm sequence (ADR 0003 fase 2).

        Cuando ``LDLTSolver`` esté implementado (backend con LDLᵀ verdadero
        de Bunch-Kaufman, p. ej. pypardiso), este método factoriza ``K_t`` y
        retorna ``factor.n_negative_pivots``. Por la ley de inercia de
        Sylvester equivale al número de autovalores negativos: 0 antes de
        bifurcación, 1 tras el primer punto crítico, etc.

        En fase 2 retorna ``None`` (LDLᵀ es placeholder; LU general no
        preserva inercia y dar un conteo aproximado sería engañoso).
        """
        return None

    def _finish(self, lambda_curr: float, step: int) -> None:
        """Registra el estado final del trazado y avisa si quedó incompleto."""
        self.lambda_final = float(lambda_curr)
        self.steps_done = int(step)
        self.reached_max_lambda = lambda_curr >= self.max_lambda - ZERO_TOL
        if not self.reached_max_lambda:
            _log.warning(
                f"{type(self).__name__}: trazado detenido en λ={lambda_curr:.4f} "
                f"< max_lambda={self.max_lambda:.4f} tras agotar max_steps="
                f"{self.max_steps}. El resultado corresponde al último paso "
                f"convergido; aumente max_steps o initial_dl para completar."
            )

    # ------------------------------------------------------------------
    # Predictor tangente (común a las variantes)
    # ------------------------------------------------------------------

    def _tangent_predictor(self, U_current, F_ext_ref, delta_U_step, step):
        """Ensambla en ``U_current``, calibra el criterio si hace falta y
        devuelve ``(du_t, sign)`` con ``du_t = K⁻¹·F_ref`` y el sentido de
        avance (para no regresar por donde se vino). Devuelve ``None`` si la
        tangente es singular (el llamador biseca la longitud de paso)."""
        K_global, F_int_global = self.assembler.assemble_non_linear_system(U_current)
        if not self.convergence.is_calibrated:
            self.convergence.calibrate(
                *default_calibration_scales(F_ext_ref, F_int_global, K_global))
        K_t_red, F_t_red, T_t, g_t = self.assembler.reduce(K_global, F_ext_ref.copy())
        try:
            du_t_red = self.corrector.solve(K_t_red, F_t_red)
        except RuntimeError:
            return None
        du_t = self.assembler.expand(du_t_red, T_t, g_t)
        sign = 1.0
        if step > 1 and np.dot(delta_U_step, du_t) < 0:
            sign = -1.0
        return du_t, sign

    def solve(self, F_ext_ref: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U_current = np.zeros(ndof)

        lambda_curr = 0.0
        step = 0
        steps_done = 0
        dl = self.dl

        delta_U_step = np.zeros(ndof)  # Historial del incremento del paso para guiar el arco

        _log.info("--- INICIANDO SOLVER NO LINEAL (MÉTODO ARC-LENGTH) ---")

        cs = self.assembler.constraint_set
        free_dofs = cs.free_dofs(ndof)

        while lambda_curr < self.max_lambda and step < self.max_steps:
            step += 1
            _log.info(f"[PASO {step}] Longitud de Arco (dl): {dl:.4e}")

            # ADR 0010 §5: hook de preparación de paso (activación de
            # discontinuidades embebidas, etc.). Evaluado con el estado
            # convergido del paso anterior para evitar chattering.
            self.assembler.prepare_all_steps(U_current)

            # --- 1. PREDICTOR ---
            pred = self._tangent_predictor(U_current, F_ext_ref, delta_U_step, step)
            if pred is None:
                _log.error("Matriz singular en predictor. Bisección de dl...")
                dl /= 2.0
                continue
            du_t, sign = pred

            dlambda = sign * dl / (np.linalg.norm(du_t) + ZERO_TOL)

            # Si el paso predictor sobrepasaría max_lambda, fijar lambda exactamente
            final_step = (sign > 0 and lambda_curr + dlambda >= self.max_lambda - ZERO_TOL)
            if final_step:
                dlambda = self.max_lambda - lambda_curr

            dU_iter = dlambda * du_t
            x0 = (U_current + dU_iter, lambda_curr + dlambda, dU_iter)

            # --- 2. CORRECTOR ITERATIVO ---
            # Último paso: lambda fijo, solo corrección de desplazamientos
            # (Newton-Raphson puro); en el resto, restricción cilíndrica.
            problem = _ArcProblem(
                self.assembler, U_current, lambda_curr, F_ext_ref, free_dofs,
                mode="newton" if final_step else "cylindrical", dl=dl,
            )
            res = self.corrector.run(
                problem, x0, check_initial=True,
                initial_delta_norm=float(np.linalg.norm(dU_iter)),
            )

            if res.converged:
                U_current, lambda_curr, delta_U_step = res.x
                _log.info(f"  -> CONVERGENCIA. (Lambda alcanzado: {lambda_curr:.4f})")
                steps_done += 1
                # Auto-ajuste de longitud de arco
                if res.n_solves < self.dl_grow_iter_threshold:
                    dl = min(dl * self.dl_grow_factor, self.dl * self.dl_max_factor)
                elif res.n_solves > self.dl_shrink_iter_threshold:
                    dl *= self.dl_shrink_factor

                if step_callback:
                    step_callback(step, U_current, lambda_curr)
                continue

            dl *= 0.5
            _log.warning(f"Bisección: reduciendo longitud de arco a {dl:.4e}")
            if dl < ARCLENGTH_MIN_DL_FACTOR * self.dl:
                raise RuntimeError("Arc-Length fracasó irreparablemente.")

        self._finish(lambda_curr, steps_done)
        return U_current
