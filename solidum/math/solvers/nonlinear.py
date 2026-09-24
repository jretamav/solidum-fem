# solidum_fem/solidum/math/solvers/nonlinear.py
"""``NonlinearSolver`` — Newton-Raphson incremental con paso adaptativo.

Estructura del bucle de Newton (un ensamblaje por iteración)
-------------------------------------------------------------
Cada iteración ensambla **una sola vez** en el iterado corriente ``U_k``:
de ese ensamblaje salen la tangente ``K_t(U_k)`` para el sistema lineal,
el residuo ``R(U_k)`` para el criterio de convergencia y el estado
*trial* de los elementos. La convergencia se evalúa **al inicio** de la
iteración con el par ``(‖R(U_k)‖, ‖δU_{k−1}‖)``: si se cumple, el estado
trial que acaba de dejar el ensamblaje es exactamente el de ``U_k`` y se
comitea sin reensamblar. Si no, se resuelve ``δU_k`` y se avanza.

Ese bucle es el :class:`~solidum.math.solvers.corrector.NewtonCorrector`
compartido (ADR 0015); este módulo aporta lo propio del control de carga
incremental: el residuo ``λ·F_ext − F_int``, la reducción por Dirichlet
con el incremento ``λ·g`` del paso, y el paso adaptativo con bisección.

La alternativa —ensamblar en ``U_k`` para resolver y otra vez en
``U_k + δU_k`` para evaluar el residuo— duplicaba el coste por iteración
(auditoría 2026-09-22) para producir la misma secuencia de decisiones,
sólo que desplazada una iteración. Con line search activo (opt-in) el
backtracking ensambla en cada ``α`` probado y devuelve ``K`` y ``F_int``
del punto aceptado, así que tampoco hay reensamblaje redundante.
"""
from __future__ import annotations

import numpy as np

from solidum.constants import (
    NEWTON_ADAPTIVE_GROWTH_FACTOR,
    NEWTON_ADAPTIVE_GROWTH_ITER_THRESHOLD,
    NEWTON_DEFAULT_MIN_DELTA_LAMBDA,
    NEWTON_LOAD_FACTOR_EPSILON,
)
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers._shared import _log, domain_is_symmetric
from solidum.math.solvers.corrector import NewtonCorrector, default_calibration_scales
from solidum.math.solvers.model_checks import ensure_statically_restrained
from solidum.registry import SolverRegistry


class _IncrementalProblem:
    """Paso de carga del Newton incremental, en el protocolo ``NewtonProblem``.

    El iterado ``x`` es el vector global ``U``; el estado es ``(K, F_int)``.
    """

    def __init__(self, assembler, F_ext_step: np.ndarray, F_ext_global: np.ndarray,
                 load_factor: float, free_dofs: np.ndarray):
        self.assembler = assembler
        self.F_ext_step = F_ext_step
        self.F_ext_global = F_ext_global
        self.load_factor = float(load_factor)
        self.free_dofs = free_dofs

    def assemble(self, U):
        return self.assembler.assemble_non_linear_system(U)

    def residual(self, U, state):
        return self.F_ext_step - state[1]

    def residual_norm(self, R):
        return float(np.linalg.norm(R[self.free_dofs]))

    def calibration_scales(self, U, state):
        # La escala de fuerza es la carga total de la corrida, no la del
        # paso: así las tolerancias no dependen del número de pasos.
        return default_calibration_scales(self.F_ext_global, state[1], state[0])

    def reference_force(self, U, state):
        return max(float(np.linalg.norm(self.F_ext_step)), float(np.linalg.norm(state[1])))

    def x_norm(self, U):
        return float(np.linalg.norm(U))

    def correction(self, U, state, R, solve):
        # Reducción con el incremento de Dirichlet ``λ·g`` del paso: en DOF
        # libres se anula y en esclavos corrige la restricción (ADR 0004).
        K_red, R_red, T_op, g_inc = self.assembler.reduce(
            state[0], R, U_current=U, load_factor=self.load_factor,
        )
        return self.assembler.expand(solve(K_red, R_red), T_op, g_inc)

    def apply(self, U, dU, alpha):
        step = alpha * dU
        return U + step, float(np.linalg.norm(step))

    def on_converged(self, U, state):
        # El estado trial es el del ensamblaje en U: coherente.
        self.assembler.commit_all_states()


@SolverRegistry.register
class NonlinearSolver:
    """Solucionador incremental-iterativo de Newton-Raphson con paso adaptativo."""

    PIPELINE_KIND = "static"

    def __init__(
        self,
        assembler,
        convergence: ConvergenceCriterion | None = None,
        *,
        max_iter: int = 20,
        num_steps: int = 10,
        adaptive: bool = True,
        min_delta_lambda: float = NEWTON_DEFAULT_MIN_DELTA_LAMBDA,
        linear_algebra: str = "auto",
        freeze_tangent_after_iter: int | None = None,
        line_search: bool = False,
    ):
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
        # Corrector compartido (ADR 0015): bucle de Newton, backend
        # algebraico con degradación a LU, Newton modificado y line search.
        self.corrector = NewtonCorrector(
            self.convergence,
            max_iter=self.max_iter,
            is_symmetric=domain_is_symmetric(assembler.domain),
            is_positive_definite=True,
            linear_algebra=self.linear_algebra,
            freeze_tangent_after_iter=self.freeze_tangent_after_iter,
            line_search=self.line_search,
            near_nullspace=assembler.near_nullspace,
        )
        # Metadatos del último análisis (los lee ``solidum.run``): pasos
        # convergidos y factor de carga alcanzado (siempre 1.0 si ``solve``
        # retorna; en caso contrario lanza una excepción tipada).
        self.steps_done: int = 0
        self.lambda_final: float = 0.0
        self.reached_max_lambda: bool = False

    def make_problem(self, F_ext_global: np.ndarray, load_factor: float) -> _IncrementalProblem:
        """Problema de Newton del paso con factor de carga ``load_factor``
        (lo usa ``solve``; expuesto para los tests de los internos)."""
        ndof = self.assembler.domain.total_dofs
        free_dofs = self.assembler.constraint_set.free_dofs(ndof)
        return _IncrementalProblem(
            self.assembler, F_ext_global * load_factor, F_ext_global, load_factor, free_dofs,
        )

    def solve(self, F_ext_global: np.ndarray, step_callback=None) -> np.ndarray:
        domain = self.assembler.domain
        ndof = domain.total_dofs
        U_current = np.zeros(ndof)

        _log.info("--- INICIANDO SOLVER NO LINEAL (CONTROL DE PASO ADAPTATIVO) ---")
        # Red de seguridad, capa 1 (ADR 0019): un mecanismo rígido no tiene
        # equilibrio estático; se rechaza antes de iterar con el movimiento
        # libre en palabras, en vez de agotar bisecciones.
        ensure_statically_restrained(self.assembler, type(self).__name__)

        load_factor = 0.0
        target_load = 1.0
        delta_lambda = 1.0 / self.num_steps
        step = 0
        n_bisections = 0  # contador de bisecciones globales del paso (ADR 0011)
        self.steps_done = 0
        self.lambda_final = 0.0
        self.reached_max_lambda = False

        while load_factor < target_load - NEWTON_LOAD_FACTOR_EPSILON:
            step += 1

            if load_factor + delta_lambda > target_load:
                delta_lambda = target_load - load_factor

            next_load_factor = load_factor + delta_lambda
            _log.info(f"[PASO {step}] Intentando Factor de Carga: {next_load_factor:.4f} (Incremento: {delta_lambda:.4f})")

            # Gancho de inicio de paso (ADR 0020, P5). Evaluado con el estado
            # convergido del paso anterior para evitar chattering por
            # predictores lineales dentro del Newton. No-op para elementos
            # estándar; lo usan los que toman decisiones discretas entre pasos.
            self.assembler.prepare_all_steps(U_current)

            # Corrector del paso. No se evalúa la convergencia en la
            # iteración 0: el iterado inicial es el convergido del paso
            # anterior y aún no incorpora el incremento de Dirichlet ``λ·g``
            # del paso nuevo (entra por ``reduce`` en la primera
            # resolución). Con control en desplazamiento su residuo en DOF
            # libres es nulo y el criterio daría convergencia sin haber
            # movido el apoyo. Todo paso hace, por tanto, al menos una
            # resolución.
            problem = self.make_problem(F_ext_global, next_load_factor)
            res = self.corrector.run(problem, U_current.copy(), check_initial=False)

            if res.converged:
                _log.info("  -> CONVERGENCIA ALCANZADA.")
                U_current = res.x
                load_factor = next_load_factor
                self.steps_done += 1

                if (self.adaptive
                        and res.n_solves < NEWTON_ADAPTIVE_GROWTH_ITER_THRESHOLD
                        and delta_lambda < (1.0 / self.num_steps)):
                    delta_lambda = min(
                        delta_lambda * NEWTON_ADAPTIVE_GROWTH_FACTOR,
                        1.0 / self.num_steps,
                    )
                    _log.info(f"  -> Acelerando el próximo incremento a {delta_lambda:.4f}")

                if step_callback:
                    step_callback(step, U_current, load_factor)
                continue

            if self.adaptive:
                delta_lambda /= 2.0
                n_bisections += 1
                _log.warning(f"NO CONVERGIÓ. Bisección: reduciendo incremento a {delta_lambda:.4f}")
                if delta_lambda < self.min_delta_lambda:
                    # Clasificar el modo y lanzar excepción tipada (ADR 0011).
                    raise res.divergence_error(
                        last_load_factor=next_load_factor,
                        n_bisections=n_bisections,
                        extra_message=(
                            f"Δλ={delta_lambda:.2e} < min={self.min_delta_lambda:.2e}; "
                            f"line_search={'on' if self.line_search else 'off'}; "
                            f"último α={res.last_alpha:.3f}."
                        ),
                    )
            else:
                raise res.divergence_error(
                    last_load_factor=next_load_factor,
                    n_bisections=n_bisections,
                    extra_message=(
                        f"Paso {step} no convergió en {self.max_iter} iter; "
                        f"adaptive=False; line_search={'on' if self.line_search else 'off'}."
                    ),
                )

        self.lambda_final = load_factor
        self.reached_max_lambda = True
        return U_current
