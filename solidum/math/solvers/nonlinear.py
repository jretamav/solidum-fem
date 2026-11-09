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

Etapas de carga (ADR 0026): ``solve(F, u0=..., constant_load=...)`` parte del
estado ``u0`` con la carga constante ``F_c`` ya aplicada; el residuo es
``F_c + λ·F − F_int`` con λ de 0 a 1, y los valores impuestos crecen de los de
``u0`` a los del modelo (:mod:`solidum.math.solvers.stages`).

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
from solidum.math.solvers.diagnostics import Termination
from solidum.math.solvers.stages import Stage, check_initial_equilibrium, make_stage
from solidum.math.solvers.model_checks import (
    check_initial_tangent,
    ensure_statically_restrained,
)
from solidum.registry import SolverRegistry


class _IncrementalProblem:
    """Paso de carga del Newton incremental, en el protocolo ``NewtonProblem``.

    El iterado ``x`` es el vector global ``U``; el estado es ``(K, F_int)``.
    """

    def __init__(self, assembler, F_ext_step: np.ndarray, F_ext_global: np.ndarray,
                 load_factor: float, free_dofs: np.ndarray, *,
                 g_start: np.ndarray | None = None):
        self.assembler = assembler
        # ``F_ext_step`` es la carga total del paso, constante incluida
        # (``F_c + λ·F``); ``F_ext_global`` la escala de calibración.
        self.F_ext_step = F_ext_step
        self.F_ext_global = F_ext_global
        self.load_factor = float(load_factor)
        self.free_dofs = free_dofs
        # Valores impuestos al empezar la etapa (ADR 0026); ``None`` = cero.
        self.g_start = g_start
        # Operador ``u = T·u_libre + g`` (ADR 0004): el equilibrio reducido es
        # ``TᵀR = 0``. Con una restricción multipunto la fuerza del esclavo
        # pasa a sus maestros y ``R[libres]`` no se anula en la solución.
        self._T = assembler.constraint_set.build(assembler.domain.total_dofs)[0]

    def assemble(self, U):
        return self.assembler.assemble_non_linear_system(U)

    def residual(self, U, state):
        return self.F_ext_step - state[1]

    def residual_norm(self, R):
        return float(np.linalg.norm(self._T.T @ R))

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
            state[0], R, U_current=U, load_factor=self.load_factor, g_start=self.g_start,
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
        # Metadatos del último análisis (los lee ``solidum.run``; ADR 0025):
        # motivo de la parada, pasos convergidos, factor de carga alcanzado
        # (1.0 si ``solve`` retorna; si lanza una excepción tipada, el del
        # último paso convergido) y su estado de desplazamientos. ``solve``
        # los reinicia al empezar.
        self.termination: Termination | None = None
        self.steps_done: int = 0
        self.lambda_final: float = 0.0
        self.reached_max_lambda: bool = False
        self.n_bisections: int = 0
        self.U_final: np.ndarray | None = None
        # Residuo relativo de equilibrio al empezar una etapa (ADR 0026);
        # ``None`` si no había nada que comprobar.
        self.initial_residual: float | None = None

    @property
    def completed(self) -> bool:
        """``True`` si el análisis llegó a la carga completa (λ = 1)."""
        return self.termination is Termination.MAX_LAMBDA

    def _fail(self, error, load_factor: float, U: np.ndarray):
        """Estado de un análisis que fracasa (``termination = "diverged"``,
        último paso convergido) y la excepción tipada, lista para lanzar."""
        self.termination = Termination.DIVERGED
        self.lambda_final = float(load_factor)
        self.U_final = np.array(U, dtype=float, copy=True)
        _log.error(
            f"{type(self).__name__}: análisis detenido ({error.mode}). Último paso "
            f"convergido: λ = {load_factor:.4f} tras {self.steps_done} pasos."
        )
        return error

    def make_problem(self, F_ext_global: np.ndarray, load_factor: float,
                     stage: Stage | None = None) -> _IncrementalProblem:
        """Problema de Newton del paso con factor de carga ``load_factor``
        (lo usa ``solve``; expuesto para los tests de los internos). Con una
        etapa (ADR 0026), la carga del paso incluye la constante y los valores
        impuestos parten de los de la etapa."""
        ndof = self.assembler.domain.total_dofs
        free_dofs = self.assembler.constraint_set.free_dofs(ndof)
        if stage is None or not stage.staged:
            return _IncrementalProblem(
                self.assembler, F_ext_global * load_factor, F_ext_global, load_factor, free_dofs,
            )
        return _IncrementalProblem(
            self.assembler, stage.F_const + F_ext_global * load_factor,
            stage.calibration_load(F_ext_global), load_factor, free_dofs,
            g_start=stage.g_start,
        )

    def solve(self, F_ext_global: np.ndarray, step_callback=None, *,
              u0: np.ndarray | None = None,
              constant_load: np.ndarray | None = None) -> np.ndarray:
        """Aplica ``F_ext_global`` de forma incremental (λ de 0 a 1).

        Parameters
        ----------
        F_ext_global
            Carga que crece con λ.
        step_callback
            ``callback(step, U, λ)`` tras cada paso convergido.
        u0, constant_load
            Etapa de carga (ADR 0026): estado de partida (los desplazamientos
            de la etapa anterior; el estado interno de los elementos lo guarda
            el ensamblador) y carga que actúa completa durante toda la etapa
            (la ya aplicada). El residuo es ``F_c + λ·F − F_int``.
        """
        # Estado del análisis anterior fuera (ADR 0025), antes de cualquier
        # comprobación que pueda lanzar.
        self.termination = None
        self.steps_done = 0
        self.lambda_final = 0.0
        self.reached_max_lambda = False
        self.n_bisections = 0
        self.U_final = None
        self.initial_residual = None
        domain = self.assembler.domain
        if domain.total_dofs == 0:
            domain.generate_equation_numbers()
        stage = make_stage(self.assembler, u0, constant_load, type(self).__name__)
        U_current = stage.U0.copy()

        _log.info("--- INICIANDO SOLVER NO LINEAL (CONTROL DE PASO ADAPTATIVO) ---")
        # Red de seguridad, capa 1 (ADR 0019): un mecanismo rígido no tiene
        # equilibrio estático; se rechaza antes de iterar con el movimiento
        # libre en palabras, en vez de agotar bisecciones.
        ensure_statically_restrained(self.assembler, type(self).__name__)
        # Capa 3 sobre la tangente inicial (ADR 0019, 2026-09-30): la primera
        # iteración factoriza K(U = 0); sus pivotes nulos delatan un mecanismo
        # interno aunque la carga no lo active (el Newton lo resolvería en
        # silencio con los pivotes perturbados de Pardiso).
        check_initial_tangent(self.assembler, self.corrector, type(self).__name__)
        # Etapa (ADR 0026): debe empezar en equilibrio con la carga constante.
        self.initial_residual = check_initial_equilibrium(
            self.assembler, stage, self.convergence, F_ext_global, type(self).__name__)

        load_factor = 0.0
        target_load = 1.0
        delta_lambda = 1.0 / self.num_steps
        step = 0
        n_bisections = 0  # contador de bisecciones globales del paso (ADR 0011)

        while load_factor < target_load - NEWTON_LOAD_FACTOR_EPSILON:
            step += 1

            if load_factor + delta_lambda > target_load:
                delta_lambda = target_load - load_factor

            next_load_factor = load_factor + delta_lambda
            _log.debug(f"[PASO {step}] intentando λ = {next_load_factor:.4f} "
                       f"(Δλ = {delta_lambda:.4f})")

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
            problem = self.make_problem(F_ext_global, next_load_factor, stage)
            res = self.corrector.run(problem, U_current.copy(), check_initial=False)

            if res.converged:
                _log.info(f"[PASO {step}] λ = {next_load_factor:.4f} | "
                          f"Δλ = {delta_lambda:.4f} | {res.n_solves} iter")
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
                    _log.debug(f"  -> Acelerando el próximo incremento a {delta_lambda:.4f}")

                if step_callback:
                    step_callback(step, U_current, load_factor)
                continue

            if self.adaptive:
                delta_lambda /= 2.0
                n_bisections += 1
                self.n_bisections = n_bisections
                _log.warning(f"[PASO {step}] sin convergencia con λ = {next_load_factor:.4f} "
                             f"({res.failure_reason()}); bisección: Δλ = {delta_lambda:.4e}")
                if delta_lambda < self.min_delta_lambda:
                    # Clasificar el modo y lanzar excepción tipada (ADR 0011).
                    raise self._fail(res.divergence_error(
                        last_load_factor=next_load_factor,
                        n_bisections=n_bisections,
                        extra_message=(
                            f"Δλ={delta_lambda:.2e} < min={self.min_delta_lambda:.2e}; "
                            f"line_search={'on' if self.line_search else 'off'}; "
                            f"último α={res.last_alpha:.3f}."
                        ),
                    ), load_factor, U_current)
            else:
                raise self._fail(res.divergence_error(
                    last_load_factor=next_load_factor,
                    n_bisections=n_bisections,
                    extra_message=(
                        f"Paso {step} no convergió en {self.max_iter} iter; "
                        f"adaptive=False; line_search={'on' if self.line_search else 'off'}."
                    ),
                ), load_factor, U_current)

        self.lambda_final = load_factor
        self.reached_max_lambda = True
        self.termination = Termination.MAX_LAMBDA
        self.U_final = np.array(U_current, dtype=float, copy=True)
        bisecciones = f", {n_bisections} bisección(es)" if n_bisections else ""
        _log.info(f"{type(self).__name__}: carga completa (termination = 'max_lambda') "
                  f"en {self.steps_done} pasos{bisecciones}.")
        return U_current
