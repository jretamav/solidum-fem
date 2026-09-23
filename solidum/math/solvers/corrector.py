"""Corrector de Newton compartido por los solvers no lineales (ADR 0015).

Los cinco solvers iterativos del proyecto —Newton-Raphson incremental,
arc-length cilíndrico y por disipación, Newton-Newmark y Newton-HHT—
ejecutan el mismo bucle de corrección:

1. ensamblar **una vez** en el iterado corriente;
2. formar el residuo y evaluar la convergencia con el par
   ``(‖R(x_k)‖, ‖δ_{k−1}‖)`` (criterio dual, ADR 0007), calibrando el
   criterio en el primer ensamblaje;
3. si converge, comitear el estado trial de ese ensamblaje;
4. si no, calcular la corrección resolviendo con el backend algebraico
   (fallback Cholesky → LU, Newton modificado con factorización congelada,
   ADR 0003), aplicar el line search opcional (descenso no monótono,
   ADR 0011) y avanzar;
5. al agotar el presupuesto, devolver el historial de residuos para que el
   solver clasifique la divergencia (ADR 0011) o biseque su paso.

Lo que cambia de un solver a otro —qué es el iterado (``U``; ``(U, λ)``;
``(u, u̇, ü)``), cómo se forma el residuo, cómo se reduce y resuelve el
sistema tangente, cómo se aplica la corrección— lo aporta un objeto
:class:`NewtonProblem` que cada solver construye por paso. El corrector
no conoce la física ni la cinemática temporal: sólo el protocolo.

Por qué un objeto por paso y no una función: el problema lleva datos del
paso (carga del paso, predictores de Newmark, longitud de arco, término
``α·X_n`` de HHT) y estado del paso (``dU`` acumulado, ``λ`` corriente)
que sólo tienen sentido dentro de él; el corrector, en cambio, conserva
entre pasos lo que sí persiste: el backend algebraico y su degradación a
LU.

Equivalencia: la migración de cada solver a este corrector no cambia
ninguna operación ni su orden; los resultados son los mismos bit a bit
salvo donde se indique en el ADR.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Protocol

import numpy as np

from solidum.constants import LINE_SEARCH_MAX_BACKTRACKS, LINE_SEARCH_RHO
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.linalg import StiffnessProperties, select_solver
from solidum.math.solvers._shared import CholeskyNotPositiveDefiniteError, _log
from solidum.math.solvers.diagnostics import SolverDivergedError, classify_divergence

SolveFn = Callable[[Any, np.ndarray], np.ndarray]


class CorrectionAborted(Exception):
    """El problema no puede producir la corrección de esta iteración
    (raíces imaginarias del arco, ``α ≈ 0`` en la restricción de
    disipación, …). El corrector abandona el paso como *no convergido*
    sin marcarlo como tangente singular."""


class NewtonProblem(Protocol):
    """Lo que un solver aporta al corrector para un paso.

    ``x`` es el iterado, opaco para el corrector (un vector, una tupla,
    …); ``state`` es lo que devuelve :meth:`assemble` (típicamente
    ``(K, F_int)``), también opaco. Todas las evaluaciones del paso se
    hacen a través de estos métodos, así que un solver nuevo escribe sólo
    su física y hereda el bucle, el backend, el line search y la
    telemetría.
    """

    def assemble(self, x) -> Any:
        """Ensambla en ``x`` (único ensamblaje por iteración) y devuelve el
        estado que consumen el resto de métodos. Deja el trial de los
        elementos evaluado en ``x``."""

    def residual(self, x, state) -> np.ndarray:
        """Residuo completo (todos los DOF) en ``x``."""

    def residual_norm(self, R: np.ndarray) -> float:
        """Norma del residuo que entra al criterio: en DOF libres
        (``‖R[libres]‖``) o reducida (``‖TᵀR‖``)."""

    def calibration_scales(self, x, state) -> tuple[float, float]:
        """``(escala de fuerza, escala de desplazamiento)`` para calibrar el
        criterio en el primer ensamblaje (ADR 0007)."""

    def reference_force(self, x, state) -> float:
        """Escala de fuerza corriente del criterio (``max(‖F_ext‖, ‖F_int‖)``)."""

    def x_norm(self, x) -> float:
        """Norma del iterado para el criterio en desplazamiento."""

    def correction(self, x, state, R: np.ndarray, solve: SolveFn) -> Any:
        """Corrección ``dx`` a partir del sistema tangente. Resuelve con
        ``solve(A, b)`` (backend del corrector, con fallback y Newton
        modificado). Lanza ``RuntimeError`` si la tangente es singular
        (lo propaga el backend) o :class:`CorrectionAborted` si la
        corrección no existe."""

    def apply(self, x, dx, alpha: float) -> tuple[Any, float]:
        """``(x + α·dx, ‖incremento‖)``: el iterado avanzado y la norma del
        incremento que entra al criterio en desplazamiento."""

    def on_converged(self, x, state) -> None:
        """Se llama con el iterado convergido y el estado de su ensamblaje
        (típicamente: comitear)."""


@dataclass
class CorrectorResult:
    """Resultado de :meth:`NewtonCorrector.run` sobre un paso."""

    converged: bool
    x: Any
    state: Any
    n_solves: int
    residual_history: list[float] = field(default_factory=list)
    delta_history: list[float] = field(default_factory=list)
    singular_tangent: bool = False
    aborted: bool = False
    last_alpha: float = 1.0
    last_residual: float = float("inf")
    last_delta: float = 0.0

    def divergence_error(self, *, last_load_factor: float, n_bisections: int = 0,
                         extra_message: str = "") -> SolverDivergedError:
        """Excepción tipada del modo de divergencia observado (ADR 0011),
        lista para lanzar."""
        err_cls = classify_divergence(
            self.residual_history, self.delta_history,
            singular_tangent_detected=self.singular_tangent,
        )
        return err_cls(
            last_residual=self.last_residual,
            last_delta=self.last_delta,
            last_load_factor=last_load_factor,
            n_bisections=n_bisections,
            extra_message=extra_message,
        )


class NewtonCorrector:
    """Bucle de corrección de Newton con backend algebraico persistente.

    Parameters
    ----------
    convergence
        Criterio dual (ADR 0007). Se calibra en el primer ensamblaje.
    max_iter
        Presupuesto de resoluciones por paso (``max_iter + 1`` evaluaciones
        del residuo).
    is_symmetric, is_positive_definite
        Propiedades declaradas de la tangente para el despachador
        (ADR 0003). Si Cholesky reporta no-positividad el corrector
        degrada a LU para el resto del análisis.
    linear_algebra
        Override del backend (``"auto"``, ``"cholesky"``, ``"lu"``).
    freeze_tangent_after_iter
        Newton modificado (ADR 0003 fase 2): factoriza fresca las primeras
        ``N`` iteraciones del paso y reusa la factorización después.
        ``None`` ⇒ Newton estándar.
    line_search
        Descenso no monótono (ADR 0011): acepta ``α = 1`` si baja el
        residuo y hace backtracking sólo cuando el paso completo lo sube.
    verbose
        Registra una línea por iteración evaluada (INFO). Los solvers
        transitorios lo desactivan y registran una línea por paso.
    """

    def __init__(self, convergence: ConvergenceCriterion, *, max_iter: int,
                 is_symmetric: bool, is_positive_definite: bool = True,
                 linear_algebra: str = "auto",
                 freeze_tangent_after_iter: int | None = None,
                 line_search: bool = False, verbose: bool = True):
        self.convergence = convergence
        self.max_iter = int(max_iter)
        self.is_symmetric = bool(is_symmetric)
        self.linear_algebra = str(linear_algebra)
        self.freeze_tangent_after_iter = freeze_tangent_after_iter
        self.line_search = bool(line_search)
        self.verbose = bool(verbose)
        self._is_pd = bool(is_positive_definite)
        self._linalg = None
        self._linalg_size: int | None = None
        self._frozen_factor = None

    # ------------------------------------------------------------------
    # Backend algebraico
    # ------------------------------------------------------------------

    @property
    def is_positive_definite(self) -> bool:
        """``False`` una vez degradado a LU."""
        return self._is_pd

    @property
    def frozen_factor(self):
        """Factorización congelada vigente (Newton modificado) o ``None``."""
        return self._frozen_factor

    def _backend(self, n: int):
        if self._linalg is None or self._linalg_size != n:
            props = StiffnessProperties(
                is_symmetric=self.is_symmetric,
                is_positive_definite=self._is_pd,
                size=n,
            )
            self._linalg = select_solver(props, override=self.linear_algebra)
            self._linalg_size = n
        return self._linalg

    def _degrade_to_lu(self, n: int):
        _log.warning("Cholesky reportó no-positividad. Degradando a LU para el resto del análisis.")
        self._is_pd = False
        self._linalg = None
        self._frozen_factor = None
        return self._backend(n)

    def solve(self, A, b: np.ndarray, iteration: int = 0) -> np.ndarray:
        """``A·x = b`` con la política del corrector: fallback Cholesky → LU
        y, con ``freeze_tangent_after_iter``, factorización congelada a
        partir de la iteración indicada. Propaga ``RuntimeError`` del
        backend (tangente singular)."""
        n = A.shape[0]
        linalg = self._backend(n)
        if self.freeze_tangent_after_iter is None:
            try:
                return linalg.solve(A, b)
            except CholeskyNotPositiveDefiniteError:
                return self._degrade_to_lu(n).solve(A, b)
        if iteration < self.freeze_tangent_after_iter or self._frozen_factor is None:
            try:
                self._frozen_factor = linalg.factorize(A)
            except CholeskyNotPositiveDefiniteError:
                self._frozen_factor = self._degrade_to_lu(n).factorize(A)
        return self._frozen_factor.solve(b)

    def reset_step(self) -> None:
        """Descarta la factorización congelada: la tangente cambia con el
        iterado convergido, así que no vale para el paso siguiente."""
        self._frozen_factor = None

    # ------------------------------------------------------------------
    # Bucle de corrección
    # ------------------------------------------------------------------

    def run(self, problem: NewtonProblem, x0, *, check_initial: bool = True,
            initial_delta_norm: float = 0.0) -> CorrectorResult:
        """Corrige ``x0`` hasta converger o agotar ``max_iter`` resoluciones.

        Parameters
        ----------
        check_initial
            Evaluar el criterio en la iteración 0 (antes de resolver). El
            Newton incremental lo desactiva: su iterado inicial es el
            convergido del paso anterior y aún no incorpora el incremento
            de Dirichlet del paso, así que todo paso hace al menos una
            resolución. Los solvers con predictor propio (arc-length,
            Newmark) sí evalúan el iterado predicho.
        initial_delta_norm
            Norma del incremento que produjo ``x0`` (el predictor), para el
            criterio en desplazamiento de la iteración 0.
        """
        x = x0
        state = problem.assemble(x)
        result = CorrectorResult(converged=False, x=x, state=state, n_solves=0)
        delta_norm = float(initial_delta_norm)
        R_norm = float("inf")

        for iteration in range(self.max_iter + 1):
            R = problem.residual(x, state)
            R_norm = problem.residual_norm(R)

            if not self.convergence.is_calibrated:
                force_scale, disp_scale = problem.calibration_scales(x, state)
                self.convergence.calibrate(force_scale, disp_scale)

            conv = None
            if iteration > 0 or check_initial:
                conv = self.convergence.evaluate(
                    residual_norm=R_norm,
                    ref_force=problem.reference_force(x, state),
                    delta_u_norm=delta_norm,
                    u_norm=problem.x_norm(x),
                )
                result.residual_history.append(R_norm)
                result.delta_history.append(delta_norm)
                result.last_residual = R_norm
                result.last_delta = delta_norm
                if self.verbose:
                    context = getattr(problem, "log_context", None)
                    ctx = context(x) if context is not None else ""
                    alpha_tag = f" | α={result.last_alpha:.3f}" if result.last_alpha < 1.0 else ""
                    _log.info(
                        f"  Iteración {result.n_solves:2d}{ctx} | "
                        f"R/tol_F: {conv.ratio_force:.4e} | "
                        f"dU/tol_d: {conv.ratio_disp:.4e}{alpha_tag}"
                    )

            if conv is not None and conv.converged:
                problem.on_converged(x, state)
                result.converged = True
                break

            if iteration == self.max_iter:
                break  # presupuesto de resoluciones agotado

            try:
                dx = problem.correction(
                    x, state, R, lambda A, b: self.solve(A, b, iteration=iteration),
                )
            except CorrectionAborted as exc:
                _log.error(str(exc))
                result.aborted = True
                break
            except RuntimeError:
                _log.error("Matriz singular detectada.")
                result.singular_tangent = True
                break
            result.n_solves += 1

            alpha, x, state, delta_norm = self.line_search_step(problem, x, dx, R_norm)
            result.last_alpha = alpha

        self.reset_step()
        result.x = x
        result.state = state
        return result

    def line_search_step(self, problem: NewtonProblem, x, dx, R_norm_current: float):
        """Avanza ``x`` con ``dx`` y ensambla en el punto aceptado.

        Devuelve ``(α, x_nuevo, estado, ‖incremento‖)``. Sin line search
        aplica ``α = 1`` y ensambla una vez. Con line search (descenso no
        monótono, variante GLL del ADR 0011) prueba ``α = 1, ρ, ρ², …``
        ensamblando en cada punto hasta que ``‖R‖`` no supere
        ``R_norm_current``; si agota los retrocesos devuelve el último
        punto probado y delega al control externo (bisección del paso).
        El estado devuelto es siempre el del ensamblaje en el punto
        aceptado: no hay reensamblaje redundante.
        """
        if not self.line_search:
            x_new, delta_norm = problem.apply(x, dx, 1.0)
            return 1.0, x_new, problem.assemble(x_new), delta_norm

        alpha = 1.0
        x_trial = state = None
        delta_norm = 0.0
        for _ in range(LINE_SEARCH_MAX_BACKTRACKS + 1):
            x_trial, delta_norm = problem.apply(x, dx, alpha)
            state = problem.assemble(x_trial)
            R_trial = problem.residual(x_trial, state)
            if problem.residual_norm(R_trial) <= R_norm_current:
                return alpha, x_trial, state, delta_norm
            alpha *= LINE_SEARCH_RHO
        return alpha, x_trial, state, delta_norm


def default_calibration_scales(F_ext: np.ndarray, F_int: np.ndarray, K) -> tuple[float, float]:
    """Escalas de calibración comunes a todos los solvers (ADR 0007):
    ``force = max(‖F_ext‖, ‖F_int‖, 1)`` y ``disp = force / max|diag K|``."""
    from solidum.math.convergence import stiffness_diag_scale
    force_scale = max(float(np.linalg.norm(F_ext)), float(np.linalg.norm(F_int)), 1.0)
    return force_scale, force_scale / stiffness_diag_scale(K)
