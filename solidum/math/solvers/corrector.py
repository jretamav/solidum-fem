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

import functools
import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Protocol

import numpy as np
import scipy.sparse as sp

from solidum.constants import (
    EQUILIBRIUM_RTOL,
    LINE_SEARCH_MAX_BACKTRACKS,
    LINE_SEARCH_RHO,
    MAX_STEP_REVISIONS,
)
from solidum.math.convergence import ConvergenceCriterion, ConvergenceState
from solidum.math.linalg import IterativeNotConvergedError, StiffnessProperties, select_solver
from solidum.math.solvers._shared import CholeskyNotPositiveDefiniteError, _log
from solidum.math.solvers.diagnostics import (
    IllPosedSystemError,
    SolverDivergedError,
    StepConstraintError,
    StepRevisionError,
    classify_divergence,
)

SolveFn = Callable[[Any, np.ndarray], np.ndarray]


@dataclass
class _Retained:
    """Factorización retenida por el corrector y una copia de la matriz que
    factorizó (ADR 0040): la copia sirve para reconocer una matriz idéntica y
    para factorizarla otra vez si el backend perdió el factor."""

    factor: Any
    matrix: Any


def _copy_matrix(A):
    """Copia de ``A`` (datos, índices y punteros si es dispersa): quien la
    pasó puede modificarla después en sitio."""
    return A.copy() if sp.issparse(A) else np.array(A, copy=True)


def _same_matrix(A, B) -> bool:
    """``True`` si ``A`` y ``B`` son la misma matriz bit a bit: misma forma y
    formato y, si son CSR/CSC, mismos punteros, índices y valores (un cero
    explícito no equivale a una entrada ausente; un NaN nunca coincide)."""
    if B is None or A.shape != B.shape:
        return False
    if sp.issparse(A) and sp.issparse(B):
        if A.format != B.format or A.format not in ("csr", "csc") or A.nnz != B.nnz:
            return False
        return (np.array_equal(A.indptr, B.indptr) and np.array_equal(A.indices, B.indices)
                and np.array_equal(A.data, B.data))
    if isinstance(A, np.ndarray) and isinstance(B, np.ndarray):
        return np.array_equal(A, B)
    return False


def _relative_residual(A, b: np.ndarray, x: np.ndarray) -> float:
    """``‖b − A·x‖/‖b‖`` (``inf`` si no es finito; 0 con ``b`` nulo)."""
    b_norm = float(np.linalg.norm(b))
    if b_norm == 0.0:
        return 0.0
    rel = float(np.linalg.norm(b - A @ x)) / b_norm
    return rel if np.isfinite(rel) else np.inf


def _finite(x: np.ndarray) -> np.ndarray:
    """Vigilancia única de la solución de toda resolución del corrector
    (fresca, reutilizada o congelada): una no finita delata una tangente
    singular o casi singular (ADR 0011)."""
    if not np.all(np.isfinite(x)):
        raise RuntimeError(
            "La solución del sistema tangente no es finita (matriz singular o "
            "casi singular).")
    return x


class CorrectionAborted(Exception):
    """El problema no puede producir la corrección de esta iteración
    (raíces imaginarias del arco, denominador nulo en la restricción del
    control indirecto de desplazamiento, …). El corrector abandona el paso
    como *no convergido*
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
        """Norma del residuo que entra al criterio: la del equilibrio
        reducido, ``‖TᵀR‖`` con ``u = T·u_libre + g`` (ADR 0004). Con sólo
        condiciones de Dirichlet coincide con ``‖R[libres]‖``; con
        restricciones multipunto ``R[libres]`` no se anula en la solución,
        porque la fuerza de cada esclavo pasa a sus maestros."""

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

    # Opcional (ADR 0041):
    # def revise(self, x) -> bool:
    #     """Revisión del paso convergido, antes de ``on_converged``: ``True``
    #     si algún elemento cambió una decisión discreta con el campo de
    #     ``x`` (``assembler.revise_all_steps``). El corrector no confirma y
    #     repite la corrección desde el iterado inicial del intento."""


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
    # Diagnóstico del backend algebraico cuando el sistema tangente no se
    # pudo resolver por una causa distinta de la singularidad (solver
    # iterativo sin convergencia, ADR 0018). Se añade al mensaje de la
    # excepción tipada para que el usuario vea la causa real.
    linear_solver_note: str = ""
    # Mayor residuo relativo ‖b − A·x‖/‖b‖ de las resoluciones del paso con
    # factorización fresca (ADR 0019). Sólo diagnóstico: el algoritmo no lo
    # usa (resolver sistemas casi singulares es legítimo cerca de un punto
    # límite), pero si el paso fracasa y es grande, la causa probable es una
    # tangente singular —p. ej. un mecanismo— y así se informa.
    max_linear_residual: float = 0.0
    aborted: bool = False
    # Motivo con que el problema abortó la corrección (``CorrectionAborted``):
    # p. ej. las raíces imaginarias de la restricción del arco (ADR 0025).
    abort_reason: str = ""
    last_alpha: float = 1.0
    last_residual: float = float("inf")
    last_delta: float = 0.0
    last_conv: ConvergenceState | None = None
    # Revisión del paso convergido (ADR 0041): veces que el intento se repitió
    # porque un elemento cambió una decisión discreta, resoluciones de las
    # rondas descartadas (``n_solves`` cuenta sólo la última) y si se agotó
    # el tope ``MAX_STEP_REVISIONS``.
    n_revisions: int = 0
    revision_solves: int = 0
    revisions_exhausted: bool = False

    def revision_note(self) -> str:
        """Fragmento para la línea INFO del paso del solver (ADR 0025): las
        revisiones del paso convergido que lo repitieron (ADR 0041)."""
        return revision_note(self.n_revisions)

    def failure_reason(self) -> str:
        """Motivo breve de un paso que no convergió, para el registro del
        solver que lo reintenta (ADR 0025)."""
        if self.revisions_exhausted:
            return (f"el estado discreto de los elementos siguió cambiando tras "
                    f"{MAX_STEP_REVISIONS} revisiones del paso convergido (ADR 0041): "
                    f"las decisiones de revise_step no se estabilizaron en el paso")
        if self.aborted:
            return self.abort_reason or "la restricción del paso no tiene solución"
        if self.linear_solver_note:
            return self.linear_solver_note
        if self.singular_tangent:
            return "tangente singular"
        return f"sin convergencia en {self.n_solves} resoluciones"

    def divergence_error(self, *, last_load_factor: float, n_bisections: int = 0,
                         extra_message: str = "") -> SolverDivergedError:
        """Excepción tipada del modo de divergencia observado (ADR 0011),
        lista para lanzar. Si el problema abortó porque la restricción del
        paso no tiene solución, :class:`StepConstraintError` con el motivo
        (ADR 0025). Si se agotaron las revisiones del paso convergido (ADR
        0041), :class:`StepRevisionError`: el Newton convergía, lo que no
        termina es el estado discreto."""
        if self.revisions_exhausted:
            return StepRevisionError(
                last_residual=self.last_residual,
                last_delta=self.last_delta,
                last_load_factor=last_load_factor,
                n_bisections=n_bisections,
                extra_message=(f"{extra_message}; " if extra_message else "") + self.failure_reason(),
            )
        if self.aborted:
            reason = self.failure_reason()
            return StepConstraintError(
                last_residual=self.last_residual,
                last_delta=self.last_delta,
                last_load_factor=last_load_factor,
                n_bisections=n_bisections,
                extra_message=(f"{extra_message}; " if extra_message else "") + reason,
            )
        bad_linear = self.max_linear_residual > EQUILIBRIUM_RTOL
        err_cls = classify_divergence(
            self.residual_history, self.delta_history,
            singular_tangent_detected=self.singular_tangent or bad_linear,
        )
        note = self.linear_solver_note
        if bad_linear and not note:
            note = (f"el sistema tangente no se pudo resolver con precisión "
                    f"(residuo lineal relativo {self.max_linear_residual:.1e}): "
                    f"la tangente es singular o casi singular")
        if note:
            extra_message = (f"{extra_message}; " if extra_message else "") + note
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
        Override del backend (``"auto"``, ``"cholesky"``, ``"pardiso"``,
        ``"lu"``, ``"iterative[:precondicionador]"``).
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
    near_nullspace
        Proveedor perezoso del casi-núcleo del sistema reducido
        (``Assembler.near_nullspace``) para el precondicionador AMG del
        backend iterativo (ADR 0018). Lo ignoran los demás backends.
    """

    def __init__(self, convergence: ConvergenceCriterion, *, max_iter: int,
                 is_symmetric: bool, is_positive_definite: bool = True,
                 linear_algebra: str = "auto",
                 freeze_tangent_after_iter: int | None = None,
                 line_search: bool = False, verbose: bool = True,
                 near_nullspace=None):
        self.convergence = convergence
        self.max_iter = int(max_iter)
        self.is_symmetric = bool(is_symmetric)
        self.linear_algebra = str(linear_algebra)
        self.freeze_tangent_after_iter = freeze_tangent_after_iter
        self.line_search = bool(line_search)
        self.verbose = bool(verbose)
        self._is_pd = bool(is_positive_definite)
        self.near_nullspace = near_nullspace
        self._linalg = None
        self._linalg_size: int | None = None
        # Factorizaciones retenidas (ADR 0040): la última fresca, que sobrevive
        # al paso y se reutiliza si la matriz es idéntica, y la congelada del
        # Newton modificado, por paso.
        self._last: _Retained | None = None
        self._frozen: _Retained | None = None
        # La cadena de análisis simbólicos (ADR 0040) se rompe hasta el
        # próximo ``release`` en cuanto una factorización perturba pivotes.
        self._chain_ok = True
        self._max_linear_residual = 0.0
        self._factor_inspector: Callable[[Any], None] | None = None

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
        return None if self._frozen is None else self._frozen.factor

    def _backend(self, n: int):
        if self._linalg is None or self._linalg_size != n:
            props = StiffnessProperties(
                is_symmetric=self.is_symmetric,
                is_positive_definite=self._is_pd,
                size=n,
                near_nullspace=self.near_nullspace,
            )
            self._linalg = select_solver(props, override=self.linear_algebra)
            self._linalg_size = n
        return self._linalg

    def _degrade_to_lu(self, n: int):
        _log.warning("Cholesky reportó no-positividad. Degradando a LU para el resto del análisis.")
        self._is_pd = False
        self._linalg = None
        self._last = None
        self._frozen = None
        return self._backend(n)

    def inspect_next_factorization(self, inspector: Callable[[Any, Any], None]) -> None:
        """Entrega a ``inspector`` la próxima factorización fresca, antes de
        usarla, y lo olvida.

        Lo usa la red de seguridad del ADR 0019
        (:func:`~solidum.math.solvers.model_checks.check_initial_tangent`)
        para buscar mecanismos internos en la tangente inicial sin factorizar
        dos veces. ``inspector(factor, A)`` recibe el objeto factorizado (con
        ``n_small_pivots``) y la matriz factorizada, o ``factor = None`` si la
        factorización directa falló por un pivote exactamente nulo; en ese
        caso el ``RuntimeError`` del backend se propaga después, como
        siempre. Si ``inspector`` lanza :class:`IllPosedSystemError`, el
        corrector no la toma por una tangente singular: se propaga hasta el
        usuario."""
        self._factor_inspector = inspector

    def _factorize(self, linalg, A, n: int, previous=None):
        """``linalg.factorize(A)`` con el fallback Cholesky → LU.
        ``previous``: factor anterior de la cadena cuyo análisis simbólico
        puede reutilizar el backend (ADR 0040); sólo llega a los backends que
        lo declaran (``reuses_analysis``), así que los demás no tienen que
        aceptar el argumento."""
        try:
            if previous is not None:
                return linalg.factorize(A, reuse_analysis_of=previous)
            return linalg.factorize(A)
        except CholeskyNotPositiveDefiniteError:
            return self._degrade_to_lu(n).factorize(A)

    def _inspected_factorize(self, linalg, A, n: int, previous=None):
        """Factoriza y entrega el factor al inspector pendiente (ver
        :meth:`inspect_next_factorization`)."""
        inspector, self._factor_inspector = self._factor_inspector, None
        try:
            factor = self._factorize(linalg, A, n, previous)
        except RuntimeError:
            inspector(None, A)       # pivote exactamente nulo
            raise
        inspector(factor, A)
        return factor

    def _current(self, entry: _Retained, linalg, n: int):
        """El factor de ``entry`` listo para resolver: si otra factorización
        lo reemplazó en el backend (``is_current`` falso), se factoriza otra
        vez **su** matriz, lo que conserva la semántica del Newton
        modificado (ADR 0040)."""
        if not getattr(entry.factor, "is_current", True):
            entry.factor = self._factorize(linalg, entry.matrix, n)
        return entry.factor

    def _fresh_factor(self, linalg, A, n: int):
        """Factorización de ``A``: la última, si ``A`` es idéntica a la matriz
        que factorizó (forma, formato, patrón y valores, bit a bit); si no,
        una nueva, que pasa a ser la última (ADR 0040).

        Con IMPL-EX o con material elástico la tangente no depende del
        iterado y se reutiliza dentro del paso (y entre pasos elásticos); con
        daño implícito, plasticidad, contacto o no linealidad geométrica la
        matriz cambia y se factoriza como siempre. La comparación no pide
        nada a materiales ni elementos."""
        last = self._last
        if last is not None and _same_matrix(A, last.matrix):
            factor = self._current(last, linalg, n)
            if self._factor_inspector is not None:
                inspector, self._factor_inspector = self._factor_inspector, None
                inspector(factor, last.matrix)
            # Un inspector puede factorizar otra matriz en el mismo backend
            # (Pardiso guarda una sola): se comprueba otra vez.
            return self._current(last, linalg, n)
        # Matriz nueva. La factorización anterior se suelta antes de calcular
        # la nueva, para no tener dos vivas a la vez (memoria pico); sólo se
        # conserva si el backend puede reutilizar su análisis simbólico, y
        # entonces la memoria del factor ya es la del handle (Pardiso).
        previous = None
        if last is not None and self._chain_ok and getattr(linalg, "reuses_analysis", False):
            previous = last.factor
        self._last = last = None
        # Se factoriza la copia, que sólo posee el corrector: el backend puede
        # conservar la matriz (Pardiso la usa en el refinamiento de cada
        # resolución), y quien pasó ``A`` puede modificarla después en sitio.
        matrix = _copy_matrix(A)
        if self._factor_inspector is not None:
            factor = self._inspected_factorize(linalg, matrix, n, previous)
        else:
            factor = self._factorize(linalg, matrix, n, previous)
        self._last = _Retained(factor, matrix)
        self._watch_chain(factor)
        return self._current(self._last, linalg, n)

    def _watch_chain(self, factor) -> None:
        """Rompe la cadena de análisis simbólicos (ADR 0040) para el resto
        del análisis si la factorización perturbó pivotes: la matriz es
        singular o casi (un mecanismo), y entonces la solución depende de la
        escala de la factorización aunque su residuo sea pequeño. Medido: una
        armadura corrotacional con mecanismo «convergía» a un colapso con el
        análisis reutilizado en vez de dar tangente singular. Sin cadena, cada
        factorización hace su propio análisis, como antes del ADR 0040. Un
        análisis bien planteado no perturba pivotes y no pierde nada."""
        if getattr(factor, "n_perturbed_pivots", 0):
            self._chain_ok = False

    def solve(self, A, b: np.ndarray, iteration: int = 0) -> np.ndarray:
        """``A·x = b`` con la política del corrector: fallback Cholesky → LU,
        factorización reutilizada si ``A`` es idéntica a la última
        factorizada (ADR 0040) y, con ``freeze_tangent_after_iter``,
        factorización congelada a partir de la iteración indicada. Propaga
        ``RuntimeError`` del backend (tangente singular) y lo lanza si la
        solución no es finita."""
        n = A.shape[0]
        linalg = self._backend(n)
        if (self.freeze_tangent_after_iter is not None and self._frozen is not None
                and iteration >= self.freeze_tangent_after_iter
                and self._frozen.matrix.shape == A.shape):
            # Tangente congelada: se resuelve con la matriz de una iteración
            # anterior a propósito, así que el residuo frente a A no significa
            # nada. Si las restricciones cambiaron el número de incógnitas
            # entre iteraciones (contacto), la congelada no sirve y se
            # factoriza la matriz nueva.
            return _finite(self._current(self._frozen, linalg, n).solve(b))
        factor = self._fresh_factor(linalg, A, n)
        if self.freeze_tangent_after_iter is not None:
            self._frozen = self._last
        x = self._solve_with_reference_fallback(factor, linalg, A, b, n)
        self._record_linear_residual(A, b, x)
        return x

    def _solve_with_reference_fallback(self, factor, linalg, A, b: np.ndarray, n: int):
        """``factor.solve(b)``; si el factor reutilizó el análisis simbólico de
        otra matriz (``shares_analysis``) y la resolución no es exacta —residuo
        relativo por encima de ``EQUILIBRIUM_RTOL``, el umbral con el que la
        red de seguridad declara casi singular una resolución (ADR 0019), o
        solución no finita—, se repite con la factorización de referencia, de
        análisis propio (ADR 0040).

        En un sistema bien condicionado el análisis reutilizado sólo cambia el
        redondeo; en uno casi singular (un mecanismo) el redondeo se amplifica
        con el condicionamiento y la solución depende de la escala de la
        factorización: medido, una armadura corrotacional con mecanismo
        «convergía» a un colapso en lugar de dar tangente singular. Con la
        de referencia, el diagnóstico es el de siempre."""
        try:
            x = factor.solve(b)
        except RuntimeError:
            if not getattr(factor, "shares_analysis", False):
                raise
            x = None
        if getattr(factor, "shares_analysis", False):
            if x is None or not (_relative_residual(A, b, x) <= EQUILIBRIUM_RTOL):
                entry = self._last
                entry.factor = None
                entry.factor = factor = self._factorize(linalg, entry.matrix, n)
                self._watch_chain(factor)
                x = factor.solve(b)
        return _finite(x)

    def _record_linear_residual(self, A, b: np.ndarray, x: np.ndarray) -> None:
        """Telemetría del ADR 0019: residuo relativo de la resolución (un
        producto matriz-vector). No interviene en el algoritmo."""
        if float(np.linalg.norm(b)) > 0.0:
            self._max_linear_residual = max(self._max_linear_residual,
                                            _relative_residual(A, b, x))

    def reset_step(self) -> None:
        """Descarta la factorización congelada: la tangente cambia con el
        iterado convergido, así que no vale para el paso siguiente. La última
        factorización se conserva: sólo se reutiliza con una matriz idéntica."""
        self._frozen = None

    def release(self) -> None:
        """Suelta las factorizaciones retenidas. Lo llama cada solver al
        terminar su ``solve()``: con LU, Cholesky o el iterativo el factor
        ocupa memoria propia, y la primera factorización del análisis
        siguiente hace su propio análisis simbólico (ADR 0040)."""
        self._last = None
        self._frozen = None
        self._chain_ok = True

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

        Revisión del paso convergido (ADR 0041): si el problema tiene
        ``revise(x)`` y devuelve ``True`` al converger, el intento no se
        confirma y la corrección se repite desde ``x0`` con la decisión
        discreta nueva, hasta que una revisión no cambia nada o se agota
        ``MAX_STEP_REVISIONS`` (entonces el intento fracasa).
        """
        x = x0
        state = problem.assemble(x)
        result = CorrectorResult(converged=False, x=x, state=state, n_solves=0)
        self._max_linear_residual = 0.0
        delta_norm = float(initial_delta_norm)
        R_norm = float("inf")
        revise = getattr(problem, "revise", None)

        iteration = 0
        while True:
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
                result.last_conv = conv
                # Detalle por iteración en DEBUG (ADR 0025): en INFO cada
                # solver escribe una línea por paso.
                if self.verbose and _log.isEnabledFor(logging.DEBUG):
                    context = getattr(problem, "log_context", None)
                    ctx = context(x) if context is not None else ""
                    alpha_tag = f" | α={result.last_alpha:.3f}" if result.last_alpha < 1.0 else ""
                    _log.debug(
                        f"  Iteración {result.n_solves:2d}{ctx} | "
                        f"R/tol_F: {conv.ratio_force:.4e} | "
                        f"dU/tol_d: {conv.ratio_disp:.4e}{alpha_tag}"
                    )

            if conv is not None and conv.converged:
                if revise is not None and revise(x):
                    # Un elemento cambió una decisión discreta con el campo
                    # convergido: el intento se repite desde ``x0`` sin
                    # confirmar (el trial se recalcula desde el estado
                    # confirmado en el primer ensamblaje).
                    result.n_revisions += 1
                    result.revision_solves += result.n_solves
                    if result.n_revisions > MAX_STEP_REVISIONS:
                        result.revisions_exhausted = True
                        _log.debug(result.failure_reason())
                        break
                    _log.debug(f"  Revisión {result.n_revisions} del paso convergido: "
                               f"cambió el estado discreto; se repite el intento.")
                    self.reset_step()
                    self._max_linear_residual = 0.0
                    x = x0
                    state = problem.assemble(x)
                    delta_norm = float(initial_delta_norm)
                    result.n_solves = 0
                    result.residual_history = []
                    result.delta_history = []
                    result.last_alpha = 1.0
                    result.last_residual = float("inf")
                    result.last_delta = 0.0
                    result.last_conv = None
                    iteration = 0
                    continue
                problem.on_converged(x, state)
                result.converged = True
                if result.n_revisions:
                    _log.debug(f"  Paso repetido {result.n_revisions} vez/veces: el campo "
                               f"convergido cambió el estado discreto de algún elemento "
                               f"({result.revision_solves} resoluciones descartadas; ADR 0041).")
                break

            if iteration == self.max_iter:
                break  # presupuesto de resoluciones agotado

            try:
                dx = problem.correction(
                    x, state, R, lambda A, b: self.solve(A, b, iteration=iteration),
                )
            # Los fracasos del paso se registran en DEBUG: el solver que llama
            # decide si el paso se recupera (bisección, WARNING con el motivo
            # de ``failure_reason``) o termina el análisis (excepción tipada).
            # ERROR queda para lo que detiene el análisis (ADR 0025).
            except CorrectionAborted as exc:
                _log.debug(str(exc))
                result.aborted = True
                result.abort_reason = str(exc)
                break
            except IterativeNotConvergedError as exc:
                # El sistema tangente no se resolvió, pero no por
                # singularidad: se registra la causa real (ADR 0018 §5).
                _log.debug(str(exc))
                result.singular_tangent = True
                result.linear_solver_note = str(exc)
                break
            except IllPosedSystemError:
                # Veredicto de la red de seguridad sobre la tangente inicial
                # (ADR 0019): un error de planteamiento del modelo, no una
                # tangente singular que bisecar.
                raise
            except RuntimeError:
                _log.debug("Matriz singular detectada.")
                result.singular_tangent = True
                break
            result.n_solves += 1

            alpha, x, state, delta_norm = self.line_search_step(problem, x, dx, R_norm)
            result.last_alpha = alpha
            iteration += 1

        self.reset_step()
        result.x = x
        result.state = state
        result.max_linear_residual = self._max_linear_residual
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


def revision_note(n_revisions: int) -> str:
    """Fragmento de la línea INFO de un paso con ``n_revisions`` revisiones
    del paso convergido (ADR 0041, ADR 0025); vacío si no hubo."""
    return f" | {n_revisions} revisión(es)" if n_revisions else ""


def releases_corrector(solve):
    """Decorador del ``solve()`` de todo solver con uno o más
    :class:`NewtonCorrector` como atributos (ADR 0040): suelta sus
    factorizaciones retenidas al empezar —la primera factorización del
    análisis hace su propio análisis simbólico, sin depender de un análisis
    anterior— y al terminar —con LU, Cholesky o el iterativo el factor ocupa
    memoria propia—. Una prueba exige el decorador en todo solver registrado
    que tenga corrector.

    Al terminar también cierra, sin confirmarlo, el intento de paso que haya
    quedado abierto (ADR 0034, ADR 0041): un trazado que acaba tras un intento
    fallido (``max_steps`` agotado después de una bisección, predictor
    singular en el último paso) no deja un intento abierto que el commit de
    posproceso de ``solidum.run`` tomaría por un paso sin confirmar."""
    @functools.wraps(solve)
    def wrapper(self, *args, **kwargs):
        _release_correctors(self)
        try:
            return solve(self, *args, **kwargs)
        finally:
            _release_correctors(self)
            assembler = getattr(self, "assembler", None)
            close = getattr(assembler, "close_step_attempt", None)
            if close is not None:
                close()
    wrapper.__releases_corrector__ = True
    return wrapper


def _release_correctors(solver) -> None:
    """Suelta todo :class:`NewtonCorrector` que el solver guarde como
    atributo (un esquema escalonado puede tener uno por campo)."""
    for value in vars(solver).values():
        if isinstance(value, NewtonCorrector):
            value.release()


def default_calibration_scales(F_ext: np.ndarray, F_int: np.ndarray, K) -> tuple[float, float]:
    """Escalas de calibración comunes a todos los solvers (ADR 0007):
    ``force = max(‖F_ext‖, ‖F_int‖, 1)`` y ``disp = force / max|diag K|``."""
    from solidum.math.convergence import stiffness_diag_scale
    force_scale = max(float(np.linalg.norm(F_ext)), float(np.linalg.norm(F_int)), 1.0)
    return force_scale, force_scale / stiffness_diag_scale(K)
