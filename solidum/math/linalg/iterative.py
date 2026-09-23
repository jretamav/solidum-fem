"""Backend iterativo de Krylov para matrices simétricas (ADR 0018).

Resuelve ``K·x = b`` sin factorizar. Su ventaja no es la velocidad —por
debajo de ~10⁵ grados de libertad los directos multihilo ganan— sino la
**memoria**: un solver directo almacena los factores ``L``/``U``, cuyo relleno
crece con la talla (pico de Pardiso ×25 sobre ``K`` a 2·10⁵ DOF, medido); el
iterativo sólo guarda ``K``, el precondicionador y unos pocos vectores. Es la
misma división que hacen ANSYS (``SPARSE`` frente a ``PCG``) y Abaqus
(*direct* frente a *iterative*): el directo por defecto, el iterativo como
opción para modelos grandes y bien condicionados. Aquí también es **sólo bajo
petición** (``linear_algebra: iterative``); el despachador nunca lo elige solo.

Métodos
-------
- **Gradiente conjugado precondicionado (CG)** para matrices simétricas
  positivas definidas. Implementado aquí, no con ``scipy.sparse.linalg.cg``,
  para detectar la **curvatura negativa** (``pᵀKp ≤ 0``): es la prueba de que
  ``K`` no es definida positiva, la misma información que da Cholesky al
  fallar, y el backend responde pasando a MINRES en vez de devolver basura.
- **MINRES** para matrices simétricas indefinidas (arc-length cerca de un
  punto límite, régimen de ablandamiento). Precondicionado con Jacobi sobre
  ``|diag K|``, porque MINRES exige precondicionador definido positivo y un
  AMG construido sobre una matriz indefinida no lo garantiza.
- Las matrices **no simétricas** no llegan aquí: el despachador las envía a
  un solver directo (ADR 0018 §2), igual que hacen los programas comerciales.

Precondicionadores
------------------
- ``"amg"``: multimalla algebraico por agregación suavizada (``pyamg``) con
  los **modos de cuerpo rígido** del modelo como casi-núcleo
  (``StiffnessProperties.near_nullspace``). Con ellos el número de iteraciones
  apenas depende de la talla (9, 12 y 11 iteraciones a 26k, 86k y 2·10⁵ DOF,
  medido); sin ellos AMG rinde peor que no precondicionar.
- ``"jacobi"``: diagonal. Barato y casi inútil en elasticidad (medido: iguala
  a no precondicionar); se ofrece para diagnóstico.
- ``"none"``: identidad.
- ``"auto"`` (default): AMG si ``pyamg`` está instalado; si no, ninguno, con
  un aviso único.

Criterio de parada
------------------
``‖b − K·x‖ ≤ rtol·‖b‖`` sobre el residuo **verdadero**, recalculado al
terminar. El residuo recursivo de CG deriva del verdadero por redondeo; si al
converger el recursivo el verdadero no cumple, CG se reinicia desde el
iterado con el residuo recalculado (*residual replacement*, hasta
``ITERATIVE_MAX_RESTARTS`` veces). Con ``rtol = 1e-10`` por defecto
(``ITERATIVE_RTOL``), cinco órdenes por debajo de la tolerancia del Newton,
la sucesión de iterados del Newton es la del Newton exacto a efectos
prácticos. Si no se alcanza, se lanza :class:`IterativeNotConvergedError`:
**nunca** se devuelve una solución que no cumple la tolerancia.

El contrato ``factorize`` / ``solve`` del ADR 0003 se conserva: "factorizar"
construye el precondicionador (la parte cara y reutilizable, que el Newton
modificado congela junto con la matriz) y cada ``solve`` es una iteración de
Krylov.
"""
from __future__ import annotations

import warnings
from typing import Callable

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from solidum.constants import (
    AMG_MAX_COARSE,
    ITERATIVE_MAX_ITER,
    ITERATIVE_MAX_RESTARTS,
    ITERATIVE_RTOL,
)
from solidum.logging import get_logger

try:  # dependencia opcional: pip install solidum-fem[iterative]
    import pyamg  # type: ignore[import-not-found]
    HAS_PYAMG = True
except ModuleNotFoundError:
    pyamg = None  # type: ignore[assignment]
    HAS_PYAMG = False

_log = get_logger("linalg.iterative")

PRECONDITIONERS: tuple[str, ...] = ("auto", "amg", "jacobi", "none")

# Piso relativo para la diagonal de Jacobi: una entrada nula o casi nula de
# ``|diag K|`` se sustituye por este múltiplo de la mayor. Mantiene el
# precondicionador definido positivo y evita que ``1/|d|`` se dispare en una
# matriz indefinida cuya diagonal pasa cerca de cero (con 1e-14 el
# precondicionador amplificaba esas filas y MINRES no convergía).
_JACOBI_FLOOR_REL = 1.0e-8


class IterativeNotConvergedError(RuntimeError):
    """El método de Krylov no alcanzó la tolerancia.

    Es ``RuntimeError`` para que los solvers la traten como cualquier fallo
    del sistema lineal (el Newton abandona el paso), pero con su propio tipo
    y mensaje, para que el diagnóstico no la confunda con una tangente
    singular: el corrector la registra aparte (ADR 0018 §5).
    """

    def __init__(self, *, method: str, preconditioner: str, iterations: int,
                 relative_residual: float, rtol: float, size: int):
        self.method = method
        self.preconditioner = preconditioner
        self.iterations = int(iterations)
        self.relative_residual = float(relative_residual)
        self.rtol = float(rtol)
        self.size = int(size)
        super().__init__(
            f"Solver iterativo ({method}, precondicionador '{preconditioner}') "
            f"sin convergencia tras {self.iterations} iteraciones sobre "
            f"{self.size} incógnitas: residuo relativo {self.relative_residual:.2e} "
            f"frente a la tolerancia {self.rtol:.0e}. Causas habituales: "
            f"matriz mal condicionada (láminas, vigas esbeltas, elementos muy "
            f"distorsionados, penalizaciones rígidas) o precondicionador "
            f"insuficiente. Salidas: linear_algebra 'auto' (solver directo) "
            f"o 'iterative:amg' si no se estaba usando (requiere pyamg)."
        )


class _NegativeCurvature(Exception):
    """CG encontró ``pᵀKp ≤ 0`` (o ``rᵀz ≤ 0``): el operador o el
    precondicionador no son definidos positivos."""


_warned: set[str] = set()


def _warn_once(key: str, message: str) -> None:
    if key in _warned:
        return
    _warned.add(key)
    warnings.warn(message, RuntimeWarning, stacklevel=4)


# ----------------------------------------------------------------------
# Precondicionadores
# ----------------------------------------------------------------------

def _identity(r: np.ndarray) -> np.ndarray:
    return r


def _jacobi(K: sp.csr_matrix) -> Callable[[np.ndarray], np.ndarray]:
    """Inversa de ``|diag K|`` con piso relativo: siempre definida positiva."""
    d = np.abs(K.diagonal())
    dmax = float(d.max()) if d.size else 0.0
    floor = _JACOBI_FLOOR_REL * dmax if dmax > 0.0 else 1.0
    inv = 1.0 / np.maximum(d, floor)
    return lambda r: inv * r


def _amg(K: sp.csr_matrix, near_nullspace: np.ndarray | None):
    """Jerarquía AMG por agregación suavizada, lista para aplicarse.

    Dos decisiones medidas (ADR 0018 §3):

    - **Nivel grueso con factorización dispersa** (``coarse_solver="splu"``),
      no con la pseudoinversa densa por defecto de ``pyamg``. Con seis modos
      rígidos cada agregado aporta seis incógnitas gruesas, así que la
      jerarquía puede detenerse con miles de incógnitas en el último nivel
      (2058 en ``Hex8 20³``); su pseudoinversa es una SVD densa O(n³) que
      costaba ~3 s frente a ~0,1 s de ``splu``.
    - **Construcción completa aquí**, no perezosa: ``pyamg`` factoriza el
      nivel grueso en la primera aplicación del precondicionador. Se fuerza
      con una aplicación de calentamiento para que todo el coste quede en
      ``factorize`` —que es lo que el Newton modificado reutiliza— y no
      aparezca escondido en la primera resolución.
    """
    ml = pyamg.smoothed_aggregation_solver(
        K, B=near_nullspace, max_coarse=AMG_MAX_COARSE, coarse_solver="splu",
    )
    M = ml.aspreconditioner(cycle="V")
    M.matvec(np.ones(K.shape[0]))
    return lambda r: M.matvec(r)


# ----------------------------------------------------------------------
# Métodos de Krylov
# ----------------------------------------------------------------------

def _pcg(K: sp.csr_matrix, r0: np.ndarray, apply_M, tol_abs: float,
         max_iter: int) -> tuple[np.ndarray, int, bool]:
    """Resuelve ``K·d = r0`` desde ``d = 0`` con CG precondicionado.

    Devuelve ``(d, iteraciones, convergió_el_recursivo)``. Lanza
    :class:`_NegativeCurvature` si ``pᵀKp ≤ 0`` o ``rᵀz ≤ 0``.
    """
    d = np.zeros_like(r0)
    r = r0.copy()
    if np.linalg.norm(r) <= tol_abs:
        return d, 0, True
    z = apply_M(r)
    rz = float(r @ z)
    if rz <= 0.0:
        raise _NegativeCurvature
    p = z.copy()
    for k in range(1, max_iter + 1):
        Kp = K @ p
        pKp = float(p @ Kp)
        if pKp <= 0.0:
            raise _NegativeCurvature
        alpha = rz / pKp
        d += alpha * p
        r -= alpha * Kp
        if np.linalg.norm(r) <= tol_abs:
            return d, k, True
        z = apply_M(r)
        rz_new = float(r @ z)
        if rz_new <= 0.0:
            raise _NegativeCurvature
        p = z + (rz_new / rz) * p
        rz = rz_new
    return d, max_iter, False


class _EmptyIterative:
    """Sistema ``0×0`` (todos los DOF prescritos): vector vacío."""

    n_negative_pivots: int | None = None

    def solve(self, b: np.ndarray) -> np.ndarray:
        return np.zeros(0)


class IterativeFactorized:
    """Precondicionador construido sobre ``K``; ``solve(b)`` itera.

    Arranca con CG. Si CG detecta curvatura negativa pasa a MINRES para el
    resto de resoluciones de esta "factorización" (la matriz es la misma,
    así que la conclusión también). Si CG agota las iteraciones sin
    converger, intenta MINRES antes de rendirse: sobre una matriz indefinida
    CG puede estancarse sin llegar a exhibir ``pᵀKp ≤ 0``.

    ``n_negative_pivots`` es ``None``: un método de Krylov no da la inercia.
    """

    n_negative_pivots: int | None = None

    def __init__(self, K: sp.csr_matrix, apply_M, preconditioner: str,
                 rtol: float, max_iter: int, method: str = "cg"):
        self._K = K
        self._apply_M = apply_M
        self.preconditioner = preconditioner
        self.rtol = float(rtol)
        self.max_iter = int(max_iter)
        self.method = method
        self._minres_M = None
        # Telemetría de la última resolución.
        self.last_iterations = 0
        self.last_relative_residual = 0.0

    # --------------------------------------------------------------
    def solve(self, b: np.ndarray) -> np.ndarray:
        b = np.ascontiguousarray(b, dtype=np.float64)
        n = self._K.shape[0]
        b_norm = float(np.linalg.norm(b))
        if b_norm == 0.0:
            self.last_iterations, self.last_relative_residual = 0, 0.0
            return np.zeros(n)

        if self.method == "cg":
            try:
                return self._solve_cg(b, b_norm)
            except _NegativeCurvature:
                _log.info(
                    "Solver iterativo: CG detectó curvatura negativa (la matriz "
                    "no es definida positiva). Se continúa con MINRES."
                )
                self.method = "minres"
            except IterativeNotConvergedError as exc:
                _log.info(
                    f"Solver iterativo: CG no convergió ({exc.iterations} it., "
                    f"residuo {exc.relative_residual:.1e}); se intenta MINRES."
                )
                try:
                    x = self._solve_minres(b, b_norm)
                except IterativeNotConvergedError:
                    raise exc from None
                # MINRES resolvió lo que CG no: la matriz es la misma en las
                # siguientes resoluciones, así que se sigue con MINRES.
                self.method = "minres"
                return x
        return self._solve_minres(b, b_norm)

    # --------------------------------------------------------------
    def _finish(self, x, b, b_norm, iterations, method, precond):
        rel = float(np.linalg.norm(b - self._K @ x)) / b_norm
        self.last_iterations, self.last_relative_residual = iterations, rel
        if rel <= self.rtol:
            return x
        raise IterativeNotConvergedError(
            method=method, preconditioner=precond, iterations=iterations,
            relative_residual=rel, rtol=self.rtol, size=self._K.shape[0],
        )

    def _solve_cg(self, b: np.ndarray, b_norm: float) -> np.ndarray:
        K = self._K
        tol_abs = self.rtol * b_norm
        x = np.zeros_like(b)
        total = 0
        for _ in range(ITERATIVE_MAX_RESTARTS + 1):
            r = b - K @ x
            if float(np.linalg.norm(r)) <= tol_abs:
                break
            d, it, _ = _pcg(K, r, self._apply_M, tol_abs, self.max_iter - total)
            x += d
            total += it
            if total >= self.max_iter:
                break
        return self._finish(x, b, b_norm, total, "CG", self.preconditioner)

    def _solve_minres(self, b: np.ndarray, b_norm: float) -> np.ndarray:
        K = self._K
        n = K.shape[0]
        if self._minres_M is None:
            jac = _jacobi(K)
            self._minres_M = spla.LinearOperator((n, n), matvec=jac, dtype=np.float64)
        tol_abs = self.rtol * b_norm
        x = np.zeros_like(b)
        total = 0
        for _ in range(ITERATIVE_MAX_RESTARTS + 1):
            r = b - K @ x
            r_norm = float(np.linalg.norm(r))
            if r_norm <= tol_abs:
                break
            counter = [0]

            def _count(_xk):
                counter[0] += 1

            # El criterio interno de MINRES es relativo a ‖A‖‖x‖ + ‖b‖ en la
            # norma del precondicionador, más laxo que el nuestro: se pide
            # una reducción del residuo de la corrección un orden más fuerte
            # y el resultado se verifica siempre contra el residuo verdadero.
            d, _info = spla.minres(
                K, r, rtol=0.1 * tol_abs / r_norm, maxiter=self.max_iter - total,
                M=self._minres_M, callback=_count,
            )
            x += d
            total += counter[0]
            if total >= self.max_iter:
                break
        return self._finish(x, b, b_norm, total, "MINRES", "jacobi |diag|")


class IterativeSolver:
    """Backend iterativo (ADR 0018): CG / MINRES con AMG, Jacobi o nada.

    Parameters
    ----------
    props
        Propiedades del sistema; de ellas sólo se usa ``near_nullspace``
        (modos de cuerpo rígido para AMG). ``None`` ⇒ AMG con el casi-núcleo
        por defecto de ``pyamg`` (el vector constante).
    preconditioner
        ``"auto"`` (AMG si ``pyamg`` está instalado, si no ninguno),
        ``"amg"``, ``"jacobi"`` o ``"none"``.
    rtol, max_iter
        Tolerancia sobre el residuo relativo verdadero y cota de iteraciones
        por resolución (``ITERATIVE_RTOL``, ``ITERATIVE_MAX_ITER``).
    """

    name = "iterative"
    # El despachador construye este backend con las propiedades del sistema
    # (necesita el casi-núcleo); el resto de backends no las reciben.
    ACCEPTS_PROPERTIES = True

    def __init__(self, props=None, *, preconditioner: str = "auto",
                 rtol: float = ITERATIVE_RTOL, max_iter: int = ITERATIVE_MAX_ITER):
        preconditioner = str(preconditioner).lower()
        if preconditioner not in PRECONDITIONERS:
            raise ValueError(
                f"IterativeSolver: precondicionador '{preconditioner}' desconocido. "
                f"Disponibles: {list(PRECONDITIONERS)}."
            )
        if preconditioner == "amg" and not HAS_PYAMG:
            raise ValueError(
                "IterativeSolver: el precondicionador 'amg' requiere pyamg, que "
                "no está instalado. Instálalo con `pip install solidum-fem[iterative]` "
                "(en Windows con Python sin binarios de pyamg, ver el anexo de la "
                "capa algebraica del manual) o usa 'iterative:none'."
            )
        if not (0.0 < float(rtol) < 1.0):
            raise ValueError(f"IterativeSolver: rtol={rtol} fuera de (0, 1).")
        if int(max_iter) < 1:
            raise ValueError(f"IterativeSolver: max_iter={max_iter} debe ser ≥ 1.")
        self.props = props
        self.preconditioner = preconditioner
        self.rtol = float(rtol)
        self.max_iter = int(max_iter)

    def _resolved_preconditioner(self) -> str:
        if self.preconditioner != "auto":
            return self.preconditioner
        if HAS_PYAMG:
            return "amg"
        _warn_once(
            "auto-no-pyamg",
            "pyamg no está instalado: el solver iterativo trabajará sin "
            "precondicionador, y en elasticidad eso cuesta del orden de cientos "
            "de iteraciones por resolución. Para activar AMG: "
            "`pip install solidum-fem[iterative]`.",
        )
        return "none"

    def _near_nullspace(self, n: int) -> np.ndarray | None:
        provider = getattr(self.props, "near_nullspace", None)
        if provider is None:
            return None
        B = provider()
        if B is None:
            return None
        B = np.ascontiguousarray(np.asarray(B, dtype=np.float64))
        if B.ndim != 2 or B.shape[0] != n or B.shape[1] == 0:
            _warn_once(
                f"nullspace-shape-{B.shape}-{n}",
                f"IterativeSolver: casi-núcleo de forma {B.shape} incompatible "
                f"con el sistema de {n} incógnitas; AMG usa su casi-núcleo por "
                f"defecto (vector constante), que en elasticidad converge peor.",
            )
            return None
        return B

    def solve(self, K: sp.spmatrix, b: np.ndarray) -> np.ndarray:
        return self.factorize(K).solve(b)

    def factorize(self, K: sp.spmatrix):
        K = K.tocsr() if not sp.isspmatrix_csr(K) else K
        if K.dtype != np.float64:
            K = K.astype(np.float64)
        n = K.shape[0]
        if n == 0:
            return _EmptyIterative()

        # Una matriz simétrica definida positiva tiene la diagonal
        # estrictamente positiva (d_i = e_iᵀ K e_i > 0). Si alguna entrada es
        # ≤ 0 la matriz no es SPD: CG no aplica y AMG no se puede construir
        # (invierte la diagonal). Se va directamente a MINRES.
        if np.any(K.diagonal() <= 0.0):
            _log.info(
                "Solver iterativo: la diagonal tiene entradas ≤ 0, así que la "
                "matriz no es definida positiva. Se usa MINRES."
            )
            return IterativeFactorized(K, _identity, "jacobi |diag|", self.rtol,
                                       self.max_iter, method="minres")

        name = self._resolved_preconditioner()
        if name == "amg":
            try:
                apply_M = _amg(K, self._near_nullspace(n))
            except Exception as exc:  # pyamg señala el fallo con errores propios
                _warn_once(
                    f"amg-setup-{type(exc).__name__}",
                    f"IterativeSolver: no se pudo construir el precondicionador "
                    f"AMG ({type(exc).__name__}: {exc}). Se usa Jacobi para esta "
                    f"matriz; la solución sigue verificándose contra la tolerancia.",
                )
                name, apply_M = "jacobi", _jacobi(K)
        elif name == "jacobi":
            apply_M = _jacobi(K)
        else:
            apply_M = _identity
        return IterativeFactorized(K, apply_M, name, self.rtol, self.max_iter)
