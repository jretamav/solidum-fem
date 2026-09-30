"""Backend Cholesky disperso (CHOLMOD vía ``scikit-sparse``).

Especializado en matrices simétricas y positivas definidas (SPD). ~2× más
rápido y ~2× menos memoria que LU para el mismo problema. Aborta si la matriz
no es SPD; los solvers degradan a LU automáticamente (ADR 0003 §5).

El import de ``sksparse.cholmod`` es responsabilidad del módulo: si la
dependencia no está instalada, este archivo lanza ``ImportError`` al
importarse y ``solidum.math.linalg.__init__`` lo absorbe silenciosamente.

Dos generaciones de la API de ``scikit-sparse``
-----------------------------------------------
La 0.5 rehízo la interfaz: ``cholesky(A)`` ya no devuelve un objeto
factorizado sino las matrices ``(R, p)``, y el objeto reutilizable pasa a
ser ``CholeskyFactor`` (``factorize``/``solve``). El extra ``cholmod``
admite ambas (``scikit-sparse>=0.4``), así que el backend también:

- **0.4.x** — ``cholesky(A)`` → ``Factor``; resolver es llamarlo,
  ``factor(b)``; ``factor.D()`` da los pivotes ``d_i`` de ``LDLᵀ`` (los
  ``L_ii²`` si la factorización es ``LLᵀ``) sin copiar el factor. Un pivote
  no positivo lanza ``CholmodNotPositiveDefiniteError`` al factorizar.
- **0.5.x** — ``CholeskyFactor(A).factorize(A)``; ``solve(b)``. La diagonal
  del factor es una vista (``factor.factor``): ``L_ii`` en ``LLᵀ`` (se
  elevan al cuadrado) o ``D_ii`` en ``LDLᵀ``. Un pivote no positivo lanza al
  factorizar, como en 0.4; pero uno *diminuto* no lanza hasta ``solve``, y
  sólo si su estimación de condición ``rcond = (min L_ii/max L_ii)²`` es
  exactamente 0. Esa comprobación se adelanta aquí a la factorización para
  que el fallback a LU actúe en el mismo sitio con las dos versiones.

Pivotes nulos (ADR 0019, capa 3; deuda #28 saldada el 2026-09-30)
-----------------------------------------------------------------
:attr:`CholeskyFactorized.n_zero_pivots` cuenta los pivotes de ``LDLᵀ`` por
debajo de ``ZERO_PIVOT_RTOL·max|K|``, la misma regla que SuperLU (la diagonal
de ``U``). Para una matriz simétrica, los pivotes de ``LDLᵀ`` son los de la
eliminación gaussiana sin pivoteo con la misma ordenación (``d_i = U_ii``);
los valores concretos dependen de la ordenación de cada backend, pero el
número de pivotes nulos de una matriz singular (su defecto de rango) no.
Pardiso aplica el umbral a la matriz escalada, así que la regla de SuperLU y
Cholesky no es invariante a un cambio de unidades y la de Pardiso sí (ADR
0019, «Lo que queda sin cubrir»). Sin este conteo, CHOLMOD factorizaba una
matriz singular con un pivote positivo diminuto (un mecanismo interno que la
carga no activa) y la solución pasaba sin aviso.
"""
from __future__ import annotations

import warnings

import numpy as np
import scipy.sparse as sp

from solidum.constants import ZERO_PIVOT_RTOL
from sksparse.cholmod import (  # type: ignore[import-not-found]
    CholmodNotPositiveDefiniteError,
)

try:  # scikit-sparse >= 0.5
    from sksparse.cholmod import (  # type: ignore[import-not-found]
        CholeskyFactor as _CholeskyFactor,
        CholmodWarning as _CholmodWarning,
    )
    _cholmod_cholesky = None
except ImportError:  # scikit-sparse 0.4.x
    from sksparse.cholmod import (  # type: ignore[import-not-found]
        cholesky as _cholmod_cholesky,
    )
    _CholeskyFactor = None
    _CholmodWarning = Warning


class CholeskyNotPositiveDefiniteError(RuntimeError):
    """``K`` no es positiva definida; señal para que el solver degrade a LU.

    Encapsula la excepción nativa de CHOLMOD para no propagar la dependencia
    externa al resto del código.
    """


class CholeskyFactorized:
    """Factorización de Cholesky retenida (CHOLMOD).

    Por construcción, una factorización de Cholesky exitosa garantiza que
    ``K`` es positiva definida; por tanto el conteo de pivotes negativos es
    siempre 0. Los pivotes **nulos** sí pueden aparecer: una matriz singular
    semidefinida puede factorizarse con un pivote positivo del orden del
    redondeo (ver :attr:`n_zero_pivots`).

    Parameters
    ----------
    solve
        ``b ↦ x``: resolución con el factor de CHOLMOD.
    pivots
        Callable sin argumentos que devuelve los pivotes ``d_i`` de ``LDLᵀ``.
        Se llama sólo si se pide :attr:`n_zero_pivots`.
    k_max
        ``max|K|``, escala del criterio de pivote nulo.
    """

    n_negative_pivots: int | None = 0

    def __init__(self, solve, pivots, k_max: float = 0.0):
        self._solve = solve
        self._pivots = pivots
        self._k_max = float(k_max)
        self._n_zero_pivots: int | None = None

    @property
    def n_zero_pivots(self) -> int:
        """Pivotes numéricamente nulos: ``|d_i| < ZERO_PIVOT_RTOL·max|K|``
        con ``d_i`` los pivotes de ``LDLᵀ`` (ADR 0019, capa 3). Se calcula
        sólo cuando se pide (el análisis estático) y una vez por
        factorización; leer la diagonal no copia el factor."""
        if self._n_zero_pivots is None:
            if self._k_max <= 0.0:
                self._n_zero_pivots = 0
            else:
                d = np.abs(np.asarray(self._pivots(), dtype=float))
                self._n_zero_pivots = int(np.sum(d < ZERO_PIVOT_RTOL * self._k_max))
        return self._n_zero_pivots

    def solve(self, b: np.ndarray) -> np.ndarray:
        return self._solve(b)


def _pivots_v05(factor) -> np.ndarray:
    """Pivotes de ``LDLᵀ`` de un ``CholeskyFactor`` (scikit-sparse 0.5)."""
    d = np.asarray(factor.factor.diagonal(), dtype=float)
    return d * d if factor.is_ll else d


def _solver_v05(factor):
    """Resolución con un ``CholeskyFactor`` (scikit-sparse 0.5).

    La 0.5 avisa en cada ``solve`` si ``rcond < n·eps`` ("Matrix is nearly
    singular"). Esa estimación, ``(min L_ii/max L_ii)²``, salta con
    tangentes legítimamente mal condicionadas —cerca de un punto límite, o
    con rigideces muy distintas entre marcos y sólidos— y no distingue una
    matriz singular; la red de seguridad del ADR 0019 decide con su propio
    criterio (``n_zero_pivots``), así que ese aviso se silencia aquí."""
    def solve(b: np.ndarray) -> np.ndarray:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Matrix is nearly singular",
                                    category=_CholmodWarning)
            return factor.solve(np.ascontiguousarray(b, dtype=np.float64))
    return solve


class CholeskySolver:
    """Cholesky disperso (``L·Lᵀ``) con permutación AMD/METIS automática.

    Solo aplicable si ``K`` es simétrica y positiva definida.
    """

    name = "cholesky"

    def solve(self, K: sp.spmatrix, b: np.ndarray) -> np.ndarray:
        return self.factorize(K).solve(b)

    def factorize(self, K: sp.spmatrix) -> CholeskyFactorized:
        K_csc = K.tocsc() if not sp.isspmatrix_csc(K) else K
        n = K_csc.shape[0]
        if n == 0:
            # Sistema vacío (todos los DOF prescritos), como en Pardiso: sin
            # incógnitas no hay nada que factorizar ni pivotes que contar.
            return CholeskyFactorized(lambda b: np.zeros(0), lambda: np.zeros(0))
        k_max = float(abs(K_csc).max()) if K_csc.nnz else 0.0
        try:
            if _CholeskyFactor is not None:
                factor = _CholeskyFactor(K_csc).factorize(K_csc)
                if factor.rcond == 0.0:
                    # 0.5 lo detectaría al resolver; se adelanta para que el
                    # fallback a LU actúe al factorizar, como con 0.4.
                    raise CholmodNotPositiveDefiniteError(
                        "Matrix is indefinite or singular to working precision.")
                return CholeskyFactorized(_solver_v05(factor),
                                          lambda: _pivots_v05(factor), k_max)
            factor = _cholmod_cholesky(K_csc)
        except CholmodNotPositiveDefiniteError as exc:
            raise CholeskyNotPositiveDefiniteError(str(exc)) from exc
        return CholeskyFactorized(factor, factor.D, k_max)
