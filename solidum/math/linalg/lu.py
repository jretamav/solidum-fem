"""Backend LU disperso (SuperLU vía ``scipy.sparse.linalg``)."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from solidum.constants import ZERO_PIVOT_RTOL
from solidum.math.linalg.base import out_of_memory_error


class LUFactorized:
    """Factorización LU retenida (``scipy.sparse.linalg.SuperLU``).

    Permite ``solve(b)`` repetido sin refactorizar — base para Newton
    modificado y reuso entre pasos en dinámica implícita.

    El conteo de pivotes negativos no se expone (``None``): SuperLU usa
    pivoteo parcial por filas que **no preserva la inercia** de la matriz,
    por lo que extraer ``sign(diag(U))`` no equivale al número de
    autovalores negativos. Para Sturm sequence verdadero ver ``LDLTSolver``.
    """

    n_negative_pivots: int | None = None

    def __init__(self, lu: spla.SuperLU, k_max: float = 0.0):
        self._lu = lu
        self._k_max = float(k_max)
        self._n_zero_pivots: int | None = None

    @property
    def n_zero_pivots(self) -> int:
        """Pivotes numéricamente nulos: ``|U_ii| < ZERO_PIVOT_RTOL·max|K|``
        (ADR 0019, capa 3). En una matriz regular es 0; > 0 delata una
        matriz singular aunque la solución parezca razonable (un mecanismo
        interno sin carga que lo active). SuperLU no lo reporta, así que se
        lee la diagonal de ``U``, lo que obliga a SciPy a copiar el factor:
        se calcula sólo cuando se pide (el análisis estático lineal) y una
        vez por factorización."""
        if self._n_zero_pivots is None:
            if self._k_max <= 0.0:
                self._n_zero_pivots = 0
            else:
                d = np.abs(self._lu.U.diagonal())
                self._n_zero_pivots = int(np.sum(d < ZERO_PIVOT_RTOL * self._k_max))
        return self._n_zero_pivots

    def solve(self, b: np.ndarray) -> np.ndarray:
        return self._lu.solve(b)


class LUSolver:
    """Factorización LU general con pivoteo parcial. Backend universal.

    No asume simetría ni positividad. Es el fallback al que se degrada cuando
    Cholesky o LDLᵀ detectan pérdida de positividad.
    """

    name = "lu"

    def solve(self, K: sp.spmatrix, b: np.ndarray) -> np.ndarray:
        # ``spsolve`` ante una matriz singular emite ``MatrixRankWarning`` y
        # devuelve ``NaN`` sin lanzar; ``splu`` sí lanza ``RuntimeError``
        # ("Factor is exactly singular"), que es lo que los solvers no
        # lineales capturan para diagnosticar tangente singular (ADR 0011).
        # Una singularidad numérica (pivotes ~1e-300) tampoco lanza, así que
        # se vigila además que la solución sea finita.
        x = self.factorize(K).solve(b)
        if not np.all(np.isfinite(x)):
            raise RuntimeError(
                "LUSolver: la solución no es finita (matriz singular o casi "
                "singular)."
            )
        return x

    def factorize(self, K: sp.spmatrix) -> LUFactorized:
        K_csc = K.tocsc() if not sp.isspmatrix_csc(K) else K
        try:
            k_max = float(abs(K_csc).max()) if K_csc.nnz else 0.0
            return LUFactorized(spla.splu(K_csc), k_max=k_max)
        except MemoryError:
            raise out_of_memory_error(K_csc.shape[0], "SuperLU") from None
        except RuntimeError as exc:
            if "memory" in str(exc).lower():
                raise out_of_memory_error(K_csc.shape[0], "SuperLU") from exc
            raise
