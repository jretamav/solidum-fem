"""Backend LU disperso (SuperLU vía ``scipy.sparse.linalg``)."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from solidum.constants import SUSPECT_PIVOT_RTOL
from solidum.math.linalg.base import out_of_memory_error, pivot_ratios


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

    def __init__(self, lu: spla.SuperLU, diag: np.ndarray | None = None):
        self._lu = lu
        self._diag = np.zeros(0) if diag is None else np.abs(np.asarray(diag, dtype=float))
        self._n_small_pivots: int | None = None

    @property
    def n_small_pivots(self) -> int:
        """Pivotes pequeños frente a la diagonal de su fila y su columna,
        ``|U_ii|/√(|K_rr|·|K_cc|) < SUSPECT_PIVOT_RTOL`` (ADR 0019, ampliación
        del 2026-10-20). Es el filtro de la red de seguridad: con alguno, el
        análisis estático busca un mecanismo interno por la energía del modo
        (``model_checks.find_internal_mechanisms``). El cociente no depende
        de las unidades. SuperLU no reporta sus pivotes, así que se lee la
        diagonal de ``U``, lo que obliga a SciPy a copiar el factor: se
        calcula sólo cuando se pide y una vez por factorización.

        ``LU = Pr·K·Pc``: el pivote ``i`` está en la fila ``argsort(perm_r)[i]``
        y la columna ``argsort(perm_c)[i]`` de ``K``."""
        if self._n_small_pivots is None:
            if self._diag.size == 0:
                self._n_small_pivots = 0
            else:
                ratios = pivot_ratios(
                    self._lu.U.diagonal(), self._diag,
                    rows=np.argsort(self._lu.perm_r), cols=np.argsort(self._lu.perm_c))
                self._n_small_pivots = int(np.sum(ratios < SUSPECT_PIVOT_RTOL))
        return self._n_small_pivots

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
            return LUFactorized(spla.splu(K_csc), diag=K_csc.diagonal())
        except MemoryError:
            raise out_of_memory_error(K_csc.shape[0], "SuperLU") from None
        except RuntimeError as exc:
            if "memory" in str(exc).lower():
                raise out_of_memory_error(K_csc.shape[0], "SuperLU") from exc
            raise
