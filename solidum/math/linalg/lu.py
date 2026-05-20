"""Backend LU disperso (SuperLU vía ``scipy.sparse.linalg``)."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


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

    def __init__(self, lu: spla.SuperLU):
        self._lu = lu

    def solve(self, b: np.ndarray) -> np.ndarray:
        return self._lu.solve(b)


class LUSolver:
    """Factorización LU general con pivoteo parcial. Backend universal.

    No asume simetría ni positividad. Es el fallback al que se degrada cuando
    Cholesky o LDLᵀ detectan pérdida de positividad.
    """

    name = "lu"

    def solve(self, K: sp.spmatrix, b: np.ndarray) -> np.ndarray:
        return spla.spsolve(K, b)

    def factorize(self, K: sp.spmatrix) -> LUFactorized:
        K_csc = K.tocsc() if not sp.isspmatrix_csc(K) else K
        return LUFactorized(spla.splu(K_csc))
