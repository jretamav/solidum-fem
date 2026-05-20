"""Backend Cholesky disperso (CHOLMOD vía ``scikit-sparse``).

Especializado en matrices simétricas y positivas definidas (SPD). ~2× más
rápido y ~2× menos memoria que LU para el mismo problema. Aborta si la matriz
no es SPD; el solver no lineal degrada a LU automáticamente.

El import de ``sksparse.cholmod`` es responsabilidad del módulo: si la
dependencia no está instalada, este archivo lanza ``ImportError`` al
importarse y ``solidum.math.linalg.__init__`` lo absorbe silenciosamente.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from sksparse.cholmod import (  # type: ignore[import-not-found]
    CholmodNotPositiveDefiniteError,
    cholesky as _cholmod_cholesky,
)


class CholeskyNotPositiveDefiniteError(RuntimeError):
    """``K`` no es positiva definida; señal para que el solver no lineal degrade a LU.

    Encapsula la excepción nativa de CHOLMOD para no propagar la dependencia
    externa al resto del código.
    """


class CholeskyFactorized:
    """Factorización ``L·Lᵀ`` retenida (CHOLMOD ``Factor``).

    Por construcción, una factorización de Cholesky exitosa garantiza que
    ``K`` es positiva definida; por tanto el conteo de pivotes negativos es
    siempre 0.
    """

    n_negative_pivots: int | None = 0

    def __init__(self, factor):
        self._factor = factor

    def solve(self, b: np.ndarray) -> np.ndarray:
        return self._factor(b)


class CholeskySolver:
    """Cholesky disperso (``L·Lᵀ``) con permutación AMD/METIS automática.

    Solo aplicable si ``K`` es simétrica y positiva definida.
    """

    name = "cholesky"

    def solve(self, K: sp.spmatrix, b: np.ndarray) -> np.ndarray:
        return self.factorize(K).solve(b)

    def factorize(self, K: sp.spmatrix) -> CholeskyFactorized:
        K_csc = K.tocsc() if not sp.isspmatrix_csc(K) else K
        try:
            factor = _cholmod_cholesky(K_csc)
        except CholmodNotPositiveDefiniteError as exc:
            raise CholeskyNotPositiveDefiniteError(str(exc)) from exc
        return CholeskyFactorized(factor)
