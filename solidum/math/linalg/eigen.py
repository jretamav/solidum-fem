"""Backend de autovalor generalizado simétrico ``K·φ = λ·M·φ`` (ADR 0009 §6).

Envuelve :func:`scipy.sparse.linalg.eigsh` (ARPACK Lanczos) con la política por
defecto que necesita el análisis modal: shift-invert centrado en ``sigma`` para
extraer los autovalores más cercanos al shift (`which="LM"` sobre el operador
desplazado equivale a `cercanos a sigma` en el espectro original).

No comparte el ``Protocol`` ``LinearAlgebraSolver`` de los backends de
``K·x = b``: la firma natural aquí es ``solve(K, M, n_modes) → (λ, φ)``, no
``solve(K, b) → x``. Mantenerlos en capas separadas evita contaminar la
abstracción de :mod:`solidum.math.linalg.base`.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from solidum.logging import get_logger

_log = get_logger("linalg.eigen")

# Shift negativo relativo a la escala espectral ``tr(K)/tr(M)`` usado como
# reintento cuando ``K`` es singular (modos rígidos, free-free) y el shift
# pedido es exactamente 0: ``K − σM`` con ``σ < 0`` es definida positiva
# y factoriza siempre; los modos rígidos salen como ``λ ≈ 0``.
_SINGULAR_RETRY_SHIFT_RTOL = 1.0e-6
# Tamaño (DOFs libres) por debajo del cual se usa el eigensolver denso.
_DENSE_FALLBACK_SIZE = 32


class EigenSolver:
    """Solver de autovalor generalizado simétrico, real y positivo (semi)definido.

    Parameters
    ----------
    sigma : float, default 0.0
        Shift para shift-invert: ARPACK factoriza ``(K − σ·M)`` y busca los
        autovalores del problema desplazado ``(K − σM)⁻¹M``, cuyos mayores
        en magnitud corresponden a los autovalores del problema original
        más cercanos a ``σ``. ``σ = 0`` recupera las frecuencias más bajas,
        que es lo relevante en ingeniería estructural.
    which : str, default "LM"
        Estrategia ARPACK. Con shift-invert el default ``"LM"`` es el
        adecuado (mayor magnitud sobre el operador desplazado = más cercano
        a σ en el original). Se permite override para diagnóstico.
    tol : float, default 1.0e-9
        Tolerancia ARPACK sobre el residuo del autovalor. ``0`` deja el
        default interno de ARPACK (≈ eps de la máquina).
    """

    name = "eigsh"

    def __init__(
        self,
        *,
        sigma: float = 0.0,
        which: str = "LM",
        tol: float = 1.0e-9,
    ) -> None:
        self.sigma = float(sigma)
        self.which = str(which)
        self.tol = float(tol)
        # ``False`` si ARPACK se detuvo sin converger todos los modos pedidos
        # (se devuelven los que convergieron). Lo consulta ``ModalSolver``.
        self.last_converged: bool = True

    def solve(
        self,
        K: sp.spmatrix,
        M: sp.spmatrix,
        n_modes: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Resuelve ``K·φ = λ·M·φ`` y devuelve ``(eigenvalues, eigenvectors)``.

        Los modos están **M-ortonormales** (``Φᵀ M Φ = I``) — ``eigsh`` lo
        garantiza al pasar ``M`` explícito y un único shift. Reordenados en
        sentido ascendente de autovalor.

        Parameters
        ----------
        K, M : scipy.sparse.spmatrix
            Matrices reducidas (tras aplicar Dirichlet). Simétricas; ``M``
            debe ser positiva definida estrictamente.
        n_modes : int
            Número de modos pedidos. Debe cumplir ``1 ≤ n_modes < K.shape[0]``
            por restricción de ARPACK.

        Returns
        -------
        eigenvalues : np.ndarray, shape (n_modes,)
            Autovalores ``λ_n = ω²_n`` en orden ascendente.
        eigenvectors : np.ndarray, shape (n_dof_red, n_modes)
            Modos ``φ_n`` en columnas, M-ortonormales.
        """
        n = K.shape[0]
        if not (1 <= n_modes <= n):
            raise ValueError(
                f"EigenSolver: n_modes={n_modes} fuera de rango. Debe cumplir "
                f"1 ≤ n_modes ≤ {n} (= número de DOFs libres)."
            )

        self.last_converged = True
        if n_modes >= n - 1 or n <= _DENSE_FALLBACK_SIZE:
            # ARPACK exige k < n y es inestable en problemas diminutos; para
            # modelos pequeños o al pedir (casi) todos los modos se resuelve el
            # problema generalizado denso. ``eigh(a, b)`` devuelve los modos
            # M-ortonormales (vᵀ·b·v = 1), la misma normalización que ARPACK.
            from scipy.linalg import eigh
            vals, vecs = eigh(np.asarray(K.todense()), np.asarray(M.todense()))
            return vals[:n_modes], vecs[:, :n_modes]
        sigma = self.sigma
        try:
            eigenvalues, eigenvectors = self._eigsh(K, M, n_modes, sigma)
        except RuntimeError as exc:
            # ``sigma = 0`` sobre ``K`` singular (estructura free-free con modos
            # rígidos): SuperLU no puede factorizar ``K − 0·M``. Que ocurra o
            # no dependía del redondeo de cada malla. Se reintenta con un
            # shift negativo pequeño, que deja el problema bien planteado y
            # conserva los modos de menor frecuencia como objetivo.
            if sigma != 0.0:
                raise
            tr_K = float(np.sum(np.abs(K.diagonal())))
            tr_M = float(np.sum(np.abs(M.diagonal())))
            sigma = -_SINGULAR_RETRY_SHIFT_RTOL * (tr_K / tr_M if tr_M > 0.0 else 1.0)
            _log.warning(
                f"EigenSolver: K − σM singular con σ=0 ({exc}); probable "
                f"estructura con modos de cuerpo rígido. Reintentando con "
                f"σ={sigma:.3e}."
            )
            eigenvalues, eigenvectors = self._eigsh(K, M, n_modes, sigma)

        order = np.argsort(eigenvalues)
        return eigenvalues[order], eigenvectors[:, order]

    def _eigsh(self, K, M, n_modes, sigma):
        try:
            return spla.eigsh(
                K, k=n_modes, M=M, sigma=sigma, which=self.which, tol=self.tol,
            )
        except spla.ArpackNoConvergence as exc:
            n_ok = int(np.size(exc.eigenvalues))
            if n_ok == 0:
                raise RuntimeError(
                    "EigenSolver: ARPACK no convergió ningún modo. Revise la "
                    "matriz de masa (¿DOF sin masa?) o aumente la tolerancia."
                ) from exc
            _log.warning(
                f"EigenSolver: ARPACK convergió sólo {n_ok} de {n_modes} modos "
                f"pedidos; se devuelven los convergidos (converged=False)."
            )
            self.last_converged = False
            return np.asarray(exc.eigenvalues), np.asarray(exc.eigenvectors)
