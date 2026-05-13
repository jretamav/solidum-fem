"""Backend de autovalor generalizado simétrico ``K·φ = λ·M·φ`` (ADR 0009 §6).

Envuelve :func:`scipy.sparse.linalg.eigsh` (ARPACK Lanczos) con la política por
defecto que necesita el análisis modal: shift-invert centrado en ``sigma`` para
extraer los autovalores más cercanos al shift (`which="LM"` sobre el operador
desplazado equivale a `cercanos a sigma` en el espectro original).

No comparte el ``Protocol`` ``LinearAlgebraSolver`` de los backends de
``K·x = b``: la firma natural aquí es ``solve(K, M, n_modes) → (λ, φ)``, no
``solve(K, b) → x``. Mantenerlos en capas separadas evita contaminar la
abstracción de :mod:`fenix.math.linalg.base`.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


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
        if not (1 <= n_modes < n):
            raise ValueError(
                f"EigenSolver: n_modes={n_modes} fuera de rango. Debe cumplir "
                f"1 ≤ n_modes < {n} (= número de DOFs libres)."
            )

        eigenvalues, eigenvectors = spla.eigsh(
            K,
            k=n_modes,
            M=M,
            sigma=self.sigma,
            which=self.which,
            tol=self.tol,
        )

        order = np.argsort(eigenvalues)
        return eigenvalues[order], eigenvectors[:, order]
