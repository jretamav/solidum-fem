"""Backend de autovalor generalizado simétrico ``K·φ = λ·M·φ`` (ADR 0009 §6).

Envuelve :func:`scipy.sparse.linalg.eigsh` (ARPACK Lanczos) con la política por
defecto que necesita el análisis modal: shift-invert centrado en ``sigma`` para
extraer los autovalores más cercanos al shift (`which="LM"` sobre el operador
desplazado equivale a `cercanos a sigma` en el espectro original).

No comparte el ``Protocol`` ``LinearAlgebraSolver`` de los backends de
``K·x = b``: la firma natural aquí es ``solve(K, M, n_modes) → (λ, φ)``, no
``solve(K, b) → x``. Mantenerlos en capas separadas evita contaminar la
abstracción de :mod:`solidum.math.linalg.base`.

Grados de libertad sin masa (ADR 0024): si ``M`` tiene diagonal nula en
algunos grados de libertad (masas sólo en ciertos nudos, elementos con
``density = 0``), ``M`` es semidefinida y ARPACK no construye la
factorización de Arnoldi. Esos grados de libertad se condensan
estáticamente, de forma exacta, y sus componentes de los modos se recuperan
después (:meth:`EigenSolver._solve_condensed`).
"""
from __future__ import annotations

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from solidum.logging import get_logger

_log = get_logger("linalg.eigen")


class SemidefiniteMassError(ValueError):
    """El problema ``K·φ = λ·M·φ`` con grados de libertad sin masa no tiene
    solución por condensación estática (ADR 0024).

    Capa algebraica: identifica los grados de libertad por su índice en las
    matrices recibidas. ``ModalSolver`` la traduce a
    :class:`~solidum.math.solvers.diagnostics.MasslessDofsError`, que los
    nombra como (nodo, dof).

    Attributes
    ----------
    cause : str
        ``"too_many_modes"`` (se piden más modos que grados de libertad con
        masa), ``"no_mass"`` (ninguno tiene masa), ``"massless_mechanism"``
        (los grados de libertad sin masa no tienen rigidez que los fije),
        ``"singular_mass"`` (la masa de los que sí la tienen es singular) o
        ``"indefinite_mass"`` (una fila con diagonal nula y fuera de la
        diagonal no nula: ``M`` no es semidefinida positiva).
    massless : np.ndarray
        Índices de los grados de libertad sin masa.
    n_mass : int
        Número de grados de libertad con masa (= frecuencias finitas).
    """

    def __init__(self, message: str, *, cause: str, massless, n_mass: int):
        super().__init__(message)
        self.cause = str(cause)
        self.massless = np.asarray(massless, dtype=np.int64)
        self.n_mass = int(n_mass)


def massless_dofs(M: sp.spmatrix) -> np.ndarray:
    """Índices de los grados de libertad sin masa: diagonal de ``M`` nula.

    Cero exacto, no relativo: en una matriz semidefinida positiva una
    diagonal nula implica fila y columna nulas, y entonces la condensación
    es exacta. Una masa pequeña pero no nula es masa (da frecuencias altas,
    no indeterminadas) y no se toca."""
    return np.flatnonzero(np.asarray(M.diagonal()) == 0.0)

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
        # Índices de los grados de libertad sin masa condensados en la
        # última resolución (vacío si ``M`` tenía masa en todos).
        self.last_massless: np.ndarray = np.empty(0, dtype=np.int64)

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
            semidefinida positiva. Si tiene grados de libertad sin masa
            (diagonal nula), se condensan (ADR 0024) y sólo hay tantas
            frecuencias finitas como grados de libertad con masa.
        n_modes : int
            Número de modos pedidos: ``1 ≤ n_modes ≤ K.shape[0]``, y no más
            que grados de libertad con masa.

        Raises
        ------
        SemidefiniteMassError
            Si hay grados de libertad sin masa y el problema no se puede
            resolver condensándolos (ver ``cause``).

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
        massless = massless_dofs(M)
        self.last_massless = massless
        if massless.size:
            return self._solve_condensed(sp.csr_matrix(K), sp.csr_matrix(M),
                                         n_modes, massless)
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

    # ------------------------------------------------------------------
    # Grados de libertad sin masa: condensación estática exacta (ADR 0024)
    # ------------------------------------------------------------------

    def _solve_condensed(self, K: sp.csr_matrix, M: sp.csr_matrix,
                         n_modes: int, zero: np.ndarray):
        """``K·φ = λ·M·φ`` con los grados de libertad ``zero`` sin masa.

        Con ``m`` los grados de libertad con masa y ``0`` los demás, la fila
        ``0`` no tiene inercia: ``φ₀ = −K₀₀⁻¹·K₀ₘ·φₘ`` (exacta). Queda
        ``K̂·φₘ = λ·Mₘₘ·φₘ`` con ``K̂ = Kₘₘ − Kₘ₀·K₀₀⁻¹·K₀ₘ``, cuyos valores
        propios son exactamente los finitos del problema original.

        - **Pocos grados de libertad con masa** (≤ ``_DENSE_FALLBACK_SIZE``)
          o casi todos los modos: ``K̂`` densa y ``eigh``.
        - **Muchos**: ARPACK sobre el problema condensado sin formar ``K̂``.
          Con ``A = K − σM``, el bloque ``mm`` de ``A⁻¹`` es
          ``(K̂ − σ·Mₘₘ)⁻¹`` (complemento de Schur: ``M`` no tiene filas
          ``0``), así que el operador de desplazamiento e inversión se
          aplica con la factorización dispersa de ``A`` del modelo completo.
        """
        n = K.shape[0]
        mass = np.setdiff1d(np.arange(n), zero)
        n_mass = int(mass.size)

        def fail(cause: str, message: str):
            return SemidefiniteMassError(message, cause=cause, massless=zero, n_mass=n_mass)

        rows = M[zero]
        if rows.nnz and np.any(rows.data != 0.0):
            raise fail("indefinite_mass",
                       "EigenSolver: la matriz de masa tiene filas con diagonal nula y "
                       "términos fuera de la diagonal no nulos: no es semidefinida "
                       "positiva.")
        if n_mass == 0:
            raise fail("no_mass", "EigenSolver: ningún grado de libertad tiene masa.")
        if n_modes > n_mass:
            raise fail("too_many_modes",
                       f"EigenSolver: se piden {n_modes} modos y sólo {n_mass} grados de "
                       f"libertad tienen masa (frecuencias finitas).")

        K_zz = K[zero][:, zero].tocsc()
        K_zm = K[zero][:, mass].tocsc()
        K_mz = K[mass][:, zero].tocsr()
        K_mm = K[mass][:, mass].tocsr()
        M_mm = M[mass][:, mass].tocsr()
        try:
            lu_zz = spla.splu(K_zz)
        except RuntimeError as exc:
            raise fail("massless_mechanism",
                       f"EigenSolver: los grados de libertad sin masa no tienen "
                       f"rigidez que los fije (K₀₀ singular: {exc}).") from exc

        _log.debug(
            f"EigenSolver: {zero.size} grado(s) de libertad sin masa condensados; "
            f"{n_mass} con masa."
        )
        if n_mass <= _DENSE_FALLBACK_SIZE or n_modes >= n_mass - 1:
            K_hat = K_mm.toarray()
            for start in range(0, n_mass, 256):          # sin K₀₀⁻¹K₀ₘ entera
                cols = slice(start, min(start + 256, n_mass))
                K_hat[:, cols] -= K_mz @ lu_zz.solve(K_zm[:, cols].toarray())
            K_hat = 0.5 * (K_hat + K_hat.T)
            try:
                vals, vecs = sla.eigh(K_hat, M_mm.toarray())
            except np.linalg.LinAlgError as exc:
                raise fail("singular_mass",
                           f"EigenSolver: la matriz de masa de los grados de libertad "
                           f"con masa es singular ({exc}).") from exc
            vals, phi_m = vals[:n_modes], vecs[:, :n_modes]
        else:
            vals, phi_m = self._condensed_arpack(K, M, K_mm, K_mz, K_zm, M_mm,
                                                 lu_zz, mass, n_modes, fail)

        phi = np.zeros((n, phi_m.shape[1]))
        phi[mass] = phi_m
        phi[zero] = -lu_zz.solve(np.asarray(K_zm @ phi_m))
        if not np.all(np.isfinite(phi)):
            raise fail("massless_mechanism",
                       "EigenSolver: los grados de libertad sin masa no tienen rigidez "
                       "que los fije (K₀₀ singular a efectos numéricos).")
        return vals, phi

    def _condensed_arpack(self, K, M, K_mm, K_mz, K_zm, M_mm, lu_zz, mass,
                          n_modes, fail):
        """ARPACK sobre ``K̂·φₘ = λ·Mₘₘ·φₘ`` con operadores (sin ``K̂``)."""
        n, n_mass = K.shape[0], mass.size

        def factor(shift):
            A = (K - shift * M) if shift != 0.0 else K
            return spla.splu(sp.csc_matrix(A))

        sigma = self.sigma
        try:
            lu_a = factor(sigma)
        except RuntimeError as exc:
            # Mismo criterio que el camino sin condensar: con σ = 0 y K
            # singular (estructura libre), desplazamiento negativo pequeño.
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
            lu_a = factor(sigma)

        def op_inv(x):
            rhs = np.zeros(n)
            rhs[mass] = np.ravel(x)
            return lu_a.solve(rhs)[mass]

        def k_hat(x):
            x = np.ravel(x)
            return K_mm @ x - K_mz @ lu_zz.solve(np.asarray(K_zm @ x))

        OP = spla.LinearOperator((n_mass, n_mass), matvec=op_inv, dtype=float)
        K_HAT = spla.LinearOperator((n_mass, n_mass), matvec=k_hat, dtype=float)
        try:
            vals, vecs = self._eigsh(K_HAT, M_mm, n_modes, sigma, OPinv=OP)
        except spla.ArpackError as exc:
            raise fail("singular_mass",
                       f"EigenSolver: ARPACK no pudo resolver el problema condensado "
                       f"({exc}); la matriz de masa de los grados de libertad con masa "
                       f"puede ser singular.") from exc
        order = np.argsort(vals)
        return vals[order], vecs[:, order]

    def _eigsh(self, K, M, n_modes, sigma, OPinv=None):
        try:
            return spla.eigsh(
                K, k=n_modes, M=M, sigma=sigma, which=self.which, tol=self.tol,
                OPinv=OPinv,
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
