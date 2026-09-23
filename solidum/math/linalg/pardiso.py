"""Backend directo Intel MKL Pardiso (vía ``pypardiso``).

Solver directo disperso **multihilo** con reordenamiento por disección
anidada (nested dissection, METIS). Resuelve el mismo sistema que
:class:`~solidum.math.linalg.lu.LUSolver` —factorización directa, sin
iteración— pero con dos ventajas que crecen con el tamaño del problema:

- **Paraleliza** la factorización sobre los núcleos disponibles; SuperLU
  (el backend de SciPy) es monohilo.
- **Reduce el relleno** (*fill-in*) con disección anidada. El relleno es lo
  que hace que el coste de un solver directo escale mucho peor que lineal:
  medido sobre ``Hex8 n³`` con SuperLU/COLAMD, ``nnz(L+U)/nnz(K)`` pasa de
  ×10 a 4 000 grados de libertad a ×25 a 28 000. El relleno domina el
  tiempo **y** la memoria.

Medición que motivó el backend (ADR 0017, Windows, 16 hilos, misma
solución a ``1e-13`` relativo):

===============  ========  ===========  ==========  =========
Malla            DOF       SuperLU      Pardiso     Speedup
===============  ========  ===========  ==========  =========
``Hex8 10³``        3 993     0,63 s      0,08 s      ×7,5
``Hex8 15³``       12 288     3,06 s      0,27 s      ×11
``Hex8 20³``       27 783    19,24 s      0,72 s      ×27
``Hex8 25³``       52 728   116,23 s      1,94 s      ×60
===============  ========  ===========  ==========  =========

El import de ``pypardiso`` es responsabilidad de este módulo: si la
dependencia no está instalada lanza ``ImportError`` al importarse y
``solidum.math.linalg.__init__`` lo absorbe, igual que con Cholesky.

Notas sobre el tipo de matriz
-----------------------------
Se usa el tipo por defecto de ``pypardiso`` (``mtype = 11``, real no
simétrica), que recibe la matriz **completa** en CSR. No se usa
``mtype = -2`` (real simétrica indefinida) aunque las matrices de Solidum
suelen ser simétricas: ese tipo espera **sólo el triángulo superior**, y
pasarle la matriz completa produce una violación de acceso en la MKL, no
una excepción de Python. Aprovechar la simetría —y con ella el conteo de
inercia de la deuda #7— exige construir el triángulo superior y validarlo
aparte; queda fuera del alcance de este backend.

Por tanto :attr:`PardisoFactorized.n_negative_pivots` es ``None``: con
``mtype = 11`` la MKL no reporta inercia (``iparm[22]`` queda en cero para
cualquier matriz, así que leerlo sería *peor* que no exponerlo).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from pypardiso import PyPardisoSolver as _PyPardisoSolver  # type: ignore[import-not-found]


class PardisoFactorized:
    """Factorización numérica de Pardiso retenida y reutilizable.

    La factorización vive dentro del handle de la MKL; cada ``solve(b)``
    ejecuta únicamente las sustituciones triangulares. Medido sobre
    ``Hex8 20³``: factorizar cuesta 0,72 s y cada ``solve`` posterior
    0,07 s — que es lo que hace rentable el Newton modificado y el reuso
    entre pasos.

    La MKL no reporta inercia con el tipo de matriz que usa este backend,
    así que ``n_negative_pivots`` es ``None`` (ver nota del módulo).
    """

    n_negative_pivots: int | None = None

    def __init__(self, solver: _PyPardisoSolver, K: sp.csr_matrix, token: int):
        self._solver = solver
        self._K = K
        self._token = token

    def solve(self, b: np.ndarray) -> np.ndarray:
        # El handle de la MKL es único por proceso y guarda **una sola**
        # factorización. Si otra la ha reemplazado, esta ya no es válida:
        # se falla de forma ruidosa en vez de devolver la solución del otro
        # sistema, que sería un error silencioso. Ocurre sólo si dos
        # factorizaciones se mantienen vivas a la vez (p. ej. M y A en un
        # solver dinámico); el patrón normal —factorizar, resolver N
        # veces— no se ve afectado.
        if self._token != PardisoSolver._token:
            raise RuntimeError(
                "PardisoFactorized: esta factorización fue invalidada por una "
                "posterior (el handle de MKL Pardiso guarda una sola). Refactoriza "
                "antes de resolver, o usa el backend 'lu' si necesitas mantener "
                "dos factorizaciones simultáneas."
            )
        x = self._solver.solve(self._K, np.ascontiguousarray(b, dtype=np.float64))
        if not np.all(np.isfinite(x)):
            raise RuntimeError(
                "PardisoSolver: la solución no es finita (matriz singular o "
                "casi singular)."
            )
        return x


class _EmptyFactorized:
    """Factorización de un sistema ``0×0``: ``solve`` devuelve el vector vacío.

    No toca la MKL, así que tampoco invalida la factorización vigente.
    """

    n_negative_pivots: int | None = 0

    def solve(self, b: np.ndarray) -> np.ndarray:
        return np.zeros(0)


class PardisoSolver:
    """Solver directo disperso multihilo (Intel MKL Pardiso).

    Aplicable a cualquier matriz real no singular: no asume simetría ni
    positividad, igual que :class:`LUSolver`, del que es sustituto directo.

    El número de hilos lo fija la MKL (``MKL_NUM_THREADS`` /
    ``OMP_NUM_THREADS``); el backend no lo impone para no pelearse con la
    configuración de hilos de Numba del ensamblaje por lotes (ADR 0014).
    """

    name = "pardiso"

    # Handle único por proceso. Construir un ``PyPardisoSolver`` localiza la
    # DLL de la MKL recorriendo el sistema de archivos con ``glob``: medido,
    # 2,5 s por instancia. Con una instancia por factorización, un Newton de
    # diez iteraciones gastaba 25 s buscando la biblioteca frente a 7 s
    # factorizando. El handle es reutilizable —``factorize(A)`` sustituye la
    # factorización previa— así que se crea una sola vez y se comparte.
    _shared_solver: _PyPardisoSolver | None = None
    # Identifica la factorización viva en el handle; ver PardisoFactorized.solve.
    _token: int = 0

    @classmethod
    def _handle(cls) -> _PyPardisoSolver:
        if cls._shared_solver is None:
            cls._shared_solver = _PyPardisoSolver()
        return cls._shared_solver

    @staticmethod
    def _as_csr(K: sp.spmatrix) -> sp.csr_matrix:
        """CSR con índices y datos en los tipos que la MKL espera.

        ``pypardiso`` exige CSR de ``float64``; un CSC o una matriz de otro
        dtype se convierte aquí y no en el punto de llamada.
        """
        K_csr = K.tocsr() if not sp.isspmatrix_csr(K) else K
        if K_csr.dtype != np.float64:
            K_csr = K_csr.astype(np.float64)
        return K_csr

    def solve(self, K: sp.spmatrix, b: np.ndarray) -> np.ndarray:
        return self.factorize(K).solve(b)

    def factorize(self, K: sp.spmatrix) -> PardisoFactorized:
        K_csr = self._as_csr(K)

        # Sistema vacío (0×0): ocurre cuando todos los DOF están prescritos
        # —un modelo de un elemento controlado por desplazamiento deja el
        # sistema reducido sin incógnitas—. SuperLU lo acepta y devuelve un
        # vector vacío; la MKL lo rechaza como singular. Se atiende aquí para
        # que el backend sea sustituto equivalente de LU también en el borde.
        if K_csr.shape[0] == 0:
            return _EmptyFactorized()

        solver = self._handle()
        try:
            solver.factorize(K_csr)
        except Exception as exc:  # la MKL señala el fallo con su propio error
            raise RuntimeError(
                f"PardisoSolver: la factorización falló ({exc}). La matriz "
                "puede ser singular."
            ) from exc
        type(self)._token += 1
        return PardisoFactorized(solver, K_csr, type(self)._token)
