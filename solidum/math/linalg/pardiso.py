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

Pivotes pequeños (ADR 0019, ampliación del 2026-10-20)
------------------------------------------------------
El handle se inicializa con ``pardisoinit`` (los mismos parámetros que la
MKL usa por omisión) más ``iparm(56) = 1``, que habilita ``pardiso_getdiag``:
tras factorizar se leen los pivotes del factor de la matriz escalada y
emparejada (MC64, ``iparm(11) = iparm(13) = 1`` por omisión), cuya diagonal
es unitaria, así que el pivote ya está relativo a la diagonal. Medido sobre
``Hex8`` con 90 000 DOF: factorizar 1,60 s sin la opción y 1,63 s con ella,
leer la diagonal 0,2 ms, la misma solución a 3e-14. El escalado MC64 no es
exactamente invariante a un cambio de unidades (el pivote mínimo de un
pórtico con zonas rígidas cambia ×10 de m a mm), lo que basta para un
filtro: decide la energía del modo (``model_checks``). Si la MKL no ofrece
``pardiso_getdiag``, el filtro es el número de pivotes perturbados
(``iparm(14)``, ``|p| < 1e-13`` en la matriz escalada).
"""
from __future__ import annotations

import ctypes

import numpy as np
import scipy.sparse as sp

from pypardiso import PyPardisoSolver as _PyPardisoSolver  # type: ignore[import-not-found]

from solidum.constants import SUSPECT_PIVOT_RTOL
from solidum.math.linalg.base import out_of_memory_error


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

    def __init__(self, solver: _PyPardisoSolver, K: sp.csr_matrix, token: int,
                 n_small_pivots: int = 0):
        self._solver = solver
        self._K = K
        self._token = token
        # Pivotes pequeños frente a la diagonal (filtro de la capa 3 del ADR
        # 0019; ver la nota del módulo). Mismo nombre que en LUFactorized.
        self.n_small_pivots = int(n_small_pivots)

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
    Sin incógnitas no hay pivotes, ni pequeños (``n_small_pivots = 0``).
    """

    n_negative_pivots: int | None = 0
    n_small_pivots: int = 0

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
    # ``pardiso_getdiag`` de la MKL, o None si no está disponible.
    _getdiag = None

    @classmethod
    def _handle(cls) -> _PyPardisoSolver:
        if cls._shared_solver is None:
            solver = _PyPardisoSolver()
            cls._getdiag = cls._enable_diagonal(solver)
            cls._shared_solver = solver
        return cls._shared_solver

    @staticmethod
    def _enable_diagonal(solver: _PyPardisoSolver):
        """Parámetros por omisión de la MKL (``pardisoinit``) más
        ``iparm(56) = 1``, que habilita ``pardiso_getdiag``. Devuelve esa
        función, o ``None`` si la MKL no la ofrece (el handle queda entonces
        con ``iparm`` a cero: los valores por omisión, como siempre)."""
        try:
            init = solver.libmkl.pardisoinit
            getdiag = solver.libmkl.pardiso_getdiag
        except AttributeError:
            return None
        init.restype = None
        getdiag.restype = None
        mtype = ctypes.c_int32(solver.mtype)
        init(solver.pt.ctypes.data_as(ctypes.c_void_p), ctypes.byref(mtype),
             solver.iparm.ctypes.data_as(ctypes.c_void_p))
        solver.set_iparm(56, 1)
        return getdiag

    @classmethod
    def _small_pivots(cls, solver: _PyPardisoSolver, n: int) -> int:
        """Pivotes del factor recién calculado por debajo de
        ``SUSPECT_PIVOT_RTOL`` frente a la diagonal. Con el escalado MC64 por
        omisión, el pivote ``df`` ya es el de la matriz escalada, de diagonal
        unitaria; sin escalado se divide por la diagonal ``da``."""
        if cls._getdiag is not None:
            df, da = np.zeros(n), np.zeros(n)
            mnum, err = ctypes.c_int32(1), ctypes.c_int32(0)
            cls._getdiag(solver.pt.ctypes.data_as(ctypes.c_void_p),
                         df.ctypes.data_as(ctypes.c_void_p),
                         da.ctypes.data_as(ctypes.c_void_p),
                         ctypes.byref(mnum), ctypes.byref(err))
            if err.value == 0:
                ratios = np.abs(df)
                if solver.get_iparm(11) == 0:
                    with np.errstate(divide="ignore", invalid="ignore"):
                        ratios = np.where(da != 0.0, ratios / np.abs(da), 0.0)
                return int(np.sum(ratios < SUSPECT_PIVOT_RTOL))
        return int(solver.get_iparm(14))

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
        except MemoryError:
            raise out_of_memory_error(K_csr.shape[0], "Pardiso") from None
        except Exception as exc:  # la MKL señala el fallo con su propio error
            # Código -2 de la MKL: memoria insuficiente.
            if getattr(exc, "value", None) == -2:
                raise out_of_memory_error(K_csr.shape[0], "Pardiso") from exc
            raise RuntimeError(
                f"PardisoSolver: la factorización falló ({exc}). La matriz "
                "puede ser singular."
            ) from exc
        type(self)._token += 1
        return PardisoFactorized(solver, K_csr, type(self)._token,
                                 n_small_pivots=self._small_pivots(solver, K_csr.shape[0]))
