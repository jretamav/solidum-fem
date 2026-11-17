"""Backend directo Intel MKL Pardiso (llamada directa a la MKL; ADR 0017 y 0040).

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

``pypardiso`` sólo localiza la biblioteca ``mkl_rt`` (su atributo público
``libmkl``); la llamada a ``pardiso``, el handle (``pt``), los parámetros
(``iparm``) y las fases son de este módulo (ADR 0040). Si la dependencia no
está instalada, importar el módulo lanza ``ImportError`` y
``solidum.math.linalg.__init__`` lo absorbe, igual que con Cholesky.

Fases de la MKL (ADR 0040)
--------------------------
Factorizar es la fase 12 (análisis + factorización numérica) o, si el
consumidor lo pide con ``factorize(K, reuse_analysis_of=factor_anterior)``,
la fase 22 sola con el análisis de ese factor. El análisis incluye el
emparejamiento y el escalado MC64, que dependen de los valores: reutilizarlo
con valores nuevos da la factorización exacta de la matriz nueva con otra
escala, y la solución cambia en el redondeo (medido: ≤ 4·10⁻¹⁵ relativo en
los desplazamientos de corridas de miles de pasos, sin crecer, como entre dos
corridas con varios hilos). Factorizar pasa de 73 a 6 ms en la viga de
Petersson con Quad8 (25 108 incógnitas). Reglas:

- Se reutiliza sólo si el handle conserva **ese** análisis (nadie analizó
  otra matriz después); si no, fase 12.
- El patrón que cuenta es el de la cadena: una matriz con menos entradas
  (las que hoy valen exactamente cero, p. ej. por la simetría de una malla
  estructurada) se escribe sobre él con ceros explícitos y reutiliza el
  análisis; una con entradas nuevas se analiza sobre la unión de los
  patrones, que sirve para las siguientes. Sin la cadena, la matriz se
  factoriza tal cual: los resultados son los de siempre, bit a bit.
- Si la fase 22 perturba algún pivote (``iparm(14) > 0``) o falla por otra
  causa que la memoria, se repite con análisis propio.
- Los pivotes que lee la red de seguridad (``n_small_pivots``) se calculan
  siempre sobre una factorización con análisis propio: en un factor con
  análisis reutilizado se calculan al pedirlos, refactorizando la misma
  matriz con análisis propio y volviendo a ligar el mismo objeto.

Resolver es la fase 33, con los índices en base 1 convertidos una vez por
factorización. Ante cualquier error de la MKL el handle olvida su análisis y
ninguna factorización anterior sigue vigente.

Comprobaciones previas
----------------------
Las mismas que hacía ``pypardiso``: matriz cuadrada, CSR de ``float64`` con
índices ordenados y ninguna fila vacía. La MKL no rechaza una fila vacía:
perturba el pivote y devuelve una solución finita y falsa; aquí es un
``RuntimeError`` (matriz singular), lo que los solvers toman por tangente
singular y la red de seguridad por un modelo mal planteado.

Notas sobre el tipo de matriz
-----------------------------
Se usa ``mtype = 11`` (real no simétrica), que recibe la matriz **completa**
en CSR. No se usa ``mtype = -2`` (real simétrica indefinida) aunque las
matrices de Solidum suelen ser simétricas: ese tipo espera **sólo el
triángulo superior**, y pasarle la matriz completa produce una violación de
acceso en la MKL, no una excepción de Python. Aprovechar la simetría —y con
ella el conteo de inercia de la deuda #7— exige construir el triángulo
superior y validarlo aparte; queda fuera del alcance de este backend.

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

# Sólo para localizar ``mkl_rt`` (ver la nota del módulo).
from pypardiso import PyPardisoSolver as _PyPardisoSolver  # type: ignore[import-not-found]

from solidum.constants import SUSPECT_PIVOT_RTOL
from solidum.math.linalg.base import out_of_memory_error

# Real no simétrica, matriz completa en CSR (ver la nota del módulo).
_MTYPE = 11
_PT_C, _PT_NP = ((ctypes.c_int64, np.int64) if ctypes.sizeof(ctypes.c_void_p) == 8
                 else (ctypes.c_int32, np.int32))
_I32P = ctypes.POINTER(ctypes.c_int32)
# Código de error de la MKL: memoria insuficiente.
_MKL_OUT_OF_MEMORY = -2


def _pattern_keys(M: sp.csr_matrix) -> np.ndarray:
    """Clave ``fila·n + columna`` de cada entrada de una CSR canónica, en
    orden ascendente."""
    n = M.shape[1]
    rows = np.repeat(np.arange(M.shape[0], dtype=np.int64), np.diff(M.indptr))
    return rows * n + M.indices.astype(np.int64)


def _on_chain_pattern(K: sp.csr_matrix, P: sp.csr_matrix) -> tuple[sp.csr_matrix, bool]:
    """``K`` llevada al patrón ``P`` de la cadena de factorizaciones (ADR 0040).

    El patrón de la matriz reducida depende de los valores: la reducción y
    las sumas de SciPy descartan las entradas que valen exactamente cero, y
    en una malla estructurada muchas se anulan por simetría hasta que la
    plasticidad o el daño la rompen. Para que eso no obligue a repetir el
    análisis simbólico, la matriz se escribe sobre el patrón de la cadena
    con ceros explícitos donde le falten entradas: la matriz es la misma.

    Devuelve ``(matriz, creció)``: con ``creció`` falso la matriz tiene el
    patrón ``P`` (``K`` misma si ya lo tenía) y el análisis de la cadena
    sirve; si ``K`` tiene entradas fuera de ``P``, la matriz va sobre la unión
    de los dos patrones, que necesita un análisis propio y sirve para las
    siguientes de la cadena."""
    if (K.nnz == P.nnz and np.array_equal(K.indptr, P.indptr)
            and np.array_equal(K.indices, P.indices)):
        return K, False
    kK, kP = _pattern_keys(K), _pattern_keys(P)
    pos = np.searchsorted(kP, kK)
    dentro = pos < kP.size
    dentro[dentro] = kP[pos[dentro]] == kK[dentro]
    if dentro.all():
        data = np.zeros(kP.size)
        data[pos] = K.data
        return sp.csr_matrix((data, P.indices, P.indptr), shape=K.shape), False
    n = K.shape[1]
    keys = np.union1d(kP, kK)
    indptr = np.zeros(K.shape[0] + 1, dtype=np.int64)
    np.cumsum(np.bincount(keys // n, minlength=K.shape[0]), out=indptr[1:])
    data = np.zeros(keys.size)
    data[np.searchsorted(keys, kK)] = K.data
    union = sp.csr_matrix((data, (keys % n).astype(np.int32), indptr.astype(np.int32)),
                          shape=K.shape)
    return union, True


class PardisoError(RuntimeError):
    """Código de error distinto de cero devuelto por la MKL en una fase."""

    def __init__(self, code: int, phase: int):
        self.value = int(code)
        self.phase = int(phase)
        super().__init__(f"la MKL Pardiso devolvió el error {self.value} en la fase {self.phase}")


class _MklPardiso:
    """La función ``pardiso`` de la MKL con su handle (``pt``) y sus
    parámetros (``iparm``). Uno por proceso: crear uno localiza la
    biblioteca, lo que cuesta ~2,5 s (``pypardiso`` recorre el disco)."""

    def __init__(self) -> None:
        lib = _PyPardisoSolver().libmkl
        try:
            fn = lib.pardiso
        except AttributeError as exc:  # pragma: no cover - depende de la MKL instalada
            raise RuntimeError(
                "PardisoSolver: la biblioteca MKL que localizó pypardiso no exporta "
                "'pardiso'. Use linear_algebra: lu o reinstale pypardiso/mkl."
            ) from exc
        fn.argtypes = [ctypes.POINTER(_PT_C), _I32P, _I32P, _I32P, _I32P, _I32P,
                       ctypes.c_void_p, _I32P, _I32P, _I32P, _I32P, _I32P, _I32P,
                       ctypes.c_void_p, ctypes.c_void_p, _I32P]
        fn.restype = None
        self._fn = fn
        self.libmkl = lib
        self.mtype = _MTYPE
        self.pt = np.zeros(64, dtype=_PT_NP)
        self.iparm = np.zeros(64, dtype=np.int32)
        self._perm = np.zeros(0, dtype=np.int32)
        self.getdiag = self._init_parameters()
        self.set_iparm(12, 0)       # A·x = b con la matriz en CSR (no la transpuesta)

    def _init_parameters(self):
        """Parámetros por omisión de la MKL (``pardisoinit``) más
        ``iparm(56) = 1``, que habilita ``pardiso_getdiag``. Devuelve esa
        función, o ``None`` si la MKL no la ofrece (``iparm`` queda entonces a
        cero: los valores por omisión)."""
        try:
            init = self.libmkl.pardisoinit
            getdiag = self.libmkl.pardiso_getdiag
        except AttributeError:
            return None
        init.restype = None
        getdiag.restype = None
        mtype = ctypes.c_int32(self.mtype)
        init(self.pt.ctypes.data_as(ctypes.c_void_p), ctypes.byref(mtype),
             self.iparm.ctypes.data_as(ctypes.c_void_p))
        self.set_iparm(56, 1)
        return getdiag

    def get_iparm(self, i: int) -> int:
        """``iparm(i)`` con la numeración de la documentación (desde 1)."""
        return int(self.iparm[i - 1])

    def set_iparm(self, i: int, value: int) -> None:
        self.iparm[i - 1] = value

    def call(self, phase: int, data: np.ndarray, ia: np.ndarray, ja: np.ndarray,
             b: np.ndarray | None = None) -> np.ndarray | None:
        """Ejecuta ``phase`` sobre la matriz CSR ``(data, ia, ja)`` en base 1.
        Con ``b`` (fase 33) devuelve la solución; lanza :class:`PardisoError`
        si la MKL devuelve un error."""
        n = ia.shape[0] - 1
        if b is None:
            # Las fases 11, 12 y 22 no leen b ni escriben x.
            b = np.zeros(n)
            nrhs = 1
        else:
            b = np.asfortranarray(b, dtype=np.float64)
            nrhs = 1 if b.ndim == 1 else b.shape[1]
        x = np.zeros_like(b)
        error = ctypes.c_int32(0)
        self._fn(self.pt.ctypes.data_as(ctypes.POINTER(_PT_C)),
                 ctypes.byref(ctypes.c_int32(1)),            # maxfct
                 ctypes.byref(ctypes.c_int32(1)),            # mnum
                 ctypes.byref(ctypes.c_int32(self.mtype)),
                 ctypes.byref(ctypes.c_int32(phase)),
                 ctypes.byref(ctypes.c_int32(n)),
                 data.ctypes.data_as(ctypes.c_void_p),
                 ia.ctypes.data_as(_I32P), ja.ctypes.data_as(_I32P),
                 self._perm.ctypes.data_as(_I32P),
                 ctypes.byref(ctypes.c_int32(nrhs)),
                 self.iparm.ctypes.data_as(_I32P),
                 ctypes.byref(ctypes.c_int32(0)),            # msglvl
                 b.ctypes.data_as(ctypes.c_void_p), x.ctypes.data_as(ctypes.c_void_p),
                 ctypes.byref(error))
        if error.value != 0:
            raise PardisoError(error.value, phase)
        return np.ascontiguousarray(x) if phase == 33 else None


class PardisoFactorized:
    """Factorización numérica de Pardiso retenida y reutilizable.

    La factorización vive dentro del handle de la MKL; cada ``solve(b)``
    ejecuta únicamente las sustituciones triangulares (fase 33). Medido sobre
    ``Hex8 20³``: factorizar cuesta 0,72 s y cada ``solve`` posterior
    0,07 s — que es lo que hace rentable el Newton modificado y el reuso
    entre pasos.

    La MKL no reporta inercia con el tipo de matriz que usa este backend,
    así que ``n_negative_pivots`` es ``None`` (ver nota del módulo).
    """

    n_negative_pivots: int | None = None

    def __init__(self, K: sp.csr_matrix, ia: np.ndarray, ja: np.ndarray, token: int,
                 analysis_id: int, n_small_pivots: int | None, n_perturbed_pivots: int = 0):
        self._K = K
        self._ia = ia
        self._ja = ja
        self._token = token
        self._analysis_id = analysis_id
        # ``None``: análisis reutilizado; se calcula al pedirlo (ver la nota
        # del módulo).
        self._n_small_pivots = n_small_pivots
        # Pivotes que la MKL perturbó al factorizar (``iparm(14)``): con
        # alguno la matriz es singular o casi, y quien encadena
        # factorizaciones deja de reutilizar el análisis (ADR 0040).
        self.n_perturbed_pivots = int(n_perturbed_pivots)

    @property
    def is_current(self) -> bool:
        """``True`` mientras el handle conserve esta factorización (ADR 0040)."""
        return self._token == PardisoSolver._token

    @property
    def shares_analysis(self) -> bool:
        """``True`` si se factorizó con el análisis de otra matriz de su cadena."""
        return self._n_small_pivots is None

    @property
    def n_small_pivots(self) -> int:
        """Pivotes pequeños frente a la diagonal (filtro de la capa 3 del ADR
        0019; ver la nota del módulo). Mismo nombre que en ``LUFactorized``.
        Siempre de una factorización con análisis propio."""
        if self._n_small_pivots is None:
            PardisoSolver._refactorize_with_own_analysis(self)
        return int(self._n_small_pivots)

    def solve(self, b: np.ndarray) -> np.ndarray:
        # El handle de la MKL es único por proceso y guarda **una sola**
        # factorización. Si otra la ha reemplazado, esta ya no es válida:
        # se falla de forma ruidosa en vez de devolver la solución del otro
        # sistema, que sería un error silencioso. Ocurre sólo si dos
        # factorizaciones se mantienen vivas a la vez (p. ej. M y A en un
        # solver dinámico); el patrón normal —factorizar, resolver N
        # veces— no se ve afectado. Quien retiene factores consulta antes
        # ``is_current`` (ADR 0040).
        if not self.is_current:
            raise RuntimeError(
                "PardisoFactorized: esta factorización fue invalidada por una "
                "posterior (el handle de MKL Pardiso guarda una sola). Refactoriza "
                "antes de resolver, o usa el backend 'lu' si necesitas mantener "
                "dos factorizaciones simultáneas."
            )
        try:
            x = PardisoSolver._handle().call(33, self._K.data, self._ia, self._ja, b)
        except PardisoError as exc:
            PardisoSolver._forget_handle_state()
            raise RuntimeError(f"PardisoSolver: la resolución falló ({exc}).") from exc
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
    is_current: bool = True

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
    # Puede reutilizar el análisis simbólico de un factor anterior de la
    # misma cadena (``factorize(K, reuse_analysis_of=...)``, ADR 0040).
    reuses_analysis = True

    # Handle único por proceso (localizar la biblioteca cuesta ~2,5 s; ver
    # ``_MklPardiso``). Es reutilizable: cada factorización sustituye a la
    # anterior.
    _shared_solver: _MklPardiso | None = None
    # Identifica la factorización viva en el handle; ver PardisoFactorized.solve.
    _token: int = 0
    # Identificador del análisis que conserva el handle (0: ninguno) y
    # contador para numerarlos (ADR 0040).
    _analysis_id: int = 0
    _analysis_seq: int = 0
    # ``pardiso_getdiag`` de la MKL, o None si no está disponible.
    _getdiag = None

    @classmethod
    def _handle(cls) -> _MklPardiso:
        if cls._shared_solver is None:
            handle = _MklPardiso()
            cls._getdiag = handle.getdiag
            cls._shared_solver = handle
        return cls._shared_solver

    @classmethod
    def _forget_handle_state(cls) -> None:
        """Tras un error de la MKL el handle queda en un estado indefinido:
        ninguna factorización anterior sigue vigente y su análisis no se
        reutiliza (ADR 0040)."""
        cls._token += 1
        cls._analysis_id = 0

    @classmethod
    def _new_analysis_id(cls) -> int:
        cls._analysis_seq += 1
        cls._analysis_id = cls._analysis_seq
        return cls._analysis_id

    @classmethod
    def _small_pivots(cls, handle: _MklPardiso, n: int) -> int:
        """Pivotes del factor recién calculado por debajo de
        ``SUSPECT_PIVOT_RTOL`` frente a la diagonal. Con el escalado MC64 por
        omisión, el pivote ``df`` ya es el de la matriz escalada, de diagonal
        unitaria; sin escalado se divide por la diagonal ``da``."""
        if cls._getdiag is not None:
            df, da = np.zeros(n), np.zeros(n)
            mnum, err = ctypes.c_int32(1), ctypes.c_int32(0)
            cls._getdiag(handle.pt.ctypes.data_as(ctypes.c_void_p),
                         df.ctypes.data_as(ctypes.c_void_p),
                         da.ctypes.data_as(ctypes.c_void_p),
                         ctypes.byref(mnum), ctypes.byref(err))
            if err.value == 0:
                ratios = np.abs(df)
                if handle.get_iparm(11) == 0:
                    with np.errstate(divide="ignore", invalid="ignore"):
                        ratios = np.where(da != 0.0, ratios / np.abs(da), 0.0)
                return int(np.sum(ratios < SUSPECT_PIVOT_RTOL))
        return handle.get_iparm(14)

    @staticmethod
    def _as_csr(K: sp.spmatrix) -> sp.csr_matrix:
        """CSR de ``float64`` con índices ordenados y sin filas vacías.

        Una matriz CSC o de otro tipo se convierte aquí y no en el punto de
        llamada. Ordenar los índices es imprescindible (la MKL da resultados
        erróneos sin orden) y se hace en sitio, sin cambiar la matriz. Una
        fila vacía hace la matriz singular, y la MKL no lo detecta (ver la
        nota del módulo)."""
        if K.shape[0] != K.shape[1]:
            raise RuntimeError(
                f"PardisoSolver: la matriz debe ser cuadrada (forma {K.shape}).")
        K_csr = K.tocsr() if not sp.isspmatrix_csr(K) else K
        if K_csr.dtype != np.float64:
            K_csr = K_csr.astype(np.float64)
        if not K_csr.has_sorted_indices:
            K_csr.sort_indices()
        if not np.diff(K_csr.indptr).all():
            raise RuntimeError(
                "PardisoSolver: la factorización falló (la matriz tiene filas "
                "vacías). La matriz es singular.")
        return K_csr

    def solve(self, K: sp.spmatrix, b: np.ndarray) -> np.ndarray:
        return self.factorize(K).solve(b)

    def factorize(self, K: sp.spmatrix, *,
                  reuse_analysis_of=None) -> PardisoFactorized | _EmptyFactorized:
        """Factoriza ``K``. Con ``reuse_analysis_of`` (factor anterior de la
        misma cadena) reutiliza su análisis si sigue en el handle y el patrón
        es idéntico (ver la nota del módulo)."""
        K_csr = self._as_csr(K)
        n = K_csr.shape[0]
        # Sistema vacío (0×0): ocurre cuando todos los DOF están prescritos
        # —un modelo de un elemento controlado por desplazamiento deja el
        # sistema reducido sin incógnitas—. SuperLU lo acepta y devuelve un
        # vector vacío; la MKL lo rechaza como singular. Se atiende aquí para
        # que el backend sea sustituto equivalente de LU también en el borde.
        if n == 0:
            return _EmptyFactorized()

        cls = type(self)
        handle = self._handle()
        prev = reuse_analysis_of
        reuse = (isinstance(prev, PardisoFactorized)
                 and cls._analysis_id != 0
                 and prev._analysis_id == cls._analysis_id
                 and prev._K.shape == K_csr.shape)
        if reuse:
            # Patrón de la cadena: una entrada que hoy vale exactamente cero
            # (simetrías de una malla estructurada que la plasticidad rompe)
            # no debe obligar a otro análisis (ver ``_on_chain_pattern``).
            K_csr, grew = _on_chain_pattern(K_csr, prev._K)
            reuse = not grew
        ia = K_csr.indptr.astype(np.int32) + 1
        ja = K_csr.indices.astype(np.int32) + 1
        try:
            if reuse:
                try:
                    handle.call(22, K_csr.data, ia, ja)
                    # Pivotes perturbados con un análisis ajeno: se repite con
                    # análisis propio.
                    reuse = handle.get_iparm(14) == 0
                except PardisoError as exc:
                    if exc.value == _MKL_OUT_OF_MEMORY:
                        raise
                    # Que una matriz se declare singular no debe depender de
                    # que se reutilizó el análisis: se repite con uno propio.
                    reuse = False
            if not reuse:
                handle.call(12, K_csr.data, ia, ja)
                cls._new_analysis_id()
        except MemoryError:
            cls._forget_handle_state()
            raise out_of_memory_error(n, "Pardiso") from None
        except PardisoError as exc:
            cls._forget_handle_state()
            if exc.value == _MKL_OUT_OF_MEMORY:
                raise out_of_memory_error(n, "Pardiso") from exc
            raise RuntimeError(
                f"PardisoSolver: la factorización falló ({exc}). La matriz "
                "puede ser singular."
            ) from exc
        cls._token += 1
        n_small = None if reuse else cls._small_pivots(handle, n)
        return PardisoFactorized(K_csr, ia, ja, cls._token, cls._analysis_id, n_small,
                                 n_perturbed_pivots=handle.get_iparm(14))

    @classmethod
    def _refactorize_with_own_analysis(cls, factor: PardisoFactorized) -> None:
        """Factoriza otra vez la matriz de ``factor`` con análisis propio y
        vuelve a ligar el mismo objeto (token y análisis nuevos), para leer
        sus pivotes. Nunca invalida ``factor``: invalida, como toda
        factorización, a los demás."""
        handle = cls._handle()
        n = factor._K.shape[0]
        try:
            handle.call(12, factor._K.data, factor._ia, factor._ja)
        except PardisoError as exc:
            cls._forget_handle_state()
            if exc.value == _MKL_OUT_OF_MEMORY:
                raise out_of_memory_error(n, "Pardiso") from exc
            raise RuntimeError(f"PardisoSolver: la factorización falló ({exc}).") from exc
        cls._token += 1
        factor._token = cls._token
        factor._analysis_id = cls._new_analysis_id()
        factor._n_small_pivots = cls._small_pivots(handle, n)
