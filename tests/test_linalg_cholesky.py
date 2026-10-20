"""Pivotes nulos del backend Cholesky (CHOLMOD) — deuda #28, ADR 0019 capa 3.

Hasta el 2026-09-30 ``CholeskyFactorized`` no exponía sus pivotes nulos, así
que la comprobación de pivotes nulos del estático lineal leía 0 con el extra
``cholmod`` instalado: si CHOLMOD factorizaba la matriz singular de un
mecanismo interno con un pivote positivo del orden del redondeo, la solución
pasaba sin aviso. Ahora el backend cuenta los pivotes de ``LDLᵀ`` pequeños
frente a la diagonal de su grado de libertad (``d_i/K_pp < SUSPECT_PIVOT_RTOL``
con la permutación del factor; desde el 2026-10-20, antes frente a
``max|K|``), que es el filtro con el que la red de seguridad busca un
mecanismo interno por la energía del modo.

``scikit-sparse`` no se puede instalar en la máquina Windows de desarrollo,
así que hay dos niveles de prueba:

1. **Simulado** (corre siempre): el módulo del backend se carga contra un
   ``sksparse.cholmod`` de pega que imita la API real —la 0.4 (``Factor``,
   ``factor(b)``, ``factor.D()``, ``factor.P()``) y la 0.5
   (``CholeskyFactor``, ``factorize``/``solve``, ``factor.factor``,
   ``is_ll``, ``rcond``, ``perm``)— con una ``LDLᵀ`` densa sin pivoteo (la
   permutación es la identidad). Imita también lo que motiva la deuda:
   ante el pivote nulo de una matriz semidefinida, el redondeo lo deja
   positivo y la factorización "tiene éxito".
2. **Real** (``skipUnless`` CHOLMOD): la armadura del ejemplo 11 sin diagonal
   con carga horizontal, con ``linear_algebra='cholesky'``, debe lanzar
   ``IllPosedSystemError``. Corre en el job ``test-cholmod`` del CI de Linux.
"""
import importlib.util
import logging
import os
import sys
import types
import unittest
import warnings
from pathlib import Path
from unittest import mock

import numpy as np
import scipy.sparse as sp

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(__file__))

from _armadura_pratt import armadura  # noqa: E402

from solidum.constants import SUSPECT_PIVOT_RTOL  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.linalg import LUSolver  # noqa: E402
from solidum.math.linalg.dispatcher import _HAS_CHOLESKY  # noqa: E402
from solidum.math.solvers import LinearSolver, NonlinearSolver  # noqa: E402
from solidum.math.solvers import corrector as corrector_module  # noqa: E402
from solidum.math.solvers import linear as linear_module  # noqa: E402
from solidum.math.solvers.diagnostics import IllPosedSystemError  # noqa: E402

BACKEND = Path(__file__).resolve().parents[1] / "solidum" / "math" / "linalg" / "cholesky.py"


def setUpModule():
    logging.disable(logging.WARNING)


def tearDownModule():
    logging.disable(logging.NOTSET)


# ----------------------------------------------------------------------
# sksparse.cholmod de pega
# ----------------------------------------------------------------------

class _NotPD(Exception):
    """Hace de ``CholmodNotPositiveDefiniteError``."""


class _CholmodWarning(UserWarning):
    pass


def _ldl_densa(K):
    """``LDLᵀ`` sin pivoteo de una matriz simétrica, como la de CHOLMOD con la
    ordenación natural. Un pivote de magnitud de redondeo (``≤ 1e-12·max|K|``)
    se deja **positivo**: es el caso de la deuda #28, un pivote nulo que
    CHOLMOD encuentra con signo positivo y acepta. Un pivote negativo de
    verdad lanza, como CHOLMOD."""
    A = K.toarray().astype(float)
    n = A.shape[0]
    escala = np.abs(A).max()
    L, d = np.eye(n), np.zeros(n)
    for j in range(n):
        d[j] = A[j, j] - np.sum(L[j, :j] ** 2 * d[:j])
        if abs(d[j]) <= 1e-12 * escala:
            d[j] = abs(d[j]) if d[j] != 0.0 else 1e-17 * escala
        elif d[j] < 0.0:
            raise _NotPD("not positive definite")
        for i in range(j + 1, n):
            L[i, j] = (A[i, j] - np.sum(L[i, :j] * L[j, :j] * d[:j])) / d[j]
    return L, d


def _resolver(L, d, b):
    y = np.linalg.solve(L, b)
    return np.linalg.solve(L.T, y / d)


def _api_04():
    """``sksparse.cholmod`` con la API 0.4: ``cholesky(A)`` → ``Factor``."""
    class Factor:
        def __init__(self, K):
            self.L, self.d = _ldl_densa(K)

        def __call__(self, b):
            return _resolver(self.L, self.d, b)

        def D(self):
            return self.d.copy()

        def P(self):
            return np.arange(self.d.size)

    mod = types.ModuleType("sksparse.cholmod")
    mod.CholmodNotPositiveDefiniteError = _NotPD
    mod.cholesky = lambda A, *a, **k: Factor(A)
    return mod


def _api_05(is_ll=True):
    """``sksparse.cholmod`` con la API 0.5: ``CholeskyFactor``. Con
    ``is_ll`` la vista del factor lleva ``L_ii = √d_i`` en la diagonal; si
    no, ``D_ii`` (forma ``LDLᵀ``)."""
    class CholeskyFactor:
        def __init__(self, A, **kwargs):
            self.A = A

        def factorize(self, A, ldl=None, beta=0.0):
            self.L, self.d = _ldl_densa(A)
            self.is_ll = is_ll
            s = np.sqrt(self.d)
            self.rcond = float((s.min() / s.max()) ** 2)
            diag = s if is_ll else self.d
            self.factor = sp.csc_matrix(np.tril(self.L, -1) + np.diag(diag))
            self.perm = np.arange(self.d.size)
            return self

        def solve(self, b, system="A"):
            if self.rcond < self.L.shape[0] * np.finfo(float).eps:
                warnings.warn(f"Matrix is nearly singular.  Results may be "
                              f"inaccurate (rcond={self.rcond:.2e}).", _CholmodWarning)
            return _resolver(self.L, self.d, b)

    mod = types.ModuleType("sksparse.cholmod")
    mod.CholmodNotPositiveDefiniteError = _NotPD
    mod.CholmodWarning = _CholmodWarning
    mod.CholeskyFactor = CholeskyFactor
    mod.cholesky = lambda A, *a, **k: (None, None)     # la 0.5 devuelve matrices
    return mod


def _backend(cholmod):
    """El módulo ``solidum/math/linalg/cholesky.py`` cargado contra
    ``cholmod`` (no toca el registro real ni ``sys.modules`` después)."""
    paquete = types.ModuleType("sksparse")
    paquete.cholmod = cholmod
    with mock.patch.dict(sys.modules, {"sksparse": paquete, "sksparse.cholmod": cholmod}):
        spec = importlib.util.spec_from_file_location(
            f"_cholesky_simulado_{id(cholmod)}", BACKEND)
        modulo = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(modulo)
    return modulo


APIS = {"0.4": _api_04, "0.5 LLᵀ": lambda: _api_05(True), "0.5 LDLᵀ": lambda: _api_05(False)}


def _K_reducida(sin_diagonal):
    dom, F = armadura(sin_diagonal=sin_diagonal, carga="horizontal")
    a = Assembler(dom)
    a.assemble_system()
    K, F_red, _, _ = a.reduce(a.K_global.tocsr(), F)
    return K, F_red


# ----------------------------------------------------------------------
# 1. Simulado
# ----------------------------------------------------------------------

class TestPivotesNulosSimulados(unittest.TestCase):

    def test_matriz_regular_sin_pivotes_nulos(self):
        K, F = _K_reducida(sin_diagonal=False)
        for api, crear in APIS.items():
            with self.subTest(api=api):
                factor = _backend(crear()).CholeskySolver().factorize(K)
                self.assertEqual(factor.n_small_pivots, 0)
                self.assertEqual(factor.n_negative_pivots, 0)
                np.testing.assert_allclose(factor.solve(F), LUSolver().solve(K, F),
                                           rtol=1e-10, atol=1e-16)

    def test_umbral_del_pivote(self):
        """``[[1, 1], [1, 1 + δ]]`` tiene pivotes ``1`` y ``δ`` y diagonal
        ``(1, 1 + δ)``: pequeño si ``δ/(1 + δ) < SUSPECT_PIVOT_RTOL``, no si
        queda por encima."""
        for api, crear in APIS.items():
            backend = _backend(crear())
            for delta, esperado in ((0.1 * SUSPECT_PIVOT_RTOL, 1), (1e3 * SUSPECT_PIVOT_RTOL, 0)):
                with self.subTest(api=api, delta=delta):
                    K = sp.csc_matrix(np.array([[1.0, 1.0], [1.0, 1.0 + delta]]))
                    self.assertEqual(backend.CholeskySolver().factorize(K).n_small_pivots,
                                     esperado)

    def test_el_filtro_no_depende_de_las_unidades(self):
        """Cambiar las unidades de un grado de libertad escala su fila y su
        columna (``S·K·S``); el cociente pivote/diagonal no cambia. Con el
        criterio anterior (``max|K|``) la matriz escalada daba un pivote
        pequeño más."""
        delta = 1e3 * SUSPECT_PIVOT_RTOL
        K = np.array([[1.0, 1.0, 0.0], [1.0, 1.0 + delta, 0.0], [0.0, 0.0, 1.0]])
        S = np.diag([1.0, 1.0, 1e9])
        for api, crear in APIS.items():
            backend = _backend(crear())
            with self.subTest(api=api):
                for M in (K, S @ K @ S):
                    self.assertEqual(
                        backend.CholeskySolver().factorize(sp.csc_matrix(M)).n_small_pivots, 0)

    def test_mecanismo_interno_da_un_pivote_nulo(self):
        """La matriz de la armadura sin diagonal: CHOLMOD la factorizaría con
        un pivote diminuto; el backend lo cuenta."""
        K, _ = _K_reducida(sin_diagonal=True)
        for api, crear in APIS.items():
            with self.subTest(api=api):
                self.assertEqual(_backend(crear()).CholeskySolver().factorize(K).n_small_pivots, 1)

    def test_no_positiva_definida_se_traduce(self):
        K = sp.csc_matrix(np.array([[1.0, 2.0], [2.0, 1.0]]))
        for api, crear in APIS.items():
            with self.subTest(api=api):
                backend = _backend(crear())
                with self.assertRaises(backend.CholeskyNotPositiveDefiniteError):
                    backend.CholeskySolver().factorize(K)

    def test_api_05_con_rcond_nulo_lanza_al_factorizar(self):
        """La 0.5 sólo detecta ``rcond == 0`` al resolver; el backend lo
        adelanta a la factorización para que el fallback a LU actúe allí."""
        backend = _backend(_api_05())
        original = backend._CholeskyFactor.factorize

        def factorize(self, A, ldl=None, beta=0.0):
            original(self, A, ldl, beta)
            self.rcond = 0.0
            return self

        with mock.patch.object(backend._CholeskyFactor, "factorize", factorize):
            with self.assertRaises(backend.CholeskyNotPositiveDefiniteError):
                backend.CholeskySolver().factorize(sp.csc_matrix(np.eye(3)))

    def test_api_05_silencia_el_aviso_de_rcond(self):
        """El aviso "Matrix is nearly singular" de la 0.5 no llega al usuario:
        la red de seguridad decide con su propio criterio."""
        K, F = _K_reducida(sin_diagonal=True)
        factor = _backend(_api_05()).CholeskySolver().factorize(K)
        with warnings.catch_warnings(record=True) as capturados:
            warnings.simplefilter("always")
            factor.solve(F)
        self.assertEqual([w for w in capturados if "nearly singular" in str(w.message)], [])

    def test_sistema_vacio(self):
        for api, crear in APIS.items():
            with self.subTest(api=api):
                factor = _backend(crear()).CholeskySolver().factorize(sp.csc_matrix((0, 0)))
                self.assertEqual(factor.n_small_pivots, 0)
                self.assertEqual(factor.solve(np.zeros(0)).shape, (0,))


class TestRedDeSeguridadConCholeskySimulado(unittest.TestCase):
    """Los solvers ven el ``n_small_pivots`` del backend Cholesky."""

    def _con_backend(self, modulo, crear_backend):
        return mock.patch.object(modulo, "select_solver",
                                 lambda props, override=None: crear_backend())

    def test_linear_solver_rechaza_el_mecanismo_no_activado(self):
        for api, crear in APIS.items():
            with self.subTest(api=api):
                backend = _backend(crear())
                dom, F = armadura(sin_diagonal=True, carga="horizontal")
                with self._con_backend(linear_module, backend.CholeskySolver):
                    with self.assertRaises(IllPosedSystemError) as cm:
                        LinearSolver(Assembler(dom)).solve(F)
                self.assertEqual(cm.exception.zero_pivots, 1)

    def test_linear_solver_resuelve_la_armadura_bien_planteada(self):
        dom, F = armadura(sin_diagonal=False, carga="vertical")
        U_lu = LinearSolver(Assembler(dom), linear_algebra="lu").solve(F)
        for api, crear in APIS.items():
            with self.subTest(api=api):
                backend = _backend(crear())
                dom, F = armadura(sin_diagonal=False, carga="vertical")
                with self._con_backend(linear_module, backend.CholeskySolver):
                    U = LinearSolver(Assembler(dom)).solve(F)
                np.testing.assert_allclose(U, U_lu, rtol=1e-9, atol=1e-15)

    def test_no_lineal_rechaza_el_mecanismo_en_la_tangente_inicial(self):
        for api, crear in APIS.items():
            with self.subTest(api=api):
                backend = _backend(crear())
                dom, F = armadura(sin_diagonal=True, carga="horizontal")
                with self._con_backend(corrector_module, backend.CholeskySolver):
                    with self.assertRaises(IllPosedSystemError) as cm:
                        NonlinearSolver(Assembler(dom), num_steps=1).solve(F)
                self.assertTrue(cm.exception.initial_tangent)


# ----------------------------------------------------------------------
# 2. Real (CI de Linux, job test-cholmod)
# ----------------------------------------------------------------------

@unittest.skipUnless(_HAS_CHOLESKY, "scikit-sparse no instalado; omitiendo CHOLMOD real")
class TestPivotesNulosConCholmodReal(unittest.TestCase):

    def test_linear_solver_cholesky_rechaza_el_mecanismo_no_activado(self):
        """Ejemplo 11, caso "sin diagonal, carga horizontal". CHOLMOD puede
        encontrar el pivote nulo con signo positivo (lo cuenta
        ``n_small_pivots``) o no positivo (degrada a LU, que lo cuenta): en
        los dos casos, ``IllPosedSystemError``."""
        dom, F = armadura(sin_diagonal=True, carga="horizontal")
        with self.assertRaises(IllPosedSystemError):
            LinearSolver(Assembler(dom), linear_algebra="cholesky").solve(F)

    def test_no_lineal_cholesky_rechaza_el_mecanismo(self):
        dom, F = armadura(sin_diagonal=True, carga="horizontal")
        with self.assertRaises(IllPosedSystemError):
            NonlinearSolver(Assembler(dom), num_steps=1, linear_algebra="cholesky").solve(F)

    def test_armadura_bien_planteada_coincide_con_lu(self):
        dom, F = armadura(sin_diagonal=False, carga="vertical")
        U_lu = LinearSolver(Assembler(dom), linear_algebra="lu").solve(F)
        dom, F = armadura(sin_diagonal=False, carga="vertical")
        U = LinearSolver(Assembler(dom), linear_algebra="cholesky").solve(F)
        np.testing.assert_allclose(U, U_lu, rtol=1e-9, atol=1e-15)

    def test_factor_real_cuenta_pivotes(self):
        from solidum.math.linalg import CholeskySolver  # type: ignore[attr-defined]
        from solidum.math.linalg.cholesky import CholeskyNotPositiveDefiniteError
        K, _ = _K_reducida(sin_diagonal=False)
        self.assertEqual(CholeskySolver().factorize(K).n_small_pivots, 0)
        K, _ = _K_reducida(sin_diagonal=True)
        try:
            n = CholeskySolver().factorize(K).n_small_pivots
        except CholeskyNotPositiveDefiniteError:
            return          # pivote no positivo: CHOLMOD lo rechaza al factorizar
        self.assertGreaterEqual(n, 1)


if __name__ == "__main__":
    unittest.main()
