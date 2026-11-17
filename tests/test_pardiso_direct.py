"""Pardiso por llamada directa a la MKL y análisis reutilizado en una cadena (ADR 0040, punto 1).

1. Las comprobaciones previas que hacía ``pypardiso`` siguen: una fila vacía
   es ``RuntimeError`` (la MKL sola devolvería una solución finita y falsa);
   una matriz no cuadrada también; índices desordenados y entrada CSC dan la
   solución correcta.
2. Ante un error de la MKL el handle olvida su análisis y ningún factor
   anterior sigue vigente.
3. Fases: sin ``reuse_analysis_of``, análisis propio (fase 12); con el factor
   anterior de la cadena, sólo la numérica (fase 22) si el handle conserva
   ese análisis y el patrón es idéntico; si otra matriz se analizó en medio o
   el patrón cambia, análisis propio. Si la fase 22 perturba pivotes, se
   repite con análisis propio.
4. La solución con análisis reutilizado es exacta al redondeo, y los pivotes
   de la red de seguridad son los de un análisis propio: el factor sigue
   vigente y resolviendo después de leerlos.
5. Con la MKL en un hilo, la factorización y la resolución (fase 12 + 33) son
   idénticas bit a bit a las de ``pypardiso`` con los mismos parámetros, que
   es como Solidum llamaba a la MKL antes del ADR 0040.
"""
import ctypes
import os
import unittest
from unittest import mock

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from solidum.math.linalg.dispatcher import _HAS_PARDISO

if _HAS_PARDISO:
    from solidum.math.linalg.pardiso import PardisoError, PardisoFactorized, PardisoSolver


def _matriz(n=60, seed=0, escala=1.0):
    """No simétrica y diagonalmente dominante, con el mismo patrón para
    cualquier ``escala`` de los valores fuera de la diagonal."""
    M = sp.random(n, n, density=0.1, random_state=seed, format="csr")
    M.data[:] = np.linspace(0.5, 1.5, M.nnz)
    return (escala * M + (2.0 * n) * sp.eye(n)).tocsr()


class _Fases:
    """Registra las fases que se ejecutan en el handle."""

    def __enter__(self):
        handle = PardisoSolver._handle()
        original, self.fases = handle.call, []

        def call(phase, *a, **k):
            self.fases.append(phase)
            return original(phase, *a, **k)
        self._patch = mock.patch.object(handle, "call", side_effect=call)
        self._patch.start()
        return self

    def __exit__(self, *exc):
        self._patch.stop()


@unittest.skipUnless(_HAS_PARDISO, "pypardiso no instalado")
class TestComprobacionesPrevias(unittest.TestCase):

    def test_una_fila_vacia_es_matriz_singular(self):
        A = _matriz().tolil()
        A[7, :] = 0.0
        A = sp.csr_matrix(A)
        A.eliminate_zeros()
        with self.assertRaisesRegex(RuntimeError, "filas vacías"):
            PardisoSolver().factorize(A)

    def test_una_matriz_no_cuadrada_es_error(self):
        with self.assertRaises(RuntimeError):
            PardisoSolver().factorize(sp.random(5, 6, density=0.5, format="csr"))

    def test_indices_desordenados_y_csc_dan_la_solucion_correcta(self):
        A = _matriz()
        b = np.arange(1.0, 61.0)
        x_ref = spla.spsolve(A.tocsc(), b)
        filas = np.repeat(np.arange(60), np.diff(A.indptr))
        orden = np.lexsort((-A.indices, filas))              # columnas descendentes
        desordenada = sp.csr_matrix((A.data[orden], A.indices[orden], A.indptr), shape=A.shape)
        desordenada.has_sorted_indices = False
        for nombre, M in (("desordenada", desordenada), ("csc", A.tocsc())):
            with self.subTest(entrada=nombre):
                x = PardisoSolver().factorize(M).solve(b)
                np.testing.assert_allclose(x, x_ref, rtol=1e-12)


@unittest.skipUnless(_HAS_PARDISO, "pypardiso no instalado")
class TestErroresDeLaMKL(unittest.TestCase):

    def test_tras_un_error_ningun_factor_anterior_sigue_vigente(self):
        f1 = PardisoSolver().factorize(_matriz())
        self.assertTrue(f1.is_current)
        handle = PardisoSolver._handle()
        with mock.patch.object(handle, "call", side_effect=PardisoError(-4, 12)):
            with self.assertRaises(RuntimeError):
                PardisoSolver().factorize(_matriz(seed=1))
        self.assertFalse(f1.is_current)
        self.assertEqual(PardisoSolver._analysis_id, 0)
        with self.assertRaisesRegex(RuntimeError, "invalidada"):
            f1.solve(np.ones(60))


@unittest.skipUnless(_HAS_PARDISO, "pypardiso no instalado")
class TestAnalisisEnUnaCadena(unittest.TestCase):

    def test_sin_cadena_analisis_propio_con_cadena_solo_la_numerica(self):
        A1, A2, A3 = _matriz(escala=1.0), _matriz(escala=3.0), _matriz(escala=5.0)
        with _Fases() as f:
            f1 = PardisoSolver().factorize(A1)
            f2 = PardisoSolver().factorize(A2, reuse_analysis_of=f1)
            f3 = PardisoSolver().factorize(A3)                     # cadena nueva
        self.assertEqual(f.fases, [12, 22, 12])
        self.assertFalse(f1.shares_analysis)
        self.assertTrue(f2.shares_analysis)
        self.assertFalse(f3.shares_analysis)

    def test_si_otra_matriz_se_analizo_en_medio_analisis_propio(self):
        A1, A2 = _matriz(escala=1.0), _matriz(escala=2.0)
        f1 = PardisoSolver().factorize(A1)
        PardisoSolver().factorize(_matriz(seed=9))                # otro consumidor
        with _Fases() as f:
            PardisoSolver().factorize(A2, reuse_analysis_of=f1)
        self.assertEqual(f.fases, [12])

    def test_con_otro_patron_analisis_propio(self):
        f1 = PardisoSolver().factorize(_matriz(seed=0))
        with _Fases() as f:
            PardisoSolver().factorize(_matriz(seed=4), reuse_analysis_of=f1)
        self.assertEqual(f.fases, [12])

    def test_con_menos_entradas_se_completa_el_patron_y_se_reutiliza(self):
        """Una entrada que pasa a valer exactamente cero (y la reducción
        descarta) no obliga a otro análisis: la matriz se escribe sobre el
        patrón de la cadena con un cero explícito, y la solución es la suya."""
        A1 = _matriz()
        A2 = A1.copy()
        fuera_de_la_diagonal = A2.indices != np.repeat(np.arange(60), np.diff(A2.indptr))
        A2.data[np.flatnonzero(fuera_de_la_diagonal)[0]] = 0.0
        A2.eliminate_zeros()
        self.assertEqual(A2.nnz, A1.nnz - 1)
        b = np.linspace(1.0, 2.0, 60)
        f1 = PardisoSolver().factorize(A1)
        with _Fases() as f:
            f2 = PardisoSolver().factorize(A2, reuse_analysis_of=f1)
        self.assertEqual(f.fases, [22])
        np.testing.assert_allclose(A2 @ f2.solve(b), b, rtol=1e-13)

    def test_con_entradas_nuevas_se_analiza_la_union_que_sirve_a_las_siguientes(self):
        A1 = _matriz()
        extra = A1.tolil()
        extra[0, 59] = 0.25                                     # fuera del patrón de A1
        A2 = extra.tocsr()
        A3 = A1.copy()
        A3.data *= 1.5
        b = np.ones(60)
        f1 = PardisoSolver().factorize(A1)
        with _Fases() as f:
            f2 = PardisoSolver().factorize(A2, reuse_analysis_of=f1)
            np.testing.assert_allclose(A2 @ f2.solve(b), b, rtol=1e-13)
            f3 = PardisoSolver().factorize(A3, reuse_analysis_of=f2)   # A1 ⊂ unión
        self.assertEqual([x for x in f.fases if x != 33], [12, 22])
        self.assertFalse(f2.shares_analysis)
        self.assertTrue(f3.shares_analysis)
        np.testing.assert_allclose(A3 @ f3.solve(b), b, rtol=1e-13)

    def test_un_error_en_la_fase_numerica_se_repite_con_analisis_propio(self):
        f1 = PardisoSolver().factorize(_matriz())
        handle = PardisoSolver._handle()
        original = handle.call

        def falla_22(phase, *a, **k):
            if phase == 22:
                raise PardisoError(-4, 22)
            return original(phase, *a, **k)
        with mock.patch.object(handle, "call", side_effect=falla_22):
            f2 = PardisoSolver().factorize(_matriz(escala=2.0), reuse_analysis_of=f1)
        self.assertFalse(f2.shares_analysis)
        b = np.ones(60)
        np.testing.assert_allclose(_matriz(escala=2.0) @ f2.solve(b), b, rtol=1e-12)

    def test_pivotes_perturbados_con_analisis_ajeno_repiten_con_analisis_propio(self):
        handle = PardisoSolver._handle()
        f1 = PardisoSolver().factorize(_matriz())
        original = handle.get_iparm
        with _Fases() as f, mock.patch.object(
                handle, "get_iparm", side_effect=lambda i: 1 if i == 14 else original(i)):
            f2 = PardisoSolver().factorize(_matriz(escala=2.0), reuse_analysis_of=f1)
        self.assertEqual(f.fases, [22, 12])
        self.assertFalse(f2.shares_analysis)

    def test_la_solucion_con_analisis_reutilizado_es_exacta_al_redondeo(self):
        A1, A2 = _matriz(escala=1.0), _matriz(escala=40.0)
        b = np.linspace(-1.0, 2.0, 60)
        f2 = PardisoSolver().factorize(A2, reuse_analysis_of=PardisoSolver().factorize(A1))
        self.assertTrue(f2.shares_analysis)
        x = f2.solve(b)
        self.assertLess(np.linalg.norm(A2 @ x - b) / np.linalg.norm(b), 1e-13)

    def test_los_pivotes_son_los_de_un_analisis_propio_y_el_factor_sigue_vigente(self):
        A1, A2 = _matriz(escala=1.0), _matriz(escala=40.0)
        propio = PardisoSolver().factorize(A2).n_small_pivots
        f2 = PardisoSolver().factorize(A2, reuse_analysis_of=PardisoSolver().factorize(A1))
        with _Fases() as f:
            n = f2.n_small_pivots
        self.assertEqual(f.fases, [12])                          # análisis propio al pedirlos
        self.assertEqual(n, propio)
        self.assertTrue(f2.is_current)
        b = np.ones(60)
        np.testing.assert_allclose(A2 @ f2.solve(b), b, rtol=1e-12)


@unittest.skipUnless(_HAS_PARDISO, "pypardiso no instalado")
class TestIgualAPypardiso(unittest.TestCase):
    """Bit a bit con la MKL en un hilo (con varios hilos la MKL no es
    determinista ni consigo misma)."""

    def test_fase_12_y_33_identicas_a_pypardiso(self):
        from pypardiso import PyPardisoSolver
        handle = PardisoSolver._handle()
        # Interfaz C (por valor); la de minúsculas es la de Fortran (puntero).
        set_local = handle.libmkl.MKL_Set_Num_Threads_Local
        set_local.argtypes = [ctypes.c_int]
        set_local.restype = ctypes.c_int
        previo = set_local(1)
        try:
            # El mismo pypardiso que usaba Solidum, con la biblioteca ya
            # localizada (evita recorrer el disco) y los mismos parámetros.
            with mock.patch.dict(os.environ, {"PYPARDISO_MKL_RT": handle.libmkl._name}):
                pyp = PyPardisoSolver()
            init = pyp.libmkl.pardisoinit
            init.restype = None
            init(pyp.pt.ctypes.data_as(ctypes.c_void_p), ctypes.byref(ctypes.c_int32(11)),
                 pyp.iparm.ctypes.data_as(ctypes.c_void_p))
            pyp.set_iparm(56, 1)
            A = _matriz(n=300, seed=7, escala=50.0)
            b = np.sin(np.arange(300.0))
            pyp.factorize(A.copy())
            x_pyp = pyp.solve(A.copy(), b)
            pyp.free_memory(everything=True)
            x = PardisoSolver().factorize(A.copy()).solve(b)
        finally:
            set_local(previo)
        self.assertTrue(np.array_equal(x, x_pyp))


if __name__ == "__main__":
    unittest.main()
