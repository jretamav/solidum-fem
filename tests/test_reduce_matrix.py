"""Reducción de matrices por el operador de las restricciones (ADR 0040, puntos 3 y 5).

1. Con sólo apoyos, ``reduce_matrix`` extrae la submatriz de los DOF libres y
   es idéntica, arreglo por arreglo, a ``Tᵀ·K·T`` en forma canónica CSR,
   también con ceros explícitos en ``K`` y con una fila libre nula; y en
   ``Assembler.reduce`` de un modelo real.
2. El camino que toma: la submatriz con sólo apoyos (y con un esclavo que el
   cierre deja sin maestros); el producto con una restricción multipunto. Las
   dos ramas dan lo mismo.
3. La salida es siempre CSR.
4. Las restricciones se leen en cada llamada: un apoyo añadido entre dos
   reducciones cambia la segunda (el contacto las cambiará entre iteraciones).
5. Ningún solver ni el ensamblador reducen una matriz por su cuenta: todo
   ``Tᵀ·X·T`` pasa por ``reduce_matrix``.
6. ``ConstraintSet.free_dofs``: sólo lectura, invalidado al añadir una
   restricción, con ``ndof`` como clave.
"""
import re
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import scipy.sparse as sp

import solidum
from solidum.bc import constraints as cmod
from solidum.bc.constraints import ConstraintSet, reduce_matrix

RAIZ = Path(__file__).resolve().parents[1]
ORIGINAL = cmod._selected_rows


def _matriz(n=30, seed=0):
    """Matriz dispersa no simétrica con ceros explícitos guardados."""
    rng = np.random.default_rng(seed)
    K = sp.random(n, n, density=0.25, random_state=seed, format="csr") + sp.eye(n, format="csr")
    K = K.tocsr()
    K.data[rng.choice(K.nnz, size=8, replace=False)] = 0.0     # ceros explícitos
    return K


def _producto(K, T):
    return sp.csr_matrix(T.T @ K @ T)


def _identicas(a, b):
    a, b = sp.csr_matrix(a), sp.csr_matrix(b)
    return (a.shape == b.shape and np.array_equal(a.indptr, b.indptr)
            and np.array_equal(a.indices, b.indices) and np.array_equal(a.data, b.data))


def _apoyos(n, dofs):
    cs = ConstraintSet()
    for d in dofs:
        cs.add_dirichlet(d, 0.0)
    return cs


class TestSubmatrizIgualAlProducto(unittest.TestCase):

    def test_con_ceros_explicitos_es_identica_arreglo_por_arreglo(self):
        K = _matriz()
        self.assertGreater(int(np.sum(K.data == 0.0)), 0)
        T, _ = _apoyos(30, [0, 1, 7, 15, 29]).build(30)
        red = reduce_matrix(K, T)
        self.assertTrue(red.has_canonical_format)
        self.assertEqual(int(np.sum(red.data == 0.0)), 0)
        self.assertTrue(_identicas(red, _producto(K, T)))

    def test_con_una_fila_libre_nula_es_identica(self):
        """Un DOF libre con rigidez nula deja una fila vacía en las dos ramas
        (la red de seguridad la detecta después como singular)."""
        K = _matriz(seed=3).tolil()
        K[4, :] = 0.0
        K[:, 4] = 0.0
        K = K.tocsr()
        T, _ = _apoyos(30, [0, 29]).build(30)
        red = reduce_matrix(K, T)
        self.assertEqual(red.indptr[4] - red.indptr[3], 0)       # libre 4 → fila 3
        self.assertTrue(_identicas(red, _producto(K, T)))

    def test_entrada_csc_y_con_duplicados_da_la_forma_canonica(self):
        K = _matriz(seed=5)
        T, _ = _apoyos(30, [2, 3]).build(30)
        coo = K.tocoo()
        duplicada = sp.csr_matrix(
            (np.r_[coo.data, 1.5], (np.r_[coo.row, 6], np.r_[coo.col, 6])), shape=K.shape)
        duplicada.has_canonical_format = False
        for entrada in (K.tocsc(), duplicada):
            with self.subTest(formato=entrada.format):
                self.assertTrue(_identicas(reduce_matrix(entrada, T), _producto(entrada, T)))

    def test_en_assembler_reduce_de_un_modelo_real(self):
        modelo = _modelo_quad4()
        a = modelo.assembler
        K, _ = a.assemble_non_linear_system(np.zeros(a.domain.total_dofs))
        T, _ = a.constraint_set.build(a.ndof)
        K_red, _, _, _ = a.reduce(K, np.zeros(a.ndof))
        self.assertEqual(K_red.format, "csr")
        self.assertTrue(_identicas(K_red, _producto(K, T)))


class TestCamino(unittest.TestCase):

    def _ruta(self, K, T):
        devueltos = []

        def espia(T_):
            r = cmod._selected_rows.__wrapped__(T_) if hasattr(cmod._selected_rows, "__wrapped__") \
                else ORIGINAL(T_)
            devueltos.append(r)
            return r
        with mock.patch.object(cmod, "_selected_rows", side_effect=espia):
            red = reduce_matrix(K, T)
        return devueltos[0], red

    def test_solo_apoyos_toma_la_submatriz(self):
        K = _matriz()
        T, _ = _apoyos(30, [0, 5]).build(30)
        filas, red = self._ruta(K, T)
        self.assertIsNotNone(filas)
        np.testing.assert_array_equal(filas, np.setdiff1d(np.arange(30), [0, 5]))
        self.assertTrue(_identicas(red, _producto(K, T)))

    def test_esclavo_sin_maestros_tras_el_cierre_toma_la_submatriz(self):
        """u_3 = 2·u_0 con u_0 apoyado: el cierre deja a 3 sin maestros (fila
        nula de T), y la submatriz sigue siendo exacta."""
        cs = _apoyos(30, [0])
        cs.add_linear(3, [0], [2.0])
        T, _ = cs.build(30)
        filas, red = self._ruta(_matriz(), T)
        self.assertIsNotNone(filas)
        self.assertTrue(_identicas(red, _producto(_matriz(), T)))

    def test_con_una_mpc_toma_el_producto(self):
        cs = _apoyos(30, [0])
        cs.add_linear(3, [4, 5], [0.5, 0.5])
        T, _ = cs.build(30)
        filas, red = self._ruta(_matriz(), T)
        self.assertIsNone(filas)
        self.assertEqual(red.format, "csr")
        self.assertTrue(_identicas(red, _producto(_matriz(), T)))


class TestRestriccionesVigentes(unittest.TestCase):

    def test_un_apoyo_anadido_entre_dos_reducciones_cambia_la_segunda(self):
        modelo = _modelo_quad4()
        a = modelo.assembler
        K, _ = a.assemble_non_linear_system(np.zeros(a.domain.total_dofs))
        n1 = a.reduce(K, np.zeros(a.ndof))[0].shape[0]
        nodo = next(n for n in a.domain.nodes.values() if not n.boundary_conditions)
        nodo.fix_dof("uy", 0.0)
        n2 = a.reduce(K, np.zeros(a.ndof))[0].shape[0]
        self.assertEqual(n2, n1 - 1)


class TestUnaSolaReduccion(unittest.TestCase):

    def test_ningun_solver_reduce_una_matriz_por_su_cuenta(self):
        patron = re.compile(r"\b(\w+)\.T\s*@\s*[\w\[\]\.]+\s*@\s*\1\b")
        archivos = sorted((RAIZ / "solidum" / "math" / "solvers").glob("*.py"))
        archivos.append(RAIZ / "solidum" / "math" / "assembly.py")
        hallados = [f"{p.name}: {m.group(0)}" for p in archivos
                    for m in patron.finditer(p.read_text(encoding="utf-8"))]
        self.assertEqual(hallados, [], "use reduce_matrix (ADR 0040)")


class TestDofsLibres(unittest.TestCase):

    def test_solo_lectura(self):
        libres = _apoyos(10, [2]).free_dofs(10)
        with self.assertRaises(ValueError):
            libres[0] = 7

    def test_anadir_una_restriccion_la_invalida(self):
        cs = _apoyos(10, [2])
        np.testing.assert_array_equal(cs.free_dofs(10), [0, 1, 3, 4, 5, 6, 7, 8, 9])
        cs.add_dirichlet(5, 1.0)
        np.testing.assert_array_equal(cs.free_dofs(10), [0, 1, 3, 4, 6, 7, 8, 9])
        cs.add_linear(7, [0], [1.0])
        np.testing.assert_array_equal(cs.free_dofs(10), [0, 1, 3, 4, 6, 8, 9])

    def test_un_dof_negativo_sin_numerar_es_error(self):
        """−1 es el número de un DOF sin numerar: marcaría otro DOF como
        restringido en silencio."""
        cs = ConstraintSet()
        with self.assertRaisesRegex(ValueError, "negativo"):
            cs.add_dirichlet(-1, 0.0)
        with self.assertRaisesRegex(ValueError, "negativo"):
            cs.add_linear(3, [-1], [1.0])
        self.assertEqual(len(cs), 0)

    def test_ndof_es_la_clave(self):
        cs = _apoyos(10, [2])
        self.assertEqual(cs.free_dofs(10).size, 9)
        self.assertEqual(cs.free_dofs(12).size, 11)
        self.assertEqual(cs.free_dofs(10).size, 9)

    def test_coincide_con_las_columnas_de_T(self):
        cs = _apoyos(10, [2, 4])
        cs.add_linear(6, [0], [0.5])
        T, _ = cs.build(10)
        libres = cs.free_dofs(10)
        self.assertEqual(T.shape[1], libres.size)
        self.assertEqual(libres.dtype, np.int64)
        np.testing.assert_array_equal(T[libres].toarray(), np.eye(libres.size))


def _modelo_quad4():
    import tempfile
    d = Path(tempfile.mkdtemp())
    ruta = d / "m.yaml"
    ruta.write_text(
        "materials:\n"
        "  - {id: 1, type: Elastic2D, E: 1000.0, nu: 0.25}\n"
        "mesh_element: Quad4\n"
        "mesh_thickness: 1.0\n"
        "mesh:\n"
        "  structured:\n"
        "    blocks:\n"
        "      - {material: 1, x: {length: 6, n: 6}, y: {length: 2, n: 2}}\n"
        "boundary_conditions_by_group:\n"
        "  xmin: {ux: 0.0, uy: 0.0}\n"
        "point_loads_by_group:\n"
        "  xmax: {uy: -1.0}\n"
        "solver: {type: LinearSolver}\n", encoding="utf-8")
    return solidum.load_yaml(ruta)


if __name__ == "__main__":
    unittest.main()
