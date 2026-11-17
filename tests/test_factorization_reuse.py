"""Factorización reutilizada por el corrector de Newton (ADR 0040, punto 2).

El corrector reutiliza la última factorización si la matriz es idéntica bit a
bit y factoriza si no. Invariantes:

1. Misma matriz (otro objeto) → una factorización; un valor distinto en un
   bit, un cero explícito frente a una entrada ausente o la matriz anterior
   modificada en sitio → otra factorización, con la solución de la matriz
   nueva.
2. Un inspector pendiente (red de seguridad) recibe el factor reutilizado.
3. La degradación a LU y ``release`` vacían lo retenido; ``reset_step``
   suelta sólo la congelada.
4. Un factor de Pardiso que otro consumidor reemplazó no se usa: se
   factoriza otra vez su matriz (la congelada del Newton modificado sigue
   siendo la de su iteración).
5. Toda resolución vigila la finitud: una solución no finita es tangente
   singular (``RuntimeError``), también en la tangente congelada.
6. En un análisis: elástico → una factorización por análisis (estático y
   arco); IMPL-EX → tantas como matrices distintas consecutivas; y con LU los
   desplazamientos son idénticos bit a bit con y sin reutilización (daño
   IMPL-EX y J2).
7. Todo solver registrado con corrector lo suelta al empezar y al terminar
   ``solve`` (``releases_corrector``).
"""
import inspect
import re
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import scipy.sparse as sp

import solidum
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.linalg import LUSolver
from solidum.math.linalg.dispatcher import _HAS_PARDISO
from solidum.math.linalg.lu import LUFactorized
from solidum.math.solvers import corrector as corrector_module
from solidum.math.solvers.corrector import NewtonCorrector
from solidum.registry import SolverRegistry

RAIZ = Path(__file__).resolve().parents[1]


def _spd(n=40, seed=0):
    M = sp.random(n, n, density=0.15, random_state=seed, format="csr")
    return (M @ M.T + n * sp.eye(n)).tocsr()


def _corrector(backend="lu", freeze=None):
    return NewtonCorrector(ConvergenceCriterion(), max_iter=5, is_symmetric=True,
                           is_positive_definite=False, linear_algebra=backend,
                           freeze_tangent_after_iter=freeze)


class _Contador:
    """Cuenta las factorizaciones de un backend sin cambiar lo que hace."""

    def __init__(self, cls):
        self.cls, self.n = cls, 0
        self.original = cls.factorize

    def __enter__(self):
        original, contador = self.original, self

        def factorize(backend, K, **kw):
            contador.n += 1
            return original(backend, K, **kw)
        self._patch = mock.patch.object(self.cls, "factorize", factorize)
        self._patch.start()
        return self

    def __exit__(self, *exc):
        self._patch.stop()


class TestMatrizIdentica(unittest.TestCase):

    def setUp(self):
        self.A = _spd()
        self.b = np.arange(1.0, 41.0)

    def test_la_misma_matriz_en_otro_objeto_se_factoriza_una_vez(self):
        c = _corrector()
        with _Contador(LUSolver) as k:
            x1 = c.solve(self.A, self.b)
            x2 = c.solve(self.A.copy(), 2 * self.b)
        self.assertEqual(k.n, 1)
        np.testing.assert_allclose(self.A @ x1, self.b, rtol=1e-12)
        np.testing.assert_allclose(self.A @ x2, 2 * self.b, rtol=1e-12)

    def test_un_bit_distinto_vuelve_a_factorizar(self):
        B = self.A.copy()
        B.data[3] = np.nextafter(B.data[3], np.inf)
        c = _corrector()
        with _Contador(LUSolver) as k:
            c.solve(self.A, self.b)
            x = c.solve(B, self.b)
        self.assertEqual(k.n, 2)
        np.testing.assert_allclose(B @ x, self.b, rtol=1e-12)

    def test_un_cero_explicito_no_equivale_a_una_entrada_ausente(self):
        sin_ceros = self.A.copy()
        i, j = next((i, j) for i in range(40) for j in range(40) if sin_ceros[i, j] == 0.0)
        coo = sin_ceros.tocoo()
        con_ceros = sp.csr_matrix((np.r_[coo.data, 0.0], (np.r_[coo.row, i], np.r_[coo.col, j])),
                                  shape=sin_ceros.shape)
        self.assertEqual(con_ceros.nnz, sin_ceros.nnz + 1)
        c = _corrector()
        with _Contador(LUSolver) as k:
            c.solve(sin_ceros, self.b)
            c.solve(con_ceros, self.b)
        self.assertEqual(k.n, 2)

    def test_la_matriz_anterior_modificada_en_sitio_vuelve_a_factorizar(self):
        A = self.A.copy()
        c = _corrector()
        with _Contador(LUSolver) as k:
            c.solve(A, self.b)
            A.data *= 2.0
            x = c.solve(A, self.b)
        self.assertEqual(k.n, 2)
        np.testing.assert_allclose(A @ x, self.b, rtol=1e-12)


class TestReglas(unittest.TestCase):

    def setUp(self):
        self.A = _spd(seed=1)
        self.b = np.ones(40)

    def test_el_inspector_pendiente_recibe_el_factor_reutilizado(self):
        c = _corrector()
        c.solve(self.A, self.b)
        recibidos = []
        c.inspect_next_factorization(lambda f, A: recibidos.append((f, A)))
        with _Contador(LUSolver) as k:
            c.solve(self.A.copy(), self.b)
            c.solve(self.A.copy(), self.b)
        self.assertEqual(k.n, 0)
        self.assertEqual(len(recibidos), 1)              # una vez, y se olvida
        self.assertIsInstance(recibidos[0][0], LUFactorized)

    def test_degradar_a_lu_y_release_vacian_lo_retenido(self):
        for vaciar in ("degradar", "release"):
            with self.subTest(vaciar=vaciar):
                c = _corrector(freeze=1)
                c.solve(self.A, self.b)
                self.assertIsNotNone(c.frozen_factor)
                if vaciar == "degradar":
                    c._degrade_to_lu(40)
                else:
                    c.release()
                self.assertIsNone(c.frozen_factor)
                with _Contador(LUSolver) as k:
                    c.solve(self.A, self.b)
                self.assertEqual(k.n, 1)

    def test_reset_step_suelta_la_congelada_y_conserva_la_ultima(self):
        c = _corrector(freeze=1)
        c.solve(self.A, self.b)
        c.reset_step()
        self.assertIsNone(c.frozen_factor)
        with _Contador(LUSolver) as k:
            c.solve(self.A.copy(), self.b)
        self.assertEqual(k.n, 0)

    def test_con_otra_forma_la_congelada_no_sirve_y_se_factoriza_la_nueva(self):
        """El contacto puede cambiar el número de incógnitas entre
        iteraciones: la tangente congelada de otra forma no se usa."""
        c = _corrector(freeze=1)
        c.solve(self.A, self.b, iteration=0)
        B = _spd(n=38, seed=4)
        x = c.solve(B, np.ones(38), iteration=1)
        np.testing.assert_allclose(B @ x, np.ones(38), rtol=1e-12)

    def test_releases_corrector_suelta_todos_los_correctores_del_solver(self):
        A, b = self.A, self.b

        class Escalonado:
            def __init__(self):
                self.mecanico, self.termico = _corrector(), _corrector()

            @corrector_module.releases_corrector
            def solve(self):
                self.mecanico.solve(A, b)
                self.termico.solve(A, b)
                return self.mecanico._last is not None and self.termico._last is not None

        s = Escalonado()
        self.assertTrue(s.solve())                           # retenidos durante el análisis
        self.assertIsNone(s.mecanico._last)
        self.assertIsNone(s.termico._last)

    def test_solucion_no_finita_es_tangente_singular_tambien_congelada(self):
        nan = mock.patch.object(LUFactorized, "solve", lambda self, b: np.full_like(b, np.nan))
        for freeze in (None, 1):
            with self.subTest(freeze=freeze):
                c = _corrector(freeze=freeze)
                with nan, self.assertRaises(RuntimeError):
                    c.solve(self.A, self.b)
        c = _corrector(freeze=1)
        c.solve(self.A, self.b)
        with nan, self.assertRaises(RuntimeError):
            c.solve(2 * self.A, self.b, iteration=1)          # rama congelada


@unittest.skipUnless(_HAS_PARDISO, "pypardiso no instalado")
class TestFactorReemplazadoEnPardiso(unittest.TestCase):

    def setUp(self):
        from solidum.math.linalg.pardiso import PardisoSolver
        self.Pardiso = PardisoSolver
        self.A, self.B = _spd(seed=2), _spd(seed=3)
        self.b = np.ones(40)

    def test_no_se_reutiliza_un_factor_que_otro_reemplazo(self):
        c = _corrector("pardiso")
        c.solve(self.A, self.b)
        self.Pardiso().factorize(self.B)                     # otro consumidor
        with _Contador(self.Pardiso) as k:
            x = c.solve(self.A.copy(), self.b)
        self.assertEqual(k.n, 1)
        np.testing.assert_allclose(self.A @ x, self.b, rtol=1e-12)

    def test_la_congelada_reemplazada_se_refactoriza_con_su_matriz(self):
        c = _corrector("pardiso", freeze=1)
        c.solve(self.A, self.b, iteration=0)
        self.Pardiso().factorize(self.B)
        x = c.solve(self.B, self.b, iteration=1)              # Newton modificado: A
        np.testing.assert_allclose(self.A @ x, self.b, rtol=1e-12)

    def test_un_inspector_que_factoriza_otra_matriz_no_deja_un_factor_invalido(self):
        """La red de seguridad puede factorizar otra matriz (``_shifted_factor``)
        en el handle único: el corrector lo comprueba después del inspector."""
        for previa in (False, True):                         # fresca y reutilizada
            with self.subTest(reutilizada=previa):
                c = _corrector("pardiso")
                if previa:
                    c.solve(self.A, self.b)
                c.inspect_next_factorization(lambda f, A: self.Pardiso().factorize(self.B))
                x = c.solve(self.A.copy(), self.b)
                np.testing.assert_allclose(self.A @ x, self.b, rtol=1e-12)

    def test_se_factoriza_la_copia_y_no_la_matriz_del_llamador(self):
        """Quien pasó la matriz la modifica en sitio y presenta después la
        original: el factor reutilizado (y el refinamiento de Pardiso, que
        relee la matriz en cada resolución) sigue siendo el de la original."""
        A = self.A.copy()
        c = _corrector("pardiso")
        x1 = c.solve(A, self.b)
        A.data *= 3.0                                        # búfer reutilizado
        x2 = c.solve(self.A.copy(), self.b)
        self.assertTrue(np.array_equal(x1, x2))


def _yaml(texto):
    d = Path(tempfile.mkdtemp())
    ruta = d / "m.yaml"
    ruta.write_text(texto, encoding="utf-8")
    return ruta


def _bloque(material, solver, desplazamiento=0.05, n=(6, 2), fuerza=None):
    """Bloque de 60×20 en tracción: desplazamiento impuesto en xmax o, con
    ``fuerza``, una carga en xmax (para el arco)."""
    if fuerza is None:
        borde, cargas = f"  xmax: {{ux: {desplazamiento}}}\n", ""
    else:
        borde, cargas = "", f"point_loads_by_group:\n  xmax: {{ux: {fuerza}}}\n"
    return (f"""\
materials:
  - {material}
mesh_element: Quad4
mesh_thickness: 10.0
mesh:
  structured:
    blocks:
      - {{material: 1, x: {{length: 60, n: {n[0]}}}, y: {{length: 20, n: {n[1]}}}}}
    groups:
      Origen: "near(x, 0) and near(y, 0)"
boundary_conditions_by_group:
  xmin: {{ux: 0.0}}
  Origen: {{uy: 0.0}}
""" + borde + cargas + f"solver: {solver}\noutput: {{log_level: ERROR}}\n")


ELASTICO = "{id: 1, type: Elastic2D, E: 30000.0, nu: 0.2}"
DANO_IMPLEX = ("{id: 1, type: IsotropicDamageCrackBand2D, E: 30000.0, nu: 0.0, f_t: 3.0, "
               "G_f: 0.1, integration: impl-ex}")
J2 = "{id: 1, type: VonMises2D, E: 70.0, nu: 0.3, sigma_y: 0.243, H: 0.135, hypothesis: plane_strain}"


class TestEnUnAnalisis(unittest.TestCase):

    def _contar(self, texto):
        with _Contador(LUSolver) as k:
            res = solidum.run_yaml(_yaml(texto))
        return k.n, res

    def test_elastico_una_factorizacion_por_analisis(self):
        casos = {
            "estático": _bloque(ELASTICO, "{type: NonlinearSolver, num_steps: 5, linear_algebra: lu}"),
            "arco": _bloque(ELASTICO, "{type: ArcLengthSolver, max_lambda: 1.0, initial_dlambda: 0.2, "
                            "linear_algebra: lu}", fuerza=1000.0),
        }
        for nombre, texto in casos.items():
            with self.subTest(solver=nombre):
                n, res = self._contar(texto)
                self.assertEqual(n, 1)
                self.assertGreater(res.num_steps, 1)

    def test_implex_factoriza_sólo_las_matrices_distintas(self):
        """Newton estándar con IMPL-EX: cada paso resuelve dos veces con la
        misma tangente; se factoriza sólo cuando la matriz cambia."""
        vistas = []
        original = NewtonCorrector.solve

        def solve(self, A, b, iteration=0):
            vistas.append(A.copy())
            return original(self, A, b, iteration)
        with mock.patch.object(NewtonCorrector, "solve", solve):
            n, _ = self._contar(_bloque(
                DANO_IMPLEX, "{type: NonlinearSolver, num_steps: 40, linear_algebra: lu}", 0.04))
        distintas = 1 + sum(not corrector_module._same_matrix(b, a) for a, b in zip(vistas, vistas[1:]))
        self.assertEqual(n, distintas)
        self.assertLess(n, len(vistas))
        self.assertGreater(n, 2)                              # el daño cambió la tangente

    def test_con_lu_los_desplazamientos_son_identicos_con_y_sin_reutilizacion(self):
        casos = {
            "daño IMPL-EX": _bloque(DANO_IMPLEX, "{type: NonlinearSolver, num_steps: 40, "
                                    "linear_algebra: lu}", 0.04),
            "J2": _bloque(J2, "{type: NonlinearSolver, num_steps: 10, linear_algebra: lu}", 0.6),
        }
        for nombre, texto in casos.items():
            with self.subTest(caso=nombre):
                U_con = solidum.run_yaml(_yaml(texto)).U
                with mock.patch.object(corrector_module, "_same_matrix", lambda A, B: False):
                    U_sin = solidum.run_yaml(_yaml(texto)).U
                self.assertTrue(np.array_equal(U_con, U_sin))


J2_DINAMICO = ("{id: 1, type: VonMises2D, E: 70.0, nu: 0.3, sigma_y: 0.243, H: 0.135, "
               "hypothesis: plane_strain, density: 1.0}")


@unittest.skipUnless(_HAS_PARDISO, "pypardiso no instalado")
class TestCadenaDelCorrectorConPardiso(unittest.TestCase):
    """La optimización principal (ADR 0040, 1b): dentro de un análisis el
    corrector encadena sus factorizaciones y Pardiso repite sólo la numérica
    (fase 22); la primera de cada análisis hace su propio análisis (fase 12),
    y la Jacobiana de Newton-Newmark no hereda el de la matriz de masa."""

    def _fases(self, texto):
        from solidum.math.linalg.pardiso import PardisoSolver
        handle = PardisoSolver._handle()
        original, fases = handle.call, []

        def call(phase, *a, **k):
            if phase in (11, 12, 22):
                fases.append(phase)
            return original(phase, *a, **k)
        with mock.patch.object(handle, "call", side_effect=call):
            solidum.run_yaml(_yaml(texto))
        return fases

    def test_j2_estatico_casi_todas_reutilizan_el_analisis(self):
        """Malla estructurada: la plasticidad rompe simetrías y aparecen
        entradas que antes valían cero; el patrón de la cadena crece unas
        pocas veces (análisis propio) y el resto reutiliza el análisis
        (medido: 4 análisis y 24 numéricas; sin la cadena, 28 análisis)."""
        texto = _bloque(J2, "{type: NonlinearSolver, num_steps: 10, linear_algebra: pardiso}", 0.6)
        for corrida in (1, 2):                                # la segunda empieza otra cadena
            with self.subTest(corrida=corrida):
                fases = self._fases(texto)
                self.assertEqual(fases[0], 12)
                self.assertGreaterEqual(fases.count(22), 3 * fases.count(12), fases)

    def test_pivotes_perturbados_rompen_la_cadena_hasta_el_proximo_analisis(self):
        """Con una matriz singular o casi (un mecanismo), la solución depende
        de la escala de la factorización aunque su residuo sea pequeño:
        medido, una armadura corrotacional con mecanismo «convergía» a un
        colapso con el análisis reutilizado. En cuanto una factorización
        perturba pivotes, el resto del análisis no reutiliza el análisis."""
        previos, perturbados = [], [1, 0, 0, 0, 0]

        class Backend:
            name, reuses_analysis = "prueba", True

            def factorize(self, K, reuse_analysis_of=None):
                previos.append(reuse_analysis_of)
                f = LUSolver().factorize(K)
                f.n_perturbed_pivots = perturbados.pop(0) if perturbados else 0
                return f

        c = _corrector()
        with mock.patch.object(corrector_module, "select_solver", lambda *a, **k: Backend()):
            for semilla in range(3):
                c.solve(_spd(seed=10 + semilla), np.ones(40))
            c.release()                                     # análisis siguiente
            for semilla in range(3):
                c.solve(_spd(seed=20 + semilla), np.ones(40))
        rotos, nuevo = previos[:3], previos[3:]
        self.assertEqual(rotos, [None, None, None])         # la 1.ª perturbó: sin cadena
        self.assertIsNone(nuevo[0])                          # cada análisis empieza la suya
        self.assertIsNotNone(nuevo[1])
        self.assertIsNotNone(nuevo[2])

    def test_newton_newmark_la_primera_jacobiana_no_hereda_la_masa(self):
        texto = _bloque(J2_DINAMICO, "{type: NewtonNewmarkSolver, t_end: 2.0, dt: 0.5, "
                        "linear_algebra: pardiso}", fuerza=30.0)
        fases = self._fases(texto)
        self.assertEqual(fases[:2], [12, 12], fases)         # M (fuera de la cadena) y la 1.ª J
        self.assertIn(22, fases)


class TestTodoSolverConCorrectorLoSuelta(unittest.TestCase):

    def test_releases_corrector_en_todo_solver_registrado_con_corrector(self):
        con_corrector = []
        for nombre in SolverRegistry.names():
            cls = SolverRegistry.get(nombre)
            fuentes = "".join(inspect.getsource(k) for k in cls.__mro__
                              if k.__module__.startswith("solidum") and k is not object)
            if "NewtonCorrector(" in fuentes:
                con_corrector.append(nombre)
                self.assertTrue(getattr(cls.solve, "__releases_corrector__", False), nombre)
        self.assertGreaterEqual(len(con_corrector), 5)

    def test_el_arco_por_disipacion_del_modulo_de_usuario_lo_suelta(self):
        """Sin importar el módulo (no se registra en esta sesión de pruebas)."""
        fuente = (RAIZ / "solidum" / "user" / "dissipation_control"
                  / "dissipation_arclength.py").read_text(encoding="utf-8")
        self.assertRegex(fuente, re.compile(r"@releases_corrector\s+def solve\(", re.S))


if __name__ == "__main__":
    unittest.main()
