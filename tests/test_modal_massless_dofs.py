"""Análisis modal con grados de libertad sin masa (ADR 0024).

Con masas sólo en algunos nudos o elementos de ``density = 0``, la matriz de
masa es semidefinida. ARPACK fallaba (``ArpackError -9999``, sin causa) y el
camino denso exigía ``M`` definida. Ahora esos grados de libertad se condensan
estáticamente, de forma exacta, y sus componentes de los modos se recuperan.

Criterios de aceptación (spec de ``ModalSolver``, bloque ``acceptance``):

1. Marco de un nivel con columnas de densidad cero (18 de sus 33 grados de
   libertad libres sin masa): ``T₁ = 0.03097 s`` y ``T₂ = 0.01295 s``,
   diferencia relativa ≤ 1e-6 contra la solución densa (valores propios de
   ``K⁻¹M``, otro algoritmo), y modos con residuo relativo
   ``‖Kφ − ω²Mφ‖/‖Kφ‖ ≤ 1e-8``.
2. Si el problema no tiene solución, un error tipado de Solidum
   (``MasslessDofsError``) que nombra los grados de libertad sin masa.
3. Los modelos con masa en todos los grados de libertad dan los mismos
   resultados que antes: siguen el camino de siempre (y la suite completa
   lo confirma).

Se cubren los dos caminos del ``EigenSolver`` condensado (denso con pocos
grados de libertad con masa, ARPACK con operadores con muchos) y las
estructuras libres (``K`` singular).
"""
import itertools
import math
import unittest
from unittest import mock

import numpy as np
import scipy.sparse as sp

from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.linalg import EigenSolver
from solidum.math.linalg.eigen import SemidefiniteMassError, massless_dofs
from solidum.math.solvers import ModalSolver
from solidum.math.solvers.diagnostics import MasslessDofsError
from solidum.sections import FiberSection, RectangularPatch

RTOL_ACCEPT = 1.0e-6
RESIDUAL_ACCEPT = 1.0e-8


def one_story_frame(n_el: int, *, column_density=0.0, beam_density=7850.0):
    """Marco de acero de un nivel: h = 3 m, L = 5 m, columnas 0.30 × 0.50 m
    (20 fibras), viga 0.40 × 0.70 m (28 fibras), 250 MPa, 3 puntos de Gauss,
    bases empotradas. Con ``column_density = 0`` sólo la viga tiene masa."""
    h, span = 3.0, 5.0

    def section(b, d, n, density):
        steel = Elastoplastic1D(E=2.0e11, sigma_y=250.0e6, H=0.0, density=density)
        return FiberSection(patches=[RectangularPatch(steel, y=(-d / 2, d / 2),
                                                      z=(-b / 2, b / 2), ny=n, nz=1)])

    dom = Domain()
    ids, eids = itertools.count(1), itertools.count(1)
    corner = {k: dom.add_node(next(ids), xy) for k, xy in
              {"bl": [0.0, 0.0], "tl": [0.0, h], "tr": [span, h], "br": [span, 0.0]}.items()}

    def member(sec, a, b):
        pts = [a]
        for i in range(1, n_el):
            t = i / n_el
            pts.append(dom.add_node(next(ids), [a.coordinates[0] + t * (b.coordinates[0] - a.coordinates[0]),
                                               a.coordinates[1] + t * (b.coordinates[1] - a.coordinates[1])]))
        pts.append(b)
        for p, q in zip(pts[:-1], pts[1:]):
            dom.add_element(Frame2DEuler(next(eids), [p, q], section=sec, n_ip=3))

    column = section(0.30, 0.50, 20, column_density)
    beam = section(0.40, 0.70, 28, beam_density)
    member(column, corner["bl"], corner["tl"])
    member(beam, corner["tl"], corner["tr"])
    member(column, corner["br"], corner["tr"])
    for k in ("bl", "br"):
        for dof in ("ux", "uy", "rz"):
            corner[k].fix_dof(dof, 0.0)
    dom.generate_equation_numbers()
    return dom


def reduced_pair(dom, lumping="consistent"):
    asm = Assembler(dom)
    asm.assemble_system()
    K_red, M_red, T = asm.reduce_pair(asm.K_global, asm.assemble_mass_matrix(lumping))
    return K_red.toarray(), M_red.toarray(), asm


def dense_periods(K: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Referencia independiente: ``μ = 1/ω²`` son los valores propios de
    ``K⁻¹M`` (no simétrico, ``numpy``); los grados de libertad sin masa dan
    ``μ = 0`` y se descartan."""
    mu = np.real(np.linalg.eigvals(np.linalg.solve(K, M)))
    mu = np.sort(mu[mu > 1e-12 * mu.max()])[::-1]
    return 2.0 * np.pi * np.sqrt(mu)


def assert_modes_ok(test, K, M, result, free):
    """Residuo relativo y M-ortonormalidad de los modos (espacio reducido)."""
    phi = result.modes[free]
    for j, omega in enumerate(result.frequencies_rad):
        Kp = K @ phi[:, j]
        res = np.linalg.norm(Kp - omega ** 2 * (M @ phi[:, j])) / np.linalg.norm(Kp)
        test.assertLessEqual(res, RESIDUAL_ACCEPT, f"modo {j + 1}")
    np.testing.assert_allclose(phi.T @ M @ phi, np.eye(phi.shape[1]), atol=1e-10)


class TestMarcoDeUnNivelSinMasaEnColumnas(unittest.TestCase):
    """Criterio 1."""

    def test_periodos_y_modos(self):
        dom = one_story_frame(4)
        K, M, asm = reduced_pair(dom)
        self.assertEqual(K.shape[0], 33)
        self.assertEqual(int(np.sum(np.diag(M) == 0.0)), 18)
        res = ModalSolver(Assembler(dom), n_modes=2).solve()
        ref = dense_periods(K, M)
        self.assertEqual(ref.size, 15)                       # 15 modos con masa
        np.testing.assert_allclose(res.periods, ref[:2], rtol=RTOL_ACCEPT)
        # Los valores del criterio, a las cifras con que se dieron.
        self.assertAlmostEqual(res.periods[0], 0.03097, delta=5e-6)
        self.assertAlmostEqual(res.periods[1], 0.01295, delta=5e-6)
        assert_modes_ok(self, K, M, res, asm.constraint_set.free_dofs(asm.ndof))

    def test_camino_arpack_con_muchos_grados_de_libertad_con_masa(self):
        """12 elementos por miembro: 39 grados de libertad con masa, más que
        el umbral del camino denso; ARPACK trabaja sobre el problema
        condensado con operadores, sin formar K̂."""
        dom = one_story_frame(12)
        for lumping in ("consistent", "lumped"):
            with self.subTest(lumping=lumping):
                K, M, asm = reduced_pair(dom, lumping)
                n_mass = int(np.sum(np.diag(M) != 0.0))
                self.assertGreater(n_mass, 32)
                with mock.patch.object(EigenSolver, "_condensed_arpack",
                                       autospec=True,
                                       side_effect=EigenSolver._condensed_arpack) as spy:
                    res = ModalSolver(Assembler(dom), n_modes=4, lumping=lumping).solve()
                self.assertEqual(spy.call_count, 1)
                np.testing.assert_allclose(res.periods, dense_periods(K, M)[:4], rtol=RTOL_ACCEPT)
                assert_modes_ok(self, K, M, res, asm.constraint_set.free_dofs(asm.ndof))

    def test_con_masa_en_todos_sigue_el_camino_de_siempre(self):
        """Criterio 3: con columnas con densidad no hay nada que condensar."""
        dom = one_story_frame(4, column_density=7850.0)
        K, M, _ = reduced_pair(dom)
        eig = EigenSolver()
        with mock.patch.object(EigenSolver, "_solve_condensed") as condensed:
            eig.solve(sp.csr_matrix(K), sp.csr_matrix(M), 2)
        condensed.assert_not_called()
        self.assertEqual(eig.last_massless.size, 0)


class TestEstructurasLibres(unittest.TestCase):
    """``K`` singular (modos de cuerpo rígido) y grados de libertad sin masa."""

    def test_dos_masas_unidas_por_una_barra_sin_masa(self):
        """Tres modos rígidos (dos traslaciones y el giro: el movimiento
        transversal relativo de una barra libre sin masa es un giro) y el
        axial, ``ω² = (EA/L)·(1/m₁ + 1/m₂)``."""
        E, A, I, length, m1, m2 = 2.0e11, 4.0e-3, 2.0e-5, 2.5, 3.0e3, 5.0e3
        dom = Domain()
        n1, n2 = dom.add_node(1, [0.0, 0.0]), dom.add_node(2, [length, 0.0])
        dom.add_element(Frame2DEuler(1, [n1, n2], material=Elastic1D(E=E, density=0.0), A=A, I=I))
        dom.add_nodal_mass(1, ux=m1, uy=m1)
        dom.add_nodal_mass(2, ux=m2, uy=m2)
        res = ModalSolver(Assembler(dom), n_modes=4).solve()
        np.testing.assert_allclose(res.frequencies_rad[:3], 0.0, atol=1e-6)
        omega_axial = math.sqrt(E * A / length * (1.0 / m1 + 1.0 / m2))
        self.assertLessEqual(abs(res.frequencies_rad[3] / omega_axial - 1.0), RTOL_ACCEPT)
        self.assertTrue(np.all(np.isinf(res.periods[:3])))

    def test_viga_libre_con_arpack_y_reintento(self):
        """Viga libre de 20 elementos sin masa, con masa en ux y uy de sus 21
        nudos (42 grados de libertad con masa: camino ARPACK). Con σ = 0,
        ``K`` es singular y se reintenta con desplazamiento negativo, como
        en el camino sin condensar. Debe coincidir con el camino denso."""
        dom = Domain()
        nodes = [dom.add_node(i + 1, [0.25 * i, 0.0]) for i in range(21)]
        mat = Elastic1D(E=2.0e11, density=0.0)
        for k in range(20):
            dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], material=mat, A=4.0e-3, I=2.0e-5))
        for i, n in enumerate(nodes):
            dom.add_nodal_mass(n.id, ux=150.0 + 3.0 * i, uy=150.0 + 3.0 * i)
        K, M, asm = reduced_pair(dom)
        K, M = sp.csr_matrix(K), sp.csr_matrix(M)
        with self.assertLogs("solidum.linalg.eigen", level="WARNING") as logs:
            lam_arpack, phi_arpack = EigenSolver().solve(K, M, 6)
        self.assertIn("Reintentando", "\n".join(logs.output))
        lam_dense, _ = EigenSolver().solve(K, M, 41)        # n_modes ≥ n_mass − 1: denso
        np.testing.assert_allclose(lam_arpack[:3], 0.0, atol=1e-6 * lam_arpack[3])
        np.testing.assert_allclose(lam_arpack[3:], lam_dense[3:6], rtol=1e-8)
        np.testing.assert_allclose(phi_arpack.T @ (M @ phi_arpack), np.eye(6), atol=1e-8)


class TestErrorTipado(unittest.TestCase):
    """Criterio 2: ``MasslessDofsError`` nombra los grados de libertad."""

    @staticmethod
    def _cantilever_with_tip_mass():
        dom = Domain()
        nodes = [dom.add_node(i + 1, [0.0, 1.5 * i]) for i in range(3)]
        mat = Elastic1D(E=2.0e11, density=0.0)
        for k in range(2):
            dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], material=mat, A=6.0e-3, I=8.5e-5))
        for dof in ("ux", "uy", "rz"):
            nodes[0].fix_dof(dof, 0.0)
        return dom, nodes

    def test_mas_modos_que_grados_de_libertad_con_masa(self):
        dom, nodes = self._cantilever_with_tip_mass()
        dom.add_nodal_mass(nodes[-1].id, ux=2.5e4)
        with self.assertRaises(MasslessDofsError) as ctx:
            ModalSolver(Assembler(dom), n_modes=2).solve()
        err = ctx.exception
        self.assertIsInstance(err, ValueError)
        self.assertIsInstance(err, SemidefiniteMassError)
        self.assertEqual(err.cause, "too_many_modes")
        self.assertEqual(err.n_mass, 1)
        self.assertEqual(len(err.massless_dofs), 5)
        for label in ("nodo 2 (ux)", "nodo 2 (uy)", "nodo 2 (rz)", "nodo 3 (uy)", "nodo 3 (rz)"):
            self.assertIn(label, err.massless_dofs)
            self.assertIn(label, str(err))
        self.assertIn("n_modes ≤ 1", str(err))

    def test_ningun_grado_de_libertad_con_masa(self):
        dom, _ = self._cantilever_with_tip_mass()
        with self.assertRaises(MasslessDofsError) as ctx:
            ModalSolver(Assembler(dom), n_modes=1).solve()
        self.assertEqual(ctx.exception.cause, "no_mass")

    def test_grados_de_libertad_sin_masa_que_forman_un_mecanismo(self):
        """Dos barras alineadas: el ``uy`` del nudo intermedio no tiene masa
        ni rigidez; su frecuencia está indeterminada."""
        dom = Domain()
        nodes = [dom.add_node(i + 1, [2.0 * i, 0.0]) for i in range(3)]
        mat = Elastic1D(E=2.0e11, density=0.0)
        for k in range(2):
            dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]], mat, A=1.0e-3))
        nodes[0].fix_dof("ux", 0.0)
        nodes[0].fix_dof("uy", 0.0)
        nodes[2].fix_dof("uy", 0.0)
        dom.add_nodal_mass(3, ux=800.0)
        with self.assertRaises(MasslessDofsError) as ctx:
            ModalSolver(Assembler(dom), n_modes=1).solve()
        self.assertEqual(ctx.exception.cause, "massless_mechanism")
        self.assertIn("nodo 2 (uy)", str(ctx.exception))


class TestEigenSolverCondensado(unittest.TestCase):
    """La capa algebraica, con matrices construidas a mano."""

    @staticmethod
    def _problem(n, n_mass, seed):
        rng = np.random.default_rng(seed)
        B = rng.standard_normal((n, n))
        K = B @ B.T + n * np.eye(n)                       # simétrica definida
        d = np.zeros(n)
        mass = np.sort(rng.choice(n, n_mass, replace=False))
        d[mass] = rng.uniform(0.5, 3.0, n_mass)
        return K, np.diag(d)

    def test_contra_la_solucion_densa(self):
        for n, n_mass, n_modes in ((40, 12, 5), (120, 50, 6)):     # denso y ARPACK
            with self.subTest(n=n, n_mass=n_mass):
                K, M = self._problem(n, n_mass, seed=n)
                eig = EigenSolver()
                lam, phi = eig.solve(sp.csr_matrix(K), sp.csr_matrix(M), n_modes)
                self.assertEqual(eig.last_massless.size, n - n_mass)
                mu = np.real(np.linalg.eigvals(np.linalg.solve(K, M)))
                ref = np.sort(1.0 / mu[mu > 1e-12 * mu.max()])[:n_modes]
                np.testing.assert_allclose(lam, ref, rtol=1e-9)
                np.testing.assert_allclose(phi.T @ M @ phi, np.eye(n_modes), atol=1e-10)
                # El residuo de los vectores lo acota la tolerancia de ARPACK
                # (tol = 1e-9 sobre los vectores de Ritz): medido, hasta 3e-10
                # según el vector inicial aleatorio. Los valores propios salen
                # mucho mejor (error cuadrático en el del vector, ~1e-14).
                residual = K @ phi - M @ phi * lam
                self.assertLess(np.linalg.norm(residual) / np.linalg.norm(K @ phi), RESIDUAL_ACCEPT)

    def test_massless_dofs_es_cero_exacto(self):
        M = sp.diags([1.0, 0.0, 1e-300, 2.0])
        np.testing.assert_array_equal(massless_dofs(M), [1])

    def test_masa_no_semidefinida(self):
        K = sp.csr_matrix(np.diag([4.0, 5.0, 6.0]))
        M = sp.csr_matrix(np.array([[0.0, 0.5, 0.0], [0.5, 1.0, 0.0], [0.0, 0.0, 1.0]]))
        with self.assertRaises(SemidefiniteMassError) as ctx:
            EigenSolver().solve(K, M, 1)
        self.assertEqual(ctx.exception.cause, "indefinite_mass")

    def test_masa_singular_mas_alla_de_los_grados_de_libertad_sin_masa(self):
        K = sp.csr_matrix(np.diag([4.0, 5.0, 6.0]))
        M = sp.csr_matrix(np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 1.0], [0.0, 1.0, 1.0]]))
        with self.assertRaises(SemidefiniteMassError) as ctx:
            EigenSolver().solve(K, M, 1)
        self.assertEqual(ctx.exception.cause, "singular_mass")


if __name__ == "__main__":
    unittest.main()
