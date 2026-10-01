"""Tests de los elementos de cable.

Materializan los criterios de `acceptance` declarados en
docs/specs/Cable2DCorot.md.
"""
import math
import os
import tempfile
import textwrap
import unittest

import numpy as np
from scipy.optimize import brentq

import solidum  # dispara autodiscover
from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.elements.cable import Cable2DCorot, Cable3DCorot
from solidum.elements.truss import Truss2D, Truss3D, Truss2DCorot, Truss3DCorot
from solidum.materials.cable_1d import CableMaterial1D
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver
from solidum.utils.yaml_parser import YamlParser


class TestCable2DCorot(unittest.TestCase):

    def _build(self, coords1, coords2, material, A=1.0):
        n1 = Node(1, list(coords1))
        n2 = Node(2, list(coords2))
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy')
            node.dofs['ux'] = 2 * k
            node.dofs['uy'] = 2 * k + 1
        return Cable2DCorot(1, [n1, n2], material, A)

    def test_acceptance_cable_tensado_como_barra_corot(self):
        """Criterio 1: tensado con material lineal coincide con barra corot."""
        E, A = 1000.0, 1.0
        mat = Elastic1D(E=E)
        cable = self._build([0.0, 0.0], [1.0, 0.0], mat, A=A)
        u_e = np.array([0.0, 0.0, 1e-3, 0.0])  # ε = 1e-3 > 0 → tensado

        K_T, F_int = cable.compute_element_state(u_e)

        # Valores esperados analíticos
        l = 1.0 + 1e-3
        N_expected = E * 1e-3 * A  # σ·A
        d = np.array([-1.0, 0.0, 1.0, 0.0])
        n = np.array([0.0, 1.0, 0.0, -1.0])
        K_M_expected = (E * A / 1.0) * np.outer(d, d)
        K_G_expected = (N_expected / l) * np.outer(n, n)

        self.assertTrue(np.allclose(K_T, K_M_expected + K_G_expected, atol=1e-10))
        self.assertTrue(np.allclose(F_int, N_expected * d, atol=1e-10))

    def test_acceptance_cable_destensado_aportacion_nula(self):
        """Criterio 2: con material unilateral + ε < 0 ⇒ K_T = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0], [1.0, 0.0], mat, A=1.0)
        u_e = np.array([0.0, 0.0, -1e-3, 0.0])  # ε < 0 → destensado

        K_T, F_int = cable.compute_element_state(u_e)

        self.assertTrue(np.allclose(K_T, 0.0, atol=1e-12))
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-12))
        self.assertEqual(cable.state.stresses_trial[0], 0.0)

    def test_acceptance_rotacion_rigida(self):
        """Criterio 3: rotación rígida ⇒ ε = 0 ⇒ σ = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0], [1.0, 0.0], mat, A=1.0)
        theta = math.pi / 4
        c, s = math.cos(theta), math.sin(theta)
        u_e = np.array([0.0, 0.0, c - 1.0, s])  # nodo 2 en (c, s)

        _, F_int = cable.compute_element_state(u_e)

        self.assertAlmostEqual(cable.state.stresses_trial[0], 0.0, places=12)
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-10))

    def test_acceptance_cruce_por_cero(self):
        """Criterio 4: ε = 0 exacto ⇒ σ = 0, E_t = 0, K_T = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0], [1.0, 0.0], mat, A=1.0)
        u_e = np.zeros(4)  # ε = 0

        K_T, F_int = cable.compute_element_state(u_e)

        self.assertEqual(cable.state.stresses_trial[0], 0.0)
        self.assertTrue(np.allclose(K_T, 0.0, atol=1e-12))
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-12))

    def test_registro_en_registry(self):
        """El elemento queda registrado vía autodiscover."""
        from solidum.registry import ElementRegistry
        self.assertIn('Cable2DCorot', ElementRegistry._items)


class TestCable3DCorot(unittest.TestCase):

    def _build(self, coords1, coords2, material, A=1.0, cls=Cable3DCorot):
        n1 = Node(1, list(coords1))
        n2 = Node(2, list(coords2))
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy'); node.add_dof('uz')
            node.dofs['ux'] = 3 * k
            node.dofs['uy'] = 3 * k + 1
            node.dofs['uz'] = 3 * k + 2
        return cls(1, [n1, n2], material, A)

    def test_acceptance_cable_tensado_como_barra_corot_3d(self):
        """Criterio 1: tensado con material lineal coincide con Truss3DCorot."""
        E, A = 1000.0, 1.0
        cable = self._build([0.0, 0.0, 0.0], [3.0, 4.0, 12.0], Elastic1D(E=E), A=A)
        truss = self._build([0.0, 0.0, 0.0], [3.0, 4.0, 12.0], Elastic1D(E=E), A=A,
                            cls=Truss3DCorot)
        axis = np.array([3.0, 4.0, 12.0]) / 13.0
        u_tiny = np.concatenate([np.zeros(3), 1e-6 * axis])

        K_c, F_c = cable.compute_element_state(u_tiny)
        K_t, F_t = truss.compute_element_state(u_tiny)

        self.assertTrue(np.allclose(K_c, K_t, rtol=1e-4))
        self.assertTrue(np.allclose(F_c, F_t, rtol=1e-4))

    def test_acceptance_cable_destensado_aportacion_nula(self):
        """Criterio 2: con material unilateral + ε < 0 ⇒ K_T = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], mat, A=1.0)
        u_e = np.array([0.0, 0.0, 0.0, -1e-3, 0.0, 0.0])  # ε < 0

        K_T, F_int = cable.compute_element_state(u_e)

        self.assertTrue(np.allclose(K_T, 0.0, atol=1e-12))
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-12))
        self.assertEqual(cable.state.stresses_trial[0], 0.0)

    def test_acceptance_rotacion_rigida_3d(self):
        """Criterio 3: rotación rígida 3D ⇒ ε = 0 ⇒ σ = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], mat, A=1.0)
        theta = math.pi / 4
        c, s = math.cos(theta), math.sin(theta)
        # R_y(theta) @ (1,0,0) = (c, 0, -s); u_2 = final - inicial
        u_e = np.array([0.0, 0.0, 0.0, c - 1.0, 0.0, -s])

        _, F_int = cable.compute_element_state(u_e)

        self.assertAlmostEqual(cable.state.stresses_trial[0], 0.0, places=12)
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-10))

    def test_acceptance_cruce_por_cero(self):
        """Criterio 4: u_e = 0 ⇒ σ = 0, E_t = 0, K_T = 0, F_int = 0."""
        mat = CableMaterial1D(E=1000.0)
        cable = self._build([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], mat, A=1.0)
        u_e = np.zeros(6)

        K_T, F_int = cable.compute_element_state(u_e)

        self.assertEqual(cable.state.stresses_trial[0], 0.0)
        self.assertTrue(np.allclose(K_T, 0.0, atol=1e-12))
        self.assertTrue(np.allclose(F_int, 0.0, atol=1e-12))

    def test_registro_en_registry(self):
        from solidum.registry import ElementRegistry
        self.assertIn('Cable3DCorot', ElementRegistry._items)


class TestUnilateralMaterialValidation(unittest.TestCase):
    """Validación declarativa material.IS_UNILATERAL vs element.ACCEPTS_UNILATERAL.

    Un material cuya rigidez tangente puede colapsar a cero (cables, gaps)
    solo es admisible en elementos preparados para manejarlo. La base
    Element atrapa la combinación incompatible en construcción.
    """

    def _make_nodes_2d(self):
        n1 = Node(1, [0.0, 0.0]); n2 = Node(2, [1.0, 0.0])
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy')
            node.dofs['ux'] = 2 * k
            node.dofs['uy'] = 2 * k + 1
        return n1, n2

    def _make_nodes_3d(self):
        n1 = Node(1, [0.0, 0.0, 0.0]); n2 = Node(2, [1.0, 0.0, 0.0])
        for k, node in enumerate((n1, n2)):
            for d in ('ux', 'uy', 'uz'):
                node.add_dof(d)
            node.dofs['ux'] = 3 * k
            node.dofs['uy'] = 3 * k + 1
            node.dofs['uz'] = 3 * k + 2
        return n1, n2

    def test_truss2d_rechaza_cable_material(self):
        n1, n2 = self._make_nodes_2d()
        mat = CableMaterial1D(E=1000.0)
        with self.assertRaisesRegex(ValueError, "unilateral"):
            Truss2D(1, [n1, n2], mat, A=1.0)

    def test_truss3d_rechaza_cable_material(self):
        n1, n2 = self._make_nodes_3d()
        mat = CableMaterial1D(E=1000.0)
        with self.assertRaisesRegex(ValueError, "unilateral"):
            Truss3D(1, [n1, n2], mat, A=1.0)

    def test_truss_corot_tambien_rechaza_cable_material(self):
        """Truss*Corot tampoco lo acepta — el material unilateral se reserva
        a Cable*Corot, donde la formulación está diseñada explícitamente
        para él."""
        n1, n2 = self._make_nodes_2d()
        with self.assertRaisesRegex(ValueError, "unilateral"):
            Truss2DCorot(1, [n1, n2], CableMaterial1D(E=1000.0), A=1.0)
        n3, n4 = self._make_nodes_3d()
        with self.assertRaisesRegex(ValueError, "unilateral"):
            Truss3DCorot(1, [n3, n4], CableMaterial1D(E=1000.0), A=1.0)

    def test_cable_corot_acepta_cable_material(self):
        """Combinación canónica — no debe lanzar."""
        n1, n2 = self._make_nodes_2d()
        Cable2DCorot(1, [n1, n2], CableMaterial1D(E=1000.0), A=1.0)
        n3, n4 = self._make_nodes_3d()
        Cable3DCorot(2, [n3, n4], CableMaterial1D(E=1000.0), A=1.0)

    def test_truss_sigue_aceptando_materiales_bilaterales(self):
        """Regresión: la nueva validación no rompe el contrato con Elastic1D."""
        n1, n2 = self._make_nodes_2d()
        Truss2D(1, [n1, n2], Elastic1D(E=1000.0), A=1.0)
        n3, n4 = self._make_nodes_3d()
        Truss3D(1, [n3, n4], Elastic1D(E=1000.0), A=1.0)


class TestPretensionLongitudSinTension(unittest.TestCase):
    """Pretensión por longitud sin tensión ``L0`` (criterios 5-10 de
    docs/specs/Cable2DCorot.md y Cable3DCorot.md, ampliación del 2026-10-01).

    Coeficientes no unitarios a propósito (cable de acero de 3 cm², vano de
    4 m, 25 kN de pretensión): un error dimensional en ``L0`` no pasa
    desapercibido.
    """

    E = 2.0e11      # Pa
    A = 3.0e-4      # m²
    L = 4.0         # m, vano completo
    T0 = 25.0e3     # N, pretensión

    @classmethod
    def _L0_para(cls, L_geo, T0):
        """Longitud sin tensión que produce T0 en la geometría L_geo:
        N = E·A·(L − L₀)/L₀ ⇒ L₀ = L/(1 + T₀/(E·A))."""
        return L_geo / (1.0 + T0 / (cls.E * cls.A))

    def _cable2d(self, coords1, coords2, L0=None, material=None):
        n1, n2 = Node(1, list(coords1)), Node(2, list(coords2))
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy')
            node.dofs['ux'] = 2 * k
            node.dofs['uy'] = 2 * k + 1
        mat = material if material is not None else CableMaterial1D(E=self.E)
        return Cable2DCorot(1, [n1, n2], mat, self.A, L0=L0)

    def _cable3d(self, coords1, coords2, L0=None):
        n1, n2 = Node(1, list(coords1)), Node(2, list(coords2))
        for k, node in enumerate((n1, n2)):
            node.add_dof('ux'); node.add_dof('uy'); node.add_dof('uz')
            node.dofs['ux'] = 3 * k
            node.dofs['uy'] = 3 * k + 1
            node.dofs['uz'] = 3 * k + 2
        return Cable3DCorot(1, [n1, n2], CableMaterial1D(E=self.E), self.A, L0=L0)

    def _w_analitica(self, P, L0_mitad):
        """Flecha central exacta de la cuerda pretensada de dos tramos:
        2·N(l)·w/l = P, con l = √(a² + w²) y N = E·A·(l − L₀)/L₀."""
        a = self.L / 2.0

        def residuo(w):
            l = math.hypot(a, w)
            N = self.E * self.A * (l - L0_mitad) / L0_mitad
            return 2.0 * N * w / l - P

        w = brentq(residuo, 0.0, a, xtol=1e-16, rtol=1e-15, maxiter=200)
        l = math.hypot(a, w)
        return w, self.E * self.A * (l - L0_mitad) / L0_mitad

    def _cuerda(self, P, tres_d=False):
        """Cuerda pretensada de dos tramos con carga P en el nodo central
        (−y en 2D, −z en 3D); extremos articulados fijos."""
        dom = Domain()
        a = self.L / 2.0
        if tres_d:
            coords = ([0.0, 0.0, 0.0], [a, 0.0, 0.0], [self.L, 0.0, 0.0])
            cls, dofs, eje = Cable3DCorot, ('ux', 'uy', 'uz'), 'uz'
        else:
            coords = ([0.0, 0.0], [a, 0.0], [self.L, 0.0])
            cls, dofs, eje = Cable2DCorot, ('ux', 'uy'), 'uy'
        n1, n2, n3 = (dom.add_node(i + 1, c) for i, c in enumerate(coords))
        L0 = self._L0_para(a, self.T0)
        c1 = cls(1, [n1, n2], CableMaterial1D(E=self.E), self.A, L0=L0)
        c2 = cls(2, [n2, n3], CableMaterial1D(E=self.E), self.A, L0=L0)
        for nodo in (n1, n3):
            for dof in dofs:
                nodo.fix_dof(dof, 0.0)
        dom.add_element(c1)
        dom.add_element(c2)
        dom.generate_equation_numbers(verbose=False)
        F = np.zeros(dom.total_dofs)
        F[n2.dofs[eje]] = -P
        conv = ConvergenceCriterion(rtol_force=1e-12, rtol_disp=1e-12)
        U = NonlinearSolver(Assembler(dom), convergence=conv, num_steps=10).solve(F)
        return U, n2, c1, L0, eje

    def test_pretension_inicial_exacta_2d(self):
        """Criterio 5: en la geometría inicial N = T₀ exacta y K_G = (T₀/L)·n⊗n."""
        L0 = self._L0_para(self.L, self.T0)
        cable = self._cable2d([0.0, 0.0], [self.L, 0.0], L0=L0)
        K, F_int = cable.compute_element_state(np.zeros(4))
        N = cable.state.stresses_trial[0] * self.A
        self.assertAlmostEqual(N / self.T0, 1.0, places=12)

        d = np.array([-1.0, 0.0, 1.0, 0.0])
        n = np.array([0.0, 1.0, 0.0, -1.0])
        K_esp = ((self.E * self.A / L0) * np.outer(d, d)
                 + (self.T0 / self.L) * np.outer(n, n))
        np.testing.assert_allclose(K, K_esp, rtol=1e-12, atol=1e-12 * np.abs(K_esp).max())
        np.testing.assert_allclose(F_int, self.T0 * d, rtol=1e-12)

        fuerzas = cable.internal_forces(np.zeros(4))
        np.testing.assert_allclose(fuerzas.components['N'], [self.T0, self.T0], rtol=1e-12)

    def test_pretension_inicial_exacta_3d_oblicuo(self):
        """Criterio 5 en 3D, con un cable oblicuo: N = T₀ y rigidez transversal T₀/L."""
        p1, p2 = np.array([0.5, -1.0, 2.0]), np.array([2.5, 1.5, 4.2])
        L_geo = float(np.linalg.norm(p2 - p1))
        L0 = self._L0_para(L_geo, self.T0)
        cable = self._cable3d(p1, p2, L0=L0)
        K, F_int = cable.compute_element_state(np.zeros(6))
        self.assertAlmostEqual(cable.state.stresses_trial[0] * self.A / self.T0, 1.0, places=12)

        e = (p2 - p1) / L_geo
        t = np.cross(e, [0.0, 0.0, 1.0])
        t /= np.linalg.norm(t)
        v = np.concatenate([np.zeros(3), t])   # desplazamiento transversal del nodo j
        self.assertAlmostEqual(v @ K @ v / (self.T0 / L_geo), 1.0, places=10)
        np.testing.assert_allclose(F_int, self.T0 * np.concatenate([-e, e]), rtol=1e-12)

    def test_cuerda_pretensada_carga_central_2d(self):
        """Criterio 6: cuerda pretensada con carga central frente a la solución
        analítica no lineal del equilibrio."""
        P = 3.0e3
        U, n2, c1, L0, eje = self._cuerda(P)
        w, N = self._w_analitica(P, L0)
        self.assertAlmostEqual(-U[n2.dofs[eje]] / w, 1.0, places=9)
        self.assertLess(abs(U[n2.dofs['ux']]), 1e-12 * w)
        N_fem = c1.internal_forces(U).components['N']
        np.testing.assert_allclose(N_fem, [N, N], rtol=1e-9)

    def test_cuerda_pretensada_limite_carga_pequena(self):
        """Criterio 6 (límite lineal): con P ≪ T₀ la flecha tiende a la de la
        cuerda tensa clásica, w = P·L/(4·T₀)."""
        P = 10.0
        U, n2, _, _, eje = self._cuerda(P)
        w_cuerda = P * self.L / (4.0 * self.T0)
        self.assertLess(abs(-U[n2.dofs[eje]] / w_cuerda - 1.0), 1e-3)

    def test_cuerda_pretensada_carga_central_3d(self):
        """Criterio 6 en 3D: misma solución analítica con la carga según −z."""
        P = 3.0e3
        U, n2, c1, L0, eje = self._cuerda(P, tres_d=True)
        w, N = self._w_analitica(P, L0)
        self.assertAlmostEqual(-U[n2.dofs[eje]] / w, 1.0, places=9)
        self.assertLess(abs(U[n2.dofs['uy']]), 1e-12 * w)
        N_fem = c1.internal_forces(U).components['N']
        np.testing.assert_allclose(N_fem, [N, N], rtol=1e-9)

    def test_cable_flojo_si_L0_mayor_que_la_geometrica(self):
        """Criterio 7: con L₀ > L y material unilateral el cable arranca flojo."""
        cable = self._cable2d([0.0, 0.0], [self.L, 0.0], L0=1.01 * self.L)
        K, F_int = cable.compute_element_state(np.zeros(4))
        self.assertTrue(np.all(K == 0.0))
        self.assertTrue(np.all(F_int == 0.0))
        self.assertEqual(cable.state.stresses_trial[0], 0.0)

    def test_L0_por_defecto_es_la_longitud_geometrica(self):
        """Criterio 8: sin L0 el elemento es idéntico al de antes de la ampliación."""
        p1, p2 = [0.3, -0.2], [3.1, 1.9]
        L_geo = math.dist(p1, p2)
        por_defecto = self._cable2d(p1, p2)
        explicito = self._cable2d(p1, p2, L0=L_geo)
        self.assertEqual(por_defecto.L0, L_geo)
        u_e = np.array([1.2e-3, -0.7e-3, 2.9e-3, 0.4e-3])
        K1, F1 = por_defecto.compute_element_state(u_e)
        K2, F2 = explicito.compute_element_state(u_e)
        np.testing.assert_array_equal(K1, K2)
        np.testing.assert_array_equal(F1, F2)

    def test_L0_invalido_se_rechaza(self):
        """Criterio 9: L0 debe ser positiva y finita."""
        for malo in (0.0, -1.0, math.inf, math.nan):
            with self.subTest(L0=malo):
                with self.assertRaisesRegex(ValueError, "L0"):
                    self._cable2d([0.0, 0.0], [self.L, 0.0], L0=malo)
                with self.assertRaisesRegex(ValueError, "L0"):
                    self._cable3d([0.0, 0.0, 0.0], [self.L, 0.0, 0.0], L0=malo)

    def test_masa_y_carga_de_cuerpo_con_longitud_sin_tension(self):
        """Criterio 10: la masa y la carga de cuerpo son las del material del
        cable descargado, ρ·A·L₀."""
        rho, g = 7850.0, 9.81
        L0 = self._L0_para(self.L, self.T0)
        cable = self._cable2d([0.0, 0.0], [self.L, 0.0], L0=L0,
                              material=CableMaterial1D(E=self.E, density=rho))
        m = rho * self.A * L0
        np.testing.assert_allclose(np.diag(cable.compute_mass_matrix("lumped")),
                                   [m / 2.0] * 4, rtol=1e-14)
        M = cable.compute_mass_matrix("consistent")
        self.assertAlmostEqual(M[0::2, 0::2].sum() / m, 1.0, places=13)
        f = cable.compute_body_load(np.array([0.0, -rho * g]))
        np.testing.assert_allclose(f, [0.0, -m * g / 2.0, 0.0, -m * g / 2.0], rtol=1e-14)

    def test_L0_desde_yaml(self):
        """El parámetro se declara en el bloque del elemento del YAML."""
        L0 = self._L0_para(self.L, self.T0)
        texto = textwrap.dedent(f"""
            nodes:
              - {{id: 1, coords: [0.0, 0.0]}}
              - {{id: 2, coords: [{self.L}, 0.0]}}
            materials:
              - {{id: 1, type: CableMaterial1D, E: {self.E!r}}}
            elements:
              - {{id: 1, type: Cable2DCorot, material: 1, nodes: [1, 2], A: {self.A!r}, L0: {L0!r}}}
            boundary_conditions_by_node:
              - {{node_id: 1, ux: 0.0, uy: 0.0}}
              - {{node_id: 2, ux: 0.0, uy: 0.0}}
            solver:
              type: NonlinearSolver
        """)
        with tempfile.TemporaryDirectory() as carpeta:
            ruta = os.path.join(carpeta, "cable.yaml")
            with open(ruta, "w", encoding="utf-8") as fh:
                fh.write(texto)
            dom = YamlParser(ruta).parse()
        self.assertEqual(dom.elements[1].L0, L0)


if __name__ == '__main__':
    unittest.main()
