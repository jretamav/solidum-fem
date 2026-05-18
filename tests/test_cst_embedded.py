"""Tests de aceptación de ``CST_Embedded2D`` (ADR 0010, fase 2).

Cubre los criterios de ``acceptance`` declarados en
``docs/specs/CST_Embedded2D.md``: verificación física (estado intacto
bit-exact con Tri3, activación Rankine, descarga elástica del bulk
post-activación, simetría KOS de la condensada, fórmula `l_d` del Cap. 6,
recuperación local del salto), tests específicos (rechazo de bulks no
elásticos, validación de JUMP_DIM, no-chattering dentro del Newton,
identificación del nodo solitario, invariancia bajo `n → −n`) y tests
arquitecturales (registro, atributos declarativos, hook `prepare_step`).
"""
from __future__ import annotations

import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic
from fenix.core.discontinuity_state import DiscontinuityState
from fenix.core.element import Element
from fenix.core.node import Node
from fenix.elements.solid_2d.embedded_cst import (
    CST_Embedded2D,
    _compute_ld,
    _max_principal_stress_2d,
)
from fenix.elements.solid_2d.tri3 import Tri3
from fenix.materials.elastic_2d import Elastic2D
from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.registry import ElementRegistry


# --- Geometría y materiales canónicos para los tests ---------------------
def _triangle_nodes(coords):
    return [Node(i, list(coords[i])) for i in range(len(coords))]


def _ref_triangle_coords():
    return np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])


def _make_bulk(E: float = 30.0e9, nu: float = 0.2):
    return Elastic2D(E=E, nu=nu, hypothesis='plane_strain')


def _make_cohesive(sigma_t0: float = 2.5e6, G_f: float = 100.0,
                   K_e: float = 1.0e13, softening: str = 'linear'):
    return CohesiveDamageIsotropic(
        sigma_t0=sigma_t0, G_f=G_f, K_e=K_e, softening=softening,
    )


def _make_pair(coords=None, *, bulk=None, cohesive=None):
    """Construye un (Tri3, CST_Embedded2D) con la misma geometría/bulk."""
    if coords is None:
        coords = _ref_triangle_coords()
    if bulk is None:
        bulk = _make_bulk()
    if cohesive is None:
        cohesive = _make_cohesive()
    nodes_a = _triangle_nodes(coords)
    nodes_b = _triangle_nodes(coords)
    tri = Tri3(1, nodes_a, bulk)
    emb = CST_Embedded2D(2, nodes_b, bulk, cohesive)
    return tri, emb


def _assign_local_dofs(element):
    """Mapea los DOFs del elemento a índices secuenciales 0..ndof-1.

    Útil en tests aislados donde no hay ``Domain`` que orqueste el mapeo
    global. Después se puede invocar ``prepare_step(U)`` con ``U`` indexada
    por estos enteros.
    """
    idx = 0
    for node in element.nodes:
        for dof_name in element.DOF_NAMES:
            node.dofs[dof_name] = idx
            idx += 1
    return idx


# =========================================================================
# Verificación — bit-exact con Tri3 en estado intacto
# =========================================================================


class TestIntactReproducesTri3(unittest.TestCase):
    """acceptance.estado_intacto_reproduce_Tri3"""

    def test_zero_displacement(self):
        tri, emb = _make_pair()
        u = np.zeros(6)
        K_t, F_t = tri.compute_element_state(u)
        K_e, F_e = emb.compute_element_state(u)
        np.testing.assert_array_equal(K_t, K_e)
        np.testing.assert_array_equal(F_t, F_e)

    def test_arbitrary_displacement(self):
        tri, emb = _make_pair()
        rng = np.random.default_rng(42)
        u = 1e-4 * rng.standard_normal(6)
        K_t, F_t = tri.compute_element_state(u)
        K_e, F_e = emb.compute_element_state(u)
        np.testing.assert_array_equal(K_t, K_e)
        np.testing.assert_array_equal(F_t, F_e)

    def test_distorted_triangle(self):
        coords = np.array([[0.1, 0.2], [1.3, 0.05], [0.4, 1.5]])
        tri, emb = _make_pair(coords)
        rng = np.random.default_rng(7)
        u = 1e-3 * rng.standard_normal(6)
        K_t, F_t = tri.compute_element_state(u)
        K_e, F_e = emb.compute_element_state(u)
        np.testing.assert_array_equal(K_t, K_e)
        np.testing.assert_array_equal(F_t, F_e)


# =========================================================================
# Activación Rankine
# =========================================================================


class TestRankineActivation(unittest.TestCase):
    """acceptance.activacion_rankine_correcta"""

    def test_uniaxial_traction_activates_with_normal_x(self):
        # σ uniaxial en x sobre el CST de referencia (0,0)-(1,0)-(0,1).
        # Plane strain con E=30e9, ν=0.2 y εxx=1e-3 da σxx ≈ 33 MPa ≫ σ_t0=1e6.
        coords = _ref_triangle_coords()
        nodes = _triangle_nodes(coords)
        emb = CST_Embedded2D(1, nodes, _make_bulk(), _make_cohesive(sigma_t0=1.0e6))
        _assign_local_dofs(emb)

        U = np.zeros(6)
        U[2] = 1.0e-3                                # ux del nodo 1 = 1e-3
        emb.prepare_step(U)
        self.assertIsNotNone(emb.discontinuity_state)
        n = emb.discontinuity_state.normal
        self.assertAlmostEqual(abs(n[0]), 1.0, delta=1e-10)
        self.assertAlmostEqual(abs(n[1]), 0.0, delta=1e-10)

    def test_below_threshold_does_not_activate(self):
        coords = _ref_triangle_coords()
        nodes = _triangle_nodes(coords)
        emb = CST_Embedded2D(1, nodes, _make_bulk(), _make_cohesive(sigma_t0=1.0e12))
        _assign_local_dofs(emb)
        emb.prepare_step(np.zeros(6))
        self.assertIsNone(emb.discontinuity_state)

    def test_activation_is_irreversible(self):
        coords = _ref_triangle_coords()
        nodes = _triangle_nodes(coords)
        emb = CST_Embedded2D(1, nodes, _make_bulk(), _make_cohesive(sigma_t0=1.0e6))
        _assign_local_dofs(emb)

        U_high = np.zeros(6)
        U_high[2] = 1.0e-3
        emb.prepare_step(U_high)
        ds_first = emb.discontinuity_state
        self.assertIsNotNone(ds_first)

        # Descargar: U → 0; la discontinuidad persiste (irreversible)
        emb.prepare_step(np.zeros(6))
        self.assertIs(emb.discontinuity_state, ds_first)


# =========================================================================
# Descarga elástica del bulk post-activación
# =========================================================================


class TestBulkElasticUnloadingAfterActivation(unittest.TestCase):
    """acceptance.bulk_descarga_elasticamente_post_activacion"""

    def test_bulk_strain_uses_net_epsilon(self):
        coords = _ref_triangle_coords()
        nodes = _triangle_nodes(coords)
        bulk = _make_bulk()
        cohesive = _make_cohesive(sigma_t0=1.0e6)
        emb = CST_Embedded2D(1, nodes, bulk, cohesive)

        # Activar manualmente con n = (1, 0)
        emb._activate(np.array([1.0, 0.0]), coords)
        ds = emb.discontinuity_state
        self.assertIsNotNone(ds)

        # Fijar jump_trial sin pasar por Newton (test directo de la cinemática)
        ds.jump_trial = np.array([2.0e-5, 0.0])  # apertura normal pequeña
        ds.jump_committed = np.copy(ds.jump_trial)

        # Aplicar un u_e con ε_xx_total y verificar σ_bulk = C_e · (ε_total − ε_enriched)
        u_e = np.array([0.0, 0.0, 1.0e-3, 0.0, 0.0, 0.0])
        # Llamar compute_element_state — sin escribir el resultado, leer estados
        emb.compute_element_state(u_e)
        # σ trial = el de la última iteración. La ε neta = ε_total − B^φ·G·[[u]].
        # Verificamos numéricamente: σ_trial debe ser C_e · ε_neta.
        from fenix.elements.solid_2d._shared import _compute_kinematics_tri3
        B, _ = _compute_kinematics_tri3(coords)
        G = np.column_stack([ds.normal, ds.tangent])
        i_star = ds.solitary_node
        B_phi = B[:, 2 * i_star:2 * i_star + 2]
        # En el Newton local el jump cambia; usamos el jump_trial final
        eps_net = B @ u_e - B_phi @ (G @ ds.jump_trial)
        sigma_expected = bulk.C @ eps_net
        np.testing.assert_allclose(
            emb.state.stresses_trial[0], sigma_expected, rtol=1e-12, atol=1e-12,
        )


# =========================================================================
# Simetría de la condensada (KOS)
# =========================================================================


class TestCondensedSymmetry(unittest.TestCase):
    """acceptance.condensacion_simetrica_KOS"""

    def test_symmetric_in_loading(self):
        coords = _ref_triangle_coords()
        nodes = _triangle_nodes(coords)
        cohesive = _make_cohesive(sigma_t0=1.0e6)
        emb = CST_Embedded2D(1, nodes, _make_bulk(), cohesive)
        emb._activate(np.array([1.0, 0.0]), coords)

        u_e = np.array([0.0, 0.0, 5.0e-4, 0.0, 0.0, 0.0])
        K_cond, _ = emb.compute_element_state(u_e)
        frob = np.linalg.norm(K_cond - K_cond.T, 'fro')
        ref = max(np.linalg.norm(K_cond, 'fro'), 1.0)
        self.assertLess(frob / ref, 1e-12)


# =========================================================================
# Fórmula cerrada de l_d (Cap. 6 Retama 2010)
# =========================================================================


class TestLdFormula(unittest.TestCase):
    """acceptance.ld_formula_cerrada"""

    def test_axis_aligned_normal_perpendicular_to_opposite_edge(self):
        # Triángulo (0,0)-(1,0)-(0,1); nodo solitario = 0; lado opuesto entre 1 y 2
        # de longitud √2 con α = 135°. Si Γ_d tiene tangente s a lo largo del
        # eje y (θ = 90°), entonces |cos(θ − α)| = |cos(−45°)| = √2/2.
        coords = _ref_triangle_coords()
        tangent = np.array([0.0, 1.0])
        l_d = _compute_ld(coords, tangent, solitary=0)
        A_e = 0.5
        h = math.sqrt(2.0) / 2.0  # altura desde (0,0) al lado entre (1,0)-(0,1)
        # Verificar h con 2A/L_opp
        L_opp = math.sqrt(2.0)
        h_check = 2.0 * A_e / L_opp
        self.assertAlmostEqual(h, h_check, places=14)
        # cos(θ − α) = cos(π/2 − 3π/4) = cos(−π/4) = √2/2
        expected = (A_e / h) * math.sqrt(2.0) / 2.0
        self.assertAlmostEqual(l_d, expected, places=14)

    def test_invariant_under_n_flip(self):
        # acceptance.ld_invariante_bajo_reflejo_de_n
        coords = np.array([[0.0, 0.0], [2.0, 0.0], [0.5, 1.5]])
        tangent = np.array([1.0, 0.5])
        tangent = tangent / np.linalg.norm(tangent)
        l_a = _compute_ld(coords, tangent, solitary=2)
        l_b = _compute_ld(coords, -tangent, solitary=2)
        self.assertAlmostEqual(l_a, l_b, places=14)


# =========================================================================
# Recuperación local del salto tras commit
# =========================================================================


class TestLocalRecovery(unittest.TestCase):
    """acceptance.recuperacion_local_del_salto

    Tras la convergencia del Newton local, R^{[[u]]} debe ser despreciable.
    """

    def test_residual_vanishes_after_loading(self):
        coords = _ref_triangle_coords()
        nodes = _triangle_nodes(coords)
        cohesive = _make_cohesive(sigma_t0=1.0e6, softening='linear')
        emb = CST_Embedded2D(1, nodes, _make_bulk(), cohesive)
        emb._activate(np.array([1.0, 0.0]), coords)

        u_e = np.array([0.0, 0.0, 2.0e-4, 0.0, 0.0, 0.0])
        emb.compute_element_state(u_e)
        ds = emb.discontinuity_state

        # Reconstruir R^{[[u]]} con el jump_trial obtenido
        from fenix.elements.solid_2d._shared import _compute_kinematics_tri3
        B, detJ = _compute_kinematics_tri3(coords)
        G = np.column_stack([ds.normal, ds.tangent])
        i_star = ds.solitary_node
        B_phi = B[:, 2 * i_star:2 * i_star + 2]
        vol = 0.5 * detJ * emb.thickness

        strain_bulk = B @ u_e - B_phi @ (G @ ds.jump_trial)
        sigma = emb.material.C @ strain_bulk
        t_local, _, _ = emb.cohesive_material.compute_traction(
            ds.jump_trial, ds.cohesive_state_committed or None,
        )
        R_jump = -G.T @ (B_phi.T @ sigma) * vol + ds.l_d * t_local
        self.assertLess(np.linalg.norm(R_jump), 1e-6)


# =========================================================================
# Validación de bulks y cohesivos
# =========================================================================


class TestConstructorValidations(unittest.TestCase):
    """acceptance.bulk_no_lineal_rechazado + jump_dim_validado"""

    def test_bulk_nonlinear_rejected(self):
        coords = _ref_triangle_coords()
        nodes = _triangle_nodes(coords)
        bulk_nl = IsotropicDamage2D(E=30.0e9, nu=0.2, kappa_0=1.0e-4, alpha=100.0)
        cohesive = _make_cohesive()
        with self.assertRaises(ValueError) as cm:
            CST_Embedded2D(1, nodes, bulk_nl, cohesive)
        self.assertIn('IsotropicDamage2D', str(cm.exception))

    def test_cohesive_must_be_CohesiveMaterial(self):
        coords = _ref_triangle_coords()
        nodes = _triangle_nodes(coords)
        bulk = _make_bulk()
        # Pasar un bulk como cohesivo → TypeError
        with self.assertRaises(TypeError):
            CST_Embedded2D(1, nodes, bulk, bulk)

    def test_jump_dim_mismatch_rejected(self):
        coords = _ref_triangle_coords()
        nodes = _triangle_nodes(coords)

        # Cohesivo falso con JUMP_DIM=3 (subclase mínima para el test)
        from fenix.core.cohesive_material import CohesiveMaterial

        class Fake3D(CohesiveMaterial):
            JUMP_DIM = 3
            def compute_traction(self, jump, state_vars=None):
                return np.zeros(3), np.zeros((3, 3)), {}

        with self.assertRaises(ValueError) as cm:
            CST_Embedded2D(1, nodes, _make_bulk(), Fake3D())
        self.assertIn('JUMP_DIM', str(cm.exception))


# =========================================================================
# No-chattering: compute_element_state no activa
# =========================================================================


class TestNoChatteringWithinNewton(unittest.TestCase):
    """acceptance.activacion_no_chattering_dentro_de_newton"""

    def test_compute_element_state_does_not_activate(self):
        coords = _ref_triangle_coords()
        nodes = _triangle_nodes(coords)
        cohesive = _make_cohesive(sigma_t0=1.0e6)
        emb = CST_Embedded2D(1, nodes, _make_bulk(), cohesive)
        # u_e que produciría σ > σ_t0 si se chequease
        u_e = np.array([0.0, 0.0, 1.0e-3, 0.0, 0.0, 0.0])
        for _ in range(5):
            emb.compute_element_state(u_e)
        self.assertIsNone(emb.discontinuity_state)


# =========================================================================
# Identificación del nodo solitario
# =========================================================================


class TestSolitaryNodeIdentification(unittest.TestCase):
    """acceptance.identificacion_nodo_solitario"""

    def test_identifies_unique_node_on_positive_side(self):
        # Triángulo (0,0)-(1,0)-(0.5, 0.866); centroide = (0.5, 0.289).
        # Con n=(0, 1), el nodo 2 (y=0.866) está claramente del lado positivo;
        # los nodos 0 y 1 (y=0) están del lado negativo.
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, math.sqrt(3) / 2.0]])
        nodes = _triangle_nodes(coords)
        emb = CST_Embedded2D(1, nodes, _make_bulk(), _make_cohesive())
        emb._activate(np.array([0.0, 1.0]), coords)
        self.assertEqual(emb.discontinuity_state.solitary_node, 2)


# =========================================================================
# Helpers — máxima principal en 2D
# =========================================================================


class TestPrincipalStressHelper(unittest.TestCase):
    def test_uniaxial_x(self):
        sigma = np.array([1.0e6, 0.0, 0.0])
        s_I, n = _max_principal_stress_2d(sigma)
        self.assertAlmostEqual(s_I, 1.0e6, places=10)
        self.assertAlmostEqual(abs(n[0]), 1.0, places=10)
        self.assertAlmostEqual(abs(n[1]), 0.0, places=10)

    def test_uniaxial_y(self):
        sigma = np.array([0.0, 1.0e6, 0.0])
        s_I, n = _max_principal_stress_2d(sigma)
        self.assertAlmostEqual(s_I, 1.0e6, places=10)
        self.assertAlmostEqual(abs(n[0]), 0.0, places=10)
        self.assertAlmostEqual(abs(n[1]), 1.0, places=10)

    def test_pure_shear(self):
        sigma = np.array([0.0, 0.0, 1.0e6])
        s_I, n = _max_principal_stress_2d(sigma)
        self.assertAlmostEqual(s_I, 1.0e6, places=10)
        # Dirección principal a 45°
        self.assertAlmostEqual(abs(n[0]), math.sqrt(2.0) / 2.0, places=10)
        self.assertAlmostEqual(abs(n[1]), math.sqrt(2.0) / 2.0, places=10)


# =========================================================================
# Arch: registro, atributos y prepare_step en base
# =========================================================================


class TestArchitecture(unittest.TestCase):
    def test_registered_in_element_registry(self):
        self.assertIn('CST_Embedded2D', ElementRegistry._items)
        self.assertIs(ElementRegistry._items['CST_Embedded2D'], CST_Embedded2D)

    def test_declarative_contract(self):
        self.assertEqual(CST_Embedded2D.DOF_NAMES, ['ux', 'uy'])
        self.assertEqual(CST_Embedded2D.STRAIN_DIM, 3)
        self.assertEqual(CST_Embedded2D.N_INTEGRATION_POINTS, 1)

    def test_element_base_prepare_step_is_noop(self):
        # Comprobamos que la base provee prepare_step (no-op) sin lanzar.
        self.assertTrue(hasattr(Element, 'prepare_step'))
        # Para un Tri3 (que no sobreescribe), debe ser invocable sin efecto.
        coords = _ref_triangle_coords()
        tri = Tri3(1, _triangle_nodes(coords), _make_bulk())
        tri.prepare_step(np.zeros(6))  # no debería lanzar nada


# =========================================================================
# Integración end-to-end: parser YAML construye el elemento correctamente
# =========================================================================


class TestYamlIntegration(unittest.TestCase):
    """Verifica que el parser YAML construye CST_Embedded2D con bulk y cohesivo."""

    def test_yaml_round_trip(self):
        import tempfile

        from fenix.utils.yaml_parser import YamlParser

        yaml_text = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [0.0, 1.0]}

materials:
  - {id: 1, type: Elastic2D, E: 30.0e9, nu: 0.2, hypothesis: plane_strain}

cohesive_materials:
  - {id: 1, type: CohesiveDamageIsotropic, sigma_t0: 2.5e6, G_f: 100.0,
     K_e: 1.0e13, softening: linear}

elements:
  - {id: 1, type: CST_Embedded2D, nodes: [1, 2, 3],
     material: 1, cohesive_material: 1}
"""
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.yaml', delete=False, encoding='utf-8',
        ) as f:
            f.write(yaml_text)
            path = f.name

        try:
            parser = YamlParser(path)
            domain = parser.parse()
            self.assertEqual(len(domain.elements), 1)
            elem = next(iter(domain.elements.values()))
            self.assertIsInstance(elem, CST_Embedded2D)
            # Bulk y cohesive_material correctamente inyectados
            self.assertEqual(type(elem.material).__name__, 'Elastic2D')
            self.assertEqual(type(elem.cohesive_material).__name__, 'CohesiveDamageIsotropic')
            self.assertIsNone(elem.discontinuity_state)  # intacto al construir
        finally:
            os.unlink(path)

    def test_yaml_rejects_unknown_cohesive_reference(self):
        import tempfile

        from fenix.utils.yaml_parser import YamlParser, YamlValidationError

        yaml_text = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [0.0, 1.0]}

materials:
  - {id: 1, type: Elastic2D, E: 30.0e9, nu: 0.2, hypothesis: plane_strain}

cohesive_materials:
  - {id: 1, type: CohesiveDamageIsotropic, sigma_t0: 2.5e6, G_f: 100.0,
     K_e: 1.0e13, softening: linear}

elements:
  - {id: 1, type: CST_Embedded2D, nodes: [1, 2, 3],
     material: 1, cohesive_material: 999}
"""
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.yaml', delete=False, encoding='utf-8',
        ) as f:
            f.write(yaml_text)
            path = f.name

        try:
            with self.assertRaises(YamlValidationError) as cm:
                YamlParser(path).parse()
            self.assertIn('cohesive_material inexistente', str(cm.exception))
        finally:
            os.unlink(path)


if __name__ == '__main__':
    unittest.main()
