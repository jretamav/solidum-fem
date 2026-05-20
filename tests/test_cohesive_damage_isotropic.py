"""Tests de aceptación de ``CohesiveDamageIsotropic`` (ADR 0010, fase 1).

Cubre los criterios de ``acceptance`` declarados en
``docs/specs/CohesiveDamageIsotropic.md``: verificación física (respuesta
elástica, inicio del daño, integral ``∫t·d[[u_n]] = G_F``, descarga,
recarga, saturación), tests específicos (simetría de la tangente,
consistencia por diferencias finitas, comportamiento en tangencial y
compresión, degeneración a elástico) y tests arquitecturales (registro,
``IS_SYMMETRIC``).
"""
from __future__ import annotations

import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.constants import DAMAGE_MAX
from solidum.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic
from solidum.registry import CohesiveMaterialRegistry


# --- Parámetros tipo hormigón para los casos de validación ---------------
SIGMA_T0 = 2.57e6      # Pa
G_F = 121.9            # N/m
# K_e suficientemente grande para que kappa_0 = sigma_t0/K_e ≪ w_c.
# Con K_e = 1e13: kappa_0 = 2.57e-7 m; w_c = 9.49e-5 m → ratio ~370.
K_E = 1.0e13


def make_linear(K_e: float = K_E) -> CohesiveDamageIsotropic:
    return CohesiveDamageIsotropic(
        sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_e, softening='linear'
    )


def make_exponential(K_e: float = K_E) -> CohesiveDamageIsotropic:
    return CohesiveDamageIsotropic(
        sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_e, softening='exponential'
    )


# =========================================================================
# Construcción y validación de parámetros
# =========================================================================


class TestConstructorValidations(unittest.TestCase):
    def test_rejects_nonpositive_sigma_t0(self):
        with self.assertRaises(ValueError):
            CohesiveDamageIsotropic(sigma_t0=0.0, G_f=G_F, K_e=K_E, softening='linear')

    def test_rejects_nonpositive_G_f(self):
        with self.assertRaises(ValueError):
            CohesiveDamageIsotropic(sigma_t0=SIGMA_T0, G_f=-1.0, K_e=K_E, softening='linear')

    def test_rejects_nonpositive_K_e(self):
        with self.assertRaises(ValueError):
            CohesiveDamageIsotropic(sigma_t0=SIGMA_T0, G_f=G_F, K_e=0.0, softening='linear')

    def test_rejects_unknown_softening(self):
        with self.assertRaises(ValueError):
            CohesiveDamageIsotropic(sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_E, softening='cubic')

    def test_rejects_exponential_with_nonpositive_H(self):
        # Forzar H = G_f - sigma_t0^2/(2 K_e) <= 0 con K_e muy pequeño
        K_low = SIGMA_T0 * SIGMA_T0 / (2.0 * G_F) * 0.5  # ⇒ H = -G_F/2 < 0
        with self.assertRaises(ValueError):
            CohesiveDamageIsotropic(sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_low, softening='exponential')

    def test_linear_derived_quantities(self):
        m = make_linear()
        self.assertAlmostEqual(m.kappa_0, SIGMA_T0 / K_E, places=20)
        self.assertAlmostEqual(m.w_c, 2.0 * G_F / SIGMA_T0, places=20)
        self.assertIsNone(m.H)

    def test_exponential_derived_quantities(self):
        m = make_exponential()
        self.assertAlmostEqual(m.kappa_0, SIGMA_T0 / K_E, places=20)
        self.assertIsNone(m.w_c)
        H_expected = G_F - 0.5 * SIGMA_T0 * (SIGMA_T0 / K_E)
        self.assertAlmostEqual(m.H, H_expected, places=10)


# =========================================================================
# Verification — comportamiento físico monotónico
# =========================================================================


class TestElasticBelowThreshold(unittest.TestCase):
    """acceptance.respuesta_elastica_bajo_umbral"""

    def test_linear(self):
        m = make_linear()
        u_n = 0.5 * m.kappa_0
        t, T, state = m.compute_traction(np.array([u_n, 0.0]))
        self.assertAlmostEqual(state['damage'], 0.0, places=15)
        self.assertAlmostEqual(t[0], K_E * u_n, delta=1e-12 * K_E * u_n)
        self.assertEqual(t[1], 0.0)
        np.testing.assert_allclose(T, np.array([[K_E, 0.0], [0.0, 0.0]]), rtol=1e-14)


class TestDamageOnset(unittest.TestCase):
    """acceptance.inicio_dano_en_umbral"""

    def test_traction_at_threshold_is_sigma_t0(self):
        m = make_linear()
        eps = m.kappa_0 * 1.0e-12
        # En κ_0 exacto: ω = 0, t_n = K_e·κ_0 = σ_t0
        t0, _, st0 = m.compute_traction(np.array([m.kappa_0, 0.0]))
        self.assertAlmostEqual(st0['damage'], 0.0, places=15)
        self.assertAlmostEqual(t0[0], SIGMA_T0, delta=1e-9 * SIGMA_T0)
        # Inmediatamente por encima: ω positivo continuo y t ≈ σ_t0
        t1, _, st1 = m.compute_traction(np.array([m.kappa_0 + eps, 0.0]))
        self.assertGreater(st1['damage'], 0.0)
        self.assertLess(st1['damage'], 1e-6)  # continuo en κ_0
        self.assertAlmostEqual(t1[0], SIGMA_T0, delta=1e-6 * SIGMA_T0)


class TestFractureEnergyLinear(unittest.TestCase):
    """acceptance.energia_fractura_softening_lineal"""

    def test_integral_equals_G_f(self):
        m = make_linear()
        # Recorrer monótonamente [[u_n]] de 0 a w_c; cuadratura trapezoidal
        n = 5000
        u_grid = np.linspace(0.0, m.w_c, n)
        t_grid = np.zeros(n)
        state = None
        for i, u in enumerate(u_grid):
            t, _, state = m.compute_traction(np.array([u, 0.0]), state)
            t_grid[i] = t[0]
        integral = np.trapezoid(t_grid, u_grid)
        self.assertAlmostEqual(integral / G_F, 1.0, delta=1e-3)


class TestFractureEnergyExponential(unittest.TestCase):
    """acceptance.energia_fractura_softening_exponencial"""

    def test_integral_equals_G_f(self):
        m = make_exponential()
        u_max = 20.0 * G_F / SIGMA_T0  # cola suficientemente larga
        n = 20000
        u_grid = np.linspace(0.0, u_max, n)
        t_grid = np.zeros(n)
        state = None
        for i, u in enumerate(u_grid):
            t, _, state = m.compute_traction(np.array([u, 0.0]), state)
            t_grid[i] = t[0]
        integral = np.trapezoid(t_grid, u_grid)
        self.assertAlmostEqual(integral / G_F, 1.0, delta=1e-3)


class TestSecantUnloading(unittest.TestCase):
    """acceptance.descarga_secante"""

    def test_unload_to_zero_gives_zero_traction_and_frozen_damage(self):
        m = make_linear()
        # Cargar hasta κ ∈ (κ_0, w_c/2)
        u_load = 0.5 * m.w_c
        _, _, state = m.compute_traction(np.array([u_load, 0.0]))
        omega_loaded = state['damage']
        self.assertGreater(omega_loaded, 0.0)

        # Descarga a 0
        t_zero, T_zero, state2 = m.compute_traction(np.array([0.0, 0.0]), state)
        self.assertAlmostEqual(t_zero[0], 0.0, delta=1e-12 * SIGMA_T0)
        self.assertAlmostEqual(state2['damage'], omega_loaded, places=14)

        # Pendiente de descarga: tangente secante (1-ω)·K_e
        expected = (1.0 - omega_loaded) * K_E
        self.assertAlmostEqual(T_zero[0, 0], expected, delta=1e-9 * K_E)


class TestReloadSameKappaNoExtraDamage(unittest.TestCase):
    """acceptance.recarga_hasta_mismo_kappa_sin_dano_extra"""

    def test_reload_does_not_grow_kappa(self):
        m = make_linear()
        # Carga, descarga, recarga al mismo nivel
        u_peak = 0.4 * m.w_c
        _, _, st1 = m.compute_traction(np.array([u_peak, 0.0]))
        kappa_peak = st1['kappa']
        omega_peak = st1['damage']

        _, _, st2 = m.compute_traction(np.array([0.0, 0.0]), st1)
        _, _, st3 = m.compute_traction(np.array([u_peak, 0.0]), st2)
        self.assertAlmostEqual(st3['kappa'], kappa_peak, places=14)
        self.assertAlmostEqual(st3['damage'], omega_peak, places=14)


class TestReloadBeyondContinuesDamaging(unittest.TestCase):
    """acceptance.recarga_mas_alla_continua_danando"""

    def test_reload_above_old_kappa_grows_damage(self):
        m = make_linear()
        u_peak = 0.3 * m.w_c
        _, _, st1 = m.compute_traction(np.array([u_peak, 0.0]))
        omega1 = st1['damage']

        # Descarga parcial y recarga más alta
        _, _, st2 = m.compute_traction(np.array([0.1 * m.w_c, 0.0]), st1)
        self.assertAlmostEqual(st2['damage'], omega1, places=14)  # no crece

        u_new = 0.6 * m.w_c
        _, _, st3 = m.compute_traction(np.array([u_new, 0.0]), st2)
        self.assertGreater(st3['damage'], omega1)
        self.assertAlmostEqual(st3['kappa'], u_new, places=14)


# =========================================================================
# Specific — propiedades algorítmicas
# =========================================================================


class TestTangentSymmetry(unittest.TestCase):
    """acceptance.tangente_simetrica_en_carga_y_descarga"""

    def _check_symmetry(self, m, state_vars, jump):
        _, T, _ = m.compute_traction(jump, state_vars)
        frob = np.linalg.norm(T - T.T, 'fro')
        ref = max(np.linalg.norm(T, 'fro'), 1.0)
        self.assertLess(frob / ref, 1e-14)

    def test_linear_modes(self):
        m = make_linear()
        # 1. Sin daño activo
        self._check_symmetry(m, None, np.array([0.5 * m.kappa_0, 0.0]))
        # 2. Carga activa
        _, _, st = m.compute_traction(np.array([0.4 * m.w_c, 0.0]))
        self._check_symmetry(m, None, np.array([0.4 * m.w_c, 0.0]))
        # 3. Descarga
        self._check_symmetry(m, st, np.array([0.1 * m.w_c, 0.0]))
        # 4. Saturado
        self._check_symmetry(m, None, np.array([m.w_c * 1.01, 0.0]))


class TestTangentFiniteDifference(unittest.TestCase):
    """acceptance.tangente_diferencia_finita_consistente"""

    def test_loading_branch_matches_finite_difference(self):
        m = make_linear()
        # Punto en carga activa, fresh (state_vars=None)
        u0 = np.array([0.4 * m.w_c, 0.0])
        t0, T_an, _ = m.compute_traction(u0)

        h = 1e-7 * m.w_c   # perturbación pequeña pero por encima de eps_machine·escala
        T_fd = np.zeros((2, 2))
        for j in range(2):
            up = u0.copy()
            up[j] += h
            tp, _, _ = m.compute_traction(up)
            T_fd[:, j] = (tp - t0) / h

        # Filtrar comparación a la columna n: la dirección tangencial es
        # exactamente 0 por construcción en Modo-I y no aporta señal.
        err = np.linalg.norm(T_fd[:, 0] - T_an[:, 0]) / max(
            np.linalg.norm(T_an[:, 0]), 1.0
        )
        self.assertLess(err, 1e-5)


class TestTangentialDoesNotDamage(unittest.TestCase):
    """acceptance.tangencial_no_dana_ni_transmite"""

    def test_pure_tangential_jump_yields_zero_traction_and_no_damage(self):
        m = make_linear()
        u_s = 10.0 * m.kappa_0     # apertura tangencial arbitrariamente grande
        t, T, state = m.compute_traction(np.array([0.0, u_s]))
        self.assertAlmostEqual(t[0], 0.0, delta=1e-14)
        self.assertAlmostEqual(t[1], 0.0, delta=1e-14)
        self.assertAlmostEqual(state['damage'], 0.0, places=15)
        # T[1,1] = 0 por construcción (Modo-I)
        self.assertEqual(T[1, 1], 0.0)


class TestCompressionDoesNotDamage(unittest.TestCase):
    """acceptance.compresion_no_dana"""

    def test_pure_compression_keeps_damage_unchanged(self):
        m = make_linear()
        # Primero dañar parcialmente
        _, _, st = m.compute_traction(np.array([0.5 * m.w_c, 0.0]))
        omega_before = st['damage']

        # Aplicar penetración pura
        t, _, st2 = m.compute_traction(np.array([-m.kappa_0, 0.0]), st)
        self.assertAlmostEqual(st2['damage'], omega_before, places=14)
        # Caveat de cierre: t_n = (1-omega)·K_e·u_n con u_n < 0
        expected = (1.0 - omega_before) * K_E * (-m.kappa_0)
        self.assertAlmostEqual(t[0], expected, delta=1e-9 * abs(expected))


class TestSaturationLinear(unittest.TestCase):
    """acceptance.saturacion_en_w_c_softening_lineal

    El cap por ``DAMAGE_MAX`` se aplica sólo a la rigidez tangente (para
    mantener el Newton no-singular). El ``ω`` reportado y ``t_n`` reflejan
    el valor físico: ``ω = 1`` exacto y ``t_n = 0`` cuando ``κ ≥ w_c``.
    """

    def test_above_w_c_saturates(self):
        m = make_linear()
        u_n = 1.01 * m.w_c
        t, T, state = m.compute_traction(np.array([u_n, 0.0]))
        self.assertEqual(state['damage'], 1.0)
        self.assertEqual(t[0], 0.0)
        self.assertAlmostEqual(T[0, 0], (1.0 - DAMAGE_MAX) * K_E,
                               delta=1e-10 * K_E)


class TestDegenerateToElastic(unittest.TestCase):
    """acceptance.degeneracion_a_elasticidad_intacta"""

    def test_huge_sigma_t0_never_activates(self):
        m = CohesiveDamageIsotropic(
            sigma_t0=1.0e20, G_f=G_F, K_e=K_E, softening='linear'
        )
        u_n = 1.0e-3  # cualquier valor "humano"
        t, T, state = m.compute_traction(np.array([u_n, 0.0]))
        self.assertAlmostEqual(state['damage'], 0.0, places=15)
        self.assertAlmostEqual(t[0], K_E * u_n, delta=1e-10 * K_E * u_n)
        self.assertAlmostEqual(T[0, 0], K_E, delta=1e-10 * K_E)


# =========================================================================
# Arch — registro y atributos de clase
# =========================================================================


class TestRegistryAndClassAttributes(unittest.TestCase):
    """acceptance.registry_kind_es_cohesive + is_symmetric_attribute"""

    def test_registered_in_cohesive_registry(self):
        self.assertIn('CohesiveDamageIsotropic', CohesiveMaterialRegistry._items)
        cls = CohesiveMaterialRegistry._items['CohesiveDamageIsotropic']
        self.assertIs(cls, CohesiveDamageIsotropic)

    def test_is_symmetric_true(self):
        self.assertTrue(CohesiveDamageIsotropic.IS_SYMMETRIC)

    def test_jump_dim_2(self):
        self.assertEqual(CohesiveDamageIsotropic.JUMP_DIM, 2)

    def test_primary_state_var_damage(self):
        self.assertEqual(CohesiveDamageIsotropic.PRIMARY_STATE_VAR, 'damage')


if __name__ == '__main__':
    unittest.main()
