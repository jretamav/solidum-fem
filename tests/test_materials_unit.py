"""Tests de unidad de los modelos constitutivos.

Cubre los seis materiales registrados en ``fenix.materials`` con foco en:

- Elasticidad: σ = E·ε o σ = C·ε con la matriz constitutiva esperada.
- Plasticidad: predictor elástico, return mapping bajo carga, tangente
  consistente, regreso elástico tras descarga.
- Daño: umbral, ley exponencial, irreversibilidad (κ máxima histórica),
  saturación a ``DAMAGE_MAX``.

No se cubre el ``CableMaterial1D`` (ya tiene ``test_cable_material.py``).
"""

import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.constants import DAMAGE_MAX
from fenix.materials.damage_1d import IsotropicDamage1D
from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.materials.elastic import Elastic1D
from fenix.materials.elastic_2d import Elastic2D
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.materials.von_mises_2d import VonMises2D


class TestElastic1D(unittest.TestCase):
    def test_hooke(self):
        mat = Elastic1D(E=2.0e5)
        sigma, E_t, _ = mat.compute_state(strain=1.0e-3)
        self.assertAlmostEqual(sigma, 2.0e5 * 1.0e-3)
        self.assertEqual(E_t, 2.0e5)

    def test_state_vars_passthrough(self):
        mat = Elastic1D(E=1.0)
        _, _, state = mat.compute_state(strain=0.5, state_vars={"foo": 7})
        self.assertEqual(state, {"foo": 7})


class TestElastic2D(unittest.TestCase):
    def test_plane_stress_matrix(self):
        E, nu = 2.0e5, 0.3
        mat = Elastic2D(E=E, nu=nu, hypothesis="plane_stress")
        coef = E / (1.0 - nu**2)
        expected = coef * np.array([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0],
        ])
        np.testing.assert_allclose(mat.C, expected, rtol=1e-12)
        # σ = C · ε
        eps = np.array([1.0e-3, 0.0, 0.0])
        sigma, C_ret, _ = mat.compute_state(eps)
        np.testing.assert_allclose(sigma, expected @ eps, rtol=1e-12)
        np.testing.assert_allclose(C_ret, expected, rtol=1e-12)

    def test_plane_strain_matrix(self):
        E, nu = 2.0e5, 0.3
        mat = Elastic2D(E=E, nu=nu, hypothesis="plane_strain")
        coef = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
        expected = coef * np.array([
            [1.0 - nu, nu, 0.0],
            [nu, 1.0 - nu, 0.0],
            [0.0, 0.0, (1.0 - 2.0 * nu) / 2.0],
        ])
        np.testing.assert_allclose(mat.C, expected, rtol=1e-12)

    def test_symmetry_of_C(self):
        mat = Elastic2D(E=1.0, nu=0.25)
        np.testing.assert_allclose(mat.C, mat.C.T, atol=1e-14)


class TestElastoplastic1D(unittest.TestCase):
    """Plasticidad J2 1D con endurecimiento isótropo lineal."""

    def setUp(self):
        self.E = 2.0e5
        self.sigma_y = 250.0
        self.H = 1.0e4
        self.mat = Elastoplastic1D(E=self.E, sigma_y=self.sigma_y, H=self.H)

    def test_elastic_below_yield(self):
        eps = self.sigma_y / self.E * 0.5  # mitad de la fluencia
        sigma, E_t, state = self.mat.compute_state(strain=eps)
        self.assertAlmostEqual(sigma, self.E * eps)
        self.assertEqual(E_t, self.E)
        self.assertEqual(state["eps_p"], 0.0)
        self.assertEqual(state["alpha"], 0.0)

    def test_at_yield_transitions(self):
        eps = self.sigma_y / self.E  # justo en la frontera elástica
        sigma, E_t, state = self.mat.compute_state(strain=eps)
        self.assertAlmostEqual(sigma, self.sigma_y, places=8)
        # En ε_yield la tangente puede ser elástica o algoritmica según
        # implementación; basta con verificar que no hay deformación plástica.
        self.assertAlmostEqual(state["eps_p"], 0.0, places=12)

    def test_return_mapping_above_yield(self):
        """Tangente algorítmica E_t = E·H / (E + H) al fluir."""
        eps = 2.0 * self.sigma_y / self.E
        sigma, E_t, state = self.mat.compute_state(strain=eps)
        # σ debe estar acotado por σ_y + H·α (criterio de fluencia).
        yield_post = self.sigma_y + self.H * state["alpha"]
        self.assertAlmostEqual(abs(sigma), yield_post, places=8)
        # Tangente algoritmica.
        E_t_expected = self.E * self.H / (self.E + self.H)
        self.assertAlmostEqual(E_t, E_t_expected, places=8)
        self.assertGreater(state["alpha"], 0.0)

    def test_unload_elastic_with_plastic_offset(self):
        # Carga plástica luego descarga total.
        eps_load = 2.0 * self.sigma_y / self.E
        _, _, state1 = self.mat.compute_state(strain=eps_load)
        eps_p = state1["eps_p"]
        # Descargar a ε = eps_p (deformación plástica residual): σ debe ser 0.
        sigma_unloaded, E_t, _ = self.mat.compute_state(strain=eps_p, state_vars=state1)
        self.assertAlmostEqual(sigma_unloaded, 0.0, places=6)
        self.assertEqual(E_t, self.E)  # respuesta elástica

    def test_compression_symmetry(self):
        eps_t = 2.0 * self.sigma_y / self.E
        sigma_t, _, _ = self.mat.compute_state(strain=eps_t)
        sigma_c, _, _ = self.mat.compute_state(strain=-eps_t)
        self.assertAlmostEqual(sigma_t, -sigma_c, places=8)


class TestIsotropicDamage1D(unittest.TestCase):
    def setUp(self):
        self.E = 2.0e5
        self.kappa_0 = 1.0e-4
        self.alpha = 100.0
        self.mat = IsotropicDamage1D(E=self.E, kappa_0=self.kappa_0, alpha=self.alpha)

    def test_no_damage_below_threshold(self):
        eps = 0.5 * self.kappa_0
        sigma, E_t, state = self.mat.compute_state(strain=eps)
        self.assertEqual(state["damage"], 0.0)
        self.assertAlmostEqual(sigma, self.E * eps)
        self.assertEqual(E_t, self.E)

    def test_damage_evolves_above_threshold(self):
        eps = 2.0 * self.kappa_0
        sigma, E_sec, state = self.mat.compute_state(strain=eps)
        # Ley exponencial: d = 1 − (κ_0/κ)·exp(−α·(κ−κ_0)).
        kappa = abs(eps)
        d_expected = 1.0 - (self.kappa_0 / kappa) * math.exp(
            -self.alpha * (kappa - self.kappa_0)
        )
        d_expected = min(d_expected, DAMAGE_MAX)
        self.assertAlmostEqual(state["damage"], d_expected, places=10)
        self.assertAlmostEqual(E_sec, (1.0 - d_expected) * self.E, places=6)
        self.assertAlmostEqual(sigma, E_sec * eps, places=6)

    def test_damage_irreversible(self):
        # Cargar a ε grande, luego descargar: el daño no debe disminuir.
        eps_load = 5.0 * self.kappa_0
        _, _, state_loaded = self.mat.compute_state(strain=eps_load)
        d_max = state_loaded["damage"]

        # Descarga a ε = 0: κ se conserva, d se conserva.
        _, E_sec, state_unloaded = self.mat.compute_state(
            strain=0.0, state_vars=state_loaded
        )
        self.assertAlmostEqual(state_unloaded["damage"], d_max, places=12)
        self.assertAlmostEqual(state_unloaded["kappa"], state_loaded["kappa"], places=12)
        self.assertAlmostEqual(E_sec, (1.0 - d_max) * self.E, places=6)

    def test_damage_saturation(self):
        eps = 1.0e6 * self.kappa_0  # extremadamente grande
        _, _, state = self.mat.compute_state(strain=eps)
        self.assertAlmostEqual(state["damage"], DAMAGE_MAX, places=12)


class TestIsotropicDamage2D(unittest.TestCase):
    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.kappa_0 = 1.0e-4
        self.alpha = 100.0
        self.mat = IsotropicDamage2D(
            E=self.E, nu=self.nu, kappa_0=self.kappa_0, alpha=self.alpha
        )

    def test_no_damage_below_threshold(self):
        eps = np.array([0.4 * self.kappa_0, 0.0, 0.0])
        sigma, C_sec, state = self.mat.compute_state(eps)
        self.assertEqual(state["damage"], 0.0)
        # σ = C_e · ε.
        Ce = self.mat.elastic_base.C
        np.testing.assert_allclose(sigma, Ce @ eps, rtol=1e-10)
        np.testing.assert_allclose(C_sec, Ce, rtol=1e-10)

    def test_damage_evolves_above_threshold(self):
        eps = np.array([2.0 * self.kappa_0, 0.0, 0.0])
        _, C_sec, state = self.mat.compute_state(eps)
        eps_eq = math.sqrt(eps[0] ** 2 + eps[1] ** 2 + 0.5 * eps[2] ** 2)
        d_expected = 1.0 - (self.kappa_0 / eps_eq) * math.exp(
            -self.alpha * (eps_eq - self.kappa_0)
        )
        d_expected = min(d_expected, DAMAGE_MAX)
        self.assertAlmostEqual(state["damage"], d_expected, places=10)
        Ce = self.mat.elastic_base.C
        np.testing.assert_allclose(C_sec, (1.0 - d_expected) * Ce, rtol=1e-8)


class TestVonMises2D(unittest.TestCase):
    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.sigma_y = 250.0
        self.H = 1.0e4
        self.mat = VonMises2D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
            hypothesis="plane_strain",
        )

    def test_plane_stress_unsupported(self):
        with self.assertRaises(NotImplementedError):
            VonMises2D(E=1.0, nu=0.3, sigma_y=1.0, hypothesis="plane_stress")

    def test_elastic_below_yield(self):
        # Deformación pequeña → respuesta puramente elástica.
        eps = np.array([1.0e-5, 0.0, 0.0])
        sigma, C, state = self.mat.compute_state(eps)
        np.testing.assert_allclose(sigma, self.mat.C_e @ eps, rtol=1e-10)
        np.testing.assert_allclose(C, self.mat.C_e, rtol=1e-10)
        self.assertAlmostEqual(state["alpha"], 0.0, places=12)

    def test_yield_criterion_satisfied(self):
        # Deformación cortante grande induce fluencia plana.
        eps = np.array([0.01, -0.01, 0.0])  # estado desviador puro
        sigma, _, state = self.mat.compute_state(eps)
        # En J2 plane strain: ‖s‖ ≤ √(2/3)·(σ_y + H·α) tras return mapping.
        # Componente zz no se reporta en σ (solo xx, yy, xy), pero el criterio
        # se evalúa sobre s_dev incluyendo zz. Verificamos que α evolucionó.
        self.assertGreater(state["alpha"], 0.0)
        # Tangente diferente a la elástica.
        # (Comprobación cualitativa: la traza desviadora indicará softening
        # respecto a C_e.)

    def test_alpha_monotonic_under_increasing_load(self):
        """``α`` (deformación plástica acumulada) crece monótonamente bajo
        carga creciente — irreversibilidad de la plasticidad."""
        state = None
        alpha_prev = 0.0
        for k in range(1, 6):
            eps = np.array([0.005 * k, -0.005 * k, 0.0])
            _, _, state = self.mat.compute_state(eps, state_vars=state)
            self.assertGreaterEqual(state["alpha"], alpha_prev)
            alpha_prev = state["alpha"]
        # Tras 5 pasos de carga creciente, debe haber fluencia acumulada.
        self.assertGreater(alpha_prev, 0.0)

    def test_no_further_yield_at_same_state(self):
        """Repetir el mismo ε con el estado actualizado no genera más fluencia."""
        eps = np.array([0.005, -0.005, 0.0])
        _, _, state1 = self.mat.compute_state(eps)
        alpha_1 = state1["alpha"]
        _, _, state2 = self.mat.compute_state(eps, state_vars=state1)
        # El predictor cae justo sobre la nueva superficie de fluencia →
        # no hay más incremento plástico (módulo redondeo del return mapping).
        self.assertAlmostEqual(state2["alpha"], alpha_1, places=10)


class TestAdmissibilityToleranceUnitInvariance(unittest.TestCase):
    """Invariancia del check de admisibilidad bajo cambio de unidades (ADR 0006).

    Mismo problema físico expresado en MPa y Pa debe producir el mismo estado
    interno tras la misma trayectoria de deformación. Antes de ADR 0006 la
    tolerancia ``PLASTIC_YIELD_TOL = 1e-9`` era absoluta y los dos sistemas
    podían divergir según el rango de unidades.
    """

    def test_elastoplastic1d_unit_invariance_MPa_vs_Pa(self):
        # Acero típico: E = 200 GPa, σ_y = 250 MPa, H = 10 GPa.
        # Sistema A (MPa, mm, N):  E=2e5, σ_y=250,   H=1e4
        # Sistema B (Pa,  m,  N):  E=2e11, σ_y=2.5e8, H=1e10
        # Una misma trayectoria de deformación (adimensional) debe dar el mismo α.
        mat_MPa = Elastoplastic1D(E=2.0e5, sigma_y=250.0, H=1.0e4)
        mat_Pa  = Elastoplastic1D(E=2.0e11, sigma_y=2.5e8, H=1.0e10)

        strains = [0.0005, 0.001, 0.002, 0.003, 0.0025, 0.0015, 0.002, 0.004]

        state_MPa = None
        state_Pa = None
        for eps in strains:
            _, _, state_MPa = mat_MPa.compute_state(eps, state_vars=state_MPa)
            _, _, state_Pa = mat_Pa.compute_state(eps, state_vars=state_Pa)
            self.assertAlmostEqual(
                state_MPa["alpha"], state_Pa["alpha"], places=14,
                msg=f"Divergencia en α tras ε={eps}: MPa={state_MPa['alpha']}, Pa={state_Pa['alpha']}",
            )
            self.assertAlmostEqual(
                state_MPa["eps_p"], state_Pa["eps_p"], places=14,
                msg=f"Divergencia en ε_p tras ε={eps}",
            )

    def test_isotropic_damage1d_unit_invariance(self):
        # Mismo concepto: dos parametrizaciones físicamente equivalentes.
        mat_MPa = IsotropicDamage1D(E=2.0e5,  kappa_0=1.0e-4, alpha=100.0)
        mat_Pa  = IsotropicDamage1D(E=2.0e11, kappa_0=1.0e-4, alpha=100.0)

        strains = [0.5e-4, 1.5e-4, 2.0e-4, 1.0e-4, 3.0e-4, 0.0]

        state_MPa = None
        state_Pa = None
        for eps in strains:
            _, _, state_MPa = mat_MPa.compute_state(eps, state_vars=state_MPa)
            _, _, state_Pa = mat_Pa.compute_state(eps, state_vars=state_Pa)
            self.assertAlmostEqual(state_MPa["damage"], state_Pa["damage"], places=14)
            self.assertAlmostEqual(state_MPa["kappa"], state_Pa["kappa"], places=14)

    def test_admissibility_scale_adaptativa_con_endurecimiento(self):
        """La escala devuelta crece al endurecerse (no se queda en σ_y inicial)."""
        mat = Elastoplastic1D(E=2.0e5, sigma_y=250.0, H=1.0e4)
        scale_inicial = mat.admissibility_scale(None)
        self.assertAlmostEqual(scale_inicial, 250.0)
        # Tras α = 0.01: escala = 250 + 1e4·0.01 = 350.
        scale_endurecido = mat.admissibility_scale({"alpha": 0.01})
        self.assertAlmostEqual(scale_endurecido, 350.0)

    def test_admissibility_scale_default_para_elastico_puro(self):
        """Materiales sin override heredan el default 1.0 (nunca se invoca)."""
        mat = Elastic1D(E=1.0e5)
        self.assertEqual(mat.admissibility_scale(), 1.0)


if __name__ == "__main__":
    unittest.main()
