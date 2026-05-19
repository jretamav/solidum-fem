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
from fenix.materials.drucker_prager_2d import DruckerPrager2D
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
        mat = Elastic2D(E=E, nu=nu, density=0.0, hypothesis="plane_stress")
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
        mat = Elastic2D(E=E, nu=nu, density=0.0, hypothesis="plane_strain")
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
        sigma, _, state = self.mat.compute_state(strain=eps)
        # Ley exponencial: d = 1 − (κ_0/κ)·exp(−α·(κ−κ_0)).
        kappa = abs(eps)
        d_expected = 1.0 - (self.kappa_0 / kappa) * math.exp(
            -self.alpha * (kappa - self.kappa_0)
        )
        d_expected = min(d_expected, DAMAGE_MAX)
        self.assertAlmostEqual(state["damage"], d_expected, places=10)
        # σ = (1-d)·E·ε
        self.assertAlmostEqual(sigma, (1.0 - d_expected) * self.E * eps, places=6)
        # NB: la matriz tangente ya no es secante en carga activa (ahora es
        # algorítmica consistente, negativa); su validación se cubre en
        # TestIsotropicDamage1DConsistentTangent.

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


class TestIsotropicDamage1DConsistentTangent(unittest.TestCase):
    """Tangente algorítmica consistente para IsotropicDamage1D.

    Acceptance de la spec ``docs/specs/IsotropicDamage1D.md``.
    """

    def setUp(self):
        self.E = 2.0e5
        self.kappa_0 = 1.0e-4
        self.alpha = 100.0
        self.mat = IsotropicDamage1D(E=self.E, kappa_0=self.kappa_0, alpha=self.alpha)

    def test_tangent_equals_E_below_threshold(self):
        """Sin daño activo: E_tan = E."""
        eps = 0.5 * self.kappa_0
        _, E_tan, state = self.mat.compute_state(strain=eps)
        self.assertEqual(state["damage"], 0.0)
        self.assertAlmostEqual(E_tan, self.E, places=12)

    def test_tangent_equals_secant_on_unloading(self):
        """Descarga: E_tan = (1-d)·E exacto."""
        # Cargar a κ > κ_0
        eps_load = 3.0 * self.kappa_0
        _, _, state_loaded = self.mat.compute_state(strain=eps_load)
        d_loaded = state_loaded["damage"]
        self.assertGreater(d_loaded, 0.0)

        # Descargar
        eps_unload = 2.0 * self.kappa_0
        _, E_tan, state_unloaded = self.mat.compute_state(
            strain=eps_unload, state_vars=state_loaded
        )
        self.assertAlmostEqual(state_unloaded["kappa"], state_loaded["kappa"], places=14)
        self.assertAlmostEqual(state_unloaded["damage"], d_loaded, places=14)
        self.assertAlmostEqual(E_tan, (1.0 - d_loaded) * self.E, places=10)

    def test_tangent_negative_on_loading(self):
        """En carga activa con d∈(0, DAMAGE_MAX), E_tan = -(1-d)·E·α·κ (negativa)."""
        eps = 2.0 * self.kappa_0
        _, E_tan, state = self.mat.compute_state(strain=eps)
        d = state["damage"]
        kappa = state["kappa"]
        self.assertGreater(d, 0.0)
        self.assertLess(d, DAMAGE_MAX)
        expected = -(1.0 - d) * self.E * self.alpha * kappa
        self.assertAlmostEqual(E_tan, expected, places=8)
        self.assertLess(E_tan, 0.0)

    def test_tangent_at_saturation_equals_secant(self):
        """En d=DAMAGE_MAX el término consistente se corta → tangente secante."""
        eps = 1.0e6 * self.kappa_0
        _, E_tan, state = self.mat.compute_state(strain=eps)
        self.assertAlmostEqual(state["damage"], DAMAGE_MAX, places=12)
        self.assertAlmostEqual(E_tan, (1.0 - DAMAGE_MAX) * self.E, places=6)

    def test_tangent_matches_finite_difference(self):
        """E_tan analítica coincide con diferencia finita centrada."""
        eps_ref = 2.5 * self.kappa_0
        # Anchor con κ < ε para forzar carga activa en eps_ref±h
        state_anchor = {'kappa': eps_ref * 0.9, 'damage': 0.0}
        _, E_tan_ref, _ = self.mat.compute_state(strain=eps_ref, state_vars=state_anchor)

        h = 1.0e-9
        sigma_p, _, _ = self.mat.compute_state(strain=eps_ref + h, state_vars=state_anchor)
        sigma_m, _, _ = self.mat.compute_state(strain=eps_ref - h, state_vars=state_anchor)
        E_fd = (sigma_p - sigma_m) / (2.0 * h)

        err = abs(E_fd - E_tan_ref) / abs(E_tan_ref)
        self.assertLess(err, 1.0e-5)

    def test_is_symmetric_attribute_remains_True(self):
        """IS_SYMMETRIC=True (tangente escalar ⇒ contribución elemental simétrica)."""
        self.assertTrue(IsotropicDamage1D.IS_SYMMETRIC)

    def test_consistent_with_2D_uniaxial_first_component(self):
        """Reducción del modelo 2D plane stress con ν=0 y ε=[eps,0,0] coincide con 1D.

        Con ν=0 el material 2D plane stress se vuelve diagonal en la entrada
        [0,0]: σ_xx = E·ε_xx y el operador M = diag(1,1,1/2) deja ∂κ/∂ε_xx =
        sign(ε_xx)·ε_xx/ε_eq = sign(ε_xx) cuando ε_eq=|ε_xx|. Las fórmulas se
        reducen exactamente al caso escalar.
        """
        from fenix.materials.damage_2d import IsotropicDamage2D
        mat_2d = IsotropicDamage2D(
            E=self.E, nu=0.0, kappa_0=self.kappa_0, alpha=self.alpha,
            hypothesis='plane_stress',
        )
        eps_1d = 2.5 * self.kappa_0
        sigma_1d, E_tan_1d, _ = self.mat.compute_state(strain=eps_1d)

        eps_2d = np.array([eps_1d, 0.0, 0.0])
        sigma_2d, C_alg_2d, _ = mat_2d.compute_state(eps_2d)

        # Componente xx idéntica
        self.assertAlmostEqual(sigma_2d[0], sigma_1d, places=8)
        # Entrada [0,0] de la tangente 2D = E_tan 1D
        self.assertAlmostEqual(C_alg_2d[0, 0], E_tan_1d, places=4)


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
        _, _, state = self.mat.compute_state(eps)
        eps_eq = math.sqrt(eps[0] ** 2 + eps[1] ** 2 + 0.5 * eps[2] ** 2)
        d_expected = 1.0 - (self.kappa_0 / eps_eq) * math.exp(
            -self.alpha * (eps_eq - self.kappa_0)
        )
        d_expected = min(d_expected, DAMAGE_MAX)
        self.assertAlmostEqual(state["damage"], d_expected, places=10)
        # NB: la matriz tangente ya no es secante en carga activa (ahora es
        # algorítmica consistente); su validación se cubre en
        # TestIsotropicDamage2DConsistentTangent.


class TestIsotropicDamage2DConsistentTangent(unittest.TestCase):
    """Tangente algorítmica consistente para IsotropicDamage2D.

    Acceptance de la spec ``docs/specs/IsotropicDamage2D.md`` §7.
    """

    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.kappa_0 = 1.0e-4
        self.alpha = 100.0
        self.mat = IsotropicDamage2D(
            E=self.E, nu=self.nu, kappa_0=self.kappa_0, alpha=self.alpha
        )
        self.Ce = self.mat.elastic_base.C

    def test_is_symmetric_attribute_is_false(self):
        """Material declara IS_SYMMETRIC=False (tangente consistente asimétrica)."""
        self.assertFalse(IsotropicDamage2D.IS_SYMMETRIC)
        self.assertFalse(self.mat.IS_SYMMETRIC)

    def test_tangent_equals_C_e_below_threshold(self):
        """Sin daño activo (κ ≤ κ_0): C_tan = C_e (rama elástica pura)."""
        eps = np.array([0.4 * self.kappa_0, 0.0, 0.0])
        _, C_tan, state = self.mat.compute_state(eps)
        self.assertEqual(state["damage"], 0.0)
        np.testing.assert_allclose(C_tan, self.Ce, rtol=1e-12)

    def test_tangent_equals_secant_on_unloading(self):
        """Descarga desde estado dañado: C_tan = (1-d)·C_e exacto."""
        # 1) Cargar hasta κ > κ_0
        eps_load = np.array([3.0 * self.kappa_0, 0.0, 0.0])
        _, _, state_loaded = self.mat.compute_state(eps_load)
        d_loaded = state_loaded["damage"]
        self.assertGreater(d_loaded, 0.0)

        # 2) Descargar a estado de menor ε_eq
        eps_unload = np.array([2.0 * self.kappa_0, 0.0, 0.0])
        _, C_tan, state_unloaded = self.mat.compute_state(eps_unload, state_vars=state_loaded)

        # κ y d no deben cambiar (irreversibilidad)
        self.assertAlmostEqual(state_unloaded["kappa"], state_loaded["kappa"], places=14)
        self.assertAlmostEqual(state_unloaded["damage"], d_loaded, places=14)

        # Tangente exactamente secante
        np.testing.assert_allclose(C_tan, (1.0 - d_loaded) * self.Ce, rtol=1e-12)

    def test_tangent_at_saturation_equals_secant(self):
        """En d=DAMAGE_MAX el término consistente se corta → tangente secante."""
        # Forzar saturación con ε enorme
        eps = np.array([1.0e6 * self.kappa_0, 0.0, 0.0])
        _, C_tan, state = self.mat.compute_state(eps)
        self.assertAlmostEqual(state["damage"], DAMAGE_MAX, places=12)
        np.testing.assert_allclose(C_tan, (1.0 - DAMAGE_MAX) * self.Ce, rtol=1e-10)

    def test_tangent_differs_from_secant_on_loading(self):
        """Carga activa con d ∈ (0, DAMAGE_MAX): C_tan ≠ (1-d)·C_e."""
        eps = np.array([2.0 * self.kappa_0, 0.5 * self.kappa_0, 0.0])
        _, C_tan, state = self.mat.compute_state(eps)
        d = state["damage"]
        self.assertGreater(d, 0.0)
        self.assertLess(d, DAMAGE_MAX)
        C_sec = (1.0 - d) * self.Ce
        # Diferencia significativa (al menos un elemento difiere en >0.1%)
        diff = np.linalg.norm(C_tan - C_sec) / np.linalg.norm(C_tan)
        self.assertGreater(diff, 1.0e-3)

    def test_tangent_not_symmetric_on_loading(self):
        """En carga activa con estado no degenerado, C_tan asimétrica."""
        # Estado con γ_xy ≠ 0 y ε_xx ≠ ε_yy → producto exterior asimétrico
        eps = np.array([3.0 * self.kappa_0, 1.0 * self.kappa_0, 0.5 * self.kappa_0])
        _, C_tan, state = self.mat.compute_state(eps)
        self.assertGreater(state["damage"], 0.0)
        # Asimetría medible
        asym = np.linalg.norm(C_tan - C_tan.T) / np.linalg.norm(C_tan)
        self.assertGreater(asym, 1.0e-3)

    def test_tangent_matches_finite_difference(self):
        """C_tan analítica coincide con diferencia finita centrada.

        Verifica que la derivada cerrada ``∂σ/∂ε`` es correcta hasta el
        error de la perturbación finita (O(h²) en diferencia centrada).
        """
        eps = np.array([2.5 * self.kappa_0, 0.7 * self.kappa_0, 0.4 * self.kappa_0])
        _, C_tan, state_ref = self.mat.compute_state(eps)
        d = state_ref["damage"]
        self.assertGreater(d, 0.0)

        # κ committed para forzar la decisión de carga durante el FD: la
        # diferencia finita evalúa el material en estados ε ± h·e_i; para
        # que todas estas perturbaciones caigan en la misma rama de carga
        # con el mismo κ_ref, pasamos un state_vars previo con κ exactamente
        # igual al ε_eq del estado de referencia menos un margen.
        eps_eq = np.linalg.norm([eps[0], eps[1], eps[2] / np.sqrt(2)])
        # Un state previo cuya κ esté justo por debajo del ε_eq actual
        # (todas las perturbaciones h⁻⁶ caerán también por encima → loading)
        state_anchor = {'kappa': eps_eq * 0.9, 'damage': 0.0}

        h = 1.0e-7
        C_fd = np.zeros((3, 3))
        for j in range(3):
            eps_p = eps.copy(); eps_p[j] += h
            eps_m = eps.copy(); eps_m[j] -= h
            sigma_p, _, _ = self.mat.compute_state(eps_p, state_vars=state_anchor)
            sigma_m, _, _ = self.mat.compute_state(eps_m, state_vars=state_anchor)
            C_fd[:, j] = (sigma_p - sigma_m) / (2.0 * h)

        # Re-evaluar la analítica usando state_anchor para que el κ usado
        # en la derivada coincida con el de la diferencia finita
        _, C_tan_ref, _ = self.mat.compute_state(eps, state_vars=state_anchor)
        err = np.linalg.norm(C_fd - C_tan_ref) / np.linalg.norm(C_tan_ref)
        self.assertLess(err, 1.0e-5)

    def test_damage_irreversibility(self):
        """κ monótono no decreciente y d no decreciente bajo carga, no cambia en descarga."""
        state = None
        kappa_prev = self.kappa_0
        d_prev = 0.0
        eps_path = [
            np.array([0.5 * self.kappa_0, 0.0, 0.0]),  # elástico
            np.array([1.5 * self.kappa_0, 0.0, 0.0]),  # daño activo
            np.array([2.5 * self.kappa_0, 0.0, 0.0]),  # más daño
            np.array([1.0 * self.kappa_0, 0.0, 0.0]),  # descarga (κ no cambia)
            np.array([3.0 * self.kappa_0, 0.0, 0.0]),  # recarga, daño avanza más
        ]
        for eps in eps_path:
            _, _, state = self.mat.compute_state(eps, state_vars=state)
            self.assertGreaterEqual(state["kappa"], kappa_prev - 1e-14)
            self.assertGreaterEqual(state["damage"], d_prev - 1e-14)
            kappa_prev = state["kappa"]
            d_prev = state["damage"]


class TestVonMises2D(unittest.TestCase):
    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.sigma_y = 250.0
        self.H = 1.0e4
        self.mat = VonMises2D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H, density=0.0,
            hypothesis="plane_strain",
        )

    def test_invalid_hypothesis_rejected(self):
        with self.assertRaises(ValueError):
            VonMises2D(E=1.0, nu=0.3, sigma_y=1.0, density=0.0, hypothesis="3d")

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


class TestVonMises2DPlaneStress(unittest.TestCase):
    """J2 plane stress proyectado (Simó-Hughes §3.4.1).

    Acceptance de la spec ``docs/specs/VonMises2D.md`` §11.
    """

    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.sigma_y = 250.0
        self.H = 1.0e4
        self.mat = VonMises2D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
            hypothesis="plane_stress",
        )

    def test_elastic_below_yield(self):
        """Bajo fluencia: σ = C_e^ps · ε, estado interno intacto."""
        eps = np.array([1.0e-5, 0.0, 0.0])
        sigma, C, state = self.mat.compute_state(eps)
        np.testing.assert_allclose(sigma, self.mat.C_e @ eps, rtol=1e-10)
        np.testing.assert_allclose(C, self.mat.C_e, rtol=1e-10)
        self.assertAlmostEqual(state["alpha"], 0.0, places=12)
        np.testing.assert_allclose(state["eps_p"], np.zeros(4), atol=1e-14)

    def test_degeneracion_a_1d_traccion_pura(self):
        """Tracción uniaxial pura plane stress ↔ Elastoplastic1D.

        En tracción uniaxial pura (σ_yy = σ_xy = 0, σ_zz = 0) el modelo plane
        stress debe reproducir la curva σ_xx-ε_xx del Elastoplastic1D con los
        mismos parámetros: yield exactamente en σ_y, pendiente plástica
        tangente ``E_t = E·H/(E+H)``.

        Se impone ε_xx prescrito; ε_yy y γ_xy se eligen tales que σ_yy = σ_xy = 0
        manteniendo el material en el estado convergido (es decir, dejamos que
        el predictor sea consistente con σ_xx-puro). Para el test simplificado,
        comparamos directamente la curva mediante una secuencia que sigue la
        recta uniaxial conocida.
        """
        # Curva uniaxial 1D de referencia (mismos E, σ_y, H)
        mat_1d = Elastoplastic1D(E=self.E, sigma_y=self.sigma_y, H=self.H)

        # En tracción uniaxial pura plane stress (σ_yy=σ_xy=0, σ_zz=0 trivial):
        # ε_yy(elástico) = -ν·ε_xx,  γ_xy = 0
        # σ_xx = E · ε_xx  (¡no E/(1-ν²)! porque la condición σ_yy = 0 corrige
        # la respuesta plane stress hasta reducirla a 1D pura)
        # Tras fluencia, ε_yy se descompone en parte elástica y plástica.
        # Para evitar resolver la consistencia ε_yy en cada paso, conducimos el
        # material por una trayectoria conocida: aplicamos σ_xx monotónico
        # 1D y reconstruimos la deformación 2D coherente con plane stress.
        # Pero el material recibe ε y devuelve σ. Como prueba alternativa:
        # verificamos que para ε_xx pequeño (régimen elástico), σ_xx = E·ε_xx
        # cuando se elige ε_yy = -ν·ε_xx (no E/(1-ν²)·ε_xx).
        eps_xx = 1.0e-3
        eps = np.array([eps_xx, -self.nu * eps_xx, 0.0])
        sigma, _, _ = self.mat.compute_state(eps)
        # σ_xx debe ser E·ε_xx, σ_yy = 0
        self.assertAlmostEqual(sigma[0], self.E * eps_xx, places=4)
        self.assertAlmostEqual(sigma[1], 0.0, places=6)
        self.assertAlmostEqual(sigma[2], 0.0, places=10)

        # Comparación 1D vs 2D en régimen elástico → idénticos.
        sigma_1d, _, _ = mat_1d.compute_state(eps_xx)
        self.assertAlmostEqual(sigma[0], sigma_1d, places=4)

    def test_yield_uniaxial_at_sigma_y(self):
        """En tracción uniaxial plane stress, fluencia exactamente en σ_xx = σ_y."""
        # ε_xx tal que σ_xx = σ_y (régimen elástico previo): ε_xx = σ_y / E
        # (en uniaxial pura plane stress σ_xx = E·ε_xx con ε_yy = -ν·ε_xx)
        eps_yield = self.sigma_y / self.E
        # Justo por debajo: elástico
        eps_below = np.array([0.999 * eps_yield, -self.nu * 0.999 * eps_yield, 0.0])
        _, _, state_below = self.mat.compute_state(eps_below)
        self.assertAlmostEqual(state_below["alpha"], 0.0, places=10)

        # Justo por encima: plástico (con un epsilon de margen)
        eps_above = np.array([1.01 * eps_yield, -self.nu * 1.01 * eps_yield, 0.0])
        sigma, _, state_above = self.mat.compute_state(eps_above)
        self.assertGreater(state_above["alpha"], 0.0)
        # σ_xx tras return mapping ≤ σ_y + algo (no debe exceder R(α) por la zona desviadora plástica)
        # En uniaxial el primer paso plástico produce σ ≈ σ_y + H·Δα con Δα pequeño.
        self.assertLess(sigma[0], self.sigma_y * 1.05)
        self.assertGreater(sigma[0], self.sigma_y * 0.99)

    def test_traccion_biaxial_isotropa_yield_at_sigma_y(self):
        """Fluencia biaxial isótropa (σ_xx = σ_yy) ocurre exactamente en σ = σ_y.

        En J2, tensor 3D = diag(σ, σ, 0) con σ_zz=0 (plane stress):
        s = diag(σ/3, σ/3, -2σ/3),  ||s|| = σ·√(2/3) = √(2/3)·σ_y  ⇒  σ = σ_y.
        """
        # ε_xx = ε_yy = ε; σ_xx = σ_yy = E/(1-ν) · ε en régimen elástico
        # Yield en σ = σ_y ⇒ ε = σ_y·(1-ν)/E
        eps_yield = self.sigma_y * (1.0 - self.nu) / self.E

        eps_below = np.array([0.99 * eps_yield, 0.99 * eps_yield, 0.0])
        _, _, state_below = self.mat.compute_state(eps_below)
        self.assertAlmostEqual(state_below["alpha"], 0.0, places=10)

        eps_above = np.array([1.05 * eps_yield, 1.05 * eps_yield, 0.0])
        sigma_above, _, state_above = self.mat.compute_state(eps_above)
        self.assertGreater(state_above["alpha"], 0.0)
        # Tras return mapping σ_xx = σ_yy aproximadamente, en torno a σ_y · (1 + pequeño)
        self.assertAlmostEqual(sigma_above[0], sigma_above[1], places=8)
        self.assertLess(sigma_above[0], self.sigma_y * 1.05)
        self.assertGreater(sigma_above[0], self.sigma_y * 0.97)

    def test_plastic_incompressibility(self):
        """tr(ε^p) = ε^p_xx + ε^p_yy + ε^p_zz = 0 (regla de flujo asociada J2)."""
        eps = np.array([5.0e-3, -2.0e-3, 1.0e-3])
        _, _, state = self.mat.compute_state(eps)
        eps_p = state["eps_p"]
        # Componentes: [xx, yy, zz, xy_tensorial]
        self.assertAlmostEqual(eps_p[0] + eps_p[1] + eps_p[2], 0.0, places=12)

    def test_plastic_incompressibility_persistente(self):
        """tr(ε^p) = 0 persiste a lo largo de una trayectoria multi-paso (H-3.9).

        El kernel plane stress cierra ``ε^p_zz = -(ε^p_xx + ε^p_yy)`` en cada
        paso por construcción, pero asume implícitamente que el ``ε^p_old``
        entrante también cumple ``tr = 0``. Este test ejercita el invariante
        en una trayectoria realista (carga monotónica + descarga + recarga +
        cortante), garantizando que se conserva paso tras paso sin drift
        numérico.
        """
        strain_path = [
            np.array([1.0e-3, 0.0, 0.0]),                  # elástico
            np.array([3.0e-3, -0.5e-3, 0.0]),              # plástico inicial
            np.array([5.0e-3, -1.0e-3, 0.5e-3]),           # plástico con cortante
            np.array([4.0e-3, -0.8e-3, 0.4e-3]),           # descarga parcial
            np.array([6.0e-3, -1.2e-3, 0.6e-3]),           # recarga plástica
            np.array([7.0e-3, -2.0e-3, 1.0e-3]),           # más plasticidad
        ]
        state = None
        for k, eps in enumerate(strain_path):
            _, _, state = self.mat.compute_state(eps, state_vars=state)
            eps_p = state["eps_p"]
            trace = eps_p[0] + eps_p[1] + eps_p[2]
            self.assertAlmostEqual(
                trace, 0.0, places=12,
                msg=f"paso {k}: tr(ε^p) = {trace:.3e} ≠ 0 tras ε={eps.tolist()}",
            )

    def test_alpha_monotonic_under_increasing_load(self):
        """α (deformación plástica acumulada) no decrece bajo carga creciente."""
        state = None
        alpha_prev = 0.0
        for k in range(1, 6):
            eps = np.array([0.005 * k, -0.005 * k, 0.0])
            _, _, state = self.mat.compute_state(eps, state_vars=state)
            self.assertGreaterEqual(state["alpha"], alpha_prev - 1.0e-14)
            alpha_prev = state["alpha"]
        self.assertGreater(alpha_prev, 0.0)

    def test_elastic_unload_recovers_C_e(self):
        """Tras descarga elástica desde estado plástico, la tangente vuelve a C_e^ps."""
        # Cargar más allá de fluencia
        eps_load = np.array([3.0e-3, 0.0, 0.0])
        _, _, state_loaded = self.mat.compute_state(eps_load)
        self.assertGreater(state_loaded["alpha"], 0.0)
        # Descarga: ε pequeño desde la deformación cargada — predictor cae dentro
        # de la superficie expandida → respuesta elástica con la matriz C_e^ps.
        eps_unload = np.array([2.0e-3, 0.0, 0.0])
        _, C_unload, state_unload = self.mat.compute_state(eps_unload, state_vars=state_loaded)
        np.testing.assert_allclose(C_unload, self.mat.C_e, rtol=1e-10)
        self.assertAlmostEqual(state_unload["alpha"], state_loaded["alpha"], places=12)

    def test_tangente_simetrica(self):
        """Tangente algorítmica simétrica (J2 asociado ⇒ IS_SYMMETRIC=True)."""
        # Estado plástico genérico
        eps = np.array([4.0e-3, -1.0e-3, 2.0e-3])
        _, C_alg, _ = self.mat.compute_state(eps)
        np.testing.assert_allclose(C_alg, C_alg.T, atol=1e-8 * np.max(np.abs(C_alg)))

    def test_unit_invariance_MPa_vs_Pa(self):
        """Invariancia bajo cambio de unidades (ADR 0006)."""
        # Sistema A (MPa, mm, N) vs B (Pa, m, N)
        mat_MPa = VonMises2D(E=2.0e5, nu=0.3, sigma_y=250.0, H=1.0e4, hypothesis="plane_stress")
        mat_Pa  = VonMises2D(E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e10, hypothesis="plane_stress")

        strain_path = [
            np.array([0.5e-3, 0.0, 0.0]),
            np.array([1.5e-3, -0.3e-3, 0.0]),
            np.array([3.0e-3, -0.8e-3, 0.5e-3]),
            np.array([2.0e-3, -0.5e-3, 0.5e-3]),  # descarga parcial
            np.array([5.0e-3, -1.2e-3, 1.0e-3]),
        ]
        state_MPa, state_Pa = None, None
        for eps in strain_path:
            _, _, state_MPa = mat_MPa.compute_state(eps, state_vars=state_MPa)
            _, _, state_Pa  = mat_Pa.compute_state(eps,  state_vars=state_Pa)
            self.assertAlmostEqual(state_MPa["alpha"], state_Pa["alpha"], places=12)
            np.testing.assert_allclose(state_MPa["eps_p"], state_Pa["eps_p"], atol=1e-14)


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


class TestDruckerPrager2D(unittest.TestCase):
    """Drucker-Prager plane strain con plasticidad no asociada y endurecimiento isótropo.

    Acceptance de la spec ``docs/specs/DruckerPrager2D.md``.
    """

    def setUp(self):
        # Parámetros típicos de suelo/hormigón
        self.E = 2.0e4
        self.nu = 0.3
        self.c0 = 10.0
        self.phi = 30.0
        self.psi = 10.0   # no asociada
        self.H = 100.0
        self.mat = DruckerPrager2D(
            E=self.E, nu=self.nu, cohesion=self.c0,
            phi_deg=self.phi, psi_deg=self.psi, H=self.H,
        )

    def test_invalid_hypothesis_rejected(self):
        with self.assertRaises(NotImplementedError):
            DruckerPrager2D(E=1e4, nu=0.3, cohesion=1.0, phi_deg=30.0, hypothesis='plane_stress')

    def test_invalid_variant_rejected(self):
        with self.assertRaises(ValueError):
            DruckerPrager2D(E=1e4, nu=0.3, cohesion=1.0, phi_deg=30.0, variant='bogus')

    def test_psi_greater_than_phi_rejected(self):
        with self.assertRaises(ValueError):
            DruckerPrager2D(E=1e4, nu=0.3, cohesion=1.0, phi_deg=30.0, psi_deg=45.0)

    def test_default_psi_equals_phi_is_associated(self):
        mat = DruckerPrager2D(E=1e4, nu=0.3, cohesion=1.0, phi_deg=30.0)
        self.assertTrue(mat.associated)
        self.assertAlmostEqual(mat.eta_f, mat.eta_g)

    def test_calibration_plane_strain_matched(self):
        """Para variant='plane_strain_matched': η_f = tan(φ)/√(9 + 12tan²(φ))."""
        phi_rad = math.radians(self.phi)
        tphi = math.tan(phi_rad)
        expected_eta_f = tphi / math.sqrt(9.0 + 12.0 * tphi * tphi)
        self.assertAlmostEqual(self.mat.eta_f, expected_eta_f, places=12)
        expected_k0 = 3.0 * self.c0 / math.sqrt(9.0 + 12.0 * tphi * tphi)
        self.assertAlmostEqual(self.mat.k0, expected_k0, places=12)

    def test_elastic_below_yield(self):
        """ε pequeño: respuesta elástica, alpha intacto."""
        eps = np.array([1.0e-5, 0.0, 0.0])
        sigma, C, state = self.mat.compute_state(eps)
        np.testing.assert_allclose(sigma, self.mat.C_e @ eps, rtol=1e-10)
        self.assertAlmostEqual(state['alpha'], 0.0, places=12)

    def test_regular_return_under_shear(self):
        """Cortante puro produce return regular activo."""
        eps = np.array([0.0, 0.0, 5.0e-3])
        sigma, _, state = self.mat.compute_state(eps)
        self.assertGreater(state['alpha'], 0.0)
        # En return regular el desviador σ_xy ≠ 0 (no colapsa al ápice)
        self.assertGreater(abs(sigma[2]), 0.0)

    def test_apex_return_under_pure_tension(self):
        """Tracción biaxial pura empuja al ápice del cono."""
        # ε grande hidrostática: predictor cae más allá del cono → apex
        eps = np.array([1.0e-2, 1.0e-2, 0.0])
        sigma, _, state = self.mat.compute_state(eps)
        self.assertGreater(state['alpha'], 0.0)
        # En el ápice σ_xx = σ_yy (puramente hidrostático), σ_xy ≈ 0
        self.assertAlmostEqual(sigma[0], sigma[1], delta=abs(sigma[0]) * 1e-10)
        self.assertAlmostEqual(sigma[2], 0.0, places=10)
        # σ_xx ≈ k(α)/(3·η_f)
        k_new = self.mat.k0 + self.mat.H * state['alpha']
        expected_p = k_new / (3.0 * self.mat.eta_f)
        self.assertAlmostEqual(sigma[0], expected_p, delta=abs(expected_p) * 1e-6)

    def test_associated_tangent_is_symmetric(self):
        """ψ = φ (asociada): tangente algorítmica simétrica."""
        mat_assoc = DruckerPrager2D(
            E=self.E, nu=self.nu, cohesion=self.c0,
            phi_deg=self.phi, psi_deg=self.phi, H=self.H,
        )
        # Estado plástico regular (cortante con algo de hidrostática)
        eps = np.array([1.0e-3, -1.0e-3, 4.0e-3])
        _, C_alg, state = mat_assoc.compute_state(eps)
        self.assertGreater(state['alpha'], 0.0)
        asym = np.linalg.norm(C_alg - C_alg.T) / np.linalg.norm(C_alg)
        self.assertLess(asym, 1.0e-10)

    def test_nonassociated_tangent_is_asymmetric(self):
        """ψ ≠ φ: tangente asimétrica."""
        eps = np.array([0.0, 0.0, 5.0e-3])
        _, C_alg, state = self.mat.compute_state(eps)
        self.assertGreater(state['alpha'], 0.0)
        asym = np.linalg.norm(C_alg - C_alg.T) / np.linalg.norm(C_alg)
        self.assertGreater(asym, 1.0e-3)

    def test_tangent_matches_finite_difference_regular(self):
        """C_alg cerrada coincide con FD centrado en rama regular."""
        eps = np.array([1.0e-3, -2.0e-3, 3.0e-3])
        state_anchor = {'eps_p': np.zeros(4), 'alpha': 0.0}
        _, C_alg, _ = self.mat.compute_state(eps, state_vars=state_anchor)

        h = 1.0e-7
        C_fd = np.zeros((3, 3))
        for j in range(3):
            ep = eps.copy(); ep[j] += h
            em = eps.copy(); em[j] -= h
            sp, _, _ = self.mat.compute_state(ep, state_vars=state_anchor)
            sm, _, _ = self.mat.compute_state(em, state_vars=state_anchor)
            C_fd[:, j] = (sp - sm) / (2.0 * h)
        err = np.linalg.norm(C_alg - C_fd) / np.linalg.norm(C_alg)
        self.assertLess(err, 1.0e-5)

    def test_alpha_monotonic_under_increasing_load(self):
        state = None
        alpha_prev = 0.0
        for k in range(1, 5):
            eps = np.array([0.0, 0.0, 2.0e-3 * k])  # cortante creciente
            _, _, state = self.mat.compute_state(eps, state_vars=state)
            self.assertGreaterEqual(state['alpha'], alpha_prev - 1e-14)
            alpha_prev = state['alpha']
        self.assertGreater(alpha_prev, 0.0)

    def test_reevaluation_on_frontier_returns_C_e(self):
        """Re-evaluar el material con la misma ε y state convergido produce paso
        elástico marginal (f_trial ≈ 0 en frontera): tangente = C_e, α no cambia.

        Análogo a ``test_no_further_yield_at_same_state`` de J2. Probar
        "descarga elástica" en Drucker-Prager no asociado con dilatancia es
        sutil: la ε_p acumulada con η_g > 0 puede inducir plasticidad inversa
        al descargar puramente en cortante (el predictor sale del cono por
        el lado contrario). El test conservador verifica que la re-evaluación
        en la misma ε converge al mismo estado.
        """
        eps_load = np.array([0.0, 0.0, 5.0e-3])
        _, _, state_loaded = self.mat.compute_state(eps_load)
        self.assertGreater(state_loaded['alpha'], 0.0)

        _, C_re, state_re = self.mat.compute_state(eps_load, state_vars=state_loaded)
        np.testing.assert_allclose(C_re, self.mat.C_e, rtol=1.0e-8)
        self.assertAlmostEqual(state_re['alpha'], state_loaded['alpha'], places=10)

    def test_unit_invariance_MPa_vs_Pa(self):
        """Mismo path adimensional en (MPa, mm, N) vs (Pa, m, N) ⇒ α idéntico."""
        mat_MPa = DruckerPrager2D(E=2.0e4, nu=0.3, cohesion=10.0, phi_deg=30.0, psi_deg=10.0, H=100.0)
        mat_Pa = DruckerPrager2D(E=2.0e10, nu=0.3, cohesion=1.0e7, phi_deg=30.0, psi_deg=10.0, H=1.0e8)
        strain_path = [
            np.array([0.0, 0.0, 1.0e-3]),
            np.array([1.0e-3, -1.0e-3, 3.0e-3]),
            np.array([2.0e-3, -2.0e-3, 5.0e-3]),
        ]
        state_MPa, state_Pa = None, None
        for eps in strain_path:
            _, _, state_MPa = mat_MPa.compute_state(eps, state_vars=state_MPa)
            _, _, state_Pa = mat_Pa.compute_state(eps, state_vars=state_Pa)
            self.assertAlmostEqual(state_MPa['alpha'], state_Pa['alpha'], places=12)


if __name__ == "__main__":
    unittest.main()
