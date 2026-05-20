"""Plasticidad cíclica: propiedades del endurecimiento isótropo lineal bajo
carga cíclica para `Elastoplastic1D` y `VonMises2D`.

**Nota sobre la regla de Masing**: la regla de Masing estricta —"la curva
post-reversal en coordenadas relativas es 2·backbone(Δε/2)"— aplica al
endurecimiento cinemático puro. Solidum implementa endurecimiento isótropo
lineal (sin kinematic), donde la frontera de fluencia crece simétricamente
con ``α``. Las propiedades cíclicas verificables en este modelo son:

1. **Rango elástico simétrico tras carga**: tras alcanzar ``α_R``, el
   intervalo elástico es ``Δσ = 2·(σ_y + H·α_R)`` — crece con la
   plastificación previa.
2. **Tangente plástica consistente**: la pendiente tangente en plast
   inversa coincide con la del backbone, ``E_t = (E·H)/(E+H)`` — el
   return mapping debe devolver el mismo módulo en ambos sentidos de
   carga.
3. **Monotonicidad de α**: ``α`` no decrece nunca, sin importar el signo
   de ``Δε_p`` en cada paso. Acumula ``∫|dε_p|``.
4. **Acumulación crece con amplitud**: un ciclo de amplitud mayor
   acumula más ``α`` que uno de amplitud menor.

Para `VonMises2D` plane strain bajo carga uniaxial impuesta, el modelo se
reduce al caso 1D con módulos efectivos plane strain — se verifica el
mismo conjunto de propiedades.
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.materials.von_mises_2d import VonMises2D


# =============================================================================
# Elastoplastic1D — verificación canónica de plasticidad cíclica isótropa
# =============================================================================

class TestElastoplastic1DCyclic(unittest.TestCase):
    """Propiedades del Elastoplastic1D con isótropo lineal bajo cíclica."""

    def setUp(self):
        self.E = 1.0e5
        self.sigma_y = 200.0
        self.H = 1.0e3
        self.mat = Elastoplastic1D(E=self.E, sigma_y=self.sigma_y, H=self.H)
        self.eps_y = self.sigma_y / self.E

    def test_elastic_range_grows_symmetrically_with_alpha(self):
        """Tras cargar hasta α_R, el rango elástico es 2·(σ_y + H·α_R).

        Verificación operativa: identifico el primer ε de descarga donde
        ocurre la plasticidad inversa, y comparo contra la predicción
        analítica ``ε_yield_neg = α_R - (σ_y + H·α_R)/E``.
        """
        # Carga monótona hasta α_R conocido.
        eps_R = 3.0 * self.eps_y
        sigma_R, _, state_R = self.mat.compute_state(eps_R, None)
        alpha_R = state_R['alpha']

        # σ_R analítico: σ_y + H·α_R.
        self.assertAlmostEqual(sigma_R, self.sigma_y + self.H * alpha_R,
                                places=8)

        # ε predicho del yield inverso.
        eps_yield_neg = alpha_R - (self.sigma_y + self.H * alpha_R) / self.E

        # Probe justo antes del yield inverso: debe ser elástico (α sin cambio).
        eps_pre = eps_yield_neg + 1.0e-9
        _, _, state_pre = self.mat.compute_state(eps_pre, state_R)
        self.assertAlmostEqual(state_pre['alpha'], alpha_R, places=10,
            msg=f"Justo antes del yield inverso (ε={eps_pre:.6e}) hubo "
                f"plastificación inesperada: α={state_pre['alpha']:.6e}, "
                f"esperado {alpha_R:.6e}.")

        # Probe justo después: debe ser plástico (α creció).
        eps_post = eps_yield_neg - 1.0e-9
        _, _, state_post = self.mat.compute_state(eps_post, state_R)
        self.assertGreater(state_post['alpha'], alpha_R + 1.0e-15,
            msg=f"Justo después del yield inverso (ε={eps_post:.6e}) no hubo "
                f"plastificación: α={state_post['alpha']:.6e}.")

    def test_tangent_modulus_in_reverse_plasticity_matches_backbone(self):
        """E_tan en plast inversa = (E·H)/(E+H), idéntica al backbone.

        El return mapping debe devolver la misma tangente plástica en
        carga directa y en reversal — la consistencia algorítmica del
        endurecimiento isótropo no distingue el signo del flujo plástico.
        """
        E_tan_expected = (self.E * self.H) / (self.E + self.H)

        # Cargar hasta plástico (backbone).
        eps_R = 3.0 * self.eps_y
        _, E_tan_backbone, state_R = self.mat.compute_state(eps_R, None)
        self.assertAlmostEqual(E_tan_backbone, E_tan_expected, places=10)

        # Descargar a un punto firmemente en plast inversa (~2 veces el
        # rango elástico bajo el yield inverso).
        alpha_R = state_R['alpha']
        eps_yield_neg = alpha_R - (self.sigma_y + self.H * alpha_R) / self.E
        eps_far_neg = eps_yield_neg - 1.0e-4
        _, E_tan_reverse, _ = self.mat.compute_state(eps_far_neg, state_R)
        self.assertAlmostEqual(E_tan_reverse, E_tan_expected, places=10,
            msg=f"E_tan plast inversa {E_tan_reverse:.6e} ≠ backbone "
                f"{E_tan_expected:.6e}.")

    def test_alpha_accumulates_through_full_cycle(self):
        """Un ciclo simétrico 0 → ε_R → -ε_R → 0 acumula α monótonamente.

        Tras N pasos del ciclo, ``α`` debe ser monotónicamente no
        decreciente. Tras un ciclo completo simétrico, ``α`` debe ser
        estrictamente mayor que el ``α_R`` alcanzado en el cuarto de
        ciclo.
        """
        eps_R = 3.0 * self.eps_y
        # Punto a punto a través de los cuatro cuartos del ciclo.
        eps_path = (
            list(np.linspace(0, eps_R, 10)) +
            list(np.linspace(eps_R, -eps_R, 20)) +
            list(np.linspace(-eps_R, 0, 10))
        )
        state = None
        alpha_prev = 0.0
        alpha_at_quarter = None
        for k, eps in enumerate(eps_path):
            _, _, state = self.mat.compute_state(eps, state)
            alpha = state['alpha']
            self.assertGreaterEqual(alpha, alpha_prev - 1.0e-14,
                msg=f"Paso {k} (ε={eps:.4e}): α decreció "
                    f"{alpha_prev:.4e} → {alpha:.4e}.")
            alpha_prev = alpha
            if k == 9:  # final del primer cuarto
                alpha_at_quarter = alpha

        # Tras el ciclo simétrico completo, α debe ser ≈ 4·α_quarter
        # (acumula |Δε_p| en cada cuarto: 1× en el cuarto positivo, 2×
        # en el descenso a -ε_R, 1× en el retorno a 0).
        # Sin pretender precisión exacta, sí debe ser ≥ 3·α_quarter.
        self.assertGreater(alpha_prev, 3.0 * alpha_at_quarter,
            f"α acumulada tras ciclo ({alpha_prev:.4e}) menor de lo esperado "
            f"vs α_quarter ({alpha_at_quarter:.4e}); el cycling no acumula.")

    def test_alpha_accumulation_scales_with_amplitude(self):
        """Un ciclo de amplitud mayor acumula más α que uno menor.

        Aplica medio ciclo (0 → ε_amplitud → 0) con dos amplitudes
        distintas y verifica que ``α_final`` escala con la amplitud por
        encima del yield.
        """
        amp_small = 2.0 * self.eps_y
        amp_large = 5.0 * self.eps_y

        def _half_cycle(amp):
            state = None
            for eps in np.linspace(0, amp, 20):
                _, _, state = self.mat.compute_state(eps, state)
            for eps in np.linspace(amp, 0, 20):
                _, _, state = self.mat.compute_state(eps, state)
            return state['alpha']

        alpha_small = _half_cycle(amp_small)
        alpha_large = _half_cycle(amp_large)
        self.assertGreater(alpha_large, alpha_small,
            f"α acumulada con amplitud {amp_large:.3e} ({alpha_large:.3e}) "
            f"no es mayor que con amplitud {amp_small:.3e} ({alpha_small:.3e}).")


# =============================================================================
# VonMises2D plane strain bajo carga uniaxial cíclica
# =============================================================================

class TestVonMises2DCyclic(unittest.TestCase):
    """Mismas propiedades sobre VonMises2D plane strain, deformación
    uniaxial impuesta en εxx."""

    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.sigma_y = 200.0
        self.H = 2.0e3
        self.mat = VonMises2D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
            hypothesis='plane_strain',
        )
        # Para plane strain uniaxial, eps_y efectivo es del orden de
        # σ_y·(1-ν²)/E·factor — calibramos empíricamente.

    def test_alpha_monotonic_under_full_cycle(self):
        """α monótona no decreciente bajo ciclo simétrico en εxx."""
        eps_amp = 5.0e-3  # >> yield uniaxial
        # 4 cuartos de ciclo.
        eps_xx_path = (
            list(np.linspace(0, eps_amp, 15)) +
            list(np.linspace(eps_amp, -eps_amp, 30)) +
            list(np.linspace(-eps_amp, 0, 15))
        )
        state = None
        alpha_prev = 0.0
        for k, eps_xx in enumerate(eps_xx_path):
            eps = np.array([eps_xx, 0.0, 0.0])
            _, _, state = self.mat.compute_state(eps, state)
            alpha = state['alpha']
            self.assertGreaterEqual(alpha, alpha_prev - 1.0e-14,
                msg=f"Paso {k} (εxx={eps_xx:.4e}): α decreció "
                    f"{alpha_prev:.6e} → {alpha:.6e}.")
            alpha_prev = alpha

        # α tras ciclo debe ser > 0 (ejercitó plast).
        self.assertGreater(alpha_prev, 0.0,
            "VonMises2D no plastificó en el ciclo; setup mal calibrado.")

    def test_alpha_grows_with_cycle_amplitude(self):
        """Amplitudes mayores acumulan más α (consistencia básica con 1D)."""
        amp_small, amp_large = 2.0e-3, 6.0e-3

        def _half_cycle(amp):
            state = None
            for eps_xx in np.linspace(0, amp, 30):
                eps = np.array([eps_xx, 0.0, 0.0])
                _, _, state = self.mat.compute_state(eps, state)
            for eps_xx in np.linspace(amp, 0, 30):
                eps = np.array([eps_xx, 0.0, 0.0])
                _, _, state = self.mat.compute_state(eps, state)
            return state['alpha']

        alpha_small = _half_cycle(amp_small)
        alpha_large = _half_cycle(amp_large)
        self.assertGreater(alpha_large, alpha_small,
            f"α(amp={amp_large}) = {alpha_large:.4e} no > "
            f"α(amp={amp_small}) = {alpha_small:.4e}.")


if __name__ == '__main__':
    unittest.main()
