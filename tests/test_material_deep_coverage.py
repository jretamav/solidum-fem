"""Profundización de cobertura física para materiales no lineales existentes.

Cubre tres regímenes no testados explícitamente en el catálogo previo:

- **DruckerPrager2D — apex, cono abierto y régimen regular**: tres
  trayectorias canónicas que ejercitan cada rama del return mapping.
- **VonMises2D plane stress — biaxial puro**: estado biaxial isótropo
  es puente físico entre uniaxial validado y multiaxial general; la
  frontera de fluencia debe ocurrir en σ = ±σ_y.
- **IsotropicDamage2D — trayectorias no proporcionales**: rotación de
  la dirección principal durante la carga; verifica que κ histórica
  responde a `ε_eq` invariante (no a una componente particular) y
  que ω no decrece bajo carga monótona en `ε_eq`.
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.materials.drucker_prager_2d import DruckerPrager2D
from fenix.materials.von_mises_2d import VonMises2D


# =============================================================================
# C.2 — DruckerPrager2D: apex (tracción), cono abierto (compresión), regular (cortante)
# =============================================================================

class TestDruckerPrager2DRegimes(unittest.TestCase):
    """Tres trayectorias que ejercitan las tres ramas físicas del DP."""

    def setUp(self):
        self.mat = DruckerPrager2D(
            E=1.0e7, nu=0.3, cohesion=1.0e3,
            phi_deg=30.0, psi_deg=10.0, H=1.0e5, hypothesis='plane_strain',
        )

    def test_isotropic_tension_returns_to_apex(self):
        """Bajo tracción isótropa creciente, el material retorna al apex:
        ``√J₂ = 0`` exacto y ``α`` crece con cada incremento."""
        state = None
        alpha_prev = 0.0
        n_apex_steps = 0
        for k in range(1, 11):
            e = 5.0e-4 * k
            eps = np.array([e, e, 0.0])
            sigma, _, state = self.mat.compute_state(eps, state)

            # Desviador 2D (la parte zz no está expuesta; en apex puro la
            # parte expuesta debe colapsar exactamente).
            sigma_h = 0.5 * (sigma[0] + sigma[1])
            sqrt_J2_2d = math.sqrt(
                (sigma[0] - sigma_h) ** 2
                + (sigma[1] - sigma_h) ** 2
                + 2.0 * sigma[2] ** 2
            )
            self.assertAlmostEqual(sqrt_J2_2d, 0.0, places=10,
                msg=f"k={k}: √J₂(2D) = {sqrt_J2_2d:.3e} ≠ 0 en apex.")

            alpha_new = state['alpha']
            self.assertGreater(alpha_new, alpha_prev - 1.0e-14)
            if alpha_new > alpha_prev + 1.0e-14:
                n_apex_steps += 1
            alpha_prev = float(alpha_new)

        self.assertGreater(n_apex_steps, 5,
            f"Sólo {n_apex_steps} pasos activos en apex.")

    def test_isotropic_compression_stays_elastic(self):
        """Bajo compresión isótropa creciente, el cono DP se abre y el
        material permanece elástico (``α = 0`` siempre).

        Esto documenta una limitación física conocida del DP "puro" (sin
        cap): no hay límite a la compresión hidrostática. La incorporación
        de un cap es modeling, no bug.
        """
        state = None
        for k in range(1, 11):
            e = -5.0e-4 * k
            eps = np.array([e, e, 0.0])
            sigma, _, state = self.mat.compute_state(eps, state)
            self.assertAlmostEqual(state['alpha'], 0.0, places=14,
                msg=f"k={k}: α={state['alpha']:.3e} ≠ 0 en compresión "
                    "isótropa — yield espurio en el cono abierto.")
            # σ debe ser lineal con ε (rama elástica): σ_xx = σ_yy.
            self.assertAlmostEqual(sigma[0], sigma[1], places=8)

    def test_pure_shear_activates_regular_branch(self):
        """Cortante puro (presión cero) entra en la rama regular del DP
        sin contribución del término friccional ``η_f·I₁``."""
        state = None
        alpha_prev = 0.0
        n_active = 0
        for k in range(1, 11):
            g = 2.0e-3 * k  # engineering shear
            eps = np.array([0.0, 0.0, g])
            sigma, _, state = self.mat.compute_state(eps, state)
            alpha_new = state['alpha']
            if alpha_new > alpha_prev + 1.0e-14:
                n_active += 1
            alpha_prev = float(alpha_new)

        self.assertGreater(n_active, 5,
            f"Cortante puro no activó plasticidad regular: "
            f"{n_active} pasos activos.")


# =============================================================================
# C.3 — VonMises2D plane stress: biaxial isótropo puro
# =============================================================================

class TestVonMises2DPlaneStressBiaxial(unittest.TestCase):
    """Estado biaxial isótropo σxx = σyy ≠ 0, σxy = 0.

    Función de fluencia plane stress (Voigt):
        f(σ) = σxx² + σyy² − σxx·σyy + 3·σxy² − σy² ≤ 0

    Bajo σxx = σyy = σ y σxy = 0:
        f = σ² − σy²    →    yield en |σ| = σy

    Idéntico al caso uniaxial. El test verifica:
    1. La frontera de fluencia bajo biaxial puro ocurre en σ ≈ σy.
    2. El flujo plástico mantiene la simetría: Δεp_xx = Δεp_yy.
    """

    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.sigma_y = 200.0
        self.H = 2.0e3  # endurecimiento pequeño
        self.mat = VonMises2D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
            hypothesis='plane_stress',
        )

    def test_biaxial_yield_at_sigma_y(self):
        """Carga biaxial isótropa monótona: el yield ocurre cuando
        σxx = σyy ≈ σy."""
        # En carga elástica biaxial isótropa con σxx = σyy = σ y plane stress:
        # σ = E·ε / (1 − ν) (estado biaxial elástico)
        # Por tanto ε_yield = σ_y·(1 − ν)/E.
        eps_yield_predicted = self.sigma_y * (1.0 - self.nu) / self.E

        state = None
        last_alpha = 0.0
        first_plastic_eps = None
        for k in range(1, 31):
            e = 0.05 * eps_yield_predicted * k  # 20 pasos hasta cerca del yield
            eps = np.array([e, e, 0.0])
            sigma, _, state = self.mat.compute_state(eps, state)
            if state['alpha'] > last_alpha + 1.0e-14 and first_plastic_eps is None:
                first_plastic_eps = e
            last_alpha = state['alpha']

        self.assertIsNotNone(first_plastic_eps,
            "Ningún paso plastificó bajo carga biaxial.")
        # El primer ε plástico debe estar cerca del predicho analítico.
        rel_err = abs(first_plastic_eps - eps_yield_predicted) / eps_yield_predicted
        self.assertLess(rel_err, 0.06,
            f"Yield biaxial en ε={first_plastic_eps:.4e}, predicho "
            f"{eps_yield_predicted:.4e} (err rel {rel_err:.2%}).")

    def test_biaxial_plastic_flow_preserves_symmetry(self):
        """Tras plastificar biaxialmente, ε_p mantiene la simetría
        ε_p_xx = ε_p_yy (la dirección del flujo es la normal exterior al
        círculo de fluencia, que en biaxial isótropo apunta radialmente)."""
        eps_yield_predicted = self.sigma_y * (1.0 - self.nu) / self.E
        # Cargar firmemente por encima del yield.
        eps = np.array([3.0 * eps_yield_predicted, 3.0 * eps_yield_predicted, 0.0])
        _, _, state = self.mat.compute_state(eps, None)

        eps_p = state['eps_p']
        # Componentes xx y yy del eps_p (índices 0, 1) deben ser iguales.
        self.assertAlmostEqual(eps_p[0], eps_p[1], places=10,
            msg=f"ε_p simetría rota: ε_p_xx={eps_p[0]:.4e}, "
                f"ε_p_yy={eps_p[1]:.4e}.")
        # Componente cortante xy (índice 3, tensorial) debe ser nula.
        self.assertAlmostEqual(eps_p[3], 0.0, places=10,
            msg=f"ε_p_xy tensorial = {eps_p[3]:.4e} ≠ 0 en carga sin cortante.")


# =============================================================================
# C.4 — IsotropicDamage2D bajo trayectoria no proporcional
# =============================================================================

class TestIsotropicDamage2DNonProportional(unittest.TestCase):
    """Trayectoria de carga no proporcional: la dirección principal de ε
    rota mientras `ε_eq` (escalar invariante histórico) crece o no.

    Para `IsotropicDamage2D` el daño depende SÓLO de la historia de
    `ε_eq = √(εᵀ·M·ε)` (escalar invariante). Bajo rotación de la dirección
    principal:
    - Si `ε_eq` no decrece, κ no debe decrecer (irreversibilidad).
    - Si `ε_eq` permanece constante, κ y ω se congelan.
    - Si `ε_eq` aumenta, κ y ω crecen monotónicamente.
    """

    def setUp(self):
        self.mat = IsotropicDamage2D(
            E=2.0e5, nu=0.3, kappa_0=1.0e-4, alpha=100.0,
            hypothesis='plane_stress',
        )

    def _eps_eq(self, eps: np.ndarray) -> float:
        """ε_eq como lo computa el material: ‖ε‖_M con M=diag(1,1,½) para
        Voigt con engineering shear."""
        return math.sqrt(eps[0] ** 2 + eps[1] ** 2 + 0.5 * eps[2] ** 2)

    def test_rotation_at_constant_eps_eq_freezes_damage(self):
        """Trayectoria circular en (εxx, εyy) con ε_eq constante: κ y ω
        no deben cambiar tras el primer punto de la trayectoria.

        Trayectoria: ε = ε_eq·(cos θ, sin θ, 0). ε_eq invariante.
        """
        eps_eq_target = 3.0 * self.mat.kappa_0
        # Anclar el estado en el primer punto θ=0.
        eps0 = np.array([eps_eq_target, 0.0, 0.0])
        _, _, state = self.mat.compute_state(eps0, None)
        kappa0 = state['kappa']
        d0 = state['damage']

        # Recorrer la circunferencia ε_eq constante en (εxx, εyy).
        for theta_deg in (30.0, 60.0, 90.0, 120.0, 180.0, 270.0, 359.0):
            theta = math.radians(theta_deg)
            eps = np.array([
                eps_eq_target * math.cos(theta),
                eps_eq_target * math.sin(theta),
                0.0,
            ])
            # Verificar que ε_eq es invariante.
            self.assertAlmostEqual(self._eps_eq(eps), eps_eq_target, places=12)
            _, _, state_new = self.mat.compute_state(eps, state)
            # κ no debe crecer (ε_eq no aumenta).
            self.assertAlmostEqual(state_new['kappa'], kappa0, places=12,
                msg=f"θ={theta_deg}°: κ cambió {kappa0:.6e} → {state_new['kappa']:.6e}.")
            self.assertAlmostEqual(state_new['damage'], d0, places=12,
                msg=f"θ={theta_deg}°: ω cambió {d0:.6e} → {state_new['damage']:.6e}.")

    def test_non_proportional_path_with_eps_eq_increasing(self):
        """Trayectoria 1: ε se mueve de (ε_a, 0, 0) a (0, ε_b, 0) pasando
        por (ε_a, ε_b, 0). Verifica que κ refleja el máximo histórico
        del ε_eq corriente, no el de las componentes individuales."""
        eps_a = 2.0 * self.mat.kappa_0
        eps_b = 3.0 * self.mat.kappa_0
        state = None

        # Paso 1: cargar en εxx puro.
        eps1 = np.array([eps_a, 0.0, 0.0])
        _, _, state = self.mat.compute_state(eps1, state)
        kappa1 = state['kappa']
        self.assertAlmostEqual(kappa1, self._eps_eq(eps1), places=10)

        # Paso 2: añadir εyy manteniendo εxx.
        eps2 = np.array([eps_a, eps_b, 0.0])
        _, _, state = self.mat.compute_state(eps2, state)
        kappa2 = state['kappa']
        # ε_eq creció (ε_eq² = εxx² + εyy²); κ debe seguirlo.
        self.assertAlmostEqual(kappa2, self._eps_eq(eps2), places=10)
        self.assertGreater(kappa2, kappa1)

        # Paso 3: retirar εxx, quedar sólo εyy.
        eps3 = np.array([0.0, eps_b, 0.0])
        _, _, state = self.mat.compute_state(eps3, state)
        kappa3 = state['kappa']
        # ε_eq ahora es eps_b < ε_eq(eps2). κ NO debe decrecer.
        self.assertAlmostEqual(kappa3, kappa2, places=12,
            msg=f"κ decreció bajo descarga parcial: {kappa2:.6e} → {kappa3:.6e}.")

    def test_principal_direction_rotation_with_growth(self):
        """ε_eq crece mientras la dirección principal rota: κ y ω crecen
        monotónicamente junto con ε_eq."""
        state = None
        kappa_prev = self.mat.kappa_0
        d_prev = 0.0
        n_active = 0

        # Rotación de la dirección principal con magnitud creciente.
        for k in range(1, 16):
            theta = math.radians(8.0 * k)  # rota 8° por paso
            magnitude = (1.5 + 0.5 * k) * self.mat.kappa_0
            eps = np.array([
                magnitude * math.cos(theta),
                magnitude * math.sin(theta),
                0.0,
            ])
            _, _, state = self.mat.compute_state(eps, state)

            self.assertGreaterEqual(state['kappa'], kappa_prev - 1.0e-14)
            self.assertGreaterEqual(state['damage'], d_prev - 1.0e-14)
            if state['damage'] > d_prev + 1.0e-14:
                n_active += 1

            kappa_prev = float(state['kappa'])
            d_prev = float(state['damage'])

        self.assertGreater(n_active, 8,
            f"Sólo {n_active} pasos con daño activo bajo rotación + crecimiento.")


if __name__ == '__main__':
    unittest.main()
