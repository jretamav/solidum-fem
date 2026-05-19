"""Consistencia termodinámica de los materiales no lineales.

Segunda ley (Clausius-Duhem para procesos isotermos): la disipación
interna debe ser **no negativa** en cada paso del proceso de carga.

Para los modelos del catálogo de Fenix, esto se traduce en dos
verificaciones operativas, fáciles de implementar y muy diagnósticas:

1. **Monotonicidad de las variables internas irreversibles** bajo carga
   monótona creciente:
   - Plasticidad: ``α`` (deformación plástica acumulada equivalente) no
     decrece, y ``‖ε_p‖`` no decrece tampoco.
   - Daño: ``κ`` (variable histórica) y ``ω`` (daño) no decrecen.

2. **Positividad de la disipación incremental por paso**:
   - Plasticidad asociada (J2): ``σ : Δε_p ≥ 0`` en cada paso plástico.
   - Daño: ``ψ_e · Δω ≥ 0`` donde ``ψ_e`` es la energía elástica
     instantánea — equivalente a ``Δω ≥ 0`` bajo carga (ψ_e ≥ 0 siempre).

Para modelos no asociados (Drucker-Prager con ``ψ ≠ φ``) la positividad
de ``σ : Δε_p`` no está garantizada por la teoría — sólo la
monotonicidad de las internas. Esto se documenta en el test.

Estas propiedades cazan bugs sutiles que no rompen la convergencia pero
producirían soluciones físicamente inadmisibles (e.g. "auto-curado" del
daño bajo carga, o un return mapping con signo invertido del flujo
plástico).
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic
from fenix.materials.damage_1d import IsotropicDamage1D
from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.materials.drucker_prager_2d import DruckerPrager2D
from fenix.materials.von_mises_2d import VonMises2D


def _voigt_inner(sigma_3: np.ndarray, eps_4: np.ndarray) -> float:
    """`σ : ε` con notación Fenix: σ es [xx, yy, xy] y ε_p es
    [xx, yy, zz, xy_tensorial].

    El producto interno tensorial es:
        σ_xx·ε_xx + σ_yy·ε_yy + σ_zz·ε_zz + 2·σ_xy·ε_xy

    Aquí σ_zz no está en el vector σ (Fenix expone sólo 3 componentes en
    plane strain/stress); su contribución se omite. Para plane strain
    con plasticidad J2 la traza desviadora de ε_p es cero (incompressible
    plastic flow), por lo que ε_p_zz = -(ε_p_xx + ε_p_yy), y σ_zz aporta;
    pero como σ_zz aproxima a ν·(σ_xx+σ_yy) en estado elástico, y el flujo
    plástico es ortogonal al desviador, el aporte de σ_zz·ε_p_zz al
    producto interno completo sigue siendo bien definido y positivo en
    asociada. Para mantener el test simple, lo aproximamos por el inner
    sobre las 3 componentes expuestas, suficiente para detectar signos
    invertidos.
    """
    return (sigma_3[0] * eps_4[0]
            + sigma_3[1] * eps_4[1]
            + 2.0 * sigma_3[2] * eps_4[3])


# =============================================================================
# Plasticidad J2 (asociada): σ:Δε_p ≥ 0 + monotonicidad de α
# =============================================================================

class TestVonMises2DDissipation(unittest.TestCase):

    def _run_monotonic_loading(self, hypothesis: str):
        mat = VonMises2D(
            E=2.0e5, nu=0.3, sigma_y=200.0, H=2.0e3, hypothesis=hypothesis,
        )
        # Trayectoria monótona creciente en ε_xx (de 0 a 5×ε_yield).
        eps_y = 200.0 / 2.0e5  # ≈ 1e-3
        n_steps = 30
        state = None
        sigma_prev = np.zeros(3)
        eps_p_prev = np.zeros(4)
        alpha_prev = 0.0
        n_plastic_steps = 0

        for k in range(1, n_steps + 1):
            eps = np.array([5.0 * eps_y * k / n_steps, 0.0, 0.0])
            sigma, _, state_new = mat.compute_state(eps, state)
            eps_p_new = state_new['eps_p']
            alpha_new = state_new['alpha']

            # Monotonicidad de α (deformación plástica acumulada).
            self.assertGreaterEqual(alpha_new, alpha_prev - 1.0e-14,
                f"{hypothesis} paso {k}: α decreció {alpha_prev:.6e} → {alpha_new:.6e}.")

            # Disipación incremental σ:Δε_p ≥ 0 cuando hay flujo plástico.
            delta_eps_p = eps_p_new - eps_p_prev
            if np.linalg.norm(delta_eps_p) > 1.0e-12:
                n_plastic_steps += 1
                sigma_avg = 0.5 * (sigma + sigma_prev)
                D = _voigt_inner(sigma_avg, delta_eps_p)
                self.assertGreaterEqual(D, -1.0e-10,
                    f"{hypothesis} paso {k}: σ:Δε_p = {D:.3e} < 0 "
                    "(disipación plástica negativa).")

            sigma_prev = sigma.copy()
            eps_p_prev = eps_p_new.copy()
            alpha_prev = float(alpha_new)
            state = state_new

        self.assertGreater(n_plastic_steps, 5,
            f"{hypothesis}: sólo {n_plastic_steps} pasos plásticos; "
            "el test no ejercitó el régimen objetivo.")

    def test_plane_strain(self):
        self._run_monotonic_loading('plane_strain')

    def test_plane_stress(self):
        self._run_monotonic_loading('plane_stress')


# =============================================================================
# Drucker-Prager (no asociado): monotonicidad de α — D ≥ 0 no garantizada
# =============================================================================

class TestDruckerPragerMonotonicity(unittest.TestCase):
    """En plasticidad NO asociada (ψ ≠ φ), la positividad de
    σ:Δε_p no está garantizada por la teoría; sólo la monotonicidad de
    la variable interna histórica ``α``."""

    def test_alpha_monotonic_under_tensile_shear_loading(self):
        mat = DruckerPrager2D(
            E=1.0e7, nu=0.3, cohesion=1.0e3,
            phi_deg=30.0, psi_deg=10.0,  # NO asociado
            H=1.0e5, hypothesis='plane_strain',
        )
        n_steps = 30
        state = None
        alpha_prev = 0.0
        n_plastic_steps = 0
        eps_scale = 1.0e3 / 1.0e7  # cohesion/E

        # Trayectoria con tracción + cortante crecientes — para DP la
        # tracción activa la rama regular más fácilmente que la compresión
        # (η_f·I_1 con I_1 > 0 hace crecer f más rápido). Magnitud
        # suficiente para entrar en plasticidad y mantenerse allí varios
        # pasos consecutivos.
        for k in range(1, n_steps + 1):
            ramp = 50.0 * k / n_steps
            eps = np.array([ramp, 0.3 * ramp, 0.5 * ramp]) * eps_scale
            _, _, state_new = mat.compute_state(eps, state)
            alpha_new = state_new['alpha']

            self.assertGreaterEqual(alpha_new, alpha_prev - 1.0e-14,
                f"DP paso {k}: α decreció {alpha_prev:.6e} → {alpha_new:.6e}.")

            if alpha_new > alpha_prev + 1.0e-14:
                n_plastic_steps += 1

            alpha_prev = float(alpha_new)
            state = state_new

        self.assertGreater(n_plastic_steps, 5,
            f"DP: sólo {n_plastic_steps} pasos plásticos.")


# =============================================================================
# Daño isótropo: κ y ω monótonas bajo carga creciente
# =============================================================================

class TestIsotropicDamageMonotonicity(unittest.TestCase):
    """Para daño: ``κ̇ ≥ 0`` y ``ω̇ ≥ 0`` bajo carga monótona creciente
    (Kuhn-Tucker + irreversibilidad)."""

    def test_damage_1d(self):
        mat = IsotropicDamage1D(E=1.0e5, kappa_0=1.0e-3, alpha=100.0)
        state = None
        kappa_prev = mat.kappa_0
        d_prev = 0.0
        n_active = 0

        for k in range(1, 21):
            eps = 3.0e-3 * k / 20.0  # creciente monotónica
            _, _, state_new = mat.compute_state(eps, state)
            kappa_new = state_new['kappa']
            d_new = state_new['damage']

            self.assertGreaterEqual(kappa_new, kappa_prev - 1.0e-14,
                f"1D paso {k}: κ decreció {kappa_prev:.6e} → {kappa_new:.6e}.")
            self.assertGreaterEqual(d_new, d_prev - 1.0e-14,
                f"1D paso {k}: ω decreció {d_prev:.6e} → {d_new:.6e}.")

            if d_new > d_prev + 1.0e-14:
                n_active += 1

            kappa_prev = float(kappa_new)
            d_prev = float(d_new)
            state = state_new

        self.assertGreater(n_active, 5,
            f"IsotropicDamage1D: sólo {n_active} pasos con daño activo.")

    def test_damage_2d(self):
        mat = IsotropicDamage2D(
            E=2.0e5, nu=0.3, kappa_0=1.0e-4, alpha=100.0,
            hypothesis='plane_stress',
        )
        state = None
        kappa_prev = mat.kappa_0
        d_prev = 0.0
        n_active = 0

        for k in range(1, 21):
            eps = np.array([5.0e-4 * k / 20.0, 0.0, 0.0])
            _, _, state_new = mat.compute_state(eps, state)
            kappa_new = state_new['kappa']
            d_new = state_new['damage']

            self.assertGreaterEqual(kappa_new, kappa_prev - 1.0e-14)
            self.assertGreaterEqual(d_new, d_prev - 1.0e-14)

            if d_new > d_prev + 1.0e-14:
                n_active += 1

            kappa_prev = float(kappa_new)
            d_prev = float(d_new)
            state = state_new

        self.assertGreater(n_active, 5)


# =============================================================================
# Cohesivo: κ y ω monótonas bajo apertura monótona
# =============================================================================

class TestCohesiveMonotonicity(unittest.TestCase):

    def _run(self, softening: str):
        sigma_t0 = 2.0e6; G_f = 100.0; K_e = 5.0e10
        mat = CohesiveDamageIsotropic(
            sigma_t0=sigma_t0, G_f=G_f, K_e=K_e, softening=softening,
        )

        state = None
        kappa_prev = mat.kappa_0
        d_prev = 0.0
        n_active = 0

        # Apertura normal creciente desde 0 hasta 5·κ_0 (lejos del cap
        # residual para que el régimen activo sea claro).
        n_steps = 20
        for k in range(1, n_steps + 1):
            u_n = 5.0 * mat.kappa_0 * k / n_steps
            jump = np.array([u_n, 0.0])
            _, _, state_new = mat.compute_traction(jump, state)
            kappa_new = state_new['kappa']
            d_new = state_new['damage']

            self.assertGreaterEqual(kappa_new, kappa_prev - 1.0e-14,
                f"{softening} paso {k}: κ decreció.")
            self.assertGreaterEqual(d_new, d_prev - 1.0e-14,
                f"{softening} paso {k}: ω decreció.")

            if d_new > d_prev + 1.0e-14:
                n_active += 1

            kappa_prev = float(kappa_new)
            d_prev = float(d_new)
            state = state_new

        self.assertGreater(n_active, 5,
            f"{softening}: sólo {n_active} pasos con daño activo.")

    def test_linear(self):
        self._run('linear')

    def test_exponential(self):
        self._run('exponential')


if __name__ == '__main__':
    unittest.main()
