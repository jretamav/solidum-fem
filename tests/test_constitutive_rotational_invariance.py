"""Invariancia rotacional de los constitutivos isótropos del catálogo.

Para un material isótropo, la relación constitutiva debe satisfacer:

    σ(R·ε·Rᵀ) = R·σ(ε)·Rᵀ        (para toda rotación R)

En notación Voigt 2D con `ε_voigt = [εxx, εyy, γxy]` (engineering shear)
y `σ_voigt = [σxx, σyy, σxy]`, esto se traduce a:

    σ(T_ε(θ) · ε) = T_σ(θ) · σ(ε)

con `T_ε` y `T_σ` distintas por la convención de engineering shear:

    T_ε(θ) = [[c², s², cs], [s², c², -cs], [-2cs, 2cs, c²-s²]]
    T_σ(θ) = [[c², s², 2cs], [s², c², -2cs], [-cs, cs, c²-s²]]

El test caza bugs en el ensamblaje del operador `B`, en el return mapping
(factores ½ olvidados en el desviador), en la imposición de hipótesis
plane strain/stress, o en transformaciones cuando algún componente
introduce una dependencia espuria de los ejes globales.

Cubre `Elastic2D` (plane strain + plane stress), `IsotropicDamage2D` y
`VonMises2D` (plane strain + plane stress).
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.materials.damage_2d import IsotropicDamage2D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.von_mises_2d import VonMises2D


def _T_eps(theta: float) -> np.ndarray:
    """Matriz de rotación para ε_voigt = [εxx, εyy, γxy] (engineering shear)."""
    c, s = math.cos(theta), math.sin(theta)
    return np.array([
        [c * c,     s * s,    c * s],
        [s * s,     c * c,   -c * s],
        [-2.0 * c * s, 2.0 * c * s, c * c - s * s],
    ])


def _T_sig(theta: float) -> np.ndarray:
    """Matriz de rotación para σ_voigt = [σxx, σyy, σxy]."""
    c, s = math.cos(theta), math.sin(theta)
    return np.array([
        [c * c,    s * s,     2.0 * c * s],
        [s * s,    c * c,    -2.0 * c * s],
        [-c * s,   c * s,     c * c - s * s],
    ])


THETAS = [math.pi / 6.0, math.pi / 4.0, math.pi / 3.0, math.pi / 2.0]


# =============================================================================
# Elastic2D
# =============================================================================

class TestElastic2DRotationalInvariance(unittest.TestCase):
    """Material puramente elástico — el test debe pasar a precisión máquina."""

    EPS = np.array([2.0e-4, -5.0e-5, 3.0e-4])  # ε no degenerada

    def _check(self, hypothesis: str):
        mat = Elastic2D(E=2.0e5, nu=0.3, hypothesis=hypothesis)
        sigma_ref, _, _ = mat.compute_state(self.EPS)
        for theta in THETAS:
            eps_rot = _T_eps(theta) @ self.EPS
            sigma_rot, _, _ = mat.compute_state(eps_rot)
            sigma_expected = _T_sig(theta) @ sigma_ref
            err = np.linalg.norm(sigma_rot - sigma_expected) / np.linalg.norm(sigma_ref)
            self.assertLess(err, 1.0e-12,
                f"Elastic2D {hypothesis}, θ={math.degrees(theta):.0f}°: "
                f"σ(T_ε·ε) − T_σ·σ(ε) tiene err rel {err:.3e}.")

    def test_plane_strain(self):
        self._check('plane_strain')

    def test_plane_stress(self):
        self._check('plane_stress')


# =============================================================================
# IsotropicDamage2D
# =============================================================================

class TestIsotropicDamage2DRotationalInvariance(unittest.TestCase):
    """El daño depende de `κ = ε_eq`, escalar invariante bajo rotación, por
    lo que la invariancia debe valer también con daño activo."""

    def test_loading_active(self):
        mat = IsotropicDamage2D(
            E=2.0e5, nu=0.3, kappa_0=1.0e-4, alpha=100.0,
        )
        eps = np.array([2.5e-4, 0.7e-4, 0.4e-4])
        # Anchor con κ ligeramente por debajo del ε_eq actual, para
        # garantizar que tanto el caso original como las rotaciones
        # caigan en carga activa (κ_eq es invariante rotacional).
        eps_eq = math.sqrt(eps[0] ** 2 + eps[1] ** 2 + 0.5 * eps[2] ** 2)
        anchor = {'kappa': eps_eq * 0.95, 'damage': 0.0}

        sigma_ref, _, _ = mat.compute_state(eps, anchor)
        for theta in THETAS:
            eps_rot = _T_eps(theta) @ eps
            sigma_rot, _, _ = mat.compute_state(eps_rot, anchor)
            sigma_expected = _T_sig(theta) @ sigma_ref
            err = np.linalg.norm(sigma_rot - sigma_expected) / np.linalg.norm(sigma_ref)
            self.assertLess(err, 1.0e-10,
                f"IsotropicDamage2D θ={math.degrees(theta):.0f}°: err rel {err:.3e}.")


# =============================================================================
# VonMises2D
# =============================================================================

class TestVonMises2DRotationalInvariance(unittest.TestCase):
    """J2 con kernels Numba para plane strain y plane stress."""

    EPS = None  # set per subtest

    def _check(self, hypothesis: str, eps: np.ndarray):
        mat = VonMises2D(
            E=2.0e5, nu=0.3, sigma_y=200.0, H=2.0e3, hypothesis=hypothesis,
        )
        # Anchor virgen (eps_p=0, alpha=0): el return mapping calcula
        # exactamente desde la rama elástica y entra en plasticidad si
        # supera el yield. Pero `eps_p` también es tensorial 2D — si el
        # anchor tuviese eps_p ≠ 0, también habría que rotarlo. Con
        # anchor virgen evitamos esa complicación: la única dependencia
        # de los ejes está en la salida.
        anchor = {'eps_p': np.zeros(4), 'alpha': 0.0}

        sigma_ref, _, _ = mat.compute_state(eps, anchor)
        for theta in THETAS:
            eps_rot = _T_eps(theta) @ eps
            sigma_rot, _, _ = mat.compute_state(eps_rot, anchor)
            sigma_expected = _T_sig(theta) @ sigma_ref
            err = np.linalg.norm(sigma_rot - sigma_expected) / max(
                np.linalg.norm(sigma_ref), 1.0e-14)
            self.assertLess(err, 1.0e-9,
                f"VonMises2D {hypothesis}, θ={math.degrees(theta):.0f}°: "
                f"err rel {err:.3e}.")

    def test_plane_strain_plastic(self):
        # ε de magnitud claramente plástica (>> ε_yield uniaxial = σ_y/E).
        self._check('plane_strain',
                    np.array([3.0e-3, -1.0e-3, 1.5e-3]))

    def test_plane_stress_plastic(self):
        self._check('plane_stress',
                    np.array([3.0e-3, -1.5e-3, 2.0e-3]))


if __name__ == '__main__':
    unittest.main()
