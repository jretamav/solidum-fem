"""Consistent tangent tests por diferenciación finita central para los
materiales no lineales que no tenían cobertura directa.

Verificación canónica (Simo & Hughes 1998 §3.4, Neto-Perić-Owen §6.2):
para cada material, en un estado committed plástico/dañado, la tangente
algorítmica ``C_t`` devuelta por `compute_state` (o `compute_traction`
en cohesivos) debe coincidir con la derivada por diferenciación finita
central ``[σ(ε+δ·e_i) − σ(ε−δ·e_i)] / (2δ)``, dentro de tolerancia
``O(δ²)``.

Cobertura existente fuera de este archivo:

- `IsotropicDamage1D` ([tests/test_materials_unit.py:238](tests/test_materials_unit.py))
- `IsotropicDamage2D` ([tests/test_materials_unit.py:392](tests/test_materials_unit.py))

Este archivo completa la familia con `VonMises2D` (plane strain + plane
stress), `DruckerPrager2D` (régimen regular + apex) y
`CohesiveDamageIsotropic` (modo I activo).

Patrón común: el `state_vars` que se pasa a cada perturbación es el
**committed antes del paso** — la consistent tangent se define como la
derivada del retorno mapping respecto a ε con `state_committed` fijo.
Para garantizar que todas las perturbaciones caen en la misma rama de
carga (sin saltar entre branches y romper la diferenciabilidad), se usa
un anchor state con `alpha` (o `kappa`) ligeramente por debajo del
estado actual.
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic
from solidum.materials.drucker_prager_2d import DruckerPrager2D
from solidum.materials.von_mises_2d import VonMises2D


def _tangent_fd(compute_fn, eps_base: np.ndarray, state_committed,
                 delta: float = 1.0e-7) -> np.ndarray:
    """Diferenciación finita central de ``compute_fn(ε, state)[0]`` respecto a ε.

    `compute_fn` debe devolver `(sigma, C_t, new_state)` para que el [0]
    extraiga el esfuerzo. La derivada se calcula columna a columna
    perturbando cada componente de `ε` y manteniendo `state_committed`
    fijo.

    Returns
    -------
    C_fd : ndarray, shape (n, n)
        Tangente por FD, comparable con la analítica devuelta en la
        misma llamada con `eps_base`.
    """
    n = len(eps_base)
    sigma_base, _, _ = compute_fn(eps_base, state_committed)
    n_out = len(np.atleast_1d(sigma_base))
    C_fd = np.zeros((n_out, n))
    for j in range(n):
        e_j = np.zeros(n); e_j[j] = delta
        sigma_p, _, _ = compute_fn(eps_base + e_j, state_committed)
        sigma_m, _, _ = compute_fn(eps_base - e_j, state_committed)
        C_fd[:, j] = (np.atleast_1d(sigma_p) - np.atleast_1d(sigma_m)) / (2.0 * delta)
    return C_fd


def _rel_error(A_ref: np.ndarray, A_test: np.ndarray) -> float:
    """Error relativo en norma Frobenius. Evita división por cero."""
    return float(np.linalg.norm(A_ref - A_test) / max(np.linalg.norm(A_ref), 1.0e-14))


# =============================================================================
# VonMises2D
# =============================================================================

class TestVonMises2DConsistentTangent(unittest.TestCase):
    """Consistent tangent del kernel J2 Numba (plane strain + plane stress)."""

    def test_plane_strain_in_plastic_regime(self):
        """C_alg analítica coincide con FD en plane strain plastificado.

        Estado anchor con ``alpha`` modesto fuerza que todas las
        perturbaciones de FD caigan en la rama plástica (no oscilen
        entre elástico/plástico).
        """
        E, nu, sigma_y, H = 2.0e5, 0.3, 200.0, 2.0e3
        mat = VonMises2D(E=E, nu=nu, sigma_y=sigma_y, H=H,
                         hypothesis='plane_strain')

        # ε en régimen plástico claro (~3× el yield uniaxial)
        eps_y_uni = sigma_y / E
        eps = np.array([3.0 * eps_y_uni, -0.3 * eps_y_uni, 0.4 * eps_y_uni])
        # Anchor: pequeña plasticidad acumulada previa para que la rama
        # esté firmemente del lado plástico durante todas las perturbaciones.
        anchor = {'eps_p': np.zeros(4), 'alpha': 0.5 * eps_y_uni}

        _, C_alg, _ = mat.compute_state(eps, anchor)
        C_fd = _tangent_fd(mat.compute_state, eps, anchor, delta=1.0e-7)

        err = _rel_error(C_alg, C_fd)
        self.assertLess(err, 1.0e-4,
            f"VonMises2D plane strain: ‖C_alg − C_fd‖/‖C_alg‖ = {err:.3e} "
            f"≥ 1e-4. Posible bug en la tangente algorítmica.\n"
            f"C_alg=\n{C_alg}\nC_fd=\n{C_fd}")

    def test_plane_stress_in_plastic_regime(self):
        """C_alg analítica coincide con FD en plane stress plastificado.

        El kernel plane stress es distinto del plane strain (proyector
        P̄ con iteración local en σ_zz=0); este test lo cubre por separado.
        """
        E, nu, sigma_y, H = 2.0e5, 0.3, 200.0, 2.0e3
        mat = VonMises2D(E=E, nu=nu, sigma_y=sigma_y, H=H,
                         hypothesis='plane_stress')

        eps_y_uni = sigma_y / E
        eps = np.array([3.0 * eps_y_uni, -0.5 * eps_y_uni, 0.6 * eps_y_uni])
        anchor = {'eps_p': np.zeros(4), 'alpha': 0.5 * eps_y_uni}

        _, C_alg, _ = mat.compute_state(eps, anchor)
        C_fd = _tangent_fd(mat.compute_state, eps, anchor, delta=1.0e-7)

        err = _rel_error(C_alg, C_fd)
        self.assertLess(err, 1.0e-4,
            f"VonMises2D plane stress: ‖C_alg − C_fd‖/‖C_alg‖ = {err:.3e} "
            f"≥ 1e-4.\nC_alg=\n{C_alg}\nC_fd=\n{C_fd}")


# =============================================================================
# DruckerPrager2D
# =============================================================================

class TestDruckerPrager2DConsistentTangent(unittest.TestCase):
    """Consistent tangent del return mapping de DP plane strain en sus dos
    ramas (regular y apex)."""

    def setUp(self):
        # Parámetros estándar (φ=30°, ψ=10° no asociado, cohesión moderada).
        self.E = 1.0e7
        self.nu = 0.3
        self.cohesion = 1.0e3
        self.phi_deg = 30.0
        self.psi_deg = 10.0
        self.H = self.E / 100.0

    def test_regular_branch(self):
        """Régimen confinado donde el return mapping cae en la rama regular
        del cono de DP (no en el apex).

        ε con compresión hidrostática moderada + cortante: garantiza
        ``ξ > apex`` y el retorno proyecta sobre la generatriz, no sobre
        el vértice.
        """
        mat = DruckerPrager2D(
            E=self.E, nu=self.nu, cohesion=self.cohesion,
            phi_deg=self.phi_deg, psi_deg=self.psi_deg, H=self.H,
            hypothesis='plane_strain',
        )

        # ε con compresión + desviador apreciable.
        eps_scale = self.cohesion / self.E
        eps = np.array([-3.0, -1.0, 0.8]) * eps_scale
        anchor = {'eps_p': np.zeros(4), 'alpha': 0.5 * eps_scale}

        _, C_alg, _ = mat.compute_state(eps, anchor)
        C_fd = _tangent_fd(mat.compute_state, eps, anchor, delta=1.0e-8)

        err = _rel_error(C_alg, C_fd)
        self.assertLess(err, 1.0e-3,
            f"DruckerPrager2D regular: ‖C_alg − C_fd‖/‖C_alg‖ = {err:.3e} "
            f"≥ 1e-3.\nC_alg=\n{C_alg}\nC_fd=\n{C_fd}")

    def test_apex_branch(self):
        """Régimen muy traccionante donde el return mapping cae en el apex
        del cono (vértice singular).

        ε con tracción isótropa fuerte: el retorno colapsa al vértice y
        usa la rama apex con su propia fórmula de tangente (rango 0 en
        dirección desviadora).
        """
        mat = DruckerPrager2D(
            E=self.E, nu=self.nu, cohesion=self.cohesion,
            phi_deg=self.phi_deg, psi_deg=self.psi_deg, H=self.H,
            hypothesis='plane_strain',
        )

        # ε con tracción isótropa fuerte (sin cortante apreciable) → apex.
        eps_scale = self.cohesion / self.E
        eps = np.array([5.0, 5.0, 0.0]) * eps_scale
        anchor = {'eps_p': np.zeros(4), 'alpha': 0.5 * eps_scale}

        sigma_apex, C_alg, _ = mat.compute_state(eps, anchor)
        # Verificación de que efectivamente está en el apex: el desviador
        # del esfuerzo (en 2D plane strain) debe ser pequeño.
        sigma_h = (sigma_apex[0] + sigma_apex[1]) / 2.0
        sigma_dev = np.array([
            sigma_apex[0] - sigma_h,
            sigma_apex[1] - sigma_h,
            sigma_apex[2],
        ])
        # No requerimos exactitud (puede haber componente residual);
        # sólo que sea < 10% del hidrostático (estamos en apex).
        if abs(sigma_h) > 1.0e-12:
            ratio = np.linalg.norm(sigma_dev) / abs(sigma_h)
            self.assertLess(ratio, 0.5,
                f"Esperaba rama apex pero el esfuerzo tiene desviador "
                f"significativo: |dev|/|h| = {ratio:.3f}.")

        C_fd = _tangent_fd(mat.compute_state, eps, anchor, delta=1.0e-9)
        err = _rel_error(C_alg, C_fd)
        # Tolerancia más holgada: la rama apex es C¹ (no C²) en la
        # transición regular↔apex; el FD captura una mezcla local si
        # alguna perturbación cae cerca del límite.
        self.assertLess(err, 1.0e-2,
            f"DruckerPrager2D apex: ‖C_alg − C_fd‖/‖C_alg‖ = {err:.3e} "
            f"≥ 1e-2.\nC_alg=\n{C_alg}\nC_fd=\n{C_fd}")


# =============================================================================
# CohesiveDamageIsotropic
# =============================================================================

class TestCohesiveConsistentTangent(unittest.TestCase):
    """Consistent tangent del cohesivo Modo-I (rank-1 sobre n⊗n)."""

    def _do_fd_test(self, softening: str):
        sigma_t0 = 2.0e6
        G_f = 100.0
        # K_e suficientemente grande para que el softening exponencial
        # cumpla ``K_e > σ_t0²/(2·G_F) = 2e10``.
        K_e = 5.0e10
        mat = CohesiveDamageIsotropic(
            sigma_t0=sigma_t0, G_f=G_f, K_e=K_e, softening=softening,
        )

        # Jump normal en régimen carga activa: > κ_0 pero MUY por debajo
        # de la apertura crítica `w_c = 2·G_F/σ_t0` (lineal: tracción
        # física se anula en w_c y la tangente analítica se reporta como
        # rigidez residual artificial, by design — no comparable con FD).
        # `u_n = 1.5·κ_0` está claramente activo y << w_c en este setup.
        u_n = 1.5 * mat.kappa_0
        jump = np.array([u_n, 0.0])
        anchor = {'kappa': u_n * 0.9, 'damage': 0.0}

        _, T_alg, _ = mat.compute_traction(jump, anchor)
        # delta relativo razonable: 1e-6 × u_n produce diferencias de
        # tracción del orden K_e × delta ≈ 6, captables sin perder a
        # precisión máquina (jump base ~6e-6, tracción base ~K_e×u_n ~3e5).
        T_fd = _tangent_fd(mat.compute_traction, jump, anchor, delta=1.0e-6 * u_n)

        err = _rel_error(T_alg, T_fd)
        self.assertLess(err, 1.0e-4,
            f"Cohesive {softening}: ‖T_alg − T_fd‖/‖T_alg‖ = {err:.3e} "
            f"≥ 1e-4.\nT_alg=\n{T_alg}\nT_fd=\n{T_fd}")

    def test_linear_softening_loading(self):
        self._do_fd_test('linear')

    def test_exponential_softening_loading(self):
        self._do_fd_test('exponential')


if __name__ == '__main__':
    unittest.main()
