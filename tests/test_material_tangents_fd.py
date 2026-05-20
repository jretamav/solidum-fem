"""Validación numérica de tangentes algorítmicas por diferencias finitas
(auditoría H-3.5).

Para materiales con return mapping iterativo (``VonMises2D`` plane strain
y plane stress, ``DruckerPrager2D`` rama regular), la tangente
algorítmica consistente ``C_alg`` debe coincidir con
``dσ/dε`` evaluado por diferencias finitas centradas sobre el mismo
``compute_state``. Coincidencia ≈ FD valida la derivación analítica
de la corrección plástica.

Los tests existentes sólo verifican convergencia indirecta del Newton
global; este archivo cierra la deuda H-3.5 con un check directo de
``C_alg ≈ (σ(ε+δε)−σ(ε−δε))/(2·δε)`` en estado plástico activo.

Rama ``apex`` de Drucker-Prager: no se incluye porque es no smooth
(transición regular↔apex no derivable en sentido clásico); su validación
queda en los tests de convergencia indirecta.
"""
from __future__ import annotations

import unittest

import numpy as np

import solidum  # autodiscover

from solidum.materials.drucker_prager_2d import DruckerPrager2D
from solidum.materials.von_mises_2d import VonMises2D


# Step de diferencia finita: suficientemente pequeño para que la FD
# atrape la pendiente sin caer en ruido numérico de máquina (con
# strains O(1e-2), δε=1e-7 da ~5 cifras de precisión).
_FD_STEP = 1.0e-7

# Tolerancia para comparar C_alg con la FD: ~5 cifras significativas
# permite la curvatura O(δε) del residual no lineal y el ruido de la
# corrección plástica del return mapping.
_TANGENT_RTOL = 5.0e-4


def _numerical_tangent_central(material, strain: np.ndarray,
                                state_vars,
                                step: float = _FD_STEP) -> np.ndarray:
    """``dσ_i/dε_j`` por diferencias finitas centradas.

    Importante: ``compute_state`` muta su ``state_vars`` (trial). Para
    evitar contaminación entre evaluaciones de perturbación, pasamos
    **copias independientes** del state cada vez. La tangente FD así
    obtenida es la "algorítmica" (con corrección plástica completa) —
    la misma que ``C_alg`` debe reproducir.
    """
    n = strain.size
    C_fd = np.zeros((n, n))
    for j in range(n):
        eps_plus = strain.copy(); eps_plus[j] += step
        eps_minus = strain.copy(); eps_minus[j] -= step
        sigma_plus, _, _ = material.compute_state(
            eps_plus, state_vars=dict(state_vars) if state_vars else None,
        )
        sigma_minus, _, _ = material.compute_state(
            eps_minus, state_vars=dict(state_vars) if state_vars else None,
        )
        C_fd[:, j] = (sigma_plus - sigma_minus) / (2.0 * step)
    return C_fd


class TestVonMises2DPlaneStrainTangent(unittest.TestCase):
    """``VonMises2D`` plane strain en carga plástica activa: ``C_alg ≈ C_fd``."""

    def test_tangent_matches_finite_difference_in_plastic_regime(self):
        mat = VonMises2D(
            E=200.0, nu=0.25, sigma_y=1.0, H=10.0,
            hypothesis="plane_strain",
        )
        # Estado plástico tras un paso uniaxial: ε > ε_y para que el
        # return mapping corrija y la tangente sea algorítmica (no elástica).
        strain = np.array([0.02, -0.005, 0.0])
        state = None  # estado virgen — el paso entra desde elástico.
        sigma, C_alg, state_new = mat.compute_state(strain, state_vars=state)

        # Sanity: debe haber plastificado (α_new > 0).
        self.assertGreater(state_new["alpha"], 0.0,
                           "el estado debe ser plástico para esta prueba")

        # Re-evaluación FD desde el mismo estado de entrada.
        C_fd = _numerical_tangent_central(mat, strain, state_vars=None)
        np.testing.assert_allclose(
            C_alg, C_fd,
            rtol=_TANGENT_RTOL,
            err_msg="VonMises2D plane strain: C_alg ≠ C_fd (regresión H-3.5)",
        )


class TestVonMises2DPlaneStressTangent(unittest.TestCase):
    """``VonMises2D`` plane stress (Simó-Hughes §3.4.1): la tangente
    cerrada incluye la corrección por ``dα/dΔγ`` no constante. La FD
    debe reproducirla."""

    def test_tangent_matches_finite_difference_in_plastic_regime(self):
        mat = VonMises2D(
            E=200.0, nu=0.25, sigma_y=1.0, H=10.0,
            hypothesis="plane_stress",
        )
        strain = np.array([0.02, -0.005, 0.0])
        sigma, C_alg, state_new = mat.compute_state(strain, state_vars=None)
        self.assertGreater(state_new["alpha"], 0.0,
                           "el estado debe ser plástico para esta prueba")

        C_fd = _numerical_tangent_central(mat, strain, state_vars=None)
        np.testing.assert_allclose(
            C_alg, C_fd,
            rtol=_TANGENT_RTOL,
            err_msg="VonMises2D plane stress: C_alg ≠ C_fd (regresión H-3.5)",
        )


class TestDruckerPrager2DRegularTangent(unittest.TestCase):
    """``DruckerPrager2D`` rama regular: ``C_alg = K v⊗v + 2G(1−β) I_dev
    + 4Gβ n̂⊗n̂ − b_g⊗b_f/A``. El término ``+4Gβ n̂⊗n̂`` captura la
    dependencia de la dirección de flujo con ε vía ``s_trial`` — su
    verificación FD es lo que la auditoría reclamaba."""

    def test_tangent_matches_finite_difference_regular_branch(self):
        # Asociado (ψ = φ) para que C_alg sea simétrica — facilita
        # comparación; el caso no asociado se cubre por sentido común
        # (mismo flujo de derivación, sólo cambia b_g ≠ b_f).
        mat = DruckerPrager2D(
            E=200.0, nu=0.25, cohesion=1.0,
            phi_deg=30.0, psi_deg=30.0, H=10.0,
        )
        # Tracción + cortante moderado para entrar en la rama regular
        # (alejarse del ápice hidrostático y de la frontera elástica).
        strain = np.array([0.03, 0.005, 0.01])
        sigma, C_alg, state_new = mat.compute_state(strain, state_vars=None)
        self.assertGreater(state_new["alpha"], 0.0,
                           "el estado debe ser plástico para esta prueba")

        C_fd = _numerical_tangent_central(mat, strain, state_vars=None)
        np.testing.assert_allclose(
            C_alg, C_fd,
            rtol=_TANGENT_RTOL,
            err_msg="DruckerPrager2D regular: C_alg ≠ C_fd (regresión H-3.5)",
        )


if __name__ == "__main__":
    unittest.main()
