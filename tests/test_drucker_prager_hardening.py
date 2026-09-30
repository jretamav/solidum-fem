"""Endurecimiento de Drucker-Prager en la parametrización de la referencia.

De Souza Neto, Perić y Owen (2008, §8.3): la cohesión se endurece con la
deformación plástica equivalente, ``c(ε̄ᵖ) = c₀ + H·ε̄ᵖ``, con ``ε̄ᵖ = ξ·γ`` y
el radio del cono ``k = ξ·c``. El parámetro ``H`` que escribe el usuario es
``dc/dε̄ᵖ`` y el return mapping trabaja con ``H_k = dk/dγ = ξ²·H``. Hasta el
2026-09-29 ``H`` era directamente ``H_k`` (hallazgo clas-6 de la auditoría;
cambio aprobado por el usuario).
"""
import math

import numpy as np
import pytest

from solidum.materials.drucker_prager_2d import DruckerPrager2D
from solidum.materials.drucker_prager_3d import DruckerPrager3D

E, NU, C0, PHI, H_C = 30.0e9, 0.2, 2.0e6, 30.0, 5.0e8


def _xi(variant):
    """ξ de la referencia para cada calibración (k₀ = ξ·c₀)."""
    phi = math.radians(PHI)
    if variant == "plane_strain_matched":
        return 3.0 / math.sqrt(9.0 + 12.0 * math.tan(phi) ** 2)
    s = math.sin(phi)
    signo = -1.0 if variant == "outer_cone" else 1.0
    return 6.0 * math.cos(phi) / (math.sqrt(3.0) * (3.0 + signo * s))


CASOS = [
    (DruckerPrager2D, "plane_strain_matched"),
    (DruckerPrager2D, "outer_cone"),
    (DruckerPrager2D, "inner_cone"),
    (DruckerPrager3D, "outer_cone"),
    (DruckerPrager3D, "inner_cone"),
]


@pytest.mark.parametrize("cls,variant", CASOS)
def test_radio_del_cono_con_cohesion_endurecida(cls, variant):
    """k(γ) = ξ·c(ε̄ᵖ) con c = c₀ + H·ε̄ᵖ y ε̄ᵖ = ξ·γ."""
    m = cls(E, NU, C0, PHI, H=H_C, variant=variant)
    xi = _xi(variant)
    assert m.k0 == pytest.approx(xi * C0, rel=1e-14)
    for gamma in (0.0, 1e-4, 3e-3):
        c = C0 + H_C * (xi * gamma)
        assert m.admissibility_scale({"alpha": gamma}) == pytest.approx(xi * c, rel=1e-14)


@pytest.mark.parametrize("cls,variant", CASOS)
def test_el_return_mapping_usa_el_modulo_convertido(cls, variant):
    """Tras un paso plástico regular, el estado cae sobre la superficie
    endurecida con H_k = ξ²·H: f(σ, γ) = √J₂ + η_f·I₁ − (k₀ + ξ²·H·γ) = 0."""
    m = cls(E, NU, C0, PHI, H=H_C, variant=variant)
    xi = _xi(variant)
    if cls is DruckerPrager2D:
        eps = np.array([4.0e-4, -1.0e-4, 6.0e-4])
    else:
        eps = np.array([4.0e-4, -1.0e-4, 0.5e-4, 6.0e-4, 1.0e-4, -2.0e-4])
    sigma, _, st = m.compute_state(eps, None)
    gamma = st["alpha"]
    assert gamma > 0.0, "el paso debe ser plástico"
    if cls is DruckerPrager2D:
        s_zz = m.out_of_plane_stress(sigma, st)
        sxx, syy, sxy = sigma
        full = np.array([[sxx, sxy, 0.0], [sxy, syy, 0.0], [0.0, 0.0, s_zz]])
    else:
        sxx, syy, szz, sxy, syz, sxz = sigma
        full = np.array([[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]])
    I1 = np.trace(full)
    dev = full - I1 / 3.0 * np.eye(3)
    sqrt_J2 = math.sqrt(0.5 * np.sum(dev * dev))
    k = m.k0 + xi ** 2 * H_C * gamma
    assert sqrt_J2 + m.eta_f * I1 - k == pytest.approx(0.0, abs=1e-8 * k)
