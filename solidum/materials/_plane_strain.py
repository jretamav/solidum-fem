"""Utilidad compartida por los materiales 2D en hipótesis *plane strain*.

En deformación plana la componente ``σ_zz`` no es nula pero tampoco forma
parte del vector Voigt 2D ``[σ_xx, σ_yy, σ_xy]`` que devuelve
``compute_state``. Para recuperarla (invariantes 3D como Von Mises en el
post-proceso) hace falta el modelo constitutivo: aquí vive la fórmula
común a los materiales isótropos con descomposición
``σ = 2G·dev(ε − ε_p) + K·tr(ε − ε_p)·I`` (elástico, J2, Drucker-Prager).

Con ``a = ε − ε_p`` y ``ε_zz = 0`` se tiene ``a_zz = −ε^p_zz`` y, llamando
``λ = K − 2G/3`` al primer parámetro de Lamé,

    σ_xx + σ_yy = 2(λ + G)·(a_xx + a_yy) + 2λ·a_zz
    σ_zz        = λ·(a_xx + a_yy) + (λ + 2G)·a_zz

de donde ``a_xx + a_yy`` se despeja de las componentes en plano conocidas.
Con ``ε_p = 0`` se recupera el resultado elástico ``σ_zz = ν(σ_xx + σ_yy)``.
"""
from __future__ import annotations


def sigma_zz_plane_strain(sigma_xx: float, sigma_yy: float,
                          eps_p_zz: float, K: float, G: float) -> float:
    """``σ_zz`` en deformación plana a partir de las componentes en plano y
    de la deformación plástica fuera del plano.

    Parameters
    ----------
    sigma_xx, sigma_yy
        Componentes normales en plano del esfuerzo (Voigt 2D).
    eps_p_zz
        Componente ``zz`` de la deformación plástica (tensorial). Cero en
        materiales sin plasticidad.
    K, G
        Módulos volumétrico y de cortante del medio.
    """
    lam = K - 2.0 * G / 3.0
    a_zz = -eps_p_zz
    a_sum = (sigma_xx + sigma_yy - 2.0 * lam * a_zz) / (2.0 * (lam + G))
    return lam * a_sum + (lam + 2.0 * G) * a_zz
