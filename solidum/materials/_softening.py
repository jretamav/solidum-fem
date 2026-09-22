"""Helpers de softening para materiales de daño escalar (auditoría H-3.4).

Solidum tiene hoy tres materiales con daño escalar y ley de softening
exponencial sobre la variable histórica ``κ``: :class:`IsotropicDamage1D`,
:class:`IsotropicDamage2D` e :class:`IsotropicDamage3D`. La fórmula y el
cap a ``DAMAGE_MAX`` son idénticos en todos; sólo difieren en el cómputo
de la deformación equivalente ``ε_eq`` y en la tangente algorítmica. Este
módulo centraliza la fórmula del daño y deja a cada material la
responsabilidad del resto.

Si en el futuro entra un cuarto modelo (Mazars con split tracción/compresión,
o un daño con softening lineal estilo cohesivo), se añade su variante aquí.

Convención: ``kappa`` y ``kappa_0`` son escalares; la función NO conoce
nada de Voigt ni de esfuerzo. El cap ``DAMAGE_MAX`` aplica directamente
al valor de ``ω`` devuelto (semántica de daño continuo — distinto del
cohesivo, donde el cap sólo afecta a la tangente; ver memoria
``feedback_damage_max_cohesivos.md``).

Compilada con Numba (ADR 0014) para que los kernels 2D/3D de daño, que
también son ``@njit``, la llamen directamente; desde Python (daño 1D) se
invoca igual que antes.
"""
from __future__ import annotations

import math

from numba import njit

from solidum.constants import DAMAGE_MAX


@njit(cache=True)
def evaluate_exponential_damage(kappa, kappa_0, alpha, cap=DAMAGE_MAX):
    """Ley de daño exponencial ``ω(κ) = 1 − (κ₀/κ)·exp(−α·(κ−κ₀))``.

    Compartida por ``IsotropicDamage1D/2D/3D`` (Simó-Ju 1987,
    Lemaitre-Chaboche 1990). Saturación en ``cap`` (default
    :data:`solidum.constants.DAMAGE_MAX`).

    Parameters
    ----------
    kappa
        Variable histórica corriente (deformación equivalente máxima).
    kappa_0
        Umbral elástico del material.
    alpha
        Velocidad de degradación (parámetro físico del material).
    cap
        Valor máximo de ``ω`` antes de saturar. Default ``DAMAGE_MAX``.

    Returns
    -------
    omega, saturated
        ``omega ∈ [0, cap]`` y ``saturated`` indicando si se alcanzó el
        cap (útil para que el llamador conmute a tangente secante en la
        rama saturada).
    """
    if kappa <= kappa_0:
        return 0.0, False
    omega = 1.0 - (kappa_0 / kappa) * math.exp(-alpha * (kappa - kappa_0))
    if omega >= cap:
        return cap, True
    return omega, False
