"""Leyes de ablandamiento de los materiales de daño escalar (auditoría H-3.4).

Dos leyes ``d(κ)`` sobre la variable histórica ``κ`` (máxima deformación
equivalente):

- **exponencial** en la deformación total, ``ω(κ) = 1 − (κ₀/κ)·exp(−α(κ−κ₀))``
  (Simó-Ju 1987): la de ``IsotropicDamage1D/2D/3D`` y, con ``α(h)`` de la banda
  de grieta, la de ``IsotropicDamageCrackBand2D/3D`` y el daño de tracción de
  ``TensionCompressionDamage2D/3D``;
- **Hordijk** en la apertura de grieta (Cornelissen, Hordijk y Reinhardt 1986;
  spec ``IsotropicDamageCrackBand2D`` §2b), sólo con banda de grieta: la
  apertura ``w = h·d·κ`` es la raíz de ``w/h + κ₀·g(w/w_c) = κ`` y
  ``d = w/(h·κ)``.

:func:`evaluate_damage_law` elige la ley por su código y devuelve el daño y su
pendiente ``dd/dκ``, que entra en la tangente consistente de cada material.
La fórmula de la exponencial y su pendiente son las de siempre, operación por
operación (los resultados con la exponencial no cambian ni en un bit).

Convención: ``kappa`` y ``kappa_0`` son escalares; las funciones NO conocen
nada de Voigt ni de esfuerzo. El tope ``cap`` aplica directamente al valor de
``ω`` devuelto (semántica de daño continuo).

Compiladas con Numba (ADR 0014) para que los kernels 2D/3D de daño, que
también son ``@njit``, las llamen directamente; desde Python (daño 1D) se
invocan igual. Caché de Numba: un cambio aquí no invalida la caché de los
núcleos que las llaman (``damage_2d``, ``damage_3d``, ``_crack_band``,
``_tension_compression``, los cuatro materiales de banda de grieta y
``fixed_smeared_crack_2d``); al
cambiar este módulo hay que tocar también ésos.
"""
from __future__ import annotations

import math

from numba import njit

from solidum.constants import DAMAGE_MAX

# Códigos de la ley (viajan en ``batch_params``).
EXPONENTIAL = 0
HORDIJK = 1
SOFTENING_CODES = {"exponential": EXPONENTIAL, "hordijk": HORDIJK}

# Ley de Hordijk: σ/f_t = g(w/w_c), g(x) = (1 + c₁³x³)·exp(−c₂x) − x·(1 + c₁³)·exp(−c₂)
# (Cornelissen, Hordijk y Reinhardt 1986, ec. 1; Hordijk 1991, ec. 6.1).
HORDIJK_C1 = 3.0
HORDIJK_C2 = 6.93
_E2 = math.exp(-HORDIJK_C2)
_C13 = HORDIJK_C1 ** 3
# Área bajo g en [0, 1], en forma cerrada:
# ∫e^{−cx} = (1 − e^{−c})/c;  ∫x³e^{−cx} = 6/c⁴ − e^{−c}(1/c + 3/c² + 6/c³ + 6/c⁴).
_c = HORDIJK_C2
_AREA = ((1.0 - _E2) / _c
         + _C13 * (6.0 / _c ** 4 - _E2 * (1.0 / _c + 3.0 / _c ** 2 + 6.0 / _c ** 3 + 6.0 / _c ** 4))
         - 0.5 * (1.0 + _C13) * _E2)
# w_c = HORDIJK_KW·G_f/f_t encierra exactamente G_f (5.13606; el de DIANA,
# Jirásek 2011 y Slobbe 2015; ATENA redondea a 5.14). Decisión del usuario,
# 2026-11-24 (Q12 de la spec).
HORDIJK_KW = 1.0 / _AREA
# |g'(0)| = c₂ + (1 + c₁³)·exp(−c₂): pendiente máxima (g es convexa en [0, 1]).
HORDIJK_SLOPE0 = HORDIJK_C2 + (1.0 + _C13) * _E2
# h_max/ℓ_ch = k_w/|g'(0)| = 0.7382 (sin retroceso en el pico).
HORDIJK_H_MAX_FACTOR = HORDIJK_KW / HORDIJK_SLOPE0
# Tope de iteraciones del Newton local (converge en ≤ 12 con h ≤ 0.999·h_max).
HORDIJK_MAX_ITER = 50
del _c


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


@njit(cache=True)
def hordijk_shape(x):
    """``g(x)`` de Hordijk, ``x = w/w_c ∈ [0, 1]``: ``g(0) = 1``, ``g(1) = 0``."""
    return (1.0 + _C13 * x * x * x) * math.exp(-HORDIJK_C2 * x) - x * (1.0 + _C13) * _E2


@njit(cache=True)
def hordijk_shape_slope(x):
    """``g'(x)``, negativa en ``[0, 1]``; ``g'(0) = −HORDIJK_SLOPE0``."""
    ex = math.exp(-HORDIJK_C2 * x)
    return (3.0 * _C13 * x * x * ex - HORDIJK_C2 * (1.0 + _C13 * x * x * x) * ex
            - (1.0 + _C13) * _E2)


@njit(cache=True)
def evaluate_hordijk_damage(kappa, kappa_0, w_c, h, cap):
    """Daño de la ley de Hordijk en la apertura (spec ``IsotropicDamageCrackBand2D`` §2b).

    La apertura ``w`` es la raíz de ``φ(w) = w/h + κ₀·g(w/w_c) − κ`` en
    ``(0, w_c)`` y ``d = w/(h·κ)``; con ``κ ≥ w_c/h`` la grieta está abierta
    del todo y ``d = cap``. ``φ`` es creciente y convexa si ``h < h_max``
    (``HORDIJK_H_MAX_FACTOR·ℓ_ch``, lo comprueba el material al recibir ``h``):
    el Newton arranca a la derecha de la raíz, en
    ``min(w_c, (κ − κ₀)/φ'(0))``, y baja de forma monótona; se detiene cuando
    el iterado deja de decrecer (en el redondeo, sin tolerancia).
    ``h ≤ 0`` (sin sembrar) da el tope, como la exponencial con ``α = ∞``.

    Returns
    -------
    omega, domega_dkappa, saturated
        El daño, su pendiente ``dd/dκ = (1/κ)·(1/(h·φ') − d)`` (0 en el tope)
        y si se alcanzó el tope.
    """
    if kappa <= kappa_0:
        return 0.0, 0.0, False
    if h <= 0.0 or kappa * h >= w_c:
        return cap, 0.0, True
    dphi0 = 1.0 / h - kappa_0 * HORDIJK_SLOPE0 / w_c
    w = w_c
    if dphi0 > 0.0:
        w_tangent = (kappa - kappa_0) / dphi0
        if w_tangent < w:
            w = w_tangent
    for _ in range(HORDIJK_MAX_ITER):
        x = w / w_c
        phi = w / h + kappa_0 * hordijk_shape(x) - kappa
        dphi = 1.0 / h + kappa_0 * hordijk_shape_slope(x) / w_c
        w_new = w - phi / dphi
        if not (w_new < w):
            break
        w = w_new
    omega = w / (h * kappa)
    if omega >= cap:
        return cap, 0.0, True
    dphi = 1.0 / h + kappa_0 * hordijk_shape_slope(w / w_c) / w_c
    return omega, (1.0 / (h * dphi) - omega) / kappa, False


@njit(cache=True)
def evaluate_damage_law(kappa, kappa_0, law, alpha, w_c, h, cap):
    """Daño ``ω(κ)``, su pendiente ``dω/dκ`` y si saturó, según la ley.

    ``law``: ``EXPONENTIAL`` (usa ``alpha``) o ``HORDIJK`` (usa la apertura
    crítica ``w_c`` y la longitud característica ``h``). Con ``κ ≤ κ₀``,
    ``(0, 0, False)``. La pendiente de la exponencial es
    ``(1 − ω)·(1/κ + α)``, la de siempre; sólo cuenta bajo el tope (el
    llamador decide la rama de la tangente).
    """
    if kappa <= kappa_0:
        return 0.0, 0.0, False
    if law == HORDIJK:
        return evaluate_hordijk_damage(kappa, kappa_0, w_c, h, cap)
    omega, saturated = evaluate_exponential_damage(kappa, kappa_0, alpha, cap)
    return omega, (1.0 - omega) * (1.0 / kappa + alpha), saturated


def parse_softening(owner: str, softening) -> int:
    """Código de la ley de ablandamiento (``'exponential'`` o ``'hordijk'``)."""
    code = SOFTENING_CODES.get(softening) if isinstance(softening, str) else None
    if code is None:
        raise ValueError(
            f"{owner}: softening={softening!r} no válida; use 'exponential' (por omisión, "
            f"en la deformación total) o 'hordijk' (en la apertura de grieta, Cornelissen, "
            f"Hordijk y Reinhardt 1986).")
    return code
