"""Regularización por banda de grieta compartida por ``IsotropicDamageCrackBand2D/3D``.

Bažant y Oh (1983): la rama de ablandamiento de cada punto de integración se
ajusta a su longitud característica ``h`` para que el material disipe
``G_f/h`` por unidad de volumen. Con la ley exponencial del daño isótropo
(``solidum.materials._softening``), en tracción uniaxial
``σ = f_t·exp(−α(κ−κ₀))`` y la energía hasta la separación completa es
``g = ½·f_t·κ₀ + f_t/α``; igualarla a ``G_f/h`` da

    α(h) = 1 / (G_f/(h·f_t) − κ₀/2),     h < h_max = 2·E·G_f/f_t².

Ver ``docs/specs/IsotropicDamageCrackBand2D.md`` §2.

Integración IMPL-EX (Oliver, Huespe y Cante 2008; spec §8, ADR 0034): con
``integration='impl-ex'`` la variable histórica se extrapola con los dos
últimos pasos convergidos, ``κ̃ = κ_n + r·(κ_n − κ_{n−1})``, el daño del paso
es ``d̃ = max(d̃_n, d(κ̃))`` y el esfuerzo ``σ = (1 − d̃)·C_e·ε`` con tangente
secante exacta; ``κ_{n+1} = max(κ_n, ε_eq)`` se actualiza implícitamente para
la extrapolación siguiente. ``r = Δτ_{n+1}/Δτ_n`` lo fija el ensamblador por
intento de paso; fuera de un paso (``r`` negativo en el kernel) el daño queda
congelado en ``d̃_n`` y el estado no cambia.

Ley de Hordijk en la apertura (spec §2b, bloque 3b): con
``softening='hordijk'`` el daño es el de ``_softening.evaluate_hordijk_damage``
con la apertura crítica ``w_c = k_w·G_f/f_t`` y la ``h`` del punto, y el límite
de ``h`` pasa a ``h_max = 0.7382·E·G_f/f_t²`` (continuidad en el pico). La
exponencial no cambia.

Por qué un módulo propio: las variantes 2D y 3D comparten la validación de
los parámetros, ``α(h)``, el límite ``h_max`` y el estado inicial por punto
(Reglas.md §1: dos casos reales); cada una conserva su núcleo compilado.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from solidum.constants import CRACK_BAND_DAMAGE_MAX
from solidum.materials._equivalent_strain import EQUIVALENT_STRAIN_CODES
from solidum.materials._softening import (
    EXPONENTIAL,
    HORDIJK,
    HORDIJK_H_MAX_FACTOR,
    HORDIJK_KW,
    evaluate_damage_law,
    parse_softening,
)

# Medidas que admiten las variantes; Rankine por omisión (decisión del
# usuario, 2026-11-11, Q1 de la spec).
CRACK_BAND_EQUIVALENT_STRAINS = tuple(EQUIVALENT_STRAIN_CODES)
DEFAULT_CRACK_BAND_EQUIVALENT_STRAIN = "rankine"

# Integración del daño (spec §7 y §8). El código viaja en ``batch_params``.
INTEGRATION_CODES = {"implicit": 0, "impl-ex": 1}
IMPLICIT, IMPL_EX = 0, 1
DEFAULT_INTEGRATION = "implicit"
# Razón de pseudo-tiempo en el kernel cuando no hay un paso en curso.
NO_STEP_RATIO = -1.0
# Ley de ablandamiento por omisión (usuario, 2026-11-24, Q13 de la spec):
# ningún resultado existente cambia; Hordijk se pide explícitamente.
DEFAULT_SOFTENING = "exponential"


@njit(cache=True)
def crack_band_alpha(h, E, f_t, G_f):
    """``α(h)`` de la banda de grieta; ``h ≤ 0`` (sin sembrar) → ``+∞``.

    Con ``+∞`` el daño salta al tope en cuanto empieza: el llamador lo
    detecta y avisa (``h`` es obligatoria en cuanto hay daño).
    """
    if h <= 0.0:
        return math.inf
    kappa_0 = f_t / E
    return 1.0 / (G_f / (h * f_t) - 0.5 * kappa_0)


@njit(cache=True)
def hordijk_critical_opening(f_t, G_f):
    """``w_c = k_w·G_f/f_t`` de la ley de Hordijk (spec §2b, Q12)."""
    return HORDIJK_KW * G_f / f_t


@njit(cache=True)
def crack_band_damage_slope(kappa, E, kappa_0, tol, cap, law, alpha, w_c, h):
    """Daño ``d(κ)`` de la ley elegida y su pendiente ``dd/dκ``, con el
    umbral de admisibilidad del padre (``E·(κ − κ_0) ≤ tol`` ⇒ ``(0, 0)``) y
    el tope ``cap`` (la pendiente no cuenta en el tope: el llamador pasa a la
    secante)."""
    if E * (kappa - kappa_0) <= tol:
        return 0.0, 0.0
    d, dd, _ = evaluate_damage_law(kappa, kappa_0, law, alpha, w_c, h, cap)
    return d, dd


@njit(cache=True)
def crack_band_damage(kappa, E, kappa_0, tol, cap, law, alpha, w_c, h):
    """Daño ``d(κ)`` de :func:`crack_band_damage_slope`."""
    d, _ = crack_band_damage_slope(kappa, E, kappa_0, tol, cap, law, alpha, w_c, h)
    return d


@njit(cache=True)
def implex_update(eps_eq, kappa_n, kappa_prev, d_n, err_n, ratio, E, kappa_0, tol, cap,
                  law, alpha, w_c, h):
    """Paso IMPL-EX de un punto (spec §8). Devuelve
    ``(d̃, κ_{n+1}, κ_prev', e)``: el daño con el que se calcula el
    esfuerzo, la variable histórica implícita, la del paso anterior para la
    extrapolación siguiente y el indicador de error ``|d(κ_{n+1}) − d̃|``.

    ``ratio < 0`` (fuera de un paso): daño congelado en ``d̃_n`` y estado
    sin cambio, así que el post-proceso y los análisis lineales sobre el
    estado de una etapa reproducen el esfuerzo del paso confirmado.
    """
    if ratio < 0.0:
        return d_n, kappa_n, kappa_prev, err_n
    kappa_tilde = kappa_n + ratio * (kappa_n - kappa_prev)
    d_tilde = crack_band_damage(kappa_tilde, E, kappa_0, tol, cap, law, alpha, w_c, h)
    if d_tilde < d_n:
        d_tilde = d_n
    kappa_new = eps_eq if eps_eq > kappa_n else kappa_n
    d_impl = crack_band_damage(kappa_new, E, kappa_0, tol, cap, law, alpha, w_c, h)
    return d_tilde, kappa_new, kappa_n, abs(d_impl - d_tilde)


def parse_integration(owner: str, integration) -> int:
    """Código de la integración (``'implicit'`` o ``'impl-ex'``)."""
    code = INTEGRATION_CODES.get(integration) if isinstance(integration, str) else None
    if code is None:
        raise ValueError(
            f"{owner}: integration={integration!r} no válida; use 'implicit' (por omisión) "
            f"o 'impl-ex' (Oliver, Huespe y Cante 2008; requiere un solver que declare su "
            f"pseudo-tiempo, ADR 0034).")
    return code


def validate_crack_band_parameters(owner: str, E, nu, f_t, G_f, damage_max, density) -> None:
    """Validación temprana de los parámetros de la variante."""
    for name, value in (("E", E), ("f_t", f_t), ("G_f", G_f)):
        try:
            ok = math.isfinite(float(value)) and float(value) > 0.0
        except (TypeError, ValueError):
            ok = False
        if not ok:
            raise ValueError(f"{owner}: {name} debe ser un número finito > 0 (recibido {value!r}).")
    if not (-1.0 < float(nu) < 0.5):
        raise ValueError(f"{owner}: nu debe estar en (−1, 0.5) (recibido {nu!r}).")
    if not (0.0 < float(damage_max) < 1.0):
        raise ValueError(f"{owner}: damage_max debe estar en (0, 1) (recibido {damage_max!r}).")
    if density is not None and density < 0.0:
        raise ValueError(f"{owner}: density={density} no puede ser negativa.")


class CrackBandMixin:
    """Lo común a las variantes 2D y 3D: ``h`` por punto, ``α(h)`` y la ley
    de ablandamiento.

    La subclase define ``E``, ``f_t``, ``G_f``, ``kappa_0`` y
    ``damage_max`` y llama a :meth:`_set_softening` antes de usar estos
    métodos.
    """

    def _set_softening(self, owner: str, softening) -> None:
        self._softening_code = parse_softening(owner, softening)
        self.w_c = (float(hordijk_critical_opening(self.f_t, self.G_f))
                    if self._softening_code == HORDIJK else 0.0)

    @property
    def softening(self) -> str:
        return "hordijk" if self._softening_code == HORDIJK else "exponential"

    @property
    def uses_characteristic_length(self) -> bool:
        return True

    @property
    def characteristic_length(self) -> float:
        """``ℓ_ch = E·G_f/f_t²`` (Hillerborg, Modéer y Petersson 1976)."""
        return self.E * self.G_f / self.f_t ** 2

    @property
    def characteristic_length_limit(self) -> float:
        """``h_max``: ``2·ℓ_ch`` con la exponencial (α > 0) y
        ``E·w_c/(f_t·|g'(0)|) = 0.7382·ℓ_ch`` con Hordijk (sin retroceso en el
        pico, spec §2b)."""
        if self._softening_code == HORDIJK:
            return HORDIJK_H_MAX_FACTOR * self.characteristic_length
        return 2.0 * self.characteristic_length

    def alpha_for(self, h: float) -> float:
        """``α(h) = 1/(G_f/(h·f_t) − κ₀/2)``."""
        return float(crack_band_alpha(float(h), self.E, self.f_t, self.G_f))

    # --- integración (spec §7 y §8; ADR 0034) ---------------------------

    @property
    def integration(self) -> str:
        return "impl-ex" if self._integration_code == IMPL_EX else "implicit"

    @property
    def uses_pseudo_time(self) -> bool:
        """IMPL-EX necesita la razón de incrementos de pseudo-tiempo."""
        return self._integration_code == IMPL_EX

    def set_pseudo_time_ratio(self, ratio) -> None:
        self._ratio = NO_STEP_RATIO if ratio is None else float(ratio)

    def _state_row(self, state_vars) -> np.ndarray:
        """Estado ``[κ, d, h, κ_prev, e]`` del punto como fila (orden de
        ``STATE_SCHEMA``), con los valores iniciales para lo que falte."""
        k0 = float(self.kappa_0)
        if state_vars is None:
            return np.array([k0, 0.0, 0.0, k0, 0.0])
        kappa = float(state_vars.get('kappa', k0))
        return np.array([kappa, float(state_vars.get('damage', 0.0)),
                         float(state_vars.get('h', 0.0)),
                         float(state_vars.get('kappa_prev', kappa)),
                         float(state_vars.get('damage_error', 0.0))])

    def _point_state(self, h: float) -> dict:
        k0 = float(self.kappa_0)
        return {'kappa': k0, 'damage': 0.0, 'h': float(h), 'kappa_prev': k0,
                'damage_error': 0.0}

    def initial_state(self) -> dict:
        """Estado sin longitud característica (``h = 0``): vale en régimen
        elástico; al empezar el daño sin ``h`` el material avisa."""
        return self._point_state(0.0)

    def initial_point_state(self, characteristic_length: float) -> dict:
        h = float(characteristic_length)
        name = type(self).__name__
        if not (math.isfinite(h) and h > 0.0):
            raise ValueError(f"{name}: longitud característica h = {h} no válida (> 0).")
        h_max = self.characteristic_length_limit
        if h >= h_max:
            params = (f"E = {self.E:.4g}, G_f = {self.G_f:.4g}, f_t = {self.f_t:.4g}, "
                      f"ℓ_ch = E·G_f/f_t² = {self.characteristic_length:.4g}")
            if self._softening_code == HORDIJK:
                raise ValueError(
                    f"{name}: la longitud característica h = {h:.4g} supera el límite de la "
                    f"ley de Hordijk h_max = {HORDIJK_H_MAX_FACTOR:.4f}·ℓ_ch = {h_max:.4g} "
                    f"({params}): la respuesta tendría retroceso en el pico (el daño saltaría "
                    f"y la energía disipada superaría G_f/h). Use elementos más pequeños.")
            raise ValueError(
                f"{name}: la longitud característica h = {h:.4g} supera el límite "
                f"h_max = 2·E·G_f/f_t² = {h_max:.4g} ({params}): la rama de ablandamiento "
                f"tendría retroceso en el propio material (la energía elástica en el pico ya "
                f"supera G_f/h). Use elementos más pequeños.")
        return self._point_state(h)

    def _missing_h(self) -> ValueError:
        return ValueError(
            f"{type(self).__name__}: el daño empezó en un punto sin longitud "
            f"característica h. La siembran los elementos sólidos (Tri3, Quad4, Tri6, "
            f"Quad8, Quad9, Hex8, Tet4, Hex20, Hex27, Tet10) al construirse; con otro "
            f"elemento este material no se admite.")

    def batch_report(self, n_flagged: int) -> None:
        """Puntos marcados por el kernel por lotes: daño sin ``h``."""
        raise self._missing_h()


__all__ = [
    "DEFAULT_INTEGRATION",
    "DEFAULT_SOFTENING",
    "EXPONENTIAL",
    "HORDIJK",
    "crack_band_damage_slope",
    "hordijk_critical_opening",
    "IMPLICIT",
    "IMPL_EX",
    "INTEGRATION_CODES",
    "NO_STEP_RATIO",
    "crack_band_damage",
    "implex_update",
    "parse_integration",
    "CRACK_BAND_DAMAGE_MAX",
    "CRACK_BAND_EQUIVALENT_STRAINS",
    "DEFAULT_CRACK_BAND_EQUIVALENT_STRAIN",
    "CrackBandMixin",
    "crack_band_alpha",
    "validate_crack_band_parameters",
]
