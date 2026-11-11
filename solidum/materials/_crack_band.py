"""Regularización por banda de grieta compartida por ``IsotropicDamageCrackBand2D/3D``.

Bažant y Oh (1983): la rama de ablandamiento de cada punto de integración se
ajusta a su longitud característica ``h`` para que el material disipe
``G_f/h`` por unidad de volumen. Con la ley exponencial del daño isótropo
(``solidum.materials._softening``), en tracción uniaxial
``σ = f_t·exp(−α(κ−κ₀))`` y la energía hasta la separación completa es
``g = ½·f_t·κ₀ + f_t/α``; igualarla a ``G_f/h`` da

    α(h) = 1 / (G_f/(h·f_t) − κ₀/2),     h < h_max = 2·E·G_f/f_t².

Ver ``docs/specs/IsotropicDamageCrackBand2D.md`` §2.

Por qué un módulo propio: las variantes 2D y 3D comparten la validación de
los parámetros, ``α(h)``, el límite ``h_max`` y el estado inicial por punto
(Reglas.md §1: dos casos reales); cada una conserva su núcleo compilado.
"""
from __future__ import annotations

import math

from numba import njit

from solidum.constants import CRACK_BAND_DAMAGE_MAX
from solidum.materials._equivalent_strain import EQUIVALENT_STRAIN_CODES

# Medidas que admiten las variantes; Rankine por omisión (decisión del
# usuario, 2026-11-11, Q1 de la spec).
CRACK_BAND_EQUIVALENT_STRAINS = tuple(EQUIVALENT_STRAIN_CODES)
DEFAULT_CRACK_BAND_EQUIVALENT_STRAIN = "rankine"


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
    """Lo común a las variantes 2D y 3D: ``h`` por punto y ``α(h)``.

    La subclase define ``E``, ``f_t``, ``G_f``, ``kappa_0`` y
    ``damage_max`` antes de usar estos métodos.
    """

    @property
    def uses_characteristic_length(self) -> bool:
        return True

    @property
    def characteristic_length_limit(self) -> float:
        """``h_max = 2·E·G_f/f_t²`` (= 2·ℓ_ch de Hillerborg)."""
        return 2.0 * self.E * self.G_f / self.f_t ** 2

    def alpha_for(self, h: float) -> float:
        """``α(h) = 1/(G_f/(h·f_t) − κ₀/2)``."""
        return float(crack_band_alpha(float(h), self.E, self.f_t, self.G_f))

    def initial_state(self) -> dict:
        """Estado sin longitud característica (``h = 0``): vale en régimen
        elástico; al empezar el daño sin ``h`` el material avisa."""
        return {'kappa': float(self.kappa_0), 'damage': 0.0, 'h': 0.0}

    def initial_point_state(self, characteristic_length: float) -> dict:
        h = float(characteristic_length)
        name = type(self).__name__
        if not (math.isfinite(h) and h > 0.0):
            raise ValueError(f"{name}: longitud característica h = {h} no válida (> 0).")
        h_max = self.characteristic_length_limit
        if h >= h_max:
            raise ValueError(
                f"{name}: la longitud característica h = {h:.4g} supera el límite "
                f"h_max = 2·E·G_f/f_t² = {h_max:.4g} (E = {self.E:.4g}, G_f = {self.G_f:.4g}, "
                f"f_t = {self.f_t:.4g}): la rama de ablandamiento tendría retroceso en el "
                f"propio material (la energía elástica en el pico ya supera G_f/h). "
                f"Use elementos más pequeños.")
        return {'kappa': float(self.kappa_0), 'damage': 0.0, 'h': h}

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
    "CRACK_BAND_DAMAGE_MAX",
    "CRACK_BAND_EQUIVALENT_STRAINS",
    "DEFAULT_CRACK_BAND_EQUIVALENT_STRAIN",
    "CrackBandMixin",
    "crack_band_alpha",
    "validate_crack_band_parameters",
]
