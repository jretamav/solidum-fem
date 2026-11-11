"""Envolvente de ablandamiento ``T(κ)`` de las leyes cohesivas de daño.

La comparten ``CohesiveDamageIsotropic`` (modo I) y
``CohesiveDamageMixedMode`` (modo mixto). Para ``κ ≤ κ_0 = σ_t0/K_e`` es la
rama elástica ``K_e·κ``; para ``κ > κ_0`` el ablandamiento se mide con la
apertura ``w = κ − κ_0`` y encierra la energía ``H = G_F − σ_t0·κ_0/2``, de
modo que el área total bajo la curva es exactamente ``G_F`` para todas las
formas (spec ``CohesiveDamageMixedMode`` §5):

- ``linear``: ``σ_t0·(1 − w/w_c')`` con ``w_c' = 2H/σ_t0``; termina en
  ``κ = 2G_F/σ_t0``;
- ``exponential``: ``σ_t0·exp(−σ_t0·w/H)``;
- ``bilinear``: quiebre en ``(w_1, r·σ_t0)`` y fin en ``w_c' = c·H/σ_t0``, con
  ``w_1 = (2 − r·c)·H/σ_t0`` (fib Model Code 2010, §5.1.5.2: ``r = 0.2``,
  ``c = 5``; Petersson 1981: ``r = 1/3``, ``c = 3.6``).

No confundir con ``solidum/materials/_softening.py``: aquél da el daño
``ω(κ)`` de los materiales continuos de daño escalar (``IsotropicDamage*``),
con tope ``DAMAGE_MAX``; aquí la envolvente se evalúa como tracción, sin tope,
porque con ``K_e`` de penalización ``1 − ω ~ κ_0/κ`` y truncar ``ω`` cambiaría
la energía disipada (spec ``CohesiveDamageIsotropic`` §12).
"""
from __future__ import annotations

import math
from typing import Optional, Tuple

SHAPES = ("linear", "exponential", "bilinear")

# fib Model Code 2010, §5.1.5.2: esfuerzo en el quiebre 0.2·f_ctm y apertura
# final 5·G_F/f_ctm.
MC2010_KNEE_STRESS_RATIO = 0.2
MC2010_FINAL_OPENING_FACTOR = 5.0


class SofteningCurve:
    """``T(κ)`` y su pendiente para una forma de ablandamiento dada.

    ``owner`` es el nombre de la ley que la usa, para los mensajes de error.
    """

    def __init__(self, sigma_t0: float, G_f: float, K_e: float, shape: str, *,
                 knee_stress_ratio: Optional[float] = None,
                 final_opening_factor: Optional[float] = None,
                 shapes: Tuple[str, ...] = SHAPES, owner: str = "SofteningCurve"):
        if shape not in shapes:
            raise ValueError(f"{owner}: softening debe ser uno de {shapes}, recibido {shape!r}.")
        self.sigma_t0 = sigma_t0
        self.kappa_0 = sigma_t0 / K_e
        self.K_e = K_e
        self.shape = shape
        self.H = G_f - 0.5 * sigma_t0 * self.kappa_0
        if self.H <= 0.0:
            raise ValueError(
                f"{owner}: se requiere K_e > σ_t0²/(2·G_F) para que la rama elástica no "
                f"consuma toda la energía (H = G_F − σ_t0·κ_0/2 > 0); con los parámetros "
                f"recibidos H = {self.H:.3e} ≤ 0."
            )

        bilinear = shape == "bilinear"
        if not bilinear and (knee_stress_ratio is not None or final_opening_factor is not None):
            raise ValueError(
                f"{owner}: knee_stress_ratio y final_opening_factor sólo valen con "
                f"softening='bilinear' (recibido softening={shape!r})."
            )
        unit = self.H / sigma_t0
        self.w_1 = None
        self.r = None
        if shape == "linear":
            self.w_end = 2.0 * unit
        elif shape == "exponential":
            self.w_end = None
        else:
            r = MC2010_KNEE_STRESS_RATIO if knee_stress_ratio is None else float(knee_stress_ratio)
            c = MC2010_FINAL_OPENING_FACTOR if final_opening_factor is None else float(final_opening_factor)
            if not (0.0 < r < 1.0) or not (c > 0.0) or not math.isfinite(c):
                raise ValueError(
                    f"{owner}: la bilineal requiere 0 < knee_stress_ratio < 1 y "
                    f"final_opening_factor > 0 (recibidos r = {r}, c = {c})."
                )
            if not (r * c < 2.0) or not (c * (1.0 + r) > 2.0):
                raise ValueError(
                    f"{owner}: la bilineal requiere r·c < 2 (quiebre después del pico) y "
                    f"c·(1 + r) > 2 (quiebre antes del final) para encerrar G_F; recibidos "
                    f"r = {r}, c = {c}."
                )
            self.r = r
            self.c = c
            self.w_1 = (2.0 - r * c) * unit
            self.w_end = c * unit
        self.kappa_end = None if self.w_end is None else self.kappa_0 + self.w_end

    def __call__(self, kappa: float) -> Tuple[float, float]:
        """``(T(κ), dT/dκ)``. Pendiente ``K_e`` en la rama elástica y la del
        ablandamiento después (negativa, nula con la grieta totalmente
        abierta)."""
        s, k0 = self.sigma_t0, self.kappa_0
        if kappa <= k0:
            return self.K_e * kappa, self.K_e
        w = kappa - k0
        if self.shape == "exponential":
            T = s * math.exp(-s * w / self.H)
            return T, -(s / self.H) * T
        if w >= self.w_end:
            return 0.0, 0.0
        if self.shape == "linear":
            slope = -s / self.w_end
            return s + slope * w, slope
        if w <= self.w_1:
            slope = -(1.0 - self.r) * s / self.w_1
            return s + slope * w, slope
        slope = -self.r * s / (self.w_end - self.w_1)
        return self.r * s + slope * (w - self.w_1), slope
