"""``CohesiveDamageMixedMode`` — ley cohesiva de modo mixto I–II con apertura
efectiva y daño escalar (material de interfaz de superficie, ADR 0033).

Apertura efectiva de Camacho y Ortiz (1996) con peso ``β`` sobre el
deslizamiento, energía de fractura independiente del modo (``G_II = G_I``,
resistencia a cortante ``β·σ_t0``) y ablandamiento lineal, exponencial o
bilineal (fib Model Code 2010). Con ``β = 0`` es ``CohesiveDamageIsotropic``.
Spec ``docs/specs/CohesiveDamageMixedMode.md``.
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np

from solidum.core.interface_material import InterfaceMaterial
from solidum.materials.interface._softening import SofteningCurve
from solidum.registry import InterfaceMaterialRegistry


@InterfaceMaterialRegistry.register
class CohesiveDamageMixedMode(InterfaceMaterial):
    """
    Daño cohesivo escalar de modo mixto con apertura efectiva.

    Opera sobre el salto en el sistema local de la superficie,
    ``(δ_n, δ_s)`` en 2D o ``(δ_n, δ_s1, δ_s2)`` en 3D, y devuelve la tracción
    en los mismos ejes::

        δ_eff = √(⟨δ_n⟩² + β²·‖δ_s‖²),        κ = máx(κ, δ_eff)
        t_n = S(κ)·δ_n  (δ_n ≥ 0),   t_n = K_e·δ_n  (δ_n < 0, sin daño)
        t_s = β²·S(κ)·δ_s,                    S(κ) = T(κ)/κ

    con ``T(κ)`` la envolvente de ablandamiento (rama elástica ``K_e·κ`` hasta
    ``κ_0 = σ_t0/K_e``). En carga activa la tracción efectiva
    ``√(t_n² + ‖t_s‖²/β²)`` sigue ``T(δ_eff)``; la iniciación es el criterio
    cuadrático ``(t_n/σ_t0)² + (‖t_s‖/(β·σ_t0))² = 1`` y toda trayectoria
    proporcional disipa ``G_F``. Sin fricción: con la grieta comprimida el
    deslizamiento también daña.

    Parameters
    ----------
    sigma_t0 : float
        Resistencia a tracción ``σ_t0`` [Pa].
    G_f : float
        Energía de fractura ``G_F`` [N/m], la misma en cualquier modo.
    K_e : float
        Rigidez de penalización normal [Pa/m]; la tangencial es ``β²·K_e``.
        Sin valor por omisión (guía: ``K_e ≈ 10·E/ℓ_c``).
    beta : float
        Peso del deslizamiento (``≥ 0``); resistencia a cortante
        ``τ_c = β·σ_t0``. Sin valor por omisión.
    softening : str
        ``'linear'``, ``'exponential'`` o ``'bilinear'``.
    knee_stress_ratio, final_opening_factor : float, optional
        Sólo con ``'bilinear'``: ``r`` (esfuerzo en el quiebre entre ``σ_t0``)
        y ``c`` (apertura final por ``σ_t0/H``). Por omisión los del Model
        Code 2010, ``0.2`` y ``5.0``.

    Notes
    -----
    Tangente consistente, simétrica en todas las ramas (potencial secante),
    con ``u = (⟨δ_n⟩, β²·δ_s)/δ_eff``::

        carga activa:          D = S·(W − u⊗u) + T'(κ)·u⊗u
        elástica / descarga:   D = S·W,           W = diag(1, β²[, β²])

    y ``D_nn = K_e``, ``D_ns = 0`` con ``δ_n < 0``. Escrita así, en modo I
    puro ``D_nn = T'`` sin restar números casi iguales.
    """
    GEOMETRY = "surface"
    JUMP_DIM = None
    PRIMARY_STATE_VAR = "damage"
    IS_SYMMETRIC = True

    def __init__(self, sigma_t0: float, G_f: float, K_e: float, beta: float, softening: str,
                 knee_stress_ratio: Optional[float] = None,
                 final_opening_factor: Optional[float] = None):
        name = type(self).__name__
        for label, value in (("sigma_t0", sigma_t0), ("G_f", G_f), ("K_e", K_e)):
            if not (math.isfinite(value) and value > 0.0):
                raise ValueError(f"{name}: {label} debe ser > 0 y finito (recibido {value}).")
        if not (math.isfinite(beta) and beta >= 0.0):
            raise ValueError(f"{name}: beta debe ser ≥ 0 y finito (recibido {beta}).")
        self.sigma_t0 = float(sigma_t0)
        self.G_f = float(G_f)
        self.K_e = float(K_e)
        self.beta = float(beta)
        self.softening = softening
        self._curve = SofteningCurve(self.sigma_t0, self.G_f, self.K_e, softening,
                                     knee_stress_ratio=knee_stress_ratio,
                                     final_opening_factor=final_opening_factor, owner=name)
        self.knee_stress_ratio = self._curve.r
        self.final_opening_factor = self._curve.c if softening == "bilinear" else None
        self.kappa_0 = self._curve.kappa_0

    @property
    def shear_strength(self) -> float:
        """Resistencia a cortante puro ``τ_c = β·σ_t0`` [Pa]."""
        return self.beta * self.sigma_t0

    def envelope(self, kappa: float):
        """``(T(κ), dT/dκ)``: tracción efectiva sobre la envolvente."""
        return self._curve(kappa)

    def compute_traction(self, jump: np.ndarray, state_vars=None):
        jump = np.asarray(jump, dtype=float)
        n = jump.size
        if jump.ndim != 1 or n not in (2, 3):
            raise ValueError(
                f"{type(self).__name__}: el salto debe tener 2 componentes (n, s) o 3 "
                f"(n, s1, s2); recibido forma {jump.shape}."
            )
        kappa_old = self.kappa_0 if state_vars is None else state_vars.get("kappa", self.kappa_0)
        b2 = self.beta * self.beta
        d_n = float(jump[0])
        d_s = jump[1:]
        d_np = max(d_n, 0.0)
        d_eff = math.sqrt(d_np * d_np + b2 * float(d_s @ d_s))

        loading = d_eff > kappa_old
        kappa = d_eff if loading else kappa_old
        T, dT = self._curve(kappa)
        S = T / kappa

        weights = np.full(n, b2)
        weights[0] = 1.0
        traction = S * weights * jump
        if loading:
            u = np.empty(n)
            u[0] = d_np / d_eff
            u[1:] = b2 * d_s / d_eff
            uu = np.outer(u, u)
            tangent = S * (np.diag(weights) - uu) + dT * uu
        else:
            tangent = S * np.diag(weights)
        if d_n < 0.0:
            # Cierre: rigidez inicial sin daño en la normal (como la ley de
            # modo I; Alfaiate, Wells y Sluys 2002, p. 667). u_n = 0, así que
            # el bloque tangencial no se acopla con la normal.
            traction[0] = self.K_e * d_n
            tangent[0, :] = 0.0
            tangent[:, 0] = 0.0
            tangent[0, 0] = self.K_e

        omega = 1.0 - S / self.K_e
        return traction, tangent, {"kappa": kappa, "damage": omega}
