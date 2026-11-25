"""Lo común a ``TensionCompressionDamage2D/3D``: estado por punto y los dos daños.

Modelo 2 de la familia de concreto (spec ``docs/specs/TensionCompressionDamage2D.md``):
dos daños escalares, ``d⁺`` de tracción y ``d⁻`` de compresión, con la ley de
esfuerzo de Faria, Oliver y Cervera (1998),

    σ = (1 − d⁺)·σ̄⁺ + (1 − d⁻)·σ̄⁻ = (1 − d⁻)·σ̄ − (d⁺ − d⁻)·σ̄⁺,

``σ̄ = C_e·ε`` el esfuerzo efectivo y ``σ̄⁺`` su parte positiva
(``_spectral_split``). Una grieta cerrada recupera la rigidez en compresión.

El daño de tracción es el del modelo 1 (``_crack_band``): Rankine sobre el
esfuerzo efectivo, ley exponencial con ``α(h)`` de la banda de grieta, tope
``damage_max`` e integración implícita o IMPL-EX (ADR 0034). Aquí vive su
actualización en un punto, la misma en 2D y en 3D (Reglas.md §1: dos casos
reales); cada dimensión conserva su núcleo compilado.

Estado por punto (orden de ``STATE_SCHEMA``)::

    [κ⁺, κ⁻, d⁺, d⁻, h, κ⁺_prev, κ⁻_prev, e]

``κ±`` las variables históricas, ``d±`` los daños con que se calculó el
esfuerzo, ``h`` la longitud característica, ``κ±_prev`` las del paso
confirmado anterior (extrapolación de IMPL-EX) y ``e`` el indicador de error
de IMPL-EX (0 en la implícita).

El daño de compresión (spec §3) es opcional: con ``f_c`` y ``G_c`` se activa
(``f_c > 0`` en el núcleo); sin ellos ``d⁻`` se arrastra del estado sin cambio
y la compresión es elástica. Medida de Drucker-Prager sobre la parte negativa
del esfuerzo efectivo, normalizada a ``|σ̄|/E`` en compresión uniaxial,

    ε⁻_eq = 3·(K·σ̄⁻_oct + τ̄⁻_oct) / ((√2 − K)·E),   K = √2·(β − 1)/(2β − 1),

con ``σ̄⁻_oct = tr σ̄⁻/3``, ``τ̄⁻_oct = √(⅔·J₂(σ̄⁻))`` del tensor 3D, y ley
parabólica de Feenstra (1993) en ``κ⁻ = max ε⁻_eq`` con la energía de
aplastamiento ``G_c/h`` en la rama descendente, ``d⁻ = 1 − σ_c(κ⁻)/(E·κ⁻)``.
``κ⁻`` empieza en ``f_c/(3E)`` (fin del tramo elástico), como ``κ⁺`` en
``κ_0``: así IMPL-EX nunca crea daño en un punto que la integración implícita
deja elástico (con ``κ⁻ = 0`` la extrapolación del crecimiento elástico
podría dar ``d̃⁻ > 0``, irreversible, a un punto que se detiene justo bajo el
umbral).
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

from solidum.constants import ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL
from solidum.materials._crack_band import (
    IMPL_EX,
    CrackBandMixin,
    crack_band_damage_slope,
    implex_update,
)

# Índices de la fila de estado.
KT, KC, DT, DC, H, KT_PREV, KC_PREV, ERR = range(8)
STATE_SCHEMA_TC = {'kappa_t': (), 'kappa_c': (), 'damage_t': (), 'damage_c': (), 'h': (),
                   'kappa_t_prev': (), 'kappa_c_prev': (), 'damage_error': ()}
# β = f_bc/f_c por omisión (Kupfer, spec §3).
DEFAULT_BETA = 1.16
SQRT2 = math.sqrt(2.0)


def drucker_prager_k(beta: float) -> float:
    """``K = √2·(β − 1)/(2β − 1)`` de la medida de compresión (spec §3)."""
    return SQRT2 * (beta - 1.0) / (2.0 * beta - 1.0)


@njit(cache=True)
def tension_damage_update(eps_eq, kappa_n, kappa_prev_n, d_n, err_n, E, kappa_0, alpha,
                          tol, cap, integration, ratio, law, w_c, h):
    """Daño de tracción ``d⁺`` de un punto (spec §2 y §8), con la ley de
    ablandamiento ``law`` del modelo 1 (exponencial con ``alpha`` o Hordijk
    con ``w_c`` y ``h``; ``_softening``).

    ``eps_eq = ⟨σ̄_I⟩₊/E`` es la medida de Rankine. Devuelve
    ``(d⁺, κ⁺, κ⁺_prev, e, dd/dκ)``, con ``dd/dκ > 0`` sólo en carga activa
    de la integración implícita con el daño entre 0 y el tope (la rama
    donde la tangente lleva el término del daño); 0 en todo otro caso.

    - Implícita: ``κ⁺ = max(κ⁺_n, ε_eq)`` y ``d⁺ = d(κ⁺)`` (Kuhn-Tucker con
      la tolerancia de admisibilidad del ADR 0006, como el modelo 1).
    - IMPL-EX: ``d̃⁺ = max(d̃⁺_n, d(κ̃⁺))`` con ``κ̃⁺`` extrapolado de los dos
      últimos pasos y ``κ⁺`` actualizada implícitamente
      (``_crack_band.implex_update``); fuera de un paso (``ratio < 0``),
      daño y estado congelados.
    """
    if integration == IMPL_EX:
        d, kappa_new, kappa_prev, err = implex_update(
            eps_eq, kappa_n, kappa_prev_n, d_n, err_n, ratio, E, kappa_0, tol, cap,
            law, alpha, w_c, h)
        return d, kappa_new, kappa_prev, err, 0.0
    loading = eps_eq > kappa_n
    kappa_new = eps_eq if loading else kappa_n
    d, dd_law = crack_band_damage_slope(kappa_new, E, kappa_0, tol, cap, law, alpha, w_c, h)
    dd = 0.0
    if loading and d > 0.0 and d < cap and eps_eq > 0.0:
        dd = dd_law
    return d, kappa_new, kappa_n, 0.0, dd


# ---------------------------------------------------------------------------
# Daño de compresión (spec §3)
# ---------------------------------------------------------------------------

@njit(cache=True)
def compression_ultimate(h, E, f_c, G_c):
    """``κ_u = κ_c + 3·G_c/(2·h·f_c)``, fin de la rama descendente.

    El área de la rama descendente, ``⅔·f_c·(κ_u − κ_c)``, es ``G_c/h``. Sin
    ``h`` sembrada (``h ≤ 0``) la rama cae de golpe (``κ_u = κ_c``): el daño
    salta al tope y el llamador avisa, como en tracción.
    """
    kappa_c = 5.0 * f_c / (3.0 * E)
    if h <= 0.0:
        return kappa_c
    return kappa_c + 1.5 * G_c / (h * f_c)


@njit(cache=True)
def compression_damage(kappa, E, f_c, kappa_u, tol, cap):
    """``(d⁻, dd⁻/dκ)`` de la parábola de Feenstra (spec §3).

    ``σ_c(κ)``: elástica hasta ``κ_c/3 = f_c/(3E)``; rama ascendente
    ``(f_c/3)·(1 + 4ξ − 2ξ²)`` hasta ``κ_c = 5f_c/(3E)`` (pendiente ``E`` al
    empezar y nula en el pico ``f_c``); rama descendente ``f_c·(1 − η²)``
    hasta ``κ_u``; nula después. ``d⁻ = 1 − σ_c/(E·κ)``, con el umbral de
    admisibilidad ``E·(κ − κ_c/3) ≤ tol`` ⇒ 0 y el tope ``cap`` (en el tope,
    ``dd/dκ = 0``).
    """
    kc3 = f_c / (3.0 * E)
    if E * (kappa - kc3) <= tol:
        return 0.0, 0.0
    kc = 5.0 * kc3
    if kappa <= kc:
        w = kc - kc3
        xi = (kappa - kc3) / w
        s = f_c / 3.0 * (1.0 + 4.0 * xi - 2.0 * xi * xi)
        ds = f_c / 3.0 * (4.0 - 4.0 * xi) / w
    elif kappa < kappa_u:
        w = kappa_u - kc
        eta = (kappa - kc) / w
        s = f_c * (1.0 - eta * eta)
        ds = -2.0 * f_c * eta / w
    else:
        return cap, 0.0
    d = 1.0 - s / (E * kappa)
    if d >= cap:
        return cap, 0.0
    if d < 0.0:
        # Redondeo junto al umbral: σ_c ≤ E·κ en toda la rama ascendente.
        d = 0.0
    return d, (s - kappa * ds) / (E * kappa * kappa)


@njit(cache=True)
def compression_damage_update(eps_eq, kappa_n, kappa_prev_n, d_n, err_n, E, f_c, kappa_u,
                              tol, cap, integration, ratio):
    """Daño de compresión ``d⁻`` de un punto (spec §3 y §8), con la misma
    estructura que :func:`tension_damage_update`. Devuelve
    ``(d⁻, κ⁻, κ⁻_prev, e, dd/dκ)``.

    IMPL-EX: ``κ̃⁻ = κ⁻_n + r·(κ⁻_n − κ⁻_{n−1})``, ``d̃⁻ = max(d̃⁻_n, d⁻(κ̃⁻))``
    y ``κ⁻`` actualizada implícitamente; fuera de un paso (``r < 0``), daño
    y estado congelados.
    """
    if integration == IMPL_EX:
        if ratio < 0.0:
            return d_n, kappa_n, kappa_prev_n, err_n, 0.0
        kappa_tilde = kappa_n + ratio * (kappa_n - kappa_prev_n)
        d_tilde, _ = compression_damage(kappa_tilde, E, f_c, kappa_u, tol, cap)
        if d_tilde < d_n:
            d_tilde = d_n
        kappa_new = eps_eq if eps_eq > kappa_n else kappa_n
        d_impl, _ = compression_damage(kappa_new, E, f_c, kappa_u, tol, cap)
        return d_tilde, kappa_new, kappa_n, abs(d_impl - d_tilde), 0.0
    loading = eps_eq > kappa_n
    kappa_new = eps_eq if loading else kappa_n
    d, dd = compression_damage(kappa_new, E, f_c, kappa_u, tol, cap)
    if not (loading and d > 0.0 and d < cap):
        dd = 0.0
    return d, kappa_new, kappa_n, 0.0, dd


@njit(cache=True)
def compression_measure(sm_xx, sm_yy, sm_zz, sm_xy, sm_yz, sm_xz, E, K):
    """``ε⁻_eq`` de Drucker-Prager sobre ``σ̄⁻`` (tensor 3D) y su gradiente.

    Devuelve ``(ε⁻_eq, g)``, ``g = ∂ε⁻_eq/∂σ̄⁻`` como tensor
    ``[g_xx, g_yy, g_zz, g_xy, g_yz, g_xz]`` (sin el factor 2 de Voigt):
    ``g = (K·I + s/τ̄_oct)/((√2 − K)·E)`` con ``s`` el desviador. Con
    ``τ̄_oct = 0`` (estado hidrostático) la parte desviadora se omite; ahí no
    hay carga activa (``ε⁻_eq ≤ 0`` con ``K ≥ 0``), así que el gradiente no se
    usa.
    """
    mean = (sm_xx + sm_yy + sm_zz) / 3.0
    sxx = sm_xx - mean
    syy = sm_yy - mean
    szz = sm_zz - mean
    J2 = 0.5 * (sxx * sxx + syy * syy + szz * szz) + sm_xy * sm_xy + sm_yz * sm_yz + sm_xz * sm_xz
    tau = math.sqrt(2.0 * J2 / 3.0)
    den = (SQRT2 - K) * E
    eps = 3.0 * (K * mean + tau) / den
    g = np.empty(6)
    if tau > 0.0:
        g[0] = (K + sxx / tau) / den
        g[1] = (K + syy / tau) / den
        g[2] = (K + szz / tau) / den
        g[3] = sm_xy / tau / den
        g[4] = sm_yz / tau / den
        g[5] = sm_xz / tau / den
    else:
        g[0] = K / den
        g[1] = K / den
        g[2] = K / den
        g[3] = 0.0
        g[4] = 0.0
        g[5] = 0.0
    return eps, g


def parse_compression(owner: str, f_c, G_c, beta, f_t: float, E: float):
    """Valida ``f_c``, ``G_c`` y ``β`` del aplastamiento (spec del modelo 2,
    §3 y ``validity``) y devuelve ``(β, f_c, G_c, κ⁻_0)``.

    Sin ``f_c`` ni ``G_c`` la compresión es elástica: ``(β, None, None, 0)``.
    ``κ⁻_0 = f_c/(3E)`` es el fin del tramo elástico de la parábola de
    Feenstra. Lo usan el modelo 2 y el concreto entre grietas del modelo 3.
    """
    try:
        beta = float(beta)
    except (TypeError, ValueError):
        beta = math.nan
    if not (math.isfinite(beta) and beta >= 1.0):
        raise ValueError(
            f"{owner}: beta = f_bc/f_c debe ser un número finito ≥ 1 (recibido {beta!r}); "
            f"1.16 por omisión (Kupfer). Con β < 1 la compresión hidrostática dañaría "
            f"el material, contra la superficie de Drucker-Prager de la spec §3.")
    if f_c is None and G_c is None:
        return beta, None, None, 0.0
    if f_c is None or G_c is None:
        raise ValueError(
            f"{owner}: f_c y G_c van juntos (aplastamiento en compresión, spec §3); "
            f"recibido f_c = {f_c!r}, G_c = {G_c!r}. Sin los dos, la compresión es "
            f"elástica.")
    for name, value in (("f_c", f_c), ("G_c", G_c)):
        try:
            ok = math.isfinite(float(value)) and float(value) > 0.0
        except (TypeError, ValueError):
            ok = False
        if not ok:
            raise ValueError(
                f"{owner}: {name} debe ser un número finito > 0 (recibido {value!r}); "
                f"f_c es la resistencia a compresión uniaxial en valor absoluto.")
    if float(f_c) <= f_t:
        raise ValueError(
            f"{owner}: f_c = {f_c} debe ser mayor que f_t = {f_t} "
            f"(spec, validity).")
    return beta, float(f_c), float(G_c), float(f_c) / (3.0 * E)


def compression_kernel_params(f_c, G_c, beta: float):
    """``(f_c, G_c, K, tol_c)`` para los núcleos; ``f_c = 0`` si no hay
    aplastamiento. ``tol_c`` es la tolerancia de admisibilidad del ADR 0006
    con la escala ``f_c/3`` (fin del tramo elástico), con la misma expresión
    que los adaptadores por lotes para que coincidan bit a bit."""
    if f_c is None:
        return 0.0, 0.0, 0.0, 0.0
    return (f_c, G_c, drucker_prager_k(beta),
            ADMISSIBILITY_TOL_ABS + ADMISSIBILITY_TOL_REL * (f_c / 3.0))


class TensionCompressionMixin(CrackBandMixin):
    """Estado por punto del modelo 2 sobre las piezas de la banda de grieta.

    La subclase define ``E``, ``f_t``, ``G_f``, ``kappa_0``, ``damage_max``,
    ``_integration_code`` y ``_ratio`` (como el modelo 1), los parámetros de
    compresión (``_set_compression``) y declara ``STATE_SCHEMA = STATE_SCHEMA_TC``.
    """

    def _set_compression(self, f_c, G_c, beta) -> None:
        """Valida y guarda ``f_c``, ``G_c`` y ``β`` (spec §3, ``validity``).

        Sin ``f_c`` ni ``G_c`` la compresión es elástica (``d⁻ ≡ 0``).
        """
        self.beta, self.f_c, self.G_c, self._kappa_c0 = parse_compression(
            type(self).__name__, f_c, G_c, beta, self.f_t, self.E)

    @property
    def has_compression_damage(self) -> bool:
        """``True`` con ``f_c`` y ``G_c`` (aplastamiento activo)."""
        return self.f_c is not None

    def _compression_params(self):
        """``(f_c, G_c, K, tol_c)`` para el núcleo (:func:`compression_kernel_params`)."""
        return compression_kernel_params(self.f_c, self.G_c, self.beta)

    def admissibility_scale(self, state_vars=None) -> float:
        """Esfuerzo umbral de tracción ``E·κ_0 = f_t`` (ADR 0006)."""
        return self.E * self.kappa_0

    def _state_row(self, state_vars) -> np.ndarray:
        """Estado del punto como fila (orden de ``STATE_SCHEMA``), con los
        valores iniciales para lo que falte."""
        row = self._initial_row(0.0)
        if state_vars is None:
            return row
        for i, name in enumerate(STATE_SCHEMA_TC):
            if name in state_vars:
                row[i] = float(state_vars[name])
        if 'kappa_t_prev' not in state_vars:
            row[KT_PREV] = row[KT]
        if 'kappa_c_prev' not in state_vars:
            row[KC_PREV] = row[KC]
        return row

    def _initial_row(self, h: float) -> np.ndarray:
        k0 = float(self.kappa_0)
        kc0 = float(self._kappa_c0)
        return np.array([k0, kc0, 0.0, 0.0, float(h), k0, kc0, 0.0])

    def _point_state(self, h: float) -> dict:
        return self._row_to_state(self._initial_row(h))

    @staticmethod
    def _row_to_state(row) -> dict:
        return {name: float(row[i]) for i, name in enumerate(STATE_SCHEMA_TC)}


__all__ = [
    "DC",
    "DEFAULT_BETA",
    "DT",
    "ERR",
    "H",
    "KC",
    "KC_PREV",
    "KT",
    "KT_PREV",
    "STATE_SCHEMA_TC",
    "TensionCompressionMixin",
    "compression_damage",
    "compression_damage_update",
    "compression_kernel_params",
    "compression_measure",
    "compression_ultimate",
    "drucker_prager_k",
    "parse_compression",
    "tension_damage_update",
]
