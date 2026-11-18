"""Lo común a ``TensionCompressionDamage2D/3D``: estado por punto y daño de tracción.

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

Bloque 1 de la implementación: sólo evoluciona ``d⁺``. ``d⁻`` se arrastra del
estado sin cambio (0 en un análisis); su ley (aplastamiento con Drucker-Prager
y la parábola de Feenstra, spec §3) es el bloque 2.
"""
from __future__ import annotations

import numpy as np
from numba import njit

from solidum.materials._crack_band import (
    IMPL_EX,
    CrackBandMixin,
    crack_band_damage,
    implex_update,
)

# Índices de la fila de estado.
KT, KC, DT, DC, H, KT_PREV, KC_PREV, ERR = range(8)
STATE_SCHEMA_TC = {'kappa_t': (), 'kappa_c': (), 'damage_t': (), 'damage_c': (), 'h': (),
                   'kappa_t_prev': (), 'kappa_c_prev': (), 'damage_error': ()}


@njit(cache=True)
def tension_damage_update(eps_eq, kappa_n, kappa_prev_n, d_n, err_n, E, kappa_0, alpha,
                          tol, cap, integration, ratio):
    """Daño de tracción ``d⁺`` de un punto (spec §2 y §8).

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
            eps_eq, kappa_n, kappa_prev_n, d_n, err_n, ratio, E, kappa_0, alpha, tol, cap)
        return d, kappa_new, kappa_prev, err, 0.0
    loading = eps_eq > kappa_n
    kappa_new = eps_eq if loading else kappa_n
    d = crack_band_damage(kappa_new, E, kappa_0, alpha, tol, cap)
    dd = 0.0
    if loading and d > 0.0 and d < cap and eps_eq > 0.0:
        dd = (1.0 - d) * (1.0 / kappa_new + alpha)
    return d, kappa_new, kappa_n, 0.0, dd


class TensionCompressionMixin(CrackBandMixin):
    """Estado por punto del modelo 2 sobre las piezas de la banda de grieta.

    La subclase define ``E``, ``f_t``, ``G_f``, ``kappa_0``, ``damage_max``,
    ``_integration_code`` y ``_ratio`` (como el modelo 1), y declara
    ``STATE_SCHEMA = STATE_SCHEMA_TC``.
    """

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
        return np.array([k0, 0.0, 0.0, 0.0, float(h), k0, 0.0, 0.0])

    def _point_state(self, h: float) -> dict:
        return self._row_to_state(self._initial_row(h))

    @staticmethod
    def _row_to_state(row) -> dict:
        return {name: float(row[i]) for i, name in enumerate(STATE_SCHEMA_TC)}


__all__ = [
    "DC",
    "DT",
    "ERR",
    "H",
    "KC",
    "KC_PREV",
    "KT",
    "KT_PREV",
    "STATE_SCHEMA_TC",
    "TensionCompressionMixin",
    "tension_damage_update",
]
