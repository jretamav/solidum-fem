# solidum_fem/solidum/materials/concrete_confinement.py
"""Funciones auxiliares del concreto confinado y de la energía de fractura.

Calculan los parámetros de ``ManderConcrete1D`` a partir del concreto sin
confinar y del refuerzo transversal; el material sólo recibe el resultado.
Ver ``docs/specs/ManderConcrete1D.md`` §5 y §6.

Unidades: las fórmulas empíricas llevan constantes con dimensión (el módulo
``5000·√f'c`` en MPa, la energía de fractura del Model Code en N/m). Cada
función recibe el valor de esas unidades en las del usuario (``MPa``,
``N_per_m``), con las del SI por omisión (esfuerzos en Pa, energía por
área en N/m).
"""
from __future__ import annotations

import math
from typing import Sequence

from solidum.logging import get_logger

_log = get_logger("materials.concrete")

# Validez declarada de la aproximación de Chang y Mander (1994): f'l/f'c ≤ 0.3.
_CHANG_MANDER_MAX_RATIO = 0.3

# Deformación de pico del concreto sin confinar que Mander (1988) propone
# cuando no hay ensayo («generally ε_co = 0.002 can be assumed», p. 1808).
DEFAULT_EPS_C0 = 0.002


def mander_elastic_modulus(fc0: float, *, MPa: float = 1.0e6) -> float:
    """Módulo tangente inicial ``E_c = 5000·√f'c`` MPa (Mander et al. 1988,
    ec. 7). ``MPa`` es el valor de 1 MPa en las unidades de esfuerzo del
    usuario (1e6 con esfuerzos en Pa)."""
    fc0 = _positive(fc0, "fc0")
    return 5000.0 * MPa * math.sqrt(fc0 / MPa)


def fib_fracture_energy(f_cm: float, *, MPa: float = 1.0e6, N_per_m: float = 1.0) -> float:
    """Energía de fractura en tracción del fib Model Code 2010:
    ``G_F = 73·f_cm^0.18`` N/m, con ``f_cm`` (resistencia media a compresión)
    en MPa. ``MPa`` y ``N_per_m`` son esas unidades en las del usuario."""
    f_cm = _positive(f_cm, "f_cm")
    return 73.0 * (f_cm / MPa) ** 0.18 * N_per_m


def confinement_effectiveness_circular(s_clear: float, d_s: float, rho_cc: float,
                                       *, spiral: bool = False) -> float:
    """Coeficiente de efectividad ``k_e`` de una sección circular (Mander et
    al. 1988): con aros, ec. 14, ``(1 − s'/2d_s)²/(1 − ρ_cc)``; con
    espiral, ec. 15, ``(1 − s'/2d_s)/(1 − ρ_cc)``. ``s_clear`` es la
    separación libre ``s'`` entre aros o espiras, ``d_s`` el diámetro de la
    espiral entre ejes de barra y ``rho_cc`` la cuantía longitudinal
    referida al núcleo."""
    s_clear = _nonnegative(s_clear, "s_clear")
    d_s = _positive(d_s, "d_s")
    rho_cc = _ratio(rho_cc, "rho_cc")
    base = 1.0 - s_clear / (2.0 * d_s)
    if base <= 0.0:
        raise ValueError(
            f"confinement_effectiveness_circular: s_clear = {s_clear} ≥ 2·d_s; "
            f"los aros están tan separados que no confinan.")
    return (base if spiral else base * base) / (1.0 - rho_cc)


def confinement_effectiveness_rectangular(b_c: float, d_c: float, s_clear: float,
                                          w_clear: Sequence[float], rho_cc: float) -> float:
    """Coeficiente de efectividad ``k_e`` de una sección rectangular con
    estribos (Mander et al. 1988, ec. 22)::

        k_e = (1 − Σ w'_i²/(6 b_c d_c)) (1 − s'/2b_c) (1 − s'/2d_c) / (1 − ρ_cc)

    ``b_c`` y ``d_c`` son las dimensiones del núcleo entre ejes del estribo
    perimetral, ``s_clear`` la separación libre entre estribos, ``w_clear``
    los anchos libres entre barras longitudinales adyacentes y ``rho_cc``
    la cuantía longitudinal referida al núcleo."""
    b_c = _positive(b_c, "b_c")
    d_c = _positive(d_c, "d_c")
    s_clear = _nonnegative(s_clear, "s_clear")
    rho_cc = _ratio(rho_cc, "rho_cc")
    w = [float(_nonnegative(x, "w_clear")) for x in w_clear]
    k = ((1.0 - sum(x * x for x in w) / (6.0 * b_c * d_c))
         * (1.0 - s_clear / (2.0 * b_c)) * (1.0 - s_clear / (2.0 * d_c)) / (1.0 - rho_cc))
    if k <= 0.0:
        raise ValueError(
            f"confinement_effectiveness_rectangular: k_e = {k:.4g} ≤ 0; revise la "
            f"geometría (anchos libres o separación demasiado grandes).")
    return k


def confined_strength_ratio(fc0: float, fl1: float, fl2: float | None = None) -> float:
    """``K = f'cc/f'c`` a partir de las presiones laterales efectivas.

    - Una presión (``fl2`` omitida; sección circular): Mander et al. (1988),
      ec. 29, ``K = −1.254 + 2.254·√(1 + 7.94 q) − 2q`` con ``q = f'l/f'c``.
    - Dos presiones (sección rectangular, iguales o no): aproximación cerrada
      de Chang y Mander (1994), ecs. 3-70a–e, p. 3-32, al ábaco de Mander
      (Fig. 4). Con presiones iguales difiere de la ec. 29 en menos de 0.003
      en ``K``; se usa también en ese caso para que la función no tenga
      saltos. Validez ``f'l/f'c ≤ 0.3``: por encima, aviso.
    """
    fc0 = _positive(fc0, "fc0")
    q1 = _nonnegative(fl1, "fl1") / fc0
    if fl2 is None:
        if q1 > _CHANG_MANDER_MAX_RATIO:
            _warn_ratio(q1)
        return -1.254 + 2.254 * math.sqrt(1.0 + 7.94 * q1) - 2.0 * q1
    q2 = _nonnegative(fl2, "fl2") / fc0
    lo, hi = min(q1, q2), max(q1, q2)
    if hi > _CHANG_MANDER_MAX_RATIO:
        _warn_ratio(hi)
    if hi == 0.0:
        return 1.0
    r = lo / hi
    x = 0.5 * (lo + hi)
    A = 6.8886 - (0.6069 + 17.275 * r) * math.exp(-4.989 * r)
    B = 4.5 / ((5.0 / A) * (0.9849 - 0.6306 * math.exp(-3.8939 * r)) - 0.1) - 5.0
    return 1.0 + A * x * (0.1 + 0.9 / (1.0 + B * x))


def mander_confined_concrete(fc0: float, *, f_yh: float, eps_su: float, k_e,
                             rho_s: float | None = None,
                             rho_x: float | None = None, rho_y: float | None = None,
                             eps_c0: float = DEFAULT_EPS_C0, E: float | None = None,
                             MPa: float = 1.0e6) -> dict:
    """Parámetros del núcleo confinado para ``ManderConcrete1D``.

    Sección circular con ``rho_s`` (cuantía volumétrica del refuerzo
    transversal), o rectangular con ``rho_x`` y ``rho_y`` (Mander et al.
    1988, ecs. 23–24). ``k_e`` es el coeficiente de efectividad o un dict
    con la geometría para calcularlo: ``{s_clear, d_s, rho_cc, spiral}``
    (circular, ecs. 14–15) o ``{b_c, d_c, s_clear, w_clear, rho_cc}``
    (rectangular, ec. 22).

    Devuelve ``{"fc", "eps_c", "eps_cu", "E"}``:

    - ``fc = K·f'c`` con ``K`` de :func:`confined_strength_ratio`;
    - ``eps_c = ε_c0·[1 + 5(K − 1)]`` (ec. 5);
    - ``eps_cu = 0.004 + 1.4·ρ_s·f_yh·ε_su/f'cc`` (balance de energía en la
      aproximación de Paulay y Priestley 1992, §3.2.2), con
      ``ρ_s = ρ_x + ρ_y`` en la sección rectangular;
    - ``E`` el dado, o ``5000·√f'c`` MPa (ec. 7).
    """
    fc0 = _positive(fc0, "fc0")
    f_yh = _positive(f_yh, "f_yh")
    eps_su = _positive(eps_su, "eps_su")
    eps_c0 = _positive(eps_c0, "eps_c0")
    circular = rho_s is not None
    rectangular = rho_x is not None or rho_y is not None
    if circular == rectangular or (rectangular and (rho_x is None or rho_y is None)):
        raise ValueError(
            "mander_confined_concrete: declare 'rho_s' (sección circular) o "
            "'rho_x' y 'rho_y' (sección rectangular), no ambos ni ninguno.")
    k_e = _effectiveness(k_e, circular)
    if circular:
        rho_s = _ratio(rho_s, "rho_s")
        fl = k_e * 0.5 * rho_s * f_yh                       # ecs. 18-19
        K = confined_strength_ratio(fc0, fl)
        rho_vol = rho_s
    else:
        rho_x = _ratio(rho_x, "rho_x")
        rho_y = _ratio(rho_y, "rho_y")
        K = confined_strength_ratio(fc0, k_e * rho_x * f_yh, k_e * rho_y * f_yh)  # ecs. 27-28
        rho_vol = rho_x + rho_y
    fcc = K * fc0
    eps_cc = eps_c0 * (1.0 + 5.0 * (K - 1.0))
    eps_cu = 0.004 + 1.4 * rho_vol * f_yh * eps_su / fcc
    if eps_cu < eps_cc:
        raise ValueError(
            f"mander_confined_concrete: la deformación última del balance de energía, "
            f"ε_cu = {eps_cu:.4g}, es menor que la de pico confinada, ε_cc = {eps_cc:.4g}: "
            f"el refuerzo transversal se rompería antes del pico. Revise f_yh, eps_su y "
            f"las cuantías.")
    if E is None:
        E = mander_elastic_modulus(fc0, MPa=MPa)
    return {"fc": fcc, "eps_c": eps_cc, "eps_cu": eps_cu, "E": float(E)}


# ---------------------------------------------------------------------------
# Validaciones
# ---------------------------------------------------------------------------

def _effectiveness(k_e, circular: bool) -> float:
    if isinstance(k_e, dict):
        geometry = dict(k_e)
        try:
            if circular:
                return confinement_effectiveness_circular(**geometry)
            return confinement_effectiveness_rectangular(**geometry)
        except TypeError as exc:
            shape = "circular" if circular else "rectangular"
            raise ValueError(
                f"mander_confined_concrete: geometría de k_e no válida para la sección "
                f"{shape} ({exc}). Circular: s_clear, d_s, rho_cc, spiral; "
                f"rectangular: b_c, d_c, s_clear, w_clear, rho_cc.") from exc
    k_e = _positive(k_e, "k_e")
    if k_e > 1.0 + 1e-12:
        raise ValueError(f"mander_confined_concrete: k_e = {k_e} > 1.")
    return k_e


def _warn_ratio(q: float) -> None:
    _log.warning(
        f"confinamiento: f'l/f'c = {q:.3g} > {_CHANG_MANDER_MAX_RATIO}, fuera del "
        f"intervalo de validez de la aproximación de Chang y Mander (1994) y del "
        f"ábaco de Mander (1988).")


def _positive(x, name: str) -> float:
    x = float(x)
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError(f"{name} debe ser un número finito > 0 (recibido {x}).")
    return x


def _nonnegative(x, name: str) -> float:
    x = float(x)
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError(f"{name} debe ser un número finito ≥ 0 (recibido {x}).")
    return x


def _ratio(x, name: str) -> float:
    x = _nonnegative(x, name)
    if x >= 1.0:
        raise ValueError(f"{name} debe ser una cuantía < 1 (recibido {x}).")
    return x
