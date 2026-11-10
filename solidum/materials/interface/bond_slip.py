"""Leyes de adherencia-deslizamiento barra–concreto (materiales de interfaz de
línea, ADR 0033). Specs ``BondSlipModelCode`` y ``BondSlipMultilinear``.

Componentes en el sistema local de la barra, normales primero: ``(n, s)`` en
2D y ``(n1, n2, s)`` en 3D, con ``s`` el deslizamiento de la barra respecto
del concreto a lo largo de su eje.

- **Tangencial**: envolvente monótona ``τ_env(|s|)`` de cada ley, con la
  misma forma al jalar y al empujar.
- **Normal**: ``t_n = k_n·Δu_n``, elástica y desacoplada (DIANA §44.5).
- **Descarga y recarga (R1)**: geometría de Eligehausen, Popov y Bertero
  (1983, Fig. 5.3b) sin degradación. Predictor ``τ* = τ_n + k0·(s − s_n)``
  acotado por

      U(s) = max(τ_fr, min(τ_env(s)·[s > 0], ℓ⁺(s)))
      L(s) = −max(τ_fr, min(τ_env(−s)·[s < 0], −ℓ⁻(s)))

  con ``ℓ±`` las rectas de pendiente ``k0`` por el deslizamiento extremo de
  cada lado y ``τ_fr = min(τ_f, τ_env(s_máx))`` la fricción.
"""
from __future__ import annotations

import math

import numpy as np

from solidum.core.interface_material import InterfaceMaterial
from solidum.registry import InterfaceMaterialRegistry

# Valores por omisión de la función auxiliar del Model Code 2010 (decisiones
# del usuario del 2026-11-10): rigidez inicial y de descarga de Eligehausen,
# Popov y Bertero (1983, §5.3.4) = 180 N/mm³.
MC2010_DEFAULT_K0 = 1.8e11          # Pa/m
_ALPHA_MC2010 = 0.4
_MPA = 1.0e6
_MM = 1.0e-3


class _BondSlipLaw(InterfaceMaterial):
    """Regla R1 común; cada ley aporta ``envelope(a) -> (τ, dτ/da)`` para
    ``a ≥ 0`` y define ``k0``, ``tau_f`` y ``k_n``."""

    GEOMETRY = "line"
    JUMP_DIM = None
    PRIMARY_STATE_VAR = "slip_max"
    IS_SYMMETRIC = True

    k0: float
    tau_f: float
    k_n: float

    def envelope(self, a: float) -> tuple[float, float]:     # pragma: no cover
        raise NotImplementedError

    # -- envolvente con signo -------------------------------------------------

    def _env_pos(self, s: float) -> tuple[float, float]:
        """``τ_env(s)·[s > 0]`` y su derivada respecto de ``s``."""
        if s <= 0.0:
            return 0.0, 0.0
        return self.envelope(s)

    # -- ley ------------------------------------------------------------------

    def compute_traction(self, jump, state_vars=None):
        jump = np.asarray(jump, dtype=np.float64).reshape(-1)
        n = jump.size
        if n not in (2, 3):
            raise ValueError(
                f"{type(self).__name__}: el desplazamiento relativo debe tener 2 (2D) "
                f"o 3 (3D) componentes; tiene {n}.")
        st = state_vars or {}
        sp, tp = st.get("s_pos", 0.0), st.get("tau_pos", 0.0)
        sn, tn = st.get("s_neg", 0.0), st.get("tau_neg", 0.0)
        s_old, tau_old = st.get("s", 0.0), st.get("tau", 0.0)

        s = float(jump[-1])
        tau, dtau = self._tangential(s, s_old, tau_old, sp, tp, sn, tn)
        if s > sp:
            sp, tp = s, tau
        if s < sn:
            sn, tn = s, tau

        traction = np.empty(n)
        traction[:-1] = self.k_n * jump[:-1]
        traction[-1] = tau
        tangent = np.zeros((n, n))
        tangent[np.arange(n - 1), np.arange(n - 1)] = self.k_n
        tangent[-1, -1] = dtau
        new_state = {"s_pos": sp, "tau_pos": tp, "s_neg": sn, "tau_neg": tn,
                     "s": s, "tau": tau, "slip_max": max(sp, -sn)}
        return traction, tangent, new_state

    def _tangential(self, s, s_old, tau_old, sp, tp, sn, tn):
        """``(τ, dτ/ds)`` de la regla R1 desde el estado confirmado."""
        k0 = self.k0
        slip_max = max(sp, -sn)
        tau_fr = min(self.tau_f, self.envelope(slip_max)[0]) if slip_max > 0.0 else 0.0

        # Cota superior: rama de fricción, envolvente positiva o recta ℓ⁺.
        env_p, denv_p = self._env_pos(s)
        if sp > 0.0:
            line_p = tp + k0 * (s - sp)
            if env_p <= line_p:
                inner_u, dinner_u = env_p, denv_p
            else:
                inner_u, dinner_u = line_p, k0
        else:
            inner_u, dinner_u = env_p, denv_p
        if tau_fr > inner_u:
            upper, dupper = tau_fr, 0.0
        else:
            upper, dupper = inner_u, dinner_u

        # Cota inferior, en espejo.
        env_m, denv_m = self._env_pos(-s)      # τ_env(−s)·[s < 0], derivada en −s
        if sn < 0.0:
            line_m = tn + k0 * (s - sn)        # ℓ⁻(s), negativa en su rama útil
            if env_m <= -line_m:
                inner_l, dinner_l = env_m, -denv_m        # d/ds de τ_env(−s)
            else:
                inner_l, dinner_l = -line_m, -k0
        else:
            inner_l, dinner_l = env_m, -denv_m
        if tau_fr > inner_l:
            lower, dlower = -tau_fr, 0.0
        else:
            lower, dlower = -inner_l, -dinner_l

        # Una cota actúa sólo si el predictor la rebasa: en un empate (el estado
        # virgen en s = 0, donde las dos cotas y la fricción valen cero) la
        # respuesta es la elástica, con tangente k0.
        trial = tau_old + k0 * (s - s_old)
        if trial > upper:
            return upper, dupper
        if trial < lower:
            return lower, dlower
        return trial, k0

    # -- validación compartida -------------------------------------------------

    @staticmethod
    def _positive(name, value, owner):
        if not (isinstance(value, (int, float)) and math.isfinite(value) and value > 0.0):
            raise ValueError(f"{owner}: {name} debe ser > 0 (recibido {value!r}).")
        return float(value)


@InterfaceMaterialRegistry.register
class BondSlipModelCode(_BondSlipLaw):
    """Adherencia-deslizamiento con la envolvente de Eligehausen, Popov y
    Bertero (1983) / fib Model Code 2010 (spec ``BondSlipModelCode``).

    ``τ_env``: tramo lineal ``k0·a`` hasta ``s0``, ``τ_max·(a/s1)^α`` hasta
    ``s1``, meseta hasta ``s2``, descenso lineal hasta ``τ_f`` en ``s3`` y
    residual ``τ_f``; ``s0 = s1·(k0·s1/τ_max)^(1/(α − 1))``.

    Parameters
    ----------
    tau_max, s1, s2, s3, alpha, tau_f
        Envolvente [Pa, m, m, m, —, Pa].
    k0 : float
        Rigidez inicial y de descarga [Pa/m], ``k0 > τ_max/s1``.
    k_n : float
        Rigidez normal [Pa/m].

    La función auxiliar :meth:`model_code_2010` da los parámetros de la
    Tabla 6.1-1 del Model Code 2010.
    """

    # Clave del YAML que traduce from_yaml (los argumentos de model_code_2010).
    YAML_KEYS = ("model_code_2010",)

    def __init__(self, tau_max: float, s1: float, s2: float, s3: float, alpha: float,
                 tau_f: float, k0: float, k_n: float):
        name = "BondSlipModelCode"
        self.tau_max = self._positive("tau_max", tau_max, name)
        self.s1 = self._positive("s1", s1, name)
        if not s2 >= s1:
            raise ValueError(f"{name}: s2 debe ser ≥ s1 (s1 = {s1!r}, s2 = {s2!r}).")
        if not s3 >= s2:
            raise ValueError(f"{name}: s3 debe ser ≥ s2 (s2 = {s2!r}, s3 = {s3!r}).")
        if not 0.0 < alpha < 1.0:
            raise ValueError(f"{name}: alpha debe estar en (0, 1) (recibido {alpha!r}).")
        if not 0.0 <= tau_f <= tau_max:
            raise ValueError(
                f"{name}: tau_f debe estar en [0, tau_max] (tau_f = {tau_f!r}, "
                f"tau_max = {tau_max!r}).")
        self.k0 = self._positive("k0", k0, name)
        if not k0 > tau_max / s1:
            raise ValueError(
                f"{name}: k0 = {k0:.4g} Pa/m debe ser > tau_max/s1 = {tau_max / s1:.4g} "
                f"Pa/m para que el tramo lineal corte la rama potencial antes de s1.")
        self.k_n = self._positive("k_n", k_n, name)
        self.s2, self.s3 = float(s2), float(s3)
        self.alpha, self.tau_f = float(alpha), float(tau_f)
        self.s0 = self.s1 * (self.k0 * self.s1 / self.tau_max) ** (1.0 / (self.alpha - 1.0))

    def envelope(self, a: float) -> tuple[float, float]:
        if a <= self.s0:
            return self.k0 * a, self.k0
        if a <= self.s1:
            tau = self.tau_max * (a / self.s1) ** self.alpha
            return tau, self.alpha * tau / a
        if a <= self.s2:
            return self.tau_max, 0.0
        if a <= self.s3:
            slope = -(self.tau_max - self.tau_f) / (self.s3 - self.s2)
            return self.tau_max + slope * (a - self.s2), slope
        return self.tau_f, 0.0

    @classmethod
    def model_code_2010(cls, fcm: float, *, bond: str = "good", failure: str = "pullout",
                        stirrups: bool = False, clear_rib_spacing: float | None = None,
                        diameter: float | None = None, E_c: float | None = None,
                        k0: float = MC2010_DEFAULT_K0,
                        k_n: float | None = None) -> "BondSlipModelCode":
        """Ley con los parámetros de la Tabla 6.1-1 del fib Model Code 2010.

        Parameters
        ----------
        fcm : float
            Resistencia media a compresión del concreto en cilindro [Pa].
        bond : {"good", "other"}
            Condición de adherencia.
        failure : {"pullout", "splitting"}
            Modo de falla: arrancamiento (concreto confinado) o hendimiento.
        stirrups : bool
            En hendimiento, con estribos (``True``) o sin confinar.
        clear_rib_spacing : float
            Separación libre entre corrugas ``c`` [m]: ``s3`` en arrancamiento,
            ``0.5·c`` en hendimiento con estribos.
        diameter, E_c : float
            ``d_b`` [m] y módulo del concreto [Pa], para ``k_n = 2·E_c/d_b``
            (Murcia-Delso y Shing 2016) si no se da ``k_n``.
        k0 : float
            Rigidez inicial [Pa/m]; 180 N/mm³ por omisión (Eligehausen et al.
            1983, §5.3.4).
        """
        if bond not in ("good", "other"):
            raise ValueError(f"model_code_2010: bond debe ser 'good' u 'other' (recibido {bond!r}).")
        if failure not in ("pullout", "splitting"):
            raise ValueError(
                f"model_code_2010: failure debe ser 'pullout' o 'splitting' (recibido {failure!r}).")
        f = cls._positive("fcm", fcm, "model_code_2010") / _MPA
        good = bond == "good"
        tau_po = (2.5 if good else 1.25) * math.sqrt(f) * _MPA
        s1_po = (1.0 if good else 1.8) * _MM
        s2_po = (2.0 if good else 3.6) * _MM
        needs_c = failure == "pullout" or stirrups
        if needs_c:
            c = cls._positive("clear_rib_spacing", clear_rib_spacing, "model_code_2010")
        if failure == "pullout":
            tau_max, s1, s2, s3 = tau_po, s1_po, s2_po, c
            tau_f = 0.40 * tau_max
        else:
            coef = {(True, False): 7.0, (True, True): 8.0,
                    (False, False): 5.0, (False, True): 5.5}[(good, bool(stirrups))]
            tau_max = coef * (f / 25.0) ** 0.25 * _MPA
            s1 = s1_po * (tau_max / tau_po) ** (1.0 / _ALPHA_MC2010)
            s2 = s1
            s3 = 0.5 * c if stirrups else 1.2 * s1
            tau_f = 0.40 * tau_max if stirrups else 0.0
        if k_n is None:
            if diameter is None or E_c is None:
                raise ValueError(
                    "model_code_2010: sin k_n hacen falta diameter (d_b) y E_c para "
                    "k_n = 2·E_c/d_b.")
            k_n = 2.0 * cls._positive("E_c", E_c, "model_code_2010") / cls._positive(
                "diameter", diameter, "model_code_2010")
        return cls(tau_max=tau_max, s1=s1, s2=s2, s3=s3, alpha=_ALPHA_MC2010,
                   tau_f=tau_f, k0=k0, k_n=k_n)

    @classmethod
    def from_yaml(cls, data: dict, resolve) -> "BondSlipModelCode":
        """Parámetros explícitos, o un bloque ``model_code_2010`` con los
        argumentos de :meth:`model_code_2010` (en unidades SI)."""
        data = dict(data)
        mc = data.pop("model_code_2010", None)
        if mc is None:
            return cls(**data)
        if data:
            raise ValueError(
                f"BondSlipModelCode: con model_code_2010 los demás parámetros van dentro "
                f"del bloque; sobran {sorted(data)}.")
        mc = dict(mc)
        fcm = mc.pop("fcm", None)
        if fcm is None:
            raise ValueError("BondSlipModelCode: model_code_2010 necesita fcm.")
        return cls.model_code_2010(fcm, **mc)


@InterfaceMaterialRegistry.register
class BondSlipMultilinear(_BondSlipLaw):
    """Adherencia-deslizamiento con envolvente multilineal dada por puntos
    (spec ``BondSlipMultilinear``): interpolación de ``(0, 0), (s_1, τ_1), …,
    (s_m, τ_m)`` y ``τ_m`` más allá. ``k0 = τ_1/s_1``; ``τ_f`` por omisión
    ``τ_m``.

    Parameters
    ----------
    points : list of (s [m], τ [Pa])
    tau_f : float, optional
    k_n : float
    """

    def __init__(self, points, k_n: float, tau_f: float | None = None):
        name = "BondSlipMultilinear"
        pts = np.asarray(points, dtype=np.float64)
        if pts.ndim != 2 or pts.shape[1] != 2 or pts.shape[0] < 1:
            raise ValueError(f"{name}: points debe ser una lista de al menos un par (s, τ).")
        s, tau = pts[:, 0], pts[:, 1]
        if not (s[0] > 0.0 and np.all(np.diff(s) > 0.0)):
            raise ValueError(f"{name}: los deslizamientos deben ser > 0 y crecientes; recibidos {s.tolist()}.")
        if np.any(tau < 0.0) or not tau[0] > 0.0:
            raise ValueError(f"{name}: las tracciones deben ser ≥ 0 y la primera > 0; recibidas {tau.tolist()}.")
        self.s_points = np.concatenate([[0.0], s])
        self.tau_points = np.concatenate([[0.0], tau])
        self.k0 = float(tau[0] / s[0])
        tau_f = float(tau[-1]) if tau_f is None else tau_f
        if not 0.0 <= tau_f <= float(tau.max()):
            raise ValueError(f"{name}: tau_f debe estar en [0, max τ_k] (recibido {tau_f!r}).")
        self.tau_f = float(tau_f)
        self.k_n = self._positive("k_n", k_n, name)

    def envelope(self, a: float) -> tuple[float, float]:
        sp, tp = self.s_points, self.tau_points
        if a >= sp[-1]:
            return float(tp[-1]), 0.0
        k = int(np.searchsorted(sp, a, side="left"))
        k = max(k, 1)
        slope = (tp[k] - tp[k - 1]) / (sp[k] - sp[k - 1])
        return float(tp[k - 1] + slope * (a - sp[k - 1])), float(slope)
