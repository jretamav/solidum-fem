# solidum_fem/solidum/materials/damage_tension_compression_2d.py
"""``TensionCompressionDamage2D`` — daño de tracción y de compresión 2D con cierre de grieta.

Modelo 2 de la familia de concreto: ley de esfuerzo de Faria, Oliver y Cervera
(1998) con dos daños escalares, ``σ = (1 − d⁺)·σ̄⁺ + (1 − d⁻)·σ̄⁻`` sobre la
partición espectral del esfuerzo efectivo ``σ̄ = C_e·ε``. Una grieta abierta
degrada sólo la parte de tracción; al cerrarse, la compresión recupera la
rigidez intacta. El daño de tracción es el de la banda de grieta del modelo 1
(Rankine, ``α(h)``, implícita o IMPL-EX). Ver
``docs/specs/TensionCompressionDamage2D.md`` y ``_tension_compression``.
"""
import math

import numpy as np
from numba import njit

from solidum.constants import ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL
from solidum.core.material import Material
from solidum.materials._crack_band import (
    CRACK_BAND_DAMAGE_MAX,
    DEFAULT_INTEGRATION,
    NO_STEP_RATIO,
    crack_band_alpha,
    parse_integration,
    validate_crack_band_parameters,
)
from solidum.materials._spectral_split import spectral_split_2d
from solidum.materials._tension_compression import (
    DC,
    DEFAULT_BETA,
    DT,
    ERR,
    H,
    KC,
    KC_PREV,
    KT,
    KT_PREV,
    STATE_SCHEMA_TC,
    TensionCompressionMixin,
    compression_damage_update,
    compression_measure,
    compression_ultimate,
    tension_damage_update,
)
from solidum.materials.elastic_2d import Elastic2D
from solidum.registry import MaterialRegistry


@njit(cache=True)
def _tc_2d_point(strain, S_old, E, kappa_0, alpha, Ce, tol, cap, integration, ratio,
                 f_c, G_c, K, nu, plane_strain, tol_c):
    """Un punto de integración, implícito o IMPL-EX (spec §3, §6–§8).

    ``Ce`` es la elasticidad plana de la hipótesis (esfuerzo plano o
    deformación plana). ``f_c > 0`` activa el daño de compresión (``G_c``,
    ``K`` de Drucker-Prager, ``tol_c`` su tolerancia de admisibilidad).
    Devuelve ``(σ, C_tan, fila de estado nueva)``.

    En deformación plana ``σ̄_zz = ν·(σ̄_xx + σ̄_yy)`` es principal y no se
    acopla con el plano en ``Q⁺``; para la medida de Rankine no supera al
    mayor del plano cuando éste es positivo (``0 ≤ ν < 0.5``, spec §6), así
    que el esfuerzo y ``C^d`` sólo necesitan el plano. La medida de
    compresión sí lo ve: ``σ̄⁻_zz = min(σ̄_zz, 0)`` entra en la traza y en
    ``J₂``, y ``∂d⁻/∂ε`` lleva el término ``g_zz·(1 − H(σ̄_zz))·ν·(C_e[0] +
    C_e[1])`` (``∂σ̄_zz/∂ε = λ·[1, 1, 0]``, spec §7).

    Caché de Numba: este núcleo se compila con ``spectral_split_2d`` y las
    funciones de ``_tension_compression`` dentro; un cambio en
    ``_spectral_split.py`` o en ``_tension_compression.py`` no invalida la
    caché de este archivo, así que al cambiarlos hay que tocar también éste
    (y el 3D).
    """
    sig_eff = Ce @ strain
    sp, Q, s1, m1 = spectral_split_2d(sig_eff)
    eps_eq = s1 / E if s1 > 0.0 else 0.0
    d_t, kt, kt_prev, err, dd = tension_damage_update(
        eps_eq, S_old[KT], S_old[KT_PREV], S_old[DT], S_old[ERR], E, kappa_0, alpha,
        tol, cap, integration, ratio)
    d_c = S_old[DC]
    kc = S_old[KC]
    kc_prev = S_old[KC_PREV]
    ddc = 0.0
    sm = sig_eff - sp
    szz_neg = False
    g = np.zeros(6)
    if f_c > 0.0:
        szz_m = 0.0
        if plane_strain:
            szz = nu * (sig_eff[0] + sig_eff[1])
            if szz <= 0.0:                  # H(0) = 0: σ̄_zz nulo cuenta como no positivo
                szz_m = szz
                szz_neg = True
        eps_c, g = compression_measure(sm[0], sm[1], szz_m, sm[2], 0.0, 0.0, E, K)
        kappa_u = compression_ultimate(S_old[H], E, f_c, G_c)
        d_c, kc, kc_prev, err_c, ddc = compression_damage_update(
            eps_c, S_old[KC], S_old[KC_PREV], S_old[DC], S_old[ERR], E, f_c, kappa_u,
            tol_c, cap, integration, ratio)
        if err_c > err:
            err = err_c
    # σ = (1 − d⁻)·σ̄ − (d⁺ − d⁻)·σ̄⁺ y, con los daños fijos,
    # C^d = (1 − d⁻)·C_e − (d⁺ − d⁻)·Q⁺·C_e (spec §7).
    sigma = (1.0 - d_c) * sig_eff - (d_t - d_c) * sp
    C_tan = (1.0 - d_c) * Ce - (d_t - d_c) * (Q @ Ce)
    if dd > 0.0:
        # Carga activa implícita: − σ̄⁺ ⊗ ∂d⁺/∂ε, ∂d⁺/∂ε = d⁺'·C_e·(n_I⊗n_I)/E
        # (en Voigt, n_I⊗n_I con el factor 2 en el cortante).
        p = m1.copy()
        p[2] *= 2.0
        v = Ce @ p
        C_tan = C_tan - (dd / E) * np.outer(sp, v)
    if ddc > 0.0:
        # − σ̄⁻ ⊗ ∂d⁻/∂ε, ∂d⁻/∂ε = d⁻'·(∂ε⁻_eq/∂σ̄⁻)·(I − Q⁺)·C_e (spec §7): el
        # gradiente respecto del vector Voigt lleva el factor 2 en el cortante.
        gv = np.empty(3)
        gv[0] = g[0]
        gv[1] = g[1]
        gv[2] = 2.0 * g[3]
        v = (gv - gv @ Q) @ Ce
        if szz_neg:
            for b in range(3):
                v[b] += g[2] * nu * (Ce[0, b] + Ce[1, b])
        C_tan = C_tan - ddc * np.outer(sm, v)
    row = np.empty(8)
    row[KT] = kt
    row[KC] = kc
    row[DT] = d_t
    row[DC] = d_c
    row[H] = S_old[H]
    row[KT_PREV] = kt_prev
    row[KC_PREV] = kc_prev
    row[ERR] = err
    return sigma, C_tan, row


@njit(cache=True)
def _tc_2d_batch(strain, S_old, S_new, params, C, sigma, C_out, flag):
    """Adaptador por lotes (ADR 0014). ``params = [E, κ_0, f_t, G_f, tol_abs,
    tol_rel, damage_max, integración, r, f_c, G_c, K, ν, deformación plana]``
    (``f_c = 0`` sin aplastamiento); ``C`` = ``C_e`` plana; estado de
    ``_tension_compression``. Daño sin ``h`` sembrada → ``flag``
    (``batch_report``)."""
    E = params[0]
    kappa_0 = params[1]
    h = S_old[H]
    alpha = crack_band_alpha(h, E, params[2], params[3])
    tol = params[4] + params[5] * (E * kappa_0)
    f_c = params[9]
    tol_c = params[4] + params[5] * (f_c / 3.0)
    sig, C_tan, row = _tc_2d_point(strain, S_old, E, kappa_0, alpha, C, tol, params[6],
                                   int(params[7]), params[8], f_c, params[10], params[11],
                                   params[12], params[13] > 0.0, tol_c)
    if h <= 0.0 and (row[DT] > 0.0 or row[ERR] > 0.0 or (f_c > 0.0 and row[DC] > 0.0)):
        flag[0] = 1
    for i in range(3):
        sigma[i] = sig[i]
    for i in range(8):
        S_new[i] = row[i]
    return C_tan


@MaterialRegistry.register
class TensionCompressionDamage2D(TensionCompressionMixin, Material):
    """
    Daño de tracción y de compresión 2D con cierre de grieta, banda de grieta e IMPL-EX.

    Parameters
    ----------
    E, nu : float
        Parámetros elásticos del medio intacto, ``0 ≤ ν < 0.5``.
    f_t : float
        Resistencia a tracción; ``κ_0 = f_t/E``.
    G_f : float
        Energía de fractura en tracción (fuerza/longitud); el punto disipa
        ``G_f/h`` por unidad de volumen en tracción uniaxial.
    hypothesis : str
        ``'plane_stress'`` (default) o ``'plane_strain'``.
    f_c, G_c : float, optional
        Resistencia a compresión uniaxial (positiva, ``f_c > f_t``) y energía
        de aplastamiento (fuerza/longitud). Van juntas: activan el daño de
        compresión ``d⁻`` (Drucker-Prager sobre la parte negativa del
        esfuerzo efectivo y parábola de Feenstra, con ``G_c/h`` en la rama
        descendente). Sin ellas la compresión es elástica.
    beta : float
        ``f_bc/f_c``, razón entre las resistencias biaxial y uniaxial en
        compresión, ``β ≥ 1``; 1.16 por omisión (Kupfer).
    damage_max : float
        Tope de los daños, en (0, 1); por omisión ``1 − 10⁻⁶``.
    integration : str
        ``'implicit'`` (default) o ``'impl-ex'`` (Oliver, Huespe y Cante
        2008; ADR 0034): daños extrapolados de los dos últimos pasos. A
        diferencia del daño isótropo, con los daños fijos el material no es
        lineal (la partición cambia con la deformación): el Newton itera con
        la tangente de daño fijo.
    density : float, optional
        Densidad (ADR 0008).

    Notes
    -----
    ``σ = (1 − d⁻)·σ̄ − (d⁺ − d⁻)·σ̄⁺`` con ``σ̄ = C_e·ε`` el esfuerzo efectivo
    del estado 3D asociado (``σ̄_zz = 0`` en esfuerzo plano; ``ε_zz = 0`` en
    deformación plana, donde ``σ̄_zz`` es principal y ``σ_zz = (1 − d±)·σ̄_zz``
    según su signo). Tracción uniaxial exacta, ``σ = (1 − d⁺)·E·ε``; una
    grieta cerrada recupera ``E`` y ``ν``. ``d⁺`` con Rankine sobre ``σ̄`` y
    ley exponencial regularizada por la longitud característica ``h`` del
    punto, ``h < 2·E·G_f/f_t²`` (modelo 1). ``d⁻`` (con ``f_c`` y ``G_c``)
    con la medida de Drucker-Prager sobre ``σ̄⁻`` del estado 3D asociado y la
    parábola de Feenstra: pico ``f_c`` en compresión uniaxial, ``β·f_c`` en
    la biaxial igual, sin daño en la hidrostática y ``G_c/h`` en la rama
    descendente, sin límite de tamaño (spec §3). Con ``ν ≠ 0`` la ley no
    deriva de un potencial y la tangente no es simétrica, tampoco con
    IMPL-EX (spec §1, §7, §9). Estado: ``kappa_t``, ``kappa_c``,
    ``damage_t``, ``damage_c``, ``h``, ``kappa_t_prev``, ``kappa_c_prev`` y
    ``damage_error`` (``_tension_compression``).
    """
    # Descarga por la secante hasta el origen (daños fijos: σ(s·ε) = s·σ(ε)):
    # ½·σ:ε es el trabajo que devuelve una descarga radial (ADR 0038).
    SECANT_UNLOADING = True
    STRAIN_DIM = 3
    PRIMARY_STATE_VAR = 'damage_t'
    IS_SYMMETRIC = False
    STATE_SCHEMA = STATE_SCHEMA_TC

    def __init__(self, E: float, nu: float, f_t: float, G_f: float,
                 hypothesis: str = 'plane_stress',
                 f_c: float | None = None, G_c: float | None = None,
                 beta: float = DEFAULT_BETA,
                 damage_max: float = CRACK_BAND_DAMAGE_MAX,
                 integration: str = DEFAULT_INTEGRATION,
                 density: float | None = None):
        owner = type(self).__name__
        validate_crack_band_parameters(owner, E, nu, f_t, G_f, damage_max, density)
        if float(nu) < 0.0:
            raise ValueError(
                f"{owner}: nu = {nu} < 0 no se admite: la ley se formula para "
                f"0 ≤ ν < 0.5 (con ν < 0, en deformación plana σ̄_zz podría superar "
                f"al mayor esfuerzo principal del plano y la medida de Rankine "
                f"dejaría de ser la del estado 3D).")
        self._integration_code = parse_integration(owner, integration)
        self._ratio = NO_STEP_RATIO
        self.E = float(E)
        self.nu = float(nu)
        self.f_t = float(f_t)
        self.G_f = float(G_f)
        self.kappa_0 = self.f_t / self.E
        self.damage_max = float(damage_max)
        self._set_compression(f_c, G_c, beta)
        self.elastic_base = Elastic2D(E, nu, hypothesis=hypothesis)
        self.density = density

    @property
    def hypothesis(self) -> str:
        return self.elastic_base.hypothesis

    def compute_state(self, strain: np.ndarray, state_vars=None):
        S_old = self._state_row(state_vars)
        h = S_old[H]
        strain = np.ascontiguousarray(strain, dtype=np.float64)
        alpha = crack_band_alpha(h, self.E, self.f_t, self.G_f)
        f_c, G_c, K, tol_c = self._compression_params()
        sigma, C_tan, row = _tc_2d_point(
            strain, S_old, self.E, self.kappa_0, alpha, self.elastic_base.C,
            self.admissibility_tol(state_vars), self.damage_max,
            self._integration_code, self._ratio, f_c, G_c, K, self.nu,
            self.hypothesis == 'plane_strain', tol_c)
        if h <= 0.0 and (row[DT] > 0.0 or row[ERR] > 0.0 or (f_c > 0.0 and row[DC] > 0.0)):
            raise self._missing_h()
        return sigma, C_tan, self._row_to_state(row)

    BATCH_KERNEL = _tc_2d_batch

    def batch_params(self) -> np.ndarray:
        f_c, G_c, K, _ = self._compression_params()
        return np.array([self.E, self.kappa_0, self.f_t, self.G_f,
                         ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL, self.damage_max,
                         float(self._integration_code), self._ratio, f_c, G_c, K, self.nu,
                         1.0 if self.hypothesis == 'plane_strain' else 0.0],
                        dtype=np.float64)

    def batch_matrix(self) -> np.ndarray:
        return self.elastic_base.C

    # --- componente fuera del plano (post-proceso) -----------------------

    def out_of_plane_stress(self, sigma, state_vars=None) -> float:
        """``σ_zz``: nula en esfuerzo plano; en deformación plana,
        ``(1 − d±)·σ̄_zz`` según el signo de ``σ̄_zz = ν·tr σ̄`` del plano.

        ``σ`` y ``σ̄`` tienen las mismas direcciones principales y
        ``σ_i = (1 − d±)·σ̄_i`` con el signo de ``σ̄_i`` (``1 − d± > 0``), así
        que ``σ̄`` se recupera de ``σ`` y de los daños sin la deformación.
        """
        if self.hypothesis != 'plane_strain':
            return 0.0
        d_t = 0.0 if state_vars is None else float(state_vars.get('damage_t', 0.0))
        d_c = 0.0 if state_vars is None else float(state_vars.get('damage_c', 0.0))
        sx, sy, sxy = (float(sigma[0]), float(sigma[1]), float(sigma[2]))
        r = math.hypot(0.5 * (sx - sy), sxy)
        mean = 0.5 * (sx + sy)
        tr_eff = 0.0
        for s in (mean + r, mean - r):
            tr_eff += s / (1.0 - d_t) if s > 0.0 else s / (1.0 - d_c)
        szz_eff = self.nu * tr_eff
        return (1.0 - d_t) * szz_eff if szz_eff > 0.0 else (1.0 - d_c) * szz_eff

    def batch_out_of_plane_stress(self, sigma, S):
        P = int(sigma.shape[0])
        if self.hypothesis != 'plane_strain':
            return np.zeros(P)
        sx, sy, sxy = sigma[:, 0], sigma[:, 1], sigma[:, 2]
        d_t, d_c = S[:, DT], S[:, DC]
        r = np.hypot(0.5 * (sx - sy), sxy)
        mean = 0.5 * (sx + sy)
        tr_eff = np.zeros(P)
        for s in (mean + r, mean - r):
            tr_eff += np.where(s > 0.0, s / (1.0 - d_t), s / (1.0 - d_c))
        szz_eff = self.nu * tr_eff
        return np.where(szz_eff > 0.0, (1.0 - d_t) * szz_eff, (1.0 - d_c) * szz_eff)
