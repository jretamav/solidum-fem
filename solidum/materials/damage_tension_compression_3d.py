# solidum_fem/solidum/materials/damage_tension_compression_3d.py
"""``TensionCompressionDamage3D`` — daño de tracción y de compresión 3D con cierre de grieta.

Versión 3D de :class:`TensionCompressionDamage2D`: la misma ley de Faria,
Oliver y Cervera (1998) con la descomposición espectral 3×3 del esfuerzo
efectivo. Voigt 6D del proyecto (ADR 0012). Ver
``docs/specs/TensionCompressionDamage3D.md``.
"""
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
from solidum.materials._spectral_split import spectral_split_3d
from solidum.materials._tension_compression import (
    DC,
    DT,
    ERR,
    H,
    KC,
    KC_PREV,
    KT,
    KT_PREV,
    STATE_SCHEMA_TC,
    TensionCompressionMixin,
    tension_damage_update,
)
from solidum.materials.elastic_3d import Elastic3D
from solidum.registry import MaterialRegistry


@njit(cache=True)
def _tc_3d_point(strain, S_old, E, kappa_0, alpha, Ce, tol, cap, integration, ratio):
    """Un punto de integración, implícito o IMPL-EX (spec 2D §6–§8 con la
    descomposición 3×3). Devuelve ``(σ, C_tan, fila de estado nueva)``.

    Caché de Numba: se compila con ``spectral_split_3d`` y
    ``tension_damage_update`` dentro; un cambio en ``_spectral_split.py`` o
    en ``_tension_compression.py`` no invalida la caché de este archivo, así
    que al cambiarlos hay que tocar también éste (y el 2D).
    """
    sig_eff = Ce @ strain
    sp, Q, s1, m1 = spectral_split_3d(sig_eff)
    eps_eq = s1 / E if s1 > 0.0 else 0.0
    d_t, kt, kt_prev, err, dd = tension_damage_update(
        eps_eq, S_old[KT], S_old[KT_PREV], S_old[DT], S_old[ERR], E, kappa_0, alpha,
        tol, cap, integration, ratio)
    d_c = S_old[DC]
    sigma = (1.0 - d_c) * sig_eff - (d_t - d_c) * sp
    C_tan = (1.0 - d_c) * Ce - (d_t - d_c) * (Q @ Ce)
    if dd > 0.0:
        # − σ̄⁺ ⊗ ∂d⁺/∂ε, ∂d⁺/∂ε = d⁺'·C_e·(n_I⊗n_I)/E (factor 2 en el cortante).
        p = m1.copy()
        p[3] *= 2.0
        p[4] *= 2.0
        p[5] *= 2.0
        v = Ce @ p
        C_tan = C_tan - (dd / E) * np.outer(sp, v)
    row = np.empty(8)
    row[KT] = kt
    row[KC] = S_old[KC]
    row[DT] = d_t
    row[DC] = d_c
    row[H] = S_old[H]
    row[KT_PREV] = kt_prev
    row[KC_PREV] = S_old[KC_PREV]
    row[ERR] = err
    return sigma, C_tan, row


@njit(cache=True)
def _tc_3d_batch(strain, S_old, S_new, params, C, sigma, C_out, flag):
    """Adaptador por lotes (ADR 0014). ``params = [E, κ_0, f_t, G_f, tol_abs,
    tol_rel, damage_max, integración, r]``; ``C`` = ``C_e``; estado de
    ``_tension_compression``. Daño sin ``h`` sembrada → ``flag``."""
    E = params[0]
    kappa_0 = params[1]
    h = S_old[H]
    alpha = crack_band_alpha(h, E, params[2], params[3])
    tol = params[4] + params[5] * (E * kappa_0)
    sig, C_tan, row = _tc_3d_point(strain, S_old, E, kappa_0, alpha, C, tol, params[6],
                                   int(params[7]), params[8])
    if h <= 0.0 and (row[DT] > 0.0 or row[ERR] > 0.0):
        flag[0] = 1
    for i in range(6):
        sigma[i] = sig[i]
    for i in range(8):
        S_new[i] = row[i]
    return C_tan


@MaterialRegistry.register
class TensionCompressionDamage3D(TensionCompressionMixin, Material):
    """
    Daño de tracción y de compresión 3D con cierre de grieta, banda de grieta e IMPL-EX.

    Parameters
    ----------
    E, nu : float
        Parámetros elásticos del medio intacto, ``0 ≤ ν < 0.5``.
    f_t : float
        Resistencia a tracción; ``κ_0 = f_t/E``.
    G_f : float
        Energía de fractura en tracción; el punto disipa ``G_f/h`` por
        unidad de volumen en tracción uniaxial.
    damage_max : float
        Tope de los daños, en (0, 1); por omisión ``1 − 10⁻⁶``.
    integration : str
        ``'implicit'`` (default) o ``'impl-ex'`` (ver
        :class:`TensionCompressionDamage2D`; ADR 0034).
    density : float, optional
        Densidad (ADR 0008).

    Notes
    -----
    Ver :class:`TensionCompressionDamage2D`; la variante 2D es la
    restricción exacta de ésta en esfuerzo plano y deformación plana.
    """
    SECANT_UNLOADING = True
    STRAIN_DIM = 6
    PRIMARY_STATE_VAR = 'damage_t'
    IS_SYMMETRIC = False
    STATE_SCHEMA = STATE_SCHEMA_TC

    def __init__(self, E: float, nu: float, f_t: float, G_f: float,
                 damage_max: float = CRACK_BAND_DAMAGE_MAX,
                 integration: str = DEFAULT_INTEGRATION,
                 density: float | None = None):
        owner = type(self).__name__
        validate_crack_band_parameters(owner, E, nu, f_t, G_f, damage_max, density)
        if float(nu) < 0.0:
            raise ValueError(
                f"{owner}: nu = {nu} < 0 no se admite: la ley se formula para "
                f"0 ≤ ν < 0.5 (spec TensionCompressionDamage2D, validity).")
        self._integration_code = parse_integration(owner, integration)
        self._ratio = NO_STEP_RATIO
        self.E = float(E)
        self.nu = float(nu)
        self.f_t = float(f_t)
        self.G_f = float(G_f)
        self.kappa_0 = self.f_t / self.E
        self.damage_max = float(damage_max)
        self.elastic_base = Elastic3D(E=self.E, nu=self.nu)
        self.density = density

    def compute_state(self, strain: np.ndarray, state_vars=None):
        S_old = self._state_row(state_vars)
        h = S_old[H]
        strain = np.ascontiguousarray(strain, dtype=np.float64)
        alpha = crack_band_alpha(h, self.E, self.f_t, self.G_f)
        sigma, C_tan, row = _tc_3d_point(
            strain, S_old, self.E, self.kappa_0, alpha, self.elastic_base.C,
            self.admissibility_tol(state_vars), self.damage_max,
            self._integration_code, self._ratio)
        if h <= 0.0 and (row[DT] > 0.0 or row[ERR] > 0.0):
            raise self._missing_h()
        return sigma, C_tan, self._row_to_state(row)

    BATCH_KERNEL = _tc_3d_batch

    def batch_params(self) -> np.ndarray:
        return np.array([self.E, self.kappa_0, self.f_t, self.G_f,
                         ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL, self.damage_max,
                         float(self._integration_code), self._ratio],
                        dtype=np.float64)

    def batch_matrix(self) -> np.ndarray:
        return self.elastic_base.C
