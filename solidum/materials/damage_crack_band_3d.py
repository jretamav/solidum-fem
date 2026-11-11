# solidum_fem/solidum/materials/damage_crack_band_3d.py
"""``IsotropicDamageCrackBand3D`` — daño isótropo 3D regularizado por banda de grieta.

Versión 3D de :class:`IsotropicDamageCrackBand2D` y variante de
:class:`IsotropicDamage3D` (que no cambia): ``α`` se calcula en cada punto
con ``G_f`` y la longitud característica ``h`` del elemento; medida de
Rankine por omisión, tope de daño propio e integración implícita o IMPL-EX
(ADR 0034). Voigt 6D del proyecto (ADR 0012).
Ver ``docs/specs/IsotropicDamageCrackBand3D.md``.
"""
import numpy as np
from numba import njit

from solidum.constants import ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL
from solidum.materials._crack_band import (
    CRACK_BAND_DAMAGE_MAX,
    CRACK_BAND_EQUIVALENT_STRAINS,
    DEFAULT_CRACK_BAND_EQUIVALENT_STRAIN,
    DEFAULT_INTEGRATION,
    IMPL_EX,
    NO_STEP_RATIO,
    CrackBandMixin,
    crack_band_alpha,
    implex_update,
    parse_integration,
    validate_crack_band_parameters,
)
from solidum.materials._equivalent_strain import (
    ENERGY,
    equivalent_strain_3d,
    parse_equivalent_strain,
)
from solidum.materials.damage_3d import IsotropicDamage3D, _compute_damage_3d
from solidum.materials.elastic_3d import Elastic3D
from solidum.registry import MaterialRegistry


@njit(cache=True)
def _crack_band_3d_point(strain, S_old, E, kappa_0, alpha, Ce, tol, cap, nu,
                         measure, k, integration, ratio):
    """Un punto de integración, implícito o IMPL-EX. Estado
    ``[κ, d, h, κ_prev, e]``. Devuelve ``(σ, C_tan, κ, d, κ_prev, e)``."""
    if integration == IMPL_EX:
        sigma_eff = Ce @ strain
        if ratio < 0.0:
            eps_eq = 0.0  # fuera de un paso: no hace falta (estado congelado)
        else:
            eps_eq, _, _ = equivalent_strain_3d(strain, sigma_eff, E, nu, measure, k)
        d, kappa_new, kappa_prev, err = implex_update(
            eps_eq, S_old[0], S_old[3], S_old[1], S_old[4], ratio,
            E, kappa_0, alpha, tol, cap)
        return (1.0 - d) * sigma_eff, (1.0 - d) * Ce, kappa_new, d, kappa_prev, err
    sigma, C_tan, kappa_new, d = _compute_damage_3d(
        strain, S_old[0], E, kappa_0, alpha, Ce, tol, cap, nu, measure, k)
    return sigma, C_tan, kappa_new, d, S_old[0], 0.0


@njit(cache=True)
def _crack_band_3d_batch(strain, S_old, S_new, params, C, sigma, C_out, flag):
    """Adaptador por lotes (ADR 0014). ``params = [E, κ_0, f_t, G_f, tol_abs,
    tol_rel, damage_max, ν, código de la medida, k, integración, r]``; estado
    ``[κ, d, h, κ_prev, e]``. Daño sin ``h`` sembrada → ``flag``
    (``batch_report``)."""
    E = params[0]
    kappa_0 = params[1]
    h = S_old[2]
    alpha = crack_band_alpha(h, E, params[2], params[3])
    tol = params[4] + params[5] * (E * kappa_0)
    sig, C_tan, kappa_new, d, kappa_prev, err = _crack_band_3d_point(
        strain, S_old, E, kappa_0, alpha, C, tol, params[6],
        params[7], int(params[8]), params[9], int(params[10]), params[11],
    )
    if h <= 0.0 and (d > 0.0 or err > 0.0):
        flag[0] = 1
    for i in range(6):
        sigma[i] = sig[i]
    S_new[0] = kappa_new
    S_new[1] = d
    S_new[2] = h
    S_new[3] = kappa_prev
    S_new[4] = err
    return C_tan


@MaterialRegistry.register
class IsotropicDamageCrackBand3D(CrackBandMixin, IsotropicDamage3D):
    """
    Daño isótropo 3D con ablandamiento exponencial regularizado por banda de grieta.

    Parameters
    ----------
    E, nu : float
        Parámetros elásticos del medio intacto.
    f_t : float
        Resistencia a tracción; ``κ_0 = f_t/E``.
    G_f : float
        Energía de fractura en modo I; el punto disipa ``G_f/h`` por unidad
        de volumen.
    equivalent_strain : str
        ``'rankine'`` (default), ``'energy'``, ``'mazars'`` o
        ``'modified_von_mises'`` (este último con ``k``).
    k : float, optional
        ``f_c/f_t``, sólo y obligatorio con ``'modified_von_mises'``.
    damage_max : float
        Tope del daño, en (0, 1); por omisión ``1 − 10⁻⁶``.
    integration : str
        ``'implicit'`` (default) o ``'impl-ex'`` (ver
        :class:`IsotropicDamageCrackBand2D`; ADR 0034).
    density : float, optional
        Densidad (ADR 0008).

    Notes
    -----
    Ver :class:`IsotropicDamageCrackBand2D`; la variante 2D es la
    restricción exacta de ésta en tensión plana y deformación plana.
    """
    STRAIN_DIM = 6
    PRIMARY_STATE_VAR = 'damage'
    IS_SYMMETRIC = False
    STATE_SCHEMA = {'kappa': (), 'damage': (), 'h': (), 'kappa_prev': (), 'damage_error': ()}

    def __init__(self, E: float, nu: float, f_t: float, G_f: float,
                 equivalent_strain: str = DEFAULT_CRACK_BAND_EQUIVALENT_STRAIN,
                 k: float | None = None,
                 damage_max: float = CRACK_BAND_DAMAGE_MAX,
                 integration: str = DEFAULT_INTEGRATION,
                 density: float | None = None):
        owner = type(self).__name__
        validate_crack_band_parameters(owner, E, nu, f_t, G_f, damage_max, density)
        self.equivalent_strain, self._measure_code, self.k = parse_equivalent_strain(
            owner, equivalent_strain, k, allowed=CRACK_BAND_EQUIVALENT_STRAINS)
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
        self.alpha = None
        self.IS_SYMMETRIC = (self._measure_code == ENERGY
                             or self._integration_code == IMPL_EX)

    def compute_state(self, strain: np.ndarray, state_vars=None):
        S_old = self._state_row(state_vars)
        h = S_old[2]
        strain = np.ascontiguousarray(strain, dtype=np.float64)
        alpha = crack_band_alpha(h, self.E, self.f_t, self.G_f)
        sigma, C_tan, kappa_new, d, kappa_prev, err = _crack_band_3d_point(
            strain, S_old, self.E, self.kappa_0, alpha,
            self.elastic_base.C, self.admissibility_tol(state_vars), self.damage_max,
            self.nu, self._measure_code, self._k_param(),
            self._integration_code, self._ratio,
        )
        if h <= 0.0 and (d > 0.0 or err > 0.0):
            raise self._missing_h()
        return sigma, C_tan, {'kappa': kappa_new, 'damage': d, 'h': h,
                              'kappa_prev': kappa_prev, 'damage_error': err}

    BATCH_KERNEL = _crack_band_3d_batch

    def batch_params(self) -> np.ndarray:
        return np.array([self.E, self.kappa_0, self.f_t, self.G_f,
                         ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL, self.damage_max,
                         self.nu, float(self._measure_code), self._k_param(),
                         float(self._integration_code), self._ratio],
                        dtype=np.float64)
