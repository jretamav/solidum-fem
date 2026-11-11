# solidum_fem/solidum/materials/damage_crack_band_3d.py
"""``IsotropicDamageCrackBand3D`` — daño isótropo 3D regularizado por banda de grieta.

Versión 3D de :class:`IsotropicDamageCrackBand2D` y variante de
:class:`IsotropicDamage3D` (que no cambia): ``α`` se calcula en cada punto
con ``G_f`` y la longitud característica ``h`` del elemento; medida de
Rankine por omisión y tope de daño propio. Voigt 6D del proyecto (ADR 0012).
Ver ``docs/specs/IsotropicDamageCrackBand3D.md``.
"""
import numpy as np
from numba import njit

from solidum.constants import ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL
from solidum.materials._crack_band import (
    CRACK_BAND_DAMAGE_MAX,
    CRACK_BAND_EQUIVALENT_STRAINS,
    DEFAULT_CRACK_BAND_EQUIVALENT_STRAIN,
    CrackBandMixin,
    crack_band_alpha,
    validate_crack_band_parameters,
)
from solidum.materials._equivalent_strain import ENERGY, parse_equivalent_strain
from solidum.materials.damage_3d import IsotropicDamage3D, _compute_damage_3d
from solidum.materials.elastic_3d import Elastic3D
from solidum.registry import MaterialRegistry


@njit(cache=True)
def _crack_band_3d_batch(strain, S_old, S_new, params, C, sigma, C_out, flag):
    """Adaptador por lotes (ADR 0014). ``params = [E, κ_0, f_t, G_f, tol_abs,
    tol_rel, damage_max, ν, código de la medida, k]``; estado ``[κ, d, h]``.
    Daño sin ``h`` sembrada → ``flag`` (``batch_report``)."""
    E = params[0]
    kappa_0 = params[1]
    h = S_old[2]
    alpha = crack_band_alpha(h, E, params[2], params[3])
    tol = params[4] + params[5] * (E * kappa_0)
    sig, C_tan, kappa_new, d = _compute_damage_3d(
        strain, S_old[0], E, kappa_0, alpha, C, tol, params[6],
        params[7], int(params[8]), params[9],
    )
    if h <= 0.0 and d > 0.0:
        flag[0] = 1
    for i in range(6):
        sigma[i] = sig[i]
    S_new[0] = kappa_new
    S_new[1] = d
    S_new[2] = h
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
    STATE_SCHEMA = {'kappa': (), 'damage': (), 'h': ()}

    def __init__(self, E: float, nu: float, f_t: float, G_f: float,
                 equivalent_strain: str = DEFAULT_CRACK_BAND_EQUIVALENT_STRAIN,
                 k: float | None = None,
                 damage_max: float = CRACK_BAND_DAMAGE_MAX,
                 density: float | None = None):
        owner = type(self).__name__
        validate_crack_band_parameters(owner, E, nu, f_t, G_f, damage_max, density)
        self.equivalent_strain, self._measure_code, self.k = parse_equivalent_strain(
            owner, equivalent_strain, k, allowed=CRACK_BAND_EQUIVALENT_STRAINS)
        self.E = float(E)
        self.nu = float(nu)
        self.f_t = float(f_t)
        self.G_f = float(G_f)
        self.kappa_0 = self.f_t / self.E
        self.damage_max = float(damage_max)
        self.elastic_base = Elastic3D(E=self.E, nu=self.nu)
        self.density = density
        self.alpha = None
        self.IS_SYMMETRIC = self._measure_code == ENERGY

    def compute_state(self, strain: np.ndarray, state_vars=None):
        if state_vars is None:
            kappa_old, h = self.kappa_0, 0.0
        else:
            kappa_old = state_vars.get('kappa', self.kappa_0)
            h = float(state_vars.get('h', 0.0))
        strain = np.ascontiguousarray(strain, dtype=np.float64)
        alpha = crack_band_alpha(h, self.E, self.f_t, self.G_f)
        sigma, C_tan, kappa_new, d = _compute_damage_3d(
            strain, float(kappa_old), self.E, self.kappa_0, alpha,
            self.elastic_base.C, self.admissibility_tol(state_vars), self.damage_max,
            self.nu, self._measure_code, self._k_param(),
        )
        if h <= 0.0 and d > 0.0:
            raise self._missing_h()
        return sigma, C_tan, {'kappa': kappa_new, 'damage': d, 'h': h}

    BATCH_KERNEL = _crack_band_3d_batch

    def batch_params(self) -> np.ndarray:
        return np.array([self.E, self.kappa_0, self.f_t, self.G_f,
                         ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL, self.damage_max,
                         self.nu, float(self._measure_code), self._k_param()],
                        dtype=np.float64)
