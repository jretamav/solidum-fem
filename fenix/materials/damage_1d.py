# fenix_fem/fenix/materials/damage_1d.py
"""``IsotropicDamage1D`` — daño escalar 1D con ablandamiento exponencial.

Réplica escalar de :class:`fenix.materials.damage_2d.IsotropicDamage2D` con
tangente algorítmica consistente en carga activa. Ver
``docs/specs/IsotropicDamage1D.md``.
"""
import math

from fenix.constants import DAMAGE_MAX
from fenix.core.material import Material
from fenix.registry import MaterialRegistry


@MaterialRegistry.register
class IsotropicDamage1D(Material):
    """
    Modelo de daño isótropo continuo 1D con ablandamiento exponencial.

    Tangente algorítmica consistente en carga activa (``E_tan = -(1-d)·E·α·κ``,
    **negativa** — refleja el softening físico); tangente secante ``(1-d)·E`` en
    descarga, sin daño activo o al saturar.

    Parameters
    ----------
    E : float
        Módulo de Young intacto (`>0`).
    kappa_0 : float
        Umbral de deformación equivalente (``|ε|``) por debajo del cual no hay
        daño. `>0`.
    alpha : float
        Velocidad de degradación exponencial. `>0`; mayor α → ablandamiento
        más abrupto.
    density : float, optional
        Densidad (ADR 0008).

    Notes
    -----
    ``IS_SYMMETRIC`` queda `True` (default): la contribución elemental
    ``B^T·E_tan·B`` es simétrica aunque pueda ser indefinida cuando
    ``E_tan < 0``. El despachador algebraico (ADR 0003) detecta no-positividad
    y degrada a LU si Cholesky falla.

    Compatible con ``Truss2D`` y ``Truss3D``. Para seguir la rama post-pico de
    softening en problemas globales se requiere ``ArcLengthSolver``.
    """
    STRAIN_DIM = 1
    PRIMARY_STATE_VAR = 'damage'

    def __init__(self, E: float, kappa_0: float, alpha: float,
                 density: float | None = None):
        if E <= 0.0:
            raise ValueError(f"IsotropicDamage1D: E debe ser > 0 (recibido {E}).")
        if kappa_0 <= 0.0:
            raise ValueError(f"IsotropicDamage1D: kappa_0 debe ser > 0 (recibido {kappa_0}).")
        if alpha <= 0.0:
            raise ValueError(f"IsotropicDamage1D: alpha debe ser > 0 (recibido {alpha}).")

        self.E = E
        self.kappa_0 = kappa_0
        self.alpha = alpha
        self.density = density

    def admissibility_scale(self, state_vars=None) -> float:
        """Esfuerzo umbral inicial ``E·κ_0`` (ADR 0006).

        Constante respecto al estado: el daño no introduce una superficie cuya
        escala evolucione — el umbral siempre se mide contra ``κ_0`` aunque
        ``κ`` histórica crezca.
        """
        return self.E * self.kappa_0

    def compute_state(self, strain: float, state_vars=None):
        kappa_old = self.kappa_0 if state_vars is None else state_vars.get('kappa', self.kappa_0)

        eps_eq = abs(strain)

        # Régimen de carga (Kuhn-Tucker)
        if eps_eq > kappa_old:
            kappa_new = eps_eq
            loading = True
        else:
            kappa_new = kappa_old
            loading = False

        f_stress = self.E * (kappa_new - self.kappa_0)
        if self.is_admissible(f_stress, state_vars):
            d = 0.0
        else:
            d = 1.0 - (self.kappa_0 / kappa_new) * math.exp(-self.alpha * (kappa_new - self.kappa_0))
            if d >= DAMAGE_MAX:
                d = DAMAGE_MAX
        saturated = (d >= DAMAGE_MAX) and f_stress > 0.0

        sigma = (1.0 - d) * self.E * strain

        # Tangente: consistente solo en carga activa con daño efectivo no saturado.
        # En 1D: E_tan = (1-d)·E - σ_eff·(∂d/∂κ)·sign(ε) = -(1-d)·E·α·κ (negativa).
        consistent_branch = loading and (d > 0.0) and (not saturated)
        if consistent_branch:
            E_tan = -(1.0 - d) * self.E * self.alpha * kappa_new
        else:
            E_tan = (1.0 - d) * self.E

        new_state = {'kappa': kappa_new, 'damage': d}
        return sigma, E_tan, new_state
