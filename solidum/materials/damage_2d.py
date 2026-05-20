# solidum_fem/solidum/materials/damage_2d.py
"""``IsotropicDamage2D`` — daño escalar isótropo 2D con ablandamiento exponencial.

Soporta tangente **algorítmica consistente** durante carga activa (recupera
convergencia cuadrática del Newton global). En descarga, por debajo del umbral
o al saturar el daño, devuelve tangente secante. Ver ``docs/specs/IsotropicDamage2D.md``.
"""
import numpy as np

from solidum.constants import DAMAGE_MAX
from solidum.core.material import Material
from solidum.materials._softening import evaluate_exponential_damage
from solidum.materials.elastic_2d import Elastic2D
from solidum.registry import MaterialRegistry


@MaterialRegistry.register
class IsotropicDamage2D(Material):
    """
    Modelo de daño isótropo continuo 2D con ablandamiento exponencial.

    Tangente algorítmica consistente en carga activa (asimétrica en general);
    tangente secante ``(1-d)·C_e`` en descarga, sin daño activo o al saturar.

    Parameters
    ----------
    E, nu : float
        Parámetros elásticos del medio intacto.
    kappa_0 : float
        Umbral inicial de deformación equivalente. ``d = 0`` mientras ``κ ≤ κ_0``.
    alpha : float
        Velocidad de degradación exponencial (``α > 0``). Mayor α → ablandamiento
        más abrupto tras el umbral.
    hypothesis : str
        ``'plane_stress'`` (default) o ``'plane_strain'``. Delegado a ``Elastic2D``
        interno.
    density : float, optional
        Densidad (ADR 0008).

    Notes
    -----
    Deformación equivalente simétrica ``ε_eq = √(ε_xx² + ε_yy² + ½γ_xy²)``: no
    distingue tracción de compresión. Para hormigón con resistencia diferente en
    tracción/compresión se necesitaría un split tipo Mazars o de Vree (out-of-scope).

    La tangente algorítmica consistente no es simétrica en estados generales de
    deformación porque el producto externo ``(C_e·ε) ⊗ (M·ε)`` con
    ``M = diag(1,1,1/2)`` no produce una matriz simétrica. Esto fuerza
    ``IS_SYMMETRIC = False`` y el despachador algebraico usa LU para el sistema
    global (ADR 0003).

    Cap ``DAMAGE_MAX`` aplica al daño escalar ``d`` y al esfuerzo nominal
    ``σ = (1-d)·C_e·ε`` — semántica de daño **continuo**. Esto contrasta con
    la familia paralela ``CohesiveMaterial`` (cohesivos traction-jump), donde
    el cap se aplica **sólo a la rigidez tangente** y no a ``ω`` ni a la
    tracción física, porque la escala del penalty ``K_e`` haría que un cap
    global desaparezca la tracción cohesiva real. Ver memoria
    ``feedback_damage_max_cohesivos.md`` para la justificación.
    """
    STRAIN_DIM = 3
    PRIMARY_STATE_VAR = 'damage'
    IS_SYMMETRIC = False  # tangente consistente asimétrica en carga activa

    def __init__(self, E: float, nu: float, kappa_0: float, alpha: float,
                 hypothesis: str = 'plane_stress', density: float | None = None):
        if E <= 0.0:
            raise ValueError(f"IsotropicDamage2D: E debe ser > 0 (recibido {E}).")
        if kappa_0 <= 0.0:
            raise ValueError(f"IsotropicDamage2D: kappa_0 debe ser > 0 (recibido {kappa_0}).")
        if alpha <= 0.0:
            raise ValueError(f"IsotropicDamage2D: alpha debe ser > 0 (recibido {alpha}).")

        self.E = E
        self.elastic_base = Elastic2D(E, nu, hypothesis=hypothesis)
        self.kappa_0 = kappa_0
        self.alpha = alpha
        self.density = density

    def admissibility_scale(self, state_vars=None) -> float:
        """Esfuerzo umbral inicial ``E·κ_0`` (ADR 0006).

        Constante respecto al estado: el daño no introduce una superficie
        cuya escala evolucione — el umbral siempre se mide contra ``κ_0``
        aunque ``κ`` histórica crezca.
        """
        return self.E * self.kappa_0

    def compute_state(self, strain: np.ndarray, state_vars=None):
        kappa_old = self.kappa_0 if state_vars is None else state_vars.get('kappa', self.kappa_0)

        # Deformación equivalente simétrica
        eps_eq = np.sqrt(strain[0]**2 + strain[1]**2 + 0.5 * strain[2]**2)

        # Régimen de carga vs descarga (Kuhn-Tucker)
        if eps_eq > kappa_old:
            kappa_new = eps_eq
            loading = True
        else:
            kappa_new = kappa_old
            loading = False

        # Check de admisibilidad en esfuerzo (ADR 0006): f = E·(κ − κ_0)
        f_stress = self.E * (kappa_new - self.kappa_0)
        if self.is_admissible(f_stress, state_vars):
            d = 0.0
        else:
            d, _ = evaluate_exponential_damage(
                kappa_new, self.kappa_0, self.alpha,
            )
        # Flag de saturación consistente al final del bloque.
        saturated = (d >= DAMAGE_MAX) and f_stress > 0.0

        Ce = self.elastic_base.C
        C_sec = (1.0 - d) * Ce
        sigma = C_sec @ strain

        # Tangente: consistente solo en carga activa con daño efectivo no saturado.
        # En todos los otros casos (descarga, sin daño, saturación), secante.
        consistent_branch = loading and (d > 0.0) and (not saturated) and (eps_eq > 0.0)
        if consistent_branch:
            # C_alg = (1-d)·C_e - (∂d/∂κ)(∂κ/∂ε)·σ_eff^T
            # ∂d/∂κ = (1-d)·(1/κ + α)
            # ∂κ/∂ε = (1/ε_eq)·M·ε  con M = diag(1, 1, 1/2)
            # σ_eff = C_e·ε
            dd_dkappa = (1.0 - d) * (1.0 / kappa_new + self.alpha)
            sigma_eff = Ce @ strain
            depseq_deps = np.array([
                strain[0] / eps_eq,
                strain[1] / eps_eq,
                0.5 * strain[2] / eps_eq,
            ])
            C_tan = C_sec - dd_dkappa * np.outer(sigma_eff, depseq_deps)
        else:
            C_tan = C_sec

        new_state = {'kappa': kappa_new, 'damage': d}
        return sigma, C_tan, new_state
