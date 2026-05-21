# solidum_fem/solidum/materials/damage_3d.py
"""``IsotropicDamage3D`` — daño escalar isótropo 3D con ablandamiento exponencial.

Soporta tangente **algorítmica consistente** durante carga activa (recupera
convergencia cuadrática del Newton global). En descarga, por debajo del umbral
o al saturar el daño, devuelve tangente secante. Convención Voigt 6D del
proyecto (ADR 0012).

Ver ``docs/specs/IsotropicDamage3D.md``.
"""
import math

import numpy as np

from solidum.constants import DAMAGE_MAX
from solidum.core.material import Material
from solidum.materials._softening import evaluate_exponential_damage
from solidum.materials.elastic_3d import Elastic3D
from solidum.registry import MaterialRegistry


@MaterialRegistry.register
class IsotropicDamage3D(Material):
    """
    Modelo de daño isótropo continuo 3D con ablandamiento exponencial.

    Tangente algorítmica consistente en carga activa (asimétrica en general);
    tangente secante ``(1-d)·C_e`` en descarga, sin daño activo o al saturar.

    Parameters
    ----------
    E, nu : float
        Parámetros elásticos del medio intacto. ``ν ∈ (-1, 0.5)`` estricto.
    kappa_0 : float
        Umbral inicial de deformación equivalente. ``d = 0`` mientras ``κ ≤ κ_0``.
    alpha : float
        Velocidad de degradación exponencial (``α > 0``). Mayor α → ablandamiento
        más abrupto tras el umbral.
    density : float, optional
        Densidad (ADR 0008). Opcional al construir; obligatoria si se ensambla
        peso propio o matriz de masa.

    Notes
    -----
    Convención Voigt 6D del proyecto (ADR 0012):
    ``ε = [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`` con ``γ_ij = 2·ε_ij``.

    Deformación equivalente simétrica
    ``ε_eq = √(ε_xx² + ε_yy² + ε_zz² + ½(γ_xy² + γ_yz² + γ_xz²))``
    (norma de Frobenius del tensor deformación 3D; extensión natural de la
    convención 2D del proyecto): no distingue tracción de compresión. Para
    hormigón con resistencia distinta en tracción/compresión se necesitaría
    un split tipo Mazars o de Vree (out-of-scope).

    La tangente algorítmica consistente no es simétrica en estados generales
    de deformación porque el producto externo ``(C_e·ε) ⊗ (M·ε)`` con
    ``M = diag(1, 1, 1, 1/2, 1/2, 1/2)`` no produce una matriz simétrica.
    Esto fuerza ``IS_SYMMETRIC = False`` y el despachador algebraico usa LU
    para el sistema global (ADR 0003).

    Cap ``DAMAGE_MAX`` aplica al daño escalar ``d`` y al esfuerzo nominal
    ``σ = (1-d)·C_e·ε`` — semántica de daño **continuo**. Esto contrasta con
    la familia paralela ``CohesiveMaterial`` (cohesivos traction-jump), donde
    el cap se aplica **sólo a la rigidez tangente** y no a ``ω`` ni a la
    tracción física (ver memoria ``feedback_damage_max_cohesivos.md``).
    """
    STRAIN_DIM = 6
    PRIMARY_STATE_VAR = 'damage'
    IS_SYMMETRIC = False  # tangente consistente asimétrica en carga activa

    def __init__(self, E: float, nu: float, kappa_0: float, alpha: float,
                 density: float | None = None):
        if E <= 0.0:
            raise ValueError(f"IsotropicDamage3D: E debe ser > 0 (recibido {E}).")
        if kappa_0 <= 0.0:
            raise ValueError(
                f"IsotropicDamage3D: kappa_0 debe ser > 0 (recibido {kappa_0})."
            )
        if alpha <= 0.0:
            raise ValueError(f"IsotropicDamage3D: alpha debe ser > 0 (recibido {alpha}).")
        if not -1.0 < nu < 0.5:
            raise ValueError(
                f"IsotropicDamage3D: nu debe estar en (-1, 0.5) (recibido {nu}). "
                f"nu = 0.5 induce singularidad (K → ∞)."
            )
        if density is not None and density < 0.0:
            raise ValueError(
                f"IsotropicDamage3D: density={density} no puede ser negativa."
            )

        self.E = float(E)
        self.nu = float(nu)
        # Base elástica reusada por composición (mismo patrón que Damage2D).
        # Elastic3D valida nu en (-1, 0.5) por su cuenta también.
        self.elastic_base = Elastic3D(E=self.E, nu=self.nu)
        self.kappa_0 = float(kappa_0)
        self.alpha = float(alpha)
        self.density = density

    def admissibility_scale(self, state_vars=None) -> float:
        """Esfuerzo umbral inicial ``E·κ_0`` (ADR 0006).

        Constante respecto al estado: el daño no introduce una superficie cuya
        escala evolucione — el umbral siempre se mide contra ``κ_0`` aunque
        ``κ`` histórica crezca. Idéntico al caso 1D/2D.
        """
        return self.E * self.kappa_0

    def compute_state(self, strain: np.ndarray, state_vars=None):
        """Devuelve ``(σ, C_tan, new_state)``.

        Variables internas en ``state_vars``:
        - ``kappa``: float — máxima ε_eq histórica (default κ_0 si state_vars is None).
        - ``damage``: float ∈ [0, DAMAGE_MAX] — variable de daño.
        """
        kappa_old = self.kappa_0 if state_vars is None else state_vars.get('kappa', self.kappa_0)

        # Deformación equivalente simétrica en Voigt 6D:
        # ε_eq² = ε_xx² + ε_yy² + ε_zz² + ½(γ_xy² + γ_yz² + γ_xz²)
        eps_eq = math.sqrt(
            strain[0] ** 2 + strain[1] ** 2 + strain[2] ** 2
            + 0.5 * (strain[3] ** 2 + strain[4] ** 2 + strain[5] ** 2)
        )

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
        saturated = (d >= DAMAGE_MAX) and f_stress > 0.0

        Ce = self.elastic_base.C
        C_sec = (1.0 - d) * Ce
        sigma = C_sec @ strain

        # Tangente: consistente sólo en carga activa con daño efectivo no saturado.
        # Otros casos (descarga, sin daño, saturación) → secante.
        consistent_branch = loading and (d > 0.0) and (not saturated) and (eps_eq > 0.0)
        if consistent_branch:
            # C_alg = (1-d)·C_e - (∂d/∂κ)·(C_e·ε) ⊗ (M·ε / ε_eq)
            # ∂d/∂κ = (1-d)·(1/κ + α)
            # ∂κ/∂ε = (1/ε_eq)·M·ε  con M = diag(1,1,1,1/2,1/2,1/2)
            # σ_eff = C_e·ε (6-vector)
            dd_dkappa = (1.0 - d) * (1.0 / kappa_new + self.alpha)
            sigma_eff = Ce @ strain
            depseq_deps = np.array([
                strain[0] / eps_eq,
                strain[1] / eps_eq,
                strain[2] / eps_eq,
                0.5 * strain[3] / eps_eq,
                0.5 * strain[4] / eps_eq,
                0.5 * strain[5] / eps_eq,
            ])
            C_tan = C_sec - dd_dkappa * np.outer(sigma_eff, depseq_deps)
        else:
            C_tan = C_sec

        new_state = {'kappa': kappa_new, 'damage': d}
        return sigma, C_tan, new_state
