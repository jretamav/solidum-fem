# solidum_fem/solidum/materials/damage_2d.py
"""``IsotropicDamage2D`` — daño escalar isótropo 2D con ablandamiento exponencial.

Soporta tangente **algorítmica consistente** durante carga activa (recupera
convergencia cuadrática del Newton global). En descarga, por debajo del umbral
o al saturar el daño, devuelve tangente secante. La deformación equivalente es
configurable (``equivalent_strain``: energética, Mazars o von Mises
modificado; ver ``solidum.materials._equivalent_strain``). Ver
``docs/specs/IsotropicDamage2D.md``.
"""
import numpy as np
from numba import njit

from solidum.constants import ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL, DAMAGE_MAX
from solidum.core.material import Material
from solidum.materials._equivalent_strain import (
    DEFAULT_EQUIVALENT_STRAIN,
    ENERGY,
    equivalent_strain_2d,
    parse_equivalent_strain,
)
from solidum.materials._softening import evaluate_exponential_damage
from solidum.materials.elastic_2d import Elastic2D
from solidum.registry import MaterialRegistry


@njit(cache=True)
def _compute_damage_2d(strain, kappa_old, E, kappa_0, alpha, Ce, tol, cap, c_zz,
                       nu, measure, k):
    """Actualización del daño isótropo 2D en un punto de Gauss.

    Devuelve ``(σ, C_tan, κ_new, d)``. ``tol`` es la tolerancia de
    admisibilidad ya evaluada (ADR 0006), ``cap`` el techo del daño
    (``DAMAGE_MAX``) y ``c_zz`` el coeficiente de la deformación fuera del
    plano, ``ε_zz = c_zz·(ε_xx + ε_yy)``: ``−ν/(1−ν)`` en tensión plana y
    ``0`` en deformación plana. ``measure`` es el código de la deformación
    equivalente (``_equivalent_strain``) y ``k = f_c/f_t`` el parámetro de
    von Mises modificado (sin uso en las otras medidas). Es el único sitio
    donde vive la física del modelo: lo llaman ``compute_state`` (camino por
    elemento) y el adaptador por lotes ``_damage_2d_batch`` (ADR 0014).
    """
    # Esfuerzo efectivo del medio intacto y deformación equivalente, con
    # ∂ε_eq/∂ε = s·v. El 2D es la restricción exacta de IsotropicDamage3D
    # con cualquiera de las tres medidas (ε_zz de arriba, γ_yz = γ_xz = 0).
    sigma_eff = Ce @ strain
    eps_eq, v, s = equivalent_strain_2d(strain, sigma_eff, E, nu, c_zz, measure, k)

    # Régimen de carga vs descarga (Kuhn-Tucker)
    if eps_eq > kappa_old:
        kappa_new = eps_eq
        loading = True
    else:
        kappa_new = kappa_old
        loading = False

    # Check de admisibilidad en esfuerzo (ADR 0006): f = E·(κ − κ_0)
    f_stress = E * (kappa_new - kappa_0)
    if f_stress <= tol:
        d = 0.0
    else:
        d, _ = evaluate_exponential_damage(kappa_new, kappa_0, alpha, cap)
    # Flag de saturación consistente al final del bloque.
    saturated = (d >= cap) and f_stress > 0.0

    C_sec = (1.0 - d) * Ce
    sigma = C_sec @ strain

    # Tangente: consistente solo en carga activa con daño efectivo no saturado.
    # En todos los otros casos (descarga, sin daño, saturación), secante.
    consistent_branch = loading and (d > 0.0) and (not saturated) and (eps_eq > 0.0)
    if consistent_branch:
        # C_alg = (1-d)·C_e − (∂d/∂κ)·(C_e·ε) ⊗ (∂ε_eq/∂ε)
        # ∂d/∂κ = (1-d)·(1/κ + α);  ∂ε_eq/∂ε = s·v
        # Con la medida energética v = C_e·ε: el producto externo es
        # simétrico bit a bit.
        dd_dkappa = (1.0 - d) * (1.0 / kappa_new + alpha)
        C_tan = C_sec - (dd_dkappa * s) * np.outer(sigma_eff, v)
    else:
        C_tan = C_sec
    return sigma, C_tan, kappa_new, d


@njit(cache=True)
def _damage_2d_batch(strain, S_old, S_new, params, C, sigma, C_out, flag):
    """Adaptador por lotes (ADR 0014). ``params = [E, κ_0, α, tol_abs,
    tol_rel, DAMAGE_MAX, c_zz, ν, código de la medida, k]``; ``C`` = ``C_e``
    del medio intacto; estado ``[κ, d]``."""
    E = params[0]
    kappa_0 = params[1]
    tol = params[3] + params[4] * (E * kappa_0)
    sig, C_tan, kappa_new, d = _compute_damage_2d(
        strain, S_old[0], E, kappa_0, params[2], C, tol, params[5], params[6],
        params[7], int(params[8]), params[9],
    )
    for i in range(3):
        sigma[i] = sig[i]
    S_new[0] = kappa_new
    S_new[1] = d
    return C_tan


@MaterialRegistry.register
class IsotropicDamage2D(Material):
    """
    Modelo de daño isótropo continuo 2D con ablandamiento exponencial.

    Tangente algorítmica consistente en carga activa; tangente secante
    ``(1-d)·C_e`` en descarga, sin daño activo o al saturar.

    Parameters
    ----------
    E, nu : float
        Parámetros elásticos del medio intacto.
    kappa_0 : float
        Umbral inicial de deformación equivalente. ``d = 0`` mientras ``κ ≤ κ_0``.
        Las tres medidas dan ``ε_eq = ε`` en tracción uniaxial de esfuerzo, así
        que ``κ_0 = f_t/E``.
    alpha : float
        Velocidad de degradación exponencial (``α > 0``). Mayor α → ablandamiento
        más abrupto tras el umbral.
    hypothesis : str
        ``'plane_stress'`` (default) o ``'plane_strain'``. Delegado a ``Elastic2D``
        interno.
    equivalent_strain : str
        Medida de deformación equivalente: ``'energy'`` (default, Simó y Ju
        1987), ``'mazars'`` (Mazars 1986) o ``'modified_von_mises'`` (de Vree,
        Brekelmans y van Gils 1995).
    k : float, optional
        ``f_c/f_t`` (≥ 1; con k < 1 el material no dañaría en tracción triaxial). Obligatorio con ``'modified_von_mises'`` y no
        admitido con las otras dos medidas.
    density : float, optional
        Densidad (ADR 0008).

    Notes
    -----
    Deformación equivalente (``_equivalent_strain``), sobre el estado 3D
    ``[ε_xx, ε_yy, ε_zz, γ_xy, 0, 0]`` con ``ε_zz = −ν/(1−ν)·(ε_xx + ε_yy)``
    en tensión plana y ``ε_zz = 0`` en deformación plana:

    - ``'energy'``: ``ε_eq = √(εᵀ·C_e·ε/E)`` con la ``C_e`` plana del material
      (en tensión plana, ``σ:ε`` del 3D con ``σ_zz = 0``). No distingue
      tracción de compresión.
    - ``'mazars'``: ``ε_eq = √(Σ⟨ε_i⟩₊²)`` sobre las tres deformaciones
      principales (las dos del plano y ``ε_zz``). En compresión uniaxial de
      esfuerzo, ``ε_eq = √2·ν·|ε|``.
    - ``'modified_von_mises'``: función de ``I₁`` y ``J₂`` con ``k``; en
      compresión uniaxial de esfuerzo, ``ε_eq = |ε|/k``.

    Con cualquiera de las tres el 2D es la restricción exacta de
    ``IsotropicDamage3D``. Hasta el 2026-09-29 la medida era la norma del
    tensor ``√(ε:ε)``, retirada por no estándar y mal calibrada
    (``ε_eq = ε·√(1+2ν²)`` en tracción uniaxial).

    Simetría de la tangente (ADR 0003): la consistente es
    ``(1-d)·C_e − (∂d/∂κ)·(C_e·ε) ⊗ ∂ε_eq/∂ε``. Con la medida energética
    ``∂ε_eq/∂ε ∝ C_e·ε`` y es **simétrica**; con Mazars y von Mises modificado
    no lo es en general. ``IS_SYMMETRIC`` es ``False`` en la clase (cota
    conservadora) y la instancia lo pone a ``True`` con ``'energy'``;
    ``domain_is_symmetric`` lo lee de la instancia.

    Cap ``DAMAGE_MAX`` aplica al daño escalar ``d`` y al esfuerzo nominal
    ``σ = (1-d)·C_e·ε`` — semántica de daño **continuo**.
    """
    STRAIN_DIM = 3
    PRIMARY_STATE_VAR = 'damage'
    # Cota conservadora de clase; la instancia la sobreescribe con 'energy'.
    IS_SYMMETRIC = False
    STATE_SCHEMA = {'kappa': (), 'damage': ()}

    def __init__(self, E: float, nu: float, kappa_0: float, alpha: float,
                 hypothesis: str = 'plane_stress',
                 equivalent_strain: str = DEFAULT_EQUIVALENT_STRAIN,
                 k: float | None = None,
                 density: float | None = None):
        if E <= 0.0:
            raise ValueError(f"IsotropicDamage2D: E debe ser > 0 (recibido {E}).")
        if kappa_0 <= 0.0:
            raise ValueError(f"IsotropicDamage2D: kappa_0 debe ser > 0 (recibido {kappa_0}).")
        if alpha <= 0.0:
            raise ValueError(f"IsotropicDamage2D: alpha debe ser > 0 (recibido {alpha}).")
        if density is not None and density < 0.0:
            raise ValueError(
                f"IsotropicDamage2D: density={density} no puede ser negativa."
            )
        self.equivalent_strain, self._measure_code, self.k = parse_equivalent_strain(
            "IsotropicDamage2D", equivalent_strain, k)

        self.E = float(E)
        self.elastic_base = Elastic2D(E, nu, hypothesis=hypothesis)
        self.kappa_0 = float(kappa_0)
        self.alpha = float(alpha)
        self.density = density
        # ε_zz = c_zz·(ε_xx + ε_yy) en la deformación equivalente.
        self.c_zz = (-float(nu) / (1.0 - float(nu))
                     if self.elastic_base.hypothesis == 'plane_stress' else 0.0)
        # Override por instancia (ADR 0003): tangente consistente simétrica
        # sólo con la medida energética.
        self.IS_SYMMETRIC = self._measure_code == ENERGY

    def admissibility_scale(self, state_vars=None) -> float:
        """Esfuerzo umbral inicial ``E·κ_0`` (ADR 0006).

        Constante respecto al estado: el daño no introduce una superficie
        cuya escala evolucione — el umbral siempre se mide contra ``κ_0``
        aunque ``κ`` histórica crezca.
        """
        return self.E * self.kappa_0

    def initial_state(self) -> dict:
        """``κ`` arranca en ``κ_0`` (es lo que ``compute_state`` asume con
        ``state_vars=None``); el daño en cero."""
        return {'kappa': float(self.kappa_0), 'damage': 0.0}

    def _k_param(self) -> float:
        """``k`` para el núcleo compilado (sin uso fuera de von Mises modificado)."""
        return 0.0 if self.k is None else float(self.k)

    def compute_state(self, strain: np.ndarray, state_vars=None):
        kappa_old = self.kappa_0 if state_vars is None else state_vars.get('kappa', self.kappa_0)
        strain = np.ascontiguousarray(strain, dtype=np.float64)
        sigma, C_tan, kappa_new, d = _compute_damage_2d(
            strain, float(kappa_old), self.E, self.kappa_0, self.alpha,
            self.elastic_base.C, self.admissibility_tol(state_vars), DAMAGE_MAX,
            self.c_zz, self.elastic_base.nu, self._measure_code, self._k_param(),
        )
        new_state = {'kappa': kappa_new, 'damage': d}
        return sigma, C_tan, new_state

    # Camino por lotes (ADR 0014)
    BATCH_KERNEL = _damage_2d_batch

    def batch_params(self) -> np.ndarray:
        return np.array([self.E, self.kappa_0, self.alpha,
                         ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL, DAMAGE_MAX,
                         self.c_zz, self.elastic_base.nu,
                         float(self._measure_code), self._k_param()],
                        dtype=np.float64)

    def batch_matrix(self) -> np.ndarray:
        return self.elastic_base.C

    def out_of_plane_stress(self, sigma, state_vars=None) -> float:
        """Delegado al elástico base: ``σ = (1−d)·C_e·ε`` escala todas las
        componentes por igual, así que ``σ_zz = ν·(σ_xx + σ_yy)`` sigue
        valiendo en plane strain con cualquier nivel de daño."""
        return self.elastic_base.out_of_plane_stress(sigma, state_vars)

    def batch_out_of_plane_stress(self, sigma, S):
        return self.elastic_base.batch_out_of_plane_stress(sigma, S)
