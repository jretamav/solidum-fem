# solidum_fem/solidum/materials/damage_2d.py
"""``IsotropicDamage2D`` — daño escalar isótropo 2D con ablandamiento exponencial.

Soporta tangente **algorítmica consistente** durante carga activa (recupera
convergencia cuadrática del Newton global). En descarga, por debajo del umbral
o al saturar el daño, devuelve tangente secante. Ver ``docs/specs/IsotropicDamage2D.md``.
"""
import math

import numpy as np
from numba import njit

from solidum.constants import ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL, DAMAGE_MAX
from solidum.core.material import Material
from solidum.materials._softening import evaluate_exponential_damage
from solidum.materials.elastic_2d import Elastic2D
from solidum.registry import MaterialRegistry


@njit(cache=True)
def _compute_damage_2d(strain, kappa_old, E, kappa_0, alpha, Ce, tol, cap, c_zz):
    """Actualización del daño isótropo 2D en un punto de Gauss.

    Devuelve ``(σ, C_tan, κ_new, d)``. ``tol`` es la tolerancia de
    admisibilidad ya evaluada (ADR 0006), ``cap`` el techo del daño
    (``DAMAGE_MAX``) y ``c_zz`` el coeficiente de la deformación fuera del
    plano, ``ε_zz = c_zz·(ε_xx + ε_yy)``: ``−ν/(1−ν)`` en tensión plana y
    ``0`` en deformación plana. Es el único sitio donde vive la física del
    modelo: lo llaman ``compute_state`` (camino por elemento) y el adaptador
    por lotes ``_damage_2d_batch`` (ADR 0014).
    """
    # Deformación equivalente: norma del tensor de deformaciones, con la
    # componente fuera del plano ε_zz, para que el 2D sea la restricción
    # exacta de IsotropicDamage3D. En tensión plana ε_zz no depende del
    # daño: σ_zz = (1−d)·(λ·tr ε + 2μ·ε_zz) = 0 ⇒ ε_zz = −ν/(1−ν)·(ε_xx + ε_yy).
    eps_zz = c_zz * (strain[0] + strain[1])
    eps_eq = math.sqrt(strain[0] ** 2 + strain[1] ** 2 + eps_zz ** 2
                       + 0.5 * strain[2] ** 2)

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
        # C_alg = (1-d)·C_e - (∂d/∂κ)(∂κ/∂ε)·σ_eff^T
        # ∂d/∂κ = (1-d)·(1/κ + α)
        # ∂κ/∂ε = (1/ε_eq)·[ε_xx + c_zz·ε_zz, ε_yy + c_zz·ε_zz, ½γ_xy]
        # σ_eff = C_e·ε
        dd_dkappa = (1.0 - d) * (1.0 / kappa_new + alpha)
        sigma_eff = Ce @ strain
        depseq_deps = np.array([
            (strain[0] + c_zz * eps_zz) / eps_eq,
            (strain[1] + c_zz * eps_zz) / eps_eq,
            0.5 * strain[2] / eps_eq,
        ])
        C_tan = C_sec - dd_dkappa * np.outer(sigma_eff, depseq_deps)
    else:
        C_tan = C_sec
    return sigma, C_tan, kappa_new, d


@njit(cache=True)
def _damage_2d_batch(strain, S_old, S_new, params, C, sigma, C_out, flag):
    """Adaptador por lotes (ADR 0014). ``params = [E, κ_0, α, tol_abs,
    tol_rel, DAMAGE_MAX, c_zz]``; ``C`` = ``C_e`` del medio intacto; estado
    ``[κ, d]``."""
    E = params[0]
    kappa_0 = params[1]
    tol = params[3] + params[4] * (E * kappa_0)
    sig, C_tan, kappa_new, d = _compute_damage_2d(
        strain, S_old[0], E, kappa_0, params[2], C, tol, params[5], params[6]
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
    Deformación equivalente = norma del tensor de deformaciones,
    ``ε_eq = √(ε_xx² + ε_yy² + ε_zz² + ½γ_xy²)``, con ``ε_zz = −ν/(1−ν)·(ε_xx +
    ε_yy)`` en tensión plana y ``ε_zz = 0`` en deformación plana: el 2D es la
    restricción exacta de ``IsotropicDamage3D`` en ambas hipótesis (hasta el
    2026-09-29 la tensión plana omitía ``ε_zz`` y retrasaba el inicio del daño
    un factor ``√(1+2ν²)/√(1+ν²)``). No distingue tracción de compresión; para
    hormigón con resistencia diferente en tracción/compresión se necesitaría
    una medida tipo Mazars o de Vree (out-of-scope).

    La tangente algorítmica consistente no es simétrica en estados generales de
    deformación porque el producto externo ``(C_e·ε) ⊗ ∂ε_eq/∂ε`` no produce
    una matriz simétrica. Esto fuerza
    ``IS_SYMMETRIC = False`` y el despachador algebraico usa LU para el sistema
    global (ADR 0003).

    Cap ``DAMAGE_MAX`` aplica al daño escalar ``d`` y al esfuerzo nominal
    ``σ = (1-d)·C_e·ε`` — semántica de daño **continuo**.
    """
    STRAIN_DIM = 3
    PRIMARY_STATE_VAR = 'damage'
    IS_SYMMETRIC = False  # tangente consistente asimétrica en carga activa
    STATE_SCHEMA = {'kappa': (), 'damage': ()}

    def __init__(self, E: float, nu: float, kappa_0: float, alpha: float,
                 hypothesis: str = 'plane_stress', density: float | None = None):
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

        self.E = E
        self.elastic_base = Elastic2D(E, nu, hypothesis=hypothesis)
        self.kappa_0 = kappa_0
        self.alpha = alpha
        self.density = density
        # ε_zz = c_zz·(ε_xx + ε_yy) en la deformación equivalente.
        self.c_zz = (-float(nu) / (1.0 - float(nu))
                     if self.elastic_base.hypothesis == 'plane_stress' else 0.0)

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

    def compute_state(self, strain: np.ndarray, state_vars=None):
        kappa_old = self.kappa_0 if state_vars is None else state_vars.get('kappa', self.kappa_0)
        strain = np.ascontiguousarray(strain, dtype=np.float64)
        sigma, C_tan, kappa_new, d = _compute_damage_2d(
            strain, float(kappa_old), self.E, self.kappa_0, self.alpha,
            self.elastic_base.C, self.admissibility_tol(state_vars), DAMAGE_MAX,
            self.c_zz,
        )
        new_state = {'kappa': kappa_new, 'damage': d}
        return sigma, C_tan, new_state

    # Camino por lotes (ADR 0014)
    BATCH_KERNEL = _damage_2d_batch

    def batch_params(self) -> np.ndarray:
        return np.array([self.E, self.kappa_0, self.alpha,
                         ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL, DAMAGE_MAX,
                         self.c_zz],
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
