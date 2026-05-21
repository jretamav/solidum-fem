# solidum_fem/solidum/materials/von_mises_3d.py
"""``VonMises3D`` — plasticidad J2 con endurecimiento isótropo lineal en 3D.

Formulación íntegramente en Voigt 6D del proyecto (ADR 0012). Sin variantes
de hipótesis: en 3D todas las componentes de ``ε`` y ``σ`` son activas y el
return mapping es **único y radial cerrado** (Simó-Hughes 1998, §3.3).
Comparado con ``VonMises2D``:

- ``plane_strain`` 2D es la restricción de este modelo bajo
  ``ε_zz = γ_yz = γ_xz = 0`` impuestos en la deformación total (`ε^p_zz`,
  `ε^p_yz`, `ε^p_xz` libres por incompresibilidad/isotropía).
- ``plane_stress`` 2D usa un algoritmo proyectado genuinamente distinto
  (Newton local sobre ``Δγ``), por lo que no comparte kernel con VM3D.

Ver ``docs/specs/VonMises3D.md`` para formulación y acceptance.
"""
import math

import numpy as np
from numba import njit

from solidum.core.material import Material
from solidum.registry import MaterialRegistry


@njit
def _compute_j2_3d(strain, eps_p_old, alpha_old, sigma_y, H, K, G, C_e, yield_tol):
    """Return mapping J2 3D radial cerrado.

    Convención Voigt 6D del proyecto (ADR 0012):

    - Entrada ``strain``: ``[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`` con
      ``γ_ij = 2·ε_ij`` *engineering*.
    - Salida ``sigma``: ``[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz]`` (cortantes
      tensoriales, sin factor 2 — es lo que el ensamblador espera).
    - Estado interno ``eps_p``: ``[ε^p_xx, ε^p_yy, ε^p_zz, ε^p_xy, ε^p_yz,
      ε^p_xz]`` con cortantes **tensoriales** (sin factor 2). Coherente con
      el slot tensorial de ``VonMises2D``; permite computar el flujo plástico
      ``ε̇^p = γ̇·N`` y la norma de Frobenius sin conversiones repetidas.

    El predictor elástico se construye como ``s_trial + p·I`` desde el
    desviador trial — equivalente a ``C_e·(ε − ε^p_n)`` pero más numérica-
    mente estable y consistente con ``ε^p_n`` acumulada en descargas o
    reevaluaciones vía ``compute_gauss_state(U_final)``.
    """
    # Deformación volumétrica (traza)
    eps_v = strain[0] + strain[1] + strain[2]

    # Desviador de la deformación 6D en convención TENSORIAL
    # (cortantes divididos por 2 al entrar; saldrán tensoriales al ensamblar).
    e_dev = np.array([
        strain[0] - eps_v / 3.0,
        strain[1] - eps_v / 3.0,
        strain[2] - eps_v / 3.0,
        strain[3] / 2.0,
        strain[4] / 2.0,
        strain[5] / 2.0
    ])

    # Predictor elástico (desviador trial, restando ε^p_n tensorial)
    e_dev_trial = e_dev - eps_p_old
    s_trial = 2.0 * G * e_dev_trial

    # Norma de Frobenius tensorial: ||s||² = Σ_diag s_ii² + 2·Σ_off s_ij²
    norm_s_trial = math.sqrt(
        s_trial[0] ** 2 + s_trial[1] ** 2 + s_trial[2] ** 2
        + 2.0 * (s_trial[3] ** 2 + s_trial[4] ** 2 + s_trial[5] ** 2)
    )

    yield_stress = sigma_y + H * alpha_old
    f_trial = norm_s_trial - math.sqrt(2.0 / 3.0) * yield_stress

    p = K * eps_v

    if f_trial <= yield_tol:
        # Predictor elástico σ = s_trial + p·I (NO C_e·strain, que ignoraría
        # ε_p — bug detectado en VM2D plane strain, 2026-05-14, aplicado
        # desde el inicio en VM3D).
        sigma = np.array([
            s_trial[0] + p,
            s_trial[1] + p,
            s_trial[2] + p,
            s_trial[3],
            s_trial[4],
            s_trial[5]
        ])
        return sigma, C_e.copy(), eps_p_old.copy(), alpha_old

    # Corrector plástico — return mapping radial cerrado
    delta_gamma = f_trial / (2.0 * G + (2.0 / 3.0) * H)
    N = s_trial / norm_s_trial   # tensorial 6D, norma de Frobenius unitaria

    s_new = s_trial - 2.0 * G * delta_gamma * N
    eps_p_new = eps_p_old + delta_gamma * N
    alpha_new = alpha_old + math.sqrt(2.0 / 3.0) * delta_gamma

    sigma = np.array([
        s_new[0] + p,
        s_new[1] + p,
        s_new[2] + p,
        s_new[3],
        s_new[4],
        s_new[5]
    ])

    # Tangente algorítmica consistente (Simó-Hughes §3.3) en Voigt 6D del proyecto.
    # En esta convención (entrada engineering, salida tensorial off-diagonal):
    #   - v = [1, 1, 1, 0, 0, 0]      (operador traza en input engineering)
    #   - I_dev con 1/2 en los cortantes (mapea engineering γ a tensorial ε_ij)
    #   - N⊗N con N tensorial — el factor 2 implícito de Frobenius queda
    #     absorbido en cómo se aplica a γ_ij = 2·ε_ij.
    beta = (2.0 * G * delta_gamma) / norm_s_trial
    gamma_factor = 1.0 / (1.0 + H / (3.0 * G)) - beta

    v = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    I_dev = np.array([
        [ 2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0],
        [-1.0 / 3.0,  2.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0],
        [-1.0 / 3.0, -1.0 / 3.0,  2.0 / 3.0, 0.0, 0.0, 0.0],
        [ 0.0,        0.0,        0.0,       0.5, 0.0, 0.0],
        [ 0.0,        0.0,        0.0,       0.0, 0.5, 0.0],
        [ 0.0,        0.0,        0.0,       0.0, 0.0, 0.5]
    ])

    N_otimes_N = np.outer(N, N)

    C_alg = (
        K * np.outer(v, v)
        + 2.0 * G * (1.0 - beta) * I_dev
        - 2.0 * G * gamma_factor * N_otimes_N
    )

    return sigma, C_alg, eps_p_new, alpha_new


@MaterialRegistry.register
class VonMises3D(Material):
    """
    Modelo de plasticidad J2 (Von Mises) 3D con endurecimiento isótropo lineal.

    Sin variantes de hipótesis: en 3D todas las componentes de ``ε`` y ``σ``
    son activas. El return mapping es radial cerrado en forma única.

    Parameters
    ----------
    E : float
        Módulo de Young (>0).
    nu : float
        Coeficiente de Poisson, ν ∈ (-1, 0.5). El límite ν → 0.5 (incompresible)
        requiere formulación mixta u-p, no soportada.
    sigma_y : float
        Esfuerzo de fluencia inicial uniaxial (>0).
    H : float, optional
        Módulo de endurecimiento isótropo lineal (≥0). ``H = 0`` ⇒ plasticidad
        perfecta. Default 0.
    density : float, optional
        Densidad (ADR 0008). Opcional al construir; obligatoria si se ensambla
        peso propio o matriz de masa.

    Notes
    -----
    Convención Voigt 6D del proyecto (ADR 0012):
    ``ε = [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`` con ``γ_ij = 2·ε_ij``.
    El estado interno ``eps_p`` se almacena en 6 componentes
    ``[ε^p_xx, ε^p_yy, ε^p_zz, ε^p_xy, ε^p_yz, ε^p_xz]`` con cortantes
    **tensoriales** (sin factor 2), mismo convenio que el slot tensorial
    de ``VonMises2D``. La parte plástica es incompresible:
    ``tr(ε^p) = ε^p_xx + ε^p_yy + ε^p_zz = 0`` exacto en cualquier estado.
    """
    STRAIN_DIM = 6
    PRIMARY_STATE_VAR = 'alpha'

    def __init__(self, E: float, nu: float, sigma_y: float, H: float = 0.0,
                 density: float | None = None):
        if E <= 0.0:
            raise ValueError(f"VonMises3D: E debe ser > 0 (recibido {E}).")
        if sigma_y <= 0.0:
            raise ValueError(f"VonMises3D: sigma_y debe ser > 0 (recibido {sigma_y}).")
        if H < 0.0:
            raise ValueError(
                f"VonMises3D: H debe ser ≥ 0 (recibido {H}). "
                f"Ablandamiento (H<0) requiere regularización no implementada."
            )
        if not -1.0 < nu < 0.5:
            raise ValueError(
                f"VonMises3D: nu debe estar en (-1, 0.5) (recibido {nu}). "
                f"nu = 0.5 induce singularidad (K → ∞)."
            )
        if density is not None and density < 0.0:
            raise ValueError(
                f"VonMises3D: density={density} no puede ser negativa."
            )

        self.E = float(E)
        self.nu = float(nu)
        self.sigma_y = float(sigma_y)
        self.H = float(H)
        self.density = density

        # Módulos elásticos
        self.K = self.E / (3.0 * (1.0 - 2.0 * self.nu))
        self.G = self.E / (2.0 * (1.0 + self.nu))

        # Matriz constitutiva elástica 6×6 (idéntica a Elastic3D)
        coef = self.E / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        shear = (1.0 - 2.0 * self.nu) / 2.0
        a = 1.0 - self.nu
        b = self.nu
        self.C_e = coef * np.array([
            [a, b, b, 0.0, 0.0, 0.0],
            [b, a, b, 0.0, 0.0, 0.0],
            [b, b, a, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, shear, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, shear, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, shear],
        ])

    def admissibility_scale(self, state_vars=None) -> float:
        """Escala característica de admisibilidad (ADR 0006).

        ``√(2/3)·(σ_y + H·α)`` — norma desviadora característica en la frontera
        de fluencia J2. Idéntica a plane strain (la dimensión de ``Δγ`` y de
        ``N`` es la misma en ambos casos; lo que cambia respecto a plane stress
        es la naturaleza del flujo plástico).
        """
        alpha = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)
        R = self.sigma_y + self.H * alpha
        return math.sqrt(2.0 / 3.0) * R

    def compute_state(self, strain: np.ndarray, state_vars=None):
        """Devuelve ``(σ, C_tangent, new_state)``.

        Variables internas en ``state_vars``:
        - ``eps_p``: ndarray(6,) — [xx, yy, zz, xy, yz, xz] (cortantes tensoriales).
        - ``alpha``: float — deformación plástica acumulada equivalente.
        """
        eps_p_old = np.zeros(6) if state_vars is None else state_vars.get('eps_p', np.zeros(6))
        alpha_old = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)

        # Tolerancia precomputada fuera del kernel @njit (ADR 0006).
        yield_tol = self.admissibility_tol({'alpha': alpha_old})

        sigma, C_alg, eps_p_new, alpha_new = _compute_j2_3d(
            strain, eps_p_old, alpha_old,
            self.sigma_y, self.H, self.K, self.G, self.C_e, yield_tol
        )

        new_state = {'eps_p': eps_p_new, 'alpha': alpha_new}
        return sigma, C_alg, new_state
