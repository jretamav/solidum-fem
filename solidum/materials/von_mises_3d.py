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

from solidum.constants import ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL
from solidum.core.material import Material
from solidum.registry import MaterialRegistry


# Constantes del tangente J2 en Voigt 6D: v = operador traza sobre la
# entrada engineering; I_dev = proyector desviador con ½ en los cortantes.
_V6 = (1.0, 1.0, 1.0, 0.0, 0.0, 0.0)
_I_DEV6 = ((2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0),
           (-1.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0),
           (-1.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0, 0.0, 0.0, 0.0),
           (0.0, 0.0, 0.0, 0.5, 0.0, 0.0),
           (0.0, 0.0, 0.0, 0.0, 0.5, 0.0),
           (0.0, 0.0, 0.0, 0.0, 0.0, 0.5))


@njit(cache=True)
def _j2_3d_core(strain, eps_p_old, alpha_old, sigma_y, H, K, G, C_e, yield_tol,
                sigma, C_out, eps_p_new):
    """Return mapping J2 3D radial cerrado **sin asignar memoria**: escribe
    ``sigma`` (6), ``C_out`` (6×6) y ``eps_p_new`` (6) in situ y devuelve
    ``alpha_new``.

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

    Misma formulación y **mismo orden de operaciones** que la versión con
    arreglos que sustituye (ADR 0014): resultados bit a bit idénticos.
    """
    # Deformación volumétrica (traza)
    eps_v = strain[0] + strain[1] + strain[2]
    third = eps_v / 3.0

    # Desviador de la deformación 6D en convención TENSORIAL (cortantes
    # divididos por 2 al entrar) menos ε^p_n tensorial: predictor elástico.
    e0 = (strain[0] - third) - eps_p_old[0]
    e1 = (strain[1] - third) - eps_p_old[1]
    e2 = (strain[2] - third) - eps_p_old[2]
    e3 = (strain[3] / 2.0) - eps_p_old[3]
    e4 = (strain[4] / 2.0) - eps_p_old[4]
    e5 = (strain[5] / 2.0) - eps_p_old[5]
    twoG = 2.0 * G
    s0 = twoG * e0
    s1 = twoG * e1
    s2 = twoG * e2
    s3 = twoG * e3
    s4 = twoG * e4
    s5 = twoG * e5

    # Norma de Frobenius tensorial: ||s||² = Σ_diag s_ii² + 2·Σ_off s_ij²
    norm_s_trial = math.sqrt(
        s0 * s0 + s1 * s1 + s2 * s2
        + 2.0 * (s3 * s3 + s4 * s4 + s5 * s5)
    )

    yield_stress = sigma_y + H * alpha_old
    f_trial = norm_s_trial - math.sqrt(2.0 / 3.0) * yield_stress

    p = K * eps_v

    if f_trial <= yield_tol:
        # Predictor elástico σ = s_trial + p·I (NO C_e·strain, que ignoraría
        # ε_p — bug detectado en VM2D plane strain, 2026-05-14, aplicado
        # desde el inicio en VM3D).
        sigma[0] = s0 + p
        sigma[1] = s1 + p
        sigma[2] = s2 + p
        sigma[3] = s3
        sigma[4] = s4
        sigma[5] = s5
        for i in range(6):
            eps_p_new[i] = eps_p_old[i]
            for j in range(6):
                C_out[i, j] = C_e[i, j]
        return alpha_old

    # Corrector plástico — return mapping radial cerrado
    delta_gamma = f_trial / (2.0 * G + (2.0 / 3.0) * H)
    # N = s_trial / ||s_trial||: tensorial 6D, norma de Frobenius unitaria
    N0 = s0 / norm_s_trial
    N1 = s1 / norm_s_trial
    N2 = s2 / norm_s_trial
    N3 = s3 / norm_s_trial
    N4 = s4 / norm_s_trial
    N5 = s5 / norm_s_trial

    c = twoG * delta_gamma
    sigma[0] = (s0 - c * N0) + p
    sigma[1] = (s1 - c * N1) + p
    sigma[2] = (s2 - c * N2) + p
    sigma[3] = s3 - c * N3
    sigma[4] = s4 - c * N4
    sigma[5] = s5 - c * N5
    eps_p_new[0] = eps_p_old[0] + delta_gamma * N0
    eps_p_new[1] = eps_p_old[1] + delta_gamma * N1
    eps_p_new[2] = eps_p_old[2] + delta_gamma * N2
    eps_p_new[3] = eps_p_old[3] + delta_gamma * N3
    eps_p_new[4] = eps_p_old[4] + delta_gamma * N4
    eps_p_new[5] = eps_p_old[5] + delta_gamma * N5
    alpha_new = alpha_old + math.sqrt(2.0 / 3.0) * delta_gamma

    # Tangente algorítmica consistente (Simó-Hughes §3.3) en Voigt 6D del proyecto.
    # En esta convención (entrada engineering, salida tensorial off-diagonal):
    #   - v = [1, 1, 1, 0, 0, 0]      (operador traza en input engineering)
    #   - I_dev con 1/2 en los cortantes (mapea engineering γ a tensorial ε_ij)
    #   - N⊗N con N tensorial — el factor 2 implícito de Frobenius queda
    #     absorbido en cómo se aplica a γ_ij = 2·ε_ij.
    # C_alg = K·v⊗v + 2G(1−β)·I_dev − 2G·γ·N⊗N
    beta = c / norm_s_trial
    gamma_factor = 1.0 / (1.0 + H / (3.0 * G)) - beta
    a2 = twoG * (1.0 - beta)
    a3 = twoG * gamma_factor
    Nv = (N0, N1, N2, N3, N4, N5)
    for i in range(6):
        for j in range(6):
            C_out[i, j] = (K * (_V6[i] * _V6[j]) + a2 * _I_DEV6[i][j]) - a3 * (Nv[i] * Nv[j])
    return alpha_new


@njit(cache=True)
def _compute_j2_3d(strain, eps_p_old, alpha_old, sigma_y, H, K, G, C_e, yield_tol):
    """Return mapping J2 3D (camino por elemento): envoltorio con asignación
    de :func:`_j2_3d_core`. Devuelve ``(σ, C_alg, ε^p_new, α_new)``."""
    sigma = np.empty(6)
    C_alg = np.empty((6, 6))
    eps_p_new = np.empty(6)
    alpha_new = _j2_3d_core(strain, eps_p_old, alpha_old, sigma_y, H, K, G, C_e,
                            yield_tol, sigma, C_alg, eps_p_new)
    return sigma, C_alg, eps_p_new, alpha_new


@njit(cache=True)
def _j2_3d_batch(strain, S_old, S_new, params, C, sigma, C_out, flag):
    """Adaptador por lotes (ADR 0014). ``params = [σ_y, H, K, G, tol_abs,
    tol_rel]``; ``C`` = ``C_e``. Tolerancia como ``admissibility_tol``. Sin
    asignaciones por punto de Gauss."""
    sigma_y = params[0]
    H = params[1]
    alpha_old = S_old[6]
    R = sigma_y + H * alpha_old
    yield_tol = params[4] + params[5] * (math.sqrt(2.0 / 3.0) * R)
    S_new[6] = _j2_3d_core(
        strain, S_old[0:6], alpha_old, sigma_y, H, params[2], params[3], C, yield_tol,
        sigma, C_out, S_new[0:6],
    )
    return C_out


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
    STATE_SCHEMA = {'eps_p': (6,), 'alpha': ()}

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

    # Camino por lotes (ADR 0014)
    BATCH_KERNEL = _j2_3d_batch

    def batch_params(self) -> np.ndarray:
        return np.array([self.sigma_y, self.H, self.K, self.G,
                         ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL], dtype=np.float64)

    def batch_matrix(self) -> np.ndarray:
        return self.C_e
