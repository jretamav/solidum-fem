# solidum_fem/solidum/materials/drucker_prager_3d.py
"""``DruckerPrager3D`` — plasticidad friccional cohesivo-friccional 3D.

Modelo Drucker-Prager (cono circular suave de Mohr-Coulomb) en Voigt 6D del
proyecto (ADR 0012). Mismas características que ``DruckerPrager2D``:

- **Plasticidad no asociada por defecto** (ángulo de dilatancia ψ ≠ ángulo de
  fricción φ); cuando ψ = φ es asociada y la tangente es simétrica.
- **Endurecimiento isótropo lineal en cohesión** ``c(α) = c_0 + H·α``.
- **Return mapping** con dos ramas cerradas: regular (cone surface) y al ápice
  (vértice hidrostático). Detección automática del régimen.

Diferencias respecto a 2D:

- Sin variantes de hipótesis (en 3D todas las componentes son activas).
- **No incluye la calibración ``plane_strain_matched``**: ese match exacto
  DP↔MC es 2D-only (en 3D MC depende del ángulo de Lode y no admite
  coincidencia circular global). Solo ``outer_cone`` (default, circunscribe
  MC en compresión) e ``inner_cone`` (inscribe MC en extensión).

Ver ``docs/specs/DruckerPrager3D.md`` para formulación y acceptance.
"""
import math

import numpy as np
from numba import njit

from solidum.core.material import Material
from solidum.materials.drucker_prager_2d import _calibrate_drucker_prager
from solidum.registry import MaterialRegistry


@njit
def _compute_drucker_prager_3d(strain, eps_p_old, alpha_old,
                               eta_f, eta_g, k0, Hk,
                               K, G, C_e, yield_tol):
    """Return mapping Drucker-Prager 3D con detección regular/apex.

    Convención Voigt 6D del proyecto (ADR 0012):

    - Entrada ``strain``: ``[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`` engineering
      (``γ_ij = 2·ε_ij``).
    - Salida ``sigma``: ``[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz]`` con cortantes
      tensoriales.
    - Estado ``eps_p`` (6,) con cortantes **tensoriales** sin factor 2.
      A diferencia de J2, ``tr(ε^p) = 3·η_g·α`` (no incompresible) — invariante
      cinemático verificable.

    Predictor elástico: ``s_trial = 2G·dev(ε - ε^p_n)``,
    ``p_trial = K·tr(ε - ε^p_n)``. La sustracción de ``ε^p_n`` antes de tomar
    el desviador es **necesaria** porque ``ε^p_n`` tiene parte volumétrica
    (dilatancia plástica). Criterio ``f = √J₂ + η_f·I₁ − k(α) ≤ 0``.

    Rama regular cerrada: ``Δγ = f_trial/(G + 9K·η_f·η_g + H_k)``; valida si
    ``√J₂^{n+1} = √J₂_trial − G·Δγ ≥ 0``. Si no, rama de ápice cerrada
    ``Δγ_apex = (I₁_trial·η_f − k(α))/(9K·η_f·η_g + H_k)`` con σ puramente
    hidrostático y tangente volumétrica reducida.
    """
    # Conversión engineering → tensorial (cortantes / 2) y sustracción de ε^p_n
    em0 = strain[0] - eps_p_old[0]
    em1 = strain[1] - eps_p_old[1]
    em2 = strain[2] - eps_p_old[2]
    em3 = strain[3] / 2.0 - eps_p_old[3]
    em4 = strain[4] / 2.0 - eps_p_old[4]
    em5 = strain[5] / 2.0 - eps_p_old[5]

    tr_eps_minus_p = em0 + em1 + em2
    third = tr_eps_minus_p / 3.0
    e_dev_trial = np.array([
        em0 - third,
        em1 - third,
        em2 - third,
        em3,
        em4,
        em5
    ])

    s_trial = 2.0 * G * e_dev_trial

    # Norma Frobenius tensorial: ||s||² = Σ_diag + 2·Σ_off
    s_norm_sq_trial = (
        s_trial[0] ** 2 + s_trial[1] ** 2 + s_trial[2] ** 2
        + 2.0 * (s_trial[3] ** 2 + s_trial[4] ** 2 + s_trial[5] ** 2)
    )
    sqrt_J2_trial = math.sqrt(0.5 * s_norm_sq_trial)

    p_trial = K * tr_eps_minus_p
    I1_trial = 3.0 * p_trial

    k_curr = k0 + Hk * alpha_old
    f_trial = sqrt_J2_trial + eta_f * I1_trial - k_curr

    if f_trial <= yield_tol:
        sigma = np.array([
            s_trial[0] + p_trial,
            s_trial[1] + p_trial,
            s_trial[2] + p_trial,
            s_trial[3],
            s_trial[4],
            s_trial[5]
        ])
        return sigma, C_e.copy(), eps_p_old.copy(), alpha_old, 0  # rama=0 elástica

    # ------------------------------------------------------------------
    # Intento 1: return regular (cone surface)
    # ------------------------------------------------------------------
    A = G + 9.0 * K * eta_f * eta_g + Hk
    delta_gamma = f_trial / A
    sqrt_J2_new = sqrt_J2_trial - G * delta_gamma

    if sqrt_J2_new >= 0.0:
        r = sqrt_J2_new / sqrt_J2_trial if sqrt_J2_trial > 0.0 else 0.0
        s_new = np.array([
            r * s_trial[0],
            r * s_trial[1],
            r * s_trial[2],
            r * s_trial[3],
            r * s_trial[4],
            r * s_trial[5]
        ])
        p_new = p_trial - 3.0 * K * eta_g * delta_gamma

        sigma = np.array([
            s_new[0] + p_new,
            s_new[1] + p_new,
            s_new[2] + p_new,
            s_new[3],
            s_new[4],
            s_new[5]
        ])

        # Δε_p = Δγ·(s_new/(2·√J2_new) + η_g·I_vol) en convención tensorial
        if sqrt_J2_new > 0.0:
            inv_2sqrtJ2 = 1.0 / (2.0 * sqrt_J2_new)
            deps_p_xx = delta_gamma * (s_new[0] * inv_2sqrtJ2 + eta_g)
            deps_p_yy = delta_gamma * (s_new[1] * inv_2sqrtJ2 + eta_g)
            deps_p_zz = delta_gamma * (s_new[2] * inv_2sqrtJ2 + eta_g)
            deps_p_xy = delta_gamma * s_new[3] * inv_2sqrtJ2
            deps_p_yz = delta_gamma * s_new[4] * inv_2sqrtJ2
            deps_p_xz = delta_gamma * s_new[5] * inv_2sqrtJ2
        else:
            # Caso límite: borde entre regular y apex
            deps_p_xx = delta_gamma * eta_g
            deps_p_yy = delta_gamma * eta_g
            deps_p_zz = delta_gamma * eta_g
            deps_p_xy = 0.0
            deps_p_yz = 0.0
            deps_p_xz = 0.0

        eps_p_new = np.array([
            eps_p_old[0] + deps_p_xx,
            eps_p_old[1] + deps_p_yy,
            eps_p_old[2] + deps_p_zz,
            eps_p_old[3] + deps_p_xy,
            eps_p_old[4] + deps_p_yz,
            eps_p_old[5] + deps_p_xz
        ])
        alpha_new = alpha_old + delta_gamma

        # Tangente algorítmica consistente — rama regular en Voigt 6D del proyecto.
        # Estructura (de Souza Neto §8.3.4 extendida a 6D):
        #     C_alg = K·v⊗v + 2G(1-β)·I_dev + 4G·β·n̂⊗n̂ − (1/A)·b_g⊗b_f
        # con I_dev 6×6 que mapea engineering input a tensorial deviator (1/2
        # en los cortantes), y v = [1,1,1,0,0,0]. n̂ y b_f, b_g en Voigt 6D
        # tensorial. Asimétrica si η_f ≠ η_g.
        beta = G * delta_gamma / sqrt_J2_trial
        inv_2sJ2t = 1.0 / (2.0 * sqrt_J2_trial)
        n_hat = np.array([
            s_trial[0] * inv_2sJ2t,
            s_trial[1] * inv_2sJ2t,
            s_trial[2] * inv_2sJ2t,
            s_trial[3] * inv_2sJ2t,
            s_trial[4] * inv_2sJ2t,
            s_trial[5] * inv_2sJ2t
        ])
        v = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])

        I_dev = np.array([
            [ 2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0],
            [-1.0 / 3.0,  2.0 / 3.0, -1.0 / 3.0, 0.0, 0.0, 0.0],
            [-1.0 / 3.0, -1.0 / 3.0,  2.0 / 3.0, 0.0, 0.0, 0.0],
            [ 0.0,        0.0,        0.0,       0.5, 0.0, 0.0],
            [ 0.0,        0.0,        0.0,       0.0, 0.5, 0.0],
            [ 0.0,        0.0,        0.0,       0.0, 0.0, 0.5]
        ])

        b_f = 2.0 * G * n_hat + 3.0 * K * eta_f * v
        b_g = 2.0 * G * n_hat + 3.0 * K * eta_g * v

        C_alg = (
            K * np.outer(v, v)
            + 2.0 * G * (1.0 - beta) * I_dev
            + 4.0 * G * beta * np.outer(n_hat, n_hat)
            - (1.0 / A) * np.outer(b_g, b_f)
        )

        return sigma, C_alg, eps_p_new, alpha_new, 1  # rama=1 regular

    # ------------------------------------------------------------------
    # Intento 2: return al ápice
    # ------------------------------------------------------------------
    A_apex = 9.0 * K * eta_f * eta_g + Hk
    if A_apex <= 0.0:
        # Degenerado (η_g=0 ∧ Hk=0): salvaguarda numérica.
        A_apex = 1.0e-30

    delta_gamma_apex = (I1_trial * eta_f - k_curr) / A_apex
    if delta_gamma_apex < 0.0:
        delta_gamma_apex = 0.0   # no debería ocurrir; salvaguarda

    alpha_new = alpha_old + delta_gamma_apex
    k_new = k0 + Hk * alpha_new
    p_new = k_new / (3.0 * eta_f) if eta_f > 1.0e-30 else p_trial

    sigma = np.array([p_new, p_new, p_new, 0.0, 0.0, 0.0])

    # ε_p tras apex: el desviador trial completo se vuelve plástico, parte
    # volumétrica añade dilatación η_g·Δγ por componente diagonal.
    third_vol = eta_g * delta_gamma_apex
    eps_p_new = np.array([
        eps_p_old[0] + e_dev_trial[0] + third_vol,
        eps_p_old[1] + e_dev_trial[1] + third_vol,
        eps_p_old[2] + e_dev_trial[2] + third_vol,
        eps_p_old[3] + e_dev_trial[3],
        eps_p_old[4] + e_dev_trial[4],
        eps_p_old[5] + e_dev_trial[5]
    ])

    # Tangente al ápice: rigidez puramente volumétrica reducida.
    v = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    if A_apex > 0.0:
        K_apex_eff = K * Hk / A_apex
    else:
        K_apex_eff = 0.0

    # Si la rigidez tangente colapsa por completo (apex con Hk=0 y η_g=0),
    # añadir un escalado mínimo de K para evitar matriz globalmente singular
    # hasta que el siguiente paso salga del ápice (mismo patrón que DP2D).
    if K_apex_eff < K * 1.0e-6:
        K_apex_eff = K * 1.0e-6
    C_alg = K_apex_eff * np.outer(v, v)

    return sigma, C_alg, eps_p_new, alpha_new, 2   # rama=2 apex


@MaterialRegistry.register
class DruckerPrager3D(Material):
    """
    Modelo Drucker-Prager 3D con plasticidad no asociada y endurecimiento
    isótropo lineal en cohesión.

    Parameters
    ----------
    E, nu : float
        Parámetros elásticos isótropos. ``ν ∈ (-1, 0.5)`` estricto.
    cohesion : float
        Cohesión inicial ``c_0`` (esfuerzo, >0).
    phi_deg : float
        Ángulo de fricción interna en grados. ``0 ≤ φ < 90``.
    psi_deg : float, optional
        Ángulo de dilatancia en grados. Si ``None`` (default), se toma
        ``ψ = φ`` (asociada). Para suelos no asociado típico ψ = 0 o ψ = φ/2.
        Debe cumplir ``0 ≤ ψ ≤ φ``.
    H : float, optional
        Módulo de endurecimiento isótropo lineal en cohesión (``≥ 0``). Default 0
        (perfectamente plástico).
    variant : str, optional
        Calibración Drucker-Prager con Mohr-Coulomb (3D):
        ``'outer_cone'`` (default, circunscribe MC en el meridiano de compresión)
        o ``'inner_cone'`` (inscribe MC en el meridiano de extensión).
        **No se acepta ``'plane_strain_matched'``** en 3D — esa calibración es
        2D-only por construcción (MC en 3D depende del ángulo de Lode).
    density : float, optional
        Densidad (ADR 0008). Opcional al construir; obligatoria si se ensambla
        peso propio o masa.

    Notes
    -----
    Convención Voigt 6D del proyecto (ADR 0012):
    ``ε = [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz]`` con ``γ_ij = 2·ε_ij``.
    El estado interno ``eps_p`` se almacena en 6 componentes con cortantes
    **tensoriales** (sin factor 2). A diferencia de J2, la parte plástica
    **no es incompresible**: ``tr(ε^p) = 3·η_g·α`` (dilatancia si ``η_g > 0``).

    Cuando ``ψ ≠ φ`` (no asociado), la tangente algorítmica consistente es
    asimétrica → ``IS_SYMMETRIC = False`` declarativo. El despachador
    algebraico (ADR 0003) elige LU para el sistema global. Cuando ``ψ = φ``
    el override por instancia activa ``IS_SYMMETRIC = True`` y el
    despachador puede usar Cholesky/LDLᵀ.
    """
    STRAIN_DIM = 6
    PRIMARY_STATE_VAR = 'alpha'
    # Default conservador a nivel de clase (caso no asociado). ``__init__``
    # sobrescribe la instancia con ``True`` cuando ``ψ = φ`` (asociada).
    IS_SYMMETRIC = False

    # Calibraciones soportadas en 3D (omitimos plane_strain_matched: 2D-only).
    _SUPPORTED_VARIANTS = ('outer_cone', 'inner_cone')

    def __init__(self, E: float, nu: float, cohesion: float, phi_deg: float,
                 psi_deg: float | None = None, H: float = 0.0,
                 variant: str = 'outer_cone',
                 density: float | None = None):
        if E <= 0.0:
            raise ValueError(f"DruckerPrager3D: E debe ser > 0 (recibido {E}).")
        if cohesion <= 0.0:
            raise ValueError(
                f"DruckerPrager3D: cohesion debe ser > 0 (recibido {cohesion})."
            )
        if not 0.0 <= phi_deg < 90.0:
            raise ValueError(
                f"DruckerPrager3D: phi_deg debe estar en [0, 90) (recibido {phi_deg})."
            )
        if psi_deg is None:
            psi_deg = phi_deg
        if not 0.0 <= psi_deg <= phi_deg:
            raise ValueError(
                f"DruckerPrager3D: psi_deg debe cumplir 0 ≤ ψ ≤ φ "
                f"(recibido ψ={psi_deg}, φ={phi_deg})."
            )
        if H < 0.0:
            raise ValueError(
                f"DruckerPrager3D: H debe ser ≥ 0 (recibido {H}). "
                f"Ablandamiento (H<0) requiere regularización no implementada."
            )
        if not -1.0 < nu < 0.5:
            raise ValueError(
                f"DruckerPrager3D: nu debe estar en (-1, 0.5) (recibido {nu})."
            )
        if variant not in self._SUPPORTED_VARIANTS:
            raise ValueError(
                f"DruckerPrager3D: variant={variant!r} no soportada en 3D. "
                f"Usar 'outer_cone' (default) o 'inner_cone'. "
                f"La calibración 'plane_strain_matched' es 2D-only por "
                f"construcción (MC en 3D depende del ángulo de Lode)."
            )
        if density is not None and density < 0.0:
            raise ValueError(
                f"DruckerPrager3D: density={density} no puede ser negativa."
            )

        self.E = float(E)
        self.nu = float(nu)
        self.cohesion_0 = float(cohesion)
        self.phi_deg = float(phi_deg)
        self.psi_deg = float(psi_deg)
        self.H = float(H)
        self.variant = variant
        self.density = density

        # Módulos elásticos
        self.K = self.E / (3.0 * (1.0 - 2.0 * self.nu))
        self.G = self.E / (2.0 * (1.0 + self.nu))

        # Matriz constitutiva elástica 6×6 (idéntica a Elastic3D / VonMises3D)
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

        phi_rad = math.radians(phi_deg)
        psi_rad = math.radians(psi_deg)
        # Reusa la calibración de DP2D (mismas fórmulas para outer/inner;
        # plane_strain_matched ya filtrado por el validador de variant).
        self.eta_f, self.k0, self.eta_g = _calibrate_drucker_prager(
            cohesion, phi_rad, psi_rad, variant
        )
        self.associated = (psi_deg == phi_deg)
        # Override por instancia (ADR 0003): tangente algorítmica simétrica
        # cuando ψ = φ (regla de flujo asociada → b_g ⊗ b_f colapsa a b_f ⊗ b_f).
        self.IS_SYMMETRIC = self.associated

    def admissibility_scale(self, state_vars=None) -> float:
        """Cohesión efectiva corriente ``k(α) = k_0 + H·α`` (ADR 0006).

        Escala con unidades de esfuerzo (idéntica a la cantidad ``f`` del criterio).
        """
        alpha = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)
        return self.k0 + self.H * alpha

    def compute_state(self, strain: np.ndarray, state_vars=None):
        """Devuelve ``(σ, C_tangent, new_state)``.

        Variables internas en ``state_vars``:
        - ``eps_p``: ndarray(6,) — [xx, yy, zz, xy, yz, xz] (cortantes tensoriales).
        - ``alpha``: float — multiplicador plástico acumulado.
        """
        eps_p_old = np.zeros(6) if state_vars is None else state_vars.get('eps_p', np.zeros(6))
        alpha_old = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)

        yield_tol = self.admissibility_tol({'alpha': alpha_old})

        sigma, C_alg, eps_p_new, alpha_new, _branch = _compute_drucker_prager_3d(
            strain, eps_p_old, alpha_old,
            self.eta_f, self.eta_g, self.k0, self.H,
            self.K, self.G, self.C_e, yield_tol,
        )

        new_state = {'eps_p': eps_p_new, 'alpha': alpha_new}
        return sigma, C_alg, new_state
