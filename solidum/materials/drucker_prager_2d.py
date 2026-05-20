# fenix_fem/solidum/materials/drucker_prager_2d.py
"""``DruckerPrager2D`` — plasticidad friccional cohesivo-friccional 2D plane strain.

Modelo Drucker-Prager (cono circular suave de Mohr-Coulomb) con:

- **Plasticidad no asociada por defecto** (ángulo de dilatancia ψ ≠ ángulo de
  fricción φ); cuando ψ = φ es asociada y la tangente es simétrica.
- **Endurecimiento isótropo lineal en cohesión** ``c(α) = c_0 + H·α``.
- **Return mapping** con dos ramas cerradas: regular (cone surface) y al ápice
  del cono (vértice hidrostático). Detección automática del régimen.

Solo soporta ``hypothesis='plane_strain'`` en esta entrega. La proyección plane
stress de Drucker-Prager es notoriamente más delicada y queda *out-of-scope*.

Ver ``docs/specs/DruckerPrager2D.md``.
"""
import math

import numpy as np
from numba import njit

from solidum.core.material import Material
from solidum.registry import MaterialRegistry


@njit
def _compute_drucker_prager_plane_strain(strain, eps_p_old, alpha_old,
                                         eta_f, eta_g, k0, Hk,
                                         K, G, C_e, yield_tol):
    """Return mapping Drucker-Prager plane strain con detección regular/apex.

    Predictor elástico desviador-volumétrico (notación 4 componentes con
    ``ε_zz = 0`` plane strain). Criterio ``f = √J₂ + η_f·I₁ − k(α) ≤ 0``.

    En la rama plástica intenta primero retorno regular cerrado
    ``Δγ = f_trial/(G + 9K·η_f·η_g + Hk)``; si ``√J₂^{n+1} < 0`` activa retorno
    al ápice con ``Δγ_apex = (I₁_trial·η_f − k(α))/(9K·η_f·η_g + Hk)``.
    """
    # Predictor elástico — descomposición desviadora 3D plane strain con
    # `ε_zz = 0`. Trabajamos directamente sobre `ε - ε_p` (que tiene parte
    # volumétrica no nula por dilatancia plástica): e_dev = dev(ε - ε_p) y
    # NO `dev(ε) - ε_p`, que difiere en la parte volumétrica del flujo.
    # eps_p tensorial 4 componentes [xx, yy, zz, xy_tens].
    eps_minus_p_xx = strain[0] - eps_p_old[0]
    eps_minus_p_yy = strain[1] - eps_p_old[1]
    eps_minus_p_zz = 0.0 - eps_p_old[2]  # ε_zz = 0 plane strain
    eps_minus_p_xy_tens = strain[2] / 2.0 - eps_p_old[3]

    tr_eps_minus_p = eps_minus_p_xx + eps_minus_p_yy + eps_minus_p_zz
    third = tr_eps_minus_p / 3.0
    e_dev_trial = np.array([
        eps_minus_p_xx - third,
        eps_minus_p_yy - third,
        eps_minus_p_zz - third,
        eps_minus_p_xy_tens
    ])

    s_trial = 2.0 * G * e_dev_trial

    s_norm_sq_trial = (s_trial[0] * s_trial[0] + s_trial[1] * s_trial[1]
                       + s_trial[2] * s_trial[2] + 2.0 * s_trial[3] * s_trial[3])
    sqrt_J2_trial = math.sqrt(0.5 * s_norm_sq_trial)

    p_trial = K * tr_eps_minus_p
    I1_trial = 3.0 * p_trial

    k_curr = k0 + Hk * alpha_old
    f_trial = sqrt_J2_trial + eta_f * I1_trial - k_curr

    if f_trial <= yield_tol:
        # Paso elástico: σ = s_trial + p_trial·I
        sigma = np.array([
            s_trial[0] + p_trial,
            s_trial[1] + p_trial,
            s_trial[3]
        ])
        return sigma, C_e.copy(), eps_p_old.copy(), alpha_old, 0  # rama=0 elastica

    # Intento 1: return regular
    A = G + 9.0 * K * eta_f * eta_g + Hk
    delta_gamma = f_trial / A
    sqrt_J2_new = sqrt_J2_trial - G * delta_gamma

    if sqrt_J2_new >= 0.0:
        # Return regular — cone surface
        r = sqrt_J2_new / sqrt_J2_trial if sqrt_J2_trial > 0.0 else 0.0
        s_new = np.array([
            r * s_trial[0],
            r * s_trial[1],
            r * s_trial[2],
            r * s_trial[3]
        ])
        p_new = p_trial - 3.0 * K * eta_g * delta_gamma

        sigma = np.array([
            s_new[0] + p_new,
            s_new[1] + p_new,
            s_new[3]
        ])

        # Actualización de ε_p tensorial:
        # Δε_p = Δγ·(s_new/(2·√J2_new) + η_g·I)
        if sqrt_J2_new > 0.0:
            inv_2sqrtJ2 = 1.0 / (2.0 * sqrt_J2_new)
            deps_p_xx = delta_gamma * (s_new[0] * inv_2sqrtJ2 + eta_g)
            deps_p_yy = delta_gamma * (s_new[1] * inv_2sqrtJ2 + eta_g)
            deps_p_zz = delta_gamma * (s_new[2] * inv_2sqrtJ2 + eta_g)
            deps_p_xy = delta_gamma * s_new[3] * inv_2sqrtJ2
        else:
            # Caso límite: borde entre regular y apex
            deps_p_xx = delta_gamma * eta_g
            deps_p_yy = delta_gamma * eta_g
            deps_p_zz = delta_gamma * eta_g
            deps_p_xy = 0.0

        eps_p_new = np.array([
            eps_p_old[0] + deps_p_xx,
            eps_p_old[1] + deps_p_yy,
            eps_p_old[2] + deps_p_zz,
            eps_p_old[3] + deps_p_xy
        ])
        alpha_new = alpha_old + delta_gamma

        # Tangente algorítmica consistente — rama regular
        # n̂ en Voigt 3 (componente xy tensorial, NO engineering):
        # n̂_voigt[i] = s_trial_voigt[i] / (2·√J2_trial) con s_trial_voigt = [s_xx, s_yy, s_xy_tens]
        # v = [1, 1, 0] (gradient I1 en Voigt 3)
        beta = G * delta_gamma / sqrt_J2_trial   # término de "rotación" desviadora
        inv_2sJ2t = 1.0 / (2.0 * sqrt_J2_trial)
        n_hat = np.array([s_trial[0] * inv_2sJ2t, s_trial[1] * inv_2sJ2t, s_trial[3] * inv_2sJ2t])
        v = np.array([1.0, 1.0, 0.0])

        # Operadores de bloque:
        I_dev_voigt = np.array([
            [ 2.0 / 3.0, -1.0 / 3.0, 0.0],
            [-1.0 / 3.0,  2.0 / 3.0, 0.0],
            [ 0.0,        0.0,        0.5]
        ])

        # b_f = 2G·n̂ + 3K·η_f·v ; b_g = 2G·n̂ + 3K·η_g·v
        b_f = 2.0 * G * n_hat + 3.0 * K * eta_f * v
        b_g = 2.0 * G * n_hat + 3.0 * K * eta_g * v

        # C_alg = K·v⊗v + 2G·(1-β)·I_dev + 4G·β·n̂⊗n̂ − (1/A)·b_g⊗b_f
        # El término +4Gβ·n̂⊗n̂ recoge la dependencia de la dirección de flujo
        # con ε a través de s_trial (no solo el cambio de magnitud); sin él la
        # tangente cierra mal en componentes de cortante alineados con el flujo
        # (verificable por diferencia finita).
        C_alg = (K * np.outer(v, v)
                 + 2.0 * G * (1.0 - beta) * I_dev_voigt
                 + 4.0 * G * beta * np.outer(n_hat, n_hat)
                 - (1.0 / A) * np.outer(b_g, b_f))

        return sigma, C_alg, eps_p_new, alpha_new, 1   # rama=1 regular

    # Intento 2: return al ápice
    A_apex = 9.0 * K * eta_f * eta_g + Hk
    if A_apex <= 0.0:
        # Caso degenerado (η_g=0 y Hk=0): el ápice no admite descenso volumétrico.
        # No debería ocurrir si los parámetros son físicamente sensatos (η_g≥0, Hk≥0
        # y al menos uno >0 para que el ápice tenga rigidez).
        A_apex = 1.0e-30

    delta_gamma_apex = (I1_trial * eta_f - k_curr) / A_apex
    if delta_gamma_apex < 0.0:
        delta_gamma_apex = 0.0   # no debería ocurrir, salvaguarda numérica

    alpha_new = alpha_old + delta_gamma_apex
    k_new = k0 + Hk * alpha_new
    p_new = k_new / (3.0 * eta_f) if eta_f > 1.0e-30 else p_trial  # σ hidrostático en el ápice

    sigma = np.array([p_new, p_new, 0.0])

    # ε_p tras apex: el desviador trial completo se vuelve plástico, parte volumétrica
    # contribuye dilatación.
    # Δε_p_dev: el desviador plástico debe anular el desviador trial → Δe_p_dev = e_dev_trial.
    # Δε_p_vol_total = 3·η_g·Δγ.
    third_vol = eta_g * delta_gamma_apex
    eps_p_new = np.array([
        eps_p_old[0] + e_dev_trial[0] + third_vol,
        eps_p_old[1] + e_dev_trial[1] + third_vol,
        eps_p_old[2] + e_dev_trial[2] + third_vol,
        eps_p_old[3] + e_dev_trial[3]
    ])

    # Tangente al ápice: rigidez puramente volumétrica reducida.
    # C_alg_apex = (K·Hk / (9K·η_f·η_g + Hk)) · v⊗v
    v = np.array([1.0, 1.0, 0.0])
    if A_apex > 0.0:
        K_apex_eff = K * Hk / A_apex
    else:
        K_apex_eff = 0.0
    C_alg = K_apex_eff * np.outer(v, v)

    # Si la rigidez tangente colapsa por completo (apex con Hk=0 y η_g=0), añadir
    # un pequeño escalado de K para evitar matriz globalmente singular hasta que
    # el próximo paso salga del ápice. Pequeño = K·1e-6.
    if K_apex_eff < K * 1.0e-6:
        C_alg = (K * 1.0e-6) * np.outer(v, v)

    return sigma, C_alg, eps_p_new, alpha_new, 2   # rama=2 apex


def _calibrate_drucker_prager(c0: float, phi_rad: float, psi_rad: float, variant: str):
    """Devuelve ``(eta_f, k0, eta_g)`` desde parámetros físicos (c, φ, ψ)."""
    tan_phi = math.tan(phi_rad)
    tan_psi = math.tan(psi_rad)
    sin_phi = math.sin(phi_rad)
    sin_psi = math.sin(psi_rad)
    cos_phi = math.cos(phi_rad)

    if variant == 'plane_strain_matched':
        denom_f = math.sqrt(9.0 + 12.0 * tan_phi * tan_phi)
        eta_f = tan_phi / denom_f
        k0 = 3.0 * c0 / denom_f
        denom_g = math.sqrt(9.0 + 12.0 * tan_psi * tan_psi)
        eta_g = tan_psi / denom_g
    elif variant == 'outer_cone':
        sqrt3 = math.sqrt(3.0)
        eta_f = 2.0 * sin_phi / (sqrt3 * (3.0 - sin_phi))
        k0 = 6.0 * c0 * cos_phi / (sqrt3 * (3.0 - sin_phi))
        eta_g = 2.0 * sin_psi / (sqrt3 * (3.0 - sin_psi))
    elif variant == 'inner_cone':
        sqrt3 = math.sqrt(3.0)
        eta_f = 2.0 * sin_phi / (sqrt3 * (3.0 + sin_phi))
        k0 = 6.0 * c0 * cos_phi / (sqrt3 * (3.0 + sin_phi))
        eta_g = 2.0 * sin_psi / (sqrt3 * (3.0 + sin_psi))
    else:
        raise ValueError(
            f"variant={variant!r} no soportada. Usar "
            f"'plane_strain_matched', 'outer_cone' o 'inner_cone'."
        )
    return eta_f, k0, eta_g


@MaterialRegistry.register
class DruckerPrager2D(Material):
    """
    Modelo Drucker-Prager 2D plane strain con plasticidad no asociada y
    endurecimiento isótropo lineal en cohesión.

    Parameters
    ----------
    E, nu : float
        Parámetros elásticos isótropos.
    cohesion : float
        Cohesión inicial ``c_0`` (esfuerzo).
    phi_deg : float
        Ángulo de fricción interna en grados. ``0 ≤ φ < 90``. Para suelos típico
        20°–40°; hormigón 30°–37°.
    psi_deg : float, optional
        Ángulo de dilatancia en grados. Si ``None`` (default), se toma ``ψ = φ``
        (asociada). Para suelos no asociado típico ψ = 0 (sin dilatación) o ψ = φ/2.
        Debe cumplir ``0 ≤ ψ ≤ φ``.
    H : float, optional
        Módulo de endurecimiento isótropo lineal en cohesión (``≥ 0``). Default 0
        (perfectamente plástico).
    hypothesis : str, optional
        Solo ``'plane_strain'`` soportado en esta entrega.
    variant : str, optional
        Calibración Drucker-Prager con Mohr-Coulomb:
        ``'plane_strain_matched'`` (default, coincide con MC en plane strain),
        ``'outer_cone'`` (circunscribe MC), ``'inner_cone'`` (inscribe MC).
    density : float, optional
        Densidad (ADR 0008).

    Notes
    -----
    Cuando ``ψ ≠ φ`` (no asociado), la tangente algorítmica consistente es
    asimétrica → ``IS_SYMMETRIC = False`` declarativo. El despachador
    algebraico (ADR 0003) elige LU para el sistema global.
    """
    STRAIN_DIM = 3
    PRIMARY_STATE_VAR = 'alpha'
    # Default conservador a nivel de clase (caso no asociado). ``__init__``
    # sobrescribe la instancia con ``True`` cuando ``ψ = φ`` (asociada,
    # tangente algorítmica simétrica). ``domain_is_symmetric`` lee del
    # instance, así que el dispatcher elige Cholesky/LDLᵀ cuando puede.
    IS_SYMMETRIC = False

    def __init__(self, E: float, nu: float, cohesion: float, phi_deg: float,
                 psi_deg: float | None = None, H: float = 0.0,
                 hypothesis: str = 'plane_strain',
                 variant: str = 'plane_strain_matched',
                 density: float | None = None):
        if hypothesis != 'plane_strain':
            raise NotImplementedError(
                f"DruckerPrager2D: hypothesis={hypothesis!r} no soportado. "
                f"Esta versión solo implementa 'plane_strain'."
            )
        if E <= 0.0:
            raise ValueError(f"DruckerPrager2D: E debe ser > 0 (recibido {E}).")
        if cohesion <= 0.0:
            raise ValueError(f"DruckerPrager2D: cohesion debe ser > 0 (recibido {cohesion}).")
        if not 0.0 <= phi_deg < 90.0:
            raise ValueError(
                f"DruckerPrager2D: phi_deg debe estar en [0, 90) (recibido {phi_deg})."
            )
        if psi_deg is None:
            psi_deg = phi_deg
        if not 0.0 <= psi_deg <= phi_deg:
            raise ValueError(
                f"DruckerPrager2D: psi_deg debe cumplir 0 ≤ ψ ≤ φ "
                f"(recibido ψ={psi_deg}, φ={phi_deg})."
            )
        if H < 0.0:
            raise ValueError(
                f"DruckerPrager2D: H debe ser ≥ 0 (recibido {H}). "
                f"Ablandamiento (H<0) requiere regularización no implementada."
            )
        if not -1.0 < nu < 0.5:
            raise ValueError(
                f"DruckerPrager2D: nu debe estar en (-1, 0.5) (recibido {nu})."
            )

        self.E = E
        self.nu = nu
        self.cohesion_0 = cohesion
        self.phi_deg = phi_deg
        self.psi_deg = psi_deg
        self.H = H
        self.variant = variant
        self.density = density

        # Módulos elásticos
        self.K = E / (3.0 * (1.0 - 2.0 * nu))
        self.G = E / (2.0 * (1.0 + nu))

        # Matriz elástica plane strain Voigt 3
        coef = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
        self.C_e = coef * np.array([
            [1.0 - nu, nu,       0.0],
            [nu,       1.0 - nu, 0.0],
            [0.0,      0.0,      (1.0 - 2.0 * nu) / 2.0]
        ])

        phi_rad = math.radians(phi_deg)
        psi_rad = math.radians(psi_deg)
        self.eta_f, self.k0, self.eta_g = _calibrate_drucker_prager(
            cohesion, phi_rad, psi_rad, variant
        )
        # ¿asociada? — útil de exponer para inspección/tests
        self.associated = (psi_deg == phi_deg)
        # Override por instancia (ADR 0003): la tangente algorítmica del
        # return mapping regular es simétrica cuando ``ψ = φ`` (regla de
        # flujo asociada) — ``b_g ⊗ b_f`` colapsa a ``b_f ⊗ b_f``.
        # ``domain_is_symmetric`` lo lee del instance y el despachador
        # elige Cholesky/LDLᵀ para problemas asociados.
        self.IS_SYMMETRIC = self.associated

    def admissibility_scale(self, state_vars=None) -> float:
        """Cohesión efectiva corriente ``k(α) = k_0 + H·α`` (ADR 0006).

        Escala con unidades de esfuerzo (idéntica a la cantidad ``f`` del criterio).
        """
        alpha = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)
        return self.k0 + self.H * alpha

    def compute_state(self, strain: np.ndarray, state_vars=None):
        eps_p_old = np.zeros(4) if state_vars is None else state_vars.get('eps_p', np.zeros(4))
        alpha_old = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)

        yield_tol = self.admissibility_tol({'alpha': alpha_old})

        sigma, C_alg, eps_p_new, alpha_new, _branch = _compute_drucker_prager_plane_strain(
            strain, eps_p_old, alpha_old,
            self.eta_f, self.eta_g, self.k0, self.H,
            self.K, self.G, self.C_e, yield_tol,
        )

        new_state = {'eps_p': eps_p_new, 'alpha': alpha_new}
        return sigma, C_alg, new_state
