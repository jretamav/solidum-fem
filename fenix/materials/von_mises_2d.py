# fenix_fem/fenix/materials/von_mises_2d.py
"""``VonMises2D`` — plasticidad J2 con endurecimiento isótropo lineal en 2D.

Soporta dos hipótesis cinemáticas mutuamente excluyentes, seleccionadas al
construir el material y resueltas con kernels Numba especializados:

- ``hypothesis='plane_strain'``: descomposición volumétrica-desviadora 3D con
  ``ε_zz = 0`` impuesto; return mapping radial cerrado (Simó-Hughes §3.3).
- ``hypothesis='plane_stress'``: algoritmo proyectado con operador ``P``;
  return mapping con Newton local escalar sobre ``Δγ`` en base de autovectores
  de ``C_e^ps·P`` (Simó-Hughes §3.4.1, Box 3.1). Tangente algorítmica
  consistente cerrada.

Ver ``docs/specs/VonMises2D.md`` para formulación y acceptance.
"""
import math

import numpy as np
from numba import njit

from fenix.core.material import Material
from fenix.registry import MaterialRegistry


# Máximo de iteraciones del Newton local del return mapping plane stress.
# La función residual f̄(Δγ) es monótona decreciente con tangente cerrada,
# así que en la práctica converge en 3-6 iteraciones; cota holgada para
# pasos plásticos extremos (carga súbita muy por encima de fluencia).
_PLANE_STRESS_MAX_LOCAL_ITER = 25


@njit
def _compute_j2_plane_strain(strain, eps_p_old, alpha_old, sigma_y, H, K, G, C_e, yield_tol):
    """Return mapping J2 plane strain. Descomposición volumétrica-desviadora 3D
    con ``ε_zz = 0`` impuesto en la deformación total; ``ε^p_zz`` evoluciona
    libremente y queda registrada en ``eps_p`` (cuarta componente del estado).
    """
    # Deformación volumétrica (traza, ε_zz = 0 por plane strain)
    eps_v = strain[0] + strain[1]

    # Desviador de la deformación 3D extendida [xx, yy, zz, xy_tensorial]
    e_dev = np.array([
        strain[0] - eps_v / 3.0,
        strain[1] - eps_v / 3.0,
        -eps_v / 3.0,
        strain[2] / 2.0
    ])

    # Predictor elástico (desviador trial)
    e_dev_trial = e_dev - eps_p_old
    s_trial = 2.0 * G * e_dev_trial
    norm_s_trial = math.sqrt(s_trial[0]**2 + s_trial[1]**2 + s_trial[2]**2 + 2.0 * s_trial[3]**2)

    yield_stress = sigma_y + H * alpha_old
    f_trial = norm_s_trial - math.sqrt(2.0 / 3.0) * yield_stress

    p = K * eps_v

    if f_trial <= yield_tol:
        # Predictor elástico σ = s_trial + p·I (NO C_e·strain, que ignoraría ε_p).
        # Diferencia material en presencia de plasticidad acumulada (descarga o
        # reevaluación post-converged): C_e·strain devolvería un esfuerzo
        # incompatible con el estado interno; s_trial + p·I respeta ε_p_old.
        sigma = np.array([
            s_trial[0] + p,
            s_trial[1] + p,
            s_trial[3]
        ])
        return sigma, C_e.copy(), eps_p_old.copy(), alpha_old

    # Corrector plástico — return mapping radial cerrado
    delta_gamma = f_trial / (2.0 * G + (2.0 / 3.0) * H)
    N = s_trial / norm_s_trial

    s_new = s_trial - 2.0 * G * delta_gamma * N
    eps_p_new = eps_p_old + delta_gamma * N
    alpha_new = alpha_old + math.sqrt(2.0 / 3.0) * delta_gamma

    sigma = np.array([
        s_new[0] + p,
        s_new[1] + p,
        s_new[3]
    ])

    # Matriz tangente algorítmica consistente
    beta = (2.0 * G * delta_gamma) / norm_s_trial
    gamma_factor = 1.0 / (1.0 + H / (3.0 * G)) - beta

    v = np.array([1.0, 1.0, 0.0])
    I_dev = np.array([
        [ 2.0/3.0, -1.0/3.0, 0.0],
        [-1.0/3.0,  2.0/3.0, 0.0],
        [ 0.0,      0.0,     0.5]
    ])

    N_voigt = np.array([N[0], N[1], N[3]])
    N_otimes_N = np.array([
        [N_voigt[0]*N_voigt[0], N_voigt[0]*N_voigt[1], N_voigt[0]*N_voigt[2]],
        [N_voigt[1]*N_voigt[0], N_voigt[1]*N_voigt[1], N_voigt[1]*N_voigt[2]],
        [N_voigt[2]*N_voigt[0], N_voigt[2]*N_voigt[1], N_voigt[2]*N_voigt[2]]
    ])

    C_alg = K * np.outer(v, v) + 2.0 * G * (1.0 - beta) * I_dev - 2.0 * G * gamma_factor * N_otimes_N

    return sigma, C_alg, eps_p_new, alpha_new


@njit
def _compute_j2_plane_stress(strain, eps_p_old, alpha_old,
                             sigma_y, H, E, nu, G, C_e_ps,
                             yield_tol, max_local_iter):
    """Return mapping J2 plane stress proyectado (Simó-Hughes §3.4.1).

    Algoritmo:
    1. Predictor elástico ``σ_trial = C_e^ps (ε - ε^p_n)``.
    2. Función de fluencia proyectada
       ``f̄ = ½ σᵀ P σ − ⅓ (σ_y + H α)²``
       con ``P = (1/3) [[2,-1,0],[-1,2,0],[0,0,6]]``.
    3. Si ``f̄_trial ≤ tol``: paso elástico.
    4. Si no: Newton local escalar sobre ``Δγ`` en base autovectores de ``C_e^ps·P``
       (autovalores ``μ₁ = E/(3(1-ν))``, ``μ₂ = μ₃ = 2G``). La función residual
       ``f̄(Δγ)`` es monótona decreciente con tangente cerrada → convergencia
       cuadrática típica en 3-6 iteraciones.
    5. Actualización: ``σ = A⁻¹ σ_trial``, ``ε^p`` por regla de flujo
       ``ε̇^p = γ̇ P σ`` con cierre de incompresibilidad
       ``ε^p_zz = -(ε^p_xx + ε^p_yy)``, y deformación plástica acumulada
       equivalente ``α = α_n + Δγ · √(2·σPσ / 3)`` (norma Frobenius del
       flujo plástico, garantiza adimensionalidad de α bajo cambio de
       unidades — la regla ``α̇ = √(2/3)·γ̇`` no es válida en plane stress
       porque ``γ̇`` aquí tiene unidades ``[1/esfuerzo]``).
    6. Tangente consistente cerrada con corrección por dα/dΔγ no constante.
    """
    # Deformación plástica plana en Voigt engineering [xx, yy, γ_xy = 2·ε_xy_tens]
    eps_p_voigt = np.array([eps_p_old[0], eps_p_old[1], 2.0 * eps_p_old[3]])

    # Predictor elástico
    sigma_trial = C_e_ps @ (strain - eps_p_voigt)

    # P·σ_trial — operador P en notación Voigt mixta (γ_xy en deformación)
    P_sigma_trial = np.empty(3)
    P_sigma_trial[0] = (2.0 * sigma_trial[0] - sigma_trial[1]) / 3.0
    P_sigma_trial[1] = (-sigma_trial[0] + 2.0 * sigma_trial[1]) / 3.0
    P_sigma_trial[2] = 2.0 * sigma_trial[2]

    half_sPs_trial = 0.5 * (
        sigma_trial[0] * P_sigma_trial[0] +
        sigma_trial[1] * P_sigma_trial[1] +
        sigma_trial[2] * P_sigma_trial[2]
    )
    R_trial = sigma_y + H * alpha_old
    f_bar_trial = half_sPs_trial - (R_trial * R_trial) / 3.0

    if f_bar_trial <= yield_tol:
        return sigma_trial.copy(), C_e_ps.copy(), eps_p_old.copy(), alpha_old

    # Autovalores de C_e^ps · P
    mu1 = E / (3.0 * (1.0 - nu))
    mu_dev = 2.0 * G  # μ₂ = μ₃ = 2G

    # Proyecciones de σ_trial sobre autovectores ortonormales:
    # v₁ = [1,1,0]/√2 (hidrostático plano), v₂ = [1,-1,0]/√2 (desviador plano), v₃ = [0,0,1]
    inv_sqrt2 = 1.0 / math.sqrt(2.0)
    a1 = (sigma_trial[0] + sigma_trial[1]) * inv_sqrt2
    a2 = (sigma_trial[0] - sigma_trial[1]) * inv_sqrt2
    a3 = sigma_trial[2]

    # Newton local sobre Δγ
    delta_gamma = 0.0

    for _ in range(max_local_iter):
        denom1 = 1.0 + delta_gamma * mu1
        denom_dev = 1.0 + delta_gamma * mu_dev
        d1_sq = denom1 * denom1
        d_dev_sq = denom_dev * denom_dev

        # ½ σᵀ P σ en base autovectores: v_i^T P v_i = 1/3, 1, 2 (i=1,2,3)
        half_sPs = 0.5 * (
            (a1 * a1) / (3.0 * d1_sq) +
            (a2 * a2) / d_dev_sq +
            2.0 * (a3 * a3) / d_dev_sq
        )
        # w := √(2 σPσ / 3) — coeficiente de Δγ en la actualización de α
        # (norma Frobenius del flujo plástico ε̇^p = γ̇ P σ).
        sPs = 2.0 * half_sPs
        w = math.sqrt(2.0 * sPs / 3.0)

        alpha_curr = alpha_old + delta_gamma * w
        R_curr = sigma_y + H * alpha_curr

        f_bar = half_sPs - (R_curr * R_curr) / 3.0

        if abs(f_bar) < yield_tol:
            break

        # Derivadas en cadena
        # d(σPσ)/dΔγ — negativa porque las componentes σ_i = a_i / (1+Δγ μ_i) decrecen
        dsPs_dg = -2.0 * (
            (a1 * a1) * mu1 / (3.0 * denom1 * d1_sq)
            + (a2 * a2) * mu_dev / (denom_dev * d_dev_sq)
            + 2.0 * (a3 * a3) * mu_dev / (denom_dev * d_dev_sq)
        )
        dhalf_dg = 0.5 * dsPs_dg
        # dw/dΔγ = (1/(2w))·(2/3)·d(σPσ)/dΔγ = dsPs_dg / (3 w)
        if w > 1.0e-30:
            dw_dg = dsPs_dg / (3.0 * w)
        else:
            dw_dg = 0.0
        dalpha_dg = w + delta_gamma * dw_dg
        df_dg = dhalf_dg - (2.0 / 3.0) * R_curr * H * dalpha_dg

        if abs(df_dg) < 1e-300:
            break

        d_delta = -f_bar / df_dg
        delta_gamma += d_delta
        if delta_gamma < 0.0:
            delta_gamma = 0.0

    # σ_n+1 — construir A⁻¹ en base canónica desde la diagonal en autovectores
    denom1 = 1.0 + delta_gamma * mu1
    denom_dev = 1.0 + delta_gamma * mu_dev
    a_aux = 0.5 / denom1 + 0.5 / denom_dev
    b_aux = 0.5 / denom1 - 0.5 / denom_dev
    c_aux = 1.0 / denom_dev

    sigma_new = np.array([
        a_aux * sigma_trial[0] + b_aux * sigma_trial[1],
        b_aux * sigma_trial[0] + a_aux * sigma_trial[1],
        c_aux * sigma_trial[2]
    ])

    # Actualización de variables internas
    P_sigma_new = np.empty(3)
    P_sigma_new[0] = (2.0 * sigma_new[0] - sigma_new[1]) / 3.0
    P_sigma_new[1] = (-sigma_new[0] + 2.0 * sigma_new[1]) / 3.0
    P_sigma_new[2] = 2.0 * sigma_new[2]

    eps_p_xx_new = eps_p_old[0] + delta_gamma * P_sigma_new[0]
    eps_p_yy_new = eps_p_old[1] + delta_gamma * P_sigma_new[1]
    # P_sigma_new[2] está en notación engineering (γ); convertir a tensorial
    eps_p_xy_new_tens = eps_p_old[3] + delta_gamma * P_sigma_new[2] / 2.0
    # Incompresibilidad plástica → ε^p_zz determinada por tr(ε^p) = 0
    eps_p_zz_new = -(eps_p_xx_new + eps_p_yy_new)

    eps_p_new = np.array([eps_p_xx_new, eps_p_yy_new, eps_p_zz_new, eps_p_xy_new_tens])

    # α final con σPσ convergido
    sPs_final = (a1 * a1) / (3.0 * denom1 * denom1) \
              + (a2 * a2) / (denom_dev * denom_dev) \
              + 2.0 * (a3 * a3) / (denom_dev * denom_dev)
    w_final = math.sqrt(2.0 * sPs_final / 3.0)
    alpha_new = alpha_old + delta_gamma * w_final

    # Tangente algorítmica consistente.
    # dσ_new = D dε - (D P σ) dΔγ          con D = A⁻¹ C_e^ps
    # ḟ̄ = (σ P D) dε - q dΔγ - (2/3) R H dα = 0
    # dα/dε    = (2 Δγ / (3 w)) σ P D
    # dα/dΔγ   = w − (2 Δγ / (3 w)) q
    # Despejando dΔγ y sustituyendo:
    # C_alg = D − (D P σ) ⊗ (σ P D) / [q + m w / (1 − m n)]
    # con m = (2/3) R H,  n = 2 Δγ / (3 w),  q = σ P D P σ
    A_inv = np.array([
        [a_aux, b_aux, 0.0],
        [b_aux, a_aux, 0.0],
        [0.0,   0.0,   c_aux]
    ])
    D = A_inv @ C_e_ps

    D_P_sigma = D @ P_sigma_new
    sigma_P_D = P_sigma_new @ D   # P simétrica ⇒ σᵀP = (Pσ)ᵀ
    q = P_sigma_new @ (D @ P_sigma_new)

    R_new = sigma_y + H * alpha_new
    m = (2.0 / 3.0) * R_new * H
    if w_final > 1.0e-30:
        n = 2.0 * delta_gamma / (3.0 * w_final)
    else:
        n = 0.0

    one_minus_mn = 1.0 - m * n
    # Forma agrupada equivalente, numéricamente estable cuando m·n → 1:
    # denom = q + m·w / (1 - m·n)
    if abs(one_minus_mn) < 1.0e-30 or H == 0.0:
        denom_beta = q  # H=0 → no hay corrección por endurecimiento
    else:
        denom_beta = q + m * w_final / one_minus_mn

    if abs(denom_beta) < 1.0e-30:
        C_alg = D.copy()
    else:
        C_alg = D - np.outer(D_P_sigma, sigma_P_D) / denom_beta

    return sigma_new, C_alg, eps_p_new, alpha_new


@MaterialRegistry.register
class VonMises2D(Material):
    """
    Modelo de plasticidad J2 (Von Mises) con endurecimiento isótropo lineal.

    Dos hipótesis cinemáticas mutuamente excluyentes, fijadas al construir y
    resueltas con kernels Numba especializados:

    - ``hypothesis='plane_strain'`` (default): ``ε_zz = 0`` impuesto.
    - ``hypothesis='plane_stress'``: ``σ_zz = 0`` impuesto mediante algoritmo
      proyectado (Simó-Hughes §3.4.1).

    Parameters
    ----------
    E : float
        Módulo de Young.
    nu : float
        Coeficiente de Poisson. Debe cumplir |ν| < 1/2; ν cercano a 0.5 induce
        locking volumétrico en plane strain con elementos de bajo orden.
    sigma_y : float
        Tensión de fluencia inicial (uniaxial).
    H : float, optional
        Módulo de endurecimiento isótropo lineal. ``H >= 0``; ``H = 0`` ⇒
        plasticidad perfecta. Default 0.
    hypothesis : str, optional
        ``'plane_strain'`` (default) o ``'plane_stress'``.
    density : float, optional
        Densidad (ADR 0008). Opcional al construir; obligatoria si se ensambla
        peso propio o matriz de masa.

    Notes
    -----
    Convención Voigt 2D: ``ε = [ε_xx, ε_yy, γ_xy]`` con ``γ_xy = 2·ε_xy``.
    El estado interno ``eps_p`` se almacena en 4 componentes
    ``[ε^p_xx, ε^p_yy, ε^p_zz, ε^p_xy_tensorial]`` (xy sin factor 2). En
    plane stress, ``ε^p_zz`` se determina por incompresibilidad plástica
    ``tr(ε^p) = 0``; en plane strain evoluciona libremente con el return
    mapping desviador.
    """
    STRAIN_DIM = 3
    PRIMARY_STATE_VAR = 'alpha'

    def __init__(self, E: float, nu: float, sigma_y: float, H: float = 0.0,
                 hypothesis: str = 'plane_strain', density: float | None = None):
        if hypothesis not in ('plane_strain', 'plane_stress'):
            raise ValueError(
                f"VonMises2D: hypothesis='{hypothesis}' no soportado. "
                f"Usar 'plane_strain' o 'plane_stress'."
            )
        if E <= 0.0:
            raise ValueError(f"VonMises2D: E debe ser > 0 (recibido {E}).")
        if sigma_y <= 0.0:
            raise ValueError(f"VonMises2D: sigma_y debe ser > 0 (recibido {sigma_y}).")
        if H < 0.0:
            raise ValueError(
                f"VonMises2D: H debe ser ≥ 0 (recibido {H}). "
                f"Ablandamiento (H<0) requiere regularización no implementada."
            )
        if not -1.0 < nu < 0.5:
            raise ValueError(
                f"VonMises2D: nu debe estar en (-1, 0.5) (recibido {nu}). "
                f"nu = 0.5 induce singularidad."
            )

        self.E = E
        self.nu = nu
        self.sigma_y = sigma_y
        self.H = H
        self.hypothesis = hypothesis
        self.density = density

        # Módulos elásticos (K solo se usa en plane strain)
        self.K = E / (3.0 * (1.0 - 2.0 * nu))
        self.G = E / (2.0 * (1.0 + nu))

        if hypothesis == 'plane_strain':
            coef = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
            self.C_e = coef * np.array([
                [1.0 - nu, nu,       0.0],
                [nu,       1.0 - nu, 0.0],
                [0.0,      0.0,      (1.0 - 2.0 * nu) / 2.0]
            ])
        else:  # 'plane_stress'
            coef_ps = E / (1.0 - nu * nu)
            self.C_e = coef_ps * np.array([
                [1.0,  nu,   0.0],
                [nu,   1.0,  0.0],
                [0.0,  0.0,  (1.0 - nu) / 2.0]
            ])

    def admissibility_scale(self, state_vars=None) -> float:
        """Escala característica del check de admisibilidad (ADR 0006).

        - Plane strain: ``√(2/3)·(σ_y + H·α)`` — norma desviadora en la frontera J2.
        - Plane stress: ``(σ_y + H·α)²/3`` — escala de la función de fluencia
          proyectada ``f̄ = ½σᵀPσ − ⅓R²`` (los dos sumandos son del mismo orden
          en la frontera).
        """
        alpha = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)
        R = self.sigma_y + self.H * alpha
        if self.hypothesis == 'plane_strain':
            return math.sqrt(2.0 / 3.0) * R
        return (R * R) / 3.0

    def compute_state(self, strain: np.ndarray, state_vars=None):
        """Devuelve ``(σ, C_tangent, new_state)`` para la hipótesis activa.

        Variables internas en ``state_vars``:
        - ``eps_p``: ndarray(4,) — [xx, yy, zz, xy_tensorial].
        - ``alpha``: float — deformación plástica acumulada equivalente.
        """
        eps_p_old = np.zeros(4) if state_vars is None else state_vars.get('eps_p', np.zeros(4))
        alpha_old = 0.0 if state_vars is None else state_vars.get('alpha', 0.0)

        # Tolerancia precomputada fuera del kernel @njit (ADR 0006).
        yield_tol = self.admissibility_tol({'alpha': alpha_old})

        if self.hypothesis == 'plane_strain':
            sigma, C_alg, eps_p_new, alpha_new = _compute_j2_plane_strain(
                strain, eps_p_old, alpha_old,
                self.sigma_y, self.H, self.K, self.G, self.C_e, yield_tol
            )
        else:  # 'plane_stress'
            sigma, C_alg, eps_p_new, alpha_new = _compute_j2_plane_stress(
                strain, eps_p_old, alpha_old,
                self.sigma_y, self.H, self.E, self.nu, self.G, self.C_e,
                yield_tol, _PLANE_STRESS_MAX_LOCAL_ITER
            )

        new_state = {'eps_p': eps_p_new, 'alpha': alpha_new}
        return sigma, C_alg, new_state
