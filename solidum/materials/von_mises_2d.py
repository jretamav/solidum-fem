# solidum_fem/solidum/materials/von_mises_2d.py
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

from solidum.constants import (
    ADMISSIBILITY_TOL_ABS,
    ADMISSIBILITY_TOL_REL,
    J2_DENOM_FLOOR,
    J2_PLANE_STRESS_MAX_LOCAL_ITER,
)
from solidum.logging import get_logger
from solidum.core.material import Material
from solidum.materials._plane_strain import sigma_zz_plane_strain
from solidum.registry import MaterialRegistry

_log = get_logger("materials")


# Re-export local del límite de iteraciones (centralizado en solidum.constants
# tras auditoría H-1.9). Conservamos el alias de módulo porque los kernels
# Numba están compilados con la referencia local.
_PLANE_STRESS_MAX_LOCAL_ITER = J2_PLANE_STRESS_MAX_LOCAL_ITER
# Piso para denominadores en la corrección de tangente del return mapping.
_DENOM_FLOOR = J2_DENOM_FLOOR


# Constantes del tangente J2 en Voigt 2D: v = gradiente de la traza,
# I_dev = proyector desviador (½ en el cortante engineering).
_V2 = (1.0, 1.0, 0.0)
_I_DEV2 = ((2.0 / 3.0, -1.0 / 3.0, 0.0),
           (-1.0 / 3.0, 2.0 / 3.0, 0.0),
           (0.0, 0.0, 0.5))


@njit(cache=True)
def _j2_plane_strain_core(strain, eps_p_old, alpha_old, sigma_y, H, K, G, C_e,
                          yield_tol, sigma, C_out, eps_p_new):
    """Return mapping J2 plane strain **sin asignar memoria**: escribe
    ``sigma`` (3), ``C_out`` (3×3) y ``eps_p_new`` (4) in situ y devuelve
    ``alpha_new``. Descomposición volumétrica-desviadora 3D con ``ε_zz = 0``
    impuesto en la deformación total; ``ε^p_zz`` evoluciona libremente y
    queda registrada en ``eps_p`` (cuarta componente del estado).

    Es la formulación de siempre escrita con escalares y **en el mismo
    orden de operaciones** que la versión con arreglos que sustituye
    (ADR 0014): los resultados coinciden bit a bit. Lo consumen
    :func:`_compute_j2_plane_strain` (camino por elemento) y el adaptador
    por lotes.
    """
    # Deformación volumétrica (traza, ε_zz = 0 por plane strain)
    eps_v = strain[0] + strain[1]
    third = eps_v / 3.0

    # Desviador de la deformación 3D extendida [xx, yy, zz, xy_tensorial]
    # menos la deformación plástica: predictor elástico.
    e0 = (strain[0] - third) - eps_p_old[0]
    e1 = (strain[1] - third) - eps_p_old[1]
    e2 = (-third) - eps_p_old[2]
    e3 = (strain[2] / 2.0) - eps_p_old[3]
    twoG = 2.0 * G
    s0 = twoG * e0
    s1 = twoG * e1
    s2 = twoG * e2
    s3 = twoG * e3
    norm_s_trial = math.sqrt(s0 * s0 + s1 * s1 + s2 * s2 + 2.0 * (s3 * s3))

    yield_stress = sigma_y + H * alpha_old
    f_trial = norm_s_trial - math.sqrt(2.0 / 3.0) * yield_stress

    p = K * eps_v

    if f_trial <= yield_tol:
        # Predictor elástico σ = s_trial + p·I (NO C_e·strain, que ignoraría ε_p).
        # Diferencia material en presencia de plasticidad acumulada (descarga o
        # reevaluación post-converged): C_e·strain devolvería un esfuerzo
        # incompatible con el estado interno; s_trial + p·I respeta ε_p_old.
        sigma[0] = s0 + p
        sigma[1] = s1 + p
        sigma[2] = s3
        for i in range(3):
            for j in range(3):
                C_out[i, j] = C_e[i, j]
        for i in range(4):
            eps_p_new[i] = eps_p_old[i]
        return alpha_old

    # Corrector plástico — return mapping radial cerrado
    delta_gamma = f_trial / (2.0 * G + (2.0 / 3.0) * H)
    N0 = s0 / norm_s_trial
    N1 = s1 / norm_s_trial
    N2 = s2 / norm_s_trial
    N3 = s3 / norm_s_trial

    c = twoG * delta_gamma
    sn0 = s0 - c * N0
    sn1 = s1 - c * N1
    sn3 = s3 - c * N3
    eps_p_new[0] = eps_p_old[0] + delta_gamma * N0
    eps_p_new[1] = eps_p_old[1] + delta_gamma * N1
    eps_p_new[2] = eps_p_old[2] + delta_gamma * N2
    eps_p_new[3] = eps_p_old[3] + delta_gamma * N3
    alpha_new = alpha_old + math.sqrt(2.0 / 3.0) * delta_gamma

    sigma[0] = sn0 + p
    sigma[1] = sn1 + p
    sigma[2] = sn3

    # Matriz tangente algorítmica consistente
    # C_alg = K·v⊗v + 2G(1−β)·I_dev − 2G·γ·N⊗N   (N en Voigt: [N0, N1, N3])
    beta = c / norm_s_trial
    gamma_factor = 1.0 / (1.0 + H / (3.0 * G)) - beta
    a2 = twoG * (1.0 - beta)
    a3 = twoG * gamma_factor
    Nv = (N0, N1, N3)
    for i in range(3):
        for j in range(3):
            C_out[i, j] = (K * (_V2[i] * _V2[j]) + a2 * _I_DEV2[i][j]) - a3 * (Nv[i] * Nv[j])
    return alpha_new


@njit(cache=True)
def _compute_j2_plane_strain(strain, eps_p_old, alpha_old, sigma_y, H, K, G, C_e, yield_tol):
    """Return mapping J2 plane strain (camino por elemento): envoltorio con
    asignación de :func:`_j2_plane_strain_core`. Devuelve
    ``(σ, C_alg, ε^p_new, α_new)``."""
    sigma = np.empty(3)
    C_alg = np.empty((3, 3))
    eps_p_new = np.empty(4)
    alpha_new = _j2_plane_strain_core(strain, eps_p_old, alpha_old, sigma_y, H, K, G, C_e,
                                      yield_tol, sigma, C_alg, eps_p_new)
    return sigma, C_alg, eps_p_new, alpha_new


@njit(cache=True)
def _compute_j2_plane_stress(strain, eps_p_old, alpha_old,
                             sigma_y, H, E, nu, G, C_e_ps,
                             tol_rel, tol_abs, max_local_iter):
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

    Tolerancia del Newton local (auditoría 2026-09-22): se recalcula en cada
    iteración con la escala corriente ``R_curr²/3`` (``tol_abs + tol_rel·R²/3``,
    ADR 0006). Con ``H > 0`` y un incremento plástico grande, ``R_curr ≫ R_old``
    y una tolerancia fija en ``α_old`` quedaba por debajo del redondeo de ``f̄``
    (~ε_mach·R_curr²), de modo que el bucle agotaba siempre ``max_local_iter``.
    El punto de partida es la cota inferior de la raíz para ``H = 0``
    (``Δγ₀ = (√(3·½σPσ_trial)/R − 1)/μ_max``): desde ``Δγ = 0`` el Newton
    necesitaba ~log₁.₅(σ_trial/R) iteraciones y con predictores lejanos
    (``σ_trial/R ≳ 3·10³``) salía sin converger y sin avisar. Devuelve un
    quinto valor ``converged_local`` que el envoltorio Python registra.
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

    if f_bar_trial <= tol_abs + tol_rel * (R_trial * R_trial) / 3.0:
        return sigma_trial.copy(), C_e_ps.copy(), eps_p_old.copy(), alpha_old, True

    # Autovalores de C_e^ps · P
    mu1 = E / (3.0 * (1.0 - nu))
    mu_dev = 2.0 * G  # μ₂ = μ₃ = 2G

    # Proyecciones de σ_trial sobre autovectores ortonormales:
    # v₁ = [1,1,0]/√2 (hidrostático plano), v₂ = [1,-1,0]/√2 (desviador plano), v₃ = [0,0,1]
    inv_sqrt2 = 1.0 / math.sqrt(2.0)
    a1 = (sigma_trial[0] + sigma_trial[1]) * inv_sqrt2
    a2 = (sigma_trial[0] - sigma_trial[1]) * inv_sqrt2
    a3 = sigma_trial[2]

    # Newton local sobre Δγ. Warm start: cota inferior de la raíz para H = 0
    # (f̄ es convexa decreciente en Δγ, así que desde la izquierda el Newton
    # es monótono). Ver docstring.
    mu_max = mu1 if mu1 > mu_dev else mu_dev
    ratio = math.sqrt(3.0 * half_sPs_trial) / R_trial
    delta_gamma = (ratio - 1.0) / mu_max if ratio > 1.0 else 0.0
    converged_local = False

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

        if abs(f_bar) < tol_abs + tol_rel * (R_curr * R_curr) / 3.0:
            converged_local = True
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
        if w > _DENOM_FLOOR:
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
    if w_final > _DENOM_FLOOR:
        n = 2.0 * delta_gamma / (3.0 * w_final)
    else:
        n = 0.0

    one_minus_mn = 1.0 - m * n
    # Forma agrupada equivalente, numéricamente estable cuando m·n → 1:
    # denom = q + m·w / (1 - m·n)
    if abs(one_minus_mn) < _DENOM_FLOOR or H == 0.0:
        denom_beta = q  # H=0 → no hay corrección por endurecimiento
    else:
        denom_beta = q + m * w_final / one_minus_mn

    if abs(denom_beta) < _DENOM_FLOOR:
        C_alg = D.copy()
    else:
        C_alg = D - np.outer(D_P_sigma, sigma_P_D) / denom_beta

    return sigma_new, C_alg, eps_p_new, alpha_new, converged_local


@njit(cache=True)
def _j2_plane_strain_batch(strain, S_old, S_new, params, C, sigma, C_out, flag):
    """Adaptador por lotes (ADR 0014) del return mapping plane strain.

    ``params = [σ_y, H, K, G, tol_abs, tol_rel]``; ``C`` = ``C_e``. La
    tolerancia de fluencia se evalúa aquí con las mismas operaciones que
    ``Material.admissibility_tol`` (ADR 0006), para que el camino por
    lotes coincida bit a bit con el camino por elemento. Sin asignaciones
    por punto de Gauss: escribe ``sigma``, ``C_out`` y la fila ``S_new``.
    """
    sigma_y = params[0]
    H = params[1]
    alpha_old = S_old[4]
    R = sigma_y + H * alpha_old
    yield_tol = params[4] + params[5] * (math.sqrt(2.0 / 3.0) * R)
    S_new[4] = _j2_plane_strain_core(
        strain, S_old[0:4], alpha_old, sigma_y, H, params[2], params[3], C, yield_tol,
        sigma, C_out, S_new[0:4],
    )
    return C_out


@njit(cache=True)
def _j2_plane_stress_batch(strain, S_old, S_new, params, C, sigma, C_out, flag):
    """Adaptador por lotes (ADR 0014) del return mapping plane stress.

    ``params = [σ_y, H, E, ν, G, tol_rel, tol_abs, max_local_iter]``;
    ``C`` = ``C_e`` plane stress. Marca ``flag = 1`` si el Newton local no
    convergió; el material lo reporta con ``batch_report``.
    """
    sig, C_alg, eps_p_new, alpha_new, local_ok = _compute_j2_plane_stress(
        strain, S_old[0:4], S_old[4],
        params[0], params[1], params[2], params[3], params[4], C,
        params[5], params[6], int(params[7]),
    )
    if not local_ok:
        flag[0] = 1
    for i in range(3):
        sigma[i] = sig[i]
    for i in range(4):
        S_new[i] = eps_p_new[i]
    S_new[4] = alpha_new
    return C_alg


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
        Esfuerzo de fluencia inicial (uniaxial).
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
    # eps_p = [ε^p_xx, ε^p_yy, ε^p_zz, ε^p_xy_tensorial]; alpha escalar.
    STATE_SCHEMA = {'eps_p': (4,), 'alpha': ()}

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
        if density is not None and density < 0.0:
            raise ValueError(
                f"VonMises2D: density={density} no puede ser negativa."
            )

        self.E = E
        self.nu = nu
        self.sigma_y = sigma_y
        self.H = H
        self.hypothesis = hypothesis
        self.density = density
        self._local_newton_warned = False

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
            sigma, C_alg, eps_p_new, alpha_new, local_ok = _compute_j2_plane_stress(
                strain, eps_p_old, alpha_old,
                self.sigma_y, self.H, self.E, self.nu, self.G, self.C_e,
                ADMISSIBILITY_TOL_REL, ADMISSIBILITY_TOL_ABS,
                _PLANE_STRESS_MAX_LOCAL_ITER,
            )
            if not local_ok and not self._local_newton_warned:
                # Aviso único por instancia: el residuo global queda
                # contaminado en este iterado; el Newton global suele
                # recuperarse (bisección), pero no debe pasar en silencio.
                self._local_newton_warned = True
                _log.warning(
                    f"VonMises2D plane stress: el Newton local del return mapping "
                    f"agotó {_PLANE_STRESS_MAX_LOCAL_ITER} iteraciones sin converger "
                    f"(predictor muy lejano de la superficie de fluencia). σ puede "
                    f"quedar fuera de la superficie en este iterado."
                )

        new_state = {'eps_p': eps_p_new, 'alpha': alpha_new}
        return sigma, C_alg, new_state

    # ------------------------------------------------------------------
    # Camino por lotes (ADR 0014)
    # ------------------------------------------------------------------

    def batch_kernel(self):
        if self.hypothesis == 'plane_strain':
            return _j2_plane_strain_batch
        return _j2_plane_stress_batch

    def batch_params(self) -> np.ndarray:
        if self.hypothesis == 'plane_strain':
            return np.array([self.sigma_y, self.H, self.K, self.G,
                             ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL], dtype=np.float64)
        return np.array([self.sigma_y, self.H, self.E, self.nu, self.G,
                         ADMISSIBILITY_TOL_REL, ADMISSIBILITY_TOL_ABS,
                         float(_PLANE_STRESS_MAX_LOCAL_ITER)], dtype=np.float64)

    def batch_matrix(self) -> np.ndarray:
        return self.C_e

    def batch_report(self, n_flagged: int) -> None:
        if not self._local_newton_warned:
            self._local_newton_warned = True
            _log.warning(
                f"VonMises2D plane stress: el Newton local del return mapping "
                f"agotó {_PLANE_STRESS_MAX_LOCAL_ITER} iteraciones sin converger "
                f"en {n_flagged} punto(s) de Gauss (predictor muy lejano de la "
                f"superficie de fluencia). σ puede quedar fuera de la superficie "
                f"en este iterado."
            )

    def out_of_plane_stress(self, sigma, state_vars=None) -> float:
        """``σ_zz`` en plane strain a partir de ``σ`` en plano y de ``ε^p_zz``
        (ver ``solidum.materials._plane_strain``); ``0`` en plane stress."""
        if self.hypothesis != 'plane_strain':
            return 0.0
        eps_p_zz = 0.0 if state_vars is None else float(
            state_vars.get('eps_p', np.zeros(4))[2]
        )
        return sigma_zz_plane_strain(sigma[0], sigma[1], eps_p_zz, self.K, self.G)

    def batch_out_of_plane_stress(self, sigma, S):
        """Versión por arreglos: ``ε^p_zz`` es la columna 2 de la fila de
        estado (``eps_p`` ocupa las columnas 0–3 del ``STATE_SCHEMA``)."""
        if self.hypothesis != 'plane_strain':
            return np.zeros(sigma.shape[0])
        return sigma_zz_plane_strain(sigma[:, 0], sigma[:, 1], S[:, 2], self.K, self.G)
