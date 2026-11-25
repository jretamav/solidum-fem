# solidum_fem/solidum/materials/fixed_smeared_crack_2d.py
"""``FixedSmearedCrack2D`` — grieta distribuida fija con retención de cortante (modelo 3).

Modelo 3 de la familia de concreto (spec ``docs/specs/FixedSmearedCrack2D.md``):
descomposición de la deformación en la del concreto entre grietas y la de la
grieta, ``ε = C_e⁻¹·σ + N·e`` (Rots 1988; de Borst et al. 2012, ecs.
6.127–6.135), con la grieta fija en la dirección principal mayor del esfuerzo
de prueba al agrietarse y una segunda grieta perpendicular. La grieta se
reparte en el ancho ``h`` del elemento (banda de grieta): apertura ``w = h·e``
y deslizamiento ``s = h·γ``.

- Ley normal (§2): ablandamiento ``f(w)`` exponencial en la apertura o de
  Hordijk (``_softening``), descarga secante hacia el origen de la grieta,
  cierre con ``w = 0``.
- Ley de cortante (§3): ``τ/s = f(w_max)/(2·w_max)``, sin parámetro; las dos
  grietas en serie, ``c_s = 2·(c_n + c_t)``. Sin recuperar rigidez al cerrar.
- Tope (§6): la flexibilidad ``c = w_max/f(w_max)`` de cada dirección no pasa
  de ``c_max = h·ω_max/((1 − ω_max)·E)`` (``damage_max``).

Bloque 1 (sin aplastamiento, §8–§9): el concreto entre grietas es elástico. En
los ejes de la grieta el cortante se desacopla del problema normal, que tiene
dos incógnitas, las aperturas ``(w_n, w_t)``:

    A·(ε_N − w/h) ∈ F(w),     A = [[C11, C12], [C12, C11]],

con ``F_i`` el grafo monótono de cada dirección (cerrada en ``w = 0``, secante
hasta ``w_max``, envolvente convexa después). Bajo el límite de ``h`` de §2 el
operador ``F(w) + A·w/h`` es estrictamente monótono y la solución es única. Se
resuelve anidado, con un problema escalar monótono en cada nivel:

1. Dentro, ``w_n`` dada ``w_t``: la rama se decide evaluando en los quiebres
   (``w = 0``, ``w = w_max``, el tope); en la secante, forma cerrada; en la
   envolvente, Newton desde la derecha, que baja de forma monótona porque la
   función es creciente y convexa (como ``_softening.evaluate_hordijk_damage``).
2. Fuera, ``w_t``: la función reducida ``R(w_t) = F_t(w_t) − σ_tt(w_n*(w_t),
   w_t)`` es creciente (minimización parcial de una función estrictamente
   convexa). Se acota entre la raíz con ``w_n`` congelada en ``w_n*(0)`` (cota
   inferior) y la raíz con ``w_n = 0`` (cota superior) y se resuelve con
   Newton protegido por bisección.

El cortante sale en forma cerrada con la historia ya actualizada,
``τ = γ_nt/(1/G + c_s/h)``. Tangente consistente (§8): la del problema normal
más la fila de cortante con ``∂τ/∂e_ii``; no simétrica en carga con cortante.

IMPL-EX (§9; ADR 0034): las aperturas máximas se extrapolan con los dos
últimos pasos confirmados y las flexibilidades quedan fijas en el paso (no
menores que las del paso anterior, guardadas como ``damage_n``/``damage_t``); el
esfuerzo sale de la ley secante en forma cerrada y la tangente es la secante,
simétrica. La historia se actualiza con la integración implícita.

Bloque 2 (aplastamiento, §4, §8b, §9; con ``f_c`` y ``G_c``): el concreto
entre grietas es el modelo 2 con ``d⁺ = 0``, ``σ_b = σ̄ − d⁻·σ̄⁻`` con
``σ̄ = C_e·ε^co`` (``_tension_compression``: Drucker-Prager sobre la parte
negativa y parábola de Feenstra con ``G_c/h``). ``d⁻`` es función explícita de
las deformaciones de la grieta a través de ``ε^co = ε − N·e``, así que el
problema local se acopla:

1. Predictor: con ``d⁻_n = 0`` se resuelve el bloque 1; si con su esfuerzo
   (el efectivo, porque el concreto es elástico) ``d⁻(κ⁻) = 0``, es la
   solución: sin aplastamiento todo es el bloque 1 bit a bit (salvo ``κ⁻``).
2. Si no, Newton semiliso sobre ``x = [w_n, w_t, s]`` (complementariedad de
   cada dirección con la función mínimo, cortante en flexibilidad) con la
   tangente algorítmica del modelo 2, en cascada (paso completo y búsqueda
   lineal desde el paso confirmado; después desde ``x = 0``). Si nada
   converge, esfuerzo y tangente no finitos: el corrector abandona el intento
   y el solver repite el paso más corto (Q5 del usuario: selección por
   continuidad, sin límite ni aviso nuevos).

Con ``d̃⁻ > 0`` IMPL-EX pierde la forma cerrada: el mismo Newton con la
historia fija y tangente consistente no simétrica.

Caché de Numba: este núcleo se compila con ``hordijk_shape`` y
``hordijk_shape_slope`` de ``_softening`` dentro; al cambiar ese módulo hay que
tocar también éste.
"""
import math

import numpy as np
from numba import njit

from solidum.constants import ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL
from solidum.core.material import Material
from solidum.materials._crack_band import (
    CRACK_BAND_DAMAGE_MAX,
    DEFAULT_INTEGRATION,
    DEFAULT_SOFTENING,
    IMPL_EX,
    NO_STEP_RATIO,
    CrackBandMixin,
    parse_integration,
    validate_crack_band_parameters,
)
from solidum.materials._softening import (
    HORDIJK,
    HORDIJK_H_MAX_FACTOR,
    hordijk_shape,
    hordijk_shape_slope,
)
from solidum.materials._spectral_split import spectral_split_2d
from solidum.materials._tension_compression import (
    DEFAULT_BETA,
    compression_damage,
    compression_kernel_params,
    compression_measure,
    compression_ultimate,
    parse_compression,
)
from solidum.materials.elastic_2d import Elastic2D
from solidum.registry import MaterialRegistry

# Índices de la fila de estado (orden de ``STATE_SCHEMA``, el del contrato de la spec).
(CRACKED, THETA, WMAX_N, WMAX_T, W_N, W_T, SLIP, DMG_N, DMG_T, KC, DC, H,
 WMAX_N_PREV, WMAX_T_PREV, KC_PREV, ERR) = range(16)
STATE_SCHEMA_FSC = {name: () for name in (
    'cracked', 'theta', 'w_max_n', 'w_max_t', 'w_n', 'w_t', 'slip', 'damage_n',
    'damage_t', 'kappa_c', 'damage_c', 'h', 'w_max_n_prev', 'w_max_t_prev',
    'kappa_c_prev', 'damage_error')}

# Ley de una dirección en el problema normal.
#   IMPLICIT_LAW: la de §2 con la historia confirmada (sin grieta si w_max = 0:
#                 cerrada hasta f_t; con grieta: cerrada, secante hasta w_max y
#                 envolvente después).
#   SECANT_LAW:   la secante fija c > 0 para toda w > 0 (la de IMPL-EX, o la
#                 rama de descarga ya decidida).
#   RIGID_LAW:    sin flexibilidad (IMPL-EX con c̃ = 0): w = 0 siempre.
#   ENVELOPE_LAW: sólo la envolvente con el tope (la rama de carga ya decidida).
IMPLICIT_LAW, SECANT_LAW, RIGID_LAW, ENVELOPE_LAW = 0, 1, 2, 3
# Rama de una dirección en la solución.
CLOSED, SECANT, LOADING, CAPPED = 0, 1, 2, 3
# Tope de iteraciones de los Newton locales (convergen en pocas; el tope sólo
# protege de un bucle infinito por redondeo).
MAX_ITER = 200
_EPS = np.finfo(np.float64).eps


# ---------------------------------------------------------------------------
# Ley de la grieta en una dirección (§2, §6)
# ---------------------------------------------------------------------------

@njit(cache=True)
def _envelope(w, f_t, G_f, law, w_c):
    """``f(w)`` y ``f'(w)`` de la ley de ablandamiento (§2), ``w ≥ 0``.

    Exponencial en la apertura, ``f_t·exp(−f_t·w/G_f)``, o la de Hordijk,
    ``f_t·g(w/w_c)``, nula desde ``w_c``. Las dos son convexas y decrecientes.
    """
    if law == HORDIJK:
        if w >= w_c:
            return 0.0, 0.0
        x = w / w_c
        return f_t * hordijk_shape(x), f_t * hordijk_shape_slope(x) / w_c
    s = f_t * math.exp(-f_t * w / G_f)
    return s, -(f_t / G_f) * s


@njit(cache=True)
def _compliance(w_max, f_t, G_f, law, w_c, c_max):
    """Flexibilidad secante ``c = w_max/f(w_max)`` con el tope ``c_max`` (§6)
    y su derivada ``dc/dw = (f − w·f')/f²`` (0 en el tope y sin grieta).

    El tope se aplica antes de dividir: con Hordijk y ``w ≥ w_c``, ``f = 0`` y
    ``w ≥ c_max·f`` da ``c_max`` sin dividir entre cero.
    """
    if w_max <= 0.0:
        return 0.0, 0.0
    f, df = _envelope(w_max, f_t, G_f, law, w_c)
    if w_max >= c_max * f:
        return c_max, 0.0
    return w_max / f, (f - w_max * df) / (f * f)


@njit(cache=True)
def _graph(w, mode, c, f_t, G_f, law, w_c, c_max):
    """``F(w)``, ``F'(w)`` y la rama de una dirección abierta (``w > 0``) en
    una rama ya decidida: la secante (``SECANT_LAW``) o la envolvente con el
    tope, ``max(f, w/c_max)`` (``ENVELOPE_LAW``)."""
    if mode == SECANT_LAW:
        return w / c, 1.0 / c, SECANT
    f, df = _envelope(w, f_t, G_f, law, w_c)
    if w >= c_max * f:
        return w / c_max, 1.0 / c_max, CAPPED
    return f, df, LOADING


@njit(cache=True)
def _envelope_root(a, b, w_lo, f_t, G_f, law, w_c, c_max):
    """Raíz ``w ≥ w_lo`` de ``E(w) + a·w = b`` sobre la envolvente con el tope,
    ``E = max(f, w/c_max)`` (§6). Devuelve ``(w, rama, E'(w))``; si
    ``E(w_lo) + a·w_lo ≥ b``, ``w_lo``.

    Con el tope, ``w = b/(a + 1/c_max)`` si ahí ``w ≥ c_max·f``. Si no, la raíz
    está antes del tope, donde ``r(w) = f(w) + a·w − b`` es creciente y convexa
    (``a > |f'(0)|``, límite de ``h``): Newton desde la derecha (la tangente en
    ``w_lo`` cae a la derecha de la raíz) baja de forma monótona hasta que deja
    de decrecer (redondeo), como ``_softening.evaluate_hordijk_damage``.
    """
    w_cap = b / (a + 1.0 / c_max)
    f, _ = _envelope(w_cap, f_t, G_f, law, w_c)
    if w_cap >= c_max * f:
        if w_cap > w_lo:
            return w_cap, CAPPED, 1.0 / c_max
        return w_lo, CAPPED, 1.0 / c_max
    f0, df0 = _envelope(w_lo, f_t, G_f, law, w_c)
    r0 = f0 + a * w_lo - b
    if r0 >= 0.0:
        return w_lo, LOADING, df0
    w = w_lo - r0 / (df0 + a)
    if not (w < w_cap):
        w = w_cap
    for _ in range(MAX_ITER):
        f, df = _envelope(w, f_t, G_f, law, w_c)
        w_new = w - (f + a * w - b) / (df + a)
        if not (w_new < w):
            break
        w = w_new
    f, df = _envelope(w, f_t, G_f, law, w_c)
    return w, LOADING, df


@njit(cache=True)
def _solve_direction(a, b, mode, w_max, c, f_t, G_f, law, w_c, c_max, tol):
    """Raíz ``w ≥ 0`` de ``F(w) + a·w ∋ b`` en una dirección. Devuelve
    ``(w, rama, F'(w))``.

    ``b`` es el esfuerzo normal del concreto con esta apertura nula y ``a`` la
    rigidez del concreto frente a la apertura (``C11/h``). La rama se decide
    evaluando en los quiebres (§8), con la tolerancia de admisibilidad ``tol``
    del ADR 0006, y no se vuelve a clasificar después:

    - ``RIGID_LAW``: cerrada.
    - Sin grieta (``IMPLICIT_LAW`` con ``w_max = 0``): cerrada si
      ``b ≤ f_t + tol``; si no, la envolvente desde 0.
    - Con grieta: cerrada si ``b ≤ tol``; secante, en forma cerrada, con
      ``SECANT_LAW`` o si el esfuerzo de prueba no supera la envolvente en
      ``w_max`` en más de ``tol`` (``b ≤ (1/c + a)·w_max + tol``); si no, la
      envolvente desde ``w_max``.

    Sin la tolerancia la decisión dependería del redondeo en el quiebre: al
    empezar un paso el punto se evalúa en su estado confirmado (``w = w_max``),
    con ``b`` igual a ``f(w_max) + a·w_max`` salvo el redondeo, y unos puntos
    caerían en la carga y otros en la descarga, con tangentes de signo
    contrario (en la barra de Quad8 eso rompía la simetría del predictor y
    llevaba a Newton a la rama localizada dentro del elemento).
    """
    if mode == RIGID_LAW:
        return 0.0, CLOSED, 0.0
    if mode == IMPLICIT_LAW and w_max <= 0.0:
        if b <= f_t + tol:
            return 0.0, CLOSED, 0.0
        return _envelope_root(a, b, 0.0, f_t, G_f, law, w_c, c_max)
    if b <= tol:
        return 0.0, CLOSED, 0.0
    if mode == SECANT_LAW or b <= (1.0 / c + a) * w_max + tol:
        return b / (1.0 / c + a), SECANT, 1.0 / c
    return _envelope_root(a, b, w_max, f_t, G_f, law, w_c, c_max)


@njit(cache=True)
def _branch_root(a, b, mode, c, w_lo, f_t, G_f, law, w_c, c_max):
    """Raíz de ``F(w) + a·w = b`` en una rama ya decidida (cotas del problema
    de la segunda dirección): la secante desde 0 o la envolvente desde ``w_lo``."""
    if mode == SECANT_LAW:
        return max(b, 0.0) / (1.0 / c + a)
    w, _, _ = _envelope_root(a, b, w_lo, f_t, G_f, law, w_c, c_max)
    return w


# ---------------------------------------------------------------------------
# Problema local en los ejes de la grieta (§8)
# ---------------------------------------------------------------------------

@njit(cache=True)
def _normal_problem(e_n, e_t, C11, C12, h, f_t, G_f, law, w_c, c_max, tol,
                    mode_n, wmax_n, c_n, mode_t, wmax_t, c_t):
    """Aperturas del problema normal de dos incógnitas, con sus ramas.

    ``e_n``, ``e_t``: deformaciones normales totales en los ejes de la grieta.
    Anidado (docstring del módulo): ``w_n*(w_t)`` por :func:`_solve_direction`
    y ``w_t`` como raíz de ``R(w_t) = F_t(w_t) + a·w_t + q·w_n* − b_t0``, con
    ``a = C11/h``, ``q = C12/h`` y ``b_t0 = C12·e_n + C11·e_t``. La rama de
    ``t`` se decide una vez, evaluando ``R`` en el quiebre ``w_t = w_max_t``
    con la tolerancia, y la raíz se busca en esa rama. Devuelve
    ``(w_n, rama_n, F_n', w_t, rama_t, F_t')``.
    """
    a = C11 / h
    q = C12 / h
    b_n0 = C11 * e_n + C12 * e_t
    b_t0 = C12 * e_n + C11 * e_t
    w_n, br_n, K_n = _solve_direction(a, b_n0, mode_n, wmax_n, c_n, f_t, G_f, law, w_c,
                                      c_max, tol)
    if mode_t == RIGID_LAW:
        return w_n, br_n, K_n, 0.0, CLOSED, 0.0
    s_tt = b_t0 - q * w_n
    virgin_t = mode_t == IMPLICIT_LAW and wmax_t <= 0.0
    if s_tt <= (f_t + tol if virgin_t else tol):
        return w_n, br_n, K_n, 0.0, CLOSED, 0.0
    t_law = SECANT_LAW if mode_t == SECANT_LAW else ENVELOPE_LAW
    w_lo = 0.0
    if mode_t == IMPLICIT_LAW and wmax_t > 0.0:
        w_at, _, _ = _solve_direction(a, b_n0 - q * wmax_t, mode_n, wmax_n, c_n, f_t, G_f,
                                      law, w_c, c_max, tol)
        if wmax_t / c_t + a * wmax_t + q * w_at - b_t0 >= -tol:
            t_law = SECANT_LAW
        else:
            w_lo = wmax_t
    # w_n*(w_t) no crece con w_t: la raíz queda entre la de w_n = w_n*(0)
    # (cota inferior) y la de w_n = 0 (cota superior).
    lo = _branch_root(a, s_tt, t_law, c_t, w_lo, f_t, G_f, law, w_c, c_max)
    hi = _branch_root(a, b_t0, t_law, c_t, w_lo, f_t, G_f, law, w_c, c_max)
    w = lo
    if hi > lo:
        for _ in range(MAX_ITER):
            w_n, br_n, K_n = _solve_direction(a, b_n0 - q * w, mode_n, wmax_n, c_n, f_t, G_f,
                                              law, w_c, c_max, tol)
            F, dF, _ = _graph(w, t_law, c_t, f_t, G_f, law, w_c, c_max)
            R = F + a * w + q * w_n - b_t0
            if R < 0.0:
                lo = w
            elif R > 0.0:
                hi = w
            else:
                break
            if hi - lo <= 4.0 * _EPS * hi:
                break
            dR = dF + a
            if br_n != CLOSED:
                dR -= q * q / (K_n + a)
            w_new = w - R / dR
            if not (lo < w_new < hi):
                w_new = 0.5 * (lo + hi)
            if abs(w_new - w) <= 4.0 * _EPS * w:
                break
            w = w_new
    w_n, br_n, K_n = _solve_direction(a, b_n0 - q * w, mode_n, wmax_n, c_n, f_t, G_f, law,
                                      w_c, c_max, tol)
    _, K_t, br_t = _graph(w, t_law, c_t, f_t, G_f, law, w_c, c_max)
    return w_n, br_n, K_n, w, br_t, K_t


@njit(cache=True)
def _crack_response(eps_loc, C11, C12, G, h, f_t, G_f, law, w_c, c_max, tol,
                    mode_n, wmax_n, c_n, mode_t, wmax_t, c_t, implicit):
    """Esfuerzo y tangente en los ejes de la grieta, ``[σ_nn, σ_tt, τ_nt]``.

    ``implicit``: la historia crece sólo en la carga sobre la envolvente
    (``w_max_i ← w_i``), el cortante usa las flexibilidades de la historia
    nueva y la tangente lleva ``∂τ/∂e_ii``; si no (IMPL-EX), las
    flexibilidades ``c_n``, ``c_t`` dadas. Devuelve ``(σ_loc, C_loc, w_n,
    w_t, s, w_max_n, w_max_t)`` con ``s = c_s·τ`` el deslizamiento.
    """
    e_n = eps_loc[0]
    e_t = eps_loc[1]
    w_n, br_n, K_n, w_t, br_t, K_t = _normal_problem(
        e_n, e_t, C11, C12, h, f_t, G_f, law, w_c, c_max, tol,
        mode_n, wmax_n, c_n, mode_t, wmax_t, c_t)
    dc_n = 0.0
    dc_t = 0.0
    wn_new = wmax_n
    wt_new = wmax_t
    cn_new = c_n
    ct_new = c_t
    if implicit:
        if br_n == LOADING or br_n == CAPPED:
            wn_new = w_n
            cn_new, dc_n = _compliance(w_n, f_t, G_f, law, w_c, c_max)
        if br_t == LOADING or br_t == CAPPED:
            wt_new = w_t
            ct_new, dc_t = _compliance(w_t, f_t, G_f, law, w_c, c_max)
    c_s = 2.0 * (cn_new + ct_new)
    D = 1.0 / G + c_s / h
    tau = eps_loc[2] / D
    sig = np.empty(3)
    sig[0] = C11 * (e_n - w_n / h) + C12 * (e_t - w_t / h)
    sig[1] = C12 * (e_n - w_n / h) + C11 * (e_t - w_t / h)
    sig[2] = tau
    # Tangente: dw = S·dε_N en las direcciones abiertas (S = M⁻¹·A_a,
    # M = diag(K_a) + A_aa/h); dσ_N = A·(dε_N − dw/h).
    a = C11 / h
    q = C12 / h
    S = np.zeros((2, 2))
    if br_n != CLOSED and br_t != CLOSED:
        m11 = K_n + a
        m22 = K_t + a
        det = m11 * m22 - q * q
        S[0, 0] = (m22 * C11 - q * C12) / det
        S[0, 1] = (m22 * C12 - q * C11) / det
        S[1, 0] = (m11 * C12 - q * C11) / det
        S[1, 1] = (m11 * C11 - q * C12) / det
    elif br_n != CLOSED:
        S[0, 0] = C11 / (K_n + a)
        S[0, 1] = C12 / (K_n + a)
    elif br_t != CLOSED:
        S[1, 0] = C12 / (K_t + a)
        S[1, 1] = C11 / (K_t + a)
    C_loc = np.zeros((3, 3))
    for j in range(2):
        C_loc[0, j] = (C11 if j == 0 else C12) - (C11 * S[0, j] + C12 * S[1, j]) / h
        C_loc[1, j] = (C12 if j == 0 else C11) - (C12 * S[0, j] + C11 * S[1, j]) / h
        # ∂τ/∂ε_N = −τ/(D·h)·Σ_i (dc_s/dw_i)·S_i, con dc_s/dw_i = 2·dc_i/dw_i
        # (sólo la carga sobre la envolvente cambia la historia; 0 en el tope).
        C_loc[2, j] = -tau / (D * h) * 2.0 * (dc_n * S[0, j] + dc_t * S[1, j])
    C_loc[2, 2] = 1.0 / D
    return sig, C_loc, w_n, w_t, c_s * tau, wn_new, wt_new


@njit(cache=True)
def _axes(theta):
    """``T`` con ``ε_loc = T·ε`` (Voigt con γ de ingeniería) y ``σ = Tᵀ·σ_loc``.

    ``n = (cos θ, sen θ)``, ``t = (−sen θ, cos θ)``; ``T = N⁻¹`` con la ``N``
    de la spec §1.
    """
    c = math.cos(theta)
    s = math.sin(theta)
    T = np.empty((3, 3))
    T[0, 0] = c * c
    T[0, 1] = s * s
    T[0, 2] = c * s
    T[1, 0] = s * s
    T[1, 1] = c * c
    T[1, 2] = -c * s
    T[2, 0] = -2.0 * c * s
    T[2, 1] = 2.0 * c * s
    T[2, 2] = c * c - s * s
    return T


@njit(cache=True)
def _cap_compliance(h, E, cap):
    """``c_max = h·ω_max/((1 − ω_max)·E)`` (§6); 0 sin ``h`` sembrada."""
    if h <= 0.0:
        return 0.0
    return h * cap / ((1.0 - cap) * E)


@njit(cache=True)
def _omega(c, h, E):
    """Daño uniaxial equivalente de la banda, ``ω = E·c/(h + E·c)`` (§6)."""
    return E * c / (h + E * c)


@njit(cache=True)
def _compliance_from_omega(omega, h, E):
    """Inversa de :func:`_omega`: ``c = h·ω/(E·(1 − ω))``."""
    return h * omega / (E * (1.0 - omega))


# ---------------------------------------------------------------------------
# Bloque 2: concreto entre grietas con aplastamiento (§4) y problema local
# acoplado (§8b, §9)
# ---------------------------------------------------------------------------

# Escalares del problema acoplado (arreglo ``P``).
(P_E, P_NU, P_H, P_K, P_G, P_FT, P_GF, P_WC, P_CMAX, P_FC, P_KDP, P_KU, P_TOLC,
 P_CAP, P_KCN, P_DFIX, P_KTOL) = range(17)
# Enteros del problema acoplado (arreglo ``I``): deformación plana, ley de
# ablandamiento, concreto implícito (1) o con d⁻ fijo (0), ley de cada dirección.
I_PS, I_LAW, I_IMPL, I_MODE_N, I_MODE_T = range(5)
# Historia de cada dirección (arreglo ``Hd``): w_max, c y la banda δ de la
# tolerancia en el quiebre w = w_max.
HD_WMAX, HD_C, HD_DELTA = 0, 2, 4
# Tope de iteraciones de cada intento del Newton acoplado (el teórico midió un
# máximo de 31 con pasos grandes) y de retrocesos de su búsqueda lineal.
MAX_COUPLED_ITER = 60
MAX_BACKTRACKS = 30


@njit(cache=True)
def _crack_N(theta):
    """``N`` de la spec §1: ``ε^cr = N·e`` (Voigt con γ de ingeniería) y
    ``[σ_nn, σ_tt, τ_nt] = Nᵀ·σ``; es la inversa de la ``T`` de :func:`_axes`."""
    c = math.cos(theta)
    s = math.sin(theta)
    N = np.empty((3, 3))
    N[0, 0] = c * c
    N[0, 1] = s * s
    N[0, 2] = -c * s
    N[1, 0] = s * s
    N[1, 1] = c * c
    N[1, 2] = c * s
    N[2, 0] = 2.0 * c * s
    N[2, 1] = -2.0 * c * s
    N[2, 2] = c * c - s * s
    return N


@njit(cache=True)
def _compression_state(sig_eff, E, nu, plane_strain, f_c, K, kappa_u, tol_c, cap, kc_n,
                       ktol):
    """Daño de compresión implícito del concreto entre grietas a partir de su
    esfuerzo efectivo (§4; el modelo 2 con ``d⁺ = 0``).

    ``κ⁻ = max(κ⁻_n, ε⁻_eq(σ̄⁻))`` con la medida de Drucker-Prager sobre la
    parte negativa del estado 3D asociado (en deformación plana
    ``σ̄_zz = ν·(σ̄_xx + σ̄_yy)`` cuenta si es negativo; ``H(0) = 0``). Devuelve
    ``(κ⁻, d⁻, dd⁻/dκ, σ̄⁻, Q⁺, g, σ̄_zz < 0)``, con ``dd⁻/dκ ≠ 0`` sólo en
    carga de compresión con ``0 < d⁻ < tope`` y ``g = ∂ε⁻_eq/∂σ̄⁻``.

    La carga es ``ε⁻_eq > κ⁻_n + ktol``. En el problema acoplado ``ktol =
    tol_c/E`` (la tolerancia de admisibilidad del ADR 0006, como la banda
    ``δ`` de las aperturas): al volver a evaluar el estado confirmado, el
    Newton acoplado reproduce ``ε⁻_eq`` sólo hasta el redondeo y, sin ella,
    ``κ⁻`` «crecía» y la tangente era la de carga en el 5–7 % de esos pasos
    (H1 de la verificación del bloque 2). Sin grieta ``ktol = 0``, como el
    modelo 2 (allí ``ε⁻_eq`` se recalcula idéntico).
    """
    sp, Q, _, _ = spectral_split_2d(sig_eff)
    sm = sig_eff - sp
    szz_m = 0.0
    szz_neg = False
    if plane_strain:
        szz = nu * (sig_eff[0] + sig_eff[1])
        if szz <= 0.0:
            szz_m = szz
            szz_neg = True
    eps_c, g = compression_measure(sm[0], sm[1], szz_m, sm[2], 0.0, 0.0, E, K)
    loading = eps_c > kc_n + ktol
    kc = eps_c if loading else kc_n
    d, dd = compression_damage(kc, E, f_c, kappa_u, tol_c, cap)
    if not (loading and d > 0.0 and d < cap):
        dd = 0.0
    return kc, d, dd, sm, Q, g, szz_neg


@njit(cache=True)
def _bulk(eps_co, Ce, P, plane_strain, implicit):
    """Concreto entre grietas (§4): ``σ_b = σ̄ − d⁻·σ̄⁻``, ``σ̄ = C_e·ε^co``.

    Devuelve ``(σ_b, C_b, κ⁻, d⁻)``. ``implicit``: ``d⁻`` evoluciona con
    ``ε^co`` (:func:`_compression_state`) y ``C_b`` es la tangente
    algorítmica del modelo 2 con ``d⁺ = 0``, ``C^d − σ̄⁻⊗∂d⁻/∂ε`` con
    ``C^d = C_e − d⁻·(I − Q⁺)·C_e`` y, en deformación plana, el término de
    ``σ̄_zz`` (§8b). Si no, ``d⁻ = d_fix`` y ``C_b = C^d`` (IMPL-EX, §9).
    """
    sig_eff = Ce @ eps_co
    E = P[P_E]
    nu = P[P_NU]
    if implicit:
        kc, d, dd, sm, Q, g, szz_neg = _compression_state(
            sig_eff, E, nu, plane_strain, P[P_FC], P[P_KDP], P[P_KU], P[P_TOLC], P[P_CAP],
            P[P_KCN], P[P_KTOL])
    else:
        sp, Q, _, _ = spectral_split_2d(sig_eff)
        sm = sig_eff - sp
        kc = P[P_KCN]
        d = P[P_DFIX]
        dd = 0.0
        g = np.zeros(6)
        szz_neg = False
    sigma = sig_eff - d * sm
    C = Ce - d * (Ce - Q @ Ce)
    if dd > 0.0:
        # −σ̄⁻ ⊗ ∂d⁻/∂ε^co, ∂d⁻/∂ε = d⁻'·(∂ε⁻_eq/∂σ̄⁻)·(I − Q⁺)·C_e (modelo 2,
        # §7): el gradiente respecto del vector Voigt lleva el factor 2 en el
        # cortante; en deformación plana, más ∂ε⁻_eq/∂σ̄_zz·ν·(C_e[0] + C_e[1]).
        gv = np.empty(3)
        gv[0] = g[0]
        gv[1] = g[1]
        gv[2] = 2.0 * g[3]
        v = (gv - gv @ Q) @ Ce
        if szz_neg:
            for b in range(3):
                v[b] += g[2] * nu * (Ce[0, b] + Ce[1, b])
        C = C - dd * np.outer(sm, v)
    return sigma, C, kc, d


@njit(cache=True)
def _dir_law(w, mode, wmax, c, delta, f_t, G_f, law, w_c, c_max):
    """``g_i(w)``, ``g_i'(w)`` de una dirección del problema acoplado (§8b) y
    si su historia crece (carga sobre la envolvente).

    El valor es continuo: secante ``w/c`` hasta ``w_max`` y la envolvente con
    el tope después; en una dirección sin grieta, ``f_t`` en ``w ≤ 0`` y la
    envolvente desde 0. La pendiente y la carga usan la banda ``w_max + δ``
    de la tolerancia del ADR 0006, ``δ = tol/(1/c + C11/h)``, la del bloque 1
    (``_solve_direction``): dentro de la banda la dirección descarga y la
    historia no crece.
    """
    if mode == SECANT_LAW:
        return w / c, 1.0 / c, False
    if wmax <= 0.0:
        if w <= 0.0:
            return f_t, 0.0, False
        F, dF, _ = _graph(w, ENVELOPE_LAW, 0.0, f_t, G_f, law, w_c, c_max)
        return F, dF, True
    if w <= wmax:
        return w / c, 1.0 / c, False
    F, dF, _ = _graph(w, ENVELOPE_LAW, 0.0, f_t, G_f, law, w_c, c_max)
    if w <= wmax + delta:
        return F, 1.0 / c, False
    return F, dF, True


@njit(cache=True)
def _coupled_eval(x, strain, N, Ce, P, I, Hd):
    """Residuo y jacobiano del problema local acoplado (§8b) en
    ``x = [w_n, w_t, s]``.

    ``ε^co = ε − N·x/h``, ``σ = σ_b(ε^co)`` y ``t = Nᵀ·σ``. Cada dirección
    normal, con la función mínimo ``r_i = min(k·w_i, g_i(w_i) − t_i)``,
    ``k = C11/h`` (Newton semiliso: la fila de la rama activa); una dirección
    rígida (IMPL-EX con ``c̃ = 0``), ``r_i = k·w_i``. Cortante, en unidades de
    esfuerzo (H2 de la revisión del bloque 2): ``r_s = (c_s·τ − s)/(c_s +
    h/G)``, con ``c_s = 2·(c_n + c_t)`` de la historia actualizada.

    Devuelve ``(r, J, σ, C_b, κ⁻, d⁻, Nᵀ·C_b, abierta, carga, c_s)``.
    """
    h = P[P_H]
    k = P[P_K]
    f_t = P[P_FT]
    G_f = P[P_GF]
    w_c = P[P_WC]
    c_max = P[P_CMAX]
    law = I[I_LAW]
    eps_co = strain - N @ (x / h)
    sig, Cb, kc, d = _bulk(eps_co, Ce, P, I[I_PS] == 1, I[I_IMPL] == 1)
    NtCb = N.T @ Cb
    tl = N.T @ sig
    A = NtCb @ N
    r = np.zeros(3)
    J = np.zeros((3, 3))
    opened = np.zeros(2, dtype=np.bool_)
    loading = np.zeros(2, dtype=np.bool_)
    c_new = np.zeros(2)
    dc = np.zeros(2)
    for i in range(2):
        mode = I[I_MODE_N + i]
        c_new[i] = Hd[HD_C + i]
        w = x[i]
        if mode == RIGID_LAW:
            r[i] = k * w
            J[i, i] = k
            continue
        g, gp, ld = _dir_law(w, mode, Hd[HD_WMAX + i], Hd[HD_C + i], Hd[HD_DELTA + i],
                             f_t, G_f, law, w_c, c_max)
        if k * w < g - tl[i]:
            r[i] = k * w
            J[i, i] = k
        else:
            opened[i] = True
            r[i] = g - tl[i]
            for j in range(3):
                J[i, j] = A[i, j] / h
            J[i, i] += gp
            if ld:
                loading[i] = True
                cc, dcc = _compliance(w, f_t, G_f, law, w_c, c_max)
                c_new[i] = cc
                dc[i] = dcc
    c_s = 2.0 * (c_new[0] + c_new[1])
    Ds = c_s + h / P[P_G]
    tau = tl[2]
    num = c_s * tau - x[2]
    r[2] = num / Ds
    for j in range(3):
        dcs = 2.0 * dc[j] if j < 2 else 0.0
        J[2, j] = (-c_s * A[2, j] / h + tau * dcs) / Ds - num * dcs / (Ds * Ds)
    J[2, 2] -= 1.0 / Ds
    return r, J, sig, Cb, kc, d, NtCb, opened, loading, c_s


@njit(cache=True)
def _solve3(J, B):
    """``X = J⁻¹·B`` (``B`` de 3×m) por eliminación de Gauss con pivote
    parcial. Devuelve ``(X, ok)``; ``ok = False`` con un pivote nulo o no
    finito."""
    A = J.copy()
    X = B.copy()
    m = X.shape[1]
    for col in range(3):
        p = col
        best = abs(A[col, col])
        for rr in range(col + 1, 3):
            if abs(A[rr, col]) > best:
                best = abs(A[rr, col])
                p = rr
        if not (best > 0.0 and math.isfinite(best)):
            return X, False
        if p != col:
            for j in range(3):
                tmp = A[col, j]
                A[col, j] = A[p, j]
                A[p, j] = tmp
            for j in range(m):
                tmp = X[col, j]
                X[col, j] = X[p, j]
                X[p, j] = tmp
        for rr in range(col + 1, 3):
            f = A[rr, col] / A[col, col]
            for j in range(col, 3):
                A[rr, j] -= f * A[col, j]
            for j in range(m):
                X[rr, j] -= f * X[col, j]
    for col in range(2, -1, -1):
        for j in range(m):
            acc = X[col, j]
            for rr in range(col + 1, 3):
                acc -= A[col, rr] * X[rr, j]
            X[col, j] = acc / A[col, col]
    return X, True


@njit(cache=True)
def _norm_inf(r):
    m = 0.0
    for v in r:
        a = abs(v)
        if not (a <= m):
            m = a
    return m


@njit(cache=True)
def _coupled_newton(x0, line_search, tol, strain, N, Ce, P, I, Hd):
    """Un intento del Newton acoplado desde ``x0``: paso completo o con
    búsqueda lineal (Armijo sobre ``‖r‖²``). Converge con ``max|r_i| ≤ tol``
    (esfuerzo) y después da un paso más de pulido, que se acepta si no
    aumenta el residuo: con la convergencia cuadrática deja la solución
    cerca del redondeo, como el algoritmo anidado del bloque 1, para que la
    tangente consistente y las diferencias finitas no vean la tolerancia.
    Devuelve ``(x, ok)``."""
    x = x0.copy()
    out = _coupled_eval(x, strain, N, Ce, P, I, Hd)
    r = out[0]
    J = out[1]
    nr = _norm_inf(r)
    for _ in range(MAX_COUPLED_ITER):
        if not math.isfinite(nr):
            return x, False
        B = np.empty((3, 1))
        for i in range(3):
            B[i, 0] = -r[i]
        X, ok = _solve3(J, B)
        if nr <= tol:
            if ok and nr > 0.0:
                xp = x + X[:, 0]
                rp = _coupled_eval(xp, strain, N, Ce, P, I, Hd)[0]
                if _norm_inf(rp) <= nr:
                    return xp, True
            return x, True
        if not ok:
            return x, False
        dx = X[:, 0]
        if not line_search:
            x = x + dx
            out = _coupled_eval(x, strain, N, Ce, P, I, Hd)
        else:
            phi0 = r @ r
            alpha = 1.0
            xt = x + dx
            out = _coupled_eval(xt, strain, N, Ce, P, I, Hd)
            for _ in range(MAX_BACKTRACKS):
                rt = out[0]
                phi = rt @ rt
                if math.isfinite(phi) and phi <= (1.0 - 1.0e-4 * alpha) * phi0:
                    break
                alpha *= 0.5
                xt = x + alpha * dx
                out = _coupled_eval(xt, strain, N, Ce, P, I, Hd)
            x = xt
        r = out[0]
        J = out[1]
        nr = _norm_inf(r)
    return x, nr <= tol


@njit(cache=True)
def _coupled_solve(x0, tol, strain, N, Ce, P, I, Hd):
    """Cascada de §8b (H1 de la revisión del bloque 2): Newton de paso completo
    desde ``x0`` (el paso confirmado), con búsqueda lineal desde ``x0`` y,
    como último recurso, los dos desde ``x = 0``. Devuelve ``(x, ok)``."""
    x, ok = _coupled_newton(x0, False, tol, strain, N, Ce, P, I, Hd)
    if ok:
        return x, True
    x, ok = _coupled_newton(x0, True, tol, strain, N, Ce, P, I, Hd)
    if ok:
        return x, True
    zero = np.zeros(3)
    x, ok = _coupled_newton(zero, False, tol, strain, N, Ce, P, I, Hd)
    if ok:
        return x, True
    return _coupled_newton(zero, True, tol, strain, N, Ce, P, I, Hd)


@njit(cache=True)
def _coupled_finish(x, strain, N, Ce, P, I, Hd):
    """Esfuerzo y tangente consistente en la solución ``x`` (§8b).

    Por el teorema de la función implícita sobre ``r(x, ε) = 0``:
    ``dx = −J⁻¹·R_ε·dε`` y ``dσ = C_b·(dε − N·dx/h)``, es decir
    ``C_tan = C_b + C_b·N·J⁻¹·R_ε/h``, con ``R_ε`` la derivada del residuo
    respecto de ``ε``: ``−(Nᵀ·C_b)_i`` en una dirección abierta, 0 en una
    cerrada y ``c_s/(c_s + h/G)·(Nᵀ·C_b)_s`` en el cortante. Es la fórmula de
    §8b con la fila de cortante en flexibilidad, válida con ``C_b`` no
    simétrica. Devuelve ``(σ, C_tan, κ⁻, d⁻, carga, ok)``.
    """
    h = P[P_H]
    r, J, sig, Cb, kc, d, NtCb, opened, loading, c_s = _coupled_eval(
        x, strain, N, Ce, P, I, Hd)
    R = np.zeros((3, 3))
    for i in range(2):
        if opened[i]:
            for j in range(3):
                R[i, j] = -NtCb[i, j]
    fs = c_s / (c_s + h / P[P_G])
    for j in range(3):
        R[2, j] = fs * NtCb[2, j]
    X, ok = _solve3(J, R)
    C_tan = Cb + (Cb @ N) @ X / h
    return sig, C_tan, kc, d, loading, ok


@njit(cache=True)
def _coupled_arrays(Ce, E, nu, h, f_t, G_f, w_c, c_max, f_c, K, kappa_u, tol_c, cap, kc_n,
                    d_fix, plane_strain, law, implicit, mode_n, mode_t, wmax_n, wmax_t,
                    c_n, c_t, tol):
    """Arreglos ``P``, ``I`` y ``Hd`` del problema acoplado, con la banda
    ``δ_i = tol/(1/c_i + C11/h)`` del quiebre ``w = w_max`` (la del bloque 1)
    y, con grieta (``mode_n`` no rígido), la tolerancia ``tol_c/E`` de la
    carga de ``κ⁻`` (:func:`_compression_state`)."""
    P = np.empty(17)
    P[P_E] = E
    P[P_NU] = nu
    P[P_H] = h
    P[P_K] = Ce[0, 0] / h
    P[P_G] = Ce[2, 2]
    P[P_FT] = f_t
    P[P_GF] = G_f
    P[P_WC] = w_c
    P[P_CMAX] = c_max
    P[P_FC] = f_c
    P[P_KDP] = K
    P[P_KU] = kappa_u
    P[P_TOLC] = tol_c
    P[P_CAP] = cap
    P[P_KCN] = kc_n
    P[P_DFIX] = d_fix
    P[P_KTOL] = 0.0 if mode_n == RIGID_LAW and mode_t == RIGID_LAW else tol_c / E
    I = np.empty(5, dtype=np.int64)
    I[I_PS] = 1 if plane_strain else 0
    I[I_LAW] = law
    I[I_IMPL] = 1 if implicit else 0
    I[I_MODE_N] = mode_n
    I[I_MODE_T] = mode_t
    Hd = np.zeros(6)
    Hd[HD_WMAX] = wmax_n
    Hd[HD_WMAX + 1] = wmax_t
    Hd[HD_C] = c_n
    Hd[HD_C + 1] = c_t
    k = Ce[0, 0] / h
    if c_n > 0.0:
        Hd[HD_DELTA] = tol / (1.0 / c_n + k)
    if c_t > 0.0:
        Hd[HD_DELTA + 1] = tol / (1.0 / c_t + k)
    return P, I, Hd


@njit(cache=True)
def _principal_max(sig):
    return 0.5 * (sig[0] + sig[1]) + math.hypot(0.5 * (sig[0] - sig[1]), sig[2])


# ---------------------------------------------------------------------------
# Integración implícita (§8) e IMPL-EX (§9)
# ---------------------------------------------------------------------------

@njit(cache=True)
def _implicit_point(strain, S_old, Ce, h, f_t, G_f, law, w_c, c_max, tol):
    """Paso implícito con la historia confirmada ``S_old``.

    Devuelve ``(σ, C_tan, agrietado, θ, w_n, w_t, s, w_max_n, w_max_t,
    falta_h)``. Sin grieta, el esfuerzo de prueba es elástico; si su principal
    mayor supera ``f_t`` (con la tolerancia de admisibilidad), la grieta se
    forma con la dirección principal de ese esfuerzo, ``θ = ½·atan2(2σ_xy,
    σ_xx − σ_yy)``, y el paso se resuelve con la grieta. ``θ`` sólo queda en
    el estado de prueba: un intento rechazado no la fija.
    """
    cracked = S_old[CRACKED] > 0.0
    theta = S_old[THETA]
    wmax_n = S_old[WMAX_N]
    wmax_t = S_old[WMAX_T]
    if not cracked:
        sig = Ce @ strain
        sx = sig[0]
        sy = sig[1]
        sxy = sig[2]
        s1 = 0.5 * (sx + sy) + math.hypot(0.5 * (sx - sy), sxy)
        if s1 - f_t <= tol:
            return sig, Ce.copy(), 0.0, theta, 0.0, 0.0, 0.0, wmax_n, wmax_t, False
        if h <= 0.0:
            return sig, Ce.copy(), 0.0, theta, 0.0, 0.0, 0.0, wmax_n, wmax_t, True
        theta = 0.5 * math.atan2(2.0 * sxy, sx - sy)
        wmax_n = 0.0
        wmax_t = 0.0
    elif h <= 0.0:
        return Ce @ strain, Ce.copy(), 1.0, theta, 0.0, 0.0, 0.0, wmax_n, wmax_t, True
    T = _axes(theta)
    c_n, _ = _compliance(wmax_n, f_t, G_f, law, w_c, c_max)
    c_t, _ = _compliance(wmax_t, f_t, G_f, law, w_c, c_max)
    sig_loc, C_loc, w_n, w_t, slip, wn_new, wt_new = _crack_response(
        T @ strain, Ce[0, 0], Ce[0, 1], Ce[2, 2], h, f_t, G_f, law, w_c, c_max, tol,
        IMPLICIT_LAW, wmax_n, c_n, IMPLICIT_LAW, wmax_t, c_t, True)
    return (T.T @ sig_loc, T.T @ C_loc @ T, 1.0, theta, w_n, w_t, slip,
            wn_new, wt_new, False)


@njit(cache=True)
def _secant_point(strain, theta, c_n, c_t, Ce, h, f_t, G_f, law, w_c, c_max, tol):
    """Ley secante con las flexibilidades fijas ``c_n``, ``c_t`` (IMPL-EX, §9).

    El problema normal tiene a lo sumo cuatro combinaciones abierta/cerrada y
    una sola consistente (convexidad de χ, §6); la dirección queda cerrada con
    ``σ_ii ≤ tol`` (el empate ``σ_ii = 0`` exacto incluido; ADR 0006).
    Tangente secante, simétrica (se simetriza para quitar el redondeo).
    Devuelve ``(σ, C, w_n, w_t, s)``.
    """
    T = _axes(theta)
    mode_n = SECANT_LAW if c_n > 0.0 else RIGID_LAW
    mode_t = SECANT_LAW if c_t > 0.0 else RIGID_LAW
    sig_loc, C_loc, w_n, w_t, slip, _, _ = _crack_response(
        T @ strain, Ce[0, 0], Ce[0, 1], Ce[2, 2], h, f_t, G_f, law, w_c, c_max, tol,
        mode_n, 0.0, c_n, mode_t, 0.0, c_t, False)
    C = T.T @ C_loc @ T
    return T.T @ sig_loc, 0.5 * (C + C.T), w_n, w_t, slip


@njit(cache=True)
def _nan_result():
    sig = np.full(3, np.nan)
    C = np.full((3, 3), np.nan)
    return sig, C


@njit(cache=True)
def _implicit_point_crush(strain, S_old, Ce, E, nu, plane_strain, h, f_t, G_f, law, w_c,
                          c_max, tol, f_c, K, kappa_u, tol_c, cap):
    """Paso implícito con aplastamiento (§8b).

    Devuelve ``(σ, C_tan, agrietado, θ, w_n, w_t, s, w_max_n, w_max_t, κ⁻,
    d⁻, falta_h, ok)``; ``ok = False`` si la cascada del Newton acoplado no
    converge (el llamador devuelve un esfuerzo no finito y el solver repite
    el paso más corto, Q5).

    1. Predictor: con ``d⁻_n = 0`` se resuelve el bloque 1
       (:func:`_implicit_point`); su esfuerzo es el efectivo del concreto
       entre grietas (elástico), y si con él ``d⁻(κ⁻) = 0`` es la solución.
    2. Si no, sin grieta: el esfuerzo de prueba es ``σ_b(ε)``; si su
       principal mayor supera ``f_t`` la grieta se forma con la dirección de
       ``σ̄ = C_e·ε`` (la misma, §2) y se resuelve acoplado desde ``x = 0``.
       Con grieta: el problema acoplado desde las deformaciones de la grieta
       del paso confirmado.
    """
    kc_n = S_old[KC]
    d_n, _ = compression_damage(kc_n, E, f_c, kappa_u, tol_c, cap)
    if d_n == 0.0:
        (sig, C, cracked, theta, w_n, w_t, slip, wn_new, wt_new,
         missing) = _implicit_point(strain, S_old, Ce, h, f_t, G_f, law, w_c, c_max, tol)
        if missing:
            return (sig, C, cracked, theta, w_n, w_t, slip, wn_new, wt_new, kc_n, 0.0,
                    True, True)
        kc, d, _, _, _, _, _ = _compression_state(sig, E, nu, plane_strain, f_c, K, kappa_u,
                                                  tol_c, cap, kc_n, 0.0)
        if d == 0.0:
            return (sig, C, cracked, theta, w_n, w_t, slip, wn_new, wt_new, kc, 0.0,
                    False, True)
    cracked = S_old[CRACKED] > 0.0
    theta = S_old[THETA]
    wmax_n = S_old[WMAX_N]
    wmax_t = S_old[WMAX_T]
    x0 = np.zeros(3)
    if not cracked:
        Pb, Ib, _ = _coupled_arrays(Ce, E, nu, 1.0, f_t, G_f, w_c, c_max, f_c, K, kappa_u,
                                    tol_c, cap, kc_n, 0.0, plane_strain, law, True,
                                    RIGID_LAW, RIGID_LAW, 0.0, 0.0, 0.0, 0.0, tol)
        sig, C, kc, d = _bulk(strain, Ce, Pb, plane_strain, True)
        if _principal_max(sig) - f_t <= tol:
            return (sig, C, 0.0, theta, 0.0, 0.0, 0.0, wmax_n, wmax_t, kc, d,
                    h <= 0.0 and d > 0.0, True)
        if h <= 0.0:
            return (sig, C, 0.0, theta, 0.0, 0.0, 0.0, wmax_n, wmax_t, kc, d, True, True)
        sb = Ce @ strain
        theta = 0.5 * math.atan2(2.0 * sb[2], sb[0] - sb[1])
        wmax_n = 0.0
        wmax_t = 0.0
    else:
        if h <= 0.0:
            sig, C = _nan_result()
            return (sig, C, 1.0, theta, 0.0, 0.0, 0.0, wmax_n, wmax_t, kc_n, d_n, True,
                    True)
        x0[0] = S_old[W_N]
        x0[1] = S_old[W_T]
        x0[2] = S_old[SLIP]
    N = _crack_N(theta)
    c_n, _ = _compliance(wmax_n, f_t, G_f, law, w_c, c_max)
    c_t, _ = _compliance(wmax_t, f_t, G_f, law, w_c, c_max)
    P, I, Hd = _coupled_arrays(Ce, E, nu, h, f_t, G_f, w_c, c_max, f_c, K, kappa_u, tol_c,
                               cap, kc_n, 0.0, plane_strain, law, True, IMPLICIT_LAW,
                               IMPLICIT_LAW, wmax_n, wmax_t, c_n, c_t, tol)
    x, ok = _coupled_solve(x0, tol, strain, N, Ce, P, I, Hd)
    sig, C = _nan_result()
    kc = kc_n
    d = d_n
    loading = np.zeros(2, dtype=np.bool_)
    if ok:
        sig, C, kc, d, loading, ok = _coupled_finish(x, strain, N, Ce, P, I, Hd)
    if not ok:
        sig, C = _nan_result()
        return (sig, C, 1.0, theta, 0.0, 0.0, 0.0, wmax_n, wmax_t, kc_n, d_n, False, False)
    wn_new = x[0] if loading[0] else wmax_n
    wt_new = x[1] if loading[1] else wmax_t
    return (sig, C, 1.0, theta, max(x[0], 0.0), max(x[1], 0.0), x[2], wn_new, wt_new, kc,
            d, False, True)


@njit(cache=True)
def _frozen_point_crush(strain, cracked, theta, c_n, c_t, d_fix, x0, Ce, E, nu,
                        plane_strain, h, f_t, G_f, law, w_c, c_max, tol):
    """Esfuerzo con la historia fija (IMPL-EX, §9): flexibilidades ``c_n``,
    ``c_t`` y daño de compresión ``d_fix``.

    Con ``d_fix = 0``, la ley secante del bloque 1 (:func:`_secant_point`, en
    forma cerrada y simétrica). Con ``d_fix > 0`` no hay forma cerrada: el
    Newton acoplado de §8b con el concreto a daño fijo y la secante de cada
    grieta; tangente consistente, no simétrica. Sin grieta, el concreto
    solo. Devuelve ``(σ, C, w_n, w_t, s, ok)``.
    """
    if not cracked:
        if d_fix == 0.0:
            return Ce @ strain, Ce.copy(), 0.0, 0.0, 0.0, True
        Pb, _, _ = _coupled_arrays(Ce, E, nu, 1.0, f_t, G_f, w_c, c_max, 0.0, 0.0, 0.0,
                                   0.0, 0.0, 0.0, d_fix, plane_strain, law, False,
                                   RIGID_LAW, RIGID_LAW, 0.0, 0.0, 0.0, 0.0, tol)
        sig, C, _, _ = _bulk(strain, Ce, Pb, plane_strain, False)
        return sig, C, 0.0, 0.0, 0.0, True
    if d_fix == 0.0:
        sig, C, w_n, w_t, slip = _secant_point(strain, theta, c_n, c_t, Ce, h, f_t, G_f, law,
                                               w_c, c_max, tol)
        return sig, C, w_n, w_t, slip, True
    mode_n = SECANT_LAW if c_n > 0.0 else RIGID_LAW
    mode_t = SECANT_LAW if c_t > 0.0 else RIGID_LAW
    N = _crack_N(theta)
    P, I, Hd = _coupled_arrays(Ce, E, nu, h, f_t, G_f, w_c, c_max, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, d_fix, plane_strain, law, False, mode_n, mode_t, 0.0,
                               0.0, c_n, c_t, tol)
    x, ok = _coupled_solve(x0, tol, strain, N, Ce, P, I, Hd)
    sig, C = _nan_result()
    if ok:
        sig, C, _, _, _, ok = _coupled_finish(x, strain, N, Ce, P, I, Hd)
    if not ok:
        sig, C = _nan_result()
        return sig, C, 0.0, 0.0, 0.0, False
    return sig, C, max(x[0], 0.0), max(x[1], 0.0), x[2], True


@njit(cache=True)
def _fsc_crush_point(strain, S_old, Ce, E, f_t, G_f, law, w_c, tol, cap, integration, ratio,
                     f_c, G_c, K, nu, plane_strain, tol_c):
    """Un punto con aplastamiento (bloque 2), implícito o IMPL-EX (§8b, §9).

    Como :func:`_fsc_b1_point`, más ``κ⁻`` (``kappa_c``), ``d⁻``
    (``damage_c``: el de la historia nueva en la implícita, el ``d̃⁻`` del
    esfuerzo en IMPL-EX) y ``kappa_c_prev``. IMPL-EX:
    ``κ̃⁻ = κ⁻_n + r·(κ⁻_n − κ⁻_{n−1})`` y ``d̃⁻ = max(d̃⁻_n, d⁻(κ̃⁻))``
    (como el modelo 2); el indicador de error incluye ``|d⁻(κ⁻) − d̃⁻|``.
    Devuelve ``(σ, C_tan, fila, falta_h)``; si el problema local no converge,
    ``σ`` y ``C_tan`` no finitos (Q5).
    """
    row = S_old.copy()
    h = S_old[H]
    c_max = _cap_compliance(h, E, cap)
    kappa_u = compression_ultimate(h, E, f_c, G_c)
    x_old = np.empty(3)
    x_old[0] = S_old[W_N]
    x_old[1] = S_old[W_T]
    x_old[2] = S_old[SLIP]
    if integration != IMPL_EX:
        (sig, C, cracked, theta, w_n, w_t, slip, wn_new, wt_new, kc, d, missing,
         ok) = _implicit_point_crush(strain, S_old, Ce, E, nu, plane_strain, h, f_t, G_f,
                                     law, w_c, c_max, tol, f_c, K, kappa_u, tol_c, cap)
        if not ok:
            return sig, C, row, False
        row[CRACKED] = cracked
        row[THETA] = theta
        row[WMAX_N] = wn_new
        row[WMAX_T] = wt_new
        row[W_N] = w_n
        row[W_T] = w_t
        row[SLIP] = slip
        if h > 0.0:
            row[DMG_N] = _omega(_compliance(wn_new, f_t, G_f, law, w_c, c_max)[0], h, E)
            row[DMG_T] = _omega(_compliance(wt_new, f_t, G_f, law, w_c, c_max)[0], h, E)
        row[KC] = kc
        row[DC] = d
        row[WMAX_N_PREV] = S_old[WMAX_N]
        row[WMAX_T_PREV] = S_old[WMAX_T]
        row[KC_PREV] = S_old[KC]
        row[ERR] = 0.0
        return sig, C, row, missing
    cracked_old = S_old[CRACKED] > 0.0
    active = cracked_old and h > 0.0
    if ratio < 0.0:
        c_n = _compliance_from_omega(S_old[DMG_N], h, E) if active else 0.0
        c_t = _compliance_from_omega(S_old[DMG_T], h, E) if active else 0.0
        sig, C, _, _, _, _ = _frozen_point_crush(strain, active, S_old[THETA], c_n, c_t,
                                                 S_old[DC], x_old, Ce, E, nu, plane_strain,
                                                 h, f_t, G_f, law, w_c, c_max, tol)
        return sig, C, row, cracked_old and h <= 0.0
    kc_n = S_old[KC]
    kappa_tilde = kc_n + ratio * (kc_n - S_old[KC_PREV])
    d_tilde, _ = compression_damage(kappa_tilde, E, f_c, kappa_u, tol_c, cap)
    if d_tilde < S_old[DC]:
        d_tilde = S_old[DC]
    om_n = S_old[DMG_N]
    om_t = S_old[DMG_T]
    c_n = 0.0
    c_t = 0.0
    if active:
        wt_n = S_old[WMAX_N] + ratio * (S_old[WMAX_N] - S_old[WMAX_N_PREV])
        wt_t = S_old[WMAX_T] + ratio * (S_old[WMAX_T] - S_old[WMAX_T_PREV])
        om_n = max(om_n, _omega(_compliance(wt_n, f_t, G_f, law, w_c, c_max)[0], h, E))
        om_t = max(om_t, _omega(_compliance(wt_t, f_t, G_f, law, w_c, c_max)[0], h, E))
        c_n = _compliance_from_omega(om_n, h, E)
        c_t = _compliance_from_omega(om_t, h, E)
    sig, C, w_n, w_t, slip, ok = _frozen_point_crush(strain, active, S_old[THETA], c_n, c_t,
                                                     d_tilde, x_old, Ce, E, nu, plane_strain,
                                                     h, f_t, G_f, law, w_c, c_max, tol)
    if not ok:
        return sig, C, row, False
    (_, _, cracked, theta, _, _, _, wn_new, wt_new, kc, d_impl, missing,
     ok) = _implicit_point_crush(strain, S_old, Ce, E, nu, plane_strain, h, f_t, G_f, law,
                                 w_c, c_max, tol, f_c, K, kappa_u, tol_c, cap)
    if not ok:
        sig, C = _nan_result()
        return sig, C, row, False
    row[CRACKED] = cracked
    row[THETA] = theta
    row[WMAX_N] = wn_new
    row[WMAX_T] = wt_new
    row[WMAX_N_PREV] = S_old[WMAX_N]
    row[WMAX_T_PREV] = S_old[WMAX_T]
    row[W_N] = w_n
    row[W_T] = w_t
    row[SLIP] = slip
    row[KC] = kc
    row[KC_PREV] = kc_n
    row[DC] = d_tilde
    err = abs(d_impl - d_tilde)
    if h > 0.0:
        row[DMG_N] = om_n
        row[DMG_T] = om_t
        err_n = abs(_omega(_compliance(wn_new, f_t, G_f, law, w_c, c_max)[0], h, E) - om_n)
        err_t = abs(_omega(_compliance(wt_new, f_t, G_f, law, w_c, c_max)[0], h, E) - om_t)
        err = max(err, err_n, err_t)
    row[ERR] = err
    return sig, C, row, missing


@njit(cache=True)
def _fsc_2d_point(strain, S_old, Ce, E, f_t, G_f, law, w_c, tol, cap, integration, ratio,
                  f_c, G_c, K, nu, plane_strain, tol_c):
    """Un punto de integración: sin aplastamiento (``f_c = 0``) el bloque 1
    tal cual (:func:`_fsc_b1_point`); con él, :func:`_fsc_crush_point`."""
    if f_c > 0.0:
        return _fsc_crush_point(strain, S_old, Ce, E, f_t, G_f, law, w_c, tol, cap,
                                integration, ratio, f_c, G_c, K, nu, plane_strain, tol_c)
    return _fsc_b1_point(strain, S_old, Ce, E, f_t, G_f, law, w_c, tol, cap, integration,
                         ratio)


@njit(cache=True)
def _fsc_b1_point(strain, S_old, Ce, E, f_t, G_f, law, w_c, tol, cap, integration, ratio):
    """Un punto de integración sin aplastamiento, implícito o IMPL-EX (spec §6–§9).

    ``Ce`` es la elasticidad plana de la hipótesis (esfuerzo plano o
    deformación plana). Devuelve ``(σ, C_tan, fila de estado nueva,
    falta_h)``; ``falta_h`` marca una grieta en un punto sin longitud
    característica sembrada.

    - Implícita: :func:`_implicit_point`. ``damage_n``, ``damage_t`` son los
      ``ω_i`` de la historia nueva.
    - IMPL-EX: ``ω̃_i = max(ω̃_i del paso anterior, ω(c(w̃_max_i)))`` con
      ``w̃_max_i = w_max_i + r·(w_max_i − w_max_i_prev)`` y ``c̃_i = c(ω̃_i)``
      (la flexibilidad es creciente en ω; la salvaguarda sobre ω la hace
      monótona también en el redondeo); el esfuerzo con la ley secante
      (:func:`_secant_point`) y la historia con la implícita. ``damage_n``,
      ``damage_t`` guardan los ``ω̃_i`` del esfuerzo y ``damage_error`` el
      indicador ``max_i |ω_i(w_max_i) − ω̃_i|``. Fuera de un paso (``r < 0``)
      el estado queda congelado y el esfuerzo es el de las ``c̃_i`` guardadas.
      Una grieta nueva entra en la extrapolación desde el paso siguiente.
    """
    row = S_old.copy()
    h = S_old[H]
    c_max = _cap_compliance(h, E, cap)
    if integration != IMPL_EX:
        sig, C, cracked, theta, w_n, w_t, slip, wn_new, wt_new, missing = _implicit_point(
            strain, S_old, Ce, h, f_t, G_f, law, w_c, c_max, tol)
        row[CRACKED] = cracked
        row[THETA] = theta
        row[WMAX_N] = wn_new
        row[WMAX_T] = wt_new
        row[W_N] = w_n
        row[W_T] = w_t
        row[SLIP] = slip
        if h > 0.0:
            row[DMG_N] = _omega(_compliance(wn_new, f_t, G_f, law, w_c, c_max)[0], h, E)
            row[DMG_T] = _omega(_compliance(wt_new, f_t, G_f, law, w_c, c_max)[0], h, E)
        row[WMAX_N_PREV] = S_old[WMAX_N]
        row[WMAX_T_PREV] = S_old[WMAX_T]
        row[ERR] = 0.0
        return sig, C, row, missing
    active = S_old[CRACKED] > 0.0 and h > 0.0
    if ratio < 0.0:
        if active:
            c_n = _compliance_from_omega(S_old[DMG_N], h, E)
            c_t = _compliance_from_omega(S_old[DMG_T], h, E)
            sig, C, _, _, _ = _secant_point(strain, S_old[THETA], c_n, c_t, Ce, h, f_t, G_f,
                                            law, w_c, c_max, tol)
            return sig, C, row, False
        return Ce @ strain, Ce.copy(), row, S_old[CRACKED] > 0.0
    w_n = 0.0
    w_t = 0.0
    slip = 0.0
    om_n = S_old[DMG_N]
    om_t = S_old[DMG_T]
    if active:
        # La salvaguarda se aplica sobre ω (guardada): ω̃ no decrece ni en el
        # redondeo de la ida y vuelta entre ω y c.
        wt_n = S_old[WMAX_N] + ratio * (S_old[WMAX_N] - S_old[WMAX_N_PREV])
        wt_t = S_old[WMAX_T] + ratio * (S_old[WMAX_T] - S_old[WMAX_T_PREV])
        om_n = max(om_n, _omega(_compliance(wt_n, f_t, G_f, law, w_c, c_max)[0], h, E))
        om_t = max(om_t, _omega(_compliance(wt_t, f_t, G_f, law, w_c, c_max)[0], h, E))
        c_n = _compliance_from_omega(om_n, h, E)
        c_t = _compliance_from_omega(om_t, h, E)
        sig, C, w_n, w_t, slip = _secant_point(strain, S_old[THETA], c_n, c_t, Ce, h, f_t,
                                               G_f, law, w_c, c_max, tol)
    else:
        sig = Ce @ strain
        C = Ce.copy()
    _, _, cracked, theta, _, _, _, wn_new, wt_new, missing = _implicit_point(
        strain, S_old, Ce, h, f_t, G_f, law, w_c, c_max, tol)
    row[CRACKED] = cracked
    row[THETA] = theta
    row[WMAX_N] = wn_new
    row[WMAX_T] = wt_new
    row[WMAX_N_PREV] = S_old[WMAX_N]
    row[WMAX_T_PREV] = S_old[WMAX_T]
    row[W_N] = w_n
    row[W_T] = w_t
    row[SLIP] = slip
    if h > 0.0:
        row[DMG_N] = om_n
        row[DMG_T] = om_t
        err_n = abs(_omega(_compliance(wn_new, f_t, G_f, law, w_c, c_max)[0], h, E) - om_n)
        err_t = abs(_omega(_compliance(wt_new, f_t, G_f, law, w_c, c_max)[0], h, E) - om_t)
        row[ERR] = max(err_n, err_t)
    return sig, C, row, missing


@njit(cache=True)
def _fsc_2d_batch(strain, S_old, S_new, params, C, sigma, C_out, flag):
    """Adaptador por lotes (ADR 0014). ``params = [E, f_t, G_f, tol_abs,
    tol_rel, damage_max, integración, r, ley de ablandamiento, w_c, f_c, G_c,
    K, ν, deformación plana]`` (``f_c = 0`` sin aplastamiento); ``C`` =
    ``C_e`` plana. Grieta o aplastamiento sin ``h`` sembrada → ``flag``
    (``batch_report``). Si el problema local no converge, ``σ`` y la
    tangente salen no finitos y el solver repite el paso (§8b, Q5)."""
    E = params[0]
    f_t = params[1]
    tol = params[3] + params[4] * f_t
    f_c = params[10]
    tol_c = params[3] + params[4] * (f_c / 3.0)
    sig, C_tan, row, missing = _fsc_2d_point(strain, S_old, C, E, f_t, params[2],
                                             int(params[8]), params[9], tol, params[5],
                                             int(params[6]), params[7], f_c, params[11],
                                             params[12], params[13], params[14] > 0.0,
                                             tol_c)
    if missing:
        flag[0] = 1
    for i in range(3):
        sigma[i] = sig[i]
    for i in range(16):
        S_new[i] = row[i]
    return C_tan


@MaterialRegistry.register
class FixedSmearedCrack2D(CrackBandMixin, Material):
    """
    Grieta distribuida fija con retención de cortante, banda de grieta e IMPL-EX.

    Parameters
    ----------
    E, nu : float
        Parámetros elásticos del concreto entre grietas, ``0 ≤ ν < 0.5``.
    f_t : float
        Resistencia a tracción: la grieta se forma cuando el esfuerzo principal
        mayor la alcanza.
    G_f : float
        Energía de fractura (fuerza/longitud); la grieta disipa ``G_f`` por
        unidad de área en tracción (``G_f/h`` por unidad de volumen).
    hypothesis : str
        ``'plane_stress'`` (default) o ``'plane_strain'``.
    softening : str
        Ley de ablandamiento en la apertura: ``'exponential'`` (default,
        ``f_t·exp(−f_t·w/G_f)``) o ``'hordijk'`` (``w_c = 5.136·G_f/f_t``).
    damage_max : float
        Tope del daño uniaxial equivalente de cada banda, en (0, 1); por
        omisión ``1 − 10⁻⁶``: la grieta conserva una rigidez residual.
    integration : str
        ``'implicit'`` (default) o ``'impl-ex'`` (Oliver, Huespe y Cante
        2008; ADR 0034): flexibilidades extrapoladas de los dos últimos pasos
        y tangente secante, simétrica; requiere un solver que declare su
        pseudo-tiempo.
    density : float, optional
        Densidad (ADR 0008).
    f_c, G_c : float, optional
        Resistencia a compresión uniaxial (positiva, ``f_c > f_t``) y energía
        de aplastamiento. Van juntas y activan el aplastamiento del concreto
        entre grietas (spec §4): el daño de compresión del modelo 2
        (Drucker-Prager sobre la parte negativa del esfuerzo efectivo y
        parábola de Feenstra con ``G_c/h``). Sin ellas la compresión es
        elástica.
    beta : float, optional
        ``f_bc/f_c`` de Drucker-Prager, ``β ≥ 1``; 1.16 por omisión (Kupfer).

    Notes
    -----
    ``ε = C_e⁻¹·σ + N·e`` con ``e = [e_nn, e_tt, γ_nt]`` las deformaciones de
    la grieta en sus ejes y ``N`` la matriz de la spec §1. La grieta queda fija
    en la dirección principal mayor del esfuerzo de prueba al agrietarse
    (``theta``); la segunda grieta es la perpendicular. Ley normal
    ``σ_nn = f(w)`` con ``w = h·e_nn``, descarga secante y cierre con
    ``w = 0``; ley de cortante ``τ = f(w_max)/(2·w_max)·s`` con ``s = h·γ_nt``
    y las dos grietas en serie. Tracción uniaxial: ``σ = f(w)`` y
    ``ε = σ/E + w/h`` con la contracción lateral sólo del concreto entre
    grietas; la energía es ``G_f/h`` exacta en esfuerzo plano y en
    deformación plana. Límite de la longitud característica (unicidad del
    problema local con las dos grietas en carga, spec §2 y §8):
    ``h < ℓ_ch/(1 + ν)`` (exponencial) o ``0.7382·ℓ_ch/(1 + ν)`` (Hordijk),
    ``ℓ_ch = E·G_f/f_t²``. A historia fija la ley deriva de un potencial y la
    descarga es secante. La tangente implícita no es simétrica en carga con
    cortante; la de IMPL-EX sí, salvo con aplastamiento. Con ``f_c`` y
    ``G_c``, la compresión paralela a una grieta abierta llega a ``f_c`` con
    cualquier apertura (sin reducción por agrietamiento transversal) y el
    problema local con concreto muy aplastado no tiene unicidad garantizada:
    se elige la solución por continuidad (spec §4, §8b). Estado: ``cracked``,
    ``theta``, ``w_max_n``, ``w_max_t`` (aperturas máximas), ``w_n``, ``w_t``,
    ``slip`` (post-proceso), ``damage_n``, ``damage_t`` (``ω`` de cada
    banda), ``kappa_c``, ``damage_c`` (``κ⁻`` y ``d⁻`` del concreto entre
    grietas; ``κ⁻`` empieza en ``f_c/(3E)``), ``h``, ``w_max_n_prev``,
    ``w_max_t_prev``, ``kappa_c_prev`` y ``damage_error`` (IMPL-EX).
    """
    # Con la historia fija la ley deriva de un potencial homogéneo de grado 2
    # (spec §6): σ(s·ε) = s·σ(ε) y ½·σ:ε es la energía recuperable (ADR 0038).
    SECANT_UNLOADING = True
    STRAIN_DIM = 3
    PRIMARY_STATE_VAR = 'damage_n'
    IS_SYMMETRIC = False
    STATE_SCHEMA = STATE_SCHEMA_FSC

    def __init__(self, E: float, nu: float, f_t: float, G_f: float,
                 hypothesis: str = 'plane_stress',
                 softening: str = DEFAULT_SOFTENING,
                 damage_max: float = CRACK_BAND_DAMAGE_MAX,
                 integration: str = DEFAULT_INTEGRATION,
                 density: float | None = None,
                 f_c: float | None = None, G_c: float | None = None,
                 beta: float | None = None):
        owner = type(self).__name__
        validate_crack_band_parameters(owner, E, nu, f_t, G_f, damage_max, density)
        if float(nu) < 0.0:
            raise ValueError(
                f"{owner}: nu = {nu} < 0 no se admite: el límite de la longitud "
                f"característica que garantiza una sola solución del problema local "
                f"(spec §8) se dedujo para 0 ≤ ν < 0.5.")
        self._integration_code = parse_integration(owner, integration)
        self._ratio = NO_STEP_RATIO
        self.E = float(E)
        self.nu = float(nu)
        self.f_t = float(f_t)
        self.G_f = float(G_f)
        self.kappa_0 = self.f_t / self.E
        self._set_softening(owner, softening)
        self.damage_max = float(damage_max)
        # Aplastamiento del concreto entre grietas (spec §4): las piezas y la
        # validación del modelo 2; sin f_c ni G_c la compresión es elástica.
        self.beta, self.f_c, self.G_c, self._kappa_c0 = parse_compression(
            owner, f_c, G_c, DEFAULT_BETA if beta is None else beta, self.f_t, self.E)
        self.elastic_base = Elastic2D(E, nu, hypothesis=hypothesis)
        self.density = density
        # IMPL-EX: tangente secante, simétrica sin aplastamiento; con él, la
        # consistente del problema con la historia fija, no simétrica (§9).
        self.IS_SYMMETRIC = self._integration_code == IMPL_EX and self.f_c is None

    @property
    def has_compression_damage(self) -> bool:
        """``True`` con ``f_c`` y ``G_c`` (aplastamiento activo, spec §4)."""
        return self.f_c is not None

    @property
    def hypothesis(self) -> str:
        return self.elastic_base.hypothesis

    # --- longitud característica (spec §2, H1) ------------------------------

    @property
    def characteristic_length_limit(self) -> float:
        """``h_max = E/((1 + ν)·|f'(0)|)``: ``ℓ_ch/(1 + ν)`` con la exponencial
        y ``0.7382·ℓ_ch/(1 + ν)`` con Hordijk (spec §2 y §8)."""
        factor = HORDIJK_H_MAX_FACTOR if self._softening_code == HORDIJK else 1.0
        return factor * self.characteristic_length / (1.0 + self.nu)

    def initial_point_state(self, characteristic_length: float) -> dict:
        h = float(characteristic_length)
        name = type(self).__name__
        if not (math.isfinite(h) and h > 0.0):
            raise ValueError(f"{name}: longitud característica h = {h} no válida (> 0).")
        h_max = self.characteristic_length_limit
        if h >= h_max:
            ley = ("0.7382·ℓ_ch/(1 + ν), ley de Hordijk" if self._softening_code == HORDIJK
                   else "ℓ_ch/(1 + ν), ley exponencial")
            raise ValueError(
                f"{name}: la longitud característica h = {h:.4g} supera el límite "
                f"h_max = {h_max:.4g} ({ley}; E = {self.E:.4g}, ν = {self.nu:.4g}, "
                f"G_f = {self.G_f:.4g}, f_t = {self.f_t:.4g}, ℓ_ch = E·G_f/f_t² = "
                f"{self.characteristic_length:.4g}): por encima, el problema local con las "
                f"dos grietas en carga puede tener varias soluciones (spec §2 y §8). Use "
                f"elementos más pequeños.")
        return self._point_state(h)

    # --- estado por punto ------------------------------------------------------

    def admissibility_scale(self, state_vars=None) -> float:
        """Esfuerzo de inicio de la grieta, ``f_t`` (ADR 0006)."""
        return self.f_t

    def _initial_row(self, h: float) -> np.ndarray:
        """Estado inicial. Con aplastamiento, ``κ⁻`` empieza en ``f_c/(3E)``
        (fin del tramo elástico), como en el modelo 2: así IMPL-EX nunca crea
        daño de compresión en un punto que la implícita deja elástico."""
        row = np.zeros(len(STATE_SCHEMA_FSC))
        row[H] = float(h)
        row[KC] = self._kappa_c0
        row[KC_PREV] = self._kappa_c0
        return row

    def _state_row(self, state_vars) -> np.ndarray:
        """Estado del punto como fila (orden de ``STATE_SCHEMA``), con los
        valores iniciales para lo que falte."""
        row = self._initial_row(0.0)
        if state_vars is None:
            return row
        for i, name in enumerate(STATE_SCHEMA_FSC):
            if name in state_vars:
                row[i] = float(state_vars[name])
        if 'w_max_n_prev' not in state_vars:
            row[WMAX_N_PREV] = row[WMAX_N]
        if 'w_max_t_prev' not in state_vars:
            row[WMAX_T_PREV] = row[WMAX_T]
        if 'kappa_c_prev' not in state_vars:
            row[KC_PREV] = row[KC]
        return row

    def _point_state(self, h: float) -> dict:
        return self._row_to_state(self._initial_row(h))

    @staticmethod
    def _row_to_state(row) -> dict:
        return {name: float(row[i]) for i, name in enumerate(STATE_SCHEMA_FSC)}

    def _missing_h(self) -> ValueError:
        return ValueError(
            f"{type(self).__name__}: la grieta empezó en un punto sin longitud "
            f"característica h. La siembran los elementos sólidos 2D (Tri3, Quad4, "
            f"Tri6, Quad8, Quad9) al construirse; con otro elemento este material no "
            f"se admite.")

    # --- respuesta -----------------------------------------------------------

    def compute_state(self, strain: np.ndarray, state_vars=None):
        S_old = self._state_row(state_vars)
        strain = np.ascontiguousarray(strain, dtype=np.float64)
        f_c, G_c, K, tol_c = compression_kernel_params(self.f_c, self.G_c, self.beta)
        sigma, C_tan, row, missing = _fsc_2d_point(
            strain, S_old, self.elastic_base.C, self.E, self.f_t, self.G_f,
            self._softening_code, self.w_c, self.admissibility_tol(state_vars),
            self.damage_max, self._integration_code, self._ratio, f_c, G_c, K, self.nu,
            self.hypothesis == 'plane_strain', tol_c)
        if missing:
            raise self._missing_h()
        return sigma, C_tan, self._row_to_state(row)

    BATCH_KERNEL = _fsc_2d_batch

    def batch_params(self) -> np.ndarray:
        f_c, G_c, K, _ = compression_kernel_params(self.f_c, self.G_c, self.beta)
        return np.array([self.E, self.f_t, self.G_f,
                         ADMISSIBILITY_TOL_ABS, ADMISSIBILITY_TOL_REL, self.damage_max,
                         float(self._integration_code), self._ratio,
                         float(self._softening_code), self.w_c, f_c, G_c, K, self.nu,
                         1.0 if self.hypothesis == 'plane_strain' else 0.0],
                        dtype=np.float64)

    def batch_matrix(self) -> np.ndarray:
        return self.elastic_base.C

    # --- componente fuera del plano (post-proceso) -------------------------

    def out_of_plane_stress(self, sigma, state_vars=None) -> float:
        """``σ_zz``: nula en esfuerzo plano; en deformación plana, la del
        concreto entre grietas (la grieta sólo está en el plano, spec §1):
        ``ν·(σ_xx + σ_yy)`` sin aplastamiento y ``σ̄_zz − d⁻·⟨σ̄_zz⟩₋`` con él.

        Con ``d⁻ > 0``, ``σ`` y ``σ̄`` tienen las mismas direcciones
        principales y ``σ_i = σ̄_i`` (``σ̄_i ≥ 0``) o ``(1 − d⁻)·σ̄_i``
        (``σ̄_i < 0``), así que ``σ̄`` se recupera de ``σ`` y de ``d⁻``, y
        ``σ̄_zz = ν·tr σ̄`` del plano (``ε^co_zz = 0``), como en el modelo 2.
        """
        if self.hypothesis != 'plane_strain':
            return 0.0
        d_c = 0.0 if state_vars is None else float(state_vars.get('damage_c', 0.0))
        if d_c == 0.0:
            return float(self.nu * (sigma[0] + sigma[1]))
        sx, sy, sxy = (float(sigma[0]), float(sigma[1]), float(sigma[2]))
        r = math.hypot(0.5 * (sx - sy), sxy)
        mean = 0.5 * (sx + sy)
        tr_eff = 0.0
        for s in (mean + r, mean - r):
            tr_eff += s if s >= 0.0 else s / (1.0 - d_c)
        szz_eff = self.nu * tr_eff
        return szz_eff if szz_eff > 0.0 else (1.0 - d_c) * szz_eff

    def batch_out_of_plane_stress(self, sigma, S):
        if self.hypothesis != 'plane_strain':
            return np.zeros(sigma.shape[0])
        plain = self.nu * (sigma[:, 0] + sigma[:, 1])
        d_c = S[:, DC]
        if not np.any(d_c > 0.0):
            return plain
        sx, sy, sxy = sigma[:, 0], sigma[:, 1], sigma[:, 2]
        r = np.hypot(0.5 * (sx - sy), sxy)
        mean = 0.5 * (sx + sy)
        tr_eff = np.zeros(sigma.shape[0])
        for s in (mean + r, mean - r):
            tr_eff += np.where(s >= 0.0, s, s / (1.0 - d_c))
        szz_eff = self.nu * tr_eff
        crushed = np.where(szz_eff > 0.0, szz_eff, (1.0 - d_c) * szz_eff)
        return np.where(d_c > 0.0, crushed, plain)
