"""``FixedSmearedCrack2D`` (modelo 3), bloque 2: aplastamiento del concreto entre grietas.

Aceptación del bloque 2 de la spec ``FixedSmearedCrack2D`` (§4, §8b, §9) en el
punto material. Cada prueba lleva en su nombre el invariante que protege.

Referencias escritas aquí, aparte del código: la parábola de Feenstra y la
medida de Drucker-Prager sobre la parte negativa del esfuerzo efectivo (spec
del modelo 2, §3), la matriz ``N`` de la spec §1 y la elasticidad plana. El
modelo 2 (``TensionCompressionDamage2D``) sirve de referencia independiente
sin grieta y con la grieta cerrada; el bloque 1 (el mismo material sin
``f_c``), mientras el concreto no se aplasta.

Concreto usual en SI: E = 30 GPa, ν = 0.2, f_t = 3 MPa, G_f = 100 N/m,
f_c = 30 MPa, G_c = 15 kN/m, β = 1.16.
"""
import math

import numpy as np
import pytest
from scipy import integrate, optimize

import solidum  # noqa: F401
from solidum.constants import CRACK_BAND_DAMAGE_MAX
from solidum.materials.damage_tension_compression_2d import TensionCompressionDamage2D
from solidum.materials.fixed_smeared_crack_2d import FixedSmearedCrack2D

E, NU, FT, GF = 30.0e9, 0.2, 3.0e6, 100.0
FC, GC, BETA = 30.0e6, 15.0e3, 1.16
H = 0.025
CAP = CRACK_BAND_DAMAGE_MAX
PS, PE = "plane_stress", "plane_strain"
KC3 = FC / (3.0 * E)
KCP = 5.0 * FC / (3.0 * E)
K_DP = math.sqrt(2.0) * (BETA - 1.0) / (2.0 * BETA - 1.0)

# Ley de Hordijk (Cornelissen, Hordijk y Reinhardt 1986), escrita aquí aparte.
C1, C2 = 3.0, 6.93


def g_hordijk(x):
    return (1.0 + C1**3 * x**3) * math.exp(-C2 * x) - x * (1.0 + C1**3) * math.exp(-C2)


KW = 1.0 / integrate.quad(g_hordijk, 0.0, 1.0, epsabs=0.0, epsrel=1e-13)[0]
WC = KW * GF / FT


def f_law(w, law="exponential"):
    if law == "hordijk":
        return 0.0 if w >= WC else FT * g_hordijk(w / WC)
    return FT * math.exp(-FT * w / GF)


def c_max(h):
    return h * CAP / ((1.0 - CAP) * E)


def compliance(w, h, law="exponential"):
    if w <= 0.0:
        return 0.0
    f = f_law(w, law)
    cm = c_max(h)
    return cm if w >= cm * f else w / f


def omega(c, h):
    return E * c / (h + E * c)


def sigma_c(kappa, h):
    """Parábola de Feenstra con G_c/h (spec del modelo 2, §3)."""
    ku = KCP + 1.5 * GC / (h * FC)
    if kappa <= KC3:
        return E * kappa
    if kappa <= KCP:
        xi = (kappa - KC3) / (KCP - KC3)
        return FC / 3.0 * (1.0 + 4.0 * xi - 2.0 * xi * xi)
    if kappa < ku:
        eta = (kappa - KCP) / (ku - KCP)
        return FC * (1.0 - eta * eta)
    return 0.0


def eps_eq(s1, s2, s3):
    """Medida de Drucker-Prager de la parte negativa de un estado principal
    3D, normalizada a |σ|/E en compresión uniaxial (modelo 2, §3)."""
    p = np.minimum(np.array([s1, s2, s3], float), 0.0)
    mean = p.sum() / 3.0
    dev = p - mean
    tau = math.sqrt(2.0 / 3.0 * 0.5 * (dev @ dev))
    return 3.0 * (K_DP * mean + tau) / ((math.sqrt(2.0) - K_DP) * E)


def N_spec(theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c * c, s * s, -c * s], [s * s, c * c, c * s],
                     [2 * c * s, -2 * c * s, c * c - s * s]])


def Ce(hyp=PS):
    if hyp == PS:
        return E / (1 - NU**2) * np.array([[1, NU, 0], [NU, 1, 0], [0, 0, (1 - NU) / 2]])
    k = E / ((1 + NU) * (1 - 2 * NU))
    return k * np.array([[1 - NU, NU, 0], [NU, 1 - NU, 0], [0, 0, (1 - 2 * NU) / 2]])


def mat(hyp=PS, crush=True, **kw):
    if crush:
        kw.setdefault("f_c", FC)
        kw.setdefault("G_c", GC)
    return FixedSmearedCrack2D(E, NU, FT, GF, hypothesis=hyp, **kw)


def model2(hyp=PS, **kw):
    return TensionCompressionDamage2D(E, NU, FT, GF, hypothesis=hyp, f_c=FC, G_c=GC, **kw)


def cracked_state(m, theta, wn, wt=0.0, h=H, law="exponential"):
    st = m.initial_point_state(h)
    st.update(cracked=1.0, theta=float(theta), w_max_n=float(wn), w_max_t=float(wt),
              w_max_n_prev=float(wn), w_max_t_prev=float(wt),
              damage_n=omega(compliance(wn, h, law), h),
              damage_t=omega(compliance(wt, h, law), h))
    return st


def principal(sig):
    return 0.5 * (sig[0] + sig[1]) + math.hypot(0.5 * (sig[0] - sig[1]), sig[2])


def voigt(t):
    """Deformación en Voigt, con γ de ingeniería."""
    return np.array([t[0, 0], t[1, 1], 2.0 * t[0, 1]])


def voigt_s(t):
    """Esfuerzo en Voigt (cortante tensorial, sin factor)."""
    return np.array([t[0, 0], t[1, 1], t[0, 1]])


def crushed_stress(kappa, h):
    """−σ en compresión uniaxial efectiva con el tope de d⁻: (1 − d⁻)·E·κ."""
    d = min(1.0 - sigma_c(kappa, h) / (E * kappa), CAP) if kappa > KC3 else 0.0
    return (1.0 - d) * E * kappa


# ---------------------------------------------------------------------------
# Contrato con aplastamiento
# ---------------------------------------------------------------------------

def test_contrato_con_aplastamiento_y_estado_inicial():
    """κ⁻ empieza en f_c/(3E), como en el modelo 2; con aplastamiento la
    tangente de IMPL-EX deja de ser simétrica (§9)."""
    m = mat()
    st = m.initial_point_state(H)
    assert st["kappa_c"] == st["kappa_c_prev"] == pytest.approx(KC3, rel=1e-15)
    assert st["damage_c"] == 0.0
    assert mat(integration="impl-ex").IS_SYMMETRIC is False
    assert mat(crush=False, integration="impl-ex").IS_SYMMETRIC is True
    assert mat(crush=False).initial_point_state(H)["kappa_c"] == 0.0


# ---------------------------------------------------------------------------
# Sin grieta: el modelo 2 con d⁺ = 0 (compression_like_model2)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
@pytest.mark.parametrize("h", [0.010, 0.050, 0.100])
@pytest.mark.parametrize("camino", ["uniaxial", "biaxial", "hidrostatica"])
def test_compresion_sin_grieta_igual_que_el_modelo_2(hyp, integration, h, camino):
    """compression_like_model2: σ, tangente y d⁻ iguales a los del modelo 2 con
    la misma integración, pasando el pico hasta la rama descendente."""
    m3 = mat(hyp, integration=integration)
    m2 = model2(hyp, integration=integration)
    if integration == "impl-ex":
        m3.set_pseudo_time_ratio(1.0)
        m2.set_pseudo_time_ratio(1.0)
    s3, s2 = m3.initial_point_state(h), m2.initial_point_state(h)
    d = {"uniaxial": [-1.0, NU, 0.0], "biaxial": [-1.0, -1.0, 0.0],
         "hidrostatica": [-1.0, -1.0, 0.0]}[camino]
    if camino == "uniaxial" and hyp == PE:
        d = [-1.0, NU / (1.0 - NU), 0.0]
    if camino == "hidrostatica" and hyp == PS:
        d = [-1.0, -0.5, 0.0]
    d_max = 0.0
    for k in range(1, 300):
        eps = np.array(d) * k * 1.0e-5
        sig3, C3, s3 = m3.compute_state(eps, s3)
        sig2, C2, s2 = m2.compute_state(eps, s2)
        np.testing.assert_allclose(sig3, sig2, rtol=1e-10, atol=1e-10 * FC)
        np.testing.assert_allclose(C3, C2, rtol=1e-10, atol=1e-10 * E)
        assert s3["damage_c"] == pytest.approx(s2["damage_c"], rel=1e-10, abs=1e-14)
        assert s3["kappa_c"] == pytest.approx(s2["kappa_c"], rel=1e-12)
        assert s3["cracked"] == 0.0
        d_max = max(d_max, s3["damage_c"])
    assert d_max > 0.3                                  # el camino aplasta


# ---------------------------------------------------------------------------
# Sin aplastamiento, el bloque 1 tal cual (compression_unchanged_below_crushing)
# ---------------------------------------------------------------------------

def _camino_aleatorio(rng, n_tramos=5, n_pasos=40, amp=1.2e-4, comp=0.0):
    """Incrementos de deformación por tramos con direcciones al azar: tracción,
    compresión (escalada por ``comp``) y cortante."""
    out = []
    eps = np.zeros(3)
    for _ in range(n_tramos):
        a = rng.uniform(0.0, math.pi)
        n = np.array([math.cos(a), math.sin(a)])
        t = np.array([-n[1], n[0]])
        tipo = rng.integers(0, 4)
        if tipo == 0:
            d = np.outer(n, n) - 0.2 * np.outer(t, t)
        elif tipo == 1:
            d = -comp * np.outer(t, t) + 0.3 * np.outer(n, n)
        elif tipo == 2:
            d = 0.5 * (np.outer(n, t) + np.outer(t, n)) + 0.3 * np.outer(n, n)
        else:
            d = np.outer(n, n) + np.outer(t, t)
        d = voigt(d) * amp * rng.uniform(0.3, 1.0) / n_pasos
        for _ in range(n_pasos):
            eps = eps + d
            out.append(eps.copy())
    return out


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_sin_aplastamiento_el_bloque_2_reproduce_el_bloque_1_bit_a_bit(hyp, integration):
    """compression_unchanged_below_crushing: mientras κ⁻ ≤ f_c/(3E) el predictor
    de §8b da el bloque 1 bit a bit (σ, tangente y estado, salvo kappa_c y
    kappa_c_prev, que con f_c se actualizan aunque no haya daño)."""
    rng = np.random.default_rng(11)
    a, b = mat(hyp, integration=integration), mat(hyp, crush=False, integration=integration)
    a.set_pseudo_time_ratio(1.0)
    b.set_pseudo_time_ratio(1.0)
    otros = [k for k in a.STATE_SCHEMA if k not in ("kappa_c", "kappa_c_prev")]
    n_agrietados = 0
    for _ in range(12):
        sa, sb = a.initial_point_state(H), b.initial_point_state(H)
        for eps in _camino_aleatorio(rng, comp=1.0):
            ra, rb = a.compute_state(eps, sa), b.compute_state(eps, sb)
            assert np.array_equal(ra[0], rb[0]) and np.array_equal(ra[1], rb[1])
            sa, sb = ra[2], rb[2]
            assert all(sa[k] == sb[k] for k in otros)
            assert sa["damage_c"] == 0.0 and sa["kappa_c"] <= KC3
        n_agrietados += sa["cracked"] > 0.0
    assert n_agrietados >= 6                            # los caminos agrietan


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("law", ["exponential", "hordijk"])
def test_problema_acoplado_reproduce_el_bloque_1_con_el_concreto_casi_elastico(hyp, law):
    """predictor_independence: con κ⁻ apenas sobre f_c/(3E) (d⁻ ≈ 10⁻¹³, el
    concreto casi elástico) el problema acoplado de §8b, otro algoritmo que el
    anidado del bloque 1, da la misma solución en estados con una y dos
    grietas, cortante, descarga, cierre y reapertura."""
    rng = np.random.default_rng(7)
    a = mat(hyp, softening=law)
    b = mat(hyp, crush=False, softening=law)
    k_eps = KC3 * (1.0 + 1e-6)
    n_casos = 0
    for _ in range(10):
        sb = b.initial_point_state(H)
        for eps in _camino_aleatorio(rng, n_pasos=30):
            sig_b, C_b, sb_new = b.compute_state(eps, sb)
            sa = dict(sb, kappa_c=k_eps, kappa_c_prev=k_eps)
            sig_a, C_a, sa_new = a.compute_state(eps, sa)
            if sb["cracked"] > 0.0:
                assert 0.0 < sa_new["damage_c"] < 1e-11       # tomó el camino acoplado
                n_casos += 1
                np.testing.assert_allclose(sig_a, sig_b, rtol=1e-10, atol=1e-10 * FT)
                for k in ("w_max_n", "w_max_t", "w_n", "w_t", "slip"):
                    assert sa_new[k] == pytest.approx(sb_new[k], rel=1e-10, abs=1e-15), k
            sb = sb_new
    assert n_casos > 500


# ---------------------------------------------------------------------------
# Aplastamiento paralelo a una grieta abierta (§4)
# ---------------------------------------------------------------------------

def _paralelo(hyp, w, h, n=3000, eps_t_max=None):
    """Grieta en n (θ = 0) con apertura máxima ``w`` y la apertura fija en
    ``w`` mientras se comprime en t: ε_n(ε_t) = w/h + (σ_nn − C12·ε_t)/C11
    (fijar ε_n cierra la grieta). Devuelve (ε_t, σ_nn, σ_tt, w_n, κ⁻)."""
    m = mat(hyp, softening="hordijk")
    st = cracked_state(m, 0.0, w, h=h, law="hordijk")
    sn = w / compliance(w, h, "hordijk") if w > 0 else 0.0
    C = Ce(hyp)
    if eps_t_max is None:
        eps_t_max = 1.5 * (KCP + 1.5 * GC / (h * FC))
    out = []
    for et in np.linspace(0.0, -eps_t_max, n + 1)[1:]:
        en = w / h + (sn - C[0, 1] * et) / C[0, 0]
        sig, _, st = m.compute_state(np.array([en, et, 0.0]), st)
        out.append((et, sig[0], sig[1], st["w_n"], st["kappa_c"]))
    return np.array(out), sn


@pytest.mark.parametrize("h", [0.010, 0.050])
@pytest.mark.parametrize("fw", [0.1, 0.5, 1.0])
def test_aplastamiento_paralelo_a_grieta_abierta_pico_f_c_con_cualquier_apertura(h, fw):
    """crushing_parallel_to_open_crack (esfuerzo plano): −σ_tt = σ_c(κ⁻) con
    κ⁻ = |E·ε_t + ν·σ_nn|/E, pico f_c, σ_nn = f(w) la del bloque 1, la
    apertura fija, y el área de la rama descendente G_c/h (cuadratura)."""
    w = fw * WC
    r, sn = _paralelo(PS, w, h)
    et, s_nn, s_tt, wn, kc = r.T
    s_bar = E * et + NU * sn                             # σ̄_tt
    kap = np.maximum.accumulate(np.maximum(-s_bar, 0.0) / E)
    d = np.array([0.0 if k <= KC3 else 1.0 - crushed_stress(k, h) / (E * k) for k in kap])
    ref = np.where(s_bar >= 0.0, -s_bar, -(1.0 - d) * s_bar)
    np.testing.assert_allclose(-s_tt, ref, rtol=1e-10, atol=1e-10 * FC)
    np.testing.assert_allclose(kc, np.maximum(kap, KC3), rtol=1e-12)
    np.testing.assert_allclose(s_nn, sn, rtol=1e-9, atol=1e-9 * FT)
    np.testing.assert_allclose(wn, w, rtol=1e-9)
    assert (-s_tt).max() == pytest.approx(FC, rel=1e-6)
    k = int(np.argmax(-s_tt))
    area = np.trapezoid(-s_tt[k:], kap[k:])
    assert area == pytest.approx(GC / h, rel=1e-3)


def _pico_deformacion_plana(sn, h):
    """Pico de −σ_tt en deformación plana con σ_nn fija (§4, H7 de la revisión
    del bloque 2): σ̄⁻ = (0, −S, ν·(σ_nn − S)), κ⁻ = ε⁻_eq y −σ_tt = σ_c(κ⁻)·S/(E·κ⁻),
    maximizado en S."""
    def menos_s(S):
        k = eps_eq(0.0, -S, NU * (sn - S))
        return -sigma_c(k, h) * S / (E * k)
    res = optimize.minimize_scalar(menos_s, bounds=(0.5 * FC, 3.0 * FC), method="bounded",
                                   options={"xatol": 1e-6 * FC})
    return -res.fun


@pytest.mark.parametrize("fw", [0.1, 1.0])
def test_aplastamiento_paralelo_en_deformacion_plana_pico_f_c_sobre_phi(fw):
    """crushing_parallel_plane_strain: el pico es f_c/φ(ρ), ρ = ν·(1 − σ_nn/|σ̄_tt|);
    1.1397·f_c con σ_nn ≈ 0 (w = w_c), algo menos con la grieta transmitiendo."""
    h = 0.025
    r, sn = _paralelo(PE, fw * WC, h, n=6000)
    pico = (-r[:, 2]).max()
    assert pico == pytest.approx(_pico_deformacion_plana(sn, h), rel=1e-6)
    if fw == 1.0:
        assert pico / FC == pytest.approx(1.1397, abs=5e-5)
    np.testing.assert_allclose(r[:, 3], fw * WC, rtol=1e-9)


# ---------------------------------------------------------------------------
# Grieta cerrada, aplastar y después agrietar, inicio tras aplastar
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("hyp", [PS, PE])
def test_grieta_cerrada_se_aplasta_como_el_modelo_2(hyp):
    """closed_crack_crushing: con la grieta cerrada, la compresión normal a ella
    sigue la curva del modelo 2 (d⁺ = 0); w_n = 0 y w_max sin cambio."""
    wmax = 0.2 * GF / FT
    m3 = mat(hyp)
    m2 = model2(hyp)
    s3 = cracked_state(m3, 0.0, wmax)
    s2 = m2.initial_point_state(H)
    lat = NU if hyp == PS else NU / (1.0 - NU)
    for k in range(1, 300):
        eps = np.array([-1.0, lat, 0.0]) * k * 1.0e-5
        sig3, _, s3 = m3.compute_state(eps, s3)
        sig2, _, s2 = m2.compute_state(eps, s2)
        np.testing.assert_allclose(sig3, sig2, rtol=1e-10, atol=1e-10 * FC)
        assert s3["w_n"] == 0.0 and s3["w_max_n"] == wmax
    assert s3["damage_c"] > 0.5


def _aplastar_y_descargar(m, ang, h=H, d_obj=0.5):
    """Compresión uniaxial de esfuerzo efectivo en la dirección ``ang`` hasta
    d⁻ ≥ d_obj en la rama descendente y descarga a deformación nula."""
    st = m.initial_point_state(h)
    n = np.array([math.cos(ang), math.sin(ang)])
    S_inv = np.linalg.inv(Ce(PS))
    s = 0.0
    while st["damage_c"] < d_obj or st["kappa_c"] < KCP:
        s += 0.002 * FC
        eps = S_inv @ voigt_s(-s * np.outer(n, n))
        _, _, st = m.compute_state(eps, st)
    sig, _, st = m.compute_state(np.zeros(3), st)
    np.testing.assert_allclose(sig, 0.0, atol=1e-9 * FC)
    return st


def test_aplastado_y_despues_traccion_forma_la_grieta_en_f_t_y_disipa_G_f():
    """crushing_then_cracking: sin acoplamiento entre los daños (§4): tras
    aplastar a d⁻ ≈ 0.5, la grieta se forma con σ_xx = f_t y θ = 0, d⁻ no
    cambia y la energía hasta el tope es D_tope (la de tension_energy)."""
    m = mat()
    st = _aplastar_y_descargar(m, 0.0)
    d_c = st["damage_c"]
    S_inv = np.linalg.inv(Ce(PS))
    for s in np.linspace(0.0, FT, 51)[1:-1]:
        _, _, st1 = m.compute_state(S_inv @ np.array([s, 0, 0]), st)
        assert st1["cracked"] == 0.0
    _, _, st1 = m.compute_state(S_inv @ np.array([FT * (1 + 1e-6), 0, 0]), st)
    assert st1["cracked"] == 1.0 and st1["theta"] == 0.0
    # energía: tracción uniaxial monótona hasta pasar el tope (como el bloque 1)
    wt = optimize.brentq(lambda w: w - c_max(H) * f_law(w), 0.0, 60 * GF / FT,
                         xtol=1e-300, rtol=4 * np.finfo(float).eps)
    area = integrate.quad(f_law, 0.0, wt, epsabs=0.0, epsrel=1e-13, limit=200)[0]
    d_tope = (area - 0.5 * f_law(wt) * wt) / H
    path = [S_inv @ np.array([s, 0, 0]) for s in np.linspace(0, FT, 51)]
    path += [S_inv @ np.array([max(f_law(w), w / c_max(H)), 0, 0]) + np.array([w / H, 0, 0])
             for w in np.linspace(0.0, 1.2 * wt, 20001)[1:]]
    W, s_prev, e_prev = 0.0, np.zeros(3), np.zeros(3)
    for eps in path:
        sig, _, st = m.compute_state(eps, st)
        W += 0.5 * (sig + s_prev) @ (eps - e_prev)
        s_prev, e_prev = sig, eps
    assert W - 0.5 * s_prev @ e_prev == pytest.approx(d_tope, rel=1e-6)
    assert st["damage_c"] == d_c


@pytest.mark.parametrize("psi", [math.pi / 2, math.pi / 2 + 0.35, math.pi / 2 - 0.35])
def test_inicio_de_grieta_tras_aplastar_con_sigma_bar_I_y_su_direccion(psi):
    """crack_initiation_after_crushing (§2): tras aplastar en x, una tracción
    uniaxial efectiva en la dirección ψ forma la grieta cuando σ̄_I = f_t, con
    θ = ψ (la dirección de σ̄)."""
    m = mat()
    st = _aplastar_y_descargar(m, 0.0)
    n = np.array([math.cos(psi), math.sin(psi)])
    S_inv = np.linalg.inv(Ce(PS))
    eps_lo = S_inv @ voigt_s(FT * (1 - 1e-7) * np.outer(n, n))
    eps_hi = S_inv @ voigt_s(FT * (1 + 1e-7) * np.outer(n, n))
    _, _, st_lo = m.compute_state(eps_lo, st)
    _, _, st_hi = m.compute_state(eps_hi, st)
    assert st_lo["cracked"] == 0.0 and st_hi["cracked"] == 1.0
    dth = (st_hi["theta"] - psi + math.pi / 2) % math.pi - math.pi / 2
    assert abs(dth) < 1e-9


# ---------------------------------------------------------------------------
# Tangente con aplastamiento (tangent_fd, secant_tangent)
# ---------------------------------------------------------------------------

def _fd(m, eps, st, rel=1e-7):
    C = np.zeros((3, 3))
    d = rel * max(np.abs(eps).max(), 1e-5)
    for j in range(3):
        ep, em = eps.copy(), eps.copy()
        ep[j] += d
        em[j] -= d
        C[:, j] = (m.compute_state(ep, st)[0] - m.compute_state(em, st)[0]) / (2.0 * d)
    return C


def _estados_con_aplastamiento(m, caso, rng, n_caminos=8):
    """Estados (ε, estado confirmado) con aplastamiento en carga en el paso
    siguiente, por caso: sin grieta, grieta cerrada, grieta abierta paralela a
    la compresión, grieta inclinada con cortante."""
    out = []
    for _ in range(n_caminos):
        st = m.initial_point_state(H)
        if caso == "sin_grieta":
            a = rng.uniform(0, math.pi)
            n = np.array([math.cos(a), math.sin(a)])
            d = voigt(-np.outer(n, n) + 0.1 * rng.uniform(-1, 1) * np.eye(2)) * 1e-5
            eps = np.zeros(3)
        else:
            theta = rng.uniform(-1.2, 1.2)
            st = cracked_state(m, theta, rng.uniform(0.05, 0.5) * GF / FT)
            n = np.array([math.cos(theta), math.sin(theta)])
            t = np.array([-n[1], n[0]])
            if caso == "grieta_cerrada":
                d = voigt(-np.outer(n, n) + 0.05 * np.outer(t, t)) * 1e-5
            elif caso == "grieta_abierta":
                d = voigt(-np.outer(t, t) + 0.25 * np.outer(n, n)) * 1e-5
            else:
                b = theta + rng.uniform(0.2, 0.6)
                n2 = np.array([math.cos(b), math.sin(b)])
                t2 = np.array([-n2[1], n2[0]])
                d = voigt(-np.outer(t2, t2) + 0.2 * np.outer(n2, n2)
                          + 0.15 * (np.outer(n2, t2) + np.outer(t2, n2))) * 1e-5
            sig_loc = N_spec(theta).T @ np.zeros(3)
            eps = N_spec(theta) @ (np.linalg.solve(Ce(m.hypothesis), sig_loc)
                                   + np.array([st["w_max_n"], 0, 0]) / H)
        for k in range(400):
            eps = eps + d
            _, _, st_new = m.compute_state(eps, st)
            if st_new["kappa_c"] > st["kappa_c"] and 0.05 < st_new["damage_c"] < 0.95:
                if k % 7 == 0:
                    out.append((eps + d, st_new))
            st = st_new
    return out


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("caso", ["sin_grieta", "grieta_cerrada", "grieta_abierta",
                                  "grieta_inclinada"])
def test_tangente_con_aplastamiento_igual_a_diferencias_finitas(hyp, caso):
    """tangent_fd del bloque 2 (§8b): |C_FD − C_tan| ≤ 10⁻⁶·max|C_FD| (10⁻⁵
    con la grieta abierta y aplastamiento en carga), lejos de los quiebres
    (dos pasos de diferencias que coinciden); no simétrica con aplastamiento."""
    rng = np.random.default_rng(3)
    m = mat(hyp)
    tol = 1e-5 if caso in ("grieta_abierta", "grieta_inclinada") else 1e-6
    n_ok, asim = 0, 0.0
    for eps, st in _estados_con_aplastamiento(m, caso, rng):
        sig, C, st_new = m.compute_state(eps, st)
        if not (st_new["kappa_c"] > st["kappa_c"]):
            continue
        C1, C2 = _fd(m, eps, st, 1e-7), _fd(m, eps, st, 2e-7)
        if np.abs(C1 - C2).max() > 1e-3 * np.abs(C1).max():
            continue                                   # junto a un quiebre
        assert np.abs(C1 - C).max() <= tol * np.abs(C1).max()
        asim = max(asim, np.abs(C - C.T).max() / np.abs(C).max())
        n_ok += 1
    assert n_ok >= 10
    assert asim > 1e-3


@pytest.mark.parametrize("hyp", [PS, PE])
def test_implex_con_daño_fijo_tangente_igual_a_diferencias_finitas_y_no_simetrica(hyp):
    """secant_tangent con d̃⁻ > 0 (§9): fuera de un paso (historia congelada) la
    tangente es la consistente del problema con la historia fija: igual a
    diferencias finitas y no simétrica."""
    rng = np.random.default_rng(5)
    m = mat(hyp, integration="impl-ex")
    n_ok, asim = 0, 0.0
    for _ in range(40):
        theta = rng.uniform(-1.2, 1.2)
        st = cracked_state(m, theta, rng.uniform(0.05, 0.5) * GF / FT)
        st.update(damage_c=rng.uniform(0.2, 0.9), kappa_c=2 * KCP, kappa_c_prev=2 * KCP)
        a = theta + rng.uniform(-0.8, 0.8)
        n = np.array([math.cos(a), math.sin(a)])
        t = np.array([-n[1], n[0]])
        # tracción a través de la grieta y compresión paralela: esfuerzos
        # principales efectivos de signos distintos (la parte no lineal, §9)
        eps = voigt(-1.5e-3 * np.outer(t, t) + rng.uniform(0.6, 1.5) * 1e-3 * np.outer(n, n)
                    + rng.uniform(-3, 3) * 1e-4 * (np.outer(n, t) + np.outer(t, n)))
        m.set_pseudo_time_ratio(None)
        sig, C, st2 = m.compute_state(eps, st)
        assert st2 == st                                # fuera de un paso, congelado
        C1, C2 = _fd(m, eps, st, 1e-7), _fd(m, eps, st, 2e-7)
        if np.abs(C1 - C2).max() > 1e-3 * np.abs(C1).max():
            continue
        assert np.abs(C1 - C).max() <= 1e-6 * np.abs(C1).max()
        asim = max(asim, np.abs(C - C.T).max() / np.abs(C).max())
        n_ok += 1
    assert n_ok >= 20 and asim > 1e-3


# ---------------------------------------------------------------------------
# IMPL-EX con aplastamiento (§9)
# ---------------------------------------------------------------------------

def _camino_paralelo_implex(n, integration, kappa_de_t, sn_w=0.3):
    """Compresión paralela a una grieta abierta (σ_nn fija), κ⁻(t) dado."""
    m = mat(PS, integration=integration)
    w = sn_w * GF / FT
    st = cracked_state(m, 0.0, w)
    sn = w / compliance(w, H)
    C = Ce(PS)
    m.set_pseudo_time_ratio(1.0)
    out = []
    for t in np.linspace(0.0, 1.0, n + 1)[1:]:
        et = -(E * kappa_de_t(t) + NU * sn) / E
        en = w / H + (sn - C[0, 1] * et) / C[0, 0]
        sig, _, st = m.compute_state(np.array([en, et, 0.0]), st)
        out.append(sig[1])
    return np.array(out)


@pytest.mark.parametrize("forma, orden", [("lineal_desde_el_umbral", 2.0),
                                           ("descarga_y_recarga", 1.0)])
def test_implex_con_grieta_abierta_paralela_convergencia_a_la_implicita(forma, orden):
    """first_order_convergence con aplastamiento: por forma del camino, como el
    modelo 2. Con la grieta abierta paralela y σ_nn fija κ⁻ es lineal en el
    paso: segundo orden (razón 4); con descarga y recarga, primer orden."""
    if forma == "lineal_desde_el_umbral":
        def kappa(t):
            return KC3 * (1.0 + 39.0 * t)
    else:
        def kappa(t):
            return KC3 * np.interp(t, [0.0, 0.4, 0.55, 1.0], [1.0, 20.0, 10.0, 40.0])
    errores = [np.abs(_camino_paralelo_implex(n, "impl-ex", kappa)
                      - _camino_paralelo_implex(n, "implicit", kappa)).max()
               for n in (100, 200, 400, 800)]
    razones = np.array(errores[:-1]) / np.array(errores[1:])
    np.testing.assert_allclose(razones, 2.0 ** orden, rtol=0.15)


def test_implex_d_tilde_no_decrece_y_el_indicador_incluye_el_aplastamiento():
    """history_monotone con aplastamiento: d̃⁻ nunca decrece con carga, descarga
    y recarga, y damage_error = |d⁻(κ⁻) − d̃⁻| en el primer paso con daño."""
    m = mat(integration="impl-ex")
    m.set_pseudo_time_ratio(1.0)
    st = m.initial_point_state(H)
    S_inv = np.linalg.inv(Ce(PS))
    d_prev = 0.0
    primero = True
    for s in list(np.linspace(0.1, 1.0, 40)) + list(np.linspace(1.0, 0.3, 15)) \
            + list(np.linspace(0.3, 0.9, 20)):
        eps = S_inv @ np.array([-s * 2.0 * FC, 0.0, 0.0])
        _, _, st_new = m.compute_state(eps, st)
        assert st_new["damage_c"] >= d_prev
        if primero and st_new["kappa_c"] > KC3 * (1 + 1e-9):
            # primer paso con carga sobre el umbral: d̃⁻ = 0 (extrapola desde el umbral)
            d_impl = 1.0 - sigma_c(st_new["kappa_c"], H) / (E * st_new["kappa_c"])
            assert st_new["damage_error"] == pytest.approx(abs(d_impl - st_new["damage_c"]),
                                                           rel=1e-9, abs=1e-15)
            primero = False
        d_prev = st_new["damage_c"]
        st = st_new
    assert d_prev > 0.3


# ---------------------------------------------------------------------------
# Paso idempotente, robustez del problema local, falla local
# ---------------------------------------------------------------------------

def test_paso_implicito_idempotente_con_aplastamiento_y_el_intento_rechazado_no_fija_nada():
    """implicit_step_idempotent con aplastamiento: el mismo intento dos veces da
    el mismo estado; un intento rechazado no fija κ⁻ ni θ."""
    m = mat()
    st = _aplastar_y_descargar(m, 0.0, d_obj=0.3)
    S_inv = np.linalg.inv(Ce(PS))
    eps = S_inv @ voigt_s(2.0 * FT * np.diag([0.2, 1.0]) - 0.5 * FC * np.diag([1.0, 0.0]))
    r1 = m.compute_state(eps, st)
    r2 = m.compute_state(eps, st)
    assert np.array_equal(r1[0], r2[0]) and r1[2] == r2[2]
    # un intento con otra deformación (rechazado) no cambia el resultado
    m.compute_state(1.7 * eps, st)
    r3 = m.compute_state(eps, st)
    assert np.array_equal(r1[0], r3[0]) and r1[2] == r3[2]
    assert st["cracked"] == 0.0


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("pasos", [60, 8])
def test_problema_local_con_aplastamiento_converge_en_caminos_al_azar(hyp, pasos):
    """local_problem_robustness_crushing: caminos al azar con aplastamiento,
    una y dos grietas, cortante, giro de direcciones e inversión de carga; h de
    10 a 50 mm. Con 60 pasos por tramo y con pasos grandes (8 por tramo, que
    ejercitan la cascada de §8b), el esfuerzo es finito en todos los pasos."""
    rng = np.random.default_rng(17)
    n_aplastados = 0
    for h in (0.010, 0.025, 0.050):
        m = mat(hyp)
        for _ in range(8):
            st = m.initial_point_state(h)
            for eps in _camino_aleatorio(rng, n_tramos=6, n_pasos=pasos, amp=3e-3, comp=1.0):
                sig, C, st = m.compute_state(eps, st)
                assert np.all(np.isfinite(sig)) and np.all(np.isfinite(C))
            n_aplastados += st["damage_c"] > 0.0
    assert n_aplastados >= 12


def _kappa_para_d(d_obj, h):
    """κ⁻ de la rama descendente con d⁻ = d_obj (bisección sobre la parábola)."""
    lo, hi = KCP, KCP + 1.5 * GC / (h * FC)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if 1.0 - sigma_c(mid, h) / (E * mid) < d_obj:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _deformacion_principal(a, b, phi):
    c, s = math.cos(phi), math.sin(phi)
    return np.array([a * c * c + b * s * s, a * s * s + b * c * c, 2.0 * (a - b) * c * s])


def _caminos_por_vertices(rng, h, n_caminos):
    """Caminos de cinco tramos entre vértices al azar (los de la verificación
    del bloque 2): tracción hasta una grieta ancha (ε hasta f_t/E + w_c/h),
    compresión hasta d⁻ = 0.9, tracción con compresión, cortante y
    bitracción, con las direcciones principales giradas entre tramos."""
    k90 = _kappa_para_d(0.9, h)
    et = FT / E + WC / h
    for _ in range(n_caminos):
        pts = [np.zeros(3)]
        phi0 = rng.uniform(0, math.pi)
        for _ in range(5):
            tipo = rng.choice(["T", "C", "TC", "S", "TT", "C"])
            u = rng.random()
            phi = phi0 if u < 0.35 else (phi0 + rng.normal(0, 0.2) if u < 0.6
                                         else rng.uniform(0, math.pi))
            if tipo == "T":
                a = rng.uniform(0.05, 1) * et
                b = rng.uniform(-0.3, 0.3) * a
            elif tipo == "C":
                a = -rng.uniform(0.2, 1) * k90
                b = rng.uniform(-0.4, 0.3) * a
            elif tipo == "TC":
                a = rng.uniform(0.05, 1) * et
                b = -rng.uniform(0.2, 1) * k90
            elif tipo == "S":
                a = rng.uniform(0.1, 1) * min(et, k90)
                b = -a
            else:
                a = rng.uniform(0.05, 1) * et
                b = rng.uniform(0.05, 1) * et
            pts.append(_deformacion_principal(a, b, phi))
        yield pts


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("law", ["exponential", "hordijk"])
def test_problema_local_con_pasos_grandes_que_necesitan_la_cascada(hyp, law):
    """local_problem_robustness_crushing, caso de pasos grandes (§8b, H1 de la
    revisión del bloque 2): tres pasos por tramo entre vértices con grietas
    anchas que descargan y aplastamiento hasta d⁻ = 0.9. Sin la cascada (o
    arrancando desde x = 0) hay pasos sin solución (5–44 por cada 13 500,
    medido en la verificación); con ella, ninguno."""
    rng = np.random.default_rng(6)
    n, n_acoplados = 0, 0
    for h in (0.010, 0.050, 0.100):
        m = mat(hyp, softening=law)
        for pts in _caminos_por_vertices(rng, h, 300):
            st = m.initial_point_state(h)
            for k in range(5):
                for j in (1, 2, 3):
                    eps = pts[k] + (pts[k + 1] - pts[k]) * j / 3.0
                    sig, C, st = m.compute_state(eps, st)
                    assert np.all(np.isfinite(sig)) and np.all(np.isfinite(C))
                    n += 1
                    n_acoplados += st["cracked"] > 0.0 and st["damage_c"] > 0.0
    assert n == 3 * 300 * 15 and n_acoplados > 3000


@pytest.mark.parametrize("hyp", [PS, PE])
def test_estado_confirmado_vuelto_a_evaluar_en_el_camino_acoplado_no_carga(hyp):
    """Al volver a evaluar el estado confirmado con la misma deformación (lo que
    hace el ensamblaje de la tangente al empezar cada paso), el camino
    acoplado no hace crecer w_max ni κ⁻ y devuelve la tangente de descarga,
    la misma que con una deformación apenas menor. Lo protegen la banda δ de
    las aperturas y la tolerancia de κ⁻ (H1 de la verificación del bloque 2:
    sin ella κ⁻ crecía por redondeo en el 5–7 % de los pasos)."""
    rng = np.random.default_rng(31)
    m = mat(hyp)
    n = 0
    for _ in range(10):
        st = m.initial_point_state(H)
        eps_prev = np.zeros(3)
        for eps in _camino_aleatorio(rng, n_tramos=6, n_pasos=60, amp=3e-3, comp=1.0):
            _, _, st_new = m.compute_state(eps, st)
            acoplado = st_new["cracked"] > 0.0 and st_new["damage_c"] > 0.0
            if acoplado and (st_new["kappa_c"] > st["kappa_c"]
                             or st_new["w_max_n"] > st["w_max_n"]):
                _, C_re, st_re = m.compute_state(eps, st_new)
                for k in ("kappa_c", "w_max_n", "w_max_t"):
                    assert st_re[k] == st_new[k], k
                _, C_desc, st_d = m.compute_state(eps - 1e-6 * (eps - eps_prev), st_new)
                assert st_d["kappa_c"] == st_new["kappa_c"]
                assert np.abs(C_re - C_desc).max() <= 1e-4 * np.abs(C_desc).max()
                n += 1
            st, eps_prev = st_new, eps
    assert n > 200


def test_esfuerzo_no_finito_si_la_deformacion_no_lo_es():
    """Si el problema local no tiene solución el material devuelve un esfuerzo
    no finito sin lanzar (el solver repite el paso más corto, §8b, Q5)."""
    m = mat()
    st = cracked_state(m, 0.3, 0.1 * GF / FT)
    st.update(damage_c=0.5, kappa_c=2 * KCP, kappa_c_prev=2 * KCP)
    sig, C, _ = m.compute_state(np.array([np.nan, 0.0, 0.0]), st)
    assert not np.all(np.isfinite(sig))


# ---------------------------------------------------------------------------
# σ_zz en deformación plana, kernel por lotes, unidades
# ---------------------------------------------------------------------------

def test_sigma_zz_en_deformacion_plana_con_aplastamiento():
    """σ_zz = σ̄_zz − d⁻·⟨σ̄_zz⟩₋ con σ̄_zz = ν·(σ̄_xx + σ̄_yy) (§4), igual en la
    versión por lotes."""
    m = mat(PE)
    st = m.initial_point_state(H)
    sigs, states = [], []
    for k in range(1, 300):
        eps = np.array([-1.0, NU / (1 - NU), 0.0]) * k * 1e-5
        sig, _, st = m.compute_state(eps, st)
        sb = Ce(PE) @ eps
        szz = NU * (sb[0] + sb[1])
        ref = szz if szz > 0 else (1 - st["damage_c"]) * szz
        assert m.out_of_plane_stress(sig, st) == pytest.approx(ref, rel=1e-10, abs=1e-6)
        sigs.append(sig)
        states.append([st[n] for n in m.STATE_SCHEMA])
    assert st["damage_c"] > 0.3
    S = np.array(states)
    batch = m.batch_out_of_plane_stress(np.array(sigs), S)
    scalar = [m.out_of_plane_stress(s, dict(zip(m.STATE_SCHEMA, row)))
              for s, row in zip(sigs, S)]
    np.testing.assert_allclose(batch, scalar, rtol=1e-14, atol=1e-6)


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_kernel_por_lotes_reproduce_compute_state_con_aplastamiento(hyp, integration):
    """batch_kernel_matches con aplastamiento: σ, tangente y fila de estado bit a bit."""
    rng = np.random.default_rng(23)
    m = mat(hyp, integration=integration)
    m.set_pseudo_time_ratio(1.0)
    names = list(m.STATE_SCHEMA)
    kern, params, Cm = m.batch_kernel(), m.batch_params(), m.batch_matrix()
    n = 0
    for _ in range(4):
        st = m.initial_point_state(H)
        for eps in _camino_aleatorio(rng, n_tramos=5, n_pasos=30, amp=3e-3, comp=1.0):
            sig, C, st_new = m.compute_state(eps, st)
            S_old = np.array([st[k] for k in names])
            S_new = np.zeros(len(names))
            sig_b = np.zeros(3)
            flag = np.zeros(1, np.int8)
            C_b = kern(eps, S_old, S_new, params, Cm, sig_b, np.zeros((3, 3)), flag)
            assert np.array_equal(sig_b, sig) and np.array_equal(C_b, C)
            assert np.array_equal(S_new, np.array([st_new[k] for k in names]))
            n += st_new["damage_c"] > 0.0
            st = st_new
    assert n > 50


@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_elemento_con_grieta_abierta_y_aplastamiento_igual_al_punto_material(integration, caplog):
    """Camino acoplado a nivel de elemento (§8b): un Quad4 de 100 mm con
    tracción en x y compresión en y proporcionales (deformación uniforme):
    la grieta se abre en x y el concreto entre grietas se aplasta en y. El
    solver termina sin bisecciones, con pocas iteraciones por paso
    (tangente consistente), y el esfuerzo es el del punto material."""
    import logging
    import re

    from solidum.elements.solid_2d import Quad4
    from solidum.math.assembly import Assembler
    from solidum.math.solvers import NonlinearSolver
    from solidum.mesh import structured_mesh

    L_el = 0.1
    ex, ey = 6.0e-4, -2.5e-3
    m = mat(integration=integration)
    mesh = structured_mesh(Quad4, (1, 1), size=(L_el, L_el), material=m, thickness=1.0)
    dom = mesh.domain
    for nd in mesh.nodes_on("xmin"):
        nd.fix_dof("ux", 0.0)
    for nd in mesh.nodes_on("xmax"):
        nd.fix_dof("ux", ex * L_el)
    for nd in mesh.nodes_on("ymin"):
        nd.fix_dof("uy", 0.0)
    for nd in mesh.nodes_on("ymax"):
        nd.fix_dof("uy", ey * L_el)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    with caplog.at_level(logging.INFO, logger="solidum"):
        U = NonlinearSolver(asm, num_steps=200).solve(np.zeros(dom.total_dofs))
    iters = [int(k) for k in re.findall(r"\| (\d+) iter", caplog.text)]
    assert "bisección" not in caplog.text and len(iters) >= 200
    assert max(iters) <= (6 if integration == "implicit" else 8)
    # punto material con el mismo camino (h del Quad4 de 100 mm: 100 mm)
    el = next(iter(dom.elements.values()))
    gp = el.state_vars[0]
    m.set_pseudo_time_ratio(1.0)
    st = m.initial_point_state(gp["h"])
    for lam in np.linspace(0.0, 1.0, 201)[1:]:
        sig, _, st = m.compute_state(np.array([ex, ey, 0.0]) * lam, st)
    assert st["cracked"] == 1.0 and st["w_n"] > 0.0 and st["damage_c"] > 0.3
    for k in ("w_max_n", "kappa_c", "damage_c"):
        assert gp[k] == pytest.approx(st[k], rel=1e-8), k
    # reacciones: deformación uniforme, σ_xx·L y σ_yy·L por unidad de espesor
    F = asm.assemble_non_linear_system(U)[1]
    Fx = sum(F[nd.dofs["ux"]] for nd in mesh.nodes_on("xmax"))
    Fy = sum(F[nd.dofs["uy"]] for nd in mesh.nodes_on("ymax"))
    assert Fx == pytest.approx(sig[0] * L_el, rel=1e-8, abs=1e-8 * FT * L_el)
    assert Fy == pytest.approx(sig[1] * L_el, rel=1e-8)


def test_invariancia_de_unidades_con_aplastamiento():
    """units_invariance: el mismo camino en Pa–m y en MPa–mm, con todos los
    parámetros físicos ≠ 1, da los mismos esfuerzos, aperturas y daños."""
    rng = np.random.default_rng(29)
    a = FixedSmearedCrack2D(E, NU, FT, GF, f_c=FC, G_c=GC, beta=1.2)
    b = FixedSmearedCrack2D(E * 1e-6, NU, FT * 1e-6, GF * 1e-3, f_c=FC * 1e-6, G_c=GC * 1e-3,
                            beta=1.2)
    sa, sb = a.initial_point_state(H), b.initial_point_state(H * 1e3)
    for eps in _camino_aleatorio(rng, n_tramos=6, n_pasos=40, amp=3e-3, comp=1.0):
        siga, _, sa = a.compute_state(eps, sa)
        sigb, _, sb = b.compute_state(eps, sb)
        np.testing.assert_allclose(sigb * 1e6, siga, rtol=1e-12, atol=1e-12 * FC)
        for k in ("w_n", "w_t", "slip", "w_max_n", "w_max_t"):
            assert sb[k] * 1e-3 == pytest.approx(sa[k], rel=1e-12, abs=1e-18), k
        for k in ("damage_n", "damage_t", "damage_c", "kappa_c"):
            assert sb[k] == pytest.approx(sa[k], rel=1e-12, abs=1e-15), k
    assert sa["damage_c"] > 0.0
