"""Ley de ablandamiento de Hordijk en la apertura de grieta (bloque 3b).

Spec ``IsotropicDamageCrackBand2D`` §2b y casos ``hordijk`` de su contrato,
para el modelo 1 (``IsotropicDamageCrackBand2D/3D``) y el daño de tracción del
modelo 2 (``TensionCompressionDamage2D/3D``) con ``softening='hordijk'``:

    σ/f_t = g(w/w_c),  g(x) = (1 + c₁³x³)·e^{−c₂x} − x·(1 + c₁³)·e^{−c₂},
    c₁ = 3, c₂ = 6.93, w_c = k_w·G_f/f_t con k_w = 1/∫₀¹g = 5.13606 (Q12),
    w = h·d·κ (banda de grieta en la apertura), h < h_max = 0.7382·ℓ_ch (Q14).

Las referencias de las pruebas (``g``, su primitiva, las constantes) se
escriben aquí en forma cerrada, sin usar las funciones del programa.

Unidades: N, mm, MPa, salvo la prueba de escalado (también SI).
"""
import math

import numpy as np
import pytest
from scipy.optimize import brentq

import solidum  # noqa: F401
from solidum.materials._softening import (
    HORDIJK_H_MAX_FACTOR,
    HORDIJK_KW,
    HORDIJK_SLOPE0,
    hordijk_shape,
    hordijk_shape_slope,
)
from solidum.materials.damage_crack_band_2d import IsotropicDamageCrackBand2D
from solidum.materials.damage_crack_band_3d import IsotropicDamageCrackBand3D
from solidum.materials.damage_tension_compression_2d import TensionCompressionDamage2D
from solidum.materials.damage_tension_compression_3d import TensionCompressionDamage3D

C1, C2 = 3.0, 6.93
# Valores de referencia de la revisión de teorico-solidum (mpmath, 40 dígitos;
# informes/teorico/verificacion_teorica_2026-11-24_hordijk_spec.md).
KW_REF = 5.1360552952511618
SLOPE0_REF = 6.957384024315071

NU = 0.2
SETS = {"concreto": (30000.0, 3.0, 0.1), "otro": (25000.0, 4.5, 0.06)}   # E, f_t, G_f
MODELS = ("CB2D", "CB3D", "TC2D", "TC3D")
DAMAGE_KEY = {"CB2D": "damage", "CB3D": "damage", "TC2D": "damage_t", "TC3D": "damage_t"}


def g_ref(x):
    return (1.0 + (C1 * x) ** 3) * math.exp(-C2 * x) - x * (1.0 + C1 ** 3) * math.exp(-C2)


def g_area(x):
    """∫₀ˣ g, en forma cerrada."""
    c = C2
    e = math.exp(-c * x)
    cubic = 6.0 / c ** 4 - e * (x ** 3 / c + 3.0 * x ** 2 / c ** 2 + 6.0 * x / c ** 3 + 6.0 / c ** 4)
    return (1.0 - e) / c + C1 ** 3 * cubic - 0.5 * x * x * (1.0 + C1 ** 3) * math.exp(-c)


def make(model, E, ft, Gf, nu=NU, **kw):
    kw.setdefault("softening", "hordijk")
    cls = {"CB2D": IsotropicDamageCrackBand2D, "CB3D": IsotropicDamageCrackBand3D,
           "TC2D": TensionCompressionDamage2D, "TC3D": TensionCompressionDamage3D}[model]
    return cls(E, nu, ft, Gf, **kw)


def uniaxial(model, e, nu=NU):
    """Tracción uniaxial de esfuerzo (esfuerzo plano en 2D)."""
    if model.endswith("2D"):
        return np.array([e, -nu * e, 0.0])
    return np.array([e, -nu * e, -nu * e, 0.0, 0.0, 0.0])


def h_max(E, ft, Gf):
    return HORDIJK_H_MAX_FACTOR * E * Gf / ft ** 2


# ---------------------------------------------------------------------------
# Constantes y parámetro
# ---------------------------------------------------------------------------

def test_constantes_de_hordijk_coinciden_con_la_referencia():
    assert HORDIJK_KW == pytest.approx(KW_REF, rel=1e-12)
    assert HORDIJK_SLOPE0 == pytest.approx(SLOPE0_REF, rel=1e-12)
    assert HORDIJK_H_MAX_FACTOR == pytest.approx(KW_REF / SLOPE0_REF, rel=1e-12)
    assert HORDIJK_H_MAX_FACTOR == pytest.approx(0.7382, abs=1e-4)
    assert hordijk_shape(0.0) == 1.0
    assert abs(hordijk_shape(1.0)) < 1e-16
    for x in np.linspace(0.0, 1.0, 11):
        assert hordijk_shape(x) == pytest.approx(g_ref(x), rel=1e-14, abs=1e-16)
    assert hordijk_shape_slope(0.0) == pytest.approx(-SLOPE0_REF, rel=1e-14)
    # k_w encierra el área unitaria: w_c = k_w·G_f/f_t disipa exactamente G_f.
    assert KW_REF * g_area(1.0) == pytest.approx(1.0, rel=1e-14)


@pytest.mark.parametrize("model", MODELS)
def test_ley_por_omision_es_la_exponencial(model):
    """Q13: sin ``softening`` la ley es la exponencial y nada cambia."""
    E, ft, Gf = SETS["concreto"]
    m = make(model, E, ft, Gf, softening="exponential")
    assert make(model, E, ft, Gf, softening="exponential").softening == "exponential"
    cls = type(m)
    default = cls(E, NU, ft, Gf)
    assert default.softening == "exponential"
    assert default.characteristic_length_limit == pytest.approx(2.0 * E * Gf / ft ** 2, rel=1e-15)
    st_a, st_b = m.initial_point_state(10.0), default.initial_point_state(10.0)
    for e in np.linspace(0.0, 40.0, 41)[1:] * ft / E:
        sa, Ca, st_a = m.compute_state(uniaxial(model, e), st_a)
        sb, Cb, st_b = default.compute_state(uniaxial(model, e), st_b)
        np.testing.assert_array_equal(sa, sb)
        np.testing.assert_array_equal(Ca, Cb)
    np.testing.assert_array_equal(m.batch_params(), default.batch_params())


@pytest.mark.parametrize("model", MODELS)
def test_softening_no_valida_se_rechaza(model):
    E, ft, Gf = SETS["concreto"]
    with pytest.raises(ValueError, match="softening"):
        make(model, E, ft, Gf, softening="lineal")


# ---------------------------------------------------------------------------
# Curva uniaxial, energía y límite de h
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("frac", [0.02, 0.1, 0.4, 0.9])
@pytest.mark.parametrize("set_name", list(SETS))
@pytest.mark.parametrize("model", MODELS)
def test_curva_uniaxial_es_la_de_hordijk(model, set_name, frac):
    """σ = f_t·g(w/w_c) con w = h·(ε − σ/E) en todo el ablandamiento
    (absoluto respecto a f_t), pico f_t en κ_0 y tope desde κ = w_c/h."""
    E, ft, Gf = SETS[set_name]
    h = frac * h_max(E, ft, Gf)
    w_c = KW_REF * Gf / ft
    k0 = ft / E
    m = make(model, E, ft, Gf)
    st = m.initial_point_state(h)
    key = DAMAGE_KEY[model]
    eps_end = 1.3 * w_c / h
    path = np.r_[np.linspace(0.0, k0, 5)[1:],
                 k0 + (w_c / h - k0) * np.linspace(0.0, 1.0, 400)[1:-1] ** 2,
                 np.linspace(w_c / h, eps_end, 5)]
    peak = 0.0
    for e in path:
        sig, _, st = m.compute_state(uniaxial(model, e), st)
        s = sig[0]
        peak = max(peak, s)
        if e <= k0:
            assert st[key] == 0.0
            assert s == pytest.approx(E * e, rel=1e-12)
        elif e < w_c / h:
            w = h * (e - s / E)
            assert abs(s - ft * g_ref(w / w_c)) <= 1e-12 * ft, (e, s)
            assert w == pytest.approx(h * st[key] * e, rel=1e-12, abs=1e-15 * w_c)
        else:
            assert st[key] == m.damage_max
    assert peak == pytest.approx(ft, rel=1e-12)


def _w_tope(E, ft, Gf, h, cap):
    """Apertura en que d = w/(h·κ(w)) llega al tope, κ(w) = w/h + κ_0·g(w/w_c)."""
    w_c = KW_REF * Gf / ft
    k0 = ft / E
    return brentq(lambda w: w / (h * (w / h + k0 * g_ref(w / w_c))) - cap,
                  0.5 * w_c, w_c, xtol=1e-16 * w_c, rtol=1e-15, maxiter=200)


@pytest.mark.parametrize("frac", [0.02, 0.4])
@pytest.mark.parametrize("model", ["CB2D", "TC2D"])
def test_energia_disipada_es_la_de_la_ley(model, frac):
    """∫σ dε − ½σε = D(w) = (1/h)(∫₀ʷσ dw − ½σw) en estados intermedios y
    D(w_tope) al final (la pérdida del tope es exacta en esa forma)."""
    E, ft, Gf = SETS["concreto"]
    h = frac * h_max(E, ft, Gf)
    w_c = KW_REF * Gf / ft
    k0 = ft / E
    m = make(model, E, ft, Gf)
    st = m.initial_point_state(h)
    eps = np.r_[np.linspace(0.0, k0, 50), k0 + (1.5 * w_c / h - k0) * np.linspace(0.0, 1.0, 400001)[1:] ** 2]
    sig = np.empty_like(eps)
    for i, e in enumerate(eps):
        s, _, st = m.compute_state(uniaxial(model, e), st)
        sig[i] = s[0]
    work = np.concatenate([[0.0], np.cumsum(0.5 * (sig[1:] + sig[:-1]) * np.diff(eps))])
    D_num = work - 0.5 * sig * eps
    for i in (len(eps) // 40, len(eps) // 4, len(eps) // 2):
        w = h * (eps[i] - sig[i] / E)
        if w >= w_c:
            continue
        D_ref = (ft * w_c * g_area(w / w_c) - 0.5 * sig[i] * w) / h
        assert D_num[i] == pytest.approx(D_ref, rel=1e-6), i
    w_t = _w_tope(E, ft, Gf, h, m.damage_max)
    s_t = ft * g_ref(w_t / w_c)
    D_final = (ft * w_c * g_area(w_t / w_c) - 0.5 * s_t * w_t) / h
    assert D_num[-1] == pytest.approx(D_final, rel=1e-6)
    assert D_final * h / Gf == pytest.approx(1.0, abs=2e-3)


@pytest.mark.parametrize("model", MODELS)
def test_h_sobre_el_limite_de_hordijk_se_rechaza(model):
    E, ft, Gf = SETS["concreto"]
    m = make(model, E, ft, Gf)
    hm = h_max(E, ft, Gf)
    assert m.characteristic_length_limit == pytest.approx(hm, rel=1e-15)
    for h in (hm, 1.01 * hm):
        with pytest.raises(ValueError, match=r"Hordijk.*h_max.*ℓ_ch"):
            m.initial_point_state(h)


@pytest.mark.parametrize("model", MODELS)
def test_cerca_del_limite_el_dano_es_continuo_y_creciente(model):
    """Con h = 0.999·h_max el daño parte de 0 en κ_0 (sin salto en el pico),
    crece con κ y el Newton local converge."""
    E, ft, Gf = SETS["concreto"]
    h = 0.999 * h_max(E, ft, Gf)
    w_c = KW_REF * Gf / ft
    k0 = ft / E
    m = make(model, E, ft, Gf)
    st = m.initial_point_state(h)
    key = DAMAGE_KEY[model]
    d_prev = 0.0
    for r in (1.0 + 1e-9, 1.0 + 1e-6, 1.0 + 1e-3, 1.01, 1.1, 1.5, 2.0, 3.0):
        e = r * k0
        sig, _, st = m.compute_state(uniaxial(model, e), st)
        d = st[key]
        assert d > d_prev
        w = h * (e - sig[0] / E)
        assert abs(sig[0] - ft * g_ref(w / w_c)) <= 1e-10 * ft
        d_prev = d
    _, _, st0 = m.compute_state(uniaxial(model, (1.0 + 1e-9) * k0), m.initial_point_state(h))
    assert st0[key] < 1e-5


@pytest.mark.parametrize("model", ["CB2D", "TC2D"])
def test_dano_depende_solo_de_kappa_sobre_kappa0_y_h_sobre_lch(model):
    """Dos materiales con E, f_t, G_f y unidades distintas (N–mm y SI) e igual
    h/ℓ_ch dan la misma d(κ/κ_0): detecta un factor olvidado en w_c ∝ G_f/f_t
    o en h."""
    mats = [(25000.0, 4.5, 0.06), (30.0e9, 3.0e6, 100.0)]
    curves = []
    for E, ft, Gf in mats:
        m = make(model, E, ft, Gf)
        st = m.initial_point_state(0.3 * E * Gf / ft ** 2)
        d = []
        for r in np.linspace(0.5, 30.0, 120):
            _, _, st = m.compute_state(uniaxial(model, r * ft / E), st)
            d.append(st[DAMAGE_KEY[model]])
        curves.append(np.array(d))
    np.testing.assert_allclose(curves[0], curves[1], rtol=0.0, atol=1e-12)
    assert curves[0].max() > 0.99


# ---------------------------------------------------------------------------
# Tangente consistente por diferencias finitas
# ---------------------------------------------------------------------------

def _fd_check(m, eps, st, step=5e-7, tol=5e-6):
    _, C, st_new = m.compute_state(eps, st)
    d = st_new.get("damage", st_new.get("damage_t"))
    assert 0.0 < d < m.damage_max, "el estado debe estar en la rama consistente"
    n = len(eps)
    C_fd = np.zeros((n, n))
    for j in range(n):
        d = step * np.abs(eps).max()
        ep, em = eps.copy(), eps.copy()
        ep[j] += d
        em[j] -= d
        C_fd[:, j] = (m.compute_state(ep, st)[0] - m.compute_state(em, st)[0]) / (2 * d)
    assert np.abs(C - C_fd).max() <= tol * np.abs(C_fd).max()


E0, FT0, GF0 = SETS["concreto"]
H0 = 25.0        # ≈ 0.1·h_max
STRAINS_2D = (np.array([3.0e-4, -0.4e-4, 1.1e-4]), np.array([2.0e-4, 1.5e-4, -0.7e-4]))
STRAINS_3D = (np.array([3.0e-4, -0.4e-4, 0.3e-4, 1.1e-4, -0.5e-4, 0.8e-4]),
              np.array([1.0e-4, 2.2e-4, -0.6e-4, -0.9e-4, 0.4e-4, 0.2e-4]))


@pytest.mark.parametrize("measure", ["rankine", "energy", "mazars", "modified_von_mises"])
@pytest.mark.parametrize("hypothesis", ["plane_stress", "plane_strain"])
def test_tangente_consistente_modelo1_2d(measure, hypothesis):
    kw = {"k": 10.0} if measure == "modified_von_mises" else {}
    m = make("CB2D", E0, FT0, GF0, hypothesis=hypothesis, equivalent_strain=measure, **kw)
    st = m.initial_point_state(H0)
    for eps in STRAINS_2D:
        _fd_check(m, eps, st)


@pytest.mark.parametrize("measure", ["rankine", "energy", "mazars", "modified_von_mises"])
def test_tangente_consistente_modelo1_3d(measure):
    kw = {"k": 10.0} if measure == "modified_von_mises" else {}
    m = make("CB3D", E0, FT0, GF0, equivalent_strain=measure, **kw)
    st = m.initial_point_state(H0)
    for eps in STRAINS_3D:
        _fd_check(m, eps, st)


@pytest.mark.parametrize("hypothesis", ["plane_stress", "plane_strain"])
def test_tangente_consistente_modelo2_2d(hypothesis):
    m = make("TC2D", E0, FT0, GF0, hypothesis=hypothesis, f_c=30.0, G_c=15.0)
    st = m.initial_point_state(H0)
    for eps in STRAINS_2D:
        _fd_check(m, eps, st)


def test_tangente_consistente_modelo2_3d():
    m = make("TC3D", E0, FT0, GF0, f_c=30.0, G_c=15.0)
    st = m.initial_point_state(H0)
    for eps in STRAINS_3D:
        _fd_check(m, eps, st)


def test_tangente_uniaxial_es_la_de_dos_resortes_en_serie():
    """dσ/dε = E·H_w/(E + H_w) con H_w = h·f_t·g'(w/w_c)/w_c (spec §2b), con
    ν = 0 para que la tangente uniaxial sea la componente xx."""
    m = make("CB2D", E0, FT0, GF0, nu=0.0)
    w_c = KW_REF * GF0 / FT0
    st = m.initial_point_state(H0)
    for frac in (0.05, 0.3, 0.7):
        e = FT0 / E0 + frac * (w_c / H0 - FT0 / E0)
        sig, C, st_new = m.compute_state(np.array([e, 0.0, 0.0]), st)
        w = H0 * (e - sig[0] / E0)
        x = w / w_c
        gp = 3 * C1 ** 3 * x * x * math.exp(-C2 * x) - C2 * (1 + (C1 * x) ** 3) * math.exp(-C2 * x) \
            - (1 + C1 ** 3) * math.exp(-C2)
        Hw = H0 * FT0 * gp / w_c
        assert C[0, 0] == pytest.approx(E0 * Hw / (E0 + Hw), rel=1e-9)
        st = st_new


# ---------------------------------------------------------------------------
# Descarga y recarga, por lotes, 3D y modelo 2
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("model", MODELS)
def test_descarga_y_recarga_por_la_secante(model):
    E, ft, Gf = SETS["concreto"]
    m = make(model, E, ft, Gf)
    key = DAMAGE_KEY[model]
    st = m.initial_point_state(H0)
    e1 = 8.0 * ft / E
    for e in np.linspace(0.0, e1, 30)[1:]:
        _, _, st = m.compute_state(uniaxial(model, e), st)
    d1, st1 = st[key], dict(st)
    assert 0.0 < d1 < 1.0
    for e in np.r_[np.linspace(e1, 0.3 * e1, 8)[1:], np.linspace(0.3 * e1, e1, 8)[1:]]:
        sig, _, st = m.compute_state(uniaxial(model, e), st)
        assert st[key] == d1
        assert sig[0] == pytest.approx((1.0 - d1) * E * e, rel=1e-12)
        if model.startswith("CB"):
            # Daño fijo: la tangente es la secante (1 − d)·C_e.
            _, C_un, _ = m.compute_state(uniaxial(model, e), st)
            np.testing.assert_allclose(C_un, (1.0 - d1) * m.elastic_base.C, rtol=1e-14, atol=0.0)
    # El mismo intento dos veces desde el mismo estado da lo mismo (trial/commit).
    a = m.compute_state(uniaxial(model, 1.2 * e1), st1)
    b = m.compute_state(uniaxial(model, 1.2 * e1), st1)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])
    assert a[2] == b[2]
    assert a[2][key] > d1


@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
@pytest.mark.parametrize("model", MODELS)
def test_kernel_por_lotes_reproduce_compute_state(model, integration):
    kw = {"f_c": 30.0, "G_c": 15.0} if model.startswith("TC") else {}
    m = make(model, E0, FT0, GF0, integration=integration, **kw)
    n = m.STRAIN_DIM
    names = list(m.STATE_SCHEMA)
    kernel = m.batch_kernel()
    C = m.batch_matrix()
    st = m.initial_point_state(H0)
    rng = np.random.default_rng(11)
    for i in range(25):
        m.set_pseudo_time_ratio(None if i % 5 == 4 else 0.5 + 0.1 * i)
        params = m.batch_params()
        eps = (rng.standard_normal(n) + np.r_[2.0, np.zeros(n - 1)]) * 1.5e-4 * (1 + 0.2 * i)
        sig_ref, C_ref, st_new = m.compute_state(eps, st)
        S_old = np.array([st[k] for k in names])
        S_new, sig, C_out, flag = (np.zeros(len(names)), np.zeros(n), np.zeros((n, n)),
                                   np.zeros(1, np.int64))
        C_t = kernel(eps, S_old, S_new, params, C, sig, C_out, flag)
        np.testing.assert_array_equal(sig, sig_ref)
        np.testing.assert_array_equal(C_t, C_ref)
        np.testing.assert_array_equal(S_new, [st_new[k] for k in names])
        assert flag[0] == 0
        st = st_new
    assert st[DAMAGE_KEY[model]] > 0.5


def _assoc_3d(eps2, nu, hypothesis):
    c_zz = -nu / (1.0 - nu) if hypothesis == "plane_stress" else 0.0
    return np.array([eps2[0], eps2[1], c_zz * (eps2[0] + eps2[1]), eps2[2], 0.0, 0.0])


@pytest.mark.parametrize("hypothesis", ["plane_stress", "plane_strain"])
@pytest.mark.parametrize("pair", [("CB2D", "CB3D"), ("TC2D", "TC3D")])
def test_restriccion_exacta_del_3d(pair, hypothesis):
    m2 = make(pair[0], E0, FT0, GF0, hypothesis=hypothesis)
    m3 = make(pair[1], E0, FT0, GF0)
    key = DAMAGE_KEY[pair[0]]
    s2, s3 = m2.initial_point_state(H0), m3.initial_point_state(H0)
    base = np.array([1.0, -0.15, 0.35])
    for f in (0.5, 1.5, 3.0, 2.0, 1.0, 4.0, 6.0, 20.0):       # carga, descarga y recarga
        e2 = f * FT0 / E0 * base
        sig2, _, s2 = m2.compute_state(e2, s2)
        sig3, _, s3 = m3.compute_state(_assoc_3d(e2, NU, hypothesis), s3)
        assert s2[key] == pytest.approx(s3[key], rel=1e-12, abs=1e-15)
        np.testing.assert_allclose(sig2, sig3[[0, 1, 3]], rtol=1e-10, atol=1e-12)
    assert s2[key] > 0.5


@pytest.mark.parametrize("dim", ["2D", "3D"])
def test_modelo2_en_traccion_uniaxial_reproduce_al_modelo1(dim):
    m1 = make("CB" + dim, E0, FT0, GF0)
    m2 = make("TC" + dim, E0, FT0, GF0)
    s1, s2 = m1.initial_point_state(H0), m2.initial_point_state(H0)
    for e in np.r_[np.linspace(0.0, 30.0, 61)[1:], np.linspace(30.0, 10.0, 5), 80.0] * FT0 / E0:
        sig1, _, s1 = m1.compute_state(uniaxial("CB" + dim, e), s1)
        sig2, _, s2 = m2.compute_state(uniaxial("TC" + dim, e), s2)
        assert s2["damage_t"] == pytest.approx(s1["damage"], rel=1e-12, abs=1e-15)
        assert sig2[0] == pytest.approx(sig1[0], rel=1e-12, abs=1e-12)


# ---------------------------------------------------------------------------
# IMPL-EX
# ---------------------------------------------------------------------------

def test_impl_ex_secante_exacta_y_disipacion_no_negativa():
    m = make("CB2D", E0, FT0, GF0, integration="impl-ex")
    st = m.initial_point_state(H0)
    Ce = m.elastic_base.C
    k0 = FT0 / E0
    camino = np.r_[np.linspace(0.0, 8.0, 41)[1:], np.linspace(8.0, 3.0, 11)[1:],
                   np.linspace(3.0, 40.0, 41)[1:]] * k0
    d_prev = 0.0
    for k, e in enumerate(camino):
        m.set_pseudo_time_ratio(2.0 if k % 7 == 0 else 1.0)
        eps = uniaxial("CB2D", e)
        sig, C, st = m.compute_state(eps, st)
        d = st["damage"]
        np.testing.assert_array_equal(C, (1.0 - d) * Ce)
        assert d >= d_prev
        assert 0.5 * eps @ Ce @ eps * (d - d_prev) >= 0.0
        d_prev = d
    assert d_prev > 0.9


@pytest.mark.parametrize("model", MODELS)
def test_impl_ex_extrapolacion_mas_alla_de_w_c_da_el_tope(model):
    """κ̃ ≥ w_c/h: el daño extrapolado es el tope (sin entrar al Newton,
    ``_softening.evaluate_hordijk_damage``), en los cuatro materiales."""
    m = make(model, E0, FT0, GF0, integration="impl-ex")
    w_c = KW_REF * GF0 / FT0
    key = DAMAGE_KEY[model]
    st = m.initial_point_state(H0)
    if model.startswith("CB"):
        st.update(kappa=0.7 * w_c / H0, kappa_prev=0.3 * w_c / H0, damage=0.99)
    else:
        st.update(kappa_t=0.7 * w_c / H0, kappa_t_prev=0.3 * w_c / H0, damage_t=0.99)
    m.set_pseudo_time_ratio(1.0)
    _, _, st_new = m.compute_state(uniaxial(model, 0.71 * w_c / H0), st)
    assert st_new[key] == m.damage_max
    assert st_new["damage_error"] > 0.0


def _camino(n, integration, kappa_de_t):
    m = make("CB2D", E0, FT0, GF0, integration=integration)
    st = m.initial_point_state(H0)
    m.set_pseudo_time_ratio(1.0)
    out = []
    for t in np.linspace(0.0, 1.0, n + 1)[1:]:
        sig, _, st = m.compute_state(uniaxial("CB2D", kappa_de_t(t)), st)
        out.append(sig[0])
    return np.array(out)


@pytest.mark.parametrize("forma, orden", [("lineal", 1.0), ("cuadratica", 2.0)])
def test_impl_ex_converge_a_la_implicita(forma, orden):
    """Primer orden con el quiebre del inicio del daño y segundo con κ(t) suave,
    como la exponencial (spec §8)."""
    k0 = FT0 / E0
    if forma == "lineal":
        def kappa(t):
            return k0 * (1.0 + 39.0 * t)
    else:
        def kappa(t):
            return k0 * (1.0 + 39.0 * t * t)
    errores = []
    for n in (100, 200, 400, 800):
        errores.append(np.abs(_camino(n, "impl-ex", kappa) - _camino(n, "implicit", kappa)).max())
    razones = np.array(errores[:-1]) / np.array(errores[1:])
    np.testing.assert_allclose(razones, 2.0 ** orden, rtol=0.15)
