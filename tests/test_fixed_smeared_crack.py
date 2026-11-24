"""``FixedSmearedCrack2D`` (modelo 3), bloque 1: punto material, IMPL-EX y tangente global.

Aceptación ``material_point`` e ``impl_ex`` de la spec ``FixedSmearedCrack2D``
(§1–§3, §6–§9) sin aplastamiento (el bloque 2, §4, no está programado). Cada
prueba lleva en su nombre el invariante que protege.

Las referencias se escriben aquí, aparte del código: la matriz ``N`` de la
spec §1 (el código usa su inversa ``T``), la ley de ablandamiento (la de
Hordijk con ``w_c`` obtenida integrando ``g`` numéricamente), la flexibilidad
con el tope, la energía disipada hasta el tope por cuadratura y, en el
problema local, la comprobación de la ley de cada dirección sobre la solución
(con una sola solución bajo el límite de ``h``, residuo nulo = la solución).

Los estados se construyen de dos maneras: (a) con la deformación que
corresponde, según la spec, a un estado de esfuerzo y de aperturas dado
(``strain_of``), y (b) imponiendo la historia directamente (``cracked_state``).
La historia fija (ley secante) se impone con IMPL-EX fuera de un paso.

Concreto usual en SI: E = 30 GPa, ν = 0.2, f_t = 3 MPa, G_f = 100 N/m
(ℓ_ch = 333 mm, límite de h 278 mm con la exponencial y 205 mm con Hordijk).
"""
import math

import numpy as np
import pytest
from scipy import integrate, optimize

import solidum  # noqa: F401
from solidum.constants import CRACK_BAND_DAMAGE_MAX
from solidum.elements.solid_2d import Quad4
from solidum.materials.fixed_smeared_crack_2d import FixedSmearedCrack2D

E, NU, FT, GF = 30.0e9, 0.2, 3.0e6, 100.0
H = 0.025
G = E / (2.0 * (1.0 + NU))
ELL = E * GF / FT**2
CAP = CRACK_BAND_DAMAGE_MAX
K0 = FT / E
PS, PE = "plane_stress", "plane_strain"
LAWS = ["exponential", "hordijk"]

# Ley de Hordijk (Cornelissen, Hordijk y Reinhardt 1986), escrita aquí aparte.
C1, C2 = 3.0, 6.93


def g_hordijk(x):
    return (1.0 + C1**3 * x**3) * math.exp(-C2 * x) - x * (1.0 + C1**3) * math.exp(-C2)


KW = 1.0 / integrate.quad(g_hordijk, 0.0, 1.0, epsabs=0.0, epsrel=1e-13)[0]
SLOPE0 = C2 + (1.0 + C1**3) * math.exp(-C2)            # |g'(0)|


def w_crit(ft=FT, gf=GF):
    return KW * gf / ft


def f_law(w, law="exponential", ft=FT, gf=GF):
    """Ley de ablandamiento de la spec §2, ``f(w)``."""
    if law == "hordijk":
        wc = w_crit(ft, gf)
        return 0.0 if w >= wc else ft * g_hordijk(w / wc)
    return ft * math.exp(-ft * w / gf)


def c_max(h, E_=E):
    return h * CAP / ((1.0 - CAP) * E_)


def compliance(w, h, law="exponential", ft=FT, gf=GF, E_=E):
    """``c = w/f(w)`` con el tope de §6."""
    if w <= 0.0:
        return 0.0
    f = f_law(w, law, ft, gf)
    cm = c_max(h, E_)
    return cm if w >= cm * f else w / f


def env(w, h, law="exponential", ft=FT, gf=GF, E_=E):
    """Envolvente con el tope, ``max(f(w), w/c_max)``."""
    return max(f_law(w, law, ft, gf), w / c_max(h, E_))


def omega(c, h, E_=E):
    return E_ * c / (h + E_ * c)


def w_tope(h, law="exponential", ft=FT, gf=GF, E_=E):
    """Raíz de ``w = c_max·f(w)`` (spec §6)."""
    hi = w_crit(ft, gf) if law == "hordijk" else 60.0 * gf / ft
    return optimize.brentq(lambda w: w - c_max(h, E_) * f_law(w, law, ft, gf), 0.0, hi,
                           xtol=1e-300, rtol=4 * np.finfo(float).eps)


def d_tope(h, law="exponential"):
    """``D_tope = (1/h)·[∫₀^{w_tope} f dw − ½·f(w_tope)·w_tope]`` (spec §6)."""
    wt = w_tope(h, law)
    area = integrate.quad(lambda w: f_law(w, law), 0.0, wt, epsabs=0.0, epsrel=1e-13,
                          limit=200)[0]
    return (area - 0.5 * f_law(wt, law) * wt) / h


def N_spec(theta):
    """Matriz ``N`` de la spec §1 (de los ejes de la grieta a los globales)."""
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c * c, s * s, -c * s], [s * s, c * c, c * s],
                     [2 * c * s, -2 * c * s, c * c - s * s]])


def Ce(hyp=PS, E_=E, nu=NU):
    if hyp == PS:
        return E_ / (1 - nu**2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])
    k = E_ / ((1 + nu) * (1 - 2 * nu))
    return k * np.array([[1 - nu, nu, 0], [nu, 1 - nu, 0], [0, 0, (1 - 2 * nu) / 2]])


def strain_of(theta, s_loc, w_n=0.0, w_t=0.0, slip=0.0, h=H, hyp=PS):
    """Deformación global de un estado con esfuerzos locales ``[σ_nn, σ_tt, τ_nt]``
    del concreto y aperturas ``w_n``, ``w_t`` y deslizamiento ``slip`` (§1):
    ``ε = N·(C_e⁻¹·σ_loc + [w_n, w_t, s]/h)``."""
    e_loc = np.linalg.solve(Ce(hyp), np.asarray(s_loc, float)) + np.array([w_n, w_t, slip]) / h
    return N_spec(theta) @ e_loc


def stress_of(theta, s_loc):
    """``σ`` global a partir de la tracción en la grieta, ``t = Nᵀ·σ``."""
    return np.linalg.solve(N_spec(theta).T, np.asarray(s_loc, float))


def local(theta, sig):
    return N_spec(theta).T @ sig


def mat(hyp=PS, nu=NU, **kw):
    return FixedSmearedCrack2D(E, nu, FT, GF, hypothesis=hyp, **kw)


def cracked_state(m, theta, wn, wt=0.0, h=H, law="exponential"):
    """Historia impuesta: grieta en ``theta`` con aperturas máximas ``wn``, ``wt``
    (y ``damage_n/t`` coherentes, que IMPL-EX usa fuera de un paso)."""
    st = m.initial_point_state(h)
    st.update(cracked=1.0, theta=float(theta), w_max_n=float(wn), w_max_t=float(wt),
              w_max_n_prev=float(wn), w_max_t_prev=float(wt),
              damage_n=omega(compliance(wn, h, law), h), damage_t=omega(compliance(wt, h, law), h))
    return st


def principal(sig):
    return 0.5 * (sig[0] + sig[1]) + math.hypot(0.5 * (sig[0] - sig[1]), sig[2])


# ---------------------------------------------------------------------------
# Parámetros, contrato y límite de h (§2, H1)
# ---------------------------------------------------------------------------

ORDEN_DEL_ESTADO = ["cracked", "theta", "w_max_n", "w_max_t", "w_n", "w_t", "slip",
                    "damage_n", "damage_t", "kappa_c", "damage_c", "h", "w_max_n_prev",
                    "w_max_t_prev", "kappa_c_prev", "damage_error"]


def test_contrato_y_estado_inicial_como_en_la_spec():
    for hyp in (PS, PE):
        m = mat(hyp)
        assert m.uses_characteristic_length
        assert m.IS_SYMMETRIC is False and m.SECANT_UNLOADING is True
        assert m.PRIMARY_STATE_VAR == "damage_n" and m.STRAIN_DIM == 3
        assert m.damage_max == CAP and m.softening == "exponential"
        assert m.integration == "implicit" and not m.uses_pseudo_time
        assert m.hypothesis == hyp
        assert list(m.STATE_SCHEMA) == ORDEN_DEL_ESTADO
        st = m.initial_point_state(H)
        assert list(st) == ORDEN_DEL_ESTADO
        assert st == {k: (H if k == "h" else 0.0) for k in ORDEN_DEL_ESTADO}
    m = mat(integration="impl-ex")
    assert m.uses_pseudo_time and m.IS_SYMMETRIC is True      # secante (§9)


def test_parametros_invalidos_se_rechazan_al_construir():
    with pytest.raises(ValueError, match="nu"):
        mat(nu=-0.1)
    with pytest.raises(ValueError, match="nu"):
        mat(nu=0.5)
    with pytest.raises(ValueError, match="f_t"):
        FixedSmearedCrack2D(E, NU, 0.0, GF)
    with pytest.raises(ValueError, match="G_f"):
        FixedSmearedCrack2D(E, NU, FT, -1.0)
    with pytest.raises(ValueError, match="damage_max"):
        mat(damage_max=1.0)
    with pytest.raises(ValueError, match="integration"):
        mat(integration="explicit")
    with pytest.raises(ValueError, match="softening"):
        mat(softening="linear")
    with pytest.raises(ValueError):
        mat("plane_strian")
    for kw in (dict(f_c=30e6, G_c=15e3), dict(beta=1.16)):
        with pytest.raises(ValueError, match="bloque 2"):
            mat(**kw)


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("law", LAWS)
def test_limite_de_h_unicidad_con_dos_grietas(law, hyp):
    """snapback_limit (§2, H1): ``h < ℓ_ch/(1+ν)`` (exponencial) o
    ``0.738·ℓ_ch/(1+ν)`` (Hordijk), igual en las dos hipótesis."""
    m = mat(hyp, softening=law)
    lim = (KW / SLOPE0 if law == "hordijk" else 1.0) * ELL / (1.0 + NU)
    assert m.characteristic_length_limit == pytest.approx(lim, rel=1e-12)
    m.initial_point_state(0.999 * lim)
    with pytest.raises(ValueError, match="supera el límite"):
        m.initial_point_state(1.001 * lim)
    for h in (float("nan"), 0.0, -1.0):
        with pytest.raises(ValueError, match="no válida"):
            m.initial_point_state(h)


def test_grieta_sin_h_avisa_y_el_lote_la_marca():
    m = mat()
    m.compute_state(np.array([0.5 * K0, 0.0, 0.0]), None)        # elástico: admitido
    eps = np.array([2.0 * K0, 0.0, 0.0])
    with pytest.raises(ValueError, match="sin longitud característica"):
        m.compute_state(eps, None)
    names = list(m.STATE_SCHEMA)
    S_old = np.array([m.initial_state()[k] for k in names])
    flag = np.zeros(1, np.int8)
    m.batch_kernel()(eps, S_old, np.zeros(len(names)), m.batch_params(), m.batch_matrix(),
                     np.zeros(3), np.zeros((3, 3)), flag)
    assert flag[0] == 1
    with pytest.raises(ValueError, match="sin longitud característica"):
        m.batch_report(1)


# ---------------------------------------------------------------------------
# Antes de agrietarse, tracción uniaxial y deformación plana (§1, §2)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("hyp", [PS, PE])
def test_elastico_antes_de_agrietarse(hyp):
    m = mat(hyp)
    st = m.initial_point_state(H)
    rng = np.random.default_rng(5)
    C = Ce(hyp)
    checked = 0
    for _ in range(60):
        eps = rng.standard_normal(3) * 2e-4
        if principal(C @ eps) >= 0.99 * FT:
            eps *= 0.98 * FT / principal(C @ eps)
        sig, Ct, st1 = m.compute_state(eps, st)
        np.testing.assert_allclose(sig, C @ eps, rtol=1e-12, atol=1e-12 * FT)
        np.testing.assert_allclose(Ct, C, rtol=1e-12)
        assert st1 == st
        checked += 1
    assert checked == 60


@pytest.mark.parametrize("theta", [0.0, math.radians(30.0)], ids=["0", "30"])
@pytest.mark.parametrize("law", LAWS)
def test_traccion_uniaxial_sigue_la_ley_y_conserva_el_angulo(law, theta):
    """uniaxial_tension_curve: con ``ε = σ/E + w/h`` en ``n`` y la contracción
    lateral ``−ν·σ/E`` del concreto entre grietas (abrir no contrae), ``σ`` es
    uniaxial en ``n`` con ``σ = f(w)``; la grieta conserva el ángulo."""
    m = mat(softening=law)
    st = m.initial_point_state(H)
    wt = w_tope(H, law)
    path = [(s, 0.0) for s in np.linspace(0.0, FT, 11)]
    path += [(env(w, H, law), w) for w in np.linspace(0.0, 1.5 * wt, 151)[1:]]
    for s, w in path:
        sig, _, st = m.compute_state(strain_of(theta, [s, 0.0, 0.0], w_n=w), st)
        np.testing.assert_allclose(sig, stress_of(theta, [s, 0.0, 0.0]), rtol=0,
                                   atol=1e-10 * FT)
        if w > 0.0:
            assert st["theta"] == pytest.approx(theta, abs=1e-12)
            assert st["w_n"] == pytest.approx(w, rel=1e-10)
            assert st["w_max_n"] == pytest.approx(w, rel=1e-10)
            assert st["w_t"] == 0.0 and st["w_max_t"] == 0.0
    assert st["damage_n"] == pytest.approx(CAP, rel=1e-12)


@pytest.mark.parametrize("law", LAWS)
def test_deformacion_plana_uniaxial_misma_ley_que_esfuerzo_plano(law):
    """plane_strain_uniaxial: con ``σ_yy = 0`` y ``ε_zz = 0``,
    ``ε_xx = (1 − ν²)·σ/E + w/h``, ``σ_zz = ν·σ_xx`` y la misma ``σ–w``."""
    m = mat(PE, softening=law)
    st = m.initial_point_state(H)
    for w in np.r_[0.0, np.geomspace(1e-9, 2.0 * w_tope(H, law), 60)]:
        s = env(w, H, law) if w > 0 else 0.97 * FT
        eps = strain_of(0.0, [s, 0.0, 0.0], w_n=w, hyp=PE)
        assert eps[0] == pytest.approx((1 - NU**2) * s / E + w / H, rel=1e-12)
        sig, _, st = m.compute_state(eps, st)
        np.testing.assert_allclose(sig, [s, 0.0, 0.0], rtol=0, atol=1e-10 * FT)
        assert m.out_of_plane_stress(sig, st) == pytest.approx(NU * s, rel=1e-10, abs=1e-10 * FT)
        assert st["w_n"] == pytest.approx(w, rel=1e-10, abs=1e-30)


# ---------------------------------------------------------------------------
# Energía (§6, §7)
# ---------------------------------------------------------------------------

def _energia_traccion(law, h, hyp, gf_model=GF, n=20000):
    """Energía disipada en tracción uniaxial monótona hasta pasar el tope, con
    los esfuerzos del modelo sobre el camino de deformación de la spec."""
    m = FixedSmearedCrack2D(E, NU, FT, gf_model, hypothesis=hyp, softening=law)
    st = m.initial_point_state(h)
    theta = 0.3
    wt = w_tope(h, law)
    path = [strain_of(theta, [s, 0, 0], h=h, hyp=hyp) for s in np.linspace(0, FT, 51)]
    path += [strain_of(theta, [env(w, h, law), 0, 0], w_n=w, h=h, hyp=hyp)
             for w in np.linspace(0.0, 1.2 * wt, n + 1)[1:]]
    W, s_prev, e_prev = 0.0, np.zeros(3), np.zeros(3)
    for eps in path:
        sig, _, st = m.compute_state(eps, st)
        W += 0.5 * (sig + s_prev) @ (eps - e_prev)
        s_prev, e_prev = sig, eps
    return W - 0.5 * s_prev @ e_prev


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("h", [0.005, 0.025, 0.100])
@pytest.mark.parametrize("law", LAWS)
def test_energia_en_traccion_es_G_f_sobre_h_hasta_el_tope(law, h, hyp):
    """tension_energy: D = D_tope de §6 a 10⁻⁶, también en deformación plana
    (sin el factor 1 − ν² del modelo 1)."""
    D = _energia_traccion(law, h, hyp)
    assert D == pytest.approx(d_tope(h, law), rel=1e-6)
    assert D * h / GF == pytest.approx(1.0, abs=3e-3)


def test_energia_en_traccion_control_negativo():
    """Con otra G_f en el material (0.1 % más) la energía se aparta de D_tope
    más que la tolerancia: la prueba anterior discrimina (Reglas §6)."""
    D = _energia_traccion("exponential", H, PS, gf_model=1.001 * GF, n=5000)
    assert abs(D / d_tope(H) - 1.0) > 1e-4


@pytest.mark.parametrize("h", [0.005, 0.025, 0.100])
@pytest.mark.parametrize("lam", [0.5, 1.0, 2.0])
@pytest.mark.parametrize("law", LAWS)
def test_energia_en_modo_mixto_es_1_mas_lambda2_medios_G_f(law, lam, h):
    """mixed_mode_energy (§7, H7): con deslizamiento proporcional a la apertura,
    ``s = λ·(w − w₁)`` desde el primer paso agrietado ``w₁``, la grieta disipa
    por unidad de área ``(1 + λ²/2)·G_f`` (con el descuento del tope y de
    ``w₁``, calculados aquí por cuadratura), igual con las tres ``h``."""
    m = mat(softening=law)
    st = m.initial_point_state(h)
    theta = -0.4
    wt = w_tope(h, law)
    # La rigidez a cortante f/(2w) cambia en la escala de w₁ al abrir: pasos
    # geométricos al principio para que la cuadratura del trabajo no pese.
    w1 = 1e-9 * wt
    ws = np.r_[np.geomspace(w1, 1e-2 * wt, 6000), np.linspace(1e-2 * wt, 1.2 * wt, 14001)[1:]]
    path = [strain_of(theta, [s, 0, 0], h=h) for s in np.linspace(0, FT, 51)]
    for w in ws:
        slip = lam * (w - w1)
        tau = slip / (2.0 * compliance(w, h, law))
        path.append(strain_of(theta, [env(w, h, law), 0.0, tau], w_n=w, slip=slip, h=h))
    W, s_prev, e_prev = 0.0, np.zeros(3), np.zeros(3)
    for eps in path:
        sig, _, st = m.compute_state(eps, st)
        W += 0.5 * (sig + s_prev) @ (eps - e_prev)
        s_prev, e_prev = sig, eps
    D = W - 0.5 * s_prev @ e_prev
    assert st["theta"] == pytest.approx(theta, abs=1e-12)
    ft_ = f_law(wt, law)
    I_n = integrate.quad(lambda w: f_law(w, law), 0.0, wt, epsabs=0, epsrel=1e-13)[0] - 0.5 * ft_ * wt
    tramos = np.r_[np.geomspace(w1, 1e-2 * wt, 8), wt]          # (w − w₁)/w varía en w₁
    I_s = (sum(integrate.quad(lambda w: f_law(w, law) * (w - w1) / (2 * w), a, b, epsabs=0,
                              epsrel=1e-12)[0] for a, b in zip(tramos[:-1], tramos[1:]))
           - ft_ * (wt - w1) ** 2 / (4 * wt))
    assert D * h == pytest.approx(I_n + lam**2 * I_s, rel=1e-5)
    assert D * h / GF == pytest.approx(1.0 + lam**2 / 2.0, rel=3e-3)


def _frozen_energy(mf, st, eps, h):
    """``Ψ(ε; c)`` con la historia de ``st`` congelada (IMPL-EX fuera de un
    paso con los daños de esa historia)."""
    fz = cracked_state(mf, st["theta"], st["w_max_n"], st["w_max_t"], h) if st["cracked"] \
        else mf.initial_point_state(h)
    return 0.5 * mf.compute_state(eps, fz)[0] @ eps


@pytest.mark.parametrize("hyp", [PS, PE])
def test_disipacion_no_negativa_en_cada_paso(hyp):
    """dissipation_non_negative: en caminos al azar con apertura, deslizamiento,
    descarga, cierre y reapertura en las dos direcciones, la energía que libera
    cada paso a deformación fija, ``Ψ(ε; c_n) − Ψ(ε; c_{n+1})``, no es negativa
    (las flexibilidades sólo crecen, también la de cortante al cerrar)."""
    m = mat(hyp)
    mf = mat(hyp, integration="impl-ex")
    mf.set_pseudo_time_ratio(None)
    rng = np.random.default_rng(11)
    theta = 0.4
    N = N_spec(theta)
    for camino in range(4):
        st = m.initial_point_state(H)
        e_loc = np.array([1.2 * K0 / (1 - NU**2), 0.0, 0.0])
        _, _, st = m.compute_state(N @ e_loc, st)
        assert st["cracked"] == 1.0
        abiertos = cerrados = 0
        for _ in range(250):
            e_loc = e_loc + rng.standard_normal(3) * np.array([1.5, 1.5, 2.0]) * K0
            e_loc = np.clip(e_loc, [-6 * K0, -6 * K0, -8 * K0], [12 * K0, 12 * K0, 8 * K0])
            eps = N @ e_loc
            _, _, st1 = m.compute_state(eps, st)
            psi_old = _frozen_energy(mf, st, eps, H)
            psi_new = _frozen_energy(mf, st1, eps, H)
            assert psi_old - psi_new >= -1e-12 * max(psi_old, 1e-30)
            assert st1["w_max_n"] >= st["w_max_n"] and st1["w_max_t"] >= st["w_max_t"]
            abiertos += st1["w_n"] > 0 and st1["w_t"] > 0
            cerrados += st1["w_n"] == 0 and st1["w_t"] == 0
            st = st1
        assert abiertos > 5 and cerrados > 5
        assert st["w_max_t"] > 0.0


# ---------------------------------------------------------------------------
# Cortante (§3) y dirección fija (§1, §10)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("law", LAWS)
@pytest.mark.parametrize("theta", [0.0, 0.4, -0.7])
@pytest.mark.parametrize("descarga", [1.0, 0.5], ids=["en_w_max", "en_descarga"])
def test_ley_de_retencion_de_cortante_y_equivalencia_con_DIANA(law, theta, descarga):
    """shear_retention_law: con la apertura fija, ``τ = f(w_max)/(2·w_max)·s``
    y ``τ/s = 1/(h·(1/G_red − 1/G))`` con ``G_red = E_red/(2(1 + ν·E_red/E))``
    y ``E_red`` el secante de la banda (§3)."""
    m = mat(softening=law)
    wmax = 0.3 * (w_crit() if law == "hordijk" else GF / FT)
    st = m.initial_point_state(H)
    _, _, st = m.compute_state(strain_of(theta, [f_law(wmax, law), 0, 0], w_n=wmax), st)
    assert st["w_max_n"] == pytest.approx(wmax, rel=1e-12)
    w = descarga * wmax
    s_nn = f_law(wmax, law) * descarga
    k_s = f_law(wmax, law) / (2.0 * wmax)
    E_red = s_nn / (s_nn / E + w / H)
    G_red = E_red / (2.0 * (1.0 + NU * E_red / E))
    assert k_s == pytest.approx(1.0 / (H * (1.0 / G_red - 1.0 / G)), rel=1e-12)
    for slip in np.linspace(-2.0 * wmax, 2.0 * wmax, 9):
        tau = k_s * slip
        sig, _, st1 = m.compute_state(strain_of(theta, [s_nn, 0.0, tau], w_n=w, slip=slip), st)
        loc = local(st1["theta"], sig)
        assert loc[2] == pytest.approx(tau, rel=1e-12, abs=1e-12 * FT)
        assert loc[0] == pytest.approx(s_nn, rel=1e-12)
        assert st1["slip"] == pytest.approx(slip, rel=1e-12, abs=1e-12 * wmax)
        assert st1["w_n"] == pytest.approx(w, rel=1e-12)
        assert st1["w_max_n"] == pytest.approx(wmax, rel=1e-12)


def test_cortante_objetivo_la_ley_tau_s_no_depende_de_h():
    """shear_objectivity: la misma ``τ–s`` con h = 5, 25 y 100 mm."""
    wmax = 1.2e-5
    curvas = []
    for h in (0.005, 0.025, 0.100):
        m = mat()
        st = m.initial_point_state(h)
        _, _, st = m.compute_state(strain_of(0.2, [f_law(wmax), 0, 0], w_n=wmax, h=h), st)
        curva = []
        for slip in np.linspace(0.0, 3.0 * wmax, 7):
            # Deformación angular total para ese deslizamiento con la τ del modelo
            # en otra h: se busca la que da ese deslizamiento.
            k_s = f_law(wmax) / (2.0 * wmax)
            eps = strain_of(0.2, [f_law(wmax), 0.0, k_s * slip], w_n=wmax, slip=slip, h=h)
            sig, _, st1 = m.compute_state(eps, st)
            curva.append((st1["slip"], local(0.2, sig)[2]))
        curvas.append(np.array(curva))
    for c in curvas[1:]:
        np.testing.assert_allclose(c, curvas[0], rtol=1e-12, atol=1e-12 * FT * 1e-6)
    np.testing.assert_allclose(curvas[0][:, 1], f_law(wmax) / (2 * wmax) * curvas[0][:, 0],
                               rtol=1e-12, atol=1e-6)


def _referencia_una_grieta(eps, theta, wmax, hyp=PS, h=H):
    """Solución independiente con la grieta ``t`` cerrada: aperturas por
    ``brentq`` en la ley normal y cortante en forma cerrada (§8)."""
    C = Ce(hyp)
    C11, C12 = C[0, 0], C[0, 1]
    e = np.linalg.solve(N_spec(theta), eps)
    b = C11 * e[0] + C12 * e[1]
    c_old = compliance(wmax, h)
    if b <= 0.0:
        w = 0.0
    elif b / (1 / c_old + C11 / h) <= wmax:
        w = b / (1 / c_old + C11 / h)
    else:
        w = optimize.brentq(lambda x: env(x, h) + C11 * x / h - b, wmax, b * h / C11,
                            xtol=1e-300, rtol=4 * np.finfo(float).eps)
    s_nn = C11 * (e[0] - w / h) + C12 * e[1]
    s_tt = C12 * (e[0] - w / h) + C11 * e[1]
    c_s = 2.0 * compliance(max(w, wmax), h)
    tau = e[2] / (1.0 / G + c_s / h)
    return stress_of(theta, [s_nn, s_tt, tau]), w, s_tt


@pytest.mark.parametrize("theta0", [0.35, -0.35], ids=["+", "-"])
@pytest.mark.parametrize("giro", [30.0, 60.0])
def test_direccion_fija_y_esfuerzo_inclinado(theta0, giro):
    """fixed_direction_and_inclined_stress: con la grieta abierta, una carga con
    los ejes principales girados no cambia θ; el tensor completo es el de la
    forma cerrada con ``N`` (§8). Se registra el caveat de §10: el esfuerzo
    principal puede superar f_t en una dirección inclinada sin grieta nueva."""
    m = mat()
    wmax = 0.02 * GF / FT          # recién pasado el pico: la grieta transmite cortante
    st = m.initial_point_state(H)
    eps0 = strain_of(theta0, [f_law(wmax), 0, 0], w_n=wmax)
    _, _, st = m.compute_state(eps0, st)
    inclinado = 0.0
    revisados = 0
    for amp in np.linspace(0.1, 6.0, 60):
        d_eps = strain_of(theta0 + math.radians(giro), [amp * FT, 0.0, 0.0], hyp=PS)
        eps = eps0 + d_eps
        ref, w, s_tt = _referencia_una_grieta(eps, theta0, wmax)
        if s_tt > FT:
            break                   # se formaría la segunda grieta (otra prueba)
        sig, _, st1 = m.compute_state(eps, st)
        assert st1["theta"] == st["theta"] and st1["w_t"] == 0.0
        np.testing.assert_allclose(sig, ref, rtol=1e-10, atol=1e-10 * FT)
        revisados += 1
        assert st1["w_n"] == pytest.approx(w, rel=1e-10, abs=1e-24)
        lam, vec = np.linalg.eigh(np.array([[sig[0], sig[2]], [sig[2], sig[1]]]))
        ang = abs(math.degrees(math.atan2(vec[1, 1], vec[0, 1])) - math.degrees(theta0)) % 180.0
        ang = min(ang, 180.0 - ang)
        if 10.0 < ang < 80.0:
            inclinado = max(inclinado, lam[1] / FT)
    # Caveat de §10 (registro): esfuerzo principal inclinado por encima de f_t
    # (medido el 2026-11-24: 1.36·f_t con el giro de 60° y 1.006·f_t con 30°).
    assert revisados >= 10
    assert inclinado > 1.0


@pytest.mark.parametrize("law", LAWS)
def test_segunda_grieta_perpendicular_y_cortante_en_serie(law):
    """second_crack_orthogonal: la grieta ``t`` se forma a θ₀ + 90° cuando
    ``σ_tt = f_t``; cada dirección sigue su ley; las dos disipan 2·D_tope; el
    cortante con las dos grietas usa ``c_s = 2·(c_n + c_t)``."""
    m = mat(softening=law)
    theta = 0.25
    st = m.initial_point_state(H)
    wt = w_tope(H, law)
    fase1 = [strain_of(theta, [s, 0, 0]) for s in np.linspace(0, FT, 21)]
    fase1 += [strain_of(theta, [env(w, H, law), 0, 0], w_n=w)
              for w in np.linspace(0, 1.2 * wt, 20001)[1:]]
    w_end = 1.2 * wt
    descarga = [strain_of(theta, [env(w_end, H, law) * x, 0, 0], w_n=w_end * x)
                for x in np.linspace(1.0, 0.0, 51)[1:]]
    fase2 = [strain_of(theta, [0, s, 0]) for s in np.linspace(0, FT, 21)]
    fase2 += [strain_of(theta, [0, env(w, H, law), 0], w_t=w)
              for w in np.linspace(0, 1.2 * wt, 20001)[1:]]
    W, s_prev, e_prev = 0.0, np.zeros(3), np.zeros(3)
    primera_t = None
    for k, eps in enumerate(fase1 + descarga + fase2):
        sig, _, st1 = m.compute_state(eps, st)
        W += 0.5 * (sig + s_prev) @ (eps - e_prev)
        s_prev, e_prev = sig, eps
        if st1["w_max_t"] > 0.0 and st["w_max_t"] == 0.0:
            primera_t = (local(theta, sig), st1)
        st = st1
        if st["cracked"]:
            assert st["theta"] == pytest.approx(theta, abs=1e-12)
    loc, st_t = primera_t
    assert loc[1] < FT and loc[1] == pytest.approx(env(st_t["w_t"], H, law), rel=1e-10)
    assert st_t["w_t"] < 1e-3 * wt            # se forma al llegar a f_t
    D = W - 0.5 * s_prev @ e_prev
    assert D == pytest.approx(2.0 * d_tope(H, law), rel=1e-6)
    # Cortante en serie con dos grietas parcialmente abiertas y cerradas.
    a, b = 0.2 * wt, 0.05 * wt
    st = cracked_state(m, theta, a, b, law=law)
    for gam in (3e-5, -7e-5):
        e_loc = np.array([-2 * K0, -1 * K0, gam])
        sig, _, st1 = m.compute_state(N_spec(theta) @ e_loc, st)
        c_s = 2.0 * (compliance(a, H, law) + compliance(b, H, law))
        assert local(theta, sig)[2] == pytest.approx(gam / (1 / G + c_s / H), rel=1e-12)
        assert st1["w_n"] == 0.0 and st1["w_t"] == 0.0


@pytest.mark.parametrize("law", LAWS)
def test_tope_rigidez_residual_sin_division_entre_cero(law):
    """cap_reached: tras el tope ``σ_nn = h·e_nn/c_max``, ``c_s ≤ 4·c_max``; con
    Hordijk y ``w ≥ w_c`` (f = 0) no hay división entre cero."""
    m = mat(softening=law)
    cm = c_max(H)
    for w in (1.5 * w_tope(H, law), 3.0 * w_crit(), 20.0 * w_crit()):
        st = m.initial_point_state(H)
        sig, C, st = m.compute_state(strain_of(0.1, [w / cm, 0, 0], w_n=w), st)
        assert np.all(np.isfinite(sig)) and np.all(np.isfinite(C))
        # σ_nn = C11·(e_nn − w/h) se resta de deformaciones ~w/h: redondeo absoluto.
        assert local(0.1, sig)[0] == pytest.approx(w / cm, rel=1e-12, abs=1e-12 * FT)
        assert st["w_n"] == pytest.approx(w, rel=1e-12)
        assert st["damage_n"] == pytest.approx(CAP, rel=1e-14)
    st = cracked_state(m, 0.1, 30 * w_crit(), 30 * w_crit(), law=law)
    gam = 4e-5
    # Cerradas con una compresión pequeña: τ (≈ 0.3 Pa) no se mezcla en la
    # rotación con esfuerzos normales de MPa.
    sig, C, _ = m.compute_state(N_spec(0.1) @ np.array([-1e-6 * K0, -1e-6 * K0, gam]), st)
    assert local(0.1, sig)[2] == pytest.approx(gam / (1 / G + 4 * cm / H), rel=1e-12)
    assert np.all(np.isfinite(C))


# ---------------------------------------------------------------------------
# Problema local (§8): robustez, idempotencia, descarga y cierre
# ---------------------------------------------------------------------------

def _comprobar_solucion_local(m, st_old, st, eps, sig, h, law, hyp):
    """Comprueba sobre la solución la ley de cada dirección (§2, §6) y el
    cortante (§8), con la historia confirmada ``st_old``."""
    theta = st["theta"]
    e = np.linalg.solve(N_spec(theta), eps)
    loc = local(theta, sig)
    C = Ce(hyp)
    C11, C12 = C[0, 0], C[0, 1]
    w = (st["w_n"], st["w_t"])
    scale = max(FT, np.abs(loc).max())
    tol = 1e-10 * scale
    assert abs(loc[0] - (C11 * (e[0] - w[0] / h) + C12 * (e[1] - w[1] / h))) <= tol
    assert abs(loc[1] - (C12 * (e[0] - w[0] / h) + C11 * (e[1] - w[1] / h))) <= tol
    agrietado = st_old["cracked"] > 0
    for i, clave in enumerate(("w_max_n", "w_max_t")):
        wm = st_old[clave] if agrietado else 0.0
        s = loc[i]
        if w[i] == 0.0:
            assert s <= (FT * (1 + 2e-10) if wm == 0.0 else 1e-12 * scale)
        elif w[i] <= wm:
            assert abs(s - w[i] / compliance(wm, h, law)) <= tol
        else:
            assert abs(s - env(w[i], h, law)) <= tol
        assert st[clave] == max(wm, w[i])
    c_s = 2.0 * (compliance(st["w_max_n"], h, law) + compliance(st["w_max_t"], h, law))
    assert abs(loc[2] - e[2] / (1 / G + c_s / h)) <= tol


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("law", LAWS)
def test_problema_local_robusto_cerca_de_los_quiebres(law, hyp):
    """local_problem_robustness: estados al azar junto a ``w = 0`` y
    ``w = w_max``, una y dos direcciones, ``h`` hasta 0.95 del límite: la
    solución cumple la ley de cada dirección y se repite idéntica."""
    m = mat(hyp, softening=law)
    rng = np.random.default_rng(19)
    lim = m.characteristic_length_limit
    escala = w_crit() if law == "hordijk" else GF / FT
    casos = {"dos_abiertas": 0, "una": 0, "cerrada": 0}
    for frac in (0.05, 0.5, 0.95):
        h = frac * lim
        for _ in range(120):
            theta = rng.uniform(-math.pi / 2, math.pi / 2)
            wm = [rng.choice([0.0, rng.uniform(0.01, 0.5) * escala]) for _ in range(2)]
            if wm[0] == 0.0 and wm[1] > 0.0:
                wm = [wm[1], 0.0]
            st_old = cracked_state(m, theta, wm[0], wm[1], h=h, law=law) if wm[0] > 0 \
                else m.initial_point_state(h)
            s_loc, w_loc = [], []
            for i in range(2):
                opcion = rng.integers(4)
                if opcion == 0:                                  # junto a w = 0
                    s_loc.append((FT if wm[i] == 0 else 0.0) * (1 + rng.uniform(-1e-3, 1e-3))
                                 - rng.uniform(0, 1e-3) * FT * (wm[i] > 0))
                    w_loc.append(0.0)
                elif opcion == 1 and wm[i] > 0:                  # junto a w_max
                    x = wm[i] * (1 + rng.uniform(-1e-6, 1e-6))
                    s_loc.append(env(x, h, law) if x > wm[i] else x / compliance(wm[i], h, law))
                    w_loc.append(x)
                elif opcion == 2:                                # en carga
                    x = max(wm[i], 1e-3 * escala) * rng.uniform(1.0, 3.0)
                    s_loc.append(env(x, h, law))
                    w_loc.append(x)
                else:                                            # comprimida
                    s_loc.append(-rng.uniform(0.1, 5.0) * FT)
                    w_loc.append(0.0)
            if wm[0] == 0.0:
                # Sin grieta: θ debe salir del esfuerzo de prueba; se fija la
                # dirección principal mayor en theta.
                s_loc = sorted(s_loc, reverse=True)
                w_loc = sorted(w_loc, reverse=True)
            eps = strain_of(theta, [s_loc[0], s_loc[1], rng.uniform(-1, 1) * FT],
                            w_n=w_loc[0], w_t=w_loc[1], h=h, hyp=hyp)
            eps = eps + rng.standard_normal(3) * 1e-4 * K0
            sig, C, st = m.compute_state(eps, st_old)
            sig2, C2, st2 = m.compute_state(eps, st_old)
            np.testing.assert_array_equal(sig, sig2)
            np.testing.assert_array_equal(C, C2)
            assert st == st2
            if st["cracked"]:
                _comprobar_solucion_local(m, st_old, st, eps, sig, h, law, hyp)
            casos["dos_abiertas" if st["w_n"] > 0 and st["w_t"] > 0 else
                  "una" if st["w_n"] > 0 or st["w_t"] > 0 else "cerrada"] += 1
    assert min(casos.values()) > 20, casos


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("law", LAWS)
def test_problema_local_misma_solucion_por_las_dos_rutas(law, hyp):
    """local_problem_robustness, «la misma solución desde arranques distintos»
    (verificación del 2026-11-24, O3): el algoritmo anidado resuelve primero
    n y después t; con θ + 90° y las historias permutadas resuelve el mismo
    problema en el otro orden. El esfuerzo global es el mismo, las aperturas
    se intercambian y el deslizamiento cambia de signo (t' = −n)."""
    m = mat(hyp, softening=law)
    rng = np.random.default_rng(31)
    escala = w_crit() if law == "hordijk" else GF / FT
    dos = 0
    for frac in (0.05, 0.5, 0.95):
        h = frac * m.characteristic_length_limit
        for _ in range(150):
            theta = rng.uniform(-1.5, 1.5)
            a, b = (rng.choice([0.0, rng.uniform(0.005, 0.5) * escala]) for _ in range(2))
            if a == 0.0 and b == 0.0:
                a = 0.1 * escala
            st1 = cracked_state(m, theta, a, b, h=h, law=law)
            st2 = cracked_state(m, theta + 0.5 * math.pi, b, a, h=h, law=law)
            e_loc = rng.uniform([-2, -2, -4], [1.5, 1.5, 4]) * escala / h
            eps = N_spec(theta) @ e_loc
            s1, _, r1 = m.compute_state(eps, st1)
            s2, _, r2 = m.compute_state(eps, st2)
            np.testing.assert_allclose(s2, s1, rtol=1e-10, atol=1e-10 * FT)
            assert r2["w_n"] == pytest.approx(r1["w_t"], rel=1e-9, abs=1e-12 * escala)
            assert r2["w_t"] == pytest.approx(r1["w_n"], rel=1e-9, abs=1e-12 * escala)
            assert r2["slip"] == pytest.approx(-r1["slip"], rel=1e-9, abs=1e-12 * escala)
            dos += r1["w_n"] > 0 and r1["w_t"] > 0
    assert dos > 20


@pytest.mark.parametrize("law", LAWS)
def test_doble_inicio_casi_equibiaxial_determinista(law):
    """Doble inicio con el esfuerzo de prueba casi equibiaxial y ``h`` = 0.95
    del límite: una sola solución, la misma en cada llamada, con las dos
    grietas abiertas y la mayor en ``n``."""
    m = mat(softening=law)
    h = 0.95 * m.characteristic_length_limit
    C = Ce(PS)
    for d in (0.0, 1e-9, 1e-6, 1e-3):
        s = 1.0498 * FT
        eps = np.linalg.solve(C, [s * (1 + d), s, 0.0])
        out = [m.compute_state(eps, m.initial_point_state(h)) for _ in range(2)]
        np.testing.assert_array_equal(out[0][0], out[1][0])
        assert out[0][2] == out[1][2]
        st = out[0][2]
        assert st["theta"] == 0.0
        assert st["w_n"] >= st["w_t"] * (1 - 1e-12) and st["w_t"] > 0.0
        _comprobar_solucion_local(m, m.initial_point_state(h), st, eps, out[0][0], h, law, PS)


@pytest.mark.parametrize("law", LAWS)
def test_estado_confirmado_vuelto_a_evaluar_no_carga(law):
    """Al empezar un paso, el punto se evalúa en su estado confirmado: la
    historia no crece y la tangente es la de descarga (secante, simétrica sin
    cortante), sin depender del redondeo en el quiebre w = w_max (decisión con
    la tolerancia de admisibilidad, ADR 0006). Sin la tolerancia, una parte de
    los puntos caía en la carga con tangente de ablandamiento (regresión: la
    barra de Quad8 se iba a la rama localizada dentro del elemento)."""
    m = mat(softening=law)
    rng = np.random.default_rng(29)
    escala = w_crit() if law == "hordijk" else GF / FT
    cargando = 0
    for _ in range(400):
        theta = rng.uniform(-1.5, 1.5)
        st0 = cracked_state(m, theta, rng.uniform(0.01, 0.4) * escala,
                            rng.choice([0.0, rng.uniform(0.01, 0.2) * escala]), law=law)
        w = st0["w_max_n"] * rng.uniform(1.0, 1.5)
        wt = st0["w_max_t"] * rng.uniform(1.0, 1.5) if st0["w_max_t"] > 0 else 0.0
        s_t = env(wt, H, law) if wt > 0 else -0.3 * FT
        eps = strain_of(theta, [env(w, H, law), s_t, 0.0], w_n=w, w_t=wt)
        _, _, st1 = m.compute_state(eps, st0)
        cargando += st1["w_max_n"] > st0["w_max_n"]
        _, C2, st2 = m.compute_state(eps, st1)
        assert st2["w_max_n"] == st1["w_max_n"] and st2["w_max_t"] == st1["w_max_t"]
        loc = N_spec(theta).T @ C2 @ N_spec(theta)
        assert loc[0, 0] > 0.0                       # secante, no ablandamiento
        np.testing.assert_allclose(C2, C2.T, rtol=0, atol=1e-12 * np.abs(C2).max())
    assert cargando >= 300


def test_paso_implicito_idempotente_y_el_intento_rechazado_no_fija_theta():
    """implicit_step_idempotent: el mismo intento dos veces da el mismo estado
    de prueba; un intento rechazado no fija θ (sólo el paso confirmado)."""
    m = mat()
    st = m.initial_point_state(H)
    rechazado = strain_of(0.7, [1.5 * FT, 0, 0], w_n=1e-6)
    _, _, st_r = m.compute_state(rechazado, st)
    assert st_r["theta"] == pytest.approx(0.7, abs=1e-12)
    bueno = strain_of(-0.2, [1.2 * FT, 0.1 * FT, 0], w_n=0.0)
    a = m.compute_state(bueno, st)
    b = m.compute_state(bueno, st)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])
    assert a[2] == b[2]
    assert a[2]["theta"] == pytest.approx(-0.2, abs=1e-12) and a[2]["cracked"] == 1.0
    assert st["cracked"] == 0.0 and st["theta"] == 0.0


@pytest.mark.parametrize("hyp", [PS, PE])
def test_descarga_secante_cierre_y_reapertura(hyp):
    """unloading_closure_reopening: descarga ``σ = f(w_max)·w/w_max``; cerrada,
    rigidez ``E`` y Poisson ``ν`` en compresión; reapertura por la secante y
    después la envolvente; ``w_max`` nunca decrece."""
    m = mat(hyp)
    theta = 0.3
    wmax = 0.25 * GF / FT
    st = m.initial_point_state(H)
    for w in np.linspace(0, wmax, 11)[1:]:
        _, _, st = m.compute_state(strain_of(theta, [f_law(w), 0, 0], w_n=w, hyp=hyp), st)
    assert st["w_max_n"] == pytest.approx(wmax, rel=1e-12)
    pasos = ([(f_law(wmax) * x, wmax * x) for x in (0.8, 0.5, 0.2)]
             + [(-s, 0.0) for s in (0.5 * FT, 3.0 * FT, 10.0 * FT)]
             + [(f_law(wmax) * x, wmax * x) for x in (0.3, 0.9)]
             + [(f_law(wmax * x), wmax * x) for x in (1.1, 1.3)])
    wm = wmax
    for s, w in pasos:
        eps = strain_of(theta, [s, 0, 0], w_n=w, hyp=hyp)
        sig, _, st1 = m.compute_state(eps, st)
        np.testing.assert_allclose(sig, stress_of(theta, [s, 0, 0]), rtol=0, atol=1e-12 * 10 * FT)
        assert st1["w_n"] == pytest.approx(w, rel=1e-12, abs=1e-30)
        assert st1["w_max_n"] >= st["w_max_n"]
        if s < 0:
            # Cerrada: compresión elástica con E y ν (deformación del concreto).
            e = np.linalg.solve(N_spec(theta), eps)
            assert e[0] == pytest.approx(s / E * ((1 - NU**2) if hyp == PE else 1.0), rel=1e-12)
        wm = max(wm, w)
        assert st1["w_max_n"] == pytest.approx(wm, rel=1e-12)
        st = st1


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("fw", [0.0, 0.5, 1.0, 2.0])
def test_compresion_paralela_a_grieta_abierta_se_conserva(hyp, fw):
    """parallel_compression_open_crack: con la grieta abierta en ``n`` y
    compresión ``ε_t`` paralela, ``σ_tt − ν·σ_nn = E·ε_t`` (esfuerzo plano) o
    ``σ_tt = E·ε_t/(1 − ν²) + ν·σ_nn/(1 − ν)`` (deformación plana), con
    cualquier apertura; tiende a la elástica completa al abrirse del todo."""
    m = mat(hyp, softening="hordijk")
    theta = -0.6
    w = max(fw * w_crit(), 1e-6 * w_crit())        # «0»: recién abierta
    st = cracked_state(m, theta, w, law="hordijk")
    C = Ce(hyp)
    e_t = -1e-3
    s_nn = env(w, H, "hordijk")
    e_n = w / H + (s_nn - C[0, 1] * e_t) / C[0, 0]
    sig, _, st1 = m.compute_state(N_spec(theta) @ np.array([e_n, e_t, 0.0]), st)
    loc = local(theta, sig)
    assert st1["w_n"] == pytest.approx(w, rel=1e-10)
    assert loc[0] == pytest.approx(s_nn, rel=1e-10, abs=1e-10 * FT)
    if hyp == PS:
        assert loc[1] - NU * loc[0] == pytest.approx(E * e_t, rel=1e-10)
    else:
        assert loc[1] == pytest.approx(E * e_t / (1 - NU**2) + NU * loc[0] / (1 - NU), rel=1e-10)


def _ciclos(seed, n_ciclos):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_ciclos):
        v = [np.array([rng.uniform(-3, 5), rng.uniform(-3, 5), rng.uniform(-4, 4)]) * K0
             for _ in range(rng.integers(3, 5))]
        out.append(v + [v[0]])
    return out


def _trabajo_del_ciclo(m, st, theta, vertices, n_lado, recuperar=False):
    N = N_spec(theta)
    W, E_max = 0.0, 0.0
    s_prev = e_prev = None
    for a, b in zip(vertices[:-1], vertices[1:]):
        for x in np.linspace(0.0, 1.0, n_lado + 1)[(0 if s_prev is None else 1):]:
            e_loc = a + x * (b - a)
            eps = N @ e_loc
            sig = m.compute_state(eps, st)[0]
            if recuperar:
                loc = local(theta, sig)
                if loc[0] <= 0.0 and loc[1] <= 0.0:              # cerrada: τ = G·γ
                    sig = stress_of(theta, [loc[0], loc[1], G * e_loc[2]])
            if s_prev is not None:
                W += 0.5 * (sig + s_prev) @ (eps - e_prev)
            E_max = max(E_max, 0.5 * sig @ eps)
            s_prev, e_prev = sig, eps
    return W, E_max


@pytest.mark.parametrize("hyp", [PS, PE])
def test_trabajo_nulo_en_ciclos_con_historia_fija(hyp):
    """zero_work_fixed_history (§6): con ν = 0.2 y la historia fija, la ley
    deriva del potencial χ y ∮σ:dε = 0 frente a la energía máxima del ciclo,
    con grietas abiertas y cerradas. Control negativo: recuperar la rigidez a
    cortante al cerrar la grieta rompe el potencial (≥ 0.1)."""
    m = mat(hyp, integration="impl-ex")
    m.set_pseudo_time_ratio(None)
    theta = 0.45
    ctrl = 0.0
    # Ciclo diseñado: desliza abierta, se cierra con el deslizamiento y vuelve.
    disenado = [np.array(v) * K0 for v in ([4, -1, 0], [4, -1, 4], [-2, -1, 4], [-2, -1, 0],
                                            [4, -1, 0])]
    for k, (wn, wt) in enumerate(((0.05 * GF / FT, 0.0), (0.08 * GF / FT, 0.02 * GF / FT))):
        st = cracked_state(m, theta, wn, wt)
        for v in _ciclos(23 + k, 6) + [disenado]:
            W, E_max = _trabajo_del_ciclo(m, st, theta, v, 3000)
            assert abs(W) <= 1e-6 * E_max
            Wc, Ec = _trabajo_del_ciclo(m, st, theta, v, 300, recuperar=True)
            ctrl = max(ctrl, abs(Wc) / Ec)
    assert ctrl >= 0.1


# ---------------------------------------------------------------------------
# Tangente consistente (§8)
# ---------------------------------------------------------------------------

def _fd(m, eps, st, rel=1e-7):
    C_fd = np.zeros((3, 3))
    d = rel * np.abs(eps).max()
    for j in range(3):
        ep, em = eps.copy(), eps.copy()
        ep[j] += d
        em[j] -= d
        C_fd[:, j] = (m.compute_state(ep, st)[0] - m.compute_state(em, st)[0]) / (2 * d)
    return C_fd


def _caso_tangente(caso, law, hyp):
    """(historia, deformación) de cada caso, lejos de los cambios de rama."""
    a = 0.15 * (w_crit() if law == "hordijk" else GF / FT)
    theta = 0.37
    if caso == "carga_normal":
        wm, w, s_tt, slip = (a, 0.0), (1.3 * a, 0.0), -0.5 * FT, 0.0
    elif caso == "carga_con_t_en_traccion":
        # Segunda dirección sin grieta, en tracción por debajo de f_t: cerrada
        # en la tangente (la envolvente empieza en f_t).
        wm, w, s_tt, slip = (a, 0.0), (1.3 * a, 0.0), 0.6 * FT, 0.4 * a
    elif caso == "carga_con_cortante":
        wm, w, s_tt, slip = (a, 0.0), (1.3 * a, 0.0), -0.5 * FT, 0.8 * a
    elif caso == "dos_grietas":
        wm, w, s_tt, slip = (a, 0.5 * a), (1.3 * a, 0.8 * a), None, 0.6 * a
    elif caso == "descarga":
        wm, w, s_tt, slip = (a, 0.4 * a), (0.6 * a, 0.2 * a), None, 0.5 * a
    else:  # cerrada
        wm, w, s_tt, slip = (a, 0.4 * a), (0.0, 0.0), None, 0.5 * a
    c_s = 2.0 * (compliance(max(wm[0], w[0]), H, law) + compliance(max(wm[1], w[1]), H, law))
    s = []
    for i in range(2):
        if w[i] == 0.0:
            s.append(s_tt if (i == 1 and s_tt is not None) else -0.7 * FT)
        elif w[i] > wm[i]:
            s.append(env(w[i], H, law))
        else:
            s.append(w[i] / compliance(wm[i], H, law))
    return wm, theta, strain_of(theta, [s[0], s[1], slip / c_s], w_n=w[0], w_t=w[1],
                                slip=slip, hyp=hyp)


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("law", LAWS)
@pytest.mark.parametrize("caso", ["carga_normal", "carga_con_cortante", "carga_con_t_en_traccion",
                                  "dos_grietas", "descarga", "cerrada"])
def test_tangente_implicita_igual_a_diferencias_finitas(caso, law, hyp):
    """tangent_fd: |C_FD − C_tan| ≤ 10⁻⁶·max|C_FD|; asimétrica en carga con
    cortante (falla si se olvida ∂τ/∂e) y simétrica en carga normal pura y a
    historia fija."""
    m = mat(hyp, softening=law)
    wm, theta, eps = _caso_tangente(caso, law, hyp)
    st = cracked_state(m, theta, wm[0], wm[1], law=law)
    sig, C, st1 = m.compute_state(eps, st)
    if caso in ("carga_normal", "carga_con_cortante", "carga_con_t_en_traccion", "dos_grietas"):
        assert st1["w_max_n"] > st["w_max_n"]
    if caso == "dos_grietas":
        assert st1["w_max_t"] > st["w_max_t"]
    C_fd = _fd(m, eps, st)
    assert np.abs(C - C_fd).max() <= 1e-6 * np.abs(C_fd).max()
    asim = np.abs(C - C.T).max() / np.abs(C).max()
    if caso in ("carga_con_cortante", "carga_con_t_en_traccion", "dos_grietas"):
        assert asim > 1e-3
    else:
        assert asim <= 1e-12


def test_tangente_en_el_paso_de_inicio_sin_d_theta():
    """En el paso de inicio θ depende de ε y la tangente no incluye ∂θ/∂ε (§8):
    con un esfuerzo de prueba sin cortante en los ejes de la grieta coincide
    con diferencias finitas en las componentes normales a la grieta."""
    m = mat()
    st = m.initial_point_state(H)
    eps = strain_of(0.0, [1.3 * FT, -0.2 * FT, 0.0], w_n=0.0)
    _, C, st1 = m.compute_state(eps, st)
    assert st1["cracked"] == 1.0
    C_fd = _fd(m, eps, st)
    np.testing.assert_allclose(C[:, :2], C_fd[:, :2], rtol=0, atol=1e-6 * np.abs(C_fd).max())


def _malla(integration, hypothesis):
    """Bloque empotrado en x = 0 con un campo de tracción y cortante no
    uniforme (N, mm) que agrieta todos los puntos en el primer paso."""
    from solidum.math.assembly import Assembler
    from solidum.mesh import structured_mesh
    E_, nu_, ft_, gf_ = 31700.0, 0.23, 3.1, 0.095
    k0 = ft_ / E_
    m = FixedSmearedCrack2D(E_, nu_, ft_, gf_, hypothesis=hypothesis, integration=integration)
    mesh = structured_mesh(Quad4, (5, 4), size=(50.0, 40.0), material=m, thickness=7.0)
    for nd in mesh.nodes_on("xmin"):
        nd.fix_dof("ux", 0.0)
        nd.fix_dof("uy", 0.0)
    dom = mesh.domain
    dom.generate_equation_numbers(verbose=False)
    U1 = np.zeros(dom.total_dofs)
    U2 = np.zeros(dom.total_dofs)
    for nd in dom.nodes.values():
        x, y = nd.coordinates[0], nd.coordinates[1]
        U1[nd.dofs["ux"]] = k0 * (3.6 * x * (1 + 0.004 * y) + 1.5 * y)
        U1[nd.dofs["uy"]] = k0 * (0.9 * x - 0.6 * y + 0.006 * x * y)
        # Segundo paso no proporcional: cortante en los ejes de las grietas.
        U2[nd.dofs["ux"]] = U1[nd.dofs["ux"]] + k0 * (2.0 * x + 3.0 * y)
        U2[nd.dofs["uy"]] = U1[nd.dofs["uy"]] + k0 * (2.5 * x + 0.8 * y - 0.01 * x * y)
    return Assembler(dom), U1, U2


def _tangente_global_frente_a_df(asm, U):
    K, _ = asm.assemble_non_linear_system(U)
    K = K.toarray() if hasattr(K, "toarray") else np.asarray(K)
    nK = K.shape[0]
    rng = np.random.default_rng(53)
    peor = 0.0
    for _ in range(3):
        dU = np.zeros_like(U)
        dU[:nK] = rng.standard_normal(nK)
        h = 1e-7 * np.abs(U).max()
        dF = (asm.assemble_non_linear_system(U + h * dU)[1]
              - asm.assemble_non_linear_system(U - h * dU)[1])[:nK] / (2.0 * h)
        peor = max(peor, np.abs(K @ dU[:nK] - dF).max() / np.abs(dF).max())
    return peor, np.abs(K - K.T).max() / np.abs(K).max()


@pytest.mark.parametrize("hyp", [PS, PE])
def test_tangente_global_implicita_ensamblada(hyp):
    """Tras un paso confirmado que agrieta todos los puntos, la tangente
    asimétrica ensamblada es la de las fuerzas internas (K·δU = δF_int)."""
    asm, U1, U2 = _malla("implicit", hyp)
    asm.prepare_all_steps(np.zeros_like(U1), pseudo_time_increment=1.0)
    asm.assemble_non_linear_system(U1)
    asm.commit_all_states()
    agrietados = [s["cracked"] for el in asm.domain.elements.values() for s in el.state_vars]
    assert min(agrietados) == 1.0                 # ningún inicio en el paso evaluado
    asm.prepare_all_steps(U1, pseudo_time_increment=1.0)
    peor, asim = _tangente_global_frente_a_df(asm, U2)
    assert peor <= 1e-6
    assert asim > 1e-3


def test_tangente_global_impl_ex_simetrica():
    asm, U1, U2 = _malla("impl-ex", PS)
    U_prev = np.zeros_like(U1)
    for k in range(1, 5):
        asm.prepare_all_steps(U_prev, pseudo_time_increment=1.0)
        asm.assemble_non_linear_system(k / 4.0 * U1)
        asm.commit_all_states()
        U_prev = k / 4.0 * U1
    asm.prepare_all_steps(U_prev, pseudo_time_increment=1.0)
    peor, asim = _tangente_global_frente_a_df(asm, U2)
    assert peor <= 1e-6
    assert asim <= 1e-12


# ---------------------------------------------------------------------------
# Unidades y kernel por lotes
# ---------------------------------------------------------------------------

def _camino_al_azar(seed, n=80):
    rng = np.random.default_rng(seed)
    base = np.array([1.0, -0.15, 0.35])
    out = []
    for f in np.r_[np.linspace(0, 6, 31)[1:], np.linspace(6, -3, 15)[1:], np.linspace(-3, 9, 25)[1:]]:
        out.append(f * K0 * base + rng.standard_normal(3) * 0.3 * K0 * (f > 0))
    return out[:n]


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
@pytest.mark.parametrize("law", LAWS)
def test_invariancia_de_unidades(hyp, integration, law):
    """units_invariance: el mismo caso en Pa–m y en MPa–mm con todos los
    parámetros ≠ 1: esfuerzos ×10⁻⁶, aperturas ×10³, daños y energías iguales."""
    kw = dict(hypothesis=hyp, integration=integration, softening=law)
    a = FixedSmearedCrack2D(31.7e9, 0.18, 2.9e6, 93.0, **kw)
    b = FixedSmearedCrack2D(31.7e3, 0.18, 2.9, 0.093, **kw)
    for m in (a, b):
        m.set_pseudo_time_ratio(1.0)
    sa, sb = a.initial_point_state(0.037), b.initial_point_state(37.0)
    Wa = Wb = 0.0
    e_prev = np.zeros(3)
    pa = pb = np.zeros(3)
    k0 = 2.9e6 / 31.7e9
    for eps in _camino_al_azar(41):
        eps = eps * k0 / K0
        siga, _, sa = a.compute_state(eps, sa)
        sigb, _, sb = b.compute_state(eps, sb)
        np.testing.assert_allclose(sigb, 1e-6 * siga, rtol=1e-12, atol=1e-12 * np.abs(sigb).max())
        for k in ("damage_n", "damage_t", "theta", "cracked", "damage_error"):
            assert sb[k] == pytest.approx(sa[k], rel=1e-12, abs=1e-14), k
        for k in ("w_max_n", "w_max_t", "w_n", "w_t", "slip"):
            assert sb[k] == pytest.approx(1e3 * sa[k], rel=1e-11, abs=1e-14), k   # mm
        Wa += 0.5 * (siga + pa) @ (eps - e_prev)
        Wb += 0.5 * (sigb + pb) @ (eps - e_prev)
        e_prev, pa, pb = eps, siga, sigb
    assert sa["damage_n"] > 0.5
    assert Wb == pytest.approx(1e-6 * Wa, rel=1e-12)


@pytest.mark.parametrize("hyp", [PS, PE])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_kernel_por_lotes_reproduce_compute_state_bit_a_bit(hyp, integration):
    m = mat(hyp, integration=integration, softening="hordijk" if hyp == PE else "exponential")
    names = list(m.STATE_SCHEMA)
    kernel = m.batch_kernel()
    C = m.batch_matrix()
    st = m.initial_point_state(H)
    for i, eps in enumerate(_camino_al_azar(3)):
        m.set_pseudo_time_ratio(None if i % 5 == 4 else 0.5 + 0.1 * i)
        params = m.batch_params()
        sig_ref, C_ref, st_new = m.compute_state(eps, st)
        S_old = np.array([st[k] for k in names])
        S_new, sig, flag = np.zeros(len(names)), np.zeros(3), np.zeros(1, np.int8)
        C_t = kernel(eps, S_old, S_new, params, C, sig, np.zeros((3, 3)), flag)
        np.testing.assert_array_equal(sig, sig_ref)
        np.testing.assert_array_equal(C_t, C_ref)
        np.testing.assert_array_equal(S_new, [st_new[k] for k in names])
        assert flag[0] == 0
        st = st_new
    assert st["damage_n"] > 0.0


def test_esfuerzo_fuera_del_plano():
    rng = np.random.default_rng(43)
    for hyp in (PS, PE):
        m = mat(hyp)
        names = list(m.STATE_SCHEMA)
        sig, S, ref = [], [], []
        st0 = cracked_state(m, 0.3, 0.1 * GF / FT, 0.0)
        for _ in range(20):
            s, _, st = m.compute_state(rng.standard_normal(3) * 3e-4, st0)
            sig.append(s)
            S.append([st[k] for k in names])
            ref.append(m.out_of_plane_stress(s, st))
            assert ref[-1] == (NU * (s[0] + s[1]) if hyp == PE else 0.0)
        np.testing.assert_allclose(m.batch_out_of_plane_stress(np.array(sig), np.array(S)), ref,
                                   rtol=1e-15, atol=0)


# ---------------------------------------------------------------------------
# IMPL-EX (§9)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("hyp", [PS, PE])
def test_implex_tangente_secante_simetrica_y_definida_positiva(hyp):
    """secant_tangent: dentro de un paso las flexibilidades extrapoladas no
    dependen de ε; la tangente es la secante (= diferencias finitas),
    simétrica y definida positiva, también con el tope."""
    m = mat(hyp, integration="impl-ex")
    m.set_pseudo_time_ratio(1.0)
    rng = np.random.default_rng(47)
    theta = 0.3
    N = N_spec(theta)
    for wn, wt in ((0.1 * GF / FT, 0.0), (0.2 * GF / FT, 0.05 * GF / FT), (1.0, 1.0)):
        st = cracked_state(m, theta, wn, wt)
        st["w_max_n_prev"], st["w_max_t_prev"] = 0.9 * wn, 0.9 * wt
        for _ in range(12):
            e_loc = rng.uniform([-4, -4, -6], [8, 8, 6]) * K0
            if min(abs(local(theta, m.compute_state(N @ e_loc, st)[0])[:2])) < 0.05 * FT:
                continue
            eps = N @ e_loc
            sig, C, _ = m.compute_state(eps, st)
            assert np.abs(C - _fd(m, eps, st)).max() <= 1e-6 * np.abs(C).max()
            np.testing.assert_array_equal(C, C.T)
            assert np.linalg.eigvalsh(C).min() > 0.0


def _camino_implex(n, integration, cortante, forma="lineal"):
    """Camino en los ejes θ = 0.3 que arranca justo en el umbral (la grieta se
    forma en el primer paso con cualquier n, como en los modelos 1 y 2):
    lineal en t, o cuadrático (velocidad nula en el umbral, sin quiebre); la
    segunda dirección no llega a f_t."""
    m = mat(integration=integration)
    m.set_pseudo_time_ratio(1.0)
    N = N_spec(0.3)
    e_c = np.linalg.solve(Ce(PS), [FT, 0.0, 0.0])
    st = m.initial_point_state(H)
    _, _, st = m.compute_state(N @ e_c, st)
    assert st["cracked"] == 0.0
    d = np.array([30.0, -3.0, 10.0 if cortante else 0.0]) * K0
    out = []
    for t in np.linspace(0.0, 1.0, n + 1)[1:]:
        sig, _, st = m.compute_state(N @ (e_c + (t if forma == "lineal" else t * t) * d), st)
        out.append(sig)
    assert st["w_max_t"] == 0.0
    return np.array(out), st


@pytest.mark.parametrize("cortante", [False, True], ids=["sin_cortante", "con_cortante"])
@pytest.mark.parametrize("forma, orden", [("lineal", 1.0), ("cuadratica", 2.0)])
def test_implex_converge_a_la_implicita(forma, orden, cortante):
    """first_order_convergence: el error máximo de σ frente a la implícita se
    reduce a la mitad con el paso (pasos Δ, Δ/2, Δ/4, Δ/8; razón 2 ± 15 %)
    con el quiebre del inicio en el camino lineal, y a la cuarta parte donde
    el camino es suave (cuadrático desde el umbral): el segundo orden es el
    de la extrapolación (sin ella, primer orden), como en los modelos 1 y 2."""
    errores = []
    for n in (100, 200, 400, 800):
        s_ex, st = _camino_implex(n, "impl-ex", cortante, forma)
        s_im, _ = _camino_implex(n, "implicit", cortante, forma)
        errores.append(np.abs(s_ex - s_im).max())
    assert st["damage_n"] > 0.9
    razones = np.array(errores[:-1]) / np.array(errores[1:])
    np.testing.assert_allclose(razones, 2.0 ** orden, rtol=0.15)


def test_implex_historia_monotona_con_extrapolacion_que_se_pasa():
    """history_monotone: carga, descarga parcial y recarga con una razón que
    hace pasarse a la extrapolación: ``ω̃`` (y la flexibilidad) nunca decrece."""
    m = mat(integration="impl-ex")
    st = m.initial_point_state(H)
    base = np.array([1.0, -0.2, 0.3])
    camino = np.r_[np.linspace(0.0, 8.0, 41)[1:], np.linspace(8.0, 3.0, 11)[1:],
                   np.linspace(3.0, 15.0, 31)[1:]] * K0
    prev = (0.0, 0.0)
    for k, e in enumerate(camino):
        m.set_pseudo_time_ratio(2.0 if k % 7 == 0 else 1.0)
        _, _, st = m.compute_state(e * base, st)
        assert st["damage_n"] >= prev[0] and st["damage_t"] >= prev[1]
        assert st["damage_error"] >= 0.0
        prev = (st["damage_n"], st["damage_t"])
    assert prev[0] > 0.9


def test_implex_grieta_nueva_entra_en_el_paso_siguiente():
    """La grieta se detecta en la etapa implícita: el paso en que se forma es
    elástico y la extrapolación la usa desde el siguiente (§9)."""
    m = mat(integration="impl-ex")
    m.set_pseudo_time_ratio(1.0)
    st = m.initial_point_state(H)
    eps = np.array([1.5 * K0, 0.0, 0.0])
    sig, C, st1 = m.compute_state(eps, st)
    np.testing.assert_allclose(sig, Ce(PS) @ eps, rtol=1e-14)
    assert st1["cracked"] == 1.0 and st1["w_max_n"] > 0.0 and st1["damage_n"] == 0.0
    assert st1["damage_error"] > 0.0
    sig2, _, st2 = m.compute_state(1.6 * eps, st1)
    assert st2["damage_n"] > 0.0 and sig2[0] < (Ce(PS) @ (1.6 * eps))[0]


def test_implex_mismo_intento_mismo_estado():
    m = mat(integration="impl-ex")
    st = m.initial_point_state(H)
    m.set_pseudo_time_ratio(1.0)
    eps0 = np.array([1.0, -0.2, 0.4])
    for f in np.linspace(0.0, 4.0, 13)[1:]:
        _, _, st = m.compute_state(f * K0 * eps0, st)
    m.set_pseudo_time_ratio(0.5)
    a = m.compute_state(4.2 * K0 * eps0, st)
    b = m.compute_state(4.2 * K0 * eps0, st)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])
    assert a[2] == b[2]


def test_implex_fuera_de_un_paso_el_estado_queda_congelado():
    m = mat(integration="impl-ex")
    st = m.initial_point_state(H)
    m.set_pseudo_time_ratio(1.0)
    for e in np.linspace(0.0, 6.0, 31)[1:] * K0:
        _, _, st = m.compute_state(np.array([e, 0.0, 0.2 * e]), st)
    m.set_pseudo_time_ratio(None)
    eps = np.array([9.0 * K0, 0.0, 1.8 * K0])
    sig, _, nuevo = m.compute_state(eps, st)
    assert nuevo == st
    # El esfuerzo es el de la ley secante con la flexibilidad guardada en
    # damage_n (la grieta t no se formó: rígida).
    assert st["damage_t"] == 0.0
    c_n = st["damage_n"] * H / (E * (1 - st["damage_n"]))
    e = np.linalg.solve(N_spec(st["theta"]), eps)
    C = Ce(PS)
    w = (C[0, 0] * e[0] + C[0, 1] * e[1]) / (1.0 / c_n + C[0, 0] / H)
    assert local(st["theta"], sig)[0] == pytest.approx(w / c_n, rel=1e-12)
