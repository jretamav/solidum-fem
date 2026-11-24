"""``TensionCompressionDamage2D/3D``: punto material e IMPL-EX.

Bloque 1 (sólo ``d⁺``, sin ``f_c``): aceptación ``material_point`` e
``impl_ex`` de la spec ``TensionCompressionDamage2D`` (§1, §2, §6–§8) en lo que
no depende del daño de compresión: tracción uniaxial exacta, deformación plana, cortante
puro, cierre de grieta, compresión paralela a una grieta abierta, descarga
secante, potencial con ν = 0, trabajo en ciclos cerrados con ν = 0.2 (contra el
prototipo independiente), energía disipada, irreversibilidad, tangentes por
diferencias finitas (también la global ensamblada con partición de signo
mixto), restricción del 3D, invariancia de unidades, kernel por lotes y la
integración IMPL-EX. Los hallazgos H1–H5 y P4 de la verificación del
2026-11-18 (``informes/teorico/verificacion_teorica_2026-11-18_modelo2_bloque1.md``)
están cubiertos.

Bloque 2 (con ``f_c`` y ``G_c``, spec §3; sección final): parámetros y estado
inicial, curva de compresión uniaxial contra la parábola de Feenstra escrita
aquí aparte (pico ``f_c`` en ``κ_c``, área ``G_c/h`` con el tope descontado en
forma cerrada), pico biaxial ``β·f_c``, compresión hidrostática sin daño,
medida de Drucker-Prager contra una evaluación independiente, irreversibilidad
de los dos daños, tangente implícita con ``d⁻`` (también con el término de
``σ̄_zz`` en deformación plana), restricción del 3D, unidades, lotes y
IMPL-EX en compresión.

Los daños fijos se imponen con la integración IMPL-EX fuera de un paso (sin
razón de pseudo-tiempo), que congela el estado: ``σ`` y la tangente son los de
los daños del estado (spec §8).

Concreto usual en SI: E = 30 GPa, ν = 0.2, f_t = 3 MPa, G_f = 100 N/m
(κ_0 = 1e-4, h_max = 0.667 m).
"""
import math

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.constants import CRACK_BAND_DAMAGE_MAX
from solidum.elements.solid_2d import Quad4
from solidum.elements.solid_3d import Hex8
from solidum.materials.damage_tension_compression_2d import TensionCompressionDamage2D
from solidum.materials.damage_tension_compression_3d import TensionCompressionDamage3D

E, NU, FT, GF = 30.0e9, 0.2, 3.0e6, 100.0
K0 = FT / E
H = 0.025
G = E / (2.0 * (1.0 + NU))


def mat2(hypothesis="plane_stress", nu=NU, E=E, **kw):
    return TensionCompressionDamage2D(E, nu, FT, GF, hypothesis=hypothesis, **kw)


def mat3(nu=NU, E=E, **kw):
    return TensionCompressionDamage3D(E, nu, FT, GF, **kw)


def mat(kind, **kw):
    """``kind`` ∈ {'plane_stress', 'plane_strain', '3d'}."""
    return mat3(**kw) if kind == "3d" else mat2(kind, **kw)


def frozen(m, d_t, d_c=0.0, h=H):
    """Estado con los daños fijos (``m`` con IMPL-EX y sin paso en curso)."""
    st = m.initial_point_state(h)
    st.update(damage_t=float(d_t), damage_c=float(d_c))
    return st


def uniaxial(kind, e, nu=NU):
    """Tracción (e > 0) o compresión uniaxial de esfuerzo; en deformación plana
    con ``σ_yy = 0`` y ``ε_zz = 0``."""
    if kind == "3d":
        return np.array([e, -nu * e, -nu * e, 0.0, 0.0, 0.0])
    if kind == "plane_strain":
        return np.array([e, -nu / (1.0 - nu) * e, 0.0])
    return np.array([e, -nu * e, 0.0])


def assoc_3d(eps2, nu, hypothesis):
    """Estado 3D asociado a un estado plano (σ̄_zz = 0 o ε_zz = 0)."""
    c_zz = -nu / (1.0 - nu) if hypothesis == "plane_stress" else 0.0
    return np.array([eps2[0], eps2[1], c_zz * (eps2[0] + eps2[1]), eps2[2], 0.0, 0.0])


def dissipated(kappa, h, E=E, ft=FT, gf=GF):
    """D(κ) del modelo 1 (spec IsotropicDamageCrackBand2D §2), por unidad de volumen."""
    k0 = ft / E
    a = 1.0 / (gf / (h * ft) - k0 / 2.0)
    x = math.exp(-a * (kappa - k0))
    return 0.5 * ft * k0 + ft / a * (1.0 - x) - 0.5 * ft * x * kappa


def exp_damage(kappa, h):
    a = 1.0 / (GF / (h * FT) - K0 / 2.0)
    return 1.0 - K0 / kappa * math.exp(-a * (kappa - K0))


# ---------------------------------------------------------------------------
# Parámetros y contrato
# ---------------------------------------------------------------------------

def test_parametros_y_estado_inicial():
    for m in (mat2(), mat2("plane_strain"), mat3()):
        assert m.uses_characteristic_length
        assert m.IS_SYMMETRIC is False
        assert m.SECANT_UNLOADING is True
        assert m.PRIMARY_STATE_VAR == "damage_t"
        assert m.damage_max == CRACK_BAND_DAMAGE_MAX
        assert m.integration == "implicit" and not m.uses_pseudo_time
        assert m.initial_point_state(H) == {
            "kappa_t": K0, "kappa_c": 0.0, "damage_t": 0.0, "damage_c": 0.0, "h": H,
            "kappa_t_prev": K0, "kappa_c_prev": 0.0, "damage_error": 0.0}
        assert m.alpha_for(H) == pytest.approx(1.0 / (GF / (H * FT) - K0 / 2.0), rel=1e-14)
    # IMPL-EX: también asimétrica (la partición cambia con la deformación, §9).
    m = mat2(integration="impl-ex")
    assert m.uses_pseudo_time and m.IS_SYMMETRIC is False


def test_parametros_invalidos():
    with pytest.raises(ValueError, match="nu"):
        mat2(nu=-0.1)
    with pytest.raises(ValueError, match="nu"):
        mat3(nu=-0.1)
    with pytest.raises(ValueError, match="nu"):
        mat2(nu=0.5)
    with pytest.raises(ValueError, match="f_t"):
        TensionCompressionDamage2D(E, NU, 0.0, GF)
    with pytest.raises(ValueError, match="G_f"):
        TensionCompressionDamage3D(E, NU, FT, -1.0)
    with pytest.raises(ValueError, match="damage_max"):
        mat2(damage_max=1.0)
    with pytest.raises(ValueError, match="integration"):
        mat3(integration="explicit")
    with pytest.raises(ValueError):
        mat2("plane_strian")


@pytest.mark.parametrize("h", [2.0 * E * GF / FT**2, float("nan"), 0.0])
def test_h_no_admisible(h):
    for m in (mat2(), mat3()):
        with pytest.raises(ValueError, match="h"):
            m.initial_point_state(h)


def test_dano_sin_h_avisa():
    m = mat2()
    m.compute_state(uniaxial("plane_stress", 0.5 * K0), None)     # elástico: admitido
    with pytest.raises(ValueError, match="sin longitud característica"):
        m.compute_state(uniaxial("plane_stress", 2.0 * K0), None)


# ---------------------------------------------------------------------------
# Ley de esfuerzo (§1)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kind", ["plane_stress", "3d"])
@pytest.mark.parametrize("d", [0.0, 0.5, 0.9])
def test_traccion_uniaxial_exacta(kind, d):
    """σ = (1 − d⁺)·E·ε con ε_⊥ = −ν·ε y las demás componentes nulas."""
    m = mat(kind, integration="impl-ex")
    st = frozen(m, d)
    e = 3.7e-4
    sig, _, _ = m.compute_state(uniaxial(kind, e), st)
    assert sig[0] == pytest.approx((1.0 - d) * E * e, rel=1e-12)
    np.testing.assert_allclose(sig[1:], 0.0, atol=1e-12 * E * e)


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
def test_traccion_uniaxial_con_dano_activo(kind):
    """Implícita: κ⁺ = ε_eq de Rankine y d⁺ de la ley exponencial (forma
    cerrada); en deformación plana ε_eq = ε/(1 − ν²)."""
    m = mat(kind)
    st = m.initial_point_state(H)
    f = 1.0 / (1.0 - NU**2) if kind == "plane_strain" else 1.0
    for e in np.linspace(0.0, 6.0 * K0, 61)[1:]:
        sig, _, st = m.compute_state(uniaxial(kind, e), st)
        kappa = max(K0, f * e)
        assert st["kappa_t"] == pytest.approx(kappa, rel=1e-13)
        d = exp_damage(kappa, H) if kappa > K0 * (1 + 1e-12) else 0.0
        assert st["damage_t"] == pytest.approx(d, rel=1e-12, abs=1e-15)
        assert sig[0] == pytest.approx((1.0 - d) * E * f * e, rel=1e-11)


@pytest.mark.parametrize("d", [0.0, 0.5, 0.9])
def test_deformacion_plana_uniaxial(d):
    """σ_yy = 0, ε_zz = 0: σ_xx = (1 − d⁺)·E·ε/(1 − ν²) y σ_zz = ν·σ_xx
    (valor independiente del 3D)."""
    m = mat2("plane_strain", integration="impl-ex")
    st = frozen(m, d)
    e = 2.9e-4
    sig, _, st1 = m.compute_state(uniaxial("plane_strain", e), st)
    sxx = (1.0 - d) * E * e / (1.0 - NU**2)
    assert sig[0] == pytest.approx(sxx, rel=1e-12)
    assert abs(sig[1]) <= 1e-12 * sxx and sig[2] == 0.0
    assert m.out_of_plane_stress(sig, st1) == pytest.approx(NU * sxx, rel=1e-12)


def test_cortante_puro_forma_cerrada():
    """σ = G·γ·[−(d⁺−d⁻)/2, −(d⁺−d⁻)/2, 1 − (d⁺+d⁻)/2] con daños fijos y, con
    daño activo, d⁺ empieza en τ = G·γ = f_t (blinda γ = 2·ε_xy)."""
    m = mat2(integration="impl-ex")
    dp, dm = 0.8, 0.25
    g = 2.3e-4
    sig, _, _ = m.compute_state(np.array([0.0, 0.0, g]), frozen(m, dp, dm))
    expected = G * g * np.array([-(dp - dm) / 2, -(dp - dm) / 2, 1.0 - (dp + dm) / 2])
    np.testing.assert_allclose(sig, expected, rtol=1e-12, atol=1e-12 * G * g)
    mi = mat2()
    g0 = FT / G
    _, _, st = mi.compute_state(np.array([0.0, 0.0, 0.999 * g0]), mi.initial_point_state(H))
    assert st["damage_t"] == 0.0
    _, _, st = mi.compute_state(np.array([0.0, 0.0, 1.001 * g0]), st)
    assert st["damage_t"] > 0.0
    assert st["kappa_t"] == pytest.approx(1.001 * g0 * G / E, rel=1e-12)


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
def test_cierre_de_grieta(kind):
    """Tracción hasta d⁺ ≥ 0.9 y después compresión: la grieta cerrada
    recupera E y ν, y la compresión no hace crecer d⁺."""
    m = mat(kind)
    st = m.initial_point_state(H)
    e = K0
    while st["damage_t"] < 0.9:
        e *= 1.05
        _, _, st = m.compute_state(uniaxial(kind, e), st)
    d = st["damage_t"]
    ec = 4.0e-4
    sig, _, st_c = m.compute_state(uniaxial(kind, -ec), st)
    f = 1.0 / (1.0 - NU**2) if kind == "plane_strain" else 1.0
    assert sig[0] == pytest.approx(-E * f * ec, rel=1e-12)
    assert st_c["damage_t"] == d and st_c["kappa_t"] == st["kappa_t"]
    # Compresión con todos los esfuerzos principales negativos: σ = C_e·ε y
    # tangente C_e (rigidez y Poisson intactos).
    n = m.STRAIN_DIM
    eps = -np.r_[3.0e-4, 2.0e-4, 1.5e-4, 0.4e-4, -0.3e-4, 0.2e-4] if n == 6 \
        else -np.r_[3.0e-4, 2.0e-4, 0.4e-4]
    sig, C, _ = m.compute_state(eps, st)
    Ce = m.elastic_base.C
    np.testing.assert_allclose(C, Ce, rtol=0, atol=1e-14 * E)
    np.testing.assert_allclose(sig, Ce @ eps, rtol=1e-13)


@pytest.mark.parametrize("en", [0.0, 1.0e-3, 2.0e-3, 4.0e-3, 5.0e-3, 6.0e-3])
def test_compresion_paralela_a_grieta_abierta(en):
    """Caveat de §1: σ_tt = E/(1 − ν²)·(ε_t + ν·ε_n) mientras sea negativa; desde
    ε_n = |ε_t|/ν (5·10⁻³ con ν = 0.2) la compresión paralela se pierde."""
    m = mat2(integration="impl-ex")
    et = -1.0e-3
    sig, _, _ = m.compute_state(np.array([en, et, 0.0]), frozen(m, m.damage_max))
    s_tt = E / (1.0 - NU**2) * (et + NU * en)
    expected = s_tt if s_tt < 0.0 else (1.0 - m.damage_max) * s_tt
    assert sig[1] == pytest.approx(expected, rel=1e-5, abs=1e-5 * E * abs(et))
    if en >= abs(et) / NU:
        assert abs(sig[1]) <= 1e-5 * E * abs(et)


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
@pytest.mark.parametrize("dp,dm", [(0.9, 0.0), (0.6, 0.3)])
def test_descarga_radial_secante(kind, dp, dm):
    """σ(s·ε) = s·σ(ε) y el trabajo devuelto en la descarga radial es ½·σ:ε."""
    m = mat(kind, integration="impl-ex")
    st = frozen(m, dp, dm)
    rng = np.random.default_rng(5)
    for _ in range(20):
        eps = rng.standard_normal(m.STRAIN_DIM) * 3.0e-4
        sig = m.compute_state(eps, st)[0]
        for s in (0.3, 2.5):
            np.testing.assert_allclose(m.compute_state(s * eps, st)[0], s * sig,
                                       rtol=1e-12, atol=1e-12 * np.abs(sig).max())
        ss = np.linspace(1.0, 0.0, 11)
        sigs = [m.compute_state(s * eps, st)[0] for s in ss]
        W = sum(0.5 * (sigs[k] + sigs[k + 1]) @ ((ss[k] - ss[k + 1]) * eps) for k in range(10))
        assert W == pytest.approx(0.5 * sig @ eps, rel=1e-12)


# ---------------------------------------------------------------------------
# Energía en ciclos cerrados con daños fijos (§1)
# ---------------------------------------------------------------------------

def _cycles(seed, cases, n_cycles):
    """Ciclos cerrados ε(t) = c + A·[cos t, sin t, sin 2t] (el generador de
    ``informes/prototipos_2026-11-18/modelo2_bloque1/ref_ciclos.py``)."""
    rng = np.random.default_rng(seed)
    out = []
    for dims, dp, dm in cases:
        for _ in range(n_cycles):
            A = rng.normal(size=(dims, 3)) * 1e-4
            c = rng.normal(size=dims) * 1e-4 * rng.uniform(0, 2)
            out.append((dims, dp, dm, A, c))
    return out


def mat_mpa(dims, nu):
    """Concreto del prototipo en N, mm, MPa (E = 30 000) con IMPL-EX, para
    imponer los daños fijos."""
    cls = TensionCompressionDamage3D if dims == 6 else TensionCompressionDamage2D
    return cls(30000.0, nu, 3.0, 0.1, integration="impl-ex")


def _cycle_work(m, dp, dm, A, c, n_pts):
    """(∮σ:dε por la regla del trapecio, máx ½·σ:ε en el ciclo)."""
    st = frozen(m, dp, dm, h=10.0)
    ts = np.linspace(0.0, 2.0 * np.pi, n_pts)
    W, psi = 0.0, 0.0
    e0 = c + A @ np.array([1.0, 0.0, 0.0])
    s0 = m.compute_state(e0, st)[0]
    psi = 0.5 * s0 @ e0
    for t in ts[1:]:
        e1 = c + A @ np.array([np.cos(t), np.sin(t), np.sin(2 * t)])
        s1 = m.compute_state(e1, st)[0]
        W += 0.5 * (s0 + s1) @ (e1 - e0)
        psi = max(psi, 0.5 * s1 @ e1)
        e0, s0 = e1, s1
    return W, psi


@pytest.mark.parametrize("dp,dm", [(0.9, 0.0), (0.9, 0.3)])
def test_potencial_con_nu_cero(dp, dm):
    """Con ν = 0 la ley deriva de un potencial: ∮σ:dε = 0 (16 001 puntos)."""
    cases = [(6, dp, dm), (3, dp, dm)]
    for dims, _, _, A, c in _cycles(11, cases, 2):
        W, psi = _cycle_work(mat_mpa(dims, 0.0), dp, dm, A, c, 16001)
        assert abs(W) <= 1e-6 * psi


# Trabajo |W|/ψ_max con ν = 0.2 calculado con el prototipo independiente de la
# ley (a) (``informes/prototipos_2026-11-18/modelo2_bloque1/ref_ciclos.py``,
# sobre ``prototipos_2026-11-17/modelo2_ley_a/check_ley_a.py``): semilla 2026,
# 20 ciclos por caso, 1001 puntos; casos (3D, 0.5, 0), (3D, 0.9, 0.3),
# (esfuerzo plano, 0.5, 0), (esfuerzo plano, 0.9, 0.3).
REF_CICLOS = [
    -3.124309835085800e-02, -2.549906746794207e-02, 2.648791663499482e-02, 3.817473674710223e-03,
    -8.433970861942823e-03, 1.896203484690645e-03, -2.070673006997850e-02, 3.383288253679560e-02,
    -8.497423465173945e-04, -4.752348926240309e-03, -3.933547939534331e-02, -3.087584164007929e-02,
    -2.513348632305987e-02, -9.988157927262833e-03, 2.577973699690617e-02, 3.848520040979180e-02,
    -1.073683016330366e-02, -3.643919098380112e-04, 1.377121695955287e-03, 9.529054253499455e-02,
    1.705765141950428e-02, 2.922609055687295e-02, -9.743744378158973e-03, 1.297730727119893e-01,
    -8.569980486184383e-04, 9.604029397838567e-04, 1.501380591322132e-02, -2.225337852557474e-01,
    5.555696187379441e-02, -1.246722732543632e-02, 5.213577225398579e-02, -3.964754493453915e-02,
    2.641038494170357e-03, 1.638119633674461e-01, 3.050901229402141e-02, -8.381581287331673e-02,
    -8.861177539470461e-02, -8.920822245662077e-02, -2.723592629090596e-02, 1.072658717461625e-01,
    -2.212292873578353e-02, 5.149938641782548e-02, 1.527829515775984e-02, -1.300689975696316e-01,
    5.715248654030142e-02, 3.760230636589119e-02, 4.742849766364039e-03, 1.888004172239924e-02,
    -3.929417689748549e-02, 2.116526050146672e-02, -1.478452796566865e-01, -3.457957251446932e-03,
    7.899234457097160e-03, 3.036031021258754e-02, 5.857583362927205e-02, 2.457934565647494e-02,
    5.491080150277897e-02, 4.556904705198019e-02, 4.275229021722378e-02, 3.913271222478978e-03,
    2.231391682458322e-02, 3.394403563157231e-02, 1.021066295971252e-02, -1.814975820911467e-01,
    2.445744530647856e-02, 5.283123693146154e-02, 2.500834306775047e-02, 3.840111420545181e-02,
    -3.202197224748843e-01, 1.105265206122865e-02, -1.095919262598269e-04, 4.674087621003856e-02,
    5.401680058332092e-03, -3.766269173407230e-02, 5.954340143324753e-02, -1.095906435387723e-01,
    1.187686425359675e-01, 7.635729446690077e-02, 1.506023053132387e-01, 6.291237176122597e-02,
]


def test_trabajo_en_ciclos_cerrados_con_poisson():
    """Defecto documentado en §1 (con ν ≠ 0 no hay potencial): trabajo de los
    dos signos, mediana de |W|/ψ_max entre 1 y 10 %, y cada ciclo igual al
    prototipo independiente."""
    cases = [(6, 0.5, 0.0), (6, 0.9, 0.3), (3, 0.5, 0.0), (3, 0.9, 0.3)]
    m3, m2 = mat_mpa(6, NU), mat_mpa(3, NU)
    r = []
    for dims, dp, dm, A, c in _cycles(2026, cases, 20):
        W, psi = _cycle_work(m3 if dims == 6 else m2, dp, dm, A, c, 1001)
        r.append(W / psi)
    r = np.array(r)
    np.testing.assert_allclose(r, REF_CICLOS, rtol=0, atol=1e-8)
    assert 0.01 < np.median(np.abs(r)) < 0.10
    assert (r > 0).any() and (r < 0).any()


# ---------------------------------------------------------------------------
# Energía disipada en tracción (§2) e irreversibilidad
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kind,h", [("plane_stress", 0.005), ("plane_stress", 0.025),
                                    ("plane_stress", 0.100), ("plane_stress", 0.400),
                                    ("3d", 0.005), ("3d", 0.400), ("plane_strain", 0.025)])
def test_energia_disipada_en_traccion(kind, h):
    """∫σ_xx dε_xx − ½σ_xx·ε_xx = D(κ) del modelo 1 en estados intermedios y
    > 0.997·G_f/h al final; en deformación plana con el lateral libre, por
    (1 − ν²) (spec §1)."""
    m = mat(kind)
    st = m.initial_point_state(h)
    f = 1.0 / (1.0 - NU**2) if kind == "plane_strain" else 1.0
    a = m.alpha_for(h)
    e_end = (K0 + 25.0 / a) / f
    path = np.concatenate([np.linspace(0.0, K0 / f, 51), np.linspace(K0 / f, e_end, 40001)[1:]])
    work, s_prev, e_prev, checked = 0.0, 0.0, 0.0, 0
    for i, e in enumerate(path):
        s, _, st = m.compute_state(uniaxial(kind, e), st)
        work += 0.5 * (s[0] + s_prev) * (e - e_prev)
        s_prev, e_prev = s[0], e
        if i % 4000 == 0 and 0.0 < st["damage_t"] < m.damage_max:
            checked += 1
            assert work - 0.5 * s[0] * e == pytest.approx(dissipated(f * e, h) / f, rel=1e-6)
    assert checked >= 3
    assert work - 0.5 * s_prev * e_prev > 0.997 * GF / h / f


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_dano_irreversible(kind, integration):
    """Carga, descarga, cierre en compresión y reapertura: d⁺ nunca decrece y
    la compresión no lo hace crecer."""
    m = mat(kind, integration=integration)
    m.set_pseudo_time_ratio(1.0)
    st = m.initial_point_state(H)
    rng = np.random.default_rng(17)
    n = m.STRAIN_DIM
    base = rng.standard_normal(n) * 0.2
    base[0] = 1.0
    path = np.r_[np.linspace(0, 4, 21)[1:], np.linspace(4, -6, 26)[1:], np.linspace(-6, 7, 31)[1:]]
    d_prev = 0.0
    for f in path:
        eps = f * K0 * base
        _, _, st = m.compute_state(eps, st)
        assert st["damage_t"] >= d_prev
        if f < 0.0 and np.all(np.linalg.eigvalsh(_tensor(m.elastic_base.C @ eps)) <= 0.0):
            assert st["damage_t"] == d_prev
        assert st["damage_c"] == 0.0
        d_prev = st["damage_t"]
    assert d_prev > 0.5


def _tensor(s):
    if len(s) == 3:
        return np.array([[s[0], s[2]], [s[2], s[1]]])
    return np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])


# ---------------------------------------------------------------------------
# Tangentes (§7)
# ---------------------------------------------------------------------------

def _fd(m, eps, st, rel=1e-6):
    n = len(eps)
    C_fd = np.zeros((n, n))
    for j in range(n):
        d = rel * np.abs(eps).max()
        ep, em = eps.copy(), eps.copy()
        ep[j] += d
        em[j] -= d
        C_fd[:, j] = (m.compute_state(ep, st)[0] - m.compute_state(em, st)[0]) / (2 * d)
    return C_fd


def _dobles(m):
    """Deformaciones con esfuerzos principales efectivos dobles del mismo signo.

    En el plano un doble es isótropo (σ̄⁺ = σ̄ o 0 con cualquier giro); en 3D,
    σ̄ = R·diag(a, a, b)·Rᵀ girado, con el doble positivo y negativo, y el
    triple.
    """
    Ce = m.elastic_base.C
    if m.STRAIN_DIM == 3:
        return [np.linalg.solve(Ce, s) for s in (np.array([6.0e6, 6.0e6, 0.0]),
                                                  np.array([-6.0e6, -6.0e6, 0.0]))]
    R, _ = np.linalg.qr(np.random.default_rng(37).standard_normal((3, 3)))
    out = []
    for lam in ((6.0e6, 6.0e6, -9.0e6), (-6.0e6, -6.0e6, 9.0e6), (3.0e6, 3.0e6, 3.0e6)):
        T = R @ np.diag(lam) @ R.T
        s = np.array([T[0, 0], T[1, 1], T[2, 2], T[0, 1], T[1, 2], T[0, 2]])
        out.append(np.linalg.solve(Ce, s))
    return out


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
def test_tangente_dano_fijo_por_diferencias_finitas(kind):
    m = mat(kind, integration="impl-ex")
    st = frozen(m, 0.7, 0.2)
    rng = np.random.default_rng(23)
    casos = [rng.standard_normal(m.STRAIN_DIM) * 3.0e-4 for _ in range(30)] + _dobles(m)
    for eps in casos:
        C = m.compute_state(eps, st)[1]
        C_fd = _fd(m, eps, st)
        assert np.abs(C - C_fd).max() <= 1e-6 * np.abs(C_fd).max()


@pytest.mark.parametrize("d_c", [0.0, 0.4])
@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
def test_tangente_implicita_en_carga_por_diferencias_finitas(kind, d_c):
    """Carga activa con d⁺: la tangente lleva −σ̄⁺⊗∂d⁺/∂ε (sin la arista de
    Rankine ni esfuerzos principales efectivos nulos, §10); también con un
    d⁻ fijo ≠ 0."""
    m = mat(kind)
    rng = np.random.default_rng(29)
    n = m.STRAIN_DIM
    checked = 0
    for _ in range(40):
        eps = rng.standard_normal(n) * 1.0e-4
        eps[0] += 3.0e-4
        s_eff = m.elastic_base.C @ eps
        lam = np.linalg.eigvalsh(_tensor(s_eff))
        if lam[-1] <= 1.2 * FT or (len(lam) > 1 and lam[-1] - lam[-2] < 0.05 * lam[-1]):
            continue
        st = m.initial_point_state(H)
        st["damage_c"] = d_c
        st["kappa_t"] = 0.97 * lam[-1] / E            # carga activa en ±δ
        st["damage_t"] = exp_damage(st["kappa_t"], H)
        _, C, st1 = m.compute_state(eps, st)
        assert st1["damage_t"] > st["damage_t"]
        C_fd = _fd(m, eps, st, rel=1e-7)
        assert np.abs(C - C_fd).max() <= 1e-6 * np.abs(C_fd).max()
        checked += 1
    assert checked >= 10


def test_escalon_nulo_en_cero():
    """H(0) = 0 (spec §7): un esfuerzo principal efectivo nulo cuenta como no
    positivo en Q⁺ (el esfuerzo no cambia; sólo la tangente sobre el quiebre)."""
    from solidum.materials._spectral_split import spectral_split_2d, spectral_split_3d
    sp, Q, _, _ = spectral_split_2d(np.array([5.0, 0.0, 0.0]))
    np.testing.assert_array_equal(sp, [5.0, 0.0, 0.0])
    np.testing.assert_array_equal(Q, np.diag([1.0, 0.0, 1.0]))
    sp, Q, _, _ = spectral_split_3d(np.array([5.0, 0.0, 0.0, 0.0, 0.0, 0.0]))
    np.testing.assert_array_equal(sp, [5.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    np.testing.assert_allclose(Q, np.diag([1.0, 0.0, 0.0, 1.0, 0.0, 1.0]), rtol=0, atol=1e-15)
    for split, n in ((spectral_split_2d, 3), (spectral_split_3d, 6)):
        sp, Q, _, _ = split(np.zeros(n))
        assert not sp.any() and not Q.any()


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
def test_simetria_de_la_tangente(kind):
    """Con daño fijo: simétrica con ν = 0 y asimétrica con ν = 0.2 (§7, §9)."""
    rng = np.random.default_rng(31)
    for nu, sym in ((0.0, True), (NU, False)):
        m = mat(kind, nu=nu, integration="impl-ex")
        st = frozen(m, 0.8, 0.1)
        asym = 0.0
        for _ in range(20):
            C = m.compute_state(rng.standard_normal(m.STRAIN_DIM) * 3e-4, st)[1]
            asym = max(asym, np.abs(C - C.T).max() / np.abs(C).max())
        if sym:
            assert asym <= 1e-14
        else:
            assert asym > 1e-2


# ---------------------------------------------------------------------------
# Restricción del 3D, unidades, lotes y post-proceso
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("hypothesis", ["plane_stress", "plane_strain"])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_restriccion_exacta_del_3d(hypothesis, integration):
    m2, m3 = mat2(hypothesis, integration=integration), mat3(integration=integration)
    for m in (m2, m3):
        m.set_pseudo_time_ratio(1.0)
    for base in (np.array([1.0, -0.15, 0.35]), np.array([-0.3, 1.0, -0.6]),
                 np.array([0.7, 0.5, 0.9])):
        s2, s3 = m2.initial_point_state(H), m3.initial_point_state(H)
        s2["damage_c"] = s3["damage_c"] = 0.3
        # Carga, descarga, compresión (σ̄_zz < 0) y recarga.
        for f in (0.5, 1.5, 3.0, 2.0, -4.0, 1.0, 4.0, 6.0, -2.0, 8.0):
            e2 = f * K0 * base
            sig2, _, s2 = m2.compute_state(e2, s2)
            sig3, _, s3 = m3.compute_state(assoc_3d(e2, NU, hypothesis), s3)
            for k in s2:
                # damage_error es la diferencia de dos daños casi iguales: absoluta.
                tol = dict(rel=0, abs=1e-12) if k == "damage_error" \
                    else dict(rel=1e-12, abs=1e-18)
                assert s2[k] == pytest.approx(s3[k], **tol), k
            scale = np.abs(sig3).max()
            np.testing.assert_allclose(sig2, sig3[[0, 1, 3]], rtol=1e-12, atol=1e-12 * scale)
            np.testing.assert_allclose(sig3[[4, 5]], 0.0, atol=1e-12 * scale)
            szz = m2.out_of_plane_stress(sig2, s2)
            assert szz == pytest.approx(sig3[2], rel=1e-9, abs=1e-12 * scale)
        assert s2["damage_t"] > 0.5


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_invariancia_de_unidades(kind, integration):
    """El mismo caso en Pa–m y en MPa–mm con todos los parámetros ≠ 1."""
    par_si = dict(E=31.7e9, nu=0.18, f_t=2.9e6, G_f=93.0)
    par_mm = dict(E=31.7e3, nu=0.18, f_t=2.9, G_f=0.093)
    h_si, h_mm = 0.037, 37.0

    def build(p):
        if kind == "3d":
            return TensionCompressionDamage3D(**p, integration=integration)
        return TensionCompressionDamage2D(**p, hypothesis=kind, integration=integration)

    a, b = build(par_si), build(par_mm)
    for m in (a, b):
        m.set_pseudo_time_ratio(1.0)
    sa, sb = a.initial_point_state(h_si), b.initial_point_state(h_mm)
    rng = np.random.default_rng(41)
    base = rng.standard_normal(a.STRAIN_DIM) * 0.3
    base[0] = 1.0
    k0 = 2.9e6 / 31.7e9
    Wa = Wb = 0.0
    ea_prev = np.zeros(a.STRAIN_DIM)
    siga_prev = sigb_prev = np.zeros(a.STRAIN_DIM)
    for f in np.r_[np.linspace(0, 6, 61)[1:], np.linspace(6, -3, 19)[1:], np.linspace(-3, 9, 25)[1:]]:
        eps = f * k0 * base
        siga, _, sa = a.compute_state(eps, sa)
        sigb, _, sb = b.compute_state(eps, sb)
        np.testing.assert_allclose(sigb, 1e-6 * siga, rtol=1e-12, atol=1e-12 * np.abs(sigb).max())
        for k in ("kappa_t", "damage_t", "damage_error"):
            assert sb[k] == pytest.approx(sa[k], rel=1e-12, abs=1e-15), k
        Wa += 0.5 * (siga + siga_prev) @ (eps - ea_prev)
        Wb += 0.5 * (sigb + sigb_prev) @ (eps - ea_prev)
        ea_prev, siga_prev, sigb_prev = eps, siga, sigb
    assert sa["damage_t"] > 0.5
    assert Wb == pytest.approx(1e-6 * Wa, rel=1e-12)


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_kernel_por_lotes_reproduce_compute_state(kind, integration):
    m = mat(kind, integration=integration)
    n = m.STRAIN_DIM
    names = list(m.STATE_SCHEMA)
    kernel = m.batch_kernel()
    C = m.batch_matrix()
    st = m.initial_point_state(H)
    st["damage_c"] = 0.2
    rng = np.random.default_rng(3)
    for i in range(30):
        # IMPL-EX: razón variable por paso y, cada cinco, «fuera de un paso».
        m.set_pseudo_time_ratio(None if i % 5 == 4 else 0.5 + 0.1 * i)
        params = m.batch_params()
        eps = (rng.standard_normal(n) + np.r_[2.0, np.zeros(n - 1)]) * 1.5e-4 * (1 if i % 7 else -1)
        sig_ref, C_ref, st_new = m.compute_state(eps, st)
        S_old = np.array([st[k] for k in names])
        S_new, sig, C_out, flag = (np.zeros(len(names)), np.zeros(n), np.zeros((n, n)),
                                   np.zeros(1, np.int8))
        C_t = kernel(eps, S_old, S_new, params, C, sig, C_out, flag)
        np.testing.assert_array_equal(sig, sig_ref)
        np.testing.assert_array_equal(C_t, C_ref)
        np.testing.assert_array_equal(S_new, [st_new[k] for k in names])
        assert flag[0] == 0
        st = st_new
    assert st["damage_t"] > 0.0


def test_kernel_por_lotes_marca_dano_sin_h():
    m = mat2()
    names = list(m.STATE_SCHEMA)
    S_old = np.array([m.initial_state()[k] for k in names])
    eps = uniaxial("plane_stress", 3.0 * K0)
    flag = np.zeros(1, np.int8)
    m.batch_kernel()(eps, S_old, np.zeros(len(names)), m.batch_params(), m.batch_matrix(),
                     np.zeros(3), np.zeros((3, 3)), flag)
    assert flag[0] == 1
    with pytest.raises(ValueError, match="sin longitud característica"):
        m.batch_report(1)


def test_esfuerzo_fuera_del_plano_por_lotes():
    rng = np.random.default_rng(43)
    for hyp in ("plane_stress", "plane_strain"):
        m = mat2(hyp, integration="impl-ex")
        names = list(m.STATE_SCHEMA)
        sig, S, ref = [], [], []
        for _ in range(25):
            st = frozen(m, rng.uniform(0, 0.99), rng.uniform(0, 0.5))
            s = m.compute_state(rng.standard_normal(3) * 3e-4, st)[0]
            sig.append(s)
            S.append([st[k] for k in names])
            ref.append(m.out_of_plane_stress(s, st))
        out = m.batch_out_of_plane_stress(np.array(sig), np.array(S))
        np.testing.assert_allclose(out, ref, rtol=1e-14, atol=0)
        if hyp == "plane_stress":
            assert not np.any(out)


# ---------------------------------------------------------------------------
# IMPL-EX (§8)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
def test_implex_tangente_con_danos_congelados(kind):
    """Dentro de un paso los daños extrapolados no dependen de ε: la tangente
    es C^d(d̃⁺, d̃⁻), igual a diferencias finitas y a la de daño fijo."""
    m = mat(kind, integration="impl-ex")
    n = m.STRAIN_DIM
    base = np.r_[1.0, 0.3 * np.ones(n - 1)] * np.r_[1.0, -0.5, 0.7, 0.2, -0.4, 0.3][:n]
    m.set_pseudo_time_ratio(1.0)
    st = m.initial_point_state(H)
    for f in np.linspace(0.0, 3.0, 16)[1:]:
        _, _, st = m.compute_state(f * K0 * base / base[0], st)
    assert st["damage_t"] > 0.0
    rng = np.random.default_rng(47)
    for _ in range(10):
        eps = 3.3 * K0 * base / base[0] + rng.standard_normal(n) * 5e-5
        sig, C, st1 = m.compute_state(eps, st)
        assert np.abs(C - _fd(m, eps, st)).max() <= 1e-6 * np.abs(C).max()
        m.set_pseudo_time_ratio(None)
        sig_f, C_f, _ = m.compute_state(eps, frozen(m, st1["damage_t"], st1["damage_c"]))
        m.set_pseudo_time_ratio(1.0)
        np.testing.assert_array_equal(C, C_f)
        np.testing.assert_array_equal(sig, sig_f)


def _camino(n, integration, kappa_de_t):
    m = mat2(integration=integration)
    st = m.initial_point_state(H)
    m.set_pseudo_time_ratio(1.0)
    out = []
    for t in np.linspace(0.0, 1.0, n + 1)[1:]:
        sig, _, st = m.compute_state(uniaxial("plane_stress", kappa_de_t(t)), st)
        out.append(sig[0])
    return np.array(out)


@pytest.mark.parametrize("forma, orden", [("lineal", 1.0), ("cuadratica", 2.0)])
def test_implex_convergencia_a_la_implicita(forma, orden):
    """Error máximo de σ frente a la implícita con pasos Δ, Δ/2, Δ/4, Δ/8:
    primer orden con el quiebre de κ(t) en el inicio del daño y segundo
    donde κ(t) es suave (como el modelo 1)."""
    if forma == "lineal":
        def kappa(t):
            return K0 * (1.0 + 39.0 * t)
    else:
        def kappa(t):
            return K0 * (1.0 + 39.0 * t * t)
    errores = [np.abs(_camino(n, "impl-ex", kappa) - _camino(n, "implicit", kappa)).max()
               for n in (100, 200, 400, 800)]
    razones = np.array(errores[:-1]) / np.array(errores[1:])
    np.testing.assert_allclose(razones, 2.0 ** orden, rtol=0.15)


def test_implex_dano_monotono():
    """Carga, descarga parcial y recarga con una razón que hace pasarse a la
    extrapolación al empezar la descarga: d̃⁺ nunca decrece."""
    m = mat2(integration="impl-ex")
    st = m.initial_point_state(H)
    camino = np.r_[np.linspace(0.0, 8.0, 41)[1:], np.linspace(8.0, 3.0, 11)[1:],
                   np.linspace(3.0, 15.0, 31)[1:]] * K0
    d_prev = 0.0
    for k, e in enumerate(camino):
        m.set_pseudo_time_ratio(2.0 if k % 7 == 0 else 1.0)
        _, _, st = m.compute_state(uniaxial("plane_stress", e), st)
        assert st["damage_t"] >= d_prev
        assert st["damage_error"] >= 0.0
        d_prev = st["damage_t"]
    assert d_prev > 0.9


def test_implex_mismo_intento_mismo_estado():
    """Repetir el intento de paso (tras una bisección) con el mismo estado
    confirmado y la misma razón da el mismo estado de prueba."""
    m = mat3(integration="impl-ex")
    st = m.initial_point_state(H)
    m.set_pseudo_time_ratio(1.0)
    eps0 = np.array([1.0, -0.2, 0.1, 0.3, -0.1, 0.2])
    for f in np.linspace(0.0, 3.0, 13)[1:]:
        _, _, st = m.compute_state(f * K0 * eps0, st)
    m.set_pseudo_time_ratio(0.5)
    a = m.compute_state(3.2 * K0 * eps0, st)
    b = m.compute_state(3.2 * K0 * eps0, st)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])
    assert a[2] == b[2]


# ---------------------------------------------------------------------------
# Tangente global en una estructura con partición de signo mixto
# ---------------------------------------------------------------------------

def _malla_cortante(cls, integration, hypothesis="plane_stress"):
    """Bloque empotrado en x = 0 con un campo de desplazamiento de cortante y
    tracción no uniforme (N, mm): esfuerzos principales efectivos de signo
    mixto y daño activo (verificación del 2026-11-18, H3)."""
    from solidum.math.assembly import Assembler
    from solidum.mesh import structured_mesh
    E_, nu_, ft_, gf_ = 31700.0, 0.23, 3.1, 0.095
    k0 = ft_ / E_
    dim = len(cls.REFERENCE_NODES[0])
    if dim == 2:
        m = TensionCompressionDamage2D(E_, nu_, ft_, gf_, hypothesis=hypothesis,
                                       integration=integration)
        mesh = structured_mesh(cls, (5, 4), size=(50.0, 40.0), material=m, thickness=7.0)
    else:
        m = TensionCompressionDamage3D(E_, nu_, ft_, gf_, integration=integration)
        mesh = structured_mesh(cls, (3, 3, 2), size=(30.0, 30.0, 20.0), material=m)
    names = ["ux", "uy", "uz"][:dim]
    for nd in mesh.nodes_on("xmin"):
        for k in names:
            nd.fix_dof(k, 0.0)
    dom = mesh.domain
    dom.generate_equation_numbers(verbose=False)
    U = np.zeros(dom.total_dofs)
    for nd in dom.nodes.values():
        x = np.asarray(nd.coordinates[:dim], dtype=float)
        u = [4.0 * k0 * x[1] + 2.0 * k0 * x[0] * (1.0 + 0.02 * x[1]),
             -3.0 * k0 * x[1] + 1.5 * k0 * x[0]]
        if dim == 3:
            u.append(k0 * x[2] - 2.0 * k0 * x[1])
        for k, nm in enumerate(names):
            U[nd.dofs[nm]] = u[k]
    return Assembler(dom), U


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


@pytest.mark.parametrize("cls,hypothesis", [(Quad4, "plane_stress"), (Quad4, "plane_strain"),
                                            (Hex8, None)], ids=["Quad4-EP", "Quad4-DP", "Hex8"])
def test_tangente_global_implicita_con_particion_mixta(cls, hypothesis):
    """Desde el estado virgen, con daño activo y esfuerzos principales de signo
    mixto: K·δU = δF_int (el ensamblaje conserva la tangente asimétrica)."""
    asm, U = _malla_cortante(cls, "implicit", hypothesis or "plane_stress")
    peor, asim = _tangente_global_frente_a_df(asm, U)
    assert peor <= 1e-6
    assert asim > 1e-2


def test_tangente_global_impl_ex_con_dano_confirmado():
    """IMPL-EX tras seis pasos confirmados (daño d̃⁺ > 0): dentro del paso la
    tangente C^d(d̃⁺) asimétrica, ensamblada, es la de las fuerzas internas."""
    asm, U_obj = _malla_cortante(Quad4, "impl-ex")
    U_prev = np.zeros_like(U_obj)
    for k in range(1, 7):
        U_k = k / 6.0 * U_obj
        asm.prepare_all_steps(U_prev, pseudo_time_increment=1.0)
        asm.assemble_non_linear_system(U_k)
        asm.commit_all_states()
        U_prev = U_k
    asm.prepare_all_steps(U_prev, pseudo_time_increment=1.0)
    peor, asim = _tangente_global_frente_a_df(asm, 1.1 * U_obj)
    assert peor <= 1e-6
    assert asim > 1e-2


def test_implex_fuera_de_un_paso_el_estado_queda_congelado():
    m = mat2(integration="impl-ex")
    st = m.initial_point_state(H)
    m.set_pseudo_time_ratio(1.0)
    for e in np.linspace(0.0, 6.0, 31)[1:] * K0:
        _, _, st = m.compute_state(uniaxial("plane_stress", e), st)
    m.set_pseudo_time_ratio(None)
    eps = uniaxial("plane_stress", 9.0 * K0)
    sig, _, nuevo = m.compute_state(eps, st)
    assert nuevo == st
    assert sig[0] == pytest.approx((1.0 - st["damage_t"]) * E * 9.0 * K0, rel=1e-12)


# ===========================================================================
# Bloque 2: daño de compresión (spec §3)
# ===========================================================================
#
# f_c = 30 MPa, G_c = 15 kN/m, β = 1.16: κ_c/3 = f_c/(3E) = 3.33e-4,
# κ_c = 5f_c/(3E) = 1.67e-3 y, con h = 25 mm, κ_u = κ_c + 0.03.

FC, GC = 30.0e6, 15.0e3
KC3 = FC / (3.0 * E)
KCP = 5.0 * FC / (3.0 * E)
BETA = 1.16
KU_H = KCP + 1.5 * GC / (H * FC)


def matc(kind, **kw):
    """Material con aplastamiento (f_c y G_c)."""
    kw.setdefault("f_c", FC)
    kw.setdefault("G_c", GC)
    return mat(kind, **kw)


def feenstra(kappa, h, fc=FC, gc=GC, E_=E):
    """σ_c(κ) de la spec §3, escrita aparte del material."""
    kc3 = fc / (3.0 * E_)
    kc = 5.0 * fc / (3.0 * E_)
    ku = kc + 3.0 * gc / (2.0 * h * fc)
    if kappa <= kc3:
        return E_ * kappa
    if kappa <= kc:
        xi = (kappa - kc3) / (kc - kc3)
        return fc / 3.0 * (1.0 + 4.0 * xi - 2.0 * xi**2)
    if kappa <= ku:
        eta = (kappa - kc) / (ku - kc)
        return fc * (1.0 - eta**2)
    return 0.0


def sigma_eff_3d(m, eps):
    """Tensor 3×3 del esfuerzo efectivo del estado 3D asociado (σ̄_zz = 0 en
    esfuerzo plano; ν·(σ̄_xx + σ̄_yy) en deformación plana)."""
    s = m.elastic_base.C @ eps
    if len(s) == 6:
        return _tensor(s)
    zz = m.nu * (s[0] + s[1]) if m.hypothesis == "plane_strain" else 0.0
    return np.array([[s[0], s[2], 0.0], [s[2], s[1], 0.0], [0.0, 0.0, zz]])


def eps_c_independiente(m, eps, beta=BETA):
    """ε⁻_eq de Drucker-Prager (spec §3) con numpy, sin el código del material."""
    lam, V = np.linalg.eigh(sigma_eff_3d(m, eps))
    neg = V @ np.diag(np.minimum(lam, 0.0)) @ V.T
    p = np.trace(neg) / 3.0
    dev = neg - p * np.eye(3)
    tau = math.sqrt(2.0 / 3.0 * 0.5 * np.sum(dev * dev))
    K = math.sqrt(2.0) * (beta - 1.0) / (2.0 * beta - 1.0)
    return 3.0 * (K * p + tau) / ((math.sqrt(2.0) - K) * E)


def test_compresion_parametros_y_estado_inicial():
    """κ⁻ empieza en f_c/(3E), fin del tramo elástico (como κ⁺ en κ_0)."""
    for m in (matc("plane_stress"), matc("plane_strain"), matc("3d")):
        assert m.has_compression_damage
        assert (m.f_c, m.G_c, m.beta) == (FC, GC, BETA)
        assert m.initial_point_state(H) == {
            "kappa_t": K0, "kappa_c": KC3, "damage_t": 0.0, "damage_c": 0.0, "h": H,
            "kappa_t_prev": K0, "kappa_c_prev": KC3, "damage_error": 0.0}
    assert not mat2().has_compression_damage and not mat3().has_compression_damage


@pytest.mark.parametrize("kind", ["plane_stress", "3d"])
def test_compresion_parametros_invalidos(kind):
    with pytest.raises(ValueError, match="van juntos"):
        mat(kind, f_c=FC)
    with pytest.raises(ValueError, match="van juntos"):
        mat(kind, G_c=GC)
    with pytest.raises(ValueError, match="f_c debe"):
        matc(kind, f_c=-FC)
    with pytest.raises(ValueError, match="G_c debe"):
        matc(kind, G_c=0.0)
    with pytest.raises(ValueError, match="mayor que f_t"):
        matc(kind, f_c=FT)
    for beta in (0.9, float("nan"), "x"):
        with pytest.raises(ValueError, match="beta"):
            matc(kind, beta=beta)


@pytest.mark.parametrize("kind", ["plane_stress", "3d"])
def test_compresion_sin_h_avisa(kind):
    m = matc(kind)
    m.compute_state(uniaxial(kind, -0.9 * KC3), None)               # elástico: admitido
    with pytest.raises(ValueError, match="sin longitud característica"):
        m.compute_state(uniaxial(kind, -2.0 * KC3), None)


@pytest.mark.parametrize("kind", ["plane_stress", "3d"])
@pytest.mark.parametrize("h", [0.01, 0.05, 0.1])
def test_curva_de_compresion_uniaxial(kind, h):
    """σ = −σ_c(κ) en carga monótona (compression_curve): elástica hasta
    f_c/3, pico f_c en κ_c = 5f_c/(3E), rama descendente hasta κ_u y, en el
    tope, la recta (1 − damage_max)·E·κ; las demás componentes, nulas; la
    compresión no hace crecer d⁺."""
    m = matc(kind)
    cap = m.damage_max
    ku = KCP + 1.5 * GC / (h * FC)
    st = m.initial_point_state(h)
    for k in np.linspace(0.0, 1.2 * ku, 3001)[1:]:
        sig, _, st = m.compute_state(uniaxial(kind, -k), st)
        ref = max(feenstra(k, h), (1.0 - cap) * E * k)
        assert abs(sig[0] + ref) <= 1e-8 * FC, k
        assert np.abs(sig[1:]).max() <= 1e-8 * FC
        assert st["damage_t"] == 0.0
    assert st["damage_c"] == cap
    sig = m.compute_state(uniaxial(kind, -KCP), m.initial_point_state(h))[0]
    assert sig[0] == pytest.approx(-FC, rel=1e-12)


@pytest.mark.parametrize("kind", ["plane_stress", "3d"])
@pytest.mark.parametrize("h", [0.01, 0.05, 0.1])
def test_area_de_la_rama_descendente(kind, h):
    """Área de la rama descendente = G_c/h (regularización por banda de grieta
    de la compresión), con el tope descontado en forma cerrada: junto a κ_u,
    donde σ_c < (1 − tope)·E·κ, el material sigue la recta del tope."""
    m = matc(kind)
    c = (1.0 - m.damage_max) * E
    w = 1.5 * GC / (h * FC)
    ku = KCP + w
    # η* donde f_c·(1 − η²) = c·(κ_c + η·w).
    eta_s = (-c * w + math.sqrt((c * w) ** 2 + 4.0 * FC * (FC - c * KCP))) / (2.0 * FC)
    ks = KCP + eta_s * w
    esperado = (GC / h - FC * w * ((1.0 - eta_s) - (1.0 - eta_s**3) / 3.0)
                + 0.5 * c * (ku**2 - ks**2))
    x, wq = np.polynomial.legendre.leggauss(6)        # exacta para cada tramo polinómico
    area = 0.0
    for a, b in ((KCP, ks), (ks, ku)):
        for xi, wi in zip(x, wq):
            k = 0.5 * (a + b) + 0.5 * (b - a) * xi
            sig = m.compute_state(uniaxial(kind, -k), m.initial_point_state(h))[0]
            area -= 0.5 * (b - a) * wi * sig[0]
    assert area == pytest.approx(esperado, rel=1e-10)
    assert esperado == pytest.approx(GC / h, rel=1e-8)      # el tope pesa < 1e-8


@pytest.mark.parametrize("kind", ["plane_stress", "3d"])
@pytest.mark.parametrize("beta", [1.0, 1.16, 1.3])
def test_pico_biaxial_igual(kind, beta):
    """Compresión biaxial igual (σ̄ = −s en x e y): ε⁻_eq = s/(β·E), así que
    el pico es β·f_c en s = β·E·κ_c (biaxial_compression_peak)."""
    m = matc(kind, beta=beta)

    def eps_de(s):
        v = np.zeros(m.STRAIN_DIM)
        v[0] = v[1] = -s
        return np.linalg.solve(m.elastic_base.C, v)

    sig, _, st = m.compute_state(eps_de(beta * E * KCP), m.initial_point_state(H))
    assert sig[0] == pytest.approx(-beta * FC, rel=1e-10)
    assert sig[1] == pytest.approx(-beta * FC, rel=1e-10)
    assert st["kappa_c"] == pytest.approx(KCP, rel=1e-12)
    st = m.initial_point_state(H)
    pico = 0.0
    for s in np.linspace(0.0, 2.0 * beta * E * KCP, 801)[1:]:
        sig, _, st = m.compute_state(eps_de(s), st)
        pico = max(pico, -sig[0])
    assert beta * FC * (1.0 - 1e-4) < pico <= beta * FC * (1.0 + 1e-12)


def test_compresion_hidrostatica_sin_dano():
    """ε⁻_eq·E/|σ̄| = −3K/(√2 − K) < 0 (−0.414 con β = 1.16): la compresión
    hidrostática nunca daña (hydrostatic_compression_no_damage)."""
    from solidum.materials._tension_compression import compression_measure
    K = math.sqrt(2.0) * (BETA - 1.0) / (2.0 * BETA - 1.0)
    eps, _ = compression_measure(-5e6, -5e6, -5e6, 0.0, 0.0, 0.0, E, K)
    assert eps * E / 5e6 == pytest.approx(-3.0 * K / (math.sqrt(2.0) - K), rel=1e-14)
    assert eps * E / 5e6 == pytest.approx(-0.414, abs=5e-4)
    for integration in ("implicit", "impl-ex"):
        m = matc("3d", integration=integration)
        m.set_pseudo_time_ratio(1.0)
        st = m.initial_point_state(H)
        for e in np.linspace(0.0, 100.0 * KCP, 101)[1:]:
            sig, _, st = m.compute_state(-e * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0]), st)
            assert st["damage_c"] == 0.0 and st["damage_t"] == 0.0
            assert st["kappa_c"] == KC3
        np.testing.assert_allclose(sig[:3], -E / (1.0 - 2.0 * NU) * e, rtol=1e-12)


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
@pytest.mark.parametrize("beta", [1.0, 1.16, 1.4])
def test_medida_de_compresion_independiente(kind, beta):
    """κ⁻ tras un paso de carga desde el estado virgen = ε⁻_eq evaluada aparte
    con numpy sobre el estado 3D asociado (traza y J₂ del tensor 3D, σ̄_zz de
    deformación plana incluida, numerador 3)."""
    m = matc(kind, beta=beta)
    rng = np.random.default_rng(59)
    n_checked = 0
    for _ in range(80):
        eps = rng.standard_normal(m.STRAIN_DIM) * 6e-4
        ref = eps_c_independiente(m, eps, beta)
        st = m.compute_state(eps, m.initial_point_state(H))[2]
        if ref > 1.001 * KC3:
            assert st["kappa_c"] == pytest.approx(ref, rel=1e-12)
            n_checked += 1
        elif ref < 0.999 * KC3:
            assert st["kappa_c"] == KC3 and st["damage_c"] == 0.0
    assert n_checked >= 10


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_danos_irreversibles_con_compresion(kind, integration):
    """Compresión, tracción, cierre y recompresión: d⁺ y d⁻ nunca decrecen;
    con todos los esfuerzos principales efectivos ≤ 0 no crece d⁺ y con todos
    ≥ 0 no crece d⁻ (damage_irreversible)."""
    m = matc(kind, integration=integration)
    m.set_pseudo_time_ratio(1.0)
    st = m.initial_point_state(H)
    rng = np.random.default_rng(61)
    n = m.STRAIN_DIM
    base = rng.standard_normal(n) * 0.2
    base[0] = 1.0
    path = np.r_[np.linspace(0, -40, 41)[1:], np.linspace(-40, 6, 47)[1:],
                 np.linspace(6, -80, 87)[1:]]
    dt_prev = dc_prev = 0.0
    for f in path:
        eps = f * K0 * base
        _, _, st = m.compute_state(eps, st)
        assert st["damage_t"] >= dt_prev and st["damage_c"] >= dc_prev
        lam = np.linalg.eigvalsh(sigma_eff_3d(m, eps))
        if lam.max() <= 0.0:
            assert st["damage_t"] == dt_prev
        if lam.min() >= 0.0:
            assert st["damage_c"] == dc_prev
        dt_prev, dc_prev = st["damage_t"], st["damage_c"]
    assert dt_prev > 0.5 and dc_prev > 0.5


def _plano(m, eps):
    """Esfuerzos principales efectivos del plano (2D) o del sólido (3D)."""
    return np.linalg.eigvalsh(_tensor(m.elastic_base.C @ eps))


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
@pytest.mark.parametrize("traccion", ["fija", "en_carga"])
def test_tangente_implicita_con_compresion_por_diferencias_finitas(kind, traccion):
    """Carga activa de d⁻ en las ramas ascendente y descendente: la tangente
    lleva −σ̄⁻⊗∂d⁻/∂ε con (I − Q⁺) y, en deformación plana, el término de
    σ̄_zz (de los dos signos); con d⁺ fijo o también en carga (tangent_fd).
    Se excluyen la arista de Rankine y los esfuerzos principales nulos (§10)."""
    m = matc(kind)
    rng = np.random.default_rng(67)
    n = m.STRAIN_DIM
    checked = 0
    signos_zz = set()
    ramas = {"ascendente": 0, "descendente": 0}
    for _ in range(600):
        eps = rng.standard_normal(n) * 2.0e-3
        ec = eps_c_independiente(m, eps)
        if not (1.5 * KC3 < ec < 0.9 * KU_H):
            continue
        lam = _plano(m, eps)
        escala = np.abs(lam).max()
        if np.abs(lam).min() < 0.02 * escala:
            continue
        if len(lam) > 1 and lam[-1] > 0.0 and lam[-2] > 0.0 and lam[-1] - lam[-2] < 0.05 * lam[-1]:
            continue
        zz = sigma_eff_3d(m, eps)[2, 2] if kind == "plane_strain" else None
        if zz is not None and abs(zz) < 0.02 * escala:
            continue
        st = m.initial_point_state(H)
        st["kappa_c"] = 0.97 * ec                     # carga activa de d⁻ en ±δ
        if traccion == "fija":
            st["kappa_t"] = 3.0 * max(lam[-1], FT) / E
        else:
            if lam[-1] <= 1.2 * FT:
                continue
            st["kappa_t"] = 0.97 * lam[-1] / E        # carga activa de d⁺ también
        st["damage_t"] = exp_damage(st["kappa_t"], H)
        _, C, st1 = m.compute_state(eps, st)
        assert st1["kappa_c"] == pytest.approx(ec, rel=1e-12)
        assert 0.0 < st1["damage_c"] < m.damage_max
        if traccion == "en_carga":
            assert st1["damage_t"] > st["damage_t"]
        C_fd = _fd(m, eps, st, rel=1e-7)
        assert np.abs(C - C_fd).max() <= 1e-6 * np.abs(C_fd).max()
        checked += 1
        ramas["ascendente" if ec < KCP else "descendente"] += 1
        if zz is not None:
            signos_zz.add(np.sign(zz))
    assert checked >= 10
    assert min(ramas.values()) >= 3, ramas
    if kind == "plane_strain":
        assert signos_zz == {-1.0, 1.0}


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
def test_tangente_implicita_en_descarga_con_dano_de_compresion(kind):
    """Implícita con d⁻ > 0 y κ⁻ del estado mayor que ε⁻_eq (descarga): la
    tangente es la de daño fijo, C^d(d⁺, d⁻), sin el término del daño
    (verificación del 2026-11-18, H2)."""
    m = matc(kind)
    ref = matc(kind, integration="impl-ex")             # C^d: IMPL-EX fuera de un paso
    rng = np.random.default_rng(89)
    n = m.STRAIN_DIM
    checked = 0
    for _ in range(200):
        eps = rng.standard_normal(n) * 2.0e-3
        ec = eps_c_independiente(m, eps)
        if not (1.5 * KC3 < ec < 0.8 * KU_H):
            continue
        lam = _plano(m, eps)
        if np.abs(lam).min() < 0.02 * np.abs(lam).max():
            continue
        st = m.initial_point_state(H)
        st["kappa_c"] = 1.2 * ec                        # descarga de d⁻
        st["damage_c"] = 1.0 - feenstra(st["kappa_c"], H) / (E * st["kappa_c"])
        st["kappa_t"] = 3.0 * max(lam[-1], FT) / E      # descarga de d⁺
        st["damage_t"] = exp_damage(st["kappa_t"], H)
        _, C, st1 = m.compute_state(eps, st)
        assert st1["damage_c"] == pytest.approx(st["damage_c"], rel=1e-12)
        assert 0.0 < st1["damage_c"] < m.damage_max
        C_d = ref.compute_state(eps, frozen(ref, st1["damage_t"], st1["damage_c"]))[1]
        np.testing.assert_allclose(C, C_d, rtol=0, atol=1e-12 * np.abs(C_d).max())
        C_fd = _fd(m, eps, st, rel=1e-7)
        assert np.abs(C - C_fd).max() <= 1e-6 * np.abs(C_fd).max()
        checked += 1
    assert checked >= 10


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
def test_tangente_implicita_con_compresion_doble(kind):
    """Esfuerzos principales efectivos negativos dobles con d⁻ en carga, en las
    ramas ascendente y descendente: en el plano, compresión biaxial igual (en
    deformación plana, con σ̄_zz < 0 distinto); en 3D, R·diag(−a, −a, b)·Rᵀ
    con b de los dos signos."""
    m = matc(kind)
    if kind == "3d":
        R, _ = np.linalg.qr(np.random.default_rng(71).standard_normal((3, 3)))
        formas = []
        for b in (0.4, -2.5):
            T = R @ np.diag([-1.0, -1.0, b]) @ R.T
            formas.append(np.array([T[0, 0], T[1, 1], T[2, 2], T[0, 1], T[1, 2], T[0, 2]]))
    else:
        formas = [np.array([-1.0, -1.0, 0.0])]
    ramas = set()
    casos = []
    for forma in formas:
        e1 = np.linalg.solve(m.elastic_base.C, forma)
        for objetivo in (2.0 * KC3, 9.0 * KC3):        # ramas ascendente y descendente
            casos.append(e1 * objetivo / eps_c_independiente(m, e1))   # ε⁻_eq es homogénea de grado 1
    for eps in casos:
        s = m.elastic_base.C @ eps
        st = m.initial_point_state(H)
        st["kappa_t"] = 3.0 * max(np.abs(s).max(), FT) / E
        st["damage_t"] = exp_damage(st["kappa_t"], H)
        ec = eps_c_independiente(m, eps)
        st["kappa_c"] = 0.97 * ec
        _, C, st1 = m.compute_state(eps, st)
        assert 0.0 < st1["damage_c"] < m.damage_max
        C_fd = _fd(m, eps, st, rel=1e-7)
        assert np.abs(C - C_fd).max() <= 1e-6 * np.abs(C_fd).max()
        ramas.add("ascendente" if ec < KCP else "descendente")
    assert ramas == {"ascendente", "descendente"}


@pytest.mark.parametrize("hypothesis", ["plane_stress", "plane_strain"])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_restriccion_exacta_del_3d_con_compresion(hypothesis, integration):
    """Historias con d⁺ y d⁻ activos y σ̄_zz de los dos signos: estado y
    esfuerzo iguales al 3D con el estado asociado (restriction_of_3d)."""
    m2 = matc(hypothesis, integration=integration)
    m3 = matc("3d", integration=integration)
    for m in (m2, m3):
        m.set_pseudo_time_ratio(1.0)
    for base in (np.array([1.0, -0.15, 0.35]), np.array([-0.3, 1.0, -0.6]),
                 np.array([0.7, 0.5, 0.9])):
        s2, s3 = m2.initial_point_state(H), m3.initial_point_state(H)
        for f in (-5.0, -12.0, -25.0, -18.0, 3.0, 6.0, -30.0, -60.0, 8.0, -90.0, 12.0):
            e2 = f * K0 * base
            sig2, _, s2 = m2.compute_state(e2, s2)
            sig3, _, s3 = m3.compute_state(assoc_3d(e2, NU, hypothesis), s3)
            for k in s2:
                tol = dict(rel=0, abs=1e-12) if k == "damage_error" \
                    else dict(rel=1e-12, abs=1e-18)
                assert s2[k] == pytest.approx(s3[k], **tol), k
            scale = np.abs(sig3).max()
            np.testing.assert_allclose(sig2, sig3[[0, 1, 3]], rtol=1e-12, atol=1e-12 * scale)
            np.testing.assert_allclose(sig3[[4, 5]], 0.0, atol=1e-12 * scale)
            szz = m2.out_of_plane_stress(sig2, s2)
            assert szz == pytest.approx(sig3[2], rel=1e-9, abs=1e-12 * scale)
        assert s2["damage_c"] > 0.1 and s2["damage_t"] > 0.1


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_invariancia_de_unidades_con_compresion(kind, integration):
    """Pa–m frente a MPa–mm con todos los parámetros ≠ 1, también f_c, G_c y β
    (units_invariance; Reglas §6)."""
    par_si = dict(E=31.7e9, nu=0.18, f_t=2.9e6, G_f=93.0, f_c=33.5e6, G_c=12.7e3, beta=1.21)
    par_mm = dict(E=31.7e3, nu=0.18, f_t=2.9, G_f=0.093, f_c=33.5, G_c=12.7, beta=1.21)
    h_si, h_mm = 0.037, 37.0

    def build(p):
        if kind == "3d":
            return TensionCompressionDamage3D(**p, integration=integration)
        return TensionCompressionDamage2D(**p, hypothesis=kind, integration=integration)

    a, b = build(par_si), build(par_mm)
    for m in (a, b):
        m.set_pseudo_time_ratio(1.0)
    sa, sb = a.initial_point_state(h_si), b.initial_point_state(h_mm)
    rng = np.random.default_rng(73)
    base = rng.standard_normal(a.STRAIN_DIM) * 0.3
    base[0] = 1.0
    k0 = 2.9e6 / 31.7e9
    Wa = Wb = 0.0
    e_prev = np.zeros(a.STRAIN_DIM)
    siga_prev = sigb_prev = np.zeros(a.STRAIN_DIM)
    for f in np.r_[np.linspace(0, 6, 31)[1:], np.linspace(6, -60, 67)[1:],
                   np.linspace(-60, 9, 35)[1:]]:
        eps = f * k0 * base
        siga, _, sa = a.compute_state(eps, sa)
        sigb, _, sb = b.compute_state(eps, sb)
        np.testing.assert_allclose(sigb, 1e-6 * siga, rtol=1e-12, atol=1e-12 * np.abs(sigb).max())
        for k in ("kappa_t", "kappa_c", "damage_t", "damage_c", "damage_error"):
            assert sb[k] == pytest.approx(sa[k], rel=1e-12, abs=1e-15), k
        Wa += 0.5 * (siga + siga_prev) @ (eps - e_prev)
        Wb += 0.5 * (sigb + sigb_prev) @ (eps - e_prev)
        e_prev, siga_prev, sigb_prev = eps, siga, sigb
    assert sa["damage_c"] > 0.3 and sa["damage_t"] > 0.3
    assert Wb == pytest.approx(1e-6 * Wa, rel=1e-12)


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
def test_kernel_por_lotes_con_compresion(kind, integration):
    """El kernel por lotes reproduce bit a bit σ, C_tan y el estado con d⁻
    activo (batch_kernel_matches)."""
    m = matc(kind, integration=integration)
    n = m.STRAIN_DIM
    names = list(m.STATE_SCHEMA)
    kernel = m.batch_kernel()
    C = m.batch_matrix()
    st = m.initial_point_state(H)
    rng = np.random.default_rng(79)
    for i in range(40):
        m.set_pseudo_time_ratio(None if i % 5 == 4 else 0.5 + 0.1 * i)
        params = m.batch_params()
        f = 1.0 if i % 3 == 0 else -2.0 - 0.2 * i
        eps = (rng.standard_normal(n) * 0.3 + np.r_[1.0, np.zeros(n - 1)]) * 1.0e-3 * f
        sig_ref, C_ref, st_new = m.compute_state(eps, st)
        S_old = np.array([st[k] for k in names])
        S_new, sig, C_out, flag = (np.zeros(len(names)), np.zeros(n), np.zeros((n, n)),
                                   np.zeros(1, np.int8))
        C_t = kernel(eps, S_old, S_new, params, C, sig, C_out, flag)
        np.testing.assert_array_equal(sig, sig_ref)
        np.testing.assert_array_equal(C_t, C_ref)
        np.testing.assert_array_equal(S_new, [st_new[k] for k in names])
        assert flag[0] == 0
        st = st_new
    assert st["damage_c"] > 0.0 and st["damage_t"] > 0.0


@pytest.mark.parametrize("kind", ["plane_stress", "3d"])
def test_kernel_por_lotes_marca_compresion_sin_h(kind):
    m = matc(kind)
    n = m.STRAIN_DIM
    names = list(m.STATE_SCHEMA)
    S_old = np.array([m.initial_state()[k] for k in names])
    flag = np.zeros(1, np.int8)
    m.batch_kernel()(uniaxial(kind, -0.9 * KC3), S_old, np.zeros(len(names)),
                     m.batch_params(), m.batch_matrix(), np.zeros(n), np.zeros((n, n)), flag)
    assert flag[0] == 0
    m.batch_kernel()(uniaxial(kind, -2.0 * KC3), S_old, np.zeros(len(names)),
                     m.batch_params(), m.batch_matrix(), np.zeros(n), np.zeros((n, n)), flag)
    assert flag[0] == 1


def _camino_c(n, integration, kappa_de_t):
    m = matc("plane_stress", integration=integration)
    st = m.initial_point_state(H)
    m.set_pseudo_time_ratio(1.0)
    out = []
    for t in np.linspace(0.0, 1.0, n + 1)[1:]:
        sig, _, st = m.compute_state(uniaxial("plane_stress", -kappa_de_t(t)), st)
        out.append(sig[0])
    return np.array(out)


@pytest.mark.parametrize("forma, orden", [("lineal_desde_el_virgen", 2.0), ("cuadratica", 2.0),
                                           ("descarga_y_recarga", 1.0)])
def test_implex_en_compresion_convergencia_a_la_implicita(forma, orden):
    """Error máximo de σ frente a la implícita con pasos Δ, Δ/2, Δ/4, Δ/8 en
    compresión, hasta la rama descendente (first_order_convergence).

    - Camino lineal desde el estado virgen: segundo orden (razón 4). La
      extrapolación es exacta desde el segundo paso; el único error es el del
      primero, (2f_c/3)·ξ₁², porque d⁻ empieza con pendiente nula
      (σ_c' = E en κ_c/3). En tracción el mismo primer paso da O(Δ).
    - Camino cuadrático: segundo orden, con el error a media historia.
    - Descarga y recarga en la rama descendente: primer orden (razón 2), por
      el quiebre de κ(t) al recargar dentro de un paso.

    Medido el 2026-11-18 (verificación de ese día, A y H5)."""
    if forma == "lineal_desde_el_virgen":
        def kappa(t):
            return KC3 * (1.0 + 39.0 * t)
    elif forma == "cuadratica":
        def kappa(t):
            return KC3 * (1.0 + 39.0 * t * t)
    else:
        def kappa(t):
            return KC3 * np.interp(t, [0.0, 0.4, 0.55, 1.0], [0.0, 20.0, 10.0, 40.0])
    errores = [np.abs(_camino_c(n, "impl-ex", kappa) - _camino_c(n, "implicit", kappa)).max()
               for n in (100, 200, 400, 800)]
    razones = np.array(errores[:-1]) / np.array(errores[1:])
    np.testing.assert_allclose(razones, 2.0 ** orden, rtol=0.15)
    assert kappa(1.0) > KCP                           # el camino pasa el pico


@pytest.mark.parametrize("kind", ["plane_stress", "3d"])
def test_implex_en_compresion_con_pasos_desiguales(kind):
    """Camino lineal en pseudo-tiempo con pasos desiguales y sus razones
    r = Δτ_{n+1}/Δτ_n: la extrapolación de κ⁻ es exacta desde el segundo paso,
    así que σ y d̃⁻ coinciden con la implícita; en el primero d̃⁻ = 0 y el
    indicador vale |d⁻(κ⁻) − d̃⁻| = d⁻(κ⁻_1) (verificación del 2026-11-18, H3)."""
    a = matc(kind, integration="impl-ex")
    b = matc(kind)
    sa, sb = a.initial_point_state(H), b.initial_point_state(H)
    taus = np.cumsum([0.0, 0.03, 0.07, 0.02, 0.11, 0.05, 0.13, 0.04, 0.09, 0.2])
    dtau_prev = None
    for k in range(1, len(taus)):
        dtau = taus[k] - taus[k - 1]
        a.set_pseudo_time_ratio(1.0 if dtau_prev is None else dtau / dtau_prev)
        dtau_prev = dtau
        kap = KC3 * (1.0 + 40.0 * taus[k])
        sig_a, _, sa = a.compute_state(uniaxial(kind, -kap), sa)
        sig_b, _, sb = b.compute_state(uniaxial(kind, -kap), sb)
        d_impl = 1.0 - feenstra(kap, H) / (E * kap)
        if k == 1:
            assert sa["damage_c"] == 0.0
            assert sa["damage_error"] == pytest.approx(d_impl, rel=1e-12)
        else:
            assert sa["damage_c"] == pytest.approx(d_impl, rel=1e-12, abs=1e-15)
            assert sa["damage_error"] < 1e-13
            np.testing.assert_allclose(sig_a, sig_b, rtol=1e-12, atol=1e-12 * FC)
    assert sa["kappa_c"] > KCP                         # el camino pasa el pico


def test_implex_dano_de_compresion_monotono():
    """Carga, descarga parcial y recarga en compresión con una razón que hace
    pasarse a la extrapolación: d̃⁻ nunca decrece (damage_monotone)."""
    m = matc("plane_stress", integration="impl-ex")
    st = m.initial_point_state(H)
    camino = np.r_[np.linspace(0.0, 30.0, 41)[1:], np.linspace(30.0, 12.0, 11)[1:],
                   np.linspace(12.0, 60.0, 31)[1:]] * KC3
    d_prev = 0.0
    for k, e in enumerate(camino):
        m.set_pseudo_time_ratio(2.0 if k % 7 == 0 else 1.0)
        _, _, st = m.compute_state(uniaxial("plane_stress", -e), st)
        assert st["damage_c"] >= d_prev
        assert st["damage_error"] >= 0.0
        d_prev = st["damage_c"]
    assert d_prev > 0.9


@pytest.mark.parametrize("kind", ["plane_stress", "plane_strain", "3d"])
def test_implex_en_compresion_tangente_y_estado_congelado(kind):
    """Dentro de un paso con d̃⁻ > 0 la tangente es C^d(d̃⁺, d̃⁻), igual a
    diferencias finitas (fixed_damage_tangent); repetir el intento da el
    mismo estado (bisection_idempotent) y fuera de un paso el estado queda
    congelado (outside_step_frozen)."""
    m = matc(kind, integration="impl-ex")
    n = m.STRAIN_DIM
    base = -np.r_[1.0, 0.3 * np.ones(n - 1)] * np.r_[1.0, -0.5, 0.7, 0.2, -0.4, 0.3][:n]
    m.set_pseudo_time_ratio(1.0)
    st = m.initial_point_state(H)
    for f in np.linspace(0.0, 12.0, 25)[1:]:
        _, _, st = m.compute_state(f * KC3 * base, st)
    assert st["damage_c"] > 0.0
    rng = np.random.default_rng(83)
    for _ in range(8):
        eps = 12.5 * KC3 * base + rng.standard_normal(n) * 2e-4
        sig, C, st1 = m.compute_state(eps, st)
        assert np.abs(C - _fd(m, eps, st)).max() <= 1e-6 * np.abs(C).max()
        again = m.compute_state(eps, st)
        np.testing.assert_array_equal(again[0], sig)
        assert again[2] == st1
    m.set_pseudo_time_ratio(None)
    sig, _, nuevo = m.compute_state(20.0 * KC3 * base, st)
    assert nuevo == st
