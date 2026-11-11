"""``IsotropicDamageCrackBand2D/3D``: punto material (spec §1–§5, aceptación
``material_point``) y siembra de la longitud característica en los sólidos
(aceptación ``elements``).

Concreto usual en SI: E = 30 GPa, f_t = 3 MPa, G_f = 100 N/m
(κ_0 = 1e-4, h_max = 2·E·G_f/f_t² = 0.667 m).
"""
import math

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.constants import CRACK_BAND_DAMAGE_MAX
from solidum.core.node import Node
from solidum.elements.solid_2d import Quad4, Quad8, Quad9, Tri3, Tri6
from solidum.elements.solid_3d import Hex8, Hex20, Hex27, Tet4, Tet10
from solidum.materials.damage_2d import IsotropicDamage2D
from solidum.materials.damage_crack_band_2d import IsotropicDamageCrackBand2D
from solidum.materials.damage_crack_band_3d import IsotropicDamageCrackBand3D
from solidum.materials.elastic_2d import Elastic2D

E, NU, FT, GF = 30.0e9, 0.2, 3.0e6, 100.0
K0 = FT / E
MEASURES = ("rankine", "energy", "mazars", "modified_von_mises")


def _kw(measure):
    return {"k": 10.0} if measure == "modified_von_mises" else {}


def mat2(measure="rankine", hypothesis="plane_stress", **kw):
    return IsotropicDamageCrackBand2D(E, NU, FT, GF, hypothesis=hypothesis,
                                      equivalent_strain=measure, **_kw(measure), **kw)


def mat3(measure="rankine", **kw):
    return IsotropicDamageCrackBand3D(E, NU, FT, GF, equivalent_strain=measure,
                                      **_kw(measure), **kw)


def uniaxial_2d(eps):
    """Tracción uniaxial de esfuerzo en tensión plana."""
    return np.array([eps, -NU * eps, 0.0])


def dissipated(kappa, h):
    """D(κ) de la spec §2: energía disipada por unidad de volumen."""
    a = 1.0 / (GF / (h * FT) - K0 / 2.0)
    x = math.exp(-a * (kappa - K0))
    return 0.5 * FT * K0 + FT / a * (1.0 - x) - 0.5 * FT * x * kappa


# ---------------------------------------------------------------------------
# Parámetros, α(h) y h_max
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("h", [0.005, 0.025, 0.100, 0.400])
def test_alpha_forma_cerrada(h):
    for m in (mat2(), mat3()):
        expected = 1.0 / (GF / (h * FT) - K0 / 2.0)
        assert m.alpha_for(h) == pytest.approx(expected, rel=1e-14)
        assert m.initial_point_state(h) == {"kappa": K0, "damage": 0.0, "h": h,
                                            "kappa_prev": K0, "damage_error": 0.0}


@pytest.mark.parametrize("h", [2.0 * E * GF / FT**2, 1.01 * 2.0 * E * GF / FT**2])
def test_h_max_rechazada(h):
    for m in (mat2(), mat3()):
        with pytest.raises(ValueError, match="h_max"):
            m.initial_point_state(h)


@pytest.mark.parametrize("h", [0.0, -0.01, float("nan"), float("inf")])
def test_h_no_valida(h):
    with pytest.raises(ValueError, match="longitud característica"):
        mat2().initial_point_state(h)


def test_parametros_invalidos():
    with pytest.raises(ValueError, match="f_t"):
        IsotropicDamageCrackBand2D(E, NU, 0.0, GF)
    with pytest.raises(ValueError, match="G_f"):
        IsotropicDamageCrackBand3D(E, NU, FT, -1.0)
    with pytest.raises(ValueError, match="damage_max"):
        IsotropicDamageCrackBand2D(E, NU, FT, GF, damage_max=1.0)
    with pytest.raises(ValueError, match="equivalent_strain"):
        IsotropicDamageCrackBand2D(E, NU, FT, GF, equivalent_strain="tresca")


def test_el_padre_no_admite_rankine():
    """Los padres validados conservan sus tres medidas."""
    with pytest.raises(ValueError, match="rankine"):
        IsotropicDamage2D(E, NU, kappa_0=K0, alpha=100.0, equivalent_strain="rankine")


def test_rankine_por_omision_y_tope_propio():
    m = IsotropicDamageCrackBand2D(E, NU, FT, GF)
    assert m.equivalent_strain == "rankine"
    assert m.damage_max == CRACK_BAND_DAMAGE_MAX == 1.0 - 1.0e-6
    assert m.uses_characteristic_length
    assert m.IS_SYMMETRIC is False
    assert mat2("energy").IS_SYMMETRIC is True


def test_dano_sin_h_avisa():
    m = mat2()
    # Régimen elástico sin h: admitido (barrido de contratos, primer paso).
    m.compute_state(uniaxial_2d(0.5 * K0), None)
    with pytest.raises(ValueError, match="sin longitud característica"):
        m.compute_state(uniaxial_2d(2.0 * K0), None)


# ---------------------------------------------------------------------------
# Pico, inicio y energía disipada
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("measure", MEASURES)
def test_pico_e_inicio(measure):
    m = mat2(measure)
    st = m.initial_point_state(0.025)
    s, _, st1 = m.compute_state(uniaxial_2d(K0), st)
    assert st1["damage"] == 0.0
    assert s[0] == pytest.approx(FT, rel=1e-12)
    peak = 0.0
    for eps in np.linspace(0.0, 5.0 * K0, 501):
        s, _, st = m.compute_state(uniaxial_2d(eps), st)
        peak = max(peak, s[0])
    assert peak == pytest.approx(FT, rel=1e-12)


@pytest.mark.parametrize("h", [0.005, 0.025, 0.100, 0.400])
def test_energia_disipada_por_volumen(h):
    """∫σ dε − ½σε = D(κ) en estados intermedios y ≈ G_f/h al final."""
    m = mat2()
    st = m.initial_point_state(h)
    a = m.alpha_for(h)
    eps_end = K0 + 25.0 / a
    path = np.concatenate([np.linspace(0.0, K0, 51), np.linspace(K0, eps_end, 40001)[1:]])
    work, s_prev, e_prev, checked = 0.0, 0.0, 0.0, 0
    for i, eps in enumerate(path):
        s, _, st = m.compute_state(uniaxial_2d(eps), st)
        work += 0.5 * (s[0] + s_prev) * (eps - e_prev)
        s_prev, e_prev = s[0], eps
        # D(κ) cerrada vale mientras el daño no llega al tope (§3).
        if i % 4000 == 0 and 0.0 < st["damage"] < m.damage_max:
            checked += 1
            D_num = work - 0.5 * s[0] * eps
            assert D_num == pytest.approx(dissipated(eps, h), rel=1e-6)
    assert checked >= 3
    assert work - 0.5 * s_prev * e_prev > 0.997 * GF / h


def test_tope_del_padre_perderia_energia():
    """Con el DAMAGE_MAX = 0.999 del padre y h = 5 mm se pierde ≈ 27 % de G_f
    (spec §3); con el tope propio, < 0.25 %."""
    h = 0.005
    out = {}
    for cap in (0.999, CRACK_BAND_DAMAGE_MAX):
        m = mat2(damage_max=cap)
        st = m.initial_point_state(h)
        eps = K0
        while True:
            eps *= 1.0005
            s, _, st = m.compute_state(uniaxial_2d(eps), st)
            if st["damage"] >= cap:
                break
        out[cap] = dissipated(st["kappa"], h) * h / GF
    assert out[0.999] == pytest.approx(0.732, abs=0.005)
    assert out[CRACK_BAND_DAMAGE_MAX] > 0.9975


# ---------------------------------------------------------------------------
# Rankine
# ---------------------------------------------------------------------------

def _plane_strain_for_stress(sx, sy, sxy, hypothesis):
    C = Elastic2D(E, NU, hypothesis=hypothesis).C
    return np.linalg.solve(C, np.array([sx, sy, sxy]))


@pytest.mark.parametrize("hypothesis", ["plane_stress", "plane_strain"])
def test_rankine_valores(hypothesis):
    m = mat2("rankine", hypothesis)
    h = 0.025
    # Tracción biaxial igual (efectiva) y cortante puro: el daño empieza en f_t.
    for stress in ((FT, FT, 0.0), (0.0, 0.0, FT)):
        if hypothesis == "plane_strain" and stress[0] > 0:
            continue    # σ̄_zz = ν·(σ̄_xx + σ̄_yy) no supera a σ̄_xx con ν < 0.5
        eps = _plane_strain_for_stress(*stress, hypothesis)
        _, _, st = m.compute_state(eps, m.initial_point_state(h))
        assert st["kappa"] == pytest.approx(K0, rel=1e-12)
        assert st["damage"] == 0.0
        _, _, st = m.compute_state(1.001 * eps, m.initial_point_state(h))
        assert st["damage"] > 0.0
    # Compresión uniaxial y biaxial: nunca daña.
    for stress in ((-100.0 * FT, 0.0, 0.0), (-100.0 * FT, -100.0 * FT, 0.0)):
        eps = _plane_strain_for_stress(*stress, hypothesis)
        _, _, st = m.compute_state(eps, m.initial_point_state(h))
        assert st["damage"] == 0.0


def test_rankine_3d_valores():
    m = mat3("rankine")
    C = m.elastic_base.C
    h = 0.025
    for stress in ([FT, FT, FT, 0, 0, 0], [0, 0, 0, FT, 0, 0], [0, 0, 0, 0, 0, FT]):
        eps = np.linalg.solve(C, np.array(stress, float))
        _, _, st = m.compute_state(eps, m.initial_point_state(h))
        assert st["kappa"] == pytest.approx(K0, rel=1e-12) and st["damage"] == 0.0
        _, _, st = m.compute_state(1.001 * eps, m.initial_point_state(h))
        assert st["damage"] > 0.0
    for stress in ([-100 * FT, 0, 0, 0, 0, 0], [-100 * FT, -100 * FT, 0, 0, 0, 0],
                   [-100 * FT] * 3 + [0, 0, 0]):
        eps = np.linalg.solve(C, np.array(stress, float))
        _, _, st = m.compute_state(eps, m.initial_point_state(h))
        assert st["damage"] == 0.0


# ---------------------------------------------------------------------------
# Tangente consistente por diferencias finitas
# ---------------------------------------------------------------------------

def _fd_check(m, eps, st):
    _, C, _ = m.compute_state(eps, st)
    n = len(eps)
    C_fd = np.zeros((n, n))
    for j in range(n):
        d = 1e-6 * max(abs(eps[j]), np.abs(eps).max())
        ep, em = eps.copy(), eps.copy()
        ep[j] += d
        em[j] -= d
        C_fd[:, j] = (m.compute_state(ep, st)[0] - m.compute_state(em, st)[0]) / (2 * d)
    assert np.abs(C - C_fd).max() <= 1e-6 * np.abs(C_fd).max()


@pytest.mark.parametrize("measure", MEASURES)
@pytest.mark.parametrize("hypothesis", ["plane_stress", "plane_strain"])
def test_tangente_consistente_2d(measure, hypothesis):
    m = mat2(measure, hypothesis)
    st = m.initial_point_state(0.025)
    for eps in (np.array([3.0e-4, -0.4e-4, 1.1e-4]), np.array([2.0e-4, 1.5e-4, -0.7e-4])):
        _fd_check(m, eps, st)


@pytest.mark.parametrize("measure", MEASURES)
def test_tangente_consistente_3d(measure):
    m = mat3(measure)
    st = m.initial_point_state(0.025)
    for eps in (np.array([3.0e-4, -0.4e-4, 0.3e-4, 1.1e-4, -0.5e-4, 0.8e-4]),
                np.array([1.0e-4, 2.2e-4, -0.6e-4, -0.9e-4, 0.4e-4, 0.2e-4])):
        _fd_check(m, eps, st)


# ---------------------------------------------------------------------------
# Restricción exacta del 3D y camino por lotes
# ---------------------------------------------------------------------------

def _assoc_3d(eps2, nu, hypothesis):
    c_zz = -nu / (1.0 - nu) if hypothesis == "plane_stress" else 0.0
    return np.array([eps2[0], eps2[1], c_zz * (eps2[0] + eps2[1]), eps2[2], 0.0, 0.0])


@pytest.mark.parametrize("measure", MEASURES)
@pytest.mark.parametrize("hypothesis", ["plane_stress", "plane_strain"])
def test_restriccion_exacta_del_3d(measure, hypothesis):
    m2, m3 = mat2(measure, hypothesis), mat3(measure)
    h = 0.025
    s2, s3 = m2.initial_point_state(h), m3.initial_point_state(h)
    base = np.array([1.0, -0.15, 0.35])
    for f in (0.5, 1.5, 3.0, 2.0, 1.0, 4.0, 6.0):       # carga, descarga y recarga
        e2 = f * K0 * base
        sig2, _, s2 = m2.compute_state(e2, s2)
        sig3, _, s3 = m3.compute_state(_assoc_3d(e2, NU, hypothesis), s3)
        assert s2["kappa"] == pytest.approx(s3["kappa"], rel=1e-12, abs=1e-18)
        assert s2["damage"] == pytest.approx(s3["damage"], rel=1e-12, abs=1e-15)
        np.testing.assert_allclose(sig2, sig3[[0, 1, 3]], rtol=1e-10, atol=1e-6)


@pytest.mark.parametrize("integration", ["implicit", "impl-ex"])
@pytest.mark.parametrize("dim", [2, 3])
def test_kernel_por_lotes_reproduce_compute_state(dim, integration):
    m = mat2(integration=integration) if dim == 2 else mat3(integration=integration)
    n = m.STRAIN_DIM
    names = list(m.STATE_SCHEMA)
    kernel = m.batch_kernel()
    C = m.batch_matrix()
    st = m.initial_point_state(0.025)
    rng = np.random.default_rng(3)
    for i in range(20):
        # IMPL-EX: razón variable por paso y, cada cinco, «fuera de un paso».
        m.set_pseudo_time_ratio(None if i % 5 == 4 else 0.5 + 0.1 * i)
        params = m.batch_params()
        eps = (rng.standard_normal(n) + np.r_[2.0, np.zeros(n - 1)]) * 1.5e-4
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


# ---------------------------------------------------------------------------
# Siembra de la longitud característica en los elementos sólidos
# ---------------------------------------------------------------------------

REF_MEASURE = {"quad": 4.0, "tri": 0.5, "hex": 8.0, "tet": 1.0 / 6.0}
SOLIDS = [
    (Tri3, 2.0 ** 0.5), (Quad4, 1.0), (Tri6, 0.5), (Quad8, 0.5), (Quad9, 0.5),
    (Hex8, 1.0), (Tet4, 6.0 ** (1.0 / 3.0)), (Hex20, 0.5), (Hex27, 0.5), (Tet10, 0.5),
]
A2 = np.array([[0.031, 0.007], [-0.004, 0.022]])                  # afín, distorsionada
A3 = np.array([[0.031, 0.007, 0.002], [-0.004, 0.022, 0.003], [0.001, -0.002, 0.027]])


def _element(cls, A, material):
    ref = np.array(cls.REFERENCE_NODES, dtype=float)
    X = ref @ A.T + 0.1
    nodes = [Node(i + 1, list(x)) for i, x in enumerate(X)]
    return cls(1, nodes, material)


@pytest.mark.parametrize("cls,factor", SOLIDS)
@pytest.mark.parametrize("shape", ["regular", "distorsionada"])
def test_h_sembrada_en_los_solidos(cls, factor, shape):
    dim = len(cls.REFERENCE_NODES[0])
    if shape == "regular":
        A = 0.025 * np.eye(dim)
    else:
        A = A2 if dim == 2 else A3
    material = mat2() if dim == 2 else mat3()
    el = _element(cls, A, material)
    assert cls.CHARACTERISTIC_LENGTH_FACTOR == pytest.approx(factor, rel=1e-15)
    measure = abs(np.linalg.det(A)) * REF_MEASURE[cls.REFERENCE_CELL]
    expected = factor * measure ** (1.0 / dim)
    hs = [v["h"] for v in el.state.vars]
    assert len(hs) == el.state.num_ip
    np.testing.assert_allclose(hs, expected, rtol=1e-12)


def test_h_sembrada_en_un_cuadrilatero_no_afin():
    """Quad4 bilineal: el área es la del polígono (fórmula del área de Gauss)."""
    X = np.array([[0.0, 0.0], [0.04, 0.005], [0.035, 0.03], [-0.003, 0.025]])
    el = Quad4(1, [Node(i + 1, list(x)) for i, x in enumerate(X)], mat2())
    x, y = X[:, 0], X[:, 1]
    area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    np.testing.assert_allclose([v["h"] for v in el.state.vars], math.sqrt(area), rtol=1e-12)


def test_sin_material_de_banda_no_se_siembra():
    X = 0.025 * np.array(Quad4.REFERENCE_NODES)
    el = Quad4(1, [Node(i + 1, list(x)) for i, x in enumerate(X)], Elastic2D(E, NU))
    assert not hasattr(el, "characteristic_lengths")


class _QuadSinH(Quad4):
    """Elemento que no proporciona h (sólo para la prueba; no registrado)."""
    CHARACTERISTIC_LENGTH_FACTOR = None


def test_elemento_sin_h_rechaza_el_material():
    X = 0.025 * np.array(Quad4.REFERENCE_NODES)
    with pytest.raises(ValueError, match="longitud característica"):
        _QuadSinH(1, [Node(i + 1, list(x)) for i, x in enumerate(X)], mat2())


def test_h_sobre_el_limite_se_rechaza_al_construir():
    """Un elemento de 1 m con un concreto de h_max = 0.667 m se rechaza al construirlo."""
    X = 0.5 * np.array(Quad4.REFERENCE_NODES)     # 1 m × 1 m
    with pytest.raises(ValueError, match="Quad4.*h_max"):
        Quad4(7, [Node(i + 1, list(x)) for i, x in enumerate(X)], mat2())
