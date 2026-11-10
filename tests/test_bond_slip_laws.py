"""Leyes de adherencia-deslizamiento ``BondSlipModelCode`` y
``BondSlipMultilinear`` (specs del mismo nombre): casos de aceptación.

Valores de Eligehausen, Popov y Bertero (1983, §5.3.1): τ1 = 13.5 MPa,
s1 = 1, s2 = 3, s3 = 10.5 mm, α = 0.4, τ_f = τ3 = 5 MPa, K = 180 N/mm³.
"""
import logging

import numpy as np
import pytest

from solidum import BondSlipModelCode, BondSlipMultilinear
from solidum.registry import InterfaceMaterialRegistry

MPA, MM = 1.0e6, 1.0e-3
T1, T3, S1, S2, S3, ALPHA, K0 = 13.5 * MPA, 5.0 * MPA, 1 * MM, 3 * MM, 10.5 * MM, 0.4, 1.8e11
KN = 3.0e12


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.WARNING)
    yield
    logging.disable(logging.NOTSET)


def epb():
    return BondSlipModelCode(tau_max=T1, s1=S1, s2=S2, s3=S3, alpha=ALPHA, tau_f=T3, k0=K0, k_n=KN)


def tau_at(law, s, state=None, dim=2):
    jump = np.zeros(dim)
    jump[-1] = s
    t, D, st = law.compute_traction(jump, state)
    return t[-1], D[-1, -1], st


def go(law, state, s):
    """Confirma el estado en ``s`` desde ``state`` (un paso)."""
    return tau_at(law, s, state)[2]


# --- envolvente --------------------------------------------------------------

def test_envolvente_y_simetria():
    law = epb()
    s0 = S1 * (K0 * S1 / T1) ** (1 / (ALPHA - 1))
    assert law.s0 == pytest.approx(s0, rel=1e-14)
    casos = {0.5 * s0: K0 * 0.5 * s0, s0: K0 * s0, 0.4 * MM: T1 * 0.4 ** ALPHA,
             S1: T1, 2 * MM: T1, 6.75 * MM: T1 - (T1 - T3) * 3.75 / 7.5, S3: T3, 20 * MM: T3}
    for s, esperado in casos.items():
        assert tau_at(law, s)[0] == pytest.approx(esperado, rel=1e-12)
        assert tau_at(law, -s)[0] == pytest.approx(-esperado, rel=1e-12)
    # continuidad en s0: la recta k0·a corta la rama potencial
    assert K0 * s0 == pytest.approx(T1 * (s0 / S1) ** ALPHA, rel=1e-12)


def test_ciclo_figura_5_3b():
    """Recorrido de la Fig. 5.3b de Eligehausen et al. (1983) sin degradación."""
    law = epb()
    st = go(law, None, 1.0 * MM)
    t = lambda s: tau_at(law, s * MM, st)[0] / MPA                     # noqa: E731
    assert t(0.95) == pytest.approx(13.5 - 180 * 0.05, rel=1e-12)       # descarga rígida
    assert t(1 - 18.5 / 180) == pytest.approx(-5.0, rel=1e-12)          # llega a −τ_f
    assert t(0.5) == pytest.approx(-5.0, rel=1e-12)                     # fricción
    assert t(-0.05) == pytest.approx(-5.0, rel=1e-12)                   # τ_env(0.05) < τ_f
    s_i = (5 / 13.5) ** 2.5                                             # corte GHI
    assert t(-s_i) == pytest.approx(-5.0, rel=1e-10)
    assert t(-0.5) == pytest.approx(-13.5 * 0.5 ** 0.4, rel=1e-12)      # envolvente opuesta
    st = go(law, st, -1.0 * MM)
    t = lambda s: tau_at(law, s * MM, st)[0] / MPA                     # noqa: E731
    assert t(-0.95) == pytest.approx(-13.5 + 9.0, rel=1e-12)
    assert t(0.0) == pytest.approx(5.0, rel=1e-12)                      # fricción +τ_f
    s_n = 1 - 8.5 / 180                                                 # punto N
    assert t(s_n - 0.01) == pytest.approx(5.0, rel=1e-12)
    assert t(0.97) == pytest.approx(13.5 + 180 * (0.97 - 1.0), rel=1e-12)   # recarga
    assert t(1.5) == pytest.approx(13.5, rel=1e-12)


def test_tangente_por_diferencias_finitas():
    law = epb()
    h = 1e-9
    estados = {
        "virgen": (None, [0.3 * law.s0, 0.4 * MM, 2 * MM, 6 * MM, 15 * MM, -0.4 * MM, -6 * MM]),
        "tras +1 mm": (go(law, None, 1 * MM), [0.97 * MM, 0.5 * MM, -0.5 * MM, -6 * MM]),
        "tras +1/−1 mm": (go(law, go(law, None, 1 * MM), -1 * MM),
                          [-0.97 * MM, 0.3 * MM, 0.97 * MM, 2 * MM]),
    }
    for nombre, (st, puntos) in estados.items():
        for s in puntos:
            for dim in (2, 3):
                _, dtau, _ = tau_at(law, s, st, dim)
                fd = (tau_at(law, s + h, st, dim)[0] - tau_at(law, s - h, st, dim)[0]) / (2 * h)
                assert dtau == pytest.approx(fd, rel=1e-6, abs=1e-3 * K0 * 1e-6), (nombre, s)


def test_tangente_inicial_es_k0():
    """En estado virgen y s = 0 las dos cotas y la fricción valen cero: la
    respuesta es la elástica, no la de fricción (si no, la interfaz no tendría
    rigidez al empezar)."""
    for law in (epb(), BondSlipMultilinear(points=[[1 * MM, 10 * MPA]], k_n=KN)):
        tau, dtau, _ = tau_at(law, 0.0)
        assert tau == 0.0 and dtau == pytest.approx(law.k0, rel=1e-14)


def test_disipacion_no_negativa():
    rng = np.random.default_rng(11)
    for law in (epb(), BondSlipMultilinear(points=[[0.1 * MM, 8 * MPA], [1 * MM, 12 * MPA],
                                                  [4 * MM, 6 * MPA]], k_n=KN)):
        for _ in range(1000 // 2):
            puntos = np.concatenate([[0.0], rng.uniform(-6, 6, 6) * MM])
            st, s_prev, tau_prev, W = None, 0.0, 0.0, 0.0
            for a, b in zip(puntos[:-1], puntos[1:]):
                for s in np.linspace(a, b, 40)[1:]:
                    tau, _, st = tau_at(law, s, st)
                    W += 0.5 * (tau + tau_prev) * (s - s_prev)
                    s_prev, tau_prev = s, tau
                    assert W - tau ** 2 / (2 * law.k0) >= -1e-9


def test_normal_desacoplada():
    law = epb()
    st = go(law, None, 1 * MM)
    for dim in (2, 3):
        jump = np.full(dim, 2.0e-5)
        jump[-1] = 0.4 * MM
        t, D, _ = law.compute_traction(jump, st)
        np.testing.assert_allclose(t[:-1], KN * jump[:-1], rtol=1e-14)
        assert np.allclose(D - np.diag(np.diag(D)), 0.0)
        np.testing.assert_allclose(np.diag(D)[:-1], KN, rtol=1e-14)


# --- Model Code 2010 -----------------------------------------------------------

def test_model_code_2010():
    fcm, c, db, Ec = 30 * MPA, 10 * MM, 16 * MM, 30e9
    f = 30.0
    po_good = 2.5 * np.sqrt(f) * MPA
    po_other = 1.25 * np.sqrt(f) * MPA
    casos = {
        ("good", "pullout", False): (po_good, 1.0 * MM, 2.0 * MM, c, 0.4 * po_good),
        ("other", "pullout", False): (po_other, 1.8 * MM, 3.6 * MM, c, 0.4 * po_other),
    }
    for (bond, stir), coef in {("good", False): 7.0, ("good", True): 8.0,
                               ("other", False): 5.0, ("other", True): 5.5}.items():
        tau = coef * (f / 25) ** 0.25 * MPA
        tpo, s1po = (po_good, 1.0 * MM) if bond == "good" else (po_other, 1.8 * MM)
        s1 = s1po * (tau / tpo) ** 2.5
        casos[(bond, "splitting", stir)] = (tau, s1, s1, 0.5 * c if stir else 1.2 * s1,
                                           0.4 * tau if stir else 0.0)
    for (bond, failure, stir), (tau, s1, s2, s3, tf) in casos.items():
        law = BondSlipModelCode.model_code_2010(fcm, bond=bond, failure=failure, stirrups=stir,
                                                clear_rib_spacing=c, diameter=db, E_c=Ec)
        got = (law.tau_max, law.s1, law.s2, law.s3, law.tau_f)
        np.testing.assert_allclose(got, (tau, s1, s2, s3, tf), rtol=1e-12, atol=0.0)
        assert law.alpha == 0.4
        assert law.k0 == pytest.approx(1.8e11, rel=1e-14)
        assert law.k_n == pytest.approx(2 * Ec / db, rel=1e-14)
    assert BondSlipModelCode.model_code_2010(fcm, clear_rib_spacing=c, k_n=5e12).k_n == 5e12


def test_rechazos():
    base = dict(tau_max=T1, s1=S1, s2=S2, s3=S3, alpha=ALPHA, tau_f=T3, k0=K0, k_n=KN)
    for cambio in ({"k0": T1 / S1}, {"s2": 0.5 * S1}, {"s3": 2 * MM}, {"alpha": 1.0},
                   {"alpha": 0.0}, {"tau_f": 2 * T1}, {"k_n": 0.0}, {"tau_max": -1.0}):
        with pytest.raises(ValueError):
            BondSlipModelCode(**{**base, **cambio})
    for kw in ({"points": [[0.0, 1.0]]}, {"points": [[1e-3, 1.0], [0.5e-3, 2.0]]},
               {"points": [[1e-3, -1.0]]}, {"points": [[1e-3, 1e6]], "tau_f": 2e6},
               {"points": []}):
        with pytest.raises(ValueError):
            BondSlipMultilinear(k_n=KN, **kw)
    with pytest.raises(ValueError, match="E_c"):
        BondSlipModelCode.model_code_2010(30 * MPA, clear_rib_spacing=0.01)
    with pytest.raises(ValueError, match="clear_rib_spacing"):
        BondSlipModelCode.model_code_2010(30 * MPA, k_n=KN)


def test_registro_y_geometria():
    for name in ("BondSlipModelCode", "BondSlipMultilinear"):
        cls = InterfaceMaterialRegistry.get(name)
        assert cls.GEOMETRY == "line" and cls.IS_SYMMETRIC
    law = epb()
    assert law.accepts("line", 2) and law.accepts("line", 3)
    assert not law.accepts("surface", 2)


# --- multilineal ----------------------------------------------------------------

def test_multilineal_por_puntos():
    law = BondSlipMultilinear(points=[[0.1 * MM, 6 * MPA], [1 * MM, 12 * MPA],
                                      [4 * MM, 3 * MPA]], k_n=KN)
    assert law.k0 == pytest.approx(6 * MPA / (0.1 * MM), rel=1e-14)
    assert law.tau_f == 3 * MPA
    for s, esperado in {0.05 * MM: 3 * MPA, 0.55 * MM: 9 * MPA, 2.5 * MM: 7.5 * MPA,
                        8 * MM: 3 * MPA}.items():
        assert tau_at(law, s)[0] == pytest.approx(esperado, rel=1e-12)
        assert tau_at(law, -s)[0] == pytest.approx(-esperado, rel=1e-12)


def test_multilineal_igual_a_model_code_con_sus_puntos():
    mc = epb()
    s = np.concatenate([np.linspace(mc.s0, S1, 400), [S2, S3]])
    pts = [[x, mc.envelope(x)[0]] for x in s]
    ml = BondSlipMultilinear(points=pts, tau_f=T3, k_n=KN)
    assert ml.k0 == pytest.approx(mc.k0, rel=1e-12)
    hist = np.concatenate([np.linspace(0, 1, 50), np.linspace(1, -1, 100),
                           np.linspace(-1, 1.5, 125)]) * MM
    st_mc = st_ml = None
    for x in hist[1:]:
        a, _, st_mc = tau_at(mc, x, st_mc)
        b, _, st_ml = tau_at(ml, x, st_ml)
        assert b == pytest.approx(a, rel=1e-3, abs=1e-3 * T1)
