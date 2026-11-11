"""Aceptación de ``CohesiveDamageMixedMode`` (spec docs/specs/CohesiveDamageMixedMode.md).

Ley cohesiva de modo mixto con apertura efectiva δ_eff = √(⟨δ_n⟩² + β²‖δ_s‖²)
(Camacho y Ortiz 1996), energía independiente del modo y ablandamiento lineal,
exponencial o bilineal (fib Model Code 2010). Los casos siguen el bloque
``acceptance`` de la spec.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from solidum.materials.interface import CohesiveDamageIsotropic, CohesiveDamageMixedMode
from solidum.registry import InterfaceMaterialRegistry

SIGMA_T0 = 3.0e6        # Pa
G_F = 100.0             # N/m
K_E = 1.0e15            # Pa/m, penalización real: κ_0/w_c ≈ 4.5e-5
BETA = 1.5
SHAPES = ("linear", "exponential", "bilinear")
W_REF = G_F / SIGMA_T0  # escala de apertura


def law(shape="linear", beta=BETA, K_e=K_E, **kw):
    return CohesiveDamageMixedMode(sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_e, beta=beta,
                                   softening=shape, **kw)


def end_opening(m):
    """Apertura efectiva desde la que la tracción es nula (o despreciable)."""
    if m.softening == "exponential":
        return 30.0 * W_REF
    return m._curve.kappa_end


def direction(theta_deg, beta=BETA, dim=2, phi_deg=40.0):
    """Salto con δ_eff = 1 en la dirección θ del plano (⟨δ_n⟩, β‖δ_s‖)."""
    th = math.radians(theta_deg)
    s = math.sin(th) / beta
    if dim == 2:
        return np.array([math.cos(th), s])
    ph = math.radians(phi_deg)
    return np.array([math.cos(th), s * math.cos(ph), s * math.sin(ph)])


def weights(dim, beta=BETA):
    w = np.full(dim, beta ** 2)
    w[0] = 1.0
    return np.diag(w)


def walk(m, path, state=None):
    """Recorre una trayectoria confirmando el estado en cada punto; devuelve
    las tracciones y el estado final."""
    tractions = []
    for jump in path:
        t, _, state = m.compute_traction(jump, state)
        tractions.append(t)
    return np.array(tractions), state


# --- verificación -------------------------------------------------------------------

@pytest.mark.parametrize("dim", [2, 3])
def test_respuesta_elastica_bajo_umbral(dim):
    m = law()
    W = weights(dim)
    for theta in (0.0, 30.0, 60.0, 90.0):
        jump = 0.5 * m.kappa_0 * direction(theta, dim=dim)
        t, D, st = m.compute_traction(jump)
        np.testing.assert_allclose(t, K_E * W @ jump, rtol=1e-12, atol=0.0)
        np.testing.assert_allclose(D, K_E * W, rtol=1e-12, atol=0.0)
        assert st["damage"] == pytest.approx(0.0, abs=1e-12)
    jump = 0.5 * m.kappa_0 * direction(60.0, dim=dim)
    jump[0] = -3.0 * m.kappa_0
    t, D, st = m.compute_traction(jump)
    np.testing.assert_allclose(t, K_E * W @ jump, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(D, K_E * W, rtol=1e-12, atol=0.0)
    assert st["kappa"] == m.kappa_0


@pytest.mark.parametrize("beta", [0.7, 1.0, 2.3])
def test_superficie_de_iniciacion(beta):
    m = law(beta=beta)
    for theta in np.linspace(0.0, 90.0, 13):
        t, _, _ = m.compute_traction(m.kappa_0 * direction(theta, beta=beta))
        g = (t[0] / SIGMA_T0) ** 2 + (t[1] / (beta * SIGMA_T0)) ** 2
        assert g == pytest.approx(1.0, rel=1e-10)
    t, _, _ = m.compute_traction(np.array([m.kappa_0, 0.0]))
    assert t[0] == pytest.approx(SIGMA_T0, rel=1e-10)
    t, _, _ = m.compute_traction(np.array([0.0, m.kappa_0 / beta]))
    assert t[1] == pytest.approx(beta * SIGMA_T0, rel=1e-10)
    assert m.shear_strength == pytest.approx(beta * SIGMA_T0, rel=1e-15)


@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("dim", [2, 3])
def test_energia_independiente_del_modo(shape, dim):
    m = law(shape)
    end = 1.001 * end_opening(m)
    # Malla con los quiebres de la envolvente como nodos: el trapecio es exacto
    # en los tramos rectos y la tolerancia absorbe el exponencial.
    marks = [0.0, m.kappa_0]
    if shape == "bilinear":
        marks.append(m.kappa_0 + m._curve.w_1)
    marks.append(end)
    lam = np.unique(np.concatenate([np.linspace(a, b, 3000) for a, b in zip(marks[:-1], marks[1:])]))
    for theta in (0.0, 30.0, 60.0, 90.0):
        d = direction(theta, dim=dim)
        path = lam[:, None] * d
        t, _ = walk(m, path)
        work = np.sum(0.5 * (t[1:] + t[:-1]) * np.diff(path, axis=0))
        assert work == pytest.approx(G_F, rel=1e-3), (shape, dim, theta)
        assert np.linalg.norm(t[-1]) <= 1e-9 * SIGMA_T0 or shape == "exponential"


def test_bilineal_model_code_2010():
    m = law("bilinear")
    k0 = m.kappa_0
    H = G_F - 0.5 * SIGMA_T0 * k0
    w1 = H / SIGMA_T0
    assert m.knee_stress_ratio == 0.2 and m.final_opening_factor == 5.0
    assert m.envelope(k0 + w1)[0] == pytest.approx(0.2 * SIGMA_T0, rel=1e-10)
    for frac in (0.25, 0.5, 0.9):
        assert m.envelope(k0 + frac * w1)[0] == pytest.approx(SIGMA_T0 * (1 - 0.8 * frac), rel=1e-10)
    for ratio in (1.5, 3.0, 4.5):
        assert m.envelope(k0 + ratio * w1)[0] == pytest.approx(SIGMA_T0 * (0.25 - 0.05 * ratio), rel=1e-10)
    assert m.envelope(k0 + 5.0 * w1) == (0.0, 0.0)
    assert m.envelope(k0 + 7.0 * w1) == (0.0, 0.0)
    # Área exacta: la curva es poligonal, el trapecio en los vértices es exacto.
    for r, c in ((0.2, 5.0), (1.0 / 3.0, 3.6), (0.5, 2.0)):
        mm = law("bilinear", knee_stress_ratio=r, final_opening_factor=c)
        kk = mm._curve
        nodes = np.array([0.0, kk.kappa_0, kk.kappa_0 + kk.w_1, kk.kappa_end])
        T = np.array([mm.envelope(k)[0] for k in nodes])
        area = np.sum(0.5 * (T[1:] + T[:-1]) * np.diff(nodes))
        assert area == pytest.approx(G_F, rel=1e-12), (r, c)
        assert mm.envelope(kk.kappa_0 + kk.w_1)[0] == pytest.approx(r * SIGMA_T0, rel=1e-10)
    pet = law("bilinear", knee_stress_ratio=1.0 / 3.0, final_opening_factor=3.6)
    Hp = G_F - 0.5 * SIGMA_T0 * pet.kappa_0
    assert pet._curve.w_1 == pytest.approx(0.8 * Hp / SIGMA_T0, rel=1e-12)


@pytest.mark.parametrize("shape", ["linear", "exponential"])
def test_limite_beta_cero(shape):
    """β = 0: tracción, tangente y estado idénticos a la ley de modo I."""
    rng = np.random.default_rng(7)
    mixed = law(shape, beta=0.0)
    mode1 = CohesiveDamageIsotropic(sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_E, softening=shape)
    for _ in range(200):
        s1 = s2 = None
        for _ in range(8):
            jump = np.array([rng.uniform(-0.5, 1.5), rng.uniform(-2.0, 2.0)]) * rng.choice(
                [mixed.kappa_0, 0.3 * W_REF, 2.0 * W_REF])
            t1, D1, s1 = mixed.compute_traction(jump, s1)
            t2, D2, s2 = mode1.compute_traction(jump, s2)
            scale_t = max(np.abs(t2).max(), 1e-30)
            np.testing.assert_allclose(t1, t2, rtol=1e-12, atol=1e-12 * scale_t)
            np.testing.assert_allclose(D1, D2, rtol=1e-12, atol=1e-12 * max(np.abs(D2).max(), 1e-30))
            assert s1["kappa"] == pytest.approx(s2["kappa"], rel=1e-14)
            assert s1["damage"] == pytest.approx(s2["damage"], rel=1e-12, abs=1e-14)


@pytest.mark.parametrize("shape", SHAPES)
def test_descarga_y_recarga_secantes(shape):
    m = law(shape)
    k1 = 0.3 * (m._curve.kappa_end if shape != "exponential" else 3.0 * W_REF)
    d45 = direction(45.0)
    _, st = walk(m, np.linspace(0.0, 1.0, 50)[1:, None] * k1 * d45)
    assert st["kappa"] == pytest.approx(k1, rel=1e-14)
    S_k1 = m.envelope(k1)[0] / k1
    W = weights(2)
    for frac in np.linspace(1.0, 0.0, 11):            # descarga al origen
        jump = frac * k1 * d45
        t, D, s = m.compute_traction(jump, st)
        np.testing.assert_allclose(t, S_k1 * W @ jump, rtol=1e-10, atol=1e-10 * SIGMA_T0)
        np.testing.assert_allclose(D, S_k1 * W, rtol=1e-10)
        assert s["kappa"] == st["kappa"]
    d80 = direction(80.0)
    t, _, s = m.compute_traction(0.9 * k1 * d80, st)     # recarga por debajo de κ
    assert s["kappa"] == st["kappa"]
    np.testing.assert_allclose(t, S_k1 * W @ (0.9 * k1 * d80), rtol=1e-10)
    t, _, s = m.compute_traction(1.2 * k1 * d80, st)     # y por encima
    assert s["kappa"] == pytest.approx(1.2 * k1, rel=1e-14)
    t_eff = math.hypot(t[0], t[1] / BETA)
    assert t_eff == pytest.approx(m.envelope(1.2 * k1)[0], rel=1e-10)


# --- específicas ----------------------------------------------------------------------

def _states(m, dim):
    """Saltos y estados confirmados que cubren todas las ramas."""
    w = end_opening(m)
    out = []
    for theta in (0.0, 30.0, 60.0, 89.9):
        for frac in (0.01, 0.05, 0.3, 0.7, 0.95):
            jump = frac * w * direction(theta, dim=dim)
            out.append((jump, {"kappa": 0.5 * frac * w}))       # carga activa
            out.append((jump, {"kappa": 2.0 * frac * w}))       # descarga
    for frac in (0.05, 0.3, 0.7):                                 # compresión y deslizamiento
        jump = frac * w * direction(90.0, dim=dim)
        jump[0] = -0.2 * frac * w
        out.append((jump, {"kappa": 0.5 * frac * w}))
        out.append((jump, {"kappa": 2.0 * frac * w}))
    out.append((0.3 * m.kappa_0 * direction(30.0, dim=dim), None))  # elástica
    return out


@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("dim", [2, 3])
def test_tangente_simetrica(shape, dim):
    m = law(shape)
    for jump, st in _states(m, dim):
        _, D, _ = m.compute_traction(jump, st)
        assert np.linalg.norm(D - D.T) <= 1e-14 * np.linalg.norm(D)
    w = 1.01 * m._curve.kappa_end if shape != "exponential" else 50 * W_REF
    _, D, _ = m.compute_traction(w * direction(30.0, dim=dim))      # totalmente abierta
    assert np.linalg.norm(D - D.T) == 0.0


@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("dim", [2, 3])
def test_tangente_diferencia_finita(shape, dim):
    m = law(shape)
    for jump, st in _states(m, dim):
        _, D, _ = m.compute_traction(jump, st)
        h = 1e-7 * np.linalg.norm(jump)
        D_fd = np.zeros_like(D)
        for j in range(dim):
            e = np.zeros(dim)
            e[j] = h
            D_fd[:, j] = (m.compute_traction(jump + e, st)[0] - m.compute_traction(jump - e, st)[0]) / (2 * h)
        scale = np.linalg.norm(D_fd)
        if scale == 0.0:                    # grieta totalmente abierta: D = 0 exacto
            assert np.linalg.norm(D) == 0.0
            continue
        err = np.linalg.norm(D - D_fd) / scale
        assert err < 1e-6, (shape, dim, jump, st, err)


@pytest.mark.parametrize("shape", SHAPES)
def test_compresion_sin_dano_normal(shape):
    m = law(shape)
    w = end_opening(m)
    for kappa in (None, {"kappa": 0.3 * w}, {"kappa": 1.5 * w}):
        for dn in (-1e-3 * m.kappa_0, -m.kappa_0, -50.0 * m.kappa_0):
            t, D, st = m.compute_traction(np.array([dn, 0.0]), kappa)
            assert t[0] == pytest.approx(K_E * dn, rel=1e-14)
            assert D[0, 0] == K_E and D[0, 1] == 0.0 and D[1, 0] == 0.0
            expected = m.kappa_0 if kappa is None else kappa["kappa"]
            assert st["kappa"] == expected
        # Continuidad de t_n en δ_n = 0.
        tiny = 1e-12 * m.kappa_0
        tp = m.compute_traction(np.array([tiny, 0.0]), kappa)[0][0]
        tm = m.compute_traction(np.array([-tiny, 0.0]), kappa)[0][0]
        assert abs(tp - tm) <= 2.0 * K_E * tiny


@pytest.mark.parametrize("shape", SHAPES)
def test_deslizamiento_bajo_compresion_dana_sin_friccion(shape):
    m = law(shape)
    end = 1.001 * end_opening(m) / BETA
    slips = np.linspace(0.0, end, 400)[1:]
    for dn in (-0.5 * m.kappa_0, -20.0 * m.kappa_0):
        st = None
        for ds in slips:
            t, _, st = m.compute_traction(np.array([dn, ds]), st)
            assert t[0] == pytest.approx(K_E * dn, rel=1e-14)
            assert st["kappa"] == pytest.approx(max(BETA * ds, m.kappa_0), rel=1e-14)
            expected = BETA * m.envelope(BETA * ds)[0] if BETA * ds > m.kappa_0 else BETA ** 2 * K_E * ds
            assert t[1] == pytest.approx(expected, rel=1e-10, abs=1e-10 * SIGMA_T0)
        if shape != "exponential":
            assert t[1] == 0.0


@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("dim", [2, 3])
def test_disipacion_no_negativa(shape, dim):
    """Trabajo menos energía secante almacenada: no negativo y no decreciente."""
    m = law(shape, K_e=1.0e13)
    rng = np.random.default_rng(11 + dim)
    w = end_opening(m)
    tol = 2e-3 * G_F
    for _ in range(1000 // 6 + 1):                 # 1000 historias entre las 6 variantes
        state, start, j_prev, t_prev = None, np.zeros(dim), np.zeros(dim), np.zeros(dim)
        work = last = 0.0
        for _ in range(6):
            target = rng.uniform(-1.0, 1.0, dim) * rng.choice([0.05, 0.4, 1.2]) * w
            target[0] = max(target[0], -3.0 * m.kappa_0)
            steps = list(np.linspace(0.0, 1.0, 41)[1:])
            if start[0] * target[0] < 0.0:
                # Nodo en el cierre δ_n = 0: ahí la rigidez normal salta a la
                # penalización K_e y el trapecio que lo cruza no es exacto.
                steps.append(start[0] / (start[0] - target[0]))
            for a in sorted(steps):
                jump = start + a * (target - start)
                t, _, state = m.compute_traction(jump, state)
                work += 0.5 * (t + t_prev) @ (jump - j_prev)
                t_prev, j_prev = t, jump
            start = target
            dissipated = work - 0.5 * t @ jump      # energía secante almacenada ½·t·δ
            assert dissipated >= -tol
            assert dissipated >= last - tol
            last = dissipated


def test_isotropia_en_el_plano_tangente():
    m = law("bilinear")
    w = end_opening(m)
    base = np.array([0.2 * w, 0.25 * w, 0.0])
    st0 = {"kappa": 0.1 * w}
    t0, D0, _ = m.compute_traction(base, st0)
    for phi in (0.3, 1.1, 2.5, -0.8):
        R = np.array([[1.0, 0.0, 0.0],
                      [0.0, math.cos(phi), -math.sin(phi)],
                      [0.0, math.sin(phi), math.cos(phi)]])
        t, D, _ = m.compute_traction(R @ base, st0)
        np.testing.assert_allclose(t, R @ t0, rtol=1e-12, atol=1e-12 * SIGMA_T0)
        np.testing.assert_allclose(D, R @ D0 @ R.T, rtol=1e-12, atol=1e-12 * np.abs(D0).max())


def test_rechazos():
    for kw in ({"sigma_t0": 0.0}, {"G_f": -1.0}, {"K_e": 0.0}, {"beta": -0.1},
               {"beta": math.inf}, {"softening": "parabolic"},
               {"K_e": SIGMA_T0 ** 2 / (2 * G_F)},
               {"softening": "linear", "knee_stress_ratio": 0.2},
               {"softening": "exponential", "final_opening_factor": 5.0},
               {"knee_stress_ratio": 1.0}, {"knee_stress_ratio": 0.0},
               {"knee_stress_ratio": 0.5, "final_opening_factor": 4.0},     # r·c = 2
               {"knee_stress_ratio": 0.2, "final_opening_factor": 1.5}):    # c(1+r) < 2
        args = dict(sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_E, beta=BETA, softening="bilinear")
        args.update(kw)
        with pytest.raises(ValueError):
            CohesiveDamageMixedMode(**args)
    m = law()
    for bad in (np.zeros(1), np.zeros(4), np.zeros((2, 2))):
        with pytest.raises(ValueError, match="componentes"):
            m.compute_traction(bad)


# --- arquitectura ---------------------------------------------------------------------

def test_registro_y_contrato():
    cls = InterfaceMaterialRegistry.get("CohesiveDamageMixedMode")
    assert cls is CohesiveDamageMixedMode
    assert cls.IS_SYMMETRIC is True
    assert cls.GEOMETRY == "surface" and cls.PRIMARY_STATE_VAR == "damage"
    m = InterfaceMaterialRegistry.create("CohesiveDamageMixedMode", sigma_t0=SIGMA_T0, G_f=G_F,
                                         K_e=K_E, beta=BETA, softening="bilinear")
    assert m.accepts("surface", 2) and m.accepts("surface", 3)
    assert not m.accepts("line", 2)
    _, _, st = m.compute_traction(np.array([1e-5, 1e-5]))
    assert set(st) == {"kappa", "damage"}
