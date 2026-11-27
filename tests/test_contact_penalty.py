"""Aceptación de ``ContactPenalty`` (spec docs/specs/ContactPenalty.md).

Contacto unilateral por penalización ``t_n = K_n·min(g_n, 0)`` con fricción de
Coulomb opcional regularizada por la adherencia ``K_t`` y retorno radial. Los
casos siguen el bloque ``acceptance`` de la spec; los parámetros son ≠ 1
(Reglas §6) para que un factor olvidado no pase como identidad.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from solidum.materials.interface import ContactPenalty
from solidum.registry import InterfaceMaterialRegistry

K_N = 3.7e4
K_T = 1.3e4
MU = 0.45
G_REF = 1.0e-3          # escala de holguras


def law(mu=MU, K_n=K_N, K_t=K_T):
    if mu == 0.0:
        return ContactPenalty(K_n=K_n)
    return ContactPenalty(K_n=K_n, mu=mu, K_t=K_t)


def walk(m, path, state=None):
    """Recorre una trayectoria confirmando el estado en cada punto."""
    out = []
    for jump in path:
        t, T, state = m.compute_traction(np.asarray(jump, dtype=float), state)
        out.append((t, T, state))
    return out


# -- normal ------------------------------------------------------------------

@pytest.mark.parametrize("dim", [2, 3])
def test_unilateral_normal(dim):
    m = law()
    for g_n in (-2.0, -1.0e-6, 1.0e-6, 3.0):
        g = np.zeros(dim)
        g[0] = g_n * G_REF
        t, T, st = m.compute_traction(g)
        if g_n < 0.0:
            assert t[0] == pytest.approx(K_N * g[0], rel=1e-14)
            assert T[0, 0] == K_N
            assert st["contact_pressure"] == pytest.approx(-K_N * g[0], rel=1e-14)
        else:
            assert np.all(t == 0.0) and np.all(T == 0.0)
            assert st["status"] == 0.0
    # Holgura nula: fuerza nula y rama cerrada (rigidez de arranque, §11).
    t, T, st = m.compute_traction(np.zeros(dim))
    assert np.all(t == 0.0)
    assert T[0, 0] == K_N
    assert st["status"] == 1.0


def test_interpenetracion_escala_con_1_sobre_Kn():
    p = 2.5e3
    for K_n in (1.0e4, 1.0e6, 1.0e8):
        m = ContactPenalty(K_n=K_n)
        g_n = -p / K_n
        t, _, st = m.compute_traction(np.array([g_n, 0.0]))
        assert st["contact_pressure"] == pytest.approx(p, rel=1e-14)


# -- fricción -----------------------------------------------------------------

def test_adherencia_y_umbral_de_coulomb():
    m = law()
    g_n = -0.8 * G_REF
    p = -K_N * g_n
    threshold = MU * p / K_T
    state = None
    for g_s in np.linspace(0.0, 5.0 * threshold, 41):
        t, T, state = m.compute_traction(np.array([g_n, g_s]), state)
        if g_s <= threshold:
            assert t[1] == pytest.approx(K_T * g_s, rel=1e-12, abs=1e-12 * MU * p)
            assert state["status"] == 1.0
        else:
            assert t[1] == pytest.approx(MU * p, rel=1e-12)
            assert state["slip"][0] == pytest.approx(g_s - threshold, rel=1e-12)
            assert state["status"] == 2.0


def _cycle_work(m, g_n, path):
    """Trabajo tangencial exacto (la tracción es lineal a trozos en g_s y la
    trayectoria pasa por los quiebres) y estado final."""
    state = None
    t_prev = np.zeros(2)
    g_prev = np.array([g_n, 0.0])
    m.compute_traction(g_prev)
    work = 0.0
    for g_s in path:
        g = np.array([g_n, g_s])
        t, _, state = m.compute_traction(g, state)
        work += 0.5 * (t[1] + t_prev[1]) * (g[1] - g_prev[1])
        t_prev, g_prev = t, g
    return work, state


def test_ciclo_de_friccion():
    """Disipación por ciclo estacionario 4·μ·p·(a − μ·p/K_t) y trabajo = ψ + D."""
    m = law()
    g_n = -0.6 * G_REF
    p = -K_N * g_n
    s_y = MU * p / K_T
    a = 3.7 * s_y
    # Trayectoria que pasa por los quiebres: 0 → a → −a → a → −a. El primer
    # tramo deja de adherir en s_y; tras una inversión, en 2·s_y.
    def leg(x0, x1, elastic=2.0 * s_y):
        brk = x0 + math.copysign(elastic, x1 - x0)
        pts = [brk] if min(x0, x1) < brk < max(x0, x1) else []
        return list(np.linspace(x0, x1, 9)[1:]) + pts
    first = sorted(leg(0.0, a, elastic=s_y))
    path_1 = first + sorted(leg(a, -a), reverse=True) + sorted(leg(-a, a))
    w1, st1 = _cycle_work(m, g_n, path_1)
    path_2 = path_1 + sorted(leg(a, -a), reverse=True) + sorted(leg(-a, a))
    w2, st2 = _cycle_work(m, g_n, path_2)
    # Ciclo completo estacionario (de a a a): sólo se disipa (ψ igual al
    # empezar y al terminar).
    assert (w2 - w1) == pytest.approx(4.0 * MU * p * (a - s_y), rel=1e-10)
    # Trabajo = ψ + D en todo el recorrido (D = μ·p·slip acumulado).
    psi = m.stored_energy(np.array([g_n, a]), st2) - 0.5 * K_N * g_n ** 2
    D = MU * p * st2["slip_accumulated"]
    assert w2 == pytest.approx(psi + D, rel=1e-10)
    assert D > 0.0


def test_paso_unico_que_cruza_ramas():
    m = law()
    # De abierta a cerrada y de adherencia a deslizamiento en un solo paso.
    g = np.array([-0.5 * G_REF, 2.0 * G_REF])
    t, T, st = m.compute_traction(g, {"slip": np.array([0.3 * G_REF])})
    p = 0.5 * G_REF * K_N
    t_trial = K_T * (2.0 - 0.3) * G_REF
    assert t_trial > MU * p
    assert t[0] == pytest.approx(-p, rel=1e-14)
    assert t[1] == pytest.approx(MU * p, rel=1e-14)
    assert st["slip"][0] == pytest.approx(2.0 * G_REF - MU * p / K_T, rel=1e-12)
    assert st["status"] == 2.0


def test_presion_decreciente_hasta_abrir():
    m = law()
    state = None
    path_n = np.linspace(-1.0, 0.2, 61) * G_REF
    prev = None
    for k, g_n in enumerate(path_n):
        g_s = 2.0 * G_REF * (1 + k / 10)          # deslizamiento sostenido
        t, _, state = m.compute_traction(np.array([g_n, g_s]), state)
        if g_n > 0:
            assert np.all(t == 0.0)
            assert state["slip"][0] == g_s
        else:
            assert abs(t[1]) == pytest.approx(MU * K_N * (-g_n), rel=1e-12, abs=1e-12)
            if prev is not None:
                assert abs(t[1]) <= prev + 1e-12     # baja de forma continua
            prev = abs(t[1])


def test_separacion_y_recontacto():
    m = law()
    g_n = -0.7 * G_REF
    st = None
    # Cerrar y deslizar.
    for g_s in np.linspace(0, 3 * G_REF, 7):
        _, _, st = m.compute_traction(np.array([g_n, g_s]), st)
    # Abrir y moverse tangencialmente.
    t, _, st = m.compute_traction(np.array([0.4 * G_REF, 5.0 * G_REF]), st)
    assert np.all(t == 0.0) and st["slip"][0] == 5.0 * G_REF
    # Volver a cerrar en la nueva posición: adherencia desde ahí.
    t, _, st = m.compute_traction(np.array([g_n, 5.0 * G_REF]), st)
    assert t[1] == 0.0
    assert st["status"] == 1.0


@pytest.mark.parametrize("theta_deg", [0.0, 37.0, 140.0])
def test_isotropia_tangencial_3d(theta_deg):
    """Trayectoria circular bajo presión constante; resultado invariante ante
    una rotación de los ejes (s1, s2) del plano tangente."""
    m = law()
    g_n = -0.9 * G_REF
    p = -K_N * g_n
    th = math.radians(theta_deg)
    rot = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    radius = 4.0 * MU * p / K_T
    st_a = st_b = None
    for phi in np.linspace(0, 2 * math.pi, 50):
        gs = radius * np.array([math.cos(phi), math.sin(phi)]) + np.array([radius, 0])
        ta, _, st_a = m.compute_traction(np.r_[g_n, gs], st_a)
        tb, _, st_b = m.compute_traction(np.r_[g_n, rot @ gs], st_b)
        np.testing.assert_allclose(tb[1:], rot @ ta[1:], rtol=1e-12, atol=1e-12 * MU * p)
        assert np.linalg.norm(ta[1:]) <= MU * p * (1 + 1e-12)
        if st_a["status"] == 2.0:
            assert np.linalg.norm(ta[1:]) == pytest.approx(MU * p, rel=1e-12)


def test_cambio_de_direccion_3d_converge_con_orden_1():
    """Deslizar en s1 y después girar 90° hacia s2 en n subpasos: la solución
    se acerca a la integración fina con error de orden 1 (spec §8)."""
    m = law()
    g_n = -0.5 * G_REF
    p = -K_N * g_n
    s_y = MU * p / K_T
    start = np.array([3 * s_y, 0.0])
    end = np.array([3 * s_y, 3 * s_y])

    def run(nsub):
        _, _, st = m.compute_traction(np.r_[g_n, start])
        t = None
        for k in range(1, nsub + 1):
            gs = start + (end - start) * k / nsub
            t, _, st = m.compute_traction(np.r_[g_n, gs], st)
        return t[1:]

    ref = run(4096)
    errs = [np.linalg.norm(run(n) - ref) for n in (8, 16, 32, 64)]
    rates = [math.log2(errs[i] / errs[i + 1]) for i in range(3)]
    assert all(0.8 < r < 1.2 for r in rates), rates


# -- tangente -----------------------------------------------------------------

def _fd_tangent(m, g, state, h):
    n = g.size
    T = np.zeros((n, n))
    for j in range(n):
        e = np.zeros(n)
        e[j] = h
        tp, _, _ = m.compute_traction(g + e, state)
        tm, _, _ = m.compute_traction(g - e, state)
        T[:, j] = (tp - tm) / (2 * h)
    return T


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("branch", ["abierta", "adherencia", "deslizamiento"])
def test_tangente_diferencia_finita(dim, branch):
    m = law()
    slip0 = np.full(dim - 1, 0.2 * G_REF)
    state = {"slip": slip0}
    if branch == "abierta":
        g = np.r_[0.5 * G_REF, np.full(dim - 1, 0.7 * G_REF)]
    elif branch == "adherencia":
        g = np.r_[-1.0 * G_REF, slip0 + 0.1 * MU * K_N * G_REF / K_T / dim]
    else:
        g = np.r_[-1.0 * G_REF, slip0 + np.linspace(3.0, 4.5, dim - 1) * MU * K_N * G_REF / K_T]
    _, T, st = m.compute_traction(g, state)
    assert st["status"] == {"abierta": 0.0, "adherencia": 1.0, "deslizamiento": 2.0}[branch]
    T_fd = _fd_tangent(m, g, state, 1e-7 * G_REF)
    np.testing.assert_allclose(T, T_fd, rtol=1e-6, atol=1e-6 * K_N)


def test_simetria_de_la_tangente():
    m = law()
    g = np.array([-G_REF, 4.0 * MU * K_N * G_REF / K_T])
    _, T, _ = m.compute_traction(g)
    assert T[1, 0] == pytest.approx(-MU * K_N, rel=1e-14)
    assert T[0, 1] == 0.0
    assert not m.IS_SYMMETRIC
    free = law(mu=0.0)
    _, T0, _ = free.compute_traction(g)
    assert np.array_equal(T0, T0.T) and free.IS_SYMMETRIC


# -- sin fricción, energía, contrato --------------------------------------------

@pytest.mark.parametrize("dim", [2, 3])
def test_sin_friccion(dim):
    m = ContactPenalty(K_n=K_N)
    g = np.r_[-G_REF, np.full(dim - 1, 7.0 * G_REF)]
    t, T, st = m.compute_traction(g)
    assert np.all(t[1:] == 0.0)
    expected = np.zeros((dim, dim))
    expected[0, 0] = K_N
    assert np.array_equal(T, expected)
    assert m.IS_SYMMETRIC
    with pytest.raises(ValueError, match="K_t sólo interviene"):
        ContactPenalty(K_n=K_N, K_t=K_T)


def test_energia_almacenada():
    m = law()
    g_n = -0.4 * G_REF
    # Trayectoria sin deslizamiento: ψ = trabajo exacto (ley lineal a trozos).
    g = np.array([g_n, 0.5 * MU * K_N * 0.4 * G_REF / K_T])
    _, _, st = m.compute_traction(g)
    work = 0.5 * K_N * g_n ** 2 + 0.5 * K_T * g[1] ** 2
    assert m.stored_energy(g, st) == pytest.approx(work, rel=1e-12)
    # Holgura abierta en U = 0: ψ = 0.
    g_open = np.array([2.0 * G_REF, 0.0])
    _, _, st = m.compute_traction(g_open)
    assert m.stored_energy(g_open, st) == 0.0


def test_escalado_con_los_parametros():
    alpha, beta = 7.3, 0.37
    g = np.array([-0.8 * G_REF, 6.0 * MU * K_N * 0.8 * G_REF / K_T])
    t1, _, _ = law().compute_traction(g)
    t2, _, _ = law(K_n=alpha * K_N, K_t=alpha * K_T).compute_traction(g / alpha)
    np.testing.assert_allclose(t2, t1, rtol=1e-12)
    t3, _, st3 = law(mu=beta * MU).compute_traction(g)
    assert st3["status"] == 2.0
    assert t3[1] == pytest.approx(beta * t1[1], rel=1e-12)


def test_prueba_y_confirmacion():
    m = law()
    committed = {"slip": np.array([0.1 * G_REF]), "slip_accumulated": 0.3 * G_REF}
    g = np.array([-G_REF, 3.0 * G_REF])
    a = m.compute_traction(g, committed)
    b = m.compute_traction(g, committed)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])
    assert a[2]["slip_accumulated"] == b[2]["slip_accumulated"]
    np.testing.assert_array_equal(a[2]["slip"], b[2]["slip"])
    assert committed["slip"][0] == 0.1 * G_REF        # el estado confirmado no cambia


def test_rechazos():
    with pytest.raises(ValueError, match="K_n"):
        ContactPenalty(K_n=0.0)
    with pytest.raises(ValueError, match="mu"):
        ContactPenalty(K_n=K_N, mu=-0.1)
    with pytest.raises(ValueError, match="hace falta K_t"):
        ContactPenalty(K_n=K_N, mu=0.3)
    with pytest.raises(ValueError, match="componentes"):
        ContactPenalty(K_n=K_N).compute_traction(np.zeros(4))


def test_registro_y_contrato():
    assert InterfaceMaterialRegistry.get("ContactPenalty") is ContactPenalty
    m = ContactPenalty(K_n=K_N)
    assert m.accepts("surface", 2) and m.accepts("surface", 3)
    assert not m.accepts("line", 2)
    assert ContactPenalty.IS_UNILATERAL and ContactPenalty.ACCEPTS_INITIAL_GAP
    assert not ContactPenalty.SECANT_UNLOADING
    assert ContactPenalty(K_n=K_N, mu=0.2, K_t=K_T).IS_SYMMETRIC is False
