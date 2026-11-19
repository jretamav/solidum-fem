"""Seguimiento de grieta del módulo de discontinuidades.

Specs ``CrackTracking`` y ``PrescribedCrackPath`` (aceptación), y la
aceptación ``tres_triangulos_traccion_hasta_separacion`` de ``CST_Embedded2D``.

Bloques:

- **Regla de la punta construida a mano**: lado del tercer nodo por la
  prolongación recta del último tramo, magnitudes de los valores nodales,
  salida por el borde libre (con placas, otro seguimiento, elementos térmicos,
  esquina y placa añadida después), propiedad «la prolongación entra siempre
  en la punta».
- **Sólo el seguimiento con un campo impuesto**: equivalencia bit a bit con el
  prototipo aprobado el 2026-11-19 en 15 mallas con trayectorias en arco
  (huellas de la secuencia de activaciones), limitación conocida en arcos y en
  rectas con campo uniforme (§6 de la spec), giro de las direcciones con media
  grieta formada y ruido en el campo.
- **Con el elemento real** y ``NonlinearSolver``: franjas en tracción con
  mallas no alineadas, mallas con astillas (``cos ψ < 0``, trayectoria
  prescrita y regla del vértice), dos grietas y un triángulo aislado hasta la
  separación.
- **Plumbing**: grupos por omisión, una evaluación por ronda, esfuerzos
  vectorizados, YAML, nombres retirados, serialización.
"""
from __future__ import annotations

import json
import logging
import pickle
import textwrap
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial import Delaunay

import solidum
from solidum import Assembler, Domain
from solidum.elements.solid_2d import Quad4, compute_kinematics_tri3
from solidum.elements.solid_2d.tri3 import Tri3
from solidum.elements.thermal import Quad4Thermal
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.math.solvers import NonlinearSolver

solidum.load_user_module("discontinuities")
from solidum.materials.interface import CohesiveDamageIsotropic  # noqa: E402
from solidum.user.discontinuities import (  # noqa: E402
    CST_Embedded2D,
    CrackTracking,
    PrescribedCrackPath,
)
from solidum.user.discontinuities import embedded_cst  # noqa: E402
from solidum.user.discontinuities.crack_tracking import (  # noqa: E402
    default_tracker,
    partition,
    polyline_signed_distance,
    root_partition,
    triangle_gradients,
)

_EDGES = ((0, 1), (1, 2), (2, 0))
E, NU, ST0, GF, TH = 30.0e9, 0.0, 3.0e6, 100.0, 0.02


# =========================================================================
# Mallas
# =========================================================================


def _ccw(pts, tris):
    out = []
    for t in tris:
        c = pts[t]
        two_a = (c[1, 0] - c[0, 0]) * (c[2, 1] - c[0, 1]) - (c[2, 0] - c[0, 0]) * (c[1, 1] - c[0, 1])
        if abs(two_a) > 1e-14:
            out.append(t if two_a > 0 else t[[0, 2, 1]])
    return np.array(out)


def _structured(L, W, nx, ny, skew=0.0):
    xs, ys = np.linspace(0, L, nx + 1), np.linspace(0, W, ny + 1)
    pts = np.array([[x + skew * (y / W - 0.5) * (L / nx) * (0 < x < L), y] for y in ys for x in xs])
    tris = []
    for j in range(ny):
        for i in range(nx):
            a = j * (nx + 1) + i
            tris += [[a, a + 1, a + nx + 2], [a, a + nx + 2, a + nx + 1]]
    return pts, _ccw(pts, np.array(tris))


def _delaunay(L, W, nx, ny, seed, jitter=0.3):
    rng = np.random.default_rng(seed)
    xs, ys = np.linspace(0, L, nx + 1), np.linspace(0, W, ny + 1)
    hx, hy = L / nx, W / ny
    pts = np.array([[x + (jitter * hx * rng.uniform(-1, 1) if 0 < x < L else 0.0),
                     y + (jitter * hy * rng.uniform(-1, 1) if 0 < y < W else 0.0)]
                    for y in ys for x in xs])
    return pts, _ccw(pts, Delaunay(pts).simplices)


def _coh(sigma_t0=ST0):
    return CohesiveDamageIsotropic(sigma_t0=sigma_t0, G_f=GF, K_e=1.0e15, softening="linear")


# =========================================================================
# Regla de la punta construida a mano (spec §1–§3)
# =========================================================================
#
# p = (a, b, r) agrietado a mano con la recta y = 0.45 (normal hacia a):
# valores a = +0.45, b = −0.55, r = −0.05. Cruces P = (1, 0.45) en a–b y
# Q = (0.1, 0.45) en a–r: el último tramo es horizontal. La punta es
# q = (a, c, b), con c = (2, y_c) a la derecha del lado a–b.

_A, _B, _R = (1.0, 0.0), (1.0, 1.0), (0.0, 0.5)


def _tip(y_c, extra=None, values=None):
    """Seguimiento con p agrietado y la punta q; devuelve (seguimiento,
    dominio, valor de punta de c). ``extra(dom, n)`` añade elementos (n: dict
    de nodos por nombre); ``values`` sustituye los valores nodales de p."""
    dom = Domain()
    n = {k: dom.add_node(i, list(x)) for i, (k, x) in
         enumerate((("a", _A), ("b", _B), ("r", _R), ("c", (2.0, y_c))), 1)}
    tracker = CrackTracking()
    bulk = Elastic2D(E=E, nu=0.2, hypothesis="plane_stress")
    p = CST_Embedded2D(1, [n["a"], n["b"], n["r"]], bulk, _coh(), crack_tracker=tracker)
    q = CST_Embedded2D(2, [n["a"], n["c"], n["b"]], bulk, _coh(), crack_tracker=tracker)
    dom.add_element(p)
    dom.add_element(q)
    if extra is not None:
        extra(dom, n, tracker)
    dom.generate_equation_numbers(verbose=False)
    p.activate(np.array([0.0, -1.0]), 0, np.array([0.5, 0.45]))
    U = np.zeros(dom.total_dofs)
    tracker.update(U)                      # adopta p; nada supera σ_t0
    if values is not None:
        tracker.cracks[0].fixed.update(values)
        tracker._built = False
        tracker.update(U)
    tracker._refresh(0)
    return tracker, dom, float(tracker._values[0][tracker._nidx[n["c"].id]])


def test_valores_nodales_de_la_grieta_formada():
    tracker, _, _ = _tip(0.3)
    np.testing.assert_allclose([tracker.cracks[0].fixed[i] for i in (1, 2, 3)], [0.45, -0.55, -0.05],
                               atol=1e-12)


@pytest.mark.parametrize("y_c,esperado", [(0.3, +0.15), (0.6, -0.15)])
def test_el_tercer_nodo_va_con_el_lado_de_la_prolongacion(y_c, esperado):
    """acceptance.regla_de_la_punta: c debajo de la prolongación y = 0.45 va
    con a (+), encima, con b (−); el valor es su distancia a la recta."""
    tracker, _, v = _tip(y_c)
    assert v == pytest.approx(esperado, abs=1e-12)
    assert tracker._tips[0] == {tracker._pos[2]}


def test_tercer_nodo_sobre_la_prolongacion_lleva_el_minimo():
    """c exactamente sobre la recta: lado positivo de d × (x − P) (el de b,
    negativo aquí), con el mínimo δ = 1e-3·h_q."""
    tracker, _, v = _tip(0.45)
    q = tracker._pos[2]
    assert v == pytest.approx(-1e-3 * np.sqrt(2.0 * tracker._area[q]), rel=1e-12)


def test_tercer_nodo_ya_guardado_no_cambia():
    """Si c ya tiene valor en la grieta (otro agrietado lo fijó), la punta no
    lo decide."""
    def extra(dom, n, tracker):
        g = dom.add_node(5, [1.0, -1.0])
        s = CST_Embedded2D(3, [n["a"], g, n["c"]], Elastic2D(E=E, nu=0.2, hypothesis="plane_stress"),
                           _coh(), crack_tracker=tracker)
        dom.add_element(s)
        s.activate(np.array([1.0, 0.0]), 2, np.array([1.5, 0.0]))   # x = 1.5 aísla c
    tracker, _, v = _tip(0.3, extra)
    assert v == tracker.cracks[0].fixed[4]          # el guardado, no la prolongación (+0.15)
    assert v == pytest.approx(-0.5, abs=1e-12)


def test_las_magnitudes_de_los_valores_forman_parte_de_la_regla():
    """acceptance.magnitudes_de_los_valores: con los valores de p reducidos a
    ±1 el último tramo deja de ser horizontal (P = (1, 0.5), Q = (0.5, 0.25))
    y c = (2, 0.6) cambia de lado."""
    _, _, v = _tip(0.6)
    assert v < 0.0
    _, _, v_signos = _tip(0.6, values={1: 1.0, 2: -1.0, 3: -1.0})
    assert v_signos > 0.0


def _x2(u, v):
    return u[0] * v[1] - u[1] * v[0]


def test_la_prolongacion_entra_siempre_en_la_punta():
    """acceptance.la_prolongacion_entra_en_la_punta (spec §2): con el lado de
    entrada cruzado y |v| ≥ δ, la cuerda de los cruces apunta hacia el
    tercer nodo: d·m > 0."""
    rng = np.random.default_rng(11)
    worst = np.inf
    for _ in range(20000):
        L, A, B, c = rng.uniform(-1, 1, (4, 2))
        if abs(_x2(A - L, B - L)) < 1e-3 or _x2(A - L, c - L) * _x2(A - L, B - L) >= 0:
            continue                                  # c al otro lado de L–A respecto de B
        vL = rng.uniform(1e-3, 1.0)
        vA, vB = -rng.uniform(1e-3, 1.0), -rng.uniform(1e-3, 1.0)
        P = L + vL / (vL - vA) * (A - L)
        Q = L + vL / (vL - vB) * (B - L)
        d = (P - Q) / np.linalg.norm(P - Q)
        e = A - L
        m = np.array([-e[1], e[0]])
        m = m if m @ (c - L) > 0 else -m
        worst = min(worst, float(d @ (m / np.linalg.norm(m))))
    assert worst > 0.0


# --- Regla del vértice en la punta, mínimo δ al activar, tercer nodo de dos
# puntas, cambio de lado de una sola evaluación (casos de la verificación del
# teórico del 2026-11-19, guiones en paso_b_impl/teorico_commit/c2.py y c8.py).


def _uniform_field(dom, phi_deg, eps=2.0e-4):
    phi = np.radians(phi_deg)
    nrm = np.array([np.cos(phi), np.sin(phi)])
    U = np.zeros(dom.total_dofs)
    for nd in dom.nodes.values():
        if "ux" in nd.dofs:
            x = np.asarray(nd.coordinates, dtype=float)[:2]
            U[nd.dofs["ux"]], U[nd.dofs["uy"]] = eps * nrm * (x @ nrm)
    return U


def _expected_partition(xs, th, normal):
    """Partición de §3 calculada aparte: (solitario, n, x_q, cos ψ)."""
    s = th >= 0.0
    lone = int(np.argmax(s)) if s.sum() == 1 else int(np.argmin(s))
    oth = [j for j in range(3) if j != lone]
    cr = [xs[lone] + th[lone] / (th[lone] - th[j]) * (xs[j] - xs[lone]) for j in oth]
    xm = 0.5 * (cr[0] + cr[1])
    n = np.asarray(normal, dtype=float) / np.linalg.norm(normal)
    n = n if n @ (xs[lone] - xm) >= 0.0 else -n
    e = xs[oth[1]] - xs[oth[0]]
    ne = np.array([-e[1], e[0]]) / np.linalg.norm(e)
    ne = ne if ne @ (xs[lone] - xs[oth[0]]) >= 0.0 else -ne
    return lone, n, xm, float(n @ ne)


def _delta(xs):
    two_a = abs((xs[1, 0] - xs[0, 0]) * (xs[2, 1] - xs[0, 1]) - (xs[2, 0] - xs[0, 0]) * (xs[1, 1] - xs[0, 1]))
    return 1e-3 * np.sqrt(two_a)


def test_regla_del_vertice_en_la_punta():
    """acceptance.regla_del_vertice_en_la_punta: con tracción a 21°, la
    partición de la punta con el valor de la prolongación (c = +0.15) da
    cos ψ = −0.075; c pasa al lado del otro nodo de entrada (−δ_q), la punta
    se activa con cos ψ = 0.829 y solitario a (no b), y c guarda la distancia
    a la recta de q con ese signo, no δ_q (§5)."""
    tracker, dom, v = _tip(0.3)
    assert v == pytest.approx(0.15, abs=1e-12)
    xs = np.array([_A, (2.0, 0.3), _B])
    phi = np.radians(21.0)
    nrm = np.array([np.cos(phi), np.sin(phi)])
    assert _expected_partition(xs, np.array([0.45, 0.15, -0.55]), nrm)[3] < 0.0
    delta = _delta(xs)
    lone, n, xm, cos_psi = _expected_partition(xs, np.array([0.45, -delta, -0.55]), nrm)
    assert lone == 0 and cos_psi > 0.0
    tracker.update(_uniform_field(dom, 21.0))
    ds = dom.elements[2].discontinuity_state
    assert ds is not None and ds.solitary_node == lone
    np.testing.assert_allclose(ds.normal, n, atol=1e-12)
    np.testing.assert_allclose(ds.centroid, xm, atol=1e-12)
    stored = tracker.cracks[0].fixed[4]
    assert stored == pytest.approx(-abs((xs[1] - xm) @ n), abs=1e-12)
    assert abs(stored) > 100 * delta
    assert tracker.unactivated_warnings == 0


def test_minimo_delta_al_activar():
    """Si la recta de la punta pasa exactamente por c, c guarda +δ_q (§1, §3):
    sin el mínimo guardaría ~1e-17 y los cruces siguientes caerían sobre c."""
    tracker, dom, _ = _tip(0.3)
    xs = np.array([_A, (2.0, 0.3), _B])
    _, _, xm, _ = _expected_partition(xs, np.array([0.45, 0.15, -0.55]), np.array([1.0, 0.0]))
    t_ = xs[1] - xm
    phi_deg = np.degrees(np.arctan2(t_[0], -t_[1])) % 180.0         # normal ⟂ (c − x_q)
    tracker.update(_uniform_field(dom, phi_deg))
    assert dom.elements[2].discontinuity_state is not None
    assert tracker.cracks[0].fixed[4] == pytest.approx(_delta(xs), rel=1e-12)


def _fan():
    """Abanico alrededor de a: agrietados p1 = (a, b, r) y p2 = (a, c2, b),
    puntas q1 = (a, r, c) y q2 = (a, c, c2) que comparten c; placas Tri3 en
    los lados exteriores r–c y c–c2 para que decida la prolongación."""
    dom = Domain()
    X = dict(a=(0.0, 0.0), b=(0.0, 1.0), r=(-1.0, 0.0), c2=(1.0, 0.0), c=(-1.0, -0.5),
             e1=(-2.0, -0.25), e2=(0.5, -1.5))
    n = {k: dom.add_node(i, list(x)) for i, (k, x) in enumerate(X.items(), 1)}
    tracker = CrackTracking()
    bulk = Elastic2D(E=E, nu=0.2, hypothesis="plane_stress")
    mk = lambda i, names: CST_Embedded2D(i, [n[k] for k in names], bulk, _coh(), crack_tracker=tracker)
    p1, p2 = mk(1, ("a", "b", "r")), mk(2, ("a", "c2", "b"))
    for el in (p1, p2, mk(3, ("a", "r", "c")), mk(4, ("a", "c", "c2"))):
        dom.add_element(el)
    dom.add_element(Tri3(10, [n["r"], n["e1"], n["c"]], Elastic2D(E=E, nu=0.2), thickness=1.0))
    dom.add_element(Tri3(11, [n["c"], n["e2"], n["c2"]], Elastic2D(E=E, nu=0.2), thickness=1.0))
    dom.generate_equation_numbers(verbose=False)
    p1.activate(np.array([0.5, -0.3]), 0, np.array([0.0, 0.5]))
    p2.activate(np.array([-0.5, -0.8]), 0, np.array([0.0, 0.5]))
    tracker.update(np.zeros(dom.total_dofs))
    tracker._refresh(0)
    return tracker, n, {k: np.array(x) for k, x in X.items()}


def test_tercer_nodo_compartido_por_dos_puntas():
    """acceptance.regla_de_la_punta (tercer nodo de dos puntas): la punta del
    primer agrietado da a c un lado y la del segundo el otro; vale la primera
    según el orden de activación (§2)."""
    tracker, n, X = _fan()
    f = tracker.cracks[0].fixed
    va, vb, vr, vc2 = f[1], f[2], f[3], f[4]

    def tip_value(P_in, Q_other, partner_far, v_far):
        d = (P_in - Q_other) / np.linalg.norm(P_in - Q_other)
        sc, sa = _x2(d, X["c"] - P_in), _x2(d, X["a"] - P_in)
        return np.sign(va if (sc >= 0) == (sa >= 0) else v_far) * abs(sc)
    Pab = X["a"] + va / (va - vb) * (X["b"] - X["a"])
    Par = X["a"] + va / (va - vr) * (X["r"] - X["a"])
    Pac2 = X["a"] + va / (va - vc2) * (X["c2"] - X["a"])
    v1 = tip_value(Par, Pab, "r", vr)          # punta q1, del agrietado p1 (primero)
    v2 = tip_value(Pac2, Pab, "c2", vc2)       # punta q2, del agrietado p2
    assert np.sign(v1) != np.sign(v2)          # el caso discrimina el orden
    assert tracker._tips[0] == {tracker._pos[3], tracker._pos[4]}
    assert tracker._values[0][tracker._nidx[n["c"].id]] == pytest.approx(v1, abs=1e-12)


def test_el_cambio_de_lado_no_sobrevive_a_la_evaluacion():
    """§5: la regla del vértice actúa también en una punta bajo σ_t0 (aquí
    cambia c con tracción a 21°, mientras sólo un elemento lejano supera su
    resistencia), pero el cambio se descarta al terminar la evaluación: en la
    siguiente, con el campo a 1°, la punta se activa igual que con una sola
    evaluación."""
    def extra(dom, n, tracker):
        far = [dom.add_node(20 + i, x) for i, x in enumerate(([10.0, 0.0], [11.0, 0.0], [10.0, 1.0]))]
        dom.add_element(CST_Embedded2D(5, far, Elastic2D(E=E, nu=0.2, hypothesis="plane_stress"),
                                       _coh(0.1 * ST0), crack_tracker=tracker))
    tr0, dom0, _ = _tip(0.3, extra)
    tr0.update(_uniform_field(dom0, 1.0))
    ref = dom0.elements[2].discontinuity_state
    tracker, dom, _ = _tip(0.3, extra)
    tracker.update(_uniform_field(dom, 21.0, eps=0.5 * ST0 * (1 - 0.2 ** 2) / E))
    assert dom.elements[2].discontinuity_state is None
    assert dom.elements[5].discontinuity_state is not None
    assert tracker._flips[0]                     # c cambió de lado en esa evaluación
    assert 4 not in tracker.cracks[0].fixed
    tracker.update(_uniform_field(dom, 1.0))
    ds = dom.elements[2].discontinuity_state
    assert ds is not None and ds.solitary_node == ref.solitary_node
    assert tracker.cracks[0].fixed[4] == tr0.cracks[0].fixed[4]



@pytest.mark.parametrize("seed,plates", [(482, False), (2413, True)])
def test_regla_del_vertice_una_vez_y_solo_si_queda_cos_positivo(seed, plates):
    """Abanicos al azar de la verificación del teórico (c6.py): dos agrietados
    alrededor de a y dos puntas q1, q2 que comparten c. En los dos casos q1 se
    activa y q2 queda con un aviso. Detecta que c cambie de lado más de una
    vez por evaluación (semilla 482: q2 se activaría y q1 no) y que cambie
    aunque no quede cos ψ > 0 (semilla 2413: ninguna se activaría, 2 avisos)."""
    rng = np.random.default_rng(seed)
    j = lambda x, y: (x + rng.uniform(-0.35, 0.35), y + rng.uniform(-0.35, 0.35))
    X = dict(a=(0.0, 0.0), b=j(0.0, 1.0), r=j(-1.0, 0.0), c2=j(1.0, 0.0), c=j(-0.3, -1.0),
             e1=(-2.5, -1.5), e2=(2.5, -1.5))
    dom = Domain()
    n = {k: dom.add_node(i, list(x)) for i, (k, x) in enumerate(X.items(), 1)}
    tracker = CrackTracking()
    bulk = Elastic2D(E=E, nu=0.2, hypothesis="plane_stress")
    mk = lambda i, names: CST_Embedded2D(i, [n[k] for k in names], bulk, _coh(), crack_tracker=tracker)
    p1, p2, q1, q2 = mk(1, ("a", "b", "r")), mk(2, ("a", "c2", "b")), mk(3, ("a", "r", "c")), mk(4, ("a", "c", "c2"))
    for el in (p1, p2, q1, q2):
        dom.add_element(el)
    if plates:
        dom.add_element(Tri3(10, [n["r"], n["e1"], n["c"]], Elastic2D(E=E, nu=0.2), thickness=1.0))
        dom.add_element(Tri3(11, [n["c"], n["e2"], n["c2"]], Elastic2D(E=E, nu=0.2), thickness=1.0))
    dom.generate_equation_numbers(verbose=False)
    xa, xb, xr, xc2 = (np.array(X[k]) for k in ("a", "b", "r", "c2"))
    w = rng.uniform(0.15, 0.85, 3)
    Pab, Par, Pac2 = xa + w[0] * (xb - xa), xa + w[1] * (xr - xa), xa + w[2] * (xc2 - xa)

    def normal(P, Q):
        t_ = Q - P
        m = np.array([-t_[1], t_[0]])
        return m if m @ (xa - P) > 0 else -m
    p1.activate(normal(Pab, Par), 0, Pab)
    p2.activate(normal(Pab, Pac2), 0, Pab)
    phi = rng.uniform(0, np.pi)
    tracker.update(_uniform_field(dom, np.degrees(phi)))
    assert q1.discontinuity_state is not None and q1.discontinuity_state.solitary_node == 0
    assert q2.discontinuity_state is None
    assert tracker.unactivated_warnings == 1
    assert tracker.side_conflicts() == []

# --- Borde libre (spec §3) -------------------------------------------------
# En la punta q = (a, c, b) los lados de salida son a–c y b–c. Con c = (2, 0.3)
# la prolongación lleva c con a (+0.15).

def _plate_on_bc(dom, n, tracker):
    """Tri3 pegado a b–c: b–c deja de ser libre y la grieta sale por a–c."""
    e = dom.add_node(6, [2.0, 1.5])
    dom.add_element(Tri3(10, [n["c"], e, n["b"]], Elastic2D(E=E, nu=0.2), thickness=1.0))


def _quad_on_ac(dom, n, cls=Quad4, **kw):
    g, h = dom.add_node(7, [1.0, -1.0]), dom.add_node(8, [2.0, -1.0])
    if cls is Quad4Thermal:
        el = cls(11, [n["a"], g, h, n["c"]], ThermalConduction(k=1.0, c=1.0, density=1.0, dim=2))
    else:
        el = cls(11, [n["a"], g, h, n["c"]], Elastic2D(E=E, nu=0.2), **kw)
    dom.add_element(el)


def test_sin_lados_libres_distintos_decide_la_prolongacion():
    """Esquina: a–c y b–c libres los dos; decide la prolongación."""
    _, _, v = _tip(0.3)
    assert v == pytest.approx(0.15, abs=1e-12)


def test_sale_por_el_unico_lado_libre():
    _, _, v = _tip(0.3, _plate_on_bc)
    assert v == pytest.approx(-0.15, abs=1e-12)         # c con b: sale por a–c


def test_una_placa_pegada_tapa_el_lado():
    def extra(dom, n, tracker):
        _plate_on_bc(dom, n, tracker)
        _quad_on_ac(dom, n)
    _, _, v = _tip(0.3, extra)
    assert v == pytest.approx(0.15, abs=1e-12)          # ningún lado libre: la prolongación


def test_otro_seguimiento_vecino_tapa_el_lado():
    """El lado común con un CST_Embedded2D de otro seguimiento no es libre
    (CST_Embedded2D declara EDGE_NODES)."""
    def extra(dom, n, tracker):
        _plate_on_bc(dom, n, tracker)
        g = dom.add_node(7, [1.0, -1.0])
        dom.add_element(CST_Embedded2D(12, [n["a"], g, n["c"]], Elastic2D(E=E, nu=0.2, hypothesis="plane_stress"),
                                       _coh(), crack_tracker=CrackTracking()))
    _, _, v = _tip(0.3, extra)
    assert v == pytest.approx(0.15, abs=1e-12)


def test_los_elementos_termicos_no_tapan_el_lado_libre():
    def extra(dom, n, tracker):
        _plate_on_bc(dom, n, tracker)
        _quad_on_ac(dom, n, cls=Quad4Thermal)
    _, _, v = _tip(0.3, extra)
    assert v == pytest.approx(-0.15, abs=1e-12)


def test_una_placa_anadida_despues_cambia_la_decision():
    """El borde libre se recalcula cuando cambian los elementos del dominio
    (Domain.topology_version), no sólo al dar de alta un miembro."""
    tracker, dom, v = _tip(0.3, _plate_on_bc)
    assert v == pytest.approx(-0.15, abs=1e-12)
    n = {k: dom.nodes[i] for k, i in (("a", 1), ("c", 4))}
    _quad_on_ac(dom, n)
    tracker._refresh(0)
    assert tracker._values[0][tracker._nidx[4]] == pytest.approx(0.15, abs=1e-12)


def test_sin_dominio_error_claro_solo_si_hace_falta():
    """Elementos que no están en ningún dominio: el seguimiento sólo necesita
    el dominio cuando una punta toca un lado sin vecino en su región."""
    from solidum.core.node import Node
    nodes = [Node(i, x) for i, x in
             enumerate((list(_A), list(_B), list(_R), [2.0, 0.3]), 1)]
    tracker = CrackTracking()
    bulk = Elastic2D(E=E, nu=0.2, hypothesis="plane_stress")
    p = CST_Embedded2D(1, [nodes[0], nodes[1], nodes[2]], bulk, _coh(), crack_tracker=tracker)
    CST_Embedded2D(2, [nodes[0], nodes[3], nodes[1]], bulk, _coh(), crack_tracker=tracker)
    for i, nd in enumerate(nodes):
        nd.dofs = {"ux": 2 * i, "uy": 2 * i + 1}
    p.activate(np.array([0.0, -1.0]), 0, np.array([0.5, 0.45]))
    tracker.update(np.zeros(8))                         # sin puntas que decidir: sin error
    with pytest.raises(RuntimeError, match="no pertenece a ningún dominio"):
        tracker._refresh(0)


# =========================================================================
# Sólo el seguimiento con un campo impuesto
# =========================================================================


def _region(pts, tris, sigma_t0, tracker):
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(pts)]
    bulk = Elastic2D(E=E, nu=0.0, hypothesis="plane_stress")
    for e, t in enumerate(tris):
        dom.add_element(CST_Embedded2D(e + 1, [nodes[k] for k in t], bulk,
                                       _coh(float(sigma_t0(pts[t].mean(axis=0)))), thickness=0.1,
                                       crack_tracker=tracker))
    dom.generate_equation_numbers(verbose=False)
    return dom, nodes


def _crossings(tracker, k=0):
    """Cruces de la grieta ``k`` con los lados de sus elementos (de los
    valores nodales), en orden de activación."""
    fixed = tracker._fixed_idx[k]
    P = []
    for eid in tracker.cracks[k].elements:
        tri = tracker._conn[tracker._pos[eid]]
        th = np.array([fixed[i] for i in tri])
        s = th >= 0.0
        for a, b in _EDGES:
            if s[a] != s[b]:
                xa, xb = tracker._X[tri[a]], tracker._X[tri[b]]
                P.append(xa + th[a] / (th[a] - th[b]) * (xb - xa))
    return np.array(P)


class _Recorder:
    """Anota cada ``CST_Embedded2D.activate`` como el registro de referencia
    del prototipo: [elemento, nodo solitario, grieta, n_x, n_y, x, y]."""

    def __init__(self, monkeypatch):
        self.log = []
        orig = embedded_cst.CST_Embedded2D.activate

        def activate(el, normal, solitary_node, point, crack=None):
            orig(el, normal, solitary_node, point, crack)
            n = np.asarray(normal, dtype=float)
            self.log.append([int(el.id), int(solitary_node), None if crack is None else int(crack),
                             *(n / np.linalg.norm(n)), *np.asarray(point, dtype=float)])
        monkeypatch.setattr(embedded_cst.CST_Embedded2D, "activate", activate)


_C_ARC, _R_ARC = np.array([-1.0, 0.5]), 1.5


def _arc(n, seed):
    """Campo radial u = c·r²·e_r alrededor de (−1, 0.5): σ_I radial, las
    trayectorias exactas son arcos; raíz en (0.5, 0.5), arco de radio 1.5.
    Una sola evaluación activa toda la grieta."""
    pts, tris = _delaunay(1.0, 1.0, n, n, seed)
    root = np.array([0.5, 0.5])
    c = 1.0e-4
    st0 = lambda x: 0.5 * E * 2.0 * c * np.hypot(*(x - _C_ARC)) * (1.0 + np.hypot(*(x - root)))
    tracker = CrackTracking(max_cracks=1)
    dom, nodes = _region(pts, tris, st0, tracker)
    U = np.zeros(dom.total_dofs)
    for nd, x in zip(nodes, pts):
        r = x - _C_ARC
        U[nd.dofs["ux"]], U[nd.dofs["uy"]] = c * np.hypot(*r) * r
    tracker.update(U)
    return tracker, 1.0 / n


# Secuencia de activaciones del prototipo aprobado (regla_continuacion.py,
# recta + borde libre, sobre 4f1b58a) en las 15 mallas de arcos; origen en el
# propio archivo. Las decisiones (elemento, nodo solitario, grieta) se exigen
# exactas; normales y puntos, a 1e-11: el redondeo de σ_I depende del núcleo
# de BLAS de la máquina (medido: hasta 3.3e-13 relativo con otro núcleo de
# OpenBLAS, revisión de arquitectura del 2026-11-19).
_REFERENCE = json.loads((Path(__file__).resolve().parent / "data" / "referencia_prototipo.json")
                        .read_text(encoding="utf-8"))
_ARC_REFERENCE = {tuple(int(v) for v in k.split("_")): rows for k, rows in _REFERENCE["arcos"].items()}


def _assert_same_sequence(log, ref):
    """Decisiones (elemento, nodo solitario, grieta) exactas; normales y
    puntos a 1e-11."""
    assert [r[:3] for r in log] == [r[:3] for r in ref], (
        "Las activaciones ya no son las de la referencia (spec CrackTracking, "
        "equivalencia_con_el_prototipo).")
    np.testing.assert_allclose(np.array([r[3:] for r in log]), np.array([r[3:] for r in ref]),
                               rtol=0, atol=1e-11)
_ARC_INCOMPLETE = {(32, 1), (32, 4), (64, 2)}


@pytest.mark.parametrize("n,seed", sorted(_ARC_REFERENCE), ids=[f"n{n}_s{s}" for n, s in sorted(_ARC_REFERENCE)])
def test_arcos_iguales_al_prototipo_aprobado(n, seed, monkeypatch):
    """acceptance.equivalencia_con_el_prototipo y trayectoria_en_arco: la
    secuencia de activaciones es la del prototipo aprobado (decisiones exactas,
    normales y puntos a 1e-11); las mismas 12 de 15 completas (limitación
    conocida, spec §6); en las completas, desviación ≤ 0.8·h."""
    rec = _Recorder(monkeypatch)
    tracker, h = _arc(n, seed)
    _assert_same_sequence(rec.log, _ARC_REFERENCE[(n, seed)])
    P = _crossings(tracker)
    complete = P[:, 1].min() < 1e-9 and P[:, 1].max() > 1.0 - 1e-9
    assert complete == ((n, seed) not in _ARC_INCOMPLETE)
    if complete:
        assert np.max(np.abs(np.hypot(*(P - _C_ARC).T) - _R_ARC)) <= 0.8 * h
        assert tracker.side_conflicts() == []


def _chord(point, t):
    s = []
    for ax in (0, 1):
        for val in (0.0, 1.0):
            if abs(t[ax]) > 1e-14:
                si = (val - point[ax]) / t[ax]
                x = point + si * t
                if -1e-9 <= x[1 - ax] <= 1 + 1e-9:
                    s.append(si)
    return max(s) - min(s)


def _straight(n, seed, phi_deg):
    """Recta en campo uniforme: tracción uniaxial según φ, raíz en el centro
    (revisión del teórico del 2026-11-19, e6.py). Devuelve (fracción de la
    cuerda recorrida, desviación máxima / h)."""
    pts, tris = _delaunay(1.0, 1.0, n, n, seed)
    phi = np.radians(phi_deg)
    nrm = np.array([np.cos(phi), np.sin(phi)])
    t = np.array([-nrm[1], nrm[0]])
    eps, root = 1.0e-4, np.array([0.5, 0.5])
    tracker = CrackTracking(max_cracks=1)
    dom, nodes = _region(pts, tris, lambda x: 0.5 * E * eps * (1.0 + 0.5 * np.hypot(*(x - root))), tracker)
    U = np.zeros(dom.total_dofs)
    for nd, x in zip(nodes, pts):
        U[nd.dofs["ux"]], U[nd.dofs["uy"]] = eps * nrm * ((x - root) @ nrm)
    tracker.update(U)
    P = _crossings(tracker)
    x0 = np.asarray(dom.elements[tracker.cracks[0].root].discontinuity_state.centroid, dtype=float)
    along = (P - x0) @ t
    return (along.max() - along.min()) / _chord(x0, t), float(np.abs((P - x0) @ nrm).max() * n)


# Rectas que no completan ≥ 98 % de la cuerda o se desvían más de h con el
# prototipo aprobado: (φ, semilla). Las de 45° con n = 16 salen al borde un
# poco antes (2–5 %), sin daño; las demás son la limitación M0 (spec §6).
_STRAIGHT_FAIL = {
    16: {(10.0, 15), (20.0, 5), (60.0, 10)} | {(45.0, s) for s in (2, 3, 4, 5, 6, 8, 9, 13, 14, 15, 17, 18)},
    32: {(0.0, 1), (0.0, 4), (10.0, 19), (30.0, 1), (45.0, 0), (75.0, 4), (75.0, 10), (75.0, 12), (75.0, 14)},
}


@pytest.mark.parametrize("n", [16, 32])
def test_recta_en_campo_uniforme_fija_la_limitacion_conocida(n):
    """acceptance.recta_campo_uniforme: 7 ángulos × 20 semillas; fallan
    exactamente las mismas que con el prototipo aprobado (limitación M0)."""
    fails = set()
    for phi in (0.0, 10.0, 20.0, 30.0, 45.0, 60.0, 75.0):
        for seed in range(20):
            frac, dev = _straight(n, seed, phi)
            if frac < 0.98 or dev > 1.0:
                fails.add((phi, seed))
    assert fails == _STRAIGHT_FAIL[n]


def test_la_parte_formada_no_cambia_y_la_nueva_sigue_la_direccion():
    """Reescritura de «rotación de direcciones»: media grieta formada con
    tracción según x (las resistencias impiden pasar de y = 0.5); con el campo
    girado 30° y mayor, los valores guardados no cambian, la grieta sigue desde
    su punta y, delante, sigue la dirección nueva (dx/dy = −tan 30°)."""
    for seed in (3, 5, 7):
        pts, tris = _delaunay(1.0, 1.0, 32, 32, seed)
        eps, h = 1.0e-4, 1.0 / 32
        st0 = lambda c: E * eps * (0.6 + 2.0 * abs(c[0] - 0.5) + (10.0 if c[1] > 0.5 else 0.0)
                                   + 0.1 * abs(c[1] - 0.25))
        tracker = CrackTracking(max_cracks=1)
        dom, nodes = _region(pts, tris, st0, tracker)

        def field(nrm, amp, x0):
            U = np.zeros(dom.total_dofs)
            for nd, x in zip(nodes, pts):
                U[nd.dofs["ux"]], U[nd.dofs["uy"]] = amp * nrm * ((x - x0) @ nrm)
            return U
        tracker.update(field(np.array([1.0, 0.0]), eps, np.array([0.5, 0.25])))
        fixed1, els1 = dict(tracker.cracks[0].fixed), list(tracker.cracks[0].elements)
        P1 = _crossings(tracker)
        assert P1[:, 1].min() < 1e-9 and 0.45 < P1[:, 1].max() < 0.55
        ang = np.radians(30.0)
        tracker.update(field(np.array([np.cos(ang), np.sin(ang)]), 20 * eps, P1[np.argmax(P1[:, 1])]))
        assert all(tracker.cracks[0].fixed[k] == v for k, v in fixed1.items())
        assert tracker.cracks[0].elements[:len(els1)] == els1
        P2 = _crossings(tracker)
        assert P2[:, 1].max() > 1.0 - 1e-9
        ahead = P2[len(P1):]
        ahead = ahead[ahead[:, 1] > P1[:, 1].max() + 3 * h]
        assert np.polyfit(ahead[:, 1], ahead[:, 0], 1)[0] == pytest.approx(-np.tan(ang), abs=0.02)
        assert tracker.side_conflicts() == [] and tracker.unactivated_warnings == 0


def test_con_ruido_en_el_campo_la_limitacion_conocida():
    """Reescritura de «con ruido la trayectoria sigue continua»: tracción
    uniaxial con ruido en los desplazamientos (desviación media de las
    direcciones ≈ 5°), todos los elementos por encima de σ_t0. Sin conflictos
    ni avisos; cruza de borde a borde en 7 de 10 mallas, y en las otras 3 la
    grieta se encuentra consigo misma (limitación de §6, como M1)."""
    complete = []
    for seed in range(10):
        pts, tris = _delaunay(1.0, 1.0, 32, 32, seed)
        eps, root = 1.0e-4, np.array([0.5, 0.5])
        tracker = CrackTracking(max_cracks=1)
        dom, nodes = _region(pts, tris, lambda c: 0.25 * E * eps * (1.0 + np.hypot(*(c - root))), tracker)
        rng = np.random.default_rng(seed)
        U = np.zeros(dom.total_dofs)
        for nd, x in zip(nodes, pts):
            U[nd.dofs["ux"]], U[nd.dofs["uy"]] = (eps * np.array([1.0, 0.0]) * (x[0] - 0.5)
                                                  + 0.1 * eps * rng.normal(0, 1, 2) / 32)
        tracker.update(U)
        assert tracker.side_conflicts() == [] and tracker.unactivated_warnings == 0
        P = _crossings(tracker)
        complete.append(bool(P[:, 1].min() < 1e-9 and P[:, 1].max() > 1.0 - 1e-9))
    assert [s for s, ok in enumerate(complete) if not ok] == [2, 3, 8]


class TestGeometria:

    def test_particion_de_la_raiz_siempre_con_cos_positivo(self):
        rng = np.random.default_rng(7)
        used_shift = 0
        for _ in range(2000):
            xs = rng.uniform(-1, 1, (3, 2))
            if abs((xs[1, 0] - xs[0, 0]) * (xs[2, 1] - xs[0, 1])
                   - (xs[2, 0] - xs[0, 0]) * (xs[1, 1] - xs[0, 1])) < 1e-3:
                continue
            a = rng.uniform(0, np.pi)
            n = np.array([np.cos(a), np.sin(a)])
            lone, nn, point = root_partition(xs, n)
            grads, _ = triangle_gradients(xs, np.array([[0, 1, 2]]))
            inward = grads[0] / np.linalg.norm(grads[0], axis=1, keepdims=True)
            assert nn @ inward[lone] > 0.0
            assert abs(abs(nn @ n) - 1.0) < 1e-14
            side = (xs - point) @ nn
            assert side[lone] > 0.0 and np.all(np.delete(side, lone) < 0.0)
            if not np.allclose(point, xs.mean(axis=0)):
                used_shift += 1
                # recta desplazada: a media distancia (proyección sobre n)
                # entre el nuevo solitario y el nodo intermedio (§3)
                s = xs @ nn
                mid = int(np.argsort(s)[1])
                assert (point - 0.5 * (xs[lone] + xs[mid])) @ nn == pytest.approx(0.0, abs=1e-12)
        assert used_shift > 0          # el caso de la otra partición aparece

    def test_particion_por_signo(self):
        xs = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        grads, _ = triangle_gradients(xs, np.array([[0, 1, 2]]))
        inward = grads[0] / np.linalg.norm(grads[0], axis=1, keepdims=True)
        lone, n, x_m, cos_psi = partition(xs, np.array([-0.5, 0.5, -0.5]), np.array([-1.0, 0.0]), inward)
        assert lone == 1 and n[0] > 0.0
        np.testing.assert_allclose(x_m, [0.5, 0.25])
        assert cos_psi == pytest.approx(1.0)
        assert partition(xs, np.array([1.0, 2.0, 3.0]), np.array([1.0, 0.0]), inward) is None

    def test_distancia_con_signo_a_la_poligonal(self):
        P = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 1.0]])
        X = np.array([[0.5, 0.3], [0.5, -0.2], [-1.0, 0.4], [3.0, 2.0 + 1e-9], [1.0, -0.1]])
        d = polyline_signed_distance(P, X)
        np.testing.assert_allclose(d[:3], [0.3, -0.2, 0.4])
        assert abs(d[3]) < 1e-8                        # prolongación del último tramo
        assert d[4] == pytest.approx(-0.1)


# =========================================================================
# Con el elemento real
# =========================================================================

L, W = 0.1, 0.05


def _strip(pts, tris, *, sigma_t0, embed, tracker=None, delta=1.6e-4, nsteps=400):
    """Franja L×W en tracción (u_x = δ a la derecha, u_x = 0 a la izquierda,
    u_y = 0 abajo). Devuelve (elementos embebidos, seguimiento, curva F-δ,
    fallos del Newton local)."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(pts)]
    bulk = Elastic2D(E=E, nu=NU, hypothesis="plane_stress")
    embedded = []
    for e, t in enumerate(tris):
        nds = [nodes[k] for k in t]
        c = pts[t].mean(axis=0)
        if embed(t, c):
            el = CST_Embedded2D(e + 1, nds, bulk, _coh(sigma_t0(c)), thickness=TH, crack_tracker=tracker)
            embedded.append(el)
        else:
            el = Tri3(e + 1, nds, bulk, thickness=TH)
        dom.add_element(el)
    right = []
    for i, p in enumerate(pts):
        if abs(p[0]) < 1e-12:
            nodes[i].fix_dof("ux", 0.0)
        if abs(p[0] - L) < 1e-12:
            nodes[i].fix_dof("ux", delta)
            right.append(nodes[i])
        if abs(p[1]) < 1e-12:
            nodes[i].fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    right = [nd.dofs["ux"] for nd in right]
    asm = Assembler(dom)
    curve = [(0.0, 0.0)]

    def record(step, U, lam):
        _, F = asm.assemble_non_linear_system(U)
        curve.append((lam * delta, float(F[right].sum())))

    level = logging.getLogger("solidum").level
    logging.getLogger("solidum").setLevel(logging.ERROR)
    try:
        NonlinearSolver(asm, num_steps=nsteps, adaptive=False, max_iter=40).solve(
            np.zeros(dom.total_dofs), step_callback=record)
    finally:
        logging.getLogger("solidum").setLevel(level)
    tracker = tracker or embedded[0].crack_tracker
    return embedded, tracker, np.array(curve), sum(e.local_newton_failures for e in embedded)


def _work(curve):
    d, F = curve[:, 0], curve[:, 1]
    return float(np.sum(0.5 * (F[1:] + F[:-1]) * np.diff(d)))


def _weak_plane(x0):
    """σ_t0 mínima en (x0, 0), creciente con |x − x0| y, poco, con y."""
    return lambda c: ST0 * (0.85 + 0.3 * c[1] / W + 3.0 * abs(c[0] - x0) / L)


def _banda(x0, half):
    return lambda t, c: abs(c[0] - x0) <= half


_FRANJAS = {
    "estructurada_16x8": lambda: _structured(L, W, 16, 8),
    "sesgada_16x8": lambda: _structured(L, W, 16, 8, skew=0.6),
    "delaunay_24x12_s1": lambda: _delaunay(L, W, 24, 12, 1),
    "delaunay_24x12_s2": lambda: _delaunay(L, W, 24, 12, 2),
}


class TestFranjasNoAlineadas:
    """acceptance.franjas_no_alineadas."""

    @pytest.mark.parametrize("malla", list(_FRANJAS))
    def test_una_grieta_coherente_que_separa_la_franja(self, malla):
        pts, tris = _FRANJAS[malla]()
        x0 = 0.0537
        els, tracker, curve, fails = _strip(pts, tris, sigma_t0=_weak_plane(x0), embed=_banda(x0, 0.02))
        F = curve[:, 1]
        assert tracker.side_conflicts() == []
        assert tracker.side_conflict_warnings == 0
        assert fails == 0
        assert abs(F[-1]) < 1.0e-3 * F.max()
        main = max(tracker.cracks, key=lambda c: len(c.elements))
        by_id = {e.id: e for e in els}
        l_d = [by_id[i].discontinuity_state.l_d for i in main.elements]
        assert min(l_d) > 0.0
        assert 0.99 <= sum(l_d) / W <= 1.012
        # Energía: G_F·W·t salvo el efecto de las normales de Rankine en Σl_d
        # (hasta ~1 % en mallas Delaunay; spec, Diálogo del 2026-11-10).
        assert 0.995 <= _work(curve) / (GF * W * TH) <= 1.012
        for e in els:
            ds = e.discontinuity_state
            if ds is not None and e.id in main.elements:
                assert ds.cohesive_state_committed.get("damage", 0.0) > 0.999


def _astillas():
    """Las 5 primeras mallas Delaunay 16×8 (jitter 0.3) en que la recta
    x = x_c deja algún elemento con cos ψ < −0.01 (escaneo del 2026-11-10)."""
    found = []
    n_x = np.array([1.0, 0.0])
    for seed in range(2, 400):
        pts, tris = _delaunay(L, W, 16, 8, seed)
        grads, _ = triangle_gradients(pts, tris)
        inward = grads / np.linalg.norm(grads, axis=2, keepdims=True)
        for xc in (0.0437, 0.0537, 0.0611):
            th = pts[:, 0] - xc
            cos_min = min((part[3] for e, t in enumerate(tris)
                           if (part := partition(pts[t], th[t], n_x, inward[e])) is not None), default=1.0)
            if cos_min < -0.01:
                found.append((seed, xc))
        if len(found) >= 5:
            return found[:5]
    raise AssertionError("no se encontraron mallas con astillas")


_ASTILLAS = _astillas()


class TestAstillas:
    """acceptance.astillas: regla del vértice con trayectoria prescrita."""

    @pytest.mark.parametrize("seed,xc", _ASTILLAS, ids=[f"s{s}_x{x}" for s, x in _ASTILLAS])
    def test_separa_con_la_regla_del_vertice(self, seed, xc, monkeypatch):
        """Además, la secuencia de activaciones es la del código anterior al
        2026-11-19 (``PrescribedCrackPath`` no cambió de comportamiento)."""
        rec = _Recorder(monkeypatch)
        pts, tris = _delaunay(L, W, 16, 8, seed)
        cut_nodes = {int(i) for t in tris if (pts[t, 0] > xc).any() and (pts[t, 0] <= xc).any()
                     for i in t}
        tracker = PrescribedCrackPath([[xc, 0.0], [xc, W]])
        els, tracker, curve, fails = _strip(
            pts, tris, sigma_t0=_weak_plane(xc), tracker=tracker, nsteps=200,
            embed=lambda t, c: bool(cut_nodes.intersection(int(i) for i in t)))
        F = curve[:, 1]
        assert fails == 0
        assert tracker.unactivated_warnings == 0
        assert tracker.side_conflicts() == []
        assert abs(F[-1]) < 1.0e-3 * F.max()
        assert all(e.discontinuity_state.l_d > 0.0 for e in els if e.discontinuity_state is not None)
        assert 0.999 <= _work(curve) / (GF * W * TH) <= 1.012
        _assert_same_sequence(rec.log, _REFERENCE["astillas"][f"{seed}_{xc}"])


def _two_planes_sigma_t0(planes):
    def sigma_t0(c):
        d = min(abs(c[0] - x0) for x0 in planes)
        return ST0 * (0.85 + 2.0 * c[1] / W + 3.0 * d / L)
    return sigma_t0


class TestDosGrietas:
    """acceptance.dos_grietas"""

    def test_dos_raices_distintas_sin_conflictos(self):
        pts, tris = _structured(L, W, 20, 8)
        planes = (0.03, 0.07)
        els, tracker, curve, fails = _strip(
            pts, tris, sigma_t0=_two_planes_sigma_t0(planes), tracker=CrackTracking(),
            delta=1.1e-5, nsteps=44, embed=lambda t, c: min(abs(c[0] - x0) for x0 in planes) <= 0.012)
        assert len(tracker.cracks) == 2
        roots = [np.mean([n.coordinates[0] for n in next(e for e in els if e.id == c.root).nodes])
                 for c in tracker.cracks]
        assert sorted(round(r, 2) for r in roots) == [0.03, 0.07]
        assert tracker.side_conflicts() == []
        assert fails == 0

    def test_una_sola_grieta_con_max_cracks_1(self):
        """acceptance.una_sola_grieta: el mismo modelo con max_cracks = 1 da
        una grieta con raíz en uno de los dos planos; ningún elemento del otro
        se agrieta, aunque sin el límite nace allí la segunda."""
        pts, tris = _structured(L, W, 20, 8)
        planes = (0.03, 0.07)
        els, tracker, curve, fails = _strip(
            pts, tris, sigma_t0=_two_planes_sigma_t0(planes), tracker=CrackTracking(max_cracks=1),
            delta=1.1e-5, nsteps=44, embed=lambda t, c: min(abs(c[0] - x0) for x0 in planes) <= 0.012)
        assert len(tracker.cracks) == 1
        by_id = {e.id: e for e in els}
        x_root = np.mean([n.coordinates[0] for n in by_id[tracker.cracks[0].root].nodes])
        assert round(x_root, 2) in planes
        other = planes[1] if round(x_root, 2) == planes[0] else planes[0]
        cerca_del_otro = [e for e in els
                          if abs(np.mean([n.coordinates[0] for n in e.nodes]) - other) <= 0.012]
        assert cerca_del_otro and all(e.discontinuity_state is None for e in cerca_del_otro)
        assert fails == 0


# =========================================================================
# Un triángulo aislado hasta la separación (aceptación de CST_Embedded2D)
# =========================================================================

# Triángulo, ψ de la partición de la raíz con n = x, y tolerancia de la
# energía disipada frente a G_F·l_d·t (None: no se separa).
_TRIANGULOS = {
    "rectangulo": (np.array([[0.0, 0.0], [0.1, 0.0], [0.0, 0.08]]), 0.0, 1.0e-3),
    "agudo": (np.array([[0.0, 0.0], [0.1, 0.02], [0.04, 0.09]]), 24.0, 3.0e-3),
    "obtuso": (np.array([[0.0, 0.0], [0.1, 0.03], [0.1542, 0.0414]]), 73.0, None),
}


@pytest.mark.parametrize("nombre", list(_TRIANGULOS))
def test_tres_triangulos_traccion_hasta_separacion(nombre):
    """acceptance.tres_triangulos_traccion_hasta_separacion (CST_Embedded2D):
    ε_xx impuesta creciente (ε_yy = γ_xy = 0). La partición de la raíz da
    cos ψ > 0 también en el obtuso, donde la recta por el centroide aislaría
    un nodo con cos ψ < 0 (con |cos| daba t_n hasta −400·σ_t0): ya no hay
    compresión espuria. Con ψ grande la formulación KOS con ley de modo I
    bloquea (el teórico, 2026-11-10: esfuerzo residual creciente con ψ): el
    obtuso (ψ = 73°) no se separa y el agudo (ψ = 24°) disipa 0.2 % de más.
    Es un elemento aislado, que no puede redistribuir; en una malla no aparece
    y la formulación se mantiene (decisión 4 del usuario, 2026-11-10)."""
    tri, psi_deg, tol = _TRIANGULOS[nombre]
    coords = tri[_ccw(tri, np.array([[0, 1, 2]]))[0]]
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
    el = CST_Embedded2D(1, nodes, Elastic2D(E=E, nu=0.2, hypothesis="plane_stress"), _coh(), thickness=TH)
    dom.add_element(el)
    dom.generate_equation_numbers(verbose=False)
    x = coords[:, 0] - coords[:, 0].min()
    B, detJ = compute_kinematics_tri3(coords)
    work, prev_power, prev_eps, tn_min = 0.0, 0.0, 0.0, np.inf
    e0 = ST0 / E
    steps = np.concatenate([np.linspace(0.0, 1.5 * e0, 3000), np.linspace(1.5 * e0, 30 * e0, 3000)[1:]])
    for eps in steps[1:]:
        u = np.zeros(6)
        u[0::2] = eps * x
        el.prepare_step(u)
        _, F = el.compute_element_state(u)
        el.commit_state()
        power = float(F[0::2] @ x)
        work += 0.5 * (power + prev_power) * (eps - prev_eps)
        prev_power, prev_eps = power, eps
        ds = el.discontinuity_state
        if ds is not None:
            t, _, _ = el.cohesive_material.compute_traction(ds.jump_committed, ds.cohesive_state_committed or None)
            tn_min = min(tn_min, t[0])
    ds = el.discontinuity_state
    assert ds is not None and ds.l_d > 0.0
    others = [i for i in range(3) if i != ds.solitary_node]
    edge = coords[others[1]] - coords[others[0]]
    cos_psi = 2.0 * ds.l_d / np.hypot(*edge)
    assert np.degrees(np.arccos(cos_psi)) == pytest.approx(psi_deg, abs=0.5)
    assert el.local_newton_failures == 0
    assert tn_min >= -1.0e-6 * ST0
    t_end, _, _ = el.cohesive_material.compute_traction(ds.jump_committed, ds.cohesive_state_committed or None)
    if tol is None:
        assert t_end[0] > 0.5 * ST0            # bloqueo KOS documentado
        return
    assert t_end[0] < 1.0e-3 * ST0
    sigma = el.state.stresses[0]
    strain_bulk = B @ u - B[:, 2 * ds.solitary_node:2 * ds.solitary_node + 2] @ (
        np.column_stack([ds.normal, ds.tangent]) @ ds.jump_committed)
    elastic = 0.5 * float(sigma @ strain_bulk) * 0.5 * detJ * TH
    assert (work - elastic) / (GF * ds.l_d * TH) == pytest.approx(1.0, abs=tol)


# =========================================================================
# Plumbing: grupos por omisión, ronda por paso, esfuerzos vectorizados, YAML
# =========================================================================


def _tri_nodes(dom, coords, start):
    return [dom.add_node(start + i, list(c)) for i, c in enumerate(coords)]


def _emb(eid, nodes, **kw):
    return CST_Embedded2D(eid, nodes, Elastic2D(E=E, nu=0.2, hypothesis="plane_stress"), _coh(), **kw)


def test_seguimiento_por_omision_por_region_conexa():
    dom = Domain()
    a = _tri_nodes(dom, [[0, 0], [1, 0], [0, 1]], 1)
    b = _tri_nodes(dom, [[5, 0], [6, 0], [5, 1]], 4)
    e1 = _emb(1, a)
    e2 = _emb(2, [a[1], dom.add_node(7, [1, 1]), a[2]])
    e3 = _emb(3, b)
    t1 = default_tracker(e1)
    assert isinstance(t1, CrackTracking)
    assert e2.crack_tracker is t1
    assert e3.crack_tracker is not None and e3.crack_tracker is not t1
    e4 = _emb(4, [a[1], b[0], b[2]])          # une las dos regiones (aún sin grietas)
    default_tracker(e4)
    assert e1.crack_tracker is e2.crack_tracker is e3.crack_tracker is e4.crack_tracker


def test_seguimiento_explicito_compartido():
    dom = Domain()
    a = _tri_nodes(dom, [[0, 0], [1, 0], [0, 1]], 1)
    tracker = CrackTracking()
    e1 = _emb(1, a, crack_tracker=tracker)
    e2 = _emb(2, [a[1], dom.add_node(4, [1, 1]), a[2]], crack_tracker=tracker)
    assert tracker.elements == (e1, e2)
    with pytest.raises(ValueError):
        PrescribedCrackPath([[0.0, 0.0]])


def test_sustituir_un_elemento_de_la_region():
    """Un ``CST_Embedded2D`` sustituido con ``Domain.replace_element`` deja
    de ser miembro del seguimiento: queda el nuevo, sin duplicar su id."""
    dom = Domain()
    a = _tri_nodes(dom, [[0, 0], [1, 0], [0, 1]], 1)
    tracker = CrackTracking()
    e1 = _emb(1, a, crack_tracker=tracker)
    e2 = _emb(2, [a[1], dom.add_node(4, [1, 1]), a[2]], crack_tracker=tracker)
    dom.add_element(e1)
    dom.add_element(e2)
    new = _emb(2, list(e2.nodes), crack_tracker=tracker)
    dom.replace_element(new)
    assert tracker.elements == (e1, new)
    tracker._build()
    assert len(tracker._conn) == 2


@pytest.mark.parametrize("valor", [0, -1, 1.5, True, "1"])
def test_max_cracks_invalido(valor):
    """acceptance.max_cracks_invalido"""
    with pytest.raises(ValueError, match="max_cracks"):
        CrackTracking(max_cracks=valor)


def _grid_region(tracker_kw=None):
    pts, tris = _structured(L, W, 4, 2)
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(pts)]
    tracker = CrackTracking(**(tracker_kw or {}))
    els = [_emb(e + 1, [nodes[k] for k in t], crack_tracker=tracker) for e, t in enumerate(tris)]
    for el in els:
        dom.add_element(el)
    dom.generate_equation_numbers(verbose=False)
    U = np.zeros(dom.total_dofs)
    for nd in nodes:
        U[nd.dofs["ux"]] = 2.0 * ST0 / E * nd.coordinates[0]
    return dom, tracker, els, U


def test_revise_step_avisa_solo_si_aumentan_los_agrietados():
    """La señal que repite el paso (ADR 0041) es el aumento del número de
    agrietados, no cualquier cambio del seguimiento; con el mismo campo, la
    segunda ronda no activa nada y devuelve False (la cascada termina)."""
    dom, tracker, els, U = _grid_region()
    assert not any(el.revise_step(np.zeros(dom.total_dofs)) for el in els)   # nada supera σ_t0
    senales = [el.revise_step(U) for el in els]
    assert senales[0] is True and not any(senales[1:])   # una evaluación por ronda
    assert any(e.discontinuity_state is not None for e in els)
    assert not any(el.revise_step(U.copy()) for el in els)   # mismo campo: nada nuevo


def test_una_actualizacion_por_ronda_e_idempotente():
    dom, tracker, els, U = _grid_region()
    calls = []
    original = tracker.update
    tracker.update = lambda U_: (calls.append(1), original(U_))[1]
    for el in els:
        el.prepare_step(U)
    assert len(calls) == 1
    cracked = [e.id for e in els if e.discontinuity_state is not None]
    assert cracked
    for el in els:                             # otra ronda con el mismo campo
        el.prepare_step(U.copy())
    assert len(calls) == 2
    assert [e.id for e in els if e.discontinuity_state is not None] == cracked
    assert tracker.side_conflicts() == []


def test_el_dominio_con_seguimiento_se_serializa():
    """Tras una evaluación que calculó el borde libre (caché con referencia
    débil al dominio), el dominio se puede serializar y el seguimiento de la
    copia sigue funcionando con el dominio de la copia."""
    dom, tracker, els, U = _grid_region()
    tracker.update(U)
    assert tracker._free_cache is not None
    back = pickle.loads(pickle.dumps(dom))
    tr2 = back.elements[1].crack_tracker
    assert tr2 is not tracker and tr2._free_cache is None
    assert [c.elements for c in tr2.cracks] == [c.elements for c in tracker.cracks]
    tr2._free_edges(0)                       # el borde libre de la copia, con su dominio
    assert tr2._free_cache[0]() is back


def test_esfuerzos_vectorizados_igual_que_el_elemento():
    pts, tris = _delaunay(L, W, 6, 3, 5)
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(pts)]
    tracker = CrackTracking()
    els = [_emb(e + 1, [nodes[k] for k in t], crack_tracker=tracker) for e, t in enumerate(tris)]
    for el in els:
        dom.add_element(el)
    dom.generate_equation_numbers(verbose=False)
    U = np.random.default_rng(0).normal(0.0, 1e-5, dom.total_dofs)
    tracker._build()
    sigma_I, normals = tracker._principal_stresses(U, np.ones(len(els), dtype=bool))
    for p, el in enumerate(tracker._els):
        s_ref, n_ref = el.principal_stress(U)
        assert sigma_I[p] == pytest.approx(s_ref, rel=1e-12)
        assert abs(normals[p] @ n_ref) == pytest.approx(1.0, abs=1e-12)


_YAML = """\
    user_modules: [discontinuities]
    nodes:
      - {id: 1, coords: [0.0, 0.0]}
      - {id: 2, coords: [1.0, 0.0]}
      - {id: 3, coords: [0.0, 1.0]}
      - {id: 4, coords: [1.0, 1.0]}
    materials:
      - {id: 1, type: Elastic2D, E: 30.0e9, nu: 0.2, hypothesis: plane_stress}
    interface_materials:
      - {id: 1, type: CohesiveDamageIsotropic, sigma_t0: 3.0e6, G_f: 100.0, K_e: 1.0e15, softening: linear}
    crack_trackers:
      - {id: 1, %s}
      - {id: 2, type: PrescribedCrackPath, points: [[0.5, 0.0], [0.5, 1.0]]}
    elements:
      - {id: 1, type: CST_Embedded2D, nodes: [1, 2, 3], material: 1, cohesive_material: 1, crack_tracker: 1}
      - {id: 2, type: CST_Embedded2D, nodes: [2, 4, 3], material: 1, cohesive_material: 1, crack_tracker: 2}
    """


def _parse(tmp_path, tracker_entry):
    from solidum.utils.yaml_parser import YamlParser
    ruta = tmp_path / "modelo.yaml"
    ruta.write_text(textwrap.dedent(_YAML % tracker_entry), encoding="utf-8")
    return YamlParser(str(ruta)).parse()


def test_yaml_crack_trackers(tmp_path):
    dom = _parse(tmp_path, "type: CrackTracking, max_cracks: 1")
    e1, e2 = dom.elements[1], dom.elements[2]
    assert isinstance(e1.crack_tracker, CrackTracking) and e1.crack_tracker.max_cracks == 1
    assert isinstance(e2.crack_tracker, PrescribedCrackPath)
    np.testing.assert_allclose(e2.crack_tracker.points, [[0.5, 0.0], [0.5, 1.0]])


@pytest.mark.parametrize("entrada,mensaje", [
    ("type: GlobalCrackTracking, max_cracks: 1", r"desconocido 'GlobalCrackTracking'.*CrackTracking"),
    ("type: CrackTracking, epsilon: 1.0e-7", r"'epsilon' no aceptado"),
])
def test_yaml_nombres_retirados(tmp_path, entrada, mensaje):
    """acceptance.nombres_retirados: el seguimiento global y su parámetro
    epsilon ya no existen; el error lo dice (con la lista de tipos o el
    parámetro desconocido)."""
    from solidum.utils.yaml_parser import YamlValidationError
    with pytest.raises(YamlValidationError, match=mensaje):
        _parse(tmp_path, entrada)
