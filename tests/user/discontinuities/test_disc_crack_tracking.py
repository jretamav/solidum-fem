"""Seguimiento de grieta del módulo de discontinuidades.

Specs ``GlobalCrackTracking`` y ``PrescribedCrackPath`` (aceptación), y la
aceptación ``tres_triangulos_traccion_hasta_separacion`` de ``CST_Embedded2D``.

Dos bloques:

- **Sólo el seguimiento** (funciones puras del problema auxiliar): trayectoria
  con un campo de direcciones conocido (arcos de círculo), sensibilidad a la
  regularización ``ε``, continuidad con ruido en las direcciones y giro de las
  direcciones con media grieta formada.
- **Con el elemento real** y ``NonlinearSolver``: franjas en tracción con
  mallas no alineadas (seguimiento global), mallas con astillas
  (``cos ψ < 0``, trayectoria prescrita y regla del vértice), dos grietas y
  un triángulo aislado hasta la separación.
"""
from __future__ import annotations

import logging
import textwrap

import numpy as np
import pytest
from scipy.spatial import Delaunay

import solidum
from solidum import Assembler, Domain
from solidum.elements.solid_2d import compute_kinematics_tri3
from solidum.elements.solid_2d.tri3 import Tri3
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.solvers import NonlinearSolver

solidum.load_user_module("discontinuities")
from solidum.user.discontinuities import (  # noqa: E402
    CST_Embedded2D,
    CohesiveDamageIsotropic,
    GlobalCrackTracking,
    PrescribedCrackPath,
)
from solidum.user.discontinuities.crack_tracking import (  # noqa: E402
    default_tracker,
    isoline_segments,
    partition,
    polyline_signed_distance,
    root_partition,
    solve_tracking_field,
    tracking_matrix,
    triangle_gradients,
)

_EDGES = ((0, 1), (1, 2), (2, 0))


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


def _edge_map(tris):
    edges = {}
    for e, tri in enumerate(tris):
        for a, b in _EDGES:
            edges.setdefault(tuple(sorted((tri[a], tri[b]))), []).append(e)
    return edges


def _locate(pts, tris, x):
    for e, t in enumerate(tris):
        a, b, c = pts[t]
        lam = np.linalg.solve(np.column_stack([b - a, c - a]), x - a)
        if lam.min() >= -1e-12 and lam.sum() <= 1 + 1e-12:
            return e
    raise ValueError(x)


# =========================================================================
# Sólo el seguimiento: problema auxiliar con campo de direcciones prescrito
# =========================================================================

_C = np.array([-1.0, 0.5])      # centro de los arcos
_R0 = 1.5                       # el arco pasa por (0.5, 0.5)


def _arc_tangents(pts, tris, noise_deg=0.0, seed=0):
    r = pts[tris].mean(axis=1) - _C
    t = np.column_stack([-r[:, 1], r[:, 0]])
    t /= np.linalg.norm(t, axis=1, keepdims=True)
    if noise_deg:
        ang = np.radians(np.random.default_rng(seed).normal(0.0, noise_deg, len(t)))
        c, s = np.cos(ang), np.sin(ang)
        t = np.column_stack([c * t[:, 0] - s * t[:, 1], s * t[:, 0] + c * t[:, 1]])
    return t


def _isoline_from(pts, tris, theta, roots):
    """Segmentos de θ = 0 de la componente conexa (por lados cruzados) que
    contiene los elementos ``roots``."""
    edges = _edge_map(tris)
    s = theta >= 0.0
    seen, stack = set(roots), list(roots)
    while stack:
        e = stack.pop()
        tri = tris[e]
        for a, b in _EDGES:
            if s[tri[a]] != s[tri[b]]:
                for q in edges[tuple(sorted((tri[a], tri[b])))]:
                    if q not in seen:
                        seen.add(q)
                        stack.append(q)
    return isoline_segments(pts, tris, theta, sorted(seen))


def _arc_path(n, *, eps=1.0e-6, noise_deg=0.0, seed=0):
    pts, tris = _delaunay(1.0, 1.0, n, n, seed=0)
    root = _locate(pts, tris, np.array([0.5, 0.5]))
    grads, area = triangle_gradients(pts, tris)
    K = tracking_matrix(grads, area, _arc_tangents(pts, tris, noise_deg, seed), len(pts), tris, eps)
    fixed = {int(i): float(np.hypot(*(pts[i] - _C)) - _R0) for i in tris[root]}
    theta = solve_tracking_field(K, np.arange(len(pts)), fixed)
    return pts, tris, 1.0 / n, _isoline_from(pts, tris, theta, [root])


def _deviation(segs):
    P = np.array([p for _, a, b in segs for p in (a, b)])
    return float(np.max(np.abs(np.hypot(*(P - _C).T) - _R0)))


def _distance_between(segs_a, segs_b):
    """Máxima distancia de los puntos de ``segs_a`` a la poligonal ``segs_b``."""
    A = np.array([p for _, a, b in segs_a for p in (a, b)])
    worst = 0.0
    for x in A:
        best = np.inf
        for _, p, q in segs_b:
            d = q - p
            w = np.clip((x - p) @ d / max(d @ d, 1e-300), 0.0, 1.0)
            best = min(best, float(np.hypot(*(x - p - w * d))))
        worst = max(worst, best)
    return worst


class TestCampoConocido:
    """acceptance.trayectoria_con_campo_conocido"""

    def test_la_isolinea_converge_al_arco(self):
        devs = []
        for n in (16, 32, 64):
            _, _, h, segs = _arc_path(n)
            devs.append(_deviation(segs))
            assert devs[-1] <= 0.02 * h, (n, devs[-1] / h)
        assert devs[0] > devs[1] > devs[2], devs

    def test_con_ruido_la_trayectoria_sigue_continua(self):
        """Con 5° de ruido por elemento la isolínea que nace de la raíz cruza
        el dominio de borde a borde (no se corta). Su desviación del arco no
        es menor que la del seguimiento local (datos del 2026-11-10 en la
        spec): el ruido equivale a una difusión transversal."""
        for seed in range(3):
            pts, tris, h, segs = _arc_path(32, noise_deg=5.0, seed=seed)
            P = np.array([p for _, a, b in segs for p in (a, b)])
            assert P[:, 1].min() < 1e-9 and P[:, 1].max() > 1.0 - 1e-9
            assert _deviation(segs) < 0.1


class TestSensibilidadEpsilon:
    """acceptance.sensibilidad_a_epsilon"""

    @pytest.mark.parametrize("n", [16, 32, 64])
    def test_la_trayectoria_apenas_depende_de_epsilon_pequeno(self, n):
        _, _, h, segs_a = _arc_path(n, eps=1.0e-6)
        _, _, _, segs_b = _arc_path(n, eps=1.0e-8)
        assert _distance_between(segs_a, segs_b) < 0.1 * h

    def test_epsilon_grande_desvia_la_trayectoria(self):
        """El regularizador difunde a través de las líneas del campo una
        distancia ~√ε·L, independiente de h: con ε = 1e-2 la isolínea se
        aparta del arco más de un elemento."""
        _, _, h, segs = _arc_path(64, eps=1.0e-2)
        assert _deviation(segs) > h


class TestGiroDeDirecciones:
    """acceptance.rotacion_de_direcciones"""

    def test_la_parte_formada_no_cambia_y_la_nueva_sigue_la_direccion(self):
        pts, tris = _delaunay(1.0, 1.0, 32, 32, seed=3)
        h = 1.0 / 32
        grads, area = triangle_gradients(pts, tris)
        cen = pts[tris].mean(axis=1)
        y0 = 0.5
        formed = [e for e, t in enumerate(tris)
                  if cen[e, 0] < 0.5 and (pts[t, 1] > y0).any() and (pts[t, 1] <= y0).any()]
        fixed = {}
        for e in formed:
            for i in tris[e]:
                d = pts[i, 1] - y0
                fixed.setdefault(int(i), np.sign(d) * max(abs(d), 1e-3 * h) if d != 0 else 1e-3 * h)
        ang = np.radians(30.0)
        t = np.tile([np.cos(ang), np.sin(ang)], (len(tris), 1))
        t[formed] = (1.0, 0.0)
        K = tracking_matrix(grads, area, t, len(pts), tris, 1.0e-6)
        theta = solve_tracking_field(K, np.arange(len(pts)), fixed)
        # La parte formada conserva su partición.
        for i, v in fixed.items():
            assert np.sign(theta[i]) == np.sign(v)
        # Por delante de la punta, la isolínea sigue la dirección nueva.
        segs = _isoline_from(pts, tris, theta, formed)
        P = np.array([p for _, a, b in segs for p in (a, b)])
        tip_x = max(p[0] for e, a, b in segs if e in formed for p in (a, b))
        ahead = P[P[:, 0] > tip_x + 2 * h]
        assert len(ahead) > 10
        assert np.max(np.abs(ahead[:, 1] - (y0 + np.tan(ang) * (ahead[:, 0] - tip_x)))) < 1.5 * h


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
            used_shift += not np.allclose(point, xs.mean(axis=0))
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

E, NU, ST0, GF, TH = 30.0e9, 0.0, 3.0e6, 100.0, 0.02
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
            coh = CohesiveDamageIsotropic(sigma_t0=sigma_t0(c), G_f=GF, K_e=1.0e15, softening="linear")
            el = CST_Embedded2D(e + 1, nds, bulk, coh, thickness=TH, crack_tracker=tracker)
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
    """acceptance.franjas_no_alineadas (seguimiento global)."""

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
    def test_separa_con_la_regla_del_vertice(self, seed, xc):
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


class TestDosGrietas:
    """acceptance.dos_grietas"""

    def test_dos_raices_distintas_sin_conflictos(self):
        pts, tris = _structured(L, W, 20, 8)
        planes = (0.03, 0.07)

        def sigma_t0(c):
            d = min(abs(c[0] - x0) for x0 in planes)
            return ST0 * (0.85 + 2.0 * c[1] / W + 3.0 * d / L)
        tracker = GlobalCrackTracking()
        els, tracker, curve, fails = _strip(
            pts, tris, sigma_t0=sigma_t0, tracker=tracker, delta=1.1e-5, nsteps=44,
            embed=lambda t, c: min(abs(c[0] - x0) for x0 in planes) <= 0.012)
        assert len(tracker.cracks) == 2
        roots = [np.mean([n.coordinates[0] for n in next(e for e in els if e.id == c.root).nodes])
                 for c in tracker.cracks]
        assert sorted(round(r, 2) for r in roots) == [0.03, 0.07]
        assert tracker.side_conflicts() == []
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
    Es la decisión pendiente del modo mixto; la prueba lo deja medido."""
    tri, psi_deg, tol = _TRIANGULOS[nombre]
    coords = tri[_ccw(tri, np.array([[0, 1, 2]]))[0]]
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
    el = CST_Embedded2D(1, nodes, Elastic2D(E=E, nu=0.2, hypothesis="plane_stress"),
                        CohesiveDamageIsotropic(sigma_t0=ST0, G_f=GF, K_e=1.0e15, softening="linear"),
                        thickness=TH)
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
    return CST_Embedded2D(eid, nodes, Elastic2D(E=E, nu=0.2, hypothesis="plane_stress"),
                          CohesiveDamageIsotropic(sigma_t0=ST0, G_f=GF, K_e=1.0e15, softening="linear"), **kw)


def test_seguimiento_por_omision_por_region_conexa():
    dom = Domain()
    a = _tri_nodes(dom, [[0, 0], [1, 0], [0, 1]], 1)
    b = _tri_nodes(dom, [[5, 0], [6, 0], [5, 1]], 4)
    e1 = _emb(1, a)
    e2 = _emb(2, [a[1], dom.add_node(7, [1, 1]), a[2]])
    e3 = _emb(3, b)
    t1 = default_tracker(e1)
    assert e2.crack_tracker is t1
    assert e3.crack_tracker is not None and e3.crack_tracker is not t1
    e4 = _emb(4, [a[1], b[0], b[2]])          # une las dos regiones (aún sin grietas)
    default_tracker(e4)
    assert e1.crack_tracker is e2.crack_tracker is e3.crack_tracker is e4.crack_tracker


def test_seguimiento_explicito_compartido():
    dom = Domain()
    a = _tri_nodes(dom, [[0, 0], [1, 0], [0, 1]], 1)
    tracker = GlobalCrackTracking(epsilon=1.0e-5)
    e1 = _emb(1, a, crack_tracker=tracker)
    e2 = _emb(2, [a[1], dom.add_node(4, [1, 1]), a[2]], crack_tracker=tracker)
    assert tracker.elements == (e1, e2)
    with pytest.raises(ValueError):
        GlobalCrackTracking(epsilon=0.0)
    with pytest.raises(ValueError):
        PrescribedCrackPath([[0.0, 0.0]])


def test_una_actualizacion_por_ronda_e_idempotente():
    pts, tris = _structured(L, W, 4, 2)
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(pts)]
    tracker = GlobalCrackTracking()
    els = [_emb(e + 1, [nodes[k] for k in t], crack_tracker=tracker) for e, t in enumerate(tris)]
    for el in els:
        dom.add_element(el)
    dom.generate_equation_numbers(verbose=False)
    U = np.zeros(dom.total_dofs)
    for nd in nodes:
        U[nd.dofs["ux"]] = 2.0 * ST0 / E * nd.coordinates[0]
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


def test_esfuerzos_vectorizados_igual_que_el_elemento():
    pts, tris = _delaunay(L, W, 6, 3, 5)
    dom = Domain()
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(pts)]
    tracker = GlobalCrackTracking()
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


def test_yaml_crack_trackers(tmp_path):
    from solidum.utils.yaml_parser import YamlParser
    texto = textwrap.dedent("""\
        user_modules: [discontinuities]
        nodes:
          - {id: 1, coords: [0.0, 0.0]}
          - {id: 2, coords: [1.0, 0.0]}
          - {id: 3, coords: [0.0, 1.0]}
          - {id: 4, coords: [1.0, 1.0]}
        materials:
          - {id: 1, type: Elastic2D, E: 30.0e9, nu: 0.2, hypothesis: plane_stress}
        cohesive_materials:
          - {id: 1, type: CohesiveDamageIsotropic, sigma_t0: 3.0e6, G_f: 100.0, K_e: 1.0e15, softening: linear}
        crack_trackers:
          - {id: 1, type: GlobalCrackTracking, epsilon: 1.0e-7}
          - {id: 2, type: PrescribedCrackPath, points: [[0.5, 0.0], [0.5, 1.0]]}
        elements:
          - {id: 1, type: CST_Embedded2D, nodes: [1, 2, 3], material: 1, cohesive_material: 1, crack_tracker: 1}
          - {id: 2, type: CST_Embedded2D, nodes: [2, 4, 3], material: 1, cohesive_material: 1, crack_tracker: 2}
        """)
    ruta = tmp_path / "modelo.yaml"
    ruta.write_text(texto, encoding="utf-8")
    dom = YamlParser(str(ruta)).parse()
    e1, e2 = dom.elements[1], dom.elements[2]
    assert isinstance(e1.crack_tracker, GlobalCrackTracking) and e1.crack_tracker.epsilon == 1.0e-7
    assert isinstance(e2.crack_tracker, PrescribedCrackPath)
    np.testing.assert_allclose(e2.crack_tracker.points, [[0.5, 0.0], [0.5, 1.0]])
