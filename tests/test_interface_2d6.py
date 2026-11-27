"""Aceptación de ``Interface2D6``, interfaz de espesor nulo 2D con caras
cuadráticas (spec ``docs/specs/Interface2D6.md``; revisión del teórico en
``informes/teorico/revision_spec_interface2d6_2026-11-26.md``).

Invariantes:
- Parcela de contacto exacta con Quad8, Quad9 y Tri6, nodal y Gauss, con una
  ley unilateral real (el signo de la normal se ve), materiales distintos y
  espesor ≠ 1; con ``split`` desde el YAML y con ``pair`` en Python.
- En un lado recto, las fuerzas de una tracción uniforme son e·L·(1/6, 1/6,
  2/3)·t con las dos reglas (detecta pesos de Gauss desordenados).
- En un lado curvo: ejes punto a punto (no los de la cuerda), pesos nodales
  de Simpson con el jacobiano del nodo (decisión Q1 del usuario) y fuerzas
  e·t·e_z × (a(−1)/3, a(+1)/3, 4a(0)/3); giro rígido sin salto ni fuerza;
  rango 6 de la rigidez con las dos reglas.
- Cara degenerada (tangente que se anula o se invierte, nodo a un cuarto):
  error al construir.
- Las pruebas de la familia: nodal desacopla pares, escalado con el espesor,
  simetría, orientación invertida en cara curva, cohesiva con holgura
  rechazada, holgura inicial, tangente por diferencias.
"""
from __future__ import annotations

import textwrap

import numpy as np
import pytest

import solidum
from solidum.core.node import Node
from solidum.elements.interface import Interface2D6
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.interface import CohesiveDamageIsotropic, ContactPenalty
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver
from solidum.mesh import pair_interfaces, structured_mesh

K_N = 4.0e5
K_T = 2.5e5
MU = 0.35
E_T = 0.37                                  # espesor ≠ 1 (Reglas §6)

RECTA = np.array([[0.2, 0.1], [1.9, 0.9], [1.05, 0.5]])         # (1, 2, 5), inclinada
CURVA = np.array([[0.0, 0.0], [2.0, 0.2], [1.0, 0.35]])          # intermedio fuera de la cuerda


# Convenio de la spec escrito a mano (no se lee de la clase, para que un cambio
# del orden en la clase lo detecten las pruebas): nodos 1, 2 = extremos de A;
# 3, 4 = extremos de B con 4 frente a 1 y 3 frente a 2; 5 = intermedio de A;
# 6 = intermedio de B. Posiciones (base 0) de A en el orden (1, 2, 5) y de su
# par en B.
POS_A = (0, 1, 4)
POS_B = (3, 2, 5)


def iface(XA, XB=None, law=None, **kw):
    """Interface2D6 con la cara A = (x_1, x_2, x_5) y la B frente a ella."""
    XB = XA if XB is None else XB
    coords = [XA[0], XA[1], XB[1], XB[0], XA[2], XB[2]]
    nodes = [Node(k + 1, list(map(float, c))) for k, c in enumerate(coords)]
    kw.setdefault("thickness", E_T)
    return Interface2D6(1, nodes, law or ContactPenalty(K_n=K_N, mu=MU, K_t=K_T), **kw)


def u_caras(du_por_par):
    """u_e con la cara A fija y cada nodo de B desplazado ``du_por_par[i]``
    (en el orden de A: 1, 2, 5)."""
    u = np.zeros(12)
    for i, b in enumerate(POS_B):
        u[2 * b:2 * b + 2] = du_por_par[i]
    return u


def tangente_de_la_cara(X, xi):
    """a(ξ) = Σ N_i'(ξ)·x_i con N = (ξ(ξ−1)/2, ξ(ξ+1)/2, 1 − ξ²)."""
    dN = np.array([xi - 0.5, xi + 0.5, -2.0 * xi])
    return dN @ X


def fd_stiffness(e, u, h=1e-9):
    snap = e.state.snapshot_trial()
    K = np.zeros((u.size, u.size))
    for j in range(u.size):
        d = np.zeros(u.size)
        d[j] = h
        _, Fp = e.compute_element_state(u + d)
        _, Fm = e.compute_element_state(u - d)
        K[:, j] = (Fp - Fm) / (2 * h)
    e.state.restore_trial(snap)
    return K


# ---------------------------------------------------------------------------
# Lado recto
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("integration", ["nodal", "gauss"])
def test_salto_uniforme_reparte_un_sexto_un_sexto_dos_tercios(integration):
    """Tracción uniforme con K_n ≠ K_t (adherencia, cara cerrada): fuerzas
    e·L·(1/6, 1/6, 2/3)·Rᵀt en B y lo opuesto en A."""
    e = iface(RECTA, integration=integration)
    n, s = e.R[0]
    g = np.array([-2e-4, 1e-5])
    _, F = e.compute_element_state(u_caras([g[0] * n + g[1] * s] * 3))
    t = np.array([K_N * g[0], K_T * g[1]])
    L = float(np.linalg.norm(RECTA[1] - RECTA[0]))
    f = t[0] * n + t[1] * s
    for frac, (a, b) in zip((1 / 6, 1 / 6, 2 / 3), zip(POS_A, POS_B)):
        np.testing.assert_allclose(F[2 * b:2 * b + 2], E_T * L * frac * f, rtol=1e-12)
        np.testing.assert_allclose(F[2 * a:2 * a + 2], -E_T * L * frac * f, rtol=1e-12)


def test_holgura_inicial_en_los_tres_puntos():
    e0 = iface(RECTA)
    n = e0.R[0, 0]
    g0 = 0.013
    e = iface(RECTA, RECTA + g0 * n, law=ContactPenalty(K_n=K_N))
    np.testing.assert_allclose(e.g0, g0, rtol=1e-12)
    _, F = e.compute_element_state(u_caras([-0.5 * g0 * n] * 3))
    np.testing.assert_allclose(F, 0.0, atol=1e-12)                  # aún abierta
    _, F = e.compute_element_state(u_caras([-1.3 * g0 * n] * 3))
    assert np.linalg.norm(F) > 0.0


def test_cohesiva_con_holgura_rechazada():
    n = iface(RECTA).R[0, 0]
    ley = CohesiveDamageIsotropic(sigma_t0=3.0, G_f=0.1, K_e=1e4, softening="linear")
    with pytest.raises(ValueError, match="no admite holgura"):
        iface(RECTA, RECTA + 0.01 * n, law=ley)


# ---------------------------------------------------------------------------
# Lado curvo
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("integration", ["nodal", "gauss"])
def test_ejes_punto_a_punto_en_cara_curva(integration):
    e = iface(CURVA, integration=integration)
    xis = ([-1.0, 1.0, 0.0] if integration == "nodal"
           else list(np.polynomial.legendre.leggauss(3)[0]))
    for g, xi in enumerate(xis):
        a = tangente_de_la_cara(CURVA, xi)
        s = a / np.linalg.norm(a)
        np.testing.assert_allclose(e.R[g], [[-s[1], s[0]], s], rtol=1e-12, atol=1e-14)


def test_pesos_de_simpson_con_el_jacobiano_del_nodo_en_cara_curva():
    """Decisión Q1 del usuario: m_i = e·w_i·J(ξ_i), w = (1/3, 1/3, 4/3); con
    un salto normal local uniforme, las fuerzas de B son
    e·t_n·e_z × (a(−1)/3, a(+1)/3, 4a(0)/3)."""
    e = iface(CURVA, law=ContactPenalty(K_n=K_N))
    a = [tangente_de_la_cara(CURVA, xi) for xi in (-1.0, 1.0, 0.0)]
    w = (1 / 3, 1 / 3, 4 / 3)
    np.testing.assert_allclose(e.measures, [E_T * wi * np.linalg.norm(ai) for wi, ai in zip(w, a)],
                               rtol=1e-13)
    delta = -3e-4
    _, F = e.compute_element_state(u_caras([delta * e.R[g, 0] for g in range(3)]))
    t_n = K_N * delta
    for wi, ai, b in zip(w, a, POS_B):
        ez_x_a = np.array([-ai[1], ai[0]])
        np.testing.assert_allclose(F[2 * b:2 * b + 2], E_T * t_n * wi * ez_x_a, rtol=1e-12)


@pytest.mark.parametrize("integration", ["nodal", "gauss"])
def test_giro_rigido_en_cara_curva(integration):
    e = iface(CURVA, integration=integration)
    X = np.array([n.coordinates for n in e.nodes], float)
    omega, c = 2e-3, np.array([0.013, -0.021])
    u = (omega * np.column_stack([-X[:, 1], X[:, 0]]) + c).ravel()
    np.testing.assert_array_equal(e.jumps(u), 0.0)
    _, F = e.compute_element_state(u)
    np.testing.assert_allclose(F, 0.0, atol=1e-12)


@pytest.mark.parametrize("integration", ["nodal", "gauss"])
def test_rango_de_la_rigidez(integration):
    """Ley lineal (adherencia, cara cerrada): 6 valores propios no nulos; con
    una regla de 2 puntos serían 4."""
    # μ alto: todos los puntos en adherencia (en deslizamiento la tangente de
    # la ley es singular en s y el rango bajaría por la ley, no por la regla).
    e = iface(CURVA, integration=integration, law=ContactPenalty(K_n=K_N, mu=5.0, K_t=K_T))
    n = e.R[0, 0]
    K, _ = e.compute_element_state(u_caras([-1e-4 * n] * 3))
    assert all(v["status"] == 1.0 for v in e.state.vars_trial)
    ev = np.abs(np.linalg.eigvals(K))
    assert np.sum(ev > 1e-8 * ev.max()) == 6


@pytest.mark.parametrize("case", ["adherencia", "deslizamiento"])
@pytest.mark.parametrize("integration", ["nodal", "gauss"])
def test_tangente_por_diferencias_en_cara_curva(case, integration):
    e = iface(CURVA, integration=integration)
    rng = np.random.default_rng(7)
    n, s = e.R[2]
    shift = -1e-3 * n + (1e-6 if case == "adherencia" else 3e-3) * s
    u = u_caras([shift] * 3) + 1e-8 * rng.standard_normal(12)
    K, _ = e.compute_element_state(u)
    np.testing.assert_allclose(K, fd_stiffness(e, u), rtol=1e-6, atol=1e-6 * K_N)


@pytest.mark.parametrize("XA,motivo", [
    (np.array([[0.0, 0.0], [1.0, 0.0], [0.25, 0.0]]), "nodo a un cuarto"),
    (np.array([[0.0, 0.0], [1.0, 0.0], [-0.2, 0.0]]), "intermedio fuera del lado"),
    (np.array([[0.0, 0.0], [1.0, 0.0], [0.125, 0.05]]), "peso consistente negativo"),
])
def test_cara_degenerada_es_error(XA, motivo):
    with pytest.raises(ValueError, match="degenerada"):
        iface(XA)


# ---------------------------------------------------------------------------
# Pruebas de la familia
# ---------------------------------------------------------------------------

def test_nodal_desacopla_pares_y_gauss_no():
    u = u_caras([np.array([0.0, -1e-3])] * 3)
    pares = [[2 * a, 2 * a + 1, 2 * b, 2 * b + 1] for a, b in zip(POS_A, POS_B)]
    K, _ = iface(RECTA, law=ContactPenalty(K_n=K_N)).compute_element_state(u)
    Kg, _ = iface(RECTA, law=ContactPenalty(K_n=K_N), integration="gauss").compute_element_state(u)
    for i in range(3):
        for j in range(3):
            if i != j:
                assert np.all(K[np.ix_(pares[i], pares[j])] == 0.0)
    assert np.any(Kg[np.ix_(pares[0], pares[2])] != 0.0)


def test_escalado_con_el_espesor():
    u = u_caras([np.array([1e-4, -1e-3])] * 3)
    K1, F1 = iface(CURVA, thickness=1.0).compute_element_state(u)
    K2, F2 = iface(CURVA, thickness=E_T).compute_element_state(u)
    np.testing.assert_allclose(K2, E_T * K1, rtol=1e-12)
    np.testing.assert_allclose(F2, E_T * F1, rtol=1e-12)


def test_simetria_sin_friccion():
    e = iface(CURVA, law=ContactPenalty(K_n=K_N), integration="gauss")
    K, _ = e.compute_element_state(u_caras([np.array([2e-4, -1e-3])] * 3))
    np.testing.assert_allclose(K, K.T, rtol=0, atol=1e-10 * np.abs(K).max())


@pytest.mark.parametrize("integration", ["nodal", "gauss"])
def test_orientacion_invertida_en_cara_curva(integration):
    """Intercambiar A y B (nodos coincidentes) da las mismas fuerzas sobre
    cada nodo físico."""
    law = ContactPenalty(K_n=K_N, mu=MU, K_t=K_T)
    e1 = iface(CURVA, law=law, integration=integration)
    n, s = e1.R[2]
    dB = -1e-3 * n + 2e-3 * s
    _, F1 = e1.compute_element_state(u_caras([dB] * 3))
    # A' = cara B recorrida al revés (2, 1, 5), B' = cara A.
    XA2 = CURVA[[1, 0, 2]]
    e2 = iface(XA2, law=law, integration=integration)
    u2 = np.zeros(12)
    for a in POS_A:
        u2[2 * a:2 * a + 2] = dB
    _, F2 = e2.compute_element_state(u2)
    # Nodo físico → posición en cada elemento: A_k de e1 ↔ B' de e2, etc.
    perm_a1 = {0: 1, 1: 0, 2: 2}      # el i-ésimo de A en e1 es el perm[i]-ésimo de A' en e2
    for i in range(3):
        a1, b1 = POS_A[i], POS_B[i]
        j = perm_a1[i]
        a2, b2 = POS_A[j], POS_B[j]
        np.testing.assert_allclose(F2[2 * b2:2 * b2 + 2], F1[2 * a1:2 * a1 + 2], rtol=1e-11,
                                   atol=1e-9)
        np.testing.assert_allclose(F2[2 * a2:2 * a2 + 2], F1[2 * b1:2 * b1 + 2], rtol=1e-11,
                                   atol=1e-9)


# ---------------------------------------------------------------------------
# Parcela de contacto (Quad8, Quad9, Tri6)
# ---------------------------------------------------------------------------

P = 35.0
E1, NU1, H1 = 2.1e4, 0.2, 0.7
E2, NU2, H2 = 0.9e4, 0.3, 1.3
KN_PARCELA = 1.0e6                    # ≈ 50·E/h


def _u_tope():
    return -P * ((1 - NU1 ** 2) * H1 / E1 + (1 - NU2 ** 2) * H2 / E2) - P / KN_PARCELA


@pytest.mark.parametrize("element", ["Quad8", "Quad9", "Tri6"])
@pytest.mark.parametrize("integration", ["nodal", "gauss"])
def test_parcela_de_contacto_desde_el_yaml(tmp_path, element, integration):
    """Dos bloques de materiales distintos unidos por ``split`` en y = H1,
    presión p arriba, deformación plana: t_n = −p, σ_yy = −p y el
    desplazamiento de la cara superior exactos (*calculado*: resortes en
    serie)."""
    ruta = tmp_path / "m.yaml"
    ruta.write_text(textwrap.dedent(f"""\
        materials:
          - {{id: 1, type: Elastic2D, E: {E1}, nu: {NU1}, hypothesis: plane_strain}}
          - {{id: 2, type: Elastic2D, E: {E2}, nu: {NU2}, hypothesis: plane_strain}}
        interface_materials:
          - {{id: 1, type: ContactPenalty, K_n: {KN_PARCELA}}}
        mesh_element: {element}
        mesh_thickness: {E_T}
        mesh:
          structured:
            blocks:
              - {{name: abajo, material: 1, x: {{length: 2, n: 4}}, y: {{length: {H1}, n: 2}}}}
              - {{name: arriba, material: 2, origin: [0, {H1}], x: {{length: 2, n: 4}},
                 y: {{length: {H2}, n: 3}}}}
            groups:
              O1: "near(x, 0) and near(y, 0)"
              O2: "near(x, 0) and near(y, {H1 + H2})"
        interfaces:
          - {{name: junta, split: abajo.ymax, material: 1, integration: {integration}}}
        boundary_conditions_by_group:
          abajo.ymin: {{uy: 0.0}}
          O1: {{ux: 0.0}}
          O2: {{ux: 0.0}}
        tractions_by_group:
          arriba.ymax: {{pressure: {P}}}
        solver: {{type: NonlinearSolver, num_steps: 1}}
        """), encoding="utf-8")
    model = solidum.load_yaml(ruta)
    res = model.solve()
    dom = model.domain
    (made,) = model.parser.applied_interfaces
    assert len(made.elements) == 4
    for eid in made.elements:
        e = dom.elements[eid]
        assert type(e) is Interface2D6
        for t in e.state.stresses:
            np.testing.assert_allclose(t, [-P, 0.0], rtol=1e-9, atol=1e-9 * P)
    for el in dom.elements.values():
        if type(el) is not Interface2D6:
            for sig in el.state.stresses:
                np.testing.assert_allclose(sig[1], -P, rtol=1e-9)
    tope = dom.physical_groups["arriba.ymax"]
    U = res.U
    for nid in tope:
        assert U[dom.nodes[nid].dofs["uy"]] == pytest.approx(_u_tope(), rel=1e-9)


def test_parcela_con_pair_y_nodos_coincidentes():
    """Dos mallas de Quad8 sin fundir, unidas con ``pair``: la misma parcela
    por desplazamiento impuesto arriba."""
    mat1 = Elastic2D(E=E1, nu=NU1, hypothesis="plane_strain")
    mat2 = Elastic2D(E=E2, nu=NU2, hypothesis="plane_strain")
    dom = structured_mesh("Quad8", lines=[np.linspace(0, 2, 5), np.linspace(0, H1, 3)],
                          material=mat1, thickness=E_T).domain
    n_a = len(dom.nodes)
    structured_mesh("Quad8", lines=[np.linspace(0, 2, 5), np.linspace(H1, H1 + H2, 4)],
                    material=mat2, thickness=E_T, domain=dom)
    y = {i: float(n.coordinates[1]) for i, n in dom.nodes.items()}
    x = {i: float(n.coordinates[0]) for i, n in dom.nodes.items()}
    dom.physical_groups = {"a": [i for i in y if abs(y[i] - H1) < 1e-12 and i <= n_a],
                           "b": [i for i in y if abs(y[i] - H1) < 1e-12 and i > n_a]}
    dom.physical_group_cells = {"a": [], "b": []}
    r = pair_interfaces(dom, "a", "b", ContactPenalty(K_n=KN_PARCELA), name="junta")
    assert len(r.elements) == 4 and all(type(dom.elements[e]) is Interface2D6 for e in r.elements)
    for i, n in dom.nodes.items():
        if abs(y[i]) < 1e-12:
            n.fix_dof("uy", 0.0)
        if abs(y[i] - (H1 + H2)) < 1e-12:
            n.fix_dof("uy", _u_tope())
        if abs(x[i]) < 1e-12 and (abs(y[i]) < 1e-12 or abs(y[i] - H1 - H2) < 1e-12):
            n.fix_dof("ux", 0.0)
    dom.generate_equation_numbers()
    asm = Assembler(dom)
    solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=1))
    for eid in r.elements:
        for t in dom.elements[eid].state.stresses:
            np.testing.assert_allclose(t, [-P, 0.0], rtol=1e-8, atol=1e-8 * P)
