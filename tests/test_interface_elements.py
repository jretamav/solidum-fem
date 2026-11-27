"""Aceptación de los elementos de interfaz de espesor nulo (specs
``Interface2D4``, ``Interface3D6``, ``Interface3D8``; ADR 0043).

Los casos siguen el bloque ``acceptance`` de las specs. Parámetros ≠ 1 (espesor,
rigideces, módulos) para que un factor olvidado no pase como identidad
(Reglas §6).
"""
from __future__ import annotations

import numpy as np
import pytest

import solidum
from solidum.core.node import Node
from solidum.elements.interface import Interface2D4, Interface3D6, Interface3D8, InterfaceElement
from solidum.elements.solid_2d import Quad4
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.interface import CohesiveDamageIsotropic, ContactPenalty
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver
from solidum.registry import ElementRegistry

K_N = 4.0e5
K_T = 2.5e5
MU = 0.35


def _nodes(coords, first=1):
    return [Node(first + i, list(map(float, c))) for i, c in enumerate(coords)]


def iface2d(xa, xb, law=None, **kw):
    """Interfaz 2D con cara A = (xa[0], xa[1]) y cara B frente a ella
    (xb[0] frente a xa[0])."""
    law = law or ContactPenalty(K_n=K_N, mu=MU, K_t=K_T)
    nodes = _nodes([xa[0], xa[1], xb[1], xb[0]])
    return Interface2D4(1, nodes, law, **kw)


def iface3d8(xa, xb=None, law=None, **kw):
    law = law or ContactPenalty(K_n=K_N, mu=MU, K_t=K_T)
    xb = xa if xb is None else xb
    return Interface3D8(1, _nodes(list(xa) + list(xb)), law, **kw)


def iface3d6(xa, xb=None, law=None, **kw):
    law = law or ContactPenalty(K_n=K_N, mu=MU, K_t=K_T)
    xb = xa if xb is None else xb
    return Interface3D6(1, _nodes(list(xa) + list(xb)), law, **kw)


def fd_stiffness(e, u, h=1e-9):
    snap = e.state.snapshot_trial()
    n = u.size
    K = np.zeros((n, n))
    for j in range(n):
        d = np.zeros(n)
        d[j] = h
        _, Fp = e.compute_element_state(u + d)
        _, Fm = e.compute_element_state(u - d)
        K[:, j] = (Fp - Fm) / (2 * h)
    e.state.restore_trial(snap)
    return K


XA_2D = np.array([[0.3, 0.2], [2.3, 1.1]])           # cara inclinada
QUAD_PLANA = np.array([[0, 0, 0], [2, 0, 0], [1.6, 1.1, 0], [0.2, 0.7, 0]], float)
QUAD_ALABEADA = np.array([[0, 0, 0], [2, 0, 0.3], [1.6, 1.1, -0.2], [0.2, 0.7, 0.25]], float)
TRI = np.array([[0.1, 0.0, 0.2], [1.9, 0.3, 0.0], [0.7, 1.4, 0.4]])


# -- cinemática ---------------------------------------------------------------------

@pytest.mark.parametrize("integration", ["nodal", "gauss"])
def test_cuerpo_rigido(integration):
    e = iface2d(XA_2D, XA_2D, integration=integration, thickness=0.7)
    X = np.vstack([XA_2D, XA_2D[::-1]])
    omega, c = 1e-6, np.array([0.4, -0.1])
    u_rot = np.concatenate([omega * np.array([-(x[1] - c[1]), x[0] - c[0]]) for x in X])
    u_tr = np.tile([0.013, -0.021], 4)
    for u in (u_tr, u_rot, u_tr + u_rot):
        np.testing.assert_array_equal(e.jumps(u), 0.0)
        _, F = e.compute_element_state(u)
        np.testing.assert_allclose(F, 0.0, atol=1e-12)


def test_salto_uniforme_y_fuerzas_nodales():
    """Cara B trasladada en la normal: tracción uniforme; fuerzas nodales
    = t·L/2·e en cada nodo con las dos reglas."""
    law = ContactPenalty(K_n=K_N)
    L = float(np.linalg.norm(XA_2D[1] - XA_2D[0]))
    s = (XA_2D[1] - XA_2D[0]) / L
    n = np.array([-s[1], s[0]])
    e_t = 0.7
    for integration in ("nodal", "gauss"):
        e = iface2d(XA_2D, XA_2D, law=law, integration=integration, thickness=e_t)
        du = -2e-4 * n                                   # B entra en A
        u = np.concatenate([[0, 0], [0, 0], du, du])
        _, F = e.compute_element_state(u)
        t_n = K_N * (-2e-4)
        np.testing.assert_allclose(e.state.stresses_trial[0], [t_n, 0.0], rtol=1e-12)
        f_node = t_n * L / 2 * e_t
        np.testing.assert_allclose(F[4:6], f_node * n, rtol=1e-12)
        np.testing.assert_allclose(F[0:2], -f_node * n, rtol=1e-12)


def test_holgura_inicial_y_cierre():
    g0 = 3e-3
    n = iface2d(XA_2D, XA_2D).R[0, 0]
    e = iface2d(XA_2D, XA_2D + g0 * n, law=ContactPenalty(K_n=K_N))   # B separada en n
    np.testing.assert_allclose(e.g0, g0, rtol=1e-12)
    for closing in (0.5, 0.99, 1.4):
        du = -closing * g0 * n
        u = np.concatenate([[0, 0], [0, 0], du, du])
        _, F = e.compute_element_state(u)
        if closing <= 1.0:
            np.testing.assert_allclose(F, 0.0, atol=1e-12)
        else:
            assert np.linalg.norm(F) > 0.0


def test_holgura_variable_en_el_elemento():
    # B inclinada respecto de A: pares sobre la normal de A con holgura
    # distinta en cada extremo.
    L = np.linalg.norm(XA_2D[1] - XA_2D[0])
    s = (XA_2D[1] - XA_2D[0]) / L
    n = np.array([-s[1], s[0]])
    gaps = np.array([1e-3, 4e-3])
    xb = XA_2D + gaps[:, None] * n
    for integration, xi in (("nodal", np.array([-1.0, 1.0])),
                            ("gauss", np.array([-1, 1]) / np.sqrt(3.0))):
        e = iface2d(XA_2D, xb, law=ContactPenalty(K_n=K_N), integration=integration)
        expected = 0.5 * (1 - xi) * gaps[0] + 0.5 * (1 + xi) * gaps[1]
        np.testing.assert_allclose(e.g0, expected, rtol=1e-12)


def test_giro_rigido_con_holgura():
    """Límite documentado (§9): con holgura, un giro rígido ω da g_s = −ω·g_0."""
    g0 = 0.07
    xa = np.array([[0.0, 0.0], [1.0, 0.0]])
    e = iface2d(xa, xa + [0.0, g0], law=ContactPenalty(K_n=K_N))
    X = np.vstack([xa, (xa + [0.0, g0])[::-1]])
    omega = 1e-4
    u = np.concatenate([omega * np.array([-x[1], x[0]]) for x in X])
    g = e.jumps(u)
    np.testing.assert_allclose(g[:, 1], -omega * g0, rtol=1e-10)
    np.testing.assert_allclose(g[:, 0], g0, rtol=1e-10)


# -- tangente y simetría --------------------------------------------------------------

@pytest.mark.parametrize("integration", ["nodal", "gauss"])
@pytest.mark.parametrize("case", ["abierta", "adherencia", "deslizamiento"])
def test_tangente_del_elemento_2d(integration, case):
    e = iface2d(XA_2D, XA_2D, integration=integration, thickness=0.7)
    s = e.R[0, 1]
    n = e.R[0, 0]
    rng = np.random.default_rng(3)
    base = {"abierta": 2e-3 * n, "adherencia": -1e-3 * n + 1e-5 * s,
            "deslizamiento": -1e-3 * n + 3e-3 * s}[case]
    u = np.concatenate([[0, 0], [0, 0], base, base]) + 1e-7 * rng.standard_normal(8)
    K, _ = e.compute_element_state(u)
    np.testing.assert_allclose(K, fd_stiffness(e, u), rtol=1e-6, atol=1e-6 * K_N)


@pytest.mark.parametrize("make", [lambda **k: iface3d8(QUAD_ALABEADA, **k),
                                  lambda **k: iface3d6(TRI, **k)])
@pytest.mark.parametrize("case", ["adherencia", "deslizamiento"])
def test_tangente_del_elemento_3d(make, case):
    e = make()
    nf = e.N_FACE_NODES
    rng = np.random.default_rng(5)
    n, s1, s2 = e.R[0]
    shift = (-1e-3 * n + (1e-5 if case == "adherencia" else 3e-3) * (0.6 * s1 + 0.8 * s2))
    u = np.concatenate([np.zeros(3 * nf), np.tile(shift, nf)]) + 1e-8 * rng.standard_normal(6 * nf)
    K, _ = e.compute_element_state(u)
    np.testing.assert_allclose(K, fd_stiffness(e, u), rtol=1e-6, atol=1e-6 * K_N)


def test_simetria_de_la_rigidez():
    e = iface2d(XA_2D, XA_2D)
    s, n = e.R[0, 1], e.R[0, 0]
    base = -1e-3 * n + 3e-3 * s
    u = np.concatenate([[0, 0], [0, 0], base, base])
    K, _ = e.compute_element_state(u)
    assert not np.allclose(K, K.T)
    e0 = iface2d(XA_2D, XA_2D, law=ContactPenalty(K_n=K_N))
    K0, _ = e0.compute_element_state(u)
    np.testing.assert_allclose(K0, K0.T, rtol=1e-14, atol=1e-14 * K_N)


def test_nodal_desacopla_pares():
    e = iface2d(XA_2D, XA_2D, law=ContactPenalty(K_n=K_N))
    u = np.concatenate([[0, 0], [0, 0], [0, -1e-3], [0, -1e-3]])
    K, _ = e.compute_element_state(u)
    # Par 1: nodos 1 (dof 0-1) y 4 (dof 6-7); par 2: nodos 2 (2-3) y 3 (4-5).
    pair1, pair2 = [0, 1, 6, 7], [2, 3, 4, 5]
    assert np.all(K[np.ix_(pair1, pair2)] == 0.0)
    eg = iface2d(XA_2D, XA_2D, law=ContactPenalty(K_n=K_N), integration="gauss")
    Kg, _ = eg.compute_element_state(u)
    assert np.any(Kg[np.ix_(pair1, pair2)] != 0.0)


def test_orientacion_invertida():
    """Nodos coincidentes: intercambiar A y B da las mismas fuerzas nodales."""
    law = ContactPenalty(K_n=K_N, mu=MU, K_t=K_T)
    e1 = iface2d(XA_2D, XA_2D, law=law)
    n, s = e1.R[0]
    dB = -1e-3 * n + 2e-3 * s
    u1 = np.concatenate([[0, 0], [0, 0], dB, dB])     # A = (1, 2), B = (3, 4) se mueve
    _, F1 = e1.compute_element_state(u1)
    # Elemento invertido: A' = cara B (orden contrario), B' = cara A.
    nodes2 = _nodes([XA_2D[1], XA_2D[0], XA_2D[0], XA_2D[1]])
    e2 = Interface2D4(2, nodes2, law)
    # Nodos de e2: A' = (B frente a 2, B frente a 1), B' = (1, 2).
    u2 = np.concatenate([dB, dB, [0, 0], [0, 0]])
    _, F2 = e2.compute_element_state(u2)
    # Fuerzas sobre cada nodo físico: e1 = [A1, A2, B2, B1]; e2 = [B2, B1, A1, A2].
    np.testing.assert_allclose(F2[[4, 5, 6, 7, 0, 1, 2, 3]], F1, rtol=1e-12, atol=1e-12)


# -- escalado, áreas, ejes --------------------------------------------------------------

def test_escalado_con_el_espesor():
    u = np.concatenate([[0, 0], [0, 0], [1e-4, -1e-3], [1e-4, -1e-3]])
    e1 = iface2d(XA_2D, XA_2D, thickness=1.0)
    e2 = iface2d(XA_2D, XA_2D, thickness=0.37)
    K1, F1 = e1.compute_element_state(u)
    K2, F2 = e2.compute_element_state(u)
    np.testing.assert_allclose(K2, 0.37 * K1, rtol=1e-12)
    np.testing.assert_allclose(F2, 0.37 * F1, rtol=1e-12)
    e1.state.commit()
    e2.state.commit()
    U = u  # numeración trivial: get_local_displacements necesita el dominio
    psi1 = sum(e1.measures[g] * e1.material.stored_energy(j, e1.state.vars[g])
               for g, j in enumerate(e1.jumps(U)))
    psi2 = sum(e2.measures[g] * e2.material.stored_energy(j, e2.state.vars[g])
               for g, j in enumerate(e2.jumps(U)))
    assert psi2 == pytest.approx(0.37 * psi1, rel=1e-12)


def test_area_tributaria_en_cara_trapezoidal():
    e = iface3d8(QUAD_PLANA)
    np.testing.assert_allclose(e.measures, [0.375, 0.44166666666666665, 0.4, 0.3333333333333333],
                               rtol=1e-12)
    eg = iface3d8(QUAD_PLANA, integration="gauss")
    assert eg.measures.sum() == pytest.approx(1.55, rel=1e-14)
    tri = iface3d6(TRI)
    area = 0.5 * np.linalg.norm(np.cross(TRI[1] - TRI[0], TRI[2] - TRI[0]))
    np.testing.assert_allclose(tri.measures, area / 3, rtol=1e-12)


def test_cara_alabeada_ejes():
    e = iface3d8(QUAD_ALABEADA, integration="gauss")
    for R in e.R:
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-14)
        assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-14)
    # La normal cambia dentro de la cara alabeada.
    assert not np.allclose(e.R[0, 0], e.R[2, 0])


# -- validación al construir --------------------------------------------------------------

def test_cohesiva_con_holgura_rechazada():
    law = CohesiveDamageIsotropic(sigma_t0=3e6, G_f=100.0, K_e=1e12, softening="linear")
    with pytest.raises(ValueError, match="no admite holgura inicial"):
        iface2d(XA_2D, XA_2D + 1e-3 * iface2d(XA_2D, XA_2D).R[0, 0], law=law)
    iface2d(XA_2D, XA_2D, law=law)          # sin holgura, se admite


def test_interpenetracion_inicial():
    L = np.linalg.norm(XA_2D[1] - XA_2D[0])
    n = iface2d(XA_2D, XA_2D).R[0, 0]
    # Bajo la tolerancia (1e-6·h): se lleva a 0.
    e = iface2d(XA_2D, XA_2D - 0.5e-6 * L * n, law=ContactPenalty(K_n=K_N))
    np.testing.assert_array_equal(e.g0, 0.0)
    with pytest.raises(ValueError, match="se interpenetran al inicio"):
        iface2d(XA_2D, XA_2D - 5e-6 * L * n, law=ContactPenalty(K_n=K_N))


def test_par_fuera_de_la_normal_rechazado():
    s = iface2d(XA_2D, XA_2D).R[0, 1]
    with pytest.raises(ValueError, match="desfase tangencial"):
        iface2d(XA_2D, XA_2D + 1e-3 * s, law=ContactPenalty(K_n=K_N))


def test_leyes_y_materiales_rechazados():
    from solidum.materials.interface import BondSlipMultilinear
    with pytest.raises(ValueError, match="superficie"):
        iface2d(XA_2D, XA_2D, law=Elastic2D(E=1.0, nu=0.2))
    line_law = BondSlipMultilinear(points=[[1e-3, 1e6]], k_n=1e12)
    with pytest.raises(ValueError, match="superficie"):
        iface2d(XA_2D, XA_2D, law=line_law)
    with pytest.raises(ValueError, match="integration"):
        iface2d(XA_2D, XA_2D, integration="lobatto")


def test_registro():
    """Todas las clases de interfaz de espesor nulo del registro (derivadas,
    no escritas a mano): sin ``STRAIN_DIM``, con ``JUMP_DIM`` 2 o 3 y su tipo
    de cara en la tabla de la base."""
    from solidum.elements.interface._shared import FACE_CELLS
    clases = [ElementRegistry.get(n) for n in ElementRegistry.names()
              if issubclass(ElementRegistry.get(n), InterfaceElement)]
    assert {c.__name__ for c in clases} >= {"Interface2D4", "Interface2D6", "Interface3D6",
                                            "Interface3D8"}
    for cls in clases:
        assert cls.STRAIN_DIM is None and cls.JUMP_DIM in (2, 3)
        assert cls.FACE_CELL in FACE_CELLS


@pytest.mark.parametrize("atributos, mensaje", [
    (dict(FACE_CELL="quadd"), "no es un tipo de cara"),
    (dict(N_FACE_NODES=3), "N_FACE_NODES"),
    (dict(B_NODES=(4, 5, 6, 6)), "deben repartir"),
    (dict(N_INTEGRATION_POINTS=9), "N_INTEGRATION_POINTS"),
    (dict(JUMP_DIM=2), "JUMP_DIM"),
])
def test_una_subclase_mal_declarada_falla_al_definirse(atributos, mensaje):
    """El contrato de la subclase se comprueba al definirla (revisión del
    arquitecto del 2026-11-26): un ``FACE_CELL`` mal escrito ya no cae en la
    regla de otra cara."""
    base = dict(DOF_NAMES=["ux", "uy", "uz"], JUMP_DIM=3, FACE_CELL="quad", N_FACE_NODES=4,
                A_NODES=(0, 1, 2, 3), B_NODES=(4, 5, 6, 7), N_INTEGRATION_POINTS=4,
                VTK_CELL_TYPE="hexahedron")
    base.update(atributos)
    with pytest.raises(TypeError, match=mensaje):
        type("Mala", (InterfaceElement,), base)
    del base["FACE_CELL"]
    type("Abstracta", (InterfaceElement,), base)    # sin cara: base intermedia, sin comprobar


# -- parcela de contacto con sólidos ------------------------------------------------------

def _bloque(dom, nid, eid, y0, H, Lx, nx, ny, mat, e_t):
    """Bloque de Quad4 en [0, Lx]×[y0, y0+H]; devuelve la rejilla de nodos."""
    grid = np.empty((ny + 1, nx + 1), dtype=object)
    for j in range(ny + 1):
        for i in range(nx + 1):
            grid[j, i] = dom.add_node(nid, [Lx * i / nx, y0 + H * j / ny])
            nid += 1
    for j in range(ny):
        for i in range(nx):
            q = [grid[j, i], grid[j, i + 1], grid[j + 1, i + 1], grid[j + 1, i]]
            dom.add_element(Quad4(eid, q, mat, thickness=e_t))
            eid += 1
    return grid, nid, eid


@pytest.mark.parametrize("integration", ["nodal", "gauss"])
@pytest.mark.parametrize("g0", [0.0, 2e-3])
def test_parcela_de_contacto_2d(integration, g0):
    """Dos bloques de materiales distintos; base con u_y = 0 y un u_x por bloque.
    Con g_0 = 0, presión p arriba; con g_0 > 0, bajada impuesta (en carga
    controlada el bloque superior flotaría). Exacto: t_n uniforme,
    u_top = −p·Σ(1−ν²)H/E − p/K_n − g_0."""
    E1, nu1, H1 = 2.1e4, 0.2, 0.7
    E2, nu2, H2 = 0.9e4, 0.3, 1.3
    Lx, nx, e_t = 2.0, 3, 0.6
    m1 = Elastic2D(E=E1, nu=nu1, hypothesis="plane_strain")
    m2 = Elastic2D(E=E2, nu=nu2, hypothesis="plane_strain")
    dom = solidum.Domain()
    g1, nid, eid = _bloque(dom, 1, 1, 0.0, H1, Lx, nx, 2, m1, e_t)
    g2, nid, eid = _bloque(dom, nid, eid, H1 + g0, H2, Lx, nx, 2, m2, e_t)
    law = ContactPenalty(K_n=K_N)
    for i in range(nx):
        # Cara A: arista superior del bloque 1 recorrida en sentido contrario
        # al de su contorno antihorario (§2), para que n apunte hacia B.
        a1, a2 = g1[-1, i], g1[-1, i + 1]
        b1, b2 = g2[0, i], g2[0, i + 1]            # b1 frente a a1, b2 frente a a2
        dom.add_element(Interface2D4(eid, [a1, a2, b2, b1], law,
                                     integration=integration, thickness=e_t))
        eid += 1
    for node in g1[0, :]:
        node.fix_dof("uy", 0.0)
    g1[0, 0].fix_dof("ux", 0.0)
    g2[0, 0].fix_dof("ux", 0.0) if g0 == 0.0 else g2[-1, 0].fix_dof("ux", 0.0)
    flex = (1 - nu1 ** 2) * H1 / E1 + (1 - nu2 ** 2) * H2 / E2 + 1.0 / K_N
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    if g0 == 0.0:
        p = 35.0
        for i in range(nx):                      # presión consistente en la cara superior
            for node in (g2[-1, i], g2[-1, i + 1]):
                F[node.dofs["uy"]] -= p * (Lx / nx) / 2 * e_t
    else:
        delta = 1.5e-3
        for node in g2[-1, :]:
            node.fix_dof("uy", -(g0 + delta))
        p = delta / flex
    asm = Assembler(dom)
    res = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=2), F_applied=F)
    U = res.U
    for el in dom.elements.values():
        if isinstance(el, Interface2D4):
            for t in el.state.stresses:
                np.testing.assert_allclose(t, [-p, 0.0], rtol=1e-9, atol=1e-9 * p)
    u_top = U[g2[-1, nx // 2].dofs["uy"]]
    assert u_top == pytest.approx(-p * flex - g0, rel=1e-9)


def test_contrato_de_posiciones_de_las_caras():
    """Cada elemento de interfaz registrado reparte sus nodos entre la cara A
    (``A_NODES``, en el orden de sus funciones de forma) y sus pares en B
    (``B_NODES``): cada posición una sola vez."""
    clases = [c for c in (ElementRegistry.get(n) for n in ElementRegistry.names())
              if isinstance(c, type) and issubclass(c, InterfaceElement)]
    assert clases
    for cls in clases:
        nf = cls.N_FACE_NODES
        assert len(cls.A_NODES) == len(cls.B_NODES) == nf, cls.__name__
        assert sorted(cls.A_NODES + cls.B_NODES) == list(range(2 * nf)), cls.__name__
