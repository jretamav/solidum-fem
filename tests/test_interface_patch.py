"""Parcelas de contacto con sólidos e interfaces (specs ``Interface2D4`` §7b,
``Interface3D6``, ``Interface3D8`` e ``Interface2D6``; verificación de
``teorico-solidum`` del 2026-11-26).

- Parcela 3D (``parcela_de_contacto_3d`` de las specs 3D): dos bloques de
  materiales distintos con malla distorsionada (caras de la interfaz planas y
  no paralelogramo en Hex8; triángulos en Tet4), bajada impuesta arriba:
  presión de contacto, esfuerzo σ_zz e interpenetración exactos, con holgura
  inicial y sin ella, nodal y Gauss.
- Cara cuadrilátera alabeada bajo presión hidrostática: con la regla nodal la
  normal y la medida de cada nodo son las del vector de área tributaria
  consistente ``∫N_i n J`` (en ``(ξ_i/3, η_i/3)``; decisión del usuario del
  2026-11-26). Con la normal de la esquina el error era de 8.7 % en los
  esfuerzos con 7.2° de alabeo.
- Lado curvo de ``Interface2D6`` bajo presión hidrostática (Simpson con el
  jacobiano del nodo): el error baja con ``1/K_n``.
- Con fricción, el modelo pide el solver no simétrico (ADR 0003).
"""
import logging

import numpy as np
import pytest

import solidum
from solidum.elements.interface import InterfaceElement
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.materials.interface import ContactPenalty
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver
from solidum.math.solvers._shared import domain_is_symmetric
from solidum.mesh import structured_mesh
from solidum.mesh.interfaces import pair_interfaces


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.WARNING)
    yield
    logging.disable(logging.NOTSET)


def _grupo(dom, nombre, ids):
    if not getattr(dom, "physical_groups", None):
        dom.physical_groups, dom.physical_group_cells = {}, {}
    dom.physical_groups[nombre] = sorted(ids)
    dom.physical_group_cells[nombre] = []


def _interfaces(dom):
    return [e for e in dom.elements.values() if isinstance(e, InterfaceElement)]


def _solidos(dom):
    return [e for e in dom.elements.values() if not isinstance(e, InterfaceElement)]


def _isostatico_3d(bloque):
    """Un nodo con ux, uy y otro con uy: quita los movimientos rígidos en el
    plano sin impedir el efecto de Poisson."""
    P = bloque[0]
    Q = max(bloque, key=lambda q: q.coordinates[0] - P.coordinates[0])
    P.fix_dof("ux", 0.0)
    P.fix_dof("uy", 0.0)
    Q.fix_dof("uy", 0.0)


# --- parcela 3D -----------------------------------------------------------------

def _parcela_3d(elemento, integracion, g0, K_n=4.0e5):
    E1, nu1, H1 = 2.1e4, 0.2, 0.7
    E2, nu2, H2 = 0.9e4, 0.3, 1.3
    Lx, Ly = 2.0, 1.5

    def mapeo(z0, H):
        def f(X):
            x, y, z = X
            zeta = (z - z0) / H
            return np.array([x + 0.11 * np.sin(1.3 * x + 0.4) * np.cos(0.9 * y + 0.2),
                             y + 0.09 * np.cos(1.1 * x) * np.sin(1.7 * y + 0.3),
                             z + 0.06 * np.sin(np.pi * zeta) * np.sin(x + 2 * y)])
        return f

    def lineas(z0, H):
        return [np.linspace(0, Lx, 4), np.linspace(0, Ly, 3), z0 + np.linspace(0, H, 3)]

    m1 = structured_mesh(elemento, lines=lineas(0.0, H1), mapping=mapeo(0.0, H1),
                         material=Elastic3D(E=E1, nu=nu1))
    dom = m1.domain
    n1 = max(dom.nodes)
    z2 = H1 + g0
    structured_mesh(elemento, lines=lineas(z2, H2), mapping=mapeo(z2, H2),
                    material=Elastic3D(E=E2, nu=nu2), domain=dom)
    _grupo(dom, "a", [i for i, n in dom.nodes.items() if abs(n.coordinates[2] - H1) < 1e-9 and i <= n1])
    _grupo(dom, "b", [i for i, n in dom.nodes.items() if abs(n.coordinates[2] - z2) < 1e-9 and i > n1])
    pair_interfaces(dom, "a", "b", ContactPenalty(K_n=K_n), name="c", integration=integracion)
    delta = 1.5e-3
    base = [n for n in dom.nodes.values() if abs(n.coordinates[2]) < 1e-9]
    tapa = [n for n in dom.nodes.values() if abs(n.coordinates[2] - (z2 + H2)) < 1e-9]
    for n in base:
        n.fix_dof("uz", 0.0)
    for n in tapa:
        n.fix_dof("uz", -(g0 + delta))
    _isostatico_3d(base)
    _isostatico_3d(tapa)
    dom.generate_equation_numbers()
    asm = Assembler(dom)
    res = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=2),
                      F_applied=np.zeros(dom.total_dofs))
    p = delta / (H1 / E1 + H2 / E2 + 1.0 / K_n)       # resortes en serie
    return dom, res, p


@pytest.mark.parametrize("g0", [0.0, 2.0e-3])
@pytest.mark.parametrize("integracion", ["nodal", "gauss"])
@pytest.mark.parametrize("elemento, clase", [("Hex8", "Interface3D8"), ("Tet4", "Interface3D6")])
def test_parcela_de_contacto_3d_presion_esfuerzo_e_interpenetracion_exactos(elemento, clase,
                                                                            integracion, g0):
    dom, res, p = _parcela_3d(elemento, integracion, g0)
    ifs = _interfaces(dom)
    assert {type(e).__name__ for e in ifs} == {clase}
    K_n = ifs[0].material.K_n
    tn = np.array([t[0] for e in ifs for t in e.state.stresses])
    ts = np.array([np.linalg.norm(t[1:]) for e in ifs for t in e.state.stresses])
    gn = np.concatenate([e.jumps(e.get_local_displacements(res.U))[:, 0] for e in ifs])
    szz = np.array([s[2] for e in _solidos(dom) for s in e.state.stresses])
    assert np.max(np.abs(tn + p)) < 1e-10 * p
    assert np.max(ts) < 1e-10 * p
    assert np.max(np.abs(gn + p / K_n)) < 1e-8 * (p / K_n)
    assert np.max(np.abs(szz + p)) < 1e-10 * p


def test_parcela_3d_hex8_tiene_caras_no_paralelogramo():
    """La parcela 3D pone a prueba los pesos nodales de caras no paralelogramo
    (con los de Newton–Cotes el error era de 4.3e-2): las esquinas de cada
    cara tienen jacobianos distintos."""
    dom, _, _ = _parcela_3d("Hex8", "nodal", 0.0)
    cociente = []
    for e in _interfaces(dom):
        XA, _ = e._face_coordinates()
        J = [e.local_axes(XA, xi)[1] for xi in ((-1., -1.), (1., -1.), (1., 1.), (-1., 1.))]
        cociente.append(max(J) / min(J))
    assert max(cociente) > 1.1


# --- cara cuadrilátera alabeada ---------------------------------------------------

def _parcela_alabeada(integracion, c=0.15, K_n=4.0e10, n=3, curvatura=0.0):
    """Dos bloques Hex8 del mismo material separados por un paraboloide
    hiperbólico ``z = H1 + c·(x − 1)(y − 0.75)`` (caras alabeadas; con
    ``curvatura``, más ``curvatura·(x − 1)²``: torsión y curvatura a la vez),
    con presión hidrostática ``p`` en todo el contorno salvo la base: σ = −p·I
    y ``t = −p·n`` en cualquier punto de la interfaz."""
    E, nu, H1, H2, Lx, Ly, p = 2.1e4, 0.25, 0.8, 1.0, 2.0, 1.5, 35.0

    def sup(x, y):
        return c * (x - Lx / 2) * (y - Ly / 2) + curvatura * (x - Lx / 2) ** 2

    def mapeo(z0, H, abajo):
        def f(X):
            x, y, z = X
            zeta = (z - z0) / H
            peso = (1 - zeta) if abajo else zeta          # 1 en la cara de la interfaz
            return np.array([x, y, z + peso * sup(x, y)])
        return f

    mat = Elastic3D(E=E, nu=nu)

    def lineas(z0, H):
        return [np.linspace(0, Lx, n + 1), np.linspace(0, Ly, n + 1), z0 + np.linspace(0, H, 3)]

    m1 = structured_mesh("Hex8", lines=lineas(0, H1), mapping=mapeo(0, H1, False), material=mat)
    dom = m1.domain
    n1 = max(dom.nodes)
    m2 = structured_mesh("Hex8", lines=lineas(H1, H2), mapping=mapeo(H1, H2, True), material=mat,
                         domain=dom)
    en = [i for i, nd in dom.nodes.items()
          if abs(nd.coordinates[2] - (H1 + sup(nd.coordinates[0], nd.coordinates[1]))) < 1e-9]
    _grupo(dom, "a", [i for i in en if i <= n1])
    _grupo(dom, "b", [i for i in en if i > n1])
    pair_interfaces(dom, "a", "b", ContactPenalty(K_n=K_n), name="c", integration=integracion)
    base = [nd for nd in dom.nodes.values() if abs(nd.coordinates[2]) < 1e-9]
    tapa = [nd for nd in dom.nodes.values() if abs(nd.coordinates[2] - (H1 + H2)) < 1e-9]
    for nd in base:
        nd.fix_dof("uz", 0.0)
    _isostatico_3d(base)
    _isostatico_3d(tapa)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)

    def presion(x, nn):
        return -p * nn

    for m, lados in ((m1, ("xmin", "xmax", "ymin", "ymax")),
                     (m2, ("xmin", "xmax", "ymin", "ymax", "zmax"))):
        for lado in lados:
            for el, k in m.boundary[lado]:
                F[el.get_global_dof_indices()] += el.compute_face_traction(k, presion)
    asm = Assembler(dom)
    solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=1), F_applied=F)
    ifs = _interfaces(dom)
    alabeo = 0.0
    for e in ifs:
        XA, _ = e._face_coordinates()
        ns = [e.local_axes(XA, q)[0][0] for q in ((-1., -1.), (1., -1.), (1., 1.), (-1., 1.))]
        alabeo = max(alabeo, max(np.degrees(np.arccos(np.clip(a @ b, -1, 1)))
                                 for a in ns for b in ns))
    tn = np.array([t[0] for e in ifs for t in e.state.stresses])
    sig = np.array([s[:3] for e in _solidos(dom) for s in e.state.stresses])
    return np.max(np.abs(tn + p)) / p, np.max(np.abs(sig + p)) / p, alabeo


@pytest.mark.parametrize("integracion", ["nodal", "gauss"])
def test_cara_alabeada_transmite_presion_uniforme(integracion):
    """Con 7.2° de alabeo por cara: antes, con la normal de la esquina, la
    regla nodal daba 1.8e-2 en la tracción y 8.7e-2 en los esfuerzos y no
    bajaba con K_n (verificación del 2026-11-26, *medido*)."""
    err_t, err_s, alabeo = _parcela_alabeada(integracion)
    assert alabeo > 7.0
    assert err_t < 1e-4
    assert err_s < 1e-6


def test_torsion_y_curvatura_residuo_nodal_de_segundo_orden():
    """Con torsión y curvatura, en cada nodo se juntan caras de normales
    distintas y la regla nodal deja un residuo en la tracción que baja ~4 veces
    al duplicar la malla (1.49e-3 → 3.75e-4, *medido* por el teórico); con la
    normal de la esquina era 17–18 % y no bajaba al refinar."""
    e3, _, _ = _parcela_alabeada("nodal", c=0.15, curvatura=0.09, n=3)
    e6, _, _ = _parcela_alabeada("nodal", c=0.15, curvatura=0.09, n=6)
    assert e3 < 3e-3
    assert e6 / e3 < 0.3


def test_cara_plana_regla_nodal_igual_al_area_tributaria():
    """En una cara plana, la medida nodal en (ξ_i/3, η_i/3) es el área
    tributaria consistente ∫N_i J de la decisión Q3 (J lineal)."""
    dom, _, _ = _parcela_3d("Hex8", "nodal", 0.0)
    x, w = np.polynomial.legendre.leggauss(4)
    esquinas = ((-1., -1.), (1., -1.), (1., 1.), (-1., 1.))
    for e in _interfaces(dom):
        XA, _ = e._face_coordinates()
        ref = np.zeros(4)
        for a, wa in zip(x, w):
            for b, wb in zip(x, w):
                N = np.array([(1 + xi * a) * (1 + eta * b) / 4 for xi, eta in esquinas])
                ref += wa * wb * e.local_axes(XA, (a, b))[1] * N
        np.testing.assert_allclose(e.measures, ref, rtol=1e-12)


# --- lado curvo de Interface2D6 -----------------------------------------------------

def _parcela_2d6_curva(integracion, K_n, flecha=0.10):
    E, nu, H1, H2, Lx, e_t, p = 2.1e4, 0.25, 0.8, 1.0, 2.0, 0.6, 35.0

    def curva(x):
        return flecha * (1 - (x - Lx / 2) ** 2 / (Lx / 2) ** 2)

    def mapeo(y0, H, abajo):
        def f(X):
            x, y = X
            zeta = (y - y0) / H
            peso = (1 - zeta) if abajo else zeta
            return np.array([x, y + peso * curva(x)])
        return f

    mat = Elastic2D(E=E, nu=nu, hypothesis="plane_strain")
    m1 = structured_mesh("Quad8", lines=[np.linspace(0, Lx, 3), np.linspace(0, H1, 3)],
                         mapping=mapeo(0.0, H1, False), material=mat, thickness=e_t)
    dom = m1.domain
    n1 = max(dom.nodes)
    m2 = structured_mesh("Quad8", lines=[np.linspace(0, Lx, 3), H1 + np.linspace(0, H2, 3)],
                         mapping=mapeo(H1, H2, True), material=mat, thickness=e_t, domain=dom)
    en = [i for i, nd in dom.nodes.items()
          if abs(nd.coordinates[1] - (H1 + curva(nd.coordinates[0]))) < 1e-9]
    _grupo(dom, "a", [i for i in en if i <= n1])
    _grupo(dom, "b", [i for i in en if i > n1])
    pair_interfaces(dom, "a", "b", ContactPenalty(K_n=K_n), name="c", integration=integracion)
    for nd in dom.nodes.values():
        if abs(nd.coordinates[1]) < 1e-9:
            nd.fix_dof("uy", 0.0)
    for y in (0.0, H1 + H2):
        eje = [nd for nd in dom.nodes.values()
               if abs(nd.coordinates[0] - Lx / 2) < 1e-9 and abs(nd.coordinates[1] - y) < 1e-9]
        eje[0].fix_dof("ux", 0.0)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)

    def presion(x, nn):
        return -p * nn

    for m, lados in ((m1, ("xmin", "xmax")), (m2, ("xmin", "xmax", "ymax"))):
        for lado in lados:
            for el, k in m.boundary[lado]:
                F[el.get_global_dof_indices()] += el.compute_edge_traction(k, presion)
    asm = Assembler(dom)
    solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=1), F_applied=F)
    ifs = _interfaces(dom)
    assert {type(e).__name__ for e in ifs} == {"Interface2D6"}
    tn = np.array([t[0] for e in ifs for t in e.state.stresses])
    return np.max(np.abs(tn + p)) / p


@pytest.mark.parametrize("integracion", ["nodal", "gauss"])
def test_lado_curvo_2d6_transmite_presion_uniforme_con_error_de_penalizacion(integracion):
    """Interfaz parabólica (flecha 0.10 en 2.0): el error es el de la
    penalización, ∝ 1/K_n (3.6e-4, 3.6e-6, 3.6e-8 con la nodal, *medido* por
    el teórico)."""
    errores = [_parcela_2d6_curva(integracion, K_n) for K_n in (4.0e5, 4.0e7, 4.0e9)]
    assert errores[2] < 1e-6
    assert errores[1] / errores[0] < 0.05 and errores[2] / errores[1] < 0.05


# --- simetría del sistema -----------------------------------------------------------

@pytest.mark.parametrize("ley, simetrico", [
    (ContactPenalty(K_n=1.0e6), True),
    (ContactPenalty(K_n=1.0e6, mu=0.3, K_t=1.0e6), False),
])
def test_con_friccion_el_modelo_pide_el_solver_no_simetrico(ley, simetrico):
    """La tangente de Coulomb no es simétrica (spec ContactPenalty §10): el
    despachador del ADR 0003 la ve a través del material del elemento."""
    m = structured_mesh("Quad4", lines=[np.linspace(0, 1, 2), np.linspace(0, 1, 3)],
                        material=Elastic2D(E=1.0e3, nu=0.2, hypothesis="plane_strain"))
    dom = m.domain
    _grupo(dom, "g", [i for i, nd in dom.nodes.items() if abs(nd.coordinates[1] - 0.5) < 1e-9])
    from solidum.mesh.interfaces import split_interfaces
    split_interfaces(dom, "g", ley, name="c")
    assert domain_is_symmetric(dom) is simetrico
