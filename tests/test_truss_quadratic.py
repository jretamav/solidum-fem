"""Armaduras de segundo orden ``Truss2D3`` y ``Truss3D3`` (specs del mismo
nombre): casos de aceptación.

Las soluciones manufacturadas generales (campo exacto con u* cúbica, órdenes
3 y 2) están en ``tests/validation/test_mms_frames.py``; aquí van las
propiedades propias del elemento cuadrático: matrices en forma cerrada,
casos analíticos de extremo a extremo, la exactitud de Tong en los extremos
frente a la superconvergencia del nodo intermedio, frecuencias axiales y
material no lineal por punto de Gauss. Parámetros físicos distintos de 1 en
todas las pruebas (E, A, ρ, L).
"""
import logging

import numpy as np
import pytest

from solidum import Assembler, Elastic1D, Elastoplastic1D, NonlinearSolver, run, run_modal
from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.elements.truss_quadratic import Truss2D3, Truss3D3, quadratic_shape_functions
from solidum.mesh import structured_mesh
from solidum.verification import ManufacturedFrame, solve_manufactured

E, A, RHO = 2.1e11, 3.0e-4, 7850.0
K_HAT = np.array([[7.0, 1.0, -8.0], [1.0, 7.0, -8.0], [-8.0, -8.0, 16.0]]) / 3.0
M_HAT = np.array([[4.0, -1.0, 2.0], [-1.0, 4.0, 2.0], [2.0, 2.0, 16.0]]) / 30.0
P0_2D, P1_2D = np.array([0.4, -0.7]), np.array([2.9, 1.1])
P0_3D, P1_3D = np.array([0.4, -0.7, 0.3]), np.array([2.9, 1.1, -1.4])


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.WARNING)
    yield
    logging.disable(logging.NOTSET)


def cerca(actual, esperado, rtol=1e-12):
    esperado = np.asarray(esperado, dtype=float)
    np.testing.assert_allclose(actual, esperado, rtol=rtol,
                               atol=rtol * max(np.abs(esperado).max(), 1e-300))


def elemento(cls, p0, p1, material=None, area=A):
    nodos = [Node(1, p0), Node(2, p1), Node(3, 0.5 * (p0 + p1))]
    return cls(1, nodos, material or Elastic1D(E, density=RHO), area)


def barra(cls, p0, p1, n, material=None):
    """Barra recta de ``n`` elementos de ``p0`` a ``p1``; nodos de los
    extremos de elemento ``1..n+1`` y los intermedios a continuación."""
    material = material or Elastic1D(E, density=RHO)
    d = Domain()
    ends = [d.add_node(k + 1, list(p0 + (p1 - p0) * k / n)) for k in range(n + 1)]
    mids = [d.add_node(n + 2 + k, list(p0 + (p1 - p0) * (k + 0.5) / n)) for k in range(n)]
    for k in range(n):
        d.add_element(cls(k + 1, [ends[k], ends[k + 1], mids[k]], material, A))
    return d, ends, mids


def axial(p0, p1):
    return (p1 - p0) / np.linalg.norm(p1 - p0)


# --- construcción y matrices ------------------------------------------------

@pytest.mark.parametrize("cls,p0,p1", [(Truss2D3, P0_2D, P1_2D), (Truss3D3, P0_3D, P1_3D)])
def test_rigidez_y_masa_en_forma_cerrada(cls, p0, p1):
    el = elemento(cls, p0, p1)
    L, e = np.linalg.norm(p1 - p0), axial(p0, p1)
    K, F = el.compute_element_state(np.zeros(3 * len(e)))
    cerca(K, np.kron(E * A / L * K_HAT, np.outer(e, e)))
    cerca(F, np.zeros_like(F), rtol=0.0)
    cerca(el.compute_mass_matrix(), np.kron(RHO * A * L * M_HAT, np.eye(len(e))))
    cerca(el.compute_mass_matrix("lumped"),
          np.kron(np.diag(RHO * A * L * np.array([1.0, 1.0, 4.0]) / 6.0), np.eye(len(e))))
    # Rango 2 (dos modos de deformación axial) y masa total ρAL por dirección.
    w = np.linalg.eigvalsh(K)
    assert np.sum(w > 1e-9 * w.max()) == 2
    assert np.sum(np.abs(w) <= 1e-9 * w.max()) == 3 * len(e) - 2
    M = el.compute_mass_matrix()
    for k in range(len(e)):
        t = np.zeros(3 * len(e)); t[k::len(e)] = 1.0
        assert t @ M @ t == pytest.approx(RHO * A * L, rel=1e-12)


def test_3d_en_el_plano_igual_que_2d():
    el2 = elemento(Truss2D3, P0_2D, P1_2D)
    el3 = elemento(Truss3D3, np.append(P0_2D, 0.0), np.append(P1_2D, 0.0))
    K2, _ = el2.compute_element_state(np.zeros(6))
    K3, _ = el3.compute_element_state(np.zeros(9))
    keep = [0, 1, 3, 4, 6, 7]
    cerca(K3[np.ix_(keep, keep)], K2)
    cerca(el3.compute_mass_matrix()[np.ix_(keep, keep)], el2.compute_mass_matrix())


def test_rechazos():
    mat = Elastic1D(E)
    a, b = Node(1, [0.0, 0.0]), Node(2, [2.0, 0.0])
    with pytest.raises(ValueError, match="3 nodos"):
        Truss2D3(1, [a, b], mat, A)
    with pytest.raises(ValueError, match="punto medio"):
        Truss2D3(1, [a, b, Node(3, [1.0, 1e-3])], mat, A)
    with pytest.raises(ValueError, match="punto medio"):
        Truss2D3(1, [a, b, Node(3, [0.7, 0.0])], mat, A)
    with pytest.raises(ValueError, match="área"):
        Truss2D3(1, [a, b, Node(3, [1.0, 0.0])], mat, 0.0)


# --- casos analíticos de extremo a extremo ----------------------------------

@pytest.mark.parametrize("cls,p0,p1", [(Truss2D3, P0_2D, P1_2D), (Truss3D3, P0_3D, P1_3D)])
@pytest.mark.parametrize("n", [1, 3])
def test_carga_puntual_en_el_extremo(cls, p0, p1, n):
    """``u(L) = F·L/(E·A)``; ``N = F`` en los puntos de Gauss y en los
    extremos. Los nodos se mueven sólo a lo largo del eje (la barra no tiene
    rigidez transversal)."""
    d, ends, mids = barra(cls, p0, p1, n)
    e, L, F = axial(p0, p1), np.linalg.norm(p1 - p0), 3.7e4
    dofs = cls.DOF_NAMES
    for name in dofs:
        ends[0].fix_dof(name, 0.0)
    k = int(np.argmax(np.abs(e)))
    for node in ends[1:] + mids:
        for j, name in enumerate(dofs):
            if j != k:
                d.add_linear_constraint((node.id, name), [(node.id, dofs[k])], [e[j] / e[k]])
    d.generate_equation_numbers()
    Fv = np.zeros(d.total_dofs)
    for j, name in enumerate(dofs):
        Fv[ends[-1].dofs[name]] = F * e[j]
    r = run(d, F_applied=Fv)
    u_end = np.array([r.U[ends[-1].dofs[name]] for name in dofs])
    cerca(u_end, F * L / (E * A) * e, rtol=1e-10)
    for el in d.elements.values():
        cerca(r.element_forces[el.id].components["N"], [F, F], rtol=1e-10)
        cerca(el.compute_internal_forces(r.U)["axial_force"], [F, F], rtol=1e-10)


def test_carga_axial_uniforme_exacta_con_un_elemento():
    """Empotrada en s = 0 y libre en s = L con carga de cuerpo axial ``b``:
    ``u(s) = b(Ls − s²/2)/E`` en todo punto y ``N(s) = bA(L − s)``."""
    p0, p1 = P0_2D, P1_2D
    e, L, b = axial(p0, p1), np.linalg.norm(p1 - p0), 2.3e5
    d, ends, mids = barra(Truss2D3, p0, p1, 1)
    for name in ("ux", "uy"):
        ends[0].fix_dof(name, 0.0)
    for node in (ends[1], mids[0]):
        d.add_linear_constraint((node.id, "uy"), [(node.id, "ux")], [e[1] / e[0]])
    d.generate_equation_numbers()
    asm = Assembler(d)
    r = run(d, assembler=asm, F_applied=asm.assemble_body_load(b * e))
    el = d.elements[1]
    s = np.linspace(0.0, 1.0, 7)
    campo, eps = el.line_field(el.get_local_displacements(r.U), s)
    cerca(campo[:, 0], b * (L * s * L - (s * L) ** 2 / 2) / E, rtol=1e-10)
    cerca(eps[:, 0], b * (L - s * L) / E, rtol=1e-10)
    cerca(r.element_forces[1].components["N"], [b * A * L, 0.0], rtol=1e-10)


def test_fuerzas_desde_nodales_iguales_a_las_del_campo():
    """``internal_forces_from_nodal`` (perturbación lineal, ADR 0031) invierte
    las proyecciones de ``F_int`` en los extremos: da el mismo ``N`` lineal."""
    el = elemento(Truss3D3, P0_3D, P1_3D)
    rng = np.random.default_rng(3)
    for _ in range(3):
        u = rng.standard_normal(9) * 1e-4
        _, F_int = el.compute_element_state(u)
        el.state.vars_trial = [None, None]
        U = u                                              # vector del elemento
        el.get_local_displacements = lambda Ug: Ug         # sin dominio
        cerca(el.internal_forces_from_nodal(F_int, U).components["N"],
              el.internal_forces(U).components["N"], rtol=1e-10)


# --- nodos: Tong en los extremos, superconvergencia en el intermedio ---------

def _errores_nodales(u, n):
    rec = lambda X: P0_2D + X[0] * (P1_2D - P0_2D)                  # noqa: E731
    m = structured_mesh("Truss2D3", (n,), mapping=rec, material=Elastic1D(E), A=A)
    sol = ManufacturedFrame.from_element({"u": u}, next(iter(m.domain.elements.values())),
                                         P0_2D)
    U = solve_manufactured(m, sol, ("xmin",))
    ends, mids = set(), set()
    for el in m.domain.elements.values():
        ends |= {el.nodes[0].id, el.nodes[1].id}
        mids.add(el.nodes[2].id)

    def err(ids):
        return max(float(np.max(np.abs(
            U[[m.domain.nodes[i].dofs[k] for k in ("ux", "uy")]]
            - sol._dof_values(m.domain.nodes[i].coordinates)))) for i in ids)
    return err(ends), err(mids), float(np.max(np.abs(U)))


def test_extremos_exactos_con_carga_polinomica():
    """Tong (1969): la función de Green (lineal con quiebre en el nodo) está
    en el espacio en los extremos de elemento, no en el nodo intermedio."""
    for n in (1, 3):
        e_ext, e_mid, escala = _errores_nodales("1e-3*(1 + s - 0.3*s**2 + 0.1*s**4)", n)
        assert e_ext < 1e-12 * escala
        assert e_mid > 1e-9 * escala           # el intermedio no es exacto
        e_ext, e_mid, escala = _errores_nodales("1e-3*(1 + s - 0.3*s**2 + 0.1*s**3)", n)
        assert max(e_ext, e_mid) < 1e-12 * escala


def test_nodo_intermedio_superconvergente():
    """En el nodo intermedio el error converge con orden 4 (medido 3.95,
    3.99, 4.00)."""
    errores = [_errores_nodales("1e-3*sin(1.3*s + 0.2)", n)[1] for n in (4, 8, 16, 32)]
    ordenes = np.log2(np.array(errores[:-1]) / np.array(errores[1:]))
    assert np.all(np.diff(errores) < 0.0)
    assert ordenes[-1] >= 4.0 - 0.1


# --- dinámica y no linealidad --------------------------------------------

def test_primera_frecuencia_axial():
    """Empotrada-libre, 4 elementos, masa consistente: ω₁ = (π/2)·√(E/ρ)/L
    con error ≤ 2e-5 (medido 1.64e-5)."""
    p0, p1 = np.array([0.0, 0.0]), np.array([3.2, 0.0])
    L = 3.2
    d, ends, mids = barra(Truss2D3, p0, p1, 4)
    ends[0].fix_dof("ux", 0.0)
    for node in ends + mids:
        node.fix_dof("uy", 0.0)
    d.generate_equation_numbers()
    r = run_modal(d, n_modes=1)
    w1 = 0.5 * np.pi * np.sqrt(E / RHO) / L
    assert abs(r.frequencies_rad[0] / w1 - 1.0) <= 2e-5
    assert r.frequencies_rad[0] > w1                        # la consistente sobreestima


def test_fluencia_en_los_dos_puntos_de_gauss():
    """Barra elastoplástica con endurecimiento lineal bajo desplazamiento
    impuesto en el extremo: deformación uniforme ``ε = δ/L`` y, en los dos
    puntos de Gauss, ``σ = σ_y + EH/(E + H)·(ε − σ_y/E)``; el estado de cada
    punto se confirma al converger. (Con plasticidad perfecta, la barra
    totalmente fluida es un mecanismo: la tangente es nula.)"""
    sy, H = 3.5e8, 0.05 * E
    mat = Elastoplastic1D(E=E, sigma_y=sy, H=H)
    p0, p1 = P0_2D, P1_2D
    e, L = axial(p0, p1), np.linalg.norm(p1 - p0)
    d, ends, mids = barra(Truss2D3, p0, p1, 2, material=mat)
    for name in ("ux", "uy"):
        ends[0].fix_dof(name, 0.0)
    delta = 3.0 * sy / E * L
    for j, name in enumerate(("ux", "uy")):
        ends[-1].fix_dof(name, delta * e[j])
    for node in ends[1:-1] + mids:
        d.add_linear_constraint((node.id, "uy"), [(node.id, "ux")], [e[1] / e[0]])
    d.generate_equation_numbers()
    asm = Assembler(d)
    r = run(d, assembler=asm, solver=NonlinearSolver(asm, num_steps=6))
    eps = delta / L
    N = A * (sy + E * H / (E + H) * (eps - sy / E))
    for el in d.elements.values():
        cerca(np.array(el.state.stresses) * A, [N, N], rtol=1e-8)
        cerca(r.element_forces[el.id].components["N"], [N, N], rtol=1e-8)
        for v in el.state.vars:
            assert v["eps_p"] == pytest.approx(eps - N / (A * E), rel=1e-8)


def test_funciones_de_forma():
    xi = np.linspace(0.0, 1.0, 11)
    N, dN = quadratic_shape_functions(xi)
    cerca(N.sum(axis=1), np.ones_like(xi))
    np.testing.assert_allclose(dN.sum(axis=1), 0.0, atol=1e-14)
    cerca(quadratic_shape_functions([0.0, 1.0, 0.5])[0], np.eye(3))
