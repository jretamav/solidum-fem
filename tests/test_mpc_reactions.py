"""Reacciones y fuerzas de restricción con restricciones multipunto (deuda #32).

Convención decidida por el usuario el 2026-09-30 (ADR 0004 §3):

- La reacción es la fuerza que un **apoyo** ejerce sobre la estructura. Un
  apoyo Dirichlet reporta su reacción completa, incluida la parte que le llega
  por nodos atados a él.
- Una MPC que impide un movimiento de sólido rígido (un rodillo oblicuo
  ``u_y = tanθ·u_x``, una relación mixta ``u_s = 2·u_m`` entre nodos
  distintos) cuenta como apoyo: su fuerza entra en las reacciones de los DOF
  de su esclavo y de sus maestros.
- Las MPC de **unión** (atado de mallas, periodicidad, unión rígida) transmiten
  una fuerza de restricción entre partes de la estructura: no son reacciones;
  van a ``SolveResult.constraint_forces``. Una MPC es de unión si una traslación
  rígida de toda la estructura, o un campo escalar uniforme, la cumple.

Formulación: con ``C`` las restricciones declaradas sin cerrar y
``r = F_int − F_applied``, los multiplicadores cumplen ``Cᵀ·μ = r``; las
reacciones son ``Σ μ_d·e_d + Σ_apoyo μ_k·C_kᵀ`` y las fuerzas de restricción
``Σ_unión μ_k·C_kᵀ``, de modo que ``F_int − F_applied = R + constraint_forces``.

Oráculo principal: atar dos nodos coincidentes con ``u_s = u_m`` en todos sus
DOF equivale a fundirlos en uno, así que las reacciones de los apoyos reales
del modelo atado deben coincidir con las del fundido. Con el código anterior
(``R = F_int − F`` en los DOF esclavos y ``reactions_by_node`` sólo en los
nodos con ``boundary_conditions``) un apoyo en un maestro perdía la parte que
le llega por el esclavo, un rodillo oblicuo no figuraba en las reacciones y
una unión dejaba fuerza en ``R`` con ``ΣR + ΣF ≠ 0``.
"""
from __future__ import annotations

import logging

import numpy as np
import pytest

import solidum
from solidum.bc.constraints import ConstraintSet
from solidum.core.domain import Domain
from solidum.elements.frame.euler import Frame2DEuler
from solidum.elements.solid_2d.quad4 import Quad4
from solidum.elements.solid_3d import Hex8
from solidum.elements.thermal.quad4_thermal import Quad4Thermal
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.materials.von_mises_2d import VonMises2D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import LinearSolver, NonlinearSolver

DUP = 1000            # desplazamiento de los ids de los nodos duplicados
TRASLACIONES = ("ux", "uy", "uz")


@pytest.fixture(autouse=True)
def _silencio():
    logging.disable(logging.CRITICAL)
    yield
    logging.disable(logging.NOTSET)


# ----------------------------------------------------------------------
# Utilidades
# ----------------------------------------------------------------------

def _resolver(dom, F, solver=None):
    asm = Assembler(dom)
    return solidum.run(dom, assembler=asm, solver=None if solver is None else solver(asm),
                       F_applied=F), asm


def _suma_por_dof(dom, v):
    """``{dof_name: Σ v}`` sobre todos los nodos."""
    tot: dict[str, float] = {}
    for node in dom.nodes.values():
        for d, k in node.dofs.items():
            tot[d] = tot.get(d, 0.0) + float(v[k])
    return tot


def _dofs_de_la_vista(dom, vista):
    return np.array(sorted(dom.nodes[n].dofs[d] for n, dd in vista.items() for d in dd),
                    dtype=np.int64)


def _equilibrio_global(dom, res, escala, rtol=1e-9, campos=TRASLACIONES + ("T",)):
    """ΣR + ΣF = 0 y Σ constraint_forces = 0 en cada traslación (o campo escalar)."""
    sr = _suma_por_dof(dom, res.R)
    sf = _suma_por_dof(dom, res.F_applied)
    sc = _suma_por_dof(dom, res.constraint_forces)
    for d in sr:
        if d in campos:
            assert abs(sr[d] + sf[d]) <= rtol * escala, (d, sr[d], sf[d])
            assert abs(sc[d]) <= rtol * escala, (d, sc[d])


def _comparar_con_fundido(dom_t, res_t, res_f, rtol):
    """Las reacciones del modelo atado coinciden con las del fundido, y ``R``
    del atado sólo tiene valores en los DOF de sus apoyos."""
    escala = max(np.abs(res_f.R).max(), np.abs(res_f.F_applied).max())
    assert res_t.reactions_by_node.keys() == res_f.reactions_by_node.keys()
    for nid, por_dof in res_f.reactions_by_node.items():
        assert res_t.reactions_by_node[nid].keys() == por_dof.keys(), nid
        for d, v in por_dof.items():
            assert res_t.reactions_by_node[nid][d] == pytest.approx(v, abs=rtol * escala), (nid, d)
    fuera = np.ones(dom_t.total_dofs, dtype=bool)
    fuera[_dofs_de_la_vista(dom_t, res_t.reactions_by_node)] = False
    assert np.abs(res_t.R[fuera]).max() <= rtol * escala
    _equilibrio_global(dom_t, res_t, escala, rtol)


def _malla_2d(dom, nx, ny, lx, ly, atado, ic, factory):
    """Malla estructurada de cuadriláteros sobre [0, lx]×[0, ly]. Con ``atado``
    la columna ``i = ic`` se duplica para los elementos de su derecha (ids
    ``DUP + id``). Devuelve ``(nid, dup)``."""
    nid, dup = {}, {}
    for j in range(ny + 1):
        for i in range(nx + 1):
            nid[i, j] = j * (nx + 1) + i + 1
            dom.add_node(nid[i, j], [lx * i / nx, ly * j / ny])
    if atado:
        for j in range(ny + 1):
            dup[ic, j] = DUP + nid[ic, j]
            dom.add_node(dup[ic, j], [lx * ic / nx, ly * j / ny])
    eid = 0
    for j in range(ny):
        for i in range(nx):
            eid += 1
            esquinas = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
            ids = [dup[c] if (atado and c[0] == ic and i >= ic) else nid[c] for c in esquinas]
            dom.add_element(factory(eid, [dom.nodes[k] for k in ids]))
    return nid, dup


def _atar(dom, dup, nid, dofs):
    for key, s in dup.items():
        for d in dofs:
            dom.add_linear_constraint((s, d), [(nid[key], d)], [1.0])


# ----------------------------------------------------------------------
# Modelos con versión atada (MPC u_s = u_m) y fundida (malla conforme)
# ----------------------------------------------------------------------

def _placa_quad4(atado, material="elastico"):
    """Placa Quad4 [0, 1]×[0, 0.5], empotrada en x = 0, con un apoyo ``uy`` en
    el nodo inferior de la interfaz (maestro en el modelo atado) y carga en el
    borde libre."""
    if material == "elastico":
        mat = Elastic2D(E=200.0e9, nu=0.3, hypothesis="plane_stress")
    else:
        mat = VonMises2D(E=200.0e9, nu=0.3, sigma_y=250.0e6, H=2.0e9,
                         hypothesis="plane_stress")
    dom = Domain()
    nx, ny, ic = 4, 2, 2
    nid, dup = _malla_2d(dom, nx, ny, 1.0, 0.5, atado, ic,
                         lambda e, n: Quad4(e, n, mat, thickness=0.2))
    for j in range(ny + 1):
        dom.nodes[nid[0, j]].fix_dof("ux", 0.0)
        dom.nodes[nid[0, j]].fix_dof("uy", 0.0)
    dom.nodes[nid[ic, 0]].fix_dof("uy", 0.0)
    if atado:
        _atar(dom, dup, nid, ("ux", "uy"))
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    for j in range(ny + 1):
        F[dom.nodes[nid[nx, j]].dofs["uy"]] = -4.0e6
        F[dom.nodes[nid[nx, j]].dofs["ux"]] = 1.0e6
    return dom, F


def _armadura(atado):
    """Armadura triangular: nodos 1 (0, 0), 2 (4, 0) y 3 (2, 2); articulación
    en 1 y rodillo ``uy`` en 2. En el modelo atado la barra 2-3 llega a un
    nodo propio (``DUP + 2``) atado al 2."""
    dom = Domain()
    mat = Elastic1D(E=200.0e9)
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [4.0, 0.0])
    n3 = dom.add_node(3, [2.0, 2.0])
    n2b = dom.add_node(DUP + 2, [4.0, 0.0]) if atado else n2
    dom.add_element(Truss2D(1, [n1, n2], mat, A=1.0e-3))
    dom.add_element(Truss2D(2, [n1, n3], mat, A=1.0e-3))
    dom.add_element(Truss2D(3, [n2b, n3], mat, A=1.0e-3))
    n1.fix_dof("ux", 0.0)
    n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    if atado:
        for d in ("ux", "uy"):
            dom.add_linear_constraint((DUP + 2, d), [(2, d)], [1.0])
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[n3.dofs["ux"]] = 1.0e4
    F[n3.dofs["uy"]] = -5.0e4
    return dom, F


def _marco(atado):
    """Viga continua Frame2DEuler de dos vanos de 3 m (nodos cada 1.5 m):
    articulación en x = 0, apoyos ``uy`` en x = 3 y x = 6. En el modelo atado
    el nodo del apoyo central se duplica para el vano derecho y se ata en
    ``ux``, ``uy`` y ``rz``."""
    dom = Domain()
    mat = Elastic1D(E=30.0e9)
    nodos = [dom.add_node(i + 1, [1.5 * i, 0.0]) for i in range(5)]
    dup = dom.add_node(DUP + 3, [3.0, 0.0]) if atado else None
    for e in range(4):
        ni, nj = nodos[e], nodos[e + 1]
        if atado and e == 2:
            ni = dup
        dom.add_element(Frame2DEuler(e + 1, [ni, nj], mat, A=0.06, I=4.5e-4))
    nodos[0].fix_dof("ux", 0.0)
    nodos[0].fix_dof("uy", 0.0)
    nodos[2].fix_dof("uy", 0.0)
    nodos[4].fix_dof("uy", 0.0)
    if atado:
        for d in ("ux", "uy", "rz"):
            dom.add_linear_constraint((DUP + 3, d), [(3, d)], [1.0])
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[nodos[1].dofs["uy"]] = -2.0e4
    F[nodos[3].dofs["uy"]] = -1.0e4
    F[nodos[3].dofs["ux"]] = 5.0e3
    F[nodos[3].dofs["rz"]] = 3.0e3
    return dom, F


def _bloque_hex8(atado):
    """Bloque Hex8 [0, 2]×[0, 1]×[0, 1] de 2×1×1 elementos, empotrado en
    x = 0, con apoyos ``uz`` en los nodos inferiores del plano x = 1 (maestros
    en el modelo atado, donde los elementos de x > 1 tienen nodos propios en
    ese plano)."""
    dom = Domain()
    mat = Elastic3D(E=70.0e9, nu=0.33)
    nx, ny, nz, ic = 2, 1, 1, 1
    nid, dup = {}, {}
    k = 0
    for kk in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                k += 1
                nid[i, j, kk] = k
                dom.add_node(k, [float(i), float(j), float(kk)])
    if atado:
        for kk in range(nz + 1):
            for j in range(ny + 1):
                dup[ic, j, kk] = DUP + nid[ic, j, kk]
                dom.add_node(dup[ic, j, kk], [float(ic), float(j), float(kk)])
    eid = 0
    for kk in range(nz):
        for j in range(ny):
            for i in range(nx):
                eid += 1
                c = [(i, j, kk), (i + 1, j, kk), (i + 1, j + 1, kk), (i, j + 1, kk),
                     (i, j, kk + 1), (i + 1, j, kk + 1), (i + 1, j + 1, kk + 1),
                     (i, j + 1, kk + 1)]
                ids = [dup[q] if (atado and q[0] == ic and i >= ic) else nid[q] for q in c]
                dom.add_element(Hex8(eid, [dom.nodes[q] for q in ids], mat))
    for j in range(ny + 1):
        for kk in range(nz + 1):
            for d in TRASLACIONES:
                dom.nodes[nid[0, j, kk]].fix_dof(d, 0.0)
    for j in range(ny + 1):
        dom.nodes[nid[ic, j, 0]].fix_dof("uz", 0.0)
    if atado:
        _atar(dom, dup, nid, TRASLACIONES)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    for j in range(ny + 1):
        for kk in range(nz + 1):
            n = dom.nodes[nid[nx, j, kk]]
            F[n.dofs["uz"]] = -1.0e6
            F[n.dofs["uy"]] = 2.0e5
    return dom, F


def _placa_termica(atado):
    """Placa Quad4Thermal con ``T = 100`` en x = 0 y ``T = 50`` en el nodo
    inferior de la columna atada (maestro en el modelo atado); calor entrante
    en el borde derecho. La «reacción» es el flujo nodal en los DOF con
    temperatura impuesta."""
    mat = ThermalConduction(k=45.0, c=460.0, density=7850.0, dim=2)
    dom = Domain()
    nx, ny, ic = 6, 3, 3
    nid, dup = _malla_2d(dom, nx, ny, 0.3, 0.15, atado, ic,
                         lambda e, n: Quad4Thermal(e, n, mat, thickness=0.02))
    for j in range(ny + 1):
        dom.nodes[nid[0, j]].fix_dof("T", 100.0)
    dom.nodes[nid[ic, 0]].fix_dof("T", 50.0)
    if atado:
        _atar(dom, dup, nid, ("T",))
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    for j in range(ny + 1):
        F[dom.nodes[nid[nx, j]].dofs["T"]] = 50.0
    return dom, F


def _estricto():
    return ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10)


ORACULOS = {
    "quad4_lineal": (lambda a: _placa_quad4(a), None, 1e-9),
    "quad4_von_mises_newton": (
        lambda a: _placa_quad4(a, "von_mises"),
        lambda asm: NonlinearSolver(asm, _estricto(), num_steps=4),
        1e-7,
    ),
    "armadura_truss2d": (_armadura, None, 1e-9),
    "marco_frame2d_euler": (_marco, None, 1e-9),
    "solido_hex8": (_bloque_hex8, None, 1e-9),
    "termico_quad4": (_placa_termica, None, 1e-9),
}


# ----------------------------------------------------------------------
# 1. Oráculo del nodo fundido
# ----------------------------------------------------------------------

@pytest.mark.parametrize("modelo", list(ORACULOS))
def test_oraculo_del_nodo_fundido(modelo):
    build, solver, rtol = ORACULOS[modelo]
    dom_f, F_f = build(False)
    res_f, _ = _resolver(dom_f, F_f, solver)
    dom_t, F_t = build(True)
    res_t, _ = _resolver(dom_t, F_t, solver)
    _comparar_con_fundido(dom_t, res_t, res_f, rtol)


def test_la_placa_von_mises_plastifica():
    """El caso no lineal del oráculo ejercita de verdad la plasticidad."""
    dom, F = _placa_quad4(True, "von_mises")
    _resolver(dom, F, lambda asm: NonlinearSolver(asm, _estricto(), num_steps=4))
    alpha = max(float(e.state.vars[g]["alpha"])
                for e in dom.elements.values() for g in range(4))
    assert alpha > 0.0


# ----------------------------------------------------------------------
# 2. Maestro con apoyo: la reacción es completa
# ----------------------------------------------------------------------

def test_maestro_con_apoyo_reaccion_completa():
    """Rodillo ``uy`` en el nodo 2 de la armadura, maestro de la barra 2-3
    atada: su reacción es la de la estática, 30 kN, con la fuerza que llega
    por el esclavo (antes daba sólo la de la barra 1-2)."""
    dom, F = _armadura(True)
    res, _ = _resolver(dom, F)
    # Momentos respecto al nodo 1: 4·R2y + 2·Fy − 2·Fx = 0.
    Fx, Fy = 1.0e4, -5.0e4
    R2y = (2.0 * Fx - 2.0 * Fy) / 4.0
    assert res.reactions_by_node[2]["uy"] == pytest.approx(R2y, rel=1e-10)
    assert res.reactions_by_node[1]["ux"] == pytest.approx(-Fx, rel=1e-10)
    assert res.reactions_by_node[1]["uy"] == pytest.approx(-Fy - R2y, rel=1e-10)
    # El esclavo no es un apoyo: ni figura en las reacciones ni tiene R.
    assert DUP + 2 not in res.reactions_by_node
    assert np.all(res.R[list(dom.nodes[DUP + 2].dofs.values())] == 0.0)


# ----------------------------------------------------------------------
# 3. Rodillo oblicuo impuesto con una MPC
# ----------------------------------------------------------------------

THETA = np.radians(30.0)


def _rodillo_oblicuo():
    """Placa Quad4 articulada en (0, 0) y con un rodillo oblicuo a 30° en
    (2, 0), impuesto con ``u_y = tan30°·u_x``; carga vertical y horizontal en
    el borde superior."""
    dom = Domain()
    mat = Elastic2D(E=30.0e9, nu=0.2)
    nx, ny = 4, 2
    nid, _ = _malla_2d(dom, nx, ny, 2.0, 1.0, False, 0,
                       lambda e, n: Quad4(e, n, mat, thickness=0.2))
    dom.nodes[nid[0, 0]].fix_dof("ux", 0.0)
    dom.nodes[nid[0, 0]].fix_dof("uy", 0.0)
    c = nid[nx, 0]
    dom.add_linear_constraint((c, "uy"), [(c, "ux")], [np.tan(THETA)])
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    for i in range(nx + 1):
        F[dom.nodes[nid[i, ny]].dofs["uy"]] = -1.0e5
    F[dom.nodes[nid[nx // 2, ny]].dofs["ux"]] = 3.0e4
    return dom, F, c


def test_rodillo_oblicuo_reaccion_normal_al_plano():
    dom, F, c = _rodillo_oblicuo()
    res, _ = _resolver(dom, F)
    assert set(res.reactions_by_node[c]) == {"ux", "uy"}
    Rx, Ry = res.reactions_by_node[c]["ux"], res.reactions_by_node[c]["uy"]
    escala = np.abs(F).max()
    # Sin componente en la dirección de rodadura (cos θ, sin θ).
    assert abs(Rx * np.cos(THETA) + Ry * np.sin(THETA)) <= 1e-9 * escala
    assert np.hypot(Rx, Ry) > 0.1 * escala
    _equilibrio_global(dom, res, escala)
    # Equilibrio de momentos respecto a la articulación (0, 0).
    M = 0.0
    for node in dom.nodes.values():
        x, y = node.coordinates[0], node.coordinates[1]
        fx = res.R[node.dofs["ux"]] + F[node.dofs["ux"]]
        fy = res.R[node.dofs["uy"]] + F[node.dofs["uy"]]
        M += x * fy - y * fx
    assert abs(M) <= 1e-9 * escala
    # Es un apoyo, no una unión.
    assert not np.any(res.constraint_forces)
    assert res.constraint_forces_by_node == {}


def test_rodillo_oblicuo_frente_a_la_estatica():
    """La estática de la placa: con la articulación en (0, 0) y el rodillo en
    (2, 0), el momento de las cargas respecto a la articulación fija la
    componente vertical de la reacción del rodillo, y la normal al plano
    fija la horizontal."""
    dom, F, c = _rodillo_oblicuo()
    res, _ = _resolver(dom, F)
    M_cargas = 0.0
    for node in dom.nodes.values():
        x, y = node.coordinates[0], node.coordinates[1]
        M_cargas += x * F[node.dofs["uy"]] - y * F[node.dofs["ux"]]
    # Reacción del rodillo: R·(−sin θ, cos θ) aplicada en (2, 0).
    Ry = -M_cargas / 2.0
    Rx = -Ry * np.tan(THETA)
    assert res.reactions_by_node[c]["uy"] == pytest.approx(Ry, rel=1e-10)
    assert res.reactions_by_node[c]["ux"] == pytest.approx(Rx, rel=1e-10)


# ----------------------------------------------------------------------
# 4. Unión: fuerzas de restricción fuera de R
# ----------------------------------------------------------------------

def test_union_fuerzas_de_restriccion_iguales_y_opuestas():
    dom, F = _placa_quad4(True)
    res, _ = _resolver(dom, F)
    escala = np.abs(F).max()
    assert res.constraint_forces_by_node, "la unión debe exponer sus fuerzas"
    for lc in dom.linear_constraints:
        (ns, d), ((nm, dm),) = lc["slave"], lc["masters"]
        fs = res.constraint_forces[dom.nodes[ns].dofs[d]]
        fm = res.constraint_forces[dom.nodes[nm].dofs[dm]]
        assert fs == pytest.approx(-fm, abs=1e-12 * escala)
        assert res.constraint_forces_by_node[ns][d] == fs
        assert res.constraint_forces_by_node[nm][dm] == fm
        # El esclavo no es apoyo: R nulo.
        assert res.R[dom.nodes[ns].dofs[d]] == 0.0
    # La interfaz transmite el cortante y el flector del voladizo.
    assert np.abs(res.constraint_forces).max() > 0.1 * escala
    _equilibrio_global(dom, res, escala)
    # Sólo hay reacciones en los apoyos Dirichlet.
    assert set(res.reactions_by_node) == {
        n.id for n in dom.nodes.values() if n.boundary_conditions}


def test_union_rigida_con_giro_es_union():
    """Unión rígida ``u_s = u_m − d_y·θ_m`` en un marco: una traslación la
    cumple (el giro no cambia con ella), así que es de unión."""
    dom = Domain()
    mat = Elastic1D(E=30.0e9)
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [2.0, 0.0])
    n3 = dom.add_node(3, [2.0, 0.5])     # satélite a 0.5 m por encima de n2
    n4 = dom.add_node(4, [4.0, 0.5])
    dom.add_element(Frame2DEuler(1, [n1, n2], mat, A=0.06, I=4.5e-4))
    dom.add_element(Frame2DEuler(2, [n3, n4], mat, A=0.06, I=4.5e-4))
    for d in ("ux", "uy", "rz"):
        n1.fix_dof(d, 0.0)
    dom.add_linear_constraint((3, "ux"), [(2, "ux"), (2, "rz")], [1.0, -0.5])
    dom.add_linear_constraint((3, "uy"), [(2, "uy")], [1.0])
    dom.add_linear_constraint((3, "rz"), [(2, "rz")], [1.0])
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[n4.dofs["uy"]] = -1.0e4
    F[n4.dofs["ux"]] = 2.0e3
    res, _ = _resolver(dom, F)
    assert set(res.reactions_by_node) == {1}
    assert set(res.constraint_forces_by_node) == {2, 3}
    _equilibrio_global(dom, res, np.abs(F).max())
    # El empotramiento reacciona contra toda la carga.
    assert res.reactions_by_node[1]["ux"] == pytest.approx(-2.0e3, rel=1e-10)
    assert res.reactions_by_node[1]["uy"] == pytest.approx(1.0e4, rel=1e-10)
    assert res.reactions_by_node[1]["rz"] == pytest.approx(
        -(4.0 * -1.0e4 - 0.5 * 2.0e3), rel=1e-10)


# ----------------------------------------------------------------------
# 5. Identidad F_int − F_applied = R + constraint_forces
# ----------------------------------------------------------------------

def _modelo_mixto():
    """Placa atada (unión, con apoyo en un maestro) más un rodillo oblicuo por
    MPC en la esquina inferior derecha (apoyo)."""
    mat = Elastic2D(E=200.0e9, nu=0.3, hypothesis="plane_stress")
    dom = Domain()
    nx, ny, ic = 4, 2, 2
    nid, dup = _malla_2d(dom, nx, ny, 1.0, 0.5, True, ic,
                         lambda e, n: Quad4(e, n, mat, thickness=0.2))
    for j in range(ny + 1):
        dom.nodes[nid[0, j]].fix_dof("ux", 0.0)
        dom.nodes[nid[0, j]].fix_dof("uy", 0.0)
    dom.nodes[nid[ic, 0]].fix_dof("uy", 0.0)
    _atar(dom, dup, nid, ("ux", "uy"))
    c = nid[nx, 0]
    dom.add_linear_constraint((c, "uy"), [(c, "ux")], [np.tan(np.radians(60.0))])
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    for j in range(ny + 1):
        F[dom.nodes[nid[nx, j]].dofs["uy"]] = -4.0e6
        F[dom.nodes[nid[nx, j]].dofs["ux"]] = 1.0e6
    return dom, F, c


@pytest.mark.parametrize("solver", [None, "newton"])
def test_identidad_residuo_igual_a_reacciones_mas_fuerzas_de_restriccion(solver):
    dom, F, c = _modelo_mixto()
    fab = None if solver is None else (lambda asm: NonlinearSolver(asm, _estricto(), num_steps=2))
    res, asm = _resolver(dom, F, fab)
    _, F_int = asm.assemble_non_linear_system(res.U)
    escala = np.abs(F).max()
    np.testing.assert_allclose(F_int - res.F_applied, res.R + res.constraint_forces,
                               rtol=0.0, atol=1e-8 * escala)
    # Cada vista contiene lo que dice: el rodillo en las reacciones, la
    # interfaz en las fuerzas de restricción.
    assert set(res.reactions_by_node[c]) == {"ux", "uy"}
    esclavos = {n for n in dom.nodes if n >= DUP}
    assert set(res.constraint_forces_by_node) == esclavos | {n - DUP for n in esclavos}
    _equilibrio_global(dom, res, escala, rtol=1e-8)


# ----------------------------------------------------------------------
# 6. Periodicidad con g ≠ 0 declarada por YAML
# ----------------------------------------------------------------------

YAML_PERIODICO = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [0.0, 1.0]}
  - {id: 4, coords: [1.0, 1.0]}
materials:
  - {id: 1, type: Elastic1D, E: 2.0e11}
elements:
  - {id: 1, type: Truss2D, material: 1, A: 1.0e-3, nodes: [1, 2]}
  - {id: 2, type: Truss2D, material: 1, A: 1.0e-3, nodes: [3, 4]}
boundary_conditions:
  - {node_id: 1, ux: 0.0, uy: 0.0}
  - {node_id: 3, ux: 0.0, uy: 0.0}
  - {node_id: 2, uy: 0.0}
  - {node_id: 4, uy: 0.0}
linear_constraints:
  # Periodicidad con salto: u_x del nodo 4 igual al del nodo 2 más g.
  - slave: {node: 4, dof: ux}
    masters:
      - {node: 2, dof: ux}
    coefficients: [1.0]
    g: __G__
point_loads:
  - {node_id: 2, ux: __P__}
solver:
  type: LinearSolver
"""


def test_periodicidad_con_salto_por_yaml(tmp_path):
    """Dos barras paralelas con ``u_x4 = u_x2 + g``: ``R`` sólo en los apoyos
    reales, y cada empotramiento reacciona contra el axil de su barra."""
    P, g, k = 1.0e5, 2.0e-4, 2.0e11 * 1.0e-3 / 1.0
    ruta = tmp_path / "periodico.yaml"
    ruta.write_text(YAML_PERIODICO.replace("__G__", repr(g)).replace("__P__", repr(P)),
                    encoding="utf-8")
    res = solidum.run_yaml(str(ruta))
    # Equilibrio: k·u2 + k·(u2 + g) = P.
    u2 = (P / k - g) / 2.0
    assert set(res.reactions_by_node) == {1, 2, 3, 4}
    assert set(res.reactions_by_node[2]) == {"uy"}
    assert set(res.reactions_by_node[4]) == {"uy"}
    assert res.reactions_by_node[1]["ux"] == pytest.approx(-k * u2, rel=1e-10)
    assert res.reactions_by_node[3]["ux"] == pytest.approx(-k * (u2 + g), rel=1e-10)
    assert (res.reactions_by_node[1]["ux"] + res.reactions_by_node[3]["ux"]
            == pytest.approx(-P, rel=1e-12))
    # La unión transmite al nodo 4 el axil de la barra 2, fuera de R.
    assert res.constraint_forces_by_node[4]["ux"] == pytest.approx(k * (u2 + g), rel=1e-10)
    assert res.constraint_forces_by_node[2]["ux"] == pytest.approx(-k * (u2 + g), rel=1e-10)
    assert abs(sum(res.R) + P) <= 1e-9 * P


# ----------------------------------------------------------------------
# 7. MPC mixta u_s = 2·u_m entre nodos distintos: cuenta como apoyo
# ----------------------------------------------------------------------

def test_mpc_mixta_cuenta_como_apoyo():
    """Dos barras colineales independientes, empotradas en sus extremos
    izquierdos, con ``u_x(s) = 2·u_x(m)`` entre sus extremos libres. Una
    traslación rígida no la cumple (1 − 2 ≠ 0): cuenta como apoyo y su fuerza
    entra en las reacciones de los dos nodos, con ``ΣR + ΣF = 0``."""
    dom = Domain()
    mat = Elastic1D(E=2.0e11)
    a = dom.add_node(1, [0.0, 0.0])
    m = dom.add_node(2, [1.0, 0.0])
    b = dom.add_node(3, [0.0, 1.0])
    s = dom.add_node(4, [1.0, 1.0])
    dom.add_element(Truss2D(1, [a, m], mat, A=1.0e-3))
    dom.add_element(Truss2D(2, [b, s], mat, A=1.0e-3))
    for n in (a, b):
        n.fix_dof("ux", 0.0)
        n.fix_dof("uy", 0.0)
    m.fix_dof("uy", 0.0)
    s.fix_dof("uy", 0.0)
    dom.add_linear_constraint((4, "ux"), [(2, "ux")], [2.0])
    dom.generate_equation_numbers(verbose=False)
    P = 1.0e5
    F = np.zeros(dom.total_dofs)
    F[m.dofs["ux"]] = P
    res, _ = _resolver(dom, F)
    assert "ux" in res.reactions_by_node[2] and "ux" in res.reactions_by_node[4]
    assert not np.any(res.constraint_forces)
    _equilibrio_global(dom, res, P)
    # k·u + 2·k·(2u) = P ⇒ u = P / (5k); la MPC toma −2·k·(2u) en el esclavo…
    k = 2.0e11 * 1.0e-3
    u = P / (5.0 * k)
    assert res.reactions_by_node[1]["ux"] == pytest.approx(-k * u, rel=1e-10)
    assert res.reactions_by_node[3]["ux"] == pytest.approx(-2.0 * k * u, rel=1e-10)
    # … y la MPC ejerce μ·(e_s − 2·e_m) con μ = k·2u (lo que empuja la barra 2).
    mu = 2.0 * k * u
    assert res.reactions_by_node[4]["ux"] == pytest.approx(mu, rel=1e-10)
    assert res.reactions_by_node[2]["ux"] == pytest.approx(-2.0 * mu, rel=1e-10)


# ----------------------------------------------------------------------
# 8. Sin MPC: R y reactions_by_node bit a bit como antes
# ----------------------------------------------------------------------

SIN_MPC = {
    "quad4": lambda: _placa_quad4(False),
    "armadura": lambda: _armadura(False),
    "marco": lambda: _marco(False),
    "hex8": lambda: _bloque_hex8(False),
    "termico": lambda: _placa_termica(False),
}


@pytest.mark.parametrize("modelo", list(SIN_MPC))
def test_sin_mpc_bit_a_bit(modelo):
    dom, F = SIN_MPC[modelo]()
    res, asm = _resolver(dom, F)
    _, F_int = asm.assemble_non_linear_system(res.U)
    bc = asm.constraint_set.slave_dofs
    R_antes = np.zeros(dom.total_dofs)
    R_antes[bc] = F_int[bc] - F[bc]
    np.testing.assert_array_equal(res.R, R_antes)
    esperado = {
        nid: {d: float(R_antes[n.dofs[d]]) for d in n.boundary_conditions}
        for nid, n in dom.nodes.items() if n.boundary_conditions
    }
    assert res.reactions_by_node == esperado
    np.testing.assert_array_equal(res.constraint_forces, np.zeros(dom.total_dofs))
    assert res.constraint_forces_by_node == {}


def test_solve_result_construido_a_mano_sigue_valiendo():
    """Los campos nuevos van al final y con valor por omisión."""
    r = solidum.SolveResult(U=np.zeros(2), F_applied=np.zeros(2), R=np.zeros(2))
    assert r.constraint_forces is None
    assert r.constraint_forces_by_node == {}


# ----------------------------------------------------------------------
# 9. Unitarios de ConstraintSet.multipliers
# ----------------------------------------------------------------------

def _C_sin_cerrar(cs, ndof):
    C = np.zeros((len(cs.declared), ndof))
    for k, c in enumerate(cs.declared):
        C[k, c.slave_dof] = 1.0
        for m, a in zip(c.masters, c.coefficients):
            C[k, m] -= a
    return C


def test_multipliers_sin_mpc_es_el_residuo():
    cs = ConstraintSet()
    cs.add_dirichlet(3, 0.0)
    cs.add_dirichlet(0, 1.5)
    r = np.array([0.1, 2.0, -3.0, 7.25, 1e-3])
    mu = cs.multipliers(r)
    assert list(mu) == [3, 0]          # orden de declaración
    assert mu == {3: 7.25, 0: 0.1}
    assert [c.slave_dof for c in cs.declared] == [3, 0]
    assert not cs.has_multipoint


def test_multipliers_cadena_hasta_un_dirichlet():
    """``u0 = 2·u1``, ``u1 = 3·u2 + 0.5``, ``u2`` prescrito: la fuerza del
    esclavo baja por la cadena hasta el apoyo."""
    cs = ConstraintSet()
    cs.add_linear(0, [1], [2.0])
    cs.add_linear(1, [2], [3.0], g=0.5)
    cs.add_dirichlet(2, 0.0)
    cs.build(4)
    r = np.array([1.0, -2.0, 4.0, 0.0])
    mu = cs.multipliers(r)
    assert mu[0] == pytest.approx(1.0)
    assert mu[1] == pytest.approx(-2.0 + 2.0 * 1.0)
    assert mu[2] == pytest.approx(4.0 + 3.0 * mu[1])
    C = _C_sin_cerrar(cs, 4)
    dofs = [c.slave_dof for c in cs.declared]
    lam = np.array([mu[d] for d in dofs])
    np.testing.assert_allclose((C.T @ lam)[dofs], r[dofs], atol=1e-14)
    # Las restricciones declaradas no las altera el cierre.
    assert cs.declared[0].masters == (1,)


def test_multipliers_varios_maestros_y_orden_de_declaracion_indiferente():
    """Una MPC con varios maestros, uno de ellos esclavo de otra MPC declarada
    después y otro con Dirichlet: ``Cᵀμ = r`` en los esclavos y, cuando ``r``
    es equilibrado en los libres, en todos los DOF."""
    cs = ConstraintSet()
    cs.add_linear(5, [1, 2, 3], [0.5, -1.5, 2.0], g=0.1)   # 1 es esclavo después
    cs.add_dirichlet(2, 0.0)
    cs.add_linear(1, [0, 4], [1.0, -0.25])
    cs.build(6)
    rng = np.random.default_rng(3)
    ndof = 6
    C = _C_sin_cerrar(cs, ndof)
    lam = rng.normal(size=len(cs.declared))
    r = C.T @ lam                     # residuo compatible: Cᵀλ en todos los DOF
    mu = cs.multipliers(r)
    got = np.array([mu[c.slave_dof] for c in cs.declared])
    np.testing.assert_allclose(got, lam, rtol=1e-13, atol=1e-13)
    dofs, arr = cs.multiplier_arrays(r)
    np.testing.assert_array_equal(dofs, [5, 2, 1])
    np.testing.assert_array_equal(arr, got)
    assert cs.has_multipoint
