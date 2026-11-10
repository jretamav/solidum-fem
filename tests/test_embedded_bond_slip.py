"""Refuerzo embebido con adherencia-deslizamiento (spec EmbeddedReinforcement
§12, casos ``bond_*``).

La extracción con el concreto fijo es la barra sobre resortes, con solución
cerrada s(x) = P/(E_s·A_s·λ)·cosh(λ(ℓ − x))/sinh(λℓ), λ = √(k·p/(E_s·A_s)).
"""
import logging
import textwrap

import numpy as np
import pytest

from solidum import (
    Assembler,
    Bond,
    BondSlipModelCode,
    BondSlipMultilinear,
    Elastic1D,
    Elastic2D,
    embed_reinforcement,
    load_yaml,
    run,
    run_modal,
)
from solidum.bc.bond_interface import BondSlipInterface
from solidum.materials.elastic_3d import Elastic3D
from solidum.mesh import structured_mesh

ES, DB = 200e9, 0.016
AB, PB = np.pi * DB ** 2 / 4, np.pi * DB
P, LB, H, T = 2.0e4, 1.0, 0.3, 0.3
EC = 30e9


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.WARNING)
    yield
    logging.disable(logging.NOTSET)


def linear_law(k, k_n=1e12):
    return BondSlipMultilinear(points=[[1.0, k * 1.0]], k_n=k_n)


def _plate(nx, ny=3, dim=2, element="Tri3"):
    if dim == 2:
        return structured_mesh(element, (nx, ny), size=(LB, H),
                               material=Elastic2D(E=EC, nu=0.2, hypothesis="plane_stress",
                                                  density=2400.0),
                               thickness=T)
    return structured_mesh(element, (nx, 2, 2), size=(LB, H, H),
                           material=Elastic3D(E=EC, nu=0.2, density=2400.0))


def pullout(k, element, nx, integration="nodal", dim=2, fixed_host=True, host="Tri3"):
    """Extracción: la barra se jala en x = 0 con P hacia −x."""
    mesh = _plate(nx, dim=dim, element=host)
    d = mesh.domain
    y = [0.0, H / 2] if dim == 2 else [0.0, H / 2, H / 2]
    q = [LB, H / 2] if dim == 2 else [LB, H / 2, H / 2]
    bar = embed_reinforcement(d, [y, q], AB, Elastic1D(E=ES, density=7850.0),
                              element=element, bond=Bond(linear_law(k), DB, integration=integration))
    names = ("ux", "uy", "uz")[:dim]
    if fixed_host:
        for nid, nd in d.nodes.items():
            if nid not in d.generated_nodes:
                for n in names:
                    nd.fix_dof(n)
    else:
        for nd in mesh.nodes_on("xmin"):
            nd.fix_dof("ux")
        nd0 = mesh.nodes_on("xmin")[0]
        for n in names[1:]:
            nd0.fix_dof(n)
        if dim == 3:
            mesh.nodes_on("xmin")[1].fix_dof("uz")
    d.generate_equation_numbers()
    F = np.zeros(d.total_dofs)
    F[d.nodes[bar.node_ids[0]].dofs["ux"]] = -P
    r = run(d, F_applied=F)
    return d, bar, r, mesh


def profile(d, bar, r):
    pts = [p for row in bar.bond_state(d, r.U) for p in row]
    x = np.array([p.position[0] for p in pts])
    o = np.argsort(x, kind="stable")
    return (x[o], -np.array([p.slip for p in pts])[o], -np.array([p.bond_stress for p in pts])[o],
            pts)


def exact_slip(x, k):
    lam = np.sqrt(k * PB / (ES * AB))
    return P / (ES * AB * lam) * np.cosh(lam * (LB - x)) / np.sinh(lam * LB)


# --- extracción contra la solución cerrada --------------------------------------

def test_extraccion_lineal_converge():
    """Medido: Truss2D3 nodal 2.7e-5 / 1.7e-6 / 1.1e-7 con nx = 5/10/20 (orden 4);
    Truss2D nodal 1.5e-2 / 3.9e-3 / 9.7e-4 (orden 2)."""
    k = 1e10
    for element, bound, order in (("Truss2D3", 5e-6, 3.9), ("Truss2D", 5e-3, 1.9)):
        errs = []
        for nx in (5, 10, 20):
            d, bar, r, _ = pullout(k, element, nx)
            x, s, _, _ = profile(d, bar, r)
            errs.append(np.max(np.abs(s - exact_slip(x, k))) / exact_slip(0.0, k))
        assert errs[1] < bound, (element, errs)
        rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
        assert np.all(rates >= order), (element, rates)


def test_sin_oscilacion_con_integracion_nodal():
    """λh = 3.5 > √6: la nodal da τ monótona; Gauss invierte su signo."""
    k = 1e12
    for integration, oscila in (("nodal", False), ("gauss", True)):
        d, bar, r, _ = pullout(k, "Truss2D", 5, integration=integration)
        _, _, tau, _ = profile(d, bar, r)
        dt = np.diff(tau)
        big = np.abs(dt) > 5e-3 * tau.max()
        flips = int(np.sum((dt[:-1] * dt[1:] < 0) & big[:-1] & big[1:]))
        assert (flips > 0 and tau.min() < -0.05 * tau.max()) == oscila, integration


def test_equilibrio_de_la_barra_y_del_concreto():
    d, bar, r, mesh = pullout(1e11, "Truss2D3", 10, fixed_host=False)
    _, _, _, pts = profile(d, bar, r)
    w = [wg for i in bar.interface_ids for wg in d.elements[i].weights]
    # fuerza de adherencia total sobre la barra = −carga aplicada
    total = sum(p.bond_force * wg for p, wg in zip(pts, w))
    assert total == pytest.approx(-P, rel=1e-10)
    normal = sum(abs(p.normal_traction[0]) * DB * wg for p, wg in zip(pts, w))
    assert normal < 1e-6 * P                                 # barra recta: sin fuerza normal
    reac = sum(r.reactions_by_node[n.id]["ux"] for n in mesh.nodes_on("xmin"))
    assert reac == pytest.approx(P, rel=1e-10)               # el concreto la devuelve


def test_traslacion_rigida_sin_fuerzas():
    d, bar, _, _ = pullout(1e11, "Truss2D3", 5)
    for iid in bar.interface_ids:
        el = d.elements[iid]
        u = np.tile([3.0e-4, -1.2e-4], len(el.nodes))
        _, F = el.compute_element_state(u)
        assert np.max(np.abs(F)) < 1e-12 * 1e11 * PB * 3e-4


def test_limite_adherencia_perfecta():
    """Con rigideces de adherencia 1e4 veces las de servicio, el sólido se
    mueve como con adherencia perfecta (MPC)."""
    def beam(bond):
        mesh = _plate(16, ny=4)
        d = mesh.domain
        host_ids = sorted(d.nodes)
        embed_reinforcement(d, [[0.0, 0.05], [LB, 0.05]], AB, Elastic1D(E=ES),
                            element="Truss2D3", bond=bond)
        for nd in mesh.nodes_on("xmin"):
            nd.fix_dof("ux"); nd.fix_dof("uy")
        d.generate_equation_numbers()
        F = np.zeros(d.total_dofs)
        for nd in mesh.nodes_on("xmax"):
            F[nd.dofs["uy"]] -= 2.0e4
        r = run(d, F_applied=F)
        return np.array([[r.U[d.nodes[i].dofs[n]] for n in ("ux", "uy")] for i in host_ids])
    U_perf = beam(None)
    U_bond = beam(Bond(linear_law(1e15, k_n=1e15), DB))
    np.testing.assert_allclose(U_bond, U_perf, rtol=0, atol=1e-4 * np.abs(U_perf).max())


def _single_interface(element="Truss2D3"):
    """Un anfitrión Quad4 fijo con una barra horizontal dentro y su interfaz."""
    mesh = structured_mesh("Quad4", (1, 1), size=(0.2, 0.1),
                           material=Elastic2D(E=EC, nu=0.2, hypothesis="plane_stress"),
                           thickness=T)
    d = mesh.domain
    law = BondSlipModelCode(tau_max=13.5e6, s1=1e-3, s2=3e-3, s3=10.5e-3, alpha=0.4,
                            tau_f=5e6, k0=1.8e11, k_n=3e12)
    bar = embed_reinforcement(d, [[0.02, 0.05], [0.18, 0.05]], AB, Elastic1D(E=ES),
                              element=element, bond=Bond(law, DB))
    return d, bar, law


def test_anclaje_con_deslizamiento_uniforme_sigue_la_ley():
    """Barra desplazada en bloque respecto del concreto fijo: la fuerza de
    extracción es τ(s)·p·ℓ con la τ de la ley en todo el ciclo de la Fig. 5.3b."""
    d, bar, law = _single_interface()
    el = d.elements[bar.interface_ids[0]]
    n_bar = len(d.elements[bar.element_ids[0]].nodes)
    ell = 0.16
    st = None
    for s in np.concatenate([np.linspace(0, 1, 21), np.linspace(1, -1, 41),
                             np.linspace(-1, 1.5, 51)])[1:] * 1e-3:
        u = np.zeros(2 * len(el.nodes))
        u[0:2 * n_bar:2] = s
        _, F = el.compute_element_state(u)
        el.commit_state()
        tau, _, st = law.compute_traction(np.array([0.0, s]), st)
        assert F[0:2 * n_bar:2].sum() == pytest.approx(tau[1] * PB * ell, rel=1e-10, abs=1e-6)
        assert -F[2 * n_bar::2].sum() == pytest.approx(tau[1] * PB * ell, rel=1e-10, abs=1e-6)


def test_quiebre_equilibrio_del_nodo():
    """Barra con un quiebre en tracción: en el nodo del quiebre la interfaz
    equilibra el cambio de dirección del axial, N1·t1 − N2·t2."""
    mesh = _plate(10, ny=4)
    d = mesh.domain
    pts = [[0.0, 0.06], [0.5, 0.15], [1.0, 0.06]]
    bar = embed_reinforcement(d, pts, AB, Elastic1D(E=ES), element="Truss2D",
                              bond=Bond(linear_law(1e11), DB))
    for nd in mesh.nodes_on("xmin"):
        nd.fix_dof("ux"); nd.fix_dof("uy")
    d.generate_equation_numbers()
    t1 = np.subtract(pts[1], pts[0]); t1 /= np.linalg.norm(t1)
    t2 = np.subtract(pts[2], pts[1]); t2 /= np.linalg.norm(t2)
    F = np.zeros(d.total_dofs)
    end = d.nodes[bar.node_ids[-1]]
    F[end.dofs["ux"]], F[end.dofs["uy"]] = 2.0e4 * t2
    r = run(d, F_applied=F)
    k_node = next(n for n in bar.node_ids
                  if np.allclose(d.nodes[n].coordinates[:2], pts[1], atol=1e-12))
    e1 = next(e for e in bar.element_ids if d.elements[e].nodes[1].id == k_node)
    e2 = next(e for e in bar.element_ids if d.elements[e].nodes[0].id == k_node)
    N1 = r.element_forces[e1].components["N"][1]
    N2 = r.element_forces[e2].components["N"][0]
    f_iface = np.zeros(2)
    for iid in bar.interface_ids:
        el = d.elements[iid]
        ids = [n.id for n in el.nodes]
        if k_node in ids:
            _, Fi = el.compute_element_state(el.get_local_displacements(r.U))
            j = ids.index(k_node)
            f_iface += Fi[2 * j:2 * j + 2]
    np.testing.assert_allclose(-f_iface, N1 * t1 - N2 * t2, rtol=0, atol=1e-8 * abs(N1))
    assert np.linalg.norm(N1 * t1 - N2 * t2) > 0.1 * abs(N1)    # el quiebre sí carga


def test_extraccion_3d():
    k = 1e10
    d, bar, r, _ = pullout(k, "Truss3D3", 10, dim=3, host="Hex8")
    x, s, _, pts = profile(d, bar, r)
    err = np.max(np.abs(s - exact_slip(x, k))) / exact_slip(0.0, k)
    assert err < 1e-4
    w = [wg for i in bar.interface_ids for wg in d.elements[i].weights]
    assert sum(p.bond_force * wg for p, wg in zip(pts, w)) == pytest.approx(-P, rel=1e-10)


def test_masa_y_modal():
    """La interfaz no tiene masa; el modal corre con la barra y el sólido."""
    d, bar, _, _ = pullout(1e11, "Truss2D3", 5, fixed_host=False)
    for i in bar.interface_ids:
        assert np.all(d.elements[i].compute_mass_matrix() == 0.0)
        assert np.all(d.elements[i].compute_mass_matrix("lumped") == 0.0)
    res = run_modal(d, n_modes=2)
    assert np.all(np.isfinite(res.frequencies_rad)) and res.frequencies_rad[0] > 0


def test_rechazos():
    from solidum import CohesiveDamageIsotropic
    with pytest.raises(ValueError, match="línea"):
        Bond(CohesiveDamageIsotropic(1e6, 100.0, 1e12, "linear"), DB)   # ley de superficie
    with pytest.raises(ValueError, match="integration"):
        Bond(linear_law(1e10), DB, integration="simpson")
    with pytest.raises(ValueError, match="diameter"):
        Bond(linear_law(1e10), 0.0)
    mesh = _plate(4)
    bar = embed_reinforcement(mesh.domain, [[0.0, 0.1], [1.0, 0.1]], AB, Elastic1D(E=ES),
                              bond=Bond(linear_law(1e10), DB))
    with pytest.raises(ValueError, match="interfaces"):
        bar.bond_forces(None)


def _yaml_text(mesh, constraint):
    d = mesh.domain
    lines = ["nodes:"]
    for nid, n in sorted(d.nodes.items()):
        lines.append(f"  - {{id: {nid}, coords: {[float(c) for c in n.coordinates]}}}")
    lines += ["materials:",
              f"  - {{id: 1, type: Elastic2D, E: {EC}, nu: 0.2, hypothesis: plane_stress}}",
              f"  - {{id: 2, type: Elastic1D, E: {ES}}}",
              "interface_materials:",
              "  - id: 1",
              "    type: BondSlipModelCode",
              f"    model_code_2010: {{fcm: 30.0e6, clear_rib_spacing: 0.01, diameter: {DB}, "
              f"E_c: {EC}}}",
              "elements:"]
    for eid, el in sorted(d.elements.items()):
        ids = [n.id for n in el.nodes]
        lines.append(f"  - {{id: {eid}, type: Tri3, nodes: {ids}, material: 1, thickness: {T}}}")
    lines += ["constraints:", f"  - {constraint}",
              "boundary_conditions_by_coord:", "  - {loc: x_min, ux: 0.0, uy: 0.0}",
              "point_loads_by_coord:", "  - {loc: x_max, uy: -2000.0}",
              "solver: {type: NonlinearSolver, num_steps: 4}"]
    return "\n".join(lines) + "\n"


def test_yaml_igual_que_python(tmp_path):
    """El bloque bond del YAML (con la ley del Model Code 2010 por su función
    auxiliar) da lo mismo que la llamada en Python."""
    pts = [[0.0, 0.05], [LB, 0.08]]
    mesh = _plate(8)
    text = _yaml_text(mesh, f"{{type: EmbeddedReinforcement, points: {pts}, area: {AB}, "
                            f"material: 2, element: Truss2D3, bond: {{law: 1, diameter: {DB}}}}}")
    path = tmp_path / "bond.yaml"
    path.write_text(text, encoding="utf-8")
    model = load_yaml(path)
    r_yaml = model.solve()
    (bar_yaml,) = model.parser.applied_constraints
    assert len(bar_yaml.interface_ids) == len(bar_yaml.element_ids) > 0

    from solidum import NonlinearSolver
    mesh_py = _plate(8)
    d = mesh_py.domain
    law = BondSlipModelCode.model_code_2010(30.0e6, clear_rib_spacing=0.01, diameter=DB, E_c=EC)
    bar_py = embed_reinforcement(d, pts, AB, Elastic1D(E=ES), element="Truss2D3",
                                 bond=Bond(law, DB))
    for nd in mesh_py.nodes_on("xmin"):
        nd.fix_dof("ux", 0.0); nd.fix_dof("uy", 0.0)
    d.generate_equation_numbers()
    F = np.zeros(d.total_dofs)
    for nd in mesh_py.nodes_on("xmax"):
        F[nd.dofs["uy"]] -= 2000.0
    asm = Assembler(d)
    r_py = run(d, assembler=asm, solver=NonlinearSolver(asm, num_steps=4), F_applied=F)
    assert bar_yaml.node_ids == bar_py.node_ids
    np.testing.assert_allclose(r_yaml.U, r_py.U, rtol=0, atol=1e-10 * np.abs(r_py.U).max())
