"""Refuerzo embebido con adherencia perfecta (spec ``EmbeddedReinforcement``,
ADR 0032).

Criterios de la spec:
- deformación uniforme con la barra a un ángulo cualquiera: ε de cada tramo =
  tᵀεt y energía de la barra ½E_sA_s(tᵀεt)²L_b (Tri3, Quad4 distorsionado,
  Hex8, Tet4);
- modos rígidos: 3 en 2D y 6 en 3D;
- la malla de la barra no cambia el resultado con anfitriones lineales;
- placa en tracción: rigidez E_c·A_g + E_s·A_s (sin descontar el área);
- fluencia de la barra en A_s·σ_y;
- viga contra la sección transformada (sin descuento);
- fuerzas de adherencia en equilibrio con la barra;
- masa de la barra en el modelo;
- YAML idéntico a Python, en 2D y en 3D.

Más la localización de puntos (inversa del mapeo isoparamétrico) en los diez
sólidos, con mallas curvas.
"""
import logging
import textwrap

import numpy as np
import pytest

from solidum import (
    Assembler,
    Elastic1D,
    Elastic2D,
    Elastoplastic1D,
    NonlinearSolver,
    embed_reinforcement,
    load_yaml,
    run,
    run_modal,
)
from solidum.materials.elastic_3d import Elastic3D
from solidum.utils.yaml_parser import YamlValidationError
from solidum.mesh import structured_mesh
from solidum.mesh.point_location import PointLocator, natural_coordinates, shape_values

EC, NU, ES = 30e9, 0.2, 200e9
SY = 420e6


@pytest.fixture(autouse=True)
def _silencio():
    logging.disable(logging.WARNING)
    yield
    logging.disable(logging.NOTSET)


def _material(dim, nu=NU, E=EC, density=None):
    if dim == 2:
        return Elastic2D(E=E, nu=nu, hypothesis="plane_stress", density=density)
    return Elastic3D(E=E, nu=nu, density=density)


def _mesh(element, divisions, size, mapping=None, nu=NU, E=EC, density=None):
    dim = len(divisions)
    kwargs = {"material": _material(dim, nu, E, density)}
    if dim == 2:
        kwargs["thickness"] = 0.2
    return structured_mesh(element, divisions, size=size, mapping=mapping, **kwargs)


def _distort_2d(X):
    x, y = X
    return np.array([x + 0.04 * np.sin(2 * np.pi * y), y + 0.03 * np.sin(np.pi * x)])


def _distort_3d(X):
    x, y, z = X
    return np.array([x + 0.03 * np.sin(np.pi * y), y + 0.03 * np.sin(np.pi * z),
                     z + 0.03 * np.sin(np.pi * x)])


# ---------------------------------------------------------------------------
# Localización de puntos
# ---------------------------------------------------------------------------

_ALL_SOLIDS = [("Tri3", 2), ("Tri6", 2), ("Quad4", 2), ("Quad8", 2), ("Quad9", 2),
               ("Tet4", 3), ("Tet10", 3), ("Hex8", 3), ("Hex20", 3), ("Hex27", 3)]


def _random_xi(cell, rng):
    if cell in ("quad", "hex"):
        return rng.uniform(-0.9, 0.9, size=2 if cell == "quad" else 3)
    w = rng.dirichlet(np.ones(3 if cell == "tri" else 4))
    return w[1:]


@pytest.mark.parametrize("element,dim", _ALL_SOLIDS, ids=[e for e, _ in _ALL_SOLIDS])
def test_inversa_del_mapeo_en_mallas_curvas(element, dim):
    """x = Σ N(ξ)·X y vuelta: ξ recuperado a 1e-9 en los diez sólidos, con
    mallas distorsionadas (lados curvos en los de orden superior)."""
    rng = np.random.default_rng(3)
    divisions, mapping = ((2, 2), _distort_2d) if dim == 2 else ((2, 1, 1), _distort_3d)
    mesh = _mesh(element, divisions, size=None, mapping=mapping)
    for el in mesh.domain.elements.values():
        X = el.get_coordinate_matrix(ndim=dim)
        for _ in range(5):
            xi = _random_xi(el.REFERENCE_CELL, rng)
            x = shape_values(el, xi) @ X
            np.testing.assert_allclose(natural_coordinates(el, x), xi, atol=1e-9)


def test_localizador_asigna_el_borde_comun_al_menor_id():
    mesh = _mesh("Quad4", (2, 1), size=(2.0, 1.0))
    locator = PointLocator(mesh.domain.elements.values(), 2)
    index, _ = locator.locate([1.0, 0.5])         # borde común de los elementos 1 y 2
    assert locator.elements[index].id == min(mesh.domain.elements)
    assert locator.locate([2.5, 0.5]) is None


# ---------------------------------------------------------------------------
# Deformación uniforme con la barra a un ángulo cualquiera
# ---------------------------------------------------------------------------

_UNIFORM = [("Tri3", 2, None), ("Quad4", 2, _distort_2d), ("Hex8", 3, None), ("Tet4", 3, None)]


@pytest.mark.parametrize("element,dim,mapping", _UNIFORM, ids=[e for e, _, _ in _UNIFORM])
def test_deformacion_uniforme_angulo_arbitrario(element, dim, mapping):
    """Con u = ε·x impuesto en todos los nodos del sólido, cada tramo de la
    barra tiene ε_b = tᵀεt y la barra almacena ½·E_s·A_s·(tᵀεt)²·L_b."""
    divisions = (4, 3) if dim == 2 else (3, 2, 2)
    mesh = _mesh(element, divisions, size=None, mapping=mapping)
    dom = mesh.domain
    if dim == 2:
        p, q = np.array([0.13, 0.21]), np.array([0.87, 0.74])
        eps = np.array([[1.0e-4, 0.35e-4], [0.35e-4, -0.6e-4]])
    else:
        p, q = np.array([0.12, 0.18, 0.21]), np.array([0.83, 0.77, 0.69])
        eps = np.array([[1.0e-4, 0.3e-4, -0.2e-4], [0.3e-4, -0.5e-4, 0.25e-4],
                        [-0.2e-4, 0.25e-4, 0.7e-4]])
    AS = 3.0e-4
    bar = embed_reinforcement(dom, [p, q], AS, Elastic1D(E=ES))
    assert len(bar.element_ids) > 3                       # cruza varios anfitriones
    names = ("ux", "uy", "uz")[:dim]
    for nid, node in dom.nodes.items():
        if nid in dom.generated_nodes:
            continue
        u = eps @ np.asarray(node.coordinates[:dim], dtype=float)
        for k, d in enumerate(names):
            node.fix_dof(d, float(u[k]))
    r = run(dom)
    t = (q - p) / np.linalg.norm(q - p)
    eps_b = float(t @ eps @ t)
    energy = 0.0
    for eid in bar.element_ids:
        el = dom.elements[eid]
        a, b = el.nodes
        ua = np.array([r.U[a.dofs[d]] for d in names])
        ub = np.array([r.U[b.dofs[d]] for d in names])
        L = np.linalg.norm(np.subtract(b.coordinates[:dim], a.coordinates[:dim]))
        assert t @ (ub - ua) / L == pytest.approx(eps_b, rel=1e-10)
        N = r.element_forces[eid].components["N"][0]
        assert N == pytest.approx(ES * AS * eps_b, rel=1e-10)
        energy += 0.5 * N * eps_b * L
    assert energy == pytest.approx(0.5 * ES * AS * eps_b ** 2 * bar.length, rel=1e-10)


# ---------------------------------------------------------------------------
# Modos rígidos
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("element,dim,n_rigid", [("Tri3", 2, 3), ("Hex8", 3, 6)])
def test_modos_rigidos(element, dim, n_rigid):
    divisions = (4, 2) if dim == 2 else (2, 2, 2)
    mesh = _mesh(element, divisions, size=None)
    dom = mesh.domain
    p = [0.1, 0.3] if dim == 2 else [0.1, 0.3, 0.4]
    q = [0.9, 0.6] if dim == 2 else [0.9, 0.6, 0.55]
    embed_reinforcement(dom, [p, q], 3.0e-4, Elastic1D(E=ES))
    dom.generate_equation_numbers()
    asm = Assembler(dom)
    asm.assemble_system()
    K_red, _, _, _ = asm.reduce(asm.K_global, np.zeros(dom.total_dofs))
    w = np.linalg.eigvalsh(K_red.toarray())
    assert int(np.sum(np.abs(w) < 1e-10 * w.max())) == n_rigid


# ---------------------------------------------------------------------------
# La malla de la barra es independiente
# ---------------------------------------------------------------------------

def _cantilever_tri3(points, subdivisions=1):
    mesh = _mesh("Tri3", (8, 2), size=(2.0, 0.4))
    dom = mesh.domain
    host_ids = sorted(dom.nodes)
    bar = embed_reinforcement(dom, points, 4.0e-4, Elastic1D(E=ES), subdivisions=subdivisions)
    for n in mesh.nodes_on("xmin"):
        n.fix_dof("ux"); n.fix_dof("uy")
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    for n in mesh.nodes_on("xmax"):
        F[n.dofs["uy"]] -= 1.0e4
    r = run(dom, F_applied=F)
    U = np.array([[r.U[dom.nodes[i].dofs[d]] for d in ("ux", "uy")] for i in host_ids])
    return U, bar


def test_malla_de_barra_independiente():
    """Con anfitriones lineales, añadir vértices colineales o subdividir los
    tramos no cambia la solución del sólido."""
    U_ref, bar_ref = _cantilever_tri3([[0.0, 0.07], [2.0, 0.07]])
    U_pts, bar_pts = _cantilever_tri3([[0.0, 0.07], [0.31, 0.07], [0.77, 0.07],
                                       [1.42, 0.07], [2.0, 0.07]])
    U_sub, bar_sub = _cantilever_tri3([[0.0, 0.07], [2.0, 0.07]], subdivisions=3)
    assert len(bar_sub.node_ids) > len(bar_pts.node_ids) > len(bar_ref.node_ids)
    np.testing.assert_allclose(U_pts, U_ref, rtol=0, atol=1e-10 * np.abs(U_ref).max())
    np.testing.assert_allclose(U_sub, U_ref, rtol=0, atol=1e-10 * np.abs(U_ref).max())


def test_barra_fuera_del_solido():
    mesh = _mesh("Quad4", (2, 1), size=(2.0, 1.0))
    with pytest.raises(ValueError, match="sale de los anfitriones|no está dentro"):
        embed_reinforcement(mesh.domain, [[0.5, 0.5], [2.5, 0.5]], 1e-4, Elastic1D(E=ES))


# ---------------------------------------------------------------------------
# Placa en tracción y fluencia de la barra
# ---------------------------------------------------------------------------

LX, LY, TH = 1.0, 0.2, 0.2
AS_PLATE = 4.0e-4


def _plate(element, E_c=EC, bar_material=None):
    mesh = _mesh(element, (5, 2), size=(LX, LY), E=E_c)
    dom = mesh.domain
    embed_reinforcement(dom, [[0.0, 0.07], [LX, 0.07]], AS_PLATE,
                        bar_material or Elastic1D(E=ES))
    return mesh, dom


@pytest.mark.parametrize("element", ["Tri3", "Quad4"])
def test_placa_traccion_uniaxial(element):
    """Rigidez axial efectiva E_c·A_g + E_s·A_s: el área de la barra no se
    descuenta del concreto (decisión del usuario, 2026-11-10)."""
    mesh, dom = _plate(element)
    delta = 1.0e-4
    for n in mesh.nodes_on("xmin"):
        n.fix_dof("ux", 0.0)
    right = mesh.nodes_on("xmax")
    for n in right:
        n.fix_dof("ux", delta)
    mesh.nodes_on("xmin")[0].fix_dof("uy", 0.0)
    r = run(dom)
    F = sum(r.reactions_by_node[n.id]["ux"] for n in right)
    assert F / (delta / LX) == pytest.approx(EC * LY * TH + ES * AS_PLATE, rel=1e-10)


def test_fluencia_de_la_barra():
    """Con un concreto de rigidez despreciable la placa cede cuando la barra
    alcanza A_s·σ_y."""
    E_soft = 1.0e3
    mesh, dom = _plate("Quad4", E_c=E_soft, bar_material=Elastoplastic1D(E=ES, sigma_y=SY))
    delta = 5 * SY / ES * LX
    for n in mesh.nodes_on("xmin"):
        n.fix_dof("ux", 0.0)
    right = mesh.nodes_on("xmax")
    for n in right:
        n.fix_dof("ux", delta)
    mesh.nodes_on("xmin")[0].fix_dof("uy", 0.0)
    dom.generate_equation_numbers()
    asm = Assembler(dom)
    r = run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=10))
    F = sum(r.reactions_by_node[n.id]["ux"] for n in right)
    host = E_soft * LY * TH * delta / LX
    assert F - host == pytest.approx(AS_PLATE * SY, rel=1e-8)


# ---------------------------------------------------------------------------
# Viga contra la sección transformada
# ---------------------------------------------------------------------------

def _beam_deflection(element, nx, ny):
    """Viga simplemente apoyada (u_y = 0 en las dos secciones extremas),
    esfuerzo plano con ν = 0, carga uniforme en la cara superior, refuerzo
    inferior. Flecha media de la sección central."""
    L, H, T, q = 4.0, 0.4, 0.2, 2.0e4
    mesh = structured_mesh(element, (nx, ny), size=(L, H),
                           material=Elastic2D(E=EC, nu=0.0, hypothesis="plane_stress"),
                           thickness=T)
    dom = mesh.domain
    embed_reinforcement(dom, [[0.0, 0.05], [L, 0.05]], 1.6e-3, Elastic1D(E=ES))
    for n in mesh.nodes_on("xmin") + mesh.nodes_on("xmax"):
        n.fix_dof("uy")
    min(mesh.nodes_on("ymin"), key=lambda n: n.coordinates[0]).fix_dof("ux")
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    for el, edge in mesh.boundary["ymax"]:
        fe = el.compute_edge_traction(edge, np.array([0.0, -q]))
        for k, nd in enumerate(el.nodes):
            F[nd.dofs["ux"]] += fe[2 * k]
            F[nd.dofs["uy"]] += fe[2 * k + 1]
    r = run(dom, F_applied=F)
    mid = [n for n in dom.nodes.values()
           if abs(n.coordinates[0] - L / 2) < 1e-9 and n.id not in dom.generated_nodes]
    return float(np.mean([r.U[n.dofs["uy"]] for n in mid]))


def _transformed_section_deflection():
    """Euler-Bernoulli con la sección transformada sin descuento más el
    cortante de Timoshenko (κ = 5/6) de la sección de concreto."""
    L, H, T, q = 4.0, 0.4, 0.2, 2.0e4
    w, Ag, ys, As = q * T, H * T, 0.05, 1.6e-3
    ybar = (EC * Ag * H / 2 + ES * As * ys) / (EC * Ag + ES * As)
    EI = EC * (T * H ** 3 / 12 + Ag * (H / 2 - ybar) ** 2) + ES * As * (ys - ybar) ** 2
    return -(5 * w * L ** 4 / (384 * EI) + w * L ** 2 / (8 * (5 / 6) * (EC / 2) * Ag))


def test_viga_seccion_transformada():
    ref = _transformed_section_deflection()
    assert _beam_deflection("Quad8", 20, 4) == pytest.approx(ref, rel=5e-3)
    errors = [abs(_beam_deflection("Quad4", nx, ny) / ref - 1)
              for nx, ny in ((20, 2), (40, 4), (80, 8))]
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 1e-2


# ---------------------------------------------------------------------------
# Fuerzas de adherencia
# ---------------------------------------------------------------------------

def test_fuerzas_de_adherencia_equilibran_la_barra():
    """En cada nodo de la barra, la fuerza de adherencia es la fuerza interna
    de sus tramos; sobre toda la barra suman cero."""
    mesh = _mesh("Quad4", (6, 2), size=(2.0, 0.4))
    dom = mesh.domain
    bar = embed_reinforcement(dom, [[0.0, 0.05], [1.9, 0.33]], 4.0e-4, Elastic1D(E=ES))
    for n in mesh.nodes_on("xmin"):
        n.fix_dof("ux"); n.fix_dof("uy")
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    for n in mesh.nodes_on("xmax"):
        F[n.dofs["uy"]] -= 1.0e4
    r = run(dom, F_applied=F)
    bond = bar.bond_forces(r)
    f_int = np.zeros_like(bond)
    for k, eid in enumerate(bar.element_ids):
        el = dom.elements[eid]
        u_e = np.array([r.U[i] for i in el.get_global_dof_indices()])
        _, fe = el.compute_element_state(u_e)
        f_int[k] += fe[:2]
        f_int[k + 1] += fe[2:]
    scale = np.abs(f_int).max()
    np.testing.assert_allclose(bond, f_int, rtol=0, atol=1e-9 * scale)
    np.testing.assert_allclose(bond.sum(axis=0), 0.0, atol=1e-9 * scale)


# ---------------------------------------------------------------------------
# Masa
# ---------------------------------------------------------------------------

def test_masa_de_la_barra():
    rho_c, rho_s, As = 2400.0, 7850.0, 4.0e-4
    mesh = _mesh("Quad4", (8, 2), size=(2.0, 0.4), density=rho_c)
    dom = mesh.domain
    points = [[0.0, 0.05], [1.2, 0.05], [2.0, 0.3]]
    bar = embed_reinforcement(dom, points, As, Elastic1D(E=ES, density=rho_s))
    for n in mesh.nodes_on("xmin"):
        n.fix_dof("ux"); n.fix_dof("uy")
    res = run_modal(dom, n_modes=2)
    expected = rho_c * 2.0 * 0.4 * 0.2 + rho_s * As * bar.length
    assert res.total_mass["ux"] == pytest.approx(expected, rel=1e-12)
    assert res.total_mass["uy"] == pytest.approx(expected, rel=1e-12)


# ---------------------------------------------------------------------------
# YAML
# ---------------------------------------------------------------------------

def _yaml_from_mesh(mesh, element, dim, constraint, extra):
    dom = mesh.domain
    lines = ["nodes:"]
    for nid, n in sorted(dom.nodes.items()):
        lines.append(f"  - {{id: {nid}, coords: {[float(c) for c in n.coordinates]}}}")
    lines.append("materials:")
    if dim == 2:
        lines.append(f"  - {{id: 1, type: Elastic2D, E: {EC}, nu: {NU}, hypothesis: plane_stress}}")
    else:
        lines.append(f"  - {{id: 1, type: Elastic3D, E: {EC}, nu: {NU}}}")
    lines.append(f"  - {{id: 2, type: Elastic1D, E: {ES}}}")
    lines.append("elements:")
    thick = ", thickness: 0.2" if dim == 2 else ""
    for eid, el in sorted(dom.elements.items()):
        ids = [n.id for n in el.nodes]
        lines.append(f"  - {{id: {eid}, type: {element}, nodes: {ids}, material: 1{thick}}}")
    lines.append("constraints:")
    lines.append(f"  - {constraint}")
    return "\n".join(lines) + "\n" + textwrap.dedent(extra)


@pytest.mark.parametrize("dim", [2, 3])
def test_yaml_igual_que_python(tmp_path, dim):
    """El bloque ``constraints`` da lo mismo que la llamada en Python. El
    apoyo por coordenada (x_min) no toma el nodo de la barra en x = 0."""
    element = "Quad4" if dim == 2 else "Hex8"
    divisions, size = ((4, 2), (2.0, 0.4)) if dim == 2 else ((4, 1, 1), (2.0, 0.4, 0.4))
    p = [0.0, 0.05] if dim == 2 else [0.0, 0.05, 0.1]
    q = [2.0, 0.3] if dim == 2 else [2.0, 0.3, 0.25]
    AS = 4.0e-4
    fixed = "{loc: x_min, ux: 0.0, uy: 0.0}" if dim == 2 else "{loc: x_min, ux: 0.0, uy: 0.0, uz: 0.0}"
    loads = "{loc: x_max, uy: -1000.0}"
    extra = f"""\
boundary_conditions_by_coord:
  - {fixed}
point_loads_by_coord:
  - {loads}
solver: {{type: LinearSolver}}
"""
    mesh = _mesh(element, divisions, size)
    text = _yaml_from_mesh(mesh, element, dim,
                           f"{{type: EmbeddedReinforcement, points: {[p, q]}, area: {AS}, material: 2}}",
                           extra)
    path = tmp_path / "modelo.yaml"
    path.write_text(text, encoding="utf-8")
    model = load_yaml(path)
    r_yaml = model.solve()
    (bar_yaml,) = model.parser.applied_constraints

    mesh_py = _mesh(element, divisions, size)
    dom = mesh_py.domain
    bar_py = embed_reinforcement(dom, [p, q], AS, Elastic1D(E=ES))
    names = ("ux", "uy", "uz")[:dim]
    for n in mesh_py.nodes_on("xmin"):
        for d in names:
            n.fix_dof(d, 0.0)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    for n in mesh_py.nodes_on("xmax"):
        F[n.dofs["uy"]] -= 1000.0
    r_py = run(dom, F_applied=F)

    assert bar_yaml.node_ids == bar_py.node_ids
    assert bar_yaml.host_ids == bar_py.host_ids
    np.testing.assert_allclose(r_yaml.U, r_py.U, rtol=0, atol=1e-12 * np.abs(r_py.U).max())


def test_yaml_tipo_de_restriccion_desconocido(tmp_path):
    mesh = _mesh("Quad4", (2, 1), (2.0, 0.4))
    text = _yaml_from_mesh(mesh, "Quad4", 2, "{type: Inexistente, points: [[0, 0.1], [1, 0.1]]}",
                           "solver: {type: LinearSolver}\n")
    path = tmp_path / "modelo.yaml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(YamlValidationError, match="Inexistente"):
        load_yaml(path)


def test_yaml_parametro_no_aceptado(tmp_path):
    mesh = _mesh("Quad4", (2, 1), (2.0, 0.4))
    text = _yaml_from_mesh(
        mesh, "Quad4", 2,
        "{type: EmbeddedReinforcement, points: [[0.1, 0.1], [1.0, 0.1]], area: 1.0e-4, "
        "material: 2, diametro: 0.016}",
        "solver: {type: LinearSolver}\n")
    path = tmp_path / "modelo.yaml"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(YamlValidationError, match="diametro"):
        load_yaml(path)
