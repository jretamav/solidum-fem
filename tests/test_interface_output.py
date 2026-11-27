"""Salida de los elementos de interfaz (ADR 0043 §7): VTK por celda, archivo
por puntos de integración e historia.

Invariantes que protegen estas pruebas:
- Cada interfaz con celda tiene ``Jump_<c>`` y ``Traction_<c>`` en ejes
  locales (media de sus puntos); ``NaN`` en las celdas de sólido, y los
  esfuerzos de sólido ``NaN`` en las de interfaz. En el pre-proceso el salto
  normal es la holgura inicial.
- Las interfaces no entran en el suavizado nodal: el esfuerzo nodal de los
  sólidos de la cara es el de sus celdas (ni la mitad ni ``NaN``).
- Un código de la ley (``STATE_CODES``: el régimen del contacto) no se
  promedia por celda; sí está en el archivo por puntos y en la historia.
- La historia lee el salto (de ``U``) y la tracción (del estado confirmado)
  por componente local y selecciona los elementos de una interfaz por su
  nombre; la adherencia barra–sólido gana la misma historia.
- La energía por zona admite ``interface: <nombre>``.
Valores con parámetros ≠ 1 (Reglas §6): presión de tres resortes en serie
``σ = δ/(2/E + 1/K_n)`` con espesor 3 (*calculado*).
"""
from __future__ import annotations

import csv
import textwrap

import numpy as np
import pytest

meshio = pytest.importorskip("meshio")

import solidum  # noqa: E402
from solidum.core.interface_material import component_names  # noqa: E402
from solidum.materials.elastic_3d import Elastic3D  # noqa: E402
from solidum.materials.interface import ContactPenalty  # noqa: E402
from solidum.mesh import structured_mesh  # noqa: E402
from solidum.mesh.interfaces import split_interfaces  # noqa: E402
from solidum.utils.vtk_exporter import VtkExporter  # noqa: E402
from solidum.utils.yaml_parser import YamlValidationError  # noqa: E402

E_B, K_N, T_B, GAP, BAJADA = 1000.0, 1.0e5, 3.0, 0.1, 0.2
SIGMA = (BAJADA - GAP) / (2.0 / E_B + 1.0 / K_N)      # presión final (calculado)

BLOQUES = f"""\
materials:
  - {{id: 1, type: Elastic2D, E: {E_B}, nu: 0.0}}
interface_materials:
  - {{id: 1, type: ContactPenalty, K_n: {K_N}}}
mesh_element: Quad4
mesh_thickness: {T_B}
mesh:
  structured:
    blocks:
      - {{name: abajo, material: 1, x: {{length: 2, n: 4}}, y: {{length: 1, n: 2}}}}
      - {{name: arriba, material: 1, origin: [0, {1 + GAP}], x: {{length: 2, n: 4}},
         y: {{length: 1, n: 2}}}}
    groups:
      Origen: "near(x, 0) and near(y, 0)"
interfaces:
  - {{name: junta, pair: [abajo.ymax, arriba.ymin], material: 1}}
boundary_conditions_by_group:
  abajo.ymin: {{uy: 0.0}}
  Origen: {{ux: 0.0}}
  arriba.ymax: {{uy: -{BAJADA}, ux: 0.0}}
solver: {{type: NonlinearSolver, num_steps: 4}}
"""


def _write(tmp_path, texto, name="m.yaml"):
    ruta = tmp_path / name
    ruta.write_text(textwrap.dedent(texto), encoding="utf-8")
    return ruta


def _history(ruta):
    with open(ruta.with_name(ruta.stem + "_historia.csv"), encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _cells(mesh, name):
    return np.concatenate(mesh.cell_data[name])


@pytest.fixture
def corrida(tmp_path):
    ruta = _write(tmp_path, BLOQUES + """\
output:
  pre_process: {export: true}
  results: {frequency: all_steps, integration_points: true}
  history:
    columns:
      - {name: gn, element: {var: jump, component: n, interface: junta}}
      - {name: gs, element: {var: jump, component: s, interface: junta, reduce: max}}
      - {name: tn, element: {var: traction, component: n, interface: junta, reduce: min}}
      - {name: p, element: {var: contact_pressure, interface: junta, reduce: max}}
      - {name: regimen, element: {var: status, interface: junta, reduce: max}}
      - {name: E_junta, energy: elastic, elements: {interface: junta}}
""")
    model = solidum.load_yaml(ruta)
    model.solve(write_output=True)
    return ruta, model


# ---------------------------------------------------------------------------
# VTK por celda
# ---------------------------------------------------------------------------

def test_celdas_de_interfaz_con_salto_y_traccion_locales(corrida):
    ruta, model = corrida
    m = meshio.read(ruta.with_name("m_final.vtu"))
    jn, tn = _cells(m, "Jump_n"), _cells(m, "Traction_n")
    is_if = ~np.isnan(jn)
    assert is_if.sum() == len(model.domain.interfaces["junta"]) == 4
    np.testing.assert_allclose(tn[is_if], -SIGMA, rtol=1e-9)
    np.testing.assert_allclose(jn[is_if], -SIGMA / K_N, rtol=1e-9)
    np.testing.assert_allclose(_cells(m, "Jump_s")[is_if], 0.0, atol=1e-12)
    np.testing.assert_allclose(_cells(m, "contact_pressure")[is_if], SIGMA, rtol=1e-9)
    # «No aplica» en las celdas del otro tipo.
    assert np.all(np.isnan(_cells(m, "Traction_n")[~is_if]))
    for name in ("Sigma_YY", "Von_Mises"):
        v = _cells(m, name)
        assert np.all(np.isnan(v[is_if])) and np.all(np.isfinite(v[~is_if]))
    np.testing.assert_allclose(_cells(m, "Sigma_YY")[~is_if], -SIGMA, rtol=1e-9)


def test_el_regimen_no_se_promedia_por_celda(corrida):
    ruta, _ = corrida
    m = meshio.read(ruta.with_name("m_final.vtu"))
    assert "status" not in m.cell_data and "contact_pressure" in m.cell_data


def test_las_interfaces_no_entran_en_el_suavizado_nodal(corrida):
    """Los nodos de la cara tocan dos sólidos y dos interfaces: con las
    interfaces en el promedio, σ_yy nodal sería NaN (o la mitad con ceros)."""
    ruta, model = corrida
    m = meshio.read(ruta.with_name("m_final.vtu"))
    np.testing.assert_allclose(m.point_data["Sigma_YY_nodal"], -SIGMA, rtol=1e-9)


def test_el_pre_proceso_muestra_la_holgura_inicial(corrida):
    ruta, _ = corrida
    m = meshio.read(ruta.with_name("m_pre_proceso.vtu"))
    jn = _cells(m, "Jump_n")
    np.testing.assert_allclose(jn[~np.isnan(jn)], GAP, rtol=1e-12)
    np.testing.assert_allclose(_cells(m, "Traction_n")[~np.isnan(jn)], 0.0)


# ---------------------------------------------------------------------------
# Archivo por puntos de integración
# ---------------------------------------------------------------------------

def test_puntos_de_las_interfaces(corrida):
    ruta, model = corrida
    m = meshio.read(ruta.with_name("m_puntos_final.vtu"))
    ids = model.domain.interfaces["junta"]
    on = np.isin(m.point_data["Element_Id"], ids)
    assert on.sum() == 2 * len(ids)                           # regla nodal: 2 por cara
    np.testing.assert_allclose(m.point_data["Traction_n"][on], -SIGMA, rtol=1e-9)
    np.testing.assert_allclose(m.point_data["status"][on], 1.0)
    np.testing.assert_allclose(m.point_data["contact_pressure"][on], SIGMA, rtol=1e-9)
    assert np.all(np.isnan(m.point_data["Sigma_YY"][on]))
    assert np.all(np.isnan(m.point_data["Jump_n"][~on]))
    # Los puntos están sobre la cara A (y = 1) y se mueven con ella: la
    # compresión del bloque de abajo, σ·1/E.
    np.testing.assert_allclose(m.points[on, 1], 1.0)
    np.testing.assert_allclose(m.point_data["Displacements"][on, 1], -SIGMA / E_B, rtol=1e-9)


# ---------------------------------------------------------------------------
# Historia
# ---------------------------------------------------------------------------

def test_historia_de_la_interfaz_paso_a_paso(corrida):
    """Bajada impuesta 0.05 por paso: abierta (holgura 0.05 y 0), cerrada con
    0.05 y 0.1 de cierre."""
    ruta, _ = corrida
    filas = _history(ruta)
    assert len(filas) == 4
    for fila in filas:
        delta = BAJADA * float(fila["factor"])
        sigma = max(delta - GAP, 0.0) / (2.0 / E_B + 1.0 / K_N)
        gn = GAP - delta if delta <= GAP else -sigma / K_N
        assert float(fila["gn"]) == pytest.approx(gn, abs=1e-12 + 1e-9 * abs(gn))
        assert float(fila["gs"]) == pytest.approx(0.0, abs=1e-12)
        assert float(fila["tn"]) == pytest.approx(-sigma, rel=1e-9, abs=1e-12)
        assert float(fila["p"]) == pytest.approx(sigma, rel=1e-9, abs=1e-12)
        assert float(fila["regimen"]) == (0.0 if delta <= GAP else 1.0)
        # ψ = ½·K_n·g_n² sobre el área de contacto 2 × 3.
        energia = 0.5 * K_N * min(gn, 0.0) ** 2 * 2.0 * T_B
        assert float(fila["E_junta"]) == pytest.approx(energia, rel=1e-9, abs=1e-12)


@pytest.mark.parametrize("columna,mensaje", [
    ("{var: jump, component: n, interface: grieta}", "interfaz 'grieta' inexistente"),
    ("{var: jump, component: s1, interface: junta}", "componente 's1'"),
    ("{var: traction, component: n, where: 'y < 1'}", "no tienen la tracción"),
    ("{var: jump, interface: junta}", "componente local"),
    ("{var: stress, component: xx, interface: junta}", "no tienen esfuerzos"),
])
def test_errores_de_las_columnas_de_interfaz(tmp_path, columna, mensaje):
    ruta = _write(tmp_path, BLOQUES + f"""\
output:
  history:
    columns:
      - {{name: c, element: {columna}}}
""")
    with pytest.raises((ValueError, YamlValidationError), match=mensaje):
        solidum.load_yaml(ruta).solve(write_output=True)


def test_element_results_admite_el_salto_y_rechaza_el_regimen(tmp_path):
    ruta = _write(tmp_path, BLOQUES + """\
output:
  results: {element_results: [Jump_n, Traction_n, contact_pressure]}
""")
    solidum.load_yaml(ruta).solve(write_output=True)
    m = meshio.read(ruta.with_name("m_final.vtu"))
    assert set(m.cell_data) == {"Jump_n", "Traction_n", "contact_pressure"}
    ruta = _write(tmp_path, BLOQUES + """\
output:
  results: {element_results: [status]}
""", name="n.yaml")
    with pytest.raises(ValueError, match="ángulos o códigos"):
        solidum.load_yaml(ruta).solve(write_output=True)


def test_la_adherencia_gana_la_misma_historia(tmp_path):
    """La adherencia barra–sólido (``BondSlipInterface``, sin celda VTK) da
    en la historia el salto y la tracción de su registro por punto."""
    ruta = _write(tmp_path, """\
materials:
  - {id: 1, type: Elastic2D, E: 3.0e10, nu: 0.2}
  - {id: 2, type: Elastic1D, E: 2.0e11}
interface_materials:
  - {id: 1, type: BondSlipMultilinear, points: [[0.001, 5.0e6]], k_n: 1.0e12}
mesh_element: Quad4
mesh_thickness: 0.1
mesh:
  structured:
    x: {length: 1.0, n: 4}
    y: {length: 0.2, n: 2}
constraints:
  - {type: EmbeddedReinforcement, points: [[0.0, 0.1], [1.0, 0.1]], area: 2.0e-4,
     material: 2, element: Truss2D, bond: {law: 1, diameter: 0.016}}
boundary_conditions_by_group:
  xmin: {ux: 0.0, uy: 0.0}
point_loads_by_coord:
  - {loc: x_max, ux: 1.0e4}
solver: {type: NonlinearSolver, num_steps: 2}
output:
  history:
    columns:
      - {name: s, element: {var: jump, component: s, reduce: max}}
      - {name: tau, element: {var: traction, component: s, reduce: max}}
""")
    model = solidum.load_yaml(ruta)
    res = model.solve(write_output=True)
    (bar,) = model.parser.applied_constraints
    pts = [p for i in bar.interface_ids
           for p in model.domain.elements[i].interface_state(res.U)]
    fila = _history(ruta)[-1]
    assert max(p.slip for p in pts) != 0.0
    assert float(fila["s"]) == pytest.approx(max(p.slip for p in pts), rel=1e-9)
    assert float(fila["tau"]) == pytest.approx(max(p.bond_stress for p in pts), rel=1e-9)


# ---------------------------------------------------------------------------
# Celdas de los demás elementos de interfaz y nombres de las componentes
# ---------------------------------------------------------------------------

def _split_3d(element):
    mat = Elastic3D(E=1000.0, nu=0.25)
    dom = structured_mesh(element, lines=[np.linspace(0, 2, 3), np.linspace(0, 1, 2),
                                          np.linspace(0, 2, 3)], material=mat).domain
    ids = sorted(i for i, n in dom.nodes.items()
                 if abs(n.coordinates[0] - 1.0) < 1e-9 and n.coordinates[2] <= 1.0 + 1e-9)
    dom.physical_groups, dom.physical_group_cells = {"grieta": ids}, {"grieta": []}
    return dom


def _split_2d(element):
    from solidum.materials.elastic_2d import Elastic2D
    dom = structured_mesh(element, lines=[np.linspace(0, 4, 5), np.linspace(0, 2, 5)],
                          material=Elastic2D(E=1000.0, nu=0.25)).domain
    ids = sorted(i for i, n in dom.nodes.items()
                 if abs(n.coordinates[0] - 2.0) < 1e-9 and n.coordinates[1] <= 1.0 + 1e-9)
    dom.physical_groups, dom.physical_group_cells = {"grieta": ids}, {"grieta": []}
    return dom


@pytest.mark.parametrize("build,element,cell,comps", [
    (_split_2d, "Quad4", "quad", ("n", "s")),
    (_split_2d, "Quad8", "quad6", ("n", "s")),
    (_split_3d, "Hex8", "hexahedron", ("n", "s1", "s2")),
    (_split_3d, "Tet4", "wedge", ("n", "s1", "s2")),
])
def test_celda_propia_y_componentes_de_cada_interfaz(tmp_path, build, element, cell, comps):
    dom = build(element)
    made = split_interfaces(dom, "grieta", ContactPenalty(K_n=1.0e5), name="g")
    assert dom.interfaces["g"] == made.elements
    dom.generate_equation_numbers()
    U = np.zeros(dom.total_dofs)
    ruta = str(tmp_path / "g.vtu")
    VtkExporter(dom).export(ruta, U)
    m = meshio.read(ruta)
    blocks = {c.type: k for k, c in enumerate(m.cells)}
    assert cell in blocks
    jn = m.cell_data["Jump_n"][blocks[cell]]
    assert np.sum(~np.isnan(jn)) == len(made.elements)
    # Ida y vuelta por meshio con los nodos en el orden del elemento (el de la
    # celda VTK, sin permutar).
    conn = m.cells[blocks[cell]].data[~np.isnan(jn)]
    for eid, row in zip(made.elements, conn):
        X = np.array([np.pad(n.coordinates, (0, 3 - len(n.coordinates)))
                      for n in dom.elements[eid].nodes], float)
        np.testing.assert_allclose(m.points[row], X)
    names = {f"{p}_{c}" for p in ("Jump", "Traction") for c in comps}
    assert names <= set(m.cell_data)
    other = {"Jump_s", "Jump_s1", "Jump_s2"} - names
    assert not other & set(m.cell_data)


def test_sin_las_tablas_de_meshio_la_celda_quad6_da_un_error_claro(tmp_path, monkeypatch):
    """Si otra versión de meshio no deja completar sus tablas, el exportador
    lo dice en vez de fallar con un ``KeyError 'quad6'`` (revisión del
    arquitecto del 2026-11-26)."""
    import meshio._common as common
    dom = _split_2d("Quad8")
    split_interfaces(dom, "grieta", ContactPenalty(K_n=1.0e5), name="g")
    dom.generate_equation_numbers()
    tabla = dict(common.num_nodes_per_cell)
    tabla.pop("quad6", None)
    monkeypatch.setattr(common, "num_nodes_per_cell", tabla)
    with pytest.raises(RuntimeError, match="no sabe escribir las celdas"):
        VtkExporter(dom).export(str(tmp_path / "g.vtu"), np.zeros(dom.total_dofs))


K_E = 1.0e4      # rigidez elástica de la cohesiva (apertura sin daño)


def test_la_celda_lleva_la_media_de_sus_puntos(tmp_path):
    """Apertura lineal a lo largo de la grieta, ``g(y) = a + b·y``, en régimen
    elástico de una cohesiva: por punto, ``g`` y ``t = K_e·g`` en los nodos
    de la cara; por celda, su media (el valor en el centro del lado)."""
    from solidum.materials.interface import CohesiveDamageIsotropic
    dom = _split_2d("Quad4")
    law = CohesiveDamageIsotropic(sigma_t0=3.0, G_f=0.1, K_e=K_E, softening="linear")
    made = split_interfaces(dom, "grieta", law, name="g")
    dom.generate_equation_numbers()
    U = np.zeros(dom.total_dofs)
    a, b = 1.0e-5, 4.0e-5
    for eid in made.elements:
        el = dom.elements[eid]
        sign = float(el.R[0][0][0])               # n = ±x
        for j in el.B_NODES:
            node = el.nodes[j]
            U[node.dofs["ux"]] = sign * (a + b * float(node.coordinates[1]))
    VtkExporter(dom).export(str(tmp_path / "g.vtu"), U)
    m = meshio.read(str(tmp_path / "g.vtu"))
    jn = _cells(m, "Jump_n")
    on = ~np.isnan(jn)
    # g de cada par (B − A) en la normal; la punta no se duplica (g = 0 allí).
    esperado = []
    for eid in made.elements:
        el = dom.elements[eid]
        n = el.R[0][0]
        g = [n @ (U[[el.nodes[jb].dofs[d] for d in ("ux", "uy")]]
                  - U[[el.nodes[ja].dofs[d] for d in ("ux", "uy")]])
             for ja, jb in zip(el.A_NODES, el.B_NODES)]
        assert g[0] != g[1]                       # apertura no uniforme en la celda
        esperado.append(np.mean(g))
    np.testing.assert_allclose(jn[on], esperado, rtol=1e-12)
    np.testing.assert_allclose(_cells(m, "Traction_n")[on], K_E * np.array(esperado),
                               rtol=1e-12)


def test_historia_de_una_interfaz_cohesiva_en_un_analisis_lineal(tmp_path):
    """El análisis lineal (sin pasos) con una cohesiva en régimen elástico:
    la historia da el salto y la tracción del resorte en serie."""
    ruta = _write(tmp_path, f"""\
materials:
  - {{id: 1, type: Elastic2D, E: {E_B}, nu: 0.0}}
interface_materials:
  - {{id: 1, type: CohesiveDamageIsotropic, sigma_t0: 3.0, G_f: 0.1, K_e: {K_E},
     softening: linear}}
mesh_element: Quad4
mesh_thickness: {T_B}
mesh:
  structured:
    x: {{length: 4, n: 4}}
    y: {{length: 1, n: 1}}
    groups:
      corte: "near(x, 2)"
interfaces:
  - {{name: g, split: corte, material: 1}}
boundary_conditions_by_group:
  xmin: {{ux: 0.0}}
  xmax: {{ux: 1.0e-4}}
  ymin: {{uy: 0.0}}
solver: {{type: LinearSolver}}
output:
  history:
    columns:
      - {{name: gn, element: {{var: jump, component: n, interface: g}}}}
      - {{name: tn, element: {{var: traction, component: n, interface: g}}}}
""")
    solidum.load_yaml(ruta).solve(write_output=True)
    (fila,) = _history(ruta)
    # Barra de dos tramos de longitud 2 con un resorte K_e entre ellos:
    # σ = δ/(4/E + 1/K_e) (calculado).
    sigma = 1.0e-4 / (4.0 / E_B + 1.0 / K_E)
    assert float(fila["tn"]) == pytest.approx(sigma, rel=1e-9)
    assert float(fila["gn"]) == pytest.approx(sigma / K_E, rel=1e-9)


def test_nombres_de_las_componentes_locales():
    assert component_names("surface", 2) == ("n", "s")
    assert component_names("surface", 3) == ("n", "s1", "s2")
    assert component_names("line", 2) == ("n", "s")
    assert component_names("line", 3) == ("n1", "n2", "s")
    with pytest.raises(ValueError):
        component_names("surface", 1)


def test_dos_interfaces_con_el_mismo_nombre_son_error():
    dom = _split_2d("Quad4")
    split_interfaces(dom, "grieta", ContactPenalty(K_n=1.0e5), name="g")
    with pytest.raises(ValueError, match="ya hay una interfaz 'g'"):
        split_interfaces(dom, "grieta", ContactPenalty(K_n=1.0e5), name="g")
