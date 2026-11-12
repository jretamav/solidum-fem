"""Mallas estructuradas por bloques (ADR 0037).

- ``axis_lines``: tramos con gradación geométrica (razón último/primero).
- ``structured_mesh``: líneas explícitas, celdas excluidas sin nodos
  huérfanos, cuatro particiones en triángulos y fusión de bloques; en todos
  los casos la malla es conforme (cada borde interior lo comparten dos
  elementos) y cubre el dominio (suma de áreas o volúmenes exacta).
- ``mesh: {structured: ...}`` en el YAML: prueba de la parcela (tracción
  uniforme) en mallas de dos bloques con gradación, para Quad4, Quad8, Tri3,
  Tri6, Hex8 y Tet4, con grupos de lado y por condición; errores de forma.
"""
import textwrap

import numpy as np
import pytest

import solidum
from solidum.core.domain import Domain
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.mesh import axis_lines, structured_mesh
from solidum.utils.expressions import ExpressionError, compile_predicate
from solidum.utils.yaml_parser import YamlValidationError

MAT2 = Elastic2D(200e3, 0.3)
MAT3 = Elastic3D(200e3, 0.3)


def _edges(domain):
    """Bordes (2D) o caras (3D) por conjunto de nodos, con cuántos elementos
    los comparten."""
    count = {}
    for el in domain.elements.values():
        attr = "FACE_NODES" if len(el.nodes[0].coordinates) == 3 else "EDGE_NODES"
        for local in getattr(el, attr):
            key = frozenset(el.nodes[j].id for j in local)
            count[key] = count.get(key, 0) + 1
    return count


def _conforme(domain) -> None:
    count = _edges(domain)
    assert max(count.values()) == 2, "un borde con más de dos elementos"
    usados = {n.id for el in domain.elements.values() for n in el.nodes}
    assert usados == set(domain.nodes), "nodos huérfanos"


def _area_tri_quad(domain) -> float:
    """Suma de áreas por la fórmula del polígono de los vértices (2D)."""
    total = 0.0
    for el in domain.elements.values():
        k = 3 if el.REFERENCE_CELL == "tri" else 4
        P = np.array([n.coordinates for n in el.nodes[:k]])
        x, y = P[:, 0], P[:, 1]
        a = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
        assert a > 0.0, f"elemento {el.id} con orientación negativa"
        total += a
    return total


# ----------------------------------------------------------------------------
# axis_lines
# ----------------------------------------------------------------------------

def test_axis_lines_uniforme_y_gradacion():
    np.testing.assert_allclose(axis_lines([(2.0, 4)], origin=1.0), [1.0, 1.5, 2.0, 2.5, 3.0])
    ln = axis_lines([(10.0, 5, 0.1)])
    h = np.diff(ln)
    assert ln[0] == 0.0 and ln[-1] == 10.0
    assert h[-1] / h[0] == pytest.approx(0.1)
    np.testing.assert_allclose(h[1:] / h[:-1], 0.1 ** 0.25)          # progresión geométrica
    ln = axis_lines([(950.0, 19, 0.1), (100.0, 40), (950.0, 19, 10.0)])
    assert ln.size == 79 and ln[-1] == pytest.approx(2000.0)
    assert np.diff(ln)[18] == pytest.approx(2.5, rel=1e-2) or np.diff(ln)[19] == pytest.approx(2.5)


@pytest.mark.parametrize("tramo", [(0.0, 2), (1.0, 0), (1.0, 2.5), (1.0, 2, -1.0), (1.0,)])
def test_axis_lines_rechaza(tramo):
    with pytest.raises(ValueError):
        axis_lines([tramo])


# ----------------------------------------------------------------------------
# structured_mesh: líneas, exclusión, diagonales, fusión
# ----------------------------------------------------------------------------

def test_lineas_explicitas_cubren_la_caja():
    lines = (axis_lines([(1.0, 3, 0.2), (2.0, 2)]), axis_lines([(1.5, 4, 3.0)], origin=-1.0))
    m = structured_mesh("Quad8", lines=lines, material=MAT2)
    assert m.divisions == (5, 4)
    assert _area_tri_quad(m.domain) == pytest.approx(3.0 * 1.5)
    xs = {round(n.coordinates[0], 12) for n in m.domain.nodes.values()}
    assert set(np.round(lines[0], 12)) <= xs
    _conforme(m.domain)


def test_lineas_y_size_a_la_vez_es_un_error():
    with pytest.raises(ValueError, match="lines"):
        structured_mesh("Quad4", lines=([0, 1], [0, 1]), size=(1, 1), material=MAT2)


def test_exclusion_de_una_ranura_sin_nodos_huerfanos():
    lines = (axis_lines([(1.0, 2), (0.5, 2), (1.0, 2)]), axis_lines([(1.0, 2), (1.0, 2)]))
    m = structured_mesh("Quad4", lines=lines, material=MAT2,
                        exclude=lambda x: 1.0 < x[0] < 1.5 and x[1] < 1.0)
    assert m.excluded == 4
    assert len(m.domain.elements) == 6 * 4 - 4
    assert _area_tri_quad(m.domain) == pytest.approx(2.5 * 2.0 - 0.5 * 1.0)
    _conforme(m.domain)
    # Los nodos del fondo de la ranura (x entre 1 y 1.5, y = 0) no existen.
    fondo = [n for n in m.domain.nodes.values()
             if 1.0 < n.coordinates[0] < 1.5 and n.coordinates[1] < 1.0]
    assert fondo == []


def test_excluir_todo_es_un_error():
    with pytest.raises(ValueError, match="todas las celdas"):
        structured_mesh("Quad4", (2, 2), material=MAT2, exclude=lambda x: True)


@pytest.mark.parametrize("element", ["Tri3", "Tri6"])
@pytest.mark.parametrize("diagonal,por_celda,vecinos_centro", [
    ("right", 2, 6), ("left", 2, 6), ("alternate", 2, 8), ("crossed", 4, 8)])
def test_diagonales(element, diagonal, por_celda, vecinos_centro):
    m = structured_mesh(element, (2, 2), size=(2.0, 2.0), material=MAT2, diagonal=diagonal)
    assert len(m.domain.elements) == 4 * por_celda
    assert _area_tri_quad(m.domain) == pytest.approx(4.0)
    _conforme(m.domain)
    centro = next(n for n in m.domain.nodes.values()
                  if np.allclose(n.coordinates, [1.0, 1.0]))
    tocan = [el for el in m.domain.elements.values() if centro in el.nodes]
    assert len(tocan) == vecinos_centro
    if diagonal == "crossed":
        # Tri3: el centro de cada celda; Tri6: el centro ya es nodo de la
        # rejilla y se añaden los 4 nodos intermedios de las semidiagonales.
        extra = 4 if element == "Tri3" else 4 * 4
        base = 9 if element == "Tri3" else 25
        assert len(m.domain.nodes) == base + extra


def test_diagonal_en_cuadrilateros_es_un_error():
    with pytest.raises(ValueError, match="triángulos"):
        structured_mesh("Quad4", (2, 2), material=MAT2, diagonal="left")
    with pytest.raises(ValueError, match="diagonal"):
        structured_mesh("Tri3", (2, 2), material=MAT2, diagonal="cruzada")


@pytest.mark.parametrize("element", ["Quad4", "Quad9", "Tri6"])
def test_fusion_de_dos_bloques(element):
    dom = Domain()
    a = structured_mesh(element, (2, 2), size=(1.0, 1.0), domain=dom, material=MAT2)
    b = structured_mesh(element, (3, 2), size=(1.5, 1.0), origin=(1.0, 0.0), domain=dom,
                        material=MAT2, merge_tol=1e-9)
    por_lado = 3 if element == "Quad4" else 5
    assert b.merged == por_lado
    assert _area_tri_quad(dom) == pytest.approx(2.5)
    _conforme(dom)
    compartidos = {n.id for n in a.nodes_on("xmax")}
    assert compartidos == {n.id for n in b.nodes_on("xmin")}


def test_sin_merge_tol_no_se_funden():
    dom = Domain()
    structured_mesh("Quad4", (1, 1), domain=dom, material=MAT2)
    b = structured_mesh("Quad4", (1, 1), origin=(1.0, 0.0), domain=dom, material=MAT2)
    assert b.merged == 0 and len(dom.nodes) == 8


def test_fusion_3d_hex_y_tet_conforme():
    for element in ("Hex8", "Tet4", "Tet10"):
        dom = Domain()
        structured_mesh(element, (2, 1, 1), size=(1.0, 1.0, 1.0), domain=dom, material=MAT3)
        structured_mesh(element, (1, 1, 1), size=(1.0, 1.0, 1.0), origin=(1.0, 0.0, 0.0),
                        domain=dom, material=MAT3, merge_tol=1e-9)
        _conforme(dom)


# ----------------------------------------------------------------------------
# Condiciones del compilador de expresiones
# ----------------------------------------------------------------------------

def test_condiciones():
    f = compile_predicate("x > 1 and not (y >= 2 or z < -1)")
    assert f(2, 0, 0) and not f(0, 0, 0) and not f(2, 3, 0) and not f(2, 0, -2)
    g = compile_predicate("near(x, 1)", functions={"near": lambda a, b: abs(a - b) < 1e-9})
    assert g(1.0, 0, 0) and not g(1.1, 0, 0)
    assert compile_predicate("0 < x < 1")(0.5, 0, 0)
    with pytest.raises(ExpressionError, match="near"):
        compile_predicate("x == 1")
    with pytest.raises(ExpressionError):
        compile_predicate("x > 1 if y else 0")


def test_las_cargas_siguen_sin_admitir_comparaciones():
    from solidum.utils.expressions import compile_expression
    with pytest.raises(ExpressionError):
        compile_expression("x > 1")


# ----------------------------------------------------------------------------
# YAML: prueba de la parcela en mallas de dos bloques con gradación
# ----------------------------------------------------------------------------

E, NU, SIGMA = 1000.0, 0.25, 3.0


def _parcela_2d(tmp_path, element, extra=""):
    ruta = tmp_path / "parcela.yaml"
    ruta.write_text(textwrap.dedent(f"""\
        materials:
          - {{id: 1, type: Elastic2D, E: {E}, nu: {NU}, hypothesis: plane_stress}}
        mesh:
          structured:
            blocks:
              - name: izq
                element: {element}
                x: [{{length: 1.0, n: 3, ratio: 0.3}}]
                y: [{{length: 0.6, n: 2}}, {{length: 0.4, n: 2, ratio: 2.0}}]
                {extra}
              - name: der
                element: {element}
                origin: [1.0, 0.0]
                x: {{length: 2.0, n: 4, ratio: 5.0}}
                y: [{{length: 0.6, n: 2}}, {{length: 0.4, n: 2, ratio: 2.0}}]
                {extra}
            groups:
              Origen: "near(x, 0) and near(y, 0)"
        boundary_conditions_by_group:
          izq.xmin: {{ux: 0.0}}
          Origen: {{uy: 0.0}}
        tractions_by_group:
          der.xmax: {{t: [{SIGMA}, 0.0]}}
        solver: {{type: LinearSolver}}
        """), encoding="utf-8")
    return ruta


def _resolver(ruta):
    model = solidum.load_yaml(ruta)
    return model, model.solve()


def _error_parcela(model, res, dim):
    """Error nodal relativo frente a la tracción uniforme exacta."""
    err = 0.0
    for n in model.domain.nodes.values():
        x = np.asarray(n.coordinates, dtype=float)
        exact = [SIGMA * x[0] / E, -NU * SIGMA * x[1] / E] + ([-NU * SIGMA * x[2] / E]
                                                             if dim == 3 else [])
        u = [res.U[n.dofs[d]] for d in ("ux", "uy", "uz")[:dim]]
        err = max(err, float(np.max(np.abs(np.subtract(u, exact)))))
    return err / (SIGMA * 3.0 / E)


@pytest.mark.parametrize("element,extra", [
    ("Quad4", ""), ("Quad8", ""), ("Quad9", ""),
    ("Tri3", "diagonal: alternate"), ("Tri6", "diagonal: crossed")])
def test_parcela_2d_por_bloques(tmp_path, element, extra):
    model, res = _resolver(_parcela_2d(tmp_path, element, extra))
    assert _error_parcela(model, res, 2) < 1e-10
    groups = model.domain.physical_groups
    assert {"izq.xmin", "izq.xmax", "der.xmin", "der.xmax", "Origen"} <= set(groups)
    assert groups["izq.xmax"] == groups["der.xmin"]           # lado común fundido


@pytest.mark.parametrize("element", ["Hex8", "Tet4"])
def test_parcela_3d_por_bloques(tmp_path, element):
    ruta = tmp_path / "parcela3d.yaml"
    ruta.write_text(textwrap.dedent(f"""\
        materials:
          - {{id: 1, type: Elastic3D, E: {E}, nu: {NU}}}
        mesh_element: {element}
        mesh:
          structured:
            blocks:
              - {{x: {{length: 1.0, n: 2, ratio: 0.5}}, y: {{length: 1.0, n: 2}},
                  z: {{length: 0.5, n: 1}}}}
              - {{origin: [1.0, 0.0, 0.0], x: {{length: 1.0, n: 2}}, y: {{length: 1.0, n: 2}},
                  z: {{length: 0.5, n: 1}}}}
            groups:
              Origen: "near(x, 0) and near(y, 0) and near(z, 0)"
              EjeY: "near(x, 0) and near(z, 0)"
        boundary_conditions_by_group:
          block1.xmin: {{ux: 0.0}}
          Origen: {{uy: 0.0, uz: 0.0}}
          EjeY: {{uz: 0.0}}
        tractions_by_group:
          block2.xmax: {{t: [{SIGMA}, 0.0, 0.0]}}
        solver: {{type: LinearSolver}}
        """), encoding="utf-8")
    model, res = _resolver(ruta)
    assert _error_parcela(model, res, 3) < 1e-10


def test_un_bloque_sin_nombre_y_su_entalla(tmp_path):
    ruta = tmp_path / "viga.yaml"
    ruta.write_text(textwrap.dedent("""\
        materials:
          - {id: 1, type: Elastic2D, E: 30000.0, nu: 0.2, hypothesis: plane_stress}
        mesh_thickness: 50.0
        mesh:
          structured:
            element: Quad4
            x: [{length: 95, n: 5, ratio: 0.2}, {length: 10, n: 10}, {length: 95, n: 5, ratio: 5}]
            y: [{length: 10, n: 5}, {length: 10, n: 2}]
            exclude: "x > 99 and x < 101 and y < 10"
            groups:
              Carga: "near(x, 100) and near(y, 20)"
        boundary_conditions_by_coord:
          - {coord: x, val: 0.0, ux: 0.0, uy: 0.0}
          - {coord: x, val: 200.0, uy: 0.0}
        point_loads_by_group:
          Carga: {uy: -10.0}
        solver: {type: LinearSolver}
        """), encoding="utf-8")
    model, res = _resolver(ruta)
    dom = model.domain
    assert set(dom.physical_groups) == {"xmin", "xmax", "ymin", "ymax", "Carga"}
    for el in dom.elements.values():
        c = np.mean([n.coordinates for n in el.nodes], axis=0)
        assert not (99 < c[0] < 101 and c[1] < 10)
    assert res.U[dom.nodes[dom.physical_groups["Carga"][0]].dofs["uy"]] < 0.0


@pytest.mark.parametrize("bloque,mensaje", [
    ({"structured": {"element": "Quad4", "y": {"length": 1, "n": 2}}}, "falta 'x'"),
    ({"structured": {"x": {"length": 0, "n": 2}, "y": {"length": 1, "n": 2}}}, "'length'"),
    ({"structured": {"x": {"length": 1, "n": 2.5}, "y": {"length": 1, "n": 2}}}, "'n'"),
    ({"structured": {"x": [{"length": 1, "n": 2, "razon": 2}], "y": {"length": 1, "n": 2}}},
     "no admitidas"),
    ({"structured": {"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2},
                     "diagonal": "cruzada"}}, "diagonal"),
    ({"structured": {"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2},
                     "exclude": "x == 1"}}, "near"),
    ({"structured": {"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2},
                     "groups": {"A": "w > 1"}}}, "'w'"),
    ({"structured": {"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2},
                     "merge_tol": -1}}, "merge_tol"),
    ({"structured": {"blocks": [{"x": {"length": 1, "n": 1}}], "element": "Quad4"}},
     "dentro de cada bloque"),
    ({"structured": {"x": {"length": 1, "n": 2}}, "geometry": "a.geo"}, "no se admiten"),
])
def test_validacion_de_forma(tmp_path, bloque, mensaje):
    import yaml
    ruta = tmp_path / "m.yaml"
    ruta.write_text(yaml.safe_dump({
        "materials": [{"id": 1, "type": "Elastic2D", "E": 1.0, "nu": 0.2}],
        "mesh_element": "Quad4", "mesh": bloque, "solver": {"type": "LinearSolver"}}),
        encoding="utf-8")
    with pytest.raises(YamlValidationError, match=mensaje):
        solidum.run_yaml(ruta)


@pytest.mark.parametrize("bloque,mensaje", [
    ({"x": {"length": 1, "n": 2}}, "2D: se esperan los ejes"),
    ({"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2}, "z": {"length": 1, "n": 1}},
     "sobran"),
    ({"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2}, "diagonal": "left"},
     "triángulos"),
    ({"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2}, "groups": {"A": "x > 5"}},
     "ningún nodo"),
    ({"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2}, "groups": {"xmin": "x < 5"}},
     "lado de un bloque"),
    ({"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2}, "element": "Cable2DCorot"},
     "no se puede mallar por bloques"),
    ({"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2}, "element": "Truss2D"},
     "1D: se esperan los ejes"),
])
def test_errores_al_mallar(tmp_path, bloque, mensaje):
    import yaml
    ruta = tmp_path / "m.yaml"
    ruta.write_text(yaml.safe_dump({
        "materials": [{"id": 1, "type": "Elastic2D", "E": 1.0, "nu": 0.2}],
        "mesh_element": "Quad4", "mesh": {"structured": bloque},
        "solver": {"type": "LinearSolver"}}), encoding="utf-8")
    with pytest.raises(YamlValidationError, match=mensaje):
        solidum.run_yaml(ruta)


def test_falta_el_elemento(tmp_path):
    import yaml
    ruta = tmp_path / "m.yaml"
    ruta.write_text(yaml.safe_dump({
        "materials": [{"id": 1, "type": "Elastic2D", "E": 1.0, "nu": 0.2}],
        "mesh": {"structured": {"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2}}},
        "solver": {"type": "LinearSolver"}}), encoding="utf-8")
    with pytest.raises(YamlValidationError, match="falta 'element'"):
        solidum.run_yaml(ruta)


def test_mesh_physical_groups_no_se_combina(tmp_path):
    import yaml
    ruta = tmp_path / "m.yaml"
    ruta.write_text(yaml.safe_dump({
        "materials": [{"id": 1, "type": "Elastic2D", "E": 1.0, "nu": 0.2}],
        "mesh_element": "Quad4",
        "mesh": {"structured": {"x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 2}}},
        "mesh_physical_groups": {"A": {"material": 1}},
        "solver": {"type": "LinearSolver"}}), encoding="utf-8")
    with pytest.raises(YamlValidationError, match="mesh_physical_groups"):
        solidum.run_yaml(ruta)


def test_mapping_desde_el_yaml(tmp_path):
    """Cuarto de corona 1 ≤ r ≤ 2 con Quad9: los nodos quedan sobre los arcos."""
    import yaml
    ruta = tmp_path / "m.yaml"
    ruta.write_text(yaml.safe_dump({
        "materials": [{"id": 1, "type": "Elastic2D", "E": 1.0, "nu": 0.2}],
        "mesh": {"structured": {
            "element": "Quad9", "x": {"length": 1, "n": 2}, "y": {"length": 1, "n": 4},
            "mapping": ["(1 + x)*cos(pi/2*y)", "(1 + x)*sin(pi/2*y)"]}},
        "boundary_conditions_by_group": {"ymin": {"uy": 0.0}, "ymax": {"ux": 0.0}},
        "solver": {"type": "LinearSolver"}}), encoding="utf-8")
    model = solidum.load_yaml(ruta)
    r = [np.hypot(*nd.coordinates) for nd in model.domain.nodes.values()]
    assert min(r) == pytest.approx(1.0) and max(r) == pytest.approx(2.0)
    interior = [nd for nd in model.domain.nodes.values()
                if nd.id in model.domain.physical_groups["xmax"]]
    np.testing.assert_allclose([np.hypot(*nd.coordinates) for nd in interior], 2.0)
