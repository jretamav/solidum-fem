"""Calidad de la malla (ADR 0037), contra valores calculados a mano.

Trapecio Quad4 de vértices (0,0), (2,0), (1,1), (0,1): det J en los vértices
es un cuarto del producto vectorial de las dos aristas que salen de cada uno,
0.5, 0.5, 0.25 y 0.25, así que la razón de jacobianos es 2; aristas 2, √2, 1
y 1 (aspecto 2); ángulos 90°, 45°, 135° y 90°.
"""
import logging
from pathlib import Path

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.core.domain import Domain
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.mesh import mesh_quality, structured_mesh
from solidum.mesh.quality import LIMITS, log_mesh_quality
from solidum.registry import ElementRegistry, ThermalMaterialRegistry

MAT2 = Elastic2D(1.0, 0.2)
MAT3 = Elastic3D(1.0, 0.2)
GMSH = Path(__file__).parent / "data" / "gmsh"


def _uno(name, coords, **kw):
    dom = Domain()
    for i, c in enumerate(coords, 1):
        dom.add_node(i, [float(v) for v in c])
    kw.setdefault("material", MAT3 if len(coords[0]) == 3 else MAT2)
    dom.add_element(ElementRegistry.create(name, element_id=7,
                                           nodes=[dom.nodes[i] for i in range(1, len(coords) + 1)],
                                           **kw))
    return dom


TRAPECIO = [(0, 0), (2, 0), (1, 1), (0, 1)]


def test_trapecio_quad4():
    q = mesh_quality(_uno("Quad4", TRAPECIO))
    assert q.jacobian_ratio[0] == pytest.approx(2.0)
    assert q.aspect_ratio[0] == pytest.approx(2.0)
    assert q.min_angle[0] == pytest.approx(45.0)
    assert q.max_angle[0] == pytest.approx(135.0)
    assert q.n_inverted == 0 and q.warnings() == []


def test_quad8_con_nodos_medios_centrados_y_desplazados():
    vert = [(0, 0), (2, 0), (3, 1), (1, 1)]                      # paralelogramo
    medios = [(1, 0), (2.5, 0.5), (2, 1), (0.5, 0.5)]
    assert mesh_quality(_uno("Quad8", vert + medios)).jacobian_ratio[0] == pytest.approx(1.0)
    curvo = [(1, 0), (2.5, 0.5), (2, 1.3), (0.5, 0.5)]           # borde superior curvo
    assert mesh_quality(_uno("Quad8", vert + curvo)).jacobian_ratio[0] > 1.5


def test_triangulo_rectangulo_y_tetraedro():
    q = mesh_quality(_uno("Tri3", [(0, 0), (1, 0), (0, 1)]))
    assert q.jacobian_ratio[0] == pytest.approx(1.0)
    assert q.aspect_ratio[0] == pytest.approx(np.sqrt(2.0))
    assert (q.min_angle[0], q.max_angle[0]) == pytest.approx((45.0, 90.0))
    q = mesh_quality(_uno("Tet4", [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]))
    assert q.jacobian_ratio[0] == pytest.approx(1.0)
    assert np.isnan(q.min_angle[0])


def test_cubo_hex8_y_hex20():
    for name, n in (("Hex8", 1), ("Hex20", 1), ("Hex27", 2)):
        m = structured_mesh(name, (n, n, n), material=MAT3)
        q = mesh_quality(m.domain)
        np.testing.assert_allclose(q.jacobian_ratio, 1.0)
        np.testing.assert_allclose(q.aspect_ratio, 1.0)


def test_elemento_invertido():
    q = mesh_quality(_uno("Quad4", TRAPECIO[::-1]))
    assert q.n_inverted == 1 and np.isinf(q.jacobian_ratio[0])
    avisos = q.warnings()
    assert any("invertido" in a and "7" in a for a in avisos)


def test_avisos_por_aspecto_y_angulo():
    largo = mesh_quality(_uno("Quad4", [(0, 0), (30, 0), (30, 1), (0, 1)]))
    assert largo.aspect_ratio[0] == pytest.approx(30.0)
    assert any("razón de aspecto > 20" in a for a in largo.warnings())
    obtuso = mesh_quality(_uno("Tri3", [(0, 0), (2, 0), (1, 0.1)]))
    assert obtuso.max_angle[0] > LIMITS["max_angle_tri"]
    assert any("ángulo" in a for a in obtuso.warnings())
    # Convexo, con ~171° en el vértice (1, 0.6).
    flecha = mesh_quality(_uno("Quad4", [(0, 0), (2, 0), (1, 0.6), (0, 1)]))
    assert flecha.max_angle[0] > LIMITS["max_angle_quad"]
    assert any("ángulo" in a for a in flecha.warnings())


def test_razon_de_jacobianos_grande():
    """Convexo casi triangular: det J = 1/4 en (0,0) y 0.02/4 en (0.51, 0.51)."""
    q = mesh_quality(_uno("Quad4", [(0, 0), (1, 0), (0.51, 0.51), (0, 1)]))
    assert q.jacobian_ratio[0] == pytest.approx(50.0)
    assert any("razón de jacobianos" in a for a in q.warnings())


def test_elemento_termico():
    mat = ThermalMaterialRegistry.create("ThermalConduction", k=1.0)
    q = mesh_quality(_uno("Quad4Thermal", TRAPECIO, material=mat))
    assert q.jacobian_ratio[0] == pytest.approx(2.0)


def test_lineas_y_cables_no_se_miden():
    from solidum.materials.elastic import Elastic1D
    dom = _uno("Truss2D", [(0, 0), (1, 0)], A=1.0, material=Elastic1D(1.0))
    assert mesh_quality(dom).element_ids.size == 0


def test_resumen_en_el_registro_al_leer_un_msh(caplog):
    from solidum.utils.gmsh_parser import GmshParser
    dom = GmshParser(str(GMSH / "quad8.msh")).parse({"element": None, "kwargs": {"material": MAT2},
                                                     "optional": {}}, {})
    with caplog.at_level(logging.INFO, logger="solidum"):
        q = log_mesh_quality(dom)
    assert q.element_ids.size == len(dom.elements)
    assert "Calidad de la malla" in caplog.text
    assert "razón de jacobianos" in caplog.text


def test_el_yaml_registra_la_calidad(tmp_path, caplog):
    ruta = tmp_path / "m.yaml"
    ruta.write_text("""\
materials:
  - {id: 1, type: Elastic2D, E: 1.0, nu: 0.2}
mesh:
  structured: {element: Quad4, x: {length: 1, n: 2}, y: {length: 1, n: 2}}
boundary_conditions_by_group:
  xmin: {ux: 0.0, uy: 0.0}
solver: {type: LinearSolver}
""", encoding="utf-8")
    with caplog.at_level(logging.INFO, logger="solidum"):
        solidum.load_yaml(ruta)
    assert "Calidad de la malla (4 elementos" in caplog.text
