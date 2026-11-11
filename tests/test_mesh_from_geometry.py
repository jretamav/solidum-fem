"""Mallado desde la geometría en el YAML (ADR 0035).

``mesh: {geometry: placa.geo, parameters: {h: ...}, order: ..., element_shape: ...}``
llama a Gmsh, escribe el ``.msh`` y lo lee. Las pruebas que mallan necesitan
Gmsh (grupo opcional ``mesh``); las de validación del bloque, no.
"""
import textwrap
from pathlib import Path

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.mesh.gmsh_geometry import validate_geometry_block
from solidum.utils.yaml_parser import YamlParser, YamlValidationError

GEO = (Path(__file__).resolve().parent / "data" / "gmsh" / "placa.geo").as_posix()

_MATERIAL = """\
    materials:
      - {id: 1, type: Elastic2D, E: 200.0, nu: 0.3, hypothesis: plane_stress}
    mesh_thickness: 0.1
    """


def _parse(tmp_path, mesh_block: str):
    ruta = tmp_path / "modelo.yaml"
    ruta.write_text(textwrap.dedent(mesh_block) + textwrap.dedent(_MATERIAL), encoding="utf-8")
    return YamlParser(str(ruta)).parse()


def _tipos(dom):
    return {type(e).__name__ for e in dom.elements.values()}


def _area(dom):
    total = 0.0
    for el in dom.elements.values():
        xy = np.array([n.coordinates for n in el.nodes])[:3 if type(el).__name__.startswith("Tri") else 4]
        x, y = xy[:, 0], xy[:, 1]
        total += 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    return total


# ---------------------------------------------------------------------------
# Validación del bloque (sin Gmsh)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bloque,mensaje", [
    ({"geometry": "a.geo", "orden": 2}, "claves no admitidas"),
    ({"geometry": "a.msh"}, "archivo .geo"),
    ({"geometry": "a.geo", "order": 3}, "se admite 1 o 2"),
    ({"geometry": "a.geo", "element_shape": "hex"}, "'tri' o 'quad'"),
    ({"geometry": "a.geo", "parameters": {"h": "fino"}}, "se esperaba un número"),
    ({"geometry": "a.geo", "serendipity": "si"}, "true o false"),
    ({"geometry": "a.geo", "output": "a.vtk"}, "archivo .msh"),
])
def test_validacion_del_bloque(bloque, mensaje):
    errores = validate_geometry_block(bloque)
    assert any(mensaje in e for e in errores), errores


def test_el_yaml_rechaza_un_bloque_mal_formado(tmp_path):
    with pytest.raises(YamlValidationError, match="se admite 1 o 2"):
        _parse(tmp_path, f"""\
            mesh:
              geometry: {GEO}
              order: 3
            """)


# ---------------------------------------------------------------------------
# Mallado con Gmsh
# ---------------------------------------------------------------------------

@pytest.fixture
def gmsh():
    return pytest.importorskip("gmsh")


def test_ruta_al_geo_con_los_valores_por_omision(tmp_path, gmsh):
    dom = _parse(tmp_path, f"mesh: {GEO}\n")
    assert _tipos(dom) == {"Tri3"}
    assert _area(dom) == pytest.approx(2.0, rel=1e-12)


def test_el_tamano_del_yaml_sobrescribe_el_del_geo(tmp_path, gmsh):
    grueso = _parse(tmp_path, f"""\
        mesh:
          geometry: {GEO}
          output: grueso.msh
        """)
    fino = _parse(tmp_path, f"""\
        mesh:
          geometry: {GEO}
          parameters: {{h: 0.125}}
          output: fino.msh
        """)
    assert len(fino.elements) > 8 * len(grueso.elements) / 2
    assert (tmp_path / "grueso.msh").is_file() and (tmp_path / "fino.msh").is_file()


@pytest.mark.parametrize("forma,orden,serendipito,principal,permitidos", [
    ("tri", 2, True, "Tri6", {"Tri6"}),
    ("quad", 1, True, "Quad4", {"Quad4", "Tri3"}),
    ("quad", 2, True, "Quad8", {"Quad8", "Tri6"}),
    ("quad", 2, False, "Quad9", {"Quad9", "Tri6"}),
])
def test_orden_y_forma(tmp_path, gmsh, forma, orden, serendipito, principal, permitidos):
    """La recombinación en cuadriláteros puede dejar algún triángulo."""
    dom = _parse(tmp_path, f"""\
        mesh:
          geometry: {GEO}
          parameters: {{h: 0.25}}
          order: {orden}
          element_shape: {forma}
          serendipity: {str(serendipito).lower()}
          output: m.msh
        """)
    tipos = _tipos(dom)
    assert principal in tipos and tipos <= permitidos
    assert _area(dom) == pytest.approx(2.0, rel=1e-12)


def test_grupos_fisicos_del_geo_llegan_al_modelo(tmp_path, gmsh):
    dom = _parse(tmp_path, f"""\
        mesh:
          geometry: {GEO}
          parameters: {{h: 0.25}}
          order: 2
          element_shape: quad
          output: m.msh
        boundary_conditions_by_group:
          - {{group_name: Izquierda, ux: 0.0, uy: 0.0}}
        """)
    assert {"Placa", "Izquierda", "Derecha"} <= set(dom.physical_groups)
    izq = [dom.nodes[i] for i in dom.physical_groups["Izquierda"]]
    assert izq and all(abs(n.coordinates[0]) < 1e-12 for n in izq)


def test_parametro_inexistente_en_el_geo_es_un_error(tmp_path, gmsh):
    with pytest.raises(ValueError, match=r"\['hh'\] no aparecen en placa.geo"):
        _parse(tmp_path, f"""\
            mesh:
              geometry: {GEO}
              parameters: {{hh: 0.1}}
            """)
