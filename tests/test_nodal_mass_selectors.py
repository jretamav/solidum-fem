"""Masas de nudo por grupo físico y por coordenada en el YAML (ADR 0024).

``nodal_masses_by_group`` y ``nodal_masses_by_coord`` seleccionan los nodos
como ``point_loads_by_group`` y ``point_loads_by_coord``: **cada** nodo
seleccionado recibe la masa indicada. Criterios:

1. Sobre una malla gmsh, la masa de grupo aumenta la diagonal de ``M`` en
   exactamente ``m`` en el grado de libertad nombrado de cada nodo del grupo,
   y no cambia nada más, con masa consistente y concentrada.
2. El selector por coordenada elige los nodos de esa coordenada; la forma dict
   de los dos bloques es igual a la forma lista; varias entradas sobre el
   mismo nodo se suman.
3. Los errores de forma y los del modelo construido se informan todos juntos
   y nombran el bloque y la entrada.
"""
import logging
import textwrap
from pathlib import Path

import numpy as np
import pytest

from solidum.math.assembly import Assembler
from solidum.utils.yaml_parser import YamlParser, YamlValidationError

KIRSCH = (Path(__file__).resolve().parents[1] / "examples" / "placa_agujero_kirsch"
          / "placa_agujero.msh")
M_NODE = 37.5       # kg
J_NODE = 2.25       # kg·m²

_BASE = f"""\
mesh: {KIRSCH.as_posix()}
mesh_material: 1
mesh_thickness: 0.01
materials:
  - {{id: 1, type: Elastic2D, E: 200.0e9, nu: 0.3, hypothesis: plane_stress, density: 7850.0}}
boundary_conditions_by_group:
  - {{group_name: "Simetria_X", ux: 0.0}}
  - {{group_name: "Simetria_Y", uy: 0.0}}
"""


@pytest.fixture(autouse=True)
def _silencio():
    logging.disable(logging.CRITICAL)
    yield
    logging.disable(logging.NOTSET)


def _model(tmp_path, extra=""):
    path = tmp_path / "modelo.yaml"
    path.write_text(_BASE + textwrap.dedent(extra), encoding="utf-8")
    parser = YamlParser(str(path))
    dom = parser.parse()
    dom.generate_equation_numbers()
    return dom


def _mass(dom, lumping="consistent"):
    return Assembler(dom).assemble_mass_matrix(lumping=lumping).toarray()


def _errors(tmp_path, extra):
    with pytest.raises(YamlValidationError) as exc:
        _model(tmp_path, extra)
    return exc.value.errors


class TestMasaPorGrupo:

    def test_cada_nodo_del_grupo_recibe_la_masa(self, tmp_path):
        dom0 = _model(tmp_path)
        dom = _model(tmp_path, f"""\
            nodal_masses_by_group:
              - {{group_name: Borde_Carga, ux: {M_NODE}}}
            """)
        group = dom.physical_groups["Borde_Carga"]
        assert len(group) > 2
        for lumping in ("consistent", "lumped"):
            expected = _mass(dom0, lumping)
            for nid in group:
                i = dom.nodes[nid].dofs["ux"]
                expected[i, i] += M_NODE
            np.testing.assert_array_equal(_mass(dom, lumping), expected)
        assert dom.nodal_masses == {(nid, "ux"): M_NODE for nid in group}


class TestMasaPorCoordenada:

    def test_selecciona_los_nodos_de_la_coordenada(self, tmp_path):
        dom = _model(tmp_path, f"""\
            nodal_masses_by_coord:
              - {{loc: x_max, uy: {M_NODE}}}
            """)
        x_max = max(n.coordinates[0] for n in dom.nodes.values() if n.dofs)
        expected = {n.id for n in dom.nodes.values()
                    if n.dofs and abs(n.coordinates[0] - x_max) < 1e-6}
        assert expected
        assert set(dom.nodal_masses) == {(nid, "uy") for nid in expected}

    def test_forma_dict_igual_que_lista_y_suma(self, tmp_path):
        lista = _model(tmp_path, f"""\
            nodal_masses_by_group:
              - {{group_name: Borde_Carga, ux: {M_NODE}, uy: {M_NODE}}}
            nodal_masses_by_coord:
              - {{coord: x, val: 0.0, uy: {J_NODE}}}
            """)
        dicc = _model(tmp_path, f"""\
            nodal_masses_by_group:
              Borde_Carga: {{ux: {M_NODE}, uy: {M_NODE}}}
            nodal_masses_by_coord:
              x_min: {{uy: {J_NODE}}}
            """)
        np.testing.assert_array_equal(_mass(lista), _mass(dicc))
        # Un nodo seleccionado dos veces acumula las dos masas.
        nid = sorted(lista.physical_groups["Borde_Carga"])[0]
        doble = _model(tmp_path, f"""\
            nodal_masses_by_group:
              - {{group_name: Borde_Carga, ux: {M_NODE}}}
            nodal_masses:
              - {{node_id: {nid}, ux: {M_NODE}}}
            """)
        assert doble.nodal_masses[(nid, "ux")] == 2 * M_NODE


class TestErrores:

    def test_forma(self, tmp_path):
        errors = _errors(tmp_path, """\
            nodal_masses_by_group:
              - {ux: 1.0}
              - {group_name: Borde_Carga, T: 1.0}
              - {group_name: Borde_Carga}
            nodal_masses_by_coord:
              - {ux: -1.0}
              - {loc: x_max, uy: .nan}
            """)
        text = "\n".join(errors)
        assert "nodal_masses_by_group[0]: falta el campo obligatorio 'group_name'" in text
        assert "nodal_masses_by_group[1]: 'T' no es un grado de libertad mecánico" in text
        assert "nodal_masses_by_group[2]: no declara ninguna masa" in text
        assert "nodal_masses_by_coord[0]: necesita 'loc' o el par 'coord'/'val'" in text
        assert "nodal_masses_by_coord[0].ux: debe ser un número finito no negativo" in text
        assert "nodal_masses_by_coord[1].uy: debe ser un número finito no negativo" in text
        assert len(errors) == 6, errors

    def test_modelo_construido(self, tmp_path):
        errors = _errors(tmp_path, f"""\
            nodal_masses_by_group:
              - {{group_name: NoExiste, ux: 1.0}}
              - {{group_name: Borde_Carga, rz: {J_NODE}}}
            nodal_masses_by_coord:
              - {{coord: x, val: 7.5, ux: 1.0}}
            """)
        text = "\n".join(errors)
        assert "nodal_masses_by_group[0]: grupo físico 'NoExiste' inexistente" in text
        assert "nodal_masses_by_group[1]: " in text
        assert "nodos no tienen el grado de libertad 'rz'" in text
        assert "nodal_masses_by_coord[0]: ningún nodo casa" in text
        assert len(errors) == 3, errors

    def test_grupo_sin_malla(self, tmp_path):
        path = tmp_path / "sin_malla.yaml"
        path.write_text(textwrap.dedent("""\
            nodes:
              - {id: 1, coords: [0.0, 0.0]}
              - {id: 2, coords: [1.0, 0.0]}
            materials:
              - {id: 1, type: Elastic1D, E: 2.0e11}
            elements:
              - {id: 1, type: Truss2D, nodes: [1, 2], material: 1, A: 1.0e-3}
            nodal_masses_by_group:
              piso: {ux: 1.0}
            """), encoding="utf-8")
        with pytest.raises(YamlValidationError, match="no tiene grupos físicos"):
            YamlParser(str(path)).parse()
