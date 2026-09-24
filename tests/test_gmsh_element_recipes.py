"""Receta de elemento en la ruta gmsh (ADR 0020, P7).

Antes la malla gmsh sólo producía Quad4 y Tri3 con materiales mecánicos: un
modelo térmico, o cualquier elemento que no fuera el de la celda, no podía
salir de gmsh. Ahora ``mesh_element`` (toda la malla) o ``element`` (un grupo
físico) eligen el tipo; sus referencias se buscan en la familia que declara
el elemento, y sus parámetros se validan contra su constructor.
"""
import textwrap
from pathlib import Path

import pytest

from solidum.utils.yaml_parser import YamlParser, YamlValidationError

EJEMPLOS = Path(__file__).resolve().parents[1] / "examples"
KIRSCH = EJEMPLOS / "placa_agujero_kirsch" / "placa_agujero.msh"   # cuadriláteros
PLACA = EJEMPLOS / "placa.msh"                                      # triángulos


def _parsear(tmp_path, malla: Path, *bloques: str):
    ruta = tmp_path / "modelo.yaml"
    texto = "".join(textwrap.dedent(b) for b in bloques)
    ruta.write_text(f"mesh: {malla.as_posix()}\n" + texto, encoding="utf-8")
    parser = YamlParser(str(ruta))
    return parser, parser.parse()


_TERMICO = """\
    thermal_materials:
      - {id: 1, type: ThermalConduction, k: 50.0}
    """
_MECANICO = """\
    materials:
      - {id: 1, type: Elastic2D, E: 2.0e11, nu: 0.3, hypothesis: plane_stress}
    """


@pytest.mark.skipif(not KIRSCH.is_file(), reason="malla de Kirsch no disponible")
class TestRecetaDeElemento:

    def test_modelo_termico_desde_gmsh(self, tmp_path):
        parser, dom = _parsear(tmp_path, KIRSCH, _TERMICO, """\
            mesh_element: Quad4Thermal
            mesh_thickness: 0.02
            """)
        tipos = {type(e).__name__ for e in dom.elements.values()}
        assert tipos == {"Quad4Thermal"}
        elem = next(iter(dom.elements.values()))
        assert elem.material is parser.thermal_materials[1]
        assert elem.thickness == 0.02

    def test_por_omision_sigue_siendo_el_elemento_de_la_celda(self, tmp_path):
        _, dom = _parsear(tmp_path, KIRSCH, _MECANICO, "mesh_quadrature: '1x1'\n")
        tipos = {type(e).__name__ for e in dom.elements.values()}
        assert tipos == {"Quad4"}
        assert next(iter(dom.elements.values())).quadrature_key == "1x1"

    def test_material_de_la_familia_equivocada(self, tmp_path):
        """El elemento térmico busca ``mesh_material`` en ``thermal_materials``."""
        with pytest.raises(YamlValidationError,
                           match=r"mesh: material=1 no existe en 'thermal_materials'"):
            _parsear(tmp_path, KIRSCH, _MECANICO, "mesh_element: Quad4Thermal\n")

    def test_parametro_no_aceptado_por_el_elemento_del_grupo(self, tmp_path):
        with pytest.raises(YamlValidationError,
                           match="parámetro 'espesor' no aceptado por 'Quad4Thermal'"):
            _parsear(tmp_path, KIRSCH, _TERMICO, """\
                mesh_element: Quad4Thermal
                mesh_physical_groups:
                  Placa: {espesor: 0.1}
                """)

    def test_parametros_sin_elemento_se_rechazan(self, tmp_path):
        """Antes un grupo con una errata (``thicknes``) se ignoraba en silencio."""
        with pytest.raises(YamlValidationError, match=r"parámetros \['thicknes'\] sin 'element'"):
            _parsear(tmp_path, KIRSCH, _MECANICO, """\
                mesh_physical_groups:
                  Placa: {thicknes: 0.1}
                """)

    def test_tipo_de_elemento_desconocido(self, tmp_path):
        with pytest.raises(YamlValidationError, match="tipo de elemento desconocido 'Quad5'"):
            _parsear(tmp_path, KIRSCH, _MECANICO, "mesh_element: Quad5\n")


@pytest.mark.skipif(not PLACA.is_file(), reason="placa.msh no disponible")
def test_elemento_incompatible_con_la_celda(tmp_path):
    """Un cuadrilátero térmico no puede ocupar una celda triangular."""
    with pytest.raises(ValueError, match="'Quad4Thermal' tiene 4 nodos y la celda gmsh 'triangle', 3"):
        _parsear(tmp_path, PLACA, _TERMICO, "mesh_element: Quad4Thermal\n")
