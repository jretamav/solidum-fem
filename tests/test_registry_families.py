"""Familias de material genéricas (ADR 0020, P1 y P2).

Un registro que declara ``YAML_SECTION`` es una familia de material: el
lector YAML construye y valida sus objetos a partir de esa sección sin
conocerla por nombre. Es lo que permite a un módulo de usuario declarar su
propia familia sin tocar el programa principal, y lo que ya usan las
familias mecánica y térmica.
"""
import textwrap

import pytest

from solidum.registry import (
    ElementRegistry,
    MaterialRegistry,
    Registry,
    ThermalMaterialRegistry,
)
from solidum.utils.yaml_parser import YamlParser, YamlValidationError


def _parsear(tmp_path, texto: str):
    ruta = tmp_path / "modelo.yaml"
    ruta.write_text(textwrap.dedent(texto), encoding="utf-8")
    parser = YamlParser(str(ruta))
    parser.parse()
    return parser


_BARRA = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
materials:
  - {id: 1, type: Elastic1D, E: 200.0e9}
elements:
  - {id: 1, type: Truss2D, nodes: [1, 2], material: 1, A: 1.0e-3}
"""

_PARED_TERMICA = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [1.0, 1.0]}
  - {id: 4, coords: [0.0, 1.0]}
thermal_materials:
  - {id: 1, type: ThermalConduction, k: 50.0}
elements:
  - {id: 1, type: Quad4Thermal, nodes: [1, 2, 3, 4], material: 1}
"""


class TestRegistro:

    def test_las_familias_del_programa_principal_declaran_su_seccion(self):
        secciones = {f.YAML_SECTION for f in Registry.families()}
        assert {"materials", "thermal_materials"} <= secciones
        assert MaterialRegistry in Registry.families()
        assert ThermalMaterialRegistry in Registry.families()

    def test_elementos_y_solvers_no_son_familias_de_material(self):
        assert ElementRegistry not in Registry.families()

    def test_cada_subclase_recibe_su_propio_almacen(self):
        """Antes, olvidar ``_items`` en una subclase compartía en silencio el
        dict de la base con todos los demás registros."""
        class Uno(Registry):
            pass

        class Dos(Registry):
            pass

        Uno.register("A", int)
        assert Uno.names() == ["A"]
        assert Dos.names() == []
        assert Registry._items == {}

    def test_dos_familias_no_pueden_declarar_la_misma_seccion(self):
        with pytest.raises(ValueError, match="thermal_materials"):
            class Duplicada(Registry):
                YAML_SECTION = "thermal_materials"
        assert all(f.__qualname__ != "TestRegistro.test_dos_familias_no_pueden_"
                   "declarar_la_misma_seccion.<locals>.Duplicada"
                   for f in Registry.families())

    def test_cada_registro_declara_su_tipo_de_spec(self):
        """``tools/spec.py`` localiza el registro por ``SPEC_KIND``, sin lista fija."""
        assert Registry.for_spec_kind("element") is ElementRegistry
        assert Registry.for_spec_kind("thermal_material") is ThermalMaterialRegistry
        assert {"element", "material", "thermal_material", "solver"} <= set(Registry.spec_kinds())
        with pytest.raises(ValueError, match="Tipo de spec 'element' ya declarado"):
            class Otra(Registry):
                SPEC_KIND = "element"

    def test_una_familia_nueva_la_recorre_el_lector_sin_tocarlo(self, tmp_path):
        """Un registro declarado fuera del lector YAML (como haría un módulo de
        usuario) aporta su sección y sus objetos."""
        class Resorte:
            def __init__(self, k: float):
                self.k = k

        class FamiliaDePrueba(Registry):
            _kind = "Resorte"
            YAML_SECTION = "springs_de_prueba"
            YAML_LABEL = "resorte"

        FamiliaDePrueba.register("Resorte", Resorte)
        try:
            parser = _parsear(tmp_path, _BARRA + "springs_de_prueba:\n  - {id: 7, type: Resorte, k: 3.0}\n")
            assert parser.springs_de_prueba[7].k == 3.0
            with pytest.raises(YamlValidationError, match="parámetro 'c' no aceptado por 'Resorte'"):
                _parsear(tmp_path, _BARRA + "springs_de_prueba:\n  - {id: 7, type: Resorte, k: 3.0, c: 1.0}\n")
        finally:
            Registry._families.remove(FamiliaDePrueba)


class TestLectorYaml:

    def test_seccion_desconocida_es_un_error(self, tmp_path):
        """Antes se ignoraba en silencio: una errata como ``boundary_condition``
        dejaba el modelo sin apoyos sin avisar."""
        with pytest.raises(YamlValidationError, match=r"Secciones desconocidas: \['boundary_condition'\]"):
            _parsear(tmp_path, _BARRA + "boundary_condition:\n  - {node_id: 1, ux: 0.0}\n")

    def test_parametro_no_aceptado_por_un_material_termico(self, tmp_path):
        """La validación de parámetros por firma, que antes sólo tenía la
        familia mecánica, vale ahora para todas."""
        with pytest.raises(YamlValidationError, match="parámetro 'kk' no aceptado por 'ThermalConduction'"):
            _parsear(tmp_path, _PARED_TERMICA.replace("k: 50.0", "k: 50.0, kk: 1.0"))

    def test_elemento_termico_con_material_mecanico_falla_al_validar(self, tmp_path):
        """Cada elemento declara a qué familia apunta su ``material`` (P3).
        Antes el id se buscaba en las dos familias y un elemento térmico
        podía recibir un material mecánico con el mismo id."""
        texto = _PARED_TERMICA.replace(
            "thermal_materials:\n  - {id: 1, type: ThermalConduction, k: 50.0}\n",
            "materials:\n  - {id: 1, type: Elastic2D, E: 1.0, nu: 0.3}\n",
        )
        with pytest.raises(YamlValidationError,
                           match=r"referencia a material térmico inexistente \(material=1\): "
                                 r"no está declarado en 'thermal_materials'"):
            _parsear(tmp_path, texto)

    def test_falta_la_referencia_obligatoria(self, tmp_path):
        with pytest.raises(YamlValidationError, match="falta el campo obligatorio 'material'"):
            _parsear(tmp_path, _BARRA.replace(", material: 1", ""))

    def test_los_objetos_se_leen_por_el_nombre_de_su_seccion(self, tmp_path):
        parser = _parsear(tmp_path, _BARRA)
        assert parser.materials[1].E == 200.0e9
        assert parser.thermal_materials == {}
        with pytest.raises(AttributeError):
            parser.seccion_inexistente
