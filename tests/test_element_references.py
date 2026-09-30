"""Referencias de un elemento a familias de material (ADR 0020, P3 y P7).

Un elemento declara en ``REFERENCE_KWARGS`` qué parámetros de su constructor
son, en el YAML, el ``id`` de un objeto de una familia de material. El lector
comprueba que el id exista en la sección de *su* familia y le pasa el objeto
construido, tanto en el bloque ``elements`` (P3) como en la receta de un grupo
físico gmsh (P7). La referencia es obligatoria si el constructor no le da
valor por defecto.

La suite principal sólo lo probaba con la referencia única ``material``. Aquí
se prueba con un elemento de juguete de dos referencias: ``material``, a
``materials``, y ``resorte``, a una familia declarada fuera del programa
principal, como la declararía un módulo de usuario. El elemento y la familia se
registran en un fixture y se retiran al terminar, para no añadir clases a los
barridos de contrato ni secciones al resto de la sesión.
"""
import textwrap
from pathlib import Path

import pytest

from solidum.elements.solid_2d.tri3 import Tri3
from solidum.registry import ElementRegistry, MaterialRegistry, Registry
from solidum.utils.yaml_parser import YamlParser, YamlValidationError

PLACA = Path(__file__).resolve().parents[1] / "examples" / "placa.msh"   # triángulos
GRUPO = "Dominio_Placa"                                                  # grupo físico 2D de placa.msh

OPCIONAL = "TriConResorte_prueba"
OBLIGATORIO = "TriConResorteObligatorio_prueba"


class _Resorte:
    def __init__(self, k: float):
        self.k = k


@pytest.fixture(scope="module")
def resortes():
    """Registra la familia ``resortes_de_prueba`` y los dos elementos de
    juguete; los retira al terminar el módulo."""

    class ResortesDePrueba(Registry):
        YAML_SECTION = "resortes_de_prueba"
        YAML_LABEL = "resorte"

    ResortesDePrueba.register("Resorte", _Resorte)

    class TriConResorte(Tri3):
        """``Tri3`` con una segunda referencia, opcional, a la familia de prueba."""
        REFERENCE_KWARGS = {"material": MaterialRegistry, "resorte": ResortesDePrueba}

        def __init__(self, element_id, nodes, material, resorte=None, thickness=1.0):
            super().__init__(element_id, nodes, material, thickness=thickness)
            self.resorte = resorte

    class TriConResorteObligatorio(TriConResorte):
        """La misma segunda referencia, ahora sin valor por defecto."""

        def __init__(self, element_id, nodes, material, resorte, thickness=1.0):
            super().__init__(element_id, nodes, material, resorte, thickness=thickness)

    elementos = {OPCIONAL: TriConResorte, OBLIGATORIO: TriConResorteObligatorio}
    for nombre, cls in elementos.items():
        ElementRegistry.register(nombre, cls)
    try:
        yield ResortesDePrueba
    finally:
        for nombre in elementos:
            ElementRegistry._items.pop(nombre, None)
        Registry._families.remove(ResortesDePrueba)


def _parsear(tmp_path, texto: str):
    ruta = tmp_path / "modelo.yaml"
    ruta.write_text(textwrap.dedent(texto), encoding="utf-8")
    parser = YamlParser(str(ruta))
    return parser, parser.parse()


_FAMILIAS = """\
materials:
  - {id: 1, type: Elastic2D, E: 2.0e11, nu: 0.3, hypothesis: plane_stress}
resortes_de_prueba:
  - {id: 7, type: Resorte, k: 3.0}
"""

_TRIANGULO = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [0.0, 1.0]}
""" + _FAMILIAS


def _en_linea(tipo: str, referencias: str) -> str:
    return (_TRIANGULO + "elements:\n"
            f"  - {{id: 1, type: {tipo}, nodes: [1, 2, 3]{referencias}}}\n")


class TestBloqueElements:
    """P3: referencias del bloque ``elements``."""

    def test_las_dos_referencias_se_resuelven(self, tmp_path, resortes):
        parser, dom = _parsear(tmp_path, _en_linea(OPCIONAL, ", material: 1, resorte: 7"))
        elem = dom.elements[1]
        assert elem.material is parser.materials[1]
        assert elem.resorte is parser.resortes_de_prueba[7]
        assert elem.resorte.k == 3.0

    def test_la_referencia_opcional_puede_faltar(self, tmp_path, resortes):
        _, dom = _parsear(tmp_path, _en_linea(OPCIONAL, ", material: 1"))
        assert dom.elements[1].resorte is None

    def test_id_inexistente_en_la_segunda_familia(self, tmp_path, resortes):
        with pytest.raises(YamlValidationError,
                           match=r"referencia a resorte inexistente \(resorte=99\): "
                                 r"no está declarado en 'resortes_de_prueba'"):
            _parsear(tmp_path, _en_linea(OPCIONAL, ", material: 1, resorte: 99"))

    def test_cada_id_se_busca_en_la_familia_de_su_referencia(self, tmp_path, resortes):
        """El id 1 es un material y el 7 un resorte: cruzados, fallan los dos."""
        with pytest.raises(YamlValidationError) as exc:
            _parsear(tmp_path, _en_linea(OPCIONAL, ", material: 7, resorte: 1"))
        texto = str(exc.value)
        assert ("referencia a material inexistente (material=7): "
                "no está declarado en 'materials'") in texto
        assert ("referencia a resorte inexistente (resorte=1): "
                "no está declarado en 'resortes_de_prueba'") in texto

    def test_falta_la_referencia_obligatoria(self, tmp_path, resortes):
        with pytest.raises(YamlValidationError, match="falta el campo obligatorio 'resorte'"):
            _parsear(tmp_path, _en_linea(OBLIGATORIO, ", material: 1"))
        parser, dom = _parsear(tmp_path, _en_linea(OBLIGATORIO, ", material: 1, resorte: 7"))
        assert dom.elements[1].resorte is parser.resortes_de_prueba[7]


def _con_malla(recetas: str) -> str:
    return f"mesh: {PLACA.as_posix()}\n" + _FAMILIAS + textwrap.dedent(recetas)


@pytest.mark.skipif(not PLACA.is_file(), reason="placa.msh no disponible")
class TestRecetaGmsh:
    """P7: referencias de la receta de un grupo físico gmsh. ``material`` sale
    de ``mesh_material`` (1 por omisión); ``resorte``, del propio grupo."""

    def test_grupo_con_las_dos_referencias(self, tmp_path, resortes):
        parser, dom = _parsear(tmp_path, _con_malla(f"""\
            mesh_physical_groups:
              {GRUPO}: {{element: {OPCIONAL}, resorte: 7}}
            """))
        elementos = list(dom.elements.values())
        assert {type(e).__name__ for e in elementos} == {"TriConResorte"}
        assert all(e.resorte is parser.resortes_de_prueba[7] for e in elementos)
        assert all(e.material is parser.materials[1] for e in elementos)

    def test_grupo_sin_la_referencia_opcional(self, tmp_path, resortes):
        _, dom = _parsear(tmp_path, _con_malla(f"""\
            mesh_physical_groups:
              {GRUPO}: {{element: {OPCIONAL}}}
            """))
        assert all(e.resorte is None for e in dom.elements.values())

    def test_el_grupo_hereda_mesh_element_y_anade_la_referencia(self, tmp_path, resortes):
        parser, dom = _parsear(tmp_path, _con_malla(f"""\
            mesh_element: {OPCIONAL}
            mesh_physical_groups:
              {GRUPO}: {{resorte: 7}}
            """))
        assert all(e.resorte is parser.resortes_de_prueba[7] for e in dom.elements.values())

    def test_id_inexistente_en_la_segunda_familia(self, tmp_path, resortes):
        with pytest.raises(YamlValidationError,
                           match=rf"mesh_physical_groups\['{GRUPO}'\]: resorte=99 no existe "
                                 r"en 'resortes_de_prueba' \(ids declarados: \[7\]\)"):
            _parsear(tmp_path, _con_malla(f"""\
                mesh_physical_groups:
                  {GRUPO}: {{element: {OPCIONAL}, resorte: 99}}
                """))

    def test_falta_la_referencia_obligatoria(self, tmp_path, resortes):
        with pytest.raises(YamlValidationError,
                           match=rf"mesh_physical_groups\['{GRUPO}'\]: falta 'resorte' "
                                 rf"para el elemento '{OBLIGATORIO}'"):
            _parsear(tmp_path, _con_malla(f"""\
                mesh_physical_groups:
                  {GRUPO}: {{element: {OBLIGATORIO}}}
                """))
