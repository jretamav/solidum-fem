"""Clasificación única de la magnitud por punto de un elemento (ADR 0043 §3).

La usan la exportación VTK, la historia de la salida, la calidad de malla, la
herramienta de specs, el camino por lotes y el barrido de contrato. Antes cada
uno la deducía por presencia de atributos, y un elemento de interfaz con
``STRAIN_DIM`` habría pasado por sólido en unos sitios y no en otros.
"""
from pathlib import Path

import pytest

import solidum  # noqa: F401 — puebla los registros
from solidum import (Bond, BondSlipMultilinear, CohesiveDamageMixedMode, Elastic1D, Elastic2D,
                     embed_reinforcement)
from solidum.bc.bond_interface import BondSlipInterface
from solidum.core.node import Node
from solidum.core.point_quantity import FLUX, JUMP, STRAIN, PointQuantity, point_quantity
from solidum.elements.interface import Interface2D4, Interface3D6, LocalJumpElement
from solidum.mesh import structured_mesh
from solidum.mesh.quality import mesh_quality
from solidum.registry import ElementRegistry
from solidum.utils.vtk_exporter import cell_type_for

from _programa_principal import nombres_programa_principal


def test_todo_elemento_del_programa_principal_tiene_un_tipo():
    for nombre in nombres_programa_principal(ElementRegistry):
        q = point_quantity(ElementRegistry.get(nombre))
        assert q.kind in (STRAIN, FLUX, JUMP) and isinstance(q.dim, int), (nombre, q)


@pytest.mark.parametrize("nombre, esperado", [
    ("Quad4", PointQuantity(STRAIN, 3)), ("Hex8", PointQuantity(STRAIN, 6)),
    ("Truss2D", PointQuantity(STRAIN, 1)), ("Frame3D", PointQuantity(STRAIN, 1)),
    ("Quad4Thermal", PointQuantity(FLUX, 2)), ("Hex8Thermal", PointQuantity(FLUX, 3)),
    ("Interface2D4", PointQuantity(JUMP, 2)), ("Interface2D6", PointQuantity(JUMP, 2)),
    ("Interface3D6", PointQuantity(JUMP, 3)),
    ("Interface3D8", PointQuantity(JUMP, 3)),
])
def test_tipo_de_cada_familia(nombre, esperado):
    assert point_quantity(ElementRegistry.get(nombre)) == esperado


def test_todo_elemento_de_salto_hereda_la_base_comun():
    """``JUMP`` es un contrato completo (revisión del arquitecto del
    2026-11-26): la historia, el VTK y las energías preguntan sólo por la
    clasificación, así que todo elemento de salto registrado —también de un
    módulo de usuario— hereda ``LocalJumpElement`` (salto, componentes,
    registro por punto, ``force_matrix``)."""
    from solidum.elements.interface import LocalJumpElement
    for nombre in ElementRegistry.names():
        cls = ElementRegistry.get(nombre)
        if point_quantity(cls).kind == JUMP:
            assert issubclass(cls, LocalJumpElement), nombre


def test_el_salto_tiene_prioridad_y_solo_valen_enteros():
    class Ambos:
        JUMP_DIM, STRAIN_DIM = 2, 3

    class SinEntero:
        STRAIN_DIM, FLUX_DIM, JUMP_DIM = None, None, True

    assert point_quantity(Ambos) == PointQuantity(JUMP, 2)
    assert point_quantity(SinEntero) == PointQuantity(None, None)
    assert point_quantity(SinEntero).attribute is None
    assert point_quantity(Ambos).attribute == "JUMP_DIM"


@pytest.mark.parametrize("dim", [2, 3])
def test_la_adherencia_es_un_salto_y_comparte_la_base(dim):
    """``BondSlipInterface`` declara el salto (no ``STRAIN_DIM``, que la haría
    pasar por sólido en la historia) y usa el bucle por punto de la base común
    (ADR 0043 §2)."""
    if dim == 2:
        mesh = structured_mesh("Tri3", (2, 1), size=(1.0, 0.3),
                               material=Elastic2D(E=30e9, nu=0.2), thickness=0.3)
        p, q = [0.0, 0.15], [1.0, 0.15]
    else:
        from solidum.materials.elastic_3d import Elastic3D
        mesh = structured_mesh("Hex8", (2, 1, 1), size=(1.0, 0.3, 0.3),
                               material=Elastic3D(E=30e9, nu=0.2))
        p, q = [0.0, 0.15, 0.15], [1.0, 0.15, 0.15]
    d = mesh.domain
    bar = embed_reinforcement(d, [p, q], 2e-4, Elastic1D(E=200e9), element="Truss2D" if dim == 2
                              else "Truss3D",
                              bond=Bond(BondSlipMultilinear(points=[[1e-3, 5e6]], k_n=1e12), 0.016))
    el = d.elements[bar.interface_ids[0]]
    assert isinstance(el, BondSlipInterface) and isinstance(el, LocalJumpElement)
    assert point_quantity(el) == PointQuantity(JUMP, dim)
    assert el.STRAIN_DIM is None
    assert type(el).compute_element_state is LocalJumpElement.compute_element_state
    assert type(el).compute_mass_matrix is LocalJumpElement.compute_mass_matrix


def _interfaz_2d():
    ley = CohesiveDamageMixedMode(sigma_t0=3e6, G_f=100.0, K_e=1e13, beta=1.0, softening="linear")
    nodos = [Node(1, [0.0, 0.0]), Node(2, [1.0, 0.0]), Node(3, [1.0, 0.0]), Node(4, [0.0, 0.0])]
    return Interface2D4(1, nodos, ley)


def test_celda_vtk_de_las_interfaces():
    """Con la tabla por número de nodos, la de 6 nodos se omitiría y la de 8
    sería un hexaedro degenerado por casualidad: la declara el elemento."""
    assert cell_type_for(_interfaz_2d()) == "quad"
    ley = CohesiveDamageMixedMode(sigma_t0=3e6, G_f=100.0, K_e=1e13, beta=1.0, softening="linear")
    cara = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    nodos = [Node(i + 1, list(x)) for i, x in enumerate(cara + cara)]
    assert cell_type_for(Interface3D6(1, nodos, ley)) == "wedge"


def test_la_calidad_de_malla_no_mide_interfaces():
    d = structured_mesh("Quad4", (1, 1), size=(1.0, 1.0), material=Elastic2D(E=1.0, nu=0.0),
                        thickness=1.0).domain
    ley = CohesiveDamageMixedMode(sigma_t0=3e6, G_f=100.0, K_e=1e13, beta=1.0, softening="linear")
    nodos = [d.add_node(100 + i, x) for i, x in
             enumerate([[0.0, 0.0], [1.0, 0.0], [1.0, 0.0], [0.0, 0.0]])]
    d.add_element(Interface2D4(99, nodos, ley))
    q = mesh_quality(d)
    assert 99 not in q.element_ids and len(q.element_ids) == 1


def test_la_spec_debe_declarar_el_campo_de_su_clase(tmp_path):
    """La herramienta de specs compara con el atributo que da la
    clasificación: una spec de interfaz que dijera ``field: displacement``
    pediría ``STRAIN_DIM`` a una clase que declara ``JUMP_DIM``."""
    from solidum.tools.spec import cross_check_with_registry, parse_spec
    src = Path(__file__).parents[1] / "docs" / "specs" / "Interface2D4.md"
    texto = src.read_text(encoding="utf-8")
    texto = texto.replace("status: draft", "status: implemented", 1)
    texto = texto.replace("field: interface", "field: displacement", 1)
    texto = texto.replace("jump_dim: 2", "strain_dim: 2", 1)
    p = tmp_path / "Interface2D4.md"
    p.write_text(texto + "\n## Implementación\n\nArchivo: x.py\n", encoding="utf-8")
    errores = cross_check_with_registry(parse_spec(p))
    assert any("pide STRAIN_DIM" in e and "declara JUMP_DIM" in e for e in errores), errores
