"""Contorno libre de una malla (``solidum.mesh.free_boundary``, ADR 0042).

Invariantes: un lado (2D) o cara (3D) es libre si lo declara un solo elemento
en ``EDGE_NODES``/``FACE_NODES``; el lado común de dos elementos de tipos
distintos no es libre; con el filtro de campo, los elementos sin esos grados
de libertad (térmicos) no tapan el contorno del sólido.
"""
from __future__ import annotations

import pytest

from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4, Tri3
from solidum.elements.solid_3d import Hex8
from solidum.elements.thermal import Quad4Thermal
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.mesh import free_boundary

M2 = Elastic2D(E=1000.0, nu=0.3)


def _nodes(dom, pts):
    return [dom.add_node(i, [float(v) for v in p]) for i, p in enumerate(pts, 1)]


def _sets(cells):
    return {frozenset(c) for c in cells}


def test_quad4_y_tri3_lado_comun_no_libre():
    dom = Domain()
    n = _nodes(dom, [(0, 0), (1, 0), (1, 1), (0, 1), (2, 0.5)])
    dom.add_element(Quad4(1, [n[0], n[1], n[2], n[3]], M2))
    dom.add_element(Tri3(2, [n[1], n[4], n[2]], M2))
    free = free_boundary(dom.elements.values(), 2)
    assert _sets(free) == {frozenset(s) for s in ((1, 2), (3, 4), (4, 1), (2, 5), (5, 3))}
    assert frozenset((2, 3)) not in _sets(free)
    # Orden local del elemento que lo tiene (orientación antihoraria).
    assert (2, 5) in free and (5, 3) in free


def test_elementos_termicos_no_tapan_el_contorno_mecanico():
    dom = Domain()
    n = _nodes(dom, [(0, 0), (1, 0), (1, 1), (0, 1)])
    dom.add_element(Quad4(1, n, M2))
    dom.add_element(Quad4Thermal(2, n, ThermalConduction(k=1.0, c=1.0, density=1.0, dim=2)))
    assert free_boundary(dom.elements.values(), 2) == []            # sin filtro, todo compartido
    assert len(free_boundary(dom.elements.values(), 2, dofs=("ux", "uy"))) == 4


def test_caras_libres_de_un_hexaedro():
    dom = Domain()
    pts = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    dom.add_element(Hex8(1, _nodes(dom, pts), Elastic3D(E=1000.0, nu=0.3)))
    faces = free_boundary(dom.elements.values(), 3)
    assert len(faces) == 6 and all(len(f) == 4 for f in faces)
    assert free_boundary(dom.elements.values(), 2) == []            # Hex8 no declara lados


def test_dimension_invalida():
    with pytest.raises(ValueError, match="dim"):
        free_boundary([], 1)
