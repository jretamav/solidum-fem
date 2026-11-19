"""Pieza P10 (ADR 0042): cada elemento conoce el dominio que lo contiene.

Invariantes que protegen estas pruebas:
- ``element.domain`` es el dominio al que lo añadió ``add_element`` o
  ``replace_element`` mientras siga vivo y el elemento siga siendo el suyo con
  ese id; si no, ``None``.
- La referencia no se serializa: tras ``pickle`` o ``copy.deepcopy`` del
  dominio, cada elemento apunta a su propio dominio (el copiado), nunca al
  original.
- Un elemento no puede estar en dos dominios vivos.
- En ``solidum/`` sólo ``Domain`` cambia ``domain.elements`` (un solo camino
  de mutación; si no, ``element.domain`` dejaría de ser fiable).
"""
from __future__ import annotations

import ast
import copy
import gc
import pickle
from pathlib import Path

import pytest

from solidum.core.domain import Domain
from solidum.core.element import Element
from solidum.elements.solid_2d import Quad4, Tri3
from solidum.materials.elastic_2d import Elastic2D

ROOT = Path(__file__).resolve().parents[1]
MAT = Elastic2D(E=1000.0, nu=0.3)


def _two_triangles():
    dom = Domain()
    for i, c in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        dom.add_node(i, [float(v) for v in c])
    n = dom.nodes
    dom.add_element(Tri3(1, [n[1], n[2], n[3]], MAT))
    dom.add_element(Tri3(2, [n[1], n[3], n[4]], MAT))
    return dom


def test_add_element_enlaza_con_el_dominio():
    dom = _two_triangles()
    assert all(el.domain is dom for el in dom.elements.values())


def test_un_elemento_suelto_no_tiene_dominio():
    dom = Domain()
    a, b, c = (dom.add_node(i, x) for i, x in ((1, [0.0, 0.0]), (2, [1.0, 0.0]), (3, [0.0, 1.0])))
    assert Tri3(1, [a, b, c], MAT).domain is None


def test_pickle_ida_y_vuelta_enlaza_con_el_dominio_nuevo():
    dom = _two_triangles()
    back = pickle.loads(pickle.dumps(dom))
    assert all(el.domain is back for el in back.elements.values())
    assert all(el.domain is dom for el in dom.elements.values())


def test_deepcopy_enlaza_la_copia_con_la_copia():
    dom = _two_triangles()
    twin = copy.deepcopy(dom)
    assert all(el.domain is twin for el in twin.elements.values())
    assert all(el.domain is dom for el in dom.elements.values())
    assert twin.elements[1] is not dom.elements[1]


def test_replace_element_enlaza_el_nuevo_y_suelta_el_viejo():
    dom = _two_triangles()
    old = dom.elements[1]
    new = Tri3(1, list(old.nodes), Elastic2D(E=2000.0, nu=0.3))
    dom.replace_element(new)
    assert dom.elements[1] is new and new.domain is dom
    assert old.domain is None
    with pytest.raises(ValueError, match="no hay ningún elemento 7"):
        dom.replace_element(Tri3(7, list(old.nodes), MAT))


def test_elemento_colado_sin_add_element_no_tiene_dominio():
    """Un elemento puesto a mano en ``elements`` no queda enlazado; el que
    había con ese id deja de reconocer el dominio como suyo."""
    dom = _two_triangles()
    old = dom.elements[2]
    intruso = Tri3(2, list(old.nodes), MAT)
    dom.elements[2] = intruso                       # a propósito: el camino prohibido
    assert intruso.domain is None
    assert old.domain is None


def test_un_elemento_no_puede_estar_en_dos_dominios_vivos():
    dom = _two_triangles()
    other = Domain()
    with pytest.raises(ValueError, match="otro dominio"):
        other.add_element(dom.elements[1])
    el = dom.elements[1]
    del dom
    gc.collect()
    other.add_element(el)                            # el primero ya no existe
    assert el.domain is other


def test_subclase_sin_registrar_tambien_se_enlaza():
    class Propio(Quad4):
        pass
    dom = Domain()
    nodes = [dom.add_node(i, [float(x), float(y)])
             for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1)]
    el = Propio(1, nodes, MAT)
    dom.add_element(el)
    assert el.domain is dom
    twin = copy.deepcopy(dom)          # local: no se puede pickle; deepcopy sí
    assert twin.elements[1].domain is twin


def test_la_referencia_al_dominio_es_un_atributo_de_clase_por_omision():
    assert Element._domain_ref is None


def _mutaciones_de_elements(path: Path) -> list[str]:
    """Escrituras en ``<algo>.elements[...]`` y llamadas que cambian el
    diccionario (``pop``, ``update``, ``setdefault``, ``clear``, ``popitem``,
    ``__setitem__``, ``__delitem__``)."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    out = []

    def es_elements(node):
        return isinstance(node, ast.Attribute) and node.attr == "elements"
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign, ast.Delete)):
            targets = node.targets if isinstance(node, (ast.Assign, ast.Delete)) else [node.target]
            for t in targets:
                if isinstance(t, ast.Subscript) and es_elements(t.value):
                    out.append(f"{path.relative_to(ROOT)}:{node.lineno}")
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr in {"pop", "update", "setdefault", "clear", "popitem",
                                       "__setitem__", "__delitem__"}
                and es_elements(node.func.value)):
            out.append(f"{path.relative_to(ROOT)}:{node.lineno}")
    return out


def test_solo_domain_cambia_el_diccionario_de_elementos():
    domain_py = ROOT / "solidum" / "core" / "domain.py"
    found = []
    for path in sorted((ROOT / "solidum").rglob("*.py")):
        if path != domain_py:
            found += _mutaciones_de_elements(path)
    assert not found, ("Escrituras directas en domain.elements; use Domain.add_element o "
                       f"Domain.replace_element (ADR 0042): {found}")


def test_el_detector_de_mutaciones_las_ve(tmp_path):
    """El control anterior no es vacuo: detecta las formas prohibidas."""
    src = tmp_path / "x.py"
    src.write_text("d.elements[1] = e\ndel d.elements[2]\nd.elements.pop(3)\nd.elements.update({})\n"
                   "self.elements = []\n", encoding="utf-8")
    global ROOT
    old = ROOT
    ROOT = tmp_path
    try:
        assert len(_mutaciones_de_elements(src)) == 4
    finally:
        ROOT = old


def test_copia_superficial_no_roba_el_enlace_al_original():
    """``copy.copy`` comparte los elementos con el original, que los conserva
    aunque la copia desaparezca."""
    dom = _two_triangles()
    shallow = copy.copy(dom)
    del shallow
    gc.collect()
    assert all(el.domain is dom for el in dom.elements.values())


def test_add_y_replace_cambian_la_huella_de_topologia():
    dom = _two_triangles()
    v = dom.topology_version
    old = dom.elements[1]
    dom.replace_element(Tri3(1, list(old.nodes), MAT))
    assert dom.topology_version == v + 1
    n = dom.nodes
    dom.add_element(Tri3(3, [n[2], n[3], n[1]], MAT))
    assert dom.topology_version == v + 2


@pytest.mark.parametrize("batch", [True, False])
def test_el_ensamblador_ve_un_elemento_sustituido(batch):
    """Tras ``replace_element`` el mismo ``Assembler`` ensambla el elemento
    nuevo (por lotes también): la huella de topología incluye
    ``Domain.topology_version``."""
    import numpy as np

    from solidum.math.assembly import Assembler

    def model():
        dom = Domain()
        pts = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]
        nd = [dom.add_node(i, [float(x), float(y)]) for i, (x, y) in enumerate(pts, 1)]
        dom.add_element(Quad4(1, [nd[0], nd[1], nd[4], nd[3]], MAT))
        dom.add_element(Quad4(2, [nd[1], nd[2], nd[5], nd[4]], MAT))
        dom.generate_equation_numbers(verbose=False)
        return dom

    def stiffness(asm, dom):
        K, _ = asm.assemble_non_linear_system(np.zeros(dom.total_dofs))
        return K.toarray() if hasattr(K, "toarray") else np.asarray(K)

    dom = model()
    asm = Assembler(dom, batch=batch)
    K1 = stiffness(asm, dom)
    old = dom.elements[2]
    dom.replace_element(Quad4(2, list(old.nodes), Elastic2D(E=5000.0, nu=0.3)))
    K2 = stiffness(asm, dom)
    K_ref = stiffness(Assembler(dom, batch=batch), dom)
    assert np.max(np.abs(K2 - K1)) > 1.0
    np.testing.assert_allclose(K2, K_ref, rtol=0, atol=1e-9 * np.max(np.abs(K_ref)))
