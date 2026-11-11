"""Lector de Gmsh: todas las celdas 2D y 3D, orden de nodos y orientación.

Mallas pequeñas versionadas en ``tests/data/gmsh/`` (generadas con
``tests/data/gmsh/generar.py``; las pruebas no necesitan Gmsh salvo la que
contrasta la tabla de referencia con el propio Gmsh):

- placa de 2 × 1 (espesor 0.1): Tri3, Tri6, Quad4, Quad8, Quad9 (las de
  cuadriláteros pueden traer algún triángulo de la recombinación), con una
  barra incrustada «Barra» y los bordes «Izquierda» y «Derecha»; Tri6 y Quad8
  también con la superficie en sentido horario;
- bloque de 2 × 1 × 1: Tet4, Tet10, Hex8, Hex20, Hex27 y un Tet10 con las
  celdas invertidas.

Prueba de la parcela: con un desplazamiento lineal ``u = A·x`` todo punto de
integración debe dar el esfuerzo ``C:ε(A)``; si el orden de los nodos o la
orientación fueran incorrectos, no lo daría. Además ``Σ u_e·F_int,e = ε:C:ε·V``
devuelve el volumen (área por espesor) de la malla.
"""
import textwrap
from pathlib import Path

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.utils.gmsh_parser import _CELL_REFERENCE, GmshParser
from solidum.utils.yaml_parser import YamlParser

DATA = Path(__file__).resolve().parent / "data" / "gmsh"
E, NU, T = 200.0, 0.3, 0.1

CASOS_2D = {
    "tri3": {"Tri3"}, "tri6": {"Tri6"}, "quad4": {"Quad4", "Tri3"}, "quad8": {"Quad8", "Tri6"},
    "quad9": {"Quad9", "Tri6"}, "tri6_horaria": {"Tri6"}, "quad8_horaria": {"Quad8"},
}
CASOS_3D = {
    "tet4": {"Tet4"}, "tet10": {"Tet10"}, "hex8": {"Hex8"}, "hex20": {"Hex20"},
    "hex27": {"Hex27"}, "tet10_invertida": {"Tet10"},
}


def _parse(nombre, material, **opt):
    p = GmshParser(str(DATA / f"{nombre}.msh"))
    spec = {"element": None, "kwargs": {"material": material},
            "optional": {"thickness": T, "quadrature": None, **opt}}
    return p, p.parse(spec)


def _parcela(dom, dim, C):
    rng = np.random.default_rng(7)
    A = rng.normal(size=(dim, dim)) * 1e-3
    if dim == 2:
        eps = np.array([A[0, 0], A[1, 1], A[0, 1] + A[1, 0]])
    else:
        eps = np.array([A[0, 0], A[1, 1], A[2, 2], A[0, 1] + A[1, 0], A[1, 2] + A[2, 1],
                        A[0, 2] + A[2, 0]])
    sig = C @ eps
    trabajo, peor = 0.0, 0.0
    for el in dom.elements.values():
        u = np.concatenate([A @ np.asarray(n.coordinates, dtype=float) for n in el.nodes])
        _, F = el.compute_element_state(u)
        trabajo += float(u @ F)
        for s in el.state.stresses_trial:
            peor = max(peor, float(np.abs(np.asarray(s) - sig).max() / np.abs(sig).max()))
    return peor, trabajo / float(eps @ sig)


@pytest.mark.parametrize("nombre", sorted(CASOS_2D))
def test_celdas_2d(nombre):
    mat = Elastic2D(E=E, nu=NU, hypothesis="plane_stress")
    parser, dom = _parse(nombre, mat)
    assert {type(e).__name__ for e in dom.elements.values()} == CASOS_2D[nombre]
    peor, vol = _parcela(dom, 2, mat.C)
    assert peor < 1e-10, nombre
    assert vol == pytest.approx(2.0 * 1.0 * T, rel=1e-10)
    assert (parser.n_reoriented > 0) == nombre.endswith("_horaria")


@pytest.mark.parametrize("nombre", sorted(CASOS_3D))
def test_celdas_3d(nombre):
    mat = Elastic3D(E=E, nu=NU)
    parser, dom = _parse(nombre, mat)
    assert {type(e).__name__ for e in dom.elements.values()} == CASOS_3D[nombre]
    assert all(len(n.coordinates) == 3 for n in dom.nodes.values())
    peor, vol = _parcela(dom, 3, mat.C)
    assert peor < 1e-10, nombre
    assert vol == pytest.approx(2.0, rel=1e-10)
    assert (parser.n_reoriented > 0) == nombre.endswith("_invertida")


@pytest.mark.parametrize("nombre,grupo,entidad", [
    ("quad8", "Derecha", "EDGE_NODES"), ("tri6", "Izquierda", "EDGE_NODES"),
    ("quad8_horaria", "Derecha", "EDGE_NODES"),
    ("hex20", "Cara_x1", "FACE_NODES"), ("tet10", "Cara_x0", "FACE_NODES"),
    ("tet10_invertida", "Cara_x0", "FACE_NODES"), ("hex27", "Cara_x0", "FACE_NODES"),
])
def test_celdas_de_contorno_casan_con_un_borde_o_cara(nombre, grupo, entidad):
    """Las cargas por grupo (ADR 0023) casan cada celda de contorno con el borde
    o la cara de un único elemento, también en mallas cuadráticas."""
    mat = Elastic3D(E=E, nu=NU) if entidad == "FACE_NODES" else Elastic2D(E=E, nu=NU)
    _, dom = _parse(nombre, mat)
    caras = {}
    for el in dom.elements.values():
        for local in getattr(el, entidad):
            caras.setdefault(frozenset(el.nodes[j].id for j in local), []).append(el.id)
    celdas = dom.physical_group_cells[grupo]
    assert celdas
    assert all(len(caras.get(frozenset(c), [])) == 1 for c in celdas)


def test_tabla_de_referencia_igual_a_gmsh():
    """``_CELL_REFERENCE`` es la de Gmsh con la permutación de meshio aplicada."""
    gmsh = pytest.importorskip("gmsh")
    from meshio.gmsh import common
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        for gtype, mtype in common._gmsh_to_meshio_type.items():
            if mtype not in _CELL_REFERENCE:
                continue
            _, dim, _, nn, loc, _ = gmsh.model.mesh.getElementProperties(gtype)
            loc = np.asarray(loc).reshape(nn, dim)
            perm = common._gmsh_to_meshio_order(mtype, np.arange(nn)[None, :])[0]
            np.testing.assert_allclose(np.array(_CELL_REFERENCE[mtype][1]), loc[perm], atol=1e-12,
                                       err_msg=mtype)
    finally:
        gmsh.finalize()


def _yaml(tmp_path, nombre, extra=""):
    ruta = tmp_path / "modelo.yaml"
    ruta.write_text(f"mesh: {(DATA / (nombre + '.msh')).as_posix()}\n" + textwrap.dedent("""\
        materials:
          - {id: 1, type: Elastic2D, E: 200.0, nu: 0.3, hypothesis: plane_stress}
          - {id: 2, type: Elastic1D, E: 2000.0}
        mesh_thickness: 0.1
        """) + textwrap.dedent(extra), encoding="utf-8")
    return YamlParser(str(ruta)).parse()


def test_barra_incrustada_como_elementos_de_linea(tmp_path):
    """Un grupo de líneas con 'element' explícito da barras que comparten los
    nodos del sólido."""
    dom = _yaml(tmp_path, "tri3", """\
        mesh_physical_groups:
          Barra: {element: Truss2D, material: 2, A: 0.01}
        """)
    barras = [e for e in dom.elements.values() if type(e).__name__ == "Truss2D"]
    assert len(barras) == 4
    nodos_solido = {n.id for e in dom.elements.values() if type(e).__name__ == "Tri3" for n in e.nodes}
    assert all(n.id in nodos_solido for b in barras for n in b.nodes)
    assert all(isinstance(b.material, Elastic1D) for b in barras)


def test_grupo_de_contorno_sin_elemento_explicito_no_crea_elementos(tmp_path):
    """El 'element' heredado de 'mesh_element' no convierte el contorno en elementos."""
    dom = _yaml(tmp_path, "tri3", """\
        mesh_element: Tri3
        mesh_physical_groups:
          Izquierda: {thickness: 0.1}
        """)
    assert {type(e).__name__ for e in dom.elements.values()} == {"Tri3"}
    assert len(dom.elements) == 32


def test_cuadratura_por_omision_de_cada_elemento(tmp_path):
    """Sin 'mesh_quadrature' cada elemento usa su regla (antes se forzaba 2x2,
    que no sirve a Tri6)."""
    dom = _yaml(tmp_path, "quad8")
    claves = {type(e).__name__: e.quadrature_key for e in dom.elements.values()}
    assert claves["Quad8"] == "3x3"
    assert claves["Tri6"] == "tri_3"
