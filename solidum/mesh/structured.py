"""Mallas estructuradas de un bloque transformado (ADR 0022).

Se malla el bloque de referencia ``[0, 1]^d`` en ``divisions`` celdas por
dirección, se colocan en cada celda los nodos del elemento a partir de sus
coordenadas naturales (``REFERENCE_NODES``), se funden los compartidos y se
lleva cada nodo a su posición física::

    X = origin + ξ·size,     x = mapping(X)   (identidad si no se da).

Con una ``mapping`` curva, los nodos intermedios de un elemento cuadrático
quedan sobre la curva: los bordes y las caras son curvos de verdad, también
en el contorno. Con una ``mapping`` suave, al refinar los elementos tienden a
paralelogramos (mallas asintóticamente afines).

Celdas (``REFERENCE_CELL`` del elemento):

- ``"quad"``: una por celda; ``"tri"``: dos por celda, separadas por la
  diagonal de ``(0, 0)`` a ``(1, 1)``;
- ``"hex"``: una por celda; ``"tet"``: seis por celda (subdivisión de Kuhn
  alrededor de la diagonal de ``(0, 0, 0)`` a ``(1, 1, 1)``; conforme entre
  celdas vecinas porque todas parten sus caras por la misma diagonal).

Los lados del bloque se llaman ``xmin``, ``xmax``, ``ymin``, ``ymax`` (y
``zmin``, ``zmax`` en 3D). :attr:`StructuredMesh.boundary` da los bordes
(2D, ``EDGE_NODES``) o caras (3D, ``FACE_NODES``) de elemento sobre cada lado
y :meth:`StructuredMesh.nodes_on`, sus nodos.

Para geometrías reales, Gmsh (``solidum.utils.gmsh_parser``); esto es para
pruebas, estudios de convergencia y modelos de geometría simple.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import permutations
from typing import Callable

import numpy as np

from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.registry import ElementRegistry

SIDES_2D = ("xmin", "xmax", "ymin", "ymax")
SIDES_3D = SIDES_2D + ("zmin", "zmax")


def _kuhn_tets():
    """Seis tetraedros de la celda unidad, todos con volumen positivo."""
    tets = []
    for a, b, _ in permutations(range(3)):
        e_a, e_b = np.eye(3)[a], np.eye(3)[b]
        v = [np.zeros(3), e_a, e_a + e_b, np.ones(3)]
        if np.linalg.det(np.array([v[1], v[2], v[3]])) < 0.0:
            v[1], v[2] = v[2], v[1]
        tets.append(tuple(tuple(p) for p in v))
    return tuple(tets)


# Vértices de las celdas por elemento en la celda unidad. En los símplices,
# en el orden de los vértices del elemento de referencia; en quad y hex sólo
# importa el origen.
_CELL_VERTICES = {
    "quad": (((0.0, 0.0),),),
    "tri": (((0.0, 0.0), (1.0, 0.0), (1.0, 1.0)),
            ((0.0, 0.0), (1.0, 1.0), (0.0, 1.0))),
    "hex": (((0.0, 0.0, 0.0),),),
    "tet": _kuhn_tets(),
}
_DIM = {"quad": 2, "tri": 2, "hex": 3, "tet": 3}

# Tolerancia para decidir que un nodo está sobre un lado del bloque, en las
# coordenadas de referencia (que son fracciones exactas de 1/divisions).
_ON_SIDE_TOL = 1e-12


def _in_cell(cell: str, vertices, xi: np.ndarray) -> np.ndarray:
    """Punto de la celda unidad que corresponde a las coordenadas naturales
    ``xi`` del elemento (``[-1, 1]^d`` o coordenadas de volumen/área)."""
    v = np.asarray(vertices)
    if cell in ("quad", "hex"):
        return v[0] + (xi + 1.0) / 2.0
    return v[0] + (v[1:] - v[0]).T @ xi


def _side(side: str, dim: int):
    sides = SIDES_2D if dim == 2 else SIDES_3D
    if side not in sides:
        raise ValueError(f"Lado {side!r} desconocido; lados: {', '.join(sides)}.")
    return "xyz".index(side[0]), (0.0 if side.endswith("min") else 1.0)


@dataclass
class StructuredMesh:
    """Malla generada por :func:`structured_mesh`.

    Attributes
    ----------
    domain : Domain
        Nodos y elementos (sin numeración de ecuaciones: se llama a
        ``generate_equation_numbers`` tras imponer las condiciones de
        contorno).
    boundary : dict[str, list[tuple[Element, int]]]
        Por lado del bloque, los pares (elemento, índice de borde en
        ``EDGE_NODES`` en 2D o de cara en ``FACE_NODES`` en 3D) que lo forman;
        listos para ``compute_edge_traction`` / ``compute_face_traction``.
    """
    domain: Domain
    divisions: tuple
    boundary: dict = field(default_factory=dict)
    reference: dict = field(default_factory=dict, repr=False)

    @property
    def dim(self) -> int:
        return len(self.divisions)

    @property
    def sides(self) -> tuple:
        return SIDES_2D if self.dim == 2 else SIDES_3D

    def nodes_on(self, side: str) -> list[Node]:
        """Nodos sobre el lado ``side`` del bloque, en orden de creación."""
        axis, value = _side(side, self.dim)
        return [self.domain.nodes[i] for i, xi in self.reference.items()
                if abs(xi[axis] - value) < _ON_SIDE_TOL]


def structured_mesh(element, divisions, *, size=None, origin=None,
                    mapping: Callable | None = None, domain: Domain | None = None,
                    **element_kwargs) -> StructuredMesh:
    """Malla estructurada de ``divisions = (nx, ny[, nz])`` celdas.

    Parameters
    ----------
    element : str o clase
        Elemento del registro que declare ``REFERENCE_CELL`` y
        ``REFERENCE_NODES``.
    size, origin :
        Caja ``origin + [0, size]`` antes de la transformación; por omisión,
        el bloque unidad.
    mapping : callable, opcional
        ``x = mapping(X)`` de la caja a la geometría física; debe conservar la
        orientación (jacobiano positivo).
    domain : Domain, opcional
        Dominio donde añadir la malla; por omisión uno nuevo. Los nodos y
        elementos se numeran a continuación de los existentes.
    **element_kwargs :
        Argumentos del constructor del elemento (``material``,
        ``thickness``, ``quadrature``...).
    """
    cls = ElementRegistry.get(element) if isinstance(element, str) else element
    cell = cls.REFERENCE_CELL
    if cls.REFERENCE_NODES is None or cell not in _CELL_VERTICES:
        raise ValueError(
            f"{cls.__name__}: el generador necesita REFERENCE_CELL en "
            f"{sorted(_CELL_VERTICES)} y REFERENCE_NODES (ADR 0022).")
    dim = _DIM[cell]
    n = tuple(int(d) for d in divisions)
    if len(n) != dim or min(n) < 1:
        raise ValueError(f"divisions={divisions}: se esperan {dim} enteros positivos.")
    ref_nodes = np.asarray(cls.REFERENCE_NODES, dtype=np.float64)
    size = np.ones(dim) if size is None else np.asarray(size, dtype=np.float64)
    origin = np.zeros(dim) if origin is None else np.asarray(origin, dtype=np.float64)
    local_entities = cls.EDGE_NODES if dim == 2 else cls.FACE_NODES
    dom = Domain() if domain is None else domain
    mesh = StructuredMesh(dom, n)
    mesh.boundary = {s: [] for s in mesh.sides}

    ids: dict[tuple, int] = {}
    next_node = max(dom.nodes, default=0) + 1
    next_elem = max(dom.elements, default=0) + 1
    for index in np.ndindex(*n):
        for vertices in _CELL_VERTICES[cell]:
            nodes = []
            for xi in ref_nodes:
                ref = (np.array(index) + _in_cell(cell, vertices, xi)) / n
                key = tuple(np.round(ref, 10))
                if key not in ids:
                    X = origin + ref * size
                    x = X if mapping is None else np.asarray(mapping(X), dtype=np.float64)
                    dom.add_node(next_node, [float(c) for c in x[:dim]])
                    ids[key] = next_node
                    mesh.reference[next_node] = ref
                    next_node += 1
                nodes.append(dom.nodes[ids[key]])
            el = cls(next_elem, nodes, **element_kwargs)
            dom.add_element(el)
            next_elem += 1
            for k, local in enumerate(local_entities):
                refs = np.array([mesh.reference[nodes[j].id] for j in local])
                for side in mesh.sides:
                    axis, value = _side(side, dim)
                    if np.all(np.abs(refs[:, axis] - value) < _ON_SIDE_TOL):
                        mesh.boundary[side].append((el, k))
    return mesh
