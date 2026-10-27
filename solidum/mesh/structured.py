"""Mallas estructuradas de un bloque transformado (ADR 0022).

Se malla el bloque de referencia ``[0, 1]^d`` en ``divisions`` celdas por
dirección, se colocan en cada celda los nodos del elemento a partir de sus
coordenadas naturales (``REFERENCE_NODES``), se funden los compartidos y se
lleva cada nodo a su posición física::

    X = origin + ξ·size,     x = mapping(X)   (identidad si no se da).

Con una ``mapping`` curva, los nodos intermedios de un elemento cuadrático
quedan sobre la curva: los bordes son curvos de verdad, también en el
contorno. Con una ``mapping`` suave, al refinar los elementos tienden a
paralelogramos (mallas asintóticamente afines).

Celdas: ``"quad"`` (una por celda) y ``"tri"`` (dos por celda, separadas por
la diagonal de ``(0, 0)`` a ``(1, 1)``). Los lados del bloque se llaman
``xmin``, ``xmax``, ``ymin`` e ``ymax``; :attr:`StructuredMesh.boundary_edges`
da los bordes de elemento sobre cada lado y :meth:`StructuredMesh.nodes_on`
sus nodos.

Para geometrías reales, Gmsh (``solidum.utils.gmsh_parser``); esto es para
pruebas, estudios de convergencia y modelos de geometría simple.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.registry import ElementRegistry

SIDES_2D = ("xmin", "xmax", "ymin", "ymax")

# Vértices de las celdas por elemento en la celda unidad, en el orden de los
# vértices de la celda de referencia del elemento.
_CELL_VERTICES = {
    "quad": (((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)),),
    "tri": (((0.0, 0.0), (1.0, 0.0), (1.0, 1.0)),
            ((0.0, 0.0), (1.0, 1.0), (0.0, 1.0))),
}

# Tolerancia para decidir que un nodo está sobre un lado del bloque, en las
# coordenadas de referencia (que son fracciones exactas de 1/divisions).
_ON_SIDE_TOL = 1e-12


def _in_cell(cell: str, vertices, xi: np.ndarray) -> np.ndarray:
    """Punto de la celda unidad que corresponde a las coordenadas naturales
    ``xi`` del elemento (cuadrado ``[-1, 1]²`` o coordenadas de área)."""
    v = np.asarray(vertices)
    if cell == "quad":
        return v[0] + (xi + 1.0) / 2.0
    return v[0] + xi[0] * (v[1] - v[0]) + xi[1] * (v[2] - v[0])


@dataclass
class StructuredMesh:
    """Malla generada por :func:`structured_mesh`.

    Attributes
    ----------
    domain : Domain
        Nodos y elementos (sin numeración de ecuaciones: se llama a
        ``generate_equation_numbers`` tras imponer las condiciones de
        contorno).
    boundary_edges : dict[str, list[tuple[Element, int]]]
        Por lado del bloque, los pares (elemento, índice de borde en
        ``EDGE_NODES``) que lo forman; listos para ``compute_edge_traction``.
    """
    domain: Domain
    divisions: tuple
    boundary_edges: dict = field(default_factory=dict)
    reference: dict = field(default_factory=dict, repr=False)

    def nodes_on(self, side: str) -> list[Node]:
        """Nodos sobre el lado ``side`` del bloque, en orden de creación."""
        axis, value = _side(side)
        return [self.domain.nodes[i] for i, xi in self.reference.items()
                if abs(xi[axis] - value) < _ON_SIDE_TOL]


def _side(side: str):
    if side not in SIDES_2D:
        raise ValueError(f"Lado {side!r} desconocido; lados: {', '.join(SIDES_2D)}.")
    return "xy".index(side[0]), (0.0 if side.endswith("min") else 1.0)


def structured_mesh(element, divisions, *, size=(1.0, 1.0), origin=(0.0, 0.0),
                    mapping: Callable | None = None, domain: Domain | None = None,
                    **element_kwargs) -> StructuredMesh:
    """Malla estructurada de ``divisions = (nx, ny)`` celdas.

    Parameters
    ----------
    element : str o clase
        Elemento del registro que declare ``REFERENCE_CELL`` y
        ``REFERENCE_NODES``.
    size, origin :
        Caja ``origin + [0, size]`` antes de la transformación.
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
    n = tuple(int(d) for d in divisions)
    if len(n) != 2 or min(n) < 1:
        raise ValueError(f"divisions={divisions}: se esperan dos enteros positivos.")
    ref_nodes = np.asarray(cls.REFERENCE_NODES, dtype=np.float64)
    size = np.asarray(size, dtype=np.float64)
    origin = np.asarray(origin, dtype=np.float64)
    dom = Domain() if domain is None else domain
    mesh = StructuredMesh(dom, n, {s: [] for s in SIDES_2D})

    ids: dict[tuple, int] = {}
    next_node = max(dom.nodes, default=0) + 1
    next_elem = max(dom.elements, default=0) + 1
    for i, j in np.ndindex(*n):
        for vertices in _CELL_VERTICES[cell]:
            nodes = []
            for xi in ref_nodes:
                ref = (np.array([i, j]) + _in_cell(cell, vertices, xi)) / n
                key = tuple(np.round(ref, 10))
                if key not in ids:
                    X = origin + ref * size
                    x = X if mapping is None else np.asarray(mapping(X), dtype=np.float64)
                    dom.add_node(next_node, [float(x[0]), float(x[1])])
                    ids[key] = next_node
                    mesh.reference[next_node] = ref
                    next_node += 1
                nodes.append(dom.nodes[ids[key]])
            el = cls(next_elem, nodes, **element_kwargs)
            dom.add_element(el)
            next_elem += 1
            for edge, local in enumerate(cls.EDGE_NODES):
                refs = np.array([mesh.reference[nodes[k].id] for k in local])
                for side in SIDES_2D:
                    axis, value = _side(side)
                    if np.all(np.abs(refs[:, axis] - value) < _ON_SIDE_TOL):
                        mesh.boundary_edges[side].append((el, edge))
    return mesh
