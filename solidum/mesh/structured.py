"""Mallas estructuradas de un bloque transformado (ADR 0022 y 0037).

Se malla el bloque de referencia ``[0, 1]^d`` en ``divisions`` celdas por
dirección, se colocan en cada celda los nodos del elemento a partir de sus
coordenadas naturales (``REFERENCE_NODES``), se funden los compartidos y se
lleva cada nodo a su posición física::

    X = origin + ξ·size,     x = mapping(X)   (identidad si no se da).

Con ``lines`` (ADR 0037), las líneas de la rejilla de cada eje se dan
explícitamente (p. ej. con :func:`axis_lines`, por tramos con gradación): la
celda ``i`` va de ``lines[i]`` a ``lines[i+1]`` y, dentro de ella, ``X`` es
lineal en ``ξ`` (los nodos intermedios de un elemento cuadrático quedan en el
centro de sus bordes).

Con una ``mapping`` curva, los nodos intermedios de un elemento cuadrático
quedan sobre la curva: los bordes y las caras son curvos de verdad, también
en el contorno. Con una ``mapping`` suave, al refinar los elementos tienden a
paralelogramos (mallas asintóticamente afines).

Celdas (``REFERENCE_CELL`` del elemento):

- ``"line"``: elementos de línea (armaduras, marcos), uno por celda; la
  ``mapping`` lleva ``[0, 1]`` a la línea física en 2D o 3D (sin ella, el
  eje ``x`` del plano);
- ``"quad"``: una por celda; ``"tri"``: dos por celda, separadas por la
  diagonal de ``(0, 0)`` a ``(1, 1)`` (``diagonal="right"``, por omisión), la
  de ``(1, 0)`` a ``(0, 1)`` (``"left"``), una u otra en tablero de ajedrez
  (``"alternate"``), o cuatro por celda alrededor de un nodo central
  (``"crossed"``);
- ``"hex"``: una por celda; ``"tet"``: seis por celda (subdivisión de Kuhn
  alrededor de la diagonal de ``(0, 0, 0)`` a ``(1, 1, 1)``; conforme entre
  celdas vecinas porque todas parten sus caras por la misma diagonal).

Los lados del bloque se llaman ``xmin``, ``xmax``, ``ymin``, ``ymax`` (y
``zmin``, ``zmax`` en 3D); una línea sólo tiene sus extremos ``xmin`` y
``xmax``. :attr:`StructuredMesh.boundary` da los bordes (2D, ``EDGE_NODES``)
o caras (3D, ``FACE_NODES``) de elemento sobre cada lado, o el nodo extremo
(índice local 0 ó 1) en una línea, y :meth:`StructuredMesh.nodes_on`, sus
nodos.

``exclude`` quita las celdas cuyo centro (en la geometría física) cumple una
condición: entallas, ranuras y huecos alineados con la rejilla. Con
``merge_tol``, un nodo que cae a esa distancia de un nodo que ya estaba en el
``domain`` se funde con él: varios bloques forman un dominio conforme si
comparten las líneas de la rejilla en sus lados comunes.

Para geometrías reales, Gmsh (``solidum.utils.gmsh_parser``); esto es para
pruebas, estudios de convergencia y modelos de geometría simple o por bloques.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import permutations
from typing import Callable

import numpy as np

from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.registry import ElementRegistry

SIDES_1D = ("xmin", "xmax")
SIDES_2D = SIDES_1D + ("ymin", "ymax")
SIDES_3D = SIDES_2D + ("zmin", "zmax")
_SIDES = {1: SIDES_1D, 2: SIDES_2D, 3: SIDES_3D}


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
    "line": (((0.0,),),),
    "quad": (((0.0, 0.0),),),
    "tri": (((0.0, 0.0), (1.0, 0.0), (1.0, 1.0)),
            ((0.0, 0.0), (1.0, 1.0), (0.0, 1.0))),
    "hex": (((0.0, 0.0, 0.0),),),
    "tet": _kuhn_tets(),
}
_DIM = {"line": 1, "quad": 2, "tri": 2, "hex": 3, "tet": 3}

# Partición de la celda en triángulos según ``diagonal``; vértices en sentido
# antihorario (jacobiano positivo).
_TRI_RIGHT = _CELL_VERTICES["tri"]
_TRI_LEFT = (((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
             ((1.0, 0.0), (1.0, 1.0), (0.0, 1.0)))
_TRI_CROSSED = (((0.0, 0.0), (1.0, 0.0), (0.5, 0.5)),
                ((1.0, 0.0), (1.0, 1.0), (0.5, 0.5)),
                ((1.0, 1.0), (0.0, 1.0), (0.5, 0.5)),
                ((0.0, 1.0), (0.0, 0.0), (0.5, 0.5)))
DIAGONALS = ("right", "left", "alternate", "crossed")


def _cell_partition(cell: str, diagonal: str, index) -> tuple:
    if cell != "tri":
        return _CELL_VERTICES[cell]
    if diagonal == "alternate":
        return _TRI_RIGHT if sum(index) % 2 == 0 else _TRI_LEFT
    return {"right": _TRI_RIGHT, "left": _TRI_LEFT, "crossed": _TRI_CROSSED}[diagonal]


def axis_lines(segments, origin: float = 0.0) -> np.ndarray:
    """Líneas de la rejilla de un eje, por tramos consecutivos desde ``origin``.

    Cada tramo es ``(length, n)`` o ``(length, n, ratio)``: ``n`` celdas que
    cubren ``length``, en progresión geométrica cuyo último tamaño es
    ``ratio`` veces el primero (``ratio = 1``, uniforme; ``ratio < 1``, se
    afinan hacia el final del tramo; ``ratio > 1``, hacia el principio).
    """
    out = [float(origin)]
    for k, seg in enumerate(segments):
        if len(seg) not in (2, 3):
            raise ValueError(f"Tramo {k}: se espera (length, n) o (length, n, ratio); "
                             f"se recibió {seg!r}.")
        length, n = float(seg[0]), seg[1]
        ratio = float(seg[2]) if len(seg) == 3 else 1.0
        if not (length > 0.0) or int(n) != n or n < 1 or not (ratio > 0.0):
            raise ValueError(f"Tramo {k}: length > 0, n entero ≥ 1 y ratio > 0; se recibió "
                             f"{seg!r}.")
        n = int(n)
        if n == 1 or ratio == 1.0:
            sizes = np.full(n, 1.0)
        else:
            sizes = ratio ** (np.arange(n) / (n - 1))
        steps = length * np.cumsum(sizes) / sizes.sum()
        steps[-1] = length
        start = out[-1]
        out.extend(start + steps)
    return np.array(out)


# Entidades de contorno de una línea: sus dos nodos extremos.
_LINE_ENDS = ((0,), (1,))

# Tolerancia para decidir que un nodo está sobre un lado del bloque, en las
# coordenadas de referencia (que son fracciones exactas de 1/divisions).
_ON_SIDE_TOL = 1e-12


def _in_cell(cell: str, vertices, xi: np.ndarray) -> np.ndarray:
    """Punto de la celda unidad que corresponde a las coordenadas naturales
    ``xi`` del elemento (``[-1, 1]^d`` o coordenadas de volumen/área)."""
    v = np.asarray(vertices)
    if cell in ("line", "quad", "hex"):
        return v[0] + (xi + 1.0) / 2.0
    return v[0] + (v[1:] - v[0]).T @ xi


def _side(side: str, dim: int):
    sides = _SIDES[dim]
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
    elements : list[Element]
        Los elementos de este bloque, en orden de creación.
    excluded, merged : int
        Celdas quitadas por ``exclude`` y nodos fundidos con los del
        ``domain`` por ``merge_tol``.
    """
    domain: Domain
    divisions: tuple
    boundary: dict = field(default_factory=dict)
    reference: dict = field(default_factory=dict, repr=False)
    elements: list = field(default_factory=list, repr=False)
    excluded: int = 0
    merged: int = 0

    @property
    def dim(self) -> int:
        return len(self.divisions)

    @property
    def sides(self) -> tuple:
        return _SIDES[self.dim]

    def nodes_on(self, side: str) -> list[Node]:
        """Nodos sobre el lado ``side`` del bloque, en orden de creación."""
        axis, value = _side(side, self.dim)
        return [self.domain.nodes[i] for i, xi in self.reference.items()
                if abs(xi[axis] - value) < _ON_SIDE_TOL]


def _grid(n, lines, size, origin, dim):
    """Divisiones y transformación ``ξ → X`` (referencia → caja) por eje."""
    if lines is None:
        if n is None:
            raise ValueError("Falta 'divisions' (o 'lines').")
        n = tuple(int(d) for d in n)
        if len(n) != dim or min(n) < 1:
            raise ValueError(f"divisions={n}: se esperan {dim} enteros positivos.")
        size = np.ones(dim) if size is None else np.asarray(size, dtype=np.float64)
        origin = np.zeros(dim) if origin is None else np.asarray(origin, dtype=np.float64)
        return n, lambda ref: origin + ref * size
    if size is not None or origin is not None:
        raise ValueError("Con 'lines' no se dan 'size' ni 'origin': las líneas ya fijan la caja.")
    if len(lines) != dim:
        raise ValueError(f"'lines' necesita {dim} ejes; se recibieron {len(lines)}.")
    axes = []
    for a, ln in enumerate(lines):
        ln = np.asarray(ln, dtype=np.float64)
        if ln.ndim != 1 or ln.size < 2 or np.any(np.diff(ln) <= 0.0):
            raise ValueError(f"lines[{a}]: se esperan al menos dos coordenadas estrictamente "
                             "crecientes.")
        axes.append(ln)
    n_lines = tuple(ln.size - 1 for ln in axes)
    if n is not None and tuple(int(d) for d in n) != n_lines:
        raise ValueError(f"divisions={tuple(n)} no coincide con las {n_lines} celdas de 'lines'.")
    knots = [np.arange(k + 1, dtype=np.float64) for k in n_lines]
    return n_lines, lambda ref: np.array(
        [np.interp(ref[a] * n_lines[a], knots[a], axes[a]) for a in range(dim)])


def structured_mesh(element, divisions=None, *, size=None, origin=None, lines=None,
                    mapping: Callable | None = None, exclude: Callable | None = None,
                    diagonal: str | None = None, merge_tol: float | None = None,
                    domain: Domain | None = None, **element_kwargs) -> StructuredMesh:
    """Malla estructurada de ``divisions = (nx, ny[, nz])`` celdas.

    Parameters
    ----------
    element : str o clase
        Elemento del registro que declare ``REFERENCE_CELL`` y
        ``REFERENCE_NODES``.
    size, origin :
        Caja ``origin + [0, size]`` antes de la transformación; por omisión,
        el bloque unidad.
    lines : secuencia de ``d`` arreglos, opcional
        Líneas de la rejilla de cada eje (estrictamente crecientes; véase
        :func:`axis_lines`). Sustituye a ``size`` y ``origin``, y fija
        ``divisions``.
    mapping : callable, opcional
        ``x = mapping(X)`` de la caja a la geometría física; debe conservar la
        orientación (jacobiano positivo). En una línea da las 2 ó 3
        coordenadas del punto.
    exclude : callable, opcional
        ``exclude(x) -> bool`` en el centro físico de cada celda: si es
        verdadero, la celda no se malla (entallas, ranuras, huecos).
    diagonal : str, opcional
        Partición de las celdas en triángulos: ``"right"`` (por omisión),
        ``"left"``, ``"alternate"`` o ``"crossed"``. Sólo con celdas ``tri``.
    merge_tol : float, opcional
        Funde cada nodo nuevo con el nodo del ``domain`` (anterior a esta
        llamada) más cercano si está a esa distancia o menos.
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
    if diagonal is not None and cell != "tri":
        raise ValueError(f"'diagonal' sólo se aplica a triángulos; {cls.__name__} tiene "
                         f"celdas {cell!r}.")
    diagonal = diagonal or "right"
    if diagonal not in DIAGONALS:
        raise ValueError(f"diagonal={diagonal!r}: se admite uno de {', '.join(DIAGONALS)}.")
    n, to_box = _grid(divisions, lines, size, origin, dim)
    ref_nodes = np.asarray(cls.REFERENCE_NODES, dtype=np.float64)
    local_entities = {1: _LINE_ENDS, 2: getattr(cls, "EDGE_NODES", None),
                      3: getattr(cls, "FACE_NODES", None)}[dim]
    dom = Domain() if domain is None else domain
    mesh = StructuredMesh(dom, n)
    mesh.boundary = {s: [] for s in mesh.sides}

    def physical(ref):
        X = to_box(ref)
        x = X if mapping is None else np.asarray(mapping(X), dtype=np.float64)
        if dim > 1:
            return x[:dim]
        return np.array([x[0], 0.0]) if mapping is None else x

    tree, old_ids = None, []
    if merge_tol is not None and dom.nodes:
        from scipy.spatial import cKDTree
        old_ids = list(dom.nodes)
        width = max(len(dom.nodes[i].coordinates) for i in old_ids)
        pts = np.zeros((len(old_ids), width))
        for r, i in enumerate(old_ids):
            c = dom.nodes[i].coordinates
            pts[r, :len(c)] = c
        tree = cKDTree(pts)

    def existing(x):
        if tree is None:
            return None
        q = np.zeros(tree.m)
        k = min(tree.m, x.size)
        q[:k] = x[:k]
        dist, r = tree.query(q)
        return old_ids[r] if dist <= merge_tol else None

    ids: dict[tuple, int] = {}
    next_node = max(dom.nodes, default=0) + 1
    next_elem = max(dom.elements, default=0) + 1
    for index in np.ndindex(*n):
        if exclude is not None and exclude(physical((np.array(index) + 0.5) / n)):
            mesh.excluded += 1
            continue
        for vertices in _cell_partition(cell, diagonal, index):
            nodes = []
            for xi in ref_nodes:
                ref = (np.array(index) + _in_cell(cell, vertices, xi)) / n
                key = tuple(np.round(ref, 10))
                if key not in ids:
                    x = physical(ref)
                    nid = existing(x)
                    if nid is None:
                        dom.add_node(next_node, [float(c) for c in x])
                        nid = next_node
                        next_node += 1
                    else:
                        mesh.merged += 1
                    ids[key] = nid
                    mesh.reference[nid] = ref
                nodes.append(dom.nodes[ids[key]])
            el = cls(next_elem, nodes, **element_kwargs)
            dom.add_element(el)
            mesh.elements.append(el)
            next_elem += 1
            for k, local in enumerate(local_entities):
                refs = np.array([mesh.reference[nodes[j].id] for j in local])
                for side in mesh.sides:
                    axis, value = _side(side, dim)
                    if np.all(np.abs(refs[:, axis] - value) < _ON_SIDE_TOL):
                        mesh.boundary[side].append((el, k))
    if not mesh.elements:
        raise ValueError("'exclude' quita todas las celdas del bloque.")
    return mesh
