"""Malla estructurada por bloques desde el YAML (ADR 0037).

``mesh: {structured: ...}`` malla uno o varios bloques con
:func:`solidum.mesh.structured.structured_mesh` y los funde en un dominio
conforme. Un bloque declara, por eje, sus tramos (``length``, ``n`` y
``ratio`` opcional de gradación, véase :func:`axis_lines`), y puede quitar
celdas (``exclude``), elegir la partición en triángulos (``diagonal``) y
transformarse (``mapping``); el resto de sus claves es la receta del elemento
(``element``, ``material``, ``thickness``, ``quadrature`` y los parámetros del
elemento), con ``mesh_element``, ``mesh_material``... por omisión::

    mesh:
      structured:
        element: Quad4
        x: [{length: 950, n: 38, ratio: 0.1}, {length: 100, n: 40},
            {length: 950, n: 38, ratio: 10}]
        y: {length: 200, n: 40}
        exclude: "abs(x - 1000) < 1 and y < 100"
        groups:
          Carga: "near(x, 1000) and near(y, 200)"

Con varios bloques, ``blocks: [...]`` y, al nivel de ``structured``,
``groups`` y ``merge_tol``. Los nodos de bloques distintos que coinciden se
funden.

Grupos (los mismos que da una malla de Gmsh: nodos y celdas de contorno, para
``boundary_conditions_by_group``, ``tractions_by_group``, ``history``...):

- los lados de cada bloque: ``xmin``, ``xmax``, ``ymin``... con un solo
  bloque sin nombre; ``<nombre>.xmin``... en otro caso (``block1``... si el
  bloque no tiene ``name``);
- los de ``groups``: los nodos cuya posición cumple la condición y los bordes
  o caras del contorno libre con todos sus nodos en el grupo. La condición
  admite ``near(a, b[, tol])``: ``|a − b| ≤ tol``, por omisión 10⁻⁶ veces la
  diagonal de la caja que envuelve la malla.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from solidum.core.domain import Domain
from solidum.logging import get_logger
from solidum.mesh.structured import DIAGONALS, axis_lines, structured_mesh
from solidum.registry import ElementRegistry
from solidum.utils.expressions import (POSITION, ExpressionError, compile_expression,
                                       compile_predicate)

_log = get_logger("mesh")

#: Claves de ``structured`` con varios bloques.
STRUCTURED_KEYS = ("blocks", "groups", "merge_tol")
#: Claves de un bloque que no son de la receta del elemento.
BLOCK_KEYS = ("name", "x", "y", "z", "origin", "mapping", "exclude", "diagonal")
_AXES = ("x", "y", "z")
_SEGMENT_KEYS = ("length", "n", "ratio")
_DIM = {"line": 1, "quad": 2, "tri": 2, "hex": 3, "tet": 3}
# Tolerancia relativa de fusión de nodos y de ``near`` (por la escala).
_REL_TOL = 1e-6


def _near_validation(a, b, tol=None):
    return abs(a - b) <= (1e-9 if tol is None else tol)


def _split(structured: dict) -> tuple[list, dict, object]:
    """(bloques, grupos, merge_tol) de la forma con ``blocks`` o de un bloque."""
    if "blocks" in structured:
        return (structured.get("blocks"), structured.get("groups") or {},
                structured.get("merge_tol"))
    block = {k: v for k, v in structured.items() if k not in ("groups", "merge_tol")}
    return [block], structured.get("groups") or {}, structured.get("merge_tol")


def _segments(spec, ctx: str, errors: list) -> list | None:
    items = [spec] if isinstance(spec, dict) else spec
    if not isinstance(items, list) or not items:
        errors.append(f"{ctx}: se esperaba un tramo {{length, n[, ratio]}} o una lista de "
                      f"tramos; se recibió {spec!r}.")
        return None
    out = []
    for k, seg in enumerate(items):
        where = f"{ctx}[{k}]" if isinstance(spec, list) else ctx
        if not isinstance(seg, dict):
            errors.append(f"{where}: se esperaba {{length, n[, ratio]}}; se recibió {seg!r}.")
            continue
        extra = sorted(set(seg) - set(_SEGMENT_KEYS))
        if extra:
            errors.append(f"{where}: claves no admitidas {extra}; se admiten "
                          f"{list(_SEGMENT_KEYS)}.")
        length, n, ratio = seg.get("length"), seg.get("n"), seg.get("ratio", 1.0)
        ok = True
        if not _is_number(length) or not length > 0:
            errors.append(f"{where}: 'length' debe ser un número > 0; se recibió {length!r}.")
            ok = False
        if isinstance(n, bool) or not isinstance(n, int) or n < 1:
            errors.append(f"{where}: 'n' debe ser un entero ≥ 1; se recibió {n!r}.")
            ok = False
        if not _is_number(ratio) or not ratio > 0:
            errors.append(f"{where}: 'ratio' debe ser un número > 0; se recibió {ratio!r}.")
            ok = False
        if ok:
            out.append((float(length), n, float(ratio)))
    return out


def _is_number(v) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def validate_structured_block(structured) -> list[str]:
    """Errores de forma de ``mesh.structured`` (todos juntos). La dimensión de
    cada bloque se comprueba al mallar, cuando se conoce su elemento."""
    errors: list[str] = []
    if not isinstance(structured, dict):
        return [f"mesh.structured: se esperaba un diccionario; se recibió {structured!r}."]
    if "blocks" in structured:
        extra = sorted(set(structured) - set(STRUCTURED_KEYS))
        if extra:
            errors.append(f"mesh.structured: con 'blocks', las claves {extra} van dentro de cada "
                          f"bloque; al nivel de 'structured' sólo {list(STRUCTURED_KEYS)}.")
    blocks, groups, merge_tol = _split(structured)
    if not isinstance(blocks, list) or not blocks:
        errors.append("mesh.structured.blocks: se esperaba una lista de bloques.")
        blocks = []
    names = []
    for k, block in enumerate(blocks):
        ctx = f"mesh.structured.blocks[{k}]" if "blocks" in structured else "mesh.structured"
        if not isinstance(block, dict):
            errors.append(f"{ctx}: se esperaba un diccionario.")
            continue
        if "name" in block:
            names.append(block["name"])
            if not isinstance(block["name"], str) or not block["name"]:
                errors.append(f"{ctx}.name: se esperaba un nombre.")
        if "x" not in block:
            errors.append(f"{ctx}: falta 'x' (los tramos del eje x).")
        for axis in _AXES:
            if axis in block:
                _segments(block[axis], f"{ctx}.{axis}", errors)
        origin = block.get("origin")
        if origin is not None and (not isinstance(origin, list)
                                   or not all(_is_number(v) for v in origin)):
            errors.append(f"{ctx}.origin: se esperaba una lista de números.")
        mapping = block.get("mapping")
        if mapping is not None:
            if not isinstance(mapping, list) or not mapping:
                errors.append(f"{ctx}.mapping: se esperaba una lista de expresiones de x, y, z.")
            else:
                for j, comp in enumerate(mapping):
                    if not _is_number(comp):
                        try:
                            compile_expression(comp, POSITION)
                        except ExpressionError as exc:
                            errors.append(f"{ctx}.mapping[{j}]: {exc}")
        if "exclude" in block:
            try:
                compile_predicate(block["exclude"], POSITION, {"near": _near_validation})
            except ExpressionError as exc:
                errors.append(f"{ctx}.exclude: {exc}")
        diagonal = block.get("diagonal")
        if diagonal is not None and diagonal not in DIAGONALS:
            errors.append(f"{ctx}.diagonal: {diagonal!r} no admitido; uno de {list(DIAGONALS)}.")
    if len(set(map(str, names))) != len(names):
        errors.append("mesh.structured: nombres de bloque repetidos.")
    if not isinstance(groups, dict):
        errors.append("mesh.structured.groups: se esperaba {nombre: condición}.")
    else:
        for name, cond in groups.items():
            try:
                compile_predicate(cond, POSITION, {"near": _near_validation})
            except ExpressionError as exc:
                errors.append(f"mesh.structured.groups[{name!r}]: {exc}")
    if merge_tol is not None and (not _is_number(merge_tol) or not merge_tol > 0):
        errors.append(f"mesh.structured.merge_tol: se esperaba un número > 0; se recibió "
                      f"{merge_tol!r}.")
    return errors


def _padded(c) -> tuple:
    return tuple(float(c[k]) if k < len(c) else 0.0 for k in range(3))


def _mapping(exprs) -> Callable | None:
    if exprs is None:
        return None
    parts = [(lambda *_, c=float(e): c) if _is_number(e) else compile_expression(e, POSITION)
             for e in exprs]
    return lambda X: np.array([f(*_padded(X)) for f in parts])


def build_structured_mesh(structured: dict, recipe: Callable, errors: list) -> Domain | None:
    """Dominio de ``mesh.structured``. ``recipe(props, ctx)`` devuelve
    ``(nombre del elemento, kwargs del constructor)`` de la receta de un
    bloque, o ``None`` tras anotar el error en ``errors``. Los errores se
    acumulan en ``errors``; si hay alguno, devuelve ``None``."""
    blocks, groups, merge_tol = _split(structured)
    single = len(blocks) == 1 and "name" not in blocks[0]
    many = "blocks" in structured

    plans = []
    for k, block in enumerate(blocks):
        ctx = f"mesh.structured.blocks[{k}]" if many else "mesh.structured"
        props = {key: v for key, v in block.items() if key not in BLOCK_KEYS}
        made = recipe(props, ctx)
        if made is None:
            continue
        element, kwargs = made
        cell = getattr(ElementRegistry.get(element), "REFERENCE_CELL", None)
        if cell not in _DIM:
            errors.append(f"{ctx}: el elemento '{element}' no se puede mallar por bloques "
                          f"(celda {cell!r}).")
            continue
        dim = _DIM[cell]
        axes = _AXES[:dim]
        missing = [a for a in axes if a not in block]
        extra = [a for a in _AXES[dim:] if a in block]
        if missing or extra:
            errors.append(f"{ctx}: el elemento '{element}' es {dim}D: se esperan los ejes "
                          f"{list(axes)}" + (f"; faltan {missing}" if missing else "")
                          + (f"; sobran {extra}" if extra else "") + ".")
            continue
        origin = block.get("origin") or [0.0] * dim
        if len(origin) != dim:
            errors.append(f"{ctx}.origin: se esperaban {dim} coordenadas.")
            continue
        lines = []
        for a, axis in enumerate(axes):
            segs = _segments(block[axis], f"{ctx}.{axis}", errors)
            if segs:
                lines.append(axis_lines(segs, origin[a]))
        if len(lines) != dim:
            continue
        name = block.get("name") or ("" if single else f"block{k + 1}")
        plans.append((ctx, name, element, kwargs, lines, block))
    if errors:
        return None

    if merge_tol is None:
        h_min = min(float(np.diff(ln).min()) for plan in plans for ln in plan[4])
        merge_tol = _REL_TOL * h_min
    dom = Domain()
    side_groups: dict[str, tuple[list, list]] = {}
    for ctx, name, element, kwargs, lines, block in plans:
        exclude = None
        if "exclude" in block:
            exclude = _block_predicate(block["exclude"], lines, block.get("mapping"))
        try:
            mesh = structured_mesh(element, lines=lines, mapping=_mapping(block.get("mapping")),
                                   exclude=exclude, diagonal=block.get("diagonal"),
                                   merge_tol=merge_tol, domain=dom, **kwargs)
        except ValueError as exc:
            errors.append(f"{ctx}: {exc}")
            continue
        entities = {1: ((0,), (1,)), 2: getattr(mesh.elements[0], "EDGE_NODES", None),
                    3: getattr(mesh.elements[0], "FACE_NODES", None)}[mesh.dim]
        for side in mesh.sides:
            gname = f"{name}.{side}" if name else side
            cells = [tuple(el.nodes[j].id for j in entities[k]) for el, k in mesh.boundary[side]]
            side_groups[gname] = (sorted({nd.id for nd in mesh.nodes_on(side)}), cells)
        _log.info(f"  -> Bloque {name or '(único)'}: {len(mesh.elements)} elementos {element}, "
                  f"{mesh.excluded} celdas excluidas, {mesh.merged} nodos fundidos con bloques "
                  "anteriores.")
    if errors:
        return None

    dom.physical_groups = {g: ids for g, (ids, _) in side_groups.items()}
    dom.physical_group_cells = {g: cells for g, (_, cells) in side_groups.items()}
    if groups:
        _predicate_groups(dom, groups, errors)
    _log.info(f"Malla estructurada: {len(dom.nodes)} nodos, {len(dom.elements)} elementos, "
              f"grupos {sorted(dom.physical_groups)}.")
    return None if errors else dom


def _block_predicate(text: str, lines, mapping_exprs) -> Callable:
    """``exclude`` de un bloque con ``near`` a la escala del bloque."""
    lo = [ln[0] for ln in lines]
    hi = [ln[-1] for ln in lines]
    scale = float(np.linalg.norm(np.subtract(hi, lo)))
    near = _make_near(_REL_TOL * scale)
    pred = compile_predicate(text, POSITION, {"near": near})
    return lambda x: pred(*_padded(x))


def _make_near(default_tol: float) -> Callable:
    def near(a, b, tol=None):
        return abs(a - b) <= (default_tol if tol is None else tol)
    return near


def _predicate_groups(dom: Domain, groups: dict, errors: list) -> None:
    coords = {i: _padded(nd.coordinates) for i, nd in dom.nodes.items()}
    pts = np.array(list(coords.values()))
    near = _make_near(_REL_TOL * float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0))))
    # Bordes o caras del contorno libre: los de un solo elemento.
    count: dict[frozenset, list] = {}
    for el in dom.elements.values():
        attr = "FACE_NODES" if len(el.nodes[0].coordinates) == 3 else "EDGE_NODES"
        for local in getattr(el, attr, None) or ():
            ids = tuple(el.nodes[j].id for j in local)
            count.setdefault(frozenset(ids), []).append(ids)
    free = [v[0] for v in count.values() if len(v) == 1]
    for name, cond in groups.items():
        if name in dom.physical_groups:
            errors.append(f"mesh.structured.groups[{name!r}]: el nombre ya es el de un lado de "
                          "un bloque.")
            continue
        pred = compile_predicate(cond, POSITION, {"near": near})
        try:
            inside = {i for i, c in coords.items() if pred(*c)}
        except (ArithmeticError, ValueError, TypeError) as exc:
            errors.append(f"mesh.structured.groups[{name!r}]: la condición {cond!r} no se puede "
                          f"evaluar: {exc}.")
            continue
        if not inside:
            errors.append(f"mesh.structured.groups[{name!r}]: ningún nodo cumple {cond!r}.")
            continue
        dom.physical_groups[name] = sorted(inside)
        dom.physical_group_cells[name] = [c for c in free if set(c) <= inside]
