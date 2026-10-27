"""Cargas distribuidas del YAML: tracción, presión, fuente y flujo (ADR 0023).

Resuelve, una vez construida la malla, las entradas del YAML en cargas de
elemento listas para ensamblar, y acumula todos los errores para mostrarlos
juntos:

- de contorno, sobre un borde (2D) o una cara (3D): la tracción mecánica
  ``t`` o la presión ``p`` (``t = −p·n``, ``p > 0`` comprime) y el flujo
  térmico ``q`` (``q > 0`` saliente). Se declaran por elemento
  (``{element, edge | face, ...}``) o por grupo físico de una malla gmsh, cuyas
  celdas de contorno son los bordes exactos que se cargan;
- de volumen: la fuente térmica ``Q`` sobre todos los elementos o una lista.

Los valores son números o expresiones (:mod:`solidum.utils.expressions`) de
``x``, ``y``, ``z`` y, en las de contorno, de la normal exterior ``nx``,
``ny``, ``nz``.
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np

from solidum.utils.expressions import ExpressionError, load_field


class BoundaryLoad(NamedTuple):
    """Carga sobre la entidad ``index`` (borde o cara) de ``element``: se
    aplica con ``getattr(element, method)(index, value)``."""
    element: object
    method: str
    index: int
    value: object


# Por familia: métodos del elemento por tipo de entidad, claves del valor y
# nombre de la carga en los mensajes.
_FAMILIES = {
    "mechanical": dict(methods={"edge": "compute_edge_traction", "face": "compute_face_traction"},
                       keys=("t", "pressure"), what="tracciones"),
    "thermal": dict(methods={"edge": "compute_edge_flux", "face": "compute_face_flux"},
                    keys=("q", "value"), what="flujos de frontera"),
}
_ENTITIES = {"edge": "EDGE_NODES", "face": "FACE_NODES"}


def _dimension(element) -> int:
    return 3 if "uz" in getattr(element, "DOF_NAMES", ()) else 2


def _value(spec: dict, family: str, element, context: str):
    """Campo de la carga de contorno de ``spec`` para ``element``."""
    if family == "thermal":
        q = spec["q"] if "q" in spec else spec["value"]
        return load_field(q, 1, normal=True, context=context)
    if "t" in spec:
        return load_field(spec["t"], _dimension(element), normal=True, context=context)
    p = load_field(spec["pressure"], 1, normal=True, context=context)
    if callable(p):
        return lambda x, n: -p(x, n) * np.asarray(n, dtype=np.float64)
    return lambda x, n: -p * np.asarray(n, dtype=np.float64)


def _value_keys_ok(spec: dict, family: str, context: str, errors: list) -> bool:
    keys = _FAMILIES[family]["keys"]
    present = [k for k in keys if k in spec]
    if family == "mechanical" and len(present) != 1:
        errors.append(f"{context}: declare exactamente uno de 't' (tracción en ejes globales) "
                      "o 'pressure' (presión normal, positiva si comprime).")
        return False
    if family == "thermal" and not present:
        errors.append(f"{context}: falta el flujo 'q'.")
        return False
    return True


def _entity_index(domain, family: str) -> dict:
    """``frozenset(ids de nodo)`` → [(elemento, tipo, índice local)] de los
    bordes y caras de los elementos que admiten la carga de ``family``."""
    methods = _FAMILIES[family]["methods"]
    index: dict = {}
    for el in domain.elements.values():
        for kind, attr in _ENTITIES.items():
            if not hasattr(el, methods[kind]):
                continue
            for k, local in enumerate(getattr(el, attr, None) or ()):
                key = frozenset(el.nodes[j].id for j in local)
                index.setdefault(key, []).append((el, kind, k))
    return index


def _group_entities(domain, group, family: str, context: str, errors: list) -> list:
    """Bordes o caras de contorno que forman el grupo físico ``group``: las
    celdas del grupo que coinciden, por sus nodos, con la entidad de un único
    elemento."""
    cells = getattr(domain, "physical_group_cells", None)
    if not cells:
        errors.append(f"{context}: el modelo no tiene grupos físicos (sólo los aporta una "
                      "malla gmsh vía 'mesh').")
        return []
    if group not in cells:
        errors.append(f"{context}: grupo físico {group!r} inexistente. Disponibles: "
                      f"{sorted(cells)}.")
        return []
    index = _entity_index(domain, family)
    found, missing, interior = [], 0, 0
    for cell in cells[group]:
        hits = index.get(frozenset(cell), [])
        if len(hits) == 1:
            found.append(hits[0])
        elif hits:
            interior += 1
        else:
            missing += 1
    if interior:
        errors.append(f"{context}: {interior} celda(s) del grupo {group!r} son bordes o caras "
                      "compartidos por dos elementos, no del contorno.")
    if missing:
        errors.append(f"{context}: {missing} celda(s) del grupo {group!r} no son un borde ni una "
                      f"cara de un elemento que admita {_FAMILIES[family]['what']} (¿es un grupo "
                      "de superficie o de volumen?).")
    return found


def _group_specs(raw, context: str, errors: list) -> list:
    """Forma dict ``{grupo: {...}}`` o lista ``[{group_name: ..., ...}]``."""
    if not raw:
        return []
    if isinstance(raw, dict):
        out = []
        for name, spec in raw.items():
            if not isinstance(spec, dict):
                errors.append(f"{context}[{name!r}]: se esperaba un diccionario con el valor de "
                              "la carga.")
                continue
            out.append({"group_name": name, **spec})
        return out
    if isinstance(raw, list) and all(isinstance(s, dict) for s in raw):
        return raw
    errors.append(f"{context}: se esperaba un diccionario {{grupo: carga}} o una lista de "
                  "diccionarios con 'group_name'.")
    return []


def resolve_boundary_loads(domain, entries, by_group, family: str, names: tuple) -> tuple:
    """Cargas de contorno de ``family`` ("mechanical" o "thermal") a partir de
    las entradas por elemento (``names[0]``) y por grupo (``names[1]``).
    Devuelve ``(cargas, errores)``."""
    errors: list = []
    loads: list = []
    fam = _FAMILIES[family]
    value_keys = set(fam["keys"])
    if entries and not isinstance(entries, list):
        entries = [entries] if isinstance(entries, dict) else None
        if entries is None:
            return [], [f"{names[0]}: se esperaba una lista de entradas."]
    for i, spec in enumerate(entries or []):
        ctx = f"{names[0]}[{i}]"
        if not isinstance(spec, dict):
            errors.append(f"{ctx}: cada entrada debe ser un diccionario.")
            continue
        unknown = set(spec) - {"element", "edge", "face"} - value_keys
        if unknown:
            errors.append(f"{ctx}: claves no admitidas {sorted(unknown)}; se admiten 'element', "
                          f"'edge' o 'face' y {sorted(value_keys)}.")
            continue
        el = domain.elements.get(spec.get("element"))
        if el is None:
            errors.append(f"{ctx}: elemento {spec.get('element')!r} inexistente.")
            continue
        kinds = [k for k in ("edge", "face") if k in spec]
        if len(kinds) != 1:
            errors.append(f"{ctx}: declare 'edge' (elementos 2D) o 'face' (elementos 3D).")
            continue
        kind = kinds[0]
        method = fam["methods"][kind]
        if not hasattr(el, method):
            errors.append(f"{ctx}: el elemento {el.id} ({type(el).__name__}) no admite "
                          f"{fam['what']} por '{kind}'.")
            continue
        n_ent = len(getattr(el, _ENTITIES[kind], None) or ())
        k = spec[kind]
        if not isinstance(k, int) or isinstance(k, bool) or not 0 <= k < n_ent:
            errors.append(f"{ctx}: {kind}={k!r} fuera de rango para {type(el).__name__} "
                          f"(0..{n_ent - 1}).")
            continue
        if not _value_keys_ok(spec, family, ctx, errors):
            continue
        try:
            loads.append(BoundaryLoad(el, method, k, _value(spec, family, el, ctx)))
        except ExpressionError as exc:
            errors.append(str(exc))
    for spec in _group_specs(by_group, names[1], errors):
        group = spec.get("group_name")
        ctx = f"{names[1]}[{group!r}]"
        unknown = set(spec) - {"group_name"} - value_keys
        if unknown:
            errors.append(f"{ctx}: claves no admitidas {sorted(unknown)}; se admiten "
                          f"{sorted(value_keys)}.")
            continue
        if not _value_keys_ok(spec, family, ctx, errors):
            continue
        entities = _group_entities(domain, group, family, ctx, errors)
        try:
            for el, kind, k in entities:
                loads.append(BoundaryLoad(el, fam["methods"][kind], k,
                                          _value(spec, family, el, ctx)))
        except ExpressionError as exc:
            errors.append(str(exc))
    return loads, errors


def resolve_body_sources(domain, sources) -> tuple:
    """Fuentes térmicas ``[{Q | value, elements?}]`` → ``[(elementos, Q)]``,
    con ``Q`` constante o campo ``Q(x)``. Devuelve ``(fuentes, errores)``."""
    errors: list = []
    out: list = []
    for i, src in enumerate(sources or []):
        ctx = f"thermal_loads.body_source[{i}]"
        if not isinstance(src, dict):
            errors.append(f"{ctx}: cada fuente debe ser un diccionario.")
            continue
        unknown = set(src) - {"Q", "value", "elements"}
        if unknown:
            errors.append(f"{ctx}: claves no admitidas {sorted(unknown)}; se admiten 'Q' y "
                          "'elements'.")
            continue
        ids = src.get("elements")
        if ids is None:
            elements = list(domain.elements.values())
        else:
            missing = [e for e in ids if e not in domain.elements]
            if missing:
                errors.append(f"{ctx}: elementos inexistentes {missing}.")
                continue
            elements = [domain.elements[e] for e in ids]
        try:
            out.append((elements, load_field(src.get("Q", src.get("value", 0.0)), 1,
                                             normal=False, context=ctx)))
        except ExpressionError as exc:
            errors.append(str(exc))
    return out, errors
