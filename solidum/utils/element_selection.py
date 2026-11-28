"""Selectores de elementos del YAML: los de la historia y las energías por zona
(ADR 0036, 0038) y los de las cargas de línea (ADR 0044).

Un selector es un diccionario con una de las claves ``id``, ``ids``, ``group``,
``where`` o ``interface``:

- ``id`` / ``ids``: los elementos de esos números;
- ``where``: los que tienen el centro donde se cumple la condición (expresión
  de ``x``, ``y``, ``z`` con comparaciones, ``and``, ``or``, ``not`` y ``near``);
- ``interface``: los elementos de la interfaz de ese nombre (ADR 0043 §7);
- ``group``: en la historia, los elementos con todos sus nodos en el grupo
  físico (:func:`select_elements`); en las cargas de línea, los que coinciden
  con una **celda** del grupo (:func:`select_by_cells`), como las tracciones
  (ADR 0023 §2).
"""
from __future__ import annotations

import numpy as np


def pad3(c) -> list:
    """Coordenadas con tres componentes (las que faltan, cero)."""
    c = [float(v) for v in c]
    return c + [0.0] * (3 - len(c))


def select_elements(dom, sel: dict, ctx: str) -> list:
    """Elementos de un selector ``id`` | ``ids`` | ``group`` (todos sus nodos
    en el grupo físico) | ``where`` (condición sobre el centro) |
    ``interface`` (los elementos de la interfaz de ese nombre, ADR 0043 §7);
    sin selector, todos."""
    if "id" in sel:
        ids = [sel["id"]]
    elif "ids" in sel:
        ids = list(sel["ids"])
    elif "interface" in sel:
        interfaces = getattr(dom, "interfaces", None) or {}
        if sel["interface"] not in interfaces:
            raise ValueError(f"{ctx}: interfaz {sel['interface']!r} inexistente. Disponibles (bloque "
                             f"'interfaces'): {sorted(interfaces)}.")
        ids = list(interfaces[sel["interface"]])
    elif "group" in sel:
        groups = getattr(dom, "physical_groups", None) or {}
        if sel["group"] not in groups:
            raise ValueError(f"{ctx}: grupo físico {sel['group']!r} inexistente. "
                             f"Disponibles: {sorted(groups)}.")
        inside = set(groups[sel["group"]])
        ids = [eid for eid, el in dom.elements.items()
               if all(nd.id in inside for nd in el.nodes)]
        if not ids:
            raise ValueError(f"{ctx}: ningún elemento tiene todos sus nodos en el grupo "
                             f"{sel['group']!r}.")
    elif "where" in sel:
        ids = _where(dom, sel["where"], ctx)
    else:
        ids = list(dom.elements)
    missing = [i for i in ids if i not in dom.elements]
    if missing:
        raise ValueError(f"{ctx}: elementos inexistentes {missing[:10]}.")
    return [dom.elements[i] for i in ids]


def _where(dom, condition: str, ctx: str) -> list:
    from solidum.utils.expressions import POSITION, compile_predicate
    pts = np.array([pad3(nd.coordinates) for nd in dom.nodes.values()])
    scale = float(np.linalg.norm(pts.max(axis=0) - pts.min(axis=0)))

    def near(a, b, tol=None):
        return abs(a - b) <= (1e-6 * scale if tol is None else tol)
    pred = compile_predicate(condition, POSITION, {"near": near})
    ids = [eid for eid, el in dom.elements.items()
           if pred(*np.mean([pad3(nd.coordinates) for nd in el.nodes], axis=0))]
    if not ids:
        raise ValueError(f"{ctx}: ningún elemento tiene el centro donde {condition!r}.")
    return ids


def select_by_cells(dom, group: str, accepts, ctx: str) -> list:
    """Elementos que coinciden por sus nodos con una celda del grupo físico
    ``group`` y cumplen ``accepts(elemento)`` (ADR 0044 §4): error si el
    grupo no existe, si ninguna celda coincide o si una coincide con dos
    elementos (se cargaría dos veces)."""
    cells = getattr(dom, "physical_group_cells", None)
    if not cells:
        raise ValueError(f"{ctx}: el modelo no tiene grupos físicos con celdas (los aportan una "
                         "malla gmsh vía 'mesh' o una malla por bloques).")
    if group not in cells:
        raise ValueError(f"{ctx}: grupo físico {group!r} inexistente. Disponibles: "
                         f"{sorted(cells)}.")
    index: dict = {}
    for el in dom.elements.values():
        index.setdefault(frozenset(nd.id for nd in el.nodes), []).append(el)
    found, double, refused = [], 0, []
    for cell in cells[group]:
        hits = index.get(frozenset(cell), [])
        ok = [el for el in hits if accepts(el)]
        if len(ok) == 1:
            found.append(ok[0])
        elif ok:
            double += 1
        elif hits:
            refused.extend(hits)
    if refused:
        kinds = sorted({type(el).__name__ for el in refused})
        raise ValueError(f"{ctx}: los elementos {[el.id for el in refused][:10]} "
                         f"({', '.join(kinds)}) del grupo {group!r} no admiten cargas de línea.")
    if double:
        raise ValueError(f"{ctx}: {double} celda(s) del grupo {group!r} coinciden con dos "
                         "elementos de línea a la vez.")
    if not found:
        raise ValueError(f"{ctx}: ninguna celda del grupo {group!r} es un elemento que admita "
                         "cargas de línea.")
    return found
