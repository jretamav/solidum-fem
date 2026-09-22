"""Post-proceso por lotes desde el dominio (ADR 0014 §9).

El exportador VTK y cualquier consumidor de ``compute_gauss_state`` sólo
ven el dominio, no el ``Assembler``. Cada familia deja en sus elementos
la referencia ``_batch_family`` al adoptar sus estados, así que desde el
dominio se puede recomponer qué elementos pertenecen a qué familia y
evaluar el post-proceso por arreglos (:meth:`Family.gauss_state`) en vez
de elemento a elemento. Los elementos sin familia (camino por elemento,
o dominio ensamblado con ``batch=False``) siguen usando su
``compute_gauss_state``: el resultado es el mismo, cambia el coste.
"""
from __future__ import annotations

import numpy as np

from solidum.math.batch.family import Family


def family_of(elem) -> Family | None:
    """Familia vigente del elemento, o ``None`` si se ensambla por elemento
    (o si el ensamblador que la creó ya la liberó)."""
    fam = getattr(elem, "_batch_family", None)
    if fam is None or not fam.adopted:
        return None
    return fam


def group_by_family(elements) -> tuple[dict[Family, list[int]], list[int]]:
    """Separa ``elements`` (secuencia) en ``{familia: [posiciones]}`` y la
    lista de posiciones sin familia. Las posiciones se refieren a
    ``elements``; dentro de cada familia se conserva el orden."""
    groups: dict[Family, list[int]] = {}
    loose: list[int] = []
    for pos, elem in enumerate(elements):
        fam = family_of(elem)
        if fam is None:
            loose.append(pos)
        else:
            groups.setdefault(fam, []).append(pos)
    return groups, loose


def _family_rows(fam: Family, members) -> np.ndarray:
    """Índice de fila en la familia de cada elemento de ``members``."""
    index = {id(e): i for i, e in enumerate(fam.elements)}
    return np.array([index[id(e)] for e in members], dtype=np.int64)


def gauss_states(domain, U: np.ndarray) -> dict[int, dict]:
    """``{elem_id: compute_gauss_state(U)}`` para todos los elementos del
    dominio que lo implementan, evaluando cada familia de una vez.

    Las entradas de los elementos en familia son vistas sobre los arreglos
    de la familia (no copias); las de los elementos sueltos son lo que
    devuelve su ``compute_gauss_state``. Un elemento sin ese método (1D)
    no aparece en el resultado.
    """
    elements = list(domain.elements.values())
    groups, loose = group_by_family(elements)
    out: dict[int, dict] = {}
    for fam, positions in groups.items():
        members = [elements[p] for p in positions]
        rows = _family_rows(fam, members)
        gs = fam.gauss_state(U)
        pg = gs["points_global"]
        for r, elem in zip(rows, members):
            entry = {"points_natural": gs["points_natural"],
                     "points_global": None if pg is None else pg[r]}
            for key in ("strain", "stress", "grad_T", "flux"):
                if key in gs:
                    entry[key] = gs[key][r]
            out[elem.id] = entry
    for p in loose:
        elem = elements[p]
        if hasattr(elem, "compute_gauss_state"):
            out[elem.id] = elem.compute_gauss_state(U)
    return out
