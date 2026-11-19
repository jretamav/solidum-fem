"""Contorno libre de una malla (ADR 0042).

Un lado (2D) o una cara (3D) es del **contorno libre** si pertenece a un solo
elemento: lo comparte con nadie más. Los lados y caras de cada elemento son los
que declara su clase en ``EDGE_NODES`` (2D) o ``FACE_NODES`` (3D); un elemento
que no los declara (barras, marcos, resortes, interfaces) no cuenta.

El filtro de campo deja contar sólo los elementos que tienen ciertos grados de
libertad: el contorno libre del sólido mecánico (``dofs=("ux", "uy")``) no lo
tapan los elementos térmicos montados sobre los mismos nodos.

Lo usan la generación de grupos de las mallas estructuradas
(:mod:`solidum.mesh.structured_blocks`) y cualquier componente que decida algo
con el borde libre del modelo entero (a través de ``Element.domain``).
"""
from __future__ import annotations

from typing import Iterable

_ATTR = {2: "EDGE_NODES", 3: "FACE_NODES"}


def free_boundary(elements: Iterable, dim: int,
                  dofs: Iterable[str] | None = None) -> list[tuple[int, ...]]:
    """Lados (``dim = 2``) o caras (``dim = 3``) del contorno libre.

    Parameters
    ----------
    elements : iterable de elementos
        Normalmente ``domain.elements.values()``: todos los del modelo, para
        que un lado compartido con un elemento de otra región no parezca libre.
    dim : {2, 3}
        Dimensión de las celdas: 2 → ``EDGE_NODES``, 3 → ``FACE_NODES``.
    dofs : iterable de str, optional
        Sólo cuentan los elementos cuyos ``DOF_NAMES`` incluyen todos estos
        grados de libertad (p. ej. ``("ux", "uy")``). Sin filtro, todos.

    Returns
    -------
    list[tuple[int, ...]]
        Ids de los nodos de cada lado o cara libre, en el orden local del
        elemento que la tiene (su orientación), en el orden en que aparecen los
        elementos.

    Raises
    ------
    ValueError
        Si ``dim`` no es 2 ni 3.
    """
    if dim not in _ATTR:
        raise ValueError(f"free_boundary: dim debe ser 2 (lados) o 3 (caras); recibido {dim!r}.")
    attr = _ATTR[dim]
    need = None if dofs is None else set(dofs)
    count: dict[frozenset, list] = {}
    for el in elements:
        if need is not None and not need <= set(getattr(el, "DOF_NAMES", ())):
            continue
        for local in getattr(el, attr, None) or ():
            ids = tuple(el.nodes[j].id for j in local)
            count.setdefault(frozenset(ids), []).append(ids)
    return [v[0] for v in count.values() if len(v) == 1]
