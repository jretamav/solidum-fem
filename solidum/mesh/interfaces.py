"""Creación de elementos de interfaz de espesor nulo desde la malla (ADR 0043 §6).

Dos formas, las dos sobre grupos de la malla (los de Gmsh o los de una malla
estructurada):

- :func:`pair_interfaces` une dos grupos **de contorno** de dos cuerpos,
  con nodos coincidentes o con holgura: cada nodo de la cara A se empareja con
  el nodo de B que está sobre la normal de A. ``n`` apunta de A hacia B (la
  normal exterior del sólido del lado A).
- :func:`split_interfaces` separa los lados (2D) o caras (3D) **interiores**
  de uno o varios grupos: duplica sus nodos y pone un elemento de interfaz en
  cada lado. Para cada nodo, los sólidos que lo tocan se separan en
  componentes conexas al cortar los lados del grupo y cada componente recibe
  su copia: la punta de una grieta no se duplica y una unión de k ramas da k
  copias.

Lados y caras del grupo: si la malla da las celdas del grupo (Gmsh), las
celdas que son lados de algún sólido; si no (los grupos de una malla
estructurada sólo llevan celdas del contorno libre), los lados con todos sus
nodos en el grupo.

Contabilidad (ADR 0043 §6): las copias son nodos del sólido; los grupos de
nodos que tenían el original reciben también las copias; las celdas de los
grupos pasan al lado dueño (las de un lado separado, a sus dos labios); se
crean los grupos ``<name>_a`` y ``<name>_b`` con las caras de cada lado. El
dominio guarda las copias en ``Domain.split_nodes``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from solidum.core.dof_names import TEMPERATURE
from solidum.core.interface_material import InterfaceMaterial
from solidum.core.node import Node
from solidum.core.point_quantity import FLUX, JUMP, STRAIN, point_quantity
from solidum.elements.interface import INTEGRATIONS, InterfaceElement
from solidum.elements.interface._shared import GAP_TOLERANCE
from solidum.logging import get_logger
from solidum.user import user_module_of

_log = get_logger("mesh.interfaces")

_SIDES = {2: "EDGE_NODES", 3: "FACE_NODES"}
_SIDE_WORD = {2: "lado", 3: "cara"}
_SHOW = 5           # ids que se muestran en un mensaje


# --- elemento de cada cara (despacho por el registro, ADR 0043 §1 y §6) -----------

def _interface_classes(main_only: bool) -> list:
    """Clases de interfaz de espesor nulo registradas (``JUMP`` y subclases de
    :class:`InterfaceElement`); con ``main_only``, sólo las del programa
    principal (un módulo de usuario se elige con ``element``)."""
    from solidum.registry import ElementRegistry
    out = []
    for name in ElementRegistry.names():
        cls = ElementRegistry.get(name)
        if not (isinstance(cls, type) and issubclass(cls, InterfaceElement)
                and point_quantity(cls).kind == JUMP):
            continue
        if main_only and user_module_of(cls) is not None:
            continue
        out.append(cls)
    return out


def _interface_class(dim: int, n_face: int, element=None):
    """Elemento de interfaz para una cara de ``n_face`` nodos en ``dim``
    dimensiones: el declarado en ``element`` (nombre registrado o clase) o, por
    omisión, la única clase del programa principal que encaja."""
    word = _SIDE_WORD[dim]
    if element is not None:
        from solidum.registry import ElementRegistry
        cls = ElementRegistry.get(element) if isinstance(element, str) else element
        if not (isinstance(cls, type) and issubclass(cls, InterfaceElement)):
            raise ValueError(f"element={element!r} no es un elemento de interfaz de espesor nulo "
                             f"(subclase de InterfaceElement).")
        if cls.JUMP_DIM != dim or cls.N_FACE_NODES != n_face:
            raise ValueError(
                f"element={cls.__name__} une {_SIDE_WORD.get(cls.JUMP_DIM, word)}s de "
                f"{cls.N_FACE_NODES} nodos en {cls.JUMP_DIM}D; los {word}s del grupo tienen "
                f"{n_face} nodos en {dim}D.")
        return cls
    fits = [c for c in _interface_classes(True) if c.JUMP_DIM == dim and c.N_FACE_NODES == n_face]
    if len(fits) == 1:
        return fits[0]
    if not fits:
        known = sorted({c.N_FACE_NODES for c in _interface_classes(True) if c.JUMP_DIM == dim})
        raise ValueError(
            f"los {word}s del grupo tienen {n_face} nodos: hay elementos de interfaz {dim}D para "
            f"{word}s de {known} nodos; los de {n_face} no existen todavía.")
    raise ValueError(f"hay varios elementos de interfaz para {word}s de {n_face} nodos "
                     f"({sorted(c.__name__ for c in fits)}): declare cuál con 'element'.")


def _check_integration(integration) -> None:
    if integration not in INTEGRATIONS:
        raise ValueError(f"integration debe ser uno de {list(INTEGRATIONS)} "
                         f"(recibido {integration!r}).")


@dataclass(frozen=True)
class InterfaceInsertion:
    """Lo que creó una entrada de interfaces."""
    name: str
    kind: str                                   # "pair" o "split"
    elements: tuple[int, ...]                   # ids de los elementos de interfaz
    copies: dict = field(default_factory=dict)  # id original → ids de sus copias (split)


# --- topología común -------------------------------------------------------------

def _coords(node, dim: int) -> np.ndarray:
    c = np.zeros(dim)
    x = np.asarray(node.coordinates, dtype=float)[:dim]
    c[:len(x)] = x
    return c


def _solids(domain) -> tuple[int, list]:
    """Dimensión del modelo y sus sólidos mecánicos con lados declarados."""
    found = {2: [], 3: []}
    for el in domain.elements.values():
        pq = point_quantity(el)
        if pq.kind != STRAIN or pq.dim not in (3, 6):
            continue
        dim = 2 if pq.dim == 3 else 3
        if getattr(el, _SIDES[dim], None):
            found[dim].append(el)
    dims = [d for d in (2, 3) if found[d]]
    if len(dims) != 1:
        raise ValueError("las interfaces necesitan un modelo de sólidos 2D o 3D con lados "
                         "declarados (EDGE_NODES o FACE_NODES)"
                         + ("; hay de los dos." if len(dims) == 2 else "; no hay ninguno."))
    return dims[0], found[dims[0]]


def _side_index(solids, dim: int) -> dict:
    """``frozenset(ids)`` → [(elemento, índice local del lado)]."""
    index: dict = {}
    for el in solids:
        for k, local in enumerate(getattr(el, _SIDES[dim])):
            index.setdefault(frozenset(el.nodes[j].id for j in local), []).append((el, k))
    return index


def _side_ids(el, k: int, dim: int) -> list[int]:
    return [el.nodes[j].id for j in getattr(el, _SIDES[dim])[k]]


def _oriented(el, k: int, dim: int) -> list[int]:
    """Nodos del lado ``k`` de ``el`` como cara A: ``n`` = normal exterior de
    ``el``. En 2D, al revés que en el contorno antihorario del elemento
    (``n = e_z × s``); en 3D, el orden de ``FACE_NODES`` (antihorario visto
    desde fuera)."""
    ids = _side_ids(el, k, dim)
    if dim == 3:
        return ids
    ids = ids[::-1]
    # Lado cuadrático (extremo, intermedio, extremo): la cara de la interfaz
    # lleva los extremos primero y el intermedio al final.
    return [ids[0], ids[2], ids[1]] if len(ids) == 3 else ids


def _face_normal(X: np.ndarray, dim: int) -> tuple[np.ndarray, float]:
    """Normal unitaria de la cara con nodos ``X`` (en su orden) y tamaño
    ``h`` (longitud, o raíz del área)."""
    if dim == 2:
        a = X[1] - X[0]
        L = float(np.linalg.norm(a))
        return np.array([-a[1], a[0]]) / L, L
    if len(X) == 3:
        c = np.cross(X[1] - X[0], X[2] - X[0])
    else:
        c = np.cross(X[2] - X[0], X[3] - X[1])
    norm = float(np.linalg.norm(c))
    return c / norm, float(np.sqrt(0.5 * norm))


def _group_faces(domain, groups: list[str], index: dict, solids, dim: int) -> list[frozenset]:
    """Lados (2D) o caras (3D) de los grupos ``groups``."""
    phys = getattr(domain, "physical_groups", None) or {}
    cells_all = getattr(domain, "physical_group_cells", None) or {}
    if not phys:
        raise ValueError("el modelo no tiene grupos (los da una malla de Gmsh o una malla "
                         "estructurada vía 'mesh').")
    nodes, cells = set(), set()
    for g in groups:
        if g not in phys:
            raise ValueError(f"grupo {g!r} inexistente. Disponibles: {sorted(phys)}.")
        nodes |= set(phys[g])
        cells |= {frozenset(c) for c in cells_all.get(g, ())}
    word = _SIDE_WORD[dim]
    candidates = [key for key in index if key <= nodes]
    given = [key for key in candidates if key in cells]
    if given:
        faces = given
    else:
        region = next((el for el in solids if all(n.id in nodes for n in el.nodes)), None)
        if region is not None:
            raise ValueError(
                f"el grupo {', '.join(map(repr, groups))} contiene todos los nodos del elemento "
                f"{region.id}: es una región, no un{'a' if dim == 3 else ''} "
                f"{'superficie' if dim == 3 else 'línea'} de {word}s.")
        faces = candidates
    if not faces:
        raise ValueError(
            f"el grupo {', '.join(map(repr, groups))} no contiene ningún {word} de los sólidos "
            f"(¿es un grupo de puntos o de otra dimensión? En Gmsh, una "
            f"{'superficie' if dim == 3 else 'curva'} interior debe estar incrustada en la "
            f"{'región' if dim == 3 else 'superficie'}: "
            f"{'Surface{…} In Volume{…}' if dim == 3 else 'Curve{…} In Surface{…}'}).")
    return sorted(faces, key=lambda f: sorted(f))


def _classes_for(faces, dim: int, element) -> dict:
    """Clase de interfaz por número de nodos de cara (errores antes de
    cambiar nada)."""
    return {nf: _interface_class(dim, nf, element) for nf in sorted({len(f) for f in faces})}


def _check_material(material, dim: int) -> None:
    if not isinstance(material, InterfaceMaterial) or not material.accepts("surface", dim):
        raise ValueError(
            f"la ley {type(material).__name__} no es un material de interfaz de superficie con "
            f"{dim} componentes (interface_materials con GEOMETRY='surface').")


def _check_new_groups(domain, name: str) -> list[str]:
    if name in getattr(domain, "interfaces", {}):
        raise ValueError(f"ya hay una interfaz {name!r}; elija otro 'name'.")
    names = [f"{name}_a", f"{name}_b"]
    taken = [n for n in names if n in (getattr(domain, "physical_groups", None) or {})]
    if taken:
        raise ValueError(f"los grupos {taken} ya existen; elija otro 'name' para la interfaz.")
    return names


def _thickness(ea, eb, dim: int) -> float:
    """Espesor de una interfaz 2D: el de los sólidos de los dos lados, que
    deben tener el mismo (spec ``Interface2D4`` §6)."""
    ta, tb = float(getattr(ea, "thickness", 1.0)), float(getattr(eb, "thickness", 1.0))
    if dim == 2 and not math.isclose(ta, tb, rel_tol=1e-12):
        raise ValueError(
            f"los sólidos de los dos lados tienen espesores distintos ({ta:g} en el elemento "
            f"{ea.id}, {tb:g} en el {eb.id}): el espesor de la interfaz no está definido.")
    return ta


def _build(cls, eid: int, nodes_a: list, nodes_b_by_a: list, material, integration: str,
           thickness: float):
    """Elemento de interfaz con la cara A y, para cada nodo de A, su par en B."""
    nodes = [None] * (2 * cls.N_FACE_NODES)
    for i in range(cls.N_FACE_NODES):
        nodes[cls.A_NODES[i]] = nodes_a[i]
        nodes[cls.B_NODES[i]] = nodes_b_by_a[i]
    kwargs = {"integration": integration}
    if cls.JUMP_DIM == 2:
        kwargs["thickness"] = thickness
    return cls(eid, nodes, material, **kwargs)


def _fmt(ids, domain, dim) -> str:
    shown = ", ".join(f"{i} {[round(float(v), 6) for v in _coords(domain.nodes[i], dim)]}"
                      for i in list(ids)[:_SHOW])
    return shown + (", …" if len(ids) > _SHOW else "")


def _add_groups(domain, names, faces_a, faces_b) -> None:
    for gname, faces in zip(names, (faces_a, faces_b)):
        domain.physical_groups[gname] = sorted({i for f in faces for i in f})
        domain.physical_group_cells[gname] = [tuple(f) for f in faces]


# --- pair ------------------------------------------------------------------------

def pair_interfaces(domain, group_a: str, group_b: str, material, *, name: str,
                    integration: str = "nodal", element=None) -> InterfaceInsertion:
    """Une con elementos de interfaz los grupos de contorno ``group_a`` y
    ``group_b`` de dos cuerpos.

    Cada nodo de la cara A se empareja con el nodo de B que está sobre la
    normal de A (desfase tangencial ≤ 10⁻⁶·h, la tolerancia del elemento);
    los nodos pueden coincidir (holgura nula) o estar separados (holgura
    inicial ``g_0``, que sólo admite una ley con ``ACCEPTS_INITIAL_GAP``). La
    normal y la medida son las de la cara A: con caras separadas, A debe ser
    la plana o la más regular (spec ``Interface2D4`` §2). En 2D el espesor es
    el de los sólidos de los dos lados, que debe ser el mismo.

    Parameters
    ----------
    domain : Domain
        Con grupos (``physical_groups`` y ``physical_group_cells``).
    group_a, group_b : str
        Grupos de contorno de la cara A y de la cara B.
    material : InterfaceMaterial
        Ley de interfaz de superficie.
    name : str
        Nombre de la interfaz: se crean los grupos ``<name>_a`` y ``<name>_b``.
    integration : {"nodal", "gauss"}
    element : str o clase, opcional
        Elemento de interfaz (nombre registrado o clase, p. ej. de un módulo de
        usuario). Por omisión, el del programa principal para esas caras.

    Returns
    -------
    InterfaceInsertion

    Raises
    ------
    ValueError
        Grupo inexistente o de otra dimensión, cara interior, nodos sin pareja,
        emparejamiento que no es uno a uno, cara de A sin cara de B enfrente,
        caras que no se miran, interpenetración inicial, holgura con una ley
        que no la admite o caras cuadráticas. No cambia el dominio si falla.
    """
    dim, solids = _solids(domain)
    _check_material(material, dim)
    _check_integration(integration)
    names = _check_new_groups(domain, name)
    index = _side_index(solids, dim)
    word = _SIDE_WORD[dim]
    faces_a = _group_faces(domain, [group_a], index, solids, dim)
    faces_b = _group_faces(domain, [group_b], index, solids, dim)
    classes = _classes_for(faces_a, dim, element)
    for faces, g in ((faces_a, group_a), (faces_b, group_b)):
        interior = [f for f in faces if len(index[f]) > 1]
        if interior:
            raise ValueError(
                f"{len(interior)} {word}(s) del grupo {g!r} son interiores (los comparten dos "
                f"elementos): pair une {word}s del contorno de dos cuerpos; un {word} interior "
                f"se separa con split.")
    if set(faces_a) & set(faces_b):
        raise ValueError(f"los grupos {group_a!r} y {group_b!r} comparten {word}s.")

    b_ids = sorted({i for f in faces_b for i in f})
    XB = np.array([_coords(domain.nodes[i], dim) for i in b_ids])
    partner: dict[int, int] = {}
    unmatched, conflicts = set(), set()
    oriented = []
    for key in faces_a:
        el, k = index[key][0]
        ids = _oriented(el, k, dim)
        X = np.array([_coords(domain.nodes[i], dim) for i in ids])
        n, h = _face_normal(X, dim)
        tol = GAP_TOLERANCE * h
        for i, x in zip(ids, X):
            d = XB - x
            dn = d @ n
            tang = np.linalg.norm(d - dn[:, None] * n, axis=1)
            cand = np.flatnonzero(tang <= tol)
            if cand.size == 0:
                unmatched.add(i)
                continue
            b = b_ids[int(cand[np.argmin(np.abs(dn[cand]))])]
            if partner.setdefault(i, b) != b:
                conflicts.add(i)
        oriented.append((el, ids, n))
    if unmatched:
        raise ValueError(
            f"{len(unmatched)} nodo(s) de {group_a!r} sin pareja en {group_b!r} sobre la normal "
            f"de su {word} (desfase tangencial > 10⁻⁶·h): {_fmt(sorted(unmatched), domain, dim)}.")
    if conflicts:
        raise ValueError(
            f"{len(conflicts)} nodo(s) de {group_a!r} tienen parejas distintas según el {word} "
            f"desde el que se buscan: {_fmt(sorted(conflicts), domain, dim)}.")
    taken: dict[int, list] = {}
    for a, b in partner.items():
        taken.setdefault(b, []).append(a)
    shared = {b: a for b, a in taken.items() if len(a) > 1}
    if shared:
        b, a = next(iter(sorted(shared.items())))
        raise ValueError(
            f"el emparejamiento no es uno a uno: {len(shared)} nodo(s) de {group_b!r} son pareja "
            f"de varios nodos de {group_a!r} (p. ej. {b} de {sorted(a)}).")

    set_b = set(faces_b)
    missing, facing = [], []
    plan = []
    for el, ids, n in oriented:
        partners = [partner[i] for i in ids]
        key_b = frozenset(partners)
        if key_b not in set_b:
            missing.append(ids)
            continue
        el_b, k_b = index[key_b][0]
        Xb = np.array([_coords(domain.nodes[i], dim) for i in _oriented(el_b, k_b, dim)])
        if float(_face_normal(Xb, dim)[0] @ n) >= 0.0:
            facing.append(ids)
            continue
        plan.append((el, el_b, ids, partners))
    if missing:
        raise ValueError(
            f"{len(missing)} {word}(s) de {group_a!r} sin {word} de {group_b!r} enfrente (los "
            f"pares de sus nodos no forman un {word} del grupo); el primero tiene los nodos "
            f"{missing[0]}.")
    if facing:
        raise ValueError(
            f"{len(facing)} {word}(s) de {group_a!r} no miran a su {word} de {group_b!r} (las dos "
            f"normales exteriores apuntan al mismo lado: ¿son los grupos de las caras que se "
            f"tocan?); el primero tiene los nodos {facing[0]}.")

    eid = max(domain.elements, default=0) + 1
    created, errors = [], []
    for el, el_b, ids, partners in plan:
        cls = classes[len(ids)]
        try:
            created.append(_build(cls, eid + len(created), [domain.nodes[i] for i in ids],
                                  [domain.nodes[i] for i in partners], material, integration,
                                  _thickness(el, el_b, dim)))
        except ValueError as exc:
            errors.append(f"{word} {ids} – {partners}: {exc}")
    if errors:
        raise ValueError(f"{len(errors)} elemento(s) de interfaz no se pueden crear; "
                         + "; ".join(errors[:_SHOW]) + (" …" if len(errors) > _SHOW else ""))
    for e in created:
        domain.add_element(e)
    domain.interfaces[name] = tuple(e.id for e in created)
    _add_groups(domain, names, [ids for _, _, ids, _ in plan], [p for _, _, _, p in plan])
    gap = max((float(np.max(e.g0)) for e in created), default=0.0)
    _log.info(f"Interfaz {name!r} (pair {group_a!r}–{group_b!r}): {len(created)} elementos "
              f"{type(created[0]).__name__}, holgura inicial máxima {gap:.6g}.")
    return InterfaceInsertion(name, "pair", tuple(e.id for e in created))


# --- split -----------------------------------------------------------------------

def _components(p: int, elements: list, index: dict, cut: set, dim: int) -> list[list]:
    """Componentes conexas de los sólidos que tocan el nodo ``p`` al cortar los
    lados de ``cut``: dos elementos están unidos si comparten un lado que
    contiene a ``p`` y no está cortado. Ordenadas por el menor id."""
    parent = {el.id: el.id for el in elements}

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for el in elements:
        for k, local in enumerate(getattr(el, _SIDES[dim])):
            key = frozenset(el.nodes[j].id for j in local)
            if p not in key or key in cut:
                continue
            for other, _ in index[key]:
                if other.id in parent and other is not el:
                    parent[root(other.id)] = root(el.id)
    comps: dict[int, list] = {}
    for el in elements:
        comps.setdefault(root(el.id), []).append(el)
    return sorted((sorted(c, key=lambda e: e.id) for c in comps.values()), key=lambda c: c[0].id)


def split_interfaces(domain, groups, material, *, name: str,
                     integration: str = "nodal", element=None) -> InterfaceInsertion:
    """Separa los lados (2D) o caras (3D) interiores de ``groups`` y pone un
    elemento de interfaz en cada uno.

    Duplica los nodos que lo necesitan (ver el módulo): la punta no se
    duplica, una unión de k ramas da k copias. Varios grupos que se cruzan
    van en la misma llamada (una lista), para que la unión se trate de una
    vez. Las copias tienen las coordenadas del original y los ids siguientes
    al mayor del dominio, en orden de nodo original. El lado A de cada cara se
    elige por su geometría (la componente mayor de ``n`` positiva), así que
    no depende de la orientación de las celdas en el archivo de malla. Si el
    dominio ya estaba numerado, se vuelve a numerar.

    Parameters
    ----------
    domain : Domain
    groups : str o lista de str
        Grupos de lados o caras interiores.
    material : InterfaceMaterial
        Ley de interfaz de superficie.
    name : str
        Nombre: se crean los grupos ``<name>_a`` y ``<name>_b`` (los dos labios).
    integration : {"nodal", "gauss"}
    element : str o clase, opcional
        Elemento de interfaz (nombre registrado o clase). Por omisión, el del
        programa principal para esas caras.

    Returns
    -------
    InterfaceInsertion
        Con ``copies``: id original → ids de todas sus copias.

    Raises
    ------
    ValueError
        Grupo inexistente o de otra dimensión, lado del contorno, caras sin
        elemento de interfaz, ningún nodo que separar, un nodo que se
        duplicaría y lo toca un elemento que no es un sólido mecánico (barra,
        marco, térmico, otra interfaz) o ya tiene un apoyo, una restricción
        multipunto o una masa (no se sabe a qué lado pertenecen), o una cara
        cuyo elemento no se puede construir. No cambia el dominio si falla.
    """
    groups = [groups] if isinstance(groups, str) else list(groups)
    dim, solids = _solids(domain)
    _check_material(material, dim)
    _check_integration(integration)
    names = _check_new_groups(domain, name)
    index = _side_index(solids, dim)
    word = _SIDE_WORD[dim]
    faces = _group_faces(domain, groups, index, solids, dim)
    classes = _classes_for(faces, dim, element)
    boundary = [f for f in faces if len(index[f]) == 1]
    if boundary:
        raise ValueError(
            f"{len(boundary)} {word}(s) del grupo están en el contorno (un solo elemento): split "
            f"separa {word}s interiores; el primero tiene los nodos {sorted(boundary[0])}.")
    shared = [f for f in faces if len(index[f]) > 2]
    if shared:
        raise ValueError(
            f"{len(shared)} {word}(s) del grupo los comparten más de dos elementos (malla no "
            f"conforme o regiones superpuestas); el primero tiene los nodos {sorted(shared[0])}.")
    for f in faces:
        (e1, _), (e2, _) = index[f]
        _thickness(e1, e2, dim)
    cut = set(faces)
    P = sorted({i for f in faces for i in f})
    touching: dict[int, list] = {p: [] for p in P}
    for el in domain.elements.values():
        for n in el.nodes:
            if n.id in touching and el not in touching[n.id]:
                touching[n.id].append(el)
    solid_ids = {el.id for el in solids}

    plan, foreign = {}, []
    for p in P:
        comps = _components(p, [e for e in touching[p] if e.id in solid_ids], index, cut, dim)
        if len(comps) > 1:
            plan[p] = comps
            foreign += [(p, e) for e in touching[p] if e.id not in solid_ids]
    if not plan:
        raise ValueError(f"el grupo no separa ningún nodo: todos sus nodos son extremos (puntas) "
                         f"de los {word}s.")
    if foreign:
        thermal = [(p, e) for p, e in foreign
                   if point_quantity(e).kind == FLUX or TEMPERATURE in getattr(e, "DOF_NAMES", ())]
        if thermal:
            p, e = thermal[0]
            raise ValueError(
                f"{len(thermal)} elemento(s) térmicos tocan nodos que se duplicarían (p. ej. "
                f"{type(e).__name__} {e.id} en el nodo {p}): separar el campo de temperatura en la "
                f"interfaz no está disponible (el contacto térmico queda fuera de este bloque).")
        p, e = foreign[0]
        raise ValueError(
            f"{len(foreign)} elemento(s) que no son sólidos mecánicos tocan nodos que se "
            f"duplicarían (p. ej. {type(e).__name__} {e.id} en el nodo {p}): no se sabe a qué "
            f"lado de la interfaz pertenecen.")
    _check_no_references(domain, sorted(plan))

    # Plan de la duplicación, sin tocar el dominio: ids de las copias y nodos
    # nuevos de cada elemento.
    next_id = max(domain.nodes) + 1
    copies: dict[int, tuple[int, ...]] = {}
    remap: dict[int, dict[int, int]] = {}
    for p in P:
        if p not in plan:
            continue
        ids = [p]
        for comp in plan[p][1:]:
            ids.append(next_id)
            for el in comp:
                remap.setdefault(el.id, {})[p] = next_id
            next_id += 1
        copies[p] = tuple(ids)
    original = {c: p for p, cs in copies.items() for c in cs}

    # Caras: A es el lado cuya normal exterior tiene positiva su componente
    # mayor. Un lado sin ningún nodo separado (sus extremos son puntas) no
    # lleva elemento: sus dos caras son los mismos nodos.
    def after(el, i):
        return remap.get(el.id, {}).get(i, i)

    inert = [f for f in faces if not f & copies.keys()]
    active = [f for f in faces if f & copies.keys()]
    faces_plan = []
    for key in active:
        (e1, k1), (e2, k2) = index[key]
        ids1 = _oriented(e1, k1, dim)
        n, _ = _face_normal(np.array([_coords(domain.nodes[i], dim) for i in ids1]), dim)
        ax = int(np.argmax(np.abs(n)))
        (ea, ka), eb = ((e1, k1), e2) if n[ax] > 0.0 else ((e2, k2), e1)
        orig_a = _oriented(ea, ka, dim)
        faces_plan.append((key, ea, eb, [after(ea, i) for i in orig_a],
                           [after(eb, i) for i in orig_a]))

    # Todos los elementos, construidos con nodos provisionales en las mismas
    # posiciones: si alguno falla (cara degenerada, desfase), nada cambió.
    def provisional(i):
        return Node(i, list(domain.nodes[original.get(i, i)].coordinates))

    eid = max(domain.elements, default=0) + 1
    errors = []
    for j, (key, ea, eb, ids_a, ids_b) in enumerate(faces_plan):
        try:
            _build(classes[len(ids_a)], eid + j, [provisional(i) for i in ids_a],
                   [provisional(i) for i in ids_b], material, integration, _thickness(ea, eb, dim))
        except ValueError as exc:
            errors.append(f"{word} {sorted(key)}: {exc}")
    if errors:
        raise ValueError(f"{len(errors)} elemento(s) de interfaz no se pueden crear; "
                         + "; ".join(errors[:_SHOW]) + (" …" if len(errors) > _SHOW else ""))

    # Mutación (todo validado).
    for p, ids in copies.items():
        for c in ids[1:]:
            domain.add_node(c, list(domain.nodes[p].coordinates))
    for el_id in sorted(remap):
        domain.reassign_element_nodes(el_id, remap[el_id])
    domain.record_split_nodes(copies)
    if inert:
        _log.warning(f"Interfaz {name!r}: {len(inert)} {word}(s) del grupo no separan ningún "
                     f"nodo (todos sus nodos son puntas) y no llevan elemento; el primero tiene "
                     f"los nodos {sorted(inert[0])}.")
    created = [_build(classes[len(ids_a)], eid + j, [domain.nodes[i] for i in ids_a],
                      [domain.nodes[i] for i in ids_b], material, integration,
                      _thickness(ea, eb, dim))
               for j, (_, ea, eb, ids_a, ids_b) in enumerate(faces_plan)]
    for e in created:
        domain.add_element(e)
    domain.interfaces[name] = tuple(e.id for e in created)
    lips_a = [ids_a for _, _, _, ids_a, _ in faces_plan]
    lips_b = [ids_b for _, _, _, _, ids_b in faces_plan]
    _remap_groups(domain, solids, copies, original, index, cut,
                  dict(zip(active, zip(lips_a, lips_b))))
    _add_groups(domain, names, lips_a, lips_b)
    domain.renumber_if_numbered()
    n_copies = sum(len(c) - 1 for c in copies.values())
    _log.info(f"Interfaz {name!r} (split {', '.join(map(repr, groups))}): {len(created)} elementos "
              f"{type(created[0]).__name__}, {n_copies} nodos nuevos.")
    return InterfaceInsertion(name, "split", tuple(e.id for e in created), dict(copies))


def _check_no_references(domain, node_ids) -> None:
    """Los nodos que se duplicarían no pueden tener ya un apoyo, una
    restricción multipunto o una masa, ni ser nodos generados (una barra
    embebida): quedarían en un solo lado sin aviso. En el YAML el bloque
    ``interfaces`` se aplica antes que todo eso; con la API de Python, cree
    las interfaces antes de declararlos."""
    ids = set(node_ids)
    constrained = {n for c in domain.linear_constraints
                   for n in [c["slave"][0], *(m[0] for m in c["masters"])]}
    found = {
        "apoyos": sorted(i for i in ids if domain.nodes[i].boundary_conditions),
        "restricciones multipunto": sorted(ids & constrained),
        "masas": sorted(ids & {nid for nid, _ in getattr(domain, "nodal_masses", {})}),
        "nodos generados (barras embebidas)": sorted(ids & set(domain.generated_nodes)),
    }
    found = {k: v for k, v in found.items() if v}
    if found:
        detail = "; ".join(f"{k} {v[:_SHOW]}{' …' if len(v) > _SHOW else ''}"
                           for k, v in found.items())
        raise ValueError(
            f"nodos que se duplicarían ya tienen {detail}: no se sabe a qué lado de la interfaz "
            f"pertenecen. Cree las interfaces antes de declararlos (en el YAML, el bloque "
            f"'interfaces' ya se aplica antes).")


def _remap_groups(domain, solids, copies, original, index, cut, lips) -> None:
    """Grupos tras duplicar nodos: los de nodos reciben las copias; las celdas
    pasan al elemento dueño (un lado separado, a sus dos labios; una celda que
    es un elemento entero, a ese elemento)."""
    whole = {}
    for el in solids:
        whole.setdefault(frozenset(original.get(n.id, n.id) for n in el.nodes), el)
    for g, ids in domain.physical_groups.items():
        if any(i in copies for i in ids):
            domain.physical_groups[g] = sorted({c for i in ids for c in copies.get(i, (i,))})
    for g, cells in domain.physical_group_cells.items():
        if not any(i in copies for c in cells for i in c):
            continue
        out = []
        for c in cells:
            key = frozenset(c)
            if not key & copies.keys():
                out.append(c)
            elif key in cut:
                out.extend(tuple(lip) for lip in lips[key])
            elif key in index or key in whole:
                el = index[key][0][0] if key in index else whole[key]
                cur = {original.get(n.id, n.id): n.id for n in el.nodes}
                out.append(tuple(cur[i] for i in c))
            elif len(c) == 1:
                out.extend((i,) for i in copies[c[0]])
            else:
                out.append(c)
        domain.physical_group_cells[g] = out
