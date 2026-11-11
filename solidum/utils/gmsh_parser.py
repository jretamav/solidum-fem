# solidum_fem/solidum/utils/gmsh_parser.py
"""Lector de mallas de Gmsh (``.msh``) a un :class:`Domain`.

Celdas admitidas (nombres de meshio) y elemento por omisión:

====================  ==========  ==================
celda                 dimensión   elemento
====================  ==========  ==================
``triangle``          2           ``Tri3``
``triangle6``         2           ``Tri6``
``quad``              2           ``Quad4``
``quad8``             2           ``Quad8``
``quad9``             2           ``Quad9``
``tetra``             3           ``Tet4``
``tetra10``           3           ``Tet10``
``hexahedron``        3           ``Hex8``
``hexahedron20``      3           ``Hex20``
``hexahedron27``      3           ``Hex27``
``line``, ``line3``   1           ninguno (contorno)
====================  ==========  ==================

Las celdas de la dimensión de la malla (la mayor presente) son elementos. Las
de dimensión menor (líneas en 2D; triángulos y cuadriláteros en 3D) sólo son
contorno, para las cargas por grupo físico (ADR 0023), salvo que su grupo
físico pida un elemento explícito (``element`` en ``mesh_physical_groups``):
así entran, p. ej., barras de refuerzo dentro de un sólido.

**Orden de los nodos.** Cada celda trae sus nodos en el orden de meshio. Se
llevan al orden del elemento de Solidum casando sus coordenadas de referencia
(``_CELL_REFERENCE``, tomadas de Gmsh con la permutación de meshio aplicada)
con las ``REFERENCE_NODES`` del elemento (ADR 0022), que usan las mismas
celdas de referencia. Un elemento sin ``REFERENCE_NODES`` recibe los nodos en
el orden de la celda.

**Orientación.** Un elemento sólido con jacobiano negativo en el centro
(celda horaria en 2D, invertida en 3D) se reordena reflejando su celda de
referencia (ξ ↔ η); uno degenerado (jacobiano nulo) es un error.
"""
from __future__ import annotations

import numpy as np

from solidum.core.domain import Domain
from solidum.logging import get_logger
from solidum.registry import ElementRegistry

try:
    import meshio
except ImportError:
    raise ImportError("Falta la librería 'meshio'. Instálala con: pip install meshio")

_log = get_logger("parsers.gmsh")

_L = ((-1.0,), (1.0,))
_TRI = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
_QUAD = ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0))
_QUAD_MID = ((0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0))
_TET = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
_HEX = ((-1.0, -1.0, -1.0), (1.0, -1.0, -1.0), (1.0, 1.0, -1.0), (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0), (1.0, -1.0, 1.0), (1.0, 1.0, 1.0), (-1.0, 1.0, 1.0))
_HEX_MID = ((0.0, -1.0, -1.0), (1.0, 0.0, -1.0), (0.0, 1.0, -1.0), (-1.0, 0.0, -1.0),
            (0.0, -1.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0), (-1.0, 0.0, 1.0),
            (-1.0, -1.0, 0.0), (1.0, -1.0, 0.0), (1.0, 1.0, 0.0), (-1.0, 1.0, 0.0))

# Celda de meshio → (celda de referencia, coordenadas de referencia de sus
# nodos en el orden de meshio, elemento por omisión). Las coordenadas salen de
# ``gmsh.model.mesh.getElementProperties`` con ``_gmsh_to_meshio_order`` de
# meshio aplicado; tests/test_gmsh_parser_celdas.py lo comprueba con Gmsh.
_CELL_REFERENCE = {
    "line": ("line", _L, None),
    "line3": ("line", _L + ((0.0,),), None),
    "triangle": ("tri", _TRI, "Tri3"),
    "triangle6": ("tri", _TRI + ((0.5, 0.0), (0.5, 0.5), (0.0, 0.5)), "Tri6"),
    "quad": ("quad", _QUAD, "Quad4"),
    "quad8": ("quad", _QUAD + _QUAD_MID, "Quad8"),
    "quad9": ("quad", _QUAD + _QUAD_MID + ((0.0, 0.0),), "Quad9"),
    "tetra": ("tet", _TET, "Tet4"),
    "tetra10": ("tet", _TET + ((0.5, 0.0, 0.0), (0.5, 0.5, 0.0), (0.0, 0.5, 0.0),
                               (0.0, 0.0, 0.5), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)), "Tet10"),
    "hexahedron": ("hex", _HEX, "Hex8"),
    "hexahedron20": ("hex", _HEX + _HEX_MID, "Hex20"),
    "hexahedron27": ("hex", _HEX + _HEX_MID + (
        (-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, -1.0), (0.0, 0.0, 1.0), (0.0, 0.0, 0.0)), "Hex27"),
}
_CELL_DIM = {"line": 1, "tri": 2, "quad": 2, "tet": 3, "hex": 3}
_N_CORNERS = {"line": 2, "tri": 3, "quad": 4, "tet": 4, "hex": 8}

# Elemento por omisión de cada celda (compatibilidad: antes sólo quad y triangle).
GMSH_TYPE_MAP = {cell: default for cell, (_, _, default) in _CELL_REFERENCE.items()
                 if default is not None}


def _match(points_from, points_to, what: str) -> list[int]:
    """Índices ``j`` tales que ``points_from[j] == points_to[k]`` para cada ``k``."""
    a = np.asarray(points_from, dtype=float)
    perm = []
    for p in np.asarray(points_to, dtype=float):
        hits = np.flatnonzero(np.all(np.abs(a - p) < 1e-9, axis=1))
        if len(hits) != 1:
            raise ValueError(what)
        perm.append(int(hits[0]))
    return perm


class GmshParser:
    """Lector de mallas generadas en Gmsh (.msh) utilizando meshio."""

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.domain = Domain()
        self._perm_cache: dict = {}
        self.n_reoriented = 0

    def parse(self, default_spec: dict, group_specs: dict | None = None) -> Domain:
        """Lee el archivo de malla y retorna un objeto Domain poblado.

        Cada *Physical Group* puede llevar su propia receta de elemento; las
        celdas sin grupo con receta usan ``default_spec``. Una receta es un
        dict (ADR 0020, P7):

        ``element``
            Tipo registrado, o ``None`` para el que corresponde a la celda
            gmsh (``GMSH_TYPE_MAP``). Las celdas de dimensión menor que la de
            la malla sólo son elementos si su grupo nombra uno.
        ``kwargs``
            Argumentos del constructor que el elemento debe aceptar
            (``material`` ya resuelto, parámetros propios del elemento).
        ``optional``
            Argumentos que se pasan sólo si el constructor los acepta
            (``thickness``, ``quadrature``: ``Tri3`` no tiene cuadratura);
            un valor ``None`` no se pasa.
        """
        from solidum.utils.yaml_parser import _constructor_kwargs
        group_specs = group_specs or {}
        _log.info(f"Leyendo malla Gmsh desde: {self.filepath} ...")
        mesh = meshio.read(self.filepath, file_format="gmsh")

        unsupported = sorted({b.type for b in mesh.cells if b.type not in _CELL_REFERENCE
                              and b.type != "vertex"})
        if unsupported:
            _log.warning(f"  -> Celdas gmsh no admitidas, se ignoran: {unsupported}")
        dims = [_CELL_DIM[_CELL_REFERENCE[b.type][0]] for b in mesh.cells if b.type in _CELL_REFERENCE]
        if not dims:
            raise RuntimeError("Error crítico: la malla no tiene celdas admitidas "
                               f"({sorted(_CELL_REFERENCE)}).")
        mesh_dim = max(dims)

        # Grupos físicos: nodos y celdas (ids de nodo de cada una). Las celdas
        # de un grupo de contorno son los bordes o caras exactos sobre los que
        # se aplica una carga distribuida (ADR 0023).
        self.domain.physical_groups = {}
        self.domain.physical_group_cells = {}
        tag_to_name: dict[int, str] = {}
        physical = mesh.cell_data.get("gmsh:physical") if hasattr(mesh, "cell_data") else None
        if hasattr(mesh, "field_data") and physical is not None:
            tag_to_name = {int(data[0]): name for name, data in mesh.field_data.items()}
            for name in tag_to_name.values():
                self.domain.physical_groups[name] = set()
                self.domain.physical_group_cells[name] = []
            for i, block in enumerate(mesh.cells):
                for j, connectivity in enumerate(block.data):
                    tag = int(physical[i][j])
                    if tag in tag_to_name:
                        name = tag_to_name[tag]
                        self.domain.physical_group_cells[name].append(
                            tuple(int(k) + 1 for k in connectivity))
                        self.domain.physical_groups[name].update(int(k) + 1 for k in connectivity)
            for name in self.domain.physical_groups:
                self.domain.physical_groups[name] = sorted(self.domain.physical_groups[name])

        for i, pt in enumerate(mesh.points):
            coords = [pt[0], pt[1], pt[2]] if mesh_dim == 3 else [pt[0], pt[1]]
            self.domain.add_node(i + 1, [float(c) for c in coords])
        _log.info(f"  -> Nodos importados: {len(self.domain.nodes)} (malla {mesh_dim}D)")

        elem_id = 1
        for i, block in enumerate(mesh.cells):
            if block.type not in _CELL_REFERENCE:
                continue
            ref_cell, ref_coords, cell_default = _CELL_REFERENCE[block.type]
            cell_dim = _CELL_DIM[ref_cell]
            tags = physical[i] if physical is not None else None
            for j, connectivity in enumerate(block.data):
                spec, group_name = default_spec, None
                if tags is not None and int(tags[j]) in tag_to_name:
                    group_name = tag_to_name[int(tags[j])]
                    spec = group_specs.get(group_name, default_spec)
                if cell_dim < mesh_dim:
                    # Contorno, salvo que su grupo pida un elemento explícito.
                    gspec = group_specs.get(group_name, {})
                    if not (gspec.get('element') and gspec.get('explicit_element', True)):
                        continue
                name = spec.get('element') or cell_default
                if name is None:
                    continue
                cls = ElementRegistry.get(name)
                nodes = self._ordered_nodes(cls, name, block.type, ref_cell, ref_coords,
                                            connectivity, group_name)
                accepted, var_kw = _constructor_kwargs(cls)
                kwargs = dict(spec.get('kwargs', {}))
                kwargs.update({k: v for k, v in spec.get('optional', {}).items()
                               if v is not None and (var_kw or k in accepted)})
                self.domain.add_element(ElementRegistry.create(
                    name, element_id=elem_id, nodes=nodes, **kwargs))
                elem_id += 1

        if len(self.domain.elements) == 0:
            raise RuntimeError("Error crítico: No se importó ningún elemento. Revisa que tu malla "
                               "tenga superficies (2D) o volúmenes (3D).")
        if self.n_reoriented:
            _log.info(f"  -> {self.n_reoriented} elemento(s) reorientado(s) (jacobiano negativo en la malla).")
        _log.info(f"  -> Elementos importados: {len(self.domain.elements)}")
        return self.domain

    # ------------------------------------------------------------------
    def _permutations(self, cls, name, cell_type, ref_cell, ref_coords, group_name):
        """(permutación meshio → elemento, reflexión del elemento, celda de
        referencia del elemento o None) para la pareja (celda, elemento)."""
        key = (cell_type, cls)
        if key in self._perm_cache:
            return self._perm_cache[key]
        ref_nodes = getattr(cls, "REFERENCE_NODES", None)
        n_cls = len(ref_nodes) if ref_nodes is not None else getattr(cls, "N_NODES", None)
        if n_cls is not None and n_cls != len(ref_coords):
            raise ValueError(
                f"Grupo {group_name!r}: el elemento '{name}' tiene {n_cls} nodos y la celda "
                f"gmsh '{cell_type}', {len(ref_coords)}.")
        if ref_nodes is None:
            out = (list(range(len(ref_coords))), None, None)
        else:
            cls_cell = getattr(cls, "REFERENCE_CELL", None)
            what = (f"Grupo {group_name!r}: los nodos de referencia del elemento '{name}' "
                    f"(celda {cls_cell!r}) no corresponden a los de la celda gmsh '{cell_type}'.")
            if cls_cell is not None and cls_cell != ref_cell:
                raise ValueError(what)
            perm = _match(ref_coords, ref_nodes, what)
            flip = None
            if _CELL_DIM[ref_cell] >= 2:
                refl = np.array(ref_nodes, dtype=float)
                refl[:, [0, 1]] = refl[:, [1, 0]]
                flip = _match(ref_nodes, refl, what)
            out = (perm, flip, ref_cell)
        self._perm_cache[key] = out
        return out

    def _ordered_nodes(self, cls, name, cell_type, ref_cell, ref_coords, connectivity, group_name):
        perm, flip, cell = self._permutations(cls, name, cell_type, ref_cell, ref_coords, group_name)
        nodes = [self.domain.get_node(int(connectivity[k]) + 1) for k in perm]
        if flip is None or cell is None:
            return nodes
        det = _center_jacobian_det(cell, [n.coordinates for n in nodes[:_N_CORNERS[cell]]])
        if abs(det) <= 1e-14 * max(1.0, _scale(nodes)) ** _CELL_DIM[cell]:
            raise ValueError(f"Elemento '{name}' degenerado (jacobiano nulo) con los nodos "
                             f"{[n.id for n in nodes]}.")
        if det < 0.0:
            self.n_reoriented += 1
            nodes = [nodes[k] for k in flip]
        return nodes


def _scale(nodes) -> float:
    X = np.array([n.coordinates for n in nodes], dtype=float)
    return float(np.ptp(X, axis=0).max()) if len(X) else 1.0


def _center_jacobian_det(cell: str, corners) -> float:
    """Determinante del jacobiano de la geometría lineal de las esquinas en el
    centro de la celda (signo = orientación)."""
    X = np.array(corners, dtype=float)
    if cell in ("tri", "tet"):
        d = 2 if cell == "tri" else 3
        J = (X[1:d + 1, :d] - X[0, :d]).T
    else:
        d = 2 if cell == "quad" else 3
        R = np.array(_QUAD if cell == "quad" else _HEX)
        J = X[:, :d].T @ R / R.shape[0]
    return float(np.linalg.det(J))
