"""Export del dominio y resultados a formato VTK (ParaView, vía ``meshio``).

Cobertura
---------
La celda VTK de cada elemento se **deriva de su contrato**, sin listas de
clases: número de nodos y dimensión espacial (``'uz' in DOF_NAMES`` en los
mecánicos, ``FLUX_DIM`` en los térmicos, ``JUMP_DIM`` en las interfaces; la
clasificación es la de ``solidum.core.point_quantity``, ADR 0043 §3).

- 2D: 3 → ``triangle``, 4 → ``quad``, 6 → ``triangle6``, 8 → ``quad8``,
  9 → ``quad9``. Cubre Tri3, Quad4, Tri6, Quad8, Quad9 y ``Quad4Thermal``
  (y cualquier elemento de un módulo de usuario con esos nodos).
- 3D: 4 → ``tetra``, 8 → ``hexahedron``, 10 → ``tetra10``,
  20 → ``hexahedron20``, 27 → ``hexahedron27``. Cubre Tet4, Hex8, Tet10,
  Hex20, Hex27 y ``Hex8Thermal``. La numeración de nodos medios del
  proyecto coincide con la de VTK (verificado contra meshio).
- Cualquier elemento de 2 nodos (trusses, cables, frames 2D/3D) → ``line``.
  No aportan campos de esfuerzo; sí sus **fuerzas internas** de extremo
  (ADR 0002), como campos por celda ``<comp>_i`` / ``<comp>_j`` (ver
  "Fuerzas internas por celda").
- Un elemento puede forzar su celda declarando ``VTK_CELL_TYPE`` como
  atributo de clase (p. ej. una cuña de 6 nodos en 3D, que la tabla
  confundiría con ``triangle6``).
- Un elemento sin celda propia por diseño (la interfaz de adherencia, que
  une los nodos de una barra con los de su anfitrión) declara
  ``VTK_EXPORT = False`` y se omite sin aviso.

Un dominio sin ninguna celda soportada lanza ``ValueError``: meshio no
puede releer un archivo sin celdas y fallar aquí es más útil que escribir
un ``.vtu`` corrupto.

Campos por celda (sólidos mecánicos)
------------------------------------
``Sigma_XX, Sigma_YY, Sigma_ZZ, Tau_XY, Tau_YZ, Tau_XZ, Von_Mises``,
promediados sobre los puntos de Gauss **evaluados en ``U``** vía
``compute_gauss_state`` (API canónica de los sólidos, ADR 0012). El
exportador no depende de qué solver produjo ``U`` ni de si el estado quedó
comiteado: es una función pura de ``(dominio, U)``; la historia interna
(plasticidad, daño) entra por el estado committed de cada punto de Gauss,
como en cualquier reevaluación.

En 2D la componente ``σ_zz`` la aporta ``material.out_of_plane_stress``:
nula en plane stress y ``≠ 0`` en plane strain, donde el Von Mises en plano
``√(σ_xx² + σ_yy² − σ_xx·σ_yy + 3τ²)`` sería el invariante equivocado. Las
componentes ``yz``/``xz`` son cero en 2D. ``Internal_State`` es el promedio
de la ``PRIMARY_STATE_VAR`` del material sobre el estado committed.

Camino por lotes (ADR 0014 §9)
------------------------------
Los elementos que un ``Assembler`` por lotes agrupó en familias se
evalúan **por familia** (``Family.gauss_state``: un kernel compilado sobre
todos sus puntos de Gauss, ``σ_zz`` y ``Internal_State`` por arreglos),
y sólo los elementos sin familia pasan por ``compute_gauss_state`` uno a
uno. El resultado es el mismo; el coste por elemento baja de decenas de
microsegundos a uno. Los nodos, las celdas y el suavizado nodal también
se calculan por arreglos.

Campos por celda (térmicos)
---------------------------
``Flux``: vector ``q`` [W/m²] promediado sobre los puntos de Gauss,
completado con ceros a 3 componentes.

Campos por celda (interfaces, ADR 0043 §7)
------------------------------------------
Los elementos de interfaz de espesor nulo tienen celda propia (``quad``,
``quad6``, ``wedge`` o ``hexahedron`` por ``VTK_CELL_TYPE``, degenerada si
las caras coinciden) con ``Jump_<c>`` y ``Traction_<c>``: el salto y la
tracción en los ejes locales de la interfaz, ``<c>`` en ``n``, ``s`` (2D) o
``n``, ``s1``, ``s2`` (3D), media de los puntos de integración de la celda.
``Jump_n`` es la holgura normal actual (con la inicial; en el pre-proceso,
la holgura inicial). Las variables internas escalares de la ley (presión de
contacto, daño…) van por su nombre, salvo los códigos (``STATE_CODES``,
p. ej. el régimen del contacto), que se leen por punto. Los campos de
esfuerzo de sólido valen ``NaN`` en estas celdas y las interfaces no entran
en el suavizado nodal. La adherencia barra–sólido no tiene celda
(``VTK_EXPORT = False``).

Fuerzas internas por celda (estructurales 1D)
---------------------------------------------
``N_i``, ``N_j``, ``V_i``, ``M_i``… (``Vy``, ``Vz``, ``T``, ``My``, ``Mz`` en
3D): fuerzas internas de extremo de ``internal_forces`` (ADR 0002) en ejes
locales, nodo i y nodo j, convención de Reglas.md §5. ``NaN`` en las celdas
sin esa componente. Con carga distribuida, pase
``element_forces=SolveResult.element_forces`` para que descuenten la carga
equivalente. ``Internal_State`` de un marco con sección por fibras promedia la
``PRIMARY_STATE_VAR`` sobre fibras y puntos de Gauss; las fibras de
compensación de las barras con ``displaces`` no entran
(``Section.material_point_states``).

Campos nodales
--------------
- ``Displacements``: vector 3D ``(ux, uy, uz)``; ``uz = 0`` en modelos 2D.
- ``Rotations``: ``(rx, ry, rz)``, sólo si algún nodo tiene DOF rotacional.
- ``Temperature``: escalar, sólo si algún nodo tiene el DOF ``T``.
- ``Supports`` / ``Supports_Rot``: 1.0 en cada componente con Dirichlet.
- ``External_Forces``: fuerzas nodales aplicadas ``(fx, fy, fz)``.
- ``Sigma_XX_nodal`` … ``Von_Mises_nodal``: suavizado nodal de orden 0
  (media simple sobre las celdas adyacentes) de los campos anteriores.
  Sólo promedia las celdas de **sólidos mecánicos**: ni las líneas ni los
  elementos térmicos ni las interfaces, tampoco cuando un térmico se superpone a un sólido con
  los mismos nodos (dominio multicampo ``(ux, uy, T)``, ADR 0016). Un modelo
  sin ningún sólido mecánico (sólo térmico, sólo estructural 1D) no exporta
  estos campos.

Restricciones multipunto (ADR 0004)
-----------------------------------
El exportador recibe el ``U`` completo, con los DOF esclavos ya expandidos
(``u = T·u_libre + g``), así que ``Displacements``, ``Temperature``, los
campos por celda (``Sigma_*``, ``Flux``…) y las fuerzas internas son los del
modelo resuelto también en los nodos esclavos. **Limitación del suavizado
nodal**: dos nodos coincidentes atados por una MPC (``u_s = u_m``, p. ej.
para unir mallas no conformes o pegar dos partes) siguen siendo nodos
distintos, y cada uno promedia sólo las celdas de su lado. En la interfaz
``Sigma_*_nodal`` del maestro y del esclavo difieren entre sí y del valor
del modelo con el nodo fundido (que, con el mismo número de celdas a cada
lado, es su media). Los valores por celda no se ven afectados; para leer el
esfuerzo en la interfaz, use los campos por celda o promedie a mano las dos
copias.
"""
from __future__ import annotations

import numpy as np

from solidum.core.dof_names import ROTATIONS as _ROTATIONAL_DOFS
from solidum.core.dof_names import TEMPERATURE as _TEMPERATURE_DOF
from solidum.core.dof_names import TRANSLATIONS as _TRANSLATIONAL_DOFS
from solidum.core.domain import Domain
from solidum.core.point_quantity import FLUX, JUMP, STRAIN, point_quantity
from solidum.core.section import Section
from solidum.logging import get_logger
from solidum.math.batch.postprocess import group_by_family

try:
    import meshio
except ImportError:
    meshio = None

_log = get_logger("exporters.vtk")

# Celdas que la tabla de tipos VTK de meshio conoce y sus tablas de dimensión
# topológica y de número de nodos no (meshio 5.3.5): sin ellas, meshio no
# crea el bloque al escribir ni al leer. VTK 30 (cuadrilátero
# cuadrático–lineal: 4 vértices y los intermedios de los lados 1-2 y 3-4) es
# la celda de Interface2D6. Celda → (dimensión, nodos).
_MESHIO_MISSING_CELLS = {"quad6": (2, 6)}


def ensure_meshio_cells() -> None:
    """Completa las tablas de meshio con las celdas que usa Solidum y le
    faltan (:data:`_MESHIO_MISSING_CELLS`). Idempotente; no cambia una
    entrada que meshio ya tenga. La llaman el exportador al importarse y el
    lector de las figuras."""
    if meshio is None:
        return
    try:
        from meshio._common import num_nodes_per_cell
        from meshio._mesh import topological_dimension
    except ImportError:         # otra versión de meshio, con otra organización
        return
    for cell, (dim, n_nodes) in _MESHIO_MISSING_CELLS.items():
        topological_dimension.setdefault(cell, dim)
        num_nodes_per_cell.setdefault(cell, n_nodes)


def _check_meshio_cells(cell_types) -> None:
    """Error claro si una celda de :data:`_MESHIO_MISSING_CELLS` no quedó en
    las tablas de meshio (otra versión con otra organización interna): sin
    esto, meshio fallaría con un ``KeyError`` sin explicación."""
    needed = [c for c in cell_types if c in _MESHIO_MISSING_CELLS]
    if not needed:
        return
    try:
        from meshio._common import num_nodes_per_cell
        from meshio._mesh import topological_dimension
        ok = all(c in num_nodes_per_cell and c in topological_dimension for c in needed)
    except ImportError:
        ok = False
    if not ok:
        raise RuntimeError(
            f"VtkExporter: esta versión de meshio ({getattr(meshio, '__version__', '?')}) no sabe "
            f"escribir las celdas {needed} (Interface2D6: celda VTK 30) y Solidum no pudo "
            f"completarle las tablas (solidum.utils.vtk_exporter.ensure_meshio_cells).")


ensure_meshio_cells()

# Celda VTK por (dimensión espacial, número de nodos). Los elementos de dos
# nodos se resuelven aparte como ``line``.
_CELL_2D = {3: "triangle", 4: "quad", 6: "triangle6", 8: "quad8", 9: "quad9"}
_CELL_3D = {4: "tetra", 8: "hexahedron", 10: "tetra10",
            20: "hexahedron20", 27: "hexahedron27"}

# Nombres y orden de las seis componentes de esfuerzo exportadas.
STRESS_COMPONENT_NAMES = ("Sigma_XX", "Sigma_YY", "Sigma_ZZ",
                          "Tau_XY", "Tau_YZ", "Tau_XZ")
#: Campos nodales que puede escribir el exportador (los de un modelo concreto
#: son un subconjunto: sin giros no hay ``Rotations``...).
NODAL_FIELD_NAMES = ("Displacements", "Rotations", "Temperature", "Supports", "Supports_Rot",
                     "External_Forces", "Nodal_Mass", "Nodal_Inertia",
                     *(f"{c}_nodal" for c in STRESS_COMPONENT_NAMES), "Von_Mises_nodal")


# ---------------------------------------------------------------------------
# Clasificación de elementos
# ---------------------------------------------------------------------------

def _is_thermal(elem) -> bool:
    return point_quantity(elem).kind == FLUX


def _spatial_dim(elem) -> int:
    q = point_quantity(elem)
    if q.kind in (FLUX, JUMP):
        return int(q.dim)
    return 3 if "uz" in elem.DOF_NAMES else 2


def cell_type_for(elem) -> str | None:
    """Celda VTK de un elemento, o ``None`` si no hay una soportada."""
    forced = getattr(elem, "VTK_CELL_TYPE", None)
    if forced:
        return str(forced)
    n = len(elem.nodes)
    if n == 2:
        return "line"
    table = _CELL_3D if _spatial_dim(elem) == 3 else _CELL_2D
    return table.get(n)


def _is_mechanical_solid(elem) -> bool:
    return point_quantity(elem).kind == STRAIN and len(elem.nodes) > 2 \
        and hasattr(elem, "compute_gauss_state")


def _is_interface_cell(elem) -> bool:
    """Elemento de salto local con su registro por punto (ADR 0043 §2):
    ``JUMP`` implica ``LocalJumpElement`` (lo exige el barrido de elementos)."""
    return point_quantity(elem).kind == JUMP


# ---------------------------------------------------------------------------
# Esfuerzos por elemento
# ---------------------------------------------------------------------------

def von_mises(s: np.ndarray) -> float:
    """Esfuerzo equivalente de Von Mises de ``[σxx, σyy, σzz, τxy, τyz, τxz]``."""
    sxx, syy, szz, txy, tyz, txz = (float(v) for v in s)
    return float(np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (txy * txy + tyz * tyz + txz * txz)
    ))


def von_mises_rows(s: np.ndarray) -> np.ndarray:
    """:func:`von_mises` fila a fila sobre ``(n, 6)`` → ``(n,)``."""
    s = np.asarray(s, dtype=float).reshape(-1, 6)
    sxx, syy, szz = s[:, 0], s[:, 1], s[:, 2]
    txy, tyz, txz = s[:, 3], s[:, 4], s[:, 5]
    return np.sqrt(
        0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
        + 3.0 * (txy * txy + tyz * tyz + txz * txz)
    )


def _six_components(avg: np.ndarray, szz: np.ndarray | None) -> np.ndarray:
    """Completa a ``(n, 6)`` un promedio ``(n, 3)`` (Voigt 2D + σ_zz) o
    devuelve el ``(n, 6)`` tal cual."""
    if avg.shape[1] == 6:
        return avg
    if avg.shape[1] == 3:
        out = np.zeros((avg.shape[0], 6))
        out[:, 0] = avg[:, 0]
        out[:, 1] = avg[:, 1]
        out[:, 2] = 0.0 if szz is None else szz
        out[:, 3] = avg[:, 2]
        return out
    raise ValueError(
        f"VtkExporter: esfuerzos con {avg.shape[1]} componentes; se esperaban "
        f"3 (Voigt 2D) ó 6 (Voigt 3D)."
    )


def element_average_stress(elem, U: np.ndarray | None) -> np.ndarray:
    """Esfuerzo promedio del elemento, shape ``(6,)`` en el orden de
    :data:`STRESS_COMPONENT_NAMES`. Camino por elemento.

    Evalúa ``compute_gauss_state(U)`` y promedia sobre los puntos de Gauss.
    En 2D completa ``σ_zz`` con ``material.out_of_plane_stress`` punto a
    punto (usa el estado committed de cada punto, como la propia
    reevaluación). Ceros si ``U`` es ``None`` o el elemento no es un sólido
    mecánico.
    """
    out = np.zeros(6)
    if U is None or not _is_mechanical_solid(elem):
        return out
    stresses = np.asarray(elem.compute_gauss_state(U)["stress"], dtype=float)
    if stresses.ndim != 2 or stresses.shape[0] == 0:
        return out
    avg = stresses.mean(axis=0)
    if avg.size == 6:
        return avg
    if avg.size == 3:
        state_vars = getattr(getattr(elem, "state", None), "vars", None)
        szz = 0.0
        for idx, sig in enumerate(stresses):
            sv = state_vars[idx] if state_vars is not None and idx < len(state_vars) else None
            szz += float(elem.material.out_of_plane_stress(sig, sv))
        szz /= stresses.shape[0]
        out[0], out[1], out[2], out[3] = avg[0], avg[1], szz, avg[2]
        return out
    raise ValueError(
        f"VtkExporter: elemento id={elem.id} devuelve esfuerzos con "
        f"{avg.size} componentes; se esperaban 3 (Voigt 2D) ó 6 (Voigt 3D)."
    )


def _element_average_flux(elem, U: np.ndarray | None) -> np.ndarray:
    out = np.zeros(3)
    if U is None:
        return out
    flux = np.asarray(elem.compute_gauss_state(U)["flux"], dtype=float)
    if flux.ndim == 2 and flux.shape[0] > 0:
        avg = flux.mean(axis=0)
        out[: avg.size] = avg
    return out


def _extract_state_scalar(sv: dict, primary_key=None) -> float:
    """Escalar representativo del dict de variables de estado de un punto.

    Prioriza ``primary_key`` (un nombre o una secuencia de nombres, en orden
    de preferencia); si no está, el primer escalar numérico.
    """
    if not isinstance(sv, dict):
        return 0.0
    keys = (primary_key,) if isinstance(primary_key, str) else tuple(primary_key or ())
    for key in keys:
        if key and key in sv:
            v = sv[key]
            if isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool):
                return float(v)
    for v in sv.values():
        if isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool):
            return float(v)
    return 0.0


def _flatten_states(entries):
    """Estados de punto material de ``state.vars``: los de un sólido son un
    dict por punto de Gauss; los de un marco con sección, una lista (un
    estado por fibra) por punto de Gauss."""
    for sv in entries:
        if isinstance(sv, (list, tuple)):
            yield from _flatten_states(sv)
        elif sv is not None:
            yield sv


def _state_var_avg(elem) -> float:
    """Promedio de la ``PRIMARY_STATE_VAR`` de los materiales del elemento
    sobre el estado committed de todos sus puntos materiales (puntos de
    Gauss; en un marco con sección, además, sus fibras)."""
    state_vars_list = getattr(getattr(elem, "state", None), "vars", None)
    if state_vars_list is None:
        return 0.0
    # Con sección, sólo los puntos materiales físicos: sin los estados
    # auxiliares (fibras de compensación de una FiberSection con displaces).
    section = getattr(elem, "section", None)
    if isinstance(section, Section):
        state_vars_list = [section.material_point_states(sv) for sv in state_vars_list]
    materials = elem.materials() if hasattr(elem, "materials") else []
    keys = [k for k in (getattr(m, "PRIMARY_STATE_VAR", None) for m in materials) if k]
    vals = [_extract_state_scalar(sv, keys) for sv in _flatten_states(state_vars_list)]
    if not vals:
        return 0.0
    return sum(vals) / len(vals)


def _is_number(v) -> bool:
    return isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool)


def _material_point_states(elem, trial: bool = False) -> list:
    """Estados (dicts) de los puntos materiales del elemento, committed o
    trial: uno por punto de Gauss en un sólido; en un marco con sección, sus
    fibras físicas (sin las de compensación)."""
    st = getattr(elem, "state", None)
    entries = getattr(st, "vars_trial" if trial else "vars", None)
    if entries is None:
        return []
    section = getattr(elem, "section", None)
    if isinstance(section, Section):
        entries = [section.material_point_states(sv) for sv in entries]
    return list(_flatten_states(entries))


def _unaveraged(material) -> set:
    """Variables internas escalares de ``material`` que el VTK por celda no
    promedia: ángulos (``STATE_ANGLES``) y códigos (``STATE_CODES``)."""
    return set(getattr(material, "STATE_ANGLES", ())) | set(getattr(material, "STATE_CODES", ()))


def _angle_names(elem) -> set:
    """Variables internas que el VTK por celda no promedia (ángulos y
    códigos) en los materiales del elemento."""
    materials = elem.materials() if hasattr(elem, "materials") else []
    return {n for m in materials for n in _unaveraged(m)}


def state_variable_names(elements) -> tuple[list[str], set]:
    """Variables internas escalares de los materiales de ``elements``, en el
    orden en que aparecen, y cuáles de ellas no se promedian por celda: los
    ángulos (``STATE_ANGLES``) y los códigos (``STATE_CODES``).

    Las escalares del ``STATE_SCHEMA`` de cada material (válidas antes del
    primer paso, cuando el estado de muchos materiales aún no existe) y,
    para los materiales sin esquema (interfaces), las que ya tenga el
    estado del elemento."""
    names: dict[str, None] = {}
    angles: set = set()
    seen: set = set()
    for elem in elements:
        if getattr(elem, "VTK_EXPORT", True) is False:
            continue
        materials = elem.materials() if hasattr(elem, "materials") else []
        key = (type(elem), tuple(id(m) for m in materials))
        if key in seen:
            continue
        seen.add(key)
        for m in materials:
            for k, shape in (getattr(m, "STATE_SCHEMA", None) or {}).items():
                if tuple(shape) == ():
                    names.setdefault(k, None)
        for sv in _material_point_states(elem):
            for k, v in sv.items():
                if _is_number(v):
                    names.setdefault(k, None)
        angles |= _angle_names(elem)
    return list(names), angles


def _state_cell_fields(elements, slot, conn_by_type: dict, groups: dict, loose: list,
                       names: list[str], angles: set) -> dict:
    """Promedio sobre los puntos materiales de cada variable interna escalar,
    por celda: ``{nombre: {tipo_de_celda: (n_celdas,)}}``. ``NaN`` en las
    celdas cuyo material no tiene esa variable; los ángulos y los códigos no
    se promedian (no figuran)."""
    wanted = [n for n in names if n not in angles]
    if not wanted:
        return {}
    fields = {n: {ct: np.full(len(conn), np.nan) for ct, conn in conn_by_type.items()}
              for n in wanted}
    for fam, positions in groups.items():
        positions = [p for p in positions if slot[p] is not None]
        if not positions or fam.is_thermal:
            continue
        rows = {id(e): r for r, e in enumerate(fam.elements)}
        fam_rows = np.array([rows[id(elements[p])] for p in positions], dtype=np.int64)
        cell_type = slot[positions[0]][0]
        idx = np.array([slot[p][1] for p in positions], dtype=np.int64)
        fam_angles = _unaveraged(fam.material)
        schema = fam.schema
        S = fam.state.S_committed.reshape(fam.n_elements, fam.n_gp, -1)
        for k, (name, shape) in enumerate(zip(schema.names, schema.shapes)):
            if shape != () or name not in fields or name in fam_angles:
                continue
            fields[name][cell_type][idx] = S[fam_rows, :, schema.offsets[k]].mean(axis=1)
    for p in loose:
        if slot[p] is None:
            continue
        elem = elements[p]
        states = _material_point_states(elem)
        if not states:
            continue
        cell_type, idx = slot[p]
        skip = _angle_names(elem)
        for name in wanted:
            if name in skip:
                continue
            vals = [float(sv[name]) for sv in states if name in sv and _is_number(sv[name])]
            if vals:
                fields[name][cell_type][idx] = sum(vals) / len(vals)
    return fields


def _internal_force_fields(elements, slot, conn_by_type: dict,
                           U: np.ndarray | None,
                           element_forces: dict | None) -> dict:
    """Fuerzas internas de extremo de los elementos estructurales 1D como
    campos por celda: ``{nombre: {tipo_de_celda: (n_celdas,)}}``.

    Por cada componente de :class:`~solidum.core.element_forces.ElementForces`
    (``N``, ``V``, ``M`` en los marcos 2D; ``N``, ``Vy``, ``Vz``, ``T``,
    ``My``, ``Mz`` en los 3D; ``N`` en barras y cables) hay dos campos,
    ``<comp>_i`` y ``<comp>_j``: el valor en el nodo i y en el j, en ejes
    locales y con la convención de Reglas.md §5. Las celdas que no tienen
    esa componente (sólidos, una barra en el campo ``M_i``…) valen ``NaN``,
    no cero: "no aplica" no es "vale cero".

    Fuente: ``element_forces`` (p. ej. ``SolveResult.element_forces``, que ya
    descuenta la carga distribuida) si se pasa; si no, ``internal_forces(U)``
    de cada elemento que lo implementa, sin carga equivalente (fuerzas
    nodales ``K·u``: exactas sin carga distribuida sobre el elemento).
    """
    from solidum.core.element import Element

    if element_forces is None and U is None:
        return {}
    fields: dict[str, dict[str, np.ndarray]] = {}
    for p, elem in enumerate(elements):
        if slot[p] is None:
            continue
        if element_forces is not None:
            forces = element_forces.get(elem.id)
        elif type(elem).internal_forces is Element.internal_forces:
            forces = None
        else:
            forces = elem.internal_forces(U)
        if forces is None:
            continue
        cell_type, idx = slot[p]
        for comp, values in forces.components.items():
            values = np.asarray(values, dtype=float).reshape(-1)
            for k, end in enumerate(("i", "j")):
                name = f"{comp}_{end}"
                per_type = fields.setdefault(name, {})
                if cell_type not in per_type:
                    per_type[cell_type] = np.full(len(conn_by_type[cell_type]), np.nan)
                per_type[cell_type][idx] = values[k]
    return fields


def _interface_cell_fields(elements, slot, conn_by_type: dict, U: np.ndarray) -> dict:
    """Salto y tracción locales de los elementos de interfaz como campos por
    celda (ADR 0043 §7): ``{nombre: {tipo_de_celda: (n_celdas,)}}`` con los
    nombres ``Jump_<c>`` y ``Traction_<c>`` y ``<c>`` las componentes del
    elemento (``n``, ``s``; ``n``, ``s1``, ``s2``). Valor de la celda: la
    media de sus puntos de integración, del registro por punto
    (``interface_points``: salto de ``U``, tracción de la ley con el estado
    confirmado). ``NaN`` en las celdas que no son de interfaz."""
    jump_fields: dict[str, dict[str, np.ndarray]] = {}
    traction_fields: dict[str, dict[str, np.ndarray]] = {}
    for p, elem in enumerate(elements):
        if slot[p] is None or not _is_interface_cell(elem):
            continue
        cell_type, idx = slot[p]
        pts = elem.interface_points(U)
        jump = np.mean([q.jump for q in pts], axis=0)
        traction = np.mean([q.traction for q in pts], axis=0)
        for k, comp in enumerate(elem.jump_components()):
            for fields, prefix, v in ((jump_fields, "Jump", jump),
                                      (traction_fields, "Traction", traction)):
                per_type = fields.setdefault(f"{prefix}_{comp}", {})
                if cell_type not in per_type:
                    per_type[cell_type] = np.full(len(conn_by_type[cell_type]), np.nan)
                per_type[cell_type][idx] = v[k]
    return {**jump_fields, **traction_fields}


def _family_cell_fields(fam, U: np.ndarray | None):
    """Campos por celda de toda una familia, por arreglos: esfuerzo medio
    ``(N, 6)``, flujo medio ``(N, 3)`` y ``Internal_State`` ``(N,)``, en
    el orden de ``fam.elements``. Con ``U = None`` los dos primeros son
    ceros (como en el camino por elemento)."""
    N = fam.n_elements
    stress = np.zeros((N, 6))
    flux = np.zeros((N, 3))
    primary = getattr(fam.material, "PRIMARY_STATE_VAR", None)
    state = np.zeros(N) if fam.is_thermal else fam.state_average(primary)
    if U is None:
        return stress, flux, state
    gs = fam.gauss_state(U)
    if fam.is_thermal:
        avg = gs["flux"].mean(axis=1)
        flux[:, : avg.shape[1]] = avg
        return stress, flux, state
    sig = gs["stress"]                       # (N, n_gp, n_sig)
    avg = sig.mean(axis=1)
    szz = None
    if sig.shape[2] == 3:
        P = N * fam.n_gp
        rows = fam.material.batch_out_of_plane_stress(
            sig.reshape(P, 3), fam.state.S_committed)
        szz = np.asarray(rows, dtype=float).reshape(N, fam.n_gp).mean(axis=1)
    stress[:] = _six_components(avg, szz)
    return stress, flux, state


# ---------------------------------------------------------------------------
# Exportador
# ---------------------------------------------------------------------------

class VtkExporter:
    """Exporta el dominio y los resultados a formato VTK para ParaView."""

    def __init__(self, domain: Domain):
        self.domain = domain

    def export(self, filepath: str, U: np.ndarray | None = None,
               F_ext: np.ndarray | None = None,
               nodal_vars: list | None = None, elem_vars: list | None = None,
               element_forces: dict | None = None):
        """Escribe ``filepath`` (``.vtu``, ``.vtk``… según extensión).

        Parameters
        ----------
        U, F_ext
            Campo global de la solución y de fuerzas aplicadas. ``None``
            exporta la malla sin deformar con campos nulos (pre-proceso).
        nodal_vars, elem_vars
            Filtro opcional de nombres de campos nodales / por celda. ``None``
            exporta todos los disponibles.
        element_forces
            ``{id de elemento: ElementForces}`` de los elementos
            estructurales 1D, p. ej. ``SolveResult.element_forces``, que
            descuenta la carga distribuida (peso propio, fuerza de cuerpo).
            ``None`` ⇒ se evalúa ``internal_forces(U)`` de cada elemento,
            exacto cuando no hay carga distribuida sobre él.
        """
        if meshio is None:
            _log.warning("La librería 'meshio' no está instalada. Omitiendo exportación VTK.")
            return

        ndof = self.domain.total_dofs
        U_eval = None if U is None else np.asarray(U, dtype=float)
        if U is None:
            U = np.zeros(ndof)
        if F_ext is None:
            F_ext = np.zeros(ndof)

        def want_nodal(name):
            return nodal_vars is None or name in nodal_vars

        def want_cell(name):
            return elem_vars is None or name in elem_vars

        # --- Puntos y campos nodales directos ---------------------------------
        node_ids = list(self.domain.nodes.keys())
        node_map = {nid: i for i, nid in enumerate(node_ids)}
        n_nodes = len(node_ids)

        points = np.zeros((n_nodes, 3))
        displacements = np.zeros((n_nodes, 3))
        rotations = np.zeros((n_nodes, 3))
        temperature = np.zeros(n_nodes)
        supports = np.zeros((n_nodes, 3))
        supports_rot = np.zeros((n_nodes, 3))
        forces = np.zeros((n_nodes, 3))
        any_rotational = False
        any_temperature = False

        for nid, node in self.domain.nodes.items():
            idx = node_map[nid]
            coords = list(node.coordinates)[:3]
            points[idx, :len(coords)] = coords

            for axis, dof_name in enumerate(_TRANSLATIONAL_DOFS):
                if dof_name in node.dofs and node.dofs[dof_name] < len(U):
                    displacements[idx, axis] = U[node.dofs[dof_name]]
                    if node.dofs[dof_name] < len(F_ext):
                        forces[idx, axis] = F_ext[node.dofs[dof_name]]
                if dof_name in node.boundary_conditions:
                    supports[idx, axis] = 1.0

            for axis, dof_name in enumerate(_ROTATIONAL_DOFS):
                if dof_name in node.dofs:
                    any_rotational = True
                    if node.dofs[dof_name] < len(U):
                        rotations[idx, axis] = U[node.dofs[dof_name]]
                    if dof_name in node.boundary_conditions:
                        supports_rot[idx, axis] = 1.0

            if _TEMPERATURE_DOF in node.dofs:
                any_temperature = True
                if node.dofs[_TEMPERATURE_DOF] < len(U):
                    temperature[idx] = U[node.dofs[_TEMPERATURE_DOF]]

        # Masas de nudo (ADR 0024): masa por traslación e inercia por giro.
        nodal_masses = getattr(self.domain, "nodal_masses", None) or {}
        node_mass = np.zeros((n_nodes, 3))
        node_inertia = np.zeros((n_nodes, 3))
        for (nid, dof_name), value in nodal_masses.items():
            if nid not in node_map:
                continue
            if dof_name in _TRANSLATIONAL_DOFS:
                node_mass[node_map[nid], _TRANSLATIONAL_DOFS.index(dof_name)] = value
            elif dof_name in _ROTATIONAL_DOFS:
                node_inertia[node_map[nid], _ROTATIONAL_DOFS.index(dof_name)] = value

        # --- Celdas agrupadas por tipo -----------------------------------------
        elements = list(self.domain.elements.values())
        conn_by_type: dict[str, list[list[int]]] = {}
        # Por celda, si es un sólido mecánico: sólo esas celdas entran en el
        # suavizado nodal de σ (un térmico superpuesto comparte tipo de celda).
        solid_by_type: dict[str, list[bool]] = {}
        slot: list[tuple[str, int] | None] = []      # (tipo, índice en su bloque)
        for elem in elements:
            if getattr(elem, "VTK_EXPORT", True) is False:
                # Elementos sin celda propia por diseño (p. ej. la interfaz de
                # adherencia): se omiten sin aviso.
                slot.append(None)
                continue
            cell_type = cell_type_for(elem)
            if cell_type is None:
                _log.warning(
                    f"VtkExporter: elemento id={getattr(elem, 'id', '?')} de tipo "
                    f"{type(elem).__name__} con {len(elem.nodes)} nodos no tiene "
                    f"celda VTK soportada; se omite. Declare VTK_CELL_TYPE en la "
                    f"clase si corresponde a una celda válida de VTK."
                )
                slot.append(None)
                continue
            block = conn_by_type.setdefault(cell_type, [])
            block.append([node_map[n.id] for n in elem.nodes])
            solid_by_type.setdefault(cell_type, []).append(_is_mechanical_solid(elem))
            slot.append((cell_type, len(block) - 1))

        if not conn_by_type:
            raise ValueError(
                "VtkExporter: ningún elemento del dominio tiene celda VTK "
                "soportada; no se escribe el archivo (meshio no puede releer "
                "una malla sin celdas)."
            )
        _check_meshio_cells(conn_by_type)

        state_by_type = {ct: np.zeros(len(c)) for ct, c in conn_by_type.items()}
        stress_by_type = {ct: np.zeros((len(c), 6)) for ct, c in conn_by_type.items()}
        flux_by_type = {ct: np.zeros((len(c), 3)) for ct, c in conn_by_type.items()}

        # --- Campos por celda: familias por arreglos, el resto por elemento ---
        groups, loose = group_by_family(elements)
        for fam, positions in groups.items():
            positions = [p for p in positions if slot[p] is not None]
            if not positions:
                continue
            stress_f, flux_f, state_f = _family_cell_fields(fam, U_eval)
            rows = {id(e): r for r, e in enumerate(fam.elements)}
            fam_rows = np.array([rows[id(elements[p])] for p in positions], dtype=np.int64)
            cell_type = slot[positions[0]][0]
            idx = np.array([slot[p][1] for p in positions], dtype=np.int64)
            stress_by_type[cell_type][idx] = stress_f[fam_rows]
            flux_by_type[cell_type][idx] = flux_f[fam_rows]
            state_by_type[cell_type][idx] = state_f[fam_rows]
        for p in loose:
            if slot[p] is None:
                continue
            elem = elements[p]
            cell_type, idx = slot[p]
            state_by_type[cell_type][idx] = _state_var_avg(elem)
            if _is_thermal(elem):
                flux_by_type[cell_type][idx] = _element_average_flux(elem, U_eval)
            elif _is_interface_cell(elem):
                # Una interfaz no tiene esfuerzo de sólido: «no aplica», no
                # «vale cero» (comparte el tipo de celda con los sólidos).
                stress_by_type[cell_type][idx] = np.nan
            else:
                stress_by_type[cell_type][idx] = element_average_stress(elem, U_eval)

        cells = []
        state_arrays, flux_arrays = [], []
        stress_arrays = {name: [] for name in STRESS_COMPONENT_NAMES}
        vm_arrays = []
        for cell_type, conn in conn_by_type.items():
            cells.append((cell_type, np.array(conn, dtype=int)))
            state_arrays.append(state_by_type[cell_type])
            s = stress_by_type[cell_type]
            for k, name in enumerate(STRESS_COMPONENT_NAMES):
                stress_arrays[name].append(s[:, k])
            vm_arrays.append(von_mises_rows(s))
            flux_arrays.append(flux_by_type[cell_type])

        # --- Campos nodales -------------------------------------------------
        point_data = {}
        if want_nodal("Displacements"):
            point_data["Displacements"] = displacements
        if any_rotational and want_nodal("Rotations"):
            point_data["Rotations"] = rotations
        if any_temperature and want_nodal("Temperature"):
            point_data["Temperature"] = temperature
        if want_nodal("Supports"):
            point_data["Supports"] = supports
        if any_rotational and want_nodal("Supports_Rot"):
            point_data["Supports_Rot"] = supports_rot
        if want_nodal("External_Forces"):
            point_data["External_Forces"] = forces
        if nodal_masses and want_nodal("Nodal_Mass"):
            point_data["Nodal_Mass"] = node_mass
        if nodal_masses and any_rotational and want_nodal("Nodal_Inertia"):
            point_data["Nodal_Inertia"] = node_inertia

        # Suavizado nodal de orden 0: media simple de los sólidos mecánicos
        # adyacentes.
        nodal = _smooth_stress_to_nodes(n_nodes, conn_by_type, stress_by_type,
                                        solid_by_type)
        if nodal is not None:
            for k, name in enumerate(STRESS_COMPONENT_NAMES):
                if want_nodal(f"{name}_nodal"):
                    point_data[f"{name}_nodal"] = nodal[:, k]
            if want_nodal("Von_Mises_nodal"):
                point_data["Von_Mises_nodal"] = von_mises_rows(nodal)

        # --- Campos por celda -----------------------------------------------
        cell_data = {}
        if want_cell("Internal_State"):
            cell_data["Internal_State"] = state_arrays
        for name in STRESS_COMPONENT_NAMES:
            if want_cell(name):
                cell_data[name] = stress_arrays[name]
        if want_cell("Von_Mises"):
            cell_data["Von_Mises"] = vm_arrays
        if any_temperature and want_cell("Flux"):
            cell_data["Flux"] = flux_arrays
        force_fields = _internal_force_fields(elements, slot, conn_by_type,
                                              U_eval, element_forces)
        # Salto y tracción de las interfaces; sin U (pre-proceso), los de
        # U = 0: el salto normal es la holgura inicial.
        interface_fields = _interface_cell_fields(elements, slot, conn_by_type, U)
        for name, per_type in {**force_fields, **interface_fields}.items():
            if want_cell(name):
                cell_data[name] = [
                    per_type.get(ct, np.full(len(conn), np.nan))
                    for ct, conn in conn_by_type.items()
                ]
        # Variables internas escalares de los materiales, por su nombre.
        state_names, angle_names = state_variable_names(elements)
        if any(want_cell(n) for n in state_names):
            state_fields = _state_cell_fields(elements, slot, conn_by_type, groups, loose,
                                              state_names, angle_names)
            for name, per_type in state_fields.items():
                if want_cell(name) and name not in cell_data:
                    cell_data[name] = [per_type[ct] for ct in conn_by_type]
        if not cell_data:
            cell_data = None

        mesh = meshio.Mesh(points, cells, point_data=point_data, cell_data=cell_data)
        mesh.write(filepath)

    def export_integration_points(self, filepath: str, U: np.ndarray,
                                  trial: bool = False) -> int:
        """Escribe los valores en cada punto de integración de los sólidos
        mecánicos y de las interfaces como una nube de puntos (celdas
        ``vertex``).

        Por punto: ``Displacements`` (``N(ξ)·u_e``, para deformar la nube en
        ParaView; en una interfaz, el de la cara A, donde está el punto),
        ``Element_Id``, ``Point`` (número del punto en su elemento, desde 1),
        las seis componentes de esfuerzo y ``Von_Mises`` (sólidos), el salto
        y la tracción locales ``Jump_<c>`` y ``Traction_<c>`` (interfaces,
        ADR 0043 §7) y **todas** las variables internas escalares del
        material, también los ángulos (``STATE_ANGLES``) y los códigos
        (``STATE_CODES``), con su nombre; ``NaN`` donde el elemento no tiene
        esa magnitud. Las posiciones son las de la geometría de referencia.
        Los valores son los del estado con que el análisis equilibró el paso:
        el confirmado (``trial=True``: el de prueba, en un análisis lineal que
        no confirma estado), los mismos que lee la historia (ADR 0038); el
        salto sale de ``U``.

        Devuelve el número de puntos escritos (0: no escribe el archivo).
        """
        if meshio is None:
            _log.warning("La librería 'meshio' no está instalada. Omitiendo exportación VTK.")
            return 0
        U = np.asarray(U, dtype=float)
        exported = [e for e in self.domain.elements.values()
                    if getattr(e, "VTK_EXPORT", True) is not False]
        elements = [e for e in exported if _is_mechanical_solid(e)]
        interfaces = [e for e in exported if _is_interface_cell(e)]
        names, _ = state_variable_names(elements + interfaces)
        blocks = []   # (posiciones, desplazamientos, σ (n, 6), ids, punto, {var: (n,)}, {campo: (n,)})
        groups, loose = group_by_family(elements)
        for fam, positions in groups.items():
            if fam.is_thermal:
                continue
            blk = _family_points(fam, [elements[p] for p in positions], U, trial, names)
            if blk is not None:
                blocks.append((*blk, {}))
        for p in loose:
            blk = _element_points(elements[p], U, trial, names)
            if blk is not None:
                blocks.append((*blk, {}))
        blocks.extend(_interface_points(el, U, trial, names) for el in interfaces)
        if not blocks:
            return 0
        pos = np.vstack([_pad_columns(b[0]) for b in blocks])
        n = pos.shape[0]
        sig = np.vstack([b[2] for b in blocks])
        point_data = {"Displacements": np.vstack([_pad_columns(b[1]) for b in blocks]),
                      "Element_Id": np.concatenate([b[3] for b in blocks]).astype(float),
                      "Point": np.concatenate([b[4] for b in blocks]).astype(float)}
        for k, comp in enumerate(STRESS_COMPONENT_NAMES):
            point_data[comp] = sig[:, k]
        point_data["Von_Mises"] = von_mises_rows(sig)
        extra = list(dict.fromkeys(f for b in blocks for f in b[6]))
        for field in extra:
            point_data[field] = np.concatenate([b[6].get(field, np.full(len(b[3]), np.nan))
                                                for b in blocks])
        for name in names:
            point_data.setdefault(name, np.concatenate([b[5][name] for b in blocks]))
        cells = [("vertex", np.arange(n, dtype=int).reshape(-1, 1))]
        meshio.Mesh(pos, cells, point_data=point_data).write(filepath)
        return n


def _pad_columns(a: np.ndarray) -> np.ndarray:
    """``(n, d)`` → ``(n, 3)`` con ceros."""
    a = np.asarray(a, dtype=float)
    out = np.zeros((a.shape[0], 3))
    out[:, :a.shape[1]] = a
    return out


def _interface_points(elem, U: np.ndarray, trial: bool, names: list[str]):
    """Bloque de la nube de puntos de un elemento de interfaz: salto de
    ``U``, tracción y variables de la ley del estado confirmado (o de
    prueba), sin esfuerzo de sólido (``NaN``)."""
    u_e = np.asarray(elem.get_local_displacements(U), dtype=float)
    X = np.asarray(elem.positions, dtype=float)
    n_gp = X.shape[0]
    side_a = getattr(elem, "side_a_displacements", None)
    disp = side_a(u_e) if side_a is not None else np.zeros_like(X)
    st = elem.state
    jumps = np.asarray(elem.jumps(u_e), dtype=float)
    tractions = np.array([np.atleast_1d(np.asarray(t, dtype=float))
                          for t in (st.stresses_trial if trial else st.stresses)])
    states = st.vars_trial if trial else st.vars
    fields = {}
    comps = elem.jump_components()
    for prefix, v in (("Jump", jumps), ("Traction", tractions)):
        for k, comp in enumerate(comps):
            fields[f"{prefix}_{comp}"] = v[:, k]
    values = {}
    for name in names:
        col = np.full(n_gp, np.nan)
        for g in range(n_gp):
            sv = states[g] if states is not None else None
            if isinstance(sv, dict) and name in sv and _is_number(sv[name]):
                col[g] = float(sv[name])
        values[name] = col
    return (X, disp, np.full((n_gp, 6), np.nan), np.full(n_gp, elem.id),
            np.arange(1, n_gp + 1), values, fields)


def _gauss_displacements(Nf: np.ndarray, u_nodes: np.ndarray) -> np.ndarray:
    """``u(ξ_g) = Σ_a N_a(ξ_g)·u_a``: ``(…, n_gp, n_dof_nodo)``."""
    return np.einsum("gn,...nd->...gd", Nf, u_nodes)


def _family_points(fam, members: list, U: np.ndarray, trial: bool, names: list[str]):
    """Bloque de la nube de puntos con los elementos ``members`` de una
    familia, por arreglos."""
    Nf = fam.gauss_shape_functions()
    if Nf is None:
        return None
    index = {id(e): i for i, e in enumerate(fam.elements)}
    rows_e = np.array([index[id(e)] for e in members], dtype=np.int64)
    n_gp = fam.n_gp
    X = fam.gauss_points_global()[rows_e].reshape(-1, fam.spatial_dim)
    dof_names = tuple(fam.element_class.DOF_NAMES)
    n_per_node = len(dof_names)
    # Sólo las traslaciones (un elemento multicampo lleva también T).
    trans = [k for k, d in enumerate(dof_names) if d in _TRANSLATIONAL_DOFS]
    u_nodes = U[fam.dof_indices[rows_e]].reshape(len(rows_e), fam.n_nodes, n_per_node)
    disp = _gauss_displacements(Nf, u_nodes[:, :, trans]).reshape(-1, len(trans))
    st = fam.state
    rows = (rows_e[:, None] * n_gp + np.arange(n_gp)[None, :]).reshape(-1)
    S = (st.S_trial if trial else st.S_committed)[rows]
    sig_v = (st.sig_trial if trial else st.sig_committed)[rows]
    szz = None
    if sig_v.shape[1] == 3:
        szz = np.asarray(fam.material.batch_out_of_plane_stress(sig_v, S), dtype=float)
    sig = _six_components(sig_v, szz)
    ids = np.repeat([e.id for e in members], n_gp)
    point = np.tile(np.arange(1, n_gp + 1), len(members))
    schema = fam.schema
    values = {}
    for name in names:
        col = np.full(len(rows), np.nan)
        if name in schema.names:
            k = schema.names.index(name)
            if schema.shapes[k] == ():
                col = S[:, schema.offsets[k]].astype(float)
        values[name] = col
    return X, disp, sig, ids, point, values


def _element_points(elem, U: np.ndarray, trial: bool, names: list[str]):
    """Bloque de la nube de puntos de un elemento sin familia."""
    pts = getattr(elem, "points", None)
    if pts is None:
        return None
    pts = np.asarray(pts, dtype=float)
    pts = pts.reshape(pts.shape[0], -1)
    Nf = elem.batch_shape_functions(pts)
    st = getattr(elem, "state", None)
    if Nf is None or st is None:
        return None
    stresses = st.stresses_trial if trial else st.stresses
    states = st.vars_trial if trial else st.vars
    n_gp = len(stresses)
    if n_gp != Nf.shape[0]:
        return None
    dim = pts.shape[1]
    coords = np.asarray(elem.get_coordinate_matrix(ndim=dim), dtype=float)
    X = Nf @ coords
    trans = [d for d in _TRANSLATIONAL_DOFS[:dim]]
    u_nodes = np.array([[U[nd.dofs[d]] if d in nd.dofs else 0.0 for d in trans]
                        for nd in elem.nodes])
    disp = Nf @ u_nodes
    sig = np.zeros((n_gp, 6))
    for g in range(n_gp):
        s = np.atleast_1d(np.asarray(stresses[g], dtype=float))
        if s.size == 3:
            sv = states[g] if states is not None else None
            szz = float(elem.material.out_of_plane_stress(s, sv))
            sig[g] = _six_components(s.reshape(1, 3), np.array([szz]))[0]
        elif s.size == 6:
            sig[g] = s
    values = {}
    for name in names:
        col = np.full(n_gp, np.nan)
        for g in range(n_gp):
            sv = states[g] if states is not None else None
            if isinstance(sv, dict) and name in sv and _is_number(sv[name]):
                col[g] = float(sv[name])
        values[name] = col
    return X, disp, sig, np.full(n_gp, elem.id), np.arange(1, n_gp + 1), values


def _smooth_stress_to_nodes(n_nodes: int, conn_by_type: dict, stress_by_type: dict,
                            solid_by_type: dict):
    """Promedio nodal de las seis componentes de σ sobre las celdas de
    sólidos mecánicos adyacentes (``solid_by_type[tipo][k]`` verdadero).

    Devuelve ``(n_nodes, 6)``, o ``None`` si el dominio no tiene ninguna
    celda de sólido mecánico. Las líneas (barras, cables, marcos) y los
    elementos térmicos no participan, tampoco cuando un térmico comparte el
    tipo de celda y los nodos con un sólido mecánico (dominio multicampo
    ``(ux, uy, T)`` del ADR 0016): su ``σ`` es nulo por construcción y
    promediarlo dividiría el esfuerzo nodal por el número de capas. Un nodo
    sin ningún sólido mecánico adyacente vale cero.
    """
    acc = np.zeros((n_nodes, 6))
    cnt = np.zeros(n_nodes, dtype=np.int64)
    has_data = False
    for cell_type, conn_list in conn_by_type.items():
        solid = np.asarray(solid_by_type.get(cell_type, ()), dtype=bool)
        if not solid.any():
            continue
        has_data = True
        conn = np.asarray(conn_list, dtype=np.int64)          # (n_cells, n_por_celda)
        s = np.asarray(stress_by_type[cell_type], dtype=float)  # (n_cells, 6)
        if not solid.all():
            conn, s = conn[solid], s[solid]
        n_per = conn.shape[1]
        np.add.at(acc, conn.ravel(), np.repeat(s, n_per, axis=0))
        cnt += np.bincount(conn.ravel(), minlength=n_nodes)
    if not has_data:
        return None
    mask = cnt > 0
    acc[mask] /= cnt[mask][:, None]
    return acc
