"""Export del dominio y resultados a formato VTK (ParaView, vía ``meshio``).

Cobertura
---------
La celda VTK de cada elemento se **deriva de su contrato**, sin listas de
clases: número de nodos y dimensión espacial (``'uz' in DOF_NAMES`` en los
mecánicos, ``FLUX_DIM`` en los térmicos).

- 2D: 3 → ``triangle``, 4 → ``quad``, 6 → ``triangle6``, 8 → ``quad8``,
  9 → ``quad9``. Cubre Tri3, Quad4, Tri6, Quad8, Quad9 y ``Quad4Thermal``
  (y cualquier elemento de un módulo de usuario con esos nodos).
- 3D: 4 → ``tetra``, 8 → ``hexahedron``, 10 → ``tetra10``,
  20 → ``hexahedron20``, 27 → ``hexahedron27``. Cubre Tet4, Hex8, Tet10,
  Hex20, Hex27 y ``Hex8Thermal``. La numeración de nodos medios del
  proyecto coincide con la de VTK (verificado contra meshio).
- Cualquier elemento de 2 nodos (trusses, cables, frames 2D/3D) → ``line``.
  No aportan campos de esfuerzo; sus fuerzas internas se consumen vía
  ``SolveResult.element_forces`` (ADR 0002).
- Un elemento puede forzar su celda declarando ``VTK_CELL_TYPE`` como
  atributo de clase (p. ej. una cuña de 6 nodos en 3D, que la tabla
  confundiría con ``triangle6``).

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

Campos nodales
--------------
- ``Displacements``: vector 3D ``(ux, uy, uz)``; ``uz = 0`` en modelos 2D.
- ``Rotations``: ``(rx, ry, rz)``, sólo si algún nodo tiene DOF rotacional.
- ``Temperature``: escalar, sólo si algún nodo tiene el DOF ``T``.
- ``Supports`` / ``Supports_Rot``: 1.0 en cada componente con Dirichlet.
- ``External_Forces``: fuerzas nodales aplicadas ``(fx, fy, fz)``.
- ``Sigma_XX_nodal`` … ``Von_Mises_nodal``: suavizado nodal de orden 0
  (media simple sobre las celdas adyacentes) de los campos anteriores.
"""
from __future__ import annotations

import numpy as np

from solidum.core.domain import Domain
from solidum.logging import get_logger
from solidum.math.batch.postprocess import group_by_family

try:
    import meshio
except ImportError:
    meshio = None

_log = get_logger("exporters.vtk")

_TRANSLATIONAL_DOFS = ("ux", "uy", "uz")
_ROTATIONAL_DOFS = ("rx", "ry", "rz")
_TEMPERATURE_DOF = "T"

# Celda VTK por (dimensión espacial, número de nodos). Los elementos de dos
# nodos se resuelven aparte como ``line``.
_CELL_2D = {3: "triangle", 4: "quad", 6: "triangle6", 8: "quad8", 9: "quad9"}
_CELL_3D = {4: "tetra", 8: "hexahedron", 10: "tetra10",
            20: "hexahedron20", 27: "hexahedron27"}

# Nombres y orden de las seis componentes de esfuerzo exportadas.
STRESS_COMPONENT_NAMES = ("Sigma_XX", "Sigma_YY", "Sigma_ZZ",
                          "Tau_XY", "Tau_YZ", "Tau_XZ")


# ---------------------------------------------------------------------------
# Clasificación de elementos
# ---------------------------------------------------------------------------

def _is_thermal(elem) -> bool:
    return hasattr(elem, "FLUX_DIM")


def _spatial_dim(elem) -> int:
    if _is_thermal(elem):
        return int(elem.FLUX_DIM)
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
    return (not _is_thermal(elem)) and len(elem.nodes) > 2 \
        and hasattr(elem, "compute_gauss_state")


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


def _extract_state_scalar(sv: dict, primary_key: str | None = None) -> float:
    """Escalar representativo del dict de variables de estado de un punto.

    Prioriza ``primary_key``; si no está, el primer escalar numérico.
    """
    if not isinstance(sv, dict):
        return 0.0
    if primary_key and primary_key in sv:
        v = sv[primary_key]
        if isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool):
            return float(v)
    for v in sv.values():
        if isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool):
            return float(v)
    return 0.0


def _state_var_avg(elem) -> float:
    """Promedio de la ``PRIMARY_STATE_VAR`` del material sobre el estado committed."""
    state_vars_list = getattr(getattr(elem, "state", None), "vars", None)
    if state_vars_list is None:
        return 0.0
    primary_key = getattr(getattr(elem, "material", None), "PRIMARY_STATE_VAR", None)
    vals = [_extract_state_scalar(sv, primary_key)
            for sv in state_vars_list if sv is not None]
    if not vals:
        return 0.0
    return sum(vals) / len(vals)


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
               nodal_vars: list | None = None, elem_vars: list | None = None):
        """Escribe ``filepath`` (``.vtu``, ``.vtk``… según extensión).

        Parameters
        ----------
        U, F_ext
            Campo global de la solución y de fuerzas aplicadas. ``None``
            exporta la malla sin deformar con campos nulos (pre-proceso).
        nodal_vars, elem_vars
            Filtro opcional de nombres de campos nodales / por celda. ``None``
            exporta todos los disponibles.
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

        # --- Celdas agrupadas por tipo -----------------------------------------
        elements = list(self.domain.elements.values())
        conn_by_type: dict[str, list[list[int]]] = {}
        slot: list[tuple[str, int] | None] = []      # (tipo, índice en su bloque)
        for elem in elements:
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
            slot.append((cell_type, len(block) - 1))

        if not conn_by_type:
            raise ValueError(
                "VtkExporter: ningún elemento del dominio tiene celda VTK "
                "soportada; no se escribe el archivo (meshio no puede releer "
                "una malla sin celdas)."
            )

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

        # Suavizado nodal de orden 0: media simple de los sólidos adyacentes.
        nodal = _smooth_stress_to_nodes(n_nodes, conn_by_type, stress_by_type)
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
        if not cell_data:
            cell_data = None

        mesh = meshio.Mesh(points, cells, point_data=point_data, cell_data=cell_data)
        mesh.write(filepath)


def _smooth_stress_to_nodes(n_nodes: int, conn_by_type: dict, stress_by_type: dict):
    """Promedio nodal de las seis componentes de σ sobre las celdas sólidas
    adyacentes. Devuelve ``(n_nodes, 6)`` o ``None`` si no hay sólidos
    mecánicos con esfuerzo no nulo (las líneas y los térmicos no participan)."""
    acc = np.zeros((n_nodes, 6))
    cnt = np.zeros(n_nodes, dtype=np.int64)
    has_data = False
    for cell_type, conn_list in conn_by_type.items():
        if cell_type == "line" or not conn_list:
            continue
        has_data = True
        conn = np.asarray(conn_list, dtype=np.int64)          # (n_cells, n_por_celda)
        s = np.asarray(stress_by_type[cell_type], dtype=float)  # (n_cells, 6)
        n_per = conn.shape[1]
        np.add.at(acc, conn.ravel(), np.repeat(s, n_per, axis=0))
        cnt += np.bincount(conn.ravel(), minlength=n_nodes)
    if not has_data:
        return None
    mask = cnt > 0
    acc[mask] /= cnt[mask][:, None]
    return acc
