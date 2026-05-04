# fenix_fem/fenix/utils/vtk_exporter.py
"""Export del dominio y resultados a formato VTK para ParaView.

Cobertura
---------
- Sólidos 2D (`Quad4`, `Tri3`): celda VTK ``quad`` / ``triangle``, con campos
  por celda Sigma_XX, Sigma_YY, Tau_XY, Von_Mises (promediados sobre puntos
  de Gauss) y la PRIMARY_STATE_VAR del material.
- Cualquier elemento de **2 nodos** (`Truss2D`, `Truss3D`, `Truss2DCorot`,
  `Truss3DCorot`, `Cable2DCorot`, `Cable3DCorot`, `Frame2DEuler`,
  `Frame2DTimoshenko`, `Frame2DEulerCorot`, `Frame3D`): celda ``line``. Los
  frames/cables no aportan campos de esfuerzos al VTK; quien quiera N/V/M
  los consume vía ``SolveResult.element_forces`` (ADR 0002).

Campos nodales
--------------
- ``Displacements``: vector 3D (ux, uy, uz) — uz=0 si el modelo es 2D.
- ``Rotations``: vector 3D (rx, ry, rz). Sólo se incluye si al menos un nodo
  tiene algún DOF rotacional (rx/ry/rz). uz y rotaciones permiten visualizar
  correctamente frames 3D y vigas corotacionales en ParaView.
- ``Supports``: vector 3D con 1.0 en cada componente con condición de
  Dirichlet en ux/uy/uz.
- ``External_Forces``: vector 3D de fuerzas nodales aplicadas (fx, fy, fz).
"""
import numpy as np
from fenix.core.domain import Domain
from fenix.elements.solid_2d import Quad4, Tri3
from fenix.logging import get_logger

try:
    import meshio
except ImportError:
    meshio = None

_log = get_logger("exporters.vtk")

_TRANSLATIONAL_DOFS = ("ux", "uy", "uz")
_ROTATIONAL_DOFS = ("rx", "ry", "rz")


def _extract_state_scalar(sv: dict, primary_key: str = None) -> float:
    """Extrae un escalar del dict de variables de estado de un punto de integración.

    Prioriza ``primary_key`` si está definido y presente en ``sv``.
    Si no, devuelve el primer valor numérico escalar encontrado en el dict.
    Retorna 0.0 si el dict está vacío o no contiene escalares numéricos.
    """
    if primary_key and primary_key in sv:
        v = sv[primary_key]
        if isinstance(v, (int, float, np.floating, np.integer)):
            return float(v)
    for v in sv.values():
        if isinstance(v, (int, float, np.floating, np.integer)) and not isinstance(v, bool):
            return float(v)
    return 0.0


def _avg_stress_from_elem(el):
    """Media de los esfuerzos comprometidos en los puntos de Gauss del elemento."""
    if hasattr(el, 'state') and hasattr(el.state, 'stresses'):
        committed = [s for s in el.state.stresses if s is not None]
        if committed:
            return np.mean(committed, axis=0)
    elif hasattr(el, 'stresses'):
        return np.mean(el.stresses, axis=0)
    return None


def _state_var_avg(elem) -> float:
    """Promedio de la PRIMARY_STATE_VAR del material sobre puntos de Gauss."""
    state_vars_list = None
    if hasattr(elem, 'state') and hasattr(elem.state, 'vars'):
        state_vars_list = elem.state.vars
    elif hasattr(elem, 'state_vars'):
        state_vars_list = elem.state_vars
    if state_vars_list is None:
        return 0.0
    primary_key = getattr(getattr(elem, 'material', None), 'PRIMARY_STATE_VAR', None)
    vals = [_extract_state_scalar(sv, primary_key)
            for sv in state_vars_list if sv is not None]
    if not vals:
        return 0.0
    return sum(vals) / len(vals)


class VtkExporter:
    """Exporta el dominio y los resultados a formato VTK para visualización en ParaView."""

    def __init__(self, domain: Domain):
        self.domain = domain

    def export(self, filepath: str, U: np.ndarray = None, F_ext: np.ndarray = None,
               nodal_vars: list = None, elem_vars: list = None):
        if meshio is None:
            _log.warning("La librería 'meshio' no está instalada. Omitiendo exportación VTK.")
            return

        if U is None:
            U = np.zeros(self.domain.total_dofs) if self.domain.total_dofs > 0 else np.zeros(0)
        if F_ext is None:
            F_ext = np.zeros(self.domain.total_dofs) if self.domain.total_dofs > 0 else np.zeros(0)

        node_ids = list(self.domain.nodes.keys())
        node_map = {nid: i for i, nid in enumerate(node_ids)}
        n_nodes = len(node_ids)

        points = np.zeros((n_nodes, 3))
        displacements = np.zeros((n_nodes, 3))
        rotations = np.zeros((n_nodes, 3))
        supports = np.zeros((n_nodes, 3))
        forces = np.zeros((n_nodes, 3))

        any_rotational = False

        for nid, node in self.domain.nodes.items():
            idx = node_map[nid]
            points[idx, 0] = node.coordinates[0]
            points[idx, 1] = node.coordinates[1]
            points[idx, 2] = node.coordinates[2] if len(node.coordinates) > 2 else 0.0

            for axis, dof_name in enumerate(_TRANSLATIONAL_DOFS):
                if dof_name in node.dofs and node.dofs[dof_name] < len(U):
                    displacements[idx, axis] = U[node.dofs[dof_name]]
                if dof_name in node.boundary_conditions:
                    supports[idx, axis] = 1.0
                if dof_name in node.dofs and node.dofs[dof_name] < len(F_ext):
                    forces[idx, axis] = F_ext[node.dofs[dof_name]]

            for axis, dof_name in enumerate(_ROTATIONAL_DOFS):
                if dof_name in node.dofs:
                    any_rotational = True
                    if node.dofs[dof_name] < len(U):
                        rotations[idx, axis] = U[node.dofs[dof_name]]

        # Conectividad por tipo de celda VTK. Cualquier elemento con 2 nodos
        # se proyecta a `line`, lo que cubre trusses, cables y frames de
        # cualquier dimensión sin enumerar las clases una a una.
        quad_conn, tri_conn, line_conn = [], [], []
        quad_state, tri_state, line_state = [], [], []
        quad_stresses, tri_stresses = [], []

        for elem in self.domain.elements.values():
            conn = [node_map[n.id] for n in elem.nodes]
            state_val = _state_var_avg(elem)

            if isinstance(elem, Quad4):
                quad_conn.append(conn)
                quad_state.append(state_val)
                s_avg = _avg_stress_from_elem(elem)
                quad_stresses.append(_voigt2d_to_components(s_avg))
            elif isinstance(elem, Tri3):
                tri_conn.append(conn)
                tri_state.append(state_val)
                s_avg = _avg_stress_from_elem(elem)
                tri_stresses.append(_voigt2d_to_components(s_avg))
            elif len(elem.nodes) == 2:
                line_conn.append(conn)
                line_state.append(state_val)
            else:
                _log.warning(
                    f"VtkExporter: elemento id={getattr(elem, 'id', '?')} de tipo "
                    f"{type(elem).__name__} con {len(elem.nodes)} nodos no tiene celda "
                    f"VTK soportada; se omite."
                )

        cells = []
        state_arrays, sxx_arrays, syy_arrays, txy_arrays, vm_arrays = [], [], [], [], []

        if quad_conn:
            cells.append(("quad", np.array(quad_conn, dtype=int)))
            state_arrays.append(np.array(quad_state, dtype=float))
            q_s = np.array(quad_stresses, dtype=float)
            sxx_arrays.append(q_s[:, 0]); syy_arrays.append(q_s[:, 1])
            txy_arrays.append(q_s[:, 2]); vm_arrays.append(q_s[:, 3])

        if tri_conn:
            cells.append(("triangle", np.array(tri_conn, dtype=int)))
            state_arrays.append(np.array(tri_state, dtype=float))
            t_s = np.array(tri_stresses, dtype=float)
            sxx_arrays.append(t_s[:, 0]); syy_arrays.append(t_s[:, 1])
            txy_arrays.append(t_s[:, 2]); vm_arrays.append(t_s[:, 3])

        if line_conn:
            cells.append(("line", np.array(line_conn, dtype=int)))
            state_arrays.append(np.array(line_state, dtype=float))
            zeros_line = np.zeros(len(line_conn))
            sxx_arrays.append(zeros_line); syy_arrays.append(zeros_line)
            txy_arrays.append(zeros_line); vm_arrays.append(zeros_line)

        point_data = {}
        if nodal_vars is None or 'Displacements' in nodal_vars:
            point_data["Displacements"] = displacements
        if any_rotational and (nodal_vars is None or 'Rotations' in nodal_vars):
            point_data["Rotations"] = rotations
        if nodal_vars is None or 'Supports' in nodal_vars:
            point_data["Supports"] = supports
        if nodal_vars is None or 'External_Forces' in nodal_vars:
            point_data["External_Forces"] = forces

        # Suavizado nodal de los esfuerzos: media simple sobre los elementos
        # adyacentes a cada nodo. Es nodal averaging de orden 0 — exacto
        # cuando el campo es constante por elemento, aproximado en general.
        # Útil para mapas continuos en ParaView sin pasar por Cell-To-Point.
        nodal_sigma = _smooth_stress_to_nodes(
            n_nodes,
            (quad_conn, quad_stresses),
            (tri_conn, tri_stresses),
        )
        if nodal_sigma is not None:
            sxx_n, syy_n, txy_n, vm_n = nodal_sigma
            if nodal_vars is None or 'Sigma_XX_nodal' in nodal_vars:
                point_data["Sigma_XX_nodal"] = sxx_n
            if nodal_vars is None or 'Sigma_YY_nodal' in nodal_vars:
                point_data["Sigma_YY_nodal"] = syy_n
            if nodal_vars is None or 'Tau_XY_nodal' in nodal_vars:
                point_data["Tau_XY_nodal"] = txy_n
            if nodal_vars is None or 'Von_Mises_nodal' in nodal_vars:
                point_data["Von_Mises_nodal"] = vm_n

        cell_data = {}
        if state_arrays and (elem_vars is None or 'Internal_State' in elem_vars):
            cell_data["Internal_State"] = state_arrays
        if state_arrays and (elem_vars is None or 'Sigma_XX' in elem_vars):
            cell_data["Sigma_XX"] = sxx_arrays
        if state_arrays and (elem_vars is None or 'Sigma_YY' in elem_vars):
            cell_data["Sigma_YY"] = syy_arrays
        if state_arrays and (elem_vars is None or 'Tau_XY' in elem_vars):
            cell_data["Tau_XY"] = txy_arrays
        if state_arrays and (elem_vars is None or 'Von_Mises' in elem_vars):
            cell_data["Von_Mises"] = vm_arrays

        if not cell_data:
            cell_data = None

        mesh = meshio.Mesh(points, cells, point_data=point_data, cell_data=cell_data)
        mesh.write(filepath)


def _smooth_stress_to_nodes(n_nodes: int, *families):
    """Promedio nodal del σ por elemento sobre los elementos adyacentes.

    families: tuplas ``(conn_list, stresses_list)`` con la conectividad y los
    [σxx, σyy, σxy, σ_VM] por celda. Devuelve ``(sxx, syy, txy, vm)`` como
    arrays nodales o ``None`` si ninguna familia aporta datos.
    """
    sxx = np.zeros(n_nodes)
    syy = np.zeros(n_nodes)
    txy = np.zeros(n_nodes)
    cnt = np.zeros(n_nodes, dtype=np.int64)
    has_data = False

    for conn_list, stress_list in families:
        if not conn_list:
            continue
        has_data = True
        for conn, sv in zip(conn_list, stress_list):
            sx, sy, sxy, _vm = sv
            for nidx in conn:
                sxx[nidx] += sx
                syy[nidx] += sy
                txy[nidx] += sxy
                cnt[nidx] += 1

    if not has_data:
        return None

    mask = cnt > 0
    sxx[mask] /= cnt[mask]
    syy[mask] /= cnt[mask]
    txy[mask] /= cnt[mask]
    vm = np.sqrt(sxx * sxx + syy * syy - sxx * syy + 3.0 * txy * txy)
    return sxx, syy, txy, vm


def _voigt2d_to_components(s_avg) -> list:
    """Empaqueta σ Voigt 2D promediada en [σxx, σyy, σxy, σ_VM]; ceros si None."""
    if s_avg is None:
        return [0.0, 0.0, 0.0, 0.0]
    sx, sy, txy = s_avg[0], s_avg[1], s_avg[2]
    vm = float(np.sqrt(sx * sx + sy * sy - sx * sy + 3.0 * txy * txy))
    return [float(sx), float(sy), float(txy), vm]
