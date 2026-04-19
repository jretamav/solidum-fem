# fenix_fem/fenix/utils/vtk_exporter.py
import numpy as np
from fenix.core.domain import Domain
from fenix.elements.solid_2d import Quad4, Tri3
from fenix.elements.structural import Truss2D, Truss3D

try:
    import meshio
except ImportError:
    meshio = None

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


class VtkExporter:
    """Exporta el dominio y los resultados a formato VTK para visualización en ParaView."""
    def __init__(self, domain: Domain):
        self.domain = domain

    def export(self, filepath: str, U: np.ndarray = None, F_ext: np.ndarray = None, nodal_vars: list = None, elem_vars: list = None):
        if meshio is None:
            print("  -> Advertencia: La librería 'meshio' no está instalada. Omita la exportación VTK.")
            return
            
        if U is None:
            U = np.zeros(self.domain.total_dofs) if self.domain.total_dofs > 0 else np.zeros(0)
        if F_ext is None:
            F_ext = np.zeros(self.domain.total_dofs) if self.domain.total_dofs > 0 else np.zeros(0)

        # 1. Extraer Nodos (Puntos) y mapear sus IDs a índices de arreglos (0, 1, 2...)
        node_ids = list(self.domain.nodes.keys())
        node_map = {nid: i for i, nid in enumerate(node_ids)}
        
        points = np.zeros((len(node_ids), 3))
        displacements = np.zeros((len(node_ids), 3))
        supports = np.zeros((len(node_ids), 3))
        forces = np.zeros((len(node_ids), 3))
        
        for nid, node in self.domain.nodes.items():
            idx = node_map[nid]
            points[idx, 0] = node.coordinates[0]
            points[idx, 1] = node.coordinates[1]
            points[idx, 2] = node.coordinates[2] if len(node.coordinates) > 2 else 0.0
            
            if len(U) > 0:
                if 'ux' in node.dofs and node.dofs['ux'] < len(U): displacements[idx, 0] = U[node.dofs['ux']]
                if 'uy' in node.dofs and node.dofs['uy'] < len(U): displacements[idx, 1] = U[node.dofs['uy']]
                
            if 'ux' in node.boundary_conditions: supports[idx, 0] = 1.0
            if 'uy' in node.boundary_conditions: supports[idx, 1] = 1.0
                
            if len(F_ext) > 0:
                if 'ux' in node.dofs and node.dofs['ux'] < len(F_ext): forces[idx, 0] = F_ext[node.dofs['ux']]
                if 'uy' in node.dofs and node.dofs['uy'] < len(F_ext): forces[idx, 1] = F_ext[node.dofs['uy']]

        # 2. Extraer Conectividad y Variables de Estado (Celdas)
        quad_conn = []
        line_conn = []
        quad_state = []
        line_state = []
        quad_stresses = []
        tri_conn = []
        tri_state = []
        tri_stresses = []
        
        for elem in self.domain.elements.values():
            conn = [node_map[n.id] for n in elem.nodes]
            
            # Extraer variable histórica (ej. 'alpha' para plasticidad o 'damage' para daño)
            # Compatible tanto con la API ElementState (elem.state.vars) como
            # con elementos que aún exponen elem.state_vars directamente.
            state_val = 0.0
            state_vars_list = None
            if hasattr(elem, 'state') and hasattr(elem.state, 'vars'):
                state_vars_list = elem.state.vars
            elif hasattr(elem, 'state_vars'):
                state_vars_list = elem.state_vars
            if state_vars_list is not None:
                primary_key = getattr(getattr(elem, 'material', None), 'PRIMARY_STATE_VAR', None)
                vals = [_extract_state_scalar(sv, primary_key)
                        for sv in state_vars_list if sv is not None]
                if vals:
                    state_val = sum(vals) / len(vals)

            # Extraer esfuerzos comprometidos (committed) desde ElementState
            def _avg_stress_from_elem(el):
                if hasattr(el, 'state') and hasattr(el.state, 'stresses'):
                    committed = [s for s in el.state.stresses if s is not None]
                    if committed:
                        return np.mean(committed, axis=0)
                elif hasattr(el, 'stresses'):
                    return np.mean(el.stresses, axis=0)
                return None

            if isinstance(elem, Quad4):
                quad_conn.append(conn)
                quad_state.append(state_val)
                s_avg = _avg_stress_from_elem(elem)
                if s_avg is not None:
                    sx, sy, txy = s_avg[0], s_avg[1], s_avg[2]
                    vm = np.sqrt(sx**2 + sy**2 - sx*sy + 3.0*txy**2)
                    quad_stresses.append([sx, sy, txy, vm])
                else:
                    quad_stresses.append([0.0, 0.0, 0.0, 0.0])

            elif isinstance(elem, Tri3):
                tri_conn.append(conn)
                tri_state.append(state_val)
                s_avg = _avg_stress_from_elem(elem)
                if s_avg is not None:
                    sx, sy, txy = s_avg[0], s_avg[1], s_avg[2]
                    vm = np.sqrt(sx**2 + sy**2 - sx*sy + 3.0*txy**2)
                    tri_stresses.append([sx, sy, txy, vm])
                else:
                    tri_stresses.append([0.0, 0.0, 0.0, 0.0])

            elif isinstance(elem, (Truss2D, Truss3D)):
                line_conn.append(conn)
                line_state.append(state_val)
            
        cells = []
        state_arrays, sxx_arrays, syy_arrays, txy_arrays, vm_arrays = [], [], [], [], []
        
        if quad_conn: 
            cells.append(("quad", np.array(quad_conn, dtype=int)))
            state_arrays.append(np.array(quad_state, dtype=float))
            q_s = np.array(quad_stresses, dtype=float)
            sxx_arrays.append(q_s[:, 0])
            syy_arrays.append(q_s[:, 1])
            txy_arrays.append(q_s[:, 2])
            vm_arrays.append(q_s[:, 3])
            
        if tri_conn:
            cells.append(("triangle", np.array(tri_conn, dtype=int)))
            state_arrays.append(np.array(tri_state, dtype=float))
            t_s = np.array(tri_stresses, dtype=float)
            sxx_arrays.append(t_s[:, 0])
            syy_arrays.append(t_s[:, 1])
            txy_arrays.append(t_s[:, 2])
            vm_arrays.append(t_s[:, 3])
            
        if line_conn: 
            cells.append(("line", np.array(line_conn, dtype=int)))
            state_arrays.append(np.array(line_state, dtype=float))
            sxx_arrays.append(np.zeros(len(line_conn)))
            syy_arrays.append(np.zeros(len(line_conn)))
            txy_arrays.append(np.zeros(len(line_conn)))
            vm_arrays.append(np.zeros(len(line_conn)))

        # Filtrado inteligente de variables basado en las peticiones del YAML
        point_data = {}
        if nodal_vars is None or 'Displacements' in nodal_vars: point_data["Displacements"] = displacements
        if nodal_vars is None or 'Supports' in nodal_vars: point_data["Supports"] = supports
        if nodal_vars is None or 'External_Forces' in nodal_vars: point_data["External_Forces"] = forces
            
        cell_data = {}
        if state_arrays and (elem_vars is None or 'Internal_State' in elem_vars): cell_data["Internal_State"] = state_arrays
        if state_arrays and (elem_vars is None or 'Sigma_XX' in elem_vars): cell_data["Sigma_XX"] = sxx_arrays
        if state_arrays and (elem_vars is None or 'Sigma_YY' in elem_vars): cell_data["Sigma_YY"] = syy_arrays
        if state_arrays and (elem_vars is None or 'Tau_XY' in elem_vars): cell_data["Tau_XY"] = txy_arrays
        if state_arrays and (elem_vars is None or 'Von_Mises' in elem_vars): cell_data["Von_Mises"] = vm_arrays
            
        if not cell_data:
            cell_data = None
        
        mesh = meshio.Mesh(points, cells, point_data=point_data, cell_data=cell_data)
        mesh.write(filepath)