# fenix_fem/fenix/utils/yaml_parser.py
import yaml
import os
import numpy as np
from fenix.core.domain import Domain
import fenix.registry_initialization
from fenix.registry import MaterialRegistry, ElementRegistry, SolverRegistry, QuadratureRegistry
import fenix.math.integration  # Asegura el registro de las cuadraturas

class YamlParser:
    """Lector automatizado de modelos estructurales desde archivos YAML."""
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.domain = Domain()
        self.materials = {}
        self.solver_config = {}
        self.point_loads = {}
        self.point_loads_by_coord = {}
        self.point_loads_by_group = {}
        self.point_loads_by_node = []
        self.output_config = {}

    def _get_quadrature(self, rule: str) -> tuple:
        """Convierte un string del YAML a una regla de cuadratura inyectable."""
        if rule == "1x1":
            print("  [!] ADVERTENCIA: Se ha seleccionado integración reducida (1x1).")
            print("      Tenga cuidado con posibles modos de energía nula (Hourglassing).")
        
        return QuadratureRegistry.get(rule)

    def parse(self) -> Domain:
        print(f"Leyendo modelo desde: {self.filepath} ...")
        with open(self.filepath, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)

        self._parse_nodes(data)
        self._parse_materials(data)
        self._parse_mesh_or_elements(data)
        self._parse_boundary_conditions(data)
        self._parse_point_loads_config(data)
        
        self.output_config = data.get('output', {})
        self.solver_config = data.get('solver', {})
        return self.domain

    def _parse_nodes(self, data: dict):
        nodes_data = data.get('nodes', [])
        if not nodes_data: return
        if isinstance(nodes_data, list):
            for node_dict in nodes_data:
                self.domain.add_node(node_dict['id'], node_dict['coords'])
        else:
            raise ValueError("El bloque 'nodes' debe ser una lista de diccionarios.")

    def _parse_materials(self, data: dict):
        for mat_data in data.get('materials', []):
            mat_id = mat_data['id']
            mat_type = mat_data['type']
            kwargs = {k: v for k, v in mat_data.items() if k not in ('id', 'type')}
            self.materials[mat_id] = MaterialRegistry.create(mat_type, **kwargs)

    def _parse_mesh_or_elements(self, data: dict):
        mesh_file = data.get('mesh', None)
        if mesh_file:
            from fenix.utils.gmsh_parser import GmshParser
            base_dir = os.path.dirname(self.filepath)
            mesh_path = os.path.join(base_dir, mesh_file)
            
            if not os.path.exists(mesh_path):
                mesh_path_sin_ext = os.path.join(base_dir, mesh_file.replace('.msh', ''))
                if os.path.exists(mesh_path_sin_ext):
                    mesh_path = mesh_path_sin_ext
            
            default_mat_id = data.get('mesh_material', 1)
            default_thickness = float(data.get('mesh_thickness', 1.0))
            default_quad_str = data.get('mesh_quadrature', '2x2')
            
            default_material = self.materials.get(default_mat_id)
            default_quadrature = self._get_quadrature(default_quad_str)
            
            physical_props = {}
            for group_name, props in data.get('mesh_physical_groups', {}).items():
                mat = self.materials.get(props.get('material', default_mat_id))
                thick = float(props.get('thickness', default_thickness))
                quad = self._get_quadrature(props.get('quadrature', default_quad_str))
                physical_props[group_name] = (mat, thick, quad)
            
            gmsh_parser = GmshParser(mesh_path)
            self.domain = gmsh_parser.parse(default_material, default_thickness, physical_props, default_quadrature)
        else:
            elements_data = data.get('elements', [])
            if isinstance(elements_data, list):
                for elem_dict in elements_data:
                    elem_id = elem_dict['id']
                    e_type = elem_dict['type']
                    mat_id = elem_dict['material']
                    node_ids = elem_dict['nodes']
                    
                    nodes = [self.domain.get_node(nid) for nid in node_ids]
                    material = self.materials[mat_id]
                    
                    kwargs = {k: v for k, v in elem_dict.items() if k not in ('id', 'type', 'material', 'nodes')}
                    self.domain.add_element(ElementRegistry.create(e_type, element_id=elem_id, nodes=nodes, material=material, **kwargs))
            elif elements_data:
                raise ValueError("El bloque 'elements' debe ser una lista de diccionarios.")

    def _parse_boundary_conditions(self, data: dict):
        bcs_data = data.get('boundary_conditions', []) or []
        bcs_by_node = data.get('boundary_conditions_by_node', []) or []
        
        for bc in bcs_data + bcs_by_node:
            node_id = bc.get('node_id')
            if node_id is None:
                raise ValueError("Falta 'node_id' en una condición de frontera.")
            node = self.domain.get_node(node_id)
            if not node: 
                raise ValueError(f"Nodo {node_id} no existe en la malla para aplicar condición de frontera.")
            for dof, value in bc.items():
                if dof != 'node_id':
                    node.fix_dof(dof, float(value))

        bcs_coord = data.get('boundary_conditions_by_coord', [])
        if isinstance(bcs_coord, dict):
            parsed_bcs = []
            for k, v in bcs_coord.items():
                if k in ['x_min', 'x_max', 'y_min', 'y_max']: v['loc'] = k
                parsed_bcs.append(v)
            bcs_coord = parsed_bcs
            
        if bcs_coord and self.domain.nodes:
            x_coords = [n.coordinates[0] for n in self.domain.nodes.values() if n.dofs]
            y_coords = [n.coordinates[1] for n in self.domain.nodes.values() if n.dofs]
            if x_coords and y_coords:
                x_min, x_max = min(x_coords), max(x_coords)
                y_min, y_max = min(y_coords), max(y_coords)
                
                for node in self.domain.nodes.values():
                    if not node.dofs: continue
                    x, y = node.coordinates
                    
                    for bcs in bcs_coord:
                        tol = float(bcs.get('tol', 1e-6))
                        loc = bcs.get('loc')
                        match = False
                        if loc == 'x_min' and abs(x - x_min) < tol: match = True
                        elif loc == 'x_max' and abs(x - x_max) < tol: match = True
                        elif loc == 'y_min' and abs(y - y_min) < tol: match = True
                        elif loc == 'y_max' and abs(y - y_max) < tol: match = True
                        elif bcs.get('coord') == 'x' and 'val' in bcs and abs(x - float(bcs['val'])) < tol: match = True
                        elif bcs.get('coord') == 'y' and 'val' in bcs and abs(y - float(bcs['val'])) < tol: match = True
                        
                        if match:
                            for dof, value in bcs.items():
                                if dof not in ['tol', 'coord', 'val', 'loc']: node.fix_dof(dof, float(value))

        bcs_group = data.get('boundary_conditions_by_group', [])
        if isinstance(bcs_group, dict):
            bcs_group = [{'group_name': k, **v} for k, v in bcs_group.items()]
            
        if bcs_group and hasattr(self.domain, 'physical_groups'):
            for bcs in bcs_group:
                group_name = bcs.get('group_name')
                if group_name and group_name in self.domain.physical_groups:
                    for node_id in self.domain.physical_groups[group_name]:
                        node = self.domain.get_node(node_id)
                        if node:
                            for dof, value in bcs.items():
                                if dof not in ['group_name']: node.fix_dof(dof, float(value))

    def _parse_point_loads_config(self, data: dict):
        self.point_loads = data.get('point_loads', [])
        self.point_loads_by_node = data.get('point_loads_by_node', [])
        
        p_loads_coord = data.get('point_loads_by_coord', [])
        self.point_loads_by_coord = list(p_loads_coord.values()) if isinstance(p_loads_coord, dict) else p_loads_coord
            
        p_loads_group = data.get('point_loads_by_group', [])
        if isinstance(p_loads_group, dict):
            self.point_loads_by_group = [{'group_name': k, **v} for k, v in p_loads_group.items()]
        else:
            self.point_loads_by_group = p_loads_group
        
    def get_external_forces(self) -> np.ndarray:
        """Construye el vector de fuerzas externas global F_ext."""
        F_ext = np.zeros(self.domain.total_dofs)
        
        for load in self.point_loads + self.point_loads_by_node:
            node_id = load.get('node_id')
            if node_id is None: continue
            node = self.domain.get_node(node_id)
            if not node: continue
            for dof, value in load.items():
                if dof != 'node_id' and dof in node.dofs:
                    F_ext[node.dofs[dof]] += float(value)
                    
        if self.point_loads_by_coord and self.domain.nodes:
            x_coords = [n.coordinates[0] for n in self.domain.nodes.values() if n.dofs]
            y_coords = [n.coordinates[1] for n in self.domain.nodes.values() if n.dofs]
            if x_coords and y_coords:
                x_min, x_max = min(x_coords), max(x_coords)
                y_min, y_max = min(y_coords), max(y_coords)
                
                for node in self.domain.nodes.values():
                    if not node.dofs: continue
                    x, y = node.coordinates
                    for loads in self.point_loads_by_coord:
                        tol = float(loads.get('tol', 1e-6))
                        loc = loads.get('loc')
                        match = False
                        if loc == 'x_min' and abs(x - x_min) < tol: match = True
                        elif loc == 'x_max' and abs(x - x_max) < tol: match = True
                        elif loc == 'y_min' and abs(y - y_min) < tol: match = True
                        elif loc == 'y_max' and abs(y - y_max) < tol: match = True
                        elif loads.get('coord') == 'x' and 'val' in loads and abs(x - float(loads['val'])) < tol: match = True
                        elif loads.get('coord') == 'y' and 'val' in loads and abs(y - float(loads['val'])) < tol: match = True
                        
                        if match:
                            for dof, value in loads.items():
                                if dof not in ['tol', 'coord', 'val', 'loc'] and dof in node.dofs:
                                    F_ext[node.dofs[dof]] += float(value)
                                    
        if self.point_loads_by_group and hasattr(self.domain, 'physical_groups'):
            for loads in self.point_loads_by_group:
                group_name = loads.get('group_name')
                if group_name and group_name in self.domain.physical_groups:
                    for node_id in self.domain.physical_groups[group_name]:
                        node = self.domain.get_node(node_id)
                        if node:
                            for dof, value in loads.items():
                                if dof not in ['group_name'] and dof in node.dofs:
                                    F_ext[node.dofs[dof]] += float(value)
                                    
        return F_ext

    def get_solver(self, assembler):
        """Construye y retorna el solver dinámicamente según la configuración YAML."""
        s_type = self.solver_config.get('type', 'LinearSolver')
        kwargs = {k: v for k, v in self.solver_config.items() if k != 'type'}
        return SolverRegistry.create(s_type, assembler=assembler, **kwargs)