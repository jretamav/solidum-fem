# fenix_fem/fenix/utils/yaml_parser.py
import yaml
import os
import numpy as np
from fenix.core.domain import Domain
from fenix.materials.von_mises_2d import VonMises2D
from fenix.materials.elastic_2d import Elastic2D
from fenix.elements.solid_2d import Quad4, Tri3
from fenix.math.solvers import LinearSolver, NonlinearSolver

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
        self.output_config = {}

    def parse(self) -> Domain:
        print(f"Leyendo modelo desde: {self.filepath} ...")
        with open(self.filepath, 'r') as f:
            data = yaml.safe_load(f)

        # 1. Instanciar Nodos
        for node_id, coords in data.get('nodes', {}).items():
            self.domain.add_node(node_id, coords)

        # 2. Instanciar Materiales
        for mat_data in data.get('materials', []):
            mat_id = mat_data['id']
            mat_type = mat_data['type']
            if mat_type == 'VonMises2D':
                self.materials[mat_id] = VonMises2D(
                    E=float(mat_data['E']),
                    nu=float(mat_data['nu']),
                    sigma_y=float(mat_data['sigma_y']),
                    H=float(mat_data.get('H', 0.0))
                )
            elif mat_type == 'Elastic2D':
                self.materials[mat_id] = Elastic2D(
                    E=float(mat_data['E']),
                    nu=float(mat_data['nu']),
                    hypothesis=mat_data.get('hypothesis', 'plane_stress')
                )

        # 3. Malla Externa (Gmsh) o Malla Manual
        mesh_file = data.get('mesh', None)
        if mesh_file:
            from fenix.utils.gmsh_parser import GmshParser
            base_dir = os.path.dirname(self.filepath)
            mesh_path = os.path.join(base_dir, mesh_file)
            
            # Evitar error si el archivo se llama "placa" en lugar de "placa.msh"
            if not os.path.exists(mesh_path):
                mesh_path_sin_ext = os.path.join(base_dir, mesh_file.replace('.msh', ''))
                if os.path.exists(mesh_path_sin_ext):
                    mesh_path = mesh_path_sin_ext
            
            mat_id = data.get('mesh_material', 1)
            thickness = float(data.get('mesh_thickness', 1.0))
            material = self.materials[mat_id]
            
            gmsh_parser = GmshParser(mesh_path)
            self.domain = gmsh_parser.parse(material, thickness)
        else:
            for elem_id, elem_data in data.get('elements', {}).items():
                e_type, mat_id, thickness = elem_data[0], elem_data[1], float(elem_data[2])
                node_ids = elem_data[3:]
                
                nodes = [self.domain.get_node(nid) for nid in node_ids]
                material = self.materials[mat_id]
                
                if e_type == 'Quad4':
                    self.domain.add_element(Quad4(elem_id, nodes, material, thickness))
                elif e_type == 'Tri3':
                    self.domain.add_element(Tri3(elem_id, nodes, material, thickness))
                else:
                    raise ValueError(f"Tipo de elemento {e_type} no soportado.")

        # 4. Aplicar Condiciones de Frontera (Desplazamientos)
        for node_id, bcs in data.get('boundary_conditions', {}).items():
            node = self.domain.get_node(node_id)
            for dof, value in bcs.items():
                node.fix_dof(dof, float(value))

        # 5. Aplicar Condiciones de Frontera por Coordenadas (Para Gmsh)
        bcs_coord = data.get('boundary_conditions_by_coord', {})
        if bcs_coord and self.domain.nodes:
            x_coords = [n.coordinates[0] for n in self.domain.nodes.values() if n.dofs]
            y_coords = [n.coordinates[1] for n in self.domain.nodes.values() if n.dofs]
            if x_coords and y_coords:
                x_min, x_max = min(x_coords), max(x_coords)
                y_min, y_max = min(y_coords), max(y_coords)
                
                for node in self.domain.nodes.values():
                    if not node.dofs: continue
                    x, y = node.coordinates
                    
                    for loc, bcs in bcs_coord.items():
                        tol = float(bcs.get('tol', 1e-6))
                        match = False
                        if loc == 'x_min' and abs(x - x_min) < tol: match = True
                        elif loc == 'x_max' and abs(x - x_max) < tol: match = True
                        elif loc == 'y_min' and abs(y - y_min) < tol: match = True
                        elif loc == 'y_max' and abs(y - y_max) < tol: match = True
                        elif bcs.get('coord') == 'x' and 'val' in bcs and abs(x - float(bcs['val'])) < tol: match = True
                        elif bcs.get('coord') == 'y' and 'val' in bcs and abs(y - float(bcs['val'])) < tol: match = True
                        
                        if match:
                            for dof, value in bcs.items():
                                if dof not in ['tol', 'coord', 'val']: node.fix_dof(dof, float(value))

        # 5.5 Aplicar Condiciones de Frontera por Grupos Físicos (Gmsh)
        bcs_group = data.get('boundary_conditions_by_group', {})
        if bcs_group and hasattr(self.domain, 'physical_groups'):
            for group_name, bcs in bcs_group.items():
                if group_name in self.domain.physical_groups:
                    for node_id in self.domain.physical_groups[group_name]:
                        node = self.domain.get_node(node_id)
                        if node:
                            for dof, value in bcs.items():
                                if dof not in ['tol', 'coord', 'val']: node.fix_dof(dof, float(value))

        # 6. Guardar configuración de cargas para ensamblarlas después
        self.point_loads = data.get('point_loads', {})
        self.point_loads_by_coord = data.get('point_loads_by_coord', {})
        self.point_loads_by_group = data.get('point_loads_by_group', {})

        self.output_config = data.get('output', {})
        self.solver_config = data.get('solver', {})
        return self.domain
        
    def get_external_forces(self) -> np.ndarray:
        """Construye el vector de fuerzas externas global F_ext."""
        F_ext = np.zeros(self.domain.total_dofs)
        
        for node_id, loads in self.point_loads.items():
            node = self.domain.get_node(node_id)
            for dof, value in loads.items():
                if dof in node.dofs:
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
                    for loc, loads in self.point_loads_by_coord.items():
                        tol = float(loads.get('tol', 1e-6))
                        match = False
                        if loc == 'x_min' and abs(x - x_min) < tol: match = True
                        elif loc == 'x_max' and abs(x - x_max) < tol: match = True
                        elif loc == 'y_min' and abs(y - y_min) < tol: match = True
                        elif loc == 'y_max' and abs(y - y_max) < tol: match = True
                        elif loads.get('coord') == 'x' and 'val' in loads and abs(x - float(loads['val'])) < tol: match = True
                        elif loads.get('coord') == 'y' and 'val' in loads and abs(y - float(loads['val'])) < tol: match = True
                        
                        if match:
                            for dof, value in loads.items():
                                if dof not in ['tol', 'coord', 'val'] and dof in node.dofs:
                                    F_ext[node.dofs[dof]] += float(value)
                                    
        if self.point_loads_by_group and hasattr(self.domain, 'physical_groups'):
            for group_name, loads in self.point_loads_by_group.items():
                if group_name in self.domain.physical_groups:
                    for node_id in self.domain.physical_groups[group_name]:
                        node = self.domain.get_node(node_id)
                        if node:
                            for dof, value in loads.items():
                                if dof not in ['tol', 'coord', 'val'] and dof in node.dofs:
                                    F_ext[node.dofs[dof]] += float(value)
                                    
        return F_ext

    def get_solver(self, assembler):
        """Construye y retorna el solver dinámicamente según la configuración YAML."""
        s_type = self.solver_config.get('type', 'LinearSolver')
        if s_type == 'NonlinearSolver':
            return NonlinearSolver(
                assembler,
                tol=float(self.solver_config.get('tol', 1e-5)),
                max_iter=int(self.solver_config.get('max_iter', 15)),
                num_steps=int(self.solver_config.get('num_steps', 10)),
                adaptive=self.solver_config.get('adaptive', True)
            )
        return LinearSolver(assembler)