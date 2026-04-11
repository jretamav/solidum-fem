# fenix_fem/fenix/utils/yaml_parser.py
import yaml
import os
import numpy as np
from fenix.core.domain import Domain
from fenix.materials.von_mises_2d import VonMises2D
from fenix.materials.elastic_2d import Elastic2D
from fenix.materials.elastic import Elastic1D
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.materials.damage_1d import IsotropicDamage1D
from fenix.elements.solid_2d import Quad4, Tri3
from fenix.elements.structural import Frame2DEuler, Frame2DTimoshenko, Truss2D, Truss3D
from fenix.math.solvers import LinearSolver, NonlinearSolver
from fenix.math.integration import GaussQuadrature

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
            return GaussQuadrature.get_points_2d_1x1()
        # Por defecto 2x2
        return GaussQuadrature.get_points_2d_2x2()

    def parse(self) -> Domain:
        print(f"Leyendo modelo desde: {self.filepath} ...")
        with open(self.filepath, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)

        # 1. Instanciar Nodos
        nodes_data = data.get('nodes', [])
        if isinstance(nodes_data, list):
            for node_dict in nodes_data:
                self.domain.add_node(node_dict['id'], node_dict['coords'])
        else:
            raise ValueError("El bloque 'nodes' debe ser una lista de diccionarios con formato moderno.")

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
            elif mat_type == 'Elastic1D':
                self.materials[mat_id] = Elastic1D(E=float(mat_data['E']))
            elif mat_type == 'Elastoplastic1D':
                self.materials[mat_id] = Elastoplastic1D(
                    E=float(mat_data['E']),
                    sigma_y=float(mat_data['sigma_y']),
                    H=float(mat_data.get('H', 0.0))
                )
            elif mat_type == 'IsotropicDamage2D':
                self.materials[mat_id] = IsotropicDamage2D(
                    E=float(mat_data['E']),
                    nu=float(mat_data['nu']),
                    kappa_0=float(mat_data['kappa_0']),
                    alpha=float(mat_data['alpha']),
                    hypothesis=mat_data.get('hypothesis', 'plane_stress')
                )
            elif mat_type == 'IsotropicDamage1D':
                self.materials[mat_id] = IsotropicDamage1D(
                    E=float(mat_data['E']),
                    kappa_0=float(mat_data['kappa_0']),
                    alpha=float(mat_data['alpha'])
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
            
            default_mat_id = data.get('mesh_material', 1)
            default_thickness = float(data.get('mesh_thickness', 1.0))
            default_quad_str = data.get('mesh_quadrature', '2x2')
            
            default_material = self.materials.get(default_mat_id)
            default_quadrature = self._get_quadrature(default_quad_str)
            
            # Extraer mapeo explícito de materiales a grupos físicos si existe
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
                # Formato moderno (Lista de diccionarios)
                for elem_dict in elements_data:
                    elem_id = elem_dict['id']
                    e_type = elem_dict['type']
                    mat_id = elem_dict['material']
                    node_ids = elem_dict['nodes']
                    
                    nodes = [self.domain.get_node(nid) for nid in node_ids]
                    material = self.materials[mat_id]
                    
                    if e_type == 'Quad4':
                        thickness = float(elem_dict.get('thickness', 1.0))
                        self.domain.add_element(Quad4(elem_id, nodes, material, thickness))
                    elif e_type == 'Tri3':
                        thickness = float(elem_dict.get('thickness', 1.0))
                        self.domain.add_element(Tri3(elem_id, nodes, material, thickness))
                    elif e_type == 'Frame2DEuler':
                        A = float(elem_dict.get('A', 0.0))
                        I = float(elem_dict.get('I', 0.0))
                        self.domain.add_element(Frame2DEuler(elem_id, nodes, material, A=A, I=I))
                    elif e_type == 'Frame2DTimoshenko':
                        A = float(elem_dict.get('A', 0.0))
                        I = float(elem_dict.get('I', 0.0))
                        As = float(elem_dict.get('As', 0.0))
                        self.domain.add_element(Frame2DTimoshenko(elem_id, nodes, material, A=A, I=I, As=As))
                    elif e_type == 'Truss2D':
                        A = float(elem_dict.get('A', 0.0))
                        self.domain.add_element(Truss2D(elem_id, nodes, material, A=A))
                    elif e_type == 'Truss3D':
                        A = float(elem_dict.get('A', 0.0))
                        self.domain.add_element(Truss3D(elem_id, nodes, material, A=A))
                    else:
                        raise ValueError(f"Tipo de elemento {e_type} no soportado.")
            elif elements_data:
                raise ValueError("El bloque 'elements' debe ser una lista de diccionarios con formato moderno.")

        # 4. Aplicar Condiciones de Frontera (Desplazamientos)
        for node_id, bcs in data.get('boundary_conditions', {}).items():
            node = self.domain.get_node(node_id)
            for dof, value in bcs.items():
                node.fix_dof(dof, float(value))
                
        for bc in data.get('boundary_conditions_by_node', []):
            node_id = bc['node_id']
            node = self.domain.get_node(node_id)
            for dof, value in bc.items():
                if dof != 'node_id':
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
        self.point_loads_by_node = data.get('point_loads_by_node', [])

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
                    
        for load in self.point_loads_by_node:
            node_id = load['node_id']
            node = self.domain.get_node(node_id)
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