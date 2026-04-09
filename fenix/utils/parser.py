# fenix_fem/fenix/utils/parser.py
import shlex
from fenix.core.domain import Domain
from fenix.materials.elastic import Elastic1D
from fenix.elements.structural import Truss2D
from fenix.math.assembly import Assembler

class FEAParser:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.mesh = Domain()
        self.assembler = Assembler(self.mesh)
        self.parsed_forces = []
        self.materials_dict = {}
        self.raw_materials = {}
        self.current_mat_id = None
        self.materials_built = False

    def parse(self):
        with open(self.filepath, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        current_block = None
        for line_number, raw_line in enumerate(lines):
            line = raw_line.split('!')[0].split('#')[0].strip()
            if not line or line.upper() in ['END', 'BATC', 'STOP', 'INTE']: continue
            
            # Normalizar MATErial,1
            tokens = shlex.split(line.replace(',', ' '))
            keyword = tokens[0].upper()
            
            if keyword in ['COORDINATES', 'COOR', 'ELEMENTS', 'ELEM', 'BOUNDARY', 'BOUN', 'FORCES', 'FORC']:
                if keyword.startswith('ELEM') and not self.materials_built:
                    self._build_materials()
                    self.materials_built = True
                current_block = keyword
                continue
                
            if keyword.startswith('MATE'):
                current_block = 'MATERIAL'
                self.current_mat_id = int(tokens[1])
                self.raw_materials[self.current_mat_id] = {'type': None, 'E': 0.0, 'area': 0.0}
                continue

            try:
                if current_block == 'MATERIAL':
                    self._parse_material_line(tokens)
                elif current_block in ['COORDINATES', 'COOR']:
                    # Formato limpio FENIX: ID X Y
                    self.mesh.add_node(int(tokens[0]), [float(tokens[1]), float(tokens[2])])
                elif current_block in ['ELEMENTS', 'ELEM']:
                    # Formato limpio FENIX: ID MAT 1 NODES 1 2
                    mat_id = int(tokens[tokens.index('MAT')+1])
                    nodes_idx = tokens.index('NODES')
                    nodes = [self.mesh.get_node(int(n)) for n in tokens[nodes_idx+1:]]
                    mat_data = self.materials_dict[mat_id]
                    self.mesh.add_element(Truss2D(int(tokens[0]), nodes, mat_data['obj'], A=mat_data['area']))
                elif current_block in ['BOUNDARY', 'BOUN']:
                    # Formato limpio FENIX: ID RestX RestY
                    node = self.mesh.get_node(int(tokens[0]))
                    if int(tokens[1]) == 1: node.fix_dof('ux', 0.0)
                    if int(tokens[2]) == 1: node.fix_dof('uy', 0.0)
                elif current_block in ['FORCES', 'FORC']:
                    # Formato limpio FENIX: ID Fx Fy
                    node_id = int(tokens[0])
                    if float(tokens[1]) != 0.0: self.parsed_forces.append((node_id, 'ux', float(tokens[1])))
                    if float(tokens[2]) != 0.0: self.parsed_forces.append((node_id, 'uy', float(tokens[2])))
            except Exception as e:
                raise SyntaxError(f"Error parseando linea {line_number + 1}: {e}")
                
        return self.mesh, self.assembler

    def _parse_material_line(self, tokens):
        key = tokens[0].upper()
        # Aquí SÍ aplicamos el formato estricto FEAP
        if key.startswith('TRUS'):
            self.raw_materials[self.current_mat_id]['type'] = 'TRUSS'
        elif key.startswith('ELAS'):
            self.raw_materials[self.current_mat_id]['E'] = float(tokens[2])
        elif key.startswith('CROS'):
            self.raw_materials[self.current_mat_id]['area'] = float(tokens[2])

    def _build_materials(self):
        for mid, data in self.raw_materials.items():
            self.materials_dict[mid] = {'obj': Elastic1D(E=data['E']), 'area': data['area']}

    def apply_parsed_forces(self):
        for node_id, dof_name, value in self.parsed_forces:
            self.assembler.apply_point_load(node_id, dof_name, value)
