# fenix_fem/fenix/utils/gmsh_parser.py
import numpy as np
from fenix.core.domain import Domain
from fenix.elements.solid_2d import Quad4, Tri3
from fenix.core.material import Material

try:
    import meshio
except ImportError:
    raise ImportError("Falta la librería 'meshio'. Instálala con: pip install meshio")

class GmshParser:
    """
    Lector de mallas generadas en Gmsh (.msh) utilizando meshio.
    """
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.domain = Domain()

    def parse(self, default_material: Material, default_thickness: float = 1.0, physical_props: dict = None, default_quadrature: tuple = None) -> Domain:
        """
        Lee el archivo de malla y retorna un objeto Domain poblado.
        Permite mapear Physical Groups a diferentes materiales y espesores.
        """
        print(f"Leyendo malla Gmsh desde: {self.filepath} ...")
        mesh = meshio.read(self.filepath, file_format="gmsh")
        
        # Extraer Grupos Físicos (Physical Groups)
        self.domain.physical_groups = {}
        if hasattr(mesh, 'field_data') and hasattr(mesh, 'cell_data') and "gmsh:physical" in mesh.cell_data:
            tag_to_name = {int(data[0]): name for name, data in mesh.field_data.items()}
            for name in tag_to_name.values():
                self.domain.physical_groups[name] = set()
                
            for i, block in enumerate(mesh.cells):
                tags = mesh.cell_data["gmsh:physical"][i]
                for j, connectivity in enumerate(block.data):
                    tag = int(tags[j])
                    if tag in tag_to_name:
                        name = tag_to_name[tag]
                        for node_idx in connectivity:
                            self.domain.physical_groups[name].add(node_idx + 1)
                            
            for name in self.domain.physical_groups:
                self.domain.physical_groups[name] = list(self.domain.physical_groups[name])

        # 1. Extraer y crear los Nodos
        # mesh.points contiene un arreglo numpy con las coordenadas (x, y, z)
        for i, pt in enumerate(mesh.points):
            node_id = i + 1  # Fenix FEM maneja IDs empezando desde 1
            # Pasamos solo X e Y para análisis 2D
            self.domain.add_node(node_id, [pt[0], pt[1]])
            
        print(f"  -> Nodos importados: {len(self.domain.nodes)}")
            
        # 2. Extraer y crear los Elementos (Filtramos solo cuadriláteros)
        elem_id_counter = 1
        for i, block in enumerate(mesh.cells):
            # Extraer tags fisicos del bloque actual si existen
            tags = mesh.cell_data["gmsh:physical"][i] if ("gmsh:physical" in mesh.cell_data) else None
            
            # Busramos el bloque de elementos tipo "quad" (cuadriláteros de 4 nodos)
            if block.type == "quad":
                for j, connectivity in enumerate(block.data):
                    # meshio devuelve índices base 0, por lo que sumamos 1
                    n1 = self.domain.get_node(connectivity[0] + 1)
                    n2 = self.domain.get_node(connectivity[1] + 1)
                    n3 = self.domain.get_node(connectivity[2] + 1)
                    n4 = self.domain.get_node(connectivity[3] + 1)
                    
                    # Resolver propiedades a asignar
                    mat, thick, quad = default_material, default_thickness, default_quadrature
                    if tags is not None and physical_props is not None:
                        tag = int(tags[j])
                        if tag in tag_to_name:
                            group_name = tag_to_name[tag]
                            if group_name in physical_props:
                                mat, thick, quad = physical_props[group_name]
                                
                    element = Quad4(elem_id_counter, [n1, n2, n3, n4], mat, thick, quadrature=quad)
                    self.domain.add_element(element)
                    elem_id_counter += 1
                    
            elif block.type == "triangle":
                for j, connectivity in enumerate(block.data):
                    n1 = self.domain.get_node(connectivity[0] + 1)
                    n2 = self.domain.get_node(connectivity[1] + 1)
                    n3 = self.domain.get_node(connectivity[2] + 1)
                    
                    mat, thick, quad = default_material, default_thickness, default_quadrature
                    if tags is not None and physical_props is not None:
                        tag = int(tags[j])
                        if tag in tag_to_name:
                            group_name = tag_to_name[tag]
                            if group_name in physical_props:
                                mat, thick, _ = physical_props[group_name]
                                
                    element = Tri3(elem_id_counter, [n1, n2, n3], mat, thick)
                    self.domain.add_element(element)
                    elem_id_counter += 1
                    
        if len(self.domain.elements) == 0:
            raise RuntimeError("Error crítico: No se importó ningún elemento. Revisa que tu malla tenga superficies.")
            
        print(f"  -> Elementos (Solid2D) importados: {len(self.domain.elements)}")
        return self.domain