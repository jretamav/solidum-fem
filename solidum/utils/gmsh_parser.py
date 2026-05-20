# fenix_fem/solidum/utils/gmsh_parser.py
from solidum.core.domain import Domain
from solidum.logging import get_logger
from solidum.registry import ElementRegistry
from solidum.core.material import Material

try:
    import meshio
except ImportError:
    raise ImportError("Falta la librería 'meshio'. Instálala con: pip install meshio")

_log = get_logger("parsers.gmsh")

# Mapeo de tipos de elemento de Gmsh/meshio a los nombres registrados en Solidum FEM.
# Esto permite extender el soporte a nuevos elementos (ej. 'triangle6': 'Tri6')
# sin modificar la lógica del parser.
GMSH_TYPE_MAP = {
    "quad": ("Quad4", [0, 1, 2, 3]),
    "triangle": ("Tri3", [0, 1, 2]),
}

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
        _log.info(f"Leyendo malla Gmsh desde: {self.filepath} ...")
        mesh = meshio.read(self.filepath, file_format="gmsh")
        
        # Extraer Grupos Físicos (Physical Groups). Inicializamos `tag_to_name`
        # fuera del condicional: si la malla no trae `field_data` pero sí
        # `cell_data["gmsh:physical"]` (caso raro pero plausible), el bucle
        # de elementos seguía usándola más abajo y disparaba `UnboundLocalError`.
        self.domain.physical_groups = {}
        tag_to_name: dict[int, str] = {}
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
            node_id = i + 1  # Solidum FEM maneja IDs empezando desde 1
            # Pasamos solo X e Y para análisis 2D
            self.domain.add_node(node_id, [pt[0], pt[1]])
            
        _log.info(f"  -> Nodos importados: {len(self.domain.nodes)}")
            
        # 2. Extraer y crear los Elementos (Filtramos solo cuadriláteros)
        elem_id_counter = 1
        for i, block in enumerate(mesh.cells):
            # Extraer tags fisicos del bloque actual si existen
            tags = mesh.cell_data["gmsh:physical"][i] if ("gmsh:physical" in mesh.cell_data) else None
            
            if block.type in GMSH_TYPE_MAP:
                fenix_name, node_map_indices = GMSH_TYPE_MAP[block.type]

                for j, connectivity in enumerate(block.data):
                    # Mapeo de nodos basado en GMSH_TYPE_MAP
                    nodes = [self.domain.get_node(connectivity[k] + 1) for k in node_map_indices]
                    
                    # Resolver propiedades a asignar
                    mat, thick, quad = default_material, default_thickness, default_quadrature
                    if tags is not None and physical_props is not None:
                        tag = int(tags[j])
                        if tag in tag_to_name:
                            group_name = tag_to_name[tag]
                            if group_name in physical_props:
                                # Tri3 no usa cuadratura, así que la ignoramos si no está
                                mat_prop, thick_prop, quad_prop = physical_props[group_name]
                                mat = mat_prop
                                thick = thick_prop
                                if fenix_name != 'Tri3':
                                    quad = quad_prop

                    # Argumentos para el constructor del elemento
                    elem_args = {
                        'element_id': elem_id_counter,
                        'nodes': nodes,
                        'material': mat,
                        'thickness': thick,
                    }
                    if fenix_name != 'Tri3':
                        elem_args['quadrature'] = quad

                    element = ElementRegistry.create(
                        fenix_name, **elem_args
                    )
                    self.domain.add_element(element)
                    elem_id_counter += 1
                    
        if len(self.domain.elements) == 0:
            raise RuntimeError("Error crítico: No se importó ningún elemento. Revisa que tu malla tenga superficies.")
            
        _log.info(f"  -> Elementos (Solid2D) importados: {len(self.domain.elements)}")
        return self.domain