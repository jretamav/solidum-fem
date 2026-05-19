# fenix_fem/fenix/utils/yaml_parser.py
import inspect
import re
import yaml
import os
import numpy as np
from fenix.core.domain import Domain
from fenix.autodiscover import initialize as _ensure_registries_initialized
from fenix.logging import get_logger
from fenix.registry import (
    CohesiveMaterialRegistry,
    ElementRegistry,
    MaterialRegistry,
    QuadratureRegistry,
    SolverRegistry,
)

_log = get_logger("parsers.yaml")

# Idempotente: garantiza que los registries estén poblados aunque se importe
# este módulo sin pasar por `fenix/__init__.py` (p. ej., en tests aislados).
_ensure_registries_initialized()


# PyYAML usa la spec YAML 1.1, que exige signo explícito en el exponente
# (`2.0e+9` sí, `2.0e9` no — se interpreta como string). Registramos el
# regex de floats de YAML 1.2 sobre SafeLoader para aceptar ambos.
_FLOAT_RESOLVER = re.compile(
    r"""^(?:
         [-+]?(?:[0-9][0-9_]*)\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\.[0-9][0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?\.(?:inf|Inf|INF)
        |\.(?:nan|NaN|NAN)
        )$""",
    re.VERBOSE,
)
yaml.SafeLoader.add_implicit_resolver(
    'tag:yaml.org,2002:float', _FLOAT_RESOLVER, list('-+0123456789.')
)


class YamlValidationError(Exception):
    """Error de validación del archivo YAML de entrada.

    Acumula todos los problemas encontrados y los presenta juntos
    para que el usuario pueda corregirlos en una sola pasada.
    """
    def __init__(self, errors: list):
        self.errors = errors
        lines = [f"  [{i+1}] {e}" for i, e in enumerate(errors)]
        super().__init__("\n\nEl archivo YAML contiene los siguientes errores:\n" + "\n".join(lines))


class YamlParser:
    """Lector automatizado de modelos estructurales desde archivos YAML."""
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.domain = Domain()
        self.materials = {}
        self.cohesive_materials = {}    # ADR 0010 — familia paralela
        self.solver_config = {}
        self.point_loads = {}
        self.point_loads_by_coord = {}
        self.point_loads_by_group = {}
        self.point_loads_by_node = []
        self.body_force = None  # np.ndarray (2,) o (3,) o None si no se especificó
        self.gravity = None     # np.ndarray (2,) o (3,) o None (ADR 0008)
        self.output_config = {}

    def _get_quadrature(self, rule: str) -> tuple:
        """Convierte un string del YAML a una regla de cuadratura inyectable."""
        if rule == "1x1":
            _log.warning("Integración reducida (1x1) seleccionada. Cuidado con modos de energía nula (hourglass).")

        return QuadratureRegistry.get(rule)

    def parse(self) -> Domain:
        _log.info(f"Leyendo modelo desde: {self.filepath} ...")
        with open(self.filepath, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict) or not data:
            raise YamlValidationError(["El archivo YAML está vacío o no tiene formato de mapa clave-valor."])

        errors = self._validate(data)
        if errors:
            raise YamlValidationError(errors)

        self._parse_nodes(data)
        self._parse_materials(data)
        self._parse_cohesive_materials(data)
        self._parse_mesh_or_elements(data)
        self._parse_boundary_conditions(data)
        self._parse_linear_constraints(data)
        self._parse_point_loads_config(data)
        self._parse_body_force(data)

        self.output_config = data.get('output', {})
        self.solver_config = data.get('solver', {})
        return self.domain

    def _validate(self, data: dict) -> list:
        """Valida la estructura del YAML y retorna una lista de errores encontrados.

        No lanza excepción: acumula todos los problemas para reportarlos juntos.
        """
        errors = []

        has_mesh = bool(data.get('mesh'))
        has_nodes = bool(data.get('nodes'))
        has_elements = bool(data.get('elements'))

        # --- Estructura top-level ---
        if not has_mesh and not has_nodes:
            errors.append("Falta el bloque 'nodes' (o 'mesh'). El modelo no tiene geometría definida.")

        if not has_mesh and has_elements and not data.get('materials'):
            errors.append("Falta el bloque 'materials'. Se definieron elementos pero no hay materiales.")

        # --- Nodos ---
        known_node_ids = set()
        nodes_data = data.get('nodes', [])
        if isinstance(nodes_data, list):
            for i, node in enumerate(nodes_data):
                ctx = f"nodes[{i}]"
                if not isinstance(node, dict):
                    errors.append(f"{ctx}: cada nodo debe ser un diccionario, se encontró {type(node).__name__}.")
                    continue
                if 'id' not in node:
                    errors.append(f"{ctx}: falta el campo obligatorio 'id'.")
                else:
                    nid = node['id']
                    if nid in known_node_ids:
                        errors.append(f"{ctx}: id de nodo duplicado ({nid}).")
                    known_node_ids.add(nid)
                if 'coords' not in node:
                    errors.append(f"{ctx} (id={node.get('id', '?')}): falta el campo obligatorio 'coords'.")
                else:
                    coords = node['coords']
                    if not isinstance(coords, (list, tuple)) or len(coords) not in (2, 3):
                        errors.append(f"{ctx} (id={node.get('id', '?')}): 'coords' debe ser una lista de 2 o 3 números, se encontró: {coords!r}.")
        elif nodes_data:
            errors.append("El bloque 'nodes' debe ser una lista de diccionarios.")

        # --- Materiales ---
        known_mat_ids = set()
        registered_materials = set(MaterialRegistry._materials.keys())
        for i, mat in enumerate(data.get('materials', [])):
            ctx = f"materials[{i}]"
            if not isinstance(mat, dict):
                errors.append(f"{ctx}: cada material debe ser un diccionario.")
                continue
            if 'id' not in mat:
                errors.append(f"{ctx}: falta el campo obligatorio 'id'.")
            else:
                known_mat_ids.add(mat['id'])
            if 'type' not in mat:
                errors.append(f"{ctx} (id={mat.get('id', '?')}): falta el campo obligatorio 'type'.")
            elif registered_materials and mat['type'] not in registered_materials:
                errors.append(
                    f"{ctx} (id={mat.get('id', '?')}): tipo de material desconocido '{mat['type']}'. "
                    f"Disponibles: {sorted(registered_materials)}."
                )

        # --- Materiales cohesivos (ADR 0010 — sección paralela) ---
        known_cohesive_ids = set()
        registered_cohesives = set(CohesiveMaterialRegistry._items.keys())
        for i, mat in enumerate(data.get('cohesive_materials', []) or []):
            ctx = f"cohesive_materials[{i}]"
            if not isinstance(mat, dict):
                errors.append(f"{ctx}: cada material cohesivo debe ser un diccionario.")
                continue
            if 'id' not in mat:
                errors.append(f"{ctx}: falta el campo obligatorio 'id'.")
            else:
                known_cohesive_ids.add(mat['id'])
            if 'type' not in mat:
                errors.append(f"{ctx} (id={mat.get('id', '?')}): falta el campo obligatorio 'type'.")
            elif registered_cohesives and mat['type'] not in registered_cohesives:
                errors.append(
                    f"{ctx} (id={mat.get('id', '?')}): tipo de cohesivo desconocido "
                    f"'{mat['type']}'. Disponibles: {sorted(registered_cohesives)}."
                )

        # --- Elementos (bloque inline, no mesh) ---
        registered_elements = set(ElementRegistry._elements.keys())
        elements_data = data.get('elements', [])
        if isinstance(elements_data, list):
            known_elem_ids = set()
            for i, elem in enumerate(elements_data):
                ctx = f"elements[{i}]"
                if not isinstance(elem, dict):
                    errors.append(f"{ctx}: cada elemento debe ser un diccionario.")
                    continue
                eid = elem.get('id')
                if eid is None:
                    errors.append(f"{ctx}: falta el campo obligatorio 'id'.")
                elif eid in known_elem_ids:
                    errors.append(f"{ctx}: id de elemento duplicado ({eid}).")
                else:
                    known_elem_ids.add(eid)

                ctx = f"elements[{i}] (id={eid or '?'})"
                if 'type' not in elem:
                    errors.append(f"{ctx}: falta el campo obligatorio 'type'.")
                elif registered_elements and elem['type'] not in registered_elements:
                    errors.append(
                        f"{ctx}: tipo de elemento desconocido '{elem['type']}'. "
                        f"Disponibles: {sorted(registered_elements)}."
                    )
                if 'material' not in elem:
                    errors.append(f"{ctx}: falta el campo obligatorio 'material'.")
                elif elem['material'] not in known_mat_ids:
                    errors.append(f"{ctx}: referencia a material inexistente (id={elem['material']}).")
                # ADR 0010 — referencia a cohesivo opcional, sólo si el elemento
                # la admite. Si se declara, validar contra `cohesive_materials`.
                if 'cohesive_material' in elem and elem['cohesive_material'] not in known_cohesive_ids:
                    errors.append(
                        f"{ctx}: referencia a cohesive_material inexistente "
                        f"(id={elem['cohesive_material']})."
                    )
                if 'nodes' not in elem:
                    errors.append(f"{ctx}: falta el campo obligatorio 'nodes'.")
                else:
                    node_refs = elem['nodes']
                    if not isinstance(node_refs, list) or len(node_refs) < 2:
                        errors.append(f"{ctx}: 'nodes' debe ser una lista con al menos 2 ids de nodo.")
                    elif known_node_ids:
                        for nref in node_refs:
                            if nref not in known_node_ids:
                                errors.append(f"{ctx}: referencia a nodo inexistente (id={nref}).")

                # Validación de kwargs: los campos extra del YAML deben existir
                # en la firma del constructor del elemento registrado.
                e_type = elem.get('type')
                if e_type in ElementRegistry._elements:
                    cls = ElementRegistry._elements[e_type]
                    sig = inspect.signature(cls.__init__)
                    accepted = {
                        name for name, p in sig.parameters.items()
                        if name != 'self' and p.kind in (
                            inspect.Parameter.POSITIONAL_OR_KEYWORD,
                            inspect.Parameter.KEYWORD_ONLY,
                        )
                    }
                    has_var_keyword = any(
                        p.kind == inspect.Parameter.VAR_KEYWORD
                        for p in sig.parameters.values()
                    )
                    if not has_var_keyword:
                        reserved = {'id', 'type', 'material', 'nodes', 'cohesive_material'}
                        extras = set(elem.keys()) - reserved
                        # 'element_id' lo inyecta el parser a partir de 'id'
                        unknown = extras - accepted
                        # Mensaje: excluir kwargs que el parser inyecta ('element_id')
                        # y los reservados YAML ('material', 'nodes', 'cohesive_material')
                        # que el usuario provee por campos dedicados.
                        advertised = sorted(accepted - {'element_id'} - reserved)
                        for kw in sorted(unknown):
                            errors.append(
                                f"{ctx}: parámetro '{kw}' no aceptado por '{e_type}'. "
                                f"Admitidos: {advertised}."
                            )

        # --- Boundary conditions (por nodo) ---
        for bc in (data.get('boundary_conditions', []) or []) + (data.get('boundary_conditions_by_node', []) or []):
            if not isinstance(bc, dict):
                continue
            nid = bc.get('node_id')
            if nid is None:
                errors.append("boundary_conditions: una entrada no tiene 'node_id'.")
            elif known_node_ids and nid not in known_node_ids:
                errors.append(f"boundary_conditions: 'node_id={nid}' no existe en el bloque 'nodes'.")

        # --- Restricciones afines lineales (MPC, ADR 0004 fase 2) ---
        for i, lc in enumerate(data.get('linear_constraints', []) or []):
            ctx = f"linear_constraints[{i}]"
            if not isinstance(lc, dict):
                errors.append(f"{ctx}: cada entrada debe ser un diccionario.")
                continue
            slave = lc.get('slave')
            if not isinstance(slave, dict) or 'node' not in slave or 'dof' not in slave:
                errors.append(f"{ctx}: 'slave' debe ser {{node: <id>, dof: <name>}}.")
            elif known_node_ids and slave['node'] not in known_node_ids:
                errors.append(f"{ctx}: 'slave.node={slave['node']}' no existe en el bloque 'nodes'.")
            masters = lc.get('masters', [])
            coefficients = lc.get('coefficients', [])
            if not isinstance(masters, list) or not isinstance(coefficients, list):
                errors.append(f"{ctx}: 'masters' y 'coefficients' deben ser listas.")
            elif len(masters) != len(coefficients):
                errors.append(
                    f"{ctx}: 'masters' ({len(masters)}) y 'coefficients' "
                    f"({len(coefficients)}) deben tener la misma longitud."
                )
            else:
                for j, m in enumerate(masters):
                    if not isinstance(m, dict) or 'node' not in m or 'dof' not in m:
                        errors.append(f"{ctx}.masters[{j}]: debe ser {{node: <id>, dof: <name>}}.")
                    elif known_node_ids and m['node'] not in known_node_ids:
                        errors.append(f"{ctx}.masters[{j}]: 'node={m['node']}' no existe en el bloque 'nodes'.")

        # --- Solver ---
        solver_cfg = data.get('solver', {})
        if solver_cfg:
            s_type = solver_cfg.get('type')
            if not s_type:
                errors.append("solver: falta el campo 'type'.")
            else:
                registered_solvers = set(SolverRegistry._solvers.keys())
                if registered_solvers and s_type not in registered_solvers:
                    errors.append(
                        f"solver: tipo desconocido '{s_type}'. "
                        f"Disponibles: {sorted(registered_solvers)}."
                    )

            # Override del backend algebraico (ADR 0003 §4) — diagnóstico, no
            # decisión de modelado. Validamos contra el registro del despachador.
            la = solver_cfg.get('linear_algebra')
            if la is not None:
                from fenix.math.linalg.dispatcher import _REGISTRY as _LA_REGISTRY
                if la != 'auto' and la not in _LA_REGISTRY:
                    errors.append(
                        f"solver: 'linear_algebra' desconocido: '{la}'. "
                        f"Disponibles: {['auto'] + sorted(_LA_REGISTRY)}."
                    )

        return errors

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

    def _parse_cohesive_materials(self, data: dict):
        """Materiales cohesivos *traction-jump* (ADR 0010, sección paralela a
        ``materials``). Se construyen con :class:`CohesiveMaterialRegistry`
        para que el contrato del parser YAML no se mezcle con los continuos.
        """
        for mat_data in data.get('cohesive_materials', []) or []:
            mat_id = mat_data['id']
            mat_type = mat_data['type']
            kwargs = {k: v for k, v in mat_data.items() if k not in ('id', 'type')}
            self.cohesive_materials[mat_id] = CohesiveMaterialRegistry.create(mat_type, **kwargs)

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
                    
                    kwargs = {k: v for k, v in elem_dict.items() if k not in ('id', 'type', 'material', 'nodes', 'cohesive_material')}
                    if 'quadrature' in kwargs and isinstance(kwargs['quadrature'], str):
                        kwargs['quadrature'] = self._get_quadrature(kwargs['quadrature'])
                    # ADR 0010 — resolver referencia a material cohesivo si está declarada.
                    if 'cohesive_material' in elem_dict:
                        kwargs['cohesive_material'] = self.cohesive_materials[elem_dict['cohesive_material']]
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

    def _parse_linear_constraints(self, data: dict):
        """Restricciones afines lineales MPC (ADR 0004 fase 2).

        Sintaxis YAML::

            linear_constraints:
              - slave:  {node: 3, dof: uy}
                masters:
                  - {node: 2, dof: uy}
                  - {node: 2, dof: rz}
                coefficients: [1.0, 0.5]
                g: 0.0          # opcional, default 0.0

        Casos de uso típicos: apoyo en plano oblicuo, periodicidad, unión
        rígida master-slave, simetrías no alineadas con ejes globales.
        """
        for lc in data.get('linear_constraints', []) or []:
            slave = lc['slave']
            masters = [(m['node'], m['dof']) for m in lc.get('masters', [])]
            coefficients = list(lc.get('coefficients', []))
            g = float(lc.get('g', 0.0))
            self.domain.add_linear_constraint(
                slave=(slave['node'], slave['dof']),
                masters=masters,
                coefficients=coefficients,
                g=g,
            )

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

    def _parse_body_force(self, data: dict):
        """Lee `body_force` y/o `gravity` del YAML (ADR 0008).

        Formas aceptadas — exclusivas entre sí::

            # Forma 1: peso propio físicamente correcto vía densidades por
            # material (recomendada para análisis estructural):
            gravity: [0.0, -9.81]            # 2D
            gravity: [0.0, 0.0, -9.81]       # 3D

            # Forma 2: fuerza de cuerpo uniforme arbitraria (útil cuando no
            # es peso propio, o para modelos monomaterial donde el usuario
            # ya precalculó ρ·g):
            body_force: [0.0, -78.5e3]       # 2D, N/m³ (acero)
            body_force: [0.0, 0.0, -78.5e3]  # 3D

        Declarar ambas en el mismo archivo es ambigüedad y lanza error.
        """
        bf = data.get('body_force')
        gv = data.get('gravity')

        if bf is not None and gv is not None:
            raise YamlValidationError([
                "No se pueden declarar 'body_force' y 'gravity' simultáneamente. "
                "Usar 'gravity' (con 'density' por material) para peso propio físicamente correcto; "
                "usar 'body_force' solo para cargas de cuerpo que no son peso propio."
            ])

        self.body_force = self._coerce_vector(bf, 'body_force') if bf is not None else None
        self.gravity = self._coerce_vector(gv, 'gravity') if gv is not None else None

    @staticmethod
    def _coerce_vector(value, name: str) -> np.ndarray:
        """Convierte una secuencia YAML a np.ndarray (2,) ó (3,)."""
        try:
            arr = np.asarray(value, dtype=float).ravel()
        except (TypeError, ValueError) as exc:
            raise YamlValidationError(
                [f"{name}: no se pudo interpretar como vector numérico: {value!r} ({exc})."]
            )
        if arr.size not in (2, 3):
            raise YamlValidationError(
                [f"{name}: longitud {arr.size}, esperado 2 ó 3 componentes."]
            )
        return arr
        
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

    def get_body_load(self, assembler) -> np.ndarray:
        """Vector global de fuerzas nodales por fuerza de cuerpo o peso propio.

        Despacha según el bloque presente en el YAML (ADR 0008):

        - ``gravity`` ⇒ ``Assembler.assemble_self_weight(g)`` — peso propio
          físicamente correcto usando ``material.density`` por elemento.
        - ``body_force`` ⇒ ``Assembler.assemble_body_load(b)`` — vector
          uniforme aplicado a todos los elementos.
        - Ninguno ⇒ vector de ceros.

        La exclusividad entre los dos bloques se valida en :meth:`_parse_body_force`.
        """
        if self.gravity is not None:
            return assembler.assemble_self_weight(self.gravity)
        if self.body_force is not None:
            return assembler.assemble_body_load(self.body_force)
        return np.zeros(self.domain.total_dofs)

    def get_solver(self, assembler):
        """Construye y retorna el solver dinámicamente según la configuración YAML."""
        s_type = self.solver_config.get('type', 'LinearSolver')
        kwargs = {k: v for k, v in self.solver_config.items() if k != 'type'}

        # Bloque `convergence:` (ADR 0007). Si está presente lo materializamos
        # como ConvergenceCriterion antes de pasar al constructor del solver.
        # Los solvers lineales (estáticos y modales) no lo admiten — se omite
        # silenciosamente. El análisis modal no es iterativo en el sentido
        # de Newton; ARPACK gestiona su propia tolerancia interna.
        if 'convergence' in kwargs:
            from fenix.math.convergence import make_convergence_from_config
            cfg = kwargs.pop('convergence')
            # Solvers lineales (sin Newton interno) no admiten `convergence`:
            # estáticos lineales, modal (ARPACK gestiona su propia tolerancia),
            # transitorios lineales (Newmark y HHT) y armónico (lineal por
            # construcción).
            if s_type not in ('LinearSolver', 'ModalSolver',
                              'NewmarkSolver', 'HHTSolver',
                              'CentralDifferenceSolver',
                              'HarmonicSolver',
                              'ResponseSpectrumSolver'):
                kwargs['convergence'] = make_convergence_from_config(cfg)

        # HarmonicSolver: si no se pasó `F_amplitude` explícito, derivar la
        # amplitud compleja del bloque estándar de cargas externas del YAML
        # (point_loads / etc.). Las cargas se interpretan como amplitudes
        # reales con fase cero; para cargas con fase, el usuario pasa
        # `F_amplitude` numpy complejo directamente desde código (YAML no
        # soporta tipos complejos nativos).
        if s_type == 'HarmonicSolver' and 'F_amplitude' not in kwargs:
            kwargs['F_amplitude'] = self.get_external_forces()

        return SolverRegistry.create(s_type, assembler=assembler, **kwargs)