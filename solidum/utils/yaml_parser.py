# solidum_fem/solidum/utils/yaml_parser.py
import inspect
import re
import yaml
import os
import numpy as np
from solidum.core.domain import Domain
from solidum.autodiscover import initialize as _ensure_registries_initialized
from solidum.logging import get_logger
from solidum.registry import (
    CohesiveMaterialRegistry,
    ElementRegistry,
    MaterialRegistry,
    QuadratureRegistry,
    SolverRegistry,
    ThermalMaterialRegistry,
)

_log = get_logger("parsers.yaml")

# Idempotente: garantiza que los registries estén poblados aunque se importe
# este módulo sin pasar por `solidum/__init__.py` (p. ej., en tests aislados).
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


class _SolidumLoader(yaml.SafeLoader):
    """``SafeLoader`` con floats YAML 1.2. Subclase propia para no mutar el
    ``SafeLoader`` global de PyYAML, que compartiría cualquier otra librería
    del proceso host (auditoría 2026-09-22)."""


_SolidumLoader.add_implicit_resolver(
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


def _constructor_kwargs(cls) -> tuple[set, bool]:
    """``(nombres aceptados, acepta **kwargs)`` del constructor de ``cls``."""
    sig = inspect.signature(cls.__init__)
    accepted = {
        name for name, p in sig.parameters.items()
        if name != 'self' and p.kind in (
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
            inspect.Parameter.KEYWORD_ONLY,
        )
    }
    has_var_keyword = any(
        p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()
    )
    return accepted, has_var_keyword


def _unknown_kwargs(cls, given: set, reserved: set) -> tuple[list, list]:
    """Claves de ``given`` que el constructor de ``cls`` no acepta, y la lista
    de admitidas para el mensaje. Vacío si el constructor toma ``**kwargs``."""
    accepted, var_kw = _constructor_kwargs(cls)
    if var_kw:
        return [], sorted(accepted - reserved)
    unknown = sorted((given - reserved) - accepted)
    return unknown, sorted(accepted - reserved)


class YamlParser:
    """Lector automatizado de modelos estructurales desde archivos YAML."""

    # Claves de control de una entrada ``*_by_coord`` (el resto son DOFs).
    _COORD_KEYS = ('tol', 'coord', 'val', 'loc')
    _COORD_LOCS = ('x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max')
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.domain = Domain()
        self.materials = {}
        self.cohesive_materials = {}    # ADR 0010 — familia paralela
        self.thermal_materials = {}     # Etapa 8 — familia paralela
        self.thermal_body_sources = []  # Etapa 8 — fuente volumétrica Q
        self.thermal_boundary_fluxes = []  # Etapa 8 — flujo prescrito q̄
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
            data = yaml.load(f, Loader=_SolidumLoader)

        if not isinstance(data, dict) or not data:
            raise YamlValidationError(["El archivo YAML está vacío o no tiene formato de mapa clave-valor."])

        errors = self._validate(data)
        if errors:
            raise YamlValidationError(errors)

        self._parse_nodes(data)
        self._parse_materials(data)
        self._parse_cohesive_materials(data)
        self._parse_thermal_materials(data)
        self._parse_mesh_or_elements(data)
        self._parse_boundary_conditions(data)
        self._parse_linear_constraints(data)
        self._parse_point_loads_config(data)
        self._parse_thermal_loads_config(data)
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

        if (not has_mesh and has_elements
                and not data.get('materials')
                and not data.get('thermal_materials')):
            errors.append(
                "Falta el bloque 'materials' (o 'thermal_materials' para un "
                "modelo térmico). Se definieron elementos pero no hay materiales."
            )

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
        registered_materials = set(MaterialRegistry.names())
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
            else:
                unknown, advertised = _unknown_kwargs(
                    MaterialRegistry.get(mat['type']), set(mat), {'id', 'type'},
                )
                for kw in unknown:
                    errors.append(
                        f"{ctx} (id={mat.get('id', '?')}): parámetro '{kw}' no aceptado "
                        f"por '{mat['type']}'. Admitidos: {advertised}."
                    )

        # --- Materiales cohesivos (ADR 0010 — sección paralela) ---
        known_cohesive_ids = set()
        registered_cohesives = set(CohesiveMaterialRegistry.names())
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

        # --- Materiales térmicos (Etapa 8 — sección paralela) ---
        known_thermal_ids = set()
        registered_thermals = set(ThermalMaterialRegistry.names())
        for i, mat in enumerate(data.get('thermal_materials', []) or []):
            ctx = f"thermal_materials[{i}]"
            if not isinstance(mat, dict):
                errors.append(f"{ctx}: cada material térmico debe ser un diccionario.")
                continue
            if 'id' not in mat:
                errors.append(f"{ctx}: falta el campo obligatorio 'id'.")
            else:
                known_thermal_ids.add(mat['id'])
            if 'type' not in mat:
                errors.append(f"{ctx} (id={mat.get('id', '?')}): falta el campo obligatorio 'type'.")
            elif registered_thermals and mat['type'] not in registered_thermals:
                errors.append(
                    f"{ctx} (id={mat.get('id', '?')}): tipo de material térmico desconocido "
                    f"'{mat['type']}'. Disponibles: {sorted(registered_thermals)}."
                )

        # --- Elementos (bloque inline, no mesh) ---
        registered_elements = set(ElementRegistry.names())
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
                elif (elem['material'] not in known_mat_ids
                      and elem['material'] not in known_thermal_ids):
                    # Un elemento térmico referencia un id de `thermal_materials`;
                    # uno mecánico, uno de `materials`. El campo es el mismo y la
                    # familia se resuelve por el bloque donde se declaró el id, de
                    # modo que el YAML no obliga al usuario a saber a qué registro
                    # pertenece cada material.
                    errors.append(
                        f"{ctx}: referencia a material inexistente (id={elem['material']}). "
                        f"Declarado ni en 'materials' ni en 'thermal_materials'."
                    )
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
                if e_type in registered_elements:
                    unknown, advertised = _unknown_kwargs(
                        ElementRegistry.get(e_type), set(elem),
                        {'id', 'type', 'material', 'nodes', 'cohesive_material', 'element_id'},
                    )
                    for kw in unknown:
                        errors.append(
                            f"{ctx}: parámetro '{kw}' no aceptado por '{e_type}'. "
                            f"Admitidos: {advertised}."
                        )

        # --- Malla externa + bloques inline: ambigüedad ---
        if has_mesh and (has_nodes or has_elements):
            errors.append(
                "'mesh' no puede combinarse con 'nodes'/'elements' inline: la malla "
                "externa sustituye por completo a la geometría del archivo."
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

        # --- Cargas puntuales (por nodo) ---
        for load in (data.get('point_loads', []) or []) + (data.get('point_loads_by_node', []) or []):
            if not isinstance(load, dict):
                errors.append("point_loads: cada entrada debe ser un diccionario.")
                continue
            nid = load.get('node_id')
            if nid is None:
                errors.append("point_loads: una entrada no tiene 'node_id'.")
            elif known_node_ids and nid not in known_node_ids:
                errors.append(f"point_loads: 'node_id={nid}' no existe en el bloque 'nodes'.")

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
                registered_solvers = set(SolverRegistry.names())
                if registered_solvers and s_type not in registered_solvers:
                    errors.append(
                        f"solver: tipo desconocido '{s_type}'. "
                        f"Disponibles: {sorted(registered_solvers)}."
                    )
                else:
                    # `convergence` se materializa aparte; `F_amplitude` y
                    # `F_func` los inyecta el parser cuando el solver los acepta.
                    unknown, advertised = _unknown_kwargs(
                        SolverRegistry.get(s_type), set(solver_cfg),
                        {'type', 'convergence', 'assembler'},
                    )
                    for kw in unknown:
                        errors.append(
                            f"solver: parámetro '{kw}' no aceptado por '{s_type}'. "
                            f"Admitidos: {advertised}."
                        )

            # Override del backend algebraico (ADR 0003 §4) — diagnóstico, no
            # decisión de modelado. Validamos contra el registro del despachador.
            la = solver_cfg.get('linear_algebra')
            if la is not None:
                from solidum.math.linalg.dispatcher import (
                    available_overrides, is_valid_override,
                )
                if not is_valid_override(la):
                    errors.append(
                        f"solver: 'linear_algebra' desconocido o no disponible: "
                        f"'{la}'. Disponibles: {available_overrides()}."
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

    def _parse_thermal_materials(self, data: dict):
        """Materiales térmicos (Etapa 8, sección paralela a ``materials``).

        Familia con registro propio, igual que los cohesivos (ADR 0010): el
        contrato es ``compute_flux(∇T)``, no ``compute_stress(ε)``, y la
        compatibilidad con el elemento la fija ``FLUX_DIM``, no ``STRAIN_DIM``.
        Separar los bloques evita que un material térmico y uno mecánico
        compartan espacio de nombres de ``type``.
        """
        for mat_data in data.get('thermal_materials', []) or []:
            mat_id = mat_data['id']
            mat_type = mat_data['type']
            kwargs = {k: v for k, v in mat_data.items() if k not in ('id', 'type')}
            self.thermal_materials[mat_id] = ThermalMaterialRegistry.create(
                mat_type, **kwargs,
            )

    def _parse_mesh_or_elements(self, data: dict):
        mesh_file = data.get('mesh', None)
        if mesh_file:
            from solidum.utils.gmsh_parser import GmshParser
            base_dir = os.path.dirname(self.filepath)
            mesh_path = os.path.join(base_dir, mesh_file)
            
            if not os.path.exists(mesh_path):
                mesh_path_sin_ext = os.path.join(base_dir, mesh_file.replace('.msh', ''))
                if os.path.exists(mesh_path_sin_ext):
                    mesh_path = mesh_path_sin_ext
            
            default_mat_id = data.get('mesh_material', 1)
            default_thickness = float(data.get('mesh_thickness', 1.0))
            default_quad_str = data.get('mesh_quadrature', '2x2')
            
            if default_mat_id not in self.materials:
                raise YamlValidationError([
                    f"mesh_material={default_mat_id!r} no existe en 'materials' "
                    f"(ids declarados: {sorted(self.materials)})."
                ])
            default_material = self.materials[default_mat_id]
            default_quadrature = self._get_quadrature(default_quad_str)
            
            physical_props = {}
            for group_name, props in data.get('mesh_physical_groups', {}).items():
                mat_id = props.get('material', default_mat_id)
                if mat_id not in self.materials:
                    raise YamlValidationError([
                        f"mesh_physical_groups[{group_name!r}]: material {mat_id!r} no "
                        f"existe en 'materials' (ids declarados: {sorted(self.materials)})."
                    ])
                mat = self.materials[mat_id]
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
                    # El id se busca en las dos familias: el bloque donde se
                    # declaró determina cuál. `materials` tiene prioridad para
                    # que un modelo puramente mecánico no cambie de comportamiento.
                    if mat_id in self.materials:
                        material = self.materials[mat_id]
                        es_termico = False
                    else:
                        material = self.thermal_materials[mat_id]
                        es_termico = True

                    kwargs = {k: v for k, v in elem_dict.items() if k not in ('id', 'type', 'material', 'nodes', 'cohesive_material')}
                    # `quadrature` viaja como CLAVE del registro: todos los
                    # elementos (2D, 3D y térmicos) la resuelven en su
                    # constructor vía ``resolve_quadrature`` y guardan
                    # ``quadrature_key`` para diagnósticos. Materializarla aquí
                    # como tupla rompía los sólidos 3D (auditoría 2026-09-22).
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

        for spec in self._normalize_coord_specs(data.get('boundary_conditions_by_coord', [])):
            for node in self._nodes_by_coord(spec, context='boundary_conditions_by_coord'):
                for dof, value in spec.items():
                    if dof not in self._COORD_KEYS:
                        node.fix_dof(dof, float(value))

        bcs_group = data.get('boundary_conditions_by_group', [])
        if isinstance(bcs_group, dict):
            bcs_group = [{'group_name': k, **v} for k, v in bcs_group.items()]
            
        for bcs in bcs_group or []:
            for node_id in self._group_node_ids(bcs.get('group_name'), 'boundary_conditions_by_group'):
                node = self.domain.get_node(node_id)
                for dof, value in bcs.items():
                    if dof != 'group_name':
                        node.fix_dof(dof, float(value))

    # ------------------------------------------------------------------
    # Selección de nodos por coordenada y por grupo físico
    # ------------------------------------------------------------------

    def _normalize_coord_specs(self, raw) -> list:
        """Forma lista de entradas ``*_by_coord``. La forma dict
        ``{x_min: {...}, ...}`` se convierte a lista inyectando ``loc``."""
        if isinstance(raw, dict):
            out = []
            for k, v in raw.items():
                v = dict(v)
                if k in self._COORD_LOCS:
                    v['loc'] = k
                out.append(v)
            return out
        return list(raw or [])

    def _nodes_by_coord(self, spec: dict, context: str) -> list:
        """Nodos (con DOFs) cuya coordenada casa con ``spec``: ``loc`` en
        ``x_min … z_max`` o ``coord``/``val``. Válido en 2D y 3D; antes
        desempaquetaba ``x, y = node.coordinates`` y reventaba con nodos 3D.
        Falla si ninguna coordenada casa (un selector que no selecciona nada
        es casi siempre un error de modelo)."""
        nodes = [n for n in self.domain.nodes.values() if n.dofs]
        axes = {'x': 0, 'y': 1, 'z': 2}

        def comp(n, ax):
            c = n.coordinates
            return float(c[ax]) if ax < len(c) else 0.0

        tol = float(spec.get('tol', 1e-6))
        loc, coord = spec.get('loc'), spec.get('coord')
        if loc is not None:
            if loc not in self._COORD_LOCS:
                raise ValueError(f"{context}: 'loc={loc}' no válido; use uno de {self._COORD_LOCS}.")
            ax = axes[loc[0]]
            vals = [comp(n, ax) for n in nodes]
            if not vals:
                return []
            target = min(vals) if loc.endswith('min') else max(vals)
        elif coord is not None and 'val' in spec:
            if coord not in axes:
                raise ValueError(f"{context}: 'coord={coord}' no válido; use x, y o z.")
            ax = axes[coord]
            target = float(spec['val'])
        else:
            raise ValueError(f"{context}: cada entrada necesita 'loc' o el par 'coord'/'val'.")
        matched = [n for n in nodes if abs(comp(n, ax) - target) < tol]
        if not matched:
            raise ValueError(f"{context}: ningún nodo casa con {spec}.")
        return matched

    def _group_node_ids(self, group_name, context: str) -> list:
        groups = getattr(self.domain, 'physical_groups', None)
        if not groups:
            raise ValueError(
                f"{context}: el modelo no tiene grupos físicos (sólo los aporta una "
                f"malla gmsh vía 'mesh')."
            )
        if group_name not in groups:
            raise ValueError(
                f"{context}: grupo físico {group_name!r} inexistente. "
                f"Disponibles: {sorted(groups)}."
            )
        return groups[group_name]

    @staticmethod
    def _add_nodal_load(F_ext, node, dof, value, context: str) -> None:
        if dof not in node.dofs:
            raise ValueError(
                f"{context}: el DOF '{dof}' no existe en el nodo {node.id} "
                f"(DOFs del nodo: {sorted(node.dofs)})."
            )
        F_ext[node.dofs[dof]] += float(value)

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

    def _parse_thermal_loads_config(self, data: dict):
        """Bloque ``thermal_loads`` (Etapa 8), con dos sub-bloques opcionales.

        Sección propia y no dentro de ``point_loads`` porque una fuente
        volumétrica y un flujo de frontera no son cargas nodales: se integran
        sobre el elemento y su borde/cara, exactamente como el peso propio no
        vive en ``point_loads``.
        """
        cfg = data.get('thermal_loads', {}) or {}
        self.thermal_body_sources = cfg.get('body_source', []) or []
        if isinstance(self.thermal_body_sources, dict):
            self.thermal_body_sources = [self.thermal_body_sources]
        self.thermal_boundary_fluxes = cfg.get('boundary_flux', []) or []
        if isinstance(self.thermal_boundary_fluxes, dict):
            self.thermal_boundary_fluxes = [self.thermal_boundary_fluxes]

    def _parse_point_loads_config(self, data: dict):
        self.point_loads = data.get('point_loads', [])
        self.point_loads_by_node = data.get('point_loads_by_node', [])

        # La forma dict conserva la clave como ``loc`` (antes se perdía y el
        # selector no casaba con ningún nodo en silencio).
        self.point_loads_by_coord = self._normalize_coord_specs(data.get('point_loads_by_coord', []))

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
        
        for load in (self.point_loads or []) + (self.point_loads_by_node or []):
            node_id = load.get('node_id')
            if node_id is None:
                raise ValueError("point_loads: una entrada no tiene 'node_id'.")
            node = self.domain.get_node(node_id)
            if node is None:
                raise ValueError(f"point_loads: el nodo {node_id} no existe en el modelo.")
            for dof, value in load.items():
                if dof != 'node_id':
                    self._add_nodal_load(F_ext, node, dof, value, 'point_loads')

        for spec in self.point_loads_by_coord:
            for node in self._nodes_by_coord(spec, context='point_loads_by_coord'):
                for dof, value in spec.items():
                    if dof not in self._COORD_KEYS:
                        self._add_nodal_load(F_ext, node, dof, value, 'point_loads_by_coord')

        for loads in self.point_loads_by_group or []:
            for node_id in self._group_node_ids(loads.get('group_name'), 'point_loads_by_group'):
                node = self.domain.get_node(node_id)
                for dof, value in loads.items():
                    if dof != 'group_name':
                        self._add_nodal_load(F_ext, node, dof, value, 'point_loads_by_group')

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

    def get_thermal_loads(self) -> np.ndarray:
        """Vector global de cargas térmicas ``F`` [W] declaradas en el YAML.

        Dos bloques, ambos opcionales:

        - ``body_source``: fuente volumétrica ``Q`` [W/m³] aplicada a un
          conjunto de elementos (o a todos). Se integra con
          ``element.compute_body_source(Q)``.
        - ``boundary_flux``: flujo prescrito ``q̄`` [W/m²] sobre bordes (2D) o
          caras (3D) de elementos concretos. **Convención de signo**
          (``Reglas.md §5``): ``q̄ > 0`` es flujo **saliente** del dominio
          (enfriamiento), de modo que el vector resultante es negativo.

        Un borde o cara sin declarar es **adiabático**, análogo al borde libre
        de tracción en mecánica: no hay que declarar el flujo nulo.

        Las cargas nodales concentradas van por el bloque estándar
        ``point_loads`` con ``T: <valor>``, que ya es genérico por nombre de
        DOF y no necesita tratamiento especial.
        """
        F = np.zeros(self.domain.total_dofs)

        for src in self.thermal_body_sources:
            Q = float(src.get('Q', src.get('value', 0.0)))
            elem_ids = src.get('elements')
            objetivo = (self.domain.elements.values() if elem_ids is None
                        else [self.domain.elements[e] for e in elem_ids])
            for elem in objetivo:
                f_e = elem.compute_body_source(Q)
                for a, node in enumerate(elem.nodes):
                    F[node.dofs["T"]] += f_e[a]

        for flux in self.thermal_boundary_fluxes:
            q_bar = float(flux.get('q', flux.get('value', 0.0)))
            elem = self.domain.elements[flux['element']]
            if 'edge' in flux:
                f_e = elem.compute_edge_flux(int(flux['edge']), q_bar)
            elif 'face' in flux:
                f_e = elem.compute_face_flux(int(flux['face']), q_bar)
            else:
                raise ValueError(
                    "boundary_flux: cada entrada necesita 'edge' (elementos 2D) "
                    "o 'face' (elementos 3D)."
                )
            for a, node in enumerate(elem.nodes):
                F[node.dofs["T"]] += f_e[a]

        return F

    def get_solver(self, assembler):
        """Construye y retorna el solver dinámicamente según la configuración YAML."""
        s_type = self.solver_config.get('type', 'LinearSolver')
        kwargs = {k: v for k, v in self.solver_config.items() if k != 'type'}

        solver_cls = SolverRegistry.get(s_type)
        accepted, _var_kw = _constructor_kwargs(solver_cls)

        # Bloque `convergence:` (ADR 0007). Se materializa como
        # ConvergenceCriterion si el constructor del solver lo acepta; si no
        # (solvers lineales, modal, armónico, espectral) se descarta CON AVISO.
        if 'convergence' in kwargs:
            from solidum.math.convergence import make_convergence_from_config
            cfg = kwargs.pop('convergence')
            if 'convergence' in accepted:
                kwargs['convergence'] = make_convergence_from_config(cfg)
            else:
                _log.warning(
                    f"solver: el bloque 'convergence' no aplica a {s_type} "
                    f"(sin Newton interno); se ignora."
                )

        # Solvers transitorios mecánicos: `F_func` es un callable que el YAML
        # no puede expresar. Las cargas del archivo (puntuales + peso propio o
        # fuerza de cuerpo) se aplican como ESCALÓN constante en el tiempo
        # desde t = 0. Antes se descartaban en silencio y el análisis corría
        # como vibración libre (auditoría 2026-09-22). Para una historia F(t)
        # arbitraria, construir el solver desde Python con su propio `F_func`.
        if (getattr(solver_cls, 'PIPELINE_KIND', None) == 'transient'
                and 'F_func' in accepted and 'F_func' not in kwargs):
            F_const = self.get_external_forces() + self.get_body_load(assembler)
            if np.any(F_const):
                _log.info(
                    f"solver {s_type}: las cargas del YAML se aplican como escalón "
                    f"constante en el tiempo (F_func = cte)."
                )
                kwargs['F_func'] = lambda t, _F=F_const: _F

        # HarmonicSolver: si no se pasó `F_amplitude` explícito, derivar la
        # amplitud compleja del bloque estándar de cargas externas del YAML
        # (point_loads / etc.). Las cargas se interpretan como amplitudes
        # reales con fase cero; para cargas con fase, el usuario pasa
        # `F_amplitude` numpy complejo directamente desde código (YAML no
        # soporta tipos complejos nativos).
        if s_type == 'HarmonicSolver' and 'F_amplitude' not in kwargs:
            kwargs['F_amplitude'] = self.get_external_forces()

        # Solver térmico transitorio: `F_func` es un callable y el YAML no
        # puede expresarlo. Se deriva el vector de carga térmica constante del
        # bloque `thermal_loads` (+ las nodales de `point_loads` con `T:`) y se
        # envuelve en una función del tiempo. Mismo patrón que `F_amplitude`
        # arriba. Para una carga variable en el tiempo, el usuario construye el
        # solver desde código y pasa su propio `F_func`.
        if (getattr(solver_cls, 'PIPELINE_KIND', None) == 'thermal_transient'
                and 'F_func' not in kwargs):
            F_termico = self.get_thermal_loads() + self.get_external_forces()
            if np.any(F_termico):
                kwargs['F_func'] = lambda t, _F=F_termico: _F

        return SolverRegistry.create(s_type, assembler=assembler, **kwargs)