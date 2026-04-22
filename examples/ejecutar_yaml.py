# fenix_fem/examples/ejecutar_yaml.py
import sys, os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import fenix
from fenix.utils.yaml_parser import YamlParser
from fenix.math.assembly import Assembler
from fenix.utils.vtk_exporter import VtkExporter


# ---------------------------------------------------------------------------
# Argumentos y validación de entrada
# ---------------------------------------------------------------------------

def _parse_args() -> str:
    """Extrae y valida la ruta del archivo YAML desde los argumentos de terminal."""
    if len(sys.argv) < 2:
        print("[!] Error: Debes especificar el archivo YAML de entrada.")
        print("Uso correcto: python ejecutar_yaml.py <ruta_al_archivo.yaml>")
        sys.exit(1)
    ruta = sys.argv[1]
    if not os.path.exists(ruta):
        print(f"[!] Error: No se encontró el archivo '{ruta}'")
        sys.exit(1)
    return ruta


# ---------------------------------------------------------------------------
# Construcción del modelo
# ---------------------------------------------------------------------------

def _build_model(ruta_yaml: str):
    """Parsea el YAML y construye parser, assembler, solver y vector de fuerzas.

    Returns
    -------
    parser, assembler, solver, F_ext
    """
    parser = YamlParser(ruta_yaml)
    mesh = parser.parse()
    mesh.generate_equation_numbers(verbose=True)

    assembler = Assembler(mesh)
    solver = parser.get_solver(assembler)
    F_ext = parser.get_external_forces()

    return parser, assembler, solver, F_ext


# ---------------------------------------------------------------------------
# Exportación VTK pre-proceso
# ---------------------------------------------------------------------------

def _export_pre_process(mesh, F_ext, ruta_base: str, n_vars, e_vars):
    """Exporta el modelo sin deformar (pre-proceso) a VTK."""
    ruta_pre = ruta_base + '_pre_proceso.vtu'
    VtkExporter(mesh).export(ruta_pre, U=None, F_ext=F_ext, nodal_vars=n_vars, elem_vars=e_vars)
    print(f"  -> Modelo de pre-proceso exportado a: {ruta_pre}")


# ---------------------------------------------------------------------------
# Callback de paso — animación y recolección de datos de texto
# ---------------------------------------------------------------------------

def _make_step_callback(mesh, F_ext, ruta_base: str, base_dir: str,
                        freq: str, n_vars, e_vars,
                        export_text: bool, text_nodes: list, text_elems: list,
                        text_n_vars: list, text_e_vars: list,
                        pvd_data: list, text_data: list, nombre_base: str):
    """Construye y retorna el callback que se invoca al final de cada paso convergido."""

    def _collect_node_data(U_step, factor) -> dict:
        result = {}
        for nid in text_nodes:
            node = mesh.get_node(nid)
            if not node:
                continue
            entry = {}
            if not text_n_vars or "Displacements" in text_n_vars:
                ux = U_step[node.dofs['ux']] if 'ux' in node.dofs else 0.0
                uy = U_step[node.dofs['uy']] if 'uy' in node.dofs else 0.0
                entry['U'] = (ux, uy)
            if not text_n_vars or "External_Forces" in text_n_vars:
                fx = (F_ext[node.dofs['ux']] * factor) if 'ux' in node.dofs else 0.0
                fy = (F_ext[node.dofs['uy']] * factor) if 'uy' in node.dofs else 0.0
                entry['F'] = (fx, fy)
            result[nid] = entry
        return result

    def _collect_elem_data() -> dict:
        result = {}
        for eid in text_elems:
            elem = mesh.elements.get(eid)
            if not elem:
                continue
            entry = {}
            if not text_e_vars or "Von_Mises" in text_e_vars:
                vm = 0.0
                if hasattr(elem, 'stresses'):
                    s_avg = np.mean(elem.stresses, axis=0)
                    sx, sy, txy = s_avg[0], s_avg[1], s_avg[2]
                    vm = np.sqrt(sx**2 + sy**2 - sx*sy + 3.0*txy**2)
                entry['Von_Mises'] = vm
            if not text_e_vars or "Internal_State" in text_e_vars:
                state_val = 0.0
                state_vars_list = None
                if hasattr(elem, 'state') and hasattr(elem.state, 'vars'):
                    state_vars_list = elem.state.vars
                elif hasattr(elem, 'state_vars'):
                    state_vars_list = elem.state_vars
                if state_vars_list and state_vars_list[0] is not None:
                    primary_key = getattr(getattr(elem, 'material', None), 'PRIMARY_STATE_VAR', None)
                    vals = [sv[primary_key] if (primary_key and primary_key in sv and
                            isinstance(sv[primary_key], (int, float)))
                            else next((v for v in sv.values()
                                       if isinstance(v, (int, float)) and not isinstance(v, bool)), 0.0)
                            for sv in state_vars_list if sv is not None]
                    if vals:
                        state_val = sum(vals) / len(vals)
                entry['Internal_State'] = state_val
            result[eid] = entry
        return result

    def exportar_paso(step, U_step, factor):
        if freq == 'all_steps':
            vtu_filename = f"{nombre_base}_{step:03d}.vtu"
            step_file = os.path.join(base_dir, vtu_filename)
            VtkExporter(mesh).export(step_file, U_step, F_ext * factor, n_vars, e_vars)
            pvd_data.append((factor, vtu_filename))

        if export_text:
            text_data.append({
                'step': step,
                'factor': factor,
                'nodes': _collect_node_data(U_step, factor),
                'elems': _collect_elem_data(),
            })

    return exportar_paso


# ---------------------------------------------------------------------------
# Ejecución del solver
# ---------------------------------------------------------------------------

def _run_solver(domain, assembler, solver, F_ext, step_callback):
    """Ejecuta el pipeline oficial ``fenix.run`` y retorna ``U``.

    Delega en ``fenix.run`` para centralizar Assembler→solver→SolveResult
    (ADR 0002); ``domain.last_result`` queda poblado tras el retorno. Para
    solvers lineales (sin ``step_callback`` nativo) sigue invocando el callback
    una sola vez post-solve, preservando la semántica de exportación del script.
    """
    has_native_callback = hasattr(solver, 'num_steps') or hasattr(solver, 'max_steps')
    result = fenix.run(
        domain,
        assembler=assembler,
        solver=solver,
        F_applied=F_ext,
        step_callback=step_callback if has_native_callback else None,
    )
    if not has_native_callback:
        step_callback(1, result.U, 1.0)
    return result.U


# ---------------------------------------------------------------------------
# Exportaciones post-proceso
# ---------------------------------------------------------------------------

def _export_pvd(pvd_data: list, ruta_base: str):
    """Genera el archivo maestro PVD para animaciones en ParaView."""
    if not pvd_data:
        return
    ruta_pvd = ruta_base + '.pvd'
    with open(ruta_pvd, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n')
        for factor, vtu_name in pvd_data:
            f.write(f'  <DataSet timestep="{factor:.6f}" group="" part="0" file="{vtu_name}"/>\n')
        f.write('</Collection>\n</VTKFile>\n')
    print(f"  -> Archivo maestro de animación (PVD) exportado a: {ruta_pvd}")


def _export_text(text_data: list, mesh, ruta_base: str, nombre_base: str,
                 text_n_vars: list, text_e_vars: list):
    """Genera el archivo de resultados en texto plano estilo FEAP."""
    if not text_data:
        return
    ruta_txt = ruta_base + '_text_export.txt'
    with open(ruta_txt, 'w') as f:
        f.write(" * * * * *  F E N I X   F E M   O U T P U T  * * * * *\n\n")
        f.write(f"     Problem Name            - - - - - - : {nombre_base}\n")
        f.write(f"     Number of Nodal Points  - - - - - - : {len(mesh.nodes)}\n")
        f.write(f"     Number of Elements  - - - - - - - - : {len(mesh.elements)}\n")
        f.write(f"     Degrees-of-Freedom/Node (Maximum) - : 2\n")
        f.write(f"     Total Degrees-of-Freedom  - - - - - : {mesh.total_dofs}\n")

        for step_info in text_data:
            f.write("\n ---------------------------------------------------\n\n")
            f.write(f"     Step:{step_info['step']:8d}    Load Factor:{step_info['factor']:14.5E}\n")

            has_disp = not text_n_vars or "Displacements" in text_n_vars
            has_force = not text_n_vars or "External_Forces" in text_n_vars
            has_vm = not text_e_vars or "Von_Mises" in text_e_vars
            has_is = not text_e_vars or "Internal_State" in text_e_vars

            if has_disp:
                f.write("\n     N o d a l   D i s p l a c e m e n t s\n\n")
                f.write("      Node       1 Displ       2 Displ\n")
                for nid, data in step_info['nodes'].items():
                    if 'U' in data:
                        f.write(f"{nid:10d}{data['U'][0]:14.4E}{data['U'][1]:14.4E}\n")

            if has_force:
                f.write("\n     N o d a l   F o r c e s\n\n")
                f.write("      Node       1 Force       2 Force\n")
                for nid, data in step_info['nodes'].items():
                    if 'F' in data:
                        f.write(f"{nid:10d}{data['F'][0]:14.4E}{data['F'][1]:14.4E}\n")

            if step_info['elems']:
                f.write("\n     E l e m e n t   V a r i a b l e s\n\n")
                header = "      Elmt"
                if has_vm: header += "     Von_Mises"
                if has_is: header += "     Int_State"
                f.write(header + "\n")
                for eid, data in step_info['elems'].items():
                    line = f"{eid:10d}"
                    if has_vm and 'Von_Mises' in data:
                        line += f"{data['Von_Mises']:14.4E}"
                    if has_is and 'Internal_State' in data:
                        line += f"{data['Internal_State']:14.4E}"
                    f.write(line + "\n")

    print(f"  -> Archivo de texto (Estilo FEAP) exportado a: {ruta_txt}")


def _export_final(mesh, U, F_ext, ruta_base: str, n_vars, e_vars):
    """Exporta el resultado final del análisis a VTK."""
    ruta_salida = ruta_base + '_final.vtu'
    VtkExporter(mesh).export(ruta_salida, U, F_ext, n_vars, e_vars)
    print(f"\n¡Análisis completado! Resultados exportados en: {os.path.dirname(ruta_base)}")


# ---------------------------------------------------------------------------
# Punto de entrada
# ---------------------------------------------------------------------------

def main():
    print("--- FENIX FEM: MOTOR DE EJECUCIÓN YAML ---")

    ruta_yaml = _parse_args()

    # 1. Construir modelo
    parser, assembler, solver, F_ext = _build_model(ruta_yaml)
    mesh = assembler.domain

    # 2. Leer configuración de salida
    out_cfg = parser.output_config
    pre_cfg = out_cfg.get('pre_process', {})
    res_cfg = out_cfg.get('results', {})
    text_cfg = out_cfg.get('text_export', {})

    base_dir = os.path.dirname(ruta_yaml)
    nombre_base = out_cfg.get('file_name', os.path.basename(ruta_yaml).replace('.yaml', ''))
    ruta_base = os.path.join(base_dir, nombre_base)

    freq = res_cfg.get('frequency', 'last_step')
    n_vars = res_cfg.get('nodal_results', None)
    e_vars = res_cfg.get('element_results', None)

    export_text = text_cfg.get('export', False)
    _nodes_in = text_cfg.get('nodes', [])
    text_nodes = list(mesh.nodes.keys()) if isinstance(_nodes_in, str) and _nodes_in.lower() == 'all' else _nodes_in
    _elems_in = text_cfg.get('elements', [])
    text_elems = list(mesh.elements.keys()) if isinstance(_elems_in, str) and _elems_in.lower() == 'all' else _elems_in
    text_n_vars = text_cfg.get('nodal_results', [])
    text_e_vars = text_cfg.get('element_results', [])

    # 3. Pre-proceso
    if pre_cfg.get('export', False):
        _export_pre_process(mesh, F_ext, ruta_base, n_vars, e_vars)

    # 4. Resolver
    pvd_data = []
    text_data = []
    step_callback = _make_step_callback(
        mesh, F_ext, ruta_base, base_dir, freq, n_vars, e_vars,
        export_text, text_nodes, text_elems, text_n_vars, text_e_vars,
        pvd_data, text_data, nombre_base,
    )
    U = _run_solver(mesh, assembler, solver, F_ext, step_callback)

    # 5. Post-proceso
    _export_pvd(pvd_data, ruta_base)
    _export_text(text_data, mesh, ruta_base, nombre_base, text_n_vars, text_e_vars)
    _export_final(mesh, U, F_ext, ruta_base, n_vars, e_vars)


if __name__ == "__main__":
    main()
