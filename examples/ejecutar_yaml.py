# fenix_fem/examples/ejecutar_yaml.py
import sys, os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fenix.utils.yaml_parser import YamlParser
from fenix.math.assembly import Assembler
from fenix.utils.vtk_exporter import VtkExporter

def main():
    print("--- FENIX FEM: MOTOR DE EJECUCIÓN YAML ---")
    
    # Exigir el archivo de entrada YAML como argumento de la terminal
    if len(sys.argv) < 2:
        print("[!] Error: Debes especificar el archivo YAML de entrada.")
        print("Uso correcto: python ejecutar_yaml.py <ruta_al_archivo.yaml>")
        return
        
    ruta_yaml = sys.argv[1]

    if not os.path.exists(ruta_yaml):
        print(f"[!] Error: No se encontró el archivo {ruta_yaml}")
        return
        
    # 1. Leer el archivo y construir el dominio
    parser = YamlParser(ruta_yaml)
    mesh = parser.parse()
    
    # 2. Generar ecuaciones y ensamblar
    mesh.generate_equation_numbers()
    assembler = Assembler(mesh)
    
    # 3. Construir el solver dinámicamente según el YAML
    solver = parser.get_solver(assembler)
    
    F_ext = parser.get_external_forces()
    
    out_cfg = parser.output_config
    pre_cfg = out_cfg.get('pre_process', {})
    res_cfg = out_cfg.get('results', {})
    
    base_dir = os.path.dirname(ruta_yaml)
    nombre_base = out_cfg.get('file_name', os.path.basename(ruta_yaml).replace('.yaml', ''))
    ruta_base = os.path.join(base_dir, nombre_base)
    freq = res_cfg.get('frequency', 'last_step')
    n_vars = res_cfg.get('nodal_results', None)
    e_vars = res_cfg.get('element_results', None)
    
    # Configuración de Sensores (Exportación a Texto/CSV)
    text_cfg = out_cfg.get('text_export', {})
    export_text = text_cfg.get('export', False)
    
    _nodes_in = text_cfg.get('nodes', [])
    text_nodes = list(mesh.nodes.keys()) if isinstance(_nodes_in, str) and _nodes_in.lower() == 'all' else _nodes_in
    
    _elems_in = text_cfg.get('elements', [])
    text_elems = list(mesh.elements.keys()) if isinstance(_elems_in, str) and _elems_in.lower() == 'all' else _elems_in
    
    text_n_vars = text_cfg.get('nodal_results', [])
    text_e_vars = text_cfg.get('element_results', [])
    
    text_data = []
    
    # 4. Exportar Pre-proceso
    if pre_cfg.get('export', False):
        ruta_pre = ruta_base + '_pre_proceso.vtu'
        VtkExporter(mesh).export(ruta_pre, U=None, F_ext=F_ext, nodal_vars=n_vars, elem_vars=e_vars)
        print(f"  -> Modelo de pre-proceso exportado a: {ruta_pre}")

    # 5. Resolver con capacidad de Animación
    pvd_data = []
    def exportar_paso(step, U_step, factor):
        if freq == 'all_steps':
            vtu_filename = f"{nombre_base}_{step:03d}.vtu"
            step_file = os.path.join(base_dir, vtu_filename)
            VtkExporter(mesh).export(step_file, U_step, F_ext * factor, n_vars, e_vars)
            pvd_data.append((factor, vtu_filename))
            
        if export_text:
            step_info = {'step': step, 'factor': factor, 'nodes': {}, 'elems': {}}
            for nid in text_nodes:
                node = mesh.get_node(nid)
                if node:
                    step_info['nodes'][nid] = {}
                    if not text_n_vars or "Displacements" in text_n_vars:
                        ux = U_step[node.dofs['ux']] if 'ux' in node.dofs else 0.0
                        uy = U_step[node.dofs['uy']] if 'uy' in node.dofs else 0.0
                        step_info['nodes'][nid]['U'] = (ux, uy)
                    if not text_n_vars or "External_Forces" in text_n_vars:
                        fx = (F_ext[node.dofs['ux']] * factor) if 'ux' in node.dofs else 0.0
                        fy = (F_ext[node.dofs['uy']] * factor) if 'uy' in node.dofs else 0.0
                        step_info['nodes'][nid]['F'] = (fx, fy)
            for eid in text_elems:
                elem = mesh.elements.get(eid)
                if elem:
                    step_info['elems'][eid] = {}
                    if not text_e_vars or "Von_Mises" in text_e_vars:
                        vm = 0.0
                        if hasattr(elem, 'stresses'):
                            s_avg = np.mean(elem.stresses, axis=0)
                            sx, sy, txy = s_avg[0], s_avg[1], s_avg[2]
                            vm = np.sqrt(sx**2 + sy**2 - sx*sy + 3.0*txy**2)
                        step_info['elems'][eid]['Von_Mises'] = vm
                    if not text_e_vars or "Internal_State" in text_e_vars:
                        state_val = 0.0
                        if hasattr(elem, 'state_vars') and elem.state_vars[0] is not None:
                            vals = [sv.get('alpha', sv.get('damage', 0.0)) for sv in elem.state_vars if sv is not None]
                            if vals: state_val = sum(vals) / len(vals)
                        step_info['elems'][eid]['Internal_State'] = state_val
            text_data.append(step_info)
            
    if hasattr(solver, 'num_steps') or hasattr(solver, 'max_steps'):
        U = solver.solve(F_ext, step_callback=exportar_paso)
    else:
        U = solver.solve(F_ext)
        exportar_paso(1, U, 1.0)
        
    # Generar archivo maestro PVD para ParaView
    if pvd_data:
        ruta_pvd = ruta_base + '.pvd'
        with open(ruta_pvd, 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n')
            for factor, vtu_name in pvd_data:
                f.write(f'  <DataSet timestep="{factor:.6f}" group="" part="0" file="{vtu_name}"/>\n')
            f.write('</Collection>\n</VTKFile>\n')
        print(f"  -> Archivo maestro de animación (PVD) exportado a: {ruta_pvd}")
        
    # Generar archivo de texto (TXT) Estilo FEAP
    if export_text and text_data:
        ruta_txt = ruta_base + '_text_export.txt'
        with open(ruta_txt, 'w') as f:
            # ENCABEZADO TIPO FEAP
            f.write(" * * * * *  F E N I X   F E M   O U T P U T  * * * * *\n\n")
            f.write(f"     Problem Name            - - - - - - : {nombre_base}\n")
            f.write(f"     Number of Nodal Points  - - - - - - : {len(mesh.nodes)}\n")
            f.write(f"     Number of Elements  - - - - - - - - : {len(mesh.elements)}\n")
            f.write(f"     Degrees-of-Freedom/Node (Maximum) - : 2\n")
            f.write(f"     Total Degrees-of-Freedom  - - - - - : {mesh.total_dofs}\n")
            
            for step_info in text_data:
                f.write("\n ---------------------------------------------------\n\n")
                f.write(f"     Step:{step_info['step']:8d}    Load Factor:{step_info['factor']:14.5E}\n")
                
                if not text_n_vars or "Displacements" in text_n_vars:
                    f.write("\n     N o d a l   D i s p l a c e m e n t s\n\n")
                    f.write("      Node       1 Displ       2 Displ\n")
                    for nid, data in step_info['nodes'].items():
                        if 'U' in data:
                            f.write(f"{nid:10d}{data['U'][0]:14.4E}{data['U'][1]:14.4E}\n")
                            
                if not text_n_vars or "External_Forces" in text_n_vars:
                    f.write("\n     N o d a l   F o r c e s\n\n")
                    f.write("      Node       1 Force       2 Force\n")
                    for nid, data in step_info['nodes'].items():
                        if 'F' in data:
                            f.write(f"{nid:10d}{data['F'][0]:14.4E}{data['F'][1]:14.4E}\n")
                            
                if step_info['elems']:
                    f.write("\n     E l e m e n t   V a r i a b l e s\n\n")
                    header = "      Elmt"
                    has_vm = not text_e_vars or "Von_Mises" in text_e_vars
                    has_is = not text_e_vars or "Internal_State" in text_e_vars
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
    
    # 6. Resumen Final
    ruta_salida = ruta_base + '_final.vtu'
    VtkExporter(mesh).export(ruta_salida, U, F_ext, n_vars, e_vars)
    print(f"\n¡Análisis completado! Resultados exportados en: {base_dir}")

if __name__ == "__main__":
    main()