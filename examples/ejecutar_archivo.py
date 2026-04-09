# fenix_fem/examples/ejecutar_archivo.py
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.utils.parser import FEAParser
from fenix.math.solvers import LinearSolver

def main():
    archivo_entrada = os.path.join(os.path.dirname(__file__), 'modelo_placa.txt')
    print(f"--- FENIX FEM: LEYENDO ARCHIVO DE ENTRADA ---")
    print(f"Archivo: {archivo_entrada}\n")
    
    # 1. Parsear el archivo (Lee datos y memoriza fuerzas)
    parser = FEAParser(archivo_entrada)
    mesh, assembler = parser.parse()
    
    # 2. Ensamblar el sistema (Esto crea K_global y F_global llenos de ceros)
    assembler.assemble_system()
    
    # 3. Aplicar Cargas y Condiciones de Frontera
    parser.apply_parsed_forces() # <--- AHORA SÍ vertemos las fuerzas memorizadas
    assembler.apply_dirichlet_bcs()
    
    # 4. Resolver
    solver = LinearSolver(assembler.K_global, assembler.F_global)
    U = solver.solve()
    
    # 5. Resultados
    print("[DESPLAZAMIENTOS]")
    for node_id, node in mesh.nodes.items():
        print(f"Nodo {node_id}: ux = {U[node.dofs['ux']]:12.6e} m, uy = {U[node.dofs['uy']]:12.6e} m")

if __name__ == "__main__":
    main()
