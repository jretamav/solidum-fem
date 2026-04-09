# fenix_fem/examples/ejecutar_armadura.py
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.utils.parser import FEAParser
from fenix.math.solvers import LinearSolver

def main():
    archivo_entrada = os.path.join(os.path.dirname(__file__), 'modelo_armadura_feap.txt')
    print("--- FENIX FEM: ANÁLISIS DE ARMADURA ---")
    print(f"Leyendo archivo: {os.path.basename(archivo_entrada)}\n")
    
    # 1. Parsear el archivo
    parser = FEAParser(archivo_entrada)
    mesh, assembler = parser.parse()
    
    # 2. Ensamblar sistema
    assembler.assemble_system()
    
    # 3. Aplicar condiciones y fuerzas memorizadas
    parser.apply_parsed_forces() 
    assembler.apply_dirichlet_bcs()
    
    # 4. Resolver
    solver = LinearSolver(assembler.K_global, assembler.F_global)
    U = solver.solve()
    
    # 5. Imprimir Desplazamientos
    print("[DESPLAZAMIENTOS]")
    for node_id in sorted(mesh.nodes.keys()):
        node = mesh.nodes[node_id]
        ux = U[node.dofs['ux']]
        uy = U[node.dofs['uy']]
        print(f"Nodo {node_id}: ux = {ux:12.6e} m, uy = {uy:12.6e} m")

if __name__ == "__main__":
    main()
