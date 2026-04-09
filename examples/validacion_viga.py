# fenix_fem/examples/validacion_viga.py
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fenix.core.domain import Domain
from fenix.elements.structural import Truss2D
from fenix.materials.elastic import Elastic1D
from fenix.math.assembly import Assembler
from fenix.math.solvers import LinearSolver

def main():
    print("--- FENIX FEM: ARMADURA 2D ---")
    mesh = Domain()
    n1 = mesh.add_node(1, [0.0, 0.0])
    n2 = mesh.add_node(2, [3.0, 0.0])
    n3 = mesh.add_node(3, [1.5, 2.0])
    
    acero = Elastic1D(E=200e9)
    A = 0.01
    
    mesh.add_element(Truss2D(1, [n1, n3], acero, A))
    mesh.add_element(Truss2D(2, [n2, n3], acero, A))
    
    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
    n2.fix_dof('ux', 0.0); n2.fix_dof('uy', 0.0)
    
    assembler = Assembler(mesh)
    assembler.assemble_system()
    assembler.apply_point_load(3, 'uy', -100000.0)
    assembler.apply_dirichlet_bcs()
    
    solver = LinearSolver(assembler.K_global, assembler.F_global)
    U = solver.solve()
    
    print("\n[DESPLAZAMIENTOS]")
    for node_id, node in mesh.nodes.items():
        print(f"Nodo {node_id}: ux = {U[node.dofs['ux']]:.6e} m, uy = {U[node.dofs['uy']]:.6e} m")

if __name__ == "__main__":
    main()
