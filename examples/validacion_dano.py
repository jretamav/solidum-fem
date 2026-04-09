# fenix_fem/examples/validacion_dano.py
import sys, os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fenix.core.domain import Domain
from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.elements.solid_2d import Quad4
from fenix.math.solvers import NonlinearSolver
from fenix.math.assembly import Assembler

def main():
    print("--- FENIX FEM: ENSAYO NO LINEAL CON MODELO DE DAÑO ---")
    mesh = Domain()
    
    # Malla de 1 elemento cuadrado (1.0 x 1.0 m)
    n1 = mesh.add_node(1, [0.0, 0.0])
    n2 = mesh.add_node(2, [1.0, 0.0])
    n3 = mesh.add_node(3, [1.0, 1.0])
    n4 = mesh.add_node(4, [0.0, 1.0])
    
    # Material con degradación: Umbral de deformación elástica kappa_0 = 0.001 (0.1%)
    # Alpha = 1000 controla la velocidad de degradación post-pico
    mat_dano = IsotropicDamage2D(E=20e9, nu=0.2, kappa_0=0.001, alpha=1000.0)
    
    mesh.add_element(Quad4(1, [n1, n2, n3, n4], mat_dano, thickness=0.1))
    
    # Empotramiento izquierdo
    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
    n4.fix_dof('ux', 0.0); n4.fix_dof('uy', 0.0)
    mesh.generate_equation_numbers()
    
    # Aplicar un desplazamiento impuesto destructivo en la cara derecha (ux = 0.003 m)
    n2.fix_dof('ux', 0.003)
    n3.fix_dof('ux', 0.003)
    n2.fix_dof('uy', 0.0)
    n3.fix_dof('uy', 0.0)
    
    F_ext = np.zeros(mesh.total_dofs)
    
    assembler = Assembler(mesh)
    solver = NonlinearSolver(assembler, tol=1e-5, max_iter=15)
    U = solver.solve(F_ext)
    
    print("\n[ESTADO INTERNO DEL MATERIAL]")
    elem = mesh.elements[1]
    for i, state in enumerate(elem.state_vars):
        d = state['damage']
        print(f"Punto de Gauss {i+1}: Nivel de Daño d = {d*100:.2f}%")

if __name__ == "__main__":
    main()
