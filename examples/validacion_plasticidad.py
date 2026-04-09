# fenix_fem/examples/validacion_plasticidad.py
import sys, os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from fenix.core.domain import Domain
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.elements.structural import Truss2D
from fenix.math.solvers import NonlinearSolver
from fenix.math.assembly import Assembler

def main():
    print("--- FENIX FEM: ENSAYO DE PLASTICIDAD 1D ---")
    mesh = Domain()
    
    # Barra simple de 1 metro
    n1 = mesh.add_node(1, [0.0, 0.0])
    n2 = mesh.add_node(2, [1.0, 0.0])
    
    # Acero elastoplástico: E = 200 GPa, Fluencia = 250 MPa, Endurecimiento = 10 GPa
    acero_plastico = Elastoplastic1D(E=200e9, sigma_y=250e6, H=10e9)
    A = 0.01  # Área transversal 0.01 m2
    
    e1 = Truss2D(1, [n1, n2], acero_plastico, A)
    mesh.add_element(e1)
    
    # Empotramos nodo 1
    n1.fix_dof('ux', 0.0)
    n1.fix_dof('uy', 0.0)
    n2.fix_dof('uy', 0.0)
    mesh.generate_equation_numbers()
    
    # La deformación de fluencia teórica es eps_y = 250e6 / 200e9 = 0.00125
    # Imponemos un desplazamiento que la supere por mucho (eps = 0.002)
    desplazamiento_impuesto = 0.002
    n2.fix_dof('ux', desplazamiento_impuesto)
    
    F_ext = np.zeros(mesh.total_dofs)
    
    assembler = Assembler(mesh)
    solver = NonlinearSolver(assembler, tol=1e-5, max_iter=10)
    
    print("\nAplicando desplazamiento en régimen inelástico...")
    U = solver.solve(F_ext)
    
    resultados = e1.compute_internal_forces(U)
    esfuerzo = resultados['stress'] / 1e6
    def_total = resultados['strain']
    def_plastica = e1.state_vars[0]['eps_p']
    
    print("\n[RESULTADOS FINALES]")
    print(f"Deformación Total Impuesta: {def_total:.5f}")
    print(f"Deformación Plástica Acumulada: {def_plastica:.5f}")
    print(f"Esfuerzo Final: {esfuerzo:.2f} MPa")
    
    # Comprobación teórica
    esfuerzo_teorico = 250.0 + (10e9 / 1e6) * def_plastica
    print(f"Esfuerzo Teórico Esperado: {esfuerzo_teorico:.2f} MPa")

if __name__ == "__main__":
    main()
