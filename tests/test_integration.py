import unittest
import numpy as np
import sys
import os

# Asegurar que podemos importar fenix
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.domain import Domain
from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.solvers import NonlinearSolver, ArcLengthSolver
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.elements.truss import Truss2D

class TestSolversIntegration(unittest.TestCase):
    """
    Pruebas de integración que validan el flujo completo:
    Dominio -> Elementos -> Materiales -> Ensamblador -> Solucionador No Lineal.
    """
    def setUp(self):
        self.domain = Domain()
        self.n1 = self.domain.add_node(1, [0.0, 0.0])
        self.n2 = self.domain.add_node(2, [1.0, 0.0]) # Longitud L = 1.0

        # Material que fluye en 150 Pa y tiene endurecimiento de 50 Pa
        self.material = Elastoplastic1D(E=100.0, sigma_y=150.0, H=50.0, density=0.0)
        # El elemento registra los DoFs en los nodos al crearse
        self.element = Truss2D(1, [self.n1, self.n2], self.material, A=1.0)

        # Empotramos el nodo 1 y evitamos movimiento en Y del nodo 2
        self.n1.fix_dof('ux', 0.0)
        self.n1.fix_dof('uy', 0.0)
        self.n2.fix_dof('uy', 0.0)
        
        self.domain.add_element(self.element)
        self.domain.generate_equation_numbers(verbose=False)
        
        self.assembler = Assembler(self.domain)
        
    def test_nonlinear_solver_elastoplastic(self):
        F_ext = np.zeros(self.domain.total_dofs)
        F_ext[self.n2.dofs['ux']] = 200.0
        
        # 4 incrementos de carga (50 N por paso)
        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=4)
        U = solver.solve(F_ext)
        
        # Validación de cinemática macroscópica
        self.assertAlmostEqual(U[self.n2.dofs['ux']], 3.0, places=4, msg="El desplazamiento total no coincide con la solución analítica.")
        
        # Validación de variables microscópicas (historial del material)
        state = self.element.state.vars[0]
        self.assertIsNotNone(state, "El elemento no propagó sus variables de estado.")
        self.assertAlmostEqual(state['eps_p'], 1.0, places=4, msg="La deformación plástica no coincide.")
        self.assertAlmostEqual(state['alpha'], 1.0, places=4)

    def test_arclength_solver_elastoplastic(self):
        F_ext = np.zeros(self.domain.total_dofs)
        F_ext[self.n2.dofs['ux']] = 200.0
        
        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = ArcLengthSolver(self.assembler, convergence=conv, max_lambda=1.0, initial_dl=0.25)
        U = solver.solve(F_ext)
        
        self.assertAlmostEqual(U[self.n2.dofs['ux']], 3.0, places=4, msg="Arc-Length falló en encontrar el equilibrio final.")
        self.assertAlmostEqual(self.element.state.vars[0]['eps_p'], 1.0, places=4)

if __name__ == '__main__':
    unittest.main()