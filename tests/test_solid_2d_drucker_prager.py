"""Tests de integración del pipeline Quad4 + DruckerPrager2D + NonlinearSolver.

Acceptance de ``docs/specs/DruckerPrager2D.md`` sección ``acceptance.integration``.

Valida que el modelo Drucker-Prager con plasticidad no asociada se integra
correctamente con el pipeline del solver. Como la tangente es asimétrica
(IS_SYMMETRIC=False), el despachador algebraico usa LU automáticamente.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.domain import Domain
from fenix.elements.solid_2d import Quad4
from fenix.materials.drucker_prager_2d import DruckerPrager2D
from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.solvers import NonlinearSolver


def _build_quad4_plane_strain_domain(material, pattern: str = 'confined_y'):
    """Cuadrado unitario con 1 Quad4. Patrones de apoyo:
    - 'confined_y': uy=0 en todos los nodos; ux=0 en lado izquierdo.
    - 'corner_fixed': nodo (0,0) totalmente fijo, (1,0) en uy, (0,1) en ux.
    """
    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [1.0, 0.0])
    n3 = domain.add_node(3, [1.0, 1.0])
    n4 = domain.add_node(4, [0.0, 1.0])

    element = Quad4(1, [n1, n2, n3, n4], material, thickness=1.0)
    domain.add_element(element)

    if pattern == 'confined_y':
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
        n2.fix_dof('uy', 0.0)
        n3.fix_dof('uy', 0.0)
        n4.fix_dof('ux', 0.0); n4.fix_dof('uy', 0.0)
    elif pattern == 'corner_fixed':
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
        n2.fix_dof('uy', 0.0)
        n4.fix_dof('ux', 0.0)
    else:
        raise ValueError(f"pattern desconocido: {pattern}")

    domain.generate_equation_numbers(verbose=False)
    return domain, element, {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4}


class TestQuad4DruckerPrager2DPipeline(unittest.TestCase):
    """Quad4 + DruckerPrager2D + NonlinearSolver — pipeline completo."""

    def setUp(self):
        self.E = 2.0e4
        self.nu = 0.3
        self.c0 = 10.0
        self.phi = 30.0
        self.psi = 10.0   # no asociada
        self.H = 100.0
        self.material = DruckerPrager2D(
            E=self.E, nu=self.nu, cohesion=self.c0,
            phi_deg=self.phi, psi_deg=self.psi, H=self.H,
            hypothesis='plane_strain',
        )

    def test_elastic_regime_below_yield(self):
        """Carga muy por debajo del yield: respuesta elástica pura."""
        domain, element, nodes = _build_quad4_plane_strain_domain(
            self.material, pattern='confined_y'
        )
        assembler = Assembler(domain)

        # Cortante mediante deformación angular pequeña: aplicar tracción
        # tangencial baja en el lado derecho. En realidad, con apoyos
        # 'confined_y' completos, la única deformación posible es ε_xx, así
        # que pasamos a tracción horizontal pequeña.
        F = np.zeros(domain.total_dofs)
        F_total = 1.0   # muy por debajo del yield
        F[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=2)
        U = solver.solve(F)

        self.assertEqual(element.state.vars[0]['alpha'], 0.0)
        self.assertGreater(U[nodes['n2'].dofs['ux']], 0.0)

    def test_plastic_regime_apex_branch_converges(self):
        """Carga uniaxial confinada empuja al ápice; solver converge."""
        domain, element, nodes = _build_quad4_plane_strain_domain(
            self.material, pattern='confined_y'
        )
        assembler = Assembler(domain)

        # Carga suficientemente alta para plastificar
        F = np.zeros(domain.total_dofs)
        F_total = 30.0
        F[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7)
        # max_iter=10: holgura porque la rama de ápice puede ser numéricamente
        # más sutil que la regular
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=4, max_iter=10)
        U = solver.solve(F)

        self.assertGreater(U[nodes['n2'].dofs['ux']], 0.0)
        self.assertGreater(element.state.vars[0]['alpha'], 0.0)

    def test_associated_flow_pipeline_converges(self):
        """Caso asociado (ψ = φ): tangente simétrica, el despachador puede usar
        Cholesky en lugar de LU. Verifica que el solver converge igualmente.
        """
        mat_assoc = DruckerPrager2D(
            E=self.E, nu=self.nu, cohesion=self.c0,
            phi_deg=self.phi, psi_deg=self.phi, H=self.H,
        )
        domain, element, nodes = _build_quad4_plane_strain_domain(
            mat_assoc, pattern='confined_y'
        )
        assembler = Assembler(domain)

        F = np.zeros(domain.total_dofs)
        F_total = 30.0
        F[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=4, max_iter=10)
        U = solver.solve(F)

        self.assertGreater(element.state.vars[0]['alpha'], 0.0)


if __name__ == '__main__':
    unittest.main()
