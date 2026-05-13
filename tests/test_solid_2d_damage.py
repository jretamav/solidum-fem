"""Tests de integración del pipeline Quad4 + IsotropicDamage2D + NonlinearSolver.

Acceptance de ``docs/specs/IsotropicDamage2D.md`` sección ``acceptance.integration``.

Hasta ahora la cobertura del modelo de daño 2D era exclusivamente unitaria
(`test_materials_unit.py · TestIsotropicDamage2D`). Estos tests validan que la
tangente algorítmica consistente recupera convergencia rápida del Newton global
en presencia de daño activo, y que el cableado solver-material-elemento no
introduce regresiones cuando el material declara `IS_SYMMETRIC=False`.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.domain import Domain
from fenix.elements.solid_2d import Quad4
from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.solvers import NonlinearSolver


def _build_uniaxial_quad4_domain(material):
    """Cuadrado unitario con 1 Quad4 y patrón de apoyos para tracción uniaxial.

    - Nodo (0,0): ux=0, uy=0
    - Nodo (1,0): uy=0
    - Nodo (0,1): ux=0
    """
    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [1.0, 0.0])
    n3 = domain.add_node(3, [1.0, 1.0])
    n4 = domain.add_node(4, [0.0, 1.0])

    element = Quad4(1, [n1, n2, n3, n4], material, thickness=1.0)
    domain.add_element(element)

    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
    n2.fix_dof('uy', 0.0)
    n4.fix_dof('ux', 0.0)

    domain.generate_equation_numbers(verbose=False)
    return domain, element, {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4}


class TestQuad4Damage2DPipeline(unittest.TestCase):
    """Quad4 + IsotropicDamage2D + NonlinearSolver — pipeline completo."""

    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.kappa_0 = 1.0e-4
        self.alpha = 100.0
        self.material = IsotropicDamage2D(
            E=self.E, nu=self.nu, kappa_0=self.kappa_0, alpha=self.alpha,
            hypothesis='plane_stress',
        )
        self.domain, self.element, self.nodes = _build_uniaxial_quad4_domain(self.material)
        self.assembler = Assembler(self.domain)

    def _apply_traction(self, F_total):
        F = np.zeros(self.domain.total_dofs)
        F[self.nodes['n2'].dofs['ux']] = 0.5 * F_total
        F[self.nodes['n3'].dofs['ux']] = 0.5 * F_total
        return F

    def test_elastic_regime_below_threshold(self):
        """Carga por debajo del umbral: sin daño, respuesta elástica pura."""
        # ε_xx target = 0.5·κ_0 ⇒ σ_xx target = E·0.5·κ_0 (plane stress uniaxial)
        eps_target = 0.5 * self.kappa_0
        F_total = self.E * eps_target  # σ_xx · área(=1)
        F_ext = self._apply_traction(F_total)

        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=2)
        U = solver.solve(F_ext)

        u_x_2 = U[self.nodes['n2'].dofs['ux']]
        self.assertAlmostEqual(u_x_2, eps_target, places=10)
        self.assertEqual(self.element.state.vars[0]['damage'], 0.0)

    def test_damage_active_converges_in_few_iter(self):
        """Carga sobre el umbral: daño activo + tangente consistente → ≤6 iter/paso.

        Con tangente secante, este mismo problema requería ~10-20 iteraciones
        por paso (convergencia lineal). Con tangente consistente algorítmica
        el Newton recupera convergencia cuadrática.
        """
        # ε_xx target = 4·κ_0 → daño moderado (no saturado)
        eps_target = 4.0 * self.kappa_0
        # Estimación grosera de σ a esa deformación: (1-d)·E·ε con d~0.3-0.5
        F_total = 0.6 * self.E * eps_target
        F_ext = self._apply_traction(F_total)

        conv = ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7)
        # max_iter=6: el test falla si el Newton necesita más
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=4, max_iter=6)
        U = solver.solve(F_ext)

        self.assertGreater(U[self.nodes['n2'].dofs['ux']], self.kappa_0)
        self.assertGreater(self.element.state.vars[0]['damage'], 0.0)

    def test_unloading_no_damage_increment(self):
        """Tras carga dañadora, descarga: κ y d no cambian."""
        # Paso 1: cargar dañando
        eps_load = 5.0 * self.kappa_0
        F_load = 0.5 * self.E * eps_load
        F_ext_load = self._apply_traction(F_load)
        conv = ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7)
        solver_load = NonlinearSolver(self.assembler, convergence=conv, num_steps=3, max_iter=8)
        solver_load.solve(F_ext_load)

        state_after_load = dict(self.element.state.vars[0])
        kappa_loaded = state_after_load['kappa']
        d_loaded = state_after_load['damage']
        self.assertGreater(d_loaded, 0.0)

        # Paso 2: descargar a un nivel inferior (no llega al umbral pero κ memoriza)
        F_unload = 0.3 * self.E * eps_load
        F_ext_unload = self._apply_traction(F_unload)
        # Nuevo solver (reinicia load_factor de 0 a 1, pero el state ya tiene memoria)
        solver_unload = NonlinearSolver(self.assembler, convergence=conv, num_steps=3, max_iter=8)
        solver_unload.solve(F_ext_unload)

        state_after_unload = self.element.state.vars[0]
        # Irreversibilidad: κ y d no cambiaron durante la descarga
        self.assertAlmostEqual(state_after_unload['kappa'], kappa_loaded, places=12)
        self.assertAlmostEqual(state_after_unload['damage'], d_loaded, places=12)

    def test_damage_progresses_under_increasing_load(self):
        """En 4 pasos de carga creciente, el daño debe avanzar paso a paso."""
        eps_target = 6.0 * self.kappa_0
        F_total = 0.5 * self.E * eps_target
        F_ext = self._apply_traction(F_total)

        damages = []

        def step_callback(step, U, load_factor):
            damages.append(self.element.state.vars[0]['damage'])

        conv = ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7)
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=4, max_iter=8)
        solver.solve(F_ext, step_callback=step_callback)

        # Daños registrados en cada paso convergido: monótonos no decrecientes
        for i in range(1, len(damages)):
            self.assertGreaterEqual(damages[i], damages[i - 1] - 1e-14)
        # Daño final > 0
        self.assertGreater(damages[-1], 0.0)


if __name__ == '__main__':
    unittest.main()
