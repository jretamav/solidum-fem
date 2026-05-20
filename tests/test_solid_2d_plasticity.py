"""Tests de integración del pipeline Quad4 + VonMises2D + NonlinearSolver.

Cierra una brecha de cobertura: hasta ahora ``test_integration`` solo cubría
plasticidad 1D (``Truss2D + Elastoplastic1D``). Estos tests validan el
flujo completo con sólidos 2D y plasticidad J2 en sus dos hipótesis
cinemáticas (plane strain, plane stress).

Acceptance de ``docs/specs/VonMises2D.md`` sección ``acceptance.integration``.
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.materials.von_mises_2d import VonMises2D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver


def _build_unit_square_domain(material, dirichlet_pattern: str):
    """Cuadrado unitario con 1 Quad4 y patrón de apoyos según ``dirichlet_pattern``.

    Devuelve ``(domain, element, nodes_dict)``.

    Patrones soportados:
    - 'confined_y': uy=0 en todos los nodos; ux=0 en lado izquierdo. Fuerza
      ``ε_yy = 0`` (plane strain transverso) y deja ε_xx libre.
    - 'free_y':     ux=0 en lado izquierdo, uy=0 solo en nodo (0,0) y (1,0).
      Permite contracción transversal por Poisson — tracción uniaxial pura.
    """
    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [1.0, 0.0])
    n3 = domain.add_node(3, [1.0, 1.0])
    n4 = domain.add_node(4, [0.0, 1.0])

    element = Quad4(1, [n1, n2, n3, n4], material, thickness=1.0)
    domain.add_element(element)

    if dirichlet_pattern == 'confined_y':
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
        n2.fix_dof('uy', 0.0)
        n3.fix_dof('uy', 0.0)
        n4.fix_dof('ux', 0.0); n4.fix_dof('uy', 0.0)
    elif dirichlet_pattern == 'free_y':
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
        n2.fix_dof('uy', 0.0)
        n4.fix_dof('ux', 0.0)
    else:
        raise ValueError(f"Patrón Dirichlet desconocido: {dirichlet_pattern}")

    domain.generate_equation_numbers(verbose=False)
    return domain, element, {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4}


class TestQuad4VonMises2DPlaneStrainPipeline(unittest.TestCase):
    """Quad4 + VonMises2D plane strain + NonlinearSolver.

    Setup: cuadrado unitario con ``ε_yy = γ_xy = 0`` impuesto por Dirichlet
    (plane strain confinado). Carga horizontal en el borde derecho. La
    respuesta debe coincidir con el material ejecutado standalone bajo la
    misma trayectoria monotónica de ε_xx.
    """

    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.sigma_y = 250.0
        self.H = 1.0e4
        self.material = VonMises2D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
            hypothesis='plane_strain',
        )
        self.domain, self.element, self.nodes = _build_unit_square_domain(
            self.material, dirichlet_pattern='confined_y'
        )
        self.assembler = Assembler(self.domain)

    def _apply_horizontal_traction(self, F_total):
        """Tracción uniforme en lado derecho: F_total repartida entre nodos 2 y 3."""
        F_ext = np.zeros(self.domain.total_dofs)
        F_ext[self.nodes['n2'].dofs['ux']] = 0.5 * F_total
        F_ext[self.nodes['n3'].dofs['ux']] = 0.5 * F_total
        return F_ext

    def test_elastic_regime(self):
        """Carga por debajo del primer yield: respuesta puramente elástica."""
        G = self.E / (2.0 * (1.0 + self.nu))
        eps_xx_yield = self.sigma_y / (2.0 * G)  # primer yield J2 plane strain confinado
        # Cargamos hasta 80% del yield
        eps_xx_target = 0.8 * eps_xx_yield
        C00 = self.E * (1.0 - self.nu) / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        sigma_xx_target = C00 * eps_xx_target
        F_total = sigma_xx_target  # área = 1·1

        F_ext = self._apply_horizontal_traction(F_total)
        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=2)
        U = solver.solve(F_ext)

        u_x_2 = U[self.nodes['n2'].dofs['ux']]
        self.assertAlmostEqual(u_x_2, eps_xx_target, places=8)
        # α debe ser exactamente 0 (régimen elástico)
        state = self.element.state.vars[0]
        self.assertAlmostEqual(state['alpha'], 0.0, places=12)

    def test_plastic_regime_matches_standalone_material(self):
        """En régimen plástico, σ_xx(ε_xx) FEM coincide con material directo."""
        G = self.E / (2.0 * (1.0 + self.nu))
        C00 = self.E * (1.0 - self.nu) / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        eps_xx_yield = self.sigma_y / (2.0 * G)
        sigma_xx_yield_elastic = C00 * eps_xx_yield
        # Carga 50% por encima del yield elástico para entrar bien en plástico
        F_total = 1.5 * sigma_xx_yield_elastic
        F_ext = self._apply_horizontal_traction(F_total)

        conv = ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7)
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=10)
        U = solver.solve(F_ext)

        u_x_final = U[self.nodes['n2'].dofs['ux']]
        self.assertGreater(u_x_final, eps_xx_yield)  # entró en plástico

        # σ_xx FEM (promedio sobre puntos de Gauss del único elemento)
        gs = self.element.compute_gauss_state(U)
        sigma_xx_fem = gs['stress'][:, 0].mean()
        alpha_fem = self.element.state.vars[0]['alpha']

        # Material standalone bajo la misma trayectoria monotónica
        # (J2 + endurecimiento lineal monotónico es path-independent en ε).
        mat_ref = VonMises2D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
            hypothesis='plane_strain',
        )
        eps_path = np.linspace(0.0, u_x_final, 50)
        state_ref = None
        sigma_ref_last = None
        for eps_xx in eps_path:
            eps = np.array([eps_xx, 0.0, 0.0])
            sigma_ref_last, _, state_ref = mat_ref.compute_state(eps, state_vars=state_ref)
        alpha_ref = state_ref['alpha']

        # Tolerancia 1e-3 abs porque ε es del orden 1e-3 y σ del orden 1e2
        self.assertAlmostEqual(sigma_xx_fem, sigma_ref_last[0], delta=abs(sigma_ref_last[0]) * 1e-3)
        self.assertAlmostEqual(alpha_fem, alpha_ref, delta=abs(alpha_ref) * 1e-3 + 1e-10)
        self.assertGreater(alpha_fem, 0.0)


class TestQuad4VonMises2DPlaneStressPipeline(unittest.TestCase):
    """Quad4 + VonMises2D plane stress + NonlinearSolver.

    Setup: cuadrado unitario con tracción uniaxial pura (lado izquierdo
    apoyado en ux, esquinas inf. apoyadas en uy). σ_yy = σ_xy = 0 emergente
    del equilibrio; ε_yy libre por Poisson. La curva σ_xx-ε_xx debe coincidir
    con ``Elastoplastic1D``.
    """

    def setUp(self):
        self.E = 2.0e5
        self.nu = 0.3
        self.sigma_y = 250.0
        self.H = 1.0e4
        self.material = VonMises2D(
            E=self.E, nu=self.nu, sigma_y=self.sigma_y, H=self.H,
            hypothesis='plane_stress',
        )
        self.domain, self.element, self.nodes = _build_unit_square_domain(
            self.material, dirichlet_pattern='free_y'
        )
        self.assembler = Assembler(self.domain)

    def _apply_horizontal_traction(self, F_total):
        F_ext = np.zeros(self.domain.total_dofs)
        F_ext[self.nodes['n2'].dofs['ux']] = 0.5 * F_total
        F_ext[self.nodes['n3'].dofs['ux']] = 0.5 * F_total
        return F_ext

    def test_elastic_uniaxial_response(self):
        """σ_xx = E · ε_xx (no E/(1-ν²)) porque σ_yy=0 emergente."""
        eps_xx_target = 0.8 * self.sigma_y / self.E
        F_total = self.E * eps_xx_target  # σ_xx = F_total/área (área=1)
        F_ext = self._apply_horizontal_traction(F_total)
        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=2)
        U = solver.solve(F_ext)

        u_x_2 = U[self.nodes['n2'].dofs['ux']]
        self.assertAlmostEqual(u_x_2, eps_xx_target, places=6)
        u_y_3 = U[self.nodes['n3'].dofs['uy']]
        # Contracción transversal por Poisson (ε_yy = -ν · ε_xx)
        self.assertAlmostEqual(u_y_3, -self.nu * eps_xx_target, places=6)
        self.assertAlmostEqual(self.element.state.vars[0]['alpha'], 0.0, places=12)

    def test_plastic_uniaxial_matches_1D(self):
        """Curva σ_xx-ε_xx FEM coincide con Elastoplastic1D bajo la misma trayectoria."""
        eps_target = 3.0 * self.sigma_y / self.E  # 3× el yield uniaxial
        F_total = 1.2 * self.sigma_y  # tracción suficientemente plástica
        F_ext = self._apply_horizontal_traction(F_total)

        conv = ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7)
        solver = NonlinearSolver(self.assembler, convergence=conv, num_steps=10)
        U = solver.solve(F_ext)

        u_x_final = U[self.nodes['n2'].dofs['ux']]
        gs = self.element.compute_gauss_state(U)
        sigma_xx_fem = gs['stress'][:, 0].mean()
        sigma_yy_fem = gs['stress'][:, 1].mean()
        alpha_fem = self.element.state.vars[0]['alpha']

        # σ_yy debe ser ≈ 0 (tracción uniaxial pura emergente)
        self.assertAlmostEqual(sigma_yy_fem, 0.0, delta=abs(sigma_xx_fem) * 1e-3)

        # Comparar contra Elastoplastic1D bajo la misma trayectoria ε_xx
        mat_1d = Elastoplastic1D(E=self.E, sigma_y=self.sigma_y, H=self.H)
        eps_path = np.linspace(0.0, u_x_final, 50)
        state_1d = None
        sigma_1d_last = 0.0
        for eps_xx in eps_path:
            sigma_1d_last, _, state_1d = mat_1d.compute_state(eps_xx, state_vars=state_1d)
        alpha_1d = state_1d['alpha']

        # σ_xx ↔ σ 1D y α ↔ α 1D
        self.assertAlmostEqual(sigma_xx_fem, sigma_1d_last, delta=abs(sigma_1d_last) * 1e-3)
        self.assertAlmostEqual(alpha_fem, alpha_1d, delta=abs(alpha_1d) * 5e-3 + 1e-10)
        self.assertGreater(alpha_fem, 0.0)


class TestQuad4VonMises2DConvergence(unittest.TestCase):
    """Convergencia del Newton global con tangente algorítmica consistente.

    En un paso plástico activo, el Newton debe converger en pocas iteraciones
    (≤ 6) con criterio de tolerancia razonable, evidencia indirecta de que
    la tangente consistente está bien derivada (de lo contrario la
    convergencia sería lineal y requeriría muchas más iteraciones).
    """

    def test_plane_strain_converges_in_few_iter(self):
        material = VonMises2D(E=2.0e5, nu=0.3, sigma_y=250.0, H=1.0e4,
                              hypothesis='plane_strain')
        domain, element, nodes = _build_unit_square_domain(material, 'confined_y')
        assembler = Assembler(domain)

        # Carga plástica significativa
        F_total = 600.0
        F_ext = np.zeros(domain.total_dofs)
        F_ext[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F_ext[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        # max_iter pequeño: el test falla si el Newton necesita muchas iteraciones
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=4, max_iter=8)
        U = solver.solve(F_ext)

        # Llegamos al final del análisis (load_factor 1.0)
        self.assertGreater(U[nodes['n2'].dofs['ux']], 0.0)
        self.assertGreater(element.state.vars[0]['alpha'], 0.0)

    def test_plane_stress_converges_in_few_iter(self):
        material = VonMises2D(E=2.0e5, nu=0.3, sigma_y=250.0, H=1.0e4,
                              hypothesis='plane_stress')
        domain, element, nodes = _build_unit_square_domain(material, 'free_y')
        assembler = Assembler(domain)

        F_total = 380.0  # > sigma_y = 250
        F_ext = np.zeros(domain.total_dofs)
        F_ext[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F_ext[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=4, max_iter=8)
        U = solver.solve(F_ext)

        self.assertGreater(U[nodes['n2'].dofs['ux']], 0.0)
        self.assertGreater(element.state.vars[0]['alpha'], 0.0)


if __name__ == '__main__':
    unittest.main()
