"""Cobertura del pipeline: sólidos 2D no-Quad4 + materiales no lineales 2D.

Cierra la deuda técnica #2 documentada en ``docs/STATUS.md``. Hasta ahora,
el catálogo declaraba compatibilidad de ``Tri3``, ``Quad8``, ``Quad9``,
``Tri6`` con ``VonMises2D``, ``DruckerPrager2D`` e ``IsotropicDamage2D``,
pero sólo ``Quad4`` estaba blindado por test de sistema
(``test_solid_2d_plasticity.py``, ``test_solid_2d_drucker_prager.py``,
``test_solid_2d_damage.py``).

Estos tests no replican los asserts físicos exhaustivos: verifican el
**cableado** — la combinación elemento + material no lineal +
``NonlinearSolver`` corre el pipeline completo, converge y produce
evidencia de actividad inelástica (``alpha > 0`` para plasticidad,
``damage > 0`` para daño).
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad8, Quad9, Tri3, Tri6
from solidum.materials.damage_2d import IsotropicDamage2D
from solidum.materials.drucker_prager_2d import DruckerPrager2D
from solidum.materials.von_mises_2d import VonMises2D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver


def _build_tri3_mesh(material):
    """Cuadrado unitario subdividido en 2 Tri3 a través de la diagonal n1-n3."""
    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [1.0, 0.0])
    n3 = domain.add_node(3, [1.0, 1.0])
    n4 = domain.add_node(4, [0.0, 1.0])
    elem_a = Tri3(1, [n1, n2, n3], material, thickness=1.0)
    elem_b = Tri3(2, [n1, n3, n4], material, thickness=1.0)
    domain.add_element(elem_a)
    domain.add_element(elem_b)
    return domain, [elem_a, elem_b], {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4}


def _build_quad8_mesh(material):
    """Cuadrado unitario con 1 Quad8 (vértices + medios de borde)."""
    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [1.0, 0.0])
    n3 = domain.add_node(3, [1.0, 1.0])
    n4 = domain.add_node(4, [0.0, 1.0])
    n5 = domain.add_node(5, [0.5, 0.0])   # medio 0-1 (lado inferior)
    n6 = domain.add_node(6, [1.0, 0.5])   # medio 1-2 (lado derecho)
    n7 = domain.add_node(7, [0.5, 1.0])   # medio 2-3 (lado superior)
    n8 = domain.add_node(8, [0.0, 0.5])   # medio 3-0 (lado izquierdo)
    elem = Quad8(1, [n1, n2, n3, n4, n5, n6, n7, n8], material, thickness=1.0)
    domain.add_element(elem)
    return domain, [elem], {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4,
                            'n5': n5, 'n6': n6, 'n7': n7, 'n8': n8}


def _build_quad9_mesh(material):
    """Cuadrado unitario con 1 Quad9 (vértices + medios + central)."""
    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [1.0, 0.0])
    n3 = domain.add_node(3, [1.0, 1.0])
    n4 = domain.add_node(4, [0.0, 1.0])
    n5 = domain.add_node(5, [0.5, 0.0])
    n6 = domain.add_node(6, [1.0, 0.5])
    n7 = domain.add_node(7, [0.5, 1.0])
    n8 = domain.add_node(8, [0.0, 0.5])
    n9 = domain.add_node(9, [0.5, 0.5])   # nodo central
    elem = Quad9(1, [n1, n2, n3, n4, n5, n6, n7, n8, n9],
                 material, thickness=1.0)
    domain.add_element(elem)
    return domain, [elem], {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4,
                            'n5': n5, 'n6': n6, 'n7': n7, 'n8': n8, 'n9': n9}


def _build_tri6_mesh(material):
    """Cuadrado unitario con 2 Tri6 compartiendo la diagonal n1-n3.

    Convención Tri6: ``EDGE_NODES = ((0, 3, 1), (1, 4, 2), (2, 5, 0))`` ⇒
    conectividad ``[v0, v1, v2, m_01, m_12, m_20]``.

    Triángulo A: vértices (n1, n2, n3); medios (n5=0.5,0), (n6=1,0.5),
    (n7=0.5,0.5).
    Triángulo B: vértices (n1, n3, n4); medios (n7 compartido), (n8=0.5,1),
    (n9=0,0.5).
    """
    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [1.0, 0.0])
    n3 = domain.add_node(3, [1.0, 1.0])
    n4 = domain.add_node(4, [0.0, 1.0])
    n5 = domain.add_node(5, [0.5, 0.0])
    n6 = domain.add_node(6, [1.0, 0.5])
    n7 = domain.add_node(7, [0.5, 0.5])    # medio de la diagonal (compartido)
    n8 = domain.add_node(8, [0.5, 1.0])
    n9 = domain.add_node(9, [0.0, 0.5])
    elem_a = Tri6(1, [n1, n2, n3, n5, n6, n7], material, thickness=1.0)
    elem_b = Tri6(2, [n1, n3, n4, n7, n8, n9], material, thickness=1.0)
    domain.add_element(elem_a)
    domain.add_element(elem_b)
    return domain, [elem_a, elem_b], {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4,
                                       'n5': n5, 'n6': n6, 'n7': n7,
                                       'n8': n8, 'n9': n9}


# Registro de builders por tipo de elemento.
MESH_BUILDERS = {
    'Tri3':  _build_tri3_mesh,
    'Quad8': _build_quad8_mesh,
    'Quad9': _build_quad9_mesh,
    'Tri6':  _build_tri6_mesh,
}


def _apply_uniaxial_dirichlet(nodes):
    """Patrón Dirichlet uniaxial robusto, válido para los 4 tipos de malla:

    - ``ux = 0`` en todos los nodos del lado izquierdo (``x = 0``).
    - ``uy = 0`` en ``(0, 0)`` y ``(1, 0)`` para fijar la rotación rígida.

    Permite contracción transversal libre por Poisson. Compatible con
    plane_strain (``ε_yy`` no impuesta a cero), plane_stress y ambos
    indistintamente.
    """
    TOL = 1e-9
    for node in nodes.values():
        x, y = node.coordinates[0], node.coordinates[1]
        if abs(x) < TOL:
            node.fix_dof('ux', 0.0)
        if abs(y) < TOL and (abs(x) < TOL or abs(x - 1.0) < TOL):
            node.fix_dof('uy', 0.0)


def _apply_right_edge_traction(domain, nodes, sigma_xx):
    """Fuerza horizontal sobre el lado derecho (``x = 1``), reparto consistente
    según número de nodos en ese lado:

    - 2 nodos (lineal): 1/2 + 1/2.
    - 3 nodos (cuadrático): 1/6 + 4/6 + 1/6.

    Asume área = 1 (lado de longitud 1, espesor 1) ⇒ ``F_total = sigma_xx``.
    """
    TOL = 1e-9
    right = [n for n in nodes.values() if abs(n.coordinates[0] - 1.0) < TOL]
    right.sort(key=lambda n: n.coordinates[1])

    F = np.zeros(domain.total_dofs)
    if len(right) == 2:
        F[right[0].dofs['ux']] = 0.5 * sigma_xx
        F[right[1].dofs['ux']] = 0.5 * sigma_xx
    elif len(right) == 3:
        F[right[0].dofs['ux']] = (1.0 / 6.0) * sigma_xx
        F[right[1].dofs['ux']] = (4.0 / 6.0) * sigma_xx
        F[right[2].dofs['ux']] = (1.0 / 6.0) * sigma_xx
    else:
        raise AssertionError(
            f"Lado derecho con {len(right)} nodos no soportado por este helper."
        )
    return F


def _max_state_value(elements, key):
    """Máximo del valor ``key`` (e.g. 'alpha', 'damage') sobre todos los
    puntos de Gauss de todos los elementos. Permite tolerar mallas con
    múltiples elementos donde la plasticidad/daño puede activarse sólo
    en uno por la asimetría de la triangulación."""
    best = 0.0
    for elem in elements:
        for gp_state in elem.state.vars:
            best = max(best, gp_state.get(key, 0.0))
    return best


# ---------------------------------------------------------------------------
# VonMises2D — plane_strain
# ---------------------------------------------------------------------------

class TestVonMises2DPlaneStrainCoverage(unittest.TestCase):
    """Pipeline Tri3/Quad8/Quad9/Tri6 + VonMises2D plane_strain converge y
    plastifica."""

    def _run(self, kind):
        mat = VonMises2D(E=2.0e5, nu=0.3, sigma_y=250.0, H=1.0e4,
                         hypothesis='plane_strain')
        domain, elements, nodes = MESH_BUILDERS[kind](mat)
        _apply_uniaxial_dirichlet(nodes)
        domain.generate_equation_numbers(verbose=False)
        # Plane strain uniaxial: σ_xx_yield ≈ E·ε_y, con ε_y dependiente del
        # estado plano. 600 N/m·m supera el yield para los parámetros dados.
        F_ext = _apply_right_edge_traction(domain, nodes, sigma_xx=600.0)

        assembler = Assembler(domain)
        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = NonlinearSolver(assembler, convergence=conv,
                                  num_steps=4, max_iter=15)
        U = solver.solve(F_ext)

        # Desplazamiento horizontal positivo en algún nodo del lado derecho
        right_dof_ux = nodes['n2'].dofs['ux']
        self.assertGreater(U[right_dof_ux], 0.0,
            f"{kind}: lado derecho no se desplazó hacia +x")
        # Plastificó al menos un punto de Gauss
        alpha = _max_state_value(elements, 'alpha')
        self.assertGreater(alpha, 0.0,
            f"{kind}: ningún punto de Gauss plastificó (alpha máx = 0)")

    def test_tri3(self):  self._run('Tri3')
    def test_quad8(self): self._run('Quad8')
    def test_quad9(self): self._run('Quad9')
    def test_tri6(self):  self._run('Tri6')


# ---------------------------------------------------------------------------
# VonMises2D — plane_stress
# ---------------------------------------------------------------------------

class TestVonMises2DPlaneStressCoverage(unittest.TestCase):
    """Pipeline Tri3/Quad8/Quad9/Tri6 + VonMises2D plane_stress converge y
    plastifica."""

    def _run(self, kind):
        mat = VonMises2D(E=2.0e5, nu=0.3, sigma_y=250.0, H=1.0e4,
                         hypothesis='plane_stress')
        domain, elements, nodes = MESH_BUILDERS[kind](mat)
        _apply_uniaxial_dirichlet(nodes)
        domain.generate_equation_numbers(verbose=False)
        # Tracción uniaxial plane stress: σ_yield = sigma_y · A = 250.
        # 320 N/m·m supera el yield uniaxial.
        F_ext = _apply_right_edge_traction(domain, nodes, sigma_xx=320.0)

        assembler = Assembler(domain)
        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = NonlinearSolver(assembler, convergence=conv,
                                  num_steps=4, max_iter=15)
        U = solver.solve(F_ext)

        self.assertGreater(U[nodes['n2'].dofs['ux']], 0.0,
            f"{kind}: lado derecho no se desplazó hacia +x")
        alpha = _max_state_value(elements, 'alpha')
        self.assertGreater(alpha, 0.0,
            f"{kind}: ningún punto de Gauss plastificó (alpha máx = 0)")

    def test_tri3(self):  self._run('Tri3')
    def test_quad8(self): self._run('Quad8')
    def test_quad9(self): self._run('Quad9')
    def test_tri6(self):  self._run('Tri6')


# ---------------------------------------------------------------------------
# DruckerPrager2D — plane_strain (única hipótesis soportada)
# ---------------------------------------------------------------------------

class TestDruckerPrager2DCoverage(unittest.TestCase):
    """Pipeline Tri3/Quad8/Quad9/Tri6 + DruckerPrager2D + Nonlinear."""

    def _run(self, kind):
        mat = DruckerPrager2D(E=2.0e4, nu=0.3, cohesion=10.0,
                              phi_deg=30.0, psi_deg=10.0, H=100.0,
                              hypothesis='plane_strain')
        domain, elements, nodes = MESH_BUILDERS[kind](mat)
        _apply_uniaxial_dirichlet(nodes)
        domain.generate_equation_numbers(verbose=False)
        # Carga moderada que plastifica con los parámetros del test Quad4
        F_ext = _apply_right_edge_traction(domain, nodes, sigma_xx=30.0)

        assembler = Assembler(domain)
        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = NonlinearSolver(assembler, convergence=conv,
                                  num_steps=4, max_iter=15)
        U = solver.solve(F_ext)

        self.assertGreater(U[nodes['n2'].dofs['ux']], 0.0,
            f"{kind}: lado derecho no se desplazó hacia +x")
        alpha = _max_state_value(elements, 'alpha')
        self.assertGreater(alpha, 0.0,
            f"{kind}: ningún punto de Gauss plastificó (alpha máx = 0)")

    def test_tri3(self):  self._run('Tri3')
    def test_quad8(self): self._run('Quad8')
    def test_quad9(self): self._run('Quad9')
    def test_tri6(self):  self._run('Tri6')


# ---------------------------------------------------------------------------
# IsotropicDamage2D — plane_stress (por consistencia con test_solid_2d_damage)
# ---------------------------------------------------------------------------

class TestIsotropicDamage2DCoverage(unittest.TestCase):
    """Pipeline Tri3/Quad8/Quad9/Tri6 + IsotropicDamage2D + Nonlinear."""

    def _run(self, kind):
        kappa_0 = 1.0e-4
        E = 2.0e5
        mat = IsotropicDamage2D(E=E, nu=0.3, kappa_0=kappa_0, alpha=100.0,
                                hypothesis='plane_stress')
        domain, elements, nodes = MESH_BUILDERS[kind](mat)
        _apply_uniaxial_dirichlet(nodes)
        domain.generate_equation_numbers(verbose=False)
        # Carga que lleva a ε_xx ~ 4·κ_0 (daño moderado, no saturado)
        # estimación grosera: σ_xx ~ 0.6·E·ε_xx por el factor (1-d)
        sigma_target = 0.6 * E * 4.0 * kappa_0
        F_ext = _apply_right_edge_traction(domain, nodes,
                                             sigma_xx=sigma_target)

        assembler = Assembler(domain)
        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = NonlinearSolver(assembler, convergence=conv,
                                  num_steps=4, max_iter=15)
        U = solver.solve(F_ext)

        self.assertGreater(U[nodes['n2'].dofs['ux']], 0.0,
            f"{kind}: lado derecho no se desplazó hacia +x")
        damage = _max_state_value(elements, 'damage')
        self.assertGreater(damage, 0.0,
            f"{kind}: ningún punto de Gauss dañado (damage máx = 0)")

    def test_tri3(self):  self._run('Tri3')
    def test_quad8(self): self._run('Quad8')
    def test_quad9(self): self._run('Quad9')
    def test_tri6(self):  self._run('Tri6')


if __name__ == '__main__':
    unittest.main()
