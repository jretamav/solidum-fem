"""Tests de la imposición de Dirichlet por eliminación directa (ADR 0004 fase 1).

Cubre los tres criterios de aceptación añadidos al ADR:

1. **Asentamiento prescrito**: viga simplemente apoyada con descenso impuesto
   en un apoyo. La solución analítica de Bernoulli-Euler para un cuerpo rígido
   girando da rotación constante y desplazamiento lineal — sin esfuerzos
   internos, sin reacciones espurias.

2. **Simetría de K_red**: tras la reducción, ``K_red − K_redᵀ`` ≈ 0 a
   redondeo cuando ``K`` es simétrica (lo es para los elementos del catálogo).

3. **Equilibrio global de reacciones**: las reacciones reportadas por
   ``SolveResult`` suman exactamente la fuerza externa aplicada (Newton 3ª).
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.bc.constraints import ConstraintSet
from solidum.core.domain import Domain
from solidum.elements.frame import Frame2DEuler
from solidum.elements.truss import Truss2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver
from solidum.results import build_solve_result


class LinearElastic1D:
    STRAIN_DIM = 1

    def __init__(self, E, nu=0.3, density: float = 0.0):
        self.E = E
        self.nu = nu

    def compute_state(self, strain, state_var=None):
        return self.E * strain, self.E, state_var


class TestConstraintSet(unittest.TestCase):
    """Validaciones de la abstracción ``ConstraintSet`` independiente del Assembler."""

    def test_dirichlet_idempotente(self):
        cs = ConstraintSet()
        cs.add_dirichlet(3, 0.5)
        cs.add_dirichlet(3, 0.5)  # no debe lanzar
        self.assertEqual(len(cs), 1)

    def test_dirichlet_inconsistente_levanta(self):
        cs = ConstraintSet()
        cs.add_dirichlet(3, 0.5)
        with self.assertRaises(ValueError):
            cs.add_dirichlet(3, 0.7)

    def test_slave_dofs_ordenados(self):
        cs = ConstraintSet()
        cs.add_dirichlet(5, 0.0)
        cs.add_dirichlet(1, 0.0)
        cs.add_dirichlet(3, 0.0)
        self.assertTrue(np.array_equal(cs.slave_dofs, [1, 3, 5]))

    def test_free_dofs_complemento(self):
        cs = ConstraintSet()
        cs.add_dirichlet(0, 0.0)
        cs.add_dirichlet(2, 0.0)
        free = cs.free_dofs(5)
        self.assertTrue(np.array_equal(free, [1, 3, 4]))


def _build_simply_supported_beam(*, settlement_at_right: float = 0.0):
    """Viga 2D Bernoulli-Euler simplemente apoyada, L=2, EI=1, EA=1.

    Apoyos: nodo 1 fijo en (ux, uy); nodo 2 fijo en uy con valor
    ``settlement_at_right`` (asentamiento si != 0). Las rotaciones quedan
    libres en ambos extremos.
    """
    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [2.0, 0.0])

    mat = LinearElastic1D(E=1.0)
    elem = Frame2DEuler(1, [n1, n2], mat, A=1.0, I=1.0)
    domain.add_element(elem)

    n1.fix_dof('ux', 0.0)
    n1.fix_dof('uy', 0.0)
    n2.fix_dof('uy', settlement_at_right)

    domain.generate_equation_numbers()
    return domain


class TestPrescribedSettlement(unittest.TestCase):
    """Asentamiento prescrito en viga simplemente apoyada — sin carga externa.

    Solución analítica: rotación rígida de la viga alrededor del apoyo izquierdo.
    El nodo derecho desciende ``δ``, ambos giros valen ``θ = δ / L`` (con signo
    según convención Reglas.md §5: nodo izquierdo gira en sentido anti-horario
    si δ < 0). Esfuerzos internos y reacciones nulos a redondeo.
    """

    def test_asentamiento_unitario_no_genera_esfuerzos(self):
        delta = -0.1  # asentamiento descendente del apoyo derecho
        L = 2.0
        domain = _build_simply_supported_beam(settlement_at_right=delta)
        assembler = Assembler(domain)
        F_ext = np.zeros(domain.total_dofs)

        U = LinearSolver(assembler).solve(F_ext)
        result = build_solve_result(domain, assembler, U, F_ext)

        n1, n2 = domain.get_node(1), domain.get_node(2)

        # Cinemática: descenso impuesto exacto, rotación rígida θ = δ/L.
        self.assertAlmostEqual(U[n2.dofs['uy']], delta, places=12)
        self.assertAlmostEqual(U[n2.dofs['ux']], 0.0, places=12)
        # Cuerpo rígido: ambos extremos rotan lo mismo.
        theta_expected = delta / L
        self.assertAlmostEqual(U[n1.dofs['rz']], theta_expected, places=10)
        self.assertAlmostEqual(U[n2.dofs['rz']], theta_expected, places=10)

        # Sin carga externa y sin deformación: reacciones cero a redondeo.
        # (La viga rota como sólido rígido, no se generan tensiones.)
        for node_id in (1, 2):
            for dof_name, val in result.reactions_by_node.get(node_id, {}).items():
                self.assertAlmostEqual(val, 0.0, places=8,
                                       msg=f"R[{node_id}][{dof_name}] = {val}")

        # Esfuerzos internos: cero a redondeo.
        ef = result.element_forces[1]
        np.testing.assert_allclose(ef.components['N'], 0.0, atol=1e-10)
        np.testing.assert_allclose(ef.components['V'], 0.0, atol=1e-10)
        np.testing.assert_allclose(ef.components['M'], 0.0, atol=1e-10)


class TestKReducedSymmetry(unittest.TestCase):
    """Tras reducir el sistema, ``K_red`` preserva la simetría de ``K``."""

    def test_simetria_armadura_2d(self):
        domain = Domain()
        n1 = domain.add_node(1, [0.0, 0.0])
        n2 = domain.add_node(2, [3.0, 4.0])
        n3 = domain.add_node(3, [6.0, 0.0])

        mat = LinearElastic1D(E=210.0e9)
        domain.add_element(Truss2D(1, [n1, n2], mat, A=1e-3))
        domain.add_element(Truss2D(2, [n2, n3], mat, A=1e-3))
        domain.add_element(Truss2D(3, [n1, n3], mat, A=1e-3))

        for dof in ('ux', 'uy'):
            n1.fix_dof(dof, 0.0)
        n3.fix_dof('uy', 0.0)

        domain.generate_equation_numbers()
        assembler = Assembler(domain)
        assembler.assemble_system()

        F = np.zeros(domain.total_dofs)
        F[n2.dofs['uy']] = -1e3

        K_red, _, _, _ = assembler.reduce(assembler.K_global.copy(), F)

        K_dense = K_red.toarray()
        # Simetría a redondeo: ‖K - Kᵀ‖∞ ≤ tol·‖K‖∞.
        asym = np.max(np.abs(K_dense - K_dense.T))
        ref = np.max(np.abs(K_dense))
        self.assertLess(asym / ref, 1e-12)


class TestReactionsBalance(unittest.TestCase):
    """Las reacciones reportadas equilibran exactamente la carga externa.

    Σ R + Σ F_ext = 0 en cada eje global; equilibrio de Newton.
    """

    def test_equilibrio_global_cantilever(self):
        domain = Domain()
        n1 = domain.add_node(1, [0.0, 0.0])
        n2 = domain.add_node(2, [1.0, 0.0])

        mat = LinearElastic1D(E=1.0)
        domain.add_element(Frame2DEuler(1, [n1, n2], mat, A=1.0, I=1.0))

        for dof in ('ux', 'uy', 'rz'):
            n1.fix_dof(dof, 0.0)

        domain.generate_equation_numbers()
        assembler = Assembler(domain)

        F = np.zeros(domain.total_dofs)
        Fy = -3.7
        F[n2.dofs['uy']] = Fy

        U = LinearSolver(assembler).solve(F.copy())
        result = build_solve_result(domain, assembler, U, F)

        # Σ Ry + Fy = 0; ux y rz solo en nodo 1.
        Ry_total = sum(rxn.get('uy', 0.0) for rxn in result.reactions_by_node.values())
        Rx_total = sum(rxn.get('ux', 0.0) for rxn in result.reactions_by_node.values())
        self.assertAlmostEqual(Ry_total + Fy, 0.0, places=10)
        self.assertAlmostEqual(Rx_total, 0.0, places=10)

        # Momento en el empotramiento equilibra el momento de Fy respecto a n1
        # (brazo = 1 m).
        Mz = result.reactions_by_node[1]['rz']
        # Convención §5 para Frame2D: comprobamos magnitud.
        self.assertAlmostEqual(abs(Mz), abs(Fy) * 1.0, places=8)


if __name__ == "__main__":
    unittest.main()
