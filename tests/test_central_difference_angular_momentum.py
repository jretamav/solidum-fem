"""Conservación del momento angular en `CentralDifferenceSolver`.

Pareja del test de momento lineal (Tanda 10): sobre un sistema libre
sin fuerzas externas ni amortiguamiento, el integrador BLM debe
conservar exactamente el momento angular respecto a un punto fijo
(Belytschko-Liu-Moran §6.2.1).

Para un sistema de N partículas con masa lumped y posiciones de
referencia ``r_i = (x_i, y_i)``:

    L(t) = Σ_i M_ii · (x_i · v_iy(t) − y_i · v_ix(t))

debe permanecer constante a precisión máquina mientras las fuerzas
internas sean centrales — para truss articulado en 2D esto se cumple
exactamente: las parejas acción-reacción de cada barra están alineadas
con el eje de la barra y producen momento neto cero respecto a
cualquier punto.

Configuración: triángulo equilátero de 3 nodos centrado en el origen
con 3 truss conectando los pares (estructura rígida en su
conectividad). Todos los DOFs libres → 6 DOFs totales, con 3 modos
rígidos (2 traslaciones + 1 rotación) + 3 modos de deformación
elástica. Excitación inicial: rotación rígida pura ``ω·k × r`` +
pequeña perturbación radial (modo de respiración) para excitar
también vibración interna. El momento angular del modo rígido
rotacional es ``L_0 = ω·I``; el modo radial no aporta a `L` (radial
contra el centro).
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.solvers import CentralDifferenceSolver


def _build_equilateral_triangle(E: float = 1.0e3, A: float = 1.0, rho: float = 1.0):
    """Triángulo equilátero centrado en el origen con 3 truss interconectados."""
    dom = Domain()
    r = 1.0
    coords = [
        (r * math.cos(0.0),               r * math.sin(0.0)),
        (r * math.cos(2.0 * math.pi / 3.0), r * math.sin(2.0 * math.pi / 3.0)),
        (r * math.cos(4.0 * math.pi / 3.0), r * math.sin(4.0 * math.pi / 3.0)),
    ]
    nodes = [dom.add_node(i + 1, list(xy)) for i, xy in enumerate(coords)]
    mat = Elastic1D(E=E, density=rho)
    # 3 truss conectando los pares (1-2, 2-3, 3-1).
    dom.add_element(Truss2D(1, [nodes[0], nodes[1]], mat, A=A))
    dom.add_element(Truss2D(2, [nodes[1], nodes[2]], mat, A=A))
    dom.add_element(Truss2D(3, [nodes[2], nodes[0]], mat, A=A))
    # Todos los DOFs libres (sistema completamente libre en plano).
    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), nodes


def _angular_momentum(M_diag: np.ndarray, v: np.ndarray, nodes) -> float:
    """L = Σ M_ii · (x_i · v_iy − y_i · v_ix), masas lumped diagonales,
    posiciones de referencia."""
    L = 0.0
    for n in nodes:
        x, y = n.coordinates[0], n.coordinates[1]
        dof_ux = n.dofs['ux']
        dof_uy = n.dofs['uy']
        L += M_diag[dof_ux] * (x * v[dof_uy] - y * v[dof_ux])
    return float(L)


class TestCentralDifferenceAngularMomentumConservation(unittest.TestCase):

    def test_lumped_mass_triangle_with_rotation_plus_breathing(self):
        """Triángulo equilátero libre con rotación rígida + respiración
        elástica: el momento angular se conserva con drift < 1e-10.

        Setup: ω = 1.0 produce un L_0 = ω · Σ m_i·|r_i|². La perturbación
        radial excita el modo de "respiración" (símetrico, no contribuye
        a L). Las fuerzas internas axiales del truss no aportan momento
        respecto al centro (parejas acción-reacción colineales con el eje
        de la barra).
        """
        E, A, rho = 1.0e3, 1.0, 1.0
        dom, asm, nodes = _build_equilateral_triangle(E=E, A=A, rho=rho)

        # Excitación inicial: rotación rígida ω + perturbación radial.
        omega = 1.0
        breath = 0.1  # pequeña velocidad radial saliente
        u0 = np.zeros(dom.total_dofs)
        v0 = np.zeros(dom.total_dofs)
        for n in nodes:
            x, y = n.coordinates[0], n.coordinates[1]
            r = math.hypot(x, y)
            # rotación + breathing
            v0[n.dofs['ux']] = -omega * y + breath * x / r
            v0[n.dofs['uy']] = omega * x + breath * y / r

        # Δt sub-crítico: ω_max ~ 2·√(E·A/(m·L²)), L = √3 (lado del
        # triángulo con r=1), m_nodal lumped = ρ·A·(2·L)·0.5 = ρ·A·L.
        L_bar = math.sqrt(3.0)
        m_nodal = rho * A * L_bar
        omega_max_est = 2.0 * math.sqrt(E * A / (m_nodal * L_bar ** 2))
        dt = 0.3 * (2.0 / omega_max_est)
        t_end = 200.0 * dt

        solver = CentralDifferenceSolver(
            asm, t_end=t_end, dt=dt,
            u0=u0, u0_dot=v0,
            F_func=None,
            lumping='lumped',
        )
        result = solver.solve()

        # Masa diagonal ground-truth (independiente del solver).
        M_global = asm.assemble_mass_matrix(lumping='lumped').toarray()
        M_diag = np.diag(M_global)

        L_history = np.array([
            _angular_momentum(M_diag, result.udot_history[:, k], nodes)
            for k in range(result.udot_history.shape[1])
        ])
        L0 = L_history[0]
        self.assertGreater(abs(L0), 0.0,
            "Momento angular inicial nulo — setup incorrecto.")

        drift = float(np.max(np.abs(L_history - L0)) / abs(L0))
        self.assertLess(drift, 1.0e-10,
            f"CentralDifference no conservó momento angular: "
            f"drift relativo = {drift:.3e}. "
            f"L0 = {L0:.4e}, L_range = [{L_history.min():.4e}, "
            f"{L_history.max():.4e}] sobre {result.udot_history.shape[1]} pasos.")


if __name__ == '__main__':
    unittest.main()
