"""``CentralDifferenceSolver(nonlinear=True)`` con material con historia.

Regresión de la auditoría 2026-09-22: el solver explícito no comiteaba el
estado de los elementos en ningún paso, así que el return mapping se
reevaluaba siempre desde el estado virgen y el oscilador elastoplástico
conservaba la energía como si el material fuera elástico no lineal.
"""
import math
import unittest

import numpy as np

import solidum  # noqa: F401 — autodiscover
from solidum.core.domain import Domain
from solidum.elements.truss import Truss2D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.assembly import Assembler
from solidum.math.solvers import CentralDifferenceSolver, NewtonNewmarkSolver

# Oscilador 1 GDL: K_red = EA/L = 25, M_red (lumped) = ρAL/2 = 1 → ω = 5.
E, RHO, A, L, OMEGA = 25.0, 2.0, 1.0, 1.0, 5.0


def _oscillator(sigma_y=5.0, H=5.0):
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    mat = Elastoplastic1D(E=E, sigma_y=sigma_y, H=H, density=RHO)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, n2.dofs["ux"]


class TestCentralDifferenceHistory(unittest.TestCase):

    def test_plasticity_accumulates_and_dissipates(self):
        T = 2.0 * math.pi / OMEGA
        dt = T / 400.0
        dom, dof = _oscillator()
        # Velocidad inicial grande: energía cinética ½·1·6² = 18 ≫ energía
        # elástica en fluencia ½·25·0.2² = 0.5, así que plastifica en ambos
        # sentidos en cada ciclo y disipa hasta que el endurecimiento lo
        # devuelve al régimen elástico.
        v0 = np.zeros(dom.total_dofs)
        v0[dof] = 6.0
        res = CentralDifferenceSolver(
            Assembler(dom), t_end=5.0 * T, dt=dt, u0_dot=v0, nonlinear=True,
        ).solve()
        u = res.u_history[dof, :]

        alpha = dom.elements[1].state.vars[0]["alpha"]
        self.assertGreater(alpha, 0.0, "el estado committed no acumula deformación plástica")
        n_per = int(round(T / dt))
        amp_first = float(np.max(np.abs(u[:n_per])))
        amp_last = float(np.max(np.abs(u[-n_per:] - np.mean(u[-n_per:]))))
        # Sin commit el material es elástico no lineal y la oscilación
        # conserva la energía: la amplitud del último ciclo igualaría la
        # del primero. Con disipación plástica cae de forma clara.
        self.assertLess(amp_last, 0.5 * amp_first)

    def test_matches_newton_newmark_in_plastic_regime(self):
        """Mismo oscilador integrado con Newton-Newmark (que sí comitea): las
        trayectorias deben coincidir dentro del error de discretización."""
        T = 2.0 * math.pi / OMEGA
        dt = T / 800.0
        dom_cd, dof = _oscillator()
        u0 = np.zeros(dom_cd.total_dofs); u0[dof] = 0.5
        u_cd = CentralDifferenceSolver(
            Assembler(dom_cd), t_end=2.0 * T, dt=dt, u0=u0, nonlinear=True,
        ).solve().u_history[dof, :]

        dom_nn, _ = _oscillator()
        u_nn = NewtonNewmarkSolver(
            Assembler(dom_nn), t_end=2.0 * T, dt=dt, u0=u0, lumping="lumped",
        ).solve().u_history[dof, :]

        err = float(np.max(np.abs(u_cd - u_nn)))
        self.assertLess(err, 0.02, f"|u_CD − u_NN|_max = {err:.3e}")
        self.assertGreater(dom_cd.elements[1].state.vars[0]["alpha"], 0.0)


if __name__ == '__main__':
    unittest.main()
