"""Conservación del momento lineal en `CentralDifferenceSolver`.

Análogo dinámico-explícito del test de conservación de energía Newmark
(Tanda 6, `test_newmark_energy_conservation.py`): para un sistema libre
sin fuerzas externas ni amortiguamiento, el integrador de diferencias
centradas debe conservar el momento lineal total exactamente
(Belytschko-Liu-Moran 2014 §6.2.1, Géradin & Rixen 2014 §7.4).

    P(t) = 1ᵀ · M · u̇(t)        debe ser independiente de t

con `1` el vector con 1 en cada DOF correspondiente al modo de
traslación rígida (en este test, todos los DOFs en x). La conservación
es exacta a precisión máquina para CD sobre cualquier M SPD (consistente
o lumped) mientras `Δt < Δt_crit`.

Configuración: cadena de 4 truss horizontales con **todos los nodos
libres en ux** (no hay apoyo axial — el sistema tiene modo rígido en
x), uy fijo en cada nodo para reducir a movimiento axial puro.
Excitación inicial: velocidad concentrada en un nodo.

Cualquier deriva creciente del momento revelaría asimetría en `K`,
asignación errónea de los DOFs del modo rígido, o un bug en la
actualización del integrador.
"""
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


def _build_free_axial_chain(n_elem: int = 4, L: float = 1.0,
                              E: float = 1.0e3, A: float = 1.0, rho: float = 1.0):
    """Cadena de `n_elem` truss horizontales con **todos los nodos libres
    en ux** (sin apoyo axial → modo rígido en x) y uy fijo en todos los
    nodos (movimiento axial puro)."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * L / n_elem, 0.0])
             for i in range(n_elem + 1)]
    mat = Elastic1D(E=E, density=rho)
    for k in range(n_elem):
        dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]], mat, A=A))
    for n in nodes:
        n.fix_dof('uy', 0.0)  # axial puro
    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), nodes


def _linear_momentum(M, v, ux_dofs):
    """P = 1ᵀ·M·v restringido a los DOFs de translación en x."""
    e_x = np.zeros(v.shape[0])
    e_x[ux_dofs] = 1.0
    return float(e_x @ (M @ v))


class TestCentralDifferenceLinearMomentumConservation(unittest.TestCase):

    def _run_and_check(self, lumping: str, tol_rel: float):
        n_elem = 4
        L, E, A, rho = 1.0, 1.0e3, 1.0, 1.0
        dom, asm, nodes = _build_free_axial_chain(n_elem, L=L, E=E, A=A, rho=rho)
        ux_dofs = [n.dofs['ux'] for n in nodes]

        # Velocidad inicial concentrada en el nodo central (excitación
        # asimétrica) → momento lineal inicial no nulo.
        v0 = np.zeros(dom.total_dofs)
        v0[nodes[n_elem // 2].dofs['ux']] = 1.0
        u0 = np.zeros(dom.total_dofs)

        # Estimación conservadora del Δt crítico: ω_max ≈ 2·√(E·A/(m·L²))
        # para cadena con masa lumped m = ρ·A·L. Δt_crit = 2/ω_max.
        m_nodal = rho * A * L / n_elem  # masa lumped por mid-nodo
        omega_max_est = 2.0 * np.sqrt(E * A / (m_nodal * (L / n_elem) ** 2))
        dt = 0.3 * (2.0 / omega_max_est)  # margen amplio
        t_end = 200.0 * dt  # muchos pasos

        solver = CentralDifferenceSolver(
            asm, t_end=t_end, dt=dt,
            u0=u0, u0_dot=v0,
            F_func=None,
            lumping=lumping,
        )
        result = solver.solve()

        # Ensamblar M ground-truth (mismo lumping).
        M_global = asm.assemble_mass_matrix(lumping=lumping).toarray()

        P_history = np.array([
            _linear_momentum(M_global, result.udot_history[:, k], ux_dofs)
            for k in range(result.udot_history.shape[1])
        ])
        P0 = P_history[0]
        self.assertGreater(abs(P0), 0.0,
            "Momento inicial nulo — setup incorrecto.")

        drift = float(np.max(np.abs(P_history - P0)) / abs(P0))
        self.assertLess(drift, tol_rel,
            f"CentralDifference ({lumping}) no conservó momento: "
            f"drift relativo = {drift:.3e} > tol = {tol_rel:.0e}. "
            f"P0 = {P0:.4e}, P_range = [{P_history.min():.4e}, "
            f"{P_history.max():.4e}] sobre {result.udot_history.shape[1]} pasos.")

    def test_lumped_mass(self):
        """Masa lumped: única discretización soportada por `CentralDifferenceSolver`
        (ADR 0009 fase 2 — el esquema explícito requiere M diagonal)."""
        self._run_and_check(lumping='lumped', tol_rel=1.0e-10)


if __name__ == '__main__':
    unittest.main()
