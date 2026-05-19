"""Conservación de energía en `NewmarkSolver` trapezoidal (β=1/4, γ=1/2).

Diagnóstico físico fuerte: la familia Newmark con `β=1/4, γ=1/2`
(regla del trapezoide / aceleración promedio) es **exactamente
conservativa** sobre sistemas lineales sin amortiguamiento (Hughes 1987
§9.1.5, Géradin & Rixen 2014 §7.1.2). Toda inconsistencia entre la
matriz de masa, la matriz de rigidez ensambladas y el integrador
temporal rompe esta propiedad — el test capta de un golpe múltiples
subsistemas.

Configuración: cadena de 3 truss horizontales empotrada en un extremo,
condición inicial de desplazamiento no-nulo y velocidad nula
(``E_kin(0)=0``, ``E_pot(0)≠0``), integración por muchos periodos sin
fuerza externa ni amortiguamiento.

Energía total a verificar:

.. code-block:: text

    E(t) = ½·u̇ᵀ·M·u̇ + ½·uᵀ·K·u

debe permanecer constante a precisión máquina relativa (ε_rel < 1e-10)
durante todo el análisis. Cualquier desviación creciente con el tiempo
indica disipación numérica o aporte espurio.

Se cubren ambas discretizaciones de masa (consistente y lumped), porque
la conservación del trapezoide es independiente del lumping mientras M
sea SPD.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.domain import Domain
from fenix.elements.truss import Truss2D
from fenix.materials.elastic import Elastic1D
from fenix.math.assembly import Assembler
from fenix.math.solvers import NewmarkSolver


def _build_axial_chain(n_elem: int = 3, L: float = 1.0,
                        E: float = 1.0e3, A: float = 1.0, rho: float = 1.0):
    """Cadena de `n_elem` truss horizontales, extremo izquierdo empotrado y
    libre del lado derecho. Movimiento axial puro (ux libre, uy fijo en
    todos los nodos libres)."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * L / n_elem, 0.0])
             for i in range(n_elem + 1)]
    mat = Elastic1D(E=E, density=rho)
    for k in range(n_elem):
        dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]], mat, A=A))

    # Empotramiento en el primer nodo.
    nodes[0].fix_dof('ux', 0.0); nodes[0].fix_dof('uy', 0.0)
    # Restringir uy en los demás para que el problema sea axial puro.
    for n in nodes[1:]:
        n.fix_dof('uy', 0.0)

    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), nodes


def _compute_energy_history(result, M, K):
    """E(t) = ½·u̇ᵀ·M·u̇ + ½·uᵀ·K·u para cada paso del integrador."""
    n_steps = result.u_history.shape[1]
    E = np.zeros(n_steps)
    for k in range(n_steps):
        u = result.u_history[:, k]
        v = result.udot_history[:, k]
        E[k] = 0.5 * v @ (M @ v) + 0.5 * u @ (K @ u)
    return E


class TestNewmarkTrapezoidalEnergyConservation(unittest.TestCase):
    """Conservación de energía del trapezoide Newmark."""

    def _run_and_check(self, lumping: str, tol_rel: float):
        n_elem = 3
        L, E, A, rho = 1.0, 1.0e3, 1.0, 1.0
        dom, asm, nodes = _build_axial_chain(n_elem, L=L, E=E, A=A, rho=rho)

        # Periodo natural aproximado del primer modo de la cadena fijo-libre
        # (modo axial fundamental): ω₁ ≈ (π/2L)·√(E/ρ) → T₁ = 4L·√(ρ/E).
        T1 = 4.0 * L * np.sqrt(rho / E)
        dt = T1 / 100.0
        n_periods = 20
        t_end = n_periods * T1

        # Condición inicial: desplazamiento triangular (perfil estático bajo
        # carga puntual en el extremo libre, no proporcional al primer modo
        # para garantizar excitación de múltiples modos).
        u0 = np.zeros(dom.total_dofs)
        for i, n in enumerate(nodes[1:], start=1):
            u0[n.dofs['ux']] = i / n_elem  # 1/3, 2/3, 1 sobre los tres libres
        u0_dot = np.zeros(dom.total_dofs)

        solver = NewmarkSolver(
            asm, t_end=t_end, dt=dt,
            beta=0.25, gamma=0.5,        # trapezoide canónico
            rayleigh=None,               # sin amortiguamiento
            u0=u0, u0_dot=u0_dot,
            F_func=None,                 # vibración libre
            lumping=lumping,
        )
        result = solver.solve()

        # Ensamblar M y K independientemente del solver para verificar la
        # energía con las matrices "ground truth".
        asm.assemble_system()
        K_global = asm.K_global.toarray()
        M_global = asm.assemble_mass_matrix(lumping=lumping).toarray()

        E_hist = _compute_energy_history(result, M_global, K_global)
        E0 = E_hist[0]
        self.assertGreater(E0, 0.0, "Energía inicial nula — setup incorrecto.")

        # Verificar conservación durante todo el análisis. El trapezoide
        # NO disipa pero acumula error de cancelación catastrófica
        # entre los términos cinético y potencial cuando se restan
        # cantidades de orden uno. Toleramos hasta tol_rel·E0.
        E_min = E_hist.min()
        E_max = E_hist.max()
        drift = max(E0 - E_min, E_max - E0) / E0
        self.assertLess(drift, tol_rel,
            f"Trapezoide ({lumping}) no conservó energía: "
            f"drift={drift:.3e} > tol={tol_rel:.0e}. "
            f"E0={E0:.4e}, E_min={E_min:.4e}, E_max={E_max:.4e} "
            f"sobre {n_periods} periodos.")

    def test_consistent_mass(self):
        """Masa consistente: conservación exacta a precisión máquina."""
        self._run_and_check(lumping='consistent', tol_rel=1.0e-10)

    def test_lumped_mass(self):
        """Masa lumped (HRZ): la conservación del trapezoide es
        independiente del lumping mientras M sea SPD (lo es por construcción).
        """
        self._run_and_check(lumping='lumped', tol_rel=1.0e-10)


if __name__ == '__main__':
    unittest.main()
