"""Snap-through del von Mises 2-bar truss — benchmark analítico para ``Truss2DCorot + ArcLengthSolver``.

Geometría canónica (Crisfield, *Non-linear Finite Element Analysis of Solids
and Structures*, Vol. 1 §9.3; tambien Wriggers §5.1.1):

    (0,0) ──────── (L, h) ──────── (2L, 0)
      |              |                |
    apoyo           apex            apoyo
                  carga -P y

Apoyos articulados en los extremos; apex con ``u_x = 0`` (simetría) y carga
vertical hacia abajo. Con ``h / L ≪ 1`` (shallow) la curva carga–desplazamiento
del apex tiene un punto límite claro: arc-length lo atraviesa, control de
carga puro no.

**Convención del test**: ``w := -u_y`` (positivo hacia abajo), ``P > 0``
descendente. El elemento ``Truss2DCorot`` (ver ``solidum/elements/truss.py``)
usa engineering strain ``ε = (L_d - L_0)/L_0``, por lo que la solución cerrada
del equilibrio en el apex es:

    L_d(w) = sqrt(L² + (h - w)²)
    N(w)   = E·A·(L_d - L_0)/L_0       (tracción positiva)
    P(w)   = 2·N(w)·(w - h)/L_d        (equilibrio vertical; sign(w-h)
                                        invierte signo tras pasar la línea
                                        de apoyos)
           = 2·E·A·(L_0 - L_d)/L_0 · (h - w)/L_d        (forma equivalente)

Casos límite verificables a mano: ``P(0) = 0`` (barras sin alargar),
``P(h) = 0`` (configuración plana), ``P(2h) = 0`` (configuración invertida
simétrica).

Este test cierra el hueco "corot + arc-length cubierto sólo en sanidad" de
``docs/validacion_matriz.md`` §"Huecos prioritarios" (Fase B, 2026-05-19).
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.truss import Truss2DCorot
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import ArcLengthSolver


def _P_analytic(E: float, A: float, L: float, h: float, L0: float, w: float) -> float:
    """Curva carga-desplazamiento cerrada del shallow von Mises 2-bar."""
    L_d = math.sqrt(L * L + (h - w) ** 2)
    return 2.0 * E * A * (L0 - L_d) / L0 * (h - w) / L_d


class TestSnapThroughVonMisesTruss(unittest.TestCase):

    def setUp(self):
        self.L = 10.0
        self.h = 1.0          # ratio h/L = 0.1 (shallow canónico)
        self.E = 1.0e5
        self.A = 1.0
        self.L0 = math.sqrt(self.L ** 2 + self.h ** 2)

        dom = Domain()
        n_left = dom.add_node(1, [0.0, 0.0])
        n_apex = dom.add_node(2, [self.L, self.h])
        n_right = dom.add_node(3, [2.0 * self.L, 0.0])
        mat = Elastic1D(E=self.E)
        dom.add_element(Truss2DCorot(1, [n_left, n_apex], mat, A=self.A))
        dom.add_element(Truss2DCorot(2, [n_apex, n_right], mat, A=self.A))

        n_left.fix_dof('ux', 0.0); n_left.fix_dof('uy', 0.0)
        n_right.fix_dof('ux', 0.0); n_right.fix_dof('uy', 0.0)
        n_apex.fix_dof('ux', 0.0)
        dom.generate_equation_numbers(verbose=False)

        self.dom = dom
        self.asm = Assembler(dom)
        self.dof_uy = n_apex.dofs['uy']

        # Magnitud de referencia de carga (descendente). El pico analítico
        # del shallow truss h/L=0.1 está en P_pico ≈ 38 N; tomar F_ref = 100 N
        # garantiza que max_lambda=1 quede muy por encima del pico y el
        # solver tenga rango para trazar pre-pico + post-pico + invertido.
        self.F_ref_mag = 10.0 * self.E * self.A * (self.h / self.L) ** 2

    def test_arc_length_curve_matches_analytical_solution(self):
        """Cada punto convergido (u_y, λ·F_ref) coincide con la solución cerrada P(w).

        Comparación punto a punto contra ``_P_analytic`` con tolerancia
        absoluta = 1e-4·F_ref_mag = 0.01 N. Es ~10⁴× el residuo nodal
        esperado por convergencia (rtol_force=1e-8 sobre F_ref≈100 N), por
        lo que el test mide la corrección de la formulación corotacional +
        ensamblaje, no el ruido del solver.
        """
        F_ext = np.zeros(self.dom.total_dofs)
        F_ext[self.dof_uy] = -self.F_ref_mag

        history: list[tuple[float, float]] = []

        def probe(_step: int, U: np.ndarray, lam: float) -> None:
            history.append((float(U[self.dof_uy]), float(lam)))

        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        solver = ArcLengthSolver(
            self.asm, convergence=conv,
            max_lambda=1.0, initial_dl=0.05, max_steps=200, max_iter=30,
        )
        solver.solve(F_ext, step_callback=probe)

        self.assertGreater(len(history), 5,
            f"Arc-length sólo registró {len(history)} pasos; insuficiente "
            "para trazar la curva.")

        abs_tol = 1.0e-4 * self.F_ref_mag

        for u_y, lam in history:
            w = -u_y
            P_real = lam * self.F_ref_mag
            P_anal = _P_analytic(self.E, self.A, self.L, self.h, self.L0, w)
            self.assertLess(
                abs(P_real - P_anal), abs_tol,
                f"En w={w:.4f}: P_arc={P_real:.4e} vs P_anal={P_anal:.4e}, "
                f"|Δ|={abs(P_real - P_anal):.2e} > tol={abs_tol:.2e}.",
            )

    def test_arc_length_traverses_limit_point_and_post_snap_region(self):
        """Los puntos trazados cubren la rama inestable (0 < w < h) y la rama
        invertida (w > h con λ negativo durante el cruce).

        Si el arc-length saltara la región problemática sin trazarla, podría
        registrar sólo configuraciones cuasi-elásticas y configuraciones
        invertidas con `w ≫ h`, sin puntos en el medio. La región inestable
        del shallow truss h/L=0.1 cae en w ∈ (0, h) — específicamente el pico
        analítico está cerca de w ≈ 0.4. Se exige al menos un paso en cada
        una de las tres regiones físicamente distintas.
        """
        F_ext = np.zeros(self.dom.total_dofs)
        F_ext[self.dof_uy] = -self.F_ref_mag

        history: list[tuple[float, float]] = []

        def probe(_step: int, U: np.ndarray, lam: float) -> None:
            history.append((float(U[self.dof_uy]), float(lam)))

        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        solver = ArcLengthSolver(
            self.asm, convergence=conv,
            max_lambda=1.0, initial_dl=0.05, max_steps=200, max_iter=30,
        )
        solver.solve(F_ext, step_callback=probe)

        ws = [-u_y for u_y, _ in history]
        lams = [lam for _, lam in history]

        n_elastic = sum(1 for w in ws if 0.0 < w < 0.2 * self.h)
        n_unstable = sum(1 for w in ws if 0.2 * self.h <= w <= 0.8 * self.h)
        n_inverted = sum(1 for w in ws if w > self.h)
        n_lambda_neg = sum(1 for lam in lams if lam < 0.0)

        self.assertGreaterEqual(n_elastic, 1,
            f"Sin pasos en la rama (cuasi)elástica (0 < w < 0.2·h): {n_elastic}.")
        self.assertGreaterEqual(n_unstable, 1,
            f"Sin pasos en la rama inestable (0.2·h ≤ w ≤ 0.8·h): {n_unstable}; "
            "arc-length saltó la región del pico.")
        self.assertGreaterEqual(n_inverted, 1,
            f"Sin pasos en la rama invertida (w > h): {n_inverted}.")
        self.assertGreaterEqual(n_lambda_neg, 1,
            f"Sin pasos con λ < 0: {n_lambda_neg}; arc-length no detectó "
            "el cambio de signo en el cruce de la línea de apoyos.")


if __name__ == '__main__':
    unittest.main()
