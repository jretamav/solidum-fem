"""Tracción uniaxial con discontinuidad embebida hasta la separación completa.

Referencia
----------
Retama Velasco, J. (2010). *Formulation and Approximation to Problems in
Solids by Embedded Discontinuity Models*, UNAM: cinemática KOS (Cap. 2 y 5),
equilibrio en Γ_d con ``l_d = A/h`` (Cap. 6, ecs. 6.24-6.25) y modelo
cohesivo de daño isótropo (Cap. 3), que sigue a Alfaiate, Wells y Sluys
(2002), *Eng. Fract. Mech.* 69, ecs. 8-17.

Concepto
--------
Probeta rectangular L×W de espesor ``t_h`` partida por la diagonal en dos
``CST_Embedded2D``. Un borde fijo en ``x`` y el opuesto con desplazamiento
impuesto ``δ``: esfuerzo uniaxial ``σ_xx`` uniforme en esfuerzo plano. Ambos
elementos se activan con la grieta vertical por el centroide; en los dos el
lado opuesto al nodo solitario es vertical, así que ``cos(θ − α) = 1``,
``l_d = A/h = W/2`` y la fila del salto da ``t_n = σ_xx`` exacto.

Solución exacta
---------------
El modo del salto separa rígidamente el nodo solitario del resto, así que el
alargamiento total es el elástico del volumen más la apertura ``w``::

    δ = f(w)·L/E + w,      F = f(w)·W·t_h

con ``f`` la envolvente de ablandamiento de la bibliografía, sin la rama
elástica de penalización (su efecto es ``O(κ_0/w_c) ~ 1e-5``):

    lineal:       f(w) = σ_t0·(1 − w/w_c),   w_c = 2·G_F/σ_t0
    exponencial:  f(w) = σ_t0·exp(−σ_t0·w/G_F)

Sin retroceso (``1 + f′·L/E > 0``) el control por desplazamiento atraviesa la
rama descendente. Con la separación completa toda la energía elástica se ha
recuperado y el trabajo externo es la energía disipada, ``G_F·W·t_h``.

Qué blinda
----------
Hasta la revisión del 2026-09-23 ningún test hacía que un solver atravesara
el ablandamiento del cohesivo, y el material devolvía en casi toda esa rama
una tangente de signo contrario a la real. Aquí intervienen a la vez la
activación de Rankine, ``l_d``, el Newton local del salto, la condensación y
la ley cohesiva, con parámetros físicos distintos de 1 (espesor incluido).
"""
from __future__ import annotations

import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from solidum.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic
from solidum.core.domain import Domain
from solidum.elements.solid_2d.embedded_cst import CST_Embedded2D
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver


E, NU = 30.0e9, 0.2
SIGMA_T0, G_F = 3.0e6, 100.0
K_E = 1.0e15                        # penalización: κ_0/w_c ≈ 4.5e-5
L, W, TH = 0.1, 0.05, 0.02          # m
DELTA_MAX = 8.0e-5                  # > w_c lineal: separación completa
N_STEPS = 75                        # el pico no cae en un paso exacto


def envelope(w: float, softening: str) -> float:
    if softening == 'linear':
        w_c = 2.0 * G_F / SIGMA_T0
        return SIGMA_T0 * max(1.0 - w / w_c, 0.0)
    return SIGMA_T0 * math.exp(-SIGMA_T0 * w / G_F)


def exact_force(delta: float, softening: str) -> float:
    """Fuerza de la solución exacta en la rama de ablandamiento: resuelve
    ``δ = f(w)·L/E + w`` por bisección (``g`` es monótona sin retroceso)."""
    g = lambda w: envelope(w, softening) * L / E + w - delta
    lo, hi = 0.0, delta
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return envelope(0.5 * (lo + hi), softening) * W * TH


def run(softening: str):
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    n3 = dom.add_node(3, [L, W])
    n4 = dom.add_node(4, [0.0, W])
    bulk = Elastic2D(E=E, nu=NU, hypothesis='plane_stress')
    coh = CohesiveDamageIsotropic(sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_E, softening=softening)
    dom.add_element(CST_Embedded2D(1, [n1, n2, n3], bulk, coh, thickness=TH))
    dom.add_element(CST_Embedded2D(2, [n1, n3, n4], bulk, coh, thickness=TH))
    # Rigidez de modo II nula (Retama 2010, p. 67): con la grieta cruzando
    # toda la probeta, la parte derecha deslizaría en y sin resistencia. Se
    # fija uy en las dos esquinas inferiores; la contracción lateral sigue
    # libre (nodos 3 y 4) y el esfuerzo es uniaxial.
    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
    n2.fix_dof('uy', 0.0)
    n4.fix_dof('ux', 0.0)
    n2.fix_dof('ux', DELTA_MAX)
    n3.fix_dof('ux', DELTA_MAX)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    solver = NonlinearSolver(asm, num_steps=N_STEPS, adaptive=False)

    right = [n2.dofs['ux'], n3.dofs['ux']]
    history = []

    def callback(step, U, lam):
        _, F_int = asm.assemble_non_linear_system(U)
        cracked = all(el.discontinuity_state is not None for el in dom.elements.values())
        history.append((lam * DELTA_MAX, float(F_int[right].sum()), cracked))

    solver.solve(np.zeros(dom.total_dofs), step_callback=callback)
    return history, dom


class TestUniaxialSofteningThroughSeparation(unittest.TestCase):

    def _check(self, softening: str, rtol: float):
        history, dom = run(softening)
        self.assertEqual(len(history), N_STEPS)

        # Rama elástica exacta antes de la activación.
        pre = [(d, F) for d, F, c in history if not c]
        self.assertGreater(len(pre), 0)
        for d, F in pre:
            self.assertAlmostEqual(F, E * d / L * W * TH, delta=1e-9 * E * d / L * W * TH)

        # Rama de ablandamiento frente a la solución exacta.
        post = [(d, F) for d, F, c in history if c]
        self.assertGreater(len(post), 50)
        F_peak = SIGMA_T0 * W * TH
        for d, F in post:
            with self.subTest(softening=softening, delta=d):
                self.assertAlmostEqual(F, exact_force(d, softening), delta=rtol * F_peak)

        # Grieta vertical por el centroide con l_d = A/h = W/2 en ambos.
        for el in dom.elements.values():
            ds = el.discontinuity_state
            self.assertAlmostEqual(abs(ds.normal[0]), 1.0, places=12)
            self.assertAlmostEqual(ds.l_d, 0.5 * W, places=14)
        return history

    def test_linear_softening(self):
        history = self._check('linear', rtol=1e-4)
        # Separación completa: fuerza nula y trabajo externo = G_F·W·t_h.
        self.assertAlmostEqual(history[-1][1], 0.0, delta=1e-9 * SIGMA_T0 * W * TH)
        d = np.array([0.0] + [h[0] for h in history])
        F = np.array([0.0] + [h[1] for h in history])
        work = float(np.sum(0.5 * (F[1:] + F[:-1]) * np.diff(d)))
        self.assertAlmostEqual(work, G_F * W * TH, delta=1e-2 * G_F * W * TH)

    def test_exponential_softening(self):
        self._check('exponential', rtol=1e-3)


if __name__ == '__main__':
    unittest.main()
