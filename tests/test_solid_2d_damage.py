"""Tests de integración del pipeline Quad4 + IsotropicDamage2D + NonlinearSolver.

Acceptance de ``docs/specs/IsotropicDamage2D.md`` sección ``acceptance.integration``.

Hasta ahora la cobertura del modelo de daño 2D era exclusivamente unitaria
(`test_materials_unit.py · TestIsotropicDamage2D`). Estos tests validan que la
tangente algorítmica consistente recupera convergencia rápida del Newton global
en presencia de daño activo, y que el cableado solver-material-elemento no
introduce regresiones con ninguna de las tres medidas de deformación
equivalente (``equivalent_strain``): la energética (por defecto, tangente
simétrica, ``IS_SYMMETRIC = True`` en la instancia), Mazars y von Mises
modificado (tangente no simétrica, ``IS_SYMMETRIC = False``).

Problema: un Quad4 unitario en tensión plana con tracción uniaxial de esfuerzo
por desplazamiento prescrito ``δ`` en el borde derecho (``u_y`` libre salvo en
el borde inferior). El campo es uniforme con ``ε = [δ, −ν·δ, 0]`` para
cualquier ``d`` (``σ_yy = (1−d)·(C_e·ε)_yy = 0``) y las tres medidas dan
``ε_eq = δ``, así que ``d(δ)`` y ``σ_xx = (1−d)·E·δ`` son analíticos.

Historia (2026-09-29): hasta esa fecha los tests con daño activo aplicaban
fuerza, ``F = 0,5·E·ε_objetivo`` o ``0,6·E·ε_objetivo`` con
``ε_objetivo = 4…6·κ₀``, por encima del pico ``σ ≈ E·κ₀`` de la ley
exponencial con ``α·κ₀ ≪ 1``: el equilibrio sólo existía con el daño
saturado en ``DAMAGE_MAX`` (rigidez residual), y el paso adaptativo ocultaba
cualquier iteración de más. Se pasaron a desplazamiento prescrito, que
recorre el ablandamiento sin saturar y admite comparación analítica.
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4
from solidum.materials.damage_2d import IsotropicDamage2D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver

MEDIDAS = (('energy', {}), ('mazars', {}), ('modified_von_mises', {'k': 10.0}))


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

    E = 2.0e5
    NU = 0.3
    KAPPA_0 = 1.0e-4
    ALPHA = 100.0

    def _modelo(self, medida='energy', **kw):
        material = IsotropicDamage2D(
            E=self.E, nu=self.NU, kappa_0=self.KAPPA_0, alpha=self.ALPHA,
            hypothesis='plane_stress', equivalent_strain=medida, **kw,
        )
        domain, element, nodes = _build_uniaxial_quad4_domain(material)
        return domain, element, nodes, Assembler(domain)

    @staticmethod
    def _prescribir(nodes, delta):
        """``u_x = δ`` en el borde derecho (nodos 2 y 3)."""
        nodes['n2'].fix_dof('ux', delta)
        nodes['n3'].fix_dof('ux', delta)

    def _d_analitico(self, delta):
        """Ley exponencial con ``κ = ε_eq = δ`` (tracción uniaxial de esfuerzo)."""
        if delta <= self.KAPPA_0:
            return 0.0
        return 1.0 - (self.KAPPA_0 / delta) * math.exp(-self.ALPHA * (delta - self.KAPPA_0))

    def _conv(self):
        return ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)

    def test_elastic_regime_below_threshold(self):
        """Carga por debajo del umbral: sin daño, respuesta elástica pura."""
        domain, element, nodes, assembler = self._modelo()
        # ε_xx target = 0.5·κ_0 ⇒ σ_xx target = E·0.5·κ_0 (plane stress uniaxial)
        eps_target = 0.5 * self.KAPPA_0
        F_ext = np.zeros(domain.total_dofs)
        F_ext[nodes['n2'].dofs['ux']] = 0.5 * self.E * eps_target
        F_ext[nodes['n3'].dofs['ux']] = 0.5 * self.E * eps_target

        solver = NonlinearSolver(assembler, convergence=self._conv(), num_steps=2)
        U = solver.solve(F_ext)

        self.assertAlmostEqual(U[nodes['n2'].dofs['ux']], eps_target, places=10)
        self.assertEqual(element.state.vars[0]['damage'], 0.0)

    def test_damage_active_converges_in_few_iter(self):
        """Ablandamiento hasta ``δ = 4·κ₀`` (``d ≈ 0,76``) en 5 pasos: con la
        tangente consistente el Newton converge en ≤ 6 iteraciones por paso
        sin bisección (``adaptive=False``: el solver lanza si no converge), y
        ``σ_xx`` coincide con ``(1−d)·E·δ``. Con las tres medidas."""
        delta = 4.0 * self.KAPPA_0
        d_esperado = self._d_analitico(delta)
        self.assertTrue(0.5 < d_esperado < 0.9)
        for medida, kw in MEDIDAS:
            with self.subTest(medida=medida):
                domain, element, nodes, assembler = self._modelo(medida, **kw)
                self._prescribir(nodes, delta)
                solver = NonlinearSolver(assembler, convergence=self._conv(), num_steps=5,
                                         max_iter=6, adaptive=False)
                U = solver.solve(np.zeros(domain.total_dofs))

                for st in element.state.vars:
                    self.assertAlmostEqual(st['kappa'] / delta, 1.0, places=10)
                    self.assertAlmostEqual(st['damage'], d_esperado, places=10)
                sigma = element.compute_gauss_state(U)['stress']
                np.testing.assert_allclose(sigma[:, 0], (1.0 - d_esperado) * self.E * delta,
                                           rtol=1e-8)
                np.testing.assert_allclose(sigma[:, 1:], 0.0, atol=1e-8 * self.E * delta)

    def test_unloading_no_damage_increment(self):
        """Tras carga dañadora, descarga: κ y d no cambian."""
        domain, element, nodes, assembler = self._modelo()
        # Paso 1: cargar dañando
        self._prescribir(nodes, 5.0 * self.KAPPA_0)
        NonlinearSolver(assembler, convergence=self._conv(), num_steps=5,
                        max_iter=8).solve(np.zeros(domain.total_dofs))

        state_after_load = dict(element.state.vars[0])
        kappa_loaded = state_after_load['kappa']
        d_loaded = state_after_load['damage']
        self.assertAlmostEqual(d_loaded, self._d_analitico(5.0 * self.KAPPA_0), places=10)

        # Paso 2: descargar a un nivel inferior. Nuevo solver: rampa de 0 a
        # δ = 3·κ₀ < κ, así que toda la historia es descarga/recarga elástica.
        self._prescribir(nodes, 3.0 * self.KAPPA_0)
        NonlinearSolver(assembler, convergence=self._conv(), num_steps=3,
                        max_iter=8).solve(np.zeros(domain.total_dofs))

        state_after_unload = element.state.vars[0]
        # Irreversibilidad: κ y d no cambiaron durante la descarga
        self.assertAlmostEqual(state_after_unload['kappa'], kappa_loaded, places=12)
        self.assertAlmostEqual(state_after_unload['damage'], d_loaded, places=12)

    def test_damage_progresses_under_increasing_load(self):
        """En 4 pasos de desplazamiento creciente hasta ``6·κ₀`` el daño
        avanza paso a paso siguiendo ``d(λ·δ)``."""
        domain, element, nodes, assembler = self._modelo()
        delta = 6.0 * self.KAPPA_0
        self._prescribir(nodes, delta)

        damages = []

        def step_callback(step, U, load_factor):
            damages.append((load_factor, element.state.vars[0]['damage']))

        solver = NonlinearSolver(assembler, convergence=self._conv(), num_steps=4, max_iter=8)
        solver.solve(np.zeros(domain.total_dofs), step_callback=step_callback)

        self.assertEqual(len(damages), 4)
        for i, (lam, d) in enumerate(damages):
            self.assertAlmostEqual(d, self._d_analitico(lam * delta), places=10)
            if i:
                self.assertGreater(d, damages[i - 1][1])


if __name__ == '__main__':
    unittest.main()
