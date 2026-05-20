"""Tests de integración end-to-end: ``CST_Embedded2D`` + solvers no lineales.

Bring-up de la fase 2 del ADR 0010 a través del pipeline completo
Domain → Assembler → solver no lineal con el hook ``prepare_step`` activo,
condensación estática elemental y activación real durante el paso.

Estos tests no son benchmarks experimentales (eso es la fase 4 del ADR
0010, Van Vliet). Su propósito es validar que el cableado del subsistema
sobrevive un paso real del solver con elementos agrietándose: que el hook
``prepare_step`` se invoca, que la activación toma efecto, que el Newton
global converge con K_cond en lugar de K_e, y que la irreversibilidad se
mantiene tras una descarga.

Setup: cuadrado unitario mallado en dos ``CST_Embedded2D`` (patch test
de tracción uniaxial). Bulk elástico (``Elastic2D`` plane strain).
Cohesivo ``CohesiveDamageIsotropic`` con softening exponencial (asintótico
a saturación, sin transición brusca como el lineal).
"""
from __future__ import annotations

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic
from solidum.core.domain import Domain
from solidum.elements.solid_2d.embedded_cst import CST_Embedded2D
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import NonlinearSolver


# --- Parámetros del problema ----------------------------------------------
E_BULK = 30.0e9
NU_BULK = 0.2
SIGMA_T0 = 2.0e6
G_F = 200.0
K_E = 1.0e13


def _build_two_cst_embedded_square():
    """Cuadrado unitario mallado en dos ``CST_Embedded2D`` con BCs de
    tracción uniaxial en x.

    Nodos::

        4 ----- 3
        |     / |
        |   /   |
        | /     |
        1 ----- 2

    Triángulos antihorarios: T1 = (1, 2, 3), T2 = (1, 3, 4).

    BCs:
        n1 = (0, 0) fix ux = uy = 0
        n2 = (1, 0) fix uy = 0
        n3 = (1, 1) libre
        n4 = (0, 1) fix ux = 0
    """
    domain = Domain()
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [1.0, 0.0])
    n3 = domain.add_node(3, [1.0, 1.0])
    n4 = domain.add_node(4, [0.0, 1.0])

    bulk = Elastic2D(E=E_BULK, nu=NU_BULK, hypothesis='plane_strain')
    cohesive = CohesiveDamageIsotropic(
        sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_E, softening='exponential',
    )

    t1 = CST_Embedded2D(1, [n1, n2, n3], bulk, cohesive)
    t2 = CST_Embedded2D(2, [n1, n3, n4], bulk, cohesive)
    domain.add_element(t1)
    domain.add_element(t2)

    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
    n2.fix_dof('uy', 0.0)
    n4.fix_dof('ux', 0.0)

    domain.generate_equation_numbers(verbose=False)

    return domain, (t1, t2), {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4}


def _apply_traction(domain, nodes, F_total):
    """Reparte ``F_total`` por igual entre nodos 2 y 3 (lado derecho)."""
    F = np.zeros(domain.total_dofs)
    F[nodes['n2'].dofs['ux']] = 0.5 * F_total
    F[nodes['n3'].dofs['ux']] = 0.5 * F_total
    return F


# =========================================================================
# Test 1 — régimen elástico, ningún elemento se activa
# =========================================================================


class TestElasticBelowThreshold(unittest.TestCase):
    """F_total por debajo del umbral de activación: respuesta lineal pura,
    ningún elemento se activa, ``prepare_step`` se invoca pero no produce
    efecto."""

    def test_no_activation(self):
        domain, elements, nodes = _build_two_cst_embedded_square()
        assembler = Assembler(domain)

        # σ_target = 0.5·σ_t0 ⇒ F_total = σ_target · área_sección (1·1)
        F_total = 0.5 * SIGMA_T0
        F_ext = _apply_traction(domain, nodes, F_total)

        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        solver = NonlinearSolver(assembler, convergence=conv, num_steps=2)
        solver.solve(F_ext)

        # Ningún elemento debió activarse.
        for elem in elements:
            self.assertIsNone(elem.discontinuity_state)


# =========================================================================
# Test 2 — el hook prepare_step se invoca correctamente por el solver
# =========================================================================


class TestPrepareStepIsCalled(unittest.TestCase):
    """Verifica el cableado del hook: el ``Assembler.prepare_all_steps``
    invoca ``Element.prepare_step`` en cada elemento del dominio, una vez
    por paso del solver, antes del primer ensamblaje del Newton.

    No depende de que el estado físico cruce el umbral (eso requiere
    afinación de pasos y manejo de softening que pertenece al benchmark
    de la fase 4, no a un test de bring-up).
    """

    def test_prepare_step_called_per_step_per_element(self):
        domain, elements, nodes = _build_two_cst_embedded_square()

        # Spy: cuenta invocaciones y registra el id del elemento por llamada.
        calls = []
        original = CST_Embedded2D.prepare_step

        def spy(self, U_committed):
            calls.append(self.id)
            original(self, U_committed)

        CST_Embedded2D.prepare_step = spy
        try:
            assembler = Assembler(domain)
            F_ext = _apply_traction(domain, nodes, 0.3 * SIGMA_T0)
            num_steps = 3
            NonlinearSolver(
                assembler,
                convergence=ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8),
                num_steps=num_steps,
            ).solve(F_ext)
        finally:
            CST_Embedded2D.prepare_step = original

        # Esperado: 2 elementos × num_steps invocaciones (sin bisección,
        # converge en cada paso de un golpe).
        self.assertEqual(len(calls), 2 * num_steps)
        # Cada elemento se invoca el mismo número de veces.
        for elem in elements:
            self.assertEqual(calls.count(elem.id), num_steps)


# =========================================================================
# Test 3 — solver lineal funciona con elementos pre-agrietados
# =========================================================================


class TestSolverWithPreCrackedElements(unittest.TestCase):
    """Pre-activar manualmente los elementos (saltando la lógica Rankine)
    y verificar que ``NonlinearSolver`` converge en régimen post-activación
    con cargas que no provocan softening severo.

    Este test cubre directamente el cableado **condensación + Newton
    global**: tras la activación, ``compute_element_state`` devuelve K_cond
    en lugar de K_e estándar; el ensamblador no debe darse cuenta. Si la
    integración está bien, el Newton global converge igual de bien que con
    elementos intactos en este régimen de carga moderada.
    """

    def test_loaded_pre_cracked_converges(self):
        domain, elements, nodes = _build_two_cst_embedded_square()

        # Pre-activar ambos elementos con n = (1, 0) (consistente con
        # tracción uniaxial en x).
        for elem in elements:
            coords = elem.get_coordinate_matrix(ndim=2)
            elem._activate(np.array([1.0, 0.0]), coords)
            self.assertIsNotNone(elem.discontinuity_state)

        assembler = Assembler(domain)

        # Carga moderada (subcrítica respecto a σ_t0): tras la activación
        # κ se inicia en κ_0 con ω = 0; la grieta abre poco, el cohesivo
        # está en rama elástica del salto. NonlinearSolver debe converger.
        F_total = 0.3 * SIGMA_T0
        F_ext = _apply_traction(domain, nodes, F_total)

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = NonlinearSolver(
            assembler, convergence=conv, num_steps=5, max_iter=20,
        )
        solver.solve(F_ext)

        # Los elementos siguen agrietados (irreversibilidad).
        for elem in elements:
            self.assertIsNotNone(elem.discontinuity_state)
            # El salto trial debe quedar consistente con el jump_committed
            # tras commit_state.
            ds = elem.discontinuity_state
            np.testing.assert_allclose(ds.jump_trial, ds.jump_committed, atol=1e-12)


# =========================================================================
# Test 4 — descarga sobre estado dañado: irreversibilidad
# =========================================================================


class TestUnloadingPreservesDamage(unittest.TestCase):
    """Pre-activar y cargar dañando; después descargar a F = 0 y verificar
    que ``ω`` y ``κ`` committed no decrecen (Kuhn-Tucker)."""

    def test_unloading_preserves_state(self):
        domain, elements, nodes = _build_two_cst_embedded_square()
        for elem in elements:
            coords = elem.get_coordinate_matrix(ndim=2)
            elem._activate(np.array([1.0, 0.0]), coords)

        assembler = Assembler(domain)
        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)

        # Fase 1: cargar para producir algo de daño cohesivo.
        F_ext_load = _apply_traction(domain, nodes, 0.5 * SIGMA_T0)
        NonlinearSolver(
            assembler, convergence=conv, num_steps=5, max_iter=20,
        ).solve(F_ext_load)

        damages = {
            e.id: e.discontinuity_state.cohesive_state_committed.get('damage', 0.0)
            for e in elements
        }
        kappas = {
            e.id: e.discontinuity_state.cohesive_state_committed.get('kappa', 0.0)
            for e in elements
        }

        # Fase 2: descargar a F = 0 (el solver arranca con U_current = 0
        # por defecto; aplicamos F_ext_unload = 0 y verificamos que la
        # corrida es trivial — ya estamos en equilibrio en U = 0 dado que
        # no hay carga, pero el daño committed se preserva).
        # Reusar el assembler y el dominio: el estado committed persiste.
        F_ext_unload = np.zeros(domain.total_dofs)
        NonlinearSolver(
            assembler, convergence=conv, num_steps=1, max_iter=20,
        ).solve(F_ext_unload)

        # Daño y κ committed: no decrecen tras la descarga.
        for elem in elements:
            d_now = elem.discontinuity_state.cohesive_state_committed.get('damage', 0.0)
            k_now = elem.discontinuity_state.cohesive_state_committed.get('kappa', 0.0)
            self.assertGreaterEqual(d_now, damages[elem.id] - 1e-15)
            self.assertGreaterEqual(k_now, kappas[elem.id] - 1e-15)
            self.assertIsNotNone(elem.discontinuity_state)


if __name__ == '__main__':
    unittest.main()
