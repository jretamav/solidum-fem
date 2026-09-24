"""Contratos del módulo de usuario de discontinuidades (ADR 0020).

El programa principal prueba sus propias clases; cada módulo de usuario trae
su barrido. Aquí:

- la batería de contrato de elementos del programa principal
  (``test_element_contract_sweep``), reutilizada y restringida a
  ``CST_Embedded2D``;
- que el elemento no va por el camino por lotes (tiene estado propio);
- la ley cohesiva: tangente consistente frente a diferencias finitas,
  monotonía de κ y ω, validación de la rama lineal;
- el aislamiento del trial del salto frente a ``compute_global_stiffness``.

Las pruebas de la ley cohesiva y del trial vivían en tests mixtos del
programa principal hasta el ADR 0020.
"""
from __future__ import annotations

import copy
import os
import sys
import unittest

import numpy as np

_TESTS = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.abspath(os.path.join(_TESTS, '..')))
sys.path.insert(0, _TESTS)

import test_element_contract_sweep as barrido  # noqa: E402

from solidum.core.node import Node  # noqa: E402
from solidum.materials.elastic_2d import Elastic2D  # noqa: E402
from solidum.math.batch import element_is_batchable  # noqa: E402
from solidum.registry import ElementRegistry  # noqa: E402
from solidum.user.discontinuities import (  # noqa: E402
    CST_Embedded2D,
    CohesiveDamageIsotropic,
    CohesiveMaterialRegistry,
)

MAT_COHESIVO = CohesiveMaterialRegistry.create(
    'CohesiveDamageIsotropic', sigma_t0=3.0e6, G_f=100.0, K_e=1.0e13,
    softening='exponential')

FABRICA = (lambda c: c(1, barrido._nodos(barrido._TRI), barrido.MAT_2D,
                       cohesive_material=MAT_COHESIVO, thickness=barrido.ESPESOR), None)


# =========================================================================
# Batería de contrato del programa principal, restringida al módulo
# =========================================================================


class _SoloElModulo:
    """Hace que la batería del programa principal recorra sólo
    ``CST_Embedded2D`` durante cada prueba."""

    def setUp(self):
        super().setUp()
        self._original = barrido.nombres_programa_principal
        barrido.nombres_programa_principal = (
            lambda reg: ['CST_Embedded2D'] if reg is ElementRegistry else self._original(reg))
        barrido.FABRICAS['CST_Embedded2D'] = FABRICA

    def tearDown(self):
        barrido.nombres_programa_principal = self._original
        barrido.FABRICAS.pop('CST_Embedded2D', None)
        super().tearDown()


class TestDiscRegistro(_SoloElModulo, barrido.TestRegistroCompleto):

    def test_todo_elemento_registrado_tiene_fabrica(self):
        self.skipTest("guardián del registro del programa principal")

    def test_el_registro_no_esta_vacio(self):
        self.assertIn('CST_Embedded2D', ElementRegistry.names())


class TestDiscContratoDeclarativo(_SoloElModulo, barrido.TestContratoDeclarativo):
    pass


class TestDiscRigidezElemental(_SoloElModulo, barrido.TestRigidezElemental):
    pass


class TestDiscMasa(_SoloElModulo, barrido.TestMasa):
    pass


class TestDiscValidacionDefensiva(_SoloElModulo, barrido.TestValidacionDefensiva):
    pass


class TestDiscCaminoPorElemento(unittest.TestCase):

    def test_no_va_por_el_camino_por_lotes(self):
        """Estado propio fuera de ``ElementState`` (el salto): camino por elemento."""
        elem = FABRICA[0](CST_Embedded2D)
        self.assertFalse(element_is_batchable(elem))


# =========================================================================
# Ley cohesiva
# =========================================================================


def _tangent_fd(compute_fn, eps_base, state_committed, delta=1.0e-7):
    """Diferenciación finita central de ``compute_fn(x, state)[0]``."""
    n = len(eps_base)
    sigma_base, _, _ = compute_fn(eps_base, state_committed)
    C_fd = np.zeros((len(np.atleast_1d(sigma_base)), n))
    for j in range(n):
        e_j = np.zeros(n)
        e_j[j] = delta
        sigma_p, _, _ = compute_fn(eps_base + e_j, state_committed)
        sigma_m, _, _ = compute_fn(eps_base - e_j, state_committed)
        C_fd[:, j] = (np.asarray(sigma_p) - np.asarray(sigma_m)) / (2.0 * delta)
    return C_fd


def _rel_error(A_ref, A_test):
    return float(np.linalg.norm(A_ref - A_test) / max(np.linalg.norm(A_ref), 1.0e-14))


class TestCohesiveConsistentTangent(unittest.TestCase):
    """Tangente consistente del cohesivo Modo I (rango 1 sobre n⊗n)."""

    def _do_fd_test(self, softening: str):
        # K_e suficientemente grande para que el ablandamiento exponencial
        # cumpla ``K_e > σ_t0²/(2·G_F) = 2e10``.
        mat = CohesiveDamageIsotropic(sigma_t0=2.0e6, G_f=100.0, K_e=5.0e10,
                                      softening=softening)
        # Carga activa: apertura > κ_0 y muy por debajo de w_c.
        u_n = 1.5 * mat.kappa_0
        jump = np.array([u_n, 0.0])
        anchor = {'kappa': u_n * 0.9, 'damage': 0.0}
        _, T_alg, _ = mat.compute_traction(jump, anchor)
        T_fd = _tangent_fd(mat.compute_traction, jump, anchor, delta=1.0e-6 * u_n)
        err = _rel_error(T_alg, T_fd)
        self.assertLess(err, 1.0e-4, f"Cohesivo {softening}: error {err:.3e}")

    def test_linear_softening_loading(self):
        self._do_fd_test('linear')

    def test_exponential_softening_loading(self):
        self._do_fd_test('exponential')


class TestCohesiveMonotonicity(unittest.TestCase):
    """κ y ω no decrecen bajo apertura monótona."""

    def _run(self, softening: str):
        mat = CohesiveDamageIsotropic(sigma_t0=2.0e6, G_f=100.0, K_e=5.0e10,
                                      softening=softening)
        state, kappa_prev, d_prev, n_active = None, mat.kappa_0, 0.0, 0
        n_steps = 20
        for k in range(1, n_steps + 1):
            jump = np.array([5.0 * mat.kappa_0 * k / n_steps, 0.0])
            _, _, state_new = mat.compute_traction(jump, state)
            kappa_new, d_new = state_new['kappa'], state_new['damage']
            self.assertGreaterEqual(kappa_new, kappa_prev - 1.0e-14, f"{softening} paso {k}: κ decreció.")
            self.assertGreaterEqual(d_new, d_prev - 1.0e-14, f"{softening} paso {k}: ω decreció.")
            if d_new > d_prev + 1.0e-14:
                n_active += 1
            kappa_prev, d_prev, state = float(kappa_new), float(d_new), state_new
        self.assertGreater(n_active, 5, f"{softening}: sólo {n_active} pasos con daño activo.")

    def test_linear(self):
        self._run('linear')

    def test_exponential(self):
        self._run('exponential')


class TestCohesiveLinearValidation(unittest.TestCase):

    def test_linear_requires_wc_greater_than_kappa0(self):
        with self.assertRaises(ValueError):
            CohesiveDamageIsotropic(K_e=10.0, sigma_t0=3.0, G_f=0.1, softening="linear")
        CohesiveDamageIsotropic(K_e=1e6, sigma_t0=3.0, G_f=0.1, softening="linear")


# =========================================================================
# Aislamiento del trial
# =========================================================================


class TestTrialIsolation(unittest.TestCase):

    def test_embebido_conserva_el_trial_del_salto(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        nodes = [Node(i, list(coords[i])) for i in range(3)]
        for k, nd in enumerate(nodes):
            nd.add_dof('ux')
            nd.add_dof('uy')
            nd.dofs['ux'] = 2 * k
            nd.dofs['uy'] = 2 * k + 1
        bulk = Elastic2D(E=30.0e9, nu=0.2, hypothesis='plane_strain')
        coh = CohesiveDamageIsotropic(sigma_t0=2.5e6, G_f=100.0, K_e=1.0e13, softening='linear')
        emb = CST_Embedded2D(1, nodes, bulk, coh)
        # Tracción uniaxial en x por encima de σ_t0 → activa la discontinuidad.
        eps = 2.0 * 2.5e6 / 30.0e9
        U = np.array([0.0, 0.0, eps, 0.0, 0.0, 0.0])
        emb.compute_element_state(U)
        emb.commit_state()
        emb.prepare_step(U)
        self.assertIsNotNone(emb.discontinuity_state)
        emb.compute_element_state(1.5 * U)
        ds = emb.discontinuity_state
        jump_trial = np.copy(ds.jump_trial)
        coh_trial = copy.deepcopy(ds.cohesive_state_trial)
        stress_trial = np.copy(emb.state.stresses_trial[0])
        emb.compute_global_stiffness()
        np.testing.assert_array_equal(ds.jump_trial, jump_trial)
        self.assertEqual(ds.cohesive_state_trial, coh_trial)
        np.testing.assert_array_equal(emb.state.stresses_trial[0], stress_trial)


if __name__ == '__main__':
    unittest.main()
