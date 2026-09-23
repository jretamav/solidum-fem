"""La rigidez en ``u = 0`` no deja huella en el estado trial (deuda #18).

``Element.compute_global_stiffness`` y ``Assembler.assemble_system`` son
evaluaciones auxiliares (rigidez inicial para análisis lineal, modal o
dinámico). Antes sobrescribían el trial con el estado en ``u = 0``: si se
llamaban a mitad de un análisis no lineal, un ``commit`` posterior
consolidaba un estado falso. Ahora el trial previo se restaura al salir,
por los dos caminos (por elemento con ``ElementState`` clásico y por
lotes con vistas sobre la familia), en sólidos, en barras 1D y en el
elemento con discontinuidad embebida (que además guarda el trial del
salto).
"""
import copy
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401
from solidum import Assembler
from solidum.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic
from solidum.core.node import Node
from solidum.elements.solid_2d.embedded_cst import CST_Embedded2D
from solidum.elements.truss import Truss2D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.registry import MaterialRegistry

from test_batch_assembly import _domain_of, _field


def _trial_alpha(dom):
    return np.array([e.state.vars_trial[i]['alpha'] for e in dom.elements.values()
                     for i in range(e.state.num_ip)])


class TestRigidezInicialNoPisaElTrial(unittest.TestCase):

    def setUp(self):
        self.mat = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9)

    def _check(self, batch):
        dom = _domain_of('Quad4', self.mat)
        U = _field(dom, 1.0e-2)
        asm = Assembler(dom, batch=batch)
        asm.assemble_system()                               # rigidez del dominio virgen
        K_virgen = asm.K_global.toarray()
        asm.assemble_non_linear_system(U)
        alpha_trial = _trial_alpha(dom)
        self.assertGreater(alpha_trial.max(), 0.0)          # plastifica
        stress_trial = [np.array(e.state.stresses_trial[0]) for e in dom.elements.values()]

        asm.assemble_system()                               # rigidez en u = 0
        np.testing.assert_array_equal(_trial_alpha(dom), alpha_trial)
        for e, s in zip(dom.elements.values(), stress_trial):
            np.testing.assert_array_equal(np.array(e.state.stresses_trial[0]), s)
        # ...y K_0 es la rigidez elástica, la misma del dominio virgen.
        np.testing.assert_array_equal(asm.K_global.toarray(), K_virgen)
        # El commit consolida el trial del ensamblaje no lineal, no el de u = 0.
        asm.commit_all_states()
        committed = np.array([e.state.vars[i]['alpha'] for e in dom.elements.values()
                              for i in range(e.state.num_ip)])
        np.testing.assert_array_equal(committed, alpha_trial)

    def test_camino_por_elemento(self):
        self._check(batch=False)

    def test_camino_por_lotes(self):
        self._check(batch=True)

    def test_compute_global_stiffness_directo_en_solido(self):
        dom = _domain_of('Quad4', self.mat)
        U = _field(dom, 1.0e-2)
        elem = next(iter(dom.elements.values()))
        u_e = elem.get_local_displacements(U)
        elem.compute_element_state(u_e)
        trial = copy.deepcopy(elem.state.vars_trial)
        K0 = elem.compute_global_stiffness()
        self.assertEqual([v['alpha'] for v in elem.state.vars_trial], [v['alpha'] for v in trial])
        self.assertGreater(max(v['alpha'] for v in trial), 0.0)
        self.assertTrue(np.all(np.isfinite(K0)))

    def test_barra_1d_con_plasticidad(self):
        n1 = Node(1, [0.0, 0.0]); n2 = Node(2, [1.0, 0.0])
        for k, nd in enumerate((n1, n2)):
            nd.add_dof('ux'); nd.add_dof('uy')
            nd.dofs['ux'] = 2 * k; nd.dofs['uy'] = 2 * k + 1
        mat = Elastoplastic1D(E=2.0e11, sigma_y=2.0e8, H=1.0e9)
        truss = Truss2D(1, [n1, n2], mat, A=1.0e-3)
        K_e_plastic, _ = truss.compute_element_state(np.array([0.0, 0.0, 5.0e-3, 0.0]))
        alpha = truss.state.vars_trial[0]['alpha']
        self.assertGreater(alpha, 0.0)
        K0 = truss.compute_global_stiffness()
        self.assertEqual(truss.state.vars_trial[0]['alpha'], alpha)
        self.assertGreater(K0[0, 0], K_e_plastic[0, 0])    # elástica > tangente plástica

    def test_embebido_conserva_el_trial_del_salto(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        nodes = [Node(i, list(coords[i])) for i in range(3)]
        for k, nd in enumerate(nodes):
            nd.add_dof('ux'); nd.add_dof('uy')
            nd.dofs['ux'] = 2 * k; nd.dofs['uy'] = 2 * k + 1
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
