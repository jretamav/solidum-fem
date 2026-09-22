"""Ensamblaje por lotes (ADR 0014, fase 3): equivalencia con el camino por elemento.

El camino por elemento es la referencia física; el camino por lotes debe
reproducirlo a precisión de máquina en **toda** combinación registrada de
elemento × material batchable. El barrido recorre ``ElementRegistry`` y
``MaterialRegistry`` completos, así que una combinación nueva entra sin
que nadie la añada a mano.

Qué se verifica
---------------
1. ``K`` y ``F_int`` idénticos (``1e-14`` relativo) en geometría perturbada
   y con un campo de desplazamientos que plastifica o daña; el estado
   trial por punto de Gauss idéntico; y lo mismo tras un commit y una
   segunda evaluación (dependencia de la historia).
2. ``assemble_system`` (rigidez en ``u = 0``) idéntica.
3. Semántica trial/commit: un paso rechazado no contamina el committed;
   ``elem.state`` sigue siendo consultable como siempre y ``invalidate``
   devuelve ``ElementState`` clásicos sin perder historia.
4. Un solve no lineal completo (Newton con plasticidad) da el mismo ``U``
   y los mismos esfuerzos de post-proceso por ambos caminos.
5. Dominios mixtos (familias + elementos sin kernel), trozos de memoria
   mínimos, y elementos térmicos (que no llevan ``ElementState``).
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401
from solidum import Assembler, Domain
from solidum.core.element_state import ElementState
from solidum.core.node import Node
from solidum.math.batch import BatchedElementState, element_is_batchable
from solidum.math.solvers import NonlinearSolver
from solidum.math.convergence import ConvergenceCriterion
from solidum.registry import ElementRegistry, MaterialRegistry, ThermalMaterialRegistry
from solidum.utils.vtk_exporter import element_average_stress

from test_element_contract_sweep import (
    _HEX20, _HEX27, _HEX_VERTICES, _QUAD, _QUAD8, _QUAD9, _TET, _TET10, _TRI, _TRI6,
)

RNG = np.random.default_rng(2026)
ESPESOR = 0.01

# Geometría de referencia y factor de tamaño de cada elemento batchable.
_GEOM = {
    'Quad4': _QUAD, 'Tri3': _TRI, 'Quad8': _QUAD8, 'Quad9': _QUAD9, 'Tri6': _TRI6,
    'Hex8': _HEX_VERTICES, 'Hex20': _HEX20, 'Hex27': _HEX27, 'Tet4': _TET, 'Tet10': _TET10,
    'Quad4Thermal': _QUAD, 'Hex8Thermal': _HEX_VERTICES,
}

# Materiales con historia o sin ella para cada dimensión Voigt. Parámetros
# elegidos para que un campo de ~1e-2 de deformación entre en régimen
# inelástico (plastifique, dañe o alcance el cono).
def _materials_for(strain_dim):
    if strain_dim == 3:
        return {
            'Elastic2D': MaterialRegistry.create('Elastic2D', E=2.0e11, nu=0.3),
            'Elastic2D-plane_strain': MaterialRegistry.create(
                'Elastic2D', E=2.0e11, nu=0.3, hypothesis='plane_strain'),
            'Orthotropic2D': MaterialRegistry.create(
                'Orthotropic2D', E1=15.0e9, E2=0.8e9, G12=0.7e9, nu12=0.35, theta=30.0),
            'VonMises2D-plane_strain': MaterialRegistry.create(
                'VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9),
            'VonMises2D-plane_stress': MaterialRegistry.create(
                'VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9, hypothesis='plane_stress'),
            'DruckerPrager2D': MaterialRegistry.create(
                'DruckerPrager2D', E=2.0e10, nu=0.2, cohesion=1.0e6, phi_deg=30.0, psi_deg=15.0),
            'IsotropicDamage2D': MaterialRegistry.create(
                'IsotropicDamage2D', E=2.0e10, nu=0.2, kappa_0=1.0e-4, alpha=200.0),
        }
    return {
        'Elastic3D': MaterialRegistry.create('Elastic3D', E=2.0e11, nu=0.3),
        'VonMises3D': MaterialRegistry.create('VonMises3D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9),
        'DruckerPrager3D': MaterialRegistry.create(
            'DruckerPrager3D', E=2.0e10, nu=0.2, cohesion=1.0e6, phi_deg=30.0, psi_deg=15.0),
        'IsotropicDamage3D': MaterialRegistry.create(
            'IsotropicDamage3D', E=2.0e10, nu=0.2, kappa_0=1.0e-4, alpha=200.0),
    }


def _make(name, nodes, material):
    cls = ElementRegistry.get(name)
    if name in ('Quad4Thermal', 'Hex8Thermal'):
        if name == 'Quad4Thermal':
            return cls(1, nodes, material, thickness=ESPESOR)
        return cls(1, nodes, material)
    if cls.STRAIN_DIM == 3:
        return cls(1, nodes, material, thickness=ESPESOR)
    return cls(1, nodes, material)


def _domain_of(name, material, n_copies=3, perturb=0.04):
    """``n_copies`` elementos del mismo tipo con nodos propios, geometría
    perturbada y desplazados a lo largo de x."""
    geom = np.asarray(_GEOM[name], dtype=float)
    d = geom.shape[1]
    span = geom[:, 0].max() - geom[:, 0].min()
    dom = Domain()
    nid = 1
    for c in range(n_copies):
        nodes = []
        for p in geom:
            q = p + perturb * span * (RNG.random(d) - 0.5)
            q[0] += c * 1.5 * span
            nodes.append(dom.add_node(nid, list(q)))
            nid += 1
        elem = _make(name, nodes, material)
        elem.id = c + 1
        dom.add_element(elem)
    dom.generate_equation_numbers()
    return dom


def _field(dom, amplitude):
    U = np.zeros(dom.total_dofs)
    for node in dom.nodes.values():
        x = np.asarray(node.coordinates, float)
        for k, dof in enumerate(node.dofs.values()):
            U[dof] = amplitude * (0.5 * x[k % len(x)] ** 2 + 0.3 * x[0] * x[-1]
                                  + 0.2 * (RNG.random() - 0.5))
    return U


def _assemble_both(dom, U, **kw):
    """Ensambla por elemento y por lotes sobre el mismo dominio.

    El primero libera sus familias al invalidar, así que el segundo parte
    exactamente del mismo estado interno."""
    ref = Assembler(dom, batch=False)
    K0, F0 = ref.assemble_non_linear_system(U)
    lot = Assembler(dom, batch=True, **kw)
    K1, F1 = lot.assemble_non_linear_system(U)
    return (ref, K0, F0), (lot, K1, F1)


def _assert_same(K0, F0, K1, F1, rtol=1e-14):
    A0 = K0.toarray()
    A1 = K1.toarray()
    escala = np.abs(A0).max()
    np.testing.assert_allclose(A1, A0, rtol=0.0, atol=rtol * escala)
    escala_f = max(np.abs(F0).max(), escala * 1e-12)
    np.testing.assert_allclose(F1, F0, rtol=0.0, atol=rtol * escala_f)


def _trial_state(dom):
    """Estado trial de todo el dominio, apilado por variable: ``{clave:
    (n_elem·n_gp, ...)}`` más ``'__stress__'``."""
    out = {}
    for elem in dom.elements.values():
        st = elem.state
        if st is None:
            continue
        for i in range(st.num_ip):
            v = st.vars_trial[i]
            if v is not None:
                for k in v:
                    out.setdefault(k, []).append(np.atleast_1d(np.array(v[k], float)))
            out.setdefault('__stress__', []).append(np.array(st.stresses_trial[i], float))
    return {k: np.stack(vals) for k, vals in out.items()}


def _assert_same_state(s0, s1, rtol=1e-13):
    """Las dos rutas ejecutan las mismas funciones compiladas por punto de
    Gauss; lo único que difiere es el orden de suma de ``B·u`` (BLAS en el
    camino por elemento, bucle explícito en el kernel), que mueve la
    deformación en el último bit y con ella el estado. Se admite ``1e-13``
    relativo a la escala de cada variable en el dominio."""
    assert set(s0) == set(s1), (sorted(s0), sorted(s1))
    for k in s0:
        escala = max(np.abs(s0[k]).max(), 1e-300)
        np.testing.assert_allclose(s1[k], s0[k], rtol=0.0, atol=rtol * escala,
                                   err_msg=f"variable {k}")


class TestEquivalenciaBarrido(unittest.TestCase):
    """Todo elemento batchable × todo material batchable de su dimensión."""

    def test_barrido_registro_completo(self):
        cubiertos = 0
        for nombre in ElementRegistry.names():
            cls = ElementRegistry.get(nombre)
            if getattr(cls, 'BATCH_KINEMATICS', None) is None:
                continue
            if nombre in ('Quad4Thermal', 'Hex8Thermal'):
                continue    # térmicos: caso propio abajo
            for etiqueta, material in _materials_for(cls.STRAIN_DIM).items():
                with self.subTest(elemento=nombre, material=etiqueta):
                    self._check(nombre, material)
                    cubiertos += 1
        self.assertGreaterEqual(cubiertos, 10 * 4)

    def _check(self, nombre, material):
        dom = _domain_of(nombre, material)
        for elem in dom.elements.values():
            self.assertTrue(element_is_batchable(elem))
        U1 = _field(dom, 1.0e-2)
        (ref, K0, F0), (lot, K1, F1) = _assemble_both(dom, U1)
        self.assertEqual(len(lot.families), 1)
        _assert_same(K0, F0, K1, F1)
        # Estado trial idéntico (mismas funciones compiladas en ambos caminos).
        lot_state = _trial_state(dom)
        lot.invalidate()
        ref2 = Assembler(dom, batch=False)
        ref2.assemble_non_linear_system(U1)
        _assert_same_state(_trial_state(dom), lot_state)
        # Historia: commit y segunda evaluación con otro campo.
        ref2.commit_all_states()
        U2 = 0.4 * U1
        K0b, F0b = ref2.assemble_non_linear_system(U2)
        ref_state = _trial_state(dom)
        lot2 = Assembler(dom, batch=True)
        # El dominio ya está comiteado en U1: el lote lo adopta tal cual.
        K1b, F1b = lot2.assemble_non_linear_system(U2)
        _assert_same(K0b, F0b, K1b, F1b)
        _assert_same_state(_trial_state(dom), ref_state)

    def test_rigidez_lineal_assemble_system(self):
        for nombre, dim in (('Quad4', 3), ('Tri6', 3), ('Hex8', 6), ('Tet10', 6)):
            with self.subTest(elemento=nombre):
                material = _materials_for(dim)['Elastic2D' if dim == 3 else 'Elastic3D']
                dom = _domain_of(nombre, material)
                ref = Assembler(dom, batch=False)
                ref.assemble_system()
                lot = Assembler(dom, batch=True)
                lot.assemble_system()
                A0 = ref.K_global.toarray()
                A1 = lot.K_global.toarray()
                np.testing.assert_allclose(A1, A0, rtol=0.0, atol=1e-14 * np.abs(A0).max())

    def test_termicos(self):
        for nombre, dim in (('Quad4Thermal', 2), ('Hex8Thermal', 3)):
            with self.subTest(elemento=nombre):
                k = np.array([[50.0, 5.0], [5.0, 30.0]]) if dim == 2 else 40.0
                mat = ThermalMaterialRegistry.create('ThermalConduction', k=k, dim=dim)
                dom = _domain_of(nombre, mat)
                T = _field(dom, 10.0)
                (ref, K0, F0), (lot, K1, F1) = _assemble_both(dom, T)
                self.assertEqual(len(lot.families), 1)
                _assert_same(K0, F0, K1, F1, rtol=1e-13)
                for elem in dom.elements.values():
                    self.assertIsNone(elem.state)


class TestSemanticaTrialCommit(unittest.TestCase):

    def setUp(self):
        self.mat = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9)
        self.dom = _domain_of('Quad4', self.mat)
        self.asm = Assembler(self.dom, batch=True)
        self.U1 = _field(self.dom, 1.0e-2)

    def test_el_estado_es_una_vista_sobre_la_familia(self):
        self.asm.assemble_non_linear_system(self.U1)
        for elem in self.dom.elements.values():
            self.assertIsInstance(elem.state, BatchedElementState)
            self.assertIsInstance(elem.state, ElementState)
            self.assertEqual(elem.state.num_ip, 4)
            self.assertIn('alpha', elem.state.vars_trial[0])
            self.assertEqual(elem.state.vars[0]['alpha'], 0.0)   # aún sin commit

    def test_paso_rechazado_no_contamina_el_committed(self):
        self.asm.assemble_non_linear_system(self.U1)
        self.asm.commit_all_states()
        elem = next(iter(self.dom.elements.values()))
        alpha_1 = elem.state.vars[0]['alpha']
        self.assertGreater(alpha_1, 0.0)
        self.asm.assemble_non_linear_system(3.0 * self.U1)      # tentativa rechazada
        self.assertGreater(elem.state.vars_trial[0]['alpha'], alpha_1)
        self.assertEqual(elem.state.vars[0]['alpha'], alpha_1)
        # La siguiente evaluación parte del committed, como en ElementState.
        self.asm.assemble_non_linear_system(self.U1)
        self.assertEqual(elem.state.vars_trial[0]['alpha'], alpha_1)

    def test_commit_de_familia_deja_trial_igual_a_committed(self):
        self.asm.assemble_non_linear_system(self.U1)
        self.asm.commit_all_states()
        for elem in self.dom.elements.values():
            for i in range(elem.state.num_ip):
                self.assertEqual(elem.state.vars[i]['alpha'], elem.state.vars_trial[i]['alpha'])
                np.testing.assert_array_equal(elem.state.stresses[i], elem.state.stresses_trial[i])

    def test_invalidate_devuelve_element_state_clasico_con_la_historia(self):
        self.asm.assemble_non_linear_system(self.U1)
        self.asm.commit_all_states()
        elem = next(iter(self.dom.elements.values()))
        alpha_1 = elem.state.vars[0]['alpha']
        self.asm.invalidate()
        self.assertEqual(type(elem.state), ElementState)
        self.assertEqual(elem.state.vars[0]['alpha'], alpha_1)
        self.assertEqual(self.asm.families, [])
        # Un nuevo ensamblaje readopta el estado sin perderlo.
        self.asm.assemble_non_linear_system(self.U1)
        self.assertEqual(elem.state.vars[0]['alpha'], alpha_1)

    def test_commit_individual_del_elemento_sigue_funcionando(self):
        self.asm.assemble_non_linear_system(self.U1)
        elems = list(self.dom.elements.values())
        elems[0].commit_state()
        self.assertGreater(elems[0].state.vars[0]['alpha'], 0.0)
        self.assertEqual(elems[1].state.vars[0]['alpha'], 0.0)

    def test_compute_gauss_state_y_post_proceso_coinciden(self):
        self.asm.assemble_non_linear_system(self.U1)
        self.asm.commit_all_states()
        gs_lot = [e.compute_gauss_state(self.U1)['stress'].copy() for e in self.dom.elements.values()]
        avg_lot = [element_average_stress(e, self.U1) for e in self.dom.elements.values()]
        self.asm.invalidate()
        gs_ref = [e.compute_gauss_state(self.U1)['stress'] for e in self.dom.elements.values()]
        avg_ref = [element_average_stress(e, self.U1) for e in self.dom.elements.values()]
        for a, b in zip(gs_lot, gs_ref):
            np.testing.assert_array_equal(a, b)
        for a, b in zip(avg_lot, avg_ref):
            np.testing.assert_array_equal(a, b)

    def test_batch_false_no_crea_familias(self):
        asm = Assembler(self.dom, batch=False)
        asm.assemble_non_linear_system(self.U1)
        self.assertEqual(asm.families, [])
        for elem in self.dom.elements.values():
            self.assertEqual(type(elem.state), ElementState)


class TestDominiosMixtosYTrozos(unittest.TestCase):

    def _quad_mesh(self, nx, ny, material, quad_name='Quad4'):
        dom = Domain()
        nid = 1
        for j in range(ny + 1):
            for i in range(nx + 1):
                dom.add_node(nid, [i / nx + 0.02 * np.sin(3 * j), j / ny + 0.02 * np.cos(2 * i)])
                nid += 1
        eid = 1
        cls = ElementRegistry.get(quad_name)
        for j in range(ny):
            for i in range(nx):
                n0 = 1 + i + j * (nx + 1)
                nodes = [dom.nodes[n0], dom.nodes[n0 + 1], dom.nodes[n0 + nx + 2], dom.nodes[n0 + nx + 1]]
                dom.add_element(cls(eid, nodes, material, thickness=ESPESOR))
                eid += 1
        return dom, eid

    def test_malla_con_nodos_compartidos_y_dos_familias(self):
        """Dos materiales J2 distintos (dos familias) en una malla conectada:
        la acumulación de F_int y de K sobre DOFs compartidos coincide."""
        m1 = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9)
        m2 = MaterialRegistry.create('VonMises2D', E=1.0e11, nu=0.25, sigma_y=1.0e8, H=0.0)
        dom, _ = self._quad_mesh(4, 3, m1)
        for k, elem in enumerate(dom.elements.values()):
            if k % 2:
                elem.material = m2
        dom.generate_equation_numbers()
        U = _field(dom, 1.0e-2)
        (ref, K0, F0), (lot, K1, F1) = _assemble_both(dom, U)
        self.assertEqual(len(lot.families), 2)
        _assert_same(K0, F0, K1, F1)

    def test_dominio_mixto_con_elemento_sin_kernel(self):
        """Quad4 (familia) + Truss2D (camino por elemento) compartiendo nodos."""
        mat = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9)
        dom, eid = self._quad_mesh(3, 2, mat)
        truss = ElementRegistry.get('Truss2D')
        mat1d = MaterialRegistry.create('Elastoplastic1D', E=2.0e11, sigma_y=2.0e8)
        dom.add_element(truss(eid, [dom.nodes[1], dom.nodes[4]], mat1d, A=1.0e-3))
        dom.add_element(truss(eid + 1, [dom.nodes[5], dom.nodes[12]], mat1d, A=1.0e-3))
        dom.generate_equation_numbers()
        U = _field(dom, 1.0e-2)
        (ref, K0, F0), (lot, K1, F1) = _assemble_both(dom, U)
        self.assertEqual(len(lot.families), 1)
        self.assertEqual(len(lot._loose), 2)
        _assert_same(K0, F0, K1, F1)
        for elem in dom.elements.values():
            if type(elem).__name__ == 'Truss2D':
                self.assertEqual(type(elem.state), ElementState)

    def test_presupuesto_de_memoria_minimo_fuerza_trozos_de_un_elemento(self):
        mat = MaterialRegistry.create('IsotropicDamage2D', E=2.0e10, nu=0.2, kappa_0=1.0e-4, alpha=200.0)
        dom, _ = self._quad_mesh(4, 4, mat)
        dom.generate_equation_numbers()
        U = _field(dom, 1.0e-2)
        (ref, K0, F0), (lot, K1, F1) = _assemble_both(dom, U, batch_memory_budget=1)
        self.assertEqual(lot.families[0].chunk, 1)
        _assert_same(K0, F0, K1, F1)

    def test_elemento_anadido_tras_ensamblar_entra_en_la_familia(self):
        mat = MaterialRegistry.create('Elastic2D', E=2.0e11, nu=0.3)
        dom, eid = self._quad_mesh(2, 2, mat)
        dom.generate_equation_numbers()
        asm = Assembler(dom, batch=True)
        asm.assemble_non_linear_system(np.zeros(dom.total_dofs))
        self.assertEqual(asm.families[0].n_elements, 4)
        n = dom.add_node(100, [2.0, 0.0])
        n2 = dom.add_node(101, [2.0, 0.5])
        dom.add_element(ElementRegistry.get('Quad4')(
            eid, [dom.nodes[3], n, n2, dom.nodes[6]], mat, thickness=ESPESOR))
        dom.generate_equation_numbers()
        K1, _ = asm.assemble_non_linear_system(np.zeros(dom.total_dofs))
        self.assertEqual(asm.families[0].n_elements, 5)
        K0, _ = Assembler(dom, batch=False).assemble_non_linear_system(np.zeros(dom.total_dofs))
        np.testing.assert_allclose(K1.toarray(), K0.toarray(), rtol=0.0, atol=1e-14 * np.abs(K0).max())

    def test_solve_no_lineal_completo_coincide(self):
        """Newton con plasticidad J2 por ambos caminos: mismo U, mismos σ."""
        def run(batch):
            mat = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=2.0e9)
            dom, _ = self._quad_mesh(4, 2, mat)
            for node in dom.nodes.values():
                if abs(node.coordinates[0]) < 1e-9 + 0.02:
                    node.fix_dof('ux'); node.fix_dof('uy')
            dom.generate_equation_numbers()
            F = np.zeros(dom.total_dofs)
            for node in dom.nodes.values():
                if node.coordinates[0] > 0.97:
                    F[node.dofs['uy']] = -4.0e7 * ESPESOR
            asm = Assembler(dom, batch=batch)
            solver = NonlinearSolver(asm, convergence=ConvergenceCriterion(rtol_force=1e-10, rtol_disp=1e-10),
                                     num_steps=4)
            U = solver.solve(F)
            asm.commit_all_states()
            sig = np.concatenate([e.compute_gauss_state(U)['stress'].ravel() for e in dom.elements.values()])
            alpha = np.array([e.state.vars[i]['alpha'] for e in dom.elements.values() for i in range(4)])
            return U, sig, alpha
        U0, s0, a0 = run(False)
        U1, s1, a1 = run(True)
        np.testing.assert_allclose(U1, U0, rtol=1e-11, atol=1e-11 * np.abs(U0).max())
        np.testing.assert_allclose(s1, s0, rtol=1e-10, atol=1e-10 * np.abs(s0).max())
        self.assertGreater(a0.max(), 0.0)
        np.testing.assert_allclose(a1, a0, rtol=1e-10, atol=1e-10 * a0.max())


if __name__ == '__main__':
    unittest.main()
