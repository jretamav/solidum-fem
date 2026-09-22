"""Post-proceso por lotes (ADR 0014 §9): equivalencia con ``compute_gauss_state``.

``Family.gauss_state(U)`` evalúa de una vez lo que ``compute_gauss_state``
devuelve elemento a elemento; el exportador VTK usa la familia cuando
existe y el camino por elemento cuando no. Lo que se verifica:

1. Barrido sobre los registros: para todo elemento batchable × materiales
   con y sin historia, tras un paso inelástico comiteado, la rodaja ``[i]``
   de ``gauss_state`` coincide con ``elements[i].compute_gauss_state(U)``
   en deformación, esfuerzo (``1e-13``: orden de suma de ``B·u``) y en la
   posición global de los puntos de Gauss (``1e-12``); igual en térmicos
   (gradiente y flujo con ``q = −k·∇T``). El trial no se toca.
2. ``gauss_states(domain, U)`` cubre familias y elementos sueltos.
3. El archivo VTK exportado con las familias vigentes es idéntico al
   exportado tras ``invalidate()`` (camino por elemento): esfuerzos,
   ``σ_zz`` de plane strain con ``ε^p_zz``, ``Internal_State``, flujo,
   suavizado nodal y Von Mises; también en dominios mixtos con barras.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401
from solidum import Assembler, Domain
from solidum.math.batch import family_of, gauss_states
from solidum.registry import ElementRegistry, MaterialRegistry, ThermalMaterialRegistry
from solidum.utils.vtk_exporter import VtkExporter

from test_batch_assembly import ESPESOR, _domain_of, _field, _materials_for

try:
    import meshio
except ImportError:      # pragma: no cover
    meshio = None


def _quad_mesh(nx, ny, material, elem_name='Quad4'):
    dom = Domain()
    nid = 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            dom.add_node(nid, [i / nx + 0.02 * np.sin(3 * j), j / ny + 0.02 * np.cos(2 * i)])
            nid += 1
    eid = 1
    cls = ElementRegistry.get(elem_name)
    for j in range(ny):
        for i in range(nx):
            n0 = 1 + i + j * (nx + 1)
            nodes = [dom.nodes[n0], dom.nodes[n0 + 1], dom.nodes[n0 + nx + 2], dom.nodes[n0 + nx + 1]]
            dom.add_element(cls(eid, nodes, material, thickness=ESPESOR))
            eid += 1
    return dom, eid


class TestGaussStatePorFamilia(unittest.TestCase):

    def _check_mechanical(self, nombre, material):
        dom = _domain_of(nombre, material)
        U1 = _field(dom, 1.0e-2)
        asm = Assembler(dom, batch=True)
        asm.assemble_non_linear_system(U1)
        asm.commit_all_states()               # historia inelástica comiteada
        U2 = 0.6 * U1
        asm.assemble_non_linear_system(1.7 * U1)   # trial distinto del committed
        fam = asm.families[0]
        trial_before = fam.state.S_trial.copy()
        gs = fam.gauss_state(U2)
        np.testing.assert_array_equal(fam.state.S_trial, trial_before)   # trial intacto
        self.assertEqual(gs['strain'].shape, (3, fam.n_gp, fam.n_sigma))
        self.assertIsNotNone(gs['points_global'])
        for i, elem in enumerate(fam.elements):
            ref = elem.compute_gauss_state(U2)
            for key in ('strain', 'stress'):
                escala = max(np.abs(ref[key]).max(), 1e-300)
                np.testing.assert_allclose(gs[key][i], ref[key], rtol=0.0, atol=1e-13 * escala,
                                           err_msg=f"{nombre} {key}")
            np.testing.assert_allclose(gs['points_natural'], ref['points_natural'], atol=1e-14)
            np.testing.assert_allclose(gs['points_global'][i], ref['points_global'], rtol=1e-12, atol=1e-12)
        # gauss_states(domain, U): mismas rodajas.
        todos = gauss_states(dom, U2)
        for i, elem in enumerate(fam.elements):
            np.testing.assert_array_equal(todos[elem.id]['stress'], gs['stress'][i])
        asm.invalidate()
        self.assertIsNone(family_of(fam.elements[0]))

    def test_barrido_mecanico(self):
        cubiertos = 0
        for nombre in ElementRegistry.names():
            cls = ElementRegistry.get(nombre)
            if getattr(cls, 'BATCH_KINEMATICS', None) is None or nombre.endswith('Thermal'):
                continue
            materiales = _materials_for(cls.STRAIN_DIM)
            for etiqueta, material in materiales.items():
                with self.subTest(elemento=nombre, material=etiqueta):
                    self._check_mechanical(nombre, material)
                    cubiertos += 1
        self.assertGreaterEqual(cubiertos, 10 * 4)

    def test_termicos(self):
        for nombre, dim in (('Quad4Thermal', 2), ('Hex8Thermal', 3)):
            with self.subTest(elemento=nombre):
                k = np.array([[50.0, 5.0], [5.0, 30.0]]) if dim == 2 else 40.0
                mat = ThermalMaterialRegistry.create('ThermalConduction', k=k, dim=dim)
                dom = _domain_of(nombre, mat)
                T = _field(dom, 10.0)
                asm = Assembler(dom, batch=True)
                asm.assemble_non_linear_system(T)
                fam = asm.families[0]
                gs = fam.gauss_state(T)
                self.assertIn('flux', gs)
                self.assertNotIn('stress', gs)
                for i, elem in enumerate(fam.elements):
                    ref = elem.compute_gauss_state(T)
                    for key in ('grad_T', 'flux'):
                        escala = max(np.abs(ref[key]).max(), 1e-300)
                        np.testing.assert_allclose(gs[key][i], ref[key], rtol=0.0, atol=1e-13 * escala)
                    np.testing.assert_allclose(gs['points_global'][i], ref['points_global'], rtol=1e-12, atol=1e-12)
                asm.invalidate()

    def test_serie_y_paralelo_bit_a_bit(self):
        mat = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9)
        dom = _domain_of('Quad4', mat, n_copies=70)
        U = _field(dom, 1.0e-2)
        out = {}
        for parallel in (False, True):
            asm = Assembler(dom, batch=True, parallel=parallel)
            asm.assemble_non_linear_system(U)
            asm.commit_all_states()
            out[parallel] = asm.families[0].gauss_state(0.5 * U)
            asm.invalidate()
        for key in ('strain', 'stress', 'points_global'):
            np.testing.assert_array_equal(out[False][key], out[True][key])

    def test_jacobiano_degenerado_lanza_con_el_id(self):
        mat = MaterialRegistry.create('Elastic2D', E=2.0e11, nu=0.3)
        dom, _ = _quad_mesh(3, 2, mat)
        dom.generate_equation_numbers()
        asm = Assembler(dom, batch=True)
        asm.assemble_system()
        malo = dom.elements[4]
        malo.nodes[1], malo.nodes[3] = malo.nodes[3], malo.nodes[1]
        fam = asm.families[0]
        fam.X[3] = malo.batch_reference_coordinates(2)      # geometría ya adoptada
        with self.assertRaises(ValueError) as ctx:
            fam.gauss_state(np.zeros(dom.total_dofs))
        self.assertIn('id=4', str(ctx.exception))


@unittest.skipIf(meshio is None, "meshio no instalado")
class TestExportadorPorLotes(unittest.TestCase):

    def setUp(self):
        self.tmp = []

    def tearDown(self):
        for p in self.tmp:
            try:
                os.remove(p)
            except OSError:
                pass

    def _path(self, name):
        p = os.path.join(os.path.dirname(__file__), name)
        self.tmp.append(p)
        return p

    def _export_both(self, dom, U, tag):
        """Exporta con las familias vigentes y tras invalidar (por elemento)."""
        asm = Assembler(dom, batch=True)
        asm.assemble_non_linear_system(U)
        asm.commit_all_states()
        asm.assemble_non_linear_system(1.5 * U)       # trial ≠ committed
        self.assertTrue(any(family_of(e) is not None for e in dom.elements.values()))
        p_lot = self._path(f"_test_vtk_lot_{tag}.vtu")
        VtkExporter(dom).export(p_lot, U=0.8 * U)
        asm.invalidate()
        self.assertTrue(all(family_of(e) is None for e in dom.elements.values()))
        p_ref = self._path(f"_test_vtk_ref_{tag}.vtu")
        VtkExporter(dom).export(p_ref, U=0.8 * U)
        return meshio.read(p_lot), meshio.read(p_ref)

    def _assert_same_mesh(self, lot, ref):
        self.assertEqual([b.type for b in lot.cells], [b.type for b in ref.cells])
        for a, b in zip(lot.cells, ref.cells):
            np.testing.assert_array_equal(a.data, b.data)
        self.assertEqual(set(lot.point_data), set(ref.point_data))
        self.assertEqual(set(lot.cell_data), set(ref.cell_data))
        for name in ref.point_data:
            escala = max(np.abs(ref.point_data[name]).max(), 1e-300)
            np.testing.assert_allclose(lot.point_data[name], ref.point_data[name],
                                       rtol=0.0, atol=1e-12 * escala, err_msg=name)
        for name in ref.cell_data:
            for a, b in zip(lot.cell_data[name], ref.cell_data[name]):
                escala = max(np.abs(b).max(), 1e-300)
                np.testing.assert_allclose(a, b, rtol=0.0, atol=1e-12 * escala, err_msg=name)

    def test_plane_strain_j2_sigma_zz_e_internal_state(self):
        mat = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9)
        dom, _ = _quad_mesh(5, 4, mat)
        dom.generate_equation_numbers()
        U = _field(dom, 1.0e-2)
        lot, ref = self._export_both(dom, U, 'j2')
        self._assert_same_mesh(lot, ref)
        self.assertGreater(np.abs(ref.cell_data['Sigma_ZZ'][0]).max(), 0.0)
        self.assertGreater(ref.cell_data['Internal_State'][0].max(), 0.0)   # α > 0

    def test_dano_y_tri6(self):
        mat = MaterialRegistry.create('IsotropicDamage2D', E=2.0e10, nu=0.2, kappa_0=1.0e-4, alpha=200.0,
                                      hypothesis='plane_strain')
        dom = _domain_of('Tri6', mat, n_copies=6)
        U = _field(dom, 1.0e-2)
        lot, ref = self._export_both(dom, U, 'dano')
        self._assert_same_mesh(lot, ref)
        self.assertGreater(ref.cell_data['Internal_State'][0].max(), 0.0)   # d > 0

    def test_hex8_j2_3d(self):
        mat = MaterialRegistry.create('VonMises3D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9)
        dom = _domain_of('Hex8', mat, n_copies=5)
        U = _field(dom, 1.0e-2)
        lot, ref = self._export_both(dom, U, 'hex8')
        self._assert_same_mesh(lot, ref)
        self.assertGreater(np.abs(ref.cell_data['Tau_YZ'][0]).max(), 0.0)

    def test_termico(self):
        mat = ThermalMaterialRegistry.create('ThermalConduction', k=np.array([[50.0, 5.0], [5.0, 30.0]]), dim=2)
        dom = _domain_of('Quad4Thermal', mat, n_copies=4)
        T = _field(dom, 10.0)
        lot, ref = self._export_both(dom, T, 'termico')
        self._assert_same_mesh(lot, ref)
        self.assertGreater(np.abs(ref.cell_data['Flux'][0]).max(), 0.0)

    def test_dominio_mixto_con_barras_y_dos_familias(self):
        m1 = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9)
        m2 = MaterialRegistry.create('Elastic2D', E=1.0e11, nu=0.25, hypothesis='plane_strain')
        dom, eid = _quad_mesh(4, 3, m1)
        for k, elem in enumerate(dom.elements.values()):
            if k % 3 == 0:
                elem.material = m2
        truss = ElementRegistry.get('Truss2D')
        mat1d = MaterialRegistry.create('Elastoplastic1D', E=2.0e11, sigma_y=2.0e8)
        dom.add_element(truss(eid, [dom.nodes[1], dom.nodes[7]], mat1d, A=1.0e-3))
        dom.generate_equation_numbers()
        U = _field(dom, 1.0e-2)
        lot, ref = self._export_both(dom, U, 'mixto')
        self._assert_same_mesh(lot, ref)
        self.assertEqual({b.type for b in ref.cells}, {'quad', 'line'})

    def test_sin_u_exporta_estado_committed_por_familia(self):
        mat = MaterialRegistry.create('VonMises2D', E=2.0e11, nu=0.3, sigma_y=2.5e8, H=1.0e9)
        dom, _ = _quad_mesh(3, 3, mat)
        dom.generate_equation_numbers()
        U = _field(dom, 1.0e-2)
        asm = Assembler(dom, batch=True)
        asm.assemble_non_linear_system(U)
        asm.commit_all_states()
        p = self._path("_test_vtk_lot_sinU.vtu")
        VtkExporter(dom).export(p, U=None)
        mesh = meshio.read(p)
        self.assertEqual(float(np.abs(mesh.cell_data['Sigma_XX'][0]).max()), 0.0)
        alpha = np.array([np.mean([e.state.vars[i]['alpha'] for i in range(4)])
                          for e in dom.elements.values()])
        np.testing.assert_allclose(mesh.cell_data['Internal_State'][0], alpha, rtol=1e-12)


if __name__ == '__main__':
    unittest.main()
