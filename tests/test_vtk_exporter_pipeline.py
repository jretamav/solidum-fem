"""Exportador VTK como función pura de ``(dominio, U)``.

Complementa ``test_vtk_exporter.py`` (que monta ``U`` a mano) con la vía
oficial ``solidum.run`` y con las familias que el exportador antiguo omitía
(sólidos 3D, térmicos) o calculaba mal (Von Mises en plane strain).
"""
import os
import unittest

import meshio
import numpy as np

import solidum
from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4
from solidum.elements.solid_3d import Hex8, Tet4
from solidum.elements.thermal import Quad4Thermal
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.utils.vtk_exporter import VtkExporter, cell_type_for


def _plate(material, thickness=0.5):
    dom = Domain()
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        dom.add_node(i, [float(x), float(y)])
    dom.add_element(Quad4(1, [dom.nodes[i] for i in (1, 2, 3, 4)], material,
                          thickness=thickness))
    dom.generate_equation_numbers()
    for k in (1, 4):
        dom.nodes[k].fix_dof('ux'); dom.nodes[k].fix_dof('uy')
    F = np.zeros(dom.total_dofs)
    F[dom.nodes[2].dofs['ux']] = 10.0
    F[dom.nodes[3].dofs['ux']] = 10.0
    return dom, F


class TestVtkExporterPipeline(unittest.TestCase):

    def setUp(self):
        self.tmp = []

    def tearDown(self):
        for p in self.tmp:
            try:
                os.remove(p)
            except OSError:
                pass

    def _tmp_path(self, name):
        path = os.path.join(os.path.dirname(__file__), name)
        self.tmp.append(path)
        return path

    def test_linear_solver_run_exporta_esfuerzos_no_nulos(self):
        """Regresión (auditoría 2026-09-22): tras ``solidum.run`` con
        ``LinearSolver`` el VTU exportaba σ = 0 porque el solver lineal no
        comiteaba y el exportador leía el estado committed."""
        dom, F = _plate(Elastic2D(E=1000.0, nu=0.3))
        res = solidum.run(dom, F_applied=F)
        path = self._tmp_path("_test_vtk_linear_run.vtu")
        VtkExporter(dom).export(path, U=res.U, F_ext=F)
        mesh = meshio.read(path)
        sxx = mesh.cell_data["Sigma_XX"][0][0]
        gs = dom.elements[1].compute_gauss_state(res.U)
        self.assertAlmostEqual(sxx, gs['stress'][:, 0].mean(), places=10)
        self.assertGreater(abs(sxx), 1.0)
        # ``run`` deja el estado committed coherente con U.
        committed = np.mean([s[0] for s in dom.elements[1].state.stresses])
        self.assertAlmostEqual(committed, sxx, places=10)

    def test_export_no_depende_del_commit(self):
        """Misma U, con o sin commit previo → mismo archivo de esfuerzos."""
        dom, F = _plate(Elastic2D(E=1000.0, nu=0.3))
        from solidum.math.assembly import Assembler
        from solidum.math.solvers import LinearSolver
        U = LinearSolver(Assembler(dom)).solve(F)       # sin commit
        p1 = self._tmp_path("_test_vtk_nocommit.vtu")
        VtkExporter(dom).export(p1, U=U)
        s1 = meshio.read(p1).cell_data["Sigma_XX"][0][0]
        gs = dom.elements[1].compute_gauss_state(U)
        self.assertAlmostEqual(s1, gs['stress'][:, 0].mean(), places=10)
        self.assertGreater(abs(s1), 1.0)

    def test_von_mises_plane_strain_incluye_sigma_zz(self):
        """Extensión uniaxial en plane strain: σ_zz = ν·(σ_xx + σ_yy) y el
        Von Mises es el invariante 3D, no la fórmula de plane stress."""
        nu = 0.3
        dom, _ = _plate(Elastic2D(E=1000.0, nu=nu, hypothesis='plane_strain'))
        eps = 1e-3
        U = np.zeros(dom.total_dofs)
        for nd in dom.nodes.values():
            U[nd.dofs['ux']] = eps * nd.coordinates[0]
        path = self._tmp_path("_test_vtk_vm_plane_strain.vtu")
        VtkExporter(dom).export(path, U=U)
        mesh = meshio.read(path)
        sxx = mesh.cell_data["Sigma_XX"][0][0]
        syy = mesh.cell_data["Sigma_YY"][0][0]
        szz = mesh.cell_data["Sigma_ZZ"][0][0]
        vm = mesh.cell_data["Von_Mises"][0][0]
        self.assertAlmostEqual(szz, nu * (sxx + syy), places=9)
        expected = np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2))
        self.assertAlmostEqual(vm, expected, places=9)
        self.assertGreater(abs(vm - np.sqrt(sxx ** 2 + syy ** 2 - sxx * syy)), 1e-3)
        # El suavizado nodal arrastra las seis componentes.
        self.assertIn("Sigma_ZZ_nodal", mesh.point_data)
        self.assertAlmostEqual(float(mesh.point_data["Sigma_ZZ_nodal"][0]), szz, places=9)

    def test_solidos_3d_y_termicos(self):
        """Hex8/Tet4 → hexahedron/tetra con seis componentes; Quad4Thermal →
        quad con campo nodal Temperature y Flux por celda."""
        dom = Domain()
        hexc = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
        for i, c in enumerate(hexc, 1):
            dom.add_node(i, [float(v) for v in c])
        m3 = Elastic3D(E=1000.0, nu=0.3)
        dom.add_element(Hex8(1, [dom.nodes[i] for i in range(1, 9)], m3))
        for i, c in enumerate([(2, 0, 0), (3, 0, 0), (2, 1, 0), (2, 0, 1)], 9):
            dom.add_node(i, [float(v) for v in c])
        dom.add_element(Tet4(2, [dom.nodes[i] for i in range(9, 13)], m3))
        for i, c in enumerate([(5, 0), (6, 0), (6, 1), (5, 1)], 13):
            dom.add_node(i, [float(v) for v in c])
        mt = ThermalConduction(k=2.0, c=1.0, density=1.0, dim=2)
        dom.add_element(Quad4Thermal(3, [dom.nodes[i] for i in range(13, 17)], mt))
        dom.generate_equation_numbers()

        U = np.zeros(dom.total_dofs)
        eps = 1e-3
        for nd in dom.nodes.values():
            if 'ux' in nd.dofs:
                U[nd.dofs['ux']] = eps * nd.coordinates[0]
            if 'T' in nd.dofs:
                U[nd.dofs['T']] = 100.0 * nd.coordinates[0]

        path = self._tmp_path("_test_vtk_3d_thermal.vtu")
        VtkExporter(dom).export(path, U=U)
        mesh = meshio.read(path)
        self.assertEqual({b.type for b in mesh.cells}, {"hexahedron", "tetra", "quad"})
        for name in ("Sigma_XX", "Sigma_YY", "Sigma_ZZ", "Tau_XY", "Tau_YZ", "Tau_XZ", "Von_Mises"):
            self.assertIn(name, mesh.cell_data)
        blocks = {b.type: k for k, b in enumerate(mesh.cells)}
        gs = dom.elements[1].compute_gauss_state(U)
        self.assertAlmostEqual(mesh.cell_data["Sigma_XX"][blocks["hexahedron"]][0],
                               gs['stress'][:, 0].mean(), places=9)
        self.assertAlmostEqual(mesh.cell_data["Sigma_ZZ"][blocks["hexahedron"]][0],
                               gs['stress'][:, 2].mean(), places=9)
        self.assertIn("Temperature", mesh.point_data)
        self.assertAlmostEqual(float(np.max(mesh.point_data["Temperature"])), 600.0, places=9)
        self.assertIn("Flux", mesh.cell_data)
        # q = −k·∇T = −2·100 en x sobre la celda térmica.
        self.assertAlmostEqual(mesh.cell_data["Flux"][blocks["quad"]][0][0], -200.0, places=9)

    def test_tabla_de_celdas_cubre_todo_el_registro_de_elementos(self):
        """Todo elemento registrado con nodos estándar obtiene celda VTK."""
        from solidum.registry import ElementRegistry
        esperado = {
            "Tri3": "triangle", "Quad4": "quad", "Tri6": "triangle6",
            "Quad8": "quad8", "Quad9": "quad9",
            "Hex8": "hexahedron", "Tet4": "tetra", "Tet10": "tetra10",
            "Hex20": "hexahedron20", "Hex27": "hexahedron27",
            "Quad4Thermal": "quad", "Hex8Thermal": "hexahedron",
        }
        for name in ElementRegistry.names():
            cls = ElementRegistry.get(name)
            if cls.__module__.startswith("solidum.user."):
                continue        # módulos de usuario: los prueban sus tests

            class _Stub:
                pass
            stub = _Stub()
            stub.DOF_NAMES = getattr(cls, "DOF_NAMES", ["ux", "uy"])
            if hasattr(cls, "FLUX_DIM"):
                stub.FLUX_DIM = cls.FLUX_DIM
            if name in esperado:
                n_nodes = {"triangle": 3, "quad": 4, "triangle6": 6, "quad8": 8,
                           "quad9": 9, "hexahedron": 8, "tetra": 4, "tetra10": 10,
                           "hexahedron20": 20, "hexahedron27": 27}[esperado[name]]
                stub.nodes = [None] * n_nodes
                self.assertEqual(cell_type_for(stub), esperado[name], name)
            else:
                stub.nodes = [None, None]
                self.assertEqual(cell_type_for(stub), "line", name)

    def test_dominio_sin_celdas_soportadas_lanza(self):
        dom = Domain()
        for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1), (2, 2)], 1):
            dom.add_node(i, [float(x), float(y)])
        elem = Quad4(1, [dom.nodes[i] for i in (1, 2, 3, 4)], Elastic2D(E=1.0, nu=0.3))
        elem.nodes = [dom.nodes[i] for i in (1, 2, 3, 4, 5)]   # 5 nodos 2D: sin celda
        dom.elements[1] = elem
        dom.generate_equation_numbers()
        with self.assertRaises(ValueError):
            VtkExporter(dom).export(self._tmp_path("_test_vtk_vacio.vtu"), U=None)


if __name__ == '__main__':
    unittest.main()
