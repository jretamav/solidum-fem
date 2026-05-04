"""Tests funcionales del exportador VTK.

Construye dominios mínimos heterogéneos (sólidos 2D + barras + frames),
exporta a un archivo temporal y vuelve a leerlo con ``meshio`` para
verificar que:
- los puntos están en 3D y conservan las coordenadas,
- cada familia de elementos genera la celda VTK correcta (quad, triangle, line),
- los desplazamientos se escriben como vectores 3D y, si hay DOFs
  rotacionales, también las rotaciones,
- los campos de esfuerzos por celda llegan al fichero.

No depende de un solver — montamos U y F a mano.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import meshio  # noqa: E402

from fenix.core.domain import Domain
from fenix.core.material import Material
from fenix.core.node import Node
from fenix.elements.solid_2d import Quad4, Tri3
from fenix.elements.truss import Truss2D, Truss3D
from fenix.elements.frame import Frame2DEuler
from fenix.elements.frame3d import Frame3D
from fenix.utils.vtk_exporter import VtkExporter


class _Elastic2D(Material):
    STRAIN_DIM = 3
    PRIMARY_STATE_VAR = None

    def __init__(self, E=1000.0, nu=0.3):
        self.E = E
        self.nu = nu
        fac = E / (1.0 - nu * nu)
        self.C = fac * np.array([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0],
        ])

    def compute_state(self, strain, sv=None):
        return self.C @ strain, self.C, sv


class _Elastic1D(Material):
    STRAIN_DIM = 1

    def __init__(self, E=1000.0, nu=0.3):
        self.E = E
        self.nu = nu

    def compute_state(self, strain, sv=None):
        return self.E * strain, self.E, sv


def _add_element_to_domain(domain, elem):
    for nd in elem.nodes:
        domain.nodes.setdefault(nd.id, nd)
    domain.elements[elem.id] = elem


def _assign_global_dofs(domain):
    """Asigna índices globales a los DOFs registrados en cada nodo."""
    counter = 0
    for nd in domain.nodes.values():
        for dof_name in list(nd.dofs.keys()):
            nd.dofs[dof_name] = counter
            counter += 1
    domain.total_dofs = counter
    return counter


class TestVtkExporterCoverage(unittest.TestCase):
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

    def test_solids_y_truss(self):
        domain = Domain()
        mat2d = _Elastic2D()
        mat1d = _Elastic1D()

        # Quad4
        n1 = Node(1, [0.0, 0.0])
        n2 = Node(2, [1.0, 0.0])
        n3 = Node(3, [1.0, 1.0])
        n4 = Node(4, [0.0, 1.0])
        quad = Quad4(1, [n1, n2, n3, n4], mat2d)
        _add_element_to_domain(domain, quad)

        # Tri3 que comparte nodos con el quad
        n5 = Node(5, [2.0, 0.0])
        tri = Tri3(2, [n2, n5, n3], mat2d)
        _add_element_to_domain(domain, tri)

        # Truss2D entre dos nodos nuevos
        n6 = Node(6, [0.0, 2.0])
        n7 = Node(7, [2.0, 2.0])
        truss = Truss2D(3, [n6, n7], mat1d, A=0.01)
        _add_element_to_domain(domain, truss)

        n6.fix_dof('ux', 0.0); n6.fix_dof('uy', 0.0)
        _assign_global_dofs(domain)
        U = np.zeros(domain.total_dofs)
        # Imponer un desplazamiento en n2.uy para verificar que llega al VTK
        U[n2.dofs['uy']] = 0.05

        path = self._tmp_path("_test_vtk_solids_truss.vtu")
        VtkExporter(domain).export(path, U=U)

        mesh = meshio.read(path)
        self.assertEqual(mesh.points.shape, (7, 3))
        cell_types = {b.type for b in mesh.cells}
        self.assertEqual(cell_types, {"quad", "triangle", "line"})

        # Displacements 3D y vector
        disp = mesh.point_data["Displacements"]
        self.assertEqual(disp.shape, (7, 3))
        # Encontrar la fila correspondiente a n2 por coordenadas (1, 0, 0)
        idx_n2 = int(np.argmin(np.linalg.norm(mesh.points - [1.0, 0.0, 0.0], axis=1)))
        self.assertAlmostEqual(disp[idx_n2, 1], 0.05, places=12)

        # Sin DOFs rotacionales → no debe aparecer "Rotations"
        self.assertNotIn("Rotations", mesh.point_data)

        # Cell data: Sigma_XX presente para quad y triangle, 0 en line
        self.assertIn("Sigma_XX", mesh.cell_data)
        # meshio agrupa por bloque de celdas; comprobamos que hay tantos
        # bloques como tipos de celda
        self.assertEqual(len(mesh.cell_data["Sigma_XX"]), len(mesh.cells))

    def test_frames_3d_y_2d_generan_lines_y_rotations(self):
        domain = Domain()
        mat1d = _Elastic1D()

        # Frame2DEuler (DOFs ux, uy, rz)
        n1 = Node(1, [0.0, 0.0])
        n2 = Node(2, [1.0, 0.0])
        f2d = Frame2DEuler(1, [n1, n2], mat1d, A=1e-3, I=1e-6)
        _add_element_to_domain(domain, f2d)

        # Frame3D (DOFs ux, uy, uz, rx, ry, rz)
        n3 = Node(3, [0.0, 0.0, 0.0])
        n4 = Node(4, [0.0, 0.0, 1.0])
        f3d = Frame3D(2, [n3, n4], mat1d, A=1e-3, Iy=1e-6, Iz=1e-6, J=1e-6)
        _add_element_to_domain(domain, f3d)

        # Truss3D para confirmar que también se exporta como line
        n5 = Node(5, [1.0, 1.0, 1.0])
        truss3d = Truss3D(3, [n3, n5], mat1d, A=1e-3)
        _add_element_to_domain(domain, truss3d)

        _assign_global_dofs(domain)
        U = np.zeros(domain.total_dofs)
        U[n2.dofs['rz']] = 0.02      # rotación en frame 2D
        U[n4.dofs['ry']] = 0.03      # rotación en frame 3D
        U[n4.dofs['uz']] = 0.04      # desplazamiento z

        path = self._tmp_path("_test_vtk_frames.vtu")
        VtkExporter(domain).export(path, U=U)

        mesh = meshio.read(path)
        cell_types = {b.type for b in mesh.cells}
        # Tres elementos, todos lineales → un solo bloque "line" con 3 celdas
        self.assertEqual(cell_types, {"line"})

        # Displacements y Rotations vectoriales 3D
        disp = mesh.point_data["Displacements"]
        rot = mesh.point_data["Rotations"]
        self.assertEqual(disp.shape[1], 3)
        self.assertEqual(rot.shape[1], 3)

        idx_n4 = int(np.argmin(np.linalg.norm(mesh.points - [0.0, 0.0, 1.0], axis=1)))
        self.assertAlmostEqual(disp[idx_n4, 2], 0.04, places=12)
        self.assertAlmostEqual(rot[idx_n4, 1], 0.03, places=12)

        idx_n2 = int(np.argmin(np.linalg.norm(mesh.points - [1.0, 0.0, 0.0], axis=1)))
        self.assertAlmostEqual(rot[idx_n2, 2], 0.02, places=12)

    def test_nodal_stress_smoothing_uniaxial(self):
        """Suavizado nodal: con σ_xx = p uniforme por elemento, Sigma_XX_nodal
        debe valer p en todos los nodos del dominio."""
        domain = Domain()
        mat = _Elastic2D(E=1000.0, nu=0.3)
        # Dos Quad4 contiguos compartiendo el borde central
        n1 = Node(1, [0.0, 0.0]); n2 = Node(2, [1.0, 0.0])
        n3 = Node(3, [1.0, 1.0]); n4 = Node(4, [0.0, 1.0])
        n5 = Node(5, [2.0, 0.0]); n6 = Node(6, [2.0, 1.0])
        q1 = Quad4(1, [n1, n2, n3, n4], mat)
        q2 = Quad4(2, [n2, n5, n6, n3], mat)
        _add_element_to_domain(domain, q1)
        _add_element_to_domain(domain, q2)
        _assign_global_dofs(domain)

        # Campo u(x,y) = (p/E)·x  ⇒  σ_xx = p, σ_yy = ν·... espera, plane
        # stress: σ_yy = 0 sólo si ε_yy = -ν·ε_xx. Tomamos ε_xx=ε* y ε_yy=-ν·ε*
        # para garantizar σ_xx = E·ε* y σ_yy = 0.
        eps = 1e-3
        nu = mat.nu
        U = np.zeros(domain.total_dofs)
        for nd in domain.nodes.values():
            x, y = nd.coordinates[0], nd.coordinates[1]
            U[nd.dofs['ux']] = eps * x
            U[nd.dofs['uy']] = -nu * eps * y

        # Aplicamos compute_element_state y commit_state para que el σ
        # commiteado refleje el campo (el exporter lee σ committed).
        for elem in domain.elements.values():
            u_e = elem.get_local_displacements(U)
            elem.compute_element_state(u_e)
            elem.commit_state()

        path = self._tmp_path("_test_vtk_nodal_smoothing.vtu")
        VtkExporter(domain).export(path, U=U)
        mesh = meshio.read(path)

        sigma_xx = mesh.point_data["Sigma_XX_nodal"]
        sigma_yy = mesh.point_data["Sigma_YY_nodal"]
        # σ_xx esperado = E·ε* (plane stress con σ_yy = 0)
        expected = mat.E * eps
        self.assertTrue(np.allclose(sigma_xx, expected, atol=1e-8))
        self.assertTrue(np.allclose(sigma_yy, 0.0, atol=1e-8))

    def test_export_sin_U_no_rompe(self):
        domain = Domain()
        mat = _Elastic2D()
        n1 = Node(1, [0.0, 0.0])
        n2 = Node(2, [1.0, 0.0])
        n3 = Node(3, [1.0, 1.0])
        n4 = Node(4, [0.0, 1.0])
        Quad4(1, [n1, n2, n3, n4], mat)
        domain.nodes = {nd.id: nd for nd in [n1, n2, n3, n4]}
        # Inyectamos manualmente el elemento sin pasar por el assembler
        for elem in []:
            pass
        # pero el VtkExporter espera elementos en domain.elements
        from fenix.elements.solid_2d import Quad4 as Q
        domain.elements[1] = Q(1, [n1, n2, n3, n4], mat)
        _assign_global_dofs(domain)
        path = self._tmp_path("_test_vtk_pre.vtu")
        VtkExporter(domain).export(path, U=None, F_ext=None)
        mesh = meshio.read(path)
        self.assertEqual(mesh.points.shape[1], 3)
        # Desplazamientos todos a cero
        self.assertTrue(np.allclose(mesh.point_data["Displacements"], 0.0))


if __name__ == '__main__':
    unittest.main()
