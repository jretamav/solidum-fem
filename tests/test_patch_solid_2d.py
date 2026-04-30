"""Patch test de MacNeal-Harder para los elementos sólidos 2D del catálogo.

Referencia
----------
MacNeal, R.H.; Harder, R.L. (1985). *A proposed standard set of problems to
test finite element accuracy*. Finite Elements in Analysis and Design, 1(1).
Recogido en NAFEMS como criterio de admisibilidad de elementos lineales.

Concepto
--------
Sobre una región rectangular dividida con al menos un nodo interior y
elementos distorsionados, se imponen desplazamientos prescritos en los
nodos del contorno correspondientes a un campo lineal:

    u(x,y) = a₀ + a₁·x + a₂·y
    v(x,y) = b₀ + b₁·x + b₂·y

Los nodos interiores quedan libres y se determinan por equilibrio. Un
elemento cumple el patch test si, tras resolver:

1. Los nodos interiores adoptan exactamente los desplazamientos del campo
   lineal (a redondeo).
2. La deformación en cada punto de integración de cada elemento es
   constante e igual a (a₁, b₂, a₂ + b₁) en notación de Voigt 2D.
3. La tensión correspondiente es uniforme en todo el dominio.

Es condición necesaria para convergencia bajo refinamiento de malla.
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.domain import Domain
from fenix.core.material import Material
from fenix.elements.solid_2d import Quad4, Tri3, _compute_kinematics, _compute_kinematics_tri3
from fenix.math.assembly import Assembler
from fenix.math.solvers import LinearSolver


class LinearElastic2DPlaneStress(Material):
    """Material elástico isótropo en tensión plana (STRAIN_DIM=3)."""

    STRAIN_DIM = 3

    def __init__(self, E: float, nu: float):
        self.E = E
        self.nu = nu
        fac = E / (1.0 - nu * nu)
        self.C = fac * np.array([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0],
        ])

    def compute_state(self, strain, state_var=None):
        return self.C @ strain, self.C, state_var


# Campo lineal de MacNeal-Harder: u = 1e-3·(x + y/2), v = 1e-3·(y + x/2).
# Genera εx = εy = 1e-3 y γxy = 1e-3 en todo el dominio.
A1, A2 = 1.0e-3, 0.5e-3
B1, B2 = 0.5e-3, 1.0e-3
EXPECTED_STRAIN = np.array([A1, B2, A2 + B1])


def _linear_field(x, y):
    return A1 * x + A2 * y, B1 * x + B2 * y


class TestPatchQuad4(unittest.TestCase):
    """Patch test estándar de MacNeal-Harder con cinco Quad4 distorsionados.

    Geometría: rectángulo (0, 0)-(0.24, 0.12), espesor 0.001.
    Cuatro nodos exteriores y cuatro nodos interiores con malla irregular.
    """

    # Coordenadas (mm escaladas a unidades arbitrarias).
    EXT_NODES = {
        1: (0.00, 0.00),
        2: (0.24, 0.00),
        3: (0.24, 0.12),
        4: (0.00, 0.12),
    }
    INT_NODES = {
        5: (0.04, 0.02),
        6: (0.18, 0.03),
        7: (0.16, 0.08),
        8: (0.08, 0.08),
    }
    # Cinco quads en sentido antihorario.
    CONNECTIVITY = {
        1: (1, 2, 6, 5),
        2: (2, 3, 7, 6),
        3: (3, 4, 8, 7),
        4: (4, 1, 5, 8),
        5: (5, 6, 7, 8),
    }

    def setUp(self):
        self.domain = Domain()
        for nid, xy in {**self.EXT_NODES, **self.INT_NODES}.items():
            self.domain.add_node(nid, list(xy))

        self.material = LinearElastic2DPlaneStress(E=1.0e6, nu=0.25)
        for eid, conn in self.CONNECTIVITY.items():
            nodes = [self.domain.get_node(nid) for nid in conn]
            self.domain.add_element(Quad4(eid, nodes, self.material, thickness=0.001))

        # Imponer el campo lineal en TODOS los DOFs de los nodos exteriores.
        for nid, (x, y) in self.EXT_NODES.items():
            u_x, u_y = _linear_field(x, y)
            node = self.domain.get_node(nid)
            node.fix_dof('ux', u_x)
            node.fix_dof('uy', u_y)

        self.domain.generate_equation_numbers()

    def test_interior_nodes_follow_linear_field(self):
        """Los nodos interiores se desplazan según el campo lineal impuesto."""
        assembler = Assembler(self.domain)
        F_ext = np.zeros(self.domain.total_dofs)
        U = LinearSolver(assembler).solve(F_ext)

        for nid, (x, y) in self.INT_NODES.items():
            u_expected, v_expected = _linear_field(x, y)
            node = self.domain.get_node(nid)
            self.assertAlmostEqual(U[node.dofs['ux']], u_expected, places=10,
                                   msg=f"ux nodo {nid}")
            self.assertAlmostEqual(U[node.dofs['uy']], v_expected, places=10,
                                   msg=f"uy nodo {nid}")

    def test_constant_strain_at_every_gauss_point(self):
        """ε y σ son constantes e iguales al valor analítico en cada punto de Gauss."""
        assembler = Assembler(self.domain)
        F_ext = np.zeros(self.domain.total_dofs)
        U = LinearSolver(assembler).solve(F_ext)

        expected_stress = self.material.C @ EXPECTED_STRAIN

        for elem in self.domain.elements.values():
            u_e = elem.get_local_displacements(U)
            coords = elem.get_coordinate_matrix(ndim=2)
            for (xi, eta) in elem.points:
                B, _ = _compute_kinematics(xi, eta, coords)
                strain = B @ u_e
                np.testing.assert_allclose(strain, EXPECTED_STRAIN, atol=1e-12,
                                           err_msg=f"ε elem {elem.id} en ({xi},{eta})")
                np.testing.assert_allclose(self.material.C @ strain,
                                           expected_stress, atol=1e-8,
                                           err_msg=f"σ elem {elem.id} en ({xi},{eta})")


class TestPatchTri3(unittest.TestCase):
    """Patch test para Tri3 sobre malla con un nodo interior y triángulos distorsionados.

    Cuadrado unitario con un nodo interior descentrado, triangulado en cuatro
    Tri3. Como CST reproduce campos lineales exactamente en cualquier triángulo
    bien formado, el test verifica además que el ensamblaje y la imposición de
    Dirichlet no introducen perturbaciones espurias.
    """

    EXT_NODES = {
        1: (0.0, 0.0),
        2: (1.0, 0.0),
        3: (1.0, 1.0),
        4: (0.0, 1.0),
    }
    INT_NODES = {
        5: (0.37, 0.42),  # interior descentrado
    }
    # Triángulos con orientación antihoraria.
    CONNECTIVITY = {
        1: (1, 2, 5),
        2: (2, 3, 5),
        3: (3, 4, 5),
        4: (4, 1, 5),
    }

    def setUp(self):
        self.domain = Domain()
        for nid, xy in {**self.EXT_NODES, **self.INT_NODES}.items():
            self.domain.add_node(nid, list(xy))

        self.material = LinearElastic2DPlaneStress(E=1.0e6, nu=0.25)
        for eid, conn in self.CONNECTIVITY.items():
            nodes = [self.domain.get_node(nid) for nid in conn]
            self.domain.add_element(Tri3(eid, nodes, self.material, thickness=1.0))

        for nid, (x, y) in self.EXT_NODES.items():
            u_x, u_y = _linear_field(x, y)
            node = self.domain.get_node(nid)
            node.fix_dof('ux', u_x)
            node.fix_dof('uy', u_y)

        self.domain.generate_equation_numbers()

    def test_interior_node_follows_linear_field(self):
        assembler = Assembler(self.domain)
        F_ext = np.zeros(self.domain.total_dofs)
        U = LinearSolver(assembler).solve(F_ext)

        for nid, (x, y) in self.INT_NODES.items():
            u_expected, v_expected = _linear_field(x, y)
            node = self.domain.get_node(nid)
            self.assertAlmostEqual(U[node.dofs['ux']], u_expected, places=10,
                                   msg=f"ux nodo {nid}")
            self.assertAlmostEqual(U[node.dofs['uy']], v_expected, places=10,
                                   msg=f"uy nodo {nid}")

    def test_constant_strain_at_every_element(self):
        assembler = Assembler(self.domain)
        F_ext = np.zeros(self.domain.total_dofs)
        U = LinearSolver(assembler).solve(F_ext)

        expected_stress = self.material.C @ EXPECTED_STRAIN

        for elem in self.domain.elements.values():
            u_e = elem.get_local_displacements(U)
            coords = elem.get_coordinate_matrix(ndim=2)
            B, _ = _compute_kinematics_tri3(coords)
            strain = B @ u_e
            np.testing.assert_allclose(strain, EXPECTED_STRAIN, atol=1e-12,
                                       err_msg=f"ε elem {elem.id}")
            np.testing.assert_allclose(self.material.C @ strain,
                                       expected_stress, atol=1e-8,
                                       err_msg=f"σ elem {elem.id}")


if __name__ == "__main__":
    unittest.main()
