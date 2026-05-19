"""Cook's membrane — benchmark cuantitativo de flexión + cortante para sólidos 2D.

Geometría trapezoidal canónica (Cook 1974, Bathe-Hughes):

    (0, 44) ──────── (48, 60)
       │                │
       │                │
    (0,  0) ──────── (48, 44)

- Lado izquierdo (x=0) totalmente empotrado.
- Lado derecho (x=48): fuerza total de cortante F=1 distribuida uniformemente
  como tracción ``t_y = 1/16`` (longitud del borde derecho = 16, espesor t=1).
- Material plane stress: E=1, ν=1/3.

**Cantidad medida**: desplazamiento vertical ``u_y`` en el punto medio del
borde derecho (48, 52). Valor de referencia convergido en la literatura
≈ 23.91 (Bathe; Belytschko; Hughes — pequeñas variaciones según fuente).

Este test cierra el hueco "Quad8/Quad9 sin benchmark cuantitativo" de la
matriz de validación (Fase B, 2026-05-19).
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.domain import Domain
from fenix.core.material import Material
from fenix.elements.solid_2d import Quad8, Quad9
from fenix.math.assembly import Assembler
from fenix.math.solvers import LinearSolver


class _PlaneStress(Material):
    """Material elástico isótropo plane stress mínimo para los tests."""
    STRAIN_DIM = 3

    def __init__(self, E: float, nu: float):
        self.E = E
        self.nu = nu
        fac = E / (1.0 - nu * nu)
        self.C = fac * np.array([
            [1.0, nu,  0.0],
            [nu,  1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0],
        ])

    def compute_state(self, strain, state_vars=None):
        return self.C @ strain, self.C, state_vars


def _cook_coords(xi: float, eta: float):
    """Mapeo (ξ,η)∈[0,1]² → trapezoid de Cook.

    Borde inferior: línea (0,0)→(48,44). Borde superior: (0,44)→(48,60).
    Ancho horizontal y crecimiento vertical son lineales con ξ.
    """
    x = 48.0 * xi
    y = 44.0 * xi + eta * (44.0 - 28.0 * xi)
    return x, y


def _build_cook_mesh(elem_cls, nx: int, ny: int, material, thickness: float = 1.0,
                     include_center: bool = False):
    """Malla estructurada de Cook con Q8 (o Q9 si ``include_center=True``).

    Devuelve ``(domain, left_corner_grid, right_edge_elements, mid_right_node)``.
    """
    domain = Domain()
    nid_grid = {}
    nid = 0
    # Grid en (2nx+1) × (2ny+1) índices paramétricos.
    for J in range(2 * ny + 1):
        for I in range(2 * nx + 1):
            is_center = (I % 2 == 1) and (J % 2 == 1)
            if is_center and not include_center:
                continue
            nid += 1
            xi = I / (2.0 * nx)
            eta = J / (2.0 * ny)
            x, y = _cook_coords(xi, eta)
            nid_grid[(I, J)] = domain.add_node(nid, [x, y])

    # Elementos: orden canónico [c1, c2, c3, c4, m12, m23, m34, m41 (, center)]
    eid = 0
    right_elements = []  # elementos cuyo borde 1 está en x=48
    for jy in range(ny):
        for ix in range(nx):
            I0 = 2 * ix
            J0 = 2 * jy
            c1 = nid_grid[(I0,   J0)]
            c2 = nid_grid[(I0+2, J0)]
            c3 = nid_grid[(I0+2, J0+2)]
            c4 = nid_grid[(I0,   J0+2)]
            m12 = nid_grid[(I0+1, J0)]
            m23 = nid_grid[(I0+2, J0+1)]
            m34 = nid_grid[(I0+1, J0+2)]
            m41 = nid_grid[(I0,   J0+1)]
            nodes = [c1, c2, c3, c4, m12, m23, m34, m41]
            if include_center:
                nodes.append(nid_grid[(I0+1, J0+1)])
            eid += 1
            elem = elem_cls(eid, nodes, material, thickness=thickness)
            domain.add_element(elem)
            if ix == nx - 1:
                right_elements.append(elem)

    # Empotramiento total en el borde izquierdo (I = 0).
    for J in range(2 * ny + 1):
        if (0, J) in nid_grid:
            n = nid_grid[(0, J)]
            n.fix_dof('ux', 0.0)
            n.fix_dof('uy', 0.0)

    # Nodo central del borde derecho (I = 2*nx, J = ny).
    mid_right = nid_grid[(2 * nx, ny)]

    domain.generate_equation_numbers(verbose=False)
    return domain, right_elements, mid_right


def _apply_right_edge_shear(domain, right_elements, total_force: float = 1.0):
    """Aplica ``F_total`` en y como tracción uniforme sobre el borde derecho.

    Borde derecho = edge 1 de cada Quad8/Quad9 rightmost.
    """
    edge_length = 16.0  # de (48,44) a (48,60)
    thickness = right_elements[0].thickness
    t_y = total_force / (edge_length * thickness)
    F_ext = np.zeros(domain.total_dofs)
    for elem in right_elements:
        f_elem = elem.compute_edge_traction(1, np.array([0.0, t_y]))
        # Ensamblar manualmente sobre los DOFs globales del elemento.
        dof_map = elem.get_global_dof_indices()
        for k_local, k_global in enumerate(dof_map):
            if k_global >= 0:
                F_ext[k_global] += f_elem[k_local]
    return F_ext


class TestCooksMembrane(unittest.TestCase):
    """Benchmark Cook 1974 — desplazamiento vertical del punto medio del
    borde derecho, contra valor convergido de referencia ≈ 23.91.
    """

    def setUp(self):
        self.material = _PlaneStress(E=1.0, nu=1.0 / 3.0)
        # Tolerancia: el valor "exacto" varía entre 23.81 y 23.96 según la
        # fuente; tomamos 23.91 (Belytschko) con tol ~0.5 (~2%).
        self.u_y_reference = 23.91
        self.tol_abs = 0.5

    def test_quad8_4x4(self):
        """Malla 4×4 Q8 (16 elementos, ~65 nodos): convergencia ≈ 0.5%."""
        domain, right_elems, mid_right = _build_cook_mesh(
            Quad8, nx=4, ny=4, material=self.material, thickness=1.0,
            include_center=False,
        )
        F = _apply_right_edge_shear(domain, right_elems, total_force=1.0)
        U = LinearSolver(Assembler(domain)).solve(F)
        u_y = U[mid_right.dofs['uy']]
        self.assertAlmostEqual(u_y, self.u_y_reference, delta=self.tol_abs,
                               msg=f"Cook Q8 4×4: u_y={u_y:.4f} (ref {self.u_y_reference})")

    def test_quad9_4x4(self):
        """Malla 4×4 Q9 (16 elementos, 81 nodos): convergencia similar a Q8."""
        domain, right_elems, mid_right = _build_cook_mesh(
            Quad9, nx=4, ny=4, material=self.material, thickness=1.0,
            include_center=True,
        )
        F = _apply_right_edge_shear(domain, right_elems, total_force=1.0)
        U = LinearSolver(Assembler(domain)).solve(F)
        u_y = U[mid_right.dofs['uy']]
        self.assertAlmostEqual(u_y, self.u_y_reference, delta=self.tol_abs,
                               msg=f"Cook Q9 4×4: u_y={u_y:.4f} (ref {self.u_y_reference})")

    def test_quad8_refinement_monotonic(self):
        """Refinamiento 2×2 → 4×4 → 6×6 reduce el error monótonamente hacia ref."""
        errors = []
        for n in (2, 4, 6):
            domain, right_elems, mid_right = _build_cook_mesh(
                Quad8, nx=n, ny=n, material=self.material, thickness=1.0,
                include_center=False,
            )
            F = _apply_right_edge_shear(domain, right_elems, total_force=1.0)
            U = LinearSolver(Assembler(domain)).solve(F)
            u_y = U[mid_right.dofs['uy']]
            errors.append(abs(u_y - self.u_y_reference))
        # Monotonía estricta del error.
        self.assertLess(errors[1], errors[0],
                        msg=f"4×4 no mejora sobre 2×2: errores={errors}")
        self.assertLess(errors[2], errors[1],
                        msg=f"6×6 no mejora sobre 4×4: errores={errors}")


if __name__ == '__main__':
    unittest.main()
