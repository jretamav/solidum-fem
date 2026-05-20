"""Patch test de MacNeal-Harder para los elementos sólidos 2D cuadráticos
del catálogo: Quad8, Quad9 y Tri6.

Extiende a la familia cuadrática el patch test que `test_patch_solid_2d.py`
cubre para Quad4 y Tri3. El test es la verificación canónica de
consistencia (Irons 1972, MacNeal-Harder 1985, Cook-Malkus-Plesha §6.4):
sobre una malla irregular, un campo lineal de desplazamientos prescrito
en el contorno debe reproducirse exactamente en los nodos interiores.

Los elementos cuadráticos pasan el patch test trivialmente cuando la base
Q2/Tri6 contiene la subbase Q1/Tri3, **pero** el ensamblaje, las
transformaciones isoparamétricas con midnodes desplazados y la imposición
de Dirichlet sobre nodos midside introducen oportunidades de bug que el
test caza. Distorsión clave: los midnodes interiores se mueven con el
corner central desplazado, garantizando que el mapeo isoparamétrico no
sea afín en ningún elemento.

Campo lineal usado (igual al de `test_patch_solid_2d.py`):

    u(x,y) = A1·x + A2·y
    v(x,y) = B1·x + B2·y

con (A1, A2, B1, B2) = (1e-3, 0.5e-3, 0.5e-3, 1e-3). Produce εxx = A1,
εyy = B2, γxy = A2 + B1 constantes en todo el dominio.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.core.material import Material
from solidum.elements.solid_2d import Quad8, Quad9, Tri6
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver


class _PlaneStress(Material):
    STRAIN_DIM = 3

    def __init__(self, E, nu):
        self.E, self.nu = E, nu
        fac = E / (1.0 - nu * nu)
        self.C = fac * np.array([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0],
        ])

    def compute_state(self, strain, state_vars=None):
        return self.C @ strain, self.C, state_vars


A1, A2 = 1.0e-3, 0.5e-3
B1, B2 = 0.5e-3, 1.0e-3
EXPECTED_STRAIN = np.array([A1, B2, A2 + B1])


def _linear_field(x, y):
    return A1 * x + A2 * y, B1 * x + B2 * y


# =============================================================================
# Quad8 / Quad9 — malla 2×2 con corner central interior desplazado
# =============================================================================

def _build_quad_higher_order(elem_cls, material, include_center: bool,
                              x_shift: float = 0.15, y_shift: float = 0.10):
    """Malla 2×2 Q8/Q9 sobre dominio [0,2]×[0,2] con corner central interior
    desplazado de (1, 1) a (1+x_shift, 1+y_shift).

    Los midnodes adyacentes al corner central también se desplazan
    (manteniéndolos en el midpoint geométrico de su arista) para que la
    malla siga siendo bien formada con el mapeo isoparamétrico no afín.

    Returns
    -------
    domain, ext_node_ids, int_node_ids, nodes_by_grid
    """
    dom = Domain()
    grid = {}
    nid = 0

    # Coordenadas regulares en grid 5x5; corner central (I=2, J=2) desplazado.
    base_x = lambda I: I * 0.5
    base_y = lambda J: J * 0.5

    for J in range(5):
        for I in range(5):
            is_center_param = (I % 2 == 1) and (J % 2 == 1)
            if is_center_param and not include_center:
                continue
            nid += 1
            x, y = base_x(I), base_y(J)
            if I == 2 and J == 2:
                x += x_shift; y += y_shift
            # Midnodes adyacentes al corner central → midpoint geométrico
            # entre sus corners reales (uno de ellos es el desplazado).
            elif (I, J) == (1, 2):
                x = 0.5 * (base_x(0) + (base_x(2) + x_shift))
                y = 0.5 * (base_y(2) + (base_y(2) + y_shift))
            elif (I, J) == (3, 2):
                x = 0.5 * ((base_x(2) + x_shift) + base_x(4))
                y = 0.5 * ((base_y(2) + y_shift) + base_y(2))
            elif (I, J) == (2, 1):
                x = 0.5 * ((base_x(2) + x_shift) + base_x(2))
                y = 0.5 * ((base_y(2) + y_shift) + base_y(0))
            elif (I, J) == (2, 3):
                x = 0.5 * ((base_x(2) + x_shift) + base_x(2))
                y = 0.5 * ((base_y(2) + y_shift) + base_y(4))
            grid[(I, J)] = dom.add_node(nid, [x, y])

    eid = 0
    for jy in range(2):
        for ix in range(2):
            I0, J0 = 2 * ix, 2 * jy
            c1 = grid[(I0, J0)]
            c2 = grid[(I0 + 2, J0)]
            c3 = grid[(I0 + 2, J0 + 2)]
            c4 = grid[(I0, J0 + 2)]
            m12 = grid[(I0 + 1, J0)]
            m23 = grid[(I0 + 2, J0 + 1)]
            m34 = grid[(I0 + 1, J0 + 2)]
            m41 = grid[(I0, J0 + 1)]
            conn = [c1, c2, c3, c4, m12, m23, m34, m41]
            if include_center:
                conn.append(grid[(I0 + 1, J0 + 1)])
            eid += 1
            dom.add_element(elem_cls(eid, conn, material, thickness=1.0))

    # Frontera externa: todos los nodos con I=0, I=4, J=0 o J=4.
    ext_ids = {(I, J) for (I, J) in grid if I in (0, 4) or J in (0, 4)}
    int_ids = set(grid.keys()) - ext_ids
    return dom, ext_ids, int_ids, grid


class _QuadHigherOrderPatchMixin:
    """Verificación común del patch test para Q8 y Q9."""

    ELEM_CLS = None
    INCLUDE_CENTER = False

    def setUp(self):
        self.material = _PlaneStress(E=1.0e6, nu=0.25)
        self.domain, self.ext_ids, self.int_ids, self.grid = (
            _build_quad_higher_order(
                self.ELEM_CLS, self.material,
                include_center=self.INCLUDE_CENTER,
            )
        )
        # Imponer el campo lineal en TODOS los DOFs de los nodos externos.
        for ij in self.ext_ids:
            node = self.grid[ij]
            x, y = node.coordinates
            u_x, u_y = _linear_field(x, y)
            node.fix_dof('ux', u_x)
            node.fix_dof('uy', u_y)
        self.domain.generate_equation_numbers(verbose=False)

    def test_interior_nodes_follow_linear_field(self):
        """Los nodos interiores (corner central + midnodes interiores + centers
        paramétricos en Q9) recuperan exactamente el campo lineal."""
        assembler = Assembler(self.domain)
        F_ext = np.zeros(self.domain.total_dofs)
        U = LinearSolver(assembler).solve(F_ext)

        for ij in self.int_ids:
            node = self.grid[ij]
            x, y = node.coordinates
            u_expected, v_expected = _linear_field(x, y)
            self.assertAlmostEqual(U[node.dofs['ux']], u_expected, places=9,
                msg=f"ux nodo interior {ij}")
            self.assertAlmostEqual(U[node.dofs['uy']], v_expected, places=9,
                msg=f"uy nodo interior {ij}")


class TestPatchQuad8(_QuadHigherOrderPatchMixin, unittest.TestCase):
    ELEM_CLS = Quad8
    INCLUDE_CENTER = False


class TestPatchQuad9(_QuadHigherOrderPatchMixin, unittest.TestCase):
    ELEM_CLS = Quad9
    INCLUDE_CENTER = True


# =============================================================================
# Tri6 — triangulación irregular de cuadrado con corner central desplazado
# =============================================================================

class TestPatchTri6(unittest.TestCase):
    """Patch test para Tri6: cuadrado [0,1]² triangulado en 4 Tri6 con
    corner central interior descentrado y midnodes desplazados acordemente.
    """

    def setUp(self):
        self.material = _PlaneStress(E=1.0e6, nu=0.25)
        dom = Domain()
        # 4 corner externos (esquinas del cuadrado) + 1 corner interior
        # descentrado + 4 midnodes externos (mitad de cada lado) + 4 midnodes
        # interiores (entre los 4 vértices externos y el corner interior).
        corners = {
            'sw': (0.0, 0.0),       # 1
            'se': (1.0, 0.0),       # 2
            'ne': (1.0, 1.0),       # 3
            'nw': (0.0, 1.0),       # 4
            'c':  (0.37, 0.42),     # 5 — interior descentrado
        }
        ids = {}
        nid = 0
        for name, xy in corners.items():
            nid += 1
            ids[name] = dom.add_node(nid, list(xy))

        # Midnodes externos: punto medio de cada arista exterior del cuadrado.
        ext_mid = {
            'se_m': (0.5, 0.0),     # entre sw-se
            'en_m': (1.0, 0.5),     # entre se-ne
            'nw_m': (0.5, 1.0),     # entre ne-nw
            'ws_m': (0.0, 0.5),     # entre nw-sw
        }
        for name, xy in ext_mid.items():
            nid += 1
            ids[name] = dom.add_node(nid, list(xy))

        # Midnodes interiores: midpoint entre cada corner exterior y el centro.
        cx, cy = corners['c']
        int_mid = {
            'sw_c': (0.5 * (0.0 + cx), 0.5 * (0.0 + cy)),
            'se_c': (0.5 * (1.0 + cx), 0.5 * (0.0 + cy)),
            'ne_c': (0.5 * (1.0 + cx), 0.5 * (1.0 + cy)),
            'nw_c': (0.5 * (0.0 + cx), 0.5 * (1.0 + cy)),
        }
        for name, xy in int_mid.items():
            nid += 1
            ids[name] = dom.add_node(nid, list(xy))

        # 4 Tri6 antihorarios. Numeración Solidum: [v0, v1, v2, m01, m12, m20].
        triangles = [
            # Inferior: sw, se, c | mids: se_m, se_c, sw_c
            ('sw', 'se', 'c',  'se_m', 'se_c', 'sw_c'),
            # Derecho:  se, ne, c | mids: en_m, ne_c, se_c
            ('se', 'ne', 'c',  'en_m', 'ne_c', 'se_c'),
            # Superior: ne, nw, c | mids: nw_m, nw_c, ne_c
            ('ne', 'nw', 'c',  'nw_m', 'nw_c', 'ne_c'),
            # Izquierdo: nw, sw, c | mids: ws_m, sw_c, nw_c
            ('nw', 'sw', 'c',  'ws_m', 'sw_c', 'nw_c'),
        ]
        eid = 0
        for tri in triangles:
            nodes = [ids[name] for name in tri]
            eid += 1
            dom.add_element(Tri6(eid, nodes, self.material, thickness=1.0))

        # Frontera externa: 4 corners exteriores + 4 midnodes externos.
        ext_names = ['sw', 'se', 'ne', 'nw',
                     'se_m', 'en_m', 'nw_m', 'ws_m']
        for name in ext_names:
            node = ids[name]
            x, y = node.coordinates
            u_x, u_y = _linear_field(x, y)
            node.fix_dof('ux', u_x)
            node.fix_dof('uy', u_y)

        dom.generate_equation_numbers(verbose=False)
        self.domain = dom
        self.ids = ids
        # Nodos interiores: corner central + 4 midnodes interiores.
        self.int_names = ['c', 'sw_c', 'se_c', 'ne_c', 'nw_c']

    def test_interior_nodes_follow_linear_field(self):
        assembler = Assembler(self.domain)
        F_ext = np.zeros(self.domain.total_dofs)
        U = LinearSolver(assembler).solve(F_ext)

        for name in self.int_names:
            node = self.ids[name]
            x, y = node.coordinates
            u_expected, v_expected = _linear_field(x, y)
            self.assertAlmostEqual(U[node.dofs['ux']], u_expected, places=9,
                msg=f"ux nodo {name}")
            self.assertAlmostEqual(U[node.dofs['uy']], v_expected, places=9,
                msg=f"uy nodo {name}")


if __name__ == '__main__':
    unittest.main()
