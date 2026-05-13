"""Tests de cargas distribuidas consistentes para elementos sólidos 2D.

Cubre Quad4 y Tri3:
- ``compute_body_load``: integral N^T b sobre el elemento.
- ``compute_edge_traction``: integral N^T t̄ sobre un borde.

Casos físicos:
1. Peso propio uniforme — la suma debe igualar ρg·A·t y el reparto es el
   consistente con las funciones de forma.
2. Tracción uniforme en un borde — la suma debe igualar t̄·L·t y el reparto
   entre los dos nodos del borde es L/2 cada uno.
3. Patch físico de tracción uniaxial: una capa de elementos en voladizo
   sometida a tracción de tracción ``(p, 0)`` reproduce ``σ_xx = p`` exacto.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.material import Material
from fenix.core.node import Node
from fenix.elements.solid_2d import Quad4, Tri3


class _Elastic2D(Material):
    """Material elástico isótropo plane stress, mínimo viable para tests."""
    STRAIN_DIM = 3

    def __init__(self, E: float, nu: float, density: float = 0.0):
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


def _make_quad4(coords, material, thickness=1.0):
    nodes = [Node(i + 1, list(c)) for i, c in enumerate(coords)]
    for k, nd in enumerate(nodes):
        nd.add_dof('ux'); nd.add_dof('uy')
        nd.dofs['ux'] = 2 * k
        nd.dofs['uy'] = 2 * k + 1
    return Quad4(1, nodes, material, thickness=thickness), nodes


def _make_tri3(coords, material, thickness=1.0):
    nodes = [Node(i + 1, list(c)) for i, c in enumerate(coords)]
    for k, nd in enumerate(nodes):
        nd.add_dof('ux'); nd.add_dof('uy')
        nd.dofs['ux'] = 2 * k
        nd.dofs['uy'] = 2 * k + 1
    return Tri3(1, nodes, material, thickness=thickness), nodes


class TestQuad4ConsistentLoads(unittest.TestCase):
    def setUp(self):
        self.mat = _Elastic2D(E=1000.0, nu=0.3)

    def test_body_load_uniform_regular_square(self):
        # Cuadrado 2x2 centrado en el origen, espesor 0.5, b = (0, -10)
        coords = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        elem, _ = _make_quad4(coords, self.mat, thickness=0.5)
        b = np.array([0.0, -10.0])
        f = elem.compute_body_load(b)
        A = 4.0
        total_y = f[1::2].sum()
        total_x = f[0::2].sum()
        self.assertAlmostEqual(total_x, 0.0, places=12)
        self.assertAlmostEqual(total_y, b[1] * A * 0.5, places=12)
        # Reparto simétrico (cuadrado regular): 1/4 por nodo
        per_node = b[1] * A * 0.5 / 4.0
        for i in range(4):
            self.assertAlmostEqual(f[2 * i + 1], per_node, places=12)

    def test_body_load_uniform_distorted(self):
        # Quad4 distorsionado: la suma debe seguir siendo b·A·t (invariante)
        coords = [(0, 0), (3, 0.2), (2.7, 2.5), (-0.3, 2.1)]
        elem, _ = _make_quad4(coords, self.mat, thickness=1.0)
        b = np.array([5.0, -7.0])
        f = elem.compute_body_load(b)
        # Área del cuadrilátero por la fórmula de shoelace
        x = np.array([c[0] for c in coords])
        y = np.array([c[1] for c in coords])
        A = 0.5 * abs(
            x[0] * (y[1] - y[3]) + x[1] * (y[2] - y[0])
            + x[2] * (y[3] - y[1]) + x[3] * (y[0] - y[2])
        )
        self.assertAlmostEqual(f[0::2].sum(), b[0] * A, places=10)
        self.assertAlmostEqual(f[1::2].sum(), b[1] * A, places=10)

    def test_edge_traction_uniform(self):
        coords = [(0, 0), (4, 0), (4, 2), (0, 2)]
        elem, _ = _make_quad4(coords, self.mat, thickness=0.5)
        # Borde 0 = (n0, n1), longitud 4, t̄ = (10, 0)
        t = np.array([10.0, 0.0])
        f = elem.compute_edge_traction(0, t)
        L = 4.0
        # Suma global: t̄·L·t
        self.assertAlmostEqual(f[0::2].sum(), t[0] * L * 0.5, places=12)
        self.assertAlmostEqual(f[1::2].sum(), 0.0, places=12)
        # Reparto: L/2·t·thick en cada uno de los 2 nodos del borde
        share = 0.5 * L * t[0] * 0.5
        self.assertAlmostEqual(f[0], share, places=12)  # nodo 0, fx
        self.assertAlmostEqual(f[2], share, places=12)  # nodo 1, fx
        # Otros nodos: cero
        self.assertAlmostEqual(f[4], 0.0, places=12)
        self.assertAlmostEqual(f[6], 0.0, places=12)

    def test_edge_index_validation(self):
        coords = [(0, 0), (1, 0), (1, 1), (0, 1)]
        elem, _ = _make_quad4(coords, self.mat)
        with self.assertRaises(ValueError):
            elem.compute_edge_traction(4, np.array([1.0, 0.0]))


class TestTri3ConsistentLoads(unittest.TestCase):
    def setUp(self):
        self.mat = _Elastic2D(E=1000.0, nu=0.3)

    def test_body_load_uniform(self):
        coords = [(0, 0), (3, 0), (0, 4)]
        elem, _ = _make_tri3(coords, self.mat, thickness=2.0)
        b = np.array([0.0, -9.81])
        f = elem.compute_body_load(b)
        A = 0.5 * 3.0 * 4.0  # área del triángulo rectángulo
        # Reparto: 1/3 por nodo (lineal)
        per_node = b[1] * A * 2.0 / 3.0
        for i in range(3):
            self.assertAlmostEqual(f[2 * i + 1], per_node, places=12)
        self.assertAlmostEqual(f[1::2].sum(), b[1] * A * 2.0, places=12)

    def test_edge_traction_uniform(self):
        coords = [(0, 0), (5, 0), (0, 3)]
        elem, _ = _make_tri3(coords, self.mat, thickness=1.0)
        # Borde 1 = (n1, n2), longitud sqrt(25 + 9) = sqrt(34)
        t = np.array([2.0, -1.0])
        f = elem.compute_edge_traction(1, t)
        L = float(np.sqrt(34.0))
        self.assertAlmostEqual(f[0::2].sum(), t[0] * L, places=12)
        self.assertAlmostEqual(f[1::2].sum(), t[1] * L, places=12)
        # Reparto entre nodos 1 y 2; el nodo 0 no recibe nada
        self.assertAlmostEqual(f[0], 0.0, places=12)
        self.assertAlmostEqual(f[1], 0.0, places=12)
        self.assertAlmostEqual(f[2], 0.5 * L * t[0], places=12)
        self.assertAlmostEqual(f[3], 0.5 * L * t[1], places=12)
        self.assertAlmostEqual(f[4], 0.5 * L * t[0], places=12)
        self.assertAlmostEqual(f[5], 0.5 * L * t[1], places=12)


class TestQuad4UniaxialTractionPatch(unittest.TestCase):
    """Patch físico: un Quad4 cuadrado con un borde libre traccionado
    uniformemente en x debe reproducir σ_xx = p con σ_yy = σ_xy = 0.

    Resolvemos el sistema 8x8 a mano (sin ensamblar globalmente): empotramos
    el borde izquierdo (ux=uy=0 en n0 y n3) y aplicamos la carga consistente
    en el borde derecho. Comparamos los desplazamientos contra la solución
    analítica u_x = p·x/E, u_y = -ν·p·y/E (plane stress).
    """

    def test_uniaxial(self):
        E = 1000.0
        nu = 0.3
        p = 5.0
        L = 2.0
        H = 1.0
        thick = 1.0

        mat = _Elastic2D(E=E, nu=nu)
        coords = [(0, 0), (L, 0), (L, H), (0, H)]
        elem, nodes = _make_quad4(coords, mat, thickness=thick)

        K_e = elem.compute_global_stiffness()
        # Borde 1 = (n1, n2), tracción (p, 0)
        F = elem.compute_edge_traction(1, np.array([p, 0.0]))

        # DOFs: 0,1 (n0), 2,3 (n1), 4,5 (n2), 6,7 (n3).
        # Para reproducir la solución uniaxial libre en Poisson:
        #   ux=0 en el borde izquierdo (n0, n3) → bloquea traslación x.
        #   uy=0 sólo en n0 → bloquea traslación y. n3 libre en y para
        #   permitir la contracción transversal -ν·p·H/E.
        # Restringidos: 0 (n0.ux), 1 (n0.uy), 6 (n3.ux). Libres: 2,3,4,5,7.
        free = [2, 3, 4, 5, 7]
        K_ff = K_e[np.ix_(free, free)]
        F_f = F[free]
        u_f = np.linalg.solve(K_ff, F_f)

        # Solución analítica plane stress: ux = p·x/E, uy = -ν·p·y/E.
        # Mapa free index → semántica:
        #   u_f[0]=n1.ux, u_f[1]=n1.uy, u_f[2]=n2.ux, u_f[3]=n2.uy, u_f[4]=n3.uy
        ux_L = p * L / E
        uy_H = -nu * p * H / E
        self.assertAlmostEqual(u_f[0], ux_L, places=10)   # n1 en (L,0)
        self.assertAlmostEqual(u_f[1], 0.0,  places=10)
        self.assertAlmostEqual(u_f[2], ux_L, places=10)   # n2 en (L,H)
        self.assertAlmostEqual(u_f[3], uy_H, places=10)
        self.assertAlmostEqual(u_f[4], uy_H, places=10)   # n3 en (0,H)

        # Reconstruir U y verificar el campo de esfuerzos en el elemento
        U = np.zeros(8)
        U[free] = u_f
        result = elem.compute_internal_forces(U)
        sigma = result['stress']
        self.assertAlmostEqual(sigma[0], p, places=8)     # σ_xx = p
        self.assertAlmostEqual(sigma[1], 0.0, places=8)   # σ_yy = 0
        self.assertAlmostEqual(sigma[2], 0.0, places=8)   # σ_xy = 0


if __name__ == '__main__':
    unittest.main()
