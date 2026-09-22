"""Regresiones del lote de elementos de la auditoría 2026-09-22.

- Tolerancia del jacobiano relativa (adimensional): mallas finas en metros
  ya no se rechazan; la orientación invertida y el colapso sí.
- ``quadrature`` acepta clave o tupla en 2D, 3D y térmicos, y rechaza una
  regla de otra familia.
- Masa consistente y capacidad térmica integradas siempre con la regla
  completa, aunque ``K`` use integración reducida.
"""
import os
import tempfile
import unittest

import numpy as np

import solidum
from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4, Quad8, Quad9, Tri6
from solidum.elements.solid_3d import Hex8, Hex20, Tet4, Tet10
from solidum.elements.thermal import Hex8Thermal, Quad4Thermal
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.registry import QuadratureRegistry

M2 = Elastic2D(E=1000.0, nu=0.3, density=1.0)
M3 = Elastic3D(E=1000.0, nu=0.3, density=1.0)
HEX = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def _hex8(h=1.0, coords=HEX, quadrature=None):
    dom = Domain()
    for i, c in enumerate(coords, 1):
        dom.add_node(i, [v * h for v in c])
    kw = {} if quadrature is None else {"quadrature": quadrature}
    e = Hex8(1, [dom.nodes[i] for i in range(1, 9)], M3, **kw)
    dom.add_element(e); dom.generate_equation_numbers()
    return e


def _quad4(h=1.0, coords=((0, 0), (1, 0), (1, 1), (0, 1)), quadrature=None, material=M2):
    dom = Domain()
    for i, c in enumerate(coords, 1):
        dom.add_node(i, [v * h for v in c])
    kw = {} if quadrature is None else {"quadrature": quadrature}
    e = Quad4(1, [dom.nodes[i] for i in range(1, 5)], material, **kw)
    dom.add_element(e); dom.generate_equation_numbers()
    return e


class TestRelativeJacobianTolerance(unittest.TestCase):

    def test_fine_meshes_in_metres_are_accepted(self):
        for h in (1e-3, 1e-4, 1e-6):
            K = _hex8(h).compute_global_stiffness()
            self.assertTrue(np.all(np.isfinite(K)), h)
        for h in (1e-4, 1e-6, 1e-8):
            K = _quad4(h).compute_global_stiffness()
            self.assertTrue(np.all(np.isfinite(K)), h)
        dom = Domain()
        for i, c in enumerate([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)], 1):
            dom.add_node(i, [v * 1e-5 for v in c])
        Tet4(1, [dom.nodes[i] for i in range(1, 5)], M3).compute_global_stiffness()

    def test_inverted_and_collapsed_elements_still_raise(self):
        inverted = [HEX[i] for i in (0, 3, 2, 1, 4, 7, 6, 5)]
        with self.assertRaises(ValueError):
            _hex8(coords=inverted).compute_global_stiffness()
        with self.assertRaises(ValueError):
            _quad4(coords=((0, 0), (0, 1), (1, 1), (1, 0))).compute_global_stiffness()
        # Colapso total en el plano: cuadrilátero de área nula.
        with self.assertRaises(ValueError):
            _quad4(coords=((0, 0), (1, 0), (2, 0), (3, 0))).compute_global_stiffness()
        # Distorsión extrema relativa (nodo sobre la diagonal): det J <= 0 en
        # algún punto de Gauss con independencia de la escala.
        for h in (1.0, 1e-4):
            with self.assertRaises(ValueError):
                _quad4(h=h, coords=((0, 0), (1, 0), (0.1, 0.1), (0, 1))).compute_global_stiffness()


class TestQuadratureArgument(unittest.TestCase):

    def test_key_and_tuple_accepted_everywhere(self):
        e_key = _hex8(quadrature="hex_3x3x3")
        e_tup = _hex8(quadrature=QuadratureRegistry.get("hex_3x3x3"))
        np.testing.assert_allclose(e_key.compute_global_stiffness(), e_tup.compute_global_stiffness())
        self.assertEqual(e_key.quadrature_key, "hex_3x3x3")
        self.assertIsNone(e_tup.quadrature_key)
        q_key = _quad4(quadrature="2x2")
        q_tup = _quad4(quadrature=QuadratureRegistry.get("2x2"))
        np.testing.assert_allclose(q_key.compute_global_stiffness(), q_tup.compute_global_stiffness())
        # Térmico con tupla.
        dom = Domain()
        for i, c in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
            dom.add_node(i, [float(v) for v in c])
        mt = ThermalConduction(k=1.0, c=1.0, density=1.0, dim=2)
        Quad4Thermal(1, [dom.nodes[i] for i in range(1, 5)], mt,
                     quadrature=QuadratureRegistry.get("2x2"))

    def test_wrong_family_is_rejected(self):
        with self.assertRaises(ValueError):
            _hex8(quadrature="2x2")
        with self.assertRaises(ValueError):
            _quad4(quadrature="tri_3")
        with self.assertRaises(ValueError):
            _quad4(quadrature="hex_2x2x2")
        dom = Domain()
        for i, c in enumerate([(0, 0), (1, 0), (0, 1), (0.5, 0), (0.5, 0.5), (0, 0.5)], 1):
            dom.add_node(i, [float(v) for v in c])
        with self.assertRaises(ValueError):
            Tri6(1, [dom.nodes[i] for i in range(1, 7)], M2, quadrature="2x2")

    def test_yaml_hex8_with_quadrature_key_runs(self):
        text = "\n".join([
            "nodes:",
            *[f"  - {{id: {i}, coords: [{c[0]}, {c[1]}, {c[2]}]}}" for i, c in enumerate(HEX, 1)],
            "materials:",
            "  - {id: 1, type: Elastic3D, E: 1000.0, nu: 0.3}",
            "elements:",
            "  - {id: 1, type: Hex8, material: 1, nodes: [1,2,3,4,5,6,7,8], quadrature: hex_3x3x3}",
            "boundary_conditions:",
            *[f"  - {{node_id: {i}, ux: 0.0, uy: 0.0, uz: 0.0}}" for i in (1, 4, 5, 8)],
            "point_loads:",
            "  - {node_id: 2, ux: 1.0}",
            "solver:",
            "  type: LinearSolver",
            "",
        ])
        fd, path = tempfile.mkstemp(suffix=".yaml", dir=os.path.dirname(__file__))
        os.close(fd)
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(text)
            res = solidum.run_yaml(path)
            self.assertGreater(float(np.max(np.abs(res.U))), 0.0)
        finally:
            os.remove(path)


class TestMassWithReducedIntegration(unittest.TestCase):

    def test_quad4_and_hex8_mass_full_rank_with_reduced_stiffness_rule(self):
        M = _quad4(quadrature="1x1").compute_mass_matrix()
        self.assertEqual(np.linalg.matrix_rank(M), 8)
        self.assertAlmostEqual(M.sum() / 2.0, 1.0, places=12)   # masa total ρ·A·t = 1
        M = _hex8(quadrature="hex_1x1x1").compute_mass_matrix()
        self.assertEqual(np.linalg.matrix_rank(M), 24)
        self.assertAlmostEqual(M.sum() / 3.0, 1.0, places=12)

    def test_quad8_quad9_mass_full_rank_with_2x2(self):
        dom = Domain()
        c8 = [(0, 0), (1, 0), (1, 1), (0, 1), (0.5, 0), (1, 0.5), (0.5, 1), (0, 0.5)]
        for i, c in enumerate(c8, 1):
            dom.add_node(i, [float(v) for v in c])
        e8 = Quad8(1, [dom.nodes[i] for i in range(1, 9)], M2, quadrature="2x2")
        dom.add_node(9, [0.5, 0.5])
        e9 = Quad9(2, [dom.nodes[i] for i in range(1, 10)], M2, quadrature="2x2")
        dom.generate_equation_numbers()
        self.assertEqual(np.linalg.matrix_rank(e8.compute_mass_matrix()), 16)
        self.assertEqual(np.linalg.matrix_rank(e9.compute_mass_matrix()), 18)

    def test_thermal_capacity_full_rank_with_reduced_rule(self):
        dom = Domain()
        for i, c in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
            dom.add_node(i, [float(v) for v in c])
        mt = ThermalConduction(k=1.0, c=2.0, density=3.0, dim=2)
        e = Quad4Thermal(1, [dom.nodes[i] for i in range(1, 5)], mt, quadrature="1x1")
        dom.add_element(e); dom.generate_equation_numbers()
        C = e.compute_capacity_matrix(lumping="consistent")
        self.assertEqual(np.linalg.matrix_rank(C), 4)
        self.assertAlmostEqual(C.sum(), 6.0, places=12)   # ρc·V = 3·2·1
        dom3 = Domain()
        for i, c in enumerate(HEX, 1):
            dom3.add_node(i, [float(v) for v in c])
        mt3 = ThermalConduction(k=1.0, c=2.0, density=3.0, dim=3)
        e3 = Hex8Thermal(1, [dom3.nodes[i] for i in range(1, 9)], mt3, quadrature="hex_1x1x1")
        dom3.add_element(e3); dom3.generate_equation_numbers()
        self.assertEqual(np.linalg.matrix_rank(e3.compute_capacity_matrix(lumping="consistent")), 8)


if __name__ == "__main__":
    unittest.main()
