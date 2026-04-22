"""Tests de los tipos de dato de fenix.results (ADR 0002)."""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.results import ElementForces, SolveResult


class TestElementForces(unittest.TestCase):
    def test_truss_ok(self):
        ef = ElementForces(kind="truss", components={"N": np.array([1.0, 1.0])})
        self.assertEqual(ef.at_node_i(), {"N": 1.0})
        self.assertEqual(ef.at_node_j(), {"N": 1.0})

    def test_frame2d_ok(self):
        comps = {
            "N": np.array([10.0, 10.0]),
            "V": np.array([2.0, -2.0]),
            "M": np.array([0.0, 5.0]),
        }
        ef = ElementForces(kind="frame2d", components=comps)
        self.assertEqual(ef.at_node_j()["M"], 5.0)

    def test_frame3d_keys(self):
        comps = {k: np.array([0.0, 0.0]) for k in ("N", "Vy", "Vz", "T", "My", "Mz")}
        ef = ElementForces(kind="frame3d", components=comps)
        self.assertEqual(set(ef.components), {"N", "Vy", "Vz", "T", "My", "Mz"})

    def test_unknown_kind_rejected(self):
        with self.assertRaises(ValueError):
            ElementForces(kind="shell", components={})  # type: ignore[arg-type]

    def test_invalid_component_key_rejected(self):
        with self.assertRaises(ValueError) as cm:
            ElementForces(kind="truss", components={"M": np.array([0.0, 0.0])})
        self.assertIn("Componentes inválidos", str(cm.exception))

    def test_wrong_shape_rejected(self):
        with self.assertRaises(ValueError):
            ElementForces(kind="truss", components={"N": np.array([1.0, 2.0, 3.0])})

    def test_frozen(self):
        ef = ElementForces(kind="truss", components={"N": np.array([1.0, 1.0])})
        with self.assertRaises(Exception):
            ef.kind = "cable"  # type: ignore[misc]


class TestSolveResult(unittest.TestCase):
    def test_defaults(self):
        n = 6
        res = SolveResult(
            U=np.zeros(n),
            F_applied=np.zeros(n),
            R=np.zeros(n),
        )
        self.assertTrue(res.converged)
        self.assertEqual(res.num_steps, 1)
        self.assertEqual(res.reactions_by_node, {})
        self.assertEqual(res.element_forces, {})

    def test_frozen(self):
        res = SolveResult(U=np.zeros(2), F_applied=np.zeros(2), R=np.zeros(2))
        with self.assertRaises(Exception):
            res.converged = False  # type: ignore[misc]


if __name__ == "__main__":
    unittest.main()
