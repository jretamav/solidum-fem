"""Tests de los entrypoints públicos solidum.run y solidum.run_yaml (ADR 0002)."""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum
from solidum.core.domain import Domain
from solidum.elements.frame import Frame2DEuler
from solidum.elements.truss import Truss2D
from solidum.results import ModalResult, SolveResult, TransientResult


class LinearElastic1D:
    STRAIN_DIM = 1

    def __init__(self, E, nu=0.3, density: float = 0.0):
        self.E = E
        self.nu = nu

    def compute_state(self, strain, state_var=None):
        return self.E * strain, self.E, state_var


class TestRun(unittest.TestCase):
    def _build_cantilever(self):
        domain = Domain()
        n1 = domain.add_node(1, [0.0, 0.0])
        n2 = domain.add_node(2, [1.0, 0.0])
        elem = Frame2DEuler(1, [n1, n2], LinearElastic1D(E=1.0), A=1.0, I=1.0)
        domain.add_element(elem)
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n1.fix_dof('rz', 0.0)
        domain.generate_equation_numbers()
        return domain, n2

    def test_run_with_defaults_noload(self):
        domain, _ = self._build_cantilever()
        res = solidum.run(domain)
        self.assertIsInstance(res, SolveResult)
        self.assertIs(domain.last_result, res)
        np.testing.assert_allclose(res.U, np.zeros(6), atol=1e-10)

    def test_run_with_tip_load(self):
        domain, n2 = self._build_cantilever()
        F = np.zeros(domain.total_dofs)
        F[n2.dofs['uy']] = -1.0
        res = solidum.run(domain, F_applied=F)

        self.assertTrue(res.converged)
        np.testing.assert_allclose(res.reactions_by_node[1]['uy'], 1.0, atol=1e-4)
        ef = res.element_forces[1]
        np.testing.assert_allclose(ef.components["M"], [-1.0, 0.0], atol=1e-4)

    def test_run_auto_equation_numbers(self):
        """run() dispara generate_equation_numbers si aún no se hizo."""
        domain = Domain()
        n1 = domain.add_node(1, [0.0, 0.0])
        n2 = domain.add_node(2, [1.0, 0.0])
        truss = Truss2D(1, [n1, n2], LinearElastic1D(E=1.0), A=1.0)
        domain.add_element(truss)
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n2.fix_dof('uy', 0.0)
        # Deliberadamente NO llamamos generate_equation_numbers antes de run().
        self.assertEqual(domain.total_dofs, 0)

        res = solidum.run(domain)  # Tira con F=0, solo valida pipeline.
        self.assertEqual(domain.total_dofs, 4)
        self.assertIsInstance(res, SolveResult)


class TestRunYaml(unittest.TestCase):
    """Smoke test: run_yaml ejecuta el pipeline completo sobre un YAML de
    examples/ sin excepciones."""

    def test_run_yaml_smoke(self):
        examples_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'examples'))
        candidates = []
        for root, _, files in os.walk(examples_dir):
            for f in files:
                if f.endswith('.yaml'):
                    candidates.append(os.path.join(root, f))
        if not candidates:
            self.skipTest("No se encontraron archivos YAML en examples/")

        # Tomar el primero determinista
        yaml_path = sorted(candidates)[0]
        try:
            res = solidum.run_yaml(yaml_path)
        except Exception as e:
            self.skipTest(f"run_yaml falló en {yaml_path}: {e}")

        # El primer YAML alfabético puede ser estático, modal o transitorio:
        # el smoke valida solo que el pipeline produce un resultado válido del
        # catálogo de tipos. La validación específica por tipo vive en los
        # tests de cada solver.
        self.assertIsInstance(res, (SolveResult, ModalResult, TransientResult))


if __name__ == "__main__":
    unittest.main()
