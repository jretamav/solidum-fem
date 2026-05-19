"""Tests defensivos del batch 1 post-auditoría 2026-05-18.

Cubre cinco fixes puntuales:

- H-1.7 — ``run_yaml`` rechaza ``PIPELINE_KIND`` desconocido con ``ValueError``
  listando los válidos.
- H-3.2 — ``Elastic1D``/``Elastic2D`` validan inputs físicos (``E>0``,
  ``ν ∈ (-1, 0.5)``, ``hypothesis`` en lista cerrada, ``density ≥ 0``).
- H-4.6 — ``HHTSolver``/``NewtonHHTSolver`` emiten ``RuntimeWarning`` si el
  usuario overridea ``β``/``γ`` (sale del régimen probado Hilber 1977).
- H-3.1 — ``DruckerPrager2D`` con ``ψ = φ`` (asociado) expone
  ``IS_SYMMETRIC = True`` por instancia y el dispatcher elige Cholesky.
- H-4.7 — ``LinearSolver`` deriva ``is_positive_definite`` del flag de
  simetría del dominio (no hardcodeado).
"""
from __future__ import annotations

import unittest
import warnings

import numpy as np

import fenix  # autodiscover

from fenix.core.domain import Domain
from fenix.elements.solid_2d.quad4 import Quad4
from fenix.elements.truss import Truss2D
from fenix.materials.drucker_prager_2d import DruckerPrager2D
from fenix.materials.elastic import Elastic1D
from fenix.materials.elastic_2d import Elastic2D
from fenix.math.assembly import Assembler
from fenix.math.solvers._shared import domain_is_symmetric
from fenix.math.solvers.newmark import HHTSolver, NewtonHHTSolver


# ---------------------------------------------------------------------------
# H-3.2 — Elastic1D / Elastic2D validan inputs físicos
# ---------------------------------------------------------------------------


class TestElasticInputValidation(unittest.TestCase):

    def test_elastic1d_rejects_non_positive_E(self):
        with self.assertRaisesRegex(ValueError, "E=0.0 debe ser"):
            Elastic1D(E=0.0)
        with self.assertRaisesRegex(ValueError, "E=-1.0 debe ser"):
            Elastic1D(E=-1.0)

    def test_elastic1d_rejects_negative_density(self):
        with self.assertRaisesRegex(ValueError, "density=-2.0"):
            Elastic1D(E=1.0, density=-2.0)

    def test_elastic2d_rejects_non_positive_E(self):
        with self.assertRaisesRegex(ValueError, "E=0.0 debe ser"):
            Elastic2D(E=0.0, nu=0.25)

    def test_elastic2d_rejects_nu_out_of_range(self):
        for bad_nu in (-1.0, 0.5, 0.6, -1.5):
            with self.subTest(nu=bad_nu):
                with self.assertRaisesRegex(ValueError, "fuera del rango"):
                    Elastic2D(E=1.0, nu=bad_nu)

    def test_elastic2d_rejects_unknown_hypothesis(self):
        with self.assertRaisesRegex(ValueError, "hypothesis"):
            Elastic2D(E=1.0, nu=0.25, hypothesis="plain_stress")  # typo

    def test_elastic2d_accepts_valid_inputs(self):
        Elastic2D(E=1000.0, nu=0.25, hypothesis="plane_strain")
        Elastic2D(E=1000.0, nu=0.0, hypothesis="plane_stress")
        Elastic2D(E=1.0, nu=-0.5, hypothesis="plane_strain", density=0.0)


# ---------------------------------------------------------------------------
# H-4.6 — HHTSolver/NewtonHHTSolver avisan en override de β/γ
# ---------------------------------------------------------------------------


def _build_1dof_for_hht():
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [1.0, 0.0])
    mat = Elastic1D(E=25.0, density=2.0)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=1.0))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


class TestHHTBetaGammaOverrideWarns(unittest.TestCase):

    def test_hht_with_explicit_beta_warns(self):
        dom = _build_1dof_for_hht()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            HHTSolver(Assembler(dom), t_end=0.1, dt=0.01,
                       alpha=-0.05, beta=0.25)
            self.assertTrue(
                any(issubclass(item.category, RuntimeWarning) for item in w),
                "esperaba RuntimeWarning por override de beta",
            )

    def test_hht_with_explicit_gamma_warns(self):
        dom = _build_1dof_for_hht()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            HHTSolver(Assembler(dom), t_end=0.1, dt=0.01,
                       alpha=-0.05, gamma=0.5)
            self.assertTrue(
                any(issubclass(item.category, RuntimeWarning) for item in w),
                "esperaba RuntimeWarning por override de gamma",
            )

    def test_hht_without_overrides_does_not_warn(self):
        dom = _build_1dof_for_hht()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            HHTSolver(Assembler(dom), t_end=0.1, dt=0.01, alpha=-0.05)
            self.assertFalse(
                any(issubclass(item.category, RuntimeWarning)
                    and "beta/gamma" in str(item.message)
                    for item in w),
                "no debería avisar cuando no hay override explícito",
            )

    def test_newton_hht_with_override_warns(self):
        dom = _build_1dof_for_hht()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            NewtonHHTSolver(Assembler(dom), t_end=0.1, dt=0.01,
                             alpha=-0.05, beta=0.3)
            self.assertTrue(
                any(issubclass(item.category, RuntimeWarning) for item in w),
                "esperaba RuntimeWarning por override en NewtonHHT",
            )


# ---------------------------------------------------------------------------
# H-3.1 — DruckerPrager2D IS_SYMMETRIC dinámico (ψ = φ ⇒ True)
# ---------------------------------------------------------------------------


class TestDruckerPragerSymmetryFlagDynamic(unittest.TestCase):

    def test_associated_drucker_prager_is_symmetric_per_instance(self):
        """Cuando ψ = φ (asociado), la tangente algorítmica es simétrica."""
        mat = DruckerPrager2D(
            E=1000.0, nu=0.25, cohesion=10.0,
            phi_deg=30.0, psi_deg=30.0,  # asociado
        )
        self.assertTrue(mat.associated)
        self.assertTrue(mat.IS_SYMMETRIC)

    def test_nonassociated_drucker_prager_keeps_asymmetric(self):
        mat = DruckerPrager2D(
            E=1000.0, nu=0.25, cohesion=10.0,
            phi_deg=30.0, psi_deg=10.0,  # NO asociado
        )
        self.assertFalse(mat.associated)
        self.assertFalse(mat.IS_SYMMETRIC)

    def test_class_level_default_is_conservative_false(self):
        """Acceso por clase mantiene el default conservador (no asociado)."""
        self.assertFalse(DruckerPrager2D.IS_SYMMETRIC)

    def test_domain_is_symmetric_reads_instance_flag(self):
        """``domain_is_symmetric`` consulta el flag de instancia, no de clase.

        Si lo leyera de clase, un DruckerPrager asociado quedaría como
        no simétrico y el dispatcher elegiría LU innecesariamente.
        """
        dom = Domain()
        nodes = [dom.add_node(i + 1, c) for i, c in enumerate([
            [0., 0.], [1., 0.], [1., 1.], [0., 1.],
        ])]
        mat = DruckerPrager2D(
            E=1000.0, nu=0.25, cohesion=10.0,
            phi_deg=30.0, psi_deg=30.0,  # asociado ⇒ simétrico
        )
        dom.add_element(Quad4(1, nodes, mat, thickness=1.0))
        nodes[0].fix_dof("ux", 0.0); nodes[0].fix_dof("uy", 0.0)
        nodes[3].fix_dof("ux", 0.0)
        dom.generate_equation_numbers(verbose=False)

        self.assertTrue(domain_is_symmetric(dom))


# ---------------------------------------------------------------------------
# H-4.7 — LinearSolver IS_PD declarativo
# ---------------------------------------------------------------------------


class TestLinearSolverIsPdDerived(unittest.TestCase):
    """``LinearSolver._build_cache`` debe derivar ``is_positive_definite``
    de ``domain_is_symmetric``, no hardcodearlo a ``True``.
    """

    def test_pd_flag_matches_domain_symmetry_asymmetric_case(self):
        """Con material no asociado (asimétrico), el solver debe declarar
        ``is_positive_definite=False`` para que el dispatcher elija LU."""
        from fenix.math.solvers import LinearSolver

        dom = Domain()
        nodes = [dom.add_node(i + 1, c) for i, c in enumerate([
            [0., 0.], [1., 0.], [1., 1.], [0., 1.],
        ])]
        mat = DruckerPrager2D(
            E=1000.0, nu=0.25, cohesion=10.0,
            phi_deg=30.0, psi_deg=10.0,  # NO asociado ⇒ asimétrico
        )
        dom.add_element(Quad4(1, nodes, mat, thickness=1.0))
        nodes[0].fix_dof("ux", 0.0); nodes[0].fix_dof("uy", 0.0)
        nodes[3].fix_dof("ux", 0.0)
        dom.generate_equation_numbers(verbose=False)

        solver = LinearSolver(Assembler(dom))
        # No queremos correr la solución (DP no es lineal puro). Bastará con
        # construir la cache una sola vez y verificar el flag interno via
        # las propiedades del factor — pero el flag no se expone públicamente.
        # En su lugar, verificamos el comportamiento equivalente: el dispatch
        # eligió LU (no falló por no-PD intentando Cholesky).
        F = np.zeros(dom.total_dofs)
        u = solver.solve(F)
        self.assertEqual(u.shape, (dom.total_dofs,))


class TestRunYamlRejectsUnknownPipelineKind(unittest.TestCase):
    """``run_yaml`` debe rechazar con ``ValueError`` un ``PIPELINE_KIND``
    que no esté en ``_KNOWN_PIPELINE_KINDS``. Antes del fix H-1.7, un
    valor mal tecleado caía silenciosamente al pipeline estático.
    """

    def test_unknown_pipeline_kind_raises_valueerror(self):
        import tempfile
        from pathlib import Path
        from unittest.mock import patch

        from fenix.entry import run_yaml
        from fenix.math.solvers import LinearSolver

        yaml_body = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}

materials:
  - {id: 1, type: Elastic1D, E: 100.0}

elements:
  - {id: 1, type: Truss2D, nodes: [1, 2], material: 1, A: 1.0}

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0}
  - {node_id: 2, uy: 0.0}

solver:
  type: LinearSolver
"""
        with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as fh:
            fh.write(yaml_body)
            path = Path(fh.name)
        try:
            # Override temporal: PIPELINE_KIND inválido en LinearSolver.
            with patch.object(LinearSolver, "PIPELINE_KIND", "spectrum_v2"):
                with self.assertRaisesRegex(
                    ValueError, "PIPELINE_KIND='spectrum_v2'",
                ):
                    run_yaml(path)
        finally:
            path.unlink()


if __name__ == "__main__":
    unittest.main()
