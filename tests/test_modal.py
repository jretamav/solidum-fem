"""Tests del análisis modal — fase 1 del ADR 0009.

Cubre:

1. Validación contra solución analítica:
   - Barra axial empotrada-libre (Truss2D): ω_n = ((2n-1)π/(2L))·√(E/ρ).
   - Viga Bernoulli-Euler biapoyada (Frame2DEuler): ω_n = (nπ/L)²·√(EI/(ρA)).
2. Propiedades algebraicas del par (K, M):
   - M-ortonormalidad de los modos: ΦᵀMΦ = I.
   - Ortogonalidad respecto a K: ΦᵀKΦ = diag(ω²).
   - Componente nula del modo en DOFs prescritos por Dirichlet.
3. Contrato del solver:
   - ``density`` faltante en cualquier material → ``ValueError`` agregado.
   - ``lumping="lumped"`` aún no implementado → ``NotImplementedError``.
   - Cableado YAML: ``solver: type: ModalSolver``.
"""
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import fenix  # autodiscover
from fenix.core.domain import Domain
from fenix.elements.frame import Frame2DEuler
from fenix.elements.truss import Truss2D
from fenix.entry import run_modal, run_yaml
from fenix.materials.elastic import Elastic1D
from fenix.math.assembly import Assembler
from fenix.math.solvers import ModalSolver
from fenix.results import ModalResult


# Parámetros físicos comunes (acero):
E = 210.0e9
RHO = 7850.0
A = 1.0e-4
L_TOTAL = 1.0
N_ELEMS = 20


def _build_axial_bar(n_elems: int = N_ELEMS) -> Domain:
    """Barra empotrada-libre discretizada en `n_elems` Truss2D alineados con x.

    Todos los DOFs ``uy`` se fijan a cero para reducir el problema al modo
    axial puro — es el caso analítico contra el que se valida.
    """
    dom = Domain()
    nodes = []
    for i in range(n_elems + 1):
        x = i * (L_TOTAL / n_elems)
        nodes.append(dom.add_node(i + 1, [x, 0.0]))

    mat = Elastic1D(E=E, density=RHO)
    for i in range(n_elems):
        dom.add_element(Truss2D(i + 1, [nodes[i], nodes[i + 1]], mat, A=A))

    # Empotramiento en x=0, restricción transversal global a uy=0.
    nodes[0].fix_dof("ux", 0.0)
    for node in nodes:
        node.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


def _build_simply_supported_beam(n_elems: int = N_ELEMS) -> Domain:
    """Viga Bernoulli-Euler simplemente apoyada (apoyos a ambos extremos
    fijos en ``ux`` y ``uy``, rotación libre)."""
    I_beam = 8.33e-10
    dom = Domain()
    nodes = []
    for i in range(n_elems + 1):
        x = i * (L_TOTAL / n_elems)
        nodes.append(dom.add_node(i + 1, [x, 0.0]))

    mat = Elastic1D(E=E, density=RHO)
    for i in range(n_elems):
        dom.add_element(Frame2DEuler(i + 1, [nodes[i], nodes[i + 1]], mat,
                                       A=A, I=I_beam))

    # Apoyos simples a ambos extremos: ux=uy=0, rz libre.
    nodes[0].fix_dof("ux", 0.0); nodes[0].fix_dof("uy", 0.0)
    nodes[-1].fix_dof("ux", 0.0); nodes[-1].fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


class TestModalAxialBar(unittest.TestCase):
    """Barra axial empotrada-libre — frecuencias propias contra analítico."""

    def test_first_three_frequencies(self):
        dom = _build_axial_bar()
        result = run_modal(dom, n_modes=5)

        self.assertIsInstance(result, ModalResult)
        self.assertEqual(result.n_modes, 5)
        self.assertTrue(result.converged)

        # ω_n = ((2n-1)π/(2L))·√(E/ρ).
        wave_speed = math.sqrt(E / RHO)
        omegas_analytic = np.array(
            [(2 * n - 1) * math.pi / (2.0 * L_TOTAL) * wave_speed
             for n in (1, 2, 3)]
        )
        np.testing.assert_allclose(
            result.frequencies_rad[:3], omegas_analytic, rtol=0.01,
            err_msg="Primeros 3 modos axiales fuera de tolerancia 1%."
        )

    def test_frequencies_strictly_ascending(self):
        dom = _build_axial_bar()
        result = run_modal(dom, n_modes=5)
        diffs = np.diff(result.frequencies_rad)
        self.assertTrue(np.all(diffs > 0.0),
                         f"Frecuencias no monótonas crecientes: {result.frequencies_rad}")

    def test_periods_are_inverse_of_frequency_hz(self):
        dom = _build_axial_bar()
        result = run_modal(dom, n_modes=3)
        np.testing.assert_allclose(
            result.periods, 1.0 / result.frequencies_hz, rtol=1e-12
        )

    def test_mode_zero_at_clamped_dof(self):
        """Componente del modo en el nodo empotrado = 0 exacto."""
        dom = _build_axial_bar()
        result = run_modal(dom, n_modes=3)
        clamped_dof = dom.nodes[1].dofs["ux"]
        np.testing.assert_allclose(
            result.modes[clamped_dof, :], 0.0, atol=1e-14
        )


class TestModalSimplySupportedBeam(unittest.TestCase):
    """Viga Bernoulli-Euler biapoyada — flexión transversal contra analítico."""

    def test_first_three_flexural_frequencies(self):
        I_beam = 8.33e-10
        dom = _build_simply_supported_beam()
        result = run_modal(dom, n_modes=5)

        # ω_n = (nπ/L)²·√(EI/(ρA)).
        beta = math.sqrt(E * I_beam / (RHO * A))
        omegas_analytic = np.array(
            [(n * math.pi / L_TOTAL) ** 2 * beta for n in (1, 2, 3)]
        )
        np.testing.assert_allclose(
            result.frequencies_rad[:3], omegas_analytic, rtol=0.02,
            err_msg="Primeros 3 modos flexionales fuera de tolerancia 2%."
        )

    def test_modes_zero_at_supports(self):
        dom = _build_simply_supported_beam()
        result = run_modal(dom, n_modes=3)
        # uy = 0 en ambos apoyos para todos los modos.
        ny0_dof = dom.nodes[1].dofs["uy"]
        nyL_dof = dom.nodes[N_ELEMS + 1].dofs["uy"]
        np.testing.assert_allclose(result.modes[ny0_dof, :], 0.0, atol=1e-14)
        np.testing.assert_allclose(result.modes[nyL_dof, :], 0.0, atol=1e-14)


class TestModalAlgebraicProperties(unittest.TestCase):
    """ΦᵀMΦ = I y ΦᵀKΦ = diag(ω²) — independientes del problema."""

    def _assemble_KM_and_modes(self, dom: Domain, n_modes: int = 4):
        asm = Assembler(dom)
        solver = ModalSolver(asm, n_modes=n_modes)
        result = solver.solve()
        # Reconstrucción de K, M globales tras solve (cacheadas en asm).
        K = asm.K_global
        M = asm.assemble_mass_matrix()
        return result, K, M

    def test_mass_orthonormality_bar(self):
        dom = _build_axial_bar()
        result, _K, M = self._assemble_KM_and_modes(dom)
        Phi = result.modes
        MtM = Phi.T @ (M @ Phi)
        np.testing.assert_allclose(MtM, np.eye(result.n_modes), atol=1e-10)

    def test_stiffness_diagonalization_bar(self):
        dom = _build_axial_bar()
        result, K, _M = self._assemble_KM_and_modes(dom)
        Phi = result.modes
        KtK = Phi.T @ (K @ Phi)
        expected = np.diag(result.frequencies_rad ** 2)
        # Tolerancia relativa al autovalor máximo del bloque.
        scale = np.max(np.abs(expected))
        np.testing.assert_allclose(KtK, expected, atol=1e-8 * scale)

    def test_mass_orthonormality_beam(self):
        dom = _build_simply_supported_beam()
        result, _K, M = self._assemble_KM_and_modes(dom)
        Phi = result.modes
        MtM = Phi.T @ (M @ Phi)
        np.testing.assert_allclose(MtM, np.eye(result.n_modes), atol=1e-10)


class TestModalSolverContract(unittest.TestCase):
    """Errores agregados y opciones no implementadas (fase 1)."""

    def test_density_missing_raises(self):
        """Material sin density → ValueError listando los materiales."""
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [1.0, 0.0])
        mat = Elastic1D(E=E)  # ¡sin density!
        dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
        n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0); n2.fix_dof("uy", 0.0)
        dom.generate_equation_numbers(verbose=False)

        with self.assertRaisesRegex(ValueError, "Elastic1D"):
            run_modal(dom, n_modes=1)

    def test_lumped_not_implemented(self):
        """``lumping='lumped'`` lanza NotImplementedError (delegado al elemento)."""
        dom = _build_axial_bar()
        asm = Assembler(dom)
        solver = ModalSolver(asm, n_modes=1, lumping="lumped")
        with self.assertRaises(NotImplementedError):
            solver.solve()

    def test_runmodal_requires_solver_or_nmodes(self):
        dom = _build_axial_bar()
        with self.assertRaisesRegex(ValueError, "solver.*n_modes"):
            run_modal(dom)


class TestModalYamlPipeline(unittest.TestCase):
    """Cableado desde YAML: ``solver: type: ModalSolver``."""

    def _write_yaml(self, content: str) -> Path:
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False,
                                          encoding="utf-8")
        f.write(content); f.close()
        return Path(f.name)

    def test_yaml_modal_axial_bar(self):
        """Pipeline completo desde YAML reproduce la frecuencia fundamental
        de la barra axial empotrada-libre."""
        # 20 elementos a lo largo de [0, 1]; nodos numerados 1..21.
        nodes_yaml = "\n".join(
            f"  - {{id: {i + 1}, coords: [{i * (L_TOTAL / N_ELEMS):.6f}, 0.0]}}"
            for i in range(N_ELEMS + 1)
        )
        elems_yaml = "\n".join(
            f"  - {{id: {i + 1}, type: Truss2D, material: 1, "
            f"nodes: [{i + 1}, {i + 2}], A: {A}}}"
            for i in range(N_ELEMS)
        )
        bcs_yaml = "  - {node_id: 1, ux: 0.0, uy: 0.0}\n"
        bcs_yaml += "\n".join(
            f"  - {{node_id: {i + 1}, uy: 0.0}}"
            for i in range(1, N_ELEMS + 1)
        )

        yaml_content = f"""
nodes:
{nodes_yaml}

materials:
  - {{id: 1, type: Elastic1D, E: {E}, density: {RHO}}}

elements:
{elems_yaml}

boundary_conditions_by_node:
{bcs_yaml}

solver:
  type: ModalSolver
  n_modes: 3
"""
        path = self._write_yaml(yaml_content)
        try:
            result = run_yaml(str(path))
        finally:
            path.unlink()

        self.assertIsInstance(result, ModalResult)
        omega_1_analytic = math.pi / (2.0 * L_TOTAL) * math.sqrt(E / RHO)
        np.testing.assert_allclose(
            result.frequencies_rad[0], omega_1_analytic, rtol=0.01
        )


if __name__ == "__main__":
    unittest.main()
