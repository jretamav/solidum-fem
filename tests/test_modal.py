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
   - ``lumping="lumped"`` corre tras ADR 0009 fase 2 (cerrada 2026-05-18);
     la validación cuantitativa vive en ``tests/test_mass_lumping.py``.
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

from _modal_fixtures import (
    A_SECTION,
    E,
    L_TOTAL,
    N_ELEMS,
    RHO,
    build_axial_bar_2d,
    build_simply_supported_beam,
)


def _build_axial_bar(n_elems: int = N_ELEMS) -> Domain:
    """Atajo histórico: barra axial Truss2D con la configuración estándar."""
    return build_axial_bar_2d(Truss2D, n_elems=n_elems)


def _build_simply_supported_beam(n_elems: int = N_ELEMS) -> Domain:
    """Atajo histórico: viga Bernoulli simplemente apoyada Frame2DEuler."""
    return build_simply_supported_beam(Frame2DEuler, n_elems=n_elems)


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
        beta = math.sqrt(E * I_beam / (RHO * A_SECTION))
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
        dom.add_element(Truss2D(1, [n1, n2], mat, A=A_SECTION))
        n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0); n2.fix_dof("uy", 0.0)
        dom.generate_equation_numbers(verbose=False)

        with self.assertRaisesRegex(ValueError, "Elastic1D"):
            run_modal(dom, n_modes=1)

    def test_lumped_runs_after_phase2(self):
        """``lumping='lumped'`` corre tras ADR 0009 fase 2 y devuelve
        frecuencias positivas. La validación cuantitativa (recuperación de
        la fundamental, diagonalidad de M) vive en ``test_mass_lumping``.
        """
        dom = _build_axial_bar()
        asm = Assembler(dom)
        result = ModalSolver(asm, n_modes=1, lumping="lumped").solve()
        self.assertGreater(result.frequencies_rad[0], 0.0)

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
            f"nodes: [{i + 1}, {i + 2}], A: {A_SECTION}}}"
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


class TestModalFreeVibration(unittest.TestCase):
    """Superposición modal para vibración libre sin amortiguamiento (ADR 0009).

    Soluciones cerradas verificables:
    - ``u₀ = φ_n``, ``u̇₀ = 0`` ⇒ ``u(t) = cos(ω_n·t)·φ_n``.
    - ``u₀ = 0``, ``u̇₀ = φ_n`` ⇒ ``u(t) = (1/ω_n)·sin(ω_n·t)·φ_n``.
    - ``t = 0`` reproduce las condiciones iniciales en el subespacio modal.
    - Modos de cuerpo rígido (``ω = 0``) integran linealmente en el tiempo.
    """

    def _bar_and_M(self):
        dom = _build_axial_bar()
        asm = Assembler(dom)
        solver = ModalSolver(asm, n_modes=5)
        result = solver.solve()
        M = asm.assemble_mass_matrix()
        return result, M

    def test_initial_displacement_first_mode_pure_cosine(self):
        result, M = self._bar_and_M()
        phi1 = result.modes[:, 0]
        omega1 = result.frequencies_rad[0]
        t = np.array([0.0, 0.25, 0.5, 0.75, 1.0]) * (2.0 * np.pi / omega1)
        u = result.free_vibration(M, u0=phi1, u0_dot=np.zeros_like(phi1), t=t)
        expected = np.outer(phi1, np.cos(omega1 * t))
        np.testing.assert_allclose(u, expected, atol=1e-10)

    def test_initial_velocity_first_mode_pure_sine(self):
        result, M = self._bar_and_M()
        phi1 = result.modes[:, 0]
        omega1 = result.frequencies_rad[0]
        t = np.linspace(0.0, np.pi / omega1, 5)
        u = result.free_vibration(M, u0=np.zeros_like(phi1), u0_dot=phi1, t=t)
        expected = np.outer(phi1 / omega1, np.sin(omega1 * t))
        np.testing.assert_allclose(u, expected, atol=1e-10)

    def test_t_zero_recovers_initial_conditions_in_modal_subspace(self):
        """Cuando ``u₀`` y ``u̇₀`` viven en el subespacio modal, ``u(0)=u₀``
        y ``u̇(0)=u̇₀`` exactos. La velocidad inicial se verifica
        analíticamente vía la proyección modal: ``u̇(0) = Φ·Φᵀ·M·u̇₀``, que
        en el subespacio coincide con ``u̇₀`` por M-ortonormalidad."""
        result, M = self._bar_and_M()
        c_disp = np.array([1.0, 0.5, -0.3, 0.1, -0.05])
        c_vel = np.array([0.0, 2.0, 0.0, -1.0, 0.5])
        u0 = result.modes @ c_disp
        u0_dot = result.modes @ c_vel
        # u(0) = u₀.
        u = result.free_vibration(M, u0, u0_dot, t=np.array([0.0]))
        np.testing.assert_allclose(u[:, 0], u0, atol=1e-10)
        # u̇(0) reconstruida por la proyección modal exacta.
        u_dot_0 = result.modes @ (result.modes.T @ (M @ u0_dot))
        np.testing.assert_allclose(u_dot_0, u0_dot, atol=1e-10)

    def test_one_dof_analytical(self):
        """Oscilador armónico simple (1 GDL): ``u(t) = u₀·cos(ωt) + (u̇₀/ω)·sin(ωt)``.

        Test sintético sin pasar por ARPACK: instanciamos un ``ModalResult``
        a mano con un único modo ``φ=1``, ``ω=5 rad/s``, ``M=1``. La
        respuesta debe coincidir con la fórmula cerrada a precisión de
        máquina."""
        omega = 5.0
        result = ModalResult(
            frequencies_rad=np.array([omega]),
            frequencies_hz=np.array([omega / (2.0 * np.pi)]),
            periods=np.array([2.0 * np.pi / omega]),
            modes=np.array([[1.0]]),
            n_modes=1,
        )
        M = np.array([[1.0]])
        u0 = np.array([2.0]); u0_dot = np.array([3.0])
        t = np.linspace(0.0, 2.0 * np.pi / omega, 11)
        u = result.free_vibration(M, u0, u0_dot, t)
        expected = u0[0] * np.cos(omega * t) + (u0_dot[0] / omega) * np.sin(omega * t)
        np.testing.assert_allclose(u[0, :], expected, atol=1e-13)

    def test_superposition_of_two_modes(self):
        result, M = self._bar_and_M()
        phi1 = result.modes[:, 0]; phi2 = result.modes[:, 1]
        w1 = result.frequencies_rad[0]; w2 = result.frequencies_rad[1]
        u0 = phi1 + 2.0 * phi2
        t = np.linspace(0.0, 1.0e-3, 7)
        u = result.free_vibration(M, u0=u0, u0_dot=np.zeros_like(u0), t=t)
        expected = (np.outer(phi1, np.cos(w1 * t))
                    + 2.0 * np.outer(phi2, np.cos(w2 * t)))
        np.testing.assert_allclose(u, expected, atol=1e-10)

    def test_rigid_body_mode_integrates_linearly_in_time(self):
        """Si ``ω_n = 0``, la coordenada modal evoluciona como ``a + b·t``.

        Construcción directa de un ModalResult con un modo único rígido
        sobre una M identidad — evita necesitar un modelo singular en eigsh.
        """
        n_dof = 3
        phi = np.array([[1.0], [0.0], [0.0]])  # un modo rígido translacional en x.
        result = ModalResult(
            frequencies_rad=np.array([0.0]),
            frequencies_hz=np.array([0.0]),
            periods=np.array([np.inf]),
            modes=phi,
            n_modes=1,
            converged=True,
        )
        M = np.eye(n_dof)
        u0 = np.array([0.0, 0.0, 0.0])
        u0_dot = np.array([1.0, 0.0, 0.0])
        t = np.array([0.0, 1.0, 2.0, 5.0])
        u = result.free_vibration(M, u0, u0_dot, t)
        # u(t) = (a + b·t) · φ_1 con a=0, b=φ_1^T·M·u̇₀=1 ⇒ u_x(t) = t.
        expected_x = t
        np.testing.assert_allclose(u[0, :], expected_x, atol=1e-12)
        np.testing.assert_allclose(u[1:, :], 0.0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
