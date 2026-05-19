"""Tests del análisis dinámico transitorio no lineal Newmark+Newton — fase 4 ADR 0009.

Cubre:

1. **Recuperación del caso lineal**: con materiales elásticos puros (o
   plásticos con σ_y enorme), ``NewtonNewmarkSolver`` reproduce exactamente
   ``NewmarkSolver``. Garantiza que el residuo se anula en una iteración y
   no se introducen errores adicionales por la maquinaria de Newton.
2. **Oscilador 1 GDL elastoplástico libre**: con u_0 más allá del yield y
   vibración libre, la respuesta muestra disipación plástica (decaimiento)
   y α > 0 al final.
3. **Convergencia rápida en régimen plástico**: con tangente algorítmica
   consistente, el Newton interno converge en pocas iteraciones por paso.
4. **Cableado YAML** con ``solver: type: NewtonNewmarkSolver``.
"""
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

import fenix  # autodiscover
from fenix.core.domain import Domain
from fenix.elements.truss import Truss2D
from fenix.entry import run_transient, run_yaml
from fenix.materials.elastic import Elastic1D
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.linalg.lu import LUSolver
from fenix.math.solvers import HHTSolver, NewmarkSolver, NewtonNewmarkSolver
from fenix.math.solvers.diagnostics import SingularTangentError
from fenix.math.solvers.newmark import NewtonHHTSolver
from fenix.results import TransientResult


# Mismos parámetros que test_newmark.py para 1 GDL:
#   K_red = EA/L = 25, M_red = ρAL·2/6 = 1, ω = 5 rad/s.
E_1DOF = 25.0
RHO_1DOF = 3.0
A_1DOF = 1.0
L_1DOF = 1.0
OMEGA_1DOF = 5.0


def _build_1dof_oscillator(material) -> Domain:
    """Truss2D 1 elemento con K_red=25, M_red=1 usando el material dado."""
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L_1DOF, 0.0])
    dom.add_element(Truss2D(1, [n1, n2], material, A=A_1DOF))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers()
    return dom


def _free_dof(dom: Domain) -> int:
    return dom.nodes[2].dofs["ux"]


class TestLinearRecovery(unittest.TestCase):
    """Con material lineal o plástico-pero-sin-yield, NewtonNewmark = Newmark."""

    def test_undamped_pure_elastic_matches_NewmarkSolver(self):
        """Material Elastic1D: el residuo de Newton se anula en 1 iter; resultado idéntico."""
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 50.0
        t_end = 2.0 * T

        # Newmark lineal
        dom_lin = _build_1dof_oscillator(Elastic1D(E=E_1DOF, density=RHO_1DOF))
        ndof = dom_lin.total_dofs
        u0_lin = np.zeros(ndof); u0_lin[_free_dof(dom_lin)] = 0.1
        sol_lin = NewmarkSolver(
            Assembler(dom_lin), t_end=t_end, dt=dt,
            u0=u0_lin, F_func=None,
        )
        res_lin = sol_lin.solve()

        # NewtonNewmark con el mismo material lineal
        dom_nl = _build_1dof_oscillator(Elastic1D(E=E_1DOF, density=RHO_1DOF))
        u0_nl = np.zeros(dom_nl.total_dofs); u0_nl[_free_dof(dom_nl)] = 0.1
        sol_nl = NewtonNewmarkSolver(
            Assembler(dom_nl), t_end=t_end, dt=dt,
            u0=u0_nl, F_func=None,
            convergence=ConvergenceCriterion(rtol_force=1e-9, rtol_disp=1e-9),
            max_iter=5,
        )
        res_nl = sol_nl.solve()

        # Historiales deben coincidir
        np.testing.assert_allclose(res_nl.u_history, res_lin.u_history, atol=1e-9)
        np.testing.assert_allclose(res_nl.udot_history, res_lin.udot_history, atol=1e-8)

    def test_step_response_matches_NewmarkSolver(self):
        """Step response con material que nunca plastifica (σ_y enorme)."""
        F0 = 1.0
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 60.0
        t_end = 1.5 * T

        # Para que Elastoplastic1D sea efectivamente elástico, σ_y >> σ_max esperado.
        # σ_max ≈ E·ε_max ≈ E·(2·F0/K) = 25·0.08 = 2 → σ_y = 100 sobra.
        mat_plast = Elastoplastic1D(E=E_1DOF, sigma_y=100.0, H=0.0, density=RHO_1DOF)
        dom = _build_1dof_oscillator(mat_plast)
        ndof = dom.total_dofs
        dof = _free_dof(dom)
        F_func = lambda t: np.array([0.0, 0.0, F0, 0.0]) if False else (
            lambda t: (lambda v: (v.__setitem__(dof, F0) or v))(np.zeros(ndof))
        )(t)

        sol = NewtonNewmarkSolver(
            Assembler(dom), t_end=t_end, dt=dt,
            F_func=F_func,
            convergence=ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8),
        )
        res = sol.solve()

        # Sin plasticidad activada: α = 0 al final.
        elem = list(dom.elements.values())[0]
        self.assertEqual(elem.state.vars[0]['alpha'], 0.0)

        # Comparar contra Newmark lineal con Elastic1D
        dom_lin = _build_1dof_oscillator(Elastic1D(E=E_1DOF, density=RHO_1DOF))
        dof_lin = _free_dof(dom_lin)
        F_lin = lambda t: (lambda v: (v.__setitem__(dof_lin, F0) or v))(np.zeros(dom_lin.total_dofs))
        sol_lin = NewmarkSolver(Assembler(dom_lin), t_end=t_end, dt=dt, F_func=F_lin)
        res_lin = sol_lin.solve()

        np.testing.assert_allclose(res.u_history, res_lin.u_history, atol=1e-7)


class TestNonlinearElastoplastic1Dof(unittest.TestCase):
    """Oscilador 1 GDL con material plástico — disipación plástica activa."""

    def test_free_vibration_with_yield_dissipates_energy(self):
        """u_0 más allá del yield: α crece y la amplitud decrece con cada ciclo."""
        # Yield en σ = σ_y = 5 → ε_y = 5/25 = 0.2 → u_y = 0.2 (con L=1, ε=u/L)
        # u_0 = 0.5 > 0.2 → primer ciclo plastifica
        sigma_y = 5.0
        H = 5.0   # endurecimiento modesto
        mat = Elastoplastic1D(E=E_1DOF, sigma_y=sigma_y, H=H, density=RHO_1DOF)
        dom = _build_1dof_oscillator(mat)
        ndof = dom.total_dofs
        dof = _free_dof(dom)
        u0 = np.zeros(ndof); u0[dof] = 0.5

        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 100.0
        t_end = 5.0 * T

        sol = NewtonNewmarkSolver(
            Assembler(dom), t_end=t_end, dt=dt,
            u0=u0,
            convergence=ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7),
            max_iter=15,
        )
        res = sol.solve()

        # α al final debe ser > 0 (hubo plasticidad)
        elem = list(dom.elements.values())[0]
        alpha_final = elem.state.vars[0]['alpha']
        self.assertGreater(alpha_final, 0.0)

        # Amplitud final < amplitud inicial (disipación)
        u_dof = res.u_history[dof, :]
        amp_inicial = abs(u_dof[0])
        # tomar amplitud en último cuarto del análisis
        amp_final = np.max(np.abs(u_dof[3 * len(u_dof) // 4 :]))
        self.assertLess(amp_final, amp_inicial)

    def test_converges_in_few_iter_in_plastic_step(self):
        """Newton interno converge típicamente en ≤6 iter con tangente algorítmica."""
        sigma_y = 5.0
        H = 5.0
        mat = Elastoplastic1D(E=E_1DOF, sigma_y=sigma_y, H=H, density=RHO_1DOF)
        dom = _build_1dof_oscillator(mat)
        u0 = np.zeros(dom.total_dofs); u0[_free_dof(dom)] = 0.5

        T = 2.0 * math.pi / OMEGA_1DOF
        # dt y t_end pequeños — pocos pasos pero ya plastifica
        dt = T / 50.0
        t_end = 0.5 * T

        # max_iter=6: si no converge en 6 iter por paso, falla
        sol = NewtonNewmarkSolver(
            Assembler(dom), t_end=t_end, dt=dt,
            u0=u0,
            convergence=ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6),
            max_iter=6,
        )
        res = sol.solve()
        # Si llega aquí, todos los pasos convergieron en ≤6 iter
        self.assertEqual(res.converged, True)


class TestYamlPipeline(unittest.TestCase):
    """Cableado YAML con ``solver.type: NewtonNewmarkSolver``."""

    def test_run_yaml_dispatches_to_transient(self):
        yaml_body = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}

materials:
  - {id: 1, type: Elastoplastic1D, E: 25.0, sigma_y: 100.0, H: 0.0, density: 3.0}

elements:
  - {id: 1, type: Truss2D, nodes: [1, 2], material: 1, A: 1.0}

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0}
  - {node_id: 2, uy: 0.0}

solver:
  type: NewtonNewmarkSolver
  t_end: 1.0
  dt: 0.05
  max_iter: 10
  convergence:
    rtol_force: 1.0e-7
    rtol_disp: 1.0e-7
"""
        with tempfile.NamedTemporaryFile('w', suffix='.yaml', delete=False) as fh:
            fh.write(yaml_body)
            path = Path(fh.name)
        try:
            result = run_yaml(path)
            self.assertIsInstance(result, TransientResult)
            self.assertGreater(result.n_steps, 0)
            self.assertTrue(result.converged)
        finally:
            path.unlink()


class TestSingularTangentDiagnostic(unittest.TestCase):
    """Regresión de la auditoría H-4.1.

    Si la tangente dinámica J = M + γΔt·C + βΔt²·K_t resulta singular en LU
    (`RuntimeError` de `scipy.sparse.linalg.splu`), `NewtonNewmarkSolver` y
    `NewtonHHTSolver` deben lanzar `SingularTangentError` (subclase tipada de
    `SolverDivergedError`, ADR 0011) — **no** `RuntimeError` plano ni otra
    subclase de divergencia.

    Antes del fix de H-4.1 el flag `singular_tangent_seen` se inicializaba
    en `False` y nunca se actualizaba; `classify_divergence` jamás recibía
    `True` y `SingularTangentError` era inaccesible desde el subsistema
    dinámico no lineal.
    """

    def _make_solver_args(self):
        mat = Elastic1D(E=E_1DOF, density=RHO_1DOF)
        dom = _build_1dof_oscillator(mat)
        return dict(
            assembler=Assembler(dom),
            t_end=0.05, dt=0.01,
            linear_algebra="lu",  # forzar LU para que el mock intercepte
        )

    @staticmethod
    def _fail_after_first_call():
        """Side-effect que deja pasar la 1ª llamada (cálculo de aceleración
        inicial con `M_red`) y lanza `RuntimeError("singular")` desde la 2ª
        en adelante (primera resolución de `J` dentro del Newton del paso 1)."""
        real_solve = LUSolver.solve
        state = {"n": 0}

        def side_effect(self_mock, K, b):
            state["n"] += 1
            if state["n"] == 1:
                return real_solve(self_mock, K, b)
            raise RuntimeError("singular factor")

        return side_effect

    def test_newton_newmark_raises_singular_tangent_on_lu_failure(self):
        sol = NewtonNewmarkSolver(**self._make_solver_args())
        with patch.object(LUSolver, "solve", autospec=True,
                           side_effect=self._fail_after_first_call()):
            with self.assertRaises(SingularTangentError):
                sol.solve()

    def test_newton_hht_raises_singular_tangent_on_lu_failure(self):
        sol = NewtonHHTSolver(**self._make_solver_args())
        with patch.object(LUSolver, "solve", autospec=True,
                           side_effect=self._fail_after_first_call()):
            with self.assertRaises(SingularTangentError):
                sol.solve()


if __name__ == '__main__':
    unittest.main()
