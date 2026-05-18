"""Tests del ``HarmonicSolver`` — fase 6 del ADR 0009.

Cubre:

1. **Función de transferencia 1 GDL** contra la fórmula analítica
   ``H(ω) = 1 / (K − ω²M + iωC)``. Amplitud y fase coinciden a precisión
   numérica (residuo del solver lineal directo).
2. **Pico de resonancia** — máximo de ``|û|`` cae en ``ω ≈ ω_n``, con la
   altura controlada por el amortiguamiento Rayleigh.
3. **Coincidencia con Newmark forzado en régimen estacionario** — un
   Newmark suficientemente largo y con dt pequeño converge a la
   amplitud del HarmonicSolver dentro de tolerancia.
4. **Barrido logarítmico** — verificación de que ``scale="log"`` produce
   espaciado geométrico.
5. **Validaciones tempranas** — ``ValueError`` por argumentos
   incoherentes.
6. **Cableado YAML** end-to-end.
"""
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import fenix  # autodiscover
from fenix.core.domain import Domain
from fenix.elements.truss import Truss2D
from fenix.entry import run_harmonic, run_transient, run_yaml
from fenix.materials.elastic import Elastic1D
from fenix.math.assembly import Assembler
from fenix.math.solvers import HarmonicSolver, NewmarkSolver
from fenix.results import HarmonicResult


# ---------------------------------------------------------------------------
# Helper: oscilador 1 GDL idéntico al de test_newmark
# ---------------------------------------------------------------------------
# Con masa CONSISTENTE: K=25, M=1, ω_n=5 rad/s.

E_1DOF = 25.0
RHO_1DOF = 3.0
A_1DOF = 1.0
L_1DOF = 1.0
K_RED = 25.0   # EA/L
M_RED = 1.0    # ρ·A·L/3 (truss2D consistente, contribución del DOF libre)
OMEGA_N = 5.0  # √(K/M)


def _build_1dof() -> Domain:
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L_1DOF, 0.0])
    mat = Elastic1D(E=E_1DOF, density=RHO_1DOF)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=A_1DOF))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


def _free_dof(dom: Domain) -> int:
    return dom.nodes[2].dofs["ux"]


# ---------------------------------------------------------------------------
# 1 GDL — función de transferencia analítica
# ---------------------------------------------------------------------------


class TestHarmonic1Dof(unittest.TestCase):
    """``H(ω) = 1 / (K − ω²M + iωC)`` reproducida a precisión numérica."""

    def test_transfer_function_no_damping(self):
        """Sin amortiguamiento, ``H(ω) = 1/(K − ω²M)``: real puro lejos de
        resonancia, signo cambia al cruzar ω_n.
        """
        dom = _build_1dof()
        ndof = dom.total_dofs
        dof = _free_dof(dom)
        F_amp = np.zeros(ndof); F_amp[dof] = 10.0  # F̂ real = 10
        # Frecuencias lejos de resonancia (evitar div por cero).
        omega = np.array([1.0, 2.0, 7.0, 10.0, 20.0])

        result = run_harmonic(
            dom, solver=HarmonicSolver(
                Assembler(dom), omega=omega, F_amplitude=F_amp,
            )
        )

        u_dof = result.u_complex[dof, :]
        H_analytic = 10.0 / (K_RED - omega**2 * M_RED)
        np.testing.assert_allclose(u_dof.real, H_analytic, rtol=1e-10)
        np.testing.assert_allclose(u_dof.imag, 0.0, atol=1e-12)

    def test_transfer_function_rayleigh_damping(self):
        """Con amortiguamiento, ``û = F̂/(K − ω²M + iωC)`` complejo."""
        dom = _build_1dof()
        ndof = dom.total_dofs
        dof = _free_dof(dom)
        F_amp = np.zeros(ndof); F_amp[dof] = 1.0
        # ξ=5% en ω_n=5 con Rayleigh masa-pura: α=2ξω, β=0 → C = α·M
        alpha = 2.0 * 0.05 * OMEGA_N
        omega = np.linspace(1.0, 10.0, 19)

        result = run_harmonic(
            dom, solver=HarmonicSolver(
                Assembler(dom), omega=omega, F_amplitude=F_amp,
                rayleigh={"alpha": alpha, "beta": 0.0},
            )
        )

        # C analítico sobre el DOF libre: c = α·M_red (M lumped o consistent
        # — Truss2D consistente da M_red = ρA L · (2/6) = 1·2/6 + ?
        # En realidad M_red para el único DOF libre es el (2,2) reducido =
        # 2·ρAL/6 = 1.0. Con Rayleigh α-mass-only: c = α·1 = α.
        c = alpha * M_RED
        H_analytic = 1.0 / (K_RED - omega**2 * M_RED + 1j * omega * c)
        u_dof = result.u_complex[dof, :]
        np.testing.assert_allclose(u_dof, H_analytic, rtol=1e-9)

    def test_resonance_peak_location(self):
        """Pico de ``|û|`` está en ``ω ≈ ω_n`` con amortiguamiento ligero."""
        dom = _build_1dof()
        ndof = dom.total_dofs
        dof = _free_dof(dom)
        F_amp = np.zeros(ndof); F_amp[dof] = 1.0
        alpha = 2.0 * 0.02 * OMEGA_N   # ξ=2% — pico estrecho pero finito
        # Barrido fino alrededor de ω_n.
        omega = np.linspace(4.0, 6.0, 201)

        result = run_harmonic(
            dom, solver=HarmonicSolver(
                Assembler(dom), omega=omega, F_amplitude=F_amp,
                rayleigh={"alpha": alpha, "beta": 0.0},
            )
        )
        amp = result.amplitude()[dof, :]
        omega_peak = omega[np.argmax(amp)]
        # Para Rayleigh α-only, el pico exacto está en ω = ω_n·√(1 − 2ξ²)
        # (oscilador amortiguado), muy cercano a ω_n = 5 con ξ=2%.
        omega_peak_analytic = OMEGA_N * math.sqrt(1.0 - 2.0 * 0.02**2)
        self.assertAlmostEqual(omega_peak, omega_peak_analytic, delta=0.01)


# ---------------------------------------------------------------------------
# Coincidencia con Newmark forzado armónico en régimen estacionario
# ---------------------------------------------------------------------------


class TestHarmonicVsNewmark(unittest.TestCase):
    """Newmark forzado por ``F(t) = F₀·cos(ωt)`` tras descartar el
    transitorio inicial reproduce la amplitud de HarmonicSolver.
    """

    def test_steady_state_amplitude_matches(self):
        omega_drive = 7.0  # rad/s — fuera de resonancia, evita divergencia
        F0 = 1.0
        alpha = 2.0 * 0.10 * OMEGA_N   # ξ=10% para que el transitorio decaiga rápido
        ndof = 4
        dof = 2

        # HarmonicSolver: una sola frecuencia.
        dom_h = _build_1dof()
        F_amp = np.zeros(ndof); F_amp[dof] = F0
        result_h = run_harmonic(
            dom_h, solver=HarmonicSolver(
                Assembler(dom_h), omega=np.array([omega_drive]),
                F_amplitude=F_amp,
                rayleigh={"alpha": alpha, "beta": 0.0},
            )
        )
        u_harmonic_amp = float(np.abs(result_h.u_complex[dof, 0]))

        # NewmarkSolver: integración suficiente para decaer el transitorio.
        # Con ξ=10% y ω_n=5: τ = 1/(ξω) = 2 s. 8τ = 16 s deja el transitorio
        # por debajo de 0.03% del estado estacionario.
        dom_n = _build_1dof()
        T_drive = 2.0 * math.pi / omega_drive
        t_end = 16.0
        dt = T_drive / 100.0

        def F_func(t, F0_=F0, w=omega_drive, dof_=dof, ndof_=ndof):
            v = np.zeros(ndof_); v[dof_] = F0_ * math.cos(w * t); return v

        result_n = run_transient(
            dom_n, t_end=t_end, dt=dt, F_func=F_func,
            rayleigh={"alpha": alpha, "beta": 0.0},
        )
        # Tomar la amplitud máxima del último ciclo (ya en régimen estacionario).
        n_last = int(round(T_drive / dt))
        u_t = result_n.u_history[dof, -n_last:]
        u_newmark_amp = float(np.max(np.abs(u_t)))

        rel_err = abs(u_newmark_amp - u_harmonic_amp) / u_harmonic_amp
        self.assertLess(rel_err, 0.02,
                          f"Newmark steady-state |u| = {u_newmark_amp:.4e} "
                          f"vs Harmonic |û| = {u_harmonic_amp:.4e}; "
                          f"rel err = {rel_err:.2%}")


# ---------------------------------------------------------------------------
# Barrido logarítmico
# ---------------------------------------------------------------------------


class TestHarmonicSweep(unittest.TestCase):
    """Espaciado del barrido — ``linear`` y ``log``."""

    def test_log_scale_geometric_spacing(self):
        dom = _build_1dof()
        solver = HarmonicSolver(
            Assembler(dom),
            omega_min=1.0, omega_max=1000.0, n_omega=4, scale="log",
        )
        np.testing.assert_allclose(solver.omega,
                                     np.geomspace(1.0, 1000.0, 4))

    def test_linear_scale(self):
        dom = _build_1dof()
        solver = HarmonicSolver(
            Assembler(dom),
            omega_min=1.0, omega_max=10.0, n_omega=10, scale="linear",
        )
        np.testing.assert_allclose(solver.omega,
                                     np.linspace(1.0, 10.0, 10))


# ---------------------------------------------------------------------------
# Validaciones tempranas
# ---------------------------------------------------------------------------


class TestHarmonicContract(unittest.TestCase):
    """``ValueError`` con mensaje informativo en argumentos incoherentes."""

    def test_missing_omega_specification(self):
        dom = _build_1dof()
        with self.assertRaisesRegex(ValueError, "omega_min.*omega_max"):
            HarmonicSolver(Assembler(dom))

    def test_inverted_range(self):
        dom = _build_1dof()
        with self.assertRaisesRegex(ValueError, "rango inv"):
            HarmonicSolver(Assembler(dom),
                             omega_min=10.0, omega_max=1.0, n_omega=5)

    def test_zero_omega_in_explicit_array(self):
        dom = _build_1dof()
        with self.assertRaisesRegex(ValueError, "positivas"):
            HarmonicSolver(Assembler(dom),
                             omega=np.array([0.0, 1.0, 2.0]))

    def test_invalid_scale(self):
        dom = _build_1dof()
        with self.assertRaisesRegex(ValueError, "scale='cubic'"):
            HarmonicSolver(Assembler(dom),
                             omega_min=1.0, omega_max=10.0, scale="cubic")


# ---------------------------------------------------------------------------
# Cableado YAML
# ---------------------------------------------------------------------------


YAML_HARM_TEMPLATE = """\
nodes:
  - {{id: 1, coords: [0.0, 0.0]}}
  - {{id: 2, coords: [1.0, 0.0]}}

materials:
  - id: 1
    type: Elastic1D
    E: 25.0
    density: 3.0

elements:
  - {{id: 1, type: Truss2D, material: 1, A: 1.0, nodes: [1, 2]}}

boundary_conditions_by_node:
  - {{node_id: 1, ux: 0.0, uy: 0.0}}
  - {{node_id: 2, uy: 0.0}}

point_loads:
  - {{node_id: 2, ux: {force}}}

solver:
  type: HarmonicSolver
  omega_min: {omega_min}
  omega_max: {omega_max}
  n_omega: {n_omega}
  scale: linear
  rayleigh:
    alpha: {alpha}
    beta: 0.0
"""


class TestHarmonicYaml(unittest.TestCase):
    """Pipeline end-to-end YAML → HarmonicSolver → HarmonicResult."""

    def test_yaml_pipeline(self):
        alpha = 2.0 * 0.05 * OMEGA_N
        yaml_content = YAML_HARM_TEMPLATE.format(
            force=1.0, omega_min=1.0, omega_max=10.0, n_omega=19, alpha=alpha,
        )
        f = tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False, encoding="utf-8"
        )
        f.write(yaml_content); f.close()
        try:
            result = run_yaml(Path(f.name))
        finally:
            Path(f.name).unlink()

        self.assertIsInstance(result, HarmonicResult)
        self.assertEqual(result.n_omega, 19)
        # Comprueba que el pico está alrededor de ω_n=5 — los valores
        # extremos del barrido (lejos de resonancia) son menores.
        amp = result.amplitude()
        dof_amp = np.max(amp, axis=0)
        idx_peak = int(np.argmax(dof_amp))
        omega_peak = result.omega[idx_peak]
        # ω_n=5; barrido en pasos de 0.5; pico debe estar en 4.5–5.5.
        self.assertGreater(omega_peak, 4.0)
        self.assertLess(omega_peak, 6.0)


if __name__ == "__main__":
    unittest.main()
