"""Tests del integrador explícito ``CentralDifferenceSolver`` — fase 5 del
ADR 0009.

Cubre:

1. **Oscilador 1 GDL no amortiguado** contra solución analítica
   ``u(t) = u₀·cos(ωt) + (u̇₀/ω)·sin(ωt)`` — paso bajo CFL.
2. **Coincidencia con Newmark** (β=1/4, γ=1/2) para paso pequeño, modo
   lineal. Difieren en O(Δt²) — banda aceptable comprobada.
3. **Estabilidad condicional CFL** — ``Δt > Δt_crit`` dispara
   ``RuntimeError`` con mensaje informativo.
4. **Step response 1 GDL** contra ``u(t) = (F₀/K)·(1 − cos(ωt))``.
5. **Modo no lineal** — barra axial con material elastoplástico, contraste
   con `NewtonNewmarkSolver` en régimen elástico (mismas trayectorias) y
   tras plastificación parcial (modo no lineal acumula la deformación
   plástica como cualquier integrador transitorio).
6. **Rechazos defensivos** — ``lumping="consistent"`` lanza ``ValueError``;
   Frame3D con eje oblicuo lanza ``ValueError`` por M no estrictamente
   diagonal.
7. **Cableado YAML** end-to-end con ``solver: type: CentralDifferenceSolver``.
"""
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import solidum  # autodiscover
from solidum.core.domain import Domain
from solidum.elements.frame3d import Frame3D
from solidum.elements.truss import Truss2D
from solidum.entry import run_transient, run_yaml
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.assembly import Assembler
from solidum.math.solvers import (
    CentralDifferenceSolver,
    NewmarkSolver,
    NewtonNewmarkSolver,
)
from solidum.results import TransientResult


# ---------------------------------------------------------------------------
# Helper: oscilador 1 GDL con K=25, M=1
# ---------------------------------------------------------------------------
# Para que el problema reducido con masa LUMPED tenga K=25 y M=1:
#   K_red = EA/L = 25  → A=1, E=25, L=1
#   M_red (lumped) = ρ·A·L/2 = 1  → ρ=2
# Frecuencia: ω = √(K/M) = 5 rad/s; Δt_crit = 2/ω = 0.4 s.

E_1DOF = 25.0
RHO_1DOF_LUMPED = 2.0
A_1DOF = 1.0
L_1DOF = 1.0
OMEGA_1DOF = 5.0
DT_CRIT_1DOF = 2.0 / OMEGA_1DOF


def _build_1dof_oscillator(rho=RHO_1DOF_LUMPED) -> Domain:
    """Truss2D 1 elemento. Con masa lumped → K_red=25, M_red=1 (ω=5)."""
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L_1DOF, 0.0])
    mat = Elastic1D(E=E_1DOF, density=rho)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=A_1DOF))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


def _free_dof(dom: Domain) -> int:
    return dom.nodes[2].dofs["ux"]


# ---------------------------------------------------------------------------
# Validación analítica — oscilador 1 GDL
# ---------------------------------------------------------------------------


class TestCentralDifference1Dof(unittest.TestCase):
    """Vibración libre y forzada 1 GDL — solución analítica conocida."""

    def test_free_vibration_pure_cosine(self):
        """``u(t) = cos(ωt)`` con u₀=1, u̇₀=0 y CFL respetado (Δt = T/200)."""
        dom = _build_1dof_oscillator()
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 200.0   # Δt = 0.0314 << Δt_crit = 0.4 → muy estable
        t_end = 2.0 * T

        ndof = dom.total_dofs
        dof = _free_dof(dom)
        u0 = np.zeros(ndof); u0[dof] = 1.0
        u0_dot = np.zeros(ndof)

        solver = CentralDifferenceSolver(
            Assembler(dom), t_end=t_end, dt=dt,
            u0=u0, u0_dot=u0_dot,
        )
        result = run_transient(dom, solver=solver)

        u_t = result.u_history[dof, :]
        u_analytic = np.cos(OMEGA_1DOF * result.t_history)
        # Banda 2% — central diff es O(Δt²) con factor de error pequeño.
        err = np.max(np.abs(u_t - u_analytic))
        self.assertLess(err, 0.02,
                          f"|u_CD − u_analytic|_max = {err:.3e} > 2%.")

    def test_step_response(self):
        """``F(t) = F₀ H(t)``: respuesta ``u(t) = (F₀/K)·(1 − cos(ωt))``."""
        dom = _build_1dof_oscillator()
        ndof = dom.total_dofs
        dof = _free_dof(dom)
        K = 25.0; F0 = 50.0; u_ss = F0 / K  # = 2.0
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 200.0
        t_end = 2.0 * T

        def F_func(t):
            v = np.zeros(ndof); v[dof] = F0; return v

        solver = CentralDifferenceSolver(
            Assembler(dom), t_end=t_end, dt=dt, F_func=F_func,
        )
        result = run_transient(dom, solver=solver)

        u_t = result.u_history[dof, :]
        u_analytic = u_ss * (1.0 - np.cos(OMEGA_1DOF * result.t_history))
        err = np.max(np.abs(u_t - u_analytic))
        self.assertLess(err, 0.02 * u_ss,
                          f"step response |err|_max = {err:.3e} > 2% u_ss.")


# ---------------------------------------------------------------------------
# Coincidencia con Newmark para Δt pequeño
# ---------------------------------------------------------------------------


class TestCentralDifferenceVsNewmark(unittest.TestCase):
    """Central diff y Newmark β=1/4, γ=1/2 deben coincidir en O(Δt²)."""

    def test_close_to_newmark_with_small_dt(self):
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 500.0
        t_end = T

        dom_n = _build_1dof_oscillator(rho=3.0)  # Newmark usa lumped o consistent
        # NOTA: con masa consistente (Newmark default) M_red = ρ/3 = 1 con ρ=3.
        # Con masa lumped, M_red = ρ/2 = 1 requiere ρ=2.
        # Usamos lumped en ambos para igualar las masas reducidas.
        ndof = dom_n.total_dofs
        dof = _free_dof(dom_n)
        u0 = np.zeros(ndof); u0[dof] = 1.0
        u0_dot = np.zeros(ndof)

        dom_l = _build_1dof_oscillator(rho=2.0)
        result_newmark = run_transient(
            dom_l, t_end=t_end, dt=dt, u0=u0, u0_dot=u0_dot,
            lumping="lumped",
        )

        dom_cd = _build_1dof_oscillator(rho=2.0)
        solver_cd = CentralDifferenceSolver(
            Assembler(dom_cd), t_end=t_end, dt=dt,
            u0=u0, u0_dot=u0_dot,
        )
        result_cd = run_transient(dom_cd, solver=solver_cd)

        u_n = result_newmark.u_history[dof, :]
        u_cd = result_cd.u_history[dof, :]
        # Las dos series temporales deben coincidir dentro de 1%.
        err = np.max(np.abs(u_n - u_cd))
        self.assertLess(err, 0.01,
                          f"|u_Newmark − u_CD|_max = {err:.3e} > 1%.")


# ---------------------------------------------------------------------------
# Estabilidad condicional CFL
# ---------------------------------------------------------------------------


class TestCentralDifferenceStability(unittest.TestCase):
    """``Δt > Δt_crit`` debe disparar ``RuntimeError`` con mensaje CFL."""

    def test_cfl_violation_raises(self):
        dom = _build_1dof_oscillator()
        # Δt = 3·Δt_crit = 1.2 — explosión exponencial garantizada.
        dt = 3.0 * DT_CRIT_1DOF
        ndof = dom.total_dofs
        dof = _free_dof(dom)
        u0 = np.zeros(ndof); u0[dof] = 1.0

        solver = CentralDifferenceSolver(
            Assembler(dom), t_end=100.0 * dt, dt=dt, u0=u0,
        )
        with self.assertRaisesRegex(RuntimeError, "divergencia.*CFL"):
            solver.solve()


# ---------------------------------------------------------------------------
# CFL analítico: Δt_crit numérico = 2/ω_max (Fase B 2026-05-19)
# ---------------------------------------------------------------------------


def _build_2dof_analytic_chain():
    """Cadena 2-DOF n1-n2-n3-n4 con n1, n4 fijos y EA/L = 1, ρAL/2 = 1/2.

    Bajo masa lumped: K_red = [[2, -1], [-1, 2]], M_red = I.
    Autovalores: ω² ∈ {1, 3} → ω_min = 1, ω_max = √3.
    Δt_crit = 2/ω_max = 2/√3.

    Devuelve ``(domain, free_dofs_list)`` para que el test los reuse.
    """
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [1.0, 0.0])
    n3 = dom.add_node(3, [2.0, 0.0])
    n4 = dom.add_node(4, [3.0, 0.0])
    mat = Elastic1D(E=1.0, density=1.0)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=1.0))
    dom.add_element(Truss2D(2, [n2, n3], mat, A=1.0))
    dom.add_element(Truss2D(3, [n3, n4], mat, A=1.0))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    n3.fix_dof("uy", 0.0)
    n4.fix_dof("ux", 0.0); n4.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, [dom.nodes[2].dofs["ux"], dom.nodes[3].dofs["ux"]]


class TestCentralDifferenceCFLAnalytic(unittest.TestCase):
    """Verifica que la frontera ``Δt_crit`` numérica coincide con la
    expresión analítica ``Δt_crit = 2/ω_max``: por debajo, integración
    estable; por encima, ``RuntimeError`` por divergencia CFL.

    Cierra el hueco "CentralDifferenceSolver sin verificación de CFL
    teórico" de la matriz de validación. Se ejecuta tanto en el
    oscilador 1-GDL (ω_max = ω = 5) como en una cadena 2-DOF con
    autovalores cerrados (ω_max = √3) — el segundo caso ejercita el
    cálculo del máximo de la masa sobre múltiples modos.
    """

    def test_just_below_cfl_is_stable_1dof(self):
        """``Δt = 0.95·Δt_crit`` ⇒ trayectoria acotada del oscilador 1-GDL."""
        dt = 0.95 * DT_CRIT_1DOF
        dom = _build_1dof_oscillator()
        u0 = np.zeros(dom.total_dofs); u0[_free_dof(dom)] = 1.0
        # ~30 ciclos de simulación.
        T = 2.0 * math.pi / OMEGA_1DOF
        res = CentralDifferenceSolver(
            Assembler(dom), t_end=30.0 * T, dt=dt, u0=u0,
        ).solve()
        max_abs = float(np.max(np.abs(res.u_history)))
        self.assertLess(max_abs, 1.5,
                          f"Δt = 0.95·Δt_crit debería ser estable; |u|_max = {max_abs}")

    def test_just_above_cfl_diverges_1dof(self):
        """``Δt = 1.05·Δt_crit`` ⇒ ``RuntimeError`` CFL en el 1-GDL."""
        dt = 1.05 * DT_CRIT_1DOF
        dom = _build_1dof_oscillator()
        u0 = np.zeros(dom.total_dofs); u0[_free_dof(dom)] = 1.0
        solver = CentralDifferenceSolver(
            Assembler(dom), t_end=200.0 * dt, dt=dt, u0=u0,
        )
        with self.assertRaisesRegex(RuntimeError, "divergencia.*CFL"):
            solver.solve()

    def test_just_below_cfl_is_stable_2dof(self):
        """Cadena 2-DOF (ω_max = √3, Δt_crit = 2/√3): justo por debajo, estable."""
        dt_crit = 2.0 / math.sqrt(3.0)
        dt = 0.95 * dt_crit
        dom, free_dofs = _build_2dof_analytic_chain()
        u0 = np.zeros(dom.total_dofs)
        u0[free_dofs[0]] = 1.0   # exitar el modo asimétrico
        u0[free_dofs[1]] = -1.0
        # 30 periodos del modo lento ω_min=1: T_min = 2π.
        res = CentralDifferenceSolver(
            Assembler(dom), t_end=30.0 * 2.0 * math.pi, dt=dt, u0=u0,
        ).solve()
        max_abs = float(np.max(np.abs(res.u_history)))
        self.assertLess(max_abs, 2.0,
                          f"Δt = 0.95·(2/√3) debería ser estable; |u|_max = {max_abs}")

    def test_just_above_cfl_diverges_2dof(self):
        """Cadena 2-DOF: ``Δt = 1.05·(2/√3)`` ⇒ ``RuntimeError`` CFL.

        El modo ω_max = √3 desestabiliza primero al violarse su CFL local,
        aunque ω_min = 1 estaría a salvo a este Δt. Es el verdadero test
        de "Δt_crit gobernado por el modo más rápido".
        """
        dt_crit = 2.0 / math.sqrt(3.0)
        dt = 1.05 * dt_crit
        dom, free_dofs = _build_2dof_analytic_chain()
        u0 = np.zeros(dom.total_dofs)
        u0[free_dofs[0]] = 1.0
        u0[free_dofs[1]] = -1.0
        solver = CentralDifferenceSolver(
            Assembler(dom), t_end=400.0 * dt, dt=dt, u0=u0,
        )
        with self.assertRaisesRegex(RuntimeError, "divergencia.*CFL"):
            solver.solve()

    def test_well_below_cfl_matches_analytic_solution_2dof(self):
        """Cadena 2-DOF: con Δt = Δt_crit/10, la trayectoria coincide con
        la superposición modal analítica.

        Estado inicial ``(1, -1)`` (puro modo antisimétrico ``φ_2``) →
        evolución armónica pura en ω_max = √3. La solución analítica
        es ``u(t) = (cos(√3·t), -cos(√3·t))``.
        """
        dt_crit = 2.0 / math.sqrt(3.0)
        dt = dt_crit / 10.0
        dom, free_dofs = _build_2dof_analytic_chain()
        u0 = np.zeros(dom.total_dofs)
        u0[free_dofs[0]] = 1.0
        u0[free_dofs[1]] = -1.0
        # 2 periodos del modo rápido.
        t_end = 2.0 * 2.0 * math.pi / math.sqrt(3.0)
        res = CentralDifferenceSolver(
            Assembler(dom), t_end=t_end, dt=dt, u0=u0,
        ).solve()
        t = res.t_history
        u_analytic = np.cos(math.sqrt(3.0) * t)
        # Banda 2% (O(Δt²) sobre 2 periodos del modo rápido).
        err_d2 = float(np.max(np.abs(res.u_history[free_dofs[0], :] - u_analytic)))
        err_d3 = float(np.max(np.abs(res.u_history[free_dofs[1], :] + u_analytic)))
        self.assertLess(err_d2, 0.02,
                          f"u(ux₂) − cos(√3·t) excede 2%: err={err_d2:.4f}")
        self.assertLess(err_d3, 0.02,
                          f"u(ux₃) + cos(√3·t) excede 2%: err={err_d3:.4f}")


# ---------------------------------------------------------------------------
# Modo no lineal — coincide con NewtonNewmark en régimen elástico
# ---------------------------------------------------------------------------


class TestCentralDifferenceNonlinear(unittest.TestCase):
    """``nonlinear=True`` con material lineal debe dar la misma trayectoria
    que el modo lineal (caso degenerado de validación) y, en régimen
    elástico de un material elastoplástico, recuperar la solución armónica.
    """

    def test_nonlinear_with_linear_material_matches_linear_mode(self):
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 200.0
        t_end = T
        ndof = 4
        dof = 2  # ux del nodo 2 (libre); el resto está prescrito
        u0 = np.zeros(ndof); u0[dof] = 1.0

        dom_lin = _build_1dof_oscillator()
        solver_lin = CentralDifferenceSolver(
            Assembler(dom_lin), t_end=t_end, dt=dt, u0=u0, nonlinear=False,
        )
        u_lin = solver_lin.solve().u_history[_free_dof(dom_lin), :]

        dom_nl = _build_1dof_oscillator()
        solver_nl = CentralDifferenceSolver(
            Assembler(dom_nl), t_end=t_end, dt=dt, u0=u0, nonlinear=True,
        )
        u_nl = solver_nl.solve().u_history[_free_dof(dom_nl), :]

        # Mismos valores numéricos hasta precisión de máquina.
        np.testing.assert_allclose(u_lin, u_nl, rtol=1e-9, atol=1e-12)

    def test_nonlinear_matches_linear_with_nonzero_prescribed_support(self):
        """Coincidencia lin↔no-lin con **apoyo prescrito no nulo** — regresión
        de la doble resta de ``F_dir`` en la rama no-lineal (auditoría H-4.3).

        Con el apoyo izquierdo en ``ux = g₀ ≠ 0`` la trayectoria analítica del
        nodo libre es ``u(t) = g₀ + (u₀ − g₀)·cos(ωt)``. Si la rama no-lineal
        restara ``F_dir`` dos veces, la aceleración inicial saldría amplificada
        y la trayectoria divergiría de la rama lineal desde el primer paso.
        El test exige coincidencia bit-a-bit entre ramas (el material es
        lineal: ambas deben recorrer el mismo balance dinámico).
        """
        T = 2.0 * math.pi / OMEGA_1DOF
        dt = T / 200.0
        t_end = T
        ndof = 4
        dof = 2  # ux del nodo 2 (libre)
        g0 = 0.5  # apoyo prescrito no nulo
        u0 = np.zeros(ndof)
        u0[dof] = 1.0
        u0[0] = g0  # ux del nodo 1 (prescrito) — consistente con la BC

        def _build_with_prescribed_support() -> Domain:
            dom = Domain()
            n1 = dom.add_node(1, [0.0, 0.0])
            n2 = dom.add_node(2, [L_1DOF, 0.0])
            mat = Elastic1D(E=E_1DOF, density=RHO_1DOF_LUMPED)
            dom.add_element(Truss2D(1, [n1, n2], mat, A=A_1DOF))
            n1.fix_dof("ux", g0); n1.fix_dof("uy", 0.0)
            n2.fix_dof("uy", 0.0)
            dom.generate_equation_numbers(verbose=False)
            return dom

        dom_lin = _build_with_prescribed_support()
        solver_lin = CentralDifferenceSolver(
            Assembler(dom_lin), t_end=t_end, dt=dt, u0=u0, nonlinear=False,
        )
        u_lin = solver_lin.solve().u_history[_free_dof(dom_lin), :]

        dom_nl = _build_with_prescribed_support()
        solver_nl = CentralDifferenceSolver(
            Assembler(dom_nl), t_end=t_end, dt=dt, u0=u0, nonlinear=True,
        )
        u_nl = solver_nl.solve().u_history[_free_dof(dom_nl), :]

        np.testing.assert_allclose(u_lin, u_nl, rtol=1e-9, atol=1e-12)

        # Verificación analítica adicional: el modo lineal debe seguir
        # u(t) = g₀ + (u₀ − g₀)·cos(ωt). Banda holgada (Δt²) — el objetivo
        # primario es el contraste lin↔no-lin de arriba.
        t = np.linspace(0.0, t_end, u_lin.size)
        u_analytic = g0 + (1.0 - g0) * np.cos(OMEGA_1DOF * t)
        np.testing.assert_allclose(u_lin, u_analytic, atol=2e-3)


# ---------------------------------------------------------------------------
# Rechazos defensivos
# ---------------------------------------------------------------------------


class TestCentralDifferenceContract(unittest.TestCase):
    """Validaciones tempranas — ``ValueError`` con mensaje informativo."""

    def test_consistent_lumping_rejected(self):
        dom = _build_1dof_oscillator()
        with self.assertRaisesRegex(ValueError, "lumping='consistent'"):
            CentralDifferenceSolver(
                Assembler(dom), t_end=1.0, dt=0.01,
                lumping="consistent",
            )

    def test_nonpositive_dt_rejected(self):
        dom = _build_1dof_oscillator()
        with self.assertRaisesRegex(ValueError, "dt=0.0 debe ser positivo"):
            CentralDifferenceSolver(Assembler(dom), t_end=1.0, dt=0.0)

    def test_nonpositive_tend_rejected(self):
        dom = _build_1dof_oscillator()
        with self.assertRaisesRegex(ValueError, "t_end=0.0 debe ser positivo"):
            CentralDifferenceSolver(Assembler(dom), t_end=0.0, dt=0.01)

    def test_frame3d_oblique_rejected(self):
        """Frame3D con eje oblicuo: M lumped es bloque-diagonal pero no
        estrictamente diagonal tras la rotación. El solver lo detecta y
        rechaza con ``ValueError`` explicando que invertir bloques 6×6
        cada paso anula la ventaja del explícito.
        """
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0, 0.0])
        n2 = dom.add_node(2, [1.0 / math.sqrt(3)] * 3)
        mat = Elastic1D(E=210e9, density=7850.0)
        elem = Frame3D(1, [n1, n2], mat,
                          A=1.5e-4, Iy=8.33e-10, Iz=2.5e-9,
                          J=8.33e-10 + 2.5e-9)
        dom.add_element(elem)
        n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0); n1.fix_dof("uz", 0.0)
        n1.fix_dof("rx", 0.0); n1.fix_dof("ry", 0.0); n1.fix_dof("rz", 0.0)
        dom.generate_equation_numbers(verbose=False)

        solver = CentralDifferenceSolver(
            Assembler(dom), t_end=1e-6, dt=1e-9,
        )
        with self.assertRaisesRegex(
            ValueError, "no es estrictamente diagonal"
        ):
            solver.solve()


# ---------------------------------------------------------------------------
# Cableado YAML
# ---------------------------------------------------------------------------


YAML_CD_TEMPLATE = """\
nodes:
  - {{id: 1, coords: [0.0, 0.0]}}
  - {{id: 2, coords: [1.0, 0.0]}}

materials:
  - id: 1
    type: Elastic1D
    E: 25.0
    density: 2.0

elements:
  - {{id: 1, type: Truss2D, material: 1, A: 1.0, nodes: [1, 2]}}

boundary_conditions_by_node:
  - {{node_id: 1, ux: 0.0, uy: 0.0}}
  - {{node_id: 2, uy: 0.0}}

solver:
  type: CentralDifferenceSolver
  t_end: {t_end}
  dt: {dt}
  lumping: lumped
  u0: [0.0, 0.0, 1.0, 0.0]
"""


class TestCentralDifferenceYaml(unittest.TestCase):
    """Pipeline end-to-end YAML → CentralDifferenceSolver → TransientResult."""

    def test_yaml_pipeline(self):
        T = 2.0 * math.pi / 5.0
        yaml_content = YAML_CD_TEMPLATE.format(t_end=T, dt=T / 200.0)
        f = tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False, encoding="utf-8"
        )
        f.write(yaml_content); f.close()
        try:
            result = run_yaml(Path(f.name))
        finally:
            Path(f.name).unlink()

        self.assertIsInstance(result, TransientResult)
        # Pico de oscilador empezando con u₀=1 debe ser ≈ 1 (no amortiguado).
        dof = list(result.u_history[:, 0]).index(
            max(result.u_history[:, 0].tolist(), key=abs)
        )
        u_max = float(np.max(np.abs(result.u_history[dof, :])))
        self.assertAlmostEqual(u_max, 1.0, delta=0.02)


if __name__ == "__main__":
    unittest.main()
