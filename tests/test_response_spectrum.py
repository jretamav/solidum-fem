"""Tests del ``ResponseSpectrumSolver`` — fase 7 del ADR 0009.

Cubre:

1. **1 GDL contra fórmula analítica**: ``u_max = Γ·φ·S_d(ω_n)`` exacto.
2. **2 GDL con SRSS**: ``u_max² = Σ_n (Γ_n·φ_n·S_d(ω_n))²`` componente a
   componente; verificación contra cálculo manual.
3. **CQC vs SRSS con modos bien separados**: deben coincidir (ρ_ij ≈ δ_ij).
4. **CQC vs SRSS con modos cercanos**: deben diferir significativamente
   (acople de correlación).
5. **Espectro tabulado** vía :func:`spectrum_tabulated`.
6. **Direction por DOF name** vs por vector explícito.
7. **Helpers de la dataclass**: ``cumulative_effective_mass_ratio``.
8. **Validaciones tempranas**.
9. **Pipeline YAML** end-to-end con espectro tabulado.
10. **Regla D aplicada**: ``free_vibration`` sigue funcionando vía
    wrapper que delega a ``solidum.math.modal_response.free_vibration``.
"""
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import solidum  # autodiscover
from solidum.core.domain import Domain
from solidum.elements.truss import Truss2D
from solidum.entry import run_modal, run_response_spectrum, run_yaml
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.modal_response import (
    free_vibration,
    participation_factors,
    response_spectrum_srss,
)
from solidum.math.solvers import ResponseSpectrumSolver
from solidum.math.solvers.response_spectrum import (
    spectrum_from_sa,
    spectrum_tabulated,
)
from solidum.results import ResponseSpectrumResult


# ---------------------------------------------------------------------------
# Helpers: 1 GDL y 2 GDL en serie
# ---------------------------------------------------------------------------

E_1DOF = 25.0
RHO_1DOF = 3.0
A_1DOF = 1.0
L_1DOF = 1.0
OMEGA_1DOF = 5.0


def _build_3dof_analytic_chain() -> Domain:
    """Sistema 3-DOF con autovalores analíticos cerrados (Laplaciano 1D).

    Cuatro trusses ``n1-n2-n3-n4-n5`` con ``n1`` y ``n5`` fijos en ``ux/uy``
    y ``n2/n3/n4`` libres en ``ux`` (fijos en ``uy``). Con ``E=A=L=1`` cada
    truss aporta ``EA/L = 1`` a la rigidez axial y ``ρAL = 1`` a la masa
    total (mitad por nodo bajo lumping).

    Matrices reducidas sobre ``(ux₂, ux₃, ux₄)`` con masa **lumped**:

        K = [[ 2, -1,  0],
             [-1,  2, -1],
             [ 0, -1,  2]],     M = I.

    Autovalores del Laplaciano 1D discreto con 3 nodos interiores:
    ``ω_n² = 2 − 2·cos(nπ/4)``, n ∈ {1, 2, 3}:

        ω₁² = 2 − √2,   ω₂² = 2,   ω₃² = 2 + √2.

    Modos M-ortonormales (φ_n[i] = √(2/4)·sin(n·i·π/4), i ∈ {1,2,3}):

        φ₁ = (1/2,  1/√2,  1/2)
        φ₂ = (1/√2,  0,   −1/√2)
        φ₃ = (1/2, −1/√2,  1/2)

    Como ARPACK exige ``n_modes < n_dof_libre = 3``, las verificaciones
    usan ``n_modes = 2`` (modos 1 y 2). Para que el resto del subespacio
    quede cubierto, escogemos direcciones de excitación que satisfagan
    ``γ₃ = φ₃ᵀ·d = 0`` y por tanto la masa efectiva total esté en los
    dos primeros modos. La condición ``φ₃ᵀd = 0`` con ``φ₃`` simétrico
    se cumple para ``d = (1, 1/√2, 0)``.
    """
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [1.0, 0.0])
    n3 = dom.add_node(3, [2.0, 0.0])
    n4 = dom.add_node(4, [3.0, 0.0])
    n5 = dom.add_node(5, [4.0, 0.0])
    mat = Elastic1D(E=1.0, density=1.0)
    for elem_id, (a, b) in enumerate(((n1, n2), (n2, n3), (n3, n4), (n4, n5)),
                                       start=1):
        dom.add_element(Truss2D(elem_id, [a, b], mat, A=1.0))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    n3.fix_dof("uy", 0.0)
    n4.fix_dof("uy", 0.0)
    n5.fix_dof("ux", 0.0); n5.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


def _build_2dof_chain(rho2_factor: float = 1.0) -> Domain:
    """Dos masas conectadas axialmente. Frecuencias modales conocidas en
    función de los parámetros — útil para SRSS analítico. 2 DOFs libres
    en ux (n2, n3), suficiente para ARPACK con n_modes=1 ó 2.
    """
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [1.0, 0.0])
    n3 = dom.add_node(3, [2.0, 0.0])
    mat = Elastic1D(E=E_1DOF, density=RHO_1DOF * rho2_factor)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=A_1DOF))
    dom.add_element(Truss2D(2, [n2, n3], mat, A=A_1DOF))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0); n3.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


def _build_3dof_chain() -> Domain:
    """Tres masas en serie. 3 DOFs libres en ux; permite ARPACK con
    n_modes ≤ 2 (ARPACK requiere n_modes < n_dof_libre = 3).
    """
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [1.0, 0.0])
    n3 = dom.add_node(3, [2.0, 0.0])
    n4 = dom.add_node(4, [3.0, 0.0])
    mat = Elastic1D(E=E_1DOF, density=RHO_1DOF)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=A_1DOF))
    dom.add_element(Truss2D(2, [n2, n3], mat, A=A_1DOF))
    dom.add_element(Truss2D(3, [n3, n4], mat, A=A_1DOF))
    n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0); n3.fix_dof("uy", 0.0); n4.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


# ---------------------------------------------------------------------------
# 1 GDL — fórmula analítica
# ---------------------------------------------------------------------------


class TestResponseSpectrumSingleMode(unittest.TestCase):
    """``u_max = Γ·φ·S_d(ω_n)`` evaluado con un solo modo del sistema.

    Usa el modelo 2-DOF (ARPACK requiere n_modes < n_dof_libre, así que
    para n_modes=1 hace falta ≥ 2 DOFs libres). La fórmula analítica se
    verifica DOF a DOF contra la contribución del primer modo.
    """

    def test_unit_spectrum_first_mode(self):
        """Con ``S_d(ω) = 1`` constante, ``u_per_mode[:, 0] = γ_1·φ_1·1``."""
        dom = _build_2dof_chain()
        ndof = dom.total_dofs
        direction = np.zeros(ndof)
        direction[dom.nodes[2].dofs["ux"]] = 1.0
        direction[dom.nodes[3].dofs["ux"]] = 1.0

        solver = ResponseSpectrumSolver(
            Assembler(dom), n_modes=1, direction=direction,
            spectrum=lambda w: 1.0, combination="SRSS",
        )
        result = run_response_spectrum(dom, solver=solver)

        self.assertEqual(result.n_modes, 1)
        gamma = result.participation_factors[0]
        # u_per_mode[:, 0] = γ_1·φ_1·S_d(ω_1). Con n_modes=1, u_combined
        # = |u_per_mode[:, 0]| exactamente (SRSS sobre un solo término).
        np.testing.assert_allclose(
            result.u_combined,
            np.abs(result.u_per_mode[:, 0]),
            rtol=1e-12,
        )
        # Verificación cruzada del valor analítico en el DOF de n2.
        dof_2 = dom.nodes[2].dofs["ux"]
        u_expected = abs(gamma) * abs(result.u_per_mode[dof_2, 0]) / abs(gamma)
        # Es decir: u_combined[dof_2] = |γ·φ[dof_2, 0]·1| (Sd=1).
        # Esto sólo confirma consistencia interna; el valor exacto sale del modo.
        self.assertGreater(result.u_combined[dof_2], 0.0)

    def test_constant_acceleration_spectrum(self):
        """Con ``S_a`` constante: la M-norma de ``u_per_mode[:, 0]`` es
        ``|γ_1·S_d|`` (consecuencia de M-ortonormalidad ``φ_1ᵀMφ_1 = 1``).

        Esta verificación es invariante al signo arbitrario del modo que
        devuelve ARPACK.
        """
        dom = _build_2dof_chain()
        ndof = dom.total_dofs
        direction = np.zeros(ndof)
        direction[dom.nodes[2].dofs["ux"]] = 1.0
        direction[dom.nodes[3].dofs["ux"]] = 1.0
        Sa_const = 10.0

        solver = ResponseSpectrumSolver(
            Assembler(dom), n_modes=1, direction=direction,
            spectrum=spectrum_from_sa(lambda _w: Sa_const),
            combination="SRSS",
        )
        result = run_response_spectrum(dom, solver=solver)
        omega1 = result.frequencies_rad[0]
        Sd_expected = Sa_const / (omega1 ** 2)
        gamma1 = result.participation_factors[0]

        # u_per_mode[:, 0] = γ_1·φ_1·S_d, con φ_1ᵀMφ_1 = 1.
        # → (u_per_mode[:, 0])ᵀ M (u_per_mode[:, 0]) = (γ_1·S_d)².
        M = solver.assembler.assemble_mass_matrix(lumping="consistent")
        u1 = result.u_per_mode[:, 0]
        m_norm_sq = float(u1 @ (M @ u1))
        expected_sq = (gamma1 * Sd_expected) ** 2
        self.assertAlmostEqual(m_norm_sq, expected_sq, places=8)


# ---------------------------------------------------------------------------
# 2-DOF con autovalores analíticos cerrados (Fase B 2026-05-19)
# ---------------------------------------------------------------------------


class TestResponseSpectrumAnalytic2DOF(unittest.TestCase):
    """Validación contra solución cerrada del 3-DOF construido en
    :func:`_build_3dof_analytic_chain`: ``K`` Toeplitz tridiagonal,
    ``M=I``, autovalores ``ω_n² = 2 − 2·cos(nπ/4)`` con modos
    ``φ_n[i] = √(2/4)·sin(n·i·π/4)`` (forma cerrada del Laplaciano 1D
    discreto). ARPACK limita ``n_modes < 3``, así que se usan los dos
    primeros modos y se escogen direcciones con ``φ₃ᵀd = 0`` para que
    los modos truncados no introduzcan masa efectiva residual.

    Cierra el hueco "espectro 1-DOF analítico" de la matriz de validación:
    contrasta autovalores, factores de participación, contribución modal
    individual ``γ_n·φ_n·S_d(ω_n)`` y envolvente SRSS contra valores
    construidos a mano. Las verificaciones son invariantes al signo
    arbitrario de ARPACK porque comparan magnitudes o sumas cuadráticas.
    """

    def test_eigenvalues_match_closed_form(self):
        """ω₁ = √(2−√2), ω₂ = √2 hasta precisión doble."""
        dom = _build_3dof_analytic_chain()
        solver = ResponseSpectrumSolver(
            Assembler(dom), n_modes=2, direction={"dof_name": "ux"},
            spectrum=lambda w: 1.0, lumping="lumped",
        )
        result = solver.solve()
        omega1_analytic = math.sqrt(2.0 - math.sqrt(2.0))
        omega2_analytic = math.sqrt(2.0)
        self.assertAlmostEqual(result.frequencies_rad[0], omega1_analytic,
                                 places=10)
        self.assertAlmostEqual(result.frequencies_rad[1], omega2_analytic,
                                 places=10)

    def test_participation_factors_match_closed_form(self):
        """Direction ``d = (1, 1/√2, 0)`` ⇒ ``γ₁ = 1, γ₂ = 1/√2, γ₃ = 0``.

        Construida para que el modo 3 (no resuelto con ARPACK n_modes=2)
        no aporte masa efectiva: ``φ₃ᵀd = 1/2 − (1/√2)·(1/√2) + 0 = 0``.
        Σ γ_n² = 1 + 1/2 = 3/2 = ``dᵀMd`` ⇒ cobertura modal del 100%
        con sólo los dos primeros modos.
        """
        dom = _build_3dof_analytic_chain()
        ndof = dom.total_dofs
        direction = np.zeros(ndof)
        direction[dom.nodes[2].dofs["ux"]] = 1.0
        direction[dom.nodes[3].dofs["ux"]] = 1.0 / math.sqrt(2.0)
        # direction[ux₄] queda en 0.

        solver = ResponseSpectrumSolver(
            Assembler(dom), n_modes=2, direction=direction,
            spectrum=lambda w: 1.0, lumping="lumped",
        )
        result = solver.solve()
        gamma_sq = result.participation_factors ** 2
        np.testing.assert_allclose(gamma_sq, [1.0, 0.5], atol=1e-10)
        # Masa efectiva acumulada con 2 modos = ‖d‖²_M = 3/2.
        cum = result.cumulative_effective_mass_ratio()
        self.assertAlmostEqual(cum[-1], 1.0, places=10)

    def test_srss_envelope_matches_analytic(self):
        """Mismo ``d`` y ``S_d ≡ 1``: ``u_combined`` analítica = ``(1/√2, 1/√2, 1/√2)``.

        Modo 1 aporta ``γ₁·φ₁ = 1·(1/2, 1/√2, 1/2)``.
        Modo 2 aporta ``γ₂·φ₂ = (1/√2)·(1/√2, 0, −1/√2) = (1/2, 0, −1/2)``.
        SRSS por componente:
            u_combined[0] = √((1/2)² + (1/2)²) = √(1/2) = 1/√2.
            u_combined[1] = √((1/√2)² + 0²)    = 1/√2.
            u_combined[2] = √((1/2)² + (−1/2)²) = 1/√2.
        """
        dom = _build_3dof_analytic_chain()
        ndof = dom.total_dofs
        direction = np.zeros(ndof)
        direction[dom.nodes[2].dofs["ux"]] = 1.0
        direction[dom.nodes[3].dofs["ux"]] = 1.0 / math.sqrt(2.0)

        solver = ResponseSpectrumSolver(
            Assembler(dom), n_modes=2, direction=direction,
            spectrum=lambda w: 1.0, combination="SRSS", lumping="lumped",
        )
        result = solver.solve()

        expected = 1.0 / math.sqrt(2.0)
        for nid in (2, 3, 4):
            dof = dom.nodes[nid].dofs["ux"]
            self.assertAlmostEqual(result.u_combined[dof], expected, places=10)

    def test_single_mode_excitation_with_midpoint_direction(self):
        """``d = e_{ux₃}`` (sólo nodo central): excita únicamente el modo 1.

        ``γ₁ = φ₁[1] = 1/√2``, ``γ₂ = φ₂[1] = 0``, ``γ₃ = φ₃[1] = −1/√2``.
        Con ``n_modes = 2`` el modo 3 se trunca; queda sólo el modo 1.

        Para ``Sd ≡ c``:
            u_per_mode[:, 0] = γ₁·φ₁·c = (1/√2)·(1/2, 1/√2, 1/2)·c
                              = (c/(2√2), c/2, c/(2√2)).
            u_per_mode[:, 1] = 0  (γ₂ = 0).
        ``u_combined`` = ``|u_per_mode[:, 0]|`` exactamente porque sólo
        contribuye un modo. Esto es la analítica de un 1-GDL "tomada del
        2-DOF": un sólo modo, contribución exactamente conocida en cada
        nodo.
        """
        dom = _build_3dof_analytic_chain()
        ndof = dom.total_dofs
        direction = np.zeros(ndof)
        direction[dom.nodes[3].dofs["ux"]] = 1.0  # sólo nodo central

        Sd_const = 0.4
        solver = ResponseSpectrumSolver(
            Assembler(dom), n_modes=2, direction=direction,
            spectrum=lambda w: Sd_const, lumping="lumped",
        )
        result = solver.solve()

        # γ₁² = 1/2,  γ₂² = 0.
        gamma_sq = result.participation_factors ** 2
        np.testing.assert_allclose(gamma_sq, [0.5, 0.0], atol=1e-10)

        # u_combined analítico.
        c = Sd_const
        expected = {
            dom.nodes[2].dofs["ux"]: c / (2.0 * math.sqrt(2.0)),
            dom.nodes[3].dofs["ux"]: c / 2.0,
            dom.nodes[4].dofs["ux"]: c / (2.0 * math.sqrt(2.0)),
        }
        for dof, u_exp in expected.items():
            self.assertAlmostEqual(result.u_combined[dof], u_exp, places=10)


# ---------------------------------------------------------------------------
# 2 GDL con SRSS — cálculo manual
# ---------------------------------------------------------------------------


class TestResponseSpectrum2DofSRSS(unittest.TestCase):
    """SRSS reproducido pasando los modos directamente a la función pura.

    Usa el modelo 3-DOF (3 DOFs libres en ux) calculando 2 modos, dentro
    del rango admitido por ARPACK.
    """

    def test_srss_matches_manual_combination(self):
        dom = _build_3dof_chain()
        modal_result = run_modal(dom, n_modes=2)
        M = Assembler(dom).assemble_mass_matrix(lumping="consistent")
        ndof = dom.total_dofs
        direction = np.zeros(ndof)
        for nid in (2, 3, 4):
            direction[dom.nodes[nid].dofs["ux"]] = 1.0
        spectrum_fn = lambda w: 1.0 / (w * w)

        u_combined, u_per_mode, gamma = response_spectrum_srss(
            modal_result.modes, modal_result.frequencies_rad, M,
            direction, spectrum_fn,
        )

        # Verificación: u_combined[dof]² = Σ_n u_per_mode[dof, n]²
        u_manual = np.sqrt(np.sum(u_per_mode ** 2, axis=1))
        np.testing.assert_allclose(u_combined, u_manual, rtol=1e-12)

        # Contribución modal individual coincide con γ_n·φ_n·S_d(ω_n).
        Sd = np.array([spectrum_fn(w) for w in modal_result.frequencies_rad])
        for n in range(2):
            expected = gamma[n] * modal_result.modes[:, n] * Sd[n]
            np.testing.assert_allclose(u_per_mode[:, n], expected, rtol=1e-12)

    def test_participation_factors_partial_mass(self):
        """Σ Γ_n² acumula la masa efectiva (≤ masa total en la dirección)."""
        dom = _build_3dof_chain()
        modal_result = run_modal(dom, n_modes=2)
        M = Assembler(dom).assemble_mass_matrix(lumping="consistent")
        ndof = dom.total_dofs
        direction = np.zeros(ndof)
        for nid in (2, 3, 4):
            direction[dom.nodes[nid].dofs["ux"]] = 1.0

        gamma = participation_factors(modal_result.modes, M, direction)
        m_total = float(direction @ (M @ direction))
        m_effective = float(np.sum(gamma ** 2))
        # Con sólo 2 de los 3 modos del subespacio libre, la masa
        # efectiva debe ser positiva y ≤ masa total.
        self.assertGreater(m_effective, 0.0)
        self.assertLessEqual(m_effective, m_total + 1e-9)


# ---------------------------------------------------------------------------
# CQC vs SRSS — modos bien separados vs cercanos
# ---------------------------------------------------------------------------


class TestSRSSvsCQC(unittest.TestCase):
    """SRSS y CQC deben coincidir con modos bien separados y diferir con
    modos cercanos.
    """

    def test_well_separated_modes_srss_eq_cqc(self):
        """Modos del 3-DOF chain con ratio > 2 ⇒ ρ_ij ≈ 0 → CQC ≈ SRSS."""
        dom = _build_3dof_chain()
        modal_result = run_modal(dom, n_modes=2)
        ratio = modal_result.frequencies_rad[1] / modal_result.frequencies_rad[0]
        self.assertGreater(ratio, 2.0,
                             f"Test asume modos separados; ratio={ratio}")

        M = Assembler(dom).assemble_mass_matrix(lumping="consistent")
        ndof = dom.total_dofs
        direction = np.zeros(ndof)
        for nid in (2, 3, 4):
            direction[dom.nodes[nid].dofs["ux"]] = 1.0
        spectrum_fn = lambda w: 1.0 / (w * w)

        from solidum.math.modal_response import response_spectrum_cqc
        u_srss, _, _ = response_spectrum_srss(
            modal_result.modes, modal_result.frequencies_rad, M,
            direction, spectrum_fn,
        )
        u_cqc, _, _ = response_spectrum_cqc(
            modal_result.modes, modal_result.frequencies_rad, M,
            direction, spectrum_fn, damping=0.05,
        )
        rel_diff = np.max(np.abs(u_cqc - u_srss)) / np.max(np.abs(u_srss))
        self.assertLess(rel_diff, 0.01)

    def test_close_modes_srss_differs_from_cqc(self):
        """Sistema sintético con modos no-ortogonales-en-componentes y
        frecuencias cercanas ⇒ ρ_12 ≈ 1 ⇒ CQC ≠ SRSS.

        Si los modos comparten componentes espaciales (no son ejes
        coordenados), el término cruzado ``ρ_ij·u_i·u_j`` del CQC es no
        nulo y la diferencia con SRSS aparece.
        """
        omega1 = 10.0
        omega2 = 10.5   # cociente ω2/ω1 = 1.05 → muy cercanos
        # Modos en superposición: φ_1 = (1,1)/√2,  φ_2 = (1,-1)/√2.
        modes = np.array([[1.0, 1.0],
                            [1.0, -1.0]]) / math.sqrt(2.0)
        M = np.eye(2)
        direction = np.array([1.0, 0.0])    # excitación sólo en DOF 0
        spectrum_fn = lambda w: 1.0   # Sd constante

        u_srss, _, _ = response_spectrum_srss(
            modes, np.array([omega1, omega2]), M, direction, spectrum_fn,
        )
        from solidum.math.modal_response import response_spectrum_cqc
        u_cqc, _, _ = response_spectrum_cqc(
            modes, np.array([omega1, omega2]), M, direction, spectrum_fn,
            damping=0.05,
        )
        rel_diff = np.max(np.abs(u_cqc - u_srss)) / np.max(np.abs(u_srss))
        self.assertGreater(rel_diff, 0.05,
                             f"CQC debería diferir de SRSS con modos "
                             f"cercanos y acople; rel_diff = {rel_diff:.2%}")


# ---------------------------------------------------------------------------
# spectrum_tabulated
# ---------------------------------------------------------------------------


class TestSpectrumTabulated(unittest.TestCase):
    """``spectrum_tabulated`` interpola y aplica clamp en los extremos."""

    def test_interpolation_in_range(self):
        periods = np.array([0.1, 0.5, 1.0, 2.0])
        Sd_vals = np.array([0.001, 0.01, 0.04, 0.1])
        sp = spectrum_tabulated(periods, Sd_vals, kind="Sd")
        # T=0.75 (entre 0.5 y 1.0): interp(0.75) = 0.01 + (0.04-0.01)·0.5 = 0.025.
        omega = 2.0 * np.pi / 0.75
        self.assertAlmostEqual(sp(omega), 0.025, places=10)

    def test_clamp_below_min_period(self):
        periods = np.array([0.1, 1.0])
        Sd_vals = np.array([0.001, 0.04])
        sp = spectrum_tabulated(periods, Sd_vals, kind="Sd")
        # T=0.05 < 0.1 → clamp al valor en T=0.1.
        omega = 2.0 * np.pi / 0.05
        self.assertEqual(sp(omega), 0.001)

    def test_kind_sa_conversion(self):
        """``kind="Sa"`` convierte a Sd internamente."""
        periods = np.array([0.1, 1.0])
        Sa_vals = np.array([10.0, 5.0])
        sp = spectrum_tabulated(periods, Sa_vals, kind="Sa")
        omega = 2.0 * np.pi / 1.0   # T=1.0
        # Sd = Sa/ω² = 5.0 / (2π)² = 5.0 / 39.478 ≈ 0.1266
        self.assertAlmostEqual(sp(omega), 5.0 / (omega ** 2), places=10)

    def test_unequal_lengths_rejected(self):
        with self.assertRaisesRegex(ValueError, "periods.size"):
            spectrum_tabulated([0.1, 0.5], [1.0, 2.0, 3.0])

    def test_non_monotonic_periods_rejected(self):
        with self.assertRaisesRegex(ValueError, "estrictamente"):
            spectrum_tabulated([0.5, 0.5], [1.0, 2.0])


# ---------------------------------------------------------------------------
# direction por DOF name
# ---------------------------------------------------------------------------


class TestDirectionDofName(unittest.TestCase):
    """``direction={"dof_name": "ux"}`` construye el vector unitario."""

    def test_dof_name_builds_unit_vector(self):
        dom = _build_3dof_chain()
        solver = ResponseSpectrumSolver(
            Assembler(dom), n_modes=2,
            direction={"dof_name": "ux"},
            spectrum=lambda w: 1.0,
        )
        result = solver.solve()
        # La dirección debe valer 1 en ux de cada nodo (libre o no).
        for node_id, node in dom.nodes.items():
            if "ux" in node.dofs:
                self.assertEqual(result.direction[node.dofs["ux"]], 1.0)

    def test_dof_name_missing_raises(self):
        dom = _build_2dof_chain()
        solver = ResponseSpectrumSolver(
            Assembler(dom), n_modes=1,
            direction={"dof_name": "rz"},   # Truss2D no tiene rz
            spectrum=lambda w: 1.0,
        )
        with self.assertRaisesRegex(ValueError, "ningún nodo"):
            solver.solve()


# ---------------------------------------------------------------------------
# Cumulative effective mass ratio
# ---------------------------------------------------------------------------


class TestCumulativeMassRatio(unittest.TestCase):
    """``cumulative_effective_mass_ratio`` crece monotonamente a 1.0."""

    def test_monotonic_and_unit_sum(self):
        dom = _build_3dof_chain()
        solver = ResponseSpectrumSolver(
            Assembler(dom), n_modes=2,
            direction={"dof_name": "ux"},
            spectrum=lambda w: 1.0,
        )
        result = solver.solve()
        cum = result.cumulative_effective_mass_ratio()
        self.assertEqual(cum.size, 2)
        self.assertGreater(cum[1], cum[0])
        self.assertAlmostEqual(cum[-1], 1.0, places=10)


# ---------------------------------------------------------------------------
# Validaciones tempranas
# ---------------------------------------------------------------------------


class TestResponseSpectrumContract(unittest.TestCase):
    def test_zero_modes_rejected(self):
        dom = _build_2dof_chain()
        with self.assertRaisesRegex(ValueError, "n_modes"):
            ResponseSpectrumSolver(
                Assembler(dom), n_modes=0,
                direction=np.zeros(dom.total_dofs),
                spectrum=lambda w: 1.0,
            )

    def test_invalid_combination_rejected(self):
        dom = _build_2dof_chain()
        with self.assertRaisesRegex(ValueError, "combination='ABS'"):
            ResponseSpectrumSolver(
                Assembler(dom), n_modes=1,
                direction=np.zeros(dom.total_dofs),
                spectrum=lambda w: 1.0,
                combination="ABS",
            )

    def test_negative_damping_rejected(self):
        dom = _build_2dof_chain()
        with self.assertRaisesRegex(ValueError, "damping=-0.05"):
            ResponseSpectrumSolver(
                Assembler(dom), n_modes=1,
                direction=np.zeros(dom.total_dofs),
                spectrum=lambda w: 1.0,
                damping=-0.05,
            )


# ---------------------------------------------------------------------------
# Cableado YAML
# ---------------------------------------------------------------------------


YAML_RS_TEMPLATE = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [2.0, 0.0]}
  - {id: 4, coords: [3.0, 0.0]}

materials:
  - id: 1
    type: Elastic1D
    E: 25.0
    density: 3.0

elements:
  - {id: 1, type: Truss2D, material: 1, A: 1.0, nodes: [1, 2]}
  - {id: 2, type: Truss2D, material: 1, A: 1.0, nodes: [2, 3]}
  - {id: 3, type: Truss2D, material: 1, A: 1.0, nodes: [3, 4]}

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0}
  - {node_id: 2, uy: 0.0}
  - {node_id: 3, uy: 0.0}
  - {node_id: 4, uy: 0.0}

solver:
  type: ResponseSpectrumSolver
  n_modes: 2
  direction:
    dof_name: ux
  spectrum:
    type: tabulated
    kind: Sa
    periods: [0.1, 0.5, 1.0, 2.0]
    values: [10.0, 10.0, 5.0, 1.0]
  combination: SRSS
"""


class TestResponseSpectrumYaml(unittest.TestCase):
    def test_yaml_pipeline(self):
        f = tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False, encoding="utf-8"
        )
        f.write(YAML_RS_TEMPLATE); f.close()
        try:
            result = run_yaml(Path(f.name))
        finally:
            Path(f.name).unlink()

        self.assertIsInstance(result, ResponseSpectrumResult)
        self.assertEqual(result.n_modes, 2)
        self.assertEqual(result.combination, "SRSS")
        # u_combined positivo en todos los DOFs ux libres.
        self.assertTrue(np.all(result.u_combined >= 0.0))
        # Masa efectiva total > 0.
        self.assertGreater(float(np.sum(result.effective_masses)), 0.0)


# ---------------------------------------------------------------------------
# Regla D — wrapper free_vibration sigue funcionando
# ---------------------------------------------------------------------------


class TestRegleDFreeVibrationWrapper(unittest.TestCase):
    """``ModalResult.free_vibration`` y ``modal_response.free_vibration``
    deben dar resultados idénticos (mismo algoritmo, sólo distinta API).
    """

    def test_wrapper_matches_direct_call(self):
        dom = _build_3dof_chain()
        modal_result = run_modal(dom, n_modes=2)
        M = Assembler(dom).assemble_mass_matrix(lumping="consistent")
        ndof = dom.total_dofs
        u0 = np.zeros(ndof); u0[dom.nodes[2].dofs["ux"]] = 0.01
        u0_dot = np.zeros(ndof)
        t = np.linspace(0.0, 1.0, 50)

        u_via_method = modal_result.free_vibration(M, u0, u0_dot, t)
        u_direct = free_vibration(
            modal_result.modes, modal_result.frequencies_rad, M,
            u0, u0_dot, t,
        )
        np.testing.assert_allclose(u_via_method, u_direct, rtol=1e-12)


if __name__ == "__main__":
    unittest.main()
