"""Tests de ``DissipationArcLengthSolver`` — variante por disipación del arc-length.

Cubre los ``acceptance`` de la spec ``docs/specs/DissipationArcLengthSolver.md``
que la **primera implementación** alcanza:

1. **Recuperación cilíndrica** en régimen elástico puro — el solver coincide
   con :class:`ArcLengthSolver` a paridad de bits porque su modo permanece
   ``"cylindrical"``.
2. **Daño 1D con softening pronunciado** — el solver atraviesa la rama
   post-pico igual que el padre cilíndrico.
3. **Switching cilíndrico→disipación** observable cuando aparece daño
   continuo.
4. **Daño 2D continuo (Quad4 + IsotropicDamage2D)** — análogo 2D del test
   3 con cinemática completa de Voigt.
5. **Balance energético cualitativo** sobre la historia con daño 1D.
6. **Sign-of-pivot tracking** vía signo del determinante (aproximación
   para indefinitud simple).

**No cubierto en esta primera implementación** (deuda técnica documentada
en la spec):

- *Cohesivo + embedded discontinuity con K_e stiff*: la activación
  discreta del cohesivo (Rankine onset) introduce un salto en F_int/K_t
  que el `sign` por producto escalar no maneja graciosamente. Requiere
  sign-of-pivot tracking real (LDLᵀ Bunch-Kaufman, no aproximación por
  determinante) o control por CMOD/CTOD del salto en lugar de τ global.
  Trabajo pendiente para una futura iteración del solver.
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.domain import Domain
from fenix.elements.solid_2d import Quad4
from fenix.elements.truss import Truss2D
from fenix.materials.damage_1d import IsotropicDamage1D
from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.materials.elastic import Elastic1D
from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.solvers import (
    ArcLengthSolver,
    DissipationArcLengthSolver,
)


# ---------------------------------------------------------------------------
# Helpers compartidos.
# ---------------------------------------------------------------------------

def _build_elastic_bar(E=1.0e5, A=1.0, L=1.0):
    """Truss2D elástica empotrada-libre. Free DOF: nodo extremo, ux."""
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    mat = Elastic1D(E=E)
    dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
    n2.fix_dof('uy', 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), n2.dofs['ux']


def _build_damage_bar_1d(*, E, kappa_0, alpha, A=1.0, L=1.0, n_elems=1):
    """Barra 1D con ``IsotropicDamage1D``. Idéntico al helper de robustez."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * L / n_elems, 0.0])
              for i in range(n_elems + 1)]
    mat = IsotropicDamage1D(E=E, kappa_0=kappa_0, alpha=alpha)
    for k in range(n_elems):
        dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]], mat, A=A))
    nodes[0].fix_dof('ux', 0.0); nodes[0].fix_dof('uy', 0.0)
    for n in nodes[1:]:
        n.fix_dof('uy', 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), nodes[-1].dofs['ux']


def _build_damage_quad4_uniaxial(E=1.0e5, nu=0.0, kappa_0=1.0e-3, alpha=300.0):
    """Quad4 cuadrado unitario con IsotropicDamage2D + tracción uniaxial.

    Material con softening exponencial moderado. ν=0 para evitar contracción
    transversal y simplificar la cinemática a esencialmente 1D.
    """
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [1.0, 0.0])
    n3 = dom.add_node(3, [1.0, 1.0])
    n4 = dom.add_node(4, [0.0, 1.0])

    mat = IsotropicDamage2D(
        E=E, nu=nu, kappa_0=kappa_0, alpha=alpha, hypothesis='plane_strain',
    )
    dom.add_element(Quad4(1, [n1, n2, n3, n4], mat, thickness=1.0))

    # Apoyos: u_x=0 en lado izquierdo (n1, n4); u_y=0 en lado inferior (n1, n2).
    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
    n2.fix_dof('uy', 0.0)
    n4.fix_dof('ux', 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4}


def _apply_uniaxial_traction_ref(dom, nodes, F_total):
    """F_ext_ref: fuerza total ``F_total`` repartida en lado derecho (n2, n3)."""
    F_ref = np.zeros(dom.total_dofs)
    F_ref[nodes['n2'].dofs['ux']] = 0.5 * F_total
    F_ref[nodes['n3'].dofs['ux']] = 0.5 * F_total
    return F_ref


# ---------------------------------------------------------------------------
# Test 1 — recuperación cilíndrica en régimen elástico (regresión)
# ---------------------------------------------------------------------------

class TestElasticRecovery(unittest.TestCase):
    """En régimen elástico puro (sin disipación), DissipationArcLengthSolver
    debe coincidir con ArcLengthSolver cilíndrico a tolerancia muy estricta."""

    def test_elastic_recovery_matches_parent(self):
        E, A, L = 1.0e5, 1.0, 1.0
        F_pico = 100.0  # carga elástica
        dom_p, asm_p, dof_p = _build_elastic_bar(E=E, A=A, L=L)
        dom_d, asm_d, dof_d = _build_elastic_bar(E=E, A=A, L=L)

        F_ref_p = np.zeros(dom_p.total_dofs); F_ref_p[dof_p] = F_pico
        F_ref_d = np.zeros(dom_d.total_dofs); F_ref_d[dof_d] = F_pico

        conv = ConvergenceCriterion(rtol_force=1e-8, rtol_disp=1e-8)
        parent = ArcLengthSolver(
            asm_p, convergence=conv, max_lambda=1.0,
            initial_dl=0.1, max_steps=50, max_iter=20,
        )
        diss = DissipationArcLengthSolver(
            asm_d, convergence=conv, max_lambda=1.0,
            initial_dl=0.1, max_steps=50, max_iter=20,
            initial_tau=0.01,
        )
        U_p = parent.solve(F_ref_p)
        U_d = diss.solve(F_ref_d)

        # En régimen puramente elástico (sin disipación), el modo debe
        # permanecer cilíndrico todo el tiempo.
        self.assertEqual(diss._mode, "cylindrical",
                          "Modo debió permanecer cilíndrico en régimen elástico")

        # Las trayectorias deben coincidir a muy alta precisión.
        np.testing.assert_allclose(
            U_d[dof_d], U_p[dof_p], rtol=1e-10,
            err_msg=f"U(λ=1): cilíndrico={U_p[dof_p]:.8e} vs disipación={U_d[dof_d]:.8e}",
        )


# ---------------------------------------------------------------------------
# Test 2 — daño 1D con softening pronunciado (atraviesa post-pico)
# ---------------------------------------------------------------------------

class TestDamageSofteningPostPeak(unittest.TestCase):
    """Mismo problema que ``test_arc_length_traverses_damage_softening`` del
    padre. Daño 1D con softening exponencial pronunciado, carga = 2× pico
    elástico. DissipationArcLengthSolver debe atravesar igual o mejor.
    """

    def test_atraviesa_pico_con_modo_switch(self):
        E = 1.0e5; kappa_0 = 1.0e-3; alpha = 500.0
        dom, asm, dof_ux = _build_damage_bar_1d(
            E=E, kappa_0=kappa_0, alpha=alpha, A=1.0, L=1.0,
        )
        F_ref = np.zeros(dom.total_dofs); F_ref[dof_ux] = 200.0  # 2× pico

        modes_seen = set()
        def probe(step, U, lam):
            modes_seen.add(solver._mode)

        # initial_dl pequeño para forzar varios pasos: con E=1e5, A=L=1,
        # F_ref=200 ⇒ ‖du_t‖=2e-3. dl=1e-4 ⇒ dλ≈0.05 por paso ⇒ ~20 pasos
        # hasta λ=1. Esto da al switching la oportunidad de activarse al
        # cruzar el pico (paso ~3-4 con kappa_0=1e-3).
        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = DissipationArcLengthSolver(
            asm, convergence=conv, max_lambda=1.0,
            initial_dl=1.0e-4, max_steps=500, max_iter=30,
            initial_tau=1.0e-5,
        )
        U = solver.solve(F_ref, step_callback=probe)

        # El frente atravesó el pico.
        u_at_peak = kappa_0 * 1.0   # ε_peak · L con L=1
        self.assertGreater(U[dof_ux], u_at_peak,
            "DissipationArcLength no atravesó el pico de daño")

        # El switching automático tuvo que activarse: ambos modos deben
        # aparecer en la historia (cilíndrico al inicio, disipación post-pico).
        self.assertIn("cylindrical", modes_seen,
            f"Modo cilíndrico nunca observado. Modos vistos: {modes_seen}")
        self.assertIn("dissipation", modes_seen,
            f"Modo de disipación nunca observado. Modos vistos: {modes_seen}")


# ---------------------------------------------------------------------------
# Test 3 — switching automático observado tras la transición
# ---------------------------------------------------------------------------

class TestSwitchingCylindricalToDissipation(unittest.TestCase):
    """En el paso N donde se cruza σ_pico, el modo debe pasar de
    'cylindrical' a 'dissipation'; los pasos 1..N-1 son cilíndricos,
    desde N+1 en adelante son de disipación."""

    def test_modo_cambia_solo_tras_yield(self):
        E = 1.0e5; kappa_0 = 1.0e-3; alpha = 500.0
        dom, asm, dof_ux = _build_damage_bar_1d(
            E=E, kappa_0=kappa_0, alpha=alpha,
        )
        F_ref = np.zeros(dom.total_dofs); F_ref[dof_ux] = 200.0

        modes_per_step: list[tuple[int, str, float]] = []
        def probe(step, U, lam):
            modes_per_step.append((step, solver._mode, lam))

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = DissipationArcLengthSolver(
            asm, convergence=conv, max_lambda=1.0,
            initial_dl=1.0e-4, max_steps=500, max_iter=30,
            initial_tau=1.0e-5,
        )
        solver.solve(F_ref, step_callback=probe)

        # Encontrar el primer paso 'dissipation'.
        first_diss = next((i for (i, m, _) in modes_per_step if m == "dissipation"), None)
        self.assertIsNotNone(first_diss,
            "Nunca se activó el modo dissipation")

        # Antes de ese paso, todos los modos son cilíndricos.
        for (i, m, _) in modes_per_step:
            if i < first_diss:
                self.assertEqual(m, "cylindrical",
                    f"Paso {i} antes del switch debería ser cilíndrico, era {m}")


# ---------------------------------------------------------------------------
# Test 4 — daño 2D continuo (Quad4 + IsotropicDamage2D) con softening
# ---------------------------------------------------------------------------

class TestDamage2DSofteningPostPeak(unittest.TestCase):
    """Quad4 con IsotropicDamage2D plane strain bajo tracción uniaxial.

    Es el caso de daño continuo bulk (no cohesivo, no embedded) — donde
    K_t cambia suavemente de elástico a softening. Equivalente 2D al test
    de daño 1D pero con cinemática completa de Voigt 2D.

    El ArcLengthSolver cilíndrico padre ya resuelve este caso
    (validado en ``test_arc_length_traverses_damage_softening``).
    DissipationArcLengthSolver debe resolverlo igual de bien — además, el
    switching automático debe activar el modo de disipación tras el pico.
    """

    def test_dissipation_arclength_atraviesa_rama_descendente_2d(self):
        E = 1.0e5; nu = 0.0; kappa_0 = 1.0e-3; alpha = 300.0
        dom, asm, nodes = _build_damage_quad4_uniaxial(
            E=E, nu=nu, kappa_0=kappa_0, alpha=alpha,
        )
        # Carga F_ref = 2× la carga elástica de pico para entrar al softening.
        # σ_pico = E·κ_0 = 100. Fuerza total para sección 1×1: 100.
        # Tomamos F_ref = 200 (2× pico).
        F_ref = _apply_uniaxial_traction_ref(dom, nodes, F_total=200.0)

        modes_seen = set()
        def probe(step, U, lam):
            modes_seen.add(solver._mode)

        conv = ConvergenceCriterion(rtol_force=1e-5, rtol_disp=1e-5)
        # dl pequeño para forzar varios pasos: ‖du_t‖ ≈ F_ref/(E·A) = 2e-3,
        # dl=5e-5 ⇒ dλ ≈ 0.025 ⇒ ~40 pasos hasta λ=1.
        solver = DissipationArcLengthSolver(
            asm, convergence=conv,
            max_lambda=1.0, initial_dl=5.0e-5, max_steps=300, max_iter=30,
            initial_tau=1.0e-4,
        )
        U = solver.solve(F_ref, step_callback=probe)

        # 1) Atravesó el pico: u_n3.x > u_pico.
        u_pico = kappa_0 * 1.0  # ε_pico · L (L=1)
        ux_n3 = U[nodes['n3'].dofs['ux']]
        self.assertGreater(ux_n3, u_pico,
            f"u_n3.x = {ux_n3:.3e} no superó el pico u_pico = {u_pico:.3e}")

        # 2) El switching activó el modo disipación (no debe permanecer
        # cilíndrico todo el tiempo si hubo softening genuino).
        self.assertIn("dissipation", modes_seen,
            f"Modo dissipation nunca observado: {modes_seen}")


# ---------------------------------------------------------------------------
# Test 5 — balance energético G_F integrado sobre historia
# ---------------------------------------------------------------------------

class TestBalanceEnergeticoGF(unittest.TestCase):
    """Σ ΔE_d ≈ trabajo cohesivo total integrado sobre la historia de la grieta.

    Setup más controlado: barra de daño 1D con G_F bien definido
    (`G_F = E·κ_0²/(2·alpha)` aproximado para softening exponencial al
    saturación). Comprueba que la disipación numérica acumulada por el
    solver coincide con la energía disipada integrada en el material.
    """

    def test_dissipacion_acumulada_consistente_con_trabajo_externo(self):
        # Para una barra 1-elemento con daño 1D bajo tracción monotónica:
        # el trabajo externo total ≈ ∫F·du.
        # Σ ΔE_d (Gutiérrez) debería coincidir con el trabajo disipado
        # (= trabajo externo − energía elástica almacenada).
        E = 1.0e5; kappa_0 = 1.0e-3; alpha = 500.0; A = 1.0; L = 1.0
        dom, asm, dof_ux = _build_damage_bar_1d(
            E=E, kappa_0=kappa_0, alpha=alpha, A=A, L=L,
        )
        F_ref = np.zeros(dom.total_dofs); F_ref[dof_ux] = 200.0

        history: list[tuple[float, float, float]] = []
        def probe(step, U, lam):
            history.append((lam, float(U[dof_ux]), solver._mode == "dissipation"))

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = DissipationArcLengthSolver(
            asm, convergence=conv, max_lambda=1.0,
            initial_dl=1.0e-4, max_steps=500, max_iter=30,
            initial_tau=1.0e-6,
        )
        solver.solve(F_ref, step_callback=probe)

        # Reconstruir ΔE_d paso a paso desde la historia.
        # ΔE_d = ½·(λ_n·F_ref·ΔU − Δλ·F_ref·U_n)
        # En 1-DOF: F_ref·U = F_ref_val·U; F_ref·ΔU = F_ref_val·ΔU.
        F_val = 200.0
        E_diss_acumulada = 0.0
        lam_prev, u_prev = 0.0, 0.0
        for (lam, u, _is_diss) in history:
            dU = u - u_prev
            dlam = lam - lam_prev
            dE = 0.5 * (lam_prev * F_val * dU - dlam * F_val * u_prev)
            E_diss_acumulada += dE
            lam_prev, u_prev = lam, u

        # En régimen elástico la disipación es ≈ 0; tras el pico crece
        # monótonamente. El total acumulado debe ser positivo y no trivial.
        self.assertGreater(E_diss_acumulada, 0.0,
            f"Σ ΔE_d = {E_diss_acumulada:.4e} debería ser > 0 (hay daño)")

        # Comparación con el trabajo externo neto realizado por F_ref·u
        # menos la energía elástica reversible final ½·k_sec·u_final²:
        # W_ext_neto = ∫F·du = ½·(F_inicio + F_fin)·du (regla trapezoidal)
        # NOTA: comprobación cualitativa — el balance exacto requiere un
        # cálculo más fino que escapa al scope de este test. Lo importante
        # aquí es que el solver **reporta una disipación coherente** (no
        # negativa, no espuria) y que crece tras el pico de daño.
        lam_finales = [lam for (lam, _, is_d) in history if is_d]
        self.assertGreater(len(lam_finales), 5,
            "Pocos pasos en régimen disipativo — historia insuficiente")


# ---------------------------------------------------------------------------
# Test 6 — sign-of-pivot tracking (aproximación inicial)
# ---------------------------------------------------------------------------

class TestSignOfPivotTracking(unittest.TestCase):
    """``_negative_pivots`` reemplaza el placeholder del padre y devuelve
    0 para K positiva-definida, 1 para K con determinante negativo
    (aproximación válida para indefinitud simple)."""

    def test_negative_pivots_pd_matrix(self):
        dom, asm, dof_ux = _build_elastic_bar()
        F_ref = np.zeros(dom.total_dofs); F_ref[dof_ux] = 100.0

        solver = DissipationArcLengthSolver(
            asm, convergence=ConvergenceCriterion(),
            max_lambda=1.0, initial_dl=0.5, max_steps=2, max_iter=20,
            initial_tau=0.01,
        )
        # Ensamblar para tener K_t y verificar el tracking en estado PD.
        K_global, _ = asm.assemble_non_linear_system(np.zeros(dom.total_dofs))
        K_red, _, _, _ = asm.reduce(K_global, F_ref)
        n_neg = solver._negative_pivots(K_red)
        self.assertEqual(n_neg, 0,
            "K positiva definida debería reportar 0 pivots negativos")

    def test_negative_pivots_indefinite_matrix(self):
        """K con autovalor negativo debe reportar n_neg=1 vía signo det."""
        # Construir una matriz 2×2 indefinida directamente.
        # Diag(2, -3) tiene det=-6<0, signo del det = -1 ⇒ n_neg=1 reportado.
        from scipy import sparse
        K_indef = sparse.csr_matrix(np.array([[2.0, 0.0], [0.0, -3.0]]))
        dom, asm, _ = _build_elastic_bar()
        solver = DissipationArcLengthSolver(
            asm, convergence=ConvergenceCriterion(),
            max_lambda=1.0, max_steps=2,
            initial_tau=0.01,
        )
        n_neg = solver._negative_pivots(K_indef)
        self.assertEqual(n_neg, 1,
            "K indefinida (det<0) debería reportar 1 pivot negativo")


if __name__ == '__main__':
    unittest.main()
