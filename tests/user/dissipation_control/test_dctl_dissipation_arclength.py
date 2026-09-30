"""Tests de ``DissipationArcLengthSolver`` — variante por disipación del arc-length.

Módulo de usuario ``dissipation_control`` (ADR 0020). Hasta el 2026-09-29 este
archivo era ``tests/test_dissipation_arclength.py``, en el programa principal.

Cubre los ``acceptance`` de la spec
``docs/user/dissipation_control/specs/DissipationArcLengthSolver.md`` que la
**primera implementación** alcanza:

1. **Recuperación cilíndrica** en régimen elástico puro — el solver coincide
   con :class:`ArcLengthSolver` a paridad de bits porque todos sus problemas
   de Newton son cilíndricos.
2. **Daño 1D con softening pronunciado** — el solver atraviesa la rama
   post-pico, resuelta con la restricción de disipación.
3. **Switching cilíndrico→disipación** en el paso siguiente al inicio del
   daño, sin oscilación.
4. **Daño 2D continuo (Quad4 + IsotropicDamage2D)** — análogo 2D del test
   3 con cinemática completa de Voigt.
5. **Balance energético** sobre la historia con daño 1D: la disipación de los
   pasos en modo de disipación coincide con la analítica del material.
6. **Daño 2D localizado con retroceso** (franja de Quad4 con
   ``IsotropicDamage2D``, medida energética, columna débil): arranque del modo
   de disipación tras el inicio del daño y balance de la disipación frente a
   la del material en los puntos de Gauss (revisión adversaria del
   2026-09-30).

Qué restricción resolvió cada paso lo dice el espía de ``_dctl_espia.py``,
que mira el modo del problema de Newton. Hasta el 2026-09-29 estos tests sólo
comprobaban el ``_mode`` que ve el callback y pasaban sin que la restricción
de disipación llegara a usarse: el predictor, con la tangente secante del
estado consolidado, devolvía cada paso a cilíndrico (spec, §5). Desde el
2026-09-30, con la tangente de carga en el predictor y el arranque por
``α_L``, la rama de ablandamiento se resuelve en modo de disipación.

**No cubierto en esta primera implementación** (deuda técnica documentada
en la spec):

- *Cohesivo + embedded discontinuity con K_e stiff* (componentes del módulo
  de usuario ``discontinuities``): la activación discreta del cohesivo
  (Rankine onset) introduce un salto en F_int/K_t que el `sign` por producto
  escalar no maneja graciosamente. Requiere sign-of-pivot tracking con una
  factorización LDLᵀ de Bunch-Kaufman (el solver no sigue hoy la inercia de
  K_t) o control por CMOD/CTOD del salto en lugar de τ global. Trabajo
  pendiente para una futura iteración del solver.

El diagnóstico ``_negative_pivots`` (signo del determinante de K_t) y sus dos
tests se retiraron el 2026-09-29: ningún solver lo llamaba.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
sys.path.insert(0, os.path.dirname(__file__))

import solidum  # noqa: E402
from solidum.constants import DAMAGE_MAX  # noqa: E402
from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.solid_2d import Quad4  # noqa: E402
from solidum.elements.truss import Truss2D  # noqa: E402
from solidum.materials.damage_1d import IsotropicDamage1D  # noqa: E402
from solidum.materials.damage_2d import IsotropicDamage2D  # noqa: E402
from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.convergence import ConvergenceCriterion  # noqa: E402
from solidum.math.solvers import ArcLengthSolver  # noqa: E402

from _dctl_espia import convergidos, espiar_problemas_de_newton  # noqa: E402

solidum.load_user_module("dissipation_control")
from solidum.user.dissipation_control import DissipationArcLengthSolver  # noqa: E402


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
        intentos = espiar_problemas_de_newton(diss)
        U_p = parent.solve(F_ref_p)
        U_d = diss.solve(F_ref_d)

        # En régimen puramente elástico (sin disipación), el modo debe
        # permanecer cilíndrico todo el tiempo: la tangente de carga coincide
        # con la secante y α_L es sólo el residuo de equilibrio.
        self.assertEqual(diss._mode, "cylindrical",
                          "Modo debió permanecer cilíndrico en régimen elástico")
        self.assertTrue(intentos)
        self.assertEqual({m for m, _ in intentos} - {"newton"}, {"cylindrical"},
                         f"problemas de Newton no cilíndricos en régimen elástico: {intentos}")

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
        u_y_modo: list[tuple[float, str]] = []   # (u, modo del problema de Newton)
        def probe(step, U, lam):
            modes_seen.add(solver._mode)
            u_y_modo.append((float(U[dof_ux]), intentos[-1][0]))

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
        intentos = espiar_problemas_de_newton(solver)
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
        # Y no sólo en el callback: los pasos posteriores al pico se
        # resolvieron con la restricción de disipación.
        post_pico = [modo for u, modo in u_y_modo if u > u_at_peak]
        self.assertGreater(convergidos(intentos, "dissipation"), 0)
        self.assertIn("dissipation", post_pico,
            f"Ningún paso post-pico en modo de disipación: {set(post_pico)}")


# ---------------------------------------------------------------------------
# Test 3 — switching automático observado tras la transición
# ---------------------------------------------------------------------------

class TestSwitchingCylindricalToDissipation(unittest.TestCase):
    """Acceptance ``switching_cilindrico_a_disipacion_en_pico``: si el daño
    empieza en el paso convergido N, los pasos 1..N son cilíndricos y desde
    N+1 en adelante son de disipación, sin oscilación. El paso N es
    cilíndrico por construcción: el arranque del modo de disipación se decide
    con la tangente de carga ``K_L`` de un paso convergido, que sólo se separa
    de la secante cuando en ese paso el daño creció."""

    def test_modo_cambia_solo_tras_yield(self):
        E = 1.0e5; kappa_0 = 1.0e-3; alpha = 500.0
        dom, asm, dof_ux = _build_damage_bar_1d(
            E=E, kappa_0=kappa_0, alpha=alpha,
        )
        F_ref = np.zeros(dom.total_dofs); F_ref[dof_ux] = 200.0
        barra = dom.elements[1]

        # (paso, _mode del solver, λ, modo del problema de Newton, daño)
        modes_per_step: list[tuple[int, str, float, str, float]] = []
        def probe(step, U, lam):
            modes_per_step.append((step, solver._mode, lam, intentos[-1][0],
                                   barra.state.vars[0]['damage']))

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = DissipationArcLengthSolver(
            asm, convergence=conv, max_lambda=1.0,
            initial_dl=1.0e-4, max_steps=500, max_iter=30,
            initial_tau=1.0e-5,
        )
        intentos = espiar_problemas_de_newton(solver)
        solver.solve(F_ref, step_callback=probe)

        # Encontrar el primer paso 'dissipation'.
        first_diss = next((i for (i, m, *_) in modes_per_step if m == "dissipation"), None)
        self.assertIsNotNone(first_diss,
            "Nunca se activó el modo dissipation")

        # Antes de ese paso, todos los modos son cilíndricos.
        for (i, m, *_) in modes_per_step:
            if i < first_diss:
                self.assertEqual(m, "cylindrical",
                    f"Paso {i} antes del switch debería ser cilíndrico, era {m}")

        # Con los problemas de Newton: cilíndricos hasta el paso N en que
        # empieza el daño (incluido) y de disipación desde el N+1, sin volver
        # a cilíndrico (el trazado acaba por max_steps con el daño lejos de la
        # saturación).
        newton = [r[3] for r in modes_per_step]
        n_dano = next(k for k, r in enumerate(modes_per_step) if r[4] > 0.0)
        self.assertEqual(set(newton[:n_dano + 1]), {"cylindrical"}, newton[:n_dano + 1])
        self.assertGreater(len(newton), n_dano + 1)
        self.assertEqual(set(newton[n_dano + 1:]), {"dissipation"},
            f"oscilación de modo tras el pico: {sorted(set(newton[n_dano + 1:]))}")


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
        u_y_modo: list[tuple[float, str]] = []   # (u_n3.x, modo del problema de Newton)
        def probe(step, U, lam):
            modes_seen.add(solver._mode)
            u_y_modo.append((float(U[nodes['n3'].dofs['ux']]), intentos[-1][0]))

        conv = ConvergenceCriterion(rtol_force=1e-5, rtol_disp=1e-5)
        # dl pequeño para forzar varios pasos: ‖du_t‖ ≈ F_ref/(E·A) = 2e-3,
        # dl=5e-5 ⇒ dλ ≈ 0.025 ⇒ ~40 pasos hasta λ=1.
        solver = DissipationArcLengthSolver(
            asm, convergence=conv,
            max_lambda=1.0, initial_dl=5.0e-5, max_steps=300, max_iter=30,
            initial_tau=1.0e-4,
        )
        intentos = espiar_problemas_de_newton(solver)
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
        # 3) Y los pasos posteriores al pico se resolvieron con la
        # restricción de disipación, no sólo el modo del callback.
        post_pico = [modo for u, modo in u_y_modo if u > u_pico]
        self.assertGreater(convergidos(intentos, "dissipation"), 0)
        self.assertIn("dissipation", post_pico,
            f"Ningún paso post-pico en modo de disipación: {set(post_pico)}")


# ---------------------------------------------------------------------------
# Test 5 — balance energético cualitativo sobre la historia de daño 1D
# ---------------------------------------------------------------------------

class TestBalanceEnergeticoCualitativo(unittest.TestCase):
    """Disipación acumulada Σ ΔE_d sobre la historia de una barra de daño
    continuo 1D (``IsotropicDamage1D``, un elemento, tracción monótona).

    No hay grieta ni ley cohesiva: el balance ``Σ ΔE_d ≈ G_F · l_d`` del caso
    cohesivo + embebido (acceptance ``balance_energetico_GF`` de la spec)
    sigue sin cubrir. Aquí se reconstruye ``ΔE_d`` paso a paso desde la
    historia ``(λ, u)`` del callback y se comprueba que la disipación
    acumulada es positiva, que el solver da varios pasos en modo de
    disipación tras el pico y (desde el 2026-09-30, cuando esos pasos
    empezaron a usar de verdad la restricción) que la disipación de los pasos
    en modo de disipación coincide con la analítica del material.
    """

    E = 1.0e5; KAPPA_0 = 1.0e-3; ALPHA = 500.0; A = 1.0; L = 1.0
    F_VAL = 200.0

    @classmethod
    def setUpClass(cls):
        dom, asm, dof_ux = _build_damage_bar_1d(
            E=cls.E, kappa_0=cls.KAPPA_0, alpha=cls.ALPHA, A=cls.A, L=cls.L,
        )
        F_ref = np.zeros(dom.total_dofs); F_ref[dof_ux] = cls.F_VAL

        # (λ, u, modo del problema de Newton del paso)
        history: list[tuple[float, float, str]] = []
        def probe(step, U, lam):
            history.append((lam, float(U[dof_ux]), intentos[-1][0]))

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = DissipationArcLengthSolver(
            asm, convergence=conv, max_lambda=1.0,
            initial_dl=1.0e-4, max_steps=500, max_iter=30,
            initial_tau=1.0e-6,
        )
        intentos = espiar_problemas_de_newton(solver)
        solver.solve(F_ref, step_callback=probe)
        cls.history = history

        # ΔE_d paso a paso desde la historia:
        # ΔE_d = ½·(λ_n·F_ref·ΔU − Δλ·F_ref·U_n); en 1 DOF, F_ref·U = F_val·u.
        cls.pasos = []   # (u_n, ΔE_d, modo)
        lam_prev, u_prev = 0.0, 0.0
        for (lam, u, modo) in history:
            dE = 0.5 * (lam_prev * cls.F_VAL * (u - u_prev) - (lam - lam_prev) * cls.F_VAL * u_prev)
            cls.pasos.append((u_prev, dE, modo))
            lam_prev, u_prev = lam, u

    @classmethod
    def _disipacion_analitica(cls, kappa):
        """Disipación acumulada del material hasta ``κ ≥ κ₀`` en carga
        monótona: ``E_d = W − Ψ``, con ``W = A·L·∫₀^κ σ dε`` y ``Ψ = ½·A·L·σ·κ``,
        ``σ = (1 − d)·E·ε = E·κ₀·exp(−α(ε − κ₀))`` en la rama de
        ablandamiento (antes del tope ``DAMAGE_MAX``). Con daño isótropo
        ``E_d`` sólo depende de ``κ`` (``Ė_d = ½·E·κ²·ḋ``), no del camino."""
        E, k0, a = cls.E, cls.KAPPA_0, cls.ALPHA
        expo = np.exp(-a * (kappa - k0))
        W = 0.5 * E * k0 ** 2 + (E * k0 / a) * (1.0 - expo)
        return cls.A * cls.L * (W - 0.5 * E * k0 * expo * kappa)

    def test_dissipacion_acumulada_consistente_con_trabajo_externo(self):
        E_diss_acumulada = sum(dE for _, dE, _ in self.pasos)

        # En régimen elástico la disipación es ≈ 0; tras el pico crece
        # monótonamente. El total acumulado debe ser positivo y no trivial.
        self.assertGreater(E_diss_acumulada, 0.0,
            f"Σ ΔE_d = {E_diss_acumulada:.4e} debería ser > 0 (hay daño)")

        # Pasos resueltos con la restricción de disipación (modo del problema
        # de Newton; hasta el 2026-09-29 se contaba el _mode del callback,
        # que no implicaba que la restricción se usara).
        pasos_disipacion = [dE for _, dE, modo in self.pasos if modo == "dissipation"]
        self.assertGreater(len(pasos_disipacion), 5,
            "Pocos pasos en régimen disipativo — historia insuficiente")
        # Cada uno disipa τ > 0: la restricción fija ΔE_d = τ.
        self.assertGreater(min(pasos_disipacion), 0.0)

    def test_disipacion_de_los_pasos_de_disipacion_igual_a_la_del_material(self):
        """Σ ΔE_d de los pasos en modo de disipación = E_d(κ_fin) − E_d(κ_ini).

        Los pasos en modo de disipación forman un tramo contiguo de carga
        monótona del único elemento (``κ = ε = u/L``) desde el estado
        consolidado previo al primero (``κ_ini``) hasta el último
        (``κ_fin``). La suma telescópica de ``ΔE_d`` es la regla del trapecio
        del trabajo externo menos el cambio de ``½·λ·F·u``. Tolerancia
        derivada: con ``rtol_force = 1e-6`` cada estado convergido cumple
        ``|λ·F − σ·A| ≤ 1e-6·λ·F``, y el error dominante está en los dos
        términos de contorno, ``≤ 2·1e-6·½·σ·A·u ≈ 1.1e-7 J`` frente a
        ``ΔE_d ≈ 2.5e-3 J`` (≈ 4.5e-5 relativo); el error del trapecio con
        pasos de ``Δε ~ 1e-7`` es de ``~1e-13 J``. Se exige 1e-4.
        """
        modos = [modo for _, _, modo in self.pasos]
        self.assertIn("dissipation", modos, "ningún paso en modo de disipación")
        primero = modos.index("dissipation")
        ultimo = len(modos) - 1 - modos[::-1].index("dissipation")
        self.assertEqual(set(modos[primero:ultimo + 1]), {"dissipation"},
                         "los pasos en modo de disipación no forman un tramo contiguo")
        kappa_ini = self.pasos[primero][0] / self.L            # u_n del primer paso
        kappa_fin = self.history[ultimo][1] / self.L
        self.assertGreater(kappa_ini, self.KAPPA_0)            # empieza tras el pico
        suma = sum(dE for _, dE, _ in self.pasos[primero:ultimo + 1])
        exacta = self._disipacion_analitica(kappa_fin) - self._disipacion_analitica(kappa_ini)
        self.assertAlmostEqual(suma / exacta, 1.0, delta=1e-4,
            msg=f"Σ ΔE_d = {suma:.6e} J frente a {exacta:.6e} J del material")


# ---------------------------------------------------------------------------
# Test 6 — daño 2D localizado con retroceso: arranque y balance energético
# ---------------------------------------------------------------------------

class TestDamage2DLocalizadoBalance(unittest.TestCase):
    """Franja de Quad4 (4×2, tensión plana, ν = 0.2) con ``IsotropicDamage2D``
    y la medida de deformación equivalente por omisión (energética), con la
    columna central un 5 % más débil. Revisión adversaria del 2026-09-30.

    Es el análogo 2D de la barra con retroceso (``test_dctl_retroceso.py``):
    el daño se localiza en la columna débil, el resto descarga
    elásticamente y el desplazamiento del extremo retrocede tras el pico
    (``(n − 1)·α·κ₀w = 1.43 > 1``). Con ``ν ≠ 0`` la columna dañada se contrae
    lateralmente distinto que sus vecinas y el campo deja de ser uniaxial:
    cinemática 2D completa.

    Se comprueba que el solver arranca en disipación en el paso siguiente al
    inicio del daño y sigue así mientras la columna se ablanda, y que la
    disipación de esos pasos coincide con la del material. Con la medida
    energética ``Y = ½·ε:C_e:ε = ½·E·ε_eq²``, así que la disipación por unidad
    de volumen sólo depende de ``κ``: ``D(κ) = ∫ ½·E·κ²·d'(κ) dκ``, la misma
    ``W − Ψ`` de tracción uniaxial que el test 5, integrada en los puntos de
    Gauss (``V_gp = t·h_x·h_y/4`` en un Quad4 rectangular con 2×2 puntos).
    """

    E = 30.0e9; NU = 0.2; K0 = 1.0e-4; WEAK = 0.95; ALPHA = 5000.0
    NX = 4; NY = 2; LX = 1.0; LY = 0.2; T = 0.1

    @classmethod
    def _modelo(cls):
        hx, hy = cls.LX / cls.NX, cls.LY / cls.NY
        dom = Domain()
        nid = {}
        for j in range(cls.NY + 1):
            for i in range(cls.NX + 1):
                nid[i, j] = dom.add_node(len(nid) + 1, [i * hx, j * hy])
        for j in range(cls.NY):
            for i in range(cls.NX):
                k0 = cls.WEAK * cls.K0 if i == cls.NX // 2 else cls.K0
                mat = IsotropicDamage2D(E=cls.E, nu=cls.NU, kappa_0=k0, alpha=cls.ALPHA)
                dom.add_element(Quad4(len(dom.elements) + 1,
                                      [nid[i, j], nid[i + 1, j], nid[i + 1, j + 1], nid[i, j + 1]],
                                      mat, thickness=cls.T))
        for j in range(cls.NY + 1):
            nid[0, j].fix_dof('ux', 0.0)
        nid[0, 0].fix_dof('uy', 0.0)
        dom.generate_equation_numbers(verbose=False)
        # Carga de referencia: el pico nominal de la columna débil, repartida
        # en el borde derecho con los pesos consistentes (½, 1, ½)/NY.
        F_peak = cls.E * cls.WEAK * cls.K0 * cls.LY * cls.T
        F = np.zeros(dom.total_dofs)
        for j in range(cls.NY + 1):
            F[nid[cls.NX, j].dofs['ux']] = (0.5 if j in (0, cls.NY) else 1.0) * F_peak / cls.NY
        return dom, Assembler(dom), F, [nid[cls.NX, j] for j in range(cls.NY + 1)], hx * hy

    @classmethod
    def setUpClass(cls):
        dom, asm, F, borde, area = cls._modelo()
        cls.F = F
        cls.V_gp = cls.T * area / 4.0
        els = list(dom.elements.values())
        cls.debil = np.array([el.material.kappa_0 < cls.K0 for el in els])
        cls.k0 = np.array([el.material.kappa_0 for el in els])

        # (λ, F·U, ‖U‖, u medio del borde, κ y d por punto de Gauss, modo)
        cls.pasos = []
        def probe(step, U, lam):
            cls.pasos.append(dict(
                lam=lam, q=float(F @ U), nU=float(np.linalg.norm(U)),
                u=float(np.mean([U[n.dofs['ux']] for n in borde])),
                kappa=np.array([[el.state.vars[g]['kappa'] for g in range(4)] for el in els]),
                d=np.array([[el.state.vars[g]['damage'] for g in range(4)] for el in els]),
                modo=intentos[-1][0]))

        # τ₀ ≈ 1/300 de la disipación de la columna débil hasta la saturación
        # (≈ 3.5 J): unos 70 pasos de ablandamiento con tau_max_factor = 5.
        cls.solver = DissipationArcLengthSolver(
            asm, max_iter=30, max_lambda=1.2, max_steps=100,
            initial_dlambda=0.02, initial_tau=1.0e-2,
        )
        intentos = espiar_problemas_de_newton(cls.solver)
        cls.solver.solve(F, step_callback=probe)

        d_max = np.array([p['d'].max() for p in cls.pasos])
        cls.inicio = int(np.flatnonzero(d_max > 0.0)[0])
        # Tramo de ablandamiento continuo: hasta el último paso sin ningún
        # punto saturado (el paso que cruza DAMAGE_MAX disipa menos que τ;
        # spec §10).
        saturados = np.flatnonzero(d_max >= DAMAGE_MAX)
        cls.fin = int(saturados[0] - 1) if len(saturados) else len(cls.pasos) - 1

    @classmethod
    def _D(cls, kappa, k0):
        """Disipación por unidad de volumen hasta ``κ`` (0 por debajo de κ₀)."""
        expo = np.exp(-cls.ALPHA * (kappa - k0))
        W = 0.5 * cls.E * k0 ** 2 + (cls.E * k0 / cls.ALPHA) * (1.0 - expo)
        return np.where(kappa > k0, W - 0.5 * cls.E * k0 * expo * kappa, 0.0)

    def _disipacion_material(self, paso):
        return self.V_gp * float(np.sum(self._D(paso['kappa'], self.k0[:, None])))

    def test_el_dano_se_localiza_y_el_extremo_retrocede(self):
        for p in self.pasos:
            self.assertFalse(np.any(p['d'][~self.debil] > 0.0),
                             "daño fuera de la columna débil")
        u = np.array([p['u'] for p in self.pasos])
        self.assertLess(u[self.inicio:self.fin + 1].min(), u[self.inicio],
                        "sin retroceso: el caso no ejercita la rama con α_L < 0")

    def test_arranca_en_disipacion_tras_el_inicio_del_dano(self):
        modos = [p['modo'] for p in self.pasos]
        self.assertEqual(set(modos[:self.inicio + 1]), {"cylindrical"}, modos[:self.inicio + 1])
        tramo = modos[self.inicio + 1:self.fin + 1]
        self.assertGreater(len(tramo), 30)
        self.assertEqual(set(tramo), {"dissipation"},
                         f"modos en la rama de ablandamiento: {sorted(set(tramo))}")

    def test_disipacion_de_los_pasos_igual_a_la_del_material(self):
        """``Σ ΔE_d`` de los pasos de disipación = ``Σ_gp V_gp·[D(κ_fin) − D(κ_ini)]``.

        ``Σ ΔE_d`` es, por telescopía, la regla del trapecio del trabajo
        externo ``∫λ d(F·U)`` menos el cambio de ``½·λ·F·U``. Cota del error,
        sin márgenes ajustados:

        - ``Ψ = ½·U·F_int`` difiere de ``½·λ·F·U`` en ``½·U·r``, con
          ``‖r‖ ≤ atol + rtol·max(λ‖F‖, ‖F_int‖)`` (criterio de convergencia) y
          ``‖F_int‖ ≤ λ·√(‖F‖² + F_total²)`` (cargas más reacciones, que suman
          ``−λ·F_total`` con el mismo signo): ``½·‖U‖·‖r‖`` en cada extremo;
        - el trapecio de ``W`` sobre el camino discreto es de segundo orden;
          la duplicación de paso estima su error en ``(T_h − T_2h)/3`` y se
          toma ``|T_h − T_2h|``.

        Medido (67 pasos de disipación, ``ΔD = 3.2307 J``): diferencia
        5.9e-4 J (1.8e-4 relativo) frente a una estimación del trapecio de
        5.3e-4 J y una cota del residuo de 8.1e-5 J; el trapecio domina, con
        pasos de ``τ`` de hasta 5e-2 J sobre una rama curva. Con el solver
        anterior al 2026-09-30 (pasos cilíndricos en la rama) hay pasos del
        tramo con ``ΔE_d < 0`` (−1.5e-11 J) y el test falla.
        """
        a, b = self.inicio, self.fin
        lam = np.array([p['lam'] for p in self.pasos])
        q = np.array([p['q'] for p in self.pasos])
        dE = 0.5 * (lam[a:b] * np.diff(q[a:b + 1]) - np.diff(lam[a:b + 1]) * q[a:b])
        self.assertGreater(dE.min(), 0.0)
        suma = float(dE.sum())
        exacta = self._disipacion_material(self.pasos[b]) - self._disipacion_material(self.pasos[a])

        conv = self.solver.convergence
        nF, F_tot = float(np.linalg.norm(self.F)), float(self.F.sum())
        cota_r = sum(0.5 * self.pasos[k]['nU']
                     * (conv.atol_force + conv.rtol_force * abs(lam[k]) * np.hypot(nF, F_tot))
                     for k in (a, b))

        def trapecio(idx):
            return float(np.sum(0.5 * (lam[idx[:-1]] + lam[idx[1:]]) * np.diff(q[idx])))
        idx = np.arange(a, b + 1)
        idx2 = idx[::2] if idx[::2][-1] == b else np.append(idx[::2], b)
        cota_t = abs(trapecio(idx) - trapecio(idx2))

        self.assertLessEqual(abs(suma - exacta), cota_r + cota_t,
            msg=f"Σ ΔE_d = {suma:.6e} J frente a {exacta:.6e} J del material "
                f"(cotas: residuo {cota_r:.2e}, trapecio {cota_t:.2e})")


if __name__ == '__main__':
    unittest.main()
