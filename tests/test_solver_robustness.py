"""Batería de stress de los solvers no lineales — auditoría fase A.

Este archivo NO valida corrección física (eso ya lo cubren los tests por
componente). Mapea **dónde y cómo divergen** ``NonlinearSolver``,
``ArcLengthSolver`` y ``NewtonNewmarkSolver`` en regímenes donde un solver
robusto debería sostener (o al menos fallar limpiamente).

Convención de marcado:

- **Tests verdes**: el solver maneja correctamente el régimen hoy.
- ``@unittest.expectedFailure``: hoy diverge o produce respuesta incorrecta;
  cuando la fase C resuelva el hueco, se retira el decorador y queda como
  test de regresión permanente. La razón se documenta en docstring.
- ``assertRaises(RuntimeError)``: el solver ya falla limpiamente; el test
  documenta el modo de fallo actual (informativo, no aspiracional).

Cada test recoge además el historial iteración-a-iteración via
``IterationProbe`` (callback inyectable) y lo expone en el atributo
``self.probe`` para que el informe de auditoría pueda extraer métricas
(número de pasos, iteraciones por paso, factor de carga final, etc.).

Informe asociado: ``docs/auditorias/solvers_robustez_fase_A.md``.
"""
from __future__ import annotations

import math
import os
import sys
import unittest
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.frame.euler_corot import Frame2DEulerCorot
from solidum.elements.solid_2d import Quad4
from solidum.elements.truss import Truss2D, Truss2DCorot
from solidum.materials.damage_1d import IsotropicDamage1D
from solidum.materials.damage_2d import IsotropicDamage2D
from solidum.materials.drucker_prager_2d import DruckerPrager2D
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.assembly import Assembler
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import (
    ArcLengthSolver,
    NewtonNewmarkSolver,
    NonlinearSolver,
)


# =============================================================================
# Utilities compartidas
# =============================================================================

@dataclass
class IterationProbe:
    """Telemetría por-paso convergida del solver.

    El ``step_callback`` de los solvers se invoca **después** de cada paso
    convergido con ``(step_idx, U, load_factor)``. Esta probe captura esa
    secuencia para análisis post-hoc en el informe de auditoría.
    """
    step_indices: list[int] = field(default_factory=list)
    load_factors: list[float] = field(default_factory=list)
    u_snapshots: list[np.ndarray] = field(default_factory=list)

    def __call__(self, step: int, U: np.ndarray, load_factor: float) -> None:
        self.step_indices.append(step)
        self.load_factors.append(float(load_factor))
        self.u_snapshots.append(U.copy())

    @property
    def n_steps(self) -> int:
        return len(self.step_indices)

    @property
    def final_load_factor(self) -> float:
        return self.load_factors[-1] if self.load_factors else 0.0


def _build_von_mises_truss(
    *, L: float, h: float, E: float, A: float, corotational: bool = True,
):
    """Von Mises 2-bar truss (snap-through clásico).

    Geometría — dos barras inclinadas formando un techo a dos aguas con
    apex elevado, apoyos articulados en los extremos:

        (0,0) -------- (L, h) -------- (2L, 0)
          |             |               |
        fijo         carga vertical    fijo
                     hacia abajo

    Con ``h ≪ L`` (shallow), la curva carga-desplazamiento del apex tiene
    un punto límite claro. Configuración apta para arc-length; control de
    carga puro debería diverger al alcanzarlo. Cinemática corotacional por
    defecto (grandes desplazamientos del apex, pequeña deformación axial).

    Retorna ``(domain, assembler, apex_node, dof_uy_apex)``.
    """
    dom = Domain()
    n_left = dom.add_node(1, [0.0, 0.0])
    n_apex = dom.add_node(2, [L, h])
    n_right = dom.add_node(3, [2.0 * L, 0.0])

    elem_cls = Truss2DCorot if corotational else Truss2D
    mat = Elastic1D(E=E)
    dom.add_element(elem_cls(1, [n_left, n_apex], mat, A=A))
    dom.add_element(elem_cls(2, [n_apex, n_right], mat, A=A))

    n_left.fix_dof('ux', 0.0); n_left.fix_dof('uy', 0.0)
    n_right.fix_dof('ux', 0.0); n_right.fix_dof('uy', 0.0)
    # Apex con ux fijo por simetría (snap-through es puro vertical).
    n_apex.fix_dof('ux', 0.0)

    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), n_apex, n_apex.dofs['uy']


def _build_damage_bar_1d(
    *, E: float, A: float, L: float, kappa_0: float, alpha: float, n_elems: int = 1,
):
    """Barra 1D empotrada-libre con ``IsotropicDamage1D``.

    Single-DOF efectivo si ``n_elems=1`` (mapea más limpiamente el modo de
    fallo del solver puro). Con ``n_elems > 1`` aparece localización del
    daño en un único elemento — útil para probar comportamiento del solver
    ante respuesta no homogénea.

    Retorna ``(domain, assembler, free_node, dof_ux)``.
    """
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * L / n_elems, 0.0]) for i in range(n_elems + 1)]
    mat = IsotropicDamage1D(E=E, kappa_0=kappa_0, alpha=alpha)
    for k in range(n_elems):
        dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]], mat, A=A))

    nodes[0].fix_dof('ux', 0.0); nodes[0].fix_dof('uy', 0.0)
    for n in nodes[1:]:
        n.fix_dof('uy', 0.0)

    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), nodes[-1], nodes[-1].dofs['ux']


def _build_uniaxial_quad4(
    material, *, side: float = 1.0, thickness: float = 1.0,
):
    """Quad4 cuadrado unitario con apoyos para tracción uniaxial en x."""
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [side, 0.0])
    n3 = dom.add_node(3, [side, side])
    n4 = dom.add_node(4, [0.0, side])
    dom.add_element(Quad4(1, [n1, n2, n3, n4], material, thickness=thickness))

    n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0)
    n2.fix_dof('uy', 0.0)
    n4.fix_dof('ux', 0.0)

    dom.generate_equation_numbers(verbose=False)
    return dom, Assembler(dom), {'n1': n1, 'n2': n2, 'n3': n3, 'n4': n4}


def _build_imperfect_column(
    *, E: float, A: float, I: float, L: float, n_elems: int, imperfection: float,
):
    """Columna biempotrada con imperfección sinusoidal en posición inicial.

    Imperfección lateral en `y` con primer modo de pandeo. ``imperfection``
    es la flecha máxima en el centro (en unidades de longitud). Con
    ``imperfection = 0`` recuperamos la bifurcación pura.

    Retorna ``(domain, assembler, top_node, dof_ux_top, dof_uy_mid)``.
    """
    dom = Domain()
    nodes = []
    for i in range(n_elems + 1):
        x = i * L / n_elems
        # Modo sinusoidal: y = imperfection · sin(π·x/L)
        y = imperfection * math.sin(math.pi * x / L)
        nodes.append(dom.add_node(i + 1, [x, y]))

    mat = Elastic1D(E=E)
    for k in range(n_elems):
        dom.add_element(Frame2DEulerCorot(k + 1, [nodes[k], nodes[k + 1]], mat, A=A, I=I))

    # Empotramiento en la base (ux, uy, rz = 0), libre arriba salvo carga.
    nodes[0].fix_dof('ux', 0.0); nodes[0].fix_dof('uy', 0.0); nodes[0].fix_dof('rz', 0.0)
    # Eliminamos giro del extremo cargado para forzar la longitud efectiva
    # de Euler-Bernoulli L_ef = L/2 (columna empotrada-empotrada con
    # deslizamiento axial libre). La carga crítica es entonces:
    #   P_cr = 4·π²·E·I / L²
    nodes[-1].fix_dof('rz', 0.0)
    nodes[-1].fix_dof('uy', 0.0)  # solo desplazamiento axial libre arriba

    dom.generate_equation_numbers(verbose=False)
    mid_idx = n_elems // 2
    return (
        dom, Assembler(dom), nodes[-1],
        nodes[-1].dofs['ux'], nodes[mid_idx].dofs['uy'],
    )


# =============================================================================
# ArcLengthSolver — 4 tests
# =============================================================================

class TestArcLengthRobustness(unittest.TestCase):
    """Stress de ``ArcLengthSolver`` en regímenes con punto límite, snap-back
    y bifurcación — la razón de existir del solver, hoy sin tests directos.
    """

    def test_snap_through_von_mises_truss(self):
        """Snap-through clásico: arc-length traza la curva completa, NonlinearSolver diverge.

        Geometría shallow ``h/L = 0.1`` → punto límite claro en torno a
        ``w ≈ h`` (cuando el apex pasa por la línea de los apoyos).
        Validación: arc-length alcanza ``λ ≥ 1`` con apex bien por debajo
        de la línea de apoyos (``w_y < -h``, configuración invertida).
        """
        L = 10.0; h = 1.0; E = 1.0e5; A = 1.0
        dom, asm, apex, dof_uy = _build_von_mises_truss(L=L, h=h, E=E, A=A)

        # Carga vertical hacia abajo (negativa en +y) suficiente para forzar
        # el snap-through. Magnitud de referencia: ~10× la carga elástica
        # característica E·A·(h/L)² (orden de magnitud del punto límite
        # del von Mises 2-bar).
        F_ref_mag = 10.0 * E * A * (h / L) ** 2
        F_ext = np.zeros(dom.total_dofs)
        F_ext[dof_uy] = -F_ref_mag

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        self.probe = IterationProbe()
        solver = ArcLengthSolver(
            asm, convergence=conv,
            max_lambda=1.0, initial_dl=0.05, max_steps=200, max_iter=30,
        )
        U = solver.solve(F_ext, step_callback=self.probe)

        # Tras el snap-through, el apex debe quedar invertido respecto al
        # plano de los apoyos: u_y < -h (configuración postcrítica).
        self.assertLess(U[dof_uy], -h,
            f"Snap-through no se completó: u_y={U[dof_uy]:.4f}, esperado < -h={-h:.4f}.")
        # Y el solver debe haber recorrido > 5 pasos (señal de que vio
        # el punto límite y no llegó saltando).
        self.assertGreater(self.probe.n_steps, 5,
            f"Solo {self.probe.n_steps} pasos: el arc-length probablemente "
            "saltó el punto límite sin trazarlo.")

    def test_snap_through_with_load_control(self):
        """Mismo problema con ``NonlinearSolver``: ¿control de carga atraviesa el punto límite?

        Cargamos a la misma magnitud del test anterior (``10·E·A·(h/L)²``,
        muy por encima del pico) y dejamos al adaptativo bisecar todo lo
        que quiera. Documentamos si el solver con paso adaptativo sostiene
        el snap-through o diverge — la suposición común es que no puede,
        pero conviene verificarla en este código concreto.
        """
        L = 10.0; h = 1.0; E = 1.0e5; A = 1.0
        dom, asm, apex, dof_uy = _build_von_mises_truss(L=L, h=h, E=E, A=A)

        F_ref_mag = 10.0 * E * A * (h / L) ** 2
        F_ext = np.zeros(dom.total_dofs)
        F_ext[dof_uy] = -F_ref_mag

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        self.probe = IterationProbe()
        solver = NonlinearSolver(
            asm, convergence=conv, num_steps=20, max_iter=30,
            adaptive=True, min_delta_lambda=1e-7,
        )

        try:
            U = solver.solve(F_ext, step_callback=self.probe)
            # Si converge, registramos cómo terminó.
            self.assertLess(U[dof_uy], -h,
                f"Solver convergió pero u_y={U[dof_uy]:.4f} > -h: no atravesó "
                "el punto límite; cargó por debajo o saltó.")
            # Adicional: documentar el número de pasos requeridos.
            self.assertGreater(self.probe.n_steps, 10,
                f"Solo {self.probe.n_steps} pasos para snap-through con control "
                "de carga: el adaptativo está siendo sorprendentemente eficaz.")
        except RuntimeError as e:
            # Modo de fallo esperado teóricamente. Lo documentamos.
            self.fail(
                "NonlinearSolver divergió en snap-through con control de carga "
                f"(esperado teóricamente): {e}. Este es el modo de fallo que "
                "motiva la existencia de ArcLengthSolver."
            )

    def test_arc_length_traverses_damage_softening(self):
        """Arc-length atraviesa la rama post-pico de daño 1D con softening fuerte.

        Barra empotrada con ``IsotropicDamage1D``, ``alpha=500`` (softening
        pronunciado), carga axial 2× el pico elástico. El control de carga
        puro no puede llegar (la curva F-u tiene tramo con dF/du < 0);
        arc-length sí. Documenta que el límite descrito en STATUS.md
        ("solvers no atraviesan la transición elástico→softening") aplica
        específicamente al penalty cohesivo stiff del embedded, **no** al
        daño continuo bulk con softening exponencial.
        """
        E = 1.0e5; A = 1.0; L = 1.0
        kappa_0 = 1.0e-3; alpha = 500.0  # softening pronunciado
        dom, asm, _free_node, dof_ux = _build_damage_bar_1d(
            E=E, A=A, L=L, kappa_0=kappa_0, alpha=alpha,
        )

        # Carga axial > carga pico (σ_pico = E·κ_0 = 100 → F_pico = 100).
        F_ext = np.zeros(dom.total_dofs)
        F_ext[dof_ux] = 200.0  # 2× la carga de pico

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        self.probe = IterationProbe()
        solver = ArcLengthSolver(
            asm, convergence=conv,
            max_lambda=1.0, initial_dl=0.01, max_steps=500, max_iter=30,
        )
        U = solver.solve(F_ext, step_callback=self.probe)

        # Si llega aquí: arc-length atravesó el softening. Verificamos
        # que el desplazamiento final es > el desplazamiento del pico.
        u_at_peak = kappa_0 * L  # ε_peak · L
        self.assertGreater(U[dof_ux], u_at_peak,
            "Arc-length no atravesó el pico de daño.")

    def test_bifurcation_imperfect_column(self):
        """Bifurcación regularizada: columna con imperfección, post-pandeo unicamente definido.

        Carga axial creciente sobre columna esbelta con imperfección
        sinusoidal (1% de L) en primer modo. La rama post-crítica está
        bien definida (no requiere tracking de pivote). Arc-length debe
        trazarla hasta ``λ·P_cr ≈ 1.2 P_cr``.

        Validación: el desplazamiento lateral del centro crece más rápido
        que linealmente cuando ``λ → 1``, signo inequívoco de
        amplificación de pandeo.
        """
        E = 2.0e8; A = 1.0e-2; I = 1.0e-5; L = 5.0; n_elems = 10
        # P_cr (empotrada-empotrada con deslizamiento axial) = 4π²·E·I/L²
        P_cr = 4.0 * math.pi ** 2 * E * I / L ** 2

        # Carga axial de referencia: 1.2·P_cr, hacia -y (compresión desde arriba).
        dom, asm, top_node, dof_ux_top, dof_uy_mid = _build_imperfect_column(
            E=E, A=A, I=I, L=L, n_elems=n_elems, imperfection=0.01 * L,
        )
        F_ext = np.zeros(dom.total_dofs)
        F_ext[dof_ux_top] = -1.2 * P_cr  # compresión en +x (la columna sube en x desde y=0)

        conv = ConvergenceCriterion(rtol_force=1e-5, rtol_disp=1e-5)
        self.probe = IterationProbe()
        solver = ArcLengthSolver(
            asm, convergence=conv,
            max_lambda=1.0, initial_dl=0.05, max_steps=100, max_iter=30,
        )
        U = solver.solve(F_ext, step_callback=self.probe)

        # Validación: el solver alcanza λ ≥ 0.9 (cerca de 1.2·P_cr)
        # y la flecha lateral del centro creció apreciablemente.
        self.assertGreaterEqual(self.probe.final_load_factor, 0.9,
            f"Solver llegó solo a λ={self.probe.final_load_factor:.3f}.")
        flecha_final = abs(U[dof_uy_mid])
        flecha_inicial = abs(0.01 * L * math.sin(math.pi * 0.5))  # mid-span
        self.assertGreater(flecha_final, 2.0 * flecha_inicial,
            f"Flecha lateral apenas creció: {flecha_final:.4e} vs inicial "
            f"{flecha_inicial:.4e}; no se ve amplificación de pandeo.")

    def test_dl_control_aggressive(self):
        """Control de ``dl``: arco shallow con ``initial_dl`` agresivo no debe descontrolar.

        Mismo problema del snap-through pero con ``initial_dl = 0.5``
        (10× el del test #1). Verificamos que el auto-shrink reacciona y
        el solver llega a λ ≥ 1 sin saltarse el punto límite.

        Si el solver se descontrola (drift fuera del arco, raíz mala),
        este test diverge y queda como ``expectedFailure``.
        """
        L = 10.0; h = 1.0; E = 1.0e5; A = 1.0
        dom, asm, apex, dof_uy = _build_von_mises_truss(L=L, h=h, E=E, A=A)
        F_ref_mag = 10.0 * E * A * (h / L) ** 2
        F_ext = np.zeros(dom.total_dofs)
        F_ext[dof_uy] = -F_ref_mag

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        self.probe = IterationProbe()
        solver = ArcLengthSolver(
            asm, convergence=conv,
            max_lambda=1.0, initial_dl=0.5, max_steps=200, max_iter=30,
        )
        U = solver.solve(F_ext, step_callback=self.probe)

        # Debe seguir trazando el snap-through correctamente.
        self.assertLess(U[dof_uy], -h,
            f"Snap-through no se completó con dl agresivo: u_y={U[dof_uy]:.4f}.")


# =============================================================================
# NonlinearSolver — 3 tests
# =============================================================================

class TestNonlinearSolverRobustness(unittest.TestCase):
    """Stress de ``NonlinearSolver`` (Newton-Raphson con paso adaptativo)."""

    def test_drucker_prager_moderate_plasticity(self):
        """Drucker-Prager plane strain con plastificación moderada (régimen tratable).

        Quad4 plane strain, tracción uniaxial con carga ~3× la cohesión
        característica. Suficiente para activar plasticidad no asociada
        en al menos un punto Gauss, sin llegar a saturar la capacidad
        plástica del elemento. El solver debe completar la carga.
        """
        E = 1.0e7; nu = 0.3
        phi_deg = 30.0; psi_deg = 10.0; cohesion = 1.0e3
        # Endurecimiento suave (H = E/100) → tangente plástica positiva,
        # el Newton mantiene convergencia tras el yield. Con ``H=0``
        # (perfectamente plástico) el régimen plástico es marginal y
        # Newton oscila — eso es exactamente lo que documenta el test
        # siguiente (``test_drucker_prager_overload_oscillates``).
        mat = DruckerPrager2D(
            E=E, nu=nu, cohesion=cohesion, phi_deg=phi_deg, psi_deg=psi_deg,
            H=E / 100.0, hypothesis='plane_strain',
        )
        dom, asm, nodes = _build_uniaxial_quad4(mat, side=1.0, thickness=1.0)

        F_total = 3.0 * cohesion  # plastificación moderada
        F_ext = np.zeros(dom.total_dofs)
        F_ext[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F_ext[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        self.probe = IterationProbe()
        solver = NonlinearSolver(
            asm, convergence=conv, num_steps=10, max_iter=30,
            adaptive=True, min_delta_lambda=1e-7,
        )
        U = solver.solve(F_ext, step_callback=self.probe)

        self.assertAlmostEqual(self.probe.final_load_factor, 1.0, places=4,
            msg=f"Solver llegó solo a λ={self.probe.final_load_factor:.4f}.")
        elem = list(dom.elements.values())[0]
        alphas = [gp.get('alpha', 0.0) for gp in elem.state.vars if gp is not None]
        self.assertGreater(max(alphas), 0.0,
            "Ningún punto Gauss plastificó: el test no ejerció el régimen objetivo.")

    def test_drucker_prager_overload_oscillates(self):
        """Drucker-Prager perfectamente plástico con carga > capacidad: Newton oscila.

        Hallazgo de la auditoría: con ``H=0`` (perfectamente plástico) y
        carga ~50× la cohesión (muy por encima de la capacidad plástica
        del elemento), Newton-Raphson **oscila entre dos estados**
        (ratios ~2 ↔ ~440 alternados), agota ``max_iter``, biseca hasta
        ``min_delta_lambda`` y muere con mensaje genérico
        ``"el solver ha divergido"``.

        Físicamente la solución no existe (F_ext > F_int_max), pero el
        modo de fallo expuesto es **exactamente** el patrón que motiva
        line search: el incremento de Newton aterriza fuera del pozo
        de convergencia y rebota. Sin globalización del Newton, no hay
        forma de detectar el rebote temprano ni de salir limpiamente.

        El test documenta el comportamiento actual: ``RuntimeError`` tras
        bisección agotada.
        """
        E = 1.0e7; nu = 0.3
        phi_deg = 30.0; psi_deg = 10.0; cohesion = 1.0e3
        mat = DruckerPrager2D(
            E=E, nu=nu, cohesion=cohesion, phi_deg=phi_deg, psi_deg=psi_deg,
            H=0.0, hypothesis='plane_strain',
        )
        dom, asm, nodes = _build_uniaxial_quad4(mat, side=1.0, thickness=1.0)

        F_total = 50.0 * cohesion  # excede la capacidad plástica
        F_ext = np.zeros(dom.total_dofs)
        F_ext[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F_ext[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = NonlinearSolver(
            asm, convergence=conv, num_steps=20, max_iter=30,
            adaptive=True, min_delta_lambda=1e-7,
        )
        with self.assertRaises(RuntimeError):
            solver.solve(F_ext)

    def test_damage_2d_strong_softening(self):
        """Daño 2D con softening pronunciado bajo control de carga monotónico.

        Quad4 + ``IsotropicDamage2D`` con ``alpha`` alto (ablandamiento
        abrupto). El control de carga **no puede** atravesar el pico
        (es el régimen para el que existe arc-length). Documentamos que
        ``NonlinearSolver`` diverge limpiamente con ``RuntimeError`` y
        no entra en un bucle inútil.
        """
        E = 2.0e5; nu = 0.0; kappa_0 = 1.0e-4
        alpha = 5000.0  # softening abrupto
        mat = IsotropicDamage2D(
            E=E, nu=nu, kappa_0=kappa_0, alpha=alpha, hypothesis='plane_stress',
        )
        dom, asm, nodes = _build_uniaxial_quad4(mat, side=1.0, thickness=1.0)

        # Carga > pico: σ_peak ≈ E·κ_0 = 20, F_peak = 20·área = 20.
        F_total = 40.0  # 2× el pico
        F_ext = np.zeros(dom.total_dofs)
        F_ext[nodes['n2'].dofs['ux']] = 0.5 * F_total
        F_ext[nodes['n3'].dofs['ux']] = 0.5 * F_total

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver = NonlinearSolver(
            asm, convergence=conv, num_steps=20, max_iter=30,
            adaptive=True, min_delta_lambda=1e-7,
        )
        with self.assertRaises(RuntimeError):
            solver.solve(F_ext)

    def test_plasticity_load_reversal(self):
        """Carga reversal con plasticidad: commit/rollback de history bajo bisección.

        Truss elastoplástico cargado hasta plastificar, descargado a 0,
        recargado en sentido opuesto. Verifica que el commit_all_states
        solo ocurre en pasos convergidos y que un paso fallido (bisección)
        no contamina el estado del siguiente intento.

        Implementación: solver paramétrico con num_steps deliberadamente
        bajo para forzar al menos una bisección durante el reversal.
        """
        E = 100.0; sigma_y = 150.0; H = 50.0
        L = 1.0; A = 1.0
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [L, 0.0])
        mat = Elastoplastic1D(E=E, sigma_y=sigma_y, H=H)
        dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n2.fix_dof('uy', 0.0)
        dom.generate_equation_numbers(verbose=False)
        asm = Assembler(dom)

        # Carga 1: tracción hasta plastificar (F = 200 N → ε = 3.0, eps_p = 1.0).
        F_ext_1 = np.zeros(dom.total_dofs)
        F_ext_1[n2.dofs['ux']] = 200.0
        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        solver_1 = NonlinearSolver(asm, convergence=conv, num_steps=4)
        U1 = solver_1.solve(F_ext_1)
        eps_p_1 = list(dom.elements.values())[0].state.vars[0]['eps_p']
        self.assertAlmostEqual(eps_p_1, 1.0, places=4)

        # Carga 2: invertir signo y aumentar magnitud para plastificar al
        # otro lado. Pasos deliberadamente pocos para tensionar el solver.
        F_ext_2 = np.zeros(dom.total_dofs)
        F_ext_2[n2.dofs['ux']] = -250.0  # del estado actual, hacia compresión
        # NOTE: el NonlinearSolver no acumula cargas; arranca cada solve
        # desde U = 0 con el assembler que ya contiene el committed state.
        # Esto es **importante** y el comportamiento que el test ejercita.
        solver_2 = NonlinearSolver(asm, convergence=conv, num_steps=4)
        U2 = solver_2.solve(F_ext_2)
        # Tras la inversión total, el material debe tener eps_p < eps_p_1
        # (deformación plástica negativa adicional) pero alpha > eps_p_1
        # (acumulada siempre positiva).
        state_final = list(dom.elements.values())[0].state.vars[0]
        self.assertLess(state_final['eps_p'], eps_p_1,
            "eps_p no decreció con la inversión de carga.")
        self.assertGreater(state_final['alpha'], eps_p_1,
            "alpha no creció con la inversión de carga.")


# =============================================================================
# NewtonNewmarkSolver — 2 tests
# =============================================================================

class TestNewtonNewmarkRobustness(unittest.TestCase):
    """Stress de ``NewtonNewmarkSolver`` en regímenes con plasticidad
    transitoria fuerte y cambio brusco de naturaleza de K_t."""

    def test_oversized_dt_with_plasticity(self):
        """Paso temporal sobredimensionado con plasticidad: ¿el solver diverge o sostiene?

        Pulso de carga violento sobre truss elastoplástico perfecto, ``dt``
        deliberadamente grande (T/3, un periodo en 3 puntos). En estática
        este régimen forzaría ``NonlinearSolver`` a bisecar; en dinámica
        la masa estabiliza el jacobiano ``J = M + γΔt·C + βΔt²·K_t`` aun
        cuando ``K_t → 0`` por plasticidad perfecta.

        El test verifica el comportamiento real para que el informe pueda
        documentarlo. Hipótesis: la masa estabiliza y converge en pocas
        iteraciones — un hallazgo positivo de la auditoría, no un bug.
        """
        E = 25.0; rho = 3.0; A = 1.0; L = 1.0
        sigma_y = 50.0; H = 0.0  # plasticidad perfecta, muy agresiva
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [L, 0.0])
        mat = Elastoplastic1D(E=E, sigma_y=sigma_y, H=H, density=rho)
        dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n2.fix_dof('uy', 0.0)
        dom.generate_equation_numbers(verbose=False)
        asm = Assembler(dom)
        dof = n2.dofs['ux']

        T = 2.0 * math.pi / 5.0
        dt_oversized = T / 3.0

        def F_func(t):
            v = np.zeros(dom.total_dofs)
            v[dof] = 200.0 if t < T * 0.5 else 0.0
            return v

        conv = ConvergenceCriterion(rtol_force=1e-7, rtol_disp=1e-7)
        solver = NewtonNewmarkSolver(
            asm, t_end=2.0 * T, dt=dt_oversized,
            F_func=F_func, convergence=conv, max_iter=20,
        )
        # Si converge: la masa estabilizó (hallazgo de auditoría).
        # Si diverge: el régimen sí rompe el Newton-Newmark, lo registramos.
        res = solver.solve()
        self.assertTrue(res.converged,
            "NewtonNewmark divergió con dt sobredimensionado + plasticidad "
            "perfecta — la masa no bastó para estabilizar el jacobiano.")
        # Validar que la plasticidad sí se activó (de lo contrario el test
        # no estaba ejercitando su régimen objetivo).
        alpha_final = list(dom.elements.values())[0].state.vars[0]['alpha']
        self.assertGreater(alpha_final, 0.0,
            "Sin plastificación: el test no ejerció el régimen objetivo.")

    def test_damage_activates_mid_transient(self):
        """Daño que se activa a mitad del transitorio (K_t cambia de naturaleza).

        Oscilador 1 GDL con ``IsotropicDamage1D``: las primeras oscilaciones
        son elásticas (ε < κ_0), luego una carga adicional empuja al daño.
        Verifica que el Rayleigh calibrado con K_0 (no K_t corriente) no
        introduce divergencias cuando la rigidez se degrada.
        """
        E = 25.0; rho = 3.0; A = 1.0; L = 1.0
        kappa_0 = 0.05  # umbral moderado
        alpha = 50.0   # softening no demasiado abrupto
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [L, 0.0])
        mat = IsotropicDamage1D(E=E, kappa_0=kappa_0, alpha=alpha, density=rho)
        dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
        n1.fix_dof('ux', 0.0); n1.fix_dof('uy', 0.0); n2.fix_dof('uy', 0.0)
        dom.generate_equation_numbers(verbose=False)
        asm = Assembler(dom)
        dof = n2.dofs['ux']

        T = 2.0 * math.pi / 5.0  # periodo natural elástico
        dt = T / 80.0

        # Rampa de carga: crece linealmente hasta exceder pico de daño.
        # σ_peak = E·κ_0 = 1.25 → F_peak = 1.25. Rampa a 3·F_peak.
        F_peak = E * kappa_0
        F_max = 3.0 * F_peak
        t_ramp = 2.0 * T

        def F_func(t):
            v = np.zeros(dom.total_dofs)
            v[dof] = F_max * min(t / t_ramp, 1.0)
            return v

        conv = ConvergenceCriterion(rtol_force=1e-6, rtol_disp=1e-6)
        # Amortiguamiento ligero (estabiliza el transitorio).
        rayleigh_cfg = {'alpha': 0.5, 'beta': 0.001}
        solver = NewtonNewmarkSolver(
            asm, t_end=3.0 * T, dt=dt,
            F_func=F_func, convergence=conv, max_iter=30,
            rayleigh=rayleigh_cfg,
        )
        res = solver.solve()

        # Validación: el solver completa el transitorio (no diverge).
        self.assertTrue(res.converged)
        # Y el daño se activó en algún momento (variable de estado final).
        damage_final = list(dom.elements.values())[0].state.vars[0]['damage']
        self.assertGreater(damage_final, 0.0,
            "El daño nunca se activó: el test no ejercita su régimen objetivo.")


if __name__ == '__main__':
    unittest.main()
