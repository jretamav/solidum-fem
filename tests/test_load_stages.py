"""Etapas de carga en los solvers estáticos incrementales (ADR 0026).

Una etapa parte del estado de la anterior (desplazamientos y estado interno
de los elementos) y mantiene constante la carga ya aplicada mientras otra
crece con λ desde 0: la gravedad que se queda y el empuje lateral.

Criterios de aceptación (specs de ``NonlinearSolver``, ``ArcLengthSolver`` e
``IndirectDisplacementSolver``):

- **P-Δ por etapas.** Voladizo corrotacional elástico (20 fibras) con
  ``P = 0.3·P_cr`` y ``H = 1e4 N``: con la gravedad en la etapa 1 y la carga
  lateral en la 2, la flecha de la cabeza coincide con la de una sola etapa
  proporcional (≤ 1e-9) y con la solución de la viga-columna con
  acortamiento axial (error < 5e-4 con 16 elementos, orden cercano a 2;
  ``δ = 2.045488e-4 m``).
- **Interacción N–M.** Columna de 20 fibras de acero perfectamente plástico
  (``Frame2DEuler``, sin P-Δ) con ``n = P/P_y = 0.3`` en la etapa 1 y empuje
  en la 2: la fuerza lateral máxima converge a ``M_p·(1 − n²)/L``.
- **Coherencia elástica.** Modelo elástico lineal: gravedad y después empuje
  da los desplazamientos de sumar los dos casos (≤ 1e-10); la reacción
  vertical suma ``P`` al terminar la gravedad y durante todo el empuje.
- **Continuación igual a una corrida.** Empujar hasta 0.03 m y continuar
  hasta 0.06 m da la misma curva que empujar hasta 0.06 m (≤ 1e-9 en la
  fuerza en los desplazamientos comunes).
- **Equilibrio inicial consistente.** Tras una etapa que plastifica, el
  residuo al empezar la siguiente, ``F_int(U_fin) − F_c``, es nulo; empezar en
  ``U = 0`` con el estado plástico se avisa.

Además: valores impuestos (apoyos con desplazamiento) y restricciones
multipunto, red de seguridad, cargas distribuidas constantes en las fuerzas
internas, y errores claros donde no hay etapas.
"""
import logging
import math
import unittest

import numpy as np

import solidum
from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler, Frame2DEulerCorot
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.convergence import ConvergenceCriterion
from solidum.math.solvers import (
    ArcLengthSolver,
    IndirectDisplacementSolver,
    LinearSolver,
    NonlinearSolver,
)
from solidum.math.solvers.diagnostics import MechanismError, Termination
from solidum.sections import FiberSection, RectangularPatch

E, L, B, H, NF, SIGMA_Y = 2.0e11, 3.0, 0.30, 0.50, 20, 250.0e6
I_FIB = B * H ** 3 / 12.0 * (1.0 - 1.0 / NF ** 2)      # inercia de 20 fibras iguales
A_SEC = B * H
P_CR = math.pi ** 2 * E * I_FIB / (4.0 * L ** 2)
P = 0.3 * P_CR
H_LAT = 1.0e4


def tight():
    return ConvergenceCriterion(rtol_force=1e-12, rtol_disp=1e-12)


def column(n, element=Frame2DEuler, fiber=None, *, density=None):
    """Voladizo vertical de ``n`` elementos con sección de 20 fibras. Devuelve
    dominio, nodo de la cabeza y las cargas de gravedad (``−P`` en ``uy``) y
    lateral (``H_LAT`` en ``ux``) en la cabeza."""
    fiber = fiber if fiber is not None else Elastic1D(E=E, density=density)
    section = FiberSection(patches=[RectangularPatch(fiber, y=(-H / 2, H / 2),
                                                     z=(-B / 2, B / 2), ny=NF, nz=1)])
    dom = Domain()
    nodes = [dom.add_node(i + 1, [0.0, i * L / n]) for i in range(n + 1)]
    for k in range(n):
        dom.add_element(element(k + 1, [nodes[k], nodes[k + 1]], section=section, n_ip=3))
    for dof in ("ux", "uy", "rz"):
        nodes[0].fix_dof(dof, 0.0)
    dom.generate_equation_numbers()
    top = nodes[-1]
    F_g = np.zeros(dom.total_dofs)
    F_g[top.dofs["uy"]] = -P
    F_h = np.zeros(dom.total_dofs)
    F_h[top.dofs["ux"]] = H_LAT
    return dom, top, F_g, F_h


def staged(dom, F_g, F_h, second, *, first=None, **run_kw):
    """Gravedad con ``NonlinearSolver`` y después ``second(asm)`` con la carga
    lateral, como etapa que continúa la primera. Devuelve los dos resultados
    y el segundo solver."""
    asm = Assembler(dom)
    s1 = first(asm) if first else NonlinearSolver(asm, num_steps=4, convergence=tight())
    r1 = solidum.run(dom, assembler=asm, solver=s1, F_applied=F_g)
    s2 = second(asm)
    r2 = solidum.run(dom, assembler=asm, solver=s2, F_applied=F_h, previous=r1, **run_kw)
    return r1, r2, s2


def beam_column_deflection(n_ip_unused=None):
    """Flecha de la viga-columna en voladizo con acortamiento axial (anexo de
    la spec): ``δ' = H/(P k')·(tan k'L' − k'L')``, con ``L' = L(1−ε)``,
    ``EI' = EI(1−ε)``, ``k' = √(P/EI')`` y ``ε = P/(EA)``."""
    eps = P / (E * A_SEC)
    L_s, EI_s = L * (1.0 - eps), E * I_FIB * (1.0 - eps)
    k = math.sqrt(P / EI_s)
    return H_LAT / (P * k) * (math.tan(k * L_s) - k * L_s)


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


class TestPDeltaPorEtapas(_Quiet):

    def test_etapas_igual_a_una_sola_etapa_proporcional(self):
        for second in (
            lambda asm: NonlinearSolver(asm, num_steps=4, convergence=tight()),
            lambda asm: ArcLengthSolver(asm, max_lambda=1.0, convergence=tight()),
        ):
            dom, top, F_g, F_h = column(8, Frame2DEulerCorot)
            r1, r2, s2 = staged(dom, F_g, F_h, second)
            dom_p, top_p, F_gp, F_hp = column(8, Frame2DEulerCorot)
            asm_p = Assembler(dom_p)
            r_p = solidum.run(dom_p, assembler=asm_p, F_applied=F_gp + F_hp,
                              solver=NonlinearSolver(asm_p, num_steps=4, convergence=tight()))
            with self.subTest(solver=type(s2).__name__):
                self.assertLessEqual(abs(r2.U[top.dofs["ux"]] / r_p.U[top_p.dofs["ux"]] - 1.0), 1e-9)
                self.assertAlmostEqual(r2.load_factor, 1.0, delta=1e-12)
                np.testing.assert_array_equal(r2.F_constant, r1.F_applied)
                np.testing.assert_allclose(r2.F_variable, F_h, rtol=0, atol=1e-12 * H_LAT)

    def test_contra_la_viga_columna_con_acortamiento(self):
        delta = beam_column_deflection()
        self.assertAlmostEqual(delta, 2.045488e-4, delta=5e-11)
        errors = {}
        for n in (8, 16):
            dom, top, F_g, F_h = column(n, Frame2DEulerCorot)
            _, r2, _ = staged(dom, F_g, F_h,
                              lambda asm: NonlinearSolver(asm, num_steps=4, convergence=tight()))
            errors[n] = abs(r2.U[top.dofs["ux"]] / delta - 1.0)
        self.assertLess(errors[16], 5e-4)
        self.assertAlmostEqual(math.log2(errors[8] / errors[16]), 2.0, delta=0.1)


class TestInteraccionNM(_Quiet):

    def test_carga_limite_con_carga_axial_constante(self):
        """``n = 0.3``: la banda central plástica es de 6 fibras exactas, así
        que la sección discreta da ``M_p·(1 − n²)`` sin error. El elemento en
        desplazamientos sobrepasa la carga límite con un error de malla de
        primer orden; con 16 elementos el exceso es del 3 % (2 veces el de
        ``n = 0``) y su extrapolación de Richardson da la carga límite."""
        n_ax = 0.3
        P_ax = n_ax * SIGMA_Y * B * H
        F_lim = SIGMA_Y * B * H ** 2 / 4.0 * (1.0 - n_ax ** 2) / L
        excess = {}
        for n in (4, 8, 16):
            steel = Elastoplastic1D(E=E, sigma_y=SIGMA_Y, H=0.0)
            dom, top, _, _ = column(n, Frame2DEuler, steel)
            F_g = np.zeros(dom.total_dofs)
            F_g[top.dofs["uy"]] = -P_ax
            F_h = np.zeros(dom.total_dofs)
            F_h[top.dofs["ux"]] = 1.0e6
            lams = []
            r1, r2, s2 = staged(
                dom, F_g, F_h,
                lambda asm: IndirectDisplacementSolver(
                    asm, control=[(top.id, "ux", 1.0)], max_lambda=1e9, initial_dl=1.5e-3,
                    control_target=0.09, max_steps=200),
                first=lambda asm: NonlinearSolver(asm, num_steps=3),
                step_callback=lambda k, U, lam: lams.append(lam))
            self.assertEqual(s2.termination, "control_target")
            self.assertLess(s2.initial_residual, 1e-12)
            # La carga axial sigue completa durante el empuje.
            self.assertAlmostEqual(r2.reactions_by_node[1]["uy"] / P_ax, 1.0, delta=1e-9)
            excess[n] = max(lams) * 1.0e6 / F_lim - 1.0
        self.assertTrue(excess[4] > excess[8] > excess[16] > 0.0)
        for coarse, fine in ((4, 8), (8, 16)):
            self.assertAlmostEqual(excess[coarse] / excess[fine], 2.0, delta=0.15)
        # Extrapolación de Richardson (primer orden): la carga límite.
        self.assertLess(abs(2.0 * excess[16] - excess[8]), 3e-3)


class TestCoherenciaElastica(_Quiet):

    def test_superposicion_y_reaccion_vertical(self):
        dom, top, F_g, F_h = column(4, Frame2DEuler)
        asm_lin = Assembler(dom)
        U_g = LinearSolver(asm_lin).solve(F_g)
        U_h = LinearSolver(asm_lin).solve(F_h)
        dom2, top2, F_g2, F_h2 = column(4, Frame2DEuler)
        asm2 = Assembler(dom2)
        asm2.assemble_system()
        K = asm2.K_global.copy()
        base_uy = dom2.nodes[1].dofs["uy"]
        reactions = []
        r1, r2, _ = staged(
            dom2, F_g2, F_h2,
            lambda asm: IndirectDisplacementSolver(asm, control=[(top2.id, "ux", 1.0)],
                                                   max_lambda=1.0, initial_dl=2.0e-5,
                                                   convergence=tight()),
            step_callback=lambda k, U, lam: reactions.append(float((K @ U)[base_uy])))
        np.testing.assert_allclose(r2.U, U_g + U_h, rtol=0, atol=1e-10 * np.max(np.abs(U_g + U_h)))
        # Flecha H·L³/(3EI) = 1.4436e-4 m y acortamiento P·L/(EA) = 5.1276e-3 m.
        self.assertAlmostEqual(r2.U[top2.dofs["ux"]] / (H_LAT * L ** 3 / (3.0 * E * I_FIB)), 1.0,
                               delta=1e-10)
        self.assertAlmostEqual(r2.U[top2.dofs["uy"]] / (-P * L / (E * A_SEC)), 1.0, delta=1e-10)
        self.assertAlmostEqual(r2.U[top2.dofs["uy"]] / -5.1276e-3, 1.0, delta=1e-5)
        self.assertAlmostEqual(r1.reactions_by_node[1]["uy"] / P, 1.0, delta=1e-10)
        self.assertGreater(len(reactions), 3)
        np.testing.assert_allclose(reactions, P, rtol=1e-10)


def plastic_cantilever(n=4):
    steel = Elastoplastic1D(E=E, sigma_y=SIGMA_Y, H=2.0e9)
    dom, top, _, F_h = column(n, Frame2DEuler, steel)
    F_h = F_h / H_LAT * 1.0e6
    return dom, top, F_h


def _push(asm, top, target, conv=None):
    return IndirectDisplacementSolver(asm, control=[(top.id, "ux", 1.0)], max_lambda=1e9,
                                      initial_dl=1.5e-3, control_target=target,
                                      convergence=conv or tight())


class TestContinuacion(_Quiet):

    def test_continuar_es_igual_a_una_corrida(self):
        dom, top, F = plastic_cantilever()
        asm = Assembler(dom)
        direct = []
        solidum.run(dom, assembler=asm, solver=_push(asm, top, 0.06), F_applied=F,
                    step_callback=lambda k, U, lam: direct.append((U[top.dofs["ux"]], lam)))

        dom2, top2, F2 = plastic_cantilever()
        asm2 = Assembler(dom2)
        staged_path = []
        r1 = solidum.run(dom2, assembler=asm2, solver=_push(asm2, top2, 0.03), F_applied=F2,
                         step_callback=lambda k, U, lam: staged_path.append((U[top2.dofs["ux"]], lam)))
        lam1 = r1.load_factor
        s2 = _push(asm2, top2, 0.06)
        r2 = solidum.run(dom2, assembler=asm2, solver=s2, F_applied=F2, previous=r1,
                         step_callback=lambda k, U, lam: staged_path.append((U[top2.dofs["ux"]], lam1 + lam)))
        direct, staged_path = np.array(direct), np.array(staged_path)
        self.assertEqual(direct.shape, staged_path.shape)
        np.testing.assert_allclose(staged_path[:, 0], direct[:, 0], rtol=1e-12)
        np.testing.assert_allclose(staged_path[:, 1], direct[:, 1], rtol=1e-9)
        self.assertLess(s2.initial_residual, 1e-12)
        self.assertAlmostEqual(r2.U[top2.dofs["ux"]], 0.06, delta=1e-15)


    def test_descargar_hasta_un_objetivo_menor(self):
        """El objetivo es absoluto y fija el sentido: tras empujar hasta
        0.06 m, una etapa con objetivo 0.045 m descarga (λ < 0), llega exacto
        y, como la descarga es elástica, con la rigidez elástica 3·E·I/L³."""
        dom, top, F = plastic_cantilever()
        asm = Assembler(dom)
        r1 = solidum.run(dom, assembler=asm, solver=_push(asm, top, 0.06), F_applied=F)
        lams = []
        s2 = _push(asm, top, 0.045)
        r2 = solidum.run(dom, assembler=asm, solver=s2, F_applied=F, previous=r1,
                         step_callback=lambda k, U, lam: lams.append(lam))
        self.assertEqual(s2.termination, Termination.CONTROL_TARGET)
        self.assertTrue(r2.converged)
        self.assertAlmostEqual(r2.U[top.dofs["ux"]], 0.045, delta=1e-15)
        self.assertTrue(np.all(np.diff([0.0] + lams) < 0.0))
        k_el = 3.0 * E * I_FIB / L ** 3
        self.assertAlmostEqual(s2.lambda_final, -k_el * 0.015 / 1.0e6,
                               delta=1e-8 * k_el * 0.015 / 1.0e6)

    def test_objetivo_igual_al_punto_de_partida(self):
        dom, top, F = plastic_cantilever()
        asm = Assembler(dom)
        r1 = solidum.run(dom, assembler=asm, solver=_push(asm, top, 0.03), F_applied=F)
        s2 = _push(asm, top, 0.03)
        logging.disable(logging.NOTSET)
        with self.assertLogs("solidum", level="WARNING") as logs:
            r2 = solidum.run(dom, assembler=asm, solver=s2, F_applied=F, previous=r1)
        self.assertTrue(any("no hay nada que trazar" in m for m in logs.output))
        self.assertEqual(s2.termination, Termination.CONTROL_TARGET)
        self.assertEqual(s2.steps_done, 0)
        np.testing.assert_array_equal(r2.U, r1.U)


class TestEquilibrioInicial(_Quiet):

    def test_tras_una_etapa_plastica(self):
        dom, top, F = plastic_cantilever()
        asm = Assembler(dom)
        r1 = solidum.run(dom, assembler=asm, solver=_push(asm, top, 0.04), F_applied=F)
        # La etapa plastificó: la rigidez secante queda bajo el 80 % de la
        # elástica 3EI/L³ (medido: 60 %).
        k_secant = r1.load_factor * 1.0e6 / 0.04
        self.assertLess(k_secant, 0.8 * 3.0 * E * I_FIB / L ** 3)
        _, F_int = asm.assemble_non_linear_system(r1.U)
        T = asm.constraint_set.build(dom.total_dofs)[0]
        residual = np.linalg.norm(T.T @ (F_int - r1.F_applied)) / np.linalg.norm(r1.F_applied)
        self.assertLess(residual, 1e-10)
        s2 = _push(asm, top, 0.05)
        solidum.run(dom, assembler=asm, solver=s2, F_applied=F, previous=r1)
        self.assertLess(s2.initial_residual, 1e-10)

    def test_empezar_en_cero_con_el_estado_plastico_se_avisa(self):
        """Lo contrario: una corrida nueva sobre el mismo modelo empieza en
        ``U = 0`` con el estado plástico comprometido, que no está en
        equilibrio ahí (fuerza interna ~7 veces la de primera fluencia)."""
        dom, top, F = plastic_cantilever()
        asm = Assembler(dom)
        solidum.run(dom, assembler=asm, solver=_push(asm, top, 0.04), F_applied=F)
        logging.disable(logging.NOTSET)
        with self.assertLogs("solidum", level="WARNING") as logs:
            s2 = _push(asm, top, 0.05)
            try:
                s2.solve(F)
            except Exception:
                pass
        self.assertGreater(s2.initial_residual, 1e-3)
        self.assertTrue([m for m in logs.output if "estado interno de un análisis anterior" in m])


class TestCompatibilidad(_Quiet):

    def test_apoyo_con_desplazamiento_impuesto(self):
        """Un asentamiento impuesto en la etapa 1 se queda en la 2; uno que
        cambia entre etapas crece con λ desde el valor de partida."""
        dom, top, F_g, F_h = column(4, Frame2DEuler)
        base = dom.nodes[1]
        base.fix_dof("uy", -0.01)
        asm = Assembler(dom)
        r1 = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=2), F_applied=F_g)
        self.assertAlmostEqual(r1.U[base.dofs["uy"]], -0.01, delta=1e-15)
        seen = []
        r2 = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=4), F_applied=F_h,
                         previous=r1, step_callback=lambda k, U, lam: seen.append(U[base.dofs["uy"]]))
        np.testing.assert_allclose(seen, -0.01, rtol=0, atol=1e-15)
        self.assertAlmostEqual(r2.U[base.dofs["uy"]], -0.01, delta=1e-15)
        base.fix_dof("uy", -0.02)            # la etapa 3 cambia el asentamiento
        path = []
        solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=4, adaptive=False),
                    F_applied=np.zeros(dom.total_dofs), previous=r2,
                    step_callback=lambda k, U, lam: path.append((lam, U[base.dofs["uy"]])))
        for lam, uy in path:
            self.assertAlmostEqual(uy, -0.01 - 0.01 * lam, delta=1e-14)

    def test_restriccion_multipunto_igual_al_nodo_fundido(self):
        """Columna partida en su nodo intermedio y atada (``u_s = u_m`` en
        ``ux``, ``uy``, ``rz``): las etapas dan lo mismo que la columna entera."""
        dom_m, top_m, F_gm, F_hm = column(4, Frame2DEulerCorot)
        _, r_merged, _ = staged(dom_m, F_gm, F_hm,
                                lambda asm: NonlinearSolver(asm, num_steps=4, convergence=tight()))

        dom = Domain()
        fiber = Elastic1D(E=E)
        section = FiberSection(patches=[RectangularPatch(fiber, y=(-H / 2, H / 2),
                                                         z=(-B / 2, B / 2), ny=NF, nz=1)])
        ys = [i * L / 4 for i in range(5)]
        nodes = [dom.add_node(i + 1, [0.0, y]) for i, y in enumerate(ys)]
        dup = dom.add_node(100, [0.0, ys[2]])
        chain = [nodes[0], nodes[1], nodes[2], dup, nodes[3], nodes[4]]
        eid = 1
        for a, b in ((0, 1), (1, 2), (3, 4), (4, 5)):
            dom.add_element(Frame2DEulerCorot(eid, [chain[a], chain[b]], section=section, n_ip=3))
            eid += 1
        for dof in ("ux", "uy", "rz"):
            nodes[0].fix_dof(dof, 0.0)
            dom.add_linear_constraint((100, dof), [(nodes[2].id, dof)], [1.0])
        dom.generate_equation_numbers()
        F_g = np.zeros(dom.total_dofs)
        F_g[nodes[4].dofs["uy"]] = -P
        F_h = np.zeros(dom.total_dofs)
        F_h[nodes[4].dofs["ux"]] = H_LAT
        _, r_tied, _ = staged(dom, F_g, F_h,
                              lambda asm: NonlinearSolver(asm, num_steps=4, convergence=tight()))
        for dof in ("ux", "uy", "rz"):
            self.assertAlmostEqual(r_tied.U[nodes[4].dofs[dof]] / r_merged.U[top_m.dofs[dof]], 1.0,
                                   delta=1e-9)

    def test_red_de_seguridad_sigue_activa(self):
        dom, top, F_g, F_h = column(2, Frame2DEuler)
        asm = Assembler(dom)
        r1 = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=2), F_applied=F_g)
        for dof in ("ux", "uy", "rz"):
            dom.nodes[1].boundary_conditions.pop(dof)
        with self.assertRaises(MechanismError):
            solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=2),
                        F_applied=F_h, previous=r1)


class TestCargaDistribuidaConstante(_Quiet):
    """Las fuerzas internas restan las cargas equivalentes con su factor: la
    constante completa y la que crece por λ (columna elástica lineal, así
    que la etapa final es la superposición)."""

    def _linear_reference(self, F_lateral_builder):
        """Una sola etapa lineal con las cargas equivalentes de las dos cargas
        sumadas a mano: el registro del ensamblador sólo guarda las de la
        última carga de cuerpo ensamblada (deuda #35)."""
        dom, top, _, _ = column(4, Frame2DEuler, density=7850.0)
        asm = Assembler(dom)
        F_sw = asm.assemble_self_weight(np.array([0.0, -9.81]))
        eq = {eid: np.array(f) for eid, f in asm.equivalent_loads.items()}
        version = asm.equivalent_loads_version
        F_lat = F_lateral_builder(asm, dom, top)
        if asm.equivalent_loads_version != version:
            for eid, f in asm.equivalent_loads.items():
                eq[eid] = eq.get(eid, 0.0) + f
        U = LinearSolver(asm).solve(F_sw + F_lat)
        return solidum.results.build_solve_result(dom, asm, U, F_sw + F_lat, equivalent_loads=eq)

    def _staged(self, F_lateral_builder):
        dom, top, _, _ = column(4, Frame2DEuler, density=7850.0)
        asm = Assembler(dom)
        F_sw = asm.assemble_self_weight(np.array([0.0, -9.81]))
        r1 = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=2, convergence=tight()),
                         F_applied=F_sw)
        F_lat = F_lateral_builder(asm, dom, top)
        return solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=3, convergence=tight()),
                           F_applied=F_lat, previous=r1)

    def _compare(self, builder):
        ref, st = self._linear_reference(builder), self._staged(builder)
        for eid, ef in ref.element_forces.items():
            for name in ("N", "V", "M"):
                a = np.asarray(ef.components[name], dtype=float)
                b = np.asarray(st.element_forces[eid].components[name], dtype=float)
                np.testing.assert_allclose(b, a, rtol=1e-8, atol=1e-8 * np.max(np.abs(a)) + 1e-6)

    def test_peso_propio_constante_y_carga_puntual_que_crece(self):
        def point(asm, dom, top):
            F = np.zeros(dom.total_dofs)
            F[top.dofs["ux"]] = 2.0e4
            return F
        self._compare(point)

    def test_peso_propio_constante_y_carga_distribuida_que_crece(self):
        """La carga que crece también es distribuida (fuerza de cuerpo
        lateral): sus cargas equivalentes, ensambladas después de la etapa 1,
        entran con λ, y las del peso propio, con 1."""
        self._compare(lambda asm, dom, top: asm.assemble_body_load(np.array([3.0e4, 0.0])))


class TestErroresClaros(_Quiet):

    def test_linear_solver_no_admite_etapas(self):
        dom, top, F_g, F_h = column(2, Frame2DEuler)
        asm = Assembler(dom)
        r1 = solidum.run(dom, assembler=asm, solver=LinearSolver(asm), F_applied=F_g)
        with self.assertRaisesRegex(ValueError, "no admite etapas"):
            solidum.run(dom, assembler=asm, solver=LinearSolver(asm), F_applied=F_h, previous=r1)

    def test_tamano_del_estado_inicial(self):
        dom, top, F_g, F_h = column(2, Frame2DEuler)
        with self.assertRaisesRegex(ValueError, "u0 tiene 3 componentes"):
            NonlinearSolver(Assembler(dom)).solve(F_h, u0=np.zeros(3))


if __name__ == "__main__":
    unittest.main()
