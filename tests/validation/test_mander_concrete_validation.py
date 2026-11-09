"""Validación de ``ManderConcrete1D`` (spec ``ManderConcrete1D``,
``acceptance.validation``).

Sección rectangular de concreto reforzado de 0.30 × 0.50 m con recubrimiento
de 4 cm: núcleo confinado (función de confinamiento de Mander con estribos
rectangulares), recubrimiento sin confinar y 3 + 3 barras de 20 mm de acero
elastoplástico con endurecimiento.

1. **Momento–curvatura** contra una integración independiente de la sección
   (4000 franjas, envolventes escritas aquí, búsqueda de la raíz de ``N``
   con bisección), sin carga axial y con ``N = −0.1·f'c·A_g``.
2. **Voladizo** de 16 ``Frame2DEuler`` hasta el pico: fuerza contra la
   integración de la curvatura de la solución independiente.
3. **Objetividad de malla**: barra a tracción con ``G_f`` (la curva
   fuerza–desplazamiento y la energía disipada no dependen de la malla y
   coinciden con la solución analítica de la banda de fisuración); voladizo
   después del pico con ``G_fc`` (8 y 16 elementos dentro del 5 %).
4. **Ciclos** de desplazamiento de amplitud creciente por etapas.
5. **Transitorio**: pulso y vibración libre; la energía disipada por el
   material, trabajo interno menos energía almacenada en las fibras, es no
   negativa y no decrece.
6. **YAML** con el confinamiento calculado desde el propio YAML: idéntico a
   Python.
"""
import logging
import math
import textwrap
import unittest

import numpy as np

import solidum
from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler
from solidum.elements.truss import Truss2D
from solidum.entry import run_yaml
from solidum.materials.concrete_1d import ManderConcrete1D
from solidum.materials.concrete_confinement import (
    mander_confined_concrete,
    mander_elastic_modulus,
)
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.solvers import (
    IndirectDisplacementSolver,
    NewtonNewmarkSolver,
    NonlinearSolver,
)
from solidum.sections import BarLayer, FiberSection, RectangularPatch
from solidum.utils.yaml_parser import YamlParser

# --- Datos de la sección (SI) ------------------------------------------------
B, H, C = 0.30, 0.50, 0.04
YB = H / 2 - C                         # eje de las barras (borde del núcleo)
FC0, FT = 30.0e6, 3.0e6
EC = mander_elastic_modulus(FC0)
EPS_T = FT / EC
EPS_TU = 10.0 * EPS_T
ES, FY, HS = 200.0e9, 420.0e6, 2.0e9
AB = math.pi * 0.020 ** 2 / 4
CONF = dict(fc0=FC0, f_yh=420.0e6, eps_su=0.09, k_e=0.75, rho_x=0.01, rho_y=0.01)
CORE = mander_confined_concrete(**CONF)
COVER = dict(fc=FC0, eps_c=0.002, eps_cu=0.004, eps_sp=0.006, E=EC)
L = 3.0


def materials(*, density=False):
    rho_c = dict(density=2400.0) if density else {}
    rho_s = dict(density=7850.0) if density else {}
    core = ManderConcrete1D(confinement=CONF, ft=FT, eps_tu=EPS_TU, **rho_c)
    cover = ManderConcrete1D(**COVER, ft=FT, eps_tu=EPS_TU, **rho_c)
    steel = Elastoplastic1D(E=ES, sigma_y=FY, H=HS, **rho_s)
    return core, cover, steel


def make_section(ny=60, *, displaces=True, density=False):
    core, cover, steel = materials(density=density)
    disp = dict(displaces=core) if displaces else {}
    return FiberSection(patches=[
        RectangularPatch(core, y=(-YB, YB), z=(-B / 2 + C, B / 2 - C), ny=ny, nz=1),
        RectangularPatch(cover, y=(YB, H / 2), z=(-B / 2, B / 2), ny=4, nz=1),
        RectangularPatch(cover, y=(-H / 2, -YB), z=(-B / 2, B / 2), ny=4, nz=1),
        RectangularPatch(cover, y=(-YB, YB), z=(-B / 2, -B / 2 + C), ny=ny, nz=1),
        RectangularPatch(cover, y=(-YB, YB), z=(B / 2 - C, B / 2), ny=ny, nz=1),
        BarLayer(steel, start=(YB, -0.1), end=(YB, 0.1), n=3, area=AB, **disp),
        BarLayer(steel, start=(-YB, -0.1), end=(-YB, 0.1), n=3, area=AB, **disp),
    ])


# --- Integración independiente de la sección ------------------------------------

def _popovics(e, fc, ec, E):
    r = E / (E - fc / ec)
    x = np.maximum(e, 0.0) / ec
    return fc * x * r / (r - 1.0 + x ** r)


def _concrete(eps, fc, ec, ecu, esp, E):
    """Envolvente monótona (σ con signo), escrita aquí sin el código de Solidum."""
    e = -eps
    comp = np.where(e <= ecu, _popovics(e, fc, ec, E), 0.0)
    if esp > ecu:
        s_cu = _popovics(np.array(ecu), fc, ec, E)
        comp = np.where((e > ecu) & (e <= esp), s_cu * (esp - e) / (esp - ecu), comp)
    tens = np.where(eps <= EPS_T, E * eps,
                    np.where(eps < EPS_TU, FT * (EPS_TU - eps) / (EPS_TU - EPS_T), 0.0))
    return np.where(eps < 0.0, -comp, tens)


def _steel(eps):
    ey = FY / ES
    Et = ES * HS / (ES + HS)
    return np.where(np.abs(eps) <= ey, ES * eps, np.sign(eps) * (FY + Et * (np.abs(eps) - ey)))


_YS = np.linspace(-H / 2, H / 2, 4001)
_YM = 0.5 * (_YS[1:] + _YS[:-1])
_DY = np.diff(_YS)
_IN_CORE = np.abs(_YM) <= YB
_W_CORE = np.where(_IN_CORE, B - 2 * C, 0.0) * _DY
_W_COVER = np.where(_IN_CORE, 2 * C, B) * _DY


def independent_NM(eps0, kappa):
    eps = eps0 - _YM * kappa
    s = (_W_CORE * _concrete(eps, CORE["fc"], CORE["eps_c"], CORE["eps_cu"], CORE["eps_cu"],
                             CORE["E"])
         + _W_COVER * _concrete(eps, FC0, 0.002, 0.004, 0.006, EC))
    N, M = s.sum(), -(_YM * s).sum()
    for yb in (YB, -YB):
        e = eps0 - yb * kappa
        f = 3 * AB * (_steel(e) - _concrete(e, CORE["fc"], CORE["eps_c"], CORE["eps_cu"],
                                            CORE["eps_cu"], CORE["E"]))
        N, M = N + f, M - yb * f
    return float(N), float(M)


def root_near(f, x0, step=1e-4, tol=1e-12):
    """Raíz de ``f`` más cercana a ``x0``: búsqueda hacia fuera y bisección
    (``N(ε₀)`` no es monótona cuando hay fibras en ablandamiento)."""
    f0 = f(x0)
    if abs(f0) < tol:
        return x0
    for _ in range(200):
        for x1 in (x0 + step, x0 - step):
            f1 = f(x1)
            if f1 * f0 <= 0.0:
                a, b, fa = x0, x1, f0
                for _ in range(200):
                    m = 0.5 * (a + b)
                    fm = f(m)
                    if abs(fm) < tol or abs(b - a) < 1e-16:
                        return m
                    if fm * fa <= 0.0:
                        b = m
                    else:
                        a, fa = m, fm
                return 0.5 * (a + b)
        step *= 1.5
    raise RuntimeError("sin raíz")


def independent_moment_curvature(N, kappas):
    eps0, out = 0.0, []
    for kap in kappas:
        eps0 = root_near(lambda e: (independent_NM(e, kap)[0] - N) / (FC0 * B * H), eps0)
        out.append(independent_NM(eps0, kap)[1])
    return np.array(out)


def solidum_moment_curvature(section, N, kappas):
    """La misma búsqueda de raíz sobre ``FiberSection.response_2d`` con los
    estados de las fibras comprometidos curvatura a curvatura."""
    states = section.initial_states()
    eps0, out = 0.0, []
    for kap in kappas:
        eps0 = root_near(lambda e: (section.response_2d([e, kap], states)[0][0] - N)
                         / (FC0 * B * H), eps0)
        s, _, states = section.response_2d([eps0, kap], states)
        out.append(s[1])
    return np.array(out)


def cantilever(n, section, *, N_axial=0.0):
    dom = Domain()
    nodes = [dom.add_node(i + 1, [0.0, i * L / n]) for i in range(n + 1)]
    for k in range(n):
        dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], section=section, n_ip=3))
    for dof in ("ux", "uy", "rz"):
        nodes[0].fix_dof(dof, 0.0)
    dom.generate_equation_numbers()
    top = nodes[-1]
    F_g = np.zeros(dom.total_dofs)
    F_g[top.dofs["uy"]] = N_axial
    F_h = np.zeros(dom.total_dofs)
    F_h[top.dofs["ux"]] = 1.0e5
    return dom, top, F_g, F_h


def push(dom, top, F_g, F_h, targets, initial_dl=2e-3):
    """Etapa de gravedad (si hay carga axial) y una etapa de desplazamiento
    controlado por cada objetivo. Devuelve los puntos (δ, F) y los solvers."""
    asm = Assembler(dom)
    previous = None
    if np.any(F_g):
        previous = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=4),
                               F_applied=F_g)
    path, solvers, lam0 = [], [], 0.0
    for target in targets:
        s = IndirectDisplacementSolver(asm, control=[(top.id, "ux", 1.0)], max_lambda=1e9,
                                       initial_dl=initial_dl, control_target=target,
                                       max_steps=3000)
        base = lam0
        r = solidum.run(dom, assembler=asm, solver=s, F_applied=F_h, previous=previous,
                        step_callback=lambda k, U, lam, base=base:
                        path.append((U[top.dofs["ux"]], (base + lam) * 1e5)))
        lam0 = base + r.load_factor
        previous = r
        solvers.append(s)
    return np.array(path), solvers


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


class TestMomentoCurvatura(_Quiet):

    def test_contra_la_integracion_independiente(self):
        section = make_section()
        for N, kap_max in ((0.0, 0.14), (-0.1 * FC0 * B * H, 0.25)):
            with self.subTest(N=N):
                kappas = np.linspace(kap_max / 140, kap_max, 140)
                M_s = solidum_moment_curvature(section, N, kappas)
                M_i = independent_moment_curvature(N, kappas[::7])
                err = np.abs(M_s[::7] - M_i) / np.abs(M_i).max()
                self.assertLess(err.max(), 2e-2)
                self.assertAlmostEqual(M_s.max() / M_i.max(), 1.0, delta=2e-2)


class TestVoladizo(_Quiet):

    def test_hasta_el_pico(self):
        """Flecha ``δ = ∫ κ(F·(L − x))·(L − x) dx`` con la curva momento–
        curvatura independiente, hasta el primer máximo del momento; se
        compara la fuerza a desplazamiento dado, que está bien condicionada
        también en la meseta posterior a la fluencia (allí ``δ(F)`` no lo
        está: con carga axial la curva es casi plana)."""
        for N in (0.0, -0.1 * FC0 * B * H):
            with self.subTest(N=N):
                kap = np.linspace(1e-4, 0.08, 400)
                M_i = independent_moment_curvature(N, kap)
                drops = np.flatnonzero(np.diff(M_i) < 0)
                ipk = int(drops[0]) if drops.size else len(M_i) - 1
                M_mon, k_mon = np.r_[0.0, M_i[:ipk + 1]], np.r_[0.0, kap[:ipk + 1]]
                x = np.linspace(0.0, L, 3001)

                def delta(F):
                    k = np.interp(F * (L - x), M_mon, k_mon)
                    return np.trapezoid(k * (L - x), x)

                F_pk = M_mon[-1] / L
                # 16 elementos: con 8, el error de discretización del elemento en
                # desplazamientos en la meseta posterior a la fluencia llega a 2.2 %
                # con carga axial; converge al refinar (no es del material).
                dom, top, F_g, F_h = cantilever(16, make_section(), N_axial=N)
                path, solvers = push(dom, top, F_g, F_h, [1.3 * delta(0.98 * F_pk)],
                                     initial_dl=5e-4)
                self.assertEqual(solvers[-1].termination, "control_target")
                F_grid = np.linspace(0.0, 0.98 * F_pk, 400)
                d_grid = np.array([delta(F) for F in F_grid])
                checked = 0
                for d, F in path:
                    if 0.02 * d_grid[-1] < d < d_grid[-1]:
                        F_ind = np.interp(d, d_grid, F_grid)
                        self.assertAlmostEqual(F / F_ind, 1.0, delta=2e-2, msg=f"δ = {d:.4g}")
                        checked += 1
                self.assertGreater(checked, 5)


class TestObjetividadBarra(_Quiet):
    """Barra de concreto simple de 0.5 m a tracción con el elemento central
    un 1 % más débil. Con ``G_f``, la banda de fisuración da la misma curva con
    cualquier malla: ``δ = σL/E`` antes del pico y
    ``δ = 2G_f/f_t − σ·(2G_f/f_t² − L/E)`` en el ablandamiento, y la energía
    disipada ``A·G_f·(1 − σ/f_t)``."""

    E, A, LB, GF = 27.4e9, 0.01, 0.5, 140.0
    FTW = 0.99 * FT

    def bar(self, n, **tension):
        strong = ManderConcrete1D(fc=FC0, eps_c=0.002, E=self.E, eps_cu=0.004, ft=FT, **tension)
        weak = ManderConcrete1D(fc=FC0, eps_c=0.002, E=self.E, eps_cu=0.004, ft=self.FTW,
                                **tension)
        dom = Domain()
        nodes = [dom.add_node(i + 1, [i * self.LB / n, 0.0]) for i in range(n + 1)]
        for k in range(n):
            dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]],
                                    weak if k == n // 2 else strong, A=self.A))
        nodes[0].fix_dof("ux", 0.0)
        for nd in nodes:
            nd.fix_dof("uy", 0.0)
        dom.generate_equation_numbers()
        end = nodes[-1]
        F = np.zeros(dom.total_dofs)
        F[end.dofs["ux"]] = self.A * FT
        return dom, end, F

    def run_bar(self, n, target, **tension):
        dom, end, F = self.bar(n, **tension)
        asm = Assembler(dom)
        path = [(0.0, 0.0)]
        s = IndirectDisplacementSolver(asm, control=[(end.id, "ux", 1.0)], max_lambda=1e9,
                                       initial_dl=target / 200, control_target=target,
                                       max_steps=5000)
        solidum.run(dom, assembler=asm, solver=s, F_applied=F,
                    step_callback=lambda k, U, lam: path.append((U[end.dofs["ux"]],
                                                                 lam * self.A * FT)))
        self.assertEqual(s.termination, "control_target")
        return np.array(path)

    def dissipated(self, path, n):
        d, f = path[:, 0], path[:, 1]
        W = float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(d)))
        sigma = f[-1] / self.A
        rest = self.LB * (1.0 - 1.0 / n)
        psi = 0.5 * sigma ** 2 / self.E * self.A * rest + 0.5 * f[-1] * (d[-1] - sigma / self.E * rest)
        return W - psi, sigma

    def test_con_energia_de_fractura(self):
        softening = lambda s: 2 * self.GF / self.FTW - s * (2 * self.GF / self.FTW ** 2 - self.LB / self.E)
        target = softening(0.05 * self.FTW)
        for n in (2, 4, 8, 16):
            with self.subTest(n=n):
                path = self.run_bar(n, target, G_f=self.GF)
                for d, F in path[1:]:
                    s = F / self.A
                    on_softening = d > self.FTW / self.E * self.LB * (1 + 1e-9)
                    ref = softening(s) if on_softening else s * self.LB / self.E
                    self.assertAlmostEqual(d, ref, delta=1e-9 * target)
                D, sigma = self.dissipated(path, n)
                # Tolerancia: la regla del trapecio en el paso que cruza el pico.
                self.assertAlmostEqual(D / (self.A * self.GF * (1 - sigma / self.FTW)), 1.0,
                                       delta=1e-4)

    def test_sin_regularizar_depende_de_la_malla(self):
        """Con ``ε_tu`` fija, la energía que disipa la grieta es el área bajo la
        curva por la longitud del elemento: con la grieta casi cerrada en
        tracción, proporcional a ``1/n``."""
        energies = {}
        for n in (2, 8):
            h = self.LB / n
            target = self.LB / self.E * 0.05 * self.FTW + 0.95 * EPS_TU * h
            path = self.run_bar(n, target, eps_tu=EPS_TU)
            energies[n], _ = self.dissipated(path, n)
        self.assertAlmostEqual(energies[2] / energies[8], 4.0, delta=0.4)


class TestCiclos(_Quiet):

    def test_ciclos_de_amplitud_creciente(self):
        N = -0.1 * FC0 * B * H
        dom, top, F_g, F_h = cantilever(8, make_section(ny=30), N_axial=N)
        targets = [0.02, -0.02, 0.04, -0.04, 0.06, -0.06]
        path, solvers = push(dom, top, F_g, F_h, targets, initial_dl=2e-3)
        self.assertTrue(all(s.termination == "control_target" for s in solvers))
        d, F = path[:, 0], path[:, 1]
        self.assertAlmostEqual(d[-1], -0.06, delta=1e-12)
        # Lazos con área: la fuerza al volver a δ = 0 no es nula (acero
        # plastificado), y el trabajo neto de los ciclos completos es positivo.
        W = np.cumsum(0.5 * (F[1:] + F[:-1]) * np.diff(d))
        self.assertGreater(W[-1], 0.0)
        # Fuerza en δ = 0 (interpolada) al volver de cada excursión: pequeña
        # tras el primer ciclo, cercano a la fluencia; del orden de la tercera
        # parte de la fuerza máxima tras los ciclos plásticos.
        residual = []
        for i in range(len(d) - 1):
            if d[i] != 0.0 and (d[i] > 0.0) != (d[i + 1] > 0.0) or d[i + 1] == 0.0 and d[i] != 0.0:
                residual.append(F[i] + (F[i + 1] - F[i]) * (0.0 - d[i]) / (d[i + 1] - d[i]))
        self.assertGreaterEqual(len(residual), 5)
        self.assertTrue(all(abs(r) > 0.15 * F.max() for r in residual[-3:]))
        self.assertAlmostEqual(-F.min() / F.max(), 1.0, delta=0.05)


class TestTransitorio(_Quiet):
    """Columna con una masa de 20 t en la cabeza: pulso lateral de medio seno
    y vibración libre (Newmark, aceleración media). El pulso agrieta el
    concreto y lleva el acero cerca de la fluencia. La energía disipada por el
    material es el trabajo de las fuerzas internas (``F_int = F_ext − M·ü`` en
    los grados de libertad libres) menos la energía almacenada en las fibras al
    final: no negativa y no decreciente."""

    def run_to(self, t_end, dt=5e-3, t_p=0.2, F0=4.0e4):
        section = make_section(ny=20, displaces=False, density=True)
        dom, top, _, _ = cantilever(6, section)
        dom.add_nodal_mass(top.id, ux=2.0e4, uy=2.0e4)
        dof = top.dofs["ux"]

        def F_func(t):
            F = np.zeros(dom.total_dofs)
            if t < t_p:
                F[dof] = F0 * math.sin(math.pi * t / t_p)
            return F

        asm = Assembler(dom)
        # Búsqueda lineal: el cierre de una grieta es un quiebre de la rigidez
        # (de la secante agrietada a la de compresión) y el Newton puro puede
        # rebotar entre las dos ramas (con dt = 5 ms, en t = 0.63 s).
        res = NewtonNewmarkSolver(asm, t_end=t_end, dt=dt, F_func=F_func,
                                  line_search=True).solve()
        self.assertEqual(res.termination, "t_end")
        M = asm.assemble_mass_matrix()
        free = [nd.dofs[k] for nd in dom.nodes.values() for k in nd.dofs
                if not nd.is_fixed(k)] if hasattr(next(iter(dom.nodes.values())), "is_fixed") \
            else [i for i in range(dom.total_dofs) if i not in
                  {dom.nodes[1].dofs[k] for k in ("ux", "uy", "rz")}]
        u, a, t = res.u_history, res.uddot_history, res.t_history
        Fint = np.array([F_func(ti)[free] for ti in t]).T - (M @ a)[free, :]
        W_int = float(np.sum(0.5 * (Fint[:, 1:] + Fint[:, :-1]) * np.diff(u[free, :], axis=1)))
        return W_int, self.stored_energy(dom, u[:, -1])

    @staticmethod
    def stored_energy(dom, U):
        """Energía libre de las fibras: ½·σ·(ε + e_p) en el concreto (spec §3–§4) y
        σ²/(2E) + ½·H·α² en el acero."""
        psi = 0.0
        for el in dom.elements.values():
            u_e = np.concatenate([U[[nd.dofs[k] for k in el.DOF_NAMES]] for nd in el.nodes])
            u_loc = el.T @ u_e
            sec = el.section
            for g in range(el.n_ip):
                e_g = el._B[g] @ u_loc
                eps = sec.strains_2d(e_g)
                for i, (st, A) in enumerate(zip(el.state.vars[g], sec.area)):
                    mat = sec._fiber_materials[i]
                    sigma = float(mat.compute_state(eps[i], st)[0])
                    if isinstance(mat, ManderConcrete1D):
                        em = st["eps_cmax"] if st else 0.0
                        ep = mat.plastic_strain(em, mat.compression_envelope(em)[0]) if em else 0.0
                        w = 0.5 * sigma * (eps[i] + ep)
                    else:
                        w = 0.5 * sigma ** 2 / mat.E + 0.5 * mat.H * (st["alpha"] if st else 0.0) ** 2
                    psi += w * A * el._wL[g]
        return psi

    def test_disipacion_del_material(self):
        D = []
        for t_end in (0.2, 0.5, 1.0):
            W_int, psi = self.run_to(t_end)
            D.append(W_int - psi)
        scale = max(abs(x) for x in D)
        self.assertGreater(D[0], -1e-4 * scale)
        self.assertGreater(D[-1], 0.0)
        for a, b in zip(D[:-1], D[1:]):
            self.assertGreater(b, a - 1e-3 * scale)


class TestObjetividadVoladizo(_Quiet):
    """Voladizo con ``N = −0.2·f'c·A_g`` empujado más allá del pico (el
    recubrimiento se aplasta). Concreto con ``eps_tu`` en tracción (rigidez
    entre grietas, que fija el refuerzo) y ``G_fc`` en compresión, calibrada
    para que con 8 elementos la rama no cambie (``λ = 1``).

    Criterio (decisión del usuario, 2026-11-09): con 8 y 16 elementos, la
    fuerza a desplazamiento dado después del pico difiere ≤ 5 % (medido:
    1.9 %). Con 4 elementos (0.75 m, más que la rótula plástica) la
    diferencia es 7.1 % (con 60 capas y hasta 0.15 m, 6.6 % frente a 16 % sin
    regularizar): el elemento en desplazamientos no concentra la curvatura en
    una zona más corta que él, y no se exige. El remedio es el elemento en
    fuerzas, pendiente (ROADMAP)."""

    N_AX = -0.2 * FC0 * B * H

    def section(self):
        h8 = (5.0 / 18.0) * L / 8                 # punto de Gauss extremo, n = 8
        g_core = ManderConcrete1D(confinement=CONF)._A_post * h8
        g_cover = ManderConcrete1D(**COVER)._A_post * h8
        core = ManderConcrete1D(confinement=CONF, ft=FT, eps_tu=EPS_TU, G_fc=g_core)
        cover = ManderConcrete1D(**COVER, ft=FT, eps_tu=EPS_TU, G_fc=g_cover)
        steel = Elastoplastic1D(E=ES, sigma_y=FY, H=HS)
        ny = 30
        return FiberSection(patches=[
            RectangularPatch(core, y=(-YB, YB), z=(-B / 2 + C, B / 2 - C), ny=ny, nz=1),
            RectangularPatch(cover, y=(YB, H / 2), z=(-B / 2, B / 2), ny=4, nz=1),
            RectangularPatch(cover, y=(-H / 2, -YB), z=(-B / 2, B / 2), ny=4, nz=1),
            RectangularPatch(cover, y=(-YB, YB), z=(-B / 2, -B / 2 + C), ny=ny, nz=1),
            RectangularPatch(cover, y=(-YB, YB), z=(B / 2 - C, B / 2), ny=ny, nz=1),
            BarLayer(steel, start=(YB, -0.1), end=(YB, 0.1), n=3, area=AB, displaces=core),
            BarLayer(steel, start=(-YB, -0.1), end=(-YB, 0.1), n=3, area=AB, displaces=core),
        ])

    def curve(self, n):
        dom, top, F_g, F_h = cantilever(n, self.section(), N_axial=self.N_AX)
        path, solvers = push(dom, top, F_g, F_h, [0.10], initial_dl=2e-3)
        self.assertEqual(solvers[-1].termination, "control_target")
        return path

    def test_ocho_y_dieciseis_elementos(self):
        p16, p8 = self.curve(16), self.curve(8)
        d = np.linspace(0.03, 0.10, 15)           # después del pico
        F16 = np.interp(d, p16[:, 0], p16[:, 1])
        F8 = np.interp(d, p8[:, 0], p8[:, 1])
        self.assertLess(np.max(np.abs(F8 / F16 - 1.0)), 5e-2)
        # Hay pico: la fuerza máxima se alcanza antes de 0.05 m y después baja.
        i_pk = int(np.argmax(p16[:, 1]))
        self.assertLess(p16[i_pk, 0], 0.05)
        self.assertLess(np.interp(0.06, p16[:, 0], p16[:, 1]), 0.97 * p16[i_pk, 1])


_YAML = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [0.0, 1.0]}
  - {id: 3, coords: [0.0, 2.0]}
  - {id: 4, coords: [0.0, 3.0]}
materials:
  - id: 1
    type: ManderConcrete1D
    confinement: {fc0: 30.0e6, f_yh: 420.0e6, eps_su: 0.09, k_e: 0.75, rho_x: 0.01, rho_y: 0.01}
    ft: 3.0e6
    G_f: 140.0
  - {id: 2, type: ManderConcrete1D, fc: 30.0e6, eps_c: 0.002, eps_cu: 0.004, eps_sp: 0.006, E: EC_, ft: 3.0e6, G_f: 140.0}
  - {id: 3, type: Elastoplastic1D, E: 2.0e11, sigma_y: 4.2e8, H: 2.0e9}
sections:
  - id: 7
    type: FiberSection
    patches:
      - {shape: rectangle, material: 1, y: [-0.21, 0.21], z: [-0.11, 0.11], ny: 30, nz: 1}
      - {shape: rectangle, material: 2, y: [0.21, 0.25], z: [-0.15, 0.15], ny: 4, nz: 1}
      - {shape: rectangle, material: 2, y: [-0.25, -0.21], z: [-0.15, 0.15], ny: 4, nz: 1}
      - {shape: rectangle, material: 2, y: [-0.21, 0.21], z: [-0.15, -0.11], ny: 30, nz: 1}
      - {shape: rectangle, material: 2, y: [-0.21, 0.21], z: [0.11, 0.15], ny: 30, nz: 1}
      - {shape: layer, material: 3, start: [0.21, -0.1], end: [0.21, 0.1], n: 3, area: AB_, displaces: 1}
      - {shape: layer, material: 3, start: [-0.21, -0.1], end: [-0.21, 0.1], n: 3, area: AB_, displaces: 1}
elements:
  - {id: 1, type: Frame2DEuler, nodes: [1, 2], section: 7}
  - {id: 2, type: Frame2DEuler, nodes: [2, 3], section: 7}
  - {id: 3, type: Frame2DEuler, nodes: [3, 4], section: 7}
boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
point_loads:
  - {node_id: 4, ux: 1.0e5}
solver:
  type: IndirectDisplacementSolver
  control:
    - {node: 4, dof: ux, coef: 1.0}
  initial_dl: 0.002
  control_target: 0.02
  max_lambda: 1.0e9
  max_steps: 500
"""


class TestYaml(_Quiet):

    def test_igual_que_python(self):
        import tempfile, os
        text = _YAML.replace("EC_", repr(EC)).replace("AB_", repr(AB))
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "concreto.yaml")
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text)
            captured = {}
            original = YamlParser.parse

            def parse(self):
                captured["dom"] = original(self)
                return captured["dom"]

            YamlParser.parse = parse
            try:
                res_y = run_yaml(path)
            finally:
                YamlParser.parse = original
        dom_y = captured["dom"]

        core = ManderConcrete1D(confinement=CONF, ft=FT, G_f=140.0)
        cover = ManderConcrete1D(**COVER, ft=FT, G_f=140.0)
        steel = Elastoplastic1D(E=ES, sigma_y=FY, H=HS)
        section = FiberSection(patches=[
            RectangularPatch(core, y=(-0.21, 0.21), z=(-0.11, 0.11), ny=30, nz=1),
            RectangularPatch(cover, y=(0.21, 0.25), z=(-0.15, 0.15), ny=4, nz=1),
            RectangularPatch(cover, y=(-0.25, -0.21), z=(-0.15, 0.15), ny=4, nz=1),
            RectangularPatch(cover, y=(-0.21, 0.21), z=(-0.15, -0.11), ny=30, nz=1),
            RectangularPatch(cover, y=(-0.21, 0.21), z=(0.11, 0.15), ny=30, nz=1),
            BarLayer(steel, start=(0.21, -0.1), end=(0.21, 0.1), n=3, area=AB, displaces=core),
            BarLayer(steel, start=(-0.21, -0.1), end=(-0.21, 0.1), n=3, area=AB, displaces=core),
        ])
        dom = Domain()
        nodes = [dom.add_node(i + 1, [0.0, float(i)]) for i in range(4)]
        for k in range(3):
            dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], section=section))
        for dof in ("ux", "uy", "rz"):
            nodes[0].fix_dof(dof, 0.0)
        dom.generate_equation_numbers()
        F = np.zeros(dom.total_dofs)
        F[nodes[3].dofs["ux"]] = 1.0e5
        asm = Assembler(dom)
        solver = IndirectDisplacementSolver(asm, control=[(4, "ux", 1.0)], initial_dl=0.002,
                                            control_target=0.02, max_lambda=1e9, max_steps=500)
        res_p = solidum.run(dom, assembler=asm, solver=solver, F_applied=F)
        for nid in (2, 3, 4):
            for k in ("ux", "uy", "rz"):
                self.assertAlmostEqual(res_y.U[dom_y.nodes[nid].dofs[k]],
                                       res_p.U[dom.nodes[nid].dofs[k]],
                                       delta=1e-12 * np.abs(res_p.U).max())
        # El concreto se agrietó (respuesta no lineal) y la regularización actuó.
        el = dom.elements[1]
        self.assertTrue(np.allclose(el.characteristic_lengths, el._gauss_w * 1.0))
        self.assertTrue(any(st and st["eps_tmax"] > EPS_T
                            for st in el.state.vars[0] if st and "eps_tmax" in st))


if __name__ == "__main__":
    unittest.main()
