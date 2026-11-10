"""Transitorios implícitos con grados de libertad sin masa (deuda #34;
decisión del usuario del 2026-11-09, ADR 0024 §6).

La aceleración inicial se toma consistente: ``üₘ = Mₘₘ⁻¹·rₘ`` en los grados
de libertad con masa y ``ü₀ = −K₀₀⁻¹·K₀ₘ·üₘ`` en los demás. Criterios:

1. Edificio de cortante de tres pisos (masas de piso, columnas sin masa: 18
   de 21 grados de libertad sin masa) con una carga lateral aplicada desde
   ``t = 0``: el desplazamiento de cada piso es el del problema 3×3
   condensado integrado con el mismo esquema (Newmark, y HHT-α escrito aquí),
   sin y con amortiguamiento de Rayleigh, a 1e-9.
2. La aceleración de los grados de libertad sin masa es la consistente,
   ``ü₀(t) = −K₀₀⁻¹·K₀ₘ·üₘ(t)``, en todo instante (a 1e-9 de la de piso).
3. Las variantes de Newton dan lo mismo que las lineales en un modelo
   elástico.
4. Una carga aplicada de golpe en un nudo sin masa se avisa.
"""
import logging
import math
import unittest

import numpy as np

from solidum import Assembler
from solidum.math.solvers import HHTSolver, NewmarkSolver, NewtonHHTSolver, NewtonNewmarkSolver
from test_nodal_mass import shear_building, shear_building_3x3
from test_transient_continuation import tight

HEIGHTS, MASSES, INERTIAS = [4.0, 3.2, 3.0], [5.2e4, 4.6e4, 3.1e4], [4.2e-4, 3.0e-4, 2.2e-4]
P = 1.0e5
RTOL = 1e-9


def _oracle(alpha, beta, gamma, C, dt, n_steps):
    """HHT-α (Newmark con α = 0) sobre el problema 3×3 condensado, carga P en
    la azotea desde t = 0."""
    K, M = shear_building_3x3(HEIGHTS, MASSES, INERTIAS)
    F = np.array([0.0, 0.0, P])
    u, v = np.zeros(3), np.zeros(3)
    a = np.linalg.solve(M, F - C @ v - K @ u)
    A = M + (1 + alpha) * gamma * dt * C + (1 + alpha) * beta * dt * dt * K
    hist_u, hist_a = [u.copy()], [a.copy()]
    for _ in range(n_steps):
        up = u + dt * v + (0.5 - beta) * dt * dt * a
        vp = v + (1 - gamma) * dt * a
        rhs = (F - (1 + alpha) * (C @ vp) + alpha * (C @ v)
               - (1 + alpha) * (K @ up) + alpha * (K @ u))
        a = np.linalg.solve(A, rhs)
        u, v = up + beta * dt * dt * a, vp + gamma * dt * a
        hist_u.append(u.copy())
        hist_a.append(a.copy())
    return np.array(hist_u).T, np.array(hist_a).T


def _rayleigh(xi):
    """Coeficientes de Rayleigh con ξ en el primer y el tercer modo."""
    if xi == 0.0:
        return 0.0, 0.0
    K, M = shear_building_3x3(HEIGHTS, MASSES, INERTIAS)
    w = np.sqrt(np.sort(np.linalg.eigvals(np.linalg.solve(M, K)).real))
    return 2 * xi * w[0] * w[2] / (w[0] + w[2]), 2 * xi / (w[0] + w[2])


class _Quiet(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)


class TestEdificioDeCortante(_Quiet):

    DT, T_END = 0.01, 1.0

    def _run(self, make, xi=0.0):
        dom, floors = shear_building(HEIGHTS, MASSES, INERTIAS)
        dom.generate_equation_numbers()
        roof = floors[-1].dofs["ux"]

        def F_func(t):
            F = np.zeros(dom.total_dofs)
            F[roof] = P
            return F
        a0, a1 = _rayleigh(xi)
        asm = Assembler(dom)
        res = make(asm, F_func, {"alpha": a0, "beta": a1}).solve()
        return dom, asm, floors, res

    def _compare(self, res, floors, alpha, beta, gamma, xi):
        K3, M3 = shear_building_3x3(HEIGHTS, MASSES, INERTIAS)
        a0, a1 = _rayleigh(xi)
        n = res.n_steps
        u_ref, a_ref = _oracle(alpha, beta, gamma, a0 * M3 + a1 * K3, self.DT, n)
        for k, floor in enumerate(floors):
            u = res.u_history[floor.dofs["ux"]]
            np.testing.assert_allclose(u, u_ref[k], rtol=0.0, atol=RTOL * np.max(np.abs(u_ref)))
            a = res.uddot_history[floor.dofs["ux"]]
            np.testing.assert_allclose(a, a_ref[k], rtol=0.0, atol=RTOL * np.max(np.abs(a_ref)))

    def test_newmark_y_hht_contra_el_problema_condensado(self):
        alpha_hht = -0.1
        beta_hht, gamma_hht = (1 - alpha_hht) ** 2 / 4, (1 - 2 * alpha_hht) / 2
        cases = {
            "Newmark": (lambda asm, F, r: NewmarkSolver(asm, t_end=self.T_END, dt=self.DT,
                                                       F_func=F, rayleigh=r),
                        (0.0, 0.25, 0.5)),
            "HHT": (lambda asm, F, r: HHTSolver(asm, t_end=self.T_END, dt=self.DT, F_func=F,
                                               alpha=alpha_hht, rayleigh=r),
                    (alpha_hht, beta_hht, gamma_hht)),
        }
        for name, (make, coefs) in cases.items():
            for xi in (0.0, 0.05):
                with self.subTest(solver=name, xi=xi):
                    dom, asm, floors, res = self._run(make, xi)
                    self.assertEqual(res.termination, "t_end")
                    self._compare(res, floors, *coefs, xi)

    def test_aceleracion_consistente_en_los_grados_sin_masa(self):
        dom, asm, floors, res = self._run(
            lambda asm, F, r: NewmarkSolver(asm, t_end=self.T_END, dt=self.DT, F_func=F))
        asm.assemble_system()
        K_red, M_red, T = asm.reduce_pair(asm.K_global, asm.assemble_mass_matrix())
        K_red, M_red = K_red.toarray(), M_red.toarray()
        free = asm.constraint_set.free_dofs(dom.total_dofs)
        zero = np.flatnonzero(np.diag(M_red) == 0.0)
        mass = np.setdiff1d(np.arange(M_red.shape[0]), zero)
        self.assertEqual((zero.size, mass.size), (18, 3))
        G = np.linalg.solve(K_red[np.ix_(zero, zero)], K_red[np.ix_(zero, mass)])
        acc = res.uddot_history[free]            # libres = maestros reducidos
        consistent = -G @ acc[mass]
        scale = np.max(np.abs(acc[mass]))
        np.testing.assert_allclose(acc[zero], consistent, rtol=0.0, atol=RTOL * scale)
        # En t = 0 la carga ya actúa: la azotea acelera P/m₃ y los grados de
        # libertad sin masa conectados a ella la siguen (no arrancan en cero).
        roof = floors[-1].dofs["ux"]
        self.assertAlmostEqual(res.uddot_history[roof, 0] / (P / MASSES[-1]), 1.0, delta=1e-12)
        self.assertGreater(np.max(np.abs(acc[zero, 0])), 1e-3 * P / MASSES[-1])

    def test_variantes_de_newton(self):
        pairs = [
            (lambda asm, F, r: NewmarkSolver(asm, t_end=0.3, dt=self.DT, F_func=F),
             lambda asm, F, r: NewtonNewmarkSolver(asm, t_end=0.3, dt=self.DT, F_func=F,
                                                   convergence=tight())),
            (lambda asm, F, r: HHTSolver(asm, t_end=0.3, dt=self.DT, F_func=F),
             lambda asm, F, r: NewtonHHTSolver(asm, t_end=0.3, dt=self.DT, F_func=F,
                                               convergence=tight())),
        ]
        for linear, newton in pairs:
            _, _, floors, r_lin = self._run(linear)
            _, _, _, r_new = self._run(newton)
            scale = np.max(np.abs(r_lin.u_history))
            np.testing.assert_allclose(r_new.u_history, r_lin.u_history, rtol=0.0,
                                       atol=1e-8 * scale)
            np.testing.assert_allclose(r_new.uddot_history, r_lin.uddot_history, rtol=0.0,
                                       atol=1e-7 * np.max(np.abs(r_lin.uddot_history)))


class TestCargaEnNudoSinMasa(unittest.TestCase):

    def test_aviso(self):
        dom, floors = shear_building(HEIGHTS, MASSES, INERTIAS)
        dom.generate_equation_numbers()
        # Un nudo intermedio de columna (sin masa) del primer entrepiso.
        mid = next(n for n in dom.nodes.values()
                   if 0.0 < n.coordinates[1] < HEIGHTS[0] and n.coordinates[0] == 0.0)
        F = np.zeros(dom.total_dofs)
        F[mid.dofs["ux"]] = P
        solver = NewmarkSolver(Assembler(dom), t_end=0.05, dt=0.01, F_func=lambda t: F)
        with self.assertLogs("solidum", level="WARNING") as logs:
            res = solver.solve()
        self.assertTrue(any("no está equilibrada" in m and f"nodo {mid.id} (ux)" in m
                            for m in logs.output))
        self.assertEqual(res.termination, "t_end")
        self.assertTrue(math.isfinite(float(np.max(np.abs(res.u_history)))))


if __name__ == "__main__":
    unittest.main()
