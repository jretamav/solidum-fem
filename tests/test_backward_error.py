"""Capa 2 de la red de seguridad: error hacia atrás normado de la solución
(ADR 0019, decisión del usuario del 2026-10-20).

Hasta el 2026-10-20 la capa 2 exigía ``‖F − K·u‖ ≤ 1e-8·‖F‖``. Un solver
directo estable garantiza un residuo del orden de ``ε·‖K‖·‖u‖``, no de
``ε·‖F‖``: con rigideces muy desproporcionadas el residuo relativo rechazaba
soluciones buenas (zonas rígidas ×1e6 en mm, ν = 0.4999999), y la propia suite
llegaba a 5.7e-9, al 57 % del umbral. Ahora se exige el error hacia atrás
normado, escalado por la diagonal,

    η = ‖r̃‖∞ / (‖K̃‖∞·‖ũ‖∞ + ‖F̃‖∞),  K̃ = D^-½·K·D^-½,  D = diag(K),

``≤ 1e-10``: mide si el solver resolvió el sistema que se le dio, sin depender
de las unidades ni del contraste. Los mecanismos internos cargados, que
delataba el residuo, los rechaza antes la capa 3 (``test_internal_mechanisms``).
"""
import logging
import os
import sys
import unittest

import numpy as np
import scipy.sparse as sp

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(__file__))

from test_internal_mechanisms import (  # noqa: E402
    DIRECT, placa_casi_incompresible, portico_con_zonas_rigidas, reducida)

from solidum.constants import BACKWARD_ERROR_RTOL  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.linalg import LUSolver  # noqa: E402
from solidum.math.solvers import LinearSolver  # noqa: E402
from solidum.math.solvers.diagnostics import IllPosedSystemError  # noqa: E402
from solidum.math.solvers.model_checks import backward_error, check_linear_solution  # noqa: E402

P = 5.0e4   # N, carga lateral en el nudo superior izquierdo


def setUpModule():
    logging.disable(logging.WARNING)


def tearDownModule():
    logging.disable(logging.NOTSET)


def portico_cargado(contraste, unidades):
    dom, _ = portico_con_zonas_rigidas(contraste, unidades)
    f = 1e3 if unidades == "mm" else 1.0
    nudo = next(n for n in dom.nodes.values()
                if np.allclose(n.coordinates, [0.0, 4.0 * f]))
    F = np.zeros(dom.total_dofs)
    F[nudo.dofs["ux"]] = P
    return dom, F, nudo


class TestContrastesExtremosSeResuelven(unittest.TestCase):

    def test_zonas_rigidas_en_m_y_en_mm(self):
        """Se resuelven, con el mismo desplazamiento en las dos unidades. El
        criterio anterior (residuo relativo ≤ 1e-8) los rechazaba: se comprueba
        aquí mismo, para que la prueba demuestre algo.

        La coincidencia entre unidades la limita el condicionamiento del
        modelo, no el solver: cada solución tiene un error relativo de hasta
        ``κ̃·ε`` (``κ̃`` el de la matriz escalada por su diagonal, el mismo en
        las dos unidades), así que la cota de la diferencia es ``2·κ̃·ε``:
        2.5e-6·2 con c = 1e6 y 2.5e-4·2 con c = 1e8. Medido: 3.3e-7 con
        Pardiso en Linux y 2e-5 en Windows; una tolerancia fija ajustada a una
        sola plataforma falló en el CI (2026-10-20). Avisar de esa pérdida de
        cifras es la carencia «mal condicionamiento sin singularidad» del
        ADR 0019."""
        for c in (1e6, 1e8):
            for la in DIRECT:
                with self.subTest(contraste=c, backend=la):
                    ux = {}
                    _, K = reducida(portico_cargado(c, "m")[0])
                    s = 1.0 / np.sqrt(np.abs(K.diagonal()))
                    kappa = np.linalg.cond(s[:, None] * K.toarray() * s[None, :])
                    rtol = 2.0 * kappa * np.finfo(float).eps
                    for unidades in ("m", "mm"):
                        dom, F, nudo = portico_cargado(c, unidades)
                        a = Assembler(dom)
                        U = LinearSolver(a, linear_algebra=la).solve(F)
                        ux[unidades] = U[nudo.dofs["ux"]]
                        if unidades == "mm":
                            _, K = reducida(dom)
                            Fr = F[a.constraint_set.free_dofs(dom.total_dofs)]
                            u = LUSolver().factorize(K).solve(Fr)
                            self.assertGreater(np.linalg.norm(Fr - K @ u) / np.linalg.norm(Fr), 1e-8)
                    self.assertLess(abs(ux["mm"] / (1e3 * ux["m"]) - 1.0), rtol)

    def test_material_casi_incompresible(self):
        for nu in (0.4999999, 0.49999999):
            for la in DIRECT:
                with self.subTest(nu=nu, backend=la):
                    dom, F = placa_casi_incompresible(nu)
                    U = LinearSolver(Assembler(dom), linear_algebra=la).solve(F)
                    self.assertTrue(np.all(np.isfinite(U)))


class TestErrorHaciaAtras(unittest.TestCase):

    def setUp(self):
        rng = np.random.default_rng(7)
        A = rng.standard_normal((30, 30))
        self.K = sp.csr_matrix(A @ A.T + 30 * np.eye(30))
        self.F = rng.standard_normal(30)
        self.u = np.linalg.solve(self.K.toarray(), self.F)

    def test_no_depende_de_las_unidades(self):
        """``S·K·S``, ``S⁻¹·u`` y ``S·F`` describen el mismo problema en otras
        unidades: el mismo ``η``. Con una solución perturbada, para que ``η``
        no sea ruido de redondeo."""
        S = np.diag(10.0 ** np.random.default_rng(3).uniform(-6, 6, 30))
        Ks = sp.csr_matrix(S @ self.K.toarray() @ S)
        u = self.u * (1 + 1e-6 * np.random.default_rng(5).standard_normal(30))
        eta = backward_error(self.K, u, self.F)
        eta_s = backward_error(Ks, np.linalg.solve(S, u), S @ self.F)
        self.assertGreater(eta, 1e-9)
        self.assertAlmostEqual(eta_s / eta, 1.0, places=6)

    def test_rechaza_una_solucion_que_no_resuelve_el_sistema(self):
        self.assertLess(check_linear_solution(self.K, self.u, self.F), BACKWARD_ERROR_RTOL)
        malo = self.u * (1 + 1e-6 * np.random.default_rng(5).standard_normal(30))
        with self.assertRaises(IllPosedSystemError) as cm:
            check_linear_solution(self.K, malo, self.F)
        self.assertGreater(cm.exception.backward_error, BACKWARD_ERROR_RTOL)
        self.assertIn("error hacia atrás", str(cm.exception))

    def test_sin_carga(self):
        self.assertEqual(backward_error(self.K, np.zeros(30), np.zeros(30)), 0.0)


if __name__ == "__main__":
    unittest.main()
