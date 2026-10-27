"""Cilindro hueco en régimen estacionario: perfil logarítmico (Carslaw y Jaeger §7.2).

Benchmark ``cilindro_hueco_logaritmico`` de las specs de ``Quad4Thermal`` y
``Hex8Thermal``, pendiente desde la Etapa 8 por falta de un mallador de
coronas; lo da ``solidum.mesh.structured_mesh`` (ADR 0022).

Corona ``a ≤ r ≤ b``, cuarto de circunferencia, con ``T(a) = T_a`` y
``T(b) = T_b`` impuestas, sin fuente y con los lados radiales adiabáticos (en
3D también las caras ``z = 0`` y ``z = 1``): por simetría el flujo es radial y

    T(r) = T_a + (T_b − T_a)·ln(r/a)/ln(b/a),

independiente de ``k``. Se resuelve con el camino de producción
(``Assembler`` y ``LinearSolver``), sin las soluciones manufacturadas.

Tolerancia de las specs: error nodal ≤ 5·10⁻³ relativo a ``|T_b − T_a|``.
Medido el 2026-10-27: 4.4·10⁻⁴, 1.2·10⁻⁴ y 3.0·10⁻⁵ con n = 4, 8 y 16
(orden 2 en los nodos), igual en 2D y en 3D.
"""
import logging
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import solidum  # noqa: F401,E402
from solidum.materials.thermal_conduction import ThermalConduction  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import LinearSolver  # noqa: E402
from solidum.mesh import structured_mesh  # noqa: E402

A_R, B_R, T_A, T_B, K = 1.0, 2.0, 100.0, 20.0, 45.0


def corona(X):
    r, t = A_R + (B_R - A_R) * X[0], X[1] * np.pi / 2.0
    return np.array([r * np.cos(t), r * np.sin(t), *X[2:]])


def analitica(x):
    r = np.hypot(x[:, 0], x[:, 1])
    return T_A + (T_B - T_A) * np.log(r / A_R) / np.log(B_R / A_R)


def error_nodal(nombre, d, n):
    """Máximo error nodal relativo a ``|T_b − T_a|``."""
    m = structured_mesh(nombre, (n,) * d, mapping=corona, material=ThermalConduction(K, dim=d))
    for nd in m.nodes_on("xmin"):
        nd.fix_dof("T", T_A)
    for nd in m.nodes_on("xmax"):
        nd.fix_dof("T", T_B)
    m.domain.generate_equation_numbers(verbose=False)
    T = LinearSolver(Assembler(m.domain)).solve(np.zeros(m.domain.total_dofs))
    nodos = list(m.domain.nodes.values())
    x = np.array([nd.coordinates for nd in nodos])
    return np.max(np.abs(T[[nd.dofs["T"] for nd in nodos]] - analitica(x))) / abs(T_B - T_A)


def setUpModule():
    logging.disable(logging.WARNING)


def tearDownModule():
    logging.disable(logging.NOTSET)


class TestCilindroHueco(unittest.TestCase):

    def test_perfil_logaritmico(self):
        for nombre, d in (("Quad4Thermal", 2), ("Hex8Thermal", 3)):
            with self.subTest(elemento=nombre):
                e = [error_nodal(nombre, d, n) for n in (4, 8, 16)]
                self.assertLess(e[0], 5e-3)
                self.assertGreater(np.log2(e[1] / e[2]), 1.9)

    def test_3d_igual_que_2d(self):
        """Con las caras ``z`` adiabáticas el campo 3D es el plano."""
        self.assertAlmostEqual(error_nodal("Hex8Thermal", 3, 4),
                               error_nodal("Quad4Thermal", 2, 4), delta=1e-12)


if __name__ == "__main__":
    unittest.main()
