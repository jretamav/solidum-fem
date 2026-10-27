"""Capa 4 de la red de seguridad: aviso de pérdida de cifras por mal
condicionamiento (ADR 0019, decisión del usuario del 2026-10-27).

Un sistema resuelto correctamente (capa 2, error hacia atrás) y no singular
(capa 3, mecanismos) puede haber perdido cifras si está mal condicionado: un
solver directo es estable hacia atrás y el residuo no lo delata. El
``LinearSolver`` estima el número de condición de la matriz escalada por la
diagonal, ``κ̃ = ‖K̃‖₁·‖K̃⁻¹‖₁``, con la factorización ya hecha (Higham y
Tisseur, ``onenormest``), lo guarda en ``SolveResult.condition_estimate`` y
avisa si supera ``CONDITION_WARNING_THRESHOLD = 1e12``.

Comprobaciones:

- el estimador coincide con ``κ₁`` exacto de la matriz densa en modelos
  pequeños y no depende de las unidades;
- avisa en un voladizo de Euler-Bernoulli de 1000 elementos (``κ̃ ≈ 1e13``,
  error real de la flecha 3e-5) y en un pórtico con zonas rígidas ×1e10
  (``κ̃ ≈ 1.3e14``); calla con zonas ×1e6 (``1.3e10``), ν = 0.4999999
  (``4.2e9``) y una placa normal;
- el caso que motivó la capa: con 10 000 elementos el voladizo da una flecha
  errónea en un 51 % (``κ̃ ≈ 1.5e17``), que antes salía sin ningún aviso;
- ``κ̃`` llega al resultado del estático lineal y del térmico estacionario, y
  es ``None`` con el backend iterativo y en el no lineal.

Medido sobre las 796 factorizaciones lineales de la suite (2026-10-27): la
mediana de ``κ̃`` es 1e3, el 99 % queda bajo 1.3e10 y el modelo legítimo más
mal condicionado (viga de MacNeal con Hex20) da 1.2e10; sólo el pórtico con
zonas rígidas ×1e8, extremo a propósito, pasa de 1e12.
"""
import logging
import os
import sys
import unittest

import numpy as np
import scipy.sparse as sp

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(__file__))

import solidum  # noqa: E402
from test_internal_mechanisms import (  # noqa: E402
    placa_casi_incompresible, portico_con_zonas_rigidas)

from solidum.constants import CONDITION_WARNING_THRESHOLD  # noqa: E402
from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.frame import Frame2DEuler  # noqa: E402
from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.materials.elastic_2d import Elastic2D  # noqa: E402
from solidum.materials.elastic_3d import Elastic3D  # noqa: E402
from solidum.materials.thermal_conduction import ThermalConduction  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import LinearSolver, NonlinearSolver  # noqa: E402
from solidum.math.solvers.model_checks import backward_error_scaling  # noqa: E402
from solidum.mesh import structured_mesh  # noqa: E402

E, A, I, L, P = 210e9, 0.01, 8.33e-6, 10.0, 1e3


def voladizo(n):
    """Voladizo Frame2DEuler de n elementos con una carga P en la punta: los
    nodos son exactos, así que el error de la flecha es la pérdida real."""
    dom = Domain()
    mat = Elastic1D(E)
    for i in range(n + 1):
        dom.add_node(i + 1, [L * i / n, 0.0])
    for i in range(n):
        dom.add_element(Frame2DEuler(i + 1, [dom.nodes[i + 1], dom.nodes[i + 2]], mat, A=A, I=I))
    for d in ("ux", "uy", "rz"):
        dom.nodes[1].fix_dof(d, 0.0)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    F[dom.nodes[n + 1].dofs["uy"]] = P
    return dom, F


def bloque(nombre, n, material):
    m = structured_mesh(nombre, (n,) * (2 if nombre in ("Quad4", "Quad4Thermal") else 3),
                        material=material)
    for nd in m.nodes_on("xmin"):
        for d in nd.dofs:
            nd.fix_dof(d, 0.0)
    m.domain.generate_equation_numbers()
    return m.domain


def portico(c, unidades="m"):
    dom, _ = portico_con_zonas_rigidas(c, unidades)
    dom.generate_equation_numbers()
    f = 1e3 if unidades == "mm" else 1.0
    nudo = next(n for n in dom.nodes.values() if np.allclose(n.coordinates, [0.0, 4.0 * f]))
    F = np.zeros(dom.total_dofs)
    F[nudo.dofs["ux"]] = 5e4
    return dom, F


def resolver(dom, F=None, **kw):
    """``(resultado, solver)`` con ``solidum.run`` y un ``LinearSolver``."""
    solver = LinearSolver(Assembler(dom), **kw)
    res = solidum.run(dom, assembler=solver.assembler, solver=solver,
                      F_applied=np.zeros(dom.total_dofs) if F is None else F)
    return res, solver


class TestEstimador(unittest.TestCase):

    def setUp(self):
        logging.disable(logging.WARNING)

    def tearDown(self):
        logging.disable(logging.NOTSET)

    def test_coincide_con_kappa_exacto(self):
        """``onenormest`` da una cota inferior; en estos modelos, exacta."""
        modelos = {
            "voladizo de 50 elementos": voladizo(50)[0],
            "placa Quad4 8x8": bloque("Quad4", 8, Elastic2D(200e3, 0.3)),
            "bloque Hex8 4x4x4": bloque("Hex8", 4, Elastic3D(200e3, 0.3)),
        }
        for nombre, dom in modelos.items():
            with self.subTest(modelo=nombre):
                _, solver = resolver(dom)
                s, _ = backward_error_scaling(solver._K_red)
                Kt = (sp.diags(s) @ solver._K_red @ sp.diags(s)).toarray()
                exacto = np.linalg.norm(Kt, 1) * np.linalg.norm(np.linalg.inv(Kt), 1)
                self.assertGreaterEqual(solver.condition_estimate, 0.5 * exacto)
                self.assertLessEqual(solver.condition_estimate, exacto * (1 + 1e-8))

    def test_no_depende_de_las_unidades(self):
        k = {u: resolver(*portico(1e6, u))[1].condition_estimate for u in ("m", "mm")}
        self.assertAlmostEqual(k["m"] / k["mm"], 1.0, places=6)

    def test_llega_al_resultado(self):
        res, solver = resolver(*voladizo(10))
        self.assertIsNotNone(res.condition_estimate)
        self.assertEqual(res.condition_estimate, solver.condition_estimate)

    def test_conduccion_estacionaria(self):
        dom = bloque("Quad4Thermal", 8, ThermalConduction(45.0))
        res, _ = resolver(dom)
        self.assertGreater(res.condition_estimate, 1.0)

    def test_sin_estimacion_con_el_iterativo(self):
        res, solver = resolver(*voladizo(10), linear_algebra="iterative")
        self.assertIsNone(res.condition_estimate)
        self.assertIsNone(solver.condition_estimate)

    def test_sin_estimacion_en_el_no_lineal(self):
        dom, F = voladizo(10)
        res = solidum.run(dom, solver=NonlinearSolver(Assembler(dom)), F_applied=F)
        self.assertIsNone(res.condition_estimate)


class TestAviso(unittest.TestCase):

    def avisa(self, dom, F):
        with self.assertLogs("solidum", level="WARNING") as cm:
            res, _ = resolver(dom, F)
        mensajes = [m for m in cm.output if "mal condicionado" in m]
        self.assertTrue(mensajes, cm.output)
        self.assertIn("cifras", mensajes[0])
        return res

    def calla(self, dom, F):
        logger = logging.getLogger("solidum")
        registros = []

        class Recolector(logging.Handler):
            def emit(self, record):
                registros.append(record.getMessage())

        h = Recolector(level=logging.WARNING)
        logger.addHandler(h)
        try:
            res, _ = resolver(dom, F)
        finally:
            logger.removeHandler(h)
        self.assertFalse([m for m in registros if "mal condicionado" in m], registros)
        return res

    def test_voladizo_de_1000_elementos(self):
        res = self.avisa(*voladizo(1000))
        self.assertGreater(res.condition_estimate, CONDITION_WARNING_THRESHOLD)

    def test_zonas_rigidas_extremas(self):
        self.avisa(*portico(1e10))

    def test_calla_en_modelos_normales_y_contrastes_moderados(self):
        casos = {
            "zonas rígidas ×1e6": portico(1e6),
            "ν = 0.4999999": placa_casi_incompresible(0.4999999),
            "voladizo de 100 elementos": voladizo(100),
        }
        for nombre, (dom, F) in casos.items():
            with self.subTest(caso=nombre):
                dom.generate_equation_numbers()
                res = self.calla(dom, F)
                self.assertLess(res.condition_estimate, CONDITION_WARNING_THRESHOLD)

    def test_el_caso_que_motivo_la_capa(self):
        """Voladizo de 10 000 elementos: la flecha sale errónea (51 % y
        κ̃ = 1.5e17 medidos el 2026-10-27; las cotas de la prueba dejan margen
        al redondeo de cada plataforma), y ahora lo avisa en vez de callar."""
        dom, F = voladizo(10_000)
        res = self.avisa(dom, F)
        flecha = res.U[dom.nodes[10_001].dofs["uy"]]
        self.assertGreater(abs(flecha / (P * L**3 / (3 * E * I)) - 1.0), 1e-3)
        self.assertGreater(res.condition_estimate, 1e15)


if __name__ == "__main__":
    unittest.main()
