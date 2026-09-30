"""Red de seguridad del ADR 0019 en ``DissipationArcLengthSolver`` (2026-09-30).

El solver del módulo ``dissipation_control`` llama a las mismas piezas
genéricas y públicas que los solvers estáticos del programa principal
(``solidum.math.solvers.model_checks``): ``ensure_statically_restrained``
antes de resolver y ``check_initial_tangent``, que inspecciona la primera
factorización del corrector, la de la tangente en ``U = 0`` (el predictor del
primer paso).

Sobre la armadura del ejemplo 11 del manual de ejemplos
(``tests/_armadura_pratt.py``), sin la diagonal del segundo panel y con una
carga horizontal que no activa el mecanismo interno:

- barras ``Truss2D`` (geométricamente lineal) ⇒ ``IllPosedSystemError``;
- barras ``Truss2DCorot`` ⇒ sólo un aviso, y el análisis sigue;
- con la diagonal, nada;
- con el backend iterativo, el aviso de que la comprobación no se hace.
"""
import contextlib
import logging
import os
import sys
import unittest

import numpy as np

_TESTS = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.abspath(os.path.join(_TESTS, '..')))
sys.path.insert(0, _TESTS)

import solidum  # noqa: E402
from solidum.elements.truss import Truss2D, Truss2DCorot  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.linalg.dispatcher import available_overrides  # noqa: E402
from solidum.math.solvers.diagnostics import IllPosedSystemError  # noqa: E402

from _armadura_pratt import armadura  # noqa: E402

solidum.load_user_module("dissipation_control")
from solidum.user.dissipation_control import DissipationArcLengthSolver  # noqa: E402

DIRECT = [b for b in ("lu", "pardiso") if b in available_overrides()]
AVISO_TANGENTE = "matriz tangente inicial (U = 0)"
AVISO_ITERATIVO = "no se comprueba si el modelo tiene un mecanismo interno"


@contextlib.contextmanager
def _avisos():
    """Mensajes WARNING del logger de los solvers durante el bloque."""
    registros = []

    class _Lista(logging.Handler):
        def emit(self, record):
            registros.append(record.getMessage())

    logger = logging.getLogger("solidum.solvers")
    handler = _Lista(level=logging.WARNING)
    previo = logging.root.manager.disable
    logging.disable(logging.NOTSET)
    logger.addHandler(handler)
    try:
        yield registros
    finally:
        logger.removeHandler(handler)
        logging.disable(previo)


def _solver(asm, la):
    return DissipationArcLengthSolver(asm, max_lambda=1.0, initial_tau=1.0,
                                      linear_algebra=la)


class TestDctlRedDeSeguridad(unittest.TestCase):

    def test_mecanismo_interno_geometricamente_lineal_se_rechaza(self):
        for la in DIRECT:
            with self.subTest(backend=la):
                dom, F = armadura(sin_diagonal=True, carga="horizontal", elemento=Truss2D)
                with self.assertRaises(IllPosedSystemError) as cm:
                    _solver(Assembler(dom), la).solve(F)
                self.assertTrue(cm.exception.initial_tangent)
                self.assertIn("DissipationArcLengthSolver", str(cm.exception))

    def test_corrotacional_solo_avisa(self):
        for la in DIRECT:
            with self.subTest(backend=la):
                dom, F = armadura(sin_diagonal=True, carga="horizontal",
                                  elemento=Truss2DCorot)
                with _avisos() as registros:
                    _solver(Assembler(dom), la).solve(F)
                aviso = [r for r in registros if AVISO_TANGENTE in r]
                self.assertEqual(len(aviso), 1, registros)
                self.assertIn("DissipationArcLengthSolver", aviso[0])
                self.assertIn("Truss2DCorot", aviso[0])

    def test_bien_planteado_no_avisa(self):
        for elemento in (Truss2D, Truss2DCorot):
            for la in DIRECT:
                with self.subTest(elemento=elemento.__name__, backend=la):
                    dom, F = armadura(sin_diagonal=False, carga="horizontal",
                                      elemento=elemento)
                    with _avisos() as registros:
                        solver = _solver(Assembler(dom), la)
                        U = solver.solve(F)
                    self.assertEqual([r for r in registros if AVISO_TANGENTE in r], [])
                    self.assertTrue(solver.reached_max_lambda)
                    self.assertTrue(np.all(np.isfinite(U)))

    def test_iterativo_avisa_una_vez(self):
        dom, F = armadura(sin_diagonal=False, carga="horizontal")
        with _avisos() as registros:
            _solver(Assembler(dom), "iterative").solve(F)
        aviso = [r for r in registros if AVISO_ITERATIVO in r]
        self.assertEqual(len(aviso), 1, registros)
        self.assertIn("DissipationArcLengthSolver", aviso[0])


if __name__ == "__main__":
    unittest.main()
