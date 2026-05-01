"""Tests del módulo de logging configurable (ADR 0005).

Verifica que el logger raíz ``fenix`` se configura una sola vez, que los
niveles se propagan a los loggers hijos y que ``set_log_level`` silencia el
progreso normal sin afectar a ``WARNING``/``ERROR``.
"""

import logging
import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import fenix
from fenix.logging import get_logger, set_log_level


class TestFenixLogger(unittest.TestCase):

    def setUp(self):
        # Restablecer al nivel por defecto antes de cada test.
        set_log_level("INFO")

    def tearDown(self):
        set_log_level("INFO")

    def test_logger_raiz_existe_con_handler_unico(self):
        root = logging.getLogger("fenix")
        self.assertGreaterEqual(len(root.handlers), 1)
        # Llamar de nuevo a get_logger no debe duplicar handlers.
        get_logger("subsistema_x")
        get_logger("subsistema_y")
        self.assertEqual(len(root.handlers), 1)

    def test_logger_no_propaga_a_root(self):
        # Evita que aplicaciones host con logging propio reciban duplicados.
        root = logging.getLogger("fenix")
        self.assertFalse(root.propagate)

    def test_set_log_level_acepta_string_y_constante(self):
        set_log_level("WARNING")
        self.assertEqual(logging.getLogger("fenix").level, logging.WARNING)
        set_log_level(logging.ERROR)
        self.assertEqual(logging.getLogger("fenix").level, logging.ERROR)

    def test_set_log_level_invalido_levanta(self):
        with self.assertRaises(ValueError):
            set_log_level("FOOBAR")

    def test_get_logger_devuelve_hijo_bajo_namespace_fenix(self):
        log = get_logger("solvers")
        self.assertEqual(log.name, "fenix.solvers")
        # Idempotente: llamar dos veces devuelve el mismo logger.
        self.assertIs(log, get_logger("solvers"))

    def test_silenciar_info_pero_dejar_warning(self):
        log = get_logger("test_subsistema")
        with self.assertLogs("fenix", level="WARNING") as captured:
            set_log_level("WARNING")
            log.info("este INFO no debe aparecer")
            log.warning("este WARNING sí debe aparecer")
        # Solo el WARNING quedó capturado.
        self.assertEqual(len(captured.records), 1)
        self.assertEqual(captured.records[0].levelname, "WARNING")

    def test_set_log_level_publico_desde_fenix(self):
        # API pública re-exportada en fenix/__init__.py.
        self.assertTrue(hasattr(fenix, "set_log_level"))
        self.assertTrue(hasattr(fenix, "get_logger"))


if __name__ == "__main__":
    unittest.main()
