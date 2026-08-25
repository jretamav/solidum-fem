"""Tests del módulo de logging configurable (ADR 0005).

Verifica que el logger raíz ``solidum`` se configura una sola vez, que los
niveles se propagan a los loggers hijos y que ``set_log_level`` silencia el
progreso normal sin afectar a ``WARNING``/``ERROR``.
"""

import logging
import os
import sys
import unittest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum
from solidum.logging import get_logger, set_log_level


class TestFenixLogger(unittest.TestCase):

    def setUp(self):
        # Restablecer al nivel por defecto antes de cada test.
        set_log_level("INFO")

    def tearDown(self):
        set_log_level("INFO")

    @staticmethod
    def _handlers_del_proyecto(root):
        """Handlers instalados por Solidum, excluyendo los de terceros.

        El logger ``solidum`` puede llevar handlers que no son suyos:
        ``_pytest.logging`` inyecta los propios (``LogCaptureHandler``,
        ``_LiveLoggingNullHandler``, ``_FileHandler``) en cuanto algún test
        de la suite engancha un handler a este logger, y una aplicación que
        embeba Solidum haría lo mismo — es un caso de uso **documentado** en
        el docstring de ``solidum/logging.py``.

        Contar todos convertiría la aserción en una afirmación sobre el
        entorno de ejecución en lugar de sobre el código del proyecto. Lo que
        interesa verificar es que ``_configure_root_once`` es idempotente:
        que Solidum instala **su** handler una sola vez.
        """
        return [h for h in root.handlers
                if type(h).__module__.split(".")[0] not in ("_pytest", "pytest")]

    def test_logger_raiz_existe_con_handler_unico(self):
        root = logging.getLogger("solidum")
        self.assertGreaterEqual(len(self._handlers_del_proyecto(root)), 1)
        # Llamar de nuevo a get_logger no debe duplicar handlers.
        get_logger("subsistema_x")
        get_logger("subsistema_y")
        self.assertEqual(len(self._handlers_del_proyecto(root)), 1)

    def test_logger_no_propaga_a_root(self):
        # Evita que aplicaciones host con logging propio reciban duplicados.
        root = logging.getLogger("solidum")
        self.assertFalse(root.propagate)

    def test_set_log_level_acepta_string_y_constante(self):
        set_log_level("WARNING")
        self.assertEqual(logging.getLogger("solidum").level, logging.WARNING)
        set_log_level(logging.ERROR)
        self.assertEqual(logging.getLogger("solidum").level, logging.ERROR)

    def test_set_log_level_invalido_levanta(self):
        with self.assertRaises(ValueError):
            set_log_level("FOOBAR")

    def test_get_logger_devuelve_hijo_bajo_namespace_fenix(self):
        log = get_logger("solvers")
        self.assertEqual(log.name, "solidum.solvers")
        # Idempotente: llamar dos veces devuelve el mismo logger.
        self.assertIs(log, get_logger("solvers"))

    def test_silenciar_info_pero_dejar_warning(self):
        log = get_logger("test_subsistema")
        with self.assertLogs("solidum", level="WARNING") as captured:
            set_log_level("WARNING")
            log.info("este INFO no debe aparecer")
            log.warning("este WARNING sí debe aparecer")
        # Solo el WARNING quedó capturado.
        self.assertEqual(len(captured.records), 1)
        self.assertEqual(captured.records[0].levelname, "WARNING")

    def test_set_log_level_publico_desde_fenix(self):
        # API pública re-exportada en solidum/__init__.py.
        self.assertTrue(hasattr(solidum, "set_log_level"))
        self.assertTrue(hasattr(solidum, "get_logger"))


if __name__ == "__main__":
    unittest.main()
