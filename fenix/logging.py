"""Logging configurable para Fenix FEM (ADR 0005).

Logger jerárquico estándar bajo el nombre ``"fenix"``. Cada submódulo pide su
logger hijo con :func:`get_logger`; el usuario final silencia o detalla la
salida con :func:`set_log_level`. El comportamiento visible por defecto es
idéntico al de los antiguos ``print(...)``: nivel ``INFO`` a stdout con
formato corto ``[NIVEL] mensaje``.

No se llama a ``logging.basicConfig()``: la configuración se aplica solo al
logger ``"fenix"`` para no contaminar aplicaciones que embeban Fenix
(frontends GUI, notebooks).

Ejemplos
--------
Silenciar el progreso normal y dejar solo advertencias y errores::

    import fenix
    fenix.set_log_level("WARNING")

Capturar el log en un panel propio (frontend GUI)::

    import logging
    handler = MiHandlerDeGUI()
    logging.getLogger("fenix").addHandler(handler)

Detallar un subsistema concreto sin tocar el resto::

    logging.getLogger("fenix.solvers").setLevel(logging.DEBUG)
"""

from __future__ import annotations

import logging
import sys
from typing import Union

_ROOT_NAME = "fenix"
_DEFAULT_LEVEL = logging.INFO
_FORMAT = "[%(levelname)s] %(message)s"

_configured = False


def _configure_root_once() -> logging.Logger:
    """Configura el logger raíz ``fenix`` una sola vez por proceso.

    Idempotente: llamadas posteriores devuelven el mismo logger sin duplicar
    handlers. ``propagate=False`` evita que los mensajes asciendan al logger
    raíz de la aplicación host.
    """
    global _configured
    root = logging.getLogger(_ROOT_NAME)
    if _configured:
        return root

    handler = logging.StreamHandler(stream=sys.stdout)
    handler.setFormatter(logging.Formatter(_FORMAT))
    root.addHandler(handler)
    root.setLevel(_DEFAULT_LEVEL)
    root.propagate = False
    _configured = True
    return root


def get_logger(name: str) -> logging.Logger:
    """Devuelve el logger hijo ``fenix.<name>`` con la configuración aplicada.

    Parameters
    ----------
    name
        Nombre del subsistema (``"solvers"``, ``"parsers.yaml"``, ``"elements"``).
        Se interpreta como sufijo bajo el namespace ``"fenix"``.
    """
    _configure_root_once()
    if name.startswith(_ROOT_NAME + ".") or name == _ROOT_NAME:
        return logging.getLogger(name)
    return logging.getLogger(f"{_ROOT_NAME}.{name}")


def set_log_level(level: Union[int, str]) -> None:
    """Ajusta el nivel del logger raíz ``fenix``.

    Acepta tanto la constante de :mod:`logging` (``logging.WARNING``) como su
    nombre (``"WARNING"``, case-insensitive).

    Niveles relevantes:
      - ``"DEBUG"``: máximo detalle (slot reservado para diagnóstico futuro).
      - ``"INFO"``: progreso normal del solver, lectura de modelos, convergencia.
      - ``"WARNING"``: integración reducida, fallback Cholesky→LU, validaciones.
      - ``"ERROR"``: matriz singular, divergencia, raíces imaginarias.
    """
    root = _configure_root_once()
    if isinstance(level, str):
        level = level.upper()
        numeric = logging.getLevelName(level)
        if not isinstance(numeric, int):
            raise ValueError(f"Nivel de log desconocido: {level!r}")
        level = numeric
    root.setLevel(level)
