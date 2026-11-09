"""Logging configurable para Solidum FEM (ADR 0005).

Logger jerárquico estándar bajo el nombre ``"solidum"``. Cada submódulo pide su
logger hijo con :func:`get_logger`; el usuario final silencia o detalla la
salida con :func:`set_log_level`. El comportamiento visible por defecto es
idéntico al de los antiguos ``print(...)``: nivel ``INFO`` a stdout con
formato corto ``[NIVEL] mensaje``.

No se llama a ``logging.basicConfig()``: la configuración se aplica solo al
logger ``"solidum"`` para no contaminar aplicaciones que embeban Solidum
(frontends GUI, notebooks).

Ejemplos
--------
Silenciar el progreso normal y dejar solo advertencias y errores::

    import solidum
    solidum.set_log_level("WARNING")

Capturar el log en un panel propio (frontend GUI)::

    import logging
    handler = MiHandlerDeGUI()
    logging.getLogger("solidum").addHandler(handler)

Detallar un subsistema concreto sin tocar el resto::

    logging.getLogger("solidum.solvers").setLevel(logging.DEBUG)
"""

from __future__ import annotations

import logging
import sys
import unicodedata
from typing import Union

_ROOT_NAME = "solidum"
_DEFAULT_LEVEL = logging.INFO
_FORMAT = "[%(levelname)s] %(message)s"

_configured = False

# Equivalentes en texto plano de los símbolos que usan los mensajes, para
# consolas cuya codificación no los tiene (cp1252 de Windows, ASCII). Las
# letras griegas se nombran a partir de Unicode (``Δ`` → ``Delta``); los
# subíndices y superíndices, por su descomposición (``₁`` → ``_1``, ``²`` →
# ``^2``); las letras acentuadas pierden el acento si hace falta (``á`` → ``a``).
_PLAIN_SYMBOLS = {
    "⇒": "=>", "→": "->", "←": "<-", "↔": "<->", "≤": "<=", "≥": ">=",
    "≈": "~", "≠": "!=", "≪": "<<", "≫": ">>", "−": "-", "‖": "||",
    "·": "*", "×": "x", "…": "...", "√": "sqrt", "∞": "inf", "∈": " in ",
    "∑": "sum", "∂": "d", "∇": "grad", "±": "+/-", "—": "--", "–": "-",
    "«": '"', "»": '"', "“": '"', "”": '"', "‘": "'", "’": "'", "°": " deg",
}


def _plain_char(ch: str, encoding: str) -> str:
    """Sustituto de ``ch`` que ``encoding`` sí puede escribir."""
    if ch in _PLAIN_SYMBOLS:
        return _PLAIN_SYMBOLS[ch]
    if unicodedata.combining(ch):
        return ""                      # tilde o acento suelto: κ̃ → kappa
    name = unicodedata.name(ch, "")
    if name.startswith("GREEK "):
        parts = name.split()
        if "LETTER" in parts:
            # "GREEK SMALL LETTER LAMDA", "... FINAL SIGMA", "... ALPHA WITH TONOS"
            words = [w for w in parts[parts.index("LETTER") + 1:] if w != "FINAL"]
            if words:
                letter = words[0].lower().replace("lamda", "lambda")
                return letter.capitalize() if "CAPITAL" in parts else letter
    decomposed = unicodedata.normalize("NFKD", ch)
    if decomposed != ch:
        plain = "".join(c for c in decomposed if not unicodedata.combining(c))
        prefix = "_" if name.startswith("SUBSCRIPT") else "^" if name.startswith("SUPERSCRIPT") else ""
        text = prefix + plain
        try:
            text.encode(encoding)
            return text
        except UnicodeEncodeError:
            return prefix + "".join(_plain_char(c, encoding) for c in plain)
    return "?"


def encodable_text(text: str, encoding: str | None) -> str:
    """``text`` tal cual si ``encoding`` puede escribirlo; si no, con cada
    carácter que no puede sustituido por su equivalente en texto plano
    (ADR 0025). Sin codificación conocida, ``text`` sin cambios."""
    if not encoding:
        return text
    try:
        text.encode(encoding)
        return text
    except (UnicodeEncodeError, LookupError):
        pass
    out = []
    for ch in text:
        try:
            ch.encode(encoding)
            out.append(ch)
        except (UnicodeEncodeError, LookupError):
            out.append(_plain_char(ch, encoding))
    return "".join(out)


class _ConsoleHandler(logging.StreamHandler):
    """``StreamHandler`` que no falla con la codificación de la consola.

    Con una salida en cp1252 (la consola de Git Bash en Windows, o cualquier
    salida encauzada con esa codificación) un mensaje con letras griegas,
    subíndices o flechas lanzaba ``UnicodeEncodeError`` dentro de
    ``logging``, que escribía una traza ``--- Logging error ---`` por mensaje
    aunque el análisis siguiera bien. Aquí el mensaje se adapta a la
    codificación del flujo antes de escribirlo (ADR 0025)."""

    def emit(self, record: logging.LogRecord) -> None:
        try:
            msg = self.format(record)
            msg = encodable_text(msg, getattr(self.stream, "encoding", None))
            self.stream.write(msg + self.terminator)
            self.flush()
        except RecursionError:
            raise
        except Exception:
            self.handleError(record)


def _configure_root_once() -> logging.Logger:
    """Configura el logger raíz ``solidum`` una sola vez por proceso.

    Idempotente: llamadas posteriores devuelven el mismo logger sin duplicar
    handlers. ``propagate=False`` evita que los mensajes asciendan al logger
    raíz de la aplicación host.
    """
    global _configured
    root = logging.getLogger(_ROOT_NAME)
    if _configured:
        return root

    handler = _ConsoleHandler(stream=sys.stdout)
    handler.setFormatter(logging.Formatter(_FORMAT))
    root.addHandler(handler)
    root.setLevel(_DEFAULT_LEVEL)
    root.propagate = False
    _configured = True
    return root


def get_logger(name: str) -> logging.Logger:
    """Devuelve el logger hijo ``solidum.<name>`` con la configuración aplicada.

    Parameters
    ----------
    name
        Nombre del subsistema (``"solvers"``, ``"parsers.yaml"``, ``"elements"``).
        Se interpreta como sufijo bajo el namespace ``"solidum"``.
    """
    _configure_root_once()
    if name.startswith(_ROOT_NAME + ".") or name == _ROOT_NAME:
        return logging.getLogger(name)
    return logging.getLogger(f"{_ROOT_NAME}.{name}")


def set_log_level(level: Union[int, str]) -> None:
    """Ajusta el nivel del logger raíz ``solidum``.

    Acepta tanto la constante de :mod:`logging` (``logging.WARNING``) como su
    nombre (``"WARNING"``, case-insensitive).

    Niveles (política del ADR 0025):
      - ``"DEBUG"``: detalle por iteración de Newton y por intento de paso.
      - ``"INFO"`` (por omisión): progreso normal; en los solvers
        incrementales, a lo sumo una línea por paso convergido y un resumen.
      - ``"WARNING"``: lo que el solver recupera por sí solo (bisección por
        tangente singular o raíces imaginarias, con su motivo), integración
        reducida, fallback Cholesky→LU, trazado detenido por ``max_steps``.
      - ``"ERROR"``: sólo lo que termina el análisis.
    """
    root = _configure_root_once()
    if isinstance(level, str):
        level = level.upper()
        numeric = logging.getLevelName(level)
        if not isinstance(numeric, int):
            raise ValueError(f"Nivel de log desconocido: {level!r}")
        level = numeric
    root.setLevel(level)
