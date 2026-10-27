"""Expresiones de las cargas del YAML (ADR 0023).

Un valor de carga del YAML es un número o una cadena con una expresión de la
posición ``x``, ``y``, ``z`` y, en las cargas de contorno, de la normal
exterior ``nx``, ``ny``, ``nz``::

    body_force: [0.0, "-78.5e3*(1 + 0.1*y)"]
    tractions_by_group:
      Agujero: {pressure: "1e3*(1 - y)"}

La expresión **no se evalúa con** ``eval`` **sobre texto arbitrario**: se
analiza con :mod:`ast` y sólo se admiten números, esas variables, las
operaciones ``+ − * / **``, el signo, la constante ``pi`` y las funciones de
:data:`FUNCTIONS`. Cualquier otra cosa (atributos, índices, llamadas a otras
funciones, nombres desconocidos) es un error al leer el archivo, con la
expresión y lo que no se admite. Así un YAML ajeno no puede ejecutar código,
y una errata (``w`` por ``y``) se detecta antes de resolver.

Las constantes enteras se convierten a reales: ``9**9**9`` desborda de
inmediato (``OverflowError``) en vez de colgar el programa calculando un
entero enorme.
"""
from __future__ import annotations

import ast
import math
from typing import Callable

import numpy as np

#: Funciones admitidas en una expresión.
FUNCTIONS = {
    "sin": math.sin, "cos": math.cos, "tan": math.tan,
    "asin": math.asin, "acos": math.acos, "atan": math.atan, "atan2": math.atan2,
    "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
    "exp": math.exp, "log": math.log, "log10": math.log10, "sqrt": math.sqrt,
    "abs": abs, "min": min, "max": max,
}
#: Constantes admitidas.
CONSTANTS = {"pi": math.pi}
#: Variables de posición y, en las cargas de contorno, de la normal exterior.
POSITION = ("x", "y", "z")
NORMAL = ("nx", "ny", "nz")

_BINARY = (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)
_UNARY = (ast.UAdd, ast.USub)


class ExpressionError(ValueError):
    """Expresión de carga no admitida o que no se puede evaluar."""


class _AsFloat(ast.NodeTransformer):
    def visit_Constant(self, node):
        return ast.copy_location(ast.Constant(float(node.value)), node)


def _check(node, allowed: set, text: str) -> None:
    """Recorre el árbol y rechaza todo lo que no esté en la lista blanca."""
    def bad(what):
        raise ExpressionError(f"Expresión {text!r}: no se admite {what}. Se admiten números, "
                              f"{', '.join(sorted(allowed))}, + - * / **, pi y las funciones "
                              f"{', '.join(sorted(FUNCTIONS))}.")

    if isinstance(node, ast.Expression):
        _check(node.body, allowed, text)
    elif isinstance(node, ast.Constant):
        if isinstance(node.value, bool) or not isinstance(node.value, (int, float)):
            bad(f"la constante {node.value!r}")
    elif isinstance(node, ast.Name):
        if node.id not in allowed and node.id not in CONSTANTS:
            bad(f"el nombre {node.id!r}")
    elif isinstance(node, ast.BinOp):
        if not isinstance(node.op, _BINARY):
            bad(f"el operador {type(node.op).__name__}")
        _check(node.left, allowed, text)
        _check(node.right, allowed, text)
    elif isinstance(node, ast.UnaryOp):
        if not isinstance(node.op, _UNARY):
            bad(f"el operador {type(node.op).__name__}")
        _check(node.operand, allowed, text)
    elif isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name) or node.func.id not in FUNCTIONS:
            bad("la llamada a " + (repr(node.func.id) if isinstance(node.func, ast.Name)
                                   else "una función que no es de la lista"))
        if node.keywords or not node.args:
            bad(f"la llamada a {node.func.id!r} con argumentos con nombre o sin argumentos")
        for arg in node.args:
            if isinstance(arg, ast.Starred):
                bad("el argumento *")
            _check(arg, allowed, text)
    else:
        bad(f"la construcción {type(node).__name__}")


def compile_expression(text: str, variables=POSITION) -> Callable[..., float]:
    """Función ``f(*variables) -> float`` de la expresión ``text``.

    Lanza :class:`ExpressionError` si la expresión no es sintácticamente
    válida o usa algo que no está en la lista blanca."""
    if not isinstance(text, str):
        raise ExpressionError(f"Se esperaba una expresión de texto; se recibió {text!r}.")
    try:
        tree = ast.parse(text.strip(), mode="eval")
    except SyntaxError as exc:
        raise ExpressionError(f"Expresión {text!r}: sintaxis no válida ({exc.msg}).") from None
    _check(tree, set(variables), text)
    body = _AsFloat().visit(tree).body
    args = ast.arguments(posonlyargs=[], args=[ast.arg(arg=v) for v in variables],
                         kwonlyargs=[], kw_defaults=[], defaults=[])
    lam = ast.fix_missing_locations(ast.Expression(ast.Lambda(args=args, body=body)))
    fn = eval(compile(lam, "<expresión YAML>", "eval"),  # noqa: S307 — árbol ya validado
              {"__builtins__": {}, **FUNCTIONS, **CONSTANTS})
    return fn


def _padded(v, n: int = 3) -> tuple:
    """Las 2 ó 3 componentes de un punto o de una normal, completadas con 0."""
    v = np.asarray(v, dtype=np.float64).reshape(-1)
    return tuple(float(v[k]) if k < v.size else 0.0 for k in range(n))


def load_field(values, n_components: int, *, normal: bool, context: str):
    """Valor de una carga del YAML: el número o el vector de números tal cual,
    o, si alguna componente es una expresión, el campo equivalente:
    ``f(x)`` (``normal=False``) o ``f(x, n)`` (``normal=True``), que devuelve
    un escalar (``n_components == 1``) o un vector.

    ``context`` encabeza los mensajes de error (``"tractions[3]"``...)."""
    scalar = n_components == 1
    items = [values] if scalar and not isinstance(values, (list, tuple)) else values
    if not isinstance(items, (list, tuple)) or len(items) != n_components:
        raise ExpressionError(f"{context}: se esperaban {n_components} componentes; se recibió "
                              f"{values!r}.")
    if not any(isinstance(v, str) for v in items):
        try:
            arr = np.array([float(v) for v in items])
        except (TypeError, ValueError):
            raise ExpressionError(f"{context}: {values!r} no es un número ni una expresión.") \
                from None
        return float(arr[0]) if scalar else arr
    names = POSITION + NORMAL if normal else POSITION
    parts = []
    for v in items:
        if isinstance(v, str):
            try:
                parts.append((compile_expression(v, names), v))
            except ExpressionError as exc:
                raise ExpressionError(f"{context}: {exc}") from None
        else:
            try:
                c = float(v)
            except (TypeError, ValueError):
                raise ExpressionError(f"{context}: {v!r} no es un número ni una expresión.") \
                    from None
            parts.append((lambda *_, c=c: c, None))

    def evaluate(x, n=None):
        args = _padded(x) + (_padded(n) if normal else ())
        out = np.empty(len(parts))
        for k, (f, text) in enumerate(parts):
            try:
                out[k] = f(*args)
            except (ArithmeticError, ValueError, TypeError) as exc:
                raise ExpressionError(f"{context}: la expresión {text!r} no se puede evaluar en "
                                      f"x = {_padded(x)}: {exc}.") from None
        return float(out[0]) if scalar else out

    if normal:
        return lambda x, n: evaluate(x, n)
    return lambda x: evaluate(x)
