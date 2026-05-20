"""Utilidades geométricas puras reutilizables entre elementos.

Funciones libres de estado, independientes de Element/Material/Solver. Sirven
de cajón único para operaciones geométricas que aparezcan duplicadas entre
módulos.
"""
from __future__ import annotations

import numpy as np


def perpendicular_projector(unit_vec: np.ndarray) -> np.ndarray:
    """Proyector ortogonal al subespacio generado por un vector unitario.

    Devuelve la matriz :math:`\\mathbf P = \\mathbf I - \\hat{\\mathbf e}\\hat{\\mathbf e}^\\top`
    de dimensión ``(n, n)`` donde ``n = len(unit_vec)``. Aplicada a un vector
    cualquiera, le quita su componente paralela a ``unit_vec`` y deja la
    perpendicular.

    Empleado por las cinemáticas corotacionales 3D (Truss3DCorot,
    Cable3DCorot) para construir la rigidez geométrica en el plano normal
    al eje del elemento.

    Parameters
    ----------
    unit_vec : np.ndarray
        Vector director (``shape=(n,)``). Se asume unitario; si no lo es,
        la matriz seguirá siendo simétrica pero ya no será un proyector
        idempotente. La normalización es responsabilidad del que llama.

    Returns
    -------
    np.ndarray
        Matriz proyectora ``(n, n)``.
    """
    return np.eye(unit_vec.size) - np.outer(unit_vec, unit_vec)
