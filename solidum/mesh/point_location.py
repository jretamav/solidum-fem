"""Localización de puntos en los elementos isoparamétricos (ADR 0032).

La inversa del mapeo isoparamétrico ``x = Σ N_i(ξ)·X_i`` se resuelve por
Newton con lo que cada elemento ya declara: su celda de referencia
(``REFERENCE_CELL``) y sus funciones de forma (``batch_shape_functions``).
Sirve así a todo sólido 2D/3D del programa principal y a los de los módulos
de usuario que declaren lo mismo, sin código por tipo de elemento.

El jacobiano ``∂x/∂ξ`` se toma por diferencias centradas de las funciones
de forma: es exacto (salvo redondeo) en los elementos de Solidum, cuyas
funciones son a lo sumo cuadráticas en cada coordenada natural, y basta a
Newton en cualquier otro.
"""
from __future__ import annotations

from typing import Iterable, Optional, Tuple

import numpy as np

# Dimensión de cada celda de referencia y su centro (punto de arranque de
# Newton). ``tri``/``tet``: coordenadas de área o volumen ξ_k ≥ 0, Σξ_k ≤ 1;
# ``quad``/``hex``: el cubo [-1, 1]^d.
CELL_DIMENSION = {"tri": 2, "quad": 2, "tet": 3, "hex": 3}
_CELL_CENTER = {
    "tri": (1.0 / 3.0, 1.0 / 3.0),
    "quad": (0.0, 0.0),
    "tet": (0.25, 0.25, 0.25),
    "hex": (0.0, 0.0, 0.0),
}
# Paso de las diferencias centradas en coordenadas naturales.
_FD_STEP = 1.0e-4
# Tolerancia de pertenencia a la celda de referencia, en coordenadas naturales.
INSIDE_TOL = 1.0e-9


def inside_reference_cell(cell: str, xi, tol: float = INSIDE_TOL) -> bool:
    """``True`` si ``ξ`` está en la celda de referencia ``cell`` (con ``tol``)."""
    xi = np.asarray(xi, dtype=float)
    if cell in ("quad", "hex"):
        return bool(np.all(np.abs(xi) <= 1.0 + tol))
    if cell in ("tri", "tet"):
        return bool(np.all(xi >= -tol) and xi.sum() <= 1.0 + tol)
    raise ValueError(f"celda de referencia desconocida: {cell!r}")


def shape_values(element, xi) -> np.ndarray:
    """Funciones de forma ``N_i(ξ)`` del elemento en un punto natural."""
    N = element.batch_shape_functions(np.asarray(xi, dtype=float).reshape(1, -1))
    if N is None:
        raise ValueError(
            f"{type(element).__name__} no declara sus funciones de forma "
            f"(BATCH_SHAPE_FUNCTIONS o batch_shape_functions): no se pueden "
            f"localizar puntos en él."
        )
    return np.asarray(N[0], dtype=float)


def is_locatable(element, ndim: int) -> bool:
    """``True`` si el elemento es un sólido isoparamétrico de dimensión
    ``ndim`` con celda de referencia y funciones de forma."""
    cell = getattr(element, "REFERENCE_CELL", None)
    if CELL_DIMENSION.get(cell) != ndim:
        return False
    try:
        shape_values(element, _CELL_CENTER[cell])
    except (ValueError, TypeError, NotImplementedError):
        return False
    return True


def natural_coordinates(element, x, *, coords: Optional[np.ndarray] = None,
                        max_iter: int = 30) -> Optional[np.ndarray]:
    """Coordenadas naturales ``ξ`` del punto físico ``x`` en ``element``.

    Newton sobre ``Σ N_i(ξ)·X_i − x = 0`` desde el centro de la celda. Devuelve
    ``None`` si no converge (un punto muy alejado de un elemento distorsionado);
    la pertenencia a la celda la decide :func:`inside_reference_cell`.
    """
    cell = element.REFERENCE_CELL
    d = CELL_DIMENSION[cell]
    X = element.get_coordinate_matrix(ndim=d) if coords is None else coords
    x = np.asarray(x, dtype=float)[:d]
    size = float(np.max(X.max(axis=0) - X.min(axis=0)))
    xi = np.array(_CELL_CENTER[cell], dtype=float)
    eye = np.eye(d)
    r = shape_values(element, xi) @ X - x
    for _ in range(max_iter):
        if np.linalg.norm(r) <= 1.0e-13 * size:
            return xi
        J = np.empty((d, d))
        for k in range(d):
            Np = shape_values(element, xi + _FD_STEP * eye[k])
            Nm = shape_values(element, xi - _FD_STEP * eye[k])
            J[:, k] = (Np - Nm) @ X / (2.0 * _FD_STEP)
        try:
            xi = xi - np.linalg.solve(J, r)
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(xi)) or np.linalg.norm(xi) > 1.0e3:
            return None
        r = shape_values(element, xi) @ X - x
    return xi if np.linalg.norm(r) <= 1.0e-10 * size else None


class PointLocator:
    """Busca el elemento que contiene un punto entre unos elementos dados.

    Los candidatos se filtran por su caja envolvente (ampliada un 5 % para
    los lados curvos de los elementos de orden superior) y se prueban en
    orden de id: un punto en el borde común de dos elementos se asigna
    siempre al de menor id.
    """

    def __init__(self, elements: Iterable, ndim: int):
        self.ndim = ndim
        self.elements = sorted(elements, key=lambda e: e.id)
        if not self.elements:
            raise ValueError("PointLocator: no hay elementos donde localizar puntos.")
        self._coords = [e.get_coordinate_matrix(ndim=ndim) for e in self.elements]
        lo = np.array([c.min(axis=0) for c in self._coords])
        hi = np.array([c.max(axis=0) for c in self._coords])
        self.sizes = np.linalg.norm(hi - lo, axis=1)
        pad = 0.05 * self.sizes[:, None]
        self._lo, self._hi = lo - pad, hi + pad

    def natural_coordinates(self, index: int, x) -> Optional[np.ndarray]:
        """``ξ`` de ``x`` en el elemento de posición ``index`` si lo contiene."""
        el = self.elements[index]
        xi = natural_coordinates(el, x, coords=self._coords[index])
        if xi is None or not inside_reference_cell(el.REFERENCE_CELL, xi):
            return None
        return xi

    def locate(self, x) -> Optional[Tuple[int, np.ndarray]]:
        """``(índice, ξ)`` del primer elemento (por id) que contiene ``x``, o
        ``None`` si ninguno lo contiene."""
        x = np.asarray(x, dtype=float)[:self.ndim]
        candidates = np.nonzero(np.all((self._lo <= x) & (x <= self._hi), axis=1))[0]
        for index in candidates:
            xi = self.natural_coordinates(int(index), x)
            if xi is not None:
                return int(index), xi
        return None
