"""Elementos de línea rectos de dos nodos: celda de referencia, campo interior y
fuerza de cuerpo variable (ADR 0022, fase 4).

Lo comparten las armaduras y los marcos. Cada elemento aporta su
interpolación en ejes locales (``_line_interpolation``), el paso de sus DOFs
globales a los locales (``_to_local``) y, si reparte la carga en ejes
locales, los ejes (``_load_axes``) y el paso inverso (``_to_global``).
:class:`LineElementMixin` da con eso:

- ``REFERENCE_CELL = "line"`` y ``REFERENCE_NODES``, que usa el generador de
  mallas estructuradas;
- :meth:`~LineElementMixin.line_field`: desplazamientos, giros y deformaciones
  generalizadas en cualquier punto del elemento, en ejes locales; con
  :meth:`~LineElementMixin.section_rigidities`, la energía de deformación
  ``½ ∫ eᵀ·D·e dx`` que la rigidez del elemento integra exactamente;
- :meth:`~LineElementMixin._field_line_load`: el vector nodal consistente de
  una fuerza de cuerpo dada como campo ``b(x)``.
"""
from __future__ import annotations

from typing import ClassVar

import numpy as np

from solidum.math.distributed_loads import field_values

# Gauss-Legendre de 4 puntos en [0, 1]: exacto hasta grado 7, es decir
# N·q con N cúbica (Hermite) y una carga q de grado ≤ 4.
_LINE_XI, _LINE_W = np.polynomial.legendre.leggauss(4)
_LINE_XI = 0.5 * (_LINE_XI + 1.0)
_LINE_W = 0.5 * _LINE_W


def hermite_functions(xi: np.ndarray, L: float):
    """Funciones cúbicas de Hermite ``v = N₁v₁ + N₂θ₁ + N₃v₂ + N₄θ₂`` en
    ``x = ξ·L``, con sus derivadas primera y segunda respecto de ``x``, cada
    una ``(n, 4)``."""
    xi = np.asarray(xi, dtype=np.float64)
    N = np.column_stack([1.0 - 3.0 * xi**2 + 2.0 * xi**3, L * (xi - 2.0 * xi**2 + xi**3),
                         3.0 * xi**2 - 2.0 * xi**3, L * (xi**3 - xi**2)])
    dN = np.column_stack([6.0 * (xi**2 - xi) / L, 1.0 - 4.0 * xi + 3.0 * xi**2,
                          6.0 * (xi - xi**2) / L, 3.0 * xi**2 - 2.0 * xi])
    d2N = np.column_stack([(12.0 * xi - 6.0) / L**2, (6.0 * xi - 4.0) / L,
                           (6.0 - 12.0 * xi) / L**2, (6.0 * xi - 2.0) / L])
    return N, dN, d2N


def linear_tangent(material) -> float:
    """Módulo ``E`` de un material 1D elástico lineal: su tangente en el
    origen."""
    return float(material.compute_state(0.0, None)[1])


class LineElementMixin:
    """Elemento de línea recto de dos nodos (ADR 0022, fase 4).

    Contrato de las subclases
    -------------------------
    ``LINE_FIELDS`` : nombres de los desplazamientos y giros locales que da
        :meth:`line_field` (``"u"`` axial, ``"v"``, ``"w"`` transversales,
        ``"theta_x"``, ``"theta_y"``, ``"theta_z"`` giros); los
        desplazamientos van primero.
    ``LINE_STRAINS`` : nombres de las deformaciones generalizadas
        (``"eps"`` axial, ``"twist"`` = θx', ``"kappa_y"``, ``"kappa_z"``,
        ``"gamma_y"``), en el orden de :meth:`section_rigidities`.
    ``_line_interpolation(xi)`` : ``(N, B)`` de forma ``(n, k, m)`` y
        ``(n, s, m)``, con ``m`` los DOFs locales: campos y deformaciones en
        los puntos ``xi ∈ [0, 1]``.
    ``_to_local(u_e)`` : DOFs locales a partir del vector global del elemento.
    ``section_rigidities()`` : rigideces de la sección por deformación.
    """
    REFERENCE_CELL: ClassVar[str] = "line"
    REFERENCE_NODES: ClassVar[tuple] = ((-1.0,), (1.0,))
    LINE_FIELDS: ClassVar[tuple] = ()
    LINE_STRAINS: ClassVar[tuple] = ()

    def line_field(self, u_e: np.ndarray, xi) -> tuple[np.ndarray, np.ndarray]:
        """Campos ``(n, len(LINE_FIELDS))`` y deformaciones generalizadas
        ``(n, len(LINE_STRAINS))`` en los puntos ``xi ∈ [0, 1]`` (``ξ = 0`` en
        el nodo i), en ejes locales, con la interpolación del elemento y su
        vector global de DOFs ``u_e``."""
        N, B = self._line_interpolation(np.atleast_1d(np.asarray(xi, dtype=np.float64)))
        u = self._to_local(np.asarray(u_e, dtype=np.float64))
        return N @ u, B @ u

    def section_rigidities(self) -> np.ndarray:
        raise NotImplementedError

    def _line_ends(self) -> np.ndarray:
        """Coordenadas ``(2, d)`` de los nodos, con ``d`` las componentes de
        desplazamiento del elemento (las que faltan, cero)."""
        d = sum(1 for n in self.DOF_NAMES if n in ("ux", "uy", "uz"))
        X = np.zeros((2, d))
        for k, node in enumerate(self.nodes):
            c = np.asarray(node.coordinates, dtype=np.float64)[:d]
            X[k, :len(c)] = c
        return X

    def _field_load_by_component(self, b, area: float) -> np.ndarray:
        """Vector nodal ``∫ N_a·A·b(x) dx`` con las funciones lineales
        ``N_a`` aplicadas a cada componente global por separado, el reparto
        de las armaduras (también la componente perpendicular al eje, que
        es fuerza nodal pura)."""
        X = self._line_ends()
        x = X[0] + _LINE_XI[:, None] * (X[1] - X[0])
        values = field_values(b, x, X.shape[1], "de cuerpo")
        N = np.column_stack([1.0 - _LINE_XI, _LINE_XI])
        return (area * self.L0 * (N * _LINE_W[:, None]).T @ values).reshape(-1)

    def _field_line_load(self, b, area: float) -> np.ndarray:
        """Vector nodal ``∫ Nᵀ·q dx`` de la carga por unidad de longitud
        ``q = A·b(x)``, proyectada sobre los ejes locales (``_load_axes``) en
        cada punto de Gauss y repartida con la interpolación de los
        desplazamientos del elemento; en ejes globales."""
        X = self._line_ends()
        x = X[0] + _LINE_XI[:, None] * (X[1] - X[0])
        values = field_values(b, x, X.shape[1], "de cuerpo")
        R = self._load_axes()
        q = area * values @ R.T
        N, _ = self._line_interpolation(_LINE_XI)
        f_local = np.einsum("g,gki,gk->i", _LINE_W * self.L0, N[:, :R.shape[0], :], q)
        return self._to_global(f_local)
