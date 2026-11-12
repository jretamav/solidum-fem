"""Calidad de la malla (ADR 0037).

Medidas por elemento de los sólidos 2D y 3D (celdas ``tri``, ``quad``,
``tet`` y ``hex``, de primer y segundo orden):

- **razón de jacobianos**: ``max det J / min det J`` en los nodos del
  elemento (vértices y nodos intermedios), con la cinemática del propio
  elemento (``BATCH_KINEMATICS``). Vale 1 en un paralelogramo, un
  paralelepípedo o un símplex de lados rectos con los nodos intermedios en el
  centro; crece con la distorsión. Si ``det J ≤ 0`` en algún nodo, el
  elemento está invertido o degenerado (``inf``);
- **razón de aspecto**: arista más larga / más corta entre vértices;
- **ángulos** de los vértices (2D): el mínimo y el máximo de cada elemento.

:func:`log_mesh_quality` resume la malla en el registro y avisa de los
elementos fuera de los límites de :data:`LIMITS`, que son los límites de
aviso por omisión de la comprobación de forma de ANSYS. Es un diagnóstico: no
detiene el análisis (un jacobiano no positivo lo detiene el ensamblaje).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit, types

from solidum.logging import get_logger
from solidum.math.batch.signatures import F64_2, F64_3, KIN_SIG

_log = get_logger("mesh")

#: Límites de aviso.
LIMITS = {
    "jacobian_ratio": 30.0,
    "aspect_ratio": 20.0,
    "max_angle_tri": 165.0,
    "max_angle_quad": 155.0,
}

_CORNERS = {"tri": 3, "quad": 4, "tet": 4, "hex": 8}
_EDGES = {
    "tri": ((0, 1), (1, 2), (2, 0)),
    "quad": ((0, 1), (1, 2), (2, 3), (3, 0)),
    "tet": ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)),
    "hex": ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7)),
}

_DETJ_SIG = types.void(types.FunctionType(KIN_SIG), F64_3, F64_2, F64_2, F64_2)


@njit(_DETJ_SIG, cache=True)
def _detj_at_points(kin, X, pts, B, out):
    for e in range(X.shape[0]):
        for p in range(pts.shape[0]):
            out[e, p] = kin(pts[p], X[e], B)


@dataclass
class MeshQuality:
    """Medidas por elemento, alineadas con ``element_ids``.

    ``jacobian_ratio`` es ``nan`` en los elementos sin cinemática compilada;
    ``min_angle`` y ``max_angle`` (grados), ``nan`` en 3D."""
    element_ids: np.ndarray
    cells: np.ndarray
    jacobian_ratio: np.ndarray
    aspect_ratio: np.ndarray
    min_angle: np.ndarray
    max_angle: np.ndarray

    @property
    def n_inverted(self) -> int:
        return int(np.sum(np.isinf(self.jacobian_ratio)))

    def summary(self) -> str:
        n = self.element_ids.size
        if n == 0:
            return "Calidad de la malla: sin elementos sólidos que medir."
        lines = [f"Calidad de la malla ({n} elementos sólidos):"]
        jr = self.jacobian_ratio[np.isfinite(self.jacobian_ratio)]
        if jr.size:
            lines.append(f"  razón de jacobianos (máx/mín de det J en los nodos): máx {jr.max():.3g}, "
                         f"media {jr.mean():.3g} (aviso > {LIMITS['jacobian_ratio']:g})")
        lines.append(f"  razón de aspecto (arista mayor/menor): máx {self.aspect_ratio.max():.3g}, "
                     f"media {self.aspect_ratio.mean():.3g} (aviso > {LIMITS['aspect_ratio']:g})")
        if np.any(np.isfinite(self.min_angle)):
            lines.append(f"  ángulos de los vértices: mín {np.nanmin(self.min_angle):.1f}°, "
                         f"máx {np.nanmax(self.max_angle):.1f}° (aviso > "
                         f"{LIMITS['max_angle_tri']:g}° en triángulos, "
                         f"> {LIMITS['max_angle_quad']:g}° en cuadriláteros)")
        lines.append(f"  elementos invertidos o degenerados (det J ≤ 0): {self.n_inverted}")
        return "\n".join(lines)

    def warnings(self) -> list[str]:
        out = []

        def report(mask, what, values, fmt):
            if np.any(mask):
                k = int(np.nanargmax(np.where(mask, values, -np.inf)))
                out.append(f"Malla: {int(mask.sum())} elemento(s) con {what} (el peor, el "
                           f"{int(self.element_ids[k])}: {fmt(values[k])}).")

        inverted = np.isinf(self.jacobian_ratio)
        if np.any(inverted):
            ids = self.element_ids[inverted]
            out.append(f"Malla: {ids.size} elemento(s) invertido(s) o degenerado(s) (det J ≤ 0 en "
                       f"algún nodo): {', '.join(map(str, ids[:10]))}"
                       + (" ..." if ids.size > 10 else "") + ".")
        jr = np.where(inverted, np.nan, self.jacobian_ratio)
        with np.errstate(invalid="ignore"):
            report(jr > LIMITS["jacobian_ratio"], "razón de jacobianos > "
                   f"{LIMITS['jacobian_ratio']:g}", jr, lambda v: f"{v:.3g}")
            report(self.aspect_ratio > LIMITS["aspect_ratio"], "razón de aspecto > "
                   f"{LIMITS['aspect_ratio']:g}", self.aspect_ratio, lambda v: f"{v:.3g}")
            limit = np.where(self.cells == "tri", LIMITS["max_angle_tri"],
                             LIMITS["max_angle_quad"])
            report(self.max_angle > limit, "un ángulo mayor que el límite (165° en triángulos, "
                   "155° en cuadriláteros)", self.max_angle, lambda v: f"{v:.1f}°")
        return out


def _corner_angles(P: np.ndarray) -> np.ndarray:
    """Ángulos (grados) en los vértices de los polígonos ``P`` (N, k, 2)."""
    prev = np.roll(P, 1, axis=1) - P
    nxt = np.roll(P, -1, axis=1) - P
    cos = np.sum(prev * nxt, axis=2) / (np.linalg.norm(prev, axis=2)
                                         * np.linalg.norm(nxt, axis=2))
    return np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))


def _detj_ratio(cls, X: np.ndarray) -> np.ndarray:
    kin = getattr(cls, "BATCH_KINEMATICS", None)
    rows = getattr(cls, "STRAIN_DIM", None) or getattr(cls, "FLUX_DIM", None)
    if kin is None or rows is None:
        return np.full(X.shape[0], np.nan)
    pts = np.ascontiguousarray(cls.REFERENCE_NODES, dtype=np.float64)
    n_dof = X.shape[1] * len(cls.DOF_NAMES)
    B = np.zeros((int(rows), n_dof))
    detj = np.zeros((X.shape[0], pts.shape[0]))
    _detj_at_points(kin, np.ascontiguousarray(X), pts, B, detj)
    lo, hi = detj.min(axis=1), detj.max(axis=1)
    with np.errstate(divide="ignore"):
        return np.where(lo > 0.0, hi / np.where(lo > 0.0, lo, 1.0), np.inf)


def mesh_quality(domain) -> MeshQuality:
    """Medidas de calidad de los sólidos 2D y 3D de ``domain``."""
    by_class: dict = {}
    for el in domain.elements.values():
        cell = getattr(el, "REFERENCE_CELL", None)
        if cell in _CORNERS and getattr(el, "REFERENCE_NODES", None) is not None:
            by_class.setdefault(type(el), []).append(el)
    parts = []
    for cls, els in by_class.items():
        cell = cls.REFERENCE_CELL
        d = 2 if cell in ("tri", "quad") else 3
        X = np.array([[list(nd.coordinates)[:d] for nd in el.nodes] for el in els],
                     dtype=np.float64)
        corners = X[:, :_CORNERS[cell]]
        lengths = np.stack([np.linalg.norm(corners[:, b] - corners[:, a], axis=1)
                            for a, b in _EDGES[cell]], axis=1)
        with np.errstate(divide="ignore"):
            aspect = lengths.max(axis=1) / lengths.min(axis=1)
        if cell in ("tri", "quad"):
            ang = _corner_angles(corners)
            amin, amax = ang.min(axis=1), ang.max(axis=1)
        else:
            amin = amax = np.full(len(els), np.nan)
        parts.append((np.array([el.id for el in els]), np.full(len(els), cell),
                      _detj_ratio(cls, X), aspect, amin, amax))
    if not parts:
        empty = np.zeros(0)
        return MeshQuality(empty.astype(int), empty.astype(str), empty, empty, empty, empty)
    cols = [np.concatenate(c) for c in zip(*parts)]
    order = np.argsort(cols[0], kind="stable")
    return MeshQuality(*(c[order] for c in cols))


def log_mesh_quality(domain) -> MeshQuality:
    """Calcula la calidad de la malla, la resume en el registro (INFO) y avisa
    (WARNING) de los elementos fuera de :data:`LIMITS`."""
    q = mesh_quality(domain)
    if q.element_ids.size:
        for line in q.summary().splitlines():
            _log.info(line)
        for w in q.warnings():
            _log.warning(w)
    return q
