# solidum_fem/solidum/core/load_vector.py
"""Vector de cargas que lleva sus cargas equivalentes por elemento (ADR 0027).

Las fuerzas internas de extremo de un elemento estructural con carga
distribuida son ``F_int − f_eq``, donde ``f_eq`` es la carga nodal equivalente
de esa carga en el elemento (ADR 0002). Para descontarla hay que saber **qué
cargas distribuidas contiene el vector aplicado y con qué factor**. Antes, el
ensamblador recordaba sólo la última carga de cuerpo ensamblada, y dos cargas
de cuerpo en una misma corrida daban fuerzas internas erróneas sin aviso
(deuda #35).

:class:`LoadVector` resuelve eso en el dato mismo: es un ``np.ndarray`` (se usa
igual que antes) que lleva un :class:`EquivalentLoads` con las cargas
equivalentes de su parte distribuida, y las combina linealmente con la
aritmética del vector:

- ``F1 + F2``, ``F1 − F2`` y ``−F``: suman, restan o cambian el signo de sus
  cargas equivalentes. Un ``np.ndarray`` simple (cargas puntuales) no aporta.
- ``a·F`` y ``F/a`` con ``a`` escalar: las escalan.
- Cualquier otra operación (producto elemento a elemento, normas,
  comparaciones, indexado, ``np.asarray``) devuelve un ``np.ndarray`` simple:
  su resultado ya no es una carga con significado.

Así ``F = F_peso + λ·F_viento`` lleva ``f_eq = f_eq,peso + λ·f_eq,viento`` en
un análisis estático, en una etapa (ADR 0026) o en la función de carga de un
transitorio, sin que el usuario haga nada.

Representación: los vectores elementales se guardan concatenados en un solo
arreglo con una disposición (:class:`EquivalentLoadLayout`, ids y offsets)
compartida por todas las cargas de un mismo ensamblador. Combinar dos cargas
con la misma disposición es una operación vectorizada sobre ese arreglo.
"""
from __future__ import annotations

import numpy as np


class EquivalentLoadLayout:
    """Orden de los elementos y tamaño de su vector en un :class:`EquivalentLoads`."""

    __slots__ = ("ids", "offsets", "_index")

    def __init__(self, element_ids, sizes):
        self.ids = np.asarray(list(element_ids), dtype=np.int64)
        sizes = np.asarray(list(sizes), dtype=np.int64)
        if self.ids.shape != sizes.shape:
            raise ValueError("EquivalentLoadLayout: ids y tamaños de distinta longitud.")
        self.offsets = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
        self._index = {int(eid): k for k, eid in enumerate(self.ids)}

    @property
    def size(self) -> int:
        return int(self.offsets[-1])

    def position(self, element_id: int) -> int | None:
        return self._index.get(int(element_id))

    def same_as(self, other: "EquivalentLoadLayout") -> bool:
        return other is self or (
            np.array_equal(self.ids, other.ids) and np.array_equal(self.offsets, other.offsets))

    def __reduce__(self):
        sizes = np.diff(self.offsets)
        return (EquivalentLoadLayout, (self.ids.tolist(), sizes.tolist()))


class EquivalentLoads:
    """Cargas nodales equivalentes por elemento de la parte distribuida de una
    carga. Inmutable: toda combinación crea un objeto nuevo."""

    __slots__ = ("layout", "values")

    def __init__(self, layout: EquivalentLoadLayout, values):
        values = np.asarray(values, dtype=float)
        if values.shape != (layout.size,):
            raise ValueError(
                f"EquivalentLoads: {values.size} valores para una disposición de "
                f"{layout.size}.")
        self.layout = layout
        self.values = values

    # --- construcción -------------------------------------------------------

    @classmethod
    def from_dict(cls, loads: dict) -> "EquivalentLoads | None":
        """Desde ``{elem_id: vector}`` (el formato de ``SolveResult.equivalent_loads``)."""
        if not loads:
            return None
        ids = list(loads)
        blocks = [np.asarray(loads[eid], dtype=float).ravel() for eid in ids]
        return cls(EquivalentLoadLayout(ids, [b.size for b in blocks]), np.concatenate(blocks))

    # --- consulta -----------------------------------------------------------

    def get(self, element_id: int) -> np.ndarray | None:
        k = self.layout.position(element_id)
        if k is None:
            return None
        return self.values[self.layout.offsets[k]:self.layout.offsets[k + 1]]

    def as_dict(self) -> dict[int, np.ndarray]:
        off = self.layout.offsets
        return {int(eid): self.values[off[k]:off[k + 1]].copy()
                for k, eid in enumerate(self.layout.ids)}

    # --- álgebra ------------------------------------------------------------

    def scaled(self, factor: float) -> "EquivalentLoads":
        return EquivalentLoads(self.layout, float(factor) * self.values)

    def plus(self, other: "EquivalentLoads", factor: float = 1.0) -> "EquivalentLoads":
        """``self + factor·other``."""
        if self.layout.same_as(other.layout):
            return EquivalentLoads(self.layout, self.values + float(factor) * other.values)
        merged = self.as_dict()
        for eid, f in other.as_dict().items():
            if eid in merged and merged[eid].shape != f.shape:
                raise ValueError(
                    f"EquivalentLoads: el elemento {eid} tiene vectores de tamaños "
                    f"distintos ({merged[eid].size} y {f.size}); ¿cargas de modelos "
                    f"distintos?")
            merged[eid] = merged.get(eid, 0.0) + float(factor) * f
        return EquivalentLoads.from_dict(merged)


def _combine(a: EquivalentLoads | None, b: EquivalentLoads | None,
             factor_b: float) -> EquivalentLoads | None:
    if b is None:
        return a
    if a is None:
        return b.scaled(factor_b)
    return a.plus(b, factor_b)


def _eq_of(x) -> EquivalentLoads | None:
    return x.equivalent_loads if isinstance(x, LoadVector) else None


def _is_scalar(x) -> bool:
    return not isinstance(x, LoadVector) and np.ndim(x) == 0


_DROP = object()   # el resultado deja de ser una carga: np.ndarray simple


class LoadVector(np.ndarray):
    """Vector global de cargas con las cargas equivalentes de su parte
    distribuida (ADR 0027). Ver la documentación del módulo.

    Parameters
    ----------
    values : array_like, shape (ndof,)
        Vector global de fuerzas nodales.
    equivalent_loads : EquivalentLoads or None
        Cargas equivalentes por elemento de su parte distribuida; ``None`` si
        no tiene (sólo cargas puntuales).
    """

    def __new__(cls, values, equivalent_loads: EquivalentLoads | None = None):
        obj = np.array(values, dtype=float).view(cls)
        obj.equivalent_loads = equivalent_loads
        return obj

    def __array_finalize__(self, obj):
        self.equivalent_loads = getattr(obj, "equivalent_loads", None)

    def __getitem__(self, key):
        # Un trozo del vector no es una carga del modelo.
        return np.ndarray.__getitem__(self.view(np.ndarray), key)

    def __array_ufunc__(self, ufunc, method, *inputs, out=None, **kwargs):
        plain = [x.view(np.ndarray) if isinstance(x, LoadVector) else x for x in inputs]
        if out is not None:
            kwargs["out"] = tuple(o.view(np.ndarray) if isinstance(o, LoadVector) else o
                                  for o in out)
        result = getattr(ufunc, method)(*plain, **kwargs)
        if result is NotImplemented or method != "__call__":
            return result
        eq = self._equivalent_of(ufunc, inputs)
        if out is not None:
            target = out[0]
            if isinstance(target, LoadVector):
                target.equivalent_loads = None if eq is _DROP else eq
                return target
            return result
        if eq is _DROP or not isinstance(result, np.ndarray) or result.ndim == 0:
            return result
        vec = result.view(LoadVector)
        vec.equivalent_loads = eq
        return vec

    @staticmethod
    def _equivalent_of(ufunc, inputs):
        if ufunc in (np.add, np.subtract) and len(inputs) == 2:
            a, b = inputs
            return _combine(_eq_of(a), _eq_of(b), 1.0 if ufunc is np.add else -1.0)
        if ufunc is np.multiply and len(inputs) == 2:
            a, b = inputs
            if isinstance(a, LoadVector) and _is_scalar(b):
                return None if a.equivalent_loads is None else a.equivalent_loads.scaled(float(b))
            if isinstance(b, LoadVector) and _is_scalar(a):
                return None if b.equivalent_loads is None else b.equivalent_loads.scaled(float(a))
            return _DROP
        if ufunc is np.true_divide and len(inputs) == 2:
            a, b = inputs
            if isinstance(a, LoadVector) and _is_scalar(b):
                return None if a.equivalent_loads is None else a.equivalent_loads.scaled(1.0 / float(b))
            return _DROP
        if ufunc is np.negative:
            eq = _eq_of(inputs[0])
            return None if eq is None else eq.scaled(-1.0)
        if ufunc is np.positive:
            return _eq_of(inputs[0])
        return _DROP

    # Copia y serialización conservan las cargas equivalentes.
    def __reduce__(self):
        constructor, args, state = super().__reduce__()
        return constructor, args, (state, self.equivalent_loads)

    def __setstate__(self, state):
        nd_state, equivalent_loads = state
        super().__setstate__(nd_state)
        self.equivalent_loads = equivalent_loads


def equivalent_loads_of(F) -> EquivalentLoads | None:
    """Cargas equivalentes de ``F`` si es un :class:`LoadVector`; ``None`` si no."""
    return F.equivalent_loads if isinstance(F, LoadVector) else None


def as_load_vector(F, equivalent_loads: dict | EquivalentLoads | None = None) -> LoadVector:
    """``F`` como :class:`LoadVector`. Si ya lo es, tal cual; si no, con las
    cargas equivalentes dadas (``{elem_id: vector}`` o :class:`EquivalentLoads`)."""
    if isinstance(F, LoadVector):
        return F
    if isinstance(equivalent_loads, dict):
        equivalent_loads = EquivalentLoads.from_dict(equivalent_loads)
    return LoadVector(F, equivalent_loads)
