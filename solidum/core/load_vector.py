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
con la misma disposición es una operación vectorizada sobre ese arreglo; con
disposiciones distintas (una carga de línea en parte de los elementos y el
peso propio en todos), un mapa de índices hacia la unión, en caché por par de
disposiciones.

Dos canales con la misma álgebra (ADR 0044 §3): el vector nodal equivalente de
cada elemento y las **muestras** de la carga a lo largo de cada elemento de
línea (:class:`LoadSamples`: sus valores en ``LINE_LOAD_POINTS``, en ejes
globales), con las que se calculan las fuerzas internas dentro del elemento.
"""
from __future__ import annotations

import weakref

import numpy as np


class EquivalentLoadLayout:
    """Orden de los elementos y tamaño de su bloque en un canal de cargas.
    ``cols`` (opcional): columnas de cada bloque (las componentes de las
    muestras de un elemento de línea)."""

    __slots__ = ("ids", "offsets", "cols", "_index", "_union_cache", "__weakref__")

    def __init__(self, element_ids, sizes, cols=None):
        self.ids = np.asarray(list(element_ids), dtype=np.int64)
        sizes = np.asarray(list(sizes), dtype=np.int64)
        if self.ids.shape != sizes.shape:
            raise ValueError("EquivalentLoadLayout: ids y tamaños de distinta longitud.")
        self.offsets = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
        self.cols = None if cols is None else np.asarray(list(cols), dtype=np.int64)
        self._index = {int(eid): k for k, eid in enumerate(self.ids)}
        # Por disposición con la que se combinó: se suelta cuando ésta muere.
        self._union_cache = weakref.WeakKeyDictionary()

    @property
    def size(self) -> int:
        return int(self.offsets[-1])

    @property
    def sizes(self) -> np.ndarray:
        return np.diff(self.offsets)

    def position(self, element_id: int) -> int | None:
        return self._index.get(int(element_id))

    def same_as(self, other: "EquivalentLoadLayout") -> bool:
        return other is self or (
            np.array_equal(self.ids, other.ids) and np.array_equal(self.offsets, other.offsets))

    def union_with(self, other: "EquivalentLoadLayout"):
        """``(unión, índices de other en la unión)``: la unión pone primero los
        bloques de ``self`` (sus índices son ``arange(self.size)``) y después
        los elementos de ``other`` que no están en ``self``. En caché por par
        de disposiciones."""
        hit = self._union_cache.get(other)
        if hit is not None:
            return hit
        pos = np.array([-1 if (p := self._index.get(int(e))) is None else p for e in other.ids],
                       dtype=np.int64)
        common = pos >= 0
        o_sizes, s_sizes = other.sizes, self.sizes
        bad = common & (s_sizes[np.where(common, pos, 0)] != o_sizes)
        if np.any(bad):
            eid = int(other.ids[np.argmax(bad)])
            raise ValueError(
                f"EquivalentLoads: el elemento {eid} tiene bloques de tamaños distintos en las "
                f"dos cargas; ¿cargas de modelos distintos?")
        new = ~common
        u_ids = np.concatenate([self.ids, other.ids[new]])
        u_sizes = np.concatenate([s_sizes, o_sizes[new]])
        u_cols = None
        if self.cols is not None or other.cols is not None:
            s_cols = self.cols if self.cols is not None else np.ones_like(s_sizes)
            o_cols = other.cols if other.cols is not None else np.ones_like(o_sizes)
            u_cols = np.concatenate([s_cols, o_cols[new]])
        union = EquivalentLoadLayout(u_ids, u_sizes, u_cols)
        # Bloque de other en la unión: el de self si es común; si no, el nuevo.
        start = np.empty(len(other.ids), dtype=np.int64)
        start[common] = self.offsets[pos[common]]
        start[new] = union.offsets[len(self.ids) + np.arange(int(new.sum()))]
        idx = np.repeat(start - other.offsets[:-1], o_sizes) + np.arange(other.size)
        self._union_cache[other] = (union, idx)
        return union, idx

    def __reduce__(self):
        sizes = np.diff(self.offsets)
        return (EquivalentLoadLayout,
                (self.ids.tolist(), sizes.tolist(),
                 None if self.cols is None else self.cols.tolist()))


def _combine_channel(la: EquivalentLoadLayout, va: np.ndarray, lb: EquivalentLoadLayout,
                     vb: np.ndarray, factor: float):
    """``a + factor·b`` de dos canales: ``(disposición, valores)``."""
    if la.same_as(lb):
        return la, va + float(factor) * vb
    union, idx = la.union_with(lb)
    out = np.zeros(union.size)
    out[:la.size] = va
    out[idx] += float(factor) * vb
    return union, out


class LineLoadSamples:
    """Carga distribuida a lo largo de un elemento de línea en un resultado
    (``SolveResult.line_load_samples``, ADR 0044 §3): sus valores ``(n, d)``
    en ejes globales en las abscisas ``points`` (``ξ ∈ [0, 1]`` desde el nodo
    i, los ``LINE_LOAD_POINTS`` del elemento). Representa la carga por el
    polinomio que interpola esos valores (spec Frame2DEuler §14)."""

    __slots__ = ("points", "values")

    def __init__(self, points, values):
        self.points = np.asarray(points, dtype=float)
        self.values = np.atleast_2d(np.asarray(values, dtype=float))
        if self.values.shape[0] != self.points.size:
            raise ValueError("LineLoadSamples: tantas filas como abscisas.")

    def at(self, xi) -> np.ndarray:
        """Carga ``(m, d)`` en las abscisas ``xi``, en ejes globales."""
        from solidum.core.element_forces import load_polynomials
        c = load_polynomials(self.points, self.values)
        xi = np.atleast_1d(np.asarray(xi, dtype=float))
        return np.vander(xi, len(self.points), increasing=True) @ c

    def __repr__(self) -> str:
        return f"LineLoadSamples(points={self.points!r}, values={self.values!r})"


class LoadSamples:
    """Muestras de la carga a lo largo de cada elemento de línea (ADR 0044 §3):
    por elemento, sus valores ``(n, d)`` en ``LINE_LOAD_POINTS`` en ejes
    globales. Inmutable: toda combinación crea un objeto nuevo."""

    __slots__ = ("layout", "values")

    def __init__(self, layout: EquivalentLoadLayout, values):
        values = np.asarray(values, dtype=float)
        if values.shape != (layout.size,) or layout.cols is None:
            raise ValueError("LoadSamples: valores o columnas incompatibles con la disposición.")
        self.layout = layout
        self.values = values

    @classmethod
    def from_dict(cls, samples: dict) -> "LoadSamples | None":
        """Desde ``{elem_id: (n, d)}``."""
        if not samples:
            return None
        ids = list(samples)
        # Arreglos (n, d) o LineLoadSamples (SolveResult.line_load_samples).
        blocks = [np.atleast_2d(np.asarray(getattr(samples[eid], "values", samples[eid]),
                                           dtype=float)) for eid in ids]
        return cls(EquivalentLoadLayout(ids, [b.size for b in blocks],
                                        [b.shape[1] for b in blocks]),
                   np.concatenate([b.ravel() for b in blocks]))

    def get(self, element_id: int) -> np.ndarray | None:
        """Muestras ``(n, d)`` del elemento, o ``None``."""
        k = self.layout.position(element_id)
        if k is None:
            return None
        off = self.layout.offsets
        return self.values[off[k]:off[k + 1]].reshape(-1, int(self.layout.cols[k]))

    def as_dict(self) -> dict[int, np.ndarray]:
        return {int(eid): self.get(int(eid)).copy() for eid in self.layout.ids}

    def scaled(self, factor: float) -> "LoadSamples":
        return LoadSamples(self.layout, float(factor) * self.values)

    def plus(self, other: "LoadSamples", factor: float = 1.0) -> "LoadSamples":
        """``self + factor·other``."""
        layout, values = _combine_channel(self.layout, self.values, other.layout,
                                          other.values, factor)
        return LoadSamples(layout, values)


def _combine_samples(a: "LoadSamples | None", b: "LoadSamples | None",
                     factor_b: float) -> "LoadSamples | None":
    if b is None:
        return a
    if a is None:
        return b.scaled(factor_b)
    return a.plus(b, factor_b)


class EquivalentLoads:
    """Cargas nodales equivalentes por elemento de la parte distribuida de una
    carga y, en ``samples``, las muestras de la carga a lo largo de los
    elementos de línea (:class:`LoadSamples`, ADR 0044). Inmutable: toda
    combinación crea un objeto nuevo."""

    __slots__ = ("layout", "values", "samples")

    def __init__(self, layout: EquivalentLoadLayout, values,
                 samples: "LoadSamples | None" = None):
        values = np.asarray(values, dtype=float)
        if values.shape != (layout.size,):
            raise ValueError(
                f"EquivalentLoads: {values.size} valores para una disposición de "
                f"{layout.size}.")
        self.layout = layout
        self.values = values
        self.samples = samples

    # --- construcción -------------------------------------------------------

    @classmethod
    def from_dict(cls, loads: dict, samples: dict | None = None) -> "EquivalentLoads | None":
        """Desde ``{elem_id: vector}`` (el formato de ``SolveResult.equivalent_loads``)
        y, opcionalmente, ``{elem_id: muestras (n, d)}``."""
        if not loads:
            return None
        ids = list(loads)
        blocks = [np.asarray(loads[eid], dtype=float).ravel() for eid in ids]
        return cls(EquivalentLoadLayout(ids, [b.size for b in blocks]), np.concatenate(blocks),
                   LoadSamples.from_dict(samples or {}))

    # --- consulta -----------------------------------------------------------

    def get(self, element_id: int) -> np.ndarray | None:
        k = self.layout.position(element_id)
        if k is None:
            return None
        return self.values[self.layout.offsets[k]:self.layout.offsets[k + 1]]

    def samples_of(self, element_id: int) -> np.ndarray | None:
        """Muestras ``(n, d)`` de la carga a lo largo del elemento, o ``None``."""
        return None if self.samples is None else self.samples.get(element_id)

    def as_dict(self) -> dict[int, np.ndarray]:
        off = self.layout.offsets
        return {int(eid): self.values[off[k]:off[k + 1]].copy()
                for k, eid in enumerate(self.layout.ids)}

    def samples_dict(self) -> dict[int, np.ndarray]:
        return {} if self.samples is None else self.samples.as_dict()

    # --- álgebra ------------------------------------------------------------

    def scaled(self, factor: float) -> "EquivalentLoads":
        return EquivalentLoads(self.layout, float(factor) * self.values,
                               None if self.samples is None else self.samples.scaled(factor))

    def plus(self, other: "EquivalentLoads", factor: float = 1.0) -> "EquivalentLoads":
        """``self + factor·other``, en los dos canales."""
        layout, values = _combine_channel(self.layout, self.values, other.layout,
                                          other.values, factor)
        return EquivalentLoads(layout, values,
                               _combine_samples(self.samples, other.samples, factor))


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


def as_load_vector(F, equivalent_loads: dict | EquivalentLoads | None = None,
                   samples: dict | None = None) -> LoadVector:
    """``F`` como :class:`LoadVector`. Si ya lo es, tal cual; si no, con las
    cargas equivalentes dadas (``{elem_id: vector}`` o :class:`EquivalentLoads`)
    y, con el diccionario, las muestras de la carga a lo largo de los elementos
    de línea (``{elem_id: (n, d)}``, como ``SolveResult.line_load_samples``)."""
    if isinstance(F, LoadVector):
        return F
    if isinstance(equivalent_loads, dict):
        equivalent_loads = EquivalentLoads.from_dict(equivalent_loads, samples)
    return LoadVector(F, equivalent_loads)
