"""Estado interno por arreglos y vistas compatibles con ``ElementState`` (ADR 0014).

``FamilyState`` guarda, para todos los puntos de Gauss de una familia,
cuatro arreglos: variables internas committed y trial ``(P, n_state)`` y
esfuerzos committed y trial ``(P, n_sigma)``, con ``P = N·n_gp``. Es la
misma semántica trial/commit de ``ElementState`` (el commit copia trial
sobre committed; un paso rechazado no hace nada porque el siguiente
ensamblaje sobrescribe el trial), pero sin diccionarios ni copias
profundas.

``BatchedElementState`` es lo que ve cada elemento de la familia en
``elem.state``: las listas ``vars``, ``vars_trial``, ``stresses`` y
``stresses_trial`` son vistas que construyen el diccionario (o el
vector) bajo demanda a partir de la fila correspondiente, y asignar en
ellas escribe la fila. Así el camino por elemento, ``compute_gauss_state``,
el exportador VTK y los tests que leen el estado siguen funcionando sin
cambios, y el coste del diccionario sólo se paga cuando alguien lo pide.
"""
from __future__ import annotations

import numpy as np

from solidum.core.element_state import ElementState
from solidum.math.batch.schema import StateSchema


class FamilyState:
    """Arreglos de estado de una familia."""

    def __init__(self, n_elem: int, n_gp: int, schema: StateSchema,
                 n_sigma: int, initial_row: np.ndarray | None = None):
        self.n_elem = int(n_elem)
        self.n_gp = int(n_gp)
        self.schema = schema
        self.n_sigma = int(n_sigma)
        P = self.n_elem * self.n_gp
        n_state = schema.n_state
        if initial_row is None:
            initial_row = np.zeros(n_state)
        self.initial_row = np.ascontiguousarray(initial_row, dtype=np.float64).reshape(n_state)
        self.S_committed = np.tile(self.initial_row, (P, 1)).astype(np.float64)
        self.S_trial = self.S_committed.copy()
        self.sig_committed = np.zeros((P, self.n_sigma))
        self.sig_trial = np.zeros((P, self.n_sigma))
        self.flags = np.zeros(P, dtype=np.int8)

    @property
    def nbytes(self) -> int:
        return (self.S_committed.nbytes + self.S_trial.nbytes
                + self.sig_committed.nbytes + self.sig_trial.nbytes
                + self.flags.nbytes)

    def commit(self) -> None:
        """Trial → committed para toda la familia (copia, no intercambio:
        tras el commit ``vars_trial`` y ``vars`` coinciden, como en
        ``ElementState``)."""
        np.copyto(self.S_committed, self.S_trial)
        np.copyto(self.sig_committed, self.sig_trial)

    def rows(self, elem_index: int) -> slice:
        r0 = elem_index * self.n_gp
        return slice(r0, r0 + self.n_gp)

    def commit_element(self, elem_index: int) -> None:
        r = self.rows(elem_index)
        self.S_committed[r] = self.S_trial[r]
        self.sig_committed[r] = self.sig_trial[r]


class _RowView:
    """Lista de longitud ``n_gp`` que lee y escribe filas de un arreglo."""

    __slots__ = ("_fs", "_attr", "_row0", "_n", "_schema")

    def __init__(self, family_state: FamilyState, attr: str, elem_index: int,
                 schema: StateSchema | None):
        self._fs = family_state
        self._attr = attr
        self._row0 = elem_index * family_state.n_gp
        self._n = family_state.n_gp
        self._schema = schema

    def _array(self) -> np.ndarray:
        return getattr(self._fs, self._attr)

    def _index(self, idx: int) -> int:
        if not isinstance(idx, (int, np.integer)):
            raise TypeError("Las vistas de estado por lotes sólo admiten índices enteros.")
        if idx < 0:
            idx += self._n
        if not 0 <= idx < self._n:
            raise IndexError(f"Índice de punto de Gauss {idx} fuera de rango (0..{self._n - 1}).")
        return self._row0 + idx

    def __len__(self) -> int:
        return self._n

    def __iter__(self):
        for i in range(self._n):
            yield self[i]

    def __getitem__(self, idx):
        row = self._array()[self._index(idx)]
        if self._schema is not None:
            return self._schema.unpack(row)
        return np.array(row, dtype=np.float64)

    def __setitem__(self, idx, value) -> None:
        row = self._array()[self._index(idx)]
        if self._schema is not None:
            self._schema.pack(value, row, self._fs.initial_row)
        else:
            if value is None:
                row[:] = 0.0
            else:
                row[:] = np.asarray(value, dtype=np.float64).reshape(-1)

    def __repr__(self) -> str:
        return f"<vista {self._attr} de {self._n} puntos de Gauss>"


class BatchedElementState(ElementState):
    """``ElementState`` cuyo almacenamiento vive en un ``FamilyState``."""

    def __init__(self, family_state: FamilyState, elem_index: int):
        # Sin super().__init__: no hay listas propias.
        self.num_ip = family_state.n_gp
        self.family_state = family_state
        self.elem_index = int(elem_index)
        schema = family_state.schema
        self.vars = _RowView(family_state, "S_committed", elem_index, schema)
        self.vars_trial = _RowView(family_state, "S_trial", elem_index, schema)
        self.stresses = _RowView(family_state, "sig_committed", elem_index, None)
        self.stresses_trial = _RowView(family_state, "sig_trial", elem_index, None)

    def commit(self) -> None:
        self.family_state.commit_element(self.elem_index)

    # ------------------------------------------------------------------
    # Conversión desde / hacia el estado por diccionarios
    # ------------------------------------------------------------------

    def adopt(self, plain: ElementState | None) -> None:
        """Copia en las filas el contenido de un ``ElementState`` clásico."""
        if plain is None:
            return
        n = min(plain.num_ip, self.num_ip)
        for i in range(n):
            self.vars[i] = plain.vars[i]
            self.vars_trial[i] = plain.vars_trial[i]
            if plain.stresses[i] is not None:
                self.stresses[i] = plain.stresses[i]
            if plain.stresses_trial[i] is not None:
                self.stresses_trial[i] = plain.stresses_trial[i]

    def to_plain(self) -> ElementState:
        """Materializa un ``ElementState`` clásico con el contenido actual."""
        plain = ElementState(self.num_ip)
        plain.vars = [self.vars[i] for i in range(self.num_ip)]
        plain.vars_trial = [self.vars_trial[i] for i in range(self.num_ip)]
        plain.stresses = [self.stresses[i] for i in range(self.num_ip)]
        plain.stresses_trial = [self.stresses_trial[i] for i in range(self.num_ip)]
        return plain
