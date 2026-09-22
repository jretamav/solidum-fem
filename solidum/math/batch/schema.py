"""Esquema de variables internas: de ``STATE_SCHEMA`` a filas de arreglo (ADR 0014).

Un material declara sus variables internas como ``{nombre: forma}``::

    STATE_SCHEMA = {"eps_p": (4,), "alpha": ()}

y este módulo deriva el almacenamiento por filas: ``eps_p`` ocupa las
columnas 0–3 y ``alpha`` la 4. El material no sabe que existe un
arreglo; su ``compute_state`` sigue recibiendo y devolviendo un ``dict``.
La traducción ``dict ↔ fila`` vive aquí y sólo se paga cuando alguien
pide el diccionario (post-proceso, camino por elemento).
"""
from __future__ import annotations

from typing import Mapping

import numpy as np


class StateSchema:
    """Disposición por columnas de un ``STATE_SCHEMA``.

    Parameters
    ----------
    schema
        ``{nombre: forma}`` con ``forma`` una tupla (``()`` para escalares).
        El orden de inserción fija el orden de las columnas.
    """

    def __init__(self, schema: Mapping[str, tuple]):
        names = []
        shapes = []
        sizes = []
        offsets = []
        total = 0
        for name, shape in schema.items():
            shape = tuple(int(s) for s in shape)
            if any(s <= 0 for s in shape):
                raise ValueError(
                    f"STATE_SCHEMA: la variable {name!r} declara forma {shape}; "
                    "todas las dimensiones deben ser positivas."
                )
            size = int(np.prod(shape)) if shape else 1
            names.append(str(name))
            shapes.append(shape)
            sizes.append(size)
            offsets.append(total)
            total += size
        self.names = tuple(names)
        self.shapes = tuple(shapes)
        self.sizes = tuple(sizes)
        self.offsets = tuple(offsets)
        self.n_state = total

    # ------------------------------------------------------------------

    def unpack(self, row: np.ndarray) -> dict | None:
        """Fila → ``dict``. ``None`` si el esquema está vacío (material sin
        historia), que es lo que hoy devuelven esos materiales."""
        if self.n_state == 0:
            return None
        out = {}
        for name, shape, size, off in zip(self.names, self.shapes,
                                          self.sizes, self.offsets):
            if shape == ():
                out[name] = float(row[off])
            else:
                out[name] = np.array(row[off:off + size], dtype=np.float64).reshape(shape)
        return out

    def pack(self, state: Mapping | None, row: np.ndarray,
             initial: np.ndarray | None = None) -> None:
        """``dict`` → fila, escrita in situ.

        ``None`` restaura la fila inicial (``initial``) o ceros: es lo que
        significa "sin estado" para un material que lo devuelve así.
        """
        if self.n_state == 0:
            return
        if state is None:
            if initial is not None:
                row[:] = initial
            else:
                row[:] = 0.0
            return
        for name, shape, size, off in zip(self.names, self.shapes,
                                          self.sizes, self.offsets):
            try:
                value = state[name]
            except KeyError as exc:
                raise ValueError(
                    f"STATE_SCHEMA: el estado devuelto por el material no "
                    f"contiene la variable {name!r} declarada en el esquema "
                    f"(claves recibidas: {sorted(state)})."
                ) from exc
            if shape == ():
                row[off] = float(value)
            else:
                arr = np.asarray(value, dtype=np.float64).reshape(-1)
                if arr.size != size:
                    raise ValueError(
                        f"STATE_SCHEMA: la variable {name!r} declara forma "
                        f"{shape} ({size} valores) pero el material devolvió "
                        f"{arr.size} valores."
                    )
                row[off:off + size] = arr

    def validate(self, state: Mapping | None) -> None:
        """Comprueba que un estado devuelto por ``compute_state`` coincide
        con el esquema en claves y formas. Lanza ``ValueError`` si no."""
        if state is None:
            if self.n_state != 0:
                raise ValueError(
                    "STATE_SCHEMA declara variables internas pero el material "
                    "devolvió None como estado."
                )
            return
        keys = set(state)
        declared = set(self.names)
        if keys != declared:
            raise ValueError(
                f"STATE_SCHEMA declara {sorted(declared)} pero el material "
                f"devolvió {sorted(keys)}."
            )
        for name, shape in zip(self.names, self.shapes):
            arr = np.asarray(state[name], dtype=np.float64)
            if shape == ():
                if arr.shape not in ((), (1,)):
                    raise ValueError(
                        f"STATE_SCHEMA: {name!r} se declara escalar pero el "
                        f"material devolvió forma {arr.shape}."
                    )
            elif arr.shape != shape:
                raise ValueError(
                    f"STATE_SCHEMA: {name!r} se declara con forma {shape} pero "
                    f"el material devolvió forma {arr.shape}."
                )
