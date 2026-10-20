"""Interfaz común y propiedades declarativas de ``K`` (ADR 0003)."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Protocol, runtime_checkable

import numpy as np
import scipy.sparse as sp


@dataclass(frozen=True)
class StiffnessProperties:
    """Propiedades declarativas de la matriz de rigidez global.

    Las usa ``select_solver`` para escoger el backend algebraico adecuado.
    Todos los flags son derivables del modelo (materiales, elementos, tipo de
    análisis) y no se piden al usuario.

    ``near_nullspace`` es un **proveedor perezoso** del casi-núcleo del
    sistema reducido (los modos de cuerpo rígido restringidos a los DOF
    libres, ``Assembler.near_nullspace``). Sólo lo consulta el precondicionador
    AMG del backend iterativo (ADR 0018); el resto de backends lo ignoran y,
    al ser perezoso, no pagan su cálculo. No interviene en la igualdad.
    """

    is_symmetric: bool
    is_positive_definite: bool
    size: int
    near_nullspace: Callable[[], np.ndarray] | None = field(
        default=None, compare=False, repr=False,
    )


@runtime_checkable
class FactorizedSolver(Protocol):
    """Factorización numérica retenida; ``solve(b)`` puede llamarse N veces.

    Habilita Newton modificado (factorización congelada durante varias
    iteraciones), dinámica implícita con paso fijo, y eigensolvers tipo
    shift-invert. La factorización es cara; reusarla amortiza el coste.

    Atributos opcionales:

    - ``n_negative_pivots`` (``int | None``): número de pivotes negativos en
      la factorización. Para LDLᵀ de matriz simétrica, equivale al número de
      autovalores negativos por la ley de inercia de Sylvester — útil para
      detectar bifurcaciones (pandeo, snap-through). ``None`` si el backend
      no lo expone.
    - ``n_small_pivots`` (``int | None``): pivotes pequeños frente a la
      diagonal, ``< SUSPECT_PIVOT_RTOL`` (ver :func:`pivot_ratios`). Es el
      filtro barato de la red de seguridad del ADR 0019 (ampliación del
      2026-10-20): con alguno, el estático lineal y la tangente inicial de
      los no lineales buscan un mecanismo interno por la energía del modo de
      menor energía, que es lo que decide. Un pivote pequeño no basta para
      declarar un mecanismo: el de un mecanismo real es ruido de redondeo que
      crece con el tamaño del modelo, y un contraste de rigidez legítimo
      (zonas rígidas, ν cerca de 1/2) también da pivotes pequeños. Lo
      exponen SuperLU, Pardiso y Cholesky; ``None`` si el backend no tiene
      pivotes que inspeccionar (el iterativo).
    """

    n_negative_pivots: int | None

    def solve(self, b: np.ndarray) -> np.ndarray:
        ...


@runtime_checkable
class LinearAlgebraSolver(Protocol):
    """Contrato uniforme para los backends que resuelven ``K·x = b``.

    Cada backend (LU, Cholesky, ...) lo implementa con su propia mecánica.
    """

    name: str

    def solve(self, K: sp.spmatrix, b: np.ndarray) -> np.ndarray:
        """Resuelve ``K·x = b`` en una sola llamada.

        Equivalente a ``factorize(K).solve(b)`` pero sin retener la factorización.
        """
        ...

    def factorize(self, K: sp.spmatrix) -> FactorizedSolver:
        """Factoriza ``K`` y retorna un objeto reutilizable para ``solve(b)``.

        Separar la fase cara (factorización) de la barata (sustitución
        triangular) habilita Newton modificado y reuso entre pasos.
        """
        ...


def pivot_ratios(pivots: np.ndarray, diag: np.ndarray, rows: np.ndarray,
                 cols: np.ndarray | None = None) -> np.ndarray:
    """Cociente de cada pivote con la diagonal de su fila y su columna,
    ``|p_i| / √(|K_rr|·|K_cc|)`` con ``r = rows[i]``, ``c = cols[i]`` (en una
    factorización simétrica, ``cols = rows`` y el cociente es ``|d_i|/K_rr``).

    Un cambio de unidades de un tipo de grado de libertad multiplica su fila
    y su columna por el mismo factor, y el cociente no cambia: es el pivote
    de la matriz escalada por su diagonal. Un pivote cuya fila o columna
    tiene diagonal nula da cociente 0 (sospechoso), porque esa matriz no se
    puede escalar."""
    cols = rows if cols is None else cols
    d = np.abs(np.asarray(diag, dtype=float))
    ref = np.sqrt(d[rows] * d[cols])
    p = np.abs(np.asarray(pivots, dtype=float))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(ref > 0.0, p / np.where(ref > 0.0, ref, 1.0), 0.0)


def out_of_memory_error(n: int, backend: str) -> MemoryError:
    """Error de memoria del solver directo con la salida para el usuario
    (ADR 0019). La memoria de un directo crece más deprisa que el modelo
    (relleno de los factores: ×14 sobre ``K`` a 26k DOF, ×25 a 200k); la de
    un iterativo es un múltiplo fijo de ``K`` (ADR 0018)."""
    return MemoryError(
        f"Memoria insuficiente para factorizar el sistema de {n} incógnitas "
        f"con el solver directo ({backend}). La memoria de un solver directo "
        f"crece más deprisa que el tamaño del modelo. Para modelos grandes "
        f"usa `linear_algebra: iterative` en el bloque `solver` (solver "
        f"iterativo con AMG, memoria proporcional al modelo; ver el anexo de "
        f"la capa algebraica del manual)."
    )

