"""Interfaz común y propiedades declarativas de ``K`` (ADR 0003)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np
import scipy.sparse as sp


@dataclass(frozen=True)
class StiffnessProperties:
    """Propiedades declarativas de la matriz de rigidez global.

    Las usa ``select_solver`` para escoger el backend algebraico adecuado.
    Todos los flags son derivables del modelo (materiales, elementos, tipo de
    análisis) y no se piden al usuario.
    """

    is_symmetric: bool
    is_positive_definite: bool
    size: int


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
