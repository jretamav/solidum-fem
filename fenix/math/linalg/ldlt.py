"""Backend LDLᵀ para matrices simétricas indefinidas (placeholder fase 2).

**Estado**: no implementado. Se activa cuando un caso real lo requiera —
típicamente al cruzar un punto de bifurcación (pandeo, snap-through) donde
``K_t`` pierde positividad pero conserva simetría.

**Por qué un placeholder y no una implementación basada en LU**: extraer
``sign(diag(U))`` de un LU general no equivale al número de autovalores
negativos. SuperLU usa pivoteo parcial por filas que no preserva la inercia
de la matriz; sólo una factorización LDLᵀ verdadera (con pivoteo simétrico
1×1/2×2 de Bunch-Kaufman) cumple la ley de inercia de Sylvester. Inventar
un Sturm sequence aproximado que diera bandera roja en bifurcaciones
incorrectas sería peor que no tenerlo.

**Backends candidatos para activarlo**:

- ``pypardiso`` (Intel MKL Pardiso): industrial, soporta LDLᵀ con conteo de
  pivotes negativos nativo. Free para uso académico (la licencia es de la
  MKL Pardiso completa, no de pypardiso). Dependencia pesada en disco pero
  el binario MKL ya viene con muchos numpy/scipy por defecto.
- ``python-mumps`` (MUMPS): equivalente open-source, multifrontal, soporta
  LDLᵀ y Sturm count. Más difícil de instalar en Windows.
- ``scipy.linalg.ldl``: LDLᵀ denso. Inviable para el rango de tallas de
  Fenix; mencionado solo para descartarlo.

Cuando llegue el caso, sustituir el cuerpo de :class:`LDLTSolver` por el
envoltorio adecuado y exponer ``n_negative_pivots`` en :class:`LDLTFactorized`.
"""
from __future__ import annotations

import warnings

import numpy as np
import scipy.sparse as sp


class LDLTNotAvailableError(RuntimeError):
    """LDLᵀ no está implementado; el despachador degrada a LU automáticamente."""


class LDLTSolver:
    """Placeholder de LDLᵀ (no implementado). El despachador degrada a LU."""

    name = "ldlt"

    _warning_emitted = False

    @classmethod
    def _warn_once(cls) -> None:
        if cls._warning_emitted:
            return
        cls._warning_emitted = True
        warnings.warn(
            "LDLTSolver no implementado en fase 2 del ADR 0003. Degradando a LU. "
            "Para activar Sturm sequence (detección de bifurcaciones), instalar "
            "pypardiso y abrir ADR para enchufar el backend.",
            RuntimeWarning,
            stacklevel=3,
        )

    def solve(self, K: sp.spmatrix, b: np.ndarray) -> np.ndarray:  # pragma: no cover
        raise LDLTNotAvailableError(
            "LDLTSolver no implementado. Use LUSolver mientras tanto."
        )

    def factorize(self, K: sp.spmatrix):  # pragma: no cover
        raise LDLTNotAvailableError(
            "LDLTSolver.factorize no implementado. Use LUSolver.factorize mientras tanto."
        )
