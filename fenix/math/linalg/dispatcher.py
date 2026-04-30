"""Despachador de solver algebraico (ADR 0003 §3).

Selecciona el backend adecuado a partir de las propiedades declarativas de
``K``. Si el usuario fuerza un backend desde el YAML (``override``), se
respeta su elección — herramienta de diagnóstico, no decisión de modelado.
"""
from __future__ import annotations

import warnings
from typing import Optional

from fenix.math.linalg.base import LinearAlgebraSolver, StiffnessProperties
from fenix.math.linalg.ldlt import LDLTSolver
from fenix.math.linalg.lu import LUSolver

# CholeskySolver es opcional: depende de ``scikit-sparse``.
try:
    from fenix.math.linalg.cholesky import CholeskySolver
    _HAS_CHOLESKY = True
except ImportError:
    _HAS_CHOLESKY = False
    CholeskySolver = None  # type: ignore[assignment,misc]


_REGISTRY: dict[str, type] = {"lu": LUSolver, "ldlt": LDLTSolver}
if _HAS_CHOLESKY:
    _REGISTRY["cholesky"] = CholeskySolver  # type: ignore[assignment]


# LDLᵀ todavía no está implementado en fase 2 (decisión documentada en
# fenix/math/linalg/ldlt.py). El despachador degrada a LU con warning.
_LDLT_AVAILABLE = False


def select_solver(
    props: StiffnessProperties,
    override: Optional[str] = None,
) -> LinearAlgebraSolver:
    """Devuelve la instancia de backend algebraico adecuada para ``K``.

    Parameters
    ----------
    props
        Propiedades declarativas de ``K`` (simetría, positividad, talla).
    override
        Nombre del backend forzado (``'lu'``, ``'cholesky'``, ...). Si es
        ``None`` o ``'auto'``, se aplica la regla automática.
    """
    if override is not None and override != "auto":
        key = override.lower()
        if key not in _REGISTRY:
            available = sorted(_REGISTRY)
            raise ValueError(
                f"Backend algebraico desconocido: '{override}'. "
                f"Disponibles: {available}."
            )
        # Override 'ldlt' con el placeholder no implementado: aviso y degrada.
        if key == "ldlt" and not _LDLT_AVAILABLE:
            LDLTSolver._warn_once()
            return _REGISTRY["lu"]()
        return _REGISTRY[key]()

    # Regla automática.
    if props.is_symmetric and props.is_positive_definite and _HAS_CHOLESKY:
        return _REGISTRY["cholesky"]()

    # Si la dependencia opcional falta y el caso era ideal para Cholesky,
    # avisamos una sola vez por sesión para que el usuario sepa por qué no
    # está obteniendo el speedup esperado.
    if props.is_symmetric and props.is_positive_definite and not _HAS_CHOLESKY:
        _warn_cholesky_unavailable_once()

    # Simétrica indefinida sería ideal para LDLᵀ (Sturm sequence en pandeo y
    # snap-through). Mientras LDLᵀ no esté implementado, LU resuelve
    # correctamente — solo se pierde el conteo de pivotes negativos como
    # diagnóstico, así que NO emitimos warning automático aquí (haría ruido
    # en cada ejecución de ArcLengthSolver). El usuario que pida el
    # diagnóstico vía override ``linear_algebra: ldlt`` sí ve el aviso.

    return _REGISTRY["lu"]()


_cholesky_warning_emitted = False


def _warn_cholesky_unavailable_once() -> None:
    global _cholesky_warning_emitted
    if _cholesky_warning_emitted:
        return
    _cholesky_warning_emitted = True
    warnings.warn(
        "scikit-sparse no está instalado: el despachador degrada a LU para una "
        "matriz que sería ideal para Cholesky. Para activar Cholesky:\n"
        "    conda install -c conda-forge scikit-sparse",
        RuntimeWarning,
        stacklevel=3,
    )
