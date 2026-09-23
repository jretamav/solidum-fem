"""Despachador de solver algebraico (ADR 0003 §3).

Selecciona el backend adecuado a partir de las propiedades declarativas de
``K``. Si el usuario fuerza un backend desde el YAML (``override``), se
respeta su elección — herramienta de diagnóstico, no decisión de modelado.
"""
from __future__ import annotations

import warnings

from solidum.math.linalg.base import LinearAlgebraSolver, StiffnessProperties
from solidum.math.linalg.iterative import (
    HAS_PYAMG as _HAS_PYAMG,
    PRECONDITIONERS as _ITERATIVE_PRECONDITIONERS,
    IterativeSolver,
)
from solidum.math.linalg.ldlt import LDLTSolver
from solidum.math.linalg.lu import LUSolver

# CholeskySolver es opcional: depende de ``scikit-sparse``.
try:
    from solidum.math.linalg.cholesky import CholeskySolver
    _HAS_CHOLESKY = True
except ImportError:
    _HAS_CHOLESKY = False
    CholeskySolver = None  # type: ignore[assignment,misc]

# PardisoSolver es opcional: depende de ``pypardiso`` (Intel MKL).
try:
    from solidum.math.linalg.pardiso import PardisoSolver
    _HAS_PARDISO = True
except ImportError:
    _HAS_PARDISO = False
    PardisoSolver = None  # type: ignore[assignment,misc]


# El backend iterativo (ADR 0018) no tiene dependencia obligatoria —sin
# ``pyamg`` funciona sin precondicionador AMG—, así que siempre se registra.
_REGISTRY: dict[str, type] = {"lu": LUSolver, "ldlt": LDLTSolver,
                              "iterative": IterativeSolver}
if _HAS_CHOLESKY:
    _REGISTRY["cholesky"] = CholeskySolver  # type: ignore[assignment]
if _HAS_PARDISO:
    _REGISTRY["pardiso"] = PardisoSolver  # type: ignore[assignment]


def _parse_override(override: str) -> tuple[str, dict]:
    """``'iterative:amg'`` → ``('iterative', {'preconditioner': 'amg'})``.

    Sólo el backend iterativo admite opción (su precondicionador); el resto
    de nombres se aceptan sin sufijo. Lanza ``ValueError`` con la lista de
    valores válidos ante cualquier otra forma."""
    key, sep, option = override.lower().partition(":")
    if key not in _REGISTRY:
        raise ValueError(
            f"Backend algebraico desconocido: '{override}'. "
            f"Disponibles: {available_overrides()}."
        )
    if not sep:
        return key, {}
    if key != "iterative" or option not in _ITERATIVE_PRECONDITIONERS:
        raise ValueError(
            f"Backend algebraico desconocido: '{override}'. Sólo 'iterative' "
            f"admite sufijo, y el precondicionador debe ser uno de "
            f"{list(_ITERATIVE_PRECONDITIONERS)}. Disponibles: "
            f"{available_overrides()}."
        )
    if option == "amg" and not _HAS_PYAMG:
        raise ValueError(
            "linear_algebra 'iterative:amg' requiere pyamg, que no está "
            "instalado (`pip install solidum-fem[iterative]`). Alternativas: "
            "'iterative:none', 'iterative:jacobi' o 'auto'."
        )
    return key, {"preconditioner": option}


def available_overrides() -> list[str]:
    """Valores válidos de ``linear_algebra`` en este entorno."""
    names = ["auto"] + sorted(_REGISTRY)
    names += [f"iterative:{p}" for p in _ITERATIVE_PRECONDITIONERS
              if p != "auto" and (p != "amg" or _HAS_PYAMG)]
    return names


def is_valid_override(override: str | None) -> bool:
    """``True`` si ``override`` es un valor aceptable de ``linear_algebra``
    (lo usa el parser YAML para fallar al leer, no al resolver)."""
    if override is None or override == "auto":
        return True
    try:
        _parse_override(str(override))
    except ValueError:
        return False
    return True


# LDLᵀ todavía no está implementado en fase 2 (decisión documentada en
# solidum/math/linalg/ldlt.py). El despachador degrada a LU con warning.
_LDLT_AVAILABLE = False


def select_solver(
    props: StiffnessProperties,
    override: str | None = None,
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
        key, options = _parse_override(override)
        if key == "iterative":
            # CG y MINRES exigen simetría. Una tangente no simétrica
            # (plasticidad no asociada, cargas seguidoras) va a un solver
            # directo, que es lo que hacen ANSYS y Abaqus con sus solvers
            # iterativos (ADR 0018 §2).
            if not props.is_symmetric:
                _warn_iterative_nonsymmetric_once()
                return select_solver(props, override=None)
            return IterativeSolver(props, **options)
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
    # está obteniendo el speedup esperado. Con Pardiso instalado el aviso
    # sobra: hay un backend mejor que LU disponible justo debajo.
    if (props.is_symmetric and props.is_positive_definite
            and not _HAS_CHOLESKY and not _HAS_PARDISO):
        _warn_cholesky_unavailable_once()

    # Pardiso ocupa el lugar de LU: mismo dominio de aplicación (real, sin
    # hipótesis de simetría ni positividad) y misma solución, pero multihilo
    # y con disección anidada. La ventaja crece con la talla — de ×7 a 4 000
    # DOF a ×60 a 53 000 (ADR 0017) — porque ataca el relleno, que es lo que
    # degrada el escalado de un solver directo.
    if _HAS_PARDISO:
        return _REGISTRY["pardiso"]()

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


_iterative_nonsym_warned = False


def _warn_iterative_nonsymmetric_once() -> None:
    global _iterative_nonsym_warned
    if _iterative_nonsym_warned:
        return
    _iterative_nonsym_warned = True
    warnings.warn(
        "linear_algebra 'iterative' pedido sobre una matriz no simétrica "
        "(plasticidad no asociada, cargas seguidoras...). CG y MINRES exigen "
        "simetría: se usa el solver directo automático para este sistema.",
        RuntimeWarning,
        stacklevel=3,
    )
