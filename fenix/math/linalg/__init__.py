"""Capa algebraica de Fenix FEM (ADR 0003).

Resolución del sistema lineal ``K·x = b`` detrás de una interfaz uniforme con
varios backends (LU, Cholesky, ...). El despachador interno selecciona el
backend adecuado a partir de propiedades declarativas de ``K``; el usuario no
elige solver algebraico desde el YAML como decisión obligatoria.
"""
from fenix.math.linalg.base import (
    FactorizedSolver,
    LinearAlgebraSolver,
    StiffnessProperties,
)
from fenix.math.linalg.dispatcher import select_solver
from fenix.math.linalg.eigen import EigenSolver
from fenix.math.linalg.ldlt import LDLTNotAvailableError, LDLTSolver
from fenix.math.linalg.lu import LUFactorized, LUSolver

# Cholesky es opcional: depende de scikit-sparse.
try:
    from fenix.math.linalg.cholesky import (  # noqa: F401
        CholeskyFactorized,
        CholeskyNotPositiveDefiniteError,
        CholeskySolver,
    )
    _HAS_CHOLESKY = True
except ImportError:
    _HAS_CHOLESKY = False

__all__ = [
    "FactorizedSolver",
    "LinearAlgebraSolver",
    "StiffnessProperties",
    "select_solver",
    "LUSolver",
    "LUFactorized",
    "LDLTSolver",
    "LDLTNotAvailableError",
    "EigenSolver",
]
if _HAS_CHOLESKY:
    __all__ += ["CholeskySolver", "CholeskyFactorized", "CholeskyNotPositiveDefiniteError"]
