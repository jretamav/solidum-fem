"""Verificación de Solidum por soluciones manufacturadas (ADR 0022).

No lo carga ``import solidum`` ni el descubrimiento automático: necesita
SymPy, que no es dependencia del programa principal
(``pip install solidum-fem[verification]``).
"""
from solidum.verification.mms import (
    ConvergenceStudy,
    ManufacturedElasticity,
    convergence_study,
    error_norms,
    interpolation_degree,
    linear_elastic_matrix,
    solve_manufactured,
)

__all__ = [
    "ConvergenceStudy",
    "ManufacturedElasticity",
    "convergence_study",
    "error_norms",
    "interpolation_degree",
    "linear_elastic_matrix",
    "solve_manufactured",
]
