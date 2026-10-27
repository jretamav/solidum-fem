"""Verificación de Solidum por soluciones manufacturadas (ADR 0022).

No lo carga ``import solidum`` ni el descubrimiento automático: necesita
SymPy, que no es dependencia del programa principal
(``pip install solidum-fem[verification]``).
"""
from solidum.verification.mms import (
    ConvergenceStudy,
    ManufacturedConduction,
    ManufacturedElasticity,
    ManufacturedFrame,
    convergence_study,
    error_norms,
    interpolation_degree,
    linear_elastic_matrix,
    solve_manufactured,
    solve_manufactured_transient,
    time_convergence_study,
)

__all__ = [
    "ConvergenceStudy",
    "ManufacturedConduction",
    "ManufacturedElasticity",
    "ManufacturedFrame",
    "convergence_study",
    "error_norms",
    "interpolation_degree",
    "linear_elastic_matrix",
    "solve_manufactured",
    "solve_manufactured_transient",
    "time_convergence_study",
]
