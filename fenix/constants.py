# fenix_fem/fenix/constants.py

"""
Constantes y tolerancias numéricas globales del motor Fenix FEM.
Centralizar estos valores evita "magic numbers" hardcodeados y 
permite un ajuste fino global para problemas de convergencia.
"""

ZERO_JACOBIAN_TOL = 1e-10  # Tolerancia para detectar un Jacobiano negativo o nulo
PLASTIC_YIELD_TOL = 1e-9   # Tolerancia para la evaluación de la función de fluencia (f <= tol)
DAMAGE_MAX = 0.999         # Límite máximo para el daño escalar (evita singularidad en la rigidez)
CONVERGENCE_TOL = 1e-5     # Tolerancia por defecto para la convergencia en solvers
ZERO_TOL = 1e-12           # Valor utilizado para evitar divisiones por cero