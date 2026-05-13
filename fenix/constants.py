# fenix_fem/fenix/constants.py

"""
Constantes y tolerancias numéricas globales del motor Fenix FEM.
Centralizar estos valores evita "magic numbers" hardcodeados y 
permite un ajuste fino global para problemas de convergencia.
"""

ZERO_JACOBIAN_TOL = 1e-10  # Tolerancia para detectar un Jacobiano negativo o nulo
DAMAGE_MAX = 0.999         # Límite máximo para el daño escalar (evita singularidad en la rigidez)
CONVERGENCE_TOL = 1e-5     # Tolerancia por defecto para la convergencia en solvers
ZERO_TOL = 1e-12           # Valor utilizado para evitar divisiones por cero

# --- Tolerancias del criterio de admisibilidad constitutiva (ADR 0006) ---
# Comparación combinada estilo ODE solver:
#     f <= ADMISSIBILITY_TOL_ABS + ADMISSIBILITY_TOL_REL · escala(estado)
# `escala` la declara cada material vía Material.admissibility_scale(state),
# en las mismas unidades que la cantidad `f` a comparar (esfuerzo).
# - REL=1e-10 ⇒ banda de seis órdenes sobre ε_machine (~1e-16) en doble.
# - ABS=1e-14 ⇒ piso cuando escala(estado) → 0 (estados degenerados o
#   modelos sin parámetro de referencia inicial). Por encima del ruido.
ADMISSIBILITY_TOL_REL = 1.0e-10
ADMISSIBILITY_TOL_ABS = 1.0e-14