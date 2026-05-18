# fenix_fem/fenix/constants.py

"""
Constantes y tolerancias numéricas globales del motor Fenix FEM.
Centralizar estos valores evita "magic numbers" hardcodeados y 
permite un ajuste fino global para problemas de convergencia.
"""

ZERO_JACOBIAN_TOL = 1e-10  # Tolerancia para detectar un Jacobiano negativo o nulo
DAMAGE_MAX = 0.999         # Límite máximo para el daño escalar (evita singularidad en la rigidez)
ZERO_TOL = 1e-12           # Pequeñas operaciones numéricas (normalizaciones, comparaciones de cierre de paso)

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

# --- Tolerancias de convergencia en solvers no lineales (ADR 0007) ---
# Patrón atol + rtol · escala separado por criterio (fuerza y desplazamiento):
#     ‖R‖   ≤ atol_force + rtol_force · max(‖F_ext‖, ‖F_int‖)
#     ‖δU‖  ≤ atol_disp  + rtol_disp  · ‖U_iter‖
# Los atol efectivos se autoderivan de la escala del problema durante
# ConvergenceCriterion.calibrate(), de modo que estas constantes son todas
# adimensionales y el código es invariante bajo cambio de unidades.
CONVERGENCE_RTOL_FORCE = 1.0e-5
CONVERGENCE_RTOL_DISP = 1.0e-5
CONVERGENCE_ATOL_FORCE_FACTOR = 1.0e-9
CONVERGENCE_ATOL_DISP_FACTOR = 1.0e-9

# --- Tolerancia mínima de longitud de arco antes de declarar fracaso ---
# Factor adimensional sobre el `initial_dl` declarado por el usuario:
# si la bisección del paso reduce ``dl`` por debajo de
# ``ARCLENGTH_MIN_DL_FACTOR · initial_dl``, ArcLengthSolver aborta con
# RuntimeError. Por debajo de ese umbral la convergencia ya no es
# físicamente significativa (ruido numérico domina la longitud de arco).
ARCLENGTH_MIN_DL_FACTOR = 1.0e-6

# --- Line search Armijo en solvers no lineales (ADR 0011) ---
# Globalización del Newton-Raphson: tras computar el incremento δU, se busca
# un factor α ∈ (0, 1] que satisfaga la condición de Armijo (suficiente
# decrecimiento del residuo):
#     ‖R(U + α·δU)‖ ≤ (1 − c₁·α) · ‖R(U)‖
# Si α=1 no satisface, se hace backtracking α ← ρ·α hasta MAX_BACKTRACKS retrocesos.
# Valores canónicos de Nocedal-Wright §3.1.
LINE_SEARCH_C1 = 1.0e-4              # Constante de Armijo (suficiente decrecimiento)
LINE_SEARCH_RHO = 0.5                # Factor de backtracking
LINE_SEARCH_MAX_BACKTRACKS = 10      # Cota superior de retrocesos por iteración