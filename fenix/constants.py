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

# --- Heurísticos del paso adaptativo del Newton-Raphson incremental ---
# Aplican en ``NonlinearSolver.solve`` con ``adaptive=True`` (default).
# El paso ``Δλ`` se acelera por ``GROWTH_FACTOR`` cuando el Newton del paso
# previo convergió en pocas iteraciones (rápido ⇒ subir el paso); se biseca
# (``Δλ ← Δλ/2``) cuando el paso diverge. ``MIN_DELTA_LAMBDA`` es el piso por
# debajo del cual el solver aborta (``LoadExceedsCapacityError``).
# El epsilon de cierre evita un paso espurio de "cola" cuando ``load_factor``
# está dentro del último ULP de ``target_load``.
NEWTON_ADAPTIVE_GROWTH_FACTOR = 1.5
NEWTON_ADAPTIVE_GROWTH_ITER_THRESHOLD = 4
NEWTON_DEFAULT_MIN_DELTA_LAMBDA = 1.0e-5
NEWTON_LOAD_FACTOR_EPSILON = 1.0e-9

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

# --- Tolerancias internas del análisis modal (ADR 0009 fase 1) ---
# Tolerancia de residuo ARPACK: criterio que ``eigsh`` usa para declarar
# converida la iteración Lanczos. Más estricto que rtol del Newton porque
# los modos se reusan en cómputos derivados (participation_factors,
# response_spectrum) que amplifican errores.
MODAL_ARPACK_TOLERANCE = 1.0e-9
# Umbral relativo para clip de autovalores espurios negativos cerca de los
# modos rígidos (``λ < MODAL_EIGENVALUE_ZERO_RTOL · λ_max`` se considera 0).
MODAL_EIGENVALUE_ZERO_RTOL = 1.0e-12

# --- Tolerancias internas del solver explícito (ADR 0009 fase 5) ---
# Umbral relativo para considerar diagonal "estrictamente" en
# ``CentralDifferenceSolver._invert_diagonal_mass``. Por encima de este,
# se rechaza el ensamble lumped (típico: Frame3D oblicuo).
LUMPED_MASS_OFF_DIAGONAL_RTOL = 1.0e-12

# --- Newton local del salto cohesivo en CST_Embedded2D (ADR 0010) ---
# Tolerancia relativa del Newton local sobre ``[[u]]`` (residuo r_jump).
# Más estricta que la del Newton global porque típicamente converge en
# 1-3 iteraciones gracias al penalty stiff K_e ≈ 1e15.
EMBEDDED_LOCAL_JUMP_RTOL = 1.0e-10
EMBEDDED_LOCAL_JUMP_MAX_ITER = 30

# --- Return mapping J2 plane stress (Simó-Hughes §3.4.1) ---
# Newton local sobre ``Δγ`` con función de fluencia proyectada
# ``f̄ = ½·σ·P·σ − R²/3``. Converge típicamente en 3-6 iteraciones por
# la tangente cerrada; cota superior conservadora.
J2_PLANE_STRESS_MAX_LOCAL_ITER = 25
# Piso para denominadores en la corrección de tangente del J2 plane stress.
# Por debajo de este, el material reporta como "casi degenerado" en la
# rama Newton local (situación patológica que tests no han forzado).
J2_DENOM_FLOOR = 1.0e-30