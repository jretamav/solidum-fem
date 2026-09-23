# solidum_fem/solidum/constants.py

"""
Constantes y tolerancias numéricas globales del motor Solidum FEM.
Centralizar estos valores evita "magic numbers" hardcodeados y 
permite un ajuste fino global para problemas de convergencia.
"""

# Tolerancia RELATIVA para detectar un jacobiano degenerado o invertido:
#     det J <= JACOBIAN_RTOL · Π_i ‖fila_i(J)‖
# El producto de las normas de las filas es la cota de Hadamard de |det J|,
# así que el cociente det J / Π‖fila_i‖ es adimensional: vale 1 en un
# elemento con filas ortogonales y tiende a 0 al colapsar. Una tolerancia
# absoluta sobre det J (que tiene unidades de L² o L³) rechazaba mallas
# finas en metros —un Hex8 de arista < 1 mm— sin que estuvieran
# distorsionadas (auditoría 2026-09-22).
JACOBIAN_RTOL = 1e-10
ZERO_JACOBIAN_TOL = JACOBIAN_RTOL  # alias histórico (mismo valor, ahora relativo)
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
# --- Ensamblaje por lotes (ADR 0014) ---
# Activa por defecto el camino por lotes en ``Assembler`` (``batch=None``
# lee esta constante). Los componentes que no declaran kernel siguen el
# camino por elemento sin que el usuario haga nada; poner ``False`` fuerza
# el camino por elemento en todo el modelo (diagnostico, comparaciones).
BATCH_ASSEMBLY_DEFAULT = True
# Presupuesto de memoria (bytes) para los temporales de un trozo de familia.
# Las matrices elementales se escriben directamente en el vector COO de la
# topologia cacheada (sin buffer aparte), asi que el unico temporal es el
# vector F_e del trozo: n_c = presupuesto / (8·n_dof). Con 64 MB, un Hex27
# (81 DOF) procesa ~100 000 elementos por trozo; en la practica un trozo
# cubre la familia entera y el limite solo actua en mallas enormes.
BATCH_MEMORY_BUDGET_BYTES = 64 * 1024 * 1024
# Variante paralela del kernel de familia (``prange`` sobre bloques de
# elementos; ADR 0014 §9). ``Assembler(parallel=None)`` lee esta constante.
# El número de hilos lo gobierna Numba (``numba.set_num_threads`` o la
# variable de entorno ``NUMBA_NUM_THREADS``); el resultado es bit a bit el
# mismo con cualquier número de hilos y que la variante serie.
BATCH_PARALLEL_DEFAULT = True

# --- Solver algebraico iterativo (ADR 0018) ---
# Tolerancia sobre el residuo relativo VERDADERO del sistema lineal,
# ``‖b − K·x‖ ≤ ITERATIVE_RTOL · ‖b‖``, recalculado al terminar (no el
# residuo recursivo de CG, que deriva del verdadero por redondeo). Cinco
# órdenes por debajo de ``CONVERGENCE_RTOL_FORCE`` (1e-5): dentro de un
# Newton es un forzado η = 1e-10 constante, de modo que la contracción del
# residuo no lineal ‖R_{k+1}‖ ≲ η‖R_k‖ + C‖R_k‖² queda dominada por el
# término cuadrático en todo el rango útil y la sucesión de iterados es la
# del Newton exacto a efectos prácticos. No es Newton inexacto: no se
# relaja con el residuo.
ITERATIVE_RTOL = 1.0e-10
# Cota de iteraciones de Krylov por resolución. Con AMG y modos de cuerpo
# rígido bastan 10-20 (medido); sin precondicionador, del orden de cientos
# a 2·10⁵ DOF. Agotarla es un fallo explícito, nunca una solución parcial.
ITERATIVE_MAX_ITER = 10_000
# Reinicios con el residuo verdadero ("residual replacement"): si al
# converger el residuo recursivo el verdadero no cumple la tolerancia, CG
# se reinicia desde el iterado actual con el residuo recalculado.
ITERATIVE_MAX_RESTARTS = 3
# Tamaño del nivel más grueso de la jerarquía AMG, resuelto de forma
# directa. 500 incógnitas mantienen barato ese solve y acortan el setup.
AMG_MAX_COARSE = 500

# --- Red de seguridad del análisis estático (ADR 0019) ---
# Residuo de equilibrio relativo máximo aceptable tras resolver un sistema
# lineal estático: ‖F − K·u‖ ≤ EQUILIBRIUM_RTOL · ‖F‖. Medido: un sistema
# bien planteado queda en ~1e-13 con cualquier solver directo y en ≤ 1e-10
# con el iterativo; un mecanismo cargado da ~1. El umbral deja cinco órdenes
# de margen hacia cada lado.
EQUILIBRIUM_RTOL = 1.0e-8
# Umbral relativo de rango para decidir si los apoyos restringen un modo de
# cuerpo rígido: valores singulares de la matriz "modo × restricción"
# (columnas normalizadas a máximo 1) por debajo de este múltiplo del mayor
# se consideran nulos ⇒ movimiento rígido libre ⇒ mecanismo.
MECHANISM_RANK_RTOL = 1.0e-8
# Un pivote de la factorización directa se considera numéricamente nulo si
# |pivote| < ZERO_PIVOT_RTOL · max|K|. Es el mismo criterio con el que MKL
# Pardiso decide perturbar un pivote (iparm(10) = 13 por defecto), aplicado
# también a SuperLU para que ambos detecten la singularidad igual. Medido:
# modelo bien planteado 3e-3; mecanismo interno 1.6e-16.
ZERO_PIVOT_RTOL = 1.0e-13
