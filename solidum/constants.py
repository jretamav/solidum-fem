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
# Factor adimensional sobre la longitud de arco del primer paso (Δl₁, la
# declarada con `initial_dl` o la derivada de `initial_dlambda`): si la
# bisección del paso reduce ``dl`` por debajo de
# ``ARCLENGTH_MIN_DL_FACTOR · Δl₁``, ArcLengthSolver aborta con
# RuntimeError. Por debajo de ese umbral la convergencia ya no es
# físicamente significativa (ruido numérico domina la longitud de arco).
ARCLENGTH_MIN_DL_FACTOR = 1.0e-6

# --- Primer paso del arc-length por omisión (2026-09-23) ---
# Fracción de la carga de referencia que lleva el primer paso en régimen
# elástico: Δl₁ = Δλ₁·‖K⁻¹·F_ref‖ (Crisfield 1991, cap. 9). Adimensional, así
# que vale igual en cualquier sistema de unidades. Sustituye al antiguo
# `initial_dl = 0.1` por omisión, una longitud que en el ejemplo 4 del manual
# de ejemplos pedía 86 veces la carga de colapso en el primer paso.
ARCLENGTH_DEFAULT_INITIAL_DLAMBDA = 0.1

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
# Umbral relativo para considerar diagonal "estrictamente" la masa reducida
# TᵀMT en ``CentralDifferenceSolver._invert_diagonal_mass``: se compara
# max|fuera de la diagonal| con este factor por max|diagonal| (relativo puro,
# independiente de las unidades desde el 2026-09-30). Por encima, se rechaza:
# Frame3D con eje oblicuo, o una restricción multipunto con varios maestros
# cuyo esclavo tiene masa.
LUMPED_MASS_OFF_DIAGONAL_RTOL = 1.0e-12

# --- Aviso de condiciones iniciales incoherentes con una MPC (ADR 0004) ---
# Tolerancia relativa (sobre el mayor valor absoluto del vector de condiciones
# iniciales o de su reconstrucción) por encima de la cual un DOF esclavo de una
# MPC con maestros se considera incoherente con ella en los solvers
# transitorios. Holgada frente al redondeo de coeficientes tecleados por el
# usuario (tan30°·u_x). Sólo gobierna un aviso de log, no un resultado.
INITIAL_CONDITION_MPC_RTOL = 1.0e-8

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
# relaja con el residuo. Validado por el usuario el 2026-09-23 (ADR 0018
# §4): fijo, deliberadamente no configurable desde el YAML.
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
# Diagnóstico del corrector de Newton (no rechaza nada): si un paso fracasa y
# alguna resolución dejó un residuo relativo ‖b − A·x‖/‖b‖ por encima de esto,
# la tangente era singular o casi singular (un sistema singular cargado no se
# puede satisfacer), y así se informa. Hasta el 2026-10-20 era también el
# criterio de la capa 2 del estático lineal; ver BACKWARD_ERROR_RTOL.
EQUILIBRIUM_RTOL = 1.0e-8
# Capa 2 del estático lineal (decisión del usuario, 2026-10-20): error hacia
# atrás normado de la solución, escalado por la diagonal,
# η = ‖r̃‖∞ / (‖K̃‖∞·‖ũ‖∞ + ‖F̃‖∞) con K̃ = D^-½·K·D^-½, ũ = D^½·u,
# F̃ = D^-½·F, r̃ = D^-½·(F − K·u), D = diag(K). Mide si el solver resolvió el
# sistema que se le dio, sin depender de las unidades ni del contraste de
# rigideces. Medido: ≤ 3.6e-16 en la suite y ≤ 5e-12 con zonas rígidas ×1e10
# en mm, donde el residuo relativo ‖F − K·u‖/‖F‖ llegaba a 1.2e-2 (y a 5.7e-9
# en la propia suite). Los mecanismos internos, cargados o no, los rechaza
# antes la capa 3 (MECHANISM_ENERGY_RTOL).
BACKWARD_ERROR_RTOL = 1.0e-10
# Umbral relativo de rango para decidir si los apoyos restringen un modo de
# cuerpo rígido: valores singulares de la matriz "modo × restricción"
# (columnas normalizadas a máximo 1) por debajo de este múltiplo del mayor
# se consideran nulos ⇒ movimiento rígido libre ⇒ mecanismo.
MECHANISM_RANK_RTOL = 1.0e-8
# Filtro de la capa 3 (ampliación del 2026-10-20): un pivote es sospechoso si
# |pivote|/√(|K_rr|·|K_cc|) < SUSPECT_PIVOT_RTOL (relativo a la diagonal, así
# que no depende de las unidades). No decide nada: con algún pivote
# sospechoso se busca el modo de menor energía y decide su energía
# (MECHANISM_ENERGY_RTOL). Medido: el pivote de un mecanismo real es ruido de
# redondeo que crece con el tamaño del modelo (1e-16 con 100 DOF, 2.4e-11 con
# 264 000); los modelos bien planteados de la suite quedan en ≥ 3.5e-8. Un
# filtro holgado sólo cuesta la comprobación en los modelos con contrastes
# de rigidez extremos, que la pasan.
SUSPECT_PIVOT_RTOL = 1.0e-8
# Un modo es un mecanismo si todos los elementos que mueve quedan sin
# deformar: energía relativa v_eᵀ·K_e·v_e / v_eᵀ·D_e·v_e < MECHANISM_ENERGY_RTOL
# en cada uno (D_e: diagonal de la K ensamblada en los DOF del elemento).
# Decisión del usuario (2026-10-20). Medido: ≤ 8e-17 en mecanismos de 17 a
# 264 000 DOF; ≥ 1e-7 con ν = 0.4999999 y 6.8e-2 en zonas rígidas de
# cualquier contraste. El umbral deja cuatro órdenes sobre el redondeo.
MECHANISM_ENERGY_RTOL = 1.0e-12
# Capa 4 del estático lineal (decisión del usuario, 2026-10-27): aviso de
# pérdida de cifras si el número de condición estimado de la matriz escalada
# por la diagonal, κ̃ = ‖K̃‖₁·‖K̃⁻¹‖₁, supera este valor; log10(κ̃) acota las
# cifras perdidas de las ~16 de la doble precisión. Medido: el error real
# queda 30 a 100 veces por debajo de la cota κ̃·ε (voladizo de 1000 elementos:
# κ̃ = 1e13, error 3e-5; zonas rígidas ×1e8: 1.3e12, 1e-5), así que con este
# umbral se avisa cuando la solución baja a unas 5 cifras. El modelo legítimo
# más mal condicionado de la suite (viga de MacNeal con Hex20) da 1.2e10.
CONDITION_WARNING_THRESHOLD = 1.0e12
