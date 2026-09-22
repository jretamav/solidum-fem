"""Firmas Numba de los kernels puntuales y del kernel de familia (ADR 0014).

Las firmas son el contrato que hace posible **un único kernel de familia
compilado una vez y cacheado en disco**: la cinemática del elemento y la
constitutiva del material entran como *funciones de primera clase*
tipadas (``types.FunctionType``), no como despachadores. Un despachador
tiene identidad propia por proceso y obligaría a recompilar el kernel
para cada pareja elemento × material en cada ejecución (medido: 3 s por
pareja); una función tipada por su firma es un tipo estructural, así que
el kernel se compila una sola vez para todas las parejas y el caché de
Numba lo reutiliza entre procesos.

Cinemática — ``KIN_SIG``
------------------------
``detJ = kin(pt, coords, B)``

- ``pt``: coordenadas naturales del punto de Gauss, ``(d,)``.
- ``coords``: coordenadas de referencia del elemento, ``(n_nodos, d)``.
- ``B``: matriz deformación–desplazamiento, ``(n_sigma, n_dof)``, que la
  función **rellena por completo** (incluidos los ceros).
- Devuelve ``det J``. **No lanza**: si el jacobiano degenera (mismo
  chequeo relativo que el camino por elemento) devuelve un valor ``≤ 0``
  y deja ``B`` sin escribir. El kernel de familia marca el punto y el
  ensamblador lanza el ``ValueError`` con el id del elemento. La razón
  es que una excepción dentro de un bucle ``prange`` de Numba se pierde,
  y el mismo contrato debe servir a las variantes serie y paralela.

Constitutiva — ``MAT_SIG``
--------------------------
``C_t = mat(strain, S_old, S_new, params, C, sigma, C_out, flag)``

- ``strain``: deformación en Voigt del proyecto, ``(n_sigma,)``.
- ``S_old`` / ``S_new``: fila de estado committed / trial, ``(n_state,)``,
  con el orden que fija ``Material.STATE_SCHEMA``. El kernel **escribe la
  fila trial completa**.
- ``params``: parámetros escalares del material, ``(n_par,)``, en el
  orden que documenta cada ``Material.batch_params``.
- ``C``: matriz del material (elástica o precalculada), ``(n_sigma, n_sigma)``.
- ``sigma``: salida, esfuerzo ``(n_sigma,)``, escrita in situ.
- ``C_out``: espacio de trabajo ``(n_sigma, n_sigma)`` donde el kernel
  puede escribir la tangente sin asignar memoria por punto de Gauss.
- ``flag``: vista de un entero (``(1,)``, ``int8``); el kernel escribe
  ``1`` para señalar una incidencia que el material reporta después con
  ``batch_report`` (p. ej. Newton local sin converger).
- Devuelve la tangente algorítmica ``(n_sigma, n_sigma)``, contigua:
  ``C`` (lineal), ``C_out`` (escrita in situ) o un arreglo propio.
"""
from numba import types

F64_1 = types.float64[::1]
F64_2 = types.float64[:, ::1]
F64_3 = types.float64[:, :, ::1]
I8_1 = types.int8[::1]
I32_1 = types.int32[::1]
I64_1 = types.int64[::1]

# reduce(data, src_ptr, src_idx, out): reducción COO → CSR por mapa inverso.
REDUCE_SIG = types.void(F64_1, I64_1, I32_1, F64_1)

KIN_SIG = types.float64(F64_1, F64_2, F64_2)
MAT_SIG = F64_2(F64_1, F64_1, F64_1, F64_1, F64_2, F64_1, F64_2, I8_1)

# gauss(kin, mat, X, u, pts, S_in, S_scratch, params, C, eps_out, sig_out, flags):
# post-proceso por familia — deformación y esfuerzo por punto de Gauss a
# partir del estado committed, sin tocar el trial (S_scratch es un espacio
# de trabajo con la forma de S_trial que se descarta).
GAUSS_SIG = types.void(
    types.FunctionType(KIN_SIG),
    types.FunctionType(MAT_SIG),
    F64_3,   # X        (N, n_nodos, d)
    F64_2,   # u        (N, n_dof)
    F64_2,   # pts      (n_gp, d)
    F64_2,   # S_in     (N·n_gp, n_state)
    F64_2,   # S_scratch(N·n_gp, n_state)
    F64_1,   # params   (n_par,)
    F64_2,   # C        (n_sigma, n_sigma)
    F64_2,   # eps_out  (N·n_gp, n_sigma)
    F64_2,   # sig_out  (N·n_gp, n_sigma)
    I8_1,    # flags    (N·n_gp,)
)

# kernel(kin, mat, X, u, pts, w, scale, S_in, S_out, params, C,
#        K_out, F_out, sig_out, flags)
FAMILY_SIG = types.void(
    types.FunctionType(KIN_SIG),
    types.FunctionType(MAT_SIG),
    F64_3,   # X      (N, n_nodos, d)
    F64_2,   # u      (N, n_dof)
    F64_2,   # pts    (n_gp, d)
    F64_1,   # w      (n_gp,)
    F64_1,   # scale  (N,)
    F64_2,   # S_in   (N·n_gp, n_state)
    F64_2,   # S_out  (N·n_gp, n_state)
    F64_1,   # params (n_par,)
    F64_2,   # C      (n_sigma, n_sigma)
    F64_3,   # K_out  (N, n_dof, n_dof)
    F64_2,   # F_out  (N, n_dof)
    F64_2,   # sig_out (N·n_gp, n_sigma)
    I8_1,    # flags  (N·n_gp,)
)
