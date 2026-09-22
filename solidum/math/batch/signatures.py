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
- Devuelve ``det J``. Lanza ``ValueError`` si el jacobiano degenera,
  exactamente como el camino por elemento.

Constitutiva — ``MAT_SIG``
--------------------------
``C_t = mat(strain, S_old, S_new, params, C, sigma, flag)``

- ``strain``: deformación en Voigt del proyecto, ``(n_sigma,)``.
- ``S_old`` / ``S_new``: fila de estado committed / trial, ``(n_state,)``,
  con el orden que fija ``Material.STATE_SCHEMA``. El kernel **escribe la
  fila trial completa**.
- ``params``: parámetros escalares del material, ``(n_par,)``, en el
  orden que documenta cada ``Material.batch_params``.
- ``C``: matriz del material (elástica o precalculada), ``(n_sigma, n_sigma)``.
- ``sigma``: salida, esfuerzo ``(n_sigma,)``, escrita in situ.
- ``flag``: vista de un entero (``(1,)``, ``int8``); el kernel escribe
  ``1`` para señalar una incidencia que el material reporta después con
  ``batch_report`` (p. ej. Newton local sin converger).
- Devuelve la tangente algorítmica ``(n_sigma, n_sigma)``, contigua.
"""
from numba import types

F64_1 = types.float64[::1]
F64_2 = types.float64[:, ::1]
F64_3 = types.float64[:, :, ::1]
I8_1 = types.int8[::1]

KIN_SIG = types.float64(F64_1, F64_2, F64_2)
MAT_SIG = F64_2(F64_1, F64_1, F64_1, F64_1, F64_2, F64_1, I8_1)

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
