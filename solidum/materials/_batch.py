"""Kernels puntuales compartidos del camino por lotes (ADR 0014).

Todo material del catálogo que participa en el ensamblaje por lotes
expone una función ``@njit`` con la firma ``MAT_SIG`` de
``solidum.math.batch.signatures``::

    C_t = kernel(strain, S_old, S_new, params, C, sigma, flag)

Este módulo contiene el kernel de los materiales **lineales sin
historia** (elásticos isótropos, ortótropo, conducción de Fourier), que
comparten la misma operación ``σ = C·ε`` con tangente constante. Los
materiales con historia definen su adaptador junto a su return mapping,
en su propio módulo, para que la física siga viviendo en un solo sitio.
"""
from numba import njit


@njit(cache=True)
def linear_material_kernel(strain, S_old, S_new, params, C, sigma, flag):
    """``σ = C·ε``, ``C_t = C``. Sin variables internas (``n_state = 0``)."""
    n = strain.shape[0]
    for a in range(n):
        s = 0.0
        for b in range(n):
            s += C[a, b] * strain[b]
        sigma[a] = s
    return C
