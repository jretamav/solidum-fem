"""Kernel de familia: bucle compilado sobre elementos y puntos de Gauss (ADR 0014).

Es el bucle exterior del ensamblaje, compilado. No contiene física
propia: la cinemática (``kin``) y la constitutiva (``mat``) llegan como
funciones tipadas (ver ``signatures.py``) y son **las mismas funciones
compiladas que usa el camino por elemento**. Lo único que hace este
kernel es recorrer, integrar y volcar en arreglos.

Aritmética con bucles explícitos en vez de ``np.dot``: para matrices de
``3×8`` a ``6×81`` la llamada a BLAS cuesta más que la operación, y
además evita asignar temporales por punto de Gauss.

Por qué recibe coordenadas de referencia y desplazamientos y no una
``B`` cacheada: principio 4 de la propuesta (nada derivado de
``x = X + u`` se cachea entre iteraciones). Cachear ``B`` duplicaría la
memoria por elemento y fijaría una hipótesis de geometría que las
formulaciones lagrangianas futuras no cumplen.
"""
import numpy as np
from numba import njit

from solidum.math.batch.signatures import FAMILY_SIG


@njit(FAMILY_SIG, cache=True)
def solid_family_kernel(kin, mat, X, u, pts, w, scale, S_in, S_out, params, C,
                        K_out, F_out, sig_out, flags):
    n_elem = X.shape[0]
    n_gp = pts.shape[0]
    n_dof = u.shape[1]
    n_sig = sig_out.shape[1]

    B = np.zeros((n_sig, n_dof))
    CB = np.zeros((n_sig, n_dof))
    strain = np.zeros(n_sig)
    sigma = np.zeros(n_sig)

    for e in range(n_elem):
        K = K_out[e]
        F = F_out[e]
        for i in range(n_dof):
            F[i] = 0.0
            for j in range(n_dof):
                K[i, j] = 0.0

        for g in range(n_gp):
            detJ = kin(pts[g], X[e], B)

            # strain = B · u_e
            for a in range(n_sig):
                s = 0.0
                for i in range(n_dof):
                    s += B[a, i] * u[e, i]
                strain[a] = s

            row = e * n_gp + g
            Ct = mat(strain, S_in[row], S_out[row], params, C, sigma,
                     flags[row:row + 1])

            dV = detJ * w[g] * scale[e]

            # CB = C_t · B
            for a in range(n_sig):
                for j in range(n_dof):
                    s = 0.0
                    for b in range(n_sig):
                        s += Ct[a, b] * B[b, j]
                    CB[a, j] = s

            # K += Bᵀ · CB · dV ;  F += Bᵀ · σ · dV
            for i in range(n_dof):
                fi = 0.0
                for a in range(n_sig):
                    fi += B[a, i] * sigma[a]
                F[i] += fi * dV
                for j in range(n_dof):
                    s = 0.0
                    for a in range(n_sig):
                        s += B[a, i] * CB[a, j]
                    K[i, j] += s * dV

            for a in range(n_sig):
                sig_out[row, a] = sigma[a]
