"""Kernels de familia: bucles compilados sobre elementos y puntos de Gauss (ADR 0014).

Son el bucle exterior del ensamblaje, compilado. No contienen física
propia: la cinemática (``kin``) y la constitutiva (``mat``) llegan como
funciones tipadas (ver ``signatures.py``) y son **las mismas funciones
compiladas que usa el camino por elemento**. Lo único que hacen estos
kernels es recorrer, integrar y volcar en arreglos.

Dos variantes con la misma firma y la misma aritmética:

- :func:`solid_family_kernel` — serie.
- :func:`solid_family_kernel_parallel` — ``prange`` sobre bloques de
  ``PARALLEL_BLOCK`` elementos (Numba ``parallel=True``). Cada bloque
  asigna su propio espacio de trabajo, así que no hay estado compartido
  entre hilos; cada elemento escribe sólo sus filas de ``K_out``,
  ``F_out``, ``S_out``, ``sig_out`` y ``flags``. El orden de las sumas
  dentro de un elemento es idéntico al de la variante serie, por lo que
  el resultado es **bit a bit el mismo** con cualquier número de hilos.

Señalización sin excepciones
----------------------------
Una excepción lanzada dentro de una región ``prange`` de Numba se pierde
(la iteración aborta en silencio y el bucle continúa). Por eso las
cinemáticas por lotes **no lanzan**: devuelven ``det J ≤ 0`` cuando el
jacobiano degenera, y el kernel marca ``flags[fila] = -1`` y salta ese
punto. Quien llama (``Family.evaluate``) convierte la marca en el
``ValueError`` de siempre, con el id del elemento. Las incidencias del
material (Newton local sin converger) usan ``flags[fila] = 1``. Ambas
variantes, serie y paralela, siguen este protocolo para comportarse igual.

Aritmética con bucles explícitos en vez de ``np.dot``: para matrices de
``3×8`` a ``6×81`` la llamada a BLAS cuesta más que la operación, y
además evita asignar temporales por punto de Gauss.

Por qué reciben coordenadas de referencia y desplazamientos y no una
``B`` cacheada: principio 4 de la propuesta (nada derivado de
``x = X + u`` se cachea entre iteraciones). Cachear ``B`` duplicaría la
memoria por elemento y fijaría una hipótesis de geometría que las
formulaciones lagrangianas futuras no cumplen.
"""
import numpy as np
from numba import njit, prange, types

from solidum.math.batch.signatures import FAMILY_SIG, GAUSS_SIG, REDUCE_SIG

# Elementos por bloque en la variante paralela. Un bloque amortiza la
# asignación del espacio de trabajo (cinco arreglos pequeños) y da al
# planificador de Numba unidades de trabajo del orden de decenas de
# microsegundos. Vive aquí y no en ``constants.py`` porque Numba congela
# las globales al compilar y el caché en disco sólo se invalida cuando
# cambia **este** archivo.
PARALLEL_BLOCK = 32

# Marca de jacobiano degenerado en ``flags`` (las incidencias del material
# usan +1).
FLAG_BAD_JACOBIAN = -1


@njit(cache=True)
def _element_integrals(kin, mat, e, X, u, pts, w, scale, S_in, S_out, params, C,
                       K, F, sig_out, flags, B, CB, strain, sigma, C_out):
    """Integra ``K_e`` y ``F_e`` del elemento ``e`` sobre sus puntos de Gauss.

    Cuerpo común a las dos variantes del kernel de familia. Los bucles
    internos recorren siempre el último índice (contiguo) de ``B``, ``CB``
    y ``K`` para que LLVM los vectorice; el orden de las sumas es fijo y
    el mismo en cada evaluación, como exige la equivalencia con el camino
    por elemento.
    """
    n_gp = pts.shape[0]
    n_dof = u.shape[1]
    n_sig = sig_out.shape[1]

    for i in range(n_dof):
        F[i] = 0.0
        for j in range(n_dof):
            K[i, j] = 0.0

    for g in range(n_gp):
        row = e * n_gp + g
        detJ = kin(pts[g], X[e], B)
        if detJ <= 0.0:
            flags[row] = FLAG_BAD_JACOBIAN
            continue

        # strain = B · u_e
        for a in range(n_sig):
            s = 0.0
            for i in range(n_dof):
                s += B[a, i] * u[e, i]
            strain[a] = s

        Ct = mat(strain, S_in[row], S_out[row], params, C, sigma, C_out,
                 flags[row:row + 1])

        dV = detJ * w[g] * scale[e]

        # CB = (C_t · B) · dV
        for a in range(n_sig):
            for j in range(n_dof):
                CB[a, j] = 0.0
            for b in range(n_sig):
                c = Ct[a, b]
                for j in range(n_dof):
                    CB[a, j] += c * B[b, j]
            for j in range(n_dof):
                CB[a, j] *= dV

        # K += Bᵀ · CB ;  F += Bᵀ · (σ · dV)
        for a in range(n_sig):
            sa = sigma[a] * dV
            for i in range(n_dof):
                bai = B[a, i]
                F[i] += bai * sa
                for j in range(n_dof):
                    K[i, j] += bai * CB[a, j]

        for a in range(n_sig):
            sig_out[row, a] = sigma[a]


@njit(FAMILY_SIG, cache=True)
def solid_family_kernel(kin, mat, X, u, pts, w, scale, S_in, S_out, params, C,
                        K_out, F_out, sig_out, flags):
    n_elem = X.shape[0]
    n_dof = u.shape[1]
    n_sig = sig_out.shape[1]

    B = np.zeros((n_sig, n_dof))
    CB = np.zeros((n_sig, n_dof))
    strain = np.zeros(n_sig)
    sigma = np.zeros(n_sig)
    C_out = np.zeros((n_sig, n_sig))

    for e in range(n_elem):
        _element_integrals(kin, mat, e, X, u, pts, w, scale, S_in, S_out, params, C,
                           K_out[e], F_out[e], sig_out, flags, B, CB, strain, sigma, C_out)


@njit(FAMILY_SIG, cache=True, parallel=True)
def solid_family_kernel_parallel(kin, mat, X, u, pts, w, scale, S_in, S_out, params, C,
                                 K_out, F_out, sig_out, flags):
    n_elem = X.shape[0]
    n_dof = u.shape[1]
    n_sig = sig_out.shape[1]
    n_blocks = (n_elem + PARALLEL_BLOCK - 1) // PARALLEL_BLOCK

    for blk in prange(n_blocks):
        B = np.zeros((n_sig, n_dof))
        CB = np.zeros((n_sig, n_dof))
        strain = np.zeros(n_sig)
        sigma = np.zeros(n_sig)
        C_out = np.zeros((n_sig, n_sig))
        e0 = blk * PARALLEL_BLOCK
        e1 = min(n_elem, e0 + PARALLEL_BLOCK)
        for e in range(e0, e1):
            _element_integrals(kin, mat, e, X, u, pts, w, scale, S_in, S_out, params, C,
                               K_out[e], F_out[e], sig_out, flags, B, CB, strain, sigma, C_out)


# ---------------------------------------------------------------------------
# Post-proceso por familia: ε y σ por punto de Gauss (ADR 0014 §9)
# ---------------------------------------------------------------------------
#
# Es ``compute_gauss_state`` sobre toda la familia: misma cinemática, misma
# constitutiva evaluada desde el estado **committed**; no integra nada y
# no escribe el trial (``S_scratch`` recibe la fila trial que el material
# produce y se descarta).

@njit(cache=True)
def _element_gauss(kin, mat, e, X, u, pts, S_in, S_scratch, params, C,
                   eps_out, sig_out, flags, B, strain, sigma, C_out):
    n_gp = pts.shape[0]
    n_dof = u.shape[1]
    n_sig = sig_out.shape[1]
    for g in range(n_gp):
        row = e * n_gp + g
        detJ = kin(pts[g], X[e], B)
        if detJ <= 0.0:
            flags[row] = FLAG_BAD_JACOBIAN
            continue
        for a in range(n_sig):
            s = 0.0
            for i in range(n_dof):
                s += B[a, i] * u[e, i]
            strain[a] = s
        mat(strain, S_in[row], S_scratch[row], params, C, sigma, C_out,
            flags[row:row + 1])
        for a in range(n_sig):
            eps_out[row, a] = strain[a]
            sig_out[row, a] = sigma[a]


@njit(GAUSS_SIG, cache=True)
def solid_family_gauss_kernel(kin, mat, X, u, pts, S_in, S_scratch, params, C,
                              eps_out, sig_out, flags):
    n_elem = X.shape[0]
    n_dof = u.shape[1]
    n_sig = sig_out.shape[1]
    B = np.zeros((n_sig, n_dof))
    strain = np.zeros(n_sig)
    sigma = np.zeros(n_sig)
    C_out = np.zeros((n_sig, n_sig))
    for e in range(n_elem):
        _element_gauss(kin, mat, e, X, u, pts, S_in, S_scratch, params, C,
                       eps_out, sig_out, flags, B, strain, sigma, C_out)


@njit(GAUSS_SIG, cache=True, parallel=True)
def solid_family_gauss_kernel_parallel(kin, mat, X, u, pts, S_in, S_scratch, params, C,
                                       eps_out, sig_out, flags):
    n_elem = X.shape[0]
    n_dof = u.shape[1]
    n_sig = sig_out.shape[1]
    n_blocks = (n_elem + PARALLEL_BLOCK - 1) // PARALLEL_BLOCK
    for blk in prange(n_blocks):
        B = np.zeros((n_sig, n_dof))
        strain = np.zeros(n_sig)
        sigma = np.zeros(n_sig)
        C_out = np.zeros((n_sig, n_sig))
        e0 = blk * PARALLEL_BLOCK
        e1 = min(n_elem, e0 + PARALLEL_BLOCK)
        for e in range(e0, e1):
            _element_gauss(kin, mat, e, X, u, pts, S_in, S_scratch, params, C,
                           eps_out, sig_out, flags, B, strain, sigma, C_out)


# ---------------------------------------------------------------------------
# Reducción COO → CSR con el mapa inverso cacheado (ver Assembler._build_csr_map)
# ---------------------------------------------------------------------------
#
# ``out[i] = Σ_k data[src_idx[k]]`` para ``k`` en ``src_ptr[i]..src_ptr[i+1]``:
# cada entrada CSR suma sus fuentes COO en orden creciente de índice COO,
# el mismo orden en que ``np.bincount`` las acumularía, así que el resultado
# es bit a bit el de ``bincount`` pero sin carreras (cada hilo escribe su
# ``out[i]``) y paralelizable.

@njit(REDUCE_SIG, cache=True)
def coo_to_csr_reduce(data, src_ptr, src_idx, out):
    for i in range(out.shape[0]):
        s = 0.0
        for k in range(src_ptr[i], src_ptr[i + 1]):
            s += data[src_idx[k]]
        out[i] = s


@njit(REDUCE_SIG, cache=True, parallel=True)
def coo_to_csr_reduce_parallel(data, src_ptr, src_idx, out):
    for i in prange(out.shape[0]):
        s = 0.0
        for k in range(src_ptr[i], src_ptr[i + 1]):
            s += data[src_idx[k]]
        out[i] = s
