"""Partición espectral de un esfuerzo en sus partes de tracción y de compresión.

Para un esfuerzo simétrico ``σ̄`` con esfuerzos principales ``σ̄_i`` y
direcciones ``n_i``,

    σ̄⁺ = Σ_i ⟨σ̄_i⟩₊ n_i⊗n_i,        σ̄⁻ = σ̄ − σ̄⁺,

y su derivada (función tensorial isótropa de un tensor; de Souza Neto, Perić
y Owen 2008, apéndice A; Faria, Oliver y Cervera 1998)

    Q⁺ = ∂σ̄⁺/∂σ̄ = Σ_i H(σ̄_i) M_ii⊗M_ii + Σ_{i<j} ½·c_ij S_ij⊗S_ij,

con ``M_ij = n_i⊗n_j``, ``S_ij = M_ij + M_ji``, ``H`` el escalón con
``H(0) = 0`` y ``c_ij = (⟨σ̄_i⟩₊ − ⟨σ̄_j⟩₊)/(σ̄_i − σ̄_j)``. La usan los
materiales de daño de tracción y compresión (``TensionCompressionDamage2D/3D``,
spec §1 y §7).

Convención Voigt de esfuerzo del proyecto: ``[σ_xx, σ_yy, σ_xy]`` en 2D y
``[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz]`` en 3D (ADR 0012). ``Q`` es la
derivada del **vector** Voigt ``σ̄⁺`` respecto del **vector** Voigt ``σ̄``:
cada columna de cortante lleva el factor 2 (variar ``σ̄_xy`` varía a la vez
``σ̄_xy`` y ``σ̄_yx`` del tensor), y ``Q = I`` cuando todos los esfuerzos
principales son positivos.

Esfuerzos principales iguales: con ``H_i = H_j`` el cociente ``c_ij`` vale
exactamente ``H_i`` (los dos positivos o los dos no positivos) y se toma así,
sin dividir; con ``H_i ≠ H_j`` los dos esfuerzos tienen signo distinto y
``σ̄_i − σ̄_j > 0``. No hace falta tolerancia: con ``c_ij = H`` en un
autoespacio repetido, la suma no depende de la base que devuelva la
descomposición dentro de él.

Las funciones devuelven además el mayor esfuerzo principal del plano (2D) o
del sólido (3D) y su proyector ``n_I⊗n_I`` en Voigt de esfuerzo, para la
medida de Rankine sin una segunda descomposición.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit


@njit(cache=True)
def spectral_split_2d(s):
    """Partición del esfuerzo plano ``s = [σ_xx, σ_yy, σ_xy]`` (forma cerrada).

    Devuelve ``(σ̄⁺, Q⁺, σ̄_1, m_1)``: la parte positiva ``(3,)``, su derivada
    ``(3, 3)``, el mayor esfuerzo principal del plano y
    ``m_1 = [n_x², n_y², n_x·n_y]`` de su dirección. La componente fuera del
    plano (``σ̄_zz`` de deformación plana, que es principal) la trata el
    material aparte: no se acopla con el plano en ``Q⁺``.
    """
    sx = s[0]
    sy = s[1]
    sxy = s[2]
    h = 0.5 * (sx - sy)
    r = math.hypot(h, sxy)              # sin desbordamiento por abajo de h²
    mean = 0.5 * (sx + sy)
    s1 = mean + r
    s2 = mean - r
    if r > 0.0:
        c2 = h / r                      # cos 2θ
        sn2 = sxy / r                   # sin 2θ
    else:
        c2 = 1.0                        # autovalor doble: cualquier base
        sn2 = 0.0
    m1 = np.empty(3)
    m2 = np.empty(3)
    s12 = np.empty(3)
    m1[0] = 0.5 * (1.0 + c2)
    m1[1] = 0.5 * (1.0 - c2)
    m1[2] = 0.5 * sn2
    m2[0] = 0.5 * (1.0 - c2)
    m2[1] = 0.5 * (1.0 + c2)
    m2[2] = -0.5 * sn2
    # S_12 = n1⊗n2 + n2⊗n1 en Voigt de esfuerzo, n1 = (cos θ, sin θ),
    # n2 = (−sin θ, cos θ).
    s12[0] = -sn2
    s12[1] = sn2
    s12[2] = c2
    H1 = 1.0 if s1 > 0.0 else 0.0
    H2 = 1.0 if s2 > 0.0 else 0.0
    sp = np.zeros(3)
    Q = np.zeros((3, 3))
    if H2 == 1.0:
        # Los dos positivos: σ̄⁺ = σ̄ exactamente.
        for a in range(3):
            sp[a] = s[a]
            Q[a, a] = 1.0
        return sp, Q, s1, m1
    if H1 == 0.0:
        # Los dos no positivos: σ̄⁺ = 0.
        return sp, Q, s1, m1
    # s1 > 0 ≥ s2: s1 − s2 ≥ s1 > 0.
    c12 = s1 / (s1 - s2)
    w = (1.0, 1.0, 2.0)
    for a in range(3):
        sp[a] = s1 * m1[a]
        for b in range(3):
            Q[a, b] = (m1[a] * m1[b] + 0.5 * c12 * s12[a] * s12[b]) * w[b]
    return sp, Q, s1, m1


@njit(cache=True)
def spectral_split_3d(s):
    """Partición del esfuerzo 3D ``s`` (Voigt 6 de esfuerzo).

    Devuelve ``(σ̄⁺, Q⁺, σ̄_I, m_I)``: la parte positiva ``(6,)``, su derivada
    ``(6, 6)``, el mayor esfuerzo principal y
    ``m_I = [n_x², n_y², n_z², n_x·n_y, n_y·n_z, n_x·n_z]`` de su dirección.
    """
    T = np.empty((3, 3))
    T[0, 0] = s[0]
    T[1, 1] = s[1]
    T[2, 2] = s[2]
    T[0, 1] = s[3]
    T[1, 0] = s[3]
    T[1, 2] = s[4]
    T[2, 1] = s[4]
    T[0, 2] = s[5]
    T[2, 0] = s[5]
    lam, V = np.linalg.eigh(T)          # autovalores en orden ascendente
    m = np.empty((3, 6))
    for i in range(3):
        n0 = V[0, i]
        n1 = V[1, i]
        n2 = V[2, i]
        m[i, 0] = n0 * n0
        m[i, 1] = n1 * n1
        m[i, 2] = n2 * n2
        m[i, 3] = n0 * n1
        m[i, 4] = n1 * n2
        m[i, 5] = n0 * n2
    m_max = m[2].copy()
    sp = np.zeros(6)
    Q = np.zeros((6, 6))
    if lam[0] > 0.0:
        # Todos positivos: σ̄⁺ = σ̄ exactamente.
        for a in range(6):
            sp[a] = s[a]
            Q[a, a] = 1.0
        return sp, Q, lam[2], m_max
    if lam[2] <= 0.0:
        # Todos no positivos: σ̄⁺ = 0.
        return sp, Q, lam[2], m_max
    w = (1.0, 1.0, 1.0, 2.0, 2.0, 2.0)
    H = np.empty(3)
    for i in range(3):
        H[i] = 1.0 if lam[i] > 0.0 else 0.0
        if H[i] == 1.0:
            for a in range(6):
                sp[a] += lam[i] * m[i, a]
                for b in range(6):
                    Q[a, b] += m[i, a] * m[i, b] * w[b]
    sij = np.empty(6)
    for i in range(3):
        for j in range(i + 1, 3):
            if H[i] == H[j]:
                c = H[i]
            else:
                pi = lam[i] if lam[i] > 0.0 else 0.0
                pj = lam[j] if lam[j] > 0.0 else 0.0
                c = (pi - pj) / (lam[i] - lam[j])
            if c == 0.0:
                continue
            a0 = V[0, i]
            a1 = V[1, i]
            a2 = V[2, i]
            b0 = V[0, j]
            b1 = V[1, j]
            b2 = V[2, j]
            sij[0] = 2.0 * a0 * b0
            sij[1] = 2.0 * a1 * b1
            sij[2] = 2.0 * a2 * b2
            sij[3] = a0 * b1 + a1 * b0
            sij[4] = a1 * b2 + a2 * b1
            sij[5] = a0 * b2 + a2 * b0
            for a in range(6):
                for b in range(6):
                    Q[a, b] += 0.5 * c * sij[a] * sij[b] * w[b]
    return sp, Q, lam[2], m_max


__all__ = ["spectral_split_2d", "spectral_split_3d"]
