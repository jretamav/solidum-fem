"""Medidas de deformación equivalente del daño isótropo escalar 2D y 3D.

``IsotropicDamage2D`` e ``IsotropicDamage3D`` gobiernan el daño con un
escalar ``ε_eq(ε)`` y la variable histórica ``κ = max(κ, ε_eq)``. Desde el
2026-09-29 la medida es configurable (parámetro ``equivalent_strain``,
aprobado por el usuario) entre tres medidas de la literatura, todas
calibradas para que ``ε_eq = ε`` en tracción uniaxial de esfuerzo, de modo
que ``κ₀ = f_t/E``:

- ``"energy"`` (Simó y Ju 1987), por defecto:
  ``ε_eq = √(εᵀ·C_e·ε / E)``. Simétrica en tracción y compresión; la
  tangente consistente que produce es simétrica.
- ``"mazars"`` (Mazars 1986): ``ε_eq = √(Σ_i ⟨ε_i⟩₊²)`` con ``ε_i`` las
  deformaciones principales del tensor 3D completo. Sólo las extensiones
  dañan: en compresión uniaxial de esfuerzo ``ε_eq = √2·ν·|ε|``.
- ``"modified_von_mises"`` (de Vree, Brekelmans y van Gils 1995): función
  de ``I₁`` y ``J₂`` de la deformación con ``k = f_c/f_t``; en compresión
  uniaxial de esfuerzo ``ε_eq = |ε|/k``.

Por qué un módulo propio: las tres medidas son las mismas en 2D y en 3D (el
2D es la restricción del 3D con ``ε_zz = c_zz·(ε_xx + ε_yy)`` y
``γ_yz = γ_xz = 0``), así que viven aquí una sola vez y los dos materiales
las consumen (Reglas.md §1: dos casos reales). La norma del tensor
``√(ε:ε)`` que usaban antes los dos materiales se retiró: no es estándar y
está mal calibrada (``ε_eq = ε·√(1+2ν²)`` en tracción uniaxial).

Contrato de las funciones compiladas
------------------------------------
Devuelven ``(ε_eq, v, s)`` con la derivada ``∂ε_eq/∂ε = s·v`` respecto de la
deformación en Voigt *engineering* del proyecto (``γ_ij = 2·ε_ij``). La
factorización en un vector y un escalar existe por la medida energética:
ahí ``v = C_e·ε`` (el esfuerzo efectivo) y ``s = 1/(E·ε_eq)``, y la tangente
``C_alg = (1−d)·C_e − (∂d/∂κ)·s·(C_e·ε)⊗v`` sale **exactamente** simétrica
en aritmética de coma flotante (``σ_i·σ_j = σ_j·σ_i``). Con ``ε_eq = 0`` la
derivada no existe y se devuelve ``s = 0``; el material no la usa, porque
la rama consistente exige ``ε_eq > κ ≥ κ₀ > 0``.

El código numérico de la medida (``ENERGY``, ``MAZARS``,
``MODIFIED_VON_MISES``) se valida al construir el material
(:func:`parse_equivalent_strain`); los núcleos compilados no lanzan
excepciones (una excepción dentro del kernel por lotes paralelo se
perdería, ADR 0014) y tratan cualquier código distinto de los dos primeros
como von Mises modificado.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit

# Códigos numéricos de la medida: viajan en ``batch_params`` (ADR 0014).
ENERGY = 0
MAZARS = 1
MODIFIED_VON_MISES = 2

EQUIVALENT_STRAIN_CODES = {
    "energy": ENERGY,
    "mazars": MAZARS,
    "modified_von_mises": MODIFIED_VON_MISES,
}
DEFAULT_EQUIVALENT_STRAIN = "energy"


def parse_equivalent_strain(owner: str, equivalent_strain, k) -> tuple[str, int, float | None]:
    """Valida la pareja ``(equivalent_strain, k)`` al construir el material.

    Devuelve ``(nombre, código, k)``, con ``k = None`` salvo en
    ``"modified_von_mises"``. Lanza ``ValueError`` con el nombre del material
    (``owner``) y del parámetro si la medida es desconocida, si falta ``k``
    con von Mises modificado, si se da ``k`` con otra medida o si ``k`` no es
    un número finito positivo.
    """
    validas = sorted(EQUIVALENT_STRAIN_CODES)
    if not isinstance(equivalent_strain, str) or equivalent_strain not in EQUIVALENT_STRAIN_CODES:
        raise ValueError(
            f"{owner}: equivalent_strain={equivalent_strain!r} no reconocida. "
            f"Valores válidos: {validas}."
        )
    code = EQUIVALENT_STRAIN_CODES[equivalent_strain]
    if code != MODIFIED_VON_MISES:
        if k is not None:
            raise ValueError(
                f"{owner}: el parámetro k (= f_c/f_t) sólo se admite con "
                f"equivalent_strain='modified_von_mises'; se recibió k={k!r} con "
                f"equivalent_strain={equivalent_strain!r}."
            )
        return equivalent_strain, code, None
    if k is None:
        raise ValueError(
            f"{owner}: equivalent_strain='modified_von_mises' exige el parámetro "
            f"k = f_c/f_t (cociente entre las resistencias uniaxiales a compresión "
            f"y a tracción), k > 0."
        )
    try:
        k_val = float(k)
    except (TypeError, ValueError):
        raise ValueError(f"{owner}: k={k!r} debe ser un número > 0.") from None
    if isinstance(k, bool) or not math.isfinite(k_val) or k_val <= 0.0:
        raise ValueError(f"{owner}: k debe ser un número finito > 0 (recibido {k!r}).")
    return equivalent_strain, code, k_val


# ---------------------------------------------------------------------------
# Núcleos compilados
# ---------------------------------------------------------------------------

@njit(cache=True)
def _modified_von_mises_3d(e, nu, k):
    """von Mises modificado sobre el estado 3D ``e`` (Voigt 6D engineering).

    ``ε_eq = (k−1)/(2k(1−2ν))·I₁ + 1/(2k)·√(b²·I₁² + c·J₂)`` con
    ``b = (k−1)/(1−2ν)`` y ``c = 12k/(1+ν)²``. Derivada:
    ``∂ε_eq/∂ε = [b + b²·I₁/S]/(2k)·∂I₁/∂ε + c/(4k·S)·∂J₂/∂ε``, con
    ``S = √(b²·I₁² + c·J₂)``, ``∂I₁/∂ε = [1, 1, 1, 0, 0, 0]`` y
    ``∂J₂/∂ε = [e_xx, e_yy, e_zz, γ_xy/2, γ_yz/2, γ_xz/2]`` (``e`` desviadora).
    """
    I1 = e[0] + e[1] + e[2]
    J2 = ((e[0] - e[1]) ** 2 + (e[1] - e[2]) ** 2 + (e[2] - e[0]) ** 2) / 6.0 \
        + 0.25 * (e[3] * e[3] + e[4] * e[4] + e[5] * e[5])
    b = (k - 1.0) / (1.0 - 2.0 * nu)
    c = 12.0 * k / ((1.0 + nu) * (1.0 + nu))
    S = math.sqrt(b * b * I1 * I1 + c * J2)
    eps_eq = (b * I1 + S) / (2.0 * k)
    if eps_eq < 0.0:
        # Sólo por redondeo (compresión hidrostática con k > 1: ε_eq = 0).
        eps_eq = 0.0
    v = np.zeros(6)
    if S > 0.0:
        g_I1 = (b + b * b * I1 / S) / (2.0 * k)
        g_J2 = c / (4.0 * k * S)
        mean = I1 / 3.0
        v[0] = g_I1 + g_J2 * (e[0] - mean)
        v[1] = g_I1 + g_J2 * (e[1] - mean)
        v[2] = g_I1 + g_J2 * (e[2] - mean)
        v[3] = g_J2 * 0.5 * e[3]
        v[4] = g_J2 * 0.5 * e[4]
        v[5] = g_J2 * 0.5 * e[5]
    return eps_eq, v


@njit(cache=True)
def _mazars_3d(e):
    """Mazars sobre el estado 3D ``e``: ``ε_eq = ‖ε₊‖`` y ``v = ε₊`` en Voigt.

    ``ε₊ = Σ_i ⟨ε_i⟩₊·n_i⊗n_i`` es la parte positiva del tensor de
    deformaciones; ``ε_eq² = Σ_i ⟨ε_i⟩₊²``. Como ``g(λ) = Σ⟨λ_i⟩₊²`` es una
    función simétrica C¹ de los valores principales, su gradiente
    ``2·Σ⟨ε_i⟩₊·n_i⊗n_i = 2·ε₊`` está bien definido también con valores
    principales repetidos (los proyectores de un autoespacio repetido llevan
    el mismo peso). En Voigt engineering la componente de ``γ_ab`` del
    gradiente es la componente tensorial ``(ε₊)_ab`` (``∂ε_i/∂γ_ab = n_a·n_b``).
    """
    v = np.zeros(6)
    T = np.empty((3, 3))
    T[0, 0] = e[0]
    T[1, 1] = e[1]
    T[2, 2] = e[2]
    T[0, 1] = 0.5 * e[3]
    T[1, 0] = 0.5 * e[3]
    T[1, 2] = 0.5 * e[4]
    T[2, 1] = 0.5 * e[4]
    T[0, 2] = 0.5 * e[5]
    T[2, 0] = 0.5 * e[5]
    w, V = np.linalg.eigh(T)          # autovalores en orden ascendente
    if w[0] >= 0.0:
        # Todas las deformaciones principales ≥ 0: ε₊ = ε.
        for a in range(3):
            v[a] = e[a]
        for a in range(3, 6):
            v[a] = 0.5 * e[a]
        s2 = e[0] * e[0] + e[1] * e[1] + e[2] * e[2] \
            + 0.5 * (e[3] * e[3] + e[4] * e[4] + e[5] * e[5])
    elif w[2] <= 0.0:
        # Todas ≤ 0: sin extensiones, ε₊ = 0.
        s2 = 0.0
    else:
        s2 = 0.0
        for i in range(3):
            lam = w[i]
            if lam > 0.0:
                s2 += lam * lam
                n0 = V[0, i]
                n1 = V[1, i]
                n2 = V[2, i]
                v[0] += lam * n0 * n0
                v[1] += lam * n1 * n1
                v[2] += lam * n2 * n2
                v[3] += lam * n0 * n1
                v[4] += lam * n1 * n2
                v[5] += lam * n0 * n2
    return math.sqrt(s2), v


@njit(cache=True)
def equivalent_strain_3d(strain, sigma_eff, E, nu, measure, k):
    """``(ε_eq, v, s)`` del estado 3D, con ``∂ε_eq/∂ε = s·v`` (Voigt 6D).

    ``sigma_eff = C_e·ε`` es el esfuerzo efectivo (lo usa la medida
    energética); ``nu`` y ``k`` los usa von Mises modificado.
    """
    if measure == ENERGY:
        w = 0.0
        for a in range(6):
            w += strain[a] * sigma_eff[a]
        eps_eq = math.sqrt(max(w, 0.0) / E)
        v = sigma_eff.copy()
        s = 1.0 / (E * eps_eq) if eps_eq > 0.0 else 0.0
        return eps_eq, v, s
    if measure == MAZARS:
        eps_eq, v = _mazars_3d(strain)
        s = 1.0 / eps_eq if eps_eq > 0.0 else 0.0
        return eps_eq, v, s
    eps_eq, v = _modified_von_mises_3d(strain, nu, k)
    return eps_eq, v, 1.0


@njit(cache=True)
def equivalent_strain_2d(strain, sigma_eff, E, nu, c_zz, measure, k):
    """``(ε_eq, v, s)`` del estado plano, con ``∂ε_eq/∂ε = s·v`` (Voigt 2D).

    El estado 3D asociado es ``[ε_xx, ε_yy, ε_zz, γ_xy, 0, 0]`` con
    ``ε_zz = c_zz·(ε_xx + ε_yy)``: ``c_zz = −ν/(1−ν)`` en tensión plana (con
    daño escalar isótropo ``σ_zz = 0`` no depende de ``d``) y ``0`` en
    deformación plana. Por la regla de la cadena,
    ``∂ε_eq/∂ε_xx = ∂ε_eq/∂ε_xx|₃D + c_zz·∂ε_eq/∂ε_zz|₃D`` (igual en ``yy``).

    - Energética: ``εᵀ·C_e·ε`` con la ``C_e`` plana del material. En tensión
      plana vale ``σ:ε`` del 3D con ``σ_zz = 0`` y en deformación plana el
      término ``σ_zz·ε_zz`` es nulo: el 2D es la restricción exacta del 3D
      sin tratar ``ε_zz`` aparte (y ``v = C_e·ε`` ya incluye la cadena).
    - Mazars: forma cerrada; ``z`` es dirección principal (``γ_xz = γ_yz = 0``)
      y las otras dos salen del bloque ``2×2`` del plano.
    - von Mises modificado: sobre el estado 3D asociado y regla de la cadena.
    """
    if measure == ENERGY:
        w = strain[0] * sigma_eff[0] + strain[1] * sigma_eff[1] + strain[2] * sigma_eff[2]
        eps_eq = math.sqrt(max(w, 0.0) / E)
        v = sigma_eff.copy()
        s = 1.0 / (E * eps_eq) if eps_eq > 0.0 else 0.0
        return eps_eq, v, s

    x = strain[0]
    y = strain[1]
    g = strain[2]
    ezz = c_zz * (x + y)
    v = np.zeros(3)

    if measure == MAZARS:
        q = 0.5 * g                     # ε_xy tensorial
        h = 0.5 * (x - y)
        r = math.sqrt(h * h + q * q)
        m = 0.5 * (x + y)
        e1 = m + r                      # principales del plano, e1 ≥ e2
        e2 = m - r
        if e2 >= 0.0:
            pxx = x
            pyy = y
            pxy = q
            s2 = x * x + y * y + 2.0 * q * q
        elif e1 <= 0.0:
            pxx = 0.0
            pyy = 0.0
            pxy = 0.0
            s2 = 0.0
        else:
            # e1 > 0 > e2 ⇒ r > 0. ε₊ del plano = e1·n1⊗n1 con
            # n1⊗n1 = (A − e2·I)/(e1 − e2); A − e2·I = [[h + r, q], [q, r − h]]
            # evaluado sin cancelación: (h + r)(r − h) = q².
            if h >= 0.0:
                a = h + r
                b = q * q / a if a > 0.0 else 0.0
            else:
                b = r - h
                a = q * q / b
            f = e1 / (2.0 * r)
            pxx = f * a
            pyy = f * b
            pxy = f * q
            s2 = e1 * e1
        zp = ezz if ezz > 0.0 else 0.0
        s2 += zp * zp
        eps_eq = math.sqrt(s2)
        v[0] = pxx + c_zz * zp
        v[1] = pyy + c_zz * zp
        v[2] = pxy
        s = 1.0 / eps_eq if eps_eq > 0.0 else 0.0
        return eps_eq, v, s

    e3 = np.zeros(6)
    e3[0] = x
    e3[1] = y
    e3[2] = ezz
    e3[3] = g
    eps_eq, v3 = _modified_von_mises_3d(e3, nu, k)
    v[0] = v3[0] + c_zz * v3[2]
    v[1] = v3[1] + c_zz * v3[2]
    v[2] = v3[3]
    return eps_eq, v, 1.0
