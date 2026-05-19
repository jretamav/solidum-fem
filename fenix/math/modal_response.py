"""Funciones de cómputo sobre modos: vibración libre y combinación modal
de respuesta espectral (ADR 0009 fases 1 y 7).

Estas funciones operan sobre **modos** y **frecuencias** (no sobre
``ModalResult`` directamente) para que sean composables y testeables sin
construir un ``Domain`` completo. Los métodos de ``ModalResult`` quedan
como wrappers delgados que delegan aquí — los detalles del algoritmo
viven en este módulo (ADR 0009 fase 7, regla D de auditoría arquitectural
2026-05-13 aplicada 2026-05-18).

Contenido:

- :func:`free_vibration` — respuesta temporal por superposición modal
  sin amortiguamiento (ADR 0009 fase 1). Algoritmo previamente en
  ``ModalResult.free_vibration``; el método queda como wrapper.
- :func:`response_spectrum_srss` — combinación SRSS (Square Root of Sum
  of Squares) de respuestas modales máximas. Apropiado cuando las
  frecuencias modales están bien separadas.
- :func:`response_spectrum_cqc` — combinación CQC (Complete Quadratic
  Combination, Der Kiureghian 1980) con coeficientes de correlación
  ``ρ_ij`` derivados del amortiguamiento. Apropiado para sistemas con
  modos cercanos en frecuencia.

Referencias:

- Chopra, A. K. (2017). *Dynamics of Structures*, cap. 13 (respuesta
  espectral; SRSS y CQC).
- Der Kiureghian, A. (1980). "Structural Response to Stationary
  Excitation". *J. Eng. Mech. Div. ASCE*, 106(EM6), 1195-1213.
- Wilson, E. L. (2002). *Three-Dimensional Static and Dynamic Analysis
  of Structures*, cap. 15 (combinación modal CQC).
"""
from __future__ import annotations

from typing import Callable

import numpy as np


def free_vibration(
    modes: np.ndarray,
    frequencies_rad: np.ndarray,
    M,
    u0: np.ndarray,
    u0_dot: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    """Respuesta temporal de vibración libre sin amortiguamiento por
    superposición modal — solución analítica del problema ``M·ü + K·u = 0``
    con condiciones iniciales ``(u(0), u̇(0))``.

    Sobre modos M-ortonormales ``φ_n``, las coordenadas modales
    ``q_n(t) = φ_nᵀ·M·u(t)`` satisfacen ``q̈_n + ω_n²·q_n = 0``, con
    solución::

        q_n(t) = a_n · cos(ω_n·t) + (b_n / ω_n) · sin(ω_n·t)

    donde ``a_n = φ_nᵀ·M·u₀`` y ``b_n = φ_nᵀ·M·u̇₀``. Para modos de cuerpo
    rígido (``ω_n = 0``) la ecuación degenera a ``q̈_n = 0`` y la solución
    es lineal en el tiempo. La respuesta completa se reconstruye como
    ``u(t) = Σ_n q_n(t)·φ_n``.

    **Limitación**: sólo reconstruye la respuesta sobre el subespacio de
    los modos calculados. Si las CI tienen contenido en modos no
    calculados, esa contribución se pierde.

    Parameters
    ----------
    modes : np.ndarray, shape (n_dof, n_modes)
        Modos en columnas, M-ortonormales.
    frequencies_rad : np.ndarray, shape (n_modes,)
        Frecuencias naturales ``ω_n`` en rad/s.
    M : scipy.sparse.spmatrix or ndarray, shape (n_dof, n_dof)
        Matriz de masa global.
    u0, u0_dot : np.ndarray, shape (n_dof,)
        Condiciones iniciales.
    t : np.ndarray, shape (n_t,)
        Vector de instantes temporales.

    Returns
    -------
    np.ndarray, shape (n_dof, n_t)
        ``u(t)`` con una columna por instante.
    """
    u0 = np.asarray(u0, dtype=float).ravel()
    u0_dot = np.asarray(u0_dot, dtype=float).ravel()
    t = np.asarray(t, dtype=float).ravel()

    a = modes.T @ (M @ u0)
    b = modes.T @ (M @ u0_dot)

    omega = np.asarray(frequencies_rad, dtype=float).ravel()
    n_modes = omega.size
    rigid = omega == 0.0
    q = np.empty((n_modes, t.size))
    if np.any(~rigid):
        w = omega[~rigid][:, None]
        tt = t[None, :]
        q[~rigid, :] = (
            a[~rigid, None] * np.cos(w * tt)
            + (b[~rigid, None] / w) * np.sin(w * tt)
        )
    if np.any(rigid):
        q[rigid, :] = a[rigid, None] + b[rigid, None] * t[None, :]

    return modes @ q


# ---------------------------------------------------------------------------
# Respuesta espectral — coeficientes modales
# ---------------------------------------------------------------------------


def participation_factors(
    modes: np.ndarray,
    M,
    direction: np.ndarray,
) -> np.ndarray:
    """Factores de participación modal ``Γ_n = φ_nᵀ·M·r``.

    Donde ``r`` es el vector de excitación rígida del modo (componente
    de aceleración del suelo proyectada a los DOFs globales). En la
    formulación clásica de análisis espectral sísmico, ``r`` corresponde
    al vector de desplazamientos rígidos resultante de un desplazamiento
    unitario del suelo en la dirección de excitación.

    Asumiendo modos M-ortonormales (``ΦᵀMΦ = I``), ``Γ_n`` es la
    proyección de ``r`` sobre cada modo. Carga modal efectiva ``m_n^*``
    (masa efectiva) ``= Γ_n²``.

    Parameters
    ----------
    modes : np.ndarray, shape (n_dof, n_modes)
        Modos M-ortonormales en columnas.
    M : scipy.sparse.spmatrix or ndarray, shape (n_dof, n_dof)
        Matriz de masa global.
    direction : np.ndarray, shape (n_dof,)
        Vector ``r`` de excitación rígida. Para excitación sísmica
        unidireccional, los DOFs traslacionales en la dirección de
        excitación valen 1, el resto 0.

    Returns
    -------
    np.ndarray, shape (n_modes,)
        Factores de participación ``Γ_n``.
    """
    direction = np.asarray(direction, dtype=float).ravel()
    return modes.T @ (M @ direction)


def response_spectrum_srss(
    modes: np.ndarray,
    frequencies_rad: np.ndarray,
    M,
    direction: np.ndarray,
    spectrum_fn: Callable[[float], float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Combinación SRSS (Square Root of Sum of Squares) de respuestas
    modales máximas.

    Para cada modo ``n``:

    .. math::

        u_n^\\max = \\Gamma_n \\cdot \\phi_n \\cdot S_d(\\omega_n)

    La respuesta combinada es

    .. math::

        u^\\max = \\sqrt{\\sum_n (u_n^\\max)^2}   \\text{(elemento a elemento)}

    Apropiado cuando las frecuencias modales están **bien separadas**
    (cociente típicamente > 1.3). Para modos cercanos, usar
    :func:`response_spectrum_cqc`.

    Parameters
    ----------
    modes : np.ndarray, shape (n_dof, n_modes)
        Modos M-ortonormales en columnas.
    frequencies_rad : np.ndarray, shape (n_modes,)
        Frecuencias naturales ``ω_n`` (rad/s). Modos de cuerpo rígido
        (``ω_n = 0``) no contribuyen al espectro de respuesta clásico —
        se filtran con un warning interno.
    M : scipy.sparse.spmatrix or ndarray
        Matriz de masa global.
    direction : np.ndarray, shape (n_dof,)
        Vector de excitación rígida (ver :func:`participation_factors`).
    spectrum_fn : callable
        Función ``S_d(ω) -> desplazamiento espectral`` (m si las
        unidades del modelo son SI). El usuario convierte aceleración
        espectral ``Sa`` a desplazamiento espectral ``Sd = Sa/ω²`` si
        su input es ``Sa``.

    Returns
    -------
    u_combined : np.ndarray, shape (n_dof,)
        Respuesta máxima combinada SRSS.
    u_per_mode : np.ndarray, shape (n_dof, n_modes)
        Contribución modal individual ``Γ_n·φ_n·S_d(ω_n)`` (con signo).
    gamma : np.ndarray, shape (n_modes,)
        Factores de participación modal.
    """
    modes = np.asarray(modes, dtype=float)
    omega = np.asarray(frequencies_rad, dtype=float).ravel()
    gamma = participation_factors(modes, M, direction)

    Sd = np.array([float(spectrum_fn(w)) if w > 0.0 else 0.0
                    for w in omega])
    # u_per_mode[:, n] = γ_n · φ_n · S_d(ω_n)
    u_per_mode = modes * (gamma * Sd)[None, :]
    u_combined = np.sqrt(np.sum(u_per_mode * u_per_mode, axis=1))
    return u_combined, u_per_mode, gamma


def _cqc_correlation_matrix(
    frequencies_rad: np.ndarray,
    damping: float,
) -> np.ndarray:
    """Matriz de correlación ``ρ_ij`` de Der Kiureghian (1980) para CQC.

    Fórmula asumiendo amortiguamiento ``ξ`` igual en todos los modos
    (caso usual en espectros normativos):

    .. math::

        \\rho_{ij} = \\frac{8 \\xi^2 (1+r) r^{3/2}}{(1-r^2)^2 + 4\\xi^2 r (1+r)^2}

    con ``r = ω_i / ω_j``. ``ρ_ii = 1``; ``ρ_ij → δ_ij`` cuando
    ``ξ → 0`` o ``r → 0`` (modos bien separados → SRSS).
    """
    omega = np.asarray(frequencies_rad, dtype=float).ravel()
    n = omega.size
    rho = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if omega[i] == 0.0 or omega[j] == 0.0:
                rho[i, j] = 1.0 if i == j else 0.0
                continue
            r = omega[i] / omega[j]
            num = 8.0 * damping**2 * (1.0 + r) * r**1.5
            den = (1.0 - r**2) ** 2 + 4.0 * damping**2 * r * (1.0 + r) ** 2
            rho[i, j] = num / den
    return rho


def response_spectrum_cqc(
    modes: np.ndarray,
    frequencies_rad: np.ndarray,
    M,
    direction: np.ndarray,
    spectrum_fn: Callable[[float], float],
    damping: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Combinación CQC (Complete Quadratic Combination, Der Kiureghian
    1980) de respuestas modales máximas.

    .. math::

        u^\\max_k = \\sqrt{\\sum_i \\sum_j \\rho_{ij} \\cdot u_{k,i}^\\max \\cdot u_{k,j}^\\max}

    con ``ρ_ij`` los coeficientes de correlación dependientes del
    amortiguamiento. Para modos bien separados ``ρ_ij ≈ δ_ij`` y CQC
    degenera a SRSS.

    Parameters
    ----------
    modes, frequencies_rad, M, direction, spectrum_fn
        Idénticos a :func:`response_spectrum_srss`.
    damping : float
        Amortiguamiento modal ``ξ`` (mismo para todos los modos). Típico
        en espectros normativos: ``0.05`` (5%).

    Returns
    -------
    u_combined : np.ndarray, shape (n_dof,)
        Respuesta máxima combinada CQC.
    u_per_mode : np.ndarray, shape (n_dof, n_modes)
        Contribución modal individual con signo.
    gamma : np.ndarray, shape (n_modes,)
        Factores de participación modal.
    """
    if damping < 0.0:
        raise ValueError(
            f"response_spectrum_cqc: damping={damping} negativo no admitido."
        )

    modes = np.asarray(modes, dtype=float)
    omega = np.asarray(frequencies_rad, dtype=float).ravel()
    gamma = participation_factors(modes, M, direction)

    Sd = np.array([float(spectrum_fn(w)) if w > 0.0 else 0.0
                    for w in omega])
    u_per_mode = modes * (gamma * Sd)[None, :]

    rho = _cqc_correlation_matrix(omega, damping)
    # u_combined_k² = Σ_i Σ_j ρ_ij · u_k,i · u_k,j
    # Vectorizado: para cada DOF k, u_combined_k² = u_k @ rho @ u_k.
    quad = np.einsum("ki,ij,kj->k", u_per_mode, rho, u_per_mode)
    quad = np.maximum(quad, 0.0)  # protección contra error de redondeo
    u_combined = np.sqrt(quad)
    return u_combined, u_per_mode, gamma
