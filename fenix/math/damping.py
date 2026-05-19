"""Modelos de amortiguamiento para análisis dinámico (ADR 0009).

Fase 3: amortiguamiento Rayleigh proporcional ``C = α·M + β·K``. Los
coeficientes ``α`` y ``β`` se calibran a partir de dos pares modales
``(ξ_n, ω_n)`` — típicamente las dos primeras frecuencias propias
obtenidas de un análisis modal previo, lo que cierra el lazo con la
fase 1.

Esquemas más generales (amortiguamiento modal por modo, Caughey series,
amortiguamiento local elemento a elemento) se difieren a fases posteriores
del ADR 0009, justificados por casos concretos.
"""
from __future__ import annotations


def rayleigh_from_modes(
    xi1: float,
    omega1: float,
    xi2: float,
    omega2: float,
) -> tuple[float, float]:
    """Calibra ``(α, β)`` de ``C = α·M + β·K`` para reproducir las razones
    de amortiguamiento objetivo ``ξ₁, ξ₂`` en las frecuencias ``ω₁, ω₂``.

    Sistema lineal 2×2 derivado de
    ``ξ_n = (α / ω_n + β · ω_n) / 2`` evaluado en ``n = 1, 2``:

    .. math::
        \\alpha = \\frac{2\\,\\omega_1\\omega_2\\,(\\xi_1\\omega_2 - \\xi_2\\omega_1)}{\\omega_2^2 - \\omega_1^2},
        \\qquad
        \\beta  = \\frac{2\\,(\\xi_2\\omega_2 - \\xi_1\\omega_1)}{\\omega_2^2 - \\omega_1^2}.

    La elección de ``ω₁, ω₂`` determina el comportamiento fuera del rango:
    entre ellos ``ξ(ω) < max(ξ₁, ξ₂)``; fuera del rango (especialmente
    para ``ω > ω₂``) ``ξ(ω)`` crece rápido (modos altos sobreamortiguados),
    propiedad útil para suavizar respuestas espurias de alta frecuencia
    sin alterar significativamente los modos bajos relevantes.

    Parameters
    ----------
    xi1, xi2 : float
        Razones de amortiguamiento objetivo (típicamente 0.02–0.05 para
        estructura civil; 0.01 para acero estructural).
    omega1, omega2 : float
        Frecuencias en las que se imponen ``xi1`` y ``xi2``, en rad/s.
        Habitualmente las dos primeras frecuencias propias del análisis
        modal. Se requiere ``omega1 != omega2``.

    Returns
    -------
    alpha : float
        Coeficiente de masa.
    beta : float
        Coeficiente de rigidez.

    Raises
    ------
    ValueError
        Si ``omega1 == omega2`` (sistema indeterminado) o si alguna
        frecuencia es no positiva (Rayleigh no define amortiguamiento
        para modos rígidos).
    """
    if omega1 <= 0.0 or omega2 <= 0.0:
        raise ValueError(
            f"rayleigh_from_modes: ω₁={omega1}, ω₂={omega2}. Las frecuencias "
            "deben ser estrictamente positivas (Rayleigh no aplica a modos rígidos)."
        )
    if omega1 == omega2:
        raise ValueError(
            f"rayleigh_from_modes: ω₁ y ω₂ coinciden ({omega1}). El sistema "
            "es indeterminado; elige dos frecuencias distintas."
        )
    denom = omega2 * omega2 - omega1 * omega1
    alpha = 2.0 * omega1 * omega2 * (xi1 * omega2 - xi2 * omega1) / denom
    beta = 2.0 * (xi2 * omega2 - xi1 * omega1) / denom
    return alpha, beta


def resolve_rayleigh_config(
    cfg: dict | None,
    *,
    source: str = "solver",
) -> tuple[float, float]:
    """Traduce el dict de Rayleigh del usuario a la tupla ``(α, β)``.

    Tres formas aceptadas:

    - ``None`` → ``(0.0, 0.0)`` (sin amortiguamiento).
    - ``{"alpha": ..., "beta": ...}`` → coeficientes directos.
    - ``{"xi1": ..., "omega1": ..., "xi2": ..., "omega2": ...}`` →
      calibración modal automática vía :func:`rayleigh_from_modes`.

    Parameters
    ----------
    cfg
        Configuración tal como llega del constructor del solver o del
        YAML. ``None`` significa "sin amortiguamiento".
    source
        Etiqueta del llamador (``"NewmarkSolver"``, ``"HarmonicSolver"``,
        ...) que se interpola en los mensajes de error. Permite que el
        usuario localice el solver al que pertenecía la configuración
        inválida.

    Returns
    -------
    alpha, beta
        Coeficientes de Rayleigh ``C = α·M + β·K``.

    Raises
    ------
    ValueError
        Si ``cfg`` no es ``None`` ni ``dict``, o si las claves del dict
        no encajan en ninguna de las dos formas aceptadas.
    """
    if cfg is None:
        return 0.0, 0.0
    if not isinstance(cfg, dict):
        raise ValueError(
            f"{source}.rayleigh: esperado dict o None, recibido "
            f"{type(cfg).__name__}."
        )
    if "alpha" in cfg and "beta" in cfg:
        return float(cfg["alpha"]), float(cfg["beta"])
    required = {"xi1", "omega1", "xi2", "omega2"}
    if required.issubset(cfg):
        return rayleigh_from_modes(
            float(cfg["xi1"]), float(cfg["omega1"]),
            float(cfg["xi2"]), float(cfg["omega2"]),
        )
    raise ValueError(
        f"{source}.rayleigh: dict admitido es {{'alpha','beta'}} o "
        f"{{'xi1','omega1','xi2','omega2'}}; recibido {set(cfg.keys())!r}."
    )


def rayleigh_xi(alpha: float, beta: float, omega: float) -> float:
    """Razón de amortiguamiento ``ξ(ω) = (α/ω + β·ω) / 2`` predicha por
    los coeficientes Rayleigh para un modo de frecuencia ``ω``.

    Útil para diagnosticar el amortiguamiento que cada modo del análisis
    modal recibirá efectivamente, dado un par ``(α, β)`` ya calibrado.
    """
    if omega <= 0.0:
        raise ValueError(
            f"rayleigh_xi: ω={omega}. La fórmula no aplica a modos rígidos."
        )
    return 0.5 * (alpha / omega + beta * omega)
