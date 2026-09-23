"""Red de seguridad del análisis estático (ADR 0019).

El usuario de Solidum no elige el solver algebraico, y no tiene por qué saber
qué hace un solver directo ante una matriz singular: SuperLU devuelve
desplazamientos de miles de kilómetros y MKL Pardiso perturba los pivotes
nulos y sigue adelante, ambos sin avisar (medido). Este módulo pone tres
capas entre ese comportamiento y el usuario:

1. **Antes de resolver** — :func:`ensure_statically_restrained`: los modos de
   cuerpo rígido que los apoyos no restringen son un mecanismo. Se detecta
   sin factorizar nada y se describe en términos del modelo ("traslación en
   la dirección x", "giro alrededor del eje z que pasa por (2, 0.5)").
2. **Después de resolver** — :func:`check_linear_solution`: el equilibrio
   ``‖F − K·u‖ ≤ EQUILIBRIUM_RTOL·‖F‖``. Atrapa lo que la capa 1 no ve:
   mecanismos internos cargados y matrices mal condicionadas.
3. **Pivotes nulos** — la misma función rechaza una factorización con
   pivotes numéricamente nulos (SuperLU y Pardiso, mismo criterio): delata
   un mecanismo interno aunque ninguna carga lo active y el residuo salga
   perfecto.

Las capas 2 y 3 son sólo para el análisis estático **lineal**. Dentro de un
Newton o un arc-length resolver un sistema casi singular es legítimo (cerca
de un punto límite) y la red es el propio Newton, que exige equilibrio real;
ahí el corrector sólo usa el residuo lineal para diagnosticar mejor.
"""
from __future__ import annotations

import numpy as np

from solidum.constants import EQUILIBRIUM_RTOL, MECHANISM_RANK_RTOL
from solidum.math.solvers.diagnostics import IllPosedSystemError, MechanismError

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
_AXIS_NAME = ("x", "y", "z")
# Una dirección se nombra por su eje si su componente dominante supera esto.
_ALIGNED = 1.0 - 1.0e-9


def _fmt(v: float) -> str:
    return f"{v:.4g}" if abs(v) >= 1e-12 else "0"


def _direction(v: np.ndarray, dims: int) -> str:
    u = v / np.linalg.norm(v)
    k = int(np.argmax(np.abs(u)))
    if abs(u[k]) >= _ALIGNED:
        return _AXIS_NAME[k]
    return "(" + ", ".join(_fmt(c) for c in u[:dims]) + ")"


def _describe_mechanical(t: np.ndarray, w: np.ndarray, center: np.ndarray,
                         length: float, dims: int) -> str:
    """Movimiento rígido infinitesimal ``u(x) = t + w × (x − c)`` en palabras."""
    wn = float(np.linalg.norm(w))
    tn = float(np.linalg.norm(t))
    if wn * length <= 1.0e-6 * max(tn, wn * length):
        return f"traslación en la dirección {_direction(t, dims)}"
    axis = w / wn
    # Punto del eje: p = c + (w × t)/|w|² (donde el movimiento es paralelo al eje).
    p = center + np.cross(w, t) / wn**2
    punto = "(" + ", ".join(_fmt(c) for c in p[:dims]) + ")"
    eje = f"eje {_direction(w, 3)}" if dims == 3 else "eje z"
    texto = f"giro alrededor del {eje} que pasa por {punto}"
    slide = float(t @ axis)
    if dims == 3 and abs(slide) > 1.0e-6 * max(tn, wn * length):
        texto += " (con deslizamiento simultáneo a lo largo del eje)"
    return texto


def _nullspace(S: np.ndarray, m: int) -> np.ndarray:
    """Base del núcleo de ``S`` (``m`` columnas); identidad si ``S`` no tiene filas."""
    if S.shape[0] == 0:
        return np.eye(m)
    # S tiene una fila por DOF restringido y sólo m ≤ 7 columnas. Una SVD
    # completa de S construiría una matriz densa (filas × filas): 200 MB con
    # 5 000 restricciones. Mismos vectores singulares derechos vía QR
    # económica: S = Q·R ⇒ SVD(S) y SVD(R) comparten Σ y Vᵀ.
    R = np.linalg.qr(S, mode="r")
    _, sv, Vt = np.linalg.svd(R, full_matrices=True)
    tol = MECHANISM_RANK_RTOL * max(1.0, float(sv.max()) if sv.size else 0.0)
    rank = int(np.sum(sv > tol))
    return Vt[rank:].T


def unrestrained_rigid_motions(assembler) -> list[str]:
    """Movimientos de cuerpo rígido que las restricciones del modelo dejan
    libres, descritos en palabras. Lista vacía ⇒ el modelo está apoyado.

    Un modo rígido ``r`` es compatible con las restricciones homogéneas si
    ``r = T·r_libre`` (ADR 0004: ``u = T·u_libre + g``), es decir, si se anula
    ``S·c = (B − T·B_libre)[esclavos]·c``. El núcleo de ``S`` son las
    combinaciones de modos que ningún apoyo impide. Se calcula por bloques
    (mecánico y cada campo escalar), que las restricciones no acoplan, para
    que ninguna descripción mezcle una traslación con una temperatura.
    """
    basis = assembler.rigid_basis()
    B = basis.B
    m = B.shape[1]
    if m == 0:
        return []

    cs = assembler.constraint_set
    ndof = assembler.ndof
    scale = np.max(np.abs(B), axis=0)
    scale[scale == 0.0] = 1.0
    Bn = B / scale

    slave = cs.slave_dofs
    if slave.size:
        T, _ = cs.build(ndof)
        free = cs.free_dofs(ndof)
        S = Bn[slave] - np.asarray(T[slave] @ Bn[free])
    else:
        S = np.zeros((0, m))

    labels = basis.labels
    dims = 3 if any((kind == "translation" and ax == "z") or (kind == "rotation" and ax != "z")
                    for kind, ax in labels) else 2
    mech = [j for j, (kind, _) in enumerate(labels) if kind != "scalar"]
    scalars = [j for j, (kind, _) in enumerate(labels) if kind == "scalar"]

    motions: list[str] = []
    if mech:
        N = _nullspace(S[:, mech], len(mech))
        for k in range(N.shape[1]):
            a = N[:, k] / scale[mech]  # coeficientes sobre las columnas originales
            t, w = np.zeros(3), np.zeros(3)
            for coef, j in zip(a, mech):
                kind, ax = labels[j]
                if kind == "translation":
                    t[_AXIS_INDEX[ax]] += coef
                else:  # columna de giro = (e × rel)/L
                    w[_AXIS_INDEX[ax]] += coef / basis.length
            motions.append(_describe_mechanical(t, w, basis.center, basis.length, dims))
    for j in scalars:
        if _nullspace(S[:, [j]], 1).shape[1]:
            name = labels[j][1]
            motions.append(
                f"el campo '{name}' no tiene ningún valor prescrito: su nivel "
                f"queda indeterminado (prescribe al menos un valor, p. ej. en "
                f"un borde)"
            )
    return motions


def ensure_statically_restrained(assembler, solver: str = "") -> None:
    """Capa 1: lanza :class:`MechanismError` si el modelo es un mecanismo
    rígido. Sólo tiene sentido en análisis **estáticos**: en uno modal o
    dinámico un modelo libre es legítimo."""
    motions = unrestrained_rigid_motions(assembler)
    if motions:
        raise MechanismError(motions, solver)


def check_linear_solution(K, u: np.ndarray, F: np.ndarray, factor=None) -> float:
    """Capas 2 y 3 del análisis estático lineal.

    Rechaza con :class:`IllPosedSystemError` una factorización con pivotes
    numéricamente nulos (si el backend lo reporta: ``n_zero_pivots``; lo
    hacen SuperLU y Pardiso con el mismo criterio) y una
    solución que no satisface ``‖F − K·u‖ ≤ EQUILIBRIUM_RTOL·‖F‖``. Devuelve
    el residuo relativo (0 si ``F = 0``)."""
    zero_pivots = int(getattr(factor, "n_zero_pivots", 0) or 0)
    if zero_pivots > 0:
        raise IllPosedSystemError(zero_pivots=zero_pivots)
    F_norm = float(np.linalg.norm(F))
    if F_norm == 0.0:
        return 0.0
    rel = float(np.linalg.norm(F - K @ u)) / F_norm
    if not np.isfinite(rel) or rel > EQUILIBRIUM_RTOL:
        raise IllPosedSystemError(relative_residual=rel, rtol=EQUILIBRIUM_RTOL)
    return rel
