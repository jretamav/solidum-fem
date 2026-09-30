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
   pivotes numéricamente nulos (SuperLU, Pardiso y Cholesky, mismo
   criterio): delata un mecanismo interno aunque ninguna carga lo active y
   el residuo salga perfecto.

Las capas 2 y 3 del estático **lineal** actúan tras resolver. Dentro de un
Newton o un arc-length resolver un sistema casi singular es legítimo (cerca
de un punto límite) y la red es el propio Newton, que exige equilibrio real;
ahí el corrector sólo usa el residuo lineal para diagnosticar mejor.

Hay una excepción en la que resolver un sistema singular nunca es legítimo:
la **tangente inicial** (``U = 0``) de un análisis estático no lineal
(:func:`check_initial_tangent`, 2026-09-30). Si ningún componente del modelo
puede ganar rigidez al deformarse —ningún elemento geométricamente no lineal
(``Element.GEOMETRICALLY_NONLINEAR``) ni material unilateral
(``Material.IS_UNILATERAL``)—, la no linealidad del material sólo reduce la
rigidez y un mecanismo en ``U = 0`` lo es durante todo el análisis: se
rechaza. Si alguno puede (una red de cables se rigidiza al tensarse), sólo
se avisa. La comprobación reutiliza la primera factorización del solver, que
es la de la tangente en ``U = 0``: no factoriza dos veces.

Un backend sin pivotes que inspeccionar (el iterativo del ADR 0018: un
método de Krylov no factoriza) no puede hacer la capa 3. En lugar de fingir
que la hace, :func:`warn_if_pivots_unavailable` lo dice una vez por análisis
y explica cómo comprobarlo.
"""
from __future__ import annotations

import numpy as np

from solidum.constants import EQUILIBRIUM_RTOL, MECHANISM_RANK_RTOL
from solidum.math.solvers._shared import _log
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


def zero_pivot_count(factor) -> int | None:
    """Pivotes numéricamente nulos de una factorización (``n_zero_pivots``,
    criterio ``|pivote| < ZERO_PIVOT_RTOL·max|K|``), o ``None`` si el backend
    no puede inspeccionarlos: el iterativo (ADR 0018) no factoriza la matriz.
    SuperLU, Pardiso y Cholesky los exponen con el mismo criterio."""
    n = getattr(factor, "n_zero_pivots", None)
    return None if n is None else int(n)


def warn_if_pivots_unavailable(factor, solver: str = "") -> bool:
    """Si ``factor`` no tiene pivotes que inspeccionar (backend iterativo),
    avisa de que la capa 3 no se hace y de cómo comprobarlo, y devuelve
    ``True``. Una vez por análisis: lo llama el solver estático sobre la
    factorización que un backend directo inspeccionaría.

    Distinguir un mecanismo (autovalor exactamente nulo) de un modelo sólo
    mal condicionado exigiría estimar el menor autovalor con una precisión
    relativa de ~1e-14, que un método de Krylov no da de forma fiable. En
    lugar de fingir la detección, el hueco deja de ser silencioso."""
    if zero_pivot_count(factor) is not None:
        return False
    quien = f"{solver}: " if solver else ""
    _log.warning(
        f"{quien}con el solver algebraico iterativo (linear_algebra: "
        f"iterative) no se comprueba si el modelo tiene un mecanismo interno "
        f"que la carga no activa: un método de Krylov no factoriza la matriz y "
        f"no tiene pivotes que inspeccionar. Si lo tuviera, la solución "
        f"cumpliría el equilibrio, pero sus desplazamientos incluirían una "
        f"parte arbitraria del mecanismo. Para comprobarlo, haz un análisis "
        f"modal del mismo modelo con densidad en los materiales (cada mecanismo "
        f"aparece como una frecuencia nula) o resuélvelo una vez con el solver "
        f"directo (linear_algebra: auto)."
    )
    return True


def check_linear_solution(K, u: np.ndarray, F: np.ndarray, factor=None) -> float:
    """Capas 2 y 3 del análisis estático lineal.

    Rechaza con :class:`IllPosedSystemError` una factorización con pivotes
    numéricamente nulos (si el backend lo reporta: ``n_zero_pivots``; lo
    hacen SuperLU, Pardiso y Cholesky con el mismo criterio) y una
    solución que no satisface ``‖F − K·u‖ ≤ EQUILIBRIUM_RTOL·‖F‖``. Devuelve
    el residuo relativo (0 si ``F = 0``)."""
    zero_pivots = zero_pivot_count(factor) or 0
    if zero_pivots > 0:
        raise IllPosedSystemError(zero_pivots=zero_pivots)
    F_norm = float(np.linalg.norm(F))
    if F_norm == 0.0:
        return 0.0
    rel = float(np.linalg.norm(F - K @ u)) / F_norm
    if not np.isfinite(rel) or rel > EQUILIBRIUM_RTOL:
        raise IllPosedSystemError(relative_residual=rel, rtol=EQUILIBRIUM_RTOL)
    return rel


def _stiffening_components(domain) -> list[str]:
    """Componentes del dominio cuya rigidez puede crecer respecto de la de
    ``U = 0``, descritos para el mensaje: elementos geométricamente no
    lineales (su rigidez geométrica crece con las fuerzas internas y la
    configuración) y materiales unilaterales (sin rigidez en ``U = 0``, la
    ganan al traccionarse). Lista vacía ⇒ la rigidez del modelo sólo puede
    bajar respecto de la inicial."""
    elements: set[str] = set()
    materials: set[str] = set()
    for elem in domain.elements.values():
        if getattr(elem, "GEOMETRICALLY_NONLINEAR", False):
            elements.add(type(elem).__name__)
        for material in elem.materials():
            if getattr(material, "IS_UNILATERAL", False):
                materials.add(type(material).__name__)
    found = []
    if elements:
        found.append(f"elementos geométricamente no lineales "
                     f"({', '.join(sorted(elements))})")
    if materials:
        found.append(f"materiales unilaterales ({', '.join(sorted(materials))})")
    return found


def _judge_initial_tangent(assembler, factor, solver: str) -> None:
    """Veredicto sobre la factorización de la tangente inicial (``U = 0``).

    ``factor is None`` significa que la factorización directa no pudo
    completarse por un pivote exactamente nulo (SuperLU: *Factor is exactly
    singular*)."""
    zero_pivots = None
    if factor is not None:
        if warn_if_pivots_unavailable(factor, solver):
            return
        zero_pivots = zero_pivot_count(factor)
        if zero_pivots == 0:
            return
    stiffening = _stiffening_components(assembler.domain)
    if not stiffening:
        raise IllPosedSystemError(zero_pivots=zero_pivots, solver=solver,
                                  initial_tangent=True)
    quien = f"{solver}: " if solver else ""
    pivotes = (f"tiene {zero_pivots} pivote(s) numéricamente nulo(s)" if zero_pivots
               else "no se pudo factorizar (tiene un pivote exactamente nulo)")
    _log.warning(
        f"{quien}la matriz tangente inicial (U = 0) {pivotes}: el modelo tiene "
        f"un mecanismo interno (una rótula o articulación de más, una barra o "
        f"elemento suelto, dos partes unidas por un solo nodo) o componentes "
        f"sin rigidez inicial. Como tiene {' y '.join(stiffening)}, su rigidez "
        f"puede crecer al deformarse y estabilizarlo (una red de cables sin "
        f"pretensado, por ejemplo, se rigidiza al tensarse), así que el "
        f"análisis sigue. Si no es el caso, localiza el mecanismo con un "
        f"análisis modal del mismo modelo con densidad en los materiales: cada "
        f"mecanismo aparece como una frecuencia nula y su forma modal muestra "
        f"qué parte se mueve."
    )


def check_initial_tangent(assembler, corrector, solver: str = "") -> None:
    """Pivotes nulos de la tangente inicial de un análisis estático **no
    lineal** (ADR 0019, 2026-09-30).

    Pide a ``corrector`` (:class:`~solidum.math.solvers.corrector.NewtonCorrector`)
    que le entregue su próxima factorización fresca antes de usarla. Llamada
    al empezar el análisis, esa factorización es la de la tangente en
    ``U = 0`` (la primera iteración del Newton, o el predictor del arco), así
    que la comprobación no factoriza dos veces. Sobre ella:

    - sin pivotes nulos: nada;
    - con pivotes nulos y ningún componente que pueda ganar rigidez al
      deformarse: :class:`IllPosedSystemError` (un mecanismo en ``U = 0`` lo es
      durante todo el análisis, porque la no linealidad del material sólo
      reduce la rigidez);
    - con pivotes nulos y elementos geométricamente no lineales o materiales
      unilaterales: aviso en el log (la rigidez geométrica, o la de un cable al
      tensarse, pueden estabilizarlo);
    - backend iterativo, sin pivotes: :func:`warn_if_pivots_unavailable`.

    Sólo para análisis **estáticos**: en uno dinámico o modal un modelo libre
    es legítimo. Pieza genérica: la llama cada solver estático no lineal justo
    después de :func:`ensure_statically_restrained`, incluidas las variantes
    de arco de fuera del programa principal."""
    corrector.inspect_next_factorization(
        lambda factor: _judge_initial_tangent(assembler, factor, solver))
