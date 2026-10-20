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
3. **Mecanismos internos** — :func:`find_internal_mechanisms`, tras
   factorizar: delata un mecanismo interno aunque ninguna carga lo active y
   el residuo salga perfecto. Desde el 2026-10-20 no decide el tamaño de un
   pivote sino la **energía del modo** (decisión del usuario): si la
   factorización tiene pivotes pequeños frente a la diagonal (filtro barato,
   ``n_small_pivots``), se calcula con ella el modo de menor energía del
   problema escalado ``K·v = λ·diag(K)·v`` y es un mecanismo si todos los
   elementos que mueve quedan sin deformar. El pivote de un mecanismo real es
   ruido de redondeo que crece con el tamaño del modelo, y un contraste de
   rigidez legítimo también da pivotes pequeños; la energía por elemento los
   separa por muchos órdenes y dice además dónde está el mecanismo.

La capa 2 del estático **lineal** actúa tras resolver y la 3 tras
factorizar. Dentro de un Newton o un arc-length resolver un sistema casi
singular es legítimo (cerca de un punto límite) y la red es el propio
Newton, que exige equilibrio real; ahí el corrector sólo usa el residuo
lineal para diagnosticar mejor.

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

from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla

from solidum.constants import (
    EQUILIBRIUM_RTOL,
    MECHANISM_ENERGY_RTOL,
    MECHANISM_RANK_RTOL,
)
from solidum.math.solvers._shared import _log
from solidum.math.solvers.diagnostics import IllPosedSystemError, MechanismError

_AXIS_INDEX = {"x": 0, "y": 1, "z": 2}
_AXIS_NAME = ("x", "y", "z")
# Una dirección se nombra por su eje si su componente dominante supera esto.
_ALIGNED = 1.0 - 1.0e-9

# Iteración inversa en bloque del modo de menor energía (capa 3). La
# separación entre el autovalor de un mecanismo (~1e-16) y el siguiente
# (≥ 1e-3 en los modelos medidos) es tan grande que dos iteraciones bastan;
# la tercera es margen.
_MODE_ITERATIONS = 3
# Un elemento no se mueve en el modo si su "masa" v_eᵀ·D_e·v_e no llega a
# esta fracción de la del elemento que más se mueve (amplitud relativa 1e-6):
# por debajo, su deformación es ruido de redondeo de la iteración.
_MOVED_RTOL = 1.0e-12
# Sólo se examinan por elementos los modos con λ < esto·MECHANISM_ENERGY_RTOL.
# Si todos los elementos que mueve un modo tienen r_e < MECHANISM_ENERGY_RTOL,
# su λ queda por debajo de ese umbral por el número de elementos que comparten
# un nodo (≤ 27 en mallas hexaédricas estructuradas): el factor es holgado.
_CANDIDATE_FACTOR = 1.0e3
# El movimiento de una parte se describe como rígido si el ajuste por
# mínimos cuadrados deja un residuo relativo por debajo de esto.
_RIGID_FIT_RTOL = 1.0e-6
_TRANSLATIONS = ("ux", "uy", "uz")


@dataclass(frozen=True)
class InternalMechanism:
    """Mecanismo interno localizado (capa 3, ampliación del 2026-10-20).

    Attributes
    ----------
    elements : tuple[int, ...]
        Ids de los elementos que el modo mueve sin deformarlos. Vacío si el
        mecanismo no se pudo localizar.
    description : str
        El mecanismo en palabras ("la parte formada por los elementos 65–128
        puede moverse sin deformarse: giro alrededor del eje z que pasa por
        (1, 1)").
    energy : float
        Mayor energía relativa ``v_eᵀ·K_e·v_e / v_eᵀ·D_e·v_e`` de los elementos
        que mueve (``< MECHANISM_ENERGY_RTOL``).
    dimension : int
        Movimientos independientes de esa parte (una parte unida al resto por
        un solo nodo en 3D puede girar alrededor de tres ejes).
    """

    elements: tuple[int, ...]
    description: str
    energy: float
    dimension: int = 1


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


def _column_space(M: np.ndarray) -> np.ndarray:
    """Base ortonormal del espacio columna de ``M`` (pocas filas y columnas)."""
    if M.size == 0:
        return M
    U, sv, _ = np.linalg.svd(M, full_matrices=False)
    tol = MECHANISM_RANK_RTOL * max(1.0, float(sv.max()) if sv.size else 0.0)
    return U[:, : int(np.sum(sv > tol))]


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
    # Un modo que sólo mueve los DOFs en el ruido de redondeo (el giro
    # alrededor de la recta de unos nodos alineados fuera del origen) es
    # nulo: normalizado por su máximo se amplificaría a un movimiento ficticio.
    negligible = scale <= MECHANISM_RANK_RTOL * max(1.0, float(scale.max()))
    if negligible.any():
        B = np.array(B, dtype=float, copy=True)
        B[:, negligible] = 0.0
    scale[negligible] = 1.0
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
        # Una combinación de modos que no mueve ningún DOF (B·c = 0) cumple
        # S·c = 0 sin que ningún apoyo la impida, pero no es un movimiento
        # del modelo: el giro alrededor de la recta que contiene todos los
        # nodos de un modelo de barras articuladas, o un giro alrededor de un
        # nodo aislado sin DOFs de giro. Se quita del núcleo antes de
        # describirlo; si B tiene rango completo, N no cambia.
        Z = _nullspace(Bn[:, mech], len(mech))
        if Z.shape[1]:
            N = _column_space(N - Z @ (Z.T @ N))
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


def small_pivot_count(factor) -> int | None:
    """Pivotes pequeños frente a la diagonal de una factorización
    (``n_small_pivots``, el filtro de la capa 3), o ``None`` si el backend no
    puede inspeccionarlos: el iterativo (ADR 0018) no factoriza la matriz.
    SuperLU, Pardiso y Cholesky los exponen."""
    n = getattr(factor, "n_small_pivots", None)
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
    if small_pivot_count(factor) is not None:
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


# ----------------------------------------------------------------------
# Capa 3: mecanismos internos por la energía del modo (2026-10-20)
# ----------------------------------------------------------------------

def _id_ranges(ids, limit: int = 12) -> str:
    """``[1, 2, 3, 7, 9, 10]`` → ``"1–3, 7, 9, 10"`` (como mucho ``limit`` grupos)."""
    ids = sorted(int(i) for i in ids)
    groups: list[tuple[int, int]] = []
    for i in ids:
        if groups and i == groups[-1][1] + 1:
            groups[-1] = (groups[-1][0], i)
        else:
            groups.append((i, i))
    parts = [str(a) if a == b else (f"{a}, {b}" if b == a + 1 else f"{a}–{b}")
             for a, b in groups]
    if len(parts) > limit:
        parts = parts[:limit] + ["…"]
    return ", ".join(parts)


def _element_blocks(assembler):
    """``(elemento, índices globales, K_e en U = 0)`` de todo el dominio y la
    diagonal de la ``K`` ensamblada (``Σ_e |diag K_e|``, sin reducir).

    Una pasada por elementos con :meth:`Element.compute_global_stiffness`, el
    coste de un ensamblaje sin lotes; sólo se paga cuando hay un modo
    candidato a mecanismo."""
    D = np.zeros(assembler.ndof)
    blocks = []
    for elem in assembler.domain.elements.values():
        idx = np.asarray(elem.get_global_dof_indices(), dtype=np.intp)
        K_e = np.asarray(elem.compute_global_stiffness(), dtype=float)
        np.add.at(D, idx, np.abs(np.diag(K_e)))
        blocks.append((elem, idx, K_e))
    return blocks, D


def _rigid_motion(nodes, v: np.ndarray):
    """Ajuste por mínimos cuadrados de ``u(x) = t + w × (x − c)`` a las
    traslaciones de ``nodes``. Devuelve ``(t, w, c, L, dims)`` si el
    movimiento es rígido (residuo relativo < ``_RIGID_FIT_RTOL``) o ``None``."""
    dims = 3 if any("uz" in n.dofs for n in nodes) else 2
    coords = np.array([np.pad(np.asarray(n.coordinates, dtype=float)[:3],
                              (0, 3 - min(3, len(n.coordinates)))) for n in nodes])
    c = coords.mean(axis=0)
    rows, rhs = [], []
    for n, x in zip(nodes, coords):
        r = x - c
        # Fila k de u = t + w × r con incógnitas (tx, ty, tz, wx, wy, wz).
        cross = ((0.0, r[2], -r[1]), (-r[2], 0.0, r[0]), (r[1], -r[0], 0.0))
        for k, d in enumerate(_TRANSLATIONS[:dims]):
            if d in n.dofs:
                row = np.zeros(6)
                row[k] = 1.0
                row[3:] = cross[k]
                rows.append(row)
                rhs.append(v[n.dofs[d]])
    A, b = np.array(rows), np.array(rhs)
    b_norm = float(np.linalg.norm(b))
    if b_norm == 0.0:
        return None
    cols = [0, 1, 2, 3, 4, 5] if dims == 3 else [0, 1, 5]
    s, *_ = np.linalg.lstsq(A[:, cols], b, rcond=None)
    if np.linalg.norm(A[:, cols] @ s - b) > _RIGID_FIT_RTOL * b_norm:
        return None
    full = np.zeros(6)
    full[cols] = s
    length = float(np.linalg.norm(coords.max(axis=0) - coords.min(axis=0))) or 1.0
    return full[:3], full[3:], c, length, dims


def _canonical_rigid(coefs: np.ndarray, length: float) -> list[tuple[np.ndarray, np.ndarray]]:
    """Base legible del subespacio de movimientos rígidos ``coefs`` (una fila
    ``[t, w]`` por modo): primero las traslaciones puras (direcciones
    ortonormales), después los giros (ejes ortonormales), cada giro sin la
    parte de traslación que ya está libre. Sin esto, dos movimientos libres
    de una misma parte (deslizar y girar) saldrían como dos giros alrededor
    de puntos arbitrarios, según la combinación que dé la iteración."""
    scaled = np.hstack([coefs[:, :3], coefs[:, 3:] * length])
    tol = 1.0e-6 * max(float(np.abs(scaled).max()), 1e-300)
    _, sw, Vt = np.linalg.svd(scaled[:, 3:].T)     # Wᵀ·a = 0 ⇒ traslación pura
    rank_w = int(np.sum(sw > tol))
    trans = (Vt[rank_w:] @ coefs)[:, :3]
    rots = Vt[:rank_w] @ coefs
    motions: list[tuple[np.ndarray, np.ndarray]] = []
    t_dirs = np.zeros((0, 3))
    if trans.shape[0]:
        _, st, Vt_t = np.linalg.svd(trans)
        t_dirs = Vt_t[: int(np.sum(st > tol))]
        motions += [(d, np.zeros(3)) for d in t_dirs]
    if rots.shape[0]:
        _, _, Vt_w = np.linalg.svd(rots[:, 3:])
        # Combinaciones de los giros con ejes ortonormales (ejes principales).
        rots = np.linalg.lstsq(rots[:, 3:].T, Vt_w[:rank_w].T, rcond=None)[0].T @ rots
        for t, w in zip(rots[:, :3], rots[:, 3:]):
            motions.append((t - t_dirs.T @ (t_dirs @ t), w))
    return motions


def _common_pivot(rotations, c: np.ndarray, length: float) -> np.ndarray | None:
    """Punto ``P`` que los tres giros ``u = t + w × (x − c)`` dejan fijo, si
    existe: una parte que puede girar alrededor de cualquier eje que pasa
    por ``P`` (unida al resto por un solo nodo). ``u(P) = 0`` es lineal en
    ``P``: ``[w]ₓ·(P − c) = −t`` para los tres giros."""
    def cross_matrix(w):
        return np.array([[0.0, -w[2], w[1]], [w[2], 0.0, -w[0]], [-w[1], w[0], 0.0]])
    A = np.vstack([cross_matrix(w) for _, w in rotations])
    b = -np.concatenate([t for t, _ in rotations])
    d, *_ = np.linalg.lstsq(A, b, rcond=None)
    scale = max(float(np.linalg.norm(b)), max(float(np.linalg.norm(w)) for _, w in rotations) * length)
    if np.linalg.norm(A @ d - b) > _RIGID_FIT_RTOL * scale:
        return None
    return c + d


def _describe_part(elements, modes: list[np.ndarray], D: np.ndarray) -> str:
    """Los modos ``modes`` (vectores completos) que mueven, sin deformarlos,
    los mismos ``elements``, en palabras: movimientos rígidos de esa parte,
    un movimiento sin deformación que no es rígido (articulaciones de más) o,
    si sólo mueven un campo escalar, una parte cuyo nivel de ese campo queda
    indeterminado."""
    ids = _id_ranges(e.id for e in elements)
    nodes = list({n.id: n for e in elements for n in e.nodes}.values())
    k = len(modes)
    formas = f" de {k} formas independientes" if k > 1 else ""
    share: dict[str, float] = {}
    for v in modes:
        for n in nodes:
            for name, eq in n.dofs.items():
                share[name] = share.get(name, 0.0) + float(D[eq] * v[eq] ** 2)
    total = sum(share.values()) or 1.0
    moving = {name for name, s in share.items() if s > 1.0e-6 * total}
    mechanical = {name for name in moving if name in _TRANSLATIONS or name in ("rx", "ry", "rz")}
    if not mechanical:
        campos = ", ".join(f"'{name}'" for name in sorted(moving))
        return (f"la parte formada por los elementos {ids} no está unida a ningún "
                f"valor prescrito del campo {campos}: su nivel queda indeterminado")
    if not mechanical & set(_TRANSLATIONS):
        giros = _id_ranges(n.id for n in nodes)
        return (f"los giros de los nodos {giros} de los elementos {ids} pueden "
                f"cambiar{formas} sin deformar ningún elemento")
    moved_nodes = [n for n in nodes if any(d in n.dofs for d in _TRANSLATIONS)]
    fits = [_rigid_motion(moved_nodes, v) for v in modes]
    if any(fit is None for fit in fits):
        return (f"los elementos {ids} pueden moverse{formas} sin deformarse, "
                f"aunque no como un sólido rígido: hay articulaciones de más "
                f"entre ellos")
    _, _, c, length, dims = fits[0]
    coefs = np.array([np.concatenate([t, w]) for t, w, *_ in fits])
    motions = _canonical_rigid(coefs, length)
    rotations = [(t, w) for t, w in motions if np.any(w)]
    pivot = _common_pivot(rotations, c, length) if dims == 3 and len(rotations) == 3 else None
    if pivot is not None:
        punto = "(" + ", ".join(_fmt(x) for x in pivot) + ")"
        textos = [_describe_mechanical(t, w, c, length, dims) for t, w in motions if not np.any(w)]
        textos.append(f"giro alrededor de cualquier eje que pasa por {punto} (como una "
                      f"rótula esférica)")
    else:
        textos = [_describe_mechanical(t, w, c, length, dims) for t, w in motions]
    return (f"la parte formada por los elementos {ids} puede moverse sin "
            f"deformarse{formas}: {'; '.join(textos)}")


def _moved_without_strain(blocks, D: np.ndarray, v: np.ndarray):
    """Índices (en ``blocks``) de los elementos que el modo ``v`` (vector
    completo) mueve y su mayor energía relativa, si no deforma ninguno de
    ellos; ``None`` si deforma alguno.

    Energía relativa por elemento ``r_e = v_eᵀ·K_e·v_e / v_eᵀ·D_e·v_e`` con
    ``D_e`` la diagonal ensamblada en sus DOF: no depende de las unidades
    (cada fila y columna de ``K`` escala como la entrada de ``D``) y no se
    anula para una barra que se mueve en transversal, cuya ``K_e`` sí tiene
    diagonal nula en esa dirección."""
    m = np.array([v[idx] @ (D[idx] * v[idx]) for _, idx, _ in blocks])
    if m.size == 0 or not m.max() > 0.0:
        return None
    moved = np.flatnonzero(m > _MOVED_RTOL * m.max())
    energy = 0.0
    for i in moved:
        _, idx, K_e = blocks[i]
        ve = v[idx]
        energy = max(energy, float(ve @ (K_e @ ve)) / m[i])
        if energy >= MECHANISM_ENERGY_RTOL:
            return None
    return tuple(int(i) for i in moved), energy


def find_internal_mechanisms(assembler, K_red, factor) -> list[InternalMechanism] | None:
    """Capa 3 (ampliación del 2026-10-20, decisión del usuario): mecanismos
    internos del sistema reducido ``K_red``, ya factorizado en ``factor``.

    1. **Filtro**: sin pivotes pequeños frente a la diagonal
       (``n_small_pivots == 0``) no hay nada que buscar y no cuesta nada. El
       backend iterativo (``None``) tampoco busca: no tiene factorización.
    2. **Modo de menor energía** del problema escalado
       ``K·v = λ·diag(K)·v``: iteración inversa en bloque con la
       factorización existente (``_MODE_ITERATIONS`` resoluciones por vector)
       y Rayleigh-Ritz. ``λ`` no depende de las unidades.
    3. **Veredicto por elementos**: un modo con ``λ`` pequeño es un mecanismo
       si todos los elementos que mueve tienen energía relativa
       ``< MECHANISM_ENERGY_RTOL``. Un contraste de rigidez grande (zonas
       rígidas, ν cerca de 1/2) da ``λ`` pequeño pero deforma algún
       elemento.

    Los modos que mueven los mismos elementos se agrupan en un mecanismo de
    ``dimension`` movimientos independientes. Devuelve los mecanismos (lista
    vacía si no hay) o ``None`` si hay pivotes pequeños pero el modo no pudo
    calcularse (solución no finita): el llamador lo trata como singular,
    igual que antes de esta ampliación.
    """
    n_small = small_pivot_count(factor)
    if not n_small:
        return []
    n = K_red.shape[0]
    D = np.abs(np.asarray(K_red.diagonal(), dtype=float))
    D = np.where(D > 0.0, D, float(D.max()) if D.size and D.max() > 0.0 else 1.0)
    sD = np.sqrt(D)
    p = min(n, n_small + 2)
    X = np.random.default_rng(0).standard_normal((n, p))
    try:
        for _ in range(_MODE_ITERATIONS):
            X = np.column_stack([factor.solve(D * X[:, j]) for j in range(p)])
            if not np.all(np.isfinite(X)):
                return None
            Q, _ = np.linalg.qr(sD[:, None] * X)   # Xᵀ·D·X = I
            X = Q / sD[:, None]
    except RuntimeError:
        return None
    A = X.T @ np.asarray(K_red @ X)
    lam, V = sla.eigh(0.5 * (A + A.T))
    candidates = [V[:, k] for k in range(p) if lam[k] < _CANDIDATE_FACTOR * MECHANISM_ENERGY_RTOL]
    if not candidates:
        return []
    T, _ = assembler.constraint_set.build(assembler.ndof)
    blocks, D_full = _element_blocks(assembler)
    groups: dict[tuple[int, ...], tuple[list[np.ndarray], float]] = {}
    for y in candidates:
        v = np.asarray(T @ (X @ y)).ravel()
        verdict = _moved_without_strain(blocks, D_full, v)
        if verdict is not None:
            moved, energy = verdict
            modes, worst = groups.get(moved, ([], 0.0))
            groups[moved] = (modes + [v], max(worst, energy))
    found = []
    for moved, (modes, energy) in groups.items():
        elements = [blocks[i][0] for i in moved]
        found.append(InternalMechanism(tuple(e.id for e in elements),
                                       _describe_part(elements, modes, D_full),
                                       energy, len(modes)))
    return found


def ensure_no_internal_mechanisms(assembler, K_red, factor, solver: str = "") -> None:
    """Capa 3 del estático lineal: :class:`IllPosedSystemError` si el sistema
    reducido tiene un mecanismo interno (ver :func:`find_internal_mechanisms`),
    con cada uno descrito y localizado. Se llama una vez por factorización,
    antes de resolver."""
    found = find_internal_mechanisms(assembler, K_red, factor)
    if found is None:
        raise IllPosedSystemError(zero_pivots=small_pivot_count(factor), solver=solver)
    if found:
        raise IllPosedSystemError(zero_pivots=sum(m.dimension for m in found), solver=solver,
                                  mechanisms=[m.description for m in found])
    if small_pivot_count(factor):
        _log.info("  Pivotes pequeños frente a la diagonal sin mecanismo interno: "
                  "el modelo tiene rigideces muy desproporcionadas, pero ningún "
                  "modo se mueve sin deformar algún elemento.")


def check_linear_solution(K, u: np.ndarray, F: np.ndarray) -> float:
    """Capa 2 del análisis estático lineal.

    Rechaza con :class:`IllPosedSystemError` una solución que no satisface
    ``‖F − K·u‖ ≤ EQUILIBRIUM_RTOL·‖F‖``. Devuelve el residuo relativo (0 si
    ``F = 0``). Los mecanismos internos que la carga no activa los rechaza
    antes :func:`ensure_no_internal_mechanisms` (capa 3)."""
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


def _judge_initial_tangent(assembler, factor, A, solver: str) -> None:
    """Veredicto sobre la factorización de la tangente inicial (``U = 0``),
    ``A`` la tangente reducida que se factorizó.

    ``factor is None`` significa que la factorización directa no pudo
    completarse por un pivote exactamente nulo (SuperLU: *Factor is exactly
    singular*)."""
    zero_pivots = None
    mechanisms = None
    if factor is not None:
        if warn_if_pivots_unavailable(factor, solver):
            return
        found = find_internal_mechanisms(assembler, A, factor)
        if found == []:
            return
        if found is None:
            zero_pivots = small_pivot_count(factor)
        else:
            zero_pivots = sum(m.dimension for m in found)
            mechanisms = [m.description for m in found]
    stiffening = _stiffening_components(assembler.domain)
    if not stiffening:
        raise IllPosedSystemError(zero_pivots=zero_pivots, solver=solver,
                                  initial_tangent=True, mechanisms=mechanisms)
    quien = f"{solver}: " if solver else ""
    if mechanisms:
        donde = "; ".join(mechanisms)
        causa = (f"es singular: el modelo tiene {zero_pivots} mecanismo(s) "
                 f"interno(s) ({donde}) o componentes sin rigidez inicial")
    elif zero_pivots:
        causa = (f"tiene {zero_pivots} pivote(s) pequeño(s) y no se pudo calcular "
                 f"su modo: el modelo tiene un mecanismo interno (una rótula o "
                 f"articulación de más, una barra o elemento suelto, dos partes "
                 f"unidas por un solo nodo) o componentes sin rigidez inicial")
    else:
        causa = ("no se pudo factorizar (tiene un pivote exactamente nulo): el "
                 "modelo tiene un mecanismo interno (una rótula o articulación "
                 "de más, una barra o elemento suelto, dos partes unidas por un "
                 "solo nodo) o componentes sin rigidez inicial")
    _log.warning(
        f"{quien}la matriz tangente inicial (U = 0) {causa}. Como tiene "
        f"{' y '.join(stiffening)}, su rigidez puede crecer al deformarse y "
        f"estabilizarlo (una red de cables sin pretensado, por ejemplo, se "
        f"rigidiza al tensarse), así que el análisis sigue. Si no es el caso, "
        f"revisa la conectividad de esos elementos; un análisis modal del mismo "
        f"modelo con densidad en los materiales también lo muestra: cada "
        f"mecanismo aparece como una frecuencia nula y su forma modal muestra "
        f"qué parte se mueve."
    )


def check_initial_tangent(assembler, corrector, solver: str = "") -> None:
    """Mecanismos internos en la tangente inicial de un análisis estático
    **no lineal** (ADR 0019, 2026-09-30; por la energía del modo desde el
    2026-10-20).

    Pide a ``corrector`` (:class:`~solidum.math.solvers.corrector.NewtonCorrector`)
    que le entregue su próxima factorización fresca antes de usarla. Llamada
    al empezar el análisis, esa factorización es la de la tangente en
    ``U = 0`` (la primera iteración del Newton, o el predictor del arco), así
    que la comprobación no factoriza dos veces. Sobre ella
    (:func:`find_internal_mechanisms`):

    - sin mecanismo interno: nada;
    - con un mecanismo y ningún componente que pueda ganar rigidez al
      deformarse: :class:`IllPosedSystemError` (un mecanismo en ``U = 0`` lo es
      durante todo el análisis, porque la no linealidad del material sólo
      reduce la rigidez);
    - con un mecanismo y elementos geométricamente no lineales o materiales
      unilaterales: aviso en el log (la rigidez geométrica, o la de un cable al
      tensarse, pueden estabilizarlo);
    - backend iterativo, sin factorización: :func:`warn_if_pivots_unavailable`.

    Sólo para análisis **estáticos**: en uno dinámico o modal un modelo libre
    es legítimo. Pieza genérica: la llama cada solver estático no lineal justo
    después de :func:`ensure_statically_restrained`, incluidas las variantes
    de arco de fuera del programa principal."""
    corrector.inspect_next_factorization(
        lambda factor, A: _judge_initial_tangent(assembler, factor, A, solver))
