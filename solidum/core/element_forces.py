"""Tipo `ElementForces` — fuerzas internas (N, V, M, T) de un elemento estructural 1D en ejes locales.

Vive en ``solidum.core`` porque es parte del contrato de :class:`Element`
(ADR 0002): cada subclase de elemento **estructural 1D** (truss, cable,
frame 2D/3D) devuelve un :class:`ElementForces` desde
``internal_forces(U)``. Mantenerlo aquí en lugar de en ``solidum.results``
evita que ``core/element.py`` importe la capa de resultados (capa más
alta), preservando la jerarquía conceptual ``core → ... → results``.

Nota terminológica: ``N``, ``V``, ``M``, ``T`` son **fuerzas internas**
(resultantes integradas del tensor de esfuerzos $\\sigma$ sobre la
sección), no esfuerzos. El término "esfuerzo" se reserva para el campo
tensorial $\\sigma(x)$ por unidad de área.

**Dominio de aplicación (ADR 0012)**: este contrato aplica **solo a
elementos estructurales 1D**. Los sólidos (2D y 3D) **no** lo exponen:
en un continuo el equivalente de "fuerza interna seccional" es el campo
tensorial ``σ(x)``, accesible vía
:meth:`Element.compute_gauss_state(U)` (devuelve ``{stress, strain,
points_natural, points_global}`` por punto de integración). Convertir
``σ(x)`` a un valor único por elemento (promedio, centroide, valor
representativo) introduce ambigüedad — cada consumidor escoge la
agregación que necesite a partir de ``compute_gauss_state``.

``solidum.results`` lo re-exporta como parte de la API pública, así
``from solidum.results import ElementForces`` sigue funcionando para
consumidores externos.

Convenciones de signo: ``Reglas.md §5``.
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from typing import Literal

import numpy as np


ElementKind = Literal["truss", "cable", "frame2d", "frame3d"]

# Claves válidas por familia. Un componente ausente del dict significa "no aplica"
# a ese tipo de elemento (p. ej. no hay M en un truss) — no se rellena con ceros.
_VALID_COMPONENTS: dict[ElementKind, frozenset[str]] = {
    "truss":   frozenset({"N"}),
    "cable":   frozenset({"N"}),
    "frame2d": frozenset({"N", "V", "M"}),
    "frame3d": frozenset({"N", "Vy", "Vz", "T", "My", "Mz"}),
}


@dataclass(frozen=True)
class ElementForces:
    """Fuerzas internas en ejes locales, evaluadas en los nodos i, j.

    Inmutable (shallow): ``frozen=True`` impide reasignar atributos; los
    arrays NumPy interiores siguen siendo mutables por contenido. Tratar
    como read-only por convención.

    Parameters
    ----------
    kind
        Familia del elemento. Determina qué claves pueden aparecer en ``components``.
    components
        Diccionario ``{nombre: array(2,)}`` con el valor en el nodo i (índice 0) y
        en el nodo j (índice 1). Las claves válidas dependen de ``kind``:

        - ``"truss"`` / ``"cable"``: ``{"N"}``.
        - ``"frame2d"``: ``{"N", "V", "M"}``. Convención de viga estructural
          (Reglas.md §5): ``V`` rota el diferencial en horario si es positivo,
          ``M`` positivo es sagging.
        - ``"frame3d"``: ``{"N", "Vy", "Vz", "T", "My", "Mz"}``. Convención
          stress-resultant / RHR pura (Reglas.md §5).
    """

    kind: ElementKind
    components: dict[str, np.ndarray]

    def __post_init__(self) -> None:
        valid = _VALID_COMPONENTS.get(self.kind)
        if valid is None:
            raise ValueError(f"ElementForces.kind desconocido: {self.kind!r}")

        extra = set(self.components) - valid
        if extra:
            raise ValueError(
                f"Componentes inválidos para kind={self.kind!r}: {sorted(extra)}. "
                f"Válidos: {sorted(valid)}"
            )

        for name, arr in self.components.items():
            if not isinstance(arr, np.ndarray) or arr.shape != (2,):
                raise ValueError(
                    f"components[{name!r}] debe ser np.ndarray shape (2,), "
                    f"recibido {type(arr).__name__} shape={getattr(arr, 'shape', None)}"
                )

    def at_node_i(self) -> dict[str, float]:
        """Valores en el nodo i como escalares."""
        return {k: float(v[0]) for k, v in self.components.items()}

    def at_node_j(self) -> dict[str, float]:
        """Valores en el nodo j como escalares."""
        return {k: float(v[1]) for k, v in self.components.items()}


# ---------------------------------------------------------------------------
# Fuerzas internas a lo largo del elemento (ADR 0044 §2)
# ---------------------------------------------------------------------------

#: Por familia, cada componente en ``ξ`` a partir de su valor en el nodo i
#: (convención pública) más términos de la carga: ``("Q", k, signo)`` suma
#: ``signo·Q_k(ξ)``, la integral de la componente local ``k`` de la carga;
#: ``("R", k, signo)``, ``signo·R_k(ξ)``, su momento respecto del corte; y
#: ``("s", c, signo)``, ``signo·s·c_i`` con ``s = ξ·l`` el brazo sobre la cuerda.
#: Specs Frame2DEuler §15, Frame3D §14, Truss2D §12 (verificadas por el
#: teórico el 2026-11-27).
_ALONG_TERMS: dict[ElementKind, dict[str, tuple]] = {
    "truss": {"N": (("Q", 0, -1.0),)},
    "frame2d": {"N": (("Q", 0, -1.0),),
                "V": (("Q", 1, 1.0),),
                "M": (("s", "V", 1.0), ("R", 1, 1.0))},
    "frame3d": {"N": (("Q", 0, -1.0),),
                "Vy": (("Q", 1, -1.0),),
                "Vz": (("Q", 2, -1.0),),
                "T": (),
                "My": (("s", "Vz", 1.0), ("R", 2, -1.0)),
                "Mz": (("s", "Vy", -1.0), ("R", 1, 1.0))},
}
#: Componentes locales de la carga que usa cada familia.
ALONG_LOAD_AXES: dict[ElementKind, int] = {"truss": 1, "frame2d": 2, "frame3d": 3}
#: Cortante cuyo cero da el extremo de cada flector (``dM/ds ∝ V``).
SHEAR_OF_MOMENT: dict[ElementKind, dict[str, str]] = {
    "truss": {}, "frame2d": {"M": "V"}, "frame3d": {"Mz": "Vy", "My": "Vz"}}


@dataclass(frozen=True)
class ElementForcesAlong:
    """Fuerzas internas en los puntos ``xi ∈ [0, 1]`` de un elemento de línea
    (``xi = 0`` en el nodo i), en ejes locales y la convención de
    :class:`ElementForces`. ``s = xi·l``: distancia sobre la cuerda (``l`` la
    longitud corriente, la de referencia en los elementos lineales)."""

    kind: ElementKind
    xi: np.ndarray
    s: np.ndarray
    components: dict[str, np.ndarray]


def load_polynomials(points, samples_local) -> np.ndarray:
    """Coeficientes ``(..., n, k)`` en potencias de ``ξ`` del polinomio que
    interpola las muestras ``(..., n, k)`` de la carga en las abscisas
    ``points`` (``n``): la carga representada ``q̃`` (spec Frame2DEuler §14)."""
    V_inv = _vandermonde_inverse(tuple(float(p) for p in np.asarray(points).ravel()))
    return np.einsum("pn,...nk->...pk", V_inv, np.asarray(samples_local, dtype=np.float64))


@functools.lru_cache(maxsize=16)
def _vandermonde_inverse(points: tuple) -> np.ndarray:
    """Inversa de la matriz de Vandermonde de las abscisas (en caché: son las
    mismas de todos los elementos de una clase)."""
    V_inv = np.linalg.inv(np.vander(np.asarray(points), len(points), increasing=True))
    V_inv.setflags(write=False)
    return V_inv


def forces_along_batch(kind: ElementKind, ends_i: dict, coeffs: np.ndarray, L0: np.ndarray,
                       l: np.ndarray, xi: np.ndarray) -> dict[str, np.ndarray]:
    """Fuerzas internas ``(E, m)`` de ``E`` elementos de una familia en las
    ``m`` abscisas ``xi``, por estática del tramo (ADR 0044 §2).

    ``ends_i``: ``{componente: (E,)}`` en el nodo i; ``coeffs``: ``(E, n, k)``
    de :func:`load_polynomials` con la carga en los ejes de las fuerzas
    internas; ``L0``, ``l``: longitudes de referencia y corriente ``(E,)``.
    La carga se integra por unidad de longitud de referencia,
    ``Q(ξ) = L0·∫₀^ξ q̃``, y su momento con los brazos sobre la cuerda,
    ``R(ξ) = L0·l·∫₀^ξ (ξ − ξ')·q̃``."""
    xi = np.asarray(xi, dtype=np.float64)
    coeffs = np.asarray(coeffs, dtype=np.float64)
    L0 = np.asarray(L0, dtype=np.float64)[:, None]
    l = np.asarray(l, dtype=np.float64)[:, None]
    n = coeffs.shape[1]
    p = np.arange(n, dtype=np.float64)
    # ∫₀^ξ ξ'^p = ξ^(p+1)/(p+1);  ∫₀^ξ (ξ − ξ')ξ'^p = ξ^(p+2)/((p+1)(p+2)).
    q_basis = xi[:, None] ** (p + 1.0) / (p + 1.0)                  # (m, n)
    r_basis = xi[:, None] ** (p + 2.0) / ((p + 1.0) * (p + 2.0))    # (m, n)
    Q = L0[:, :, None] * np.einsum("mn,enk->emk", q_basis, coeffs)  # (E, m, k)
    R = (L0 * l)[:, :, None] * np.einsum("mn,enk->emk", r_basis, coeffs)
    s = l * xi[None, :]
    out = {}
    for comp, terms in _ALONG_TERMS[kind].items():
        value = np.repeat(np.asarray(ends_i[comp], dtype=np.float64)[:, None], len(xi), axis=1)
        for what, key, sign in terms:
            if what == "Q":
                value = value + sign * Q[:, :, key]
            elif what == "R":
                value = value + sign * R[:, :, key]
            else:
                value = value + sign * s * np.asarray(ends_i[key], dtype=np.float64)[:, None]
        out[comp] = value
    return out


def load_moment_integral(points, samples_local, xi) -> np.ndarray:
    """``(m, k)``: ``∫₀^ξ (ξ − ξ')·q̃(ξ') dξ'`` en las abscisas ``xi`` para cada
    componente local de la carga representada ``q̃`` (las muestras ``(n, k)``
    en ``points``), sin los factores de longitud. Lo usa el término
    geométrico de los diagramas de una perturbación (spec Frame2DEulerCorot
    §11)."""
    xi = np.asarray(xi, dtype=np.float64)
    coeffs = load_polynomials(points, samples_local)                 # (n, k)
    p = np.arange(coeffs.shape[0], dtype=np.float64)
    r_basis = xi[:, None] ** (p + 2.0) / ((p + 1.0) * (p + 2.0))    # (m, n)
    return r_basis @ coeffs


def shear_zeros(kind: ElementKind, ends_i: dict, coeffs: np.ndarray, L0: float) -> dict[str, np.ndarray]:
    """Abscisas ``ξ ∈ (0, 1)`` donde se anula cada cortante de un elemento
    (donde su flector es extremo): ``{flector: ξ}``. ``coeffs`` ``(n, k)``."""
    out = {}
    for moment, shear in SHEAR_OF_MOMENT[kind].items():
        (_, key, sign), = _ALONG_TERMS[kind][shear]
        c = np.asarray(coeffs, dtype=np.float64)[:, key]
        # V(ξ) = V_i + sign·L0·Σ c_p ξ^(p+1)/(p+1).
        poly = np.concatenate([[float(ends_i[shear])],
                               sign * L0 * c / np.arange(1, len(c) + 1)])
        scale = np.max(np.abs(poly))
        if scale == 0.0:
            out[moment] = np.empty(0)
            continue
        trimmed = np.trim_zeros(np.where(np.abs(poly) > 1e-14 * scale, poly, 0.0), "b")
        if len(trimmed) < 2:
            out[moment] = np.empty(0)
            continue
        roots = np.polynomial.polynomial.polyroots(trimmed)
        real = roots[np.abs(roots.imag) <= 1e-10].real
        out[moment] = np.sort(real[(real > 1e-12) & (real < 1.0 - 1e-12)])
    return out


def forces_along(forces: "ElementForces", points, samples_local, L0: float, l: float,
                 xi) -> ElementForcesAlong:
    """Fuerzas internas de un elemento en ``xi`` a partir de sus fuerzas de
    extremo ``forces`` (las del nodo i), las muestras de la carga en los ejes
    de las fuerzas internas (``(n, k)`` en ``points``; ``None`` sin carga) y
    las longitudes ``L0``, ``l``."""
    kind = forces.kind
    if kind not in _ALONG_TERMS:
        raise ValueError(f"forces_along: la familia {kind!r} no da fuerzas internas a lo "
                         "largo del elemento.")
    xi = np.atleast_1d(np.asarray(xi, dtype=np.float64))
    if np.any((xi < 0.0) | (xi > 1.0)):
        raise ValueError("forces_along: xi debe estar en [0, 1] (fracción de la longitud "
                         "desde el nodo i).")
    k = ALONG_LOAD_AXES[kind]
    n = len(points)
    if samples_local is None:
        coeffs = np.zeros((1, n, k))
    else:
        coeffs = load_polynomials(points, np.asarray(samples_local)[:, :k])[None]
    ends_i = {c: np.array([v[0]]) for c, v in forces.components.items()}
    comps = forces_along_batch(kind, ends_i, coeffs, np.array([L0]), np.array([l]), xi)
    return ElementForcesAlong(kind=kind, xi=xi, s=float(l) * xi,
                              components={c: v[0] for c, v in comps.items()})


def axial_end_forces(F_net: np.ndarray, direction: np.ndarray) -> tuple[float, float]:
    """``(N_i, N_j)`` de un elemento de dos nodos a partir de sus fuerzas
    nodales netas ``F_net = F_int − f_eq`` (globales, layout ``[F_i, F_j]``)
    y de la dirección unitaria ``i → j``. Cara ``−x`` en ``i``:
    ``N_i = −F_i·d``; cara ``+x`` en ``j``: ``N_j = +F_j·d``. Compartido por
    trusses y cables (Reglas.md §5, tracción positiva)."""
    direction = np.asarray(direction, dtype=float)
    k = direction.size
    F_net = np.asarray(F_net, dtype=float)
    return -float(F_net[:k] @ direction), float(F_net[k:] @ direction)
