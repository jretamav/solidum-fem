"""Cargas distribuidas consistentes: fuerza de cuerpo y carga de borde (ADR 0022).

Una carga distribuida es un **valor constante** o un **campo**, una función de
la posición. Los elementos la integran en sus puntos de Gauss sobre la
geometría real (isoparamétrica):

- fuerza de cuerpo ``f_i = ∫ N_i b(x) dV``, con ``b(x)`` evaluado en cada
  punto: :func:`body_load`;
- carga de borde 2D ``f_k = ∫_Γe N_k t(x, n) ds``, con ``n`` la normal
  exterior unitaria del borde **discreto** en el punto: :func:`edge_load`;
- carga de cara 3D ``f_k = ∫_Γe N_k t(x, n) dA``, igual sobre la cara:
  :func:`face_load`.

Que el campo de borde reciba la normal tiene dos usos. Una presión es
``t(x, n) = −p·n``, sin un método aparte. Y una solución manufacturada
``u*`` se impone con ``t = σ(u*)·n`` sobre la geometría de la malla, de modo
que ``u*`` es la solución exacta del problema discreto en su propio dominio
(``solidum.verification``).

El borde se integra con su interpolación de Lagrange (2 nodos: lineal; 3
nodos en el orden ``(inicio, medio, fin)``: cuadrática), Gauss de 3 puntos y
el jacobiano ``|dx/ds|``. Es exacto para tracción constante en un borde recto,
con el nodo medio centrado o no, y la resultante de una presión uniforme lo es
sobre cualquier borde cuadrático (``∫ n ds`` es la cuerda girada). Hasta el
2026-10-27 los bordes cuadráticos repartían la cuerda con pesos fijos 1-4-1,
exactos sólo con el nodo medio centrado sobre un borde recto (deuda #13).
"""
from __future__ import annotations

import numpy as np

# Gauss-Legendre de 3 puntos en [-1, 1]: exacta hasta grado 5, de sobra para
# N_k·|dx/ds| con N_k cuadrática y borde recto (grado 3).
_EDGE_POINTS, _EDGE_WEIGHTS = np.polynomial.legendre.leggauss(3)


def _lagrange_1d(n_nodes: int, s: np.ndarray):
    """Funciones de forma de Lagrange en ``s ∈ [-1, 1]`` y sus derivadas,
    ``(n_g, n_nodes)``, con los nodos en el orden de ``EDGE_NODES``:
    ``(-1, 1)`` o ``(-1, 0, 1)``."""
    if n_nodes == 2:
        N = np.column_stack([(1.0 - s) / 2.0, (1.0 + s) / 2.0])
        dN = np.column_stack([np.full_like(s, -0.5), np.full_like(s, 0.5)])
    elif n_nodes == 3:
        N = np.column_stack([s * (s - 1.0) / 2.0, 1.0 - s * s, s * (s + 1.0) / 2.0])
        dN = np.column_stack([s - 0.5, -2.0 * s, s + 0.5])
    else:
        raise ValueError(f"Borde de {n_nodes} nodos no admitido (2 ó 3).")
    return N, dN


def field_values(value, points, n_components: int, name: str, *extra) -> np.ndarray:
    """Valores ``(n_g, n_components)`` de una carga en los puntos ``points``.

    Una constante se repite; un campo se evalúa como ``value(x, *extra)`` en
    cada punto, con ``extra`` los argumentos adicionales por punto (la normal
    de un borde)."""
    n_g = len(points)
    if not callable(value):
        arr = np.asarray(value, dtype=np.float64).reshape(-1)
        if arr.size != n_components:
            raise ValueError(f"La carga {name} tiene {arr.size} componentes; "
                             f"se esperaban {n_components}.")
        return np.broadcast_to(arr, (n_g, n_components))
    out = np.empty((n_g, n_components))
    for g in range(n_g):
        v = np.asarray(value(points[g], *(e[g] for e in extra)), dtype=np.float64).reshape(-1)
        if v.size != n_components:
            raise ValueError(f"El campo {name} devolvió {v.size} componentes en "
                             f"x = {points[g]}; se esperaban {n_components}.")
        out[g] = v
    return out


def body_load(N: np.ndarray, dV: np.ndarray, x: np.ndarray, b,
              n_components: int | None = None) -> np.ndarray:
    """``f_i = Σ_g N_i(ξ_g)·b(x_g)·dV_g``.

    Parameters
    ----------
    N : (n_g, n) funciones de forma en los puntos de Gauss.
    dV : (n_g,) diferencial ``detJ·w`` (por el espesor en los sólidos planos).
    x : (n_g, d) posición global de los puntos.
    b : constante o campo ``b(x)``.
    n_components : componentes de ``b``; por omisión ``d`` (carga vectorial).

    Returns
    -------
    (n·n_components,) en el orden ``[f_1x, f_1y, ..., f_nx, f_ny]``.
    """
    m = x.shape[1] if n_components is None else n_components
    values = field_values(b, x, m, "de cuerpo")
    return ((N * dV[:, None]).T @ values).reshape(-1)


def edge_load(X: np.ndarray, t, n_components: int = 2) -> np.ndarray:
    """``f_k = ∫ N_k·t(x, n) ds`` sobre un borde 2D de 2 ó 3 nodos.

    ``X`` ``(k, 2)`` son las coordenadas de los nodos del borde en el orden de
    ``EDGE_NODES``, que recorre el contorno del elemento en sentido
    antihorario; la normal exterior es la tangente girada −90°. ``t`` es una
    constante o un campo ``t(x, n)``. Devuelve ``(k·n_components,)``.
    """
    X = np.asarray(X, dtype=np.float64)[:, :2]
    N, dN = _lagrange_1d(len(X), _EDGE_POINTS)
    g = dN @ X
    ds = np.hypot(g[:, 0], g[:, 1])
    if np.any(ds <= 0.0):
        raise ValueError("Borde degenerado: nodos coincidentes.")
    normal = np.column_stack([g[:, 1], -g[:, 0]]) / ds[:, None]
    values = field_values(t, N @ X, n_components, "de borde", normal)
    return ((N * (ds * _EDGE_WEIGHTS)[:, None]).T @ values).reshape(-1)


def face_load(N: np.ndarray, dN: np.ndarray, w: np.ndarray, X: np.ndarray, t,
              n_components: int = 3) -> np.ndarray:
    """``f_k = ∫ N_k·t(x, n) dA`` sobre una cara 3D.

    ``N`` ``(n_g, k)`` y ``dN`` ``(n_g, 2, k)`` son las funciones de forma de
    la cara y sus derivadas respecto de sus parámetros ``(s, r)`` en los
    puntos de la regla de pesos ``w``; ``X`` ``(k, 3)``, las coordenadas de
    los nodos de la cara en el orden de ``FACE_NODES``, antihorario visto
    desde fuera, así que ``∂x/∂s × ∂x/∂r`` apunta hacia fuera. ``t`` es una
    constante o un campo ``t(x, n)``. Devuelve ``(k·n_components,)``.
    """
    X = np.asarray(X, dtype=np.float64)[:, :3]
    c = np.cross(dN[:, 0, :] @ X, dN[:, 1, :] @ X)
    dA = np.linalg.norm(c, axis=1)
    if np.any(dA <= 0.0):
        raise ValueError("Cara degenerada: nodos coincidentes o alineados.")
    values = field_values(t, N @ X, n_components, "de cara", c / dA[:, None])
    return ((N * (dA * np.asarray(w))[:, None]).T @ values).reshape(-1)
