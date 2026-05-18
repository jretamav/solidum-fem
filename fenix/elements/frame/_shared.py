"""Utilidades compartidas por los elementos del paquete ``fenix.elements.frame``.

- :func:`build_geometry_2d` — longitud, cosenos directores y matriz de
  transformación ``T`` 6×6 para un elemento 2D de dos nodos. Antes
  estaba duplicado bit-a-bit en ``Frame2DEuler._build_geometry`` y
  ``Frame2DTimoshenko._build_geometry``.
- :func:`_frame2d_consistent_body_load` — vector nodal consistente con
  carga distribuida uniforme; idéntico para Euler, Timoshenko (Hermite
  blended bajo carga uniforme) y EulerCorot (config de referencia).
- :func:`_frame2d_consistent_mass_local` — matriz de masa 6×6 en ejes
  locales: axial lineal + Hermitiana cúbica + inercia rotacional de
  sección (ADR 0009).
- :func:`_frame2d_forces_from_local` — traduce el vector local
  ``[Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]`` a ``ElementForces`` en la
  convención de viga estructural de Reglas.md §5.
"""
import math
from typing import List

import numpy as np

from fenix.core.node import Node
from fenix.logging import get_logger
from fenix.results import ElementForces

# Logger compartido por todos los elementos del paquete `frame`. Los
# módulos consumidores hacen `from fenix.elements.frame._shared import _log`
# en vez de declarar uno propio.
_log = get_logger("elements.frame")


def build_geometry_2d(nodes: List[Node]):
    """Calcula longitud, cosenos directores y matriz ``T`` 6×6.

    Para vigas 2D de dos nodos con DOFs ``[ux, uy, rz]`` por nodo. La
    matriz ``T`` mapea ejes globales a locales: filas 0-2 corresponden al
    nodo i, filas 3-5 al nodo j; las rotaciones ``rz`` no se rotan.

    Parameters
    ----------
    nodes
        Lista de exactamente 2 :class:`Node` extremos del elemento.

    Returns
    -------
    L0, c, s, T
        ``L0`` longitud (en configuración de referencia para corotacionales),
        ``c`` y ``s`` cosenos directores del eje local x, y ``T`` la matriz
        6×6 global → local.

    Raises
    ------
    ValueError
        Si la longitud del elemento es cero.
    """
    coords1 = nodes[0].coordinates
    coords2 = nodes[1].coordinates
    dx, dy = coords2[0] - coords1[0], coords2[1] - coords1[1]
    L0 = math.sqrt(dx**2 + dy**2)
    if L0 == 0.0:
        raise ValueError("La longitud del elemento no puede ser cero.")
    c = dx / L0
    s = dy / L0
    T = np.array([
        [ c,  s, 0, 0, 0, 0],
        [-s,  c, 0, 0, 0, 0],
        [ 0,  0, 1, 0, 0, 0],
        [ 0,  0, 0, c, s, 0],
        [ 0,  0, 0,-s, c, 0],
        [ 0,  0, 0, 0, 0, 1],
    ])
    return L0, c, s, T


def _frame2d_consistent_body_load(b: np.ndarray, A: float, L: float,
                                   T: np.ndarray) -> np.ndarray:
    """Vector nodal consistente con ∫NᵀbA dx para una viga 2D de 2 nodos.

    Compartido por ``Frame2DEuler``, ``Frame2DTimoshenko`` y
    ``Frame2DEulerCorot`` porque la fórmula es idéntica: shape lineal en
    axial (mitad-mitad) y shape Hermite cúbico en transversal (qL/2 en
    fuerza, ±qL²/12 en momento). Para Timoshenko con interpolación
    Hermite-blended el resultado coincide exactamente con Euler-Bernoulli
    en el caso de carga uniforme, que es lo que cubre este método.

    Parameters
    ----------
    b : np.ndarray, shape (2,)
        Fuerza de cuerpo por unidad de volumen en ejes globales.
    A : float
        Área de la sección transversal.
    L : float
        Longitud del elemento (en configuración de referencia para
        formulaciones corotacionales).
    T : np.ndarray, shape (6, 6)
        Matriz de transformación local↔global del elemento.

    Returns
    -------
    np.ndarray, shape (6,)
        ``[Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]`` en ejes globales.
    """
    b = np.asarray(b, dtype=float).reshape(2)
    # Carga distribuida por unidad de longitud (en globales).
    q_global = A * b
    # Transformar a ejes locales: T[0:2, 0:2] es la sub-rotación 2D global→local.
    R = T[0:2, 0:2]
    q_local = R @ q_global  # (q_axial, q_transversal)
    qa, qt = q_local[0], q_local[1]
    f_local = np.array([
        0.5 * qa * L,
        0.5 * qt * L,
        qt * L * L / 12.0,
        0.5 * qa * L,
        0.5 * qt * L,
        -qt * L * L / 12.0,
    ])
    return T.T @ f_local


def _frame2d_consistent_mass_local(rho: float, A: float, I: float,
                                     L: float) -> np.ndarray:
    """Matriz de masa consistente 6×6 en ejes locales del elemento 2D (ADR 0009).

    Tres contribuciones:

    - **Axial** (ux_i, ux_j): masa traslacional lineal ``(ρAL/6)·[[2,1],[1,2]]``.
    - **Flexión transversal** (uy_i, rz_i, uy_j, rz_j): masa Hermitiana cúbica
      ``(ρAL/420)·H``, idéntica para Euler-Bernoulli y Timoshenko (la
      densidad lineal ``ρA`` no depende ni de ``I`` ni de ``As/Φ``).
    - **Inercia rotacional propia de sección** (rz_i, rz_j): masa rotacional
      ``(ρIL/6)·[[2,1],[1,2]]`` que recoge la energía cinética asociada a la
      rotación de la sección alrededor del nodo. Estándar en vigas peraltadas
      (Timoshenko) y mejora pequeña pero correcta en Bernoulli. Sin esta
      contribución, M sería rotacionalmente "singular" en ausencia de inercia
      a través de la deflexión Hermitiana — patrón consistente con la
      rigidez ``K_local`` que sí tiene términos rotacionales puros.

    Layout local: ``[ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]``.

    Parameters
    ----------
    rho : float
        Densidad del material.
    A : float
        Área de la sección.
    I : float
        Momento de inercia de la sección respecto al eje perpendicular al
        plano de flexión (``Iz`` en convención 3D). Mismo parámetro usado
        en la rigidez.
    L : float
        Longitud del elemento (en configuración de referencia para
        formulaciones corotacionales — Lagrangeano total, ADR 0009).
    """
    c_ax = rho * A * L / 6.0
    c_tr = rho * A * L / 420.0
    c_rot = rho * I * L / 6.0
    L2 = L * L
    M = np.zeros((6, 6))
    # Axial (ux_i, ux_j): masa de barra lineal 1D.
    M[0, 0] = 2.0 * c_ax;  M[0, 3] = 1.0 * c_ax
    M[3, 0] = 1.0 * c_ax;  M[3, 3] = 2.0 * c_ax
    # Flexional (uy_i, rz_i, uy_j, rz_j): masa Hermitiana cúbica.
    H = c_tr * np.array([
        [156.0,   22.0 * L,   54.0,  -13.0 * L],
        [22.0 * L,  4.0 * L2,  13.0 * L, -3.0 * L2],
        [54.0,    13.0 * L,  156.0,  -22.0 * L],
        [-13.0 * L, -3.0 * L2, -22.0 * L,  4.0 * L2],
    ])
    idx = [1, 2, 4, 5]
    for i, gi in enumerate(idx):
        for j, gj in enumerate(idx):
            M[gi, gj] = H[i, j]
    # Inercia rotacional propia de sección (rz_i, rz_j): ρI·L masa rotacional.
    M[2, 2] += 2.0 * c_rot;  M[2, 5] += 1.0 * c_rot
    M[5, 2] += 1.0 * c_rot;  M[5, 5] += 2.0 * c_rot
    return M


def _frame2d_lumped_mass_local(rho: float, A: float, I: float,
                                 L: float) -> np.ndarray:
    """Matriz de masa lumped 6×6 en ejes locales del frame 2D (ADR 0009 fase 2).

    Esquema **nodal directo**: la masa se distribuye en partes iguales entre
    los dos nodos del elemento, con la misma magnitud en cada dirección
    traslacional y una inercia rotacional proporcional al momento de
    inercia geométrico de la sección:

    .. math::

        m_t = \\rho A L / 2 \\quad \\text{(por DOF traslacional)}

        m_r = \\rho I L / 2 \\quad \\text{(por DOF rotacional)}

    Layout local: ``[ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]``.

    Propiedades:

    - Preserva la masa total ``ρAL`` exactamente (suma sobre nodos en
      cada dirección).
    - Preserva la inercia rotacional propia total ``ρIL`` (suma sobre
      nodos).
    - **Diagonal en global** tras la rotación ``T^T·M·T``: el bloque
      traslacional 2×2 por nodo es ``m_t·I_2``, invariante a rotación 2D;
      el bloque rotacional 1×1 es escalar. Esta diagonalidad en global es
      la condición que vuelve trivial la inversión ``M⁻¹`` requerida por
      diferencias centradas (ADR 0009 fase 5).

    Frente al lumping HRZ canónico (Hinton-Rock-Zienkiewicz 1976) aplicado
    a la diagonal consistente: el HRZ daría masas distintas en los DOFs
    ``ux_local`` (lineal) y ``uy_local`` (Hermitiana) — diagonal en local
    pero **no** en global tras rotar. El esquema nodal directo sacrifica
    la "fidelidad a la diagonal consistente" a cambio de la diagonalidad
    global, que es la propiedad operacional relevante para integradores
    explícitos. Es el esquema lumped estándar en frames (SAP2000, OpenSees,
    Cook-Malkus-Plesha §11.4).

    Parameters
    ----------
    rho, A, I, L
        Mismos parámetros físicos que :func:`_frame2d_consistent_mass_local`.

    Returns
    -------
    np.ndarray, shape (6, 6)
        Diagonal, positiva definida.
    """
    m_t = rho * A * L / 2.0
    m_r = rho * I * L / 2.0
    return np.diag([m_t, m_t, m_r, m_t, m_t, m_r])


def _frame2d_forces_from_local(F_local: np.ndarray) -> ElementForces:
    """Traduce un vector F_local (6,) en convención stress-resultant interna a
    ``ElementForces`` en convención de viga estructural (Reglas.md §5).

    F_local tiene semántica ``K·u = F_ext``: fuerzas/momentos externos sobre los
    nodos del elemento en ejes locales, layout ``[Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]``.

    Mapeo §5:
      - Nodo i (cara −x): N = −F_local[0], V = +F_local[1], M = −F_local[2]
      - Nodo j (cara +x): N = +F_local[3], V = −F_local[4], M = +F_local[5]

    Verificado con cantilever (carga en punta → M hogging) y cantilever con
    momento en punta (M constante sagging).
    """
    return ElementForces(
        kind="frame2d",
        components={
            "N": np.array([-F_local[0],  F_local[3]]),
            "V": np.array([ F_local[1], -F_local[4]]),
            "M": np.array([-F_local[2],  F_local[5]]),
        },
    )
