"""Utilidades compartidas por los elementos del paquete ``solidum.elements.frame``.

- :func:`build_geometry_2d` — longitud, cosenos directores y matriz de
  transformación ``T`` 6×6 para un elemento 2D de dos nodos. Antes
  estaba duplicado bit-a-bit en ``Frame2DEuler._build_geometry`` y
  ``Frame2DTimoshenko._build_geometry``.
- :func:`_frame2d_consistent_body_load` — vector nodal consistente con
  carga distribuida uniforme; idéntico para Euler, Timoshenko (Hermite
  blended bajo carga uniforme) y EulerCorot (config de referencia).
- :func:`hermite_consistent_mass` y :func:`hrz_rotational_mass` — bloques
  de masa de la flexión con los cúbicos de Hermite (consistente 4×4 y la
  masa rotacional concentrada HRZ), compartidos con ``Frame3D``.
- :func:`_frame2d_consistent_mass_local` y :func:`_frame2d_lumped_mass_local`
  — matrices de masa 6×6 en ejes locales (ADR 0009): Euler-Bernoulli sin
  inercia rotatoria por defecto; Timoshenko pasa su inercia rotatoria ``ρI``.
- :func:`_frame2d_forces_from_local` — traduce el vector local
  ``[Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]`` a ``ElementForces`` en la
  convención de viga estructural de Reglas.md §5.
- :class:`Frame2DLine` — campo interior y fuerza de cuerpo como campo de los
  marcos 2D (ADR 0022, fase 4).
"""
import math
from typing import List

import numpy as np

from solidum.core.node import Node
from solidum.elements._line import LineElementMixin, hermite_functions, linear_tangent
from solidum.logging import get_logger
from solidum.results import ElementForces

# Logger compartido por todos los elementos del paquete `frame`. Los
# módulos consumidores hacen `from solidum.elements.frame._shared import _log`
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


class Frame2DLine(LineElementMixin):
    """Campo interior de los marcos 2D (ADR 0022, fase 4) en ejes locales:
    axial lineal y flexión con los cúbicos de Hermite (Euler-Bernoulli,
    ``θ = v'`` y ``κ = v''``); el ``Frame2DTimoshenko`` sustituye la
    flexión por su propia interpolación. Con él, la fuerza de cuerpo como
    campo se reparte con las mismas funciones que dan la rigidez."""
    LINE_FIELDS = ("u", "v", "theta_z")
    LINE_STRAINS = ("eps", "kappa_z")

    def _line_transform(self) -> np.ndarray:
        """``T`` 6×6 global → local de la geometría de referencia."""
        return self.T

    def _to_local(self, u_e):
        return self._line_transform() @ u_e

    def _to_global(self, f_local):
        return self._line_transform().T @ f_local

    def _load_axes(self):
        return self._line_transform()[0:2, 0:2]

    def _line_interpolation(self, xi):
        L = self.L0
        n = len(xi)
        N = np.zeros((n, 3, 6))
        B = np.zeros((n, 2, 6))
        H, dH, d2H = hermite_functions(xi, L)
        N[:, 0, 0], N[:, 0, 3] = 1.0 - xi, xi
        N[:, 1, [1, 2, 4, 5]] = H
        N[:, 2, [1, 2, 4, 5]] = dH
        B[:, 0, 0], B[:, 0, 3] = -1.0 / L, 1.0 / L
        B[:, 1, [1, 2, 4, 5]] = d2H
        return N, B

    def section_rigidities(self) -> np.ndarray:
        """``[EA, EI]``; sólo con material elástico lineal (con sección por
        fibras la rigidez no es una constante)."""
        if getattr(self, "section", None) is not None:
            raise ValueError(f"{type(self).__name__}(id={self.id}): con sección por fibras "
                             "las rigideces de la sección no son constantes.")
        E = linear_tangent(self.material)
        return np.array([E * self.A, E * self.I])


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


def hermite_consistent_mass(rhoA: float, L: float) -> np.ndarray:
    """Masa consistente 4×4 de la flexión, ``∫₀ᴸ ρA NᵀN dx`` con los cúbicos
    de Hermite, en ``(v₁, θ₁, v₂, θ₂)`` con ``θ = dv/dx`` (plano ``xy``
    local; Cook-Malkus-Plesha, masa consistente de viga)::

        ρAL/420 · [[ 156,   22L,   54,  −13L],
                   [ 22L,   4L²,  13L,  −3L²],
                   [  54,   13L,  156,  −22L],
                   [−13L,  −3L², −22L,   4L²]]

    Es la energía cinética de la traslación transversal ``v(x)`` y nada más:
    no contiene inercia rotatoria de la sección. Reproduce exactamente la
    energía cinética de los movimientos rígidos (traslación y giro), que la
    interpolación cúbica contiene.

    Parameters
    ----------
    rhoA : float
        Densidad lineal ``ρA`` (masa por unidad de longitud).
    L : float
        Longitud del elemento.
    """
    L2 = L * L
    return (rhoA * L / 420.0) * np.array([
        [156.0,      22.0 * L,   54.0,      -13.0 * L],
        [22.0 * L,    4.0 * L2,  13.0 * L,   -3.0 * L2],
        [54.0,       13.0 * L,  156.0,      -22.0 * L],
        [-13.0 * L,  -3.0 * L2, -22.0 * L,    4.0 * L2],
    ])


def hrz_rotational_mass(rhoA: float, L: float) -> float:
    """Masa rotacional concentrada por nodo de la flexión, ``m_θ = ρAL³/78``.

    Lumping HRZ (Hinton, Rock y Zienkiewicz 1976; Cook-Malkus-Plesha):
    cada término de la diagonal consistente de :func:`hermite_consistent_mass`
    se escala por el factor que hace que los términos traslacionales sumen
    la masa ``ρAL``. La diagonal es ``ρAL/420·(156, 4L², 156, 4L²)`` y el
    factor ``α = ρAL/(312·ρAL/420) = 420/312``, de modo que ``m_v = ρAL/2``
    y ``m_θ = 4L²·ρAL/420·420/312 = ρAL³/78``.
    """
    return rhoA * L**3 / 78.0


def _frame2d_consistent_mass_local(rhoA: float, L: float,
                                   rotary_inertia: float | None = None) -> np.ndarray:
    """Matriz de masa consistente 6×6 en ejes locales del marco 2D (ADR 0009).

    - **Axial** ``(ux_i, ux_j)``: masa de barra lineal ``(ρAL/6)·[[2,1],[1,2]]``.
    - **Flexión** ``(uy_i, rz_i, uy_j, rz_j)``: masa hermítica
      :func:`hermite_consistent_mass`.
    - **Inercia rotatoria de la sección** ``(rz_i, rz_j)``, sólo si se pasa
      ``rotary_inertia = ρI``: ``(ρIL/6)·[[2,1],[1,2]]``, la energía cinética
      del giro de la sección interpolado linealmente. Es cinemática de
      Timoshenko (giro independiente de la pendiente) y la usa sólo
      ``Frame2DTimoshenko``; Euler-Bernoulli la desprecia por hipótesis y
      no pasa el argumento.

    Layout local: ``[ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]``.

    Parameters
    ----------
    rhoA : float
        Densidad lineal ``ρA`` (``ρ·A`` con material; ``Σρ_i A_i`` con
        sección).
    L : float
        Longitud del elemento (en la configuración de referencia en las
        formulaciones corrotacionales, ADR 0009).
    rotary_inertia : float, optional
        ``ρI`` por unidad de longitud de la inercia rotatoria de la sección;
        ``None`` (defecto) la desprecia (Euler-Bernoulli).
    """
    c_ax = rhoA * L / 6.0
    M = np.zeros((6, 6))
    # Axial (ux_i, ux_j): masa de barra lineal 1D.
    M[0, 0] = 2.0 * c_ax;  M[0, 3] = 1.0 * c_ax
    M[3, 0] = 1.0 * c_ax;  M[3, 3] = 2.0 * c_ax
    # Flexión (uy_i, rz_i, uy_j, rz_j): masa hermítica.
    idx = [1, 2, 4, 5]
    M[np.ix_(idx, idx)] = hermite_consistent_mass(rhoA, L)
    if rotary_inertia is not None:
        # Inercia rotatoria de la sección (Timoshenko): ρI·L/6·[[2,1],[1,2]].
        c_rot = rotary_inertia * L / 6.0
        M[2, 2] += 2.0 * c_rot;  M[2, 5] += 1.0 * c_rot
        M[5, 2] += 1.0 * c_rot;  M[5, 5] += 2.0 * c_rot
    return M


def _frame2d_lumped_mass_local(rhoA: float, L: float,
                               rotary_inertia: float | None = None) -> np.ndarray:
    """Matriz de masa concentrada 6×6 en ejes locales del marco 2D (ADR 0009).

    Masa traslacional ``m_t = ρAL/2`` por nodo en ``ux`` y en ``uy`` (el HRZ
    de la axial lineal y de la flexión hermítica da el mismo valor). Masa
    rotacional por nodo:

    - Euler-Bernoulli (``rotary_inertia=None``, defecto): la HRZ de la
      flexión, ``m_θ = ρAL³/78`` (:func:`hrz_rotational_mass`). Sin inercia
      rotatoria de la sección, coherente con la masa consistente.
    - Timoshenko (``rotary_inertia = ρI``): ``m_r = ρIL/2``, la inercia
      rotatoria de la sección repartida entre los nodos.

    Layout local: ``[ux_i, uy_i, rz_i, ux_j, uy_j, rz_j]``.

    Propiedades:

    - Conserva la masa total ``ρAL`` en cada dirección traslacional.
    - **Diagonal también en globales** tras ``Tᵀ·M·T``: el bloque
      traslacional por nodo es ``m_t·I₂``, invariante ante la rotación 2D,
      y el rotacional es escalar. Es lo que vuelve trivial ``M⁻¹`` en
      diferencias centradas (ADR 0009 fase 5).

    Parameters
    ----------
    rhoA, L, rotary_inertia
        Como en :func:`_frame2d_consistent_mass_local`.

    Returns
    -------
    np.ndarray, shape (6, 6)
        Diagonal, definida positiva.
    """
    m_t = rhoA * L / 2.0
    if rotary_inertia is None:
        m_r = hrz_rotational_mass(rhoA, L)
    else:
        m_r = rotary_inertia * L / 2.0
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
