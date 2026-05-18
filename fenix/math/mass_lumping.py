"""Algoritmos de mass lumping (ADR 0009 fase 2).

Implementa el esquema HRZ (Hinton-Rock-Zienkiewicz 1976) para sólidos y
elementos isoparamétricos. Para vigas/marcos cada elemento usa **lumping
nodal directo** con fórmula cerrada (masa traslacional uniforme por nodo +
inercia rotacional proporcional a la inercia geométrica de sección): es el
único esquema que sobrevive a la rotación local→global manteniendo
``M_lumped`` diagonal — propiedad indispensable para diferencias centradas
(ADR 0009 fase 5). Cuando el bloque rotacional no es escalar (Frame3D con
``Iy ≠ Iz``), la rotación deja un bloque 3×3 lleno en los DOFs rotacionales;
es una limitación documentada del lumped en frames 3D y es estándar en la
literatura (Cook-Malkus-Plesha §11.4).

El esquema HRZ canónico:

.. math::

    \\alpha = \\frac{D \\cdot m_\\text{total}}{\\sum_{i \\in \\text{trasl}} M_{ii}^\\text{cons}}

    M^\\text{lumped}_{ii} = \\alpha \\cdot M^\\text{cons}_{ii}
    \\quad \\forall i \\text{ (incluyendo rotacionales)}

donde ``D`` es el número de direcciones traslacionales (2 en 2D, 3 en 3D) y
``m_total = ρ·V_e`` es la masa total del elemento. Preserva la masa total
sumada sobre las ``D`` direcciones traslacionales, mantiene la
proporcionalidad de la diagonal consistente (evita masas nulas o negativas
en elementos de orden superior — Tri6, Quad8, Quad9), y deja los DOFs
rotacionales con inercia positiva (``M`` no singular para modal/dinámico).

Referencias:

- Hinton, E., Rock, T., Zienkiewicz, O. C. (1976). *A note on mass lumping
  and related processes in the finite element method*. Earthquake Eng.
  Struct. Dyn., 4(3), 245-249.
- Bathe, K.-J. (2014). *Finite Element Procedures*, §9.2.4.
- Cook, R. D., Malkus, D. S., Plesha, M. E. (2002). *Concepts and
  Applications of Finite Element Analysis*, §11.4.
"""
from __future__ import annotations

import numpy as np


def lump_hrz(
    M_consistent: np.ndarray,
    *,
    total_mass: float,
    n_translational_dirs: int,
    translational_dofs: np.ndarray | None = None,
) -> np.ndarray:
    """Aplica el esquema HRZ a una matriz consistente y devuelve la lumped.

    Parameters
    ----------
    M_consistent : (n, n) ndarray
        Matriz de masa consistente del elemento, ya integrada y simétrica.
    total_mass : float
        Masa física total ``ρ·V_e`` del elemento (un escalar). Se usa para
        normalizar el factor de escala α de modo que la suma de los DOFs
        traslacionales lumped iguale ``n_translational_dirs · total_mass``.
    n_translational_dirs : int
        Número de direcciones traslacionales del elemento (2 en 2D, 3 en 3D).
    translational_dofs : ndarray of int or None, optional
        Índices de los DOFs traslacionales. Si ``None``, asume que **todos**
        los DOFs son traslacionales (caso usual en sólidos isoparamétricos).

    Returns
    -------
    (n, n) ndarray
        Matriz diagonal. Devuelta como ``np.diag(...)`` ndarray densa por
        consistencia con el contrato de ``Element.compute_mass_matrix``.

    Raises
    ------
    ValueError
        Si la suma diagonal de los DOFs traslacionales no es positiva
        (matriz consistente mal construida o cuadratura insuficiente que
        deja modos nulos espurios).
    """
    diag_M = np.diag(M_consistent).copy()
    n = len(diag_M)
    if translational_dofs is None:
        translational_idx = np.arange(n)
    else:
        translational_idx = np.asarray(translational_dofs, dtype=int)

    sum_translational = float(np.sum(diag_M[translational_idx]))
    if sum_translational <= 0.0:
        raise ValueError(
            "lump_hrz: la suma diagonal de los DOFs traslacionales no es "
            f"positiva ({sum_translational:.3e}). La matriz consistente está "
            "mal construida o la cuadratura no integra exactamente el producto "
            "Nᵀ·N (modos nulos espurios en la diagonal)."
        )

    alpha = (n_translational_dirs * float(total_mass)) / sum_translational
    return np.diag(alpha * diag_M)
