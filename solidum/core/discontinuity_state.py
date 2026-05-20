# fenix_fem/solidum/core/discontinuity_state.py
"""Estado de una discontinuidad embebida (ADR 0010 §4).

Estructura paralela a :class:`ElementState` con semántica trial/commit, pero
específica de elementos con discontinuidad interior. Recoge la geometría de
``Γ_d`` (normal, tangente, centroide, nodo solitario, longitud efectiva),
los DOFs enriquecidos ``[[u]]`` y el estado interno del material cohesivo
asociado.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field

import numpy as np


@dataclass
class DiscontinuityState:
    """Estado de un elemento con discontinuidad activa.

    Attributes
    ----------
    normal : np.ndarray (2,)
        Normal unitaria a ``Γ_d``, fijada al activar (perpendicular a ``σ_I``
        en criterio Rankine). Apunta hacia el nodo solitario (``Ω⁺``).
    tangent : np.ndarray (2,)
        Tangente ``s = (−n_y, n_x)`` (regla de la mano derecha con ``+z``
        saliendo del plano, consistente con Reglas.md §5).
    centroid : np.ndarray (2,)
        Punto por el que pasa ``Γ_d`` (en CST, el centroide del triángulo).
    solitary_node : int
        Índice (0, 1 o 2 en CST) del nodo que queda en ``Ω⁺``. Los otros dos
        forman ``Ω⁻``. La función de soporte ``φ`` vale 1 en este nodo y 0
        en los otros.
    l_d : float
        Longitud efectiva de ``Γ_d`` dentro del elemento. Para CST KOS,
        ``l_d = (A_e/h)·cos(θ − α)`` por el Cap. 6 de Retama (2010).
    jump_committed : np.ndarray (2,)
        Salto ``[[u]] = ([[u_n]], [[u_s]])`` en frame local, al cierre del
        último paso convergido. Inicializado a cero al activar.
    jump_trial : np.ndarray (2,)
        Salto ``[[u]]`` dentro de la iteración Newton del paso actual.
    cohesive_state_committed : dict
        Estado interno del material cohesivo al cierre del último paso
        convergido (típicamente ``{'kappa': ..., 'damage': ...}``).
    cohesive_state_trial : dict
        Estado interno del cohesivo en la iteración Newton actual.

    Notes
    -----
    La activación es **irreversible**: una vez creada una ``DiscontinuityState``
    no se elimina aunque el siguiente paso descargue el elemento. Se mantiene
    para que el material cohesivo conserve su historial ``κ`` (damage no
    decrece, Kuhn-Tucker).
    """

    normal: np.ndarray
    tangent: np.ndarray
    centroid: np.ndarray
    solitary_node: int
    l_d: float
    jump_committed: np.ndarray = field(default_factory=lambda: np.zeros(2))
    jump_trial: np.ndarray = field(default_factory=lambda: np.zeros(2))
    cohesive_state_committed: dict = field(default_factory=dict)
    cohesive_state_trial: dict = field(default_factory=dict)

    def commit(self) -> None:
        """Fija el estado trial (jump + cohesivo) como nuevo committed.

        Se invoca tras la convergencia del Newton global, en
        :meth:`Element.commit_state`. Usa copia profunda para que committed
        y trial queden completamente independientes.
        """
        self.jump_committed = np.copy(self.jump_trial)
        self.cohesive_state_committed = copy.deepcopy(self.cohesive_state_trial)
