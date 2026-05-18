# fenix_fem/fenix/core/cohesive_material.py
"""Clase base abstracta para materiales cohesivos *traction-jump* (ADR 0010).

Familia paralela a ``Material``: opera sobre el salto ``[[u]]`` y devuelve
tracciones ``t`` sobre la superficie de discontinuidad ``Γ_d``, no sobre
deformaciones y esfuerzos del bulk.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar


class CohesiveMaterial(ABC):
    """Clase base de los materiales cohesivos de Fenix FEM (ADR 0010).

    Contrato de subclases
    ---------------------
    JUMP_DIM : ClassVar[int]
        Dimensión del vector de salto ``[[u]]`` que el material espera:
            2 → ``[[u]] = ([[u_n]], [[u_s]])``, frame local 2D
            3 → ``[[u]] = ([[u_n]], [[u_s1]], [[u_s2]])``, frame local 3D
        Los elementos con discontinuidad embebida validarán al construirse
        que ``cohesive_material.JUMP_DIM`` coincide con la dimensión que su
        formulación requiere.

    PRIMARY_STATE_VAR : ClassVar[str | None]
        Nombre de la clave principal del dict de estado para exportación
        al post-proceso (ej. ``'damage'``, ``'kappa'``). ``None`` si el
        material es elástico puro sin variables internas relevantes.

    IS_SYMMETRIC : ClassVar[bool]
        ``True`` si la contribución del material a la matriz tangente del
        elemento es simétrica. Lo es típicamente en daño isótropo Modo-I;
        modo mixto I-II y plasticidad cohesiva no asociada introducen
        asimetría. El despachador algebraico (ADR 0003) lo agrega aguas
        arriba para elegir Cholesky/LDLᵀ/LU.

    Notas sobre unidades
    --------------------
    Los cohesivos trabajan con parámetros físicos cuyas unidades difieren
    de los continuos: ``σ_t0`` (Pa, resistencia a tracción), ``G_F`` (N/m,
    energía de fractura), ``K_e`` (Pa/m, rigidez del salto / penalty). No
    tienen ``E``, ``ν`` ni densidad propia — la inercia es del bulk.
    """

    JUMP_DIM: ClassVar[int]
    PRIMARY_STATE_VAR: ClassVar[str | None] = None
    IS_SYMMETRIC: ClassVar[bool] = True

    @abstractmethod
    def compute_traction(self, jump, state_vars=None):
        """Calcula la tracción, la tangente algorítmica y el nuevo estado.

        Parameters
        ----------
        jump : np.ndarray
            Salto ``[[u]] ∈ ℝ^{JUMP_DIM}`` en el frame local ``(n, s)`` de
            ``Γ_d`` (en 3D ``(n, s1, s2)``). La transformación entre ejes
            globales y locales es responsabilidad del elemento consumidor.
        state_vars : dict | None
            Estado interno *committed* del paso anterior. ``None`` indica
            estado virgen (material intacto, sin historial).

        Returns
        -------
        (traction, tangent, new_state_vars)
            ``traction``: ``np.ndarray`` de tamaño ``JUMP_DIM``.
            ``tangent``: ``np.ndarray`` ``(JUMP_DIM, JUMP_DIM)`` — tangente
            algorítmica ``∂t/∂[[u]]`` en el frame local.
            ``new_state_vars``: dict con el estado *trial* tras este paso;
            el elemento decide cuándo *commit*.
        """
        pass
