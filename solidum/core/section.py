# solidum_fem/solidum/core/section.py
"""Contrato de las secciones de marco (ADR 0021).

Una **sección** es la ley constitutiva de un marco a nivel de sección
transversal: recibe las deformaciones generalizadas del eje en un punto de
Gauss y devuelve las fuerzas internas de la sección y su tangente. Es a un
marco lo que un material es a un sólido, pero con otro contrato (un material
recibe una deformación y devuelve un esfuerzo; una sección recibe
``e = [ε₀, κ]`` y devuelve ``N`` y ``M``), y por eso vive en su propia
familia (``SectionRegistry``) y no en la de materiales.

Magnitudes en los ejes locales del elemento, referidas al **eje del
elemento** (la línea de nodos), con la convención de Reglas.md §5:

- 2D: ``e = [ε₀, κ_z]``, ``s = [N, M_z]``.
- 3D: ``e = [ε₀, κ_y, κ_z]``, ``s = [N, M_y, M_z]``.

``N``, ``M_y`` y ``M_z`` son **fuerzas internas** (resultantes integradas del
esfuerzo ``σ_xx`` sobre la sección), no esfuerzos. ``M_z > 0`` tracciona las
fibras con ``y < 0``; ``M_y > 0``, las fibras con ``z > 0``.

**Sin estado.** Una sección se declara una vez y la comparten muchos
elementos. La historia (plasticidad, daño) de cada punto de la sección en
cada punto de Gauss vive en el elemento: el elemento guarda la lista que
devuelve :meth:`Section.initial_states`, se la pasa a ``response_*`` y guarda
la lista nueva como estado *trial*, con la misma semántica trial/commit que
el estado de un material en un sólido.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Tuple

import numpy as np


class Section(ABC):
    """Clase base abstracta de las secciones de marco (ADR 0021).

    Atributos
    ---------
    GJ : float | None
        Rigidez torsional elástica ``G·J`` de la sección. Opcional: ``None``
        si no se declara. La exige ``Frame3D`` (torsión elástica desacoplada,
        ``T = GJ·θₓ'``); los marcos 2D no la usan.

    Contrato
    --------
    response_2d(e, states) -> (s, k_s, states')
        Fuerzas internas ``s = [N, M_z]`` (2,), tangente ``k_s = ∂s/∂e``
        (2, 2) y estados nuevos (trial) para ``e = [ε₀, κ_z]``.
    response_3d(e, states) -> (s, k_s, states')
        Ídem en 3D con ``e = [ε₀, κ_y, κ_z]``, ``s = [N, M_y, M_z]``.
    strains_2d(e), strains_3d(e)
        Deformación axial de cada punto material de la sección para ``e``.
    initial_states()
        Lista de estados iniciales (uno por punto material de la sección)
        que el elemento guarda en cada punto de Gauss.
    materials()
        Materiales distintos que usa la sección, en orden de aparición.
    mass_properties()
        ``{"rhoA", "rhoIy", "rhoIz"}``: densidad lineal y momentos de
        inercia másicos por unidad de longitud respecto del eje del
        elemento. Lanza ``ValueError`` si falta la densidad de algún
        material.
    total_area()
        Área bruta ``∫ dA`` de la sección: la usa la carga de una fuerza de
        cuerpo uniforme ``b`` (fuerza por unidad de volumen), cuya carga por
        unidad de longitud es ``A·b``.
    """

    GJ: float | None = None

    @abstractmethod
    def response_2d(self, e, states) -> Tuple[np.ndarray, np.ndarray, list]:
        """``(s, k_s, states')`` en 2D; ``states`` no se modifica."""

    @abstractmethod
    def response_3d(self, e, states) -> Tuple[np.ndarray, np.ndarray, list]:
        """``(s, k_s, states')`` en 3D; ``states`` no se modifica."""

    @abstractmethod
    def strains_2d(self, e) -> np.ndarray:
        """Deformación axial de cada punto material para ``e = [ε₀, κ_z]``."""

    @abstractmethod
    def strains_3d(self, e) -> np.ndarray:
        """Deformación axial de cada punto material para ``e = [ε₀, κ_y, κ_z]``."""

    @abstractmethod
    def initial_states(self) -> list:
        """Estados iniciales, uno por punto material (lista nueva en cada llamada)."""

    @abstractmethod
    def materials(self) -> list:
        """Materiales distintos de la sección, en orden de aparición."""

    @abstractmethod
    def mass_properties(self) -> dict:
        """``{"rhoA", "rhoIy", "rhoIz"}`` por unidad de longitud."""

    @abstractmethod
    def total_area(self) -> float:
        """Área bruta de la sección."""
