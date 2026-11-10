"""Materiales de interfaz: leyes tracción–desplazamiento relativo (ADR 0033).

Familia paralela a ``Material``: no relaciona el esfuerzo con la deformación
en Voigt, sino la **tracción** ``t`` sobre una interfaz con el
**desplazamiento relativo** ``Δu`` entre sus dos lados, ambos en el sistema
local de la interfaz. La usan los elementos que tienen una interfaz: una
superficie de discontinuidad (grieta, junta) o la línea de una barra embebida
en un sólido (adherencia-deslizamiento).

Orden de las componentes (normales primero, tangenciales después):

- ``GEOMETRY = "surface"`` (superficie: grieta, junta; codimensión 1): una
  normal ``n`` y ``d − 1`` tangenciales: ``(n, s)`` en 2D, ``(n, s1, s2)``
  en 3D.
- ``GEOMETRY = "line"`` (línea de una barra dentro de un sólido): ``d − 1``
  normales y una tangencial a lo largo de la barra: ``(n, s)`` en 2D,
  ``(n1, n2, s)`` en 3D.

En 2D las dos geometrías tienen el mismo orden. El paso entre ejes globales y
locales y el ancho sobre el que actúa cada componente (en una barra, el
perímetro para la tangencial) son del elemento, no de la ley.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar

GEOMETRIES = ("surface", "line")


class InterfaceMaterial(ABC):
    """Clase base de los materiales de interfaz (ADR 0033).

    Contrato de las subclases
    -------------------------
    GEOMETRY : ClassVar[str]
        ``"surface"`` o ``"line"``: fija qué componentes del desplazamiento
        relativo son normales y cuáles tangenciales (ver el módulo).
    JUMP_DIM : ClassVar[int | None]
        Número de componentes que admite la ley; ``None`` si admite las de
        2D y las de 3D. El elemento consumidor lo valida al construirse.
    PRIMARY_STATE_VAR : ClassVar[str | None]
        Variable interna principal para el post-proceso (``'damage'``,
        ``'slip_max'``…); ``None`` sin variables internas.
    IS_SYMMETRIC : ClassVar[bool]
        ``True`` si la tangente ``∂t/∂Δu`` es simétrica (ADR 0003).

    Unidades: la tracción es fuerza por unidad de área (Pa) y las rigideces de
    la interfaz, fuerza por unidad de área y de longitud (Pa/m). ``density``
    vale 0.0: la inercia es de los elementos que la rodean.
    """

    GEOMETRY: ClassVar[str]
    JUMP_DIM: ClassVar[int | None] = None
    # Sin masa por diseño (ADR 0008: density = 0.0 declarada): la inercia es
    # de los elementos que rodean la interfaz.
    density: ClassVar[float] = 0.0
    PRIMARY_STATE_VAR: ClassVar[str | None] = None
    IS_SYMMETRIC: ClassVar[bool] = True

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        geometry = cls.__dict__.get("GEOMETRY")
        if geometry is not None and geometry not in GEOMETRIES:
            raise TypeError(
                f"{cls.__name__}.GEOMETRY = {geometry!r}: debe ser uno de {GEOMETRIES}.")

    def accepts(self, geometry: str, n_components: int) -> bool:
        """``True`` si la ley sirve a una interfaz de esa geometría con ese
        número de componentes."""
        return (self.GEOMETRY == geometry
                and (self.JUMP_DIM is None or self.JUMP_DIM == n_components))

    @abstractmethod
    def compute_traction(self, jump, state_vars=None):
        """Tracción, tangente algorítmica y estado de prueba.

        Parameters
        ----------
        jump : np.ndarray
            Desplazamiento relativo ``Δu`` en el sistema local de la interfaz,
            en el orden de ``GEOMETRY``.
        state_vars : dict | None
            Estado confirmado del paso anterior; ``None`` en estado virgen.

        Returns
        -------
        (traction, tangent, new_state_vars)
            ``traction`` con la forma de ``jump``; ``tangent`` la matriz
            ``∂t/∂Δu`` en el mismo sistema; ``new_state_vars`` el estado de
            prueba, que el elemento confirma al converger.
        """
