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

import numpy as np

GEOMETRIES = ("surface", "line")

_COMPONENT_NAMES = {
    ("surface", 2): ("n", "s"), ("surface", 3): ("n", "s1", "s2"),
    ("line", 2): ("n", "s"), ("line", 3): ("n1", "n2", "s"),
}


def all_component_names() -> tuple[str, ...]:
    """Todos los nombres de componentes locales de la familia, sin repetir
    (los que admite ``component`` en la historia)."""
    return tuple(dict.fromkeys(c for names in _COMPONENT_NAMES.values() for c in names))


def component_names(geometry: str, n_components: int) -> tuple[str, ...]:
    """Nombres de las componentes locales del salto y de la tracción en el
    orden de ``geometry`` (ver el módulo): los que usan la salida VTK
    (``Jump_<c>``, ``Traction_<c>``) y la historia (``component: <c>``)."""
    try:
        return _COMPONENT_NAMES[(geometry, int(n_components))]
    except KeyError:
        raise ValueError(f"interfaz {geometry!r} con {n_components} componentes: se admiten "
                         f"{sorted(_COMPONENT_NAMES)}.") from None


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
    SECANT_UNLOADING : ClassVar[bool]
        ``True`` si la ley descarga por la secante hasta el origen (cohesivas
        de daño): ``½·t·Δu`` es la energía elástica recuperable (ADR 0038).
        Default ``False``.
    IS_UNILATERAL : ClassVar[bool]
        ``True`` si la ley no tiene rigidez en un sentido (contacto abierto):
        la red de seguridad del ADR 0019 avisa en vez de rechazar una
        tangente inicial singular, y los análisis lineal y modal la
        rechazan. Default ``False``.
    ACCEPTS_INITIAL_GAP : ClassVar[bool]
        ``True`` si la ley admite una holgura inicial ``g_0 ≠ 0`` (ADR 0043
        §4): su primera componente normal es entonces la holgura actual
        ``g_0 + Δu_n``. Una ley que no lo declara (cohesivas, adherencia)
        sólo se admite con ``g_0 = 0``, y el elemento lo comprueba al
        construirse. Default ``False``.
    STATE_SCHEMA : ClassVar[dict | None]
        Variables internas ``{nombre: forma}``: ``()`` escalar (la historia
        de la salida la registra), ``(None,)`` vector cuyo tamaño depende de
        la dimensión. ``None``: no declarado.
    STATE_CODES : ClassVar[tuple[str, ...]]
        Variables internas escalares que son un código y no una magnitud (el
        régimen de un punto de contacto: abierto, adherido, deslizando). El
        VTK por celda no las promedia —la media de dos regímenes no es un
        régimen—; se leen en el archivo por puntos de integración y en la
        historia. Default ``()``.

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
    SECANT_UNLOADING: ClassVar[bool] = False
    IS_UNILATERAL: ClassVar[bool] = False
    ACCEPTS_INITIAL_GAP: ClassVar[bool] = False
    STATE_SCHEMA: ClassVar[dict | None] = None
    STATE_CODES: ClassVar[tuple[str, ...]] = ()

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

    def stored_energy(self, jump, state_vars=None) -> float:
        """Energía almacenada por unidad de área ``ψ`` en el punto (ADR 0043 §4).

        ``state_vars`` es el estado que devolvió :meth:`compute_traction` para
        ese mismo ``jump`` (el del paso convergido). Por omisión ``½·t·Δu``,
        que sólo es la energía recuperable si la ley descarga por la secante
        (``SECANT_UNLOADING``); una ley que no lo cumple (fricción,
        plasticidad) o que admite holgura inicial la redefine.
        """
        if not self.SECANT_UNLOADING:
            raise NotImplementedError(
                f"{type(self).__name__} no descarga por la secante y no define "
                f"stored_energy(): su energía almacenada no se conoce.")
        traction, _, _ = self.compute_traction(jump, state_vars)
        return 0.5 * float(np.dot(traction, np.asarray(jump, dtype=np.float64)))
