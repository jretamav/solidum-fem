"""Magnitud por punto de integración de un elemento (ADR 0043 §3).

Cada elemento entrega a su material, en cada punto de integración, una de tres
magnitudes, y lo declara con un atributo de su contrato:

- ``STRAIN_DIM`` (1, 3 o 6): deformación en Voigt (o escalar en las barras) y
  esfuerzo; sólidos, barras, cables y marcos.
- ``FLUX_DIM`` (2 o 3): gradiente de temperatura y flujo; elementos térmicos.
- ``JUMP_DIM`` (2 o 3): salto (desplazamiento relativo) y tracción en ejes
  locales de una interfaz; elementos de interfaz.

La clasificación vive aquí y no en cada consumidor (exportación VTK, historia
de la salida, calidad de malla, herramienta de specs, camino por lotes): antes
cada uno la deducía por presencia de atributos y un elemento nuevo podía pasar
por sólido en unos sitios y no en otros. Prioridad: ``JUMP_DIM``, después
``FLUX_DIM``, después ``STRAIN_DIM``; un atributo vale sólo si es un entero.
"""
from __future__ import annotations

from dataclasses import dataclass

#: Deformación/esfuerzo (``STRAIN_DIM``).
STRAIN = "strain"
#: Gradiente de temperatura/flujo (``FLUX_DIM``).
FLUX = "flux"
#: Salto/tracción de interfaz (``JUMP_DIM``).
JUMP = "jump"

_ATTRIBUTES = (("JUMP_DIM", JUMP), ("FLUX_DIM", FLUX), ("STRAIN_DIM", STRAIN))


@dataclass(frozen=True)
class PointQuantity:
    """Tipo (``STRAIN``, ``FLUX``, ``JUMP`` o ``None``) y número de
    componentes de la magnitud por punto de un elemento."""
    kind: str | None
    dim: int | None

    @property
    def attribute(self) -> str | None:
        """Nombre del atributo del contrato que la declara."""
        return next((a for a, k in _ATTRIBUTES if k == self.kind), None)


_NONE = PointQuantity(None, None)


def point_quantity(element) -> PointQuantity:
    """Magnitud por punto de ``element`` (instancia o clase).

    Un elemento sin ninguno de los tres atributos como entero (p. ej. una
    clase base abstracta) da ``PointQuantity(None, None)``.
    """
    for attr, kind in _ATTRIBUTES:
        value = getattr(element, attr, None)
        if isinstance(value, int) and not isinstance(value, bool):
            return PointQuantity(kind, value)
    return _NONE
