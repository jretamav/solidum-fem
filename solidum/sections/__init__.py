"""Secciones de marco (ADR 0021).

Familia de secciones (``SectionRegistry``, sección YAML ``sections``): leyes
de sección que consumen los marcos Euler-Bernoulli (``Frame2DEuler``,
``Frame2DEulerCorot``, ``Frame3D``) a través de su parámetro ``section``.
Cada archivo del paquete se registra solo con ``@SectionRegistry.register``
(autodescubrimiento, ``solidum/autodiscover.py``).

- ``fiber`` — :class:`FiberSection` y sus primitivas de discretización
  (:class:`RectangularPatch`, :class:`CircularPatch`, :class:`BarLayer`,
  :class:`Fiber`).
"""
from solidum.sections.fiber import (
    BarLayer,
    CircularPatch,
    Fiber,
    FiberSection,
    RectangularPatch,
)

__all__ = [
    "BarLayer",
    "CircularPatch",
    "Fiber",
    "FiberSection",
    "RectangularPatch",
]
