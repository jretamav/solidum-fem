"""Paquete de elementos de marco/viga 2D.

Históricamente un único archivo, dividido por clase para abaratar las
extensiones futuras (Reglas.md §1). Cada formulación vive en su módulo:

- ``euler`` — :class:`Frame2DEuler` (Euler-Bernoulli, vigas esbeltas).
- ``timoshenko`` — :class:`Frame2DTimoshenko` (deformación por cortante,
  vigas peraltadas).
- ``euler_corot`` — :class:`Frame2DEulerCorot` (Updated Lagrangian, grandes
  desplazamientos y rotaciones rígidas).

Las tres clases son independientes (heredan de :class:`Element`, no entre
sí). Conviven en este paquete por familia temática y porque comparten
unos pocos helpers libres en :mod:`_shared`:
:func:`build_geometry_2d` (longitud + cosenos directores + matriz ``T``
6×6) y los integradores nodales de masa, carga de cuerpo y traducción
``F_local → ElementForces``.

Los reexports a nivel de paquete preservan los imports históricos
``from fenix.elements.frame import Frame2DEuler`` sin alterar al
consumidor. La importación del paquete dispara el
``@ElementRegistry.register`` de cada submódulo como efecto colateral,
idéntico al comportamiento del archivo monolítico previo.
"""
from fenix.elements.frame.euler import Frame2DEuler
from fenix.elements.frame.euler_corot import Frame2DEulerCorot
from fenix.elements.frame.timoshenko import Frame2DTimoshenko

__all__ = [
    "Frame2DEuler",
    "Frame2DEulerCorot",
    "Frame2DTimoshenko",
]
