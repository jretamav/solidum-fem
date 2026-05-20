"""Paquete de elementos sólidos 3D (ADR 0012, Etapa 7).

Componentes:

- ``hex8`` — :class:`Hex8` (hexaedro trilineal isoparamétrico).
- ``tet4`` — :class:`Tet4` (tetraedro lineal CST 3D).

Los helpers libres (funciones de forma, derivadas, kinematics Numba,
expansor escalar→bloque 3D) viven en :mod:`_shared`. Convención Voigt 6D
del proyecto fijada en ADR 0012.
"""
from solidum.elements.solid_3d.hex8 import Hex8
from solidum.elements.solid_3d.tet4 import Tet4

__all__ = [
    "Hex8",
    "Tet4",
]
