"""Paquete de elementos sólidos 3D (ADR 0012, Etapa 7 + sub-etapa A.ter).

Componentes:

- ``hex8`` — :class:`Hex8` (hexaedro trilineal isoparamétrico).
- ``tet4`` — :class:`Tet4` (tetraedro lineal CST 3D).
- ``hex20`` — :class:`Hex20` (hexaedro serendípito 3D de orden 2, sub-etapa A.ter).

Los helpers libres (funciones de forma, derivadas, kinematics Numba,
kinematics genérico de orden superior 3D, expansor escalar→bloque 3D)
viven en :mod:`_shared`. Convención Voigt 6D del proyecto fijada en
ADR 0012.
"""
from solidum.elements.solid_3d.hex8 import Hex8
from solidum.elements.solid_3d.hex20 import Hex20
from solidum.elements.solid_3d.hex27 import Hex27
from solidum.elements.solid_3d.tet4 import Tet4
from solidum.elements.solid_3d.tet10 import Tet10

__all__ = [
    "Hex8",
    "Hex20",
    "Hex27",
    "Tet4",
    "Tet10",
]
