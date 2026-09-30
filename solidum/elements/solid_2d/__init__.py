"""Paquete de elementos sólidos 2D (deformación/esfuerzo plano).

Históricamente un único archivo de 929 líneas, dividido por elemento para
abaratar la incorporación de nuevas formulaciones (Reglas.md §1). Cada
clase vive en su módulo:

- ``quad4`` — :class:`Quad4` (cuadrilátero bilineal isoparamétrico).
- ``tri3`` — :class:`Tri3` (triángulo lineal CST).
- ``quad8`` — :class:`Quad8` (cuadrilátero serendípito).
- ``quad9`` — :class:`Quad9` (cuadrilátero Lagrangiano).
- ``tri6`` — :class:`Tri6` (triángulo cuadrático completo P₂).

Los helpers libres (funciones de forma, derivadas, kinematics Numba para
los lineales, ``_expand_scalar_mass``, ``_quadratic_edge_traction`` y la
base interna ``_HigherOrderSolid2D`` compartida por Quad8, Quad9 y Tri6)
viven en :mod:`_shared`.

API pública de cinemática (ADR 0020, P6): :func:`compute_kinematics_tri3`
(matriz ``B`` y ``det J`` del triángulo lineal) y :func:`compute_integrands`
(contribuciones ``Bᵀ·C·B·dV`` y ``Bᵀ·σ·dV`` de un punto de Gauss, con
``dV = det J · w · espesor``). Las usan Tri3 y Quad4, y un elemento de
usuario las importa desde aquí, no desde :mod:`_shared`, que es interno y
puede reorganizarse.

El paquete reexporta además algunos helpers privados que importan tests del
repositorio (``test_patch_solid_2d``, ``test_higher_order_solid_2d``); no
forman parte del contrato y quedan fuera de ``__all__``.
"""
# Privados: sólo para tests existentes (test_patch_solid_2d,
# test_higher_order_solid_2d). compute_kinematics_tri3 y compute_integrands
# son API pública (ADR 0020, P6) y figuran en __all__.
from solidum.elements.solid_2d._shared import (
    _compute_kinematics,
    compute_integrands,
    compute_kinematics_tri3,
    _dN_quad8,
    _dN_quad9,
    _dN_tri6,
    _N_quad8,
    _N_quad9,
    _N_tri6,
)
from solidum.elements.solid_2d.quad4 import Quad4
from solidum.elements.solid_2d.quad8 import Quad8
from solidum.elements.solid_2d.quad9 import Quad9
from solidum.elements.solid_2d.tri3 import Tri3
from solidum.elements.solid_2d.tri6 import Tri6

__all__ = [
    "Quad4",
    "Quad8",
    "Quad9",
    "Tri3",
    "Tri6",
    "compute_integrands",
    "compute_kinematics_tri3",
]
