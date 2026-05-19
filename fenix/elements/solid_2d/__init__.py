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
viven en :mod:`_shared`. Algunos tests del repositorio importan estos
helpers privados directamente (``test_patch_solid_2d``,
``test_higher_order_solid_2d``); el paquete los reexporta para preservar
esos imports sin cambios.
"""
# Helpers privados expuestos para tests existentes (test_patch_solid_2d,
# test_higher_order_solid_2d). Reexports estables; mantenerlos al partir
# evita un re-cableado de tests sin valor.
from fenix.elements.solid_2d._shared import (
    _compute_kinematics,
    _compute_kinematics_tri3,
    _dN_quad8,
    _dN_quad9,
    _dN_tri6,
    _N_quad8,
    _N_quad9,
    _N_tri6,
)
from fenix.elements.solid_2d.quad4 import Quad4
from fenix.elements.solid_2d.quad8 import Quad8
from fenix.elements.solid_2d.quad9 import Quad9
from fenix.elements.solid_2d.tri3 import Tri3
from fenix.elements.solid_2d.tri6 import Tri6

__all__ = [
    "Quad4",
    "Quad8",
    "Quad9",
    "Tri3",
    "Tri6",
]
