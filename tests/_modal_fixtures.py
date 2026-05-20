"""Fixtures compartidos por los tests de análisis modal y dinámico.

Centraliza la construcción de modelos canónicos que aparecen una y otra
vez en ``tests/test_modal.py`` y ``tests/test_modal_catalog.py``:

- :func:`build_axial_bar_2d` — barra alineada con ``x`` discretizada en
  ``n_elems`` elementos 1D (Truss2D, Truss2DCorot, Cable2DCorot, …);
  ``uy = 0`` global y ``ux = 0`` en el primer nodo.
- :func:`build_axial_bar_3d` — análogo en 3D para Truss3D, Truss3DCorot,
  Cable3DCorot, …; restringe ``uy = uz = 0`` globalmente.
- :func:`build_simply_supported_beam` — viga biapoyada para Frame2DEuler,
  Frame2DEulerCorot y Frame2DTimoshenko; ``ux = uy = 0`` en ambos
  extremos, ``rz`` libre.

Las constantes físicas (acero) coinciden con las usadas históricamente
por los tests modales para no romper expectativas de tolerancia.
"""
from __future__ import annotations

from typing import Any

import numpy as np  # noqa: F401  — disponible para tests que importen este módulo.

from solidum.core.domain import Domain
from solidum.materials.elastic import Elastic1D


# Constantes físicas comunes (acero) — heredadas de los tests originales.
E = 210.0e9
RHO = 7850.0
A_SECTION = 1.0e-4
L_TOTAL = 1.0
N_ELEMS = 20


def build_axial_bar_2d(elem_cls: type, *, n_elems: int = N_ELEMS,
                        **elem_kwargs: Any) -> Domain:
    """Barra 2D empotrada-libre alineada con ``x``, restringida a modo axial.

    Parameters
    ----------
    elem_cls
        Clase del elemento (Truss2D, Truss2DCorot, Cable2DCorot, …).
    n_elems
        Número de elementos a lo largo de ``L_TOTAL``.
    **elem_kwargs
        Argumentos extra del constructor del elemento. ``A=A_SECTION`` es
        el default razonable cuando no se pasa.

    Returns
    -------
    Domain
        Listo para resolver: ``uy = 0`` global, ``ux = 0`` en el primer
        nodo. ``generate_equation_numbers`` ya invocado.
    """
    elem_kwargs.setdefault("A", A_SECTION)
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * L_TOTAL / n_elems, 0.0])
             for i in range(n_elems + 1)]
    mat = Elastic1D(E=E, density=RHO)
    for i in range(n_elems):
        dom.add_element(elem_cls(i + 1, [nodes[i], nodes[i + 1]],
                                  mat, **elem_kwargs))
    nodes[0].fix_dof("ux", 0.0)
    for n in nodes:
        n.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


def build_axial_bar_3d(elem_cls: type, *, n_elems: int = N_ELEMS,
                        **elem_kwargs: Any) -> Domain:
    """Barra 3D empotrada-libre alineada con ``x``, restringida a modo axial.

    Mismo patrón que :func:`build_axial_bar_2d` pero en 3D: ``uy = uz = 0``
    en todos los nodos, ``ux = 0`` en el primero.
    """
    elem_kwargs.setdefault("A", A_SECTION)
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * L_TOTAL / n_elems, 0.0, 0.0])
             for i in range(n_elems + 1)]
    mat = Elastic1D(E=E, density=RHO)
    for i in range(n_elems):
        dom.add_element(elem_cls(i + 1, [nodes[i], nodes[i + 1]],
                                  mat, **elem_kwargs))
    nodes[0].fix_dof("ux", 0.0)
    for n in nodes:
        n.fix_dof("uy", 0.0); n.fix_dof("uz", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom


def build_simply_supported_beam(elem_cls: type, *, n_elems: int = N_ELEMS,
                                  **elem_kwargs: Any) -> Domain:
    """Viga 2D simplemente apoyada para tests modales / dinámicos.

    Parameters
    ----------
    elem_cls
        Clase del elemento (Frame2DEuler, Frame2DEulerCorot,
        Frame2DTimoshenko, …).
    n_elems
        Número de elementos a lo largo de ``L_TOTAL``.
    **elem_kwargs
        Argumentos extra del constructor del elemento. ``A=A_SECTION`` y
        ``I=8.33e-10`` son defaults razonables cuando no se pasan.

    Returns
    -------
    Domain
        Listo para resolver: ``ux = uy = 0`` en ambos apoyos extremos,
        ``rz`` libre. ``generate_equation_numbers`` ya invocado.
    """
    elem_kwargs.setdefault("A", A_SECTION)
    elem_kwargs.setdefault("I", 8.33e-10)
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * L_TOTAL / n_elems, 0.0])
             for i in range(n_elems + 1)]
    mat = Elastic1D(E=E, density=RHO)
    for i in range(n_elems):
        dom.add_element(elem_cls(i + 1, [nodes[i], nodes[i + 1]],
                                  mat, **elem_kwargs))
    nodes[0].fix_dof("ux", 0.0); nodes[0].fix_dof("uy", 0.0)
    nodes[-1].fix_dof("ux", 0.0); nodes[-1].fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom
