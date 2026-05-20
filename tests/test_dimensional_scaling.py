"""Tests de escalamiento dimensional (regresión H-2.2).

La auditoría 2026-05-18 destacó que ~27 ocurrencias de ``thickness=1.0`` en
los tests dejaban un agujero: el factor ``thickness`` dentro de
``_compute_integrands`` se cancela contra 1.0 y una omisión accidental
pasaría inadvertida. Memoria del proyecto:
``feedback_tests_con_coeficientes_unitarios.md`` — *los coeficientes
unitarios esconden bugs dimensionales*. El bug histórico del ``thickness``
en ``CST_Embedded2D`` se descubrió precisamente por un benchmark con
``thickness=0.1``.

Este archivo añade tests fail-fast que comparan ``K_e`` con varios valores
de ``thickness`` y exigen el escalamiento lineal ``K_e ∝ thickness``.
Cubre los cinco sólidos 2D del catálogo (``Tri3``, ``Quad4``, ``Tri6``,
``Quad8``, ``Quad9``).
"""
from __future__ import annotations

import unittest

import numpy as np

import solidum  # autodiscover

from solidum.core.domain import Domain
from solidum.elements.solid_2d.quad4 import Quad4
from solidum.elements.solid_2d.quad8 import Quad8
from solidum.elements.solid_2d.quad9 import Quad9
from solidum.elements.solid_2d.tri3 import Tri3
from solidum.elements.solid_2d.tri6 import Tri6
from solidum.materials.elastic_2d import Elastic2D


# Geometría unitaria, mismas reglas en todos los elementos.
_QUAD4 = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]
_TRI3 = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
_TRI6 = [
    [0.0, 0.0], [1.0, 0.0], [0.0, 1.0],
    [0.5, 0.0], [0.5, 0.5], [0.0, 0.5],
]
_QUAD8 = [
    [0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0],
    [0.5, 0.0], [1.0, 0.5], [0.5, 1.0], [0.0, 0.5],
]
_QUAD9 = _QUAD8 + [[0.5, 0.5]]

# (ElementClass, coords). Cubre los cinco sólidos 2D del catálogo.
_SOLID_2D_CASES = [
    (Tri3, _TRI3),
    (Quad4, _QUAD4),
    (Tri6, _TRI6),
    (Quad8, _QUAD8),
    (Quad9, _QUAD9),
]


def _build_single_element(elem_cls, coords, thickness: float):
    """Dominio con un elemento aislado y el primer nodo fijado para que
    ``generate_equation_numbers`` no se queje de DOFs colgantes.
    """
    dom = Domain()
    nodes = [dom.add_node(i + 1, c) for i, c in enumerate(coords)]
    mat = Elastic2D(E=1000.0, nu=0.25, hypothesis="plane_strain")
    elem = elem_cls(1, nodes, mat, thickness=thickness)
    dom.add_element(elem)
    nodes[0].fix_dof("ux", 0.0)
    nodes[0].fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return elem


class TestStiffnessScalesWithThickness(unittest.TestCase):
    """``K_e ∝ thickness`` para los 5 sólidos 2D — fail-fast contra omisión
    accidental del factor en ``_compute_integrands`` o equivalente.
    """

    def test_K_e_scales_linearly_with_thickness(self):
        for elem_cls, coords in _SOLID_2D_CASES:
            with self.subTest(elem=elem_cls.__name__):
                K_unit = _build_single_element(
                    elem_cls, coords, thickness=1.0,
                ).compute_global_stiffness()
                K_thick = _build_single_element(
                    elem_cls, coords, thickness=2.5,
                ).compute_global_stiffness()
                K_thin = _build_single_element(
                    elem_cls, coords, thickness=0.3,
                ).compute_global_stiffness()

                # Norma de Frobenius como invariante intrínseco —
                # captura escalamiento uniforme sin sensibilidad al ruido
                # numérico en los ceros estructurales del array. Una
                # omisión del factor en _compute_integrands rompería
                # la relación lineal de la norma.
                norm_unit = np.linalg.norm(K_unit)
                self.assertGreater(norm_unit, 0.0)
                np.testing.assert_allclose(
                    np.linalg.norm(K_thick), 2.5 * norm_unit, rtol=1e-12,
                    err_msg=f"{elem_cls.__name__}: ‖K_e‖ no escala con thickness=2.5",
                )
                np.testing.assert_allclose(
                    np.linalg.norm(K_thin), 0.3 * norm_unit, rtol=1e-12,
                    err_msg=f"{elem_cls.__name__}: ‖K_e‖ no escala con thickness=0.3",
                )


if __name__ == "__main__":
    unittest.main()
