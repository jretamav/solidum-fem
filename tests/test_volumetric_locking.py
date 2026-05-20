"""Locking volumétrico en sólidos 2D plane strain con ν → 0.5 — Fase D del plan
de validación: documentación de limitación arquitectural conocida (no bug).

El elemento `Quad4` con integración 2×2 sobre material plane strain
incompresible (ν → 0.5) no puede satisfacer la restricción `tr(ε) = 0`
punto a punto sin sobre-rigidizarse. El resultado es que la respuesta del
elemento a modos de flexión colapsa cuando ν → 0.5: aparece una rigidez
artificial que reduce el desplazamiento del extremo libre muy por debajo
del valor físicamente correcto. El `Quad8`, al disponer de una base
polinómica cuadrática, mitiga sustancialmente el locking aunque no lo
elimina.

Este test **demuestra y blinda** el límite descrito en
`docs/STATUS.md` §"Limitaciones declaradas" (entrada *Locking volumétrico*).
No intenta corregir el locking — mitigaciones B-bar/F-bar están diferidas
hasta que aparezca un caso de uso real con materiales incompresibles
(plasticidad J2 perfectamente plástica, hiperelasticidad, etc.).

Configuración del benchmark:

    Cantilever rectangular L×h con L/h = 8 (esbelto, dominado por flexión).
    Empotramiento total en x=0; carga vertical descendente P en el nodo
    central del borde derecho x=L. Material Elastic2D plane strain con
    E=1e5 y ν barrido en {0.3, 0.4999}. Espesor unitario.

Cantidad medida: desplazamiento vertical del nodo cargado. Se compara
el **ratio** ``u_y(ν=0.4999) / u_y(ν=0.3)`` para cada elemento — un valor
cercano a 1 indicaría ausencia de locking (la elasticidad lineal de viga
es prácticamente independiente de ν), mientras que un valor << 1
documenta el colapso por locking.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4, Quad8
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver


_L = 8.0     # luz de la viga
_h = 1.0     # canto
_t = 1.0     # espesor (plane strain unitario)
_nx = 4      # elementos en x
_ny = 1      # elementos en y (slenderness L/h = 8 con 4×1 ya es suficiente
             # para que el modo de flexión domine).
_P = 1.0     # carga puntual descendente


def _build_cantilever_quad4(material):
    """Malla 4×1 Quad4 sobre cantilever L=8, h=1."""
    dom = Domain()
    grid = {}
    nid = 0
    for J in range(_ny + 1):
        for I in range(_nx + 1):
            nid += 1
            x = I * _L / _nx
            y = J * _h / _ny
            grid[(I, J)] = dom.add_node(nid, [x, y])

    eid = 0
    for jy in range(_ny):
        for ix in range(_nx):
            eid += 1
            c1 = grid[(ix,     jy    )]
            c2 = grid[(ix + 1, jy    )]
            c3 = grid[(ix + 1, jy + 1)]
            c4 = grid[(ix,     jy + 1)]
            dom.add_element(Quad4(eid, [c1, c2, c3, c4], material, thickness=_t))

    for J in range(_ny + 1):
        n = grid[(0, J)]
        n.fix_dof('ux', 0.0); n.fix_dof('uy', 0.0)

    # Nodo central del borde derecho. Con _ny=1 hay dos nodos en x=L (J=0 y J=1);
    # el "central" se toma como J=ny//2 = 0. Para ny par, sería J=ny/2. Aquí
    # cualquiera de los dos sirve para una carga puntual sobre el borde derecho;
    # tomamos J=0 (esquina inferior derecha) consistentemente entre Q4 y Q8.
    tip_node = grid[(_nx, 0)]
    dom.generate_equation_numbers(verbose=False)
    return dom, tip_node


def _build_cantilever_quad8(material):
    """Malla 4×1 Quad8 sobre cantilever L=8, h=1. Grid (2*_nx+1) × (2*_ny+1) sin centers."""
    dom = Domain()
    grid = {}
    nid = 0
    for J in range(2 * _ny + 1):
        for I in range(2 * _nx + 1):
            is_center = (I % 2 == 1) and (J % 2 == 1)
            if is_center:
                continue
            nid += 1
            x = I * _L / (2.0 * _nx)
            y = J * _h / (2.0 * _ny)
            grid[(I, J)] = dom.add_node(nid, [x, y])

    eid = 0
    for jy in range(_ny):
        for ix in range(_nx):
            I0, J0 = 2 * ix, 2 * jy
            c1 = grid[(I0,     J0    )]
            c2 = grid[(I0 + 2, J0    )]
            c3 = grid[(I0 + 2, J0 + 2)]
            c4 = grid[(I0,     J0 + 2)]
            m12 = grid[(I0 + 1, J0    )]
            m23 = grid[(I0 + 2, J0 + 1)]
            m34 = grid[(I0 + 1, J0 + 2)]
            m41 = grid[(I0,     J0 + 1)]
            eid += 1
            dom.add_element(Quad8(eid, [c1, c2, c3, c4, m12, m23, m34, m41],
                                  material, thickness=_t))

    for J in range(2 * _ny + 1):
        if (0, J) in grid:
            n = grid[(0, J)]
            n.fix_dof('ux', 0.0); n.fix_dof('uy', 0.0)

    tip_node = grid[(2 * _nx, 0)]
    dom.generate_equation_numbers(verbose=False)
    return dom, tip_node


def _tip_deflection(builder, nu: float) -> float:
    """Carga el cantilever con `nu` y devuelve |u_y| en el nodo del extremo."""
    mat = Elastic2D(E=1.0e5, nu=nu, hypothesis='plane_strain')
    dom, tip = builder(mat)
    F = np.zeros(dom.total_dofs)
    F[tip.dofs['uy']] = -_P
    U = LinearSolver(Assembler(dom)).solve(F)
    return abs(float(U[tip.dofs['uy']]))


class TestVolumetricLocking(unittest.TestCase):
    """Documenta el locking volumétrico de Quad4 plane strain ν → 0.5 y la
    mitigación parcial del Quad8, vía ratios de desplazamiento del extremo
    de un cantilever esbelto."""

    def test_quad4_locks_severely(self):
        """Quad4 plane strain pierde >50% de deflexión al pasar ν=0.3 → 0.4999.

        La elasticidad de viga esbelta predice una dependencia suave en ν
        (el módulo efectivo plane strain cambia por factor 1/(1-ν²), del
        orden ~30% en este rango). Un colapso del ratio por debajo de 0.5
        sólo se explica por locking volumétrico del Q4.
        """
        u_03 = _tip_deflection(_build_cantilever_quad4, nu=0.3)
        u_05 = _tip_deflection(_build_cantilever_quad4, nu=0.4999)
        ratio = u_05 / u_03
        self.assertLess(ratio, 0.5,
            f"Q4 ratio u_y(ν=0.4999)/u_y(ν=0.3) = {ratio:.3f} >= 0.5; "
            f"el Quad4 no exhibe el locking severo esperado (u_03={u_03:.4e}, "
            f"u_05={u_05:.4e}). Si esto falla, la formulación del Quad4 ha "
            "cambiado — revisar si se introdujo alguna mitigación.")

    def test_quad8_locks_less_than_quad4(self):
        """Quad8 plane strain lockea menos que Quad4 en el mismo régimen ν → 0.5.

        La base cuadrática del Q8 dispone de modos internos que pueden
        aproximar la restricción de incompresibilidad sin sobre-rigidizar.
        Comparación relativa entre ambos elementos (mismo benchmark) en
        lugar de umbrales absolutos: ``ratio_Q8 > ratio_Q4`` es la métrica
        físicamente significativa — el Quad8 lockea, pero menos.
        """
        ratio_q4 = (_tip_deflection(_build_cantilever_quad4, nu=0.4999) /
                    _tip_deflection(_build_cantilever_quad4, nu=0.3))
        ratio_q8 = (_tip_deflection(_build_cantilever_quad8, nu=0.4999) /
                    _tip_deflection(_build_cantilever_quad8, nu=0.3))
        self.assertGreater(ratio_q8, ratio_q4,
            f"Q8 lockea más o igual que Q4: ratio_Q8={ratio_q8:.3f}, "
            f"ratio_Q4={ratio_q4:.3f}. La mitigación esperada del Q8 no "
            "aparece — revisar formulación.")
        # Margen mínimo: el Q8 debe preservar al menos 1.5× más respuesta
        # que el Q4 en el límite incompresible (criterio amplio, robusto
        # a cambios de malla).
        self.assertGreater(ratio_q8, 1.5 * ratio_q4,
            f"Q8 sólo mejora marginalmente sobre Q4: ratio_Q8={ratio_q8:.3f}, "
            f"ratio_Q4={ratio_q4:.3f}; mitigación esperada >1.5×.")

    def test_quad4_underestimates_quad8_in_incompressible_limit(self):
        """Con ν=0.4999 el Quad4 da una deflexión muy inferior al Quad8.

        Hace explícita la diferencia entre los dos elementos en el régimen
        problemático: el Q4 colapsa, el Q8 mantiene rango físico. Si esta
        relación se invierte o desaparece, hay regresión.
        """
        u_q4 = _tip_deflection(_build_cantilever_quad4, nu=0.4999)
        u_q8 = _tip_deflection(_build_cantilever_quad8, nu=0.4999)
        self.assertLess(u_q4 / u_q8, 0.5,
            f"En ν=0.4999: u_Q4={u_q4:.4e}, u_Q8={u_q8:.4e}; el ratio "
            f"{u_q4 / u_q8:.3f} no muestra el locking diferencial del Q4.")


if __name__ == '__main__':
    unittest.main()
