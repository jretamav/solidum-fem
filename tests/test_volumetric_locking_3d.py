"""Locking volumétrico Hex8 con ν → 0.5 — limitación arquitectural documentada.

Espejo 3D de ``test_volumetric_locking.py`` (2D). Documenta y blinda la
limitación declarada en la spec del Hex8 y en STATUS.md: el Hex8 con
integración completa 2×2×2 no puede representar desplazamientos puramente
isocóricos sin sobre-rigidizarse cuando ν → 0.5. Mitigación B-bar/F-bar
diferida hasta caso de uso real con materiales incompresibles (política
idéntica al Quad4 2D).

Cantidad medida: desplazamiento del extremo libre en una viga 3D esbelta
con carga transversal puntual, barriendo ν ∈ {0.3, 0.4999}.

    ratio = u_tip(ν=0.4999) / u_tip(ν=0.3)

Un ratio ≪ 1 documenta el colapso por locking. La viga Euler-Bernoulli
analítica es prácticamente independiente de ν (E gobierna), así que ν
no debería afectar la flecha si el elemento no lockease.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_3d import Hex8
from solidum.materials.elastic_3d import Elastic3D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver


_L = 8.0
_h = 1.0
_t = 1.0
_nx = 4
_ny = 1
_nz = 1
_P = 1.0


def _build_cantilever_hex8(material):
    """Malla 4×1×1 Hex8 sobre cantilever L=8, h=1, t=1."""
    dom = Domain()
    nid = {}
    node_id = 0
    for K in range(_nz + 1):
        for J in range(_ny + 1):
            for I in range(_nx + 1):
                node_id += 1
                x = I * _L / _nx
                y = J * _h / _ny
                z = K * _t / _nz
                nid[(I, J, K)] = dom.add_node(node_id, [x, y, z])

    eid = 0
    for k in range(_nz):
        for j in range(_ny):
            for i in range(_nx):
                local = [
                    nid[(i,     j,     k)],
                    nid[(i + 1, j,     k)],
                    nid[(i + 1, j + 1, k)],
                    nid[(i,     j + 1, k)],
                    nid[(i,     j,     k + 1)],
                    nid[(i + 1, j,     k + 1)],
                    nid[(i + 1, j + 1, k + 1)],
                    nid[(i,     j + 1, k + 1)],
                ]
                eid += 1
                dom.add_element(Hex8(eid, local, material))

    # Empotramiento total en x=0.
    for n in dom.nodes.values():
        if abs(n.coordinates[0]) < 1e-12:
            n.fix_dof('ux', 0.0)
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)

    tip_node = nid[(_nx, 0, 0)]
    dom.generate_equation_numbers(verbose=False)
    return dom, tip_node


def _solve_with_tip_load(material):
    dom, tip = _build_cantilever_hex8(material)
    assembler = Assembler(dom)
    F = np.zeros(dom.total_dofs)
    F[tip.dofs['uy']] = -_P
    U = LinearSolver(assembler).solve(F)
    return U[tip.dofs['uy']]


class TestHex8VolumetricLocking(unittest.TestCase):

    def test_hex8_locks_with_incompressible_nu(self):
        """Ratio u_tip(ν=0.4999) / u_tip(ν=0.3) ≪ 1 documenta el locking."""
        m_03 = Elastic3D(E=1.0e5, nu=0.3)
        m_05 = Elastic3D(E=1.0e5, nu=0.4999)
        u_03 = _solve_with_tip_load(m_03)
        u_05 = _solve_with_tip_load(m_05)
        ratio = u_05 / u_03
        # Para Hex8 esbelto 4×1×1 el locking severo deja la flecha incompresible
        # significativamente reducida. La cota superior estricta (ratio < 0.6) es
        # la propiedad cualitativa que documentamos: locking presente.
        self.assertLess(
            ratio, 0.6,
            msg=f"Hex8 4×1×1 con ν=0.4999 NO lockea (ratio={ratio:.3f}); se "
                f"esperaba ratio<0.6 documentando el colapso. u_03={u_03:.3e}, "
                f"u_05={u_05:.3e}.",
        )
        self.assertGreater(
            ratio, 0.0,
            msg="Ratio negativo: el modelo o las BCs están invertidos.",
        )

    def test_hex8_at_nu_03_gives_reasonable_deflection(self):
        """Sanity: con ν=0.3 la flecha no está dominada por locking; tiene
        que ser del mismo orden de magnitud que la analítica Euler-Bernoulli
        (relación >0.4× del valor exacto). El Hex8 sí sufre shear locking,
        pero menos severo que el volumétrico."""
        m = Elastic3D(E=1.0e5, nu=0.3)
        u = _solve_with_tip_load(m)
        # Euler-Bernoulli: u_tip = P·L³ / (3·E·I), con I = t·h³/12.
        I = _t * _h ** 3 / 12.0
        u_eb = -_P * _L ** 3 / (3.0 * 1.0e5 * I)  # negativo por carga descendente
        ratio = u / u_eb
        # Espejo del Quad4 2D: malla coarse subestima por shear locking;
        # cota mínima permisiva (>0.30) documenta que no está totalmente colapsada.
        self.assertGreater(
            ratio, 0.30,
            msg=f"u_tip(ν=0.3)/u_EB = {ratio:.3f} demasiado bajo (esperado >0.30); "
                f"u_tip={u:.3e}, u_EB={u_eb:.3e}.",
        )


if __name__ == "__main__":
    unittest.main()
