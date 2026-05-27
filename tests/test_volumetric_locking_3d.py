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
from solidum.elements.solid_3d import Hex8, Hex20, Hex27
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


def _build_cantilever_hex20(material):
    """Malla nx×1×1 de Hex20 (4 elementos) sobre cantilever L=8, h=1, t=1.

    Espejo de ``_build_cantilever_hex8`` pero con elementos cuadráticos
    serendípitos. Para Hex20 se necesitan 2*nx+1 índices por dirección
    (vértices en pares, medios en impares) y sólo se allocan posiciones
    con como mucho un índice impar (vértice + arista, sin centros de
    cara ni interior).
    """
    dom = Domain()
    nIx, nIy, nIz = 2 * _nx + 1, 2 * _ny + 1, 2 * _nz + 1
    nid = {}
    node_id = 0
    for K in range(nIz):
        for J in range(nIy):
            for I in range(nIx):
                n_odd = (I % 2) + (J % 2) + (K % 2)
                if n_odd <= 1:
                    node_id += 1
                    x = I * _L / (2 * _nx)
                    y = J * _h / (2 * _ny)
                    z = K * _t / (2 * _nz)
                    nid[(I, J, K)] = dom.add_node(node_id, [x, y, z])

    eid = 0
    for k in range(_nz):
        for j in range(_ny):
            for i in range(_nx):
                ic, jc, kc = 2 * i, 2 * j, 2 * k
                vertices = [
                    nid[(ic,     jc,     kc)],
                    nid[(ic + 2, jc,     kc)],
                    nid[(ic + 2, jc + 2, kc)],
                    nid[(ic,     jc + 2, kc)],
                    nid[(ic,     jc,     kc + 2)],
                    nid[(ic + 2, jc,     kc + 2)],
                    nid[(ic + 2, jc + 2, kc + 2)],
                    nid[(ic,     jc + 2, kc + 2)],
                ]
                midsides = [
                    nid[(ic + 1, jc,     kc)],       # 8
                    nid[(ic + 2, jc + 1, kc)],       # 9
                    nid[(ic + 1, jc + 2, kc)],       # 10
                    nid[(ic,     jc + 1, kc)],       # 11
                    nid[(ic + 1, jc,     kc + 2)],   # 12
                    nid[(ic + 2, jc + 1, kc + 2)],   # 13
                    nid[(ic + 1, jc + 2, kc + 2)],   # 14
                    nid[(ic,     jc + 1, kc + 2)],   # 15
                    nid[(ic,     jc,     kc + 1)],   # 16
                    nid[(ic + 2, jc,     kc + 1)],   # 17
                    nid[(ic + 2, jc + 2, kc + 1)],   # 18
                    nid[(ic,     jc + 2, kc + 1)],   # 19
                ]
                eid += 1
                dom.add_element(Hex20(eid, vertices + midsides, material))

    for n in dom.nodes.values():
        if abs(n.coordinates[0]) < 1e-12:
            n.fix_dof('ux', 0.0)
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)

    tip_node = nid[(2 * _nx, 0, 0)]
    dom.generate_equation_numbers(verbose=False)
    return dom, tip_node


def _solve_with_tip_load_hex20(material):
    dom, tip = _build_cantilever_hex20(material)
    assembler = Assembler(dom)
    F = np.zeros(dom.total_dofs)
    F[tip.dofs['uy']] = -_P
    U = LinearSolver(assembler).solve(F)
    return U[tip.dofs['uy']]


class TestHex20VolumetricLocking(unittest.TestCase):
    """Locking volumétrico Hex20 — limitación atenuada respecto al Hex8.

    El Hex20 mantiene locking volumétrico (el espacio serendípito
    triquadrático tampoco contiene desplazamientos puramente isocóricos
    suficientes), pero significativamente atenuado: para la misma viga
    esbelta 4×1×1 que en el Hex8, el ratio u(ν=0.4999)/u(ν=0.3) sube de
    < 0.6 (Hex8) a ≈ 0.8 (Hex20). La mitigación es parcial, **no
    completa** — requeriría B-bar/F-bar para eliminarla por completo
    (fuera de alcance, política idéntica al 2D).
    """

    def test_hex20_locks_less_than_hex8(self):
        """El locking del Hex20 es atenuado respecto al Hex8."""
        m_03 = Elastic3D(E=1.0e5, nu=0.3)
        m_05 = Elastic3D(E=1.0e5, nu=0.4999)
        u_03_h20 = _solve_with_tip_load_hex20(m_03)
        u_05_h20 = _solve_with_tip_load_hex20(m_05)
        ratio_h20 = u_05_h20 / u_03_h20

        u_03_h8 = _solve_with_tip_load(m_03)
        u_05_h8 = _solve_with_tip_load(m_05)
        ratio_h8 = u_05_h8 / u_03_h8

        # Documentado: Hex20 mantiene > 70% de la flecha al pasar
        # ν=0.3→0.4999, mientras que Hex8 cae bajo 60%.
        self.assertGreater(
            ratio_h20, 0.70,
            msg=f"Hex20 con ν=0.4999 lockea más de lo esperado: "
                f"ratio={ratio_h20:.3f}; se esperaba > 0.70.",
        )
        # La mitigación respecto al Hex8 debe ser clara.
        self.assertGreater(
            ratio_h20, ratio_h8,
            msg=f"Hex20 no atenúa el locking respecto al Hex8: "
                f"ratio_h20={ratio_h20:.3f}, ratio_h8={ratio_h8:.3f}.",
        )

    def test_hex20_still_locks_some(self):
        """Aún con Hex20 el locking volumétrico no desaparece — sigue
        siendo necesaria B-bar/F-bar para eliminarlo completamente.

        Esta cota superior documenta la **limitación que permanece**:
        no podemos pretender que el Hex20 está libre de locking
        volumétrico — es una mitigación parcial, no una cura.
        """
        m_03 = Elastic3D(E=1.0e5, nu=0.3)
        m_05 = Elastic3D(E=1.0e5, nu=0.4999)
        u_03 = _solve_with_tip_load_hex20(m_03)
        u_05 = _solve_with_tip_load_hex20(m_05)
        ratio = u_05 / u_03
        # Si esto fallara hacia arriba (ratio > 0.95), significaría que
        # el Hex20 está libre de locking volumétrico — algo que iría en
        # contra de la teoría establecida (Cook §13.6). Probablemente
        # señalaría un bug que está enmascarando la incompresibilidad.
        self.assertLess(
            ratio, 0.95,
            msg=f"Hex20 con ν=0.4999 NO muestra locking volumétrico "
                f"(ratio={ratio:.3f} ≥ 0.95); revisar formulación, podría "
                f"haber un bug que enmascara la incompresibilidad.",
        )


def _build_cantilever_hex27(material):
    """Malla nx×1×1 de Hex27 sobre cantilever L=8, h=1, t=1.

    Espejo de ``_build_cantilever_hex20`` pero allocando **todas** las
    posiciones (incluyendo face centers y body centers).
    """
    dom = Domain()
    nIx, nIy, nIz = 2 * _nx + 1, 2 * _ny + 1, 2 * _nz + 1
    nid = {}
    node_id = 0
    for K in range(nIz):
        for J in range(nIy):
            for I in range(nIx):
                node_id += 1
                x = I * _L / (2 * _nx)
                y = J * _h / (2 * _ny)
                z = K * _t / (2 * _nz)
                nid[(I, J, K)] = dom.add_node(node_id, [x, y, z])

    eid = 0
    for k in range(_nz):
        for j in range(_ny):
            for i in range(_nx):
                ic, jc, kc = 2 * i, 2 * j, 2 * k
                vertices = [
                    nid[(ic,     jc,     kc)],
                    nid[(ic + 2, jc,     kc)],
                    nid[(ic + 2, jc + 2, kc)],
                    nid[(ic,     jc + 2, kc)],
                    nid[(ic,     jc,     kc + 2)],
                    nid[(ic + 2, jc,     kc + 2)],
                    nid[(ic + 2, jc + 2, kc + 2)],
                    nid[(ic,     jc + 2, kc + 2)],
                ]
                midsides = [
                    nid[(ic + 1, jc,     kc)],
                    nid[(ic + 2, jc + 1, kc)],
                    nid[(ic + 1, jc + 2, kc)],
                    nid[(ic,     jc + 1, kc)],
                    nid[(ic + 1, jc,     kc + 2)],
                    nid[(ic + 2, jc + 1, kc + 2)],
                    nid[(ic + 1, jc + 2, kc + 2)],
                    nid[(ic,     jc + 1, kc + 2)],
                    nid[(ic,     jc,     kc + 1)],
                    nid[(ic + 2, jc,     kc + 1)],
                    nid[(ic + 2, jc + 2, kc + 1)],
                    nid[(ic,     jc + 2, kc + 1)],
                ]
                face_centers = [
                    nid[(ic,     jc + 1, kc + 1)],
                    nid[(ic + 2, jc + 1, kc + 1)],
                    nid[(ic + 1, jc,     kc + 1)],
                    nid[(ic + 1, jc + 2, kc + 1)],
                    nid[(ic + 1, jc + 1, kc)],
                    nid[(ic + 1, jc + 1, kc + 2)],
                ]
                body_center = [nid[(ic + 1, jc + 1, kc + 1)]]
                eid += 1
                dom.add_element(
                    Hex27(eid, vertices + midsides + face_centers + body_center,
                          material)
                )

    for n in dom.nodes.values():
        if abs(n.coordinates[0]) < 1e-12:
            n.fix_dof('ux', 0.0)
            n.fix_dof('uy', 0.0)
            n.fix_dof('uz', 0.0)

    tip_node = nid[(2 * _nx, 0, 0)]
    dom.generate_equation_numbers(verbose=False)
    return dom, tip_node


def _solve_with_tip_load_hex27(material):
    dom, tip = _build_cantilever_hex27(material)
    assembler = Assembler(dom)
    F = np.zeros(dom.total_dofs)
    F[tip.dofs['uy']] = -_P
    U = LinearSolver(assembler).solve(F)
    return U[tip.dofs['uy']]


class TestHex27VolumetricLocking(unittest.TestCase):
    """Locking volumétrico Hex27 — limitación comparable al Hex20.

    El espacio Lagrangiano completo añade flexibilidad pero no elimina
    la incompresibilidad (Cook §13.6). Esperamos comportamiento similar
    al Hex20 (ratio u(ν=0.4999)/u(ν=0.3) ≈ 0.8), mucho mejor que el Hex8
    (< 0.6) pero todavía con locking apreciable.
    """

    def test_hex27_locks_less_than_hex8(self):
        m_03 = Elastic3D(E=1.0e5, nu=0.3)
        m_05 = Elastic3D(E=1.0e5, nu=0.4999)
        u_03_h27 = _solve_with_tip_load_hex27(m_03)
        u_05_h27 = _solve_with_tip_load_hex27(m_05)
        ratio_h27 = u_05_h27 / u_03_h27

        u_03_h8 = _solve_with_tip_load(m_03)
        u_05_h8 = _solve_with_tip_load(m_05)
        ratio_h8 = u_05_h8 / u_03_h8

        self.assertGreater(
            ratio_h27, 0.70,
            msg=f"Hex27 con ν=0.4999 lockea más de lo esperado: "
                f"ratio={ratio_h27:.3f}; se esperaba > 0.70.",
        )
        self.assertGreater(
            ratio_h27, ratio_h8,
            msg=f"Hex27 no atenúa el locking respecto al Hex8: "
                f"ratio_h27={ratio_h27:.3f}, ratio_h8={ratio_h8:.3f}.",
        )

    def test_hex27_still_locks_some(self):
        """El espacio Q_2 lagrangiano completo no elimina el locking.

        Si esto fallara hacia arriba (ratio > 0.95), indicaría que el
        espacio Lagrangiano completo está eliminando la incompresibilidad,
        lo cual contradice Cook §13.6 — probablemente sería un bug.
        """
        m_03 = Elastic3D(E=1.0e5, nu=0.3)
        m_05 = Elastic3D(E=1.0e5, nu=0.4999)
        u_03 = _solve_with_tip_load_hex27(m_03)
        u_05 = _solve_with_tip_load_hex27(m_05)
        ratio = u_05 / u_03
        self.assertLess(
            ratio, 0.95,
            msg=f"Hex27 con ν=0.4999 NO muestra locking volumétrico "
                f"(ratio={ratio:.3f} ≥ 0.95); revisar formulación.",
        )


if __name__ == "__main__":
    unittest.main()
