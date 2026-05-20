"""Tests de mass lumping — fase 2 del ADR 0009.

Cubre:

1. **Helper ``lump_hrz``** — unidad pura del algoritmo HRZ
   (Hinton-Rock-Zienkiewicz 1976).
2. **Matrices elementales con ``lumping="lumped"``**:
   - Diagonalidad estricta (los off-diagonals son cero).
   - Masa total preservada: ``Σ diag traslacional == D · ρ·V_e``.
   - Entradas positivas (M no singular, indispensable para modal/dinámico).
   - Invariancia a la orientación del elemento (rotación rígida no
     cambia la masa total ni rompe la diagonalidad en frames 2D).
3. **Recuperación modal con lumped** — frecuencia fundamental de la
   barra axial empotrada-libre dentro de la tolerancia esperada del
   lumping ``O((L/N)²)``.
"""
import math
import unittest

import numpy as np

import solidum  # autodiscover
from solidum.core.domain import Domain
from solidum.elements.cable import Cable2DCorot, Cable3DCorot
from solidum.elements.frame import (
    Frame2DEuler,
    Frame2DEulerCorot,
    Frame2DTimoshenko,
)
from solidum.elements.frame3d import Frame3D
from solidum.elements.solid_2d import Quad4, Quad8, Quad9, Tri3, Tri6
from solidum.elements.truss import Truss2D, Truss3D
from solidum.entry import run_modal
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.mass_lumping import lump_hrz
from solidum.math.solvers import ModalSolver

from _modal_fixtures import (
    A_SECTION,
    E,
    L_TOTAL,
    N_ELEMS,
    RHO,
    build_axial_bar_2d,
)


# ---------------------------------------------------------------------------
# Helper lump_hrz — tests unitarios
# ---------------------------------------------------------------------------


class TestLumpHrzHelper(unittest.TestCase):
    """``lump_hrz`` debe satisfacer el contrato HRZ canónico (Bathe §9.2.4)."""

    def test_diagonal_output(self):
        """La salida es diagonal: off-diagonals son cero."""
        M_c = np.array([[2.0, 0.5, 0.3],
                        [0.5, 2.0, 0.5],
                        [0.3, 0.5, 2.0]])
        M_l = lump_hrz(M_c, total_mass=1.0, n_translational_dirs=1)
        np.testing.assert_allclose(M_l - np.diag(np.diag(M_l)), 0.0)

    def test_preserves_total_mass(self):
        """Σ diagonal traslacional == D · m_total tras lumping HRZ."""
        # 4×4 con diagonal arbitraria; todos traslacionales (D=2 direcciones).
        M_c = np.diag([2.0, 2.0, 2.0, 2.0])
        m_total = 5.0
        M_l = lump_hrz(M_c, total_mass=m_total, n_translational_dirs=2)
        self.assertAlmostEqual(np.sum(np.diag(M_l)), 2 * m_total)

    def test_translational_subset(self):
        """Cuando hay DOFs rotacionales, ``translational_dofs`` los excluye
        del cómputo del factor α — pero α se aplica también a los
        rotacionales (su diagonal se escala, no se elimina).
        """
        # Frame ficticio: 2 trasl + 1 rot + 2 trasl + 1 rot
        diag = np.array([0.3, 0.3, 0.05, 0.3, 0.3, 0.05])
        M_c = np.diag(diag)
        m_total = 1.0
        trasl = np.array([0, 1, 3, 4])
        M_l = lump_hrz(M_c, total_mass=m_total, n_translational_dirs=2,
                        translational_dofs=trasl)
        # α = D · m_total / Σ_trasl diag = 2·1 / (0.3·4) = 1/0.6 = 5/3.
        alpha_expected = 2.0 * m_total / (0.3 * 4)
        np.testing.assert_allclose(np.diag(M_l), alpha_expected * diag)
        # Y la masa total traslacional == D · m_total.
        self.assertAlmostEqual(np.sum(np.diag(M_l)[trasl]), 2.0 * m_total)

    def test_raises_on_zero_diagonal(self):
        """Suma diagonal traslacional ≤ 0 dispara ``ValueError`` con mensaje
        diagnóstico (modo nulo espurio o cuadratura insuficiente).
        """
        M_c = np.zeros((3, 3))
        with self.assertRaisesRegex(ValueError, "no es positiva"):
            lump_hrz(M_c, total_mass=1.0, n_translational_dirs=1)


# ---------------------------------------------------------------------------
# Truss / Cable — fórmula nodal directa
# ---------------------------------------------------------------------------


def _build_truss_segment(elem_cls, *, ndim, length=2.0, A=1.5e-4,
                          rho=7850.0, angle=0.0):
    """Construye un solo elemento barra orientado a ``angle`` (en el plano xy
    para ndim=2; en el plano xy del 3D para ndim=3). Devuelve el elemento."""
    dom = Domain()
    if ndim == 2:
        x1, y1 = 0.0, 0.0
        x2, y2 = length * math.cos(angle), length * math.sin(angle)
        n1 = dom.add_node(1, [x1, y1])
        n2 = dom.add_node(2, [x2, y2])
    else:
        x1, y1, z1 = 0.0, 0.0, 0.0
        x2 = length * math.cos(angle)
        y2 = length * math.sin(angle)
        z2 = 0.0
        n1 = dom.add_node(1, [x1, y1, z1])
        n2 = dom.add_node(2, [x2, y2, z2])
    mat = Elastic1D(E=210e9, density=rho)
    elem = elem_cls(1, [n1, n2], mat, A=A)
    dom.add_element(elem)
    return elem, length


class TestLumpedTrussCable(unittest.TestCase):
    """Truss/Cable: ``lumping="lumped"`` = ``ρAL/2`` por DOF traslacional."""

    def _check_truss_like(self, elem_cls, ndim):
        for angle in (0.0, math.pi / 6, math.pi / 3, math.pi / 2):
            elem, length = _build_truss_segment(elem_cls, ndim=ndim,
                                                  angle=angle)
            M_l = elem.compute_mass_matrix(lumping="lumped")
            n_dofs = M_l.shape[0]
            m_e = elem.material.density * elem.A * length
            # Diagonalidad estricta.
            np.testing.assert_allclose(M_l - np.diag(np.diag(M_l)), 0.0,
                                         atol=1e-14)
            # Cada entrada = ρAL/2.
            np.testing.assert_allclose(np.diag(M_l), 0.5 * m_e,
                                         err_msg=f"{elem_cls.__name__} @ {angle}")
            # Masa total por dirección = ρAL.
            n_dirs = ndim
            self.assertAlmostEqual(np.sum(np.diag(M_l)) / n_dirs, m_e)

    def test_truss2d(self):
        self._check_truss_like(Truss2D, ndim=2)

    def test_truss3d(self):
        self._check_truss_like(Truss3D, ndim=3)

    def test_cable2d(self):
        self._check_truss_like(Cable2DCorot, ndim=2)

    def test_cable3d(self):
        self._check_truss_like(Cable3DCorot, ndim=3)


# ---------------------------------------------------------------------------
# Frame 2D — Euler, Timoshenko, EulerCorot
# ---------------------------------------------------------------------------


def _build_frame2d_segment(elem_cls, *, length=2.0, A=1.5e-4, I=8.33e-10,
                            rho=7850.0, angle=0.0):
    dom = Domain()
    x2, y2 = length * math.cos(angle), length * math.sin(angle)
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [x2, y2])
    mat = Elastic1D(E=210e9, density=rho)
    elem = elem_cls(1, [n1, n2], mat, A=A, I=I)
    dom.add_element(elem)
    return elem, length


class TestLumpedFrame2D(unittest.TestCase):
    """Frame 2D: ``lumping="lumped"`` = ``ρAL/2`` traslacional + ``ρIL/2``
    rotacional por DOF nodal. Diagonal en global para cualquier orientación.
    """

    def _check_frame2d(self, elem_cls):
        for angle in (0.0, math.pi / 6, math.pi / 4, math.pi / 3, math.pi / 2):
            elem, length = _build_frame2d_segment(elem_cls, angle=angle)
            M_l = elem.compute_mass_matrix(lumping="lumped")
            m_e = elem.material.density * elem.A * length
            m_r_total = elem.material.density * elem.I * length
            # Diagonalidad estricta tras rotación.
            off = M_l - np.diag(np.diag(M_l))
            self.assertLess(np.max(np.abs(off)), 1e-12,
                             f"{elem_cls.__name__} no diagonal @ angle={angle}")
            # Masa traslacional total preservada (por dirección).
            diag = np.diag(M_l)
            trasl_idx = [0, 1, 3, 4]  # ux1, uy1, ux2, uy2
            rot_idx = [2, 5]            # θz1, θz2
            self.assertAlmostEqual(
                np.sum(diag[trasl_idx]) / 2.0, m_e,
                msg=f"{elem_cls.__name__}: masa traslacional ≠ ρAL "
                    f"@ angle={angle}"
            )
            # Inercia rotacional total = ρIL.
            self.assertAlmostEqual(
                np.sum(diag[rot_idx]), m_r_total,
                msg=f"{elem_cls.__name__}: inercia rot ≠ ρIL @ angle={angle}"
            )
            # Entradas positivas.
            self.assertTrue(np.all(diag > 0.0))

    def test_frame2d_euler(self):
        self._check_frame2d(Frame2DEuler)

    def test_frame2d_timoshenko(self):
        # Timoshenko necesita E, ν: ajusto el constructor.
        for angle in (0.0, math.pi / 4):
            dom = Domain()
            n1 = dom.add_node(1, [0.0, 0.0])
            n2 = dom.add_node(2, [2.0 * math.cos(angle), 2.0 * math.sin(angle)])
            mat = Elastic1D(E=210e9, density=7850.0)
            elem = Frame2DTimoshenko(1, [n1, n2], mat,
                                       A=1.5e-4, I=8.33e-10, As=1.0e-4)
            dom.add_element(elem)
            M_l = elem.compute_mass_matrix(lumping="lumped")
            m_e = 7850.0 * 1.5e-4 * 2.0
            # Diagonalidad.
            off = M_l - np.diag(np.diag(M_l))
            self.assertLess(np.max(np.abs(off)), 1e-12)
            # Masa traslacional preservada.
            self.assertAlmostEqual(
                (M_l[0, 0] + M_l[1, 1] + M_l[3, 3] + M_l[4, 4]) / 2.0, m_e
            )

    def test_frame2d_euler_corot(self):
        self._check_frame2d(Frame2DEulerCorot)


# ---------------------------------------------------------------------------
# Frame 3D — Iy=Iz da diagonal en global; Iy≠Iz da bloque-diagonal
# ---------------------------------------------------------------------------


class TestLumpedFrame3D(unittest.TestCase):
    """Frame3D lumped: diagonal en local **siempre**; en global la
    diagonalidad sólo se preserva cuando el eje del elemento coincide con
    un eje global (la rotación T no mezcla DOFs rotacionales locales con
    globales). En orientación oblicua queda **bloque-diagonal** por nodo:
    bloque traslacional 3×3 diagonal (m_t·I_3 es invariante a SO(3)) +
    bloque rotacional 3×3 lleno (porque ``ρ·Jp·L`` ≠ ``ρ·Iy·L`` ≠
    ``ρ·Iz·L`` en general — la torsional usa el momento polar y las dos
    flexionales el geométrico, valores siempre distintos para sección
    real). Es la limitación estándar del lumped en frames 3D (Cook-
    Malkus-Plesha §11.4).
    """

    def _build(self, *, Iy, Iz, A=1.5e-4, length=2.0, axis="x"):
        dom = Domain()
        if axis == "x":
            n1 = dom.add_node(1, [0.0, 0.0, 0.0])
            n2 = dom.add_node(2, [length, 0.0, 0.0])
        elif axis == "skew":
            n1 = dom.add_node(1, [0.0, 0.0, 0.0])
            n2 = dom.add_node(2, [length / math.sqrt(3),
                                    length / math.sqrt(3),
                                    length / math.sqrt(3)])
        mat = Elastic1D(E=210e9, density=7850.0)
        elem = Frame3D(1, [n1, n2], mat,
                         A=A, Iy=Iy, Iz=Iz, J=Iy + Iz)
        dom.add_element(elem)
        return elem, length

    def test_diagonal_when_aligned_with_global_axis(self):
        """Eje X global: M_lumped es completamente diagonal."""
        for Iy, Iz in ((8.33e-10, 8.33e-10), (8.33e-10, 2.5e-9)):
            elem, _ = self._build(Iy=Iy, Iz=Iz, axis="x")
            M_l = elem.compute_mass_matrix(lumping="lumped")
            off = M_l - np.diag(np.diag(M_l))
            self.assertLess(np.max(np.abs(off)), 1e-12,
                             f"Iy={Iy}, Iz={Iz}, axis=x: esperaba diagonal.")

    def test_block_diagonal_when_skew(self):
        """Orientación oblicua: traslacional 3×3 diagonal, rotacional 3×3
        lleno, sin acoplamiento entre nodos distintos.
        """
        elem, length = self._build(Iy=8.33e-10, Iz=2.5e-9, axis="skew")
        M_l = elem.compute_mass_matrix(lumping="lumped")
        # Bloque cruzado nodo i ↔ nodo j debe ser cero.
        cross = M_l[0:6, 6:12]
        self.assertLess(np.max(np.abs(cross)), 1e-12,
                         "Lumped 3D no debe acoplar nodos distintos.")
        # Bloque traslacional diagonal (3×3 cada nodo).
        for n_off in (0, 6):
            trans_block = M_l[n_off:n_off + 3, n_off:n_off + 3]
            off = trans_block - np.diag(np.diag(trans_block))
            self.assertLess(np.max(np.abs(off)), 1e-12,
                             "Bloque traslacional 3×3 debe ser diagonal.")
        # Masa traslacional total = 2·ρAL.
        rho_a_l = 7850.0 * 1.5e-4 * length
        trasl_idx = [0, 1, 2, 6, 7, 8]
        self.assertAlmostEqual(np.sum(np.diag(M_l)[trasl_idx]) / 3.0,
                                 rho_a_l, places=10)


# ---------------------------------------------------------------------------
# Sólidos 2D — HRZ
# ---------------------------------------------------------------------------


def _build_tri3():
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [2.0, 0.0])
    n3 = dom.add_node(3, [0.0, 1.0])
    mat = Elastic2D(E=210e9, nu=0.3, density=7850.0)
    elem = Tri3(1, [n1, n2, n3], mat, thickness=0.05)
    dom.add_element(elem)
    return elem, 0.5 * 2.0 * 1.0 * 0.05  # área * espesor


def _build_quad4():
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [2.0, 0.0])
    n3 = dom.add_node(3, [2.0, 1.0])
    n4 = dom.add_node(4, [0.0, 1.0])
    mat = Elastic2D(E=210e9, nu=0.3, density=7850.0)
    elem = Quad4(1, [n1, n2, n3, n4], mat, thickness=0.05)
    dom.add_element(elem)
    return elem, 2.0 * 1.0 * 0.05


def _build_tri6():
    dom = Domain()
    # Triángulo recto vértices (0,0),(2,0),(0,1); puntos medios
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [2.0, 0.0])
    n3 = dom.add_node(3, [0.0, 1.0])
    n4 = dom.add_node(4, [1.0, 0.0])    # medio 1-2
    n5 = dom.add_node(5, [1.0, 0.5])    # medio 2-3
    n6 = dom.add_node(6, [0.0, 0.5])    # medio 3-1
    mat = Elastic2D(E=210e9, nu=0.3, density=7850.0)
    elem = Tri6(1, [n1, n2, n3, n4, n5, n6], mat, thickness=0.05)
    dom.add_element(elem)
    return elem, 0.5 * 2.0 * 1.0 * 0.05


def _build_quad8():
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [2.0, 0.0])
    n3 = dom.add_node(3, [2.0, 1.0])
    n4 = dom.add_node(4, [0.0, 1.0])
    n5 = dom.add_node(5, [1.0, 0.0])
    n6 = dom.add_node(6, [2.0, 0.5])
    n7 = dom.add_node(7, [1.0, 1.0])
    n8 = dom.add_node(8, [0.0, 0.5])
    mat = Elastic2D(E=210e9, nu=0.3, density=7850.0)
    elem = Quad8(1, [n1, n2, n3, n4, n5, n6, n7, n8], mat, thickness=0.05)
    dom.add_element(elem)
    return elem, 2.0 * 1.0 * 0.05


def _build_quad9():
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [2.0, 0.0])
    n3 = dom.add_node(3, [2.0, 1.0])
    n4 = dom.add_node(4, [0.0, 1.0])
    n5 = dom.add_node(5, [1.0, 0.0])
    n6 = dom.add_node(6, [2.0, 0.5])
    n7 = dom.add_node(7, [1.0, 1.0])
    n8 = dom.add_node(8, [0.0, 0.5])
    n9 = dom.add_node(9, [1.0, 0.5])
    mat = Elastic2D(E=210e9, nu=0.3, density=7850.0)
    elem = Quad9(1, [n1, n2, n3, n4, n5, n6, n7, n8, n9], mat,
                   thickness=0.05)
    dom.add_element(elem)
    return elem, 2.0 * 1.0 * 0.05


class TestLumpedSolid2D(unittest.TestCase):
    """Sólidos 2D: HRZ → diagonal, masa total preservada, entradas positivas."""

    def _check_solid(self, build_fn, rho=7850.0):
        elem, volume = build_fn()
        M_l = elem.compute_mass_matrix(lumping="lumped")
        m_total = rho * volume
        # Diagonalidad.
        off = M_l - np.diag(np.diag(M_l))
        self.assertLess(np.max(np.abs(off)), 1e-10,
                         f"{type(elem).__name__}: lumped no diagonal.")
        # Masa total por dirección.
        diag = np.diag(M_l)
        self.assertAlmostEqual(np.sum(diag) / 2.0, m_total,
                                  msg=f"{type(elem).__name__}: masa total no preservada")
        # Entradas positivas (M no singular).
        self.assertTrue(np.all(diag > 0.0),
                          f"{type(elem).__name__}: entradas no positivas {diag}")

    def test_tri3(self):
        self._check_solid(_build_tri3)

    def test_quad4(self):
        self._check_solid(_build_quad4)

    def test_tri6(self):
        self._check_solid(_build_tri6)

    def test_quad8(self):
        self._check_solid(_build_quad8)

    def test_quad9(self):
        self._check_solid(_build_quad9)


# ---------------------------------------------------------------------------
# Recuperación modal — barra axial con lumped
# ---------------------------------------------------------------------------


class TestModalLumpedRecovery(unittest.TestCase):
    """Modal con masa lumped recupera la frecuencia fundamental dentro de
    la tolerancia esperada del lumping en una malla razonable.
    """

    def test_axial_bar_first_frequency(self):
        """Barra axial empotrada-libre, N=20 elementos lineales.

        ω_exacto = (π/2L)·√(E/ρ). Con N elementos y lumping nodal:
        - Consistente: error ~0.1%.
        - Lumped (HRZ): error ~0.1–0.5% (mismo orden, signo opuesto).

        Tolerancia 2% es holgada y robusta a la elección de N.
        """
        dom = build_axial_bar_2d(Truss2D, n_elems=N_ELEMS)
        asm = Assembler(dom)
        omega_lumped = ModalSolver(asm, n_modes=1,
                                      lumping="lumped").solve().frequencies_rad[0]
        wave_speed = math.sqrt(E / RHO)
        omega_exact = (math.pi / (2.0 * L_TOTAL)) * wave_speed
        rel_err = abs(omega_lumped - omega_exact) / omega_exact
        self.assertLess(rel_err, 0.02,
                          f"lumped fundamental error {rel_err:.2%} > 2%")


if __name__ == "__main__":
    unittest.main()
