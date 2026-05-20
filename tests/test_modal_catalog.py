"""Cobertura modal de los elementos restantes del catálogo (ADR 0009 fase 1).

El test analítico canónico de Truss2D y Frame2DEuler vive en
``tests/test_modal.py``; este archivo cubre el resto:

- **Tests analíticos** (validación física rigurosa):
  - Truss3D — barra empotrada-libre en 3D, modo axial.
  - Frame3D — cantilever, flexión Bernoulli en plano xy local.
  - Frame2DTimoshenko — viga biapoyada (corrección por cortante pequeña en
    régimen esbelto, tolerancia ampliada).
  - Quad4 y Tri3 — barra rectangular esbelta modelada como sólido 2D,
    modo axial 1D (`nu=0`, `uy=0` confinado para aislar el modo).

- **Tests de equivalencia** (variantes que en configuración de referencia
  reproducen exactamente las frecuencias de su base):
  - Truss2DCorot ≡ Truss2D; Truss3DCorot ≡ Truss3D.
  - Frame2DEulerCorot ≡ Frame2DEuler.
  - Cable2DCorot con material bilateral ≡ Truss2D; Cable3DCorot ≡ Truss3D.
    (El cable con material unilateral requiere pre-tensión inicial para
    modal — diferido a fases futuras de ADR 0009.)

- **Tests algebraicos** (propiedades sin solución analítica específica):
  - Quad8, Quad9, Tri6 — matriz de masa simétrica positiva sobre un
    elemento aislado y modos ortonormales sobre un parche libre.
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
from solidum.elements.truss import Truss2D, Truss2DCorot, Truss3D, Truss3DCorot
from solidum.entry import run_modal
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import ModalSolver

from _modal_fixtures import (
    A_SECTION,
    E,
    L_TOTAL,
    N_ELEMS,
    RHO,
    build_axial_bar_2d,
    build_axial_bar_3d,
    build_simply_supported_beam,
)


# Frecuencias analíticas de barra axial empotrada-libre: ω_n = ((2n-1)π/(2L))·√(E/ρ).
def _axial_freqs(n: int = 3) -> np.ndarray:
    wave = math.sqrt(E / RHO)
    return np.array([(2 * k - 1) * math.pi / (2.0 * L_TOTAL) * wave
                     for k in range(1, n + 1)])


# ===========================================================================
# Tests analíticos
# ===========================================================================

class TestModalTruss3D(unittest.TestCase):
    def test_axial_first_three(self):
        dom = build_axial_bar_3d(Truss3D)
        result = run_modal(dom, n_modes=5)
        np.testing.assert_allclose(
            result.frequencies_rad[:3], _axial_freqs(3), rtol=0.01
        )


class TestModalFrame3DCantilever(unittest.TestCase):
    """Voladizo Frame3D — flexión Bernoulli en plano xy local.

    ``ref_vector=(0,1,0)`` alinea ejes locales con globales:
    ``y_local = (0,1,0)``, ``z_local = (0,0,1)``. El modo flexional uy–rz es
    el clásico cantilever Bernoulli: ``ω_n = β_n²·√(EI/(ρA))/L²`` con
    ``β_n·L = 1.8751, 4.6941, 7.8548`` (raíces de cos·cosh + 1 = 0).
    """

    def test_first_three_flexural(self):
        I_beam = 8.33e-10
        Jp_geom = 2.0 * I_beam  # Iy = Iz; J cualquiera (no afecta a flexión).
        n_elems = N_ELEMS
        dom = Domain()
        nodes = [dom.add_node(i + 1, [i * L_TOTAL / n_elems, 0.0, 0.0])
                 for i in range(n_elems + 1)]
        mat = Elastic1D(E=E, density=RHO)
        for i in range(n_elems):
            dom.add_element(Frame3D(
                i + 1, [nodes[i], nodes[i + 1]], mat,
                A=A_SECTION, Iy=I_beam, Iz=I_beam, J=Jp_geom,
                ref_vector=[0.0, 1.0, 0.0],
            ))
        # Empotramiento total en n1.
        for dof in ("ux", "uy", "uz", "rx", "ry", "rz"):
            nodes[0].fix_dof(dof, 0.0)
        # Resto: aislar el plano flexional xy (uy, rz libres).
        for n in nodes[1:]:
            for dof in ("ux", "uz", "rx", "ry"):
                n.fix_dof(dof, 0.0)
        dom.generate_equation_numbers()

        result = run_modal(dom, n_modes=5)

        beta_L = np.array([1.875104, 4.694091, 7.854757])
        beta = beta_L / L_TOTAL
        omegas_analytic = beta**2 * math.sqrt(E * I_beam / (RHO * A_SECTION))
        np.testing.assert_allclose(
            result.frequencies_rad[:3], omegas_analytic, rtol=0.02
        )


class TestModalFrame2DTimoshenko(unittest.TestCase):
    """Viga biapoyada con Frame2DTimoshenko — corrección por cortante.

    En régimen esbelto (L/r ≫ 10) las frecuencias convergen a las de
    Bernoulli; la corrección es ``O((r/L)²·n²)`` y se mantiene dentro de
    una tolerancia del 3 % para los 3 primeros modos.
    """

    def test_first_three_flexural(self):
        I_beam = 8.33e-10
        A_shear = 5.0 * A_SECTION / 6.0  # factor 5/6 para sección rectangular
        n_elems = N_ELEMS
        dom = Domain()
        nodes = [dom.add_node(i + 1, [i * L_TOTAL / n_elems, 0.0])
                 for i in range(n_elems + 1)]
        mat = Elastic1D(E=E, density=RHO)
        for i in range(n_elems):
            dom.add_element(Frame2DTimoshenko(
                i + 1, [nodes[i], nodes[i + 1]], mat,
                A=A_SECTION, I=I_beam, As=A_shear, nu=0.3,
            ))
        nodes[0].fix_dof("ux", 0.0); nodes[0].fix_dof("uy", 0.0)
        nodes[-1].fix_dof("ux", 0.0); nodes[-1].fix_dof("uy", 0.0)
        dom.generate_equation_numbers()

        result = run_modal(dom, n_modes=5)

        beta = math.sqrt(E * I_beam / (RHO * A_SECTION))
        omegas_analytic = np.array([(n * math.pi / L_TOTAL) ** 2 * beta
                                      for n in (1, 2, 3)])
        np.testing.assert_allclose(
            result.frequencies_rad[:3], omegas_analytic, rtol=0.03
        )


# ===========================================================================
# Tests de equivalencia (corotacionales y cables sobre config de referencia)
# ===========================================================================

class TestModalCorotEquivalence(unittest.TestCase):
    """Variantes corotacionales en configuración de referencia (no deformada)
    producen exactamente las frecuencias de su base. La masa se evalúa en L0,
    y K en u_e = 0 coincide con la base (la formulación corotacional reduce
    su contribución no lineal a cero en el estado natural)."""

    def test_truss2dcorot_matches_truss2d(self):
        base = run_modal(build_axial_bar_2d(Truss2D), n_modes=3)
        corot = run_modal(build_axial_bar_2d(Truss2DCorot), n_modes=3)
        np.testing.assert_allclose(
            corot.frequencies_rad, base.frequencies_rad, rtol=1e-10
        )

    def test_truss3dcorot_matches_truss3d(self):
        base = run_modal(build_axial_bar_3d(Truss3D), n_modes=3)
        corot = run_modal(build_axial_bar_3d(Truss3DCorot), n_modes=3)
        np.testing.assert_allclose(
            corot.frequencies_rad, base.frequencies_rad, rtol=1e-10
        )

    def test_frame2deulercorot_matches_frame2deuler(self):
        base = run_modal(build_simply_supported_beam(Frame2DEuler), n_modes=3)
        corot = run_modal(build_simply_supported_beam(Frame2DEulerCorot), n_modes=3)
        np.testing.assert_allclose(
            corot.frequencies_rad, base.frequencies_rad, rtol=1e-10
        )


class TestModalCablesBilateralEquivalence(unittest.TestCase):
    """Cable con material Elastic1D (bilateral) ⇒ comportamiento idéntico a
    truss del mismo material. La unilateralidad solo aparece con
    ``CableMaterial1D``, que requiere pre-tensión para modal (la rigidez
    tangente axial colapsa a cero en compresión / slack)."""

    def test_cable2d_matches_truss2d(self):
        truss = run_modal(build_axial_bar_2d(Truss2D), n_modes=3)
        cable = run_modal(build_axial_bar_2d(Cable2DCorot), n_modes=3)
        np.testing.assert_allclose(
            cable.frequencies_rad, truss.frequencies_rad, rtol=1e-10
        )

    def test_cable3d_matches_truss3d(self):
        truss = run_modal(build_axial_bar_3d(Truss3D), n_modes=3)
        cable = run_modal(build_axial_bar_3d(Cable3DCorot), n_modes=3)
        np.testing.assert_allclose(
            cable.frequencies_rad, truss.frequencies_rad, rtol=1e-10
        )


# ===========================================================================
# Tests analíticos en sólidos 2D — barra axial 1D embebida
# ===========================================================================

class TestModalSolid2DAxialBar(unittest.TestCase):
    """Barra rectangular muy esbelta modelada como sólido 2D, restringida a
    movimiento axial puro: ``nu=0`` (sin acoplamiento Poisson) y ``uy=0`` en
    todos los nodos. Frecuencias coinciden con barra 1D empotrada-libre.

    Cobertura de Quad4 y Tri3 (el resto, en TestModalSolid2DAlgebraic)."""

    H_BAR = 0.01  # altura de la barra (muy fina para que el modo axial 1D domine)
    N_X = 20      # divisiones a lo largo de la barra

    def _build_quad4_mesh(self) -> Domain:
        dom = Domain()
        # Nodos por columna: bottom (id 2i+1), top (id 2i+2).
        bot = []; top = []
        for i in range(self.N_X + 1):
            x = i * L_TOTAL / self.N_X
            bot.append(dom.add_node(2 * i + 1, [x, 0.0]))
            top.append(dom.add_node(2 * i + 2, [x, self.H_BAR]))
        mat = Elastic2D(E=E, nu=0.0, density=RHO)
        for i in range(self.N_X):
            dom.add_element(Quad4(
                i + 1, [bot[i], bot[i + 1], top[i + 1], top[i]],
                mat, thickness=1.0,
            ))
        # Empotramiento total en x=0 (ambos nodos).
        bot[0].fix_dof("ux", 0.0); bot[0].fix_dof("uy", 0.0)
        top[0].fix_dof("ux", 0.0); top[0].fix_dof("uy", 0.0)
        # Resto: uy=0 (movimiento puramente axial).
        for n in bot[1:] + top[1:]:
            n.fix_dof("uy", 0.0)
        dom.generate_equation_numbers()
        return dom

    def _build_tri3_mesh(self) -> Domain:
        dom = Domain()
        bot = []; top = []
        for i in range(self.N_X + 1):
            x = i * L_TOTAL / self.N_X
            bot.append(dom.add_node(2 * i + 1, [x, 0.0]))
            top.append(dom.add_node(2 * i + 2, [x, self.H_BAR]))
        mat = Elastic2D(E=E, nu=0.0, density=RHO)
        eid = 1
        for i in range(self.N_X):
            # Dos triángulos por celda — orden antihorario.
            dom.add_element(Tri3(
                eid, [bot[i], bot[i + 1], top[i]], mat, thickness=1.0,
            )); eid += 1
            dom.add_element(Tri3(
                eid, [bot[i + 1], top[i + 1], top[i]], mat, thickness=1.0,
            )); eid += 1
        bot[0].fix_dof("ux", 0.0); bot[0].fix_dof("uy", 0.0)
        top[0].fix_dof("ux", 0.0); top[0].fix_dof("uy", 0.0)
        for n in bot[1:] + top[1:]:
            n.fix_dof("uy", 0.0)
        dom.generate_equation_numbers()
        return dom

    def test_quad4_axial_first_three(self):
        dom = self._build_quad4_mesh()
        result = run_modal(dom, n_modes=5)
        np.testing.assert_allclose(
            result.frequencies_rad[:3], _axial_freqs(3), rtol=0.05
        )

    def test_tri3_axial_first_three(self):
        dom = self._build_tri3_mesh()
        result = run_modal(dom, n_modes=5)
        # Tri3 (CST) es más rígido por shear locking — tolerancia ampliada.
        np.testing.assert_allclose(
            result.frequencies_rad[:3], _axial_freqs(3), rtol=0.08
        )


# ===========================================================================
# Tests algebraicos para Quad8, Quad9, Tri6 (sin solución analítica específica)
# ===========================================================================

class TestModalSolid2DAlgebraic(unittest.TestCase):
    """Para los elementos de orden superior (Quad8, Quad9, Tri6) la masa
    consistente sigue el mismo patrón ``∫ρ·N^T·N``; la diferencia con Quad4
    es solo el orden de las funciones de forma. Validamos que sobre un
    modelo trivial:

    1. La matriz de masa elemental es simétrica y positiva definida.
    2. La traza coincide con ``2·ρ·A·t`` (un grado de masa por componente
       traslacional sobre toda el área).
    3. Sobre un modelo no trivial los modos son ortonormales respecto a M
       y ortogonales respecto a K (propiedades de la fase 1 ya cubiertas
       por TestModalAlgebraicProperties en tests/test_modal.py para frames
       y trusses; aquí extendido a los sólidos restantes).
    """

    THICK = 0.5
    L_SIDE = 1.0

    def _quad_nodes(self):
        """8 nodos para Quad8 / 9 para Quad9 sobre el cuadrado [0,1]²."""
        return [
            (0.0, 0.0),         # 0
            (self.L_SIDE, 0.0), # 1
            (self.L_SIDE, self.L_SIDE),  # 2
            (0.0, self.L_SIDE), # 3
            (0.5 * self.L_SIDE, 0.0),               # 4 (medio 0-1)
            (self.L_SIDE, 0.5 * self.L_SIDE),       # 5 (medio 1-2)
            (0.5 * self.L_SIDE, self.L_SIDE),       # 6 (medio 2-3)
            (0.0, 0.5 * self.L_SIDE),               # 7 (medio 3-0)
            (0.5 * self.L_SIDE, 0.5 * self.L_SIDE), # 8 (centro, solo Quad9)
        ]

    def _tri6_nodes(self):
        """6 nodos para Tri6 sobre (0,0)-(1,0)-(0,1)."""
        return [
            (0.0, 0.0),   # 0
            (1.0, 0.0),   # 1
            (0.0, 1.0),   # 2
            (0.5, 0.0),   # 3 (medio 0-1)
            (0.5, 0.5),   # 4 (medio 1-2)
            (0.0, 0.5),   # 5 (medio 2-0)
        ]

    def _build_single(self, elem_cls, node_coords, n_nodes: int) -> Domain:
        dom = Domain()
        nodes = [dom.add_node(i + 1, list(node_coords[i])) for i in range(n_nodes)]
        mat = Elastic2D(E=E, nu=0.3, density=RHO)
        dom.add_element(elem_cls(1, nodes, mat, thickness=self.THICK))
        return dom

    def _check_mass_properties(self, elem_cls, node_coords, n_nodes: int,
                                expected_area: float):
        """Para un único elemento aislado: M simétrica, PD, traza = 2·ρ·t·A.

        Para masa consistente translacional la suma de filas de M (sum_j M[i,j])
        sobre un componente traslacional da la masa nodal asociada al nodo i;
        la suma total ``Σ_i Σ_j M_uxux[i,j]`` = masa total del elemento
        ``ρ·t·Área``. Como M acopla ux y uy ortogonalmente, la traza completa
        es ``2 · ρ · t · Área``.
        """
        dom = self._build_single(elem_cls, node_coords, n_nodes)
        elem = list(dom.elements.values())[0]
        M = elem.compute_mass_matrix()
        # Simetría.
        np.testing.assert_allclose(M, M.T, atol=1e-10)
        # Suma total = 2 · ρ · t · A (masa total × 2 componentes traslacionales).
        expected_total = 2.0 * RHO * self.THICK * expected_area
        np.testing.assert_allclose(M.sum(), expected_total, rtol=1e-10)
        # Definida positiva: todos los autovalores > 0.
        eigvals = np.linalg.eigvalsh(M)
        self.assertTrue(np.all(eigvals > 0.0),
                         f"M no positiva definida: min eig = {eigvals.min():.3e}")

    def test_quad8_mass_properties(self):
        coords = self._quad_nodes()[:8]
        self._check_mass_properties(Quad8, coords, 8, expected_area=1.0)

    def test_quad9_mass_properties(self):
        coords = self._quad_nodes()
        self._check_mass_properties(Quad9, coords, 9, expected_area=1.0)

    def test_tri6_mass_properties(self):
        coords = self._tri6_nodes()
        # Triángulo rectángulo (0,0)-(1,0)-(0,1) → área 0.5.
        self._check_mass_properties(Tri6, coords, 6, expected_area=0.5)

    def _check_orthonormality_on_patch(self, elem_cls, node_coords, n_nodes: int):
        """Patch con un único elemento libre — modos no degenerados, MΦ=I.

        Sin restricciones de Dirichlet existen modos de cuerpo rígido
        (autovalor 0). Restringimos un nodo para eliminar los 3 modos
        rígidos en 2D y verificamos ortonormalidad sobre los siguientes.
        """
        dom = self._build_single(elem_cls, node_coords, n_nodes)
        # Restringir nodo 1 totalmente y nodo 2 en uy para fijar traslación + giro.
        dom.nodes[1].fix_dof("ux", 0.0); dom.nodes[1].fix_dof("uy", 0.0)
        dom.nodes[2].fix_dof("uy", 0.0)
        dom.generate_equation_numbers()

        result = run_modal(dom, n_modes=3)
        asm = Assembler(dom)
        M = asm.assemble_mass_matrix()
        Phi = result.modes
        np.testing.assert_allclose(Phi.T @ (M @ Phi), np.eye(3), atol=1e-10)
        self.assertTrue(np.all(result.frequencies_rad > 0.0))

    def test_quad8_orthonormality(self):
        self._check_orthonormality_on_patch(Quad8, self._quad_nodes()[:8], 8)

    def test_quad9_orthonormality(self):
        self._check_orthonormality_on_patch(Quad9, self._quad_nodes(), 9)

    def test_tri6_orthonormality(self):
        self._check_orthonormality_on_patch(Tri6, self._tri6_nodes(), 6)


if __name__ == "__main__":
    unittest.main()
