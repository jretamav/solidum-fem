"""Participación modal en el resultado del análisis modal y movimientos de
cuerpo rígido (``solidum.math.rigid_body``).

Criterios de aceptación:

1. Edificio de cortante de tres pisos (masas de piso, columnas sin masa): los
   factores de participación ``Γₙ = φₙᵀ·m·1`` y las masas efectivas ``Γₙ²``
   son los del problema 3×3 resuelto a mano (diferencia relativa ≤ 1e-6); la
   suma de las tres masas efectivas es la masa de los pisos; ``condensed_dofs``
   nombra exactamente los grados de libertad libres sin masa.
2. Completitud: con todos los modos, ``Σ Γₙ² = bᵀ·M_LL⁻¹·b`` con
   ``b = (M·r)_L`` sobre los grados de libertad libres ``L``; con masa
   concentrada, la masa de los nudos libres. Con masa consistente y
   concentrada, en ``ux`` y ``uy``.
3. Las masas efectivas del modal son las del espectral con los mismos modos.
4. Un movimiento de cuerpo rígido no deforma: ``K·r = 0`` sin apoyos, para
   las traslaciones y los giros alrededor de un punto cualquiera, en marcos 2D
   y 3D con barras inclinadas, armaduras 3D y sólidos 2D y 3D.
5. Un giro alrededor de un eje a distancia ``d`` de una masa ``m`` con
   inercia ``J`` da ``rᵀ·M·r = m·d² + J`` (teorema de los ejes paralelos).
6. Los errores de uso se informan con un mensaje claro.
"""
import math
import unittest

import numpy as np
import scipy.linalg as sla

from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler
from solidum.elements.frame3d import Frame3D
from solidum.elements.solid_2d import Quad4
from solidum.elements.solid_3d.hex8 import Hex8
from solidum.elements.truss import Truss3D
from solidum.materials.elastic import Elastic1D
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.elastic_3d import Elastic3D
from solidum.math.modal_response import participation_factors
from solidum.math.rigid_body import dof_unit_vector, rigid_body_vector
from solidum.math.solvers import ModalSolver, ResponseSpectrumSolver
from test_nodal_mass import shear_building, shear_building_3x3

E, RHO, A, I = 2.0e11, 7850.0, 6.0e-3, 8.5e-5
RTOL_ACCEPT = 1.0e-6


def portal(density=RHO):
    """Marco de dos columnas y una viga, con una diagonal inclinada."""
    dom = Domain()
    mat = Elastic1D(E=E, density=density)
    pts = {1: (0.0, 0.0), 2: (0.0, 3.2), 3: (2.5, 3.2), 4: (5.0, 3.2), 5: (5.0, 0.0)}
    nodes = {k: dom.add_node(k, list(v)) for k, v in pts.items()}
    for eid, (a, b) in enumerate([(1, 2), (2, 3), (3, 4), (4, 5), (1, 3)], start=1):
        dom.add_element(Frame2DEuler(eid, [nodes[a], nodes[b]], material=mat, A=A, I=I))
    for n in (nodes[1], nodes[5]):
        for dof in ("ux", "uy", "rz"):
            n.fix_dof(dof, 0.0)
    dom.generate_equation_numbers()
    return dom


class TestEdificioDeCortante(unittest.TestCase):
    """Criterio 1."""

    def test_participacion_contra_el_problema_3x3(self):
        heights, masses = [4.0, 3.2, 3.0], [5.2e4, 4.6e4, 3.1e4]
        inertias = [4.2e-4, 3.0e-4, 2.2e-4]
        K3, M3 = shear_building_3x3(heights, masses, inertias)
        lam, vec = sla.eigh(K3, M3)
        gamma3 = vec.T @ M3 @ np.ones(3)
        dom, floors = shear_building(heights, masses, inertias)
        res = ModalSolver(Assembler(dom), n_modes=3).solve()

        self.assertEqual(set(res.participation_factors), {"ux", "uy"})
        self.assertEqual(set(res.effective_masses), {"ux", "uy"})
        roof = floors[-1].dofs["ux"]
        for n in range(3):
            # El signo de Γ sigue al del modo: Γₙ·φₙ no depende de él.
            self.assertLessEqual(
                abs(res.participation_factors["ux"][n] * res.modes[roof, n]
                    / (gamma3[n] * vec[2, n]) - 1.0), RTOL_ACCEPT)
        np.testing.assert_allclose(res.effective_masses["ux"], gamma3 ** 2, rtol=RTOL_ACCEPT)
        self.assertAlmostEqual(np.sum(res.effective_masses["ux"]) / sum(masses), 1.0, delta=1e-10)
        # Sin masa en uy: los pisos no se mueven en vertical.
        np.testing.assert_allclose(res.effective_masses["uy"], 0.0, atol=1e-20)

        # Libres sin masa: todos los de los nudos intermedios de las columnas.
        # Los nudos de piso sólo tienen ux libre, con masa (el derecho es
        # esclavo del izquierdo).
        floor_y = {f.coordinates[1] for f in floors}
        expected = {(n.id, d) for n in dom.nodes.values() for d in n.dofs
                    if d not in n.boundary_conditions and n.coordinates[1] not in floor_y}
        self.assertEqual(set(res.condensed_dofs), expected)
        self.assertEqual(len(res.condensed_dofs), len(expected))


class TestCompletitud(unittest.TestCase):
    """Criterios 2 y 3."""

    def test_suma_de_todas_las_masas_efectivas(self):
        for lumping in ("consistent", "lumped"):
            dom = portal()
            asm = Assembler(dom)
            M = asm.assemble_mass_matrix(lumping=lumping).toarray()
            free = np.array([i for i in range(dom.total_dofs)
                             if not any(i == n.dofs[d] for n in dom.nodes.values()
                                        for d in n.boundary_conditions)])
            M_LL = M[np.ix_(free, free)]
            n_mass = int(np.count_nonzero(np.diag(M_LL)))
            res = ModalSolver(asm, n_modes=n_mass, lumping=lumping).solve()
            for name in ("ux", "uy"):
                with self.subTest(lumping=lumping, direction=name):
                    r = dof_unit_vector(dom, name)
                    b = (M @ r)[free]
                    if lumping == "lumped":
                        expected = float(r[free] @ M_LL @ r[free])
                    else:
                        expected = float(b @ np.linalg.solve(M_LL, b))
                    total = float(np.sum(res.effective_masses[name]))
                    self.assertAlmostEqual(total / expected, 1.0, delta=1e-9)
                    np.testing.assert_allclose(
                        res.participation_factors[name],
                        participation_factors(res.modes, M, r), rtol=1e-12, atol=0.0)
            # La masa concentrada de los marcos lleva inercia de giro: todos
            # los grados de libertad libres tienen masa.
            self.assertEqual(res.condensed_dofs, ())

    def test_igual_que_el_espectral(self):
        dom = portal()
        modal = ModalSolver(Assembler(dom), n_modes=4).solve()
        for name in ("ux", "uy"):
            spec = ResponseSpectrumSolver(
                Assembler(dom), n_modes=4, direction={"dof_name": name},
                spectrum={"type": "constant_acceleration", "Sa": 1.0}).solve()
            np.testing.assert_allclose(modal.effective_masses[name], spec.effective_masses,
                                       rtol=1e-10, atol=1e-12 * spec.total_mass)


def _free_frame3d():
    dom = Domain()
    mat = Elastic1D(E=E, density=RHO)
    pts = [(0.0, 0.0, 0.0), (3.0, 0.5, 1.0), (4.0, 3.0, 2.5), (1.0, 2.0, 4.0)]
    nodes = [dom.add_node(i + 1, list(p)) for i, p in enumerate(pts)]
    for k, (a, b) in enumerate([(0, 1), (1, 2), (2, 3), (3, 0)], start=1):
        dom.add_element(Frame3D(k, [nodes[a], nodes[b]], mat, A=A, Iy=I, Iz=0.7 * I,
                                J=1.3 * I, nu=0.3))
    return dom


def _free_frame2d():
    dom = Domain()
    mat = Elastic1D(E=E, density=RHO)
    pts = [(0.0, 0.0), (3.0, 1.2), (4.5, 4.0), (0.5, 3.0)]
    nodes = [dom.add_node(i + 1, list(p)) for i, p in enumerate(pts)]
    for k, (a, b) in enumerate([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)], start=1):
        dom.add_element(Frame2DEuler(k, [nodes[a], nodes[b]], material=mat, A=A, I=I))
    return dom


def _free_truss3d():
    dom = Domain()
    mat = Elastic1D(E=E, density=RHO)
    pts = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.3), (0.4, 1.8, 0.0), (0.7, 0.6, 1.9)]
    nodes = [dom.add_node(i + 1, list(p)) for i, p in enumerate(pts)]
    for k, (a, b) in enumerate([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)], start=1):
        dom.add_element(Truss3D(k, [nodes[a], nodes[b]], mat, A=A))
    return dom


def _free_quad4():
    dom = Domain()
    mat = Elastic2D(E=E, nu=0.3, density=RHO)
    pts = [(0.0, 0.0), (2.0, 0.2), (2.3, 1.7), (0.1, 1.4)]
    nodes = [dom.add_node(i + 1, list(p)) for i, p in enumerate(pts)]
    dom.add_element(Quad4(1, nodes, mat, thickness=0.2))
    return dom


def _free_hex8():
    dom = Domain()
    mat = Elastic3D(E=E, nu=0.3, density=RHO)
    coords = [(0.0, 0.0, 0.0), (1.2, 0.1, 0.0), (1.3, 1.1, 0.1), (0.0, 1.0, 0.0),
              (0.1, 0.0, 0.9), (1.1, 0.0, 1.0), (1.2, 1.2, 1.1), (0.0, 1.0, 1.0)]
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
    dom.add_element(Hex8(1, nodes, mat))
    return dom


class TestMovimientoDeCuerpoRigido(unittest.TestCase):
    """Criterios 4 y 5."""

    def test_no_deforma(self):
        cases = {
            "Frame2DEuler": (_free_frame2d, ("ux", "uy", "rz"), (0.7, -1.3)),
            "Quad4": (_free_quad4, ("ux", "uy", "rz"), (-2.0, 0.4)),
            "Frame3D": (_free_frame3d, ("ux", "uy", "uz", "rx", "ry", "rz"), (0.3, -1.1, 2.2)),
            "Truss3D": (_free_truss3d, ("ux", "uy", "uz", "rx", "ry", "rz"), (1.5, 0.2, -0.8)),
            "Hex8": (_free_hex8, ("ux", "uy", "uz", "rx", "ry", "rz"), (-0.5, 2.0, 0.6)),
        }
        for label, (build, motions, about) in cases.items():
            dom = build()
            dom.generate_equation_numbers()
            asm = Assembler(dom)
            asm.assemble_system()
            K = asm.K_global.toarray()
            scale = np.max(np.abs(K))
            for motion in motions:
                with self.subTest(element=label, motion=motion):
                    kw = {"about": about} if motion.startswith("r") else {}
                    r = rigid_body_vector(dom, motion, **kw)
                    self.assertGreater(np.max(np.abs(r)), 0.0)
                    self.assertLessEqual(np.max(np.abs(K @ r)),
                                         1e-11 * scale * np.max(np.abs(r)))

    def test_ejes_paralelos(self):
        m, J, x, y = 3.0e3, 450.0, 2.5, 1.0
        dom = Domain()
        mat = Elastic1D(E=E, density=0.0)
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [x, y])
        dom.add_element(Frame2DEuler(1, [n1, n2], material=mat, A=A, I=I))
        dom.add_nodal_mass(2, ux=m, uy=m, rz=J)
        dom.generate_equation_numbers()
        M = Assembler(dom).assemble_mass_matrix().toarray()
        c = (-1.0, 4.0)
        r = rigid_body_vector(dom, "rz", about=c)
        d2 = (x - c[0]) ** 2 + (y - c[1]) ** 2
        self.assertAlmostEqual(float(r @ M @ r) / (m * d2 + J), 1.0, delta=1e-14)
        # El giro mueve el punto perpendicularmente al radio: e_z × (x − c).
        self.assertAlmostEqual(r[n2.dofs["ux"]], -(y - c[1]), delta=1e-15)
        self.assertAlmostEqual(r[n2.dofs["uy"]], x - c[0], delta=1e-15)
        self.assertEqual(r[n2.dofs["rz"]], 1.0)


class TestErrores(unittest.TestCase):
    """Criterio 6."""

    def setUp(self):
        self.dom = _free_quad4()
        self.dom.generate_equation_numbers()

    def test_mensajes(self):
        cases = [
            (dict(motion="rz"), "necesita el punto"),
            (dict(motion="ux", about=(0.0, 0.0)), "no tiene punto de giro"),
            (dict(motion="rx", about=(0.0, 0.0)), "sólo existe el giro 'rz'"),
            (dict(motion="rz", about=(0.0, 0.0, 0.0)), "debe tener 2 coordenadas"),
            (dict(motion="theta", about=None), "desconocido"),
            (dict(motion="uz"), "ningún nodo tiene el grado de libertad 'uz'"),
        ]
        for kw, msg in cases:
            with self.subTest(**{k: str(v) for k, v in kw.items()}):
                motion = kw.pop("motion")
                with self.assertRaisesRegex(ValueError, msg):
                    rigid_body_vector(self.dom, motion, **kw)


if __name__ == "__main__":
    unittest.main()
