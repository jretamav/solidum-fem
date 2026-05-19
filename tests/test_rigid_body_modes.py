"""Rigid body modes explícitos para todos los elementos del catálogo.

Verifica las propiedades canónicas de cualquier elemento finito bien
formado:

1. **Traslación rígida** (`u(x) = (a, b)` constante en todos los nodos)
   produce deformación nula en cada punto de integración y fuerza
   interna nula en cada DOF.
2. **Rotación rígida infinitesimal** (`u(x) = (-θ·y, θ·x)`) produce
   deformación nula (la parte antisimétrica del gradiente no contribuye
   a `ε = sym(∇u)`) y fuerza interna nula.
3. **Rank-deficiency de `K`**: la matriz de rigidez elemental tiene
   exactamente 3 modos cero en 2D (2 traslaciones + 1 rotación), o
   6 en 3D.

Para elementos corotacionales (`Truss2DCorot`, `Frame2DEulerCorot`), se
verifica adicionalmente que una **rotación rígida GRANDE** (90°)
produce deformación axial nula — propiedad definitoria de la
formulación corotacional.

Estos comportamientos están cubiertos hoy implícitamente por otros
tests (patch tests, validación analítica de cantilevers, etc.), pero
un test directo es más diagnóstico y se ejecuta en aislamiento sobre
un único elemento sin acoplamiento del solver.
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.domain import Domain
from fenix.core.material import Material
from fenix.elements.frame.euler_corot import Frame2DEulerCorot
from fenix.elements.solid_2d import Quad4, Quad8, Quad9, Tri3, Tri6
from fenix.elements.solid_3d import Hex8, Tet4
from fenix.elements.truss import Truss2D, Truss2DCorot
from fenix.materials.elastic import Elastic1D
from fenix.materials.elastic_3d import Elastic3D


class _PlaneStress(Material):
    STRAIN_DIM = 3

    def __init__(self, E, nu):
        self.E, self.nu = E, nu
        fac = E / (1.0 - nu * nu)
        self.C = fac * np.array([
            [1.0, nu, 0.0],
            [nu, 1.0, 0.0],
            [0.0, 0.0, (1.0 - nu) / 2.0],
        ])

    def compute_state(self, strain, state_vars=None):
        return self.C @ strain, self.C, state_vars


def _rigid_translation_u_e(elem, ax: float, ay: float) -> np.ndarray:
    """u_e de traslación rígida (ax, ay) leyendo `elem.DOF_NAMES`."""
    dof_names = elem.DOF_NAMES
    n_dofs = len(dof_names)
    n_nodes = len(elem.nodes)
    u_e = np.zeros(n_nodes * n_dofs)
    for i in range(n_nodes):
        for j, name in enumerate(dof_names):
            if name == 'ux':
                u_e[i * n_dofs + j] = ax
            elif name == 'uy':
                u_e[i * n_dofs + j] = ay
            # rz / uz / etc: 0 en una traslación pura.
    return u_e


def _rigid_rotation_u_e(elem, theta: float, center=(0.0, 0.0)) -> np.ndarray:
    """u_e de rotación rígida infinitesimal alrededor de `center`.

    Para nodos con `rz`, asigna además `θ` al DOF rotacional para que la
    rotación sea consistente con la cinemática de viga.
    """
    dof_names = elem.DOF_NAMES
    n_dofs = len(dof_names)
    n_nodes = len(elem.nodes)
    u_e = np.zeros(n_nodes * n_dofs)
    cx, cy = center
    for i, node in enumerate(elem.nodes):
        x, y = node.coordinates[0], node.coordinates[1]
        for j, name in enumerate(dof_names):
            if name == 'ux':
                u_e[i * n_dofs + j] = -theta * (y - cy)
            elif name == 'uy':
                u_e[i * n_dofs + j] = theta * (x - cx)
            elif name == 'rz':
                u_e[i * n_dofs + j] = theta
    return u_e


def _build_quad4(material):
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [1.0, 0.0])
    n3 = dom.add_node(3, [1.0, 1.0])
    n4 = dom.add_node(4, [0.0, 1.0])
    elem = Quad4(1, [n1, n2, n3, n4], material, thickness=1.0)
    dom.add_element(elem)
    return elem


def _build_quad8(material):
    dom = Domain()
    nodes = []
    coords = [
        (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0),   # corners
        (0.5, 0.0), (1.0, 0.5), (0.5, 1.0), (0.0, 0.5),   # midnodes
    ]
    for i, xy in enumerate(coords, start=1):
        nodes.append(dom.add_node(i, list(xy)))
    elem = Quad8(1, nodes, material, thickness=1.0)
    dom.add_element(elem)
    return elem


def _build_quad9(material):
    dom = Domain()
    nodes = []
    coords = [
        (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0),   # corners
        (0.5, 0.0), (1.0, 0.5), (0.5, 1.0), (0.0, 0.5),   # midnodes
        (0.5, 0.5),                                        # center
    ]
    for i, xy in enumerate(coords, start=1):
        nodes.append(dom.add_node(i, list(xy)))
    elem = Quad9(1, nodes, material, thickness=1.0)
    dom.add_element(elem)
    return elem


def _build_tri3(material):
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [1.0, 0.0])
    n3 = dom.add_node(3, [0.3, 1.0])
    elem = Tri3(1, [n1, n2, n3], material, thickness=1.0)
    dom.add_element(elem)
    return elem


def _build_tri6(material):
    dom = Domain()
    # 3 vértices distorsionados + 3 midnodes (en midpoint geométrico)
    v1, v2, v3 = (0.0, 0.0), (1.0, 0.0), (0.3, 1.0)
    m12 = (0.5 * (v1[0] + v2[0]), 0.5 * (v1[1] + v2[1]))
    m23 = (0.5 * (v2[0] + v3[0]), 0.5 * (v2[1] + v3[1]))
    m31 = (0.5 * (v3[0] + v1[0]), 0.5 * (v3[1] + v1[1]))
    coords = [v1, v2, v3, m12, m23, m31]
    nodes = []
    for i, xy in enumerate(coords, start=1):
        nodes.append(dom.add_node(i, list(xy)))
    elem = Tri6(1, nodes, material, thickness=1.0)
    dom.add_element(elem)
    return elem


# =============================================================================
# Sólidos 2D — translation + rotation + rank-deficiency
# =============================================================================

class _RigidBodyModesSolid2DMixin:
    """Verificación común para Q4/Q8/Q9/Tri3/Tri6 (small strain isótropos)."""

    BUILDER = None

    def setUp(self):
        self.material = _PlaneStress(E=2.0e5, nu=0.3)
        self.elem = self.BUILDER(self.material)

    def test_rigid_translation_produces_zero_internal_force(self):
        u_e = _rigid_translation_u_e(self.elem, ax=0.5, ay=-0.3)
        _, F_int = self.elem.compute_element_state(u_e)
        np.testing.assert_allclose(F_int, 0.0, atol=1.0e-9,
            err_msg=f"{type(self.elem).__name__}: traslación rígida produce "
                    f"fuerzas internas no nulas: max|F_int|={np.max(np.abs(F_int)):.3e}")

    def test_rigid_rotation_produces_zero_internal_force(self):
        u_e = _rigid_rotation_u_e(self.elem, theta=1.0e-4, center=(0.5, 0.5))
        _, F_int = self.elem.compute_element_state(u_e)
        # Para rotación infinitesimal, ε = sym(∇u) = 0 exactamente; F_int
        # debería ser nulo a precisión máquina (no hay términos cuadráticos
        # en small strain). Tolerancia laxa para absorber redondeo en B.
        np.testing.assert_allclose(F_int, 0.0, atol=1.0e-6,
            err_msg=f"{type(self.elem).__name__}: rotación rígida infinitesimal "
                    f"produce fuerzas internas no nulas: "
                    f"max|F_int|={np.max(np.abs(F_int)):.3e}")

    def test_stiffness_has_three_zero_modes(self):
        """K elemental tiene exactamente 3 autovalores ≈ 0 en 2D (2 trans + 1 rot)."""
        K, _ = self.elem.compute_element_state(np.zeros(self._dof_count()))
        eigvals = np.sort(np.abs(np.linalg.eigvalsh(0.5 * (K + K.T))))
        # Escala de comparación: máximo autovalor de K.
        scale = float(eigvals[-1])
        threshold = 1.0e-9 * scale
        n_zeros = int(np.sum(eigvals < threshold))
        self.assertEqual(n_zeros, 3,
            f"{type(self.elem).__name__}: K tiene {n_zeros} modos ≈ 0, "
            f"esperaba 3. eigvals[:5]={eigvals[:5]}, threshold={threshold:.3e}")

    def _dof_count(self) -> int:
        return len(self.elem.nodes) * len(self.elem.DOF_NAMES)


class TestRBMQuad4(_RigidBodyModesSolid2DMixin, unittest.TestCase):
    BUILDER = staticmethod(_build_quad4)


class TestRBMQuad8(_RigidBodyModesSolid2DMixin, unittest.TestCase):
    BUILDER = staticmethod(_build_quad8)


class TestRBMQuad9(_RigidBodyModesSolid2DMixin, unittest.TestCase):
    BUILDER = staticmethod(_build_quad9)


class TestRBMTri3(_RigidBodyModesSolid2DMixin, unittest.TestCase):
    BUILDER = staticmethod(_build_tri3)


class TestRBMTri6(_RigidBodyModesSolid2DMixin, unittest.TestCase):
    BUILDER = staticmethod(_build_tri6)


# =============================================================================
# Sólidos 3D — translation + rotation (x3 ejes) + rank-deficiency (6 modos)
# =============================================================================

def _build_hex8(material):
    dom = Domain()
    coords = [
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    ]
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
    elem = Hex8(1, nodes, material)
    dom.add_element(elem)
    return elem


def _build_tet4(material):
    dom = Domain()
    coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
    elem = Tet4(1, nodes, material)
    dom.add_element(elem)
    return elem


def _rigid_translation_u_e_3d(elem, ax, ay, az):
    n_nodes = len(elem.nodes)
    u_e = np.zeros(n_nodes * 3)
    for i in range(n_nodes):
        u_e[3 * i]     = ax
        u_e[3 * i + 1] = ay
        u_e[3 * i + 2] = az
    return u_e


def _rigid_rotation_u_e_3d(elem, omega, center=(0.0, 0.0, 0.0)):
    """Rotación rígida infinitesimal: u = ω × (x - c) con ω vector 3D pequeño."""
    omega = np.asarray(omega, dtype=np.float64)
    center = np.asarray(center, dtype=np.float64)
    n_nodes = len(elem.nodes)
    u_e = np.zeros(n_nodes * 3)
    for i, node in enumerate(elem.nodes):
        x = np.asarray(node.coordinates[:3], dtype=np.float64)
        u_e[3 * i:3 * i + 3] = np.cross(omega, x - center)
    return u_e


class _RigidBodyModesSolid3DMixin:
    """Verificación común para Hex8 y Tet4."""

    BUILDER = None

    def setUp(self):
        self.material = Elastic3D(E=2.0e5, nu=0.3)
        self.elem = self.BUILDER(self.material)

    def test_rigid_translation_produces_zero_internal_force(self):
        u_e = _rigid_translation_u_e_3d(self.elem, ax=0.5, ay=-0.3, az=0.7)
        _, F_int = self.elem.compute_element_state(u_e)
        np.testing.assert_allclose(
            F_int, 0.0, atol=1.0e-9,
            err_msg=f"{type(self.elem).__name__}: traslación rígida produce "
                    f"max|F_int|={np.max(np.abs(F_int)):.3e}",
        )

    def test_rigid_rotation_x_produces_zero_internal_force(self):
        u_e = _rigid_rotation_u_e_3d(
            self.elem, omega=(1.0e-4, 0.0, 0.0), center=(0.5, 0.5, 0.5)
        )
        _, F_int = self.elem.compute_element_state(u_e)
        np.testing.assert_allclose(F_int, 0.0, atol=1.0e-6)

    def test_rigid_rotation_y_produces_zero_internal_force(self):
        u_e = _rigid_rotation_u_e_3d(
            self.elem, omega=(0.0, 1.0e-4, 0.0), center=(0.5, 0.5, 0.5)
        )
        _, F_int = self.elem.compute_element_state(u_e)
        np.testing.assert_allclose(F_int, 0.0, atol=1.0e-6)

    def test_rigid_rotation_z_produces_zero_internal_force(self):
        u_e = _rigid_rotation_u_e_3d(
            self.elem, omega=(0.0, 0.0, 1.0e-4), center=(0.5, 0.5, 0.5)
        )
        _, F_int = self.elem.compute_element_state(u_e)
        np.testing.assert_allclose(F_int, 0.0, atol=1.0e-6)

    def test_stiffness_has_six_zero_modes(self):
        """K elemental tiene exactamente 6 autovalores ≈ 0 en 3D (3 trans + 3 rot)."""
        n_dofs = len(self.elem.nodes) * 3
        K, _ = self.elem.compute_element_state(np.zeros(n_dofs))
        eigvals = np.sort(np.abs(np.linalg.eigvalsh(0.5 * (K + K.T))))
        scale = float(eigvals[-1])
        threshold = 1.0e-9 * scale
        n_zeros = int(np.sum(eigvals < threshold))
        self.assertEqual(
            n_zeros, 6,
            f"{type(self.elem).__name__}: K tiene {n_zeros} modos ≈ 0, "
            f"esperaba 6. eigvals[:8]={eigvals[:8]}, threshold={threshold:.3e}",
        )


class TestRBMHex8(_RigidBodyModesSolid3DMixin, unittest.TestCase):
    BUILDER = staticmethod(_build_hex8)


class TestRBMTet4(_RigidBodyModesSolid3DMixin, unittest.TestCase):
    BUILDER = staticmethod(_build_tet4)


# =============================================================================
# Corotacionales — translation + large rotation
# =============================================================================

class TestRBMTruss2DCorot(unittest.TestCase):
    """`Truss2DCorot`: además del small-rotation, verifica rotación rígida grande."""

    def setUp(self):
        self.material = Elastic1D(E=2.0e5)
        dom = Domain()
        # Barra horizontal de longitud L=1 con extremo en el origen.
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [1.0, 0.0])
        self.elem = Truss2DCorot(1, [n1, n2], self.material, A=1.0)
        dom.add_element(self.elem)

    def test_rigid_translation_produces_zero_axial_force(self):
        u_e = _rigid_translation_u_e(self.elem, ax=0.7, ay=-0.4)
        _, F_int = self.elem.compute_element_state(u_e)
        np.testing.assert_allclose(F_int, 0.0, atol=1.0e-9,
            err_msg=f"Traslación rígida produce N no nulo: "
                    f"max|F_int|={np.max(np.abs(F_int)):.3e}")

    def test_large_rigid_rotation_produces_zero_axial_force(self):
        """Rotación rígida de 90° alrededor del nodo 1.

        Configuración corriente: nodo 1 en (0,0), nodo 2 en (0,1). La
        longitud actual sigue siendo L=1 (no estiramiento), así que el
        corotacional debe captar ε = 0 y N = 0 sin importar la magnitud
        de la rotación.
        """
        theta = math.pi / 2.0  # 90°
        # Nodo 1 fijo en (0,0). Nodo 2 inicial (1,0) → final (cos θ, sin θ) = (0,1).
        u_e = np.array([
            0.0, 0.0,                                   # nodo 1
            math.cos(theta) - 1.0, math.sin(theta) - 0.0,  # nodo 2: desde (1,0) a (0,1)
        ])
        _, F_int = self.elem.compute_element_state(u_e)
        np.testing.assert_allclose(F_int, 0.0, atol=1.0e-9,
            err_msg=f"Rotación rígida 90° produce N no nulo: "
                    f"max|F_int|={np.max(np.abs(F_int)):.3e}; corot no "
                    "está extrayendo la rotación rígida correctamente.")


class TestRBMFrame2DEulerCorot(unittest.TestCase):
    """`Frame2DEulerCorot`: traslación + rotación rígida grande con DOFs `rz`."""

    def setUp(self):
        self.material = Elastic1D(E=2.0e5)
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [1.0, 0.0])
        self.elem = Frame2DEulerCorot(
            1, [n1, n2], self.material, A=1.0, I=1.0e-4,
        )
        dom.add_element(self.elem)

    def test_rigid_translation_produces_zero_forces(self):
        u_e = _rigid_translation_u_e(self.elem, ax=0.3, ay=0.5)
        _, F_int = self.elem.compute_element_state(u_e)
        np.testing.assert_allclose(F_int, 0.0, atol=1.0e-9,
            err_msg=f"Frame2DEulerCorot traslación: max|F_int|={np.max(np.abs(F_int)):.3e}")

    def test_large_rigid_rotation_produces_zero_forces(self):
        """Rotación rígida 90° alrededor del nodo 1, con `rz = π/2` en ambos
        nodos (la sección rota solidariamente con el frame)."""
        theta = math.pi / 2.0
        # Posiciones finales: nodo 1 en (0,0), nodo 2 en (0, 1).
        # DOFs: [ux1, uy1, rz1, ux2, uy2, rz2].
        u_e = np.array([
            0.0, 0.0, theta,                                   # nodo 1: rota in situ
            math.cos(theta) - 1.0, math.sin(theta), theta,     # nodo 2
        ])
        _, F_int = self.elem.compute_element_state(u_e)
        np.testing.assert_allclose(F_int, 0.0, atol=1.0e-7,
            err_msg=f"Frame2DEulerCorot rotación 90°: "
                    f"max|F_int|={np.max(np.abs(F_int)):.3e}; corot no "
                    "extrae la rotación rígida.")


if __name__ == '__main__':
    unittest.main()
