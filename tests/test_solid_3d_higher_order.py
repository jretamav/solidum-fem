"""Tests unitarios de la familia 3D de orden superior (ADR 0012, sub-etapa A.ter).

Cubre los elementos sólidos 3D de orden cuadrático que entran en A.ter
sobre el subsistema 3D lineal de la Etapa 7:

- ``Hex20`` (sub-fase 1, 2026-05-27).
- ``Hex27`` (sub-fase 2, 2026-05-27 — dispara la centralización en
  :class:`_HigherOrderSolid3D`).
- ``Tet10`` se añadirá en sub-fase 3 con su propio mixin (caras
  triangulares Tri6).

Espejo 3D de ``tests/test_higher_order_solid_2d.py`` (familia
Quad8/Quad9/Tri6).

Para evitar duplicación, los tests comunes a `Hex20` y `Hex27` viven en
:class:`_HexQuadraticElementMixin`. Los tests específicos al espacio
de aproximación (patch test triquadrático completo, conteo exacto de
modos de hourglass) van en clases concretas.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.core.domain import Domain
from solidum.elements.solid_3d import Hex20, Hex27, Tet10
from solidum.materials.elastic_3d import Elastic3D


# =============================================================================
# Helpers
# =============================================================================

HEX20_UNIT_REF = [
    # Vértices (orden Hex8)
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    # Medios de aristas de la cara inferior (z=0): 0-1, 1-2, 2-3, 3-0
    (0.5, 0.0, 0.0), (1.0, 0.5, 0.0), (0.5, 1.0, 0.0), (0.0, 0.5, 0.0),
    # Medios de aristas de la cara superior (z=1): 4-5, 5-6, 6-7, 7-4
    (0.5, 0.0, 1.0), (1.0, 0.5, 1.0), (0.5, 1.0, 1.0), (0.0, 0.5, 1.0),
    # Medios de aristas verticales: 0-4, 1-5, 2-6, 3-7
    (0.0, 0.0, 0.5), (1.0, 0.0, 0.5), (1.0, 1.0, 0.5), (0.0, 1.0, 0.5),
]

HEX27_UNIT_REF = HEX20_UNIT_REF + [
    # Centros de cara: -x, +x, -y, +y, -z, +z
    (0.0, 0.5, 0.5), (1.0, 0.5, 0.5),
    (0.5, 0.0, 0.5), (0.5, 1.0, 0.5),
    (0.5, 0.5, 0.0), (0.5, 0.5, 1.0),
    # Centro de cuerpo
    (0.5, 0.5, 0.5),
]


def _build_hex(element_cls, ref_coords, material, coords=None,
               quadrature: str | None = None):
    dom = Domain()
    if coords is None:
        coords = ref_coords
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
    kwargs = {} if quadrature is None else {"quadrature": quadrature}
    elem = element_cls(1, nodes, material, **kwargs)
    dom.add_element(elem)
    return elem, dom


def _apply_linear_field(elem, A: np.ndarray) -> np.ndarray:
    """u_e correspondiente a ``u_i = A_ij·x_j``."""
    n_nodes = len(elem.nodes)
    u_e = np.zeros(n_nodes * 3)
    for i, n in enumerate(elem.nodes):
        x = np.asarray(n.coordinates[:3], dtype=np.float64)
        u_e[3 * i:3 * i + 3] = A @ x
    return u_e


def _assign_sequential_dofs(elem):
    """Asigna índices globales secuenciales a los DOFs del elemento."""
    for i, n in enumerate(elem.nodes):
        n.dofs['ux'] = 3 * i
        n.dofs['uy'] = 3 * i + 1
        n.dofs['uz'] = 3 * i + 2


# =============================================================================
# Mixin compartido por Hex20 y Hex27
# =============================================================================

class _HexQuadraticElementMixin:
    """Tests comunes a los hexaedros 3D cuadráticos (Hex20, Hex27)."""

    # Cada subclase declara:
    ELEMENT_CLS = None       # clase del elemento (Hex20 o Hex27)
    REF_COORDS = ()          # coordenadas de referencia del cubo unitario
    N_NODES = 0
    N_DOFS = 0               # = 3 * N_NODES
    N_GAUSS_DEFAULT = 27     # 3×3×3 default para ambos
    REDUCED_QUADRATURE = "hex_2x2x2"

    def _build(self, material=None, coords=None, quadrature: str | None = None):
        if material is None:
            material = self.material
        return _build_hex(self.ELEMENT_CLS, self.REF_COORDS, material,
                          coords=coords, quadrature=quadrature)

    def setUp(self):
        self.material = Elastic3D(E=2.0e5, nu=0.3)
        self.elem, self.dom = self._build()

    def test_dimensiones_y_dofs(self):
        self.assertEqual(self.elem.STRAIN_DIM, 6)
        self.assertEqual(self.elem.DOF_NAMES, ['ux', 'uy', 'uz'])
        self.assertEqual(self.elem.N_INTEGRATION_POINTS, self.N_GAUSS_DEFAULT)
        K = self.elem.compute_global_stiffness()
        self.assertEqual(K.shape, (self.N_DOFS, self.N_DOFS))

    def test_simetria_K(self):
        K, _ = self.elem.compute_element_state(np.zeros(self.N_DOFS))
        max_K = float(np.max(np.abs(K)))
        np.testing.assert_allclose(K, K.T, atol=1.0e-9 * max_K)

    def test_patch_test_lineal_strain_constante(self):
        eps0 = 1.0e-3
        nu = self.material.nu
        A = np.diag([eps0, -nu * eps0, -nu * eps0])
        u_e = _apply_linear_field(self.elem, A)
        _, F_int = self.elem.compute_element_state(u_e)
        F_nodal = F_int.reshape(self.N_NODES, 3)
        np.testing.assert_allclose(F_nodal.sum(axis=0), 0.0, atol=1.0e-9)

        _assign_sequential_dofs(self.elem)
        gs = self.elem.compute_gauss_state(u_e)
        eps_expected = np.array([eps0, -nu * eps0, -nu * eps0, 0.0, 0.0, 0.0])
        for k in range(self.elem.N_INTEGRATION_POINTS):
            np.testing.assert_allclose(
                gs['strain'][k], eps_expected, atol=1.0e-12,
                err_msg=f"strain incorrecto en Gauss point {k}",
            )

    def test_patch_test_cuadratico_simple(self):
        """``u_x = c·x²`` debe dar ``ε_xx(x) = 2c·x`` exacto en todos los Gauss.

        Cuadrático simple en una sola variable — reproducido tanto por
        Hex20 serendípito como por Hex27 lagrangiano.
        """
        c = 1.0e-3
        u_e = np.zeros(self.N_DOFS)
        for i, n in enumerate(self.elem.nodes):
            x = n.coordinates[0]
            u_e[3 * i] = c * x * x

        _assign_sequential_dofs(self.elem)
        gs = self.elem.compute_gauss_state(u_e)
        for k in range(self.elem.N_INTEGRATION_POINTS):
            x_g = gs['points_global'][k, 0]
            self.assertAlmostEqual(
                gs['strain'][k, 0], 2.0 * c * x_g, places=12,
                msg=f"ε_xx en Gauss point {k} (x_g={x_g})",
            )
            for j in (1, 2, 3, 4, 5):
                self.assertAlmostEqual(gs['strain'][k, j], 0.0, places=12)

    def test_jacobiano_negativo_aborta(self):
        bad = list(self.REF_COORDS)
        swap_z = {0.0: 1.0, 1.0: 0.0}
        bad_swapped = []
        for c in bad:
            x, y, z = c
            new_z = swap_z.get(z, z)  # mids con z=0.5 invariantes
            bad_swapped.append((x, y, new_z))
        elem_bad, _ = self._build(coords=bad_swapped)
        with self.assertRaises(ValueError):
            elem_bad.compute_element_state(np.zeros(self.N_DOFS))

    def test_body_load_uniforme_balance(self):
        f = self.elem.compute_body_load(np.array([0.0, 0.0, -1.0]))
        self.assertAlmostEqual(f[0::3].sum(), 0.0, places=10)
        self.assertAlmostEqual(f[1::3].sum(), 0.0, places=10)
        self.assertAlmostEqual(f[2::3].sum(), -1.0, places=10)

    def test_face_traction_balance(self):
        """Tracción uniforme en cara +ξ debe sumar la presión total."""
        p = 7.5
        t_vec = np.array([p, 0.0, 0.0])
        f = self.elem.compute_face_traction(3, t_vec)
        self.assertAlmostEqual(f[0::3].sum(), p, places=10)
        self.assertAlmostEqual(f[1::3].sum(), 0.0, places=10)
        self.assertAlmostEqual(f[2::3].sum(), 0.0, places=10)
        # Nodos fuera de la cara reciben cero.
        face_local = set(self.elem.FACE_NODES[3])
        for i in range(self.N_NODES):
            if i not in face_local:
                mag = float(np.linalg.norm(f[3 * i:3 * i + 3]))
                self.assertAlmostEqual(
                    mag, 0.0, places=12,
                    msg=f"nodo {i} fuera de face 3 recibe carga: {mag}",
                )

    def test_face_traction_indice_fuera_de_rango(self):
        with self.assertRaises(ValueError):
            self.elem.compute_face_traction(6, np.array([1.0, 0.0, 0.0]))
        with self.assertRaises(ValueError):
            self.elem.compute_face_traction(-1, np.array([1.0, 0.0, 0.0]))

    def test_mass_matrix_consistent_total(self):
        m = Elastic3D(E=1.0, nu=0.3, density=7850.0)
        elem, _ = self._build(material=m)
        M = elem.compute_mass_matrix()
        self.assertEqual(M.shape, (self.N_DOFS, self.N_DOFS))
        np.testing.assert_allclose(M, M.T, atol=1.0e-9)
        self.assertAlmostEqual(M.sum(), 3.0 * 7850.0 * 1.0, places=4)

    def test_mass_matrix_lumped_HRZ(self):
        m = Elastic3D(E=1.0, nu=0.3, density=7850.0)
        elem, _ = self._build(material=m)
        M = elem.compute_mass_matrix(lumping="lumped")
        off = M - np.diag(np.diag(M))
        np.testing.assert_allclose(off, 0.0, atol=1.0e-12)
        self.assertAlmostEqual(np.diag(M).sum(), 3.0 * 7850.0 * 1.0, places=4)
        self.assertTrue(
            float(np.min(np.diag(M))) > 0.0,
            msg=f"HRZ produjo masa no positiva: min={np.min(np.diag(M))}",
        )

    def test_mass_matrix_independent_of_K_quadrature(self):
        m = Elastic3D(E=1.0, nu=0.3, density=7850.0)
        elem_full, _ = self._build(material=m, quadrature="hex_3x3x3")
        elem_red, _ = self._build(material=m, quadrature=self.REDUCED_QUADRATURE)
        M_full = elem_full.compute_mass_matrix()
        M_red = elem_red.compute_mass_matrix()
        np.testing.assert_allclose(M_full, M_red, atol=1.0e-12)

    def test_gauss_state_shapes(self):
        _assign_sequential_dofs(self.elem)
        gs = self.elem.compute_gauss_state(np.zeros(self.N_DOFS))
        n_g = self.elem.N_INTEGRATION_POINTS
        self.assertEqual(gs['strain'].shape, (n_g, 6))
        self.assertEqual(gs['stress'].shape, (n_g, 6))
        self.assertEqual(gs['points_natural'].shape, (n_g, 3))
        self.assertEqual(gs['points_global'].shape, (n_g, 3))
        np.testing.assert_allclose(gs['strain'], 0.0, atol=1.0e-14)
        np.testing.assert_allclose(gs['stress'], 0.0, atol=1.0e-14)


# =============================================================================
# Hex20 — concreto
# =============================================================================

class TestHex20Element(_HexQuadraticElementMixin, unittest.TestCase):
    ELEMENT_CLS = Hex20
    REF_COORDS = HEX20_UNIT_REF
    N_NODES = 20
    N_DOFS = 60

    def test_face_traction_serendipity_share(self):
        """Reparto vértice/medio sobre cara plana cuadrada serendípita.

        Específico al Hex20 — el Hex27 con Quad9 face tiene reparto
        diferente (incluye el nodo central de cara).
        """
        p = 7.5
        t_vec = np.array([p, 0.0, 0.0])
        f = self.elem.compute_face_traction(3, t_vec)
        for v in (1, 2, 5, 6):
            self.assertAlmostEqual(
                f[3 * v], -p / 12.0, places=10,
                msg=f"vértice {v} reparto incorrecto: {f[3*v]}",
            )
        for m_node in (9, 13, 17, 18):
            self.assertAlmostEqual(
                f[3 * m_node], 4.0 * p / 12.0, places=10,
                msg=f"medio {m_node} reparto incorrecto: {f[3*m_node]}",
            )

    def test_reduced_integration_introduces_hourglass_hex20(self):
        """Hex20 con hex_2x2x2: 12 modos null (6 rígidos + 6 hourglass)."""
        elem_red, _ = self._build(quadrature="hex_2x2x2")
        K, _ = elem_red.compute_element_state(np.zeros(60))
        eigs = np.sort(np.abs(np.linalg.eigvalsh(0.5 * (K + K.T))))
        scale = float(eigs[-1])
        threshold = 1.0e-6 * scale
        n_zeros = int(np.sum(eigs < threshold))
        self.assertEqual(
            n_zeros, 12,
            f"Hex20 reducido: {n_zeros} modos null (esperaba 12 = 6+6).",
        )


# =============================================================================
# Hex27 — concreto
# =============================================================================

class TestHex27Element(_HexQuadraticElementMixin, unittest.TestCase):
    ELEMENT_CLS = Hex27
    REF_COORDS = HEX27_UNIT_REF
    N_NODES = 27
    N_DOFS = 81

    def test_patch_test_triquadratico_completo(self):
        """``u_x = c·x²·y²·z²`` — capacidad distintiva del Hex27.

        El espacio Q_2 lagrangiano completo contiene los términos puros
        triquadráticos que el Hex20 serendípito **no** captura. Este test
        falla para Hex20 (overshoot O(c)) y pasa exacto para Hex27.
        """
        c = 1.0e-3
        u_e = np.zeros(81)
        for i, n in enumerate(self.elem.nodes):
            x, y, z = n.coordinates
            u_e[3 * i] = c * (x ** 2) * (y ** 2) * (z ** 2)

        _assign_sequential_dofs(self.elem)
        gs = self.elem.compute_gauss_state(u_e)
        for k in range(self.elem.N_INTEGRATION_POINTS):
            x_g, y_g, z_g = gs['points_global'][k]
            # ε_xx = ∂u_x/∂x = 2c·x·y²·z²
            exp_eps_xx = 2.0 * c * x_g * (y_g ** 2) * (z_g ** 2)
            self.assertAlmostEqual(
                gs['strain'][k, 0], exp_eps_xx, places=12,
                msg=f"ε_xx triquadrático incorrecto en Gauss point {k} "
                    f"(x_g={x_g}, y_g={y_g}, z_g={z_g})",
            )

    def test_reduced_integration_introduces_hourglass_hex27(self):
        """Hex27 con hex_2x2x2: 33 modos null (6 rígidos + 27 hourglass).

        Es el caso más problemático del catálogo 3D: el Lagrangiano
        completo es subintegrado severamente por 2×2×2.
        """
        elem_red, _ = self._build(quadrature="hex_2x2x2")
        K, _ = elem_red.compute_element_state(np.zeros(81))
        eigs = np.sort(np.abs(np.linalg.eigvalsh(0.5 * (K + K.T))))
        scale = float(eigs[-1])
        threshold = 1.0e-6 * scale
        n_zeros = int(np.sum(eigs < threshold))
        self.assertEqual(
            n_zeros, 33,
            f"Hex27 reducido: {n_zeros} modos null (esperaba 33 = 6+27).",
        )


# =============================================================================
# Tet10 — concreto (cara triangular Tri6, cuadratura tet_4 + tet_15 mass)
# =============================================================================

TET10_REF_COORDS = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
    (0.5, 0.0, 0.0),    # 4: medio 0-1
    (0.5, 0.5, 0.0),    # 5: medio 1-2
    (0.0, 0.5, 0.0),    # 6: medio 2-0
    (0.0, 0.0, 0.5),    # 7: medio 0-3
    (0.5, 0.0, 0.5),    # 8: medio 1-3
    (0.0, 0.5, 0.5),    # 9: medio 2-3
]


class TestTet10Element(unittest.TestCase):
    """Tests del Tet10. Geometría triangular incompatible con la mixin
    hexaédrica (volumen 1/6 vs 1, cara triangular vs cuadrilateral), así
    que sus tests viven aquí en una clase concreta — paritarios en
    estructura con `_HexQuadraticElementMixin` pero adaptados al tet.
    """

    REF_VOLUME = 1.0 / 6.0

    def setUp(self):
        self.material = Elastic3D(E=2.0e5, nu=0.3)
        self.dom = Domain()
        nodes = [self.dom.add_node(i + 1, list(c)) for i, c in enumerate(TET10_REF_COORDS)]
        self.elem = Tet10(1, nodes, self.material)
        self.dom.add_element(self.elem)

    def test_dimensiones_y_dofs(self):
        self.assertEqual(self.elem.STRAIN_DIM, 6)
        self.assertEqual(self.elem.DOF_NAMES, ['ux', 'uy', 'uz'])
        self.assertEqual(self.elem.N_INTEGRATION_POINTS, 4)
        K = self.elem.compute_global_stiffness()
        self.assertEqual(K.shape, (30, 30))

    def test_simetria_K(self):
        K, _ = self.elem.compute_element_state(np.zeros(30))
        max_K = float(np.max(np.abs(K)))
        np.testing.assert_allclose(K, K.T, atol=1.0e-9 * max_K)

    def test_patch_test_lineal_strain_constante(self):
        eps0 = 1.0e-3
        nu = self.material.nu
        A = np.diag([eps0, -nu * eps0, -nu * eps0])
        u_e = _apply_linear_field(self.elem, A)
        _, F_int = self.elem.compute_element_state(u_e)
        F_nodal = F_int.reshape(10, 3)
        np.testing.assert_allclose(F_nodal.sum(axis=0), 0.0, atol=1.0e-9)

        _assign_sequential_dofs(self.elem)
        gs = self.elem.compute_gauss_state(u_e)
        eps_expected = np.array([eps0, -nu * eps0, -nu * eps0, 0.0, 0.0, 0.0])
        for k in range(self.elem.N_INTEGRATION_POINTS):
            np.testing.assert_allclose(gs['strain'][k], eps_expected, atol=1.0e-12)

    def test_patch_test_cuadratico(self):
        """``u_x = c·x²`` da ``ε_xx(x) = 2c·x`` exacto en todos los Gauss."""
        c = 1.0e-3
        u_e = np.zeros(30)
        for i, n in enumerate(self.elem.nodes):
            x = n.coordinates[0]
            u_e[3 * i] = c * x * x

        _assign_sequential_dofs(self.elem)
        gs = self.elem.compute_gauss_state(u_e)
        for k in range(self.elem.N_INTEGRATION_POINTS):
            x_g = gs['points_global'][k, 0]
            self.assertAlmostEqual(gs['strain'][k, 0], 2.0 * c * x_g, places=12)

    def test_body_load_uniforme_balance(self):
        """``b = (0,0,-1)`` sobre tet ref con V=1/6 ⇒ Σf_z = -1/6."""
        f = self.elem.compute_body_load(np.array([0.0, 0.0, -1.0]))
        self.assertAlmostEqual(f[0::3].sum(), 0.0, places=10)
        self.assertAlmostEqual(f[1::3].sum(), 0.0, places=10)
        self.assertAlmostEqual(f[2::3].sum(), -self.REF_VOLUME, places=10)

    def test_face_traction_balance(self):
        """Tracción uniforme sobre cara 0 (opuesta al vértice 0)."""
        p = 5.0
        f = self.elem.compute_face_traction(0, np.array([p, 0.0, 0.0]))
        # Cara 0 = triángulo (1,2,3) en (1,0,0)-(0,1,0)-(0,0,1).
        # Área = sqrt(3)/2.
        A_face = np.sqrt(3.0) / 2.0
        self.assertAlmostEqual(f[0::3].sum(), p * A_face, places=10)
        self.assertAlmostEqual(f[1::3].sum(), 0.0, places=10)
        # Nodos fuera de la cara (0, 4, 6, 7) reciben cero.
        face_local = set(self.elem.FACE_NODES[0])
        for i in range(10):
            if i not in face_local:
                mag = float(np.linalg.norm(f[3 * i:3 * i + 3]))
                self.assertAlmostEqual(
                    mag, 0.0, places=12,
                    msg=f"nodo {i} fuera de face 0 recibe carga: {mag}",
                )

    def test_face_traction_indice_fuera_de_rango(self):
        with self.assertRaises(ValueError):
            self.elem.compute_face_traction(4, np.array([1.0, 0.0, 0.0]))
        with self.assertRaises(ValueError):
            self.elem.compute_face_traction(-1, np.array([1.0, 0.0, 0.0]))

    def test_mass_matrix_consistent_with_tet15(self):
        """Σ M = 3·ρ·V_e (conservación de masa total).

        La cuadratura `_MASS_QUADRATURE = "tet_15"` integra exactamente el
        producto cuadrático×cuadrático; la suma total debe coincidir con
        3·ρ·V a precisión máquina.
        """
        m = Elastic3D(E=1.0, nu=0.3, density=7850.0)
        dom = Domain()
        nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(TET10_REF_COORDS)]
        elem = Tet10(1, nodes, m)
        M = elem.compute_mass_matrix()
        self.assertEqual(M.shape, (30, 30))
        np.testing.assert_allclose(M, M.T, atol=1.0e-9)
        self.assertAlmostEqual(M.sum(), 3.0 * 7850.0 * self.REF_VOLUME, places=6)

    def test_mass_matrix_lumped_HRZ(self):
        m = Elastic3D(E=1.0, nu=0.3, density=7850.0)
        dom = Domain()
        nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(TET10_REF_COORDS)]
        elem = Tet10(1, nodes, m)
        M = elem.compute_mass_matrix(lumping="lumped")
        off = M - np.diag(np.diag(M))
        np.testing.assert_allclose(off, 0.0, atol=1.0e-12)
        self.assertAlmostEqual(np.diag(M).sum(), 3.0 * 7850.0 * self.REF_VOLUME, places=6)
        self.assertTrue(
            float(np.min(np.diag(M))) > 0.0,
            msg=f"HRZ Tet10 produjo masa no positiva: min={np.min(np.diag(M))}",
        )

    def test_gauss_state_shapes(self):
        _assign_sequential_dofs(self.elem)
        gs = self.elem.compute_gauss_state(np.zeros(30))
        self.assertEqual(gs['strain'].shape, (4, 6))
        self.assertEqual(gs['stress'].shape, (4, 6))
        self.assertEqual(gs['points_natural'].shape, (4, 3))
        self.assertEqual(gs['points_global'].shape, (4, 3))
        np.testing.assert_allclose(gs['strain'], 0.0, atol=1.0e-14)


if __name__ == "__main__":
    unittest.main()
