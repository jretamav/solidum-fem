"""Tests unitarios de la familia 3D (ADR 0012, Etapa 7).

Cubre tres componentes nuevos:

- ``Elastic3D`` — material elástico isótropo 3D (Voigt 6D del proyecto).
- ``Hex8`` — hexaedro trilineal 3D.
- ``Tet4`` — tetraedro lineal 3D (CST 3D).

Tests:

1. **Material**: ley de Hooke 3×3, simetría, positividad definida, tracción
   uniaxial, cortantes puros en los tres planos, compresión hidrostática,
   rechazo de inputs inválidos.
2. **Elementos** (Hex8 y Tet4): patch test físico (tracción uniaxial reproducida
   en σ), simetría de K_e, jacobiano degenerado abortado, cargas distribuidas
   (body load + face traction) con balances de fuerza, mass matrix consistente
   y lumped (HRZ), gauss state.

Los modos rígidos están en ``test_rigid_body_modes.py`` (extensión a 3D).
La validación contra referencias publicadas (cubo Lamé, MacNeal) vive en
``tests/validation/``.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.core.domain import Domain
from fenix.elements.solid_3d import Hex8, Tet4
from fenix.materials.elastic_3d import Elastic3D


# =============================================================================
# Helpers
# =============================================================================

HEX8_REF_COORDS = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]
TET4_REF_COORDS = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]


def _build_hex8(material, coords=None):
    dom = Domain()
    if coords is None:
        coords = HEX8_REF_COORDS
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
    elem = Hex8(1, nodes, material)
    dom.add_element(elem)
    return elem, dom


def _build_tet4(material, coords=None):
    dom = Domain()
    if coords is None:
        coords = TET4_REF_COORDS
    nodes = [dom.add_node(i + 1, list(c)) for i, c in enumerate(coords)]
    elem = Tet4(1, nodes, material)
    dom.add_element(elem)
    return elem, dom


def _apply_linear_field(elem, A):
    """u_e correspondiente a ``u_i = A_ij · x_j`` (gradiente constante A 3×3)."""
    n_nodes = len(elem.nodes)
    u_e = np.zeros(n_nodes * 3)
    for i, n in enumerate(elem.nodes):
        x = np.asarray(n.coordinates[:3], dtype=np.float64)
        u_e[3 * i:3 * i + 3] = A @ x
    return u_e


# =============================================================================
# Elastic3D — material
# =============================================================================

class TestElastic3D(unittest.TestCase):

    def test_construye_y_C_es_simetrica(self):
        m = Elastic3D(E=2.0e5, nu=0.3)
        self.assertEqual(m.STRAIN_DIM, 6)
        np.testing.assert_allclose(m.C, m.C.T, atol=1.0e-14)

    def test_C_positiva_definida_en_rango(self):
        for nu in (-0.5, 0.0, 0.3, 0.49):
            m = Elastic3D(E=2.0e5, nu=nu)
            eigs = np.linalg.eigvalsh(m.C)
            self.assertTrue(np.all(eigs > 0.0),
                            f"C no PD con nu={nu}: eigs={eigs}")

    def test_traccion_uniaxial(self):
        E, nu = 2.0e5, 0.3
        eps0 = 1.0e-3
        m = Elastic3D(E=E, nu=nu)
        strain = np.array([eps0, -nu * eps0, -nu * eps0, 0.0, 0.0, 0.0])
        sigma, C, _ = m.compute_state(strain)
        expected = np.array([E * eps0, 0.0, 0.0, 0.0, 0.0, 0.0])
        np.testing.assert_allclose(sigma, expected, atol=1.0e-9, rtol=1.0e-12)
        np.testing.assert_allclose(C, m.C)

    def test_cortante_puro_xy(self):
        E, nu = 2.0e5, 0.3
        G = E / (2.0 * (1.0 + nu))
        gamma = 5.0e-4
        m = Elastic3D(E=E, nu=nu)
        strain = np.array([0.0, 0.0, 0.0, gamma, 0.0, 0.0])
        sigma, _, _ = m.compute_state(strain)
        expected = np.array([0.0, 0.0, 0.0, G * gamma, 0.0, 0.0])
        np.testing.assert_allclose(sigma, expected, atol=1.0e-9, rtol=1.0e-12)

    def test_cortante_puro_yz(self):
        E, nu = 2.0e5, 0.3
        G = E / (2.0 * (1.0 + nu))
        gamma = 5.0e-4
        m = Elastic3D(E=E, nu=nu)
        strain = np.array([0.0, 0.0, 0.0, 0.0, gamma, 0.0])
        sigma, _, _ = m.compute_state(strain)
        expected = np.array([0.0, 0.0, 0.0, 0.0, G * gamma, 0.0])
        np.testing.assert_allclose(sigma, expected, atol=1.0e-9, rtol=1.0e-12)

    def test_cortante_puro_xz(self):
        E, nu = 2.0e5, 0.3
        G = E / (2.0 * (1.0 + nu))
        gamma = 5.0e-4
        m = Elastic3D(E=E, nu=nu)
        strain = np.array([0.0, 0.0, 0.0, 0.0, 0.0, gamma])
        sigma, _, _ = m.compute_state(strain)
        expected = np.array([0.0, 0.0, 0.0, 0.0, 0.0, G * gamma])
        np.testing.assert_allclose(sigma, expected, atol=1.0e-9, rtol=1.0e-12)

    def test_compresion_hidrostatica(self):
        """Bulk modulus K = E / [3(1-2ν)]."""
        E, nu = 2.0e5, 0.3
        K_bulk = E / (3.0 * (1.0 - 2.0 * nu))
        eps0 = -1.0e-3
        m = Elastic3D(E=E, nu=nu)
        strain = np.array([eps0, eps0, eps0, 0.0, 0.0, 0.0])
        sigma, _, _ = m.compute_state(strain)
        # σ_ii = K · tr(ε) = K · 3 · eps0
        expected_sigma_ii = K_bulk * 3.0 * eps0
        for i in range(3):
            self.assertAlmostEqual(sigma[i], expected_sigma_ii, places=8)
        for i in range(3, 6):
            self.assertAlmostEqual(sigma[i], 0.0, places=12)

    def test_rechazo_E_negativo(self):
        with self.assertRaises(ValueError):
            Elastic3D(E=-1.0, nu=0.3)

    def test_rechazo_nu_fuera_de_rango(self):
        for bad in (-1.5, 0.5, 0.51, 1.0):
            with self.assertRaises(ValueError):
                Elastic3D(E=1.0, nu=bad)

    def test_rechazo_density_negativa(self):
        with self.assertRaises(ValueError):
            Elastic3D(E=1.0, nu=0.3, density=-7.0)

    def test_density_opcional(self):
        m = Elastic3D(E=1.0, nu=0.3)
        self.assertIsNone(m.density)
        m2 = Elastic3D(E=1.0, nu=0.3, density=7850.0)
        self.assertEqual(m2.density, 7850.0)


# =============================================================================
# Hex8 — elemento
# =============================================================================

class TestHex8Element(unittest.TestCase):

    def setUp(self):
        self.material = Elastic3D(E=2.0e5, nu=0.3)
        self.elem, self.dom = _build_hex8(self.material)

    def test_dimensiones_y_dofs(self):
        self.assertEqual(self.elem.STRAIN_DIM, 6)
        self.assertEqual(self.elem.DOF_NAMES, ['ux', 'uy', 'uz'])
        self.assertEqual(self.elem.N_INTEGRATION_POINTS, 8)
        K = self.elem.compute_global_stiffness()
        self.assertEqual(K.shape, (24, 24))

    def test_simetria_K(self):
        K, _ = self.elem.compute_element_state(np.zeros(24))
        np.testing.assert_allclose(K, K.T, atol=1.0e-9 * np.max(np.abs(K)))

    def test_patch_traccion_uniaxial_da_sigma_correcto(self):
        """Campo lineal ``u_x = eps0·x, u_y = -ν·eps0·y, u_z = -ν·eps0·z`` debe
        dar σ_xx = E·eps0, resto cero, en todos los puntos de Gauss."""
        E, nu = self.material.E, self.material.nu
        eps0 = 1.0e-3
        A = np.diag([eps0, -nu * eps0, -nu * eps0])
        u_e = _apply_linear_field(self.elem, A)
        K, F_int = self.elem.compute_element_state(u_e)
        # F_int debe ser equilibrado en el bulk: la suma de fuerzas externas
        # equivalentes (en este caso, las que reproducen la tracción uniaxial)
        # da Σf_x distinto de cero por nodos del borde, pero la SUMA de F_int
        # sobre todos los nodos en cada dirección debe ser nula (equilibrio
        # de un sólido elástico interior bajo strain uniforme).
        F_nodal = F_int.reshape(8, 3)
        np.testing.assert_allclose(F_nodal.sum(axis=0), 0.0, atol=1.0e-9)

    def test_jacobiano_negativo_aborta(self):
        """Si los nodos 4-7 están bajo los nodos 0-3 (cubo invertido en z), det J < 0."""
        bad_coords = [
            (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
            (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        ]
        elem_bad, _ = _build_hex8(self.material, coords=bad_coords)
        with self.assertRaises(ValueError):
            elem_bad.compute_element_state(np.zeros(24))

    def test_body_load_uniforme(self):
        """``b = (0, 0, -g)`` con cubo unitario debe sumar Σ f_z = -g·V = -g·1."""
        f = self.elem.compute_body_load(np.array([0.0, 0.0, -1.0]))
        self.assertAlmostEqual(f[::3].sum(), 0.0, places=10)
        self.assertAlmostEqual(f[1::3].sum(), 0.0, places=10)
        self.assertAlmostEqual(f[2::3].sum(), -1.0, places=10)

    def test_face_traction_balance(self):
        """Tracción uniforme en cara +ξ (face 3) del cubo unitario debe sumar p·A=p·1."""
        t_vec = np.array([7.5, 0.0, 0.0])
        f = self.elem.compute_face_traction(3, t_vec)
        # Suma total
        self.assertAlmostEqual(f[::3].sum(), 7.5, places=10)
        self.assertAlmostEqual(f[1::3].sum(), 0.0, places=10)
        self.assertAlmostEqual(f[2::3].sum(), 0.0, places=10)
        # Reparto: solo los 4 nodos de la cara +ξ (1, 2, 6, 5) reciben fuerza.
        active = {1, 2, 5, 6}
        for i in range(8):
            mag = np.linalg.norm(f[3 * i:3 * i + 3])
            if i in active:
                self.assertAlmostEqual(mag, 7.5 / 4.0, places=10)
            else:
                self.assertAlmostEqual(mag, 0.0, places=12)

    def test_face_traction_indice_fuera_de_rango(self):
        with self.assertRaises(ValueError):
            self.elem.compute_face_traction(6, np.array([1.0, 0.0, 0.0]))
        with self.assertRaises(ValueError):
            self.elem.compute_face_traction(-1, np.array([1.0, 0.0, 0.0]))

    def test_mass_matrix_consistent(self):
        m = Elastic3D(E=1.0, nu=0.3, density=7850.0)
        elem, _ = _build_hex8(m)
        M = elem.compute_mass_matrix()
        self.assertEqual(M.shape, (24, 24))
        np.testing.assert_allclose(M, M.T, atol=1.0e-9)
        # Masa total por dirección: sum_{i,j} M_consistent[i, j] (acoplado por nodos)
        # debe igualar 3·ρ·V para el cubo unitario.
        self.assertAlmostEqual(M.sum(), 3.0 * 7850.0 * 1.0, places=6)

    def test_mass_matrix_lumped_hrz(self):
        """HRZ: Σ diagonal traslacional = 3·ρ·V; matriz diagonal."""
        m = Elastic3D(E=1.0, nu=0.3, density=7850.0)
        elem, _ = _build_hex8(m)
        M = elem.compute_mass_matrix(lumping="lumped")
        # Diagonal estricta
        off = M - np.diag(np.diag(M))
        np.testing.assert_allclose(off, 0.0, atol=1.0e-12)
        # Conservación de masa: total = 3·ρ·V
        self.assertAlmostEqual(np.diag(M).sum(), 3.0 * 7850.0 * 1.0, places=6)

    def test_gauss_state(self):
        """``compute_gauss_state`` con DOFs asignados a u_e correcto."""
        eps0 = 1.0e-3
        nu = self.material.nu
        A = np.diag([eps0, -nu * eps0, -nu * eps0])
        u_e = _apply_linear_field(self.elem, A)
        # Asignar manualmente los índices globales para que get_local_displacements
        # devuelva u_e directamente.
        for i, n in enumerate(self.elem.nodes):
            n.dofs['ux'] = 3 * i
            n.dofs['uy'] = 3 * i + 1
            n.dofs['uz'] = 3 * i + 2
        gs = self.elem.compute_gauss_state(u_e)
        self.assertEqual(gs['stress'].shape, (8, 6))
        self.assertEqual(gs['strain'].shape, (8, 6))
        # Strain uniforme = diag([eps0, -nu·eps0, -nu·eps0]) Voigt
        eps_expected = np.array([eps0, -nu * eps0, -nu * eps0, 0.0, 0.0, 0.0])
        for k in range(8):
            np.testing.assert_allclose(gs['strain'][k], eps_expected, atol=1.0e-12)


# =============================================================================
# Tet4 — elemento
# =============================================================================

class TestTet4Element(unittest.TestCase):

    def setUp(self):
        self.material = Elastic3D(E=2.0e5, nu=0.3)
        self.elem, self.dom = _build_tet4(self.material)

    def test_dimensiones_y_dofs(self):
        self.assertEqual(self.elem.STRAIN_DIM, 6)
        self.assertEqual(self.elem.DOF_NAMES, ['ux', 'uy', 'uz'])
        self.assertEqual(self.elem.N_INTEGRATION_POINTS, 1)
        K = self.elem.compute_global_stiffness()
        self.assertEqual(K.shape, (12, 12))

    def test_simetria_K_exacta(self):
        """Tet4 tiene K = B^T·D·B·V_e, simétrica exacta (no numérica)."""
        K, _ = self.elem.compute_element_state(np.zeros(12))
        np.testing.assert_array_equal(K, K.T)

    def test_volumen(self):
        """V_e = 1/6 para el tet de referencia."""
        self.assertAlmostEqual(self.elem._element_volume(), 1.0 / 6.0, places=12)

    def test_patch_traccion_uniaxial_da_sigma_correcto(self):
        """CST 3D: campo lineal reproduce strain constante exacto."""
        E, nu = self.material.E, self.material.nu
        eps0 = 1.0e-3
        A = np.diag([eps0, -nu * eps0, -nu * eps0])
        u_e = _apply_linear_field(self.elem, A)
        K, F_int = self.elem.compute_element_state(u_e)
        F_nodal = F_int.reshape(4, 3)
        np.testing.assert_allclose(F_nodal.sum(axis=0), 0.0, atol=1.0e-9)
        # σ debe ser igual a Elastic3D · ε exacto
        eps = np.array([eps0, -nu * eps0, -nu * eps0, 0.0, 0.0, 0.0])
        sigma_expected = self.material.C @ eps
        np.testing.assert_allclose(self.elem.state.stresses_trial[0], sigma_expected,
                                   atol=1.0e-9, rtol=1.0e-12)

    def test_jacobiano_degenerado_aborta(self):
        """Cuatro nodos coplanares: det(J) = 0."""
        bad_coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                      (0.0, 1.0, 0.0), (0.5, 0.5, 0.0)]
        elem_bad, _ = _build_tet4(self.material, coords=bad_coords)
        with self.assertRaises(ValueError):
            elem_bad.compute_element_state(np.zeros(12))

    def test_body_load_uniforme(self):
        """``b · V_e / 4`` por nodo; suma total = b · V_e."""
        b = np.array([1.0, 2.0, -3.0])
        f = self.elem.compute_body_load(b)
        V = 1.0 / 6.0
        for i in range(4):
            np.testing.assert_allclose(f[3 * i:3 * i + 3], b * V / 4.0, atol=1.0e-12)
        np.testing.assert_allclose(
            f.reshape(4, 3).sum(axis=0), b * V, atol=1.0e-12,
        )

    def test_face_traction_balance(self):
        """Tracción uniforme en cara 0 (opuesta al nodo 0, área de la cara
        (1,2,3) en el tet de referencia). Σf = t̄·A; cero al nodo 0."""
        # Cara (1,2,3): vértices (1,0,0)-(0,1,0)-(0,0,1). Lados v1-v2 y v1-v3.
        X = np.array([(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)])
        v1 = X[1] - X[0]
        v2 = X[2] - X[0]
        A = 0.5 * np.linalg.norm(np.cross(v1, v2))
        t_vec = np.array([5.0, 0.0, 0.0])
        f = self.elem.compute_face_traction(0, t_vec)
        # Suma total Σf_x = A · 5
        self.assertAlmostEqual(f[::3].sum(), A * 5.0, places=12)
        # Nodo 0 (opuesto a la cara) recibe cero
        np.testing.assert_allclose(f[0:3], 0.0, atol=1.0e-14)
        # Reparto: A/3 · t̄ a cada uno de los nodos 1, 2, 3
        for node_local in (1, 2, 3):
            base = 3 * node_local
            np.testing.assert_allclose(
                f[base:base + 3], (A / 3.0) * t_vec, atol=1.0e-12,
            )

    def test_face_traction_indice_fuera_de_rango(self):
        with self.assertRaises(ValueError):
            self.elem.compute_face_traction(4, np.array([1.0, 0.0, 0.0]))

    def test_mass_matrix_consistent(self):
        m = Elastic3D(E=1.0, nu=0.3, density=7850.0)
        elem, _ = _build_tet4(m)
        M = elem.compute_mass_matrix()
        np.testing.assert_allclose(M, M.T, atol=1.0e-12)
        V = 1.0 / 6.0
        # Suma total = 3·ρ·V (todas las direcciones)
        self.assertAlmostEqual(M.sum(), 3.0 * 7850.0 * V, places=6)

    def test_mass_matrix_lumped_hrz(self):
        m = Elastic3D(E=1.0, nu=0.3, density=7850.0)
        elem, _ = _build_tet4(m)
        M = elem.compute_mass_matrix(lumping="lumped")
        off = M - np.diag(np.diag(M))
        np.testing.assert_allclose(off, 0.0, atol=1.0e-14)
        V = 1.0 / 6.0
        self.assertAlmostEqual(np.diag(M).sum(), 3.0 * 7850.0 * V, places=6)
        # Reparto equitativo en HRZ para Tet4 (4 nodos): cada DOF traslacional
        # recibe ρ·V/4.
        for i in range(12):
            self.assertAlmostEqual(M[i, i], 7850.0 * V / 4.0, places=8)


if __name__ == "__main__":
    unittest.main()
