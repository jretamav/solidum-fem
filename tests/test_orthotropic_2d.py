"""Elasticidad ortótropa 2D — `Orthotropic2D` (plane stress).

Cubre los criterios de `acceptance` de docs/specs/Orthotropic2D.md:

- Degeneración al isótropo (el test de regresión más valioso: el isótropo es
  el caso particular del ortótropo).
- Ausencia de acoplamiento en ejes del material; aparición del acoplamiento
  tracción-cortante al rotar.
- Invariancia bajo rotación simultánea de material y deformación, que blinda
  el factor 2 del cortante engineering (ver ADR 0013 y §12 de la spec).
- Admisibilidad ortótropa, que **no** es la isótropa: nu12 > 0.5 es legítimo.

Órdenes de magnitud de bambú tomados de literatura general (E1/E2 ~ 20),
no de una especie concreta — ver §Diálogo de la spec.
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.orthotropic_2d import Orthotropic2D


# Datos representativos de bambú (literatura general, no especie concreta).
BAMBU = dict(E1=15.0e9, E2=0.8e9, G12=0.7e9, nu12=0.35)


def _T_eps(theta_deg: float) -> np.ndarray:
    """Transformación de deformaciones Voigt engineering (gamma = 2*eps)."""
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    return np.array([
        [c * c,          s * s,          c * s],
        [s * s,          c * c,         -c * s],
        [-2.0 * c * s,   2.0 * c * s,    c * c - s * s],
    ])


def _T_sig(theta_deg: float) -> np.ndarray:
    """Transformación de esfuerzos Voigt.

    **No coincide** con `_T_eps`: difieren en factores 2 en el bloque cortante
    por la convención engineering. Confundirlas es el error clásico.
    """
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    return np.array([
        [c * c,      s * s,     2.0 * c * s],
        [s * s,      c * c,    -2.0 * c * s],
        [-c * s,     c * s,     c * c - s * s],
    ])


def _acoplamiento_relativo(C: np.ndarray) -> float:
    """Términos de acoplamiento normalizados por la escala de C.

    La comparación debe ser **adimensional**: C tiene magnitudes ~1e10 Pa, de
    modo que un cero exacto en aritmética de punto flotante deja residuos de
    ~1e-6 en valor absoluto que son 1e-16 en relativo. Comparar contra cero
    absoluto convertiría el test en una medida de la magnitud de E1, no de la
    física. Patrón canónico del proyecto (ADR 0006).
    """
    escala = np.abs(C).max()
    return max(abs(C[0, 2]), abs(C[1, 2])) / escala


class TestDegeneracionAlIsotropo(unittest.TestCase):
    """El isótropo es el caso particular del ortótropo."""

    def test_C_coincide_con_elastic2d_plane_stress(self):
        E, nu = 210.0e9, 0.3
        G = E / (2.0 * (1.0 + nu))

        orto = Orthotropic2D(E1=E, E2=E, G12=G, nu12=nu, theta=0.0)
        iso = Elastic2D(E=E, nu=nu, hypothesis='plane_stress')

        np.testing.assert_allclose(orto.C, iso.C, rtol=1e-14, atol=0.0)

    def test_degeneracion_se_mantiene_bajo_rotacion(self):
        """Un isótropo rotado sigue siendo el mismo isótropo."""
        E, nu = 210.0e9, 0.3
        G = E / (2.0 * (1.0 + nu))
        iso = Elastic2D(E=E, nu=nu, hypothesis='plane_stress')

        for theta in (0.0, 17.0, 45.0, 90.0, 123.4):
            with self.subTest(theta=theta):
                orto = Orthotropic2D(E1=E, E2=E, G12=G, nu12=nu, theta=theta)
                escala = np.abs(iso.C).max()
                self.assertLess(np.abs(orto.C - iso.C).max() / escala, 1e-14)

    def test_reciprocidad_degenera_a_nu_unico(self):
        E, nu = 210.0e9, 0.3
        G = E / (2.0 * (1.0 + nu))
        orto = Orthotropic2D(E1=E, E2=E, G12=G, nu12=nu)
        self.assertAlmostEqual(orto.nu21, nu, places=14)


class TestEjesDelMaterial(unittest.TestCase):
    """En ejes del material no hay acoplamiento."""

    def test_bloque_cortante_nulo_en_theta_0(self):
        m = Orthotropic2D(**BAMBU, theta=0.0)
        self.assertAlmostEqual(m.C[0, 2], 0.0, places=14)
        self.assertAlmostEqual(m.C[1, 2], 0.0, places=14)

    def test_bloque_cortante_nulo_en_theta_90(self):
        m = Orthotropic2D(**BAMBU, theta=90.0)
        self.assertLess(_acoplamiento_relativo(m.C), 1e-15)

    def test_theta_90_intercambia_E1_y_E2(self):
        """A 90 grados los papeles de las dos direcciones se invierten."""
        m0 = Orthotropic2D(**BAMBU, theta=0.0)
        m90 = Orthotropic2D(**BAMBU, theta=90.0)
        escala = np.abs(m0.C).max()
        self.assertLess(abs(m0.C[0, 0] - m90.C[1, 1]) / escala, 1e-15)
        self.assertLess(abs(m0.C[1, 1] - m90.C[0, 0]) / escala, 1e-15)

    def test_traccion_pura_no_distorsiona_en_ejes(self):
        m = Orthotropic2D(**BAMBU, theta=0.0)
        eps = np.array([1.0e-4, 0.0, 0.0])
        sigma, _, _ = m.compute_state(eps)
        self.assertLess(abs(sigma[2]) / np.abs(sigma).max(), 1e-15)

    def test_rigidez_longitudinal_domina(self):
        """E1 >> E2 debe reflejarse en la constitutiva."""
        m = Orthotropic2D(**BAMBU, theta=0.0)
        self.assertGreater(m.C[0, 0], 10.0 * m.C[1, 1])


class TestAcoplamientoTraccionCortante(unittest.TestCase):
    """Al rotar aparece acoplamiento inexistente en isótropos."""

    def test_theta_45_genera_cortante_bajo_traccion_pura(self):
        m = Orthotropic2D(**BAMBU, theta=45.0)
        eps = np.array([1.0e-4, 0.0, 0.0])
        sigma, _, _ = m.compute_state(eps)
        self.assertGreater(abs(sigma[2]), 0.0)

    def test_C_glob_contra_calculo_independiente(self):
        """Verifica T^T C T contra una construcción hecha aparte en el test."""
        theta = 30.0
        m = Orthotropic2D(**BAMBU, theta=theta)

        nu21 = BAMBU['nu12'] * BAMBU['E2'] / BAMBU['E1']
        d = 1.0 - BAMBU['nu12'] * nu21
        C_mat = np.array([
            [BAMBU['E1'] / d,                 BAMBU['nu12'] * BAMBU['E2'] / d, 0.0],
            [BAMBU['nu12'] * BAMBU['E2'] / d, BAMBU['E2'] / d,                 0.0],
            [0.0,                             0.0,               BAMBU['G12']],
        ])
        T = _T_eps(theta)
        ref = T.T @ C_mat @ T
        self.assertLess(np.abs(m.C - ref).max() / np.abs(ref).max(), 1e-14)

    def test_acoplamiento_se_anula_en_ejes_principales(self):
        for theta in (0.0, 90.0, 180.0):
            with self.subTest(theta=theta):
                m = Orthotropic2D(**BAMBU, theta=theta)
                self.assertLess(_acoplamiento_relativo(m.C), 1e-15)


class TestInvarianciaRotacional(unittest.TestCase):
    """Rotar material y deformación juntos no cambia la física.

    Es el test que blinda el factor 2 del cortante engineering: si se
    confundieran T_eps y T_sig, este test falla.
    """

    def test_sigma_rota_consistentemente(self):
        eps = np.array([3.0e-4, -1.0e-4, 2.0e-4])

        for theta in (0.0, 23.0, 45.0, 67.5, 90.0):
            with self.subTest(theta=theta):
                m0 = Orthotropic2D(**BAMBU, theta=0.0)
                sigma0, _, _ = m0.compute_state(eps)

                # Material rotado +theta, deformación rotada +theta.
                m_rot = Orthotropic2D(**BAMBU, theta=theta)
                eps_rot = np.linalg.inv(_T_eps(theta)) @ eps
                sigma_rot, _, _ = m_rot.compute_state(eps_rot)

                ref = np.linalg.inv(_T_sig(theta)) @ sigma0
                self.assertLess(
                    np.abs(sigma_rot - ref).max() / np.abs(sigma0).max(),
                    1e-12,
                )

    def test_energia_de_deformacion_invariante(self):
        eps = np.array([3.0e-4, -1.0e-4, 2.0e-4])
        m0 = Orthotropic2D(**BAMBU, theta=0.0)
        W0 = 0.5 * eps @ (m0.C @ eps)

        for theta in (15.0, 45.0, 80.0):
            with self.subTest(theta=theta):
                m = Orthotropic2D(**BAMBU, theta=theta)
                eps_rot = np.linalg.inv(_T_eps(theta)) @ eps
                W = 0.5 * eps_rot @ (m.C @ eps_rot)
                self.assertAlmostEqual(W / W0, 1.0, places=10)


class TestPropiedadesDeC(unittest.TestCase):

    def test_simetrica_para_todo_theta(self):
        for theta in np.linspace(0.0, 180.0, 19):
            with self.subTest(theta=theta):
                m = Orthotropic2D(**BAMBU, theta=float(theta))
                self.assertLess(
                    np.abs(m.C - m.C.T).max() / np.abs(m.C).max(), 1e-15)

    def test_definida_positiva_para_todo_theta(self):
        for theta in np.linspace(0.0, 180.0, 19):
            with self.subTest(theta=theta):
                m = Orthotropic2D(**BAMBU, theta=float(theta))
                self.assertGreater(np.linalg.eigvalsh(m.C).min(), 0.0)

    def test_reciprocidad(self):
        m = Orthotropic2D(**BAMBU)
        self.assertAlmostEqual(m.nu12 / m.E1, m.nu21 / m.E2, places=20)

    def test_nu21_pequeno_en_material_muy_ortotropo(self):
        """E1 >> E2 obliga a nu21 << nu12. No es un error de medición."""
        m = Orthotropic2D(**BAMBU)
        self.assertLess(m.nu21, 0.03)

    def test_tangente_constante_sin_historia(self):
        m = Orthotropic2D(**BAMBU, theta=30.0)
        _, C1, _ = m.compute_state(np.array([1.0e-4, 0.0, 0.0]))
        _, C2, _ = m.compute_state(np.array([0.0, -5.0e-4, 3.0e-4]))
        np.testing.assert_allclose(C1, C2, rtol=1e-14, atol=0.0)


class TestAdmisibilidad(unittest.TestCase):
    """El criterio ortótropo NO es el isótropo."""

    def test_nu12_mayor_que_medio_es_admisible(self):
        """Con E1/E2 = 20 el límite es sqrt(20) = 4.47, no 0.5.

        Blinda que no se copió la validación de `Elastic2D`.
        """
        m = Orthotropic2D(E1=15.0e9, E2=0.75e9, G12=0.7e9, nu12=0.6)
        self.assertGreater(np.linalg.eigvalsh(m.C).min(), 0.0)

    def test_rechaza_nu12_sobre_el_limite(self):
        limite = math.sqrt(15.0e9 / 0.8e9)
        with self.assertRaises(ValueError) as ctx:
            Orthotropic2D(E1=15.0e9, E2=0.8e9, G12=0.7e9, nu12=limite + 0.1)
        self.assertIn("nu12", str(ctx.exception))

    def test_rechaza_modulos_no_positivos(self):
        for campo in ("E1", "E2", "G12"):
            with self.subTest(campo=campo):
                datos = dict(BAMBU)
                datos[campo] = 0.0
                with self.assertRaises(ValueError) as ctx:
                    Orthotropic2D(**datos)
                self.assertIn(campo, str(ctx.exception))

    def test_rechaza_densidad_negativa(self):
        with self.assertRaises(ValueError):
            Orthotropic2D(**BAMBU, density=-1.0)

    def test_density_opcional(self):
        """ADR 0008: opcional al construir, no obligatoria."""
        self.assertIsNone(Orthotropic2D(**BAMBU).density)


class TestContratoMaterial(unittest.TestCase):

    def test_declaraciones_de_clase(self):
        self.assertEqual(Orthotropic2D.STRAIN_DIM, 3)
        self.assertIsNone(Orthotropic2D.PRIMARY_STATE_VAR)
        self.assertTrue(Orthotropic2D.IS_SYMMETRIC)

    def test_registrado_en_el_registry(self):
        import solidum  # noqa: F401  (dispara el autodiscover)
        from solidum.registry import MaterialRegistry
        self.assertIn("Orthotropic2D", MaterialRegistry.names())

    def test_compute_state_preserva_state_vars(self):
        m = Orthotropic2D(**BAMBU)
        _, _, sv = m.compute_state(np.array([1.0e-4, 0.0, 0.0]), None)
        self.assertIsNone(sv)


if __name__ == "__main__":
    unittest.main()
