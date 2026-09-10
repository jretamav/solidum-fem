"""Verificación de la infraestructura del estudio (Fase 0).

Las métricas del estudio deciden sus conclusiones, así que deben estar
blindadas antes de producir un solo resultado. Este archivo comprueba que se
comportan como el diseño (`../DISENO_DEL_ESTUDIO.md` §7) especifica.

No es parte de la suite del proyecto: `pytest paper/bambu/estudio/`.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from _comun import (  # noqa: E402
    CENTRO, RANGO_ANISOTROPIA, RANGO_G_NORMALIZADO, RANGO_NU,
    UMBRAL_ADMISIBLE, UMBRAL_RESERVAS, ParametrosBambu, clasificar,
    constantes_aparentes, isotropo_I_L, isotropo_I_T, m1_error_modulo,
    m2_error_constitutiva, m3_acoplamiento_omitido, m4_error_desplazamiento,
    m5_error_esfuerzos,
)


class TestParametros(unittest.TestCase):

    def test_reciprocidad(self):
        p = ParametrosBambu()
        self.assertAlmostEqual(p.nu_LT / p.E_L, p.nu_TL / p.E_T, places=20)

    def test_anisotropia_define_E_T(self):
        p = ParametrosBambu(E_L=20.0e9, anisotropia=8.0)
        self.assertAlmostEqual(p.E_L / p.E_T, 8.0, places=12)

    def test_todo_el_espacio_de_parametros_es_admisible(self):
        """Ningún punto del barrido viola |nu12| < sqrt(E1/E2).

        Si fallara, el barrido tendría agujeros y habría que acotar los rangos
        del §5 — mejor descubrirlo aquí que a mitad de la Fase 1.
        """
        for a in np.linspace(*RANGO_ANISOTROPIA, 15):
            for nu in np.linspace(*RANGO_NU, 9):
                with self.subTest(anisotropia=a, nu=nu):
                    p = ParametrosBambu(anisotropia=float(a), nu_LT=float(nu))
                    self.assertTrue(p.es_admisible())
                    p.ortotropo(0.0)   # no debe lanzar

    def test_centro_dentro_de_los_rangos(self):
        self.assertTrue(
            RANGO_ANISOTROPIA[0] <= CENTRO['anisotropia'] <= RANGO_ANISOTROPIA[1])
        self.assertTrue(
            RANGO_G_NORMALIZADO[0] <= CENTRO['g_normalizado'] <= RANGO_G_NORMALIZADO[1])
        self.assertTrue(RANGO_NU[0] <= CENTRO['nu_LT'] <= RANGO_NU[1])


class TestMetricasConstitutivas(unittest.TestCase):

    def setUp(self):
        self.p = ParametrosBambu()

    def test_M1_nula_para_I_L_en_el_eje(self):
        """I-L calibra E = E_L, así que a theta=0 acierta el módulo aparente
        por construcción. Si esto fallara, la idealización estaría mal
        definida."""
        e = m1_error_modulo(isotropo_I_L(self.p), self.p.ortotropo(0.0))
        self.assertLess(e, 1e-14)

    def test_M1_nula_para_I_T_a_noventa_grados(self):
        """Simétricamente, I-T acierta cuando la fibra es perpendicular."""
        e = m1_error_modulo(isotropo_I_T(self.p), self.p.ortotropo(90.0))
        self.assertLess(e, 1e-14)

    def test_M2_no_nula_aunque_M1_lo_sea(self):
        """Acertar E_x NO implica acertar la constitutiva.

        Es la base de la hipótesis H3: un isótropo calibrado para reproducir
        rigidez en una dirección sigue errando el resto de la matriz.
        """
        orto = self.p.ortotropo(0.0)
        iso = isotropo_I_L(self.p)
        self.assertLess(m1_error_modulo(iso, orto), 1e-14)
        self.assertGreater(m2_error_constitutiva(iso, orto), 0.1)

    def test_M2_es_simetrica_en_su_definicion(self):
        """La norma de Frobenius de la diferencia no depende del orden salvo
        por el denominador; se comprueba que el valor es finito y positivo."""
        v = m2_error_constitutiva(isotropo_I_L(self.p), self.p.ortotropo(30.0))
        self.assertGreater(v, 0.0)
        self.assertTrue(np.isfinite(v))

    def test_M3_nula_en_ejes_principales(self):
        for theta in (0.0, 90.0):
            with self.subTest(theta=theta):
                self.assertLess(m3_acoplamiento_omitido(
                    self.p.ortotropo(theta)), 1e-14)

    def test_M3_no_nula_fuera_de_ejes(self):
        for theta in (15.0, 30.0, 45.0, 60.0, 75.0):
            with self.subTest(theta=theta):
                self.assertGreater(m3_acoplamiento_omitido(
                    self.p.ortotropo(theta)), 1e-3)

    def test_M3_nula_para_material_isotropo(self):
        """Un isótropo no acopla, gire lo que gire: confirma que M3 mide el
        fenómeno y no un artefacto de la rotación."""
        p_iso = ParametrosBambu(anisotropia=1.0, nu_LT=0.3,
                                g_normalizado=1.0 / (2.0 * (1.0 + 0.3)))
        for theta in (0.0, 23.0, 45.0, 71.0):
            with self.subTest(theta=theta):
                self.assertLess(m3_acoplamiento_omitido(
                    p_iso.ortotropo(theta)), 1e-14)

    def test_metricas_invariantes_de_escala(self):
        """Cambiar las unidades de E_L no puede cambiar un error adimensional.

        Blinda que el barrido pueda parametrizarse por relaciones en vez de por
        valores absolutos (§5 del diseño).
        """
        base = ParametrosBambu(E_L=15.0e9)
        escalado = ParametrosBambu(E_L=15.0e3)   # mismas relaciones, otra escala
        for theta in (0.0, 30.0, 60.0):
            with self.subTest(theta=theta):
                e1 = m1_error_modulo(isotropo_I_L(base), base.ortotropo(theta))
                e2 = m1_error_modulo(isotropo_I_L(escalado),
                                     escalado.ortotropo(theta))
                self.assertAlmostEqual(e1, e2, places=10)


class TestMetricasDeCampo(unittest.TestCase):

    def test_M4_nula_para_campos_identicos(self):
        u = np.array([1.0, -2.0, 0.5])
        self.assertAlmostEqual(m4_error_desplazamiento(u, u), 0.0, places=15)

    def test_M4_es_relativa(self):
        """Escalar ambos campos no cambia el error relativo."""
        a, b = np.array([1.0, 2.0]), np.array([1.1, 2.2])
        self.assertAlmostEqual(
            m4_error_desplazamiento(a, b),
            m4_error_desplazamiento(1000 * a, 1000 * b), places=12)

    def test_M5_nula_para_esfuerzos_identicos(self):
        s = np.array([[1.0, 2.0, 0.3], [0.5, -1.0, 0.2]])
        self.assertAlmostEqual(m5_error_esfuerzos(s, s), 0.0, places=15)

    def test_M5_no_depende_del_numero_de_puntos(self):
        """Duplicar los puntos de Gauss con el mismo campo no cambia M5: es lo
        que garantiza que la métrica compare mallas distintas."""
        s_iso = np.array([[1.1, 2.2, 0.33]])
        s_ort = np.array([[1.0, 2.0, 0.30]])
        e_1 = m5_error_esfuerzos(s_iso, s_ort)
        e_2 = m5_error_esfuerzos(np.repeat(s_iso, 4, axis=0),
                                 np.repeat(s_ort, 4, axis=0))
        self.assertAlmostEqual(e_1, e_2, places=12)


class TestCriterioDeAdmisibilidad(unittest.TestCase):
    """Los umbrales del §7.3 están declarados antes de ver resultados."""

    def test_clasificacion(self):
        self.assertEqual(clasificar(0.01), "admisible")
        self.assertEqual(clasificar(0.10), "con reservas")
        self.assertEqual(clasificar(0.50), "inadmisible")

    def test_fronteras(self):
        self.assertEqual(clasificar(UMBRAL_ADMISIBLE - 1e-9), "admisible")
        self.assertEqual(clasificar(UMBRAL_ADMISIBLE), "con reservas")
        self.assertEqual(clasificar(UMBRAL_RESERVAS), "inadmisible")


if __name__ == '__main__':
    unittest.main()
