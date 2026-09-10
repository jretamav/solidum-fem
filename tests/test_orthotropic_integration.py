"""Integración de `Orthotropic2D` con el resto del sistema.

Complementa los otros dos frentes ortótropos:

- `tests/test_orthotropic_2d.py` — verificación de la constitutiva.
- `tests/validation/test_off_axis_orthotropic.py` — validación externa contra
  Jones (1999) §2.8, incluidos modelos FEM con los cinco elementos 2D.
- `tests/test_patch_orthotropic_2d.py` — patch test con constitutiva llena.

Lo que queda fuera de esos tres y cubre este archivo: los caminos por los que
el material se usa **desde fuera** —YAML, matriz de masa, contrato genérico—
y que ninguno de ellos ejercita.

Por qué importa cada bloque
---------------------------
1. **Regresión del ejemplo publicado.** `examples/bambu_ortotropo/README.md`
   publica una tabla de siete filas con desplazamientos y módulos aparentes.
   `tests/test_examples_yaml.py` ejecuta ese YAML pero sólo comprueba que no
   lance excepción — los valores no están anclados. Sin este bloque, un cambio
   en el material dejaría la documentación mintiendo en silencio.

   *Alcance honesto*: las cifras se generaron con este mismo código, así que
   esto es **regresión**, no validación — detecta cambios respecto al estado
   validado, no errores que ya estuvieran presentes. Quien valida contra una
   referencia externa es el benchmark de Jones.

2. **Densidad y masa.** El material acepta `density` y el YAML de bambú la
   declara (700 kg/m³), pero ninguna matriz de masa se construye nunca con él.
   La masa no depende de la ortotropía —es un escalar por el volumen— y
   precisamente por eso conviene comprobar que la anisotropía de `C` no se
   filtra a un sitio donde no pinta nada.

3. **Contrato genérico.** ADR 0013 §2 declara `theta`-en-el-material un atajo
   desechable y promete que migrarlo al elemento será **aditivo**: la firma
   pasaría a `compute_state(strain, state_vars=None, orientation=None)` sin que
   ningún test existente cambie. Ese bloque fija hoy lo que esa migración debe
   preservar, para que la promesa del ADR sea verificable y no sólo declarada.
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum
from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.elements.solid_2d import Quad4
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.orthotropic_2d import Orthotropic2D
from solidum.registry import MaterialRegistry
from solidum.results import SolveResult

_EJEMPLO = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    '..', 'examples', 'bambu_ortotropo', 'tira_45_grados.yaml')

# Material del ejemplo publicado.
BAMBU_EJEMPLO = dict(E1=15.0e9, E2=0.8e9, G12=0.7e9, nu12=0.35)


# =============================================================================
# 1. Regresión del ejemplo publicado
# =============================================================================

class TestRegresionEjemploBambu(unittest.TestCase):
    """Ancla las cifras que `examples/bambu_ortotropo/README.md` publica.

    La tabla del README es material citable —se apoya en ella la lectura
    física del acoplamiento— así que debe romperse ruidosamente si el material
    cambia de comportamiento.
    """

    # (theta, ux, uy) del nodo 25 (esquina superior derecha), en metros.
    # Reproducidos con `barrido_orientacion.py` sobre el estado validado.
    TABLA = (
        (0.0,   0.000667, -0.000233),
        (15.0,  0.001500, -0.003374),
        (30.0,  0.003747, -0.005621),
        (45.0,  0.006746, -0.006313),
        (60.0,  0.009664, -0.005338),
        (75.0,  0.011748, -0.003091),
        (90.0,  0.012500, -0.000233),
    )

    SIGMA = 10.0e6      # tracción aplicada [Pa]
    L = 1.0             # longitud de la tira [m]

    @staticmethod
    def _correr(theta_deg):
        """Ejecuta el YAML del ejemplo con la orientación pedida.

        Sustituye `theta` en memoria como hace `barrido_orientacion.py`, sin
        tocar el archivo en disco.
        """
        import re
        import tempfile

        with open(_EJEMPLO, encoding='utf-8') as f:
            texto = f.read()
        texto = re.sub(r'theta: [0-9.]+', f'theta: {theta_deg}', texto)

        tmp = tempfile.NamedTemporaryFile(
            mode='w', suffix='.yaml', delete=False, encoding='utf-8')
        try:
            tmp.write(texto)
            tmp.close()
            resultado = solidum.run_yaml(tmp.name)
        finally:
            os.unlink(tmp.name)

        return resultado

    def test_el_yaml_del_ejemplo_corre(self):
        r = self._correr(45.0)
        self.assertIsInstance(r, SolveResult)

    def test_tabla_de_desplazamientos_del_readme(self):
        """Los ux, uy publicados, a la precisión con que están tabulados."""
        for theta, ux_ref, uy_ref in self.TABLA:
            with self.subTest(theta=theta):
                U = self._correr(theta).U
                ux, uy = U[-2], U[-1]
                # El README tabula en mm con 3-4 decimales: 1e-6 m de holgura
                # es más fino que la cifra publicada y no es frágil.
                self.assertAlmostEqual(ux, ux_ref, delta=1.0e-6,
                                       msg=f"ux(θ={theta}) fuera de tabla")
                self.assertAlmostEqual(uy, uy_ref, delta=1.0e-6,
                                       msg=f"uy(θ={theta}) fuera de tabla")

    def test_modulo_aparente_cae_por_debajo_del_diez_por_ciento(self):
        """La cifra citable del README y de la spec: a 45° queda < 10% de E1.

        Es contraintuitivo —E1/E2 = 18.75 sugeriría no bajar del 5%— y es
        justamente la conclusión física que el ejemplo existe para mostrar.
        """
        ux = self._correr(45.0).U[-2]
        E_x = self.SIGMA / (ux / self.L)
        self.assertLess(E_x / BAMBU_EJEMPLO['E1'], 0.10)
        # Y no es un colapso numérico: sigue siendo un módulo razonable.
        self.assertGreater(E_x / BAMBU_EJEMPLO['E1'], 0.05)

    def test_poisson_puro_en_ejes_principales(self):
        """En θ=0° y 90° el cociente uy/ux es exactamente −ν12 y −ν21.

        Fuera de eje ese cociente deja de ser Poisson: es acoplamiento. Aquí
        se comprueba el caso limpio, que además verifica la reciprocidad
        propagada hasta el modelo completo.
        """
        nu12 = BAMBU_EJEMPLO['nu12']
        nu21 = nu12 * BAMBU_EJEMPLO['E2'] / BAMBU_EJEMPLO['E1']

        U0 = self._correr(0.0).U
        self.assertAlmostEqual(U0[-1] / U0[-2], -nu12, delta=1.0e-3)

        U90 = self._correr(90.0).U
        self.assertAlmostEqual(U90[-1] / U90[-2], -nu21, delta=1.0e-3)

    def test_fuera_de_eje_el_cociente_excede_cualquier_poisson(self):
        """A 15° el cociente uy/ux supera 2 — imposible para contracción
        transversal, y por tanto prueba directa de acoplamiento."""
        U = self._correr(15.0).U
        self.assertGreater(abs(U[-1] / U[-2]), 2.0)


# =============================================================================
# 2. Densidad y matriz de masa
# =============================================================================

class TestMasaConMaterialOrtotropo(unittest.TestCase):
    """La masa es escalar: la anisotropía de C no debe filtrarse a ella."""

    L, T, RHO = 2.0, 0.05, 700.0

    def _elemento(self, theta):
        nodes = [Node(1, [0.0, 0.0]), Node(2, [self.L, 0.0]),
                 Node(3, [self.L, self.L]), Node(4, [0.0, self.L])]
        mat = Orthotropic2D(**BAMBU_EJEMPLO, theta=theta, density=self.RHO)
        return Quad4(1, nodes, mat, thickness=self.T)

    def test_masa_total_consistente(self):
        """Suma de la masa consistente = ρ·A·t, para cualquier orientación."""
        esperada = self.RHO * self.L * self.L * self.T
        for theta in (0.0, 30.0, 45.0, 90.0):
            with self.subTest(theta=theta):
                M = self._elemento(theta).compute_mass_matrix('consistent')
                # Cada DOF traslacional aparece en la suma; en 2D la masa
                # total se recupera sumando una de las dos direcciones.
                total = M[0::2, 0::2].sum()
                self.assertAlmostEqual(total / esperada, 1.0, places=12)

    def test_masa_independiente_de_la_orientacion(self):
        """Rotar la fibra no mueve un gramo: la masa no ve la constitutiva."""
        M_ref = self._elemento(0.0).compute_mass_matrix('consistent')
        for theta in (15.0, 45.0, 90.0, 123.0):
            with self.subTest(theta=theta):
                M = self._elemento(theta).compute_mass_matrix('consistent')
                np.testing.assert_allclose(M, M_ref, rtol=1e-14, atol=0.0)

    def test_masa_lumped_conserva_el_total(self):
        esperada = self.RHO * self.L * self.L * self.T
        M = self._elemento(45.0).compute_mass_matrix('lumped')
        total = np.diag(M)[0::2].sum()
        self.assertAlmostEqual(total / esperada, 1.0, places=12)

    def test_densidad_es_opcional(self):
        """ADR 0008: sólo la exigen peso propio y masa, no la rigidez."""
        mat = Orthotropic2D(**BAMBU_EJEMPLO, theta=30.0)
        self.assertIsNone(mat.density)
        nodes = [Node(1, [0.0, 0.0]), Node(2, [1.0, 0.0]),
                 Node(3, [1.0, 1.0]), Node(4, [0.0, 1.0])]
        K = Quad4(1, nodes, mat, thickness=0.01).compute_global_stiffness()
        self.assertTrue(np.all(np.isfinite(K)))


# =============================================================================
# 3. Contrato genérico del material
# =============================================================================

class TestContratoGenerico(unittest.TestCase):
    """Lo que la migración del ADR 0013 a orientación-en-el-elemento debe
    preservar. Fija hoy el contrato para que la promesa de "migración aditiva"
    sea verificable y no sólo declarada."""

    def setUp(self):
        self.mat = Orthotropic2D(**BAMBU_EJEMPLO, theta=37.0)

    def test_atributos_declarativos(self):
        self.assertEqual(self.mat.STRAIN_DIM, 3)
        self.assertIsNone(self.mat.PRIMARY_STATE_VAR)
        self.assertTrue(self.mat.IS_SYMMETRIC)

    def test_registrado_por_nombre(self):
        self.assertIn('Orthotropic2D', MaterialRegistry.names())
        creado = MaterialRegistry.create('Orthotropic2D', **BAMBU_EJEMPLO,
                                         theta=37.0)
        np.testing.assert_allclose(creado.C, self.mat.C, rtol=1e-15)

    def test_compute_state_devuelve_la_terna(self):
        eps = np.array([1.0e-3, -2.0e-4, 5.0e-4])
        sigma, C, sv = self.mat.compute_state(eps, None)

        self.assertEqual(sigma.shape, (3,))
        self.assertEqual(C.shape, (3, 3))
        self.assertIsNone(sv)
        np.testing.assert_allclose(sigma, self.mat.C @ eps, rtol=1e-15)

    def test_state_vars_pasa_sin_tocarse(self):
        """Material sin historia: lo que entra como estado, sale igual."""
        centinela = object()
        _, _, sv = self.mat.compute_state(np.zeros(3), centinela)
        self.assertIs(sv, centinela)

    def test_tangente_constante_y_simetrica(self):
        """`IS_SYMMETRIC = True` debe ser cierto, no sólo declarado.

        El despachador algebraico (ADR 0003) elige Cholesky apoyándose en esta
        bandera: si C dejara de ser simétrica, la elección sería inválida.
        """
        C_ref = None
        for eps in (np.zeros(3), np.array([1e-3, 0.0, 0.0]),
                    np.array([-4e-3, 7e-4, -2e-3])):
            _, C, _ = self.mat.compute_state(eps, None)
            np.testing.assert_allclose(C, C.T, rtol=1e-14)
            if C_ref is None:
                C_ref = C
            else:
                np.testing.assert_allclose(C, C_ref, rtol=1e-15)

    def test_tangente_definida_positiva(self):
        """Barrido de orientaciones: C debe seguir siendo definida positiva.

        Es la condición de admisibilidad de la spec §12 y lo que garantiza que
        la rigidez global sea factorizable.
        """
        for theta in range(0, 180, 7):
            with self.subTest(theta=theta):
                mat = Orthotropic2D(**BAMBU_EJEMPLO, theta=float(theta))
                autovalores = np.linalg.eigvalsh(mat.C)
                self.assertGreater(autovalores.min(), 0.0)

    def test_theta_es_periodico_en_ciento_ochenta_grados(self):
        """Un eje material no tiene sentido: θ y θ+180° son el mismo material.

        Blinda que la rotación se aplica con funciones del ángulo doble, como
        corresponde a un tensor de segundo orden, y no con una convención que
        distinga sentidos.
        """
        for theta in (0.0, 23.0, 45.0, 71.0):
            with self.subTest(theta=theta):
                a = Orthotropic2D(**BAMBU_EJEMPLO, theta=theta)
                b = Orthotropic2D(**BAMBU_EJEMPLO, theta=theta + 180.0)
                escala = np.abs(a.C).max()
                np.testing.assert_allclose(b.C, a.C, rtol=0.0,
                                           atol=1e-12 * escala)


class TestValidacionDesdeYAML(unittest.TestCase):
    """El validador ortótropo también protege la entrada por YAML.

    ADR 0013 §4 documenta la trampa: heredar el criterio isótropo
    (−1 < ν < 0.5) rechazaría datos de bambú perfectamente válidos. Aquí se
    comprueba por el camino que recorre un usuario real.
    """

    def test_rechaza_nu12_inadmisible(self):
        """|ν12| ≥ sqrt(E1/E2) hace C singular o indefinida."""
        with self.assertRaises(ValueError) as ctx:
            MaterialRegistry.create('Orthotropic2D', E1=15.0e9, E2=0.8e9,
                                    G12=0.7e9, nu12=5.0)
        mensaje = str(ctx.exception)
        self.assertIn('nu12', mensaje)
        self.assertIn('sqrt(E1/E2)', mensaje)

    def test_acepta_nu12_mayor_que_medio(self):
        """El límite isótropo de 0.5 NO aplica: con E1/E2 = 18.75 el techo
        real es 4.33. Un ν12 de 0.6 es legítimo en bambú."""
        mat = MaterialRegistry.create('Orthotropic2D', E1=15.0e9, E2=0.8e9,
                                      G12=0.7e9, nu12=0.6)
        self.assertGreater(np.linalg.eigvalsh(mat.C).min(), 0.0)

    def test_rechaza_modulos_no_positivos(self):
        for parametro in ('E1', 'E2', 'G12'):
            with self.subTest(parametro=parametro):
                kwargs = dict(E1=15.0e9, E2=0.8e9, G12=0.7e9, nu12=0.35)
                kwargs[parametro] = -1.0
                with self.assertRaises(ValueError) as ctx:
                    MaterialRegistry.create('Orthotropic2D', **kwargs)
                self.assertIn(parametro, str(ctx.exception))


class TestDegeneracionIsotropaEndToEnd(unittest.TestCase):
    """Con datos isótropos, un modelo ortótropo debe dar lo mismo que Elastic2D.

    Es el test de regresión más valioso de la entrega según ADR 0013: el
    isótropo es el caso particular verificable del ortótropo, y aquí se
    comprueba sobre un modelo completo, no sólo sobre la matriz C.
    """

    E, NU = 210.0e9, 0.3

    def _rigidez(self, material):
        nodes = [Node(1, [0.0, 0.0]), Node(2, [1.3, 0.1]),
                 Node(3, [1.2, 1.4]), Node(4, [0.1, 1.1])]
        return Quad4(1, nodes, material, thickness=0.02).compute_global_stiffness()

    def test_rigidez_coincide_con_elastic2d(self):
        G = self.E / (2.0 * (1.0 + self.NU))
        orto = Orthotropic2D(E1=self.E, E2=self.E, G12=G, nu12=self.NU)
        iso = Elastic2D(E=self.E, nu=self.NU, hypothesis='plane_stress')

        K_orto = self._rigidez(orto)
        K_iso = self._rigidez(iso)
        escala = np.abs(K_iso).max()
        np.testing.assert_allclose(K_orto, K_iso, rtol=0.0, atol=1e-10 * escala)

    def test_rigidez_isotropa_no_depende_de_theta(self):
        """Un isótropo rotado sigue siendo el mismo isótropo, y su elemento
        también: la rotación no debe dejar residuo en la rigidez."""
        G = self.E / (2.0 * (1.0 + self.NU))
        K_ref = self._rigidez(
            Orthotropic2D(E1=self.E, E2=self.E, G12=G, nu12=self.NU))
        escala = np.abs(K_ref).max()

        for theta in (17.0, 45.0, 90.0, 133.0):
            with self.subTest(theta=theta):
                K = self._rigidez(Orthotropic2D(
                    E1=self.E, E2=self.E, G12=G, nu12=self.NU, theta=theta))
                np.testing.assert_allclose(K, K_ref, rtol=0.0,
                                           atol=1e-9 * escala)


if __name__ == '__main__':
    unittest.main()
