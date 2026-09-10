"""Tests end-to-end de regresión sobre `examples/*.yaml`.

Carga cada YAML del directorio `examples/`, lo ejecuta vía
`solidum.run_yaml(...)`, y compara una cantidad característica del
resultado contra el valor cacheado. Blinda el pipeline completo
(parser YAML → dispatcher por `PIPELINE_KIND` → solver → serialización
del resultado) contra regresiones silenciosas.

Valores de referencia: fijados ejecutando cada ejemplo el 2026-05-19 y
registrando una cantidad estable (típicamente `|u|_∞`, `ω_1`, o el
máximo del envelope). Cualquier desviación >0.1 % en un cambio futuro
indica regresión real — no perfeccionar la tolerancia hasta que aparezca
una causa específica que lo justifique.

Cada ejemplo está en su propio test method para localizar el fallo en
caso de regresión sin necesidad de ejecutar todos.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum
from solidum.results import (
    HarmonicResult,
    ModalResult,
    ResponseSpectrumResult,
    SolveResult,
    TransientResult,
)


_EXAMPLES_DIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'examples',
))


def _path(name: str) -> str:
    return os.path.join(_EXAMPLES_DIR, name)


def _malla_faltante(ruta_yaml: str) -> bool:
    """True si el YAML declara una malla `.msh` que no está en disco.

    Los `.geo` de gmsh sí se versionan; los `.msh` que producen, no todos.
    Distinguir "falta una dependencia externa declarada" de "el ejemplo está
    roto" evita que el test se vuelva ruido y acabe silenciado.
    """
    with open(ruta_yaml, encoding='utf-8') as f:
        for linea in f:
            limpia = linea.split('#')[0].strip()
            if limpia.startswith('mesh:'):
                nombre = limpia.split(":", 1)[1].strip().strip("\"'")
                if nombre.endswith('.msh'):
                    destino = os.path.join(
                        os.path.dirname(ruta_yaml), nombre)
                    return not os.path.isfile(destino)
    return False


class TestExamplesYAMLRegression(unittest.TestCase):
    """Regresión end-to-end de cada YAML en `examples/`."""

    def test_modelo_elastico_2d(self):
        """Estático lineal 2D: `LinearSolver` sobre Quad4 plane stress."""
        r = solidum.run_yaml(_path('modelo_elastico_2d.yaml'))
        self.assertIsInstance(r, SolveResult)
        self.assertAlmostEqual(float(np.max(np.abs(r.U))), 3.000000e-03,
                               delta=1.0e-8)

    def test_modelo_marco(self):
        """Marco 2D no lineal: `NonlinearSolver` sobre Frame2DEulerCorot."""
        r = solidum.run_yaml(_path('modelo_marco.yaml'))
        self.assertIsInstance(r, SolveResult)
        self.assertAlmostEqual(float(np.max(np.abs(r.U))), 1.921236e-01,
                               delta=1.0e-5)

    def test_modelo_plasticidad(self):
        """Plasticidad J2 plane strain sobre Quad4 con paso adaptativo."""
        r = solidum.run_yaml(_path('modelo_plasticidad.yaml'))
        self.assertIsInstance(r, SolveResult)
        self.assertAlmostEqual(float(np.max(np.abs(r.U))), 2.000000e-03,
                               delta=1.0e-8)

    def test_modelo_modal(self):
        """Modal: `ModalSolver` ARPACK shift-invert."""
        r = solidum.run_yaml(_path('modelo_modal.yaml'))
        self.assertIsInstance(r, ModalResult)
        # frecuencia angular del primer modo (rad/s)
        omega1 = float(r.frequencies_rad[0])
        self.assertAlmostEqual(omega1, 7979.703815, delta=1.0e-3)

    def test_modelo_dinamico_plastico(self):
        """Newton-Newmark sobre material plástico, transitorio con desplazamiento prescrito."""
        r = solidum.run_yaml(_path('modelo_dinamico_plastico.yaml'))
        self.assertIsInstance(r, TransientResult)
        self.assertEqual(r.u_history.shape[1], 505)
        self.assertAlmostEqual(float(r.t_history[-1]), 6.3, delta=1.0e-6)
        self.assertAlmostEqual(float(np.max(np.abs(r.u_history))), 0.5,
                               delta=1.0e-8)

    def test_modelo_central_difference(self):
        """Central difference explícito sobre transitorio lineal."""
        r = solidum.run_yaml(_path('modelo_central_difference.yaml'))
        self.assertIsInstance(r, TransientResult)
        self.assertEqual(r.u_history.shape[1], 401)
        self.assertAlmostEqual(float(r.t_history[-1]), 2.5133, delta=1.0e-3)
        self.assertAlmostEqual(float(np.max(np.abs(r.u_history))), 1.0,
                               delta=1.0e-8)

    def test_modelo_harmonic(self):
        """Barrido armónico complejo en frecuencia."""
        r = solidum.run_yaml(_path('modelo_harmonic.yaml'))
        self.assertIsInstance(r, HarmonicResult)
        self.assertEqual(len(r.omega), 91)
        self.assertAlmostEqual(float(r.omega[0]), 1.0, delta=1.0e-12)
        self.assertAlmostEqual(float(r.omega[-1]), 10.0, delta=1.0e-12)
        amp = r.amplitude()
        self.assertEqual(amp.shape, (4, 91))
        self.assertAlmostEqual(float(np.max(np.abs(amp))), 0.4,
                               delta=1.0e-8)

    def test_modelo_response_spectrum(self):
        """Combinación modal espectral SRSS/CQC."""
        r = solidum.run_yaml(_path('modelo_response_spectrum.yaml'))
        self.assertIsInstance(r, ResponseSpectrumResult)
        self.assertEqual(r.u_combined.shape, (8,))
        self.assertAlmostEqual(
            float(np.max(np.abs(r.u_combined))), 0.2749997,
            delta=1.0e-5,
        )

    def test_modelo_placa(self):
        """Placa rectangular con malla cuadrilátera estructurada."""
        r = solidum.run_yaml(_path('modelo_placa.yaml'))
        self.assertIsInstance(r, SolveResult)
        self.assertAlmostEqual(float(np.max(np.abs(r.U))), 5.070000e-04,
                               delta=1.0e-9)

    def test_placa_gmsh(self):
        """Placa con malla importada desde Gmsh (.msh)."""
        r = solidum.run_yaml(_path('placa_gmsh.yaml'))
        self.assertIsInstance(r, SolveResult)
        self.assertAlmostEqual(float(np.max(np.abs(r.U))), 2.000000e-03,
                               delta=1.0e-8)


class TestTodosLosEjemplosEjecutan(unittest.TestCase):
    """Descubrimiento **recursivo**: ningún YAML de `examples/` queda fuera.

    Los tests de arriba fijan cifras concretas para los ejemplos históricos,
    pero los enumeran a mano: un YAML nuevo —o uno colocado en un
    subdirectorio, como pide la convención de `examples/README.md`— no
    entraría en la regresión y podría romperse en silencio.

    Este test no verifica valores; verifica que **todo** ejemplo del
    directorio sigue ejecutándose de extremo a extremo. Es el invariante que
    hace segura la organización en carpetas.
    """

    def test_todos_los_yaml_se_ejecutan(self):
        yamls = sorted(
            os.path.join(raiz, f)
            for raiz, _, archivos in os.walk(_EXAMPLES_DIR)
            for f in archivos
            if f.endswith('.yaml')
        )
        self.assertGreater(len(yamls), 0, "No se encontró ningún YAML")

        for ruta in yamls:
            rel = os.path.relpath(ruta, _EXAMPLES_DIR)
            with self.subTest(ejemplo=rel):
                if _malla_faltante(ruta):
                    # Caso legítimo: los ejemplos de la charla consumen un
                    # `.msh` que NO está versionado —sólo el `.geo` que lo
                    # genera—, así que requieren correr gmsh antes. No es un
                    # ejemplo roto; es una dependencia externa declarada.
                    self.skipTest(
                        f"{rel} necesita una malla .msh no versionada "
                        f"(generar con gmsh desde el .geo correspondiente)"
                    )
                resultado = solidum.run_yaml(ruta)
                self.assertIsInstance(
                    resultado,
                    (SolveResult, ModalResult, TransientResult,
                     HarmonicResult, ResponseSpectrumResult),
                    f"{rel} devolvió un tipo de resultado inesperado",
                )

    def test_los_ejemplos_compuestos_tienen_readme(self):
        """La convención de `examples/README.md`: un ejemplo en carpeta propia
        se documenta con su propio README, o deja de ser navegable."""
        for entrada in sorted(os.listdir(_EXAMPLES_DIR)):
            ruta = os.path.join(_EXAMPLES_DIR, entrada)
            if not os.path.isdir(ruta) or entrada.startswith(('_', '.')):
                continue
            with self.subTest(carpeta=entrada):
                self.assertTrue(
                    os.path.isfile(os.path.join(ruta, 'README.md')),
                    f"examples/{entrada}/ no tiene README.md "
                    f"(ver la convención en examples/README.md)",
                )


if __name__ == '__main__':
    unittest.main()
