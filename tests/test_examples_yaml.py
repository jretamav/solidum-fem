"""Tests end-to-end de regresión sobre `examples/*.yaml`.

Carga cada YAML del directorio `examples/`, lo ejecuta vía
`fenix.run_yaml(...)`, y compara una cantidad característica del
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

import fenix
from fenix.results import (
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


class TestExamplesYAMLRegression(unittest.TestCase):
    """Regresión end-to-end de cada YAML en `examples/`."""

    def test_modelo_elastico_2d(self):
        """Estático lineal 2D: `LinearSolver` sobre Quad4 plane stress."""
        r = fenix.run_yaml(_path('modelo_elastico_2d.yaml'))
        self.assertIsInstance(r, SolveResult)
        self.assertAlmostEqual(float(np.max(np.abs(r.U))), 3.000000e-03,
                               delta=1.0e-8)

    def test_modelo_marco(self):
        """Marco 2D no lineal: `NonlinearSolver` sobre Frame2DEulerCorot."""
        r = fenix.run_yaml(_path('modelo_marco.yaml'))
        self.assertIsInstance(r, SolveResult)
        self.assertAlmostEqual(float(np.max(np.abs(r.U))), 1.921236e-01,
                               delta=1.0e-5)

    def test_modelo_plasticidad(self):
        """Plasticidad J2 plane strain sobre Quad4 con paso adaptativo."""
        r = fenix.run_yaml(_path('modelo_plasticidad.yaml'))
        self.assertIsInstance(r, SolveResult)
        self.assertAlmostEqual(float(np.max(np.abs(r.U))), 2.000000e-03,
                               delta=1.0e-8)

    def test_modelo_modal(self):
        """Modal: `ModalSolver` ARPACK shift-invert."""
        r = fenix.run_yaml(_path('modelo_modal.yaml'))
        self.assertIsInstance(r, ModalResult)
        # frecuencia angular del primer modo (rad/s)
        omega1 = float(r.frequencies_rad[0])
        self.assertAlmostEqual(omega1, 7979.703815, delta=1.0e-3)

    def test_modelo_dinamico_plastico(self):
        """Newton-Newmark sobre material plástico, transitorio con desplazamiento prescrito."""
        r = fenix.run_yaml(_path('modelo_dinamico_plastico.yaml'))
        self.assertIsInstance(r, TransientResult)
        self.assertEqual(r.u_history.shape[1], 505)
        self.assertAlmostEqual(float(r.t_history[-1]), 6.3, delta=1.0e-6)
        self.assertAlmostEqual(float(np.max(np.abs(r.u_history))), 0.5,
                               delta=1.0e-8)

    def test_modelo_central_difference(self):
        """Central difference explícito sobre transitorio lineal."""
        r = fenix.run_yaml(_path('modelo_central_difference.yaml'))
        self.assertIsInstance(r, TransientResult)
        self.assertEqual(r.u_history.shape[1], 401)
        self.assertAlmostEqual(float(r.t_history[-1]), 2.5133, delta=1.0e-3)
        self.assertAlmostEqual(float(np.max(np.abs(r.u_history))), 1.0,
                               delta=1.0e-8)

    def test_modelo_harmonic(self):
        """Barrido armónico complejo en frecuencia."""
        r = fenix.run_yaml(_path('modelo_harmonic.yaml'))
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
        r = fenix.run_yaml(_path('modelo_response_spectrum.yaml'))
        self.assertIsInstance(r, ResponseSpectrumResult)
        self.assertEqual(r.u_combined.shape, (8,))
        self.assertAlmostEqual(
            float(np.max(np.abs(r.u_combined))), 0.2749997,
            delta=1.0e-5,
        )

    def test_modelo_placa(self):
        """Placa rectangular con malla cuadrilátera estructurada."""
        r = fenix.run_yaml(_path('modelo_placa.yaml'))
        self.assertIsInstance(r, SolveResult)
        self.assertAlmostEqual(float(np.max(np.abs(r.U))), 5.070000e-04,
                               delta=1.0e-9)

    def test_placa_gmsh(self):
        """Placa con malla importada desde Gmsh (.msh)."""
        r = fenix.run_yaml(_path('placa_gmsh.yaml'))
        self.assertIsInstance(r, SolveResult)
        self.assertAlmostEqual(float(np.max(np.abs(r.U))), 2.000000e-03,
                               delta=1.0e-8)


if __name__ == '__main__':
    unittest.main()
