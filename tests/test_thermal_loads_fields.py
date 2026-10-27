"""Fuente y flujo de frontera variables en los elementos térmicos (ADR 0022, fase 3).

``compute_body_source`` acepta un campo ``Q(x)``; ``compute_edge_flux`` (2D)
y ``compute_face_flux`` (3D), un campo ``q̄(x, n)`` con ``n`` la normal
exterior del borde o de la cara discretos. Comprobaciones analíticas sobre
elementos de nodos perturbados al azar:

- fuente lineal: ``Σ f_i = ∫Q dΩ`` y ``Σ f_i·x_i = ∫Q·x dΩ``
  (``Σ N_i x_i = x``), contra una cuadratura de Gauss de 8 puntos por
  dirección sobre la misma geometría, en lo que la regla de la capacidad
  integra exactamente; la fuente no depende de la regla de la conductividad
  (1×1 o 2×2);
- flujo ``q̄ = g·n`` con ``g`` constante: en un borde, ``−(g·n)·L·t/2`` en cada
  nodo con ``n`` la normal exterior; en una cara alabeada, resultante
  ``−g·½(d₁ × d₂)``, con ``d₁``, ``d₂`` las diagonales (vector área de un
  cuadrilátero de lados rectos);
- flujo lineal sobre un borde: ``−L·t·(2q_a + q_c)/6`` y ``−L·t·(q_a + 2q_c)/6``.

Signo (Reglas §5): ``q̄ > 0`` sale del dominio y el vector nodal es negativo.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401,E402
from solidum.core.node import Node  # noqa: E402
from solidum.elements.thermal import Hex8Thermal, Quad4Thermal  # noqa: E402
from solidum.materials.thermal_conduction import ThermalConduction  # noqa: E402

MAT = {2: ThermalConduction([[3.0, 1.0], [1.0, 2.0]]),
       3: ThermalConduction([[3.0, 0.5, 0.2], [0.5, 2.0, 0.4], [0.2, 0.4, 1.5]])}
ESPESOR = 0.3


def cerca(actual, esperado, rtol=1e-13):
    """``assert_allclose`` con tolerancia absoluta escalada (puede haber
    componentes nulas, cuyo redondeo cambia con la plataforma)."""
    esperado = np.asarray(esperado, dtype=float)
    np.testing.assert_allclose(actual, esperado, rtol=rtol,
                               atol=rtol * max(np.abs(esperado).max(), 1e-300))


def elemento(cls, semilla=0, perturbacion=0.15, **kw):
    """Elemento suelto: sus nodos de referencia llevados a ``[0, 1]^d``,
    deformados por una matriz y movidos al azar hasta ``perturbacion``."""
    xi = (np.asarray(cls.REFERENCE_NODES, dtype=float) + 1.0) / 2.0
    d = xi.shape[1]
    M = np.array([[2.0, 0.3, 0.1], [0.4, 1.5, -0.2], [0.0, 0.3, 1.2]])[:d, :d]
    X = xi @ M + np.random.default_rng(semilla).uniform(-perturbacion, perturbacion, xi.shape)
    nodos = [Node(i + 1, list(x)) for i, x in enumerate(X)]
    if d == 2:
        return cls(1, nodos, MAT[2], thickness=ESPESOR, **kw)
    return cls(1, nodos, MAT[3], **kw)


def integral(el, f):
    """``∫ f(x) dΩ`` (por el espesor en 2D) con Gauss de 8 puntos por
    dirección sobre la geometría isoparamétrica del elemento."""
    d = el.FLUX_DIM
    a, w = np.polynomial.legendre.leggauss(8)
    malla = np.meshgrid(*([a] * d), indexing="ij")
    pts = np.column_stack([m.ravel() for m in malla])
    W = np.prod(np.meshgrid(*([w] * d), indexing="ij"), axis=0).ravel()
    coords = el.get_coordinate_matrix(ndim=d)
    total = 0.0
    for p, wp in zip(pts, W):
        x = el._shape_functions(p) @ coords
        total = total + f(x) * el._gradient(p, coords)[1] * wp * el.thickness
    return total


def Q_lineal(x):
    return 5.0 + 2.0 * x[0] - 3.0 * x[1] + (1.5 * x[2] if len(x) > 2 else 0.0)


class TestFuenteVariable(unittest.TestCase):

    def test_fuente_lineal_exacta(self):
        """En el hexaedro distorsionado ``detJ`` es de grado 2 por dirección:
        2×2×2 integra ``Q·detJ`` (grado 3) pero no ``Q·x·detJ`` (grado 4), así
        que el momento se comprueba sobre un paralelepípedo."""
        for cls in (Quad4Thermal, Hex8Thermal):
            for semilla in range(3):
                with self.subTest(elemento=cls.__name__, semilla=semilla):
                    el = elemento(cls, semilla)
                    cerca(el.compute_body_source(Q_lineal).sum(), integral(el, Q_lineal))
                    if cls is Hex8Thermal:
                        el = elemento(cls, semilla, perturbacion=0.0)
                    X = el.get_coordinate_matrix(ndim=el.FLUX_DIM)
                    cerca(el.compute_body_source(Q_lineal) @ X,
                          integral(el, lambda x: Q_lineal(x) * x))

    def test_constante_y_campo_coinciden(self):
        for cls in (Quad4Thermal, Hex8Thermal):
            with self.subTest(elemento=cls.__name__):
                el = elemento(cls)
                cerca(el.compute_body_source(lambda x: 7.0), el.compute_body_source(7.0))
                cerca(el.compute_body_source(7.0).sum(), integral(el, lambda x: 7.0))

    def test_no_depende_de_la_regla_de_la_conductividad(self):
        for cls, reducida in ((Quad4Thermal, "1x1"), (Hex8Thermal, "hex_1x1x1")):
            with self.subTest(elemento=cls.__name__):
                cerca(elemento(cls, quadrature=reducida).compute_body_source(Q_lineal),
                      elemento(cls).compute_body_source(Q_lineal))

    def test_campo_con_componentes_de_mas(self):
        with self.assertRaises(ValueError):
            elemento(Quad4Thermal).compute_body_source(lambda x: [1.0, 2.0])


class TestFlujoDeBorde(unittest.TestCase):

    def test_flujo_normal_con_la_normal_exterior(self):
        el = elemento(Quad4Thermal)
        X = el.get_coordinate_matrix(ndim=2)
        centro = X.mean(axis=0)
        g = np.array([4.0, -1.5])
        for e, (a, c) in enumerate(el.EDGE_NODES):
            with self.subTest(borde=e):
                vistas = []
                f = el.compute_edge_flux(e, lambda x, n: vistas.append(n) or g @ n)
                d = X[c] - X[a]
                L = np.linalg.norm(d)
                n = np.array([d[1], -d[0]]) / L
                self.assertGreater(n @ ((X[a] + X[c]) / 2 - centro), 0.0)
                cerca(np.array(vistas), np.tile(n, (len(vistas), 1)))
                esperado = np.zeros(4)
                esperado[[a, c]] = -(g @ n) * L * ESPESOR / 2
                cerca(f, esperado)

    def test_flujo_lineal_a_lo_largo_del_borde(self):
        el = elemento(Quad4Thermal, semilla=1)
        X = el.get_coordinate_matrix(ndim=2)
        for e, (a, c) in enumerate(el.EDGE_NODES):
            with self.subTest(borde=e):
                f = el.compute_edge_flux(e, lambda x, n: Q_lineal(x))
                L = np.linalg.norm(X[c] - X[a])
                qa, qc = Q_lineal(X[a]), Q_lineal(X[c])
                esperado = np.zeros(4)
                esperado[a] = -L * ESPESOR * (2 * qa + qc) / 6
                esperado[c] = -L * ESPESOR * (qa + 2 * qc) / 6
                cerca(f, esperado)

    def test_constante_sin_cambio(self):
        el = elemento(Quad4Thermal)
        for e in range(4):
            cerca(el.compute_edge_flux(e, 2.5), el.compute_edge_flux(e, lambda x, n: 2.5))


class TestFlujoDeCara(unittest.TestCase):

    def test_resultante_con_el_vector_area(self):
        """Caras alabeadas (nodos perturbados): ``Σ f = −g·½(d₁ × d₂)``."""
        g = np.array([4.0, -1.5, 2.0])
        for semilla in range(3):
            el = elemento(Hex8Thermal, semilla)
            X = el.get_coordinate_matrix(ndim=3)
            for k, cara in enumerate(el.FACE_NODES):
                with self.subTest(semilla=semilla, cara=k):
                    P = X[list(cara)]
                    area = 0.5 * np.cross(P[2] - P[0], P[3] - P[1])
                    f = el.compute_face_flux(k, lambda x, n: g @ n)
                    cerca(f.sum(), -g @ area)
                    fuera = [i for i in range(8) if i not in cara]
                    self.assertTrue(np.all(f[fuera] == 0.0))

    def test_normal_exterior(self):
        el = elemento(Hex8Thermal)
        X = el.get_coordinate_matrix(ndim=3)
        centro = X.mean(axis=0)
        for k, cara in enumerate(el.FACE_NODES):
            with self.subTest(cara=k):
                vistas = []
                el.compute_face_flux(k, lambda x, n: vistas.append((x, n)) or 0.0)
                for x, n in vistas:
                    self.assertAlmostEqual(np.linalg.norm(n), 1.0, places=14)
                    self.assertGreater(n @ (X[list(cara)].mean(axis=0) - centro), 0.0)

    def test_constante_sin_cambio(self):
        el = elemento(Hex8Thermal, semilla=2)
        for k in range(6):
            cerca(el.compute_face_flux(k, 2.5), el.compute_face_flux(k, lambda x, n: 2.5))


if __name__ == "__main__":
    unittest.main()
