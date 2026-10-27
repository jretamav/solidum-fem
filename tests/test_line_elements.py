"""Elementos de línea: celda de referencia, campo interior y fuerza de cuerpo
como campo (ADR 0022, fase 4).

Las armaduras y los marcos declaran ``REFERENCE_CELL = "line"`` y exponen su
campo interior en ejes locales (``line_field``) con las rigideces de la
sección (``section_rigidities``). Se contrasta con lo que no puede mentir:

- **energía**: ``∫ eᵀ·D·e dx`` con el campo interior es ``uᵀ·K·u`` con la
  rigidez del elemento, para ``u`` al azar y un elemento inclinado (en 3D,
  orientado al azar). Comprueba a la vez la interpolación, las rigideces y el
  paso a ejes locales; en el ``Frame2DTimoshenko``, que las funciones de
  Friedman y Kosmatka dan su rigidez con la corrección ``Φ``;
- **extremos**: en ``ξ = 0`` y ``ξ = 1`` el campo vale los DOFs locales;
- **carga**: un campo constante da lo mismo que la constante; una carga
  lineal en ejes locales da las fuerzas de empotramiento clásicas
  (``L(7q₁ + 3q₂)/20``, ``L²(3q₁ + 2q₂)/60``...), con los signos de la
  regla de la mano derecha en el plano ``xz`` del ``Frame3D``.
"""
import logging
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401,E402
from solidum.core.node import Node  # noqa: E402
from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.mesh import structured_mesh  # noqa: E402
from solidum.registry import ElementRegistry  # noqa: E402

MAT = Elastic1D(210.0)
X2 = ([0.3, 0.2], [1.5, 1.1])
X3 = ([0.3, 0.2, -0.4], [1.5, 1.1, 0.6])
KWARGS = {
    "Truss2D": dict(A=2.0), "Truss2DCorot": dict(A=2.0), "Truss3D": dict(A=2.0),
    "Truss3DCorot": dict(A=2.0),
    "Frame2DEuler": dict(A=2.0, I=0.7), "Frame2DEulerCorot": dict(A=2.0, I=0.7),
    "Frame2DTimoshenko": dict(A=2.0, I=0.7, As=1.6, nu=0.25),
    "Frame3D": dict(A=2.0, Iy=0.7, Iz=0.4, J=0.9, nu=0.25),
}
DOFS_DE_LINEA = {"ux", "uy", "uz", "rx", "ry", "rz"}


def cerca(actual, esperado, rtol=1e-12):
    """``assert_allclose`` con tolerancia absoluta escalada (puede haber
    componentes nulas, cuyo redondeo cambia con la plataforma)."""
    esperado = np.asarray(esperado, dtype=float)
    np.testing.assert_allclose(actual, esperado, rtol=rtol,
                               atol=rtol * max(np.abs(esperado).max(), 1e-300))


def de_linea():
    """Elementos de línea del programa principal: dos nodos y sólo DOFs de
    desplazamiento y giro (los módulos de usuario pueden estar cargados)."""
    out = []
    for n in ElementRegistry.names():
        cls = ElementRegistry.get(n)
        if (cls.__module__.startswith("solidum.elements.")
                and set(cls.DOF_NAMES) <= DOFS_DE_LINEA and "ux" in cls.DOF_NAMES
                and getattr(cls, "STRAIN_DIM", None) == 1):
            out.append(cls)
    return out


def lineales():
    return [c for c in de_linea() if not c.GEOMETRICALLY_NONLINEAR
            and c.__name__ in KWARGS]


def elemento(cls):
    X = X3 if "uz" in cls.DOF_NAMES else X2
    return cls(1, [Node(1, X[0]), Node(2, X[1])], MAT, **KWARGS[cls.__name__])


def setUpModule():
    logging.disable(logging.WARNING)


def tearDownModule():
    logging.disable(logging.NOTSET)


class TestContrato(unittest.TestCase):

    def test_todo_elemento_de_linea_declara_su_celda(self):
        """Ningún elemento de línea queda fuera del generador ni de la
        verificación por soluciones manufacturadas (salvo los cables, que no
        tienen campo interior)."""
        nombres = {c.__name__ for c in de_linea()}
        self.assertGreaterEqual(nombres, {"Truss2D", "Truss3D", "Frame2DEuler",
                                          "Frame2DTimoshenko", "Frame3D"})
        for cls in de_linea():
            if cls.__name__.startswith("Cable"):
                continue
            with self.subTest(elemento=cls.__name__):
                self.assertEqual(cls.REFERENCE_CELL, "line")
                self.assertTrue(cls.LINE_FIELDS)
                self.assertTrue(cls.LINE_STRAINS)
                self.assertIn(cls.__name__, KWARGS, "falta en KWARGS de esta prueba")

    def test_energia_del_campo_interior_es_la_de_la_rigidez(self):
        xg, wg = np.polynomial.legendre.leggauss(6)
        xg, wg = (xg + 1) / 2, wg / 2
        rng = np.random.default_rng(7)
        for cls in lineales():
            with self.subTest(elemento=cls.__name__):
                el = elemento(cls)
                n = 2 * len(cls.DOF_NAMES)
                K, _ = el.compute_element_state(np.zeros(n))
                D = el.section_rigidities()
                self.assertEqual(len(D), len(cls.LINE_STRAINS))
                for _ in range(3):
                    u = rng.standard_normal(n)
                    _, e = el.line_field(u, xg)
                    energia = el.L0 * (wg @ np.einsum("gi,i,gi->g", e, D, e))
                    cerca(energia, u @ K @ u)

    def test_extremos_valen_los_dofs_locales(self):
        rng = np.random.default_rng(8)
        for cls in lineales():
            with self.subTest(elemento=cls.__name__):
                el = elemento(cls)
                u = rng.standard_normal(2 * len(cls.DOF_NAMES))
                campo, _ = el.line_field(u, [0.0, 1.0])
                cerca(campo.reshape(-1), el._to_local(u))


class TestCargaComoCampo(unittest.TestCase):

    def test_campo_constante_igual_que_la_constante(self):
        for cls in de_linea():
            if cls.__name__.startswith("Cable"):
                continue
            with self.subTest(elemento=cls.__name__):
                el = elemento(cls)
                b = np.array([0.3, -1.2, 0.5])[:2 + ("uz" in cls.DOF_NAMES)]
                cerca(el.compute_body_load(lambda x: b), el.compute_body_load(b))

    def test_armadura_carga_lineal(self):
        """Cada componente global con funciones lineales:
        ``f_i = A·L·(2b_i + b_j)/6``."""
        for cls in (ElementRegistry.get("Truss2D"), ElementRegistry.get("Truss3D")):
            with self.subTest(elemento=cls.__name__):
                el = elemento(cls)
                X = el._line_ends()
                g = np.array([0.4, -0.7, 1.1])[:X.shape[1]]
                b = lambda x: 2.0 + g * (x @ g)
                bi, bj = b(X[0]), b(X[1])
                L, A = el.L0, el.A
                cerca(el.compute_body_load(b),
                      np.concatenate([A * L * (2 * bi + bj) / 6, A * L * (bi + 2 * bj) / 6]))

    def empotramiento(self, L, q1, q2):
        """Fuerzas y momentos de empotramiento de una carga transversal lineal
        de ``q1`` a ``q2`` (Hermite)."""
        return (L * (7 * q1 + 3 * q2) / 20, L**2 * (3 * q1 + 2 * q2) / 60,
                L * (3 * q1 + 7 * q2) / 20, -L**2 * (2 * q1 + 3 * q2) / 60)

    def carga_local_lineal(self, el, c1, c2):
        """Campo ``b(x)`` cuya carga por unidad de longitud en ejes locales
        varía linealmente de ``c1`` a ``c2`` (por componente local)."""
        X = el._line_ends()
        R = el._load_axes()
        e = (X[1] - X[0]) / el.L0

        def b(x):
            s = (x - X[0]) @ e / el.L0
            return R.T @ ((1 - s) * np.asarray(c1) + s * np.asarray(c2)) / el._axial_area()
        return b

    def test_marco_2d_carga_lineal(self):
        el = elemento(ElementRegistry.get("Frame2DEuler"))
        p1, p2, q1, q2 = 1.5, -0.5, 2.0, -3.0
        f = el.T @ el.compute_body_load(self.carga_local_lineal(el, [p1, q1], [p2, q2]))
        F1, M1, F2, M2 = self.empotramiento(el.L0, q1, q2)
        L = el.L0
        cerca(f, [L * (2 * p1 + p2) / 6, F1, M1, L * (p1 + 2 * p2) / 6, F2, M2])

    def test_marco_3d_carga_lineal(self):
        """Plano ``xy`` como en 2D; plano ``xz`` con los momentos ``M_y`` de
        signo opuesto (``w' = −θ_y``)."""
        el = elemento(ElementRegistry.get("Frame3D"))
        c1, c2 = np.array([1.5, 2.0, -1.0]), np.array([-0.5, -3.0, 4.0])
        f = el.T @ el.compute_body_load(self.carga_local_lineal(el, c1, c2))
        L = el.L0
        Fy1, Mz1, Fy2, Mz2 = self.empotramiento(L, c1[1], c2[1])
        Fz1, My1, Fz2, My2 = self.empotramiento(L, c1[2], c2[2])
        cerca(f, [L * (2 * c1[0] + c2[0]) / 6, Fy1, Fz1, 0.0, -My1, Mz1,
                  L * (c1[0] + 2 * c2[0]) / 6, Fy2, Fz2, 0.0, -My2, Mz2])

    def test_timoshenko_carga_uniforme_independiente_de_phi(self):
        """Con carga uniforme las funciones de Friedman y Kosmatka dan
        ``qL/2`` y ``±qL²/12`` para cualquier ``Φ``, como Hermite."""
        cls = ElementRegistry.get("Frame2DTimoshenko")
        for As in (1e-3, 1.6, 1e6):
            with self.subTest(As=As):
                el = cls(1, [Node(1, X2[0]), Node(2, X2[1])], MAT, A=2.0, I=0.7, As=As, nu=0.25)
                b = np.array([0.3, -1.2])
                cerca(el.compute_body_load(lambda x: b), el.compute_body_load(b))

    def test_timoshenko_esbelto_reparte_como_euler(self):
        """Con ``Φ → 0`` (``A_s`` enorme) la carga lineal se reparte como en
        Euler-Bernoulli."""
        cls = ElementRegistry.get("Frame2DTimoshenko")
        el = cls(1, [Node(1, X2[0]), Node(2, X2[1])], MAT, A=2.0, I=0.7, As=1e12, nu=0.25)
        f = el.T @ el.compute_body_load(self.carga_local_lineal(el, [0.0, 2.0], [0.0, -3.0]))
        F1, M1, F2, M2 = self.empotramiento(el.L0, 2.0, -3.0)
        cerca(f, [0.0, F1, M1, 0.0, F2, M2], rtol=1e-9)

    def test_cables_rechazan_el_campo(self):
        for nombre in ("Cable2DCorot", "Cable3DCorot"):
            cls = ElementRegistry.get(nombre)
            X = X3 if "uz" in cls.DOF_NAMES else X2
            el = cls(1, [Node(1, X[0]), Node(2, X[1])], MAT, A=2.0)
            with self.subTest(elemento=nombre), self.assertRaisesRegex(ValueError, "campo"):
                el.compute_body_load(lambda x: np.zeros(len(X[0])))


class TestMallaDeLinea(unittest.TestCase):

    def test_linea_inclinada_en_3d(self):
        X0, X1 = np.array([0.5, -0.3, 0.2]), np.array([1.7, 0.9, 1.0])
        m = structured_mesh("Frame3D", (4,), mapping=lambda X: X0 + X[0] * (X1 - X0),
                            material=MAT, **KWARGS["Frame3D"])
        self.assertEqual((len(m.domain.nodes), len(m.domain.elements)), (5, 4))
        self.assertEqual(m.sides, ("xmin", "xmax"))
        cerca(m.nodes_on("xmax")[0].coordinates, X1)
        cerca(m.nodes_on("xmin")[0].coordinates, X0)
        (el, k), = m.boundary["xmax"]
        self.assertIs(el.nodes[k], m.nodes_on("xmax")[0])
        for el in m.domain.elements.values():
            self.assertAlmostEqual(el.L0, np.linalg.norm(X1 - X0) / 4, places=14)

    def test_sin_transformacion_en_el_eje_x(self):
        m = structured_mesh("Truss2D", (3,), size=(2.0,), material=MAT, A=1.0)
        cerca(np.array([n.coordinates for n in m.domain.nodes.values()]),
              [[0, 0], [2 / 3, 0], [4 / 3, 0], [2, 0]])


if __name__ == "__main__":
    unittest.main()
