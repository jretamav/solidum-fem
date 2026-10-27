"""Cargas distribuidas variables en el espacio de los sólidos 2D (ADR 0022).

``compute_body_load`` acepta un campo ``b(x)`` y ``compute_edge_traction`` un
campo ``t(x, n)``, con ``n`` la normal exterior del borde. La carga de borde
se integra sobre la geometría real del borde: hasta el 2026-10-27 los bordes
cuadráticos repartían la cuerda con pesos fijos 1-4-1, exactos sólo con el
nodo medio centrado sobre un borde recto (deuda #13). Las comprobaciones son
analíticas:

- borde recto con el nodo medio en ``a_m`` (posición sobre la cuerda de
  longitud ``L``): con ``|dx/ds| = L/2 + s·(L − 2·a_m)``, la tracción
  uniforme ``t`` da ``f = t·(L/6 − (L − 2a_m)/3, 2L/3, L/6 + (L − 2a_m)/3)``;
  con el nodo centrado, el 1-4-1 de siempre;
- presión uniforme sobre cualquier borde: ``∫ n ds`` es la cuerda girada −90°,
  así que la resultante de ``−p·n`` es ``−p·(Δy, −Δx)`` aunque el borde sea
  curvo;
- tracción lineal en un borde de dos nodos: ``f_a = L(2t_a + t_c)/6``,
  ``f_c = L(t_a + 2t_c)/6``;
- fuerza de cuerpo lineal en un elemento afín: resultante ``b(centroide)·A·e``.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401,E402
from solidum.core.node import Node  # noqa: E402
from solidum.materials.elastic_2d import Elastic2D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.mesh import structured_mesh  # noqa: E402
from solidum.registry import ElementRegistry  # noqa: E402

MAT = Elastic2D(200e3, 0.3)
E = 0.2   # espesor
REF = ("Quad4", "Quad8", "Quad9", "Tri3", "Tri6")


def elemento(nombre):
    """Elemento suelto con sus nodos de referencia llevados a un
    paralelogramo."""
    cls = ElementRegistry.get(nombre)
    xi = np.asarray(cls.REFERENCE_NODES, dtype=float)
    if cls.REFERENCE_CELL == "quad":
        xi = (xi + 1.0) / 2.0
    X = xi @ np.array([[2.0, 0.3], [0.4, 1.5]]) + np.array([1.0, -0.5])
    nodos = [Node(i + 1, list(x)) for i, x in enumerate(X)]
    return cls(1, nodos, MAT, thickness=E)


def nodal(f, k):
    return f[2 * k:2 * k + 2]


def cerca(actual, esperado, rtol=1e-13):
    """``assert_allclose`` con tolerancia absoluta escalada a la carga: hay
    componentes exactamente nulas (los vértices del Tri6 con b constante,
    ``∫N_vértice dA = 0``) cuyo redondeo cambia con la plataforma, y con
    tolerancia sólo relativa la prueba falla en Linux (CI, 2026-10-27)."""
    esperado = np.asarray(esperado)
    np.testing.assert_allclose(actual, esperado, rtol=rtol,
                               atol=rtol * np.abs(esperado).max())


class TestFuerzaDeCuerpo(unittest.TestCase):

    def test_campo_constante_igual_a_la_constante(self):
        for nombre in REF:
            with self.subTest(elemento=nombre):
                el = elemento(nombre)
                b = np.array([3.0, -7.0])
                cerca(el.compute_body_load(lambda x: b), el.compute_body_load(b))

    def test_campo_lineal(self):
        """Resultante ``∫ b dV = b(centroide)·A·e`` en un elemento afín."""
        g = np.array([[2.0, -1.0], [0.5, 3.0]])
        campo = lambda x: g @ x + np.array([1.0, 0.0])  # noqa: E731
        for nombre in REF:
            with self.subTest(elemento=nombre):
                el = elemento(nombre)
                X = el.get_coordinate_matrix(ndim=2)
                vert = X[:3] if nombre.startswith("Tri") else X[:4]
                area = 2.0 * 1.5 - 0.3 * 0.4
                if nombre.startswith("Tri"):
                    area /= 2
                c = vert.mean(axis=0)
                f = el.compute_body_load(campo).reshape(-1, 2).sum(axis=0)
                cerca(f, campo(c) * area * E)

    def test_componentes_erroneas(self):
        el = elemento("Quad4")
        with self.assertRaisesRegex(ValueError, "componentes"):
            el.compute_body_load(lambda x: [1.0, 2.0, 3.0])
        with self.assertRaisesRegex(ValueError, "componentes"):
            el.compute_body_load([1.0, 2.0, 3.0])

    def test_ensamblador_acepta_un_campo(self):
        m = structured_mesh("Quad8", (2, 3), mapping=lambda X: X + 0.05 * np.sin(3 * X[::-1]),
                            material=MAT, thickness=E)
        m.domain.generate_equation_numbers(verbose=False)
        a = Assembler(m.domain)
        F = a.assemble_body_load(lambda x: np.array([x[1], -2.0]))
        G = np.zeros_like(F)
        for el in m.domain.elements.values():
            G[el.get_global_dof_indices()] += el.compute_body_load(lambda x: np.array([x[1], -2.0]))
        cerca(F, G)


class TestCargaDeBorde(unittest.TestCase):

    def test_nodo_medio_descentrado(self):
        """Deuda #13: borde recto inferior de longitud ``L`` con el nodo medio
        en ``a_m = 0.3·L``; el 1-4-1 antiguo daba (1/6, 4/6, 1/6)·L."""
        for nombre in ("Quad8", "Quad9", "Tri6"):
            with self.subTest(elemento=nombre):
                a, m, c = ElementRegistry.get(nombre).EDGE_NODES[0]
                el = elemento(nombre)
                X = el.get_coordinate_matrix(ndim=2)
                el.nodes[m].coordinates[:2] = list(X[a] + 0.3 * (X[c] - X[a]))
                L = np.linalg.norm(X[c] - X[a])
                t = np.array([2.0, -1.0])
                f = el.compute_edge_traction(0, t)
                pesos = (L / 6 - (L - 0.6 * L) / 3, 2 * L / 3, L / 6 + (L - 0.6 * L) / 3)
                for k, w in zip((a, m, c), pesos):
                    cerca(nodal(f, k), w * t * E)
                el.nodes[m].coordinates[:2] = list((X[a] + X[c]) / 2)
                f = el.compute_edge_traction(0, t)
                for k, w in zip((a, m, c), (L / 6, 2 * L / 3, L / 6)):
                    cerca(nodal(f, k), w * t * E)

    def test_presion_sobre_borde_curvo(self):
        """Resultante de ``−p·n`` = ``−p·(Δy, −Δx)`` sobre un arco: exacta."""
        p = 5.0
        for nombre in ("Quad8", "Quad9", "Tri6"):
            for e in range(len(ElementRegistry.get(nombre).EDGE_NODES)):
                with self.subTest(elemento=nombre, borde=e):
                    el = elemento(nombre)
                    a, m, c = el.EDGE_NODES[e]
                    X = el.get_coordinate_matrix(ndim=2)
                    d = X[c] - X[a]
                    bombeo = (X[a] + X[c]) / 2 + 0.15 * np.array([d[1], -d[0]])
                    el.nodes[m].coordinates[:2] = list(bombeo)
                    f = el.compute_edge_traction(e, lambda x, n: -p * n)
                    np.testing.assert_allclose(f.reshape(-1, 2).sum(axis=0),
                                               -p * np.array([d[1], -d[0]]) * E, atol=1e-12)

    def test_normal_exterior(self):
        """Con ``t = n`` sobre el borde derecho de un cuadrado unidad, cada nodo
        recibe su parte de ``(1, 0)·L``."""
        for nombre in REF:
            with self.subTest(elemento=nombre):
                m = structured_mesh(nombre, (1, 1), material=MAT, thickness=E)
                for lado, n in (("xmax", (1, 0)), ("xmin", (-1, 0)),
                                ("ymax", (0, 1)), ("ymin", (0, -1))):
                    total = sum(el.compute_edge_traction(e, lambda x, nn: nn).reshape(-1, 2).sum(axis=0)
                                for el, e in m.boundary[lado])
                    np.testing.assert_allclose(total, np.array(n) * E, atol=1e-14)

    def test_traccion_lineal_en_borde_de_dos_nodos(self):
        for nombre in ("Quad4", "Tri3"):
            with self.subTest(elemento=nombre):
                el = elemento(nombre)
                a, c = el.EDGE_NODES[0]
                X = el.get_coordinate_matrix(ndim=2)
                L = np.linalg.norm(X[c] - X[a])
                campo = lambda x, n: np.array([x[0], 3.0 - x[1]])  # noqa: E731
                ta, tc = campo(X[a], None), campo(X[c], None)
                f = el.compute_edge_traction(0, campo)
                cerca(nodal(f, a), L * (2 * ta + tc) / 6 * E)
                cerca(nodal(f, c), L * (ta + 2 * tc) / 6 * E)

    def test_borde_fuera_de_rango(self):
        for nombre in REF:
            with self.subTest(elemento=nombre):
                el = elemento(nombre)
                with self.assertRaises(ValueError):
                    el.compute_edge_traction(len(el.EDGE_NODES), np.zeros(2))


if __name__ == "__main__":
    unittest.main()
