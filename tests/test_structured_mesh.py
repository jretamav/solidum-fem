"""Generador de mallas estructuradas ``solidum.mesh.structured_mesh`` (ADR 0022).

Se comprueba lo que una malla debe cumplir con independencia del elemento:
cuántos nodos y elementos tiene (los compartidos se funden), que cubre el
dominio sin huecos ni solapes (la suma de las áreas de los elementos es el
área del dominio), que todos los jacobianos son positivos (orientación
antihoraria), que los bordes de contorno están donde dicen y que la
transformación lleva los nodos intermedios sobre la geometría curva.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401,E402
from solidum.core.domain import Domain  # noqa: E402
from solidum.materials.elastic_2d import Elastic2D  # noqa: E402
from solidum.mesh import SIDES_2D, structured_mesh  # noqa: E402
from solidum.registry import QuadratureRegistry  # noqa: E402

MAT = Elastic2D(200e3, 0.3)
N = 3

# Nodos y elementos de una malla de N×N celdas.
CONTEO = {
    "Quad4": ((N + 1) ** 2, N * N),
    "Quad8": ((2 * N + 1) ** 2 - N * N, N * N),
    "Quad9": ((2 * N + 1) ** 2, N * N),
    "Tri3": ((N + 1) ** 2, 2 * N * N),
    "Tri6": ((2 * N + 1) ** 2, 2 * N * N),
}


def corona(X):
    r, t = 1.0 + X[0], X[1] * np.pi / 2
    return np.array([r * np.cos(t), r * np.sin(t)])


def areas(mesh):
    """Área de cada elemento con la regla de su carga (exacta para detJ
    polinómico de grado ≤ 2·p) y el menor jacobiano encontrado."""
    total, minimo = 0.0, np.inf
    for el in mesh.domain.elements.values():
        pts, wts = QuadratureRegistry.get(el._LOAD_QUADRATURE)
        pts = np.asarray(pts, dtype=float).reshape(len(wts), 2)
        coords = np.ascontiguousarray(el.get_coordinate_matrix(ndim=2))
        B = np.empty((3, 2 * len(el.nodes)))
        detJ = np.array([type(el).BATCH_KINEMATICS(p, coords, B) for p in pts])
        total += detJ @ np.asarray(wts)
        minimo = min(minimo, detJ.min())
    return total, minimo


class TestMallaEstructurada(unittest.TestCase):

    def test_conteo_de_nodos_y_elementos(self):
        for nombre, (nn, ne) in CONTEO.items():
            with self.subTest(elemento=nombre):
                m = structured_mesh(nombre, (N, N), material=MAT)
                self.assertEqual(len(m.domain.nodes), nn)
                self.assertEqual(len(m.domain.elements), ne)

    def test_sin_nodos_duplicados(self):
        for nombre in CONTEO:
            with self.subTest(elemento=nombre):
                m = structured_mesh(nombre, (N, N), mapping=corona, material=MAT)
                X = np.array([n.coordinates for n in m.domain.nodes.values()])
                d = np.linalg.norm(X[:, None] - X[None], axis=-1) + np.eye(len(X))
                self.assertGreater(d.min(), 1e-6)

    def test_cubre_la_caja_con_jacobianos_positivos(self):
        for nombre in CONTEO:
            with self.subTest(elemento=nombre):
                m = structured_mesh(nombre, (N, 2), size=(2.0, 0.5), origin=(1.0, -1.0),
                                    material=MAT, thickness=0.1)
                area, minimo = areas(m)
                self.assertAlmostEqual(area, 1.0, places=12)
                self.assertGreater(minimo, 0.0)

    def test_corona_curva(self):
        """Los nodos (también los intermedios) caen sobre las circunferencias, y
        el área converge a la de la corona: con orden 2 los elementos lineales
        (cuerdas) y con orden 4 los cuadráticos (arcos parabólicos)."""
        exacta = 3 * np.pi / 4
        for nombre in CONTEO:
            with self.subTest(elemento=nombre):
                m = structured_mesh(nombre, (N, N), mapping=corona, material=MAT)
                for lado, r in (("xmin", 1.0), ("xmax", 2.0)):
                    for nodo in m.nodes_on(lado):
                        self.assertAlmostEqual(np.hypot(*nodo.coordinates[:2]), r, places=12)
                e = [abs(areas(structured_mesh(nombre, (k, k), mapping=corona,
                                               material=MAT))[0] - exacta) for k in (4, 8)]
                orden = np.log2(e[0] / e[1])
                self.assertGreater(orden, 1.9 if nombre in ("Quad4", "Tri3") else 3.9)

    def test_bordes_de_contorno(self):
        for nombre in CONTEO:
            with self.subTest(elemento=nombre):
                m = structured_mesh(nombre, (N, N + 1), size=(3.0, 2.0), material=MAT)
                for lado in SIDES_2D:
                    eje, valor = "xy".index(lado[0]), (0.0 if lado.endswith("min") else
                                                       (3.0, 2.0)["xy".index(lado[0])])
                    bordes = m.boundary_edges[lado]
                    self.assertEqual(len(bordes), N + 1 if eje == 0 else N)
                    for el, e in bordes:
                        for k in el.EDGE_NODES[e]:
                            self.assertAlmostEqual(el.nodes[k].coordinates[eje], valor, places=12)
                    nodos = {n.id for el, e in bordes for n in (el.nodes[k] for k in el.EDGE_NODES[e])}
                    self.assertEqual(nodos, {n.id for n in m.nodes_on(lado)})

    def test_anade_a_un_dominio_existente(self):
        dom = Domain()
        dom.add_node(7, [10.0, 10.0])
        m = structured_mesh("Quad4", (1, 1), domain=dom, material=MAT)
        self.assertIs(m.domain, dom)
        self.assertEqual(sorted(dom.nodes), [7, 8, 9, 10, 11])

    def test_rechazos(self):
        with self.assertRaises(ValueError):
            structured_mesh("Truss2D", (2, 2), material=MAT)
        with self.assertRaises(ValueError):
            structured_mesh("Quad4", (0, 2), material=MAT)
        with self.assertRaises(ValueError):
            structured_mesh("Quad4", (1, 1), material=MAT).nodes_on("izquierda")


if __name__ == "__main__":
    unittest.main()
