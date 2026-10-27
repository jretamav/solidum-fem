"""Generador de mallas estructuradas ``solidum.mesh.structured_mesh`` (ADR 0022).

Se comprueba lo que una malla debe cumplir con independencia del elemento,
en 2D y en 3D: cuántos nodos y elementos tiene (los compartidos se funden),
que cubre el dominio sin huecos ni solapes (la suma de las áreas o volúmenes
de los elementos es la del dominio), que todos los jacobianos son positivos
(orientación correcta, también en los seis tetraedros de Kuhn de cada cubo),
que los bordes y caras de contorno están donde dicen y que la transformación
lleva los nodos intermedios sobre la geometría curva.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401,E402
from solidum.core.domain import Domain  # noqa: E402
from solidum.materials.elastic_2d import Elastic2D  # noqa: E402
from solidum.materials.elastic_3d import Elastic3D  # noqa: E402
from solidum.mesh import SIDES_2D, SIDES_3D, structured_mesh  # noqa: E402
from solidum.registry import ElementRegistry, QuadratureRegistry  # noqa: E402

MAT = {2: Elastic2D(200e3, 0.3), 3: Elastic3D(200e3, 0.3)}


def conteo(nombre, n):
    """Nodos y elementos de una malla de n celdas por dirección."""
    c = n + 1                     # vértices por dirección
    q = 2 * n + 1                 # vértices y nodos medios por dirección
    return {
        "Quad4": (c ** 2, n ** 2), "Quad8": (q ** 2 - n ** 2, n ** 2), "Quad9": (q ** 2, n ** 2),
        "Tri3": (c ** 2, 2 * n ** 2), "Tri6": (q ** 2, 2 * n ** 2),
        "Hex8": (c ** 3, n ** 3), "Hex20": (c ** 3 + 3 * n * c ** 2, n ** 3), "Hex27": (q ** 3, n ** 3),
        "Tet4": (c ** 3, 6 * n ** 3), "Tet10": (q ** 3, 6 * n ** 3),
    }[nombre]


ELEMENTOS = ("Quad4", "Quad8", "Quad9", "Tri3", "Tri6", "Hex8", "Hex20", "Hex27", "Tet4", "Tet10")
LINEALES = ("Quad4", "Tri3", "Hex8", "Tet4")


def dim(nombre):
    return 2 if ElementRegistry.get(nombre).REFERENCE_CELL in ("quad", "tri") else 3


def malla(nombre, n, **kw):
    d = dim(nombre)
    return structured_mesh(nombre, (n,) * d if isinstance(n, int) else n, material=MAT[d], **kw)


def corona(X):
    """Cuarto de corona 1 ≤ r ≤ 2 en 2D; en 3D, sector de cilindro hueco de
    altura 1."""
    r, t = 1.0 + X[0], X[1] * np.pi / 2
    return np.array([r * np.cos(t), r * np.sin(t), *X[2:]])


def medidas(mesh):
    """Área o volumen total con la regla de la carga de cada elemento y el
    menor jacobiano encontrado."""
    d = mesh.dim
    total, minimo = 0.0, np.inf
    for el in mesh.domain.elements.values():
        pts, wts = QuadratureRegistry.get(el._LOAD_QUADRATURE)
        pts = np.asarray(pts, dtype=float).reshape(len(wts), d)
        coords = np.ascontiguousarray(el.get_coordinate_matrix(ndim=d))
        B = np.empty((3 if d == 2 else 6, d * len(el.nodes)))
        detJ = np.array([type(el).BATCH_KINEMATICS(p, coords, B) for p in pts])
        total += detJ @ np.asarray(wts)
        minimo = min(minimo, detJ.min())
    return total, minimo


class TestMallaEstructurada(unittest.TestCase):

    def test_conteo_de_nodos_y_elementos(self):
        for nombre in ELEMENTOS:
            n = 3 if dim(nombre) == 2 else 2
            with self.subTest(elemento=nombre):
                m = malla(nombre, n)
                self.assertEqual((len(m.domain.nodes), len(m.domain.elements)), conteo(nombre, n))

    def test_sin_nodos_duplicados(self):
        for nombre in ELEMENTOS:
            with self.subTest(elemento=nombre):
                m = malla(nombre, 2, mapping=corona)
                X = np.array([n.coordinates for n in m.domain.nodes.values()])
                d = np.linalg.norm(X[:, None] - X[None], axis=-1) + np.eye(len(X))
                self.assertGreater(d.min(), 1e-6)

    def test_cubre_la_caja_con_jacobianos_positivos(self):
        for nombre in ELEMENTOS:
            with self.subTest(elemento=nombre):
                if dim(nombre) == 2:
                    m = malla(nombre, (3, 2), size=(2.0, 0.5), origin=(1.0, -1.0), thickness=0.1)
                else:
                    m = malla(nombre, (2, 3, 1), size=(2.0, 0.5, 4.0), origin=(1.0, -1.0, 0.5))
                medida, minimo = medidas(m)
                self.assertAlmostEqual(medida, 1.0 if dim(nombre) == 2 else 4.0, places=12)
                self.assertGreater(minimo, 0.0)

    def test_corona_curva(self):
        """Los nodos (también los intermedios) caen sobre las circunferencias, y
        el área o el volumen converge al exacto, ``3π/4``: con orden 2 los
        elementos lineales (cuerdas) y con orden 4 los cuadráticos (arcos
        parabólicos)."""
        for nombre in ELEMENTOS:
            with self.subTest(elemento=nombre):
                m = malla(nombre, 2, mapping=corona)
                for lado, r in (("xmin", 1.0), ("xmax", 2.0)):
                    for nodo in m.nodes_on(lado):
                        self.assertAlmostEqual(np.hypot(*nodo.coordinates[:2]), r, places=12)
                niveles = (4, 8) if dim(nombre) == 2 else (2, 4)
                e = [abs(medidas(malla(nombre, k, mapping=corona))[0] - 3 * np.pi / 4)
                     for k in niveles]
                self.assertGreater(np.log2(e[0] / e[1]), 1.9 if nombre in LINEALES else 3.9)

    def test_contorno(self):
        """Cada lado está formado por las entidades de contorno que dice
        (bordes en 2D, caras en 3D), y sus nodos son los de ``nodes_on``."""
        for nombre in ELEMENTOS:
            with self.subTest(elemento=nombre):
                d = dim(nombre)
                n = (3, 4) if d == 2 else (2, 3, 1)
                size = (3.0, 2.0) if d == 2 else (3.0, 2.0, 1.5)
                m = malla(nombre, n, size=size)
                cls = ElementRegistry.get(nombre)
                locales = cls.EDGE_NODES if d == 2 else cls.FACE_NODES
                por_celda = 2 if cls.REFERENCE_CELL == "tet" else 1
                for lado in (SIDES_2D if d == 2 else SIDES_3D):
                    eje = "xyz".index(lado[0])
                    valor = 0.0 if lado.endswith("min") else size[eje]
                    entidades = m.boundary[lado]
                    celdas = int(np.prod([n[k] for k in range(d) if k != eje]))
                    self.assertEqual(len(entidades), por_celda * celdas)
                    for el, k in entidades:
                        for j in locales[k]:
                            self.assertAlmostEqual(el.nodes[j].coordinates[eje], valor, places=12)
                    nodos = {el.nodes[j].id for el, k in entidades for j in locales[k]}
                    self.assertEqual(nodos, {nd.id for nd in m.nodes_on(lado)})

    def test_anade_a_un_dominio_existente(self):
        dom = Domain()
        dom.add_node(7, [10.0, 10.0])
        m = structured_mesh("Quad4", (1, 1), domain=dom, material=MAT[2])
        self.assertIs(m.domain, dom)
        self.assertEqual(sorted(dom.nodes), [7, 8, 9, 10, 11])

    def test_rechazos(self):
        with self.assertRaises(ValueError):
            structured_mesh("Truss2D", (2, 2), material=MAT[2])
        with self.assertRaises(ValueError):
            structured_mesh("Quad4", (0, 2), material=MAT[2])
        with self.assertRaises(ValueError):
            structured_mesh("Hex8", (2, 2), material=MAT[3])
        with self.assertRaises(ValueError):
            structured_mesh("Quad4", (1, 1), material=MAT[2]).nodes_on("zmin")


if __name__ == "__main__":
    unittest.main()
