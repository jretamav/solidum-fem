"""Cargas distribuidas variables en el espacio de los sólidos 3D (ADR 0022, fase 2).

``compute_body_load`` acepta un campo ``b(x)`` y ``compute_face_traction`` un
campo ``t(x, n)``, con ``n`` la normal exterior de la cara discreta. Las
comprobaciones son analíticas:

- fuerza de cuerpo lineal en un elemento afín: resultante ``b(centroide)·V``;
- **vector área**: con ``t = n``, la suma de las fuerzas nodales es
  ``∫ n dA``, que sobre cualquier cara, plana o curva, vale
  ``½ ∮ x × dx`` a lo largo de su contorno recorrido en sentido antihorario
  visto desde fuera. Se comprueba con las caras de elementos de nodos
  perturbados al azar (caras alabeadas y bordes curvos): verifica a la vez la
  integral sobre la geometría real y el sentido de la normal; una presión
  uniforme tiene resultante ``−p`` por el vector área;
- tracción lineal sobre una cara del Tet4:
  ``f_k = A·(2t_k + Σ_{j≠k} t_j)/12`` (``∫ N_i N_j dA = A(1 + δ_ij)/12``).
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401,E402
from solidum.core.node import Node  # noqa: E402
from solidum.materials.elastic_3d import Elastic3D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.mesh import SIDES_3D, structured_mesh  # noqa: E402
from solidum.registry import ElementRegistry  # noqa: E402

MAT = Elastic3D(200e3, 0.3)
REF = ("Hex8", "Hex20", "Hex27", "Tet4", "Tet10")
M = np.array([[2.0, 0.3, 0.1], [0.4, 1.5, -0.2], [0.0, 0.3, 1.2]])


def elemento(nombre, perturbacion=0.0, semilla=0):
    """Elemento suelto: sus nodos de referencia llevados a un paralelepípedo
    (``x = ξ·M + c``) y, si se pide, movidos al azar."""
    cls = ElementRegistry.get(nombre)
    xi = np.asarray(cls.REFERENCE_NODES, dtype=float)
    if cls.REFERENCE_CELL == "hex":
        xi = (xi + 1.0) / 2.0
    X = xi @ M + np.array([1.0, -0.5, 0.2])
    X += perturbacion * np.random.default_rng(semilla).uniform(-1.0, 1.0, X.shape)
    return cls(1, [Node(i + 1, list(x)) for i, x in enumerate(X)], MAT)


def volumen(nombre):
    return abs(np.linalg.det(M)) / (1.0 if ElementRegistry.get(nombre).REFERENCE_CELL == "hex" else 6.0)


def cerca(actual, esperado, rtol=1e-13):
    """``assert_allclose`` con tolerancia absoluta escalada (puede haber
    componentes nulas, cuyo redondeo cambia con la plataforma)."""
    esperado = np.asarray(esperado)
    np.testing.assert_allclose(actual, esperado, rtol=rtol, atol=rtol * np.abs(esperado).max())


def bordes_de_cara(n_nodos):
    """Bordes del contorno de una cara en el orden de ``FACE_NODES``, como
    índices locales de la cara ``(inicio, [medio,] fin)``."""
    if n_nodos in (3, 4):
        return [(k, (k + 1) % n_nodos) for k in range(n_nodos)]
    vertices = 3 if n_nodos == 6 else 4
    return [(k, vertices + k, (k + 1) % vertices) for k in range(vertices)]


def vector_area(X_cara):
    """``½ ∮ x × dx`` sobre el contorno de la cara (bordes de Lagrange de 2 ó
    3 nodos; Gauss de 3 puntos, exacto para ``x × x'`` de grado 3)."""
    s, w = np.polynomial.legendre.leggauss(3)
    total = np.zeros(3)
    for borde in bordes_de_cara(len(X_cara)):
        X = X_cara[list(borde)]
        if len(borde) == 2:
            N = np.column_stack([(1 - s) / 2, (1 + s) / 2])
            dN = np.column_stack([-0.5 * np.ones(3), 0.5 * np.ones(3)])
        else:
            N = np.column_stack([s * (s - 1) / 2, 1 - s * s, s * (s + 1) / 2])
            dN = np.column_stack([s - 0.5, -2 * s, s + 0.5])
        total += 0.5 * (w[:, None] * np.cross(N @ X, dN @ X)).sum(axis=0)
    return total


class TestFuerzaDeCuerpo(unittest.TestCase):

    def test_campo_constante_igual_a_la_constante(self):
        b = np.array([3.0, -7.0, 2.0])
        for nombre in REF:
            with self.subTest(elemento=nombre):
                el = elemento(nombre)
                cerca(el.compute_body_load(lambda x: b), el.compute_body_load(b))

    def test_campo_lineal(self):
        """Resultante ``∫ b dV = b(centroide)·V`` en un elemento afín."""
        g = np.array([[2.0, -1.0, 0.5], [0.5, 3.0, 0.0], [-1.0, 0.2, 1.0]])
        campo = lambda x: g @ x + np.array([1.0, 0.0, -2.0])  # noqa: E731
        for nombre in REF:
            with self.subTest(elemento=nombre):
                el = elemento(nombre)
                X = el.get_coordinate_matrix(ndim=3)
                vertices = X[:8] if ElementRegistry.get(nombre).REFERENCE_CELL == "hex" else X[:4]
                f = el.compute_body_load(campo).reshape(-1, 3).sum(axis=0)
                cerca(f, campo(vertices.mean(axis=0)) * volumen(nombre))

    def test_ensamblador_acepta_un_campo(self):
        m = structured_mesh("Tet10", (2, 1, 1), material=MAT)
        m.domain.generate_equation_numbers(verbose=False)
        campo = lambda x: np.array([x[1], -2.0, x[0] * x[2]])  # noqa: E731
        F = Assembler(m.domain).assemble_body_load(campo)
        G = np.zeros_like(F)
        for el in m.domain.elements.values():
            G[el.get_global_dof_indices()] += el.compute_body_load(campo)
        cerca(F, G)


class TestCargaDeCara(unittest.TestCase):

    def test_vector_area_en_caras_alabeadas(self):
        for nombre in REF:
            el = elemento(nombre, perturbacion=0.08, semilla=3)
            X = el.get_coordinate_matrix(ndim=3)
            for f, cara in enumerate(el.FACE_NODES):
                with self.subTest(elemento=nombre, cara=f):
                    total = el.compute_face_traction(f, lambda x, n: n).reshape(-1, 3).sum(axis=0)
                    cerca(total, vector_area(X[list(cara)]), rtol=1e-12)

    def test_presion_uniforme(self):
        p = 5.0
        for nombre in REF:
            el = elemento(nombre, perturbacion=0.08, semilla=4)
            X = el.get_coordinate_matrix(ndim=3)
            for f, cara in enumerate(el.FACE_NODES):
                with self.subTest(elemento=nombre, cara=f):
                    total = el.compute_face_traction(f, lambda x, n: -p * n).reshape(-1, 3).sum(axis=0)
                    cerca(total, -p * vector_area(X[list(cara)]), rtol=1e-12)

    def test_normal_exterior_en_un_cubo(self):
        for nombre in REF:
            with self.subTest(elemento=nombre):
                m = structured_mesh(nombre, (1, 1, 1), material=MAT)
                for lado in SIDES_3D:
                    n = np.zeros(3)
                    n["xyz".index(lado[0])] = -1.0 if lado.endswith("min") else 1.0
                    total = sum(el.compute_face_traction(k, lambda x, nn: nn).reshape(-1, 3).sum(axis=0)
                                for el, k in m.boundary[lado])
                    cerca(total, n)

    def test_traccion_lineal_en_cara_del_tet4(self):
        el = elemento("Tet4", perturbacion=0.1, semilla=5)
        X = el.get_coordinate_matrix(ndim=3)
        campo = lambda x, n: np.array([x[0], 3.0 - x[1], x[0] + 2 * x[2]])  # noqa: E731
        for f, cara in enumerate(el.FACE_NODES):
            with self.subTest(cara=f):
                Xc = X[list(cara)]
                A = 0.5 * np.linalg.norm(np.cross(Xc[1] - Xc[0], Xc[2] - Xc[0]))
                t = np.array([campo(x, None) for x in Xc])
                f_el = el.compute_face_traction(f, campo)
                for k, i in enumerate(cara):
                    cerca(f_el[3 * i:3 * i + 3], A * (t[k] + t.sum(axis=0)) / 12)

    def test_cara_fuera_de_rango(self):
        for nombre in REF:
            with self.subTest(elemento=nombre):
                el = elemento(nombre)
                with self.assertRaises(ValueError):
                    el.compute_face_traction(len(el.FACE_NODES), np.zeros(3))


if __name__ == "__main__":
    unittest.main()
