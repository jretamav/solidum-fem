"""Contrato ``REFERENCE_CELL`` / ``REFERENCE_NODES`` sobre todo el registro (ADR 0022).

Con esos dos datos el generador de mallas estructuradas (``solidum.mesh``)
coloca los nodos de cualquier elemento y la verificación por soluciones
manufacturadas deduce su orden de convergencia. Un dato mal declarado daría
mallas con los nodos cambiados de sitio, así que se contrasta con lo único
que no puede mentir: las funciones de forma del propio elemento valen
``N_i(ξ_j) = δ_ij`` en los nodos declarados.

Además, los bordes (``EDGE_NODES``, 2D) y las caras (``FACE_NODES``, 3D)
deben estar sobre el contorno de la celda de referencia, recorridos en
sentido antihorario visto desde fuera, y con sus nodos donde los esperan las
funciones de forma del borde o de la cara: es lo que suponen la carga de
borde o de cara (la normal exterior sale de esa orientación) y el generador
(entidades sobre el contorno del bloque).

Y ningún sólido del programa principal puede quedarse sin declararlos:
quedaría fuera, sin aviso, del barrido de soluciones manufacturadas. Los
elementos de los módulos de usuario (ADR 0020) los declaran si quieren usar
el generador o la verificación; su barrido es suyo.
"""
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum  # noqa: F401,E402  — dispara el autodiscover que puebla el registro
from solidum.materials.elastic_2d import Elastic2D  # noqa: E402
from solidum.materials.elastic_3d import Elastic3D  # noqa: E402
from solidum.mesh import structured_mesh  # noqa: E402
from solidum.registry import ElementRegistry  # noqa: E402

CELDAS = {"quad": 2, "tri": 2, "hex": 3, "tet": 3}
MAT = {2: Elastic2D(1.0, 0.25), 3: Elastic3D(1.0, 0.25)}
DOFS = {2: ["ux", "uy"], 3: ["ux", "uy", "uz"]}
STRAIN = {2: 3, 3: 6}


def declarados():
    return [ElementRegistry.get(n) for n in ElementRegistry.names()
            if ElementRegistry.get(n).REFERENCE_NODES is not None]


def solidos(dim):
    """Sólidos del programa principal (los módulos de usuario pueden estar
    cargados por otras pruebas de la misma sesión)."""
    return [ElementRegistry.get(n) for n in ElementRegistry.names()
            if ElementRegistry.get(n).__module__.startswith("solidum.elements.")
            and list(getattr(ElementRegistry.get(n), "DOF_NAMES", [])) == DOFS[dim]
            and getattr(ElementRegistry.get(n), "STRAIN_DIM", None) == STRAIN[dim]]


def un_elemento(cls):
    d = CELDAS[cls.REFERENCE_CELL]
    m = structured_mesh(cls, (1,) * d, material=MAT[d])
    return next(iter(m.domain.elements.values()))


def centro_de_celda(cls):
    nodos = np.asarray(cls.REFERENCE_NODES, dtype=float)
    if cls.REFERENCE_CELL in ("quad", "hex"):
        return np.zeros(nodos.shape[1])
    return nodos[:nodos.shape[1] + 1].mean(axis=0)


def sobre_el_contorno(cls, X):
    """``True`` si todos los puntos ``X`` están sobre un mismo lado o cara de
    la celda de referencia."""
    if cls.REFERENCE_CELL in ("quad", "hex"):
        return any(np.allclose(X[:, k], s) for k in range(X.shape[1]) for s in (-1.0, 1.0))
    return (any(np.allclose(X[:, k], 0.0) for k in range(X.shape[1]))
            or np.allclose(X.sum(axis=1), 1.0))


class TestNodosDeReferencia(unittest.TestCase):

    def test_hay_elementos_declarados(self):
        self.assertGreaterEqual({c.__name__ for c in declarados()},
                                {"Quad4", "Quad8", "Quad9", "Tri3", "Tri6",
                                 "Hex8", "Hex20", "Hex27", "Tet4", "Tet10"})

    def test_todo_solido_declara_su_celda(self):
        for d in (2, 3):
            for cls in solidos(d):
                with self.subTest(elemento=cls.__name__):
                    self.assertEqual(CELDAS.get(cls.REFERENCE_CELL), d)
                    self.assertIsNotNone(cls.REFERENCE_NODES)

    def test_funciones_de_forma_valen_delta_en_los_nodos(self):
        for cls in declarados():
            if cls.REFERENCE_CELL not in CELDAS:
                continue
            with self.subTest(elemento=cls.__name__):
                nodos = np.asarray(cls.REFERENCE_NODES, dtype=float)
                N = un_elemento(cls).batch_shape_functions(nodos)
                np.testing.assert_allclose(N, np.eye(len(nodos)), atol=1e-12)

    def test_bordes_antihorarios_sobre_los_lados_de_la_celda(self):
        for cls in declarados():
            if CELDAS.get(cls.REFERENCE_CELL) != 2:
                continue
            nodos = np.asarray(cls.REFERENCE_NODES, dtype=float)
            centro = centro_de_celda(cls)
            for e, borde in enumerate(cls.EDGE_NODES):
                with self.subTest(elemento=cls.__name__, borde=e):
                    a, c = nodos[borde[0]], nodos[borde[-1]]
                    self.assertTrue(sobre_el_contorno(cls, nodos[list(borde)]))
                    d = c - a
                    normal = np.array([d[1], -d[0]])
                    self.assertGreater(normal @ ((a + c) / 2 - centro), 0.0,
                                       "la normal tangente·(−90°) debe salir de la celda")
                    if len(borde) == 3:
                        np.testing.assert_allclose(nodos[borde[1]], (a + c) / 2, atol=1e-14)

    def test_caras_antihorarias_vistas_desde_fuera(self):
        """En la celda de referencia las caras son planas: la transformación
        que dan las funciones de forma de la cara sobre sus nodos debe ser la
        lineal (triángulo) o bilineal (cuadrilátero) de sus vértices, lo que
        sólo ocurre con los nodos intermedios en su sitio; y
        ``∂x/∂s × ∂x/∂r`` debe salir de la celda."""
        rng = np.random.default_rng(1)
        for cls in declarados():
            if CELDAS.get(cls.REFERENCE_CELL) != 3:
                continue
            nodos = np.asarray(cls.REFERENCE_NODES, dtype=float)
            centro = centro_de_celda(cls)
            for f, cara in enumerate(cls.FACE_NODES):
                with self.subTest(elemento=cls.__name__, cara=f):
                    X = nodos[list(cara)]
                    self.assertTrue(sobre_el_contorno(cls, X))
                    tri = cls.REFERENCE_CELL == "tet"
                    st = rng.dirichlet(np.ones(3), 6)[:, 1:] if tri else rng.uniform(-1, 1, (6, 2))
                    for s, r in st:
                        x = np.asarray(cls._FACE_N_FN(s, r)) @ X
                        if tri:
                            v = X[0] + s * (X[1] - X[0]) + r * (X[2] - X[0])
                        else:
                            w = np.array([(1 - s) * (1 - r), (1 + s) * (1 - r),
                                          (1 + s) * (1 + r), (1 - s) * (1 + r)]) / 4
                            v = w @ X[:4]
                        np.testing.assert_allclose(x, v, atol=1e-13)
                        dN = np.asarray(cls._FACE_DN_FN(s, r))
                        normal = np.cross(dN[0] @ X, dN[1] @ X)
                        self.assertGreater(normal @ (x - centro), 0.0,
                                           "∂x/∂s × ∂x/∂r debe salir de la celda")


if __name__ == "__main__":
    unittest.main()
