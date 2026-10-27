"""Contrato ``REFERENCE_CELL`` / ``REFERENCE_NODES`` sobre todo el registro (ADR 0022).

Con esos dos datos el generador de mallas estructuradas (``solidum.mesh``)
coloca los nodos de cualquier elemento y la verificación por soluciones
manufacturadas deduce su orden de convergencia. Un dato mal declarado daría
mallas con los nodos cambiados de sitio, así que se contrasta con lo único
que no puede mentir: las funciones de forma del propio elemento valen
``N_i(ξ_j) = δ_ij`` en los nodos declarados.

Además, los bordes (``EDGE_NODES``) deben ser lados de la celda de
referencia recorridos en sentido antihorario, con el nodo medio en el centro:
es lo que suponen la carga de borde (normal exterior = tangente girada −90°)
y el generador (bordes sobre el contorno del bloque).

Y ningún sólido plano del programa principal puede quedarse sin declararlos:
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
from solidum.mesh import structured_mesh  # noqa: E402
from solidum.registry import ElementRegistry  # noqa: E402

CELDAS_2D = {"quad", "tri"}
MAT = Elastic2D(1.0, 0.25)


def declarados():
    return [ElementRegistry.get(n) for n in ElementRegistry.names()
            if ElementRegistry.get(n).REFERENCE_NODES is not None]


def solidos_planos():
    """Sólidos planos del programa principal (los módulos de usuario pueden
    estar cargados por otras pruebas de la misma sesión)."""
    return [ElementRegistry.get(n) for n in ElementRegistry.names()
            if ElementRegistry.get(n).__module__.startswith("solidum.elements.")
            and list(getattr(ElementRegistry.get(n), "DOF_NAMES", [])) == ["ux", "uy"]
            and getattr(ElementRegistry.get(n), "STRAIN_DIM", None) == 3]


def un_elemento(cls):
    return next(iter(structured_mesh(cls, (1, 1), material=MAT).domain.elements.values()))


class TestNodosDeReferencia(unittest.TestCase):

    def test_hay_elementos_declarados(self):
        self.assertGreaterEqual({c.__name__ for c in declarados()},
                                {"Quad4", "Quad8", "Quad9", "Tri3", "Tri6"})

    def test_todo_solido_plano_declara_su_celda(self):
        for cls in solidos_planos():
            with self.subTest(elemento=cls.__name__):
                self.assertIn(cls.REFERENCE_CELL, CELDAS_2D)
                self.assertIsNotNone(cls.REFERENCE_NODES)

    def test_funciones_de_forma_valen_delta_en_los_nodos(self):
        for cls in declarados():
            if cls.REFERENCE_CELL not in CELDAS_2D:
                continue
            with self.subTest(elemento=cls.__name__):
                nodos = np.asarray(cls.REFERENCE_NODES, dtype=float)
                N = un_elemento(cls).batch_shape_functions(nodos)
                np.testing.assert_allclose(N, np.eye(len(nodos)), atol=1e-12)

    def test_bordes_antihorarios_sobre_los_lados_de_la_celda(self):
        for cls in declarados():
            if cls.REFERENCE_CELL not in CELDAS_2D:
                continue
            nodos = np.asarray(cls.REFERENCE_NODES, dtype=float)
            centro = nodos[:3].mean(axis=0) if cls.REFERENCE_CELL == "tri" else np.zeros(2)
            for e, borde in enumerate(cls.EDGE_NODES):
                with self.subTest(elemento=cls.__name__, borde=e):
                    a, c = nodos[borde[0]], nodos[borde[-1]]
                    d = c - a
                    normal = np.array([d[1], -d[0]])
                    self.assertGreater(normal @ ((a + c) / 2 - centro), 0.0,
                                       "la normal tangente·(−90°) debe salir de la celda")
                    if len(borde) == 3:
                        np.testing.assert_allclose(nodos[borde[1]], (a + c) / 2, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
