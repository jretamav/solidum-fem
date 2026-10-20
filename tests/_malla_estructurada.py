"""Mallas de prueba por coordenadas: nodos compartidos fundidos y cajas
estructuradas de Quad4/Hex8 (o de sus gemelos térmicos).

Lo usan las pruebas de la red de seguridad del ADR 0019
(``test_safety_net_initial_tangent.py``, ``test_internal_mechanisms.py``).
"""
import numpy as np

from solidum.core.domain import Domain
from solidum.registry import ElementRegistry

# Orden local de los nodos en el cuadrado / cubo unidad de referencia.
QUAD4 = ((0, 0), (1, 0), (1, 1), (0, 1))
HEX8 = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))


class Malla:
    """Nodos por coordenadas (los compartidos se funden) y mallas
    estructuradas de Quad4/Hex8 sobre cajas."""

    def __init__(self):
        self.dom = Domain()
        self._ids = {}

    def nodo(self, x):
        clave = tuple(round(float(v), 9) for v in x)
        if clave not in self._ids:
            nid = len(self.dom.nodes) + 1
            while nid in self.dom.nodes:
                nid += 1
            self._ids[clave] = nid
            self.dom.add_node(nid, [float(v) for v in x])
        return self.dom.nodes[self._ids[clave]]

    def elemento(self, tipo, coords, *args, **kwargs):
        el = ElementRegistry.get(tipo)(len(self.dom.elements) + 1,
                                       [self.nodo(c) for c in coords], *args, **kwargs)
        self.dom.add_element(el)
        return el

    def caja(self, tipo, origen, tam, n, material, **kwargs):
        ref = QUAD4 if len(n) == 2 else HEX8
        h = np.asarray(tam, float) / np.asarray(n, float)
        for celda in np.ndindex(*n):
            o = np.asarray(origen, float) + np.asarray(celda) * h
            self.elemento(tipo, [o + np.asarray(r) * h for r in ref], material, **kwargs)

    def nodos(self, condicion):
        return [self.dom.nodes[i] for x, i in self._ids.items() if condicion(np.asarray(x))]

    def cerrar(self, cargas=()):
        self.dom.generate_equation_numbers(verbose=False)
        F = np.zeros(self.dom.total_dofs)
        for nodo, dof, valor in cargas:
            F[nodo.dofs[dof]] += valor
        return self.dom, F


def empotrar(nodos, dofs):
    for nodo in nodos:
        for dof in dofs:
            nodo.fix_dof(dof, 0.0)
