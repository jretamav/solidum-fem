"""Armadura Pratt del ejemplo 11 del manual de ejemplos, para las pruebas de
la red de seguridad del análisis estático (ADR 0019).

Cuatro paneles de 3 m × 3 m, pasador en el nudo inferior izquierdo y rodillo
vertical en el derecho (``examples/modelo_mal_planteado/run.py``). Con
``sin_diagonal=True`` falta la diagonal del segundo panel: un cuadrilátero
articulado, mecanismo **interno** (los apoyos impiden todo movimiento rígido,
así que la capa 1 no lo ve). La carga ``"horizontal"`` (P en el rodillo) no
lo activa: la cuerda inferior la lleva al pasador y el sistema es singular
pero consistente; la ``"vertical"`` (P en los tres nudos inferiores
interiores) lo activa.

Se reescribe aquí, en vez de importar el script del ejemplo, para que las
pruebas de la red no dependan de un archivo de documentación.
"""
from __future__ import annotations

import numpy as np

from solidum.core.domain import Domain
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D

PANEL, ALTO, N_PANELES = 3.0, 3.0, 4     # [m]
E, AREA, P = 200.0e9, 20.0e-4, 100.0e3   # acero, 20 cm², 100 kN


def armadura(sin_diagonal: bool = False, carga: str = "horizontal",
             elemento=Truss2D, material=None):
    """``(dominio, F)`` de la armadura. ``elemento`` es la clase de barra
    (``Truss2D``, ``Truss2DCorot``…); ``material``, el de todas las barras
    (acero elástico por omisión)."""
    s = N_PANELES + 1
    xy = [(PANEL * i, 0.0) for i in range(s)] + [(PANEL * i, ALTO) for i in range(s)]
    barras = ([(i, i + 1) for i in range(N_PANELES)]               # cuerda inferior
              + [(s + i, s + i + 1) for i in range(N_PANELES)]     # cuerda superior
              + [(i, s + i) for i in range(s)])                    # montantes
    for i in range(N_PANELES):
        if sin_diagonal and i == 1:
            continue                                               # falta esta diagonal
        barras.append((s + i, i + 1) if i < N_PANELES // 2 else (i, s + i + 1))

    dom = Domain()
    mat = material if material is not None else Elastic1D(E=E)
    nodos = [dom.add_node(k + 1, list(p)) for k, p in enumerate(xy)]
    for a, b in barras:
        dom.add_element(elemento(len(dom.elements) + 1, [nodos[a], nodos[b]], mat, A=AREA))
    inferiores = nodos[:s]
    inferiores[0].fix_dof("ux", 0.0)
    inferiores[0].fix_dof("uy", 0.0)
    inferiores[-1].fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)

    F = np.zeros(dom.total_dofs)
    if carga == "vertical":
        for nodo in inferiores[1:-1]:
            F[nodo.dofs["uy"]] = -P
    elif carga == "horizontal":
        F[inferiores[-1].dofs["ux"]] = P
    else:
        raise ValueError(f"carga desconocida: {carga!r}")
    return dom, F


def nodo_de_control(dom):
    """Nudo inferior derecho (el del rodillo), cuyo ``ux`` crece con la carga
    horizontal: sirve de magnitud controlada al control indirecto."""
    return dom.nodes[N_PANELES + 1]
