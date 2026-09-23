"""Malla del ejemplo 2: cuarto de placa con agujero circular (Gmsh).

Genera ``placa_agujero.msh`` (versionado junto al ejemplo, para que resolverlo
no requiera Gmsh). Sólo hace falta volver a ejecutarlo si cambia la geometría
o la densidad de malla::

    python examples/placa_agujero_kirsch/malla.py

Geometría (cuarto superior derecho por doble simetría): agujero de radio
``a`` centrado en el origen, placa de ``W × H``. Malla **estructurada** de
cuadriláteros (Quad4) en cuatro parches transfinitos:

- un **anillo polar** ``a ≤ r ≤ R1``, partido a 45°: líneas radiales rectas y
  circunferencias concéntricas, con elementos que crecen en progresión
  geométrica desde el agujero. Ahí están los gradientes de esfuerzo fuertes.
- la **transición** del círculo ``r = R1`` al contorno cuadrado, partida por la
  diagonal hacia la esquina ``(W, H)``.

Una primera versión sin anillo (dos parches del agujero al contorno) daba
líneas "circunferenciales" quebradas sobre la diagonal, porque la diagonal es
más larga que los bordes y reparte sus nodos a otras distancias del agujero:
el error en los puntos de Gauss vecinos llegaba al 11 % de σ∞.

Grupos físicos: ``Placa`` (superficies), ``Simetria_X`` (x = 0),
``Simetria_Y`` (y = 0), ``Borde_Carga`` (x = W), ``Agujero``.
"""
from __future__ import annotations

import math
from pathlib import Path

AQUI = Path(__file__).resolve().parent

A = 0.01           # radio del agujero [m]
W = H = 20 * A     # cuarto de placa: 0.2 m × 0.2 m
R1 = 5 * A         # radio exterior del anillo polar
N_ARCO = 24        # elementos en cada octavo de circunferencia
N_ANILLO = 36      # elementos radiales en el anillo (a → R1)
Q_ANILLO = 1.06    # razón de tamaños radiales en el anillo
N_EXT = 16         # elementos radiales en la transición (R1 → contorno)
Q_EXT = 1.15


def generar(ruta: Path = AQUI / "placa_agujero.msh") -> Path:
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("placa_agujero")
    geo = gmsh.model.geo
    c = math.cos(math.pi / 4)

    o = geo.addPoint(0, 0, 0)
    pa, pm, pe = (geo.addPoint(A, 0, 0), geo.addPoint(A * c, A * c, 0),
                  geo.addPoint(0, A, 0))                       # agujero
    qa, qm, qe = (geo.addPoint(R1, 0, 0), geo.addPoint(R1 * c, R1 * c, 0),
                  geo.addPoint(0, R1, 0))                      # r = R1
    pb, pc, pd = geo.addPoint(W, 0, 0), geo.addPoint(W, H, 0), geo.addPoint(0, H, 0)

    agujero = [geo.addCircleArc(pa, o, pm), geo.addCircleArc(pm, o, pe)]
    circulo = [geo.addCircleArc(qa, o, qm), geo.addCircleArc(qm, o, qe)]
    # Radiales, siempre del centro hacia fuera (la progresión crece hacia fuera).
    ra_in, rm_in, re_in = geo.addLine(pa, qa), geo.addLine(pm, qm), geo.addLine(pe, qe)
    ra_ex, rm_ex, re_ex = geo.addLine(qa, pb), geo.addLine(qm, pc), geo.addLine(qe, pd)
    bc, cd = geo.addLine(pb, pc), geo.addLine(pc, pd)

    parches = [
        (geo.addPlaneSurface([geo.addCurveLoop([ra_in, circulo[0], -rm_in, -agujero[0]])]),
         [pa, qa, qm, pm]),
        (geo.addPlaneSurface([geo.addCurveLoop([rm_in, circulo[1], -re_in, -agujero[1]])]),
         [pm, qm, qe, pe]),
        (geo.addPlaneSurface([geo.addCurveLoop([ra_ex, bc, -rm_ex, -circulo[0]])]),
         [qa, pb, pc, qm]),
        (geo.addPlaneSurface([geo.addCurveLoop([rm_ex, cd, -re_ex, -circulo[1]])]),
         [qm, pc, pd, qe]),
    ]

    for curva in agujero + circulo + [bc, cd]:
        geo.mesh.setTransfiniteCurve(curva, N_ARCO + 1)
    for curva in (ra_in, rm_in, re_in):
        geo.mesh.setTransfiniteCurve(curva, N_ANILLO + 1, "Progression", Q_ANILLO)
    for curva in (ra_ex, rm_ex, re_ex):
        geo.mesh.setTransfiniteCurve(curva, N_EXT + 1, "Progression", Q_EXT)
    for superficie, esquinas in parches:
        geo.mesh.setTransfiniteSurface(superficie, cornerTags=esquinas)
        geo.mesh.setRecombine(2, superficie)
    geo.synchronize()

    for dim, tags, nombre in ((2, [s for s, _ in parches], "Placa"),
                              (1, [re_in, re_ex], "Simetria_X"),
                              (1, [ra_in, ra_ex], "Simetria_Y"),
                              (1, [bc], "Borde_Carga"),
                              (1, agujero, "Agujero")):
        gmsh.model.setPhysicalName(dim, gmsh.model.addPhysicalGroup(dim, tags), nombre)

    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(str(ruta))
    gmsh.finalize()
    return ruta


if __name__ == "__main__":
    print(generar())
