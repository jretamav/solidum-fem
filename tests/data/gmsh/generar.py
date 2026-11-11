"""Genera las mallas pequeñas de tests/test_gmsh_parser_celdas.py con Gmsh.

Uso (sólo para regenerarlas; las pruebas leen los .msh versionados y no
necesitan Gmsh): python tests/data/gmsh/generar.py

2D: placa de 2 × 1 con la línea interior «Barra» (y = 0.5) incrustada, grupos
«Placa», «Izquierda» (x = 0), «Derecha» (x = 2) y «Barra». Una variante con la
superficie en sentido horario da celdas con jacobiano negativo.
3D: bloque de 2 × 1 × 1, grupos «Bloque», «Cara_x0» y «Cara_x1». Una variante
de Tet10 con las celdas invertidas se escribe con meshio.
"""
from pathlib import Path

import gmsh
import meshio
import numpy as np

HERE = Path(__file__).resolve().parent


def _setup(name):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.model.add(name)


def placa(name, quad=False, order=1, incomplete=False, horaria=False):
    _setup(name)
    g = gmsh.model.geo
    p = [g.addPoint(x, y, 0, 0.5) for x, y in [(0, 0), (2, 0), (2, 1), (0, 1), (0, 0.5), (2, 0.5)]]
    # contorno partido en y = 0.5 para incrustar la barra de lado a lado
    l = [g.addLine(p[0], p[1]), g.addLine(p[1], p[5]), g.addLine(p[5], p[2]), g.addLine(p[2], p[3]),
         g.addLine(p[3], p[4]), g.addLine(p[4], p[0])]
    barra = g.addLine(p[4], p[5])
    loop = [l[0], l[1], l[2], l[3], l[4], l[5]]
    if horaria:
        loop = [-c for c in reversed(loop)]
    s = g.addPlaneSurface([g.addCurveLoop(loop)])
    g.synchronize()
    gmsh.model.mesh.embed(1, [barra], 2, s)
    gmsh.model.addPhysicalGroup(2, [s], name="Placa")
    gmsh.model.addPhysicalGroup(1, [l[4], l[5]], name="Izquierda")
    gmsh.model.addPhysicalGroup(1, [l[1], l[2]], name="Derecha")
    gmsh.model.addPhysicalGroup(1, [barra], name="Barra")
    if quad:
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.model.mesh.generate(2)
    if order == 2:
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1 if incomplete else 0)
        gmsh.model.mesh.setOrder(2)
    gmsh.write(str(HERE / f"{name}.msh"))
    gmsh.finalize()


def bloque(name, hexa=False, order=1, incomplete=False):
    _setup(name)
    g = gmsh.model.occ
    v = g.addBox(0, 0, 0, 2, 1, 1)
    g.synchronize()
    caras = gmsh.model.getBoundary([(3, v)], oriented=False)
    x0 = [t for d, t in caras if abs(gmsh.model.occ.getCenterOfMass(d, t)[0]) < 1e-9]
    x1 = [t for d, t in caras if abs(gmsh.model.occ.getCenterOfMass(d, t)[0] - 2) < 1e-9]
    gmsh.model.addPhysicalGroup(3, [v], name="Bloque")
    gmsh.model.addPhysicalGroup(2, x0, name="Cara_x0")
    gmsh.model.addPhysicalGroup(2, x1, name="Cara_x1")
    if hexa:
        for d, t in gmsh.model.getEntities(1):
            gmsh.model.mesh.setTransfiniteCurve(t, 3)
        for d, t in gmsh.model.getEntities(2):
            gmsh.model.mesh.setTransfiniteSurface(t)
            gmsh.model.mesh.setRecombine(2, t)
        gmsh.model.mesh.setTransfiniteVolume(v)
    else:
        gmsh.option.setNumber("Mesh.MeshSizeMin", 1.0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 1.0)
    gmsh.model.mesh.generate(3)
    if order == 2:
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1 if incomplete else 0)
        gmsh.model.mesh.setOrder(2)
    gmsh.write(str(HERE / f"{name}.msh"))
    gmsh.finalize()


def invertir_tet10():
    m = meshio.read(HERE / "tet10.msh")
    swap = [0, 2, 1, 3, 6, 5, 4, 7, 9, 8]   # ξ ↔ η en el orden de meshio
    cells = []
    for b in m.cells:
        cells.append((b.type, b.data[:, swap] if b.type == "tetra10" else b.data))
    out = meshio.Mesh(m.points, cells, cell_data=m.cell_data, field_data=m.field_data)
    meshio.write(HERE / "tet10_invertida.msh", out, file_format="gmsh22", binary=False)


if __name__ == "__main__":
    placa("tri3")
    placa("tri6", order=2)
    placa("quad4", quad=True)
    placa("quad8", quad=True, order=2, incomplete=True)
    placa("quad9", quad=True, order=2)
    placa("tri6_horaria", order=2, horaria=True)
    placa("quad8_horaria", quad=True, order=2, incomplete=True, horaria=True)
    bloque("tet4")
    bloque("tet10", order=2)
    bloque("hex8", hexa=True)
    bloque("hex20", hexa=True, order=2, incomplete=True)
    bloque("hex27", hexa=True, order=2)
    invertir_tet10()
    for f in sorted(HERE.glob("*.msh")):
        print(f.name, f.stat().st_size)
