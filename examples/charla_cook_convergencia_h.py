"""Estudio de convergencia h para Cook's membrane.

Corre tres mallas (gruesa, media, fina) con identicos parametros materiales y
de carga, registra la curva F-u de cada una, y las superpone en un solo PNG.
Demuestra que el cortante ultimo converge a un valor asintotico — sustituto
estandar de la "solucion exacta" inexistente en este problema.

Usa la API Python de gmsh para construir las mallas en un archivo temporal,
que se borra al terminar. Solo deja como artefacto el PNG final.
"""
import os
import sys
import tempfile

import gmsh
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import solidum
from solidum.math.assembly import Assembler
from solidum.utils.yaml_parser import YamlParser


HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "charla_cook_convergencia_h.png")

# (etiqueta, N_horiz, N_vert)
MALLAS = [
    ("gruesa", 12, 10),
    ("media",  24, 20),
    ("fina",   48, 40),
]

YAML_TEMPLATE = """\
mesh: {msh_path}
mesh_material: 1
mesh_thickness: 1.0

materials:
  - id: 1
    type: VonMises2D
    E: 70.0
    nu: 0.3333
    sigma_y: 0.243
    H: 0.0
    hypothesis: plane_strain
    density: 1.0

boundary_conditions_by_group:
  - {{group_name: "Apoyo_Izquierdo", ux: 0.0, uy: 0.0}}
  - {{group_name: "Carga_Derecha",   uy: 10.0}}

solver:
  type: NonlinearSolver
  num_steps: 100
  max_iter: 25
  adaptive: true
  convergence:
    rtol_force: 1.0e-6
    rtol_disp:  1.0e-6
"""

X_LEFT = 0.0
X_RIGHT = 48.0
Y_MID = 52.0
TOL = 1e-6


def generar_msh_cook(msh_path: str, n_horiz: int, n_vert: int) -> None:
    """Construye y guarda una malla transfinita del trapecio de Cook."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("cook")

    p1 = gmsh.model.geo.addPoint( 0.0,  0.0, 0)
    p2 = gmsh.model.geo.addPoint(48.0, 44.0, 0)
    p3 = gmsh.model.geo.addPoint(48.0, 60.0, 0)
    p4 = gmsh.model.geo.addPoint( 0.0, 44.0, 0)

    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)

    cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([cl])

    gmsh.model.geo.mesh.setTransfiniteCurve(l1, n_horiz + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l3, n_horiz + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l2, n_vert + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l4, n_vert + 1)
    gmsh.model.geo.mesh.setTransfiniteSurface(surf)
    gmsh.model.geo.mesh.setRecombine(2, surf)

    gmsh.model.geo.synchronize()

    gmsh.model.addPhysicalGroup(1, [l4], tag=10, name="Apoyo_Izquierdo")
    gmsh.model.addPhysicalGroup(1, [l2], tag=11, name="Carga_Derecha")
    gmsh.model.addPhysicalGroup(2, [surf], tag=12, name="Dominio_Cook")

    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(msh_path)
    gmsh.finalize()


resultados = {}

with tempfile.TemporaryDirectory() as tmp:
    for nombre, nh, nv in MALLAS:
        n_elem = nh * nv
        print(f"\n=== Malla {nombre}: {nh}x{nv} = {n_elem} Quad4 ===")

        msh_path = os.path.join(tmp, f"cook_{nombre}.msh")
        yml_path = os.path.join(tmp, f"cook_{nombre}.yaml")

        generar_msh_cook(msh_path, nh, nv)

        with open(yml_path, "w") as f:
            f.write(YAML_TEMPLATE.format(msh_path=msh_path.replace("\\", "/")))

        parser = YamlParser(yml_path)
        domain = parser.parse()
        domain.generate_equation_numbers()
        assembler = Assembler(domain)
        solver = parser.get_solver(assembler)
        F_ext = parser.get_external_forces() + parser.get_body_load(assembler)

        left_uy = []
        mid_node = None
        mid_dist = float("inf")
        for node in domain.nodes.values():
            x, y = node.coordinates[0], node.coordinates[1]
            if "uy" not in node.dofs:
                continue
            if abs(x - X_LEFT) < TOL:
                left_uy.append(node.dofs["uy"])
            elif abs(x - X_RIGHT) < TOL:
                d = abs(y - Y_MID)
                if d < mid_dist:
                    mid_dist = d
                    mid_node = node

        mid_uy = mid_node.dofs["uy"]
        history = []

        def cb(step, U, lam,
               _h=history, _m=mid_uy, _l=left_uy, _a=assembler):
            _, Fi = _a.assemble_non_linear_system(U)
            _h.append((float(U[_m]), -float(np.sum(Fi[_l]))))

        solidum.run(
            domain, assembler=assembler, solver=solver,
            F_applied=F_ext, step_callback=cb,
        )

        resultados[nombre] = {
            "n_elem": n_elem,
            "u_mid": np.array([h[0] for h in history]),
            "V":     np.array([h[1] for h in history]),
            "V_final": history[-1][1],
        }
        print(f"  Cortante final: V = {resultados[nombre]['V_final']:.4f}")

print("\n--- Resumen de convergencia ---")
print(f"{'Malla':10s} {'Elementos':>10s} {'V_final':>10s}")
for nombre, _, _ in MALLAS:
    r = resultados[nombre]
    print(f"{nombre:10s} {r['n_elem']:>10d} {r['V_final']:>10.4f}")

fig, ax = plt.subplots(figsize=(7.8, 5.0), dpi=130)

colors = ["#9bc2e6", "#1f4e79", "#003366"]
styles = ["o--", "s-", "^-"]

for (nombre, _, _), color, style in zip(MALLAS, colors, styles):
    r = resultados[nombre]
    ax.plot(r["u_mid"], r["V"], style, lw=1.6, ms=4, color=color,
            label=f"{nombre.capitalize()} ({r['n_elem']} Quad4) — "
                  f"$V_{{\\mathrm{{final}}}}={r['V_final']:.3f}$")

ax.set_xlabel(r"Desplazamiento vertical del punto medio  $u_y^{(48,52)}$")
ax.set_ylabel(r"Cortante resultante en el apoyo  $V$")
ax.set_title("Cook's membrane — convergencia h\n"
             r"J2 plane strain, plasticidad perfecta ($H=0$),  "
             r"$E=70$, $\nu=1/3$, $\sigma_y=0.243$")
ax.grid(True, alpha=0.3)
ax.legend(loc="lower right", fontsize=9)
plt.tight_layout()
plt.savefig(OUT, bbox_inches="tight")
print(f"\nGrafica guardada en: {OUT}")
