"""Curva fuerza-desplazamiento para Cook's membrane.

Eje horizontal: desplazamiento vertical del **punto medio del borde derecho**
(referencia clasica en literatura: Simo-Armero 1992, Wriggers, Bathe).
Eje vertical: cortante total transmitido por el apoyo izquierdo (= resultante
vertical de las reacciones).
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import solidum
from solidum.math.assembly import Assembler
from solidum.utils.yaml_parser import YamlParser


HERE = os.path.dirname(__file__)
YAML = os.path.join(HERE, "charla_cook.yaml")
TOL = 1e-6

# Geometria
X_LEFT = 0.0
X_RIGHT = 48.0
Y_MID_RIGHT = 52.0   # punto medio del borde derecho: (48, (44+60)/2 = 52)

parser = YamlParser(YAML)
domain = parser.parse()
domain.generate_equation_numbers()
assembler = Assembler(domain)
solver = parser.get_solver(assembler)
F_ext = parser.get_external_forces() + parser.get_body_load(assembler)

# DOFs uy del borde izquierdo (apoyo) y DOFs uy del borde derecho (carga).
# El punto medio del borde derecho lo identificamos por proximidad a y=52.
left_dofs_uy = []
right_dofs_uy = []
mid_right_node = None
mid_right_dist = float("inf")
prescribed_uy = None

for node in domain.nodes.values():
    x, y = node.coordinates[0], node.coordinates[1]
    if "uy" not in node.dofs:
        continue
    if abs(x - X_LEFT) < TOL:
        left_dofs_uy.append(node.dofs["uy"])
    elif abs(x - X_RIGHT) < TOL:
        right_dofs_uy.append(node.dofs["uy"])
        if prescribed_uy is None:
            prescribed_uy = node.boundary_conditions.get("uy", 0.0)
        # candidato a punto medio del borde derecho
        d = abs(y - Y_MID_RIGHT)
        if d < mid_right_dist:
            mid_right_dist = d
            mid_right_node = node

assert left_dofs_uy and right_dofs_uy, "No se identificaron los bordes"
assert mid_right_node is not None, "No se encontro el punto medio del borde derecho"

mid_uy_dof = mid_right_node.dofs["uy"]
print(f"Nodos en borde izquierdo (apoyo): {len(left_dofs_uy)}")
print(f"Nodos en borde derecho  (carga): {len(right_dofs_uy)}")
print(f"Punto medio borde derecho: nodo en y={mid_right_node.coordinates[1]:.3f} "
      f"(referencia y=52)")
print(f"Desplazamiento vertical prescrito final: {prescribed_uy:.3f}\n")

history = []   # (u_mid_aplicado, F_apoyo, F_carga)

def callback(step, U, load_factor):
    _, F_int = assembler.assemble_non_linear_system(U)
    R_left = float(np.sum(F_int[left_dofs_uy]))
    R_right = float(np.sum(F_int[right_dofs_uy]))
    u_mid = float(U[mid_uy_dof])
    history.append((u_mid, R_left, R_right))
    print(f"  paso {step:3d}  lambda={load_factor:.3f}  "
          f"u_mid={u_mid:.4f}  F_apoyo={R_left:+.4f}  F_carga={R_right:+.4f}")

print("Re-ejecutando para registrar la curva F-u...\n")
solidum.run(domain, assembler=assembler, solver=solver,
            F_applied=F_ext, step_callback=callback)

u_arr = np.array([h[0] for h in history])
F_apoyo = np.array([h[1] for h in history])
F_carga = np.array([h[2] for h in history])

imbalance = F_apoyo[-1] + F_carga[-1]
print(f"\nEquilibrio en y al final: R_izq + R_der = {imbalance:.3e}  (~0)")

fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=130)
ax.plot(u_arr, np.abs(F_carga), "o-", lw=2.0, ms=4, color="#1f4e79",
        label="Solidum FEM (480 Quad4)")
ax.set_xlabel(r"Desplazamiento vertical del punto medio  $u_y^{(48,52)}$")
ax.set_ylabel(r"Cortante resultante en el apoyo  $V$")
ax.set_title("Cook's membrane — J2 plane strain, plasticidad perfecta ($H=0$)\n"
             r"$E=70$,  $\nu=1/3$,  $\sigma_y=0.243$")
ax.grid(True, alpha=0.3)
ax.legend(loc="lower right")
plt.tight_layout()

out_png = os.path.join(HERE, "charla_cook_curva_Fu.png")
plt.savefig(out_png)
print(f"\nCurva F-u guardada en: {out_png}")
