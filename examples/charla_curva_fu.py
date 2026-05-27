"""Curva fuerza-desplazamiento para la charla.

Re-ejecuta `charla_placa_perfecta.yaml` con un callback que en cada paso
convergido suma las fuerzas internas en los DOFs ux del borde izquierdo
(reaccion en el apoyo) y del borde derecho (traccion en la cara cargada).

Genera un PNG con la curva F-u.
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
YAML = os.path.join(HERE, "charla_placa_perfecta.yaml")
TOL = 1e-6
L = 2.0          # lado de la placa
THICKNESS = 1.0  # espesor unitario declarado en el YAML

parser = YamlParser(YAML)
domain = parser.parse()
domain.generate_equation_numbers()
assembler = Assembler(domain)
solver = parser.get_solver(assembler)
F_ext = parser.get_external_forces() + parser.get_body_load(assembler)

# Identificacion de DOFs por coordenada (sin depender de los physical groups).
right_dofs = []
left_dofs = []
prescribed_ux = None
for node in domain.nodes.values():
    x = node.coordinates[0]
    if "ux" not in node.dofs:
        continue
    if abs(x - L / 2) < TOL:
        right_dofs.append(node.dofs["ux"])
        if prescribed_ux is None:
            prescribed_ux = node.boundary_conditions.get("ux", 0.0)
    elif abs(x + L / 2) < TOL:
        left_dofs.append(node.dofs["ux"])

assert right_dofs and left_dofs, "No se identificaron los nodos del borde derecho/izquierdo"
print(f"Nodos en borde derecho (carga):  {len(right_dofs)}")
print(f"Nodos en borde izquierdo (apoyo): {len(left_dofs)}")
print(f"Desplazamiento prescrito final:   {prescribed_ux*1e3:.3f} mm\n")

history = []  # (u_aplicado [m], R_izq [N/m de espesor], R_der [N/m])

def callback(step, U, load_factor):
    # assemble_non_linear_system es idempotente sobre el estado ya commited:
    # recalcula F_int para la U dada, asi obtenemos las reacciones del paso.
    _, F_int = assembler.assemble_non_linear_system(U)
    R_left = float(np.sum(F_int[left_dofs]))
    R_right = float(np.sum(F_int[right_dofs]))
    u_now = prescribed_ux * load_factor
    history.append((u_now, R_left, R_right))
    print(f"  paso {step:2d}  lambda={load_factor:.3f}  u={u_now*1e3:.4f} mm  "
          f"F_apoyo={R_left/1e3:+.2f} kN/m  F_carga={R_right/1e3:+.2f} kN/m")

print("Re-ejecutando con callback para registrar la curva F-u...\n")
solidum.run(domain, assembler=assembler, solver=solver,
            F_applied=F_ext, step_callback=callback)

u_arr = np.array([h[0] for h in history]) * 1e3       # mm
F_right = np.array([h[2] for h in history]) / 1e3     # kN/m de espesor
F_left = np.array([h[1] for h in history]) / 1e3

# Diagnostico de equilibrio (debe ser ~0)
imbalance = F_left[-1] + F_right[-1]
print(f"\nEquilibrio en x al final: R_izq + R_der = {imbalance:.3e} kN/m  (~0)")

# Limite teorico de colapso plastico (J2 plane strain, criterio de tubo):
#   sigma_y_pl_strain = sigma_y * 2/sqrt(3) = 288.7 MPa  (sin agujero)
#   con agujero, el area resistente neta es (L - 2r)*t = (2 - 0.8)*1 = 1.2 m
#   carga limite "neta": ~ 288.7 MPa * 1.2 m * 1 m = 346.4 kN/m  (cota superior)
sigma_y = 250.0e6
collapse_ref = (sigma_y * 2.0 / np.sqrt(3.0)) * (L - 0.8) * THICKNESS / 1e3  # kN/m
print(f"Cota de colapso teorica (net-section, J2 plane strain): {collapse_ref:.1f} kN/m")

fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=130)
ax.plot(u_arr, np.abs(F_right), "o-", lw=2.0, ms=5, color="#1f4e79",
        label="Solidum FEM (601 Quad4)")
ax.axhline(collapse_ref, ls="--", color="#c00000", lw=1.3,
           label=f"Cota neta J2 plane strain ({collapse_ref:.0f} kN/m)")
ax.set_xlabel(r"Desplazamiento prescrito  $u_x$  (mm)")
ax.set_ylabel(r"Reaccion resultante en x / espesor  (kN/m)")
ax.set_title("Placa con agujero — J2 plane strain, plasticidad perfecta ($H=0$)\n"
             r"$E=200$ GPa,  $\nu=0.3$,  $\sigma_y=250$ MPa")
ax.grid(True, alpha=0.3)
ax.legend(loc="lower right")
plt.tight_layout()

out_png = os.path.join(HERE, "charla_placa_perfecta_curva_Fu.png")
plt.savefig(out_png)
print(f"\nCurva F-u guardada en: {out_png}")
