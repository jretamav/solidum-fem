"""Figura de definicion del problema para la charla.

Esquema ingenieril: placa cuadrada con agujero central, empotramiento en el
borde izquierdo, desplazamiento prescrito en el borde derecho, cotas y datos
del material.
"""
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch
from matplotlib.lines import Line2D


HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "charla_definicion_problema.png")

# Parametros del problema
L = 2.0          # lado de la placa (m)
r = 0.4          # radio del agujero (m)
u_pres = 4.0e-3  # desplazamiento prescrito (m)

fig, ax = plt.subplots(figsize=(10.5, 7.5), dpi=140)

# --- Placa con agujero (relleno claro + contorno) -------------------------
plate = Rectangle((-L/2, -L/2), L, L, facecolor="#e8eef5",
                  edgecolor="#1f4e79", lw=2.0, zorder=2)
ax.add_patch(plate)
hole = Circle((0, 0), r, facecolor="white", edgecolor="#1f4e79",
              lw=2.0, zorder=3)
ax.add_patch(hole)

# --- Empotramiento en borde izquierdo: hatching oblicuo -------------------
n_hatch = 12
hatch_len = 0.13
y_hatch = np.linspace(-L/2, L/2, n_hatch)
for yh in y_hatch:
    ax.plot([-L/2 - hatch_len, -L/2],
            [yh + hatch_len, yh],
            color="#1f4e79", lw=1.4, zorder=4)
# Linea solida que une los inicios del hatching (pared del apoyo)
ax.plot([-L/2 - hatch_len, -L/2 - hatch_len],
        [-L/2 + hatch_len, L/2 + hatch_len],
        color="#1f4e79", lw=2.0, zorder=4)

ax.text(-L/2 - 0.22, 0.0, r"$u_x = u_y = 0$",
        ha="right", va="center", fontsize=12, color="#1f4e79",
        rotation=90)

# --- Desplazamiento prescrito en borde derecho: flechas hacia +x ---------
n_arrows = 7
y_arrows = np.linspace(-L/2 + 0.15, L/2 - 0.15, n_arrows)
arrow_len = 0.32
for ya in y_arrows:
    ax.add_patch(FancyArrowPatch(
        (L/2, ya), (L/2 + arrow_len, ya),
        arrowstyle="-|>", mutation_scale=14,
        color="#c00000", lw=1.8, zorder=5,
    ))

ax.text(L/2 + arrow_len + 0.15, 0.0,
        rf"$\bar u_x = {u_pres*1e3:.0f}\,$mm",
        ha="left", va="center", fontsize=13, color="#c00000",
        fontweight="bold")
ax.text(L/2 + arrow_len + 0.15, -0.22, r"$u_y$ libre",
        ha="left", va="center", fontsize=10, color="#c00000")

# --- Cotas L y r ---------------------------------------------------------
# Cota L horizontal (debajo de la placa)
y_cota = -L/2 - 0.42
ax.add_patch(FancyArrowPatch((-L/2, y_cota), (L/2, y_cota),
             arrowstyle="<|-|>", mutation_scale=12,
             color="#444444", lw=1.2, zorder=4))
ax.plot([-L/2, -L/2], [y_cota - 0.05, y_cota + 0.05], color="#444444", lw=1.2)
ax.plot([L/2, L/2], [y_cota - 0.05, y_cota + 0.05], color="#444444", lw=1.2)
ax.text(0, y_cota - 0.13, rf"$L = {L:.1f}$ m",
        ha="center", va="top", fontsize=12, color="#444444")

# Cota r (radio del agujero)
ax.add_patch(FancyArrowPatch((0, 0), (r * np.cos(np.deg2rad(35)),
                                       r * np.sin(np.deg2rad(35))),
             arrowstyle="-|>", mutation_scale=11,
             color="#444444", lw=1.2, zorder=6))
ax.text(r * np.cos(np.deg2rad(35)) * 0.55,
        r * np.sin(np.deg2rad(35)) * 0.55 + 0.06,
        rf"$r = {r:.1f}$ m",
        ha="center", va="bottom", fontsize=11, color="#444444")

# --- Ejes x, y ------------------------------------------------------------
ax_origin = (-L/2 - 0.55, -L/2 - 0.85)
ax.add_patch(FancyArrowPatch(ax_origin, (ax_origin[0] + 0.35, ax_origin[1]),
             arrowstyle="-|>", mutation_scale=12, color="black", lw=1.5))
ax.add_patch(FancyArrowPatch(ax_origin, (ax_origin[0], ax_origin[1] + 0.35),
             arrowstyle="-|>", mutation_scale=12, color="black", lw=1.5))
ax.text(ax_origin[0] + 0.40, ax_origin[1], "$x$",
        ha="left", va="center", fontsize=13)
ax.text(ax_origin[0], ax_origin[1] + 0.40, "$y$",
        ha="center", va="bottom", fontsize=13)

# --- Caja de propiedades del material ------------------------------------
props_text = (
    "Material:  J2 plane strain (Von Mises 2D)\n"
    "Plasticidad perfecta  ($H = 0$)\n"
    "$E = 200$ GPa\n"
    r"$\nu = 0.3$" + "\n"
    "$\\sigma_y = 250$ MPa\n"
    "Espesor  $t = 1.0$ m\n"
    "\n"
    "Malla:  601 Quad4 (gmsh)\n"
    "Solver: Newton-Raphson, 20 pasos"
)
ax.text(L/2 + 1.35, L/2 + 0.05, props_text,
        ha="left", va="top", fontsize=10.5,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#fffbe6",
                  edgecolor="#aaa", lw=1.0))

# --- Titulo ---------------------------------------------------------------
ax.set_title("Placa cuadrada con agujero central en traccion uniaxial\n"
             "Definicion del problema",
             fontsize=14, pad=14)

# --- Limpieza de ejes y proporcion ---------------------------------------
ax.set_aspect("equal")
ax.set_xlim(-L/2 - 1.1, L/2 + 3.6)
ax.set_ylim(-L/2 - 1.3, L/2 + 0.7)
ax.axis("off")

plt.tight_layout()
plt.savefig(OUT, bbox_inches="tight")
print(f"Figura de definicion guardada en: {OUT}")
