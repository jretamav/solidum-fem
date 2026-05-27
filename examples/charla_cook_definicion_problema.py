"""Figura de definicion del problema de Cook's membrane."""
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch


HERE = os.path.dirname(__file__)
OUT = os.path.join(HERE, "charla_cook_definicion_problema.png")

# Vertices del trapecio
P1 = np.array([0.0,  0.0])    # esquina inferior izquierda
P2 = np.array([48.0, 44.0])   # esquina inferior derecha
P3 = np.array([48.0, 60.0])   # esquina superior derecha
P4 = np.array([0.0,  44.0])   # esquina superior izquierda

u_pres = 10.0   # desplazamiento prescrito vertical en el borde derecho

fig, ax = plt.subplots(figsize=(10.5, 7.5), dpi=140)

# --- Trapecio relleno + contorno ---
trapecio = Polygon([P1, P2, P3, P4], closed=True,
                   facecolor="#e8eef5", edgecolor="#1f4e79",
                   lw=2.0, zorder=2)
ax.add_patch(trapecio)

# --- Empotramiento en borde izquierdo: hatching oblicuo ---
n_hatch = 10
hatch_len = 3.0
y_hatch = np.linspace(P1[1], P4[1], n_hatch)
for yh in y_hatch:
    ax.plot([P1[0] - hatch_len, P1[0]],
            [yh + hatch_len, yh],
            color="#1f4e79", lw=1.4, zorder=4)
ax.plot([P1[0] - hatch_len, P1[0] - hatch_len],
        [P1[1] + hatch_len, P4[1] + hatch_len],
        color="#1f4e79", lw=2.0, zorder=4)

ax.text(P1[0] - 6.0, (P1[1] + P4[1]) / 2.0,
        r"$u_x = u_y = 0$",
        ha="right", va="center", fontsize=12,
        color="#1f4e79", rotation=90)

# --- Desplazamiento prescrito en borde derecho: flechas verticales ---
n_arrows = 6
y_arrows = np.linspace(P2[1] + 1.5, P3[1] - 1.5, n_arrows)
arrow_len = 8.0
for ya in y_arrows:
    ax.add_patch(FancyArrowPatch(
        (P2[0], ya), (P2[0], ya + arrow_len),
        arrowstyle="-|>", mutation_scale=14,
        color="#c00000", lw=1.8, zorder=5,
    ))

ax.text(P2[0] + 6.0, (P2[1] + P3[1]) / 2.0 + arrow_len / 2.0,
        rf"$\bar u_y = {u_pres:.0f}$",
        ha="left", va="center", fontsize=13,
        color="#c00000", fontweight="bold")
ax.text(P2[0] + 6.0, (P2[1] + P3[1]) / 2.0 + arrow_len / 2.0 - 5.0,
        r"$u_x$ libre",
        ha="left", va="center", fontsize=10, color="#c00000")

# --- Cotas ---
# Cota horizontal 48 (debajo del borde inferior)
y_cota_h = -10.0
ax.add_patch(FancyArrowPatch((P1[0], y_cota_h), (P2[0], y_cota_h),
             arrowstyle="<|-|>", mutation_scale=12,
             color="#444444", lw=1.2, zorder=4))
ax.plot([P1[0], P1[0]], [y_cota_h - 1.2, y_cota_h + 1.2], color="#444444", lw=1.2)
ax.plot([P2[0], P2[0]], [y_cota_h - 1.2, y_cota_h + 1.2], color="#444444", lw=1.2)
ax.text(24.0, y_cota_h - 3.0, "48",
        ha="center", va="top", fontsize=12, color="#444444")

# Cota altura izquierda 44 (a la izquierda del empotramiento)
x_cota_izq = -18.0
ax.add_patch(FancyArrowPatch((x_cota_izq, P1[1]), (x_cota_izq, P4[1]),
             arrowstyle="<|-|>", mutation_scale=12,
             color="#444444", lw=1.2, zorder=4))
ax.plot([x_cota_izq - 1.2, x_cota_izq + 1.2], [P1[1], P1[1]], color="#444444", lw=1.2)
ax.plot([x_cota_izq - 1.2, x_cota_izq + 1.2], [P4[1], P4[1]], color="#444444", lw=1.2)
ax.text(x_cota_izq - 1.5, (P1[1] + P4[1]) / 2.0, "44",
        ha="right", va="center", fontsize=12, color="#444444")

# Cota altura derecha 16 (a la derecha del borde cargado)
x_cota_der = 70.0
ax.add_patch(FancyArrowPatch((x_cota_der, P2[1]), (x_cota_der, P3[1]),
             arrowstyle="<|-|>", mutation_scale=12,
             color="#444444", lw=1.2, zorder=4))
ax.plot([x_cota_der - 1.2, x_cota_der + 1.2], [P2[1], P2[1]], color="#444444", lw=1.2)
ax.plot([x_cota_der - 1.2, x_cota_der + 1.2], [P3[1], P3[1]], color="#444444", lw=1.2)
ax.text(x_cota_der + 2.0, (P2[1] + P3[1]) / 2.0, "16",
        ha="left", va="center", fontsize=12, color="#444444")

# --- Punto de control: medio del borde derecho ---
mid_right = (P2 + P3) / 2.0
ax.plot(mid_right[0], mid_right[1], "o", color="#1f4e79",
        ms=8, zorder=6)
ax.annotate(r"$(48, 52)$" + "\npunto de control",
            xy=(mid_right[0], mid_right[1]),
            xytext=(mid_right[0] + 12.0, mid_right[1] - 12.0),
            fontsize=10, color="#1f4e79",
            arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.0))

# --- Ejes globales ---
ax_origin = (-25.0, -15.0)
ax.add_patch(FancyArrowPatch(ax_origin, (ax_origin[0] + 7.0, ax_origin[1]),
             arrowstyle="-|>", mutation_scale=12, color="black", lw=1.5))
ax.add_patch(FancyArrowPatch(ax_origin, (ax_origin[0], ax_origin[1] + 7.0),
             arrowstyle="-|>", mutation_scale=12, color="black", lw=1.5))
ax.text(ax_origin[0] + 8.5, ax_origin[1], "$x$", fontsize=13,
        ha="left", va="center")
ax.text(ax_origin[0], ax_origin[1] + 8.5, "$y$", fontsize=13,
        ha="center", va="bottom")

# --- Caja de propiedades ---
props_text = (
    "Material:  J2 plane strain (Von Mises 2D)\n"
    "Plasticidad perfecta  ($H = 0$)\n"
    "$E = 70$\n"
    r"$\nu = 1/3$" + "\n"
    "$\\sigma_y = 0.243$\n"
    "Espesor  $t = 1$\n"
    "\n"
    "Malla:  480 Quad4 (transfinita 24x20)\n"
    "Solver: Newton-Raphson, 100 pasos\n"
    "\n"
    "Ref:  Simo & Armero (1992)"
)
ax.text(85.0, 60.0, props_text,
        ha="left", va="top", fontsize=10.5,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#fffbe6",
                  edgecolor="#aaa", lw=1.0))

# --- Titulo ---
ax.set_title("Cook's membrane — definicion del problema",
             fontsize=14, pad=14)

# Limpieza
ax.set_aspect("equal")
ax.set_xlim(-30.0, 160.0)
ax.set_ylim(-20.0, 80.0)
ax.axis("off")

plt.tight_layout()
plt.savefig(OUT, bbox_inches="tight")
print(f"Figura guardada en: {OUT}")
