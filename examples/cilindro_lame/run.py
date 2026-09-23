"""Ejemplo 3 del manual de ejemplos — cilindro de pared gruesa (Lamé).

Estudio de convergencia con el API de Python: un cuarto de corona circular
bajo presión interna, en deformación plana, resuelto con cuatro elementos
(Tri3, Quad4, Tri6, Quad8) en mallas cada vez más finas. Se mide el error en
desplazamientos y en esfuerzos contra la solución de Lamé y se estima el orden
de convergencia de cada elemento.

Uso::

    python examples/cilindro_lame/run.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))

from ejemplos_comun import (  # noqa: E402
    COLOR_ANALITICO, COLOR_FE, COLOR_SECUNDARIO, COLOR_TERCIARIO,
    estilo_figuras, guardar_figura, guardar_resultados,
)

from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.solid_2d import Quad4, Quad8, Tri3, Tri6  # noqa: E402
from solidum.materials.elastic_2d import Elastic2D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.batch.postprocess import gauss_states  # noqa: E402
from solidum.math.solvers import LinearSolver  # noqa: E402

RI, RE = 0.1, 0.2          # radios interior y exterior [m]
P = 10.0e6                 # presión interna [Pa]
E, NU = 200.0e9, 0.3       # acero, deformación plana

# Mallas de n × n celdas (radial × circunferencial). Los cuadráticos llegan a
# la precisión de máquina antes, así que no necesitan las más finas.
MALLAS = {"Tri3": (2, 4, 8, 16, 32, 64), "Quad4": (2, 4, 8, 16, 32, 64),
          "Tri6": (2, 4, 8, 16, 32), "Quad8": (2, 4, 8, 16, 32)}
CLASES = {"Tri3": Tri3, "Quad4": Quad4, "Tri6": Tri6, "Quad8": Quad8}

# Órdenes teóricos (Strang y Fix; Hughes §4): con interpolación de grado p,
# el desplazamiento converge como h^(p+1) y el esfuerzo como h^p. Lo que el
# capítulo afirma: el esfuerzo converge con el orden teórico (±0.15) y el
# desplazamiento al menos con él (déficit ≤ 0.15). En desplazamientos se
# permite superarlo: en mallas regulares los valores nodales superconvergen
# (medido 2026-09-23: Quad8 da orden ≈ 4, no 3).
TOLERANCIAS = {
    **{f"{e}.deficit_orden_u": 0.15 for e in ("Tri3", "Quad4", "Tri6", "Quad8")},
    **{f"{e}.error_orden_s": 0.15 for e in ("Tri3", "Quad4", "Tri6", "Quad8")},
}
ORDEN_TEORICO = {"Tri3": (2, 1), "Quad4": (2, 1), "Tri6": (3, 2), "Quad8": (3, 2)}


def lame(r):
    """σrr, σθθ y u_r de Lamé en deformación plana."""
    r = np.asarray(r, dtype=float)
    a = P * RI**2 / (RE**2 - RI**2)
    b = P * RI**2 * RE**2 / (RE**2 - RI**2)
    ur = (1 + NU) / E * ((1 - 2 * NU) * a * r + b / r)
    return a - b / r**2, a + b / r**2, ur


def modelo(nombre: str, n: int, material=None, cuadratura: str | None = None):
    """Cuarto de corona con n × n celdas de ``nombre``. Devuelve el dominio,
    las aristas cargadas ``[(elemento, arista)]`` y los nodos. Por omisión,
    acero elástico en deformación plana y la cuadratura propia de cada
    elemento (el ejemplo 4 pasa un material plástico y, para el Quad8,
    integración reducida ``"2x2"``)."""
    cls = CLASES[nombre]
    cuadratico = nombre in ("Tri6", "Quad8")
    k = 2 if cuadratico else 1          # nodos por celda en cada dirección
    ni = nj = k * n + 1
    if material is None:
        material = Elastic2D(E=E, nu=NU, hypothesis="plane_strain")
    domain = Domain()
    nodo = {}
    for j in range(nj):
        th = 0.5 * np.pi * j / (nj - 1)
        for i in range(ni):
            if nombre == "Quad8" and i % 2 == 1 and j % 2 == 1:
                continue                    # Quad8 no tiene nodo central
            r = RI + (RE - RI) * i / (ni - 1)
            nodo[i, j] = domain.add_node(len(nodo) + 1, [r * np.cos(th), r * np.sin(th)])

    aristas = []
    for cj in range(n):
        for ci in range(n):
            i, j = k * ci, k * cj
            c1, c2, c3, c4 = nodo[i, j], nodo[i + k, j], nodo[i + k, j + k], nodo[i, j + k]
            eid = len(domain.elements) + 1
            if nombre == "Quad4":
                el = [cls(eid, [c1, c2, c3, c4], material, thickness=1.0,
                          quadrature=cuadratura)]
                cara = (el[0], 3)                          # arista c4 → c1, en r = RI
            elif nombre == "Quad8":
                m = [nodo[i + 1, j], nodo[i + 2, j + 1], nodo[i + 1, j + 2], nodo[i, j + 1]]
                el = [cls(eid, [c1, c2, c3, c4] + m, material, thickness=1.0,
                          quadrature=cuadratura)]
                cara = (el[0], 3)
            elif nombre == "Tri3":
                el = [cls(eid, [c1, c2, c3], material, thickness=1.0),
                      cls(eid + 1, [c1, c3, c4], material, thickness=1.0)]
                cara = (el[1], 2)                          # arista c4 → c1
            else:                                          # Tri6
                cen = nodo[i + 1, j + 1]
                el = [cls(eid, [c1, c2, c3, nodo[i + 1, j], nodo[i + 2, j + 1], cen],
                          material, thickness=1.0),
                      cls(eid + 1, [c1, c3, c4, cen, nodo[i + 1, j + 2], nodo[i, j + 1]],
                          material, thickness=1.0)]
                cara = (el[1], 2)
            for e in el:
                domain.add_element(e)
            if ci == 0:
                aristas.append(cara)

    for i in range(ni):                    # simetría
        if (i, 0) in nodo:
            nodo[i, 0].fix_dof("uy", 0.0)
        if (i, nj - 1) in nodo:
            nodo[i, nj - 1].fix_dof("ux", 0.0)
    domain.generate_equation_numbers()
    return domain, aristas, nodo


def carga_de_presion(domain, aristas, p: float) -> np.ndarray:
    """Vector de cargas consistente con una presión ``p`` sobre las aristas.

    La presión empuja contra la cara: t = −p n, con n la normal exterior del
    sólido. Se integra ∫ Nᵀ t dΓ sobre la geometría discretizada de cada
    arista (recta o parabólica), con tres puntos de Gauss: n varía a lo largo
    de una arista curva, así que no sirve una tracción constante.
    """
    F = np.zeros(domain.total_dofs)
    xg, wg = np.polynomial.legendre.leggauss(3)
    for el, arista in aristas:
        loc = el.EDGE_NODES[arista]
        X = np.array([el.nodes[q].coordinates[:2] for q in loc], dtype=float)
        dofs = el.get_global_dof_indices()
        for s, w in zip(xg, wg):
            if len(loc) == 2:
                N = np.array([(1 - s) / 2, (1 + s) / 2])
                dN = np.array([-0.5, 0.5])
            else:                                          # (extremo, medio, extremo)
                N = np.array([s * (s - 1) / 2, 1 - s**2, s * (s + 1) / 2])
                dN = np.array([s - 0.5, -2 * s, s + 0.5])
            dx, dy = dN @ X
            # Recorriendo la arista en el sentido del elemento (antihorario), la
            # normal exterior multiplicada por |J| es (dy, −dx).
            t = -p * np.array([dy, -dx])
            for a, q in enumerate(loc):
                for d in range(2):
                    g = dofs[2 * q + d]
                    if g >= 0:
                        F[g] += w * N[a] * t[d] * el.thickness
    return F


def resolver(nombre: str, n: int) -> dict:
    domain, aristas, nodo = modelo(nombre, n)
    F = carga_de_presion(domain, aristas, P)
    U = LinearSolver(Assembler(domain)).solve(F)

    # Desplazamiento radial en todos los nodos.
    xy = np.array([nd.coordinates for nd in nodo.values()])
    ux = np.array([U[nd.dofs["ux"]] for nd in nodo.values()])
    uy = np.array([U[nd.dofs["uy"]] for nd in nodo.values()])
    r = np.hypot(xy[:, 0], xy[:, 1])
    ur = (ux * xy[:, 0] + uy * xy[:, 1]) / r
    ur_ref = lame(r)[2]
    e_u = float(np.sqrt(np.mean((ur - ur_ref) ** 2) / np.mean(ur_ref**2)))

    # Esfuerzos en los puntos de Gauss, en polares.
    est = gauss_states(domain, U)
    pg = np.vstack([g["points_global"] for g in est.values()])
    sg = np.vstack([g["stress"] for g in est.values()])
    rg = np.hypot(pg[:, 0], pg[:, 1])
    c, s = pg[:, 0] / rg, pg[:, 1] / rg
    srr = sg[:, 0] * c**2 + sg[:, 1] * s**2 + 2 * sg[:, 2] * c * s
    stt = sg[:, 0] * s**2 + sg[:, 1] * c**2 - 2 * sg[:, 2] * c * s
    srr_ref, stt_ref, _ = lame(rg)
    e_s = float(np.sqrt(np.mean((srr - srr_ref) ** 2 + (stt - stt_ref) ** 2)
                        / np.mean(srr_ref**2 + stt_ref**2)))
    return {"n": n, "h": (RE - RI) / n, "gdl": int(domain.total_dofs),
            "error_u": e_u, "error_s": e_s,
            "r_gauss": rg, "srr": srr, "stt": stt}


def _orden(h, e) -> float:
    """Pendiente log-log entre las dos mallas más finas."""
    return float(np.log(e[-2] / e[-1]) / np.log(h[-2] / h[-1]))


def calcular() -> dict:
    out = {"Ri": RI, "Re": RE, "p": P, "E": E, "nu": NU}
    _, stt_a, ur_a = lame(RI)
    out["sigma_tt_Ri_MPa"] = stt_a / 1e6
    out["u_r_Ri_um"] = ur_a * 1e6
    for nombre, mallas in MALLAS.items():
        filas = [resolver(nombre, n) for n in mallas]
        h = np.array([f["h"] for f in filas])
        eu = np.array([f["error_u"] for f in filas])
        es = np.array([f["error_s"] for f in filas])
        pu, ps = ORDEN_TEORICO[nombre]
        orden_u, orden_s = _orden(h, eu), _orden(h, es)
        out[nombre] = {
            "n": list(mallas), "h": h, "gdl": [f["gdl"] for f in filas],
            "error_u": eu, "error_s": es,
            "orden_u": orden_u, "orden_s": orden_s,
            "orden_u_teorico": pu, "orden_s_teorico": ps,
            "deficit_orden_u": max(0.0, pu - orden_u), "error_orden_s": abs(orden_s - ps),
            "error_s_n8_pct": 100 * es[list(mallas).index(8)],
            "error_u_n8_pct": 100 * eu[list(mallas).index(8)],
            "gdl_n8": filas[list(mallas).index(8)]["gdl"],
            "n_final": mallas[-1], "gdl_final": filas[-1]["gdl"],
            "error_s_final_pct": 100 * es[-1],
        }
    # Perfiles radiales para la figura (malla de 4 × 4 celdas).
    for nombre in ("Quad4", "Quad8"):
        f = resolver(nombre, 4)
        out[f"perfil_{nombre}"] = {"r": f["r_gauss"], "srr": f["srr"], "stt": f["stt"]}
    return out


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    # --- Esquema: malla de 6 × 6 Quad4, presión interna y apoyos de simetría
    domain, aristas, _ = modelo("Quad4", 6)
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    polys = [np.array([nd.coordinates for nd in el.nodes]) for el in domain.elements.values()]
    ax.add_collection(PolyCollection(polys, facecolor="#e8eef6", edgecolor="0.4", lw=0.5))
    for th in np.linspace(0.05, 0.5 * np.pi - 0.05, 9):
        d = np.array([np.cos(th), np.sin(th)])
        ax.annotate("", xy=RI * d, xytext=(RI - 0.035) * d,
                    arrowprops=dict(arrowstyle="-|>", color=COLOR_SECUNDARIO, lw=1.2))
    ax.text(0.035, 0.035, r"$p$", color=COLOR_SECUNDARIO, fontsize=12)
    for x in np.linspace(RI, RE, 5):
        ax.plot(x, -0.006, "^", color=COLOR_FE, ms=6)
        ax.plot(-0.006, x, ">", color=COLOR_FE, ms=6)
    ax.text(0.15, -0.022, r"$u_y = 0$", ha="center", color=COLOR_FE)
    ax.text(-0.024, 0.15, r"$u_x = 0$", rotation=90, va="center", color=COLOR_FE)
    ax.text(0.055, 0.005, r"$R_i$", fontsize=9)
    ax.text(0.19, 0.12, r"$R_e$", fontsize=9)
    ax.set_xlim(-0.035, 0.215)
    ax.set_ylim(-0.035, 0.215)
    ax.set_aspect("equal")
    ax.axis("off")
    guardar_figura(fig, carpeta, "fig_esquema")

    # --- Convergencia
    estilos = {"Tri3": ("v", COLOR_TERCIARIO, "--"), "Quad4": ("s", COLOR_FE, "--"),
               "Tri6": ("^", COLOR_SECUNDARIO, "-"), "Quad8": ("o", COLOR_ANALITICO, "-")}
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.5), sharex=True)
    for nombre, (mk, col, ls) in estilos.items():
        d = res[nombre]
        h = np.asarray(d["h"]) * 1e3
        axs[0].loglog(h, d["error_u"], mk + ls, color=col, ms=4.5, lw=1,
                      label=f"{nombre} ({d['orden_u']:.2f})")
        axs[1].loglog(h, d["error_s"], mk + ls, color=col, ms=4.5, lw=1,
                      label=f"{nombre} ({d['orden_s']:.2f})")
    axs[0].set_ylabel(r"error relativo en $u_r$")
    axs[1].set_ylabel(r"error relativo en $\sigma$")
    for ax in axs:
        ax.set_xlabel(r"tamaño radial del elemento $h$ [mm]")
        ax.legend(title="elemento (orden medido)", fontsize=8, title_fontsize=8)
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_convergencia")

    # --- Perfiles radiales con la malla de 4 × 4
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    rr = np.linspace(RI, RE, 200)
    srr, stt, _ = lame(rr)
    ax.plot(rr * 1e3, srr / 1e6, color=COLOR_ANALITICO, lw=1.3, label="Lamé")
    ax.plot(rr * 1e3, stt / 1e6, color=COLOR_ANALITICO, lw=1.3)
    for nombre, mk, col in (("Quad4", "s", COLOR_FE), ("Quad8", "o", COLOR_SECUNDARIO)):
        d = res[f"perfil_{nombre}"]
        r = np.asarray(d["r"]) * 1e3
        ax.plot(r, np.asarray(d["srr"]) / 1e6, mk, color=col, ms=3.5, mfc="none",
                label=f"{nombre} 4×4")
        ax.plot(r, np.asarray(d["stt"]) / 1e6, mk, color=col, ms=3.5, mfc="none")
    ax.text(150, 11.0, r"$\sigma_{\theta\theta}$")
    ax.text(150, -4.6, r"$\sigma_{rr}$")
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set_xlabel(r"$r$ [mm]")
    ax.set_ylabel("esfuerzo [MPa]")
    ax.legend(loc="center right")
    guardar_figura(fig, carpeta, "fig_perfiles")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    for nombre in MALLAS:
        d = res[nombre]
        print(f"{nombre:6s} orden u = {d['orden_u']:.2f} (teórico {d['orden_u_teorico']})"
              f"   orden σ = {d['orden_s']:.2f} (teórico {d['orden_s_teorico']})"
              f"   error σ 8×8 = {d['error_s_n8_pct']:.2f} %")
    return res


if __name__ == "__main__":
    main()
