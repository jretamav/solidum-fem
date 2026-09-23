"""Ejemplo 2 del manual de ejemplos — placa con agujero circular (Kirsch).

Cuarto de placa por doble simetría, malla estructurada de Quad4 generada con
Gmsh (``malla.py``), carga por desplazamiento impuesto. Se compara con la
solución de Kirsch (1898) para un agujero en una placa infinita sometida a
tracción uniaxial σ∞:

- esfuerzos en todos los puntos de Gauss contra el campo de Kirsch evaluado en
  el mismo punto;
- esfuerzo circunferencial en el borde del agujero, σθθ = σ∞ (1 − 2 cos 2θ);
- factor de concentración de esfuerzos K_t = σxx(0, a) / σ∞ = 3.

Uso::

    python examples/placa_agujero_kirsch/run.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))

from ejemplos_comun import (  # noqa: E402
    COLOR_ANALITICO, COLOR_FE, COLOR_SECUNDARIO, error_relativo,
    esfuerzos_nodales, estilo_figuras, guardar_figura, guardar_resultados,
    resolver_yaml_estatico,
)

from solidum.math.batch.postprocess import gauss_states  # noqa: E402

# Datos (deben coincidir con malla.py y modelo.yaml).
A = 0.01
W = H = 0.2
ESPESOR = 0.01
R_CERCANO = 2 * A      # zona de gradiente fuerte junto al agujero
R_LEJANO = 5 * A       # más allá, el error lo domina el modelo, no la malla

# Cotas de lo que el capítulo afirma (errores relativos a σ∞ salvo K_t).
# Estudio de refinamiento (2026-09-23, mallas con 1/2, 1, 2 y 4 veces la
# densidad de ésta): K_t = 3.035, 3.026, 3.017, 3.011 → converge hacia 3; el
# error en los puntos de Gauss con r ≤ 2a baja a la mitad con cada
# refinamiento (primer orden), y el de r > 5a se estanca en ~1 %: es la
# diferencia entre el desplazamiento uniforme impuesto y la tracción uniforme
# en el infinito que supone Kirsch.
TOLERANCIAS = {
    "error_Kt": 0.015,
    "error_borde_max": 0.05,
    "error_gauss_cercano": 0.08,
    "rms_gauss_cercano": 0.025,
    "error_gauss_lejano": 0.02,
}


def kirsch(x, y, s_inf: float):
    """Campo de Kirsch en cartesianas: (σxx, σyy, σxy) en (x, y)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    r2 = x**2 + y**2
    th = np.arctan2(y, x)
    q, q2 = A**2 / r2, A**4 / r2**2
    c2, s2 = np.cos(2 * th), np.sin(2 * th)
    srr = 0.5 * s_inf * (1 - q) + 0.5 * s_inf * (1 - 4 * q + 3 * q2) * c2
    stt = 0.5 * s_inf * (1 + q) - 0.5 * s_inf * (1 + 3 * q2) * c2
    srt = -0.5 * s_inf * (1 + 2 * q - 3 * q2) * s2
    c, s = np.cos(th), np.sin(th)
    sxx = srr * c**2 + stt * s**2 - 2 * srt * s * c
    syy = srr * s**2 + stt * c**2 + 2 * srt * s * c
    sxy = (srr - stt) * s * c + srt * (c**2 - s**2)
    return sxx, syy, sxy


def calcular() -> dict:
    domain, _, res = resolver_yaml_estatico(AQUI / "modelo.yaml")

    # Esfuerzo nominal lejano: reacción total en el borde cargado / área.
    grupo = domain.physical_groups
    rx = sum(res.reactions_by_node[n]["ux"] for n in grupo["Borde_Carga"])
    s_inf = rx / (H * ESPESOR)

    # Puntos de Gauss: Solidum frente a Kirsch en el mismo punto.
    estados = gauss_states(domain, res.U)
    pts = np.vstack([gs["points_global"] for gs in estados.values()])
    sig = np.vstack([gs["stress"] for gs in estados.values()])
    sk = np.column_stack(kirsch(pts[:, 0], pts[:, 1], s_inf))
    err_pg = np.abs(sig - sk).max(axis=1) / s_inf
    r_pg = np.hypot(pts[:, 0], pts[:, 1])
    cercano = r_pg <= R_CERCANO
    lejano = r_pg > R_LEJANO

    # Esfuerzos nodales recuperados: línea x = 0 y borde del agujero.
    nodal = esfuerzos_nodales(domain, estados)
    coords = {nid: np.asarray(domain.nodes[nid].coordinates) for nid in nodal}

    eje = sorted(grupo["Simetria_X"], key=lambda n: coords[n][1])
    y_eje = np.array([coords[n][1] for n in eje])
    sxx_eje = np.array([nodal[n][0] for n in eje])

    borde = sorted(grupo["Agujero"], key=lambda n: np.arctan2(coords[n][1], coords[n][0]))
    th = np.array([np.arctan2(coords[n][1], coords[n][0]) for n in borde])
    s_b = np.array([nodal[n] for n in borde])
    c, s = np.cos(th), np.sin(th)
    stt_b = s_b[:, 0] * s**2 + s_b[:, 1] * c**2 - 2 * s_b[:, 2] * s * c
    stt_k = s_inf * (1 - 2 * np.cos(2 * th))

    kt = sxx_eje[0] / s_inf        # nodo (0, a): primero de la línea x = 0
    return {
        "a": A, "W": W, "H": H, "espesor": ESPESOR,
        "n_elementos": len(domain.elements), "n_nodos": len(domain.nodes),
        "n_puntos_gauss": int(len(pts)), "n_puntos_cercanos": int(cercano.sum()),
        "reaccion_total": rx, "sigma_inf": s_inf, "sigma_inf_MPa": s_inf / 1e6,
        "Kt": kt, "error_Kt": error_relativo(kt, 3.0),
        "error_Kt_pct": 100 * (kt - 3.0) / 3.0,
        "sigma_min_borde": float(stt_b.min() / s_inf),
        "error_borde_max": float(np.abs(stt_b - stt_k).max() / s_inf),
        "error_borde_max_pct": float(100 * np.abs(stt_b - stt_k).max() / s_inf),
        "error_gauss_cercano": float(err_pg[cercano].max()),
        "error_gauss_cercano_pct": float(100 * err_pg[cercano].max()),
        "rms_gauss_cercano": float(np.sqrt(np.mean(err_pg[cercano] ** 2))),
        "rms_gauss_cercano_pct": float(100 * np.sqrt(np.mean(err_pg[cercano] ** 2))),
        "error_gauss_lejano": float(err_pg[lejano].max()),
        "error_gauss_lejano_pct": float(100 * err_pg[lejano].max()),
        "primer_elemento_radial_sobre_a": float(
            (np.sort(np.hypot(*np.array([domain.nodes[n].coordinates
                                         for n in grupo["Simetria_Y"]]).T))[1] - A) / A),
        "eje_y": y_eje, "eje_sxx": sxx_eje,
        "borde_theta": th, "borde_stt": stt_b,
    }


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    s_inf = res["sigma_inf"]

    # --- Malla y condiciones de contorno
    import meshio
    m = meshio.read(AQUI / "placa_agujero.msh")
    quads = np.vstack([b.data for b in m.cells if b.type == "quad"])
    xy = m.points[:, :2] / A
    fig, (ax, det) = plt.subplots(1, 2, figsize=(7.0, 3.6),
                                  gridspec_kw={"width_ratios": [1.15, 1]})
    ax.add_collection(PolyCollection(xy[quads], facecolor="none", edgecolor="0.55", lw=0.25))
    ax.plot([0, 0], [1, 20], color=COLOR_FE, lw=2.2)
    ax.plot([1, 20], [0, 0], color=COLOR_FE, lw=2.2)
    ax.plot([20, 20], [0, 20], color=COLOR_SECUNDARIO, lw=2.2)
    ax.text(-0.8, 10, r"$u_x = 0$", rotation=90, ha="right", va="center", color=COLOR_FE)
    ax.text(10, -0.8, r"$u_y = 0$", ha="center", va="top", color=COLOR_FE)
    ax.text(20.8, 10, r"$u_x = \bar u$", rotation=90, ha="left", va="center",
            color=COLOR_SECUNDARIO)
    ax.add_patch(plt.Rectangle((0, 0), 2.2, 2.2, fill=False, ls="--", lw=0.8, color="0.2"))
    ax.set_xlim(-3, 23.5)
    ax.set_ylim(-3, 21)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x/a$")
    ax.set_ylabel(r"$y/a$")
    ax.grid(False)
    det.add_collection(PolyCollection(xy[quads], facecolor="none", edgecolor="0.35", lw=0.3))
    det.set_xlim(0, 2.2)
    det.set_ylim(0, 2.2)
    det.set_aspect("equal")
    det.set_xlabel(r"$x/a$  (detalle)")
    det.grid(False)
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_malla")

    # --- σxx en la línea x = 0 y σθθ en el borde del agujero
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 3.3))
    yy = np.linspace(1, 20, 400)
    a1.plot(yy, 1 + 0.5 / yy**2 + 1.5 / yy**4, color=COLOR_ANALITICO, lw=1.3, label="Kirsch")
    y_n = np.asarray(res["eje_y"]) / A
    a1.plot(y_n, np.asarray(res["eje_sxx"]) / s_inf, "o", ms=3.2, color=COLOR_FE,
            label="Solidum (nodos)")
    a1.set_xscale("log")
    a1.set_xlabel(r"$y/a$  (sobre $x = 0$)")
    a1.set_ylabel(r"$\sigma_{xx}/\sigma_\infty$")
    a1.legend()
    tt = np.linspace(0, np.pi / 2, 200)
    a2.plot(np.degrees(tt), 1 - 2 * np.cos(2 * tt), color=COLOR_ANALITICO, lw=1.3)
    a2.plot(np.degrees(res["borde_theta"]), np.asarray(res["borde_stt"]) / s_inf, "o",
            ms=3.2, color=COLOR_FE)
    a2.set_xlabel(r"$\theta$ [°]")
    a2.set_ylabel(r"$\sigma_{\theta\theta}/\sigma_\infty$ en $r = a$")
    a2.set_xticks([0, 30, 60, 90])
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_kirsch")

    # --- Campo σxx
    domain, _, sol = resolver_yaml_estatico(AQUI / "modelo.yaml")
    estados = gauss_states(domain, sol.U)
    polys, vals = [], []
    for eid, gs in estados.items():
        el = domain.elements[eid]
        polys.append(np.array([n.coordinates for n in el.nodes]) / A)
        vals.append(np.mean(np.asarray(gs["stress"])[:, 0]) / s_inf)
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    pc = PolyCollection(polys, array=np.array(vals), cmap="viridis", edgecolor="face", lw=0.4)
    ax.add_collection(pc)
    ax.set_xlim(0, 6)
    ax.set_ylim(0, 6)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x/a$")
    ax.set_ylabel(r"$y/a$")
    ax.grid(False)
    fig.colorbar(pc, ax=ax, label=r"$\sigma_{xx}/\sigma_\infty$ (promedio por elemento)")
    guardar_figura(fig, carpeta, "fig_campo")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    print(f"σ∞ = {res['sigma_inf_MPa']:.3f} MPa   K_t = {res['Kt']:.4f}  (error {res['error_Kt']:.2%})")
    print(f"Borde del agujero: error máx {res['error_borde_max']:.2%} de σ∞")
    print(f"Puntos de Gauss r ≤ 2a: máx {res['error_gauss_cercano']:.2%}, rms "
          f"{res['rms_gauss_cercano']:.2%};  r > 5a: máx {res['error_gauss_lejano']:.2%}")
    return res


if __name__ == "__main__":
    main()
