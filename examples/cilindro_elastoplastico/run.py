"""Ejemplo 4 del manual de ejemplos — cilindro elastoplástico (J2).

El cilindro del ejemplo 3, ahora de material elastoplástico perfecto de von
Mises, en deformación plana, bajo presión interna creciente. Dos referencias
exactas:

- presión de primera fluencia (en la cara interior, con σzz = ν(σrr + σθθ)):
  p_e = σy (b² − a²) / √(3 b⁴ + (1 − 2ν)² a⁴);
- presión de colapso (análisis límite; no depende de E ni de ν):
  p_lim = (2/√3) σy ln(b/a).

La curva presión-desplazamiento se traza con el método de longitud de arco,
que sigue la meseta de colapso donde el control de carga no puede pasar. Se
comparan cinco elementos: los lineales con integración completa muestran
bloqueo volumétrico en el régimen plástico, que es isocórico.

Uso::

    python examples/cilindro_elastoplastico/run.py
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))

from ejemplos_comun import (  # noqa: E402
    COLOR_ANALITICO, COLOR_FE, COLOR_SECUNDARIO, COLOR_TERCIARIO,
    estilo_figuras, guardar_figura, guardar_resultados,
)

import solidum  # noqa: E402
from solidum.materials.von_mises_2d import VonMises2D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import ArcLengthSolver  # noqa: E402

# Geometría, malla y carga de presión: las del ejemplo 3.
_spec = importlib.util.spec_from_file_location("cilindro_lame", AQUI.parent / "cilindro_lame" / "run.py")
lame_ej = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lame_ej)
modelo, carga_de_presion = lame_ej.modelo, lame_ej.carga_de_presion
RI, RE, E, NU = lame_ej.RI, lame_ej.RE, lame_ej.E, lame_ej.NU

SY = 240.0e6                                    # límite elástico [Pa]
P_E = SY * (RE**2 - RI**2) / np.sqrt(3 * RE**4 + (1 - 2 * NU) ** 2 * RI**4)
P_LIM = 2 / np.sqrt(3) * SY * np.log(RE / RI)
U_B_ELASTICO = lame_ej.lame(RE)[2] / lame_ej.P   # u_r(b) por unidad de presión

# (etiqueta, elemento, celdas n × n, cuadratura). Mallas con un número de
# grados de libertad parecido: 8 × 8 para los lineales, 4 × 4 para los cuadráticos.
CASOS = [("Tri3", "Tri3", 8, None), ("Quad4", "Quad4", 8, None),
         ("Tri6", "Tri6", 4, None), ("Quad8", "Quad8", 4, None),
         ("Quad8R", "Quad8", 4, "2x2")]
DL = 2.0e-5          # longitud de arco inicial [m]: del orden del desplazamiento elástico
PASOS = 150
U_LECTURA = (0.5e-3, 1.0e-3)   # desplazamientos exteriores donde se lee la carga

TOLERANCIAS = {
    **{f"{c}.error_flexibilidad": 0.001 for c in ("Tri6", "Quad8", "Quad8R")},
    "Quad4.error_flexibilidad": 0.01,
    "Tri3.error_flexibilidad": 0.03,
    **{f"{c}.error_colapso": 0.001 for c in ("Tri6", "Quad8", "Quad8R")},
    "Tri3.error_colapso": 0.01,
    # Bloqueo: el Quad4 sigue ganando carga entre 0.5 y 1 mm (sin meseta), con
    # un esfuerzo medio muy fuera del intervalo exacto; los que no se bloquean
    # quedan dentro salvo el error de discretización.
    "Quad4.crecimiento_meseta": (0.005, 0.05),
    **{f"{c}.crecimiento_meseta": (-0.002, 0.002) for c in ("Tri6", "Quad8", "Quad8R")},
    "Quad4.medio_fuera_MPa": (50.0, 1000.0),
    **{f"{c}.medio_fuera_MPa": 10.0 for c in ("Tri6", "Quad8", "Quad8R")},
    # El paso de arco excesivo aterriza en λ = λ_max con desplazamientos absurdos.
    "paso_excesivo.lambda": (1.49, 1.51),
}


def trazar(elemento: str, n: int, cuadratura, dl: float = DL, pasos: int = PASOS):
    """Curva (λ, u_b) con p = λ·p_lim, u_b = desplazamiento radial exterior.

    Devuelve también el esfuerzo medio de los puntos de Gauss en el primer
    paso con u_b ≥ 1 mm (``None`` si no se llega): un punto fijo del recorrido,
    comparable entre elementos."""
    material = VonMises2D(E=E, nu=NU, sigma_y=SY, H=0.0, hypothesis="plane_strain")
    domain, aristas, nodo = modelo(elemento, n, material=material, cuadratura=cuadratura)
    F = carga_de_presion(domain, aristas, P_LIM)          # λ = 1 ⇔ presión de colapso
    assembler = Assembler(domain)
    solver = ArcLengthSolver(assembler, max_iter=30, max_lambda=1.5,
                             initial_dl=dl, max_steps=pasos)
    exterior = nodo[max(i for i, _ in nodo), 0]           # (b, 0)
    curva, medio = [], []

    def al_converger(paso, U, lam):
        # El solver llama aquí con el estado ya consolidado (commit).
        curva.append((lam, U[exterior.dofs["ux"]]))
        if not medio and curva[-1][1] >= U_LECTURA[1]:
            medio.append(presion_media(domain))

    solidum.run(domain, assembler=assembler, solver=solver, F_applied=F,
                step_callback=al_converger)
    lam, ub = np.array(curva).T
    return domain, lam, ub, (medio[0] if medio else None)


def presion_media(domain) -> np.ndarray:
    """Esfuerzo medio (σxx + σyy + σzz)/3 en cada punto de Gauss, en el
    último estado convergido."""
    out = []
    for el in domain.elements.values():
        for k, sig in enumerate(el.state.stresses):
            szz = el.material.out_of_plane_stress(sig, el.state.vars[k])
            out.append((sig[0] + sig[1] + szz) / 3)
    return np.array(out)


def calcular() -> dict:
    # En el colapso σzz = (σrr + σθθ)/2 y σθθ − σrr = 2k, así que el esfuerzo
    # medio es σrr + k, entre −p_lim + k (cara interior) y k (exterior).
    k = SY / np.sqrt(3)
    medio_lo, medio_hi = -P_LIM + k, k
    out = {"sigma_y_MPa": SY / 1e6, "p_e_MPa": P_E / 1e6, "p_lim_MPa": P_LIM / 1e6,
           "p_e_sobre_p_lim": P_E / P_LIM, "dl": DL, "pasos": PASOS,
           "u_lectura_mm": [1e3 * u for u in U_LECTURA],
           "medio_exacto_min_MPa": medio_lo / 1e6, "medio_exacto_max_MPa": medio_hi / 1e6}
    for etiqueta, elemento, n, cuad in CASOS:
        domain, lam, ub, pm = trazar(elemento, n, cuad)
        flex = ub[0] / (lam[0] * P_LIM)                    # primer paso: elástico
        l05, l10 = (float(np.interp(u, ub, lam)) for u in U_LECTURA)
        out[etiqueta] = {
            "elemento": elemento, "n": n, "cuadratura": cuad or "completa",
            "gdl": int(domain.total_dofs),
            "lambda": lam, "u_b_mm": 1e3 * ub,
            "error_flexibilidad": abs(flex / U_B_ELASTICO - 1),
            "lambda_05": l05, "lambda_10": l10,
            "error_colapso": abs(l10 - 1), "crecimiento_meseta": l10 - l05,
            "p_10_MPa": l10 * P_LIM / 1e6,
            "medio_min_MPa": pm.min() / 1e6, "medio_max_MPa": pm.max() / 1e6,
            # Cuánto se sale el esfuerzo medio del intervalo exacto [MPa].
            "medio_fuera_MPa": max(medio_lo - pm.min(), pm.max() - medio_hi, 0.0) / 1e6,
        }
    # Advertencia del capítulo: paso de arco inicial por omisión (0.1 m).
    _, lam, ub, _ = trazar("Quad8", 4, None, dl=0.1, pasos=1)
    out["paso_excesivo"] = {"dl": 0.1, "lambda": float(lam[-1]), "u_b_m": float(ub[-1]),
                            "u_b_sobre_elastico": float(ub[-1] / (U_B_ELASTICO * P_LIM))}
    return out


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt

    estilos = {"Tri3": ("v", COLOR_TERCIARIO), "Quad4": ("s", COLOR_FE),
               "Tri6": ("^", COLOR_SECUNDARIO), "Quad8": ("o", COLOR_ANALITICO),
               "Quad8R": ("D", "#7570b3")}
    fig, (ax, det) = plt.subplots(1, 2, figsize=(7.2, 3.6))
    ub_el = np.linspace(0, U_B_ELASTICO * P_E * 1e3, 10)
    ax.plot(ub_el, ub_el / (U_B_ELASTICO * P_LIM * 1e3), color="0.6", lw=4, alpha=0.5,
            label="elástico (Lamé)")
    for a in (ax, det):
        a.axhline(1.0, color=COLOR_ANALITICO, lw=1.0, ls="--")
    ax.axhline(res["p_e_sobre_p_lim"], color="0.5", lw=0.8, ls=":")
    ax.text(0.59, 1.02, r"$p_{\mathrm{lim}}$", ha="right", va="bottom")
    ax.text(0.59, res["p_e_sobre_p_lim"] + 0.02, r"$p_e$", ha="right", color="0.4")
    etiquetas = {"Quad8R": "Quad8 (2×2)"}
    for c, (mk, col) in estilos.items():
        d = res[c]
        ax.plot(d["u_b_mm"], d["lambda"], "-", color=col, lw=1.0, marker=mk, ms=3,
                markevery=3)
        det.plot(d["u_b_mm"], d["lambda"], "-", color=col, lw=1.2, marker=mk, ms=3,
                 markevery=12, label=etiquetas.get(c, c))
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel(r"$u_r(b)$ [mm]")
    ax.set_ylabel(r"$p / p_{\mathrm{lim}}$")
    ax.legend(loc="lower right", fontsize=8)
    det.set_xlim(0.3, 1.65)
    det.set_ylim(0.985, 1.03)
    det.set_xlabel(r"$u_r(b)$ [mm]  (detalle de la meseta)")
    det.legend(loc="upper left", fontsize=8, ncol=2)
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_curvas")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    print(f"p_e = {res['p_e_MPa']:.2f} MPa   p_lim = {res['p_lim_MPa']:.2f} MPa")
    for c, *_ in CASOS:
        d = res[c]
        print(f"{c:7s} gdl={d['gdl']:4d}  flexibilidad {d['error_flexibilidad']:.2e}  "
              f"λ(0.5 mm)={d['lambda_05']:.4f}  λ(1 mm)={d['lambda_10']:.4f}")
    print(f"dl = 0.1: λ = {res['paso_excesivo']['lambda']:.3f}, u_b = {res['paso_excesivo']['u_b_m']:.2f} m")
    return res


if __name__ == "__main__":
    main()
