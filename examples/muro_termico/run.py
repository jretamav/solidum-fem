"""Ejemplo 9 del manual de ejemplos — conducción térmica transitoria en un muro.

Muro de concreto de espesor ``L``, inicialmente a ``T₀``. En ``t = 0`` su cara
``x = 0`` pasa bruscamente a ``T_s``; la cara ``x = L`` es adiabática. La
solución exacta es la serie de Carslaw y Jaeger (1959, §3.3):

    (T − T_s)/(T₀ − T_s) = Σ 4/((2n+1)π) · sin(λ_n x) · exp(−α λ_n² t),
    λ_n = (2n + 1)·π/(2L),   α = k/(ρ c)

Se integra con el método θ (``ThetaMethodSolver``): Euler implícito (θ = 1,
primer orden, L-estable) y Crank-Nicolson (θ = 1/2, segundo orden, no
L-estable), con capacidad concentrada o consistente.

Uso::

    python examples/muro_termico/run.py
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
from solidum.elements.thermal import Quad4Thermal  # noqa: E402
from solidum.materials.thermal_conduction import ThermalConduction  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import ThetaMethodSolver  # noqa: E402

K, RHO, CP = 1.4, 2400.0, 900.0      # concreto: W/(m·K), kg/m³, J/(kg·K)
ALPHA = K / (RHO * CP)               # difusividad [m²/s]
L, N = 0.20, 40                      # espesor [m] y elementos
H = L / N
T0, TS = 20.0, 80.0                  # temperatura inicial y de la cara [°C]
HORA = 3600.0
T_FIN = 4 * HORA                     # instante de comparación

PASOS_CONVERGENCIA = (1200.0, 600.0, 300.0, 150.0, 75.0)   # Δt [s]

TOLERANCIAS = {
    # Euler implícito: nunca sale del rango de los datos, orden 1.
    "euler_600.violacion": 1e-10,
    "orden_euler": (0.9, 1.1),
    # Crank-Nicolson: orden 2 con pasos pequeños; con Δt = 600 s oscila y
    # sale del rango (el nodo junto a la cara supera T_s).
    "orden_cn": (1.8, 2.2),
    "cn_600.T_max": (95.0, 115.0),
    "cn_300.T_max": (85.0, 100.0),
    "cn_150.violacion": 1e-10,
    # Capacidad consistente: por debajo de Δt ≈ h²/(6α) la temperatura baja
    # de T₀ antes de subir; la concentrada nunca.
    "consistente_1.T_min": (18.5, 19.5),
    "concentrada_1.violacion": 1e-10,
    "consistente_8.violacion": 1e-10,
}


def muro():
    """Franja de ``N`` elementos ``Quad4Thermal`` a lo ancho del muro; los
    bordes sin condición son adiabáticos. Temperatura ``T_s`` en ``x = 0``."""
    dom = Domain()
    concreto = ThermalConduction(k=K, c=CP, density=RHO)
    nodos = {}
    for j in range(2):
        for i in range(N + 1):
            nodos[(i, j)] = dom.add_node(len(nodos) + 1, [i * H, j * H])
    for i in range(N):
        dom.add_element(Quad4Thermal(i + 1, [nodos[(i, 0)], nodos[(i + 1, 0)],
                                             nodos[(i + 1, 1)], nodos[(i, 1)]], concreto))
    for j in range(2):
        nodos[(0, j)].fix_dof("T", TS)
    dom.generate_equation_numbers(verbose=False)
    return dom, nodos


def exacta(x, t, terminos: int = 400):
    """Serie de Carslaw y Jaeger para la cara ``x = 0`` a ``T_s`` y la cara
    ``x = L`` adiabática."""
    n = np.arange(terminos)
    lam = (2 * n + 1) * np.pi / (2 * L)
    v = np.sum(4 / ((2 * n + 1) * np.pi) * np.sin(np.outer(np.atleast_1d(x), lam))
               * np.exp(-ALPHA * lam**2 * t), axis=1)
    return TS + (T0 - TS) * v


def integrar(theta: float, dt: float, t_fin: float, capacidad: str = "lumped"):
    """Temperatura en la línea ``y = 0`` en todos los instantes."""
    dom, nodos = muro()
    solver = ThetaMethodSolver(Assembler(dom), dt=dt, n_steps=int(round(t_fin / dt)),
                               T_initial=T0, theta=theta, lumping=capacidad)
    resultado = solver.solve()
    filas = [nodos[(i, 0)].dofs["T"] for i in range(N + 1)]
    x = np.array([nodos[(i, 0)].coordinates[0] for i in range(N + 1)])
    return np.asarray(resultado.t_history), np.asarray(resultado.T_history)[filas], x


def _medir(t, T, x) -> dict:
    """Error final frente a la serie y salida del rango [T₀, T_s]."""
    fuera = max(0.0, T0 - T[:, 1:].min(), T[:, 1:].max() - TS)
    return {"error_final": float(np.max(np.abs(T[:, -1] - exacta(x, t[-1]))) / (TS - T0)),
            "violacion": float(fuera / (TS - T0)),
            "T_min": float(T[:, 1:].min()), "T_max": float(T[:, 1:].max())}


def calcular() -> dict:
    out = {"L_cm": 100 * L, "N": N, "h_mm": 1e3 * H, "k": K, "rho": RHO, "cp": CP,
           "alpha": ALPHA, "T0": T0, "Ts": TS, "t_fin_h": T_FIN / HORA,
           "escala_h": L * L / ALPHA / HORA, "h2_alpha_s": H * H / ALPHA,
           "umbral_consistente_s": H * H / (6 * ALPHA)}

    # Convergencia en t = 4 h. Frente a la serie exacta el error incluye el
    # espacial de la malla; el temporal se aísla frente a una referencia en
    # la misma malla con Δt = 1 s (Crank-Nicolson).
    _, T_ref, _ = integrar(0.5, 1.0, T_FIN)
    for theta, nombre in ((1.0, "euler"), (0.5, "cn")):
        errores, temporales = [], []
        for dt in PASOS_CONVERGENCIA:
            t, T, x = integrar(theta, dt, T_FIN)
            errores.append(_medir(t, T, x)["error_final"])
            temporales.append(float(np.max(np.abs(T[:, -1] - T_ref[:, -1])) / (TS - T0)))
        out[f"errores_{nombre}"] = errores
        out[f"temporales_{nombre}"] = temporales
        # Pendiente log-log del error temporal entre los dos pasos más pequeños.
        out[f"orden_{nombre}"] = float(np.log(temporales[-2] / temporales[-1]) / np.log(2))
        out[f"{nombre}_por_paso"] = {str(int(dt)): {"error": e, "temporal": et}
                                     for dt, e, et in zip(PASOS_CONVERGENCIA, errores, temporales)}
    out["error_espacial"] = out["errores_cn"][-1]

    # Paso de 10 min con los dos esquemas.
    for theta, clave in ((1.0, "euler_600"), (0.5, "cn_600")):
        t, T, x = integrar(theta, 600.0, T_FIN)
        out[clave] = _medir(t, T, x)
        out[clave]["historia_1"] = {"t_h": t / HORA, "T": T[1]}      # nodo a h de la cara
    t = np.linspace(0, T_FIN, 400)
    out["exacta_1"] = {"t_h": t / HORA, "T": np.array([exacta(H, ti)[0] for ti in t])}

    # Umbral de Crank-Nicolson: con 5 min todavía sale del rango; con 2.5 min no.
    for dt in (300.0, 150.0):
        t, T, x = integrar(0.5, dt, 3 * dt)
        out[f"cn_{int(dt)}"] = {**_medir(t, T, x), "pasos_por_h2_alpha": dt / (H * H / ALPHA)}

    # Capacidad consistente frente a concentrada con pasos muy pequeños.
    for dt in (1.0, 8.0):
        for capacidad, nombre in (("consistent", "consistente"), ("lumped", "concentrada")):
            t, T, x = integrar(1.0, dt, 30 * dt, capacidad)
            out[f"{nombre}_{int(dt)}"] = _medir(t, T, x)

    # Perfiles con Euler implícito y Δt = 60 s.
    perfiles = {}
    for horas in (0.25, 1.0, 4.0, 12.0):
        t, T, x = integrar(1.0, 60.0, horas * HORA)
        perfiles[f"{horas:g}"] = {"x_cm": 100 * x, "T": T[:, -1],
                                  "T_exacta": exacta(x, t[-1])}
    out["perfiles"] = perfiles
    return out


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt

    fig, (ax, bx, cx) = plt.subplots(1, 3, figsize=(7.4, 2.9))
    colores = [COLOR_FE, COLOR_SECUNDARIO, COLOR_TERCIARIO, "#7570b3"]
    for c, (horas, p) in zip(colores, res["perfiles"].items()):
        ax.plot(p["x_cm"], p["T_exacta"], color=COLOR_ANALITICO, lw=3.0, alpha=0.25)
        ax.plot(p["x_cm"], p["T"], "o", color=c, ms=2.5, markevery=2, label=f"{horas} h")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("T [°C]")
    ax.legend(fontsize=7, title="Euler, Δt = 1 min", title_fontsize=7,
              loc="upper right", bbox_to_anchor=(1.02, 0.86), ncol=2, columnspacing=0.6)

    e = res["exacta_1"]
    bx.plot(e["t_h"], e["T"], color=COLOR_ANALITICO, lw=3.0, alpha=0.25, label="exacta")
    for clave, c, etiqueta in (("euler_600", COLOR_FE, "Euler"), ("cn_600", COLOR_SECUNDARIO, "Crank-Nicolson")):
        h = res[clave]["historia_1"]
        bx.plot(h["t_h"], h["T"], "o-", color=c, ms=2.5, lw=0.8, label=etiqueta)
    bx.axhline(res["Ts"], color="0.5", lw=0.6, ls=":")
    bx.set_xlim(0, 2)
    bx.set_xlabel("t [h]")
    bx.set_ylabel(f"T a {res['h_mm']:.0f} mm de la cara [°C]")
    bx.legend(fontsize=7, title="Δt = 10 min", title_fontsize=7)

    dts = np.array([1200, 600, 300, 150, 75])
    cx.loglog(dts, res["temporales_euler"], "o-", color=COLOR_FE, ms=4, label="Euler")
    cx.loglog(dts, res["temporales_cn"], "s-", color=COLOR_SECUNDARIO, ms=4, label="Crank-Nicolson")
    cx.axhline(res["error_espacial"], color="0.5", lw=0.8, ls="--")
    cx.text(80, res["error_espacial"] * 1.3, "error de la malla", fontsize=7, color="0.4")
    cx.set_xlabel("Δt [s]")
    cx.set_ylabel("error temporal a las 4 h")
    cx.legend(fontsize=7)
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_muro")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    print(f"α = {res['alpha']:.3e} m²/s, L²/α = {res['escala_h']:.1f} h, h²/α = {res['h2_alpha_s']:.1f} s")
    print("temporal Euler:", res["temporales_euler"], f" orden {res['orden_euler']:.2f}")
    print("temporal CN   :", res["temporales_cn"], f" orden {res['orden_cn']:.2f}")
    print(f"error espacial (CN, Δt = 75 s, frente a la serie): {res['error_espacial']:.1e}")
    for clave in ("euler_600", "cn_600", "consistente_1", "concentrada_1", "consistente_8", "concentrada_8"):
        d = res[clave]
        print(f"{clave:14s} T ∈ [{d['T_min']:.3f}, {d['T_max']:.3f}]  fuera de rango {d['violacion']:.3f}")
    return res


if __name__ == "__main__":
    main()
