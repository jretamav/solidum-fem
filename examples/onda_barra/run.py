"""Ejemplo 8 del manual de ejemplos — propagación de una onda en una barra.

Barra de acero empotrada en ``x = 0`` y libre en ``x = L``, con una fuerza
en escalón ``F₀`` aplicada en el extremo libre en ``t = 0``. La solución
exacta (d'Alembert) es una onda de esfuerzo que viaja a ``c = √(E/ρ)`` y se
refleja en los extremos:

- el extremo libre oscila en diente de sierra entre 0 y ``2·F₀L/(EA)``, con
  periodo ``4L/c``;
- el esfuerzo en el empotramiento es una onda cuadrada: 0 hasta ``t = L/c``,
  ``2σ₀`` hasta ``3L/c``, 0 hasta ``5L/c``… (``σ₀ = F₀/A``).

Se comparan las diferencias centradas explícitas (masa concentrada) y los
métodos implícitos de Newmark (trapezoidal) y HHT-α (masa consistente).

Uso::

    python examples/onda_barra/run.py
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
from solidum.elements.truss import Truss2D  # noqa: E402
from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import CentralDifferenceSolver, HHTSolver, NewmarkSolver  # noqa: E402

L, N = 1.0, 100            # longitud [m] y elementos
A = 1.0e-4                 # sección [m²] (1 cm²)
E, RHO = 200.0e9, 7850.0   # acero
F0 = 1.0e3                 # fuerza en escalón [N]

C = np.sqrt(E / RHO)       # velocidad de la onda
H = L / N                  # tamaño de elemento
SIGMA0 = F0 / A
U_EST = F0 * L / (E * A)   # desplazamiento estático del extremo
PASO_CFL = H / C           # tiempo que tarda la onda en cruzar un elemento

# (clave, solver, Δt / (h/c), opciones)
CASOS = [("explicito_1", CentralDifferenceSolver, 1.0, {}),
         ("explicito_05", CentralDifferenceSolver, 0.5, {}),
         ("newmark_1", NewmarkSolver, 1.0, {}),
         ("hht_1", HHTSolver, 1.0, {"alpha": -0.3}),
         ("newmark_4", NewmarkSolver, 4.0, {}),
         ("hht_4", HHTSolver, 4.0, {"alpha": -0.3})]

TOLERANCIAS = {
    # Explícito con masa concentrada y Δt = h/c: exacto en los nodos.
    "explicito_1.error_u": 1e-10,
    "explicito_1.sobrepaso": 1e-10,
    # Con la mitad del paso aparece la dispersión.
    "explicito_05.sobrepaso": (0.2, 0.3),
    "explicito_05.error_u": (0.01, 0.03),
    # Implícitos: el mismo sobrepaso; HHT lo reduce poco.
    "newmark_1.sobrepaso": (0.2, 0.3),
    "hht_reduce_sobrepaso": (0.005, 0.05),
    # ... y amortigua la cola de oscilaciones más deprisa.
    "hht_sobre_newmark_cola": (0.0, 0.9),
    # Con Δt cuatro veces mayor siguen estables, con más error.
    "newmark_4.error_u": (0.04, 0.08),
    # Límite de estabilidad del explícito: 2/ω_max = h/c.
    "dt_critico_sobre_h_c": (1.0, 1.0001),
    "inestable.diverge": (1, 1),
}


def barra():
    """Barra de ``N`` elementos ``Truss2D`` empotrada en ``x = 0``. Devuelve el
    dominio, los nodos y el vector de la fuerza en el extremo libre."""
    dom = Domain()
    acero = Elastic1D(E=E, density=RHO)
    nodos = [dom.add_node(i + 1, [i * H, 0.0]) for i in range(N + 1)]
    for k in range(N):
        dom.add_element(Truss2D(k + 1, [nodos[k], nodos[k + 1]], acero, A=A))
    nodos[0].fix_dof("ux", 0.0)
    for nodo in nodos:
        nodo.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[nodos[-1].dofs["ux"]] = F0
    return dom, nodos, F


def exacta(t):
    """Desplazamiento del extremo libre y esfuerzo en el empotramiento."""
    tau = np.mod(t, 4 * L / C)
    u = np.where(tau <= 2 * L / C, U_EST * C * tau / L, U_EST * (4 - C * tau / L))
    sigma = np.where((tau > L / C) & (tau < 3 * L / C), 2 * SIGMA0, 0.0)
    return u, sigma


def integrar(solver_cls, fraccion_cfl: float, **opciones):
    """Integra una vuelta completa de la onda (``4L/c``) con Δt = fracción·h/c."""
    dom, nodos, F = barra()
    solver = solver_cls(Assembler(dom), 4 * L / C, fraccion_cfl * PASO_CFL,
                        F_func=lambda t: F, **opciones)
    resultado = solver.solve()
    t = np.asarray(resultado.t_history)
    U = np.asarray(resultado.u_history)
    u_libre = U[nodos[-1].dofs["ux"]]
    # Esfuerzo del primer elemento: E·(u₁ − u₀)/h, con u₀ = 0.
    sigma_empotramiento = E * U[nodos[1].dofs["ux"]] / H
    return t, u_libre, sigma_empotramiento


def calcular() -> dict:
    out = {"L_m": L, "N": N, "h_cm": 100 * H, "A_cm2": A * 1e4, "E_GPa": E / 1e9,
           "rho": RHO, "c_m_s": C, "F0_kN": F0 / 1e3, "sigma0_MPa": SIGMA0 / 1e6,
           "u_est_um": 1e6 * U_EST, "t_cruce_us": 1e6 * L / C,
           "paso_cfl_us": 1e6 * PASO_CFL}
    historias = {}
    for clave, solver_cls, fraccion, opciones in CASOS:
        t, u, sigma = integrar(solver_cls, fraccion, **opciones)
        u_ex, sigma_ex = exacta(t)
        # Oscilación de la cola, lejos de los frentes: 1.5 < t·c/L < 2.8.
        cola = (t * C / L > 1.5) & (t * C / L < 2.8)
        out[clave] = {"fraccion_cfl": fraccion, "pasos": len(t) - 1,
                      "error_u": float(np.max(np.abs(u - u_ex)) / U_EST),
                      "sobrepaso": float(sigma.max() / (2 * SIGMA0) - 1),
                      "sobrepaso_pct": float(100 * (sigma.max() / (2 * SIGMA0) - 1)),
                      "error_u_pct": float(100 * np.max(np.abs(u - u_ex)) / U_EST),
                      "cola_rms": float(np.sqrt(np.mean(((sigma - sigma_ex)[cola]
                                                         / (2 * SIGMA0))**2)))}
        historias[clave] = (t * C / L, sigma / (2 * SIGMA0))
    out["hht_reduce_sobrepaso"] = out["newmark_1"]["sobrepaso"] - out["hht_1"]["sobrepaso"]
    out["hht_sobre_newmark_cola"] = out["hht_1"]["cola_rms"] / out["newmark_1"]["cola_rms"]

    # Límite de estabilidad del explícito: 2/ω_max del problema con masa
    # concentrada, frente al tiempo de cruce de un elemento.
    dom, nodos, _ = barra()
    asm = Assembler(dom)
    K, _ = asm.assemble_non_linear_system(np.zeros(dom.total_dofs))
    M = asm.assemble_mass_matrix("lumped")
    libres = [n.dofs["ux"] for n in nodos[1:]]
    k_diag = K.tocsr()[libres][:, libres].toarray()
    m_diag = M.tocsr()[libres][:, libres].diagonal()
    w_max = np.sqrt(np.max(np.linalg.eigvalsh(k_diag / np.sqrt(np.outer(m_diag, m_diag)))))
    out["dt_critico_sobre_h_c"] = float(2 / w_max / PASO_CFL)
    try:
        integrar(CentralDifferenceSolver, 1.005)
        out["inestable"] = {"diverge": 0, "fraccion_cfl": 1.005}
    except RuntimeError:
        out["inestable"] = {"diverge": 1, "fraccion_cfl": 1.005}

    out["historias"] = {k: {"t": v[0], "sigma": v[1]} for k, v in historias.items()
                        if k in ("explicito_1", "explicito_05", "newmark_1", "hht_1")}
    t = np.linspace(0, 4 * L / C, 2001)
    out["historias"]["exacta"] = {"t": t * C / L, "sigma": exacta(t)[1] / (2 * SIGMA0)}
    return out


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt

    h = res["historias"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True)
    for eje in (ax, bx):
        eje.plot(h["exacta"]["t"], h["exacta"]["sigma"], color=COLOR_ANALITICO, lw=3.5,
                 alpha=0.25, label="exacta")
        eje.set_xlabel(r"$t\,c/L$")
        eje.set_xlim(0, 4)
    ax.plot(h["explicito_05"]["t"], h["explicito_05"]["sigma"], color=COLOR_FE, lw=0.9,
            label=r"$\Delta t = 0.5\,h/c$")
    ax.plot(h["explicito_1"]["t"], h["explicito_1"]["sigma"], color=COLOR_SECUNDARIO, lw=1.2,
            label=r"$\Delta t = h/c$")
    ax.set_title("diferencias centradas, masa concentrada", fontsize=9)
    bx.plot(h["newmark_1"]["t"], h["newmark_1"]["sigma"], color=COLOR_FE, lw=0.9,
            label="Newmark trapezoidal")
    bx.plot(h["hht_1"]["t"], h["hht_1"]["sigma"], color=COLOR_TERCIARIO, lw=0.9,
            label=r"HHT, $\alpha = -0.3$")
    bx.set_title(r"implícitos, masa consistente, $\Delta t = h/c$", fontsize=9)
    ax.set_ylabel(r"$\sigma(0, t) / 2\sigma_0$")
    ax.set_ylim(-0.4, 1.4)
    for eje in (ax, bx):
        eje.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_onda")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    print(f"c = {res['c_m_s']:.1f} m/s, h/c = {res['paso_cfl_us']:.3f} µs, "
          f"Δt_crít = {res['dt_critico_sobre_h_c']:.6f}·h/c; con 1.005·h/c diverge: {res['inestable']['diverge']}")
    for clave, *_ in CASOS:
        d = res[clave]
        print(f"{clave:13s} Δt = {d['fraccion_cfl']}·h/c  pasos {d['pasos']:4d}  "
              f"error u_L {d['error_u']:.2e}  sobrepaso σ {d['sobrepaso']:+.3f}")
    return res


if __name__ == "__main__":
    main()
