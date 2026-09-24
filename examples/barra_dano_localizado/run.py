"""Ejemplo 5 del manual de ejemplos — barra con daño: localización y retroceso.

Barra a tracción con daño isótropo y ablandamiento exponencial
(``IsotropicDamage1D``), empotrada en un extremo y cargada en el otro. El
elemento central tiene una resistencia un 5 % menor: el daño se localiza en
él y el resto de la barra se descarga elásticamente. Con el esfuerzo uniforme,
la respuesta exacta depende sólo de la deformación ``ε_w`` del elemento débil:

    σ(ε_w) = (1 − d)·E·ε_w,   d = 1 − (κ₀w/ε_w)·exp(−α(ε_w − κ₀w))
    F = A·σ,                   u = h·ε_w + (L − h)·σ/E

Dos barras con el mismo tamaño de elemento ``h``: en la corta la curva baja
sin más; en la larga retrocede (snap-back), porque libera más energía elástica
de la que el elemento dañado puede disipar. Se comparan tres controles: de
carga (se detiene en el pico), de desplazamiento del extremo (no pasa el
retroceso) e indirecto sobre el alargamiento del elemento débil (sigue las
dos curvas exactas).

Uso::

    python examples/barra_dano_localizado/run.py
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

from solidum.constants import DAMAGE_MAX  # noqa: E402
from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.truss import Truss2D  # noqa: E402
from solidum.materials.damage_1d import IsotropicDamage1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import (  # noqa: E402
    ArcLengthSolver, IndirectDisplacementSolver, NonlinearSolver,
)

E = 30.0e9            # módulo de Young [Pa]
A = 0.01              # sección [m²]
KAPPA_0 = 1.0e-4      # deformación umbral del daño
ALPHA = 5000.0        # velocidad de ablandamiento
DEBIL = 0.95          # el elemento central resiste un 5 % menos
H = 0.1               # tamaño de elemento [m], igual en las dos barras
BARRAS = {"corta": 0.2, "larga": 1.0}   # longitudes [m]

K0W = DEBIL * KAPPA_0
F_PICO = A * E * K0W                      # carga máxima: fluye el elemento débil
DLAMBDA = 0.02        # paso del control indirecto: 2 % de la carga de pico
DLAMBDA_GRANDE = 0.1  # paso de la advertencia (mayor que el margen del 5 %)

TOLERANCIAS = {
    # El control de carga converge hasta el pico y no más allá.
    **{f"{b}.carga.F_max_sobre_pico": (0.99, 1.0 + 1e-9) for b in BARRAS},
    # Control de desplazamiento: exacto en la corta; en la larga se detiene
    # en el pico, al empezar el retroceso.
    "corta.desplazamiento.error_F": 1e-12,
    "larga.desplazamiento.u_final_sobre_pico": (0.99, 1.01),
    # Control indirecto: exacto en las dos, con un solo elemento dañado.
    **{f"{b}.indirecto.error_F": 1e-12 for b in BARRAS},
    **{f"{b}.indirecto.error_u": 1e-12 for b in BARRAS},
    **{f"{b}.indirecto.danados": (1, 1) for b in BARRAS},
    "larga.indirecto.error_u_min": 1e-3,
    # El arco cilíndrico no pasa del pico de la barra larga: no retrocede.
    "larga.cilindrico.u_min_sobre_pico": (0.99, 100.0),
    # Pasos grandes: otro equilibrio, con todos los elementos dañados.
    "salto.danados": (10, 10),
    "salto.error_u": (0.1, 100.0),
    "salto_desplazamiento.danados": (2, 2),
}


def exacta(eps_w, L: float):
    """Solución exacta ``(u, F)`` con el daño localizado en el elemento débil,
    en función de su deformación ``ε_w`` (incluye la rama elástica y el tope
    ``DAMAGE_MAX`` del material)."""
    eps_w = np.asarray(eps_w, dtype=float)
    d = np.where(eps_w > K0W, 1.0 - (K0W / eps_w) * np.exp(-ALPHA * (eps_w - K0W)), 0.0)
    sigma = (1.0 - np.minimum(d, DAMAGE_MAX)) * E * eps_w
    return H * eps_w + (L - H) * sigma / E, A * sigma


def barra(L: float):
    """Barra de ``L/H`` elementos ``Truss2D`` con el central debilitado,
    empotrada en ``x = 0``. Devuelve el dominio, el nodo del extremo libre y
    el elemento débil; cada control añade su carga o su apoyo."""
    n = round(L / H)
    dom = Domain()
    nodos = [dom.add_node(i + 1, [i * H, 0.0]) for i in range(n + 1)]
    for k in range(n):
        umbral = K0W if k == n // 2 else KAPPA_0
        material = IsotropicDamage1D(E=E, kappa_0=umbral, alpha=ALPHA)
        dom.add_element(Truss2D(k + 1, [nodos[k], nodos[k + 1]], material, A=A))
    nodos[0].fix_dof("ux", 0.0)
    for nodo in nodos:
        nodo.fix_dof("uy", 0.0)
    extremo, debil = nodos[-1], dom.elements[n // 2 + 1]
    return dom, extremo, debil


def _registro(dom, extremo, debil, curva: list, factor):
    """Callback que guarda (u, F, ε_w, elementos dañados) en cada paso."""
    a, b = debil.nodes

    def al_converger(paso, U, lam):
        danados = sum(1 for el in dom.elements.values() if el.state.vars[0]["damage"] > 0.0)
        eps_w = (U[b.dofs["ux"]] - U[a.dofs["ux"]]) / H
        curva.append((U[extremo.dofs["ux"]], factor(U, lam), eps_w, danados))
    return al_converger


def control_de_carga(L: float, F_max: float = 1.2 * F_PICO, pasos: int = 24):
    """Newton con control de carga hasta 1.2 veces la de pico."""
    dom, extremo, debil = barra(L)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[extremo.dofs["ux"]] = F_max
    solver = NonlinearSolver(Assembler(dom), num_steps=pasos, max_iter=30)
    curva: list = []
    try:
        solver.solve(F, step_callback=_registro(dom, extremo, debil, curva,
                                                lambda U, lam: lam * F_max))
        error = None
    except Exception as exc:                  # el solver no pasa el pico
        error = type(exc).__name__
    return np.array(curva), error


def control_de_desplazamiento(L: float, u_max: float, pasos: int):
    """Newton con el desplazamiento del extremo impuesto; la fuerza es la
    reacción en el extremo."""
    dom, extremo, debil = barra(L)
    extremo.fix_dof("ux", u_max)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    solver = NonlinearSolver(asm, num_steps=pasos, max_iter=30, adaptive=False)
    fila = extremo.dofs["ux"]
    curva: list = []
    try:
        solver.solve(np.zeros(dom.total_dofs), step_callback=_registro(
            dom, extremo, debil, curva,
            lambda U, lam: asm.assemble_non_linear_system(U)[1][fila]))
        error = None
    except Exception as exc:
        error = type(exc).__name__
    return np.array(curva), error


def control_indirecto(L: float, dlambda: float = DLAMBDA, pasos: int = 500):
    """Control indirecto sobre el alargamiento del elemento débil."""
    dom, extremo, debil = barra(L)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[extremo.dofs["ux"]] = F_PICO
    a, b = debil.nodes
    solver = IndirectDisplacementSolver(
        Assembler(dom), control=[(b, "ux", 1.0), (a, "ux", -1.0)],
        initial_dlambda=dlambda, max_lambda=1.2, max_steps=pasos, max_iter=30)
    curva: list = []
    solver.solve(F, step_callback=_registro(dom, extremo, debil, curva,
                                            lambda U, lam: lam * F_PICO))
    return np.array(curva)


def arco_cilindrico(L: float, dlambda: float = DLAMBDA, pasos: int = 500):
    """ArcLengthSolver (cilíndrico) con el mismo primer paso, sin crecimiento."""
    dom, extremo, debil = barra(L)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[extremo.dofs["ux"]] = F_PICO
    solver = ArcLengthSolver(Assembler(dom), initial_dlambda=dlambda, dl_max_factor=1.0,
                             max_lambda=1.2, max_steps=pasos, max_iter=30)
    curva: list = []
    try:
        solver.solve(F, step_callback=_registro(dom, extremo, debil, curva,
                                                lambda U, lam: lam * F_PICO))
        error = None
    except Exception as exc:
        error = type(exc).__name__
    return np.array(curva), error


def _errores(curva: np.ndarray, L: float) -> dict:
    """Error máximo frente a la solución exacta con el mismo ε_w."""
    u_ex, F_ex = exacta(curva[:, 2], L)
    return {"error_F": float(np.max(np.abs(curva[:, 1] - F_ex)) / F_PICO),
            "error_u": float(np.max(np.abs(curva[:, 0] - u_ex)) / (L * K0W))}


def calcular() -> dict:
    out = {"E_GPa": E / 1e9, "A_cm2": A * 1e4, "kappa_0": KAPPA_0, "alpha": ALPHA,
           "alpha_kappa_0w": ALPHA * K0W, "sigma_t_MPa": E * KAPPA_0 / 1e6,
           "sigma_t_debil_MPa": E * K0W / 1e6, "debil_pct": 100 * (1 - DEBIL),
           "h_m": H, "F_pico_kN": F_PICO / 1e3, "dlambda": DLAMBDA,
           "dlambda_pct": 100 * DLAMBDA, "dlambda_grande_pct": 100 * DLAMBDA_GRANDE}
    # Energía disipada por el elemento débil hasta la rotura (sin el tope):
    # A·h·(½·E·κ₀w² + E·κ₀w/α), igual en las dos barras.
    out["energia_J"] = A * H * (0.5 * E * K0W**2 + E * K0W / ALPHA)
    # Curva exacta hasta donde el daño alcanza el tope DAMAGE_MAX del material
    # (más allá queda una rigidez residual (1 − DAMAGE_MAX)·E que ningún
    # trazado del ejemplo recorre).
    eps = K0W * np.geomspace(1.0, 60.0, 200000)
    eps = eps[(K0W / eps) * np.exp(-ALPHA * (eps - K0W)) > 1.0 - DAMAGE_MAX]
    for nombre, L in BARRAS.items():
        n = round(L / H)
        u_pico = L * K0W
        u_ex, F_ex = exacta(eps, L)
        d = {"L_m": L, "n": n, "u_pico_um": 1e6 * u_pico,
             # Energía elástica almacenada en el pico, toda la barra a σ = E·κ₀w.
             "energia_elastica_pico_J": F_PICO**2 * L / (2 * E * A),
             "criterio": (n - 1) * ALPHA * K0W,
             "u_min_exacto_sobre_pico": float(u_ex.min() / u_pico)}
        curva, err = control_de_carga(L)
        d["carga"] = {"F_max_sobre_pico": float(curva[:, 1].max() / F_PICO), "excepcion": err}
        u_max = (6.0 if nombre == "corta" else 2.0) * u_pico
        curva, err = control_de_desplazamiento(L, u_max, pasos=600)
        d["desplazamiento"] = {**_errores(curva, L), "excepcion": err,
                               "u_final_sobre_pico": float(curva[-1, 0] / u_pico),
                               "u": curva[:, 0] / u_pico, "F": curva[:, 1] / F_PICO}
        curva = control_indirecto(L)
        post = curva[curva[:, 3] >= 1]
        d["indirecto"] = {**_errores(curva, L), "pasos": len(curva),
                          "danados": int(curva[:, 3].max()),
                          "u_min_sobre_pico": float(post[:, 0].min() / u_pico),
                          "F_final_sobre_pico": float(curva[-1, 1] / F_PICO),
                          "u": curva[:, 0] / u_pico, "F": curva[:, 1] / F_PICO}
        d["indirecto"]["error_u_min"] = abs(d["indirecto"]["u_min_sobre_pico"]
                                            - d["u_min_exacto_sobre_pico"])
        d["exacta"] = {"u": u_ex[::200] / u_pico, "F": F_ex[::200] / F_PICO}
        out[nombre] = d
    # El arco cilíndrico en la barra larga: desplazamiento mínimo alcanzado
    # tras el inicio del daño (1 si no llega a dañarse).
    L = BARRAS["larga"]
    curva, err = arco_cilindrico(L)
    post = curva[curva[:, 3] >= 1]
    out["larga"]["cilindrico"] = {
        "excepcion": err, "danados": int(curva[:, 3].max()),
        "u_min_sobre_pico": float(post[:, 0].min() / (L * K0W)) if len(post) else 1.0}
    # Advertencia: el mismo control indirecto con pasos del 10 % de la carga.
    L = BARRAS["larga"]
    curva = control_indirecto(L, dlambda=DLAMBDA_GRANDE, pasos=100)
    out["salto"] = {**_errores(curva, L), "danados": int(curva[:, 3].max()),
                    "u": curva[:, 0] / (L * K0W), "F": curva[:, 1] / F_PICO}
    # Y el control de desplazamiento de la barra corta con 12 pasos.
    L = BARRAS["corta"]
    curva, _ = control_de_desplazamiento(L, 6.0 * L * K0W, pasos=12)
    out["salto_desplazamiento"] = {"danados": int(curva[:, 3].max()), "pasos": 12}
    return out


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt

    fig, ejes = plt.subplots(1, 2, figsize=(7.2, 3.5), sharey=True)
    for ax, nombre in zip(ejes, BARRAS):
        d = res[nombre]
        ax.plot(d["exacta"]["u"], d["exacta"]["F"], color=COLOR_ANALITICO, lw=3.5,
                alpha=0.25, label="exacta")
        ax.plot(d["indirecto"]["u"], d["indirecto"]["F"], "-", color=COLOR_SECUNDARIO,
                lw=1.2, label="control indirecto")
        dd = d["desplazamiento"]
        ax.plot(dd["u"], dd["F"], "o", color=COLOR_FE, ms=2.5, markevery=6,
                label="control de desplazamiento")
        ax.axhline(res[nombre]["carga"]["F_max_sobre_pico"], color=COLOR_TERCIARIO,
                   lw=0.9, ls="--", label="límite del control de carga")
        ax.set_title(f"barra {nombre} (L = {d['L_m']:g} m, {d['n']} elementos)", fontsize=9)
        ax.set_xlabel(r"$u / u_{\mathrm{pico}}$")
        ax.set_xlim(0, 4.5 if nombre == "corta" else 1.6)
        ax.set_ylim(0, 1.08)
    ejes[0].set_ylabel(r"$F / F_{\mathrm{pico}}$")
    ejes[1].legend(loc="center right", fontsize=7.5)
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_curvas")

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    d = res["larga"]
    ax.plot(d["exacta"]["u"], d["exacta"]["F"], color=COLOR_ANALITICO, lw=3.5, alpha=0.25,
            label="exacta")
    ax.plot(d["indirecto"]["u"], d["indirecto"]["F"], "-", color=COLOR_SECUNDARIO, lw=1.0,
            label=f"pasos del {res['dlambda_pct']:.0f} %")
    s = res["salto"]
    ax.plot(s["u"], s["F"], "s-", color=COLOR_FE, ms=3, lw=0.8,
            label=f"pasos del {res['dlambda_grande_pct']:.0f} %")
    ax.set_xlim(0, 4.0)
    ax.set_ylim(0, 1.08)
    ax.set_xlabel(r"$u / u_{\mathrm{pico}}$")
    ax.set_ylabel(r"$F / F_{\mathrm{pico}}$")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_salto")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    print(f"F_pico = {res['F_pico_kN']:.1f} kN   energía disipada = {res['energia_J']:.3f} J")
    for nombre in BARRAS:
        d = res[nombre]
        print(f"{nombre}: criterio de retroceso {d['criterio']:.3f}  "
              f"carga hasta {d['carga']['F_max_sobre_pico']:.4f} F_pico ({d['carga']['excepcion']})  "
              f"desplazamiento hasta u = {d['desplazamiento']['u_final_sobre_pico']:.3f} u_pico "
              f"({d['desplazamiento']['excepcion']}), error F {d['desplazamiento']['error_F']:.1e}  "
              f"indirecto: error F {d['indirecto']['error_F']:.1e}, u_min {d['indirecto']['u_min_sobre_pico']:.4f} "
              f"(exacto {d['u_min_exacto_sobre_pico']:.4f}), {d['indirecto']['danados']} dañado(s)")
    print(f"pasos del {res['dlambda_grande_pct']:.0f} %: {res['salto']['danados']} elementos dañados, "
          f"error u {res['salto']['error_u']:.2f}; desplazamiento corta con 12 pasos: "
          f"{res['salto_desplazamiento']['danados']} dañados")
    return res


if __name__ == "__main__":
    main()
