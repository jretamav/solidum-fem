"""Ejemplo 1 del manual de ejemplos — voladizo de marco 2D.

Parte A (``modelo.yaml``): voladizo Euler-Bernoulli con carga puntual en el
extremo. Se comparan desplazamientos, giros, reacciones y fuerzas internas con
la solución analítica y se lee la convención de signos del API público
(Reglas.md §5): ``V > 0`` tiende a girar el diferencial en sentido horario,
``M > 0`` tracciona la fibra inferior.

Parte B (API de Python): cociente de flechas Timoshenko / Euler-Bernoulli en
función de la esbeltez ``L/h`` para una sección rectangular, contra
``1 + 3EI/(κGAL²)``.

Uso::

    python examples/voladizo_marco/run.py

escribe ``resultados.json`` y las figuras en esta carpeta.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))

from ejemplos_comun import (  # noqa: E402
    COLOR_ANALITICO, COLOR_FE, COLOR_SECUNDARIO, error_relativo,
    estilo_figuras, guardar_figura, guardar_resultados, resolver_yaml_estatico,
)

from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.frame import Frame2DEuler, Frame2DTimoshenko  # noqa: E402
from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import LinearSolver  # noqa: E402

# Datos de la parte A (deben coincidir con modelo.yaml).
L = 4.0
B = H = 0.1
E = 200.0e9
P = 1000.0
I = B * H**3 / 12.0

# Parte B: sección rectangular b × h, ν = 0.3, κ = 5/6.
NU = 0.3
KAPPA = 5.0 / 6.0
ESBELTECES = (1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 100, 1000)

# Tolerancias de lo que el manual afirma (las comprueba
# tests/test_examples_manual.py). Parte A: la hermítica cúbica es exacta con
# cargas nodales; parte B: la rigidez de Timoshenko con Φ es la exacta de la
# viga de Timoshenko, así que un elemento basta a cualquier esbeltez.
TOLERANCIAS = {
    "A.error_max": 1e-10,
    "B.error_max": 1e-10,
}


def _analitica_voladizo(x):
    """Flecha, giro, cortante y flector del voladizo con P hacia abajo en x = L.

    v(x) = −P x² (3L − x) / (6EI),  θ(x) = v'(x) = −P x (2L − x) / (2EI),
    V(x) = +P (convención de viga: el tramo derecho empuja hacia abajo, lo
    que tiende a girar el diferencial en sentido horario),
    M(x) = −P (L − x) (flector negativo: tracción en la fibra superior).
    """
    x = np.asarray(x, dtype=float)
    v = -P * x**2 * (3 * L - x) / (6 * E * I)
    th = -P * x * (2 * L - x) / (2 * E * I)
    V = np.full_like(x, P)
    M = -P * (L - x)
    return v, th, V, M


def _parte_a() -> dict:
    domain, _, res = resolver_yaml_estatico(AQUI / "modelo.yaml")
    nodos = sorted(domain.nodes.values(), key=lambda n: n.coordinates[0])
    x = np.array([n.coordinates[0] for n in nodos])
    v = np.array([res.U[n.dofs["uy"]] for n in nodos])
    th = np.array([res.U[n.dofs["rz"]] for n in nodos])
    v_a, th_a, _, _ = _analitica_voladizo(x)

    # Fuerzas internas en ambos extremos de cada elemento (ejes locales =
    # globales porque la viga es horizontal).
    x_ef, V_ef, M_ef = [], [], []
    for eid in sorted(res.element_forces):
        el = domain.elements[eid]
        xi, xj = el.nodes[0].coordinates[0], el.nodes[1].coordinates[0]
        f = res.element_forces[eid].components
        x_ef += [xi, xj]
        V_ef += list(f["V"])
        M_ef += list(f["M"])
    x_ef, V_ef, M_ef = map(np.array, (x_ef, V_ef, M_ef))
    _, _, V_a, M_a = _analitica_voladizo(x_ef)

    reac = res.reactions_by_node[1]
    delta = -v[-1]
    delta_a = P * L**3 / (3 * E * I)
    theta = th[-1]
    theta_a = -P * L**2 / (2 * E * I)

    errores = {
        "flecha": error_relativo(delta, delta_a),
        "giro": error_relativo(theta, theta_a),
        "Ry": error_relativo(reac["uy"], P),
        "Mz": error_relativo(reac["rz"], P * L),
        "M_empotramiento": error_relativo(M_ef[0], -P * L),
        "V": float(np.max(np.abs(V_ef - V_a)) / P),
        "M": float(np.max(np.abs(M_ef - M_a)) / (P * L)),
        "v_nodal": float(np.max(np.abs(v - v_a)) / delta_a),
    }
    return {
        "L": L, "b": B, "h": H, "E": E, "I": I, "P": P,
        "flecha_fe": delta, "flecha_analitica": delta_a,
        "flecha_fe_mm": 1e3 * delta, "flecha_analitica_mm": 1e3 * delta_a,
        "giro_fe": theta, "giro_analitico": theta_a,
        "reaccion_Ry": reac["uy"], "reaccion_Rx": reac["ux"], "reaccion_Mz": reac["rz"],
        "M_empotramiento": M_ef[0], "V_empotramiento": V_ef[0],
        "errores": errores,
        "error_max": max(errores.values()),
        "nodos_x": x, "nodos_v": v, "nodos_giro": th,
        "ef_x": x_ef, "ef_V": V_ef, "ef_M": M_ef,
    }


def _flecha_punta(cls, L_h: float) -> float:
    """Flecha en la punta de un voladizo de un solo elemento, sección b × h = 0.1 × 0.2 m."""
    h = 0.2
    b = 0.1
    Lb = L_h * h
    A = b * h
    Ib = b * h**3 / 12.0
    domain = Domain()
    mat = Elastic1D(E=E)
    n1 = domain.add_node(1, [0.0, 0.0])
    n2 = domain.add_node(2, [Lb, 0.0])
    if cls is Frame2DTimoshenko:
        el = Frame2DTimoshenko(1, [n1, n2], mat, A=A, I=Ib, As=KAPPA * A, nu=NU)
    else:
        el = Frame2DEuler(1, [n1, n2], mat, A=A, I=Ib)
    domain.add_element(el)
    for d in ("ux", "uy", "rz"):
        n1.fix_dof(d, 0.0)
    domain.generate_equation_numbers()
    F = np.zeros(domain.total_dofs)
    F[n2.dofs["uy"]] = -P
    U = LinearSolver(Assembler(domain)).solve(F)
    return -U[n2.dofs["uy"]]


def _parte_b() -> dict:
    lh = np.array(ESBELTECES, dtype=float)
    cociente = np.array([_flecha_punta(Frame2DTimoshenko, s) / _flecha_punta(Frame2DEuler, s)
                         for s in lh])
    # 3EI/(κGAL²) con I = bh³/12, A = bh, G = E/(2(1+ν)):  (1+ν)/(2κ) · (h/L)²
    c = (1 + NU) / (2 * KAPPA)
    analitico = 1 + c * lh**-2
    err = np.abs(cociente - analitico) / analitico
    return {
        "nu": NU, "kappa": KAPPA, "coef": c,
        "L_h": lh, "cociente_fe": cociente, "cociente_analitico": analitico,
        "cociente_Lh1": cociente[0], "cociente_Lh10": cociente[list(lh).index(10)],
        "aporte_cortante_Lh10_pct": 100 * (cociente[list(lh).index(10)] - 1),
        "aporte_cortante_Lh5_pct": 100 * (cociente[list(lh).index(5)] - 1),
        "exceso_Lh1000": cociente[list(lh).index(1000)] - 1,
        "error_max": float(err.max()),
    }


def _cota(error: float) -> int:
    """k tal que error < 10^-k. El capítulo cita la cota y no el error: al
    nivel del redondeo, la cifra exacta cambia entre ejecuciones y máquinas."""
    return int(np.floor(-np.log10(error)))


def calcular() -> dict:
    a, b = _parte_a(), _parte_b()
    a["cota_orden"] = _cota(a["error_max"])
    b["cota_orden"] = _cota(b["error_max"])
    return {"A": a, "B": b}


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, Rectangle

    a = res["A"]
    xs = np.linspace(0, L, 200)
    v_a, _, V_a, M_a = _analitica_voladizo(xs)

    # --- Esquema del problema
    fig, ax = plt.subplots(figsize=(6.4, 1.9))
    ax.plot([0, L], [0, 0], color=COLOR_ANALITICO, lw=4, solid_capstyle="butt")
    ax.add_patch(Rectangle((-0.18, -0.45), 0.18, 0.9, facecolor="none",
                           edgecolor=COLOR_ANALITICO, hatch="////", lw=1.0))
    ax.add_patch(FancyArrowPatch((L, 0.75), (L, 0.06), arrowstyle="-|>",
                                 mutation_scale=16, color=COLOR_SECUNDARIO, lw=1.8))
    ax.text(L + 0.08, 0.5, r"$P = 1$ kN", color=COLOR_SECUNDARIO, va="center")
    ax.annotate("", xy=(0, -0.62), xytext=(L, -0.62),
                arrowprops=dict(arrowstyle="<->", color="0.35", lw=0.9))
    ax.text(L / 2, -0.78, r"$L = 4$ m", ha="center", va="top", color="0.25")
    ax.text(0.25, 0.2, r"$E = 200$ GPa,  $b = h = 0.1$ m", color="0.25", fontsize=9)
    for xn in (0, 1, 2, 3, 4):
        ax.plot(xn, 0, "o", ms=5, mfc="white", mec=COLOR_FE, zorder=3)
    ax.set_xlim(-0.5, L + 1.1)
    ax.set_ylim(-1.05, 0.95)
    ax.set_aspect("equal")
    ax.axis("off")
    guardar_figura(fig, carpeta, "fig_esquema")

    # --- Deformada, cortante y flector
    fig, axs = plt.subplots(3, 1, figsize=(6.4, 5.4), sharex=True)
    axs[0].plot(xs, 1e3 * v_a, color=COLOR_ANALITICO, lw=1.3, label="analítica")
    axs[0].plot(a["nodos_x"], 1e3 * np.asarray(a["nodos_v"]), "o", color=COLOR_FE,
                label="Solidum (nodos)")
    axs[0].set_ylabel(r"$u_y$ [mm]")
    axs[0].legend(loc="lower left")
    axs[1].plot(xs, 1e-3 * V_a, color=COLOR_ANALITICO, lw=1.3)
    axs[1].plot(a["ef_x"], 1e-3 * np.asarray(a["ef_V"]), "s", color=COLOR_FE, ms=5)
    axs[1].set_ylabel(r"$V$ [kN]")
    axs[1].set_ylim(0, 1.5)
    axs[2].plot(xs, 1e-3 * M_a, color=COLOR_ANALITICO, lw=1.3)
    axs[2].plot(a["ef_x"], 1e-3 * np.asarray(a["ef_M"]), "s", color=COLOR_FE, ms=5)
    axs[2].axhline(0, color="0.5", lw=0.6)
    axs[2].set_ylabel(r"$M$ [kN·m]")
    axs[2].set_xlabel(r"$x$ [m]")
    fig.align_ylabels(axs)
    guardar_figura(fig, carpeta, "fig_diagramas")

    # --- Timoshenko / Euler-Bernoulli
    b = res["B"]
    fig, ax = plt.subplots()
    lh = np.geomspace(1, 1000, 300)
    ax.semilogx(lh, 1 + b["coef"] * lh**-2, color=COLOR_ANALITICO, lw=1.3,
                label=r"$1 + 3EI/(\kappa G A L^2)$")
    ax.semilogx(b["L_h"], b["cociente_fe"], "o", color=COLOR_FE,
                label="Solidum (1 elemento)")
    ax.axhline(1.0, color="0.5", lw=0.6)
    ax.set_xlabel(r"esbeltez $L/h$")
    ax.set_ylabel(r"$\delta_\mathrm{Timoshenko}\,/\,\delta_\mathrm{Euler}$")
    ax.legend()
    guardar_figura(fig, carpeta, "fig_timoshenko")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    print(f"Flecha: {res['A']['flecha_fe_mm']:.4f} mm  (analítica {res['A']['flecha_analitica_mm']:.4f} mm)")
    print(f"Error máximo parte A: {res['A']['error_max']:.2e}  parte B: {res['B']['error_max']:.2e}")
    return res


if __name__ == "__main__":
    main()
