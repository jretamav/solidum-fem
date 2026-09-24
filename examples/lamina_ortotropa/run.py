"""Ejemplo 10 del manual de ejemplos — lámina ortótropa fuera de ejes.

Una lámina de carbono-epoxi con las fibras a un ángulo θ del eje de carga.
Bajo tracción uniaxial libre el estado es uniforme y la teoría clásica de
láminas (Jones, *Mechanics of Composite Materials*, cap. 2) da en forma
cerrada el módulo aparente, el Poisson aparente y el coeficiente de
influencia mutua η_xy,x, el acoplamiento tracción-cortante que el material
isótropo no tiene. Con la flexibilidad rotada S̄ = T⁻¹·S·T⁻ᵀ (T, la matriz
de transformación de las deformaciones de ingeniería):

    E_x = 1/S̄₁₁,   ν_xy = −S̄₁₂/S̄₁₁,   η_xy,x = S̄₁₆/S̄₁₁

Después se sujeta la misma tira con mordazas rígidas, que impiden la
distorsión: la tira se deforma en S y el módulo que se mediría crece.

Uso::

    python examples/lamina_ortotropa/run.py
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
from solidum.elements.solid_2d import Quad8  # noqa: E402
from solidum.materials.orthotropic_2d import Orthotropic2D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import LinearSolver  # noqa: E402

# Carbono-epoxi unidireccional (valores típicos de un T300/5208).
E1, E2, G12, NU12 = 181.0e9, 10.3e9, 7.17e9, 0.28
W = 0.02                   # ancho de la probeta [m]
SIGMA = 10.0e6             # tracción aplicada en la tira libre [Pa]
ANGULOS = (0.0, 10.0, 20.0, 30.0, 45.0, 60.0, 75.0, 90.0)
ESBELTECES = (2.0, 4.0, 8.0)   # L/w de la tira con mordazas

TOLERANCIAS = {
    # Tira libre: estado uniforme, el Quad8 es exacto.
    **{f"libre.{int(a)}.error_E": 1e-9 for a in ANGULOS},
    **{f"libre.{int(a)}.error_nu": 1e-9 for a in ANGULOS},
    **{f"libre.{int(a)}.error_eta": 1e-9 for a in ANGULOS},
    # A 45° la rigidez cae por debajo del 10 % de E₁ (y no a un tercio).
    "E45_sobre_E1": (0.05, 0.10),
    # Mordazas: más rígida cuanto más corta, y más a 30° que a 45°.
    "mordazas.45.2.E_sobre_Ex": (1.1, 1.3),
    "mordazas.45.8.E_sobre_Ex": (1.0, 1.03),
    "mordazas.30.2.E_sobre_Ex": (1.4, 1.8),
    # Teoremas de energía: E_x ≤ E_aparente ≤ Q̄₁₁, es decir, la fracción del
    # rigidizamiento máximo posible queda entre 0 y 1.
    **{f"mordazas.{t}.{int(r)}.fraccion_de_cota": (0.0, 1.0) for t in (30, 45) for r in ESBELTECES},
    # Las esquinas sujetas concentran esfuerzos y convergen despacio: la malla
    # de 16 elementos por ancho cambia a lo sumo un 1 % respecto de la de 8.
    **{f"mordazas.{t}.{int(r)}.cambio_malla": 0.012 for t in (30, 45) for r in ESBELTECES},
}


def flexibilidad_rotada(theta_deg: float) -> np.ndarray:
    """S̄ en ejes de carga (Voigt de ingeniería), independiente de Solidum."""
    S = np.array([[1 / E1, -NU12 / E1, 0.0], [-NU12 / E1, 1 / E2, 0.0], [0.0, 0.0, 1 / G12]])
    t = np.radians(theta_deg)
    c, s = np.cos(t), np.sin(t)
    T = np.array([[c * c, s * s, c * s], [s * s, c * c, -c * s],
                  [-2 * c * s, 2 * c * s, c * c - s * s]])     # ε₁₂ = T·ε_xy
    Ti = np.linalg.inv(T)
    return Ti @ S @ Ti.T


def jones(theta_deg: float) -> dict:
    """Constantes aparentes de la lámina (Jones, cap. 2). E_x se comprueba
    también con la fórmula explícita en cosenos y senos."""
    Sb = flexibilidad_rotada(theta_deg)
    t = np.radians(theta_deg)
    c, s = np.cos(t), np.sin(t)
    E_x = 1 / (c**4 / E1 + (1 / G12 - 2 * NU12 / E1) * c**2 * s**2 + s**4 / E2)
    assert abs(E_x * Sb[0, 0] - 1) < 1e-12        # las dos vías coinciden
    return {"E": E_x, "nu": -Sb[0, 1] / Sb[0, 0], "eta": Sb[0, 2] / Sb[0, 0]}


def tira(theta: float, L: float, nx: int, ny: int):
    """Tira ``L × W`` de ``nx × ny`` elementos ``Quad8`` con las fibras a
    ``theta`` grados del eje ``x``."""
    dom = Domain()
    lamina = Orthotropic2D(E1=E1, E2=E2, G12=G12, nu12=NU12, theta=theta)
    nodos = {}
    for j in range(2 * ny + 1):
        for i in range(2 * nx + 1):
            if i % 2 and j % 2:
                continue                       # centro: no es nodo del Quad8
            nodos[(i, j)] = dom.add_node(len(nodos) + 1, [L * i / (2 * nx), W * j / (2 * ny)])
    elementos = []
    for j in range(ny):
        for i in range(nx):
            a, b = 2 * i, 2 * j
            conectividad = [nodos[(a, b)], nodos[(a + 2, b)], nodos[(a + 2, b + 2)], nodos[(a, b + 2)],
                            nodos[(a + 1, b)], nodos[(a + 2, b + 1)], nodos[(a + 1, b + 2)], nodos[(a, b + 1)]]
            el = Quad8(len(dom.elements) + 1, conectividad, lamina)
            dom.add_element(el)
            elementos.append((el, i, j))
    return dom, nodos, elementos


def traccion_libre(theta: float, L: float = 2 * W, nx: int = 4, ny: int = 2) -> dict:
    """Tracción ``SIGMA`` en el borde derecho; el izquierdo sólo impide ``u_x``
    (y ``u_y`` en una esquina), así que la tira puede estrecharse y
    distorsionarse libremente: el estado es uniforme."""
    dom, nodos, elementos = tira(theta, L, nx, ny)
    for (i, j), nodo in nodos.items():
        if i == 0:
            nodo.fix_dof("ux", 0.0)
    nodos[(0, 0)].fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    for el, i, j in elementos:
        if i == nx - 1:                       # arista derecha del Quad8: índice 1
            f = el.compute_edge_traction(1, np.array([SIGMA, 0.0]))
            F[el.get_global_dof_indices()] += f
    U = LinearSolver(Assembler(dom)).solve(F)

    def u(i, j, d):
        return U[nodos[(i, j)].dofs[d]]
    eps_x = u(2 * nx, 0, "ux") / L
    eps_y = u(0, 2 * ny, "uy") / W
    gamma = u(2 * nx, 0, "uy") / L            # u_x = 0 en x = 0 ⇒ γ = ∂u_y/∂x
    return {"E": SIGMA / eps_x, "nu": -eps_y / eps_x, "eta": gamma / eps_x,
            "deformada": _deformada(nodos, U, nx, ny)}


def traccion_mordazas(theta: float, esbeltez: float, nx_por_w: int = 16, ny: int = 16) -> dict:
    """Mordazas rígidas: en los dos extremos ``u_y = 0`` y ``u_x`` uniforme
    (0 y δ). El módulo aparente es la fuerza de reacción entre ``w·δ/L``."""
    L = esbeltez * W
    nx = int(nx_por_w * esbeltez)
    dom, nodos, _ = tira(theta, L, nx, ny)
    delta = 1.0e-4 * L
    derecha = []
    for (i, j), nodo in nodos.items():
        if i == 0:
            nodo.fix_dof("ux", 0.0)
            nodo.fix_dof("uy", 0.0)
        elif i == 2 * nx:
            nodo.fix_dof("ux", delta)
            nodo.fix_dof("uy", 0.0)
            derecha.append(nodo)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    U = LinearSolver(asm).solve(np.zeros(dom.total_dofs))
    _, F_int = asm.assemble_non_linear_system(U)
    P = sum(F_int[nodo.dofs["ux"]] for nodo in derecha)
    return {"E": (P / W) / (delta / L), "deformada": _deformada(nodos, U, nx, ny)}


def _deformada(nodos, U, nx, ny) -> dict:
    """Contorno de la tira (sin y con desplazamientos), para la figura."""
    borde = ([(i, 0) for i in range(2 * nx + 1)] + [(2 * nx, j) for j in range(1, 2 * ny + 1)]
             + [(i, 2 * ny) for i in range(2 * nx - 1, -1, -1)] + [(0, j) for j in range(2 * ny - 1, -1, -1)])
    xy = np.array([nodos[k].coordinates for k in borde])
    uv = np.array([[U[nodos[k].dofs["ux"]], U[nodos[k].dofs["uy"]]] for k in borde])
    return {"x": xy[:, 0], "y": xy[:, 1], "ux": uv[:, 0], "uy": uv[:, 1]}


def calcular() -> dict:
    out = {"E1_GPa": E1 / 1e9, "E2_GPa": E2 / 1e9, "G12_GPa": G12 / 1e9, "nu12": NU12,
           "nu21": NU12 * E2 / E1, "W_mm": 1e3 * W, "sigma_MPa": SIGMA / 1e6, "libre": {}, "mordazas": {}}
    for a in ANGULOS:
        ref, fe = jones(a), traccion_libre(a)
        out["libre"][str(int(a))] = {
            "E_GPa": fe["E"] / 1e9, "E_jones_GPa": ref["E"] / 1e9, "nu": fe["nu"],
            "nu_jones": ref["nu"], "eta": fe["eta"], "eta_jones": ref["eta"],
            "error_E": abs(fe["E"] / ref["E"] - 1), "error_nu": abs(fe["nu"] - ref["nu"]),
            "error_eta": abs(fe["eta"] - ref["eta"])}
    out["E45_sobre_E1"] = jones(45.0)["E"] / E1
    out["E45_pct"] = 100 * out["E45_sobre_E1"]
    out["deformada_libre"] = traccion_libre(45.0, L=4 * W, nx=8, ny=2)["deformada"]
    out["cota"] = {}
    for theta in (30, 45):
        out["mordazas"][str(theta)] = {}
        E_x = jones(theta)["E"]
        # Cota superior: la tira entera sin alargamiento transversal ni
        # distorsión, E = Q̄₁₁ (rigidez reducida de tensión plana, Q̄ = S̄⁻¹).
        cota = np.linalg.inv(flexibilidad_rotada(theta))[0, 0] / E_x
        out["cota"][str(theta)] = cota
        for r in ESBELTECES:
            gruesa = traccion_mordazas(theta, r, nx_por_w=8, ny=8)
            fina = traccion_mordazas(theta, r)
            out["mordazas"][str(theta)][str(int(r))] = {
                "E_sobre_Ex": fina["E"] / E_x, "exceso_pct": 100 * (fina["E"] / E_x - 1),
                "fraccion_de_cota": (fina["E"] / E_x - 1) / (cota - 1),
                "cambio_malla": abs(fina["E"] / gruesa["E"] - 1)}
            if theta == 45 and r == 4.0:
                out["deformada_mordazas"] = fina["deformada"]
    # Curva continua de Jones para la figura.
    th = np.linspace(0, 90, 181)
    out["curva"] = {"theta": th, "E_sobre_E1": [jones(t)["E"] / E1 for t in th],
                    "eta": [jones(t)["eta"] for t in th]}
    return out


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.3))
    c = res["curva"]
    ax.plot(c["theta"], c["E_sobre_E1"], color=COLOR_ANALITICO, lw=1.2, label=r"$E_x/E_1$ (Jones)")
    ax2 = ax.twinx()
    ax2.plot(c["theta"], c["eta"], color=COLOR_SECUNDARIO, lw=1.2, ls="--", label=r"$\eta_{xy,x}$ (Jones)")
    for a, d in res["libre"].items():
        ax.plot(float(a), d["E_GPa"] / res["E1_GPa"], "o", color=COLOR_FE, ms=4)
        ax2.plot(float(a), d["eta"], "s", color=COLOR_SECUNDARIO, ms=3.5)
    ax.set_xlabel(r"ángulo de las fibras $\theta$ [°]")
    ax.set_ylabel(r"$E_x / E_1$")
    ax2.set_ylabel(r"$\eta_{xy,x}$", color=COLOR_SECUNDARIO)
    ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
    lineas = ax.get_legend_handles_labels()[0] + ax2.get_legend_handles_labels()[0]
    ax.legend(lineas, [l.get_label() for l in lineas], fontsize=7.5, loc="center right")
    ax2.grid(False)
    ax2.spines["right"].set_visible(True)

    # Las dos deformadas se amplifican hasta el mismo alargamiento aparente
    # (un 6 % de la longitud): lo que se compara es la forma, no la magnitud.
    w = res["W_mm"] / 1e3
    for clave, color, etiqueta, desplazar in (("deformada_libre", COLOR_FE, "libre", 0.0),
                                              ("deformada_mordazas", COLOR_SECUNDARIO, "con mordazas", -1.7)):
        d = res[clave]
        x0, y0 = np.asarray(d["x"]) / w, np.asarray(d["y"]) / w + desplazar
        ux, uy = np.asarray(d["ux"]) / w, np.asarray(d["uy"]) / w
        escala = 0.06 * x0.max() / ux.max()
        bx.plot(x0, y0, color="0.7", lw=0.8)
        bx.plot(x0 + escala * ux, y0 + escala * uy, color=color, lw=1.4, label=etiqueta)
    bx.set_aspect("equal")
    bx.set_xlabel(r"$x / w$")
    bx.set_yticks([])
    bx.legend(fontsize=7.5, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.18))
    bx.spines["left"].set_visible(False)
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_lamina")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    print(f"ν21 = {res['nu21']:.4f};  E_x(45°)/E1 = {res['E45_sobre_E1']:.4f}")
    for a, d in res["libre"].items():
        print(f"θ = {a:>2}°: E_x = {d['E_GPa']:.3f} GPa ({d['E_jones_GPa']:.3f}), ν = {d['nu']:.4f} "
              f"({d['nu_jones']:.4f}), η = {d['eta']:+.4f} ({d['eta_jones']:+.4f})")
    for t, filas in res["mordazas"].items():
        for r, d in filas.items():
            print(f"mordazas θ = {t}°, L/w = {r}: E/E_x = {d['E_sobre_Ex']:.4f} (cota {res['cota'][t]:.3f}), "
                  f"cambio de malla {d['cambio_malla']:.1e}")
    return res


if __name__ == "__main__":
    main()
