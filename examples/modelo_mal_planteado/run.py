"""Ejemplo 11 del manual de ejemplos — un modelo mal planteado.

Una armadura Pratt de cuatro paneles, primero bien apoyada y después con tres
errores de modelado habituales: dos rodillos, un rodillo mal orientado y una
diagonal de menos (con y sin una carga que active el mecanismo). En cada caso
se compara lo que dice Solidum (red de seguridad del análisis estático, ADR
0019) con lo que dice la estática sola. La matriz de equilibrio de los nudos
no usa rigideces. Su rango da el número de mecanismos, ``2j − rango``, y el
núcleo de su traspuesta, que es la matriz de compatibilidad, da su forma.

Uso::

    python examples/modelo_mal_planteado/run.py
"""
from __future__ import annotations

import logging
import sys
import textwrap
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))

from ejemplos_comun import (  # noqa: E402
    COLOR_FE, COLOR_SECUNDARIO,
    estilo_figuras, guardar_figura, guardar_resultados,
)

from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.truss import Truss2D  # noqa: E402
from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.linalg.dispatcher import available_overrides  # noqa: E402
from solidum.math.solvers import LinearSolver, ModalSolver, NonlinearSolver  # noqa: E402
from solidum.math.solvers.diagnostics import IllPosedSystemError, MechanismError  # noqa: E402

PANEL, ALTO, N_PANELES = 3.0, 3.0, 4     # [m]
E, A, P = 200.0e9, 20.0e-4, 100.0e3      # acero, 20 cm², 100 kN
RHO = 7850.0                             # sólo para el análisis modal [kg/m³]

# Apoyos en los nudos inferiores extremos: (nudo, desplazamiento impedido).
APOYOS = {
    "correcto": [(0, "ux"), (0, "uy"), (-1, "uy")],        # pasador y rodillo vertical
    "dos_rodillos": [(0, "uy"), (-1, "uy")],
    "rodillo_girado": [(0, "ux"), (0, "uy"), (-1, "ux")],  # el rodillo impide u_x
}
# Cada caso: apoyos, si falta la diagonal del segundo panel, y la carga.
CASOS = {
    "correcto": ("correcto", False, "vertical"),
    "dos_rodillos": ("dos_rodillos", False, "vertical"),
    "rodillo_girado": ("rodillo_girado", False, "vertical"),
    "sin_diagonal": ("correcto", True, "vertical"),
    "sin_diagonal_H": ("correcto", True, "horizontal"),
}
MAL_PLANTEADOS = [c for c in CASOS if c != "correcto"]

TOLERANCIAS = {
    # El modelo bien apoyado coincide con la estática.
    "correcto.error_N": 1e-10,
    "correcto.error_flecha": 1e-10,
    "casos.correcto.rechazado": (0, 0),
    "casos.correcto.mecanismos": (0, 0),
    # Cada modelo mal planteado es un mecanismo, y Solidum lo rechaza.
    **{f"casos.{c}.mecanismos": (1, 1) for c in MAL_PLANTEADOS},
    **{f"casos.{c}.rechazado": (1, 1) for c in MAL_PLANTEADOS},
    # El rodillo girado cumple la cuenta de Maxwell; los otros no.
    "casos.rodillo_girado.deficit": (0, 0),
    "casos.dos_rodillos.deficit": (1, 1),
    "casos.sin_diagonal.deficit": (1, 1),
    # Con estas cargas hay equilibrio en unos casos y en otros no.
    "casos.correcto.con_equilibrio": (1, 1),
    "casos.dos_rodillos.con_equilibrio": (1, 1),
    "casos.rodillo_girado.con_equilibrio": (0, 0),
    "casos.sin_diagonal.con_equilibrio": (0, 0),
    "casos.sin_diagonal_H.con_equilibrio": (1, 1),
    # Las fuerzas extremas: la diagonal de apoyo, 1.5·√2·P, y la cuerda
    # superior central, −M/h = −2P.
    "correcto.error_N_max": 1e-10,
    "correcto.error_N_min": 1e-10,
    # El área en cm² en vez de m²: nada cambia salvo la flecha.
    "unidades.cambio_N": 1e-10,
    "unidades.razon_flecha": (0.99e-4, 1.01e-4),
    # Todos los solvers directos disponibles dan el mismo diagnóstico.
    "backends_discrepantes": 0,
    # El análisis modal localiza el mecanismo interno: una frecuencia nula,
    # con la forma del mecanismo de la estática.
    "modal.frecuencias_nulas": (1, 1),
    "modal.coseno_con_estatica": (0.999, 1.0 + 1e-9),
    # Fuera del alcance de la red: el iterativo y el no lineal aceptan el
    # mecanismo interno que la carga no activa.
    "limites.acepta_iterativo": (1, 1),
    "limites.acepta_no_lineal": (1, 1),
}


def geometria(sin_diagonal: bool = False):
    """Nudos y barras de la armadura Pratt: cuerdas, montantes y diagonales
    que bajan hacia el centro. Los nudos inferiores son ``0 … N_PANELES`` y
    los superiores les siguen en el mismo orden."""
    s = N_PANELES + 1
    nudos = [(PANEL * i, 0.0) for i in range(s)] + [(PANEL * i, ALTO) for i in range(s)]
    barras = ([(i, i + 1) for i in range(N_PANELES)]              # cuerda inferior
              + [(s + i, s + i + 1) for i in range(N_PANELES)]    # cuerda superior
              + [(i, s + i) for i in range(s)])                   # montantes
    for i in range(N_PANELES):
        if sin_diagonal and i == 1:
            continue                                              # falta esta diagonal
        barras.append((s + i, i + 1) if i < N_PANELES // 2 else (i, s + i + 1))
    return np.array(nudos), barras


def armadura(caso: str, area: float = A, densidad: float | None = None):
    """Modelo de Solidum del caso: dominio, nudos y vector de cargas."""
    apoyos, sin_diagonal, carga = CASOS[caso]
    xy, barras = geometria(sin_diagonal)
    dom = Domain()
    acero = Elastic1D(E=E, density=densidad)
    nodos = [dom.add_node(k + 1, list(p)) for k, p in enumerate(xy)]
    for a, b in barras:
        dom.add_element(Truss2D(len(dom.elements) + 1, [nodos[a], nodos[b]], acero, A=area))
    inferiores = nodos[:N_PANELES + 1]
    for k, dof in APOYOS[apoyos]:
        inferiores[k].fix_dof(dof, 0.0)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    if carga == "vertical":                    # P hacia abajo en los nudos inferiores interiores
        for nodo in inferiores[1:-1]:
            F[nodo.dofs["uy"]] = -P
    else:                                      # P horizontal en el rodillo
        F[inferiores[-1].dofs["ux"]] = P
    return dom, nodos, F


def analizar(caso: str, area: float = A, linear_algebra: str = "auto") -> dict:
    """Resuelve el caso. Si el modelo está mal planteado, Solidum lo rechaza
    con una excepción que explica el problema; aquí se captura para verla."""
    dom, nodos, F = armadura(caso, area)
    try:
        U = LinearSolver(Assembler(dom), linear_algebra=linear_algebra).solve(F)
    except MechanismError as e:            # antes de resolver: mecanismo rígido
        return {"resultado": "MechanismError", "movimiento": e.motions[0],
                "n_movimientos": len(e.motions), "mensaje": str(e)}
    except IllPosedSystemError as e:       # después de factorizar: mecanismo interno
        return {"resultado": "IllPosedSystemError", "pivotes_nulos": e.zero_pivots or 0,
                "mensaje": str(e)}
    N = [el.internal_forces(U).components["N"][0] for el in dom.elements.values()]
    return {"resultado": "solución", "N": np.array(N),
            "flecha": -U[nodos[N_PANELES // 2].dofs["uy"]]}


def buscar_mecanismo(caso: str, n_modos: int = 4):
    """Análisis modal del mismo modelo, con densidad: cada frecuencia nula es
    un mecanismo, y su forma modal muestra qué parte se mueve sin deformarse."""
    dom, nodos, _ = armadura(caso, densidad=RHO)
    modal = ModalSolver(Assembler(dom), n_modes=n_modos).solve()
    f = modal.frequencies_hz
    nulas = f < 1e-6 * f.max()
    modos = [np.array([modal.modes[nodo.dofs[d], k] for nodo in nodos for d in ("ux", "uy")])
             for k in np.flatnonzero(nulas)]
    return f, modos


def no_lineal(caso: str) -> str:
    """El mismo caso con ``NonlinearSolver``: ¿lo acepta?"""
    dom, _, F = armadura(caso)
    try:
        NonlinearSolver(Assembler(dom)).solve(F)
    except Exception as e:                 # noqa: BLE001 — sólo se registra el tipo
        return type(e).__name__
    return "solución"


def estatica(caso: str) -> dict:
    """Método de los nudos en forma matricial, sin rigideces: ``A·s = −F``,
    con ``s`` las fuerzas axiales (tracción positiva) y las reacciones."""
    apoyos, sin_diagonal, carga = CASOS[caso]
    xy, barras = geometria(sin_diagonal)
    n2 = 2 * len(xy)
    columnas = []
    for a, b in barras:
        e = (xy[b] - xy[a]) / np.linalg.norm(xy[b] - xy[a])
        c = np.zeros(n2)
        c[2 * a:2 * a + 2], c[2 * b:2 * b + 2] = e, -e     # la barra tira de sus dos nudos
        columnas.append(c)
    for k, dof in APOYOS[apoyos]:
        c = np.zeros(n2)
        c[2 * (k % (N_PANELES + 1)) + (dof == "uy")] = 1.0
        columnas.append(c)
    A_eq = np.column_stack(columnas)

    def cargas(nudos_y_valores):
        F = np.zeros(n2)
        for i, fx, fy in nudos_y_valores:
            F[2 * i], F[2 * i + 1] = fx, fy
        return F

    if carga == "vertical":
        F = cargas([(i, 0.0, -P) for i in range(1, N_PANELES)])
    else:
        F = cargas([(N_PANELES, P, 0.0)])
    s = np.linalg.lstsq(A_eq, -F, rcond=None)[0]
    # Núcleo de la compatibilidad Aᵀ: desplazamientos sin alargar ninguna
    # barra ni mover ningún apoyo, es decir, los mecanismos.
    _, sv, Vt = np.linalg.svd(A_eq.T)
    rango = int(np.sum(sv > 1e-10 * sv[0]))
    modos = Vt[rango:]
    # Trabajos virtuales: fuerzas n de una carga unidad en el centro.
    n_u = np.linalg.lstsq(A_eq, -cargas([(N_PANELES // 2, 0.0, -1.0)]), rcond=None)[0]
    L = np.array([np.linalg.norm(xy[b] - xy[a]) for a, b in barras])
    m = len(barras)
    return {
        "barras": m, "reacciones": len(APOYOS[apoyos]), "dos_j": n2,
        "deficit": n2 - m - len(APOYOS[apoyos]), "mecanismos": n2 - rango,
        "con_equilibrio": int(np.linalg.norm(A_eq @ s + F) <= 1e-10 * np.linalg.norm(F)),
        "N": s[:m], "flecha": float(np.sum(s[:m] * n_u[:m] * L) / (E * A)),
        "modo": modos[0] if len(modos) else np.zeros(n2), "xy": xy, "conexiones": barras,
    }


def calcular() -> dict:
    logging.disable(logging.WARNING)      # los rechazos también se registran en el log
    try:
        casos, detalle = {}, {}
        for caso in CASOS:
            ref, fe = estatica(caso), analizar(caso)
            casos[caso] = {
                "dos_j": ref["dos_j"], "m_mas_r": ref["barras"] + ref["reacciones"],
                "deficit": ref["deficit"], "mecanismos": ref["mecanismos"],
                "con_equilibrio": ref["con_equilibrio"],
                "equilibrio_txt": "sí" if ref["con_equilibrio"] else "no",
                "mensaje": _partir(fe.get("mensaje", "")),
                "rechazado": int(fe["resultado"] != "solución"), "resultado": fe["resultado"],
                "movimiento": fe.get("movimiento", ""), "pivotes_nulos": fe.get("pivotes_nulos", 0)}
            detalle[caso] = (ref, fe)
        ref, fe = detalle["correcto"]
        escala = np.abs(ref["N"]).max()
        correcto = {
            "error_N": float(np.abs(fe["N"] - ref["N"]).max() / escala),
            "error_flecha": abs(fe["flecha"] / ref["flecha"] - 1),
            "flecha_mm": 1e3 * fe["flecha"], "N_max_kN": 1e-3 * fe["N"].max(),
            "N_min_kN": 1e-3 * fe["N"].min(), "N": fe["N"],
            "error_N_max": abs(fe["N"].max() / (1.5 * np.sqrt(2) * P) - 1),
            "error_N_min": abs(fe["N"].min() / (-2 * P) - 1)}
        mal = analizar("correcto", area=A * 1e4)          # 20 "m²" en vez de 20 cm²
        unidades = {"cambio_N": float(np.abs(mal["N"] - fe["N"]).max() / escala),
                    "razon_flecha": mal["flecha"] / fe["flecha"],
                    "flecha_um": 1e6 * mal["flecha"]}
        f, modos_nulos = buscar_mecanismo("sin_diagonal")
        d = detalle["sin_diagonal"][0]["modo"]
        modal = {"frecuencias_nulas": len(modos_nulos), "f1_hz": float(f[len(modos_nulos)]),
                 "coseno_con_estatica": float(abs(modos_nulos[0] @ d)
                                              / (np.linalg.norm(modos_nulos[0]) * np.linalg.norm(d)))}
        limites = {
            "acepta_iterativo": int(analizar("sin_diagonal_H", linear_algebra="iterative")["resultado"]
                                    == "solución"),
            "acepta_no_lineal": int(no_lineal("sin_diagonal_H") == "solución")}
        # El diagnóstico no depende del solver directo que elija el despachador.
        directos = [b for b in ("lu", "pardiso") if b in available_overrides()]
        discrepantes = sum(
            analizar(c, linear_algebra=b)["resultado"] != casos[c]["resultado"]
            for b in directos for c in CASOS)
    finally:
        logging.disable(logging.NOTSET)
    modos = {c: {"d": detalle[c][0]["modo"]} for c in CASOS}
    return {"P_kN": P / 1e3, "panel_m": PANEL, "alto_m": ALTO, "E_GPa": E / 1e9,
            "A_cm2": A * 1e4, "n_nudos": casos["correcto"]["dos_j"] // 2,
            "n_barras": casos["correcto"]["m_mas_r"] - len(APOYOS["correcto"]),
            "casos": casos, "correcto": correcto, "unidades": unidades,
            "backends_discrepantes": discrepantes, "modos": modos, "modal": modal,
            "limites": limites,
            "xy": detalle["correcto"][0]["xy"],
            "conexiones": [list(b) for b in detalle["correcto"][0]["conexiones"]],
            "conexiones_sin_diagonal": [list(b) for b in detalle["sin_diagonal"][0]["conexiones"]]}


def _partir(mensaje: str, ancho: int = 84) -> str:
    """El mensaje de error en líneas cortas, para el listado del manual."""
    return "\n".join(textwrap.fill(linea, ancho, subsequent_indent="    ")
                     for linea in mensaje.splitlines())


def _apoyos(ax, xy, tipo):
    """Pasador: triángulo lleno; rodillo: círculo del lado que impide."""
    ultimo = N_PANELES
    dofs = {}
    for k, dof in APOYOS[tipo]:
        dofs.setdefault(k % (ultimo + 1), set()).add(dof)
    for k, d in dofs.items():
        x, y = xy[k]
        if d == {"ux", "uy"}:
            ax.plot(x, y - 0.35, marker="^", ms=8, color="k")
        elif d == {"uy"}:
            ax.plot(x, y - 0.35, marker="o", ms=7, mfc="white", mec="k")
        else:
            ax.plot(x + 0.4, y, marker="o", ms=7, mfc="white", mec="k")


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt

    xy = np.asarray(res["xy"])
    fig, ejes = plt.subplots(2, 2, figsize=(7.2, 4.4))
    paneles = (("correcto", "(a) bien apoyada"), ("dos_rodillos", "(b) dos rodillos"),
               ("rodillo_girado", "(c) rodillo horizontal"), ("sin_diagonal", "(d) sin diagonal"))
    for ax, (caso, titulo) in zip(ejes.flat, paneles):
        conexiones = res["conexiones_sin_diagonal"] if caso == "sin_diagonal" else res["conexiones"]
        if caso == "correcto":
            N = np.asarray(res["correcto"]["N"])
            for (a, b), n in zip(conexiones, N):
                color = COLOR_FE if n > 1e-6 * abs(N).max() else (
                    COLOR_SECUNDARIO if n < -1e-6 * abs(N).max() else "0.7")
                ax.plot(*xy[[a, b]].T, color=color, lw=0.8 + 2.5 * abs(n) / abs(N).max())
            ax.plot([], [], color=COLOR_FE, lw=2, label="tracción")
            ax.plot([], [], color=COLOR_SECUNDARIO, lw=2, label="compresión")
            ax.legend(fontsize=7, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02), frameon=False)
        else:
            d = np.asarray(res["modos"][caso]["d"]).reshape(-1, 2)
            d = d * 0.2 * PANEL / np.abs(d).max()
            if d.flat[np.argmax(np.abs(d))] < 0:
                d = -d                                    # signo del mecanismo: arbitrario
            for a, b in conexiones:
                ax.plot(*xy[[a, b]].T, color="0.75", lw=0.8)
                ax.plot(*(xy + d)[[a, b]].T, color=COLOR_SECUNDARIO, lw=1.4)
        _apoyos(ax, xy, CASOS[caso][0])
        ax.set_title(titulo, fontsize=9, loc="left")
        ax.set_aspect("equal")
        ax.set_xlim(-1.0, N_PANELES * PANEL + 1.3)
        ax.set_ylim(-1.3, ALTO + 1.3)
        ax.axis("off")
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_armadura")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    for caso, d in res["casos"].items():
        print(f"{caso:15s} 2j = {d['dos_j']}, m + r = {d['m_mas_r']}, mecanismos = {d['mecanismos']}, "
              f"equilibrio = {d['con_equilibrio']} → {d['resultado']} {d['movimiento']}"
              f"{' pivotes nulos: ' + str(d['pivotes_nulos']) if d['pivotes_nulos'] else ''}")
    c = res["correcto"]
    print(f"correcto: error N = {c['error_N']:.1e}, flecha = {c['flecha_mm']:.3f} mm "
          f"(error {c['error_flecha']:.1e}), N en [{c['N_min_kN']:.1f}, {c['N_max_kN']:.1f}] kN")
    u = res["unidades"]
    print(f"área ×10⁴: cambio de N = {u['cambio_N']:.1e}, flecha ×{u['razon_flecha']:.2e}")
    print(f"backends discrepantes: {res['backends_discrepantes']}")
    return res


if __name__ == "__main__":
    main()
