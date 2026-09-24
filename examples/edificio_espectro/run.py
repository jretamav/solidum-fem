"""Ejemplo 7 del manual de ejemplos — análisis modal y espectro de respuesta de un edificio.

Marco de cinco niveles modelado con ``Frame2DEuler``: columnas sin masa y
losas que llevan la masa de cada piso. Con losas rígidas y columnas
inextensibles el marco es un **edificio de cortante**: una cadena de masas y
resortes iguales, empotrada en la base, con solución cerrada:

    k = 2·12·E·I/h³ por piso (dos columnas biempotradas)
    ω_j = 2·√(k/m)·sin[(2j − 1)·π / (2(2n + 1))]
    φ_j(i) = sin[(2j − 1)·i·π / (2n + 1)],   i = 1..n

El espectro de respuesta se combina además a mano con esos modos analíticos,
independientemente de Solidum, para comparar desplazamientos, masas
efectivas y cortante basal.

Uso::

    python examples/edificio_espectro/run.py
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
from solidum.elements.frame import Frame2DEuler  # noqa: E402
from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import ModalSolver, ResponseSpectrumSolver  # noqa: E402

N = 5                     # niveles
H, CLARO = 3.0, 6.0       # altura de entrepiso y claro [m]
E = 25.0e9                # concreto [Pa]
LADO = 0.40               # columnas cuadradas de 40 cm
A_COL, I_COL = LADO**2, LADO**4 / 12
M_PISO = 30.0e3           # masa de cada nivel [kg]
RIGIDEZ = 1.0e4           # losas y columnas axiales, relativas a la columna
G = 9.81

K_PISO = 2 * 12 * E * I_COL / H**3       # rigidez lateral de entrepiso

# Espectro de diseño simplificado (aceleración espectral, ξ = 5 %): rampa
# hasta la meseta, meseta y caída 1/T.
A0 = 0.30 * G                            # aceleración del terreno
T_B, T_C = 0.10, 0.50                    # inicio y fin de la meseta [s]


def Sa(T: float) -> float:
    """Aceleración espectral [m/s²] del espectro de diseño simplificado."""
    if T < T_B:
        return A0 * (1 + 1.5 * T / T_B)
    if T <= T_C:
        return 2.5 * A0
    return 2.5 * A0 * T_C / T


TOLERANCIAS = {
    # Frecuencias frente a la solución cerrada: la diferencia es la rigidez
    # finita de losas y columnas, y baja con ella.
    "error_w_max": 5e-4,
    "error_w_max_rigido": 5e-5,
    # Masas efectivas, desplazamientos SRSS y cortante basal frente a la
    # combinación hecha a mano con los modos analíticos.
    "error_masa_efectiva": 1e-3,
    "error_u_srss": 1e-3,
    "error_cortante": 1e-3,
    # Combinar desplazamientos y luego restar subestima la deriva del último
    # entrepiso (un 5 %).
    "deriva_mal_sobre_bien": (0.90, 0.97),
    # Formas modales frente a las analíticas.
    "error_formas": 1e-3,
    # CQC y SRSS casi coinciden: los modos están bien separados.
    "cqc_sobre_srss_max": (0.99, 1.02),
}


def marco(rigidez: float = RIGIDEZ):
    """Marco de ``N`` niveles y un vano. Columnas ``Frame2DEuler`` sin masa;
    cada losa es una viga ``rigidez`` veces más rígida que una columna, con
    la masa del piso. Las columnas tienen área ``rigidez`` veces la real
    (inextensibles). Devuelve el dominio y los nodos por ``(nivel, lado)``."""
    dom = Domain()
    columna = Elastic1D(E=E, density=0.0)
    A_losa = rigidez * A_COL
    losa = Elastic1D(E=E, density=M_PISO / (A_losa * CLARO))
    nodos = {}
    for i in range(N + 1):
        for lado, x in enumerate((0.0, CLARO)):
            nodos[(i, lado)] = dom.add_node(len(nodos) + 1, [x, i * H])
    for i in range(N):
        for lado in (0, 1):
            dom.add_element(Frame2DEuler(len(dom.elements) + 1,
                                         [nodos[(i, lado)], nodos[(i + 1, lado)]],
                                         columna, A=rigidez * A_COL, I=I_COL))
        dom.add_element(Frame2DEuler(len(dom.elements) + 1,
                                     [nodos[(i + 1, 0)], nodos[(i + 1, 1)]],
                                     losa, A=A_losa, I=rigidez * I_COL))
    for lado in (0, 1):
        for dof in ("ux", "uy", "rz"):
            nodos[(0, lado)].fix_dof(dof, 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, nodos


def cerrada():
    """Frecuencias, modos, masas efectivas y respuesta espectral del edificio
    de cortante, con los modos analíticos (combinación hecha a mano)."""
    j = np.arange(1, N + 1)
    w = 2 * np.sqrt(K_PISO / M_PISO) * np.sin((2 * j - 1) * np.pi / (2 * (2 * N + 1)))
    i = np.arange(1, N + 1)
    phi = np.sin(np.outer(i, 2 * j - 1) * np.pi / (2 * N + 1))       # (piso, modo)
    L = M_PISO * phi.sum(axis=0)
    Mn = M_PISO * (phi**2).sum(axis=0)
    gamma = L / Mn
    masa_efectiva = L**2 / Mn
    sa = np.array([Sa(2 * np.pi / wj) for wj in w])
    u_modal = phi * gamma * sa / w**2                                # (piso, modo)
    cortante = masa_efectiva * sa
    return {"w": w, "T": 2 * np.pi / w, "phi": phi, "masa_efectiva": masa_efectiva,
            "Sa": sa, "u_modal": u_modal, "u_srss": np.sqrt((u_modal**2).sum(axis=1)),
            "V_srss": float(np.sqrt((cortante**2).sum())), "V_modal": cortante}


def espectral(combinacion: str = "SRSS"):
    """Análisis espectral con Solidum; devuelve el resultado y los nodos."""
    dom, nodos = marco()
    solver = ResponseSpectrumSolver(
        Assembler(dom), n_modes=N, direction={"dof_name": "ux"},
        spectrum=lambda w: Sa(2 * np.pi / w) / w**2,       # S_d = S_a/ω²
        combination=combinacion, damping=0.05)
    return solver.solve(), nodos


def calcular() -> dict:
    ref = cerrada()
    out = {"N": N, "H_m": H, "claro_m": CLARO, "E_GPa": E / 1e9, "lado_cm": 100 * LADO,
           "m_t": M_PISO / 1e3, "k_MN_m": K_PISO / 1e6, "rigidez": RIGIDEZ,
           "a0_g": A0 / G, "T_B": T_B, "T_C": T_C,
           "T_exacto": ref["T"], "masa_efectiva_pct_exacta": 100 * ref["masa_efectiva"] / (N * M_PISO),
           "Sa_g": ref["Sa"] / G, "V_srss_kN": ref["V_srss"] / 1e3,
           "V_srss_sobre_W": ref["V_srss"] / (N * M_PISO * G)}

    # Modal: frecuencias frente a la solución cerrada, con la rigidez del
    # ejemplo y con losas diez veces más rígidas.
    for clave, rigidez in (("error_w_max", RIGIDEZ), ("error_w_max_rigido", 10 * RIGIDEZ)):
        dom, nodos = marco(rigidez)
        modal = ModalSolver(Assembler(dom), n_modes=N).solve()
        orden = np.argsort(np.asarray(modal.frequencies_rad))
        w = np.asarray(modal.frequencies_rad)[orden]
        out[clave] = float(np.max(np.abs(w / ref["w"] - 1)))
        if rigidez == RIGIDEZ:
            out["T_fe"] = 2 * np.pi / w
            # Forma modal en el lado izquierdo, normalizada al último piso.
            forma = []
            for k in range(N):
                phi = np.asarray(modal.modes)[:, orden[k]]
                u = np.array([phi[nodos[(i, 0)].dofs["ux"]] for i in range(1, N + 1)])
                forma.append(u / u[-1])
            out["formas_fe"] = np.array(forma).T
            out["formas_exactas"] = ref["phi"] / ref["phi"][-1]
            out["error_formas"] = float(np.max(np.abs(out["formas_fe"] - out["formas_exactas"])))

    # Espectro con SRSS y con CQC.
    srss, nodos = espectral("SRSS")
    cqc, _ = espectral("CQC")
    fila = [nodos[(i, 0)].dofs["ux"] for i in range(1, N + 1)]
    u_srss = np.asarray(srss.u_combined)[fila]
    u_cqc = np.asarray(cqc.u_combined)[fila]
    masa_ef = np.asarray(srss.effective_masses)
    orden = np.argsort(np.asarray(srss.frequencies_rad))
    masa_ef = masa_ef[orden]
    out["u_srss_mm"] = 1e3 * u_srss
    out["u_srss_exacto_mm"] = 1e3 * ref["u_srss"]
    out["error_u_srss"] = float(np.max(np.abs(u_srss / ref["u_srss"] - 1)))
    out["error_masa_efectiva"] = float(np.max(np.abs(masa_ef - ref["masa_efectiva"])) / (N * M_PISO))
    out["masa_efectiva_acumulada_pct"] = 100 * np.cumsum(masa_ef) / (N * M_PISO)
    # Cortante basal: masa efectiva por aceleración espectral de cada modo, SRSS.
    w_s = np.asarray(srss.frequencies_rad)[orden]
    V = np.sqrt(sum((m * Sa(2 * np.pi / w))**2 for m, w in zip(masa_ef, w_s)))
    out["V_fe_kN"] = V / 1e3
    out["error_cortante"] = abs(V / ref["V_srss"] - 1)
    out["cqc_sobre_srss_max"] = float(np.max(u_cqc / u_srss))
    out["cqc_sobre_srss_min"] = float(np.min(u_cqc / u_srss))

    # Deriva de entrepiso: se combina la deriva de cada modo (bien), no la
    # diferencia de los desplazamientos ya combinados (mal).
    por_modo = np.asarray(srss.u_per_mode)[fila][:, orden]            # (piso, modo)
    derivas_modales = np.diff(np.vstack([np.zeros(N), por_modo]), axis=0)
    deriva_bien = np.sqrt((derivas_modales**2).sum(axis=1))
    deriva_mal = np.diff(np.concatenate([[0.0], u_srss]))
    out["deriva_bien_mm"] = 1e3 * deriva_bien
    out["deriva_mal_mm"] = 1e3 * deriva_mal
    out["deriva_mal_sobre_bien"] = float(deriva_mal[-1] / deriva_bien[-1])
    out["deriva_bien_pct_h"] = 100 * deriva_bien / H
    # Vistas escalares por modo y por piso, para las tablas del capítulo.
    out["modo"] = {str(k + 1): {
        "T_exacto": float(ref["T"][k]), "T_fe": float(out["T_fe"][k]),
        "masa_pct": float(100 * masa_ef[k] / (N * M_PISO)),
        "masa_acum_pct": float(out["masa_efectiva_acumulada_pct"][k]),
        "Sa_g": float(ref["Sa"][k] / G)} for k in range(N)}
    out["piso"] = {str(i + 1): {
        "u_mm": float(1e3 * u_srss[i]), "u_exacto_mm": float(1e3 * ref["u_srss"][i]),
        "deriva_bien_mm": float(1e3 * deriva_bien[i]),
        "deriva_mal_mm": float(1e3 * deriva_mal[i])} for i in range(N)}
    return out


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw={"width_ratios": [1.25, 1]})
    niveles = np.arange(0, N + 1)
    colores = [COLOR_FE, COLOR_SECUNDARIO, COLOR_TERCIARIO, "#7570b3", "#e7298a"]
    for k in range(3):
        fe = np.concatenate([[0.0], np.asarray(res["formas_fe"])[:, k]])
        ex = np.concatenate([[0.0], np.asarray(res["formas_exactas"])[:, k]])
        ax.plot(ex, niveles, color=colores[k], lw=3.0, alpha=0.3)
        ax.plot(fe, niveles, "o-", color=colores[k], ms=4, lw=1.0,
                label=f"modo {k + 1}, T = {res['T_fe'][k]:.3f} s")
    ax.axvline(0, color="0.5", lw=0.6)
    ax.set_xlabel("forma modal (1 en el último nivel)")
    ax.set_ylabel("nivel")
    ax.legend(fontsize=7.5, loc="lower left")
    T = np.linspace(0.01, 1.2, 300)
    bx.plot(T, [Sa(t) / G for t in T], color=COLOR_ANALITICO, lw=1.2)
    for k, t in enumerate(res["T_fe"]):
        bx.plot(t, Sa(t) / G, "o", color=colores[k], ms=5)
        bx.annotate(str(k + 1), (t, Sa(t) / G), textcoords="offset points",
                    xytext=(4, -11 if k > 1 else 4), fontsize=8, color=colores[k])
    bx.set_ylim(0, 0.85)
    bx.set_xlabel("periodo T [s]")
    bx.set_ylabel(r"$S_a / g$")
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_modos")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    print("T exactos  :", np.round(res["T_exacto"], 4))
    print("T Solidum  :", np.round(res["T_fe"], 4))
    print(f"error ω máx {res['error_w_max']:.1e} (losas ×10 más rígidas: {res['error_w_max_rigido']:.1e}); "
          f"error formas {res['error_formas']:.1e}")
    print("masa efectiva acumulada %:", np.round(res["masa_efectiva_acumulada_pct"], 2))
    print("u SRSS [mm]:", np.round(res["u_srss_mm"], 3), f" error {res['error_u_srss']:.1e}")
    print(f"cortante basal {res['V_fe_kN']:.1f} kN ({res['V_srss_sobre_W']:.3f} W), error {res['error_cortante']:.1e}")
    print(f"CQC/SRSS en [{res['cqc_sobre_srss_min']:.4f}, {res['cqc_sobre_srss_max']:.4f}]")
    print("deriva bien [mm]:", np.round(res["deriva_bien_mm"], 3))
    print("deriva mal  [mm]:", np.round(res["deriva_mal_mm"], 3), f" último piso: {res['deriva_mal_sobre_bien']:.3f}")
    return res


if __name__ == "__main__":
    main()
