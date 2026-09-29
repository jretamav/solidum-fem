"""Ensayo de tracción de Van Vliet con grieta prescrita (Retama 2010, §8.1).

Modelo compartido por el test de validación (``test_disc_van_vliet.py``) y por
la figura de ``docs/user/discontinuities/validacion/``. Datos de la tesis
(PDF pp. 65-68):

- Probeta en hueso de perro «type C»: D = 200 mm, entallas circulares de
  radio 145 mm y 40 mm de profundidad (cuello de 120 mm), tramos rectos de
  0.25 D, espesor 100 mm, esfuerzo plano.
- Concreto: E = 39.8 GPa, ν = 0.2, σ_t0 = 2.57 MPa, G_F = 0.1219 N/mm.
- Dos placas de acero de dos capas de elementos; carga puntual con
  excentricidad e = D/50 y giro libre de las placas.
- Grieta prescrita: una banda horizontal de elementos embebidos a media
  altura; el resto, CST estándar.
- Deformación δ: desplazamiento vertical relativo entre dos puntos del eje
  separados L_s = 0.6 D.

Supuestos propios (la tesis no los da): E del acero = 210 GPa, ν = 0.3; las
rótulas no se desplazan lateralmente (la reacción horizontal resulta nula);
la carga se aplica como desplazamiento impuesto en la rótula superior.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

import solidum  # noqa: F401
from solidum.core.domain import Domain
from solidum.elements.solid_2d import Tri3
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver
from solidum.user.discontinuities import CST_Embedded2D, CohesiveDamageIsotropic

D, RADIO, PROFUNDIDAD, ESPESOR = 0.200, 0.145, 0.040, 0.100      # [m]
E_C, NU_C = 39.8e9, 0.2
SIGMA_T0, G_F = 2.57e6, 121.9                                     # [Pa], [N/m]
E_ACERO, NU_ACERO = 210.0e9, 0.3
EXCENTRICIDAD = D / 50
L_S = 0.6 * D
H_BANDA, H_CAPA = 0.008, 0.008
K_E = 1.0e14                        # penalización [Pa/m]: 1e13 da lo mismo
_XC = D / 2 - PROFUNDIDAD + RADIO   # centro del arco de la entalla derecha

# Valores leídos de la fig. 8.3 de la tesis (±0.3 kN en el pico, ±0.5 kN en
# las ramas, ±1 µm en δ).
TESIS = {
    "exponential": {"P_max": 27.0, "d_pico": 7.0, 20: 19.1, 50: 8.3, 100: 3.1},
    "linear": {"P_max": 27.4, "d_pico": 7.5, 20: 22.4, 50: 8.5, 100: 2.4},
}


def semiancho(y: float) -> float:
    return _XC - np.sqrt(RADIO**2 - y**2) if abs(y) <= D / 2 else D / 2


def modelo(softening: str, n_lado: int = 7):
    """Malla estructurada mapeada al contorno de las entallas. Devuelve el
    dominio, la rótula superior, la inferior y los dos puntos de medida."""
    pos = [H_BANDA / 2 + 0.008 * k for k in range(8)] + [0.070 + 0.010 * k for k in range(9)]
    ys = sorted([-y for y in pos] + pos)          # sin nodo en y = 0: la banda la cruza
    paso = (1 - 0.04) / n_lado
    xi_pos = [0.04 + paso * k for k in range(n_lado + 1)]      # ξ = ±0.04: rótulas en x = −e
    xis = sorted([-x for x in xi_pos] + [0.0] + xi_pos)

    dom = Domain()
    concreto = Elastic2D(E=E_C, nu=NU_C, hypothesis="plane_stress")
    acero = Elastic2D(E=E_ACERO, nu=NU_ACERO, hypothesis="plane_stress")
    cohesivo = CohesiveDamageIsotropic(sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_E, softening=softening)
    nodos: dict = {}

    def nodo(x, y):
        clave = (round(x, 9), round(y, 9))
        if clave not in nodos:
            nodos[clave] = dom.add_node(len(nodos) + 1, [x, y])
        return nodos[clave]

    def triangulo(a, b, c, tipo):
        eid = len(dom.elements) + 1
        if tipo == "banda":
            dom.add_element(CST_Embedded2D(eid, [a, b, c], concreto, cohesivo, thickness=ESPESOR))
        else:
            dom.add_element(Tri3(eid, [a, b, c], acero if tipo == "placa" else concreto,
                                 thickness=ESPESOR))

    def fila(xs0, y0, xs1, y1, tipo):
        for i in range(len(xs0) - 1):
            a, b, c, d = nodo(xs0[i], y0), nodo(xs0[i + 1], y0), nodo(xs1[i + 1], y1), nodo(xs1[i], y1)
            # Diagonales simétricas respecto del eje: rectángulos casi exactos
            # en el cuello, así que la banda tiene triángulos rectángulos.
            if xs0[i] + xs0[i + 1] < 0:
                triangulo(a, b, c, tipo)
                triangulo(a, c, d, tipo)
            else:
                triangulo(a, b, d, tipo)
                triangulo(b, c, d, tipo)

    for j in range(len(ys) - 1):
        y0, y1 = ys[j], ys[j + 1]
        fila([x * semiancho(y0) for x in xis], y0, [x * semiancho(y1) for x in xis], y1,
             "banda" if y0 < 0 < y1 else "probeta")

    xs_placa = sorted([-0.6 * D] + [x * D / 2 for x in xis] + [0.6 * D])
    borde = 0.75 * D
    for capa in range(2):
        for signo in (1, -1):
            y0, y1 = borde + capa * H_CAPA, borde + (capa + 1) * H_CAPA
            y0, y1 = (y0, y1) if signo > 0 else (-y1, -y0)
            fila(xs_placa, y0, xs_placa, y1, "placa")

    tope = borde + 2 * H_CAPA
    return (dom, nodo(-EXCENTRICIDAD, tope), nodo(-EXCENTRICIDAD, -tope),
            (nodo(0.0, L_S / 2), nodo(0.0, -L_S / 2)))


def ensayo(softening: str, u_max: float = 0.12e-3, paso: float = 0.5e-6):
    """Curva P-δ con desplazamiento impuesto en la rótula superior. El paso
    de 0.5 µm está convergido (0.25 y 0.125 µm cambian el pico < 0.1 kN);
    con 1 µm el pico sube ~0.9 kN por el sobrepaso de la activación, que se
    evalúa con el estado del paso anterior."""
    dom, arriba, abajo, (m1, m2) = modelo(softening)
    abajo.fix_dof("ux", 0.0)
    abajo.fix_dof("uy", 0.0)
    arriba.fix_dof("ux", 0.0)
    arriba.fix_dof("uy", u_max)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    curva = [(0.0, 0.0)]

    def registrar(step, U, lam):
        _, F_int = asm.assemble_non_linear_system(U)
        curva.append((U[m1.dofs["uy"]] - U[m2.dofs["uy"]], F_int[arriba.dofs["uy"]]))

    nivel = logging.getLogger("solidum").level
    logging.getLogger("solidum").setLevel(logging.ERROR)
    try:
        NonlinearSolver(asm, num_steps=int(round(u_max / paso))).solve(
            np.zeros(dom.total_dofs), step_callback=registrar)
    finally:
        logging.getLogger("solidum").setLevel(nivel)
    return np.array(curva)


def resumen(curva: np.ndarray) -> dict:
    d_um, P_kN = 1e6 * curva[:, 0], 1e-3 * curva[:, 1]
    i = int(np.argmax(P_kN))
    out = {"P_max": float(P_kN[i]), "d_pico": float(d_um[i])}
    for d in (20, 50, 100):
        out[d] = float(np.interp(d, d_um, P_kN))
    return out


def figura(destino: Path) -> None:
    """Curvas completas hasta 200 µm frente a los valores leídos de la tesis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.2))
    for ley, color, marca, nombre in (("exponential", "#1f5fa8", "o", "exponencial"),
                                      ("linear", "#d95f02", "s", "lineal")):
        c = ensayo(ley, u_max=0.25e-3)
        ax.plot(1e6 * c[:, 0], 1e-3 * c[:, 1], color=color, lw=1.6, label=f"Solidum, {nombre}")
        t = TESIS[ley]
        ax.plot([t["d_pico"], 20, 50, 100], [t["P_max"], t[20], t[50], t[100]], marca,
                mfc="white", mec=color, label=f"Tesis, {nombre} (leído de la fig. 8.3)")
    ax.plot([6.6, 20, 50, 100, 200], [27.5, 14.7, 8.3, 3.7, 1.3], "k^", ms=5,
            label="Experimental Van Vliet (leído)")
    ax.set_xlabel("deformación δ en L_s = 120 mm [µm]")
    ax.set_ylabel("P [kN]")
    ax.set_xlim(0, 210)
    ax.set_ylim(0, 30)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(destino, dpi=130)


if __name__ == "__main__":
    raiz = Path(__file__).resolve().parents[3]
    figura(raiz / "docs" / "user" / "discontinuities" / "validacion" / "van_vliet_fase1.png")
