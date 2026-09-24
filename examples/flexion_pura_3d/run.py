"""Ejemplo 6 del manual de ejemplos — flexión pura de un prisma 3D: Hex8 frente a Hex20.

Prisma de sección cuadrada en flexión pura, con la solución exacta de
Saint-Venant (curvatura anticlástica incluida). Con curvatura κ,
``σ_xx = E·κ·z`` y el resto de componentes nulas; el desplazamiento es

    u = κ·x·z,   v = −ν·κ·y·z,   w = −κ·x²/2 + (ν·κ/2)·(y² − z²)

y el momento, ``M = E·κ·I``. Se impone la curvatura con el desplazamiento
axial del extremo (``u = κ·L·z``) y se mide el momento de reacción.

El Hex20 contiene los polinomios cuadráticos completos y reproduce la
solución exacta con un solo elemento. El Hex8 se bloquea por cortante: dentro
de cada elemento su flecha es lineal en ``x`` y aparece una distorsión
angular parásita, ``σ_xz = G·κ·a/(2√3)`` en los puntos de Gauss, cuya energía
lo hace demasiado rígido.

Uso::

    python examples/flexion_pura_3d/run.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))

from ejemplos_comun import (  # noqa: E402
    COLOR_ANALITICO, COLOR_FE, COLOR_SECUNDARIO,
    estilo_figuras, guardar_figura, guardar_resultados,
)

from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.solid_3d import Hex8, Hex20  # noqa: E402
from solidum.materials.elastic_3d import Elastic3D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import LinearSolver  # noqa: E402

L, B = 1.0, 0.1          # longitud y lado de la sección cuadrada [m]
E, NU = 200.0e9, 0.3     # acero [Pa]
KAPPA = 0.01             # curvatura impuesta [1/m]: σ_max = E·κ·B/2 = 100 MPa
I = B**4 / 12
G = E / (2 * (1 + NU))
SIGMA_MAX = E * KAPPA * B / 2

# (elemento, elementos a lo largo, elementos por lado de la sección)
CASOS = [("Hex8", 5, 2), ("Hex8", 10, 2), ("Hex8", 20, 2), ("Hex8", 40, 2),
         ("Hex8", 40, 4), ("Hex8", 80, 4), ("Hex20", 1, 1), ("Hex20", 5, 1)]


def _clave(elemento, nx, n):
    return f"{elemento}_{nx}x{n}x{n}"


TOLERANCIAS = {
    # El Hex20 es exacto: momento, desplazamientos y esfuerzos al redondeo.
    **{f"{_clave(*c)}.error_M": 1e-10 for c in CASOS if c[0] == "Hex20"},
    **{f"{_clave(*c)}.error_u": 1e-10 for c in CASOS if c[0] == "Hex20"},
    **{f"{_clave(*c)}.sxz_sobre_smax": 1e-10 for c in CASOS if c[0] == "Hex20"},
    # El Hex8 es demasiado rígido: 2.6 veces en la malla más gruesa.
    "Hex8_5x2x2.M_sobre_exacto": (2.5, 2.7),
    "Hex8_80x4x4.M_sobre_exacto": (1.005, 1.03),
    # El cortante parásito en los puntos de Gauss es el de la fórmula.
    **{f"{_clave(*c)}.error_sxz_teorico": 1e-9 for c in CASOS if c[0] == "Hex8"},
    # El resto del exceso depende sólo de la sección (anticlástica): igual
    # para todas las longitudes de elemento con la misma sección.
    **{f"Hex8_{nx}x2x2.exceso_resto": (0.035, 0.042) for nx in (5, 10, 20, 40)},
    **{f"Hex8_{nx}x4x4.exceso_resto": (0.008, 0.011) for nx in (40, 80)},
    # Con ν = 0 el exceso es exactamente el del cortante parásito.
    "nu_cero.error_formula": 1e-9,
}


def exacta(x, y, z, nu=NU):
    """Desplazamiento exacto de la flexión pura (Saint-Venant)."""
    return (KAPPA * x * z, -nu * KAPPA * y * z,
            -KAPPA * x * x / 2 + nu * KAPPA / 2 * (y * y - z * z))


def malla(elemento: str, nx: int, n: int, nu: float = NU):
    """Prisma [0, L] × [−B/2, B/2]² con ``nx × n × n`` elementos ``Hex8`` o
    ``Hex20``. Los nodos se indexan en una rejilla (i, j, k); el Hex20 usa
    los vértices y los medios de arista de la rejilla de paso mitad."""
    dom = Domain()
    material = Elastic3D(E=E, nu=nu)
    q = 2 if elemento == "Hex20" else 1
    nodos = {}
    for k in range(q * n + 1):
        for j in range(q * n + 1):
            for i in range(q * nx + 1):
                if q == 2 and (i % 2) + (j % 2) + (k % 2) > 1:
                    continue             # centros de cara e interior: no son nodos del Hex20
                xyz = [L * i / (q * nx), -B / 2 + B * j / (q * n), -B / 2 + B * k / (q * n)]
                nodos[(i, j, k)] = dom.add_node(len(nodos) + 1, xyz)
    for k in range(n):
        for j in range(n):
            for i in range(nx):
                a, b, c = q * i, q * j, q * k
                vertices = [nodos[(a, b, c)], nodos[(a + q, b, c)], nodos[(a + q, b + q, c)],
                            nodos[(a, b + q, c)], nodos[(a, b, c + q)], nodos[(a + q, b, c + q)],
                            nodos[(a + q, b + q, c + q)], nodos[(a, b + q, c + q)]]
                eid = len(dom.elements) + 1
                if elemento == "Hex8":
                    dom.add_element(Hex8(eid, vertices, material))
                    continue
                medios = [nodos[(a + 1, b, c)], nodos[(a + 2, b + 1, c)], nodos[(a + 1, b + 2, c)],
                          nodos[(a, b + 1, c)], nodos[(a + 1, b, c + 2)], nodos[(a + 2, b + 1, c + 2)],
                          nodos[(a + 1, b + 2, c + 2)], nodos[(a, b + 1, c + 2)], nodos[(a, b, c + 1)],
                          nodos[(a + 2, b, c + 1)], nodos[(a + 2, b + 2, c + 1)], nodos[(a, b + 2, c + 1)]]
                dom.add_element(Hex20(eid, vertices + medios, material))
    return dom, nodos, q


def resolver(elemento: str, nx: int, n: int, nu: float = NU) -> dict:
    """Curvatura impuesta con ``u`` en los extremos; sólido rígido fijado en
    dos esquinas de la raíz con los valores exactos."""
    dom, nodos, q = malla(elemento, nx, n, nu)
    extremo = []
    for nodo in dom.nodes.values():
        x, y, z = nodo.coordinates
        if abs(x) < 1e-12:
            nodo.fix_dof("ux", 0.0)
        elif abs(x - L) < 1e-12:
            nodo.fix_dof("ux", KAPPA * L * z)
            extremo.append(nodo)
    esquina_a, esquina_b = nodos[(0, 0, 0)], nodos[(0, q * n, 0)]
    esquina_a.fix_dof("uy", exacta(*esquina_a.coordinates, nu)[1])
    esquina_a.fix_dof("uz", exacta(*esquina_a.coordinates, nu)[2])
    esquina_b.fix_dof("uz", exacta(*esquina_b.coordinates, nu)[2])

    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    U = LinearSolver(asm).solve(np.zeros(dom.total_dofs))
    _, F_int = asm.assemble_non_linear_system(U)
    # Momento de reacción en el extremo: Σ R_x·z (R = fuerza interna nodal).
    M = sum(F_int[nodo.dofs["ux"]] * nodo.coordinates[2] for nodo in extremo)

    error_u = max(abs(U[nodo.dofs[d]] - exacta(*nodo.coordinates, nu)[k])
                  for nodo in dom.nodes.values() for k, d in enumerate(("ux", "uy", "uz")))
    sxz = 0.0
    for el in dom.elements.values():
        estado = el.compute_gauss_state(U)
        sxz = max(sxz, float(np.max(np.abs(np.asarray(estado["stress"])[:, 5]))))
    a = L / nx                                    # longitud del elemento
    g = E / (2 * (1 + nu))
    sxz_teorico = g * KAPPA * a / (2 * np.sqrt(3)) if elemento == "Hex8" else 0.0
    # Energía del cortante parásito frente a la de flexión: (a/h)²/(2(1+ν)).
    cortante = (a / B) ** 2 / (2 * (1 + nu)) if elemento == "Hex8" else 0.0
    return {"elemento": elemento, "nx": nx, "n": n, "gdl": int(dom.total_dofs),
            "a_sobre_h": a / B, "M_sobre_exacto": M / (E * KAPPA * I),
            "exceso_cortante": cortante,
            "exceso_resto": M / (E * KAPPA * I) - 1 - cortante,
            "error_M": abs(M / (E * KAPPA * I) - 1),
            "error_M_pct": 100 * abs(M / (E * KAPPA * I) - 1),
            "error_u": error_u / (KAPPA * L * L / 2),
            "sxz_sobre_smax": sxz / SIGMA_MAX,
            "sxz_teorico_sobre_smax": sxz_teorico / SIGMA_MAX,
            "error_sxz_teorico": abs(sxz - sxz_teorico) / SIGMA_MAX}


def calcular() -> dict:
    out = {"L_m": L, "B_cm": 100 * B, "E_GPa": E / 1e9, "nu": NU, "kappa": KAPPA,
           "sigma_max_MPa": SIGMA_MAX / 1e6, "M_kNm": E * KAPPA * I / 1e3,
           "flecha_mm": 1e3 * KAPPA * L * L / 2,
           "factor_sxz": 1 / (2 * np.sqrt(3) * (1 + NU))}
    for caso in CASOS:
        out[_clave(*caso)] = resolver(*caso)
    # Con ν = 0 no hay curvatura anticlástica: el exceso del Hex8 es sólo el
    # del cortante parásito, 1 + (a/h)²/2 exacto.
    d = resolver("Hex8", 5, 2, nu=0.0)
    out["nu_cero"] = {"M_sobre_exacto": d["M_sobre_exacto"],
                      "error_formula": abs(d["M_sobre_exacto"] - 1 - d["exceso_cortante"])}
    return out


def figuras(res: dict, carpeta: Path) -> None:
    estilo_figuras()
    import matplotlib.pyplot as plt

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for n, marca, color in ((2, "o", COLOR_FE), (4, "s", COLOR_SECUNDARIO)):
        d = [res[_clave(e, nx, m)] for e, nx, m in CASOS if e == "Hex8" and m == n]
        ax.plot([x["gdl"] for x in d], [x["M_sobre_exacto"] for x in d], marca + "-",
                color=color, ms=4, label=f"Hex8, {n}×{n} en la sección")
        bx.plot([x["a_sobre_h"] for x in d], [x["sxz_sobre_smax"] for x in d], marca,
                color=color, ms=5, label=f"Hex8, {n}×{n}")
    d = [res[_clave(e, nx, m)] for e, nx, m in CASOS if e == "Hex20"]
    ax.plot([x["gdl"] for x in d], [x["M_sobre_exacto"] for x in d], "^", color=COLOR_ANALITICO,
            ms=6, label="Hex20, 1×1 en la sección")
    ax.axhline(1.0, color="0.5", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("grados de libertad")
    ax.set_ylabel(r"$M / M_{\mathrm{exacto}}$ (misma curvatura)")
    ax.legend(fontsize=7.5)
    r = np.linspace(0, 2.1, 50)
    bx.plot(r, res["factor_sxz"] * r, color=COLOR_ANALITICO, lw=1.0,
            label=r"$G\kappa a/(2\sqrt{3})$")
    bx.set_xlabel(r"longitud del elemento / canto, $a/h$")
    bx.set_ylabel(r"$\sigma_{xz}$ parásito $/\ \sigma_{\max}$")
    bx.legend(fontsize=7.5)
    fig.tight_layout()
    guardar_figura(fig, carpeta, "fig_bloqueo")


def main() -> dict:
    res = calcular()
    guardar_resultados(AQUI, res)
    figuras(res, AQUI)
    print(f"M exacto = {res['M_kNm']:.3f} kN·m   σ_max = {res['sigma_max_MPa']:.0f} MPa")
    for caso in CASOS:
        d = res[_clave(*caso)]
        print(f"{_clave(*caso):14s} gdl={d['gdl']:5d}  M/M_ex={d['M_sobre_exacto']:.4f}  "
              f"error u={d['error_u']:.1e}  σxz/σmax={d['sxz_sobre_smax']:.4f} "
              f"(teórico {d['sxz_teorico_sobre_smax']:.4f})")
    return res


if __name__ == "__main__":
    main()
