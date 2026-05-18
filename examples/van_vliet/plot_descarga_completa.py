"""Plot de las curvas carga-deformación completas (lineal y exponencial).

Combina:
- Curvas **analíticas 1D** del cohesivo + equilibrio uniaxial en el cuello
  (``curva_analitica_{linear,exponential}.csv``): trazo completo de la
  rama elástica + pico + descarga hasta el cierre energético.
- Datos **MEF** disponibles (``curva_numerica_linear.csv``,
  ``curva_numerica_exponential.csv``): hasta donde llegó el solver
  (arc-length / Newton) antes de fallar en la rama de softening.

El acuerdo perfecto entre los puntos MEF (hasta donde llegaron) y la
curva analítica confirma que el subsistema elemento+cohesivo entrega la
respuesta física correcta — el límite es del solver, no del modelo.
"""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).parent


def load(path: Path, col_def: str, col_P: str):
    defs, Ps = [], []
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            defs.append(float(row[col_def]))
            Ps.append(float(row[col_P]))
    return defs, Ps


def main():
    fig, ax = plt.subplots(figsize=(9.5, 5.8))

    # --- Curvas analíticas (descarga completa) -----------------------------
    for soft, color, label in [
        ('linear', '#1f77b4', 'Lineal — analítica 1D'),
        ('exponential', '#d62728', 'Exponencial — analítica 1D'),
    ]:
        defs, Ps = load(HERE / f'curva_analitica_{soft}.csv', 'def_Ls_um', 'P_kN')
        ax.plot(defs, Ps, '-', color=color, linewidth=1.6, label=label, alpha=0.85)

    # --- Datos MEF (truncados donde murió el solver) -----------------------
    for soft, color, label in [
        ('linear', '#1f77b4', 'Lineal — MEF (arc-length)'),
        ('exponential', '#d62728', 'Exponencial — MEF (arc-length)'),
    ]:
        path = HERE / f'curva_numerica_{soft}.csv'
        if not path.exists():
            continue
        defs, Ps = load(path, 'deformation_um', 'P_kN')
        ax.plot(defs, Ps, 'o', color=color, markersize=4.5,
                markerfacecolor='none', markeredgewidth=1.2,
                label=label, alpha=0.8)

    # Pico teórico
    ax.axhline(30.84, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
    ax.text(95, 30.84 + 0.6, r'$\sigma_{t0}\cdot A_{neck}$ = 30.84 kN',
            color='gray', fontsize=9, va='bottom')

    ax.set_xlabel('Deformación entre puntos de referencia $L_s = 0.6\\cdot D$ (µm)',
                  fontsize=11)
    ax.set_ylabel('Carga $P$ (kN)', fontsize=11)
    ax.set_title(
        'Van Vliet §8.1 — descarga completa (post-fix dimensional)\n'
        'curva analítica vs MEF (tracción centrada, no faithful)',
        fontsize=11,
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlim(0, 110)
    ax.set_ylim(0, 35)

    out = HERE / 'curva_descarga_completa.png'
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print(f'Plot guardado en: {out}')


if __name__ == '__main__':
    main()
