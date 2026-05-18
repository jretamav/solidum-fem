"""Plot diagnóstico de la curva carga-deformación del ensayo Van Vliet.

Lee ``curva_numerica_linear.csv`` (resultado del ensayo de tracción
centrada con softening lineal, post-fix del bug del thickness) y la
dibuja con anotaciones de pico y rama de softening.

Este ensayo **no es Van Vliet faithful** — la carga está centrada
(``x = 0`` en lugar de ``x = +e = 4 mm``) porque la malla con ``n_x = 21``
no tiene nodo en la línea de excentricidad. Sirve como diagnóstico del
subsistema cohesivo + condensación + arc-length.
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # backend no interactivo (headless)
import matplotlib.pyplot as plt

HERE = Path(__file__).parent


def main():
    csv_path = HERE / 'curva_numerica_linear.csv'
    steps, lambdas, P_kN, def_um = [], [], [], []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            steps.append(int(row['step']))
            lambdas.append(float(row['load_factor']))
            P_kN.append(float(row['P_kN']))
            def_um.append(float(row['deformation_um']))

    # Identificar el pico
    i_peak = max(range(len(P_kN)), key=lambda i: P_kN[i])
    P_peak = P_kN[i_peak]
    def_peak = def_um[i_peak]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.plot(def_um, P_kN, '-o', color='#1f77b4', linewidth=1.3,
            markersize=3.5, label='Numérico (softening lineal)')

    # Pico
    ax.axvline(def_peak, color='gray', linestyle=':', linewidth=0.7, alpha=0.6)
    ax.axhline(P_peak, color='gray', linestyle=':', linewidth=0.7, alpha=0.6)
    ax.plot([def_peak], [P_peak], 'o', color='#d62728', markersize=8,
            zorder=5, label=f'Pico: P={P_peak:.2f} kN, def={def_peak:.2f} µm')

    # Referencia teórica: sigma_t0 * A_neck = 30.84 kN
    P_theo = 2.57e6 * (0.120 * 0.100) / 1000  # kN
    ax.axhline(P_theo, color='#2ca02c', linestyle='--', linewidth=1.2,
               alpha=0.7, label=f'Teórico σ_t0·A_neck = {P_theo:.2f} kN')

    # Sombrear rama de softening
    if i_peak < len(P_kN) - 1:
        ax.axvspan(def_peak, max(def_um), color='red', alpha=0.06,
                   label='Rama de softening')

    ax.set_xlabel('Deformación entre puntos de referencia L_s = 0.6·D (µm)',
                  fontsize=11)
    ax.set_ylabel('Carga P (kN)', fontsize=11)
    ax.set_title(
        'Van Vliet §8.1 — diagnóstico de subsistema embedded discontinuity\n'
        '(tracción centrada, no faithful; post-fix dimensional thickness)',
        fontsize=11,
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower right', fontsize=9)
    ax.set_xlim(0, max(def_um) * 1.05)
    ax.set_ylim(0, P_peak * 1.15)

    out = HERE / 'curva_diagnostico_linear.png'
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print(f'Plot guardado en: {out}')


if __name__ == '__main__':
    main()
