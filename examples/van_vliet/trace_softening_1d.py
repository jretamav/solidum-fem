"""Curva carga-deformación analítica 1D del problema Van Vliet diagnóstico.

Recorre directamente el material cohesivo barriendo ``jump_n`` desde 0
hasta el cierre energético, y para cada estado compone el equilibrio 1D
(tracción uniforme en el cuello + apertura puntual de la grieta) para
obtener ``(P, def_Ls)``. **No depende de ningún solver** — sortea las
limitaciones del arc-length / Newton en la rama de softening y muestra
la respuesta completa que el subsistema elemento+cohesivo entrega *si*
se le condujera por desplazamiento de la grieta.

Para cada ``jump_n``:
- ``t = (1-ω)·K_e·jump_n`` (cohesivo, salida directa).
- ``σ_neck = t`` (continuidad de tracción normal en Γ_d).
- ``P = σ_neck · A_neck``.
- Deformación elástica del gauge ``L_s = 0.6·D``: integral de ``σ(z)/E``
  sobre ``z ∈ [-L_s/2, L_s/2]`` (todo dentro del cuello curvo).
- ``def_Ls = elastic_def + jump_n`` (apertura de Γ_d se suma).
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from examples.van_vliet._mesh import (
    D, R, DEPTH, E_CONCRETE, SIGMA_T0, G_F, K_E_COHESIVE,
    BODY_MIN_WIDTH, X_C_ARC, L_S, half_width,
)
from fenix.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic


HERE = Path(__file__).parent

A_NECK = BODY_MIN_WIDTH * DEPTH   # 0.120 m × 0.100 m = 0.012 m²


def gauge_elastic_def(P: float, n_quad: int = 200) -> float:
    """Deformación elástica de L_s = 0.6·D bajo la carga ``P`` (axial uniforme).

    Integra ``ε(z) = σ(z)/E = P/(2·half_width(z)·depth·E)`` sobre el
    gauge ``z ∈ [-L_s/2, +L_s/2]`` (todo dentro del cuello curvo).
    """
    z_a, z_b = -L_S / 2.0, +L_S / 2.0
    zs = np.linspace(z_a, z_b, n_quad)
    eps = np.array([P / (2.0 * half_width(z) * DEPTH * E_CONCRETE) for z in zs])
    return float(np.trapezoid(eps, zs))


def trace(softening: str, n_points: int = 400):
    coh = CohesiveDamageIsotropic(
        sigma_t0=SIGMA_T0, G_f=G_F, K_e=K_E_COHESIVE, softening=softening,
    )

    # Sweep de jump_n: desde 0 hasta una apertura "muy abierta".
    # En lineal w_c es exacto; en exponencial usamos ~10·G_F/σ_t0 para llegar
    # a la cola asintótica.
    if softening == 'linear':
        jump_max = coh.w_c * 1.02
    else:
        jump_max = 10.0 * G_F / SIGMA_T0

    # Densificamos al inicio (rama elástica + zona del pico) y aflojamos en la cola.
    elastic_part = np.linspace(0.0, coh.kappa_0 * 1.5, n_points // 4)
    post_peak = np.linspace(coh.kappa_0 * 1.5, jump_max, 3 * n_points // 4)
    jumps = np.unique(np.concatenate([elastic_part, post_peak]))

    records = []
    for j in jumps:
        jump_vec = np.array([j, 0.0])
        t_local, _, state = coh.compute_traction(jump_vec, state_vars=None)
        sigma_neck = float(t_local[0])
        P = sigma_neck * A_NECK
        if P < 0.0:
            break
        elastic_def = gauge_elastic_def(P) if P > 0.0 else 0.0
        total_def = elastic_def + j  # apertura de la grieta se suma directamente
        records.append((j, sigma_neck, P, elastic_def, total_def,
                        state['damage'], state['kappa']))
    return records


def write_csv(records, path: Path):
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['jump_n_m', 'sigma_neck_MPa', 'P_kN',
                         'elastic_def_um', 'def_Ls_um', 'damage', 'kappa'])
        for j, s, P, el, tot, d, k in records:
            writer.writerow([
                f'{j:.4e}', f'{s*1e-6:.4f}', f'{P*1e-3:.4f}',
                f'{el*1e6:.4f}', f'{tot*1e6:.4f}', f'{d:.4f}', f'{k:.4e}',
            ])


if __name__ == '__main__':
    print(f'A_neck = {A_NECK*1e6:.1f} mm² = {A_NECK:.4f} m²')
    print(f'σ_t0 · A_neck = {SIGMA_T0 * A_NECK / 1000:.2f} kN (pico teórico)')
    print()
    for soft in ('linear', 'exponential'):
        recs = trace(soft)
        out = HERE / f'curva_analitica_{soft}.csv'
        write_csv(recs, out)
        i_peak = max(range(len(recs)), key=lambda i: recs[i][2])
        print(f'[{soft:11s}] {len(recs):4d} puntos | '
              f'pico P={recs[i_peak][2]*1e-3:.2f} kN @ '
              f'def={recs[i_peak][4]*1e6:.2f} µm | '
              f'cola P={recs[-1][2]*1e-3:.3f} kN @ def={recs[-1][4]*1e6:.1f} µm')
        print(f'             → {out}')
