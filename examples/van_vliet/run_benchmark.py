"""Benchmark Van Vliet (§8.1, Retama 2010) — fase 4 del ADR 0010.

Ejecuta la simulación del dog-bone Van Vliet con ``CST_Embedded2D``
distribuido + crack prescrita pre-activada + ``ArcLengthSolver``.

Uso::

    python -m examples.van_vliet.run_benchmark --softening linear
    python -m examples.van_vliet.run_benchmark --softening exponential

Salida: CSV con la curva carga-deformación numérica para comparativa
contra ``curvas_referencia.csv`` (digitalización de fig 8.3 de la tesis).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np

# Permitir ejecución directa desde la raíz del repo
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.solvers import ArcLengthSolver

from examples.van_vliet._mesh import build_van_vliet_domain


HERE = Path(__file__).parent

# Carga de referencia: 35 kN, ligeramente por encima del pico esperado
# (~27 kN). El arc-length escala esta referencia con su factor λ. El
# pico físico debe alcanzarse en λ ≈ 27/35 ≈ 0.77.
P_REF = 35.0e3  # N


def run(softening: str, out_csv: str | None = None, verbose: bool = True):
    """Construye el dominio, corre arc-length, devuelve la curva.

    Returns
    -------
    list of (step, lambda, P_N, deformation_m)
    """
    if verbose:
        print(f'[Van Vliet §8.1] Construyendo dominio con softening={softening}...')
    domain, info = build_van_vliet_domain(softening=softening)

    if verbose:
        print(f'  Nodos: {len(domain.nodes)}, elementos: {len(domain.elements)}, '
              f'DOFs: {domain.total_dofs}')
        print(f'  Crack elements (pre-activados): {len(info["crack_elements"])}')
        print(f'  load_node: id={info["load_node"].id} '
              f'coords={info["load_node"].coordinates}')
        print(f'  fix_node:  id={info["fix_node"].id}')
        print(f'  ref_top:   coords={info["ref_top_node"].coordinates}')
        print(f'  ref_bot:   coords={info["ref_bot_node"].coordinates}')

    assembler = Assembler(domain)

    # Vector de fuerza externa de referencia: P_REF en uy del load_node.
    F_ext_ref = np.zeros(domain.total_dofs)
    F_ext_ref[info['load_node'].dofs['uy']] = P_REF

    # Captura de la curva por callback.
    records = []
    ref_top_uy_dof = info['ref_top_node'].dofs['uy']
    ref_bot_uy_dof = info['ref_bot_node'].dofs['uy']

    def step_callback(step, U_current, load_factor):
        P_total = load_factor * P_REF
        defo = U_current[ref_top_uy_dof] - U_current[ref_bot_uy_dof]
        records.append((step, load_factor, P_total, defo))
        if verbose:
            print(f'    Paso {step:3d}: lambda={load_factor:.4f}  '
                  f'P={P_total*1e-3:6.2f} kN  def={defo*1e6:7.2f} um')

    # ArcLengthSolver: control por longitud de arco. Adecuado para curvas
    # con softening (snap-back potencial).
    # rtol_disp aflojado a 1e-2: con K_e cohesivo de 1e15 N/m³ la
    # calibración disp_scale = force/K_diag colapsa a valores muy chicos;
    # el criterio dual saturaría la tolerancia de desplazamiento. La
    # fuerza converge con margen amplio (1e-4) y eso es el indicador
    # real de equilibrio en este problema.
    conv = ConvergenceCriterion(rtol_force=1e-4, rtol_disp=1e-2)
    solver = ArcLengthSolver(
        assembler,
        convergence=conv,
        max_iter=30,
        max_lambda=1.5,
        initial_dl=5.0e-6,   # ~ 5 µm de norma de incremento de desplazamiento
        max_steps=200,
        dl_grow_factor=1.3,
        dl_max_factor=4.0,
    )

    if verbose:
        print(f'[Van Vliet §8.1] Resolviendo...')

    try:
        solver.solve(F_ext_ref, step_callback=step_callback)
        status = 'OK'
    except RuntimeError as exc:
        status = f'STOPPED: {exc}'

    if verbose:
        print(f'[Van Vliet §8.1] Solver terminó: {status}')
        print(f'  Pasos convergidos: {len(records)}')

    # Escribir CSV
    if out_csv is None:
        out_csv = HERE / f'curva_numerica_{softening}.csv'

    with open(out_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['step', 'load_factor', 'P_kN', 'deformation_um'])
        for step, lam, P, defo in records:
            writer.writerow([step, f'{lam:.6f}', f'{P*1e-3:.4f}', f'{defo*1e6:.4f}'])
    if verbose:
        print(f'  Curva escrita en: {out_csv}')

    return records


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--softening', choices=['linear', 'exponential'], default='linear')
    p.add_argument('--out', default=None, help='CSV de salida (opcional)')
    p.add_argument('--quiet', action='store_true')
    args = p.parse_args()
    run(softening=args.softening, out_csv=args.out, verbose=not args.quiet)
