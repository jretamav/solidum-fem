"""Van Vliet diagnóstico — control por desplazamiento del nodo de carga.

A diferencia de ``run_benchmark.py`` (control por fuerza vía
``ArcLengthSolver``), aquí se prescribe el desplazamiento ``uy`` del nodo
de carga y se lee la reacción como ``P``. Esto sortea el snap-back del
arc-length cilíndrico en la rama de softening y permite trazar la curva
carga-deformación completa hasta el residuo de saturación.

Mismo modelo (geometría, materiales, crack pre-activado, BCs salvo la
carga) que ``run_benchmark.py``.
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from fenix.math.assembly import Assembler
from fenix.math.convergence import ConvergenceCriterion
from fenix.math.solvers import NonlinearSolver

from examples.van_vliet._mesh import build_van_vliet_domain


HERE = Path(__file__).parent

# Desplazamiento objetivo en el nodo de carga (uy). Sweep de 0 a U_TARGET
# en ``num_steps`` pasos. Para Van Vliet con linear softening y
# `w_c ≈ 95 µm`, ~50 µm de uy en el load_node basta para trazar buena
# parte de la rama post-pico.
U_TARGET = 50.0e-6  # m

NUM_STEPS = 100


def run(softening: str, out_csv: str | None = None, verbose: bool = True):
    if verbose:
        print(f'[Van Vliet disp-ctrl] Construyendo dominio con softening={softening}...')
    domain, info = build_van_vliet_domain(softening=softening)

    # Cambiar BC del load_node: prescribir uy = U_TARGET en lugar de
    # dejarla libre. Quedaba "libre + carga aplicada"; ahora queda
    # "prescrita + sin carga directa". El factor de carga del solver
    # escala el valor prescrito de 0 a U_TARGET.
    load_node = info['load_node']
    load_node.fix_dof('uy', U_TARGET)

    if verbose:
        print(f'  Nodos: {len(domain.nodes)}, elementos: {len(domain.elements)}')
        print(f'  Prescrito uy({load_node.id}) = {U_TARGET*1e6:.2f} µm')

    assembler = Assembler(domain)

    # F_ext nula: la carga "entra" por la prescripción del desplazamiento.
    F_ext = np.zeros(domain.total_dofs)

    load_uy_dof = load_node.dofs['uy']
    ref_top_uy_dof = info['ref_top_node'].dofs['uy']
    ref_bot_uy_dof = info['ref_bot_node'].dofs['uy']

    records = []

    def step_callback(step, U_current, load_factor):
        # Reacción = F_int en el DOF prescrito tras converger. La obtenemos
        # re-ensamblando el sistema interno (no implica iteración extra
        # significativa porque K se reusa internamente en el solver).
        _, F_int = assembler.assemble_non_linear_system(U_current)
        P_reaction = F_int[load_uy_dof]  # N
        uy_load = U_current[load_uy_dof]
        defo = U_current[ref_top_uy_dof] - U_current[ref_bot_uy_dof]
        records.append((step, load_factor, P_reaction, uy_load, defo))
        if verbose:
            print(f'    Paso {step:3d}: lam={load_factor:.4f} '
                  f'uy_load={uy_load*1e6:6.2f} um  '
                  f'P={P_reaction*1e-3:7.3f} kN  def_Ls={defo*1e6:7.2f} um')

    conv = ConvergenceCriterion(rtol_force=1e-4, rtol_disp=1e-2)
    solver = NonlinearSolver(
        assembler,
        convergence=conv,
        max_iter=30,
        num_steps=NUM_STEPS,
        adaptive=True,
        min_delta_lambda=1e-5,
    )

    if verbose:
        print(f'[Van Vliet disp-ctrl] Resolviendo...')

    try:
        solver.solve(F_ext, step_callback=step_callback)
        status = 'OK'
    except RuntimeError as exc:
        status = f'STOPPED: {exc}'

    if verbose:
        print(f'[Van Vliet disp-ctrl] Solver terminó: {status}')
        print(f'  Pasos convergidos: {len(records)}')

    if out_csv is None:
        out_csv = HERE / f'curva_dispctrl_{softening}.csv'

    with open(out_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['step', 'load_factor', 'P_kN', 'uy_load_um', 'def_Ls_um'])
        for step, lam, P, uy, defo in records:
            writer.writerow([
                step, f'{lam:.6f}', f'{P*1e-3:.4f}',
                f'{uy*1e6:.4f}', f'{defo*1e6:.4f}',
            ])
    if verbose:
        print(f'  Curva escrita en: {out_csv}')

    return records


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--softening', choices=['linear', 'exponential'], default='linear')
    p.add_argument('--out', default=None)
    p.add_argument('--quiet', action='store_true')
    args = p.parse_args()
    run(softening=args.softening, out_csv=args.out, verbose=not args.quiet)
