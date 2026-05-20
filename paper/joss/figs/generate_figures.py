"""Genera la figura mosaico de validación del paper JOSS.

Reutiliza los helpers de ``tests/validation`` para no duplicar la física;
cada subplot corre el mismo benchmark que los tests ejecutan en CI y
añade el plot del campo o la curva de convergencia correspondiente.

Salida: ``paper/joss/figs/validation_mosaic.png`` (DPI 150).

Uso:
    python paper/joss/figs/generate_figures.py
"""
import math
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Permite importar tanto de ``solidum`` como de ``tests.validation``.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.solid_2d import Quad4, Quad8  # noqa: E402
from solidum.materials.elastic_2d import Elastic2D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.math.solvers import LinearSolver  # noqa: E402

from tests.validation.test_lame_cylinder import (  # noqa: E402
    R_INNER as LAME_RI,
    R_OUTER as LAME_RE,
    _build_polar_mesh as _build_lame_mesh,
    _solve_and_collect_gauss_stress as _solve_lame,
    _lame_stress_polar,
)
from tests.validation.test_nafems_le1 import (  # noqa: E402
    SIGMA_YY_AT_D_REF,
    _apply_outer_pressure as _apply_nafems_pressure,
    _build_elliptic_mesh as _build_nafems_mesh,
    _sigma_yy_at_nearest_gauss,
    E_YOUNG as E_NAFEMS,
    NU as NU_NAFEMS,
)
from tests.validation.test_slender_cantilever import (  # noqa: E402
    U_TIP_TIMOSHENKO,
    _apply_tip_shear,
    _build_solid_cantilever,
    P_TIP,
    E_YOUNG as E_MAC,
    NU as NU_MAC,
)
from tests.validation.test_hill_cylinder_j2 import (  # noqa: E402
    K_SHEAR,
    R_INNER as HILL_RI,
    R_OUTER as HILL_RE,
    _hill_stress_polar,
    _hill_transition_radius,
    _solve_pressure as _solve_hill,
    _gather_gauss_state as _gather_hill,
)


# ---------------------------------------------------------------------------
# Subplots
# ---------------------------------------------------------------------------

def plot_lame(ax) -> None:
    """(a) Cilindro grueso de Lamé — σ_rr y σ_θθ vs r."""
    material = Elastic2D(E=1000.0, nu=0.3, hypothesis='plane_strain')
    domain, *_ = _build_lame_mesh(Quad8, nr=4, nt=4, material=material)
    _, _, srr_g, stt_g, rg, _ = _solve_lame(domain)

    r_ref = np.linspace(LAME_RI, LAME_RE, 200)
    srr_ref, stt_ref = _lame_stress_polar(r_ref)

    ax.plot(r_ref, srr_ref, 'b-', lw=1.5, label=r'$\sigma_{rr}$ (Lamé)')
    ax.plot(r_ref, stt_ref, 'r-', lw=1.5, label=r'$\sigma_{\theta\theta}$ (Lamé)')
    ax.scatter(rg, srr_g, c='b', s=14, marker='o',
               edgecolors='white', linewidths=0.5, label='FE Gauss', zorder=3)
    ax.scatter(rg, stt_g, c='r', s=14, marker='s',
               edgecolors='white', linewidths=0.5, zorder=3)
    ax.axhline(0.0, color='gray', lw=0.5, ls=':')
    ax.set_xlabel(r'$r$')
    ax.set_ylabel(r'Stress')
    ax.set_title('(a) Lamé thick cylinder')
    ax.legend(loc='center right', fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3)


def plot_nafems_le1(ax) -> None:
    """(b) NAFEMS LE1 — convergencia h de σ_yy(D) hacia 92.7."""
    mesh_sizes = [4, 8, 16, 32]
    sigma_q4 = []
    for n in mesh_sizes:
        material = Elastic2D(E=E_NAFEMS, nu=NU_NAFEMS, hypothesis='plane_stress')
        domain, outer_edges, _, _, _, point_D = _build_nafems_mesh(
            Quad4, nr=n, nt=n, material=material,
        )
        F = _apply_nafems_pressure(domain, outer_edges)
        U = LinearSolver(Assembler(domain)).solve(F)
        sigma_q4.append(_sigma_yy_at_nearest_gauss(point_D, U, domain))

    mesh_sizes_q8 = [2, 4, 8]
    sigma_q8 = []
    for n in mesh_sizes_q8:
        material = Elastic2D(E=E_NAFEMS, nu=NU_NAFEMS, hypothesis='plane_stress')
        domain, outer_edges, _, _, _, point_D = _build_nafems_mesh(
            Quad8, nr=n, nt=n, material=material,
        )
        F = _apply_nafems_pressure(domain, outer_edges)
        U = LinearSolver(Assembler(domain)).solve(F)
        sigma_q8.append(_sigma_yy_at_nearest_gauss(point_D, U, domain))

    ax.axhline(SIGMA_YY_AT_D_REF, color='k', lw=1.0, ls='--',
               label=fr'NAFEMS ref = {SIGMA_YY_AT_D_REF}')
    ax.plot(mesh_sizes, sigma_q4, 'o-', color='C0', lw=1.5,
            markersize=7, label='Quad4')
    ax.plot(mesh_sizes_q8, sigma_q8, 's-', color='C1', lw=1.5,
            markersize=7, label='Quad8')
    ax.set_xscale('log', base=2)
    ax.set_xlabel(r'mesh size $n_r = n_t$')
    ax.set_ylabel(r'$\sigma_{yy}$ at point $D$')
    ax.set_title('(b) NAFEMS LE1 — convergence to benchmark')
    ax.legend(loc='lower right', fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3, which='both')


def plot_macneal(ax) -> None:
    """(c) MacNeal-Harder slender cantilever — shear locking de Q4 vs Q8."""
    material = Elastic2D(E=E_MAC, nu=NU_MAC, hypothesis='plane_stress')

    cases = [
        ('Quad4 30×1\n(shear locking)', Quad4, 30, 1, 'C0', 'o'),
        ('Quad4 60×4', Quad4, 60, 4, 'C1', 's'),
        ('Quad4 120×8', Quad4, 120, 8, 'C2', '^'),
        ('Quad8 12×1', Quad8, 12, 1, 'C3', 'D'),
    ]
    labels, u_ratios = [], []
    for label, cls, nx, ny, _, _ in cases:
        material = Elastic2D(E=E_MAC, nu=NU_MAC, hypothesis='plane_stress')
        domain, mid_right, right_elems = _build_solid_cantilever(
            cls, nx=nx, ny=ny, material=material,
        )
        F = _apply_tip_shear(domain, right_elems, total_force=-P_TIP)
        U = LinearSolver(Assembler(domain)).solve(F)
        u_tip = abs(U[mid_right.dofs['uy']])
        labels.append(label)
        u_ratios.append(u_tip / U_TIP_TIMOSHENKO)

    colors = [c[4] for c in cases]
    x = np.arange(len(labels))
    ax.bar(x, u_ratios, color=colors, edgecolor='black', linewidth=0.5)
    ax.axhline(1.0, color='k', ls='--', lw=1.0,
               label='Timoshenko exact')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel(r'$u_{\rm tip} / u_{\rm Timoshenko}$')
    ax.set_title('(c) MacNeal-Harder slender beam')
    ax.set_ylim(0.0, 1.15)
    ax.legend(loc='lower right', fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3, axis='y')


def plot_hill(ax) -> None:
    """(d) Cilindro de Hill J2 — σ_rr y σ_θθ vs r en régimen elastoplástico."""
    p = 0.70
    r_c = _hill_transition_radius(p)
    domain, U, *_ = _solve_hill(p, num_steps=10)
    r_g, srr_g, stt_g, _ = _gather_hill(domain, U)

    r_ref = np.linspace(HILL_RI, HILL_RE, 300)
    srr_ref = np.array([_hill_stress_polar(r, p, r_c)[0] for r in r_ref])
    stt_ref = np.array([_hill_stress_polar(r, p, r_c)[1] for r in r_ref])

    ax.axvspan(HILL_RI, r_c, alpha=0.12, color='red', label='plastic zone')
    ax.axvspan(r_c, HILL_RE, alpha=0.10, color='blue', label='elastic zone')
    ax.plot(r_ref, srr_ref, 'b-', lw=1.5, label=r'$\sigma_{rr}$ (Hill)')
    ax.plot(r_ref, stt_ref, 'r-', lw=1.5, label=r'$\sigma_{\theta\theta}$ (Hill)')
    ax.scatter(r_g, srr_g, c='b', s=10, marker='o',
               edgecolors='white', linewidths=0.4, zorder=3)
    ax.scatter(r_g, stt_g, c='r', s=10, marker='s',
               edgecolors='white', linewidths=0.4, zorder=3)
    ax.axvline(r_c, color='k', ls=':', lw=0.8)
    ax.annotate(fr'$r_c \approx {r_c:.2f}$',
                xy=(r_c, 0.0), xycoords='data',
                xytext=(6, 8), textcoords='offset points',
                fontsize=8, ha='left', va='bottom')
    ax.set_xlabel(r'$r$')
    ax.set_ylabel(r'Stress')
    ax.set_title(fr'(d) Hill $J_2$ cylinder, $p={p}$')
    ax.legend(loc='lower right', fontsize=7.5, framealpha=0.9, ncol=2)
    ax.grid(True, alpha=0.3)


# ---------------------------------------------------------------------------
# Mosaico
# ---------------------------------------------------------------------------

def make_mosaic(output_path: str) -> None:
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 9.0,
        'axes.titlesize': 10.0,
        'axes.labelsize': 9.0,
        'legend.fontsize': 8.0,
    })
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.0), constrained_layout=True)
    plot_lame(axes[0, 0])
    plot_nafems_le1(axes[0, 1])
    plot_macneal(axes[1, 0])
    plot_hill(axes[1, 1])
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Wrote: {output_path}")


if __name__ == '__main__':
    out = os.path.join(os.path.dirname(__file__), 'validation_mosaic.png')
    make_mosaic(out)
