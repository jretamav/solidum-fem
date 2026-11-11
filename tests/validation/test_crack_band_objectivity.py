"""Objetividad de malla de la banda de grieta (``IsotropicDamageCrackBand2D/3D``).

Aceptación ``mesh_objectivity`` de las specs: barra en tracción con una
columna (capa) de elementos debilitada (f_t·0.95) en el centro, elementos
cuadrados (cúbicos) de lado Δ = 20, 10, 5 mm, desplazamiento impuesto en un
extremo. La deformación se localiza en la columna debilitada, que disipa
``G_f/h`` por unidad de volumen; con ``h = Δ`` la energía disipada al final
es ``G_f·A_grieta`` con cualquier malla y la rama posterior al pico es la
misma.

- ν = 0: D/(G_f·A) = 1 a ±1 % (Quad4, Tri3, Hex8, Tet4) y la carga en
  u = 0.05 mm igual entre mallas a ±1 %.
- ν = 0.2 (Quad4): D/(G_f·A) = 1 − ν² a ±1 %. Las columnas vecinas, que se
  descargan, impiden la contracción lateral de la banda: la banda trabaja con
  ε_yy ≈ 0, Rankine ve ε_eq = ε_xx/(1−ν²) y el trabajo de σ_xx en la banda
  es (1−ν²)·G_f/h. Es física del modelo continuo (igual en todas las mallas),
  no un error de discretización.
- Triángulos y tetraedros: con el ``√A`` / ``∛V`` de Abaqus la columna
  disipaba 1.41 / 1.81·G_f·A (los dos triángulos o seis tetraedros de un
  cuadrado o cubo forman la banda); con √(2A) (Rots 1988) y ∛(6V), 1.00
  (decisión del usuario 2026-11-11, Q2).

Unidades: N, mm, MPa. E = 30 000 MPa, f_t = 3 MPa, G_f = 0.1 N/mm.
"""
import logging

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4, Tri3
from solidum.elements.solid_3d import Hex8, Tet4
from solidum.materials.damage_crack_band_2d import IsotropicDamageCrackBand2D
from solidum.materials.damage_crack_band_3d import IsotropicDamageCrackBand3D
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver

E, FT, GF = 30000.0, 3.0, 0.1
L, H, T = 100.0, 20.0, 10.0
U_MAX, STEPS = 0.4, 400
KUHN = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]


def _solve(dom, right_nodes, u_max=U_MAX, steps=STEPS):
    dom.generate_equation_numbers(verbose=False)
    right = [nd.dofs["ux"] for nd in right_nodes]
    asm = Assembler(dom)
    curve = [(0.0, 0.0)]

    def cb(step, U, lam):
        _, F = asm.assemble_non_linear_system(U)
        curve.append((U[right[0]], F[right].sum()))

    level = logging.getLogger("solidum").level
    logging.getLogger("solidum").setLevel(logging.ERROR)
    try:
        NonlinearSolver(asm, num_steps=steps).solve(np.zeros(dom.total_dofs), step_callback=cb)
    finally:
        logging.getLogger("solidum").setLevel(level)
    c = np.array(curve)
    W = np.trapezoid(c[:, 1], c[:, 0])
    return c, W - 0.5 * c[-1, 1] * c[-1, 0]


def bar_2d(kind, nx, nu):
    dx = L / nx
    ny = int(round(H / dx))
    dom = Domain()
    mat = IsotropicDamageCrackBand2D(E, nu, FT, GF)
    weak = IsotropicDamageCrackBand2D(E, nu, 0.95 * FT, GF)
    nodes = {(i, j): None for j in range(ny + 1) for i in range(nx + 1)}
    for k, (i, j) in enumerate(nodes):
        nodes[(i, j)] = dom.add_node(k + 1, [i * dx, j * H / ny])
    eid = 0
    for j in range(ny):
        for i in range(nx):
            m = weak if i == nx // 2 else mat
            a, b, c, d = nodes[(i, j)], nodes[(i + 1, j)], nodes[(i + 1, j + 1)], nodes[(i, j + 1)]
            if kind == "Quad4":
                eid += 1
                dom.add_element(Quad4(eid, [a, b, c, d], m, thickness=T))
            else:
                for tri in ([a, b, c], [a, c, d]):
                    eid += 1
                    dom.add_element(Tri3(eid, tri, m, thickness=T))
    for j in range(ny + 1):
        nodes[(0, j)].fix_dof("ux", 0.0)
        nodes[(nx, j)].fix_dof("ux", U_MAX)
    nodes[(0, 0)].fix_dof("uy", 0.0)
    return _solve(dom, [nodes[(nx, j)] for j in range(ny + 1)]), GF * H * T


def bar_3d(kind, nx):
    dx = L / nx
    n = int(round(H / dx))
    dom = Domain()
    mat = IsotropicDamageCrackBand3D(E, 0.0, FT, GF)
    weak = IsotropicDamageCrackBand3D(E, 0.0, 0.95 * FT, GF)
    nodes = {}
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(nx + 1):
                nodes[(i, j, k)] = dom.add_node(len(nodes) + 1, [i * dx, j * dx, k * dx])
    eid = 0
    for k in range(n):
        for j in range(n):
            for i in range(nx):
                m = weak if i == nx // 2 else mat
                c = [nodes[(i, j, k)], nodes[(i + 1, j, k)], nodes[(i + 1, j + 1, k)],
                     nodes[(i, j + 1, k)], nodes[(i, j, k + 1)], nodes[(i + 1, j, k + 1)],
                     nodes[(i + 1, j + 1, k + 1)], nodes[(i, j + 1, k + 1)]]
                if kind == "Hex8":
                    eid += 1
                    dom.add_element(Hex8(eid, c, m))
                    continue
                for t in KUHN:
                    tn = [c[q] for q in t]
                    X = np.array([nd.coordinates for nd in tn])
                    if np.linalg.det(np.array([X[1] - X[0], X[2] - X[0], X[3] - X[0]])) < 0:
                        tn = [tn[0], tn[2], tn[1], tn[3]]
                    eid += 1
                    dom.add_element(Tet4(eid, tn, m))
    for k in range(n + 1):
        for j in range(n + 1):
            nodes[(0, j, k)].fix_dof("ux", 0.0)
            nodes[(nx, j, k)].fix_dof("ux", U_MAX)
    nodes[(0, 0, 0)].fix_dof("uy", 0.0)
    nodes[(0, 0, 0)].fix_dof("uz", 0.0)
    nodes[(0, n, 0)].fix_dof("uz", 0.0)
    right = [nodes[(nx, j, k)] for k in range(n + 1) for j in range(n + 1)]
    return _solve(dom, right), GF * H * H


def _check(results, expected_ratio):
    p05 = []
    for nx, ((curve, D), GA) in results.items():
        assert D / GA == pytest.approx(expected_ratio, rel=1e-2), nx
        # El pico es f_t·0.95·A; con Δu = 1 µm la curva lo muestrea por debajo.
        P_peak = 0.95 * FT * GA / GF
        assert 0.97 * P_peak < curve[:, 1].max() <= P_peak * (1 + 1e-9)
        p05.append(np.interp(0.05, curve[:, 0], curve[:, 1]))
    assert max(p05) / min(p05) - 1.0 < 1e-2


@pytest.mark.parametrize("kind", ["Quad4", "Tri3"])
def test_barra_2d_energia_objetiva(kind):
    _check({nx: bar_2d(kind, nx, 0.0) for nx in (5, 10, 20)}, 1.0)


def test_barra_2d_con_poisson_disipa_1_menos_nu2():
    nu = 0.2
    _check({nx: bar_2d("Quad4", nx, nu) for nx in (5, 10, 20)}, 1.0 - nu**2)


@pytest.mark.parametrize("kind", ["Hex8", "Tet4"])
def test_barra_3d_energia_objetiva(kind):
    _check({nx: bar_3d(kind, nx) for nx in (5, 10)}, 1.0)
