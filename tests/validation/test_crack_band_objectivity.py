"""Objetividad de malla de la banda de grieta (``IsotropicDamageCrackBand2D/3D``).

Aceptación ``mesh_objectivity`` de las specs: barra en tracción con una
columna (capa) de elementos debilitada (f_t·0.95) en el centro, celdas
cuadradas (cúbicas) de lado Δ = 20, 10, 5 mm, desplazamiento impuesto en un
extremo. La deformación se localiza en la columna debilitada, que disipa
``G_f/h`` por unidad de volumen; con ``h = Δ`` la energía disipada al final
es ``G_f·A_grieta`` con cualquier malla y la rama posterior al pico es la
misma.

- ν = 0: D/(G_f·A) = 1 a ±1 % en los diez sólidos y la carga en u = 0.05 mm
  igual entre mallas a ±1 %.
- ν = 0.2 (Quad4): D/(G_f·A) = 1 − ν² a ±1 %. Las columnas vecinas, que se
  descargan, impiden la contracción lateral de la banda: la banda trabaja con
  ε_yy ≈ 0, Rankine ve ε_eq = ε_xx/(1−ν²) y el trabajo de σ_xx en la banda
  es (1−ν²)·G_f/h. Es física del modelo continuo (igual en todas las mallas),
  no un error de discretización.
- Factor de ``h = c·Ω_e^(1/d)`` (regla Q2 de la spec: si la barra se aparta
  más de un 5 % de G_f·A, se corrige el factor):
  - triángulos y tetraedros: con el ``√A`` / ``∛V`` de Abaqus la columna
    disipaba 1.41 / 1.81·G_f·A (los dos triángulos o seis tetraedros de un
    cuadrado o cubo forman la banda); con √(2A) (Rots 1988) y ∛(6V), 1.00
    (decisión del usuario 2026-11-11, Q2);
  - cuadráticos: con el ½ de Abaqus disipaban 1.97 (``Quad8``, ``Quad9``,
    ``Hex20``, ``Hex27``), 2.68 (``Tri6``) y 3.23 (``Tet10``)·G_f·A: la
    deformación se localiza en el elemento entero, no en media fila de puntos
    de integración. Usan el factor del elemento lineal de su celda y dan 1.00
    (usuario, 2026-11-11, tras el diagnóstico de la viga SEN).

Unidades: N, mm, MPa. E = 30 000 MPa, f_t = 3 MPa, G_f = 0.1 N/mm.
"""
import logging

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.elements.solid_2d import Quad4, Quad8, Quad9, Tri3, Tri6
from solidum.elements.solid_3d import Hex8, Hex20, Hex27, Tet4, Tet10
from solidum.materials.damage_crack_band_2d import IsotropicDamageCrackBand2D
from solidum.materials.damage_crack_band_3d import IsotropicDamageCrackBand3D
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver
from solidum.mesh import structured_mesh

E, FT, GF = 30000.0, 3.0, 0.1
L, H, T = 100.0, 20.0, 10.0
U_MAX, STEPS = 0.4, 400


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


def bar(cls, nx, nu=0.0):
    """Barra de L × H (× H) con la columna central de celdas debilitada."""
    dim = len(cls.REFERENCE_NODES[0])
    dx = L / nx
    n = int(round(H / dx))
    if dim == 2:
        mat, weak = (IsotropicDamageCrackBand2D(E, nu, f, GF) for f in (FT, 0.95 * FT))
        kw = {"thickness": T}
        mesh = structured_mesh(cls, (nx, n), size=(L, H), material=mat, **kw)
        GA = GF * H * T
    else:
        mat, weak = (IsotropicDamageCrackBand3D(E, nu, f, GF) for f in (FT, 0.95 * FT))
        kw = {}
        mesh = structured_mesh(cls, (nx, n, n), size=(L, H, H), material=mat)
        GA = GF * H * H
    dom = mesh.domain
    x_weak = (nx // 2 + 0.5) * dx
    for eid, el in list(dom.elements.items()):
        xc = np.mean([nd.coordinates[0] for nd in el.nodes])
        if abs(xc - x_weak) < 0.5 * dx:
            dom.elements[eid] = cls(eid, el.nodes, weak, **kw)
    for nd in mesh.nodes_on("xmin"):
        nd.fix_dof("ux", 0.0)
    right = mesh.nodes_on("xmax")
    for nd in right:
        nd.fix_dof("ux", U_MAX)
    corners = {tuple(np.round(nd.coordinates, 9)): nd for nd in mesh.nodes_on("xmin")}
    origin = corners[(0.0,) * dim]
    origin.fix_dof("uy", 0.0)
    if dim == 3:
        origin.fix_dof("uz", 0.0)
        corners[(0.0, H, 0.0)].fix_dof("uz", 0.0)
    return _solve(dom, right), GA


def _check(results, expected_ratio):
    p05 = []
    for nx, ((curve, D), GA) in results.items():
        assert D / GA == pytest.approx(expected_ratio, rel=1e-2), nx
        # El pico es f_t·0.95·A; con Δu = 1 µm la curva lo muestrea por debajo.
        P_peak = 0.95 * FT * GA / GF
        assert 0.97 * P_peak < curve[:, 1].max() <= P_peak * (1 + 1e-9)
        p05.append(np.interp(0.05, curve[:, 0], curve[:, 1]))
    assert max(p05) / min(p05) - 1.0 < 1e-2


@pytest.mark.parametrize("cls", [Quad4, Tri3, Quad8, Quad9, Tri6], ids=lambda c: c.__name__)
def test_barra_2d_energia_objetiva(cls):
    _check({nx: bar(cls, nx) for nx in (5, 10, 20)}, 1.0)


def test_barra_2d_con_poisson_disipa_1_menos_nu2():
    nu = 0.2
    _check({nx: bar(Quad4, nx, nu) for nx in (5, 10, 20)}, 1.0 - nu**2)


@pytest.mark.parametrize("cls", [Hex8, Tet4, Hex20, Hex27, Tet10], ids=lambda c: c.__name__)
def test_barra_3d_energia_objetiva(cls):
    _check({nx: bar(cls, nx) for nx in (5, 10)}, 1.0)
