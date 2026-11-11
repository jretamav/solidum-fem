"""Ensayo de tracción de Van Vliet con fisuración distribuida (banda de grieta).

Aceptación ``structural: van_vliet_tension`` de la spec
``IsotropicDamageCrackBand2D``: la probeta en hueso de perro de la tesis
(Retama 2010, §8.1; geometría y datos como en
``tests/user/discontinuities/_van_vliet.py``), ahora con TODO el concreto de
``IsotropicDamageCrackBand2D`` (Rankine, ``Tri3``): la grieta no se prescribe,
se forma donde el esfuerzo alcanza ``f_t``.

- Dos mallas (filas de 8 y 4 mm junto al cuello). Carga máxima y energía
  disipada en ``u = 0.12 mm`` iguales entre mallas a ±3 %.
- Al final (``u = 1.2 mm``, la carga ya es < 0.5 % del pico), energía
  disipada = ``(1−ν²)·G_f·A`` a ±3 %, con ``A`` el ligamento del cuello: la
  banda queda confinada lateralmente y disipa ``(1−ν²)·G_f`` por unidad de
  área (decisión del usuario 2026-11-11: documentado, sin corregir).
- IMPL-EX con paso fino (0.0625 µm): carga máxima y energía a ±2 % de la
  implícita. Con 0.5 µm la extrapolación sobrepasa el pico un 7.7 %; el
  error baja con el paso (5.1, 2.9 y 0.8 % con 0.25, 0.125 y 0.0625 µm).

Supuesto propio: las placas de acero son dos capas de 25 mm. Con dos capas de
8 mm (las de la prueba de la discontinuidad embebida, donde sólo puede
agrietarse la banda del cuello) la carga puntual flexiona las placas y el
concreto junto a ellas llega a ``σ_I`` 1.7 veces mayor que en el cuello:
con fisuración en todo el concreto la probeta se rompía junto a las placas.
La curva de la tesis (pico 27.0 kN, discontinuidad embebida) sólo es
referencia: aquí 26.5–26.9 kN.

Unidades: N, m, Pa.
"""
from __future__ import annotations

import logging

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.core.domain import Domain
from solidum.elements.solid_2d import Tri3
from solidum.materials.damage_crack_band_2d import IsotropicDamageCrackBand2D
from solidum.materials.elastic_2d import Elastic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver

D, RADIO, PROFUNDIDAD, ESPESOR = 0.200, 0.145, 0.040, 0.100
E_C, NU_C = 39.8e9, 0.2
F_T, G_F = 2.57e6, 121.9
E_ACERO, NU_ACERO = 210.0e9, 0.3
CAPA_PLACA = 0.025
EXCENTRICIDAD = D / 50
CUELLO = 0.120
GFA = G_F * CUELLO * ESPESOR
_XC = D / 2 - PROFUNDIDAD + RADIO


def semiancho(y: float) -> float:
    return _XC - np.sqrt(RADIO**2 - y**2) if abs(y) <= D / 2 else D / 2


def modelo(fino: bool, integration: str):
    """Malla estructurada mapeada al contorno de las entallas, con una fila
    de elementos centrada en el cuello. Devuelve el dominio y las rótulas
    superior e inferior."""
    if fino:
        hb = 0.004
        pos = ([hb / 2 + hb * k for k in range(15)] + [0.060, 0.065]
               + [0.070 + 0.005 * k for k in range(18)])
        n_lado = 14
    else:
        hb = 0.008
        pos = [hb / 2 + hb * k for k in range(8)] + [0.070 + 0.010 * k for k in range(9)]
        n_lado = 7
    ys = sorted([-y for y in pos] + pos)
    paso = (1 - 0.04) / n_lado
    xi_pos = [0.04 + paso * k for k in range(n_lado + 1)]
    xis = sorted([-x for x in xi_pos] + [0.0] + xi_pos)

    dom = Domain()
    concreto = IsotropicDamageCrackBand2D(E_C, NU_C, F_T, G_F, integration=integration)
    acero = Elastic2D(E=E_ACERO, nu=NU_ACERO, hypothesis="plane_stress")
    nodos: dict = {}

    def nodo(x, y):
        clave = (round(x, 9), round(y, 9))
        if clave not in nodos:
            nodos[clave] = dom.add_node(len(nodos) + 1, [x, y])
        return nodos[clave]

    def triangulo(a, b, c, material):
        dom.add_element(Tri3(len(dom.elements) + 1, [a, b, c], material, thickness=ESPESOR))

    def fila(xs0, y0, xs1, y1, material):
        for i in range(len(xs0) - 1):
            a, b, c, d = nodo(xs0[i], y0), nodo(xs0[i + 1], y0), nodo(xs1[i + 1], y1), nodo(xs1[i], y1)
            if xs0[i] + xs0[i + 1] < 0:
                triangulo(a, b, c, material)
                triangulo(a, c, d, material)
            else:
                triangulo(a, b, d, material)
                triangulo(b, c, d, material)

    for j in range(len(ys) - 1):
        y0, y1 = ys[j], ys[j + 1]
        fila([x * semiancho(y0) for x in xis], y0, [x * semiancho(y1) for x in xis], y1, concreto)
    xs_placa = sorted([-0.6 * D] + [x * D / 2 for x in xis] + [0.6 * D])
    borde = 0.75 * D
    for capa in range(2):
        for signo in (1, -1):
            y0, y1 = borde + capa * CAPA_PLACA, borde + (capa + 1) * CAPA_PLACA
            y0, y1 = (y0, y1) if signo > 0 else (-y1, -y0)
            fila(xs_placa, y0, xs_placa, y1, acero)
    tope = borde + 2 * CAPA_PLACA
    return dom, nodo(-EXCENTRICIDAD, tope), nodo(-EXCENTRICIDAD, -tope)


def ensayo(fino: bool, integration: str, u_max: float, paso: float) -> dict:
    """Desplazamiento impuesto en la rótula superior; devuelve la carga
    máxima [kN] y la energía disipada [J] = W − ½·F·u (descarga secante)."""
    dom, arriba, abajo = modelo(fino, integration)
    abajo.fix_dof("ux", 0.0)
    abajo.fix_dof("uy", 0.0)
    arriba.fix_dof("ux", 0.0)
    arriba.fix_dof("uy", u_max)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    curva = [(0.0, 0.0)]

    def registrar(step, U, lam):
        _, F_int = asm.assemble_non_linear_system(U)
        curva.append((U[arriba.dofs["uy"]], F_int[arriba.dofs["uy"]]))

    nivel = logging.getLogger("solidum").level
    logging.getLogger("solidum").setLevel(logging.ERROR)
    try:
        NonlinearSolver(asm, num_steps=int(round(u_max / paso))).solve(
            np.zeros(dom.total_dofs), step_callback=registrar)
    finally:
        logging.getLogger("solidum").setLevel(nivel)
    c = np.array(curva)
    W = np.trapezoid(c[:, 1], c[:, 0])
    return {"P_max": c[:, 1].max() / 1e3, "P_fin": c[-1, 1] / 1e3,
            "D": W - 0.5 * c[-1, 1] * c[-1, 0]}


@pytest.fixture(scope="module")
def implicita():
    return {fino: ensayo(fino, "implicit", 0.12e-3, 1.0e-6) for fino in (False, True)}


def test_carga_maxima_y_energia_entre_mallas(implicita):
    g, f = implicita[False], implicita[True]
    assert f["P_max"] == pytest.approx(g["P_max"], rel=0.03)
    assert f["D"] == pytest.approx(g["D"], rel=0.03)
    assert 25.0 < g["P_max"] < 28.5          # tesis (embebida): 27.0 kN, sólo referencia


def test_energia_al_final_igual_a_la_del_ligamento():
    r = ensayo(False, "implicit", 1.2e-3, 2.0e-6)
    assert r["P_fin"] < 0.005 * r["P_max"]
    assert r["D"] / ((1.0 - NU_C**2) * GFA) == pytest.approx(1.0, abs=0.03)


def test_impl_ex_con_paso_fino_igual_a_la_implicita(implicita):
    r = ensayo(False, "impl-ex", 0.12e-3, 0.0625e-6)
    ref = implicita[False]
    assert r["P_max"] == pytest.approx(ref["P_max"], rel=0.02)
    assert r["D"] == pytest.approx(ref["D"], rel=0.02)
