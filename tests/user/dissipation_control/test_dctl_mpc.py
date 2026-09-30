"""``DissipationArcLengthSolver`` con restricciones multipunto (deuda #31).

El solver hereda de ``ArcProblem`` la norma del residuo reducido ``‖TᵀR‖``
(ADR 0004); hasta el 2026-09-30 medía ``‖R[libres]‖``, que con una MPC no se
anula en la solución, y no convergía. Oráculo del nodo fundido: la barra con
daño localizado de ``test_dctl_retroceso.py`` con un nodo duplicado y atado por
una MPC (``ux_dup = ux``) recorre paso a paso el mismo retroceso que la barra
de malla conforme, incluido el tramo resuelto con la restricción de disipación.
Con MPC homogéneas los productos ``F_ref·ΔU`` de la restricción valen lo mismo
en los vectores completos que en los reducidos (``F·T·u = (TᵀF)·u``).
"""
from __future__ import annotations

import logging
import os
import sys

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
sys.path.insert(0, os.path.dirname(__file__))

import solidum  # noqa: E402
from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.truss import Truss2D  # noqa: E402
from solidum.materials.damage_1d import IsotropicDamage1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402

from _dctl_espia import convergidos, espiar_problemas_de_newton  # noqa: E402

solidum.load_user_module("dissipation_control")
from solidum.user.dissipation_control import DissipationArcLengthSolver  # noqa: E402

E, A, L = 30.0e9, 0.01, 1.0
K0, ALPHA, WEAK = 1.0e-4, 5000.0, 0.95
N_ELEM = 10
H_ELEM = L / N_ELEM
F_PEAK = A * E * WEAK * K0
NODO_ATADO = 3  # índice del nodo entre los elementos 3 y 4


def _barra(atada: bool):
    """Barra de 10 ``Truss2D`` con el elemento central debilitado. Atada: el
    elemento 4 empieza en un nodo propio coincidente con el nodo 3, cuyo ``ux``
    es esclavo del de ese nodo."""
    dom = Domain()
    nodos = [dom.add_node(i + 1, [i * H_ELEM, 0.0]) for i in range(N_ELEM + 1)]
    duplicado = dom.add_node(100, [NODO_ATADO * H_ELEM, 0.0]) if atada else None
    for k in range(N_ELEM):
        kappa_0 = WEAK * K0 if k == N_ELEM // 2 else K0
        a, b = nodos[k], nodos[k + 1]
        if atada and k == NODO_ATADO:
            a = duplicado
        dom.add_element(Truss2D(k + 1, [a, b], IsotropicDamage1D(E=E, kappa_0=kappa_0, alpha=ALPHA), A=A))
    nodos[0].fix_dof("ux", 0.0)
    for n in nodos:
        n.fix_dof("uy", 0.0)
    if atada:
        duplicado.fix_dof("uy", 0.0)
        dom.add_linear_constraint((duplicado.id, "ux"), [(nodos[NODO_ATADO].id, "ux")], [1.0])
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[nodos[-1].dofs["ux"]] = F_PEAK
    return dom, F, nodos[-1]


def _traza(atada: bool):
    dom, F, extremo = _barra(atada)
    solver = DissipationArcLengthSolver(Assembler(dom), max_iter=30, max_lambda=1.2, max_steps=120,
                                        initial_dlambda=0.02, initial_tau=1.0e-3)
    intentos = espiar_problemas_de_newton(solver)
    historia = []
    logging.disable(logging.CRITICAL)
    try:
        solver.solve(F, step_callback=lambda k, U, lam: historia.append((U[extremo.dofs["ux"]], lam)))
    finally:
        logging.disable(logging.NOTSET)
    return np.array(historia), convergidos(intentos, "dissipation")


def test_barra_atada_sigue_el_mismo_retroceso_que_la_fundida():
    fundida, n_dis_fundida = _traza(False)
    atada, n_dis_atada = _traza(True)
    assert n_dis_fundida > 50, "el control debe ejercitar la restricción de disipación"
    assert n_dis_atada == n_dis_fundida
    assert len(atada) == len(fundida)
    escala = np.abs(fundida).max(axis=0)
    np.testing.assert_allclose(atada / escala, fundida / escala, rtol=0.0, atol=1e-9)
    # El trazado retrocede: el desplazamiento del extremo baja tras el pico.
    i_pico = int(np.argmax(fundida[:, 1]))
    assert fundida[i_pico:, 0].min() < 0.7 * fundida[i_pico, 0]
