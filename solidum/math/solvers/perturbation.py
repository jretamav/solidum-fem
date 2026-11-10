# solidum_fem/solidum/math/solvers/perturbation.py
"""Perturbación lineal estática sobre el estado de una etapa (ADR 0031).

Un ``LinearSolver`` que continúa una etapa (``solidum.run(..., previous=…)``,
o una etapa ``stages`` del YAML) no avanza el estado: calcula la respuesta
lineal a las cargas nuevas con la rigidez tangente del estado de la etapa
(decisión del usuario del 2026-11-09),

    K_t(u_b)·ΔU = ΔF,        U = u_b + ΔU,

con los valores impuestos que cambien respecto de la etapa (``Δg``). Las
reacciones y las fuerzas internas son lineales en ``ΔU``: las del estado base
más ``K_t·ΔU``. El estado interno de los elementos no cambia, así que varios
casos lineales sobre la misma etapa son superponibles y la etapa siguiente
continúa la etapa base. Es el caso de carga lineal sobre el estado P-Δ de
ETABS y SAP2000 y el paso de perturbación lineal de Abaqus.
"""
from __future__ import annotations

import numpy as np

from solidum.core.load_vector import as_load_vector, equivalent_loads_of
from solidum.math.linalg import StiffnessProperties, select_solver
from solidum.math.solvers._shared import _log, domain_is_symmetric
from solidum.math.solvers.stages import make_stage


def _tangent_and_force(elem, u_e: np.ndarray):
    """``(K_e, f_e)`` del elemento en ``u_e`` sin dejar huella en su estado
    trial: es la evaluación de un estado ya comprometido, no un iterado."""
    state = getattr(elem, "state", None)
    if state is not None and hasattr(state, "snapshot_trial"):
        snapshot = state.snapshot_trial()
        try:
            return elem.compute_element_state(u_e)
        finally:
            state.restore_trial(snapshot)
    return elem.compute_element_state(u_e)


def linear_perturbation(domain, assembler, solver, F_applied, previous):
    """``SolveResult`` de la perturbación lineal de ``F_applied`` sobre el
    estado de ``previous`` (ver el docstring del módulo)."""
    from solidum.results import build_solve_result

    name = type(solver).__name__
    stage = make_stage(assembler, previous.U, previous.F_applied, name)
    u_b = stage.U0
    K_t, F_int_b = assembler.assemble_non_linear_system(u_b)
    dF = np.asarray(F_applied, dtype=float).reshape(-1)

    # Δg: los valores impuestos del modelo menos los del estado base.
    K_red, F_red, T, dg = assembler.reduce(K_t, dF, g_start=stage.g_start)
    props = StiffnessProperties(
        is_symmetric=domain_is_symmetric(domain), is_positive_definite=False,
        size=K_red.shape[0], near_nullspace=assembler.near_nullspace)
    lin = select_solver(props, override=getattr(solver, "linear_algebra", "auto"))
    du_red = lin.solve(K_red.tocsr(), F_red)
    dU = assembler.expand(du_red, T, dg)
    U = u_b + dU
    _log.info(f"{name}: perturbación lineal sobre el estado de la etapa anterior "
              f"(ADR 0031), rigidez tangente; el estado no cambia.")

    # Cargas: la constante de la etapa (con sus cargas equivalentes) más la nueva.
    F_constant = as_load_vector(previous.F_applied, previous.equivalent_loads)
    F_new = F_applied if equivalent_loads_of(F_applied) is not None else as_load_vector(dF, {})
    F_total = F_constant + F_new
    carried = equivalent_loads_of(F_total)
    eq_loads = {} if carried is None else carried.as_dict()

    # Fuerzas lineales: las del estado base más K_t·ΔU.
    F_int_lin = F_int_b + K_t @ dU
    element_forces = {}
    for elem_id, elem in domain.elements.items():
        if not _defines_from_nodal(elem):
            continue
        idx = elem.get_global_dof_indices()
        K_e, f_e = _tangent_and_force(elem, u_b[idx])
        ef = elem.internal_forces_from_nodal(f_e + K_e @ dU[idx], u_b, eq_loads.get(elem_id))
        if ef is not None:
            element_forces[elem_id] = ef

    return build_solve_result(
        domain, assembler, U, F_total,
        F_constant=F_constant, equivalent_loads=eq_loads,
        F_int=F_int_lin, element_forces=element_forces, perturbation=True,
        condition_estimate=None,
    )


def _defines_from_nodal(elem) -> bool:
    """``True`` si la clase del elemento redefine la conversión de fuerzas
    nodales a fuerzas internas (los elementos estructurales 1D)."""
    from solidum.core.element import Element
    return type(elem).internal_forces_from_nodal is not Element.internal_forces_from_nodal
