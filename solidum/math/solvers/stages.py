# solidum_fem/solidum/math/solvers/stages.py
"""Etapas de carga de los solvers estáticos incrementales (ADR 0026).

Una **etapa** parte de un estado de equilibrio ya alcanzado ``U₀`` (con el
estado interno comprometido de los elementos, que vive en el ensamblador) y
de una carga constante ``F_c`` que actúa completa durante toda la etapa,
mientras otra carga de referencia crece con el factor ``λ``, que vale 0 al
empezar:

    R(U, λ) = F_c + λ·F_ref − F_int(U).

Los valores impuestos (Dirichlet y término independiente de las restricciones
multipunto) crecen de lo que valían en ``U₀`` al valor declarado en el modelo,
``g(λ) = g₀ + λ·(g − g₀)``: un apoyo que no cambia entre etapas queda
constante y uno que cambia crece con la etapa. Sin estado inicial, ``g₀ = 0``
y todo es como en una corrida desde cero.

Es la operación de los pasos de carga de Abaqus y ANSYS y del ``loadConst``
de OpenSees: aplicar primero la gravedad, que se queda, y empujar después.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np

from solidum.math.solvers._shared import _log
from solidum.math.solvers.corrector import default_calibration_scales


@dataclass
class Stage:
    """Datos de partida de una etapa (ADR 0026).

    Attributes
    ----------
    U0 : np.ndarray
        Desplazamientos de partida (ceros si la etapa empieza desde cero).
    F_const : np.ndarray
        Carga constante de la etapa (ceros si no hay).
    g_start : np.ndarray or None
        Parte independiente de los valores impuestos en ``U0``,
        ``g₀ = U₀ − T·U₀[libres]``; ``None`` si es nula (``g(λ) = λ·g``).
    staged : bool
        ``True`` si se dio estado inicial o carga constante.
    """

    U0: np.ndarray
    F_const: np.ndarray
    g_start: np.ndarray | None
    staged: bool

    def calibration_load(self, F_ref: np.ndarray) -> np.ndarray:
        """Carga con la que se calibra el criterio de convergencia: la de la
        etapa completa, constante más referencia."""
        return self.F_const + F_ref if self.staged else F_ref


def make_stage(assembler, u0, constant_load, solver_name: str) -> Stage:
    """Valida y prepara el estado de partida de una etapa."""
    ndof = assembler.domain.total_dofs

    def as_vector(value, name):
        arr = np.asarray(value, dtype=float).reshape(-1)
        if arr.size != ndof:
            raise ValueError(f"{solver_name}: {name} tiene {arr.size} componentes y el "
                             f"modelo {ndof} grados de libertad.")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{solver_name}: {name} tiene valores no finitos.")
        return arr.copy()

    U0 = np.zeros(ndof) if u0 is None else as_vector(u0, "u0")
    F_const = np.zeros(ndof) if constant_load is None else as_vector(constant_load, "constant_load")
    g_start = None
    if u0 is not None:
        cs = assembler.constraint_set
        T, _ = cs.build(ndof)
        g0 = U0 - T @ U0[cs.free_dofs(ndof)]
        if np.any(g0 != 0.0):
            g_start = np.asarray(g0, dtype=float)
    return Stage(U0=U0, F_const=F_const, g_start=g_start,
                 staged=(u0 is not None or constant_load is not None))


def check_initial_equilibrium(assembler, stage: Stage, convergence, F_ref: np.ndarray,
                              solver_name: str) -> float | None:
    """Residuo relativo ``‖Tᵀ(F_c − F_int(U₀))‖ / ‖F‖`` al empezar.

    Se comprueba en una etapa (``u0`` o ``constant_load``) y en una corrida
    desde ``U = 0`` sobre un modelo que conserva el estado interno de un
    análisis anterior. En los dos casos la etapa debe empezar en equilibrio:
    si no, el estado de los elementos y los desplazamientos no corresponden y
    se avisa. Con el criterio de fuerza del solver (ADR 0007), calibrado con
    la carga de la etapa en una copia, para no alterar el del solver.
    Devuelve ``None`` si no hay nada que comprobar.
    """
    if not stage.staged and not getattr(assembler.domain, "committed_history", False):
        return None
    K, F_int = assembler.assemble_non_linear_system(stage.U0)
    criterion = copy.deepcopy(convergence)
    if not criterion.is_calibrated:
        criterion.calibrate(*default_calibration_scales(stage.calibration_load(F_ref), F_int, K))
    T = assembler.constraint_set.build(assembler.domain.total_dofs)[0]
    residual = float(np.linalg.norm(T.T @ (stage.F_const - F_int)))
    ref = max(float(np.linalg.norm(stage.F_const)), float(np.linalg.norm(F_int)))
    state = criterion.evaluate(residual_norm=residual, ref_force=ref, delta_u_norm=0.0,
                               u_norm=float(np.linalg.norm(stage.U0)))
    ratio = residual / ref if ref > 0.0 else 0.0
    if state.ratio_force > 1.0:
        if stage.staged:
            _log.warning(
                f"{solver_name}: la etapa empieza fuera de equilibrio "
                f"(‖F_c − F_int(u0)‖/‖F‖ = {ratio:.2e}). Lo habitual es continuar "
                f"desde el resultado de la etapa anterior: u0 = sus desplazamientos "
                f"y constant_load = su carga aplicada, con el mismo ensamblador "
                f"(guarda el estado interno de los elementos). solidum.run(..., "
                f"previous=resultado) lo hace así."
            )
        else:
            _log.warning(
                f"{solver_name}: el modelo conserva el estado interno de un análisis "
                f"anterior, pero esta corrida empieza en U = 0 sin carga, y ese estado "
                f"no está en equilibrio ahí (‖F_int(0)‖/‖F‖ = {ratio:.2e}). Para "
                f"continuar el análisis anterior, pase u0 y constant_load (o use "
                f"solidum.run(..., previous=resultado)); para un análisis nuevo, "
                f"construya el modelo de nuevo."
            )
    return ratio
