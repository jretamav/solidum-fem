"""Entrypoints públicos de Fenix FEM (ADR 0002).

Encapsulan el pipeline completo ``Assembler → solver.solve → SolveResult`` para
consumidores externos (GUIs, scripts). Retornan ``SolveResult`` inmutable y lo
cachean en ``domain.last_result``.

- :func:`run`: recibe un ``Domain`` ya construido (y opcionalmente solver y
  cargas) y resuelve. Es la vía para modelos construidos programáticamente.
- :func:`run_modal`: análoga para análisis modal (ADR 0009). Devuelve
  :class:`ModalResult` en vez de :class:`SolveResult` — son análisis de
  naturaleza distinta y no se comparte el pipeline.
- :func:`run_yaml`: parsea un archivo YAML, arma el modelo y despacha a
  :func:`run` o :func:`run_modal` según el ``type`` de solver declarado.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Optional, Union

import numpy as np

from fenix.core.domain import Domain
from fenix.math.assembly import Assembler
from fenix.math.solvers import LinearSolver, ModalSolver
from fenix.results import ModalResult, SolveResult, build_solve_result


StepCallback = Callable[[int, np.ndarray, float], None]


def _invoke_solve(solver: Any, F_applied: np.ndarray,
                   step_callback: Optional[StepCallback]) -> np.ndarray:
    """Llama ``solver.solve`` pasando ``step_callback`` solo si el solver lo acepta.

    ``LinearSolver.solve`` no lo admite; los no-lineales sí. Detectar por firma
    en tiempo de llamada evita forzar el callback donde no aplica.
    """
    if step_callback is None:
        return solver.solve(F_applied)
    try:
        return solver.solve(F_applied, step_callback=step_callback)
    except TypeError:
        return solver.solve(F_applied)


def run(
    domain: Domain,
    *,
    assembler: Optional[Assembler] = None,
    solver: Optional[Any] = None,
    F_applied: Optional[np.ndarray] = None,
    step_callback: Optional[StepCallback] = None,
) -> SolveResult:
    """Ejecuta el pipeline completo y retorna el resultado agregado.

    Parameters
    ----------
    domain
        Dominio con nodos, elementos y BCs ya configurados. Si ``total_dofs==0``
        se invoca ``generate_equation_numbers`` automáticamente.
    assembler
        Assembler opcional. Si es ``None`` se construye uno nuevo sobre ``domain``.
    solver
        Instancia de solver ya vinculada al assembler. Si es ``None`` se usa un
        ``LinearSolver`` por defecto.
    F_applied
        Vector global ``(n_dof,)`` de cargas aplicadas. Si es ``None`` se usa
        ceros (útil para problemas controlados por desplazamientos).
    step_callback
        Callback opcional ``(step, U, load_factor) -> None`` para solvers
        incrementales. Se ignora en solvers que no lo admiten.

    Returns
    -------
    SolveResult
        Resultado inmutable. También queda en ``domain.last_result``.
    """
    if domain.total_dofs == 0:
        domain.generate_equation_numbers()
    if assembler is None:
        assembler = Assembler(domain)
    if solver is None:
        solver = LinearSolver(assembler)
    if F_applied is None:
        F_applied = np.zeros(domain.total_dofs)

    U = _invoke_solve(solver, F_applied, step_callback)

    result = build_solve_result(domain, assembler, U, F_applied)
    domain.last_result = result
    return result


def run_modal(
    domain: Domain,
    *,
    assembler: Optional[Assembler] = None,
    solver: Optional[ModalSolver] = None,
    n_modes: Optional[int] = None,
) -> ModalResult:
    """Ejecuta un análisis modal y retorna el resultado.

    Parameters
    ----------
    domain
        Dominio con nodos, elementos y BCs ya configurados. Si
        ``total_dofs == 0`` se invoca ``generate_equation_numbers``
        automáticamente.
    assembler
        Assembler opcional. Si es ``None`` se construye sobre ``domain``.
    solver
        Instancia de ``ModalSolver`` ya vinculada al assembler. Si es
        ``None`` se construye con ``ModalSolver(assembler, n_modes=n_modes)``
        y los demás parámetros por defecto.
    n_modes
        Número de modos a calcular cuando ``solver is None``. Requerido en
        ese caso; ignorado si ``solver`` se pasa explícitamente.

    Returns
    -------
    ModalResult
        Resultado inmutable. Queda también en ``domain.last_result``.
    """
    if domain.total_dofs == 0:
        domain.generate_equation_numbers()
    if assembler is None:
        assembler = Assembler(domain)
    if solver is None:
        if n_modes is None:
            raise ValueError(
                "run_modal: pasar `solver` o `n_modes` (no ambos None)."
            )
        solver = ModalSolver(assembler, n_modes=n_modes)

    result = solver.solve()
    domain.last_result = result
    return result


def run_yaml(
    path: str | Path,
    *,
    step_callback: Optional[StepCallback] = None,
) -> Union[SolveResult, ModalResult]:
    """Parsea un archivo YAML, arma el modelo completo y lo resuelve.

    Despacha según el tipo de solver declarado en el YAML:

    - ``ModalSolver`` (ADR 0009) → :func:`run_modal`, retorna ``ModalResult``.
    - cualquier otro → :func:`run`, retorna ``SolveResult``.

    Returns
    -------
    SolveResult | ModalResult
        Resultado agregado. Inspeccionar el tipo concreto si el consumidor
        necesita distinguir (e.g. ``isinstance(result, ModalResult)``).
    """
    # Import diferido para evitar costes si el consumidor no usa YAML.
    from fenix.utils.yaml_parser import YamlParser

    parser = YamlParser(str(path))
    domain = parser.parse()
    # Debe numerarse antes de get_external_forces, que indexa por DOF.
    domain.generate_equation_numbers()
    assembler = Assembler(domain)
    solver = parser.get_solver(assembler)

    # ModalSolver tiene firma propia `solve()` sin F_applied y devuelve
    # ModalResult: no encaja en el pipeline genérico de `run`. Se despacha
    # a `run_modal`, que también cachea el resultado en `domain.last_result`.
    if isinstance(solver, ModalSolver):
        return run_modal(domain, assembler=assembler, solver=solver)

    F_ext = parser.get_external_forces() + parser.get_body_load(assembler)
    return run(
        domain,
        assembler=assembler,
        solver=solver,
        F_applied=F_ext,
        step_callback=step_callback,
    )
