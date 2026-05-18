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
from typing import Any, Callable

import numpy as np

from fenix.core.domain import Domain
from fenix.math.assembly import Assembler
from fenix.math.solvers import LinearSolver, ModalSolver, NewmarkSolver  # noqa: F401 — NewmarkSolver es el default de run_transient
from fenix.results import ModalResult, SolveResult, TransientResult, build_solve_result


StepCallback = Callable[[int, np.ndarray, float], None]


def _invoke_solve(solver: Any, F_applied: np.ndarray,
                   step_callback: StepCallback | None) -> np.ndarray:
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
    assembler: Assembler | None = None,
    solver: Any | None = None,
    F_applied: np.ndarray | None = None,
    step_callback: StepCallback | None = None,
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
    assembler: Assembler | None = None,
    solver: ModalSolver | None = None,
    n_modes: int | None = None,
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


def run_transient(
    domain: Domain,
    *,
    assembler: Assembler | None = None,
    solver: Any | None = None,
    t_end: float | None = None,
    dt: float | None = None,
    **solver_kwargs,
) -> TransientResult:
    """Ejecuta un análisis dinámico transitorio y retorna el resultado.

    Acepta cualquier solver con ``PIPELINE_KIND="transient"`` y firma
    ``solve() -> TransientResult`` (NewmarkSolver, NewtonNewmarkSolver,
    HHTSolver, NewtonHHTSolver, CentralDifferenceSolver). Si ``solver is
    None`` construye un ``NewmarkSolver`` con los kwargs siguientes.

    Parameters
    ----------
    domain
        Dominio con nodos, elementos y BCs ya configurados.
    assembler
        Assembler opcional. Si es ``None`` se construye sobre ``domain``.
    solver
        Instancia de solver transitorio ya vinculada al assembler.
    t_end, dt
        Tiempo final y paso temporal. Requeridos cuando ``solver is None``.
    **solver_kwargs
        Resto de parámetros del constructor del NewmarkSolver default.

    Returns
    -------
    TransientResult
        Historiales temporales de ``u, u̇, ü``. Queda también en
        ``domain.last_result``.
    """
    if domain.total_dofs == 0:
        domain.generate_equation_numbers()
    if assembler is None:
        assembler = Assembler(domain)
    if solver is None:
        if t_end is None or dt is None:
            raise ValueError(
                "run_transient: pasar `solver` o (`t_end` y `dt`)."
            )
        solver = NewmarkSolver(assembler, t_end=t_end, dt=dt, **solver_kwargs)

    result = solver.solve()
    domain.last_result = result
    return result


def run_yaml(
    path: str | Path,
    *,
    step_callback: StepCallback | None = None,
) -> SolveResult | ModalResult | TransientResult:
    """Parsea un archivo YAML, arma el modelo completo y lo resuelve.

    Despacha según el atributo de clase ``PIPELINE_KIND`` del solver
    declarado en el YAML:

    - ``"modal"`` → :func:`run_modal`, retorna ``ModalResult``.
    - ``"transient"`` → :func:`run_transient`, retorna ``TransientResult``.
    - ``"static"`` (default) → :func:`run`, retorna ``SolveResult``.

    El despacho declarativo (atributo de clase en lugar de cadena de
    ``isinstance``) permite añadir nuevos solvers no clásicos sin tocar
    este archivo: el solver declara su ``PIPELINE_KIND`` y queda
    automáticamente despachado al pipeline correcto.

    Returns
    -------
    SolveResult | ModalResult | TransientResult
        Resultado agregado. Inspeccionar el tipo concreto si el consumidor
        necesita distinguir (e.g. ``isinstance(result, TransientResult)``).
    """
    # Import diferido para evitar costes si el consumidor no usa YAML.
    from fenix.utils.yaml_parser import YamlParser

    parser = YamlParser(str(path))
    domain = parser.parse()
    # Debe numerarse antes de get_external_forces, que indexa por DOF.
    domain.generate_equation_numbers()
    assembler = Assembler(domain)
    solver = parser.get_solver(assembler)

    pipeline_kind = getattr(solver, "PIPELINE_KIND", "static")

    if pipeline_kind == "modal":
        return run_modal(domain, assembler=assembler, solver=solver)
    if pipeline_kind == "transient":
        return run_transient(domain, assembler=assembler, solver=solver)

    # "static" o desconocido — cae al pipeline genérico de fuerzas externas.
    F_ext = parser.get_external_forces() + parser.get_body_load(assembler)
    return run(
        domain,
        assembler=assembler,
        solver=solver,
        F_applied=F_ext,
        step_callback=step_callback,
    )
