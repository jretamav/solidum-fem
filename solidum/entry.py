"""Entrypoints públicos de Solidum FEM (ADR 0002).

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

from solidum.core.domain import Domain
from solidum.math.assembly import Assembler
from solidum.math.solvers import (  # noqa: F401 — NewmarkSolver default de run_transient
    HarmonicSolver,
    LinearSolver,
    ModalSolver,
    NewmarkSolver,
    ResponseSpectrumSolver,
)
from solidum.logging import get_logger
from solidum.results import (
    HarmonicResult,
    ModalResult,
    ResponseSpectrumResult,
    SolveResult,
    TransientResult,
    build_solve_result,
)

_log = get_logger("solidum.entry")

StepCallback = Callable[[int, np.ndarray, float], None]


# Valores válidos para el atributo declarativo PIPELINE_KIND de los solvers.
# Cualquier solver con un valor fuera de este conjunto será rechazado en
# ``run_yaml`` con ValueError (regla C del ADR 0009 + fail-fast del proyecto).
_KNOWN_PIPELINE_KINDS = frozenset({
    "static", "modal", "transient", "harmonic", "spectrum",
})


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


def run_harmonic(
    domain: Domain,
    *,
    assembler: Assembler | None = None,
    solver: HarmonicSolver | None = None,
    **solver_kwargs,
) -> HarmonicResult:
    """Ejecuta un análisis de respuesta armónica y retorna el resultado.

    Acepta cualquier solver con ``PIPELINE_KIND="harmonic"`` y firma
    ``solve() -> HarmonicResult``. Si ``solver is None`` construye un
    :class:`HarmonicSolver` con los kwargs siguientes.

    Parameters
    ----------
    domain
        Dominio con nodos, elementos y BCs ya configurados.
    assembler
        Assembler opcional. Si es ``None`` se construye sobre ``domain``.
    solver
        Instancia ya vinculada al assembler. Si es ``None``, se construye
        con los kwargs siguientes.
    **solver_kwargs
        Parámetros del constructor de :class:`HarmonicSolver`
        (``omega`` / ``omega_min`` / ``omega_max`` / ``n_omega`` /
        ``scale`` / ``F_amplitude`` / ``rayleigh`` / ``lumping``).

    Returns
    -------
    HarmonicResult
        Amplitud compleja del desplazamiento para cada frecuencia del
        barrido. Queda también en ``domain.last_result``.
    """
    if domain.total_dofs == 0:
        domain.generate_equation_numbers()
    if assembler is None:
        assembler = Assembler(domain)
    if solver is None:
        solver = HarmonicSolver(assembler, **solver_kwargs)

    result = solver.solve()
    domain.last_result = result
    return result


def run_response_spectrum(
    domain: Domain,
    *,
    assembler: Assembler | None = None,
    solver: ResponseSpectrumSolver | None = None,
    **solver_kwargs,
) -> ResponseSpectrumResult:
    """Ejecuta un análisis de respuesta espectral y retorna el resultado.

    Acepta cualquier solver con ``PIPELINE_KIND="spectrum"`` y firma
    ``solve() -> ResponseSpectrumResult``. Si ``solver is None``
    construye un :class:`ResponseSpectrumSolver` con los kwargs
    siguientes.

    Parameters
    ----------
    domain
        Dominio con nodos, elementos y BCs ya configurados.
    assembler
        Assembler opcional. Si es ``None`` se construye sobre ``domain``.
    solver
        Instancia ya vinculada al assembler. Si es ``None``, se construye
        con los kwargs siguientes.
    **solver_kwargs
        Parámetros del constructor de :class:`ResponseSpectrumSolver`
        (``n_modes``, ``direction``, ``spectrum``, ``combination``,
        ``damping``, ``sigma``, ``lumping``).

    Returns
    -------
    ResponseSpectrumResult
        Respuesta máxima envolvente. Queda también en
        ``domain.last_result``.
    """
    if domain.total_dofs == 0:
        domain.generate_equation_numbers()
    if assembler is None:
        assembler = Assembler(domain)
    if solver is None:
        solver = ResponseSpectrumSolver(assembler, **solver_kwargs)

    result = solver.solve()
    domain.last_result = result
    return result


def run_yaml(
    path: str | Path,
    *,
    step_callback: StepCallback | None = None,
) -> SolveResult | ModalResult | TransientResult | HarmonicResult | ResponseSpectrumResult:
    """Parsea un archivo YAML, arma el modelo completo y lo resuelve.

    Despacha según el atributo de clase ``PIPELINE_KIND`` del solver
    declarado en el YAML:

    - ``"modal"`` → :func:`run_modal`, retorna ``ModalResult``.
    - ``"transient"`` → :func:`run_transient`, retorna ``TransientResult``.
    - ``"harmonic"`` → :func:`run_harmonic`, retorna ``HarmonicResult``.
    - ``"spectrum"`` → :func:`run_response_spectrum`, retorna ``ResponseSpectrumResult``.
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
    from solidum.utils.yaml_parser import YamlParser

    parser = YamlParser(str(path))
    domain = parser.parse()
    # Debe numerarse antes de get_external_forces, que indexa por DOF.
    domain.generate_equation_numbers()
    assembler = Assembler(domain)
    solver = parser.get_solver(assembler)

    pipeline_kind = getattr(solver, "PIPELINE_KIND", "static")
    if pipeline_kind not in _KNOWN_PIPELINE_KINDS:
        raise ValueError(
            f"run_yaml: PIPELINE_KIND={pipeline_kind!r} no reconocido en "
            f"{type(solver).__name__}. Valores válidos: "
            f"{sorted(_KNOWN_PIPELINE_KINDS)}."
        )

    # Aviso del usuario: ``step_callback`` sólo lo consume el pipeline estático
    # (``run`` → ``NonlinearSolver.solve(..., step_callback=...)`` o
    # ``ArcLengthSolver.solve(...)``). Pasarlo a ``run_yaml`` con un solver
    # modal/transitorio/armónico/espectral es no-op silencioso, lo cual
    # confunde al usuario que espera ver progreso. Aviso explícito (H-1.8).
    if step_callback is not None and pipeline_kind != "static":
        _log.warning(
            "run_yaml: step_callback ignorado — solver %s declara "
            "PIPELINE_KIND=%r, no 'static'. step_callback sólo afecta "
            "al pipeline estático (NonlinearSolver/ArcLengthSolver).",
            type(solver).__name__, pipeline_kind,
        )

    if pipeline_kind == "modal":
        return run_modal(domain, assembler=assembler, solver=solver)
    if pipeline_kind == "transient":
        return run_transient(domain, assembler=assembler, solver=solver)
    if pipeline_kind == "harmonic":
        return run_harmonic(domain, assembler=assembler, solver=solver)
    if pipeline_kind == "spectrum":
        return run_response_spectrum(
            domain, assembler=assembler, solver=solver,
        )

    # "static" — pipeline genérico con fuerzas externas y peso propio.
    F_ext = parser.get_external_forces() + parser.get_body_load(assembler)
    return run(
        domain,
        assembler=assembler,
        solver=solver,
        F_applied=F_ext,
        step_callback=step_callback,
    )
