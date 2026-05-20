# solidum_fem/solidum/math/solvers/__init__.py
"""Paquete de solvers no lineales y dinámicos.

Históricamente un único archivo, dividido por solver para abaratar las
extensiones futuras (Reglas.md §1). Cada solver vive en su módulo:

- ``linear`` — :class:`LinearSolver` (un paso, sistema algebraico lineal).
- ``nonlinear`` — :class:`NonlinearSolver` (Newton-Raphson con paso adaptativo).
- ``arclength`` — :class:`ArcLengthSolver` (Crisfield cilíndrico, post-crítico).
- ``modal`` — :class:`ModalSolver` (problema de valores característicos
  ``K·φ = ω²·M·φ``; ADR 0009 fase 1).
- ``newmark`` — :class:`NewmarkSolver` (integración Newmark-β lineal para
  ``M·ü + C·u̇ + K·u = F(t)``; ADR 0009 fase 3) y :class:`NewtonNewmarkSolver`
  (variante no lineal con Newton dentro de cada paso; ADR 0009 fase 4).
  Variantes HHT-α (Hilber-Hughes-Taylor 1977) con disipación numérica
  controlada: :class:`HHTSolver` (lineal) y :class:`NewtonHHTSolver`
  (no lineal).
- ``central_difference`` — :class:`CentralDifferenceSolver` (integración
  explícita por diferencias centradas, lineal y no lineal en una sola
  clase con parámetro ``nonlinear``; ADR 0009 fase 5). Requiere masa
  lumped diagonal.
- ``harmonic`` — :class:`HarmonicSolver` (respuesta forzada armónica en
  el dominio de la frecuencia, ADR 0009 fase 6). Aritmética compleja con
  barrido sobre ``ω``. ``PIPELINE_KIND = "harmonic"`` despacha a
  ``solidum.entry.run_harmonic``.
- ``response_spectrum`` — :class:`ResponseSpectrumSolver` (combinación
  modal SRSS/CQC sobre un espectro de respuesta, ADR 0009 fase 7).
  Orquesta ModalSolver interno y delega la combinación en
  ``solidum.math.modal_response``. ``PIPELINE_KIND = "spectrum"`` despacha
  a ``solidum.entry.run_response_spectrum``.

Los reexports a nivel de paquete preservan los imports históricos
``from solidum.math.solvers import LinearSolver`` sin alterar al consumidor.
La importación del paquete dispara el ``@SolverRegistry.register`` de cada
submódulo como efecto colateral, idéntico al comportamiento del archivo
monolítico previo.

Despacho de pipeline declarativo (atributo de clase ``PIPELINE_KIND``)
elimina la cadena de ``isinstance`` en ``solidum.entry.run_yaml``: cada
solver declara su pipeline (``"static"``, ``"modal"``, ``"transient"``) y
``run_yaml`` despacha por valor literal. Nuevos solvers no clásicos no
requieren tocar ``entry.py`` (regla C de auditoría arquitectural 2026-05-13,
aplicada 2026-05-18 al añadir :class:`CentralDifferenceSolver`).
"""
from solidum.math.solvers._shared import CholeskyNotPositiveDefiniteError
from solidum.math.solvers.linear import LinearSolver
from solidum.math.solvers.nonlinear import NonlinearSolver
from solidum.math.solvers.arclength import ArcLengthSolver
from solidum.math.solvers.dissipation_arclength import DissipationArcLengthSolver
from solidum.math.solvers.modal import ModalSolver
from solidum.math.solvers.newmark import (
    HHTSolver,
    NewmarkSolver,
    NewtonHHTSolver,
    NewtonNewmarkSolver,
)
from solidum.math.solvers.central_difference import CentralDifferenceSolver
from solidum.math.solvers.harmonic import HarmonicSolver
from solidum.math.solvers.response_spectrum import ResponseSpectrumSolver

__all__ = [
    "CholeskyNotPositiveDefiniteError",
    "LinearSolver",
    "NonlinearSolver",
    "ArcLengthSolver",
    "DissipationArcLengthSolver",
    "ModalSolver",
    "NewmarkSolver",
    "NewtonNewmarkSolver",
    "HHTSolver",
    "NewtonHHTSolver",
    "CentralDifferenceSolver",
    "HarmonicSolver",
    "ResponseSpectrumSolver",
]
