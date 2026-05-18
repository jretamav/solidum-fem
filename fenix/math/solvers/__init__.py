# fenix_fem/fenix/math/solvers/__init__.py
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

Los reexports a nivel de paquete preservan los imports históricos
``from fenix.math.solvers import LinearSolver`` sin alterar al consumidor.
La importación del paquete dispara el ``@SolverRegistry.register`` de cada
submódulo como efecto colateral, idéntico al comportamiento del archivo
monolítico previo.
"""
from fenix.math.solvers._shared import CholeskyNotPositiveDefiniteError
from fenix.math.solvers.linear import LinearSolver
from fenix.math.solvers.nonlinear import NonlinearSolver
from fenix.math.solvers.arclength import ArcLengthSolver
from fenix.math.solvers.modal import ModalSolver
from fenix.math.solvers.newmark import (
    HHTSolver,
    NewmarkSolver,
    NewtonHHTSolver,
    NewtonNewmarkSolver,
)

__all__ = [
    "CholeskyNotPositiveDefiniteError",
    "LinearSolver",
    "NonlinearSolver",
    "ArcLengthSolver",
    "ModalSolver",
    "NewmarkSolver",
    "NewtonNewmarkSolver",
    "HHTSolver",
    "NewtonHHTSolver",
]
