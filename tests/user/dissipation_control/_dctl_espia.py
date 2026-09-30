"""Espía de los problemas de Newton de ``DissipationArcLengthSolver``.

Compartido por ``test_dctl_dissipation_arclength.py`` y ``test_dctl_retroceso.py``.

El atributo ``_mode`` del solver, que ven los callbacks, no dice con qué
restricción se resolvió cada paso: hasta el 2026-09-30 el modo pasaba a
``"dissipation"`` tras el commit y el predictor siguiente lo devolvía a
``"cylindrical"`` antes de construir el problema de Newton, así que los tests
veían el modo de disipación sin que la restricción de disipación llegara a
usarse. El espía mira el problema que de verdad se resuelve: envuelve
``solver.corrector.run`` (el corrector compartido, atributo público de
``ArcLengthSolver``) y anota el ``mode`` del problema (atributo público de
``solidum.math.solvers.ArcProblem``) y si el paso convergió.
"""
from __future__ import annotations


def espiar_problemas_de_newton(solver) -> list[tuple[str, bool]]:
    """Instala el espía en ``solver`` y devuelve la lista que irá llenando.

    Cada intento de paso añade ``(modo, convergió)``, con ``modo`` en
    ``{"cylindrical", "dissipation", "newton"}`` (``"newton"``: llegada exacta
    a ``max_lambda``). El callback de paso se llama después del ``run`` que
    convergió, así que el último elemento de la lista es el del paso que el
    callback está viendo.
    """
    intentos: list[tuple[str, bool]] = []
    run = solver.corrector.run

    def run_espiado(problem, x0, **kwargs):
        res = run(problem, x0, **kwargs)
        intentos.append((problem.mode, bool(res.converged)))
        return res

    solver.corrector.run = run_espiado
    return intentos


def convergidos(intentos, modo: str) -> int:
    """Número de problemas de Newton en ``modo`` que convergieron."""
    return sum(1 for m, ok in intentos if ok and m == modo)
