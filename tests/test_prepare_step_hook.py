"""Gancho de inicio de paso ``Element.prepare_step`` (ADR 0020, P5).

Es parte del contrato genérico de un elemento, como la tarea de inicio de
paso de un elemento de usuario de FEAP: los solvers que avanzan por pasos lo
invocan en cada elemento, una vez por intento de paso, con el campo
**convergido** del paso anterior. Un elemento que toma decisiones discretas
entre pasos (activar una grieta, un contacto) lo sobreescribe; en el resto es
no-op.

La prueba usa un elemento de juguete sin registrar (un ``Truss2D`` que anota
cada llamada), para no depender de ninguna formulación concreta ni añadir
clases a los barridos de contrato.
"""
import numpy as np
import pytest

from solidum.core.domain import Domain
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.solvers import (
    ArcLengthSolver,
    DissipationArcLengthSolver,
    IndirectDisplacementSolver,
    LinearSolver,
    NonlinearSolver,
)
from solidum.math.solvers.newmark import NewtonHHTSolver, NewtonNewmarkSolver

E, A, L, P = 200.0e9, 1.0e-3, 1.0, 1.0e3


class _BarraEspia(Truss2D):
    """Anota el ``U_committed`` que recibe en cada inicio de paso."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.llamadas = []

    def prepare_step(self, U_committed):
        self.llamadas.append(np.array(U_committed, copy=True))


def _barra():
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    barra = _BarraEspia(1, [n1, n2], Elastic1D(E=E, density=7850.0), A=A)
    dom.add_element(barra)
    n1.fix_dof("ux", 0.0)
    n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = P
    return dom, barra, n2, F


def _estaticos():
    return {
        "NonlinearSolver": lambda asm, n2: NonlinearSolver(asm, num_steps=3),
        "ArcLengthSolver": lambda asm, n2: ArcLengthSolver(asm, initial_dlambda=0.4),
        "IndirectDisplacementSolver": lambda asm, n2: IndirectDisplacementSolver(
            asm, initial_dlambda=0.4, control=[(n2, "ux", 1.0)]),
        "DissipationArcLengthSolver": lambda asm, n2: DissipationArcLengthSolver(
            asm, initial_dlambda=0.4, initial_tau=1.0e-3),
    }


@pytest.mark.parametrize("nombre", list(_estaticos()))
def test_los_solvers_estaticos_por_pasos_llaman_al_gancho(nombre):
    dom, barra, n2, F = _barra()
    _estaticos()[nombre](Assembler(dom), n2).solve(F)
    assert len(barra.llamadas) >= 2, f"{nombre}: {len(barra.llamadas)} llamadas"
    # El primer paso parte del estado inicial; los siguientes, del convergido.
    np.testing.assert_array_equal(barra.llamadas[0], 0.0)
    assert barra.llamadas[1][n2.dofs["ux"]] > 0.0


@pytest.mark.parametrize("solver_cls", [NewtonNewmarkSolver, NewtonHHTSolver])
def test_los_solvers_dinamicos_no_lineales_llaman_al_gancho(solver_cls):
    dom, barra, n2, F = _barra()
    dt, pasos = 1.0e-5, 4
    solver_cls(Assembler(dom), pasos * dt, dt, F_func=lambda t: F).solve()
    assert len(barra.llamadas) == pasos
    np.testing.assert_array_equal(barra.llamadas[0], 0.0)


def test_el_estatico_lineal_no_llama_al_gancho():
    """Un análisis lineal no avanza por pasos: no hay decisiones entre pasos."""
    dom, barra, n2, F = _barra()
    LinearSolver(Assembler(dom)).solve(F)
    assert barra.llamadas == []


def test_la_base_es_no_op():
    dom, barra, n2, F = _barra()
    Truss2D.prepare_step(barra, np.zeros(dom.total_dofs))      # no lanza
