"""Gancho de revisión del paso convergido ``Element.revise_step`` (ADR 0041, P9).

Al converger un intento de paso y antes de confirmarlo, el corrector
compartido pregunta a los elementos si el campo convergido cambia una
decisión discreta suya (activar una grieta, un contacto). Si alguno cambia,
el intento no se confirma y se repite desde su iterado inicial, en cascada,
hasta que ninguno cambia. Es la comprobación de fin de paso de un elemento de
usuario de FEAP.

Las pruebas usan barras de juguete sin registrar: una ``Truss2D`` que se
«rompe» (su módulo cae a una fracción, para siempre) cuando, en una revisión,
su deformación supera un umbral. Así se comprueba lo que el gancho promete sin
depender de ninguna formulación:

- el paso que cruza el umbral se confirma ya roto (sin sobrepaso de un paso);
- entre la repetición y la original no se confirma nada;
- la cascada (al romperse una barra se carga la otra, que se rompe en el mismo
  paso) se resuelve en un solo paso;
- los seis solvers que avanzan por pasos con el corrector compartido la hacen;
- un elemento que viola el contrato (cambia siempre) agota el tope y el
  intento fracasa con un motivo que lo dice;
- sin elementos que redefinan el gancho, nada cambia.
"""
import numpy as np
import pytest

from solidum.core.domain import Domain
from solidum.core.element import Element
from solidum.elements.truss import Truss2D
from solidum.materials.elastic import Elastic1D
from solidum.math.assembly import Assembler
from solidum.math.solvers import (
    ArcLengthSolver,
    IndirectDisplacementSolver,
    NonlinearSolver,
)
from solidum.constants import MAX_STEP_REVISIONS
from solidum.math.solvers.newmark import NewtonHHTSolver, NewtonNewmarkSolver

E, A, L = 200.0e9, 1.0e-3, 1.0
EPS_ROTURA = 1.0e-3          # deformación a la que se rompe la barra débil
FRACCION = 0.25              # módulo que le queda a una barra rota


class _BarraFragil(Truss2D):
    """Barra que, al revisarse con ``ε > eps_rotura``, conserva sólo
    ``FRACCION`` de su rigidez (irreversible). Anota las revisiones, las
    roturas y los commits."""

    def __init__(self, *args, eps_rotura=EPS_ROTURA, **kwargs):
        super().__init__(*args, **kwargs)
        self.eps_rotura = eps_rotura
        self.rota = False
        self.revisiones = 0
        self.commits = 0
        self.rotura_en_commit = None    # nº de commits confirmados antes de romperse

    def _strain(self, U):
        u = self.get_local_displacements(U)
        return (u[2] - u[0]) / L

    def compute_element_state(self, u_e):
        K_e, F_int_e = super().compute_element_state(u_e)
        if self.rota:
            return FRACCION * K_e, FRACCION * F_int_e
        return K_e, F_int_e

    def revise_step(self, U_trial):
        self.revisiones += 1
        if not self.rota and self._strain(U_trial) > self.eps_rotura:
            self.rota = True
            self.rotura_en_commit = self.commits
            return True
        return False

    def commit_state(self):
        super().commit_state()
        self.commits += 1


class _BarraQueSiempreCambia(Truss2D):
    """Viola el contrato: dice haber cambiado en cada revisión."""

    def revise_step(self, U_trial):
        return True


def _barra(cls=_BarraFragil, n_barras=1, **kw):
    """``n_barras`` barras en paralelo entre el nodo 1 (fijo) y el 2 (libre
    en x). La primera es la débil (umbral ``EPS_ROTURA``); las demás rompen
    con 1.5 veces esa deformación."""
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    barras = []
    for i in range(n_barras):
        extra = {} if cls is not _BarraFragil else {
            "eps_rotura": EPS_ROTURA * (1.0 if i == 0 else 1.5)}
        b = cls(i + 1, [n1, n2], Elastic1D(E=E, density=7850.0), A=A, **extra, **kw)
        dom.add_element(b)
        barras.append(b)
    n1.fix_dof("ux", 0.0)
    n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    return dom, barras, n2


def _desplazamiento(dom, n2, delta):
    n2.fix_dof("ux", delta)
    dom.generate_equation_numbers(verbose=False)
    return np.zeros(dom.total_dofs)


# ----------------------------------------------------------------------
# Lo esencial: el paso que cruza el umbral se confirma ya roto
# ----------------------------------------------------------------------

def test_el_paso_que_cruza_el_umbral_se_confirma_ya_roto():
    """Desplazamiento impuesto en 10 pasos hasta 2.05·ε_rotura·L. El paso 5
    (ε = 1.025·ε_rotura) cruza el umbral: se confirma con la barra rota, con
    la fuerza del módulo reducido. Con la decisión al inicio del paso
    siguiente se confirmaría intacta, con la fuerza 1.025·E·A·ε_rotura."""
    dom, (barra,), n2 = _barra()
    delta = 2.05 * EPS_ROTURA * L
    F = _desplazamiento(dom, n2, delta)
    asm = Assembler(dom)
    fuerzas = []

    def cb(step, U, lam):
        _, F_int = asm.assemble_non_linear_system(U)
        fuerzas.append((lam * delta / L, float(F_int[n2.dofs["ux"]]), barra.rota))

    NonlinearSolver(asm, num_steps=10, adaptive=False).solve(F, step_callback=cb)
    eps, fuerza, rota = map(np.array, zip(*fuerzas))
    k = int(np.argmax(eps > EPS_ROTURA))                # primer paso por encima
    assert rota[k] and not rota[k - 1]
    assert barra.rotura_en_commit == k                  # rota antes de confirmar el paso k
    np.testing.assert_allclose(fuerza[k], FRACCION * E * A * eps[k], rtol=1e-12)
    # Ningún estado confirmado intacto por encima del umbral.
    assert np.all(eps[~rota] <= EPS_ROTURA)


def test_la_repeticion_no_confirma_el_intento_descartado():
    dom, (barra,), n2 = _barra()
    F = _desplazamiento(dom, n2, 2.05 * EPS_ROTURA * L)
    solver = NonlinearSolver(Assembler(dom), num_steps=10, adaptive=False)
    solver.solve(F)
    # Un commit por paso: el intento descartado no se confirmó.
    assert barra.commits == 10
    # Una revisión por paso más la del intento que se repitió.
    assert barra.revisiones == 11


def test_cascada_en_un_solo_paso():
    """Dos barras en paralelo con carga de fuerza. Al romperse la débil, la
    otra toma su parte y supera su propio umbral en el mismo paso: las dos
    roturas se resuelven antes de confirmar ese paso."""
    dom, barras, n2 = _barra(n_barras=2)
    # Fuerza que lleva a la débil justo por encima de su umbral en el paso 4 de 5.
    P_rotura = 2.0 * E * A * EPS_ROTURA
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = 1.2 * P_rotura
    asm = Assembler(dom)
    estados = []

    def cb(step, U, lam):
        estados.append(tuple(b.rota for b in barras))

    solver = NonlinearSolver(asm, num_steps=5, adaptive=False)
    solver.solve(F, step_callback=cb)
    k = estados.index((True, True))
    assert estados[k - 1] == (False, False)             # las dos en el mismo paso
    assert barras[0].rotura_en_commit == barras[1].rotura_en_commit == k
    assert solver.n_revisions == 2                      # dos rondas que cambiaron
    assert barras[0].commits == 5                       # un commit por paso


# ----------------------------------------------------------------------
# Los solvers que avanzan por pasos con el corrector compartido
# ----------------------------------------------------------------------

def _estaticos():
    return {
        "NonlinearSolver": lambda asm, n2: NonlinearSolver(asm, num_steps=10, adaptive=False),
        "ArcLengthSolver": lambda asm, n2: ArcLengthSolver(asm, initial_dlambda=0.1, max_steps=200),
        "IndirectDisplacementSolver": lambda asm, n2: IndirectDisplacementSolver(
            asm, initial_dlambda=0.1, max_steps=200, control=[(n2, "ux", 1.0)]),
    }


@pytest.mark.parametrize("nombre", list(_estaticos()))
def test_los_solvers_estaticos_revisan_antes_de_confirmar(nombre):
    """Carga de fuerza hasta 1.5 veces la de rotura de la barra débil (con
    dos barras en paralelo). Con cualquier solver estático, todo estado
    confirmado con la débil intacta está por debajo de su umbral."""
    dom, barras, n2 = _barra(n_barras=2)
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = 0.75 * 2.0 * E * A * EPS_ROTURA * 2.0
    asm = Assembler(dom)
    confirmados = []

    def cb(step, U, lam):
        confirmados.append((barras[0]._strain(U), barras[0].rota))

    _estaticos()[nombre](asm, n2).solve(F, step_callback=cb)
    assert barras[0].rota, nombre
    assert all(eps <= EPS_ROTURA * (1 + 1e-12) for eps, rota in confirmados if not rota), nombre
    assert barras[0].revisiones > barras[0].commits, nombre     # hubo una repetición


@pytest.mark.parametrize("solver_cls", [NewtonNewmarkSolver, NewtonHHTSolver])
def test_los_solvers_dinamicos_revisan_antes_de_confirmar(solver_cls):
    # Período de la barra ≈ 0.9 ms: en 0.5 ms la carga súbita lleva el
    # desplazamiento por encima del de rotura (ω·t ≈ 3.5).
    dom, (barra,), n2 = _barra()
    dt, pasos = 5.0e-5, 10
    P = 3.0 * E * A * EPS_ROTURA
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = P
    solver_cls(Assembler(dom), pasos * dt, dt, F_func=lambda t: F).solve()
    assert barra.revisiones > barra.commits
    assert barra.commits == pasos


# ----------------------------------------------------------------------
# Contrato
# ----------------------------------------------------------------------

def test_un_elemento_que_siempre_cambia_agota_el_tope():
    """El tope protege de un elemento que viola el contrato (cambios no
    finitos): el intento fracasa; sin bisección, el análisis se detiene con
    un mensaje que nombra la revisión."""
    dom, _, n2 = _barra(cls=_BarraQueSiempreCambia)
    F = _desplazamiento(dom, n2, 1.0e-4)
    with pytest.raises(Exception) as exc:
        NonlinearSolver(Assembler(dom), num_steps=2, adaptive=False).solve(F)
    assert str(MAX_STEP_REVISIONS) in str(exc.value)
    assert "revis" in str(exc.value)


def test_sin_elementos_que_revisen_no_hay_revisiones():
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    dom.add_element(Truss2D(1, [n1, n2], Elastic1D(E=E, density=7850.0), A=A))
    n1.fix_dof("ux", 0.0)
    n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    assert asm.revise_all_steps(np.zeros(dom.total_dofs)) is False
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = 1.0e3
    solver = NonlinearSolver(asm, num_steps=3)
    U = solver.solve(F)
    assert U[n2.dofs["ux"]] == pytest.approx(1.0e3 * L / (E * A), rel=1e-12)


def test_la_base_es_no_op():
    dom, (barra,), n2 = _barra()
    assert Element.revise_step(barra, np.zeros(dom.total_dofs)) is False


# ----------------------------------------------------------------------
# Pruebas de la revisión de arquitectura (informes/revision_arquitectura_
# 2026-11-18_adr0041.md)
# ----------------------------------------------------------------------

class _BarraQueAnotaRevisiones(_BarraFragil):
    """Anota el desplazamiento del extremo con que se revisa cada vez."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.u_revisados = []

    def revise_step(self, U_converged):
        self.u_revisados.append(float(self.get_local_displacements(U_converged)[2]))
        return super().revise_step(U_converged)


@pytest.mark.parametrize("objetivo", ["max_lambda", "control_target"])
def test_el_arco_no_revisa_un_punto_que_rebasa_el_final(objetivo):
    """Un paso del arco que rebasa ``max_lambda`` o ``control_target`` no se
    confirma (llegada exacta): tampoco se revisa. Ninguna revisión ve un
    desplazamiento mayor que el final, y la llegada sí se revisa."""
    dom, (barra,), n2 = _barra(cls=_BarraQueAnotaRevisiones, eps_rotura=10.0)  # no se rompe
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = 1.0e3
    asm = Assembler(dom)
    u_final_exacto = 1.0e3 * L / (E * A)
    if objetivo == "max_lambda":
        solver = ArcLengthSolver(asm, initial_dlambda=0.3, max_lambda=1.0, max_steps=50)
    else:
        solver = ArcLengthSolver(asm, initial_dlambda=0.3, max_lambda=10.0, max_steps=50,
                                 control=[(n2, "ux", 1.0)], control_target=u_final_exacto)
    U = solver.solve(F)
    u_final = U[n2.dofs["ux"]]
    assert u_final == pytest.approx(u_final_exacto, rel=1e-9)
    assert max(barra.u_revisados) <= u_final * (1.0 + 1e-12)
    assert barra.u_revisados[-1] == pytest.approx(u_final, rel=1e-12)   # la llegada se revisa
    assert barra.revisiones == barra.commits                           # una por punto confirmado


class _BarraRotaDeAntemano(_BarraFragil):
    """Se rompe al inicio del paso ``paso_rotura`` (en ``prepare_step``), sin
    revisión: la referencia de «la decisión tomada de antemano»."""

    def __init__(self, *args, paso_rotura, **kwargs):
        super().__init__(*args, **kwargs)
        self.paso_rotura = paso_rotura
        self.pasos = 0

    def prepare_step(self, U_committed):
        self.pasos += 1
        if self.pasos == self.paso_rotura:
            self.rota = True

    def revise_step(self, U_converged):
        self.revisiones += 1
        return False


def _paralelo_con_fuerza(cls_debil, **kw):
    """Barra débil (que se rompe) en paralelo con una elástica, con fuerza
    impuesta en el nodo libre (un GDL libre): con 10 pasos, la débil cruza su
    umbral en el paso 5 (ε = 1.025·ε_rotura)."""
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    debil = cls_debil(1, [n1, n2], Elastic1D(E=E, density=7850.0), A=A, **kw)
    dom.add_element(debil)
    dom.add_element(Truss2D(2, [n1, n2], Elastic1D(E=E, density=7850.0), A=A))
    n1.fix_dof("ux", 0.0)
    n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = 4.1 * E * A * EPS_ROTURA
    return dom, debil, n2, F


@pytest.mark.parametrize("opciones", [
    {}, {"freeze_tangent_after_iter": 0}, {"line_search": True}],
    ids=["newton", "newton_modificado", "line_search"])
def test_el_paso_repetido_es_el_de_la_decision_tomada_de_antemano(opciones):
    """Reiniciar desde el iterado inicial del intento, sin confirmar, da bit a
    bit el mismo paso que si la decisión se hubiera tomado antes de empezarlo
    (backend LU, control de fuerza con un GDL libre): la ronda descartada no
    deja rastro (factorización congelada, historiales, estado de prueba)."""
    resultados = []
    for cls, kw in ((_BarraFragil, {}), (_BarraRotaDeAntemano, {"paso_rotura": 5})):
        dom, debil, n2, F = _paralelo_con_fuerza(cls, **kw)
        asm = Assembler(dom)
        hist = []
        NonlinearSolver(asm, num_steps=10, adaptive=False, linear_algebra="lu", **opciones).solve(
            F, step_callback=lambda s, U, lam: hist.append(U.copy()))
        resultados.append((hist, debil.rotura_en_commit if cls is _BarraFragil else None))
    (h1, k), (h2, _) = resultados
    assert k == 4                                   # se rompió al revisar el paso 5
    assert len(h1) == len(h2) == 10
    for U1, U2 in zip(h1, h2):
        np.testing.assert_array_equal(U1, U2)


def test_tras_una_revision_el_intento_se_repite_desde_su_iterado_inicial():
    """Espía estructural del reinicio: el primer ensamblaje después de una
    revisión que cambió algo es en el iterado inicial del intento (el campo
    confirmado del paso anterior en el Newton incremental), bit a bit, y no en
    el convergido descartado."""
    dom, debil, n2, F = _paralelo_con_fuerza(_BarraFragil)
    asm = Assembler(dom)
    eventos = []
    ensamblar, confirmar, revisar = asm.assemble_non_linear_system, asm.commit_all_states, debil.revise_step

    def espia_ensamblar(U):
        eventos.append(("A", np.array(U, copy=True)))
        return ensamblar(U)

    def espia_confirmar():
        eventos.append(("C", None))
        return confirmar()

    def espia_revisar(U):
        cambio = revisar(U)
        eventos.append(("R", cambio))
        return cambio

    asm.assemble_non_linear_system = espia_ensamblar
    asm.commit_all_states = espia_confirmar
    debil.revise_step = espia_revisar
    NonlinearSolver(asm, num_steps=10, adaptive=False).solve(F)
    x0 = np.zeros(dom.total_dofs)       # iterado inicial del paso en curso
    ultimo = None
    reinicios = 0
    esperar_x0 = False
    for tipo, dato in eventos:
        if tipo == "A":
            if esperar_x0:
                np.testing.assert_array_equal(dato, x0)
                reinicios += 1
                esperar_x0 = False
            ultimo = dato
        elif tipo == "R" and dato:
            assert not np.array_equal(ultimo, x0)       # el convergido no es x0
            esperar_x0 = True
        elif tipo == "C":
            x0 = ultimo                                  # el confirmado: x0 del siguiente
    assert reinicios == 1


def test_sin_elementos_que_revisen_el_resultado_no_cambia():
    """Bit a bit y con el mismo número de ensamblajes que un corrector sin
    revisión (el problema sin el método ``revise``)."""
    from solidum.math.solvers import nonlinear

    def corrida():
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [L, 0.0])
        dom.add_element(Truss2D(1, [n1, n2], Elastic1D(E=E, density=7850.0), A=A))
        n1.fix_dof("ux", 0.0)
        n1.fix_dof("uy", 0.0)
        n2.fix_dof("uy", 0.0)
        dom.generate_equation_numbers(verbose=False)
        asm = Assembler(dom)
        cuenta = [0]
        original = asm.assemble_non_linear_system

        def contado(U):
            cuenta[0] += 1
            return original(U)
        asm.assemble_non_linear_system = contado
        F = np.zeros(dom.total_dofs)
        F[n2.dofs["ux"]] = 1.0e3
        return NonlinearSolver(asm, num_steps=4).solve(F), cuenta[0]

    U_con, n_con = corrida()
    revise = nonlinear._IncrementalProblem.revise
    try:
        del nonlinear._IncrementalProblem.revise
        U_sin, n_sin = corrida()
    finally:
        nonlinear._IncrementalProblem.revise = revise
    np.testing.assert_array_equal(U_con, U_sin)
    assert n_con == n_sin


def test_el_tope_lanza_un_error_tipado():
    from solidum.math.solvers.diagnostics import StepRevisionError
    dom, _, n2 = _barra(cls=_BarraQueSiempreCambia)
    F = _desplazamiento(dom, n2, 1.0e-4)
    with pytest.raises(StepRevisionError, match="revise_step"):
        NonlinearSolver(Assembler(dom), num_steps=2, adaptive=False).solve(F)


def test_el_tope_con_biseccion_sigue_siendo_tipado():
    """Con paso adaptativo, el intento que agota el tope se biseca como
    cualquier fallo; al llegar al paso mínimo, el error sigue siendo el de la
    revisión (no un modo de divergencia del Newton, que sí convergía)."""
    from solidum.math.solvers.diagnostics import StepRevisionError
    dom, _, n2 = _barra(cls=_BarraQueSiempreCambia)
    F = _desplazamiento(dom, n2, 1.0e-4)
    with pytest.raises(StepRevisionError) as exc:
        NonlinearSolver(Assembler(dom), num_steps=2, adaptive=True, min_delta_lambda=1e-2).solve(F)
    assert exc.value.n_bisections > 0


class _BarraQueRevisaYFallaUnIntento(Truss2D):
    """Revisa (sin cambiar nunca) y hace fracasar el primer intento del
    paso 2 con una fuerza espuria que cambia de signo en cada ensamblaje."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pasos = 0
        self.signo = 1.0

    def prepare_step(self, U_committed):
        self.pasos += 1

    def revise_step(self, U_converged):
        return False

    def compute_element_state(self, u_e):
        K_e, F_int_e = super().compute_element_state(u_e)
        if self.pasos == 2:
            self.signo = -self.signo
            F_int_e = F_int_e + self.signo * 1.0e3 * np.array([0.0, 0.0, 1.0, 0.0])
        return K_e, F_int_e


def test_un_trazado_que_acaba_tras_un_intento_fallido_no_deja_el_intento_abierto():
    """Hallazgo 1 de la revisión del diff: el arco agota ``max_steps`` justo
    después de un intento fallido. El commit de posproceso de ``solidum.run``
    (ensamblar en U y confirmar) no debe tomarlo por un paso sin revisión."""
    dom, (barra,), n2 = _barra(cls=_BarraQueRevisaYFallaUnIntento)
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = 1.0e3
    asm = Assembler(dom)
    solver = ArcLengthSolver(asm, initial_dlambda=0.4, max_steps=2)
    U = solver.solve(F)
    assert solver.n_bisections == 1                     # el paso 2 fracasó y se redujo
    asm.assemble_non_linear_system(U)                   # lo que hace solidum.run al final
    asm.commit_all_states()                             # no lanza


def test_un_solver_que_confirma_sin_revisar_lo_rechaza_el_ensamblador():
    """Guarda del commit: con elementos que revisan, confirmar un intento de
    paso cuya evaluación no se revisó es un error (el solver no hace la
    revisión que piden; sus decisiones llegarían tarde o nunca)."""
    from solidum.math.solvers.central_difference import CentralDifferenceSolver
    dom, (barra,), n2 = _barra()
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = 1.0e3
    with pytest.raises(ValueError, match="ADR 0041"):
        CentralDifferenceSolver(Assembler(dom), 4.0e-6, 1.0e-6, nonlinear=True,
                                F_func=lambda t: F).solve()


class _BarraQueCambiaUnApoyo(Truss2D):
    def revise_step(self, U_converged):
        self.nodes[1].fix_dof("uy", 1.0e-6)
        return True


def test_una_revision_que_cambia_restricciones_es_un_error():
    dom, _, n2 = _barra(cls=_BarraQueCambiaUnApoyo)
    F = _desplazamiento(dom, n2, 1.0e-4)
    with pytest.raises(ValueError, match="restricciones o la topología"):
        NonlinearSolver(Assembler(dom), num_steps=2, adaptive=False).solve(F)


def test_un_elemento_que_revisa_no_va_por_lotes():
    """P8 del ADR 0020: un elemento que redefine un gancho de paso guarda
    estado que el núcleo compilado no ve; va por el camino por elemento."""
    from solidum.elements.solid_2d.tri3 import Tri3
    from solidum.materials.elastic_2d import Elastic2D
    from solidum.math.batch.family import element_is_batchable

    class _Tri3QueRevisa(Tri3):
        def revise_step(self, U_converged):
            return False

    class _Tri3QuePrepara(Tri3):
        def prepare_step(self, U_committed):
            pass

    dom = Domain()
    n = [dom.add_node(i + 1, c) for i, c in enumerate([[0, 0], [1, 0], [0, 1]])]
    mat = Elastic2D(E=1.0e9, nu=0.2, hypothesis="plane_stress")
    assert element_is_batchable(Tri3(1, n, mat, thickness=0.1))
    assert not element_is_batchable(_Tri3QueRevisa(2, n, mat, thickness=0.1))
    assert not element_is_batchable(_Tri3QuePrepara(3, n, mat, thickness=0.1))
