"""Integración IMPL-EX de ``IsotropicDamageCrackBand2D/3D`` y pseudo-tiempo (ADR 0034).

Aceptación ``impl_ex`` de la spec ``IsotropicDamageCrackBand2D`` (§8) y la
pieza genérica del ADR 0034:

- punto material: tangente secante exacta, simétrica y definida positiva;
  convergencia a la integración implícita de primer orden en un quiebre de
  ``κ(t)`` (inicio del daño) y de segundo orden donde ``κ(t)`` es suave;
  disipación nunca negativa con carga, descarga y recarga; fuera de un paso
  el daño queda congelado y el estado no cambia;
- ensamblador: razón ``r = Δτ_{n+1}/Δτ_n`` entregada al material y a la
  familia por lotes, idempotencia tras una bisección, rechazo de un solver
  que no declara su pseudo-tiempo;
- global: en una barra que se ablanda, el Newton converge en una iteración
  por paso, la energía disipada es ``G_f·A`` en todas las mallas (±1 %) y el
  post-proceso reproduce el esfuerzo del paso confirmado; los solvers de
  arco, de Newmark y de diferencias centrales declaran su pseudo-tiempo.

Unidades: N, mm, MPa. E = 30 000 MPa, f_t = 3 MPa, G_f = 0.1 N/mm.
"""
import logging

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.core.domain import Domain
from solidum.elements.solid_2d import Quad4
from solidum.materials.damage_crack_band_2d import IsotropicDamageCrackBand2D
from solidum.materials.damage_crack_band_3d import IsotropicDamageCrackBand3D
from solidum.math.assembly import Assembler
from solidum.math.solvers import (
    ArcLengthSolver,
    CentralDifferenceSolver,
    NewtonNewmarkSolver,
    NonlinearSolver,
)

E, NU, FT, GF = 30000.0, 0.2, 3.0, 0.1
K0 = FT / E
H = 10.0


def mat2(integration="impl-ex", **kw):
    return IsotropicDamageCrackBand2D(E, NU, FT, GF, integration=integration, **kw)


def uniaxial(e):
    """Tracción uniaxial de esfuerzo en tensión plana."""
    return np.array([e, -NU * e, 0.0])


# ---------------------------------------------------------------------------
# Punto material
# ---------------------------------------------------------------------------

def test_parametro_integration():
    assert mat2("implicit").uses_pseudo_time is False
    m = mat2()
    assert m.uses_pseudo_time is True and m.integration == "impl-ex"
    # IMPL-EX: tangente secante, simétrica con cualquier medida.
    assert m.IS_SYMMETRIC is True
    assert mat2("implicit").IS_SYMMETRIC is False          # Rankine implícito
    with pytest.raises(ValueError, match="integration"):
        mat2("explicit")
    m3 = IsotropicDamageCrackBand3D(E, NU, FT, GF, integration="impl-ex")
    assert m3.uses_pseudo_time and m3.IS_SYMMETRIC


@pytest.mark.parametrize("dim", [2, 3])
def test_tangente_secante_exacta(dim):
    m = mat2() if dim == 2 else IsotropicDamageCrackBand3D(E, NU, FT, GF, integration="impl-ex")
    n = m.STRAIN_DIM
    st = m.initial_point_state(H)
    rng = np.random.default_rng(7)
    m.set_pseudo_time_ratio(1.0)
    for _ in range(30):
        eps = (rng.standard_normal(n) + np.r_[3.0, np.zeros(n - 1)]) * 2.0e-4
        sig, C, st = m.compute_state(eps, st)
        d = st["damage"]
        np.testing.assert_array_equal(C, (1.0 - d) * m.elastic_base.C)
        np.testing.assert_allclose(sig, C @ eps, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(C, C.T, rtol=0, atol=0)
        assert np.linalg.eigvalsh(C).min() > 0.0
    assert st["damage"] > 0.5


def _camino(n, integration, kappa_de_t):
    m = mat2(integration)
    st = m.initial_point_state(H)
    m.set_pseudo_time_ratio(1.0)
    out = []
    for t in np.linspace(0.0, 1.0, n + 1)[1:]:
        sig, _, st = m.compute_state(uniaxial(kappa_de_t(t)), st)
        out.append(sig[0])
    return np.array(out)


@pytest.mark.parametrize("forma, orden", [("lineal", 1.0), ("cuadratica", 2.0)])
def test_convergencia_a_la_implicita(forma, orden):
    """Error máximo de σ frente a la implícita con pasos Δ, Δ/2, Δ/4, Δ/8 de
    pseudo-tiempo. Camino que arranca en κ_0: lineal (quiebre de κ(t) en el
    inicio del daño → primer orden) o cuadrático (κ(t) suave → segundo)."""
    if forma == "lineal":
        def kappa(t):
            return K0 * (1.0 + 39.0 * t)
    else:
        def kappa(t):
            return K0 * (1.0 + 39.0 * t * t)
    errores = []
    for n in (100, 200, 400, 800):
        errores.append(np.abs(_camino(n, "impl-ex", kappa) - _camino(n, "implicit", kappa)).max())
    razones = np.array(errores[:-1]) / np.array(errores[1:])
    np.testing.assert_allclose(razones, 2.0 ** orden, rtol=0.15)


def test_disipacion_nunca_negativa():
    """Carga hasta el ablandamiento, descarga parcial y recarga, con una
    razón que hace pasarse a la extrapolación al empezar la descarga."""
    m = mat2()
    st = m.initial_point_state(H)
    Ce = m.elastic_base.C
    camino = np.r_[np.linspace(0.0, 8.0, 41)[1:], np.linspace(8.0, 3.0, 11)[1:],
                   np.linspace(3.0, 15.0, 31)[1:]] * K0
    d_prev = 0.0
    for k, e in enumerate(camino):
        m.set_pseudo_time_ratio(2.0 if k % 7 == 0 else 1.0)
        eps = uniaxial(e)
        _, _, st = m.compute_state(eps, st)
        d = st["damage"]
        assert d >= d_prev
        psi0 = 0.5 * eps @ Ce @ eps
        assert psi0 * (d - d_prev) >= 0.0
        assert st["damage_error"] >= 0.0
        d_prev = d
    assert d_prev > 0.9


def test_fuera_de_un_paso_el_dano_queda_congelado():
    m = mat2()
    st = m.initial_point_state(H)
    m.set_pseudo_time_ratio(1.0)
    for e in np.linspace(0.0, 6.0, 31)[1:] * K0:
        _, _, st = m.compute_state(uniaxial(e), st)
    m.set_pseudo_time_ratio(None)
    eps = uniaxial(9.0 * K0)
    sig, C, nuevo = m.compute_state(eps, st)
    assert nuevo == st
    np.testing.assert_array_equal(C, (1.0 - st["damage"]) * m.elastic_base.C)
    np.testing.assert_allclose(sig, C @ eps, rtol=1e-14, atol=1e-14)


def test_extrapolacion_con_la_razon():
    """κ̃ = κ_n + r·(κ_n − κ_{n−1}): con r = 0 el daño del paso es d(κ_n)."""
    m = mat2()
    st = m.initial_point_state(H)
    m.set_pseudo_time_ratio(1.0)
    for e in (1.5, 2.0, 2.5):
        _, _, st = m.compute_state(uniaxial(e * K0), st)
    ref = mat2("implicit")
    st_ref = ref.initial_point_state(H)
    _, _, st_ref = ref.compute_state(uniaxial(2.5 * K0), st_ref)
    m.set_pseudo_time_ratio(0.0)
    _, _, nuevo = m.compute_state(uniaxial(2.6 * K0), st)
    assert nuevo["damage"] == pytest.approx(st_ref["damage"], rel=1e-14)
    assert nuevo["kappa_prev"] == st["kappa"]
    assert nuevo["kappa"] == pytest.approx(2.6 * K0, rel=1e-14)


# ---------------------------------------------------------------------------
# Ensamblador (ADR 0034)
# ---------------------------------------------------------------------------

def _un_elemento(integration="impl-ex", batch=None):
    dom = Domain()
    pts = [(0, 0), (H, 0), (H, H), (0, H)]
    nodes = [dom.add_node(i + 1, list(p)) for i, p in enumerate(pts)]
    mat = mat2(integration)
    dom.add_element(Quad4(1, nodes, mat, thickness=1.0))
    nodes[0].fix_dof("ux", 0.0)
    nodes[0].fix_dof("uy", 0.0)
    nodes[3].fix_dof("ux", 0.0)
    for nd in (nodes[1], nodes[2]):
        nd.fix_dof("ux", 0.0)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom) if batch is None else Assembler(dom, batch=batch)
    return dom, asm, mat


@pytest.mark.parametrize("batch", [True, False])
def test_razon_entregada_por_el_ensamblador(batch):
    dom, asm, mat = _un_elemento(batch=batch)
    U = np.zeros(dom.total_dofs)
    asm.assemble_non_linear_system(U)           # construye las familias
    fams = [f for f in asm.families if f.material is mat]
    assert (len(fams) == 1) == batch

    def razon():
        r = mat.batch_params()[-1]
        for f in fams:
            assert f.mat_params[-1] == r
        return r

    assert razon() == -1.0                       # antes del primer paso
    asm.prepare_all_steps(U, pseudo_time_increment=0.1)
    assert razon() == 1.0 and asm.pseudo_time_increment == 0.1
    asm.prepare_all_steps(U, pseudo_time_increment=0.05)   # bisección, sin referencia aún
    assert razon() == 1.0
    asm.assemble_non_linear_system(U)
    asm.commit_all_states()                      # Δτ_n = 0.05
    assert razon() == -1.0 and asm.pseudo_time_increment is None
    asm.prepare_all_steps(U, pseudo_time_increment=0.1)
    assert razon() == pytest.approx(2.0)
    asm.prepare_all_steps(U, pseudo_time_increment=0.025)  # bisección
    assert razon() == pytest.approx(0.5)
    asm.prepare_all_steps(U)
    asm.set_pseudo_time_increment(0.05)          # tras el predictor (arco)
    assert razon() == pytest.approx(1.0)
    with pytest.raises(ValueError, match="pseudo-tiempo"):
        asm.set_pseudo_time_increment(0.0)


def test_mismo_intento_tras_una_biseccion():
    dom, asm, mat = _un_elemento()
    U_step = np.zeros(dom.total_dofs)
    asm.prepare_all_steps(U_step, pseudo_time_increment=1.0)
    for nd in dom.nodes.values():
        if nd.coordinates[1] == H:
            U_step[nd.dofs["uy"]] = 4.0 * K0 * H
    asm.assemble_non_linear_system(U_step)
    asm.commit_all_states()
    U_next = U_step * 1.5
    asm.prepare_all_steps(U_step, pseudo_time_increment=0.8)
    K1, F1 = asm.assemble_non_linear_system(U_next)
    asm.prepare_all_steps(U_step, pseudo_time_increment=0.4)
    asm.assemble_non_linear_system(U_next)
    asm.prepare_all_steps(U_step, pseudo_time_increment=0.8)
    K2, F2 = asm.assemble_non_linear_system(U_next)
    np.testing.assert_array_equal(F1, F2)
    np.testing.assert_array_equal(K1.toarray(), K2.toarray())


def test_solver_sin_pseudo_tiempo_se_rechaza():
    dom, asm, _ = _un_elemento()
    U = np.zeros(dom.total_dofs)
    asm.prepare_all_steps(U)                     # intento sin incremento declarado
    asm.assemble_non_linear_system(U)
    with pytest.raises(ValueError, match="IMPL-EX"):
        asm.commit_all_states()
    # Con la implícita no hace falta el pseudo-tiempo.
    dom, asm, _ = _un_elemento("implicit")
    asm.prepare_all_steps(U)
    asm.assemble_non_linear_system(U)
    asm.commit_all_states()
    # Fuera de un paso (análisis lineal, commit final de solidum.run): sin error.
    dom, asm, _ = _un_elemento()
    asm.assemble_non_linear_system(U)
    asm.commit_all_states()


# ---------------------------------------------------------------------------
# Barra que se ablanda (global)
# ---------------------------------------------------------------------------

L_BAR, H_BAR, T_BAR = 100.0, 20.0, 10.0


def barra(nx, integration, nu=0.0, density=None):
    """Barra en tracción con la columna central debilitada (f_t·0.95);
    desplazamiento impuesto en el extremo derecho."""
    dx = L_BAR / nx
    ny = int(round(H_BAR / dx))
    dom = Domain()
    mat = IsotropicDamageCrackBand2D(E, nu, FT, GF, integration=integration, density=density)
    weak = IsotropicDamageCrackBand2D(E, nu, 0.95 * FT, GF, integration=integration,
                                      density=density)
    nodes = {}
    for j in range(ny + 1):
        for i in range(nx + 1):
            nodes[(i, j)] = dom.add_node(len(nodes) + 1, [i * dx, j * H_BAR / ny])
    for j in range(ny):
        for i in range(nx):
            m = weak if i == nx // 2 else mat
            dom.add_element(Quad4(len(dom.elements) + 1,
                                  [nodes[(i, j)], nodes[(i + 1, j)], nodes[(i + 1, j + 1)],
                                   nodes[(i, j + 1)]], m, thickness=T_BAR))
    for j in range(ny + 1):
        nodes[(0, j)].fix_dof("ux", 0.0)
    nodes[(0, 0)].fix_dof("uy", 0.0)
    right = [nodes[(nx, j)] for j in range(ny + 1)]
    return dom, right


def _con_desplazamiento(dom, right, u_max, steps):
    for nd in right:
        nd.fix_dof("ux", u_max)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    dofs = [nd.dofs["ux"] for nd in right]
    curva = [(0.0, 0.0)]

    def cb(step, U, lam):
        _, F = asm.assemble_non_linear_system(U)
        curva.append((U[dofs[0]], F[dofs].sum()))

    return asm, cb, curva


class _Iteraciones(logging.Handler):
    def __init__(self):
        super().__init__(logging.INFO)
        self.iters = []

    def emit(self, record):
        msg = record.getMessage()
        if msg.startswith("[PASO") and "iter" in msg:
            self.iters.append(int(msg.rsplit("|", 1)[1].split()[0]))


def _resolver(asm, cb, steps, n):
    log = logging.getLogger("solidum")
    nivel = log.level
    contador = _Iteraciones()
    log.addHandler(contador)
    log.setLevel(logging.INFO)
    try:
        NonlinearSolver(asm, num_steps=steps).solve(np.zeros(n), step_callback=cb)
    finally:
        log.removeHandler(contador)
        log.setLevel(nivel)
    return contador.iters


def test_newton_de_una_iteracion_en_el_ablandamiento():
    """Con IMPL-EX el problema de cada paso es lineal: el Newton converge en
    una iteración, más la de comprobación."""
    dom, right = barra(10, "impl-ex")
    asm, cb, curva = _con_desplazamiento(dom, right, 0.4, 400)
    iters = _resolver(asm, cb, 400, dom.total_dofs)
    # Una resolución y la de comprobación (corrección nula) en cada paso.
    assert len(iters) == 400 and set(iters) == {2}
    # La implícita necesita más (y bisecciones en el pico).
    dom, right = barra(10, "implicit")
    asm, cb, _ = _con_desplazamiento(dom, right, 0.4, 400)
    assert max(_resolver(asm, cb, 400, dom.total_dofs)) > 2


@pytest.mark.parametrize("nx", [5, 10, 20])
def test_objetividad_y_energia_con_impl_ex(nx):
    """D = G_f·A a ±1 % y la carga en u = 0.05 mm igual a la implícita a
    ±0.5 % (paso de 0.125 µm; con pasos más gruesos la extrapolación sobrepasa
    el pico y disipa de más: 1.01–1.03·G_f·A con 0.5 µm)."""
    dom, right = barra(nx, "impl-ex")
    asm, cb, curva = _con_desplazamiento(dom, right, 0.4, 3200)
    _resolver(asm, cb, 3200, dom.total_dofs)
    c = np.array(curva)
    D = np.trapezoid(c[:, 1], c[:, 0]) - 0.5 * c[-1, 1] * c[-1, 0]
    assert D / (GF * H_BAR * T_BAR) == pytest.approx(1.0, abs=0.01)
    dom_i, right_i = barra(nx, "implicit")
    asm_i, cb_i, curva_i = _con_desplazamiento(dom_i, right_i, 0.4, 400)
    _resolver(asm_i, cb_i, 400, dom_i.total_dofs)
    ci = np.array(curva_i)
    assert np.interp(0.05, c[:, 0], c[:, 1]) == pytest.approx(
        np.interp(0.05, ci[:, 0], ci[:, 1]), rel=5e-3)


def test_post_proceso_reproduce_el_esfuerzo_confirmado():
    dom, right = barra(5, "impl-ex")
    asm, cb, _ = _con_desplazamiento(dom, right, 0.05, 100)
    U = NonlinearSolver(asm, num_steps=100).solve(np.zeros(dom.total_dofs))
    hubo_dano = False
    for e in dom.elements.values():
        gs = e.compute_gauss_state(U)
        for ip, sig in enumerate(np.asarray(gs["stress"])):
            st = e.state.vars[ip]
            eps = np.asarray(gs["strain"])[ip]
            np.testing.assert_allclose(
                sig, (1.0 - st["damage"]) * (e.material.elastic_base.C @ eps), rtol=1e-12, atol=1e-12)
            hubo_dano |= st["damage"] > 0.5
    assert hubo_dano


def test_arco_declara_su_pseudo_tiempo():
    """Barra con fuerza en el extremo y longitud de arco: llega al
    ablandamiento sin el error de pseudo-tiempo no declarado."""
    dom, right = barra(5, "impl-ex")
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    F = np.zeros(dom.total_dofs)
    for nd in right:
        F[nd.dofs["ux"]] = FT * H_BAR * T_BAR / len(right)
    log = logging.getLogger("solidum")
    nivel = log.level
    log.setLevel(logging.ERROR)
    try:
        solver = ArcLengthSolver(asm, control=[(right[0].id, "ux", 1.0)], control_target=0.05,
                                 max_lambda=10.0)
        U = solver.solve(F)
    finally:
        log.setLevel(nivel)
    assert U[right[0].dofs["ux"]] == pytest.approx(0.05, rel=1e-6)
    assert max(st["damage"] for e in dom.elements.values() for st in e.state.vars) > 0.5


@pytest.mark.parametrize("solver", ["newton_newmark", "central_difference"])
def test_dinamica_declara_su_pseudo_tiempo(solver):
    """Barra con masa y fuerza creciente en el extremo: los solvers dinámicos
    declaran Δt como pseudo-tiempo y el daño evoluciona con IMPL-EX."""
    dom, right = barra(5, "impl-ex", density=2.4e-9)       # t/mm³
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    dofs = [nd.dofs["ux"] for nd in right]
    P = 1.5 * FT * H_BAR * T_BAR / len(right)

    def F(t):
        f = np.zeros(dom.total_dofs)
        f[dofs] = P * min(t / 1.0e-4, 1.0)
        return f

    log = logging.getLogger("solidum")
    nivel = log.level
    log.setLevel(logging.ERROR)
    try:
        if solver == "newton_newmark":
            NewtonNewmarkSolver(asm, t_end=1.5e-4, dt=1.0e-6, F_func=F).solve()
        else:
            CentralDifferenceSolver(asm, t_end=1.5e-4, dt=1.0e-6, F_func=F, nonlinear=True).solve()
    finally:
        log.setLevel(nivel)
    assert max(st["damage"] for e in dom.elements.values() for st in e.state.vars) > 0.5
