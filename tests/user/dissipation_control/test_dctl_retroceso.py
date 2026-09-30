"""Retroceso (*snap-back*) de una barra con daño localizado trazado con
``DissipationArcLengthSolver`` (módulo de usuario ``dissipation_control``).

Mismo caso de validación que ``tests/test_indirect_displacement.py``
(``IndirectDisplacementSolver``), replicado aquí. Barra de longitud ``L`` y
sección ``A`` en 10 elementos ``Truss2D`` con ``IsotropicDamage1D``
(ablandamiento exponencial), empotrada en ``x = 0`` y cargada en el extremo; el
elemento central tiene un umbral ``κ₀`` un 5 % menor, el daño se localiza en él
y el resto se descarga elásticamente. Con la deformación ``ε_w`` del elemento
débil como parámetro, la respuesta exacta es::

    σ(ε_w) = (1 − d(ε_w))·E·ε_w,     d = 1 − (κ₀w/ε_w)·exp(−α(ε_w − κ₀w))
    F = A·σ,                          u = h·ε_w + (L − h)·σ/E,   h = L/n

(con el tope ``DAMAGE_MAX`` sobre ``d``). Con ``(n − 1)·α·κ₀w = 4.275 > 1`` hay
retroceso: el desplazamiento del extremo baja hasta 0.616·u_pico.

Se comprueba, con las tolerancias del test de ``IndirectDisplacementSolver``,
que los puntos posteriores al pico caen sobre la curva exacta y que el
desplazamiento retrocede hasta el mínimo exacto; y, con el espía de
``_dctl_espia.py``, que la rama de ablandamiento se resuelve de verdad con la
restricción de disipación (el modo del problema de Newton, no sólo el
``_mode`` que ve el callback).

Historia. Hasta el 2026-09-29 el solver no seguía el retroceso (deuda #27 de
``docs/STATUS.md``) y los dos tests de la curva iban con ``xfail(strict=True)``:
el predictor del modo de disipación se ensamblaba en el estado consolidado,
donde ``IsotropicDamage1D`` devuelve la secante ``K_s`` y el denominador de la
restricción queda en ``α = −½·F·K_s⁻¹·r_n ≈ 0`` (sólo el residuo de
equilibrio), y el arco cilíndrico se atascaba en la esquina del pico sin dar un
paso con ``ΔE_d`` sobre el umbral de conmutación. El 2026-09-30 el usuario
aprobó el cambio de formulación (spec, §5): (a) el predictor usa la tangente de
carga ``K_L`` del último iterado convergido y (b) el modo de disipación
arranca cuando ``|α_L| > 1e-3·½·|λ_n·F·δu_L|``. Con ello los tests pasan sin
marca y sin tocar tolerancias.

Revisión adversaria del 2026-09-30: todos los puntos con daño, también los de
la rama de rigidez residual tras la saturación, sobre la solución exacta en su
propia ``ε_w``; y el trazado, idéntico con el modelo en otras unidades
(``test_el_trazado_no_depende_de_las_unidades``).

Límite conocido (``test_limite_conocido_con_un_primer_paso_grande``): con
``initial_dlambda = 0.1`` el paso cilíndrico que cruza el pico converge a otro
equilibrio, con los diez elementos ablandándose. No es un fallo del solver sino
la no unicidad del problema incremental con ablandamiento; el test lo deja
fijado como caso conocido.
"""
from __future__ import annotations

import logging
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
sys.path.insert(0, os.path.dirname(__file__))

import solidum  # noqa: E402
from solidum.constants import DAMAGE_MAX  # noqa: E402
from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.truss import Truss2D  # noqa: E402
from solidum.materials.damage_1d import IsotropicDamage1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402

from _dctl_espia import espiar_problemas_de_newton  # noqa: E402

solidum.load_user_module("dissipation_control")
from solidum.user.dissipation_control import DissipationArcLengthSolver  # noqa: E402

E, A, L = 30.0e9, 0.01, 1.0
K0, ALPHA, WEAK = 1.0e-4, 5000.0, 0.95
K0W = WEAK * K0
F_PEAK = A * E * K0W
N_ELEM = 10
H_ELEM = L / N_ELEM

# Tramo muestreado de la curva exacta: ε_w ∈ [κ₀w, 60·κ₀w].
EW_MAX = 60.0 * K0W

# τ₀ ≈ E_d,total/700. La disipación total del elemento débil es
# (½·E·κ₀w² + E·κ₀w/α)·A·h ≈ 0.705 J; con tau_max_factor = 5 el paso llega a
# τ = 5e-3 J, así que la rama de ablandamiento lleva del orden de 140 pasos
# (el test pide más de 100 puntos post-pico).
INITIAL_TAU = 1.0e-3

# Columnas de la historia por paso convergido.
U_END, FORCE, N_DAMAGED, EPS_W, D_W = range(5)


def _bar(f_long: float = 1.0, f_fuerza: float = 1.0):
    """Barra con el elemento central debilitado; carga de referencia = pico.

    ``f_long`` y ``f_fuerza`` cambian las unidades de longitud y de fuerza
    (número en las unidades nuevas = número en SI × factor): ``E`` escala con
    ``f_fuerza/f_long²``, ``A`` con ``f_long²`` y las longitudes con
    ``f_long``; las deformaciones, ``κ₀`` y ``α`` no cambian. Con los factores
    por omisión el modelo es el de SI, bit a bit."""
    E_u = E * f_fuerza / f_long ** 2
    A_u = A * f_long ** 2
    h_u = H_ELEM * f_long
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * h_u, 0.0]) for i in range(N_ELEM + 1)]
    for k in range(N_ELEM):
        kappa_0 = K0W if k == N_ELEM // 2 else K0
        mat = IsotropicDamage1D(E=E_u, kappa_0=kappa_0, alpha=ALPHA)
        dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]], mat, A=A_u))
    nodes[0].fix_dof('ux', 0.0)
    for node in nodes:
        node.fix_dof('uy', 0.0)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[nodes[-1].dofs['ux']] = F_PEAK * f_fuerza
    weak = dom.elements[N_ELEM // 2 + 1]
    return dom, Assembler(dom), F, nodes[-1], weak


def _exact_curve():
    """(u, F) exactos en la rama de ablandamiento, densamente muestreados."""
    ew = K0W * np.concatenate([np.linspace(1.0, 3.0, 20000),
                               np.geomspace(3.0, EW_MAX / K0W, 20000)])
    d = np.minimum(1.0 - (K0W / ew) * np.exp(-ALPHA * (ew - K0W)), DAMAGE_MAX)
    sigma = (1.0 - d) * E * ew
    return H_ELEM * ew + (L - H_ELEM) * sigma / E, A * sigma


def _trazar(initial_dlambda: float, max_steps: int, *, f_long: float = 1.0,
            f_fuerza: float = 1.0) -> SimpleNamespace:
    """Traza la barra y devuelve la historia por paso convergido.

    - ``H``: columnas ``(u_extremo, F, n_dañados, ε_w, d_w)`` (ver ``U_END`` …
      ``D_W``), con ``d_w`` el daño del elemento débil, siempre en SI aunque
      el modelo se resuelva en otras unidades (``f_long``, ``f_fuerza``: ver
      ``_bar``; ``τ₀`` escala con la energía, ``f_fuerza·f_long``);
    - ``modos``: modo del problema de Newton que produjo cada paso;
    - ``intentos``: ``(modo, convergió)`` de todos los intentos de paso,
      fallidos incluidos;
    - ``error``: el ``RuntimeError`` si el solver aborta (longitud de paso
      mínima); la historia llega entonces hasta el último paso convergido.
    """
    dom, asm, F, end, weak = _bar(f_long, f_fuerza)
    a, b = weak.nodes
    solver = DissipationArcLengthSolver(
        asm, max_iter=30, max_lambda=1.2, max_steps=max_steps,
        initial_dlambda=initial_dlambda, initial_tau=INITIAL_TAU * f_fuerza * f_long,
    )
    intentos = espiar_problemas_de_newton(solver)
    hist, modos = [], []

    def cb(step, U, lam):
        damaged = sum(1 for el in dom.elements.values() if el.state.vars[0]['damage'] > 0.0)
        hist.append((U[end.dofs['ux']] / f_long, lam * F_PEAK, damaged,
                     (U[b.dofs['ux']] - U[a.dofs['ux']]) / (H_ELEM * f_long),
                     weak.state.vars[0]['damage']))
        modos.append(intentos[-1][0])

    logging.disable(logging.CRITICAL)
    error = None
    try:
        solver.solve(F, step_callback=cb)
    except RuntimeError as exc:
        error = exc
    finally:
        logging.disable(logging.NOTSET)
    return SimpleNamespace(H=np.array(hist), modos=np.array(modos), intentos=list(intentos),
                           error=error)


@pytest.fixture(scope="module")
def traza():
    """Trazado con el primer paso del caso de validación (Δλ₁ = 0.02)."""
    return _trazar(initial_dlambda=0.02, max_steps=400)


def _post_peak(H, F_ex):
    """Puntos con daño en la rama de ablandamiento, dentro del tramo muestreado.

    Mismo filtro que el test de ``IndirectDisplacementSolver`` (daño y
    ``F > 1.5·min F_exacta``) más ``ε_w ≤ 60·κ₀w``: el paso de disipación que
    cruza la saturación (``d = DAMAGE_MAX``) salta por la rama de rigidez
    residual ``(1 − DAMAGE_MAX)·E`` hasta ``ε_w ≈ 46·κ₀w``, el solver sigue
    después en modo cilíndrico y sale del tramo en que está muestreada la
    curva exacta.
    """
    return H[(H[:, N_DAMAGED] >= 1) & (H[:, FORCE] > 1.5 * F_ex.min())
             & (H[:, EPS_W] <= EW_MAX)]


def test_solo_se_dana_el_elemento_debil(traza):
    assert int(traza.H[:, N_DAMAGED].max()) == 1


def test_puntos_post_pico_sobre_la_curva_exacta(traza):
    H, error = traza.H, traza.error
    u_ex, F_ex = _exact_curve()
    u_pk = u_ex[0]
    post = _post_peak(H, F_ex)
    assert len(post) > 100, f"{len(post)} puntos post-pico (solver: {error})"
    dist = [np.min(np.hypot((u_ex - u) / u_pk, (F_ex - f) / F_PEAK))
            for u, f in post[:, [U_END, FORCE]]]
    assert max(dist) < 1e-3


def test_todos_los_puntos_con_dano_sobre_la_rama_exacta_en_su_eps_w(traza):
    """Complemento del anterior sin el filtro ``ε_w ≤ 60·κ₀w`` (revisión
    adversaria del 2026-09-30): cada punto con daño se compara con la
    solución exacta evaluada en su propia ``ε_w``, sin muestreo, así que
    entran también los ~190 pasos cilíndricos de la rama de rigidez residual
    tras el salto de saturación (hasta ``ε_w ≈ 230·κ₀w``), que el filtro deja
    fuera del test anterior.

    Tolerancia: la del test de ``IndirectDisplacementSolver`` (1e-3). Cota
    del criterio de convergencia (``rtol_force = 1e-5``): la fuerza axial de
    cada elemento difiere de ``λ·F`` en ``√n·‖R‖ ≤ √10·1e-5·λ·F``; en la rama
    residual (rigidez ``(1 − DAMAGE_MAX)·E``) eso da ``δε_w ≤ √10·1e-5·ε_w`` y
    ``δu/u_pico ≤ h·√10·1e-5·ε_w/(κ₀w·L) ≤ 7.2e-4`` con ``ε_w = 229·κ₀w``; en
    la rama de ablandamiento, 6.7e-5 (ver el test del primer paso grande).
    Medido: 7.0e-6 en ``u`` y 7.7e-6 en ``F``.
    """
    H = traza.H[traza.H[:, N_DAMAGED] >= 1]
    ew = H[:, EPS_W]
    d = np.minimum(1.0 - (K0W / ew) * np.exp(-ALPHA * (ew - K0W)), DAMAGE_MAX)
    sigma = (1.0 - d) * E * ew
    u_ex = H_ELEM * ew + (L - H_ELEM) * sigma / E
    u_pico = K0W * L
    assert H[:, EPS_W].max() > 100.0 * K0W       # la rama residual está incluida
    np.testing.assert_array_less(np.abs(H[:, U_END] - u_ex) / u_pico, 1e-3)
    np.testing.assert_array_less(np.abs(H[:, FORCE] - A * sigma) / F_PEAK, 1e-3)


def test_el_desplazamiento_retrocede_hasta_el_minimo_exacto(traza):
    H, error = traza.H, traza.error
    u_ex, F_ex = _exact_curve()
    u_pk = u_ex[0]
    u_min_exact = u_ex[F_ex < 0.999 * F_PEAK].min()
    u_min = H[H[:, N_DAMAGED] >= 1, U_END].min()
    assert u_min < 0.7 * u_pk, f"u_min/u_pico = {u_min / u_pk:.5f} (solver: {error})"
    assert abs(u_min / u_pk - u_min_exact / u_pk) < 2e-3


def test_la_rama_de_ablandamiento_se_resuelve_con_la_restriccion_de_disipacion(traza):
    """Los pasos de la rama de ablandamiento salen de problemas de Newton en
    modo de disipación, no sólo del ``_mode`` que ve el callback.

    Secuencia esperada (acceptance ``switching_cilindrico_a_disipacion_en_pico``):

    - rama elástica y paso en que empieza el daño: cilíndricos. El arranque
      del modo de disipación necesita un paso convergido con daño en
      crecimiento, cuya tangente de carga ``K_L`` se separe de la secante;
    - desde el paso siguiente y mientras el elemento débil se ablanda
      (``0 < d_w < DAMAGE_MAX``): todos en modo de disipación, sin oscilar;
    - tras la saturación del daño la tangente de carga vuelve a ser la
      secante (``α_L ≈ 0``) y el solver vuelve al modo cilíndrico.
    """
    H, modos = traza.H, traza.modos
    d_w = H[:, D_W]
    ablandando = np.flatnonzero((d_w > 0.0) & (d_w < DAMAGE_MAX))
    saturados = np.flatnonzero(d_w >= DAMAGE_MAX)
    assert len(ablandando) > 100, f"{len(ablandando)} pasos con el elemento débil ablandándose"
    inicio = ablandando[0]
    assert set(modos[:inicio + 1]) == {"cylindrical"}, modos[:inicio + 1]
    assert set(modos[ablandando[1:]]) == {"dissipation"}, (
        f"modos en la rama de ablandamiento: {sorted(set(modos[ablandando[1:]]))}")
    assert np.array_equal(ablandando, np.arange(inicio, ablandando[-1] + 1)), (
        "la rama de ablandamiento del elemento débil no es un tramo contiguo")
    # El paso que cruza la saturación es aún de disipación; los siguientes, no.
    assert len(saturados) > 1
    assert saturados[0] == ablandando[-1] + 1
    assert set(modos[saturados[1:]]) == {"cylindrical"}, sorted(set(modos[saturados[1:]]))


# ---------------------------------------------------------------------------
# Invariancia frente a las unidades
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("f_long, f_fuerza", [
    pytest.param(1.0e3, 1.0, id="mm"),           # energías × 1e3
    pytest.param(1.0e-6, 1.0e-6, id="micro"),    # energías × 1e-12
])
def test_el_trazado_no_depende_de_las_unidades(traza, f_long, f_fuerza):
    """El mismo trazado con el modelo en otras unidades (revisión adversaria
    del 2026-09-30).

    Todo lo que el solver compara es adimensional: el arranque por
    ``|α_L|/(½·|λ_n·F·δu_L|)``, la permanencia por ``|α|/escala``, la vuelta
    por ``ΔE_d/(‖F‖·‖U‖)``, la guarda de ``α`` nulo (relativa a los dos
    términos de ``α``) y el criterio de convergencia. Los dos modelos hacen
    las mismas operaciones con números reescalados por potencias de diez, así
    que la secuencia de intentos (modo y convergencia) debe ser idéntica y los
    puntos deben coincidir salvo redondeo.

    ``mm`` es el cambio de unidades de un factor 1e3. ``micro`` lleva las
    energías a ~1e-12 (el orden de una probeta de micras en SI): hasta esta
    revisión, la guarda ``|α| < ZERO_TOL`` en unidades de energía devolvía el
    paso de ablandamiento a cilíndrico en la esquina del pico y el solver
    abortaba en el paso 24.

    Tolerancia: la única operación no covariante que queda en el camino es el
    ``+ ZERO_TOL`` del ``Δλ`` del predictor cilíndrico (la misma fórmula que
    ``ArcLengthSolver``), que en ``micro`` (``‖K⁻¹·F‖ ≈ 1e-10``) mueve un 1 %
    el punto de partida de Newton; el punto convergido lo fijan la
    restricción y el equilibrio (medido: 1.5e-10 en ``micro``, 4e-13 en
    ``mm``). 1e-8 separa ese redondeo de un umbral dependiente de las
    unidades, que aparecería como otra secuencia de modos o al nivel de la
    tolerancia de convergencia (``rtol_force = 1e-5``). Ochenta intentos
    cubren la rama elástica, los cilíndricos fallidos de la esquina, el
    arranque del modo de disipación y unos cincuenta pasos de ablandamiento.
    """
    otra = _trazar(initial_dlambda=0.02, max_steps=80, f_long=f_long, f_fuerza=f_fuerza)
    assert otra.error is None, f"el solver abortó: {otra.error}"
    n = len(otra.intentos)
    assert otra.intentos == traza.intentos[:n]
    assert "dissipation" in set(otra.modos)
    k = len(otra.H)
    u_pico = K0W * L
    np.testing.assert_array_less(np.abs(otra.H[:, U_END] - traza.H[:k, U_END]) / u_pico, 1e-8)
    np.testing.assert_array_less(np.abs(otra.H[:, FORCE] - traza.H[:k, FORCE]) / F_PEAK, 1e-8)


# ---------------------------------------------------------------------------
# Límite conocido: primer paso grande
# ---------------------------------------------------------------------------

def _eps_softening(F, kappa_0):
    """Deformación de un elemento en su rama de ablandamiento bajo la fuerza
    ``F``: ``F/A = E·κ₀·exp(−α(ε − κ₀))`` (``σ = (1 − d)·E·ε`` con el ``d``
    exponencial de ``IsotropicDamage1D``, antes del tope ``DAMAGE_MAX``)."""
    return kappa_0 - np.log(F / (A * E * kappa_0)) / ALPHA


@pytest.fixture(scope="module")
def traza_paso_grande():
    """Trazado con Δλ₁ = 0.1. Treinta pasos bastan: el salto de rama ocurre
    en el quinto."""
    return _trazar(initial_dlambda=0.1, max_steps=30)


def test_limite_conocido_con_un_primer_paso_grande(traza_paso_grande):
    """Con ``initial_dlambda = 0.1`` el trazado no sigue la rama del retroceso.

    Qué pasa: los pasos cilíndricos crecen (×1.5 hasta 5·Δl₁) y el que cruza
    el pico lleva el predictor elástico, de deformación uniforme, por encima
    del umbral de los diez elementos a la vez (el margen entre ``κ₀w`` y
    ``κ₀`` es del 5 %). Newton converge entonces a un equilibrio con los diez
    elementos en carga activa, ablandándose en serie. Es un equilibrio
    legítimo: con ablandamiento el problema incremental no tiene solución
    única (a una misma fuerza, cada elemento puede cargar o descargar), y el
    arco no tiene cómo elegir la rama localizada. El salto ocurre en un paso
    cilíndrico, antes de que el criterio de arranque del modo de disipación
    pueda actuar: éste necesita un paso convergido con daño en crecimiento.
    Desde ahí el solver sigue, ya en modo de disipación, esa otra rama.

    No es un fallo del solver sino una limitación documentada (spec §10): el
    primer paso debe ser lo bastante pequeño para que el paso que cruza el
    pico no rebase el umbral de los elementos que deben descargarse. Si una
    mejora futura lo evitara, este test fallará y habrá que actualizarlo.
    """
    H, modos, error = traza_paso_grande.H, traza_paso_grande.modos, traza_paso_grande.error
    assert error is None, f"el solver abortó: {error}"
    con_dano = np.flatnonzero(H[:, N_DAMAGED] >= 1)
    salto = con_dano[0]
    # El paso que cruza el pico daña los diez elementos, y es cilíndrico.
    assert int(H[salto, N_DAMAGED]) == N_ELEM
    assert modos[salto] == "cylindrical"
    assert set(modos[salto + 1:]) == {"dissipation"}
    assert len(modos) - salto - 1 > 10
    # Los puntos están sobre la rama con los diez elementos ablandándose:
    #   u(F) = h·ε(F; κ₀w) + (L − h)·ε(F; κ₀).
    # Tolerancia derivada del criterio de convergencia por omisión
    # (rtol_force = 1e-5 sobre ‖R‖ en los DOF libres): la fuerza axial de cada
    # elemento difiere de λ·F en la suma de los residuos nodales desde él
    # hasta el extremo, así que |σ_i·A − λ·F| ≤ √n·‖R‖ ≤ √10·1e-5·λ·F; en
    # ablandamiento |dσ/dε| = α·σ, de modo que δε_i ≤ √10·1e-5/α y
    # δu ≤ L·√10·1e-5/α, es decir, δu/u_pico ≤ √10·1e-5/(α·κ₀w) = 6.7e-5
    # (u_pico = κ₀w·L). 1e-4 deja 1.5× de margen (medido: 4.5e-13).
    F = H[salto:, FORCE]
    u_rama = (H_ELEM * _eps_softening(F, K0W)
              + (L - H_ELEM) * _eps_softening(F, K0))
    u_pico = K0W * L
    np.testing.assert_array_less(np.abs(H[salto:, U_END] - u_rama) / u_pico, 1e-4)
    # Y no es la rama localizada: el desplazamiento no retrocede.
    assert H[salto:, U_END].min() > u_pico
