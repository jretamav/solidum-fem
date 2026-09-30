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
desplazamiento retrocede hasta el mínimo exacto.

Estado (2026-09-29): **no lo sigue** (deuda #27 de ``docs/STATUS.md``); los dos
tests van con ``xfail(strict=True)``. El cambio del predictor a Δλ = τ/α no
basta, porque ese predictor no llega a ejecutarse:

1. El predictor del modo de disipación se ensambla en el estado consolidado,
   donde ``IsotropicDamage1D`` no está en carga activa (ε_eq = κ) y devuelve la
   rigidez secante ``K_s``. Como en daño con descarga secante
   ``F_int(U_n) = K_s·U_n``, el denominador de la restricción resulta
   ``α = ½(λ_n·F·K_s⁻¹F − F·U_n) = −½·F·K_s⁻¹·r_n`` (no confundir con el
   ``ALPHA`` del material), es decir, sólo el residuo de equilibrio: la salvaguarda
   ``|α| < 1e-6·escala`` devuelve el paso al modo cilíndrico. En esta barra (y
   en la de un elemento de ``test_dctl_dissipation_arclength.py``) ningún
   problema de Newton se construye en modo de disipación.
2. El pico es una esquina de la curva (la tangente del elemento débil salta de
   ``E`` a ``−ALPHA·κ₀w·E = −0.475·E``, con ``ALPHA`` el parámetro de
   ablandamiento del material) y el retroceso empieza en él. El arco cilíndrico se
   atasca allí: biseca hasta pasos de ~1e-10 m, cuya disipación
   (ΔE_d ≤ 1.1e-12 J) no supera el umbral de conmutación
   (1e-12·‖F‖‖U‖ ≈ 5.3e-12 J), y aborta al llegar a la longitud de arco mínima.
   Con pasos iniciales mayores (Δλ₁ = 0.1) salta el retroceso y converge a
   otro equilibrio con los diez elementos dañados.

Un prototipo de diagnóstico (fuera del repositorio), con la tangente del último
iterado convergido en el predictor y la conmutación decidida por ese α, sigue
el retroceso con estas mismas tolerancias (distancia máxima 3.2e-4,
u_min/u_pico = 0.616376). Es un cambio de formulación pendiente de decisión del
usuario. Cuando se resuelva, ``strict=True`` hará fallar estos tests hasta que
se retire la marca.
"""
from __future__ import annotations

import logging
import os
import sys

import numpy as np
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))

import solidum  # noqa: E402
from solidum.constants import DAMAGE_MAX  # noqa: E402
from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.truss import Truss2D  # noqa: E402
from solidum.materials.damage_1d import IsotropicDamage1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402

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

XFAIL_27 = pytest.mark.xfail(
    strict=True,
    reason=(
        "Deuda #27: el predictor del modo de disipación usa la tangente del estado "
        "consolidado, que en daño con descarga secante es la secante (α ≈ 0 ⇒ vuelta "
        "al modo cilíndrico), y el arco cilíndrico se atasca en la esquina del pico "
        "sin dar un paso con ΔE_d sobre el umbral de conmutación. Cambio de "
        "formulación pendiente del usuario (2026-09-29)."
    ),
)


def _bar():
    """Barra con el elemento central debilitado; carga de referencia = pico."""
    dom = Domain()
    nodes = [dom.add_node(i + 1, [i * H_ELEM, 0.0]) for i in range(N_ELEM + 1)]
    for k in range(N_ELEM):
        kappa_0 = K0W if k == N_ELEM // 2 else K0
        mat = IsotropicDamage1D(E=E, kappa_0=kappa_0, alpha=ALPHA)
        dom.add_element(Truss2D(k + 1, [nodes[k], nodes[k + 1]], mat, A=A))
    nodes[0].fix_dof('ux', 0.0)
    for node in nodes:
        node.fix_dof('uy', 0.0)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[nodes[-1].dofs['ux']] = F_PEAK
    weak = dom.elements[N_ELEM // 2 + 1]
    return dom, Assembler(dom), F, nodes[-1], weak


def _exact_curve():
    """(u, F) exactos en la rama de ablandamiento, densamente muestreados."""
    ew = K0W * np.concatenate([np.linspace(1.0, 3.0, 20000),
                               np.geomspace(3.0, EW_MAX / K0W, 20000)])
    d = np.minimum(1.0 - (K0W / ew) * np.exp(-ALPHA * (ew - K0W)), DAMAGE_MAX)
    sigma = (1.0 - d) * E * ew
    return H_ELEM * ew + (L - H_ELEM) * sigma / E, A * sigma


@pytest.fixture(scope="module")
def traza():
    """Historia ``(u_extremo, F, n_dañados, ε_w)`` por paso convergido.

    Si el solver aborta (``RuntimeError`` por longitud de paso mínima), la
    historia llega hasta el último paso convergido y el error se guarda para
    los mensajes de los tests.
    """
    dom, asm, F, end, weak = _bar()
    a, b = weak.nodes
    solver = DissipationArcLengthSolver(
        asm, max_iter=30, max_lambda=1.2, max_steps=400,
        initial_dlambda=0.02, initial_tau=INITIAL_TAU,
    )
    hist = []

    def cb(step, U, lam):
        damaged = sum(1 for el in dom.elements.values() if el.state.vars[0]['damage'] > 0.0)
        hist.append((U[end.dofs['ux']], lam * F_PEAK, damaged,
                     (U[b.dofs['ux']] - U[a.dofs['ux']]) / H_ELEM))

    logging.disable(logging.CRITICAL)
    error = None
    try:
        solver.solve(F, step_callback=cb)
    except RuntimeError as exc:
        error = exc
    finally:
        logging.disable(logging.NOTSET)
    return np.array(hist), error


def _post_peak(H, F_ex):
    """Puntos con daño en la rama de ablandamiento, dentro del tramo muestreado.

    Mismo filtro que el test de ``IndirectDisplacementSolver`` (daño y
    ``F > 1.5·min F_exacta``) más ``ε_w ≤ 60·κ₀w``: pasada la saturación
    (``d = DAMAGE_MAX``) este solver sigue en modo cilíndrico por la rama de
    rigidez residual ``(1 − DAMAGE_MAX)·E`` y sale del tramo en que está
    muestreada la curva exacta.
    """
    return H[(H[:, 2] >= 1) & (H[:, 1] > 1.5 * F_ex.min()) & (H[:, 3] <= EW_MAX)]


def test_solo_se_dana_el_elemento_debil(traza):
    H, _ = traza
    assert int(H[:, 2].max()) == 1


@XFAIL_27
def test_puntos_post_pico_sobre_la_curva_exacta(traza):
    H, error = traza
    u_ex, F_ex = _exact_curve()
    u_pk = u_ex[0]
    post = _post_peak(H, F_ex)
    assert len(post) > 100, f"{len(post)} puntos post-pico (solver: {error})"
    dist = [np.min(np.hypot((u_ex - u) / u_pk, (F_ex - f) / F_PEAK)) for u, f in post[:, :2]]
    assert max(dist) < 1e-3


@XFAIL_27
def test_el_desplazamiento_retrocede_hasta_el_minimo_exacto(traza):
    H, error = traza
    u_ex, F_ex = _exact_curve()
    u_pk = u_ex[0]
    u_min_exact = u_ex[F_ex < 0.999 * F_PEAK].min()
    u_min = H[H[:, 2] >= 1, 0].min()
    assert u_min < 0.7 * u_pk, f"u_min/u_pico = {u_min / u_pk:.5f} (solver: {error})"
    assert abs(u_min / u_pk - u_min_exact / u_pk) < 2e-3
