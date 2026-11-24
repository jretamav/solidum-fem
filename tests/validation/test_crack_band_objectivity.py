"""Objetividad de malla de la banda de grieta (``IsotropicDamageCrackBand2D/3D``).

Aceptación ``mesh_objectivity`` de las specs: barra en tracción con una
columna (capa) de elementos debilitada (f_t·0.95) en el centro, celdas
cuadradas (cúbicas) de lado Δ = 20, 10, 5 mm, desplazamiento impuesto en un
extremo. La deformación se localiza en la columna debilitada, que disipa
``G_f/h`` por unidad de volumen; con ``h = Δ`` la energía disipada al final
es ``G_f·A_grieta`` con cualquier malla y la rama posterior al pico es la
misma.

- ν = 0: D/(G_f·A) = 1 a ±1 % en los diez sólidos y la carga en u = 0.05 mm
  igual entre mallas a ±1 %.
- ν = 0.2 (Quad4): D/(G_f·A) = 1 − ν² a ±1 %. Las columnas vecinas, que se
  descargan, impiden la contracción lateral de la banda: la banda trabaja con
  ε_yy ≈ 0, Rankine ve ε_eq = ε_xx/(1−ν²) y el trabajo de σ_xx en la banda
  es (1−ν²)·G_f/h. Es física del modelo continuo (igual en todas las mallas),
  no un error de discretización.
- Factor de ``h = c·Ω_e^(1/d)`` (regla Q2 de la spec: si la barra se aparta
  más de un 5 % de G_f·A, se corrige el factor):
  - triángulos y tetraedros: con el ``√A`` / ``∛V`` de Abaqus la columna
    disipaba 1.41 / 1.81·G_f·A (los dos triángulos o seis tetraedros de un
    cuadrado o cubo forman la banda); con √(2A) (Rots 1988) y ∛(6V), 1.00
    (decisión del usuario 2026-11-11, Q2);
  - cuadráticos: con el ½ de Abaqus disipaban 1.97 (``Quad8``, ``Quad9``,
    ``Hex20``, ``Hex27``), 2.68 (``Tri6``) y 3.23 (``Tet10``)·G_f·A: con el
    esfuerzo uniforme de la barra la deformación se localiza en el elemento
    entero. Usan el ancho del elemento lineal de su celda (usuario,
    2026-11-11).
- Factor de localización ``α`` (``CRACK_BAND_LOCALIZATION``, usuario,
  2026-11-12): con gradiente, en ``Quad8`` y ``Quad9`` sólo se ablandan dos
  de sus tres filas de puntos de Gauss 3×3 (Jirásek y Bauer 2012, §4; medido
  0.72 en la viga de Petersson) y ``h = (13/18)·√A``. En esta barra, de
  esfuerzo uniforme, el elemento se ablanda entero y disipa ``G_f·A/α =
  1.385·G_f·A``: el caso en que ``α`` sobrestima la energía, documentado aquí
  a propósito. El caso que ``α`` corrige es el de flexión: viga entallada
  (``test_viga_entallada_*``), donde ``Quad8`` disipa lo mismo que ``Quad4``.
  ``Tri6`` y los cuadráticos 3D, sin factor establecido (``α = 1``), siguen
  dando 1.00 en la barra.

Modelo 2, ``TensionCompressionDamage2D/3D`` (aceptación ``bar_in_tension`` de
su spec): la misma barra con los diez sólidos. En tracción uniaxial su ley es
exacta, ``σ = (1 − d⁺)·E·ε`` (spec §1), y en la banda con Poisson los dos
esfuerzos principales efectivos son de tracción, así que las energías son las
del modelo 1: ``G_f·A``, ``G_f·A/α`` en ``Quad8``/``Quad9`` y ``(1 − ν²)·G_f·A``
con ν = 0.2; también con IMPL-EX y pasos de 0.125 µm (con 1 µm IMPL-EX disipa
de más, 1.04–1.09·G_f·A según la malla, igual en los dos modelos).

Unidades: N, mm, MPa. E = 30 000 MPa, f_t = 3 MPa, G_f = 0.1 N/mm.
"""
import logging

import numpy as np
import pytest

import solidum  # noqa: F401
from solidum.elements.solid_2d import Quad4, Quad8, Quad9, Tri3, Tri6
from solidum.elements.solid_3d import Hex8, Hex20, Hex27, Tet4, Tet10
from solidum.materials.damage_crack_band_2d import IsotropicDamageCrackBand2D
from solidum.materials.damage_crack_band_3d import IsotropicDamageCrackBand3D
from solidum.materials.damage_tension_compression_2d import TensionCompressionDamage2D
from solidum.materials.damage_tension_compression_3d import TensionCompressionDamage3D
from solidum.math.assembly import Assembler
from solidum.math.solvers import NonlinearSolver
from solidum.mesh import structured_mesh

E, FT, GF = 30000.0, 3.0, 0.1
L, H, T = 100.0, 20.0, 10.0
U_MAX, STEPS = 0.4, 400
# Materiales (2D, 3D) de cada modelo con banda de grieta.
MODELO_1 = (IsotropicDamageCrackBand2D, IsotropicDamageCrackBand3D)
MODELO_2 = (TensionCompressionDamage2D, TensionCompressionDamage3D)


def _solve(dom, right_nodes, u_max=U_MAX, steps=STEPS):
    dom.generate_equation_numbers(verbose=False)
    right = [nd.dofs["ux"] for nd in right_nodes]
    asm = Assembler(dom)
    curve = [(0.0, 0.0)]

    def cb(step, U, lam):
        _, F = asm.assemble_non_linear_system(U)
        curve.append((U[right[0]], F[right].sum()))

    level = logging.getLogger("solidum").level
    logging.getLogger("solidum").setLevel(logging.ERROR)
    try:
        NonlinearSolver(asm, num_steps=steps).solve(np.zeros(dom.total_dofs), step_callback=cb)
    finally:
        logging.getLogger("solidum").setLevel(level)
    c = np.array(curve)
    W = np.trapezoid(c[:, 1], c[:, 0])
    return c, W - 0.5 * c[-1, 1] * c[-1, 0]


def bar(cls, nx, nu=0.0, modelo=MODELO_1, integration="implicit", steps=STEPS):
    """Barra de L × H (× H) con la columna central de celdas debilitada."""
    dim = len(cls.REFERENCE_NODES[0])
    dx = L / nx
    n = int(round(H / dx))
    if dim == 2:
        mat, weak = (modelo[0](E, nu, f, GF, integration=integration) for f in (FT, 0.95 * FT))
        kw = {"thickness": T}
        mesh = structured_mesh(cls, (nx, n), size=(L, H), material=mat, **kw)
        GA = GF * H * T
    else:
        mat, weak = (modelo[1](E, nu, f, GF, integration=integration) for f in (FT, 0.95 * FT))
        kw = {}
        mesh = structured_mesh(cls, (nx, n, n), size=(L, H, H), material=mat)
        GA = GF * H * H
    dom = mesh.domain
    x_weak = (nx // 2 + 0.5) * dx
    for eid, el in list(dom.elements.items()):
        xc = np.mean([nd.coordinates[0] for nd in el.nodes])
        if abs(xc - x_weak) < 0.5 * dx:
            dom.replace_element(cls(eid, el.nodes, weak, **kw))
    for nd in mesh.nodes_on("xmin"):
        nd.fix_dof("ux", 0.0)
    right = mesh.nodes_on("xmax")
    for nd in right:
        nd.fix_dof("ux", U_MAX)
    corners = {tuple(np.round(nd.coordinates, 9)): nd for nd in mesh.nodes_on("xmin")}
    origin = corners[(0.0,) * dim]
    origin.fix_dof("uy", 0.0)
    if dim == 3:
        origin.fix_dof("uz", 0.0)
        corners[(0.0, H, 0.0)].fix_dof("uz", 0.0)
    return _solve(dom, right, steps=steps), GA


def _check(results, expected_ratio):
    p05 = []
    for nx, ((curve, D), GA) in results.items():
        assert D / GA == pytest.approx(expected_ratio, rel=1e-2), nx
        # El pico es f_t·0.95·A; con Δu = 1 µm la curva lo muestrea por debajo.
        P_peak = 0.95 * FT * GA / GF
        assert 0.97 * P_peak < curve[:, 1].max() <= P_peak * (1 + 1e-9)
        p05.append(np.interp(0.05, curve[:, 0], curve[:, 1]))
    assert max(p05) / min(p05) - 1.0 < 1e-2


@pytest.mark.parametrize("cls", [Quad4, Tri3, Tri6], ids=lambda c: c.__name__)
def test_barra_2d_energia_objetiva(cls):
    _check({nx: bar(cls, nx) for nx in (5, 10, 20)}, 1.0)


@pytest.mark.parametrize("cls", [Quad8, Quad9], ids=lambda c: c.__name__)
def test_barra_2d_cuadraticos_con_factor_disipan_1_sobre_alfa(cls):
    """Esfuerzo uniforme: el elemento se ablanda entero y α = 13/18 da
    G_f·A/α (objetivo entre mallas, pero 1.385 veces G_f·A)."""
    _check({nx: bar(cls, nx) for nx in (5, 10, 20)}, 18.0 / 13.0)


def test_barra_2d_con_poisson_disipa_1_menos_nu2():
    nu = 0.2
    _check({nx: bar(Quad4, nx, nu) for nx in (5, 10, 20)}, 1.0 - nu**2)


@pytest.mark.parametrize("cls", [Hex8, Tet4, Hex20, Hex27, Tet10], ids=lambda c: c.__name__)
def test_barra_3d_energia_objetiva(cls):
    _check({nx: bar(cls, nx) for nx in (5, 10)}, 1.0)


# ---------------------------------------------------------------------------
# Modelo 2: TensionCompressionDamage2D/3D (bar_in_tension de su spec)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cls", [Quad4, Tri3, Tri6], ids=lambda c: c.__name__)
def test_modelo2_barra_2d_energia_objetiva(cls):
    _check({nx: bar(cls, nx, modelo=MODELO_2) for nx in (5, 10, 20)}, 1.0)


@pytest.mark.parametrize("cls", [Quad8, Quad9], ids=lambda c: c.__name__)
def test_modelo2_barra_2d_cuadraticos_con_factor_disipan_1_sobre_alfa(cls):
    _check({nx: bar(cls, nx, modelo=MODELO_2) for nx in (5, 10, 20)}, 18.0 / 13.0)


def test_modelo2_barra_2d_con_poisson_disipa_1_menos_nu2():
    nu = 0.2
    _check({nx: bar(Quad4, nx, nu, modelo=MODELO_2) for nx in (5, 10, 20)}, 1.0 - nu**2)


@pytest.mark.parametrize("cls", [Hex8, Tet4, Hex20, Hex27, Tet10], ids=lambda c: c.__name__)
def test_modelo2_barra_3d_energia_objetiva(cls):
    _check({nx: bar(cls, nx, modelo=MODELO_2) for nx in (5, 10)}, 1.0)


def test_modelo2_barra_impl_ex():
    """IMPL-EX con Newton completo dentro del paso (spec §8), con los criterios
    y el paso de 0.125 µm de la prueba del modelo 1
    (``test_damage_crack_band_implex``): D = G_f·A a ±1 % y la carga en
    u = 0.05 mm igual a la implícita a ±0.5 %. El pico no se compara con
    f_t·A: la extrapolación retrasa el daño un paso y lo rebasa un poco. Con
    pasos más gruesos disipa de más, igual en los dos modelos (abajo)."""
    for nx in (5, 10, 20):
        (c, D), GA = bar(Quad4, nx, modelo=MODELO_2, integration="impl-ex", steps=3200)
        assert D / GA == pytest.approx(1.0, abs=0.01), nx
        (ci, _), _ = bar(Quad4, nx, modelo=MODELO_2)
        assert np.interp(0.05, c[:, 0], c[:, 1]) == pytest.approx(
            np.interp(0.05, ci[:, 0], ci[:, 1]), rel=5e-3), nx


def test_modelo2_barra_impl_ex_igual_al_modelo_1():
    """Con ν = 0 la barra es uniaxial: los dos modelos dan la misma curva
    también con el paso grueso (1 µm). Con ese paso IMPL-EX rebasa el pico un
    13–15 % y disipa 1.036, 1.062 y 1.094·G_f·A con 5, 10 y 20 columnas: el
    exceso crece al refinar (reporte de validación del 2026-11-18,
    ``informes/reportes/reporte_modelo2_barra_2026-11-18.pdf``; hipótesis: el error
    depende de Δu/h, no de Δu)."""
    (c1, D1), _ = bar(Quad4, 5, modelo=MODELO_1, integration="impl-ex")
    (c2, D2), _ = bar(Quad4, 5, modelo=MODELO_2, integration="impl-ex")
    np.testing.assert_allclose(c2, c1, rtol=1e-9, atol=1e-9 * np.abs(c1).max())
    assert D2 == pytest.approx(D1, rel=1e-9)


# ---------------------------------------------------------------------------
# Modelo 2 en compresión (bar_in_compression de su spec)
# ---------------------------------------------------------------------------
#
# Barra de longitud fija L = 200 mm y sección cuadrada de 40 mm (celdas
# cuadradas de 40, 20 y 10 mm), ν = 0, columna central con f_c·0.95,
# acortamiento impuesto. Con ν = 0 el esfuerzo es uniforme a lo largo de la
# barra, σ = F/A, así que el trabajo sobre la columna de la banda es
# ∫F·dδ_banda, con δ_banda el acortamiento de la columna. Desde el pico, la
# columna recorre la rama descendente de la parábola de Feenstra, de área
# G_c/h: el trabajo debe ser G_c·A con cualquier malla. Sin retroceso: con
# h ≤ 40 mm la pendiente de la rama, a lo sumo 2·f_c/(κ_u − κ_c) = 80·h
# (MPa), cumple (L − h)·80·h/E < h. El trabajo externo TOTAL no es objetivo
# (1.116, 1.081, 1.064·G_c·A con h = 40, 20, 10 mm): además de G_c·A tiene
# dos partes, la energía de la banda antes del pico (0.079, 0.039, 0.020·G_c·A,
# proporcional a h) y el daño de la rama ascendente en el resto de la barra
# (0.037, 0.042, 0.044·G_c·A, proporcional a L − h); con h → 0 el total
# tiende a 1.046·G_c·A (spec §10; verificación y validación del 2026-11-18).

FC, GC = 30.0, 15.0                    # MPa, N/mm
LC, HC = 200.0, 40.0
UC_MAX, STEPS_C = 1.0, 1000


def bar_compression(cls, nx, integration="implicit", steps=STEPS_C, u_max=UC_MAX):
    """(curva [(acortamiento, fuerza)], trabajo sobre la banda desde su pico, G_c·A)."""
    dim = len(cls.REFERENCE_NODES[0])
    dx = LC / nx
    n = int(round(HC / dx))
    if dim == 2:
        mat, weak = (TensionCompressionDamage2D(E, 0.0, FT, GF, f_c=f, G_c=GC,
                                                integration=integration)
                     for f in (FC, 0.95 * FC))
        kw = {"thickness": T}
        mesh = structured_mesh(cls, (nx, n), size=(LC, HC), material=mat, **kw)
        A = HC * T
    else:
        mat, weak = (TensionCompressionDamage3D(E, 0.0, FT, GF, f_c=f, G_c=GC,
                                                integration=integration)
                     for f in (FC, 0.95 * FC))
        kw = {}
        mesh = structured_mesh(cls, (nx, n, n), size=(LC, HC, HC), material=mat)
        A = HC * HC
    dom = mesh.domain
    x_weak = (nx // 2 + 0.5) * dx
    for eid, el in list(dom.elements.items()):
        xc = np.mean([nd.coordinates[0] for nd in el.nodes])
        if abs(xc - x_weak) < 0.5 * dx:
            dom.replace_element(cls(eid, el.nodes, weak, **kw))
    for nd in mesh.nodes_on("xmin"):
        nd.fix_dof("ux", 0.0)
    right = mesh.nodes_on("xmax")
    for nd in right:
        nd.fix_dof("ux", -u_max)
    corners = {tuple(np.round(nd.coordinates, 9)): nd for nd in mesh.nodes_on("xmin")}
    origin = corners[(0.0,) * dim]
    origin.fix_dof("uy", 0.0)
    if dim == 3:
        origin.fix_dof("uz", 0.0)
        corners[(0.0, HC, 0.0)].fix_dof("uz", 0.0)
    dom.generate_equation_numbers(verbose=False)
    faces = [[nd for nd in dom.nodes.values() if abs(nd.coordinates[0] - x) < 1e-9]
             for x in (x_weak - 0.5 * dx, x_weak + 0.5 * dx)]
    dofs_r = [nd.dofs["ux"] for nd in right]
    dofs_faces = [[nd.dofs["ux"] for nd in f] for f in faces]
    asm = Assembler(dom)
    rows = [(0.0, 0.0, 0.0)]

    def cb(step, U, lam):
        _, F = asm.assemble_non_linear_system(U)
        band = np.mean(U[dofs_faces[0]]) - np.mean(U[dofs_faces[1]])   # acortamiento > 0
        rows.append((-U[dofs_r[0]], -F[dofs_r].sum(), band))

    level = logging.getLogger("solidum").level
    logging.getLogger("solidum").setLevel(logging.ERROR)
    try:
        NonlinearSolver(asm, num_steps=steps).solve(np.zeros(dom.total_dofs), step_callback=cb)
    finally:
        logging.getLogger("solidum").setLevel(level)
    c = np.array(rows)
    k = int(np.argmax(c[:, 1]))
    W_band = np.trapezoid(c[k:, 1], c[k:, 2])
    return c, W_band, GC * A


def _check_compression(cls, nxs, expected_ratio, **kw):
    W = []
    for nx in nxs:
        c, Wb, GA = bar_compression(cls, nx, **kw)
        assert Wb / GA == pytest.approx(expected_ratio, rel=0.02), nx
        A = GA / GC
        assert c[:, 1].max() == pytest.approx(0.95 * FC * A, rel=1e-3), nx
        assert c[-1, 1] < 1e-3 * c[:, 1].max(), nx
        W.append(Wb)
    assert max(W) / min(W) - 1.0 < 0.02


@pytest.mark.parametrize("cls", [Quad4, Tri3, Tri6, Hex8, Tet4], ids=lambda c: c.__name__)
def test_modelo2_barra_en_compresion_objetiva(cls):
    """Trabajo sobre la banda desde su pico = G_c·A a ±2 % e igual entre
    mallas a ±2 %; pico 0.95·f_c·A y carga final ≈ 0 (rama completa)."""
    dim = len(cls.REFERENCE_NODES[0])
    _check_compression(cls, (5, 10, 20) if dim == 2 else (5, 10), 1.0)


@pytest.mark.parametrize("cls", [Hex20, Tet10], ids=lambda c: c.__name__)
def test_modelo2_barra_en_compresion_3d_cuadraticos(cls):
    """Sólidos 3D cuadráticos, sin factor de localización (α = 1): G_c·A.
    Hex27 da lo mismo (medido el 2026-11-18) pero tarda 80 s; queda fuera de
    la suite por tiempo, cubierto por Hex20 (misma h y el mismo α)."""
    _check_compression(cls, (5, 10), 1.0)


@pytest.mark.parametrize("cls", [Quad8, Quad9], ids=lambda c: c.__name__)
def test_modelo2_barra_en_compresion_cuadraticos_con_factor_disipan_1_sobre_alfa(cls):
    """Quad8/Quad9 siembran h = (13/18)·√A, y el material la usa también en
    κ_u: con el esfuerzo uniforme de la barra el elemento se ablanda entero y
    la banda entrega G_c·A/α = 1.385·G_c·A, objetivo entre mallas. Es el
    mismo comportamiento que en tracción (arriba), documentado sin decidir
    para la compresión (verificación del 2026-11-18, H1)."""
    _check_compression(cls, (5, 10, 20), 18.0 / 13.0, steps=1500, u_max=1.5)
#
# Flexión en cuatro puntos, H = 100, espesor 50, malla estructurada de celdas
# de 10 mm con un número impar de columnas (41): la entalla es la celda
# central (ancho 10, profundidad 50) y la banda sube por la columna central,
# alineada con la malla. Cargas por desplazamiento impuesto en dos nodos de
# la cara superior a 55 mm del eje; IMPL-EX con 1000 pasos hasta δ = 2 mm
# (carga final ≈ 1 % del pico). Concreto de Petersson (E 30 000, ν 0.2,
# f_t 3.33, G_f 0.124). Medido el 2026-11-12 (D/G_f·A en δ = 2 mm):
# Quad4 1.015, Quad8 sin α 0.773, con α 1.029; Tri3 1.038, Tri6 sin α 0.790,
# con ⅔ 1.168 (la fracción de Tri6 varía de un triángulo a otro: sin factor).

BEAM_E, BEAM_NU, BEAM_FT, BEAM_GF = 30000.0, 0.2, 3.33, 0.124
BEAM_H, BEAM_T, BEAM_DX = 100.0, 50.0, 10.0
_beams: dict = {}


def notched_beam(cls):
    """(curva [(δ, P)], D/(G_f·A_ligamento)) de la viga entallada; en caché."""
    if cls in _beams:
        return _beams[cls]
    dx = BEAM_DX
    nx, ny = 41, int(round(BEAM_H / dx))
    Lb = nx * dx
    xc = 0.5 * Lb
    mat = IsotropicDamageCrackBand2D(BEAM_E, BEAM_NU, BEAM_FT, BEAM_GF, integration="impl-ex")
    kw = {"diagonal": "right"} if cls.REFERENCE_CELL == "tri" else {}
    mesh = structured_mesh(cls, (nx, ny), size=(Lb, BEAM_H), material=mat, thickness=BEAM_T,
                           exclude=lambda x: abs(x[0] - xc) < 0.5 * dx and x[1] < 0.5 * BEAM_H,
                           **kw)
    dom = mesh.domain
    pts = {tuple(np.round(nd.coordinates[:2], 6)): nd for nd in dom.nodes.values()}
    pts[(0.0, 0.0)].fix_dof("ux", 0.0)
    pts[(0.0, 0.0)].fix_dof("uy", 0.0)
    pts[(round(Lb, 6), 0.0)].fix_dof("uy", 0.0)
    loads = [pts[(round(xc - 5.5 * dx, 6), BEAM_H)], pts[(round(xc + 5.5 * dx, 6), BEAM_H)]]
    for nd in loads:
        nd.fix_dof("uy", -2.0)
    dom.generate_equation_numbers(verbose=False)
    asm = Assembler(dom)
    dofs = [nd.dofs["uy"] for nd in loads]
    curve = [(0.0, 0.0)]

    def cb(step, U, lam):
        _, F = asm.assemble_non_linear_system(U)
        curve.append((-U[dofs[0]], -F[dofs].sum()))

    level = logging.getLogger("solidum").level
    logging.getLogger("solidum").setLevel(logging.ERROR)
    try:
        NonlinearSolver(asm, num_steps=1000, freeze_tangent_after_iter=1).solve(
            np.zeros(dom.total_dofs), step_callback=cb)
    finally:
        logging.getLogger("solidum").setLevel(level)
    c = np.array(curve)
    D = np.trapezoid(c[:, 1], c[:, 0]) - 0.5 * c[-1, 1] * c[-1, 0]
    _beams[cls] = (c, D / (BEAM_GF * 0.5 * BEAM_H * BEAM_T))
    return _beams[cls]


def test_viga_entallada_lineales_disipan_G_f():
    for cls in (Quad4, Tri3):
        curve, ratio = notched_beam(cls)
        assert ratio == pytest.approx(1.0, abs=0.06), cls.__name__
        assert curve[-1, 1] < 0.02 * curve[:, 1].max()      # ensayo casi completo


def test_viga_entallada_quad8_con_alfa_disipa_como_quad4():
    """α = 13/18: Quad8 disipa lo mismo que Quad4 (sin α, un 24 % menos) y su
    pico queda dentro de un 8 %."""
    c4, r4 = notched_beam(Quad4)
    c8, r8 = notched_beam(Quad8)
    assert r8 == pytest.approx(r4, rel=0.05)
    assert c8[:, 1].max() == pytest.approx(c4[:, 1].max(), rel=0.08)


def test_viga_entallada_tri6_sin_factor_disipa_de_menos():
    """Tri6 sin factor (decisión del usuario 2026-11-12): disipa menos que Tri3
    (0.76 de su energía), del lado seguro; con ⅔ disiparía de más (1.13)."""
    _, r3 = notched_beam(Tri3)
    _, r6 = notched_beam(Tri6)
    assert 0.65 < r6 / r3 < 0.9
