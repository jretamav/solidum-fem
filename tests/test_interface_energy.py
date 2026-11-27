"""Energías de la historia con elementos de interfaz (ADR 0038 y ADR 0043 §5).

La energía almacenada es ``E = Σ_e E_e``: ``½·u_e·f_e`` en los sólidos y la de
la ley (``Σ_g m_g·ψ_g``) en las interfaces. Con holgura inicial o con fricción
``½·u_e·f_e`` no es la energía almacenada: con el resorte unilateral del
arquitecto (holgura 1, cierre 2) daba una disipación de −0.5.

Modelo: un bloque Quad4 apoyado, a través de un ``Interface2D4``, sobre una
base fija; la cara superior del bloque baja y se desliza con desplazamientos
impuestos proporcionales. Con ν = 0 el bloque es un resorte en serie con la
interfaz, así que las energías tienen forma cerrada.
"""
import csv
import logging

import numpy as np
import pytest

from solidum import Assembler, ContactPenalty, Elastic2D, Interface2D4, NonlinearSolver, Quad4, run
from solidum.core.domain import Domain
from solidum.utils.output_writer import OutputWriter


@pytest.fixture(autouse=True)
def _quiet():
    logging.disable(logging.WARNING)
    yield
    logging.disable(logging.NOTSET)


def _leer(path):
    with open(path, newline="", encoding="utf-8") as f:
        filas = list(csv.reader(f))
    return {c: np.array([float(r[k]) for r in filas[1:]]) for k, c in enumerate(filas[0])}


def _modelo(law, gap, uy_top, ux_top, E_bloque=1.0e9):
    """Base fija (nodos 1, 2 en y = 0), interfaz con la cara A en y = 0 y la
    cara B a ``gap`` encima (nodos 3, 4), bloque Quad4 de 1×1 encima con la
    cara superior (nodos 5, 6) desplazada ``(ux_top, uy_top)``. Con ν = 0 el
    bloque se comprime uniformemente: un resorte de rigidez ``E_bloque`` en
    serie con la interfaz (``K_n·A``, A = 1)."""
    d = Domain()
    for nid, x in ((1, [0.0, 0.0]), (2, [1.0, 0.0]), (3, [1.0, gap]), (4, [0.0, gap]),
                   (5, [1.0, gap + 1.0]), (6, [0.0, gap + 1.0])):
        d.add_node(nid, x)
    n = d.nodes
    d.add_element(Interface2D4(1, [n[1], n[2], n[3], n[4]], law))
    bloque = Elastic2D(E=E_bloque, nu=0.0, hypothesis="plane_stress")
    d.add_element(Quad4(2, [n[4], n[3], n[5], n[6]], bloque, thickness=1.0))
    for nid in (1, 2):
        n[nid].fix_dof("ux")
        n[nid].fix_dof("uy")
    for nid in (5, 6):
        n[nid].fix_dof("ux", ux_top)
        n[nid].fix_dof("uy", uy_top)
    d.generate_equation_numbers()
    return d


def _historia(tmp_path, d, columnas, pasos):
    asm = Assembler(d)
    w = OutputWriter(d, asm, {"file_name": "contacto", "history": {"columns": columnas}},
                     tmp_path / "contacto.yaml", np.zeros(d.total_dofs), static=True)
    run(d, assembler=asm, solver=NonlinearSolver(asm, num_steps=pasos),
        F_applied=np.zeros(d.total_dofs), step_callback=w.step)
    return _leer(tmp_path / "contacto_historia.csv")


ENERGIAS = [{"name": "W", "energy": "external_work"}, {"name": "E", "energy": "elastic"},
            {"name": "D", "energy": "dissipated"}]


# Holgura g_0 = 1 y bajada 2 de la cara superior: la interfaz se cierra con la
# bajada 1 y el resto, δ = 1, se reparte entre la interfaz (K_n·A = 1) y el
# bloque (E·A/h = 1e3) en serie. Bloque de E = 1e3 y no más rígido: mientras la
# holgura está abierta el bloque baja sin fuerza, y con E = 1e9 el redondeo de
# su residuo (≈ 1e-9 N) queda en el piso absoluto del criterio de convergencia
# (1e-9·1 N sin fuerza de referencia, ADR 0007): converge o no según el paso
# (*medido* el 2026-11-26).
K_SERIE = 1.0 * 1.0e3 / (1.0 + 1.0e3)
E_ALMACENADA = 0.5 * K_SERIE * 1.0 ** 2


def _holgura():
    return _modelo(ContactPenalty(K_n=1.0), gap=1.0, uy_top=-2.0, ux_top=0.0, E_bloque=1.0e3)


def test_holgura_inicial_sin_disipacion_espuria(tmp_path):
    """Contraejemplo del arquitecto: holgura g_0 = 1, bajada 2, K_n·A = 1.
    Antes: E = ½·U·F_int ≈ 1.0 y D ≈ −0.5 (con el bloque rígido). Ahora
    E = ψ + energía del bloque = ½·k·δ² ≈ 0.5 y D = 0 (los pasos caen en el
    cierre: sin error del trapecio)."""
    h = _historia(tmp_path, _holgura(), ENERGIAS, pasos=4)
    assert h["E"][-1] == pytest.approx(E_ALMACENADA, rel=1e-9)
    assert h["W"][-1] == pytest.approx(E_ALMACENADA, rel=1e-9)
    assert abs(h["D"][-1]) < 1e-9 * E_ALMACENADA


@pytest.mark.parametrize("pasos", [3, 30])
def test_cierre_dentro_de_un_paso_disipacion_del_trapecio(tmp_path, pasos):
    """Con el cierre dentro de un paso, el trapecio sobrestima el trabajo de
    ese paso y aparece una disipación positiva que baja con el paso (nota del
    manual, decisión del usuario del 2026-11-26). Con 3 pasos (bajadas 2/3,
    4/3, 2) el cierre, en la bajada 1, cae a mitad del segundo: fuerza de 0 a
    k/3; trapecio ½·(k/3)·(2/3) = k/9, exacto ½·k·(1/3)² = k/18: disipación
    aparente k/18. Con 30 pasos el cierre cae justo en el paso 15: sin error."""
    h = _historia(tmp_path, _holgura(), ENERGIAS, pasos=pasos)
    assert h["E"][-1] == pytest.approx(E_ALMACENADA, rel=1e-9)   # no depende del paso
    esperada = K_SERIE / 18.0 if pasos == 3 else 0.0
    assert h["D"][-1] == pytest.approx(esperada, abs=1e-9 * E_ALMACENADA)


def test_friccion_trabajo_es_almacenada_mas_disipada(tmp_path):
    """Cierre c y deslizamiento δ = c proporcionales, μ = 0.3, K_t = K_n: la
    interfaz desliza desde el inicio (K_t·δ > μ·K_n·c). Forma cerrada por
    unidad de área: W = ½K_n c² + ½μK_n c δ, ψ = ½K_n c² + ½(μK_n c)²/K_t,
    D = ½μK_n c (δ − μK_n c/K_t). Con ½·u·f la disipación saldría 0."""
    Kn, Kt, mu, c, delta = 1.0e3, 1.0e3, 0.3, 1.0e-3, 1.0e-3
    d = _modelo(ContactPenalty(K_n=Kn, mu=mu, K_t=Kt), gap=0.0, uy_top=-c, ux_top=delta)
    h = _historia(tmp_path, d, ENERGIAS, pasos=5)
    W = 0.5 * Kn * c * c + 0.5 * mu * Kn * c * delta
    psi = 0.5 * Kn * c * c + 0.5 * (mu * Kn * c) ** 2 / Kt
    D = 0.5 * mu * Kn * c * (delta - mu * Kn * c / Kt)
    assert h["W"][-1] == pytest.approx(W, rel=1e-4)
    assert h["E"][-1] == pytest.approx(psi, rel=1e-4)
    assert h["D"][-1] == pytest.approx(D, rel=1e-4)
    assert np.all(np.diff(h["D"]) >= -1e-12 * W)     # la disipación no decrece


def test_friccion_con_parametros_no_unitarios_escala_con_el_area(tmp_path):
    """La misma forma cerrada con ancho L = 1.7, espesor e = 0.6, K_t = 0.4·K_n
    y μ = 0.35 (Reglas §6: con valores unitarios un factor de área olvidado es
    la identidad). Las tres energías escalan con el área A = L·e."""
    Kn, Kt, mu, c, delta, L, e = 1.0e3, 4.0e2, 0.35, 1.0e-3, 1.0e-3, 1.7, 0.6
    d = Domain()
    for nid, x in ((1, [0.0, 0.0]), (2, [L, 0.0]), (3, [L, 0.0]), (4, [0.0, 0.0]),
                   (5, [L, 1.0]), (6, [0.0, 1.0])):
        d.add_node(nid, x)
    n = d.nodes
    d.add_element(Interface2D4(1, [n[1], n[2], n[3], n[4]],
                               ContactPenalty(K_n=Kn, mu=mu, K_t=Kt), thickness=e))
    d.add_element(Quad4(2, [n[4], n[3], n[5], n[6]],
                        Elastic2D(E=1.0e9, nu=0.0, hypothesis="plane_stress"), thickness=e))
    for nid in (1, 2):
        n[nid].fix_dof("ux")
        n[nid].fix_dof("uy")
    for nid in (5, 6):
        n[nid].fix_dof("ux", delta)
        n[nid].fix_dof("uy", -c)
    d.generate_equation_numbers()
    assert Kt * delta > mu * Kn * c          # desliza desde el inicio
    h = _historia(tmp_path, d, ENERGIAS, pasos=5)
    A = L * e
    W = A * (0.5 * Kn * c * c + 0.5 * mu * Kn * c * delta)
    psi = A * (0.5 * Kn * c * c + 0.5 * (mu * Kn * c) ** 2 / Kt)
    D = A * 0.5 * mu * Kn * c * (delta - mu * Kn * c / Kt)
    assert h["W"][-1] == pytest.approx(W, rel=1e-4)
    assert h["E"][-1] == pytest.approx(psi, rel=1e-4)
    assert h["D"][-1] == pytest.approx(D, rel=1e-4)


@pytest.mark.parametrize("defecto, error", [("sin_nodal_forces", TypeError),
                                            ("desplazamientos_grandes", ValueError)])
def test_la_energia_se_comprueba_al_crear_la_columna(tmp_path, defecto, error):
    """Revisión del arquitecto del 2026-11-26: quien entrega su energía debe
    dar también ``nodal_forces``, y la ``G_e`` fija de las interfaces no vale
    con desplazamientos grandes; las dos cosas fallan al crear la columna, no
    en el primer paso."""
    d = _modelo(ContactPenalty(K_n=1.0), gap=0.0, uy_top=-1.0, ux_top=0.0)
    el = d.elements[1]
    if defecto == "sin_nodal_forces":
        el.nodal_forces = None
    else:
        el.GEOMETRICALLY_NONLINEAR = True
    asm = Assembler(d)
    with pytest.raises(error):
        OutputWriter(d, asm, {"file_name": "c", "history": {"columns": ENERGIAS}},
                     tmp_path / "c.yaml", np.zeros(d.total_dofs), static=True)


def test_energias_por_zona_con_interfaces(tmp_path):
    """La interfaz como zona: su energía almacenada es la de la ley y la
    disipación por fricción está en ella; con el bloque, la partición suma
    las energías del modelo."""
    Kn, mu = 1.0e3, 0.3
    d = _modelo(ContactPenalty(K_n=Kn, mu=mu, K_t=Kn), gap=0.0, uy_top=-1e-3, ux_top=1e-3)
    zonas = [{"name": f"{n}{z}", "energy": e, "elements": {"ids": [i]}}
             for z, i in (("i", 1), ("b", 2))
             for n, e in (("W", "external_work"), ("E", "elastic"), ("D", "dissipated"))]
    h = _historia(tmp_path, d, ENERGIAS + zonas, pasos=5)
    for k in ("W", "E", "D"):
        np.testing.assert_allclose(h[k + "i"] + h[k + "b"], h[k], rtol=1e-8, atol=1e-15)
    np.testing.assert_allclose(h["Di"], h["D"], rtol=1e-6)
    assert np.all(np.abs(h["Db"]) <= 1e-9 * h["W"])


def test_ley_sin_energia_conocida_es_un_error(tmp_path):
    """Una interfaz con una ley que no descarga por la secante ni define su
    energía: 'elastic' y 'dissipated' no se pueden registrar."""
    from solidum.core.interface_material import InterfaceMaterial

    class SinEnergia(InterfaceMaterial):
        GEOMETRY = "surface"
        STATE_SCHEMA = {"x": ()}

        def compute_traction(self, jump, state_vars=None):
            return 1e3 * np.asarray(jump), 1e3 * np.eye(len(jump)), {"x": 0.0}

    d = _modelo(SinEnergia(), gap=0.0, uy_top=-1e-3, ux_top=0.0)
    asm = Assembler(d)
    with pytest.raises(ValueError, match="SinEnergia"):
        OutputWriter(d, asm, {"file_name": "c", "history": {"columns": ENERGIAS}},
                     tmp_path / "c.yaml", np.zeros(d.total_dofs), static=True)
    # El trabajo externo sí.
    h = _historia(tmp_path, d, ENERGIAS[:1], pasos=2)
    assert h["W"][-1] > 0.0
