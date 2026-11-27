"""Apoyos y desplazamientos impuestos en las etapas del YAML (ADR 0026 §7) e
historia de un análisis por etapas (ADR 0036, 0038).

Invariantes que protegen estas pruebas:
- El valor que declara una etapa es el **total** al final de la etapa
  (decisión del usuario, 2026-11-27): crece en línea recta con λ desde el que
  había al empezarla. Un grado de libertad libre que una etapa prescribe
  parte de su valor actual, sin salto.
- Lo que impone una etapa sigue en las siguientes mientras ninguna lo cambie.
- Con etapas, la historia lleva la columna ``stage`` y registra reacciones y
  energías: la reacción con la carga del paso, ``F_c + λ·F_ref``; el trabajo
  externo se acumula de una etapa a otra.
- Las fuerzas internas de un marco con carga distribuida descuentan, en cada
  etapa, la carga equivalente constante de las anteriores más λ veces la de
  la etapa.
- Errores claros: apoyos dentro de ``loads``, en una etapa en el tiempo o de
  perturbación lineal (decisiones físicas pendientes), nodo, grupo o DOF
  inexistente (al leer), reacción de un grado de libertad que nunca se
  prescribe, columna llamada ``stage``.

Valores exactos con parámetros ≠ 1 (Reglas §6), *calculados*:
- Fricción: dos bloques que se mueven como sólidos rígidos (todos sus nodos
  prescritos) con una interfaz ``ContactPenalty`` (K_n = K_t = 1000, μ = 0.5,
  espesor 3, A = 6). Penetración 0.01 → p = 10, μpA = 30, pA = 60. Ciclo
  tangencial 0 → 0.02 → 0 con el cambio adherencia/deslizamiento en los
  bordes de paso (u = μp/K_t = 0.005): W = 1.125, E = 0.375 (normal 0.3 y
  tangencial 0.075), D = μpA·(0.015 + 0.01) = 0.75.
- Barra elástica k = EA/L = 300: carga P = 60 en la etapa 1, desplazamiento
  impuesto 0.5 en la 2 y 0 en la 3; R = k·u − P, W = E = ½·k·u², D = 0.
- Voladizo de marcos con peso propio q en la etapa 1 y carga de punta P en la
  2: M en el empotramiento −q·L²/2 − λ·P·L.
"""
from __future__ import annotations

import csv
import textwrap

import numpy as np
import pytest

import solidum
from solidum.utils.output_writer import validate_output_block
from solidum.utils.yaml_parser import YamlValidationError


def _write(tmp_path, texto, name="m.yaml"):
    ruta = tmp_path / name
    ruta.write_text(textwrap.dedent(texto), encoding="utf-8")
    return ruta


def _history(ruta):
    """Columnas del CSV de la historia; una celda vacía («no aplica») es NaN."""
    with open(ruta.with_name(ruta.stem + "_historia.csv"), newline="", encoding="utf-8") as f:
        filas = list(csv.reader(f))
    datos = np.array([[float(x) if x != "" else np.nan for x in fila] for fila in filas[1:]])
    return {c: datos[:, k] for k, c in enumerate(filas[0])}


# ---------------------------------------------------------------------------
# Fricción en tres etapas: cierre, deslizamiento y retorno
# ---------------------------------------------------------------------------

K, MU, T_B, GAP, PEN = 1000.0, 0.5, 3.0, 0.1, 0.01
AREA = 2.0 * T_B
P_N = K * PEN                  # presión de contacto
T_MAX = MU * P_N               # tracción tangencial de deslizamiento

FRICCION = f"""\
materials:
  - {{id: 1, type: Elastic2D, E: 1000.0, nu: 0.0}}
interface_materials:
  - {{id: 1, type: ContactPenalty, K_n: {K}, K_t: {K}, mu: {MU}}}
mesh_element: Quad4
mesh_thickness: {T_B}
mesh:
  structured:
    blocks:
      - {{name: abajo, material: 1, x: {{length: 2, n: 2}}, y: {{length: 1, n: 1}}}}
      - {{name: arriba, material: 1, origin: [0, {1 + GAP}], x: {{length: 2, n: 2}},
         y: {{length: 1, n: 1}}}}
    groups:
      Abajo: "y < {1 + GAP / 2}"
      Arriba: "y > {1 + GAP / 2}"
interfaces:
  - {{name: junta, pair: [abajo.ymax, arriba.ymin], material: 1}}
boundary_conditions_by_group:
  Abajo: {{ux: 0.0, uy: 0.0}}
  Arriba: {{ux: 0.0, uy: 0.0}}
stages:
  - name: acerca
    boundary_conditions_by_group:
      Arriba: {{uy: {-GAP}}}
    solver: {{type: NonlinearSolver, num_steps: 1, adaptive: false}}
  - name: aprieta
    boundary_conditions_by_group:
      Arriba: {{uy: {-(GAP + PEN)}}}
    solver: {{type: NonlinearSolver, num_steps: 1, adaptive: false}}
  - name: desliza
    boundary_conditions_by_group:
      Arriba: {{ux: 0.02}}
    solver: {{type: NonlinearSolver, num_steps: 4, adaptive: false}}
  - name: vuelve
    boundary_conditions_by_group:
      Arriba: {{ux: 0.0}}
    solver: {{type: NonlinearSolver, num_steps: 4, adaptive: false}}
output:
  log_level: WARNING
  history:
    columns:
      - {{name: u, displacement: {{group: Arriba, dof: ux}}}}
      - {{name: v, displacement: {{group: Arriba, dof: uy}}}}
      - {{name: Rx, reaction: {{group: Arriba, dof: ux}}}}
      - {{name: Ry, reaction: {{group: Arriba, dof: uy}}}}
      - {{name: W, energy: external_work}}
      - {{name: E, energy: elastic}}
      - {{name: D, energy: dissipated}}
"""


@pytest.fixture(scope="module")
def friccion(tmp_path_factory):
    ruta = _write(tmp_path_factory.mktemp("friccion"), FRICCION)
    res = solidum.run_yaml(ruta, write_output=True)
    return res, _history(ruta)


def test_la_historia_con_etapas_lleva_el_numero_de_etapa(friccion):
    _, h = friccion
    np.testing.assert_array_equal(h["stage"], [1, 2, 3, 3, 3, 3, 4, 4, 4, 4])
    np.testing.assert_array_equal(h["step"], np.arange(1, 11))
    np.testing.assert_allclose(h["factor"], [1, 1, .25, .5, .75, 1, .25, .5, .75, 1], rtol=0,
                               atol=1e-15)


def test_el_valor_de_una_etapa_es_total_y_crece_desde_el_anterior(friccion):
    """``ux: 0.02`` y después ``ux: 0.0``: la ida y la vuelta en línea recta
    con λ desde el valor con que empieza cada etapa."""
    _, h = friccion
    lam = h["factor"]
    st = h["stage"]
    np.testing.assert_allclose(h["u"][st == 3], 0.02 * lam[st == 3], rtol=0, atol=1e-15)
    np.testing.assert_allclose(h["u"][st == 4], 0.02 * (1 - lam[st == 4]), rtol=0, atol=1e-15)


def test_lo_que_impone_una_etapa_sigue_en_las_siguientes(friccion):
    """``uy`` lo cambian sólo las etapas 1 y 2; en la 3 y la 4 no se declara
    y se queda; ``ux`` sigue en 0 durante las dos primeras."""
    res, h = friccion
    st = h["stage"]
    np.testing.assert_allclose(h["v"][st >= 2], -(GAP + PEN), rtol=0, atol=1e-15)
    np.testing.assert_allclose(h["u"][st <= 2], 0.0, rtol=0, atol=0)
    assert res.names == ("acerca", "aprieta", "desliza", "vuelve")


def test_reacciones_de_contacto_y_friccion_exactas(friccion):
    _, h = friccion
    st = h["stage"]
    np.testing.assert_allclose(h["Ry"][st >= 2], -P_N * AREA, rtol=1e-12)
    # Adherencia hasta μp/K_t = 0.005 (borde del primer paso) y después
    # deslizamiento con μpA; a la vuelta, descarga elástica y deslizamiento
    # inverso.
    exacto = np.array([0, 0, 1, 1, 1, 1, 0, -1, -1, -1]) * T_MAX * AREA
    np.testing.assert_allclose(h["Rx"], exacto, rtol=0, atol=1e-10 * T_MAX * AREA)


def test_energias_exactas_a_traves_de_las_etapas(friccion):
    _, h = friccion
    W = 0.5 * P_N * PEN * AREA                                 # cierre
    E_n = W
    W_pasos = [W]
    for du, t0, t1 in ((0.005, 0, T_MAX), (0.005, T_MAX, T_MAX), (0.005, T_MAX, T_MAX),
                       (0.005, T_MAX, T_MAX), (-0.005, T_MAX, 0), (-0.005, 0, -T_MAX),
                       (-0.005, -T_MAX, -T_MAX), (-0.005, -T_MAX, -T_MAX)):
        W += 0.5 * (t0 + t1) * du * AREA
        W_pasos.append(W)
    E_t = 0.5 * K * (T_MAX / K) ** 2 * AREA
    np.testing.assert_allclose(h["W"], [0.0] + W_pasos, rtol=0, atol=1e-12)
    np.testing.assert_allclose(h["W"][-1], 1.125, rtol=1e-12)
    np.testing.assert_allclose(h["E"][-1], E_n + E_t, rtol=1e-12)
    np.testing.assert_allclose(h["D"][-1], T_MAX * AREA * (0.015 + 0.01), rtol=1e-12)
    np.testing.assert_allclose(h["D"][-1], 0.75, rtol=1e-12)
    assert np.all(np.diff(h["D"]) >= -1e-12)


# ---------------------------------------------------------------------------
# Barra elástica: carga, desplazamiento impuesto sobre un gdl libre y vuelta
# ---------------------------------------------------------------------------

E_B, A_B, L_B, P_B = 200.0, 3.0, 2.0, 60.0
K_B = E_B * A_B / L_B
U1, U2 = P_B / K_B, 0.5

BARRA = f"""\
nodes:
  - {{id: 1, coords: [0.0, 0.0]}}
  - {{id: 2, coords: [{L_B}, 0.0]}}
materials:
  - {{id: 1, type: Elastic1D, E: {E_B}}}
elements:
  - {{id: 1, type: Truss2D, material: 1, nodes: [1, 2], A: {A_B}}}
boundary_conditions_by_node:
  - {{node_id: 1, ux: 0.0, uy: 0.0}}
  - {{node_id: 2, uy: 0.0}}
"""

BARRA_ETAPAS = BARRA + f"""\
stages:
  - name: carga
    loads:
      point_loads_by_node:
        - {{node_id: 2, ux: {P_B}}}
    solver: {{type: NonlinearSolver, num_steps: 2, adaptive: false}}
  - name: impone
    boundary_conditions_by_node:
      - {{node_id: 2, ux: {U2}}}
    solver: {{type: NonlinearSolver, num_steps: 3, adaptive: false}}
  - name: vuelve
    boundary_conditions_by_node:
      - {{node_id: 2, ux: 0.0}}
    solver: {{type: NonlinearSolver, num_steps: 2, adaptive: false}}
"""

HISTORIA_BARRA = """\
output:
  log_level: WARNING
  history:
    columns:
      - {name: u, displacement: {node: 2, dof: ux}}
      - {name: R, reaction: {node: 2, dof: ux}}
      - {name: W, energy: external_work}
      - {name: E, energy: elastic}
      - {name: D, energy: dissipated}
"""


@pytest.fixture(scope="module")
def barra(tmp_path_factory):
    ruta = _write(tmp_path_factory.mktemp("barra"), BARRA_ETAPAS + HISTORIA_BARRA)
    solidum.run_yaml(ruta, write_output=True)
    return _history(ruta)


def test_prescribir_un_gdl_libre_parte_de_su_valor_actual(barra):
    h = barra
    lam, st = h["factor"], h["stage"]
    np.testing.assert_allclose(h["u"][st == 1], U1 * lam[st == 1], rtol=1e-12)
    np.testing.assert_allclose(h["u"][st == 2], U1 + lam[st == 2] * (U2 - U1), rtol=1e-12)
    np.testing.assert_allclose(h["u"][st == 3], U2 * (1 - lam[st == 3]), rtol=0, atol=1e-14)


def test_reaccion_con_la_carga_constante_de_la_etapa(barra):
    """En la etapa 1 el gdl está libre: no hay apoyo y la celda queda vacía
    (no el residuo del equilibrio); en la 2 y la 3 es ``k·u − P``, con la
    carga de la etapa 1 aún aplicada."""
    h = barra
    st = h["stage"]
    assert np.all(np.isnan(h["R"][st == 1]))
    np.testing.assert_allclose(h["R"][st >= 2], K_B * h["u"][st >= 2] - P_B, rtol=0,
                               atol=1e-9 * P_B)
    np.testing.assert_allclose(h["R"][-1], -P_B, rtol=1e-12)


def test_carga_y_descarga_elastica_no_disipa(barra):
    h = barra
    np.testing.assert_allclose(h["E"], 0.5 * K_B * h["u"] ** 2, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(h["W"], h["E"], rtol=0, atol=1e-12 * K_B * U2 ** 2)
    np.testing.assert_allclose(h["D"], 0.0, rtol=0, atol=1e-12 * K_B * U2 ** 2)
    assert abs(h["W"][-1]) < 1e-12 * K_B * U2 ** 2


def test_mismo_resultado_que_las_etapas_desde_python(tmp_path):
    """El YAML hace lo mismo que la API de Python con ``fix_dof`` entre etapas
    (ADR 0026 §2)."""
    from solidum import Assembler
    from solidum.math.solvers import NonlinearSolver
    ruta = _write(tmp_path, BARRA_ETAPAS)
    res = solidum.run_yaml(ruta)
    model = solidum.load_yaml(_write(tmp_path, BARRA + "solver: {type: LinearSolver}\n", "p.yaml"))
    dom = model.domain
    asm = Assembler(dom)
    n2 = dom.nodes[2]
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = P_B
    r1 = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=2, adaptive=False),
                     F_applied=F)
    n2.fix_dof("ux", U2)
    r2 = solidum.run(dom, assembler=asm, solver=NonlinearSolver(asm, num_steps=3, adaptive=False),
                     F_applied=np.zeros(dom.total_dofs), previous=r1)
    np.testing.assert_allclose(res["impone"].U, r2.U, rtol=0, atol=1e-15)
    np.testing.assert_allclose(res["impone"].R, r2.R, rtol=0, atol=1e-12)


def test_cada_solve_empieza_con_los_apoyos_del_nivel_raiz(tmp_path):
    """Las etapas cambian los apoyos del dominio; un segundo ``solve()`` del
    mismo modelo empieza otra vez con los del nivel raíz (el nodo 2 libre en
    la etapa 1) y da lo mismo. Al terminar, el dominio tiene los de la última
    etapa, los del estado final."""
    model = solidum.load_yaml(_write(tmp_path, BARRA_ETAPAS))
    n2 = model.domain.nodes[2]
    assert "ux" not in n2.boundary_conditions
    r1 = model.solve()
    assert n2.boundary_conditions["ux"] == 0.0
    seen = []
    r2 = model.solve(step_callback=lambda k, U, lam: seen.append(U[n2.dofs["ux"]]))
    np.testing.assert_allclose(seen[:2], [U1 / 2, U1], rtol=1e-12)
    for name in r1.names:
        np.testing.assert_allclose(r2[name].U, r1[name].U, rtol=0, atol=1e-14)


# ---------------------------------------------------------------------------
# Fuerzas internas de marco con carga distribuida en etapas
# ---------------------------------------------------------------------------

VOLADIZO = """\
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [2.0, 0.0]}
  - {id: 4, coords: [3.0, 0.0]}
  - {id: 5, coords: [4.0, 0.0]}
materials:
  - {id: 1, type: Elastic1D, E: 210.0e9, density: 0.0}
elements:
  - {id: 1, type: Frame2DEuler, material: 1, nodes: [1, 2], A: 1.0e-2, I: 8.33e-6}
  - {id: 2, type: Frame2DEuler, material: 1, nodes: [2, 3], A: 1.0e-2, I: 8.33e-6}
  - {id: 3, type: Frame2DEuler, material: 1, nodes: [3, 4], A: 1.0e-2, I: 8.33e-6}
  - {id: 4, type: Frame2DEuler, material: 1, nodes: [4, 5], A: 1.0e-2, I: 8.33e-6}
boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
stages:
  - name: peso
    loads:
      body_force: [0.0, -78.5e3]
    solver: {type: NonlinearSolver, num_steps: 2, adaptive: false}
  - name: punta
    loads:
      point_loads_by_node:
        - {node_id: 5, uy: -2000.0}
    solver: {type: NonlinearSolver, num_steps: 2, adaptive: false}
output:
  log_level: WARNING
  history:
    columns:
      - {name: M0, element: {id: 1, var: M, end: i}}
      - {name: V0, element: {id: 1, var: V, end: i}}
      - {name: M3, element: {id: 4, var: M, end: i}}
      - {name: W, energy: external_work}
      - {name: E, energy: elastic}
"""
Q, LV, P_V = 78.5e3 * 1.0e-2, 4.0, 2000.0


def test_fuerzas_internas_del_marco_por_etapas(tmp_path):
    ruta = _write(tmp_path, VOLADIZO, "voladizo.yaml")
    solidum.run_yaml(ruta, write_output=True)
    h = _history(ruta)
    lam, st = h["factor"], h["stage"]
    q = np.where(st == 1, lam, 1.0) * Q                # peso: crece en 1, constante en 2
    p = np.where(st == 2, lam, 0.0) * P_V              # punta: sólo en la 2
    np.testing.assert_allclose(h["M0"], -q * LV ** 2 / 2 - p * LV, rtol=1e-9)
    np.testing.assert_allclose(h["V0"], q * LV + p, rtol=1e-9)
    np.testing.assert_allclose(h["M3"], -q * 1.0 ** 2 / 2 - p * 1.0, rtol=1e-9)
    # Elástico lineal: W = E aunque la carga no sea proporcional.
    np.testing.assert_allclose(h["W"], h["E"], rtol=1e-9)


# ---------------------------------------------------------------------------
# Errores
# ---------------------------------------------------------------------------

def _errores_al_leer(tmp_path, texto) -> str:
    with pytest.raises((YamlValidationError, ValueError)) as exc:
        solidum.load_yaml(_write(tmp_path, texto))
    return str(exc.value)


def _barra_con(etapas: str) -> str:
    return BARRA + textwrap.dedent(etapas)


def test_apoyos_dentro_de_loads_es_un_error(tmp_path):
    msg = _errores_al_leer(tmp_path, _barra_con("""\
        stages:
          - name: mal
            loads:
              boundary_conditions_by_node:
                - {node_id: 2, ux: 0.1}
            solver: {type: NonlinearSolver}
        """))
    assert "no son cargas" in msg and "stages[0].loads" in msg


def test_nodo_inexistente_en_una_etapa_se_informa_al_leer(tmp_path):
    msg = _errores_al_leer(tmp_path, _barra_con("""\
        stages:
          - name: mal
            boundary_conditions_by_node:
              - {node_id: 9, ux: 0.1}
            solver: {type: NonlinearSolver}
        """))
    assert "stages[0].boundary_conditions_by_node" in msg and "node_id=9" in msg


def test_dof_inexistente_en_una_etapa_se_informa_al_leer(tmp_path):
    msg = _errores_al_leer(tmp_path, _barra_con("""\
        stages:
          - name: mal
            boundary_conditions_by_node:
              - {node_id: 2, rz: 0.0}
            solver: {type: NonlinearSolver}
        """))
    assert "stages[0].boundary_conditions" in msg and "'rz'" in msg


def test_grupo_inexistente_en_una_etapa_se_informa_al_leer(tmp_path):
    texto = FRICCION.replace("    boundary_conditions_by_group:\n      Arriba: {ux: 0.02}",
                             "    boundary_conditions_by_group:\n      Nada: {ux: 0.02}")
    assert texto != FRICCION
    msg = _errores_al_leer(tmp_path, texto)
    assert "stages[2].boundary_conditions_by_group" in msg and "'Nada'" in msg


@pytest.mark.parametrize("segunda,que", [
    ("solver: {type: LinearSolver}", "perturbación lineal"),
    ("solver: {type: NewmarkSolver, t_end: 0.1, dt: 0.01}", "'transient'"),
])
def test_apoyos_en_una_etapa_que_no_cambia_el_estado_estatico_es_un_error(tmp_path, segunda, que):
    """Se informa al leer, antes de correr las etapas anteriores."""
    texto = BARRA.replace("{id: 1, type: Elastic1D, E: 200.0}",
                          "{id: 1, type: Elastic1D, E: 200.0, density: 1.0}") + f"""\
stages:
  - name: carga
    loads:
      point_loads_by_node:
        - {{node_id: 2, ux: {P_B}}}
    solver: {{type: NonlinearSolver, num_steps: 2}}
  - name: segunda
    boundary_conditions_by_node:
      - {{node_id: 2, ux: 0.1}}
    {segunda}
"""
    msg = _errores_al_leer(tmp_path, texto)
    assert "stages[1]" in msg and que in msg and "decisión física pendiente" in msg


def test_un_linear_solver_como_primera_etapa_si_admite_apoyos(tmp_path):
    """Sin etapa anterior no es perturbación: es un análisis lineal más."""
    texto = BARRA + """\
stages:
  - name: impone
    boundary_conditions_by_node:
      - {node_id: 2, ux: 0.5}
    solver: {type: LinearSolver}
"""
    res = solidum.run_yaml(_write(tmp_path, texto))
    np.testing.assert_allclose(res.final.R[res.final.R != 0], [-K_B * 0.5, K_B * 0.5], rtol=1e-12)


def test_apoyo_de_etapa_sobre_un_esclavo_es_un_error_al_leer(tmp_path):
    texto = BARRA.replace("  - {node_id: 2, uy: 0.0}\n", "") + """\
linear_constraints:
  - {slave: {node: 2, dof: uy}, masters: [{node: 1, dof: uy}], coefficients: [1.0]}
stages:
  - name: mal
    boundary_conditions_by_node:
      - {node_id: 2, uy: 0.1}
    solver: {type: NonlinearSolver}
"""
    msg = _errores_al_leer(tmp_path, texto)
    assert "stages[0]" in msg and "esclavos" in msg and "(2, 'uy')" in msg


def test_fuerzas_de_marco_con_carga_distribuida_y_etapa_en_el_tiempo_es_un_error(tmp_path):
    """La carga equivalente de una etapa estática sigue en la etapa en el
    tiempo, donde la historia no la conoce: se rechaza al crear la salida, no
    se escribe un momento sin descontarla."""
    texto = VOLADIZO.replace("density: 0.0", "density: 7850.0").replace(
        "  - name: punta\n    loads:\n      point_loads_by_node:\n        - {node_id: 5, uy: -2000.0}\n"
        "    solver: {type: NonlinearSolver, num_steps: 2, adaptive: false}\n",
        "  - name: vibra\n    loads:\n      point_loads_by_node:\n        - {node_id: 5, uy: -2000.0}\n"
        "    solver: {type: NewmarkSolver, t_end: 0.01, dt: 0.001}\n").replace(
        "      - {name: W, energy: external_work}\n      - {name: E, energy: elastic}\n", "")
    assert "NewmarkSolver" in texto and "energy" not in texto
    with pytest.raises(ValueError) as exc:
        solidum.run_yaml(_write(tmp_path, texto, "voladizo.yaml"), write_output=True)
    assert "si ninguna es en el tiempo" in str(exc.value)


def test_reaccion_de_un_gdl_que_ninguna_etapa_prescribe_es_un_error(tmp_path):
    texto = BARRA + """\
stages:
  - name: carga
    loads:
      point_loads_by_node:
        - {node_id: 2, ux: 60.0}
    solver: {type: NonlinearSolver}
output:
  history:
    columns:
      - {name: R, reaction: {node: 2, dof: ux}}
"""
    with pytest.raises(ValueError) as exc:
        solidum.run_yaml(_write(tmp_path, texto), write_output=True)
    assert "en el nivel raíz o en alguna etapa" in str(exc.value)


def test_stage_es_un_nombre_reservado_de_la_historia():
    errores = validate_output_block(
        {"history": {"columns": [{"name": "stage", "displacement": {"node": 1, "dof": "ux"}}]}})
    assert any("'stage'" in e for e in errores)
