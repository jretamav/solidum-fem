"""Energías y magnitudes de elemento en la historia de la salida (ADR 0038).

1. Barra de la prueba de objetividad de la banda de grieta, por el camino del
   YAML (malla de tres bloques, la columna central debilitada al 95 %), Quad4
   de Δ = 10 mm, ν = 0: la energía disipada al separarse es G_f·A con ±1 %;
   antes del pico D = 0 y E = W; ``½·U·F_int`` coincide con ``½·P·u``; la
   media de σ_xx en la columna debilitada es P/A en cada paso y la de ε_xx
   en el bloque sano, P/(A·E); el daño máximo llega a la saturación.
2. Voladizo de cuatro Frame2DEuler con peso propio: M y V en el empotramiento
   son −q·L²/2 y +q·L (convención de viga, Reglas §5) en cada paso, con la
   carga equivalente descontada.
3. Errores: energía elástica con plasticidad, variable inexistente,
   componente inválida, energía en un análisis no estático.
"""
import csv
import textwrap

import numpy as np
import pytest

import solidum
from solidum.utils.output_writer import validate_output_block

E, FT, GF = 30000.0, 3.0, 0.1
L, H, T = 100.0, 20.0, 10.0
A = H * T


def _leer(path):
    with open(path, newline="", encoding="utf-8") as f:
        filas = list(csv.reader(f))
    cab = filas[0]
    datos = np.array(filas[1:], dtype=float)
    return {c: datos[:, k] for k, c in enumerate(cab)}


BARRA = f"""\
materials:
  - {{id: 1, type: IsotropicDamageCrackBand2D, E: {E}, nu: 0.0, f_t: {FT}, G_f: {GF}}}
  - {{id: 2, type: IsotropicDamageCrackBand2D, E: {E}, nu: 0.0, f_t: {0.95 * FT}, G_f: {GF}}}
mesh_element: Quad4
mesh_thickness: {T}
mesh:
  structured:
    blocks:
      - {{name: a, material: 1, x: {{length: 50, n: 5}}, y: {{length: {H}, n: 2}}}}
      - {{name: w, material: 2, origin: [50, 0], x: {{length: 10, n: 1}}, y: {{length: {H}, n: 2}}}}
      - {{name: b, material: 1, origin: [60, 0], x: {{length: 40, n: 4}}, y: {{length: {H}, n: 2}}}}
    groups:
      Origen: "near(x, 0) and near(y, 0)"
boundary_conditions_by_group:
  a.xmin: {{ux: 0.0}}
  Origen: {{uy: 0.0}}
  b.xmax: {{ux: 0.4}}
solver: {{type: NonlinearSolver, num_steps: 400}}
output:
  file_name: barra
  log_level: WARNING
  history:
    columns:
      - {{name: u, displacement: {{group: b.xmax, dof: ux}}}}
      - {{name: P, reaction: {{group: b.xmax, dof: ux}}}}
      - {{name: W, energy: external_work}}
      - {{name: Ee, energy: elastic}}
      - {{name: D, energy: dissipated}}
      - {{name: d_max, element: {{var: damage, reduce: max}}}}
      - {{name: sxx_w, element: {{where: "x > 50 and x < 60", var: stress, component: xx}}}}
      - {{name: exx_a, element: {{where: "x < 50", var: strain, component: xx}}}}
      - {{name: mises_w, element: {{where: "x > 50 and x < 60", var: stress, component: mises}}}}
      - {{name: Wa, energy: external_work, elements: {{where: "x < 50"}}}}
      - {{name: Ww, energy: external_work, elements: {{where: "x > 50 and x < 60"}}}}
      - {{name: Wb, energy: external_work, elements: {{where: "x > 60"}}}}
      - {{name: Ea, energy: elastic, elements: {{where: "x < 50"}}}}
      - {{name: Ewb, energy: elastic, elements: {{where: "x > 50"}}}}
      - {{name: Da, energy: dissipated, elements: {{where: "x < 50"}}}}
      - {{name: Dw, energy: dissipated, elements: {{where: "x > 50 and x < 60"}}}}
      - {{name: Dwb, energy: dissipated, elements: {{where: "x > 50"}}}}
"""


@pytest.fixture(scope="module")
def barra(tmp_path_factory):
    d = tmp_path_factory.mktemp("barra")
    ruta = d / "barra.yaml"
    ruta.write_text(BARRA, encoding="utf-8")
    solidum.run_yaml(ruta, write_output=True)
    return _leer(d / "barra_historia.csv")


def test_energia_disipada_es_gf_por_area(barra):
    assert barra["D"][-1] / (GF * A) == pytest.approx(1.0, abs=0.01)
    np.testing.assert_allclose(barra["D"], barra["W"] - barra["Ee"], atol=1e-9)


def test_antes_del_pico_no_hay_disipacion(barra):
    pico = int(np.argmax(barra["P"]))
    pre = slice(0, max(pico - 5, 1))
    np.testing.assert_allclose(barra["Ee"][pre], barra["W"][pre], rtol=1e-8)
    assert np.max(np.abs(barra["D"][pre])) < 1e-8 * barra["W"][pico]
    assert np.all(barra["d_max"][pre] == 0.0)


def test_energia_elastica_es_medio_p_por_u(barra):
    np.testing.assert_allclose(barra["Ee"], 0.5 * barra["P"] * barra["u"],
                               rtol=1e-6, atol=1e-9 * GF * A)
    # El trabajo es el área bajo P–u (trapecio entre pasos).
    W = np.concatenate([[0.0], np.cumsum(0.5 * (barra["P"][1:] + barra["P"][:-1])
                                         * np.diff(barra["u"]))]) + 0.5 * barra["P"][0] * barra["u"][0]
    np.testing.assert_allclose(barra["W"], W, rtol=1e-6, atol=1e-9 * GF * A)


def test_esfuerzo_y_deformacion_de_elemento(barra):
    np.testing.assert_allclose(barra["sxx_w"] * A, barra["P"], rtol=1e-6, atol=1e-9 * FT * A)
    np.testing.assert_allclose(barra["exx_a"], barra["P"] / (A * E), rtol=1e-6,
                               atol=1e-12 * FT / E)
    # Tracción uniaxial (ν = 0, esfuerzo plano): von Mises = |σ_xx|.
    np.testing.assert_allclose(barra["mises_w"], np.abs(barra["sxx_w"]), rtol=1e-9, atol=1e-12)
    assert barra["d_max"][-1] > 0.99


def test_las_energias_por_zona_suman_las_del_modelo(barra):
    """Con una partición del modelo (espesor 10 ≠ 1), la suma de las energías
    de las zonas es la del modelo completo en cada paso (hasta el redondeo
    del CSV, 10 cifras)."""
    tol = dict(rtol=1e-8, atol=1e-8 * GF * A)
    np.testing.assert_allclose(barra["Wa"] + barra["Ww"] + barra["Wb"], barra["W"], **tol)
    np.testing.assert_allclose(barra["Ea"] + barra["Ewb"], barra["Ee"], **tol)
    np.testing.assert_allclose(barra["Da"] + barra["Dwb"], barra["D"], **tol)


def test_la_disipacion_esta_en_la_columna_que_se_agrieta(barra):
    """La columna debilitada disipa G_f·A; los bloques que descargan
    elásticamente, nada."""
    assert barra["Dw"][-1] / (GF * A) == pytest.approx(1.0, abs=0.01)
    assert np.max(np.abs(barra["Da"])) < 1e-6 * GF * A
    np.testing.assert_allclose(barra["Dwb"], barra["Dw"], rtol=1e-6, atol=1e-6 * GF * A)


def test_energia_por_zona_en_un_analisis_lineal(tmp_path):
    """El análisis lineal no confirma estado: la zona lee el de prueba."""
    ruta = _placa(tmp_path, ["{name: W, energy: external_work}",
                             "{name: W1, energy: external_work, elements: {where: \"x < 0.5\"}}",
                             "{name: W2, energy: external_work, elements: {where: \"x > 0.5\"}}",
                             "{name: E2, energy: elastic, elements: {where: \"x > 0.5\"}}"],
                  material="{id: 1, type: Elastic2D, E: 200000.0, nu: 0.3}")
    ruta.write_text(ruta.read_text(encoding="utf-8").replace(
        "solver: {type: NonlinearSolver, num_steps: 2}", "solver: {type: LinearSolver}"),
        encoding="utf-8")
    solidum.run_yaml(ruta, write_output=True)
    h = _leer(tmp_path / "placa_historia.csv")
    assert h["W"][0] > 0.0
    np.testing.assert_allclose(h["W1"] + h["W2"], h["W"], rtol=1e-8)
    np.testing.assert_allclose(h["E2"], h["W2"], rtol=1e-8)    # elástico: E = W


# ----------------------------------------------------------------------------
# Marco con peso propio
# ----------------------------------------------------------------------------

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
body_force: [0.0, -78.5e3]
"""
Q, LV = 78.5e3 * 1.0e-2, 4.0


@pytest.mark.parametrize("solver,pasos", [("solver: {type: LinearSolver}\n", 1),
                                          ("solver: {type: NonlinearSolver, num_steps: 4}\n", 4)])
def test_fuerzas_internas_del_marco(tmp_path, solver, pasos):
    ruta = tmp_path / "voladizo.yaml"
    ruta.write_text(VOLADIZO + solver + textwrap.dedent("""\
        output:
          file_name: voladizo
          history:
            columns:
              - {name: M0, element: {id: 1, var: M, end: i}}
              - {name: V0, element: {id: 1, var: V, end: i}}
              - {name: M3, element: {id: 4, var: M, end: i}}
              - {name: Mfin, element: {id: 4, var: M, end: j}}
              - {name: W, energy: external_work}
              - {name: Ee, energy: elastic}
        """), encoding="utf-8")
    solidum.run_yaml(ruta, write_output=True)
    h = _leer(tmp_path / "voladizo_historia.csv")
    assert len(h["step"]) == pasos
    f = h["factor"]
    np.testing.assert_allclose(h["M0"], -f * Q * LV ** 2 / 2, rtol=1e-9)
    np.testing.assert_allclose(h["V0"], f * Q * LV, rtol=1e-9)
    np.testing.assert_allclose(h["M3"], -f * Q * 1.0 ** 2 / 2, rtol=1e-9)
    np.testing.assert_allclose(h["Mfin"], 0.0, atol=1e-9 * Q * LV ** 2)
    # Elástico lineal con carga proporcional: W = E (trapecio exacto).
    np.testing.assert_allclose(h["W"], h["Ee"], rtol=1e-9)


def test_energia_por_zona_con_marcos_es_un_error(tmp_path):
    ruta = tmp_path / "voladizo.yaml"
    ruta.write_text(VOLADIZO + "solver: {type: LinearSolver}\n" + textwrap.dedent("""\
        output:
          file_name: voladizo
          history:
            columns:
              - {name: W, energy: external_work, elements: {ids: [1, 2]}}
        """), encoding="utf-8")
    with pytest.raises(ValueError, match="sólidos isoparamétricos"):
        solidum.run_yaml(ruta, write_output=True)


# ----------------------------------------------------------------------------
# SECANT_UNLOADING de las leyes cohesivas (los materiales de volumen, en el
# barrido de contratos)
# ----------------------------------------------------------------------------

@pytest.mark.parametrize("nombre", ["CohesiveDamageIsotropic", "CohesiveDamageMixedMode"])
def test_cohesivas_descargan_por_la_secante(nombre):
    from solidum.registry import InterfaceMaterialRegistry
    sigma_t0, G_f, K_e = 3.0e6, 100.0, 1.0e15
    extra = {"beta": 1.5} if nombre == "CohesiveDamageMixedMode" else {}
    ley = InterfaceMaterialRegistry.create(nombre, sigma_t0=sigma_t0, G_f=G_f, K_e=K_e,
                                           softening="linear", **extra)
    assert ley.SECANT_UNLOADING
    salto = np.array([0.3 * G_f / sigma_t0, 0.1 * G_f / sigma_t0])   # en ablandamiento
    t1, _, st1 = ley.compute_traction(salto)
    assert 0.0 < np.linalg.norm(t1) < sigma_t0
    for frac in (0.5, 0.0):
        t, _, _ = ley.compute_traction(frac * salto, st1)
        np.testing.assert_allclose(t, frac * np.asarray(t1), rtol=1e-9, atol=1e-6)


# ----------------------------------------------------------------------------
# Errores
# ----------------------------------------------------------------------------

PLACA = """\
materials:
  - {id: 1, type: VonMises2D, E: 200000.0, nu: 0.3, sigma_y: 250.0}
mesh_element: Quad4
mesh:
  structured: {x: {length: 1, n: 2}, y: {length: 1, n: 2}}
boundary_conditions_by_group:
  xmin: {ux: 0.0, uy: 0.0}
  xmax: {ux: 0.003}
solver: {type: NonlinearSolver, num_steps: 2}
"""


def _placa(tmp_path, columnas, material=None):
    ruta = tmp_path / "placa.yaml"
    texto = PLACA if material is None else PLACA.replace(
        "{id: 1, type: VonMises2D, E: 200000.0, nu: 0.3, sigma_y: 250.0}", material)
    ruta.write_text(texto + "output:\n  file_name: placa\n  history:\n    columns:\n"
                    + "".join(f"      - {c}\n" for c in columnas), encoding="utf-8")
    return ruta


def test_energia_elastica_con_plasticidad_es_un_error(tmp_path):
    with pytest.raises(ValueError, match="VonMises2D"):
        solidum.run_yaml(_placa(tmp_path, ["{name: E, energy: elastic}"]), write_output=True)


def test_trabajo_externo_y_deformacion_plastica_si(tmp_path):
    ruta = _placa(tmp_path, ["{name: W, energy: external_work}",
                             "{name: a, element: {var: alpha, reduce: max}}",
                             "{name: s, element: {id: 1, var: stress, component: max_principal}}"])
    solidum.run_yaml(ruta, write_output=True)
    h = _leer(tmp_path / "placa_historia.csv")
    assert h["W"][-1] > h["W"][0] > 0.0
    assert h["a"][-1] > 0.0                      # fluyó
    assert np.all(np.isfinite(h["s"]))


def test_energia_disipada_por_zona_con_plasticidad_es_un_error(tmp_path):
    with pytest.raises(ValueError, match="VonMises2D"):
        solidum.run_yaml(_placa(tmp_path, ["{name: D, energy: dissipated, elements: {id: 1}}"]),
                         write_output=True)


@pytest.mark.parametrize("columna,mensaje", [
    ("{name: x, element: {var: damage}}", "variable interna escalar 'damage'"),
    ("{name: x, element: {id: 1, var: eps_p}}", "variable interna escalar 'eps_p'"),
    ("{name: x, element: {id: 1, var: N, end: i}}", "fuerza interna 'N'"),
    ("{name: x, element: {id: 99, var: stress, component: xx}}", "inexistentes"),
    ("{name: x, element: {id: 1, var: stress, component: yz}}", "no existe con 3 componentes"),
    ("{name: x, element: {where: \"x > 5\", var: stress, component: xx}}", "ningún elemento"),
])
def test_errores_al_resolver(tmp_path, columna, mensaje):
    with pytest.raises(ValueError, match=mensaje):
        solidum.run_yaml(_placa(tmp_path, [columna]), write_output=True)


@pytest.mark.parametrize("bloque,mensaje", [
    ({"history": {"columns": [{"name": "E", "energy": "cinetica"}]}}, "no admitida"),
    ({"history": {"columns": [{"name": "s", "element": {"var": "stress"}}]}}, "component"),
    ({"history": {"columns": [{"name": "s", "element": {"var": "M"}}]}}, "extremo"),
    ({"history": {"columns": [{"name": "s", "element": {"var": "damage", "id": 1, "ids": [2]}}]}},
     "como mucho uno"),
    ({"history": {"columns": [{"name": "s", "element": {"var": "damage", "reduce": "sum"}}]}},
     "'mean', 'min' o 'max'"),
    ({"history": {"columns": [{"name": "s", "element": {"var": "damage", "where": "x == 1"}}]}},
     "near"),
    ({"history": {"columns": [{"name": "s", "element": {"var": "damage", "component": "xx"}}]}},
     "sólo con var: stress, strain, jump o traction"),
    ({"history": {"columns": [{"name": "D", "energy": "dissipated", "elements": "banda"}]}},
     "un diccionario con uno de"),
    ({"history": {"columns": [{"name": "D", "energy": "dissipated",
                               "elements": {"id": 1, "group": "g"}}]}}, "exactamente uno de"),
    ({"history": {"columns": [{"name": "D", "energy": "dissipated",
                               "elements": {"where": "x == 1"}}]}}, "near"),
    ({"history": {"columns": [{"name": "u", "displacement": {"node": 1, "dof": "ux"},
                               "elements": {"id": 1}}]}}, "no admitidas"),
])
def test_validacion(bloque, mensaje):
    errores = validate_output_block(bloque)
    assert any(mensaje in e for e in errores), errores


def test_energia_en_un_transitorio_es_un_error(tmp_path):
    ruta = tmp_path / "din.yaml"
    ruta.write_text(textwrap.dedent("""\
        nodes:
          - {id: 1, coords: [0.0, 0.0]}
          - {id: 2, coords: [1.0, 0.0]}
        materials:
          - {id: 1, type: Elastic1D, E: 1000.0, density: 1.0}
        elements:
          - {id: 1, type: Truss2D, material: 1, nodes: [1, 2], A: 1.0}
        boundary_conditions:
          - {node_id: 1, ux: 0.0, uy: 0.0}
          - {node_id: 2, uy: 0.0}
        point_loads:
          - {node_id: 2, ux: 1.0}
        solver: {type: NewmarkSolver, t_end: 0.01, dt: 0.001}
        output:
          history:
            columns:
              - {name: W, energy: external_work}
        """), encoding="utf-8")
    with pytest.raises(ValueError, match="estáticos"):
        solidum.run_yaml(ruta, write_output=True)


def test_la_zona_usa_la_regla_de_integracion_de_cada_elemento():
    """Dos Quad4 con reglas distintas del mismo número de puntos (Gauss 2×2 y
    Newton–Cotes en las esquinas): ``G_e·σ_e`` de la zona es la fuerza
    interna de cada elemento con su propia regla."""
    from solidum.core.domain import Domain
    from solidum.registry import ElementRegistry, MaterialRegistry
    from solidum.utils.output_writer import _ZoneGroup
    mat = MaterialRegistry.create("Elastic2D", E=1000.0, nu=0.25)
    dom = Domain()
    coords = [(0, 0), (1, 0), (2, 0), (0, 1), (1.1, 1.2), (2, 1)]
    for k, c in enumerate(coords, start=1):
        dom.add_node(k, list(map(float, c)))
    Q4 = ElementRegistry.get("Quad4")
    esquinas = (np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]]), np.ones(4))
    n = dom.nodes
    e1 = Q4(1, [n[1], n[2], n[5], n[4]], mat, thickness=3.0)
    e2 = Q4(2, [n[2], n[3], n[6], n[5]], mat, thickness=3.0, quadrature=esquinas)
    dom.add_element(e1)
    dom.add_element(e2)
    dom.generate_equation_numbers()
    U = np.random.default_rng(0).normal(size=dom.total_dofs) * 1e-3
    f_ref = [el.compute_element_state(U[el.get_global_dof_indices()])[1] for el in (e1, e2)]
    grupo = _ZoneGroup([e1, e2])
    f = np.einsum("mij,mj->mi", grupo.G, grupo.stresses(trial=True))
    for m in range(2):
        np.testing.assert_allclose(f[m], f_ref[m], rtol=1e-12, atol=1e-12 * np.abs(f_ref[m]).max())
