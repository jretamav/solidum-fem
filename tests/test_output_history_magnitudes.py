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
     "sólo con var: stress o strain"),
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
