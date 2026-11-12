"""Salida de resultados del YAML (ADR 0036): historia y frecuencia de VTK.

Barra de dos tramos (Truss2D, EA = 2000, L = 2) con el extremo derecho
desplazado u = 0.01 en 10 pasos: la reacción del extremo es P = EA/L·u en
cada paso y la del apoyo, −P.
"""
import csv
import textwrap
from pathlib import Path

import numpy as np
import pytest

import solidum
from solidum.utils.output_writer import validate_output_block
from solidum.utils.yaml_parser import YamlValidationError

EA, L, U_END, STEPS = 2000.0, 2.0, 0.01, 10


def _modelo(tmp_path, output: str, solver: str = None) -> Path:
    ruta = tmp_path / "barra.yaml"
    ruta.write_text(textwrap.dedent("""\
        nodes:
          - {id: 1, coords: [0.0, 0.0]}
          - {id: 2, coords: [1.0, 0.0]}
          - {id: 3, coords: [2.0, 0.0]}
        materials:
          - {id: 1, type: Elastic1D, E: 2000.0}
        elements:
          - {id: 1, type: Truss2D, material: 1, nodes: [1, 2], A: 1.0}
          - {id: 2, type: Truss2D, material: 1, nodes: [2, 3], A: 1.0}
        boundary_conditions:
          - {node_id: 1, ux: 0.0, uy: 0.0}
          - {node_id: 2, uy: 0.0}
          - {node_id: 3, ux: 0.01, uy: 0.0}
        """) + textwrap.dedent(solver or """\
        solver:
          type: NonlinearSolver
          num_steps: 10
        """) + textwrap.dedent(output), encoding="utf-8")
    return ruta


_HIST = """\
    output:
      file_name: barra
      history:
        columns:
          - {name: u, displacement: {node: 3, dof: ux}}
          - {name: P, reaction: {node: 3, dof: ux}}
          - {name: R1, reaction: {at: [0.0, 0.0], dof: ux}}
          - {name: u_media, displacement: {nodes: [1, 2, 3], dof: ux, reduce: mean}}
    """


def _leer(path):
    with open(path, newline="", encoding="utf-8") as f:
        filas = list(csv.reader(f))
    return filas[0], np.array(filas[1:], dtype=float)


def test_historia_reaccion_y_desplazamiento(tmp_path):
    ruta = _modelo(tmp_path, _HIST)
    solidum.run_yaml(ruta, write_output=True)
    cab, datos = _leer(tmp_path / "barra_historia.csv")
    assert cab == ["step", "factor", "u", "P", "R1", "u_media"]
    assert len(datos) == STEPS
    np.testing.assert_allclose(datos[:, 1], np.arange(1, STEPS + 1) / STEPS)
    np.testing.assert_allclose(datos[:, 2], datos[:, 1] * U_END, rtol=1e-9)
    np.testing.assert_allclose(datos[:, 3], EA / L * datos[:, 2], rtol=1e-8)
    np.testing.assert_allclose(datos[:, 4], -datos[:, 3], rtol=1e-8)
    np.testing.assert_allclose(datos[:, 5], datos[:, 2] / 2.0, rtol=1e-9)
    assert (tmp_path / "barra_final.vtu").is_file()


def test_sin_write_output_no_se_escribe_nada(tmp_path):
    ruta = _modelo(tmp_path, _HIST)
    solidum.run_yaml(ruta)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["barra.yaml"]


@pytest.mark.parametrize("freq,pasos", [
    ('"all_steps"', list(range(1, 11))),
    ("{every: 4}", [4, 8]),
    ("{at: [0.25, 0.5, 1.0]}", [3, 5, 10]),
    ('"last_step"', []),
])
def test_frecuencia_de_vtk(tmp_path, freq, pasos):
    ruta = _modelo(tmp_path, f"""\
        output:
          file_name: barra
          results:
            frequency: {freq}
        """)
    solidum.run_yaml(ruta, write_output=True)
    escritos = sorted(int(p.stem.split("_")[-1]) for p in tmp_path.glob("barra_[0-9]*.vtu"))
    assert escritos == pasos
    assert (tmp_path / "barra_final.vtu").is_file()
    assert (tmp_path / "barra.pvd").is_file() == (len(pasos) >= 1)


def test_solver_lineal_escribe_una_fila(tmp_path):
    ruta = _modelo(tmp_path, _HIST, solver="solver: {type: LinearSolver}\n")
    solidum.run_yaml(ruta, write_output=True)
    _, datos = _leer(tmp_path / "barra_historia.csv")
    assert datos.shape[0] == 1
    assert datos[0, 3] == pytest.approx(EA / L * U_END, rel=1e-9)


def test_python_m_solidum(tmp_path):
    from solidum.__main__ import main
    ruta = _modelo(tmp_path, _HIST)
    assert main([str(ruta)]) == 0
    assert (tmp_path / "barra_historia.csv").is_file()


@pytest.mark.parametrize("bloque,mensaje", [
    ({"results": {"frequency": "a_veces"}}, "no admitido"),
    ({"results": {"frequency": {"every": 0}}}, "entero ≥ 1"),
    ({"results": {"frequency": {"at": "fin"}}}, "lista de números"),
    ({"history": {"columns": []}}, "se esperaba {columns"),
    ({"history": {"columns": [{"name": "P", "reaction": {"node": 3}}]}}, "falta el grado de libertad"),
    ({"history": {"columns": [{"name": "P", "fuerza": {"node": 3, "dof": "ux"}}]}}, "exactamente una magnitud"),
    ({"history": {"columns": [{"name": "u", "displacement": {"node": 3, "group": "A", "dof": "ux"}}]}},
     "exactamente uno de"),
    ({"history": {"columns": [{"name": "step", "displacement": {"node": 3, "dof": "ux"}}]}}, "nombre único"),
    ({"history": {"columns": [{"name": "u", "displacement": {"node": 3, "dof": "ux", "reduce": "sum"}}]}},
     "'mean', 'min' o 'max'"),
])
def test_validacion(bloque, mensaje):
    errores = validate_output_block(bloque)
    assert any(mensaje in e for e in errores), errores


def test_el_yaml_rechaza_una_frecuencia_invalida(tmp_path):
    ruta = _modelo(tmp_path, """\
        output:
          results: {frequency: siempre}
        """)
    with pytest.raises(YamlValidationError, match="frequency"):
        solidum.run_yaml(ruta)


def test_reaccion_en_un_nodo_libre_es_un_error(tmp_path):
    ruta = _modelo(tmp_path, """\
        output:
          history:
            columns:
              - {name: P, reaction: {node: 2, dof: ux}}
        """)
    with pytest.raises(ValueError, match="necesita 'ux' prescrito"):
        solidum.run_yaml(ruta, write_output=True)
