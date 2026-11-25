"""Un análisis que se detiene conserva lo calculado (ADR 0036, ampliación del
2026-11-25).

Barra de la prueba de objetividad de la banda de grieta (tres bloques, la
columna central debilitada) con una **carga** en el extremo mayor que su
resistencia: el Newton no puede pasar del pico y el análisis se detiene. Antes
sólo quedaba la historia (escrita fila a fila); ahora se escriben además el
VTK final, el ``.pvd`` y las figuras con el último paso convergido, y el
error original sigue su curso. ``python -m solidum`` lo explica y termina con
el código 3.
"""
import csv
import textwrap

import numpy as np
import pytest

import solidum
from solidum.math.solvers.diagnostics import SolverDivergedError

pytest.importorskip("matplotlib")

MODELO = """\
materials:
  - {id: 1, type: IsotropicDamageCrackBand2D, E: 30000.0, nu: 0.0, f_t: 3.0, G_f: 0.1}
  - {id: 2, type: IsotropicDamageCrackBand2D, E: 30000.0, nu: 0.0, f_t: 2.85, G_f: 0.1}
mesh_element: Quad4
mesh_thickness: 10.0
mesh:
  structured:
    blocks:
      - {name: a, material: 1, x: {length: 50, n: 5}, y: {length: 20, n: 2}}
      - {name: w, material: 2, origin: [50, 0], x: {length: 10, n: 1}, y: {length: 20, n: 2}}
      - {name: b, material: 1, origin: [60, 0], x: {length: 40, n: 4}, y: {length: 20, n: 2}}
    groups:
      Origen: "near(x, 0) and near(y, 0)"
boundary_conditions_by_group:
  a.xmin: {ux: 0.0}
  Origen: {uy: 0.0}
point_loads_by_group:
  b.xmax: {ux: 400.0}
solver: {type: NonlinearSolver, num_steps: 20}
output:
  file_name: det
  log_level: ERROR
  results:
    frequency: {at: [0.2, 1.0]}
    integration_points: true
  history:
    columns:
      - {name: u, displacement: {group: b.xmax, dof: ux}}
  plots:
    - {type: deformed, file: deformada, at: [final], field: Internal_State}
    - {type: curve, file: curva, x: u, y: factor}
"""


def _modelo(tmp_path):
    ruta = tmp_path / "det.yaml"
    ruta.write_text(MODELO, encoding="utf-8")
    return ruta


def _historia(path):
    with open(path, newline="", encoding="utf-8") as f:
        filas = list(csv.reader(f))
    return {c: np.array([float(r[k]) for r in filas[1:]]) for k, c in enumerate(filas[0])}


def test_el_analisis_detenido_escribe_el_ultimo_paso_convergido(tmp_path):
    import meshio
    ruta = _modelo(tmp_path)
    with pytest.raises(SolverDivergedError):
        solidum.run_yaml(ruta, write_output=True)
    h = _historia(tmp_path / "det_historia.csv")
    assert 0.0 < h["factor"][-1] < 1.0                   # se detuvo antes del final
    for nombre in ("det_final.vtu", "det.pvd", "det_puntos_final.vtu", "det_puntos.pvd",
                   "deformada.png", "curva.png"):
        assert (tmp_path / nombre).is_file(), nombre
    # El VTK final es el del último paso de la historia.
    m = meshio.read(tmp_path / "det_final.vtu")
    ux = np.asarray(m.point_data["Displacements"])[:, 0]
    extremo = np.isclose(m.points[:, 0], 100.0)
    assert np.mean(ux[extremo]) == pytest.approx(h["u"][-1], rel=1e-9)


def test_python_m_solidum_explica_la_detencion(tmp_path, capsys):
    from solidum.__main__ import main
    ruta = _modelo(tmp_path)
    assert main([str(ruta)]) == 3
    err = capsys.readouterr().err
    assert "El análisis se detuvo" in err and "Resultados escritos hasta el último paso" in err
    assert (tmp_path / "det_final.vtu").is_file()


def test_sin_pasos_convergidos_la_terminal_no_anuncia_resultados(tmp_path, capsys):
    from solidum.__main__ import main
    ruta = tmp_path / "det.yaml"
    ruta.write_text(MODELO.replace("  a.xmin: {ux: 0.0}\n", ""), encoding="utf-8")
    assert main([str(ruta)]) == 3
    err = capsys.readouterr().err
    assert "El análisis se detuvo" in err and "Resultados escritos" not in err


def test_un_error_antes_del_primer_paso_no_escribe_resultados(tmp_path):
    """Sin ningún paso convergido (aquí, un mecanismo: la barra sin apoyo en x)
    no hay estado que escribir: sólo el error."""
    from solidum.math.solvers.diagnostics import MechanismError
    ruta = tmp_path / "det.yaml"
    ruta.write_text(MODELO.replace("  a.xmin: {ux: 0.0}\n", ""), encoding="utf-8")
    with pytest.raises(MechanismError):
        solidum.run_yaml(ruta, write_output=True)
    assert not (tmp_path / "det_final.vtu").exists()


def test_el_fallo_del_cierre_no_tapa_el_error_del_analisis(tmp_path, monkeypatch):
    from solidum.utils import output_writer as ow

    def roto(self, U_final=None):
        raise OSError("disco lleno")
    monkeypatch.setattr(ow.OutputWriter, "finish", roto)
    with pytest.raises(SolverDivergedError):
        solidum.run_yaml(_modelo(tmp_path), write_output=True)


def test_el_resultado_normal_no_cambia(tmp_path):
    """Con desplazamiento impuesto el mismo modelo termina y escribe lo de siempre."""
    ruta = tmp_path / "det.yaml"
    ruta.write_text(textwrap.dedent(MODELO).replace(
        "point_loads_by_group:\n  b.xmax: {ux: 400.0}\n", "").replace(
        "  Origen: {uy: 0.0}\n", "  Origen: {uy: 0.0}\n  b.xmax: {ux: 0.4}\n"), encoding="utf-8")
    solidum.run_yaml(ruta, write_output=True)
    h = _historia(tmp_path / "det_historia.csv")
    assert h["factor"][-1] == pytest.approx(1.0)
    assert (tmp_path / "det_final.vtu").is_file()
