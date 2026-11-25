"""Figuras de los resultados desde el YAML (ADR 0039).

Las figuras se dibujan con lo que Solidum escribe (VTK por instante con su
``.pvd``, ``_final.vtu`` e historia CSV): al terminar ``run_yaml(...,
write_output=True)`` y, sin volver a calcular, con ``python -m solidum.post``.
Se comprueba que se escriben imágenes válidas del tamaño esperado, que la
regeneración no resuelve de nuevo, la comparación de curvas de dos modelos,
la validación y los mensajes de error.
"""
import textwrap
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("matplotlib")

import solidum  # noqa: E402
from solidum.post import ResultFiles, make_plots, validate_plots_block  # noqa: E402
from solidum.post.__main__ import main as post_main  # noqa: E402

MODELO = """\
materials:
  - {{id: 1, type: IsotropicDamageCrackBand2D, E: 30000.0, nu: 0.0, f_t: 3.0, G_f: 0.1}}
  - {{id: 2, type: IsotropicDamageCrackBand2D, E: 30000.0, nu: 0.0, f_t: 2.85, G_f: 0.1}}
mesh_element: {element}
mesh_thickness: 10.0
mesh:
  structured:
    blocks:
      - {{name: a, material: 1, x: {{length: 50, n: 5}}, y: {{length: 20, n: 2}}}}
      - {{name: w, material: 2, origin: [50, 0], x: {{length: 10, n: 1}}, y: {{length: 20, n: 2}}}}
      - {{name: b, material: 1, origin: [60, 0], x: {{length: 40, n: 4}}, y: {{length: 20, n: 2}}}}
    groups:
      Origen: "near(x, 0) and near(y, 0)"
boundary_conditions_by_group:
  a.xmin: {{ux: 0.0}}
  Origen: {{uy: 0.0}}
  b.xmax: {{ux: 0.4}}
solver: {{type: NonlinearSolver, num_steps: 100}}
output:
  file_name: {name}
  log_level: WARNING
  results:
    frequency: {{at: [0.05, 1.0]}}
  history:
    columns:
      - {{name: u, displacement: {{group: b.xmax, dof: ux}}}}
      - {{name: P, reaction: {{group: b.xmax, dof: ux}}}}
      - {{name: D, energy: dissipated}}
  plots:
{plots}
"""

PLOTS = """\
    - {type: mesh, file: malla}
    - {type: deformed, file: deformada.png, at: [0.05, 1.0], field: Internal_State,
       field_label: "daño", label_column: u}
    - {type: deformed, file: zoom.pdf, at: [final], field: displacement, zoom: [40, 70, 0, 20],
       scale: 20}
    - {type: curve, file: curva, x: u, y: [P, D]}
"""


def _modelo(d: Path, name="barra", element="Quad4", plots=PLOTS) -> Path:
    ruta = d / f"{name}.yaml"
    ruta.write_text(MODELO.format(name=name, element=element, plots=textwrap.indent(
        textwrap.dedent(plots), "    ").rstrip()), encoding="utf-8")
    return ruta


def _png_size(path: Path):
    from matplotlib.image import imread
    return imread(path).shape


@pytest.fixture(scope="module")
def resuelto(tmp_path_factory):
    d = tmp_path_factory.mktemp("figuras")
    ruta = _modelo(d)
    solidum.run_yaml(ruta, write_output=True)
    return d, ruta


def test_se_escriben_las_figuras_al_resolver(resuelto):
    d, _ = resuelto
    for name in ("malla.png", "deformada.png", "zoom.pdf", "curva.png"):
        assert (d / name).is_file() and (d / name).stat().st_size > 2000, name
    h, w, _ = _png_size(d / "malla.png")
    assert w > 2 * h                       # barra de 100 × 20: figura apaisada
    h2, w2, _ = _png_size(d / "deformada.png")
    assert h2 > h                          # dos paneles
    assert (d / "barra.pvd").is_file()


def test_post_regenera_sin_resolver(resuelto, monkeypatch):
    d, ruta = resuelto
    for p in d.glob("*.png"):
        p.unlink()

    def no_resolver(*a, **k):
        raise AssertionError("python -m solidum.post no debe resolver el modelo")
    monkeypatch.setattr(solidum, "run_yaml", no_resolver)
    assert post_main([str(ruta)]) == 0
    assert (d / "malla.png").is_file() and (d / "curva.png").is_file()


def test_comparar_curvas_de_dos_modelos(tmp_path):
    a = _modelo(tmp_path, "q4", "Quad4", "    - {type: mesh}\n")
    b = _modelo(tmp_path, "q8", "Quad8", "    - {type: mesh}\n")
    solidum.run_yaml(a, write_output=True)
    solidum.run_yaml(b, write_output=True)
    out = make_plots(ResultFiles.of_yaml(a), [
        {"type": "curve", "file": "comparacion", "x": "u", "y": "P", "compare": ["q8.yaml"],
         "labels": ["Quad4", "Quad8"]}])
    assert out[0].name == "comparacion.png" and out[0].stat().st_size > 2000
    assert (tmp_path / "q4_mesh_1.png").is_file()      # nombre por omisión


def test_errores_de_las_figuras(resuelto):
    d, ruta = resuelto
    files = ResultFiles.of_yaml(ruta)
    with pytest.raises(ValueError, match="inexistente en el VTK"):
        make_plots(files, [{"type": "deformed", "field": "Temperatura"}])
    with pytest.raises(ValueError, match="no tiene la columna"):
        make_plots(files, [{"type": "curve", "x": "u", "y": "Q"}])
    with pytest.raises(ValueError, match="label_column"):
        make_plots(files, [{"type": "deformed", "at": [1.0], "label_column": "v"}])


def test_un_campo_que_no_existe_no_tira_el_analisis(tmp_path):
    ruta = _modelo(tmp_path, plots="    - {type: deformed, field: Nada}\n")
    solidum.run_yaml(ruta, write_output=True)          # sin excepción
    assert (tmp_path / "barra_historia.csv").is_file()
    assert post_main([str(ruta)]) == 3


def test_sin_instantes_es_un_error_claro(tmp_path):
    ruta = _modelo(tmp_path, plots="    - {type: mesh}\n")
    texto = ruta.read_text(encoding="utf-8").replace("    frequency: {at: [0.05, 1.0]}\n", "")
    ruta.write_text(texto.replace("  results:\n", ""), encoding="utf-8")
    solidum.run_yaml(ruta, write_output=True)
    with pytest.raises(ValueError, match="frequency"):
        make_plots(ResultFiles.of_yaml(ruta), [{"type": "deformed", "at": [0.5]}])
    # 'final' siempre existe.
    out = make_plots(ResultFiles.of_yaml(ruta), [{"type": "deformed", "at": ["final"]}])
    assert out[0].is_file()


def test_malla_3d_no_disponible(tmp_path):
    ruta = tmp_path / "cubo.yaml"
    ruta.write_text(textwrap.dedent("""\
        materials:
          - {id: 1, type: Elastic3D, E: 1.0, nu: 0.2}
        mesh:
          structured: {element: Hex8, x: {length: 1, n: 1}, y: {length: 1, n: 1}, z: {length: 1, n: 1}}
        boundary_conditions_by_group:
          xmin: {ux: 0.0, uy: 0.0, uz: 0.0}
        point_loads_by_group:
          xmax: {ux: 1.0}
        solver: {type: LinearSolver}
        output:
          file_name: cubo
          plots:
            - {type: mesh}
        """), encoding="utf-8")
    solidum.run_yaml(ruta, write_output=True)
    with pytest.raises(ValueError, match="3D"):
        make_plots(ResultFiles.of_yaml(ruta), [{"type": "mesh"}])


@pytest.mark.parametrize("plots,mensaje", [
    ([], "lista de figuras"),
    ([{"type": "contorno"}], "type: mesh"),
    ([{"type": "mesh", "field": "x"}], "no admitidas"),
    ([{"type": "mesh", "file": "a.jpg"}], "no admitido"),
    ([{"type": "mesh", "file": "a"}, {"type": "curve", "file": "a", "x": "u", "y": "P"}], "repetido"),
    ([{"type": "deformed", "at": "final"}], "lista de factores"),
    ([{"type": "deformed", "scale": -1}], "'auto' o un número"),
    ([{"type": "deformed", "edges": "si"}], "'auto', true o false"),
    ([{"type": "deformed", "zoom": [0, 1]}], "xmin, xmax, ymin, ymax"),
    ([{"type": "curve", "y": "P"}], "eje x"),
    ([{"type": "curve", "x": "u", "y": 3}], "lista de columnas"),
    ([{"type": "curve", "x": "u", "y": "P", "compare": "b.yaml"}], "lista de archivos"),
    ([{"type": "curve", "x": "u", "y": "P", "xlim": [1, 0]}], "mínimo < máximo"),
    ([{"type": "curve", "x": "u", "y": "P", "yscale": "logaritmica"}], "'linear' o 'log'"),
    ([{"type": "curve", "x": "u", "y": "P", "reference": {"file": "a.csv"}}],
     "lista de curvas o bandas"),
    ([{"type": "curve", "x": "u", "y": "P", "reference": [{"label": "ensayo"}]}],
     "'file' (un CSV con encabezado) o 'points'"),
    ([{"type": "curve", "x": "u", "y": "P", "reference": [{"file": "a.csv", "y": "P"}]}],
     "columna del eje x"),
    ([{"type": "curve", "x": "u", "y": "P",
       "reference": [{"file": "a.csv", "x": "d", "y": ["a", "b", "c"]}]}], "[mínimo, máximo]"),
    ([{"type": "curve", "x": "u", "y": "P", "reference": [{"points": [[0, 0]]}]}],
     "al menos dos puntos"),
    ([{"type": "curve", "x": "u", "y": "P", "reference": [{"points": [[0, 0], [1, 2, 3]]}]}],
     "[x, y_min, y_max]"),
    ([{"type": "curve", "x": "u", "y": "P", "reference": [{"points": [[0, 0], [1, 1], 5]}]}],
     "[x, y_min, y_max]"),
    ([{"type": "curve", "x": "u", "y": "P",
       "reference": [{"points": [[0, 0], [1, 1]], "style": "puntos"}]}], "'line' o 'markers'"),
    ([{"type": "curve", "x": "u", "y": "P",
       "reference": [{"points": [[0, 0], [1, 1]], "scale": 2}]}], "[factor x, factor y]"),
])
def test_validacion(plots, mensaje):
    errores = validate_plots_block(plots)
    assert any(mensaje in e for e in errores), errores


def test_el_yaml_valida_las_figuras(tmp_path):
    from solidum.utils.yaml_parser import YamlValidationError
    ruta = _modelo(tmp_path, plots="    - {type: contorno}\n")
    with pytest.raises(YamlValidationError, match="output.plots"):
        solidum.run_yaml(ruta)


def test_paleta_de_curvas_en_orden_fijo(resuelto):
    """Las dos curvas de la figura llevan los dos primeros colores de la paleta."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from solidum.post import plots as P
    d, ruta = resuelto
    created = []
    real = plt.subplots

    def spy(*a, **k):
        fig, ax = real(*a, **k)
        created.append(ax)
        return fig, ax
    try:
        plt.subplots = spy
        P.plot_curve(ResultFiles.of_yaml(ruta), {"type": "curve", "x": "u", "y": ["P", "D"]},
                     d / "colores.png")
    finally:
        plt.subplots = real
    colors = [matplotlib.colors.to_hex(line.get_color()) for line in created[0].get_lines()]
    assert colors == [P._SERIES[0], P._SERIES[1]]
    np.testing.assert_equal(len(created[0].get_legend().get_texts()), 2)


# ----------------------------------------------------------------------------
# Referencias externas de las curvas (ensayo, otra publicación)
# ----------------------------------------------------------------------------

def _ejes_de(spec, files, out):
    """Dibuja ``spec`` y devuelve sus ejes (espía de ``plt.subplots``)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from solidum.post import plots as P
    created = []
    real = plt.subplots

    def spy(*a, **k):
        fig, ax = real(*a, **k)
        created.append(ax)
        return fig, ax
    try:
        plt.subplots = spy
        P.plot_curve(files, spec, out)
    finally:
        plt.subplots = real
    return created[0]


def test_curva_con_referencias_de_un_csv_y_en_linea(resuelto):
    """Curva y banda de un CSV (separador ';', unidades convertidas con
    ``scale``) y una curva en línea: en tinta neutra, debajo de las de
    Solidum, que conservan la paleta; límites y escala de los ejes."""
    import matplotlib
    from solidum.post import plots as P
    d, ruta = resuelto
    (d / "ensayo.csv").write_text("d_mm;P_kN;Pmin_kN;Pmax_kN\n0;0;0;0\n0.01;0.5;0.45;0.55\n"
                                  "0.2;0.1;;0.2\n", encoding="utf-8")
    spec = {"type": "curve", "x": "u", "y": "P", "xlim": [0.0, 0.3], "yscale": "log",
            "reference": [
                {"file": "ensayo.csv", "x": "d_mm", "y": "P_kN", "scale": [1, 1000],
                 "label": "ensayo", "style": "markers"},
                {"file": "ensayo.csv", "x": "d_mm", "y": ["Pmin_kN", "Pmax_kN"],
                 "scale": [1, 1000], "label": "envolvente"},
                {"points": [[0.0, 1.0], [0.3, 500.0]], "label": "otra"}]}
    ax = _ejes_de(spec, ResultFiles.of_yaml(ruta), d / "referencias.png")
    lines = ax.get_lines()
    assert [ln.get_label() for ln in lines] == ["ensayo", "otra", "barra"]
    np.testing.assert_allclose(lines[0].get_ydata(), [0.0, 500.0, 100.0])   # kN → N
    assert matplotlib.colors.to_hex(lines[-1].get_color()) == P._SERIES[0]
    assert all(matplotlib.colors.to_hex(ln.get_color()) == P._INK_2 for ln in lines[:-1])
    assert all(ln.get_zorder() < lines[-1].get_zorder() for ln in lines[:-1])
    assert len(ax.collections) == 1                      # la banda
    assert ax.get_xlim() == (0.0, 0.3) and ax.get_yscale() == "log"
    assert len(ax.get_legend().get_texts()) == 4


def test_errores_de_las_referencias(resuelto):
    d, ruta = resuelto
    files = ResultFiles.of_yaml(ruta)
    (d / "malo.csv").write_text("d,P\n0,0\n0.1,mucho\n", encoding="utf-8")
    (d / "bueno.csv").write_text("d,P\n0,0\n0.1,5\n", encoding="utf-8")
    base = {"type": "curve", "x": "u", "y": "P"}
    with pytest.raises(FileNotFoundError, match="no_existe.csv"):
        make_plots(files, [{**base, "reference": [{"file": "no_existe.csv", "x": "d", "y": "P"}]}])
    with pytest.raises(ValueError, match="no tiene las columnas"):
        make_plots(files, [{**base, "reference": [{"file": "bueno.csv", "x": "d", "y": "Q"}]}])
    with pytest.raises(ValueError, match="línea 3"):
        make_plots(files, [{**base, "reference": [{"file": "malo.csv", "x": "d", "y": "P"}]}])


def test_celdas_sin_el_campo_en_gris(tmp_path):
    """Un campo que sólo tiene un material (NaN en las demás celdas) se dibuja
    con esas celdas en gris, no transparentes."""
    texto = MODELO.format(name="mixto", element="Quad4", plots="    - {type: mesh}")
    texto = texto.replace("{id: 1, type: IsotropicDamageCrackBand2D, E: 30000.0, nu: 0.0, "
                          "f_t: 3.0, G_f: 0.1}", "{id: 1, type: Elastic2D, E: 30000.0, nu: 0.0}")
    ruta = tmp_path / "mixto.yaml"
    ruta.write_text(texto, encoding="utf-8")
    solidum.run_yaml(ruta, write_output=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from solidum.post import plots as P
    created = []
    real = plt.subplots

    def spy(*a, **k):
        fig, axes = real(*a, **k)
        created.append(axes)
        return fig, axes
    try:
        plt.subplots = spy
        P.plot_deformed(ResultFiles.of_yaml(ruta), {"type": "deformed", "field": "damage"},
                        tmp_path / "d.png")
    finally:
        plt.subplots = real
    pc = created[0][0, 0].collections[0]
    assert np.isnan(np.asarray(pc.get_array(), dtype=float)).sum() == 18   # 18 elásticas
    assert matplotlib.colors.to_hex(pc.get_cmap().get_bad()) == P._NO_FIELD
