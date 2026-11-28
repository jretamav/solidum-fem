"""Cargas por unidad de longitud en barras y marcos y fuerzas internas a lo
largo del elemento (ADR 0044; specs Frame2DEuler §14–§15, Frame3D §13–§14,
Frame2DTimoshenko §14, Truss2D §12).

Parámetros ≠ 1 y elementos oblicuos (Reglas §6): con valores unitarios un
factor olvidado es la identidad. Los valores esperados son las formas
cerradas de la spec (verificadas por el teórico el 2026-11-27).
"""
from __future__ import annotations

import copy
import pickle
import textwrap

import numpy as np
import pytest

from solidum import Assembler, Domain
from solidum.core.element import Element
from solidum.core.load_vector import LoadVector
from solidum.elements._line import LineElementMixin
from solidum.elements.frame import Frame2DEuler, Frame2DTimoshenko
from solidum.elements.frame3d import Frame3D
from solidum.elements.truss import Truss2D
from solidum.entry import load_yaml
from solidum.materials.elastic import Elastic1D
from solidum.registry import ElementRegistry
from solidum.utils.yaml_parser import YamlValidationError

E, A, I = 2.1e8, 0.013, 2.1e-4
L, ANG = 3.7, np.deg2rad(37.0)
X0 = np.array([0.4, -0.2])
T = np.array([np.cos(ANG), np.sin(ANG)])
W, Q0 = 5.3e3, 7.1e3


def frame(cls=Frame2DEuler, **kw):
    dom = Domain()
    n1 = dom.add_node(1, list(X0))
    n2 = dom.add_node(2, list(X0 + L * T))
    mat = Elastic1D(E=E, density=7850.0)
    if cls is Frame2DTimoshenko:
        kw.setdefault("nu", 0.3)
    el = cls(1, [n1, n2], mat, A=A, I=I, **kw)
    dom.add_element(el)
    dom.generate_equation_numbers()
    return el, np.zeros(dom.total_dofs)


def triangular_local(x):
    """De 0 en i a −Q0 en j, en y' (hacia −y')."""
    s = (np.asarray(x, dtype=float) - X0) @ T
    return [0.0, -Q0 * s / L]


def clamped(el, U, samples):
    """Fuerzas de una viga biempotrada de un elemento (U = 0): −f_eq."""
    f = el.line_load_vector(samples)
    ef = el.internal_forces_from_nodal(np.zeros_like(f), U, equivalent_load=f)
    return f, ef


# ---------------------------------------------------------------------------
# Vector nodal consistente (spec Frame2DEuler §14)
# ---------------------------------------------------------------------------

class TestFormasCerradas:

    def test_uniforme_y_triangular_en_locales(self):
        el, U = frame()
        f = el.compute_line_load([0.0, W], local=True)
        np.testing.assert_allclose(el.T @ f, [0, W * L / 2, W * L**2 / 12,
                                              0, W * L / 2, -W * L**2 / 12], rtol=1e-12)
        s = el.line_load_samples(lambda x: [0.0, Q0 * ((np.asarray(x) - X0) @ T) / L], local=True)
        np.testing.assert_allclose(el.T @ el.line_load_vector(s),
                                   [0, 3 * Q0 * L / 20, Q0 * L**2 / 30,
                                    0, 7 * Q0 * L / 20, -Q0 * L**2 / 20], rtol=1e-12, atol=1e-9)

    def test_global_igual_a_local_y_a_fuerza_de_cuerpo(self):
        el, _ = frame()
        q_local = np.array([1.3e3, -W])
        q_global = el._load_axes().T @ q_local
        f_loc = el.compute_line_load(q_local, local=True)
        np.testing.assert_allclose(el.compute_line_load(q_global), f_loc, rtol=1e-12, atol=1e-9)
        np.testing.assert_allclose(el.compute_body_load(q_global / A), f_loc, rtol=1e-12, atol=1e-9)

    def test_proyectada_por_longitud_horizontal(self):
        el, _ = frame()
        qy = -W
        f = el.compute_line_load([0.0, qy], projected=True)
        assert f[1] + f[4] == pytest.approx(qy * L * abs(T[0]), rel=1e-12)
        qx = 2.7e3
        f = el.compute_line_load([qx, 0.0], projected=True)
        assert f[0] + f[3] == pytest.approx(qx * L * abs(T[1]), rel=1e-12)

    def test_timoshenko_triangular_depende_de_phi(self):
        el, U = frame(Frame2DTimoshenko, As=0.011)
        phi = el._phi()
        assert phi > 1e-3
        s = el.line_load_samples(lambda x: [0.0, Q0 * ((np.asarray(x) - X0) @ T) / L], local=True)
        f = el.T @ el.line_load_vector(s)
        k = 1.0 + phi
        np.testing.assert_allclose(
            f, [0, Q0 * L * (9 + 10 * phi) / (60 * k), Q0 * L**2 * (4 + 5 * phi) / (120 * k),
                0, Q0 * L * (21 + 20 * phi) / (60 * k), -Q0 * L**2 * (6 + 5 * phi) / (120 * k)],
            rtol=1e-12, atol=1e-9)

    def test_frame3d_uniforme_en_y_y_z_locales(self):
        dom = Domain()
        a = np.array([0.3, -1.1, 0.7])
        b = a + np.array([1.1, 2.0, 0.7])
        n1, n2 = dom.add_node(1, list(a)), dom.add_node(2, list(b))
        el = Frame3D(1, [n1, n2], Elastic1D(E=E), A=A, Iy=I, Iz=1.7 * I, J=0.9 * I, nu=0.3,
                     ref_vector=[0.2, 0.1, 1.0])
        qy, qz = 1.7e3, -2.9e3
        Le = el.L0
        f = el.T @ el.compute_line_load([0.0, qy, qz], local=True)
        np.testing.assert_allclose(f[[1, 5, 7, 11]], [qy * Le / 2, qy * Le**2 / 12,
                                                     qy * Le / 2, -qy * Le**2 / 12], rtol=1e-12)
        np.testing.assert_allclose(f[[2, 4, 8, 10]], [qz * Le / 2, -qz * Le**2 / 12,
                                                     qz * Le / 2, qz * Le**2 / 12], rtol=1e-12)
        g = el.compute_line_load([0.0, 0.0, qz], projected=True)
        assert g[2] + g[8] == pytest.approx(qz * np.hypot(1.1, 2.0), rel=1e-12)


# ---------------------------------------------------------------------------
# Fuerzas internas a lo largo del elemento (spec Frame2DEuler §15)
# ---------------------------------------------------------------------------

class TestDiagramas2D:

    def test_biempotrada_uniforme(self):
        el, U = frame()
        s = el.line_load_samples([0.0, -W], local=True)
        _, ef = clamped(el, U, s)
        al = el.internal_forces_along(ef, U, [0.0, 0.5, 1.0], s)
        np.testing.assert_allclose(al.components["M"], np.array([-1, 0.5, -1]) * W * L**2 / 12,
                                   rtol=1e-12)
        np.testing.assert_allclose(al.components["V"], [W * L / 2, 0.0, -W * L / 2],
                                   rtol=1e-12, atol=1e-9)

    def test_biempotrada_triangular_y_cierre_en_j(self):
        el, U = frame()
        s = el.line_load_samples(triangular_local, local=True)
        _, ef = clamped(el, U, s)
        al = el.internal_forces_along(ef, U, [0.0, 1.0], s)
        np.testing.assert_allclose(al.components["M"], [-Q0 * L**2 / 30, -Q0 * L**2 / 20],
                                   rtol=1e-11)
        np.testing.assert_allclose(al.components["V"], [3 * Q0 * L / 20, -7 * Q0 * L / 20],
                                   rtol=1e-11)
        for c in ("N", "V", "M"):
            assert al.components[c][1] == pytest.approx(ef.components[c][1], abs=1e-9 * Q0 * L**2)

    def test_carga_no_polinomica_cierra_con_la_representada(self):
        el, U = frame()

        def q(x):
            xi = ((np.asarray(x) - X0) @ T) / L
            return [0.0, -(2.3 * np.exp(1.3 * xi) + 5 * np.sin(5 * xi) - 1.1) * 1e3]
        s = el.line_load_samples(q, local=True)
        _, ef = clamped(el, U, s)
        al = el.internal_forces_along(ef, U, [1.0], s)
        assert al.components["M"][0] == pytest.approx(ef.components["M"][1], abs=1e-9 * 1e4)

    def test_derivada_del_flector_es_el_cortante(self):
        el, U = frame()
        s = el.line_load_samples(triangular_local, local=True)
        _, ef = clamped(el, U, s)
        xi = np.linspace(0.05, 0.95, 7)
        h = 1e-6
        al = el.internal_forces_along(ef, U, xi, s)
        dM = (el.internal_forces_along(ef, U, xi + h, s).components["M"]
              - el.internal_forces_along(ef, U, xi - h, s).components["M"]) / (2 * h * L)
        np.testing.assert_allclose(dM, al.components["V"], rtol=1e-6, atol=1e-3)


YAML_APOYADA = """
nodes:
  - {{id: 1, coords: [0.0, 0.0]}}
  - {{id: 2, coords: [{L}, 0.0]}}
materials:
  - {{id: 1, type: Elastic1D, E: 2.1e8}}
elements:
  - {{id: 1, type: Frame2DEuler, material: 1, nodes: [1, 2], A: 0.013, I: 2.1e-4}}
boundary_conditions_by_node:
  - {{node_id: 1, ux: 0.0, uy: 0.0}}
  - {{node_id: 2, uy: 0.0}}
{loads}
solver:
  type: LinearSolver
{output}
"""


def write(tmp_path, text, name="m.yaml"):
    p = tmp_path / name
    p.write_text(textwrap.dedent(text), encoding="utf-8")
    return p


class TestApoyadaYaml:

    def test_triangular_extremo_en_L_raiz_de_3(self, tmp_path):
        p = write(tmp_path, YAML_APOYADA.format(
            L=L, loads=f"line_loads:\n  - {{elements: {{id: 1}}, q: [0.0, \"-{Q0}*x/{L}\"]}}",
            output=""))
        m = load_yaml(p)
        r = m.solve()
        el = m.domain.elements[1]
        xs = 1 / np.sqrt(3)
        al = r.internal_forces_along(el, [0.0, xs, 1.0])
        np.testing.assert_allclose(al.components["V"], [Q0 * L / 6, 0.0, -Q0 * L / 3],
                                   rtol=1e-9, atol=1e-6)
        assert al.components["M"][1] == pytest.approx(Q0 * L**2 / (9 * np.sqrt(3)), rel=1e-10)

    def test_diagrama_historia_y_vtk_dicen_lo_mismo(self, tmp_path):
        import meshio
        out = textwrap.dedent("""
            output:
              results: {diagrams: {stations: 5}}
              history:
                columns:
                  - {name: M_centro, element: {id: 1, var: M, at: 0.5}}
                  - {name: V_i, element: {id: 1, var: V, end: i}}
            """)
        p = write(tmp_path, YAML_APOYADA.format(
            L=L, loads=f"line_loads:\n  - {{elements: {{id: 1}}, q: [0.0, -{W}]}}", output=out))
        load_yaml(p).solve(write_output=True)
        hist = np.genfromtxt(tmp_path / "m_historia.csv", delimiter=",", names=True)
        assert float(hist["M_centro"]) == pytest.approx(W * L**2 / 8, rel=1e-10)
        assert float(hist["V_i"]) == pytest.approx(W * L / 2, rel=1e-10)
        d = np.genfromtxt(tmp_path / "m_diagramas_final.csv", delimiter=",", names=True)
        assert d["V"][0] == pytest.approx(W * L / 2, rel=1e-10)
        assert np.max(d["M"]) == pytest.approx(W * L**2 / 8, rel=1e-10)
        vtk = meshio.read(tmp_path / "m_final.vtu")
        assert float(vtk.cell_data["V_i"][0][0]) == pytest.approx(W * L / 2, rel=1e-10)


class TestArchivoDeDiagramas:

    def test_incluye_el_cero_del_cortante(self, tmp_path):
        """Con carga triangular el flector máximo cae en ξ = 1/√3, entre las
        estaciones uniformes: el archivo lo incluye exacto."""
        p = write(tmp_path, YAML_APOYADA.format(
            L=L, loads=f"line_loads:\n  - {{elements: {{id: 1}}, q: [0.0, \"-{Q0}*x/{L}\"]}}",
            output="output:\n  results: {diagrams: {stations: 4}}\n"))
        load_yaml(p).solve(write_output=True)
        d = np.genfromtxt(tmp_path / "m_diagramas_final.csv", delimiter=",", names=True)
        assert len(d) == 5
        k = int(np.argmax(d["M"]))
        assert d["xi"][k] == pytest.approx(1 / np.sqrt(3), abs=1e-9)
        assert d["M"][k] == pytest.approx(Q0 * L**2 / (9 * np.sqrt(3)), rel=1e-9)
        assert abs(d["V"][k]) < 1e-6 * Q0 * L
        np.testing.assert_allclose(d["x"], d["xi"] * L, atol=1e-12)

    def test_cero_del_cortante_en_una_estacion_uniforme_sin_fila_repetida(self, tmp_path):
        """Con carga uniforme y estaciones impares el cero del cortante cae
        en ξ = 1/2, que ya es estación: una fila, no dos (validación del
        2026-11-27). Y el axil nulo de la viga horizontal se escribe 0, no
        «-0»."""
        p = write(tmp_path, YAML_APOYADA.format(
            L=L, loads=f"line_loads:\n  - {{elements: {{id: 1}}, q: [0.0, -{W}]}}",
            output="output:\n  results: {diagrams: {stations: 5}}\n"))
        load_yaml(p).solve(write_output=True)
        text = (tmp_path / "m_diagramas_final.csv").read_text(encoding="utf-8")
        d = np.genfromtxt(tmp_path / "m_diagramas_final.csv", delimiter=",", names=True)
        assert len(d) == 5
        np.testing.assert_allclose(d["xi"], np.linspace(0.0, 1.0, 5), atol=1e-12)
        assert d["M"][2] == pytest.approx(W * L**2 / 8, rel=1e-10)
        assert ",-0," not in text and not any(line.endswith(",-0") for line in text.splitlines())


YAML_CONTINUA = """
nodes:
  - {{id: 1, coords: [0.0, 0.0]}}
  - {{id: 2, coords: [5.0, 0.0]}}
  - {{id: 3, coords: [12.0, 0.0]}}
{extra_nodes}
materials:
  - {{id: 1, type: Elastic1D, E: 2.1e8}}
elements:
{elements}
boundary_conditions_by_node:
  - {{node_id: 1, ux: 0.0, uy: 0.0}}
  - {{node_id: 2, uy: 0.0}}
  - {{node_id: 3, uy: 0.0}}
line_loads:
  - {{elements: {{ids: {ids}}}, q: [0.0, -10.0e3]}}
solver:
  type: LinearSolver
output:
  results: {{diagrams: true}}
"""


def continua(tmp_path, n_por_metro=False):
    """Viga continua de 5 + 7 m con 10 kN/m (caso v1/v2 de la validación del
    2026-11-27): M_B = −q(L1³ + L2³)/[8(L1 + L2)] = −48750, máximos de claro
    11628.125 y 39300.064 (ecuación de los tres momentos)."""
    if not n_por_metro:
        extra, ids = "", [1, 2]
        elements = ("  - {id: 1, type: Frame2DEuler, material: 1, nodes: [1, 2], A: 0.013, I: 2.1e-4}\n"
                    "  - {id: 2, type: Frame2DEuler, material: 1, nodes: [2, 3], A: 0.013, I: 2.1e-4}")
    else:
        xs = [x for x in range(13) if x not in (0, 5, 12)]
        extra = "\n".join(f"  - {{id: {10 + x}, coords: [{float(x)}, 0.0]}}" for x in xs)
        nid = {0: 1, 5: 2, 12: 3, **{x: 10 + x for x in xs}}
        ids = list(range(1, 13))
        elements = "\n".join(
            f"  - {{id: {x + 1}, type: Frame2DEuler, material: 1, nodes: [{nid[x]}, {nid[x + 1]}],"
            " A: 0.013, I: 2.1e-4}" for x in range(12))
    p = write(tmp_path, YAML_CONTINUA.format(extra_nodes=extra, elements=elements, ids=ids))
    load_yaml(p).solve(write_output=True)
    return p


def labels_of_diagram(p, var):
    """Dibuja la figura ``diagram`` de ``var`` y devuelve sus rótulos (texto y
    caja en pantalla)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from solidum.post import plots as P
    from solidum.post import ResultFiles
    created = []
    real = plt.subplots

    def spy(*a, **k):
        fig, axes = real(*a, **k)
        created.append(axes)
        return fig, axes
    try:
        plt.subplots = spy
        P.plot_diagram(ResultFiles.of_yaml(p), {"type": "diagram", "var": var},
                       p.parent / f"d_{var}.png")
    finally:
        plt.subplots = real
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    ax = created[0][0, 0]
    canvas = FigureCanvasAgg(ax.figure)   # la figura cerrada perdió su lienzo
    canvas.draw()
    renderer = canvas.get_renderer()
    return [(t.get_text(), t.get_window_extent(renderer)) for t in ax.texts]


class TestFiguraDiagram:

    def test_rotula_el_maximo_del_claro_aunque_sea_menor_que_el_del_apoyo(self, tmp_path):
        """Manual de usuario (figura ``diagram``): se rotula el flector máximo
        del claro, donde el cortante se anula. Con un elemento por claro de la
        viga continua, los máximos 11628 y 39300 son menores que el flector
        del apoyo (48750) y antes no se rotulaban (validación 2026-11-27)."""
        values = [float(t) for t, _ in labels_of_diagram(continua(tmp_path), "M")]
        for expected in (48750.0, 11628.125, 39300.064):
            assert any(v == pytest.approx(expected, rel=1e-4) for v in values), (expected, values)

    def test_rotula_el_extremo_de_menor_valor_absoluto(self, tmp_path):
        """Marco portal de la validación del 2026-11-27: en la columna con
        viento el flector es negativo en todo el tramo y su extremo donde el
        cortante se anula (−334.64 en s = 0.632 m) es el menor valor absoluto
        del tramo: también se rotula (extremo del flector, no de |M|)."""
        text = textwrap.dedent("""
            nodes:
              - {id: 1, coords: [0.0, 0.0]}
              - {id: 2, coords: [0.0, 4.0]}
              - {id: 3, coords: [3.0, 4.0]}
              - {id: 4, coords: [6.0, 4.0]}
              - {id: 5, coords: [6.0, 0.0]}
            materials:
              - {id: 1, type: Elastic1D, E: 25.0e9}
            elements:
              - {id: 1, type: Frame2DEuler, material: 1, nodes: [1, 2], A: 0.09, I: 6.75e-4}
              - {id: 2, type: Frame2DEuler, material: 1, nodes: [2, 3], A: 0.12, I: 1.6e-3}
              - {id: 3, type: Frame2DEuler, material: 1, nodes: [4, 3], A: 0.12, I: 1.6e-3}
              - {id: 4, type: Frame2DEuler, material: 1, nodes: [5, 4], A: 0.09, I: 6.75e-4}
            boundary_conditions_by_node:
              - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}
              - {node_id: 5, ux: 0.0, uy: 0.0, rz: 0.0}
            line_loads:
              - {elements: {ids: [2, 3]}, q: [0.0, -12.0e3]}
              - {elements: {id: 1}, q: [3.0e3, 0.0]}
            solver:
              type: LinearSolver
            output:
              results: {diagrams: true}
            """)
        p = write(tmp_path, text)
        load_yaml(p).solve(write_output=True)
        values = [float(t) for t, _ in labels_of_diagram(p, "M")]
        for expected in (334.641, 933.879, 33095.640):
            assert any(v == pytest.approx(expected, rel=1e-3) for v in values), (expected, values)

    def test_rotulos_sin_encimarse(self, tmp_path):
        """Con doce elementos de 1 m los rótulos de nudos vecinos y del máximo
        del claro se encimaban (validación 2026-11-27): ninguna caja toca a
        otra, y el máximo del claro, el mayor de su zona, queda."""
        for var in ("M", "V"):
            labels = labels_of_diagram(continua(tmp_path, n_por_metro=True), var)
            assert len(labels) >= 3
            boxes = [b for _, b in labels]
            for i in range(len(boxes)):
                for j in range(i + 1, len(boxes)):
                    assert not boxes[i].overlaps(boxes[j]), (labels[i][0], labels[j][0])
            if var == "M":
                values = [float(t) for t, _ in labels]
                assert any(v == pytest.approx(11628.125, rel=1e-4) for v in values), values

    def test_flector_del_lado_de_la_traccion_sin_importar_la_numeracion(self):
        """Q1 del usuario (2026-11-27): el flector se dibuja hacia la fibra
        en tracción; N y V, hacia unos ejes de dibujo fijos. Permutar los
        nodos invierte t, y' y el signo de M: el dibujo no cambia."""
        from solidum.post.plots import _drawing_normal
        t = np.array([0.6, 0.8])
        # M > 0 tracciona y' < 0: se dibuja hacia −y'.
        np.testing.assert_allclose(_drawing_normal(t, "M"), [0.8, -0.6])
        for var, sign in (("M", -1.0), ("V", 1.0), ("N", 1.0)):
            # Valor del elemento permutado: −M, V y N iguales (spec §15).
            a = 1.0 * _drawing_normal(t, var)
            b = sign * _drawing_normal(-t, var)
            np.testing.assert_allclose(a, b, atol=1e-15)
        # Barra horizontal: V positivo hacia arriba; vertical: hacia la izquierda.
        np.testing.assert_allclose(_drawing_normal(np.array([-1.0, 0.0]), "V"), [0.0, 1.0])
        np.testing.assert_allclose(_drawing_normal(np.array([0.0, -1.0]), "V"), [-1.0, 0.0])

    def test_figura_desde_el_yaml_y_validacion(self, tmp_path):
        from solidum.post.plots import validate_plots_block
        assert validate_plots_block([{"type": "diagram"}])
        assert validate_plots_block([{"type": "diagram", "var": "M", "scale": -1}])
        assert not validate_plots_block([{"type": "diagram", "var": "V", "at": ["final"]}])
        out = "output:\n  plots:\n    - {type: diagram, var: M, file: dm.png}\n"
        p = write(tmp_path, YAML_APOYADA.format(
            L=L, loads=f"line_loads:\n  - {{elements: {{id: 1}}, q: [0.0, -{W}]}}", output=out))
        load_yaml(p).solve(write_output=True)
        assert (tmp_path / "dm.png").stat().st_size > 1000
        assert (tmp_path / "m_diagramas_final.csv").is_file()
        assert (tmp_path / "m_diagramas_indice.csv").is_file()

    def test_diagramas_rechazados_en_el_tiempo(self, tmp_path):
        text = YAML_APOYADA.format(L=L, loads="", output="output:\n  results: {diagrams: true}\n")
        text = text.replace("materials:\n  - {id: 1, type: Elastic1D, E: 2.1e8}",
                            "materials:\n  - {id: 1, type: Elastic1D, E: 2.1e8, density: 7850.0}")
        text = text.replace("solver:\n  type: LinearSolver",
                            "solver:\n  type: NewmarkSolver\n  dt: 1.0e-3\n  t_end: 2.0e-3")
        with pytest.raises(ValueError, match="sólo se escriben en análisis estáticos"):
            load_yaml(write(tmp_path, text)).solve(write_output=True)


class TestVtkDescuentaLaCarga:

    def test_voladizo_con_fuerza_de_cuerpo(self, tmp_path):
        """El VTK escribía las fuerzas nodales K·u (−6018.3 en el caso del
        arquitecto) en vez de las de extremo; ahora las del SolveResult."""
        import meshio
        text = textwrap.dedent(f"""
            nodes:
              - {{id: 1, coords: [0.0, 0.0]}}
              - {{id: 2, coords: [{L / 2}, 0.0]}}
              - {{id: 3, coords: [{L}, 0.0]}}
            materials:
              - {{id: 1, type: Elastic1D, E: 2.1e8}}
            elements:
              - {{id: 1, type: Frame2DEuler, material: 1, nodes: [1, 2], A: {A}, I: {I}}}
              - {{id: 2, type: Frame2DEuler, material: 1, nodes: [2, 3], A: {A}, I: {I}}}
            boundary_conditions_by_node:
              - {{node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}}
            body_force: [0.0, {-W / A}]
            solver:
              type: LinearSolver
            output:
              results: {{frequency: last_step}}
            """)
        p = write(tmp_path, text)
        r = load_yaml(p).solve(write_output=True)
        vtk = meshio.read(tmp_path / "m_final.vtu")
        assert float(vtk.cell_data["M_i"][0][0]) == pytest.approx(-W * L**2 / 2, rel=1e-10)
        assert float(r.element_forces[1].components["M"][0]) == pytest.approx(-W * L**2 / 2,
                                                                                rel=1e-10)


# ---------------------------------------------------------------------------
# Armaduras (spec Truss2D §12)
# ---------------------------------------------------------------------------

class TestArmaduras:

    def test_barra_colgada_y_carga_transversal(self):
        dom = Domain()
        n1, n2 = dom.add_node(1, [0.0, 0.0]), dom.add_node(2, [0.0, -L])
        el = Truss2D(1, [n1, n2], Elastic1D(E=E), A=A)
        dom.add_element(el)
        dom.generate_equation_numbers()
        U = np.zeros(dom.total_dofs)
        s = el.line_load_samples([0.0, -W])          # a lo largo del eje, hacia abajo
        f = el.line_load_vector(s)
        np.testing.assert_allclose(f, [0, -W * L / 2, 0, -W * L / 2], rtol=1e-12)
        # Colgada del nodo 1: el nodo 1 le da +W·L hacia arriba, el 2 nada.
        P = np.array([0.0, W * L, 0.0, 0.0])
        ef = el.internal_forces_from_nodal(P + f, U, equivalent_load=f)
        al = el.internal_forces_along(ef, U, [0.0, 0.5, 1.0], s)
        np.testing.assert_allclose(al.components["N"], [W * L, W * L / 2, 0.0], atol=1e-8)
        tri = el.line_load_samples(lambda x: [Q0 * (-x[1]) / L, 0.0])
        g = el.line_load_vector(tri)
        np.testing.assert_allclose(g[[0, 2]], [Q0 * L / 6, Q0 * L / 3], rtol=1e-12)


class TestCorrotacional:
    """Giro rígido de 90° (F_int = 0, λ = 1, cos φ = 0), spec Frame2DEulerCorot
    §11: con carga uniforme el diagrama cierra; con carga lineal queda
    ΔM = Δq⊥·L0²·(λ·cos φ − 1)/60 = −Δq⊥·L0²/60 (caso 11 del teórico)."""

    def setup_method(self):
        from solidum.elements.frame import Frame2DEulerCorot
        dom = Domain()
        n1, n2 = dom.add_node(1, [0.0, 0.0]), dom.add_node(2, [L, 0.0])
        self.el = Frame2DEulerCorot(1, [n1, n2], Elastic1D(E=E), A=A, I=I)
        dom.add_element(self.el)
        dom.generate_equation_numbers()
        U = np.zeros(dom.total_dofs)
        U[n2.dofs["ux"]], U[n2.dofs["uy"]] = -L, L
        U[n1.dofs["rz"]] = U[n2.dofs["rz"]] = np.pi / 2
        self.U = U

    def diagram(self, q):
        s = self.el.line_load_samples(q)
        f = self.el.line_load_vector(s)
        ef = self.el.internal_forces(self.U, equivalent_load=f)
        return ef, self.el.internal_forces_along(ef, self.U, [0.0, 1.0], s)

    def test_sin_fuerza_interna_y_cierre_con_carga_uniforme(self):
        _, F_int = self.el._evaluate_without_trace(self.el.get_local_displacements(self.U))
        assert np.max(np.abs(F_int)) < 1e-6
        ef, al = self.diagram([0.0, -W])
        for c in ("N", "V", "M"):
            assert al.components[c][1] == pytest.approx(ef.components[c][1], abs=1e-9 * W * L**2)

    def test_cierre_exacto_con_giro_y_alargamiento(self):
        """Giro de 30°, alargamiento del 1 % (λ ≠ 1) y flexión (giros nodales
        distintos del de la cuerda, cortante en F_int): con carga uniforme el
        cierre es exacto si los brazos se miden con la longitud corriente."""
        from solidum.elements.frame import Frame2DEulerCorot
        dom = Domain()
        n1, n2 = dom.add_node(1, [0.0, 0.0]), dom.add_node(2, [L, 0.0])
        el = Frame2DEulerCorot(1, [n1, n2], Elastic1D(E=E), A=A, I=I)
        dom.add_element(el)
        dom.generate_equation_numbers()
        phi, lam = np.deg2rad(30.0), 1.01
        U = np.zeros(dom.total_dofs)
        U[n2.dofs["ux"]] = L * lam * np.cos(phi) - L
        U[n2.dofs["uy"]] = L * lam * np.sin(phi)
        U[n1.dofs["rz"]], U[n2.dofs["rz"]] = phi + 0.01, phi - 0.02
        s = el.line_load_samples([0.0, -W])
        f = el.line_load_vector(s)
        ef = el.internal_forces(U, equivalent_load=f)
        _, F_int = el._evaluate_without_trace(el.get_local_displacements(U))
        assert max(abs(F_int[2]), abs(F_int[5])) > 1e-3 * W * L**2   # hay flexión en F_int
        al = el.internal_forces_along(ef, U, [1.0], s)
        scale = max(abs(ef.components["M"][0]), W * L**2)
        assert al.components["M"][0] == pytest.approx(ef.components["M"][1], abs=1e-10 * scale)

    def test_error_de_cierre_con_carga_lineal(self):
        ef, al = self.diagram(lambda x: [0.0, -Q0 * x[0] / L])
        dq = -Q0                                  # variación de q⊥ inicial de i a j
        assert al.components["M"][1] - ef.components["M"][1] == pytest.approx(
            -dq * L**2 / 60, rel=1e-6)
        for c in ("N", "V"):
            assert al.components[c][1] == pytest.approx(ef.components[c][1], abs=1e-8 * Q0 * L)


class TestOrientacion:

    def test_permutar_los_nodos_cambia_el_signo_del_flector(self, tmp_path):
        base = YAML_APOYADA.format(
            L=L, loads=f"line_loads:\n  - {{elements: {{id: 1}}, q: [0.0, \"-{Q0}*x/{L}\"]}}",
            output="")
        rev = base.replace("nodes: [1, 2], A:", "nodes: [2, 1], A:")
        m2 = load_yaml(write(tmp_path, rev, "r.yaml"))
        r2 = m2.solve()
        m1 = load_yaml(write(tmp_path, base))
        r1 = m1.solve()
        xi = np.linspace(0.0, 1.0, 7)
        a = r1.internal_forces_along(m1.domain.elements[1], xi).components
        b = r2.internal_forces_along(m2.domain.elements[1], 1.0 - xi).components
        np.testing.assert_allclose(b["M"], -a["M"], rtol=1e-10, atol=1e-6)
        np.testing.assert_allclose(b["V"], a["V"], rtol=1e-10, atol=1e-6)
        np.testing.assert_allclose(b["N"], a["N"], atol=1e-6)

    def test_voladizo_3d_flector_my(self, tmp_path):
        text = textwrap.dedent(f"""
            nodes:
              - {{id: 1, coords: [0.0, 0.0, 0.0]}}
              - {{id: 2, coords: [{L}, 0.0, 0.0]}}
            materials:
              - {{id: 1, type: Elastic1D, E: 2.1e8}}
            elements:
              - {{id: 1, type: Frame3D, material: 1, nodes: [1, 2], A: {A}, Iy: {I},
                 Iz: {1.7 * I}, J: {0.9 * I}, nu: 0.3}}
            boundary_conditions_by_node:
              - {{node_id: 1, ux: 0.0, uy: 0.0, uz: 0.0, rx: 0.0, ry: 0.0, rz: 0.0}}
            line_loads:
              - {{elements: {{id: 1}}, q_local: [0.0, {-W}, {-Q0}]}}
            solver:
              type: LinearSolver
            """)
        m = load_yaml(write(tmp_path, text))
        r = m.solve()
        xi = np.linspace(0.0, 1.0, 5)
        c = r.internal_forces_along(m.domain.elements[1], xi).components
        s = xi * L
        np.testing.assert_allclose(c["My"], Q0 * (L - s) ** 2 / 2, rtol=1e-9, atol=1e-6)
        np.testing.assert_allclose(c["Vz"], -Q0 * (L - s), rtol=1e-9, atol=1e-6)
        # Plano xy: igual que el 2D (M_2D = Mz, V_2D = −Vy); voladizo con q_y' = −W.
        np.testing.assert_allclose(c["Mz"], -W * (L - s) ** 2 / 2, rtol=1e-9, atol=1e-6)
        np.testing.assert_allclose(c["Vy"], -W * (L - s), rtol=1e-9, atol=1e-6)


# ---------------------------------------------------------------------------
# Casos de la verificación del teórico (2026-11-27)
# ---------------------------------------------------------------------------

def _section(mat2):
    from solidum.sections import FiberSection, RectangularPatch
    return FiberSection(patches=[RectangularPatch(mat2, y=(-0.25, 0.25), z=(-0.15, 0.15),
                                                  ny=8, nz=4)], GJ=3.1e6)


def _line_elements():
    """Una instancia de cada elemento de línea registrado, oblicua, con
    material y, si lo admite, con sección por fibras (ρA ≠ ρ·A)."""
    from solidum.elements.frame import Frame2DEulerCorot
    from solidum.elements.truss import Truss2DCorot, Truss3D, Truss3DCorot
    from solidum.elements.truss_quadratic import Truss2D3, Truss3D3
    mat = Elastic1D(E=E, density=7850.0)
    mat2 = Elastic1D(E=0.3 * E, density=2400.0)
    a2, b2 = np.array([0.4, -0.2]), np.array([0.4, -0.2]) + L * T
    a3, b3 = np.array([0.3, -1.1, 0.7]), np.array([1.4, 0.9, 1.4])
    out = []

    def build(cls, coords, **kw):
        dom = Domain()
        nodes = [dom.add_node(k + 1, list(c)) for k, c in enumerate(coords)]
        el = cls(1, nodes, **kw)
        dom.add_element(el)
        dom.generate_equation_numbers()
        out.append(el)
    for cls in (Frame2DEuler, Frame2DEulerCorot):
        build(cls, [a2, b2], material=mat, A=A, I=I)
        build(cls, [a2, b2], section=_section(mat2))
    build(Frame2DTimoshenko, [a2, b2], material=mat, A=A, I=I, As=0.011, nu=0.3)
    build(Frame3D, [a3, b3], material=mat, A=A, Iy=I, Iz=1.7 * I, J=0.9 * I, nu=0.3,
          ref_vector=[0.2, 0.1, 1.0])
    build(Frame3D, [a3, b3], section=_section(mat2), ref_vector=[0.2, 0.1, 1.0])
    for cls in (Truss2D, Truss2DCorot):
        build(cls, [a2, b2], material=mat, A=A)
    for cls in (Truss3D, Truss3DCorot):
        build(cls, [a3, b3], material=mat, A=A)
    build(Truss2D3, [a2, b2, (a2 + b2) / 2], material=mat, A=A)
    build(Truss3D3, [a3, b3, (a3 + b3) / 2], material=mat, A=A)
    return out


class TestBarridoDeElementos:

    def test_todas_las_clases_de_linea_estan_en_el_barrido(self):
        assert {type(el) for el in _line_elements()} == set(_line_classes())

    @pytest.mark.parametrize("k", range(len(_line_elements())))
    def test_vector_igual_a_muestras_con_material_y_seccion(self, k):
        el = _line_elements()[k]
        d = 3 if "uz" in el.DOF_NAMES else 2
        b = np.array([1.3, -2.7, 0.9])[:d]
        field = (lambda x: [0.5 * x[0], -1.0 - x[1]]) if d == 2 else (
            lambda x: [0.5 * x[0], -1.0 - x[1], 0.3 * x[2]])
        g = np.array([0.0, -9.81, 0.0])[:d] if d == 2 else np.array([0.0, 0.0, -9.81])
        for load in (b, field):
            ref = el.compute_body_load(load)
            np.testing.assert_allclose(el.line_load_vector(el.body_load_samples(load)), ref,
                                       rtol=1e-13, atol=1e-12 * np.max(np.abs(ref)))
        ref = el.compute_self_weight(g)
        np.testing.assert_allclose(el.line_load_vector(el.self_weight_samples(g)), ref,
                                   rtol=1e-13, atol=1e-12 * np.max(np.abs(ref)))

    @pytest.mark.parametrize("k", range(len(_line_elements())))
    def test_extremos_y_derivadas_con_desplazamientos_al_azar(self, k):
        el = _line_elements()[k]
        if not el.INTERNAL_FORCES_ALONG:
            pytest.skip("sin fuerzas internas a lo largo del elemento")
        d = 3 if "uz" in el.DOF_NAMES else 2
        rng = np.random.default_rng(7 + k)
        n = len(el.nodes[0].dofs)
        U = np.zeros(n * len(el.nodes))
        U[:] = 1e-3 * rng.standard_normal(U.size)
        # En los corrotacionales el cierre del flector es exacto sólo con carga
        # uniforme (spec Frame2DEulerCorot §11: ΔM = Δq⊥·L0²·(λcosφ − 1)/60).
        curved = 0.0 if getattr(el, "GEOMETRICALLY_NONLINEAR", False) else 0.3
        q = lambda x: list(1e3 * (np.array([0.7, -1.3, 0.4])[:d]           # noqa: E731
                                  + curved * np.asarray(x)[:d] ** 2))
        s = el.line_load_samples(q)
        f = el.line_load_vector(s)
        ef = el.internal_forces(U, equivalent_load=f)
        al = el.internal_forces_along(ef, U, [0.0, 1.0], s)
        for c, v in ef.components.items():
            scale = 1e-9 * max(1.0, np.max(np.abs(v)), np.max(np.abs(al.components[c])))
            np.testing.assert_allclose(al.components[c], v, atol=scale, err_msg=c)
        pairs = {"frame2d": (("M", "V", 1.0),), "frame3d": (("Mz", "Vy", -1.0), ("My", "Vz", 1.0))}
        R, L0, l = el.internal_force_axes(U)
        xi = np.linspace(0.1, 0.9, 5)
        h = 1e-6
        for moment, shear, sign in pairs.get(ef.kind, ()):
            plus = el.internal_forces_along(ef, U, xi + h, s).components[moment]
            minus = el.internal_forces_along(ef, U, xi - h, s).components[moment]
            dM = (plus - minus) / (2 * h * l)
            V = el.internal_forces_along(ef, U, xi, s).components[shear]
            np.testing.assert_allclose(dM, sign * V, rtol=1e-5, atol=1e-6 * np.max(np.abs(V)))


class TestCasosDelTeorico:

    def test_timoshenko_biempotrada_triangular(self):
        el, U = frame(Frame2DTimoshenko, As=0.011)
        phi, k = el._phi(), 1.0 + el._phi()
        s = el.line_load_samples(triangular_local, local=True)
        _, ef = clamped(el, U, s)
        al = el.internal_forces_along(ef, U, [0.0, 1.0], s)
        np.testing.assert_allclose(al.components["M"],
                                   [-Q0 * L**2 * (4 + 5 * phi) / (120 * k),
                                    -Q0 * L**2 * (6 + 5 * phi) / (120 * k)], rtol=1e-11)
        assert al.components["V"][0] == pytest.approx(Q0 * L * (9 + 10 * phi) / (60 * k), rel=1e-11)

    def test_voladizo_triangular_2d(self, tmp_path):
        text = YAML_APOYADA.format(
            L=L, loads=f"line_loads:\n  - {{elements: {{id: 1}}, q: [0.0, \"-{Q0}*x/{L}\"]}}",
            output="").replace("  - {node_id: 1, ux: 0.0, uy: 0.0}\n  - {node_id: 2, uy: 0.0}",
                               "  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}")
        m = load_yaml(write(tmp_path, text))
        r = m.solve()
        xi = np.linspace(0, 1, 6)
        sx = xi * L
        c = r.internal_forces_along(m.domain.elements[1], xi).components
        np.testing.assert_allclose(c["M"], -Q0 * (2 * L**3 - 3 * L**2 * sx + sx**3) / (6 * L),
                                   rtol=1e-9, atol=1e-6)
        assert c["V"][0] == pytest.approx(Q0 * L / 2, rel=1e-9)

    def test_proyectada_en_un_elemento_paralelo_es_nula(self):
        dom = Domain()
        n1, n2 = dom.add_node(1, [0.2, 0.0]), dom.add_node(2, [0.2, L])
        el = Frame2DEuler(1, [n1, n2], Elastic1D(E=E), A=A, I=I)
        assert np.max(np.abs(el.compute_line_load([0.0, -W], projected=True))) == 0.0

    def test_corrotacional_cierre_de_n_y_v_con_alargamiento(self):
        from solidum.elements.frame import Frame2DEulerCorot
        dom = Domain()
        n1, n2 = dom.add_node(1, [0.0, 0.0]), dom.add_node(2, [L, 0.0])
        el = Frame2DEulerCorot(1, [n1, n2], Elastic1D(E=E), A=A, I=I)
        dom.add_element(el)
        dom.generate_equation_numbers()
        phi, lam = np.deg2rad(52.0), 1.013
        U = np.zeros(dom.total_dofs)
        U[n2.dofs["ux"]], U[n2.dofs["uy"]] = L * lam * np.cos(phi) - L, L * lam * np.sin(phi)
        U[n1.dofs["rz"]], U[n2.dofs["rz"]] = phi + 0.02, phi - 0.03
        s = el.line_load_samples(lambda x: [0.3 * W, -W * (1 + x[0] / L)])
        f = el.line_load_vector(s)
        ef = el.internal_forces(U, equivalent_load=f)
        al = el.internal_forces_along(ef, U, [1.0], s)
        for c in ("N", "V"):
            assert al.components[c][0] == pytest.approx(ef.components[c][1], abs=1e-9 * W * L)

    def test_armadura_corrotacional_con_giro_grande(self):
        from solidum.elements.truss import Truss2DCorot
        dom = Domain()
        n1, n2 = dom.add_node(1, [0.0, 0.0]), dom.add_node(2, [L, 0.0])
        el = Truss2DCorot(1, [n1, n2], Elastic1D(E=E), A=A)
        dom.add_element(el)
        dom.generate_equation_numbers()
        U = np.zeros(dom.total_dofs)
        phi = np.deg2rad(70.0)
        U[n2.dofs["ux"]], U[n2.dofs["uy"]] = 1.01 * L * np.cos(phi) - L, 1.01 * L * np.sin(phi)
        s = el.line_load_samples([0.0, -W])
        f = el.line_load_vector(s)
        ef = el.internal_forces(U, equivalent_load=f)
        R, _, l = el.internal_force_axes(U)
        np.testing.assert_allclose(R[0], [np.cos(phi), np.sin(phi)], atol=1e-14)
        assert l == pytest.approx(1.01 * L, rel=1e-12)
        al = el.internal_forces_along(ef, U, [1.0], s)
        assert al.components["N"][0] == pytest.approx(ef.components["N"][1], abs=1e-9 * W * L)

    def test_escalado_y_unidades(self, tmp_path):
        """N, m frente a kN, mm: M en N·m = kN·mm; N y V en N = 1e3·kN. Y la
        respuesta es lineal en la carga (Reglas §6)."""
        def model(scale_len, scale_force, factor):
            Lm = 4.2 * scale_len
            Em = 2.1e11 * scale_force / scale_len**2
            q = factor * 7.3e3 * scale_force / scale_len
            return textwrap.dedent(f"""
                nodes:
                  - {{id: 1, coords: [0.0, 0.0]}}
                  - {{id: 2, coords: [{0.6 * Lm}, {0.8 * Lm}]}}
                  - {{id: 3, coords: [{1.6 * Lm}, {0.8 * Lm}]}}
                materials:
                  - {{id: 1, type: Elastic1D, E: {Em}}}
                elements:
                  - {{id: 1, type: Frame2DEuler, material: 1, nodes: [1, 2], A: {0.013 * scale_len**2}, I: {2.1e-4 * scale_len**4}}}
                  - {{id: 2, type: Frame2DEuler, material: 1, nodes: [2, 3], A: {0.013 * scale_len**2}, I: {2.1e-4 * scale_len**4}}}
                boundary_conditions_by_node:
                  - {{node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}}
                  - {{node_id: 3, ux: 0.0, uy: 0.0}}
                line_loads:
                  - {{elements: {{id: 1}}, q: [0.0, {-q}], projected: true}}
                  - {{elements: {{id: 2}}, q_local: [0.0, "{-q}*(x - {0.6 * Lm})/{Lm}"]}}
                solver:
                  type: LinearSolver
                """)
        xi = np.linspace(0, 1, 7)

        def run(name, *args):
            m = load_yaml(write(tmp_path, model(*args), name))
            r = m.solve()
            return [r.internal_forces_along(m.domain.elements[e], xi).components for e in (1, 2)]
        a = run("si.yaml", 1.0, 1.0, 1.0)
        b = run("kn_mm.yaml", 1e3, 1e-3, 1.0)
        c = run("doble.yaml", 1.0, 1.0, 2.5)
        for e in range(2):
            Ms = np.max(np.abs(a[e]["M"]))
            np.testing.assert_allclose(b[e]["M"], a[e]["M"], rtol=1e-10, atol=1e-10 * Ms)
            for comp in ("N", "V"):
                np.testing.assert_allclose(1e3 * b[e][comp], a[e][comp], rtol=1e-10,
                                           atol=1e-10 * np.max(np.abs(a[e][comp])))
            for comp in ("N", "V", "M"):
                np.testing.assert_allclose(c[e][comp], 2.5 * a[e][comp], rtol=1e-10,
                                           atol=1e-10 * np.max(np.abs(a[e][comp])))

    def test_seccion_plastificada_isostatica(self, tmp_path):
        """Caso 10: la viga apoyada con FiberSection plastificada (qL²/8 =
        0.93·Mp) cumple la estática en todo x, también a medio cargar."""
        Lb, bw, hh, sy = 3.7, 0.15, 0.30, 250.0e6
        Mp = sy * bw * hh**2 / 4
        q = 0.93 * Mp * 8 / Lb**2
        nel = 6
        nodes = "\n".join(f"  - {{id: {k + 1}, coords: [{k * Lb / nel!r}, 0.0]}}"
                          for k in range(nel + 1))
        els = "\n".join(f"  - {{id: {k + 1}, type: Frame2DEuler, nodes: [{k + 1}, {k + 2}], "
                        f"section: 1, n_ip: 5}}" for k in range(nel))
        text = f"""
nodes:
{nodes}
materials:
  - {{id: 1, type: Elastoplastic1D, E: 2.0e11, sigma_y: {sy}, H: 2.0e9}}
sections:
  - id: 1
    type: FiberSection
    patches:
      - {{shape: rectangle, material: 1, y: [{-hh / 2}, {hh / 2}], z: [{-bw / 2}, {bw / 2}], ny: 40, nz: 1}}
elements:
{els}
boundary_conditions_by_node:
  - {{node_id: 1, ux: 0.0, uy: 0.0}}
  - {{node_id: {nel + 1}, uy: 0.0}}
line_loads:
  - {{elements: {{ids: {list(range(1, nel + 1))}}}, q: [0.0, {-q}]}}
solver:
  type: NonlinearSolver
  num_steps: 10
output:
  results: {{diagrams: {{stations: 5, at: [0.5, final]}}}}
"""
        p = tmp_path / "plast.yaml"
        p.write_text(text, encoding="utf-8")
        m = load_yaml(p)
        r = m.solve(write_output=True)
        mid = m.domain.nodes[nel // 2 + 1]
        EI = 2.0e11 * bw * hh**3 / 12
        assert r.U[mid.dofs["uy"]] / (-5 * q * Lb**4 / (384 * EI)) > 1.1     # plastificó
        for name, fac in (("plast_diagramas_final.csv", 1.0),):
            d = np.genfromtxt(tmp_path / name, delimiter=",", names=True)
            np.testing.assert_allclose(d["M"], fac * q * d["x"] * (Lb - d["x"]) / 2,
                                       rtol=1e-6, atol=1e-6 * q * Lb**2)
            np.testing.assert_allclose(d["V"], fac * q * (Lb / 2 - d["x"]), atol=1e-6 * q * Lb)
        idx = (tmp_path / "plast_diagramas_indice.csv").read_text(encoding="utf-8").splitlines()
        fac, fname = idx[1].split(",")
        d = np.genfromtxt(tmp_path / fname, delimiter=",", names=True)
        np.testing.assert_allclose(d["M"], float(fac) * q * d["x"] * (Lb - d["x"]) / 2,
                                   rtol=1e-6, atol=1e-6 * q * Lb**2)

    def test_archivo_de_diagramas_3d_con_cero_del_cortante(self, tmp_path):
        text = textwrap.dedent(f"""
            nodes:
              - {{id: 1, coords: [0.0, 0.0, 0.0]}}
              - {{id: 2, coords: [{L}, 0.0, 0.0]}}
            materials:
              - {{id: 1, type: Elastic1D, E: 2.1e8}}
            elements:
              - {{id: 1, type: Frame3D, material: 1, nodes: [1, 2], A: {A}, Iy: {I},
                 Iz: {1.7 * I}, J: {0.9 * I}, nu: 0.3}}
            boundary_conditions_by_node:
              - {{node_id: 1, ux: 0.0, uy: 0.0, uz: 0.0, rx: 0.0}}
              - {{node_id: 2, uy: 0.0, uz: 0.0}}
            line_loads:
              - {{elements: {{id: 1}}, q_local: [0.0, 0.0, "-{Q0}*x/{L}"]}}
            solver:
              type: LinearSolver
            output:
              results: {{diagrams: {{stations: 4}}}}
            """)
        load_yaml(write(tmp_path, text)).solve(write_output=True)
        d = np.genfromtxt(tmp_path / "m_diagramas_final.csv", delimiter=",", names=True)
        k = int(np.argmin(d["My"]))
        assert d["xi"][k] == pytest.approx(1 / np.sqrt(3), abs=1e-9)
        assert d["My"][k] == pytest.approx(-Q0 * L**2 / (9 * np.sqrt(3)), rel=1e-9)


# ---------------------------------------------------------------------------
# Hallazgos de la revisión del diff del arquitecto (2026-11-27)
# ---------------------------------------------------------------------------

class TestRevisionDelDiff:

    def _cantilever(self):
        from solidum.math.solvers import LinearSolver
        dom = Domain()
        nodes = [dom.add_node(k + 1, [k * L / 2, 0.0]) for k in range(3)]
        mat = Elastic1D(E=E)
        for k in range(2):
            dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], mat, A=A, I=I))
        for dof in ("ux", "uy", "rz"):
            nodes[0].fix_dof(dof, 0.0)
        dom.generate_equation_numbers()
        asm = Assembler(dom)
        return dom, asm, LinearSolver(asm)

    def test_asarray_da_error_en_vez_de_un_diagrama_sin_la_carga(self):
        """h. 1: la carga de línea pasa por el registro heredado: sin sus
        equivalentes (np.asarray), el resultado las recupera para las fuerzas
        de extremo y el diagrama da error en vez de salir sin la carga."""
        import solidum
        dom, asm, solver = self._cantilever()
        F = np.asarray(asm.assemble_line_loads([(el, [0.0, -W], False, False)
                                                for el in dom.elements.values()]))
        r = solidum.run(dom, assembler=asm, solver=solver, F_applied=F)
        assert float(r.element_forces[1].components["M"][0]) == pytest.approx(-W * L**2 / 2,
                                                                             rel=1e-10)
        with pytest.raises(ValueError, match="no lleva su carga a lo largo"):
            r.internal_forces_along(dom.elements[1], [0.5])

    def test_vtk_de_una_etapa_estatica_antes_de_un_dinamico(self, tmp_path):
        """h. 2: en un análisis mixto, el VTK de la etapa estática descuenta
        la carga; en los instantes en el tiempo con carga distribuida, NaN."""
        import meshio
        stages = textwrap.dedent(f"""
            stages:
              - {{name: carga, loads: {{line_loads: [{{elements: {{id: 1}}, q: [0.0, -{W}]}}]}},
                 solver: {{type: NonlinearSolver, num_steps: 1}}}}
              - {{name: vibra, loads: {{}}, solver: {{type: NewmarkSolver, dt: 1.0e-3, t_end: 2.0e-3}}}}
            output:
              results: {{frequency: all_steps}}
            """)
        text = YAML_APOYADA.format(L=L, loads="", output="").replace(
            "solver:\n  type: LinearSolver\n", stages).replace(
            "materials:\n  - {id: 1, type: Elastic1D, E: 2.1e8}",
            "materials:\n  - {id: 1, type: Elastic1D, E: 2.1e8, density: 7850.0}")
        r = load_yaml(write(tmp_path, text)).solve(write_output=True)
        first = r.results[0]
        vtk = meshio.read(tmp_path / "m_0001.vtu")
        assert float(vtk.cell_data["V_i"][0][0]) == pytest.approx(
            float(first.element_forces[1].components["V"][0]), rel=1e-10)
        assert float(vtk.cell_data["V_i"][0][0]) == pytest.approx(W * L / 2, rel=1e-10)
        later = meshio.read(tmp_path / "m_0002.vtu")
        assert np.isnan(float(later.cell_data["M_i"][0][0]))

    def test_la_cache_de_la_union_se_suelta(self):
        """h. 3: la caché de la fusión no retiene las disposiciones muertas."""
        import gc
        from solidum.core.load_vector import EquivalentLoadLayout
        a = EquivalentLoadLayout([1, 2, 3], [6, 6, 6])
        for _ in range(5):
            b = EquivalentLoadLayout([2, 4], [6, 6])
            a.union_with(b)
        del b
        gc.collect()
        assert len(a._union_cache) == 0

    def test_grupo_por_celdas_y_where(self):
        """h. 4: `group` toma las celdas exactas; una celda de un elemento que
        no admite la carga (un cable) da error, como con `ids`."""
        from solidum.elements.cable import Cable2DCorot
        from solidum.materials.cable_1d import CableMaterial1D
        from solidum.utils.yaml_loads import resolve_line_loads
        dom = Domain()
        n = [dom.add_node(k + 1, [float(k), 0.0]) for k in range(4)]
        mat = Elastic1D(E=E)
        dom.add_element(Frame2DEuler(1, [n[0], n[1]], mat, A=A, I=I))
        dom.add_element(Frame2DEuler(2, [n[1], n[2]], mat, A=A, I=I))
        dom.add_element(Cable2DCorot(3, [n[2], n[3]], CableMaterial1D(E=E), A=A))
        dom.physical_groups = {"Viga": [1, 2, 3], "Mixto": [2, 3, 4]}
        dom.physical_group_cells = {"Viga": [(1, 2), (2, 3)], "Mixto": [(2, 3), (3, 4)]}
        loads, errors = resolve_line_loads(dom, [{"elements": {"group": "Viga"}, "q": [0, -1.0]}])
        assert not errors and sorted(el.id for el, *_ in loads) == [1, 2]
        _, errors = resolve_line_loads(dom, [{"elements": {"group": "Mixto"}, "q": [0, -1.0]}])
        assert errors and "Cable2DCorot" in errors[0]
        _, errors = resolve_line_loads(dom, [{"elements": {"group": "Nada"}, "q": [0, -1.0]}])
        assert errors and "inexistente" in errors[0]
        loads, errors = resolve_line_loads(dom, [{"elements": {"where": "x < 1.5"}, "q": [0, -1.0]}])
        assert not errors and [el.id for el, *_ in loads] == [1]

    def test_muestras_del_resultado_con_sus_abscisas(self, tmp_path):
        """h. 5: SolveResult.line_load_samples lleva las abscisas y evalúa la
        carga representada en cualquier punto."""
        from solidum.core.load_vector import LineLoadSamples
        p = write(tmp_path, YAML_APOYADA.format(
            L=L, loads=f"line_loads:\n  - {{elements: {{id: 1}}, q: [0.0, \"-{Q0}*x/{L}\"]}}",
            output=""))
        r = load_yaml(p).solve()
        s = r.line_load_samples[1]
        assert isinstance(s, LineLoadSamples)
        np.testing.assert_allclose(s.at([0.0, 0.5, 1.0])[:, 1], [0.0, -Q0 / 2, -Q0],
                                   atol=1e-9 * Q0)

    def test_regresion_transitorio_con_fuerza_de_cuerpo(self, tmp_path):
        """h. 9: el rechazo del transitorio con carga distribuida, con la
        fuerza de cuerpo (con el código anterior escribía −286 en vez de −352,
        sin aviso)."""
        text = YAML_APOYADA.format(
            L=L, loads=f"body_force: [0.0, {-W / A}]",
            output="output:\n  history:\n    columns:\n"
                   "      - {name: M, element: {id: 1, var: M, end: i}}\n")
        text = text.replace("materials:\n  - {id: 1, type: Elastic1D, E: 2.1e8}",
                            "materials:\n  - {id: 1, type: Elastic1D, E: 2.1e8, density: 7850.0}")
        text = text.replace("solver:\n  type: LinearSolver",
                            "solver:\n  type: NewmarkSolver\n  dt: 1.0e-3\n  t_end: 2.0e-3")
        with pytest.raises(ValueError, match="sólo se registran en análisis estáticos"):
            load_yaml(write(tmp_path, text)).solve(write_output=True)

    def test_regresion_tracciones_de_una_etapa_al_leer(self, tmp_path):
        """h. 9: una tracción inválida de la etapa 2 se informa al leer, con
        su prefijo (con el código anterior, tras resolver la etapa 1)."""
        stage = textwrap.dedent("""
            stages:
              - {name: a, loads: {point_loads_by_node: [{node_id: 2, uy: -1.0}]},
                 solver: {type: LinearSolver}}
              - {name: b, loads: {tractions: [{element: 99, edge: 0, t: [1.0, 0.0]}]},
                 solver: {type: LinearSolver}}
            """)
        p = write(tmp_path, YAML_APOYADA.format(L=L, loads="", output="").replace(
            "solver:\n  type: LinearSolver\n", stage))
        with pytest.raises(YamlValidationError, match=r"stages\[1\]\.loads\.tractions"):
            load_yaml(p)


# ---------------------------------------------------------------------------
# Contrato y arquitectura (ADR 0044 §1–§3)
# ---------------------------------------------------------------------------

def _registered():
    return [ElementRegistry.get(n) for n in ElementRegistry.names()]


def _line_classes():
    return [c for c in _registered() if issubclass(c, LineElementMixin)]


class TestContrato:

    def test_todo_elemento_de_linea_declara_sus_capacidades(self):
        assert _line_classes()
        for cls in _line_classes():
            for flag in ("LINE_LOADS", "LINE_LOAD_LOCAL_AXES", "INTERNAL_FORCES_ALONG",
                         "LINE_LOAD_BY_COMPONENT"):
                assert isinstance(getattr(cls, flag), bool), (cls.__name__, flag)
            assert len(cls.LINE_LOAD_POINTS) == len(cls.LINE_LOAD_WEIGHTS)

    def test_los_demas_elementos_no_reciben_cargas_de_linea(self):
        for cls in _registered():
            if not issubclass(cls, LineElementMixin):
                assert not getattr(cls, "LINE_LOADS", False), cls.__name__

    def test_bit_a_bit_la_fuerza_de_cuerpo_como_campo_en_los_marcos(self):
        """Mismas operaciones que antes del ADR 0044 (expresión copiada)."""
        from solidum.elements._line import _LINE_W, _LINE_XI
        el, _ = frame()

        def b(x):
            return [3.0 * x[0], -7.0 + x[1] ** 2]
        vals = np.array([b(x) for x in el._sample_points()])
        R = el._load_axes()
        q = A * vals @ R.T
        N, _ = el._line_interpolation(_LINE_XI)
        old = el._to_global(np.einsum("g,gki,gk->i", _LINE_W * el.L0, N[:, :2, :], q))
        assert np.array_equal(old, el.compute_body_load(b))


class TestCanales:

    def setup_method(self):
        dom = Domain()
        nodes = [dom.add_node(k + 1, [k * 1.3, 0.2 * k]) for k in range(4)]
        mat = Elastic1D(E=E, density=7850.0)
        for k in range(3):
            dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], mat, A=A, I=I))
        dom.generate_equation_numbers()
        self.dom = dom
        self.asm = Assembler(dom)
        self.Fg = self.asm.assemble_self_weight([0.0, -9.81])
        els = list(dom.elements.values())
        self.Fl = self.asm.assemble_line_loads([(els[1], [0.0, -W], True, False)])

    def test_combinar_disposiciones_distintas(self):
        F = self.Fg + 2.5 * self.Fl
        eq = F.equivalent_loads
        for eid in (1, 2, 3):
            expected = self.Fg.equivalent_loads.samples_of(eid).copy()
            if eid == 2:
                expected += 2.5 * self.Fl.equivalent_loads.samples_of(2)
            np.testing.assert_allclose(eq.samples_of(eid), expected, rtol=1e-14)
            exp_f = self.Fg.equivalent_loads.get(eid) + (2.5 * self.Fl.equivalent_loads.get(2)
                                                         if eid == 2 else 0.0)
            np.testing.assert_allclose(eq.get(eid), exp_f, rtol=1e-14)

    def test_copia_serializacion_y_asarray(self):
        F = self.Fg + self.Fl
        for G in (copy.deepcopy(F), pickle.loads(pickle.dumps(F))):
            np.testing.assert_array_equal(G.equivalent_loads.samples_of(2),
                                          F.equivalent_loads.samples_of(2))
        assert not isinstance(np.asarray(F), LoadVector)


class TestYamlErrores:

    def test_errores_juntos_y_de_las_etapas_al_leer(self, tmp_path):
        loads = textwrap.dedent("""
            line_loads:
              - {elements: {id: 1}, q: [0.0, -1.0], q_local: [0.0, 1.0]}
              - {elements: {id: 1}, q_local: [0.0, -1.0], projected: true}
              - {elements: {id: 7}, q: [0.0, -1.0]}
              - {q: [0.0, -1.0]}
              - {elements: {id: 1}, q: [0.0, -1.0, 2.0]}
            """)
        p = write(tmp_path, YAML_APOYADA.format(L=L, loads=loads, output=""))
        with pytest.raises(YamlValidationError) as exc:
            load_yaml(p)
        assert len(exc.value.errors) == 5
        stage = textwrap.dedent("""
            stages:
              - {name: a, loads: {line_loads: [{elements: {id: 1}, q: [0.0, -1.0]}]},
                 solver: {type: LinearSolver}}
              - {name: b, loads: {line_loads: [{elements: {id: 9}, q: [0.0, -1.0]}]},
                 solver: {type: LinearSolver}}
            """)
        p = write(tmp_path, YAML_APOYADA.format(L=L, loads="", output="").replace(
            "solver:\n  type: LinearSolver\n", stage), "e.yaml")
        with pytest.raises(YamlValidationError, match=r"stages\[1\]\.loads\.line_loads"):
            load_yaml(p)

    def test_transitorio_rechaza_fuerzas_de_marco_con_carga(self, tmp_path):
        text = YAML_APOYADA.format(
            L=L, loads=f"line_loads:\n  - {{elements: {{id: 1}}, q: [0.0, -{W}]}}",
            output="output:\n  history:\n    columns:\n"
                   "      - {name: M, element: {id: 1, var: M, end: i}}\n")
        text = text.replace("materials:\n  - {id: 1, type: Elastic1D, E: 2.1e8}",
                            "materials:\n  - {id: 1, type: Elastic1D, E: 2.1e8, density: 7850.0}")
        text = text.replace("solver:\n  type: LinearSolver",
                            "solver:\n  type: NewmarkSolver\n  dt: 1.0e-3\n  t_end: 2.0e-3")
        p = write(tmp_path, text)
        with pytest.raises(ValueError, match="sólo se registran en análisis estáticos"):
            load_yaml(p).solve(write_output=True)
