"""Cargas distribuidas y variables en el YAML (ADR 0023; deudas #23 y #24).

- **Expresiones**: sólo números, las variables de la carga, aritmética, ``pi``
  y una lista de funciones; todo lo demás se rechaza al leer el archivo.
- **Igual que el API de Python**: las cargas declaradas en el YAML dan
  exactamente el vector que dan ``compute_edge_traction`` /
  ``compute_face_traction`` / ``compute_body_source`` / ``compute_edge_flux``
  con los mismos campos escritos en Python.
- **Resultantes analíticas** con la malla gmsh del ejemplo de Kirsch (cuarto
  de placa de 0.2 m con agujero de radio a = 0.01 m, espesor 0.01 m): una
  tracción uniforme ``σ`` en ``Borde_Carga`` (x = W) da ``ΣFx = σ·H·e``; una
  presión ``p`` en ``Agujero`` (``p > 0`` comprime) empuja la placa hacia
  fuera del agujero, ``ΣF = p·a·e·(1, 1)`` (``∫ n ds`` sobre el arco es la
  cuerda girada, también en el polígono discreto).
- **Extremo a extremo**: placa a tracción uniforme por ``tractions`` en una
  malla distorsionada, resuelta con ``solidum.run_yaml``: el campo de
  desplazamientos lineal exacto.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pytest

import solidum
from solidum.core.node import Node
from solidum.elements.solid_3d.hex8 import Hex8
from solidum.materials.elastic_3d import Elastic3D
from solidum.utils.expressions import ExpressionError, compile_expression, load_field
from solidum.utils.yaml_parser import YamlParser, YamlValidationError

KIRSCH = (pathlib.Path(__file__).resolve().parents[1] / "examples" / "placa_agujero_kirsch"
          / "placa_agujero.msh")
A, W, H, ESPESOR = 0.01, 0.2, 0.2, 0.01

# Dos Quad4 distorsionados (el lado común 2-5 inclinado), 1-2-5-4 y 2-3-6-5,
# con el contorno rectangular: libre de tracción bajo un esfuerzo uniforme σxx.
NODOS = ("nodes:\n"
         "  - {id: 1, coords: [0.0, 0.0]}\n  - {id: 2, coords: [1.1, 0.0]}\n"
         "  - {id: 3, coords: [2.0, 0.0]}\n  - {id: 4, coords: [0.0, 1.0]}\n"
         "  - {id: 5, coords: [0.8, 1.0]}\n  - {id: 6, coords: [2.0, 1.0]}\n")
MECANICO = (NODOS +
            "materials:\n  - {id: 1, type: Elastic2D, E: 200.0, nu: 0.25, hypothesis: plane_stress}\n"
            "elements:\n"
            "  - {id: 1, type: Quad4, nodes: [1, 2, 5, 4], material: 1, thickness: 0.5}\n"
            "  - {id: 2, type: Quad4, nodes: [2, 3, 6, 5], material: 1, thickness: 0.5}\n")
TERMICO = (NODOS +
           "thermal_materials:\n  - {id: 1, type: ThermalConduction, k: 2.0}\n"
           "elements:\n"
           "  - {id: 1, type: Quad4Thermal, nodes: [1, 2, 5, 4], material: 1}\n"
           "  - {id: 2, type: Quad4Thermal, nodes: [2, 3, 6, 5], material: 1}\n")


@pytest.fixture
def escribir(tmp_path):
    def _escribir(texto: str, nombre: str = "modelo.yaml") -> str:
        ruta = tmp_path / nombre
        ruta.write_text(texto, encoding="utf-8")
        return str(ruta)
    return _escribir


def leer(ruta):
    parser = YamlParser(ruta)
    parser.parse()
    parser.domain.generate_equation_numbers()
    return parser


def cerca(actual, esperado, rtol=1e-12):
    esperado = np.asarray(esperado, dtype=float)
    np.testing.assert_allclose(actual, esperado, rtol=rtol,
                               atol=rtol * max(np.abs(esperado).max(), 1e-300))


# ----------------------------------------------------------------------
# Expresiones
# ----------------------------------------------------------------------

class TestExpresiones:

    def test_aritmetica_y_funciones(self):
        f = compile_expression("1e3*(1 - y/2) + sin(pi*x)**2 - abs(-z) + max(x, 2)")
        assert f(0.5, 1.0, -0.25) == pytest.approx(500.0 + 1.0 - 0.25 + 2.0, rel=1e-15)

    @pytest.mark.parametrize("texto", [
        '__import__("os")', "x.real", "w + 1", "(lambda: 1)()", "x[0]", 'open("f")',
        '"a"', "x if y else 2", "sin(x=1)", "().__class__", "x == 1", "[x]", "True",
    ])
    def test_rechaza_lo_que_no_esta_en_la_lista(self, texto):
        with pytest.raises(ExpressionError, match="no se admite"):
            compile_expression(texto)

    def test_sintaxis_invalida(self):
        with pytest.raises(ExpressionError, match="sintaxis"):
            compile_expression("1 +* x")

    def test_potencia_enorme_desborda_en_vez_de_colgar(self):
        with pytest.raises(OverflowError):
            compile_expression("9**9**9**9")(0.0, 0.0, 0.0)

    def test_la_normal_solo_en_las_cargas_de_contorno(self):
        assert load_field("nx", 1, normal=True, context="c")([0, 0], [0.6, 0.8]) == 0.6
        with pytest.raises(ExpressionError, match="'nx'"):
            load_field("nx", 1, normal=False, context="body_force")

    def test_error_de_evaluacion_con_contexto(self):
        Q = load_field("log(x)", 1, normal=False, context="thermal_loads.body_source[0]")
        with pytest.raises(ExpressionError, match=r"body_source\[0\].*log\(x\)"):
            Q([-1.0, 0.0])


# ----------------------------------------------------------------------
# Igual que el API de Python
# ----------------------------------------------------------------------

class TestComoElAPI:

    def test_traccion_y_presion_por_elemento(self, escribir):
        p = leer(escribir(MECANICO + (
            "tractions:\n"
            "  - {element: 1, edge: 3, t: [\"1e3*y\", \"-2e2*x*ny\"]}\n"
            "  - {element: 2, edge: 1, pressure: \"5e2*(1 + y)\"}\n"
            "  - {element: 2, edge: 0, t: [10.0, -3.0]}\n")))
        el1, el2 = p.domain.elements[1], p.domain.elements[2]
        F = np.zeros(p.domain.total_dofs)
        F[el1.get_global_dof_indices()] += el1.compute_edge_traction(
            3, lambda x, n: np.array([1e3 * x[1], -2e2 * x[0] * n[1]]))
        F[el2.get_global_dof_indices()] += el2.compute_edge_traction(
            1, lambda x, n: -5e2 * (1 + x[1]) * np.asarray(n))
        F[el2.get_global_dof_indices()] += el2.compute_edge_traction(0, np.array([10.0, -3.0]))
        cerca(p.get_external_forces(), F)

    def test_fuerza_de_cuerpo_como_campo(self, escribir):
        p = leer(escribir(MECANICO + "body_force: [0.0, \"-78.5e3*(1 + 0.1*y)\"]\n"))
        assembler = solidum.math.assembly.Assembler(p.domain)
        b = lambda x: np.array([0.0, -78.5e3 * (1 + 0.1 * x[1])])
        cerca(p.get_body_load(assembler), assembler.assemble_body_load(b))

    def test_cargas_termicas(self, escribir):
        p = leer(escribir(TERMICO + (
            "thermal_loads:\n"
            "  body_source: {Q: \"5e3*exp(-x)\", elements: [2]}\n"
            "  boundary_flux:\n"
            "    - {element: 1, edge: 3, q: \"-200*(1 + y)*nx\"}\n"
            "    - {element: 2, edge: 1, q: 50.0}\n")))
        el1, el2 = p.domain.elements[1], p.domain.elements[2]
        F = np.zeros(p.domain.total_dofs)
        F[el2.get_global_dof_indices()] += el2.compute_body_source(lambda x: 5e3 * np.exp(-x[0]))
        F[el1.get_global_dof_indices()] += el1.compute_edge_flux(
            3, lambda x, n: -200 * (1 + x[1]) * n[0])
        F[el2.get_global_dof_indices()] += el2.compute_edge_flux(1, 50.0)
        cerca(p.get_thermal_loads(), F)

    def test_presion_en_una_cara_3d(self, escribir):
        """Hexaedro alabeado: la resultante de una presión uniforme es
        ``−p·½(d₁ × d₂)`` en la cara (vector área)."""
        X = np.array([[0, 0, 0], [1, 0, 0.1], [1.1, 1, 0], [0, 1, -0.1],
                      [0, 0.1, 1], [1, 0, 1.2], [1, 1, 1], [-0.1, 1, 0.9]])
        nodos = "".join(f"  - {{id: {i + 1}, coords: {list(map(float, x))}}}\n"
                        for i, x in enumerate(X))
        p = leer(escribir(
            "nodes:\n" + nodos +
            "materials:\n  - {id: 1, type: Elastic3D, E: 200.0, nu: 0.25}\n"
            "elements:\n  - {id: 1, type: Hex8, nodes: [1, 2, 3, 4, 5, 6, 7, 8], material: 1}\n"
            "tractions:\n  - {element: 1, face: 1, pressure: 3.0}\n"))
        el = p.domain.elements[1]
        P = X[list(Hex8.FACE_NODES[1])]
        area = 0.5 * np.cross(P[2] - P[0], P[3] - P[1])
        cerca(p.get_external_forces().reshape(-1, 3).sum(axis=0), -3.0 * area)
        cerca(p.get_external_forces()[el.get_global_dof_indices()],
              el.compute_face_traction(1, lambda x, n: -3.0 * np.asarray(n)))


# ----------------------------------------------------------------------
# Por grupo físico (malla gmsh de Kirsch)
# ----------------------------------------------------------------------

def kirsch(cargas: str) -> str:
    return (f"mesh: {KIRSCH.as_posix()}\nmesh_material: 1\nmesh_thickness: {ESPESOR}\n"
            "materials:\n  - {id: 1, type: Elastic2D, E: 200.0e9, nu: 0.3, "
            "hypothesis: plane_stress}\n" + cargas)


@pytest.mark.skipif(not KIRSCH.is_file(), reason="malla de Kirsch no disponible")
class TestPorGrupo:

    def resultante(self, p):
        F = p.get_external_forces()
        dom = p.domain
        ux = [n.dofs["ux"] for n in dom.nodes.values()]
        uy = [n.dofs["uy"] for n in dom.nodes.values()]
        return np.array([F[ux].sum(), F[uy].sum()])

    def test_traccion_uniforme_en_el_borde(self, escribir):
        p = leer(escribir(kirsch("tractions_by_group:\n  Borde_Carga: {t: [1.0e6, 0.0]}\n")))
        cerca(self.resultante(p), [1.0e6 * H * ESPESOR, 0.0], rtol=1e-12)

    def test_presion_en_el_agujero_empuja_hacia_fuera(self, escribir):
        p = leer(escribir(kirsch("tractions_by_group:\n  Agujero: {pressure: 2.0e6}\n")))
        cerca(self.resultante(p), [2.0e6 * A * ESPESOR] * 2, rtol=1e-12)

    def test_forma_de_lista(self, escribir):
        p = leer(escribir(kirsch(
            "tractions_by_group:\n  - {group_name: Borde_Carga, t: [1.0e6, 0.0]}\n")))
        cerca(self.resultante(p), [1.0e6 * H * ESPESOR, 0.0], rtol=1e-12)

    def test_flujo_termico_por_grupo(self, escribir):
        p = leer(escribir(
            f"mesh: {KIRSCH.as_posix()}\nmesh_element: Quad4Thermal\nmesh_material: 1\n"
            "thermal_materials:\n  - {id: 1, type: ThermalConduction, k: 45.0}\n"
            "thermal_loads:\n  boundary_flux_by_group:\n    Borde_Carga: {q: -300.0}\n"))
        assert p.get_thermal_loads().sum() == pytest.approx(300.0 * H, rel=1e-12)

    def test_grupo_de_superficie_rechazado(self, escribir):
        with pytest.raises(YamlValidationError, match="no son un borde"):
            leer(escribir(kirsch("tractions_by_group:\n  Placa: {t: [1.0, 0.0]}\n")))

    def test_grupo_inexistente(self, escribir):
        with pytest.raises(YamlValidationError, match="inexistente.*Borde_Carga"):
            leer(escribir(kirsch("tractions_by_group:\n  Borde: {t: [1.0, 0.0]}\n")))


# ----------------------------------------------------------------------
# Extremo a extremo
# ----------------------------------------------------------------------

def test_placa_a_traccion_uniforme(escribir):
    """Malla distorsionada de dos Quad4, ``σ`` en ``x = 2`` (borde 1 del
    elemento 2), ``ux = 0`` en ``x = 0`` y ``uy = 0`` en el nodo 1: el campo
    es el lineal exacto ``ux = σx/E``, ``uy = −νσy/E`` (esfuerzo plano)."""
    sigma, E, nu = 3.0, 200.0, 0.25
    ruta = escribir(MECANICO + (
        "boundary_conditions:\n"
        "  - {node_id: 1, ux: 0.0, uy: 0.0}\n  - {node_id: 4, ux: 0.0}\n"
        f"tractions:\n  - {{element: 2, edge: 1, t: [{sigma}, 0.0]}}\n"
        "solver:\n  type: LinearSolver\n"))
    res = solidum.run_yaml(ruta)
    p = leer(ruta)
    for node in p.domain.nodes.values():
        x, y = node.coordinates[:2]
        assert res.U[node.dofs["ux"]] == pytest.approx(sigma * x / E, abs=1e-12)
        assert res.U[node.dofs["uy"]] == pytest.approx(-nu * sigma * y / E, abs=1e-12)


# ----------------------------------------------------------------------
# Errores al leer el archivo
# ----------------------------------------------------------------------

@pytest.mark.parametrize("bloque, mensaje", [
    ("tractions:\n  - {element: 9, edge: 0, t: [1, 0]}\n", "inexistente"),
    ("tractions:\n  - {element: 1, edge: 4, t: [1, 0]}\n", "fuera de rango"),
    ("tractions:\n  - {element: 1, edge: 0}\n", "exactamente uno"),
    ("tractions:\n  - {element: 1, edge: 0, t: [1, 0], pressure: 2}\n", "exactamente uno"),
    ("tractions:\n  - {element: 1, t: [1, 0]}\n", "'edge'.*'face'"),
    ("tractions:\n  - {element: 1, edge: 0, t: [1, 0, 0]}\n", "2 componentes"),
    ("tractions:\n  - {element: 1, edge: 0, presion: 2}\n", "claves no admitidas"),
    ("tractions:\n  - {element: 1, edge: 0, pressure: \"2*w\"}\n", r"tractions\[0\].*'w'"),
    ("tractions_by_group:\n  Borde: {t: [1, 0]}\n", "grupos físicos"),
    ("body_force: [0.0, \"-9.81*rho\"]\n", "'rho'"),
])
def test_errores(escribir, bloque, mensaje):
    with pytest.raises(YamlValidationError, match=mensaje):
        leer(escribir(MECANICO + bloque))


def test_elemento_sin_tracciones(escribir):
    ruta = escribir(
        "nodes:\n  - {id: 1, coords: [0.0, 0.0]}\n  - {id: 2, coords: [1.0, 0.0]}\n"
        "materials:\n  - {id: 1, type: Elastic1D, E: 200.0}\n"
        "elements:\n  - {id: 1, type: Frame2DEuler, nodes: [1, 2], material: 1, A: 1.0, I: 1.0}\n"
        "tractions:\n  - {element: 1, edge: 0, t: [1, 0]}\n")
    with pytest.raises(YamlValidationError, match="no admite tracciones"):
        leer(ruta)


def test_sub_bloque_termico_desconocido(escribir):
    with pytest.raises(YamlValidationError, match="sub-bloques desconocidos"):
        leer(escribir(TERMICO + "thermal_loads:\n  boundary_fluxes: []\n"))
