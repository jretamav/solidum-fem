"""Restricciones multipunto (MPC, ADR 0004) declaradas en YAML.

Dos bloques:

1. **Misma solución que la API Python.** Un modelo con MPC escrito en YAML
   (``linear_constraints``) da el mismo resultado, bit a bit, que el mismo
   modelo construido con ``Domain.add_linear_constraint``: en el problema
   mecánico (``LinearSolver``) y en el térmico (estacionario con
   ``LinearSolver`` y transitorio con ``ThetaMethodSolver``), con ``g = 0`` y
   ``g ≠ 0``. El modelo es una placa partida en dos por una columna de nodos
   duplicada y atada con ``u_s = u_m + g`` (con ``g = 0`` equivale a fundir
   los nodos, el oráculo de ``tests/test_mpc_nonlinear_static.py``).

2. **Validación temprana** (auditoría de las MPC, 2026-09-30, TE-7 y EST-9).
   Antes sólo se comprobaban la forma y la existencia de los nodos, y sólo
   con bloque ``nodes``; un DOF inexistente, un esclavo con apoyo, un ciclo
   o un coeficiente no numérico fallaban después, uno a uno y con números de
   ecuación internos, y una clave con errata (``G: 5``) se ignoraba y el
   modelo corría con ``g = 0``. Ahora todos esos errores se acumulan en un
   ``YamlValidationError`` que nombra la entrada y el par (nodo, dof), también
   con malla gmsh; y el ``Assembler`` traduce los errores del
   ``ConstraintSet`` de la API Python a (nodo, dof).
"""
from __future__ import annotations

import logging
import re
import textwrap
from pathlib import Path

import numpy as np
import pytest

import solidum
from solidum.core.domain import Domain
from solidum.elements.solid_2d.quad4 import Quad4
from solidum.elements.thermal.quad4_thermal import Quad4Thermal
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.math.assembly import Assembler
from solidum.utils.yaml_parser import YamlParser, YamlValidationError

K_COND, C_HEAT, RHO = 45.0, 460.0, 7850.0
E, NU, ESP = 200.0e9, 0.3, 0.02
NX, NY, L, H, IC = 4, 2, 0.4, 0.2, 2
DUP = 1000                    # desplazamiento de ids de los nodos duplicados

KIRSCH = (Path(__file__).resolve().parents[1] / "examples" / "placa_agujero_kirsch"
          / "placa_agujero.msh")


@pytest.fixture(autouse=True)
def _silencio():
    logging.disable(logging.CRITICAL)
    yield
    logging.disable(logging.NOTSET)


# ----------------------------------------------------------------------
# Modelos: placa partida por la columna i = IC, atada con MPC
# ----------------------------------------------------------------------

def _placa(fisica: str, g: float = 0.0):
    """Placa [0, L]×[0, H] de Quad4 (``'mecanica'``) o Quad4Thermal
    (``'termica'``). La columna ``i = IC`` se duplica para los elementos de su
    derecha, y cada nodo duplicado se ata a su original con ``u_s = u_m + g``
    (``ux`` y ``uy``, o ``T``). Mecánica: empotrada en ``x = 0``, carga vertical
    en el borde derecho. Térmica: ``T = 100`` en ``x = 0``, calor entrante en
    el borde derecho."""
    dom = Domain()
    if fisica == 'mecanica':
        mat, dofs = Elastic2D(E=E, nu=NU, hypothesis='plane_stress'), ('ux', 'uy')
    else:
        mat, dofs = ThermalConduction(k=K_COND, c=C_HEAT, density=RHO, dim=2), ('T',)
    nid, dup = {}, {}
    for j in range(NY + 1):
        for i in range(NX + 1):
            nid[i, j] = j * (NX + 1) + i + 1
            dom.add_node(nid[i, j], [L * i / NX, H * j / NY])
    for j in range(NY + 1):
        dup[IC, j] = DUP + nid[IC, j]
        dom.add_node(dup[IC, j], [L * IC / NX, H * j / NY])
    eid = 0
    for j in range(NY):
        for i in range(NX):
            eid += 1
            esquinas = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
            ids = [dup[a, b] if (a == IC and i >= IC) else nid[a, b] for a, b in esquinas]
            nodos = [dom.nodes[k] for k in ids]
            elem = (Quad4(eid, nodos, mat, thickness=ESP) if fisica == 'mecanica'
                    else Quad4Thermal(eid, nodos, mat, thickness=ESP))
            dom.add_element(elem)
    for j in range(NY + 1):
        for d in dofs:
            dom.nodes[nid[0, j]].fix_dof(d, 0.0 if fisica == 'mecanica' else 100.0)
    for j in range(NY + 1):
        for d in dofs:
            dom.add_linear_constraint((dup[IC, j], d), [(nid[IC, j], d)], [1.0], g=g)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    for j in range(NY + 1):
        n = dom.nodes[nid[NX, j]]
        if fisica == 'mecanica':
            F[n.dofs['uy']] = -1.0e4
        else:
            F[n.dofs['T']] = 50.0
    return dom, F


def _yaml_de(dom: Domain, F: np.ndarray, fisica: str, solver: str) -> str:
    """El mismo modelo que ``dom`` escrito en YAML (mismos ids y mismo orden
    de inserción de los nodos ⇒ misma numeración de ecuaciones)."""
    lineas = ["nodes:"]
    for n, node in dom.nodes.items():
        lineas.append(f"  - {{id: {n}, coords: [{node.coordinates[0]!r}, {node.coordinates[1]!r}]}}")
    if fisica == 'mecanica':
        lineas += ["materials:",
                   f"  - {{id: 1, type: Elastic2D, E: {E!r}, nu: {NU!r}, hypothesis: plane_stress}}"]
        tipo = 'Quad4'
    else:
        lineas += ["thermal_materials:",
                   f"  - {{id: 1, type: ThermalConduction, k: {K_COND!r}, c: {C_HEAT!r}, "
                   f"density: {RHO!r}}}"]
        tipo = 'Quad4Thermal'
    lineas.append("elements:")
    for e, el in dom.elements.items():
        lineas.append(f"  - {{id: {e}, type: {tipo}, nodes: {[x.id for x in el.nodes]}, "
                      f"material: 1, thickness: {ESP!r}}}")
    lineas.append("boundary_conditions:")
    for n, node in dom.nodes.items():
        for d, v in node.boundary_conditions.items():
            lineas.append(f"  - {{node_id: {n}, {d}: {v!r}}}")
    lineas.append("linear_constraints:")
    for lc in dom.linear_constraints:
        maestros = ", ".join(f"{{node: {m[0]}, dof: {m[1]}}}" for m in lc['masters'])
        lineas += [f"  - slave: {{node: {lc['slave'][0]}, dof: {lc['slave'][1]}}}",
                   f"    masters: [{maestros}]",
                   f"    coefficients: {lc['coefficients']!r}",
                   f"    g: {lc['g']!r}"]
    lineas.append("point_loads:")
    for n, node in dom.nodes.items():
        for d, idx in node.dofs.items():
            if F[idx]:
                lineas.append(f"  - {{node_id: {n}, {d}: {float(F[idx])!r}}}")
    lineas.append(solver)
    return "\n".join(lineas) + "\n"


def _escribir(tmp_path, texto: str, nombre: str = "modelo.yaml") -> str:
    ruta = tmp_path / nombre
    ruta.write_text(texto, encoding="utf-8")
    return str(ruta)


LINEAL = "solver:\n  type: LinearSolver"
THETA = ("solver:\n  type: ThetaMethodSolver\n  dt: 5.0\n  n_steps: 8\n"
         "  T_initial: 20.0\n  theta: 0.5")


# ----------------------------------------------------------------------
# 1. YAML y API Python dan el mismo resultado (TE-9, parte YAML)
# ----------------------------------------------------------------------

class TestMismoResultadoQueLaApiPython:

    @pytest.mark.parametrize("fisica", ["mecanica", "termica"])
    @pytest.mark.parametrize("g", [0.0, 3.0e-4])
    def test_estacionario(self, tmp_path, fisica, g):
        g = g if fisica == 'mecanica' else g * 1.0e4        # 3 K en la térmica
        dom, F = _placa(fisica, g)
        U_py = solidum.run(dom, F_applied=F).U
        U_yaml = solidum.run_yaml(_escribir(tmp_path, _yaml_de(dom, F, fisica, LINEAL))).U
        np.testing.assert_array_equal(U_yaml, U_py)
        # La restricción se cumple en la solución: u_s − u_m = g.
        d = 'uy' if fisica == 'mecanica' else 'T'
        s, m = dom.nodes[DUP + NX + 1 + IC + 1], dom.nodes[NX + 1 + IC + 1]
        assert U_yaml[s.dofs[d]] - U_yaml[m.dofs[d]] == pytest.approx(g, abs=1e-12)

    @pytest.mark.parametrize("g", [0.0, 3.0])
    def test_transitorio_termico(self, tmp_path, g):
        dom, F = _placa('termica', g)
        py = solidum.run_thermal_transient(dom, dt=5.0, n_steps=8, T_initial=20.0,
                                           theta=0.5, F_func=lambda t: F)
        dom2, F2 = _placa('termica', g)
        yml = solidum.run_yaml(_escribir(tmp_path, _yaml_de(dom2, F2, 'termica', THETA)))
        np.testing.assert_array_equal(yml.T_history, py.T_history)

    def test_con_g_nulo_equivale_a_fundir_los_nodos(self, tmp_path):
        """Oráculo del nodo fundido por la vía YAML: la placa atada con
        ``g = 0`` reproduce la placa conforme (la columna sin duplicar)."""
        dom, F = _placa('termica')
        U_mpc = solidum.run_yaml(_escribir(tmp_path, _yaml_de(dom, F, 'termica', LINEAL))).U
        fundido = Domain()
        mat = ThermalConduction(k=K_COND, c=C_HEAT, density=RHO, dim=2)
        for n, node in dom.nodes.items():
            if n < DUP:
                fundido.add_node(n, list(node.coordinates))
        for e, el in dom.elements.items():
            ids = [x.id - DUP if x.id > DUP else x.id for x in el.nodes]
            fundido.add_element(Quad4Thermal(e, [fundido.nodes[k] for k in ids], mat,
                                             thickness=ESP))
        for n, node in dom.nodes.items():
            for d, v in node.boundary_conditions.items():
                fundido.nodes[n].fix_dof(d, v)
        fundido.generate_equation_numbers()
        Ff = np.zeros(fundido.total_dofs)
        for n, node in dom.nodes.items():
            if F[node.dofs['T']]:
                Ff[fundido.nodes[n].dofs['T']] = F[node.dofs['T']]
        U_f = solidum.run(fundido, F_applied=Ff).U
        for n, node in dom.nodes.items():
            src = n - DUP if n > DUP else n
            assert U_mpc[node.dofs['T']] == pytest.approx(
                U_f[fundido.nodes[src].dofs['T']], abs=1e-10)


# ----------------------------------------------------------------------
# 2. Validación temprana del bloque linear_constraints (TE-7)
# ----------------------------------------------------------------------

def _modelo_termico_con(bloque: str, bc_extra: str = "") -> str:
    """Placa térmica de la sección 1 con ``bloque`` como ``linear_constraints``."""
    dom, F = _placa('termica')
    texto = _yaml_de(dom, F, 'termica', LINEAL)
    cabeza, resto = texto.split("linear_constraints:\n", 1)
    resto = resto[resto.index("point_loads:"):]
    return cabeza + bc_extra + bloque + resto


S, M = DUP + IC + 1, IC + 1          # esclavo y maestro de la fila j = 0
M2 = NX + 1 + IC + 1                 # maestro de la fila j = 1


def _lc(slave, masters, coefs, extra=""):
    ms = ", ".join(f"{{node: {n}, dof: {d}}}" for n, d in masters)
    return (f"  - slave: {{node: {slave[0]}, dof: {slave[1]}}}\n"
            f"    masters: [{ms}]\n    coefficients: {coefs}\n{extra}")


def _errores(tmp_path, texto: str) -> list[str]:
    with pytest.raises(YamlValidationError) as info:
        YamlParser(_escribir(tmp_path, texto)).parse()
    return info.value.errors


CASOS = {
    "esclavo_dof_inexistente": (
        _lc((S, 'ux'), [(M, 'T')], [1.0]),
        rf"linear_constraints\[0\]: el esclavo \(nodo {S}, 'ux'\) no existe: el nodo {S} "
        r"tiene los grados de libertad \['T'\]"),
    "maestro_dof_inexistente": (
        _lc((S, 'T'), [(M, 'Tz')], [1.0]),
        rf"linear_constraints\[0\]: el maestro \(nodo {M}, 'Tz'\) no existe"),
    "maestro_nodo_inexistente": (
        _lc((S, 'T'), [(99999, 'T')], [1.0]),
        r"linear_constraints\[0\]\.masters\[0\]: 'node=99999' no existe"),
    "esclavo_con_apoyo": (
        _lc((1, 'T'), [(2, 'T')], [1.0]),
        r"linear_constraints\[0\]: el esclavo \(nodo 1, 'T'\) tiene además un apoyo "
        r"\(T = 100\.0\).*Quita el apoyo del esclavo"),
    "ciclo": (
        _lc((S, 'T'), [(M, 'T')], [1.0]) + _lc((M, 'T'), [(S, 'T')], [1.0]),
        rf"restricción cíclica \(nodo {S}, 'T'\) ← \(nodo {M}, 'T'\) ← \(nodo {S}, 'T'\) "
        r"\(entradas \[0, 1\]\)"),
    "ciclo_largo": (
        _lc((S, 'T'), [(M, 'T')], [1.0]) + _lc((M, 'T'), [(M2, 'T')], [1.0])
        + _lc((M2, 'T'), [(S, 'T')], [1.0]),
        r"restricción cíclica .*\(entradas \[0, 1, 2\]\)"),
    "autorreferencia": (
        _lc((S, 'T'), [(S, 'T')], [1.0]),
        rf"el esclavo \(nodo {S}, 'T'\) figura entre sus propios maestros"),
    "esclavo_redeclarado_distinto": (
        _lc((S, 'T'), [(M, 'T')], [1.0]) + _lc((S, 'T'), [(M, 'T')], [2.0]),
        rf"linear_constraints\[1\]: el esclavo \(nodo {S}, 'T'\) ya está declarado en "
        r"linear_constraints\[0\]"),
    "coeficiente_no_numerico": (
        _lc((S, 'T'), [(M, 'T')], "[uno]"),
        r"linear_constraints\[0\]\.coefficients\[0\]: debe ser un número finito, "
        r"se encontró 'uno'"),
    "coeficiente_booleano": (
        _lc((S, 'T'), [(M, 'T')], "[true]"),
        r"coefficients\[0\]: debe ser un número finito, se encontró True"),
    "g_no_numerico": (
        _lc((S, 'T'), [(M, 'T')], [1.0], "    g: cinco\n"),
        r"linear_constraints\[0\]: 'g' debe ser un número finito, se encontró 'cinco'"),
    "clave_desconocida": (
        _lc((S, 'T'), [(M, 'T')], [1.0], "    G: 5.0\n"),
        r"linear_constraints\[0\]: claves desconocidas \['G'\]"),
    "clave_desconocida_en_maestro": (
        f"  - slave: {{node: {S}, dof: T}}\n"
        f"    masters: [{{node: {M}, dof: T, coef: 2.0}}]\n    coefficients: [1.0]\n",
        r"linear_constraints\[0\]\.masters\[0\]: claves desconocidas \['coef'\]"),
    "sin_esclavo": (
        f"  - masters: [{{node: {M}, dof: T}}]\n    coefficients: [1.0]\n",
        r"linear_constraints\[0\]: falta el campo obligatorio 'slave'"),
    "dof_no_es_nombre": (
        _lc((S, 3), [(M, 'T')], [1.0]),
        r"linear_constraints\[0\]\.slave: 'dof' debe ser el nombre"),
}


class TestValidacion:

    @pytest.mark.parametrize("caso", sorted(CASOS))
    def test_error_claro_al_leer(self, tmp_path, caso):
        bloque, patron = CASOS[caso]
        errores = _errores(tmp_path, _modelo_termico_con("linear_constraints:\n" + bloque))
        assert any(re.search(patron, e) for e in errores), errores

    def test_bloque_que_no_es_lista(self, tmp_path):
        errores = _errores(tmp_path, _modelo_termico_con(
            f"linear_constraints:\n  slave: {{node: {S}, dof: T}}\n"))
        assert errores == [e for e in errores if "'linear_constraints' debe ser una lista" in e]
        assert len(errores) == 1

    def test_errores_acumulados(self, tmp_path):
        """Varios problemas del modelo construido se informan juntos."""
        bloque = (_lc((S, 'ux'), [(M, 'T')], [1.0]) + _lc((1, 'T'), [(2, 'T')], [1.0])
                  + _lc((M, 'T'), [(M2, 'T')], [1.0]) + _lc((M2, 'T'), [(M, 'T')], [1.0]))
        errores = _errores(tmp_path, _modelo_termico_con("linear_constraints:\n" + bloque))
        assert len(errores) == 3, errores
        assert "no existe" in errores[0]
        assert "tiene además un apoyo" in errores[1]
        assert "restricción cíclica" in errores[2]

    def test_redeclaracion_identica_se_admite(self, tmp_path):
        """El ``ConstraintSet`` admite la misma restricción dos veces
        (idempotente); el lector también."""
        bloque = _lc((S, 'T'), [(M, 'T')], [1.0]) + _lc((S, 'T'), [(M, 'T')], [1.0])
        dom = YamlParser(_escribir(tmp_path, _modelo_termico_con(
            "linear_constraints:\n" + bloque))).parse()
        assert len(dom.linear_constraints) == 2

    def test_maestro_con_apoyo_se_admite(self, tmp_path):
        """Un maestro prescrito es legítimo: el cierre transitivo transmite su
        valor al esclavo."""
        dom = YamlParser(_escribir(tmp_path, _modelo_termico_con(
            "linear_constraints:\n" + _lc((S, 'T'), [(1, 'T')], [1.0], "    g: 5\n")))).parse()
        assert dom.linear_constraints[0]['g'] == 5.0
        dom.generate_equation_numbers()
        r = solidum.run(dom, F_applied=np.zeros(dom.total_dofs))
        assert r.U[dom.nodes[S].dofs['T']] == pytest.approx(105.0)


@pytest.mark.skipif(not KIRSCH.is_file(), reason="malla de Kirsch no disponible")
class TestValidacionConMallaGmsh:
    """Sin bloque ``nodes`` la validación de forma no conoce los nodos; la
    del modelo construido sí."""

    _BASE = """\
        thermal_materials:
          - {id: 1, type: ThermalConduction, k: 50.0}
        mesh_element: Quad4Thermal
        mesh_thickness: 0.02
        """

    def _ids(self, tmp_path):
        ruta = tmp_path / "base.yaml"
        ruta.write_text(f"mesh: {KIRSCH.as_posix()}\n" + textwrap.dedent(self._BASE),
                        encoding="utf-8")
        dom = YamlParser(str(ruta)).parse()
        return sorted(n for n, node in dom.nodes.items() if node.dofs)

    def test_errores_con_malla_gmsh(self, tmp_path):
        a, b = self._ids(tmp_path)[:2]
        texto = (f"mesh: {KIRSCH.as_posix()}\n" + textwrap.dedent(self._BASE)
                 + f"boundary_conditions_by_node:\n  - {{node_id: {a}, T: 20.0}}\n"
                 + "linear_constraints:\n"
                 + _lc((10 ** 9, 'T'), [(b, 'T')], [1.0])
                 + _lc((b, 'ux'), [(a, 'T')], [1.0])
                 + _lc((a, 'T'), [(b, 'T')], [1.0]))
        errores = _errores(tmp_path, texto)
        assert len(errores) == 3, errores
        assert f"no hay ningún nodo {10 ** 9} en el modelo" in errores[0]
        assert f"el esclavo (nodo {b}, 'ux') no existe" in errores[1]
        assert f"el esclavo (nodo {a}, 'T') tiene además un apoyo" in errores[2]

    def test_restriccion_valida_con_malla_gmsh(self, tmp_path):
        a, b = self._ids(tmp_path)[:2]
        texto = (f"mesh: {KIRSCH.as_posix()}\n" + textwrap.dedent(self._BASE)
                 + "linear_constraints:\n" + _lc((b, 'T'), [(a, 'T')], [1.0]))
        dom = YamlParser(_escribir(tmp_path, texto)).parse()
        assert dom.linear_constraints[0]['slave'] == (b, 'T')


# ----------------------------------------------------------------------
# 3. Mensajes del Assembler con la API Python (EST-9)
# ----------------------------------------------------------------------

class TestMensajesDelAssembler:

    def test_apoyo_redundante_en_el_esclavo(self):
        """Apoyo en los dos nodos coincidentes de una interfaz atada: el
        mensaje nombra esclavo y maestro por (nodo, dof) y dice qué hacer."""
        dom, _ = _placa('mecanica')
        dom.nodes[IC + 1].fix_dof('uy', 0.0)
        dom.nodes[S].fix_dof('uy', 0.0)
        with pytest.raises(ValueError) as info:
            Assembler(dom).constraint_set
        msg = str(info.value)
        assert f"esclavo (nodo {S}, 'uy')" in msg
        assert f"maestros [(nodo {M}, 'uy')]" in msg
        assert "Quita el apoyo del esclavo" in msg
        assert f"DOF {dom.nodes[S].dofs['uy']} (nodo {S}, 'uy')" in msg

    def test_ciclo_con_la_api_python(self):
        dom, _ = _placa('termica')
        dom.linear_constraints.clear()
        dom.add_linear_constraint((S, 'T'), [(M, 'T')], [1.0])
        dom.add_linear_constraint((M, 'T'), [(S, 'T')], [1.0])
        with pytest.raises(ValueError, match=rf"\(nodo ({S}|{M}), 'T'\)"):
            Assembler(dom).constraint_set

    def test_ciclo_de_tres_con_la_api_python(self):
        """Con tres o más esclavos el cierre transitivo del ``ConstraintSet``
        no llega a la autorreferencia: rota hasta agotar las iteraciones y
        sólo dice «posible ciclo», sin grados de libertad. El ``Assembler``
        nombra el ciclo por (nodo, dof)."""
        dom, _ = _placa('termica')
        dom.linear_constraints.clear()
        dom.add_linear_constraint((S, 'T'), [(M, 'T')], [1.0])
        dom.add_linear_constraint((M, 'T'), [(M2, 'T')], [1.0])
        dom.add_linear_constraint((M2, 'T'), [(S, 'T')], [1.0])
        with pytest.raises(ValueError) as info:
            Assembler(dom).constraint_set
        msg = str(info.value)
        assert "cíclicas" in msg
        for n in (S, M, M2):
            assert f"(nodo {n}, 'T')" in msg

    def test_autorreferencia_de_un_esclavo_con_apoyo(self):
        """El mensaje describe el error que ocurre (el esclavo entre sus
        maestros), no el apoyo del esclavo."""
        dom, _ = _placa('termica')
        dom.linear_constraints.clear()
        dom.add_linear_constraint((1, 'T'), [(1, 'T')], [1.0])
        with pytest.raises(ValueError) as info:
            Assembler(dom).constraint_set
        msg = str(info.value)
        assert "el esclavo figura entre sus propios maestros" in msg
        assert "Quita el apoyo" not in msg
