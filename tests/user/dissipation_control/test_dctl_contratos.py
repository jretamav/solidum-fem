"""Contratos del módulo de usuario ``dissipation_control`` (ADR 0020).

El programa principal prueba sus propios solvers; cada módulo de usuario trae
su barrido. Aquí, sobre los solvers que registra el módulo:

- el contrato transversal que los hace despachables desde YAML
  (``PIPELINE_KIND``, firma del constructor y de ``solve``), con la tabla de
  familias de pipeline de ``test_solver_contract_sweep``;
- la declaración del primer paso heredada de ``ArcLengthSolver``
  (``initial_dlambda`` e ``initial_dl``, excluyentes) y la validación de
  ``initial_tau``;
- el gancho de inicio de paso ``Element.prepare_step`` (ADR 0020, P5): una
  llamada por intento de paso con el campo convergido, repetida con el mismo
  campo tras una bisección;
- la carga desde YAML: con ``user_modules: [dissipation_control]`` el solver
  se construye y resuelve; sin esa clave, el error lo dice y lista los
  módulos disponibles.

Las pruebas del gancho y del primer paso con este solver vivían en tests del
programa principal (``test_prepare_step_hook.py``,
``test_arclength_paso_inicial.py``) hasta que salió al módulo el 2026-09-29.
El elemento espía del gancho se repite aquí, en pequeño, para no depender del
estado de aquel archivo.
"""
from __future__ import annotations

import inspect
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

import numpy as np

_TESTS = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.abspath(os.path.join(_TESTS, '..')))
sys.path.insert(0, _TESTS)

import test_solver_contract_sweep as barrido  # noqa: E402

import solidum  # noqa: E402
from solidum.core.domain import Domain  # noqa: E402
from solidum.elements.truss import Truss2D  # noqa: E402
from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.math.assembly import Assembler  # noqa: E402
from solidum.registry import SolverRegistry  # noqa: E402

solidum.load_user_module("dissipation_control")
from solidum.user.dissipation_control import DissipationArcLengthSolver  # noqa: E402

ROOT = Path(_TESTS).parent
PAQUETE = "solidum.user.dissipation_control"


def solvers_del_modulo() -> list:
    """Solvers registrados cuya clase pertenece al módulo."""
    return sorted(n for n in SolverRegistry.names()
                  if (SolverRegistry.get(n).__module__ + ".").startswith(PAQUETE + "."))


# =========================================================================
# Contrato transversal de solver (el del barrido del programa principal)
# =========================================================================

class TestDctlRegistro(unittest.TestCase):

    def test_el_modulo_registra_su_solver(self):
        """Sin esto, los barridos de abajo pasarían por vacuidad."""
        self.assertEqual(solvers_del_modulo(), ["DissipationArcLengthSolver"])
        self.assertIs(SolverRegistry.get("DissipationArcLengthSolver"),
                      DissipationArcLengthSolver)


class TestDctlContratoDeclarativo(unittest.TestCase):

    def test_pipeline_kind_despachable(self):
        """``PARAMETROS_ESPERADOS`` cubre exactamente las familias canónicas
        de ``entry.py`` (lo comprueba el barrido del programa principal)."""
        for nombre in solvers_del_modulo():
            with self.subTest(solver=nombre):
                kind = getattr(SolverRegistry.get(nombre), 'PIPELINE_KIND', None)
                self.assertIn(kind, barrido.PARAMETROS_ESPERADOS)

    def test_primer_parametro_es_el_assembler(self):
        for nombre in solvers_del_modulo():
            with self.subTest(solver=nombre):
                params = [p for p in barrido._firma_init(nombre).parameters if p != 'self']
                self.assertEqual(params[0], 'assembler')

    def test_parametros_coherentes_con_la_familia(self):
        for nombre in solvers_del_modulo():
            with self.subTest(solver=nombre):
                clase = SolverRegistry.get(nombre)
                esperados = barrido.PARAMETROS_ESPERADOS[clase.PIPELINE_KIND]
                faltantes = esperados - barrido._parametros_requeridos(nombre)
                self.assertFalse(faltantes, f"{nombre} no pide {sorted(faltantes)}")
                if clase.PIPELINE_KIND == 'static':
                    intrusos = barrido._parametros_requeridos(nombre) & {'dt', 't_end', 'n_steps'}
                    self.assertFalse(intrusos, f"{nombre} es estático y exige {sorted(intrusos)}")


class TestDctlContratoDeSolve(unittest.TestCase):

    def test_los_estaticos_reciben_las_cargas_y_el_callback_es_opcional(self):
        for nombre in solvers_del_modulo():
            clase = SolverRegistry.get(nombre)
            with self.subTest(solver=nombre):
                firma = inspect.signature(clase.solve)
                params = [p for p in firma.parameters if p != 'self']
                if clase.PIPELINE_KIND == 'static':
                    self.assertIs(firma.parameters[params[0]].default,
                                  inspect.Parameter.empty,
                                  f"{nombre}.solve: las cargas deben ser obligatorias")
                if 'step_callback' in firma.parameters:
                    self.assertIsNot(firma.parameters['step_callback'].default,
                                     inspect.Parameter.empty)


# =========================================================================
# Primer paso y parámetros propios
# =========================================================================

def _barra(cls=Truss2D, E=200.0e9, A=1.0e-3, L=1.0, P=1.0e3):
    """Barra elástica empotrada-libre cargada en ``ux`` del extremo."""
    dom = Domain()
    n1 = dom.add_node(1, [0.0, 0.0])
    n2 = dom.add_node(2, [L, 0.0])
    barra = cls(1, [n1, n2], Elastic1D(E=E, density=7850.0), A=A)
    dom.add_element(barra)
    n1.fix_dof("ux", 0.0)
    n1.fix_dof("uy", 0.0)
    n2.fix_dof("uy", 0.0)
    dom.generate_equation_numbers(verbose=False)
    F = np.zeros(dom.total_dofs)
    F[n2.dofs["ux"]] = P
    return dom, barra, n2, F


class TestDctlParametros(unittest.TestCase):

    def test_declaracion_invalida_del_primer_paso(self):
        """``initial_dl`` e ``initial_dlambda`` son excluyentes, como en el padre."""
        dom, _, _, _ = _barra()
        with self.assertRaises(ValueError):
            DissipationArcLengthSolver(Assembler(dom), initial_tau=0.01, initial_dl=0.1,
                                       initial_dlambda=0.1)

    def test_initial_tau_obligatorio_y_positivo(self):
        dom, _, _, _ = _barra()
        with self.assertRaises(TypeError):
            DissipationArcLengthSolver(Assembler(dom))
        for tau in (0.0, -1.0e-3):
            with self.subTest(initial_tau=tau), self.assertRaises(ValueError):
                DissipationArcLengthSolver(Assembler(dom), initial_tau=tau)


# =========================================================================
# Gancho de inicio de paso (ADR 0020, P5)
# =========================================================================

class _BarraEspia(Truss2D):
    """Anota el ``U_committed`` que recibe en cada inicio de paso."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.llamadas = []

    def prepare_step(self, U_committed):
        self.llamadas.append(np.array(U_committed, copy=True))


class _BarraQueFallaUnIntento(_BarraEspia):
    """Hace fracasar el primer intento del segundo paso con una fuerza
    interna espuria ``±P`` que cambia de signo en cada ensamblaje: el Newton
    agota sus iteraciones, el solver biseca y repite el paso."""

    P = 1.0e3

    def prepare_step(self, U_committed):
        super().prepare_step(U_committed)
        self.envenenado = len(self.llamadas) == 2
        self.signo = 1.0

    def compute_element_state(self, u_e):
        K_e, F_int_e = super().compute_element_state(u_e)
        if getattr(self, "envenenado", False):
            self.signo = -self.signo
            F_int_e = F_int_e + self.signo * self.P * np.array([0.0, 0.0, 1.0, 0.0])
        return K_e, F_int_e


def _solver(asm):
    return DissipationArcLengthSolver(asm, initial_dlambda=0.4, initial_tau=1.0e-3)


class TestDctlGanchoDeInicioDePaso(unittest.TestCase):

    def test_llama_al_gancho_con_el_campo_convergido(self):
        dom, barra, n2, F = _barra(_BarraEspia)
        _solver(Assembler(dom)).solve(F)
        self.assertGreaterEqual(len(barra.llamadas), 2)
        np.testing.assert_array_equal(barra.llamadas[0], 0.0)
        self.assertGreater(barra.llamadas[1][n2.dofs["ux"]], 0.0)

    def test_tras_una_biseccion_el_gancho_se_repite_con_el_mismo_campo(self):
        E, A, L, P = 200.0e9, 1.0e-3, 1.0, _BarraQueFallaUnIntento.P
        dom, barra, n2, F = _barra(_BarraQueFallaUnIntento, E=E, A=A, L=L, P=P)
        U = _solver(Assembler(dom)).solve(F)
        self.assertGreaterEqual(len(barra.llamadas), 3)
        self.assertGreater(barra.llamadas[1][n2.dofs["ux"]], 0.0)
        np.testing.assert_array_equal(barra.llamadas[2], barra.llamadas[1])
        np.testing.assert_allclose(U[n2.dofs["ux"]], P * L / (E * A), rtol=1.0e-6)


# =========================================================================
# Carga desde YAML
# =========================================================================

_BARRA_YAML = textwrap.dedent("""\
    nodes:
      - {id: 1, coords: [0.0, 0.0]}
      - {id: 2, coords: [1.0, 0.0]}
    materials:
      - {id: 1, type: Elastic1D, E: 1000.0}
    elements:
      - {id: 1, type: Truss2D, material: 1, A: 1.0, nodes: [1, 2]}
    boundary_conditions_by_node:
      - {node_id: 1, ux: 0.0, uy: 0.0}
      - {node_id: 2, uy: 0.0}
    point_loads_by_node:
      - {node_id: 2, ux: 10.0}
    solver: {type: DissipationArcLengthSolver, initial_dlambda: 0.2, initial_tau: 1.0e-3}
    """)


class TestDctlYaml(unittest.TestCase):

    def test_con_user_modules_se_construye_y_resuelve(self):
        with tempfile.TemporaryDirectory() as d:
            ruta = Path(d) / "modelo.yaml"
            ruta.write_text("user_modules: [dissipation_control]\n" + _BARRA_YAML,
                            encoding="utf-8")
            resultado = solidum.run_yaml(ruta)
        # Barra elástica: el modo sigue cilíndrico y llega a λ = 1 con u = PL/(EA).
        self.assertTrue(resultado.converged)
        np.testing.assert_allclose(resultado.U.max(), 10.0 * 1.0 / (1000.0 * 1.0), rtol=1.0e-8)

    def test_sin_user_modules_el_error_lo_dice(self):
        """En un proceso limpio (en esta sesión el módulo ya está cargado):
        ``import solidum`` no lo carga, y un YAML que pide su solver sin
        declararlo falla con una pista que nombra ``user_modules`` y los
        módulos disponibles."""
        script = textwrap.dedent('''\
            import sys
            from pathlib import Path

            import solidum
            from solidum.utils.yaml_parser import YamlParser, YamlValidationError

            ruta = Path(sys.argv[1])
            modelo = Path(sys.argv[2]).read_text(encoding="utf-8")
            assert "dissipation_control" in solidum.available_user_modules()
            assert not any(m.startswith("solidum.user.dissipation_control")
                           for m in sys.modules), "import solidum cargó el módulo"

            ruta.write_text(modelo, encoding="utf-8")
            try:
                YamlParser(str(ruta)).parse()
                raise AssertionError("sin user_modules el solver del módulo debía fallar")
            except YamlValidationError as exc:
                texto = str(exc)
                assert "solver: tipo desconocido 'DissipationArcLengthSolver'" in texto, texto
                assert "decláralo en 'user_modules'" in texto, texto
                assert "dissipation_control" in texto, texto

            ruta.write_text("user_modules: [dissipation_control]\\n" + modelo, encoding="utf-8")
            parser = YamlParser(str(ruta))
            dom = parser.parse()
            dom.generate_equation_numbers(verbose=False)
            from solidum.math.assembly import Assembler
            solver = parser.get_solver(Assembler(dom))
            assert type(solver).__module__.startswith("solidum.user.dissipation_control."), \\
                type(solver).__module__
            assert solver.initial_tau == 1.0e-3 and solver.initial_dlambda == 0.2
            print("OK")
            ''')
        with tempfile.TemporaryDirectory() as d:
            base = Path(d) / "modelo_base.yaml"
            base.write_text(_BARRA_YAML, encoding="utf-8")
            r = subprocess.run(
                [sys.executable, "-c", script, str(Path(d) / "modelo.yaml"), str(base)],
                cwd=ROOT, capture_output=True, text=True, encoding="utf-8", timeout=300,
            )
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertTrue(r.stdout.strip().endswith("OK"), r.stdout)


if __name__ == '__main__':
    unittest.main()
