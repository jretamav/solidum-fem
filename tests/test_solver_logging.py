"""Registro de los solvers estáticos incrementales (ADR 0025).

Política de niveles:

- **(a)** Con el nivel por omisión (``INFO``), un empuje de 30 pasos escribe a
  lo sumo una línea por paso y ninguna línea de iteración; con ``DEBUG`` se
  recupera el detalle por iteración.
- **(b)** Un evento que el solver recupera por sí solo (tangente singular
  resuelta bisecando) se registra como ``WARNING`` con su motivo, nunca como
  ``ERROR``; ``ERROR`` queda para lo que termina el análisis.

Consola en cp1252: correr cualquier solver con
``PYTHONIOENCODING=cp1252`` no produce trazas ``--- Logging error ---``; los
caracteres que la consola no tiene se escriben en texto plano.

Las tangentes singulares se inyectan de forma determinista (``mock``): en los
modelos reales su aparición depende de la plataforma (el redondeo decide si un
mecanismo plástico deja la tangente exactamente singular).
"""
import io
import logging
import os
import subprocess
import sys
import textwrap
import unittest
from contextlib import contextmanager
from unittest import mock

import numpy as np

from solidum import Assembler, Domain
from solidum.elements.frame import Frame2DEuler
from solidum.logging import _ConsoleHandler, encodable_text
from solidum.materials.elastic import Elastic1D
from solidum.materials.plastic_1d import Elastoplastic1D
from solidum.math.solvers import ArcLengthSolver, IndirectDisplacementSolver, NonlinearSolver
from solidum.math.solvers.arclength import _ArcProblem
from solidum.math.solvers.diagnostics import SolverDivergedError
from solidum.math.solvers.nonlinear import _IncrementalProblem
from solidum.sections import FiberSection, RectangularPatch

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
E, L, B, H = 2.0e11, 3.0, 0.30, 0.50


def plastic_cantilever(n_el=4):
    """Voladizo de 20 fibras de acero con endurecimiento del 1 %, carga
    lateral de referencia de 1 MN en la cabeza."""
    steel = Elastoplastic1D(E=E, sigma_y=250.0e6, H=2.0e9)
    section = FiberSection(patches=[RectangularPatch(steel, y=(-H / 2, H / 2),
                                                     z=(-B / 2, B / 2), ny=20, nz=1)])
    dom = Domain()
    nodes = [dom.add_node(i + 1, [0.0, i * L / n_el]) for i in range(n_el + 1)]
    for k in range(n_el):
        dom.add_element(Frame2DEuler(k + 1, [nodes[k], nodes[k + 1]], section=section, n_ip=3))
    for dof in ("ux", "uy", "rz"):
        nodes[0].fix_dof(dof, 0.0)
    dom.generate_equation_numbers()
    F = np.zeros(dom.total_dofs)
    F[nodes[-1].dofs["ux"]] = 1.0e6
    return dom, nodes[-1], F


@contextmanager
def captured(level):
    """Registros de ``solidum`` con el nivel dado, sin escribirlos en pantalla."""
    root = logging.getLogger("solidum")
    records = []

    class Collect(logging.Handler):
        def emit(self, record):
            records.append(record)

    saved_handlers, saved_level = list(root.handlers), root.level
    root.handlers = [Collect(level=logging.DEBUG)]
    root.setLevel(level)
    logging.disable(logging.NOTSET)
    try:
        yield records
    finally:
        root.handlers = saved_handlers
        root.setLevel(saved_level)


def _messages(records, levelno):
    return [r.getMessage() for r in records if r.levelno == levelno]


def _push(steps=30, **kw):
    dom, tip, F = plastic_cantilever()
    kw.setdefault("max_steps", steps)
    solver = IndirectDisplacementSolver(Assembler(dom), control=[(tip.id, "ux", 1.0)],
                                        max_lambda=1.0e9, initial_dl=1.5e-3, **kw)
    solver.solve(F)
    return solver


class TestVerbosidad(unittest.TestCase):
    """(a) Una línea por paso en INFO; el detalle, en DEBUG."""

    def test_empuje_de_30_pasos_en_info(self):
        with captured(logging.INFO) as records:
            solver = _push(30)
        self.assertEqual(solver.steps_done, 30)
        info = _messages(records, logging.INFO)
        self.assertFalse([m for m in info if "Iteración" in m])
        self.assertEqual(len([m for m in info if m.startswith("[PASO")]), 30)
        # Inicio, primer paso y una línea por paso; nada más.
        self.assertLessEqual(len(info), 30 + 3)

    def test_debug_recupera_el_detalle_por_iteracion(self):
        with captured(logging.DEBUG) as records:
            _push(10)
        debug = _messages(records, logging.DEBUG)
        self.assertGreaterEqual(len([m for m in debug if "Iteración" in m]), 10)

    def test_newton_incremental_en_info(self):
        dom, tip, F = plastic_cantilever()
        with captured(logging.INFO) as records:
            NonlinearSolver(Assembler(dom), num_steps=20).solve(1.2 * F)
        info = _messages(records, logging.INFO)
        self.assertFalse([m for m in info if "Iteración" in m])
        steps = [m for m in info if m.startswith("[PASO")]
        self.assertGreaterEqual(len(steps), 20)
        self.assertLessEqual(len(info), len(steps) + 2)     # inicio y resumen


def _fail_once(original, at_call):
    """Envoltura que lanza ``RuntimeError`` (tangente singular) sólo en la
    llamada número ``at_call`` y delega en ``original`` en las demás."""
    calls = {"n": 0}

    def wrapper(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == at_call:
            raise RuntimeError("Factor is exactly singular")
        return original(*args, **kwargs)
    return wrapper, calls


class TestEventoRecuperadoNoEsError(unittest.TestCase):
    """(b) WARNING con el motivo para lo recuperado; ERROR sólo al detenerse."""

    def _assert_recovered(self, records, solver, needle):
        self.assertTrue(solver.completed)
        self.assertEqual(_messages(records, logging.ERROR), [])
        warnings = _messages(records, logging.WARNING)
        self.assertTrue([w for w in warnings if needle in w and "bisección" in w], warnings)
        self.assertGreaterEqual(solver.n_bisections, 1)

    def test_tangente_singular_en_el_corrector_del_arco(self):
        wrapper, calls = _fail_once(_ArcProblem.correction, at_call=5)
        with mock.patch.object(_ArcProblem, "correction", autospec=True, side_effect=wrapper):
            with captured(logging.INFO) as records:
                solver = _push(200, control_target=3.0e-2)
        self.assertGreaterEqual(calls["n"], 5)
        self._assert_recovered(records, solver, "tangente singular")
        summary = [m for m in _messages(records, logging.INFO) if "trazado completo" in m]
        self.assertTrue(summary and "bisección" in summary[0])

    def test_tangente_singular_en_el_predictor_del_arco(self):
        wrapper, calls = _fail_once(ArcLengthSolver._tangent_predictor, at_call=4)

        def returns_none(*args, **kwargs):
            try:
                return wrapper(*args, **kwargs)
            except RuntimeError:
                return None          # así señala el predictor una tangente singular
        with mock.patch.object(ArcLengthSolver, "_tangent_predictor", autospec=True,
                               side_effect=returns_none):
            with captured(logging.INFO) as records:
                solver = _push(200, control_target=1.0e-2)
        self.assertGreaterEqual(calls["n"], 4)
        self._assert_recovered(records, solver, "tangente singular en el predictor")

    def test_tangente_singular_en_el_newton_incremental(self):
        dom, tip, F = plastic_cantilever()
        wrapper, calls = _fail_once(_IncrementalProblem.correction, at_call=3)
        with mock.patch.object(_IncrementalProblem, "correction", autospec=True, side_effect=wrapper):
            with captured(logging.INFO) as records:
                solver = NonlinearSolver(Assembler(dom), num_steps=8)
                solver.solve(1.2 * F)
        self._assert_recovered(records, solver, "tangente singular")

    def test_error_solo_al_detenerse(self):
        """El arco que no pasa un retroceso: un único ERROR, al final."""
        from tests.test_trace_termination import damaged_bar
        dom, F = damaged_bar()
        solver = ArcLengthSolver(Assembler(dom), max_iter=30, max_lambda=1.2,
                                 max_steps=400, initial_dlambda=0.02)
        with captured(logging.INFO) as records:
            with self.assertRaises(SolverDivergedError):
                solver.solve(F)
        errors = _messages(records, logging.ERROR)
        self.assertEqual(len(errors), 1)
        self.assertIn("análisis detenido", errors[0])
        warnings = _messages(records, logging.WARNING)
        self.assertTrue([w for w in warnings if "Raíces imaginarias" in w])


class TestConsolaSinUnicode(unittest.TestCase):
    """El registro se adapta a la codificación de la consola."""

    SAMPLE = ("Primer paso: Δλ₁ = 0.5 ⇒ Δl₁ = 3.0e-05 | α=0.5 | ω₁ = 12.6 rad/s | "
              "κ̃ ≤ 1e12 | ‖R‖ | σ ς | kg·m² | áéñ")

    def test_texto_plano_en_cp1252_y_ascii(self):
        cp = encodable_text(self.SAMPLE, "cp1252")
        self.assertIn("Deltalambda_1 = 0.5 => Deltal_1", cp)
        self.assertIn("alpha=0.5 | omega_1", cp)
        self.assertIn("kappa <= 1e12 | ||R|| | sigma sigma | kg·m² | áéñ", cp)   # cp1252 tiene ·, ², á
        cp.encode("cp1252")
        asc = encodable_text(self.SAMPLE, "ascii")
        self.assertIn("kg*m^2 | aen", asc)
        asc.encode("ascii")

    def test_sin_cambios_si_la_codificacion_lo_admite(self):
        self.assertEqual(encodable_text(self.SAMPLE, "utf-8"), self.SAMPLE)
        self.assertEqual(encodable_text(self.SAMPLE, None), self.SAMPLE)

    def test_el_manejador_escribe_en_un_flujo_cp1252(self):
        raw = io.BytesIO()
        stream = io.TextIOWrapper(raw, encoding="cp1252", errors="strict")
        handler = _ConsoleHandler(stream=stream)
        handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        record = logging.LogRecord("solidum.t", logging.INFO, __file__, 1, self.SAMPLE, None, None)
        with mock.patch.object(handler, "handleError", side_effect=AssertionError("Logging error")):
            handler.emit(record)
        stream.flush()
        self.assertIn("[INFO] Primer paso: Deltalambda_1", raw.getvalue().decode("cp1252"))

    def test_solvers_en_una_consola_cp1252_y_ascii(self):
        """Lineal, Newton, arco, control indirecto, modal y HHT en DEBUG (el
        nivel con más mensajes) en un subproceso con la consola en cp1252 y en
        ASCII: ninguna traza de registro."""
        code = textwrap.dedent(f"""
            import sys
            sys.path.insert(0, {ROOT!r})
            import numpy as np
            import solidum
            from solidum import Assembler
            from solidum.math.solvers import (ArcLengthSolver, HHTSolver, IndirectDisplacementSolver,
                                              LinearSolver, ModalSolver, NonlinearSolver)
            from tests.test_solver_logging import plastic_cantilever
            solidum.set_log_level("DEBUG")

            def model():
                dom, tip, F = plastic_cantilever()
                for el in dom.elements.values():
                    for m in el.materials():
                        m.density = 7850.0
                return dom, tip, F

            dom, tip, F = model(); LinearSolver(Assembler(dom)).solve(0.1 * F)
            dom, tip, F = model(); NonlinearSolver(Assembler(dom), num_steps=4).solve(1.1 * F)
            dom, tip, F = model(); ArcLengthSolver(Assembler(dom), max_lambda=1.1).solve(F)
            dom, tip, F = model()
            IndirectDisplacementSolver(Assembler(dom), control=[(tip.id, "ux", 1.0)], max_lambda=1e9,
                                       control_target=2.0e-2, initial_dl=4.0e-3).solve(F)
            dom, tip, F = model(); ModalSolver(Assembler(dom), n_modes=2).solve()
            dom, tip, F = model()
            HHTSolver(Assembler(dom), alpha=-0.05, t_end=0.01, dt=0.005,
                      F_func=lambda t: 0.1 * F).solve()
            print("FIN")
        """)
        for encoding in ("cp1252", "ascii"):
            with self.subTest(encoding=encoding):
                out = subprocess.run(
                    [sys.executable, "-c", code], capture_output=True, cwd=ROOT, timeout=600,
                    env={**os.environ, "PYTHONIOENCODING": encoding})
                stdout = out.stdout.decode(encoding, errors="replace")
                stderr = out.stderr.decode("utf-8", errors="replace")
                self.assertEqual(out.returncode, 0, stderr[-2000:])
                self.assertNotIn("Logging error", stderr)
                self.assertNotIn("UnicodeEncodeError", stderr)
                self.assertIn("Deltalambda_1", stdout)
                self.assertIn("FIN", stdout)


if __name__ == "__main__":
    unittest.main()
