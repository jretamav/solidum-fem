"""Red de seguridad del análisis estático (ADR 0019).

Antes de esta red, un modelo mal apoyado resuelto con ``LinearSolver``
devolvía **en silencio** desplazamientos de miles de kilómetros con SuperLU
y de decenas de kilómetros con Pardiso (que perturba los pivotes nulos y
sigue). Estas pruebas fijan el comportamiento nuevo:

1. **Capa 1 — mecanismos rígidos, antes de resolver**: se detectan con
   cualquier backend y en los tres solvers estáticos, y el mensaje describe
   el movimiento libre (traslación, giro por un punto, campo escalar sin
   valor prescrito). Las restricciones lineales (MPC) cuentan como apoyo.
   **No** se aplica a análisis dinámicos, donde un modelo libre es legítimo.
2. **Capas 2 y 3 — equilibrio y pivotes, después de resolver** (sólo
   ``LinearSolver``): un mecanismo interno se rechaza, cargado o no, con
   SuperLU y con Pardiso.
3. **Diagnóstico del no lineal**: un mecanismo interno se reporta como
   tangente singular, no como "modo no clasificado".
4. **Memoria**: un solver directo sin memoria sugiere el iterativo.
"""
import logging
import os
import sys
import unittest
import warnings
from unittest import mock

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum
from solidum import Assembler, Domain
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.math.linalg import LUSolver
from solidum.math.linalg.dispatcher import _HAS_PARDISO
from solidum.math.solvers.arclength import ArcLengthSolver
from solidum.math.solvers.diagnostics import (
    IllPosedSystemError,
    MechanismError,
    SingularTangentError,
)
from solidum.math.solvers.linear import LinearSolver
from solidum.math.solvers.model_checks import unrestrained_rigid_motions
from solidum.math.solvers.nonlinear import NonlinearSolver
from solidum.registry import ElementRegistry, MaterialRegistry

EL = ElementRegistry.get
Elastic1D = MaterialRegistry.get("Elastic1D")
Elastic2D = MaterialRegistry.get("Elastic2D")
Elastic3D = MaterialRegistry.get("Elastic3D")

DIRECT = ("lu",) + (("pardiso",) if _HAS_PARDISO else ())


def setUpModule():
    logging.disable(logging.WARNING)
    warnings.simplefilter("ignore", RuntimeWarning)


def tearDownModule():
    logging.disable(logging.NOTSET)
    warnings.resetwarnings()


# ----------------------------------------------------------------------
# Modelos
# ----------------------------------------------------------------------

def plate(supports, nx=6, ny=2, lx=4.0, ly=1.0):
    """Placa Quad4 ``[0, lx] × [0, ly]``. ``supports(domain, left_nodes)``
    impone los apoyos. Devuelve ``(domain, nodo de la esquina superior
    derecha)``."""
    d = Domain()
    nid = 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            d.add_node(nid, [lx * i / nx, ly * j / ny]); nid += 1
    m = Elastic2D(E=210e9, nu=0.3)
    for j in range(ny):
        for i in range(nx):
            n0 = 1 + i + j * (nx + 1)
            d.add_element(EL("Quad4")(len(d.elements) + 1,
                                      [d.nodes[n0], d.nodes[n0 + 1], d.nodes[n0 + nx + 2], d.nodes[n0 + nx + 1]], m))
    left = [d.nodes[1 + j * (nx + 1)] for j in range(ny + 1)]
    supports(d, left)
    d.generate_equation_numbers()
    return d, d.nodes[(nx + 1) * (ny + 1)]


def clamp(d, left):
    for n in left:
        n.fix_dof("ux", 0.0); n.fix_dof("uy", 0.0)


def rollers_uy(d, left):
    for n in left:
        n.fix_dof("uy", 0.0)


def single_pin(d, left):
    left[0].fix_dof("ux", 0.0); left[0].fix_dof("uy", 0.0)


def hinged_blocks():
    """Bloque A ``[0,2]×[0,1]`` empotrado en ``x = 0``; bloque B
    ``[2,4]×[1,2]`` unido a A por **un solo nodo**, ``(2, 1)``: B puede girar
    respecto a A. Mecanismo interno: los apoyos restringen todo movimiento
    rígido del conjunto, así que la capa 1 no lo ve."""
    d = Domain(); m = Elastic2D(E=210e9, nu=0.3); ids = {}

    def node(x, y):
        k = (round(x, 9), round(y, 9))
        if k not in ids:
            ids[k] = len(ids) + 1; d.add_node(ids[k], [x, y])
        return d.nodes[ids[k]]

    def block(x0, y0):
        for j in range(2):
            for i in range(4):
                xs = [x0 + i * 0.5, x0 + (i + 1) * 0.5]; ys = [y0 + j * 0.5, y0 + (j + 1) * 0.5]
                d.add_element(EL("Quad4")(len(d.elements) + 1,
                                          [node(xs[0], ys[0]), node(xs[1], ys[0]), node(xs[1], ys[1]), node(xs[0], ys[1])], m))
    block(0.0, 0.0); block(2.0, 1.0)
    for (x, _y), k in ids.items():
        if abs(x) < 1e-9:
            d.nodes[k].fix_dof("ux", 0.0); d.nodes[k].fix_dof("uy", 0.0)
    d.generate_equation_numbers()
    return d, ids


def run_static(d, loads, solver_factory):
    a = Assembler(d)
    a.assemble_system()
    for node_id, dof, value in loads:
        a.apply_point_load(node_id, dof, value)
    return solidum.run(d, assembler=a, solver=solver_factory(a), F_applied=a.F_global.copy())


# ----------------------------------------------------------------------
# 1. Capa 1: mecanismos rígidos
# ----------------------------------------------------------------------

class TestRigidMechanismDetection(unittest.TestCase):

    def _motions(self, d):
        a = Assembler(d)
        a.assemble_system()
        return unrestrained_rigid_motions(a)

    def test_well_supported_models_pass(self):
        d, _ = plate(clamp)
        self.assertEqual(self._motions(d), [])

    def test_rollers_on_a_vertical_line_leave_translation_and_rotation(self):
        """Impedir sólo ``uy`` en ``x = 0`` deja libres la traslación en x y
        el giro alrededor de cualquier punto de esa recta (mueve esos nodos
        en horizontal, que está permitido)."""
        d, _ = plate(rollers_uy)
        motions = self._motions(d)
        self.assertEqual(len(motions), 2)
        text = " | ".join(motions)
        self.assertIn("traslación en la dirección x", text)
        self.assertIn("giro alrededor del eje z que pasa por (0,", text)

    def test_single_pin_leaves_rotation_about_the_pin(self):
        d, _ = plate(single_pin)
        self.assertEqual(self._motions(d), ["giro alrededor del eje z que pasa por (0, 0)"])

    def test_free_3d_solid_has_six_motions(self):
        d = Domain(); nid = 1; ids = {}
        for k in range(2):
            for j in range(2):
                for i in range(2):
                    d.add_node(nid, [float(i), float(j), float(k)]); ids[(i, j, k)] = nid; nid += 1
        c = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
        d.add_element(EL("Hex8")(1, [d.nodes[ids[q]] for q in c], Elastic3D(E=1e9, nu=0.3)))
        d.generate_equation_numbers()
        motions = self._motions(d)
        self.assertEqual(len(motions), 6)
        self.assertEqual(sum("traslación" in m for m in motions), 3)
        self.assertEqual(sum("giro" in m for m in motions), 3)
        for q in c[:4]:                            # empotrar la cara z = 0
            for dof in ("ux", "uy", "uz"):
                d.nodes[ids[q]].fix_dof(dof, 0.0)
        self.assertEqual(self._motions(d), [])

    def test_frame_on_rollers_can_slide(self):
        d = Domain()
        for i, p in enumerate([(0, 0), (0, 3), (4, 3), (4, 0)], 1):
            d.add_node(i, list(p))
        m = Elastic1D(E=200e9)
        for e, (a_, b_) in enumerate([(1, 2), (2, 3), (3, 4)], 1):
            d.add_element(EL("Frame2DEuler")(e, [d.nodes[a_], d.nodes[b_]], m, A=1e-2, I=1e-4))
        d.nodes[1].fix_dof("uy", 0.0); d.nodes[4].fix_dof("uy", 0.0)
        d.generate_equation_numbers()
        self.assertEqual(self._motions(d), ["traslación en la dirección x"])
        d.nodes[1].fix_dof("ux", 0.0)              # un apoyo fijo lo resuelve
        self.assertEqual(self._motions(d), [])

    def test_linear_constraints_count_as_supports(self):
        """Una restricción lineal que ata ``ux`` a un DOF prescrito restringe;
        una que ata dos DOF libres entre sí, no."""
        # Atar ux de las dos esquinas izquierdas entre sí (dos DOF libres): el
        # giro alrededor de un punto de x = 0 las movería distinto en x, así
        # que queda impedido; la traslación en x, no (se mueven igual).
        d, _ = plate(rollers_uy)
        n_bl, n_tl = d.nodes[1], d.nodes[1 + 2 * 7]
        d.add_linear_constraint(slave=(n_bl.id, "ux"), masters=[(n_tl.id, "ux")], coefficients=[1.0])
        self.assertEqual(self._motions(d), ["traslación en la dirección x"])

        # Atar ux de cada esquina a su propio uy prescrito (ux = uy = 0): las
        # dos esquinas quedan fijas en x a través de la restricción lineal, y
        # con los rodillos en y ya no queda movimiento libre.
        d2, _ = plate(rollers_uy)
        for n in (d2.nodes[1], d2.nodes[1 + 2 * 7]):
            d2.add_linear_constraint(slave=(n.id, "ux"), masters=[(n.id, "uy")], coefficients=[1.0])
        self.assertEqual(self._motions(d2), [])

        # Con una sola esquina fijada así queda el giro alrededor de ella: el
        # resto del borde sólo tiene impedido uy, y ese giro lo mueve en x.
        d3, _ = plate(rollers_uy)
        d3.add_linear_constraint(slave=(d3.nodes[1].id, "ux"), masters=[(d3.nodes[1].id, "uy")],
                                 coefficients=[1.0])
        self.assertEqual(self._motions(d3), ["giro alrededor del eje z que pasa por (0, 0)"])

    def test_steady_thermal_without_dirichlet(self):
        d = Domain()
        for i, p in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
            d.add_node(i, list(p))
        d.add_element(EL("Quad4Thermal")(1, [d.nodes[k] for k in (1, 2, 3, 4)], ThermalConduction(k=50.0)))
        d.generate_equation_numbers()
        motions = self._motions(d)
        self.assertEqual(len(motions), 1)
        self.assertIn("el campo 'T' no tiene ningún valor prescrito", motions[0])
        d.nodes[1].fix_dof("T", 20.0)
        self.assertEqual(self._motions(d), [])

    def test_mixed_domain_reports_only_the_free_field(self):
        d = Domain()
        for i, p in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
            d.add_node(i, list(p))
        ns = [d.nodes[k] for k in (1, 2, 3, 4)]
        d.add_element(EL("Quad4")(1, ns, Elastic2D(E=210e9, nu=0.3)))
        d.add_element(EL("Quad4Thermal")(2, ns, ThermalConduction(k=50.0)))
        for n in (d.nodes[1], d.nodes[4]):
            n.fix_dof("ux", 0.0); n.fix_dof("uy", 0.0)
        d.generate_equation_numbers()
        motions = self._motions(d)
        self.assertEqual(len(motions), 1)
        self.assertIn("'T'", motions[0])

    def test_all_static_solvers_refuse_with_any_backend(self):
        loads = lambda tip: [(tip.id, "ux", 1e3), (tip.id, "uy", -1e3)]
        factories = {
            "LinearSolver": lambda la: (lambda a: LinearSolver(a, linear_algebra=la)),
            "NonlinearSolver": lambda la: (lambda a: NonlinearSolver(a, num_steps=1, linear_algebra=la)),
            "ArcLengthSolver": lambda la: (lambda a: ArcLengthSolver(a, max_lambda=1.0, linear_algebra=la)),
        }
        for name, make in factories.items():
            for la in DIRECT + ("iterative:none",):
                with self.subTest(solver=name, backend=la):
                    d, tip = plate(rollers_uy)
                    with self.assertRaises(MechanismError) as cm:
                        run_static(d, loads(tip), make(la))
                    self.assertIn(name, str(cm.exception))
                    self.assertEqual(len(cm.exception.motions), 2)

    def test_dynamics_of_a_free_body_is_legitimate(self):
        """En dinámica un cuerpo libre es un problema bien planteado (la masa
        lo regulariza): la capa 1 no debe intervenir."""
        from solidum.math.solvers.newmark import NewmarkSolver
        d = Domain(); nid = 1
        for j in range(3):
            for i in range(5):
                d.add_node(nid, [i * 0.5, j * 0.5]); nid += 1
        m = Elastic2D(E=210e9, nu=0.3, density=7850.0)
        for j in range(2):
            for i in range(4):
                n0 = 1 + i + j * 5
                d.add_element(EL("Quad4")(len(d.elements) + 1,
                                          [d.nodes[n0], d.nodes[n0 + 1], d.nodes[n0 + 6], d.nodes[n0 + 5]], m))
        d.generate_equation_numbers()
        a = Assembler(d)
        F = np.zeros(d.total_dofs); F[d.nodes[15].dofs["ux"]] = 1e3
        s = NewmarkSolver(a, t_end=3e-4, dt=1e-4, F_func=lambda t: F)
        r = solidum.entry.run_transient(d, assembler=a, solver=s)
        self.assertTrue(np.all(np.isfinite(r.u_history)))


# ----------------------------------------------------------------------
# 2. Capas 2 y 3: equilibrio y pivotes (LinearSolver)
# ----------------------------------------------------------------------

class TestLinearSolutionChecks(unittest.TestCase):

    def test_internal_mechanism_loaded(self):
        d, ids = hinged_blocks()
        for la in DIRECT:
            with self.subTest(backend=la):
                with self.assertRaises(IllPosedSystemError):
                    run_static(d, [(ids[(4.0, 2.0)], "uy", -1e3)],
                               lambda a: LinearSolver(a, linear_algebra=la))

    def test_internal_mechanism_not_activated_by_loads(self):
        """Carga sólo en el bloque empotrado: el sistema es singular pero
        consistente y el residuo sale perfecto; lo delatan los pivotes
        nulos, con el mismo criterio en SuperLU y en Pardiso."""
        d, ids = hinged_blocks()
        for la in DIRECT:
            with self.subTest(backend=la):
                with self.assertRaises(IllPosedSystemError) as cm:
                    run_static(d, [(ids[(1.0, 1.0)], "uy", -1e3)],
                               lambda a: LinearSolver(a, linear_algebra=la))
                self.assertGreaterEqual(cm.exception.zero_pivots, 1)

    def test_well_posed_model_passes_every_check(self):
        for la in DIRECT + ("iterative",):
            with self.subTest(backend=la):
                d, tip = plate(clamp)
                r = run_static(d, [(tip.id, "uy", -1e3)], lambda a: LinearSolver(a, linear_algebra=la))
                self.assertTrue(np.all(np.isfinite(r.U)))
                self.assertLess(np.abs(r.U).max(), 1e-3)

    def test_lu_zero_pivots_zero_on_regular_matrix(self):
        d, _ = plate(clamp)
        a = Assembler(d); a.assemble_system()
        K, _, _, _ = a.reduce(a.K_global.tocsr(), np.zeros(d.total_dofs))
        self.assertEqual(LUSolver().factorize(K).n_zero_pivots, 0)


# ----------------------------------------------------------------------
# 3. Diagnóstico del no lineal
# ----------------------------------------------------------------------

class TestNonlinearDiagnosis(unittest.TestCase):

    def test_internal_mechanism_reported_as_singular_tangent(self):
        """Antes: ``UnknownDivergenceError`` ("modo no clasificado"). Ahora el
        corrector registra que el sistema tangente no se resolvió con
        precisión y la divergencia se clasifica como tangente singular."""
        for la in DIRECT:
            with self.subTest(backend=la):
                d, ids = hinged_blocks()
                with self.assertRaises(SingularTangentError) as cm:
                    run_static(d, [(ids[(4.0, 2.0)], "uy", -1e3)],
                               lambda a: NonlinearSolver(a, num_steps=1, linear_algebra=la))
                self.assertIn("mecanismo", str(cm.exception))


# ----------------------------------------------------------------------
# 4. Memoria
# ----------------------------------------------------------------------

class TestOutOfMemoryMessage(unittest.TestCase):

    def test_superlu_out_of_memory_suggests_iterative(self):
        d, _ = plate(clamp)
        a = Assembler(d); a.assemble_system()
        K, _, _, _ = a.reduce(a.K_global.tocsr(), np.zeros(d.total_dofs))
        with mock.patch("solidum.math.linalg.lu.spla.splu", side_effect=MemoryError):
            with self.assertRaises(MemoryError) as cm:
                LUSolver().factorize(K)
        self.assertIn("linear_algebra: iterative", str(cm.exception))

    @unittest.skipUnless(_HAS_PARDISO, "pypardiso no instalado")
    def test_pardiso_out_of_memory_suggests_iterative(self):
        from pypardiso.pardiso_wrapper import PyPardisoError
        from solidum.math.linalg import PardisoSolver  # type: ignore[attr-defined]
        d, _ = plate(clamp)
        a = Assembler(d); a.assemble_system()
        K, _, _, _ = a.reduce(a.K_global.tocsr(), np.zeros(d.total_dofs))
        handle = PardisoSolver._handle()
        with mock.patch.object(handle, "factorize", side_effect=PyPardisoError(-2)):
            with self.assertRaises(MemoryError) as cm:
                PardisoSolver().factorize(K)
        self.assertIn("linear_algebra: iterative", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
