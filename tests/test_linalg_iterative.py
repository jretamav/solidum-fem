"""Backend algebraico iterativo (ADR 0018).

Cubre, en este orden:

1. **Modos de cuerpo rígido**: son el núcleo exacto de ``K`` en modelos sin
   apoyos, en todas las familias del catálogo (sólidos 2D/3D, armaduras,
   marcos 2D/3D, térmico y dominio mixto), y el ensamblador los restringe a
   los DOF libres y refresca la caché si cambian los apoyos.
2. **Backend**: equivalencia con el directo dentro de la tolerancia, residuo
   verdadero bajo ``rtol``, factorización reutilizable, curvatura negativa →
   MINRES, diagonal no positiva → MINRES directo, no convergencia → error
   explícito (nunca una solución parcial), bordes ``0×0`` y ``b = 0``.
3. **Despachador y YAML**: sintaxis ``iterative[:precondicionador]``, matriz
   no simétrica → solver directo con aviso, valores inválidos rechazados.
4. **Corrector**: un iterativo sin convergencia se reporta con su causa, no
   como tangente singular a secas.
5. **Pipelines**: ``LinearSolver``, ``NonlinearSolver`` (J2),
   ``ArcLengthSolver``, ``NewmarkSolver``, ``NewtonNewmarkSolver`` y
   ``ThetaMethodSolver`` dan con ``linear_algebra='iterative'`` el mismo
   resultado que con el directo.

Los casos con AMG se saltan si ``pyamg`` no está instalado; el resto corre
siempre (CG sin precondicionar o con Jacobi no tiene dependencias).
"""
import logging
import os
import sys
import tempfile
import unittest
import warnings
from unittest import mock

import numpy as np
import scipy.sparse as sp

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import solidum
from solidum import Assembler, Domain
from solidum.constants import AMG_MAX_COARSE, ITERATIVE_RTOL
from solidum.materials.thermal_conduction import ThermalConduction
from solidum.math.linalg import (
    IterativeNotConvergedError,
    IterativeSolver,
    LUSolver,
    StiffnessProperties,
    rigid_body_modes,
    select_solver,
)
from solidum.math.linalg.dispatcher import available_overrides, is_valid_override
from solidum.math.linalg.iterative import HAS_PYAMG, IterativeFactorized
from solidum.registry import ElementRegistry, MaterialRegistry

EL = ElementRegistry.get
Elastic1D = MaterialRegistry.get("Elastic1D")
Elastic2D = MaterialRegistry.get("Elastic2D")
Elastic3D = MaterialRegistry.get("Elastic3D")
VonMises2D = MaterialRegistry.get("VonMises2D")


# Los modelos de prueba emiten avisos de construcción irrelevantes aquí
# (``nu`` por defecto en marcos, barra vertical sin ``ref_vector``). Se
# silencian SÓLO mientras corre este módulo: un ``logging.disable`` a nivel
# de módulo se ejecutaría al recolectar y apagaría el logging de toda la
# sesión, rompiendo las pruebas de otros archivos que verifican avisos.
def setUpModule():
    logging.disable(logging.WARNING)


def tearDownModule():
    logging.disable(logging.NOTSET)


# ----------------------------------------------------------------------
# Modelos de prueba
# ----------------------------------------------------------------------

def quad4_grid(nx, ny, material=None, *, lx=4.0, ly=1.0, fix_left=True, density=None):
    """Placa ``nx × ny`` de Quad4 ligeramente distorsionada (ningún elemento
    es un rectángulo exacto), empotrada en ``x = 0`` si ``fix_left``."""
    if material is None:
        material = Elastic2D(E=210e9, nu=0.3, density=density)
    d = Domain()
    nid = 1
    for j in range(ny + 1):
        for i in range(nx + 1):
            x = lx * i / nx + 0.02 * lx / nx * np.sin(1.7 * j)
            y = ly * j / ny + 0.02 * ly / ny * np.cos(1.3 * i)
            d.add_node(nid, [x, y]); nid += 1
    eid = 1
    for j in range(ny):
        for i in range(nx):
            n0 = 1 + i + j * (nx + 1)
            nodes = [d.nodes[n0], d.nodes[n0 + 1], d.nodes[n0 + nx + 2], d.nodes[n0 + nx + 1]]
            d.add_element(EL("Quad4")(eid, nodes, material)); eid += 1
    if fix_left:
        for j in range(ny + 1):
            node = d.nodes[1 + j * (nx + 1)]
            node.fix_dof("ux", 0.0); node.fix_dof("uy", 0.0)
    d.generate_equation_numbers()
    return d


def right_edge_nodes(d, nx, ny):
    return [d.nodes[(nx + 1) * (j + 1)] for j in range(ny + 1)]


def reduced_stiffness(d):
    a = Assembler(d)
    a.assemble_system()
    K, _, _, _ = a.reduce(a.K_global.tocsr(), np.zeros(d.total_dofs))
    return a, K.tocsr()


def free_floating_models():
    """Un modelo sin apoyos por familia del catálogo, con su número de modos
    rígidos esperado."""
    out = []

    d = quad4_grid(3, 2, fix_left=False)
    out.append(("Quad4", d, 3))

    d = Domain(); nid = 1; ids = {}
    for k in range(2):
        for j in range(2):
            for i in range(3):
                d.add_node(nid, [i + 0.1 * j, j + 0.05 * k, k + 0.1 * i]); ids[(i, j, k)] = nid; nid += 1
    m3 = Elastic3D(E=210e9, nu=0.3)
    for e, i in enumerate(range(2), 1):
        c = [(i, 0, 0), (i + 1, 0, 0), (i + 1, 1, 0), (i, 1, 0),
             (i, 0, 1), (i + 1, 0, 1), (i + 1, 1, 1), (i, 1, 1)]
        d.add_element(EL("Hex8")(e, [d.nodes[ids[q]] for q in c], m3))
    out.append(("Hex8", d, 6))

    m1 = Elastic1D(E=200e9)
    d = Domain()
    for i, p in enumerate([(0, 0), (3, 0), (1.5, 2)], 1):
        d.add_node(i, list(p))
    for e, (a_, b_) in enumerate([(1, 2), (2, 3), (3, 1)], 1):
        d.add_element(EL("Truss2D")(e, [d.nodes[a_], d.nodes[b_]], m1, A=1e-3))
    out.append(("Truss2D", d, 3))

    d = Domain()
    for i, p in enumerate([(0, 0, 0), (2, 0, 0), (0, 2, 0), (0.3, 0.4, 2)], 1):
        d.add_node(i, list(p))
    for e, (a_, b_) in enumerate([(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)], 1):
        d.add_element(EL("Truss3D")(e, [d.nodes[a_], d.nodes[b_]], m1, A=1e-3))
    out.append(("Truss3D", d, 6))

    for name, kw in (("Frame2DEuler", dict(A=1e-2, I=1e-4)),
                     ("Frame2DTimoshenko", dict(A=1e-2, I=1e-4, As=8e-3, nu=0.3))):
        d = Domain()
        for i, p in enumerate([(0, 0), (0, 3), (4, 3.5), (4, 0)], 1):
            d.add_node(i, list(p))
        for e, (a_, b_) in enumerate([(1, 2), (2, 3), (3, 4)], 1):
            d.add_element(EL(name)(e, [d.nodes[a_], d.nodes[b_]], m1, **kw))
        out.append((name, d, 3))

    d = Domain()
    for i, p in enumerate([(0, 0, 0), (0.2, 0, 3), (4, 0.5, 3), (4, 2, 0)], 1):
        d.add_node(i, list(p))
    for e, (a_, b_) in enumerate([(1, 2), (2, 3), (3, 4)], 1):
        d.add_element(EL("Frame3D")(e, [d.nodes[a_], d.nodes[b_]], m1,
                                    A=1e-2, Iy=1e-4, Iz=2e-4, J=5e-5, nu=0.3,
                                    ref_vector=[0.0, 1.0, 0.0] if e != 2 else [0.0, 0.0, 1.0]))
    out.append(("Frame3D", d, 6))

    d = Domain()
    for i, p in enumerate([(0, 0), (1, 0), (1.1, 1), (0, 1)], 1):
        d.add_node(i, list(p))
    d.add_element(EL("Quad4Thermal")(1, [d.nodes[k] for k in (1, 2, 3, 4)], ThermalConduction(k=50.0)))
    out.append(("Quad4Thermal", d, 1))

    d = Domain()
    for i, p in enumerate([(0, 0), (1, 0), (1.1, 1), (0, 1)], 1):
        d.add_node(i, list(p))
    ns = [d.nodes[k] for k in (1, 2, 3, 4)]
    d.add_element(EL("Quad4")(1, ns, Elastic2D(E=210e9, nu=0.3)))
    d.add_element(EL("Quad4Thermal")(2, ns, ThermalConduction(k=50.0)))
    out.append(("Quad4 + Quad4Thermal", d, 4))
    return out


def _rel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


# ----------------------------------------------------------------------
# 1. Modos de cuerpo rígido
# ----------------------------------------------------------------------

class TestRigidBodyModes(unittest.TestCase):

    def test_modes_are_exact_nullspace_of_free_floating_models(self):
        """``K·B = 0`` a precisión de máquina en cada familia: son los modos
        rígidos exactos, no una aproximación."""
        for label, d, n_modes in free_floating_models():
            with self.subTest(modelo=label):
                if d.total_dofs == 0:
                    d.generate_equation_numbers()
                a = Assembler(d)
                a.assemble_system()
                K = a.K_global.tocsr()
                B = rigid_body_modes(d)
                self.assertEqual(B.shape, (d.total_dofs, n_modes))
                self.assertEqual(np.linalg.matrix_rank(B), n_modes)
                num = np.linalg.norm(K @ B, axis=0)
                den = abs(K).max() * np.linalg.norm(B, axis=0)
                self.assertLess(float(np.max(num / den)), 1e-13)

    def test_mixed_domain_blocks_do_not_mix(self):
        """En un dominio u + T los modos mecánicos valen cero en ``T`` y el
        modo constante térmico vale cero en ``u``."""
        _, d, _ = free_floating_models()[-1]
        B = rigid_body_modes(d)
        T_rows = [n.dofs["T"] for n in d.nodes.values()]
        u_rows = [n.dofs[k] for n in d.nodes.values() for k in ("ux", "uy")]
        for j in range(B.shape[1]):
            on_T = np.linalg.norm(B[T_rows, j]) > 0
            on_u = np.linalg.norm(B[u_rows, j]) > 0
            self.assertNotEqual(on_T, on_u, f"el modo {j} mezcla campos")

    def test_assembler_restricts_to_free_dofs_and_refreshes_cache(self):
        d = quad4_grid(6, 2)
        a = Assembler(d)
        a.assemble_system()
        free = a.constraint_set.free_dofs(d.total_dofs)
        B = a.near_nullspace()
        np.testing.assert_array_equal(B, rigid_body_modes(d)[free])
        self.assertIs(a.near_nullspace(), B)  # cacheado

        d.nodes[2].fix_dof("uy", 0.0)          # nuevo apoyo
        B2 = a.near_nullspace()
        self.assertEqual(B2.shape[0], B.shape[0] - 1)


# ----------------------------------------------------------------------
# 2. Backend
# ----------------------------------------------------------------------

class TestIterativeBackend(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # 30 × 10 Quad4 empotrada: 660 incógnitas, por encima de
        # AMG_MAX_COARSE, así que AMG construye de verdad una jerarquía
        # (por debajo sería un solve directo y CG convergería en 1 it.).
        cls.domain = quad4_grid(30, 10)
        cls.assembler, cls.K = reduced_stiffness(cls.domain)
        assert cls.K.shape[0] > AMG_MAX_COARSE
        rng = np.random.default_rng(0)
        cls.b = rng.standard_normal(cls.K.shape[0])
        cls.x_ref = LUSolver().solve(cls.K, cls.b)
        cls.props = StiffnessProperties(
            is_symmetric=True, is_positive_definite=True, size=cls.K.shape[0],
            near_nullspace=cls.assembler.near_nullspace,
        )

    def _check(self, preconditioner):
        f = IterativeSolver(self.props, preconditioner=preconditioner).factorize(self.K)
        x = f.solve(self.b)
        true_res = np.linalg.norm(self.b - self.K @ x) / np.linalg.norm(self.b)
        self.assertLessEqual(true_res, ITERATIVE_RTOL)
        self.assertAlmostEqual(f.last_relative_residual, true_res, delta=1e-14)
        self.assertLess(_rel(x, self.x_ref), 1e-7)
        self.assertEqual(f.method, "cg")
        return f

    def test_cg_without_preconditioner(self):
        self._check("none")

    def test_cg_jacobi(self):
        self._check("jacobi")

    @unittest.skipUnless(HAS_PYAMG, "pyamg no instalado")
    def test_cg_amg_with_rigid_body_modes_needs_few_iterations(self):
        f_amg = self._check("amg")
        f_none = IterativeSolver(self.props, preconditioner="none").factorize(self.K)
        f_none.solve(self.b)
        self.assertLess(f_amg.last_iterations, 40)
        self.assertLess(f_amg.last_iterations, f_none.last_iterations / 5)

    @unittest.skipUnless(HAS_PYAMG, "pyamg no instalado")
    def test_amg_without_nullspace_still_correct(self):
        """Sin proveedor de modos AMG usa su casi-núcleo por defecto: converge
        peor, pero la solución sigue cumpliendo la tolerancia."""
        props = StiffnessProperties(True, True, self.K.shape[0])
        x = IterativeSolver(props, preconditioner="amg").solve(self.K, self.b)
        self.assertLess(_rel(x, self.x_ref), 1e-7)

    def test_factorization_is_reusable(self):
        f = IterativeSolver(self.props, preconditioner="jacobi").factorize(self.K)
        for scale in (1.0, -3.0, 0.5):
            np.testing.assert_allclose(f.solve(scale * self.b), scale * self.x_ref,
                                       rtol=1e-7, atol=1e-7 * np.abs(self.x_ref).max())

    def test_negative_curvature_switches_to_minres(self):
        """Matriz simétrica con dos autovalores negativos (una tangente que
        acaba de cruzar un punto límite): CG detecta ``pᵀKp ≤ 0`` y el
        backend termina con MINRES, con la solución correcta."""
        lam = np.linalg.eigvalsh(self.K.toarray())[:3]
        self.assertGreater(lam[2] - lam[1], 1e-3 * lam[2])  # autovalores distintos
        Ks = (self.K - 0.5 * (lam[1] + lam[2]) * sp.eye(self.K.shape[0])).tocsr()
        x_ref = LUSolver().solve(Ks, self.b)
        for pc in ("none",) + (("amg",) if HAS_PYAMG else ()):
            with self.subTest(precondicionador=pc):
                f = IterativeSolver(self.props, preconditioner=pc).factorize(Ks)
                x = f.solve(self.b)
                self.assertEqual(f.method, "minres")
                self.assertLess(_rel(x, x_ref), 1e-6)
                self.assertLessEqual(f.last_relative_residual, ITERATIVE_RTOL)

    def test_nonpositive_diagonal_goes_straight_to_minres(self):
        f = IterativeSolver(self.props).factorize((-self.K).tocsr())
        self.assertEqual(f.method, "minres")
        x = f.solve(self.b)
        self.assertLess(_rel(x, -self.x_ref), 1e-7)

    def test_non_convergence_raises_never_returns_partial(self):
        f = IterativeSolver(self.props, preconditioner="none", max_iter=3).factorize(self.K)
        with self.assertRaises(IterativeNotConvergedError) as cm:
            f.solve(self.b)
        err = cm.exception
        self.assertGreater(err.relative_residual, err.rtol)
        self.assertEqual(err.size, self.K.shape[0])
        self.assertIn("Solver iterativo", str(err))
        self.assertIsInstance(err, RuntimeError)  # los solvers la tratan como fallo del sistema

    def test_edge_cases(self):
        self.assertEqual(IterativeSolver().factorize(sp.csr_matrix((0, 0))).solve(np.zeros(0)).shape, (0,))
        np.testing.assert_array_equal(
            IterativeSolver(self.props).factorize(self.K).solve(np.zeros(self.K.shape[0])),
            np.zeros(self.K.shape[0]))

    def test_invalid_arguments(self):
        with self.assertRaises(ValueError):
            IterativeSolver(preconditioner="ilu")
        with self.assertRaises(ValueError):
            IterativeSolver(rtol=0.0)
        with self.assertRaises(ValueError):
            IterativeSolver(max_iter=0)

    def test_amg_requires_pyamg(self):
        with mock.patch("solidum.math.linalg.iterative.HAS_PYAMG", False):
            with self.assertRaises(ValueError):
                IterativeSolver(preconditioner="amg")


# ----------------------------------------------------------------------
# 3. Despachador y YAML
# ----------------------------------------------------------------------

class TestDispatcherIterative(unittest.TestCase):

    def test_override_syntax(self):
        p = StiffnessProperties(True, True, 10)
        s = select_solver(p, override="iterative")
        self.assertIsInstance(s, IterativeSolver)
        self.assertEqual(s.preconditioner, "auto")
        s = select_solver(p, override="iterative:jacobi")
        self.assertEqual(s.preconditioner, "jacobi")
        self.assertIs(s.props, p)  # recibe las propiedades (casi-núcleo)

    def test_auto_never_selects_iterative(self):
        for sym, pd in ((True, True), (True, False), (False, False)):
            p = StiffnessProperties(sym, pd, 10)
            self.assertNotIsInstance(select_solver(p), IterativeSolver)

    def test_nonsymmetric_goes_to_direct_with_warning(self):
        p = StiffnessProperties(is_symmetric=False, is_positive_definite=False, size=10)
        with mock.patch("solidum.math.linalg.dispatcher._iterative_nonsym_warned", False):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                s = select_solver(p, override="iterative")
        self.assertNotIsInstance(s, IterativeSolver)
        self.assertTrue(any("no simétrica" in str(x.message) for x in w))

    def test_invalid_overrides(self):
        p = StiffnessProperties(True, True, 10)
        for bad in ("iterative:ilu", "pardiso:amg", "lu:jacobi", "cg"):
            with self.subTest(override=bad):
                self.assertFalse(is_valid_override(bad))
                with self.assertRaises(ValueError):
                    select_solver(p, override=bad)
        for good in ("auto", None, "lu", "iterative", "iterative:none", "iterative:jacobi"):
            self.assertTrue(is_valid_override(good), good)
        self.assertIn("iterative:jacobi", available_overrides())
        self.assertEqual("iterative:amg" in available_overrides(), HAS_PYAMG)

    def test_yaml_accepts_and_rejects(self):
        from solidum.utils.yaml_parser import YamlParser, YamlValidationError
        src = os.path.join(os.path.dirname(__file__), "..", "examples", "modelo_elastico_2d.yaml")
        with open(src, encoding="utf-8") as f:
            base = f.read()
        self.assertIn("\nsolver:\n", base)
        for value, ok in (("iterative:jacobi", True), ("iterative:foo", False)):
            with self.subTest(linear_algebra=value):
                text = base.replace("\nsolver:\n", f"\nsolver:\n  linear_algebra: {value}\n", 1)
                with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False,
                                                 encoding="utf-8") as tmp:
                    tmp.write(text)
                try:
                    if ok:
                        YamlParser(tmp.name).parse()
                    else:
                        with self.assertRaises(YamlValidationError):
                            YamlParser(tmp.name).parse()
                finally:
                    os.unlink(tmp.name)


# ----------------------------------------------------------------------
# 4. Corrector: diagnóstico honesto
# ----------------------------------------------------------------------

class TestCorrectorReportsIterativeFailure(unittest.TestCase):

    def test_divergence_error_carries_iterative_cause(self):
        from solidum.math.solvers.diagnostics import SolverDivergedError
        from solidum.math.solvers.nonlinear import NonlinearSolver

        d = quad4_grid(6, 2)
        a = Assembler(d)
        a.assemble_system()
        for n in right_edge_nodes(d, 6, 2):
            a.apply_point_load(n.id, "uy", -1.0e4)
        F = a.F_global.copy()

        def _fail(self, b):
            raise IterativeNotConvergedError(method="CG", preconditioner="none", iterations=7,
                                             relative_residual=1e-3, rtol=1e-10, size=10)

        solver = NonlinearSolver(a, num_steps=1, linear_algebra="iterative:none")
        with mock.patch.object(IterativeFactorized, "solve", _fail):
            with self.assertRaises(SolverDivergedError) as cm:
                solver.solve(F)
        self.assertIn("Solver iterativo", str(cm.exception))


# ----------------------------------------------------------------------
# 5. Pipelines completos: mismo resultado que el directo
# ----------------------------------------------------------------------

class TestPipelinesMatchDirect(unittest.TestCase):

    CHOICE = "iterative"  # precondicionador automático (AMG si hay pyamg)

    def _static(self, make_solver, material=None, load=-2.0e5, nx=20, ny=8):
        out = {}
        for la in ("auto", self.CHOICE):
            d = quad4_grid(nx, ny, material)
            a = Assembler(d)
            a.assemble_system()
            for n in right_edge_nodes(d, nx, ny):
                a.apply_point_load(n.id, "uy", load)
            F = a.F_global.copy()
            r = solidum.run(d, assembler=a, solver=make_solver(a, la), F_applied=F)
            self.assertTrue(r.converged)
            out[la] = r.U.copy()
        return _rel(out[self.CHOICE], out["auto"])

    def test_linear_solver(self):
        from solidum.math.solvers.linear import LinearSolver
        self.assertLess(self._static(lambda a, la: LinearSolver(a, linear_algebra=la)), 1e-8)

    def test_nonlinear_solver_j2_plasticity(self):
        from solidum.math.solvers.nonlinear import NonlinearSolver
        mat = VonMises2D(E=210e9, nu=0.3, sigma_y=250e6, H=2e9)
        diff = self._static(lambda a, la: NonlinearSolver(a, num_steps=4, linear_algebra=la),
                            material=mat, load=-1.5e6, nx=12, ny=4)
        self.assertLess(diff, 1e-6)

    def test_arclength_solver(self):
        from solidum.math.solvers.arclength import ArcLengthSolver
        diff = self._static(lambda a, la: ArcLengthSolver(a, max_lambda=1.0, initial_dl=0.5,
                                                          linear_algebra=la), nx=10, ny=4)
        self.assertLess(diff, 1e-6)

    def _newmark(self, cls, **kw):
        out = {}
        for la in ("auto", self.CHOICE):
            d = quad4_grid(12, 4, density=7850.0)
            a = Assembler(d)
            a.assemble_system()
            F = np.zeros(d.total_dofs)
            for n in right_edge_nodes(d, 12, 4):
                F[n.dofs["uy"]] = -1.0e4
            s = cls(a, t_end=5e-4, dt=1e-4, F_func=lambda t: F, linear_algebra=la, **kw)
            r = solidum.entry.run_transient(d, assembler=a, solver=s)
            out[la] = np.asarray(r.u_history)
        return _rel(out[self.CHOICE], out["auto"])

    def test_newmark_linear(self):
        from solidum.math.solvers.newmark import NewmarkSolver
        self.assertLess(self._newmark(NewmarkSolver), 1e-8)

    def test_newton_newmark(self):
        from solidum.math.solvers.newmark import NewtonNewmarkSolver
        self.assertLess(self._newmark(NewtonNewmarkSolver), 1e-6)

    def test_theta_method_thermal(self):
        from solidum.math.solvers.theta_method import ThetaMethodSolver
        out = {}
        for la in ("auto", self.CHOICE):
            d = Domain()
            nx, ny = 20, 8
            nid = 1
            for j in range(ny + 1):
                for i in range(nx + 1):
                    d.add_node(nid, [i / nx, 0.4 * j / ny]); nid += 1
            mat = ThermalConduction(k=50.0, c=460.0, density=7850.0)
            for j in range(ny):
                for i in range(nx):
                    n0 = 1 + i + j * (nx + 1)
                    nodes = [d.nodes[n0], d.nodes[n0 + 1], d.nodes[n0 + nx + 2], d.nodes[n0 + nx + 1]]
                    d.add_element(EL("Quad4Thermal")(len(d.elements) + 1, nodes, mat))
            for j in range(ny + 1):
                d.nodes[1 + j * (nx + 1)].fix_dof("T", 100.0)
            d.generate_equation_numbers()
            a = Assembler(d)
            s = ThetaMethodSolver(a, dt=50.0, n_steps=5, T_initial=0.0, linear_algebra=la)
            r = solidum.entry.run_thermal_transient(d, assembler=a, solver=s)
            out[la] = np.asarray(r.T_history)
        self.assertLess(_rel(out[self.CHOICE], out["auto"]), 1e-8)


if __name__ == "__main__":
    unittest.main()
