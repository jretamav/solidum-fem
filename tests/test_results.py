"""Tests de los tipos de dato de solidum.results (ADR 0002)."""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.results import ElementForces, SolveResult


class TestElementForces(unittest.TestCase):
    def test_truss_ok(self):
        ef = ElementForces(kind="truss", components={"N": np.array([1.0, 1.0])})
        self.assertEqual(ef.at_node_i(), {"N": 1.0})
        self.assertEqual(ef.at_node_j(), {"N": 1.0})

    def test_frame2d_ok(self):
        comps = {
            "N": np.array([10.0, 10.0]),
            "V": np.array([2.0, -2.0]),
            "M": np.array([0.0, 5.0]),
        }
        ef = ElementForces(kind="frame2d", components=comps)
        self.assertEqual(ef.at_node_j()["M"], 5.0)

    def test_frame3d_keys(self):
        comps = {k: np.array([0.0, 0.0]) for k in ("N", "Vy", "Vz", "T", "My", "Mz")}
        ef = ElementForces(kind="frame3d", components=comps)
        self.assertEqual(set(ef.components), {"N", "Vy", "Vz", "T", "My", "Mz"})

    def test_unknown_kind_rejected(self):
        with self.assertRaises(ValueError):
            ElementForces(kind="shell", components={})  # type: ignore[arg-type]

    def test_invalid_component_key_rejected(self):
        with self.assertRaises(ValueError) as cm:
            ElementForces(kind="truss", components={"M": np.array([0.0, 0.0])})
        self.assertIn("Componentes inválidos", str(cm.exception))

    def test_wrong_shape_rejected(self):
        with self.assertRaises(ValueError):
            ElementForces(kind="truss", components={"N": np.array([1.0, 2.0, 3.0])})

    def test_frozen(self):
        ef = ElementForces(kind="truss", components={"N": np.array([1.0, 1.0])})
        with self.assertRaises(Exception):
            ef.kind = "cable"  # type: ignore[misc]


class TestSolveResult(unittest.TestCase):
    def test_defaults(self):
        n = 6
        res = SolveResult(
            U=np.zeros(n),
            F_applied=np.zeros(n),
            R=np.zeros(n),
        )
        self.assertTrue(res.converged)
        self.assertEqual(res.num_steps, 1)
        self.assertEqual(res.reactions_by_node, {})
        self.assertEqual(res.element_forces, {})

    def test_frozen(self):
        res = SolveResult(U=np.zeros(2), F_applied=np.zeros(2), R=np.zeros(2))
        with self.assertRaises(Exception):
            res.converged = False  # type: ignore[misc]


class TestTransientResultInternalForcesHistory(unittest.TestCase):
    """``TransientResult.internal_forces_history(domain)`` debe devolver
    la historia de ``ElementForces`` para cada elemento que implemente el
    contrato ADR 0002 (regresión auditoría H-1.3).

    Caso de validación: oscilador 1 GDL no amortiguado con material
    elástico lineal. El esfuerzo axial en cada paso ``k`` debe coincidir
    con ``E·A·(u_2(t_k) − u_1(t_k))/L``, y la trayectoria de la fuerza
    axial debe replicar la del desplazamiento del DOF libre porque la
    relación es lineal y los apoyos son fijos.
    """

    def test_truss_internal_forces_track_displacement(self):
        import math

        import solidum  # autodiscover
        from solidum.core.domain import Domain
        from solidum.elements.truss import Truss2D
        from solidum.entry import run_transient
        from solidum.materials.elastic import Elastic1D

        # Parámetros 1 GDL: ω = 5 rad/s con K=25, M=1.
        E, A, L = 25.0, 1.0, 1.0
        rho_lumped = 2.0
        omega = 5.0

        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [L, 0.0])
        mat = Elastic1D(E=E, density=rho_lumped)
        dom.add_element(Truss2D(1, [n1, n2], mat, A=A))
        n1.fix_dof("ux", 0.0); n1.fix_dof("uy", 0.0)
        n2.fix_dof("uy", 0.0)
        dom.generate_equation_numbers(verbose=False)

        T = 2.0 * math.pi / omega
        dt = T / 200.0
        u0 = np.zeros(dom.total_dofs); u0[n2.dofs["ux"]] = 0.1
        result = run_transient(dom, t_end=T, dt=dt, u0=u0)

        history = result.internal_forces_history(dom)

        # Truss2D registrado con elem_id=1.
        self.assertIn(1, history)
        per_step = history[1]
        # n_steps + 1 entradas (incluye t=0).
        self.assertEqual(len(per_step), result.u_history.shape[1])

        # ElementForces type para Truss tiene componente "N".
        self.assertEqual(per_step[0].kind, "truss")
        self.assertIn("N", per_step[0].components)

        # Validación física: para Truss2D, N = E·A·(u_2x − u_1x)/L₀,
        # con el primer nodo fijado en x=0 ⇒ N(t) = E·A·u_2x(t)/L₀.
        u2x = result.u_history[n2.dofs["ux"], :]
        N_expected = E * A * u2x / L
        N_history = np.array([ef.components["N"][0] for ef in per_step])
        np.testing.assert_allclose(N_history, N_expected, rtol=1e-10)


if __name__ == "__main__":
    unittest.main()
