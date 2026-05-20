"""Tests de restricciones afines lineales (MPC), ADR 0004 fase 2.

Cubre:

1. **Cierre transitivo** de cadenas master-slave (`a ← b`, `b ← c` ⇒ `a ← c`).
2. **Detección de ciclos** (`a ← b`, `b ← a`).
3. **Apoyo en plano oblicuo**: rodillo a 45° en armadura — comparación con la
   misma armadura rotada y empotrada en eje global.
4. **Periodicidad**: dos barras paralelas con igualdad de desplazamientos en
   los nodos de un extremo, bajo carga en uno solo. La carga se reparte por
   simetría.
"""

import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from solidum.bc.constraints import ConstraintSet
from solidum.core.domain import Domain
from solidum.elements.frame import Frame2DEuler
from solidum.elements.truss import Truss2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver
from solidum.results import build_solve_result


class DummyMaterial:
    STRAIN_DIM = 1

    def __init__(self, E):
        self.E = E

    def compute_state(self, strain, state_var=None):
        return self.E * strain, self.E, state_var


class TestTransitiveClosure(unittest.TestCase):
    """Cierre transitivo de cadenas master-slave."""

    def test_chain_a_b_c(self):
        """``a ← b``, ``b ← c`` ⇒ tras cierre, ``a`` y ``b`` apuntan a ``c``."""
        cs = ConstraintSet()
        cs.add_linear(0, masters=[1], coefficients=[1.0])  # a = b
        cs.add_linear(1, masters=[2], coefficients=[1.0])  # b = c
        ndof = 3

        T, g = cs.build(ndof)

        # Tras cierre: a y b dependen ambos de c (única columna libre).
        self.assertEqual(T.shape, (3, 1))
        T_dense = T.toarray()
        np.testing.assert_array_equal(T_dense, [[1.0], [1.0], [1.0]])
        np.testing.assert_array_equal(g, [0.0, 0.0, 0.0])

    def test_chain_with_offsets(self):
        """``a ← 2·b + 1``, ``b ← 3·c + 0.5`` ⇒ ``a ← 6·c + 2``."""
        cs = ConstraintSet()
        cs.add_linear(0, masters=[1], coefficients=[2.0], g=1.0)
        cs.add_linear(1, masters=[2], coefficients=[3.0], g=0.5)
        T, g = cs.build(ndof=3)

        T_dense = T.toarray()
        # Columna única corresponde al DOF libre 2.
        np.testing.assert_array_equal(T_dense[:, 0], [6.0, 3.0, 1.0])
        np.testing.assert_allclose(g, [2.0, 0.5, 0.0])

    def test_cycle_detected(self):
        """``a ← b`` y ``b ← a`` debe lanzar al construir."""
        cs = ConstraintSet()
        cs.add_linear(0, masters=[1], coefficients=[1.0])
        cs.add_linear(1, masters=[0], coefficients=[1.0])
        with self.assertRaises(ValueError):
            cs.build(ndof=2)

    def test_self_reference_direct(self):
        """``add_linear`` rechaza maestro = esclavo de forma directa."""
        cs = ConstraintSet()
        with self.assertRaises(ValueError):
            cs.add_linear(0, masters=[0], coefficients=[1.0])

    def test_inconsistent_redeclare(self):
        """Reclarar el mismo esclavo con coeficientes distintos lanza."""
        cs = ConstraintSet()
        cs.add_linear(0, masters=[1], coefficients=[1.0])
        with self.assertRaises(ValueError):
            cs.add_linear(0, masters=[1], coefficients=[2.0])


class TestObliqueRoller(unittest.TestCase):
    """Apoyo en plano oblicuo en el extremo libre de una barra horizontal.

    Geometría: n1=(0,0) empotrado, n2=(1,0) con rodillo a 45°. La restricción
    ``u_y = u_x`` permite el movimiento a lo largo de la dirección
    ``[1, 1]/√2`` y bloquea la perpendicular ``[1, −1]/√2``.

    Caso 1 — carga horizontal pura: la barra resiste solo en su eje;
    ``u_x = PL/EA``, ``u_y`` se ajusta a ``u_x`` por la restricción.

    Caso 2 — carga perpendicular al plano del rodillo (``[P, −P]``): la
    restricción la absorbe íntegramente y los desplazamientos son nulos.
    """

    def _build_horizontal_bar_with_oblique_roller(self):
        domain = Domain()
        n1 = domain.add_node(1, [0.0, 0.0])
        n2 = domain.add_node(2, [1.0, 0.0])

        E, A = 1.0e6, 1.0e-3
        mat = DummyMaterial(E=E)
        domain.add_element(Truss2D(1, [n1, n2], mat, A=A))

        n1.fix_dof("ux", 0.0)
        n1.fix_dof("uy", 0.0)

        # Rodillo a 45°: u_y = u_x (esclavo = u_y de n2, maestro = u_x de n2).
        domain.add_linear_constraint(
            slave=(2, "uy"),
            masters=[(2, "ux")],
            coefficients=[1.0],
            g=0.0,
        )
        domain.generate_equation_numbers()
        return domain, n2, E, A

    def test_horizontal_load_extends_bar(self):
        domain, n2, E, A = self._build_horizontal_bar_with_oblique_roller()
        L, P = 1.0, 1.0
        F = np.zeros(domain.total_dofs)
        F[n2.dofs["ux"]] = P

        U = LinearSolver(Assembler(domain)).solve(F)

        expected = P * L / (E * A)
        self.assertAlmostEqual(U[n2.dofs["ux"]], expected, places=10)
        # Restricción cumplida a redondeo.
        self.assertAlmostEqual(U[n2.dofs["uy"]], U[n2.dofs["ux"]], places=12)

    def test_load_perpendicular_to_roller_plane_absorbed(self):
        domain, n2, _, _ = self._build_horizontal_bar_with_oblique_roller()
        F = np.zeros(domain.total_dofs)
        P = 1.0
        # Carga en dirección [1, -1] (perpendicular al plano del rodillo).
        F[n2.dofs["ux"]] = P
        F[n2.dofs["uy"]] = -P

        U = LinearSolver(Assembler(domain)).solve(F)

        # La restricción absorbe toda la carga: desplazamientos nulos.
        self.assertAlmostEqual(U[n2.dofs["ux"]], 0.0, places=12)
        self.assertAlmostEqual(U[n2.dofs["uy"]], 0.0, places=12)


class TestPeriodicity(unittest.TestCase):
    """Periodicidad ``u_x_n2 = u_x_n4`` en barras paralelas — la carga
    se distribuye por simetría.

    Geometría:
        n1=(0,0) — n2=(1,0)   bar 1 (área A)
        n3=(0,1) — n4=(1,1)   bar 2 (área A)

    n1 y n3 empotrados (ux,uy). En n4 también uy=0.
    Carga horizontal P aplicada en n2; restricción periódica
    ``u_x_n2 = u_x_n4``.

    Sin restricción la barra 1 toma toda P; con la restricción ambas se
    estiran lo mismo y reparten la carga ⇒ ``u_x_n2 = P·L / (2·E·A)``.
    """

    def test_load_split_by_periodicity(self):
        domain = Domain()
        n1 = domain.add_node(1, [0.0, 0.0])
        n2 = domain.add_node(2, [1.0, 0.0])
        n3 = domain.add_node(3, [0.0, 1.0])
        n4 = domain.add_node(4, [1.0, 1.0])

        E, A = 1.0e6, 1.0e-3
        mat = DummyMaterial(E=E)
        domain.add_element(Truss2D(1, [n1, n2], mat, A=A))
        domain.add_element(Truss2D(2, [n3, n4], mat, A=A))

        # Empotramientos.
        for dof in ("ux", "uy"):
            n1.fix_dof(dof, 0.0)
            n3.fix_dof(dof, 0.0)
        # n2 y n4 no resisten carga vertical (truss horizontal): pinear uy.
        n2.fix_dof("uy", 0.0)
        n4.fix_dof("uy", 0.0)

        # Periodicidad: u_x de n4 igual al de n2 (sin offset).
        domain.add_linear_constraint(
            slave=(4, "ux"),
            masters=[(2, "ux")],
            coefficients=[1.0],
            g=0.0,
        )

        domain.generate_equation_numbers()

        P = 1.0
        F = np.zeros(domain.total_dofs)
        F[n2.dofs["ux"]] = P
        # Sin carga en n4 — la carga viaja vía la restricción.

        assembler = Assembler(domain)
        U = LinearSolver(assembler).solve(F)

        u_x2 = U[n2.dofs["ux"]]
        u_x4 = U[n4.dofs["ux"]]
        # Restricción exacta a redondeo.
        self.assertAlmostEqual(u_x2, u_x4, places=12)

        # Reparto: cada barra toma P/2.
        L = 1.0
        expected = (P / 2.0) * L / (E * A)
        self.assertAlmostEqual(u_x2, expected, places=10)

        # Reacciones: n1 y n3 deben absorber P/2 cada uno (compresión axial).
        result = build_solve_result(domain, assembler, U, F)
        self.assertAlmostEqual(result.reactions_by_node[1]["ux"], -P / 2.0, places=10)
        self.assertAlmostEqual(result.reactions_by_node[3]["ux"], -P / 2.0, places=10)


class TestRigidLink(unittest.TestCase):
    """Unión rígida master-slave entre dos nodos en 2D.

    La cinemática del cuerpo rígido en el plano fija el desplazamiento de
    un nodo esclavo en función del desplazamiento y la rotación de un nodo
    maestro. Para un offset ``r = (r_x, r_y) = pos_slave − pos_master``:

        u_x_slave = u_x_master − r_y · θ_master
        u_y_slave = u_y_master + r_x · θ_master
        θ_slave   = θ_master

    Caso de prueba: cantilever Bernoulli-Euler de longitud ``L`` con un
    nodo satélite a offset ``(0, 0.5)`` rígidamente unido al extremo libre.
    Bajo momento puro ``M`` aplicado al nodo maestro, la solución analítica
    del cantilever es ``θ = ML/EI``, ``v = ML²/(2EI)`` en el extremo. El
    nodo satélite debe seguir la cinemática rígida.
    """

    def test_rigid_link_with_rotation(self):
        L, EI, EA, M = 1.0, 1.0, 1.0, 1.0
        offset_x, offset_y = 0.0, 0.5

        domain = Domain()
        n1 = domain.add_node(1, [0.0, 0.0])  # empotrado
        n2 = domain.add_node(2, [L, 0.0])    # extremo libre — maestro
        n3 = domain.add_node(3, [L + offset_x, offset_y])  # satélite — esclavo

        # El frame solo une n1-n2; n3 no tiene elemento, sus DOFs los fija
        # íntegramente la unión rígida — los registramos a mano porque
        # ``_register_dofs`` solo lo hacen los elementos.
        mat = DummyMaterial(E=EA)
        domain.add_element(Frame2DEuler(1, [n1, n2], mat, A=EA, I=EI))
        for dof in ("ux", "uy", "rz"):
            n3.add_dof(dof)

        # Empotramiento total en n1.
        for dof in ("ux", "uy", "rz"):
            n1.fix_dof(dof, 0.0)

        # Unión rígida n3 ← n2:
        #   u_x_n3 = u_x_n2 − r_y · rz_n2
        #   u_y_n3 = u_y_n2 + r_x · rz_n2
        #   rz_n3  = rz_n2
        domain.add_linear_constraint(
            slave=(3, "ux"),
            masters=[(2, "ux"), (2, "rz")],
            coefficients=[1.0, -offset_y],
            g=0.0,
        )
        domain.add_linear_constraint(
            slave=(3, "uy"),
            masters=[(2, "uy"), (2, "rz")],
            coefficients=[1.0, offset_x],
            g=0.0,
        )
        domain.add_linear_constraint(
            slave=(3, "rz"),
            masters=[(2, "rz")],
            coefficients=[1.0],
            g=0.0,
        )

        domain.generate_equation_numbers()
        F = np.zeros(domain.total_dofs)
        F[n2.dofs["rz"]] = M  # momento puro en el extremo libre

        U = LinearSolver(Assembler(domain)).solve(F)

        # Solución analítica del cantilever Bernoulli-Euler bajo momento puro.
        theta_n2 = M * L / EI
        v_n2 = M * L**2 / (2.0 * EI)
        u_x_n2 = 0.0  # sin carga axial

        self.assertAlmostEqual(U[n2.dofs["ux"]], u_x_n2, places=10)
        self.assertAlmostEqual(U[n2.dofs["uy"]], v_n2, places=10)
        self.assertAlmostEqual(U[n2.dofs["rz"]], theta_n2, places=10)

        # Cinemática rígida en el satélite.
        expected_ux_n3 = u_x_n2 - offset_y * theta_n2
        expected_uy_n3 = v_n2 + offset_x * theta_n2
        expected_rz_n3 = theta_n2
        self.assertAlmostEqual(U[n3.dofs["ux"]], expected_ux_n3, places=12)
        self.assertAlmostEqual(U[n3.dofs["uy"]], expected_uy_n3, places=12)
        self.assertAlmostEqual(U[n3.dofs["rz"]], expected_rz_n3, places=12)


class TestYamlLinearConstraints(unittest.TestCase):
    """El parser YAML reconoce el bloque ``linear_constraints`` y construye
    las MPC equivalentes a la API programática."""

    def test_yaml_periodicity(self):
        import os
        import tempfile

        yaml_text = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [0.0, 1.0]}
  - {id: 4, coords: [1.0, 1.0]}

materials:
  - {id: 1, type: Elastic1D, E: 1.0e+6, density: 0.0}

elements:
  - {id: 1, type: Truss2D, material: 1, nodes: [1, 2], A: 1.0e-3}
  - {id: 2, type: Truss2D, material: 1, nodes: [3, 4], A: 1.0e-3}

boundary_conditions:
  - {node_id: 1, ux: 0.0, uy: 0.0}
  - {node_id: 3, ux: 0.0, uy: 0.0}
  - {node_id: 2, uy: 0.0}
  - {node_id: 4, uy: 0.0}

linear_constraints:
  - slave: {node: 4, dof: ux}
    masters:
      - {node: 2, dof: ux}
    coefficients: [1.0]

point_loads:
  - {node_id: 2, ux: 1.0}

solver:
  type: LinearSolver
"""
        from solidum.utils.yaml_parser import YamlParser

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False, encoding="utf-8"
        ) as fh:
            fh.write(yaml_text)
            path = fh.name
        try:
            parser = YamlParser(path)
            domain = parser.parse()
            self.assertEqual(len(domain.linear_constraints), 1)
            self.assertEqual(domain.linear_constraints[0]["slave"], (4, "ux"))
            self.assertEqual(domain.linear_constraints[0]["masters"], [(2, "ux")])
            self.assertEqual(domain.linear_constraints[0]["coefficients"], [1.0])
            self.assertEqual(domain.linear_constraints[0]["g"], 0.0)
        finally:
            os.unlink(path)


if __name__ == "__main__":
    unittest.main()
