"""Regresiones del lote YAML de la auditoría 2026-09-22.

- ``*_by_coord`` funciona con nodos 3D e incluye ``z_min``/``z_max``.
- Cargas puntuales con nodo o DOF inexistente fallan con ``ValueError``
  (antes se ignoraban en silencio).
- Solvers transitorios reciben las cargas del YAML como escalón constante.
- ``mesh_material`` inexistente, ``mesh`` + ``nodes`` y kwargs desconocidos
  de materiales/solvers se rechazan en validación.
- El ``SafeLoader`` global de PyYAML no se modifica.
"""
import os
import tempfile
import unittest

import numpy as np
import yaml

import solidum
from solidum.utils.yaml_parser import YamlParser, YamlValidationError

HEX = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]


def _hex_yaml(extra: str) -> str:
    return "\n".join([
        "nodes:",
        *[f"  - {{id: {i}, coords: [{c[0]}, {c[1]}, {c[2]}]}}" for i, c in enumerate(HEX, 1)],
        "materials:",
        "  - {id: 1, type: Elastic3D, E: 1000.0, nu: 0.3, density: 1.0}",
        "elements:",
        "  - {id: 1, type: Hex8, material: 1, nodes: [1,2,3,4,5,6,7,8]}",
        extra,
        "",
    ])


BAR_2D = "\n".join([
    "nodes:",
    "  - {id: 1, coords: [0.0, 0.0]}",
    "  - {id: 2, coords: [1.0, 0.0]}",
    "materials:",
    "  - {id: 1, type: Elastic1D, E: 25.0, density: 2.0}",
    "elements:",
    "  - {id: 1, type: Truss2D, material: 1, A: 1.0, nodes: [1, 2]}",
    "boundary_conditions:",
    "  - {node_id: 1, ux: 0.0, uy: 0.0}",
    "  - {node_id: 2, uy: 0.0}",
])


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.paths = []

    def tearDown(self):
        for p in self.paths:
            try:
                os.remove(p)
            except OSError:
                pass

    def _write(self, text):
        fd, path = tempfile.mkstemp(suffix=".yaml", dir=os.path.dirname(__file__))
        os.close(fd)
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
        self.paths.append(path)
        return path


class TestByCoord3D(_Tmp):

    def test_x_min_and_z_max_on_hex(self):
        path = self._write(_hex_yaml("\n".join([
            "boundary_conditions_by_coord:",
            "  - {loc: x_min, ux: 0.0, uy: 0.0, uz: 0.0}",
            "point_loads_by_coord:",
            "  - {loc: z_max, uz: 1.0}",
            "solver:",
            "  type: LinearSolver",
        ])))
        res = solidum.run_yaml(path)
        self.assertGreater(float(np.max(res.U)), 0.0)
        dom = res.U.size
        self.assertEqual(dom, 24)
        # Los 4 nodos de z = 1 reciben la carga.
        self.assertAlmostEqual(float(res.F_applied.sum()), 4.0, places=12)

    def test_dict_form_keeps_loc_for_loads(self):
        path = self._write(_hex_yaml("\n".join([
            "boundary_conditions_by_coord:",
            "  x_min: {ux: 0.0, uy: 0.0, uz: 0.0}",
            "point_loads_by_coord:",
            "  x_max: {ux: 2.0}",
            "solver:",
            "  type: LinearSolver",
        ])))
        res = solidum.run_yaml(path)
        self.assertAlmostEqual(float(res.F_applied.sum()), 8.0, places=12)

    def test_selector_without_match_raises(self):
        path = self._write(_hex_yaml("\n".join([
            "boundary_conditions_by_coord:",
            "  - {coord: z, val: 5.0, ux: 0.0}",
        ])))
        with self.assertRaises(ValueError):
            YamlParser(path).parse()


class TestStrictPointLoads(_Tmp):

    def test_unknown_dof_raises(self):
        path = self._write(BAR_2D + "\npoint_loads:\n  - {node_id: 2, fy: -1.0}\n")
        parser = YamlParser(path)
        parser.parse(); parser.domain.generate_equation_numbers()
        with self.assertRaises(ValueError):
            parser.get_external_forces()

    def test_unknown_node_is_rejected_in_validation(self):
        path = self._write(BAR_2D + "\npoint_loads:\n  - {node_id: 99, ux: -1.0}\n")
        with self.assertRaises(YamlValidationError):
            YamlParser(path).parse()

    def test_group_without_mesh_raises(self):
        path = self._write(BAR_2D + "\npoint_loads_by_group:\n  - {group_name: borde, ux: 1.0}\n")
        parser = YamlParser(path)
        parser.parse(); parser.domain.generate_equation_numbers()
        with self.assertRaises(ValueError):
            parser.get_external_forces()


class TestTransientLoadsFromYaml(_Tmp):

    def test_newmark_receives_step_load(self):
        path = self._write(BAR_2D + "\n".join([
            "",
            "point_loads:",
            "  - {node_id: 2, ux: 10.0}",
            "solver:",
            "  type: NewmarkSolver",
            "  t_end: 1.2566",
            "  dt: 0.006283",
            "  lumping: lumped",
            "",
        ]))
        res = solidum.run_yaml(path)
        dof = 2  # ux del nodo 2 (numeración por nodos: ux1, uy1, ux2, uy2)
        u = res.u_history[dof, :]
        # Respuesta al escalón: u = (F/K)(1 − cos ωt), K = 25, F = 10 → máx 0.8.
        self.assertGreater(float(np.max(u)), 0.75)
        self.assertLess(float(np.max(u)), 0.85)

    def test_explicit_f_func_is_not_overridden(self):
        path = self._write(BAR_2D + "\npoint_loads:\n  - {node_id: 2, ux: 10.0}\nsolver:\n  type: NewmarkSolver\n  t_end: 0.1\n  dt: 0.01\n")
        from solidum.math.assembly import Assembler
        parser = YamlParser(path)
        dom = parser.parse(); dom.generate_equation_numbers()
        solver = parser.get_solver(Assembler(dom))
        self.assertIsNotNone(solver.F_func)
        np.testing.assert_allclose(solver.F_func(0.37)[2], 10.0)


class TestValidation(_Tmp):

    def test_unknown_material_and_solver_kwargs(self):
        path = self._write(BAR_2D.replace("E: 25.0, density: 2.0", "E: 25.0, densidad: 2.0")
                           + "\nsolver:\n  type: NonlinearSolver\n  max_iters: 5\n")
        with self.assertRaises(YamlValidationError) as cm:
            YamlParser(path).parse()
        msg = str(cm.exception)
        self.assertIn("densidad", msg)
        self.assertIn("max_iters", msg)

    def test_mesh_with_inline_nodes_is_rejected(self):
        path = self._write(BAR_2D + "\nmesh: placa.msh\nmesh_material: 1\n")
        with self.assertRaises(YamlValidationError):
            YamlParser(path).parse()

    def test_mesh_material_missing_is_rejected(self):
        msh = os.path.join(os.path.dirname(__file__), "..", "examples", "placa.msh")
        if not os.path.isfile(msh):
            self.skipTest("placa.msh no disponible")
        path = self._write("\n".join([
            "materials:",
            "  - {id: 1, type: Elastic2D, E: 1.0, nu: 0.3}",
            f"mesh: {os.path.relpath(msh, os.path.dirname(__file__))}",
            "mesh_material: 7",
            "",
        ]))
        with self.assertRaises(YamlValidationError):
            YamlParser(path).parse()

    def test_convergence_on_linear_solver_warns_instead_of_silently_dropping(self):
        import logging
        path = self._write(BAR_2D + "\nsolver:\n  type: LinearSolver\n  convergence: {rtol_force: 1e-8}\n")
        from solidum.math.assembly import Assembler
        parser = YamlParser(path)
        dom = parser.parse(); dom.generate_equation_numbers()
        with self.assertLogs("solidum.parsers.yaml", level=logging.WARNING):
            parser.get_solver(Assembler(dom))

    def test_global_safe_loader_untouched(self):
        # El resolver de floats 1.2 vive en un Loader propio; el global de
        # PyYAML sigue interpretando '1e5' como cadena (YAML 1.1).
        self.assertIsInstance(yaml.safe_load("v: 1e5")["v"], str)
        path = self._write(BAR_2D.replace("E: 25.0", "E: 25e0"))
        parser = YamlParser(path); parser.parse()
        self.assertEqual(parser.materials[1].E, 25.0)


if __name__ == "__main__":
    unittest.main()
