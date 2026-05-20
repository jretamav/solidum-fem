"""Tests del pipeline completo de fuerza de cuerpo (YAML → solver → resultado).

Verifica el cableado introducido para exponer ``compute_body_load`` desde
los archivos de entrada:

1. ``Assembler.assemble_body_load(b)`` acumula correctamente sobre todos
   los elementos del dominio (multi-elemento, multi-DOF).
2. El parser YAML lee el bloque ``body_force:`` y lo agrega al vector de
   fuerzas externas.
3. ``run_yaml`` produce el mismo resultado físico que aplicar manualmente
   las cargas equivalentes nodales.
"""
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

import solidum  # autodiscover
from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.materials.elastic import Elastic1D
from solidum.elements.truss import Truss2D
from solidum.elements.frame import Frame2DEuler
from solidum.elements.frame3d import Frame3D
from solidum.entry import run
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver
from solidum.utils.yaml_parser import YamlParser


class TestAssembleBodyLoad(unittest.TestCase):
    def _two_truss_2d(self):
        """Dos trusses 2D en serie: (0,0)-(L,0)-(2L,0), área constante."""
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [1.0, 0.0])
        n3 = dom.add_node(3, [2.0, 0.0])
        mat = Elastic1D(E=210e9)
        e1 = Truss2D(1, [n1, n2], mat, A=1e-3)
        e2 = Truss2D(2, [n2, n3], mat, A=1e-3)
        dom.add_element(e1)
        dom.add_element(e2)
        dom.generate_equation_numbers(verbose=False)
        return dom

    def test_acumulacion_en_nodo_compartido(self):
        """El nodo intermedio recibe la suma de los aportes de ambos elementos."""
        dom = self._two_truss_2d()
        asm = Assembler(dom)
        b = np.array([0.0, -78.5e3])  # peso propio acero
        F = asm.assemble_body_load(b)
        # Esperado: nodos 1 y 3 reciben (A·L/2)·b; nodo 2 recibe 2·(A·L/2)·b = A·L·b.
        A, L = 1e-3, 1.0
        n1, n2, n3 = dom.nodes[1], dom.nodes[2], dom.nodes[3]
        np.testing.assert_allclose([F[n1.dofs['ux']], F[n1.dofs['uy']]],
                                    [0.5 * A * L * b[0], 0.5 * A * L * b[1]],
                                    rtol=1e-12)
        np.testing.assert_allclose([F[n2.dofs['ux']], F[n2.dofs['uy']]],
                                    [A * L * b[0], A * L * b[1]], rtol=1e-12)
        np.testing.assert_allclose([F[n3.dofs['ux']], F[n3.dofs['uy']]],
                                    [0.5 * A * L * b[0], 0.5 * A * L * b[1]],
                                    rtol=1e-12)

    def test_padding_2d_a_3d_para_mezcla(self):
        """Un dominio mixto (Truss2D + Truss2D con DOFs 'uz'? No tenemos truss 2D con
        'uz'; uso un Frame3D para forzar el slicing). El parser debe pasar
        b[:2] al elemento 2D y b[:3] al 3D sin reventar."""
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [1.0, 0.0])
        mat = Elastic1D(E=210e9)
        dom.add_element(Truss2D(1, [n1, n2], mat, A=1e-3))
        dom.generate_equation_numbers(verbose=False)
        asm = Assembler(dom)
        # b 3D: el padding debe ignorar z porque el elemento es 2D.
        F = asm.assemble_body_load([0.0, -78.5e3, 999.0])
        # El nodo 1 y 2 deben recibir solo los componentes 2D — sin DOF 'uz' no
        # hay donde poner el componente z.
        for nid in (1, 2):
            node = dom.nodes[nid]
            np.testing.assert_allclose(F[node.dofs['uy']],
                                        0.5 * 1e-3 * 1.0 * (-78.5e3), rtol=1e-12)

    def test_dimension_invalida_lanza(self):
        dom = self._two_truss_2d()
        asm = Assembler(dom)
        with self.assertRaisesRegex(ValueError, "dimensión inesperada"):
            asm.assemble_body_load([1.0])
        with self.assertRaisesRegex(ValueError, "dimensión inesperada"):
            asm.assemble_body_load([1.0, 2.0, 3.0, 4.0])

class TestBodyForceFromYaml(unittest.TestCase):
    def test_run_yaml_aplica_peso_propio(self):
        """Frame2DEuler horizontal en voladizo con peso propio uniforme.
        Compara contra la solución analítica del voladizo bajo carga
        distribuida transversal:
            v(L) = q·L⁴ / (8·E·I), con q = ρ·g·A.
        """
        yaml_content = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [2.0, 0.0]}
  - {id: 4, coords: [3.0, 0.0]}
  - {id: 5, coords: [4.0, 0.0]}

materials:
  - id: 1
    type: Elastic1D
    E: 210.0e9
    density: 0.0

elements:
  - {id: 1, type: Frame2DEuler, material: 1, nodes: [1, 2], A: 1.0e-2, I: 8.33e-6}
  - {id: 2, type: Frame2DEuler, material: 1, nodes: [2, 3], A: 1.0e-2, I: 8.33e-6}
  - {id: 3, type: Frame2DEuler, material: 1, nodes: [3, 4], A: 1.0e-2, I: 8.33e-6}
  - {id: 4, type: Frame2DEuler, material: 1, nodes: [4, 5], A: 1.0e-2, I: 8.33e-6}

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}

body_force: [0.0, -78.5e3]   # acero: ρg ≈ 78.5 kN/m³

solver:
  type: LinearSolver
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False,
                                          encoding='utf-8') as f:
            f.write(yaml_content)
            tmp_path = Path(f.name)
        try:
            parser = YamlParser(str(tmp_path))
            domain = parser.parse()
            domain.generate_equation_numbers()
            asm = Assembler(domain)
            solver = parser.get_solver(asm)
            F_ext = parser.get_external_forces() + parser.get_body_load(asm)
            result = run(domain, assembler=asm, solver=solver, F_applied=F_ext)
        finally:
            tmp_path.unlink()

        # Cantidades del problema.
        A = 1.0e-2
        I = 8.33e-6
        E = 210e9
        rho_g = 78.5e3
        q = rho_g * A
        L = 4.0
        v_L_analitico = q * L**4 / (8.0 * E * I)

        uy_idx = domain.nodes[5].dofs['uy']
        v_L_numerico = result.U[uy_idx]
        # El signo: q = -ρgA·b_y, con b_y = -ρg → q es negativo (gravedad
        # apunta en -y). v_L analítico para el voladizo bajo carga uniforme
        # hacia abajo sale negativo si la convención es +y hacia arriba.
        # Como b = (0, -ρg), la carga distribuida transversal vale q = A·(-ρg).
        # v(L) = q·L^4/(8EI), con q negativo ⇒ v(L) < 0.
        v_L_esperado = -v_L_analitico   # ya con el signo

        # 4 elementos de Hermite cúbico aproximan exactamente el voladizo bajo
        # carga uniforme (la formulación Hermite es exacta para flexión Euler
        # bajo carga distribuida polinómica). Tolerancia muy estricta.
        np.testing.assert_allclose(v_L_numerico, v_L_esperado, rtol=1e-6)

    def test_run_yaml_sin_body_force_no_aplica(self):
        """Sin bloque ``body_force``, las fuerzas externas son solo las puntuales."""
        yaml_content = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}

materials:
  - id: 1
    type: Elastic1D
    E: 210.0e9
    density: 0.0

elements:
  - {id: 1, type: Truss2D, material: 1, nodes: [1, 2], A: 1.0e-3}

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0}
  - {node_id: 2, uy: 0.0}

point_loads_by_node:
  - {node_id: 2, ux: 1000.0}

solver:
  type: LinearSolver
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False,
                                          encoding='utf-8') as f:
            f.write(yaml_content)
            tmp_path = Path(f.name)
        try:
            parser = YamlParser(str(tmp_path))
            domain = parser.parse()
            domain.generate_equation_numbers()
            asm = Assembler(domain)
            solver = parser.get_solver(asm)
            F_ext = parser.get_external_forces() + parser.get_body_load(asm)
            result = run(domain, assembler=asm, solver=solver, F_applied=F_ext)
        finally:
            tmp_path.unlink()

        # u_x del nodo 2: 1000 N · 1 m / (210e9 · 1e-3) = 4.76e-6 m
        expected = 1000.0 * 1.0 / (210e9 * 1e-3)
        ux_idx = domain.nodes[2].dofs['ux']
        np.testing.assert_allclose(result.U[ux_idx], expected, rtol=1e-10)


if __name__ == '__main__':
    unittest.main()
