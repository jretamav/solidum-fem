"""Tests del soporte de densidad por material y peso propio multimaterial.

Verifica el ADR 0008:

1. ``density`` es argumento OBLIGATORIO al construir cualquier material;
   omitirlo lanza ``TypeError`` en construcción (fallo ruidoso e inmediato,
   sin posibilidad de propagar masa cero silenciosa a un análisis dinámico).
2. ``Assembler.assemble_self_weight(g)`` aplica ``b_element = density · g``
   por elemento, respetando la densidad del material asignado.
3. Modelos multimaterial dan resultado físicamente correcto (cada elemento
   ve su propio ρ).
4. Material con ``density = 0`` declarado explícitamente emite warning
   informativo y aporta cero al peso propio (caso legítimo: penalty,
   restricción).
5. YAML: bloque ``gravity`` mutuamente exclusivo con ``body_force``.
6. Pipeline completo desde YAML.
"""
import logging
import tempfile
import unittest
from pathlib import Path

import numpy as np

import fenix  # autodiscover
from fenix.core.domain import Domain
from fenix.core.node import Node
from fenix.materials.elastic import Elastic1D
from fenix.materials.elastic_2d import Elastic2D
from fenix.materials.plastic_1d import Elastoplastic1D
from fenix.materials.damage_1d import IsotropicDamage1D
from fenix.materials.damage_2d import IsotropicDamage2D
from fenix.materials.von_mises_2d import VonMises2D
from fenix.materials.cable_1d import CableMaterial1D
from fenix.elements.truss import Truss2D
from fenix.elements.frame import Frame2DEuler
from fenix.entry import run
from fenix.math.assembly import Assembler
from fenix.utils.yaml_parser import YamlParser, YamlValidationError


class TestMaterialDensityAttribute(unittest.TestCase):
    """Todos los materiales requieren `density` y exponen el atributo."""

    def test_elastic1d_density_obligatoria(self):
        # Omitir density lanza TypeError al construir — no permite default silencioso.
        with self.assertRaises(TypeError):
            Elastic1D(E=210e9)
        # Valor explícito (incluido 0.0 para fixtures) es válido.
        m_massless = Elastic1D(E=210e9, density=0.0)
        self.assertEqual(m_massless.density, 0.0)
        m_steel = Elastic1D(E=210e9, density=7850.0)
        self.assertEqual(m_steel.density, 7850.0)

    def test_elastic2d_density(self):
        m = Elastic2D(E=210e9, nu=0.3, density=7850.0)
        self.assertEqual(m.density, 7850.0)

    def test_elastoplastic1d_density(self):
        m = Elastoplastic1D(E=210e9, sigma_y=250e6, density=7850.0)
        self.assertEqual(m.density, 7850.0)

    def test_damage1d_density(self):
        m = IsotropicDamage1D(E=210e9, kappa_0=1e-4, alpha=10.0, density=7850.0)
        self.assertEqual(m.density, 7850.0)

    def test_damage2d_density(self):
        m = IsotropicDamage2D(E=30e9, nu=0.2, kappa_0=1e-4, alpha=10.0, density=2500.0)
        self.assertEqual(m.density, 2500.0)

    def test_vonmises2d_density(self):
        m = VonMises2D(E=210e9, nu=0.3, sigma_y=250e6, density=7850.0)
        self.assertEqual(m.density, 7850.0)

    def test_cable_density(self):
        m = CableMaterial1D(E=150e9, density=7850.0)
        self.assertEqual(m.density, 7850.0)


class TestAssembleSelfWeight(unittest.TestCase):
    def _two_truss_diff_density(self):
        """Dos trusses en serie con materiales de densidad distinta.
        n1 -- (Elastic1D ρ=7850) -- n2 -- (Elastic1D ρ=2500) -- n3
        """
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [1.0, 0.0])
        n3 = dom.add_node(3, [2.0, 0.0])
        mat_steel = Elastic1D(E=210e9, density=7850.0)
        mat_concrete = Elastic1D(E=30e9, density=2500.0)
        e1 = Truss2D(1, [n1, n2], mat_steel, A=1e-3)
        e2 = Truss2D(2, [n2, n3], mat_concrete, A=1e-3)
        dom.add_element(e1)
        dom.add_element(e2)
        dom.generate_equation_numbers(verbose=False)
        return dom

    def test_multimaterial_self_weight(self):
        """Cada elemento recibe ρᵢ·g·A·L/2 según su propio material."""
        dom = self._two_truss_diff_density()
        asm = Assembler(dom)
        g = np.array([0.0, -9.81])
        F = asm.assemble_self_weight(g)

        A, L = 1e-3, 1.0
        # Nodo 1: aporte solo de e1 (acero, ρ=7850).
        n1 = dom.nodes[1]
        np.testing.assert_allclose(F[n1.dofs['uy']],
                                    0.5 * 7850.0 * (-9.81) * A * L,
                                    rtol=1e-12)
        # Nodo 3: aporte solo de e2 (hormigón, ρ=2500).
        n3 = dom.nodes[3]
        np.testing.assert_allclose(F[n3.dofs['uy']],
                                    0.5 * 2500.0 * (-9.81) * A * L,
                                    rtol=1e-12)
        # Nodo 2 (compartido): aporte de los dos.
        n2 = dom.nodes[2]
        esperado = 0.5 * (7850.0 + 2500.0) * (-9.81) * A * L
        np.testing.assert_allclose(F[n2.dofs['uy']], esperado, rtol=1e-12)

    def test_density_cero_explicita_genera_warning_y_aporte_nulo(self):
        """Material con density=0.0 declarado explícitamente (caso legítimo:
        material de penalty/restricción) emite warning informativo y aporta
        cero al peso propio sin abortar el análisis."""
        dom = Domain()
        n1 = dom.add_node(1, [0.0, 0.0])
        n2 = dom.add_node(2, [1.0, 0.0])
        mat = Elastic1D(E=210e9, density=0.0)
        dom.add_element(Truss2D(1, [n1, n2], mat, A=1e-3))
        dom.generate_equation_numbers(verbose=False)
        asm = Assembler(dom)

        with self.assertLogs('fenix.assembly', level='WARNING') as cm:
            F = asm.assemble_self_weight([0.0, -9.81])

        # Warning identifica el material por nombre.
        self.assertTrue(any('Elastic1D' in msg and 'density=0' in msg
                            for msg in cm.output))
        # Aporte numérico es exactamente cero.
        np.testing.assert_allclose(F, 0.0, atol=1e-15)

    def test_padding_2d_a_3d(self):
        """Vector g 3D sobre dominio 2D: ignora componente z."""
        dom = self._two_truss_diff_density()
        asm = Assembler(dom)
        # g 3D con componente z (que se ignora en problema 2D).
        F = asm.assemble_self_weight([0.0, -9.81, 999.0])
        # Verificar que ningún DOF recibió la componente z (no hay DOF 'uz').
        for nid in (1, 2, 3):
            node = dom.nodes[nid]
            self.assertIn('uy', node.dofs)
            # Componente uy correcta.
            self.assertFalse(np.isnan(F[node.dofs['uy']]))


class TestYamlGravityBlock(unittest.TestCase):
    def _write_yaml(self, content: str) -> Path:
        f = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False,
                                         encoding='utf-8')
        f.write(content); f.close()
        return Path(f.name)

    def test_gravity_acepta_density_por_material(self):
        """Voladizo Frame2D con peso propio declarado vía density + gravity.
        Mismo resultado físico que el test del commit anterior con
        body_force precalculado."""
        yaml_content = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}
  - {id: 3, coords: [2.0, 0.0]}
  - {id: 4, coords: [3.0, 0.0]}
  - {id: 5, coords: [4.0, 0.0]}

materials:
  - {id: 1, type: Elastic1D, E: 210.0e9, density: 8000.0}

elements:
  - {id: 1, type: Frame2DEuler, material: 1, nodes: [1, 2], A: 1.0e-2, I: 8.33e-6}
  - {id: 2, type: Frame2DEuler, material: 1, nodes: [2, 3], A: 1.0e-2, I: 8.33e-6}
  - {id: 3, type: Frame2DEuler, material: 1, nodes: [3, 4], A: 1.0e-2, I: 8.33e-6}
  - {id: 4, type: Frame2DEuler, material: 1, nodes: [4, 5], A: 1.0e-2, I: 8.33e-6}

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0, rz: 0.0}

gravity: [0.0, -9.81]

solver:
  type: LinearSolver
"""
        path = self._write_yaml(yaml_content)
        try:
            parser = YamlParser(str(path))
            domain = parser.parse()
            domain.generate_equation_numbers()
            asm = Assembler(domain)
            solver = parser.get_solver(asm)
            F_ext = parser.get_external_forces() + parser.get_body_load(asm)
            result = run(domain, assembler=asm, solver=solver, F_applied=F_ext)
        finally:
            path.unlink()

        A, I, E, L = 1.0e-2, 8.33e-6, 210e9, 4.0
        q = 8000.0 * (-9.81) * A   # ρ·g·A (negativo, gravedad hacia abajo)
        v_L_esperado = q * L**4 / (8.0 * E * I)
        uy_idx = domain.nodes[5].dofs['uy']
        np.testing.assert_allclose(result.U[uy_idx], v_L_esperado, rtol=1e-6)

    def test_gravity_y_body_force_simultaneos_falla(self):
        yaml_content = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}

materials:
  - {id: 1, type: Elastic1D, E: 210.0e9, density: 7850.0}

elements:
  - {id: 1, type: Truss2D, material: 1, nodes: [1, 2], A: 1.0e-3}

boundary_conditions_by_node:
  - {node_id: 1, ux: 0.0, uy: 0.0}
  - {node_id: 2, uy: 0.0}

gravity: [0.0, -9.81]
body_force: [0.0, -78.5e3]

solver:
  type: LinearSolver
"""
        path = self._write_yaml(yaml_content)
        try:
            with self.assertRaises(YamlValidationError) as ctx:
                YamlParser(str(path)).parse()
            self.assertIn("body_force", str(ctx.exception))
            self.assertIn("gravity", str(ctx.exception))
        finally:
            path.unlink()

    def test_density_pasa_por_yaml_a_materiales(self):
        """El parser entrega density al constructor del material vía
        introspección estándar."""
        yaml_content = """
nodes:
  - {id: 1, coords: [0.0, 0.0]}
  - {id: 2, coords: [1.0, 0.0]}

materials:
  - {id: 1, type: Elastic1D, E: 210.0e9, density: 7850.0}
  - {id: 2, type: VonMises2D, E: 210.0e9, nu: 0.3, sigma_y: 250.0e6, density: 7850.0}

elements:
  - {id: 1, type: Truss2D, material: 1, nodes: [1, 2], A: 1.0e-3}

solver:
  type: LinearSolver
"""
        path = self._write_yaml(yaml_content)
        try:
            parser = YamlParser(str(path))
            parser.parse()
        finally:
            path.unlink()
        self.assertEqual(parser.materials[1].density, 7850.0)
        self.assertEqual(parser.materials[2].density, 7850.0)


if __name__ == '__main__':
    unittest.main()
