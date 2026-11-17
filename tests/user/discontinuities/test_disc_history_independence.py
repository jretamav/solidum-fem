"""Los resultados de la discontinuidad embebida no dependen de lo que pide la salida (ADR 0040, punto 4).

Antes del ADR 0040 la historia reensamblaba en ``U`` tras el commit para leer
las reacciones, y eso modificaba el salto de prueba (``jump_trial``) del
Newton local de ``CST_Embedded2D``, del que arranca el paso siguiente: los
desplazamientos podían depender de si el YAML pedía una columna de reacción.
Ahora la historia toma la ``F_int`` del paso confirmado sin reensamblar.

Caso: la probeta de tracción uniaxial de dos ``CST_Embedded2D`` de
``test_disc_embedded_uniaxial_softening`` (Retama 2010), por el camino del
YAML, hasta la separación completa; ``U`` de todos los pasos, con y sin la
columna de reacción, idéntica bit a bit.

Es una salvaguarda, no una prueba de regresión: con este modelo el código
anterior también daba ``U`` idéntica (medido el 2026-11-17), porque el Newton
local converge exactamente desde el salto confirmado. La dependencia sólo
aparecía cuando el Newton local no convergía (deuda #12), que este caso no
ejercita. Protege que la historia no vuelva a reensamblar.
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

import solidum

MODELO = """\
nodes:
  - {{id: 1, coords: [0.0, 0.0]}}
  - {{id: 2, coords: [0.1, 0.0]}}
  - {{id: 3, coords: [0.1, 0.05]}}
  - {{id: 4, coords: [0.0, 0.05]}}
materials:
  - {{id: 1, type: Elastic2D, E: 30.0e9, nu: 0.2, hypothesis: plane_stress}}
user_modules: [discontinuities]
interface_materials:
  - {{id: 1, type: CohesiveDamageIsotropic, sigma_t0: 3.0e6, G_f: 100.0, K_e: 1.0e15,
     softening: linear}}
elements:
  - {{id: 1, type: CST_Embedded2D, nodes: [1, 2, 3], material: 1, cohesive_material: 1,
     thickness: 0.02}}
  - {{id: 2, type: CST_Embedded2D, nodes: [1, 3, 4], material: 1, cohesive_material: 1,
     thickness: 0.02}}
boundary_conditions_by_node:
  - {{node_id: 1, ux: 0.0, uy: 0.0}}
  - {{node_id: 2, ux: 8.0e-5, uy: 0.0}}
  - {{node_id: 3, ux: 8.0e-5}}
  - {{node_id: 4, ux: 0.0}}
solver: {{type: NonlinearSolver, num_steps: 75, adaptive: false}}
output:
  file_name: probeta
  log_level: ERROR
{historia}"""

HISTORIA = """\
  history:
    columns:
      - {name: P, reaction: {nodes: [2, 3], dof: ux}}
"""


def _desplazamientos(con_historia: bool) -> np.ndarray:
    d = Path(tempfile.mkdtemp())
    ruta = d / "probeta.yaml"
    ruta.write_text(MODELO.format(historia=HISTORIA if con_historia else ""), encoding="utf-8")
    pasos = []
    solidum.load_yaml(ruta).solve(step_callback=lambda s, U, lam: pasos.append(np.array(U)),
                                  write_output=True)
    return np.array(pasos)


class TestResultadosIndependientesDeLaSalida(unittest.TestCase):

    def test_U_identica_con_y_sin_columna_de_reaccion(self):
        con, sin = _desplazamientos(True), _desplazamientos(False)
        self.assertEqual(con.shape, sin.shape)
        self.assertGreater(np.max(np.abs(con)), 0.0)
        self.assertTrue(np.array_equal(con, sin))


if __name__ == "__main__":
    unittest.main()
