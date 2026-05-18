"""Validación experimental de la aportación del Cap. 6 de Retama (2010).

Verifica numéricamente que la longitud efectiva ``l_d = (A_e/h)·cos(θ − α)``
del Cap. 6 cura el *stress locking* que aparecería con la alternativa ingenua
``l_d_ingenuo = A_e/h``. La tesis no incluyó este experimento; lo aporta
Fenix como fase 3 del ADR 0010 (decisión arquitectural cerrada — la fórmula
correcta es la única implementada en producción; los tests prueban la
**alternativa rechazada** mediante monkey-patch para documentar el motivo).

Tres tests:

1. **Aritmético oblicuo** — para una geometría con ``Γ_d`` no paralela al
   lado opuesto al nodo solitario, ``l_d_correcto`` y ``l_d_ingenuo``
   difieren por un factor exacto ``|cos(θ − α)|``.

2. **Aritmético paralelo** — cuando ``Γ_d`` es paralela al lado opuesto,
   ``cos(θ − α) = ±1`` y ambas fórmulas coinciden.

3. **Stress locking físico** — con ``Γ_d`` oblicua, dada una misma carga
   aplicada, la apertura ``[[u_n]]`` resultante del Newton local es **mayor**
   con ``l_d_correcto`` que con ``l_d_ingenuo``: la versión ingenua "infla"
   la rigidez cohesiva (``l_d·T_coh``) y subestima la apertura — comportamiento
   más rígido de lo físico, *stress locking*.
"""
from __future__ import annotations

import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from fenix.cohesive_materials.damage_isotropic import CohesiveDamageIsotropic
from fenix.core.node import Node
from fenix.elements.solid_2d.embedded_cst import CST_Embedded2D, _compute_ld
from fenix.materials.elastic_2d import Elastic2D


# --- Helpers ---------------------------------------------------------------

def _triangle_nodes(coords):
    return [Node(i, list(coords[i])) for i in range(len(coords))]


def _make_elem(coords, *, sigma_t0=2.5e6, softening='linear'):
    bulk = Elastic2D(E=30.0e9, nu=0.2, hypothesis='plane_strain')
    cohesive = CohesiveDamageIsotropic(
        sigma_t0=sigma_t0, G_f=100.0, K_e=1.0e13, softening=softening,
    )
    return CST_Embedded2D(1, _triangle_nodes(coords), bulk, cohesive)


def _assign_local_dofs(element):
    idx = 0
    for node in element.nodes:
        for dof_name in element.DOF_NAMES:
            node.dofs[dof_name] = idx
            idx += 1


# =========================================================================
# 1. Aritmético — Γ_d oblicua
# =========================================================================


class TestArithmeticObliqueCase(unittest.TestCase):
    """Geometría de prueba: triángulo rectángulo (0,0), (1,0), (0,1).

    Activado con ``n = (sin 60°, −cos 60°) = (√3/2, −1/2)`` ⇒ tangente
    ``s = (1/2, √3/2)`` ⇒ ``θ = 60°``. El nodo solitario resulta el 1.
    Lado opuesto al nodo 1: lado 0–2 (vertical, ``α = 90°``).

    Resultado analítico: ``cos(θ − α) = cos(−30°) = √3/2``. Como ``A_e = 1/2``
    y ``h = 1`` (distancia de (1,0) al eje y), ``l_d_correcto = √3/4`` y
    ``l_d_ingenuo = 1/2``, con razón ``√3/2 ≈ 0.866``.
    """

    def test_l_d_correcto_vs_ingenuo(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        tangent = np.array([0.5, math.sqrt(3) / 2.0])  # θ = 60°
        l_d_correcto = _compute_ld(coords, tangent, solitary=1)

        A_e = 0.5
        h = 1.0                                   # distancia de (1,0) al eje y
        l_d_ingenuo = A_e / h

        expected_correcto = (A_e / h) * (math.sqrt(3) / 2.0)
        self.assertAlmostEqual(l_d_correcto, expected_correcto, places=14)
        self.assertAlmostEqual(l_d_ingenuo, 0.5, places=14)
        # Razón exacta cos(θ - α) = √3/2
        self.assertAlmostEqual(
            l_d_correcto / l_d_ingenuo, math.sqrt(3) / 2.0, places=14,
        )

    def test_activation_path_produces_same_l_d(self):
        """La activación vía ``_activate`` (con la convención de orientar n
        hacia el solitario) produce el mismo l_d que la función auxiliar."""
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        emb = _make_elem(coords)
        # Forzar n input que provoque solitario = 1 sin inversión
        n_in = np.array([math.sqrt(3) / 2.0, -0.5])  # apunta hacia nodo 1
        emb._activate(n_in, coords)
        ds = emb.discontinuity_state
        self.assertEqual(ds.solitary_node, 1)
        self.assertAlmostEqual(ds.l_d, (0.5) * (math.sqrt(3) / 2.0), places=14)


# =========================================================================
# 2. Aritmético — Γ_d paralela al lado opuesto
# =========================================================================


class TestArithmeticParallelCase(unittest.TestCase):
    """Γ_d horizontal sobre el mismo triángulo: tangente ``s = (−1, 0)`` ⇒
    ``θ = 180°`` (post-inversión por la convención de n hacia el solitario).
    Lado opuesto al nodo 2 = lado 0–1, ``α = 0``. ``cos(θ − α) = −1``,
    ``|cos| = 1`` ⇒ ``l_d_correcto = l_d_ingenuo = A_e/h``.
    """

    def test_l_d_coincide(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        emb = _make_elem(coords)
        emb._activate(np.array([0.0, 1.0]), coords)  # n hacia (0,1)
        ds = emb.discontinuity_state
        self.assertEqual(ds.solitary_node, 2)

        # Lado opuesto al nodo 2 = lado 0-1 (horizontal), α = 0.
        # Tangente post-inversión = (-1, 0), θ = 180°. |cos(180°)| = 1.
        A_e = 0.5
        h = 1.0     # distancia de (0,1) al eje x
        l_d_ingenuo = A_e / h
        self.assertAlmostEqual(ds.l_d, l_d_ingenuo, places=14)


# =========================================================================
# 3. Stress locking físico — l_d ingenuo subestima la apertura
# =========================================================================


class TestNaiveLdCausesStressLocking(unittest.TestCase):
    """Para una misma carga sobre el elemento, ``l_d_ingenuo > l_d_correcto``
    en orientación oblicua. La rigidez cohesiva proyectada al sistema local
    del salto, ``l_d · T_coh``, es **mayor** con la versión ingenua. En
    equilibrio (R^{[[u]]} = 0), la apertura ``[[u_n]]`` que produce el
    cohesivo para balancear el bulk es **menor** con la versión ingenua.
    Eso es el *stress locking* que la corrección del Cap. 6 evita.

    Se controla el experimento por *monkey-patch* de ``ds.l_d``: la
    implementación de producción siempre usa ``l_d_correcto`` (decisión
    cerrada en ADR 0010 §6). El test simula la alternativa rechazada para
    documentar su impacto cuantitativo.
    """

    def test_naive_ld_underestimates_jump(self):
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        # Softening exponencial: ω → 1 asintóticamente. Si fuera lineal, una
        # carga suficiente para activar lleva fácilmente a κ ≥ w_c (saturación
        # exacta, t = 0), donde l_d no aparece en R^{[[u]]} y ambas versiones
        # producirían la misma apertura — el locking sólo se manifiesta en
        # rama activa de softening con t ≠ 0.
        emb = _make_elem(coords, sigma_t0=1.0e6, softening='exponential')
        _assign_local_dofs(emb)

        # Activar con orientación oblicua (la del primer test aritmético).
        n_in = np.array([math.sqrt(3) / 2.0, -0.5])
        emb._activate(n_in, coords)
        ds = emb.discontinuity_state
        self.assertEqual(ds.solitary_node, 1)
        l_d_correcto = ds.l_d

        # Carga: u_x del nodo 1 = 1e-4 (provoca apertura post-activación
        # pero en rama de softening exponencial activa, sin saturar).
        u_e = np.array([0.0, 0.0, 1.0e-4, 0.0, 0.0, 0.0])

        # --- Caso 1: con l_d correcto (producción) ----------------------
        emb.compute_element_state(u_e)
        jump_correcto = np.copy(ds.jump_trial)

        # --- Caso 2: l_d ingenuo (monkey-patch) -------------------------
        # Reset del estado interno del elemento + cohesivo:
        ds.jump_trial = np.zeros(2)
        ds.jump_committed = np.zeros(2)
        ds.cohesive_state_trial = {}
        ds.cohesive_state_committed = {}
        emb.state.vars_trial = [None]
        emb.state.vars = [None]
        emb.state.stresses_trial = [np.zeros(3)]
        emb.state.stresses = [np.zeros(3)]
        # Inyectar l_d ingenuo (= A_e/h sin el factor cos(θ−α))
        ds.l_d = 0.5

        emb.compute_element_state(u_e)
        jump_ingenuo = np.copy(ds.jump_trial)

        # --- Verificaciones del locking ---------------------------------
        # La componente normal de la apertura debe ser positiva (la grieta abre).
        self.assertGreater(jump_correcto[0], 0.0)
        self.assertGreater(jump_ingenuo[0], 0.0)

        # La versión correcta produce MAYOR apertura normal:
        self.assertGreater(jump_correcto[0], jump_ingenuo[0])

        # Diferencia relativa significativa (no ruido numérico).
        rel_diff = (jump_correcto[0] - jump_ingenuo[0]) / jump_correcto[0]
        self.assertGreater(rel_diff, 0.05)

        # Comprobación adicional: l_d_correcto < l_d_ingenuo en oblicuo.
        self.assertLess(l_d_correcto, 0.5)


if __name__ == '__main__':
    unittest.main()
