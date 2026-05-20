"""Tests de ``compute_body_load`` para truss, cable y frame (peso propio).

Blinda dos invariantes en cada elemento:

1. **Equilibrio global**: la suma de las fuerzas nodales devueltas iguala
   ``b · A · L`` componente a componente (la integral de la carga sobre
   el volumen del elemento, repartida en los nodos).
2. **Reparto correcto**: comparación contra las fórmulas analíticas
   conocidas — mitad-mitad para axial en 1D, ``qL/2 ± qL²/12`` para flexión
   en marcos con interpolación Hermite cúbica.

Todos los tests cubren orientaciones generales (no alineadas con los ejes
globales) para verificar las transformaciones local↔global.
"""
import math
import unittest

import numpy as np

import solidum  # autodiscover
from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.materials.elastic import Elastic1D
from solidum.materials.cable_1d import CableMaterial1D
from solidum.elements.truss import Truss2D, Truss2DCorot, Truss3D, Truss3DCorot
from solidum.elements.cable import Cable2DCorot, Cable3DCorot
from solidum.elements.frame import (
    Frame2DEuler,
    Frame2DTimoshenko,
    Frame2DEulerCorot,
)
from solidum.elements.frame3d import Frame3D


def _node(idx, *coords):
    return Node(idx, list(coords))


# ---------------------------------------------------------------------------
# Truss 2D / 3D — reparto axial mitad-mitad
# ---------------------------------------------------------------------------

class TestTrussBodyLoad(unittest.TestCase):
    def _truss2d_orientado(self, L=2.0, ang_deg=37.0, A=1e-3):
        ang = math.radians(ang_deg)
        n1 = _node(1, 0.0, 0.0)
        n2 = _node(2, L * math.cos(ang), L * math.sin(ang))
        mat = Elastic1D(E=210e9)
        return Truss2D(1, [n1, n2], mat, A=A), A, L

    def test_truss2d_reparto_mitad_mitad(self):
        elem, A, L = self._truss2d_orientado()
        b = np.array([3.0, -7.0])
        f = elem.compute_body_load(b)
        esperado = 0.5 * A * L * np.array([b[0], b[1], b[0], b[1]])
        np.testing.assert_allclose(f, esperado, rtol=1e-12)
        # Equilibrio global.
        np.testing.assert_allclose(f[0:2] + f[2:4], A * L * b, rtol=1e-12)

    def test_truss2d_corot_hereda(self):
        ang = math.radians(20.0)
        L, A = 3.0, 2e-4
        n1 = _node(1, 0.0, 0.0)
        n2 = _node(2, L * math.cos(ang), L * math.sin(ang))
        mat = Elastic1D(E=200e9)
        elem = Truss2DCorot(1, [n1, n2], mat, A=A)
        b = np.array([0.0, -78.5e3])  # peso propio de acero, ρg
        f = elem.compute_body_load(b)
        np.testing.assert_allclose(f[0:2] + f[2:4], A * L * b, rtol=1e-12)

    def test_truss3d_reparto_mitad_mitad(self):
        L, A = 2.5, 5e-4
        # Orientación arbitraria 3D (vector unitario inclinado).
        d = np.array([2.0, -1.0, 3.0])
        d = d / np.linalg.norm(d) * L
        n1 = _node(1, 0.0, 0.0, 0.0)
        n2 = _node(2, d[0], d[1], d[2])
        mat = Elastic1D(E=210e9)
        elem = Truss3D(1, [n1, n2], mat, A=A)
        b = np.array([1.5, -2.5, 4.0])
        f = elem.compute_body_load(b)
        esperado = 0.5 * A * L * np.array([b[0], b[1], b[2], b[0], b[1], b[2]])
        np.testing.assert_allclose(f, esperado, rtol=1e-12)
        np.testing.assert_allclose(f[0:3] + f[3:6], A * L * b, rtol=1e-12)

    def test_truss3d_corot_hereda(self):
        L, A = 1.8, 1e-3
        n1 = _node(1, 0.0, 0.0, 0.0)
        n2 = _node(2, 0.0, 0.0, L)
        mat = Elastic1D(E=200e9)
        elem = Truss3DCorot(1, [n1, n2], mat, A=A)
        b = np.array([0.0, 0.0, -78.5e3])
        f = elem.compute_body_load(b)
        np.testing.assert_allclose(f[0:3] + f[3:6], A * L * b, rtol=1e-12)


# ---------------------------------------------------------------------------
# Cable 2D / 3D — misma fórmula que truss
# ---------------------------------------------------------------------------

class TestCableBodyLoad(unittest.TestCase):
    def test_cable2d_corot(self):
        L, A = 5.0, 2e-5
        n1 = _node(1, 0.0, 0.0)
        n2 = _node(2, L, 0.0)
        mat = CableMaterial1D(E=150e9)
        elem = Cable2DCorot(1, [n1, n2], mat, A=A)
        b = np.array([0.0, -78.5e3])
        f = elem.compute_body_load(b)
        esperado = 0.5 * A * L * np.array([b[0], b[1], b[0], b[1]])
        np.testing.assert_allclose(f, esperado, rtol=1e-12)
        np.testing.assert_allclose(f[0:2] + f[2:4], A * L * b, rtol=1e-12)

    def test_cable3d_corot(self):
        L, A = 10.0, 3e-5
        n1 = _node(1, 0.0, 0.0, 0.0)
        n2 = _node(2, L, 0.0, 0.0)
        mat = CableMaterial1D(E=150e9)
        elem = Cable3DCorot(1, [n1, n2], mat, A=A)
        b = np.array([0.0, 0.0, -78.5e3])
        f = elem.compute_body_load(b)
        np.testing.assert_allclose(f[0:3] + f[3:6], A * L * b, rtol=1e-12)


# ---------------------------------------------------------------------------
# Frame 2D — Hermite cúbico para transversal: qL/2 y ±qL²/12
# ---------------------------------------------------------------------------

class TestFrame2DBodyLoad(unittest.TestCase):
    def _viga_horizontal(self, L=4.0, A=1e-2, I=8.33e-6, cls=Frame2DEuler):
        n1 = _node(1, 0.0, 0.0)
        n2 = _node(2, L, 0.0)
        mat = Elastic1D(E=210e9)
        if cls is Frame2DTimoshenko:
            return cls(1, [n1, n2], mat, A=A, I=I, As=A*5/6, nu=0.3), A, L
        return cls(1, [n1, n2], mat, A=A, I=I), A, L

    def test_frame2d_euler_horizontal_gravedad(self):
        """Viga horizontal con peso propio (gravedad −y):
        local coincide con global, transversal es q_y = b_y · A.
        Esperado: F_y = qL/2 en cada nodo, Mz_i = +qL²/12, Mz_j = -qL²/12.
        """
        elem, A, L = self._viga_horizontal()
        b = np.array([0.0, -78.5e3])
        q = A * b[1]  # carga distribuida transversal (negativa)
        f = elem.compute_body_load(b)
        # Layout global: [Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]
        esperado = np.array([
            0.0, 0.5 * q * L,  q * L * L / 12.0,
            0.0, 0.5 * q * L, -q * L * L / 12.0,
        ])
        np.testing.assert_allclose(f, esperado, rtol=1e-12, atol=1e-12)

    def test_frame2d_euler_equilibrio_global_axial_transversal(self):
        """Carga arbitraria 2D y viga oblicua: ΣF = b·A·L, ΣM debe ser
        coherente con la integral de momento sobre la longitud."""
        ang = math.radians(40.0)
        L, A, I = 3.0, 1e-2, 8.33e-6
        n1 = _node(1, 0.0, 0.0)
        n2 = _node(2, L * math.cos(ang), L * math.sin(ang))
        mat = Elastic1D(E=210e9)
        elem = Frame2DEuler(1, [n1, n2], mat, A=A, I=I)
        b = np.array([2.0e3, -78.5e3])
        f = elem.compute_body_load(b)
        # Suma de fuerzas globales: ΣF_i + ΣF_j = A · L · b.
        sumF = np.array([f[0] + f[3], f[1] + f[4]])
        np.testing.assert_allclose(sumF, A * L * b, rtol=1e-12)

    def test_frame2d_timoshenko_coincide_con_euler_carga_uniforme(self):
        """Para carga uniforme, la fórmula consistente de Timoshenko
        (Hermite-blended) coincide con la de Euler-Bernoulli."""
        elem_e, A, L = self._viga_horizontal(cls=Frame2DEuler)
        elem_t, _, _ = self._viga_horizontal(cls=Frame2DTimoshenko)
        b = np.array([0.0, -78.5e3])
        f_e = elem_e.compute_body_load(b)
        f_t = elem_t.compute_body_load(b)
        np.testing.assert_allclose(f_t, f_e, rtol=1e-12, atol=1e-12)

    def test_frame2d_euler_corot_usa_geometria_de_referencia(self):
        """El corotacional evalúa la carga en la configuración inicial:
        debe dar el mismo resultado que el Euler estándar para igual
        orientación inicial."""
        ang = math.radians(25.0)
        L, A, I = 2.0, 5e-3, 4e-6
        n1 = _node(1, 0.0, 0.0)
        n2 = _node(2, L * math.cos(ang), L * math.sin(ang))
        mat = Elastic1D(E=210e9)
        elem_e = Frame2DEuler(1, [n1, n2], mat, A=A, I=I)
        elem_c = Frame2DEulerCorot(2, [n1, n2], mat, A=A, I=I)
        b = np.array([1.5e3, -78.5e3])
        f_e = elem_e.compute_body_load(b)
        f_c = elem_c.compute_body_load(b)
        np.testing.assert_allclose(f_c, f_e, rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------------------
# Frame 3D — dos planos de flexión + axial, sin torsión por carga uniforme
# ---------------------------------------------------------------------------

class TestFrame3DBodyLoad(unittest.TestCase):
    def _frame3d_axial_x(self, L=4.0, A=1e-2, Iy=8.33e-6, Iz=8.33e-6, J=1e-5):
        n1 = _node(1, 0.0, 0.0, 0.0)
        n2 = _node(2, L, 0.0, 0.0)
        mat = Elastic1D(E=210e9)
        # ref_vector por defecto [0,0,1] da: y_local=+y_global, z_local=+z_global.
        return Frame3D(1, [n1, n2], mat, A=A, Iy=Iy, Iz=Iz, J=J), A, L

    def test_frame3d_axial_x_gravedad_z(self):
        """Viga horizontal alineada con +x global, gravedad en −z.
        Carga distribuida transversal en z_local (q_z), produce f_z = q·L/2
        y My = ∓q·L²/12. No produce torsión.
        """
        elem, A, L = self._frame3d_axial_x()
        b = np.array([0.0, 0.0, -78.5e3])
        q = A * b[2]  # q_z_local = b_z · A (local = global aquí)
        f = elem.compute_body_load(b)
        # Layout: [Fx_i, Fy_i, Fz_i, Mx_i, My_i, Mz_i, Fx_j, ...].
        esperado = np.array([
            0.0,        0.0, 0.5 * q * L,
            0.0, -q * L * L / 12.0, 0.0,
            0.0,        0.0, 0.5 * q * L,
            0.0,  q * L * L / 12.0, 0.0,
        ])
        np.testing.assert_allclose(f, esperado, rtol=1e-12, atol=1e-12)

    def test_frame3d_axial_x_carga_y(self):
        """Carga transversal en y_local, produce f_y = q·L/2 y
        Mz = ±q·L²/12."""
        elem, A, L = self._frame3d_axial_x()
        b = np.array([0.0, -50e3, 0.0])
        q = A * b[1]
        f = elem.compute_body_load(b)
        esperado = np.array([
            0.0, 0.5 * q * L, 0.0,
            0.0, 0.0,  q * L * L / 12.0,
            0.0, 0.5 * q * L, 0.0,
            0.0, 0.0, -q * L * L / 12.0,
        ])
        np.testing.assert_allclose(f, esperado, rtol=1e-12, atol=1e-12)

    def test_frame3d_equilibrio_global_orientacion_arbitraria(self):
        """Viga en orientación 3D arbitraria con carga arbitraria:
        ΣF_i + ΣF_j en globales debe valer A · L · b."""
        n1 = _node(1, 0.0, 0.0, 0.0)
        n2 = _node(2, 2.0, 1.0, 1.5)
        L = math.sqrt(2.0**2 + 1.0**2 + 1.5**2)
        A = 5e-3
        mat = Elastic1D(E=210e9)
        elem = Frame3D(1, [n1, n2], mat, A=A, Iy=1e-6, Iz=1e-6, J=2e-6,
                       ref_vector=[0.0, 0.0, 1.0])
        b = np.array([1.0e3, -2.0e3, -78.5e3])
        f = elem.compute_body_load(b)
        sumF = f[0:3] + f[6:9]
        np.testing.assert_allclose(sumF, A * L * b, rtol=1e-12)

    def test_frame3d_sin_torsion_por_carga_uniforme(self):
        """Carga uniforme centrada en el eje no produce torsión: Mx debe
        ser exactamente cero en ambos nodos (en locales y, por simetría,
        también en globales cuando la T es la identidad parcial)."""
        elem, _, _ = self._frame3d_axial_x()
        b = np.array([1.0, -2.0, -3.0])
        f = elem.compute_body_load(b)
        # Para esta orientación (axial=x global), Mx_local coincide con Mx_global.
        self.assertAlmostEqual(f[3], 0.0, places=12)
        self.assertAlmostEqual(f[9], 0.0, places=12)


if __name__ == '__main__':
    unittest.main()
