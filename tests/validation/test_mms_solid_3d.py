"""Verificación por soluciones manufacturadas de los sólidos 3D (ADR 0022, fase 2).

Mismo método y mismo criterio que ``test_mms_solid_2d.py``: un desplazamiento
suave y no polinómico, ``b* = −∇·σ*`` y ``t* = σ*·n`` derivados con SymPy, y
el orden de convergencia observado al refinar,

    ‖u* − u_h‖_L2 ~ h^(p+1),        ‖u* − u_h‖_E ~ h^p.

**Criterio de aceptación** (decisión del usuario, 2026-10-27): el orden
observado entre los dos últimos niveles es ≥ teórico − 0.1, en L² y en
energía, y el error baja en cada refinamiento.

Geometrías (mallas asintóticamente afines):

- cubo distorsionado: cada coordenada se desplaza ``a·sin(π·X_i)`` por un
  factor de las otras dos, de modo que las caras siguen planas con los nodos
  deslizados sobre ellas (nodos intermedios descentrados en el contorno) y
  los elementos interiores son curvos; Dirichlet en ``x = 0``, ``y = 0`` y
  ``z = 0``, tracción en las otras tres caras;
- sector de cilindro hueco ``1 ≤ r ≤ 2``, ``0 ≤ θ ≤ π/2``, ``0 ≤ z ≤ 1``:
  tracción sobre las dos caras cilíndricas y la superior, Dirichlet en las
  dos caras radiales y la inferior.

Niveles: n = 2, 4, 8 para los cuadráticos y 4, 8, 16 para los lineales. El
Tet4 converge al orden asintótico desde abajo y despacio (L²: 1.50, 1.75,
1.91 en el cubo; medido el 2026-10-27): con n = 8 todavía no cumpliría el
criterio. El barrido recorre el registro: todo sólido 3D que se registre
queda verificado sin escribir nada aquí.

La prueba de la parcela (``u*`` lineal) es exacta en una malla afín con
todas las reglas. En las mallas curvas la dejan de pasar el Tet10 con su
regla por omisión (``tet_4``) y el Hex20 con 2×2×2: ninguna integra
exactamente ``Bᵀ·detJ`` sobre un elemento de lados curvos. Queda un error de
consistencia que baja como h^3.3 y que no estropea la convergencia; con la
regla completa la parcela es exacta (``SOLO_AFIN``).
"""
import logging
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import solidum  # noqa: F401,E402
from solidum.materials.elastic_3d import Elastic3D  # noqa: E402
from solidum.mesh import structured_mesh  # noqa: E402
from solidum.registry import ElementRegistry  # noqa: E402
from solidum.verification import (  # noqa: E402
    ManufacturedElasticity,
    convergence_study,
    interpolation_degree,
)

U_ESTRELLA = ("1e-3*sin(1.1*x + 0.5*y + 0.3*z + 0.2)*exp(0.3*y)",
              "1e-3*cos(0.7*x - 1.2*y + 0.4*z)*(1 + 0.2*x)",
              "1e-3*sin(0.9*x + 0.6*y - 1.1*z)*(1 + 0.1*y)")
U_LINEAL = ("1e-3*(0.5 + 2*x - y + 0.3*z)", "1e-3*(-0.3 + x + 3*y - z)",
            "1e-3*(0.1 - 0.4*x + 0.2*y + 1.5*z)")
A = 0.08
TOLERANCIA_ORDEN = 0.1
MAT = Elastic3D(200e3, 0.3)


def distorsion(X):
    x, y, z = X
    return np.array([x + A * np.sin(np.pi * x) * (1.0 - y / 2.0) * (1.0 - z / 3.0),
                     y + A * np.sin(np.pi * y) * (1.0 - z / 2.0) * (1.0 - x / 3.0),
                     z + A * np.sin(np.pi * z) * (1.0 - x / 2.0) * (1.0 - y / 3.0)])


def cilindro(X):
    r, t = 1.0 + X[0], X[1] * np.pi / 2.0
    return np.array([r * np.cos(t), r * np.sin(t), X[2]])


GEOMETRIAS = {
    "cubo distorsionado": (distorsion, ("xmin", "ymin", "zmin")),
    "cilindro hueco": (cilindro, ("ymin", "ymax", "zmin")),
}


def solidos_3d():
    """Elementos del registro con celda de referencia 3D y DOFs (ux, uy, uz)."""
    return [n for n in ElementRegistry.names()
            if ElementRegistry.get(n).REFERENCE_CELL in ("hex", "tet")
            and list(ElementRegistry.get(n).DOF_NAMES) == ["ux", "uy", "uz"]
            and not ElementRegistry.get(n).GEOMETRICALLY_NONLINEAR]


def variantes():
    """(elemento, kwargs): la regla por omisión de cada elemento, más el Hex20
    con 2×2×2 (integración reducida, como el Quad8 2×2 en 2D)."""
    return [(n, {}) for n in solidos_3d()] + [("Hex20", {"quadrature": "hex_2x2x2"})]


def grado(nombre):
    el = next(iter(structured_mesh(nombre, (1, 1, 1), material=MAT).domain.elements.values()))
    return interpolation_degree(el)


def setUpModule():
    logging.disable(logging.WARNING)


def tearDownModule():
    logging.disable(logging.NOTSET)


class TestOrdenDeConvergencia(unittest.TestCase):

    def test_hay_elementos_que_barrer(self):
        self.assertGreaterEqual(set(solidos_3d()), {"Hex8", "Hex20", "Hex27", "Tet4", "Tet10"})

    def test_orden_teorico(self):
        sol = ManufacturedElasticity.from_material(U_ESTRELLA, MAT)
        for nombre, kwargs in variantes():
            niveles = (4, 8, 16) if grado(nombre) == 1 else (2, 4, 8)
            for geo, (phi, dirichlet) in GEOMETRIAS.items():
                with self.subTest(elemento=nombre, **kwargs, geometria=geo):
                    st = convergence_study(nombre, sol, niveles, mapping=phi,
                                           dirichlet=dirichlet, material=MAT, **kwargs)
                    for norma in ("L2", "energy"):
                        self.assertTrue(np.all(np.diff(st.errors[norma]) < 0.0), st.table())
                        self.assertGreaterEqual(st.orders(norma)[-1],
                                                st.expected_order(norma) - TOLERANCIA_ORDEN,
                                                st.table())


# Reglas de rigidez que no integran exactamente ``Bᵀ·detJ`` sobre un elemento
# de lados curvos (el integrando es de grado mayor que el que integran): con
# ellas la prueba de la parcela sólo es exacta en mallas de lados rectos. En
# las curvas queda un error de consistencia que baja como h^3.3 (medido el
# 2026-10-27) y no estropea el orden de convergencia, que el barrido comprueba
# con estas mismas reglas. Con la regla completa (el valor) la parcela es
# exacta también en las mallas curvas.
SOLO_AFIN = {("Tet10", "tet_4"): "tet_15", ("Hex20", "hex_2x2x2"): "hex_3x3x3"}


def afin(X):
    return np.array([[1.0, 0.3, 0.1], [0.2, 1.2, -0.2], [0.0, 0.25, 0.9]]) @ X


class TestPruebaDeLaParcela(unittest.TestCase):
    """Con ``u*`` lineal todo elemento isoparamétrico es exacto: el error es
    redondeo. Separa un fallo de las cargas o del contorno (no exacto aquí)
    de uno de la aproximación."""

    def error(self, nombre, phi, dirichlet, **kwargs):
        sol = ManufacturedElasticity.from_material(U_LINEAL, MAT)
        st = convergence_study(nombre, sol, (2,), mapping=phi, dirichlet=dirichlet,
                               material=MAT, **kwargs)
        return st.errors["L2"][0], st.errors["energy"][0]

    def test_exacta_en_malla_afin(self):
        for nombre, kwargs in variantes():
            with self.subTest(elemento=nombre, **kwargs):
                e_l2, e_en = self.error(nombre, afin, ("xmin", "ymin", "zmin"), **kwargs)
                self.assertLess(e_l2, 1e-14)
                self.assertLess(e_en, 1e-9)

    def test_mallas_curvas(self):
        for nombre, kwargs in variantes():
            regla = kwargs.get("quadrature") or next(iter(structured_mesh(
                nombre, (1, 1, 1), material=MAT).domain.elements.values())).quadrature_key
            for geo, (phi, dirichlet) in GEOMETRIAS.items():
                with self.subTest(elemento=nombre, regla=regla, geometria=geo):
                    e_l2, _ = self.error(nombre, phi, dirichlet, **kwargs)
                    completa = SOLO_AFIN.get((nombre, regla))
                    if completa is None:
                        self.assertLess(e_l2, 1e-14)
                    else:
                        self.assertGreater(e_l2, 1e-10)
                        self.assertLess(self.error(nombre, phi, dirichlet, quadrature=completa)[0],
                                        1e-14)


class TestSolucionManufacturada(unittest.TestCase):

    def test_cortantes_en_el_orden_voigt_del_proyecto(self):
        """``u* = (0, 0, x² + 3y²)``: ``γ_yz = 6y`` y ``γ_xz = 2x`` en las
        posiciones 4 y 5 del Voigt del proyecto, ``σ_xz = 2Gx``,
        ``σ_yz = 6Gy`` y ``b*_z = −(∂σ_xz/∂x + ∂σ_yz/∂y) = −8G``. Con
        ``γ_yz`` y ``γ_xz`` permutados saldría ``b*_z = 0``."""
        G = MAT.C[3, 3]
        sol = ManufacturedElasticity.from_material(("0", "0", "x**2 + 3*y**2"), MAT)
        x = np.array([0.7, -1.3, 0.4])
        np.testing.assert_allclose(sol.strain(x), [0, 0, 0, 0, 6 * x[1], 2 * x[0]], atol=1e-15)
        np.testing.assert_allclose(sol.body_force(x), [0.0, 0.0, -8.0 * G], rtol=1e-14, atol=1e-9)
        np.testing.assert_allclose(sol.traction(x, [1.0, 0.0, 0.0]), [0.0, 0.0, 2 * G * x[0]],
                                   rtol=1e-14, atol=1e-9)
        np.testing.assert_allclose(sol.traction(x, [0.0, 1.0, 0.0]), [0.0, 0.0, 6 * G * x[1]],
                                   rtol=1e-14, atol=1e-9)

    def test_grado_de_interpolacion(self):
        for nombre, p in {"Hex8": 1, "Tet4": 1, "Hex20": 2, "Hex27": 2, "Tet10": 2}.items():
            with self.subTest(elemento=nombre):
                self.assertEqual(grado(nombre), p)


if __name__ == "__main__":
    unittest.main()
