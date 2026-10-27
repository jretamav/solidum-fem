"""Verificación por soluciones manufacturadas de los sólidos planos (ADR 0022).

Se inventa un desplazamiento suave y no polinómico (ningún elemento lo
reproduce), SymPy deriva la fuerza de cuerpo ``b* = −∇·σ*`` y la tracción
``t* = σ*·n``, y se resuelve en mallas cada vez más finas. Con ``p`` el grado
del polinomio completo que reproduce el elemento, la teoría da

    ‖u* − u_h‖_L2 ~ h^(p+1),        ‖u* − u_h‖_E ~ h^p.

Una sola prueba ejerce a la vez lo que las soluciones analíticas comprueban
por separado: malla distorsionada con elementos de lados curvos, fuerza de
cuerpo y tracción variables, Dirichlet no homogéneo, contorno curvo y
material anisótropo girado. Un error de un factor 2, de un signo o de una
cuadratura en cualquiera de esas piezas rompe el orden.

**Criterio de aceptación** (decisión del usuario, 2026-10-27): el orden
observado entre los dos últimos niveles es ≥ teórico − 0.1, en L² y en
energía, y el error baja en cada refinamiento. Es de un solo lado para no
penalizar una superconvergencia.

Geometrías (mallas asintóticamente afines, así que también los serendípitos
deben dar el orden completo; Arnold, Boffi y Falk 2002):

- cuadrado distorsionado por ``x = X + a·sin(πX)·(1 − Y/2)``,
  ``y = Y + a·sin(πY)·(1 − X/2)``: lados rectos con los nodos deslizados a lo
  largo de ellos (nodos medios descentrados en el contorno, el caso de la
  deuda #13) y elementos interiores curvos; Dirichlet en ``x = 0`` e
  ``y = 0``, tracción en los otros dos lados;
- cuarto de corona ``1 ≤ r ≤ 2``: tracción sobre los dos arcos y Dirichlet
  en los lados rectos.

El barrido recorre el registro: todo sólido plano de desplazamientos que se
registre queda verificado sin escribir nada aquí
(``tests/test_reference_nodes_contract.py`` impide que uno se quede fuera por
no declarar su celda). Con el reparto antiguo 1-4-1 de la carga de borde
(deuda #13), Quad8 y Tri6 caían de orden 3 a orden 2 en L² en las dos
geometrías (medido el 2026-10-27).
"""
import logging
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import solidum  # noqa: F401,E402
from solidum.materials.elastic_2d import Elastic2D  # noqa: E402
from solidum.materials.orthotropic_2d import Orthotropic2D  # noqa: E402
from solidum.mesh import structured_mesh  # noqa: E402
from solidum.registry import ElementRegistry  # noqa: E402
from solidum.verification import (  # noqa: E402
    ManufacturedElasticity,
    convergence_study,
    interpolation_degree,
    linear_elastic_matrix,
)

U_ESTRELLA = ("1e-3*sin(1.1*x + 0.5*y + 0.2)*exp(0.3*y)",
              "1e-3*cos(0.7*x - 1.2*y)*(1 + 0.2*x)")
U_LINEAL = ("1e-3*(0.5 + 2*x - y)", "1e-3*(-0.3 + x + 3*y)")
A = 0.08
TOLERANCIA_ORDEN = 0.1

MATERIALES = {
    "deformación plana": Elastic2D(200e3, 0.3, "plane_strain"),
    "esfuerzo plano": Elastic2D(200e3, 0.3, "plane_stress"),
    "ortótropo a 30°": Orthotropic2D(150e3, 50e3, 30e3, 0.25, theta=np.radians(30.0)),
}


def distorsion(X):
    x, y = X
    return np.array([x + A * np.sin(np.pi * x) * (1.0 - y / 2.0),
                     y + A * np.sin(np.pi * y) * (1.0 - x / 2.0)])


def corona(X):
    r, t = 1.0 + X[0], X[1] * np.pi / 2.0
    return np.array([r * np.cos(t), r * np.sin(t)])


GEOMETRIAS = {
    "cuadrado distorsionado": (distorsion, ("xmin", "ymin")),
    "corona": (corona, ("ymin", "ymax")),
}


def solidos_planos():
    """Elementos del registro con celda de referencia plana y DOFs (ux, uy)."""
    return [n for n in ElementRegistry.names()
            if ElementRegistry.get(n).REFERENCE_CELL in ("quad", "tri")
            and list(ElementRegistry.get(n).DOF_NAMES) == ["ux", "uy"]
            and not ElementRegistry.get(n).GEOMETRICALLY_NONLINEAR]


def variantes():
    """(elemento, kwargs): la regla por omisión de cada elemento, más el Quad8
    con 2×2, la mitigación del bloqueo que recomienda la deuda #6."""
    casos = [(n, {}) for n in solidos_planos()]
    return casos + [("Quad8", {"quadrature": "2x2"})]


def grado(nombre):
    mat = MATERIALES["deformación plana"]
    el = next(iter(structured_mesh(nombre, (1, 1), material=mat).domain.elements.values()))
    return interpolation_degree(el)


def setUpModule():
    logging.disable(logging.WARNING)


def tearDownModule():
    logging.disable(logging.NOTSET)


class TestOrdenDeConvergencia(unittest.TestCase):

    def test_hay_elementos_que_barrer(self):
        self.assertGreaterEqual(set(solidos_planos()), {"Quad4", "Quad8", "Quad9", "Tri3", "Tri6"})

    def test_orden_teorico(self):
        for nombre, kwargs in variantes():
            p = grado(nombre)
            niveles = (4, 8, 16, 32) if p == 1 else (2, 4, 8, 16)
            for geo, (phi, dirichlet) in GEOMETRIAS.items():
                for nmat, mat in MATERIALES.items():
                    if geo == "corona" and nmat != "deformación plana":
                        continue
                    with self.subTest(elemento=nombre, **kwargs, geometria=geo, material=nmat):
                        sol = ManufacturedElasticity.from_material(U_ESTRELLA, mat)
                        st = convergence_study(nombre, sol, niveles, mapping=phi,
                                               dirichlet=dirichlet, material=mat, **kwargs)
                        for norma in ("L2", "energy"):
                            self.assertTrue(np.all(np.diff(st.errors[norma]) < 0.0), st.table())
                            self.assertGreaterEqual(st.orders(norma)[-1],
                                                    st.expected_order(norma) - TOLERANCIA_ORDEN,
                                                    st.table())


class TestPruebaDeLaParcela(unittest.TestCase):
    """Con ``u*`` lineal todo elemento isoparamétrico es exacto, también con
    lados curvos: el error es redondeo. Separa un fallo de las cargas o del
    contorno (no exacto aquí) de uno de la aproximación."""

    def test_campo_lineal_exacto(self):
        mat = MATERIALES["ortótropo a 30°"]
        sol = ManufacturedElasticity.from_material(U_LINEAL, mat)
        for nombre, kwargs in variantes():
            for geo, (phi, dirichlet) in GEOMETRIAS.items():
                with self.subTest(elemento=nombre, **kwargs, geometria=geo):
                    st = convergence_study(nombre, sol, (3,), mapping=phi,
                                           dirichlet=dirichlet, material=mat, **kwargs)
                    self.assertLess(st.errors["L2"][0], 1e-14)
                    self.assertLess(st.errors["energy"][0], 1e-9)


class _NoLineal:
    STRAIN_DIM = 3

    def compute_state(self, strain, state_vars=None):
        return np.asarray(strain) + 1e3 * np.asarray(strain) ** 2, np.eye(3), state_vars


class TestSolucionManufacturada(unittest.TestCase):

    def test_fuerza_de_cuerpo_y_traccion(self):
        """``u* = (x², x·y)``: ``ε* = (2x, x, y)`` y ``b* = −(2C₁₁ + C₁₂ + C₃₃, 0)``
        con ``C`` isótropa (``C₁₃ = C₂₃ = 0``)."""
        mat = MATERIALES["esfuerzo plano"]
        C = mat.C
        sol = ManufacturedElasticity.from_material(("x**2", "x*y"), mat)
        x = np.array([0.7, -1.3])
        np.testing.assert_allclose(sol.strain(x), [2 * x[0], x[0], x[1]], rtol=1e-14)
        np.testing.assert_allclose(sol.body_force(x),
                                   [-(2 * C[0, 0] + C[0, 1] + C[2, 2]), 0.0], rtol=1e-14, atol=1e-9)
        s = C @ sol.strain(x)
        n = np.array([0.6, 0.8])
        np.testing.assert_allclose(sol.traction(x, n),
                                   [s[0] * n[0] + s[2] * n[1], s[2] * n[0] + s[1] * n[1]], rtol=1e-14)

    def test_grado_de_interpolacion(self):
        esperado = {"Quad4": 1, "Tri3": 1, "Quad8": 2, "Quad9": 2, "Tri6": 2}
        for nombre, p in esperado.items():
            with self.subTest(elemento=nombre):
                self.assertEqual(grado(nombre), p)

    def test_rechaza_material_no_lineal(self):
        with self.assertRaisesRegex(ValueError, "elástico lineal"):
            linear_elastic_matrix(_NoLineal())

    def test_rechaza_material_distinto_del_de_la_solucion(self):
        sol = ManufacturedElasticity.from_material(U_ESTRELLA, MATERIALES["esfuerzo plano"])
        with self.assertRaisesRegex(ValueError, "no es la de la solución"):
            convergence_study("Quad4", sol, (2,), material=MATERIALES["deformación plana"])


if __name__ == "__main__":
    unittest.main()
