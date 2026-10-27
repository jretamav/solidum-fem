"""Verificación por soluciones manufacturadas de barras y marcos (ADR 0022, fase 4).

Línea recta inclinada (en 3D, orientada al azar respecto de los ejes
globales), con la solución exacta impuesta en el nodo inicial y, en el
final, sus fuerzas internas como cargas nodales; la carga de línea
``q = −(fuerza interna)'`` entra como fuerza de cuerpo ``b(x)``.

**Exactitud nodal** (Tong 1969). En barras y vigas de sección constante las
funciones de forma son soluciones de la ecuación sin carga (lineales en
axial y torsión, cúbicas de Hermite en flexión, las de Friedman y Kosmatka
en el Timoshenko), así que los desplazamientos y giros nodales son exactos
si la carga se integra exactamente: con campos polinómicos cuya carga es de
grado ≤ 4 (lo que integra la regla de 4 puntos), el error nodal es redondeo,
con 1 y con 3 elementos. Cualquier error de signo o de factor en la rigidez,
en la carga de línea o en las cargas de extremo la rompe; y un reparto de la
carga del Timoshenko con Hermite deja un error nodal de 12 %, 2.8 % y 0.7 %
con 2, 4 y 8 elementos (medido el 2026-10-27), por eso usa sus propias
funciones.

**Orden entre nodos**, con campos suaves no polinómicos; criterio de
aceptación del usuario (≥ teórico − 0.1 entre los dos últimos niveles y error
decreciente), componente a componente:

- axial y armaduras: ``u`` orden 2 en L², ``ε`` orden 1 en energía;
- Navier-Bernoulli (Hermite): ``v``, ``w`` orden 4; giros, 3; curvaturas, 2;
- Timoshenko: ``θ`` orden 3, ``κ`` 2, ``γ`` 1 (constante por elemento) y
  ``v`` orden 2. Al refinar, ``Φ = 12EI/(G·A_s·L²)`` crece como ``1/h²`` y la
  ``v`` del elemento tiende a lineal más una burbuja que sale de los giros:
  aproxima ``v''`` con ``θ'``, que difiere en ``γ' ≠ 0``. Con cortante casi
  rígido (``Φ ≪ 1``) se comporta como Navier-Bernoulli, orden 4;
- torsión del ``Frame3D``: ``θx`` lineal (no hay torsor distribuido), exacta.
"""
import logging
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import solidum  # noqa: F401,E402
from solidum.materials.elastic import Elastic1D  # noqa: E402
from solidum.mesh import structured_mesh  # noqa: E402
from solidum.registry import ElementRegistry  # noqa: E402
from solidum.verification import (  # noqa: E402
    ManufacturedFrame,
    convergence_study,
    error_norms,
    solve_manufactured,
)

MAT = Elastic1D(210.0)
TOLERANCIA_ORDEN = 0.1
LINEA_2D = (np.array([0.5, -0.3]), np.array([2.1, 0.9]))
LINEA_3D = (np.array([0.5, -0.3, 0.2]), np.array([1.7, 0.9, 1.0]))

U_P = "1e-3*(1 + s - 0.3*s**2 + 0.1*s**4)"
V_P = "1e-3*(0.5 - s + 0.4*s**3 - 0.05*s**6)"
W_P = "1e-3*(-0.2 + 0.7*s - 0.3*s**2 + 0.04*s**5)"
PSI_P = "1e-3*(0.5 - s + 0.4*s**3 - 0.02*s**7)"
U_S = "1e-3*sin(1.3*s + 0.2)"
V_S = "1e-3*cos(1.7*s)*exp(0.3*s)"
W_S = "1e-3*sin(0.9*s - 0.4)*(1 + 0.2*s)"
TX = "1e-3*(0.3 - 0.8*s)"

# Por elemento: argumentos, línea, campos polinómicos y suaves, órdenes
# esperados por componente y componentes exactas.
CASOS = {
    "Truss2D": dict(kw=dict(A=2.0), linea=LINEA_2D, poli={"u": U_P}, suave={"u": U_S},
                    ordenes={"L2:u": 2, "energy:eps": 1}),
    "Truss3D": dict(kw=dict(A=2.0), linea=LINEA_3D, poli={"u": U_P}, suave={"u": U_S},
                    ordenes={"L2:u": 2, "energy:eps": 1}),
    "Frame2DEuler": dict(
        kw=dict(A=2.0, I=0.3), linea=LINEA_2D, poli={"u": U_P, "v": V_P},
        suave={"u": U_S, "v": V_S},
        ordenes={"L2:u": 2, "L2:v": 4, "L2:theta_z": 3, "energy:eps": 1, "energy:kappa_z": 2}),
    "Frame2DTimoshenko": dict(
        kw=dict(A=2.0, I=0.3, As=1.6, nu=0.25), linea=LINEA_2D, poli={"u": U_P, "psi": PSI_P},
        suave={"u": U_S, "psi": V_S},
        ordenes={"L2:u": 2, "L2:v": 2, "L2:theta_z": 3, "energy:eps": 1, "energy:kappa_z": 2,
                 "energy:gamma_y": 1}),
    "Frame3D": dict(
        kw=dict(A=2.0, Iy=0.3, Iz=0.2, J=0.4, nu=0.25), linea=LINEA_3D,
        poli={"u": U_P, "v": V_P, "w": W_P, "theta_x": TX},
        suave={"u": U_S, "v": V_S, "w": W_S, "theta_x": TX},
        ordenes={"L2:u": 2, "L2:v": 4, "L2:w": 4, "L2:theta_y": 3, "L2:theta_z": 3,
                 "energy:eps": 1, "energy:kappa_y": 2, "energy:kappa_z": 2},
        exactas=("L2:theta_x", "energy:twist")),
}


def recta(linea):
    X0, X1 = linea
    return lambda X: X0 + X[0] * (X1 - X0)


def lineales():
    """Elementos de línea lineales del registro (programa principal)."""
    return [n for n in ElementRegistry.names()
            if ElementRegistry.get(n).REFERENCE_CELL == "line"
            and ElementRegistry.get(n).__module__.startswith("solidum.elements.")
            and not ElementRegistry.get(n).GEOMETRICALLY_NONLINEAR]


def solucion(nombre, campos, **cambios):
    caso = CASOS[nombre]
    kw = {**caso["kw"], **cambios}
    m = structured_mesh(nombre, (1,), mapping=recta(caso["linea"]), material=MAT, **kw)
    return ManufacturedFrame.from_element(campos, next(iter(m.domain.elements.values())),
                                          caso["linea"][0]), kw


def setUpModule():
    logging.disable(logging.WARNING)


def tearDownModule():
    logging.disable(logging.NOTSET)


class TestExactitudNodal(unittest.TestCase):

    def test_todo_elemento_lineal_tiene_su_caso(self):
        self.assertGreaterEqual(set(lineales()), set(CASOS))
        self.assertEqual(set(lineales()) - set(CASOS), set(),
                         "elemento de línea sin caso: añadir sus campos y órdenes a CASOS")

    def test_campos_polinomicos(self):
        for nombre in CASOS:
            sol, kw = solucion(nombre, CASOS[nombre]["poli"])
            for n in (1, 3):
                with self.subTest(elemento=nombre, elementos=n):
                    m = structured_mesh(nombre, (n,), mapping=recta(CASOS[nombre]["linea"]),
                                        material=MAT, **kw)
                    U = solve_manufactured(m, sol, ("xmin",))
                    e = error_norms(m.domain, U, sol)
                    self.assertLess(e["nodal"], 1e-12 * np.max(np.abs(U)))


class TestOrdenEntreNodos(unittest.TestCase):

    def test_orden_teorico(self):
        for nombre, caso in CASOS.items():
            with self.subTest(elemento=nombre):
                sol, kw = solucion(nombre, caso["suave"])
                st = convergence_study(nombre, sol, (4, 8, 16, 32), mapping=recta(caso["linea"]),
                                       dirichlet=("xmin",), expected=caso["ordenes"],
                                       material=MAT, **kw)
                for norma in caso["ordenes"]:
                    self.assertTrue(np.all(np.diff(st.errors[norma]) < 0.0), st.table())
                    self.assertGreaterEqual(st.orders(norma)[-1],
                                            st.expected_order(norma) - TOLERANCIA_ORDEN,
                                            st.table())
                # Los nodos, exactos salvo la cuadratura de la carga.
                self.assertLess(st.errors["nodal"][-1], 1e-12)
                for norma in caso.get("exactas", ()):
                    self.assertLess(st.errors[norma][-1], 1e-12, norma)

    def test_timoshenko_esbelto_como_navier_bernoulli(self):
        """Con ``A_s`` enorme (``Φ ≪ 1`` en todos los niveles) la flecha
        converge con orden 4, como con Hermite."""
        sol, kw = solucion("Frame2DTimoshenko", CASOS["Frame2DTimoshenko"]["suave"], As=1e5)
        st = convergence_study("Frame2DTimoshenko", sol, (2, 4, 8),
                               mapping=recta(LINEA_2D), dirichlet=("xmin",),
                               expected={"L2:v": 4}, material=MAT, **kw)
        self.assertGreaterEqual(st.orders("L2:v")[-1], 4 - TOLERANCIA_ORDEN, st.table())


class TestSolucionManufacturada(unittest.TestCase):

    def test_voladizo_con_carga_en_el_extremo(self):
        """``v = P/(6EI)·(3L·s² − s³)``: la flecha de un voladizo con una
        fuerza ``P`` en ``+y`` en el extremo libre. Sin carga de línea; en el
        extremo, ``Vy = P`` y ``Mz = 0``; en el empotramiento,
        ``Mz = P·L``. En el plano ``xz``, ``w`` igual con ``P`` en ``+z``:
        ``Vz = P``, ``My(0) = −P·L``."""
        P, L = 3.0, 2.0
        sol, _ = solucion("Frame3D", {})
        EIy, EIz = sol.rigidities["kappa_y"], sol.rigidities["kappa_z"]
        sol, _ = solucion("Frame3D", {"v": f"{P}/(6*{EIz})*(3*{L}*s**2 - s**3)",
                                      "w": f"{P}/(6*{EIy})*(3*{L}*s**2 - s**3)"})
        X0, X1 = LINEA_3D
        e = (X1 - X0) / np.linalg.norm(X1 - X0)
        punto = lambda s: X0 + s * e
        np.testing.assert_allclose(sol._eval("line_load", punto(0.7)), 0.0, atol=1e-12)
        np.testing.assert_allclose(sol._eval("end_forces", punto(L)), [0, P, P, 0, 0, 0],
                                   atol=1e-12)
        np.testing.assert_allclose(sol._eval("end_forces", punto(0.0)),
                                   [0, P, P, 0, -P * L, P * L], atol=1e-12)

    def test_timoshenko_sin_momento_distribuido(self):
        """De ``ψ`` salen ``θ = ψ'`` y ``v = ψ − (EI/GA_s)·ψ''``:
        ``Mz' = −Vy`` sin momento distribuido y ``q_y = EI·ψ''''``."""
        sol, _ = solucion("Frame2DTimoshenko", {"psi": "s**5"})
        EI, GA = sol.rigidities["kappa_z"], sol.rigidities["gamma_y"]
        X0, X1 = LINEA_2D
        e = (X1 - X0) / np.linalg.norm(X1 - X0)
        x = X0 + 0.8 * e
        s = 0.8
        np.testing.assert_allclose(sol.local(["v", "theta_z", "kappa_z", "gamma_y"], x[None])[0],
                                   [s**5 - EI / GA * 20 * s**3, 5 * s**4, 20 * s**3,
                                    -EI / GA * 60 * s**2], rtol=1e-13)
        N, Vy, _, _, _, Mz = sol._eval("end_forces", x)
        self.assertAlmostEqual(Vy, -EI * 60 * s**2, places=9)
        self.assertAlmostEqual(sol._eval("line_load", x)[1], EI * 120 * s, places=9)

    def test_rechazos(self):
        with self.assertRaisesRegex(ValueError, "no admitidos"):
            solucion("Frame2DTimoshenko", {"v": "s**2"})
        with self.assertRaisesRegex(ValueError, "lineal"):
            solucion("Frame3D", {"theta_x": "s**2"})
        sol, kw = solucion("Frame2DEuler", CASOS["Frame2DEuler"]["suave"])
        with self.assertRaisesRegex(ValueError, "rigideces"):
            convergence_study("Frame2DEuler", sol, (2,), mapping=recta(LINEA_2D),
                              dirichlet=("xmin",), material=MAT, A=2.0, I=0.5)
        with self.assertRaisesRegex(ValueError, "ejes"):
            convergence_study("Frame2DEuler", sol, (2,), mapping=recta(LINEA_3D[::-1]),
                              dirichlet=("xmin",), material=MAT, **kw)


if __name__ == "__main__":
    unittest.main()
