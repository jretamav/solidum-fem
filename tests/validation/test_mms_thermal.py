"""Verificación por soluciones manufacturadas de la conducción de calor (ADR 0022, fase 3).

Mismo método que ``test_mms_solid_2d.py``: una temperatura suave y no
polinómica ``T*``, y SymPy deriva la fuente y el flujo de frontera que la
hacen exacta con conductividad **anisótropa** (tensor lleno):

    q* = −k·∇T*,     Q* = ρc·∂T*/∂t + ∇·q*,     q̄*(x, n) = q*·n  (> 0 saliente).

**Criterio de aceptación** (decisión del usuario, 2026-10-27): el orden
observado entre los dos últimos niveles es ≥ teórico − 0.1, en L² y en
energía (``‖e‖_E² = ∫ ∇e·k·∇e dΩ``), y el error baja en cada refinamiento.

Estacionario: ``‖e‖_L2 ~ h^(p+1)`` y ``‖e‖_E ~ h^p`` en las geometrías de los
sólidos (cuadrado y cubo distorsionados, corona y cilindro hueco), con
Dirichlet en unos lados y flujo variable en los demás.

Transitorio (método θ), con los dos errores separados:

- **en el espacio**: Crank-Nicolson con ``Δt = h/2`` hasta ``t = 0.5``; el
  error en el tiempo, ``O(Δt²) = O(h²)``, no altera los órdenes del espacio.
  Con la capacidad concentrada (la de omisión) y con la consistente;
- **en el tiempo**: ``T*`` lineal en ``x`` con cualquier dependencia de
  ``t``, sobre una malla fija. La semidiscretización en el espacio es exacta
  (el elemento reproduce ``T*`` y la capacidad consistente y la fuente usan la
  misma regla), así que el error en ``t_end`` es sólo el del esquema: orden 2
  con Crank-Nicolson y 1 con θ = 1, 2/3 y 0 (éste con ``Δt`` dentro de su
  límite de estabilidad, ``Δt_crit ≈ 0.0055`` en la malla 4×4, medido el
  2026-10-27);
- **parcela en el tiempo**: con ``T*`` lineal en ``x`` y en ``t`` el método θ es
  exacto para todo θ, porque la solución exacta cumple la ecuación del paso
  (``(T_{n+1} − T_n)/Δt = Ṫ`` y ``θF_{n+1} + (1 − θ)F_n = F(t_n + θΔt)``). Un
  error en la ponderación de la carga o del Dirichlet variable en el tiempo
  (incluido el término ``C_fp·(ḡ_{n+1} − ḡ_n)``) rompería la igualdad.

El barrido recorre el registro: todo elemento térmico que se registre queda
verificado sin escribir nada aquí.
"""
import logging
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import solidum  # noqa: F401,E402
from solidum.materials.thermal_conduction import ThermalConduction  # noqa: E402
from solidum.mesh import structured_mesh  # noqa: E402
from solidum.registry import ElementRegistry  # noqa: E402
from solidum.verification import (  # noqa: E402
    ManufacturedConduction,
    convergence_study,
    interpolation_degree,
    time_convergence_study,
)

T_ESTRELLA = {
    2: "sin(1.1*x + 0.5*y + 0.2)*exp(0.3*y) + 0.5*cos(0.7*x - 1.2*y)",
    3: "sin(1.1*x + 0.5*y + 0.3*z + 0.2)*exp(0.3*y) + 0.5*cos(0.7*x - 1.2*y + 0.4*z)",
}
T_TRANSITORIA = {
    2: "(sin(1.1*x + 0.5*y + 0.2)*exp(0.3*y) + 0.5*cos(0.7*x - 1.2*y))*exp(-0.5*t)"
       " + 0.3*x*y*sin(2*t)",
    3: "(sin(1.1*x + 0.5*y + 0.3*z + 0.2)*exp(0.3*y) + 0.5*cos(0.7*x - 1.2*y + 0.4*z))"
       "*exp(-0.5*t) + 0.3*x*y*z*sin(2*t)",
}
T_LINEAL = {2: "1 + 2*x - 0.7*y", 3: "1 + 2*x - 0.7*y + 0.4*z"}
T_LINEAL_EN_X = {
    2: "(1 + 2*x - 0.7*y)*cos(3*t) + (0.5 - x + 1.5*y)*exp(-t) + 2*sin(4*t)",
    3: "(1 + 2*x - 0.7*y + 0.4*z)*cos(3*t) + (0.5 - x + 1.5*y - z)*exp(-t) + 2*sin(4*t)",
}
T_LINEAL_EN_X_Y_T = {2: "1 + 2*x - 0.7*y + (3 - x + 0.5*y)*t",
                     3: "1 + 2*x - 0.7*y + 0.4*z + (3 - x + 0.5*y - 0.2*z)*t"}

MAT = {2: ThermalConduction([[3.0, 1.0], [1.0, 2.0]], c=2.0, density=1.5),
       3: ThermalConduction([[3.0, 0.5, 0.2], [0.5, 2.0, 0.4], [0.2, 0.4, 1.5]],
                            c=2.0, density=1.5)}
A = 0.08
TOLERANCIA_ORDEN = 0.1


def distorsion(X):
    if len(X) == 2:
        x, y = X
        return np.array([x + A * np.sin(np.pi * x) * (1.0 - y / 2.0),
                         y + A * np.sin(np.pi * y) * (1.0 - x / 2.0)])
    x, y, z = X
    return np.array([x + A * np.sin(np.pi * x) * (1.0 - y / 2.0) * (1.0 - z / 3.0),
                     y + A * np.sin(np.pi * y) * (1.0 - z / 2.0) * (1.0 - x / 3.0),
                     z + A * np.sin(np.pi * z) * (1.0 - x / 2.0) * (1.0 - y / 3.0)])


def curva(X):
    """Cuarto de corona ``1 ≤ r ≤ 2`` en 2D; sector de cilindro hueco de
    altura 1 en 3D."""
    r, t = 1.0 + X[0], X[1] * np.pi / 2.0
    return np.array([r * np.cos(t), r * np.sin(t), *X[2:]])


GEOMETRIAS = {
    2: {"cuadrado distorsionado": (distorsion, ("xmin", "ymin")),
        "corona": (curva, ("ymin", "ymax"))},
    3: {"cubo distorsionado": (distorsion, ("xmin", "ymin", "zmin")),
        "cilindro hueco": (curva, ("ymin", "ymax", "zmin"))},
}


def termicos():
    """(nombre, dimensión) de los elementos térmicos del registro con celda de
    referencia."""
    dims = {"quad": 2, "tri": 2, "hex": 3, "tet": 3}
    return [(n, dims[ElementRegistry.get(n).REFERENCE_CELL]) for n in ElementRegistry.names()
            if list(ElementRegistry.get(n).DOF_NAMES) == ["T"]
            and ElementRegistry.get(n).REFERENCE_CELL in dims]


def grado(nombre, d):
    el = next(iter(structured_mesh(nombre, (1,) * d, material=MAT[d]).domain.elements.values()))
    return interpolation_degree(el)


def niveles(nombre, d):
    return (4, 8, 16) if grado(nombre, d) == 1 else (2, 4, 8)


def setUpModule():
    logging.disable(logging.WARNING)


def tearDownModule():
    logging.disable(logging.NOTSET)


class _Criterio(unittest.TestCase):

    def cumple(self, st):
        for norma in ("L2", "energy"):
            self.assertTrue(np.all(np.diff(st.errors[norma]) < 0.0), st.table())
            self.assertGreaterEqual(st.orders(norma)[-1],
                                    st.expected_order(norma) - TOLERANCIA_ORDEN, st.table())


class TestEstacionario(_Criterio):

    def test_hay_elementos_que_barrer(self):
        self.assertGreaterEqual({n for n, _ in termicos()}, {"Quad4Thermal", "Hex8Thermal"})

    def test_orden_teorico(self):
        for nombre, d in termicos():
            sol = ManufacturedConduction.from_material(T_ESTRELLA[d], MAT[d])
            for geo, (phi, dirichlet) in GEOMETRIAS[d].items():
                with self.subTest(elemento=nombre, geometria=geo):
                    self.cumple(convergence_study(nombre, sol, niveles(nombre, d), mapping=phi,
                                                  dirichlet=dirichlet, material=MAT[d]))

    def test_parcela(self):
        """``T*`` lineal: todo elemento isoparamétrico es exacto, también con
        lados curvos; el error es redondeo."""
        for nombre, d in termicos():
            sol = ManufacturedConduction.from_material(T_LINEAL[d], MAT[d])
            for geo, (phi, dirichlet) in GEOMETRIAS[d].items():
                with self.subTest(elemento=nombre, geometria=geo):
                    st = convergence_study(nombre, sol, (2,), mapping=phi, dirichlet=dirichlet,
                                           material=MAT[d])
                    self.assertLess(st.errors["L2"][0], 1e-14)
                    self.assertLess(st.errors["energy"][0], 1e-13)


class TestTransitorioEnElEspacio(_Criterio):

    def test_orden_teorico_con_crank_nicolson(self):
        for nombre, d in termicos():
            sol = ManufacturedConduction.from_material(T_TRANSITORIA[d], MAT[d])
            n = niveles(nombre, d)
            # En 3D sólo el cubo y la capacidad de omisión: el resto lo cubren
            # el 2D y el estacionario 3D, y cada nivel 3D cuesta segundos.
            casos = ([(geo, lump) for geo in GEOMETRIAS[d] for lump in ("lumped", "consistent")]
                     if d == 2 else [("cubo distorsionado", "lumped")])
            for geo, lumping in casos:
                phi, dirichlet = GEOMETRIAS[d][geo]
                with self.subTest(elemento=nombre, geometria=geo, capacidad=lumping):
                    self.cumple(convergence_study(
                        nombre, sol, n, mapping=phi, dirichlet=dirichlet, material=MAT[d],
                        t_end=0.5, n_steps=n, theta=0.5, lumping=lumping))


class TestTransitorioEnElTiempo(_Criterio):

    PASOS = {1.0: (16, 32, 64, 128), 2.0 / 3.0: (16, 32, 64, 128), 0.5: (8, 16, 32, 64)}

    def test_orden_del_metodo_theta(self):
        for nombre, d in termicos():
            sol = ManufacturedConduction.from_material(T_LINEAL_EN_X[d], MAT[d])
            for theta, pasos in self.PASOS.items():
                if d == 3 and theta == 2.0 / 3.0:
                    continue
                with self.subTest(elemento=nombre, theta=theta):
                    st = time_convergence_study(nombre, sol, 4 if d == 2 else 2, pasos,
                                                t_end=1.0, theta=theta, mapping=distorsion,
                                                material=MAT[d])
                    self.assertEqual(st.expected_order("L2"), 2 if theta == 0.5 else 1)
                    self.cumple(st)

    def test_explicito_dentro_de_su_limite(self):
        for nombre, d in termicos():
            if d != 2:
                continue
            with self.subTest(elemento=nombre):
                sol = ManufacturedConduction.from_material(T_LINEAL_EN_X[d], MAT[d])
                self.cumple(time_convergence_study(nombre, sol, 4, (64, 128, 256), t_end=0.25,
                                                   theta=0.0, mapping=distorsion,
                                                   material=MAT[d]))

    def test_parcela_en_el_tiempo(self):
        for nombre, d in termicos():
            sol = ManufacturedConduction.from_material(T_LINEAL_EN_X_Y_T[d], MAT[d])
            for theta in (0.0, 0.5, 2.0 / 3.0, 1.0):
                with self.subTest(elemento=nombre, theta=theta):
                    pasos = 64 if theta == 0.0 else 4
                    st = time_convergence_study(nombre, sol, 3 if d == 2 else 2, (pasos,),
                                                t_end=0.25, theta=theta, mapping=distorsion,
                                                material=MAT[d])
                    self.assertLess(st.errors["L2"][0], 1e-13)
                    self.assertLess(st.errors["energy"][0], 1e-12)


class TestSolucionManufacturada(unittest.TestCase):

    def test_signos_de_la_fuente_y_del_flujo(self):
        """``T* = −(x² + y²)`` con ``k = diag(3, 2)``: un máximo en el origen
        sostenido por una fuente ``Q* = 2·3 + 2·2 = 10 > 0``, con el calor
        saliendo por el contorno: ``q* = (6x, 4y)`` y ``q̄* > 0`` con la normal
        exterior."""
        sol = ManufacturedConduction("-(x**2 + y**2)", np.diag([3.0, 2.0]))
        x = np.array([1.0, 0.5])
        self.assertAlmostEqual(sol.source(x), 10.0, places=14)
        np.testing.assert_allclose(sol.flux(x), [6.0, 2.0], rtol=1e-15)
        self.assertAlmostEqual(sol.boundary_flux(x, [1.0, 0.0]), 6.0, places=14)
        self.assertAlmostEqual(sol.boundary_flux(x, [0.0, -1.0]), -2.0, places=14)

    def test_termino_cruzado_de_la_conductividad(self):
        """``T* = x·y``: ``∇T* = (y, x)`` y ``Q* = ∇·q* = −2·k_xy``; con ``k``
        diagonal sería nula."""
        sol = ManufacturedConduction.from_material("x*y", MAT[2])
        x = np.array([0.3, -0.8])
        self.assertAlmostEqual(sol.source(x), -2.0, places=14)
        np.testing.assert_allclose(sol.flux(x), -MAT[2].conductivity @ [x[1], x[0]], rtol=1e-15)

    def test_fuente_transitoria(self):
        """``T* = t``: el cuerpo se calienta uniformemente con ``Q* = ρc``."""
        sol = ManufacturedConduction.from_material("t", MAT[3])
        self.assertTrue(sol.time_dependent)
        self.assertAlmostEqual(sol.source(np.zeros(3), 0.7), 3.0, places=14)
        self.assertAlmostEqual(sol.temperature(np.zeros(3), 0.7), 0.7, places=15)

    def test_grado_de_interpolacion(self):
        for nombre, p in {"Quad4Thermal": 1, "Hex8Thermal": 1}.items():
            with self.subTest(elemento=nombre):
                self.assertEqual(grado(nombre, 2 if nombre.startswith("Quad") else 3), p)

    def test_rechazos(self):
        with self.assertRaisesRegex(ValueError, "rho_c"):
            ManufacturedConduction("x*t", np.eye(2))
        with self.assertRaisesRegex(ValueError, "Símbolos desconocidos"):
            ManufacturedConduction("x*z", np.eye(2))
        with self.assertRaisesRegex(ValueError, "lineal en x"):
            time_convergence_study("Quad4Thermal",
                                   ManufacturedConduction.from_material(T_TRANSITORIA[2], MAT[2]),
                                   2, (4,), t_end=0.1, material=MAT[2])
        with self.assertRaisesRegex(ValueError, "conductividad del material"):
            convergence_study("Quad4Thermal", ManufacturedConduction(T_ESTRELLA[2], np.eye(2)),
                              (2,), material=MAT[2])
        with self.assertRaisesRegex(ValueError, "t_end"):
            convergence_study("Quad4Thermal",
                              ManufacturedConduction.from_material(T_TRANSITORIA[2], MAT[2]),
                              (2,), material=MAT[2])

    def test_divergencia_del_explicito(self):
        """Por encima de su límite el explícito diverge, y el estudio lo dice
        en vez de medir un error sin sentido."""
        sol = ManufacturedConduction.from_material(T_LINEAL_EN_X[2], MAT[2])
        with self.assertRaisesRegex(RuntimeError, "divergió"):
            time_convergence_study("Quad4Thermal", sol, 4, (16,), t_end=1.0, theta=0.0,
                                   mapping=distorsion, material=MAT[2])


if __name__ == "__main__":
    unittest.main()
