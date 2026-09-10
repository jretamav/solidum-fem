"""Lámina ortótropa fuera de eje — constantes aparentes vs solución cerrada.

Referencia
----------
Jones R.M. (1999). *Mechanics of Composite Materials*, 2ª ed. Taylor &
Francis. §2.6 (transformación de la constitutiva), §2.8 (constantes
aparentes de una lámina cargada fuera de sus ejes principales).

Definición
----------
Una lámina ortótropa cuyos ejes principales (1, 2) forman un ángulo θ con
los ejes de carga (x, y) responde con **constantes aparentes** que dependen
de θ. Bajo tracción uniaxial pura σ_xx, con σ_yy = σ_xy = 0, la teoría
clásica da tres cantidades medibles con forma cerrada:

**Módulo aparente** (Jones ec. 2.85):

    1/E_x(θ) = c⁴/E₁ + (1/G₁₂ − 2ν₁₂/E₁)·c²s² + s⁴/E₂

**Poisson aparente** (Jones ec. 2.86):

    ν_xy(θ) = E_x(θ)·[ ν₁₂(c⁴ + s⁴)/E₁ − (1/E₁ + 1/E₂ − 1/G₁₂)·c²s² ]

**Coeficiente de influencia mutua** — el acoplamiento tracción-cortante,
que no tiene análogo isótropo (Jones §2.8, "coefficient of mutual
influence" η_xy,x). Se mide como la distorsión angular generada por una
tracción normal pura:

    η_xy,x(θ) = γ_xy / ε_xx    bajo σ_xx puro

con c = cos θ, s = sin θ.

Por qué este benchmark
----------------------
Los tests unitarios de `tests/test_orthotropic_2d.py` son **verificación**:
comprueban que el código implementa las ecuaciones de la spec (degeneración
al isótropo, invariancia rotacional, simetría, admisibilidad). Todos
contrastan el código contra la propia formulación del proyecto.

Este archivo es **validación externa**: contrasta la formulación contra
resultados publicados independientes. Es la distinción que `Reglas.md §6`
exige para toda formulación nueva.

Las tres cantidades son sensibles a errores distintos:

1. `E_x(θ)` valida la **rotación de la constitutiva** y, en particular, el
   término cruzado `(1/G₁₂ − 2ν₁₂/E₁)`, que es donde G₁₂ deja de ser
   derivable de E y ν — el rasgo que distingue al ortótropo del isótropo.
2. `ν_xy(θ)` valida la **reciprocidad** y el acoplamiento entre las dos
   direcciones normales.
3. `η_xy,x(θ)` valida el **acoplamiento tracción-cortante**, y es el que
   caza un factor 2 mal puesto en el cortante *engineering* — el error
   clásico al programar la transformación de Voigt (ver spec §12).

Cobertura
---------
- Constantes aparentes contra forma cerrada, barriendo θ ∈ [0°, 90°].
- Consistencia FEM: el mismo E_x medido sobre un modelo `Quad4` real
  traccionado, no sólo sobre la matriz constitutiva.
- Casos límite conocidos: en θ = 0° y 90° las constantes degeneran a las
  del material y el acoplamiento se anula.
- Degeneración isótropa: con datos isótropos, E_x(θ) = E para todo θ.

Datos
-----
Bambú, órdenes de literatura general (E₁/E₂ ≈ 19). No corresponden a una
especie concreta — ver §Diálogo de docs/specs/Orthotropic2D.md.
"""
import math
import os
import sys
import unittest

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from solidum.core.node import Node
from solidum.elements.solid_2d.quad4 import Quad4
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.orthotropic_2d import Orthotropic2D


BAMBU = dict(E1=15.0e9, E2=0.8e9, G12=0.7e9, nu12=0.35)

# Barrido de ángulos. Incluye los dos ejes principales, el punto de máximo
# acoplamiento (45°) y ángulos intermedios sin simetría particular.
ANGULOS = (0.0, 5.0, 15.0, 30.0, 45.0, 60.0, 75.0, 85.0, 90.0)


# ----------------------------------------------------------------------
# Soluciones cerradas de Jones (1999) §2.8 — escritas aquí de forma
# independiente de la implementación, a partir de las ecuaciones publicadas.
# ----------------------------------------------------------------------

def _E_x_analitico(theta_deg: float, E1, E2, G12, nu12) -> float:
    """Módulo aparente en la dirección de carga. Jones ec. 2.85."""
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    inv = (c**4 / E1
           + (1.0 / G12 - 2.0 * nu12 / E1) * c**2 * s**2
           + s**4 / E2)
    return 1.0 / inv


def _nu_xy_analitico(theta_deg: float, E1, E2, G12, nu12) -> float:
    """Coeficiente de Poisson aparente. Jones ec. 2.86."""
    t = math.radians(theta_deg)
    c, s = math.cos(t), math.sin(t)
    Ex = _E_x_analitico(theta_deg, E1, E2, G12, nu12)
    return Ex * (nu12 * (c**4 + s**4) / E1
                 - (1.0 / E1 + 1.0 / E2 - 1.0 / G12) * c**2 * s**2)


# ----------------------------------------------------------------------
# Medición sobre el material: constantes aparentes desde la flexibilidad
# ----------------------------------------------------------------------

def _constantes_aparentes(material) -> tuple:
    """Extrae (E_x, nu_xy, eta_xy_x) invirtiendo C bajo tracción uniaxial.

    Bajo σ = [1, 0, 0] (tracción unitaria en x, resto nulo), la deformación
    es ε = S·σ con S = C⁻¹, luego:

        E_x      = 1 / S[0,0]
        nu_xy    = −S[1,0] / S[0,0]
        eta_xy,x =  S[2,0] / S[0,0]

    Se impone el estado de **esfuerzo** y no el de deformación porque las
    constantes aparentes de Jones se definen así: un ensayo de tracción
    uniaxial deja libres las demás componentes.
    """
    S = np.linalg.inv(material.C)
    E_x = 1.0 / S[0, 0]
    nu_xy = -S[1, 0] / S[0, 0]
    eta = S[2, 0] / S[0, 0]
    return E_x, nu_xy, eta


class TestModuloAparente(unittest.TestCase):
    """E_x(θ) contra Jones ec. 2.85."""

    def test_barrido_de_angulos(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                m = Orthotropic2D(**BAMBU, theta=theta)
                E_x, _, _ = _constantes_aparentes(m)
                ref = _E_x_analitico(theta, **BAMBU)
                self.assertLess(abs(E_x - ref) / ref, 1e-12)

    def test_degenera_a_E1_en_eje_de_fibra(self):
        m = Orthotropic2D(**BAMBU, theta=0.0)
        E_x, _, _ = _constantes_aparentes(m)
        self.assertLess(abs(E_x - BAMBU['E1']) / BAMBU['E1'], 1e-12)

    def test_degenera_a_E2_a_noventa_grados(self):
        m = Orthotropic2D(**BAMBU, theta=90.0)
        E_x, _, _ = _constantes_aparentes(m)
        self.assertLess(abs(E_x - BAMBU['E2']) / BAMBU['E2'], 1e-12)

    def test_decrecimiento_monotono_de_cero_a_noventa(self):
        """Con E1 >> E2 el módulo aparente cae monótonamente al girar."""
        valores = []
        for theta in ANGULOS:
            m = Orthotropic2D(**BAMBU, theta=theta)
            E_x, _, _ = _constantes_aparentes(m)
            valores.append(E_x)
        for previo, actual in zip(valores, valores[1:]):
            self.assertLess(actual, previo)

    def test_caida_de_rigidez_fuera_de_eje(self):
        """Cifras citables de la severidad ortótropa del bambú.

        La caída de E_x con θ es **mucho más brusca de lo que sugiere el
        cociente E1/E2**, y ésa es la consecuencia práctica que hace de la
        orientación de fibra un dato de primer orden y no un detalle:

            θ =  0°  → 15.00 GPa   (1.000 · E1)
            θ = 15°  →  6.67 GPa   (0.444 · E1)   ya menos de la mitad
            θ = 30°  →  2.67 GPa   (0.178 · E1)
            θ = 45°  →  1.48 GPa   (0.099 · E1)   menos del 10 %
            θ = 90°  →  0.80 GPa   (0.053 · E1)   = E2

        A 45° la rigidez cae por debajo de la décima parte pese a que
        E1/E2 = 18.75: el mínimo NO se alcanza interpolando entre E1 y E2,
        porque el término cruzado (1/G12 − 2ν12/E1) domina en la zona
        intermedia. Con G12 pequeño —como en bambú— ese término hunde la
        curva muy por debajo de la media de los dos módulos.
        """
        esperado = {0.0: 1.0, 15.0: 0.4444, 30.0: 0.1779,
                    45.0: 0.0988, 90.0: 0.0533}
        E_0 = _constantes_aparentes(Orthotropic2D(**BAMBU, theta=0.0))[0]
        for theta, ratio_ref in esperado.items():
            with self.subTest(theta=theta):
                E_x = _constantes_aparentes(
                    Orthotropic2D(**BAMBU, theta=theta))[0]
                self.assertLess(abs(E_x / E_0 - ratio_ref) / ratio_ref, 1e-3)

    def test_minimo_no_esta_entre_E1_y_E2_interpolados(self):
        """A 45° E_x cae por debajo de E2·(algo): el término cruzado manda.

        Blinda que la fórmula no degeneró a una interpolación ingenua entre
        los dos módulos principales, que es el error conceptual típico.
        """
        E_45 = _constantes_aparentes(Orthotropic2D(**BAMBU, theta=45.0))[0]
        media_armonica = 2.0 / (1.0 / BAMBU['E1'] + 1.0 / BAMBU['E2'])
        self.assertLess(E_45, media_armonica)


class TestPoissonAparente(unittest.TestCase):
    """nu_xy(θ) contra Jones ec. 2.86."""

    def test_barrido_de_angulos(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                m = Orthotropic2D(**BAMBU, theta=theta)
                _, nu_xy, _ = _constantes_aparentes(m)
                ref = _nu_xy_analitico(theta, **BAMBU)
                escala = max(abs(ref), 1e-3)
                self.assertLess(abs(nu_xy - ref) / escala, 1e-10)

    def test_degenera_a_nu12_en_eje(self):
        m = Orthotropic2D(**BAMBU, theta=0.0)
        _, nu_xy, _ = _constantes_aparentes(m)
        self.assertLess(abs(nu_xy - BAMBU['nu12']) / BAMBU['nu12'], 1e-12)

    def test_degenera_a_nu21_a_noventa(self):
        """A 90° el Poisson aparente es el menor, nu21 = nu12*E2/E1."""
        nu21 = BAMBU['nu12'] * BAMBU['E2'] / BAMBU['E1']
        m = Orthotropic2D(**BAMBU, theta=90.0)
        _, nu_xy, _ = _constantes_aparentes(m)
        self.assertLess(abs(nu_xy - nu21) / nu21, 1e-12)


class TestAcoplamientoTraccionCortante(unittest.TestCase):
    """eta_xy,x(θ) — el coeficiente sin análogo isótropo.

    Es el test que caza un factor 2 mal puesto en el cortante engineering.
    """

    def test_nulo_en_ejes_principales(self):
        for theta in (0.0, 90.0):
            with self.subTest(theta=theta):
                m = Orthotropic2D(**BAMBU, theta=theta)
                _, _, eta = _constantes_aparentes(m)
                self.assertLess(abs(eta), 1e-14)

    def test_no_nulo_fuera_de_ejes(self):
        for theta in (15.0, 30.0, 45.0, 60.0, 75.0):
            with self.subTest(theta=theta):
                m = Orthotropic2D(**BAMBU, theta=theta)
                _, _, eta = _constantes_aparentes(m)
                self.assertGreater(abs(eta), 1e-3)

    def test_antisimetria_respecto_a_noventa_grados(self):
        """eta(θ) y eta(−θ) tienen signos opuestos: la distorsión cambia
        de sentido al espejar la orientación de la fibra."""
        for theta in (15.0, 30.0, 45.0):
            with self.subTest(theta=theta):
                _, _, eta_pos = _constantes_aparentes(
                    Orthotropic2D(**BAMBU, theta=theta))
                _, _, eta_neg = _constantes_aparentes(
                    Orthotropic2D(**BAMBU, theta=-theta))
                escala = abs(eta_pos)
                self.assertLess(abs(eta_pos + eta_neg) / escala, 1e-12)

    def test_nulo_para_material_isotropo(self):
        """Ningún isótropo acopla tracción con cortante, gire lo que gire."""
        E, nu = 210.0e9, 0.3
        G = E / (2.0 * (1.0 + nu))
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                m = Orthotropic2D(E1=E, E2=E, G12=G, nu12=nu, theta=theta)
                _, _, eta = _constantes_aparentes(m)
                self.assertLess(abs(eta), 1e-14)


class TestConsistenciaFEM(unittest.TestCase):
    """El módulo aparente medido sobre un modelo Quad4 real.

    Valida no sólo la constitutiva sino el camino completo
    material -> elemento -> sistema -> desplazamientos.
    """

    @staticmethod
    def _tirar_de_una_placa(theta_deg: float) -> float:
        """Placa unitaria traccionada en x; devuelve el E_x medido.

        Un solo Quad4 con carga uniforme reproduce un estado de esfuerzo
        homogéneo exactamente (el Quad4 es exacto para campos lineales), así
        que la comparación con la solución cerrada no arrastra error de
        discretización — el residuo es sólo aritmético.
        """
        mat = Orthotropic2D(**BAMBU, theta=theta_deg)
        L, t = 1.0, 0.01
        nodes = [Node(1, [0.0, 0.0]), Node(2, [L, 0.0]),
                 Node(3, [L, L]), Node(4, [0.0, L])]
        el = Quad4(1, nodes, mat, thickness=t)

        K = el.compute_global_stiffness()

        # DOFs: [ux1, uy1, ux2, uy2, ux3, uy3, ux4, uy4]
        # Tracción uniaxial: bordes libres salvo restricciones de sólido
        # rígido. ux = 0 en el borde izquierdo (nodos 1, 4); uy = 0 en el
        # nodo 1 para bloquear la traslación vertical. El resto libre, de
        # modo que la contracción de Poisson y la distorsión por
        # acoplamiento pueden desarrollarse sin coartar el estado.
        fijos = [0, 6, 1]           # ux1, ux4, uy1
        libres = [i for i in range(8) if i not in fijos]

        F = np.zeros(8)
        sigma_x = 1.0e6
        F[2] = F[4] = sigma_x * t * L / 2.0   # ux2, ux3

        U = np.zeros(8)
        U[libres] = np.linalg.solve(
            K[np.ix_(libres, libres)], F[libres])

        eps_xx = (U[2] - U[0]) / L            # alargamiento en x
        return sigma_x / eps_xx

    def test_E_x_fem_contra_analitico(self):
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                E_fem = self._tirar_de_una_placa(theta)
                ref = _E_x_analitico(theta, **BAMBU)
                self.assertLess(abs(E_fem - ref) / ref, 1e-10)


class TestDegeneracionIsotropa(unittest.TestCase):
    """Un isótropo no tiene dirección privilegiada: E_x es constante."""

    def test_E_x_independiente_de_theta(self):
        E, nu = 210.0e9, 0.3
        G = E / (2.0 * (1.0 + nu))
        for theta in ANGULOS:
            with self.subTest(theta=theta):
                m = Orthotropic2D(E1=E, E2=E, G12=G, nu12=nu, theta=theta)
                E_x, nu_xy, _ = _constantes_aparentes(m)
                self.assertLess(abs(E_x - E) / E, 1e-12)
                self.assertLess(abs(nu_xy - nu) / nu, 1e-12)

    def test_coincide_con_elastic2d_validado(self):
        """Cierra la cadena de validación: Elastic2D ya está validado contra
        Lamé y NAFEMS LE1, así que reproducirlo hereda esa validación."""
        E, nu = 210.0e9, 0.3
        G = E / (2.0 * (1.0 + nu))
        orto = Orthotropic2D(E1=E, E2=E, G12=G, nu12=nu, theta=37.0)
        iso = Elastic2D(E=E, nu=nu, hypothesis='plane_stress')
        escala = np.abs(iso.C).max()
        self.assertLess(np.abs(orto.C - iso.C).max() / escala, 1e-14)


if __name__ == "__main__":
    unittest.main()
