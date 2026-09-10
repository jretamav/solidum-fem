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

from solidum.core.domain import Domain
from solidum.core.node import Node
from solidum.elements.solid_2d import Quad4, Quad8, Quad9, Tri3, Tri6
from solidum.materials.elastic_2d import Elastic2D
from solidum.materials.orthotropic_2d import Orthotropic2D
from solidum.math.assembly import Assembler
from solidum.math.solvers import LinearSolver


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



def _malla_cuadrada(elem_cls, material, n, L, t):
    """Malla n x n del cuadrado [0,L]^2 con el elemento pedido.

    Devuelve (dominio, nodos_borde_izquierdo, nodos_borde_derecho). Soporta
    los cinco elementos sólidos 2D: los lineales usan la rejilla de corners;
    los cuadráticos añaden midnodes (y centro paramétrico en Quad9); los
    triángulos parten cada celda en dos.
    """
    cuadratico = elem_cls in (Quad8, Quad9, Tri6)
    paso = 2 if cuadratico else 1
    n_lin = n * paso + 1           # nodos por lado en la rejilla fina
    h = L / (n * paso)

    dom = Domain()
    grid = {}
    nid = 0
    for J in range(n_lin):
        for I in range(n_lin):
            # Quad8 es serendipity: no lleva centro paramétrico. Tri6 sí usa
            # ese nodo, aunque no como centro: es el midnode de la diagonal
            # que parte la celda en dos triángulos.
            if cuadratico and elem_cls is Quad8 and I % 2 == 1 and J % 2 == 1:
                continue
            nid += 1
            grid[(I, J)] = dom.add_node(nid, [I * h, J * h])

    eid = 0
    for jy in range(n):
        for ix in range(n):
            I0, J0 = ix * paso, jy * paso
            if elem_cls is Quad4:
                conn = [grid[(I0, J0)], grid[(I0 + 1, J0)],
                        grid[(I0 + 1, J0 + 1)], grid[(I0, J0 + 1)]]
                eid += 1
                dom.add_element(elem_cls(eid, conn, material, thickness=t))
            elif elem_cls is Tri3:
                c = [grid[(I0, J0)], grid[(I0 + 1, J0)],
                     grid[(I0 + 1, J0 + 1)], grid[(I0, J0 + 1)]]
                for tri in ((c[0], c[1], c[2]), (c[0], c[2], c[3])):
                    eid += 1
                    dom.add_element(elem_cls(eid, list(tri), material,
                                             thickness=t))
            elif elem_cls in (Quad8, Quad9):
                conn = [grid[(I0, J0)], grid[(I0 + 2, J0)],
                        grid[(I0 + 2, J0 + 2)], grid[(I0, J0 + 2)],
                        grid[(I0 + 1, J0)], grid[(I0 + 2, J0 + 1)],
                        grid[(I0 + 1, J0 + 2)], grid[(I0, J0 + 1)]]
                if elem_cls is Quad9:
                    conn.append(grid[(I0 + 1, J0 + 1)])
                eid += 1
                dom.add_element(elem_cls(eid, conn, material, thickness=t))
            else:  # Tri6 — diagonal del cuadrante, [v0,v1,v2,m01,m12,m20]
                sw, se = grid[(I0, J0)], grid[(I0 + 2, J0)]
                ne, nw = grid[(I0 + 2, J0 + 2)], grid[(I0, J0 + 2)]
                s_m, e_m = grid[(I0 + 1, J0)], grid[(I0 + 2, J0 + 1)]
                n_m, w_m = grid[(I0 + 1, J0 + 2)], grid[(I0, J0 + 1)]
                diag = grid[(I0 + 1, J0 + 1)]
                for tri in ((sw, se, ne, s_m, e_m, diag),
                            (sw, ne, nw, diag, n_m, w_m)):
                    eid += 1
                    dom.add_element(elem_cls(eid, list(tri), material,
                                             thickness=t))

    izq = [grid[(0, J)] for J in range(n_lin) if (0, J) in grid]
    der = [grid[(n_lin - 1, J)] for J in range(n_lin) if (n_lin - 1, J) in grid]
    return dom, izq, der


def _cargas_de_borde(elem_cls, nodos_der, L, t, sigma):
    """Reparte una tracción uniforme en cargas nodales consistentes.

    Para elementos lineales el reparto es trapezoidal (1/2 en los extremos).
    Para los cuadráticos, la carga consistente de una arista de tres nodos bajo
    presión uniforme es (1/6, 4/6, 1/6) por arista, que acumulada sobre aristas
    contiguas da el patrón 1:4:2:4:...:4:1.
    """
    F_total = sigma * L * t
    m = len(nodos_der)
    if elem_cls in (Quad8, Quad9, Tri6):
        pesos = np.zeros(m)
        for a in range(0, m - 1, 2):       # una arista por cada par de tramos
            pesos[a] += 1.0
            pesos[a + 1] += 4.0
            pesos[a + 2] += 1.0
        pesos /= pesos.sum()
    else:
        pesos = np.ones(m)
        pesos[0] = pesos[-1] = 0.5
        pesos /= pesos.sum()
    return F_total * pesos


ELEMENTOS_2D = (Quad4, Tri3, Quad8, Quad9, Tri6)

# Geometría de la probeta y nivel de tracción del ensayo virtual.
L_PROBETA = 1.0
T_PROBETA = 0.01
SIGMA_ENSAYO = 1.0e6


def _constantes_aparentes_fem(elem_cls, theta, n=2):
    """Mide (E_x, nu_xy, eta_xy_x) sobre un modelo FEM completo.

    Reproduce numéricamente el ensayo que define las constantes aparentes de
    Jones 2.8: tracción uniaxial pura, con las demás componentes de esfuerzo
    libres. Las restricciones son las mínimas para eliminar el sólido rígido
    sin coartar ni la contracción de Poisson ni la distorsión por acoplamiento
    — empotrar un borde entero destruiría justamente el efecto a medir.

    Bajo tracción uniaxial el campo es homogéneo, así que las tres constantes
    se leen de los desplazamientos de las esquinas:

        eps_xx   = u_x(esquina inferior derecha) / L
        eps_yy   = u_y(esquina superior izquierda) / L
        gamma_xy = u_y(esquina inferior derecha) / L

    El borde izquierdo tiene ux = 0, de modo que la rotación de sólido rígido
    está eliminada y el u_y del borde derecho es distorsión angular pura.

    Todos los elementos del catálogo reproducen exactamente un campo lineal,
    así que no hay error de discretización: la comparación con la forma cerrada
    es aritmética pura y admite tolerancia estricta.
    """
    material = Orthotropic2D(**BAMBU, theta=theta)
    L, t, sigma = L_PROBETA, T_PROBETA, SIGMA_ENSAYO

    dom, izq, der = _malla_cuadrada(elem_cls, material, n, L, t)
    cargas = _cargas_de_borde(elem_cls, der, L, t, sigma)

    for node in izq:
        node.fix_dof('ux', 0.0)
    izq[0].fix_dof('uy', 0.0)      # sólo una esquina: fija la traslación
    dom.generate_equation_numbers(verbose=False)

    F = np.zeros(dom.total_dofs)
    for node, valor in zip(der, cargas):
        F[node.dofs['ux']] += valor
    U = LinearSolver(Assembler(dom)).solve(F)

    inf_der, sup_izq = der[0], izq[-1]

    eps_xx = U[inf_der.dofs['ux']] / L
    eps_yy = U[sup_izq.dofs['uy']] / L
    gamma_xy = U[inf_der.dofs['uy']] / L

    return sigma / eps_xx, -eps_yy / eps_xx, gamma_xy / eps_xx


class TestConsistenciaFEM(unittest.TestCase):
    """Las tres constantes de Jones medidas sobre modelos FEM reales.

    Amplía la verificación a los **cinco** elementos sólidos 2D del catálogo.
    Cada uno tiene su propia matriz B, su propia cuadratura y su propio mapeo
    isoparamétrico; que los cinco reproduzcan la misma solución cerrada es lo
    que descarta un error localizado en un elemento concreto.

    Mide las tres constantes —no sólo E_x— porque son sensibles a errores
    distintos: `eta_xy_x` en particular es la que detecta el factor 2 del
    cortante engineering, y hasta ahora sólo se medía sobre la matriz C.
    """

    def test_E_x_contra_jones(self):
        for elem_cls in ELEMENTOS_2D:
            for theta in ANGULOS:
                with self.subTest(elemento=elem_cls.__name__, theta=theta):
                    E_fem, _, _ = _constantes_aparentes_fem(elem_cls, theta)
                    ref = _E_x_analitico(theta, **BAMBU)
                    self.assertLess(abs(E_fem - ref) / ref, 1e-9)

    def test_nu_xy_contra_jones(self):
        for elem_cls in ELEMENTOS_2D:
            for theta in ANGULOS:
                with self.subTest(elemento=elem_cls.__name__, theta=theta):
                    _, nu_fem, _ = _constantes_aparentes_fem(elem_cls, theta)
                    ref = _nu_xy_analitico(theta, **BAMBU)
                    escala = max(abs(ref), 1e-3)
                    self.assertLess(abs(nu_fem - ref) / escala, 1e-8)

    def test_eta_acoplamiento_contra_la_constitutiva(self):
        """El acoplamiento medido sobre el modelo coincide con el del material.

        `_constantes_aparentes` lo obtiene invirtiendo C; aquí se obtiene de
        los desplazamientos de un modelo con malla, ensamblaje y solver. Que
        ambos caminos coincidan blinda la recuperación del cortante engineering
        a lo largo de todo el pipeline, no sólo en la matriz constitutiva.
        """
        for elem_cls in ELEMENTOS_2D:
            for theta in ANGULOS:
                with self.subTest(elemento=elem_cls.__name__, theta=theta):
                    _, _, eta_fem = _constantes_aparentes_fem(elem_cls, theta)
                    _, _, eta_mat = _constantes_aparentes(
                        Orthotropic2D(**BAMBU, theta=theta))
                    escala = max(abs(eta_mat), 1e-3)
                    self.assertLess(abs(eta_fem - eta_mat) / escala, 1e-8)

    def test_acoplamiento_no_nulo_fuera_de_ejes_en_el_modelo(self):
        """Guardia: fuera de eje el modelo FEM debe distorsionarse de verdad.

        Sin esta comprobación, un fallo que anulase el acoplamiento haría que
        los tests anteriores comparasen dos ceros y siguieran verdes.
        """
        for elem_cls in ELEMENTOS_2D:
            for theta in (15.0, 30.0, 45.0, 60.0, 75.0):
                with self.subTest(elemento=elem_cls.__name__, theta=theta):
                    _, _, eta = _constantes_aparentes_fem(elem_cls, theta)
                    self.assertGreater(
                        abs(eta), 1e-3,
                        "el modelo no distorsiona: acoplamiento perdido")

    def test_sin_acoplamiento_en_ejes_principales_en_el_modelo(self):
        """En theta = 0 y 90 la probeta se alarga y contrae, pero no se tuerce."""
        for elem_cls in ELEMENTOS_2D:
            for theta in (0.0, 90.0):
                with self.subTest(elemento=elem_cls.__name__, theta=theta):
                    _, _, eta = _constantes_aparentes_fem(elem_cls, theta)
                    self.assertLess(abs(eta), 1e-9)



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
