# fenix_fem/fenix/math/integration.py
import numpy as np
from fenix.registry import QuadratureRegistry

class GaussQuadrature:
    @staticmethod
    def get_points_2d_1x1():
        points = [(0.0, 0.0)]
        weights = [4.0]
        return points, weights

    @staticmethod
    def get_points_2d_2x2():
        pt = 1.0 / np.sqrt(3.0)
        points = [(-pt, -pt), ( pt, -pt), ( pt,  pt), (-pt,  pt)]
        weights = [1.0, 1.0, 1.0, 1.0]
        return points, weights

    @staticmethod
    def get_points_2d_3x3():
        """Gauss-Legendre 3×3 (9 puntos). Exacta para polinomios hasta grado 5;
        es la cuadratura por defecto de los Quad8/Quad9."""
        a = np.sqrt(3.0 / 5.0)
        xs = [-a, 0.0, a]
        ws = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
        points, weights = [], []
        for i, xi in enumerate(xs):
            for j, eta in enumerate(xs):
                points.append((xi, eta))
                weights.append(ws[i] * ws[j])
        return points, weights

    @staticmethod
    def get_points_tri_3():
        """Cuadratura triangular de 3 puntos (puntos medios de los lados),
        exacta para polinomios hasta grado 2; estándar en Tri6.

        Coordenadas baricéntricas (ξ, η) en el triángulo de referencia
        (0,0)-(1,0)-(0,1); peso 1/6 en cada uno.
        """
        points = [(0.5, 0.0), (0.5, 0.5), (0.0, 0.5)]
        weights = [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]
        return points, weights

    # ------------------------------------------------------------------
    # Cuadraturas 3D (ADR 0012 — sólidos 3D).
    # ------------------------------------------------------------------

    @staticmethod
    def get_points_hex_1x1x1():
        """Gauss-Legendre 1×1×1 sobre el cubo de referencia [-1,1]^3.

        Único punto central, peso 8.0. **Aviso**: introduce 12 modos de
        hourglass en Hex8; sin estabilización Flanagan-Belytschko (no
        implementada) los modos espurios pueden contaminar la respuesta.
        Disponible para experimentación; el default del Hex8 es 2×2×2.
        """
        points = [(0.0, 0.0, 0.0)]
        weights = [8.0]
        return points, weights

    @staticmethod
    def get_points_hex_2x2x2():
        """Gauss-Legendre 2×2×2 (8 puntos) sobre el cubo de referencia
        [-1,1]^3. Exacta para polinomios hasta grado 3 en cada dirección;
        cuadratura por defecto del Hex8."""
        pt = 1.0 / np.sqrt(3.0)
        xs = [-pt, pt]
        ws = [1.0, 1.0]
        points, weights = [], []
        for i, xi in enumerate(xs):
            for j, eta in enumerate(xs):
                for k, zeta in enumerate(xs):
                    points.append((xi, eta, zeta))
                    weights.append(ws[i] * ws[j] * ws[k])
        return points, weights

    @staticmethod
    def get_points_hex_3x3x3():
        """Gauss-Legendre 3×3×3 (27 puntos) sobre el cubo de referencia.

        Exacta para polinomios hasta grado 5 en cada dirección. Útil para
        materiales no lineales severos donde la 2×2×2 podría subintegrar el
        integrando real, y para Hex20/Hex27 cuando entren.
        """
        a = np.sqrt(3.0 / 5.0)
        xs = [-a, 0.0, a]
        ws = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
        points, weights = [], []
        for i, xi in enumerate(xs):
            for j, eta in enumerate(xs):
                for k, zeta in enumerate(xs):
                    points.append((xi, eta, zeta))
                    weights.append(ws[i] * ws[j] * ws[k])
        return points, weights

    @staticmethod
    def get_points_tet_1():
        """Cuadratura tetraédrica de 1 punto sobre el tetraedro de referencia
        con vértices (0,0,0), (1,0,0), (0,1,0), (0,0,1).

        Punto en el baricentro (1/4, 1/4, 1/4) con peso 1/6 (= volumen del
        tetraedro de referencia). Exacta para polinomios lineales en (ξ,η,ζ);
        suficiente para Tet4 (CST 3D) donde B y σ son constantes.
        """
        points = [(0.25, 0.25, 0.25)]
        weights = [1.0 / 6.0]
        return points, weights

    @staticmethod
    def get_points_tri_6():
        """Dunavant 6 puntos, orden 4 exacto sobre el triángulo de referencia.

        Necesaria para integrar exactamente productos cuadrático×cuadrático
        (orden 4 sobre ξ, η) — caso de la masa consistente del Tri6. La
        cuadratura ``tri_3`` (orden 2) la subintegraría, produciendo modos
        nulos espurios en la masa.

        Coordenadas baricéntricas: dos ternas simétricas
        (a, a, 1−2a) y (b, b, 1−2b) con sus 3 permutaciones cada una.
        Referencia: Dunavant (1985), "High degree efficient symmetrical
        Gaussian quadrature rules for the triangle", IJNME.
        """
        a1 = 0.445948490915965
        w1 = 0.111690794839005
        a2 = 0.091576213509771
        w2 = 0.054975871827661
        # En el triángulo de referencia (0,0)-(1,0)-(0,1) con
        # L1 = 1 − ξ − η, L2 = ξ, L3 = η: cada terna (L1, L2, L3) mapea a
        # (ξ, η) = (L2, L3). Sumas de pesos = 0.5 (= área del triángulo).
        triples_1 = [(a1, a1, 1.0 - 2.0 * a1),
                     (a1, 1.0 - 2.0 * a1, a1),
                     (1.0 - 2.0 * a1, a1, a1)]
        triples_2 = [(a2, a2, 1.0 - 2.0 * a2),
                     (a2, 1.0 - 2.0 * a2, a2),
                     (1.0 - 2.0 * a2, a2, a2)]
        points = [(t[1], t[2]) for t in triples_1] + \
                 [(t[1], t[2]) for t in triples_2]
        weights = [w1] * 3 + [w2] * 3
        return points, weights

# Registrar reglas automáticamente
QuadratureRegistry.register("1x1", *GaussQuadrature.get_points_2d_1x1())
QuadratureRegistry.register("2x2", *GaussQuadrature.get_points_2d_2x2())
QuadratureRegistry.register("3x3", *GaussQuadrature.get_points_2d_3x3())
QuadratureRegistry.register("tri_3", *GaussQuadrature.get_points_tri_3())
QuadratureRegistry.register("tri_6", *GaussQuadrature.get_points_tri_6())
QuadratureRegistry.register("hex_1x1x1", *GaussQuadrature.get_points_hex_1x1x1())
QuadratureRegistry.register("hex_2x2x2", *GaussQuadrature.get_points_hex_2x2x2())
QuadratureRegistry.register("hex_3x3x3", *GaussQuadrature.get_points_hex_3x3x3())
QuadratureRegistry.register("tet_1", *GaussQuadrature.get_points_tet_1())
