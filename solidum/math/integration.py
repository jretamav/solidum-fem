# solidum_fem/solidum/math/integration.py
import numpy as np
from solidum.registry import QuadratureRegistry

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
    def get_points_tet_4():
        """Cuadratura tetraédrica de 4 puntos (Stroud T3:2-1) sobre el
        tetraedro de referencia con vértices (0,0,0), (1,0,0), (0,1,0),
        (0,0,1).

        Cuatro puntos simétricos en coordenadas baricéntricas
        ``(b, a, a, a)`` y permutaciones, donde
        ``a = (5 − √5)/20 ≈ 0.13820`` y ``b = (5 + 3√5)/20 ≈ 0.58541``.
        Peso ``1/24`` en cada uno (suma 4·1/24 = 1/6 = volumen de
        referencia). **Exacta para polinomios hasta grado 2 en (ξ, η, ζ)**
        — suficiente para integrar ``B^T D B`` del Tet10 sobre Jacobiano
        constante (el integrando es a lo sumo de grado 2 con ``B`` lineal
        en coordenadas baricéntricas).

        Para la masa consistente del Tet10 el integrando ``ρ N^T N`` es
        de grado 4 y esta regla **subintegra**; el `Tet10` declara una
        cuadratura de masa específica suficientemente exacta
        (`_MASS_QUADRATURE`).

        Referencia: Stroud A.H. (1971). *Approximate calculation of
        multiple integrals*, sec. 8.8.
        """
        import numpy as _np
        a = (5.0 - _np.sqrt(5.0)) / 20.0
        b = (5.0 + 3.0 * _np.sqrt(5.0)) / 20.0
        # Coords naturales (ξ, η, ζ) = (L_2, L_3, L_4); L_1 = 1 − ξ − η − ζ.
        # Cada punto tiene una L cerca de ``b`` y las otras tres cerca de ``a``.
        points = [(a, a, a),  # L_1 = b, otros L = a
                  (b, a, a),  # L_2 = b
                  (a, b, a),  # L_3 = b
                  (a, a, b)]  # L_4 = b
        weights = [1.0 / 24.0] * 4
        return points, weights

    @staticmethod
    def get_points_tet_15():
        """Cuadratura tetraédrica de 15 puntos (Keast 15-point) sobre el
        tetraedro de referencia.

        **Orden de exactitud 5** — exacta para polinomios hasta grado 5
        en (ξ, η, ζ). Necesaria para integrar exactamente productos
        cuadrático×cuadrático (orden 4) — caso de la matriz de masa
        consistente del Tet10. La cuadratura ``tet_4`` (orden 2) la
        subintegraría, produciendo M aproximada.

        Pesos todos positivos excepto el del centroide (negativo, común
        en cuadraturas tetraédricas de alto orden — no afecta la
        positividad de la masa total ni de la masa lumped HRZ porque la
        ponderación global ``ρ ∫ N_i N_j dV`` queda con valores netos
        positivos).

        Referencia: Keast P. (1986). "Moderate-degree tetrahedral
        quadrature formulas", *CMAME* 55(3), 339-348 (Tabla 8).
        """
        import numpy as _np
        sqrt15 = _np.sqrt(15.0)
        # Suma de pesos = 1/6 (volumen del tetraedro de referencia).
        points = []
        weights = []

        # Centroide.
        points.append((0.25, 0.25, 0.25))
        weights.append(16.0 / 135.0 / 6.0)  # = 16/810

        # Grupo A: (b, a, a, a) y permutaciones; a, w_a:
        a_A = (7.0 - sqrt15) / 34.0
        b_A = 1.0 - 3.0 * a_A
        w_A = (2665.0 + 14.0 * sqrt15) / 226800.0
        for perm in [(b_A, a_A, a_A), (a_A, b_A, a_A),
                     (a_A, a_A, b_A), (a_A, a_A, a_A)]:
            points.append(perm)
            weights.append(w_A)

        # Grupo B: (b, a, a, a) y permutaciones; a, w_b (segundo grupo):
        a_B = (7.0 + sqrt15) / 34.0
        b_B = 1.0 - 3.0 * a_B
        w_B = (2665.0 - 14.0 * sqrt15) / 226800.0
        for perm in [(b_B, a_B, a_B), (a_B, b_B, a_B),
                     (a_B, a_B, b_B), (a_B, a_B, a_B)]:
            points.append(perm)
            weights.append(w_B)

        # Grupo C: (c, c, d, d) y las 6 permutaciones únicas; c, w_c:
        c_C = (10.0 - 2.0 * sqrt15) / 40.0
        d_C = (10.0 + 2.0 * sqrt15) / 40.0
        w_C = 20.0 / 378.0 / 6.0  # = 20/2268
        # Las 6 permutaciones de (c, c, d, d) en (L1, L2, L3, L4):
        # mapeo a (ξ, η, ζ) = (L2, L3, L4) cuando L1 = 1 - L2 - L3 - L4.
        # Lista de las 6 distribuciones de "qué dos coords valen c":
        for perm_natural in [
            (c_C, c_C, d_C),  # L1=d, L2=c, L3=c, L4=d
            (c_C, d_C, c_C),  # L1=d, L2=c, L3=d, L4=c
            (d_C, c_C, c_C),  # L1=d, L2=d, L3=c, L4=c
            (c_C, d_C, d_C),  # L1=c, L2=c, L3=d, L4=d
            (d_C, c_C, d_C),  # L1=c, L2=d, L3=c, L4=d
            (d_C, d_C, c_C),  # L1=c, L2=d, L3=d, L4=c
        ]:
            points.append(perm_natural)
            weights.append(w_C)

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
QuadratureRegistry.register("tet_4", *GaussQuadrature.get_points_tet_4())
QuadratureRegistry.register("tet_15", *GaussQuadrature.get_points_tet_15())
