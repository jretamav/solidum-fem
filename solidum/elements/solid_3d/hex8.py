"""``Hex8`` — hexaedro trilineal 3D isoparamétrico (ADR 0012, Etapa 7)."""
from typing import List

import numpy as np

from solidum.core.element import Element, validate_lumping_kwarg
from solidum.core.material import Material
from solidum.core.node import Node
from solidum.elements.solid_3d._shared import (
    _compute_integrands_3d,
    _compute_kinematics_hex8,
    _det_jacobian_hex8,
    _expand_scalar_mass_3d,
    _shape_functions_hex8,
)
from solidum.math.mass_lumping import lump_hrz
from solidum.registry import ElementRegistry, QuadratureRegistry


@ElementRegistry.register
class Hex8(Element):
    """Elemento hexaédrico trilineal 3D isoparamétrico.

    Parameters
    ----------
    element_id : int
    nodes : List[Node]
        8 nodos en orden VTK_HEXAHEDRON (0-3 cara inferior antihorario
        visto desde dentro, 4-7 cara superior antihorario visto desde fuera).
    material : Material
        Material 3D (STRAIN_DIM=6).
    quadrature : str, optional
        Nombre de cuadratura desde ``QuadratureRegistry``. Default ``hex_2x2x2``
        (8 puntos). Alternativas: ``hex_3x3x3`` (27 puntos, no lineales severos),
        ``hex_1x1x1`` (1 punto, reducido — propenso a hourglass).

    Notes
    -----
    Locking volumétrico con ν → 0.5 y hourglass con integración reducida son
    limitaciones declaradas en la spec; sin mitigación implementada.
    """

    DOF_NAMES = ['ux', 'uy', 'uz']
    STRAIN_DIM = 6
    N_INTEGRATION_POINTS = 8  # default Gauss 2×2×2

    # ADR 0012 — 6 caras con normal saliente. Cada tupla son los 4 nodos
    # locales de la cara en orden tal que (a-b)×(c-b) apunte hacia fuera.
    FACE_NODES = (
        (0, 3, 2, 1),  # 0: −ζ (inferior)
        (4, 5, 6, 7),  # 1: +ζ (superior)
        (0, 1, 5, 4),  # 2: −η (frontal)
        (1, 2, 6, 5),  # 3: +ξ (derecha)
        (2, 3, 7, 6),  # 4: +η (trasera)
        (3, 0, 4, 7),  # 5: −ξ (izquierda)
    )

    def __init__(self, element_id: int, nodes: List[Node], material: Material,
                 quadrature: str = "hex_2x2x2"):
        pts, ws = QuadratureRegistry.get(quadrature)
        self.points = pts
        self.weights = ws
        self.N_INTEGRATION_POINTS = len(self.points)
        super().__init__(element_id, nodes, material)

    # ------------------------------------------------------------------
    # Contrato principal
    # ------------------------------------------------------------------

    def compute_element_state(self, u_e: np.ndarray):
        K_e = np.zeros((24, 24))
        F_int_e = np.zeros(24)

        coords = self.get_coordinate_matrix(ndim=3)

        for idx, (p, w) in enumerate(zip(self.points, self.weights)):
            xi, eta, zeta = p
            B, detJ = _compute_kinematics_hex8(xi, eta, zeta, coords)
            strain = B @ u_e

            sigma, C_tangent, new_state = self.material.compute_state(
                strain, self.state.vars[idx]
            )
            self.state.vars_trial[idx] = new_state
            self.state.stresses_trial[idx] = sigma

            K_contrib, F_contrib = _compute_integrands_3d(
                B, C_tangent, sigma, detJ, w
            )
            K_e += K_contrib
            F_int_e += F_contrib

        return K_e, F_int_e

    # ------------------------------------------------------------------
    # Salida por punto de Gauss (ADR 0012: sólidos no exponen internal_forces)
    # ------------------------------------------------------------------

    def compute_internal_forces(self, U_global: np.ndarray) -> dict:
        """Promedio del campo σ y ε sobre los puntos de Gauss (utilidad rápida).

        Para visualizar el campo no promediado úsese ``compute_gauss_state``.
        """
        gs = self.compute_gauss_state(U_global)
        return {
            'stress': gs['stress'].mean(axis=0),
            'strain': gs['strain'].mean(axis=0),
        }

    def compute_gauss_state(self, U_global: np.ndarray) -> dict:
        """Estado en cada punto de Gauss del Hex8.

        Returns
        -------
        dict con claves:
        - ``points_natural`` : ndarray (n_g, 3) con (ξ, η, ζ).
        - ``points_global``  : ndarray (n_g, 3) con (x, y, z) en globales.
        - ``strain`` : ndarray (n_g, 6) — Voigt 3D [ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz].
        - ``stress`` : ndarray (n_g, 6) — Voigt 3D [σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz].
        """
        u_e = self.get_local_displacements(U_global)
        coords = self.get_coordinate_matrix(ndim=3)
        n_g = self.N_INTEGRATION_POINTS

        nat = np.asarray(self.points, dtype=np.float64).reshape(n_g, 3)
        glb = np.zeros((n_g, 3))
        eps = np.zeros((n_g, 6))
        sig = np.zeros((n_g, 6))

        for idx, (xi, eta, zeta) in enumerate(self.points):
            N = _shape_functions_hex8(xi, eta, zeta)
            glb[idx] = N @ coords
            B, _ = _compute_kinematics_hex8(xi, eta, zeta, coords)
            strain = B @ u_e
            stress, _, _ = self.material.compute_state(strain, self.state.vars[idx])
            eps[idx] = strain
            sig[idx] = stress

        return {
            'points_natural': nat,
            'points_global': glb,
            'strain': eps,
            'stress': sig,
        }

    # ------------------------------------------------------------------
    # Cargas distribuidas consistentes
    # ------------------------------------------------------------------

    def compute_body_load(self, b: np.ndarray) -> np.ndarray:
        """Vector de cargas nodales consistente con una fuerza de cuerpo uniforme.

        Integra ``∫ N^T b dV`` sobre el elemento usando la cuadratura del
        propio elemento. ``b`` se da en coordenadas globales ``(b_x, b_y, b_z)``.

        Parameters
        ----------
        b : ndarray (3,)
            Fuerza de cuerpo por unidad de volumen (e.g. ``ρg`` para peso
            propio: ``b = (0, 0, -ρg)``).

        Returns
        -------
        ndarray (24,)
            Cargas equivalentes en los DOFs locales del elemento.
        """
        b = np.asarray(b, dtype=np.float64).reshape(3)
        coords = self.get_coordinate_matrix(ndim=3)
        f = np.zeros(24)
        for (xi, eta, zeta), w in zip(self.points, self.weights):
            N = _shape_functions_hex8(xi, eta, zeta)
            detJ = _det_jacobian_hex8(xi, eta, zeta, coords)
            if detJ <= 0.0:
                raise ValueError(f"Jacobiano no positivo en Hex8 id={self.id}.")
            factor = detJ * w
            for i in range(8):
                f[3 * i]     += N[i] * b[0] * factor
                f[3 * i + 1] += N[i] * b[1] * factor
                f[3 * i + 2] += N[i] * b[2] * factor
        return f

    def compute_face_traction(self, face: int, t_vec: np.ndarray) -> np.ndarray:
        """Vector de cargas nodales consistente con una tracción uniforme en una cara.

        Integra ``∫_Γ N^T t̄ dS`` sobre la cara cuadrilateral indicada usando
        cuadratura Gauss 2×2 sobre el cuadrado de referencia [-1,1]² con
        funciones de forma bilineales. Exacto para cara plana y tracción
        constante. ``t_vec`` se asume **constante** sobre la cara.

        Parameters
        ----------
        face : int
            Índice de cara según ``FACE_NODES`` (0..5). Numeración con
            normal saliente fijada por ADR 0012.
        t_vec : ndarray (3,)
            Tracción uniforme en globales sobre la cara.

        Returns
        -------
        ndarray (24,)
            Cargas equivalentes en los DOFs locales del elemento (no nulas
            solo en los 4 nodos de la cara).
        """
        if face not in range(6):
            raise ValueError(f"face={face} fuera de rango para Hex8 (0..5).")
        t_vec = np.asarray(t_vec, dtype=np.float64).reshape(3)

        face_local = self.FACE_NODES[face]
        X = np.array(
            [self.nodes[i].coordinates[:3] for i in face_local],
            dtype=np.float64,
        )  # 4×3

        # Cuadratura 2×2 sobre el cuadrado [-1,1]² con N_i = (1+s_i·s)(1+t_i·t)/4.
        gp = 1.0 / np.sqrt(3.0)
        s_pts = (-gp, gp, gp, -gp)
        t_pts = (-gp, -gp, gp, gp)
        s_sign = (-1.0, 1.0, 1.0, -1.0)
        t_sign = (-1.0, -1.0, 1.0, 1.0)

        f = np.zeros(24)
        for sp, tp in zip(s_pts, t_pts):
            N = np.zeros(4)
            dN_ds = np.zeros(4)
            dN_dt = np.zeros(4)
            for i in range(4):
                N[i] = 0.25 * (1.0 + s_sign[i] * sp) * (1.0 + t_sign[i] * tp)
                dN_ds[i] = 0.25 * s_sign[i] * (1.0 + t_sign[i] * tp)
                dN_dt[i] = 0.25 * t_sign[i] * (1.0 + s_sign[i] * sp)
            g_s = dN_ds @ X
            g_t = dN_dt @ X
            dA = float(np.linalg.norm(np.cross(g_s, g_t)))
            # Peso del Gauss 2×2 cada uno = 1·1; la jacobiana superficial
            # da el factor de área (g_s, g_t ya incluyen el dx/ds, dx/dt).
            for i, node_local in enumerate(face_local):
                base = 3 * node_local
                f[base]     += N[i] * t_vec[0] * dA
                f[base + 1] += N[i] * t_vec[1] * dA
                f[base + 2] += N[i] * t_vec[2] * dA
        return f

    # ------------------------------------------------------------------
    # Masa elemental (ADR 0009)
    # ------------------------------------------------------------------

    def compute_mass_matrix(self, lumping: str = "consistent") -> np.ndarray:
        """Masa del Hex8 por cuadratura del propio elemento.

        **Consistente** (default): ``M_e = ∫ρ·N^T·N dV``. Para Hex8 trilineal
        con cuadratura 2×2×2 la integral es exacta.

        **Lumped**: HRZ canónico (ADR 0009 fase 2). Para Hex8 regular coincide
        con row-sum por simetría — ``ρ·V_e/8`` en cada DOF traslacional.

        Returns
        -------
        np.ndarray, shape (24, 24)
        """
        rho = self.material.density
        coords = self.get_coordinate_matrix(ndim=3)
        M_s = np.zeros((8, 8))
        m_total = 0.0
        for (xi, eta, zeta), w in zip(self.points, self.weights):
            N = _shape_functions_hex8(xi, eta, zeta)
            detJ = _det_jacobian_hex8(xi, eta, zeta, coords)
            if detJ <= 0.0:
                raise ValueError(f"Jacobiano no positivo en Hex8 id={self.id}.")
            weight = detJ * w
            M_s += rho * np.outer(N, N) * weight
            m_total += rho * weight
        M_consistent = _expand_scalar_mass_3d(M_s)
        validate_lumping_kwarg(lumping, type(self).__name__)
        if lumping == "lumped":
            return lump_hrz(M_consistent, total_mass=m_total, n_translational_dirs=3)
        return M_consistent
